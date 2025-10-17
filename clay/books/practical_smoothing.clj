(ns books.practical-smoothing
  (:require [clojisr.v1.r :as rr]
            [fastmath.dev.ggplot :as gg]
            [scicloj.kindly.v4.kind :as kind]
            [tablecloth.api :as tc]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.ml.regression :as reg]
            [fastmath.matrix :as mat]
            [fastmath.stats :as stats]
            [fastmath.grid :as grid])
  (:import [java.util Calendar Locale]
           [org.apache.commons.math3.linear RealMatrix]))

;; # Practical Smoothing {.unnumbered}

;; Reimplementation and study of examples from **Practical Smoothing: The Joys of P-splines** by Paul Eilers and Brian Marx.

;; [https://psplines.bitbucket.io](https://psplines.bitbucket.io)

;; ## R packages

(rr/require-r '[base]
              '[ggplot2]
              '[colorspace]
              '[JOPS]
              '[AGD]
              '[MASS]
              '[boot]
              '[rpart]
              '[datasets]
              '[graphics]
              '[spam]
              '[grDevices]
              '[metR]
              '[gamlss.data]
              '[fds])

;; ## Common functions

(defn b-spline-base
  "Builds B-spline base function."
  ([{:keys [xl xr nseg degree] :or {xl 0.0 xr 1.0 nseg 15 degree 3}}] (b-spline-base xl xr nseg degree))
  ([xl xr nseg degree] (reg/b-spline-transformer xl xr nseg degree)))

;; Following creates a row of coordinates of B-spline bases (5+2 knots, degree 2)

(let [bfn (b-spline-base 0 1 5 2)]
  (bfn [0.25]))

(defn bbase-matrix
  "Calculates B-spline coordinates for given base function and points."
  [bsfn xs]
  (mat/rows->RealMatrix (map (comp bsfn v/vec->Vec) xs)))

(defn diffs-matrix
  "Creates differences (penalty) matrix for smoothing."
  [size degree]
  (mat/differences (mat/eye size true) degree))

(defn hdiffs-matrix
  "Creates harmonic differences (penalty) matrix for smoothing"
  [size ^double phi]
  (let [^RealMatrix D (diffs-matrix size 2)
        phi' (m/* -2.0 (m/cos phi))]
    (dotimes [n (mat/nrow D)]
      (.setEntry D n (m/inc n) phi'))
    D))

(defn range->str
  "Convert range values to formatted string (0-leading)"
  ([mx] (range->str 0 mx))
  ([mn mx] (range->str mn mx 1))
  ([mn mx step] (map (partial format "%02d") (range mn mx step))))

(defn b-spline-solve
  "Solve p-spline equation for given base B, weights W, differences D, lambda and ys"
  ([B ys] (b-spline-solve B nil nil nil ys))
  ([B D lambda ys] (b-spline-solve B nil D lambda ys))
  ([B W ys] (b-spline-solve B W nil nil ys))
  ([B W D lambda ys]
   (let [W (mat/diagonal (or W (repeat (count ys) 1.0)))
         D (or D (mat/zero (mat/ncol B) true))
         left (mat/add (mat/tmulm B (mat/mulm W B))
                       (mat/muls (mat/tmulm D D) (or lambda 1.0)))
         right (mat/mulv (mat/transpose B) (mat/mulv W (v/vec->RealVector ys)))]
     (-> (mat/solve left right)
         (v/vec->array)))))

(defn roughness
  "Roughness of coefficients"
  ([a d] (roughness a (diffs-matrix (count a) d) d))
  ([a D d]
   (let [n (count a)
         diffs (mat/mulv D (v/vec->array a))]
     (m/sqrt (m// (v/sum (v/sq diffs)) (- n d))))))

(defn fit-p-spline
  "Fit p-spline GLM"
  ([family data tr lambda] (fit-p-spline family data tr lambda 2))
  ([family data tr lambda order] (fit-p-spline family data tr nil lambda order))
  ([family data tr weights lambda order]
   (reg/glm (:y data) (:x data) {:intercept? false
                                 :weights weights
                                 :family family
                                 :transformer tr
                                 :augmentation :diffs
                                 :augmentation-params {:lambda lambda :order order}})))

(def fit-poisson (partial fit-p-spline :poisson))
(def fit-gaussian (partial fit-p-spline :gaussian))
(def fit-binomial (partial fit-p-spline :binomial))

(defn HFS
  "Find optimal lambda using HFS algorithm

  Arguments:

  * data - a dataset with x and y columns
  * fitter - fitting function, eg: `fit-gaussian`
  * tr - b-spline transformer

  Parameters:

  * `:order` - differences order (default: 2)
  * `:init-lambda` - initial lambda (default 1.0)
  * `:dispersion` - a sig2, `:mean-deviance` (default) or `pearson`
  * `:max-iters` - maximum number of iterations (default: 40)
  * `:eps` - absolute tolerance for lambda change (default: 1.0e-10)
  * `:verbose?` - show progress? (default: false)

  Returns a map containing:

  * `:fit` - fitted object
  * `:lambda` - optimal lambda
  * `:lambdas` - consecutive lambdas
  * `:iters` - number of iters"
  ([data fitter tr] (HFS data fitter tr {}))
  ([data fitter tr {:keys [order dispersion init-lambda eps max-iters weights verbose?]
                    :or {order 2 dispersion :mean-deviance init-lambda 1.0
                         eps 1.0e-10 max-iters 40 verbose? false}}]
   (loop [lambda init-lambda
          i 0
          diff ##Inf
          fit nil
          lambdas [init-lambda]]
     (if (or (> i max-iters)
             (< diff eps))
       {:fit fit :lambda lambda :iters i :lambdas lambdas}
       (let [fit (fitter data tr weights lambda order)
             ed (:effective-dimension fit)
             sig2 (get-in fit [:dispersions dispersion]) ;; mean-deviance works too
             tau2 (/ (v/sum (v/sq (v/differences (:beta fit) order))) (- ed order))
             lanew (/ sig2 tau2)]
         (when verbose? (println "iter=" i ", lambda=" lanew))
         (recur lanew (inc i) (abs (- lambda lanew)) fit (conj lambdas lanew)))))))

(defn box-product
  "See p. 71"
  [bx by]
  (for [x bx y by] (* x y)))

;; ## Introduction

;; Source file `f-ps-show.R`.

;; Illustration of:
;;
;; - B-spline basis and coefficients as peaks
;; - P-spline smoothing

(defn ps-show-create-data
  "Creates example dataset."
  ([] (ps-show-create-data nil))
  ([{:keys [xl xr n] :or {xl 0.0 xr 1.0 n 150}}]
   (let [x (m/slice-range xl xr n)
         y (-> (v/mult x 1.3)
               (v/shift 0.3)
               (v/sin)
               (v/add (repeatedly n #(r/grand 0.15)))
               (v/shift 0.3))]
     (tc/dataset {:x x :y y :id "15"}))))

(defn ps-show-peaks
  "Calculates peaks positions."
  [{:keys [nseg degree]}]
  (let [dk (/ 1.0 nseg)]
    (v/shift (v/mult (range 1 (inc (+ nseg degree))) dk)
             (- (* (inc degree) (/ dk 2.0))))))

(def ps-show-data (ps-show-create-data))

(gg/->image (gg/line (:x ps-show-data) (:y ps-show-data)))

(defn ps-show-datasets
  "Builds data for a chart."
  [data {:keys [xl xr ng lambda] :or {xl 0.0 xr 1.0 ng 500 lambda 0.1} :as config}]
  (let [bsfn (b-spline-base config)
        B (bbase-matrix bsfn (:x data))
        nb (mat/ncol B)
        xg (m/slice-range xl xr ng)
        Bg (bbase-matrix bsfn xg)
        xa (ps-show-peaks config)
        D (diffs-matrix nb 2)
        a (b-spline-solve B D lambda (:y data))
        z (v/vec->array (mat/mulv Bg (v/vec->RealVector a)))]
    {:Zf (tc/dataset {:x xg :y z :id "01"})
     :C (tc/dataset {:x xa :y a :id "01"})
     :Bf (tc/dataset {:x (mapcat identity (repeat nb xg))
                      :y (mapcat v/vec->seq (mat/cols (mat/mulm Bg (mat/diagonal a))))
                      :id (mapcat identity (map (partial repeat ng)
                                                (range->str 1 (inc nb))))})
     :nb nb}))

(defn ps-show-plot
  "Creates f-ps-show.R plot."
  [data config]
  (let [{:keys [Zf C Bf nb]} (ps-show-datasets data config)
        pal (colorspace/rainbow_hcl nb :c 80 :l 80 :start 10 :end 350)]
    (rr/r+ (ggplot2/ggplot Bf (ggplot2/aes :x :x :y :y :group :id :color :id))
           (ggplot2/geom_line :size 0.7)
           (ggplot2/geom_hline :yintercept 0.0 :size 0.3)
           (ggplot2/geom_line :data data :color "grey60")
           (ggplot2/geom_point :data data :color "grey60" :size 0.9)
           (ggplot2/geom_line :data Zf :size 1 :color "blue")
           (ggplot2/geom_point :data C :color pal :size 2.5)
           (ggplot2/xlab "")
           (ggplot2/ylab "")
           (ggplot2/ggtitle "Figure 1.1")
           (ggplot2/theme_light)
           (ggplot2/theme :legend.position "none")
           (ggplot2/scale_x_continuous :breaks (range 0.0 1.01 0.2))
           (ggplot2/scale_color_manual :values pal))))

(def ps-show-config {:n 150
                   :nseg 15
                   :degree 3
                   :ng 500
                   :lambda 0.1
                   :xl 0.0
                   :xr 1.0})

(gg/->image (ps-show-plot ps-show-data ps-show-config)
            {:width 800 :height 600})

;; ## Bases, Penalties, and Likelihoods

;; ### Linear and Polynomial Regression

;; Source file: `f-air-wind.R`

;; Linear and polynomial regression comparison.

(rr/data 'airquality 'JOPS)

(def airquality-data (-> (tc/dataset {:x (rr/r->clj (rr/r$ airquality 'Wind))
                                    :y (rr/r->clj (rr/r$ airquality 'Ozone))})
                       (tc/drop-missing)))

(defn air-wind-models
  [data]
  {:poly-lm (reg/lm (:y data) (:x data) {:transformer (reg/polynomial-transformer (:x data) 2)})
   :lin-lm (reg/lm (:y data) (:x data))})

(defn air-wind-plot
  [data]
  (let [{:keys [poly-lm lin-lm]} (air-wind-models data)
        xs (m/slice-range (:x data) 50)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point)
           (ggplot2/geom_smooth :data (tc/dataset {:x xs :y (map poly-lm xs)}) :colour "red" :se false)
           (ggplot2/geom_smooth :data (tc/dataset {:x xs :y (map lin-lm xs)}) :colour "blue" :se false :lty 2)
           (ggplot2/theme_light)
           (ggplot2/xlab "Wind speed (mph)")
           (ggplot2/ylab "Ozone concetntraion (ppb)")
           (ggplot2/ggtitle "Figure 2.1: New York air quality"))))

(gg/->image (air-wind-plot airquality-data))

;; ---

;; Source file: `f-motpol1.R`

;; Two polynomial (degree 9) fits for two datasets:

;; - blue line - whole dataset
;; - red line - times > 5ms

(rr/data 'mcycle 'MASS)

(def mcycle-data (tc/dataset {:times (rr/r->clj (rr/r$ mcycle 'times))
                            :accel (rr/r->clj (rr/r$ mcycle 'accel))}))

(defn motpol1-models
  [data]
  (let [mc1 (tc/select-rows data (fn [row] (> (:times row) 5)))
        lm1 (reg/lm (:accel mc1) (:times mc1)
                    {:transformer (reg/polynomial-transformer (:times mc1) 9)})
        xg1 (m/slice-range (:times mc1) 100)
        mc2 (tc/select-rows data (fn [row] (pos? (:times row))))
        lm2 (reg/lm (:accel mc2) (:times mc2)
                    {:transformer (reg/polynomial-transformer (:times mc2) 9)})
        xg2 (m/slice-range (:times mc2) 100)]
    {:F1 (tc/dataset {:times xg1 :accel (map lm1 xg1)})
     :F2 (tc/dataset {:times xg2 :accel (map lm2 xg2)})}))

(defn motpol1-plot
  [data]
  (let [{:keys [F1 F2]} (motpol1-models data)]
    (rr/r+ (ggplot2/ggplot mcycle (ggplot2/aes :x :times :y :accel))
           (ggplot2/geom_point)
           (ggplot2/xlab "Time (ms)")
           (ggplot2/ylab "Acceleration (g)")
           (ggplot2/ggtitle "Figure 2.2: Polynomial fits to motorcycle helmet data")
           (ggplot2/geom_hline :yintercept 0.0)
           (ggplot2/geom_line :data F1 :color "red" :size 1 :lty 6)
           (ggplot2/geom_line :data F2 :color "blue" :size 1 :lty 1)
           (ggplot2/theme_light))))

(gg/->image (motpol1-plot mcycle-data))

;; ### B-splines

;; Source file: `f-persp.R`

;; B-splines in two perspectives: vertical and plotted on the top of each other

(defn persp-base
  ([] (persp-base nil))
  ([{:keys [xmin xmax ng deg1 ndx]
     :or {xmin 0.0 xmax 4.0 ng 500 deg1 3 ndx 4}}]
   (let [bbase-fn (b-spline-base xmin xmax ndx deg1)
         xg (m/slice-range xmin xmax ng)
         Bg (bbase-matrix bbase-fn xg)
         nb1 (mat/ncol Bg)
         Bsc1 (mat/add Bg (mat/outer (repeat ng 1.0) (range 1 (inc nb1))))]
     {:Bg Bg :Bsc1 Bsc1 :nb1 nb1 :xg xg})))

(defn perps-get-row
  [B k]
  (let [bk2 (v/vec->seq (mat/row B k))
        bk1 (map-indexed (fn [v i] (+ v i 1)) bk2)]
    [bk1 (map (fn [v] (if (< v 0.001) ##NaN v)) bk2)]))

(defn persp-data
  ([] (persp-data nil))
  ([{:keys [ng k] :or {ng 500 k 159} :as config}]
   (let [{:keys [Bg Bsc1 nb1 xg] :as bbase} (persp-base config)
         x (mapcat identity (repeat nb1 xg))
         id (mapcat identity (map (partial repeat ng) (range->str 1 (inc nb1))))
         Bf1 (tc/dataset {:x x :id id :y (mapcat v/vec->seq (mat/cols Bsc1))})
         Bf2 (-> (tc/dataset {:x x :id id :y (mapcat v/vec->seq (mat/cols Bg))})
                 (tc/drop-rows (fn [row] (< (:y row) 0.0001))))
         [bk1 bk2] (perps-get-row Bg k)
         xk (nth xg k)
         factor (range->str 1 (inc nb1))
         Fk1 (tc/dataset {:x xk :y bk1 :id factor})
         Fk2 (-> (tc/dataset {:x xk :y bk2 :id factor})
                 (tc/drop-missing))]
     (assoc bbase :Bf1 Bf1 :Bf2 Bf2 :Fk1 Fk1 :Fk2 Fk2 :xk xk))))

(defn persp-plot-stub
  [B F xk nb1]
  (rr/r+ (ggplot2/ggplot B (ggplot2/aes :x :x :y :y :group :id :color :id))
         (ggplot2/geom_line :size 1.2)
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)
         (ggplot2/geom_vline :xintercept xk :lty 2)
         (ggplot2/geom_point (ggplot2/aes :color :id) :size 3 :data F)
         (ggplot2/theme :legend.position "none")
         (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (inc nb1) :start 10 :end 350))))

(defn persp-plot
  ([] (persp-plot nil))
  ([config]
   (let [{:keys [Bf1 Fk1 Bf2 Fk2 nb1 xk]} (persp-data config)]
     [(rr/r+ (persp-plot-stub Bf1 Fk1 xk nb1)
             (ggplot2/ggtitle "Figure 2.3: Perspective view"))
      (rr/r+ (persp-plot-stub Bf2 Fk2 xk nb1)
             (ggplot2/ggtitle "Figure 2.3: Columns of a B-spline basis")
             (ggplot2/geom_hline :yintercept 0 :size 0.3))])))

(let [[p1 p2] (persp-plot)]
  (kind/table
   [[(gg/->image p1)
     (gg/->image p2)]]))

;; ---

;; Source file: `f-B-lin-quad.R`

;; Linear and quadratic B-spline basis

(defn b-lin-quad-data
  [xg degree]
  (let [ng (count xg)
        bbase (b-spline-base 0.0 4.0 4 degree)
        Bf (bbase-matrix bbase xg)
        nb (mat/ncol Bf)]
    [nb (-> (tc/dataset {:x (mapcat identity (repeat nb xg))
                         :y (mapcat v/vec->seq (mat/cols Bf))
                         :id (mapcat identity (map (partial repeat ng) (range->str 1 (inc nb))))})
            (tc/drop-rows (fn [row] (< (:y row) 0.0001))))]))

(defn b-lin-quad-plot-stub
  [Bf nb]
  (rr/r+ (ggplot2/ggplot Bf (ggplot2/aes :x :x :y :y :group :id :color :id))
         (ggplot2/geom_line :size 1.2)
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)
         (ggplot2/theme :legend.position "none")
         (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (inc nb) :start 10 :end 350))))

(defn b-lin-quad-plot
  []
  (let [xg (m/slice-range 0.0 4.0 500)
        [nb1 Bf1] (b-lin-quad-data xg 1)
        [nb2 Bf2] (b-lin-quad-data xg 2)]
    [(rr/r+ (b-lin-quad-plot-stub Bf1 nb1)
            (ggplot2/ggtitle "Figure 2.4: Linear B-splines"))
     (rr/r+ (b-lin-quad-plot-stub Bf2 nb2)
            (ggplot2/ggtitle "Figure 2.4: Quadratic B-splines"))]))

(let [[p1 p2] (b-lin-quad-plot)]
  (kind/table
   [[(gg/->image p1)
     (gg/->image p2)]]))

;; ---

;; Source file: `f-mot-bsp.R`

;; Two B-spline (segments: 5, degree: 3) fits for two datasets:

;; - blue line - whole dataset
;; - red line - times > 5ms

(def mcycle2-data (tc/rename-columns mcycle-data {:times :x :accel :y}))

(defn mot-bsp-data
  [data]
  (let [x (:x data)
        y (:y data)
        [xl xr] (stats/extent x)
        base (b-spline-base xl xr 5 3)
        B (bbase-matrix base x)
        a (b-spline-solve B y)
        asel (b-spline-solve B (map #(if (> % 5) 1.0 0.0) x) y)
        xg (m/slice-range x 1000)
        Bg (bbase-matrix base xg)]
    {:Zf1 (tc/dataset {:x xg :y (v/vec->seq (mat/mulv Bg a))})
     :Zf2 (-> (tc/dataset {:x xg :y (v/vec->seq (mat/mulv Bg asel))})
              (tc/select-rows (fn [row] (> (:x row) 5))))}))

(defn mot-bsp-plot
  [data]
  (let [{:keys [Zf1 Zf2]} (mot-bsp-data data)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :size 1)
           (ggplot2/geom_hline :yintercept 0 :size 0.3)
           (ggplot2/geom_line :data Zf1 :size 1 :color "blue")
           (ggplot2/geom_line :data Zf2 :size 1 :color "red" :lty 6)
           (ggplot2/xlab "Time (ms")
           (ggplot2/ylab "Acceleration (g)")
           (ggplot2/ggtitle "Figure 2.5: Motorcycle helmet impact data")
           (ggplot2/ylim [-150 100])
           (ggplot2/theme_light))))

(gg/->image (mot-bsp-plot mcycle2-data))

;; ---

;; Source file: `f-mot-bsize.R`

;; Small (10) and large (20) basis

(defn mot-bsize-calculate
  [y x xg nseg]
  (let [base (b-spline-base 0.0 60.0 nseg 3)
        B (bbase-matrix base x)
        Bg (bbase-matrix base xg)
        a (b-spline-solve B y)]
    (v/vec->seq (mat/mulv Bg a))))

(defn mot-bsize-data
  [data]
  (let [x (:x data)
        y (:y data)
        xg (m/slice-range x 1000)]
    {:Z1 (tc/dataset {:x xg :y (mot-bsize-calculate y x xg 10)})
     :Z2 (tc/dataset {:x xg :y (mot-bsize-calculate y x xg 20)})}))

(defn mot-bsize-plot
  [data]
  (let [{:keys [Z1 Z2]} (mot-bsize-data data)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :size 1)
           (ggplot2/geom_hline :yintercept 0 :size 0.3)
           (ggplot2/geom_line :data Z1 :size 1.2 :color "red" :linetype 2)
           (ggplot2/geom_line :data Z1 :size 0.6 :color "red" :linetype 1)
           (ggplot2/geom_line :data Z2 :size 1 :color "blue" :linetype 1)
           (ggplot2/xlab "Time (ms")
           (ggplot2/ylab "Acceleration (g)")
           (ggplot2/ggtitle "Figure 2.6: Motorcycle helmet impact data")
           (ggplot2/ylim [-200 100])
           (ggplot2/theme_light))))

(gg/->image (mot-bsize-plot mcycle2-data))

;; ---

;; Source file `f-bsize.R`

;; Small (5) and large bases (15)

(defn bsize-create-data
  "Creates example dataset."
  []
  (let [x (m/slice-range 0 1 150)
        y (-> (v/mult x 1.2)
              (v/shift 0.3)
              (v/sin)
              (v/add (repeatedly 150 #(r/grand 0.1)))
              (v/shift 0.3))]
    (tc/dataset {:x x :y y :id "05"})))

(defn bsize-datasets
  [data xg ndx]
  (let [base (b-spline-base 0 1 ndx 3)
        B (bbase-matrix base (:x data))
        nb (mat/ncol B)
        Bg (bbase-matrix base xg)
        a (b-spline-solve B (:y data))
        z (mat/mulv Bg a)
        Bsc (mat/mulm Bg (mat/diagonal a))]
    {:nb nb
     :Zf (tc/dataset {:x xg :y z :id "01"})
     :Bf (-> (tc/dataset {:x (mapcat identity (repeat nb xg))
                          :y (mapcat v/vec->seq (mat/cols Bsc))
                          :id (mapcat identity (map (partial repeat 500) (range->str 1 (inc nb))))})
             (tc/drop-rows (fn [row] (< (abs (:y row)) 0.0001))))}))

(defn bsize-plot
  [data xg ndx title]
  (let [{:keys [nb Zf Bf]} (bsize-datasets data xg ndx)]
    (rr/r+ (ggplot2/ggplot Bf (ggplot2/aes :x :x :y :y :group :id :color :id))
           (ggplot2/geom_line :size 0.7)
           (ggplot2/ggtitle (str "Figure 2.7: " title))
           (ggplot2/geom_hline :yintercept 0 :size 0.3)
           (ggplot2/geom_line :data Zf :size 1 :color "blue")
           (ggplot2/geom_point :data data :color "grey60" :size 0.8)
           (ggplot2/xlab "")
           (ggplot2/ylab "")
           (ggplot2/theme_light)
           (ggplot2/theme :legend.position "none")
           (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (inc nb) :start 10 :end 350)))))

(defn bsize-plots
  [data]
  (let [xg (m/slice-range 0 1 500)]
    [(bsize-plot data xg 5 "Small basis")
     (bsize-plot data xg 15 "Large basis")]))

(kind/table
 (let [[p1 p2] (bsize-plots (bsize-create-data))]
   [[(gg/->image p1)
     (gg/->image p2)]]))

;; ---

;; Source file: `f-brough2.R`

;; Roughness for different coefficients `a`

(defn brough2-data
  [id n]
  (case (int id)
    1 (repeatedly n rand)
    2 (v/add (v/mult (brough2-data 1 n) 0.2)
             (v/mult (v/sin (v/mult (brough2-data 3 n) 2)) 0.8))
    3 (v/div (range 1 (inc n)) n)
    4 (repeat n 1)))

(defn brough2-datasets
  [B u id n knots groups]
  (let [a (brough2-data id n)
        r (format "r=%3.2f" (roughness a 1))
        Bsc (mat/mulm B (mat/diagonal a))]
    {:r r
     :Zf (tc/dataset {:x u :y (v/vec->array (mat/mulv B (v/vec->array a))) :id "01"})
     :Bf (-> (tc/dataset {:x (mapcat identity (repeat n u))
                          :y (mapcat v/vec->seq (mat/cols Bsc))
                          :id groups})
             (tc/drop-rows (fn [row] (< (:y row) 0.0001))))
     :Fa (tc/dataset {:x knots :y a :id "01"})}))

(defn brough2-plot
  [n {:keys [r Zf Bf Fa]}]
  (rr/r+ (ggplot2/ggplot Bf (ggplot2/aes :x :x :y :y :group :id :color :id))
         (ggplot2/geom_line :size 0.7)
         (ggplot2/ggtitle r)
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/geom_line :data Zf :size 1 :col "blue")
         (ggplot2/geom_point :data Fa :color "red" :size 2 :shape 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)
         (ggplot2/theme :legend.position "none")
         (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (inc n) :start 10 :end 350))))

(defn brough2-plots
  [nseg]
  (let [u (m/slice-range 0 1 200)
        base (b-spline-base 0 1 nseg 3)
        B (bbase-matrix base u)
        n (mat/ncol B)
        knots (v/div (v/shift (range 1 (inc n)) -2) nseg)
        groups (mapcat identity (map (partial repeat 200) (range->str 1 (inc n))))]
    (for [id (range 4)
          :let [data (brough2-datasets B u (inc id) n knots groups)]]
      (brough2-plot n data))))

(kind/table
 (let [[p1 p2 p3 p4] (map gg/->image (brough2-plots 10))]
   [[p1 p2]
    [p3 p4]]))

;; ### Penalized Least Squares

;; Source files: `f-d1pen.R`, `f-d2pen.R`

;; Roughness and stddev of residuals

(defn dxpen-datasets
  [pre u B Bu D degree lambda y knots]
  (let [a (b-spline-solve B D lambda y)
        z (mat/mulv Bu a)
        mu (mat/mulv B a)
        s (v/sqrt (stats/mean (v/sq (v/sub y mu))))
        r (roughness a D degree)]
    {:title (str pre "λ=" lambda " | s=" (m/approx s 2) " | r=" (m/approx r 2))
     :Zf (tc/dataset {:x u :y z :id "01"})
     :Fa (tc/dataset {:x knots :y a :id "01"})}))

(defn dxpen-plot
  [data {:keys [title Zf Fa]}]
  (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y) :ylim [0 1.5])
         (ggplot2/geom_point :color "grey40")
         (ggplot2/geom_line :data Zf :color "blue" :size 1)
         (ggplot2/geom_point :data Fa :color "red" :size 2 :shape 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/ggtitle title)
         (ggplot2/theme_light)))

(defn dxpen-plots
  [pre nseg degree lambdas]
  (let [x (repeatedly 50 rand)
        y (-> (v/add (v/sin (v/mult x 2.5))
                     (repeatedly 50 #(r/grand 0.1)))
              (v/shift 0.2))
        u (m/slice-range 0 1 200)
        base (b-spline-base 0 1 nseg 3)
        Bu (bbase-matrix base u)
        B (bbase-matrix base x)
        n (mat/ncol B)
        D (diffs-matrix n degree)
        knots (v/div (v/shift (range 1 (inc n)) -2) nseg)
        dataset (tc/dataset {:x x :y y})]
    (for [lambda lambdas
          :let [data (dxpen-datasets pre u B Bu D degree lambda y knots)]]
      (dxpen-plot dataset data))))

;; Differences degree = 1

(kind/table
 (let [[p1 p2 p3 p4] (map gg/->image (dxpen-plots "Figure 2.9: " 20 1 [0.1 1 10 100]))]
   [[p1 p2] [p3 p4]]))

;; Differences degree = 2

(kind/table
 (let [[p1 p2 p3 p4] (map gg/->image (dxpen-plots "Figure 2.10: " 20 2 [0.1 5 500 10000]))]
   [[p1 p2] [p3 p4]]))

;; ### Interpolation and Extrapolation

;; Source file: `f-extrapol1.R`

(defn extrapol1-data
  []
  (let [x (->> (repeatedly 50 rand)
               (filter (fn [x] (or (< 0.2 x 0.4)
                                  (< 0.6 x 0.8)))))
        y (-> (v/add (v/sin (v/mult x 2.5))
                     (repeatedly 50 #(r/grand 0.05)))
              (v/shift 0.2))]
    (tc/dataset {:x x :y y})))

(defn extrapol1-datasets
  [xg B Bg nb y degree knots]
  (let [D (diffs-matrix nb degree)
        a (b-spline-solve B D 1 y)
        z (mat/mulv Bg a)]
    {:Fa (tc/dataset {:x knots :y a})
     :Zf (tc/dataset {:x xg :y z})}))

(def extrapol1-titles {1 "First" 2 "Second"})

(defn extrapol1-plot
  [degree data {:keys [Fa Zf]}]
  (rr/r+ (ggplot2/ggplot :mapping (ggplot2/aes :x :x :y :y) :ylim [nil 1.5])
         (ggplot2/geom_line :data Zf :size 0.6 :col "blue")
         (ggplot2/ggtitle (str "Figure 2.11: " (extrapol1-titles degree) " differences"))
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/geom_point :data data :color "darkgrey")
         (ggplot2/geom_point :data Fa :color "red" :shape 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)))

(defn extrapol1-plots
  [ndx degrees]
  (let [data (extrapol1-data)
        base (b-spline-base 0 1 ndx 3)
        B (bbase-matrix base (:x data))
        nb (mat/ncol B)
        xg (m/slice-range 0 1 500)
        Bg (bbase-matrix base xg)
        knots (v/div (v/shift (range 1 (inc nb)) -2) ndx)]
    (for [degree degrees
          :let [datasets (extrapol1-datasets xg B Bg nb (:y data) degree knots)]]
      (extrapol1-plot degree data datasets))))

(kind/table
 [(map gg/->image (extrapol1-plots 25 [1 2]))])

;; ### Derivatives

;; Source file: `f-slope-height.R`

(rr/data 'boys7482 'AGD)

(def boys7482-data-hgt (-> (tc/dataset {:x (rr/r->clj (rr/r$ boys7482 'age))
                                      :y (rr/r->clj (rr/r$ boys7482 'hgt))})
                         (tc/drop-missing)
                         (tc/select-rows (fn [row] (< 0 (:x row) 20)))
                         (tc/shuffle)
                         (tc/head 1000)))

(defn slope-height-derivative
  [xg n mn mx nseg bdeg deriv-deg]
  (let [base (b-spline-base mn mx nseg (dec bdeg))]
    (-> (bbase-matrix base xg)
        (mat/mulm (diffs-matrix n deriv-deg))
        (mat/muls (/ nseg (- mx mn))))))

(defn slope-height-dataset
  [data nseg bdeg]
  (let [[mn mx] (stats/extent (:x data))
        base (b-spline-base mn mx nseg bdeg)
        B (bbase-matrix base (:x data))
        n (mat/ncol B)
        xg (range mn (+ mx m/EPSILON) 0.1)
        Bg (bbase-matrix base xg)
        B1 (slope-height-derivative xg n mn mx nseg bdeg 1)
        D (diffs-matrix n 2)
        a (b-spline-solve B D 100 (:y data))
        z (mat/mulv Bg a)
        g (mat/mulv B1 a)]
    (tc/dataset {:x xg :z z :g g})))

(defn slope-height-plot-1
  [data dataset]
  (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point :color "darkgrey" :size 1)
         (ggplot2/xlab "Age")
         (ggplot2/ylab "Height (cm)")
         (ggplot2/xlim [0 20])
         (ggplot2/geom_line :data dataset (ggplot2/aes :x :x :y :z) :size 1 :color "blue")
         (ggplot2/ggtitle "Figure 2.12: Heights of Dutch boys")
         (ggplot2/theme_light)))

(defn slope-height-plot-2
  [dataset]
  (rr/r+ (ggplot2/ggplot dataset (ggplot2/aes :x :x :y :g))
         (ggplot2/xlab "Age")
         (ggplot2/ylab "Growth speed (cm/y)")
         (ggplot2/ylim [0 25])
         (ggplot2/xlim [0 20])
         (ggplot2/geom_line  :size 1 :color "blue")
         (ggplot2/ggtitle "Figure 2.12: Growth speed of Dutch boys")
         (ggplot2/theme_light)))

(defn slope-height-plots
  [data nseg bdeg]
  (let [dataset (slope-height-dataset data nseg bdeg)]
    [(slope-height-plot-1 data dataset)
     (slope-height-plot-2 dataset)]))

(kind/table
 [(map gg/->image (slope-height-plots boys7482-data-hgt 50 3))])

;; ### The Effective Dimension

;; Source file: `f-Peffdim.R`

(defn peffdim-fit
  [data bbase lambin order]
  (map (fn [^double lambda]
         (:effective-dimension (fit-gaussian data bbase lambda order))) lambin))

(defn peffdim-fit-orders
  [data bbase lambin orders]
  (pmap (partial peffdim-fit data bbase lambin) orders))

(defn peffdim-dataset
  [data bbase]
  (let [llambs (range -5 5.01 0.2)
        lambin (map (fn [^double v] (m/pow 10.0 v)) llambs)
        [ed1 ed2 ed3] (peffdim-fit-orders data bbase lambin [1 2 3])]
    (tc/concat (tc/dataset {:llambs llambs :ED ed1 :order "1"})
               (tc/dataset {:llambs llambs :ED ed2 :order "2"})
               (tc/dataset {:llambs llambs :ED ed3 :order "3"}))))

(defn peffdim-plot
  [data]
  (let [bbase (reg/b-spline-transformer (:x data) 20 3)]
    (rr/r+ (ggplot2/ggplot (peffdim-dataset data bbase)
                           (ggplot2/aes :x :llambs :y :ED :color :order :linetype :order))
           (ggplot2/geom_line :size 1.5)
           (ggplot2/labs :x "log10(λ)")
           (ggplot2/geom_hline :yintercept 0)
           (ggplot2/ggtitle "Figure 2.13: Effective dimensions, across penalty order")
           (ggplot2/theme_light))))

(gg/->image (peffdim-plot mcycle2-data))

;; ### Standard Errors

;; Source file: `f-se.R`

(defn se-dataset
  [data]
  (let [[xl xr] (stats/extent (:x data))
        res (fit-gaussian data (reg/b-spline-transformer xl xr 20 3) 0.8 2)]
    (-> (tc/dataset (map (fn [^double x] (-> (select-keys (res x true) [:fit :stderr])
                                            (assoc :x x)))
                         (m/slice-range xl xr 100)))
        (tc/map-columns :selow [:fit :stderr] (fn [^double f ^double s] (- f (* 2 s))))
        (tc/map-columns :seup [:fit :stderr] (fn [^double f ^double s] (+ f (* 2 s)))))))

(defn se-plot
  [data]
  (let [dataset (se-dataset data)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :color "grey50")
           (ggplot2/geom_line :data dataset (ggplot2/aes :x :x :y :fit) :color "blue" :size 1)
           (ggplot2/geom_line :data dataset (ggplot2/aes :x :x :y :selow) :color "red" :size 0.8 :linetype 2)
           (ggplot2/geom_line :data dataset (ggplot2/aes :x :x :y :seup) :color "red" :size 0.8 :linetype 2)
           (ggplot2/ggtitle "Figure 2.14: P-spline fit to motorcycle helmet data")
           (ggplot2/xlab "Time (ms)")
           (ggplot2/ylab "Acceleration (g)")
           (ggplot2/theme_light))))

(gg/->image (se-plot mcycle2-data))

;; ### P-splines as a Parametric Model

;; Source file: `f-bcoeff.R`

(def bcoeff-data
  (let [x (repeatedly 50 rand)
        n (repeatedly 50 r/grand)]
    (tc/dataset {:x x
                 :y (-> x (v/mult 1.2) (v/shift 0.3) (v/sin) (v/shift 0.3) (v/add (v/mult n 0.1)))})))

(defn bcoeff-base
  [data degree]
  (let [xg (m/slice-range 0 1 500)
        bbase (b-spline-base 0 1 10 degree)
        B (bbase-matrix bbase (:x data))
        Bg (bbase-matrix bbase xg)
        nb (mat/ncol B)
        knots (v/div (v/shift (range 1 (inc nb)) (if (m/one? degree) -1 -2)) 10)
        D (mat/differences (mat/eye nb true) 2)
        ids (mapcat identity (map (partial repeat 500) (range->str 1 (inc nb))))]
    {:xg xg :bbase bbase :B B :Bg Bg :nb nb :knots knots :D D :ids ids}))

(defn bcoeff-datasets
  [data {:keys [xg B Bg nb knots D ids]} lambda]
  (let [a (b-spline-solve B D lambda (:y data))
        z (v/vec->seq (mat/mulv Bg a))
        Bsc (mat/mulm Bg (mat/diagonal a))]
    {:A (tc/dataset {:x knots :y a :id "01"})
     :Zf (tc/dataset {:x xg :y z :id "01"})
     :Bf (-> (tc/dataset {:x (mapcat identity (repeat nb xg))
                          :y (mapcat v/vec->seq (mat/cols Bsc))
                          :id ids})
             (tc/drop-rows (fn [row] (< (:y row) 1.0e-4))))
     :nb nb}))

(defn bcoeff-plot
  [{:keys [Bf Zf A nb]} title]
  (rr/r+ (ggplot2/ggplot Bf (ggplot2/aes :x :x :y :y :group :id :color :id))
         (ggplot2/geom_line :size 0.7)
         (ggplot2/ggtitle title)
         (ggplot2/geom_hline :yintercept 0.0 :size 0.3)
         (ggplot2/geom_line :data Zf :size 1 :color "blue")
         (ggplot2/geom_point :data A :color "red" :size 1.5)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)
         (ggplot2/theme :legend.position "none")
         (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (inc nb) :start 10 :end 350))))

(defn bcoeff-plots
  [data degree no]
  (let [base (bcoeff-base data degree)
        dss (map (partial bcoeff-datasets data base) [0.1 3.0])]
    (->> ["Wiggly curve" "Smooth curve"]
         (map (fn [t] (str "Figure " no ": " t)))
         (map bcoeff-plot dss))))

(kind/table
 [(map gg/->image (bcoeff-plots bcoeff-data 3 "2.15"))])

;; ---

;; Source file: `f-bcoeff-lin.R`

;; Same as above, with `degree=1`

(kind/table
 [(map gg/->image (bcoeff-plots bcoeff-data 1 "2.16"))])

;; ### Equivalent Kernels

;; Source file: `f-eff-kernel.R`

(defn eff-kernels-dataset
  [xs E D2 lambda ids]
  (let [P (mat/muls D2 lambda)
        H (mat/inverse (mat/add E P))
        hs (map (partial mat/col H) (m/slice-range 0 200 5))]
    (tc/dataset {:x (mapcat identity (repeat 5 xs))
                 :y (mapcat v/vec->seq hs)
                 :id ids})))

(defn eff-kernels-plot
  [dataset lambda]
  (rr/r+ (ggplot2/ggplot dataset (ggplot2/aes :x :x :y :y :color :id))
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/geom_line)
         (ggplot2/ggtitle (str "Figure 2.17: λ=" lambda))
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/theme_light)
         (ggplot2/theme :legend.position "none")))

(defn eff-kernels-plots
  [lambdas]
  (let [E (mat/eye 201 true)
        D (mat/differences E 2)
        D2 (mat/tmulm D D)
        ids (mapcat identity (map (partial repeat 201) (range->str 1 6)))
        xs (range 1 202)]
    (for [lambda lambdas
          :let [ds (eff-kernels-dataset xs E D2 lambda ids)]]
      (eff-kernels-plot ds lambda))))

(kind/table
 (->> (eff-kernels-plots [1.0 100.0 1e4 1e6])
      (map gg/->image)
      (partition 2)))

;; ### Smoothing of a Non-normal Response

;; Based on GLM

;; #### Poisson Smoothing

;; Source file: `f-coal-smooth.R`

(rr/data 'coal 'boot)

(def coal-histogram-data (let [hist (->> (rr/r$ 'coal 'date)
                                       (rr/r->clj)
                                       (map int)
                                       (frequencies))
                             year (range 1851 1963)]
                         (-> (tc/dataset {:Year year :Count (map #(get hist % 0) year)})
                             (tc/rename-columns {:Year :x :Count :y}))))

(defn coal-smooth-datasets
  [data]
  (let [[xl xr] (stats/extent (:x data))
        bbase (reg/b-spline-transformer xl xr 20 3)
        fit1 (fit-poisson data bbase 1)
        fit2 (fit-poisson data bbase 100)
        n (count (:beta fit1))
        knots (v/shift (v/mult (v/shift (range 1 (inc n)) -2) (/ (- xr xl) 20)) xl)
        grid (m/slice-range xl xr 100)
        p1 (map #(fit1 % true) grid)
        p2 (map #(fit2 % true) grid)]
    {:F2 (tc/dataset {:xg grid
                      :yg1 (map :fit p1) :yg2 (map :fit p2)
                      :eta1 (map :link p1) :eta2 (map :link p2)})
     :F3 (tc/dataset {:knots knots :pcoef1 (:beta fit1) :pcoef2 (:beta fit2)})}))

(defn coal-smooth-plot-1
  [F1 F2]
  (rr/r+ (ggplot2/ggplot)
         (ggplot2/geom_point :data F1 (ggplot2/aes :x :x :y :y) :size 1.5 :color "darkgrey")
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :yg1) :data F2 :size 1 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :yg2) :data F2 :size 1 :color "red" :lty 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "Accident count")
         (ggplot2/ggtitle "Figure 2.18: Poisson P-spline fits")
         (ggplot2/theme_light)))

(defn coal-smooth-plot-2
  [F2 F3]
  (rr/r+ (ggplot2/ggplot)
         (ggplot2/geom_point :data F3 (ggplot2/aes :x :knots :y :pcoef1) :size 1.5 :pch 25 :color "blue")
         (ggplot2/geom_point :data F3 (ggplot2/aes :x :knots :y :pcoef2) :size 1.5 :pch 24 :color "red")
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :eta1) :data F2 :size 1 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :eta2) :data F2 :size 1 :color "red" :lty 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "")
         (ggplot2/ggtitle "Figure 2.19: Linear predictor")
         (ggplot2/theme_light)))

(defn coal-smooth-plots
  [data]
  (let [{:keys [F2 F3]} (coal-smooth-datasets data)]
    [(coal-smooth-plot-1 data F2)
     (coal-smooth-plot-2 F2 F3)]))

(kind/table
 [(map gg/->image (coal-smooth-plots coal-histogram-data))])

;; #### Binomial Smoothing

;; Source file: `f-kypho-smooth.R`

(rr/data 'kyphosis 'rpart)

(def kyphosis-data (-> (rr/r->clj kyphosis)
                     (tc/map-columns :y :Kyphosis (fn [c] (if (= :present c) 1 0)))
                     (tc/select-columns [:Age :y])
                     (tc/rename-columns {:Age :x})))

(defn kypho-smooth-datasets
  [data]
  (let [[xl xr] (stats/extent (:x data))
        bbase (reg/b-spline-transformer xl xr 20 3)
        fit1 (fit-binomial data bbase 1)
        fit2 (fit-binomial data bbase 100)
        n (count (:beta fit1))
        knots (v/shift (v/mult (v/shift (range 1 (inc n)) -2) (/ (- xr xl) 20)) xl)
        grid (m/slice-range xl xr 100)
        p1 (map #(fit1 % true) grid)
        p2 (map #(fit2 % true) grid)]
    {:F2 (tc/dataset {:xg grid
                      :yg1 (map :fit p1) :yg2 (map :fit p2)
                      :eta1 (map :link p1) :eta2 (map :link p2)})
     :F3 (tc/dataset {:knots knots :pcoef1 (:beta fit1) :pcoef2 (:beta fit2)})}))

(defn kypho-smooth-plot-1
  [F1 F2]
  (rr/r+ (ggplot2/ggplot)
         (ggplot2/geom_point :data F1 (ggplot2/aes :x :x :y :y) :size 1.5 :color "darkgrey")
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :yg1) :data F2 :size 1 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :yg2) :data F2 :size 1 :color "red" :lty 1)
         (ggplot2/xlab "Age (months)")
         (ggplot2/ylab "P(presence)")
         (ggplot2/ggtitle "Figure 2.19: Binomial P-spline fits")
         (ggplot2/theme_light)))

(defn kypho-smooth-plot-2
  [F2 F3]
  (rr/r+ (ggplot2/ggplot)
         (ggplot2/geom_point :data F3 (ggplot2/aes :x :knots :y :pcoef1) :size 1.5 :pch 25 :color "blue")
         (ggplot2/geom_point :data F3 (ggplot2/aes :x :knots :y :pcoef2) :size 1.5 :pch 24 :color "red")
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :eta1) :data F2 :size 1 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :x :xg :y :eta2) :data F2 :size 1 :color "red" :lty 1)
         (ggplot2/xlab "Age (months)")
         (ggplot2/ylab "logit(p)")
         (ggplot2/ggtitle "Figure 2.19: Linear predictor")
         (ggplot2/theme_light)))

(defn kypho-smooth-plots
  [data]
  (let [{:keys [F2 F3]} (kypho-smooth-datasets data)]
    [(kypho-smooth-plot-1 data F2)
     (kypho-smooth-plot-2 F2 F3)]))

(kind/table
 [(map gg/->image (kypho-smooth-plots kyphosis-data))])

;; ## Optimal Smoothing in Action

;; ### Cross-Validation

;; Source file: `f-cv-plot.R`

(defn f-cv-plot-dataset
  [xg data tr lambda]
  (let [fit (fit-gaussian data tr lambda)]
    {:F2 (tc/dataset {:xg xg :yg (map fit xg)})
     :title (str "λ=" lambda
                 " | CV=" (m/approx (:cv fit) 1)
                 " | ED=" (m/round (:effective-dimension fit)))}))

(defn f-cv-plot-plots
  [data tr]
  (let [xg (m/slice-range (:x data) 100)
        lambdas (map (fn [ll] (m/pow 10.0 ll)) (range -2 3))]
    (for [lambda lambdas
          :let [{:keys [F2 title]} (f-cv-plot-dataset xg data tr lambda)]]
      (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
             (ggplot2/geom_point :size 1.5 :color "darkgrey")
             (ggplot2/geom_line (ggplot2/aes :x :xg :y :yg) :data F2 :size 1 :color "blue" :lty 1)
             (ggplot2/xlab "Time (ms)")
             (ggplot2/ylab "Acceleration (g)")
             (ggplot2/ggtitle (str "Figure 3.1: " title))
             (ggplot2/theme_light)))))

(defn f-cv-plot-cv-dataset
  [data tr]
  (let [llambs (range -4 2.01 0.1)
        lambin (map (partial m/pow 10.0) llambs)
        cvs (pmap (comp :cv (partial fit-gaussian data tr)) lambin)]
    (tc/dataset {:x llambs :y cvs})))

(defn f-cv-plot-cv-plot
  [data tr]
  (let [F3 (f-cv-plot-cv-dataset data tr)]
    (rr/r+ (ggplot2/ggplot :data F3 (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_line :size 1 :color "blue" :lty 1)
           (ggplot2/geom_point :size 1 :color "blue")
           (ggplot2/xlab "log10(λ)")
           (ggplot2/ylab "LOOCV")
           (ggplot2/ggtitle "Cross-validation profile")
           (ggplot2/theme_light))))

(defn f-cv-plot-all-plots
  [data]
  (let [tr (reg/b-spline-transformer (:x data) 50 3)
        plots (vec (f-cv-plot-plots data tr))]
    (conj plots (f-cv-plot-cv-plot data tr))))

(kind/table
 (->> (f-cv-plot-all-plots mcycle2-data)
      (map gg/->image)
      (partition 2)))

;; ### Akaike’s Information Criterion

;; Source file: `f-kyphopt.R`

(defn kyphopt-dataset
  [data tr]
  (let [llamb (range -2 3.1 0.25)]
    (-> (for [l llamb
              :let [lambda (m/pow 10.0 l)
                    aic2 (-> (fit-binomial data tr lambda 2) :ll :aic-dev)
                    aic3 (-> (fit-binomial data tr lambda 3) :ll :aic-dev)]]
          {:x l :lambda lambda :y2 aic2 :y3 aic3})
        (tc/dataset))))

(defn kyphopt-aic-plot
  [F1 selector title]
  (rr/r+ (ggplot2/ggplot F1 (ggplot2/aes :x :x :y selector))
         (ggplot2/geom_point :size 1.5 :color "blue")
         (ggplot2/xlab "log10(λ)")
         (ggplot2/ylab "AIC")
         (ggplot2/ggtitle (str  "Figure 3.2: " title))
         (ggplot2/theme_light)))

(defn kyphopt-lambda-dataset
  [data tr F1]
  (let [lambda2 (-> (tc/order-by F1 :y2)
                    (tc/rows :as-maps)
                    (first)
                    (:lambda))
        lambda3 (-> (tc/order-by F1 :y3)
                    (tc/rows :as-maps)
                    (first)
                    (:lambda))
        fit2 (fit-binomial data tr lambda2 2)
        fit3 (fit-binomial data tr lambda3 3)
        xg (m/slice-range (:x data) 100)]
    (-> (for [x xg] {:x x :y2 (fit2 x) :y3 (fit3 x)})
        (tc/dataset))))

(defn kyphopt-lambda-plot
  [data F3 selector title]
  (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point :size 1.5 :color "blue")
         (ggplot2/geom_line :data F3 :color "red" :size 1 (ggplot2/aes :x :x :y selector))
         (ggplot2/xlab "Age (months)")
         (ggplot2/ylab "(P(Kyphosis)")
         (ggplot2/ggtitle (str  "Figure 3.2: " title))
         (ggplot2/theme_light)))

(defn kyphopt-all-plots
  [data]
  (let [tr (reg/b-spline-transformer (:x data) 10 3)
        F1 (kyphopt-dataset data tr)
        F3 (kyphopt-lambda-dataset data tr F1)]
    [(kyphopt-aic-plot F1 :y2 "Second-order penalty")
     (kyphopt-aic-plot F1 :y3 "Third-order penalty")
     (kyphopt-lambda-plot data F3 :y2 "Second-order penalty")
     (kyphopt-lambda-plot data F3 :y3 "Third-order penalty")]))

(kind/table
 (->> (map gg/->image (kyphopt-all-plots kyphosis-data))
      (partition 2)))

;; ### Density estimation

;; Source file: `f-geyseropt.R`

(rr/data 'faithful 'datasets)

(def faithful-data (rr/r->clj faithful))

(defn geyseropt-histogram
  [data bw]
  (let [breaks (range 1.0 (+ 6.0 (* 0.1 bw)) bw)
        hist (:bins-maps (stats/histogram (:eruptions data) breaks))]
    (-> (tc/dataset (map #(select-keys % [:mid :count]) hist))
        (tc/rename-columns {:mid :x :count :y}))))

(defn geyseropt-aic-datasets
  [data xg lambdas bw]
  (let [hist (geyseropt-histogram data bw)
        tr (reg/b-spline-transformer (:x hist) 20 3)
        aics (vec (pmap (comp :aic-dev :ll (partial fit-poisson hist tr)) lambdas))
        ka (apply min-key aics (range (count aics)))
        fit (fit-poisson hist tr (lambdas ka))]
    {:hist hist
     :F1 (tc/dataset {:x (v/log10 lambdas) :y aics})
     :Fit (tc/dataset {:x xg :y (map fit xg)})}))

(defn geyseropt-aic-plot
  [F1 bw]
  (rr/r+ (ggplot2/ggplot F1 (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point)
         (ggplot2/ggtitle (str "Figure 3.3: Bin width " bw))
         (ggplot2/xlab "log10(lambda)")
         (ggplot2/ylab "AIC")
         (ggplot2/theme_light)))

(defn geyseropt-hist-plot
  [hist Fit]
  (rr/r+ (ggplot2/ggplot hist (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_bar :stat "identity" :fill "wheat3")
         (ggplot2/geom_hline :yintercept 0)
         (ggplot2/xlab "Eruption length (min)")
         (ggplot2/ylab "Frequency")
         (ggplot2/ggtitle "Figure 3.3: Old Faithful geyser")
         (ggplot2/geom_line :data Fit :color "steelblue" :size 1)
         (ggplot2/theme_light)))

(defn geyseropt-plots
  [data bw]
  (let [xg (m/slice-range 1.0 6.0 100)
        lambdas (mapv (partial m/pow 10.0) (range -3.0 0.01 0.1))
        {:keys [hist F1 Fit]} (geyseropt-aic-datasets data xg lambdas bw)]
    [(geyseropt-aic-plot F1 bw)
     (geyseropt-hist-plot hist Fit)]))


(kind/table
 [(map gg/->image (geyseropt-plots faithful-data 0.02))
  (map gg/->image (geyseropt-plots faithful-data 0.05))])

;; ---

;; Source file: `f-suicide-opt.R`

(rr/data 'Suicide 'JOPS)

;; We need to use `R` version of histogram, since it slightly differently solves ties comparing to `fastmath` version.

(defn suicide-opt-histogram
  [xl xr]
  (let [h (graphics/hist Suicide :breaks (base/seq xl xr :by 10) :plot false)]
    (tc/dataset {:x (rr/r->clj (rr/r$ h 'mids))
                 :y (rr/r->clj (rr/r$ h 'counts))})))

(defn suicide-opt-aic-datasets
  [xg lambdas xl xr]
  (let [hist (suicide-opt-histogram xl xr)
        tr (reg/b-spline-transformer (:x hist) 20 3)
        aics (vec (pmap (comp :aic-dev :ll (partial fit-poisson hist tr)) lambdas))
        ka (apply min-key aics (range (count aics)))
        fit (fit-poisson hist tr (lambdas ka))]
    {:hist hist
     :F1 (tc/dataset {:x (v/log10 lambdas) :y aics})
     :Fit (tc/dataset {:x xg :y (map fit xg)})}))

(defn suicide-opt-aic-plot
  [F1 title]
  (rr/r+ (ggplot2/ggplot F1 (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point)
         (ggplot2/ggtitle (str "Figure 3.4: " title))
         (ggplot2/xlab "log10(lambda)")
         (ggplot2/ylab "AIC")
         (ggplot2/theme_light)))

(defn suicide-opt-hist-plot
  [hist Fit]
  (rr/r+ (ggplot2/ggplot hist (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_bar :stat "identity" :fill "wheat3")
         (ggplot2/geom_hline :yintercept 0)
         (ggplot2/xlab "Treatment length (months)")
         (ggplot2/ylab "Frequency")
         (ggplot2/ggtitle "Figure 3.4: Suicides")
         (ggplot2/geom_line :data Fit :color "steelblue" :size 1)
         (ggplot2/theme_light)))

(defn suicide-opt-plots
  [title xl xr]
  (let [xg (m/slice-range xl xr 100)
        lambdas (mapv (partial m/pow 10.0) (range -3 4.01 0.2))
        {:keys [hist F1 Fit]} (suicide-opt-aic-datasets xg lambdas xl xr)]
    [(suicide-opt-aic-plot F1 title)
     (suicide-opt-hist-plot hist Fit)]))

(kind/table
 [(map gg/->image (suicide-opt-plots "Boundary at zero" 0 800))
  (map gg/->image (suicide-opt-plots "Boundary at -100" -100 800))])

;; ### Mixed Models

;; Only HFS algorithm presented

(defn mot-lme-mix-dataset
  [data order]
  (let [[xl xr] (stats/extent (:x data))
        tr (reg/b-spline-transformer xl xr 20 3)
        fit (:fit (HFS data fit-gaussian tr {:order order}))
        xg (m/slice-range xl xr 200)]
    (tc/dataset {:x xg :y (map fit xg)})))

(defn mot-lme-mix-plot
  ([data] (mot-lme-mix-plot data 2))
  ([data order]
   (let [F2 (mot-lme-mix-dataset mcycle2-data order)]
     (rr/r+ (ggplot2/ggplot :mapping (ggplot2/aes :x :x :y :y))
            (ggplot2/geom_hline :yintercept 0 :size 0.3)
            (ggplot2/geom_point :data data :size 1.5 :color "darkgrey")
            (ggplot2/geom_line :data F2 :size 1 :color "blue")
            (ggplot2/xlab "Time (ms)")
            (ggplot2/ylab "Acceleration (g)")
            (ggplot2/ylim [-150 100])
            (ggplot2/ggtitle "Motorcycle data fit with HFS algorithm")
            (ggplot2/labs :color "")
            (ggplot2/guides :linetype false)
            (ggplot2/theme_light)))))

(gg/->image (mot-lme-mix-plot mcycle2-data))

;; ---

;; Source file: `f-HFS-convergence.R`

;; To get similar result we have to use seeded RNG.

(def HFS-convergence-rng (r/rng :mersenne 8))

(def HFS-convergence-data (tc/dataset {:x (m/slice-range 100)
                                     :r (repeatedly 100 #(r/grandom HFS-convergence-rng))}))

(defn HFS-convergence-dataset
  [data tr phi]
  (let [sphi (str phi)
        y (v/add (v/sin (v/mult (:x data) 10.0))
                 (v/mult (:r data) phi))
        res (-> (HFS (tc/add-column data :y y) fit-gaussian tr)
                :lambdas v/log10 v/differences v/abs v/log10)        ]
    (tc/dataset {:it (range 21)
                 :dla res
                 :phi sphi})))

(defn HFS-convergence-phis-dataset
  [data phis]
  (let [tr (reg/b-spline-transformer (:x data) 20 3)]
    (->> (map (partial HFS-convergence-dataset data tr) phis)
         (reduce tc/concat))))

(defn HFS-convergence-plot
  [data phis]
  (let [L (HFS-convergence-phis-dataset data phis)]
    (rr/r+ (ggplot2/ggplot L (ggplot2/aes :x :it :y :dla :group :phi :color :phi))
           (ggplot2/geom_point)
           (ggplot2/xlab "Iteration")
           (ggplot2/ylab "Δlog10(λ)")
           (ggplot2/ggtitle "Figure 3.6: Convergence of the HFS algorithm")
           (ggplot2/labs :color "ϕ")
           (ggplot2/theme_light))))

(gg/->image (HFS-convergence-plot HFS-convergence-data [0.2 0.4 0.8]))

;; ---

;; Source file: `f-regula-falsi.R`

(defn regula-falsi-u
  [data tr llambda]
  (let [fit (fit-gaussian data tr (m/pow 10.0 llambda))
        ed (:effective-dimension fit)
        sig2 (get-in fit [:dispersions :mean-deviance])
        tau2 (/ (v/sum (v/sq (v/differences (:beta fit) 2))) (- ed 2))]
    (- (+ (m/log10 tau2) llambda)
       (m/log10 sig2))))

(defn regula-falsi-loop
  [data tr ll1 ll2]
  (loop [i 0
         [u1 u2] [(regula-falsi-u data tr ll1) (regula-falsi-u data tr ll2)]
         [ll1 ll2] [ll1 ll2]
         llambdas []]
    (if (or (> i 20)
            (< (abs (- ll2 ll1)) 1.0e-15))
      (v/log10 (v/abs (v/differences llambdas)))
      (let [dla (/ (* (- u1) (- ll2 ll1)) (- u2 u1))
            ll (+ ll1 dla)
            u (regula-falsi-u data tr ll)]
        (if (pos? (* u u1))
          (recur (inc i) [u u2] [ll ll2] (conj llambdas ll))
          (recur (inc i) [u1 u] [ll1 ll] (conj llambdas ll)))))))

(defn regula-falsi-dataset
  [data tr phi]
  (let [sphi (str phi)
        y (v/add (v/sin (v/mult (:x data) 10.0))
                 (v/mult (:r data) phi))
        res (regula-falsi-loop (tc/add-column data :y y) tr -2 3)]
    (tc/dataset {:it (range 21)
                 :dla res
                 :phi sphi})))

(defn regula-falsi-phis-dataset
  [data phis]
  (let [tr (reg/b-spline-transformer (:x data) 20 3)]
    (->> (map (partial regula-falsi-dataset data tr) phis)
         (reduce tc/concat)
         (tc/drop-missing))))

(defn regula-falsi-plot
  [data phis]
  (let [L (regula-falsi-phis-dataset data phis)]
    (rr/r+ (ggplot2/ggplot L (ggplot2/aes :x :it :y :dla :group :phi :color :phi))
           (ggplot2/geom_point)
           (ggplot2/xlab "Iteration")
           (ggplot2/ylab "Δlog10(λ)")
           (ggplot2/xlim [0 21])
           (ggplot2/ggtitle "Figure 3.7: Fast convergence of the regula falsi algorithm")
           (ggplot2/labs :color "ϕ")
           (ggplot2/theme_light))))

(gg/->image (regula-falsi-plot HFS-convergence-data [0.2 0.4 0.8]))

;; ### Bayesian P-splines

;; skipping

;; ### Dangers of Automatic Smoothing

;; Source file: `f-heyser-with-DP.R`

(rr/data 'geyser 'MASS)

(def geyser-with-DP-data (let [u (v/round (v/mult (rr/r->clj (rr/r$ geyser 'duration)) 60))]
                         (-> (->> (stats/histogram u (v/shift (range 40 351 5) 0.5))
                                  (:bins-maps)
                                  (map #(select-keys % [:mid :count]))
                                  (tc/dataset))
                             (tc/rename-columns {:mid :x :count :y}))))

(defn geyser-with-DP-data-aics
  [data tr weights llambdas]
  (map (fn [llambda]
         (let [lambda (m/pow 10.0 llambda)
               fit (fit-poisson data tr weights lambda 2)]
           [(get-in fit [:ll :aic-dev]) lambda])) llambdas))


(defn geyser-with-DP-dataset
  [data weights llambdas]
  (let [[xl xr] (stats/extent (:x data))
        tr (reg/b-spline-transformer xl xr 20 3)
        lambda (->> (geyser-with-DP-data-aics data tr weights llambdas)
                    (sort-by first)
                    (first)
                    (second))
        fit (fit-poisson data tr weights lambda 2)
        xs (m/slice-range xl xr 100)]
    (tc/dataset {:x xs :y (map fit xs)})))

(defn geyser-with-DP-plot
  ([data llambdas] (geyser-with-DP-plot data nil llambdas))
  ([data weights llambdas]
   (let [fit (geyser-with-DP-dataset data weights llambdas)]
     (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
            (ggplot2/geom_bar :stat "identity" :fill "wheat3")
            (ggplot2/geom_hline :yintercept 0)
            (ggplot2/geom_line :data fit :color "steelblue" :size 1)
            (ggplot2/xlab "Eruption length (seconds)")
            (ggplot2/ylab "Frequency")
            (ggplot2/ggtitle "Figure 3.13: Digit preference in Old Faithful geyser data")
            (ggplot2/theme_light)))))

(gg/->image (geyser-with-DP-plot geyser-with-DP-data (range -3 1.01 0.1)))

;; ---

;; Source file: `f-geyser-weights.R`

;; Set weights to `0.0` for outliers

(defn geyser-weights-plot
  [data llambdas]
  (let [weights (map #(if (> % 30) 0.0 1.0) (:y data))]
    (geyser-with-DP-plot data weights llambdas)))

(gg/->image (geyser-weights-plot geyser-with-DP-data (range -3 1.01 0.1)))

;; ---

;; Source file: `f-Greece-AIC.R`

(rr/data 'Greece_deaths 'JOPS)

(def greece-deaths-data (-> (rr/r->clj 'Greece_deaths)
                          (tc/select-rows (fn [row] (and (>= (:Age row) 40)
                                                        (< (:Age row) 85))))
                          (tc/rename-columns {:Age :x :Males :y})
                          (tc/drop-columns :Females)))


#_(kind/table greece-deaths-data {:use-datatables true})

(defn greece-AIC-optimal-lambda
  [pairs]
  (reduce (fn [t1 t2] (if (< (last t1) (last t2)) t1 t2 )) pairs))

(defn greece-AIC-data-for-nseg
  [data tr llambdas]
  (let [aics (for [llambda llambdas
                   :let [fit (fit-poisson data tr (m/pow 10.0 llambda) 2)]]
               [llambda (get-in fit [:ll :aic-dev])])
        [optimal-lambda] (greece-AIC-optimal-lambda aics)
        fit (fit-poisson data tr (m/pow 10.0 optimal-lambda) 2)]
    {:aics (map second aics) :mu (:fitted fit)}))

(defn greece-AIC-datasets
  [data llambdas]
  (let [nseg1 (greece-AIC-data-for-nseg data (reg/b-spline-transformer (:x data) 10 3) llambdas)
        nseg2 (greece-AIC-data-for-nseg data (reg/b-spline-transformer (:x data) 40 3) llambdas)]
    {:F1 (tc/dataset {:llambdas llambdas :aic1 (:aics nseg1) :aic2 (:aics nseg2)})
     :F2 (tc/add-columns data {:mu1 (:mu nseg1) :mu2 (:mu nseg2)})}))

(defn greece-AIC-plot-aic
  [F1 selector title]
  (rr/r+ (ggplot2/ggplot F1)
         (ggplot2/geom_point (ggplot2/aes :x :llambdas :y selector) :size 0.7 :color "blue")
         (ggplot2/ggtitle (str "Figure 3.15: " title " basis segments"))
         (ggplot2/xlab "log10(lambda)")
         (ggplot2/ylab "AIC")
         (ggplot2/theme_light)))

(defn greece-AIC-plot-mu
  [F2 selector title]
  (rr/r+ (ggplot2/ggplot F2)
         (ggplot2/geom_bar (ggplot2/aes :x :x :y :y) :stat "identity" :width 0.5 :fill "wheat3")
         (ggplot2/geom_line (ggplot2/aes :x :x :y selector) :color "blue")
         (ggplot2/ggtitle (str "Figure 3.15: " title " basis segments"))
         (ggplot2/xlab "Age")
         (ggplot2/ylab "Counts")
         (ggplot2/theme_light)))

(defn greece-AIC-plots
  [data]
  (let [{:keys [F1 F2]} (greece-AIC-datasets data (range -4 0.01 0.1))]
    [(greece-AIC-plot-aic F1 :aic1 "10")
     (greece-AIC-plot-mu F2 :mu1 "10")
     (greece-AIC-plot-aic F1 :aic2 "40")
     (greece-AIC-plot-mu F2 :mu2 "40")]))

(kind/table
 (->> (greece-AIC-plots greece-deaths-data)
      (map gg/->image)
      (partition 2)))

;; ---

;; Source file: `f-Greece-QL.R`

(defn greece-QL-dataset
  [data]
  (let [mu1 (-> (HFS data fit-poisson (reg/b-spline-transformer (:x data) 10 3)) :fit :fitted)
        mu2 (-> (HFS data fit-poisson (reg/b-spline-transformer (:x data) 40 3)) :fit :fitted)]
    (tc/add-columns data {:mu1 mu1 :mu2 mu2})))

(defn greece-QL-plot
  [data]
  (let [F2 (greece-QL-dataset data)]
    (rr/r+ (ggplot2/ggplot F2)
           (ggplot2/geom_bar (ggplot2/aes :x :x :y :y) :stat "identity" :fill "wheat3" :width 0.5)
           (ggplot2/geom_line (ggplot2/aes :x :x :y :mu2) :color "pink" :size 1.5)
           (ggplot2/geom_line (ggplot2/aes :x :x :y :mu1) :color "blue" :size 1.5 :linetype 2)
           (ggplot2/xlab "Age")
           (ggplot2/ylab "Counts")
           (ggplot2/ggtitle "Figure 3.16: 10 and 40 basis segments")
           (ggplot2/theme_light))))

(gg/->image (greece-QL-plot greece-deaths-data))

;; ### L- and V-curves

;; Source file: `f-wood-surf.R`

(rr/data 'Woodsurf 'spam)

(def woodsurf-data (let [d (rr/r->clj Woodsurf)]
                   (tc/add-column d :x (range 1 (inc (count (:y d)))))))

(defn wood-surf-dataset
  [data]
  (let [m (tc/row-count data)
        E (mat/eye m true)
        D (mat/differences E 2)
        P (mat/muls (mat/tmulm D D) 2800000)
        z (mat/solve (mat/add E P) (v/vec->RealVector (:y data)))]
    (tc/add-column data :z (v/vec->array z))))

(defn wood-surf-plot
  [data]
  (let [F1 (wood-surf-dataset woodsurf-data)]
    (rr/r+ (ggplot2/ggplot F1)
           (ggplot2/geom_point (ggplot2/aes :x :x :y :y) :color "darkgrey")
           (ggplot2/geom_line (ggplot2/aes :x :x :y :z) :color "blue" :size 1 :lty 1)
           (ggplot2/geom_line (ggplot2/aes :x :x :y :y) :color "darkgrey")
           (ggplot2/xlab "Position")
           (ggplot2/ylab "Height (unknown units)")
           (ggplot2/ggtitle "Figure 3.17: Wood surface")
           (ggplot2/theme_light))))

(gg/->image (wood-surf-plot woodsurf-data))

;; ---

;; Source file: `f-L-and-V-wood.R`

(defn l-and-v-wood-matrices
  [data]
  (let [m (tc/row-count data)
        E (mat/eye m true)
        D (mat/differences E 2)]
    [E (mat/tmulm D D)]))

(defn l-and-v-wood-fits-pens
  [E D2 y lambdas]
  (let [vy (v/vec->RealVector y)]
    (for [lambda lambdas
          :let [z (mat/solve (mat/add E (mat/muls D2 lambda)) vy)
                fit (m/log10 (v/sum (v/sq (v/sub vy z))))
                pen (m/log10 (v/sum (v/sq (v/differences z 2))))]]
      [fit pen])))

(defn l-and-v-wood-datasets
  [data lambdas]
  (let [[E D2] (l-and-v-wood-matrices data)
        [fits pens] (apply map vector (l-and-v-wood-fits-pens E D2 (:y data) lambdas))
        dfits (v/differences fits)
        dpens (v/differences pens)
        v (vec (v/sqrt (v/add (v/sq dfits) (v/sq dpens))))
        mlambdas (v/sqrt (v/emult (rest lambdas) (butlast lambdas)))
        k (apply min-key v (range (count v)))
        lambda (nth mlambdas k)
        [fitopt penopt] (first (l-and-v-wood-fits-pens E D2 (:y data) [lambda]))]
    {:F1 (tc/dataset {:x fits :y pens})
     :F2 (tc/dataset {:x (v/log10 mlambdas) :y v})
     :lamopt (m/log10 lambda)
     :fitopt fitopt
     :penopt penopt
     :nla (count lambdas)}))

(defn l-and-v-wood-L-plot
  [{:keys [F1 fitopt penopt nla]}]
  (rr/r+ (ggplot2/ggplot F1 (ggplot2/aes :x :x :y :y))
         (ggplot2/ggtitle "Figure 3.18: L-curve")
         (ggplot2/geom_vline :xintercept fitopt :lty 2 :color "darkgrey")
         (ggplot2/geom_hline :yintercept penopt :lty 2 :color "darkgrey")
         (ggplot2/xlab "ϕ")
         (ggplot2/ylab "ψ")
         (ggplot2/geom_point :size 0.3 :color (grDevices/rainbow nla))
         (ggplot2/xlim [0 8])
         (ggplot2/ylim [-6 6])
         (ggplot2/theme_light)))

(defn l-and-v-wood-V-plot
  [{:keys [F2 lamopt nla]}]
  (rr/r+ (ggplot2/ggplot F2 (ggplot2/aes :x :x :y :y))
         (ggplot2/ggtitle "Figure 3.18: V-curve")
         (ggplot2/geom_vline :xintercept lamopt :lty 2 :color "darkgrey")
         (ggplot2/xlab "log10(distance)")
         (ggplot2/ylab "v")
         (ggplot2/geom_point :size 0.5 :color (grDevices/rainbow (dec nla)))
         (ggplot2/theme_light)))

(defn l-and-v-wood-plots
  [data]
  (let [d (l-and-v-wood-datasets data (map (partial m/pow 10.0) (range 0 8.01 0.1)))]
    [(l-and-v-wood-L-plot d)
     (l-and-v-wood-V-plot d)]))

(kind/table
 [(map gg/->image (l-and-v-wood-plots woodsurf-data))])

;; ### Transformation of the Independent Variable

;; Source file: `f-BMI-lin.R`

(def boys7482-data-bmi (-> (tc/dataset {:x (rr/r->clj (rr/r$ boys7482 'age))
                                      :y (rr/r->clj (rr/r$ boys7482 'bmi))})
                         (tc/drop-missing)
                         (tc/select-rows (fn [row] (< 10 (:y row) 30)))))

(defn bmi-dataset
  [data init-lambda]
  (let [[xl xr] (stats/extent (:x data))
        tr (reg/b-spline-transformer xl xr 50 3)
        xg (range xl (+ xr m/EPSILON) 0.01)
        fit (:fit (HFS data fit-gaussian tr {:init-lambda init-lambda :eps 1.0e-3 :verbose? true}))]
    (tc/dataset {:x xg :y (map fit xg)})))

(defn bmi-plot
  [data init-lambda]
  (let [Z (bmi-dataset data init-lambda)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :color "darkgrey" :size 0.4)
           (ggplot2/geom_line :data Z :size 1 :color "blue")
           (ggplot2/xlab "Age")
           (ggplot2/ylab "BMI")
           (ggplot2/ggtitle "Figure 3.19")
           (ggplot2/theme_light))))

(gg/->image (bmi-plot boys7482-data-bmi 10.0))

;; ---

;; Source file: `f-BMI-sqrt.R`

(gg/->image (bmi-plot (tc/sqrt boys7482-data-bmi :x :x) 100.0))

;; ## Multidimensional Smoothing

;; ### Generalized Additive Models

;; Source file: `f-ethanol-gam.R`

(rr/data 'ethanol 'JOPS)

(def ethanol-data (rr/r->clj ethanol))

(defn ethanol-gam-datasets
  [data]
  (let [b1 (reg/b-spline-transformer 7.5 18 10 3)
        b2 (reg/b-spline-transformer 0.5 1.25 10 3)
        tr (fn [[^double x1 ^double x2]] (concat (b1 [x1]) (b2 [x2])))
        n1 (count (b1 [7.5]))
        n2 (count (b2 [0.5]))
        D1 (mat/differences (mat/eye n1) 2)
        D2 (mat/differences (mat/eye n2) 2)
        P (mat/add (mat/block-diagonal [(mat/tmulm D1 D1)
                                        (mat/muls (mat/tmulm D2 D2) 0.1)])
                   (mat/diagonal (repeat (+ n1 n2) 1.0e-6)))
        fit (reg/pglm (:NOx data) (-> data
                                      (tc/select-columns [:C :E])
                                      (tc/rows))
                      P {:transformer tr})
        diff (map < (:NOx data) (:fitted fit))
        cgrid (m/slice-range 7.5 18 50)
        egrid (m/slice-range 0.5 1.25 50)]
    {:dens (tc/dataset (for [x cgrid y egrid] {:x x :y y :NOx (fit [x y])}))
     :Fc (tc/dataset {:cgrid cgrid :zc (map (fn [x] (fit [x ##-Inf])) cgrid)})
     :Fe (tc/dataset {:egrid egrid :ze (map (fn [x] (fit [##-Inf x])) egrid)})
     :Fmod (-> data
               (tc/add-column :f1 (map (fn [x] (fit [x ##-Inf])) (:C data)))
               (tc/add-column :f2 (map (fn [x] (fit [##-Inf x])) (:E data)))
               (tc/add-column :cols (map #(if %1 "blue" "yellow") diff))
               (tc/add-column :pchs (map #(if %1 "+" "-") diff)))}))

(defn ethanol-gam-plot1
  [{:keys [Fmod]}]
  (rr/r+ (ggplot2/ggplot Fmod (ggplot2/aes :x :C :y :E))
         (ggplot2/geom_point :color "darkgrey")
         (ggplot2/xlab "Compression ratio (C)")
         (ggplot2/ylab "Equivalence ratio (E)")
         (ggplot2/ggtitle "Figure 4.1: Experiment design")
         (ggplot2/theme_light)))

(defn ethanol-gam-plot2
  [{:keys [Fmod Fe]}]
  (rr/r+ (ggplot2/ggplot Fmod (ggplot2/aes :x :E :y '(- NOx f1)))
         (ggplot2/geom_point :color "darkgrey")
         (ggplot2/geom_line (ggplot2/aes :x :egrid :y :ze) :data Fe :size 1 :color "blue")
         (ggplot2/xlab "Equivalence ratio")
         (ggplot2/ylab "Partial residuals")
         (ggplot2/ggtitle "Figure 4.1: Partial response")
         (ggplot2/theme_light)))

(defn ethanol-gam-plot3
  [{:keys [Fmod Fc]}]
  (rr/r+ (ggplot2/ggplot Fmod (ggplot2/aes :x :C :y '(- NOx f2)))
         (ggplot2/geom_point :color "darkgrey")
         (ggplot2/geom_line (ggplot2/aes :x :cgrid :y :zc) :data Fc :size 1 :color "blue")
         (ggplot2/xlab "Compression ratio")
         (ggplot2/ylab "Partial residuals")
         (ggplot2/ggtitle "Figure 4.1: Partial response")
         (ggplot2/theme_light)))

(defn ethanol-gam-plot4
  [{:keys [Fmod]}]
  (rr/r+ (ggplot2/ggplot Fmod (ggplot2/aes :x '(+ f1 f2) :y :NOx))
         (ggplot2/geom_point :color "darkgrey")
         (ggplot2/geom_abline :slope 1 :intercept 0 :color "blue" :size 1)
         (ggplot2/xlab "Fitted NOx")
         (ggplot2/ylab "Observed")
         (ggplot2/ggtitle "Figure 4.1: Compare fit to data")
         (ggplot2/theme_light)))

(defn ethanol-gam-plots
  [data]
  (let [dss (ethanol-gam-datasets data)]
    [(ethanol-gam-plot1 dss)
     (ethanol-gam-plot2 dss)
     (ethanol-gam-plot3 dss)
     (ethanol-gam-plot4 dss)]))

(kind/table
 (->> (ethanol-gam-plots ethanol-data)
      (map gg/->image)
      (partition 2)))

;; ---

;; Source file: `f-ethanol-gam-surf.R`

(defn ethanol-gam-surf-plot
  [data]
  (let [{:keys [dens Fmod]} (ethanol-gam-datasets data)]
    (rr/r+ (ggplot2/ggplot dens (ggplot2/aes :x :x :y :y :fill :NOx :z :NOx))
           (ggplot2/geom_raster :show.legend true)
           (ggplot2/scale_fill_gradientn :colours (grDevices/terrain-colors 100))
           (ggplot2/geom_contour :data dens :color "blue" :show.legend true)
           (metR/geom_text_contour :color "blue" :size 3)
           (ggplot2/xlab "Compression ratio")
           (ggplot2/ylab "Equivalence ratio")
           (ggplot2/ggtitle "Figure 4.2: GAM for NOx emission, ethanol data")
           (ggplot2/geom_point :data Fmod (ggplot2/aes :x :C :y :E) :shape (:pchs Fmod) :size 5)
           (ggplot2/theme_light)
           (ggplot2/theme :panel.grid.major (ggplot2/element_blank)
                          :panel.grid.minor (ggplot2/element_blank)))))

(gg/->image (ethanol-gam-surf-plot ethanol-data))

;; ### Varying Coefficient Models

;; Source file: `f-vcm4up.R`

(rr/data 'Disks 'JOPS)

(def disks-data (-> (rr/r->clj Disks)
                  (tc/map-rows (fn [{:keys [Year Month PriceDG]}]
                                 {:Month (if (= Year 2000) (+ Month 12) Month)
                                  :Price (/ PriceDG 2.2)}))))

(def disks-data-groups (-> disks-data
                         (tc/group-by :Month)
                         (tc/groups->map)))

(def vcm4up-months (->> (Locale. "en")
                      (.getDisplayNames (Calendar/getInstance) Calendar/MONTH Calendar/LONG)
                      (into {})
                      (map (comp vec reverse))
                      (into {})))

(defn vcm4up-slopes
  [[month ds]]
  (let [xg (m/slice-range (:Size ds) 100)
        fit (reg/lm (:Price ds) (:Size ds))]
    [month {:slope (-> fit :beta first (m/approx 1))
            :ds ds
            :month month
            :month-name (vcm4up-months (mod (dec month) 12))
            :fit (tc/dataset {:x xg :y (map fit xg)})}]))

(defn vcm4up-plot
  [{:keys [slope fit ds month-name]}]
  (rr/r+ (ggplot2/ggplot ds (ggplot2/aes :x :Size :y :Price))
         (ggplot2/ylim [100 500])
         (ggplot2/xlim [5 35])
         (ggplot2/geom_point :color "grey40")
         (ggplot2/geom_line :data fit :mapping (ggplot2/aes :x :x :y :y) :color "blue")
         (ggplot2/xlab "Hard drive size (GB)")
         (ggplot2/ylab "Price (Euro)")
         (ggplot2/ggtitle (str "Figure 4.3: " month-name " slope = " slope))
         (ggplot2/theme_light)))

(defn vvm4up-plots
  [data months]
  (let [dss (into {} (map vcm4up-slopes data))]
    (map (comp vcm4up-plot dss) months)))

(kind/table
 (->> (vvm4up-plots disks-data-groups [2 5 9 12])
      (map gg/->image)
      (partition 2)))

;; ---

;; Source file: `f-vcmsmooth.R`

(defn vcmsmooth-vara
  [xtxinv v0 v]
  (let [vv (v/vec->RealVector (concat v0 (v/vec->seq v)))]
    (v/dot vv (mat/mulv xtxinv vv))))

(defn vcmsmooth-dataset
  [data]
  (let [tr-month (reg/b-spline-transformer (:Month data) 10 3)
        B (bbase-matrix tr-month (:Month data))
        n (count (tr-month [1]))
        D (mat/differences (mat/eye n) 2)
        P (mat/muls (mat/tmulm D D) 100.0)
        C1 (mat/mulm (mat/diagonal (:Size data)) B)
        C (mat/bind-cols B C1)
        Q (mat/add (mat/block-diagonal P P)
                   (mat/muls (mat/eye (mat/ncol C)) 1.0e-6))
        fit (reg/pglm (:Price data) (mat/rows C) Q)
        a2 (drop n (:beta fit)) ;; only slope coeffs
        
        xg (m/slice-range 2 13 100)
        Bg (bbase-matrix tr-month xg)
        trendgrid (v/vec->seq (mat/mulv Bg (v/vec->RealVector a2))) ;; slope only

        vara (map (partial vcmsmooth-vara (:xtxinv fit) (repeat n 0.0)) (mat/rows Bg)) ;; skip intercept
        sea2 (v/mult (v/sqrt vara) (* 2.0 (m/sqrt (:dispersion fit))))]
    (tc/dataset {:xg xg :trendgrid trendgrid
                 :se2up (v/add trendgrid sea2) :se2lo (v/sub trendgrid sea2)})))

(defn vcmsmooth-slope
  [[month ds]]
  (let [fit (reg/lm (:Price ds) (:Size ds))]
    {:slope (-> fit :beta first) :month month}))


(defn vcmsmooth-plot
  [data groups]
  (let [F2dat (vcmsmooth-dataset data)
        slopeframe (tc/dataset (map vcmsmooth-slope groups))]
    (rr/r+ (ggplot2/ggplot slopeframe (ggplot2/aes :x :month :y :slope))
           (ggplot2/geom_point :size 2.2)
           (ggplot2/scale_x_continuous :breaks (range 2 14)
                                       :labels ["Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                                "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"])
           (ggplot2/geom_line (ggplot2/aes :x :xg :y :trendgrid) :data F2dat :size 1 :color "blue")
           (ggplot2/geom_line (ggplot2/aes :x :xg :y :se2up) :data F2dat :size 1 :lty 2 :color "red")
           (ggplot2/geom_line (ggplot2/aes :x :xg :y :se2lo) :data F2dat :size 1 :lty 2 :color "red")
           (ggplot2/geom_hline :yintercept 0 :color "darkgrey")
           (ggplot2/ylab "Slope: Euro/GB")
           (ggplot2/xlab "Month (1999-2000)")
           (ggplot2/ylim [0 30])
           (ggplot2/theme_light))))

(gg/->image (vcmsmooth-plot disks-data disks-data-groups))

;; ---

;; Source file: `f-polio-vcm.R`

(rr/data 'polio 'gamlss.data)

(def polio-data (let [y (rr/r->clj (base/as-numeric polio))
                    x (->> (range 1 (inc (count y)))
                           (map (fn [m] (+ 1970 (/ (- m 0.5) 12)))))]
                (tc/dataset {:y y :x x })))

(defn polio-vcm-penalty-Q
  [K P lambdas]
  (-> (map mat/muls (repeat P) lambdas)
      (mat/block-diagonal)
      (mat/add K)))

(defn polio-vcm-model
  [data]
  (let [tr (reg/b-spline-transformer (:x data) 10 3)
        B (bbase-matrix tr (:x data))
        n (mat/ncol B)
        D (mat/differences (mat/eye n) 2)
        P (mat/tmulm D D)
        om (v/mult (seq (:x data)) m/TWO_PI)
        cs (v/cos om)
        sn (v/sin om)
        cs2 (v/cos (v/mult om 2))
        sn2 (v/sin (v/mult om 2))
        C (mat/bind-cols B
                         (mat/mulm (mat/diagonal cs) B)
                         (mat/mulm (mat/diagonal sn) B)
                         (mat/mulm (mat/diagonal cs2) B)
                         (mat/mulm (mat/diagonal sn2) B))
        K (mat/diagonal (repeat (* 5 n) 1.0e-10))]
    {:P P :C C :n n :K K :cs cs :sn sn :cs2 cs2 :sn2 sn2 :B B}))

;; We need effective dimension per each term

(defn pre-hat
  [xtxinv xss weights]
  (let [sqrt-weights (map m/sqrt weights)
        whx (mat/mulm (mat/diagonal sqrt-weights) xss)]
    (->> (mat/rows whx)
         (map (fn [row] (map v/vec->seq [row (mat/mulv xtxinv row)]))))))

(defn polio-vcm-fit
  [data]
  (let [{:keys [P C n K] :as model} (polio-vcm-model data)
        design-m (mat/rows C)]
    (loop [i 0
           lambdas [1 1 1 1 1]
           dla ##Inf
           fit nil]
      (if (or (= i 20) (< dla 1.0e-10))
        (assoc model
               :as (partition n (:beta fit))
               :mu (:fitted fit))
        (let [nfit (reg/pglm (:y data) design-m (polio-vcm-penalty-Q K P lambdas) {:family :poisson})
              a (->> (partition n (:beta nfit))
                     (map (comp v/sum v/sq)))
              h (->> (pre-hat (:xtxinv nfit) C (get-in nfit [:weights :weights]))
                     (map (fn [[l r]]
                            (map v/dot (partition n l) (partition n r))))
                     (apply map vector)
                     (map v/sum))
              nlambdas (v/ediv h a)]
          (recur (inc i) nlambdas (v/mx (v/abs (v/sub (v/log10 lambdas) (v/log10 nlambdas)))) nfit))))))

(defn polio-vcm-dataset
  [data]
  (let [{:keys [B as cs sn cs2 sn2 mu]} (polio-vcm-fit data)
        [fit1 fit2 fit3 fit4 fit5] (map (comp v/vec->seq (partial mat/mulv B) v/vec->RealVector) as)]
    (tc/add-columns data {:mu mu :trend fit1
                          :seas1 (v/add (v/emult fit2 cs) (v/emult fit3 sn))
                          :seas2 (v/add (v/emult fit4 cs2) (v/emult fit5 sn2))
                          :amp1 (v/sqrt (v/add (v/sq fit2) (v/sq fit3)))
                          :amp2 (v/sqrt (v/add (v/sq fit4) (v/sq fit5)))})))

(defn polio-vcm-plot1
  [ds]
  (rr/r+ (ggplot2/ggplot ds (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point :size 0.8)
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/geom_line (ggplot2/aes :y :mu) :size 0.7 :color "blue" :lty 1)
         (ggplot2/xlab "")
         (ggplot2/ylab "Counts per month")
         (ggplot2/ggtitle "Figure 4.5: USA polio cases")
         (ggplot2/ylim [0 15])
         (ggplot2/theme_light)))

(defn polio-vcm-plot2
  [ds]
  (rr/r+ (ggplot2/ggplot ds (ggplot2/aes :x :x :y :trend))
         (ggplot2/geom_line :color "green4" :size 0.5 :lty 1)
         (ggplot2/geom_hline :yintercept 0 :size 0.3)
         (ggplot2/geom_line (ggplot2/aes :y :seas1) :size 0.5 :color "blue" :lty 1) 
         (ggplot2/geom_line (ggplot2/aes :y :seas2) :size 0.5 :color "red" :lty 1) 
         (ggplot2/geom_line (ggplot2/aes :y :amp1) :size 0.4 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :y '(- amp1)) :size 0.4 :color "blue" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :y :amp2) :size 0.4 :color "red" :lty 2)
         (ggplot2/geom_line (ggplot2/aes :y '(- amp2)) :size 0.4 :color "red" :lty 2)
         (ggplot2/xlab "Year")
         (ggplot2/ylab "")
         (ggplot2/ggtitle "Figure 4.5: Linear predictors")
         (ggplot2/theme_light)))


(defn polio-vcm-plots
  [data]
  (let [ds (polio-vcm-dataset data)]
    [(polio-vcm-plot1 ds)
     (polio-vcm-plot2 ds)]))

(kind/table
 (map (comp vector #(gg/->image % {:height 400 :width 800})) (polio-vcm-plots polio-data)))

;; ### Tensor Product Models

;; skipped

;; ### Two-Dimensional Penalties

;; Source file: `f-ethanol-psp2.R`

(defn ethanol-psp2-base
  [data]
  (let [trnx (reg/b-spline-transformer 7 19 10 3)
        trny (reg/b-spline-transformer 0.5 1.25 10 3)
        Bx (bbase-matrix trnx (:C data))
        By (bbase-matrix trny (:E data))
        nx (mat/ncol Bx)
        ny (mat/ncol By)
        B1 (mat/kronecker (mat/cols->RealMatrix (repeat ny [1])) Bx)
        B2 (mat/kronecker By (mat/cols->RealMatrix (repeat nx [1])))
        B (mat/emulm B1 B2)]
    {:B B :nx nx :ny ny :trnx trnx :trny trny}))

(defn ethanol-psp2-penalty
  [{:keys [nx ny] :as base} l1 l2]
  (let [Ix (mat/eye nx true)
        Iy (mat/eye ny true)
        Dx (mat/differences Ix 2)
        Dy (mat/differences Iy 2)
        Px (mat/kronecker Iy (mat/add (mat/tmulm Dx Dx) (mat/muls Ix 1.0e-10)))
        Py (mat/kronecker (mat/add (mat/tmulm Dy Dy) (mat/muls Iy 1.0e-10)) Ix)]
    (assoc base :D (mat/add (mat/muls Px l1) (mat/muls Py l2)))))

(defn ethanol-psp2-dataset
  [data l1 l2]
  (let [{:keys [B D trnx trny]} (-> data ethanol-psp2-base (ethanol-psp2-penalty l1 l2))
        u (m/slice-range 7 19 50)
        v (m/slice-range 0.5 1.25 50)
        fit (reg/pglm (:NOx data) (mat/rows B) D)
        a (:beta fit)
        diff (map < (:NOx data) (:fitted fit))]
    {:dens (tc/dataset (for [x u y v
                             :let [bx (trnx [x])
                                   by (trny [y])]]
                         {:x x :y y :NOx (v/dot a (box-product by bx))}))
     :Fmod (tc/add-column data :pchs (map #(if %1 "+" "-") diff))}))

(defn ethanol-psp2-surf-plot
  [data l1 l2]
  (let [{:keys [dens Fmod]} (ethanol-psp2-dataset data l1 l2)]
    (rr/r+ (ggplot2/ggplot dens (ggplot2/aes :x :x :y :y :fill :NOx :z :NOx))
           (ggplot2/geom_raster :show.legend true)
           (ggplot2/scale_fill_gradientn :colours (grDevices/terrain-colors 100))
           (ggplot2/geom_contour :data dens :color "blue" :show.legend true)
           (metR/geom_text_contour :color "blue" :size 3)
           (ggplot2/xlab "Compression ratio")
           (ggplot2/ylab "Equivalence ratio")
           (ggplot2/ggtitle "Figure 4.8: 2D P-splines for NOx emission, ethanol data")
           (ggplot2/geom_point :data Fmod (ggplot2/aes :x :C :y :E) :shape (:pchs Fmod) :size 5 :color "blue")
           (ggplot2/theme_light)
           (ggplot2/theme :panel.grid.major (ggplot2/element_blank)
                          :panel.grid.minor (ggplot2/element_blank)))))

(gg/->image (ethanol-psp2-surf-plot ethanol-data 1.0 0.1))

;; ### Interpolation and Extrapolation

;; Source file: `f-psp-corner.R`

(def rng2 (r/rng :mersenne 1))

(def psp-corner-data (let [[x y] (split-at 200 (repeatedly 400 #(r/drandom rng2 -1 1)))
                         z (v/add (v/sq (v/exp (v/sub (v/sub (v/sq x)) (v/sq y))))
                                  (repeatedly 200 #(r/grandom rng2 0.1)))]
                     (tc/dataset {:x x :y y :z z})))

(def psp-corner-data2 (tc/select-rows psp-corner-data (fn [{:keys [x y]}] (or (pos? x) (pos? y)))))

(defn psp-corner-base
  [data]
  (let [trnx (reg/b-spline-transformer -1 1 10 3)
        trny (reg/b-spline-transformer -1 1 20 3)
        Bx (bbase-matrix trnx (:x data))
        By (bbase-matrix trny (:y data))
        nx (mat/ncol Bx)
        ny (mat/ncol By)
        B1 (mat/kronecker (mat/cols->RealMatrix (repeat ny [1])) Bx)
        B2 (mat/kronecker By (mat/cols->RealMatrix (repeat nx [1])))
        B (mat/emulm B1 B2)]
    {:B B :nx nx :ny ny :trnx trnx :trny trny}))

(defn psp-corner-penalty
  [{:keys [nx ny] :as base}]
  (let [Ix (mat/eye nx true)
        Iy (mat/eye ny true)
        Dx (mat/differences Ix 2)
        Dy (mat/differences Iy 2)
        Px (mat/kronecker Iy (mat/add (mat/tmulm Dx Dx) (mat/muls Ix 1.0e-10)))
        Py (mat/kronecker (mat/add (mat/tmulm Dy Dy) (mat/muls Iy 1.0e-10)) Ix)]
    (assoc base :D (mat/add Px Py))))

(defn psp-corner-dataset
  [data]
  (let [{:keys [B D trnx trny]} (-> data psp-corner-base psp-corner-penalty)
        u (m/slice-range -1 1 50)
        v (m/slice-range -1 1 50)
        fit (reg/pglm (:z data) (mat/rows B) D)
        a (:beta fit)
        diff (map < (:z data) (:fitted fit))]
    {:dens (tc/dataset (for [x u y v
                             :let [bx (trnx [x])
                                   by (trny [y])]]
                         {:x x :y y :z (v/dot a (box-product by bx))}))
     :Fmod (tc/add-column data :pchs (map #(if %1 "+" "-") diff))}))

(defn psp-corner-surf-plot
  [data]
  (let [{:keys [dens Fmod]} (psp-corner-dataset data)]
    (rr/r+ (ggplot2/ggplot dens (ggplot2/aes :x :x :y :y :fill :z :z :z))
           (ggplot2/geom_raster :show.legend true)
           (ggplot2/scale_fill_gradientn :colours (grDevices/terrain-colors 100))
           (ggplot2/geom_contour :data dens :color "blue" :show.legend true)
           (metR/geom_text_contour :color "blue" :size 3)
           (ggplot2/xlab "")
           (ggplot2/ylab "")
           (ggplot2/ggtitle "Figure 4.9")
           (ggplot2/geom_point :data Fmod (ggplot2/aes :x :x :y :y) :shape (:pchs Fmod) :size 5 :color "blue")
           (ggplot2/theme_light)
           (ggplot2/theme :panel.grid.major (ggplot2/element_blank)
                          :panel.grid.minor (ggplot2/element_blank)))))

(kind/table
 [[(gg/->image (psp-corner-surf-plot psp-corner-data))
   (gg/->image (psp-corner-surf-plot psp-corner-data2))]])

;; ### Generalized Two-Dimensional Smoothing

;; Source file: `f-geyser-aniso.R`

(defn hist2d
  ([xs ys] (hist2d xs ys 100))
  ([xs ys b] (hist2d xs ys (stats/extent xs false) (stats/extent ys false) b))
  ([xs ys [xl xr] [yl yr] b]
   (let [xf (m/make-norm xl xr 0.0 1.0)
         yf (m/make-norm yl yr 0.0 1.0)
         xfr (m/make-norm 0.0 1.0 xl xr)
         yfr (m/make-norm 0.0 1.0 yl yr)
         s (/ 1.0 b)
         g (grid/grid :square s)
         r (range b)
         cells (into {} (for [x r y r
                              :let [[mx my] (grid/cell->mid g x y)]]
                          [[(xfr mx) (yfr my)] 0]))
         fq (->> (for [[x y] (map vector xs ys)
                       :let [[mx my] (grid/coords->mid g (xf x) (yf y))]]
                   [(xfr mx) (yfr my)])
                 (reduce (fn [c k] (update c k inc)) cells)
                 #_(frequencies))
         kfq (keys fq)]
     (tc/dataset {:x (map first kfq) :y (map second kfq) :z (vals fq)}))))

(def faithful-hist2d-data (hist2d (faithful-data :eruptions) (faithful-data :waiting) [1 6] [40 100] 100))

(defn geyser-base
  [data]
  (let [trnx (reg/b-spline-transformer 1 6 10 3)
        trny (reg/b-spline-transformer 40 100 10 3)
        Bx (bbase-matrix trnx (:x data))
        By (bbase-matrix trny (:y data))
        nx (mat/ncol Bx)
        ny (mat/ncol By)
        B1 (mat/kronecker (mat/cols->RealMatrix (repeat ny [1])) Bx)
        B2 (mat/kronecker By (mat/cols->RealMatrix (repeat nx [1])))
        B (mat/emulm B1 B2)]
    {:B B :nx nx :ny ny :trnx trnx :trny trny}))

(defn geyser-penalty
  [nx ny l1 l2]
  (let [Ix (mat/eye nx true)
        Iy (mat/eye ny true)
        Dx (mat/differences Ix 2)
        Dy (mat/differences Iy 2)
        Px (mat/kronecker Iy (mat/add (mat/tmulm Dx Dx) (mat/muls Ix 1.0e-4)))
        Py (mat/kronecker (mat/add (mat/tmulm Dy Dy) (mat/muls Iy 1.0e-4)) Ix)]
    (mat/add (mat/muls Px l1) (mat/muls Py l2))))

;; This is too slow, need to refactor to optimize fit only.

(defn geyser-aniso-optimal
  [z xy nx ny]
  (fn [^double l1 ^double l2]
    (let [pl1 (m/pow 10.0 l1)
          pl2 (m/pow 10.0 l2)
          D (geyser-penalty nx ny pl1 pl2)
          fit (reg/pglm z xy D {:family :poisson})]
      (get-in fit [:ll :aic-dev]))))

(defn geyser-dataset
  [data l1 l2]
  (let [{:keys [B nx ny trnx trny]} (geyser-base data)
        ;; slow!
        #_#_optl1l2 (opt/minimize :powell (geyser-aniso-optimal (:z data) (mat/rows B) nx ny)
                                  {:bounds [[-5 5] [-5 5]]
                                   :initial [0 0]})
        #_ (println optl1l2 (map (partial m/pow 10) optl1l2))
        D (geyser-penalty nx ny l1 l2)
        u (m/slice-range 1 6 100)
        v (m/slice-range 40 100 100)
        fit (reg/pglm (:z data) (mat/rows B) D {:family :poisson :epsilon 1.0e-4})
        a (:beta fit)]
    (tc/dataset (for [x u y v
                      :let [bx (trnx [x])
                            by (trny [y])]]
                  {:x x :y y :z (v/sqrt (v/exp (v/dot a (box-product by bx))))}))))

(defn geyser-plot
  [data title l1 l2]
  (let [dens (geyser-dataset data l1 l2)]
    (rr/r+ (ggplot2/ggplot dens (ggplot2/aes :x :x :y :y :fill :z :z :z))
           (ggplot2/geom_raster :show.legend false)
           (ggplot2/scale_fill_gradient :high "darkgreen" :low "white")
           (ggplot2/ggtitle title)
           (ggplot2/xlab "Eruption length (min)")
           (ggplot2/ylab "Waiting time (min)")
           (ggplot2/geom_contour :color "steelblue" :show.legend true)
           (ggplot2/theme_light))))

;; No optimization (too slow), just used inferred data from R

(kind/table
 [[(gg/->image (geyser-plot faithful-hist2d-data "Figure 4.12: Old Faithful, isotropic smoothing (square root of density)"
                            (m/pow 10 -2.249608) (m/pow 10 -0.5806777)))
   (let [l (m/pow 10 -1.476649)]
     (gg/->image (geyser-plot faithful-hist2d-data "Figure 4.13: Old Faithful, anisotropic smoothing (square root of density)" l l)))]])

;; ### Nested Bases and PS-ANOVA

;; skipped

;; ## Smoothing of Scale and Shape

;; ### Quantile Smoothing

;; Only iterative algorithm (third image)

;; Source file: `mot-median-iter.R`

(defn mot-median-dataset
  [data]
  (let [[xlo xhi] (stats/extent (:x data) false)
        tr (reg/b-spline-transformer xlo xhi 50 3)
        B (bbase-matrix tr (:x data))
        xg (m/slice-range xlo xhi 1000)
        Bg (bbase-matrix tr xg)
        n (mat/ncol B)
        D (mat/differences (mat/eye n true) 2)
        P (mat/tmulm D D)
        fit (reg/pglm (:y data) (mat/rows B) P {:mode {:method :quantile :tau 0.5}})]
    (tc/dataset {:x xg :y (map fit (mat/rows Bg))})))

(defn mot-median-plot
  [data]
  (let [ds (mot-median-dataset data)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :size 1)
           (ggplot2/geom_hline :yintercept 0 :size 0.3)
           (ggplot2/geom_line :data ds :size 0.5 :color "red" :lty 1)
           (ggplot2/geom_line :data ds :size 1.5 :color "red" :lty 2)
           (ggplot2/xlab "Time (ms)")
           (ggplot2/ylab "Acceleration (g)")
           (ggplot2/ggtitle "Figure 5.3: Motorcycle helmet impact data")
           (ggplot2/ylim [-150 100])
           (ggplot2/theme_light))))

(gg/->image (mot-median-plot mcycle2-data))

;; ---

;; Source file: `f-quant-curves-1000.R`

(defn boys7482-data-wgt-fn [size] (let [ds (-> (rr/r->clj boys7482)
                                            (tc/select-columns [:age :wgt])
                                            (tc/select-rows (fn [{:keys [age]}] (> age 5)))
                                            (tc/drop-missing)
                                            (tc/rename-columns {:age :x :wgt :y}))
                                     n (tc/row-count ds)]
                                 (tc/select-rows ds (repeatedly size #(r/irand n)))))

(def boys7482-data-wgt-100 (boys7482-data-wgt-fn 100))
(def boys7482-data-wgt-1000 (boys7482-data-wgt-fn 1000))

(defn quant-curves-fit
  [ys B P xg Bg pp]
  (let [fit (reg/pglm ys (mat/rows B) P {:mode {:method :quantile :tau pp}})]
    (tc/dataset {:x xg :y (map fit (mat/rows Bg)) :tau (str (m/approx pp))})))

(defn quant-curves-dataset
  [data pps]
  (let [[xl xr] (stats/extent (:x data) false)
        tr (reg/b-spline-transformer xl xr 20 3)
        B (bbase-matrix tr (:x data))
        xg (m/slice-range xl xr 100)
        Bg (bbase-matrix tr xg)
        D (mat/differences (mat/eye (mat/ncol B) true) 2)
        P (mat/tmulm D D)]
    (->> pps
         (map (partial quant-curves-fit (:y data) B P xg Bg))
         (reduce tc/concat))))

(defn quant-curves-plot
  [data title pps]
  (let [ds (quant-curves-dataset data pps)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :color "grey75" :size 0.9)
           (ggplot2/geom_line :data ds (ggplot2/aes :x :x :y :y :group :tau :color :tau) :size 1)
           (ggplot2/labs :color '(bquote tau))
           (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (count pps) :start 10 :end 350))
           (ggplot2/xlab "Age (yr)")
           (ggplot2/ylab "Weight (kg)")
           (ggplot2/ggtitle title)
           (ggplot2/theme_light))))


(gg/->image (quant-curves-plot boys7482-data-wgt-1000
                               "Figure 5.4: Quantile curves for weights of 1000 Dutch boys"
                               (range 0.05 0.951 0.1)))

;; ---

;; Source file: `f-quant-curves-100.R`

(gg/->image (quant-curves-plot boys7482-data-wgt-100
                               "Figure 5.5: Quantile curves for weights of 100 Dutch boys"
                               [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95]))

;; ---

;; Source file: `f-quant-sheet-100.R`

;; skipped

;; ### Expectile Smoothing

(defn exp-curves-fit
  [ys B P xg Bg pp]
  (let [fit (reg/pglm ys (mat/rows B) P {:mode {:method :expectile :tau pp}})]
    (tc/dataset {:x xg :y (map fit (mat/rows Bg)) :tau (str (m/approx pp 5))})))

(defn exp-curves-dataset
  [data pps lambda]
  (let [tr (reg/b-spline-transformer 5 22 20 3)
        B (bbase-matrix tr (:x data))
        xg (m/slice-range 5 22 100)
        Bg (bbase-matrix tr xg)
        D (mat/differences (mat/eye (mat/ncol B) true) 2)
        P (mat/muls (mat/tmulm D D) lambda)]
    (->> pps
         (map (partial exp-curves-fit (:y data) B P xg Bg))
         (reduce tc/concat))))

(defn exp-curves-plot
  [data title pps lambda]
  (let [ds (exp-curves-dataset data pps lambda)]
    (rr/r+ (ggplot2/ggplot data (ggplot2/aes :x :x :y :y))
           (ggplot2/geom_point :color "grey75" :size 0.9)
           (ggplot2/geom_line :data ds (ggplot2/aes :x :x :y :y :group :tau :color :tau) :size 1)
           (ggplot2/labs :color '(bquote tau))
           (ggplot2/scale_color_manual :values (colorspace/rainbow_hcl (count pps) :start 10 :end 350))
           (ggplot2/xlab "Age (yr)")
           (ggplot2/ylab "Weight (kg)")
           (ggplot2/ggtitle title)
           (ggplot2/theme_light))))

(kind/table
 [[(gg/->image (exp-curves-plot boys7482-data-wgt-100
                                "Figure 5.7: Expectile curves for weights of 100 Dutch boys (λ=1)"
                                [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95] 1.0))
   (gg/->image (exp-curves-plot boys7482-data-wgt-100
                                "Figure 5.7: Expectile curves for weights of 100 Dutch boys (λ=10)"
                                [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95] 10.0))]])

;; ---

;; Source file: `f-exp-curves-1000.R`

(gg/->image (exp-curves-plot boys7482-data-wgt-1000
                             "Figure 5.9: Expectile curves for weights of 1000 Dutch boys (λ=10)"
                             [0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.5, 0.8, 0.9, 0.97, 0.98, 0.99, 0.997, 0.999] 10.0))

;; ---

;; Source file: `f-exp-curves-1000-dens.R`

(defn exp-curves-1000-dens-fit-15
  [ys B P x15 pp]
  (let [fit (reg/pglm ys (mat/rows B) P {:mode {:method :expectile :tau pp}})]
    (fit x15)))

(defn exp-curves-1000-dens-ze
  [data pps]
  (let [tr (reg/b-spline-transformer 5 22 20 3)
        B (bbase-matrix tr (:x data))
        D (mat/differences (mat/eye (mat/ncol B) true) 2)
        P (mat/muls (mat/tmulm D D) 10.0)]
    (->> pps
         (map (partial exp-curves-1000-dens-fit-15 (:y data) B P (tr [15]))))))

(exp-curves-1000-dens-ze boys7482-data-wgt-1000 [0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.5, 0.8, 0.9, 0.97, 0.98, 0.99, 0.997, 0.999])

;; skipped the rest

;; ### Baseline Estimation

;; Source file: `f-baseline.R`

(rr/data 'indiumoxide 'JOPS)

(def indiumoxide-data (-> (rr/r->clj indiumoxide)
                        (tc/rename-columns {1 :x 2 :y})
                        (tc/select-rows (fn [{:keys [x]}] (m/< 47 x 65)))))

(defn baseline-dataset
  [data]
  (let [tr (reg/b-spline-transformer (:x data) 40 3)
        B (bbase-matrix tr (:x data))
        D (mat/differences (mat/eye (mat/ncol B) true) 2)
        P (mat/muls (mat/tmulm D D) 1000)
        fit (reg/pglm (:y data) (mat/rows B) P {:mode {:method :expectile :tau 0.04}})]
    (-> (tc/add-column data :z (:fitted fit))
        (tc/- :r [:y :z]))))

(defn baseline-plot1
  [dat]
  (rr/r+ (ggplot2/ggplot dat (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_line :color "darkgrey" :size 0.2)
         (ggplot2/geom_line (ggplot2/aes :y :z) :color "blue" :size 0.6)
         (ggplot2/ggtitle "Figure 5.16: X-ray diffraction scan of indium oxide with baseline")
         (ggplot2/xlab "Diffraction angle")
         (ggplot2/ylab "Counts")
         (ggplot2/theme_light)
         (ggplot2/ylim [-25 750])))

(defn baseline-plot2
  [dat]
  (rr/r+ (ggplot2/ggplot dat (ggplot2/aes :x :x :y :r))
         (ggplot2/geom_line :color "darkgrey" :size 0.2)
         (ggplot2/ggtitle "Figure 5.16: Baseline subtracted")
         (ggplot2/xlab "Diffraction angle")
         (ggplot2/ylab "Counts")
         (ggplot2/theme_light)
         (ggplot2/ylim [-25 750])))

(defn baseline-plots
  [data]
  (let [dat (baseline-dataset data)]
    [(baseline-plot1 dat)
     (baseline-plot2 dat)]))

(kind/table
 [(map gg/->image (baseline-plots indiumoxide-data))])

;; ## Complex Counts and Composite Links

;; skipped

;; ## Signal Regression

;; Source file: `f-psrcoef.R`

(rr/data 'nirc 'fds)

(def nirc-psr-data (let [X (-> (rr/r$ nirc 'y)
                             (rr/r->clj)
                             (tc/select-rows (range 49 650))
                             (tc/drop-columns [:$row.names])
                             (tc/rows)
                             (mat/rows->RealMatrix))
                       iindex (-> (rr/r$ nirc 'x)
                                  (rr/r->clj)
                                  (subvec 49 650))
                       diindex (rest iindex)
                       y (-> (rr/bra fds/labc 1 (range 1 41))
                             (rr/r->clj)
                             (->> (keep-indexed #(when (not= %1 22) %2))))
                       dX (-> (mat/differences X)
                              (mat/cols)
                              (->> (keep-indexed #(when (not= %1 22) %2)))
                              (mat/rows->RealMatrix))]
                   {:y y :dX dX :diindex diindex}))


;; skipped the rest, possible return in future

;; ## Special Subjects

;; ### Harmonic Smoothing

;; Source file: `f-varstar.R`

(rr/data 'Varstar 'JOPS)

(def varstar-data (let [ds (-> (tc/dataset {:x (rr/r->clj (rr/r$ Varstar 'V1))
                                          :y (rr/r->clj (rr/r$ Varstar 'V2))})
                             (tc/select-rows (fn [{:keys [x y]}] (and (< y 16)
                                                                     (or (< x 3000) (> x 3300))))))
                      m (stats/mean (:y ds))]
                  (tc/map-columns ds :y (fn [x] (- x m)))))

(defn varstar-dataset
  [data]
  (let [xl 2400
        xr 4200
        nseg 100
        tr (reg/b-spline-transformer xl xr 100 3)
        B (bbase-matrix tr (:x data))
        n (mat/ncol B)
        D (diffs-matrix n 2)
        D2 (mat/muls (mat/tmulm D D) 0.1)
        phi (/ (* m/TWO_PI (/ (- xr xl) nseg)) 110.0)
        Dp (hdiffs-matrix n phi)
        Dp2 (mat/muls (mat/tmulm Dp Dp) 0.1)
        fit1 (reg/pglm (:y data) (mat/rows B) D2)
        fit2 (reg/pglm (:y data) (mat/rows B) Dp2)
        xg (m/slice-range xl xr 500)
        Bg (bbase-matrix tr xg)]
    (tc/dataset {:x xg :z (map fit1 (mat/rows Bg)) :zh (map fit2 (mat/rows Bg))})))

(defn varstar-plot
  [data ds selector title]
  (rr/r+ (ggplot2/ggplot :data data (ggplot2/aes :x :x :y :y))
         (ggplot2/geom_point :color "tomato" :size 0.9)
         (ggplot2/geom_line (ggplot2/aes :x :x :y selector) :data ds :color "blue")
         (ggplot2/theme_light)
         (ggplot2/ylim [-1.5 1.5])
         (ggplot2/ggtitle (str "Figure 8.2: " title))
         (ggplot2/xlab "Day")
         (ggplot2/ylab "Centered magnitude")))

(defn varstar-plots
  [data]
  (let [ds (varstar-dataset data)]
    [(varstar-plot data ds :z "Variable star, standard penalty, d = 2")
     (varstar-plot data ds :zh "Variable star, harmonic penalty")]))

(kind/table
 [(map gg/->image (varstar-plots varstar-data))])

;; ### Circular Smoothing

;; Source file: `f-periodic.R`

;; TODO

;; 
