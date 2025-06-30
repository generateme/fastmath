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
            [fastmath.grid :as g]))

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
              '[rpart])

;; ## Common functions

(defn b-spline-base
  "Builds B-spline base function."
  ([{:keys [xl xr nseg degree] :or {xl 0.0 xr 1.0 nseg 15 degree 3}}] (b-spline-base xl xr nseg degree))
  ([xl xr nseg degree] (reg/b-spline-transformer xl xr nseg degree)))

;; Following creates a row of coordinates of B-spline bases (5+2 knots, degree 2)

(let [bfn (b-spline-base 0 1 5 2)]
  (bfn 0.25))

(defn bbase-matrix
  "Calculates B-spline coordinates for given base function and points."
  [bsfn xs]
  (mat/rows->RealMatrix (map bsfn xs)))

(defn diffs-matrix
  "Creates differences (penalty) matrix for smoothing."
  [size degree]
  (mat/differences (mat/eye size) degree))

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

(def boys7482-data (-> (tc/dataset {:x (rr/r->clj (rr/r$ boys7482 'age))
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
 [(map gg/->image (slope-height-plots boys7482-data 50 3))])

;; ### The Effective Dimension

(defn peffdim-fit
  [data bbase lambin order]
  (map (fn [^double lambda]
         (:effective-dimension (reg/glm (:y data) (:x data)
                                        {:intercept? false
                                         :transformer bbase
                                         :augmentation :diffs                                         
                                         :augmentation-param {:lambda lambda :order order}}))) lambin))

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

(defn se-dataset
  [data]
  (let [[xl xr] (stats/extent (:x data))
        res (reg/glm (:y data) (:x data)
                     {:intercept? false
                      :transformer (reg/b-spline-transformer xl xr 20 3)
                      :augmentation :diffs
                      :augmentation-param {:lambda 0.8 :order 2}})]
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

;; ### Equivalent Kernels

;; ### Smoothing of a Non-normal Response

;; #### Poisson Smoothing

(rr/data 'coal 'boot)

(let [tr (reg/b-spline-transformer 0 4 20 2)
      fit (reg/glm [7 3 0 1 1] [0 1 2 3 4]
                   {:intercept? false
                    :family :poisson
                    :transformer tr
                    :augmentation :diffs
                    :augmentation-param {:lambda 1 :order 2}})]
  (:ll fit))

;; #### Binomial Smoothing

(rr/data 'kyphosis 'rpart)

(def kyphosis-data (-> (rr/r->clj kyphosis)
                     (tc/map-columns :y :Kyphosis (fn [c] (if (= :present c) 1 0)))))

(let [tr (reg/b-spline-transformer (:Age kyphosis-data) 20 3)
      fit1 (reg/glm (:y kyphosis-data) (:Age kyphosis-data)
                    {:intercept? false
                     :family :binomial
                     :transformer tr
                     :augmentation :diffs
                     :augmentation-param {:lambda 1 :order 2}})]
  (fit1 100 true))

(stats/extent (:Age kyphosis-data))
