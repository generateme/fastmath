(ns fastmath.utils.gen-images
  "Generate images from examples attached to metadata."
  (:require [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.transform :as t]
            [fastmath.interpolation :as i]
            [cljplot.build :as b]
            [cljplot.core :as cljplot]
            [clojure2d.color :as clr]
            [clojure2d.core :as c2d]))

(r/set-seed! r/default-rng 1234)

(def ^:const bg-color 0x30426a)
(def ^:const fg-color 0xb2bfdc)

(defn save-chart
  "Save chart under given path and name"
  ([c prefix n suff]
   (cljplot/save c (str "docs/images/" prefix "/" n suff)))
  ([c prefix n]
   (save-chart c prefix n ".jpg")))

(defn symbol->fn
  "Convert symbol to function"
  [s] (eval `(fn [x#] (~s x#))))

(defn symbol->fn2
  "Convert symbol to 2 parameter function"
  [s] (eval `(fn [x# y#] (~s x# y#))))

;; core

(defn function-chart
  ([f] (function-chart f nil))
  ([f d]
   (cljplot/xy-chart {:width 400 :height 200 :background bg-color}
                     (b/series [:vline 0 {:color [60 100 120]}]
                               [:hline 0 {:color [60 100 120]}]
                               [:function f {:domain (or d [-3.1 3.1])
                                             :color :white
                                             :samples 300}])
                     (b/update-scale :x :ticks 5)
                     (b/update-scale :y :ticks 5)
                     (b/add-axes :bottom {:ticks {:color fg-color}
                                          :line {:color fg-color}})
                     (b/add-axes :left {:ticks {:color fg-color}
                                        :line {:color fg-color}}))))

(defn function2d-chart
  [f d]
  (cljplot/xy-chart {:width 300 :height 300 :background bg-color}
                    (b/series [:function-2d f (merge d {:gradient (clr/gradient [bg-color :white])})])
                    (b/update-scale :x :ticks 5)
                    (b/update-scale :y :ticks 5)
                    (b/add-axes :bottom {:ticks {:color fg-color}
                                         :line {:color fg-color}})
                    (b/add-axes :left {:ticks {:color fg-color}
                                       :line {:color fg-color}})))


(doseq [f `(m/abs m/acos m/acosh m/acot m/acoth m/acsc m/acsch m/asec m/asech m/asin m/asinh m/atan m/atanh
                  m/cbrt m/ceil m/cos m/cosh m/cot m/coth m/csc m/csch
                  [m/digamma [0 5]]
                  m/erf m/erfc m/exp m/expm1
                  m/floor m/frac
                  [m/gamma [0 5]]
                  m/haversine [m/high-2-exp [0.01 10]]
                  m/iabs m/inv-erf m/inv-erfc [m/inv-gamma-1pm1 [-0.5 1.5]] m/itrunc
                  m/ln m/log m/log-gamma [m/log-gamma-1p [-0.5 1.5]] m/log10 m/log1p m/log1pexp m/log2 [m/low-2-exp [0.01 10]]
                  m/logit
                  m/pow3 m/pow2 
                  [m/qsqrt [0.01 5]] m/qsin [m/qlog [0.01 5]] m/qexp m/qcos
                  [m/rqsqrt [0.01 5]] m/round-up-pow2 m/round m/rint 
                  m/sin m/sqrt m/sq m/sinh m/sinc m/signum m/sigmoid m/sgn m/sfrac m/sech m/sec m/safe-sqrt 
                  m/trunc [m/trigamma [0 2]] m/tanh m/tan)]
  (let [[f d] (if (vector? f) f [f nil])]
    (save-chart (function-chart (symbol->fn f) d) "m" (name f) ".png")))

(doseq [f `([m/bessel-j {:x [0 5] :y [0 5]}]
            [m/log-beta {:x [0 5] :y [0 5]}]
            m/atan2 m/hypot m/hypot-sqrt)]
  (let [[f d] (if (vector? f) f [f nil])]
    (save-chart (function2d-chart (symbol->fn2 f) (or d {})) "m" (name f) ".png")))


(doseq [i `(m/cos-interpolation m/lerp m/wrap m/smooth-interpolation m/quad-interpolation)]
  (save-chart (function-chart (eval `(fn [x#] (~i 0.0 1.0 x#))) [0 1]) "m" (name i) ".png"))


(save-chart (function2d-chart (fn [x y] (m/erf x y)) {}) "m" "erf2" ".png")

;; random

(defn scatter-chart
  [d]
  (cljplot/xy-chart {:width 300 :height 300 :background bg-color}
                    (b/series [:scatter d {:color fg-color}])
                    (b/update-scale :x :ticks 5)
                    (b/update-scale :y :ticks 5)
                    (b/add-axes :bottom {:ticks {:color fg-color}
                                         :line {:color fg-color}})
                    (b/add-axes :left {:ticks {:color fg-color}
                                       :line {:color fg-color}})))

(doseq [nm r/sequence-generators-list]
  (save-chart (scatter-chart (take 500 (r/sequence-generator nm 2))) "r" (name nm) ".jpg"))

(doseq [nm r/sequence-generators-list]
  (save-chart (scatter-chart (take 500 (r/jittered-sequence-generator nm 2 0.5))) "r" (str "j" (name nm)) ".jpg"))

;; noise

(doseq [[nm noise-fn] [[:noise r/noise]
                       [:vnoise r/vnoise]
                       [:simplex r/simplex]
                       [:random1 (r/random-noise-fn {:seed 1234 :generator :fbm :interpolation :linear})]
                       [:random2 (r/random-noise-fn {:seed 1234 :generator :billow})]
                       [:random3 (r/random-noise-fn {:seed 1234 :generator :ridgemulti :octaves 2})]
                       [:fbm (r/fbm-noise {:seed 1234})]
                       [:billow (r/billow-noise {:seed 1234})]
                       [:ridgedmulti (r/ridgedmulti-noise {:seed 1234})]
                       [:single (r/single-noise {:seed 1234})]]]
  (save-chart (function2d-chart noise-fn {:x [-2 2] :y [-2 2]}) "n" (name nm) ".jpg"))

(save-chart (function2d-chart r/discrete-noise {:x [-10 10] :y [-10 10]}) "n" "discrete_noise" ".jpg")
(save-chart (function2d-chart (r/warp-noise-fn) {:x [-2 2] :y [-2 2]}) "n" "warp" ".jpg")

;; distributions

(def pal (reverse (clr/palette-presets :category20)))

(defn distr-chart
  ([f d pnames ps] (distr-chart f d pnames ps nil))
  ([f d pnames ps domain]
   (let [ff (case f
              :cdf r/cdf
              :pdf r/pdf
              :icdf r/icdf)
         fs (map-indexed (fn [id p]
                           [:function
                            (partial ff (r/distribution d (zipmap pnames p)))
                            {:domain (or domain [-3.1 3.1])
                             :color (nth pal id)
                             :samples 300}]) ps)]
     (cljplot/xy-chart {:width 400 :height 230 :background bg-color}
                       (-> (b/series [:vline 0 {:color [60 100 120]}]
                                     [:hline 0 {:color [60 100 120]}])
                           (b/add-series fs))
                       (b/update-scale :x :ticks 5)
                       (b/update-scale :y :ticks 5)
                       (b/add-axes :bottom {:ticks {:color fg-color}
                                            :line {:color fg-color}})
                       (b/add-axes :left {:ticks {:color fg-color}
                                          :line {:color fg-color}})
                       (b/add-label :bottom (name d) {:color fg-color})
                       (b/add-label :left (name f) {:color fg-color})))))

(defn save-distr-charts
  ([distr params ps] (save-distr-charts distr params ps nil))
  ([distr params ps domain]
   (binding [c2d/*jpeg-image-quality* 0.75]
     (save-chart (distr-chart :pdf distr params ps domain) "d" (str "pdf-" (name distr)) ".jpg")
     (save-chart (distr-chart :cdf distr params ps domain) "d" (str "cdf-" (name distr)) ".jpg")
     (save-chart (distr-chart :icdf distr params ps [0 0.99]) "d" (str "icdf-" (name distr)) ".jpg"))))

#_(defn show-distr-charts
    ([distr params ps] (show-distr-charts distr params ps nil))
    ([distr params ps domain]
     (binding [c2d/*jpeg-image-quality* 0.75]
       (cljplot/show (pdf-chart :pdf distr params ps domain))
       (cljplot/show (pdf-chart :cdf distr params ps domain))
       (cljplot/show (pdf-chart :icdf distr params ps [0 0.99])))))

(def empirical-data (sort (repeatedly 1000 #(r/grand (r/randval 0.3 -1 1) 1))))

(save-distr-charts :beta [:alpha :beta] [[5 2] [2 5] [0.5 0.5] [5 5]] [0.01 0.99])
(save-distr-charts :cauchy [:median :scale] [[0 1] [1 2] [-1 0.5]])
(save-distr-charts :chi-squared [:degrees-of-freedom] [[2] [3] [5]] [0 5])
(save-distr-charts :exponential [:mean] [[0.1] [0.5] [1] [2]] [0 1])
(save-distr-charts :f [:numerator-degrees-of-freedom :denominator-degrees-of-freedom] [[1 1] [2 2] [5 2] [100 100]] [0 5])
(save-distr-charts :gamma [:shape :scale] [[2 2] [10 0.5] [5 0.5]] [0 10])
(save-distr-charts :gumbel [:mu :beta] [[1 2] [0.5 1] [3 4]] [-5 10])
(save-distr-charts :laplace [:mu :beta] [[-1 2] [0.5 0.5] [0 1]])
(save-distr-charts :levy [:mu :c] [[-1 2] [0.5 0.5] [0 1]] [-1 5])
(save-distr-charts :logistic [:mu :s] [[-1 2] [0.5 0.5] [0 1]])
(save-distr-charts :log-normal [:scale :shape] [[1 1] [0.5 0.5] [2 2]] [0 5])
(save-distr-charts :nakagami [:mu :omega] [[1 1] [0.5 0.5] [2 2]] [0 5])
(save-distr-charts :normal [:mu :sd] [[0 1] [1 2] [-1 0.5]])
(save-distr-charts :pareto [:scale :shape] [[2 1] [0.5 0.5] [1 3]] [0 5])
(save-distr-charts :t [:degrees-of-freedom] [[1] [5] [0.5]])
(save-distr-charts :triangular [:a :c :b] [[-1 0 1] [-3 -1 3] [0.5 1 3]])
(save-distr-charts :uniform-real [:lower :upper] [[-1 1] [-3 2] [1.5 3]])
(save-distr-charts :weibull [:alpha :beta] [[2 1] [5 1] [1 1]] [0 5])
(save-distr-charts :empirical [:data :bin-count] [[empirical-data 1000]])
(save-distr-charts :enumerated-real [:data] [[empirical-data]])

(save-distr-charts :negative-binomial [:r :p] [[20 0.5] [10 0.9] [10 0.2]] [0 50])
(save-distr-charts :bernoulli [:p] [[0.5] [0.9] [0.2]] [0 2])
(save-distr-charts :enumerated-int [:data] [[(map int empirical-data)]])
(save-distr-charts :binomial [:trials :p] [[20 0.5] [10 0.9] [10 0.2]] [0 20])
(save-distr-charts :geometric [:p] [[0.5] [0.9] [0.2]] [0 7])
(save-distr-charts :hypergeometric [:population-size :number-of-successes :sample-size] [[100 50 25] [50 10 5] [50 40 30]] [0 40])
(save-distr-charts :pascal [:r :p] [[20 0.5] [10 0.9] [10 0.2]] [0 50])
(save-distr-charts :poisson [:p] [[0.5] [0.9] [0.2]] [0 7])
(save-distr-charts :uniform-int [:lower :upper] [[-1 1] [-3 2] [1.5 3]] [-4 4])
(save-distr-charts :zipf [:number-of-elements :exponent] [[100 3] [50 0.5] [10 1]] [0 7])

(save-distr-charts :anderson-darling [:n] [[1] [3]] [0 3])
(save-distr-charts :inverse-gamma [:alpha :beta] [[2 1] [1 2] [2 2]] [0 5])
(save-distr-charts :chi [:nu] [[1] [2] [3]] [0 3])
(save-distr-charts :chi-squared-noncentral [:nu :lambda] [[1 1] [1 3] [2 0.5]] [0 5])
(save-distr-charts :erlang [:k :lambda] [[1 1] [2 1] [2 2]] [0 3])
(save-distr-charts :fatigue-life [:mu :beta :gamma] [[0 1 1] [1 2 3] [0 2 2]] [0 3])
(save-distr-charts :folded-normal [:mu :sigma] [[0 1] [1 2] [1 0.5]] [0 5])
(save-distr-charts :frechet [:alpha :beta :delta] [[1 1 0] [2 1 -1] [0.5 0.5 2]] [-1 5])
(save-distr-charts :hyperbolic-secant [:mu :sigma] [[0 1] [1 2] [-1 0.5]])
(save-distr-charts :inverse-gaussian [:mu :lambda] [[2 1] [1 2] [2 2]] [0 5])
(save-distr-charts :hypoexponential-equal [:n :k :h] [[1 1 1] [2 2 2] [2 2 3]] [0 5])

(save-distr-charts :johnson-sb [:gamma :delta :xi :lambda] [[0 1 0 1] [1 1 -2 2] [-2 2 1 1]])
(save-distr-charts :johnson-sl [:gamma :delta :xi :lambda] [[0 1 0 1] [1 1 -2 2] [-2 2 1 1]] [-3 5])
(save-distr-charts :johnson-su [:gamma :delta :xi :lambda] [[0 1 0 1] [1 1 -2 2] [-2 2 1 1]] [-6 5])

(save-distr-charts :kolmogorov-smirnov [:n] [[1] [2] [3]] [0 3])
(save-distr-charts :kolmogorov-smirnov+ [:n] [[1] [2] [3]] [0 3])

(save-distr-charts :log-logistic [:alpha :beta] [[3 1] [1 3] [2 2]] [0 5])
(save-distr-charts :pearson-6 [:alpha1 :alpha2 :beta] [[1 1 1] [0.5 2 2] [3 3 0.5]] [0 5])
(save-distr-charts :power [:a :b :c] [[0 1 2] [0 2 3] [1 3 2]] [0 5])
(save-distr-charts :rayleigh [:a :beta] [[0 1] [2 0.5] [-1 2]] [-3 5])

(save-distr-charts :watson-g [:n] [[2] [40]] [0 2])
(save-distr-charts :watson-u [:n] [[2] [40]] [0 0.5])

(save-distr-charts :hypoexponential [:lambdas] [[[1.0]] [[2 3 4]] [[0.5 0.1 2 5]]] [0 5])

(save-distr-charts :reciprocal-sqrt [:a] [[0.5] [2] [3]] [0 5])

(save-distr-charts :continuous-distribution [:data] [[empirical-data]])
(save-distr-charts :real-discrete-distribution [:data] [[empirical-data]])
(save-distr-charts :integer-discrete-distribution [:data] [[empirical-data]])

(save-distr-charts :half-cauchy [:scale] [[1] [2] [0.5]] [0 5])

(save-chart (function2d-chart (fn [x y] (r/pdf (r/distribution :multi-normal) [x y])) {:x [-3.1 3.1] :y [-3.1 3.1]})
            "d" "multi-normal" ".jpg")
(save-chart (function2d-chart (fn [x y] (r/pdf (r/distribution :multi-normal {:covariances [[1 -1] [-1 2]]}) [x y])) {:x [-3.1 3.1] :y [-3.1 3.1]})
            "d" "multi-normal2" ".jpg")
(save-chart (function2d-chart (fn [x y] (r/pdf (r/distribution :dirichlet {:alpha [2 0.8]}) [x y])) {:x [0 1] :y [0 1]})
            "d" "dirichlet" ".jpg")

#_(comment (do
             (defn generate-math-graph
               "Generate math graphs for function"
               [canvas f a b]
               (let [mid-x (m/norm 0 a b 0 200)
                     mid-y (- 200 mid-x)]
                 (-> canvas
                     (set-color 60 100 120)
                     (set-stroke 0.8)
                     (line 0 mid-y 200 mid-y)
                     (line mid-x 0 mid-x 200)
                     (set-color :white 130) 
                     (set-stroke 1.5))
                 (dotimes [x 400]
                   (try
                     (let [xx (m/norm x 0 400 a b) 
                           y (m/norm (f xx) b a 0 200)]
                       (point canvas (* 0.5 x) y))
                     (catch Exception e ))))
               canvas)

             (doseq [s m/single-list]
               (let [n (name s)
                     f (symbol->fn s)
                     c (with-canvas-> (canvas 200 200)
                         (set-background bg-color)
                         (generate-math-graph f (- m/PI) m/PI))]
                 (save-canvas c "m" n ".png")))

             (doseq [s m/interp-list]
               (let [n (name s)
                     f (eval `(fn [x#] (~s 0.0 1.0 x#)))
                     c (with-canvas-> (canvas 200 200)
                         (set-background bg-color)
                         (generate-math-graph f -0.1 1.1))]
                 (save-canvas c "m" n ".png")))

             ;; complex

             (defn generate-complex-graph
               "Generate graph for complex fn."
               [canvas f]
               (let [w (width canvas)
                     h (height canvas)]
                 ;; (set-stroke canvas 1.5)
                 (set-color canvas :white 20)
                 (dotimes [x w]
                   (dotimes [y h]
                     (let [xx (m/norm x 0 w (- m/PI) m/PI)
                           yy (m/norm y 0 h (- m/PI) m/PI)
                           res (f (c/complex xx yy))
                           resx (m/norm (res 0) (- m/PI) m/PI 0 w)
                           resy (m/norm (res 1) (- m/PI) m/PI 0 h)]
                       (rect canvas resx resy 1 1)))))
               canvas)

             (defn generate-complex-graph2
               "Generate graph for complex fn."
               [canvas f]
               (let [w (width canvas)
                     h (height canvas)]
                 (dotimes [x w]
                   (dotimes [y h]
                     (let [xx (m/norm x 0 w (- m/PI) m/PI)
                           yy (m/norm y 0 h (- m/PI) m/PI)
                           res (f (c/complex xx yy))
                           arg (m/norm (c/arg res) (- m/PI) m/PI 0 255)
                           mag (m/cnorm (c/abs res) (- m/PI) m/PI 0.0 255.0)]
                       (set-color canvas (from-HSB (v/vec4 arg 255 mag 255)))
                       (rect canvas x y 1 1)))))
               canvas)


             (binding [*jpeg-image-quality* 0.9]
               (doseq [s c/fn-list]
                 (let [n (name s)
                       f (symbol->fn s)
                       can (canvas 200 200 :low)
                       c (with-canvas [c can]
                           (set-background c bg-color)
                           (generate-complex-graph2 c f)
                           (if-not (= n "identity")
                             (generate-complex-graph c f)
                             c))]
                   (save-canvas c "c" n))))

             ;; random/noise

             (def background-vec (to-color bg-color))
             (def white (to-color :white))

             (defn icolor 
               [t]
               (v/interpolate background-vec white t))

             ;; discrete noise
             (let [canvas (canvas 200 200)]
               (with-canvas [c canvas]
                 (set-background c bg-color)
                 (dotimes [x 180]
                   (dotimes [y 180]
                     (let [n (r/discrete-noise x y)]
                       (set-color c (icolor n))
                       (point c (+ x 10) (+ y 10))))))
               (save-canvas canvas "n" "discrete_noise"))

             (defn draw-and-save-noise
               "Draw noise and save to the file"
               [f nm]
               (let [canvas (canvas 200 200)]
                 (with-canvas [c canvas] 
                   (set-background c bg-color)
                   (dotimes [x 180]
                     (dotimes [y 180]
                       (let [xx (/ x 50.0)
                             yy (/ y 50.0)
                             n (f xx yy)]
                         (set-color c (icolor n))
                         (point c (+ x 10) (+ y 10))))))
                 (save-canvas canvas "n" nm)))

             ;; noise

             (doseq [s `(r/noise r/vnoise r/simplex)]
               (let [nm (name s)
                     f (symbol->fn2 s)]
                 (draw-and-save-noise f nm)))

             ;; combined noise (fbm...)

             (draw-and-save-noise (r/make-fbm-noise {:interpolation :linear
                                                     :noise-type :value}) "fbm")

             (draw-and-save-noise (r/make-ridgedmulti-noise {:octaves 3
                                                             :lacunarity 2.1
                                                             :gain 0.7
                                                             :noise-type :simplex}) "ridgedmulti")

             (draw-and-save-noise (r/make-billow-noise {:seed 12345
                                                        :interpolation :linear}) "billow")

             (draw-and-save-noise (r/make-single-noise {:interpolation :linear}) "single")

             (draw-and-save-noise (r/make-random-noise-fn) "random1")
             (draw-and-save-noise (r/make-random-noise-fn) "random2")
             (draw-and-save-noise (r/make-random-noise-fn) "random3")


             ;; random distributions
             ;; random series

             (defn draw-and-save-random-2d
               "Draw noise and save to the file"
               [f nm]
               (let [canvas (canvas 200 200)]
                 (with-canvas [c canvas]
                   (set-background c bg-color)
                   (set-color c :white 200)
                   (set-stroke c 2.0)
                   (doseq [^fastmath.vector.Vec2 res (take 1000 f)]
                     (let [x (* (.x res) 200)
                           y (* (.y res) 200)]
                       (point c x y))))
                 (save-canvas canvas "r" nm)))

             (draw-and-save-random-2d ((r/make-sequence-generator :sobol 2)) "sobol")
             (draw-and-save-random-2d ((r/make-sequence-generator :halton 2)) "halton")
             (draw-and-save-random-2d ((r/make-sequence-generator :default 2)) "default")
             (draw-and-save-random-2d (map #(v/add (v/mult % 0.1) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :gaussian 2))) "gaussian")
             (draw-and-save-random-2d (map #(v/add (v/mult % 0.4) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :sphere 2))) "sphere")

;;;; interpolations

             (defn real-fn
               ""
               [x]
               (let [a (m/norm x 0 200 0 (* 2.2 m/PI))
                     f (m/sin (* a (* 0.5 (m/cos (inc a)))))]
                 (m/norm f -1.0 1.0 20 180)))

             (def xs [20 50 58 66 100 121 140 150 160 170 180])
             (def ys (map real-fn xs))

             (def real-fn-path (for [x (range 200)]
                                 (v/vec2 x (real-fn x))))

             (defn draw-interpolation
               "Draw interpolation line"
               [interp-name interp]
               (let [canvas (canvas 200 200)
                     p (remove nil? (for [x (range 200)]
                                      (try
                                        (v/vec2 x (interp x))
                                        (catch Exception e nil))))]
                 (with-canvas [c canvas]
                   (set-background c bg-color)
                   (set-color c :cyan)
                   (doseq [[x y] (map vector xs ys)]
                     (ellipse c x y 8 8 true))
                   (set-color c 255 100 100 200)
                   (path c real-fn-path)
                   (set-color c :white 200)
                   (set-stroke c 2.0)
                   (path c p))
                 (save-canvas canvas "i" (name interp-name))))

             (draw-interpolation :shepard (i/shepard-interpolator xs ys))
             (draw-interpolation :shepard-5 (i/shepard-interpolator xs ys 5))
             (draw-interpolation :shepard-09 (i/shepard-interpolator xs ys 0.9))

             (draw-interpolation :linear-smile (i/linear-smile-interpolator xs ys))
             (draw-interpolation :kriging (i/kriging-spline-interpolator xs ys))
             (draw-interpolation :cubic-spline (i/cubic-spline-interpolator xs ys))
             (draw-interpolation :spline (i/spline-interpolator xs ys))
             (draw-interpolation :neville (i/neville-interpolator xs ys))
             (draw-interpolation :divided-diff (i/divided-difference-interpolator xs ys))
             (draw-interpolation :akima (i/akima-spline-interpolator xs ys))
             (draw-interpolation :linear (i/linear-interpolator xs ys))

             (draw-interpolation :loess (i/loess-interpolator xs ys))
             (draw-interpolation :loess-07-4 (i/loess-interpolator xs ys 0.7 4))
             (draw-interpolation :loess-02-1 (i/loess-interpolator xs ys 0.2 1))

             (draw-interpolation :rbf-gaussian (i/rbf-interpolator xs ys (i/rbf :gaussian 120)))
             (draw-interpolation :rbf-imultiquadratic (i/rbf-interpolator xs ys (i/rbf :inverse-multiquadratic 80)))
             (draw-interpolation :rbf-multiquadratic (i/rbf-interpolator xs ys (i/rbf :multiquadratic 120)))
             (draw-interpolation :rbf-thinplate (i/rbf-interpolator xs ys (i/rbf :thinplate 80)))

             (draw-interpolation :microsphere (i/microsphere-projection-interpolator xs ys 8 0.9 0.000001 0 1.5 false 1))

             ;; draw 2d interpolation

             (def xs [20 50 58 66 100 121 140 150 160 170 180])
             (def ys [20 30 58 66 90  121 140 152 170 172 180])

             (def vs (map (fn [_] (repeatedly (count xs) #(r/irand 256))) ys))

             (defn draw-interpolation-2d
               "Draw interpolation map"
               [interp-name interp]
               (let [canvas (canvas 200 200)]
                 (with-canvas [c canvas]
                   (set-background c bg-color)
                   (dotimes [x 200]
                     (dotimes [y 200]
                       (try
                         (set-color c :white (interp x y))
                         (rect c x y 1 1)
                         (catch Exception e nil))))
                   (set-color c 255 100 100 200)
                   (doseq [x xs
                           y ys]
                     (ellipse c x y 5 5 true)))
                 (save-canvas canvas "i" (name interp-name))))

             (draw-interpolation-2d :bicubic (i/bicubic-interpolator xs ys vs))
             (draw-interpolation-2d :piecewise-bicubic (i/piecewise-bicubic-interpolator xs ys vs))
             (draw-interpolation-2d :microsphere-2d (i/microsphere-2d-projection-interpolator xs ys vs 32 0.99 0.00001 1 2 true 1))
             (draw-interpolation-2d :bilinear (i/bilinear-interpolator xs ys vs))
             (draw-interpolation-2d :bicubic-smile (i/bicubic-smile-interpolator xs ys vs))
             (draw-interpolation-2d :cubic-2d (i/cubic-2d-interpolator xs ys vs))

             ;;

             (defn draw-distribution
               "Draw real distribution."
               [a b s distr-name & distr]
               (let [canvas (canvas 200 200)
                     p (for [d distr]
                         (for [x (range 1 200)]
                           (let [xx (m/norm x 0 200 a b)
                                 y (- 200 (* s (r/pdf d xx)))]
                             (v/vec2 x y))))]
                 (with-canvas [c canvas]
                   (set-background c bg-color)
                   (doseq [[p col] (map vector p [:white
                                                  (color 102 190 141)
                                                  (color 258 164 88)
                                                  (color 244 77 81)])] 
                     (set-color c col)
                     (path c p)))
                 (save-canvas canvas "d" (name distr-name))))

             (draw-distribution 0 3 200 :levy
                                (r/distribution :levy {:c 0.5})
                                (r/distribution :levy {:c 1})
                                (r/distribution :levy {:c 2}))

             (draw-distribution 0.0001 1 50 :beta
                                (r/distribution :beta {:alpha 2 :beta 5})
                                (r/distribution :beta {:alpha 0.5 :beta 0.5})
                                (r/distribution :beta {:alpha 2.5 :beta 2}))

             (draw-distribution -4 4 250 :cauchy
                                (r/distribution :cauchy {})
                                (r/distribution :cauchy {:mean 0.5 :scale 3.0})
                                (r/distribution :cauchy {:scale 0.5}))

             (draw-distribution 0 8 250 :chi-squared
                                (r/distribution :chi-squared {:degrees-of-freedom 1})
                                (r/distribution :chi-squared {:degrees-of-freedom 5})
                                (r/distribution :chi-squared {:degrees-of-freedom 9}))

             (draw-distribution -1.1 1.1 60 :empirical
                                (r/distribution :empirical {:bin-count 10
                                                            :data (repeatedly 10000 #(m/sin (r/drand 0 m/TWO_PI)))})
                                (r/distribution :empirical {:bin-count 15
                                                            :data (repeatedly 10000 #(* (r/drand -1 1) (r/drand -1 1) (r/drand -1 1)))}))

             (draw-distribution 0.0001 2.5 50 :weibull
                                (r/distribution :weibull)
                                (r/distribution :weibull {:alpha 0.5 :beta 1.0})
                                (r/distribution :weibull {:alpha 5.0 :beta 1.0}))

             (draw-distribution 0 3 100 :exponential
                                (r/distribution :exponential)
                                (r/distribution :exponential {:mean 0.5})
                                (r/distribution :exponential {:mean 1.5}))

             (draw-distribution 0 3 200 :f
                                (r/distribution :f)
                                (r/distribution :f {:denominator-degrees-of-freedom 5 :numerator-degrees-of-freedom 5})
                                (r/distribution :f {:denominator-degrees-of-freedom 5 :numerator-degrees-of-freedom 1}))

             (draw-distribution 0 20 300 :gamma
                                (r/distribution :gamma)
                                (r/distribution :gamma {:shape 1 :scale 2})
                                (r/distribution :gamma {:shape 2 :scale 0.5}))

             (draw-distribution 0 23 600 :binomial
                                (r/distribution :binomial)
                                (r/distribution :binomial {:trials 20 :p 0.9})
                                (r/distribution :binomial {:trials 10 :p 0.4}))

             (draw-distribution 0 10 300 :geometric
                                (r/distribution :geometric)
                                (r/distribution :geometric {:p 0.9})
                                (r/distribution :geometric {:p 0.1}))

             (draw-distribution 0 31 600 :hypergeometric
                                (r/distribution :hypergeometric)
                                (r/distribution :hypergeometric {:population-size 500 :number-of-successes 100 :sample-size 25})
                                (r/distribution :hypergeometric {:population-size 500 :number-of-successes 400 :sample-size 30}))

             (draw-distribution 0 32 300 :pascal
                                (r/distribution :pascal)
                                (r/distribution :pascal {:r 30 :p 0.6})
                                (r/distribution :pascal {:r 1 :p 0.5}))

             (draw-distribution 0 6 200 :poisson
                                (r/distribution :poisson)
                                (r/distribution :poisson {:p 0.1})
                                (r/distribution :poisson {:p 0.9}))

             (draw-distribution 0 9 150 :uniform-int
                                (r/distribution :uniform-int {:lower 2 :upper 5})
                                (r/distribution :uniform-int {:lower 1 :upper 2})
                                (r/distribution :uniform-int {:lower 7 :upper 7}))

             (draw-distribution 0 10 160 :zipf
                                (r/distribution :zipf)
                                (r/distribution :zipf {:exponent 1 :number-of-elements 10})
                                (r/distribution :zipf {:exponent 10}))

             (draw-distribution -5 10 460 :gumbel
                                (r/distribution :gumbel)
                                (r/distribution :gumbel {:mu 0.5 :beta 2.0})
                                (r/distribution :gumbel {:mu 3.0 :beta 1.0}))

             (draw-distribution -10 10 300 :laplace
                                (r/distribution :laplace)
                                (r/distribution :laplace {:mu 0.5 :beta 1.0})
                                (r/distribution :laplace {:mu 3.0 :beta 2.5}))

             (draw-distribution -10 10 460 :logistic
                                (r/distribution :logistic)
                                (r/distribution :logistic {:mu 0.5 :s 2.0})
                                (r/distribution :logistic {:mu 3.0 :s 1.0}))

             (draw-distribution -0.1 10 460 :log-normal
                                (r/distribution :log-normal)
                                (r/distribution :log-normal {:scale 1.0 :shape 0.5})
                                (r/distribution :log-normal {:scale 0.5 :shape 1.0}))

             (draw-distribution 0 3 100 :nakagami
                                (r/distribution :nakagami)
                                (r/distribution :nakagami {:mu 5.0 :omega 1.0})
                                (r/distribution :nakagami {:mu 0.5 :omega 1.0}))

             (draw-distribution -5 5 200 :normal
                                (r/distribution :normal)
                                (r/distribution :normal {:mu 0.0 :sd 0.5})
                                (r/distribution :normal {:mu 1.0 :sd 2.0}))

             (draw-distribution 0 5 40 :pareto
                                (r/distribution :pareto)
                                (r/distribution :pareto {:scale 2.0 :shape 0.5})
                                (r/distribution :pareto {:scale 0.5 :shape 2.0}))

             (draw-distribution -4 4 500 :t
                                (r/distribution :t)
                                (r/distribution :t {:degrees-of-freedom 3})
                                (r/distribution :t {:degrees-of-freedom 5}))

             (draw-distribution -2 2 100 :triangular
                                (r/distribution :triangular)
                                (r/distribution :triangular {:a -1.5 :b 1.0 :c 1.5})
                                (r/distribution :triangular {:a -1.0 :b -0.8 :c 1.0}))

             (draw-distribution 0 9 150 :uniform-real
                                (r/distribution :uniform-real {:lower 2 :upper 5})
                                (r/distribution :uniform-real {:lower 0.5 :upper 2.2})
                                (r/distribution :uniform-real {:lower 6.1 :upper 7.1}))

             (draw-distribution -10 10 860 :enumerated-real
                                (r/distribution :enumerated-real {:data (repeatedly 10000 #(int (* 10 (m/sin (r/drand 0 m/TWO_PI)))))}))


             (draw-distribution -10 10 360 :enumerated-int
                                (r/distribution :enumerated-int {:data [-4.0 0 3 4 5]
                                                                 :probabilities [0.1 0.5 0.2 0.05 0.15]}))


             ;; denoise

             (let [c (canvas 512 200)
                   s (for [x (range 512)]
                       (let [xx (m/norm x 0 512 0 m/TWO_PI)]
                         (+ (m/sin xx) (r/randval 0.0 (r/grand 0 0.1)))))
                   ss (t/denoise (t/transformer :packet :daubechies-5) s false)
                   p1 (map v/vec2 (range 512) (map #(+ 100 (* 70.0 %)) s))
                   p2 (map v/vec2 (range 512) (map #(+ 100 (* 70.0 %)) ss))] 
               (with-canvas-> c
                 (set-background bg-color)
                 (set-color 180 200 200)
                 (path p1)
                 (set-stroke 1.8)
                 (set-color :black)
                 (path p2))
               (save-canvas c "t" "denoise"))

             ;; compressed

             (let [c (canvas 512 200)
                   s (for [x (range 512)]
                       (let [xx (m/norm x 0 512 0 m/TWO_PI)]
                         (+ (m/sin xx) (r/randval 0.0 (r/grand 0 0.1)))))
                   tr (t/transformer :packet :haar)
                   ss (t/reverse-1d tr (t/compress 5 (t/forward-1d tr s)))
                   p1 (map v/vec2 (range 512) (map #(+ 100 (* 70.0 %)) s))
                   p2 (map v/vec2 (range 512) (map #(+ 100 (* 70.0 %)) ss))] 
               (with-canvas-> c
                 (set-background bg-color)
                 (set-color 180 200 200)
                 (path p1)
                 (set-stroke 2)
                 (set-color :black)
                 (path p2))
               (save-canvas c "t" "compress"))




;;;


;;;;; Add images

             (def single-list `(sin cos tan cot sec csc asin acos atan acot asec acsc
                                    sinh cosh tanh coth sech csch asinh acosh atanh acoth asech acsch
                                    qsin qcos exp log log10 ln log1p sqrt cbrt qexp qsqrt rqsqrt
                                    erf erfc inv-erf inv-erfc sinc log2 qlog
                                    sq pow2 pow3 safe-sqrt floor ceil round rint abs iabs trunc
                                    frac sfrac low-2-exp high-2-exp round-up-pow2 next-double prev-double
                                    signum sgn sigmoid
                                    gamma log-gamma digamma log-gamma-1p trigamma inv-gamma-1pm1))

             (def interp-list `(quad-interpolation smooth-interpolation wrap lerp cos-interpolation))

             (def fn-list (concat single-list interp-list))

             (defmacro ^:private add-image-examples
               []
               `(do
                  ~@(for [x fn-list]
                      `(add-examples ~x
                                     (example-image ~(str "Plot of " (name x)) ~(str "images/m/" (name x) ".png"))))))

             (add-image-examples)


             ;;

             (def fn-list `(atan asin acos log exp csc sec tanh tan sinh sin cosh cos sqrt sq sqrt1z reciprocal identity))

             (defmacro ^:private add-image-examples
               []
               `(do
                  ~@(for [x fn-list]
                      `(add-examples ~x
                                     (example-image ~(str "Plot of " (name x)) ~(str "images/c/" (name x) ".jpg"))))))

             (add-image-examples)
             ))


