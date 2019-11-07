(ns fastmath.utils.gen-images
  "Generate images from examples attached to metadata."
  (:require [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.kernel :as k]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.transform :as t]
            [fastmath.interpolation :as i]
            [fastmath.distance :as d]
            [fastmath.fields :as fields]
            [fastmath.easings :as e]
            [fastmath.grid :as g]
            [fastmath.optimization :as opt]
            [cljplot.build :as b]
            [cljplot.core :as cljplot]
            [clojure2d.color :as clr]
            [clojure2d.core :as c2d]
            [clojure2d.extra.utils :as utls]))

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

;; kernels

(doseq [rbf k/rbf-list]
  (save-chart (function-chart (k/rbf rbf)) "k" (str "rbf_" (name rbf)) ".png"))

(doseq [ks k/kernels-list]
  (save-chart (function2d-chart (let [k (k/kernel ks)]
                                  (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" (str "k_" (name ks)) ".jpg"))

(def density-data (repeatedly 200 r/grand))

(defn density-chart
  [k]
  (cljplot/xy-chart {:width 400 :height 200 :background bg-color}
                    (b/series [:vline 0 {:color [60 100 120]}]
                              [:hline 0 {:color [60 100 120]}]
                              [:histogram density-data {:type :lollipops :density? true :palette [[200 200 220]] :bins 35}]
                              [:density density-data {:kernel-type k :color (clr/set-alpha fg-color 180) :area? true
                                                      :kernel-bandwidth 0.25}])
                    (b/update-scale :x :ticks 5)
                    (b/update-scale :y :ticks 5)
                    (b/add-axes :bottom {:ticks {:color fg-color}
                                         :line {:color fg-color}})
                    (b/add-axes :left {:ticks {:color fg-color}
                                       :line {:color fg-color}})))

(doseq [kd k/kernel-density-list]
  (save-chart (density-chart kd) "k" (str "d_" (name kd)) ".jpg"))

(save-chart (function2d-chart (let [k (k/approx (k/kernel :gaussian) 1)]
                                (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" "approx" ".jpg")

(save-chart (function2d-chart (let [k (k/cpd->pd (k/kernel :periodic))]
                                (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" "cpdpd" ".jpg")

(save-chart (function2d-chart (let [k (k/exp (k/kernel :dirichlet) 5.0)]
                                (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" "exp" ".jpg")


(defn ci-chart
  []
  (let [d (k/kernel-density-ci :epanechnikov density-data 0.5)
        r (range -3.0 3.0 0.05)
        top (map #(vector % (second (d %))) r)
        bottom (map #(vector % (last (d %))) r)]
    (cljplot/xy-chart {:width 400 :height 200 :background bg-color}
                      (b/series [:vline 0 {:color [60 100 120]}]
                                [:hline 0 {:color [60 100 120]}]
                                [:histogram density-data {:type :lollipops :density? true :palette [[200 200 220]] :bins 35}]
                                [:ci [top bottom] {:color (clr/set-alpha (clr/darken (clr/darken fg-color)) 200)}]
                                [:function (comp first d) {:domain [-3.0 3.0] :color fg-color}])
                      (b/update-scale :x :ticks 5)
                      (b/update-scale :y :ticks 5)
                      (b/add-axes :bottom {:ticks {:color fg-color}
                                           :line {:color fg-color}})
                      (b/add-axes :left {:ticks {:color fg-color}
                                         :line {:color fg-color}}))))


(save-chart (ci-chart) "k" "ci" ".jpg")

(save-chart (function2d-chart (let [k1 (k/kernel :periodic)
                                    k2 (k/kernel :laplacian)
                                    k (k/mult k1 k2)]
                                (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" "mult" ".jpg")

(save-chart (function2d-chart (let [k1 (k/kernel :periodic)
                                    k2 (k/kernel :laplacian)
                                    k (k/wadd [0.2 0.8] [k1 k2])]
                                (fn [x y] (k [x] [y]))) {:x [-3 3] :y [-3 3]}) "k" "wadd" ".jpg")

(save-chart (function-chart (k/rbf :thin-plate 2 1)) "k" "thin-plate" ".jpg")

;; easings

(doseq [[e ef] e/easings-list]
  (save-chart (function-chart ef [0.0 1.0]) "e" (name e) ".png"))

(save-chart (function-chart (e/out e/sin-in) [0 1]) "e" "out" ".png")
(save-chart (function-chart (e/in-out e/sin-out) [0 1]) "e" "in-out" ".png")
(save-chart (function-chart (e/reflect e/elastic-in-out 0.2) [0 1]) "e" "reflect" ".png")

;; complex

(defn complex-chart
  [f]
  (cljplot/xy-chart {:width 300 :height 300 :background bg-color}
                    (b/series [:complex f])
                    (b/update-scale :x :ticks 5)
                    (b/update-scale :y :ticks 5)
                    (b/add-axes :bottom {:ticks {:color fg-color}
                                         :line {:color fg-color}})
                    (b/add-axes :left {:ticks {:color fg-color}
                                       :line {:color fg-color}})))

(doseq [c `(c/acos c/asin c/atan c/cos c/cosh c/csc c/exp c/log c/reciprocal
                   c/sec c/sin c/sinh c/sq c/sqrt c/sqrt1z c/tan c/tanh)]
  (save-chart (complex-chart (symbol->fn c)) "c" (name c) ".jpg"))

(save-chart (complex-chart identity) "c" "identity" ".jpg")

;; interpolation

(defn ifun
  ^double [^double x]
  (m/sin (* x (* 0.5 (m/cos (inc x))))))

(defn interpolation-chart
  ([inter] (interpolation-chart inter false))
  ([inter r?] (apply interpolation-chart inter (if r? [0.69 6.22] [0 7])))
  ([inter mn mx]
   (let [xs [0.69 1.73 2.0 2.28 3.46 4.18 4.84 5.18 5.53 5.87 6.22]
         ys (map ifun xs)]
     (cljplot/xy-chart {:width 450 :height 250 :background bg-color}
                       (b/series [:function ifun {:domain [0 7] :color fg-color :stroke {:dash [5 5]}}]
                                 [:function (inter xs ys) {:domain [mn mx] :samples 300 :color :white}]
                                 [:scatter (map vector xs ys) {:color :lightgoldenrodyellow}])
                       (b/update-scale :x :ticks 5)
                       (b/update-scale :y :ticks 5)
                       (b/add-axes :bottom {:ticks {:color fg-color}
                                            :line {:color fg-color}})
                       (b/add-axes :left {:ticks {:color fg-color}
                                          :line {:color fg-color}})))))

(save-chart (function-chart ifun [0 7]) "i" "1d" ".png")

(save-chart (interpolation-chart i/akima-spline true) "i" "akima" ".png")
(save-chart (interpolation-chart i/divided-difference true) "i" "divided-difference" ".png")
(save-chart (interpolation-chart i/linear true) "i" "linear" ".png")

(save-chart (interpolation-chart i/loess true) "i" "loess" ".png")
(save-chart (interpolation-chart (partial i/loess 0.7 2) true) "i" "loess2" ".png")
(save-chart (interpolation-chart (partial i/loess 0.2 1) true) "i" "loess1" ".png")

(save-chart (interpolation-chart i/neville true) "i" "neville" ".png")
(save-chart (interpolation-chart i/spline true) "i" "spline" ".png")

(save-chart (interpolation-chart (partial i/microsphere-projection 6 0.1 0.1 0.1 1.5 false 0.01) true) "i" "microsphere" ".png")

(save-chart (interpolation-chart i/cubic-spline) "i" "cubic-spline" ".png")
(save-chart (interpolation-chart i/kriging-spline) "i" "kriging-spline" ".png")
(save-chart (interpolation-chart i/linear-smile) "i" "linear-smile" ".png")

(save-chart (interpolation-chart i/rbf) "i" "rbf" ".png")
(save-chart (interpolation-chart (partial i/rbf (k/rbf :mattern-c0))) "i" "rbf1" ".png")
(save-chart (interpolation-chart (partial i/rbf (k/rbf :gaussian) true)) "i" "rbf2" ".png")
(save-chart (interpolation-chart (partial i/rbf (k/rbf :truncated-power 3 0.3))) "i" "rbf3" ".png")
(save-chart (interpolation-chart (partial i/rbf (k/rbf :wendland-53))) "i" "rbf4" ".png")

(save-chart (interpolation-chart i/shepard) "i" "shepard" ".png")
(save-chart (interpolation-chart (partial i/shepard 0.9)) "i" "shepard1" ".png")

(save-chart (interpolation-chart i/step) "i" "step" ".png")
(save-chart (interpolation-chart i/step-after) "i" "step-after" ".png")
(save-chart (interpolation-chart i/step-before) "i" "step-before" ".png")
(save-chart (interpolation-chart i/monotone) "i" "monotone" ".png")

(defn ifun2d
  [x y]
  (m/sin (* (/ (- x 100.0) 10.0) (m/cos (/ y 20.0)))))

(defn interpolation2d-chart
  [f]
  (let [xs [20 50 58 66 100 121 140 150 160 170 180]
        ys [20 30 58 66 90  121 140 152 170 172 180] 
        vs (partition (count ys) (for [x xs y ys] (ifun2d x y)))]
    (cljplot/xy-chart {:width 400 :height 400 :background bg-color}
                      (b/series [:function-2d (f xs ys vs) {:x [20 180] :y [20 180] :gradient (clr/gradient [bg-color :white])}] 
                                [:scatter (for [x xs y ys] [x y]) {:color :lightgoldenrodyellow :margins {:x [0 0] :y [0 0]}}])
                      ;; (b/update-scale :x :ticks 5)
                      ;; (b/update-scale :y :ticks 5)
                      (b/add-axes :bottom {:ticks {:color fg-color}
                                           :line {:color fg-color}})
                      (b/add-axes :left {:ticks {:color fg-color}
                                         :line {:color fg-color}}))))

(save-chart (interpolation2d-chart i/bicubic) "i" "bicubic" ".jpg")
(save-chart (interpolation2d-chart i/piecewise-bicubic) "i" "piecewise-bicubic" ".jpg")
(save-chart (interpolation2d-chart i/bilinear) "i" "bilinear" ".jpg")
(save-chart (interpolation2d-chart i/bicubic-smile) "i" "bicubic-smile" ".jpg")
(save-chart (interpolation2d-chart i/cubic-2d) "i" "cubic-2d" ".jpg")
(save-chart (interpolation2d-chart (partial i/microsphere-2d-projection 10 0.5 0.0001 0.5 1.5 false 0.1)) "i" "microsphere-2d" ".jpg")

(save-chart (function2d-chart ifun2d {:x [0 200]
                                      :y [0 200]}) "i" "2d" ".jpg")

;; grids

(doseq [gs g/cell-names]
  (let [grid (g/grid gs 40)
        cells (distinct (map (partial g/coords->mid grid) (take 100 (map #(v/add (v/mult % 100) (v/vec2 150 150)) (r/sequence-generator :sphere 2)))))]
    (c2d/save (c2d/with-canvas [c (c2d/canvas 300 300 :highest)]
                (-> (c2d/set-background c bg-color)
                    (c2d/set-color fg-color)
                    (c2d/set-stroke 2))
                (doseq [[x y] cells]
                  (c2d/grid-cell c grid x y true))
                (c2d/set-stroke c 1)
                (doseq [[x y :as mid] cells
                        :let [[ax ay] (g/coords->anchor grid mid)]]
                  (-> c
                      (c2d/set-color 250 100 100)
                      (c2d/crect x y 3 3)
                      (c2d/set-color 100 250 250)
                      (c2d/ellipse ax ay 7 7 true)))
                c) (str "docs/images/g/" (name grid) ".jpg"))))

;; optimization

(defn opt-1d-chart
  ([f d pts]
   (cljplot/xy-chart {:width 400 :height 200 :background bg-color}
                     (b/series [:vline 0 {:color [60 100 120]}]
                               [:hline 0 {:color [60 100 120]}]
                               [:function f {:domain (or d [-3.1 3.1])
                                             :color :white
                                             :samples 300}]
                               [:scatter pts {:color :red :size 8}])
                     (b/update-scale :x :ticks 5)
                     (b/update-scale :y :ticks 5)
                     (b/add-axes :bottom {:ticks {:color fg-color}
                                          :line {:color fg-color}})
                     (b/add-axes :left {:ticks {:color fg-color}
                                        :line {:color fg-color}}))))


(defn opt-2d-chart
  [f d pts]
  (cljplot/xy-chart {:width 300 :height 300 :background bg-color}
                    (b/series [:contour-2d f (merge d {:palette (clr/resample 100 [bg-color :white])
                                                       :contours 100})]
                              [:scatter pts {:color (clr/color :red 140) :size 8}])
                    (b/update-scale :x :ticks 5)
                    (b/update-scale :y :ticks 5)
                    (b/add-axes :bottom {:ticks {:color fg-color}
                                         :line {:color fg-color}})
                    (b/add-axes :left {:ticks {:color fg-color}
                                       :line {:color fg-color}})))


(doseq [ooo [:powell :nelder-mead :multidirectional-simplex :cmaes :gradient :brent]]
  (let [bounds [[-5.0 5.0]]
        f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))
        o1 (opt/minimize ooo f {:bounds bounds})
        o2 (opt/maximize ooo f {:bounds bounds})]
    (save-chart (opt-1d-chart f (first bounds)
                              (map (juxt ffirst second) [o1 o2])) "o" (str (name ooo) "-1d") ".png")))

(doseq [ooo [:powell :nelder-mead :multidirectional-simplex :cmaes :gradient :bobyqa]]
  (let [bounds [[-5.0 5.0] [-5.0 5.0]]
        f (fn [x y] (+ (m/sq (+ (* x x) y -11.0))
                      (m/sq (+ x (* y y) -7.0)))) ;; Himmelblau's function
        o1 (opt/minimize ooo f {:bounds bounds})]
    (save-chart (opt-2d-chart f {:x (first bounds)
                                 :y (second bounds)} [(first o1)]) "o" (str (name ooo) "-2d") ".jpg")))

(let [bounds [[-5.0 5.0] [-5.0 5.0]]
      f (fn [x y] (+ (m/sq (+ (* x x) y -11.0))
                    (m/sq (+ x (* y y) -7.0)))) ;; Himmelblau's function
      bo (nth (opt/bayesian-optimization (fn [x y] (- (f x y))) {:bounds bounds
                                                                :init-points 5
                                                                :utility-function-type :poi}) 20)]
  #_(cljplot/show (opt-2d-chart f {:x (first bounds)
                                   :y (second bounds)} (:xs bo)))
  (save-chart (opt-2d-chart f {:x (first bounds)
                               :y (second bounds)} (:xs bo)) "o" "bo" ".jpg"))

#_(cljplot/show (interpolation2d-chart i/piecewise-bicubic))
