^{:kindly/hide-code true}
(ns interpolation
  (:require [scicloj.kindly.v4.kind :as kind]
            [fastmath.interpolation :as i]
            [fastmath.interpolation.linear :as linear]
            [fastmath.core :as m]
            [fastmath.dev.ggplot :as ggplot]
            [clojisr.v1.r :as R]
            [fastmath.random :as r]
            [fastmath.interpolation.cubic :as cubic]
            [fastmath.calculus :as calc]
            [fastmath.interpolation.acm :as acm]
            [fastmath.interpolation.ssj :as ssj]
            [fastmath.interpolation.step :as step]
            [fastmath.interpolation.shepard :as shepard]
            [fastmath.kernel :as kernel]
            [fastmath.kernel.vector :as kvector]
            [fastmath.interpolation.rbf :as rbf]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.optimization :as opt]
            [fastmath.interpolation.kriging :as kriging]))

;; # Interpolation

;; Interpolation namespace defines the unified API for various interpolation methods. Most of them also extrapolates. Methods include:

;; * 1d interpolation
;; * 2d interpolation on irregular grid points
;; * Multivariate interpolation
;; * Kernel based interpolation
;; * Smoothing

;; All methods are accessible from `fastmath.interpolation` namespace via a multimethod `interpolation`. Additionally each method is implemented as a regular function in the dedicated namespace. `interpolation` returns an interpolant function

;; Both examples below are equivalent:

(require '[fastmath.interpolation :as i]
         '[fastmath.interpolation.linear :as linear])

(def i1 (i/interpolation :linear [1 2 3] [4 5 -1]))
(def i2 (linear/linear [1 2 3] [4 5 -1]))

{:i1 (i1 2.5)
 :i2 (i2 2.5)}

;; List of all possible methods:

(sort (keys (methods i/interpolation)))

;; The following functions and samples will be used as a target to illustrate usage of described method.

;; ## 1d target

;; $$f(x)=\sin\left(\frac{x\cos(x+1)}{2}\right)$$

(defn target-1d [x] (m/sin (* 0.5 x (m/cos (inc x)))))

(target-1d 4.0)

;; Points used in interpolation

(def xs1 [0.5 0.69 1.73 2.0 2.28 3.46 3.5 4.18 4.84 5.18 5.53 5.87 6.22 6.5])
(def ys1 (map target-1d xs1))

^{:kindly/hide-code true}
(-> (ggplot/function+scatter target-1d xs1 ys1 {:x [0.0 7.0]})
    #_(ggplot/title "Target function with sampled points")
    (ggplot/->file))

(-> (ggplot/function+scatter (ssj/b-spline xs1 ys1) xs1 ys1 {:x [0.5 6.5]})
    (ggplot/->file))

;; ## 2d target

;; $$f(x,y)=\sin\left(\frac{x-100}{10}\cos\left(\frac{y}{20}\right)\right)+\frac{x}{100}+\left(\frac{y-100}{100}\right)^2+1$$

(defn target-2d [[x y]] (m/+ 1.0 (m/sin (* (/ (- x 100.0) 10.0) (m/cos (/ y 20.0))))
                          (m// x 100.0)
                          (m/sq (m// (m/- y 100.0) 100.0))))

(target-2d [20 20])

;; ### Grid

;; Points for grid interpolation

(def xs2 [20 25 30 35 40 50 58 66 100 121 140 150 160 170 180])
(def ys2 [20 30 58 66 90  121 140 152 170     180])
(def zss (for [x xs2]
         (for [y ys2]
           (target-2d [x y]))))

^{:kindly/hide-code true}
(let [points (->> (for [x xs2 y ys2] [x y]) (apply map vector))]
  (-> (ggplot/function2d+scatter target-2d (first points) (second points))
      (ggplot/title "Target 2d function with grid points")
      (ggplot/->image)))

;; ### Random

;; For non-grid methods random samples will be used

(def uniform-seed-44 (r/rng :uniform 44))

(def xss (repeatedly 300 #(vector (r/drandom uniform-seed-44 20 180)
                                (r/drandom uniform-seed-44 20 180))))
(def ys3 (map target-2d xss))

^{:kindly/hide-code true}
(-> (ggplot/function2d+scatter target-2d (map first xss) (map second xss))
    (ggplot/title "Target 2d function with random points")
    (ggplot/->image))

;; ### Error

;; The error between original function and interpolant is defined as:

;; $$error_{1d}(f,g)=\|f-g\|=\sqrt{\int_{0.5}^{6.5}|f(x)-g(x)|^2\,dx}$$

(require '[fastmath.calculus :as calc])

(defn error-1d
  [interpolant]
  (m/sqrt (calc/integrate (fn [^double x] (m/sq (m/- (target-1d x) (interpolant x)))) 0.5 6.5)))

(error-1d (linear/linear xs1 ys1))

;; For 2d case the following formula will be used:

;; $$error_{2d}(f,g)=\|f-g\|=\sqrt{\int_{20}^{180}\int_{20}^{180}|f(x,y)-g(x,y)|^2\,dx dy}$$

(defn error-2d
  [interpolant]
  (m/sqrt (calc/cubature (fn [xy] (m/sq (m/- (target-2d xy) (interpolant xy))))
                         [20.0 20.0]
                         [180.0 180.0])))

(error-2d (linear/bilinear xs2 ys2 zss))

;; ## 1d

^{:kindly/hide-code true}
(defn chart-1d
  [f title]
  (-> (ggplot/function+scatter [["interp" f]
                                ["orig" target-1d]] xs1 ys1 {:x [0.0 7.5]})
      (ggplot/ylim [-1.2 1.2])
      (ggplot/title title)
      (ggplot/->image)))

;; ### Linear

;; Linear piecewise interpolation and extrapolation. Extrapolation uses a slope from the boundaries.
;; See more on [Wikipedia](https://en.wikipedia.org/wiki/Linear_interpolation)

(require '[fastmath.interpolation.linear :as linear])

(def linear (linear/linear xs1 ys1))

(linear 4.0)

(error-1d linear)

^{:kindly/hide-code true}
(chart-1d linear "Linear interpolation")

;; ### Cubic

;; Natural cubic spline (second derivatives at boundary points have value $0$) interpolation and extrapolation.
;; See more on [Wikipedia](https://en.wikipedia.org/wiki/Spline_interpolation)

(require '[fastmath.interpolation.cubic :as cubic])

(def cubic (cubic/cubic xs1 ys1))

(cubic 4.0)

(error-1d cubic)

^{:kindly/hide-code true}
(chart-1d cubic "Cubic spline interpolation")

;; ### Akima

;; See more on [Wikipedia](https://en.wikipedia.org/wiki/Akima_spline)

(require '[fastmath.interpolation.acm :as acm])

(def akima (acm/akima xs1 ys1))

(akima 4.0)

(error-1d akima)

^{:kindly/hide-code true}
(chart-1d akima "Akima interpolation")

;; ### Neville

;; See more on [Wikipedia](https://en.wikipedia.org/wiki/Neville%27s_algorithm)

(require '[fastmath.interpolation.acm :as acm])

(def neville (acm/neville xs1 ys1))
(neville 4.0)
(error-1d neville)

^{:kindly/hide-code true}
(chart-1d neville "Neville interpolation")

;; ### Barycentric

;; Rational interpolation as described in [Numerical Recipes ch. 3.4](https://numerical.recipes/book.html). The `order` (default $1$) parameter contols number of points used to calculate weights. Higher order means better accuracy.

(require '[fastmath.interpolation.barycentric :as barycentric])

(defn barycentric
  ([] (barycentric/barycentric xs1 ys1))
  ([order] (barycentric/barycentric xs1 ys1 {:order order})))

((barycentric) 4.0)

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["order" "error" "barrycentric(4.0)" "error at 4.0"]
 :row-vectors (for [i (range 6)
                    :let [v ((barycentric i) 4.0)]]
                [i (error-1d (barycentric i))
                 v (m/abs (- v (target-1d 4.0)))])}

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (barycentric 0) "Barycentric, order=0")
                (chart-1d (barycentric 1) "Barycentric, order=1")]
               [(chart-1d (barycentric 2) "Barycentric, order=2")
                (chart-1d (barycentric 3) "Barycentric, order=3")]
               [(chart-1d (barycentric 4) "Barycentric, order=4")
                (chart-1d (barycentric 5) "Barycentric, order=5")]]}

;; ### B-spline

;;

(require '[fastmath.interpolation.ssj :as ssj])

(defn b-spline
  ([] (ssj/b-spline xs1 ys1))
  ([degree] (b-spline degree nil))
  ([degree hp1] (ssj/b-spline xs1 ys1 {:degree degree :hp1 hp1})))

((b-spline) 4.0)

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (b-spline 1) "B-spline, degree=1")
                (chart-1d (b-spline 2) "B-spline, degree=2")]
               [(chart-1d (b-spline 3) "B-spline, degree=3")
                (chart-1d (b-spline 5) "B-spline, degree=5")]
               [(chart-1d (b-spline 3 7) "B-spline, degree=3, hp1=7")
                (chart-1d (b-spline 5 7) "B-spline, degree=5, hp1=7")]]}

;; ### Divided difference

(require '[fastmath.interpolation.acm :as acm])

(def divided-difference (acm/divided-difference xs1 ys1))

(divided-difference 4.0)

(error-1d divided-difference)

^{:kindly/hide-code true}
(chart-1d divided-difference "Divided difference interpolation")

;; ### Polynomial

(require '[fastmath.interpolation.ssj :as ssj])

(def polynomial (ssj/polynomial xs1 ys1))

(polynomial 4.0)

(error-1d polynomial)

^{:kindly/hide-code true}
(chart-1d polynomial "Polynomial interpolation")


;; ### Monotone

(require '[fastmath.interpolation.monotone :as monotone])

(def monotone (monotone/monotone xs1 ys1))

(monotone 4.0)

(error-1d monotone)

^{:kindly/hide-code true}
(chart-1d monotone "Monotone interpolation")

;; ### Step

(require '[fastmath.interpolation.step :as step])

(defn step
  ([] (step/step xs1 ys1))
  ([point] (step/step xs1 ys1 {:point point})))
(def step-before (step/step-before xs1 ys1))
(def step-after (step/step-after xs1 ys1))

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["method" "error" "value at 4.0"]
 :row-vectors [["step-before" (error-1d step-before) (step-before 4.0)]
               ["step-after" (error-1d step-after) (step-after 4.0)]
               ["step" (error-1d (step)) ((step) 4.0)]
               ["step (point=0.55)" (error-1d (step 0.55)) ((step 0.55) 4.0)]
               ["step (point=0.25)" (error-1d (step 0.25)) ((step 0.25) 4.0)]
               ["step (point=0.75)" (error-1d (step 0.75)) ((step 0.75) 4.0)]]}


^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d step-before "Step before")
                (chart-1d step-after "Step after")]
               [(chart-1d (step) "Step (point=0.5, default)")
                (chart-1d (step 0.55) "Step (point=0.55)")]
               [(chart-1d (step 0.25) "Step (point=0.25)")
                (chart-1d (step 0.75) "Step (point=0.75)")]]}

;; ### Loess

(require '[fastmath.interpolation.acm :as acm])

(defn loess
  ([] (acm/loess xs1 ys1))
  ([bandwidth] (acm/loess xs1 ys1 {:bandwidth bandwidth})))

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (loess) "Loess (default)")
                (chart-1d (loess 0.8) "Loess (bw=0.8)")]
               [(chart-1d (loess 0.2) "Loess (bw=0.2)")
                (chart-1d (loess 1.0) "Loess (bw=1.0)")]]}

;; ### Cubic smoothing

(require '[fastmath.interpolation.ssj :as ssj])

(defn cubic-smoothing
  ([] (ssj/cubic-smoothing xs1 ys1))
  ([rho] (ssj/cubic-smoothing xs1 ys1 {:rho rho})))

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (cubic-smoothing)     "Cubic smoothing (default, rho=1.0)")
                (chart-1d (cubic-smoothing 0.8) "Cubic smoothing (rho=0.8)")]
               [(chart-1d (cubic-smoothing 0.2) "Cubic smoothing (rho=0.2)")
                (chart-1d (cubic-smoothing 0.0) "Cubic smoothing (rho=0.0)")]]}

;; ## 2d grid

^{:kindly/hide-code true}
(defn chart-2d
  ([f title]
   (-> (ggplot/function2d f {:steps 100 :x [20 180] :y [20 180]})
       (ggplot/title title)
       (ggplot/->image))))

;; ### Bilinear

(require '[fastmath.interpolation.linear :as linear])

(def bilinear (linear/bilinear xs2 ys2 zss))

(chart-2d bilinear "Bilinear")

(error-2d bilinear)

;; ### Bicubic

(require '[fastmath.interpolation.acm :as acm])

(def bicubic (acm/bicubic xs2 ys2 zss))

(chart-2d bicubic "Bicubic")

(error-2d bicubic)

;; ### Cubic 2d

(require '[fastmath.interpolation.cubic :as cubic])

(def cubic-2d (cubic/cubic-2d xs2 ys2 zss))

(chart-2d cubic-2d "Cubic 2d")

(error-2d cubic-2d)

;; ## Multivariate and kernel based

;; ### Microsphere projection

(require '[fastmath.interpolation.acm :as acm])

^{:kindly/kind :kind/table :kindly/hide-code true}
[[(chart-1d (acm/microsphere-projection xs1 ys1) "Microsphere projection (1d)")
  (chart-2d (acm/microsphere-projection xss ys3) "Microsphere projection (2d)")]]

(error-1d (acm/microsphere-projection xs1 ys1))
(error-2d (acm/microsphere-projection xss ys3))

;; ### Shepard

(require '[fastmath.interpolation.shepard :as shepard])

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (shepard/shepard xs1 ys1) "Shepard (1d)")
                (chart-2d (shepard/shepard xss ys3) "Shepard (2d)")]
               [(chart-1d (shepard/shepard xs1 ys1 {:p 3}) "Shepard (1d, p=3)")
                (chart-2d (shepard/shepard xss ys3 {:p 3}) "Shepard (2d, p=3)")]]}

(defn chart-2d->file
  ([f]
   (-> (ggplot/function2d f {:steps 100 :x [20 180] :y [20 180]})
       (ggplot/->file))))

;; ### Radial Basis Function

(require '[fastmath.interpolation.rbf :as rbf])

(defn chart-f [f title] (-> (ggplot/function f {:x [-5 5]})
                            (ggplot/title title)
                            (ggplot/->image)))

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-f (kernel/rbf :gaussian) "Gaussian")
                (chart-f (kernel/rbf :matern-c2) "Matern C2, 3/2")]
               [(chart-f (kernel/rbf :gaussians-laguerre-22) "Gaussians-Laguerre (dim=2,deg=2)")
                (chart-f (kernel/rbf :thin-plate) "Thin plate")]]}

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussian)) "Gaussian (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1})) "Gaussian (2d, shape=0.1)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :matern-c2)) "Matern C2, 3/2 (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :matern-c2 {:shape 0.15}))
                          "Matern C2, 3/2 (2d, shape=0.15)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussians-laguerre-22)) "Gaussians-Laguerre (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussians-laguerre-22 {:shape 0.07}))
                          "Gaussians-Laguerre (2d, dim=2, deg=2, shape=0.07)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :thin-plate)) "Thin plate (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :thin-plate)) "Thin plate (2d)")]]}

;; #### Polynomial term

(defn polynomial-terms-1d [^double x]
  [1.0 x (m/sq x)])

(defn polynomial-terms-2d [[^double x ^double y]]
  [1.0 x y (m/* x y) (m/sq x) (m/sq y)])

^{:kindly/kind :kind/table :kindly/hide-code true}
{:column-names ["" ""]
 :row-vectors [[(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussian)
                                   {:polynomial-terms polynomial-terms-1d})
                          "Gaussian (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1})
                                   {:polynomial-terms polynomial-terms-2d})
                          "Gaussian (2d, shape=0.1)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :matern-c2)
                                   {:polynomial-terms polynomial-terms-1d})
                          "Matern C2, 3/2 (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :matern-c2 {:shape 0.15})
                                   {:polynomial-terms polynomial-terms-2d})
                          "Matern C2, 3/2 (2d, shape=0.15)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussians-laguerre-22)
                                   {:polynomial-terms polynomial-terms-1d})
                          "Gaussians-Laguerre (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussians-laguerre-22 {:shape 0.07})
                                   {:polynomial-terms polynomial-terms-2d})
                          "Gaussians-Laguerre (2d, dim=2, deg=2, shape=0.07)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :thin-plate) {:polynomial-terms polynomial-terms-1d})
                          "Thin plate (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :thin-plate) {:polynomial-terms polynomial-terms-2d})
                          "Thin plate (2d)")]]}


(error-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1})))
(error-2d (rbf/rbf xss ys3 (kernel/rbf :matern-c2 {:shape 0.15})))
(error-2d (rbf/rbf xss ys3 (kernel/rbf :gaussians-laguerre-22 {:shape 0.07})))
(error-2d (rbf/rbf xss ys3 (kernel/rbf :thin-plate)))

;; #### Smoothing

^{:kindly/kind :kind/table :kindly/hide-code true}
{:row-vectors [[(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussian) {:lambda 0.0
                                                                   :polynomial-terms polynomial-terms-1d})
                          "Gaussian (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1}) {:lambda 0.0
                                                                                :polynomial-terms polynomial-terms-2d}) "Gaussian (2d, shape=0.1)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussian) {:lambda 1.0
                                                                   :polynomial-terms polynomial-terms-1d})
                          "Gaussian (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1}) {:lambda 1.0
                                                                                :polynomial-terms polynomial-terms-2d}) "Gaussian (2d, shape=0.1)")]
               [(chart-1d (rbf/rbf xs1 ys1 (kernel/rbf :gaussian) {:lambda 5.0
                                                                   :polynomial-terms polynomial-terms-1d})
                          "Gaussian (1d)")
                (chart-2d (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1}) {:lambda 5.0
                                                                                :polynomial-terms polynomial-terms-2d}) "Gaussian (2d, shape=0.1)")]]}

#_(defn chart-f->file [f] (-> (ggplot/function f {:x [0 7]})                         
                              (ggplot/->file)))


#_(chart-2d->file (rbf/rbf xss ys3 (kernel/rbf :gaussian {:shape 0.1}) {:lambda 20}))

#_(opt/minimize :lbfgsb (fn [sh] (error-2d (rbf/rbf xss ys3 (keranel/rbf :gaussian {:shape sh}))))
                {:bounds [[0.01 0.2]]})

;; ### Kriging

;; #### Variograms

(require '[fastmath.kernel.variogram :as variogram]
         '[fastmath.interpolation.kriging :as kriging])

(defn svar-image [f emp title]
  (let [x (map :h emp)
        y (map :gamma emp)]
    (-> (ggplot/function+scatter f x y)
        (ggplot/title title)
        (ggplot/ylim [0 nil])
        (ggplot/->file))))

(def empirical-matheron-1d (variogram/empirical xs1 ys1))

(def empirical-matheron (variogram/empirical xss ys3 {:size 20}))
empirical-matheron

(let [x (map :h empirical-matheron)
      y (map :gamma empirical-matheron)]
  (-> (ggplot/function+scatter (linear/linear x y) x y)
      (ggplot/title "Classical (Matheron) variogram estimation"                    )
      (ggplot/ylim [0 nil])
      (ggplot/->image)))

(def empirical-cressie (variogram/empirical xss ys3 {:estimator :cressie :size 20}))

(let [x (map :h empirical-cressie)
      y (map :gamma empirical-cressie)]
  (-> (ggplot/function+scatter (linear/linear x y) x y)
      (ggplot/title "Cressie variogram estimation")
      (ggplot/ylim [0 nil])
      (ggplot/->image)))

(def empirical-highly-robust (variogram/empirical xss ys3 {:estimator :highly-robust :size 20
                                                         :remove-outliers? true}))
empirical-highly-robust

(def empirical-quantile (variogram/empirical xss ys3 {:estimator :quantile :size 50
                                                    :quantile 0.92}))



(def empirical-M-robust (variogram/empirical xss ys3 {:estimator :m-robust :size 50}))

(let [x (map :h empirical-highly-robust)
      y (map :gamma empirical-highly-robust)]
  (-> (ggplot/function+scatter (linear/linear x y) x y)
      (ggplot/title "Highly robust (quantiles) variogram estimation")
      (ggplot/ylim [0 nil])
      (ggplot/->image)))

;; #### Semi-variograms

(def variogram-linear  (variogram/fit empirical-quantile :linear))
(def variogram-gaussian (variogram/fit empirical-highly-robust :gaussian))
(def variogram-pentaspherical (variogram/fit empirical-highly-robust :pentaspherical))
(def variogram-rbf-wendland-2-3 (variogram/fit empirical-highly-robust (kernel/rbf :wendland {:s 2 :k 3})))

(def variogram-superspherical (variogram/fit empirical-matheron :superspherical))
(def variogram-superspherical-1d (variogram/fit empirical-matheron-1d :tplstable {:order 1.9 :defaults {:beta 14.0}}))

(((variogram/->superspherical 1.0) {:nugget 0.1 :psill 0.5 :range 1.0}) 0.4)

(ggplot/->file (ggplot/function variogram-superspherical-1d {:x [0 5]}))

(variogram-superspherical-1d)

^{:kindly/kind :kind/table :kindly/hide-code true}
{:row-vectors [[(svar-image variogram-linear empirical-quantile "Linear")
                (svar-image variogram-gaussian empirical-highly-robust "Linear")]
               [(svar-image variogram-pentaspherical empirical-highly-robust "Linear")
                (svar-image variogram-rbf-wendland-2-3 empirical-highly-robust "Linear")]]}

(def kriging-linear (kriging/kriging xss ys3 variogram-linear))
(def kriging-gaussian (kriging/kriging xss ys3 variogram-gaussian))
(def kriging-pentaspherical (kriging/kriging xss ys3 variogram-pentaspherical))
(def kriging-rbf-wendland-2-3 (kriging/kriging xss ys3 variogram-rbf-wendland-2-3))

(def kriging-superspherical (kriging/kriging xss ys3 variogram-superspherical))

(-> (ggplot/function2d kriging-superspherical {:x [20 180] :y [20 180] :steps 100})
    (ggplot/->file))

(error-2d kriging-linear)

^{:kindly/kind :kind/table :kindly/hide-code true}
{:row-vectors [[(-> (ggplot/function2d kriging-linear {:x [20 180] :y [20 180] :steps 100})
                    (ggplot/->file))
                (-> (ggplot/function2d kriging-gaussian {:x [20 180] :y [20 180] :steps 100})
                    (ggplot/->file))]
               [(-> (ggplot/function2d kriging-pentaspherical {:x [20 180] :y [20 180] :steps 100})
                    (ggplot/->file))
                (-> (ggplot/function2d kriging-rbf-wendland-2-3 {:x [20 180] :y [20 180] :steps 100})
                    (ggplot/->file))]]}

(defn target [^double n ^double s ^double r]
  (let [l (variogram/gaussian {:nugget n :sill s :range r})
        k (kriging/kriging xss ys3 l)]
    (error-2d k)))

(def vl (variogram/linear {:nugget 0.03 :sill 0.5 :range 14.0}))
(svar-image vl empirical-cressie "Linear")

(error-2d (kriging/kriging xss ys3 vl))




(opt/scan-and-minimize :lbfgsb target {:bounds [[0.0 0.2]
                                                [0.0 1.0]
                                                [10.0 100.0]]})

(-> (ggplot/function2d (kriging/kriging xss ys3 vl) {:x [20 180] :y [20 180] :steps 100})
    (ggplot/->file))

;; #### Smoothing

(-> (ggplot/function2d (kriging/kriging xss ys3 variogram-pentaspherical {:error 5})
                       {:x [20 180] :y [20 180] :steps 100})
    (ggplot/->image))

;; ### Gaussian processes

(let [k (kvector/triangular 2)]
  (-> (ggplot/function2d (fn [[x y]] (k x y)) {:x [-6 6] :y [-6 6]}) ;; tst
      (ggplot/->image))
  )
