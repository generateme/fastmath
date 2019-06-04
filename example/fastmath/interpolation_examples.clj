(ns fastmath.interpolation-examples
  (:require [fastmath.interpolation :refer :all]
            [incanter.charts :as c]
            [fastmath.core :as m]
            [metadoc.examples :refer :all]
            [fastmath.kernel :as k]
            [incanter.core :as i]))

(defsnippet fastmath.interpolation save-graph
  "Save graph"
  (let [fname (str "images/i/" (first opts) ".png")
        [x1 x2] params]
    (incanter.core/save (c/function-plot f x1 x2 :y-label "value") (str "docs/" fname) :width 250 :height 250)
    fname))

(defsnippet fastmath.interpolation save-interpolation
  "Save interpolation graph"
  (let [xs [0.69 1.73 2.0 2.28 3.46 4.18 4.84 5.18 5.53 5.87 6.22]
        ft (fn [^double x] (m/sin (* x (* 0.5 (m/cos (inc x))))))
        ys (map ft xs)
        [x1 x2] params
        fname (str "images/i/" (first opts) ".png")]
    (incanter.core/save (doto (c/function-plot ft 0 7 :y-label "value")
                          (c/set-theme-bw)
                          (c/add-points xs ys)
                          (c/add-function (f xs ys) x1 x2)) (str "docs/" fname) :width 600 :height 300)
    fname))

(add-examples akima-spline
  (example-snippet "Akima plot" save-interpolation :image akima-spline 0.69 6.22))

(add-examples linear-smile
  (example-snippet "Linear (SMILE version) plot" save-interpolation :image linear-smile 0 7))

(add-examples linear
  (example-snippet "Linear (Apache Commons version) plot" save-interpolation :image linear 0.69 6.22))

(add-examples shepard
  (example-snippet "Shepard plot" save-interpolation :image shepard 0 7)
  (example-snippet "Shepard plot, p=5" save-interpolation :image (fn [xs ys] (shepard 5 xs ys)) 0 7)
  (example-snippet "Shepard plot, p=0.9" save-interpolation :image (fn [xs ys] (shepard 0.9 xs ys)) 0 7))

(add-examples kriging-spline
  (example "Usage" {:test-value -0.07}
    (let [interpolator (kriging-spline [2 5 9 10 11] [0.4 1.0 -1.0 -0.5 0.0])]
      (m/approx (interpolator 7.0))))
  (example-snippet "Kriging spline plot" save-interpolation :image kriging-spline 0 7))

(add-examples cubic-spline
  (example-snippet "Cubic spline plot" save-interpolation :image cubic-spline 0 7))

(add-examples spline
  (example-snippet "Spline plot" save-interpolation :image spline 0.69 6.22))

(add-examples neville
  (example-snippet "Neville plot" save-interpolation :image neville 0.65 6.5))

(add-examples divided-difference
  (example-snippet "Divided Difference plot" save-interpolation :image divided-difference 0.65 6.5))

(add-examples step-after
  (example-snippet "Step function plot" save-interpolation :image step-after 0 7))

(add-examples step-before
  (example-snippet "Step function plot" save-interpolation :image step-before 0 7))

(add-examples loess
  (example-snippet "Loess plot" save-interpolation :image loess 0.69 6.22)
  (example-snippet "Loess plot, bandwidth=0.7, robustness-iters=4" save-interpolation :image (fn [xs ys] (loess 0.7 4 xs ys)) 0.69 6.22)
  (example-snippet "Loess plot, bandwidth=0.2, robustness-iters=1" save-interpolation :image (fn [xs ys] (loess 0.2 1 xs ys)) 0.69 6.22))

(add-examples microsphere-projection
  (example-snippet "Microsphere projection plot" save-interpolation :image (fn [xs ys] (microsphere-projection 4 0.5 0.0000001 0.1 2.5 false 0.1 xs ys)) 0 7))

(add-examples cubic-2d
  (example "Usage" {:test-value 4.68}
    (let [interpolator (cubic-2d [2 5 9] [2 3 10] [[4 0 2]
                                                   [-1 2 -2]
                                                   [-2 0 1]])]
      (m/approx (interpolator 5.0 5.0))))
  (example "Array layout"
    (let [intrp (cubic-2d [2 5] [1 6] [[-1 -2]
                                       [3 4]])]
      [(intrp 2 1)
       (intrp 2 6)
       (intrp 5 1)
       (intrp 5 6)])))

(add-examples rbf
  (example-snippet "RBF - Linear" save-interpolation :image (fn [xs ys] (rbf (k/rbf :linear) xs ys)) 0 7)
  (example-snippet "RBF - Gaussian" save-interpolation :image (fn [xs ys] (rbf (k/rbf :gaussian) xs ys)) 0 7)
  (example-snippet "RBF - Multiquadratic" save-interpolation :image (fn [xs ys] (rbf (k/rbf :multiquadratic) xs ys)) 0 7)
  (example-snippet "RBF - Inverse ultiquadratic" save-interpolation :image (fn [xs ys] (rbf (k/rbf :inverse-multiquadratic) xs ys)) 0 7)
  ;; (example-snippet "RBF - Inverse quadratic" save-interpolation :image (fn [xs ys] (rbf (k/rbf :inverse-quadratic) xs ys)) 0 7)
  (example-snippet "RBF - Thinplate" save-interpolation :image (fn [xs ys] (rbf (k/rbf :thin-plate) xs ys)) 0 7)
  ;; (example-snippet "RBF - Polyharmonic 3" save-interpolation :image (fn [xs ys] (rbf (k/rbf :polyharmonic 3) xs ys)) 0 7)
  ;; (example-snippet "RBF - Polyharmonic 4" save-interpolation :image (fn [xs ys] (rbf (k/rbf :polyharmonic 4) xs ys)) 0 7)
  ;; (example-snippet "RBF - Wendland" save-interpolation :image (fn [xs ys] (rbf (k/rbf :wendland) xs ys)) 0 7)
  ;; (example-snippet "RBF - Wu" save-interpolation :image (fn [xs ys] (rbf (k/rbf :wu) xs ys)) 0 7)
  )


(add-examples interpolators-1d-list
  (example "List of names" (keys interpolators-1d-list)))

(add-examples interpolators-2d-list
  (example "List of names" (keys interpolators-2d-list)))

(comment let [f (fn [xs ys] (rbf xs ys (k/rbf :wu 2)))
              xs [0.69 1.73 2.0 2.28 3.46 4.18 4.84 5.18 5.53 5.87 6.22]
              ft (fn [^double x] (m/sin (* x (* 0.5 (m/cos (inc x))))))
              ys (map ft xs)
              [x1 x2] [0.69 6.22]]
         (incanter.core/view (doto (c/function-plot ft 0 7 :y-label "value")
                               (c/set-theme-bw)
                               (c/add-points xs ys)
                               (c/add-function (f xs ys) x1 x2)) :width 600 :height 300))


