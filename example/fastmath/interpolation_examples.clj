(ns fastmath.interpolation-examples
  (:require [fastmath.interpolation :refer :all]
            [incanter.charts :as c]
            [fastmath.core :as m]
            [metadoc.examples :refer :all]
            [incanter.core :as i]))

(defsnippet fastmath.interpolation save-graph
  "Save graph"
  (let [fname (str "images/i/" (first opts) ".png")
        [x1 x2] params]
    (incanter.core/save (c/function-plot f x1 x2 :y-label "value") (str "docs/" fname) :width 250 :height 250)
    fname))

(add-examples rbf
  (example-snippet "Linear" save-graph :image (rbf :linear) -2.0 2.0)
  (example-snippet "Gaussian" save-graph :image (rbf :gaussian) -2.0 2.0)
  (example-snippet "Multiquadratic" save-graph :image (rbf :multiquadratic) -3.0 3.0)
  (example-snippet "Inverse multiquadratic" save-graph :image (rbf :inverse-multiquadratic) -3.0 3.0)
  (example-snippet "Inverse quadratic" save-graph :image (rbf :inverse-quadratic) -3.0 3.0)
  (example-snippet "Thinplate" save-graph :image (rbf :thinplate) -1.0 3.0)
  (example-snippet "Polyharmonic 3" save-graph :image (rbf :polyharmonic 3) -3.0 3.0)
  (example-snippet "Polyharmonic 4" save-graph :image (rbf :polyharmonic 4) -3.0 3.0)
  (example-snippet "Wendland" save-graph :image (rbf :wendland) -3.0 3.0)
  (example-snippet "Wu" save-graph :image (rbf :wu) -3.0 3.0))


(comment incanter.core/view (c/function-plot (rbf :polyharmonic 2) -3.0 3.0 :y-label "value") :width 250 :height 250)

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
  (example-snippet "Shepard plot, p=5" save-interpolation :image (fn [xs ys] (shepard xs ys 5)) 0 7)
  (example-snippet "Shepard plot, p=0.9" save-interpolation :image (fn [xs ys] (shepard xs ys 0.9)) 0 7))

(add-examples kriging-spline
  (example-snippet "Kriging spline plot" save-interpolation :image kriging-spline 0 7))

(add-examples cubic-spline
  (example-snippet "Cubic spline plot" save-interpolation :image cubic-spline 0 7))

(add-examples spline
  (example-snippet "Spline plot" save-interpolation :image spline 0.69 6.22))

(add-examples neville
  (example-snippet "Neville plot" save-interpolation :image neville 0.65 6.5))

(add-examples divided-difference
  (example-snippet "Divided Difference plot" save-interpolation :image divided-difference 0.65 6.5))

(add-examples loess
  (example-snippet "Loess plot" save-interpolation :image loess 0.69 6.22)
  (example-snippet "Loess plot, bandwidth=0.7, robustness-iters=4" save-interpolation :image (fn [xs ys] (loess xs ys 0.7 4)) 0.69 6.22)
  (example-snippet "Loess plot, bandwidth=0.2, robustness-iters=1" save-interpolation :image (fn [xs ys] (loess xs ys 0.2 1)) 0.69 6.22))

(add-examples microsphere-projection
  (example-snippet "Microsphere projection plot" save-interpolation :image (fn [xs ys] (microsphere-projection xs ys 4 0.5 0.0000001 0.1 2.5 false 0.1)) 0 7))

(add-examples rbf-interpolator
  (example-snippet "RBF - Linear" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :linear))) 0 7)
  (example-snippet "RBF - Gaussian" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :gaussian))) 0 7)
  (example-snippet "RBF - Multiquadratic" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :multiquadratic))) 0 7)
  (example-snippet "RBF - Inverse ultiquadratic" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :inverse-multiquadratic))) 0 7)
  (example-snippet "RBF - Inverse quadratic" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :inverse-quadratic))) 0 7)
  (example-snippet "RBF - Thinplate" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :thinplate))) 0 7)
  (example-snippet "RBF - Polyharmonic 3" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :polyharmonic 3))) 0 7)
  (example-snippet "RBF - Polyharmonic 4" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :polyharmonic 4))) 0 7)
  (example-snippet "RBF - Wendland" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :wendland))) 0 7)
  (example-snippet "RBF - Wu" save-interpolation :image (fn [xs ys] (rbf-interpolator xs ys (rbf :wu))) 0 7))


(comment let [f (fn [xs ys] (rbf-interpolator xs ys (rbf :wu 2)))
              xs [0.69 1.73 2.0 2.28 3.46 4.18 4.84 5.18 5.53 5.87 6.22]
              ft (fn [^double x] (m/sin (* x (* 0.5 (m/cos (inc x))))))
              ys (map ft xs)
              [x1 x2] [0.69 6.22]      ]
         (incanter.core/view (doto (c/function-plot ft 0 7 :y-label "value")
                               (c/set-theme-bw)
                               (c/add-points xs ys)
                               (c/add-function (f xs ys) x1 x2)) :width 600 :height 300))


