(ns fastmath.interpolation-examples
  (:require [fastmath.interpolation :refer :all]
            [fastmath.core :as m]
            [metadoc.examples :refer :all]
            [fastmath.kernel :as k]))

(defsnippet fastmath.interpolation interpolate
  "1d interpolation"
  (let [fun (fn [x] (m/sin (* x (* 0.5 (m/cos (inc x))))))
        xs [0.69 1.73 2.0 2.28 3.46 4.18 4.84 5.18 5.53 5.87 6.22]
        ys (map fun xs)
        interpolator ((apply partial f params) xs ys)]
    (interpolator 5.0)))

(defsnippet fastmath.interpolation interpolate2d
  "2d interpolation"
  (let [fun (fn [x y] (m/sin (* (/ (- x 100.0) 10.0) (m/cos (/ y 20.0)))))
        xs [20 50 58 66 100 121 140 150 160 170 180]
        ys [20 30 58 66 90  121 140 152 170 172 180] 
        vs (partition (count ys) (for [x xs y ys] (fun x y)))
        interpolator ((apply partial f params) xs ys vs)]
    (interpolator 105 155)))

(add-examples akima-spline
  (example-snippet "Usage" interpolate :simple akima-spline)
  (example-image "Akima spline plot" "images/i/akima.png"))

(add-examples divided-difference
  (example-snippet "Usage" interpolate :simple divided-difference)
  (example-image "Divided difference plot" "images/i/divided-difference.png"))

(add-examples linear-smile
  (example-snippet "Usage" interpolate :simple linear-smile)
  (example-image "Linear (Smile) plot" "images/i/linear-smile.png"))

(add-examples linear
  (example-snippet "Usage" interpolate :simple linear)
  (example-image "Linear (Apache) plot" "images/i/linear.png"))

(add-examples loess
  (example-snippet "Usage" interpolate :simple loess)
  (example-image "Loess plot" "images/i/loess.png")
  (example-snippet "Usage (0.7, 2.0)" interpolate :simple loess 0.7 2)
  (example-image "Loess (0.7, 2.0) plot" "images/i/loess2.png")
  (example-snippet "Usage (0.2, 1.0)" interpolate :simple loess 0.2 1)
  (example-image "Loess (0.2, 1.0) plot" "images/i/loess1.png"))

(add-examples neville
  (example-snippet "Usage" interpolate :simple neville)
  (example-image "Neville plot" "images/i/neville.png"))

(add-examples spline
  (example-snippet "Usage" interpolate :simple spline)
  (example-image "Spline plot" "images/i/spline.png"))

(add-examples cubic-spline
  (example-snippet "Usage" interpolate :simple cubic-spline)
  (example-image "Cubic spline plot" "images/i/cubic-spline.png"))

(add-examples kriging-spline
  (example-snippet "Usage" interpolate :simple kriging-spline)
  (example-image "Kriging spline plot" "images/i/kriging-spline.png"))

(add-examples step
  (example-snippet "Usage" interpolate :simple step)
  (example-image "Step plot" "images/i/step.png"))

(add-examples step-after
  (example-snippet "Usage" interpolate :simple step-after)
  (example-image "Step (after) plot" "images/i/step-after.png"))

(add-examples step-before
  (example-snippet "Usage" interpolate :simple step-before)
  (example-image "Step (before) plot" "images/i/step-before.png"))

(add-examples monotone
  (example-snippet "Usage" interpolate :simple monotone)
  (example-image "Monotone plot" "images/i/monotone.png"))

(add-examples polynomial
  (example-snippet "Usage" interpolate :simple polynomial)
  (example-image "Polynomial plot" "images/i/polynomial.png"))

(add-examples b-spline
  (example-snippet "Usage" interpolate :simple b-spline)
  (example-image "B-Spline plot" "images/i/bspline1.png")
  (example-image "B-Spline plot (degree=1)" "images/i/bspline2.png")
  (example-image "B-Spline plot (with knots)" "images/i/bspline3.png"))

(add-examples b-spline-interp
  (example-snippet "Usage" interpolate :simple b-spline-interp)
  (example-image "B-Spline interpolation plot" "images/i/bsplinei1.png")
  (example-image "B-Spline interpolation plot (degree=5)" "images/i/bsplinei2.png")
  (example-image "B-Spline interpolation plot (degree=3, h=6)" "images/i/bsplinei3.png"))

(add-examples microsphere-projection
  (example-snippet "Usage" interpolate :simple microsphere-projection 6 0.1 0.1 0.1 1.5 false 0.01)
  (example-image "Microsphere plot" "images/i/microsphere.png"))

(add-examples rbf
  (example-snippet "Usage" interpolate :simple rbf)
  (example-image "Rbf plot" "images/i/rbf.png")
  (example-snippet "Usage (mattern-c0 kernel)" interpolate :simple rbf (k/rbf :mattern-c0))
  (example-image "Rbf (mattern-c0 kernel) plot" "images/i/rbf1.png")
  (example-snippet "Usage (gaussian kernel, normalized)" interpolate :simple rbf (k/rbf :gaussian) true)
  (example-image "Rbf (gaussian kernel, normalized) plot" "images/i/rbf2.png")
  (example-snippet "Usage (truncated-power kernel)" interpolate :simple rbf (k/rbf :truncated-power 3 0.3))
  (example-image "Rbf (truncated-power kernel) plot" "images/i/rbf3.png")
  (example-snippet "Usage (wendland-53 kernel)" interpolate :simple rbf (k/rbf :wendland-53))
  (example-image "Rbf (wendland-53 kernel) plot" "images/i/rbf4.png"))

(add-examples shepard
  (example-snippet "Usage" interpolate :simple shepard)
  (example-image "Shepard plot" "images/i/shepard.png")
  (example-snippet "Usage (0.9)" interpolate :simple shepard 0.9)
  (example-image "Shepard (0.9) plot" "images/i/shepard1.png"))

;; 2d

(add-examples bicubic
  (example-snippet "Usage" interpolate2d :simple bicubic)
  (example-image "Bicubic plot" "images/i/bicubic.jpg"))

(add-examples piecewise-bicubic
  (example-snippet "Usage" interpolate2d :simple piecewise-bicubic)
  (example-image "Piecewise bicubic plot" "images/i/piecewise-bicubic.jpg"))

(add-examples bilinear
  (example-snippet "Usage" interpolate2d :simple bilinear)
  (example-image "Bilinear plot" "images/i/bilinear.jpg"))

(add-examples bicubic-smile
  (example-snippet "Usage" interpolate2d :simple bicubic-smile)
  (example-image "Bicubic (Smile) plot" "images/i/bicubic-smile.jpg"))

(add-examples cubic-2d
  (example-snippet "Usage" interpolate2d :simple cubic-2d)
  (example-image "Cubic-2d plot" "images/i/cubic-2d.jpg"))

(add-examples microsphere-2d-projection
  (example-snippet "Usage" interpolate2d :simple microsphere-2d-projection 10 0.5 0.0001 0.5 1.5 false 0.1)
  (example-image "Microsphere 2d plot" "images/i/microsphere-2d.jpg"))

(add-examples interpolators-1d-list
  (example "List of names" (keys interpolators-1d-list)))

(add-examples interpolators-2d-list
  (example "List of names" (keys interpolators-2d-list)))


