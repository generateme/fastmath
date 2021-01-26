(ns fastmath.optimization-examples
  (:require [fastmath.optimization :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.core :as m]))

(add-examples minimize
  (example "1d function"
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))]
      {:powell (minimize :powell f {:bounds bounds})
       :brent (minimize :brent f {:bounds bounds})
       :brent-with-initial-point (minimize :brent f {:bounds bounds :initial [2.0]})}))
  (example "2d function"
    (let [bounds [[-5.0 5.0] [-5.0 5.0]]
          ;; Himmelblau's function
          f (fn [x y] (+ (m/sq (+ (* x x) y -11.0))
                        (m/sq (+ x (* y y) -7.0))))]
      {:bobyqa (minimize :bobyqa f {:bounds bounds})
       :gradient (minimize :gradient f {:bounds bounds})
       :bfgs (minimize :bfgs f {:bounds bounds})}))
  (example "With stats"
    (minimize :gradient #(m/sin %) {:bounds [[-5 5]]
                                    :stats? true})))

(add-examples scan-and-minimize
  (example "1d function"
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))]
      {:powell (scan-and-minimize :powell f {:bounds bounds})
       :brent (scan-and-minimize :brent f {:bounds bounds})
       :bfgs (scan-and-minimize :bfgs f {:bounds bounds})}))
  (example "2d function"
    (let [bounds [[-5.0 5.0] [-5.0 5.0]]
          ;; Himmelblau's function
          f (fn [x y] (+ (m/sq (+ (* x x) y -11.0))
                        (m/sq (+ x (* y y) -7.0))))]
      {:bobyqa (scan-and-minimize :bobyqa f {:bounds bounds})
       :gradient (scan-and-minimize :gradient f {:bounds bounds})
       :bfgs (scan-and-minimize :bfgs f {:bounds bounds})})))

(add-examples maximize
  (example
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))]
      {:powell (maximize :powell f {:bounds bounds})
       :brent (maximize :brent f {:bounds bounds})
       :bfgs (maximize :bfgs f {:bounds bounds})})))

(add-examples maximizer
  (example
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))
          optimizer (maximizer :cmaes f {:bounds bounds})]
      {:optimizer optimizer
       :run-1 (optimizer)
       :run-2 (optimizer [4.5])
       :run-3 (optimizer [-4.5])})))

(add-examples minimizer
  (example
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))
          optimizer (minimizer :brent f {:bounds bounds})]
      {:optimizer optimizer
       :run-1 (optimizer)
       :run-2 (optimizer [4.5])
       :run-3 (optimizer [-4.5])})))

(add-examples scan-and-maximize
  (example
    (let [bounds [[-5.0 5.0]]
          f (fn [x] (+ (* 0.2 (m/sin (* 10.0 x))) (/ (+ 6.0 (- (* x x) (* 5.0 x))) (inc (* x x)))))]
      {:powell (scan-and-maximize :powell f {:bounds bounds})
       :brent (scan-and-maximize :brent f {:bounds bounds})
       :bfgs (scan-and-maximize :bfgs f {:bounds bounds})})))

(add-examples bayesian-optimization
  (example
    (let [bounds [[-5.0 5.0] [-5.0 5.0]]
          ;; Himmelblau's function
          f (fn [x y] (+ (m/sq (+ (* x x) y -11.0))
                        (m/sq (+ x (* y y) -7.0))))]
      (nth (bayesian-optimization (fn [x y] (- (f x y))) {:bounds bounds
                                                         :init-points 5
                                                         :utility-function-type :poi}) 10)))
  (example-image "Bayesian optimization points" "images/o/bo.jpg"))

(defmacro add-1d-images
  []
  `(add-examples minimize
     ~@(for [x [:powell :nelder-mead :multidirectional-simplex :cmaes :gradient :brent]]
         `(example-image ~(str "min/max of f using `" x "` optimizer") ~(str "images/o/" (name x) "-1d.png")))))

(add-1d-images)

(defmacro add-2d-images
  []
  `(add-examples minimize
     ~@(for [x [:powell :nelder-mead :multidirectional-simplex :cmaes :gradient :bobyqa]]
         `(example-image ~(str "min/max of f using `" x "` optimizer") ~(str "images/o/" (name x) "-2d.jpg")))))

(add-2d-images)
