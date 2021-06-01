(ns fastmath.vector-examples
  (:require [fastmath.vector :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.core :as m]))

(add-examples vec2
  (example-session "Usage" (vec2 0.5 -0.5) (vec2))
  (example "Destructuring"
    (let [[^double x ^double y] (vec2 4.3 2.2)]
      (+ x y)))
  (example "As function"
    [((vec2 11 22) 0) ((vec2 11 22) 1)])
  (example-session "As sequence"
    (map #(* 2.0 ^double %) (vec2 1 2))
    (reduce clojure.core/+ (vec2 6 4))
    (cons 6.0 (vec2 3 4))
    (first (vec2 3 4))
    (second (vec2 3 4))
    (last (vec2 3 4))
    (nth (vec2 3 4) 1)
    (filter clojure.core/neg? (vec2 -1 2)))
  (example "Count" (count (vec2 4 3))))

(add-examples vec3
  (example-session "Usage"
    (vec3)
    (vec3 0.5 -0.5 1.0)
    (let [v (vec2 1 2)]
      (vec3 v -1.0))))

(add-examples vec4
  (example-session "Usage"
    (vec4)
    (vec4 0.5 -0.5 1.0 -1.0)
    (let [v (vec2 1 2)]
      (vec4 v -1.0 0.1))
    (let [v (vec3 0 1 2)]
      (vec4 v 0.1))))

(add-examples array-vec
  (example-session "Usage"
    (array-vec [1 2 3 4 5 6 7])
    (array-vec (range 0.0 1.0 0.25))
    (array-vec 11))
  (example-session "Operations"
    (nth (array-vec [9 8 7 6]) 2)
    (count (array-vec (range 0.1 1.0 0.05)))
    (seq (array-vec [1 2]))))

(add-examples ediv
  (example "Usage" (ediv (vec2 1 3) (vec2 2 4))))

(add-examples average-vectors
  (example-session "Usage"
    (average-vectors [[1 2] [0 1] [3 4] [1 2] [4 -1]])
    (average-vectors (vec2 0 0) [(vec2 1 1) (vec2 1 1) (vec2 1 1)])))

(add-examples dist
  (example-session "Usage"
    (dist (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-sq
  (example-session "Usage"
    (dist-sq (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-sq [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-abs
  (example-session "Usage"
    (dist-abs (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-abs [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-cheb
  (example-session "Usage"
    (dist-cheb (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-cheb [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-discrete
  (example-session "Usage"
    (dist-discrete (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-discrete [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-canberra
  (example-session "Usage"
    (dist-canberra (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-canberra [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-emd
  (example-session "Usage"
    (dist-emd (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-emd [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples sim-cos
  (example-session "Usage"
    (sim-cos (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (sim-cos [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples dist-ang
  (example-session "Usage"
    (dist-ang (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
    (dist-ang [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5])))

(add-examples distances
  (example "List of distances" (sort (keys distances))))

(add-examples normalize
  (example "Usage" (normalize (vec2 1.0 -1.0))))

(add-examples set-mag
  (example-session "Usage"
    (set-mag (vec2 0.22 0.22) (m/sqrt 2.0))
    (set-mag (vec2 1.0 1.0) 0.0)))

(add-examples limit
  (example-session "Usage"
    (limit (vec3 1.0 1.0 1.0) 1.0)
    (limit (vec3 1.0 1.0 1.0) 2.0)))

(add-examples angle-between
  (example-session "Usage"
    (m/degrees (angle-between (vec3 1.0 0.0 0.0) (vec3 0.0 1.0 0.0)))
    (m/degrees (angle-between (vec (repeatedly 50 rand)) (vec (repeatedly 50 rand))))))

(add-examples relative-angle-between
  (example-session "Usage"
    (m/degrees (relative-angle-between (vec3 1.0 0.0 0.0) (vec3 0.0 1.0 0.0)))
    (m/degrees (relative-angle-between (vec (repeatedly 50 rand)) (vec (repeatedly 50 rand))))))

(add-examples aligned?
  (example-session "Usage" {:test-values [true false true]}
    (aligned? (vec2 1.0 1.0) (vec2 2.0 2.000001))
    (aligned? (vec2 1.0 1.0) (vec2 2.0 2.00001))
    (aligned? (vec2 1.0 1.0) (vec2 2.0 2.00001) 0.001)))

(add-examples faceforward
  (example-session "Usage"
    (faceforward (vec2 1.0 1.0) (vec2 -1.0 -3.0))
    (faceforward (vec2 1.0 1.0) (vec2 1.0 0.0))))

(add-examples generate-vec2 (example-session "Usage" (generate-vec2 (constantly 2)) (generate-vec2 rand (constantly 1))))
(add-examples generate-vec3 (example-session "Usage" (generate-vec3 rand) (generate-vec3 rand (constantly 1) (constantly 2))))
(add-examples generate-vec4 (example-session "Usage" (generate-vec4 rand) (generate-vec4 rand rand (constantly 1) (constantly 2))))

(add-examples array->vec2 (example (array->vec2 (double-array [11 22 33 44 55]))))
(add-examples array->vec3 (example (array->vec3 (double-array [11 22 33 44 55]))))
(add-examples array->vec4 (example (array->vec4 (double-array [11 22 33 44 55]))))

(add-examples seq->vec2 (example (seq->vec2 [11 22 33 44 55])))
(add-examples seq->vec3 (example (seq->vec3 (lazy-seq '(1 2)))))
(add-examples seq->vec4 (example (seq->vec4 (double-array [11 22 33 44 55]))))

(add-examples mag (example "Length of the vector" {:test-value (m/sqrt 2.0)} (mag (vec2 1 1))))
(add-examples magsq (example "Length of the vector, squared" {:test-value 10.0} (magsq [1 2 1 2])))
(add-examples abs (example "Usage" (abs (array-vec [-1 2 -2 1]))))
(add-examples add (example "Usage" (add (vec2 1 2) (vec2 -1 -2))))
(add-examples applyf (example "Usage" (applyf [5 3] (fn [v] (m/sin v)))))
(add-examples fmap (example "Usage" (fmap [5 3] (fn [v] (m/sin v)))))
(add-examples approx (example-session "Usage" (approx (vec2 (m/sin -1) (m/cos 1))) (approx (vec2 (m/sin -1) (m/cos 1)) 6)))
(add-examples div
  (example "Usage" (div [5 4 3 5] 4.0))
  (example "Applied to vector is reciprocal" (div [1 2 4 -5])))
(add-examples mult (example "Usage" (mult (vec4 5 4 3 5) 4.0)))
(add-examples dot (example "Usage" (dot (vec4 1 0 0 0) (vec4 -1 0 -1 0))))
(add-examples econstrain (example "Usage" (econstrain (vec3 2 0 -2) -1 1)))
(add-examples emn (example "Usage" (emn [-1 2] [1 -2])))
(add-examples emx (example "Usage" (emx [-1 2] [1 -2])))
(add-examples emult (example "Usage" (emult (vec3 1 2 3) (vec3 9 9 9))))
(add-examples is-near-zero?
  (example-session "Usage" {:test-values [true false true]}
    (is-near-zero? [0 0.0000001])
    (is-near-zero? [0 0.000001])
    (is-near-zero? [0 0.000001] 0.0001)))
(add-examples is-zero? (example-session "Usage" (is-zero? [0 0.0000001]) (is-zero? [0 0.0])))
(add-examples mn (example "Usage" (mn (vec4 -1 -2 3 4))))
(add-examples mx (example "Usage" (mx (vec4 -1 -2 3 4))))
(add-examples permute (example "Usage" (permute (vec4 1 2 3 4) (vec4 0 3 2 1))))
(add-examples reciprocal (example "Usage" (reciprocal (vec3 1 2 5))))
(add-examples sub (example "Usage" (sub (vec2 1 2) (vec2 -1 -2))))
(add-examples sum (example "Usage" (sum (vec (range 1000)))))

(add-examples transform
  (example-session "Usage"
    (transform (vec2 1 1) (vec2 -1 -1) (vec2 1.0 0.0) (vec2 0.0 1.0))
    (transform (vec3 1 1 1) (vec3 -1 -1 0) (vec3 1.0 0.0 1.0) (vec3 0.0 1.0 0.0) (vec3 0.0 1.0 1.0))))

(add-examples rotate
  (example-session "Usage"
    (rotate (vec2 1 2) (m/radians 90))
    (rotate (vec3 1 2 3) (m/radians 90) 0 0)))

(add-examples perpendicular
  (example-session "Usage"
    (perpendicular (vec2 -4 0))
    (perpendicular (vec3 1.0 0.0 0.0) (vec3 0.0 -1.0 0.0))))

(add-examples maxdim
  (example-session "Usage"
    (let [v (vec (repeatedly 100 (fn [] (- (int (rand-int 200)) 100))))
          mdim (maxdim v)]
      [mdim (v mdim)])
    (maxdim (vec3 1 2 3))))

(add-examples mindim
  (example-session "Usage"
    (let [v (vec (repeatedly 100 (fn [] (- (int (rand-int 200)) 100))))
          mdim (mindim v)]
      [mdim (v mdim)])
    (mindim (vec3 1 2 3))))

(add-examples heading
  (example-session "Usage"
    (m/degrees (heading (vec2 1.0 1.0)))
    (m/degrees (heading (vec3 1.0 0.0 1.0)))
    (m/degrees (heading (vec4 1.0 -1.0 1.0 -1.0)))
    (heading [1 2 3 4 5 6 7 8 9])
    (heading (array-vec [1 2 3 4 5 6 7 8 9]))))

(add-examples from-polar
  (example-session "Usage"
    (from-polar (vec2 1.0 (m/radians 90)))
    (from-polar (vec3 1.0 (m/radians 90) (m/radians 45)))))

(add-examples to-polar
  (example-session "Usage"
    (to-polar (vec2 1.0 1.0))
    (to-polar (vec3 1.0 0.0 1.0))))

(add-examples interpolate
  (example-session "Usage"
    (interpolate (vec2 -1 -1) (vec2 1 1) 0.25)
    (interpolate (vec2 -1 -1) (vec2 1 1) 0.25 m/smooth-interpolation)))

(add-examples einterpolate
  (example-session "Usage"
    (einterpolate (vec2 -1 -1) (vec2 1 1) (vec2 0.25 0.75))
    (einterpolate (vec2 -1 -1) (vec2 1 1) (vec2 0.25 0.75) m/smooth-interpolation)))

(add-examples axis-rotate
  (example-session "Usage"
    (axis-rotate (vec3 1.0 0.0 0.0) m/HALF_PI (vec3 0.0 1.0 0.0))
    (axis-rotate (vec3 1.0 0.0 0.0) m/HALF_PI (vec3 0.0 1.0 0.0) (vec3 1.0 1.0 1.0))))

(add-examples base-from
  (example-session "Usage"
    (base-from (vec3 0.1 1.0 -1.0))
    (base-from (vec2 1.0 0.0))))

(add-examples cross
  (example-session "Usage"
    (cross (vec3 1 2 3)
           (vec3 4 3 2))
    (cross (vec2 1 2)
           (vec2 -1 2))))

(add-examples vec->Vec
  (example-session "Check types"
    (type (vec->Vec [1 2 3]))
    (type (vec->Vec (vec2 1 2)))
    (type (vec->Vec (vec3 1 2 3)))
    (type (vec->Vec (vec4 1 2 3 4)))
    (type (vec->Vec (array-vec 1)))))

(add-examples as-vec
  (example-session "Create vectors"
    (as-vec [1 2 3])
    (as-vec [3 2 3] (repeat -10))
    (as-vec (vec3 1 2 3) (repeat -10))
    (type (as-vec (array-vec [1 1 11 2 3]) (repeat -10)))))

(add-examples clamp
  (example-session "Clamp"
    (clamp [-1 -2 -3 1 2 3])
    (clamp [-1 -2 -3 1 2 3] -2.5 2.5)))

(add-examples make-vector
  (example-session "Usage"
    (make-vector 2)
    (make-vector 2 (cycle [1 2 3]))
    (make-vector 3)
    (make-vector 3 (cycle [1 2 3]))
    (make-vector 4)
    (make-vector 4 (cycle [1 2 3]))
    (make-vector 6)
    (make-vector 6 (cycle [1 2 3]))
    (type (make-vector 6))))

(add-examples nonzero-count
  (example (nonzero-count [1 2 3 -1 0 2 0 0])))

(add-examples zero-count
  (example (zero-count [1 2 3 -1 0 2 0 0])))

;;

(defmacro math-ops
  [ops]
  `(do ~@(for [op ops]
           `(add-examples ~op (example (~op (vec4 0.5 -1.5 2.1 0.0)))))))

(math-ops [sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh
           cot sec csc acot asec acsc coth sech csch acoth asech csch
           sq safe-sqrt sqrt cbrt exp expm1 log log10 log2 ln log1p
           radians degrees sinc sigmoid
           floor ceil round rint trunc frac sfrac signum sgn])
