;; https://github.com/MartinWeigel/Quaternion/blob/master/TestQuaternion.c
(ns fastmath.quaternion-test
  (:require [fastmath.quaternion :as sut]
            [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.complex :as c]
            [fastmath.matrix :as mat]
            [fastmath.vector :as v]
            [fastmath.core :as m]))

(t/deftest rotation-quaternion
  (t/is (v/delta-eq (sut/rotation-quaternion m/HALF_PI (v/vec3 1 0 0))
                    (sut/quaternion 0.7071 0.7071 0 0) 1.0e-4)))

(t/deftest multiplication
  (t/are [a2 b2 c2 d2 ra rb rc rd]
      (v/delta-eq (sut/mult (sut/quaternion 1 0 1 0)
                            (sut/quaternion a2 b2 c2 d2))
                  (sut/quaternion ra rb rc rd) 1.0e-4)
    1 0.5 0.5 0.75 0.5 1.25 1.5 0.25
    1 0 1 0 0 0 2 0
    2 1 0.1 0.1 1.9 1.1 2.1 -0.9))

(t/deftest rotate
  (t/are [x y z a b c d rx ry rz]
      (v/delta-eq (v/vec3 rx ry rz)
                  (sut/rotate (v/vec3 x y z) (sut/quaternion a b c d)) 1.0e-4)
    5.1 6.8 -5.3 1 0 0 0 5.1 6.8 -5.3
    1 0 0 0.5 0.5 0.5 0.5 0 1 0
    1 0 0 0.6532815, -0.270598, 0.270598, 0.6532815 0 0.7071 -0.7071
    1 0 0 0.7071, 0, 0, 0.7071 0 1 0))

(t/deftest euler

  (t/are [r p y a b c d] (v/delta-eq (sut/from-euler (v/fmap (v/vec3 r p y) m/radians))
                                     (sut/quaternion a b c d) 1.0e-4)
    165 63 122 0.5070333, 0.3501829, 0.7724199, -0.1538071
    90 0 90 0.5 0.5 0.5 0.5
    90 0 180 0 0 0.7071 0.7071
    90 0 0 0.7071 0.7071 0 0)
  
  (t/are [r p y a b c d] (v/delta-eq (sut/to-euler (sut/quaternion a b c d))
                                     (v/fmap (v/vec3 r p y) m/radians) 1.0e-4)
    165 63 122 0.5070333, 0.3501829, 0.7724199, -0.1538071
    90 0 90 0.5 0.5 0.5 0.5
    90 0 180 0 0 0.7071 0.7071
    90 0 0 0.7071 0.7071 0 0))

(t/deftest slerp
  (let [q1 (sut/quaternion 0.6532815, -0.270598, 0.270598, 0.6532815)
        q2 (sut/quaternion 0.5 0.5 0.5 0.5)]
    (t/is (= q1 (sut/slerp q1 q2 0.0)))
    (t/is (= q2 (sut/slerp q1 q2 1.0)))
    (t/is (v/delta-eq (sut/slerp q1 q2 0.62)
                      (sut/quaternion 0.6119266025696755
                                      0.22069444274723088
                                      0.4498729015909088
                                      0.6119266025696755) 1.0e-4))))

;; ---------------------------------------------------------------------------
;; audit 2026-10-27: reference data (Quaternions.jl, scipy, pyquaternion), contracts and edge cases

(def ^:private reference (delay (edn/read-string (slurp (io/resource "quaternion/quaternion_reference.edn")))))
(def ^:private rotation-reference (delay (edn/read-string (slurp (io/resource "quaternion/quaternion_rotation_reference.edn")))))

(defn- parts
  [^fastmath.vector.Vec4 q]
  [(.x q) (.y q) (.z q) (.w q)])

(defn- parts3
  [^fastmath.vector.Vec3 x]
  [(.x x) (.y x) (.z x)])

(defn- max-diff
  ^double [a b]
  (reduce m/max 0.0 (map (fn [x y] (m/abs (- (double x) (double y)))) a b)))

(defn- rotation-error
  "Distance of two quaternions up to the sign, q and -q are the same rotation."
  ^double [a b]
  (m/min (max-diff a b) (max-diff a (map - b))))

(defn- normwise-error
  "Error of the quaternion `q` against the reference parts, relative to the magnitude of the reference, infinite parts must be equal."
  ^double [expected q]
  (let [r (parts q)]
    (cond
      (some m/nan? r) ##Inf
      (some true? (map (fn [a b] (or (and (m/inf? a) (not= a b)) (and (m/inf? b) (not (m/inf? a))))) expected r)) ##Inf
      :else (let [fe (map #(if (m/inf? %) 0.0 %) expected)
                  d (map (fn [a b] (if (m/inf? a) 0.0 (- b a))) expected r)
                  h (fn [xs] (m/sqrt (reduce + (map #(* % %) xs))))]
              (/ (h d) (m/max (h fe) 1.0e-320))))))

(def ^:private functions
  {:exp sut/exp :log sut/log :sqrt sut/sqrt
   :sin sut/sin :cos sut/cos :tan sut/tan :sinh sut/sinh :cosh sut/cosh :tanh sut/tanh
   :asin sut/asin :acos sut/acos :atan sut/atan :asinh sut/asinh :acosh sut/acosh :atanh sut/atanh
   :sec sut/sec :csc sut/csc :cot sut/cot :sech sut/sech :csch sut/csch :coth sut/coth})

(t/deftest reference-functions
  ;; Quaternions.jl, 257 points with a non zero vector part (including 1e-9 and 1e5 sized vectors):
  ;; the maximum normwise error is 1.2e-15 (coth), 9e-16 (tanh, sec), below 7e-16 for the others
  (let [rf @reference
        points (:points rf)]
    (doseq [[f g] functions
            :let [bad (filter (fn [[p e]] (> (normwise-error e (g (apply sut/quaternion p))) 5.0e-15))
                              (map vector points (get rf f)))]]
      (t/is (empty? bad) (str f ": " (count bad) " points above 5e-15, the first " (first bad))))))

(t/deftest reference-arithmetic
  (let [rf @reference
        args (:pair-args rf)]
    (doseq [[k g tol] [[:mult sut/mult 1.0e-15] [:div sut/div 3.0e-15] [:pow sut/pow 1.0e-14]]
            :let [bad (filter (fn [[[a b] e]] (> (normwise-error e (g (apply sut/quaternion a) (apply sut/quaternion b))) tol))
                              (map vector args (get rf k)))]]
      (t/is (empty? bad) (str k ": " (count bad) " pairs above " tol ", the first " (first bad))))
    (let [bad (filter (fn [[p e]] (> (normwise-error e (sut/reciprocal (apply sut/quaternion p))) 3.0e-15))
                      (map vector (:points rf) (:inv rf)))]
      (t/is (empty? bad) "reciprocal against inv"))
    (let [bad (filter (fn [[p e]] (> (normwise-error e (sut/conjugate (apply sut/quaternion p))) 0.0))
                      (map vector (:points rf) (:conj rf)))]
      (t/is (empty? bad) "conjugate"))))

(t/deftest derived-functions
  ;; secant, ... , inverse secant ... have no counterparts in Quaternions.jl: identities in the plane of the quaternion
  (doseq [p (take-nth 11 (:points @reference))
          :let [q (apply sut/quaternion p)
                nq (sut/normalize q)
                small (sut/scale nq 0.4)]]
    (t/is (sut/delta-eq (sut/sec q) (sut/reciprocal (sut/cos q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/csc q) (sut/reciprocal (sut/sin q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/cot q) (sut/reciprocal (sut/tan q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/sech q) (sut/reciprocal (sut/cosh q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/csch q) (sut/reciprocal (sut/sinh q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/coth q) (sut/reciprocal (sut/tanh q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/asec q) (sut/acos (sut/reciprocal q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/acsc q) (sut/asin (sut/reciprocal q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/acot q) (sut/atan (sut/reciprocal q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/asech q) (sut/acosh (sut/reciprocal q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/acsch q) (sut/asinh (sut/reciprocal q)) 1.0e-12))
    (t/is (sut/delta-eq (sut/acoth q) (sut/atanh (sut/reciprocal q)) 1.0e-12))
    ;; the principal values invert the functions in the neighbourhood of zero
    (t/is (sut/delta-eq small (sut/asin (sut/sin small)) 1.0e-13))
    (t/is (sut/delta-eq small (sut/atan (sut/tan small)) 1.0e-13))
    (t/is (sut/delta-eq small (sut/asinh (sut/sinh small)) 1.0e-13))
    (t/is (sut/delta-eq small (sut/atanh (sut/tanh small)) 1.0e-13))
    (t/is (sut/delta-eq (if (m/neg? (sut/re small)) (sut/neg small) small) (sut/acos (sut/cos small)) 1.0e-13)
          "acos(cos q) = q for the scalar part in [0, pi], -q for a negative one")
    (when (< (sut/norm q) 10.0)
      (t/is (sut/delta-eq sut/ONE (sut/add (sut/sq (sut/sin q)) (sut/sq (sut/cos q))) 1.0e-9) "sin^2 + cos^2 = 1 in the plane of q"))
    (t/is (sut/delta-eq q (sut/exp (sut/log q)) (* 1.0e-12 (m/max 1.0 (sut/norm q)))))
    (t/is (sut/delta-eq q (sut/sq (sut/sqrt q)) (* 1.0e-12 (m/max 1.0 (sut/norm q)))))))

(t/deftest generated-function-range
  (t/testing "exp does not give NaN for a large exponent"
    (t/is (= (v/vec4 ##Inf ##Inf 0.0 0.0) (sut/exp (sut/quaternion 1000.0 1.0 0.0 0.0))) "was [Inf Inf NaN NaN]")
    (t/is (not (sut/nan? (sut/exp (sut/quaternion 1000.0 1.0 2.0 3.0)))))
    (t/is (not (sut/nan? (sut/sin (sut/quaternion 1.0 1.0e200 1.0e200 0.0)))) "huge vector part"))
  (t/testing "the magnitude of the vector part does not overflow"
    (let [r (sut/log (sut/quaternion 1.0 1.0e200 1.0e200 1.0e200))]
      (t/is (m/delta-eq (+ (* 200.0 (m/log 10.0)) (* 0.5 (m/log 3.0))) (.x r) 1.0e-12))
      (t/is (m/delta-eq m/HALF_PI (m/hypot (.y r) (.z r) (.w r)) 1.0e-12))))
  (t/testing "a zero vector part: the i axis"
    (t/is (= (v/vec4 0.0 2.0 0.0 0.0) (sut/sqrt (sut/quaternion -4.0 0.0 0.0 0.0))))
    (t/is (= (v/vec4 (m/log 2.0) m/PI 0.0 0.0) (sut/log (sut/quaternion -2.0 0.0 0.0 0.0))))
    (t/is (= (v/vec4 2.0 0.0 0.0 0.0) (sut/sqrt (sut/quaternion 4.0 0.0 0.0 0.0))))
    (t/is (= (v/vec4 ##-Inf 0.0 0.0 0.0) (sut/log sut/ZERO)))
    (t/is (= sut/ONE (sut/exp sut/ZERO))))
  (t/testing "docstrings"
    (t/is (every? #(> (count (:doc (meta (ns-resolve 'fastmath.quaternion (symbol (name %)))))) 100) (keys functions)))))

(t/deftest arithmetic-contracts
  (let [a (sut/quaternion 1.0 2.0 3.0 4.0)
        b (sut/quaternion -0.5 0.25 2.0 -1.5)
        d (sut/quaternion 0.3 -1.0 0.0 2.0)]
    (t/testing "multiplication table"
      (t/is (= sut/K (sut/mult sut/I sut/J)))
      (t/is (= sut/I (sut/mult sut/J sut/K)))
      (t/is (= sut/J (sut/mult sut/K sut/I)))
      (t/is (= sut/-K (sut/mult sut/J sut/I)))
      (t/is (= (sut/neg sut/ONE) (sut/mult sut/I sut/I)))
      (t/is (= (sut/neg sut/ONE) (sut/mult (sut/mult sut/I sut/J) sut/K)))
      (t/is (not= (sut/mult a b) (sut/mult b a)) "not commutative"))
    (t/testing "variadic"
      (t/is (= sut/ZERO (sut/add)))
      (t/is (= sut/ONE (sut/mult)))
      (t/is (= a (sut/add a)))
      (t/is (= a (sut/mult a)))
      (t/is (= (sut/neg a) (sut/sub a)))
      (t/is (= (sut/add (sut/add a b) d) (sut/add a b d)))
      (t/is (= (sut/sub (sut/sub a b) d) (sut/sub a b d)))
      (t/is (= (sut/mult (sut/mult a b) d) (sut/mult a b d)) "left to right")
      (t/is (not= (sut/mult a b d) (sut/mult d b a)))
      (t/is (= (sut/mult (sut/mult (sut/mult a b) d) a) (sut/mult a b d a))))
    (t/testing "numbers and complex numbers"
      (t/is (= (sut/scale a 2.0) (sut/mult a 2.0)))
      (t/is (= (sut/scale a 2.0) (sut/mult 2.0 a)))
      (t/is (= (sut/add a (sut/quaternion 3.0)) (sut/add a 3.0)))
      (t/is (= (sut/adds a 3.0) (sut/add a 3.0)))
      (t/is (= (sut/sub a (sut/quaternion 3.0 1.0 0.0 0.0)) (sut/sub a (c/complex 3.0 1.0))))
      (t/is (= (sut/mult a sut/I) (sut/mult a (c/complex 0.0 1.0))))
      (t/is (sut/delta-eq (sut/pow a (sut/quaternion 2.0)) (sut/pow a 2.0) 1.0e-15))
      (t/is (sut/delta-eq (sut/sq a) (sut/pow a 2.0) 1.0e-12))
      (t/is (sut/delta-eq (sut/div a (sut/quaternion 2.0)) (sut/div a 2.0) 1.0e-15))
      (t/is (sut/delta-eq (sut/logb a (sut/quaternion 2.0)) (sut/logb a 2.0) 1.0e-15))
      (t/is (= 5.0 (sut/dot (sut/quaternion 1.0 2.0 0.0 0.0) (sut/quaternion 1.0 2.0 0.0 0.0))))
      (t/is (= 3.0 (sut/dot 3.0 1.0)))
      (t/is (= (sut/mult 2.0 3.0) (sut/quaternion 6.0)))
      (t/is (= a (sut/ensure-quaternion a))))
    (t/testing "division and reciprocal"
      (t/is (sut/delta-eq sut/ONE (sut/mult a (sut/reciprocal a)) 1.0e-15))
      (t/is (sut/delta-eq sut/ONE (sut/mult (sut/reciprocal a) a) 1.0e-15))
      (t/is (= (sut/reciprocal a) (sut/inverse a)))
      (t/is (sut/delta-eq a (sut/mult (sut/div a b) b) 1.0e-14))
      (t/is (sut/delta-eq (sut/conjugate (sut/normalize a)) (sut/inverse (sut/normalize a)) 1.0e-15) "unit quaternion: the inverse is the conjugate")
      (t/is (sut/nan? (sut/reciprocal sut/ZERO)))
      (t/is (sut/zero? (sut/reciprocal (sut/quaternion ##Inf 1.0 2.0 3.0))) "the reciprocal of an infinite quaternion is zero")
      (doseq [s [1.0e-200 1.0e200 1.0e-300]
              :let [q (sut/scale a s)]]
        (t/is (sut/delta-eq sut/ONE (sut/mult q (sut/reciprocal q)) 1.0e-14) (str "q * q^-1 = 1 for the scale " s))))
    (t/testing "norm"
      (t/is (= 5.0 (sut/norm (sut/quaternion 3.0 4.0 0.0 0.0))))
      (t/is (= 2.0e200 (sut/norm (sut/quaternion 1.0e200 1.0e200 1.0e200 1.0e200))) "no overflow")
      (t/is (= 2.0e-200 (sut/norm (sut/quaternion 1.0e-200 1.0e-200 1.0e-200 1.0e-200))) "no underflow")
      (t/is (= (m/sqrt 30.0) (sut/norm a)))
      (t/is (m/delta-eq 1.0 (sut/norm (sut/normalize a)) 1.0e-15))
      (t/is (= sut/ZERO (sut/normalize sut/ZERO)))
      (t/is (= 5.0 (sut/norm (c/complex 3.0 4.0))) "a complex number as a quaternion")
      (t/is (= 3.0 (+ 1.0 (sut/norm (sut/quaternion 2.0))))))
    (t/testing "parts and predicates"
      (t/is (= 1.0 (sut/scalar a) (sut/re a)))
      (t/is (= (v/vec3 2.0 3.0 4.0) (sut/vector a)))
      (t/is (= [2.0 3.0 4.0] [(sut/im-i a) (sut/im-j a) (sut/im-k a)]))
      (t/is (sut/real? (sut/quaternion 3.0)))
      (t/is (not (sut/real? a)))
      (t/is (sut/imaginary? (sut/quaternion 0.0 1.0 2.0 3.0)))
      (t/is (sut/imaginary? sut/ZERO))
      (t/is (sut/zero? sut/ZERO))
      (t/is (not (sut/zero? sut/I)))
      (t/is (sut/inf? (sut/quaternion 1.0 ##Inf 0.0 0.0)))
      (t/is (sut/nan? (sut/quaternion 1.0 ##NaN 0.0 0.0)))
      (t/is (sut/invalid? (sut/quaternion 1.0 ##NaN 0.0 0.0)))
      (t/is (sut/valid? a))
      (t/is (not (sut/valid? (sut/quaternion ##Inf 1.0 0.0 0.0))))
      (t/is (= sut/I (sut/complex->quaternion (c/complex 0.0 1.0))))
      (t/is (= a (sut/quaternion 1.0 [2.0 3.0 4.0])))
      (t/is (= (sut/quaternion 5.0) (sut/quaternion 5.0 0.0 0.0 0.0)))
      (t/is (= (v/vec4 1.0 -2.0 -3.0 -4.0) (sut/conjugate a)))
      (t/is (= 1.0 (sut/qsgn a)))
      (t/is (= -1.0 (sut/qsgn (sut/neg a))))
      (t/is (= 0.0 (sut/qsgn sut/ZERO)))
      (t/is (= 0.0 (sut/qsgn 0.0 0.0 0.0 0.0)))
      (t/is (= 1.0 (sut/qsgn 2.0 -1.0 0.0 0.0)))
      (t/is (m/delta-eq (m/atan2 (m/sqrt 29.0) 1.0) (sut/arg a) 1.0e-15))
      (t/is (= m/PI (sut/arg (sut/quaternion -1.0 0.0 0.0 0.0)))))))

(t/deftest rotation-reference-conversions
  ;; scipy: 60 random rotations, away from the gimbal lock, maximum errors 2.7e-15 (to-euler), 8.9e-16 (to-angles), the others below 6e-16
  (let [rr @rotation-reference]
    (doseq [[qq e a mt aa] (map vector (:q rr) (:euler rr) (:angles rr) (:matrix rr) (:axis-angle rr))
            :let [q (apply sut/quaternion qq)]]
      (t/is (< (rotation-error qq (parts (sut/from-euler (nth e 0) (nth e 1) (nth e 2)))) 1.0e-14) "from-euler")
      (t/is (< (max-diff e (parts3 (sut/to-euler q))) 1.0e-13) "to-euler")
      (t/is (< (rotation-error qq (parts (sut/from-euler (v/vec3 (nth e 0) (nth e 1) (nth e 2))))) 1.0e-14) "from-euler, vector")
      (t/is (< (rotation-error qq (parts (sut/from-angles (nth a 0) (nth a 1) (nth a 2)))) 1.0e-14) "from-angles")
      (t/is (< (max-diff a (parts3 (sut/to-angles q))) 1.0e-13) "to-angles")
      (t/is (< (max-diff mt (let [^fastmath.matrix.Mat3x3 mm (sut/to-rotation-matrix q)]
                              [(.a00 mm) (.a01 mm) (.a02 mm) (.a10 mm) (.a11 mm) (.a12 mm) (.a20 mm) (.a21 mm) (.a22 mm)]))
               1.0e-14) "to-rotation-matrix")
      (t/is (< (rotation-error qq (parts (sut/from-rotation-matrix (apply mat/mat3x3 mt)))) 1.0e-14) "from-rotation-matrix")
      (t/is (< (rotation-error qq (parts (sut/from-axis-angle (first aa) (apply v/vec3 (rest aa))))) 1.0e-14) "from-axis-angle")
      (t/is (< (rotation-error qq (parts (sut/rotation-quaternion (first aa) (apply v/vec3 (rest aa))))) 1.0e-14) "rotation-quaternion")
      (let [r (sut/to-axis-angle q)]
        (t/is (< (max-diff aa (cons (:angle r) (parts3 (:axis r)))) 1.0e-13) "to-axis-angle")))))

(t/deftest rotation-consistency
  (let [rr @rotation-reference]
    (doseq [qq (take 20 (:q rr))
            :let [q (apply sut/quaternion qq)
                  mm (sut/to-rotation-matrix q)
                  x [0.3 -1.2 2.5]
                  by-matrix (v/vec3 (+ (* (.a00 mm) 0.3) (* (.a01 mm) -1.2) (* (.a02 mm) 2.5))
                                    (+ (* (.a10 mm) 0.3) (* (.a11 mm) -1.2) (* (.a12 mm) 2.5))
                                    (+ (* (.a20 mm) 0.3) (* (.a21 mm) -1.2) (* (.a22 mm) 2.5)))]]
      (t/is (v/delta-eq by-matrix (sut/rotate x q) 1.0e-13) "rotate = matrix times vector")
      (t/is (v/delta-eq by-matrix (sut/rotate x (sut/scale q 7.0)) 1.0e-13) "a non unit rotation quaternion is normalized")
      (t/is (v/delta-eq by-matrix (sut/rotate x (sut/neg q)) 1.0e-13) "q and -q")
      (t/is (v/delta-eq (v/vec3 0.3 -1.2 2.5) (sut/rotate by-matrix (sut/inverse q)) 1.0e-13) "the inverse rotates back"))
    (t/testing "the angles functions are the rotation matrix of the matrix namespace"
      (doseq [[x y z] [[0.3 -0.5 1.2] [1.0 0.2 -2.0] [-2.5 1.4 0.1]]
              :let [^fastmath.matrix.Mat3x3 a (sut/to-rotation-matrix (sut/from-angles x y z))
                    ^fastmath.matrix.Mat3x3 b (mat/rotation-matrix-3d x y z)]]
        (t/is (< (max-diff [(.a00 a) (.a01 a) (.a02 a) (.a10 a) (.a11 a) (.a12 a) (.a20 a) (.a21 a) (.a22 a)]
                           [(.a00 b) (.a01 b) (.a02 b) (.a10 b) (.a11 b) (.a12 b) (.a20 b) (.a21 b) (.a22 b)]) 1.0e-14))))
    (t/testing "rotate around the axis"
      (t/is (v/delta-eq (v/vec3 0.0 1.0 0.0) (sut/rotate (v/vec3 1.0 0.0 0.0) m/HALF_PI (v/vec3 0.0 0.0 1.0)) 1.0e-15))
      (t/is (v/delta-eq (v/axis-rotate (v/vec3 1.0 2.0 3.0) 0.7 (v/vec3 0.3 -0.2 1.0))
                        (sut/rotate (v/vec3 1.0 2.0 3.0) 0.7 (v/vec3 0.3 -0.2 1.0)) 1.0e-13))
      (t/is (v/delta-eq (v/vec3 1.0 2.0 3.0) (sut/rotate [1 2 3] sut/ZERO) 1.0e-15) "the zero quaternion does not rotate")
      (t/is (v/delta-eq (v/vec3 1.0 2.0 3.0) (sut/rotate [1 2 3] sut/ONE) 1.0e-15)))))

(t/deftest rotation-quaternion-contract
  (t/is (= 1.0 (sut/norm (sut/rotation-quaternion 1.234 (v/vec3 3.0 -4.0 12.0)))))
  (t/is (m/delta-eq 1.0 (sut/norm (sut/rotation-quaternion 1.234 (v/vec3 3.0 -4.0 12.0))) 1.0e-15))
  (t/is (= sut/ONE (sut/rotation-quaternion 0.0 (v/vec3 0.0 0.0 5.0))))
  (t/is (sut/delta-eq (sut/neg sut/ONE) (sut/rotation-quaternion m/TWO_PI (v/vec3 0.0 0.0 5.0)) 1.0e-15) "a full turn is -1")
  (t/is (sut/delta-eq sut/K (sut/rotation-quaternion m/PI (v/vec3 0.0 0.0 5.0)) 1.0e-15))
  (t/is (= (sut/from-axis-angle 0.4 (v/vec3 0.0 1.0 0.0)) (sut/rotation-quaternion 0.4 (v/vec3 0.0 1.0 0.0))))
  (t/testing "the zero axis is rejected"
    (let [e (try (sut/rotation-quaternion 1.0 (v/vec3 0.0 0.0 0.0)) nil (catch clojure.lang.ExceptionInfo e e))]
      (t/is (some? e))
      (t/is (= :u (:option (ex-data e)))))
    (t/is (thrown? clojure.lang.ExceptionInfo (sut/from-axis-angle 1.0 (v/vec3 0.0 0.0 0.0))))))

(t/deftest axis-angle
  (let [axis (v/normalize (v/vec3 1.0 -2.0 0.5))]
    (t/testing "to-axis-angle"
      (let [r (sut/to-axis-angle (sut/rotation-quaternion 2.0 axis))]
        (t/is (m/delta-eq 2.0 (:angle r) 1.0e-14))
        (t/is (v/delta-eq axis (:axis r) 1.0e-14)))
      (let [r (sut/to-axis-angle (sut/neg (sut/rotation-quaternion 2.0 axis)))]
        (t/is (m/delta-eq 2.0 (:angle r) 1.0e-14) "-q is the same rotation")
        (t/is (v/delta-eq axis (:axis r) 1.0e-14)))
      (let [r (sut/to-axis-angle (sut/rotation-quaternion 4.0 axis))]
        (t/is (m/delta-eq (- m/TWO_PI 4.0) (:angle r) 1.0e-14) "the shortest angle")
        (t/is (v/delta-eq (v/mult axis -1.0) (:axis r) 1.0e-14)))
      (t/is (m/delta-eq 1.0 (:angle (sut/to-axis-angle (sut/scale (sut/rotation-quaternion 1.0 axis) 5.0))) 1.0e-14) "normalized")
      (t/is (= {:axis (v/vec3 1.0 0.0 0.0) :angle 0.0} (sut/to-axis-angle sut/ONE)) "no rotation")
      (t/is (= {:axis (v/vec3 1.0 0.0 0.0) :angle 0.0} (sut/to-axis-angle (sut/neg sut/ONE))))
      (t/is (m/delta-eq m/PI (:angle (sut/to-axis-angle sut/I)) 1.0e-15)))))

(t/deftest from-two-vectors
  (doseq [[u w] [[[1.0 0.0 0.0] [0.0 1.0 0.0]]
                 [[1.0 2.0 3.0] [-3.0 0.5 2.0]]
                 [[0.0 0.0 5.0] [0.0 0.0 2.0]]
                 [[1.0 0.0 0.0] [-1.0 0.0 0.0]]
                 [[0.2 0.5 -0.1] [-0.2 -0.5 0.1]]
                 [[0.0 1.0 0.0] [0.0 -1.0 0.0]]
                 [[1.0 1.0 1.0] [-1.0 -1.0 -1.0000000001]]
                 [[1.0e-100 2.0e-100 0.0] [3.0e100 -1.0e100 4.0e100]]]
          :let [u (apply v/vec3 u)
                w (apply v/vec3 w)
                q (sut/from-two-vectors u w)]]
    (t/is (m/delta-eq 1.0 (sut/norm q) 1.0e-14) "unit quaternion")
    (t/is (v/delta-eq (v/normalize w) (v/normalize (sut/rotate u q)) 1.0e-7)
          (str "rotates the direction of " u " to " w))
    (t/is (v/delta-eq (v/normalize w) (sut/rotate (v/normalize u) q) 1.0e-7)))
  (t/testing "the shortest rotation"
    (let [q (sut/from-two-vectors (v/vec3 1.0 0.0 0.0) (v/vec3 0.0 1.0 0.0))]
      (t/is (sut/delta-eq (sut/rotation-quaternion m/HALF_PI (v/vec3 0.0 0.0 1.0)) q 1.0e-15)))
    (t/is (= sut/ONE (sut/from-two-vectors (v/vec3 1.0 0.0 0.0) (v/vec3 5.0 0.0 0.0))) "parallel: no rotation")
    (t/is (m/delta-eq m/PI (:angle (sut/to-axis-angle (sut/from-two-vectors (v/vec3 1.0 0.0 0.0) (v/vec3 -2.0 0.0 0.0)))) 1.0e-12)
          "opposite: a half turn"))
  (t/testing "zero vectors"
    (t/is (= :u (:option (ex-data (try (sut/from-two-vectors (v/vec3 0.0 0.0 0.0) (v/vec3 1.0 0.0 0.0)) (catch clojure.lang.ExceptionInfo e e))))))
    (t/is (= :v (:option (ex-data (try (sut/from-two-vectors (v/vec3 1.0 0.0 0.0) (v/vec3 0.0 0.0 0.0)) (catch clojure.lang.ExceptionInfo e e))))))))

(t/deftest slerp-contract
  (let [rr @rotation-reference
        z (v/vec3 0.0 0.0 1.0)]
    (t/testing "pyquaternion: the same rotation (the sign follows the other convention for a negative dot product)"
      (doseq [[[a b t] e] (map vector (:slerp-args rr) (:slerp rr))
              :let [r (parts (sut/slerp (apply sut/quaternion a) (apply sut/quaternion b) t))]]
        (t/is (< (rotation-error e r) 1.0e-13) (str "slerp " t))
        (t/is (>= (reduce + (map * a r)) 0.0) "the path starts at q1 and goes the short way")))
    (t/testing "endpoints are returned unchanged"
      (let [a (sut/rotation-quaternion 0.3 z)
            b (sut/rotation-quaternion 1.9 z)]
        (t/is (= a (sut/slerp a b 0.0)))
        (t/is (= b (sut/slerp a b 1.0)))
        (t/is (= a (sut/slerp a b -5.0)) "t is constrained")
        (t/is (= b (sut/slerp a b 7.0)))))
    (t/testing "constant angular velocity"
      (let [a (sut/rotation-quaternion 0.3 z)
            b (sut/rotation-quaternion 1.9 z)]
        (doseq [t [0.1 0.25 0.5 0.8]]
          (t/is (sut/delta-eq (sut/rotation-quaternion (+ 0.3 (* t 1.6)) z) (sut/slerp a b t) 1.0e-14)))))
    (t/testing "the shortest arc"
      (let [a sut/ONE
            b (sut/neg (sut/rotation-quaternion 0.2 z))]
        (t/is (sut/delta-eq (sut/rotation-quaternion 0.1 z) (sut/slerp a b 0.5) 1.0e-15) "was a rotation by 2pi - 0.1")
        (t/is (sut/delta-eq (sut/rotation-quaternion 0.05 z) (sut/slerp a b 0.25) 1.0e-15))))
    (t/testing "the result is a unit quaternion"
      (doseq [t [0.01 0.3 0.77 0.999]]
        (t/is (m/delta-eq 1.0 (sut/norm (sut/slerp (sut/rotation-quaternion 0.3 z) (sut/rotation-quaternion 2.4 (v/vec3 1.0 1.0 0.0)) t)) 1.0e-14))))
    (t/testing "almost parallel quaternions"
      (let [a (sut/rotation-quaternion 0.3 z)
            b (sut/rotation-quaternion (+ 0.3 1.0e-6) z)]
        (t/is (sut/delta-eq (sut/rotation-quaternion (+ 0.3 5.0e-7) z) (sut/slerp a b 0.5) 1.0e-12))
        (t/is (sut/delta-eq a (sut/slerp a a 0.4) 1.0e-15) "identical: the same quaternion")))
    (t/testing "zero quaternions and numbers"
      (t/is (= sut/ZERO (sut/slerp sut/ZERO sut/ZERO 0.5)))
      (t/is (sut/delta-eq (sut/rotation-quaternion m/HALF_PI (v/vec3 1.0 0.0 0.0)) (sut/slerp sut/ONE sut/I 0.5) 1.0e-15) "i is the half turn around x: the middle is pi/2 around x")
      (t/is (sut/delta-eq (sut/quaternion (m/cos 0.5) (m/sin 0.5) 0.0 0.0) (sut/slerp 1.0 (c/complex 0.0 1.0) (/ 1.0 m/PI)) 1.0e-12) "complex arguments")))
  (t/testing "nlerp"
    (let [a (sut/rotation-quaternion 0.3 (v/vec3 0.0 0.0 1.0))
          b (sut/neg (sut/rotation-quaternion 0.5 (v/vec3 0.0 0.0 1.0)))]
      (t/is (sut/delta-eq (sut/rotation-quaternion 0.4 (v/vec3 0.0 0.0 1.0)) (sut/nlerp a b 0.5) 1.0e-15) "the shortest arc")
      (t/is (m/delta-eq 1.0 (sut/norm (sut/nlerp a b 0.3)) 1.0e-15))
      (t/is (sut/delta-eq a (sut/nlerp a b 0.0) 1.0e-15))
      (t/is (sut/delta-eq a (sut/nlerp a b -3.0) 1.0e-15) "t is constrained")
      (t/is (sut/delta-eq (sut/neg b) (sut/nlerp a b 5.0) 1.0e-15) "the end of the shortest arc is -b"))))

(t/deftest gimbal-lock
  (let [vectors [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.3 -1.2 2.5]]
        same-rotation? (fn [q1 q2 tol] (every? #(v/delta-eq (sut/rotate % q1) (sut/rotate % q2) tol) vectors))]
    ;; measured rotation errors: 2.4e-9 close to the lock (the dropped roll is of the order of cos(pitch)), 1.2e-9 at 1e-6 from it
    (t/testing "to-euler at pitch +-pi/2 and close to it: the rotation is preserved"
      (doseq [pitch [m/HALF_PI (- m/HALF_PI) (- m/HALF_PI 1.0e-12) (- m/HALF_PI 1.0e-9) (- m/HALF_PI 1.0e-6)
                     (+ (- m/HALF_PI) 1.0e-12) (+ (- m/HALF_PI) 1.0e-9) (+ (- m/HALF_PI) 1.0e-6)]
              roll [-2.0 0.3 1.0]
              yaw [-1.0 0.4 2.5]
              :let [q (sut/from-euler roll pitch yaw)
                    e (sut/to-euler q)
                    tol (if (< (m/abs (- (m/abs pitch) m/HALF_PI)) 1.0e-7) 1.0e-7 1.0e-8)]]
        (t/is (same-rotation? q (sut/from-euler e) tol) (str "euler " roll " " pitch " " yaw " -> " e))
        (t/is (<= (m/abs (.y e)) (+ m/HALF_PI 1.0e-15)))))
    (t/testing "at the lock the roll is zero and the yaw takes the whole rotation"
      (let [e (sut/to-euler (sut/from-euler 0.7 m/HALF_PI 0.2))]
        (t/is (= 0.0 (.x e)))
        (t/is (m/delta-eq m/HALF_PI (.y e) 1.0e-15))
        (t/is (m/delta-eq -0.5 (.z e) 1.0e-12) "yaw - roll"))
      (let [e (sut/to-euler (sut/from-euler 0.7 (- m/HALF_PI) 0.2))]
        (t/is (= 0.0 (.x e)))
        (t/is (m/delta-eq 0.9 (.z e) 1.0e-12) "yaw + roll")))
    (t/testing "to-angles at y = +-pi/2 and close to it"
      (doseq [y [m/HALF_PI (- m/HALF_PI) (- m/HALF_PI 1.0e-12) (- m/HALF_PI 1.0e-9) (- m/HALF_PI 1.0e-6)
                 (+ (- m/HALF_PI) 1.0e-12) (+ (- m/HALF_PI) 1.0e-9) (+ (- m/HALF_PI) 1.0e-6)]
              x [-2.0 0.3 1.0]
              z [-1.0 0.4 2.5]
              :let [q (sut/from-angles x y z)
                    e (sut/to-angles q)
                    tol (if (< (m/abs (- (m/abs y) m/HALF_PI)) 1.0e-7) 1.0e-7 1.0e-8)]]
        (t/is (same-rotation? q (sut/from-angles e) tol) (str "angles " x " " y " " z " -> " e))
        (t/is (<= (m/abs (.y e)) (+ m/HALF_PI 1.0e-15))))
      (let [e (sut/to-angles (sut/from-angles 0.7 m/HALF_PI 0.2))]
        (t/is (= 0.0 (.x e)))
        (t/is (m/delta-eq 0.9 (.z e) 1.0e-12) "x + z"))
      (let [e (sut/to-angles (sut/from-angles 0.7 (- m/HALF_PI) 0.2))]
        (t/is (= 0.0 (.x e)))
        (t/is (m/delta-eq -0.5 (.z e) 1.0e-12) "z - x")))
    (t/testing "far from the lock the round trip is exact"
      (doseq [[r p y] [[0.3 0.5 0.4] [-2.0 -1.0 3.0] [1.0 1.4 -1.0]]]
        (t/is (v/delta-eq (v/vec3 r p y) (sut/to-euler (sut/from-euler r p y)) 1.0e-14))
        (t/is (v/delta-eq (v/vec3 r p y) (sut/to-angles (sut/from-angles r p y)) 1.0e-14))))))

(t/deftest negative-units
  (t/is (= (v/vec4 0.0 -1.0 0.0 0.0) sut/-I))
  (t/is (= (v/vec4 0.0 0.0 -1.0 0.0) sut/-J))
  (t/is (= (v/vec4 0.0 0.0 0.0 -1.0) sut/-K))
  (t/is (= sut/ONE (sut/mult sut/I sut/-I)))
  (t/is (= sut/ONE (sut/mult sut/-J sut/J)))
  (t/is (= sut/-K (sut/mult sut/-I sut/J)))
  (t/is (= sut/-K (sut/mult sut/-J sut/-I)) "(-j)(-i) = ji = -k"))


(t/deftest test-phase-edge-cases
  ;; found in the test phase 2026-10-27 (BR-1, BR-9, BR-10, BR-11)
  (t/testing "subnormal vectors (BR-1)"
    (t/is (sut/delta-eq (v/vec4 0.7071067811865475 0.0 0.0 0.7071067811865475)
                        (sut/from-two-vectors (v/vec3 1.0e-320 0.0 0.0) (v/vec3 0.0 1.0e-320 0.0)) 1.0e-15))
    (t/is (sut/delta-eq (sut/rotation-quaternion 1.0 (v/vec3 0.0 0.0 1.0))
                        (sut/rotation-quaternion 1.0 (v/vec3 0.0 0.0 1.0e-320)) 1.0e-15))
    (t/is (m/delta-eq 1.0 (sut/norm (sut/normalize (sut/quaternion 1.0e-320 1.0e-320 0.0 0.0))) 1.0e-3) "few significant bits of a subnormal number"))
  (t/testing "reciprocal of a subnormal quaternion (BR-9)"
    (let [r (sut/reciprocal (sut/quaternion 1.0e-320 0.0 0.0 0.0))]
      (t/is (= ##Inf (.x r)) "1e320 is not representable")
      (t/is (not (sut/nan? r)) "the vector part is a signed zero, not NaN")
      (t/is (sut/zero? (sut/quaternion 0.0 (.y r) (.z r) (.w r)))))
    (let [q (sut/quaternion 1.0e-308 2.0e-308 -3.0e-308 4.0e-308)
          p (sut/mult q (sut/reciprocal q))]
      (t/is (sut/delta-eq sut/ONE p 1.0e-8) "subnormal inputs with a representable reciprocal, a product close to one"))
    (t/is (sut/nan? (sut/reciprocal sut/ZERO)) "the zero quaternion stays NaN"))
  (t/testing "NaN as the interpolation parameter (BR-10)"
    (t/is (sut/nan? (sut/slerp sut/ONE sut/I ##NaN)))
    (t/is (sut/nan? (sut/nlerp sut/ONE sut/I ##NaN)))
    (t/is (= sut/I (sut/slerp sut/ONE sut/I 1.0)))
    (t/is (= sut/I (sut/slerp sut/ONE sut/I ##Inf)) "infinity is constrained to 1"))
  (t/testing "axes and vectors which are not finite (BR-11)"
    (doseq [bad [(v/vec3 ##Inf 0.0 0.0) (v/vec3 0.0 ##-Inf 1.0) (v/vec3 ##NaN 0.0 1.0)]]
      (t/is (= :u (:option (ex-data (try (sut/rotation-quaternion 1.0 bad) (catch clojure.lang.ExceptionInfo e e))))))
      (t/is (= :u (:option (ex-data (try (sut/from-two-vectors bad (v/vec3 1.0 0.0 0.0)) (catch clojure.lang.ExceptionInfo e e))))))
      (t/is (= :v (:option (ex-data (try (sut/from-two-vectors (v/vec3 1.0 0.0 0.0) bad) (catch clojure.lang.ExceptionInfo e e)))))))))

(t/deftest normalize-longer-than-max
  ;; RV-17: the norm overflows, the versor does not
  (t/is (sut/delta-eq (v/vec4 0.5 0.5 -0.5 0.5) (sut/normalize (sut/quaternion 1.7e308 1.7e308 -1.7e308 1.7e308)) 1.0e-15))
  (t/is (= ##Inf (sut/norm (sut/quaternion 1.7e308 1.7e308 -1.7e308 1.7e308))) "the length itself is not representable")
  (t/is (m/delta-eq 1.0 (sut/norm (sut/normalize (sut/quaternion 1.7e308 1.7e308 1.7e308 1.7e308))) 1.0e-15)))

(t/deftest directions-longer-than-max
  ;; BR-12: a finite vector is a direction also when its length overflows
  (let [long-axis (v/vec3 1.7e308 1.7e308 1.7e308)
        unit-axis (v/vec3 1.0 1.0 1.0)]
    (t/is (sut/delta-eq (sut/rotation-quaternion 1.0 unit-axis) (sut/rotation-quaternion 1.0 long-axis) 1.0e-15))
    (t/is (sut/delta-eq (sut/from-axis-angle 1.0 unit-axis) (sut/from-axis-angle 1.0 long-axis) 1.0e-15))
    (let [q (sut/from-two-vectors long-axis (v/vec3 -1.7e308 0.0 0.0))]
      (t/is (m/delta-eq 1.0 (sut/norm q) 1.0e-15))
      (t/is (v/delta-eq (v/vec3 -1.0 0.0 0.0) (sut/rotate (v/normalize long-axis) q) 1.0e-14)
            "rotates the direction of the first vector to the direction of the second"))
    (t/is (sut/delta-eq (sut/from-two-vectors (v/vec3 1.0 1.0 1.0) (v/vec3 0.0 1.0 0.0))
                        (sut/from-two-vectors long-axis (v/vec3 0.0 1.7e308 0.0)) 1.0e-15)))
  (t/testing "still rejected: zero and non finite"
    (doseq [bad [(v/vec3 0.0 0.0 0.0) (v/vec3 ##Inf 1.0 0.0) (v/vec3 1.7e308 ##-Inf 0.0) (v/vec3 ##NaN 1.0 1.0)]]
      (t/is (= :u (:option (ex-data (try (sut/rotation-quaternion 1.0 bad) (catch clojure.lang.ExceptionInfo e e)))))))))

