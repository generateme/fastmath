(ns fastmath.solver.closed-form-test
  (:require [fastmath.core :as m]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.solver.closed-form :as cf])
  (:import [fastmath.vector Vec2 Vec3]))

(defn- close?
  "Relative (or absolute near zero) closeness of two vectors, NaN equal to NaN."
  [expected actual tol]
  (and (= (count expected) (count actual))
       (every? true? (map (fn [e a]
                            (let [e (double e) a (double a)]
                              (cond (and (m/nan? e) (m/nan? a)) true
                                    (or (m/nan? e) (m/nan? a)) false
                                    (== e a) true
                                    :else (<= (m/abs (- e a)) (* tol (max 1.0 (m/abs e)))))))
                          expected actual))))

(defn- rel-close? [e a tol] (or (== e a) (<= (m/abs (- e a)) (* tol (m/abs e)))))

(t/deftest quadratic-basic
  (t/is (instance? Vec2 (cf/quadratic 1 0 -4)))
  (t/is (= [-2.0 2.0] (vec (cf/quadratic 1 0 -4))))
  (t/is (= [-1.0 -1.0] (vec (cf/quadratic 1 2 1))) "double root twice")
  (t/is (= [1.0 2.0] (vec (cf/quadratic 1 -3 2))))
  (t/is (= [1.0 2.0] (vec (cf/quadratic -1 3 -2))) "negative leading coefficient")
  (t/is (= [1.0 2.0] (vec (cf/quadratic 0.5 -1.5 1.0))))
  (t/is (close? [1.0e-8 1.0e8] (vec (cf/quadratic 1 -1e8 1)) 1e-15) "no cancellation in the small root")
  (t/is (close? [-1.618033988749895 0.6180339887498949] (vec (cf/quadratic 1 1 -1)) 1e-15))
  (t/is (= [-0.5 3.0] (vec (cf/quadratic 2 -5 -3))) "rational roots")
  (t/is (= [1.0 2.0] (vec (cf/quadratic 1.0 -3.0 2.0))) "doubles")
  (t/is (= [1.0 2.0] (vec (cf/quadratic 1N -3 2))) "other number types"))

(t/deftest quadratic-boundaries
  ;; zeros
  (t/is (= [0.0 0.0] (vec (cf/quadratic 1 0 0))) "x^2 = 0")
  (t/is (= [0 0] (map #(m/double-bits %) (cf/quadratic 1 0 0))) "no negative zero")
  (t/is (= [-1.0 0.0] (vec (cf/quadratic 1 1 0))))
  (t/is (= [0.0 1.0] (vec (cf/quadratic -1 1 0))) "sorted")
  (t/is (= [0.0 0.0] (vec (cf/quadratic 5 -0.0 0.0))) "negative zero coefficients")
  (t/is (= [-2.0 2.0] (vec (cf/quadratic 1 0.0 -4))))
  ;; no real roots
  (t/is (every? m/nan? (cf/quadratic 1 0 1)) "x^2 + 1")
  (t/is (every? m/nan? (cf/quadratic 1 1 1)))
  (t/is (every? m/nan? (cf/quadratic 1e-300 1e-300 1e-300)) "no real roots at a tiny scale")
  (t/is (every? m/nan? (cf/quadratic 1e300 1e300 1e300)))
  ;; degree lower than two
  (t/is (close? [-2.0 ##NaN] (vec (cf/quadratic 0 2 4)) 0.0) "linear")
  (t/is (close? [0.0 ##NaN] (vec (cf/quadratic 0 2 0)) 0.0) "linear, root zero")
  (t/is (every? m/nan? (cf/quadratic 0 0 0)) "0 = 0")
  (t/is (every? m/nan? (cf/quadratic 0 0 4)) "4 = 0")
  (t/is (close? [-1e-300 ##NaN] (vec (cf/quadratic 0 1e300 1.0)) 1e-12) "linear: no overflow in the ratio")
  ;; non finite
  (doseq [bad [##NaN ##Inf ##-Inf]]
    (t/is (every? m/nan? (cf/quadratic bad 1 1)) (str "a = " bad))
    (t/is (every? m/nan? (cf/quadratic 1 bad 1)) (str "b = " bad))
    (t/is (every? m/nan? (cf/quadratic 1 1 bad)) (str "c = " bad)))
  ;; the discriminant is exactly zero / just below
  (t/is (= [0.5 0.5] (vec (cf/quadratic 4 -4 1))))
  (t/is (every? m/nan? (cf/quadratic 1 (m/prev-double 2.0) 1)) "discriminant just below zero is not real")
  (t/is (every? #(not (m/nan? %)) (cf/quadratic 1 (m/next-double 2.0) 1)) "just above zero is real"))

(t/deftest quadratic-scale
  (doseq [k [-1000 -500 -100 -1 1 100 500 1000 1020]
          :let [s (m/scalb 1.0 k)]]
    (t/is (= [1.0 2.0] (vec (cf/quadratic s (* -3.0 s) (* 2.0 s)))) (str "scaled by 2^" k)))
  ;; the coefficients of a very different size
  (t/is (close? [-1.0 1.0] (vec (cf/quadratic 1e200 0 -1e200)) 1e-15))
  (t/is (close? [-1e200 -1e-200] (vec (cf/quadratic 1 1e200 1)) 1e-15))
  (t/is (close? [-1.618033988749895 0.6180339887498949] (vec (cf/quadratic 1e-200 1e-200 -1e-200)) 1e-15))
  (t/is (close? [-1.618033988749895 0.6180339887498949] (vec (cf/quadratic 1e200 1e200 -1e200)) 1e-15))
  (t/is (close? [-1.618033988749895 0.6180339887498949] (vec (cf/quadratic Double/MIN_VALUE Double/MIN_VALUE (- Double/MIN_VALUE))) 1e-15) "subnormal coefficients")
  (t/is (close? [-1e100 1e100] (vec (cf/quadratic 1e-100 0 -1e100)) 1e-12) "roots of 1e100")
  (t/is (close? [-1e-100 1e-100] (vec (cf/quadratic 1e100 0 -1e-100)) 1e-12) "roots of 1e-100")
  ;; roots which are not representable
  (let [[r0 r1] (vec (cf/quadratic 1e-300 1e300 1e-300))]
    (t/is (== ##-Inf r0) "a root beyond the range")
    (t/is (zero? r1) "a root below the range"))
  (let [[r0 r1] (vec (cf/quadratic 1 1e300 1e-300))]
    (t/is (== -1e300 r0))
    (t/is (zero? r1) "the root -1e-600 underflows"))
  ;; invariance: the same polynomial multiplied by a power of two gives the same roots bit by bit
  (doseq [[a b c] [[1.0 -3.0 2.0] [0.1 0.3 -0.7] [3e5 -1e-3 -1e-9]]
          k [-300 -50 50 300]
          :let [s (m/scalb 1.0 k)]]
    (t/is (= (vec (cf/quadratic a b c)) (vec (cf/quadratic (* a s) (* b s) (* c s)))) (str [a b c] " 2^" k))))

(t/deftest quadratic-accuracy-random
  ;; roots r1 r2 of very different sizes: the backward stable formula keeps the relative accuracy of both
  (let [rnd (java.util.Random. 20261026)]
    (dotimes [_ 2000]
      (let [e1 (- (* 20.0 (.nextDouble rnd)) 10.0)
            e2 (- (* 20.0 (.nextDouble rnd)) 10.0)
            r1 (* (if (.nextBoolean rnd) 1 -1) (m/pow 10.0 e1))
            r2 (* (if (.nextBoolean rnd) 1 -1) (m/pow 10.0 e2))
            ;; the exact coefficients (BigDecimal) rounded to doubles: the roots move by a relative 1e-16 at most for the product and the sum
            a 1.0
            b (- (+ r1 r2))
            c (* r1 r2)
            [x0 x1] (vec (cf/quadratic a b c))
            [e0 e1'] (sort [r1 r2])]
        ;; b and c carry rounding errors of a relative 1e-16, which is a relative error of at most 1e-16 * cond of the roots
        (t/is (and (rel-close? e0 x0 1e-9) (rel-close? e1' x1 1e-9)) (str [r1 r2] " -> " [x0 x1]))))))

(t/deftest cubic-basic
  (t/is (instance? Vec3 (cf/cubic 1 -6 11 -6)))
  (t/is (close? [1.0 2.0 3.0] (vec (cf/cubic 1 -6 11 -6)) 1e-14))
  (t/is (close? [1.0 2.0 3.0] (vec (cf/cubic -2 12 -22 12)) 1e-14) "negative leading coefficient")
  (t/is (close? [2.0 ##NaN ##NaN] (vec (cf/cubic 1 0 0 -8)) 1e-14) "one real root")
  (t/is (close? [-1.0 ##NaN ##NaN] (vec (cf/cubic 1 0 0 1)) 1e-14))
  (t/is (close? [-2.0 1.0 1.0] (vec (cf/cubic 1 0 -3 2)) 1e-7) "double root")
  (t/is (close? [-1.0 -1.0 2.0] (vec (cf/cubic 1 0 -3 -2)) 1e-7) "double root, the other side")
  (t/is (close? [1.0 1.0 1.0] (vec (cf/cubic 1 -3 3 -1)) 1e-4) "a triple root is known to the cube root of the rounding error")
  (t/is (close? [0.0 0.0 0.0] (vec (cf/cubic 1 0 0 0)) 0.0) "x^3")
  (t/is (close? [-0.5 0.25 2.0] (vec (cf/cubic 8 -14 -5 2)) 1e-14) "rational roots")
  (t/is (close? [1.0 2.0 3.0] (vec (cf/cubic 1N -6N 11N -6N)) 1e-14) "other number types")
  (t/is (apply <= (cf/cubic 1 -6 11 -6)) "ascending"))

(t/deftest cubic-boundaries
  ;; a zero constant term: an exact zero root
  (t/is (= [-1.0 0.0 1.0] (vec (cf/cubic 1 0 -1 0))) "x^3 - x")
  (t/is (= [-1.0 0.0 0.0] (vec (cf/cubic 1 1 0 0))) "x^3 + x^2")
  (t/is (close? [0.0 ##NaN ##NaN] (vec (cf/cubic 1 0 1 0)) 0.0) "x^3 + x: one real root 0")
  (t/is (close? [0.0 2.0 3.0] (vec (cf/cubic 1 -5 6 0)) 1e-14))
  (t/is (close? [0.0 (/ 2.0 3.0) 3e20] (vec (cf/cubic 1e-20 -3 2 0)) 1e-14) "a tiny leading coefficient and d = 0")
  ;; degree lower than three
  (t/is (close? [-2.0 2.0 ##NaN] (vec (cf/cubic 0 1 0 -4)) 1e-15) "quadratic")
  (t/is (every? m/nan? (cf/cubic 0 1 0 1)) "quadratic without real roots: no exception")
  (t/is (every? m/nan? (cf/cubic 0 1 2 3)))
  (t/is (close? [2.0 ##NaN ##NaN] (vec (cf/cubic 0 0 2 -4)) 1e-15) "linear")
  (t/is (every? m/nan? (cf/cubic 0 0 0 1)))
  (t/is (every? m/nan? (cf/cubic 0 0 0 0)))
  ;; negative zero
  (t/is (close? [1.0 2.0 3.0] (vec (cf/cubic 1 -6 11 -6.0)) 1e-14))
  (t/is (close? [-1.0 1.0 ##NaN] (vec (cf/cubic -0.0 1 0 -1)) 1e-15) "a = -0.0")
  ;; non finite
  (doseq [bad [##NaN ##Inf ##-Inf]]
    (t/is (every? m/nan? (cf/cubic bad 1 1 1)) (str "a = " bad))
    (t/is (every? m/nan? (cf/cubic 1 bad 1 1)))
    (t/is (every? m/nan? (cf/cubic 1 1 bad 1)))
    (t/is (every? m/nan? (cf/cubic 1 1 1 bad))))
  ;; the arc cosine of the trigonometric formula exceeded 1 by rounding: the real root must not be lost
  (let [coefs [1.0 -0.15435828008735172 -0.02225630683100868 -6.02300324876155E-4]
        r (vec (cf/cubic 1.0 -0.15435828008735172 -0.02225630683100868 -6.02300324876155E-4))
        x (first (remove m/nan? r))]
    (t/is (some? x) "a real root is returned")
    (t/is (some #(rel-close? 0.25211322 % 1e-7) (remove m/nan? r)) "the root 0.25211322 is among the roots")
    (doseq [x (remove m/nan? r)
            :let [p (reduce (fn [acc k] (+ (* acc x) k)) 0.0 coefs)]]
      (t/is (< (m/abs p) 1e-15) "a small residual")))
  ;; a root which is real to the rounding error
  (let [r (vec (cf/cubic 1 (- -5.0 1e-10) (+ 8.0 2e-10) -4.0))]
    (t/is (= 3 (count (remove m/nan? r))) "three real roots 1, 2, 2 + 1e-10 (double within the rounding)")
    (t/is (close? [1.0 2.0 2.0] r 1e-6))))

(t/deftest cubic-scale
  (doseq [k [-1000 -300 -50 50 300 1000]
          :let [s (m/scalb 1.0 k)]]
    (t/is (= (vec (cf/cubic 1 -6 11 -6)) (vec (cf/cubic s (* -6.0 s) (* 11.0 s) (* -6.0 s)))) (str "coefficients scaled by 2^" k))
    (t/is (close? (vec (cf/cubic 2 -3 -1 5)) (vec (cf/cubic (* 2 s) (* -3 s) (* -1 s) (* 5 s))) 0.0) (str "one real root, scaled by 2^" k)))
  ;; the roots s, 2s, 3s
  (doseq [s [1e-100 1e-60 1e-10 1e10 1e60 1e100]]
    (t/is (close? [s (* 2 s) (* 3 s)] (vec (cf/cubic 1.0 (* -6.0 s) (* 11.0 s s) (* -6.0 s s s))) 1e-13) (str "roots of the size " s)))
  ;; roots of very different sizes
  (doseq [[r0 r1 r2] [[1e-6 1.0 1e6] [1e-8 1.0 1e8] [-1e-6 1.0 1e6] [-1e5 -1e-3 1e2]]
          :let [rs (vec (cf/cubic 1.0 (- (+ r0 r1 r2)) (+ (* r0 r1) (* r0 r2) (* r1 r2)) (- (* r0 r1 r2))))]]
    (t/is (every? true? (map #(rel-close? %1 %2 1e-9) (sort [r0 r1 r2]) rs)) (str [r0 r1 r2] " -> " rs)))
  ;; roots which cannot be represented
  (t/is (every? #(not (m/nan? %)) (remove m/nan? (cf/cubic 1e-320 1 -3 2))) "a tiny leading coefficient: no NaN in the roots found")
  (let [rs (vec (cf/cubic 1e-300 -3 2 1e-300))]
    (t/is (rel-close? -5.0e-301 (rs 0) 1e-12) "a tiny root of a tiny constant term")
    (t/is (rel-close? (/ 2.0 3.0) (rs 1) 1e-12))
    (t/is (rel-close? 3e300 (rs 2) 1e-12) "a huge root of a tiny leading coefficient")))

(t/deftest cubic-random-residuals
  ;; every returned root has a small componentwise backward error
  (let [rnd (java.util.Random. 20261027)
        eps 2.220446049250313E-16]
    (dotimes [_ 3000]
      (let [co (vec (repeatedly 4 #(* (if (.nextBoolean rnd) 1 -1) (m/pow 10.0 (- (* 8.0 (.nextDouble rnd)) 4.0)))))
            r (vec (cf/cubic (co 0) (co 1) (co 2) (co 3)))]
        (doseq [x r :when (not (m/nan? x))]
          (let [p (reduce (fn [acc c] (+ (* acc x) c)) 0.0 co)
                s (reduce (fn [acc c] (+ (* acc (m/abs x)) (m/abs c))) 0.0 co)]
            (t/is (<= (/ (m/abs p) s) (* 200 eps)) (str co " x=" x))))
        (t/is (apply <= (remove m/nan? r)) "ascending")
        (t/is (seq (remove m/nan? r)) (str "a cubic has a real root: " co))))))

;; roots of random polynomials computed with mpmath from the exact double coefficients (120 digits),
;; generated by utils/fastmath/dev/generate_solver_closed_form_reference.py
(def reference (delay (edn/read-string (slurp (io/resource "solver/closed_form_reference.edn")))))

(t/deftest reference-roots
  (doseq [[kind solve tol] [[:quadratic cf/quadratic 1e-15] [:cubic cf/cubic 1e-12]]
          {:keys [coefs roots]} (get @reference kind)
          :let [got (vec (remove m/nan? (apply solve coefs)))]]
    (t/is (= (count roots) (count got)) (str kind " " coefs " roots " roots " got " got))
    (when (= (count roots) (count got))
      ;; relative error per root; observed maximum 3.1e-16 (quadratic) and 5.2e-14 (cubic, the conditioning of the polynomial), tolerances 1e-15 and 1e-12
      (t/is (every? true? (map #(rel-close? %1 %2 tol) roots got)) (str kind " " coefs " roots " roots " got " got)))))
