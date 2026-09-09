(ns fastmath.special-test
  (:require [fastmath.special :as sut]
            [clojure.test :as t]
            [clojure.java.io :as io]
            [clojure.edn :as edn]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]))

;; Reference values for the Bessel/Airy comparisons below were computed once
;; from R (`base` and `Bessel` packages) and are stored in
;; `test/resources/special/reference.edn`. See
;; `utils/fastmath/dev/special_ref_gen.clj` for the generator (requires the
;; `:dev` profile and a working R + `Bessel` package installation) if the
;; reference data ever needs to be regenerated or extended.
;;
;; Reference source per input: R `base::bessel*` (vectorised, x >= 0) or
;; `Bessel::Bessel*` (per element, handles negative arguments/orders).

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "special/reference.edn")))))

(defn- blk
  "Fetch a precomputed comparison block `{:order? :arg :ref}` from the reference data."
  [test-key block-key]
  (get-in @reference [test-key block-key]))

(def ^:private ABS 1.0e-9)

(defn- check1
  "Compare `(f arg)` against a single-argument reference block."
  ([f b] (check1 f b 1.0e-11))
  ([f b rel] (let [{:keys [arg ref]} b] (v/edelta-eq (mapv f arg) ref ABS rel))))

(defn- check2
  "Compare `(f order arg)` against a two-argument reference block."
  ([f b] (check2 f b 1.0e-10))
  ([f b rel] (let [{:keys [order arg ref]} b] (v/edelta-eq (mapv f order arg) ref ABS rel))))

(def xl [1e12, 5e12, 1e13, 5e13, 1e14, 5e14, 1e15, 5e15, 1e16, 5e16, 1e17, 5e17, 1e18, 5e18, 1e19, 5e19, 1e20, 1e22, 1e25, 1e30, 1e40])

(t/deftest bessel-J0
  (t/is (m/nan? (sut/bessel-J0 ##NaN)))
  (t/is (m/one? (sut/bessel-J0 0.0)))
  (t/is (m/zero? (sut/bessel-J0 ##Inf)))
  (t/is (m/zero? (sut/bessel-J0 ##-Inf)))
  ;; positive & negative
  (t/is (check1 sut/bessel-J0 (blk :bessel-J0 :pos)))
  (t/is (check1 sut/bessel-J0 (blk :bessel-J0 :neg)))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-J0 xl) [1.016712505004068e-7
                                             -2.1276975389854557e-7
                                             1.192648473966565e-7
                                             -7.094408384899425e-8
                                             -6.698265203680451e-8
                                             -3.4394958970536735e-8
                                             6.156638646885022e-9
                                             -1.0644174809027939e-8
                                             8.661427680921681e-10
                                             -1.8463226228625315e-9
                                             -2.4087235483673835e-9
                                             6.690193890936787e-11
                                             -4.934387036790141e-10
                                             -6.814178646566312e-11
                                             -2.3228731060101177e-10
                                             -4.7522631032150424e-11
                                             6.6980090407034224e-12
                                             -1.856105106510822e-12
                                             1.1543496219672643e-13
                                             -5.589003016686146e-16
                                             -6.538288347442135e-22] 1.0e-8 1.0e-8)))

(t/deftest bessel-J1
  (t/is (m/nan? (sut/bessel-J1 ##NaN)))
  (t/is (m/zero? (sut/bessel-J1 0.0)))
  (t/is (m/zero? (sut/bessel-J0 ##Inf)))
  (t/is (m/zero? (sut/bessel-J0 ##-Inf)))
  ;; positive & negative
  (t/is (check1 sut/bessel-J1 (blk :bessel-J1 :pos)))
  (t/is (check1 sut/bessel-J1 (blk :bessel-J1 :neg)))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-J1 xl) [-7.913802683850442e-7
                                             2.8644892441665137e-7
                                             -2.2234629165382475e-7
                                             -8.774583986821614e-8
                                             4.3353454877231535e-8
                                             -9.498754768402967e-9
                                             2.4468665123771328e-8
                                             -3.745063031294846e-9
                                             7.931694266803266e-9
                                             -3.0534387532186807e-9
                                             7.511648229358563e-10
                                             -1.1263941030142736e-9
                                             -6.270071914094412e-10
                                             -3.5025797836849376e-10
                                             -9.85118397478559e-11
                                             1.0234253752537956e-10
                                             -7.95068198242545e-11
                                             -7.759951744073064e-12
                                             -2.2435852276969216e-13
                                             5.694297368089524e-16
                                             7.952011386537066e-21] 1.0e-8 1.0e-8)))

(t/deftest bessel-J
  (t/is (m/nan? (sut/bessel-J 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-J -0.5 ##NaN)))
  (t/is (every? m/zero? (map #(sut/bessel-J % 0.0) (range 0.2 100.1 0.1))))
  (t/is (every? m/zero? (map #(sut/bessel-J % 0.0) (range -100 0))))
  ;; integer orders
  (t/is (check2 sut/bessel-J (blk :bessel-J :int-xx)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :int-x)))
  ;; fractional orders
  (t/is (check2 sut/bessel-J (blk :bessel-J :frac-xx)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :frac-x)))
  ;; large orders
  (t/is (check2 sut/bessel-J (blk :bessel-J :large-xx) 1.0e-9))
  (t/is (check2 sut/bessel-J (blk :bessel-J :large-x) 1.0e-9))
  ;; order sweeps, fixed arg
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-pos-015)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-pos-21)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-pos-421)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-pos-1421)))
  ;; negative order sweeps, fixed arg
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-neg-015)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-neg-21)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-neg-421)))
  (t/is (check2 sut/bessel-J (blk :bessel-J :vs-neg-1421)))
  ;; negative order, negative arg
  (t/is (check2 sut/bessel-J (blk :bessel-J :neg-ord)))
  ;; explicit (order, arg) pairs
  (t/is (check2 sut/bessel-J (blk :bessel-J :are))))

;;;;

(t/deftest bessel-Y0
  (t/is (m/nan? (sut/bessel-Y0 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-Y0 0.0)))
  (t/is (m/zero? (sut/bessel-Y0 ##Inf)))
  (t/is (m/nan? (sut/bessel-Y0 ##-Inf)))
  ;; positive
  (t/is (check1 sut/bessel-Y0 (blk :bessel-Y0 :pos)))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-Y0 xl) [ -7.91380268385095e-7
                                             2.8644892441667265e-7
                                             -2.223462916538307e-7
                                             -8.774583986821542e-8
                                             4.335345487723187e-8
                                             -9.498754768402932e-9
                                             2.446866512377132e-8
                                             -3.745063031294845e-9
                                             7.931694266803265e-9
                                             -3.0534387532186807e-9
                                             7.511648229358563e-10
                                             -1.1263941030142736e-9
                                             -6.270071914094412e-10
                                             -3.5025797836849376e-10
                                             -9.85118397478559e-11
                                             1.0234253752537956e-10
                                             -7.95068198242545e-11
                                             -7.759951744073064e-12
                                             -2.2435852276969216e-13
                                             5.694297368089524e-16
                                             7.952011386537066e-21] 1.0e-8 1.0e-8)))

(t/deftest bessel-Y1
  (t/is (m/nan? (sut/bessel-Y1 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-Y1 0.0)))
  (t/is (m/zero? (sut/bessel-Y1 ##Inf)))
  (t/is (m/nan? (sut/bessel-Y1 ##-Inf)))
  ;; positive
  (t/is (check1 sut/bessel-Y1 (blk :bessel-Y1 :pos)))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-Y1 xl) [-1.0167125050080249e-7
                                             2.127697538985742e-7
                                             -1.1926484739666762e-7
                                             7.094408384899338e-8
                                             6.698265203680473e-8
                                             3.439495897053673e-8
                                             -6.15663864688501e-9
                                             1.0644174809027939e-8
                                             -8.661427680921677e-10
                                             1.8463226228625315e-9
                                             2.4087235483673835e-9
                                             -6.690193890936787e-11
                                             4.934387036790141e-10
                                             6.814178646566312e-11
                                             2.3228731060101177e-10
                                             4.7522631032150424e-11
                                             -6.6980090407034224e-12
                                             1.856105106510822e-12
                                             -1.1543496219672643e-13
                                             5.589003016686146e-16
                                             6.538288347442135e-22] 1.0e-8 1.0e-8)))

(t/deftest bessel-Y
  (t/is (m/nan? (sut/bessel-Y 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-Y -0.5 ##NaN)))
  ;; integer orders
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :int-xx)))
  ;; fractional orders
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :frac-xx)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :frac-x)))
  ;; large orders
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :large-xx) 1.0e-9))
  ;; order sweeps, fixed arg
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-pos-015)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-pos-21)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-pos-421)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-pos-1421)))
  ;; negative order sweeps, fixed arg
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-neg-015)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-neg-21)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-neg-421)))
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :vs-neg-1421)))
  ;; explicit (order, arg) pairs
  (t/is (check2 sut/bessel-Y (blk :bessel-Y :are))))

;;

(t/deftest bessel-I0
  (t/is (m/nan? (sut/bessel-I0 ##NaN)))
  (t/is (m/one? (sut/bessel-I0 0.0)))
  (t/is (m/nan? (sut/bessel-I0 ##Inf)))
  (t/is (m/nan? (sut/bessel-I0 ##-Inf)))
  ;; positive & negative
  (t/is (check1 sut/bessel-I0 (blk :bessel-I0 :pos)))
  (t/is (check1 sut/bessel-I0 (blk :bessel-I0 :neg))))

(t/deftest bessel-I1
  (t/is (m/nan? (sut/bessel-I1 ##NaN)))
  (t/is (m/zero? (sut/bessel-I1 0.0)))
  (t/is (m/nan? (sut/bessel-I1 ##Inf)))
  (t/is (m/nan? (sut/bessel-I1 ##-Inf)))
  ;; positive & negative
  (t/is (check1 sut/bessel-I1 (blk :bessel-I1 :pos)))
  (t/is (check1 sut/bessel-I1 (blk :bessel-I1 :neg))))

;;

(t/deftest bessel-I
  (t/is (m/nan? (sut/bessel-I 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-I -0.5 ##NaN)))
  ;; general orders
  (t/is (check2 sut/bessel-I (blk :bessel-I :gen-xx)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :gen-x)))
  ;; fractional orders
  (t/is (check2 sut/bessel-I (blk :bessel-I :frac-xx)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :frac-x)))
  ;; order sweeps, fixed arg
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-pos-015)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-pos-21)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-pos-421)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-pos-1421)))
  ;; negative order sweeps, fixed arg
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-neg-015)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-neg-21)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-neg-421)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :vs-neg-1421)))
  ;; negative order, negative & positive arg
  (t/is (check2 sut/bessel-I (blk :bessel-I :negord-xx)))
  (t/is (check2 sut/bessel-I (blk :bessel-I :negord-axx)))
  ;; explicit (order, arg) pairs
  (t/is (check2 sut/bessel-I (blk :bessel-I :are))))

;;

(t/deftest bessel-K0
  (t/is (m/nan? (sut/bessel-K0 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-K0 0.0)))
  (t/is (m/zero? (sut/bessel-K0 ##Inf)))
  (t/is (m/nan? (sut/bessel-K0 ##-Inf)))
  ;; positive
  (t/is (check1 sut/bessel-K0 (blk :bessel-K0 :pos))))

(t/deftest bessel-K1
  (t/is (m/nan? (sut/bessel-K1 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-K1 0.0)))
  (t/is (m/zero? (sut/bessel-K1 ##Inf)))
  (t/is (m/nan? (sut/bessel-K1 ##-Inf)))
  ;; positive
  (t/is (check1 sut/bessel-K1 (blk :bessel-K1 :pos))))

(t/deftest bessel-K
  (t/is (m/nan? (sut/bessel-K 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-K -0.5 ##NaN)))
  ;; order/arg sweep
  (t/is (check2 sut/bessel-K (blk :bessel-K :nu-sweep)))
  ;; fractional orders
  (t/is (check2 sut/bessel-K (blk :bessel-K :frac-xx)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :frac-x)))
  ;; order sweeps, fixed arg
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-pos-015)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-pos-21)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-pos-421)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-pos-1421)))
  ;; negative order sweeps, fixed arg
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-neg-015)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-neg-21)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-neg-421)))
  (t/is (check2 sut/bessel-K (blk :bessel-K :vs-neg-1421)))
  ;; negative order, positive arg
  (t/is (check2 sut/bessel-K (blk :bessel-K :negord-axx)))
  ;; explicit (order, arg) pairs
  (t/is (check2 sut/bessel-K (blk :bessel-K :are))))

(t/deftest bessel-K-half
  (doseq [o (range 1 100 2)
          x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          :let [xx (* o x)
                oh (* o 0.5)]]
    (t/is (m/delta-eq (sut/bessel-K oh xx) (sut/bessel-K-half-odd o xx) 1.0e-13 1.0e-13))))

;;

;; `digamma`/`trigamma` reference values below were computed independently
;; with Python `scipy.special.digamma`/`polygamma(1, x)` and are stored in
;; `test/resources/special/digamma_trigamma_reference.edn`. See
;; `utils/fastmath/dev/generate_digamma_trigamma_reference.py` for the
;; generator. In addition to that external check, both functions are also
;; validated against their recurrence and reflection identities over dense
;; grids (input grids are offset by a fractional amount to avoid landing on
;; the non-positive-integer poles).

(def ^:private gamma-reference
  (delay (edn/read-string (slurp (io/resource "special/digamma_trigamma_reference.edn")))))

(def ^:private gamma-recurrence-xs
  (concat (range 0.05 10.0 0.1) (range -9.95 -0.05 0.1) [1.0e3 1.0e6 -1000.37 -1000000.37]))

(def ^:private gamma-reflection-xs
  (concat (range 0.05 5.0 0.1) (range -4.95 -0.05 0.1)))

(t/deftest digamma
  (t/is (m/nan? (sut/digamma ##NaN)))
  (t/is (m/nan? (sut/digamma ##-Inf)))
  (t/is (m/pos-inf? (sut/digamma ##Inf)))
  (t/is (m/neg-inf? (sut/digamma 0.0)))
  (t/is (m/neg-inf? (sut/digamma -0.0)))
  (t/testing "closed-form special values"
    (t/is (m/delta-eq (m/- m/GAMMA) (sut/digamma 1.0) 1.0e-14))
    (t/is (m/delta-eq (m/- 1.0 m/GAMMA) (sut/digamma 2.0) 1.0e-14))
    (t/is (m/delta-eq (m/- (m/- m/GAMMA) (m/* 2.0 m/LN2)) (sut/digamma 0.5) 1.0e-14))
    (t/is (m/delta-eq (m/+ 2.0 (m/- (m/- m/GAMMA) (m/* 2.0 m/LN2))) (sut/digamma 1.5) 1.0e-13)))
  (t/testing "recurrence: digamma(x+1) = digamma(x) + 1/x"
    (doseq [x gamma-recurrence-xs]
      (t/is (m/delta-eq (m/+ (sut/digamma x) (m// x)) (sut/digamma (m/inc x)) 1.0e-7))))
  (t/testing "reflection: digamma(1-x) - digamma(x) = PI*cot(PI*x)"
    (doseq [x gamma-reflection-xs]
      (t/is (m/delta-eq (m/* m/PI (m/cot (m/* m/PI x)))
                        (m/- (sut/digamma (m/- 1.0 x)) (sut/digamma x)) 1.0e-8))))
  (t/testing "vs scipy"
    (t/is (check1 sut/digamma {:arg (:arg @gamma-reference) :ref (:digamma @gamma-reference)} 1.0e-8))))

(t/deftest trigamma
  (t/is (m/nan? (sut/trigamma ##NaN)))
  (t/is (m/nan? (sut/trigamma ##-Inf)))
  (t/is (m/zero? (sut/trigamma ##Inf)))
  (t/is (m/pos-inf? (sut/trigamma 0.0)))
  (t/is (m/pos-inf? (sut/trigamma -0.0)))
  (t/testing "closed-form special values"
    (t/is (m/delta-eq (m// (m/sq m/PI) 6.0) (sut/trigamma 1.0) 1.0e-14))
    (t/is (m/delta-eq (m/- (m// (m/sq m/PI) 6.0) 1.0) (sut/trigamma 2.0) 1.0e-14))
    (t/is (m/delta-eq (m// (m/sq m/PI) 2.0) (sut/trigamma 0.5) 1.0e-13))
    (t/is (m/delta-eq (m/- (m// (m/sq m/PI) 2.0) 4.0) (sut/trigamma 1.5) 1.0e-13)))
  (t/testing "recurrence: trigamma(x+1) = trigamma(x) - 1/x^2"
    (doseq [x gamma-recurrence-xs]
      (t/is (m/delta-eq (m/- (sut/trigamma x) (m// (m/sq x))) (sut/trigamma (m/inc x)) 1.0e-6))))
  (t/testing "reflection: trigamma(x) + trigamma(1-x) = (PI/sin(PI*x))^2"
    (doseq [x gamma-reflection-xs]
      (t/is (m/delta-eq (m/sq (m// m/PI (m/sin (m/* m/PI x))))
                        (m/+ (sut/trigamma x) (sut/trigamma (m/- 1.0 x))) 1.0e-9))))
  (t/testing "vs scipy"
    (t/is (check1 sut/trigamma {:arg (:arg @gamma-reference) :ref (:trigamma @gamma-reference)} 1.0e-7))))

;;

;; `beta`/`regularized-beta`/`incomplete-beta` reference values below were
;; computed independently with Python `mpmath` (`mpmath.beta`,
;; `mpmath.betainc`, which support analytic continuation for negative
;; parameters) and are stored in `test/resources/special/beta_reference.edn`.
;; See `utils/fastmath/dev/generate_beta_reference.py` for the generator.
;;
;; Domain notes for negative arguments (see also `beta`'s docstring):
;;  - `beta`: a genuine (non-removable) pole occurs when `p` or `q` itself is
;;    a non-positive integer. When only `p+q` is a non-positive integer, the
;;    singularity is removable and the correct value is `0.0`.
;;  - `regularized-beta`: a genuine pole occurs whenever `a+b` is a
;;    non-positive integer (even if `a`,`b` individually are not) - unlike
;;    `beta`, this is *not* removable, since it corresponds to dividing a
;;    finite incomplete-beta value by a vanishing `beta(a,b)`.
;;  - `incomplete-beta` (unnormalized) has no such issue: its negative-domain
;;    formula never divides by `beta(a,b)`, so it stays finite even where
;;    `regularized-beta` is singular.
;; Reference grids avoid true poles (`p`, `q`, `p+q` / `a`, `b`, `a+b` being
;; non-positive integers) accordingly.

(def ^:private beta-reference
  (delay (edn/read-string (slurp (io/resource "special/beta_reference.edn")))))

(t/deftest beta
  (t/is (m/nan? (sut/beta ##NaN 1.0)))
  (t/is (m/nan? (sut/beta 1.0 ##NaN)))
  (t/testing "positive domain"
    (t/is (m/delta-eq 1.0 (sut/beta 1.0 1.0)))
    (t/is (m/delta-eq m/PI (sut/beta 0.5 0.5) 1.0e-14))
    (t/is (m/delta-eq (m// 1.0 12.0) (sut/beta 2.0 3.0) 1.0e-14)))
  (t/testing "negative arguments: removable singularity (p+q non-positive integer,
              p and q themselves finite) evaluates to 0.0"
    (t/is (m/zero? (sut/beta -0.5 -0.5)))
    (t/is (m/zero? (sut/beta -1.5 -2.5)))
    (t/is (m/zero? (sut/beta -0.25 -3.75)))
    (t/is (m/zero? (sut/beta -10.5 9.5))))
  (t/testing "negative arguments: true pole (p or q itself a non-positive integer) is NaN"
    (t/is (m/nan? (sut/beta -1.0 2.0)))
    (t/is (m/nan? (sut/beta 2.0 -1.0)))
    (t/is (m/nan? (sut/beta -2.0 -3.0)))
    (t/is (m/nan? (sut/beta -1.0 1.0))))
  (t/testing "symmetry: beta(p,q) = beta(q,p)"
    (let [ps (get-in @beta-reference [:beta :p])
          qs (get-in @beta-reference [:beta :q])]
      (t/is (v/edelta-eq (mapv sut/beta ps qs) (mapv sut/beta qs ps) 1.0e-9 1.0e-10))))
  (t/testing "recurrence: beta(p+1,q) = beta(p,q) * p/(p+q)"
    (doseq [[p q] (map vector (get-in @beta-reference [:beta :p]) (get-in @beta-reference [:beta :q]))]
      (t/is (m/delta-eq (sut/beta (m/inc p) q) (m/* (sut/beta p q) (m// p (m/+ p q))) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath (negative/mixed arguments)"
    (t/is (check2 sut/beta {:order (get-in @beta-reference [:beta :p])
                            :arg (get-in @beta-reference [:beta :q])
                            :ref (get-in @beta-reference [:beta :ref])} 1.0e-10))))

;; `log-beta` extends over the same domain as `beta` (see its docstring): it
;; is `(log (beta p q))`, so it reuses the same `beta` reference grid, split
;; by the sign of the reference `beta` value (log undefined for negative
;; beta -> NaN).
(t/deftest log-beta
  (t/is (m/nan? (sut/log-beta ##NaN 1.0)))
  (t/is (m/nan? (sut/log-beta 1.0 ##NaN)))
  (t/testing "positive domain"
    (t/is (m/delta-eq (m/log (m// 1.0 12.0)) (sut/log-beta 2.0 3.0) 1.0e-14))
    (t/testing "independent identity: log-beta(p,q) = log-gamma(p)+log-gamma(q)-log-gamma(p+q)"
      (doseq [[p q] [[2.0 3.0] [0.5 7.3] [100.0 50.0] [1.0 1.0] [0.1 0.2]]]
        (t/is (m/delta-eq (sut/log-beta p q)
                          (m/- (m/+ (sut/log-gamma p) (sut/log-gamma q)) (sut/log-gamma (m/+ p q)))
                          1.0e-12 1.0e-12)))))
  (t/testing "negative arguments: removable singularity (p+q non-positive integer,
              p and q themselves finite) evaluates to -Infinity (log of beta's 0.0)"
    (t/is (m/neg-inf? (sut/log-beta -0.5 -0.5)))
    (t/is (m/neg-inf? (sut/log-beta -1.5 -2.5)))
    (t/is (m/neg-inf? (sut/log-beta -0.25 -3.75)))
    (t/is (m/neg-inf? (sut/log-beta -10.5 9.5))))
  (t/testing "negative arguments: true pole (p or q itself a non-positive integer) is NaN"
    (t/is (m/nan? (sut/log-beta -1.0 2.0)))
    (t/is (m/nan? (sut/log-beta 2.0 -1.0)))
    (t/is (m/nan? (sut/log-beta -2.0 -3.0)))
    (t/is (m/nan? (sut/log-beta -1.0 1.0))))
  (t/testing "negative arguments: beta(p,q) < 0 has no real logarithm -> NaN"
    (t/is (m/nan? (sut/log-beta -0.5 2.0))) ;; beta(-0.5,2.0) = -4.0
    (t/is (m/nan? (sut/log-beta 3.2 -3.7)))) ;; negative beta value, see beta-reference
  (t/testing "matches beta over the same reference grid"
    (let [ps (get-in @beta-reference [:beta :p])
          qs (get-in @beta-reference [:beta :q])
          refs (get-in @beta-reference [:beta :ref])]
      (doseq [[p q r] (map vector ps qs refs)]
        (if (m/pos? r)
          (t/is (m/delta-eq (m/log r) (sut/log-beta p q) 1.0e-9 1.0e-10))
          (t/is (m/nan? (sut/log-beta p q))))))))

(t/deftest regularized-beta
  (t/is (m/nan? (sut/regularized-beta ##NaN 1.0 2.0)))
  (t/is (m/nan? (sut/regularized-beta 0.5 ##NaN 2.0)))
  (t/testing "positive domain"
    (t/is (m/zero? (sut/regularized-beta 0.0 2.0 3.0)))
    (t/is (m/one? (sut/regularized-beta 1.0 2.0 3.0)))
    (t/is (m/delta-eq 0.4 (sut/regularized-beta 0.4 1.0 1.0) 1.0e-14)))
  (t/testing "negative arguments: true (non-removable) pole when a+b is a non-positive integer"
    (t/is (m/nan? (sut/regularized-beta 0.3 -0.5 -1.5))))
  (t/testing "reflection: I_x(a,b) = 1 - I_1-x(b,a)"
    (let [xs (get-in @beta-reference [:incbeta :x])
          as (get-in @beta-reference [:incbeta :a])
          bs (get-in @beta-reference [:incbeta :b])]
      (doseq [[x a b] (map vector xs as bs)]
        (t/is (m/delta-eq (sut/regularized-beta x a b)
                          (m/- 1.0 (sut/regularized-beta (m/- 1.0 x) b a)) 1.0e-6 1.0e-6)))))
  (t/testing "vs mpmath (negative/mixed arguments)"
    (let [{:keys [x a b reg]} (:incbeta @beta-reference)]
      (t/is (v/edelta-eq (mapv sut/regularized-beta x a b) reg 1.0e-9 1.0e-8)))))

(t/deftest incomplete-beta
  (t/is (m/nan? (sut/incomplete-beta ##NaN 1.0 2.0)))
  (t/is (m/nan? (sut/incomplete-beta 0.5 ##NaN 2.0)))
  (t/testing "relation to regularized-beta and beta: incomplete-beta = regularized-beta * beta(a,b)"
    (let [xs (get-in @beta-reference [:incbeta :x])
          as (get-in @beta-reference [:incbeta :a])
          bs (get-in @beta-reference [:incbeta :b])]
      (t/is (v/edelta-eq (mapv sut/incomplete-beta xs as bs)
                         (map (fn [x a b] (m/* (sut/regularized-beta x a b) (sut/beta a b))) xs as bs)
                         1.0e-6 1.0e-6))))
  (t/testing "negative arguments: stays finite where regularized-beta is singular (a+b non-positive integer)"
    (t/is (m/delta-eq -0.24939187455541445 (sut/incomplete-beta 0.3 -0.5 -1.5) 1.0e-9)))
  (t/testing "vs mpmath (negative/mixed arguments)"
    (let [{ib-ref :inc :keys [x a b]} (:incbeta @beta-reference)]
      (t/is (v/edelta-eq (mapv sut/incomplete-beta x a b) ib-ref 1.0e-9 1.0e-8)))))

;;

;; Abramowitz and Stegun p.511
(t/deftest kummers-m
  (t/are [a b x res] (m/delta-eq res (sut/kummers-M a b x))
    0.3 0.2 -0.1 0.8578490
    -0.1 0.2 0.1 (* 0.8578490 (m/pow m/E 0.1))
    17 16 1 2.8881744
    -1 16 -1 1.0625
    -1.3 0.2 0.1 0.3582123
    -1.3 1.2 0.1 0.8924108
    -0.3 1.2 0.1 0.9745952
    1 1 1 m/E
    2 2 2 (m/exp 2)
    0.3 0.4 0.5 (/ 1.724128 (/ 0.7 0.6))
    ;; p.533
    -1 1 9 -8
    -1 0.6 9 -14
    0 1 9 1
    1 1 9 (m/exp 9))
  (t/are [a b x res acc] (m/delta-eq res (sut/kummers-M a b x) acc)
    0.9 0.1 10 1227235 1
    -52.5 0.1 1 -16.34 0.2)
  (t/is (m/pos-inf? (sut/kummers-M 1 0 1)))
  (t/is (m/neg-inf? (sut/kummers-M -1 0 1)))
  (t/is (m/pos-inf? (sut/kummers-M -1 0 -1)))
  (t/is (m/neg-inf? (sut/kummers-M 1 0 -1)))
  (t/is (m/nan? (sut/kummers-M -1 -1 2))) ;; example 4
  (t/testing "neg b and pos a is nan" (t/is (m/nan? (sut/kummers-M 10 -10 2))))
  (t/testing "neg b and neg a and a<b is nan"(t/is (m/nan? (sut/kummers-M -20 -10 2))))) 

(t/deftest whittaker-m
  (t/is (m/delta-eq 1.10622 (sut/whittaker-M 0 -0.4 1) 1.0e-5)))

(t/deftest besselk
  (t/are [order ress] (v/delta-eq ress (mapv (partial sut/bessel-K-half-odd order) [0.5 1 1.33 2.5 5]))
    1 [1.075047603 0.461068504 0.287423621 0.065065943 0.003776613]
    3 [3.225142810 0.922137009 0.503531608 0.091092320 0.004531936]
    5 [20.425904466  3.227479531  1.423209202  0.174376728  0.006495775]
    7 [207.48418748  17.05953466   5.85394214   0.43984578   0.01102771]
    9 [2.925204529e+03 1.226442222e+02 3.223343101e+01 1.405944900e+00 2.193457048e-02]))

;;

(t/deftest airy
  (t/is (check1 sut/airy-Ai (blk :airy :ai-small) 1.0e-9))
  (t/is (check1 sut/airy-Bi (blk :airy :bi-small) 1.0e-9))
  (t/is (check1 sut/airy-Ai' (blk :airy :ai'-small) 1.0e-9))
  (t/is (check1 sut/airy-Bi' (blk :airy :bi'-small) 1.0e-9))
  (t/is (check1 sut/airy-Ai (blk :airy :ai-large) 1.0e-9))
  (t/is (check1 sut/airy-Bi (blk :airy :bi-large) 1.0e-9))
  (t/is (check1 sut/airy-Ai' (blk :airy :ai'-large) 1.0e-9))
  (t/is (check1 sut/airy-Bi' (blk :airy :bi'-large) 1.0e-9)))

;;

(t/deftest lambert-w
  (doseq [x (repeatedly 1000 #(r/drand -1.0 500.0))]
    (t/is (v/delta-eq (sut/lambert-W (m/* x (m/exp x))) x)))
  (doseq [x (repeatedly 1000 #(r/drand -500.0 -1.0))]
    (t/is (v/delta-eq (sut/lambert-W-1(m/* x (m/exp x))) x))))

;;

;; owens_t

(def hvec [0.0625, 6.5, 7.0, 4.78125, 2.0, 1.0, 0.0625, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.125, 0.125, 0.125, 0.125, 0.0078125
         , 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0625, 0.5, 0.9, 2.5, 7.33, 0.6, 1.6, 2.33, 2.33])
(def avec [0.25, 0.4375, 0.96875, 0.0625, 0.5, 0.9999975, 0.999999125, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 10, 100
         , 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999
         , 0.99999])
(def cvec [0.0389119302347013668966224771378, 2.00057730485083154100907167685e-11, 6.399062719389853083219914429e-13
         , 1.06329748046874638058307112826e-7, 0.00862507798552150713113488319155, 0.0667418089782285927715589822405
         , 0.1246894855262192
         , 0.04306469112078537, 0.06674188216570097, 0.0784681869930841, 0.0792995047488726, 0.06448860284750375, 0.1066710629614485
         , 0.1415806036539784, 0.1510840430760184, 0.07134663382271778, 0.1201285306350883, 0.1666128410939293, 0.1847501847929859
         , 0.07317273327500386, 0.1237630544953746, 0.1737438887583106, 0.1951190307092811, 0.07378938035365545
         , 0.1249951430754052, 0.1761984774738108, 0.1987772386442824, 0.2340886964802671, 0.2479460829231492
         , 0.1246895548850743676554299881345328280176736760739893903915691894
         , 0.1066710629614484543187382775527753945264849005582264731161129477
         , 0.0750909978020473015080760056431386286348318447478899039422181015
         , 0.0030855526911589942124216949767707430484322201889568086810922629
         , 5.7538182971139187466647478665179637676531179007295252996453e-14, 0.0995191725772188724714794470740785702586033387786949658229016920
         , 0.0258981646643923680014142514442989928165349517076730515952020227
         , 0.0049025023268168675126146823752680242063832053551244071400100690
         , 0.0049024988349089450612896251009169062698683918433614542387524648])

(t/deftest owens_t
  (doseq [[h a res] (map vector hvec avec cvec)]
    (t/is (m/delta-eq res (sut/owens-t h a)))))
