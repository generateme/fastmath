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

(defn- check3
  "Compare `(f a b x)` against a three-argument reference block."
  ([f blk] (check3 f blk 1.0e-10))
  ([f blk rel] (let [{:keys [a b x ref]} blk]
                (v/edelta-eq (mapv f a b x) ref ABS rel))))

(defn- check4
  "Compare `(f a b c x)` against a four-argument reference block."
  ([f blk] (check4 f blk 1.0e-10))
  ([f blk rel] (let [{:keys [a b c x ref]} blk]
                (v/edelta-eq (mapv f a b c x) ref ABS rel))))

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

;; `zeta`/`eta`/`dirichlet-beta`/`xi` reference values below were computed
;; independently with Python `mpmath` (`mpmath.zeta`, `mpmath.altzeta`;
;; `dirichlet-beta`/`xi` via their standard defining formulas using mpmath's
;; own gamma/zeta) and are stored in `test/resources/special/zeta_reference.edn`.
;; See `utils/fastmath/dev/generate_zeta_reference.py` for the generator.
;;
;; KNOWN LIMITATION (see also `zeta`'s docstring): the 2-arity Hurwitz
;; `(zeta s z)` becomes unreliable for sufficiently negative `s` (roughly
;; `s < -5` for `z ≲ 2`, more for larger `z`) - its `:zeta2` reference grid
;; below is restricted to a domain verified to be accurate. The 1-arity
;; Riemann `(zeta s)` has no such restriction (`:zeta1` covers a wide range).

(def ^:private zeta-reference
  (delay (edn/read-string (slurp (io/resource "special/zeta_reference.edn")))))

(t/deftest zeta
  (t/testing "1-arity (Riemann): edge cases"
    (t/is (m/nan? (sut/zeta ##NaN)))
    (t/is (m/nan? (sut/zeta ##-Inf)))
    (t/is (m/nan? (sut/zeta 1.0))) ;; pole
    (t/is (m/one? (sut/zeta ##Inf))))
  (t/testing "1-arity: closed-form special values"
    (t/is (m/delta-eq -0.5 (sut/zeta 0.0)))
    (t/is (m/delta-eq (m// (m/pow m/PI 2.0) 6.0) (sut/zeta 2.0) 1.0e-14))
    (t/is (m/delta-eq (m// (m/pow m/PI 4.0) 90.0) (sut/zeta 4.0) 1.0e-13))
    (t/is (m/delta-eq (m// -1.0 12.0) (sut/zeta -1.0) 1.0e-14))
    (t/is (m/delta-eq (m// 1.0 120.0) (sut/zeta -3.0) 1.0e-14))
    (t/is (m/delta-eq (m// -1.0 252.0) (sut/zeta -5.0) 1.0e-14)))
  (t/testing "1-arity: trivial zeros at negative even integers, exact"
    (doseq [n [2 4 20 50 100 200 1000 10000]]
      (t/is (m/zero? (sut/zeta (m/- (double n)))))))
  (t/testing "vs mpmath, wide domain (positive and negative, both fixed via
              the log-space reflection rewrite)"
    (t/is (check1 sut/zeta {:arg (get-in @zeta-reference [:zeta1 :arg])
                            :ref (get-in @zeta-reference [:zeta1 :ref])} 1.0e-10)))
  (t/testing "2-arity (Hurwitz): shortcuts and edge cases"
    (t/is (m/delta-eq (sut/zeta 3.3) (sut/zeta 3.3 0.0) 1.0e-14))
    (t/is (m/delta-eq (sut/zeta 3.3) (sut/zeta 3.3 1.0) 1.0e-14))
    (t/is (m/delta-eq (sut/trigamma 2.5) (sut/zeta 2.0 2.5) 1.0e-14))
    (t/is (m/nan? (sut/zeta ##NaN 1.5)))
    (t/is (m/nan? (sut/zeta 2.0 ##NaN)))
    (t/is (m/nan? (sut/zeta ##-Inf 1.5)))
    (t/is (m/nan? (sut/zeta ##Inf -0.5)))
    (t/is (m/zero? (sut/zeta ##Inf 2.0)))
    (t/is (m/pos-inf? (sut/zeta ##Inf 0.5))))
  (t/testing "2-arity: closed forms via Bernoulli polynomials"
    (doseq [z [0.3 1.7 5.2]]
      (t/is (m/delta-eq (m/- 0.5 z) (sut/zeta 0.0 z) 1.0e-13))
      (t/is (m/delta-eq (m/* -0.5 (m/+ (m/* z z) (m/- z) (m// 1.0 6.0))) (sut/zeta -1.0 z) 1.0e-12))))
  (t/testing "2-arity: recurrence zeta(s,z+1) = zeta(s,z) - z^-s (safe domain)"
    (doseq [s [0.3 2.0 3.0 10.5] z [0.3 1.5 3.7]]
      (t/is (m/delta-eq (sut/zeta s (m/inc z)) (m/- (sut/zeta s z) (m/pow z (m/- s))) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath, restricted safe domain"
    (t/is (check2 sut/zeta {:order (get-in @zeta-reference [:zeta2 :s])
                            :arg (get-in @zeta-reference [:zeta2 :z])
                            :ref (get-in @zeta-reference [:zeta2 :ref])} 1.0e-9))))

(t/deftest eta
  (t/is (m/nan? (sut/eta ##NaN)))
  (t/is (m/nan? (sut/eta ##-Inf)))
  (t/is (m/one? (sut/eta ##Inf)))
  (t/testing "closed-form special values"
    (t/is (m/delta-eq 0.5 (sut/eta 0.0)))
    (t/is (m/delta-eq m/LN2 (sut/eta 1.0) 1.0e-14))
    (t/is (m/delta-eq (m// (m/pow m/PI 2.0) 12.0) (sut/eta 2.0) 1.0e-14))
    (t/is (m/delta-eq 0.25 (sut/eta -1.0) 1.0e-14)))
  (t/testing "relation to zeta: eta(s) = (1 - 2^(1-s)) * zeta(s)"
    (doseq [s (get-in @zeta-reference [:zeta1 :arg])
            :when (not (m/one? s))]
      (t/is (m/delta-eq (sut/eta s) (m/* (m/- 1.0 (m/pow 2.0 (m/- 1.0 s))) (sut/zeta s)) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath (altzeta), wide domain"
    (t/is (check1 sut/eta {:arg (get-in @zeta-reference [:zeta1 :arg])
                           :ref (get-in @zeta-reference [:eta :ref])} 1.0e-10))))

(t/deftest dirichlet-beta
  (t/is (m/nan? (sut/dirichlet-beta ##NaN)))
  (t/testing "closed-form special values"
    (t/is (m/delta-eq (m// m/PI 4.0) (sut/dirichlet-beta 1.0) 1.0e-14))
    (t/is (m/delta-eq 0.5 (sut/dirichlet-beta 0.0)))
    (t/is (m/near-zero? (sut/dirichlet-beta -1.0) 1.0e-9))
    (t/is (m/delta-eq -0.5 (sut/dirichlet-beta -2.0) 1.0e-9))
    (t/is (m/near-zero? (sut/dirichlet-beta -3.0) 1.0e-9))
    (t/is (m/delta-eq 2.5 (sut/dirichlet-beta -4.0) 1.0e-9)))
  (t/testing "converges to 1.0 as x -> +Inf, including for large finite x where
              the general formula would otherwise hit a 0 * (Inf - Inf) failure
              mode (`4^-x` underflows to 0.0 before `4^x` in zeta(x,0.25) overflows)"
    (t/is (m/one? (sut/dirichlet-beta 100.0)))
    (t/is (m/one? (sut/dirichlet-beta 600.0)))
    (t/is (m/one? (sut/dirichlet-beta 1000.0)))
    (t/is (m/one? (sut/dirichlet-beta ##Inf))))
  (t/testing "-Inf has no well-defined limit (oscillating growth), stays NaN"
    (t/is (m/nan? (sut/dirichlet-beta ##-Inf))))
  (t/testing "vs mpmath, wide positive / restricted negative domain"
    (t/is (check1 sut/dirichlet-beta {:arg (get-in @zeta-reference [:dbeta :arg])
                                      :ref (get-in @zeta-reference [:dbeta :ref])} 1.0e-9))))

(t/deftest xi
  (t/is (m/nan? (sut/xi ##NaN)))
  (t/is (m/pos-inf? (sut/xi ##Inf)))
  (t/is (m/pos-inf? (sut/xi ##-Inf)))
  (t/testing "closed-form special values"
    (t/is (m/delta-eq 0.5 (sut/xi 0.0)))
    (t/is (m/delta-eq 0.5 (sut/xi 1.0))))
  (t/testing "functional equation: xi(s) = xi(1-s), wide domain (fixed via the
              log-space rewrite)"
    (doseq [s (get-in @zeta-reference [:xi :arg])]
      (t/is (m/delta-eq (sut/xi s) (sut/xi (m/- 1.0 s)) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath, wide domain"
    (t/is (check1 sut/xi {:arg (get-in @zeta-reference [:xi :arg])
                          :ref (get-in @zeta-reference [:xi :ref])} 1.0e-10))))

;; Reference values for `polygamma` were computed independently with
;; `mpmath.polygamma` and are stored in
;; `test/resources/special/polygamma_reference.edn`. See
;; `utils/fastmath/dev/generate_polygamma_reference.py` for the generator,
;; and for why `x<=0` is restricted to small orders/moderate `|x|` there
;; (a known limitation of the internal `cotderiv` helper, see also
;; `polygamma`'s own docstring).

(def ^:private polygamma-reference
  (delay (edn/read-string (slurp (io/resource "special/polygamma_reference.edn")))))

(t/deftest polygamma
  (t/testing "negative order -> NaN"
    (t/is (m/nan? (sut/polygamma -1 1.0)))
    (t/is (m/nan? (sut/polygamma -5 2.0))))
  (t/testing "order 0/1 delegate to digamma/trigamma"
    (doseq [x [0.5 1.5 2.5 10.0 -0.5 -2.5]]
      (t/is (m/delta-eq (sut/polygamma 0 x) (sut/digamma x) 1.0e-12))
      (t/is (m/delta-eq (sut/polygamma 1 x) (sut/trigamma x) 1.0e-12))))
  (t/testing "NaN/Inf x"
    (t/is (m/nan? (sut/polygamma 3 ##NaN)))
    (t/is (m/zero? (sut/polygamma 3 ##Inf)))
    (t/is (m/nan? (sut/polygamma 3 ##-Inf))))
  (t/testing "pole at x=0 (and -0.0): +-Inf, sign depending on order parity"
    (t/is (m/neg-inf? (sut/polygamma 2 0.0)))
    (t/is (m/pos-inf? (sut/polygamma 3 0.0)))
    (t/is (m/neg-inf? (sut/polygamma 2 -0.0)))
    (t/is (m/pos-inf? (sut/polygamma 3 -0.0))))
  (t/testing "poles at negative integers: large magnitude (a floating-point
              approximation of a true pole via `cotderiv`'s cot(pi*z),
              not exactly +-Inf)"
    (t/is (> (Math/abs (sut/polygamma 2 -1.0)) 1.0e40))
    (t/is (> (Math/abs (sut/polygamma 3 -2.0)) 1.0e40)))
  (t/testing "closed form at x=1: psi^(m)(1) = (-1)^(m+1) * m! * zeta(m+1)"
    (doseq [m [2 3 4 5 10 50 100 150]]
      (let [expected (* (if (odd? m) 1.0 -1.0)
                         (Math/exp (sut/log-gamma (inc m)))
                         (sut/zeta (inc m)))]
        (t/is (m/delta-eq (sut/polygamma m 1.0) expected 1.0e-9 1.0e-9)))))
  (t/testing "recurrence: psi^(m)(x+1) = psi^(m)(x) + (-1)^m * m!/x^(m+1)"
    (doseq [m (range 2 11) x [0.7 1.3 2.5 5.5]]
      (let [rhs (+ (sut/polygamma m x)
                   (* (if (even? m) 1.0 -1.0) (Math/exp (sut/log-gamma (inc m)))
                      (Math/pow x (- (inc m)))))]
        (t/is (m/delta-eq (sut/polygamma m (inc x)) rhs 1.0e-8 1.0e-8)))))
  (t/testing "sign pattern for x>0: negative for even order, positive for odd"
    (doseq [m (range 2 11) x [0.5 1.0 2.5 10.0]]
      (let [v (sut/polygamma m x)]
        (t/is (if (even? m) (neg? v) (pos? v))))))
  (t/testing "decays to 0 for large x"
    (t/is (m/delta-eq 0.0 (sut/polygamma 2 1e12) 1.0e-9))
    (t/is (m/delta-eq 0.0 (sut/polygamma 3 1e12) 1.0e-9)))
  (t/testing "vs mpmath, x>0, wide order range (up to m=300, after the
              log-space Gamma-overflow fix)"
    (t/is (check2 sut/polygamma (:pos @polygamma-reference) 1.0e-8)))
  (t/testing "vs mpmath, x<=0, small orders (see namespace comment above for
              why larger orders/|x| are excluded here)"
    (t/is (check2 sut/polygamma (:neg @polygamma-reference) 1.0e-5))))

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
  ;; a=b=-1 is a removable 0/0 coincidence in the series (both Pochhammer
  ;; symbols vanish at the same term), not a genuine pole: the mathematically
  ;; correct value is the truncated series 1 + x (confirmed against mpmath),
  ;; not NaN nor exp(x) -- see fastmath.special-test/hypergeometric-1F1.
  (t/is (m/delta-eq 3.0 (sut/kummers-M -1 -1 2)))
  (t/testing "neg b and pos a is nan" (t/is (m/nan? (sut/kummers-M 10 -10 2))))
  (t/testing "neg b and neg a and a<b is nan"(t/is (m/nan? (sut/kummers-M -20 -10 2))))) 

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
    (t/is (v/delta-eq (sut/lambert-W-1 (m/* x (m/exp x))) x))))

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

;; `Si` reference values below were computed independently with Python
;; `mpmath` (`mpmath.si`, mpmath's own name for the standard sine integral -
;; not to be confused with fastmath's lowercase `si`, which is the shifted
;; quantity `Si(x) - pi/2` and is tested below purely via that relation) and
;; are stored in `test/resources/special/si_reference.edn`. See
;; `utils/fastmath/dev/generate_si_reference.py` for the generator. The grid
;; is dense near 0 and around each branch switch of the piecewise
;; implementation (`x = +-6`, `x = +-12`, and the `x*x` double-overflow
;; threshold around `+-1.34e154`) and sparse-but-wide out to `+-1e155`.
;;
;; A handful of points spanning all branches were additionally cross-checked
;; independently against R (`pracma::Si`) and Julia (`SpecialFunctions.sinint`);
;; both agree with the `mpmath` values here to double precision (see
;; `si-cross-check` below).

(def ^:private si-reference
  (delay (edn/read-string (slurp (io/resource "special/si_reference.edn")))))

(t/deftest Si
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Si ##NaN)))
    (t/is (m/delta-eq m/HALF_PI (sut/Si ##Inf)))
    (t/is (m/delta-eq m/-HALF_PI (sut/Si ##-Inf)))
    (t/is (m/zero? (sut/Si 0.0))))
  (t/testing "oddness: Si(-x) = -Si(x), over the full reference domain"
    (doseq [x (get-in @si-reference [:Si :arg])]
      (t/is (m/delta-eq (sut/Si (m/- x)) (m/- (sut/Si x)) 1.0e-12 1.0e-12))))
  (t/testing "continuity across the internal branch switches (x=6, x=12)"
    (doseq [edge [6.0 12.0]]
      (t/is (m/delta-eq (sut/Si (m/- edge 1.0e-9)) (sut/Si (m/+ edge 1.0e-9)) 1.0e-8 1.0e-8))
      (t/is (m/delta-eq (sut/Si (m/- (m/- edge) 1.0e-9)) (sut/Si (m/+ (m/- edge) 1.0e-9)) 1.0e-8 1.0e-8))))
  (t/testing "vs mpmath, wide domain across all branches (small/medium/large/overflow)"
    (t/is (check1 sut/Si {:arg (get-in @si-reference [:Si :arg])
                          :ref (get-in @si-reference [:Si :ref])} 1.0e-10)))
  (t/testing "vs R (pracma::Si) and Julia (SpecialFunctions.sinint), spot check across all branches"
    (let [si-cross-check {:arg [0.5 1.0 3.0 6.0 6.0001 8.0 12.0 12.0001 50.0 1000.0 1.0e8]
                          :ref [0.4931074180430667 0.946083070367183 1.8486525279994683
                                1.4246875512805066 1.4246828951944845 1.5741868217069421
                                1.5049712415263734 1.5049667704556322 1.551617072485936
                                1.5702331219687713 1.5707963304287473]}]
      (t/is (check1 sut/Si si-cross-check 1.0e-10)))))

(t/deftest si
  (t/is (m/nan? (sut/si ##NaN)))
  (t/testing "relation to Si: si(x) = Si(x) - pi/2, over the full reference domain"
    (doseq [x (get-in @si-reference [:Si :arg])]
      (t/is (m/delta-eq (sut/si x) (m/- (sut/Si x) m/HALF_PI) 1.0e-14 1.0e-14))))
  (t/testing "value at 0 and at the infinities (exact, since it's a plain
              double subtraction of Si's own exact limits)"
    (t/is (m/delta-eq m/-HALF_PI (sut/si 0.0)))
    (t/is (m/zero? (sut/si ##Inf)))
    (t/is (m/delta-eq (m/- m/PI) (sut/si ##-Inf))))
  (t/testing "converges to 0 as x -> +Inf, to -pi as x -> -Inf (large finite x)"
    (t/is (m/delta-eq 0.0 (sut/si 1.0e8) 1.0e-8))
    (t/is (m/delta-eq (m/- m/PI) (sut/si -1.0e8) 1.0e-8)))
  (t/testing "not an odd function (unlike Si): si(-x) != -si(x) in general"
    (t/is (not (m/delta-eq (sut/si -1.0) (m/- (sut/si 1.0)) 1.0e-6 1.0e-6)))))

;; `Ci` reference values below (`x >= 0` only, matching fastmath's own domain
;; restriction) were computed independently with Python `mpmath`
;; (`mpmath.ci`). `Cin` reference values (all real `x`, including negative and
;; zero) were computed independently from mpmath's own Euler-Mascheroni
;; constant, log and ci -- NOT by calling fastmath's own formula. Both are
;; stored in `test/resources/special/ci_reference.edn`. See
;; `utils/fastmath/dev/generate_ci_reference.py` for the generator. As with
;; `Si`, the grids are dense near the branch switches of the piecewise
;; implementation (`x = 3, 6, 12` and the `x*x` double-overflow threshold
;; around `1.34e154`).
;;
;; A handful of points spanning all branches were additionally cross-checked
;; independently against R (`pracma::Ci`) and Julia
;; (`SpecialFunctions.cosint`); both agree with the `mpmath` values here to
;; double precision (see `ci-cross-check` below).

(def ^:private ci-reference
  (delay (edn/read-string (slurp (io/resource "special/ci_reference.edn")))))

(t/deftest Ci
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Ci ##NaN)))
    (t/is (m/neg-inf? (sut/Ci 0.0)))
    (t/is (m/zero? (sut/Ci ##Inf)))
    (t/is (thrown? AssertionError (sut/Ci -1.0)))
    (t/is (thrown? AssertionError (sut/Ci ##-Inf))))
  (t/testing "continuity across the internal branch switches (x=3, x=6, x=12)"
    (doseq [edge [3.0 6.0 12.0]]
      (t/is (m/delta-eq (sut/Ci (m/- edge 1.0e-9)) (sut/Ci (m/+ edge 1.0e-9)) 1.0e-8 1.0e-8))))
  (t/testing "vs mpmath, wide domain across all branches (small/medium/large/overflow)"
    (t/is (check1 sut/Ci {:arg (get-in @ci-reference [:Ci :arg])
                          :ref (get-in @ci-reference [:Ci :ref])} 1.0e-9)))
  (t/testing "vs R (pracma::Ci) and Julia (SpecialFunctions.cosint), spot check across all branches"
    (let [ci-cross-check {:arg [0.5 1.0 3.0 3.0001 6.0 6.0001 8.0 12.0 12.0001 50.0 1000.0]
                          :ref [-0.1777840788066129 0.33740392290096816 0.11962978600800032
                                0.1195967865729574 -0.06805724389324713 -0.06804124095567483
                                0.12243388253200956 -0.04978000688411367 -0.04977297457353216
                                -0.005628386324116306 0.0008263155110906822]}]
      (t/is (check1 sut/Ci ci-cross-check 1.0e-9)))))

(t/deftest Cin
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Cin ##NaN)))
    (t/is (m/zero? (sut/Cin 0.0)))
    (t/is (m/pos-inf? (sut/Cin ##Inf)))
    (t/is (m/pos-inf? (sut/Cin ##-Inf))))
  (t/testing "evenness: Cin(-x) = Cin(x), over the full reference domain
              (unlike Ci, Cin is defined and finite for negative x too)"
    (doseq [x (get-in @ci-reference [:Cin :arg])]
      (t/is (m/delta-eq (sut/Cin (m/- x)) (sut/Cin x) 1.0e-9 1.0e-9))))
  (t/testing "relation to Ci: Cin(x) = gamma + log(x) - Ci(x), for x > 0"
    (doseq [x (get-in @ci-reference [:Ci :arg])
            :when (m/pos? x)]
      (t/is (m/delta-eq (sut/Cin x) (m/- (m/+ m/GAMMA (m/log x)) (sut/Ci x)) 1.0e-9 1.0e-9))))
  (t/testing "non-negative everywhere"
    (doseq [x (get-in @ci-reference [:Cin :arg])]
      (t/is (m/>= (sut/Cin x) -1.0e-9))))
  (t/testing "vs mpmath, wide domain (all real x, incl. negative and 0)"
    (t/is (check1 sut/Cin {:arg (get-in @ci-reference [:Cin :arg])
                           :ref (get-in @ci-reference [:Cin :ref])} 1.0e-8))))

;; Reference values for `E0`/`E1`/`Ein`/`En`/`Ei`/`li`/`Li` below were computed
;; independently with Python `mpmath` (`mpmath.expint` for `E0`/`E1`/`En`,
;; `mpmath.ei` for `Ei`, `mpmath.li` for `li`/`Li`; `Ein` computed directly
;; from mpmath's own Euler-Mascheroni constant/log/expint for `x > 0` and from
;; the alternating entire power series for `x <= 0`, NOT via fastmath's own
;; formula) and are stored in `test/resources/special/ei_reference.edn`. See
;; `utils/fastmath/dev/generate_ei_reference.py` for the generator.
;;
;; Spot checks for `E1`, `En` (orders 2 and 5) and `Ei` were additionally
;; cross-validated independently against R (`expint` package's `expint_E1`,
;; `expint_En`, `expint_Ei`) and Julia (`SpecialFunctions.expint`/`expinti`);
;; all three agree with the `mpmath` values to double precision (see
;; `e-cross-check` below). `li`/`Li`/`Ein`/`E0` have no direct equivalents in
;; those two libraries, so they are validated against `mpmath` only.
;;
;; Three genuine numerical bugs were found and fixed while building this
;; reference data (all independently confirmed against `mpmath`):
;;  - `E0`: `exp(-x)` alone overflowed to `##Inf` for `x` roughly in
;;    `(-716.35, -709.78)`, even though the true ratio `exp(-x)/x` still fits
;;    in double range there; fixed via a log-space fallback.
;;  - `Ein`: was `##NaN` for every negative `x` (delegated entirely to `E1`,
;;    which is undefined there); now extended to `x >= -2.15` via its own
;;    entire power series (still `##NaN` below that, a documented remaining
;;    limitation).
;;  - `Ei`: a "near its real zero" branch (`|x - 0.3725| < 0.3`) using a
;;    13-term Taylor series was inaccurate by up to ~1% inside its own claimed
;;    domain (insufficient terms for that radius); removed entirely in favor
;;    of the already-present, already-accurate general Taylor branches, which
;;    turned out to cover that whole region correctly (including exactly at
;;    the zero) with no replacement needed. Also fixed the same premature-
;;    overflow issue as `E0` (the `x > 710.0 -> Inf` shortcut was replaced
;;    with an exact log-space computation, extending the accurate range to
;;    the true double-precision limit).
;;  - `En`: `(En n 0.0)` for `0 < n < 1` returned a finite (wrong) value
;;    instead of the mathematically correct `##Inf` (the defining integral
;;    itself diverges there, same as for `n <= 0`).

(def ^:private ei-reference
  (delay (edn/read-string (slurp (io/resource "special/ei_reference.edn")))))

(t/deftest E0
  (t/testing "edge cases"
    (t/is (m/nan? (sut/E0 ##NaN)))
    (t/is (m/pos-inf? (sut/E0 0.0)))
    (t/is (m/zero? (sut/E0 ##Inf)))
    (t/is (m/neg-inf? (sut/E0 ##-Inf))))
  (t/testing "relation to En: E0(x) = En(0, x)"
    (doseq [x (get-in @ei-reference [:E0 :arg])]
      (t/is (m/delta-eq (sut/E0 x) (sut/En 0.0 x) 1.0e-9 1.0e-9))))
  (t/testing "no premature overflow for very negative x where the true ratio
              exp(-x)/x still fits in double range"
    (t/is (m/valid-double? (sut/E0 -710.0)))
    (t/is (m/valid-double? (sut/E0 -716.0)))
    (t/is (m/neg-inf? (sut/E0 -717.0))))
  (t/testing "vs mpmath, wide domain (all real x)"
    (t/is (check1 sut/E0 {:arg (get-in @ei-reference [:E0 :arg])
                          :ref (get-in @ei-reference [:E0 :ref])} 1.0e-9))))

(t/deftest E1
  (t/testing "edge cases"
    (t/is (m/nan? (sut/E1 ##NaN)))
    (t/is (m/pos-inf? (sut/E1 0.0)))
    (t/is (m/zero? (sut/E1 ##Inf)))
    (t/is (m/nan? (sut/E1 -1.0))))
  (t/testing "relation to En: E1(x) = En(1, x)"
    (doseq [x (get-in @ei-reference [:E1 :arg])]
      (t/is (m/delta-eq (sut/E1 x) (sut/En 1.0 x) 1.0e-9 1.0e-9))))
  (t/testing "continuity across the internal branch switches (looser near the
              small-x end, where E1's derivative ~ -1/x is steep enough that
              even a 2e-9 step in x amplifies to a much larger step in value)"
    (doseq [edge [0.0044 0.053 0.6 2.15 4.0 10.0 20.0 200.0]]
      (t/is (m/delta-eq (sut/E1 (m/- edge 1.0e-9)) (sut/E1 (m/+ edge 1.0e-9)) 1.0e-6 1.0e-6))))
  (t/testing "vs mpmath, wide domain (x >= 0)"
    (t/is (check1 sut/E1 {:arg (get-in @ei-reference [:E1 :arg])
                          :ref (get-in @ei-reference [:E1 :ref])} 1.0e-9)))
  (t/testing "vs R (expint::expint_E1) and Julia (SpecialFunctions.expint), spot check"
    (let [e1-cross-check {:arg [0.001 0.1 1.0 2.15 4.0 10.0 20.0 200.0 500.0]
                          :ref [6.331539364136149 1.8229239584193906 0.21938393439552029
                                0.0398034560480185 0.0037793524098489058 4.1569689296853246e-6
                                9.8355252906498815e-11 6.8852261063076359e-90 1.4220767822536385e-220]}]
      (t/is (check1 sut/E1 e1-cross-check 1.0e-9)))))

(t/deftest Ein
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Ein ##NaN)))
    (t/is (m/zero? (sut/Ein 0.0)))
    (t/is (m/pos-inf? (sut/Ein ##Inf)))
    (t/is (m/nan? (sut/Ein ##-Inf)))
    (t/is (m/nan? (sut/Ein -5.0)))
    (t/is (m/valid-double? (sut/Ein -2.15))))
  (t/testing "relation to E1: Ein(x) = E1(x) + log(x) + gamma, for x > 0"
    (doseq [x (get-in @ei-reference [:E1 :arg])]
      (t/is (m/delta-eq (sut/Ein x) (m/+ (sut/E1 x) (m/log x) m/GAMMA) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath, wide domain (x >= -2.15)"
    (t/is (check1 sut/Ein {:arg (get-in @ei-reference [:Ein :arg])
                           :ref (get-in @ei-reference [:Ein :ref])} 1.0e-8))))

(t/deftest En
  (t/testing "edge cases"
    (t/is (m/nan? (sut/En ##NaN 1.0)))
    (t/is (m/nan? (sut/En 1.0 ##NaN)))
    (t/is (m/pos-inf? (sut/En 1.0 0.0)))
    (t/is (m/pos-inf? (sut/En 0.5 0.0)))
    (t/is (m/pos-inf? (sut/En -0.5 0.0)))
    (t/is (m/delta-eq 2.0 (sut/En 1.5 0.0)))
    (t/is (m/delta-eq 1.0 (sut/En 2.0 0.0))))
  (t/testing "negative x: NaN for positive integer or fractional order, real for non-positive integer order"
    (t/is (m/nan? (sut/En 1.0 -1.0)))
    (t/is (m/nan? (sut/En 2.0 -1.0)))
    (t/is (m/nan? (sut/En 0.5 -1.0)))
    (t/is (m/valid-double? (sut/En 0.0 -1.0)))
    (t/is (m/valid-double? (sut/En -1.0 -1.0)))
    (t/is (m/valid-double? (sut/En -2.0 -1.0))))
  (t/testing "vs mpmath, grid of integer and fractional orders"
    (t/is (check2 sut/En {:order (get-in @ei-reference [:En :order])
                          :arg (get-in @ei-reference [:En :arg])
                          :ref (get-in @ei-reference [:En :ref])} 1.0e-6)))
  (t/testing "vs R (expint::expint_En) and Julia (SpecialFunctions.expint), spot check for orders 2 and 5"
    (let [xs [0.001 0.1 1.0 2.15 4.0 10.0 20.0 200.0 500.0]
          en2-cross-check {:order (repeat 9 2.0) :arg xs
                           :ref [0.99266896046923891 0.72254502219402039 0.14849550677592205
                                 0.0309067272702572 0.0031982292493385554 3.8302404656316078e-6
                                 9.4048564308581459e-11 6.8513054752104106e-90 1.419249547309342e-220]}
          en5-cross-check {:order (repeat 9 5.0) :arg xs
                           :ref [0.24966691650035003 0.21901595224028042 0.070454237461720332
                                 0.017887851745349667 0.0021555113535254578 3.089728914253682e-6
                                 8.3071305994177173e-11 6.7515104974077901e-90 1.4108347621366221e-220]}]
      (t/is (check2 sut/En en2-cross-check 1.0e-9))
      (t/is (check2 sut/En en5-cross-check 1.0e-9)))))

(t/deftest Ei
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Ei ##NaN)))
    (t/is (m/neg-inf? (sut/Ei 0.0)))
    (t/is (m/pos-inf? (sut/Ei ##Inf)))
    (t/is (m/delta-eq -0.0 (sut/Ei ##-Inf))))
  (t/testing "relation to E1: Ei(x) = -E1(-x), for x < 0"
    (doseq [x (get-in @ei-reference [:E1 :arg])
            :when (m/pos? x)]
      (t/is (m/delta-eq (sut/Ei (m/- x)) (m/- (sut/E1 x)) 1.0e-9 1.0e-9))))
  (t/testing "no premature overflow for very large x where the true value still fits in double range"
    (t/is (m/valid-double? (sut/Ei 710.0)))
    (t/is (m/valid-double? (sut/Ei 716.0)))
    (t/is (m/pos-inf? (sut/Ei 717.0))))
  (t/testing "accurate through the whole 'near its real zero' region (a former
              inaccurate special-cased branch, now removed)"
    (doseq [x [0.073 0.083 0.093 0.103 0.3725074107813666 0.5 0.6]]
      (t/is (m/valid-double? (sut/Ei x)))))
  (t/testing "vs mpmath, wide domain (all real x != 0)"
    (t/is (check1 sut/Ei {:arg (get-in @ei-reference [:Ei :arg])
                          :ref (get-in @ei-reference [:Ei :ref])} 1.0e-8)))
  (t/testing "vs R (expint::expint_Ei) and Julia (SpecialFunctions.expinti), spot check"
    (let [ei-cross-check {:arg [-5.0 -1.0 -0.1 0.1 0.6 1.0 2.15 4.0 10.0 20.0 200.0 500.0]
                          :ref [-0.0011482955912753257 -0.21938393439552029 -1.8229239584193906
                                -1.6228128139692763 0.76988128993735938 1.8951178163559366
                                5.5302550090155798 19.630874470056217 2492.2289762418782
                                25615652.664056588 3.6312352331593567e84 2.8128213978862945e214]}]
      (t/is (check1 sut/Ei ei-cross-check 1.0e-8)))))

(t/deftest li
  (t/testing "edge cases"
    (t/is (m/nan? (sut/li ##NaN)))
    (t/is (m/nan? (sut/li -1.0)))
    (t/is (m/delta-eq -0.0 (sut/li 0.0)))
    (t/is (m/neg-inf? (sut/li 1.0))))
  (t/testing "relation to Ei: li(x) = Ei(ln(x))"
    (doseq [x (get-in @ei-reference [:li :arg])]
      (t/is (m/delta-eq (sut/li x) (sut/Ei (m/log x)) 1.0e-9 1.0e-9))))
  (t/testing "single positive real zero at the Ramanujan-Soldner constant"
    (t/is (m/near-zero? (sut/li 1.45136923488338105028396848589) 1.0e-9)))
  (t/testing "vs mpmath, wide domain (x > 0, x != 1)"
    (t/is (check1 sut/li {:arg (get-in @ei-reference [:li :arg])
                          :ref (get-in @ei-reference [:li :ref])} 1.0e-9))))

(t/deftest Li
  (t/testing "edge cases"
    (t/is (m/nan? (sut/Li ##NaN)))
    (t/is (m/nan? (sut/Li -1.0)))
    (t/is (m/neg-inf? (sut/Li 1.0)))
    (t/is (m/near-zero? (sut/Li 2.0) 1.0e-9)))
  (t/testing "relation to li: Li(x) = li(x) - li(2)"
    (doseq [x (get-in @ei-reference [:li :arg])]
      (t/is (m/delta-eq (sut/Li x) (m/- (sut/li x) (sut/li 2.0)) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath, wide domain (x > 0, x != 1)"
    (t/is (check1 sut/Li {:arg (get-in @ei-reference [:Li :arg])
                          :ref (get-in @ei-reference [:Li :ref])} 1.0e-9))))

;; Reference values for the low-order generalized hypergeometric functions
;; below were computed once from `mpmath` (`mpmath.hyper` for `0F0`/`1F0`,
;; `mpmath.hyp0f1` for `0F1`) and are stored in
;; `test/resources/special/hyp_low_order_reference.edn`. See
;; `utils/fastmath/dev/generate_hyp_low_order_reference.py` for the
;; generator.
;;
;; Spot checks were additionally cross-validated against R's `hypergeo`
;; package (`genhypergeo`), which agrees with `mpmath` to double precision on
;; the (non-oscillatory) points checked below; it was not used for `1F0`
;; points outside its own domain of support (it returns `NA` for some
;; negative-order/negative-argument combinations that `mpmath` and fastmath
;; both handle fine), nor for `0F1` at large negative `x` (oscillatory
;; Bessel-`J` regime, where its series summation loses several digits
;; relative to `mpmath`).
;;
;; No numerical bugs were found in `0F0`, `1F0` or `0F1`; both are thin,
;; already-correct wrappers (`exp`, `(1-x)^(-a)`, and the Bessel-`I`/`J`
;; relation, respectively).

(def ^:private hyp-reference
  (delay (edn/read-string (slurp (io/resource "special/hyp_low_order_reference.edn")))))

(t/deftest hypergeometric-0F0
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-0F0 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-0F0 0.0)))
    (t/is (m/pos-inf? (sut/hypergeometric-0F0 ##Inf)))
    (t/is (m/zero? (sut/hypergeometric-0F0 ##-Inf))))
  (t/testing "identity: 0F0(x) = exp(x)"
    (doseq [x (get-in @hyp-reference [:0F0 :arg])]
      (t/is (m/delta-eq (sut/hypergeometric-0F0 x) (m/exp x) 1.0e-12 1.0e-12))))
  (t/testing "vs mpmath, wide domain incl. near double over/underflow boundaries"
    (t/is (check1 sut/hypergeometric-0F0 {:arg (get-in @hyp-reference [:0F0 :arg])
                                          :ref (get-in @hyp-reference [:0F0 :ref])} 1.0e-12)))
  (t/testing "vs R (hypergeo::genhypergeo), spot check"
    (t/is (check1 sut/hypergeometric-0F0 {:arg [-5.0 -1.0 0.0 1.0 5.0 10.0]
                                          :ref [0.0067379469990846378 0.36787944117144245
                                                1.0 2.7182818284590455
                                                148.4131591025766 22026.465794806714]}
                  1.0e-9))))

(t/deftest hypergeometric-1F0
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-1F0 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-1F0 0.5 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-1F0 0.0 0.5)))
    (t/is (m/one? (sut/hypergeometric-1F0 0.0 1.0)))
    (t/is (m/pos-inf? (sut/hypergeometric-1F0 2.0 1.0)))
    (t/is (m/zero? (sut/hypergeometric-1F0 -2.0 1.0))))
  (t/testing "negative-integer order: series truncates to a polynomial, real for x >= 1 too"
    (t/is (m/delta-eq -1.0 (sut/hypergeometric-1F0 -1.0 2.0)))
    (t/is (m/delta-eq 0.0 (sut/hypergeometric-1F0 -2.0 1.0)))
    (t/is (m/delta-eq 1.0 (sut/hypergeometric-1F0 -2.0 2.0))))
  (t/testing "x >= 1, non-integer order: complex off the real line, NaN"
    (t/is (m/nan? (sut/hypergeometric-1F0 0.5 2.0)))
    (t/is (m/nan? (sut/hypergeometric-1F0 -2.5 1.5))))
  (t/testing "x >= 1, integer (non-negative) order: still on the real line"
    (t/is (m/delta-eq 1.0 (sut/hypergeometric-1F0 2.0 2.0))))
  (t/testing "identity: 1F0(a, x) = (1-x)^(-a)"
    (doseq [[a x] (map vector (get-in @hyp-reference [:1F0 :order])
                       (get-in @hyp-reference [:1F0 :arg]))]
      (t/is (m/delta-eq (sut/hypergeometric-1F0 a x) (m/pow (m/- 1.0 x) (m/- a)) 1.0e-9 1.0e-9))))
  (t/testing "vs mpmath, grid of orders (incl. negative integers, x < 1 elsewhere)"
    (t/is (check2 sut/hypergeometric-1F0 {:order (get-in @hyp-reference [:1F0 :order])
                                          :arg (get-in @hyp-reference [:1F0 :arg])
                                          :ref (get-in @hyp-reference [:1F0 :ref])} 1.0e-9)))
  (t/testing "vs R (hypergeo::genhypergeo), spot check"
    (t/is (check2 sut/hypergeometric-1F0 {:order [2.0] :arg [0.5] :ref [3.9999999999999991]} 1.0e-9))
    (t/is (check2 sut/hypergeometric-1F0 {:order [-2.0] :arg [1.0] :ref [0.0]} 1.0e-9))))

(t/deftest hypergeometric-0F1
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-0F1 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F1 1.5 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-0F1 1.5 0.0)))
    (t/is (m/one? (sut/hypergeometric-0F1 0.0 0.0)))
    (t/is (m/one? (sut/hypergeometric-0F1 -3.0 0.0))))
  (t/testing "poles at non-positive integer order, for any nonzero x"
    (t/is (m/nan? (sut/hypergeometric-0F1 0.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F1 -1.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F1 -2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F1 -3.0 -0.5))))
  (t/testing "negative non-integer order: finite, real"
    (t/is (m/valid-double? (sut/hypergeometric-0F1 -1.5 0.5)))
    (t/is (m/valid-double? (sut/hypergeometric-0F1 -2.5 -0.5))))
  (t/testing "relation to Bessel I / J for x > 0 / x < 0"
    (doseq [[a x] (map vector (get-in @hyp-reference [:0F1 :order])
                       (get-in @hyp-reference [:0F1 :arg]))
            :when (and (m/pos? a) (not (m/zero? x)))]
      (let [xx (m/* 2.0 (m/sqrt (m/abs x)))
            a- (m/dec a)
            expected (if (m/pos? x)
                       (m// (m/* (sut/bessel-I a- xx) (sut/gamma a)) (m/pow (m/* 0.5 xx) a-))
                       (m// (m/* (sut/bessel-J a- xx) (sut/gamma a)) (m/pow (m/* 0.5 xx) a-)))]
        (t/is (m/delta-eq (sut/hypergeometric-0F1 a x) expected 1.0e-9 1.0e-9))))
    ;; x = 0.0 is handled directly (1.0), the identity above needs x != 0
    (t/is (m/one? (sut/hypergeometric-0F1 2.5 0.0))))
  (t/testing "vs mpmath, grid of orders (positive and negative non-integer), wide x range"
    (t/is (check2 sut/hypergeometric-0F1 {:order (get-in @hyp-reference [:0F1 :order])
                                          :arg (get-in @hyp-reference [:0F1 :arg])
                                          :ref (get-in @hyp-reference [:0F1 :ref])} 1.0e-7)))
  (t/testing "vs R (hypergeo::genhypergeo), spot check"
    (t/is (check2 sut/hypergeometric-0F1 {:order [1.5 -1.5 0.5]
                                          :arg [0.5 0.5 10.0]
                                          :ref [1.3682988720085907 0.89370818366377003
                                                279.05568512996319]}
                  1.0e-9))))

;; Reference values for `hypergeometric-1F1` (Kummer's confluent
;; hypergeometric function M) below were computed once from `mpmath`
;; (`mpmath.hyp1f1`) and are stored in
;; `test/resources/special/hyp1f1_reference.edn`, split into blocks
;; `:generic` (a, b non-integer/positive, b never a non-positive integer),
;; `:overflow` (a probe for large |x| against genuine double-precision
;; overflow), `:equal` (a = b, both the exp(x) sub-case and the
;; truncated-series sub-case), `:neg-int-a` (a a non-positive integer, b
;; generic) and `:neg-int-both-safe` (both a, b non-positive integers,
;; a >= b, a != b). See `utils/fastmath/dev/generate_hyp1f1_reference.py`
;; for the generator.
;;
;; Spot checks were additionally cross-validated against R's `gsl` package
;; (`hyperg_1F1`), which agrees with `mpmath` to double precision on every
;; point checked below EXCEPT the a = b = (negative integer) coincidence,
;; where `gsl` returns `exp(x)` unconditionally (e.g. `hyperg_1F1(-3,-3,2)`
;; = `exp(2)` = 7.389..., but the correct value, confirmed independently by
;; hand-derivation and by mpmath's own direct series summation, is the
;; truncated series 1+2+2+4/3 = 6.333...) -- so `gsl`, like fastmath's own
;; `kummers-M` before this session's fix, does not handle that specific
;; edge case correctly; it is therefore not used as a cross-check there.
;;
;; Two genuine numerical bugs were found and fixed in `kummers-M` while
;; building this reference data (both confirmed against `mpmath`):
;;  - `a = b` on a negative integer previously returned `##NaN` (for
;;    `a = b = -1` specifically) or `exp(x)` (for every other negative
;;    integer `a = b`); both were wrong. The series has a removable `0/0`
;;    coincidence at the term where both Pochhammer symbols vanish
;;    together; the correct value is the truncated exponential series
;;    `sum_{n=0}^{|a|} x^n/n!`.
;;  - The positive-`x` series loop returned `##NaN` once its running sum
;;    first overflowed to `##Inf` (`Inf - Inf` in the convergence-check
;;    arithmetic), instead of `##Inf`, for `x` large enough that the true
;;    value genuinely exceeds double range (e.g. `M(2,3,1000)`, true value
;;    ~3.94e431).

(def ^:private hyp1f1-reference
  (delay (edn/read-string (slurp (io/resource "special/hyp1f1_reference.edn")))))

(t/deftest hypergeometric-1F1
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-1F1 ##NaN 2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-1F1 2.0 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-1F1 2.0 2.0 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-1F1 2.0 3.0 0.0)))
    (t/is (m/one? (sut/hypergeometric-1F1 0.0 3.0 0.5)))
    (t/is (m/one? (sut/hypergeometric-1F1 0.0 -3.0 0.5))))
  (t/testing "poles: b non-positive integer, a positive or a < b (both integers)"
    (t/is (m/nan? (sut/hypergeometric-1F1 2.0 -3.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-1F1 -5.0 -3.0 0.5))))
  (t/testing "b = 0.0, a != 0.0: signed infinity"
    (t/is (m/pos-inf? (sut/hypergeometric-1F1 1.0 0.0 1.0)))
    (t/is (m/neg-inf? (sut/hypergeometric-1F1 -1.0 0.0 1.0)))
    (t/is (m/pos-inf? (sut/hypergeometric-1F1 -1.0 0.0 -1.0)))
    (t/is (m/neg-inf? (sut/hypergeometric-1F1 1.0 0.0 -1.0))))
  (t/testing "a = b on a negative integer: removable 0/0 coincidence, truncated
              exponential series (fixed this session), not NaN nor exp(x)"
    (t/is (m/delta-eq 1.5 (sut/hypergeometric-1F1 -1.0 -1.0 0.5)))
    (t/is (m/delta-eq 1.625 (sut/hypergeometric-1F1 -2.0 -2.0 0.5)))
    (t/is (m/delta-eq 1.6458333333333333 (sut/hypergeometric-1F1 -3.0 -3.0 0.5))))
  (t/testing "a = b elsewhere (zero, positive, non-integer): exp(x)"
    (t/is (m/delta-eq (m/exp 0.5) (sut/hypergeometric-1F1 0.5 0.5 0.5)))
    (t/is (m/delta-eq (m/exp 0.5) (sut/hypergeometric-1F1 2.0 2.0 0.5))))
  (t/testing "no premature NaN for large positive x where the true value
              overflows double range (fixed this session)"
    (t/is (m/pos-inf? (sut/hypergeometric-1F1 2.0 3.0 1000.0)))
    (t/is (m/valid-double? (sut/hypergeometric-1F1 2.0 3.0 700.0))))
  (t/testing "vs mpmath, generic (a, b non-integer/positive, b never a non-positive integer)"
    (t/is (check3 sut/hypergeometric-1F1 {:a (get-in @hyp1f1-reference [:generic :a])
                                          :b (get-in @hyp1f1-reference [:generic :b])
                                          :x (get-in @hyp1f1-reference [:generic :x])
                                          :ref (get-in @hyp1f1-reference [:generic :ref])} 1.0e-6)))
  (t/testing "vs mpmath, large |x| overflow probe"
    (t/is (check3 sut/hypergeometric-1F1 {:a (get-in @hyp1f1-reference [:overflow :a])
                                          :b (get-in @hyp1f1-reference [:overflow :b])
                                          :x (get-in @hyp1f1-reference [:overflow :x])
                                          :ref (get-in @hyp1f1-reference [:overflow :ref])} 1.0e-6)))
  (t/testing "vs mpmath, a = b (both sub-cases)"
    (t/is (check3 sut/hypergeometric-1F1 {:a (get-in @hyp1f1-reference [:equal :a])
                                          :b (get-in @hyp1f1-reference [:equal :b])
                                          :x (get-in @hyp1f1-reference [:equal :x])
                                          :ref (get-in @hyp1f1-reference [:equal :ref])} 1.0e-9)))
  (t/testing "vs mpmath, a non-positive integer, b generic (always terminates)"
    (t/is (check3 sut/hypergeometric-1F1 {:a (get-in @hyp1f1-reference [:neg-int-a :a])
                                          :b (get-in @hyp1f1-reference [:neg-int-a :b])
                                          :x (get-in @hyp1f1-reference [:neg-int-a :x])
                                          :ref (get-in @hyp1f1-reference [:neg-int-a :ref])} 1.0e-9)))
  (t/testing "vs mpmath, both a, b non-positive integers, a >= b, a != b (safe termination)"
    (t/is (check3 sut/hypergeometric-1F1 {:a (get-in @hyp1f1-reference [:neg-int-both-safe :a])
                                          :b (get-in @hyp1f1-reference [:neg-int-both-safe :b])
                                          :x (get-in @hyp1f1-reference [:neg-int-both-safe :x])
                                          :ref (get-in @hyp1f1-reference [:neg-int-both-safe :ref])} 1.0e-9)))
  (t/testing "vs R (gsl::hyperg_1F1), spot check"
    (t/is (check3 sut/hypergeometric-1F1 {:a [2.0 0.5 5.0 -1.0 -2.0 1.5]
                                          :b [3.0 3.0 0.5 2.0 -5.0 2.5]
                                          :x [5.0 -10.0 -1.0 3.0 3.0 0.5]
                                          :ref [47.572210912824517 0.44148780381255504
                                                -0.23303836033449427 -0.5 2.6500000000000004
                                                1.3612908263697017]}
                  1.0e-7))))

;; Reference values for `hypergeometric-0F2` below were computed once from
;; `mpmath` (`mpmath.hyper([], [a, b], x)`) and are stored in
;; `test/resources/special/hyp0f2_reference.edn`, split into blocks
;; `:generic` (full a, b grid incl. negative, |x| <= 200 -- see precision
;; note below), `:negative-x-wide` (a, b both positive only, |x| up to
;; 5000) and `:wide-positive` (full a, b grid, x positive only, up to 1e7).
;; See `utils/fastmath/dev/generate_hyp0f2_reference.py` for the generator.
;; Poles (a or b a non-positive integer, x != 0) are not in the bulk
;; reference (mpmath raises there, matching fastmath's NaN); covered by
;; dedicated edge-case assertions below instead.
;;
;; Spot checks were additionally cross-validated against R's `hypergeo`
;; package (`genhypergeo`), which agrees with `mpmath` to double precision
;; on every point checked.
;;
;; One genuine numerical bug was found and fixed in `hypergeometric-0F2`
;; while building this reference data: when `a` or `b` is a non-positive
;; integer (a genuine pole -- confirmed via mpmath raising "pole in
;; hypergeometric series"), it previously returned inconsistent, sometimes
;; silently wrong results depending on the sign of `x`: `##Inf`/`##-Inf` for
;; positive `x` (via `maclaurin-0F2`), but for negative `x` (via
;; `weniger-0F2`) it often returned a large but finite, WRONG number instead
;; (e.g. `(hypergeometric-0F2 -1.0 2.5 -0.5)` was `-1.02e17`, not a signal
;; of the pole at all). Fixed by adding an explicit pole check in
;; `hypergeometric-0F2` itself, ahead of both code paths.
;;
;; A separate, NOT fixed (by explicit user choice, documented as a known
;; limitation) numerical-stability issue was also found: `weniger-0F2`
;; (negative `x`) progressively loses accuracy as `|x|` grows, and how fast
;; depends heavily on `a`/`b`. When `a` and `b` are BOTH positive, it stays
;; accurate (~1e-10 relative or better) out to at least `|x| = 5000`, only
;; breaking down catastrophically (wrong order of magnitude, even wrong
;; sign) beyond roughly `|x| ~ 20000-50000` (e.g. `a=1.5, b=2.5`: `x=-20000`
;; is already ~1e-5 relative error; `x=-50000` gives a positive ~7.6e16
;; where the true value is a negative ~-4.0e17). When `a` or `b` is
;; negative (even far from any actual pole), the safe range shrinks
;; dramatically -- e.g. at `a=-1.7, b=-2.7`, `x=-5000` is already off by
;; more than 100% (wrong sign); `|x| <= 200` was empirically confirmed safe
;; (~1e-8 relative or better) across the whole a, b grid used here. This is
;; deliberately NOT tested beyond those confirmed-safe boundaries.

(def ^:private hyp0f2-reference
  (delay (edn/read-string (slurp (io/resource "special/hyp0f2_reference.edn")))))

(t/deftest hypergeometric-0F2
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-0F2 ##NaN 2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 2.0 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 2.0 2.0 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-0F2 2.0 3.0 0.0)))
    (t/is (m/one? (sut/hypergeometric-0F2 0.0 2.5 0.0)))
    (t/is (m/one? (sut/hypergeometric-0F2 -3.0 -5.0 0.0))))
  (t/testing "poles: a or b a non-positive integer, x != 0.0 (fixed this session:
              previously NaN/Inf/-Inf/silently-wrong-finite depending on sign of x)"
    (t/is (m/nan? (sut/hypergeometric-0F2 0.0 2.5 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 -1.0 2.5 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 -1.0 2.5 -0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 1.5 -2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 -1.0 -2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-0F2 -1.0 -2.0 -0.5))))
  (t/testing "no NaN gap at the genuine double-precision overflow boundary for large positive x"
    (t/is (m/valid-double? (sut/hypergeometric-0F2 1.5 2.5 1.0e7)))
    (t/is (m/pos-inf? (sut/hypergeometric-0F2 1.5 2.5 1.5e7)))
    (t/is (m/pos-inf? (sut/hypergeometric-0F2 1.5 2.5 1.0e9))))
  (t/testing "vs mpmath, generic (a, b non-integer or positive, |x| <= 200)"
    (t/is (check3 sut/hypergeometric-0F2 {:a (get-in @hyp0f2-reference [:generic :a])
                                          :b (get-in @hyp0f2-reference [:generic :b])
                                          :x (get-in @hyp0f2-reference [:generic :x])
                                          :ref (get-in @hyp0f2-reference [:generic :ref])} 1.0e-7)))
  (t/testing "vs mpmath, negative x extended range, a and b both positive"
    (t/is (check3 sut/hypergeometric-0F2 {:a (get-in @hyp0f2-reference [:negative-x-wide :a])
                                          :b (get-in @hyp0f2-reference [:negative-x-wide :b])
                                          :x (get-in @hyp0f2-reference [:negative-x-wide :x])
                                          :ref (get-in @hyp0f2-reference [:negative-x-wide :ref])} 1.0e-6)))
  (t/testing "vs mpmath, positive x extended range (up to 1e7)"
    (t/is (check3 sut/hypergeometric-0F2 {:a (get-in @hyp0f2-reference [:wide-positive :a])
                                          :b (get-in @hyp0f2-reference [:wide-positive :b])
                                          :x (get-in @hyp0f2-reference [:wide-positive :x])
                                          :ref (get-in @hyp0f2-reference [:wide-positive :ref])} 1.0e-9)))
  (t/testing "vs R (hypergeo::genhypergeo), spot check"
    (t/is (check3 sut/hypergeometric-0F2 {:a [1.5 1.5 0.5 1.5]
                                          :b [2.5 2.5 0.5 2.5]
                                          :x [0.5 -0.5 -100.0 100.0]
                                          :ref [1.1371833737326889 0.87043608108448023
                                                532.81415522956274 1089.2355638194263]}
                  1.0e-9))))

;; Reference values for `whittaker-M` below were computed once from mpmath
;; (`mpmath.whitm`) and are stored in
;; `test/resources/special/whittaker_m_reference.edn` (kappa, mu grid incl.
;; negative values, x in (0, 300]; see `utils/fastmath/dev/
;; generate_whittaker_m_reference.py` for the generator). No R equivalent
;; was found (`gsl`, `Bessel` packages checked; no Whittaker function
;; available), so this is mpmath-only.
;;
;; Two genuine numerical bugs, both confirmed against mpmath, were found and
;; fixed this session:
;;  - `kummers-M`'s `b = 0.0` branch was checked before its `x = 0.0`
;;    branch, so `kummers-M(a, 0.0, 0.0)` returned a signed `##Inf` instead
;;    of the mathematically correct `1.0` (mpmath: `hyp1f1(2,0,0) = 1.0`),
;;    inconsistent with the x=0-always-wins convention already established
;;    and tested for `hypergeometric-0F1`/`0F2`. Fixed by reordering.
;;  - `whittaker-M(kappa, -0.5, 0.0)` returned `##NaN` (from `mu+0.5 = 0.0`
;;    times `log(0) = ##-Inf` = `##NaN` in its own formula) instead of the
;;    correct limiting value `1.0` (`x^(mu+1/2) = x^0 = 1` identically at
;;    `mu = -0.5`, no actual singularity). Fixed with an explicit x=0,
;;    mu=-0.5 short-circuit.
;;
;; A separate, NOT fixed (by explicit user choice, documented as a known
;; limitation) numerical issue was also found: `whittaker-M` can incorrectly
;; return `##Inf` when `kummers-M`'s own intermediate value individually
;; overflows double range, even though the final SCALED product
;; (`kummers-M`'s huge value times the tiny `z^2` prefactor) would be a
;; valid finite double -- e.g. `whittaker-M(-3.0, -2.3, 700.0)` returns
;; `##Inf`, but the true value is `~9.3e159`, well within double range. A
;; real fix would need a log-space/rescaled variant of `kummers-M` (a
;; nontrivial addition, out of scope here). Confirmed via the wider probe
;; grid used while building this reference: unaffected for x <= 300 (the
;; domain used below), affected for roughly 7% of points at x in {700,
;; 1000} in that wider probe.

(def ^:private whittaker-m-reference
  (delay (edn/read-string (slurp (io/resource "special/whittaker_m_reference.edn")))))

(t/deftest whittaker-M
  (t/testing "edge cases"
    (t/is (m/nan? (sut/whittaker-M ##NaN 0.5 1.0)))
    (t/is (m/nan? (sut/whittaker-M 1.0 ##NaN 1.0)))
    (t/is (m/nan? (sut/whittaker-M 1.0 0.5 ##NaN)))
    (t/is (m/nan? (sut/whittaker-M 1.0 0.5 -1.0))))
  (t/testing "x -> 0+ limit depends on the sign of mu+0.5"
    (t/is (m/zero? (sut/whittaker-M 1.0 0.5 0.0)))
    (t/is (m/zero? (sut/whittaker-M 1.0 -0.2 0.0)))
    (t/is (m/pos-inf? (sut/whittaker-M 1.0 -0.8 0.0)))
    ;; mu = -0.5 exactly: fixed this session, was NaN, correct limit is 1.0
    ;; for any kappa
    (t/is (m/delta-eq 1.0 (sut/whittaker-M 1.0 -0.5 0.0)))
    (t/is (m/delta-eq 1.0 (sut/whittaker-M 2.0 -0.5 0.0)))
    (t/is (m/delta-eq 1.0 (sut/whittaker-M -3.5 -0.5 0.0))))
  (t/testing "known limitation (not fixed, documented): intermediate overflow
              of kummers-M's own value can make whittaker-M wrongly return
              ##Inf even where the true, properly-scaled value is finite"
    (t/is (m/pos-inf? (sut/whittaker-M -3.0 -2.3 700.0))) ;; true value ~9.3e159
    (t/is (m/valid-double? (sut/whittaker-M -3.0 -2.3 300.0))))
  (t/testing "vs mpmath, grid of kappa, mu (incl. negative), x in (0, 300]"
    (t/is (check3 sut/whittaker-M {:a (get-in @whittaker-m-reference [:kappa])
                                   :b (get-in @whittaker-m-reference [:mu])
                                   :x (get-in @whittaker-m-reference [:x])
                                   :ref (get-in @whittaker-m-reference [:ref])} 1.0e-9))))

;; Reference values for `hypergeometric-2F0` below were computed once from
;; mpmath (`mpmath.hyper([a, b], [], x)`) and are stored in
;; `test/resources/special/hyp2f0_reference.edn`, split into `:terminating`
;; (a or b a non-positive integer, both signs of x), `:generic-negative-safe`
;; (|x| <= 0.5) and `:generic-negative-moderate` (0.5 < |x| <= 5.0). See
;; `utils/fastmath/dev/generate_hyp2f0_reference.py` for the generator.
;;
;; `hypergeometric-2F0`'s underlying series diverges for every x != 0
;; (radius of convergence 0) unless a or b is a non-positive integer, in
;; which case it is an exact finite polynomial; otherwise it is necessarily
;; computed via resummation of the divergent series.
;;
;; One genuine numerical bug was found and fixed in the terminating case
;; this session: `weniger-2F0`'s general resummation algorithm occasionally
;; returned `##NaN` at specific coincidental (a, b, x), e.g.
;; `a=-2.0, b=1.0, x=-1.0` exactly (while x=-0.999 or -1.001 both correctly
;; gave values close to the true 5.0) -- despite the true value being a
;; perfectly well-defined finite polynomial there. Fixed by adding a direct
;; polynomial-evaluation branch in `hypergeometric-2F0` itself whenever a or
;; b is a non-positive integer (mirroring the analogous `hypergeometric-1F1`
;; fix earlier this session), which also improved precision there (from
;; occasional outright failure to ~1e-13 relative).
;;
;; A separate, NOT fixed (by explicit user choice, documented as a known
;; limitation) numerical-stability issue was also found in the generic
;; (non-terminating) case: for x < 0, `weniger-2F0`'s resummation agrees
;; with mpmath's to ~1e-15 relative at |x| <= 0.1, but the agreement
;; degrades progressively as |x| grows -- ~1e-6 at |x| <= 0.5, ~1e-2 at
;; |x| <= 5, over 100% (wrong order of magnitude) by |x| ~ 1000. This is a
;; genuine limitation of resumming this particular divergent series, not a
;; simple bug; this generator therefore does not extend the generic-case
;; reference domain past |x| = 5.0.
;;
;; Separately, for x > 0 in the generic case, mpmath's own resummation
;; returns a genuinely COMPLEX number (confirmed e.g. at
;; a=1.5, b=2.5, x=0.5), whose real part does not match fastmath's real
;; result there -- this is expected, not a bug: resumming a divergent
;; series is inherently branch-dependent once x crosses to the other side
;; of the series' natural boundary, and different resummation conventions
;; (Borel summation vs. the Pade/Weniger-type acceleration used here) can
;; legitimately disagree. x > 0 in the generic case is therefore not tested
;; against any external reference at all.

(def ^:private hyp2f0-reference
  (delay (edn/read-string (slurp (io/resource "special/hyp2f0_reference.edn")))))

(t/deftest hypergeometric-2F0
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-2F0 ##NaN 2.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-2F0 2.0 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-2F0 2.0 2.0 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-2F0 2.0 3.0 0.0)))
    (t/is (m/one? (sut/hypergeometric-2F0 -3.0 2.5 0.0))))
  (t/testing "terminating case (a or b a non-positive integer): exact finite
              polynomial, valid for any sign of x, incl. the specific point
              that used to be ##NaN before this session's fix"
    (t/is (m/delta-eq 5.0 (sut/hypergeometric-2F0 -2.0 1.0 -1.0)))
    (t/is (m/delta-eq 37.0 (sut/hypergeometric-2F0 3.0 -2.0 2.0)))
    (t/is (m/delta-eq 181.0 (sut/hypergeometric-2F0 -3.0 -2.0 5.0)))
    (t/is (m/one? (sut/hypergeometric-2F0 0.0 2.5 0.5))))
  (t/testing "vs mpmath, terminating case, both signs of x"
    (t/is (check3 sut/hypergeometric-2F0 {:a (get-in @hyp2f0-reference [:terminating :a])
                                          :b (get-in @hyp2f0-reference [:terminating :b])
                                          :x (get-in @hyp2f0-reference [:terminating :x])
                                          :ref (get-in @hyp2f0-reference [:terminating :ref])} 1.0e-9)))
  (t/testing "vs mpmath, generic case, x < 0, |x| <= 0.5 (high precision)"
    (t/is (check3 sut/hypergeometric-2F0 {:a (get-in @hyp2f0-reference [:generic-negative-safe :a])
                                          :b (get-in @hyp2f0-reference [:generic-negative-safe :b])
                                          :x (get-in @hyp2f0-reference [:generic-negative-safe :x])
                                          :ref (get-in @hyp2f0-reference [:generic-negative-safe :ref])} 1.0e-5)))
  (t/testing "vs mpmath, generic case, x < 0, 0.5 < |x| <= 5.0 (known-limitation, looser tolerance)"
    (t/is (check3 sut/hypergeometric-2F0 {:a (get-in @hyp2f0-reference [:generic-negative-moderate :a])
                                          :b (get-in @hyp2f0-reference [:generic-negative-moderate :b])
                                          :x (get-in @hyp2f0-reference [:generic-negative-moderate :x])
                                          :ref (get-in @hyp2f0-reference [:generic-negative-moderate :ref])} 5.0e-2))))

;; Reference values for `tricomis-U` below were computed once from mpmath
;; (`mpmath.hyperu`) and are stored in
;; `test/resources/special/tricomis_u_reference.edn` (a, b grid incl. the
;; a=b and a=b-1 internal special-cased coincidences, x >= 0.01 only -- see
;; below for why x=0 and very small x are excluded from the bulk
;; reference). See `utils/fastmath/dev/generate_tricomis_u_reference.py`
;; for the generator.
;;
;; x = 0.0 is NOT covered by the bulk reference: mpmath's own direct
;; `hyperu(a,b,0)` evaluation for b >= 1 was found UNRELIABLE (e.g. at
;; a=-0.5, b=3.0 it literally returns +inf, while evaluating at a sequence
;; of very small positive x clearly trends to -inf instead). The x=0
;; boundary values used below (in "edge cases") were instead derived from
;; those small-x limit trends by hand, and, for a a non-positive integer,
;; from an independently confirmed closed form.
;;
;; Three genuine numerical bugs, all confirmed against mpmath (via small-x
;; limit trends for the x=0 cases, since direct mpmath x=0 evaluation is
;; unreliable there, see above), were found and fixed this session:
;;  - `tricomis-U(a, b, 0.0)` for `b >= 1.0` unconditionally returned
;;    `##NaN`, even though `U` genuinely diverges there for generic `a`.
;;    Fixed to return a signed `##Inf`, with sign = sign(gamma(a)) (since
;;    gamma(b-1) is always positive for b>1, and the same sign rule was
;;    confirmed, separately, to also hold for the logarithmic divergence at
;;    b=1 exactly).
;;  - Digging further: when `a` is ALSO a non-positive integer (`gamma(a)`
;;    itself a pole), `U(a,b,0)` does NOT actually diverge -- confirmed via
;;    small-x limit trends to be finite, equal to `(-1)^n (b)_n` (`a=-n`,
;;    `(b)_n` the rising Pochhammer symbol) for both `b<1` and `b>=1`,
;;    unifying with (and superseding) the pre-existing `b<1` gamma-ratio
;;    formula for that specific sub-case.
;;  - The general branch (`a != b`, `a != b-1`) called `weniger-2F0`
;;    DIRECTLY instead of through the already-fixed `hypergeometric-2F0`
;;    wrapper, so it was still exposed to that function's own
;;    (already-fixed-elsewhere) terminating-case `##NaN` glitch, e.g.
;;    `tricomis-U(1.0, 5.0, 2.0)` was `##NaN` instead of `2.375`. Fixed by
;;    calling `hypergeometric-2F0` instead.
;;  - The `a = b` branch computed `exp(x) * upper-incomplete-gamma(1-a,x)`;
;;    `exp(x)` alone overflows to `##Inf` for `x > ~710`, giving `##NaN` (or,
;;    right at the boundary, a silently wrong `##Inf`) even though the true
;;    value is a tiny, perfectly finite double there (confirmed:
;;    `U(3,3,1000)` true value `~9.97e-10`). Found via `whittaker-W`'s own
;;    testing (its x grid reaches x=1000, exposing this for the a=b case).
;;    Fixed by falling back to the general asymptotic-series formula
;;    (`x^(-a) * hypergeometric-2F0(a,1,-1/x)`) whenever `exp(x)` overflows
;;    -- confirmed to already agree with the `a=b` formula to ~13
;;    significant digits well before the overflow boundary (e.g. at
;;    `x=300`).
;;
;; A separate, NOT a new bug: the general branch computes
;; `x^(-a) * hypergeometric-2F0(a, 1+a-b, -1/x)`; for small x, `-1/x` is a
;; large-magnitude negative argument, surfacing `hypergeometric-2F0`'s own
;; already-documented precision limitation there (see that function's test
;; section). Observed max relative error ~1.9e-3 for x >= 0.01 in the grid
;; used here; smaller x is not covered by the bulk reference.
;;
;; x < 0.0 is outside `tricomis-U`'s validated domain (see docstring) and
;; is not tested against any external reference.

(def ^:private tricomis-u-reference
  (delay (edn/read-string (slurp (io/resource "special/tricomis_u_reference.edn")))))

(t/deftest tricomis-U
  (t/testing "edge cases"
    (t/is (m/nan? (sut/tricomis-U ##NaN 2.0 0.5)))
    (t/is (m/nan? (sut/tricomis-U 2.0 ##NaN 0.5)))
    (t/is (m/nan? (sut/tricomis-U 2.0 2.0 ##NaN)))
    (t/is (m/nan? (sut/tricomis-U 1.5 2.5 -0.5))))
  (t/testing "x = 0.0, b < 1.0: closed form, pre-existing, unaffected"
    (t/is (m/delta-eq 1.772453850905516 (sut/tricomis-U 1.5 0.5 0.0))))
  (t/testing "x = 0.0, b >= 1.0, a not a non-positive integer: genuinely
              diverges, signed by gamma(a) (fixed this session, was ##NaN)"
    (t/is (m/pos-inf? (sut/tricomis-U 1.5 2.5 0.0)))
    (t/is (m/pos-inf? (sut/tricomis-U 1.5 1.0 0.0)))
    (t/is (m/neg-inf? (sut/tricomis-U -0.5 3.0 0.0)))
    (t/is (m/neg-inf? (sut/tricomis-U -0.5 1.0 0.0))))
  (t/testing "x = 0.0, a a non-positive integer: divergence is removable,
              finite (-1)^n (b)_n (fixed this session, was ##NaN or, with
              the first fix alone, wrongly ##Inf)"
    (t/is (m/delta-eq -24.0 (sut/tricomis-U -3.0 2.0 0.0)))
    (t/is (m/delta-eq 12.0 (sut/tricomis-U -2.0 3.0 0.0)))
    (t/is (m/delta-eq 59.0625 (sut/tricomis-U -4.0 1.5 0.0)))
    (t/is (m/delta-eq -5.0 (sut/tricomis-U -1.0 5.0 0.0)))
    (t/is (m/delta-eq 0.75 (sut/tricomis-U -2.0 0.5 0.0))))
  (t/testing "general branch used to route around hypergeometric-2F0's fix
              (fixed this session): now consistent"
    (t/is (m/delta-eq 2.375 (sut/tricomis-U 1.0 5.0 2.0))))
  (t/testing "a = b branch: no premature NaN/Inf for large x where exp(x)
              alone overflows but the true (tiny) value fits in double
              range (fixed this session)"
    (t/is (m/valid-double? (sut/tricomis-U 3.0 3.0 700.0)))
    (t/is (m/delta-eq 2.78225113631986E-9 (sut/tricomis-U 3.0 3.0 710.0) 1.0e-9 1.0e-9))
    (t/is (m/delta-eq 9.970119403575002E-10 (sut/tricomis-U 3.0 3.0 1000.0) 1.0e-9 1.0e-9)))
  (t/testing "vs mpmath, wide (a, b) grid incl. a=b and a=b-1 coincidences, x >= 0.01"
    (t/is (check3 sut/tricomis-U {:a (get-in @tricomis-u-reference [:a])
                                  :b (get-in @tricomis-u-reference [:b])
                                  :x (get-in @tricomis-u-reference [:x])
                                  :ref (get-in @tricomis-u-reference [:ref])} 5.0e-3))))

;; Reference values for `whittaker-W` below were computed once from mpmath
;; (`mpmath.whitw`) and are stored in
;; `test/resources/special/whittaker_w_reference.edn` (kappa, mu grid incl.
;; negative values, x in [0.1, 1000]; see `utils/fastmath/dev/
;; generate_whittaker_w_reference.py` for the generator).
;;
;; x = 0.0 is a separate boundary case, not covered by the bulk reference
;; (its correct value depends discontinuously on the sign of mu+0.5,
;; mirroring whittaker-M); covered by dedicated edge-case assertions
;; instead.
;;
;; One genuine numerical bug, directly analogous to the whittaker-M fix
;; earlier this session, was found and fixed here: whittaker-W(kappa, -0.5,
;; 0.0) returned ##NaN (from mu+0.5=0.0 times log(0)=##-Inf = ##NaN in its
;; own formula) instead of the correct limiting value, tricomis-U's own
;; value at that point (x^(mu+1/2) = x^0 = 1 identically at mu=-0.5, no
;; actual singularity). Fixed with an explicit x=0, mu=-0.5 short-circuit,
;; confirmed via small-x mpmath trends to match tricomis-U's already-correct
;; x=0 handling exactly (e.g. kappa=1 -> 0.0, kappa=2 -> 0.0, kappa=-3.5 ->
;; ~0.08597), for three different kappa spanning tricomis-U's different x=0
;; sub-branches.
;;
;; A SEPARATE genuine bug was found and fixed in `tricomis-U` itself (not
;; whittaker-W) while testing this -- see that function's own test section
;; for details (an exp(x) overflow in its a=b branch, surfaced by
;; whittaker-W's x grid reaching x=1000).
;;
;; Unlike whittaker-M, no intermediate-overflow limitation was found here:
;; tricomis-U is the decaying/bounded solution, so it does not tend to
;; overflow at large x the way kummers-M does. A DIFFERENT, NOT fixed
;; (documented) limitation was found instead: for small x, tricomis-U's
;; general branch feeds a large-magnitude negative argument (-1/x) to
;; hypergeometric-2F0, surfacing that function's own already-documented
;; precision limitation. Through the mu, kappa -> a, b transform used here,
;; this shows up at a larger x threshold than for tricomis-U tested
;; directly (the mu, kappa grid used here reaches larger |a|, |b|, e.g.
;; kappa=5.0, mu=-3.0 gives a=-7.5, b=-5.0): x=0.01 still shows ~95%
;; relative error for some points in this grid, dropping to ~1e-7 by
;; x=0.1 (the domain used below) and to machine precision by x=0.2.

(def ^:private whittaker-w-reference
  (delay (edn/read-string (slurp (io/resource "special/whittaker_w_reference.edn")))))

(t/deftest whittaker-W
  (t/testing "edge cases"
    (t/is (m/nan? (sut/whittaker-W ##NaN 0.5 1.0)))
    (t/is (m/nan? (sut/whittaker-W 1.0 ##NaN 1.0)))
    (t/is (m/nan? (sut/whittaker-W 1.0 0.5 ##NaN)))
    (t/is (m/nan? (sut/whittaker-W 1.0 0.5 -1.0))))
  (t/testing "x -> 0+ limit depends on the sign of mu+0.5"
    (t/is (m/zero? (sut/whittaker-W 1.0 0.5 0.0)))
    (t/is (m/pos-inf? (sut/whittaker-W 1.0 -0.8 0.0))))
  (t/testing "mu = -0.5 exactly: fixed this session, was ##NaN; correct
              limit is tricomis-U's own value at that point (varies by
              kappa, unlike whittaker-M's uniform 1.0)"
    (t/is (m/delta-eq 0.0 (sut/whittaker-W 1.0 -0.5 0.0) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq 0.0 (sut/whittaker-W 2.0 -0.5 0.0) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq 0.08597174606442001 (sut/whittaker-W -3.5 -0.5 0.0))))
  (t/testing "vs mpmath, grid of kappa, mu (incl. negative), x in [0.1, 1000]"
    (t/is (check3 sut/whittaker-W {:a (get-in @whittaker-w-reference [:kappa])
                                   :b (get-in @whittaker-w-reference [:mu])
                                   :x (get-in @whittaker-w-reference [:x])
                                   :ref (get-in @whittaker-w-reference [:ref])} 1.0e-5))))

;; Reference values for `hypergeometric-2F1` (Gauss's hypergeometric
;; function) below were computed once from mpmath (`mpmath.hyp2f1`) and are
;; stored in `test/resources/special/hyp2f1_reference.edn`, split into
;; `:generic` (a, b non-positive-integer-free, x in [-5, 1)), `:terminating`
;; (a or b a non-positive integer, x incl. past 1), `:x1-boundary` (x=1.0
;; exactly, c-a-b>0 sub-case only -- the c-a-b<=0 diverging sub-case is
;; covered by dedicated edge-case assertions instead, see below) and
;; `:near-integer-diff` (b-a very close to an integer, both signs, x incl.
;; past 1 -- the region that contained this session's bugs, see below). See
;; `utils/fastmath/dev/generate_hyp2f1_reference.py` for the generator. For
;; x > 1 with generic (non-terminating) a, b, c, mpmath returns a genuinely
;; complex value (a real branch cut, not a bug); not included in the bulk
;; reference.
;;
;; TWO genuine bugs were found and fixed in `hypergeometric-2F1` this
;; session, both confirmed against mpmath:
;;  - A genuine INFINITE LOOP (confirmed hung for several minutes, force-
;;    killed, not just a wrong value) for real x routing to the internal
;;    `inf-2F1` helper (roughly `|x| >= 1.39`) whenever invoked with
;;    `b < a` (root cause: an internal loop bounded by `(== n m)`, `n` only
;;    ever increasing from 0, never terminates for negative `m =
;;    round(b-a)`). This can happen even though the top-level function
;;    already normalizes `a <= b` once, because `general-2F1`'s own
;;    `c-a-b < 0` transformation can produce an internal recursive call
;;    with `b < a` again, e.g. `hypergeometric-2F1(1.5, 2.5, 3.5, 5.0)`
;;    hung (internally reaching `inf-2F1(2.0, 1.0, 3.5, 5.0)`). Fixed by
;;    swapping a, b (exact, by 2F1's own a<->b symmetry) at the start of
;;    `inf-2F1` whenever `b < a`, restoring `m >= 0`.
;;    NOTE: a first attempt at a blanket "b-a near an integer -> ##NaN"
;;    guard (before finding the true root cause above) was tried and
;;    reverted after it caused 69 test regressions in `regularized-beta`
;;    (which legitimately routes through `inf-2F1` with near-integer
;;    differences that converge fine) -- "b-a near an integer" alone does
;;    NOT predict a hang; the real condition is specifically `b < a`.
;;  - `a` or `b` a non-positive integer (an exact finite polynomial, valid
;;    for any x) was only handled when BOTH were non-positive integers AND
;;    `|x| < 0.72`, giving `##NaN` for every other combination (e.g. only
;;    `a` non-positive, or `|x| >= 0.72`, or `x > 1`) even though the true
;;    value is perfectly finite and real (e.g.
;;    `hypergeometric-2F1(-3.0, 2.5, 3.5, 5.0)` was `##NaN`, true value
;;    `-24.8658...`). Fixed with a direct polynomial evaluation, added
;;    ahead of every other special-cased branch.

(def ^:private hyp2f1-reference
  (delay (edn/read-string (slurp (io/resource "special/hyp2f1_reference.edn")))))

(t/deftest hypergeometric-2F1
  (t/testing "edge cases"
    (t/is (m/nan? (sut/hypergeometric-2F1 ##NaN 2.0 3.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-2F1 2.0 ##NaN 3.0 0.5)))
    (t/is (m/nan? (sut/hypergeometric-2F1 2.0 3.0 ##NaN 0.5)))
    (t/is (m/nan? (sut/hypergeometric-2F1 2.0 3.0 4.0 ##NaN)))
    (t/is (m/one? (sut/hypergeometric-2F1 2.0 3.0 4.0 0.0)))
    (t/is (m/one? (sut/hypergeometric-2F1 0.0 3.0 4.0 0.5))))
  (t/testing "x = 1.0: Gauss's summation theorem"
    (t/is (m/delta-eq 2.5 (sut/hypergeometric-2F1 1.0 1.5 3.5 1.0)))
    (t/is (m/pos-inf? (sut/hypergeometric-2F1 1.0 1.5 2.5 1.0))) ;; c-a-b=0
    (t/is (m/pos-inf? (sut/hypergeometric-2F1 1.5 2.5 3.0 1.0)))) ;; c-a-b<0
  (t/testing "x > 1.0, generic (non-terminating) a, b: leaves the real line, ##NaN"
    (t/is (m/nan? (sut/hypergeometric-2F1 1.5 2.5 3.5 5.0)))
    (t/is (m/nan? (sut/hypergeometric-2F1 1.5 2.5 3.5 1.4))))
  (t/testing "pole: c a non-positive integer, not terminated by a or b first"
    (t/is (m/nan? (sut/hypergeometric-2F1 1.5 2.5 -3.0 0.5))))
  (t/testing "terminating (a or b a non-positive integer): exact finite
              polynomial, valid past x=1 too (fixed this session, was
              ##NaN outside |x|<0.72 or when only one of a,b qualified)"
    (t/is (m/delta-eq -24.86580086580087 (sut/hypergeometric-2F1 -3.0 2.5 3.5 5.0)))
    (t/is (m/delta-eq 19.095238095238095 (sut/hypergeometric-2F1 -3.0 -2.0 3.5 5.0)))
    (t/is (m/delta-eq 6.78125 (sut/hypergeometric-2F1 -2.0 2.5 -3.0 1.5)))
    (t/is (m/nan? (sut/hypergeometric-2F1 -5.0 2.5 -3.0 1.5)))) ;; pole: n=5 > |c|=3
  (t/testing "no infinite loop, and correct value (##NaN, genuinely complex
              per mpmath), where b < a used to hang internally (fixed this
              session)"
    (t/is (m/nan? (sut/hypergeometric-2F1 1.45 2.5 3.5 5.0)))
    (t/is (m/nan? (sut/hypergeometric-2F1 1.5 2.5 3.5 1.4))))
  (t/testing "vs mpmath, generic (a, b non-positive-integer-free), x in [-5, 1)"
    (t/is (check4 sut/hypergeometric-2F1 {:a (get-in @hyp2f1-reference [:generic :a])
                                          :b (get-in @hyp2f1-reference [:generic :b])
                                          :c (get-in @hyp2f1-reference [:generic :c])
                                          :x (get-in @hyp2f1-reference [:generic :x])
                                          :ref (get-in @hyp2f1-reference [:generic :ref])} 1.0e-8)))
  (t/testing "vs mpmath, terminating case, x incl. past 1"
    (t/is (check4 sut/hypergeometric-2F1 {:a (get-in @hyp2f1-reference [:terminating :a])
                                          :b (get-in @hyp2f1-reference [:terminating :b])
                                          :c (get-in @hyp2f1-reference [:terminating :c])
                                          :x (get-in @hyp2f1-reference [:terminating :x])
                                          :ref (get-in @hyp2f1-reference [:terminating :ref])} 1.0e-9)))
  (t/testing "vs mpmath, x=1.0 boundary, c-a-b>0 sub-case"
    (t/is (check4 sut/hypergeometric-2F1 {:a (get-in @hyp2f1-reference [:x1-boundary :a])
                                          :b (get-in @hyp2f1-reference [:x1-boundary :b])
                                          :c (get-in @hyp2f1-reference [:x1-boundary :c])
                                          :x (get-in @hyp2f1-reference [:x1-boundary :x])
                                          :ref (get-in @hyp2f1-reference [:x1-boundary :ref])} 1.0e-9)))
  (t/testing "vs mpmath, b-a near an integer (both signs), x incl. past 1
              (the region that contained this session's bugs)"
    (t/is (check4 sut/hypergeometric-2F1 {:a (get-in @hyp2f1-reference [:near-integer-diff :a])
                                          :b (get-in @hyp2f1-reference [:near-integer-diff :b])
                                          :c (get-in @hyp2f1-reference [:near-integer-diff :c])
                                          :x (get-in @hyp2f1-reference [:near-integer-diff :x])
                                          :ref (get-in @hyp2f1-reference [:near-integer-diff :ref])} 1.0e-9))))
