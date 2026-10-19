(ns fastmath.kernel.variogram-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel.rbf :as rbf]
            [fastmath.kernel.variogram :as vg]
            [fastmath.optimization :as optim]
            [fastmath.random :as r]))

;; Reference data: `utils/fastmath/dev/generate_variogram_reference.py` (mpmath, scikit-gstat, numpy; file
;; `variogram_reference.edn`) and `generate_variogram_gstat_reference.R` (gstat; file
;; `variogram_gstat_reference.edn`).

(def ^:private reference (delay (edn/read-string (slurp (io/resource "kernel/variogram_reference.edn")))))
(def ^:private gstat-reference (delay (edn/read-string (slurp (io/resource "kernel/variogram_gstat_reference.edn")))))

(defn- max-err
  "Maximum absolute difference of two sequences of numbers; NaN if any of the pair is NaN."
  ^double [as bs]
  (reduce (fn [^double acc [a b]]
            (let [e (Math/abs (- (double a) (double b)))]
              (if (or (Double/isNaN e) (> e acc)) e acc)))
          0.0 (map vector as bs)))

(defn- ex-info-of
  "Data of the `ex-info` thrown by `f` (nil when nothing or another exception is thrown)."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e (assoc (ex-data e) ::message (.getMessage e)))
       (catch Throwable _ nil)))

(defn- invalid-input?
  "True when `f` throws `ex-info` with `:type :invalid-input` and the given data."
  [f expected]
  (let [d (ex-info-of f)]
    (and (= :invalid-input (:type d))
         (= expected (select-keys d (keys expected))))))

(def ^:private base {:nugget 0.3 :psill 2.0 :range 5.0})

(def ^:private constructors
  {:power vg/power :exponential vg/exponential :expower vg/expower :gaussian vg/gaussian
   :rational vg/rational :spherical vg/spherical :cubic vg/cubic :pentaspherical vg/pentaspherical
   :circular vg/circular :linear vg/linear :tpower vg/tpower :hole vg/hole :cauchy vg/cauchy
   :bessel vg/bessel :matern vg/matern :hyperspherical vg/hyperspherical :superspherical vg/superspherical
   :tplstable vg/tplstable})

(def ^:private beta-models #{:power :expower :rational :tpower :cauchy :tplstable})

(def ^:private bounded-models #{:spherical :cubic :pentaspherical :circular :linear :tpower
                                :hyperspherical :superspherical})

(defn- valid-params
  [model]
  (cond-> base (beta-models model) (assoc :beta 1.0)))

(deftest model-reference-test
  (testing "gstat (Sph, Exp, Gau, Pen, Cir, Lin, Wav), identical parametrisation"
    ;; measured maximum error 4.4e-16 for spherical, exponential, gaussian, pentaspherical, linear.
    ;; circular and hole (Wav): 1.7e-13, growing linearly with the lag; the gstat C code uses a truncated
    ;; constant (pi), the error against 80-digit mpmath is 1e-15 for fastmath and 1.7e-13 for gstat
    (doseq [{:keys [model nugget psill range lag gamma]} (:models @gstat-reference)]
      (let [f ((constructors model) {:nugget nugget :psill psill :range range})
            tol (if (#{:circular :hole} model) 3e-13 1e-14)]
        (is (< (max-err (map f lag) gamma) tol) (str model " " nugget " " range)))))
  (testing "mpmath, textbook formulas"
    ;; measured maximum relative error over all the cases: 2.3e-15 (power at the lag 10 * range)
    (doseq [{:keys [model nugget psill range beta lag gamma]} (:models @reference)]
      (let [f ((constructors model) (cond-> {:nugget nugget :psill psill :range range} beta (assoc :beta beta)))]
        (is (every? true? (map #(m/delta-eq %1 %2 1e-14 1e-14) (map f lag) gamma)) (str model " " beta " " nugget))))))

(deftest model-boundaries-test
  (doseq [[model ctor] constructors
          :let [f (ctor (valid-params model))
                {:keys [nugget range]} base]]
    (testing (str model)
      (is (== 0.0 (f 0.0)) "lag 0 is exactly 0")
      (is (== 0.0 (f -0.0)) "negative zero lag")
      (is (m/delta-eq nugget (f 1e-300) 1e-12) "tiny lag: the nugget")
      (is (m/delta-eq nugget (f 1e-12) 1e-9) "lag 1e-12")
      (is (m/nan? (f ##NaN)) "NaN lag")
      (is (every? #(and (Double/isFinite (double %)) (>= (double %) (- nugget 1e-12)))
                  (map f [1e-300 1e-9 0.1 1.0 range 3.0 (* 1.0001 range) (* 3 range) 1e6])))))
  (testing "bounded models are exactly constant beyond the range and reach the sill at the range"
    (doseq [model bounded-models
            :let [f ((constructors model) (valid-params model))
                  sill (+ (:nugget base) (:psill base))]]
      (is (m/delta-eq sill (f (:range base)) 1e-13) (str model " at the range"))
      (is (== sill (f (* 1.0000001 (:range base)))) (str model " just beyond the range"))
      (is (== sill (f (* 100 (:range base)))) (str model " far beyond the range"))
      (is (== sill (f ##Inf)) (str model " infinite lag"))))
  (testing "infinite lag of the sill models: the sill; of power: infinity"
    (doseq [model [:exponential :gaussian :expower :rational :cauchy :hole :bessel :matern :tplstable]
            :let [f ((constructors model) (valid-params model))]]
      (is (m/delta-eq 2.3 (f ##Inf) 1e-12) (str model)))
    (is (== ##Inf ((vg/power (valid-params :power)) ##Inf)))))

(deftest model-identities-test
  (let [p base
        xs [1e-6 0.01 0.1 0.5 1.0 2.0 5.0 9.0]]
    (testing "linear is the superspherical of order 0 and a power with beta 1 within the range"
      (is (< (max-err (map (vg/superspherical (assoc p :order 0.0)) xs) (map (vg/linear p) xs)) 1e-14))
      (is (< (max-err (map (vg/power (assoc p :beta 1.0)) [0.1 1.0 4.0]) (map (vg/linear p) [0.1 1.0 4.0])) 1e-14)))
    (testing "circular, spherical as hyper- and superspherical"
      (is (< (max-err (map (vg/hyperspherical (assoc p :order 2.0)) xs) (map (vg/circular p) xs)) 1e-14))
      (is (< (max-err (map (vg/hyperspherical (assoc p :order 3.0)) xs) (map (vg/spherical p) xs)) 1e-14))
      (is (< (max-err (map (vg/superspherical (assoc p :order 0.5)) xs) (map (vg/circular p) xs)) 1e-14))
      (is (< (max-err (map (vg/superspherical p) xs) (map (vg/spherical p) xs)) 1e-14) "default order 1.0"))
    (testing "exponential is expower beta 1, gaussian is expower beta 2 and matern order 0.5"
      (is (< (max-err (map (vg/expower (assoc p :beta 1.0)) xs) (map (vg/exponential p) xs)) 1e-15))
      (is (< (max-err (map (vg/expower (assoc p :beta 2.0)) xs) (map (vg/gaussian p) xs)) 1e-15))
      ;; measured 1.6e-15
      (is (< (max-err (map (vg/matern (assoc p :order 0.5)) xs) (map (vg/exponential p) xs)) 1e-14)))
    (testing "tpower beta 1 is linear"
      (is (< (max-err (map (vg/tpower (assoc p :beta 1.0)) xs) (map (vg/linear p) xs)) 1e-14)))))

(deftest exponential-doc-test
  (testing "range is the scale: 95% of the structured part at about 3 * range"
    (let [f (vg/exponential {:nugget 0.0 :psill 1.0 :range 4.0})]
      (is (m/delta-eq (- 1.0 (Math/exp -3.0)) (f 12.0) 1e-15))
      (is (m/delta-eq 0.95 (f (* 4.0 2.995732273553991)) 1e-14)))))

(deftest hole-test
  (testing "normalised sinc: the first time at the sill is x = range, then oscillation"
    (let [f (vg/hole {:nugget 0.1 :psill 2.0 :range 3.0})]
      (is (m/delta-eq 2.1 (f 3.0) 1e-15))
      (is (m/delta-eq 0.1 (f 1e-9) 1e-12))
      ;; sin(pi x / r) / (pi x / r) at x = 4.5, r = 3
      (is (m/delta-eq (+ 0.1 (* 2.0 (- 1.0 (/ (Math/sin (* Math/PI 1.5)) (* Math/PI 1.5))))) (f 4.5) 1e-14))
      (is (m/delta-eq 2.1 (f ##Inf) 1e-15) "tends to the sill"))))

(deftest rbf->variogram-test
  (testing "kernel is evaluated at x / range: the RBF scale multiplies the range"
    (let [vg-scale (vg/rbf->variogram (rbf/gaussian {:scale 2.0}))
          f (vg-scale {:nugget 0.2 :psill 1.5 :range 3.0})
          g (vg/gaussian {:nugget 0.2 :psill 1.5 :range 6.0})]
      (is (< (max-err (map f [0.0 0.1 1.0 3.0 6.0 9.0 20.0]) (map g [0.0 0.1 1.0 3.0 6.0 9.0 20.0])) 1e-15))))
  (testing "other kernels"
    (let [f ((vg/rbf->variogram (rbf/truncated-power {:k 1.0})) {:nugget 0.0 :psill 1.0 :range 4.0})]
      (is (m/delta-eq 0.5 (f 2.0) 1e-15))
      (is (m/delta-eq 1.0 (f 5.0) 1e-15))
      (is (== 0.0 (f 0.0)))))
  (testing "parameters are checked"
    (is (invalid-input? #((vg/rbf->variogram (rbf/gaussian)) {:nugget 0.0 :psill 1.0 :range 0.0})
                        {:model :rbf :property :range}))))

(deftest parameter-validation-test
  (doseq [[model ctor] constructors
          :let [ok (valid-params model)]]
    (testing (str model)
      (is (fn? (ctor ok)) "valid values are accepted")
      (is (fn? (ctor (assoc ok :nugget 0.0 :psill 0.0))) "nugget 0 and psill 0 are valid")
      (doseq [[k vs] {:range [0.0 -1.0 ##NaN ##Inf ##-Inf]
                      :psill [-1.0 -1e-300 ##NaN ##Inf]
                      :nugget [-1.0 -1e-300 ##NaN ##-Inf]}
              v vs]
        (is (invalid-input? #(ctor (assoc ok k v)) {:model model :property k}) (str k " " v)))
      (doseq [k (cond-> [:nugget :psill :range] (beta-models model) (conj :beta))]
        (is (invalid-input? #(ctor (dissoc ok k)) {:model model :property k}) (str "missing " k))
        (is (invalid-input? #(ctor (assoc ok k nil)) {:model model :property k}) (str "nil " k))
        (is (invalid-input? #(ctor (assoc ok k "1")) {:model model :property k}) (str "string " k)))
      (is (invalid-input? #(ctor nil) {:model model}) "nil parameter map")
      (is (invalid-input? #(ctor {}) {:model model}) "empty parameter map")))
  (testing "beta intervals"
    (doseq [model [:power :expower :tplstable]]
      (is (fn? ((constructors model) (assoc base :beta 2.0))) (str model " beta 2 is valid"))
      (is (fn? ((constructors model) (assoc base :beta 1e-10))) (str model " tiny beta is valid"))
      (doseq [beta [0.0 -0.5 2.0000001 3.0 ##NaN ##Inf]]
        (is (invalid-input? #((constructors model) (assoc base :beta beta)) {:model model :property :beta}) (str model " " beta))))
    (doseq [model [:rational :tpower :cauchy]]
      (is (fn? ((constructors model) (assoc base :beta 100.0))) (str model " large beta is valid"))
      (doseq [beta [0.0 -1.0 -1e-300 ##NaN ##Inf]]
        (is (invalid-input? #((constructors model) (assoc base :beta beta)) {:model model :property :beta}) (str model " " beta)))))
  (testing "orders and the Hurst exponent"
    (is (invalid-input? #(vg/matern (assoc base :order 0.0)) {:model :matern :property :order}))
    (is (invalid-input? #(vg/matern (assoc base :order -1.0)) {:model :matern :property :order}))
    (is (invalid-input? #(vg/bessel (assoc base :order -0.5)) {:model :bessel :property :order}))
    (is (fn? (vg/bessel (assoc base :order 0.0))))
    (is (invalid-input? #(vg/hyperspherical (assoc base :order 0.0)) {:model :hyperspherical :property :order}))
    (is (invalid-input? #(vg/superspherical (assoc base :order -1.0)) {:model :superspherical :property :order}))
    (is (invalid-input? #(vg/tplstable (assoc base :beta 1.0 :H 0.0)) {:model :tplstable :property :H}))
    (is (invalid-input? #(vg/tplstable (assoc base :beta 1.0 :H 1.0)) {:model :tplstable :property :H}))
    (is (fn? (vg/tplstable (assoc base :beta 1.0 :H 0.999)))))
  (testing "the message names the property and the value; the check is at construction, not per evaluation"
    (let [d (ex-info-of #(vg/exponential (assoc base :range -2.0)))]
      (is (= -2.0 (:value d)))
      (is (re-find #":range" (::message d)))
      (is (re-find #"-2.0" (::message d))))
    (is (fn? (vg/exponential base)))
    (is (m/delta-eq 0.3 ((vg/exponential base) 1e-300) 1e-12)))
  (testing "long values are accepted"
    (is (m/delta-eq ((vg/exponential {:nugget 0.0 :psill 1.0 :range 10.0}) 5.0)
                    ((vg/exponential {:nugget 0 :psill 1 :range 10}) 5) 1e-15))))

(deftest ->bessel-test
  (testing "default order 1.0, equal to bessel"
    (let [p base
          xs [1e-9 0.1 1.0 2.0 5.0 9.0 20.0]]
      (is (< (max-err (map ((vg/->bessel) p) xs) (map (vg/bessel p) xs)) 1e-15))
      (is (< (max-err (map ((vg/->bessel) p) xs) (map ((vg/->bessel 1.0) p) xs)) 1e-15))
      (is (> (max-err (map ((vg/->bessel) p) xs) (map ((vg/->bessel 0.0) p) xs)) 1e-3) "order 0 is not the default")))
  (testing "bessel reference (mpmath), orders 0 to 300, lags 1e-300 to 1000"
    ;; measured maximum error 8.1e-14 (order 300), 1.9e-14 at order 200 and 1.1e-16 below order 100
    (let [{:keys [order lag gamma]} (:bessel @reference)
          values (map (fn [o x] (((vg/->bessel o) {:nugget 0.0 :psill 1.0 :range 1.0}) x)) order lag)]
      (is (= (count lag) (count values) 256))
      (is (< (max-err values gamma) 1e-12))))
  (testing "no NaN and no infinity for orders up to 200 and lags from 1e-300 to 3 * range"
    (doseq [order [0.0 0.5 1.0 2.5 5.0 10.0 20.0 50.0 100.0 150.0 170.0 171.0 200.0]
            range [1e-3 1.0 250.0]
            :let [f ((vg/->bessel order) {:nugget 0.1 :psill 2.0 :range range})]]
      (doseq [rel [1e-300 1e-200 1e-100 1e-50 1e-20 1e-10 1e-5 1e-3 0.1 0.5 1.0 2.0 3.0]]
        (let [v (f (* rel range))]
          (is (and (Double/isFinite v) (>= v 0.1)) (str order " " range " " rel " " v))))))
  (testing "tiny lag tends to the nugget"
    (is (== 0.1 ((vg/bessel {:nugget 0.1 :psill 2.0 :range 1.0 :order 200.0}) 1e-300))))
  (testing "order 0.5: closed form 1 - sin(d) / d"
    (let [f (vg/bessel {:nugget 0.0 :psill 1.0 :range 2.0 :order 0.5})]
      (doseq [x [0.1 1.0 3.0 7.0]]
        (is (m/delta-eq (- 1.0 (/ (Math/sin (/ x 2.0)) (/ x 2.0))) (f x) 1e-14))))))

(deftest ->matern-test
  (testing "matern reference (mpmath), orders 0.1 to 200, lags 1e-300 to 100"
    ;; measured maximum error 4.4e-15 (orders below 20 and 20..22: old branch, orders from 20 up the log-domain expansion)
    (let [{:keys [order lag gamma]} (:matern @reference)
          values (map (fn [o x] (((vg/->matern o) {:nugget 0.0 :psill 1.0 :range 1.0}) x)) order lag)]
      (is (= 713 (count values)))
      (is (< (max-err values gamma) 1e-12))))
  (testing "no NaN and no infinity for orders up to 200 and lags from 1e-300 to 3 * range"
    (doseq [order [0.1 0.5 1.0 1.5 2.0 2.5 5.0 10.0 15.0 19.9 20.0 25.0 50.0 100.0 150.0 170.0 171.0 200.0]
            range [1e-3 1.0 250.0]
            :let [f ((vg/->matern order) {:nugget 0.1 :psill 2.0 :range range})]]
      (doseq [rel [1e-300 1e-200 1e-100 1e-50 1e-20 1e-10 1e-5 1e-3 0.1 0.5 1.0 2.0 3.0]]
        (let [v (f (* rel range))]
          (is (and (Double/isFinite v) (<= 0.1 v 2.1)) (str order " " range " " rel " " v))))))
  (testing "tiny lag tends to the nugget, also for large order"
    (doseq [order [1.5 5.0 10.0 20.0 100.0 200.0]]
      (is (== 0.1 ((vg/matern {:nugget 0.1 :psill 2.0 :range 1.0 :order order}) 1e-300)) (str order))))
  (testing "the default order is 1.0"
    (is (== ((vg/matern base) 1.5) ((vg/matern (assoc base :order 1.0)) 1.5)))
    (is (== (((vg/->matern) base) 1.5) (((vg/->matern 1.0) base) 1.5)))
    (is (fn? ((vg/->matern) base))))
  (testing "order 0.5 is the exponential model with the same range"
    (is (< (max-err (map (vg/matern (assoc base :order 0.5)) [0.1 1.0 5.0 20.0]) (map (vg/exponential base) [0.1 1.0 5.0 20.0])) 1e-14)))
  (testing "gaussian limit needs the range * sqrt(2)"
    (let [xs (range 0.05 3.0 0.05)
          p {:nugget 0.0 :psill 1.0 :range 1.0}
          f (vg/matern (assoc p :order 200.0))]
      ;; measured order 200: 2.0e-3 against the gaussian with the range sqrt(2), 0.249 against the gaussian with the range 1
      (is (< (max-err (map f xs) (map (vg/gaussian (assoc p :range (Math/sqrt 2.0))) xs)) 5e-3))
      (is (> (max-err (map f xs) (map (vg/gaussian p) xs)) 0.2)))))

(deftest hyperspherical-test
  (testing "hyperspherical reference (mpmath) inside the documented limits"
    ;; measured maximum error 3.5e-13 (odd dimension 31), 2.5e-14 for the other dimensions up to 60
    (let [{:keys [order lag gamma]} (:hyperspherical @reference)
          values (map (fn [o h] (((vg/->hyperspherical o) {:nugget 0.0 :psill 1.0 :range 1.0}) h)) order lag)]
      (is (= 304 (count values)))
      (is (< (max-err values gamma) 1e-12))))
  (testing "superspherical reference (mpmath) inside the documented limits"
    ;; measured maximum error 3.5e-13 (order 15)
    (let [{:keys [order lag gamma]} (:superspherical @reference)
          values (map (fn [o h] (((vg/->superspherical o) {:nugget 0.0 :psill 1.0 :range 1.0}) h)) order lag)]
      (is (= 176 (count values)))
      (is (< (max-err values gamma) 1e-12))))
  (testing "limits are documented"
    (is (re-find #"dimension 60" (:doc (meta #'vg/->hyperspherical))))
    (is (re-find #"d = -15" (:doc (meta #'vg/->hyperspherical))))
    (is (re-find #"Accuracy limits" (:doc (meta #'vg/->superspherical))))))

(deftest tplstable-test
  (testing "tplstable reference (mpmath expint), integer orders 1 + 2H/beta and H 0.8, beta 1.9"
    ;; measured maximum error 1.8e-15
    (doseq [{:keys [H beta nugget psill range lag gamma]} (:tplstable @reference)
            :let [f (vg/tplstable {:H H :beta beta :nugget nugget :psill psill :range range})]]
      (is (< (max-err (map f lag) gamma) 1e-12) (str H " " beta))))
  (testing "tplstable reference, fractional orders (needs the fixed fastmath.special/En: H 0.5, beta 1.2 was off by 5e-4)"
    ;; measured maximum error 6.3e-16 with the fixed En
    (doseq [{:keys [H beta nugget psill range lag gamma]} (:tplstable-fractional @reference)
            :let [f (vg/tplstable {:H H :beta beta :nugget nugget :psill psill :range range})]]
      (is (< (max-err (map f lag) gamma) 1e-12) (str H " " beta))))
  (testing "H defaults to 0.5, ->tplstable gives the same"
    (let [p (assoc base :beta 1.0)]
      (is (== ((vg/tplstable p) 2.0) ((vg/tplstable (assoc p :H 0.5)) 2.0)))
      (is (== ((vg/tplstable p) 2.0) (((vg/->tplstable) p) 2.0))))))

(deftest variogram->covariance-test
  (let [p {:nugget 0.1 :psill 2.0 :range 3.0}
        sill (+ 0.1 2.0)]
    (testing "exponential: C(h) = psill exp(-h / range) for h > 0 and the variance at zero"
      (let [c (vg/variogram->covariance (vg/exponential p) sill)]
        (is (== sill (c 0.0)))
        (doseq [h [1e-9 0.1 1.0 3.0 10.0 100.0]]
          (is (m/delta-eq (* 2.0 (Math/exp (- (/ h 3.0)))) (c h) 1e-15) (str h)))
        (is (m/delta-eq 0.0 (c ##Inf) 1e-15))
        (is (m/nan? (c ##NaN)))))
    (testing "bounded model: covariance is 0 beyond the range and psill - ... inside"
      (let [c (vg/variogram->covariance (vg/spherical p) sill)]
        (doseq [h [1.0 2.0 3.0]]
          (let [u (/ h 3.0)]
            (is (m/delta-eq (* 2.0 (- 1.0 (- (* 1.5 u) (* 0.5 u u u)))) (c h) 1e-15))))
        (is (== 0.0 (c 3.0000001)))
        (is (== 0.0 (c 100.0)))))
    (testing "gaussian, sill - gamma equals the correlation"
      (let [c (vg/variogram->covariance (vg/gaussian p) sill)]
        (is (m/delta-eq (* 2.0 (Math/exp -1.0)) (c 3.0) 1e-15))))
    (testing "invalid input"
      (is (invalid-input? #(vg/variogram->covariance (vg/gaussian p) -1.0) {:property :sill}))
      (is (invalid-input? #(vg/variogram->covariance (vg/gaussian p) ##NaN) {:property :sill}))
      (is (invalid-input? #(vg/variogram->covariance (vg/gaussian p) ##Inf) {:property :sill}))
      (is (invalid-input? #(vg/variogram->covariance (vg/gaussian p) nil) {:property :sill}))
      (is (invalid-input? #(vg/variogram->covariance nil 1.0) {:property :semivariogram}))
      (is (fn? (vg/variogram->covariance (vg/gaussian p) 0.0))))
    (testing "the sill of a fit"
      (let [res (vg/fit-params [{:h 1.0 :gamma 0.5 :n 10} {:h 2.0 :gamma 1.0 :n 10} {:h 3.0 :gamma 1.4 :n 10}]
                               :exponential {:rng (r/rng :isaac 1)})
            {:keys [nugget psill]} (:params res)
            c (vg/variogram->covariance (:semivariogram res) (+ nugget psill))]
        (is (m/delta-eq (+ nugget psill) (c 0.0) 1e-15))
        (is (m/delta-eq (- (+ nugget psill) ((:semivariogram res) 2.0)) (c 2.0) 1e-15))))))

;; ------------------------------------------------------------------ estimators

(defn- gamma-of
  "Half of the estimate of the variogram."
  ^double [estimator diffs]
  (* 0.5 (double (estimator (count diffs) diffs))))

(deftest estimators-test
  (let [cases (:estimators @reference)]
    (testing "matheron vs scikit-gstat" ;; measured 2.2e-16 relative
      (doseq [{:keys [diffs matheron]} cases]
        (is (m/delta-eq matheron (gamma-of vg/matheron-estimator diffs) 1e-14 1e-14))))
    (testing "cressie vs scikit-gstat (Hawkins-Cressie constant with the last term)" ;; measured 1.3e-15 relative
      (doseq [{:keys [diffs cressie n]} cases]
        (is (m/delta-eq cressie (gamma-of vg/cressie-estimator diffs) 1e-14 1e-14) (str n))))
    (testing "dowd vs scikit-gstat: 2.198109 * median(abs d)^2" ;; measured 2.8e-16 relative
      (doseq [{:keys [diffs dowd n]} cases]
        (is (m/delta-eq dowd (gamma-of vg/dowd-estimator diffs) 1e-14 1e-14) (str n))))
    (testing "genton vs the same estimator with the Weibull quantile of numpy"
      (doseq [{:keys [diffs genton-weibull n]} cases
              :when genton-weibull]
        (is (m/delta-eq genton-weibull (gamma-of vg/highly-robust-estimator diffs) 1e-12 1e-12) (str n))))
    (testing "m-robust vs an independent numpy implementation"
      (doseq [{:keys [diffs m-robust n]} cases
              :when m-robust]
        (is (m/delta-eq m-robust (gamma-of vg/robust-m-estimator diffs) 1e-12 1e-12) (str n)))))
  (testing "dowd is the squared median of the absolute differences, not the median of the squares"
    ;; even number of pairs: |d| = 1 1 3 3, median 2, (median |d|)^2 = 4 while median(d^2) = 5
    (is (m/delta-eq (* 2.198109338317728 4.0) (vg/dowd-estimator 4 [1.0 -1.0 3.0 -3.0]) 1e-13))
    (is (m/delta-eq (* 2.198109338317728 9.0) (vg/dowd-estimator 3 [1.0 -3.0 3.0]) 1e-13))
    (is (m/delta-eq (* 2.198109338317728 0.0) (vg/dowd-estimator 1 [0.0]) 1e-13))
    (is (m/nan? (vg/dowd-estimator 0 []))))
  (testing "the constant is the reciprocal of the chi-squared(1) median"
    (is (m/delta-eq 2.198109338317728 (/ 1.0 (r/icdf (r/distribution :chi-squared) 0.5)) 1e-12)))
  (testing "cressie constant, N = 1: (|d|^0.5)^4 / (0.457 + 0.494 + 0.045)"
    (is (m/delta-eq (/ 4.0 0.996) (vg/cressie-estimator 1 [2.0]) 1e-14)))
  (testing "matheron, empty and single"
    (is (m/nan? (vg/matheron-estimator 0 [])))
    (is (m/delta-eq 9.0 (vg/matheron-estimator 1 [-3.0]) 1e-15)))
  (testing "quantile estimator"
    (let [e (vg/->quantile-estimator 0.5)]
      (is (m/delta-eq (/ 4.0 0.45493642311957283) (e 3 [1.0 -2.0 3.0]) 1e-9)))
    (is (invalid-input? #(vg/->quantile-estimator 0.0) {:property :quantile}))
    (is (invalid-input? #(vg/->quantile-estimator 1.0) {:property :quantile}))
    (is (invalid-input? #(vg/->quantile-estimator ##NaN) {:property :quantile}))
    (is (invalid-input? #(vg/->quantile-estimator -0.1) {:property :quantile})))
  (testing "genton and m-robust are not defined for a single pair"
    (is (m/nan? (vg/highly-robust-estimator 1 [2.0])))
    (is (m/nan? (vg/robust-m-estimator 1 [2.0])))))

;; ------------------------------------------------------------------ empirical

(deftest empirical-gstat-test
  (testing "classical estimator equals gstat variogram with identical bin edges, 6 data sets"
    ;; measured: counts and mean lags identical, maximum gamma difference 4.4e-15 (relative 5e-15)
    (doseq [{:keys [name xs ys cutoff bins n h gamma]} (:empirical @gstat-reference)
            :let [res (vg/empirical xs ys {:cutoff cutoff :bins bins})]]
      (testing (str name)
        (is (= n (map :n res)))
        (is (< (max-err h (map :h res)) 1e-13))
        (is (< (max-err gamma (map :gamma res)) 1e-13))
        (is (every? integer? (map :n res)))
        (is (= (* (count xs) (dec (count xs)) 0.5) (double (reduce + 0 (map :n (vg/empirical xs ys {:cutoff 1e9 :bins bins})))))
            "all pairs are binned for a large cutoff")))))

(deftest empirical-test
  (let [xs [[0 0] [1 0] [2 0] [3 0]]
        ys [1.0 2.0 5.0 4.0]]
    (testing "n is the number of pairs in the bin, h the mean lag, gamma the classical estimate"
      (let [res (vg/empirical xs ys {:cutoff 5.0 :bins 3})]
        ;; d = 1: 3 pairs, diffs -1 -3 1; d = 2: 2 pairs, diffs -4 -2; d = 3: 1 pair, diff -3
        (is (= [3 2 1] (map :n res)))
        (is (= [1.0 2.0 3.0] (map :h res)))
        (is (m/delta-eq (/ 11.0 6.0) (:gamma (first res)) 1e-14))
        (is (m/delta-eq (/ 20.0 4.0) (:gamma (second res)) 1e-14))
        (is (m/delta-eq (/ 9.0 2.0) (:gamma (nth res 2)) 1e-14))
        (is (vector? res))))
    (testing "scalar locations"
      (is (= (map :n (vg/empirical [0.0 1.0 2.0 3.0] ys {:cutoff 5.0 :bins 3}))
             [3 2 1])))
    (testing "one bin takes all the pairs, bins are on (0, largest kept distance]"
      (let [res (vg/empirical xs ys {:cutoff 5.0 :bins 1})]
        (is (= [6] (map :n res)))
        (is (m/delta-eq (/ 10.0 6.0) (:h (first res)) 1e-14))))
    (testing "cutoff keeps only the near pairs; the lag exactly at the cutoff is kept"
      (is (= [3] (map :n (vg/empirical xs ys {:cutoff 1.0}))))
      (is (= [3 2] (map :n (vg/empirical xs ys {:cutoff 2.0 :bins 2}))))
      (is (= [3] (map :n (vg/empirical xs ys {:cutoff 1.5 :bins 4})))))
    (testing "default cutoff is the bounding box diagonal / 3"
      (is (= [3] (map :n (vg/empirical xs ys)))))
    (testing ":diagonal-den and :distance"
      (is (= [3 2] (map :n (vg/empirical xs ys {:diagonal-den 1.5 :bins 3}))))
      (is (= [3 2 1] (map :n (vg/empirical xs ys {:cutoff 5.0 :bins 3 :distance (fn [a b] (Math/abs (- (double (first a)) (double (first b)))))}))))
      (is (= [6] (map :n (vg/empirical xs ys {:cutoff 5.0 :bins 1 :distance (fn [_ _] 1.0)})))))
    (testing "nil options are the same as missing options"
      (is (= (vg/empirical xs ys {:cutoff 5.0 :bins 3})
             (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator nil :quantile nil :distance nil :remove-outliers? nil
                                  :diagonal-den nil})))
      (is (= (vg/empirical xs ys {:cutoff 2.0}) (vg/empirical xs ys {:cutoff 2.0 :bins nil}))))
    (testing "custom estimator and the quantile estimator"
      (is (= [2.0 2.0 2.0] (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator (fn [n _] (double (* 4 (/ n n))))}))))
      (is (= (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator :quantile :quantile 0.5}))
             (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator (vg/->quantile-estimator 0.5)})))))
    (testing "estimator keywords"
      (doseq [[k f] {:classical vg/matheron-estimator :matheron vg/matheron-estimator :cressie vg/cressie-estimator
                     :dowd vg/dowd-estimator}]
        (let [res (vg/empirical xs ys {:cutoff 5.0 :bins 1 :estimator k})
              diffs (for [i (range 4) j (range (inc i) 4)] (- (double (ys i)) (double (ys j))))]
          (is (m/delta-eq (* 0.5 (double (f 6 diffs))) (:gamma (first res)) 1e-13) (str k))))
      (is (= (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 1 :estimator :genton}))
             (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 1 :estimator :highly-robust}))))))
  (testing "zero-distance pairs (duplicate locations) are dropped"
    (let [res (vg/empirical [[0 0] [0 0] [1 0] [2 0]] [1.0 2.0 3.0 4.0] {:cutoff 5.0})]
      (is (= [3 2] (map :n res)))
      (is (= [1.0 2.0] (map :h res)))
      (is (m/delta-eq 1.0 (:gamma (first res)) 1e-15))
      (is (m/delta-eq 3.25 (:gamma (second res)) 1e-15))))
  (testing "robust estimators drop the bins with a non-finite estimate"
    ;; bins of 3, 2 and 1 pair: the single pair has no estimate
    (let [xs [[0 0] [1 0] [2 0] [3 0]]
          ys [1.0 2.0 5.0 4.0]]
      (is (= [3 2] (map :n (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator :genton}))))
      (is (every? #(Double/isFinite (double %)) (map :gamma (vg/empirical xs ys {:cutoff 5.0 :bins 3 :estimator :genton}))))
      (is (= [3 2] (map :n (vg/empirical xs [0.0 1.0 3.0 7.0] {:cutoff 5.0 :bins 3 :estimator :m-robust})))))
    (testing "a bin which is not defined for the data of the bin (ties) is dropped, the others stay"
      (is (= [] (try (vec (vg/empirical [[0 0] [1 0] [2 0] [3 0]] [1.0 2.0 3.0 4.0]
                                        {:cutoff 5.0 :bins 3 :estimator :m-robust}))
                     (catch clojure.lang.ExceptionInfo _ [])))))
    (testing "NaN in the observations is not dropped"
      (is (some #(m/nan? (:gamma %)) (vg/empirical [[0 0] [1 0] [2 0] [3 0]] [1.0 ##NaN 5.0 4.0] {:cutoff 5.0 :bins 3})))
      (is (some #(m/nan? (:gamma %)) (vg/empirical [[0 0] [1 0] [2 0] [3 0]] [1.0 ##NaN 5.0 4.0]
                                                   {:cutoff 5.0 :bins 3 :estimator :genton})))))
  (testing "outliers are removed from the points"
    (let [xs (mapv (fn [i] [(double i) 0.0]) (range 10))
          ys [1.0 1.1 0.9 1.0 100.0 1.2 0.8 1.0 1.1 0.9]
          all (vg/empirical xs ys {:cutoff 20.0 :bins 1})
          fence (vg/empirical xs ys {:cutoff 20.0 :bins 1 :remove-outliers? true})
          mad (vg/empirical xs ys {:cutoff 20.0 :bins 1 :remove-outliers? :mad})]
      (is (= [45] (map :n all)))
      (is (= [36] (map :n fence)))
      (is (= [36] (map :n mad)))
      (is (< (:gamma (first fence)) 0.1))
      (is (> (:gamma (first all)) 100.0))
      (is (= [[0 [0.0] 1.0]] (vg/maybe-remove-outliers [[0 [0.0] 1.0]] [1.0] false)))
      (is (= [[0 [0.0] 1.0]] (vg/maybe-remove-outliers [[0 [0.0] 1.0]] [1.0] nil)))))
  (testing "errors"
    (let [xs [[0 0] [1 0] [2 0]]
          ys [1.0 2.0 3.0]]
      (is (invalid-input? #(vg/empirical [] []) {:property :xs}) "empty")
      (is (invalid-input? #(vg/empirical [[0 0]] [1.0]) {:property :xs}) "one point")
      (is (invalid-input? #(vg/empirical [[0 0] [0 0]] [1.0 2.0]) {:property :cutoff}) "one location")
      (is (invalid-input? #(vg/empirical [[0 0] [0 0] [0 0]] [1.0 2.0 3.0] {:cutoff 5.0}) {:property :cutoff}) "all duplicates")
      (is (invalid-input? #(vg/empirical xs ys {:cutoff 0.5}) {:property :cutoff}) "cutoff too small")
      (is (invalid-input? #(vg/empirical xs ys {:cutoff 0.0}) {:property :cutoff}) "cutoff 0")
      (is (invalid-input? #(vg/empirical xs ys {:cutoff -1.0}) {:property :cutoff}) "negative cutoff")
      (is (invalid-input? #(vg/empirical xs ys {:cutoff ##NaN}) {:property :cutoff}) "NaN cutoff")
      (is (invalid-input? #(vg/empirical xs [1.0 2.0]) {:property :ys}) "unequal lengths")
      (is (invalid-input? #(vg/empirical [[0 0] [1 0]] [1.0 2.0 3.0]) {:property :ys}) "unequal lengths, more ys")
      (is (invalid-input? #(vg/empirical xs ys {:bins 0}) {:property :bins}) ":bins 0")
      (is (invalid-input? #(vg/empirical xs ys {:bins -3}) {:property :bins}) ":bins negative")
      (is (invalid-input? #(vg/empirical xs ys {:estimator :nope}) {:property :estimator}) "unknown estimator")
      (is (invalid-input? #(vg/empirical xs ys {:estimator "classical"}) {:property :estimator}) "string estimator")
      (is (invalid-input? #(vg/empirical xs ys {:estimator :quantile :quantile 1.5}) {:property :quantile}) "quantile level")
      (is (invalid-input? #(vg/empirical [[0 0] [1 0]] [1.0 2.0] {:cutoff 5.0 :estimator :genton}) {:property :estimator})
          "no bin with an estimate")
      (let [d (ex-info-of #(vg/empirical xs ys {:estimator :nope}))]
        (is (some #{:classical} (:expected d)) "the message lists the accepted values")
        (is (re-find #":nope" (::message d)))))))

;; ------------------------------------------------------------------ fit

(defn- synthetic
  "Empirical semivariogram from a model on lags 2, 5, ..., 38."
  [f]
  (vec (for [h (range 2.0 40.0 3.0)] {:h h :gamma (f h) :n 50})))

(deftest fit-tpower-test
  (let [truth {:nugget 0.1 :psill 2.0 :range 25.0 :beta 1.5}
        emp (synthetic (vg/tpower truth))
        t0 (System/nanoTime)
        res (vg/fit-params emp :tpower {:rng (r/rng :isaac 1)})]
    (testing "tpower is fitted by tpower (the loss of power was 21.27)"
      ;; measured: loss 2.8e-19, parameters within 2e-9
      (is (< (:loss res) 1e-12))
      (doseq [[k v] truth]
        (is (m/delta-eq v (k (:params res)) 1e-6) (str k)))
      (is (< (max-err (map (:semivariogram res) [1.0 10.0 24.0 25.0 30.0 100.0])
                      (map (vg/tpower truth) [1.0 10.0 24.0 25.0 30.0 100.0]))
             1e-6))
      (is (m/delta-eq 2.1 ((:semivariogram res) 1000.0) 1e-6) "constant beyond the range, unlike power")
      (is (m/delta-eq ((vg/tpower (:params res)) 12.0) (((:semivariogram-creator res) (:params res)) 12.0) 1e-15))
      (is (> (Math/abs (- ((vg/power (:params res)) 12.0) ((vg/tpower (:params res)) 12.0))) 1e-3)
          "the creator is not power"))
    (testing "fast"
      (is (< (/ (- (System/nanoTime) t0) 1e9) 20.0)))
    (testing "fit returns the function"
      (is (m/delta-eq ((vg/tpower truth) 10.0) ((vg/fit emp :tpower {:rng (r/rng :isaac 1)}) 10.0) 1e-6)))))

(deftest fit-bounds-test
  (testing "the default beta bounds are the valid intervals: power family models are fitted without errors"
    (let [truth {:nugget 0.1 :psill 1.0 :range 30.0}]
      (doseq [model [:power :expower :rational :cauchy :tplstable :tpower]
              :let [emp (synthetic (vg/exponential truth))
                    res (vg/fit-params emp model {:rng (r/rng :isaac 3) :N 50})
                    beta (:beta (:params res))]]
        (is (some? beta) (str model))
        (is (<= 0.0 beta 30.0) (str model " " beta))
        (when (#{:power :expower :tplstable} model) (is (<= beta 2.0) (str model " " beta)))
        (is (fn? (:semivariogram res))))))
  (testing "the optimum at the upper bound of beta (a numerical gradient steps over it, the model would reject beta > 2)"
    (let [emp (synthetic (vg/gaussian {:nugget 0.1 :psill 1.0 :range 20.0}))
          res (vg/fit-params emp :expower {:rng (r/rng :isaac 11) :N 60})]
      (is (m/delta-eq 2.0 (:beta (:params res)) 1e-6))
      (is (< (:loss res) 1e-10))))
  (testing "explicit bounds are used"
    (let [emp (synthetic (vg/gaussian {:nugget 0.0 :psill 1.0 :range 100.0}))
          narrow (vg/fit-params emp :gaussian {:rng (r/rng :isaac 2) :bounds {:range [1.0 50.0]}})
          wide (vg/fit-params emp :gaussian {:rng (r/rng :isaac 2) :bounds {:range [1.0 500.0]}})]
      (is (<= (:range (:params narrow)) 50.0))
      (is (< (:loss wide) (:loss narrow)))))
  (testing "defaults fix parameters; when all are fixed there is no optimisation"
    (let [emp (synthetic (vg/gaussian {:nugget 0.0 :psill 1.0 :range 20.0}))
          res (vg/fit-params emp :gaussian {:rng (r/rng :isaac 2) :defaults {:nugget 0.0}})]
      (is (== 0.0 (:nugget (:params res))))
      (is (m/delta-eq 20.0 (:range (:params res)) 1e-5)))
    (let [res (vg/fit-params [{:h 1.0 :gamma 0.5 :n 5}] :exponential {:defaults {:nugget 0.1 :psill 1.0 :range 2.0}})]
      (is (= {:nugget 0.1 :psill 1.0 :range 2.0} (:params res))))
    (let [res (vg/fit-params [{:h 1.0 :gamma 0.5 :n 5} {:h 2.0 :gamma 0.9 :n 5}] :power
                             {:defaults {:nugget 0.1 :psill 1.0 :range 2.0 :beta 1.5}})]
      (is (= 1.5 (:beta (:params res))))
      (is (m/delta-eq ((vg/power {:nugget 0.1 :psill 1.0 :range 2.0 :beta 1.5}) 1.0) ((:semivariogram res) 1.0) 1e-15))))
  (testing "unknown keys in :defaults and :bounds"
    (let [emp (synthetic (vg/gaussian {:nugget 0.0 :psill 1.0 :range 20.0}))]
      (is (invalid-input? #(vg/fit-params emp :gaussian {:defaults {:sill 1.0}}) {:property :defaults}))
      (is (invalid-input? #(vg/fit-params emp :gaussian {:bounds {:rng [1.0 2.0]}}) {:property :bounds})))))

(deftest fit-input-test
  (let [emp (synthetic (vg/gaussian {:nugget 0.1 :psill 1.0 :range 20.0}))]
    (testing "NaN in gamma throws at once (it hung the fit for minutes)"
      (let [t0 (System/nanoTime)]
        (is (invalid-input? #(vg/fit-params (assoc-in emp [3 :gamma] ##NaN) :gaussian) {:property :gamma}))
        (is (invalid-input? #(vg/fit (assoc-in emp [0 :gamma] ##Inf) :gaussian) {:property :gamma}))
        (is (invalid-input? #(vg/fit (assoc-in emp [0 :h] ##NaN) :gaussian) {:property :h}))
        (is (invalid-input? #(vg/fit (assoc-in emp [0 :n] nil) :gaussian) {:property :n}))
        (is (invalid-input? #(vg/fit (assoc-in emp [0 :n] ##Inf) :gaussian) {:property :n}))
        (is (< (/ (- (System/nanoTime) t0) 1e9) 1.0))))
    (testing "NaN in the observations reaches the fit as NaN gamma"
      (is (invalid-input? #(vg/fit (vg/empirical [[0 0] [1 0] [2 0] [3 0]] [1.0 ##NaN 5.0 4.0] {:cutoff 5.0 :bins 3}) :gaussian)
                          {:property :gamma})))
    (testing "empty and all-zero empirical semivariogram"
      (is (invalid-input? #(vg/fit [] :gaussian) {:property :empirical-semivariogram}))
      (is (invalid-input? #(vg/fit nil :gaussian) {:property :empirical-semivariogram}))
      (is (invalid-input? #(vg/fit [{:h 1.0 :gamma 0.0 :n 3} {:h 2.0 :gamma 0.0 :n 3}] :gaussian) {:property :gamma})))
    (testing "unknown model, estimation and weights list the accepted values"
      (let [d (ex-info-of #(vg/fit emp :nope))]
        (is (= :invalid-input (:type d)))
        (is (= :model (:property d)))
        (is (some #{:gaussian} (:expected d)))
        (is (re-find #":nope" (::message d))))
      (is (invalid-input? #(vg/fit emp nil) {:property :model}))
      (is (invalid-input? #(vg/fit emp "gaussian") {:property :model}))
      (let [d (ex-info-of #(vg/fit emp :gaussian {:estimation :nope}))]
        (is (= :estimation (:property d)))
        (is (some #{:huber} (:expected d))))
      (let [d (ex-info-of #(vg/fit emp :gaussian {:weights :nope}))]
        (is (= :weights (:property d)))
        (is (re-find #":ngg2" (:expected d))))
      (is (invalid-input? #(vg/fit emp :gaussian {:weights [1.0 2.0]}) {:property :weights}))
      (is (invalid-input? #(vg/fit emp :gaussian {:weights (repeat 14 ##NaN)}) {:property :weights}))
      (is (invalid-input? #(vg/fit emp :gaussian {:weights (concat (repeat 12 1.0) [nil])}) {:property :weights})))
    (testing "valid options: all estimations, all weights"
      (doseq [estimation [:sq :abs :soft :huber :cauchy :atan]]
        (is (fn? (vg/fit emp :gaussian {:estimation estimation :rng (r/rng :isaac 4) :N 30})) (str estimation)))
      (doseq [weights [nil false :n :nh :nh+ :nhh :nhh+ :nexp :nsqrt :ngg :ngg2 (repeat 13 2.0) (fn [^double h] (+ 1.0 h))]]
        (is (fn? (vg/fit emp :gaussian {:weights weights :rng (r/rng :isaac 4) :N 30})) (str weights))))
    (testing "nil options are missing options"
      (let [a (vg/fit-params emp :gaussian {:rng (r/rng :isaac 9)})
            b (vg/fit-params emp :gaussian {:rng (r/rng :isaac 9) :estimation nil :defaults nil
                                            :bounds nil :parameter nil :N nil})]
        (is (= (:params a) (:params b)))))
    (testing "except :weights: nil means the uniform weights"
      (is (every? #(== 1.0 %) (:weights (vg/fit-params emp :gaussian {:weights nil :rng (r/rng :isaac 9)})))))
    (testing "a function is a kernel; a model constructor is not accepted"
      (let [res (vg/fit-params emp (rbf/gaussian) {:rng (r/rng :isaac 5)})]
        (is (< (:loss res) 1e-6))
        (is (m/delta-eq 20.0 (:range (:params res)) 1e-4)))
      (is (some? (ex-info-of #(vg/fit emp vg/gaussian)))))))

(deftest fit-weights-test
  (let [truth {:nugget 0.1 :psill 1.0 :range 20.0}
        emp (mapv (fn [e] (update e :gamma #(+ % (* 0.03 (Math/sin (:h e)))))) (synthetic (vg/gaussian truth)))
        eps m/EPSILON]
    (testing ":ngg2 weights are n * gamma / model^3 with non-truncated n * gamma"
      ;; measured maximum relative difference 2.2e-16
      (let [res (vg/fit-params emp :gaussian {:weights :ngg2 :rng (r/rng :isaac 6)})
            sv (:semivariogram res)
            expected (map (fn [{:keys [h gamma n]}] (/ (* n gamma) (max eps (m/cb (sv h))))) emp)]
        (is (every? true? (map #(m/delta-eq %1 %2 0.0 1e-14) expected (:weights res))))
        (is (> (count (filter #(not (zero? (- (* (:n %) (:gamma %)) (Math/floor (* (:n %) (:gamma %)))))) emp)) 8)
            "n * gamma is not an integer, so a truncation would show")))
    (testing ":ngg weights are n / model^2"
      (let [res (vg/fit-params emp :gaussian {:weights :ngg :rng (r/rng :isaac 6)})
            sv (:semivariogram res)
            expected (map (fn [{:keys [h n]}] (/ n (max eps (m/sq (sv h))))) emp)]
        (is (every? true? (map #(m/delta-eq %1 %2 0.0 1e-14) expected (:weights res))))))
    (testing "other weights"
      (let [ns (map :n emp)]
        (is (= ns (:weights (vg/fit-params emp :gaussian {:weights :n :defaults {:nugget 0.1 :psill 1.0 :range 20.0}}))))
        (is (every? #(== 1.0 %) (:weights (vg/fit-params emp :gaussian {:weights nil :defaults {:nugget 0.1 :psill 1.0 :range 20.0}}))))
        (is (= (repeat 13 2.0) (:weights (vg/fit-params emp :gaussian {:weights (repeat 13 2.0)
                                                                      :defaults {:nugget 0.1 :psill 1.0 :range 20.0}}))))
        (is (every? true? (map #(m/delta-eq (/ %1 %2) %3 1e-14)
                               ns (map :h emp)
                               (:weights (vg/fit-params emp :gaussian {:weights :nh :defaults {:nugget 0.1 :psill 1.0 :range 20.0}})))))))
    (testing ":ngg2 on integer-free n: the old truncation changed the weights by up to 12%"
      (let [res (vg/fit-params emp :gaussian {:weights :ngg2 :defaults {:nugget 0.1 :psill 1.0 :range 20.0}})
            w (first (:weights res))
            sv (:semivariogram res)
            e (first emp)
            truncated (/ (Math/floor (* (:n e) (:gamma e))) (max eps (m/cb (sv (:h e)))))]
        (is (not (m/delta-eq w truncated 0.0 1e-4)))))))

(deftest fit-parameter-test
  (let [emp (synthetic (vg/gaussian {:nugget 0.1 :psill 1.0 :range 20.0}))]
    (testing ":bessel default order is 1.0, the same as bessel and ->bessel"
      (let [res (vg/fit-params emp :bessel {:rng (r/rng :isaac 7) :N 30})
            p (:params res)
            creator (:semivariogram-creator res)]
        (is (< (max-err (map (creator p) [1.0 5.0 10.0 30.0]) (map (vg/bessel p) [1.0 5.0 10.0 30.0])) 1e-15))
        (is (> (max-err (map (creator p) [5.0 10.0]) (map ((vg/->bessel 0.0) p) [5.0 10.0])) 1e-3))))
    (testing ":parameter is the order / H"
      (let [res (vg/fit-params emp :matern {:rng (r/rng :isaac 7) :N 30 :parameter 2.5})]
        (is (m/delta-eq ((vg/matern (assoc (:params res) :order 2.5)) 7.0) ((:semivariogram res) 7.0) 1e-15)))
      (let [res (vg/fit-params emp :bessel {:rng (r/rng :isaac 7) :N 30 :parameter 0.5})]
        (is (m/delta-eq ((vg/bessel (assoc (:params res) :order 0.5)) 7.0) ((:semivariogram res) 7.0) 1e-15)))
      (is (invalid-input? #(vg/fit-params emp :matern {:parameter -1.0}) {:model :matern :property :order})))))

(deftest fit-rng-test
  (let [emp (vec (map (fn [e] (update e :gamma #(+ % (* 0.05 (Math/cos (* 3.0 (:h e))))))) (synthetic (vg/spherical {:nugget 0.2 :psill 1.0 :range 25.0}))))]
    (testing "the same seeded :rng gives bitwise identical parameters, for several models"
      (doseq [model [:spherical :gaussian :matern :hole :linear :cauchy :power :tplstable]]
        (let [a (vg/fit-params emp model {:rng (r/rng :mersenne 123)})
              b (vg/fit-params emp model {:rng (r/rng :mersenne 123)})]
          (is (= (:params a) (:params b)) (str model))
          (is (== (:loss a) (:loss b)) (str model))
          (is (= (:weights a) (:weights b))))))
    (testing "different seeds are different runs of the same optimum"
      (let [a (vg/fit-params emp :spherical {:rng (r/rng :mersenne 1)})
            b (vg/fit-params emp :spherical {:rng (r/rng :mersenne 2)})]
        (is (m/delta-eq (:loss a) (:loss b) 1e-6))))
    (testing ":rng and :N are passed to scan-and-minimize"
      (let [seen (atom nil)
            orig optim/scan-and-minimize
            rng (r/rng :mersenne 5)]
        (with-redefs [optim/scan-and-minimize (fn [method f opts] (reset! seen opts) (orig method f opts))]
          (vg/fit emp :spherical {:rng rng :N 17}))
        (is (identical? rng (:rng @seen)))
        (is (== 17 (:N @seen)))
        (is (contains? @seen :bounds))
        (reset! seen nil)
        (with-redefs [optim/scan-and-minimize (fn [method f opts] (reset! seen opts) (orig method f opts))]
          (vg/fit emp :spherical))
        (is (not (contains? @seen :rng)) "without :rng the optimiser draws a new generator")
        (is (not (contains? @seen :N)))))))

(deftest fit-kriging-api-test
  (testing "the calls used by interpolation.kriging keep working"
    (let [xss (mapv (fn [i] [(double (mod (* i 7) 11)) (double (mod (* i 5) 13))]) (range 30))
          ys (mapv (fn [[x y]] (+ (Math/sin (/ x 3.0)) (* 0.1 y))) xss)
          emp (vg/empirical xss ys)
          f (vg/fit emp :gaussian)]
      (is (seq emp))
      (is (every? (fn [e] (and (integer? (:n e)) (Double/isFinite (:h e)) (Double/isFinite (:gamma e)))) emp))
      (is (== 0.0 (f 0.0)))
      (is (Double/isFinite (f 3.0))))))
