(ns fastmath.kernel.window-test
  "Tests of `fastmath.kernel.window`: N contract, `:normalize?` contract, parameter domains, oracle
  comparison (mpmath, scipy, Octave) of the windows changed by the kernel audit.
  Properties (integrals, means, Princen-Bradley) and identities are in `fastmath.kernel.window2-test`."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel.window :as win]))

;; All windows defined in `fastmath.kernel.window` (the registry of `fastmath.kernel` is tested in
;; `fastmath.kernel.registry-test`); a new window is added here by hand.
(def ^:private window-fns
  {:rectangular05 win/rectangular05
   :rectangular win/rectangular
   :triangular win/triangular
   :parzen win/parzen
   :b-spline win/b-spline
   :welch win/welch
   :connes win/connes
   :parzen-algebraic win/parzen-algebraic
   :singla-singh win/singla-singh
   :sinc win/sinc
   :fejer win/fejer
   :de-la-vallee-poussin win/de-la-vallee-poussin
   :lanczos win/lanczos
   :hamming win/hamming
   :hamming-exact win/hamming-exact
   :hann win/hann
   :raised-cosine win/raised-cosine
   :webster-hamming win/webster-hamming
   :power-of-cosine win/power-of-cosine
   :raised-power-of-cosine win/raised-power-of-cosine
   :parzen-cosine win/parzen-cosine
   :bohman win/bohman
   :trapezoid win/trapezoid
   :tukey win/tukey
   :bartlett-hann win/bartlett-hann
   :blackman-harris-family win/blackman-harris-family
   :blackman win/blackman
   :blackman-exact win/blackman-exact
   :blackman-harris win/blackman-harris
   :blackman-harris-61db win/blackman-harris-61db
   :blackman-harris-67db win/blackman-harris-67db
   :blackman-harris-74db win/blackman-harris-74db
   :blackman-harris-92db win/blackman-harris-92db
   :nutall-3-1st win/nutall-3-1st
   :nutall-3-3rd win/nutall-3-3rd
   :blackman-nutall win/blackman-nutall
   :nutall-1st win/nutall-1st
   :nutall-3rd win/nutall-3rd
   :nutall-5th win/nutall-5th
   :mottaghi-kashtiban-shayesteh win/mottaghi-kashtiban-shayesteh
   :exponential win/exponential
   :hanning-poisson win/hanning-poisson
   :gaussian win/gaussian
   :parzen-exponential win/parzen-exponential
   :dolph-chebyshev win/dolph-chebyshev
   :taylor win/taylor
   :cauchy win/cauchy
   :parzen-geometric win/parzen-geometric
   :kaiser-bessel win/kaiser-bessel
   :cosh win/cosh
   :avci-nacaroglu win/avci-nacaroglu
   :knab win/knab
   :ultraspherical win/ultraspherical
   :saramaki win/saramaki
   :legendre win/legendre
   :bessel-I1 win/bessel-I1
   :shayesteh-kashtiban win/shayesteh-kashtiban
   :kaiser-bessel-derived win/kaiser-bessel-derived
   :vorbis win/vorbis
   :flat-top win/flat-top
   :flat-top-3 win/flat-top-3
   ;; additions (WP 5.1 - 5.4), tested in detail in `fastmath.kernel.window-additions-test`
   :cosine win/cosine
   :kaiser win/kaiser
   :general-gaussian win/general-gaussian
   :flat-top-hp win/flat-top-hp
   :flat-top-srs win/flat-top-srs
   :flat-top-bk win/flat-top-bk
   :hft70 win/hft70
   :hft90d win/hft90d
   :hft95 win/hft95
   :hft116d win/hft116d
   :hft144d win/hft144d
   :hft169d win/hft169d
   :hft196d win/hft196d
   :hft223d win/hft223d
   :hft248d win/hft248d})

(defn- ex-data-of
  "`ex-data` of the `ExceptionInfo` thrown by `(f)`, nil when nothing is thrown; other exceptions are not caught."
  [f]
  (try (f) nil (catch clojure.lang.ExceptionInfo e (or (ex-data e) {:no-data true}))))

(defn- finite-seq?
  [xs]
  (every? #(m/valid-double? (double %)) xs))

(defn- max-err
  "Maximum of `abs(a - b) / max(1, abs(b))`, `##Inf` for different lengths."
  ^double [a b]
  (let [a (vec a)
        b (vec b)]
    (if (= (count a) (count b))
      (reduce m/max 0.0 (map (fn [x y] (m// (m/abs (m/- (double x) (double y))) (m/max 1.0 (m/abs (double y))))) a b))
      ##Inf)))

(deftest registry-list
  (testing "the hand-written list covers all windows of the namespace"
    (is (= 76 (count window-fns)))
    (doseq [[k f] window-fns]
      (is (identical? f (deref (ns-resolve 'fastmath.kernel.window (symbol (name k))))) (str k)))
    (let [public-windows (->> (ns-publics 'fastmath.kernel.window)
                              (filter (fn [[_ v]] (some #{'[N opts]} (:arglists (meta v)))))
                              (map (comp keyword name key))
                              set)]
      (is (= public-windows (set (keys window-fns)))))))

;; ---------------------------------------------------------------------------------------------
;; WP 3.1: N contract

(deftest n-contract
  (doseq [[k w] window-fns
          call [(fn [N] (w N)) (fn [N] (w N nil)) (fn [N] (w N {}))]]
    (testing (str k)
      (testing "N = -1 throws ex-info"
        (let [d (ex-data-of #(call -1))]
          (is (= :invalid-input (:type d)))
          (is (= :N (:property d)))
          (is (= k (:function d)))
          (is (= -1 (:value d)))
          (is (some? (:expected d)))))
      (testing "N = 0 is empty, N = 1 is (1.0)"
        (is (= () (call 0)))
        (is (= [1.0] (vec (call 1)))))
      (testing "length, finite values and symmetry for N in 2 3 8 9"
        (doseq [N [2 3 8 9]]
          (let [r (vec (call N))]
            (is (= N (count r)) (str N))
            (is (finite-seq? r) (str N))
            (is (< (max-err r (reverse r)) 1.0e-12) (str N)))))))
  (testing "invalid N values of every kind throw ex-info"
    (doseq [[k w] window-fns
            N [3.7 nil "x" :a ##NaN ##Inf -0.5 [3]]]
      (let [d (ex-data-of #(w N {}))]
        (is (= :invalid-input (:type d)) (str k " " (pr-str N)))
        (is (= :N (:property d)) (str k " " (pr-str N))))))
  (testing "an integer valued double is accepted"
    (is (= (win/hann 4) (win/hann 4.0)))
    (is (= 0 (count (win/hann 0.0))))))

(deftest n-contract-generic-helpers
  (testing "sample-window"
    (is (= () (win/sample-window m/sin 0)))
    (is (= [1.0] (vec (win/sample-window m/cos 1))))
    (is (= :N (:property (ex-data-of #(win/sample-window m/sin -1)))))
    (is (= :N (:property (ex-data-of #(win/sample-window m/sin 2.5))))))
  (testing "continuous->window-fn and discrete->window-fn"
    (let [cw (win/continuous->window-fn (fn ^double [^double x] (m/- 1.0 (m/* 2.0 (m/abs x)))))
          dw (win/discrete->window-fn (fn [N] (map inc (range N))))]
      (is (= () (cw 0)))
      (is (= () (cw 0 {:normalize? :LInf})))
      (is (= [1.0] (vec (cw 1))))
      (is (= :N (:property (ex-data-of #(cw -1)))))
      (is (= () (dw 0)))
      (is (= [1.0] (vec (dw 1))))
      (is (= [0.5 1.0] (vec (dw 2))))
      (is (= :N (:property (ex-data-of #(dw -1))))))))

;; ---------------------------------------------------------------------------------------------
;; WP 3.2: :normalize? contract

(deftest normalize-coefficients-contract
  (testing "values"
    (is (= [0.5 1.0] (vec (win/normalize-coefficients [1.0 2.0] true))))
    (is (= [0.5 1.0] (vec (win/normalize-coefficients [1.0 2.0] :LInf))))
    (is (= [0.5 1.0] (vec (win/normalize-coefficients [1.0 2.0] :Linf))))
    (is (< (max-err [(/ -1.0 3.0) (/ 2.0 3.0)] (win/normalize-coefficients [-1.0 2.0] :L1)) 1.0e-15)
        ":L1 divides by the sum of the absolute values")
    (is (< (max-err [0.6 0.8] (win/normalize-coefficients [3.0 4.0] :L2)) 1.0e-15))
    (is (< (max-err [(/ 2.0 3.0) (/ 4.0 3.0)] (win/normalize-coefficients [1.0 2.0] :N)) 1.0e-15)
        ":N is :L1 times the number of coefficients")
    (is (= [0.5 1.0] (vec (win/normalize-coefficients [1.0 2.0] 2))))
    (is (= [1.0 2.0] (vec (win/normalize-coefficients [1.0 2.0] false))))
    (is (= [-1.0 -2.0] (vec (win/normalize-coefficients [1.0 2.0] -1.0))) "a negative divisor is a finite non-zero number"))
  (testing "invalid values throw ex-info"
    (doseq [bad [nil :bogus "x" 0 0.0 ##NaN ##Inf ##-Inf :linf :l1 [1] {}]]
      (let [d (ex-data-of #(win/normalize-coefficients [1.0 2.0] bad))]
        (is (= :invalid-input (:type d)) (pr-str bad))
        (is (= :normalize? (:property d)) (pr-str bad)))))
  (testing "empty input never throws for a valid mode"
    (doseq [mode [true false :LInf :Linf :L1 :L2 :N 2.0]]
      (is (= [] (vec (win/normalize-coefficients [] mode))) (pr-str mode))
      (is (= [] (vec (win/normalize-coefficients () mode))) (pr-str mode)))))

(deftest normalize-option-contract
  (doseq [[k w] window-fns]
    (testing (str k)
      (testing "invalid values throw, also for N = 0"
        (doseq [bad [:bogus "x" 0 0.0 ##NaN ##Inf :linf [1]]
                N [0 1 5]]
          (let [d (ex-data-of #(w N {:normalize? bad}))]
            (is (= :invalid-input (:type d)) (str N " " (pr-str bad)))
            (is (= :normalize? (:property d)) (str N " " (pr-str bad))))))
      (testing "valid modes accepted for N = 0, 1 and 5"
        (doseq [mode [true false :LInf :Linf :L1 :L2 :N 2.0]]
          (is (= () (w 0 {:normalize? mode})) (pr-str mode))
          (is (= 1 (count (w 1 {:normalize? mode}))) (pr-str mode))
          (is (= 5 (count (w 5 {:normalize? mode}))) (pr-str mode))))
      (testing "nil is the default of the window"
        (is (= (vec (w 8)) (vec (w 8 {:normalize? nil})))))
      (testing ":Linf is an alias of :LInf, the modes scale the same shape"
        (is (= (vec (w 8 {:normalize? :LInf})) (vec (w 8 {:normalize? :Linf}))))
        (let [raw (vec (w 9 {:normalize? false}))
              l1 (vec (w 9 {:normalize? :L1}))
              nn (vec (w 9 {:normalize? :N}))
              n2 (vec (w 9 {:normalize? 4.0}))]
          (is (< (max-err l1 (map #(/ % (reduce + (map m/abs raw))) raw)) 1.0e-12))
          (is (< (max-err nn (map #(* 9.0 %) l1)) 1.0e-12))
          (is (< (max-err n2 (map #(/ % 4.0) raw)) 1.0e-12)))))))

(deftest normalize-option-values
  (testing "documented examples"
    (is (= () (win/hann 0 {:normalize? :LInf})))
    (is (= [0.0 0.5 1.0 0.5 0.0] (vec (win/hann 5 {:normalize? nil}))))
    (is (= (win/hann 7) (win/hann 7 {:normalize? true})))
    (is (= :normalize? (:property (ex-data-of #(win/hann 5 {:normalize? :bogus}))))))
  (testing "true is the division by f(0) for continuous model windows (even N: maximum below 1)"
    (is (< (apply max (win/hann 4)) 1.0))
    (is (m/delta-eq 0.75 (apply max (win/hann 4)) 1.0e-14))
    (is (= 1.0 (apply max (win/hann 5)))))
  (testing "true is :LInf for discrete windows"
    (is (= 1.0 (apply max (win/dolph-chebyshev 8))))
    (is (= 1.0 (apply max (win/saramaki 8)))))
  (testing "triangular and parzen discrete default paths honour every mode"
    (doseq [[f opts] [[win/triangular {}] [win/triangular {:shift 3}] [win/parzen {}]]
            N [4 5 8]]
      (let [base (vec (f N (assoc opts :normalize? false)))
            sum (reduce + base)]
        (is (= base (vec (f N opts))) "peak 1 formula: true and false agree")
        (is (= base (vec (f N (assoc opts :normalize? true)))))
        (is (< (max-err (f N (assoc opts :normalize? :L1)) (map #(/ % sum) base)) 1.0e-14) (str N))
        (is (m/delta-eq 1.0 (m/sqrt (reduce + (map m/sq (f N (assoc opts :normalize? :L2))))) 1.0e-14))
        (is (< (max-err (f N (assoc opts :normalize? :N)) (map #(* N (/ % sum)) base)) 1.0e-14))
        (is (< (max-err (f N (assoc opts :normalize? 4.0)) (map #(/ % 4.0) base)) 1.0e-14))
        (is (< (max-err (f N (assoc opts :normalize? :LInf)) (map #(/ % (apply max base)) base)) 1.0e-14))))
    (is (= [0.25 0.75 0.75 0.25] (vec (win/triangular 4 {:normalize? false}))))
    (is (m/delta-eq 1.0 (reduce + (win/triangular 4 {:normalize? :L1})) 1.0e-14))
    (is (not= (vec (win/parzen 8)) (vec (win/parzen 8 {:normalize? :L1}))))))

(deftest discrete-selection-options
  (testing ":discrete? must be a boolean, :shift an integer >= 0 (nil is the default)"
    (is (= (vec (win/parzen 8)) (vec (win/parzen 8 {:discrete? nil}))))
    (is (= :discrete? (:property (ex-data-of #(win/parzen 8 {:discrete? "x"})))))
    (is (= (vec (win/triangular 8)) (vec (win/triangular 8 {:shift nil}))))
    (is (= :shift (:property (ex-data-of #(win/triangular 8 {:shift 1.5})))))))

;; ---------------------------------------------------------------------------------------------
;; WP 3.3: parameter domains.
;; window, parameter, accepted values (at and inside the bound), rejected values (just outside).
;;
;; Domain table:
;;   triangular :shift int >= 0                      parzen :discrete? boolean
;;   b-spline :order int [1, 1000]                   singla-singh :order int [0, 1000]
;;   connes :alpha (0, Inf)                          parzen-algebraic :gamma (0, 1], :u (0, Inf)
;;   lanczos :L [0, Inf)                             raised-cosine :alpha [0.5, 1]
;;   webster-hamming :v (-0.5, 1000]                 power-of-cosine :m [0, 1000]
;;   raised-power-of-cosine :alpha [0, 1], :m [0, 1000]
;;   parzen-cosine :gamma [0, 1], :m (0, Inf)        trapezoid :alpha [0, 0.5]
;;   tukey :alpha [0, Inf) (alpha > 1 is hann)
;;   :alpha (0, Inf): exponential, hanning-poisson, gaussian, cauchy, kaiser-bessel, cosh,
;;                    avci-nacaroglu, knab, bessel-I1, kaiser-bessel-derived, parzen-exponential,
;;                    parzen-geometric (also :r (0, Inf))
;;   :level finite, |level| <= 5000 (sign ignored): dolph-chebyshev, saramaki, legendre, taylor,
;;                    ultraspherical (taylor: a level giving a non-positive middle, e.g. 0, throws)
;;   ultraspherical :alpha (-1, 150]                 taylor :n int [2, 400]
;;   kaiser :beta [0, 1e300]                         general-gaussian :p (0, 1000], :sig (0, Inf)
;;   blackman-harris-family :coeffs non-empty sequence of finite numbers, first > 0
;;   mottaghi-kashtiban-shayesteh, flat-top, flat-top-3 and the fixed windows have no parameter.
;; An open bound is rejected, a closed bound accepted; NaN and +-Inf are rejected everywhere.

(def ^:private alpha-ok [1.0e-3 0.5 1.0 10.0 100.0])
(def ^:private alpha-bad [0 0.0 -1 -1.0e-9 ##NaN ##Inf ##-Inf])
(def ^:private level-ok [-5000.0 -100.0 -50.0 -30.0 30.0 50.0 5000.0])
(def ^:private level-bad [-5000.1 5000.1 -1.0e6 ##NaN ##Inf ##-Inf "x"])

(def ^:private domain-table
  (vec (concat
        [[:triangular :shift [0 1 2 7] [-1 0.5 1.5 ##NaN ##Inf "x"]]
         [:b-spline :order [1 2 3 60 1000 4.0] [0 -1 1001 2.5 ##NaN ##Inf "x"]]
         [:singla-singh :order [0 1 20 1000 3.0] [-1 1001 1.5 ##NaN 1.0e10 ##Inf]]
         [:connes :alpha alpha-ok alpha-bad]
         [:parzen-algebraic :gamma [1.0e-3 0.5 1.0] [0 -0.1 1.0000001 2 ##NaN]]
         [:parzen-algebraic :u [1.0e-3 1.0 3.0 100.0] alpha-bad]
         [:lanczos :L [0 0.5 1.0 3.0 100.0] [-1.0e-4 -1 ##Inf ##NaN]]
         [:raised-cosine :alpha [0.5 0.75 1.0] [0.4999999 1.0000001 0 1.5 ##NaN]]
         [:webster-hamming :v [-0.49 -0.4 0 1.0 1000.0] [-0.5 -1 1000.001 ##NaN ##Inf]]
         [:power-of-cosine :m [0 0.5 1.0 2.0 1000.0] [-0.001 -1 1000.001 ##NaN ##Inf]]
         [:raised-power-of-cosine :alpha [0 0.05 0.5 1.0] [-0.001 1.001 ##NaN]]
         [:raised-power-of-cosine :m [0 1.0 1000.0] [-0.001 1000.001 ##NaN]]
         [:parzen-cosine :gamma [0 0.5 1.0] [-0.001 1.001 ##NaN]]
         [:parzen-cosine :m [1.0e-3 1.0 2.0 50.0] alpha-bad]
         [:trapezoid :alpha [0 0.1 0.25 0.5] [-0.001 0.5000001 1 ##NaN]]
         [:tukey :alpha [0 0.1 1.0 1.5 100.0] [-0.001 ##NaN ##Inf]]
         [:parzen-exponential :r [0.25 0.5 1.0 5.0] alpha-bad]
         [:parzen-geometric :r [0.25 0.5 1.0 5.0] alpha-bad]
         [:ultraspherical :alpha [-0.99 -0.5 0 0.5 2.0 150.0] [-1 -1.0001 150.001 ##NaN ##Inf]]
         [:taylor :n [2 4 400 8.0] [1 0 401 2.5 ##NaN]]
         [:kaiser :beta [0 1.0e-300 0.5 3.0 100.0 1.0e5 1.0e300] [-1.0e-9 -1 1.0e301 ##NaN ##Inf ##-Inf "x"]]
         [:general-gaussian :p [1.0e-300 1.0e-3 0.5 1.0 10.0 1000.0] [0 0.0 -1 1000.001 ##NaN ##Inf ##-Inf "x"]]
         [:general-gaussian :sig [1.0e-3 0.05 0.5 1.0 100.0 1.0e100] alpha-bad]
         [:blackman-harris-family :coeffs [[1.0] [0.5 0.5] [1.0 1.5 0.7 -0.3] '(1 2)]
          [[] [0] [-1 1] [##NaN 1] [1 ##Inf] [1 "x"] 1 "x" {}]]]
        (for [k [:exponential :hanning-poisson :gaussian :cauchy :kaiser-bessel :cosh :avci-nacaroglu
                 :knab :bessel-I1 :kaiser-bessel-derived :parzen-exponential :parzen-geometric]]
          [k :alpha alpha-ok alpha-bad])
        (for [k [:dolph-chebyshev :saramaki :legendre :ultraspherical]]
          [k :level level-ok level-bad])
        [[:taylor :level (remove #{-5000.0 5000.0} (conj level-ok -20.0)) level-bad]
         [:taylor :level [-5000.0 5000.0] level-bad]])))

(deftest parameter-domains
  (doseq [[k param ok bad] domain-table
          :let [w (window-fns k)]]
    (testing (str k " " param)
      (doseq [v ok]
        (let [r (vec (w 9 {param v}))]
          (is (= 9 (count r)) (pr-str v))
          (is (finite-seq? r) (str "no NaN / Inf for " (pr-str v)))))
      (doseq [v bad]
        (let [d (ex-data-of #(w 9 {param v}))]
          (is (= :invalid-input (:type d)) (pr-str v))
          (is (= param (:property d)) (pr-str v))
          (is (= k (:function d)) (pr-str v))
          (is (some? (:expected d)) (pr-str v))
          (when-not (and (number? v) (m/nan? (double v)))
            (is (= v (:value d)) (pr-str v)))))
      (testing "nil is a missing parameter"
        (is (= (vec (w 9)) (vec (w 9 {param nil}))))))))

(deftest parameter-domain-details
  (testing "the sign of :level is ignored"
    (doseq [k [:dolph-chebyshev :saramaki :legendre :ultraspherical :taylor]]
      (is (= (vec ((window-fns k) 9 {:level -50.0})) (vec ((window-fns k) 9 {:level 50.0}))) (str k))))
  (testing "taylor: a level with a non-positive window middle throws on :level"
    (is (= :level (:property (ex-data-of #(win/taylor 9 {:level 0.0}))))))
  (testing "raised-cosine: the old code clipped a negative lobe and returned a window for alpha 0.25"
    (is (= :alpha (:property (ex-data-of #(win/raised-cosine 9 {:alpha 0.25}))))))
  (testing "kaiser-bessel-derived: an alpha with no non-zero kernel sample throws instead of NaN"
    (is (= :alpha (:property (ex-data-of #(win/kaiser-bessel-derived 9 {:alpha 1.0e5}))))))
  (testing "defaults"
    (is (= (vec (win/taylor 9)) (vec (win/taylor 9 {:level -50.0 :n 4}))))
    (is (= (vec (win/tukey 9)) (vec (win/tukey 9 {:alpha 0.5}))))
    (is (= (vec (win/ultraspherical 9)) (vec (win/ultraspherical 9 {:level -50.0 :alpha 2.0}))))
    (is (= (vec (win/singla-singh 9)) (vec (win/singla-singh 9 {:order 1}))))
    (is (= (vec (win/b-spline 9)) (vec (win/b-spline 9 {:order 3}))))))

;; ---------------------------------------------------------------------------------------------
;; Oracle comparison: test/resources/kernel/window_reference.edn (utils/fastmath/dev/generate_window_reference.py,
;; mpmath 60 digits and scipy) and window_octave_reference.edn (generate_window_octave.m).

(def ^:private reference (delay (edn/read-string (slurp (io/resource "kernel/window_reference.edn")))))
(def ^:private octave-reference (delay (edn/read-string (slurp (io/resource "kernel/window_octave_reference.edn")))))

(defn- case-error
  "Maximum error (see `max-err`) of the window of the reference case `c`; the ends are skipped for `:interior`.
  A `:periodic` case is the symmetric window of `N + 1` samples without the last one."
  ^double [c]
  (let [w (window-fns (:window c))
        a (vec (w (if (:periodic c) (inc (long (:N c))) (:N c)) (assoc (:opts c) :normalize? (:normalize? c))))
        a (if (:periodic c) (pop a) a)
        e (vec (:expected c))
        trim (fn [v] (if (:interior c) (subvec v 1 (dec (count v))) v))]
    (max-err (trim a) (trim e))))

(defn- check-group
  "Asserts that all cases of `group` (of `ref`) are within `tol`; `n` is the expected number of cases."
  [ref group tol n]
  (let [cs (filter #(and (= group (:group %)) (window-fns (:window %))) (:cases @ref))]
    (is (= n (count cs)) (str group " number of cases"))
    (doseq [c cs]
      (let [e (case-error c)]
        (is (<= e tol) (str (:id c) " error " e))))))

;; Observed maximum errors (relative to max(1, |expected|)) are in the comments.

(deftest b-spline-reference
  ;; WP 4.1, W1-01, W1-05: raw scale order * B(order x), stable up to order 1000. Observed 8.9e-16.
  (check-group reference :b-spline 1.0e-13 48)
  (testing "no NaN for large orders (the old truncated power sum gave NaN from order 101)"
    (doseq [order [101 150 200 500 1000]]
      (let [w (vec (win/b-spline 9 {:order order :normalize? false}))]
        (is (finite-seq? w) (str order))
        (is (< (max-err (win/b-spline 9 {:order order}) (map #(/ % (nth w 4)) w)) 1.0e-12))))
    (is (m/delta-eq 1.0 (nth (vec (win/b-spline 9 {:order 200})) 4) 1.0e-14))))

(deftest tukey-reference
  ;; WP 4.1, W1-02: plateau 2 / (2 - alpha). Observed 4.4e-16.
  (check-group reference :tukey 1.0e-14 24)
  (is (m/delta-eq (/ 2.0 1.5) (nth (vec (win/tukey 9 {:alpha 0.5 :normalize? false})) 4) 1.0e-14)))

(deftest webster-hamming-reference
  ;; WP 4.1, W1-03: A uses Gamma((4 + v) / 2). Observed 3.7e-16 (v < 0: the two singular ends are skipped).
  (check-group reference :webster-hamming 1.0e-13 12)
  (testing "v = 0 raw peak equals hamming-exact (the old code gave twice as much)"
    (is (m/delta-eq (nth (vec (win/hamming-exact 9 {:normalize? false})) 4)
                    (nth (vec (win/webster-hamming 9 {:v 0.0 :normalize? false})) 4)
                    1.0e-14))))

(deftest quadrature-constants-reference
  ;; WP 4.3, W2-10: avci-nacaroglu, parzen-exponential and parzen-geometric (r < 1) raw scale. Observed 7.9e-16.
  ;; The old gk-quadrature constants were off by 1e-9.
  (check-group reference :quadrature 1.0e-13 26))

(deftest knab-reference
  ;; WP 4.4, W2-03: the ends have the limit value pi alpha / sinh(pi alpha). Observed 6.2e-16.
  (check-group reference :knab 1.0e-13 24)
  (testing "documented values"
    (is (m/delta-eq 0.0460 (first (win/knab 5 {:normalize? false})) 1.0e-4))
    (is (not= [0.0 0.0] (vec (win/knab 2))))
    (is (m/delta-eq (/ (* m/PI 2.0) (m/sinh (* m/PI 2.0))) (first (win/knab 2)) 1.0e-14)))
  (testing "knab-fn and avci-nacaroglu-fn are private"
    (is (not (contains? (ns-publics 'fastmath.kernel.window) 'knab-fn)))
    (is (not (contains? (ns-publics 'fastmath.kernel.window) 'avci-nacaroglu-fn)))))

(deftest low-sidelobe-removed
  ;; WP 4.5, W2-01: the window overflowed the stack for every N.
  (is (not (contains? (ns-publics 'fastmath.kernel.window) 'low-sidelobe)))
  (is (not (contains? (ns-publics 'fastmath.kernel.window) 'low-sidelobe-continuous))))

(deftest vorbis-reference
  ;; WP 4.6, W2-09: sin(pi/2 sin^2(pi (n + 1/2) / N)). Observed 4.4e-16.
  (check-group reference :vorbis 1.0e-14 16)
  (testing "default :normalize? is false"
    (is (= (vec (win/vorbis 8)) (vec (win/vorbis 8 {:normalize? false}))))
    (is (= 1.0 (apply max (win/vorbis 9 {:normalize? true})))))
  (testing "vorbis-continuous is kept"
    (is (m/delta-eq 1.660592492619163 (win/vorbis-continuous 0.0) 1.0e-14))
    (is (< (win/vorbis-continuous 0.5) 1.0e-30))))

(deftest blackman-harris-family-reference
  ;; WP 4.7, W2-02: the family is not clipped; mpmath cosine sum. Observed 5.0e-16.
  (check-group reference :family 1.0e-14 24)
  ;; scipy.signal.windows.flattop (general_cosine) on the same coefficients. Observed 7.2e-16.
  (check-group reference :family-scipy 1.0e-14 5)
  (testing "the negative lobes of the flat-top coefficients are kept (the old code clipped them to 0)"
    (let [flat [0.21557895 0.41663158 0.277263158 0.083578947 0.006947368]
          w (vec (win/blackman-harris-family 65 {:coeffs flat}))]
      (is (neg? (apply min w)))
      (is (m/delta-eq -0.0004 (first w) 1.0e-4))))   ; scipy flattop end value -4.21e-4
  (testing "->blackman-harris: the clip? flag; the named positive windows clip"
    (let [flat [0.21557895 0.41663158 0.277263158 0.083578947 0.006947368]
          xs (range -0.5 0.5001 0.01)
          unclipped (win/->blackman-harris flat false)
          clipped (win/->blackman-harris flat)
          clipped-explicit (win/->blackman-harris flat true)]
      (is (neg? (apply min (map unclipped xs))))
      (is (= 0.0 (apply min (map clipped xs))))
      (is (= (map clipped xs) (map clipped-explicit xs)))
      (is (every? #(>= (win/blackman-harris-92db-continuous %) 0.0) xs))))
  (testing "empty :coeffs and non-positive a0 throw instead of NPE / NaN"
    (is (= :coeffs (:property (ex-data-of #(win/blackman-harris-family 5 {:coeffs []})))))
    (is (= :coeffs (:property (ex-data-of #(win/blackman-harris-family 5 {:coeffs [0.0 1.0]})))))))

(deftest flat-top-reference
  ;; WP 4.7, W2-08: flat-top-3 = FTNI of Heinzel et al. 2002 (0.2810639 0.5208972 0.1980399) and flat-top (the
  ;; 5-term scipy `flattop` set) keep negative lobes. Observed 1.4e-15.
  (check-group reference :flat-top-3 1.0e-14 20)
  (is (neg? (apply min (win/flat-top 65))))
  (is (neg? (apply min (win/flat-top-3 65)))))

(deftest mottaghi-kashtiban-shayesteh-reference
  ;; WP 4.8, W2-13: unclipped published formula for N >= 2, (1.0) for N = 1. Observed 4.1e-16.
  (check-group reference :mks 1.0e-14 20)
  (is (= [1.0] (vec (win/mottaghi-kashtiban-shayesteh 1))))
  (testing "N 2..4 have negative end values, N >= 5 is non-negative"
    (doseq [N [2 3 4]]
      (is (neg? (first (win/mottaghi-kashtiban-shayesteh N))) (str N)))
    (doseq [N [5 6 8 9 16 65 256]]
      (is (every? #(>= % 0.0) (win/mottaghi-kashtiban-shayesteh N)) (str N)))))

(deftest singla-singh-reference
  ;; WP 4.2, W1-04: 2 (1 - I_t(m + 1, m + 1)) by the regularized beta function. Observed 1.2e-14 (order 100)
  ;; for orders 0..1000; the old alternating sum had error 1.0 from order 21.
  (check-group reference :singla-singh 1.0e-12 64)
  (testing "an integer valued double order is accepted, a huge order throws immediately"
    (is (= (vec (win/singla-singh 9 {:order 3})) (vec (win/singla-singh 9 {:order 3.0}))))
    (let [t0 (System/nanoTime)
          d (ex-data-of #(win/singla-singh 9 {:order 1.0e10}))]
      (is (= :order (:property d)))
      (is (< (/ (- (System/nanoTime) t0) 1.0e9) 1.0)))))

(deftest taylor-reference
  ;; WP 4.9, W2-06: Doerry eq. 208-212 on the fastmath grid (endpoints included, differs from scipy). Observed 5.6e-16.
  (check-group reference :taylor 1.0e-14 60)
  (testing "level sign is ignored and :n is the number of nearly constant sidelobes"
    (is (= (vec (win/taylor 9 {:level 50.0 :n 4})) (vec (win/taylor 9 {:level -50.0 :n 4}))))
    (is (not= (vec (win/taylor 9 {:n 2})) (vec (win/taylor 9 {:n 4}))))))

(deftest discrete-windows-reference
  ;; WP 4.1, W2-11: dolph-chebyshev, ultraspherical, saramaki, legendre, shayesteh-kashtiban from the spectrum
  ;; definition in mpmath, raw scale = mean 1. Observed 1.6e-14 (dolph-chebyshev N 33 level -80, the shared
  ;; accuracy limit 1e-16 N^2 of T_(N-1) near 1).
  (check-group reference :discrete 1.0e-13 78))

(deftest scipy-reference
  ;; Observed maxima: chebwin 5.6e-14 (N 65, level 40: the accuracy limit shared by fastmath and scipy,
  ;; Findings W2-05), kaiser 7.8e-16, kaiser_bessel_derived 1.1e-16, gaussian 2.2e-16, tukey 3.2e-15.
  (check-group reference :scipy 1.0e-13 64))

(deftest octave-reference-test
  ;; Octave 11 signal 1.4.7: chebwin (observed 5.6e-14), nuttallwin, blackmanharris, kaiser (kaiser-bessel and
  ;; kaiser), gausswin and flattopwin (flat-top-bk, N 1..40, 63..65, 100, 127..129, 255, 256, symmetric and periodic).
  (let [cs (filter #(window-fns (:window %)) (:cases @octave-reference))]
    (is (= 168 (count cs)) "all cases are used")
    (doseq [c (remove #(#{:flat-top-bk :kaiser} (:window %)) cs)]
      (let [e (case-error c)]
        (is (<= e 1.0e-13) (str (:id c) " error " e))))
    ;; Observed maxima: flat-top-bk 1.2e-15 (N 256 and N 255 periodic: Octave computes cos(2 pi j k / N) with a rounded
    ;; argument and is itself 1.0e-15 away from a 50 digit evaluation, fastmath 5.6e-16), kaiser 1.2e-15.
    (doseq [c (filter #(= :flat-top-bk (:window %)) cs)]
      (let [e (case-error c)]
        (is (<= e 2.0e-15) (str (:id c) " error " e))))
    (doseq [c (filter #(= :kaiser (:window %)) cs)]
      (let [e (case-error c)]
        (is (<= e 1.0e-14) (str (:id c) " error " e))))))
