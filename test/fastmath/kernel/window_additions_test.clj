(ns fastmath.kernel.window-additions-test
  "Tests of the windows added by the kernel audit (WP 5.1 - 5.4: `cosine`, `kaiser`, `general-gaussian` and the
  twelve flat-top sets), of the achieved sidelobe levels documented for the `:level` windows (WP 6.3) and the
  scripted documentation check of every public var of `fastmath.kernel.window` (WP 6.1, 6.2).
  The N / `:normalize?` / parameter contracts of the new windows are in `fastmath.kernel.window-test`
  (they are listed in `window-fns` there)."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel.window :as win]
            [fastmath.special :as spec]
            [fastmath.transform :as t]))

(def ^:private additions
  {:exponential win/exponential     ; not an addition: the scipy mapping stated in its docstring
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

(def ^:private flat-top-sets [:flat-top-hp :flat-top-srs :flat-top-bk :hft70 :hft90d :hft95 :hft116d :hft144d
                              :hft169d :hft196d :hft223d :hft248d])

(defn- max-err
  "Maximum of `abs(a - b) / max(1, abs(b))`, `##Inf` for different lengths."
  ^double [a b]
  (let [a (vec a)
        b (vec b)]
    (if (= (count a) (count b))
      (reduce m/max 0.0 (map (fn [x y] (m// (m/abs (m/- (double x) (double y))) (m/max 1.0 (m/abs (double y))))) a b))
      ##Inf)))

(defn- ex-data-of
  [f]
  (try (f) nil (catch clojure.lang.ExceptionInfo e (or (ex-data e) {:no-data true}))))

(defn- simpson
  "Composite Simpson rule of `f` on `[a, b]` with `n` (even) intervals."
  ^double [f ^double a ^double b ^long n]
  (let [h (/ (- b a) n)]
    (* (/ h 3.0)
       (+ (double (f a)) (double (f b))
          (reduce + 0.0 (map (fn [^long i] (* (if (odd? i) 4.0 2.0) (double (f (+ a (* i h)))))) (range 1 n)))))))

;; ---------------------------------------------------------------------------------------------
;; Reference data: window_reference.edn (generate_window_reference.py: mpmath 60 digits, scipy) and
;; window_octave_reference.edn (generate_window_octave.m: Octave signal 1.4.7)

(def ^:private reference (delay (edn/read-string (slurp (io/resource "kernel/window_reference.edn")))))
(def ^:private octave-reference (delay (edn/read-string (slurp (io/resource "kernel/window_octave_reference.edn")))))
(def ^:private psl-reference (delay (edn/read-string (slurp (io/resource "kernel/flattop_psl_reference.edn")))))

(defn- case-error
  "Maximum error of the window of the reference case `c`; a `:periodic` case is the symmetric window of `N + 1`
  samples without the last one."
  ^double [c]
  (let [w (additions (:window c))
        N (:N c)
        a (vec (w (if (:periodic c) (inc (long N)) N) (assoc (:opts c) :normalize? (:normalize? c))))]
    (max-err (if (:periodic c) (pop a) a) (:expected c))))

(defn- group-errors
  [ref group]
  (let [cs (filter #(and (= group (:group %)) (additions (:window %))) (:cases @ref))]
    {:n (count cs)
     :max (reduce m/max 0.0 (map case-error cs))
     :worst (:id (apply max-key case-error cs))}))

(defn- check-group
  "All `n` cases of `group` are within `tol`."
  [group tol n]
  (let [{mx :max worst :worst cnt :n} (group-errors reference group)]
    (is (= n cnt) (str group " number of cases"))
    (is (<= mx tol) (str group " error " mx " (" worst ")"))))

(deftest cosine-reference
  ;; WP 5.1. mpmath sin(pi (n + 1/2) / N), observed 3.9e-16; scipy.signal.windows.cosine symmetric and periodic
  ;; (sym=False), observed 4.4e-16.
  (check-group :cosine 1.0e-15 22)
  (check-group :cosine-scipy 1.0e-15 21))

(deftest kaiser-reference
  ;; WP 5.2. scipy.signal.windows.kaiser (symmetric and periodic), observed 5.3e-15 (beta 100, the Bessel function
  ;; I0 of both libraries differs by a few ulps at large arguments); raw scale from the mpmath definition, observed
  ;; 1.5e-15.
  (check-group :kaiser-scipy 1.0e-14 96)
  (check-group :kaiser-mp 1.0e-14 48)
  (testing "Octave kaiser(N, beta), observed 1.2e-15"
    (let [cs (filter #(= :kaiser (:window %)) (:cases @octave-reference))]
      (is (= 24 (count cs)))
      (doseq [c cs]
        (is (<= (case-error c) 1.0e-14) (:id c))))))

(deftest general-gaussian-reference
  ;; WP 5.3. scipy.signal.windows.general_gaussian(N, p, sig * (N - 1)) (periodic: sig * N), observed 1.1e-15;
  ;; raw scale (integral 1) from the mpmath definition, observed 4.4e-16.
  (check-group :general-gaussian-scipy 2.0e-15 90)
  (check-group :general-gaussian-mp 2.0e-15 112))

(deftest documented-scipy-mapping-reference
  ;; The mapping stated in the docstring of `exponential`: scipy exponential(N, tau = (N - 1) / (2 alpha)), observed 4.4e-16.
  (check-group :exponential-scipy 1.0e-15 18))

(defn- same?
  "Windows `a` and `b` (functions of `[N opts]`) agree for N in 9 16 in both scales within `tol`."
  [a opts-a b opts-b tol]
  (every? (fn [[N normalize?]]
            (< (max-err (a N (assoc opts-a :normalize? normalize?)) (b N (assoc opts-b :normalize? normalize?))) tol))
          (for [N [9 16] n [true false]] [N n])))

(deftest documented-identities
  ;; Identities stated in the docstrings. Observed maxima: below 2e-15.
  (testing "cos^4 and cos^6 are the Nuttall windows with the third and fifth derivative continuous"
    (is (same? win/power-of-cosine {:m 4.0} win/nutall-3-3rd {} 1.0e-13))
    (is (same? win/power-of-cosine {:m 6.0} win/nutall-5th {} 1.0e-13)))
  (testing "welch is parzen-algebraic gamma 1 u 2, power-of-cosine m 0 and lanczos L 0 are rectangular"
    (is (same? win/welch {} win/parzen-algebraic {:gamma 1.0 :u 2.0} 1.0e-14))
    (is (same? win/power-of-cosine {:m 0.0} win/rectangular {} 1.0e-14))
    (is (same? win/lanczos {:L 0.0} win/rectangular {} 1.0e-14)))
  (testing "triangular shift 2 end values: 1/N for even N, 2/(N+1) for odd N"
    (is (m/delta-eq (/ 1.0 8.0) (first (win/triangular 8)) 1.0e-15))
    (is (m/delta-eq (/ 2.0 10.0) (first (win/triangular 9)) 1.0e-15)))
  (testing "the end values stated for the exponential family"
    (is (m/delta-eq (m/exp -1.5) (first (win/exponential 9 {:alpha 1.5})) 1.0e-14))
    (is (m/delta-eq (m/exp (* -0.5 4.0)) (first (win/gaussian 9 {:alpha 2.0})) 1.0e-14))
    (is (m/delta-eq (m/exp (- (* m/PI 2.0))) (first (win/avci-nacaroglu 9 {:alpha 2.0})) 1.0e-14))
    (is (m/delta-eq (/ 1.0 (spec/bessel-I0 (* m/PI 1.5))) (first (win/kaiser-bessel 9 {:alpha 1.5})) 1.0e-14)))
  (testing "the Nuttall four-term 5th derivative window has the sidelobes of cos^6 and the stated ends"
    (is (m/delta-eq 0.0 (first (win/nutall-5th 9)) 1.0e-15)))
  (testing "mottaghi-kashtiban-shayesteh has the stated coefficients (here N = 9)"
    (is (< (max-err (win/mottaghi-kashtiban-shayesteh 9)
                    (win/blackman-harris-family 9 {:coeffs [(- 0.5363 (/ 0.14 8.0)) (- 0.996 (- 0.5363 (/ 0.14 8.0))) 0.0 0.004]}))
           1.0e-14))))

(deftest flat-top-sets-reference
  ;; WP 5.4. mpmath cosine sums with the positive magnitudes of the sets, N 2 3 8 9 64, raw (integral 1) and default
  ;; scale, observed 1.7e-15 (hft196d, hft223d raw: the values are up to 10, the error is relative to max(1, |y|));
  ;; scipy general_cosine, symmetric and periodic, observed 8.9e-16.
  (check-group :flat-top-sets 4.0e-15 120)
  (check-group :flat-top-sets-scipy 2.0e-15 72))

(deftest flat-top-bk-is-octave-flattopwin
  ;; Octave 11 flattopwin(N, 'symmetric' / 'periodic'), N 1 .. 40, 63 .. 65, 100, 127 .. 129, 255, 256: observed
  ;; maximum 1.2e-15 (N 256 symmetric, N 255 periodic). Octave itself is 1.0e-15 away from a 50 digit evaluation
  ;; (it computes cos(2 pi j / N k) with a rounded argument), fastmath 5.6e-16 (N 256).
  (let [cs (filter #(= :flat-top-bk (:window %)) (:cases @octave-reference))]
    (is (= 97 (count cs)))
    (is (= (set (concat (range 1 41) [63 64 65 100 127 128 129 255 256])) (set (map :N cs))))
    (doseq [c cs]
      (is (<= (case-error c) 2.0e-15) (:id c))))
  (testing "N = 1 is (1.0) and the values of Octave flattopwin(2)"
    (is (= [1.0] (vec (win/flat-top-bk 1))))
    (is (m/delta-eq (/ 0.0042 4.6402) (first (win/flat-top-bk 2)) 1.0e-15))))

;; ---------------------------------------------------------------------------------------------
;; cosine, kaiser, general-gaussian: identities and details

(deftest cosine-details
  (testing "the default is the formula, true is the maximum 1"
    (is (= (vec (win/cosine 8)) (vec (win/cosine 8 {:normalize? false}))))
    (is (m/delta-eq (m/cos (/ m/PI 16.0)) (apply max (win/cosine 8)) 1.0e-15) "even N: cos(pi / 2N)")
    (is (= 1.0 (apply max (win/cosine 9))) "odd N")
    (is (= 1.0 (apply max (win/cosine 8 {:normalize? true}))))
    (is (every? pos? (win/cosine 64)) "no zero ends, the first value is sin(pi / 2N)")
    (is (m/delta-eq (m/sin (/ m/PI 128.0)) (first (win/cosine 64)) 1.0e-15)))
  (testing "the other normalisations"
    (is (m/delta-eq 1.0 (reduce + (win/cosine 9 {:normalize? :L1})) 1.0e-14))
    (is (m/delta-eq 1.0 (m/sqrt (reduce + (map m/sq (win/cosine 9 {:normalize? :L2})))) 1.0e-14)))
  (testing "cosine-continuous: integral 1, the discrete window samples it at the half-sample points"
    (is (m/delta-eq 1.0 (simpson win/cosine-continuous -0.5 0.5 2000) 1.0e-12))
    (is (< (max-err (win/cosine 8) (map #(/ (win/cosine-continuous (- (/ (+ % 0.5) 8.0) 0.5)) m/HALF_PI) (range 8)))
           1.0e-15))
    (is (< (max-err (win/power-of-cosine 9 {:m 1.0 :normalize? false})
                    (map win/cosine-continuous (m/sample identity -0.5 0.5 9)))
           1.0e-14)
        "sampled on the inclusive grid it is power-of-cosine m 1, not the cosine window")
    (is (> (max-err (win/cosine 9) (win/power-of-cosine 9)) 0.01))))

(deftest kaiser-details
  (testing "kaiser is kaiser-bessel with alpha = beta / pi, both scales"
    (doseq [beta [0.5 m/PI 8.6 30.0]
            N [8 9 64]
            normalize? [true false]]
      (is (< (max-err (win/kaiser N {:beta beta :normalize? normalize?})
                      (win/kaiser-bessel N {:alpha (/ beta m/PI) :normalize? normalize?}))
             1.0e-15)
          (str beta " " N " " normalize?))))
  (testing "defaults: beta = pi is kaiser-bessel with alpha 1"
    (is (< (max-err (win/kaiser 9) (win/kaiser-bessel 9)) 1.0e-15))
    (is (< (max-err (win/kaiser 16 {:beta nil}) (win/kaiser-bessel 16 {:alpha 1.0})) 1.0e-15)))
  (testing "beta = 0 is the rectangular window, raw and default"
    (is (= (vec (win/rectangular 9)) (vec (win/kaiser 9 {:beta 0}))))
    (is (= (vec (win/rectangular 9 {:normalize? false})) (vec (win/kaiser 9 {:beta 0.0 :normalize? false})))))
  (testing "the raw window has the integral 1 and the end value 1 / I0(beta) of the middle value"
    (is (m/delta-eq 1.0 (simpson (win/->kaiser 8.6) -0.5 0.5 2000) 1.0e-12))
    (is (m/delta-eq (/ 1.0 (spec/bessel-I0 8.6)) (first (win/kaiser 9 {:beta 8.6})) 1.0e-15)))
  (testing "a huge beta gives a finite window with a single non-zero value"
    (is (= [0.0 0.0 1.0 0.0 0.0] (vec (win/kaiser 5 {:beta 1.0e10}))))))

(deftest general-gaussian-details
  (testing "documented defaults"
    (is (= (vec (win/general-gaussian 9)) (vec (win/general-gaussian 9 {:p 1.5 :sig 0.25})))))
  (testing "p = 1 is gaussian with alpha = 1 / (2 sig), both scales"
    (doseq [alpha [0.5 1.0 2.5 4.0]
            N [8 9 33]
            normalize? [true false]]
      (is (< (max-err (win/general-gaussian N {:p 1.0 :sig (/ 1.0 (* 2.0 alpha)) :normalize? normalize?})
                      (win/gaussian N {:alpha alpha :normalize? normalize?}))
             2.0e-14)
          (str alpha " " N " " normalize?))))
  (testing "parzen-exponential with r = 2 p and alpha = 2^(-1/(2p)) / (2 sig) is the same family"
    (doseq [[p sig] [[1.5 0.3] [0.7 0.2] [3.0 0.4]]
            normalize? [true false]]
      (is (< (max-err (win/general-gaussian 17 {:p p :sig sig :normalize? normalize?})
                      (win/parzen-exponential 17 {:r (* 2.0 p) :alpha (/ (m/pow 2.0 (/ -1.0 (* 2.0 p))) (* 2.0 sig))
                                                  :normalize? normalize?}))
             1.0e-13)
          (str p " " sig " " normalize?))))
  (testing "the integral of the raw continuous function is 1"
    (doseq [[p sig] [[1.5 0.25] [4.0 0.3] [1.0 2.0]]]
      (is (m/delta-eq 1.0 (simpson (win/->general-gaussian p sig) -0.5 0.5 4000) 1.0e-10) (str p " " sig))))
  (testing "a sig so small that no sample is non-zero throws, a tiny but representable one does not"
    (is (= :sig (:property (ex-data-of #(win/general-gaussian 9 {:p 1.0 :sig 1.0e-320})))))
    (is (every? #(m/valid-double? %) (win/general-gaussian 9 {:p 1.0 :sig 1.0e-200})))))

;; ---------------------------------------------------------------------------------------------
;; Flat-top sets: peak sidelobe level, equivalent noise bandwidth and end values against
;; flattop_psl_reference.edn (generate_flattop_psl.py: exact transfer function in mpmath, float64 FFT)

(def ^:private fft (t/transformer :real :fft))

(defn- periodic-window
  "The periodic window of `n` samples as a double array: the symmetric window of `n + 1` samples without the last."
  ^doubles [w n]
  (let [^doubles a (m/seq->double-array (w (inc n) {:normalize? false}))]
    (java.util.Arrays/copyOf a (int n))))

(defn- dft-psl
  "Peak sidelobe level in dB of the periodic window `w` of `n` samples: float64 FFT with zero padding to `pad`, the
  maximum of the magnitude over `[f0, f-max]` bins relative to the magnitude at 0."
  [w n pad f0 f-max]
  (let [^doubles a (periodic-window w n)
        ^doubles x (java.util.Arrays/copyOf a (int pad))
        ^doubles out (t/forward-1d fft x)
        mag (fn ^double [^long k]
              (if (zero? k)
                (m/abs (aget out 0))
                (Math/hypot (aget out (* 2 k)) (aget out (inc (* 2 k))))))
        scale (/ (double pad) (double n))
        k0 (long (m/ceil (* f0 (/ (double pad) (double n)))))
        k1 (long (m/floor (* f-max scale)))
        mx (loop [k k0
                  mx 0.0]
             (if (<= k k1) (recur (inc k) (m/max mx (mag k))) mx))]
    (* 20.0 (m/log10 (/ mx (mag 0))))))

(deftest flat-top-sets-peak-sidelobe-level
  ;; Periodic window n = 8192, FFT 2^18, sidelobes from the first zero to 100 bins.
  ;; Observed: |fastmath - mpmath transfer function| <= 0.025 dB (hft95: float64 rounding of the samples limits
  ;; the DFT, not the coefficients; hft196d 0.013, hft248d 0.001), |fastmath - Heinzel Table 2| <= 0.06 dB
  ;; (the table has one decimal). Double precision samples resolve even the -248 dB sidelobes of hft248d.
  (let [{:keys [n-dft pad f-max cases]} @psl-reference]
    (is (= 12 (count cases)))
    (is (= (set flat-top-sets) (set (map :window cases))))
    (doseq [{:keys [window first-zero psl-db heinzel]} cases
            :let [w (additions window)
                  psl (dft-psl w n-dft pad first-zero f-max)]]
      (is (<= (m/abs (- psl psl-db)) 0.03) (str window " dft " psl " transfer function " psl-db))
      (when heinzel
        (is (<= (m/abs (- psl (:psl-db heinzel))) 0.06) (str window " dft " psl " Heinzel " (:psl-db heinzel)))))))

(deftest flat-top-sets-coefficients-and-bandwidth
  (doseq [{:keys [window coeffs nenbw end-value heinzel]} (:cases @psl-reference)
          :let [w (additions window)]]
    (testing (str window)
      (testing "the window is the cosine sum of the coefficients of the reference (not clipped)"
        (doseq [N [9 16]
                normalize? [true false]]
          (is (< (max-err (w N {:normalize? normalize?})
                          (win/blackman-harris-family N {:coeffs coeffs :normalize? normalize?}))
                 1.0e-14)
              (str N " " normalize?))))
      (testing "normalized equivalent noise bandwidth N sum(w^2) / sum(w)^2 of the periodic window"
        (let [^doubles a (periodic-window w 64)
              s (reduce + (seq a))
              s2 (reduce + (map m/sq (seq a)))]
          (is (m/delta-eq nenbw (/ (* 64.0 s2) (* s s)) 1.0e-12))
          (when heinzel
            ;; Heinzel prints 4 decimals of coefficients with 6-12 digits: observed |difference| <= 1.1e-4
            (is (<= (m/abs (- nenbw (:nenbw heinzel))) 2.0e-4)))))
      (testing "end values: zero for the D sets (continuous, differentiable), a small jump for the others"
        (let [ends (let [r (vec (w 9))] [(first r) (peek r)])]
          (is (every? #(<= (m/abs (- % end-value)) 1.0e-15) ends) (str ends " " end-value))
          (when (or (#{:flat-top-srs} window) (str/ends-with? (name window) "d"))
            (is (every? #(<= (m/abs %) 1.0e-15) ends) "end value 0 +- 1e-15"))))
      (testing "negative lobes are kept"
        (is (neg? (apply min (w 65))))))))

(deftest hft169d-coefficient-correction
  ;; The printed coefficient 0.33673420 (Heinzel eq. 182) gives an end value 3.3e-9 (relative) and PSL -167.9 dB,
  ;; with 0.33673418 the end value is 0 and the PSL -169.52 dB (report: -169.5 dB at 10.41 bins).
  (let [c (first (filter #(= :hft169d (:window %)) (:cases @psl-reference)))
        printed (:printed c)]
    (is (= 0.33673418 (nth (:coeffs c) 4)))
    (is (= 0.3367342 (nth (:coeffs printed) 4)))
    (is (> (:end-value printed) 1.0e-9))
    (is (< (:psl-db c) (- (:psl-db printed) 1.5)))
    (is (> (m/abs (first (win/blackman-harris-family 9 {:coeffs (:coeffs printed)}))) 1.0e-9)
        "the printed set is not continuous at the ends")
    (is (< (m/abs (first (win/hft169d 9))) 1.0e-15))))

;; ---------------------------------------------------------------------------------------------
;; WP 6.3: achieved peak sidelobe levels of the :level windows (documented in the docstrings).
;; PSL = maximum of the magnitude of the DFT (zero padding) from the first local minimum on, relative to DC.

(defn- level-psl
  "Peak sidelobe level in dB of the window values `w`: FFT with zero padding to 2^17, from the first local minimum of
  the magnitude to 100 bins."
  ^double [w]
  (let [n (count w)
        pad (bit-shift-left 1 17)
        ^doubles x (java.util.Arrays/copyOf ^doubles (m/seq->double-array w) (int pad))
        ^doubles out (t/forward-1d fft x)
        mag (fn ^double [^long k]
              (if (zero? k) (m/abs (aget out 0)) (Math/hypot (aget out (* 2 k)) (aget out (inc (* 2 k))))))
        k0 (loop [k 1] (if (<= (mag (inc k)) (mag k)) (recur (inc k)) k))
        k1 (min (dec (quot pad 2)) (* 100 (quot pad n)))
        mx (loop [k k0 mx 0.0] (if (<= k k1) (recur (inc k) (m/max mx (mag k))) mx))]
    (* 20.0 (m/log10 (/ mx (mag 0))))))

(deftest level-windows-achieved-sidelobe-level
  ;; N = 513 (measured also at 2049: the levels differ by 0.1 dB at most). The requested :level is Doerry's design
  ;; parameter; only dolph-chebyshev achieves it. Measured (dB, N = 513):
  ;;   level      dolph   saramaki  ultraspherical(a 2)  legendre   taylor(n 4)
  ;;   -30        -30.00  -30.96    -33.81               -33.61     -30.47
  ;;   -50        -50.00  -47.16    -47.15               -53.17     -44.66
  ;;   -80        -80.00  -73.50    -70.35               -83.57     -55.95
  (doseq [[level dolph sara ultra legendre taylor] [[-30 -30.00 -30.96 -33.81 -33.61 -30.47]
                                                    [-50 -50.00 -47.16 -47.15 -53.17 -44.66]
                                                    [-80 -80.00 -73.50 -70.35 -83.57 -55.95]]
          :let [N 513]]
    (is (<= (m/abs (- dolph (level-psl (win/dolph-chebyshev N {:level level})))) 0.05) (str "dolph " level))
    (is (<= (m/abs (- sara (level-psl (win/saramaki N {:level level})))) 0.05) (str "saramaki " level))
    (is (<= (m/abs (- ultra (level-psl (win/ultraspherical N {:level level})))) 0.05) (str "ultraspherical " level))
    (is (<= (m/abs (- legendre (level-psl (win/legendre N {:level level})))) 0.05) (str "legendre " level))
    (is (<= (m/abs (- taylor (level-psl (win/taylor N {:level level})))) 0.05) (str "taylor " level))))

;; ---------------------------------------------------------------------------------------------
;; WP 6.1, 6.2: documentation of every public var

(def ^:private window-names
  "Public vars which are windows: functions of `[N opts]`."
  (->> (ns-publics 'fastmath.kernel.window)
       (filter (fn [[_ v]] (some #{'[N opts]} (:arglists (meta v)))))
       (map key)
       set))

(defn- variants
  "Constructors `->name` and continuous functions `name-continuous`."
  [s]
  (let [n (name s)]
    (or (str/starts-with? n "->") (str/ends-with? n "-continuous"))))

(deftest public-vars-have-docstrings
  (let [publics (ns-publics 'fastmath.kernel.window)]
    (is (< 150 (count publics)) "the namespace is the one under test")
    (testing "every public var has a non-empty docstring"
      (doseq [[s v] publics]
        (is (and (string? (:doc (meta v))) (not (str/blank? (:doc (meta v))))) (str s))))
    (testing "the namespace has a docstring which names the main functions"
      (let [doc (:doc (meta (the-ns 'fastmath.kernel.window)))]
        (is (string? doc))
        (is (< 1000 (count doc)))
        (doseq [s ['sample-window 'normalize-coefficients 'continuous->window-fn 'discrete->window-fn 'hann 'kaiser
                   'hft248d 'dolph-chebyshev]]
          (is (str/includes? doc (str s)) (str "namespace docstring mentions " s)))))
    (testing "registered windows (functions of N and options): Parameters, Returns and See also"
      (is (<= 76 (count window-names)))
      (doseq [s (sort window-names)
              :let [doc (:doc (meta (publics s)))]]
        (is (str/includes? doc "Parameters") (str s " Parameters"))
        (is (str/includes? doc "Returns") (str s " Returns"))
        (is (str/includes? doc "See also") (str s " See also"))
        (is (str/includes? doc "normalize?") (str s " :normalize?"))))
    (testing "the other functions with arguments (not constructors and continuous functions) state what they return"
      (doseq [[s v] (sort-by key publics)
              :when (and (not (window-names s)) (not (variants s))
                         (some seq (:arglists (meta v)))
                         (not (:macro (meta v))))
              :let [doc (:doc (meta v))]]
        (is (str/includes? doc "Returns") (str s " Returns"))))
    (testing "constructors and continuous functions: a one line docstring is enough, but not empty"
      (doseq [[s v] publics
              :when (variants s)]
        (is (<= 20 (count (:doc (meta v)))) (str s))))
    (testing "docstring style: no LaTeX, no double quotes, no backslashes, no headings, no code blocks"
      (doseq [[s v] publics
              :let [doc (:doc (meta v))]]
        (is (not (str/includes? doc "$")) (str s " LaTeX"))
        (is (not (str/includes? doc "\"")) (str s " double quote"))
        (is (not (str/includes? doc "\\")) (str s " backslash"))
        (is (not (str/includes? doc "```")) (str s " code block"))
        (is (not (re-find #"(?m)^\s*#+ " doc)) (str s " heading")))
      (let [doc (:doc (meta (the-ns 'fastmath.kernel.window)))]
        (is (not (str/includes? doc "\"")))
        (is (not (str/includes? doc "\\")))
        (is (not (str/includes? doc "```")))))))
