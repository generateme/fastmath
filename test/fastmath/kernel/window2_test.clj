(ns fastmath.kernel.window2-test
  "Property tests and identities of `fastmath.kernel.window`: the raw scale (integral 1 of the continuous
  functions, mean 1 of the discrete windows), Princen-Bradley, identities between windows."
  (:require [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel.window :as win]))

(defn- max-err
  "Maximum of `abs(a - b) / max(1, abs(b))`, `##Inf` for different lengths."
  ^double [a b]
  (let [a (vec a)
        b (vec b)]
    (if (= (count a) (count b))
      (reduce m/max 0.0 (map (fn [x y] (m// (m/abs (m/- (double x) (double y))) (m/max 1.0 (m/abs (double y))))) a b))
      ##Inf)))

;; ---------------------------------------------------------------------------------------------
;; WP 4.1, 7.3: the integral of every continuous model window is 1.
;; Independent integrator: tanh-sinh (double exponential) on [0, 1/2] split at the kinks (the functions are even).

(defn- tanh-sinh
  ^double [f ^double a ^double b]
  (let [h 0.0625
        hw (* 0.5 (- b a))
        mid (* 0.5 (+ a b))]
    (* hw h
       (double
        (reduce + 0.0
                (for [k (range -64 65)
                      :let [t (* k h)
                            u (* 0.5 Math/PI (Math/sinh t))
                            x (Math/tanh u)
                            c (Math/cosh u)]
                      :when (< (Math/abs x) 1.0)]
                  (* 0.5 Math/PI (Math/cosh t) (/ 1.0 (* c c)) (double (f (+ mid (* hw x)))))))))))

(defn- integral
  "Integral of the even function `f` over `[-0.5, 0.5]`; `breaks` are the points in `(0, 0.5)` where `f` is not smooth."
  ^double [f breaks]
  (let [pts (concat [0.0] (sort (filter #(< 0.0 % 0.5) breaks)) [0.5])]
    (* 2.0 (double (reduce + 0.0 (map (fn [[a b]] (tanh-sinh f a b)) (partition 2 1 pts)))))))

(def ^:private flat-top-coeffs [0.21557895 0.41663158 0.277263158 0.083578947 0.006947368])

(defn- mks-continuous
  [N]
  (let [a0 (- 0.5363 (/ 0.14 (dec N)))]
    (win/->blackman-harris [a0 (- 0.996 a0) 0.0 0.004] false)))

(defn- continuous-windows
  "Vector of `[label continuous-function breaks]`."
  []
  (vec (concat
        [["rectangular05" win/rectangular05-continuous []]
         ["rectangular" win/rectangular-continuous []]
         ["triangular" win/triangular-continuous []]
         ["parzen" win/parzen-continuous [0.25]]
         ["welch" win/welch-continuous []]
         ["sinc" win/sinc-continuous []]
         ["fejer" win/fejer-continuous []]
         ["de-la-vallee-poussin" win/de-la-vallee-poussin-continuous []]
         ["hamming" win/hamming-continuous []]
         ["hamming-exact" win/hamming-exact-continuous []]
         ["hann" win/hann-continuous []]
         ["bohman" win/bohman-continuous []]
         ["bartlett-hann" win/bartlett-hann-continuous []]
         ["blackman" win/blackman-continuous []]
         ["blackman-exact" win/blackman-exact-continuous []]
         ["blackman-harris" win/blackman-harris-continuous []]
         ["blackman-harris-61db" win/blackman-harris-61db-continuous []]
         ["blackman-harris-67db" win/blackman-harris-67db-continuous []]
         ["blackman-harris-74db" win/blackman-harris-74db-continuous []]
         ["blackman-harris-92db" win/blackman-harris-92db-continuous []]
         ["nutall-3-1st" win/nutall-3-1st-continuous []]
         ["nutall-3-3rd" win/nutall-3-3rd-continuous []]
         ["blackman-nutall" win/blackman-nutall-continuous []]
         ["nutall-1st" win/nutall-1st-continuous []]
         ["nutall-3rd" win/nutall-3rd-continuous []]
         ["nutall-5th" win/nutall-5th-continuous []]
         ["flat-top" win/flat-top-continuous []]
         ["flat-top-3" win/flat-top-3-continuous []]
         ["vorbis" win/vorbis-continuous []]
         ["family flat-top unclipped" (win/->blackman-harris flat-top-coeffs false) []]
         ["family hann" (win/->blackman-harris [0.5 0.5] false) []]
         ["mks N=8" (mks-continuous 8) []]
         ["mks N=64" (mks-continuous 64) []]]
        (for [order [1 2 3 4 5 8 20 60]]
          [(str "b-spline " order) (win/->b-spline order) (map #(- 0.5 (/ % order)) (range (inc order)))])
        (for [a [0.5 1.0 2.0]] [(str "connes " a) (win/->connes a) [(/ a 2.0)]])
        (for [[g u] [[1.0 3.0] [0.5 2.0] [0.3 5.0] [1.0 1.0]]]
          [(str "parzen-algebraic " g " " u) (win/->parzen-algebraic g u) []])
        (for [o [0 1 2 5 20 100]] [(str "singla-singh " o) (win/->singla-singh o) []])
        (for [L [0.0 0.5 1.0 3.0 6.0]] [(str "lanczos " L) (win/->lanczos L) []])
        (for [a [0.5 0.75 1.0]] [(str "raised-cosine " a) (win/->raised-cosine a) []])
        (for [v [-0.4 0.0 1.0 5.0 50.0]] [(str "webster-hamming " v) (win/->webster-hamming v) []])
        (for [mm [0.0 0.5 1.0 2.0 5.0]] [(str "power-of-cosine " mm) (win/->power-of-cosine mm) []])
        (for [[mm a] [[1.0 0.05] [2.0 0.5] [1.0 0.0] [1.0 1.0]]]
          [(str "raised-power-of-cosine " mm " " a) (win/->raised-power-of-cosine mm a) []])
        (for [[mm g] [[2.0 1.0] [2.0 0.5] [1.0 1.0] [3.0 0.7]]]
          [(str "parzen-cosine " mm " " g) (win/->parzen-cosine mm g) []])
        (for [a [0.0 0.1 0.25 0.5]] [(str "trapezoid " a) (win/->trapezoid a) [a]])
        (for [a [0.0 0.1 0.5 1.0 1.5]] [(str "tukey " a) (win/->tukey a) [(- 0.5 (* 0.5 (min a 1.0)))]])
        (for [a [0.1 1.0 5.0]] [(str "exponential " a) (win/->exponential a) []])
        (for [a [0.1 1.0 5.0]] [(str "hanning-poisson " a) (win/->hanning-poisson a) []])
        (for [a [0.1 1.0 5.0]] [(str "gaussian " a) (win/->gaussian a) []])
        (for [[a r] [[1.0 1.0] [3.0 0.5] [1.0 0.25] [2.0 2.0] [1.0 5.0]]]
          [(str "parzen-exponential " a " " r) (win/->parzen-exponential a r) []])
        (for [a [0.1 1.0 5.0]] [(str "cauchy " a) (win/->cauchy a) []])
        (for [[a r] [[1.0 1.0] [1.0 0.5] [2.0 0.8] [1.0 3.0] [0.5 0.25]]]
          [(str "parzen-geometric " a " " r) (win/->parzen-geometric a r) []])
        (for [a [0.1 1.0 3.0 10.0]] [(str "kaiser-bessel " a) (win/->kaiser-bessel a) []])
        (for [a [0.1 2.0 5.0]] [(str "cosh " a) (win/->cosh a) []])
        (for [a [0.1 0.5 2.0 5.0]] [(str "avci-nacaroglu " a) (win/->avci-nacaroglu a) []])
        (for [a [0.1 0.5 2.0 5.0]] [(str "knab " a) (win/->knab a) []])
        (for [a [0.1 2.0 5.0]] [(str "bessel-I1 " a) (win/->bessel-I1 a) []])
        (for [[level n] [[-50.0 4] [-30.0 2] [-80.0 8] [-35.0 3]]]
          [(str "taylor " level " " n) (win/->taylor level n) []]))))

(def ^:private integral-tolerances
  "Per label tolerance of the raw integral, otherwise 1e-12. Observed maxima are in comments."
  {"webster-hamming -0.4" 1.0e-9})   ; the end singularity cos^-0.4: the tanh-sinh nodes round to x = 1/2

(deftest continuous-windows-integral-one
  ;; Observed maximum |integral - 1| over the 133 rows: 1.1e-15 (1.2e-11 for webster-hamming -0.4, the error of the
  ;; integrator at the cos^-0.4 end singularity). The old code gave 1/order for b-spline, (4 - 2 alpha) / alpha for
  ;; tukey and 2.754 for webster-hamming -0.4.
  (doseq [[label f breaks] (continuous-windows)]
    (let [z (integral f breaks)]
      (is (< (m/abs (- z 1.0)) (get integral-tolerances label 1.0e-12)) (str label " integral " z)))))

(deftest raw-windows-sample-the-continuous-function
  (testing "the raw window is the continuous function sampled at -0.5 ... 0.5 (default alpha / order)"
    (doseq [[w f] [[win/b-spline (win/->b-spline 3)] [win/tukey (win/->tukey 0.5)]
                   [win/webster-hamming (win/->webster-hamming 1.0)] [win/knab (win/->knab 2.0)]
                   [win/avci-nacaroglu (win/->avci-nacaroglu 2.0)] [win/singla-singh (win/->singla-singh 1)]
                   [win/parzen-exponential (win/->parzen-exponential 1.0 1.0)]]]
      (is (< (max-err (w 9 {:normalize? false}) (map f (m/sample identity -0.5 0.5 9))) 1.0e-14)))))

;; ---------------------------------------------------------------------------------------------
;; WP 4.1: discrete windows have mean 1 (Doerry).

(deftest discrete-windows-mean-one
  (doseq [[label w opts] (concat
                          (for [l [-30.0 -50.0 -100.0]] ["dolph-chebyshev" win/dolph-chebyshev {:level l}])
                          (for [a [-0.5 0.5 2.0 3.0] l [-30.0 -70.0]] ["ultraspherical" win/ultraspherical {:level l :alpha a}])
                          (for [l [-30.0 -80.0]] ["saramaki" win/saramaki {:level l}])
                          (for [l [-30.0 -50.0 -100.0]] ["legendre" win/legendre {:level l}])
                          [["shayesteh-kashtiban" win/shayesteh-kashtiban {}]])
          N [1 2 3 8 9 64 65]]
    (let [raw (vec (w N (assoc opts :normalize? false)))
          mean (/ (reduce + raw) N)]
      (is (= N (count raw)))
      ;; observed |mean - 1| <= 4e-15
      (is (< (m/abs (- mean 1.0)) 1.0e-12) (str label " " opts " N=" N " mean " mean)))))

(deftest shayesteh-kashtiban-raw-scale
  ;; WP 4.1: the old code had sum 1 (mean 1/N); Doerry eq. 277 has mean 1.
  (is (m/delta-eq 64.0 (reduce + (win/shayesteh-kashtiban 64 {:normalize? false})) 1.0e-10)))

(deftest kaiser-bessel-derived-raw-scale
  ;; The raw window is the power complementary one, not the unit mean scale (documented).
  (let [w (vec (win/kaiser-bessel-derived 64))]
    (is (= w (vec (win/kaiser-bessel-derived 64 {:normalize? false}))) "default is raw")
    (is (< (/ (reduce + w) 64.0) 0.9))
    (is (= 1.0 (apply max (win/kaiser-bessel-derived 64 {:normalize? true}))))))

;; ---------------------------------------------------------------------------------------------
;; WP 4.6: Princen-Bradley, w[n]^2 + w[n + N/2]^2 = 1.

(defn- princen-bradley
  ^double [w]
  (let [w (vec w)
        h (quot (count w) 2)]
    (reduce m/max 0.0 (map (fn [i] (m/abs (- (+ (m/sq (w i)) (m/sq (w (+ i h)))) 1.0))) (range h)))))

(deftest vorbis-princen-bradley
  ;; Observed maximum 2.2e-16 over even N 2..512 (the old fastmath vorbis: .15 at N 4, .0097 at N 256).
  (doseq [N (range 2 513 2)]
    (is (<= (princen-bradley (win/vorbis N)) 1.0e-15) (str N))))

(deftest kaiser-bessel-derived-princen-bradley
  ;; Observed maximum 4.4e-16.
  (doseq [alpha [0.5 1.0 2.0 4.0]
          N [2 4 8 16 64 256]]
    (is (<= (princen-bradley (win/kaiser-bessel-derived N {:alpha alpha})) 1.0e-15) (str alpha " " N))))

;; ---------------------------------------------------------------------------------------------
;; WP 7.3, W2-20: identities between windows (default and raw scale, N = 9 and 16).

(defn- same?
  "Windows `a` and `b` (functions of `[N opts]`) agree for N in 9 16 in both scales within `tol`."
  [a opts-a b opts-b tol]
  (every? (fn [[N normalize?]]
            (< (max-err (a N (assoc opts-a :normalize? normalize?)) (b N (assoc opts-b :normalize? normalize?))) tol))
          (for [N [9 16] n [true false]] [N n])))

(deftest identities-continuous
  ;; Observed maxima: b-spline 4 / parzen 1.2e-15, lanczos L 4 / de-la-vallee-poussin 7.7e-16, parzen-exponential r 2 /
  ;; gaussian 3.6e-16; all others below 1e-15 (the quadrature constants of the old code were only 1e-9 accurate).
  (testing "parzen-exponential r 1 = exponential, r 2 = gaussian with alpha sqrt(2)"
    (doseq [a [0.5 1.0 3.0]]
      (is (same? win/parzen-exponential {:alpha a :r 1.0} win/exponential {:alpha a} 1.0e-14) (str a))
      (is (same? win/parzen-exponential {:alpha a :r 2.0} win/gaussian {:alpha (* a m/SQRT2)} 1.0e-14) (str a))))
  (testing "parzen-geometric r 2 = cauchy"
    (doseq [a [0.5 1.0 3.0]]
      (is (same? win/parzen-geometric {:alpha a :r 2.0} win/cauchy {:alpha a} 1.0e-14) (str a))))
  (testing "b-spline 1 / 2 / 4 = rectangular / triangular shift 0 / continuous parzen"
    (is (same? win/b-spline {:order 1} win/rectangular {} 1.0e-14))
    (is (same? win/b-spline {:order 2} win/triangular {:shift 0} 1.0e-14))
    (is (same? win/b-spline {:order 4} win/parzen {:discrete? false} 1.0e-14)))
  (testing "singla-singh 0 = triangular shift 0, trapezoid 0 = triangular, trapezoid 0.5 = rectangular"
    (is (same? win/singla-singh {:order 0} win/triangular {:shift 0} 1.0e-14))
    (is (same? win/trapezoid {:alpha 0.0} win/triangular {:shift 0} 1.0e-14))
    (is (same? win/trapezoid {:alpha 0.5} win/rectangular {} 1.0e-14)))
  (testing "power-of-cosine m 2 = hann, raised-power-of-cosine alpha 0 = power-of-cosine"
    (is (same? win/power-of-cosine {:m 2.0} win/hann {} 1.0e-14))
    (is (same? win/raised-power-of-cosine {:alpha 0.0 :m 2.5} win/power-of-cosine {:m 2.5} 1.0e-14)))
  (testing "lanczos L 1 / 2 / 4 = sinc / fejer / de-la-vallee-poussin"
    (is (same? win/lanczos {:L 1.0} win/sinc {} 1.0e-14))
    (is (same? win/lanczos {:L 2.0} win/fejer {} 1.0e-14))
    (is (same? win/lanczos {:L 4.0} win/de-la-vallee-poussin {} 1.0e-14)))
  (testing "tukey 0 = rectangular, tukey 1 and larger = hann"
    (is (same? win/tukey {:alpha 0.0} win/rectangular {} 1.0e-14))
    (is (same? win/tukey {:alpha 1.0} win/hann {} 1.0e-14))
    (is (same? win/tukey {:alpha 1.7} win/hann {} 1.0e-14)))
  (testing "raised-cosine = hann, hamming, hamming-exact; webster-hamming v 0 = hamming-exact"
    (is (same? win/raised-cosine {:alpha 0.5} win/hann {} 1.0e-14))
    (is (same? win/raised-cosine {:alpha 0.54} win/hamming {} 1.0e-14))
    (is (same? win/raised-cosine {:alpha (/ 25.0 46.0)} win/hamming-exact {} 1.0e-14))
    (is (same? win/webster-hamming {:v 0.0} win/hamming-exact {} 1.0e-14)))
  (testing "blackman-harris-family with the coefficients of a named window"
    (is (same? win/blackman-harris-family {:coeffs [0.5 0.5]} win/hann {} 1.0e-14))
    (is (same? win/blackman-harris-family {:coeffs [0.42 0.5 0.08]} win/blackman {} 1.0e-14))
    (is (same? win/blackman-harris-family {:coeffs [0.35875 0.48829 0.14128 0.01168]} win/blackman-harris-92db {} 1.0e-14))
    (is (same? win/blackman-harris-family {:coeffs [0.4243801 0.4973406 0.0782793]} win/blackman-harris {} 1.0e-14))))

(deftest identities-discrete
  ;; Observed 0 (saramaki / ultraspherical) and below 1e-15 (legendre).
  (testing "saramaki = ultraspherical alpha 1, ultraspherical alpha 0 = dolph-chebyshev"
    (doseq [l [-30.0 -50.0 -80.0]]
      (is (same? win/saramaki {:level l} win/ultraspherical {:level l :alpha 1.0} 1.0e-14) (str l))
      (is (same? win/ultraspherical {:level l :alpha 0.0} win/dolph-chebyshev {:level l} 1.0e-14) (str l))))
  (testing "legendre L = ultraspherical alpha 1/2 with the level 1.0754 abs(L) + 1.7388 (Doerry eq. 266)"
    (doseq [l [-30.0 -50.0 -80.0]]
      (is (same? win/legendre {:level l} win/ultraspherical {:level (- (+ (* 1.0754 (m/abs l)) 1.7388)) :alpha 0.5} 1.0e-14)
          (str l)))))

(deftest family-flat-top-is-not-clipped
  ;; WP 4.7: raw integral 1 with negative lobes, the family with the flat-top coefficients is flat-top.
  (is (same? win/blackman-harris-family {:coeffs flat-top-coeffs} win/flat-top {} 1.0e-14))
  (is (neg? (apply min (win/blackman-harris-family 65 {:coeffs flat-top-coeffs :normalize? false})))))
