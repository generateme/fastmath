(ns fastmath.kernel.density-test
  "Tests of `fastmath.kernel.density`: kernels, `kde-data`, exact and binned evaluation, validation, confidence intervals.
  The oracle is `test/resources/kernel/kde_reference.edn` (utils/fastmath/dev/generate_kde_reference.py): exact direct sums
  in mpmath (40 digits, no cut-off) and the constants of the kernels by mpmath quadrature."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.kernel.density :as d]
            [fastmath.stats :as stats]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "kernel/kde_reference.edn")))))

(def ^:private kernels
  [:uniform :triangular :epanechnikov :quartic :triweight :tricube :gaussian :cosine
   :logistic :sigmoid :silverman :wigner :cauchy :laplace])

(def ^:private compact-kernels
  #{:uniform :triangular :epanechnikov :quartic :triweight :tricube :cosine :wigner})

(defn- close?
  "Relative and absolute error `abs(a - b) <= tol * max(1, abs(b))`, `##Inf` and `##NaN` must be equal."
  [a b tol]
  (let [a (double a)
        b (double b)]
    (cond
      (m/nan? b) (m/nan? a)
      (m/inf? b) (== a b)
      :else (<= (m/abs (- a b)) (* tol (max 1.0 (m/abs b)))))))

(defn- max-rel-error
  "Maximum of `abs(f(x) - y)` over the points divided by the maximum of `abs(y)`."
  [f xs ys]
  (/ (reduce max 0.0 (map (fn [x y] (m/abs (- (double (f x)) (double y)))) xs ys))
     (reduce max (map #(m/abs (double %)) ys))))

(defn- invalid-input?
  "True when `thunk` throws `ex-info` with `:type :invalid-input` and the `property` (and `function`) given."
  [property thunk]
  (try (thunk) false
       (catch clojure.lang.ExceptionInfo e
         (let [data (ex-data e)]
           (and (= :invalid-input (:type data))
                (= property (:property data))
                (contains? data :function)
                (contains? data :value)
                (contains? data :expected))))))

;; ---------------------------------------------------------------------------------------------
;; kernels and kde-data

(deftest kernel-functions-test
  (let [{:keys [x values]} (:kernel-values @reference)]
    (doseq [kernel kernels
            :let [f (get-in d/kde-data [kernel :kernel])]]
      (testing (str "values of the kernel " kernel)
        (is (every? true? (map (fn [xi yi] (close? (f xi) yi 1.0e-15)) x (get values kernel)))
            (str kernel)))
      (testing (str "symmetry of the kernel " kernel)
        (is (every? #(= (f %) (f (- %))) [0.1 0.5 0.9999 1.0 1.5 3.0 10.0 100.0 1.0e300]))))))

(deftest kernel-edges-test
  (testing "compact kernels: closed support [-1, 1], zero outside, values just inside"
    (doseq [kernel compact-kernels
            :let [f (get-in d/kde-data [kernel :kernel])]]
      (is (zero? (f (Math/nextUp 1.0))) (str kernel))
      (is (zero? (f (- (Math/nextUp 1.0)))) (str kernel))
      (is (>= (f (- 1.0 1.0e-9)) 0.0) (str kernel))))
  (testing "silverman at infinity is 0, not NaN (KE-18)"
    (is (zero? (d/silverman ##Inf)))
    (is (zero? (d/silverman ##-Inf))))
  (testing "no NaN of any kernel at large arguments"
    (doseq [kernel kernels
            :let [f (get-in d/kde-data [kernel :kernel])]
            x [1.0e3 1.0e100 ##Inf ##-Inf]]
      (is (not (m/nan? (f x))) (str kernel " " x))))
  (testing "(1 - |x|)(1 + |x|) form: no cancellation next to the edge of the support"
    (is (close? (d/epanechnikov (- 1.0 1.0e-12)) (* 0.75 2.0e-12) 1.0e-9))
    (is (close? (d/quartic (- 1.0 1.0e-12)) (* 0.9375 2.0e-12 2.0e-12) 1.0e-9))
    (is (close? (d/triweight (- 1.0 1.0e-12)) (* 1.09375 (Math/pow 2.0e-12 3)) 1.0e-9))))

(deftest kde-data-constants-test
  (let [consts (:kernel-constants @reference)]
    (is (= (set kernels) (set (keys d/kde-data))))
    (doseq [kernel kernels
            :let [data (d/kde-data kernel)
                  ref (consts kernel)]]
      (testing (str "constants of " kernel)
        (is (= (:kernel data) (get-in d/kde-data [kernel :kernel])))
        (is (close? (:k2 data) (:k2 ref) 1.0e-14))
        (is (close? (:x2k data) (:x2k ref) 1.0e-14))
        (is (close? (:delta0 data) (:delta0 ref) 1.0e-14))
        (is (close? (:nrd-scale data) (:nrd-scale ref) 1.0e-14))
        (is (close? (:efficiency data) (:efficiency ref) 1.0e-14))))
    (testing "gaussian is the unit of the canonical bandwidth"
      (is (close? (:nrd-scale (d/kde-data :gaussian)) 1.0 1.0e-15)))
    (testing "undefined canonical bandwidth: silverman (x2k = 0), cauchy (x2k = Inf)"
      (is (= ##Inf (:nrd-scale (d/kde-data :silverman))))
      (is (zero? (:nrd-scale (d/kde-data :cauchy)))))))

(deftest kde-data-radius-test
  (let [consts (:kernel-constants @reference)]
    (testing "compact kernels: radius = extent = 1"
      (doseq [kernel compact-kernels]
        (is (= 1.0 (:radius (d/kde-data kernel)) (:extent (d/kde-data kernel))) (str kernel))))
    (testing "radius is where |K| < 1e-10: not below the mpmath crossing, at most 0.15 above"
      (doseq [kernel [:gaussian :logistic :sigmoid :silverman :laplace]
              :let [crossing (:radius (consts kernel))
                    radius (:radius (d/kde-data kernel))]]
        (is (<= crossing radius (+ crossing 0.15)) (str kernel " " crossing " " radius))
        (is (< (m/abs ((get-in d/kde-data [kernel :kernel]) radius)) 1.0e-10) (str kernel))
        (is (= radius (:extent (d/kde-data kernel))))))
    (testing "cauchy: every point is summed (radius Inf), the extent is 1000 (KE-07, KE-15)"
      (is (= ##Inf (:radius (d/kde-data :cauchy))))
      (is (= 1000.0 (:extent (d/kde-data :cauchy)))))))

;; ---------------------------------------------------------------------------------------------
;; exact evaluation

(deftest kernel-density-anchors-test
  (testing "KE-01: the density in a gap of the data is not lost (the forward scan returned 0.0 here)"
    (doseq [{:keys [kernel data h x density id]} (:anchors @reference)]
      (is (close? ((d/kernel-density kernel data {:bandwidth h}) x) density 1.0e-12) (str id))))
  (testing "hand computed"
    (is (close? ((d/kernel-density :epanechnikov [0.0 10.0] {:bandwidth 1.0}) 9.5) 0.28125 1.0e-15))
    (is (close? ((d/kernel-density :uniform [0.0 10.0] {:bandwidth 1.0}) 9.5) 0.25 1.0e-15))
    (is (close? ((d/kernel-density :gaussian [0.0 100.0] {:bandwidth 1.0}) 98.0)
                (* 0.5 (/ (Math/exp -2.0) (Math/sqrt (* 2.0 Math/PI)))) 1.0e-14))))

(deftest kernel-density-reference-test
  (testing "every kernel, 13 datasets (gaps, ties, tails, offset 1e6, single point), 240 points each"
    (doseq [{:keys [kernel data h x density dataset]} (:cases @reference)
            :let [f (d/kernel-density kernel data {:bandwidth h})
                  err (max-rel-error f x density)
                  ;; the sum omits the kernel beyond 1.01 * h * radius where |K| < 1e-10: 3e-10 measured for the tails
                  tol (cond (compact-kernels kernel) 1.0e-12
                            (= kernel :cauchy) 1.0e-12
                            :else 5.0e-10)]]
      (is (<= err tol) (str kernel " " dataset " " err)))))

(deftest kernel-density-properties-test
  (let [data [0.0 0.5 0.6 3.0 3.1 9.0]]
    (testing "the density is a probability density: the integral over [mn, mx] is 1 up to the cut-off mass"
      (doseq [kernel kernels
              :let [{:keys [kde mn mx]} (d/kernel-density+ kernel data {:bandwidth 0.7})
                    n 400000
                    dx (/ (- mx mn) n)
                    integral (* dx (+ (* 0.5 (+ (kde mn) (kde mx)))
                                      (reduce + (map #(kde (+ mn (* % dx))) (range 1 n)))))]]
        (is (<= (- 1.0 (if (= kernel :cauchy) 4.0e-4 1.0e-4)) integral (+ 1.0 1.0e-4)) (str kernel " " integral))))
    (testing "mn and mx: the data extended by 1.01 * h * extent"
      (doseq [kernel kernels
              :let [{:keys [mn mx factor h]} (d/kernel-density+ kernel data {:bandwidth 0.7})
                    reach (* 1.01 0.7 (:extent (d/kde-data kernel)))]]
        (is (close? mn (- 0.0 reach) 1.0e-12) (str kernel))
        (is (close? mx (+ 9.0 reach) 1.0e-12) (str kernel))
        (is (close? factor (/ 1.0 (* 6 0.7)) 1.0e-15))
        (is (= 0.7 h))))
    (testing "the density is zero away from the data only for compact kernels"
      (is (zero? ((d/kernel-density :epanechnikov data {:bandwidth 0.7}) 5.0)))
      (is (pos? ((d/kernel-density :gaussian data {:bandwidth 0.7}) 5.0)))
      (is (pos? ((d/kernel-density :cauchy data {:bandwidth 0.7}) 1.0e6))))
    (testing "the order of the data and a repeated call do not matter, the input is not modified"
      (let [a (double-array [3.0 1.0 2.0 2.0 5.0])
            f1 (d/kernel-density :gaussian a {:bandwidth 0.5})
            f2 (d/kernel-density :gaussian [5.0 2.0 2.0 3.0 1.0] {:bandwidth 0.5})]
        (is (= [3.0 1.0 2.0 2.0 5.0] (vec a)))
        (is (= (f1 2.2) (f2 2.2) (f1 2.2)))))
    (testing "NaN in, NaN out; infinity gives 0 (infinite support kernels too)"
      (is (m/nan? ((d/kernel-density :gaussian data {:bandwidth 0.5}) ##NaN)))
      (is (zero? ((d/kernel-density :gaussian data {:bandwidth 0.5}) ##Inf)))
      (is (zero? ((d/kernel-density :epanechnikov data {:bandwidth 0.5}) ##-Inf)))
      (is (zero? ((d/kernel-density :cauchy data {:bandwidth 0.5}) ##Inf))))
    (testing "single point, one data value repeated"
      (is (close? ((d/kernel-density :gaussian [2.0] {:bandwidth 0.5}) 2.0) (/ (d/gaussian 0.0) 0.5) 1.0e-15))
      (is (close? ((d/kernel-density :gaussian [2.0 2.0 2.0] {:bandwidth 0.5}) 2.0) (/ (d/gaussian 0.0) 0.5) 1.0e-15)))))

(deftest custom-kernel-test
  (let [phi (fn [x] (* 0.1 (d/gaussian (/ x 10.0))))
        data [0.0 1.0 2.0 50.0]
        exact (fn [x] (/ (reduce + (map #(phi (/ (- x %) 1.0)) data)) (* 4 1.0)))]
    (testing "KE-08: a custom kernel without options is summed over every point"
      (let [f (d/kernel-density phi data {:bandwidth 1.0})]
        (doseq [x [-30.0 0.0 1.0 25.0 49.0 50.0 80.0]]
          (is (close? (f x) (exact x) 1.0e-12) (str x)))))
    (testing "radius and extent options"
      (let [{:keys [mn mx kde]} (d/kernel-density+ phi data {:bandwidth 1.0 :radius 100.0 :extent 20.0})]
        (is (close? mn (- 0.0 (* 1.01 20.0)) 1.0e-12))
        (is (close? mx (+ 50.0 (* 1.01 20.0)) 1.0e-12))
        (is (close? (kde 25.0) (exact 25.0) 1.0e-12)))
      (let [{:keys [mn mx]} (d/kernel-density+ phi data {:bandwidth 1.0})]
        (is (close? mn (- 0.0 (* 1.01 5.0)) 1.0e-12))
        (is (close? mx (+ 50.0 (* 1.01 5.0)) 1.0e-12))))
    (testing "a kernel function which is not a primitive function"
      (let [boxed (fn [x] (Math/exp (- (* x x))))]
        (is (close? ((d/kernel-density boxed [0.0 1.0] {:bandwidth 1.0}) 0.5)
                    (/ (+ (Math/exp -0.25) (Math/exp -0.25)) 2.0) 1.0e-15))))
    (testing "invalid radius and extent"
      (doseq [v [0 -1 ##NaN "a"]]
        (is (invalid-input? :radius #(d/kernel-density phi data {:bandwidth 1.0 :radius v})) (str v)))
      (doseq [v [0 -1 ##Inf ##NaN]]
        (is (invalid-input? :extent #(d/kernel-density phi data {:bandwidth 1.0 :extent v})) (str v))))
    (testing "the options are ignored for a predefined kernel"
      (is (= ((d/kernel-density :epanechnikov data {:bandwidth 1.0}) 1.0)
             ((d/kernel-density :epanechnikov data {:bandwidth 1.0 :radius 50.0 :extent 50.0}) 1.0))))))

;; ---------------------------------------------------------------------------------------------
;; binned evaluation

(deftest binned-density-test
  (let [cases (filter #(#{:n30-normal :n200-normal :bimodal-wide-gap :heavy-tail-t2 :offset-1e6 :ties200} (:dataset %))
                      (:cases @reference))
        worst (fn [kernel divisor]
                (reduce max
                        (for [{:keys [data h x density]} (filter #(= kernel (:kernel %)) cases)]
                          (max-rel-error (d/kernel-density kernel data {:bandwidth h :binned? divisor}) x density))))]
    (testing "KE-02: no exception at any point, 14 kernels (the old code threw ArrayIndexOutOfBounds)"
      (doseq [{:keys [kernel data h x]} cases
              :let [f (d/kernel-density kernel data {:bandwidth h :binned? true})]]
        (is (every? #(not (m/nan? (f %))) x) (str kernel))))
    (testing "second order convergence for smooth kernels: error 5e-3 at divisor 5, 2e-4 at divisor 40 (measured 3e-3, 4e-5),
              at least 8 times smaller from divisor 5 to 20"
      (doseq [kernel [:gaussian :logistic :sigmoid :silverman]
              :let [e5 (worst kernel 5)
                    e20 (worst kernel 20)
                    e40 (worst kernel 40)]]
        (is (<= e5 5.0e-3) (str kernel " " e5))
        (is (<= e40 2.0e-4) (str kernel " " e40))
        (is (<= (* 8.0 e20) e5) (str kernel " " e5 " " e20))))
    (testing "cauchy: smooth but the peak is narrow (the effective bandwidth is h / 2): 2e-2 at divisor 5, 6e-4 at divisor 40 (measured 1.9e-2, 5.3e-4),
              at least 8 times smaller from divisor 5 to 20"
      (let [e5 (worst :cauchy 5)
            e20 (worst :cauchy 20)
            e40 (worst :cauchy 40)]
        (is (<= e5 2.0e-2) (str e5))
        (is (<= e40 6.0e-4) (str e40))
        (is (<= (* 8.0 e20) e5) (str e5 " " e20))))
    (testing "kernels with kinks (the support edge, the origin) and jumps: the error decreases with the divisor,
              1e-1 at divisor 5 (measured 7e-2 uniform excluded), 2e-2 at divisor 40 (measured 1.2e-2)"
      (doseq [kernel [:triangular :epanechnikov :quartic :triweight :tricube :cosine :laplace]
              :let [e5 (worst kernel 5)
                    e40 (worst kernel 40)]]
        (is (<= e5 1.0e-1) (str kernel " " e5))
        (is (<= e40 2.0e-2) (str kernel " " e40))
        (is (< (* 2.0 e40) e5) (str kernel))))
    (testing "uniform and wigner (jump, square root edge): converge slowly, only bounded here (measured 0.27 and 0.13 at divisor 5)"
      (doseq [kernel [:uniform :wigner]]
        (is (<= (worst kernel 5) 0.4) (str kernel))
        (is (<= (worst kernel 40) (worst kernel 5)) (str kernel))))
    (testing "binned? true is the divisor 5; the factor and the support are those of the exact estimator"
      (let [data [0.0 1.0 2.5 4.0 4.1]
            exact (d/kernel-density+ :gaussian data {:bandwidth 0.8})
            binned (d/kernel-density+ :gaussian data {:bandwidth 0.8 :binned? true})]
        (is (= (:factor exact) (:factor binned)))
        (is (= (:mn exact) (:mn binned)))
        (is (= (:mx exact) (:mx binned)))
        (is (= ((:kde binned) 2.0) ((d/kernel-density :gaussian data {:bandwidth 0.8 :binned? 5}) 2.0)))
        (is (= ((:kde binned) 2.0) ((d/kernel-density :gaussian data {:bandwidth 0.8 :binned? 5.0}) 2.0)))))
    (testing "the mass is preserved by the linear binning: the estimate of a single point is exact on the node"
      (is (close? ((d/kernel-density :gaussian [0.3] {:bandwidth 1.0 :binned? true}) 0.3) (d/gaussian 0.0) 1.0e-12))
      (is (close? ((d/kernel-density :gaussian [1.0 1.0 1.0] {:bandwidth 1.0 :binned? true}) 1.0) (d/gaussian 0.0) 1.0e-12)))
    (testing "a datum on the last node, equal data, a huge divisor"
      (is (close? ((d/kernel-density :gaussian [0.0 1.0] {:bandwidth 1.0 :binned? 5}) 1.0)
                  ((d/kernel-density :gaussian [0.0 1.0] {:bandwidth 1.0}) 1.0) 1.0e-12))
      (is (close? ((d/kernel-density :gaussian [0.0 1.0 2.0] {:bandwidth 1.0 :binned? 1000}) 1.0)
                  ((d/kernel-density :gaussian [0.0 1.0 2.0] {:bandwidth 1.0}) 1.0) 1.0e-3)))))

(deftest binned-invalid-test
  (testing "KE-03: a non positive or NaN divisor throws at once (the old code hung)"
    (doseq [v [-1 0 0.0 -5.5 ##NaN ##Inf ##-Inf "a" :x]]
      (is (invalid-input? :binned? #(d/kernel-density :gaussian [1.0 2.0 3.0] {:bandwidth 1.0 :binned? v})) (str v))))
  (testing "nil and false: no binning"
    (is (= ((d/kernel-density :gaussian [1.0 2.0 3.0] {:bandwidth 1.0 :binned? nil}) 1.5)
           ((d/kernel-density :gaussian [1.0 2.0 3.0] {:bandwidth 1.0 :binned? false}) 1.5)
           ((d/kernel-density :gaussian [1.0 2.0 3.0] {:bandwidth 1.0}) 1.5))))
  (testing "more than 2^24 grid nodes"
    (is (invalid-input? :binned? #((d/kernel-density :gaussian [0.0 1.0e9] {:bandwidth 1.0 :binned? 5}) 0.0)))
    (is (invalid-input? :binned? #((d/kernel-density :gaussian [0.0 1.0e5] {:bandwidth 1.0 :binned? 1000}) 0.0)))
    (is (number? ((d/kernel-density :gaussian [0.0 1.0e5] {:bandwidth 1.0 :binned? 5}) 0.0))))
  (testing "stats/modes with :binned? <= 0 throws instead of hanging"
    (is (thrown? clojure.lang.ExceptionInfo
                 (stats/modes [1.0 2.0 2.0 3.0 7.0] :kde {:binned? -1})))))

;; ---------------------------------------------------------------------------------------------
;; validation

(deftest validation-test
  (let [good [1.0 2.0 4.0 8.0]]
    (testing "data: not empty, finite numbers"
      (doseq [bad [[] nil '() [1.0 ##NaN] [1.0 ##Inf] [##-Inf] ["a" 1.0] [nil 1.0]]]
        (is (invalid-input? :data #(d/kernel-density :gaussian bad {:bandwidth 1.0})) (pr-str bad))
        (is (invalid-input? :data #(d/kernel-density+ :gaussian bad {:bandwidth 1.0})) (pr-str bad))
        (is (invalid-input? :data #(d/kernel-density-ci :gaussian bad {:bandwidth 1.0})) (pr-str bad))))
    (testing "kernel: unknown keyword, anything but a keyword or a function"
      (doseq [bad [:foo nil 1 "gaussian" 'gaussian [:gaussian]]]
        (is (invalid-input? :kernel #(d/kernel-density bad good {:bandwidth 1.0})) (pr-str bad))))
    (testing "bandwidth: unknown keyword and not a positive finite number"
      (is (invalid-input? :bandwidth #(d/kernel-density :gaussian good {:bandwidth :foo})))
      (doseq [bad [0 0.0 -1 ##NaN ##Inf ##-Inf "a" [1.0]]]
        (is (invalid-input? :bandwidth #(d/kernel-density :gaussian good {:bandwidth bad})) (pr-str bad))))
    (testing "ex-data: the invalid value and the message with the valid keywords"
      (let [e (try (d/kernel-density :foo good) (catch clojure.lang.ExceptionInfo e e))]
        (is (= :foo (:value (ex-data e))))
        (is (re-find #":gaussian" (ex-message e))))
      (let [e (try (d/kernel-density :gaussian good {:bandwidth :foo}) (catch clojure.lang.ExceptionInfo e e))]
        (is (re-find #":lscv" (ex-message e)))))
    (testing "nil option = missing option"
      (is (= ((d/kernel-density :gaussian good) 3.0)
             ((d/kernel-density :gaussian good {:bandwidth nil :binned? nil :bandwidth-options nil :radius nil :extent nil}) 3.0)
             ((d/kernel-density :gaussian good nil) 3.0))))
    (testing "the default bandwidth is :nrd"
      (is (= (:h (d/kernel-density+ :gaussian good)) (d/bandwidth :gaussian good :nrd))))
    (testing "bandwidth-options are forwarded to the bandwidth selector"
      (is (= (:h (d/kernel-density+ :gaussian [1.0 2.0 4.0 8.0 3.0 5.0 9.0 2.5]
                                    {:bandwidth :lscv :bandwidth-options {:bounds [0.5 2.0]}}))
             (d/bandwidth :gaussian [1.0 2.0 4.0 8.0 3.0 5.0 9.0 2.5] :lscv {:bounds [0.5 2.0]}))))
    (testing "a single point needs a bandwidth"
      (is (number? ((d/kernel-density :gaussian [1.0] {:bandwidth 1.0}) 1.0)))
      (is (invalid-input? :data #(d/kernel-density :gaussian [1.0]))))))

;; ---------------------------------------------------------------------------------------------
;; confidence intervals

(deftest kernel-density-ci-test
  (let [data [1.0 2.0 2.5 3.0 4.0 7.0]
        h 0.8
        za 1.959963984540054
        n (count data)]
    (testing "the formula: f -+ z * sqrt(k2 f / (n h))"
      (doseq [kernel kernels
              :let [f (d/kernel-density kernel data {:bandwidth h})
                    ci (d/kernel-density-ci kernel data {:bandwidth h})
                    k2 (:k2 (d/kde-data kernel))]
              x [0.0 2.0 3.5 7.5]
              :let [fx (f x)
                    band (* za (Math/sqrt (/ (* k2 (max fx 0.0)) (* n h))))
                    [e lo hi] (ci x)]]
        (is (close? e fx 1.0e-15) (str kernel))
        (is (close? lo (- fx band) 1.0e-12) (str kernel))
        (is (close? hi (+ fx band) 1.0e-12) (str kernel))))
    (testing "alpha"
      (let [[f lo hi] ((d/kernel-density-ci :gaussian data {:bandwidth h :alpha 0.1}) 3.0)]
        (is (close? (- hi f) (* 1.6448536269514722 (Math/sqrt (/ (* (:k2 (d/kde-data :gaussian)) f) (* n h)))) 1.0e-12))
        (is (close? (- f lo) (- hi f) 1.0e-12))))
    (testing "KE-17: the negative density of silverman gives no NaN"
      (let [ci (d/kernel-density-ci :silverman data {:bandwidth 0.3})
            negative-x (first (filter #(neg? ((d/kernel-density :silverman data {:bandwidth 0.3}) %))
                                      (range -10.0 20.0 0.05)))]
        (is (some? negative-x))
        (is (every? #(not (m/nan? %)) (ci negative-x)))
        (let [[f lo hi] (ci negative-x)]
          (is (neg? f))
          (is (= lo hi f)))))
    (testing "KE-11: invalid kernel and alpha throw ex-info, also with assertions off"
      (binding [*assert* false]
        (is (invalid-input? :kernel #(d/kernel-density-ci :foo data {:bandwidth h})))
        (is (invalid-input? :kernel #(d/kernel-density-ci (fn [x] x) data {:bandwidth h}))))
      (doseq [alpha [0 0.0 1 1.0 -0.1 1.5 ##NaN "a"]]
        (is (invalid-input? :alpha #(d/kernel-density-ci :gaussian data {:bandwidth h :alpha alpha})) (pr-str alpha))))
    (testing "nil alpha is the default 0.05"
      (is (= ((d/kernel-density-ci :gaussian data {:bandwidth h :alpha nil}) 3.0)
             ((d/kernel-density-ci :gaussian data {:bandwidth h}) 3.0))))))

;; ---------------------------------------------------------------------------------------------
;; wrappers in fastmath.kernel

(deftest kernel-wrappers-test
  (let [data [1.0 2.0 2.5 3.0 4.0 7.0 7.5 3.3]]
    (testing "KE-16: a keyword or a number as the third argument"
      (is (= (:h (k/kernel-density :gaussian data :nrd0 true)) (k/bandwidth :gaussian data :nrd0)))
      (is (= (:h (k/kernel-density :gaussian data 0.7 true)) 0.7))
      (is (= ((k/kernel-density :gaussian data :nrd0) 3.0)
             ((k/kernel-density :gaussian data {:bandwidth :nrd0}) 3.0)))
      (is (= ((k/kernel-density-ci :gaussian data :nrd0) 3.0)
             ((k/kernel-density-ci :gaussian data {:bandwidth :nrd0}) 3.0))))
    (testing "bandwidth-options through the wrapper"
      (is (= (k/bandwidth :gaussian data :lscv {:bounds [0.5 2.0]})
             (:h (k/kernel-density :gaussian data {:bandwidth :lscv :bandwidth-options {:bounds [0.5 2.0]}} true)))))
    (testing "the wrappers throw ex-info for invalid input"
      (is (invalid-input? :kernel #(k/kernel-density :foo data)))
      (is (invalid-input? :binned? #(k/kernel-density :gaussian data {:binned? -1})))
      (is (invalid-input? :alpha #(k/kernel-density-ci :gaussian data {:alpha 2}))))))

;; ---------------------------------------------------------------------------------------------
;; documentation

(deftest documentation-test
  (testing "every public var has a docstring in the house style (no double quotes, backslashes, code blocks, headings)"
    (doseq [[sym var] (ns-publics 'fastmath.kernel.density)
            :let [doc (:doc (meta var))]]
      (is (some? doc) (str sym))
      (is (not (re-find #"[\"\\]|```|(?m)^#" doc)) (str sym))
      (when (fn? @var)
        (is (re-find #"Returns" doc) (str sym " Returns"))
        (is (re-find #"See also" doc) (str sym " See also")))))
  (testing "the namespace has a docstring naming the functions"
    (let [doc (:doc (meta (the-ns 'fastmath.kernel.density)))]
      (is (every? #(re-find (re-pattern (str "\\[\\[" % "\\]\\]")) doc)
                  ["kernel-density" "kernel-density+" "kernel-density-ci" "bandwidth" "kde-data"]))))
  (testing "public API: functions in kde-data and the three estimators plus bandwidth"
    (is (= #{'kde-data 'bandwidth 'kernel-density 'kernel-density+ 'kernel-density-ci
             'uniform 'gaussian 'triangular 'epanechnikov 'quartic 'triweight 'tricube 'cosine 'logistic
             'sigmoid 'silverman 'laplace 'wigner 'cauchy}
           (set (keys (ns-publics 'fastmath.kernel.density)))))))
