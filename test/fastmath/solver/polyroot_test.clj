(ns fastmath.solver.polyroot-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.polynomials :as poly]
            [fastmath.solver.polyroot :as pr])
  (:import [clojure.lang ExceptionInfo]
           [fastmath.vector Vec2]))

(defn- thrown-data [thunk] (try (thunk) ::no-throw (catch ExceptionInfo e (ex-data e))))

(defn- pairs [roots] (mapv (fn [^Vec2 z] [(.x z) (.y z)]) roots))

(defn- near?
  "Two complex numbers `[re im]` equal within the relative tolerance."
  [[ar ai] [br bi] tol]
  (<= (m/hypot (- ar br) (- ai bi)) (* tol (max (m/hypot br bi) Double/MIN_NORMAL))))

(defn- match-all?
  "Every expected root has a computed root within `tol`, and the numbers of roots are equal."
  [expected computed tol]
  (and (= (count expected) (count computed))
       (loop [exp expected free (vec computed)]
         (or (empty? exp)
             (let [e (first exp)
                   best (when (seq free) (apply min-key (fn [c] (m/hypot (- (c 0) (e 0)) (- (c 1) (e 1)))) free))]
               (and best (near? best e tol)
                    (recur (rest exp) (let [i (.indexOf ^java.util.List free best)] (into (subvec free 0 i) (subvec free (inc i)))))))))))

(defn- max-error
  "The largest relative distance of an expected root to its nearest computed root."
  [expected computed]
  (reduce max 0.0 (map (fn [e] (apply min (map (fn [c] (/ (m/hypot (- (c 0) (e 0)) (- (c 1) (e 1))) (m/hypot (e 0) (e 1)))) computed))) expected)))

(t/deftest sorted-roots-keep-all
  (let [sorted-roots #'pr/sorted-roots
        root (fn [re im] (pr/->Root re im 0.0 1.0 true))
        out (sorted-roots [(root 1.0 -2.0) (root 3.0 0.0)] true)]
    (t/is (= 2 (count out)) "a lower half plane root without a partner is not dropped")
    (t/is (= 3 (count (sorted-roots [(root 1.0 2.0) (root 1.0 -2.0) (root 5.0 -1.0)] true))))))

(t/deftest basic
  (let [rs (pairs (pr/polyroot [-6 11 -6 1]))]
    (t/is (= 3 (count rs)))
    (t/is (every? #(zero? (% 1)) rs) "real roots have exact zero imaginary parts")
    (t/is (every? true? (map #(<= (m/abs (- (%1 0) %2)) 1e-14) rs [1.0 2.0 3.0])))
    (t/is (instance? Vec2 (first (pr/polyroot [-6 11 -6 1]))))
    (t/is (vector? (pr/polyroot [-6 11 -6 1]))))
  (t/is (match-all? [[0.0 1.0] [0.0 -1.0]] (pairs (pr/polyroot [1 0 1])) 1e-14) "z^2 + 1")
  (t/is (match-all? [[2.0 0.0]] (pairs (pr/polyroot [-4 2])) 1e-15) "linear")
  (t/is (match-all? [[0.0 -1.0]] (pairs (pr/polyroot [[0 1] 1])) 1e-15) "i + z = 0 -> z = -i, complex coefficients")
  (t/is (match-all? [[1.0 0.0] [2.0 0.0]] (pairs (pr/polyroot [2 -3 1])) 1e-14) "quadratic")
  (t/is (match-all? [[1.5 0.5] [1.5 -0.5]] (pairs (pr/polyroot [2.5 -3 1])) 1e-14) "quadratic, complex roots"))

(t/deftest input-forms
  (let [expected (pairs (pr/polyroot [-6 11 -6 1]))]
    (t/is (= expected (pairs (pr/polyroot (double-array [-6 11 -6 1])))) "array")
    (t/is (= expected (pairs (pr/polyroot (list -6 11 -6 1)))) "list")
    (t/is (= expected (pairs (pr/polyroot [-6N 11 -6.0 1/1]))) "other number types")
    (t/is (= expected (pairs (pr/polyroot (poly/polynomial [-6 11 -6 1])))) "polynomial object")
    (t/is (= expected (pairs (pr/polyroot (poly/ratio-polynomial [-6 11 -6 1])))) "ratio polynomial object")
    (t/is (= expected (pairs (pr/polyroot [[-6 0] (Vec2. 11.0 0.0) -6 [1.0 0.0]]))) "pairs and Vec2")
    (t/is (= expected (pairs (pr/polyroot [-6 11 -6 1 0 0]))) "leading zero coefficients are dropped")))

(t/deftest boundaries
  (t/is (= [] (pr/polyroot [5])) "a constant has no roots")
  (t/is (= [] (pr/polyroot [[0 3]])) "a complex constant")
  (t/is (= [[0.0 0.0]] (pairs (pr/polyroot [0 7]))) "7 z")
  (t/is (= [[0.0 0.0] [0.0 0.0] [0.0 0.0]] (pairs (pr/polyroot [0 0 0 2]))) "2 z^3")
  (t/is (= [[0.0 0.0] [0.0 0.0] [0.5 0.0] [1.0 0.0]] (pairs (pr/polyroot [0 0 1 -3 2]))) "trailing zeros are exact zero roots")
  (t/is (= [] (pr/polyroot [3 0 0 0])) "leading zeros of a constant")
  ;; the zero polynomial and invalid input
  (t/is (= :coefficients (:option (thrown-data #(pr/polyroot [0 0 0])))))
  (t/is (= :coefficients (:option (thrown-data #(pr/polyroot [])))))
  (t/is (= :coefficients (:option (thrown-data #(pr/polyroot nil)))))
  (t/is (= :coefficients (:option (thrown-data #(pr/polyroot 5)))))
  (doseq [bad [##NaN ##Inf ##-Inf "1" :a nil [1] [1 2 3] [1 ##NaN]]]
    (let [d (thrown-data #(pr/polyroot [1 bad 1]))]
      (t/is (= :coefficients (:option d)) (str (pr/polyroot [1 1]) " invalid " (pr-str bad)))
      (t/is (= 1 (:index d)) "the index of the coefficient")))
  (t/is (= :max-iters (:option (thrown-data #(pr/polyroot [1 2 1] {:max-iters 0})))))
  (t/is (= :max-iters (:option (thrown-data #(pr/polyroot [1 2 1] {:max-iters 2.5})))))
  ;; the error data never holds the coefficients
  (let [coefs (vec (concat (repeat 100000 1.0) [##NaN]))
        d (thrown-data #(pr/polyroot coefs))]
    (t/is (= 100000 (:index d)))
    (t/is (< (count (pr-str d)) 300) "data is small")))

(defn- from-roots
  "Coefficients (ascending) of the monic polynomial with the given complex roots `[re im]`."
  [roots]
  (let [mul (fn [co [rr ri]]
              ;; co * (z - r), co ascending as [re im]
              (let [shifted (cons [0.0 0.0] co)
                    scaled (concat (map (fn [[a b]] [(- (+ (* a rr) (- (* b ri)))) (- (+ (* a ri) (* b rr)))]) co) [[0.0 0.0]])]
                (mapv (fn [[a b] [c d]] [(+ a c) (+ b d)]) shifted scaled)))]
    (reduce mul [[1.0 0.0]] roots)))

(t/deftest ordering
  ;; real coefficients: real roots ascending, then the pairs of complex conjugates by the real part, the positive imaginary part first
  (let [co (map first (from-roots [[3.0 0.0] [-2.0 0.0] [1.0 2.0] [1.0 -2.0] [-1.0 1.0] [-1.0 -1.0] [0.5 0.0]]))
        rs (pairs (pr/polyroot co))]
    (t/is (= 7 (count rs)))
    (t/is (every? zero? (map #(% 1) (take 3 rs))) "three real roots first")
    (t/is (apply < (map first (take 3 rs))) "ascending")
    (t/is (every? true? (map #(<= (m/abs (- (%1 0) %2)) 1e-12) (take 3 rs) [-2.0 0.5 3.0])))
    (let [[a b c d] (drop 3 rs)]
      (t/is (and (pos? (a 1)) (neg? (b 1)) (== (a 0) (b 0)) (== (a 1) (- (b 1)))) "an exact conjugate pair, positive first")
      (t/is (and (pos? (c 1)) (neg? (d 1)) (== (c 0) (d 0)) (== (c 1) (- (d 1)))))
      (t/is (< (a 0) (c 0)) "pairs by the real part: -1 before 1")))
  ;; complex coefficients: by the real part, then the imaginary
  (let [co (from-roots [[2.0 1.0] [-1.0 3.0] [-1.0 -3.0] [0.0 0.0] [5.0 -1.0]])
        rs (pairs (pr/polyroot co))]
    (t/is (= rs (sort-by (fn [[re im]] [re im]) rs)))
    (t/is (match-all? [[2.0 1.0] [-1.0 3.0] [-1.0 -3.0] [0.0 0.0] [5.0 -1.0]] rs 1e-10))))

(t/deftest roots-of-unity
  (doseq [n [2 3 4 5 8 16 31 64]
          :let [rs (pairs (pr/polyroot (concat [-1] (repeat (dec n) 0) [1])))]]
    (t/is (= n (count rs)))
    (t/is (every? (fn [[re im]] (<= (m/abs (- (m/hypot re im) 1.0)) 1e-13)) rs) (str n " roots on the unit circle"))
    (t/is (match-all? (for [k (range n)] [(m/cos (/ (* 2 m/PI k) n)) (m/sin (/ (* 2 m/PI k) n))]) rs 1e-12) (str "z^" n " - 1"))))

(t/deftest special-polynomials
  ;; Wilkinson polynomial of degree 10: exact coefficients, roots 1..10
  (let [w10 (map first (from-roots (map #(vector (double %) 0.0) (range 1 11))))
        rs (pairs (pr/polyroot w10))]
    (t/is (match-all? (map #(vector (double %) 0.0) (range 1 11)) rs 1e-6) "Wilkinson 10, roots 1..10")
    (t/is (every? #(zero? (% 1)) rs) "real roots are real"))
  ;; Chebyshev T_12: cos((2k - 1) pi / 24)
  (let [t (loop [n 1 a [1] b [0 1]]
            (if (== n 12) b
                (recur (inc n) b (let [x2b (concat [0] (map #(* 2 %) b))
                                       a (concat a (repeat (- (count x2b) (count a)) 0))]
                                   (mapv - x2b a)))))
        rs (pairs (pr/polyroot t))]
    (t/is (match-all? (for [k (range 1 13)] [(m/cos (/ (* (- (* 2 k) 1) m/PI) 24)) 0.0]) rs 1e-12) "Chebyshev T12"))
  ;; a multiple root: the accuracy is the fifth root of the machine epsilon
  (let [rs (pairs (pr/polyroot [-1 5 -10 10 -5 1]))]
    (t/is (= 5 (count rs)))
    (t/is (every? (fn [[re im]] (<= (m/hypot (- re 1.0) im) 5e-3)) rs) "(z-1)^5"))
  ;; a cluster of close roots
  (let [roots (map #(vector (+ 1.0 (* 1e-3 %)) 0.0) (range 5))
        {rs :roots be :backward-errors} (pr/polyroot (map first (from-roots roots)) {:info? true})]
    ;; the roots 1e-3 apart are known to 1e-16 / (1e-3)^4 relative: the backward error is the measure of the method
    (t/is (match-all? roots (pairs rs) 2e-3) "five roots 1e-3 apart, to the conditioning")
    ;; the stopping rule weights the terms by 3.8 k + 1: up to 20 eps for the degree 5, observed 8.5 eps
    (t/is (every? #(<= % (* 20 2.220446049250313E-16)) be) (str "backward errors at the rounding level " (apply max be))))
  ;; very different sizes of the roots, the reversed evaluation of the roots outside of the unit circle
  (let [rs (pairs (pr/polyroot [-1.0 (+ 1e100 1.0) (- (+ 1e100 1.0)) 1.0]))]
    (t/is (match-all? [[1e100 0.0] [1.0 0.0] [1e-100 0.0]] rs 1e-12) "roots 1e100, 1, 1e-100"))
  (t/is (match-all? [[1e100 0.0] [-1e100 0.0]] (pairs (pr/polyroot [-1e200 0 1])) 1e-13) "huge roots")
  (t/is (match-all? [[1e-100 0.0] [-1e-100 0.0]] (pairs (pr/polyroot [-1e-200 0 1])) 1e-13) "tiny roots")
  ;; roots of the size 1e60 give the coefficients up to 1e240
  (let [roots [[3e60 0.0] [-3e60 0.0] [1e60 2e60] [1e60 -2e60]]]
    (t/is (match-all? roots (pairs (pr/polyroot (map first (from-roots roots)))) 1e-10) "coefficients up to 1e240")))

(t/deftest scale-invariance
  ;; roots of c p(z): the same roots when the coefficients are multiplied by a constant of any size
  (let [co [-6 11 -6 1]
        base (pairs (pr/polyroot co))]
    (doseq [k [-300 -100 100 300]
            :let [s (m/scalb 1.0 k)]]
      (t/is (match-all? base (pairs (pr/polyroot (map #(* s %) co))) 1e-14) (str "scaled by 2^" k)))))

(t/deftest backward-error
  (let [{:keys [roots backward-errors condition-numbers converged? status]}
        (pr/polyroot (map first (from-roots (map #(vector (double %) 0.0) (range 1 21)))) {:info? true})]
    (t/is (= :converged status))
    (t/is (= 20 (count roots)))
    (t/is (every? true? converged?))
    ;; Wilkinson 20 is so ill conditioned that its roots are not the integers, but each root is exact for a coefficient perturbation of the rounding size
    (t/is (every? #(<= % (* 100 2.220446049250313E-16)) backward-errors) (str "backward errors " (apply max backward-errors)))
    (t/is (every? pos? condition-numbers))))

(t/deftest info
  (let [info (pr/polyroot [-6 11 -6 1] {:info? true})]
    (t/is (= #{:roots :backward-errors :condition-numbers :converged? :iterations :status} (set (keys info))))
    (t/is (= :converged (:status info)))
    (t/is (= 3 (count (:roots info)) (count (:backward-errors info)) (count (:condition-numbers info)) (count (:converged? info))))
    (t/is (every? #(<= 0.0 % 1e-15) (:backward-errors info)))
    (t/is (<= (apply max (:backward-errors info)) 2.220446049250313E-16)))
  (let [info (pr/polyroot [0 0 1 -3 2] {:info? true})]
    (t/is (= [##Inf ##Inf] (take 2 (:condition-numbers info))) "a zero root has an infinite condition number")
    (t/is (= [0.0 0.0] (take 2 (:backward-errors info)))))
  (t/is (= {:status :converged :roots []} (select-keys (pr/polyroot [5] {:info? true}) [:status :roots])))
  ;; the condition number of a root of z^2 - 1 is 1: a relative change of the coefficients by e changes the root by about e
  (t/is (m/delta-eq 1.0 (first (:condition-numbers (pr/polyroot [-1 0 1] {:info? true}))) 1e-12)))

(t/deftest limits
  (let [co (map first (from-roots (map #(vector (double %) 0.0) (range 1 16))))
        d (thrown-data #(pr/polyroot co {:max-iters 1}))]
    (t/is (= :max-evals (:status d)))
    (t/is (pos? (:unconverged d)))
    (t/is (= 15 (:degree d)))
    (let [info (pr/polyroot co {:max-iters 1 :info? true})]
      (t/is (= :max-evals (:status info)))
      (t/is (= 15 (count (:roots info))))
      (t/is (some false? (:converged? info))))))

(t/deftest deterministic
  (let [co [1 2 3 4 5 6 7]]
    (t/is (= (pairs (pr/polyroot co)) (pairs (pr/polyroot co))))
    (t/is (= 1 (count (distinct (pmap (fn [_] (pairs (pr/polyroot co))) (range 32))))) "thread safe")))

(def reference (delay (edn/read-string (slurp (io/resource "solver/polyroot_reference.edn")))))

;; mpmath roots (150 digits) of random polynomials with well separated roots,
;; generated by utils/fastmath/dev/generate_solver_polyroot_reference.py
(t/deftest reference-roots
  (let [worst (atom 0.0)]
    (doseq [{:keys [kind coefs roots]} @reference
            :let [got (pairs (pr/polyroot coefs))
                  err (max-error roots got)]]
      (swap! worst max err)
      (t/is (= (count roots) (count got)) (str kind " " (count coefs) " coefficients"))
      (t/is (<= err 1e-9) (str kind " degree " (dec (count coefs)) " error " err))
      (when (= kind :real)
        (let [reals (filter #(zero? (% 1)) got)
              expected-reals (filter #(zero? (% 1)) roots)]
          (t/is (= (count expected-reals) (count reals)) "the same number of real roots"))))
    ;; the observed maximum of the relative error over the 160 polynomials is checked against the tolerance 1e-9
    (t/is (<= @worst 1e-9) (str "worst error " @worst))))
