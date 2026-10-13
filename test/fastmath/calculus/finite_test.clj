(ns fastmath.calculus.finite-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.calculus.finite :as sut]))

(set! *warn-on-reflection* true)

;; helpers

(defn- thrown-data
  "ex-data of the `ExceptionInfo` thrown by `thunk`, nil when nothing is thrown."
  [thunk]
  (try (thunk) nil
       (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- invalid?
  "Is the `ex-info` thrown with the given option and value and an expectation text?"
  [thunk option value]
  (let [d (thrown-data thunk)]
    (and (= option (:option d))
         (= value (:value d))
         (string? (:expected d)))))

(defn- invalid-option?
  [thunk option]
  (= option (:option (thrown-data thunk))))

(defn- close?
  "Absolute error when the expected value is below 1, relative otherwise."
  [^double a ^double b ^double tol]
  (<= (Math/abs (- a b)) (* tol (Math/max 1.0 (Math/abs b)))))

(defn- all-close?
  [as bs tol]
  (and (= (count as) (count bs))
       (every? true? (map #(close? %1 %2 tol) as bs))))

(defn- rel-err
  ^double [^double a ^double b]
  (if (zero? b) (Math/abs a) (/ (Math/abs (- a b)) (Math/abs b))))

;; ------------------------------------------------------------------ coefficients

(def ^:private fd-reference
  (delay (edn/read-string (slurp (io/resource "calculus/fd_coeffs_reference.edn")))))

(t/deftest fd-coeffs-reference
  (t/testing "168 combinations (n 1-7, acc 1-8, 3 kinds) against sympy finite_diff_weights"
    (t/is (= 168 (count @fd-reference)))
    (doseq [{:keys [n acc kind offsets coeffs]} @fd-reference]
      (let [[o c] (sut/fd-coeffs-for-accuracy n acc kind)]
        (t/is (= offsets (vec o)) (str "offsets " [n acc kind]))
        (t/is (= (count coeffs) (count c)))
        ;; measured max relative difference: 0.0 (the weights are the nearest doubles of the exact ratios)
        (t/is (every? true? (map (fn [^double a ^double b]
                                   (if (zero? b)
                                     (zero? a)
                                     (<= (Math/abs (- a b)) (* 1.0e-15 (Math/abs b)))))
                                 c coeffs))
              (str "coeffs " [n acc kind]))))))

(t/deftest fd-coeffs-arities
  (t/is (= (sut/fd-coeffs-for-accuracy 1 2 :central) (sut/fd-coeffs-for-accuracy 1)))
  (t/is (= (sut/fd-coeffs-for-accuracy 2 2 :central) (sut/fd-coeffs-for-accuracy 2)))
  (t/is (= (sut/fd-coeffs-for-accuracy 3 4 :central) (sut/fd-coeffs-for-accuracy 3 4))))

(t/deftest fd-coeffs-known-values
  (t/is (= [[-1 0 1] [-0.5 0.0 0.5]] (let [[o c] (sut/fd-coeffs 1 [-1 0 1])] [o (vec c)])))
  (t/is (= [1.0 -2.0 1.0] (vec (second (sut/fd-coeffs 2 [-1 0 1])))))
  (t/is (= [-1.5 2.0 -0.5] (vec (second (sut/fd-coeffs 1 [0 1 2])))))
  (t/is (= [-1.0 3.0 -3.0 1.0] (vec (second (sut/fd-coeffs 3 [0 1 2 3])))))
  (t/testing "n = 0 is the interpolation: weight 1 at the offset 0"
    (t/is (= [0.0 1.0 0.0] (vec (second (sut/fd-coeffs 0 [-1 0 1]))))))
  (t/testing "single offset"
    (t/is (= [1.0] (vec (second (sut/fd-coeffs 0 [0])))))
    (t/is (= [0.0] (vec (second (sut/fd-coeffs 1 [0]))))))
  (t/testing "fewer offsets than n + 1: the derivative is not determined, zero weights"
    (t/is (= [0.0 0.0] (vec (second (sut/fd-coeffs 2 [0 1]))))))
  (t/testing "non-integer and unordered offsets"
    (t/is (all-close? (second (sut/fd-coeffs 1 [-0.5 0.5])) [-1.0 1.0] 1.0e-15))
    (t/is (all-close? (second (sut/fd-coeffs 2 [1 -1 0])) [1.0 1.0 -2.0] 1.0e-15))
    (t/is (all-close? (second (sut/fd-coeffs 1 [0 1/2 1])) [-3.0 4.0 -1.0] 1.0e-15)))
  (t/testing "offsets are returned as given"
    (t/is (= [1 -1 0] (first (sut/fd-coeffs 2 [1 -1 0]))))))

(t/deftest fd-coeffs-moment-conditions
  (t/testing "high orders: sum(c_i * x_i^k) = k! if k = n else 0, for k below the number of nodes"
    (doseq [[n offsets] [[9 (range -6 7)] [10 (range 0 14)] [12 (range -8 9)]]
            :let [[_ c] (sut/fd-coeffs n offsets)]
            k (range (count offsets))
            :let [moment (reduce + (map (fn [x ^double w] (* w (Math/pow (double x) k))) offsets c))
                  scale (reduce + (map (fn [x ^double w] (Math/abs (* w (Math/pow (double x) k)))) offsets c))
                  expected (if (= k n) (m/factorial n) 0.0)]]
      ;; rounding of the sum of terms of size `scale`
      (t/is (<= (Math/abs (- moment expected)) (* 1.0e-13 (Math/max 1.0 (double scale))))
            (str "n=" n " k=" k)))))

;; ------------------------------------------------------------------ validation

(t/deftest validation
  (t/testing "derivative"
    (t/is (invalid? #(sut/derivative identity -1) :n -1))
    (t/is (invalid? #(sut/derivative identity 1.5) :n 1.5))
    (t/is (invalid? #(sut/derivative identity 1 {:acc 0}) :acc 0))
    (t/is (invalid? #(sut/derivative identity 1 {:acc -2}) :acc -2))
    (t/is (invalid? #(sut/derivative identity 1 {:acc 2.5}) :acc 2.5))
    (t/is (invalid? #(sut/derivative identity 1 {:method :bogus}) :method :bogus))
    (t/is (invalid? #(sut/derivative identity 1 {:h -1.0}) :h -1.0))
    (t/is (invalid-option? #(sut/derivative identity 1 {:h ##NaN}) :h))
    (t/is (invalid-option? #(sut/derivative identity 1 {:h ##Inf}) :h))
    (t/is (invalid? #(sut/derivative identity 1 {:h "a"}) :h "a")))
  (t/testing "derivative accepts the boundary values"
    (t/is (fn? (sut/derivative identity 0)))
    (t/is (fn? (sut/derivative identity 1 {:acc 1})))
    (t/is (fn? (sut/derivative identity 1 {:h 0})))
    (t/is (fn? (sut/derivative identity 1 {:h 0.0 :acc 1 :method :forward}))))
  (t/testing "fd-coeffs-for-accuracy"
    (t/is (invalid? #(sut/fd-coeffs-for-accuracy -1) :n -1))
    (t/is (invalid? #(sut/fd-coeffs-for-accuracy 1 0) :acc 0))
    (t/is (invalid? #(sut/fd-coeffs-for-accuracy 1 2 :bogus) :method :bogus)))
  (t/testing "fd-coeffs"
    (t/is (invalid? #(sut/fd-coeffs -1 [0 1]) :n -1))
    (t/is (invalid? #(sut/fd-coeffs 1 []) :offsets []))
    (t/is (invalid? #(sut/fd-coeffs 1 [1 2 1]) :offsets [1 2 1]))
    (t/is (invalid? #(sut/fd-coeffs 1 [0.5 1/2]) :offsets [0.5 1/2]))
    (t/is (invalid-option? #(sut/fd-coeffs 1 [0 ##Inf]) :offsets)))
  (t/testing "gradient"
    (t/is (invalid? #(sut/gradient identity {:acc 0}) :acc 0))
    (t/is (invalid? #(sut/gradient identity {:acc 2.0}) :acc 2.0))
    (t/is (invalid? #(sut/gradient identity {:h -1.0}) :h -1.0))
    (t/is (invalid? #((sut/gradient (fn [_] 0.0)) [1.0] -1.0) :h -1.0)))
  (t/testing "jacobian"
    (t/is (invalid? #(sut/jacobian identity {:acc 0}) :acc 0))
    (t/is (invalid? #(sut/jacobian identity {:h -1.0}) :h -1.0)))
  (t/testing "hessian, hessian-diagonal, laplacian"
    (t/is (invalid? #(sut/hessian identity {:h -1.0}) :h -1.0))
    (t/is (invalid-option? #(sut/hessian identity {:h ##NaN}) :h))
    (t/is (invalid? #(sut/hessian-diagonal identity {:h -1.0}) :h -1.0))
    (t/is (invalid? #(sut/laplacian identity {:h -1.0}) :h -1.0))
    (t/is (invalid? #((sut/hessian (fn [_] 0.0)) [1.0] -2.0) :h -2.0))))

(t/deftest extrapolate-validation
  (t/is (invalid? #(sut/extrapolate identity {:contract 1.0}) :contract 1.0))
  (t/is (invalid? #(sut/extrapolate identity {:contract 2.0}) :contract 2.0))
  (t/is (invalid? #(sut/extrapolate identity {:contract 0.0}) :contract 0.0))
  (t/is (invalid? #(sut/extrapolate identity {:contract -0.5}) :contract -0.5))
  (t/is (invalid-option? #(sut/extrapolate identity {:contract ##NaN}) :contract))
  (t/is (invalid? #(sut/extrapolate identity {:max-evals 0}) :max-evals 0))
  (t/is (invalid? #(sut/extrapolate identity {:max-evals -3}) :max-evals -3))
  (t/is (invalid? #(sut/extrapolate identity {:max-evals 2.5}) :max-evals 2.5))
  (t/is (invalid? #(sut/extrapolate identity {:init-h 0.0}) :init-h 0.0))
  (t/is (invalid? #(sut/extrapolate identity {:init-h 0}) :init-h 0))
  (t/is (invalid-option? #(sut/extrapolate identity {:init-h ##NaN}) :init-h))
  (t/is (invalid-option? #(sut/extrapolate identity {:init-h ##Inf}) :init-h))
  (t/testing "valid options are accepted"
    (t/is (fn? (sut/extrapolate identity {:contract 0.9 :max-evals 1 :init-h -0.5})))
    (t/is (fn? (sut/extrapolate identity)))
    (t/is (fn? (sut/extrapolate identity nil)))))

;; ------------------------------------------------------------------ derivative

(defn- exp-fn ^double [^double x] (Math/exp x))
(defn- sin-fn ^double [^double x] (Math/sin x))

(defn- exact-derivative
  ^double [kind ^long n ^double x]
  (case kind
    :exp (Math/exp x)
    :sin (Math/sin (+ x (* n m/HALF_PI)))))

(t/deftest derivative-of-order-zero
  (t/is (identical? exp-fn (sut/derivative exp-fn 0)))
  (t/is (identical? exp-fn (sut/derivative exp-fn 0 {:acc 4}))))

(t/deftest derivative-arities
  (let [d (sut/derivative exp-fn 2)
        dh (sut/derivative exp-fn 2 {:h 1.0e-3})
        h 0.01]
    (t/testing "1-arity, 2-arity with the step, 3-arity with the step and its power"
      (t/is (number? (d 0.5)))
      (t/is (= (d 0.5 h) (d 0.5 h (m/fpow h 2))))
      (t/is (= (d 0.5 h) (dh 0.5 h))))
    (t/testing "explicit :h option is an absolute step"
      (t/is (= (dh 0.5) (d 0.5 1.0e-3)))
      (t/is (= (dh 5.0) (d 5.0 1.0e-3)))
      (t/is (= (dh 30.0) (d 30.0 1.0e-3))))
    (t/testing "zero step is the automatic one"
      (t/is (= (d 0.5) (d 0.5 0.0)))
      (t/is (= (d 5.0) (d 5.0 0.0)))
      (t/is (= (dh 0.5) ((sut/derivative exp-fn 2 {:h 0.0}) 0.5 1.0e-3)))
      (t/is (= (d 0.5) ((sut/derivative exp-fn 2 {:h 0}) 0.5))))
    (t/testing "the 2-arity function is what extrapolate needs"
      (t/is (close? ((sut/extrapolate (sut/derivative exp-fn 1) {:power 2.0}) 0.5)
                    (Math/exp 0.5) 1.0e-9)))))

(t/deftest derivative-stencil-consistency
  (t/testing "every inline stencil (n 1-6, acc 1-6) and the general one (n 7-8, acc 7-8) equals the sum of fd-coeffs"
    (doseq [n (range 1 9)
            acc (range 1 9)
            method [:central :forward :backward]
            :let [h 0.05
                  x 0.3
                  [offsets coeffs] (sut/fd-coeffs-for-accuracy n acc method)
                  terms (map (fn [o ^double c] (* c (Math/exp (+ x (* ^long o h))))) offsets coeffs)
                  expected (/ (reduce + terms) (m/fpow h n))
                  scale (/ (reduce + (map #(Math/abs ^double %) terms)) (m/fpow h n))
                  d (sut/derivative exp-fn n {:acc acc :method method})]]
      ;; both are sums of the same terms, only the order of the additions differs
      (t/is (<= (Math/abs (- ^double (d x h) expected)) (* 1.0e-13 scale))
            (str "n=" n " acc=" acc " " method)))))

(t/deftest derivative-accuracy-analytic
  (t/testing "all orders with acc 4 and an explicit step"
    (doseq [kind [:exp :sin]
            n (range 1 7)
            x [-2.0 0.0 0.3 1.7]
            :let [f (if (= kind :exp) exp-fn sin-fn)
                  d (sut/derivative f n {:acc 4 :h (if (> n 3) 0.1 0.02)})
                  exact (exact-derivative kind n x)]]
      (t/is (close? (d x) exact 1.0e-3) (str kind " n=" n " x=" x))))
  (t/testing "default step, central acc 2 and 4, n 1-3, relative errors measured below 1.0e-6"
    (doseq [kind [:exp :sin]
            n [1 2 3]
            acc [2 4]
            x [-3.0 -0.5 0.0 0.3 2.0 5.0]
            :let [f (if (= kind :exp) exp-fn sin-fn)
                  d (sut/derivative f n {:acc acc})
                  exact (exact-derivative kind n x)]]
      (t/is (close? (d x) exact 1.0e-5) (str kind " n=" n " acc=" acc " x=" x))))
  (t/testing "forward and backward methods"
    (doseq [method [:forward :backward]
            n [1 2 3]
            :let [d (sut/derivative exp-fn n {:acc 3 :method method})]]
      (t/is (close? (d 0.5) (Math/exp 0.5) 1.0e-3) (str method " n=" n))
      (t/is (close? (d -4.0) (Math/exp -4.0) 1.0e-3) (str method " n=" n))))
  (t/testing "polynomials are exact up to the order of the stencil"
    (let [cubic (fn ^double [^double x] (+ (* x x x) (* -2.0 x x) x 3.0))]
      (t/is (close? ((sut/derivative cubic 1 {:acc 4 :h 0.1}) 2.0) 5.0 1.0e-10))
      (t/is (close? ((sut/derivative cubic 2 {:acc 2 :h 0.1}) 2.0) 8.0 1.0e-9))
      (t/is (close? ((sut/derivative cubic 3 {:acc 2 :h 0.5}) 2.0) 6.0 1.0e-9))
      (t/is (close? ((sut/derivative cubic 3 {:acc 1 :h 0.5 :method :forward}) 2.0) 6.0 1.0e-9)))))

(t/deftest derivative-order-of-accuracy
  (t/testing "the error shrinks as h^acc: ratio of the errors for h and h/2 is close to 2^acc"
    (let [ratio (fn [n acc method h]
                  (let [e (fn [h] (Math/abs (- ^double ((sut/derivative exp-fn n {:acc acc :h h :method method}) 0.3)
                                                (Math/exp 0.3))))]
                    (/ (e h) (e (/ h 2.0)))))]
      (t/are [n acc method h lo hi] (<= lo (ratio n acc method h) hi)
        1 2 :central 0.1 3.9 4.1
        2 2 :central 0.1 3.9 4.1
        3 2 :central 0.2 3.9 4.1
        1 4 :central 0.1 15.5 16.5
        4 3 :central 0.2 15.5 16.5
        1 1 :forward 0.05 1.9 2.1
        1 2 :forward 0.1 3.9 4.4
        1 2 :backward 0.1 3.7 4.2
        ;; the former LU coefficients gave a negative observed order here
        5 7 :forward 0.2 100.0 400.0))))

(t/deftest derivative-high-orders
  (t/testing "n = 7, 8 (general stencil) are not garbage with the automatic step"
    ;; measured: n=8 acc=2: 0.99736, acc=4: 1.00076, n=7 acc=8: 0.99999998 (exact 1)
    (t/is (close? ((sut/derivative exp-fn 8) 0.0) 1.0 1.0e-2))
    (t/is (close? ((sut/derivative exp-fn 8 {:acc 4}) 0.0) 1.0 1.0e-2))
    (t/is (close? ((sut/derivative exp-fn 7 {:acc 8}) 0.0) 1.0 1.0e-6))))

(def ^:private old-default-errors
  "Absolute errors of the former default step (eps^(1/(n+2))) for acc 4 and 6 (central), [function x n acc error].
  Measured with the code before the audit."
  [[:exp 0.3 1 4 3.7718939083219993E-11]
   [:exp 0.3 1 6 7.036815574679167E-12]
   [:exp 0.3 2 4 2.683746536114029E-8]
   [:exp 0.3 2 6 2.4697195621925516E-8]
   [:exp 0.3 3 4 5.441561909957926E-7]
   [:exp 0.3 3 6 1.9329276072399537E-5]
   [:exp 0.3 4 4 5.308172108198228E-5]
   [:exp 0.3 4 6 4.1841072250692335E-5]
   [:exp 0.7 1 4 2.5717650231626976E-11]
   [:exp 0.7 1 6 8.346701108052912E-11]
   [:exp 0.7 2 4 3.8003273417075434E-8]
   [:exp 0.7 2 6 4.393017416504108E-8]
   [:exp 0.7 3 4 4.942254152595638E-7]
   [:exp 0.7 3 6 3.0172089121371926E-5]
   [:exp 0.7 4 4 1.6772561648492612E-4]
   [:exp 0.7 4 6 2.0407263545063614E-4]
   [:exp 5.0 1 4 4.916302032142994E-9]
   [:exp 5.0 1 6 1.3324608971743146E-8]
   [:exp 5.0 2 4 1.474333004125583E-6]
   [:exp 5.0 2 6 1.1160848032432114E-6]
   [:exp 5.0 3 4 7.596593801224572E-4]
   [:exp 5.0 3 6 0.003193869443634867]
   [:exp 5.0 4 4 0.00808601248775176]
   [:exp 5.0 4 6 0.010485211483796775]
   [:sin 0.3 1 4 4.1701087027945505E-12]
   [:sin 0.3 1 6 7.41295913542217E-13]
   [:sin 0.3 2 4 4.812803655607212E-9]
   [:sin 0.3 2 6 4.689326593698695E-9]
   [:sin 0.3 3 4 2.6297906963179685E-8]
   [:sin 0.3 3 6 4.9290879445829106E-6]
   [:sin 0.3 4 4 9.594448170668901E-6]
   [:sin 0.3 4 6 1.2973819768724582E-5]
   [:sin 0.7 1 4 4.807265696626928E-13]
   [:sin 0.7 1 6 1.5279000287193867E-11]
   [:sin 0.7 2 4 2.428900414130908E-8]
   [:sin 0.7 2 6 1.1776658537243634E-8]
   [:sin 0.7 3 4 9.152481087459918E-7]
   [:sin 0.7 3 6 1.0657381026346613E-5]
   [:sin 0.7 4 4 3.1960388147012075E-5]
   [:sin 0.7 4 6 1.784212458166845E-5]
   [:sin 5.0 1 4 1.2193468457155632E-11]
   [:sin 5.0 1 6 4.772149342358034E-11]
   [:sin 5.0 2 4 2.7134956948060562E-8]
   [:sin 5.0 2 6 1.2646977798347336E-8]
   [:sin 5.0 3 4 1.5353354136937725E-6]
   [:sin 5.0 3 6 1.2059545142439898E-5]
   [:sin 5.0 4 4 2.200872111557395E-6]
   [:sin 5.0 4 6 4.1533979238117524E-6]])

(t/deftest derivative-default-step-dominance
  (t/testing "acc 4 and 6: the new default step is better than the former one in all 48 cells (n <= 4, x in 0.3 0.7 5)"
    (doseq [[kind x n acc old-error] old-default-errors
            :let [f (if (= kind :exp) exp-fn sin-fn)
                  d (sut/derivative f n {:acc acc})
                  new-error (Math/abs (- ^double (d x) (exact-derivative kind n x)))]]
      (t/is (<= new-error old-error) (str kind " x=" x " n=" n " acc=" acc " new " new-error " old " old-error))))
  (t/testing "acc 2: up to |x| = 10 the default step is the former one, eps^(1/(n+2)), so the error is the same in all 24 cells"
    (doseq [kind [:exp :sin]
            n (range 1 5)
            x [0.3 0.7 5.0 10.0]
            :let [f (if (= kind :exp) exp-fn sin-fn)
                  former (Math/pow m/MACHINE-EPSILON (/ 1.0 (+ n 2.0)))]]
      (t/is (= ((sut/derivative f n) x) ((sut/derivative f n) x former)) (str kind " x=" x " n=" n))))
  (t/testing "above |x| = 10 the step grows with |x|: 50 is 5 times the step at 10"
    (let [d (sut/derivative exp-fn 1)
          h (Math/pow m/MACHINE-EPSILON (/ 1.0 3.0))]
      (t/is (= (d 50.0) (d 50.0 (* 5.0 h))))
      (t/is (= (d 10.0) (d 10.0 h)))))
  (t/testing "higher accuracy is not worse than acc 2 at the default step (n = 1)"
    (doseq [x [0.3 0.7 5.0]
            acc [4 6 8]]
      (t/is (<= (rel-err ((sut/derivative exp-fn 1 {:acc acc}) x) (Math/exp x))
                (+ 1.0e-10 (rel-err ((sut/derivative exp-fn 1 {:acc 2}) x) (Math/exp x))))
            (str "x=" x " acc=" acc)))))

(t/deftest derivative-huge-and-special-arguments
  (let [sq (fn ^double [^double x] (* x x))]
    (t/testing "x = 1e8 (the step grows with |x| above 10); measured relative error 1.3e-11"
      (t/is (< (rel-err ((sut/derivative sq 1) 1.0e8) 2.0e8) 1.0e-6))
      (t/is (< (rel-err ((sut/derivative sq 1) -1.0e8) -2.0e8) 1.0e-6))
      (t/is (< (rel-err ((sut/derivative sq 1) 1.0e12) 2.0e12) 1.0e-6))
      (t/is (< (rel-err ((sut/derivative sq 2) 1.0e8) 2.0) 1.0e-3)))
    (t/testing "an explicit step is absolute: h below the spacing of doubles at a huge x gives zero"
      (t/is (= 0.0 ((sut/derivative sq 1 {:h 1.0e-9}) 1.0e8))))
    (t/testing "zero, negative argument"
      (t/is (close? ((sut/derivative sin-fn 1) 0.0) 1.0 1.0e-9))
      (t/is (close? ((sut/derivative sin-fn 2) 0.0) 0.0 1.0e-6))
      (t/is (close? ((sut/derivative (fn ^double [^double x] (* x x x)) 3) 0.0) 6.0 1.0e-6))
      (t/is (close? ((sut/derivative exp-fn 1) -30.0) (Math/exp -30.0) 1.0e-9)))
    (t/testing "NaN and infinite arguments give NaN"
      (t/is (m/invalid-double? ((sut/derivative sin-fn 1) ##NaN)))
      (t/is (m/invalid-double? ((sut/derivative sin-fn 1) ##Inf))))))

(t/deftest derivative-extrapolated
  (doseq [n [1 2 3]
          x [-1.3 0.0 0.5 2.0]]
    (t/is (close? ((sut/derivative sin-fn n {:acc 4 :extrapolate? true}) x) (exact-derivative :sin n x) 1.0e-6)
          (str "n=" n " x=" x)))
  (t/testing "options of extrapolate"
    (t/is (close? ((sut/derivative exp-fn 1 {:extrapolate? {:init-h 0.1 :contract 0.5}}) 1.0) (Math/exp 1.0) 1.0e-8))
    (t/is (invalid? #(sut/derivative exp-fn 1 {:extrapolate? {:contract 1.0}}) :contract 1.0))))

;; ------------------------------------------------------------------ gradient

(defn- grad-fn [[x y z]]
  (+ (* (Math/exp x) (Math/sin y)) (* x y z z) (/ 1.0 (+ 1.0 (* z z)))))

(defn- grad-exact [[x y z]]
  [(+ (* (Math/exp x) (Math/sin y)) (* y z z))
   (+ (* (Math/exp x) (Math/cos y)) (* x z z))
   (+ (* 2.0 x y z) (/ (* -2.0 z) (Math/pow (+ 1.0 (* z z)) 2.0)))])

(t/deftest gradient-analytic
  (t/testing "acc 1-8 with the automatic step; measured max error 6.3e-10 (acc 2, f of size 1e2)"
    (doseq [acc (range 1 9)
            point [[0.3 0.7 -1.2] [0.0 0.0 0.0] [-5.0 2.0 3.0] [1.0 1.0 1.0]]]
      (t/is (all-close? ((sut/gradient grad-fn {:acc acc}) point) (grad-exact point) 1.0e-9)
            (str "acc=" acc " " point))))
  (t/testing "default acc"
    (t/is (all-close? ((sut/gradient grad-fn) [0.3 0.7 -1.2]) (grad-exact [0.3 0.7 -1.2]) 1.0e-9)))
  (t/testing "explicit absolute step"
    (t/is (all-close? ((sut/gradient grad-fn {:h 1.0e-5}) [0.3 0.7 -1.2]) (grad-exact [0.3 0.7 -1.2]) 1.0e-8))
    (t/is (all-close? ((sut/gradient grad-fn {:h 1.0e-3 :acc 4}) [0.3 0.7 -1.2]) (grad-exact [0.3 0.7 -1.2]) 1.0e-9))
    (t/is (all-close? ((sut/gradient grad-fn {:h 1.0e-3 :acc 8}) [0.3 0.7 -1.2]) (grad-exact [0.3 0.7 -1.2]) 1.0e-11)))
  (t/testing "huge components get a scaled step"
    (let [g (sut/gradient (fn [[x]] (* x x)) {:acc 2})]
      (t/is (all-close? (g [1.0e8]) [2.0e8] 1.0e-6))))
  (t/testing "the odd accuracy equals the next even one"
    (t/is (= ((sut/gradient grad-fn {:acc 1}) [0.3 0.7 -1.2]) ((sut/gradient grad-fn {:acc 2}) [0.3 0.7 -1.2])))
    (t/is (= ((sut/gradient grad-fn {:acc 3}) [0.3 0.7 -1.2]) ((sut/gradient grad-fn {:acc 4}) [0.3 0.7 -1.2])))))

(t/deftest gradient-shape-and-arities
  (let [g (sut/gradient grad-fn)
        point [0.3 0.7 -1.2]
        res (g point)]
    (t/testing "returns a vector"
      (t/is (vector? res))
      (t/is (= 3 (count res)))
      (t/is (every? double? res)))
    (t/testing "any seqable input and a double array"
      (t/is (= res (g '(0.3 0.7 -1.2))))
      (t/is (= res (g (double-array point))))
      (t/is (= res (g (map identity point))))
      (t/is (every? double? (g [1 2 3]))))
    (t/testing "(g v h) uses an absolute step, 0 is the automatic one"
      (t/is (= ((sut/gradient grad-fn {:h 1.0e-4}) point) (g point 1.0e-4)))
      (t/is (= res (g point 0.0))))
    (t/testing "f gets a vector"
      (let [seen (atom nil)]
        ((sut/gradient (fn [v] (reset! seen v) 0.0)) (double-array [1.0 2.0]))
        (t/is (vector? @seen))))))

(t/deftest gradient-eager
  (t/testing "f is called when the gradient is called, not when the result is consumed"
    (let [calls (atom 0)
          g (sut/gradient (fn [[x y]] (swap! calls inc) (* x y)))]
      (t/is (zero? @calls))
      (g [1.0 2.0])
      (t/is (= 4 @calls))
      (reset! calls 0)
      ((sut/gradient (fn [[x y]] (swap! calls inc) (* x y)) {:acc 4}) [1.0 2.0])
      (t/is (= 8 @calls)))))

(t/deftest gradient-edge-cases
  (t/is (= [] ((sut/gradient (fn [_] 1.0)) [])))
  (t/is (= [] ((sut/gradient (fn [_] 1.0)) nil)))
  (t/testing "one variable equals derivative"
    (let [g ((sut/gradient (fn [[x]] (Math/sin x))) [0.4])]
      (t/is (close? (first g) (Math/cos 0.4) 1.0e-9)))))

;; ------------------------------------------------------------------ jacobian

(t/deftest jacobian-analytic
  (t/testing "R^2 -> R^3"
    (let [j ((sut/jacobian (fn [[x y]] [(* x y) (Math/sin x) (+ x y)])) [1.0 2.0])]
      ;; measured max error 2.9e-11
      (t/is (vector? j))
      (t/is (every? vector? j))
      (t/is (= [3 2] [(count j) (count (first j))]))
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-8) j [[2.0 1.0] [(Math/cos 1.0) 0.0] [1.0 1.0]])))))
  (t/testing "R^3 -> R^2 with acc 4"
    (let [j ((sut/jacobian (fn [[x y z]] [(+ (* x x) (* y z)) (* (Math/exp x) (Math/sin z))]) {:acc 4}) [0.3 -0.5 0.8])
          e (Math/exp 0.3)]
      ;; measured max error 6.7e-14
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-9)
                               j
                               [[0.6 0.8 -0.5] [(* e (Math/sin 0.8)) 0.0 (* e (Math/cos 0.8))]])))))
  (t/testing "polar to cartesian, f returns a seq / double array"
    (let [r 2.0 t 0.7
          expected [[(Math/cos t) (* (- r) (Math/sin t))] [(Math/sin t) (* r (Math/cos t))]]]
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-8)
                               ((sut/jacobian (fn [[r t]] (list (* r (Math/cos t)) (* r (Math/sin t))))) [r t])
                               expected)))
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-8)
                               ((sut/jacobian (fn [[r t]] (double-array [(* r (Math/cos t)) (* r (Math/sin t))]))) [r t])
                               expected)))))
  (t/testing "identity map has the identity Jacobian, a scalar-valued function in a vector is a one-row Jacobian"
    (t/is (every? true? (map #(all-close? %1 %2 1.0e-9)
                             ((sut/jacobian identity) [3.0 -4.0 5.0])
                             [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]])))
    (let [j ((sut/jacobian (fn [[x y]] [(* x y)])) [3.0 4.0])]
      (t/is (= 1 (count j)))
      (t/is (all-close? (first j) [4.0 3.0] 1.0e-8))))
  (t/testing "the Jacobian of a gradient-like map is the Hessian"
    (let [h ((sut/jacobian (fn [[x y]] [(* 2.0 x y) (+ (* x x) (* 3.0 y y))]) {:acc 4}) [1.5 -0.5])]
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-9) h [[-1.0 3.0] [3.0 -3.0]]))))))

(t/deftest jacobian-arities-and-edge-cases
  (let [f (fn [[x y]] [(* x y) (+ x y)])
        j (sut/jacobian f)]
    (t/is (= (j [1.0 2.0]) (j (double-array [1.0 2.0]))))
    (t/is (= (j [1.0 2.0]) (j '(1.0 2.0) 0.0)))
    (t/is (= ((sut/jacobian f {:h 1.0e-4}) [1.0 2.0]) (j [1.0 2.0] 1.0e-4))))
  (t/testing "empty point: m empty rows"
    (t/is (= [[] []] ((sut/jacobian (fn [_] [1.0 2.0])) []))))
  (t/testing "eager"
    (let [calls (atom 0)]
      ((sut/jacobian (fn [[x y]] (swap! calls inc) [x y])) [1.0 2.0])
      (t/is (= 4 @calls))))
  (t/testing "wrong shapes"
    (t/is (invalid-option? #((sut/jacobian (fn [_] 1.0)) [1.0 2.0]) :f))
    (t/is (invalid-option? #((sut/jacobian (fn [[x _]] (if (> x 1.0) [1.0] [1.0 2.0]))) [1.0 2.0]) :f))))

;; ------------------------------------------------------------------ hessian

(t/deftest hessian-automatic-step
  (t/testing "exp(xy) at (0.5, 1.0): relative error below 1e-6 (measured 6.9e-8)"
    (let [f (fn [[x y]] (Math/exp (* x y)))
          e (Math/exp 0.5)
          expected [[e (* 1.5 e)] [(* 1.5 e) (* 0.25 e)]]
          h ((sut/hessian f) [0.5 1.0])]
      (doseq [i (range 2)
              j (range 2)]
        (t/is (< (rel-err (get-in h [i j]) (get-in expected [i j])) 1.0e-6) (str i j)))))
  (t/testing "the step scales with the point: x^2 y at (1e4, 2e3)"
    (let [h ((sut/hessian (fn [[x y]] (* x x y))) [1.0e4 2.0e3])]
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-6) h [[4000.0 20000.0] [20000.0 0.0]])))))
  (t/testing "quadratic form, the Hessian is exact to rounding"
    (let [h ((sut/hessian (fn [[x y z]] (+ (* x x) (* 2.0 x y) (* 3.0 y y) (* y z) (* 4.0 z z)))) [0.3 -0.2 0.9])]
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-6) h [[2.0 2.0 0.0] [2.0 6.0 1.0] [0.0 1.0 8.0]]))))))

(t/deftest hessian-shape-and-arities
  (let [f (fn [[x y]] (+ (* x x x) (* x y y)))
        hs (sut/hessian f)
        h (hs [1.0 2.0])]
    (t/testing "vector of vectors, symmetric"
      (t/is (vector? h))
      (t/is (every? vector? h))
      (t/is (= h (apply mapv vector h)))
      (t/is (every? double? (flatten h))))
    (t/testing "values"
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-5) h [[6.0 4.0] [4.0 2.0]]))))
    (t/testing "arities and explicit step"
      (t/is (= h (hs (double-array [1.0 2.0]))))
      (t/is (= h (hs '(1.0 2.0) 0.0)))
      (t/is (= ((sut/hessian f {:h 1.0e-3}) [1.0 2.0]) (hs [1.0 2.0] 1.0e-3)))
      (t/is (every? true? (map #(all-close? %1 %2 1.0e-6) (hs [1.0 2.0] 1.0e-3) [[6.0 4.0] [4.0 2.0]]))))
    (t/testing "mixed term is divided by 4 h_i h_j, also for an explicit step"
      (t/is (close? (get-in (hs [1.0 2.0] 0.25) [0 1]) 4.0 1.0e-12)))))

(t/deftest hessian-edge-cases
  (t/is (= [] ((sut/hessian (fn [_] 1.0)) [])))
  (t/is (every? true? (map #(all-close? %1 %2 1.0e-6)
                           ((sut/hessian (fn [[x]] (Math/sin x))) [0.5])
                           [[(- (Math/sin 0.5))]])))
  (t/testing "f is called eagerly"
    (let [calls (atom 0)]
      ((sut/hessian (fn [[x y]] (swap! calls inc) (* x y))) [1.0 2.0])
      ;; f(v), 2 for each diagonal term and 4 for the mixed one
      (t/is (= 9 @calls)))))

(t/deftest hessian-diagonal-test
  (let [f (fn [[x y]] (+ (* x x x) (* 3.0 y y y y)))
        hd (sut/hessian-diagonal f)]
    (t/testing "returns double[]"
      (t/is (= (Class/forName "[D") (class (hd [1.0 2.0])))))
    (t/testing "automatic step"
      (t/is (all-close? (vec (hd [1.0 2.0])) [6.0 144.0] 1.0e-6))
      (t/is (all-close? (vec (hd (double-array [1.0 2.0]))) [6.0 144.0] 1.0e-6)))
    (t/testing "explicit step is absolute (the former behaviour used by acm): the second difference of a cubic is exact"
      (let [hd (sut/hessian-diagonal (fn [[x]] (* x x x)) {:h 0.1})]
        (t/is (close? (aget ^doubles (hd [1.0]) 0) 6.0 1.0e-10)))
      (t/is (= (vec ((sut/hessian-diagonal f {:h 1.0e-3}) [1.0 2.0])) (vec (hd [1.0 2.0] 1.0e-3))))
      (t/is (= (vec (hd [1.0 2.0])) (vec (hd [1.0 2.0] 0.0)))))
    (t/testing "agrees with the diagonal of the hessian"
      (t/is (all-close? (vec (hd [1.0 2.0])) (map-indexed (fn [i row] (row i)) ((sut/hessian f) [1.0 2.0])) 1.0e-9)))
    (t/testing "empty point"
      (t/is (= 0 (alength ^doubles ((sut/hessian-diagonal (fn [_] 0.0))  [])))))))

;; ------------------------------------------------------------------ laplacian

(t/deftest laplacian-analytic
  (t/testing "x^2 + y^2 + z^2 = 6 (measured error 8e-9)"
    (t/is (close? ((sut/laplacian (fn [[x y z]] (+ (* x x) (* y y) (* z z)))) [0.3 0.4 0.5]) 6.0 1.0e-7)))
  (t/testing "sin x sin y sin z = -3 f (measured error 6e-10)"
    (let [f (fn [[x y z]] (* (Math/sin x) (Math/sin y) (Math/sin z)))
          p [0.3 0.4 0.5]]
      (t/is (close? ((sut/laplacian f) p) (* -3.0 ^double (f p)) 1.0e-8))))
  (t/testing "exp(x + 2y) = 5 f (measured relative error 1.3e-9)"
    (t/is (< (rel-err ((sut/laplacian (fn [[x y]] (Math/exp (+ x (* 2.0 y))))) [0.3 0.2]) (* 5.0 (Math/exp 0.7))) 1.0e-8)))
  (t/testing "harmonic function x^2 - y^2 and the point far from the origin"
    (t/is (close? ((sut/laplacian (fn [[x y]] (- (* x x) (* y y)))) [3.0 -2.0]) 0.0 1.0e-7))
    (t/is (close? ((sut/laplacian (fn [[x y]] (+ (* x x) (* 2.0 y y)))) [1.0e5 -2.0e5]) 6.0 1.0e-3)))
  (t/testing "explicit step, arities"
    (let [l (sut/laplacian (fn [[x y]] (+ (* x x x) (* y y))))]
      (t/is (close? (l [1.0 1.0] 0.1) 8.0 1.0e-9))
      (t/is (= (l [1.0 1.0]) (l [1.0 1.0] 0.0)))
      (t/is (= (l [1.0 1.0] 1.0e-3) ((sut/laplacian (fn [[x y]] (+ (* x x x) (* y y))) {:h 1.0e-3}) [1.0 1.0])))
      (t/is (double? (l [1.0 1.0])))))
  (t/testing "empty point"
    (t/is (= 0.0 ((sut/laplacian (fn [_] 1.0)) [])))))

;; ------------------------------------------------------------------ extrapolate

(t/deftest extrapolate-behaviour
  (t/testing "forward difference extrapolated to h -> 0 (measured error 1.1e-13)"
    (let [g (fn ^double [^double x ^double h] (/ (- (Math/exp (+ x h)) (Math/exp x)) h))]
      (t/is (close? ((sut/extrapolate g) 0.5) (Math/exp 0.5) 1.0e-10))))
  (t/testing "limits at infinity: g(0, 1/h) for +Inf and g(0, -1/h) for -Inf"
    (let [g (fn ^double [^double _ ^double a] (+ 1.0 (/ 1.0 (* a a))))]
      (t/is (close? ((sut/extrapolate g {:power 2}) ##Inf) 1.0 1.0e-10))
      (t/is (close? ((sut/extrapolate g {:power 2}) ##-Inf) 1.0 1.0e-10)))
    (let [seen (atom [])
          g (fn ^double [^double x ^double a] (swap! seen conj [x a]) 1.0)]
      ((sut/extrapolate g {:max-evals 1 :init-h 0.25}) ##Inf)
      ((sut/extrapolate g {:max-evals 1 :init-h 0.25}) ##-Inf)
      (t/is (= [[0.0 4.0] [0.0 -4.0]] @seen))))
  (t/testing "NaN gives NaN"
    (t/is (m/invalid-double? ((sut/extrapolate (fn ^double [^double _ ^double _] 1.0)) ##NaN))))
  (t/testing ":max-evals 1 evaluates g once, at :init-h"
    (let [calls (atom 0)
          g (fn ^double [^double x ^double h] (swap! calls inc) (+ x h))]
      (t/is (= 1.25 ((sut/extrapolate g {:max-evals 1 :init-h 0.25}) 1.0)))
      (t/is (= 1 @calls))))
  (t/testing ":max-evals bounds the number of evaluations"
    (let [calls (atom 0)
          g (fn ^double [^double x ^double h] (swap! calls inc) (+ x (* 1.0e-3 (Math/sin (* 50.0 h)))))]
      ((sut/extrapolate g {:max-evals 5 :abs 0.0 :rel 0.0 :tol 1.0e10}) 1.0)
      (t/is (<= @calls 5))))
  (t/testing "negative :init-h is accepted"
    (let [g (fn ^double [^double x ^double h] (/ (- (Math/exp (+ x h)) (Math/exp x)) h))]
      (t/is (close? ((sut/extrapolate g {:init-h -0.5}) 0.5) (Math/exp 0.5) 1.0e-9)))))

;; ------------------------------------------------------------------ extrapolate against Richardson.jl

;; Reference values: Richardson.jl, `utils/fastmath/dev/generate_richardson_reference.jl`.
;; The stopping rule of the port tests the previous error estimate, so the results agree up to the tolerance
;; of the extrapolation: measured differences are at most 2.8e-13 in 15 of the 16 cases and 2.6e-9 in the noisy one.

(def ^:private richardson-functions
  {"central-sin" (fn ^double [^double x ^double h] (/ (- (Math/sin (+ x h)) (Math/sin (- x h))) (* 2.0 h)))
   "forward-exp" (fn ^double [^double x ^double h] (/ (- (Math/exp (+ x h)) (Math/exp x)) h))
   "sinc" (fn ^double [^double x ^double h] (/ (Math/sin (+ x h)) (+ x h)))
   "euler" (fn ^double [^double _ ^double h] (Math/pow (+ 1.0 h) (/ 1.0 h)))
   "central-atan" (fn ^double [^double x ^double h] (/ (- (Math/atan (+ x h)) (Math/atan (- x h))) (* 2.0 h)))
   "central-sqrt" (fn ^double [^double x ^double h] (/ (- (Math/sqrt (+ x h)) (Math/sqrt (- x h))) (* 2.0 h)))
   "noisy" (fn ^double [^double _ ^double h] (+ 1.0 (* h h) (* 1.0e-9 (Math/sin (/ 1.0 h)))))})

(def ^:private richardson-reference
  (delay (edn/read-string (slurp (io/resource "calculus/richardson_reference.edn")))))

(t/deftest extrapolate-reference
  (doseq [{:keys [id x init-h power abs rel max-evals value]} @richardson-reference]
    (let [result ((sut/extrapolate (richardson-functions id) {:contract 0.5 :power power :init-h init-h
                                                              :abs abs :rel rel :max-evals max-evals :tol 2.0}) x)]
      (t/is (m/delta-eq value result 1.0e-8)
            (str id " x " x " h " init-h " abs " abs " rel " rel " max-evals " max-evals)))))
