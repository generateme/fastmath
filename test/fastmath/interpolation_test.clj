(ns fastmath.interpolation-test
  (:require [fastmath.interpolation :as sut]
            [fastmath.kernel.variogram :as vg]
            [clojure.test :as t]
            [fastmath.core :as m]))

;; `extrapolation` wraps a 1D interpolator with a boundary policy controlling
;; behaviour outside [start, end]. Inside the domain (inclusive of the
;; boundaries) it must always delegate to the wrapped interpolator, whatever
;; the policy is.

(def ^:private xs [0.0 1.0 2.0 3.0])
(def ^:private ys [0.0 10.0 20.0 5.0])
(def ^:private lin (sut/linear xs ys))

;; a second, structurally different interpolator used to confirm `extrapolation`
;; is not specific to `linear`
(def ^:private cub (sut/cubic xs ys))

(t/deftest extrapolation-default-is-skip
  (t/testing "3-arg arity defaults to :skip - returns the interpolator unchanged"
    (let [e (sut/extrapolation lin 0.0 3.0)]
      (t/is (identical? lin e))
      (t/are [x] (m/delta-eq (lin x) (e x))
        -10.0 -1.0 0.0 1.5 3.0 4.0 10.0)))
  (t/testing "explicit :skip behaves the same as the default"
    (let [e (sut/extrapolation lin :skip 0.0 3.0)]
      (t/is (identical? lin e))
      (t/is (m/delta-eq (lin -1.0) (e -1.0)))
      (t/is (m/delta-eq (lin 4.0) (e 4.0))))))

(t/deftest extrapolation-inside-domain-always-delegates
  (t/testing "every policy delegates to the interpolator inside [start,end], boundaries included"
    (doseq [method [:skip :constant :zero :error 99.0 [-1.0 -2.0] {:left -1.0 :right -2.0}]]
      (let [e (sut/extrapolation lin method 0.0 3.0)]
        (t/are [x] (m/delta-eq (lin x) (e x))
          0.0 0.5 1.0 1.5 2.0 2.5 3.0)))))

(t/deftest extrapolation-constant
  (let [e (sut/extrapolation lin :constant 0.0 3.0)
        sv (lin 0.0)
        ev (lin 3.0)]
    (t/testing "outside domain clamps to the boundary value of the interpolator"
      (t/is (m/delta-eq sv (e -1.0)))
      (t/is (m/delta-eq sv (e -100.0)))
      (t/is (m/delta-eq ev (e 4.0)))
      (t/is (m/delta-eq ev (e 100.0))))
    (t/testing "boundaries and inside still delegate to the interpolator"
      (t/is (m/delta-eq sv (e 0.0)))
      (t/is (m/delta-eq ev (e 3.0)))
      (t/is (m/delta-eq (lin 1.5) (e 1.5))))))

(t/deftest extrapolation-zero
  (let [e (sut/extrapolation lin :zero 0.0 3.0)]
    (t/testing "outside domain returns exact zero"
      (t/are [x] (m/zero? (e x))
        -1.0 -100.0 4.0 100.0))
    (t/testing "boundaries and inside still delegate to the interpolator"
      (t/is (m/delta-eq (lin 0.0) (e 0.0)))
      (t/is (m/delta-eq (lin 3.0) (e 3.0)))
      (t/is (m/delta-eq (lin 1.5) (e 1.5))))))

(t/deftest extrapolation-number
  (t/testing "a double constant is returned outside the domain"
    (let [e (sut/extrapolation lin 42.0 0.0 3.0)]
      (t/are [x v] (m/delta-eq v (e x))
        -1.0 42.0
        4.0  42.0
        0.0  (lin 0.0)
        3.0  (lin 3.0)
        1.5  (lin 1.5))))
  (t/testing "a long constant is accepted and coerced to double"
    (let [e (sut/extrapolation lin 7 0.0 3.0)]
      (t/is (m/delta-eq 7.0 (e -1.0)))
      (t/is (m/delta-eq 7.0 (e 4.0))))))

(t/deftest extrapolation-left-right-vector
  (let [e (sut/extrapolation lin [-1.0 -2.0] 0.0 3.0)]
    (t/are [x v] (m/delta-eq v (e x))
      -1.0   -1.0
      -100.0 -1.0
      4.0    -2.0
      100.0  -2.0
      0.0    (lin 0.0)
      3.0    (lin 3.0)
      1.5    (lin 1.5))))

(t/deftest extrapolation-left-right-map
  (t/testing "the {:left :right} map form is equivalent to the two-element vector form"
    (let [ev (sut/extrapolation lin [-1.0 -2.0] 0.0 3.0)
          em (sut/extrapolation lin {:left -1.0 :right -2.0} 0.0 3.0)]
      (t/are [x] (m/delta-eq (ev x) (em x))
        -5.0 -1.0 0.0 1.5 3.0 4.0 10.0))))

(t/deftest extrapolation-error
  (let [e (sut/extrapolation lin :error 0.0 3.0)]
    (t/testing "throws outside [start,end]"
      (t/is (thrown? IndexOutOfBoundsException (e -0.0001)))
      (t/is (thrown? IndexOutOfBoundsException (e 3.0001)))
      (t/is (thrown? IndexOutOfBoundsException (e -100.0)))
      (t/is (thrown? IndexOutOfBoundsException (e 100.0))))
    (t/testing "does not throw at the boundaries or inside the domain"
      (t/is (m/delta-eq (lin 0.0) (e 0.0)))
      (t/is (m/delta-eq (lin 3.0) (e 3.0)))
      (t/is (m/delta-eq (lin 1.5) (e 1.5))))
    (t/testing "exception message reports x and the [start,end] range"
      (let [msg (try (e -1.0) (catch IndexOutOfBoundsException ex (.getMessage ex)))]
        (t/is (re-find #"-1\.0" msg))
        (t/is (re-find #"0\.0" msg))
        (t/is (re-find #"3\.0" msg))))))

(t/deftest extrapolation-works-with-any-1d-interpolator
  (t/testing "wrapping a different interpolator (cubic) behaves consistently"
    (let [e (sut/extrapolation cub :zero 0.0 3.0)]
      (t/is (m/zero? (e -1.0)))
      (t/is (m/zero? (e 4.0)))
      (t/is (m/delta-eq (cub 1.5) (e 1.5)))
      (t/is (m/delta-eq (cub 0.0) (e 0.0)))
      (t/is (m/delta-eq (cub 3.0) (e 3.0))))))

;; extrapolation: invalid configuration and NaN

(defn- property
  "Returns the :property of the ex-info thrown by `f`, nil when it does not throw one."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :extrapolation (:method d))) (:property d))))))

(t/deftest extrapolation-invalid-configuration
  (t/testing "unknown keyword, nil and values of other types throw"
    (t/is (= :method (property #(sut/extrapolation lin :constnat 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin nil 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin "zero" 0.0 3.0)))))
  (t/testing "sides: a vector needs two numbers, a map needs :left and :right"
    (t/is (= :method (property #(sut/extrapolation lin [1.0] 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin [] 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin [1.0 2.0 3.0] 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin [1.0 nil] 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin [1.0 "a"] 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin {:left 1.0} 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin {} 0.0 3.0))))
    (t/is (= :method (property #(sut/extrapolation lin {:right 1.0} 0.0 3.0)))))
  (t/testing "the domain has to be a non-empty interval"
    (t/is (= :domain (property #(sut/extrapolation lin :zero 3.0 0.0))))
    (t/is (= :domain (property #(sut/extrapolation lin :zero 1.0 1.0))))
    (t/is (= :domain (property #(sut/extrapolation lin :zero ##NaN 1.0))))
    (t/is (= :domain (property #(sut/extrapolation lin :zero 0.0 ##NaN))))
    (t/is (= :domain (property #(sut/extrapolation lin 0.0 3.0 0.0)))))
  (t/testing "valid configurations still work, also with infinite boundary"
    (t/is (fn? (sut/extrapolation lin :zero 0.0 3.0)))
    (t/is (= 0.0 ((sut/extrapolation lin :zero 0.0 ##Inf) -1.0)))
    (t/is (= 0.0 ((sut/extrapolation lin :zero ##-Inf 3.0) 4.0)))))

(t/deftest extrapolation-nan
  (t/testing "NaN gives NaN with every policy but :error"
    (doseq [method [:constant :zero 7.0 [1.0 2.0] {:left 1.0 :right 2.0}]]
      (t/is (Double/isNaN ((sut/extrapolation lin method 0.0 3.0) ##NaN)) (str method))))
  (t/testing ":error throws for NaN"
    (t/is (thrown? IndexOutOfBoundsException ((sut/extrapolation lin :error 0.0 3.0) ##NaN))))
  (t/testing ":skip passes NaN to the interpolator"
    (t/is (Double/isNaN ((sut/extrapolation lin :skip 0.0 3.0) ##NaN))))
  (t/testing "infinite arguments are outside the domain"
    (let [e (sut/extrapolation lin [1.0 2.0] 0.0 3.0)]
      (t/is (= 2.0 (e ##Inf)))
      (t/is (= 1.0 (e ##-Inf))))))

;; dispatch: the multimethod and the named function give the same interpolator

(def ^:private dxs (mapv double (range 15)))
(def ^:private dys (mapv #(+ (Math/sin (* 0.5 %)) (* 0.1 %)) dxs))
(def ^:private query (concat [-1.5 0.0 0.25 3.5 7.0 13.9 14.0 16.0] (range 0.3 14.0 1.7)))

(defn- same-values?
  [f g]
  (every? (fn [x] (m/delta-eq (double (f x)) (double (g x)) 1.0e-12)) query))

(t/deftest dispatch-1d
  (t/testing "every 1D keyword equals its named function"
    (doseq [[kw f] {:linear sut/linear :cubic sut/cubic :monotone sut/monotone :pchip sut/pchip
                    :akima sut/akima :neville sut/neville :divided-difference sut/divided-difference
                    :polynomial sut/polynomial :sprague sut/sprague :step-before sut/step-before
                    :step-after sut/step-after :step sut/step :b-spline sut/b-spline
                    :barycentric sut/barycentric :loess sut/loess :cubic-smoothing sut/cubic-smoothing
                    :pava sut/pava :cir sut/cir}]
      (t/is (same-values? (sut/interpolation kw dxs dys) (f dxs dys)) (str kw))))
  (t/testing "options reach the implementation"
    (t/is (same-values? (sut/interpolation :akima dxs dys {:method :makima}) (sut/akima dxs dys {:method :makima})))
    (t/is (not (same-values? (sut/interpolation :akima dxs dys {:method :makima}) (sut/akima dxs dys))))
    (t/is (same-values? (sut/interpolation :b-spline dxs dys {:degree 5}) (sut/b-spline dxs dys {:degree 5})))
    (t/is (not (same-values? (sut/interpolation :b-spline dxs dys {:degree 5}) (sut/b-spline dxs dys))))
    (t/is (same-values? (sut/interpolation :barycentric dxs dys {:order 3}) (sut/barycentric dxs dys {:order 3})))
    (t/is (same-values? (sut/interpolation :step dxs dys {:point 0.2}) (sut/step dxs dys {:point 0.2})))
    (t/is (same-values? (sut/interpolation :loess dxs dys {:bandwidth 0.6}) (sut/loess dxs dys {:bandwidth 0.6})))
    (t/is (same-values? (sut/interpolation :cubic-smoothing dxs dys {:rho 0.5}) (sut/cubic-smoothing dxs dys {:rho 0.5}))))
  (t/testing "the interpolators of interpolating methods pass through the nodes"
    (doseq [kw [:linear :cubic :monotone :pchip :akima :neville :polynomial :sprague :step-before :step-after :b-spline :barycentric]]
      (let [f (sut/interpolation kw dxs dys)]
        (t/is (every? (fn [[x y]] (m/delta-eq y (f x) 1.0e-9)) (map vector dxs dys)) (str kw)))))
  (t/testing "unknown keyword"
    (t/is (thrown? IllegalArgumentException (sut/interpolation :no-such-method dxs dys)))))

(t/deftest dispatch-isotonic-method
  (let [ys [1.0 3.0 2.0 4.0 3.5 6.0 7.0]
        xs [0.0 1.0 2.0 3.0 4.0 5.0 6.0]
        by-kw (sut/pava xs ys {:method :cubic})
        by-fn (sut/pava xs ys {:method sut/cubic})
        by-var (sut/pava xs ys {:method #'sut/cubic})
        by-default (sut/pava xs ys)
        by-linear (sut/pava xs ys {:method :linear})]
    (t/testing "a keyword, a function and a var select the method"
      (t/is (same-values? by-kw by-fn))
      (t/is (same-values? by-kw by-var))
      (t/is (same-values? by-default by-linear))
      (t/is (not (same-values? by-kw by-linear))))
    (t/testing "cir with function-like methods"
      (t/is (same-values? (sut/cir xs ys {:method :cubic}) (sut/cir xs ys {:method #'sut/cubic}))))
    (t/testing ":desc is the same as :decreasing"
      (t/is (same-values? (sut/pava xs (reverse ys) {:order :desc}) (sut/pava xs (reverse ys) {:order :decreasing})))
      (t/is (same-values? (sut/cir xs (reverse ys) {:order :desc}) (sut/cir xs (reverse ys) {:order :decreasing}))))))

(t/deftest dispatch-grid-and-scattered
  (let [gx [0.0 1.0 2.5 3.0 4.0]
        gy [0.0 0.5 1.5 2.0]
        vss (mapv (fn [x] (mapv (fn [y] (+ (* x x 0.1) (Math/sin y))) gy)) gx)
        points (map (fn [[x y]] [x y]) [[0.5 0.5] [3.5 1.0] [2.0 1.2] [-0.5 2.5]])]
    (doseq [[kw f] {:bilinear sut/bilinear :bicubic sut/bicubic :cubic-2d sut/cubic-2d}]
      (let [a (sut/interpolation kw gx gy vss)
            b (f gx gy vss)]
        (t/is (every? (fn [[x y]] (and (m/delta-eq (double (a x y)) (double (b x y)) 1.0e-12)
                                       (m/delta-eq (double (a [x y])) (double (b x y)) 1.0e-12))) points) (str kw))
        (t/is (every? (fn [[xi yj v]] (m/delta-eq v (double (a xi yj)) 1.0e-10))
                      (for [[i xi] (map-indexed vector gx) [j yj] (map-indexed vector gy)]
                        [xi yj (get-in vss [i j])])) (str kw " at the nodes")))))
  (let [xss (for [i (range 5) j (range 5)] [(* 0.5 i) (* 0.7 j)])
        ys (map (fn [[x y]] (+ (Math/sin x) (* 0.5 y))) xss)
        qs [[0.3 0.4] [1.1 1.9] [1.9 2.5]]
        variogram (vg/spherical {:nugget 0.0 :psill 1.0 :range 3.0})
        same? (fn [a b] (every? (fn [q] (m/delta-eq (double (a q)) (double (b q)) 1.0e-12)) qs))]
    (t/is (same? (sut/interpolation :shepard xss ys) (sut/shepard xss ys)))
    (t/is (same? (sut/interpolation :shepard xss ys {:p 3.0}) (sut/shepard xss ys {:p 3.0})))
    (t/is (same? (sut/interpolation :rbf xss ys) (sut/rbf xss ys)))
    (t/is (same? (sut/interpolation :kriging xss ys variogram) (sut/kriging xss ys variogram)))
    (t/is (same? (sut/interpolation :gp xss ys) (sut/gp xss ys)))
    (t/is (same? (sut/interpolation :microsphere-projection xss ys) (sut/microsphere-projection xss ys)))))
