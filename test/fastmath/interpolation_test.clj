(ns fastmath.interpolation-test
  (:require [fastmath.interpolation :as sut]
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
