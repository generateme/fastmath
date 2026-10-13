(ns fastmath.calculus.cubature-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.special :as special]
            [fastmath.vector :as v]
            [fastmath.calculus.cubature :as sut]))

;; Reference values: HCubature.jl, `utils/fastmath/dev/generate_cubature_reference.jl`.

(def ^:private integrands
  {"sincosexp3" (fn ^double [v] (m/* (m/sin (v 0)) (m/cos (v 1)) (m/exp (v 2))))
   "gauss2" (fn ^double [v] (m/exp (m/- (m/+ (m/sq (v 0)) (m/sq (v 1))))))
   "peak2" (fn ^double [v] (m// (m/+ 1.0e-2 (m/sq (m/- (v 0) 0.3)) (m/sq (m/- (v 1) 0.6)))))
   "abs3" (fn ^double [v] (m/abs (m/- (m/+ (m/sq (v 0)) (m/sq (v 1)) (m/sq (v 2))) 0.125)))
   "smalldim3" (fn ^double [v] (m/+ (m/* 1.0e-3 (m/sin (m/* 9.0 (v 2)))) 1.0))
   "anis4" (fn ^double [v] (m/* (m/exp (m/* 10.0 (v 0))) (m/inc (m/* 0.01 (v 1) (v 2) (v 3)))))
   "poly2" (fn ^double [v] (m/* (m/sq (v 0)) (v 1)))
   "sincos2" (fn ^double [v] (m/* (m/sin (v 0)) (m/cos (v 1))))
   "expsum5" (fn ^double [v] (m/exp (m/- (m/+ (m/sq (v 0)) (m/sq (v 1)) (m/sq (v 2)) (m/sq (v 3)) (m/sq (v 4))))))})

(def ^:private reference (delay (edn/read-string (slurp (io/resource "calculus/cubature_reference.edn")))))
(def ^:private axis-reference (delay (edn/read-string (slurp (io/resource "calculus/cubature_axis_reference.edn")))))

(defn- counted
  "Returns `[f' calls]` where `calls` is an atom counting the calls of `f'`."
  [f]
  (let [calls (atom 0)]
    [(fn ^double [v] (swap! calls inc) (f v)) calls]))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest reference-cases
  ;; measured: :result within 1e-15 relative, :error within a few eps * |I| (cancellation in I - I')
  (doseq [{:keys [id lower upper rel abs initdiv result error evaluations]} @reference]
    (let [[f calls] (counted (integrands id))
          res (sut/cubature f lower upper {:rel rel :abs abs :initdiv initdiv :max-iters 1000000 :info? true})]
      (t/testing (str id " " lower upper " rel " rel " initdiv " initdiv)
        (t/is (m/delta-eq result (:result res) 1.0e-12 1.0e-12) "result")
        (t/is (m/delta-eq error (:error res) (m/* 1.0e-13 (m/max 1.0 (m/abs result))) 1.0e-5) "error")
        (t/is (= evaluations (:evaluations res)) "evaluations of the reference")
        (t/is (= @calls (:evaluations res)) "evaluations are the calls of the integrand")
        (t/is (false? (:fail? res)))))))

(t/deftest split-axis
  (let [gm @#'sut/genz-malik
        integrate-gm @#'sut/integrate-gm
        fs {:sin (fn [k amp freq _] (fn [v] (m/* amp (m/sin (m/* freq (v k))))))
            :const (fn [_ amp _ _] (fn [_] amp))
            :sum (fn [_ amp freq _] (fn [v] (m/* amp (m/+ (m/sin (m/* freq (v 0))) (m/sin (m/* freq (v 1)))))))}]
    (doseq [{:keys [kind lower upper k amp freq axis]} @axis-reference]
      (let [n (count lower)
            {:keys [p w w']} (gm n)
            [gp1 gp2 gp3 gp4] p
            box (integrate-gm ((fs kind) k amp freq nil) lower upper n gp1 gp2 gp3 gp4 w w')]
        (t/is (= axis (.kdiv ^fastmath.calculus.cubature.CubatureBox box))
              (str (name kind) " n " n " k " k " amp " amp " freq " freq " box " upper))))))

(t/deftest limits-and-failure
  (let [hard (integrands "peak2")
        easy (integrands "poly2")
        gm-evals 17]
    (t/testing ":fail? is present in every branch"
      (t/is (= false (:fail? (sut/cubature easy [0 0] [1 1] {:initdiv 1 :info? true}))))
      (t/is (= false (:fail? (sut/cubature hard [0 0] [1 1] {:rel 1.0e-3 :info? true})))))
    (t/testing ":max-iters bounds the rounds including the initial one"
      (doseq [mi [1 2 3 10]]
        (let [res (sut/cubature hard [0 0] [1 1] {:rel 1.0e-12 :initdiv 1 :max-iters mi :info? true})]
          (t/is (= mi (:iterations res)))
          (t/is (= :max-iters (:fail? res)))
          (t/is (= (m/* gm-evals (m/inc (m/* 2 (m/dec mi)))) (:evaluations res))))))
    (t/testing "convergence in the last allowed round wins over the limit"
      (t/is (= false (:fail? (sut/cubature easy [0 0] [1 1] {:initdiv 1 :max-iters 1 :info? true})))))
    (t/testing ":max-evals is never exceeded, unless the initial round alone is above it"
      (doseq [me [100 500 2000]]
        (let [res (sut/cubature hard [0 0] [1 1] {:rel 1.0e-12 :initdiv 1 :max-evals me :info? true})]
          (t/is (<= (:evaluations res) me))
          (t/is (> (+ (:evaluations res) (* 2 gm-evals)) me) "no more room for another split")
          (t/is (= :max-evals (:fail? res)))))
      (let [res (sut/cubature hard [0 0] [1 1] {:rel 1.0e-12 :max-evals 10 :info? true})]
        (t/is (= 68 (:evaluations res)) "four initial boxes of 17 evaluations")
        (t/is (= :max-evals (:fail? res)))))
    (t/testing "without :info? a number is returned"
      (t/is (number? (sut/cubature hard [0 0] [1 1] {:max-iters 2}))))
    (t/testing "a NaN integrand never converges and is refined until the limit"
      (let [res (sut/cubature (fn [_] ##NaN) [0 0] [1 1] {:max-iters 5 :info? true})]
        (t/is (m/nan? (:result res)))
        (t/is (= :max-iters (:fail? res)))))))

(t/deftest tied-errors
  (t/testing "boxes with equal (zero) errors are all kept"
    (doseq [[k dims boxes] [[2 2 4] [3 2 9] [2 3 8]]]
      (let [res (sut/cubature (fn [_] 1.0) (repeat dims 0.0) (repeat dims 1.0) {:initdiv k :info? true})]
        (t/is (= boxes (:subdivisions res)))
        (t/is (m/delta-eq 1.0 (:result res) 1.0e-14)))))
  (t/testing "a symmetric integrand has equal errors in all the mirrored boxes"
    (t/is (m/delta-eq m/PI (sut/cubature (integrands "gauss2") [##-Inf ##-Inf] [##Inf ##Inf] {:rel 1.0e-9 :max-iters 10000}) 1.0e-8))))

(t/deftest dimensions
  (t/testing "one dimension is integrated by the Gauss-Kronrod engine; f still gets a vector"
    (let [seen (atom nil)
          res (sut/cubature (fn [v] (reset! seen v) (m/sq (v 0))) [0] [1] {:info? true})]
      (t/is (vector? @seen))
      (t/is (m/delta-eq (m// 1.0 3.0) (:result res) 1.0e-12))
      (t/is (= #{:result :error :evaluations :iterations :subdivisions :fail?} (set (keys res)))))
    (t/is (m/delta-eq 1.0 (sut/cubature (fn [v] (m/exp (m/- (v 0)))) [0] [##Inf] {:rel 1.0e-10}) 1.0e-9))
    (t/is (m/delta-eq (m/- (m/- 1.0 (m// m/E))) (sut/cubature (fn [v] (m/exp (m/- (v 0)))) [1] [0] {:rel 1.0e-10}) 1.0e-9) "reversed")
    (t/is (m/delta-eq m/PI (sut/cubature (fn [v] (m// (m/inc (m/sq (v 0))))) [##-Inf] [##Inf] {:rel 1.0e-10}) 1.0e-9))
    (t/testing "the options are passed on"
      (t/is (= 3 (:subdivisions (sut/cubature (fn [v] (v 0)) [0.0] [1.0] {:initdiv 3 :info? true}))))
      (t/is (= :max-iters (:fail? (sut/cubature (fn [v] (m/abs (m/sin (m/* 20.0 (v 0))))) [0.0] [1.0] {:rel 1.0e-14 :max-iters 3 :info? true}))))))
  (t/testing "empty bounds and different lengths throw ex-info"
    (doseq [[lower upper] [[[] []] [[0 0] [1]] [[0] [1 1]] [nil nil]]]
      (t/is (some? (ex-data-of #(sut/cubature (fn [_] 1.0) lower upper))))))
  (t/testing ":initdiv must be an integer >= 1"
    (doseq [bad [0 -1 1.5 "a" nil]]
      (t/is (= :initdiv (:option (ex-data-of #(sut/cubature (fn [_] 1.0) [0 0] [1 1] {:initdiv bad})))))
      (t/is (= :initdiv (:option (ex-data-of #(sut/cubature (fn [_] 1.0) [0] [1] {:initdiv bad})))))))
  (t/testing "equal bounds give zero"
    (t/is (= 0.0 (sut/cubature (fn [_] 1.0) [1 1] [1 1])))
    (t/is (= 0.0 (sut/cubature (fn [_] 1.0) [0 1] [1 1])))
    (t/is (= 0.0 (sut/cubature (fn [_] 1.0) [2] [2])))))

(t/deftest orientation
  (let [f (integrands "sincos2")
        straight (sut/cubature f [0.0 -1.0] [2.0 3.0] {:rel 1.0e-9})]
    (t/testing "every reversed dimension changes the sign"
      (t/is (m/delta-eq straight (sut/cubature f [0.0 -1.0] [2.0 3.0] {:rel 1.0e-9}) 1.0e-12))
      (t/is (m/delta-eq (m/- straight) (sut/cubature f [2.0 -1.0] [0.0 3.0] {:rel 1.0e-9}) 1.0e-8))
      (t/is (m/delta-eq (m/- straight) (sut/cubature f [0.0 3.0] [2.0 -1.0] {:rel 1.0e-9}) 1.0e-8))
      (t/is (m/delta-eq straight (sut/cubature f [2.0 3.0] [0.0 -1.0] {:rel 1.0e-9}) 1.0e-8)))
    (t/testing "the info map holds the oriented result and a positive error"
      (let [res (sut/cubature f [2.0 -1.0] [0.0 3.0] {:rel 1.0e-9 :info? true})]
        (t/is (m/delta-eq (m/- straight) (:result res) 1.0e-8))
        (t/is (m/pos? (:error res)))))
    (t/testing "infinite bounds in reverse order"
      (let [g (integrands "gauss2")]
        (t/is (m/delta-eq m/PI (sut/cubature g [##-Inf ##-Inf] [##Inf ##Inf] {:rel 1.0e-9 :max-iters 10000}) 1.0e-8))
        (t/is (m/delta-eq m/PI (sut/cubature g [##Inf ##Inf] [##-Inf ##-Inf] {:rel 1.0e-9 :max-iters 10000}) 1.0e-8))
        (t/is (m/delta-eq (m/- m/PI) (sut/cubature g [##Inf ##-Inf] [##-Inf ##Inf] {:rel 1.0e-9 :max-iters 10000}) 1.0e-8))))))

(t/deftest improper-integrals
  (t/testing "quarter plane and half plane of the standard bivariate Gauss"
    (let [g (integrands "gauss2")
          opts {:rel 1.0e-9 :max-iters 10000}]
      (t/is (m/delta-eq (m// m/PI 4.0) (sut/cubature g [0 0] [##Inf ##Inf] opts) 1.0e-8))
      (t/is (m/delta-eq (m// m/PI 2.0) (sut/cubature g [##-Inf 0] [##Inf ##Inf] opts) 1.0e-8))
      (t/is (m/delta-eq (m// m/PI 4.0) (sut/cubature g [##-Inf ##-Inf] [0 0] opts) 1.0e-8)))))

;; https://github.com/JohannesBuchner/cuba/blob/master/demo/demo-c.c
;; + wolfram alpha results

(def ^:private fun1-res 0.664669679781377268445214)
(defn- fun1 [[^double x ^double y ^double z]] (m/* (m/sin x) (m/cos y) (m/exp z)))

(def ^:private fun1+-res 3.019450507398802024611853)
(defn- fun1+ [[^double x ^double y ^double z]] (m/+ (m/sin x) (m/cos y) (m/exp z)))

(def ^:private fun2-res 5.26852)
(defn- fun2 [[^double x ^double y ^double z]] (m/* (m// (m/+ (m/sq (m/+ x y)) 0.003))
                                                (m/cos y) (m/exp z)))

(def ^:private fun3-res 0.307807)
(defn- fun3 [[^double x ^double y ^double z]] (m// (m/- 3.75 (m/cos (m/* m/PI x))
                                                     (m/cos (m/* m/PI y))
                                                     (m/cos (m/* m/PI z)))))

(def ^:private fun4-res 0.877314)
(defn- fun4 [v] (m/abs (m/- (v/magsq v) 0.125)))

(def ^:private fun5-res 0.416538385886638169609660243601001780010941391101)
(defn- fun5 [v] (m/exp (m/- (v/magsq v))))

(def ^:private fun6-res 1.202056902995100878845769442835673986997509412384)
(defn- fun6 [[^double x ^double y ^double z]] (m// 1.0 (m/- 1.0 (m/* x y z) -1.0e-10)))

(def ^:private fun7-res 0.709615885636934)
(defn- fun7 [[^double x ^double y ^double z]] (m/sqrt (m/abs (m/- x y z))))

(def ^:private fun8-res 0.891213)
(defn- fun8 [[^double x ^double y ^double z]] (m/exp (m/- (m/* x y z))))

(def ^:private fun9-res 0.0801856)
(defn- fun9 [[^double x ^double y ^double z]] (m// (m/sq x)
                                                (m/+ 5.0 (m/cos (m/+ 1.0 x y z)))))

(def ^:private fun11-res 0.5235987755982988)
(defn- fun11 [v] (if (m/< (v/magsq v) 1.0) 1.0 0.0))

(def ^:private normal-pdf (partial r/pdf (r/distribution :multi-normal {:means [0 0 0]
                                                                         :covariances [[1 0 0]
                                                                                       [0 0.5 0]
                                                                                       [0 0 2]]})))

(defn- gaussian [[^double _phi ^double r]] (m/* r (m/exp (m/* -0.5 r r))))

;; https://tutorial.math.lamar.edu/classes/calciii/TISphericalCoords.aspx
(defn- polar-3d
  [[^double ro ^double _theta ^double phi]]
  (m/* 16.0 ro ro ro (m/sin phi) (m/cos phi)))

(t/deftest multivariate-integrals
  (t/testing "unit cube; the default initial division is enough"
    (t/are [f res] (m/delta-eq res (sut/cubature f [0 0 0] [1 1 1] {:rel 1.0e-5 :abs 1.0e-5 :max-iters 20000}) 1.0e-3)
      fun1+ fun1+-res
      fun1 fun1-res
      fun2 fun2-res
      fun3 fun3-res
      fun4 fun4-res
      fun5 fun5-res
      fun6 fun6-res
      fun7 fun7-res
      fun8 fun8-res
      fun9 fun9-res))
  (t/testing "discontinuous integrand: the boxes along the surface converge slowly (measured error 2.1e-3 after 20000 rounds)"
    (t/is (m/delta-eq fun11-res (sut/cubature fun11 [0 0 0] [1 1 1] {:rel 1.0e-5 :abs 1.0e-5 :max-iters 20000}) 5.0e-3)))
  (t/testing "improper integrals"
    (t/are [lower upper res] (m/delta-eq res (sut/cubature normal-pdf lower upper {:rel 1.0e-6 :abs 1.0e-6 :max-iters 20000}) 1.0e-3)
      [##-Inf ##-Inf ##-Inf] [##Inf ##Inf ##Inf] 1.0
      [##-Inf ##-Inf ##-Inf] [0 0 0] 0.125
      [0 0 0] [##Inf ##Inf ##Inf] 0.125
      [##-Inf 0 ##-Inf] [0 ##Inf ##Inf] 0.25)
    (t/is (m/delta-eq m/TWO_PI (sut/cubature gaussian [0.0 0.0] [m/TWO_PI ##Inf] {:max-iters 20000}) 1.0e-3))
    (t/is (m/delta-eq (m/* 4.0 m/PI)
                      (sut/cubature polar-3d [0 0 0] [1.0 m/TWO_PI m/HALF_PI] {:max-iters 20000}) 1.0e-3))))

;; https://github.com/stevengj/cubature/blob/master/test.c

(defn- ftest0 ^double [v] (v/prod (v/cos v)))
(defn- ftest1 ^double [v] (let [[^double val ^double scale] (reduce (fn [[^double val ^double scale] ^double x]
                                                                   (if (m/pos? x)
                                                                     [(m/+ val (m/sq (m// (m/- 1.0 x) x)))
                                                                      (m/* scale (m// m/M_2_SQRTPI
                                                                                      (m/sq x)))]
                                                                     (reduced [0.0 0.0]))) [0.0 1.0] v)]
                         (m/* (m/exp (m/- val)) scale)))

(def ^:private ^{:const true :tag 'double} radius 0.50124145262344534123412)
(def ^:private ^{:const true :tag 'double} radius2 (m/sq radius))

(defn- ftest2 ^double [v] (if (m/< (v/magsq v) radius2) 1.0 0.0))
(defn- nball-V
  ^double [^long dims ^double radius]
  (let [hdims (m/* 0.5 dims)]
    (m/* (m/fpow radius dims)
         (m// (m/pow m/PI hdims)
              (special/gamma (m/inc hdims))))))

(def ^:private ^{:const true :tag 'double} fdata 0.1)

(defn- ftest3 ^double [v] (v/prod (v/mult v 2.0)))
(defn- ftest4 ^double [^long dim v] (let [sum (m/- (v/sum (v/sq (v/shift v -0.5))))]
                                   (m/* (m/fpow (m// m/M_2_SQRTPI (m/* 2.0 fdata)) dim)
                                        (m/exp (m// sum (m/sq fdata))))))
(defn- ftest7 ^double [^long dim v] (let [p (m// dim)
                                       prod (m/fpow (m/inc p) dim)]
                                   (m/* prod (v/prod (map (fn [^double x] (m/pow x p)) v)))))

(t/deftest stevengj-battery
  (let [mn (repeat 3 0.0)
        mx (repeat 3 1.0)
        dims (count mn)]
    (t/are [f res] (m/delta-eq res (sut/cubature f mn mx {:max-iters 20000 :rel 1.0e-3 :abs 1.0e-3}) 1.0e-3)
      ftest0 (v/prod (v/sin mx))
      ftest1 1.0
      ftest2 (m// (nball-V dims radius) (m/fpow 2 dims))
      ftest3 1.0
      (partial ftest4 dims) 1.0
      (partial ftest7 dims) 1.0)))

(t/deftest scipy-example
  ;; https://docs.scipy.org/doc/scipy/tutorial/integrate.html
  (t/is (v/delta-eq (map (fn [^double n]
                           (sut/cubature (fn [[^double t ^double x]] (m// (m/exp (m/* -1.0 x t)) (m/pow t n))) [1 0] [##Inf ##Inf]
                                         {:max-iters 5000}))
                         (range 2 10))
                    (map (fn [^double n] (m// n)) (range 2 10)) 1.0e-5)))
