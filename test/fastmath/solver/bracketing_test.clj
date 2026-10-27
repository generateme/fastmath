(ns fastmath.solver.bracketing-test
  "One battery which every bracketing solver `fastmath.solver.*` has to pass."
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.solver.brent :as brent]
            [fastmath.solver.bisection :as bisection]
            [fastmath.solver.regula-falsi :as regula-falsi]
            [fastmath.solver.illinois :as illinois]
            [fastmath.solver.pegasus :as pegasus]
            [fastmath.solver.ridders :as ridders]
            [fastmath.solver.muller :as muller]
            [fastmath.solver.muller2 :as muller2]
            [fastmath.solver.secant :as secant]
            [fastmath.solver.toms748 :as toms748]
            [fastmath.solver.newton :as newton]
            [fastmath.solver.halley :as halley])
  (:import [clojure.lang ExceptionInfo]))

(def solvers
  {:brent brent/brent
   :bisection bisection/bisection
   ;; the regula falsi converges only linearly, a convex function needs thousands of evaluations
   :regula-falsi (fn rf
                   ([f lo hi] (rf f lo hi nil))
                   ([f lo hi opts] (regula-falsi/regula-falsi f lo hi (merge {:max-iters 100000} opts))))
   :illinois illinois/illinois
   :pegasus pegasus/pegasus
   :ridders ridders/ridders
   :muller muller/muller
   :muller2 muller2/muller2
   :secant secant/secant
   :toms748 toms748/toms748
   ;; the safeguarded derivative methods with finite difference derivatives
   :newton newton/newton
   :halley halley/halley})

;; solvers which keep the correct root at huge scales (the others fail cleanly or are not tested there)
(def large-scale-solvers #{:brent :bisection})

;; cases a solver is not expected to solve: the regula falsi converges only linearly, and not at all in practice on a flat multiple root,
;; a pole or the extremely convex exp x - 1e10 (more than 1e5 evaluations)
(def excluded
  {:regula-falsi #{"exp x - 1e10" "(x-1)^3" "x^3 at zero" "pole"}
   ;; the methods without a bracket fail cleanly (:diverged), never with a wrong root: saturated or exploding values, poles and jumps
   :muller2 #{"exp x - 1e10" "1/x - 2, pole outside" "pole"}
   :secant #{"tanh steep" "(x-1)^2 (x-3)" "exp x - 1e10" "1/x - 2, pole outside" "pole" "jump" "plateau"}})

(defn- skip? [nm label] (boolean ((get excluded nm #{}) label)))

(defn- thrown-data [thunk] (try (thunk) ::no-throw (catch ExceptionInfo e (ex-data e))))

(defn- err [^double r ^double root] (m/abs (- r root)))

(defn- cube [x] (* x x x))

(def smooth
  "[name f lower upper root]"
  [["cos x - x" (fn [^double x] (- (m/cos x) x)) 0.0 1.0 0.7390851332151607]
   ["x^3 - 2" (fn [^double x] (- (cube x) 2.0)) 0.0 3.0 1.2599210498948732]
   ["tanh steep" (fn [^double x] (m/tanh (* 100.0 (- x 0.3)))) 0.0 1.0 0.3]
   ["(x-1)^2 (x-3)" (fn [^double x] (* (- x 1.0) (- x 1.0) (- x 3.0))) 2.0 5.0 3.0]
   ["x^3 + x" (fn [^double x] (+ (cube x) x)) -1.0 2.0 0.0]
   ["exp x - 1e10" (fn [^double x] (- (m/exp x) 1e10)) 0.0 100.0 23.025850929940457]
   ["1/x - 2, pole outside" (fn [^double x] (- (/ 1.0 x) 2.0)) 0.1 5.0 0.5]
   ["atan" (fn [^double x] (- (m/atan x) 1.0)) 0.0 10.0 (m/tan 1.0)]])

(t/deftest smooth-functions
  (doseq [[nm solve] solvers
          [fname f lo hi root] smooth
          :when (not (skip? nm fname))]
    (let [r (solve f lo hi {:absolute-accuracy 1e-12 :relative-accuracy 0.0})]
      (t/is (<= (err r root) (* 10.0 1e-12)) (str nm " " fname " got " r))
      (t/is (<= (err (solve f hi lo {:absolute-accuracy 1e-12 :relative-accuracy 0.0}) root) 1e-11) (str nm " " fname " reversed bounds")))
    (let [r (solve f lo hi)]
      (t/is (<= (err r root) (* 1e-13 (max 1.0 (m/abs root)))) (str nm " " fname " at the default tolerance, got " r)))))

(t/deftest root-at-an-end
  (doseq [[nm solve] solvers
          :let [f (fn [^double x] (- x 1.0))]]
    (t/is (== 1.0 (solve f 1.0 3.0)) (str nm " lower end"))
    (t/is (== 1.0 (solve f 0.0 1.0)) (str nm " upper end"))
    (t/is (<= (:evaluations (solve f 1.0 3.0 {:info? true})) 2) (str nm " costs at most two evaluations"))
    (t/is (== 0.0 (solve (fn [^double x] x) -1.0 0.0)) (str nm " zero as the upper end"))))

(t/deftest multiple-and-tiny-values
  (doseq [[nm solve] solvers]
    ;; an odd multiplicity: the value is flat near the root, the bracket decides
    (when-not (skip? nm "(x-1)^3")
      (let [r (solve (fn [^double x] (cube (- x 1.0))) 0.0 3.0 {:absolute-accuracy 1e-12 :relative-accuracy 0.0})]
        (t/is (<= (err r 1.0) 1e-11) (str nm " (x-1)^3, got " r))))
    ;; values below 1e-15 do not stop the search
    (let [r (solve (fn [^double x] (* 1e-16 (- x 0.3))) 0.0 1.0)]
      (t/is (<= (err r 0.3) 1e-14) (str nm " tiny values, got " r)))
    (let [r (solve (fn [^double x] (- x 1e-20)) -1.0 1.0)]
      (t/is (<= (err r 1e-20) 1e-33) (str nm " a root far below the interval scale, got " r)))
    ;; the default relative tolerance scales with the root
    (let [r (solve (fn [^double x] (- x 1.234e10)) 0.0 1e11)]
      (t/is (<= (err r 1.234e10) 1e-5) (str nm " a root of 1e10, got " r)))
    (let [r (solve (fn [^double x] (- x -1.234e-10)) -1.0 1.0)]
      (t/is (<= (err r -1.234e-10) 1e-25) (str nm " a small negative root, got " r)))
    ;; a jump: the sign change is the root
    (when-not (skip? nm "jump")
      (let [r (solve (fn [^double x] (m/signum (- x 0.3))) 0.0 1.0)]
        (t/is (<= (err r 0.3) 1e-12) (str nm " jump, got " r))))
    ;; the root at zero is found to the absolute tolerance
    (when-not (skip? nm "x^3 at zero")
      (let [r (solve (fn [^double x] (cube x)) -1.0 2.0 {:absolute-accuracy 1e-10})]
        (t/is (<= (m/abs r) 1e-9) (str nm " root at zero within 10 tolerances, got " r))))))

(t/deftest scale
  (doseq [[nm solve] solvers
          :when (large-scale-solvers nm)]
    (t/is (<= (err (solve (fn [^double x] (- x 1.2345e290)) 0.0 1e291) 1.2345e290) 1e276) (str nm " root 1e290"))
    (t/is (<= (err (solve (fn [^double x] (- x -1.2345e290)) -1e291 0.0) -1.2345e290) 1e276) (str nm " root -1e290"))
    (t/is (<= (err (solve (fn [^double x] (- x 1e-200)) 0.0 1.0) 1e-200) 1e-213) (str nm " root 1e-200"))
    (t/is (m/delta-eq 0.5 (solve (fn [^double x] (- x 0.5)) -1.7e308 1.7e308 {:max-iters 5000}) 1e-9) (str nm " the widest interval, no overflow"))))

(t/deftest failures
  (doseq [[nm solve] solvers]
    (let [d (thrown-data #(solve (fn [^double x] (+ 1.0 (* x x))) -1.0 2.0))]
      (t/is (= :no-bracket (:status d)) (str nm " no sign change"))
      (t/is (= 2 (:evaluations d))))
    (t/is (= :nan (:status (thrown-data #(solve (fn [^double x] (if (m/> x 0.7) ##NaN (- x 0.5))) 0.0 1.0)))) (str nm " NaN is reported"))
    (t/is (= :nan (:status (thrown-data #(solve (fn [^double _] ##NaN) 0.0 1.0)))) (str nm " NaN at the ends"))
    (t/is (= :invalid-result (:status (thrown-data #(solve (fn [^double _] nil) 0.0 1.0)))) (str nm " nil result"))
    (t/is (= :invalid-result (:status (thrown-data #(solve (fn [^double _] "a") 0.0 1.0)))) (str nm " string result"))
    (let [d (thrown-data #(solve (fn [^double x] (- (m/cos x) x)) 0.0 1.0 {:max-iters 3 :absolute-accuracy 0.0 :relative-accuracy 0.0}))]
      (t/is (= :max-evals (:status d)) (str nm " the limit of evaluations"))
      (t/is (= 3 (:evaluations d))))
    (let [info (solve (fn [^double x] (- (m/cos x) x)) 0.0 1.0 {:max-iters 3 :info? true})]
      (t/is (= :max-evals (:status info)) (str nm " :info? returns the status"))
      (t/is (m/nan? (:root info))))
    (t/is (= :lower-bound (:option (thrown-data #(solve (fn [^double x] x) ##-Inf 1.0)))) (str nm " infinite bound"))
    (t/is (= :max-iters (:option (thrown-data #(solve (fn [^double x] x) -1.0 1.0 {:max-iters 50.0})))) (str nm " a double limit is rejected"))
    (t/is (thrown? ExceptionInfo (solve (fn [^double _] (throw (ex-info "f failed" {}))) 0.0 1.0)) (str nm " an exception of f propagates"))
    (t/is (= "f failed" (try (solve (fn [^double _] (throw (ex-info "f failed" {}))) 0.0 1.0) (catch ExceptionInfo e (ex-message e)))) (str nm " unchanged"))))

(t/deftest infinite-values
  ;; an infinite value is a legal value of the function: the sign decides
  (doseq [[nm solve] solvers]
    (when-not (skip? nm "pole")
      (let [r (solve (fn [^double x] (/ 1.0 (- x 0.3))) 0.0 1.0 {:absolute-accuracy 1e-12 :relative-accuracy 0.0})]
        (t/is (<= (err r 0.3) 1e-11) (str nm " a pole is a sign change, got " r))))
    (when-not (skip? nm "plateau")
      (let [r (solve (fn [^double x] (if (m/< x 0.4) ##-Inf 1.0)) 0.0 1.0 {:absolute-accuracy 1e-12 :relative-accuracy 0.0})]
        (t/is (<= (err r 0.4) 1e-11) (str nm " infinite plateau, got " r))))))

(t/deftest tolerance
  (doseq [[nm solve] solvers]
    (let [info (solve (fn [^double x] (- x (/ 1.0 3.0))) 0.0 1.0 {:absolute-accuracy 1e-3 :relative-accuracy 0.0 :info? true})]
      (t/is (<= (err (:root info) (/ 1.0 3.0)) 1e-3) (str nm " a loose tolerance is respected"))
      (t/is (<= (:evaluations info) 20) (str nm " a loose tolerance saves evaluations")))
    (t/is (<= (err (solve (fn [^double x] (- x (/ 1.0 3.0))) 0.0 1.0 {:absolute-accuracy 0.0 :relative-accuracy 0.0}) (/ 1.0 3.0)) 1e-15)
          (str nm " zero tolerance: the resolution of doubles"))))

(t/deftest initial-value
  (doseq [[nm solve] solvers
          :let [f (fn [^double x] (- x 0.3))]]
    (t/is (<= (err (solve f 0.0 1.0 {:initial-value 0.31}) 0.3) 1e-14) (str nm " accepts :initial-value (used or ignored)"))
    (t/is (<= (err (solve f 0.0 1.0 {:initial-value 0.3}) 0.3) 1e-15) (str nm " the start is the root"))))

(t/deftest sequences
  ;; repeated calls are independent and deterministic
  (doseq [[nm solve] solvers
          :let [f (fn [^double x] (- (m/cos x) x))]]
    (t/is (== (solve f 0.0 1.0) (solve f 0.0 1.0)) (str nm " deterministic"))
    (t/is (= [(solve f 0.0 1.0)] (distinct (pmap (fn [_] (solve f 0.0 1.0)) (range 64)))) (str nm " thread safe"))))
