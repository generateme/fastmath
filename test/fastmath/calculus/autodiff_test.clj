(ns fastmath.calculus.autodiff-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.calculus.autodiff :as ad]
            [fastmath.calculus.autodiff.dual :as dual]
            [fastmath.calculus.autodiff.derivatives :as derivatives])
  (:import [fastmath.calculus.autodiff.dual Dual]))

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

(defn- tag ^long [^Dual d] (.tag d))
(defn- partials [^Dual d] (vec (.partials d)))

(defn- dual-is?
  "Is `d` a dual with the given tag, value and partials (structural equality, `clojure.core/=`)?"
  [d t v p]
  (and (instance? Dual d)
       (= d (ad/dual t v p))))

(defn- nan? [x] (and (number? x) (m/nan? (double x))))

(def ^:private D ad/derivative)

;; dual numbers

(t/deftest api-test
  (t/testing "the public API of the facade"
    (t/is (= #{"*" "+" "-" "/" "<" "<=" "==" ">" ">=" "cos" "derivative" "dual" "exp" "gradient" "hessian"
               "jacobian" "log" "pow" "primal" "scalar" "sin"}
             (set (map str (keys (ns-publics 'fastmath.calculus.autodiff))))))
    (t/is (every? #(:no-doc (meta (find-ns %)))
                  ['fastmath.calculus.autodiff.dual 'fastmath.calculus.autodiff.ops 'fastmath.calculus.autodiff.derivatives])
          "internal namespaces are excluded from the documentation"))
  (t/testing "invalid input throws (caller error, the exception type is not specified)"
    (t/is (thrown? Exception (ad/+ nil 1)))
    (t/is (thrown? Exception (ad/sin nil)))
    (t/is (thrown? Exception (ad/* (ad/dual 1 1 [1]) "a")))))

(t/deftest dual-test
  (t/testing "construction and coercion"
    (t/is (dual-is? (ad/dual 1 2 [1]) 1 2.0 [1.0]))
    (t/is (dual-is? (ad/dual 0 1/2 [3/4 1]) 0 0.5 [0.75 1.0]) "tag 0 is valid, ratios are coerced")
    (t/is (dual-is? (ad/dual 2 (ad/dual 1 1 [1]) [(ad/dual 1 0 [1])]) 2 (ad/dual 1 1.0 [1.0]) [(ad/dual 1 0.0 [1.0])])
          "nested duals with older tags")
    (t/is (= 3 (count (partials (ad/dual 1 0 [1 2 3])))) "more than one partial"))
  (t/testing "invalid arguments"
    (t/is (invalid? #(ad/dual -1 1 [1]) :tag -1))
    (t/is (invalid? #(ad/dual 1 1 []) :partials []))
    (t/is (invalid? #(ad/dual 1 (ad/dual 1 1 [1]) [1]) :tag 1) "equal nested tag")
    (t/is (invalid? #(ad/dual 1 (ad/dual 2 1 [1]) [1]) :tag 1) "newer nested tag in the value")
    (t/is (invalid? #(ad/dual 1 1 [(ad/dual 5 1 [1])]) :tag 1) "newer nested tag in the partials"))
  (t/testing "a dual is not a number"
    (let [d (ad/dual 1 2 [1])]
      (t/is (not (instance? Number d)))
      (t/is (thrown? ClassCastException (double d)))
      (t/is (= d d))))
  (t/testing "structural equality and hash (clojure.core/=)"
    (let [d (ad/dual 2 (ad/dual 1 3 [1]) [0 (ad/dual 1 -1 [2])])
          same (ad/dual 2 (ad/dual 1 3 [1]) [0 (ad/dual 1 -1 [2])])
          nan (ad/dual 1 ##NaN [##NaN])]
      (t/is (= d same) "same tag, value and partials")
      (t/is (= (hash d) (hash same)) "hasheq consistent with =")
      (t/is (= (.hashCode ^Object d) (.hashCode ^Object same)) "hashCode consistent with equals")
      (t/is (.equals ^Object d same))
      (t/is (not= d (ad/dual 3 (ad/dual 1 3 [1]) [0 (ad/dual 1 -1 [2])])) "different tag")
      (t/is (not= d (ad/dual 2 (ad/dual 1 3 [1]) [0 (ad/dual 1 -1 [3])])) "different nested partial")
      (t/is (not= (ad/dual 1 1 [1]) (ad/dual 1 1 [1 0])) "different number of partials")
      (t/is (not= (ad/dual 1 1 [0]) 1.0) "a dual is never = to a number")
      (t/is (not= 1.0 (ad/dual 1 1 [0])))
      (t/is (not= (ad/dual 1 1 [1]) nil))
      (t/is (= nan nan) "reflexive for the same object")
      (t/is (not= nan (ad/dual 1 ##NaN [##NaN])) "NaN leaves compare as clojure.core/=: a different NaN object is not equal")
      (t/is (= (ad/dual 1 0.0 [1 -0.0]) (ad/dual 1 -0.0 [1 0.0])) "0.0 equals -0.0, as clojure.core/=")
      (t/is (= (hash (ad/dual 1 0.0 [1 -0.0])) (hash (ad/dual 1 -0.0 [1 0.0]))) "and their hashes are equal")
      (t/is (= (hash d) (hash d) (.hashCode ^Object d)) "cached hash is stable, hashCode is hasheq")
      (t/is (= (hash d) (hash (ad/dual 2 (ad/dual 1 3 [1]) [0 (ad/dual 1 -1 [2])]))) "a fresh dual has the same hash")
      (t/is (= 1 (count (hash-set d same))) "set")
      (t/is (= :v (get {d :v} same)) "map key")
      (t/is (= :v (get {nan :v} nan)) "NaN dual as a map key, same object")))
  (t/testing "more than 32 partials"
    (let [p33 (range 33)
          d (ad/dual 1 0.5 p33)
          same (ad/dual 1 0.5 (map double p33))
          other (ad/dual 1 0.5 (concat (range 32) [99]))]
      (t/is (= d same))
      (t/is (= (hash d) (hash same)))
      (t/is (not= d other) "the last partial decides")
      (t/is (= 33 (count (partials d))))
      (t/is (= (str "#dual[1 0.5 " (pr-str (mapv double p33)) "]") (pr-str d)))
      (t/is (dual-is? (ad/* d 2) 1 1.0 (map #(* 2.0 %) p33)))
      (t/is (ad/< (ad/dual 1 0.5 (repeat 33 0)) d) "comparison over 33 partials")))
  (t/testing "primal and scalar"
    (let [d (ad/dual 2 (ad/dual 1 3 [1]) [1])]
      (t/is (instance? Dual (ad/primal d)))
      (t/is (= 3.0 (ad/primal (ad/primal d))))
      (t/is (= 3.0 (ad/scalar d)))
      (t/is (= 2.0 (ad/primal 2)) "a number is coerced")
      (t/is (= 0.5 (ad/scalar 1/2)))))
  (t/testing "printing"
    (t/is (= "#dual[2 #dual[1 3.0 [1.0]] [0.0]]" (pr-str (ad/dual 2 (ad/dual 1 3 [1]) [0])))))
  (t/testing "tags are unique across threads"
    (let [fresh-tag dual/fresh-tag
          tags (doall (pmap (fn [_] (fresh-tag)) (range 10000)))]
      (t/is (= 10000 (count (distinct tags)))))))

;; arithmetic

(t/deftest arithmetic-test
  (t/testing "arities and coercion of numbers"
    (t/is (= 0.0 (ad/+)))
    (t/is (= 1.0 (ad/*)))
    (t/is (= 3.0 (ad/+ 3)))
    (t/is (= 0.5 (ad/* 1/2)))
    (t/is (= -3.0 (ad/- 3)))
    (t/is (= 0.25 (ad// 4)))
    (t/is (= 10.0 (ad/+ 1 2 3 4)))
    (t/is (= 24.0 (ad/* 1 2 3 4)))
    (t/is (= -8.0 (ad/- 1 2 3 4)))
    (t/is (= 0.125 (ad// 1 2 4)))
    (t/is (= 3.5 (ad/+ 3 1/2)) "long and ratio")
    (t/is (= ##Inf (ad// 1 0)))
    (t/is (= ##-Inf (ad// -1 0)))
    (t/is (nan? (ad// 0 0))))
  (let [x (ad/dual 1 3 [1])
        y (ad/dual 1 2 [5])]
    (t/testing "dual and number, both orders"
      (t/is (dual-is? (ad/+ x 2) 1 5.0 [1.0]))
      (t/is (dual-is? (ad/+ 2 x) 1 5.0 [1.0]))
      (t/is (dual-is? (ad/- x 2) 1 1.0 [1.0]))
      (t/is (dual-is? (ad/- 2 x) 1 -1.0 [-1.0]))
      (t/is (dual-is? (ad/* x 2) 1 6.0 [2.0]))
      (t/is (dual-is? (ad/* 2 x) 1 6.0 [2.0]))
      (t/is (dual-is? (ad// x 2) 1 1.5 [0.5]))
      (t/is (dual-is? (ad// 6 x) 1 2.0 [(/ -6.0 9.0)]) "d(6/x) = -6/x^2"))
    (t/testing "two duals with the same tag"
      (t/is (dual-is? (ad/+ x y) 1 5.0 [6.0]))
      (t/is (dual-is? (ad/- x y) 1 1.0 [-4.0]))
      (t/is (dual-is? (ad/* x y) 1 6.0 [17.0]) "1*2 + 5*3")
      (t/is (dual-is? (ad// x y) 1 1.5 [(- 0.5 (* 5.0 0.75))]) "dx/y - dy*x/y^2")
      (t/is (dual-is? (ad/- x) 1 -3.0 [-1.0]))
      (t/is (dual-is? (ad// x) 1 (/ 3.0) [(/ -1.0 9.0)]))
      (t/is (dual-is? (ad/+ x y x) 1 8.0 [7.0]) "variadic")
      (t/is (dual-is? (ad/* x y 2) 1 12.0 [34.0])))
    (t/testing "duals with different tags: the newer tag is the outer layer"
      (let [z (ad/dual 2 10 [1])]
        (t/is (dual-is? (ad/+ x z) 2 (ad/dual 1 13.0 [1.0]) [1.0]))
        (t/is (dual-is? (ad/+ z x) 2 (ad/dual 1 13.0 [1.0]) [1.0]))
        (t/is (dual-is? (ad/* x z) 2 (ad/dual 1 30.0 [10.0]) [x]))
        (t/is (dual-is? (ad/* z x) 2 (ad/dual 1 30.0 [10.0]) [x])))))
  (t/testing "several partials"
    (let [x (ad/dual 1 2 [1 0])
          y (ad/dual 1 3 [0 1])]
      (t/is (dual-is? (ad/* x y) 1 6.0 [3.0 2.0]))))
  (t/testing "same tag, different numbers of partials: caller error"
    (let [x (ad/dual 1 1 [1 2])
          y (ad/dual 1 1 [1])
          mismatch? (fn [thunk] (= {:count-x 2 :count-y 1} (thrown-data thunk)))]
      (t/is (mismatch? #(ad/+ x y)))
      (t/is (mismatch? #(ad/- x y)))
      (t/is (mismatch? #(ad/* x y)))
      (t/is (mismatch? #(ad// x y)))
      (t/is (mismatch? #(ad/< x y)))
      (t/is (dual-is? (ad/+ x (ad/dual 2 1 [1])) 2 (ad/dual 1 2.0 [1.0 2.0]) [1.0]) "different tags may differ"))))

;; elementary functions

(t/deftest elementary-test
  (t/testing "numbers"
    (t/is (= (m/sin 0.5) (ad/sin 0.5)))
    (t/is (= (m/cos 0.5) (ad/cos 1/2)))
    (t/is (= 1.0 (ad/exp 0)))
    (t/is (= 0.0 (ad/log 1)))
    (t/is (nan? (ad/log -1)))
    (t/is (= ##-Inf (ad/log 0))))
  (let [x (ad/dual 1 0.5 [2])]
    (t/testing "chain rule"
      (t/is (dual-is? (ad/sin x) 1 (m/sin 0.5) [(* 2.0 (m/cos 0.5))]))
      (t/is (dual-is? (ad/cos x) 1 (m/cos 0.5) [(* 2.0 (- (m/sin 0.5)))]))
      (t/is (dual-is? (ad/exp x) 1 (m/exp 0.5) [(* 2.0 (m/exp 0.5))]))
      (t/is (dual-is? (ad/log x) 1 (m/log 0.5) [4.0]))))
  (t/testing "log at the domain boundary"
    (t/is (dual-is? (ad/log (ad/dual 1 0 [1])) 1 ##-Inf [##Inf]))
    (let [d (ad/log (ad/dual 1 -1 [1]))]
      (t/is (nan? (ad/primal d)))
      (t/is (every? nan? (partials d)) "outside the domain the derivative is NaN, as the value")))
  (t/testing "log outside the domain: NaN derivatives of every order"
    (doseq [x [-1.0 -2.0 -1.0e308 ##-Inf]]
      (t/is (nan? ((D ad/log) x)) (str x))
      (t/is (nan? ((D ad/log 2) x)) (str x " second order")))
    (t/is (nan? ((D #(ad/log (ad/log %))) 0.5)) "log of a negative intermediate value")
    (t/is (nan? ((D #(ad/log (ad/log %)) 2) 0.5)))
    (t/is (= ##-Inf ((D ad/log) -0.0)) "-0.0 is on the domain boundary, as ForwardDiff")
    (t/is (= 1.0 ((D ad/log) 1.0)) "inside the domain unchanged")))

;; pow

(defn- d-pow
  "Derivative of `f` at `x` and the value."
  [f x]
  [((D f) x) (f x)])

(t/deftest pow-test
  (t/testing "numbers"
    (t/is (= 8.0 (ad/pow 2 3)))
    (t/is (nan? (ad/pow -2 0.5))))
  (t/testing "variable base, constant exponent: y x^(y-1)"
    (t/is (= [0.0 0.0] (d-pow #(ad/pow % 2) 0)))
    (t/is (= [0.0 1.0] (d-pow #(ad/pow % 0) 0)) "y = 0 at x = 0")
    (t/is (= [0.0 1.0] (d-pow #(ad/pow % 0) -3)) "y = 0 at a negative x")
    (t/is (= [1.0 0.0] (d-pow #(ad/pow % 1) 0)))
    (t/is (= [##Inf 0.0] (d-pow #(ad/pow % 0.5) 0)))
    (t/is (= [0.0 0.0] (d-pow #(ad/pow % 1.5) 0)))
    (t/is (= [12.0 -8.0] (d-pow #(ad/pow % 3) -2)) "negative base, integer exponent")
    (t/is (every? nan? (d-pow #(ad/pow % 2.5) -2)) "negative base, fractional exponent")
    (t/is (= 2.0 ((D #(ad/pow % 2) 2) 0)) "second derivative of x^2 at 0"))
  (t/testing "constant base, variable exponent: x^y log(x)"
    (t/is (= [(* 2.0 (m/log 2.0)) 2.0] (d-pow #(ad/pow 2 %) 1)))
    (t/is (= [0.0 0.0] (d-pow #(ad/pow 0 %) 1)) "x = 0, y > 0")
    (t/is (= [##-Inf 1.0] (d-pow #(ad/pow 0 %) 0)) "x = 0, y = 0")
    (t/is (= [##-Inf ##Inf] (d-pow #(ad/pow 0 %) -1)) "x = 0, y < 0")
    (let [[d v] (d-pow #(ad/pow -2 %) 2)]
      (t/is (nan? d) "negative base")
      (t/is (= 4.0 v))))
  (t/testing "variable base and exponent"
    (t/is (= [1.0 1.0] (d-pow #(ad/pow % %) 1)))
    (t/is (every? nan? (d-pow #(ad/pow % %) -0.5)))
    (let [[d v] (d-pow #(ad/pow % %) 0)]
      (t/is (nan? d) "0^0, as in ForwardDiff")
      (t/is (= 1.0 v)))
    (let [[d v] (d-pow #(ad/pow % %) -1)]
      (t/is (nan? d) "the derivative with respect to the exponent needs log(-1)")
      (t/is (= -1.0 v)))
    (let [e9log3 (* 9.0 (m/log 3.0))]
      (t/is (= [(* 2.0 e9log3) 9.0] (d-pow #(ad/pow 3 (ad/* 2 %)) 1)) "3^(2x)")))
  (t/testing "a dual with zero partials is a constant"
    (let [x (ad/dual 1 -2 [1])
          y (ad/dual 1 3 [0])
          r (ad/pow x y)]
      (t/is (dual-is? r 1 -8.0 [12.0]) "negative base, exponent with zero partials, no log(-2)"))
    (let [x (ad/dual 1 0 [0])
          y (ad/dual 1 0.5 [1])
          r (ad/pow x y)]
      (t/is (dual-is? r 1 0.0 [0.0]) "zero base with zero partials, no 0 * Inf"))
    (let [x (ad/dual 1 0 [0])
          r (ad/pow x 0.5)]
      (t/is (dual-is? r 1 0.0 [0.0]) "x^0.5 at 0 with zero partials is not NaN"))
    (let [x (ad/dual 1 -2 [1])
          y (ad/dual 1 0 [0])]
      (t/is (dual-is? (ad/pow x y) 1 1.0 [0.0]) "exponent zero with zero partials"))))

;; comparisons

(t/deftest comparison-test
  (t/testing "numbers"
    (t/is (ad/< 1 2))
    (t/is (not (ad/< 2 2)))
    (t/is (ad/<= 2 2))
    (t/is (ad/> 3 2))
    (t/is (ad/>= 2 2))
    (t/is (ad/== 2 2.0))
    (t/is (ad/== 0.0 -0.0))
    (t/is (ad/< 1) "one argument")
    (t/is (ad/== ##NaN) "one argument, as clojure.core")
    (t/is (ad/< 1 2 3 4))
    (t/is (not (ad/< 1 2 2 4)))
    (t/is (ad/<= 1 2 2 4))
    (t/is (ad/== 1 1 1))
    (t/is (not (ad/== 1 1 2)))
    (t/is (ad/> 4 3 1/2)))
  (t/testing "NaN is unordered"
    (t/is (not-any? (fn [f] (f ##NaN 1.0)) [ad/< ad/<= ad/> ad/>= ad/==]))
    (t/is (not-any? (fn [f] (f (ad/dual 1 ##NaN [1]) 1.0)) [ad/< ad/<= ad/> ad/>= ad/==])))
  (t/testing "lexicographic on value and partials"
    (let [x (ad/dual 1 1 [1])]
      (t/is (not (ad/== x 1.0)) "a dual with nonzero partials is not equal to its value")
      (t/is (ad/== (ad/dual 1 1 [0]) 1.0) "zero partials")
      (t/is (ad/> x 1.0))
      (t/is (ad/< 1.0 x))
      (t/is (ad/< x 1.5) "values decide first")
      (t/is (ad/< (ad/dual 1 1 [1]) (ad/dual 1 1 [2])))
      (t/is (ad/<= (ad/dual 1 1 [1 2]) (ad/dual 1 1 [1 2])))
      (t/is (ad/< (ad/dual 1 1 [1 2]) (ad/dual 1 1 [1 3])) "second partial decides")
      (t/is (ad/< (ad/dual 1 1 [-1]) 1.0))
      (t/is (ad/> (ad/dual 2 1 [1]) (ad/dual 1 1 [0])) "the older dual is a constant of the newer layer")
      (t/is (ad/< (ad/dual 2 1 [1]) (ad/dual 1 1 [5])) "values first: 1 < 1 + 5e (older perturbation)")))
  (t/testing "branches follow the perturbation"
    (let [abs (fn [x] (if (ad/< x 0) (ad/- x) x))]
      (t/is (= 1.0 ((D abs) 0)))
      (t/is (= 1.0 ((D #(abs (ad/- %))) 0)) "one-sided derivative in the direction of the perturbation, as ForwardDiff")
      (t/is (= -1.0 ((D abs) -1))))
    (let [sinc (fn [x] (if (ad/== x 0) 1.0 (ad// (ad/sin x) x)))
          sinc-scalar (fn [x] (if (zero? (ad/scalar x)) 1.0 (ad// (ad/sin x) x)))]
      (t/is (nan? ((D sinc) 0)) "the dual at 0 is not equal to 0, sin(0)/0 is evaluated (ForwardDiff gives NaN too)")
      (t/is (= 0.0 ((D sinc-scalar) 0)) "branch on the value: the derivative of the constant branch"))))

;; derivative

(def ^:private reference-functions
  "Functions of the reference file, written as in utils/fastmath/dev/generate_autodiff_reference.jl."
  {:sin ad/sin
   :cos ad/cos
   :exp ad/exp
   :exp-sin (fn [x] (ad/* (ad/exp x) (ad/sin x)))
   :poly (fn [x] (ad/+ (ad/* 3.0 (ad/pow x 4.0)) (ad/* -2.0 (ad/pow x 3.0)) x -5.0))
   :rational (fn [x] (ad// (ad/+ (ad/pow x 2.0) 1.0) (ad/+ (ad/cos x) 2.0)))
   :log ad/log
   :log-pow (fn [x] (ad/log (ad/+ (ad/pow x 2.5) 1.0)))
   :const-pow (fn [x] (ad/pow 2.0 x))
   :pow-pow (fn [x] (ad/pow x x))
   :pow-sin (fn [x] (ad/pow (ad/+ (ad/sin x) 2.0) (ad/cos x)))
   :reciprocal (fn [x] (ad/- (ad// x) (ad// x (ad/+ (ad/pow x 2.0) 3.0))))})

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "calculus/autodiff_reference.edn")))))

(t/deftest derivative-test
  (t/testing "order"
    (t/is (identical? ad/sin (D ad/sin 0)) "n = 0 is the function itself")
    (t/is (= (m/cos 1.0) ((D ad/sin 1) 1.0)))
    (t/is (= ((D ad/sin) 1.0) ((D ad/sin 1) 1.0)) "default order is 1")
    (t/is (= (m/cos 1.0) ((apply D [ad/sin 1.7]) 1.0)) "fractional order is truncated")
    (t/is (invalid? #(D ad/sin -1) :n -1))
    (t/is (= (m/sin 1.0) ((D ad/sin 4) 1.0)) "fourth order"))
  (t/testing "points"
    (t/is (= (m/cos 1.0) ((D ad/sin) 1)) "long point")
    (t/is (= (m/cos 0.5) ((D ad/sin) 1/2)) "ratio point"))
  (t/testing "constant functions"
    (t/is (= 0.0 ((D (constantly 3)) 1.0)) "long result")
    (t/is (= 0.0 ((D (constantly 3.0)) 1.0)))
    (t/is (= 0.0 ((D (fn [x] (ad/- x x))) 1.0)) "x - x has the zero derivative")
    (t/is (= 0.0 ((D (constantly 3.0) 3) 1.0)) "higher order"))
  (t/testing "bulk comparison with ForwardDiff"
    (doseq [[k f] reference-functions
            :let [{:keys [x] :as r} (@reference k)]
            n [0 1 2 3]
            :let [d (D f n)]
            [xx expected] (map vector x (r (keyword (str "d" n))))]
      ;; measured: max relative error 4.4e-16, max absolute error 1.1e-13 (reciprocal, n=3, value 740)
      (t/is (m/delta-eq expected (double (d xx)) 1.0e-14 1.0e-14) (str k " n=" n " x=" xx)))))

;; perturbation confusion

(t/deftest perturbation-confusion-test
  (t/testing "Siskind and Pearlmutter 2005: d/dx (x * d/dy (x + y)) at 1 is 1, not 2"
    (t/is (= 1.0 ((D (fn [x] (ad/* x ((D (fn [y] (ad/+ x y))) 1)))) 1)))
    (let [constant-one (fn [x] ((D (fn [y] (ad/+ x y))) 1))]
      (t/is (= ((D (fn [x] (ad/* x (constant-one x)))) 1)
               ((D (fn [x] (ad/* x 1))) 1))
            "referential transparency")))
  (t/testing "inner closure over the outer variable"
    (let [g (fn [x] ((D (fn [y] (ad/* x y))) 2.0))]
      (t/is (= 6.0 ((D (fn [x] (ad/* x (g x)))) 3.0)) "d/dx x^2 = 2x")))
  (t/testing "one derivative function nested in itself"
    (let [dsin (D ad/sin)]
      (t/is (= (- (m/sin 1.0)) ((D dsin) 1.0)))
      (t/is (= (- (m/cos 1.0)) ((D (D dsin)) 1.0)))
      (t/is (m/delta-eq (* (m/sin (m/cos 1.0)) (m/sin 1.0)) ((D (fn [x] (dsin (dsin x)))) 1.0) 1.0e-15)
            "d/dx cos(cos(x)) = sin(cos(x)) sin(x)")))
  (t/testing "mixed nesting: d/dx (x * d/dy (x * y^2)) = d/dx (x * 2xy) at y = 3"
    ;; x * 2 * x * 3 = 6x^2, derivative 12x
    (t/is (= 24.0 ((D (fn [x] (ad/* x ((D (fn [y] (ad/* x (ad/pow y 2)))) 3)))) 2.0))))
  (t/testing "derivatives in parallel threads"
    (let [d2 (D (fn [x] (ad/* (ad/sin x) (ad/exp x))) 2)
          xs (range -2.0 2.0 0.01)]
      (t/is (= (map d2 xs) (doall (pmap d2 xs)))))))

;; errors

(t/deftest error-test
  (t/testing "a function returning a function"
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"returning a function"
                            ((D (fn [x] (fn [y] (ad/* x y)))) 1.0))))
  (t/testing "a function returning a vector"
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"neither a number nor a dual"
                            ((D (fn [x] [x x])) 1.0))))
  (t/testing "a dual escaping through an atom"
    (let [a (atom nil)
          g (D (fn [y] (reset! a y) y))
          data (thrown-data #((D (fn [x] (g x) (ad/* x @a))) 1.0))]
      (t/is (some? data))
      (t/is (< (:tag data) (:result-tag data)))))
  (t/testing "a dual escaping through a lazy sequence realized later"
    (let [ys (atom nil)
          g (D (fn [y] (reset! ys (map (fn [k] (ad/* y k)) [1 2 3])) y))]
      (t/is (some? (thrown-data #((D (fn [x] (g x) (ad/+ x (first @ys)))) 1.0))))))
  (t/testing "limitation: a stale dual of a finished, unrelated derivative call is not detected"
    (let [leak (atom nil)]
      ((D (fn [y] (reset! leak y) y)) 2.0)
      ;; the stale dual is older than the new tag, so it is treated as a constant
      (t/is (= 1.0 ((D (fn [x] (ad/+ x @leak))) 1.0)))))
  (t/testing "limitation: an escaped dual used only in a comparison is not detected"
    (let [leak (atom nil)
          g (D (fn [y] (reset! leak y) y))]
      ;; the correct derivative is 1 (x < x is false), the escaped dual takes the other branch
      (t/is (= 10.0 ((D (fn [x] (g x) (if (ad/< x @leak) (ad/* 10 x) x))) 1.0))))))

;; gradient, jacobian, hessian

(def ^:private multivariate-reference
  (delay (edn/read-string (slurp (io/resource "calculus/autodiff_multivariate_reference.edn")))))

(def ^:private scalar-functions
  "Functions of the multivariate reference file, written as in utils/fastmath/dev/generate_autodiff_multivariate_reference.jl."
  {:rosen (fn [v] (reduce ad/+ 0.0 (map (fn [a b] (ad/+ (ad/* 100.0 (ad/pow (ad/- b (ad/pow a 2.0)) 2.0))
                                                         (ad/pow (ad/- 1.0 a) 2.0))) v (rest v))))
   :dense (fn [v] (let [s (reduce ad/+ 0.0 v)
                        q (reduce ad/+ 0.0 (map (fn [a] (ad/* a a)) v))]
                    (reduce ad/+ 0.0 (map (fn [a] (ad/* (ad/sin (ad/* a s)) (ad/exp (ad// q 100.0)))) v))))})

(def ^:private vector-functions
  {:vec32 (fn [[x y z]] [(ad/* x y (ad/sin z)) (ad// (ad/exp x) (ad/+ 1.0 (ad/pow y 2.0) (ad/pow z 2.0)))])
   :vec23 (fn [[x y]] [(ad// x y) (ad/log (ad/+ x (ad/pow y 2.0))) (ad/pow x y)])})

(defn- close-to-reference?
  "Elementwise comparison of nested vectors with the reference."
  [expected actual]
  (if (sequential? expected)
    (and (sequential? actual)
         (= (count expected) (count actual))
         (every? true? (map close-to-reference? expected actual)))
    ;; measured: max error relative to max(1, |expected|) 8.8e-16 (gradient, hessian), 2.2e-16 (jacobian)
    (m/delta-eq (double expected) (double actual) 1.0e-14 1.0e-14)))

(defn- matrix-max-abs ^double [h] (reduce m/max 0.0 (map m/abs (flatten h))))

(t/deftest chunk-policy-test
  (t/testing "default chunk size: balanced chunks of at most 8"
    (t/is (= [0 1 8 5 8 6 8 8] (map @#'derivatives/default-chunk-size [0 1 8 9 16 17 50 100]))))
  (t/testing "options"
    (doseq [api [ad/gradient ad/jacobian ad/hessian]]
      (t/is (invalid? #(api identity {:chunk-size 0}) :chunk-size 0))
      (t/is (invalid? #(api identity {:chunk-size -1}) :chunk-size -1))
      (t/is (fn? (api identity {:chunk-size 1})) "smallest valid chunk size")
      (t/is (invalid? #(api identity {:chunk-size ##Inf}) :chunk-size ##Inf))
      (t/is (invalid? #(api identity {:chunk-size ##-Inf}) :chunk-size ##-Inf))
      (t/is (let [data (thrown-data #(api identity {:chunk-size ##NaN}))]
              (and (= :chunk-size (:option data)) (nan? (:value data)))) "NaN")
      (t/is (invalid? #(api identity {:chunk-size 0.5}) :chunk-size 0.5))
      (t/is (thrown? ClassCastException (api identity {:chunk-size false})) "false is not a missing option")
      (t/is (thrown? ClassCastException (api identity {:chunk-size "2"}))))
    (let [f (:rosen scalar-functions)
          x [0.1 0.2 0.3]
          r ((ad/gradient f) x)]
      (t/is (= r ((ad/gradient f {}) x) ((ad/gradient f nil) x) ((ad/gradient f {:chunk-size nil}) x))
            "no options, empty options and nil options are the same")
      (t/is (= r ((ad/gradient f {:chunk-size 100}) x)) "chunk size larger than the number of inputs is clamped")
      (t/is (= r ((ad/gradient f {:chunk-size 1.0e30}) x) ((ad/gradient f {:chunk-size Long/MAX_VALUE}) x))
            "chunk sizes beyond the long range are clamped")
      (t/is (= r ((ad/gradient f {:chunk-size 1.9}) x) ((ad/gradient f {:chunk-size 3/2}) x)) "fractional chunk sizes are truncated")
      (t/is (= ((ad/hessian f) x) ((ad/hessian f {}) x) ((ad/hessian f nil) x)))
      (t/is (= ((ad/jacobian (:vec32 vector-functions)) x) ((ad/jacobian (:vec32 vector-functions) nil) x)))))
  (t/testing "passes: the chunk inputs are duals of one fresh tag, the others stay numbers"
    (let [calls (atom [])
          f (fn [v] (swap! calls conj v) (reduce ad/+ v))
          x [1.0 2.0 3.0 4.0 5.0]]
      ((ad/gradient f {:chunk-size 2}) x)
      (let [passes @calls
            dual-positions (map (fn [v] (keep-indexed (fn [i e] (when (instance? Dual e) i)) v)) passes)
            tags (map (fn [v] (distinct (map tag (filter #(instance? Dual %) v)))) passes)]
        (t/is (= [[0 1] [2 3] [4]] dual-positions) "chunks of 2, the last one shorter")
        (t/is (every? #(= 1 (count %)) tags) "one tag per pass")
        (t/is (apply < (map first tags)) "a fresh, newer tag for every pass")
        (t/is (every? true? (map (fn [v] (every? (fn [[e xe]] (if (instance? Dual e) (= xe (ad/primal e)) (= xe e))) (map vector v x))) passes))
              "values unchanged, non-chunk inputs are plain numbers")
        (t/is (= [[1.0 0.0] [0.0 1.0]] (map partials (filter #(instance? Dual %) (first passes)))) "basis partials"))))
  (t/testing "documented: a branch at a tie of a comparison depends on the chunks"
    (let [min-branch (fn [[a b]] (if (ad/< a b) a b))]
      ;; one pass: a = 1 + e1, b = 1 + e2, a > b by the partials, b is taken
      (t/is (= [0.0 1.0] ((ad/gradient min-branch {:chunk-size 2}) [1.0 1.0])))
      ;; two passes: a > b (constant) in the first, a (constant) < b in the second, the constant is taken twice
      (t/is (= [0.0 0.0] ((ad/gradient min-branch {:chunk-size 1}) [1.0 1.0]))))))

(t/deftest gradient-test
  (t/testing "ForwardDiff reference, all chunk sizes give identical results"
    (doseq [[id f] scalar-functions
            {:keys [x value gradient]} (get-in @multivariate-reference [:scalar id])
            :let [n (count x)
                  default ((ad/gradient f) x)]]
      (t/is (close-to-reference? gradient default) (str id " n=" n))
      (t/is (m/delta-eq value (double (f x)) 1.0e-14 1.0e-14))
      (doseq [k [1 3 n (inc n)]]
        (t/is (= default ((ad/gradient f {:chunk-size k}) x)) (str id " n=" n " chunk " k)))))
  (t/testing "points"
    (let [f (fn [[x y]] (ad/* x (ad/pow y 2.0)))
          r [4.0 16.0]]
      (t/is (= r ((ad/gradient f) [4 2])) "longs")
      (t/is (= r ((ad/gradient f) (double-array [4 2]))) "double[]")
      (t/is (= r ((ad/gradient f) '(4.0 2.0))) "list")
      (t/is (= [1.0] ((ad/gradient (fn [[x]] x)) [3])) "one input")))
  (t/testing "empty point"
    (let [g (ad/gradient (fn [_] (throw (Exception. "f must not be called"))))]
      (t/is (= [] (g [])))
      (t/is (= [] (g nil)))
      (t/is (= [] (g (double-array 0))))))
  (t/testing "sizes around the chunk threshold: one pass (n=8) and two unequal passes (n=9) agree with one input per pass"
    (doseq [n [8 9]
            [id f] scalar-functions
            :let [x (vec (map (fn [i] (- (* 0.3 i) 0.8)) (range 1 (inc n))))]]
      (t/is (= ((ad/gradient f {:chunk-size 1}) x) ((ad/gradient f) x)) (str id " n=" n))
      (t/is (= ((ad/hessian f {:chunk-size 1}) x) ((ad/hessian f) x)) (str id " n=" n))))
  (t/testing "more than 32 inputs in one pass agree with one input per pass"
    (let [x (vec (map (fn [i] (- (* 0.05 i) 0.8)) (range 40)))]
      (doseq [[id f] scalar-functions]
        (t/is (= ((ad/gradient f {:chunk-size 1}) x) ((ad/gradient f {:chunk-size 40}) x)) (str id " gradient")))
      (t/is (= ((ad/jacobian (fn [v] [(reduce ad/* v) (reduce ad/+ v)]) {:chunk-size 1}) x)
               ((ad/jacobian (fn [v] [(reduce ad/* v) (reduce ad/+ v)]) {:chunk-size 40}) x))
            "jacobian")))
  (t/testing "constant function and unused inputs"
    (t/is (= [0.0 0.0] ((ad/gradient (constantly 5)) [1 2])))
    (t/is (= [1.0 0.0 0.0] ((ad/gradient (fn [v] (v 0)) {:chunk-size 1}) [1 2 3]))))
  (t/testing "errors"
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"neither a number nor a dual"
                            ((ad/gradient (fn [v] v)) [1 2])))
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"returning a function"
                            ((ad/gradient (fn [_] identity)) [1 2])))))

(t/deftest jacobian-test
  (t/testing "ForwardDiff reference, all chunk sizes give identical results"
    (doseq [[id f] vector-functions
            {:keys [x value jacobian]} (get-in @multivariate-reference [:vector id])
            :let [n (count x)
                  default ((ad/jacobian f) x)]]
      (t/is (close-to-reference? jacobian default) (str id))
      (t/is (close-to-reference? value (f x)))
      (t/is (= (count value) (count default)) "m rows")
      (t/is (every? #(= n (count %)) default) "n columns")
      (doseq [k [1 2 n]]
        (t/is (= default ((ad/jacobian f {:chunk-size k}) x)) (str id " chunk " k)))))
  (t/testing "numbers and constant outputs give zero rows"
    (t/is (= [[0.0 0.0] [1.0 0.0]] ((ad/jacobian (fn [[x _]] [1 x])) [5 6]))))
  (t/testing "empty point: f is called once"
    (let [calls (atom 0)
          j (ad/jacobian (fn [_] (swap! calls inc) [1.0 2.0 3.0]))]
      (t/is (= [[] [] []] (j [])))
      (t/is (= 1 @calls)))
    (t/is (= [] ((ad/jacobian (constantly [])) [1 2])) "no outputs"))
  (t/testing "errors"
    (doseq [[r msg] [[1.0 "number"] [(double-array [1.0]) "double[]"] [identity "function"]]]
      (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"sequential of numbers or duals"
                              ((ad/jacobian (constantly r)) [1 2])) msg))
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"Jacobian of a function returning a sequential with an element which is neither"
                            ((ad/jacobian (constantly [:a])) [1 2])) "an element is not a number")
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"Jacobian of a function returning a sequential with an element which is neither"
                            ((ad/jacobian (constantly [1.0 identity])) [1 2])) "an element is a function")
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"Jacobian of a function returning a sequential with an element which is neither"
                            ((ad/jacobian (constantly [nil])) [])) "an element is nil, empty point")
    (let [f (fn [v] (if (instance? Dual (v 0)) [(v 0)] [(v 0) (v 1)]))]
      (t/is (= {:expected 1 :actual 2} (thrown-data #((ad/jacobian f {:chunk-size 1}) [1 2])))
            "the number of values changes between the passes"))))

(t/deftest hessian-test
  (t/testing "ForwardDiff reference, all chunk sizes give identical results"
    (doseq [[id f] scalar-functions
            {:keys [x hessian]} (get-in @multivariate-reference [:scalar id])
            :let [n (count x)
                  h ((ad/hessian f) x)]]
      (t/is (close-to-reference? hessian h) (str id " n=" n))
      (t/is (every? true? (for [i (range n) j (range n)]
                            (m/<= (m/abs (m/- (double (get-in h [i j])) (double (get-in h [j i]))))
                                  (m/* 1.0e-14 (matrix-max-abs h)))))
            "symmetric up to the rounding")
      (doseq [k [1 3 n]]
        (t/is (= h ((ad/hessian f {:chunk-size k}) x)) (str id " n=" n " chunk " k)))))
  (t/testing "closed form"
    (t/is (= [[6.0 3.0] [3.0 0.0]] ((ad/hessian (fn [[x y]] (ad/* (ad/pow x 3.0) y))) [1 1]))))
  (t/testing "empty point"
    (t/is (= [] ((ad/hessian (fn [_] (throw (Exception. "f must not be called")))) []))))
  (t/testing "errors"
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"neither a number nor a dual"
                            ((ad/hessian (fn [v] v)) [1 2])))))

(t/deftest multivariate-nesting-test
  (t/testing "gradient inside derivative: d/dx of dg/dx (x, 2) for g = x^2 y^3 is 2 y^3"
    (t/is (= 16.0 ((D (fn [x] (nth ((ad/gradient (fn [[a b]] (ad/* (ad/pow a 2.0) (ad/pow b 3.0)))) [x 2.0]) 0))) 1.5))))
  (t/testing "the gradient inside a derivative returns duals"
    (let [seen (atom nil)]
      ((D (fn [x] (let [g ((ad/gradient (fn [[a b]] (ad/* a b))) [x 2.0])] (reset! seen g) (g 0)))) 1.0)
      (t/is (every? #(instance? Dual %) @seen))))
  (t/testing "derivative inside gradient"
    (t/is (= [3.0 2.0] ((ad/gradient (fn [[a b]] (ad/* a ((D (fn [t] (ad/* t b))) 1.0)))) [2.0 3.0]))))
  (t/testing "hessian inside derivative: d/dx H00 of x^3 y = 6y"
    (t/is (= 6.0 ((D (fn [x] (get-in ((ad/hessian (fn [[a b]] (ad/* (ad/pow a 3.0) b))) [x 1.0]) [0 0]))) 2.0))))
  (t/testing "a gradient function nested in itself gives the hessian rows"
    (let [f (:dense scalar-functions)
          x [0.1 -0.4 0.7]
          g (ad/gradient f)
          h ((ad/hessian f) x)]
      (doseq [i (range 3)]
        (t/is (= (h i) ((ad/gradient (fn [v] ((g v) i))) x))))))
  (t/testing "Siskind-Pearlmutter with gradients: d/dx (x * dg/dy (x + y)) at 1 is 1"
    (t/is (= 1.0 ((D (fn [x] (ad/* x (first ((ad/gradient (fn [[y]] (ad/+ x y))) [1.0]))))) 1.0)))))

(t/deftest multivariate-escape-test
  (t/testing "a dual of a finished gradient pass leaking into an enclosing derivative throws"
    (let [leak (atom nil)
          g (ad/gradient (fn [v] (reset! leak (v 0)) (v 0)))]
      (t/is (some? (thrown-data #((D (fn [x] (g [x]) (ad/+ x @leak))) 1.0))))))
  (t/testing "limitation: a dual of an earlier pass of the same call is older, it is treated as a constant"
    (let [first-seen (atom nil)
          f (fn [v] (when-not @first-seen (reset! first-seen (v 0))) (ad/+ (v 1) @first-seen))]
      ;; second pass: v1 + (stale dual of v0, a constant of the newer layer)
      (t/is (= [1.0 1.0] ((ad/gradient f {:chunk-size 1}) [1.0 2.0]))))))
