(ns fastmath.calculus.autodiff-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.calculus.autodiff :as ad])
  (:import [fastmath.calculus.autodiff Dual]))

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
(defn- partials [^Dual d] (.partials d))

(defn- same?
  "Structural equality of numbers and duals (tags, values and partials)."
  [a b]
  (if (and (instance? Dual a) (instance? Dual b))
    (and (= (tag a) (tag b))
         (same? (ad/primal a) (ad/primal b))
         (= (count (partials a)) (count (partials b)))
         (every? true? (map same? (partials a) (partials b))))
    (= a b)))

(defn- dual-is?
  "Is `d` a dual with the given tag, value and partials?"
  [d t v p]
  (and (instance? Dual d)
       (same? d (ad/dual t v p))))

(defn- nan? [x] (and (number? x) (m/nan? (double x))))

(def ^:private D ad/derivative)

;; dual numbers

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
      (t/is (not (= d (ad/dual 1 2 [1]))) "clojure.core/= is identity")
      (t/is (= d d))))
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
    (let [fresh-tag @#'ad/fresh-tag
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
      (t/is (= [-1.0] (partials d)) "the derivative 1/x is defined, the value is not"))))

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
    (t/is (ad/= 2 2.0))
    (t/is (ad/= 0.0 -0.0))
    (t/is (ad/< 1) "one argument")
    (t/is (ad/= ##NaN) "one argument, as clojure.core")
    (t/is (ad/< 1 2 3 4))
    (t/is (not (ad/< 1 2 2 4)))
    (t/is (ad/<= 1 2 2 4))
    (t/is (ad/= 1 1 1))
    (t/is (not (ad/= 1 1 2)))
    (t/is (ad/> 4 3 1/2)))
  (t/testing "NaN is unordered"
    (t/is (not-any? (fn [f] (f ##NaN 1.0)) [ad/< ad/<= ad/> ad/>= ad/=]))
    (t/is (not-any? (fn [f] (f (ad/dual 1 ##NaN [1]) 1.0)) [ad/< ad/<= ad/> ad/>= ad/=])))
  (t/testing "lexicographic on value and partials"
    (let [x (ad/dual 1 1 [1])]
      (t/is (not (ad/= x 1.0)) "a dual with nonzero partials is not equal to its value")
      (t/is (ad/= (ad/dual 1 1 [0]) 1.0) "zero partials")
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
    (let [sinc (fn [x] (if (ad/= x 0) 1.0 (ad// (ad/sin x) x)))
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
