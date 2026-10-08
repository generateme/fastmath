(ns fastmath.optimization.problems-test
  (:require [fastmath.optimization.problems :as sut]
            [fastmath.optimization :as opt]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t]))

;; Reference values.
;;
;; Univariate problems: the collection of one dimensional global optimization problems by Gavana,
;; https://infinity77.net/global_optimization/test_functions_1d.html (published minima), where the values of
;; the collection were not at hand the minimum was found with a dense grid search (python, numpy, 4-5 million
;; points): problem10 (7.978665, -7.916727), problem21 (4.7954, -9.508350), problem22 (-1 at 5 pi / 2).
;; problem03 and problem08 are sums of k = 1..5, they used to have a sixth term and the published minima
;; -12.03125 and -14.508 were not reproduced.
;;
;; Multivariate problems: https://www.sfu.ca/~ssurjano/optimization.html (published minima and domains),
;; Himmelblau: Wikipedia.
;;
;; Derivatives and gradients are compared with central finite differences of the functions.

(def ^:private univariate
  ;; number -> [f* [x* ...]] documented in the docstrings
  {"02" [-1.899599 [5.145735]]
   "03" [-12.031249 [-6.774576 -0.491391 5.791794]]
   "04" [-3.850450 [2.868034]]
   "05" [-1.489073 [0.966090]]
   "06" [-0.824239 [0.679560]]
   "07" [-1.601308 [5.199780]]
   "08" [-14.508008 [-7.083506 -0.800321 5.482864]]
   "09" [-1.905961 [17.039199]]
   "10" [-7.916727 [7.978666]]
   "11" [-1.5 [(/ (* 2.0 Math/PI) 3.0) (/ (* 4.0 Math/PI) 3.0)]]
   "12" [-1.0 [Math/PI (/ (* 3.0 Math/PI) 2.0)]]
   "13" [-1.587401 [0.707107]]
   "14" [-0.788685 [0.224885]]
   "15" [-0.035534 [2.414214]]
   "18" [0.0 [2.0]]
   "20" [-0.063491 [1.195137]]
   "21" [-9.508350 [4.795400]]
   "22" [-1.0 [(/ (* 5.0 Math/PI) 2.0)]]})

(defn- var-of [prefix n] (deref (ns-resolve 'fastmath.optimization.problems (symbol (str prefix n)))))

(t/deftest public-api
  (let [publics (ns-publics 'fastmath.optimization.problems)]
    (t/is (seq publics))
    #_(t/is (empty? (remove (comp :doc meta val) publics)) "every public var is documented")
    (t/is (:doc (meta (the-ns 'fastmath.optimization.problems))))
    ;; names with the typo are gone, the new names are there
    (t/is (nil? (ns-resolve 'fastmath.optimization.problems '->auckley)))
    (t/is (nil? (ns-resolve 'fastmath.optimization.problems 'auckley-bounds)))
    (doseq [s '[->ackley ackley-bounds rosenbrock-gradient himmelblau-gradient beale-gradient sphere-gradient]]
      (t/is (some? (ns-resolve 'fastmath.optimization.problems s)) (str s)))))

;; univariate

(t/deftest univariate-bounds
  (doseq [n    (keys univariate)
          :let [bounds ((var-of "problem" (str n "-bounds")))]]
    (t/is (= 1 (count bounds)) n)
    (t/is (vector? bounds) n)
    (let [[lo hi] (first bounds)]
      (t/is (and (double? lo) (double? hi) (< lo hi)) n)
      ;; the published minimizers are inside of the domain
      (doseq [x (second (univariate n))]
        (t/is (<= lo x hi) (str n " " x))))))

(t/deftest univariate-minima
  (doseq [[n [fmin xs]] univariate
          :let [f (var-of "problem" n)
                [[lo hi]] ((var-of "problem" (str n "-bounds")))]]
    ;; the documented minimizers give the documented minimum ...
    (doseq [x xs]
      (t/is (m/delta-eq fmin (double (f x)) 1.0e-5) (str n " at " x)))
    ;; ... and no point of a dense grid is lower (a global minimum)
    (let [N 100000
          lowest (reduce (fn [^double acc ^long i] (m/min acc (double (f (m/+ lo (m/* (m// i (double N)) (m/- hi lo)))))))
                         Double/POSITIVE_INFINITY (range (inc N)))]
      (t/is (>= lowest (- fmin 1.0e-5)) (str n " grid minimum " lowest))
      (t/is (<= lowest (+ fmin 1.0e-4)) (str n " grid minimum " lowest)))
    ;; the function returns a double, also for integers
    (t/is (double? (f (long (Math/ceil lo)))) n)))

(t/deftest univariate-derivatives
  (doseq [n (keys univariate)
          :let [f (var-of "problem" n)
                df (var-of "dproblem" n)
                [[lo hi]] ((var-of "problem" (str n "-bounds")))
                ;; interior points, not the kink of problem18 at x = 3
                xs (remove #(and (= n "18") (< (Math/abs (- % 3.0)) 1.0e-3))
                           (map #(+ lo (* % (- hi lo))) [0.05 0.17 0.31 0.5 0.62 0.77 0.93]))]
          x xs
          :let [h 1.0e-6
                fd (/ (- (f (+ x h)) (f (- x h))) (* 2.0 h))
                d (df [x])]]
    (t/is (vector? d) n)
    (t/is (= 1 (count d)) n)
    (t/is (double? (first d)) n)
    ;; central differences: error about h^2 f''' + eps f / h
    (t/is (<= (Math/abs (- (first d) fd)) (* 1.0e-5 (+ 1.0 (Math/abs (first d))))) (str n " at " x ": " (first d) " vs " fd)))
  ;; the derivative vanishes at the minima inside of the domains (the documented digits limit the accuracy)
  (doseq [[n [_ xs]] univariate
          x xs
          :let [d (first ((var-of "dproblem" n) [x]))]]
    ;; 1e-2: the minimizers are given with 6 digits and the steepest function (problem05, f'' about 500) has
    ;; derivative 2e-3 at the rounded minimizer
    (t/is (m/delta-eq 0.0 d 1.0e-2) (str n " at " x ": " d)))
  ;; the derivative accepts any sequence with a number
  (t/is (= (sut/dproblem02 [3.0]) (sut/dproblem02 '(3.0)) (sut/dproblem02 (double-array [3.0])))))

(t/deftest problem18-both-branches
  ;; (x - 2)^2 up to x = 3, then 2 ln(x - 2) + 1, continuous with value 1 at x = 3
  (t/is (m/delta-eq 1.0 (sut/problem18 3.0) 1.0e-12))
  (t/is (m/delta-eq 1.0 (sut/problem18 3.000000001) 1.0e-8))
  (t/is (m/delta-eq 4.0 (sut/problem18 0.0)))
  (t/is (m/delta-eq (+ 1.0 (* 2.0 (Math/log 2.0))) (sut/problem18 4.0)))
  (t/is (= [2.0] (sut/dproblem18 [3.0])))
  (t/is (= [1.0] (sut/dproblem18 [4.0])) "2 / (x - 2)"))

;; multivariate

(t/deftest multivariate-bounds
  (doseq [[nm bounds expected-count] [["bukin-no-6" (sut/bukin-no-6-bounds) 2] ["cross-in-tray" (sut/cross-in-tray-bounds) 2]
                                      ["drop-wave" (sut/drop-wave-bounds) 2] ["egg-holder" (sut/egg-holder-bounds) 2]
                                      ["himmelblau" (sut/himmelblau-bounds) 2] ["beale" (sut/beale-bounds) 2]
                                      ["rosenbrock" (sut/rosenbrock-bounds 4) 4] ["sphere" (sut/sphere-bounds 3) 3]
                                      ["ackley" (sut/ackley-bounds 5) 5]]]
    (t/is (= expected-count (count bounds)) nm)
    (t/is (every? (fn [[lo hi]] (and (double? lo) (double? hi) (< lo hi))) bounds) nm))
  ;; published domains
  (t/is (= [[-4.5 4.5] [-4.5 4.5]] (sut/beale-bounds)) "the second range used to be [4.5 4.5]")
  (t/is (= [[-32.768 32.768] [-32.768 32.768]] (sut/ackley-bounds 2)) "usually [-32.768, 32.768]")
  (t/is (= [[-5.12 5.12] [-5.12 5.12] [-5.12 5.12]] (sut/sphere-bounds 3)))
  (t/is (= [[-5.0 10.0] [-5.0 10.0]] (sut/rosenbrock-bounds 2)))
  (t/is (= [[-15.0 -5.0] [-3.0 3.0]] (sut/bukin-no-6-bounds)))
  (t/is (= [[-10.0 10.0] [-10.0 10.0]] (sut/cross-in-tray-bounds)))
  (t/is (= [[-5.12 5.12] [-5.12 5.12]] (sut/drop-wave-bounds)))
  (t/is (= [[-512.0 512.0] [-512.0 512.0]] (sut/egg-holder-bounds)))
  (t/is (= [[-5.0 5.0] [-5.0 5.0]] (sut/himmelblau-bounds)))
  ;; any number of dimensions, one is the boundary
  (doseq [f [sut/ackley-bounds sut/rosenbrock-bounds sut/sphere-bounds]]
    (t/is (= 1 (count (f 1))))
    (t/is (empty? (f 0)))
    (t/is (= 50 (count (f 50))))))

(t/deftest multivariate-minima
  ;; the minimizers are inside of the domains
  (let [inside? (fn [bounds pt] (every? true? (map (fn [[lo hi] x] (<= lo x hi)) bounds pt)))]
    (t/is (inside? (sut/himmelblau-bounds) [3.0 2.0]))
    (t/is (inside? (sut/beale-bounds) [3.0 0.5]))
    (t/is (inside? (sut/bukin-no-6-bounds) [-10.0 1.0]))
    (t/is (inside? (sut/egg-holder-bounds) [512.0 404.2319]))
    (t/is (inside? (sut/rosenbrock-bounds 5) [1.0 1.0 1.0 1.0 1.0])))
  ;; himmelblau: four minima
  (doseq [pt [[3.0 2.0] [-2.805118 3.131312] [-3.779310 -3.283186] [3.584428 -1.848126]]]
    (t/is (m/delta-eq 0.0 (sut/himmelblau pt) 1.0e-5) (pr-str pt)))
  (t/is (m/delta-eq 0.0 (sut/himmelblau [3 2])) "integers")
  ;; beale
  (t/is (= 0.0 (sut/beale [3.0 0.5])))
  (t/is (> (sut/beale [4.5 4.5]) 1.0e4) "peaks at the corners")
  ;; bukin 6
  (t/is (m/delta-eq 0.0 (sut/bukin-no-6 [-10.0 1.0])))
  (t/is (> (sut/bukin-no-6 [-5.0 3.0]) 1.0))
  ;; cross-in-tray: four minima
  (doseq [pt [[1.34941 1.34941] [-1.34941 1.34941] [1.34941 -1.34941] [-1.34941 -1.34941]]]
    (t/is (m/delta-eq -2.06261 (sut/cross-in-tray pt) 1.0e-5) (pr-str pt)))
  ;; drop-wave
  (t/is (m/delta-eq -1.0 (sut/drop-wave [0.0 0.0])))
  (t/is (> (sut/drop-wave [3.0 3.0]) -1.0))
  (t/is (m/delta-eq -1.0 (sut/drop-wave (double-array [0.0 0.0]))))
  ;; egg-holder
  (t/is (m/delta-eq -959.6407 (sut/egg-holder [512.0 404.2319]) 1.0e-4))
  ;; rosenbrock: minimum in the valley for any number of dimensions
  (doseq [n [2 3 10]]
    (t/is (= 0.0 (sut/rosenbrock (repeat n 1.0))) (str n)))
  ;; sphere
  (t/is (= 0.0 (sut/sphere [0.0 0.0 0.0])))
  (t/is (= 14.0 (sut/sphere [1.0 2.0 3.0])))
  (t/is (= 25.0 (sut/sphere [3 4])))
  (t/is (= 0.0 (sut/sphere [])) "no coordinates"))

(t/deftest rosenbrock-short-points
  (t/is (= 0.0 (sut/rosenbrock [])))
  (t/is (= 0.0 (sut/rosenbrock [5.0])) "no consecutive pair")
  (t/is (= 101.0 (sut/rosenbrock [0.0 1.0])))
  (t/is (= 4.0 (sut/rosenbrock [-1.0 1.0])) "(x1 - 1)^2 = 4, x2 = x1^2")
  (t/is (m/delta-eq 24.2 (sut/rosenbrock [-1.2 1.0]) 1.0e-9) "the classical start point")
  ;; the sum goes over consecutive pairs
  (t/is (= (+ (sut/rosenbrock [2.0 3.0]) (sut/rosenbrock [3.0 4.0])) (sut/rosenbrock [2.0 3.0 4.0]))))

(t/deftest ackley
  (let [f (sut/->ackley)]
    (doseq [n [1 2 5 20]]
      (t/is (m/delta-eq 0.0 (f (repeat n 0.0)) 1.0e-12) (str n " dimensions")))
    (t/is (double? (f [0.0])))
    ;; a point away from the minimum has a higher value, the function is symmetric
    (t/is (> (f [1.0 1.0]) 3.0))
    (t/is (m/delta-eq (f [1.0 -2.0 0.5]) (f [-1.0 2.0 -0.5])))
    ;; reference: the formula of the library of simulation experiments evaluated with numpy
    ;; -a exp(-b sqrt(mean(x^2))) - exp(mean(cos(c x))) + a + e; the function used exp(-b sum(x^2)) before:
    ;; wrong for more than one non-zero coordinate or more than one dimension
    (doseq [[pt expected] [[[1 1] 3.6253849384403627] [[2 0 0 0] 3.6253849384403627] [[0.5 -1.5 2.5] 8.13725728226161]
                           [[0.3 -0.7] 4.0262342249673075] [[10.0] 17.293294335267746]
                           [[-32.768 32.768 1.0] 21.118865107019342] [[1 2 3 4 5] 9.697286414061548]]]
      (t/is (m/delta-eq expected (f pt) 1.0e-9) (pr-str pt)))
    ;; explicit defaults and other parameters
    (t/is (= (f [0.3 -0.7]) ((sut/->ackley {:a 20.0 :b 0.2 :c m/TWO_PI}) [0.3 -0.7])))
    (t/is (m/delta-eq 3.1499641909857066 ((sut/->ackley {:a 10.0 :b 0.3 :c 3.0}) [0.3 -0.7]) 1.0e-9))
    (t/is (m/delta-eq 0.0 ((sut/->ackley {:a 10.0 :b 0.3 :c 3.0}) [0.0 0.0]) 1.0e-12) "the minimum does not depend on the parameters")))

;; gradients

(defn- numeric-gradient
  [f pt]
  (let [h 1.0e-6]
    (mapv (fn [i]
            (let [p+ (assoc (vec pt) i (+ (nth pt i) h))
                  p- (assoc (vec pt) i (- (nth pt i) h))]
              (/ (- (f p+) (f p-)) (* 2.0 h))))
          (range (count pt)))))

(defn- gradient-close? [expected actual]
  (and (= (count expected) (count actual))
       (every? true? (map (fn [e a] (<= (Math/abs (- e a)) (* 1.0e-5 (+ 1.0 (Math/abs e))))) expected actual))))

(t/deftest gradients-vs-finite-differences
  (doseq [[nm f g pts]
          [["himmelblau" sut/himmelblau sut/himmelblau-gradient
            [[0.0 0.0] [1.0 1.0] [3.0 2.0] [-2.0 3.5] [4.9 -4.9] [0.1 -0.3] [-3.7 -3.2]]]
           ["beale" sut/beale sut/beale-gradient
            [[0.0 0.0] [1.0 1.0] [3.0 0.5] [-2.0 3.5] [4.5 4.5] [-4.5 -4.5] [0.3 -0.7] [2.0 -1.0]]]
           ["rosenbrock 2" sut/rosenbrock sut/rosenbrock-gradient
            [[0.0 0.0] [1.0 1.0] [-1.2 1.0] [3.0 -2.0] [10.0 -5.0] [0.5 0.25]]]
           ["rosenbrock 3" sut/rosenbrock sut/rosenbrock-gradient
            [[0.0 0.0 0.0] [1.0 1.0 1.0] [-1.2 1.0 0.5] [2.0 -3.0 4.0] [0.1 0.2 0.3]]]
           ["rosenbrock 7" sut/rosenbrock sut/rosenbrock-gradient
            [[1.0 2.0 3.0 4.0 5.0 6.0 7.0] [-1.0 0.5 -0.5 2.0 -2.0 0.1 0.0] [0.9 0.9 0.9 0.9 0.9 0.9 0.9]]]
           ["sphere" sut/sphere sut/sphere-gradient
            [[0.0 0.0] [1.0 2.0 3.0] [-5.0 5.0 0.1 0.0]]]]
          pt pts]
    (let [gv (g pt)]
      (t/is (every? number? gv) (str nm (pr-str pt)))
      (t/is (gradient-close? (numeric-gradient f pt) gv) (str nm (pr-str pt) " " gv " vs " (numeric-gradient f pt))))))

(t/deftest gradient-forms
  ;; every gradient is a vector of doubles with a value for every coordinate
  (doseq [[g pt] [[sut/himmelblau-gradient [1.0 2.0]] [sut/beale-gradient [1.0 2.0]] [sut/rosenbrock-gradient [1.0 2.0 3.0]]]]
    (let [r (g pt)]
      (t/is (vector? r))
      (t/is (= (count pt) (count r)))
      (t/is (every? double? r))))
  ;; integers, lazy sequences and arrays are accepted
  (t/is (= (sut/himmelblau-gradient [1.0 2.0]) (sut/himmelblau-gradient [1 2]) (sut/himmelblau-gradient (map double [1 2])) (sut/himmelblau-gradient (double-array [1 2]))))
  (t/is (= (sut/beale-gradient [1.0 2.0]) (sut/beale-gradient [1 2]) (sut/beale-gradient (double-array [1 2]))))
  (t/is (= (sut/rosenbrock-gradient [1.0 2.0 3.0]) (sut/rosenbrock-gradient [1 2 3]) (sut/rosenbrock-gradient (map double [1 2 3])) (sut/rosenbrock-gradient (double-array [1 2 3]))))
  ;; sphere: the form of the argument
  (t/is (= [2.0 4.0] (sut/sphere-gradient [1.0 2.0])))
  (t/is (= [2.0 4.0] (vec (sut/sphere-gradient (double-array [1.0 2.0])))))
  ;; short points: nothing to differentiate
  (t/is (= [] (sut/rosenbrock-gradient [])))
  (t/is (= [0.0] (sut/rosenbrock-gradient [5.0])))
  ;; closed form of the rosenbrock gradient at (0, 0): (-2, 0), and at (1, 2): (-400 (2 - 1) + 0, 200)
  (t/is (= [-2.0 0.0] (sut/rosenbrock-gradient [0.0 0.0])))
  (t/is (= [-400.0 200.0] (sut/rosenbrock-gradient [1.0 2.0]))))

(t/deftest gradients-vanish-at-minima
  (doseq [pt [[3.0 2.0] [-2.805118 3.131312] [-3.779310 -3.283186] [3.584428 -1.848126]]]
    (t/is (every? #(m/delta-eq 0.0 % 1.0e-3) (sut/himmelblau-gradient pt)) (pr-str pt)))
  (t/is (= [0.0 0.0] (sut/beale-gradient [3.0 0.5])))
  (t/is (= [0.0 0.0] (sut/sphere-gradient [0.0 0.0])))
  (doseq [n [2 3 10]]
    (t/is (= (repeat n 0.0) (sut/rosenbrock-gradient (repeat n 1.0))) (str n))))

;; the problems with the optimizers

(t/deftest problems-with-optimizers
  (let [run (fn [method f bounds initial gradient]
              (opt/minimize method f (cond-> {:bounds bounds :initial initial}
                                       gradient (assoc :gradient gradient))))]
    ;; the gradients are accepted by the optimizers which use them
    (doseq [method [:lbfgsb :conjugate-gradient]]
      (let [[pt val] (run method sut/himmelblau (sut/himmelblau-bounds) [1.0 1.0] sut/himmelblau-gradient)]
        (t/is (m/delta-eq 0.0 val 1.0e-3) (str method))
        (t/is (some #(v/delta-eq (vec pt) % 1.0e-3) [[3.0 2.0] [-2.805118 3.131312] [-3.779310 -3.283186] [3.584428 -1.848126]]) (str method)))
      (let [[pt val] (run method sut/rosenbrock (sut/rosenbrock-bounds 5) [-1.2 1.0 -1.2 1.0 -1.2] sut/rosenbrock-gradient)]
        (t/is (m/delta-eq 0.0 val 1.0e-3) (str method))
        (t/is (v/delta-eq (vec pt) [1.0 1.0 1.0 1.0 1.0] 1.0e-2) (str method))))
    (let [[pt val] (run :lbfgsb sut/beale (sut/beale-bounds) [2.0 0.0] sut/beale-gradient)]
      (t/is (m/delta-eq 0.0 val 1.0e-4))
      (t/is (v/delta-eq (vec pt) [3.0 0.5] 1.0e-3)))
    (let [[pt val] (run :lbfgsb sut/sphere (sut/sphere-bounds 4) [1.0 -2.0 3.0 -4.0] sut/sphere-gradient)]
      (t/is (m/delta-eq 0.0 val 1.0e-9))
      (t/is (v/delta-eq (vec pt) [0.0 0.0 0.0 0.0] 1.0e-4)))
    ;; with and without the gradient the optimizer ends in the same place
    (t/is (v/delta-eq (first (run :lbfgsb sut/himmelblau (sut/himmelblau-bounds) [1.0 1.0] nil))
                      (first (run :lbfgsb sut/himmelblau (sut/himmelblau-bounds) [1.0 1.0] sut/himmelblau-gradient)) 1.0e-5))
    ;; global problems on their domains, scanned from many points
    (let [[pt val] (opt/scan-and-minimize :lbfgsb (sut/->ackley) {:bounds (sut/ackley-bounds 2) :N 2000 :n 20 :gradient-h 1.0e-6})]
      (t/is (some? pt))
      (t/is (< val 1.0e-3)))
    (let [[_ val] (opt/scan-and-minimize :lbfgsb sut/egg-holder {:bounds (sut/egg-holder-bounds) :N 2000 :n 50})]
      (t/is (< val -950.0)))
    (let [[_ val] (opt/scan-and-minimize :lbfgsb sut/cross-in-tray {:bounds (sut/cross-in-tray-bounds) :N 500 :n 20})]
      (t/is (m/delta-eq -2.06261 val 1.0e-3)))
    ;; univariate problems with the optimizers: the global minimum of problem03 (the published one, not the lower value of the sum with six terms)
    (let [[pt val] (opt/scan-and-minimize :brent sut/problem03 {:bounds (sut/problem03-bounds) :N 400 :n 40})]
      (t/is (m/delta-eq -12.031249 val 1.0e-5))
      (t/is (some #(m/delta-eq % pt 1.0e-3) [-6.774576 -0.491391 5.791794])))
    (let [[_ val] (opt/scan-and-minimize :brent sut/problem08 {:bounds (sut/problem08-bounds) :N 400 :n 40})]
      (t/is (m/delta-eq -14.508008 val 1.0e-5)))
    (let [[pt val] (opt/minimize :lbfgsb (fn [[x]] (sut/problem02 x)) {:bounds (sut/problem02-bounds) :gradient sut/dproblem02})]
      (t/is (m/delta-eq 5.145735 (first pt) 1.0e-4))
      (t/is (m/delta-eq -1.899599 val 1.0e-5)))
    (let [[pt val] (opt/minimize :lbfgsb sut/problem02 {:bounds (sut/problem02-bounds) :gradient sut/dproblem02 :vector-arg? false})]
      (t/is (m/delta-eq 5.145735 (first pt) 1.0e-4))
      (t/is (m/delta-eq -1.899599 val 1.0e-5)))))
