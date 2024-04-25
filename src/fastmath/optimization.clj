(ns fastmath.optimization
  "Optimization.

  Namespace provides various optimization methods.

  * Brent (1d functions)
  * Bobyqa (2d+ functions)
  * Powell
  * Nelder-Mead
  * Multidirectional simplex
  * CMAES
  * Gradient
  * L-BFGS-B
  * Bayesian Optimization (see below)
  * Linear optimization

  All optimizers require bounds.

  ## Optimizers

  To optimize functions call one of the following functions:

  * [[minimize]] or [[maximize]] - to perform actual optimization
  * [[scan-and-minimize]] or [[scan-and-maximize]] - functions find initial point using brute force and then perform optimization paralelly for best initialization points. Brute force scan is done using jitter low discrepancy sequence generator.

  You can also create optimizer (function which performs optimization) by calling [[minimizer]] or [[maximizer]]. Optimizer accepts initial point.

  All above accept:

  * one of the optimization method, ie: `:brent`, `:bobyqa`, `:nelder-mead`, `:multidirectional-simplex`, `:cmaes`, `:gradient`, `:bfgs` and `:lbfgsb`
  * function to optimize
  * parameters as a map

  For parameters meaning refer [Optim package](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html?org/apache/commons/math3/optim/package-summary.html)
  
  ### Common parameters

  * `:bounds` (obligatory) - search ranges for each dimensions as a seqence of [low high] pairs
  * `:initial` - initial point other then mid of the bounds as vector
  * `:max-evals` - maximum number of function evaluations
  * `:max-iters` - maximum number of algorithm interations
  * `:bounded?` - should optimizer force to keep search within bounds (some algorithms go outside desired ranges)
  * `:stats?` - return number of iterations and evaluations along with result
  * `:rel` and `:abs` - relative and absolute accepted errors

  For `scan-and-...` functions additionally you can provide:

  * `:N` - number of brute force iterations
  * `:n` - fraction of N which are used as initial points to parallel optimization
  * `:jitter` - jitter factor for sequence generator (for scanning domain)
  
  ### Specific parameters

  * BOBYQA - `:number-of-points`, `:initial-radius`, `:stopping-radius`
  * Nelder-Mead - `:rho`, `:khi`, `:gamma`, `:sigma`, `:side-length`
  * Multidirectional simples - `:khi`, `:gamma`, `:side-length`
  * CMAES - `:check-feasable-count`, `:diagonal-only`, `:stop-fitness`, `:active-cma?`, `:population-size`
  * Gradient - `:bracketing-range`, `:formula` (`:polak-ribiere` or `:fletcher-reeves`), `:gradient-h` (finite differentiation step, default: `0.01`) 

  ## Bayesian Optimization

  Bayesian optimizer can be used for optimizing expensive to evaluate black box functions. Refer this [article](http://krasserm.github.io/2018/03/21/bayesian-optimization/) or this [article](https://nextjournal.com/a/LKqpdDdxiggRyHhqDG5FH?token=Ss1Qq3MzHWN8ZyEt9UC1ZZ)

  ## Linear optimization "
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.kernel :as k]
            [fastmath.interpolation.gp :as gp]
            [fastmath.optimization.lbfgsb :as lbfgsb])
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunction ObjectiveFunctionGradient]
           [org.apache.commons.math3.optim.univariate SearchInterval BrentOptimizer UnivariateObjectiveFunction UnivariatePointValuePair]
           [org.apache.commons.math3.optim BaseOptimizer OptimizationData MaxEval MaxIter SimpleBounds SimpleValueChecker InitialGuess PointValuePair]
           [org.apache.commons.math3.optim.linear LinearObjectiveFunction LinearConstraint
            Relationship LinearConstraintSet SimplexSolver NonNegativeConstraint PivotSelectionRule]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction MultivariateVectorFunction]
           [org.apache.commons.math3.optim.nonlinear.scalar MultivariateFunctionMappingAdapter]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer PowellOptimizer NelderMeadSimplex SimplexOptimizer MultiDirectionalSimplex CMAESOptimizer CMAESOptimizer$PopulationSize CMAESOptimizer$Sigma]
           [org.apache.commons.math3.optim.nonlinear.scalar.gradient NonLinearConjugateGradientOptimizer NonLinearConjugateGradientOptimizer$Formula]
           [org.generateme.lbfgsb LBFGSB]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def ^:private univariate-set #{:brent})
(def ^:private multivariate-set #{:bobyqa :powell :nelder-mead :multidirectional-simplex :cmaes :gradient :bfgs :lbfgsb})
(def ^:private unbounded-set #{:powell :nelder-mead :multidirectional-simplex :gradient})

(defn- brent
  [{:keys [^double rel ^double abs]
    :or {rel 1.0e-6  abs 1.0e-10}}]
  (BrentOptimizer. rel abs))

(defn- bobyqa
  [{:keys [number-of-points ^int dim initial-radius stopping-radius]
    :or {initial-radius BOBYQAOptimizer/DEFAULT_INITIAL_RADIUS
         stopping-radius BOBYQAOptimizer/DEFAULT_STOPPING_RADIUS}}]
  (let [number-of-points (or number-of-points (/ (+ 3 (* 3 dim)) 2))]
    (BOBYQAOptimizer. number-of-points initial-radius stopping-radius)))

(defn- powell
  [{:keys [^double rel ^double abs]
    :or {rel 1.0e-6  abs 1.0e-10}}]
  (PowellOptimizer. rel abs))

(defn- nelder-mead
  [{:keys [^int dim ^double rho ^double khi ^double gamma ^double sigma ^double side-length]
    :or {rho 1.0 khi 2.0 gamma 0.5 sigma 0.5 side-length 1.0}}]
  (NelderMeadSimplex. dim side-length rho khi gamma sigma))

(defn- multidirectional-simplex
  [{:keys [^int dim ^double khi ^double gamma ^double side-length]
    :or {khi 2.0 gamma 0.5 side-length 1.0}}]
  (MultiDirectionalSimplex. dim side-length khi gamma))

(defn- cmaes
  [{:keys [^double rel ^double abs active-cma?
           ^int max-iters ^int check-feasable-count ^int diagonal-only
           ^double stop-fitness]
    :or {rel 1.0e-6
         abs 1.0e-10
         active-cma? true
         max-iters Integer/MAX_VALUE
         check-feasable-count 0
         diagonal-only 0
         stop-fitness 1.0e-6}}]
  (let [checker (SimpleValueChecker. rel abs)]
    (CMAESOptimizer. max-iters stop-fitness (boolean active-cma?) diagonal-only
                     check-feasable-count r/default-rng false checker)))

(defn- simplex
  [{:keys [^double rel ^double abs]
    :or {rel 1.0e-6  abs 1.0e-10}}]
  (SimplexOptimizer. rel abs))

(defn- non-linear-gradient
  [{:keys [^double rel ^double abs ^double bracketing-range formula]
    :or {rel 1.0e-8 abs 1.0e-8 bracketing-range 1.0e-8 formula :polak-ribiere}}]
  (let [checker (SimpleValueChecker. rel abs)]
    (NonLinearConjugateGradientOptimizer. (if (= formula :polak-ribiere)
                                            NonLinearConjugateGradientOptimizer$Formula/POLAK_RIBIERE
                                            NonLinearConjugateGradientOptimizer$Formula/FLETCHER_REEVES)
                                          checker, rel abs bracketing-range)))

(def ^:private optimizers
  {:brent brent
   :bobyqa bobyqa
   :powell powell
   :nelder-mead simplex
   :multidirectional-simplex simplex
   :cmaes cmaes
   :gradient non-linear-gradient})

(defn- wrap-univariate-function ^UnivariateFunction [f] (reify UnivariateFunction (value [_ x] (f x))))
(defn- wrap-univariate-objective-function [f] (UnivariateObjectiveFunction. (wrap-univariate-function f)))

(defn- multivariate-function [f]
  (reify
    MultivariateFunction
    (value [_ xs] (apply f xs))))

(defn- negative-multivariate-function [f]
  (reify
    MultivariateFunction
    (value [_ xs] (- ^double (apply f xs)))))


(defn- wrap-multivariate-function [f] (ObjectiveFunction. (multivariate-function f)))

(defn- wrap-objective-function [f] (ObjectiveFunction. f))

(defn- finite-differences
  [^MultivariateFunction f ^doubles xs ^double step]
  (let [step2 (* 2.0 step)]
    (double-array (map-indexed (fn [^long id ^double x]
                                 (let [a (aclone ^doubles xs)
                                       v1 (do (aset a id (+ x step))
                                              (.value f a))
                                       v2 (do (aset a id (- x step))
                                              (.value f a))]
                                   (/ (- v1 v2) step2))) xs))))

(defn- multivariate-gradient
  "Calculate gradient numerically."
  [^MultivariateFunction f ^double step] 
  (reify MultivariateVectorFunction
    (value [_ xs]
      (finite-differences f xs step))))

(defn- wrap-objective-function-gradient [f ^double step]
  (ObjectiveFunctionGradient. (multivariate-gradient f step)))

(defn- wrap-evals [evals] (if evals (MaxEval. evals) (MaxEval/unlimited)))
(defn- wrap-iters [iters] (if iters (MaxIter. iters) (MaxIter/unlimited)))

(defn- infer-lo-high
  [lo high]
  [(or lo m/EPSILON)
   (or high (- 1.0 m/EPSILON))])

(defn- search-interval
  [[lo high] init]
  (let [[lo high] (infer-lo-high lo high)]
    (if init
      (SearchInterval. lo high (if (sequential? init) (first init) init))
      (SearchInterval. lo high))))

(defn- multi-bounds
  [bounds]
  (SimpleBounds. (double-array (map first bounds))
                 (double-array (map second bounds))))

(defn- wrap-bounds
  [method bounds]
  (when (and bounds (multivariate-set method))
    (multi-bounds bounds)))

(defn- mid-point
  [bounds]
  (map (fn [[^double l ^double h]] (* 0.5 (+ l h))) bounds))

(defn- initial-guess
  [bounds initial ^MultivariateFunctionMappingAdapter mfma]
  (let [guess (double-array (if initial
                              (if (sequential? initial) initial [initial])
                              (mid-point bounds)))]
    (InitialGuess. (if mfma (.boundedToUnbounded mfma guess) guess))))

(defn- wrap-initial
  [method bounds initial mfma]
  (if (univariate-set method)
    (search-interval bounds initial)
    (initial-guess bounds initial mfma)))

(defn- wrap-function
  [method f]
  (if (univariate-set method)
    (wrap-univariate-objective-function f)
    (wrap-multivariate-function f)))

(defn- parse-result
  ([res] (parse-result nil res))
  ([^MultivariateFunctionMappingAdapter mfma res]
   (condp instance? res
     UnivariatePointValuePair (let [^UnivariatePointValuePair res res]
                                [(list (.getPoint res)) (.getValue res)])
     PointValuePair (let [^PointValuePair res res]
                      [(seq (if mfma
                              (.unboundedToBounded mfma (.getPointRef res))
                              (.getPointRef res))) (.getValue res)])
     res)))

(defn- find-dimensions
  ^long [bounds]
  (if (sequential? (first bounds)) (count bounds) 1))

(defn- fix-brent-bounds
  [method bounds]
  (if (and (= method :brent)
           (sequential? (first bounds)))
    (first bounds) bounds))

(defn- maybe-stats?
  [stats? ^BaseOptimizer optimizer res]
  (if-not stats?
    res
    {:result res
     :evaluations (.getEvaluations optimizer)
     :iterations (.getIterations optimizer)}))

(defn- optimizer
  [method f {:keys [max-evals max-iters goal bounds stats? population-size bounded? gradient-h]
             :or {gradient-h 0.0001}
             :as config}]
  (if (= method :lbfgsb)
    (lbfgsb/lbfgsb-fn f (assoc config :bounded? true))
    (do
      (assert (not (nil? bounds)) "Provide bounds")
      (let [bounds (fix-brent-bounds method bounds)
            dim (find-dimensions bounds)
            config (assoc config :dim dim :bounds bounds)
            
            ^SimpleBounds b (wrap-bounds method bounds)

            bounded? (and bounded? (unbounded-set method))

            mfma (when bounded?
                   (MultivariateFunctionMappingAdapter. (multivariate-function f) (.getLower b) (.getUpper b)))
            
            ;; create initial optimization data
            base-opt-data [(wrap-evals max-evals)
                           (wrap-iters max-iters)
                           (if (= goal :maximize) GoalType/MAXIMIZE GoalType/MINIMIZE)
                           (if bounded?
                             (wrap-objective-function mfma)
                             (wrap-function method f))]
            
            ;; powell and simplex methods do not accept bounds
            base-opt-data (if (and b (not (unbounded-set method))) (conj base-opt-data b) base-opt-data)

            ;; simplex methods should have also specific siumplex algorithms, also for cmaes we add additional stuff
            base-opt-data (case method
                            :nelder-mead (conj base-opt-data (nelder-mead config))
                            :multidirectional-simplex (conj base-opt-data (multidirectional-simplex config))
                            ;; when function is wrapped to bounding adapter, we need to use it to calculate gradient
                            :gradient (conj base-opt-data (wrap-objective-function-gradient (or mfma (multivariate-function f))
                                                                                            gradient-h))
                            :cmaes (conj base-opt-data
                                         (CMAESOptimizer$PopulationSize. (or population-size (int (+ 4.5 (* 3.0 (m/ln dim))))))
                                         (CMAESOptimizer$Sigma. (double-array (map #(* 0.75 (- ^double %2 ^double %1)) (.getLower b) (.getUpper b)))))
                            base-opt-data)
            
            builder (optimizers method)        
            ^BaseOptimizer optimizer (builder config)]
        
        (fn local-optimizer
          ([] (local-optimizer nil))
          ([init] (->> (conj base-opt-data (wrap-initial method bounds init mfma))
                       (into-array OptimizationData)
                       (.optimize optimizer)
                       (parse-result mfma)
                       (maybe-stats? stats? optimizer))))))))

(defn minimizer
  "Create optimizer which minimizes function.

  Returns function which performs optimization for optionally given initial point."
  [method f config] (optimizer method f (assoc config :goal :minimize)))

(defn maximizer
  "Create optimizer which maximizes function.

  Returns function which performs optimization for optionally given initial point."
  [method f config] (optimizer method f (assoc config :goal :maximize)))

(defmacro ^:private with-optimizer [opt m f c] `((~opt ~m ~f ~c) (:initial ~c)))

#_(defn optimize [method f config] (with-optimizer optimizer method f config))
(defn minimize
  "Minimize given function.

  Parameters: optimization method, function and configuration."
  [method f config] (with-optimizer minimizer method f config))

(defn maximize
  "Maximize given function.

  Parameters: optimization method, function and configuration."
  [method f config] (with-optimizer maximizer method f config))

(defn- goal-comparator
  [goal]
  (if (= goal :minimize)
    (fn [^double a ^double b] (< a b))
    (fn [^double a ^double b] (> a b))))

(defn- generate-points
  [method f bounds goal N n jitter]
  (let [bounds (fix-brent-bounds method bounds)
        dim (find-dimensions bounds)
        [lo high inter genf] (if (= method :brent)
                               [(first bounds) (second bounds) m/lerp f]
                               [(mapv first bounds) (mapv second bounds) v/einterpolate (partial apply f)])
        N (int N)
        N (max (m/fpow 3 dim) N)
        n (double n)
        nbest (max 10 (long (if (> n 1.0) n (m/floor (* n N)))))
        gen (r/jittered-sequence-generator (if (< dim 15) :r2 :sobol) dim jitter)]
    (->> (take N (if (and (not= method :brent)
                          (m/one? dim)) (map vector gen) gen))
         (map #(let [p (inter lo high %)]
                 [(genf p) p]))
         (sort-by first (goal-comparator goal))
         (take nbest)
         (map second))))

(defn- scan-and-
  "For cheap functions, scan domain and bruteforcely search for set of minimal values, then use part of this values as initial points for parallel optimization.

  Additional parameters in config:

  * N - number of total grid points
  * n - fraction of total points N) used for optimization (default: 0.05, minimum 10)"
  [optimizer-fn goal method f {:keys [bounds ^int N ^double n ^double jitter]
                               :or {N 100 n 0.05 jitter 0.25}
                               :as config}]
  (assert (not (nil? bounds)) "Provide search bounds.")
  (let [goal (or goal (get config :goal :minimize))
        samples (generate-points method f bounds goal N n jitter)
        opt (repeatedly (count samples) #(optimizer-fn method f config))
        ^int tk (get config :take 1)
        taker (if (> tk 1) (partial take tk) first)]
    (->> (pmap #(%1 %2) opt samples)
         (sort-by second (goal-comparator goal))
         (taker))))

#_(def scan-and-optimize (partial scan-and- optimizer))
(def scan-and-minimize (partial scan-and- minimizer :minimize))
(def scan-and-maximize (partial scan-and- maximizer :maximize))

;; bo

(defmulti ^:private utility-function (fn [t & _] t))

(defmethod utility-function :default [_ p]
  (utility-function :ucb p))

(defmethod utility-function :ucb
  [_ ^double kappa]
  (fn [gp x _]
    (let [[^double mean ^double stddev] (gp/predict gp x true)]
      (+ mean (* kappa stddev)))))

(defmethod utility-function :ei
  [_ ^double xi]
  (fn [gp x ^double y-max]
    (let [[^double mean ^double stddev] (gp/predict gp x true)
          diff (- mean y-max xi)
          z (/ diff stddev)]
      (+ (* diff ^double (r/cdf r/default-normal z))
         (* stddev ^double (r/pdf r/default-normal z))))))

(defmethod utility-function :poi
  [_ ^double xi]
  (fn [gp x ^double y-max]
    (let [[^double mean ^double stddev] (gp/predict gp x true)]
      (r/cdf r/default-normal (/ (- mean y-max xi) stddev)))))

(defn- gen-sequence
  [init-points bounds jitter]
  (let [dims (count bounds)
        int-fn (if (== dims 1)
                 #(vector (m/lerp (ffirst bounds) (second (first bounds)) %))
                 #(v/einterpolate (mapv first bounds) (mapv second bounds) %))]
    (->> (r/jittered-sequence-generator (if (< dims 15) :r2 :sobol) dims jitter)
         (take init-points)
         (map int-fn))))

(defn- initial-values
  [f init-points bounds jitter]
  (let [pts (if (sequential? init-points)
              init-points
              (gen-sequence init-points bounds jitter))]
    [pts (map f pts)]))

(defn- bayesian-step-fn
  [f util-fn warm-up bounds gp jitter optimizer optimizer-params]
  (let [params (merge optimizer-params {:N warm-up :n 0.02 :bounded? true
                                        :bounds bounds :jitter jitter})]
    (fn [{:keys [x ^double y xs ys]}]
      (let [curr-gp (gp xs ys)
            curr-util (fn [& r] (util-fn curr-gp r y))
            bx (first (scan-and-maximize optimizer curr-util params))
            ^double by (f bx)
            nxs (conj xs bx)
            nys (conj ys by)]
        {:x (if (> by y) bx x)
         :y (if (> by y) by y)
         :util-fn curr-util
         :gp curr-gp
         :xs nxs
         :ys nys
         :util-best bx}))))

(defn bayesian-optimization
  "Bayesian optimizer

  Parameters are:

  * `:warm-up` - number of brute force iterations to find maximum of utility function
  * `:init-points` - number of initial evaluation before bayesian optimization starts. Points are selected using jittered low discrepancy sequence generator (see: [[jittered-sequence-generator]]
  * `:bounds` - bounds for each dimension
  * `:utility-function-type` - one of `:ei`, `:poi` or `:ucb`
  * `:utility-param` - parameter for utility function (kappa for `ucb` and xi for `ei` and `poi`)
  * `:kernel` - kernel, default `:mattern-52`, see [[fastmath.kernel]]
  * `:kscale` - scaling factor for kernel
  * `:jitter` - jitter factor for sequence generator (used to find initial points)
  * `:noise` - noise (lambda) factor for gaussian process
  * `:optimizer` - name of optimizer (used to optimized utility function)
  * `:optimizer-params` - optional parameters for optimizer
  * `:normalize?` - normalize data in gaussian process?

  Returns lazy sequence with consecutive executions. Each step consist:

  * `:x` - maximum `x`
  * `:y` - value
  * `:xs` - list of all visited x's
  * `:ys` - list of values for every visited x
  * `:gp` - current gaussian process regression instance
  * `:util-fn` - current utility function
  * `:util-best` - best x in utility function"
  [f {:keys [warm-up init-points bounds utility-function-type utility-param kernel kscale jitter noise optimizer optimizer-params normalize?]
      :or {kscale 1.0
           kernel (k/kernel :mattern-52)
           warm-up (* ^int (count bounds) 1000)
           init-points 3
           utility-function-type :ucb
           utility-param (if (#{:ei :poi} utility-function-type) 0.001 2.576)
           jitter 0.25
           normalize? true
           noise 1.0e-8}}]
  (let [optimizer (or optimizer (if (m/one? (count bounds)) :cmaes :lbfgsb))
        f (partial apply f)
        [xs ys] (initial-values f init-points bounds jitter)
        [maxx maxy] (first (sort-by second clojure.core/> (map vector xs ys)))
        util-fn (utility-function utility-function-type utility-param)
        gp #(gp/gaussian-process %1 %2 {:normalize? normalize? :kernel kernel :kscale kscale :noise noise})
        step-fn (bayesian-step-fn f util-fn warm-up bounds gp jitter optimizer optimizer-params)]
    (rest (iterate step-fn {:x maxx
                            :y maxy
                            :xs xs
                            :ys ys}))))

;; linear optimization

(defn- constraint-relations
  ^Relationship [relation]
  (cond
    (#{:>= '>= :geq} relation) Relationship/GEQ
    (#{:<= '<= :leq} relation) Relationship/LEQ
    :else Relationship/EQ))

(defn- build-constraint
  ^LinearConstraint [[left relation right]]
  (let [relationship (constraint-relations relation)]
    (if (number? right)
      (LinearConstraint. (m/seq->double-array left) relationship (double right))
      (LinearConstraint. (m/seq->double-array (butlast left))
                         (double (last left))
                         relationship
                         (m/seq->double-array (butlast right))
                         (double (last right))))))

(defn linear-optimization
  "Solves a linear problem.

   Target is defined as a vector of coefficients and constant as the last value:
   `[a1 a2 a3 ... c]` means `f(x1,x2,x3) = a1*x1 + a2*x2 + a3*x3 + ... +  c`

   Constraints are defined as a sequence of one of the following triplets:

   * `[a1 a2 a3 ...] R n` - which means `a1*x1+a2*x2+a3*x3+... R n` 
   * `[a1 a2 a3 ... ca] R [b1 b2 b3 ... cb]` - which means `a1*x1+a2*x2+a3*x3+...+ca R b1*x1+b2*x2+b3*x3+...+cb`

   `R` is a relationship and can be one of `<=`, `>=` or `=` as symbol or keyword. Also `:leq`, `:geq` and `:eq` are valid.

  Function returns pair of optimal point and function value. If `stat?` option is set to true, returns also information about number of iterations.  

  Possible options:
  
  * `:goal` - `:minimize` (default) or `:maximize`
  * `:rule` - pivot selection rule, `:dantzig` (default) or `:bland`
  * `:max-iter` - maximum number of iterations, maximum integer by default
  * `:non-negative?` - allow non-negative variables only, default: `false` 
  * `:epsilon` - convergence value, default: `1.0e-6`:
  * `:max-ulps` - floating point comparisons, default: `10` ulp
  * `:cut-off` - pivot elements smaller than cut-off are treated as zero, default: `1.0e-10`

  ```clojure
  (linear-optimization [-1 4 0] [[-3 1] :<= 6
                                 [-1 -2] :>= -4
                                 [0 1] :>= -3])
  ;; => [(9.999999999999995 -3.0) -21.999999999999993]
  ```"
  ([target constraints] (linear-optimization target constraints {}))
  ([target constraints {:keys [goal ^double epsilon ^int max-ulps ^double cut-off
                               rule non-negative? ^int max-iter stats?]
                        :or {goal :minimize epsilon 1.0e-6 max-ulps 10 cut-off 1.0e-10
                             rule :dantzig non-negative? false max-iter Integer/MAX_VALUE}}]
   (let [goal (if (= goal :minimize) GoalType/MINIMIZE GoalType/MAXIMIZE)
         rule (if (= rule :dantzig) PivotSelectionRule/DANTZIG PivotSelectionRule/BLAND)
         max-iter (MaxIter. max-iter)
         non-negative? (NonNegativeConstraint. non-negative?)
         target (LinearObjectiveFunction. (m/seq->double-array (butlast target))
                                          (double (last target)))
         constraints (->> constraints
                          (partition 3)
                          ^java.util.Collection (map build-constraint)
                          (LinearConstraintSet.))
         ^BaseOptimizer solver (SimplexSolver. epsilon max-ulps cut-off)]
     (->> [goal rule max-iter non-negative? target constraints]
          (into-array OptimizationData)
          (.optimize solver)
          (parse-result nil)
          (maybe-stats? stats? solver)))))

(m/unuse-primitive-operators)


#_(let [f (fn [^double x ^double y] (inc (- (- (* x x)) (m/sq (dec y)))))
        bounds [[-4 4] [-3 3]]
        bo (bayesian-optimization f {:bounds bounds
                                     :utility-function-type :ei
                                     :utility-param 0.1
                                     :optimizer :powell})]
    (println (f 0 1))
    (last (take 30 (map (juxt :x :y) bo)))    )

#_(let [f (fn [^double x] (- (+ (/ (m/sin (* 10 m/PI x)) (+ x x)) (m/pow (dec x) 4))))
        bounds [[-0.2 0.2]]
        bo (bayesian-optimization f {:bounds bounds
                                     :utility-function-type :ucb
                                     ;; :utility-param 0
                                     :optimizer :lbfgsb})]
    (take 3 (drop 30 (map (juxt :x :y) bo))))
;; => ([(0.14580050823621077) 2.867142408546129] [(0.14564645826782044) 2.868128502656343] [(0.14550225115172954) 2.868981775657352])
;; => ([(0.14973586005947528) 2.8164430833342324] [(0.14973586005947528) 2.8164430833342324] [(0.14973586005947528) 2.8164430833342324])

;; tests

#_(defn target-1d ^double [^double x]
    (+ (/ (m/sin (* 10 m/PI x)) (+ x x)) (m/pow (dec x) 4)))

#_(defn target-2d-schwefel
    ^double [^double x ^double y]
    (- 418.9829
       (+ (* x (m/sin (m/sqrt (m/abs x))))
          (* y (m/sin (m/sqrt (m/abs y)))))))

#_(defn target-2d-booth
    ^double [^double x ^double y]
    (+ (m/sq (+ x y y -7))
       (m/sq (+ x x y -5))))

#_(time (let [f (minimizer :lbfgsb target-2d-booth {:bounds [[-10 10]
                                                             [-10 10]]
                                                    :tolerance 1.0e-10})]
          (f (v/generate-vec2 #(r/drand -9 9)))))

#_(scan-and-minimize :lbfgsb target-2d-schwefel {:bounds [[-500 500]
                                                          [-500 500]] :N 1000 :bounded? true})

#_(scan-and-minimize :gradient target-1d {:bounds [[-0.5 2.5]] :gradient-step 0.1})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#_(do

    (defn bfn6 ^double [^double x ^double y] (+ (* 100 (m/sqrt (m/abs (- y (* 0.01 x x)))))
                                             (* 0.01 (m/abs (+ 10.0 x)))))

    (defn h ^double [^double x ^double y] (+ (m/sq (+ (* x x) y -11))
                                          (m/sq (+ x (* y y) -7))))

    (defn d5 [a b c d e] (reduce #(+ ^double %1 (m/sq %2)) 0.0 [a b c d e]))

    (time (scan-and-maximize :bobyqa bfn6 {:bounds [[-15 -3] [15 3]] :N 100 :n 0.2}))

    (minimize :powell bfn6 {:bounds [[-15 -3] [15 3]] :initial [-13.217532309719662 0.9280007414397415]})

    (time (scan-and-optimize :powell #(m/cos %) {:bounds [-3 3] :initial -2 :N 100 :n 0.2 :goal :maximize})))


;; => [1.9210981963566007 -5.751481824637489 0.3304425131054902]

#_(defn rosenbrock
    [& vs]
    (reduce (fn [^double fx [^double xi ^double xi+1]]
              (let [t1 (- 1.0 xi)
                    t2 (* 10.0 (- xi+1 (* xi xi)))]
                (+ fx (* t1 t1) (* t2 t2)))) 0.0 (partition 2 1 vs)))


#_(minimize :lbfgsb rosenbrock {:bounds (repeat 20 [-5 10])
                                :init [2 -4 2 4 -2] :m 50 :N 10 :n 1
                                :max-iters 1000})



#_(comment (defn hump
             [^double x ^double y]
             (let [x2 (* x x)
                   x4 (* x2 x2)]
               (+ (* 2.0 x2)
                  (* -1.05 x4)
                  (* x4 x2 m/SIXTH)
                  (* x y)
                  (* y y))))

           (defn hump-grad
             [^double x ^double y]
             (let [x2 (* x x)
                   x4 (* x2 x2)]
               [(+ (* 4.0 x)
                   (* -4.2 x2 x)
                   (* x4 x)
                   y)
                (+ x (* 2.0 y))]))

           (minimize :lbfgsb hump {:bounds [[-5 5.1] [-5 5.1]]
                                   :gradient-f hump-grad :N 1000
                                   :weak-wolfe? false
                                   :stats? true}))

;; sudoku solver

#_(defn sudoku-idx ^long [^long v ^long r ^long c] (+ (* 81 v) (* 9 r) c))

#_(defn sudoku-solver
    [input]
    (let [zeros (vec (repeat (* 9 9 9) 0))
          c1 (for [r (range 9)
                   c (range 9)]
               (reduce (fn [z v]
                         (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
          c2 (for [v (range 9)
                   c (range 9)]
               (reduce (fn [z r]
                         (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
          c3 (for [v (range 9)
                   r (range 9)]
               (reduce (fn [z c]
                         (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
          c4 (for [v (range 9)
                   p (range 3)
                   q (range 3)
                   :let [rc (for [r (range (* 3 p) (* 3 (inc p)))
                                  c (range (* 3 q) (* 3 (inc q)))]
                              [r c])]]
               (reduce (fn [z [r c]]
                         (assoc z (sudoku-idx v r c) 1)) zeros rc))
          c5 (for [[v r c] input]
               (assoc zeros (sudoku-idx (dec v) r c) 1))
          all (reduce (fn [buff c]
                        (conj buff c :eq 1)) [] (mapcat identity [c1 c2 c3 c4 c5]))
          res (vec (first (linear-optimization (conj zeros 0) all {:non-negative? true})))
          solution (vec (repeat 81 0))]
      (->> (for [r (range 9)
                 c (range 9)
                 v (range 9)
                 :let [id (sudoku-idx v r c)]
                 :when (m/one? (m/round (res id)))]
             [v r c])
           (reduce (fn [s [v r c]]
                     (assoc s (+ c (* r 9)) (inc v))) solution)
           (partition 9))))

;; v - value 1-9
;; r - row 0-8
;; c - column 0-8
;;
;;               v r c
#_(sudoku-solver [[8 0 2]
                  [4 0 4]
                  [2 0 5]
                  [3 0 8]

                  [6 1 1]
                  [1 1 7]
                  
                  [3 2 1]
                  [7 2 5]

                  [5 3 5]
                  [3 3 6]

                  [4 4 1]
                  [9 4 3]
                  [1 4 5]
                  [7 4 7]

                  [5 5 2]
                  [2 5 3]

                  [3 6 3]
                  [4 6 7]

                  [5 7 1]
                  [6 7 7]

                  [7 8 0]
                  [6 8 3]
                  [1 8 4]
                  [9 8 6]])

;; => ((9 7 8 1 4 2 6 5 3)
;;     (2 6 4 5 9 3 7 1 8)
;;     (5 3 1 8 6 7 4 2 9)
;;     (1 2 7 4 8 5 3 9 6)
;;     (8 4 6 9 3 1 5 7 2)
;;     (3 9 5 2 7 6 1 8 4)
;;     (6 1 9 3 5 8 2 4 7)
;;     (4 5 3 7 2 9 8 6 1)
;;     (7 8 2 6 1 4 9 3 5))
