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
  * Bayesian Optimization (see below)

  All optimizers require bounds.

  ## Optimizers

  To optimize functions call one of the following functions:

  * [[minimize]] or [[maximize]] - to perform actual optimization
  * [[scan-and-minimize]] or [[scan-and-maximize]] - functions find initial point using brute force and then perform optimization paralelly for best initialization points. Brute force scan is done using jitter low discrepancy sequence generator.

  You can also create optimizer (function which performs optimization) by calling [[minimizer]] or [[maximizer]]. Optimizer accepts initial point.

  All above accept:

  * one of the optimization method, ie: `:brent`, `:bobyqa`, `:nelder-mead`, `:multidirectional-simplex`, `:cmaes`, `:gradient`
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

  Bayesian optimizer can be used for optimizing expensive to evaluate black box functions. Refer this [article](http://krasserm.github.io/2018/03/21/bayesian-optimization/) or this [article](https://nextjournal.com/a/LKqpdDdxiggRyHhqDG5FH?token=Ss1Qq3MzHWN8ZyEt9UC1ZZ)"
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunction ObjectiveFunctionGradient]
           [org.apache.commons.math3.optim.univariate SearchInterval BrentOptimizer UnivariateObjectiveFunction UnivariatePointValuePair]
           [org.apache.commons.math3.optim BaseOptimizer OptimizationData MaxEval MaxIter SimpleBounds SimpleValueChecker InitialGuess PointValuePair]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction MultivariateVectorFunction]
           [org.apache.commons.math3.optim.nonlinear.scalar MultivariateFunctionMappingAdapter]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer PowellOptimizer NelderMeadSimplex SimplexOptimizer MultiDirectionalSimplex CMAESOptimizer CMAESOptimizer$PopulationSize CMAESOptimizer$Sigma]
           [org.apache.commons.math3.optim.nonlinear.scalar.gradient NonLinearConjugateGradientOptimizer NonLinearConjugateGradientOptimizer$Formula]
           [smile.math BFGS DifferentiableMultivariateFunction])
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.kernel :as k]
            [fastmath.gp :as gp]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def ^:private univariate-set #{:brent})
(def ^:private multivariate-set #{:bobyqa :powell :nelder-mead :multidirectional-simplex :cmaes :gradient :bfgs})
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
  (double-array (map-indexed (fn [^long id ^double x]
                               (let [x+ (+ x step)
                                     x- (- x step)
                                     step2 (+ (- x+ x) (- x x-))
                                     a (aclone ^doubles xs)
                                     v1 (do (aset a id (+ x step))
                                            (.value f a))
                                     v2 (do (aset a id (- x step))
                                            (.value f a))]
                                 (/ (- v1 v2) step2))) xs)))

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
  [^MultivariateFunctionMappingAdapter mfma res]
  (condp instance? res
    UnivariatePointValuePair (let [^UnivariatePointValuePair res res]
                               [(list (.getPoint res)) (.getValue res)])
    PointValuePair (let [^PointValuePair res res]
                     [(seq (if mfma
                             (.unboundedToBounded mfma (.getPointRef res))
                             (.getPointRef res))) (.getValue res)])
    res))

(defn- find-dimensions
  ^long [bounds]
  (if (sequential? (first bounds)) (count bounds) 1))

(defn- fix-brent-bounds
  [method bounds]
  (if (and (= method :brent)
           (sequential? (first bounds)))
    (first bounds) bounds))

(defn- smile-fn
  [^MultivariateFunction mf ^double gradient-h]
  (reify DifferentiableMultivariateFunction
    (f [_ xs] (.value mf xs))
    (g [_ xs gout]
      (let [^doubles gradient (finite-differences mf xs gradient-h)]
        (System/arraycopy gradient 0 gout 0 (alength gradient))
        (.value mf xs)))))

(defn- bfgs
  [f {:keys [^int max-iters ^int m ^double tolerance goal bounds bounded? ^double gradient-h]
      :or {max-iters 1000 m 5 tolerance 1.0e-8 goal :minimize bounded? true gradient-h 0.0001}}]
  (let [dim (find-dimensions bounds)
        mf (if (= goal :minimize)
             (multivariate-function f)
             (negative-multivariate-function f))
        ^DifferentiableMultivariateFunction sf (smile-fn mf gradient-h)]
    (if bounded?
      (let [l (double-array (map first bounds))
            u (double-array (map second bounds))]
        (fn smile-bounded
          ([] (smile-bounded nil))
          ([init]
           (let [i (double-array (or init (mid-point bounds)))
                 res (BFGS/minimize sf m i l u tolerance max-iters)]
             [(seq i) (if (= goal :minimize) res (- res))]))))
      (fn smile-unbounded
        ([] (smile-unbounded nil))
        ([init]
         (let [i (double-array (or init (repeat dim 0.0)))
               res (BFGS/minimize sf m i tolerance max-iters)]
           [(seq i) (if (= goal :minimize) res (- res))]))))))

(defn- optimizer
  [method f {:keys [max-evals max-iters goal bounds stats? population-size bounded? gradient-h]
             :or {gradient-h 0.0001}
             :as config}]
  (if (= method :bfgs)
    (bfgs f config)
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
          ([init] (let [res (->> (conj base-opt-data (wrap-initial method bounds init mfma))
                                 (into-array OptimizationData)
                                 (.optimize optimizer)
                                 (parse-result mfma))]
                    (if-not stats?
                      res
                      {:result res
                       :evaluations (.getEvaluations optimizer)
                       :iterations (.getIterations optimizer)}))))))))

(defn minimizer
  "Create optimizer which minimizes function.

  Returns function which performs optimization for optionally given initial point."
  [method f config] (optimizer method f (assoc config :goal :minimize)))

(defn maximizer
  "Create optimizer which maximizer function.

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
        N (max 10 ^int N)
        n (max 4 (m/floor (* ^double n N)))
        gen (r/jittered-sequence-generator (if (< dim 5) :r2 :sobol) dim jitter)]
    (->> (take N (if (and (not= method :brent)
                          (m/one? dim)) (map vector gen) gen))
         (map #(let [p (inter lo high %)]
                 [(genf p) p]))
         (sort-by first (goal-comparator goal))
         (take n)
         (map second))))

(defn- scan-and-
  "For cheap functions, scan domain and bruteforcely search for set of minimal values, then use part of this values as initial points for parallel optimization.

  Additional parameters in config:

  * N - number of points to scan per dimension (default: 100, minumum 10)
  * n - fraction of total points N) used for optimization (default: 0.05, minimum 10)"
  [optimizer goal method f {:keys [bounds ^int N ^double n ^double jitter]
                            :or {N 100 n 0.05 jitter 0.25}
                            :as config}]
  (assert (not (nil? bounds)) "Provide search bounds.")
  (let [goal (or goal (get config :goal :minimize))
        samples (generate-points method f bounds goal N n jitter)
        opt (repeatedly (count samples) #(optimizer method f config))
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
  (let [optimizer (or optimizer (if (m/one? (count bounds)) :cmaes :bfgs))
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

;; root finding methods


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
                                     :optimizer :bfgs})]
    (take 3 (drop 30 (map (juxt :x :y) bo))))

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

#_(time (let [f (minimizer :bfgs target-2d-booth {:bounds [[-10 10]
                                                           [-10 10]]})]
          (f (v/generate-vec2 #(r/drand -9 9)))))

#_(scan-and-minimize :cmaes target-2d-schweel {:bounds [[-500 500]
                                                        [-500 500]] :N 100 :bounded? true})

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
