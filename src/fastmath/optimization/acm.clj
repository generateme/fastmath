(ns fastmath.optimization.acm
  "Apache Commons Math optimizers"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.finite :as finite])
  (:import [org.apache.commons.math3.optim.univariate BrentOptimizer UnivariatePointValuePair SearchInterval UnivariateObjectiveFunction BracketFinder]
           [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunctionGradient ObjectiveFunction]
           [org.apache.commons.math3.optim MaxEval MaxIter OptimizationData SimpleBounds InitialGuess BaseOptimizer PointValuePair SimpleValueChecker]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer CMAESOptimizer CMAESOptimizer$PopulationSize CMAESOptimizer$Sigma NelderMeadSimplex MultiDirectionalSimplex SimplexOptimizer PowellOptimizer]
           [org.apache.commons.math3.optim.nonlinear.scalar.gradient NonLinearConjugateGradientOptimizer NonLinearConjugateGradientOptimizer$Formula NonLinearConjugateGradientOptimizer$IdentityPreconditioner Preconditioner]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction MultivariateVectorFunction]

           [org.apache.commons.math3.random RandomGenerator JDKRandomGenerator]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; common

(defn goal-type
  [goal]
  (if-not goal
    GoalType/MINIMIZE
    (case goal
      :minimize GoalType/MINIMIZE
      :maximize GoalType/MAXIMIZE)))

(defn max-eval [max-evals] (MaxEval. (or max-evals 10000)))
(defn max-iter [max-iters] (MaxIter. (or max-iters 10000)))

(defn- optimization-data ^"[Lorg.apache.commons.math3.optim.OptimizationData;" [data] (into-array OptimizationData data))

;; univariate

(defn- search-interval
  ([^double low ^double high]
   (SearchInterval. low high))
  ([^double low ^double high init]
   (if init
     (SearchInterval. low high (if (sequential? init) (first init) init))
     (search-interval low high))))

(defn- univariate-function
  ^UnivariateFunction [f vector-arg?]
  (if vector-arg?
    (reify UnivariateFunction (value [_ x] (f [x])))
    (reify UnivariateFunction (value [_ x] (f x)))))

(defn- univariate-objective-function ^UnivariateObjectiveFunction [^UnivariateFunction f] (UnivariateObjectiveFunction. f))

(defn- process-bounds
  [bounds find-bracket ^UnivariateFunction uf goal]
  (let [bounds (if (sequential? (first bounds)) (first bounds) bounds)]
    (if find-bracket
      (let [{:keys [^double grow-limit ^long max-evals]
             :or {grow-limit 100.0 max-evals 500}} (if (map? find-bracket) find-bracket {})
            [^double x ^double y] bounds
            ^BracketFinder bf (BracketFinder. grow-limit max-evals)]
        (.search bf uf goal x y)
        [(.getLo bf) (.getHi bf)])
      bounds)))

(defn- parse-result
  [^UnivariatePointValuePair res vector-arg?]
  (if vector-arg?
    [[(.getPoint res)] (.getValue res)]
    [(.getPoint res) (.getValue res)]))

(defn brent-data
  [f {:keys [bounds initial max-evals max-iters ^double rel ^double abs find-bracket vector-arg? goal]
      :or {rel 1.0e-8 abs 1.0e-10}}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [goal (goal-type goal)
        uf (univariate-function f vector-arg?)
        uof (univariate-objective-function uf)
        [^double lo ^double high] (process-bounds bounds find-bracket uf goal)
        si (search-interval lo high initial)
        me (max-eval max-evals)
        mi (max-iter max-iters)]
    {:goal goal
     :univariate-objective-function uof
     :search-interval si
     :max-eval me
     :max-iter mi
     :rel rel :abs abs
     :vector-arg? vector-arg?}))

(defn update-univariate-initial
  [{:keys [^SearchInterval search-interval] :as brent-data} initial]
  (let [initial (if (sequential? initial) (first initial) initial)]
    (assoc brent-data :search-interval (SearchInterval. (.getMin search-interval)
                                                        (.getMax search-interval)
                                                        initial))))

(def ^:private univariate-extract-optimization-data (juxt :goal :univariate-objective-function :search-interval :max-eval :max-iter))

(defn brent
  ([f opts] (brent (brent-data f opts)))
  ([{:keys [^double rel ^double abs vector-arg?] :as brent-data}]
   (let [^BrentOptimizer bo (BrentOptimizer. rel abs)]
     (-> (univariate-extract-optimization-data brent-data)
         (optimization-data)
         (->> (.optimize bo))
         (parse-result vector-arg?)))))

;; multivariate

(defn- multivariate-function
  ^MultivariateFunction [f vector-arg?]
  (if vector-arg?
    (reify MultivariateFunction (value [_ x] (f x)))
    (reify MultivariateFunction (value [_ x] (apply f x)))))

(defn- objective-function ^ObjectiveFunction [^MultivariateFunction mf] (ObjectiveFunction. mf))

(defn- multivariate-bounds
  ^SimpleBounds [bounds]
  (when bounds
    (SimpleBounds. (double-array (map first bounds))
                   (double-array (map second bounds)))))

(defn- multivariate-initial
  ^InitialGuess [^SimpleBounds bounds initial]
  (if initial
    (InitialGuess. (m/seq->double-array initial))
    (if bounds
      (InitialGuess. (v/interpolate (.getLower bounds) (.getUpper bounds) 0.5))
      (throw (ex-info "No initial point provided or inferred" {:bounds bounds :initial initial})))))

(defn- multivariate-common-data
  [f {:keys [bounds initial max-evals max-iters vector-arg? goal stats?]
      :or {vector-arg? true stats? false}}]
  (let [bounds (multivariate-bounds bounds)
        initial (multivariate-initial bounds initial)
        opts {:goal (goal-type goal)
              :max-eval (max-eval max-evals)
              :max-iter (max-iter max-iters)
              :vector-arg? vector-arg?
              :objective-function (objective-function (multivariate-function f vector-arg?))
              :stats? stats?
              :initial initial}]
    (cond-> opts
      bounds (assoc :bounds bounds))))

(defn update-multivariate-initial
  [data initial]
  (assoc data :initial (InitialGuess. (m/seq->double-array initial))))

(def ^:private multivariate-optimization-data-keys #{:goal :max-eval :max-iter :objective-function :bounds :initial})

(defn multivariate-optimize
  [^BaseOptimizer opt data stats?]
  (let [^PointValuePair res (->> (optimization-data data)
                                 (.optimize opt))
        point (vec (.getPointRef res))
        value (.getValue res)]
    (if stats?
      {:point point
       :value value
       :evaluations (.getEvaluations opt)
       :iterations (.getIterations opt)}
      [point value])))

;;

(defn- bobyqa-radius
  ^double [bounds]
  (* 0.5 (v/mx (map (fn [[^double l ^double u]] (m/- u l)) bounds))))

(defn bobyqa-data
  [f {:keys [bounds number-of-points initial-radius ^double stopping-radius]
      :or {initial-radius :inferred
           stopping-radius BOBYQAOptimizer/DEFAULT_STOPPING_RADIUS}
      :as opts}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [n (count bounds)]
    (when (m/< n 2) (throw (ex-info "Number of dimensions should be equal or greater than 2" {:n n :bounds bounds})))
    (let [number-of-points (long (or number-of-points (m/inc (m/* 2 n))))]
      (assoc (multivariate-common-data f opts)
             :initial-radius (condp = initial-radius
                               :default BOBYQAOptimizer/DEFAULT_INITIAL_RADIUS
                               :inferred (bobyqa-radius bounds))
             :stopping-radius stopping-radius
             :number-of-points number-of-points))))

(defn bobyqa
  ([f opts] (bobyqa (bobyqa-data f opts)))
  ([{:keys [^long number-of-points ^double initial-radius ^double stopping-radius stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys)
         ^BOBYQAOptimizer bo (BOBYQAOptimizer. number-of-points initial-radius stopping-radius)]
     (multivariate-optimize bo opt-data stats?))))

;;

(defn- bounds->steps
  ^doubles [bounds ^double length]
  (double-array (map (fn [[^double low ^double high]]
                       (m/* length (m/- high low))) bounds)))

(defn cmaes-data
  [f {:keys [bounds ^double stop-fitness active-cma? ^long diagonal-only ^long check-feasible-count rng ^double rel ^double abs
             population-size ^double sigma]
      :or {stop-fitness 1.0e-10 active-cma? true diagonal-only 0 check-feasible-count 0 rng (JDKRandomGenerator.)
           rel 1.0e-10 abs 1.0e-10 sigma 0.2}
      :as opts}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [dims (count bounds)
        population-size (CMAESOptimizer$PopulationSize. (long (or population-size (m/ceil (m/+ 4 (m/* 3 (m/log dims)))))))
        sigma (CMAESOptimizer$Sigma. (bounds->steps bounds sigma))
        checker (SimpleValueChecker. rel abs)]
    (assoc (multivariate-common-data f opts)
           :stop-fitness stop-fitness
           :active-cma? (boolean active-cma?)
           :diagonal-only diagonal-only
           :check-feasible-count check-feasible-count
           :rng rng
           :population-size population-size
           :sigma sigma
           :checker checker)))

(def ^:private multivariate-optimization-data-keys-cmaes (conj multivariate-optimization-data-keys :population-size :sigma))

(defn cmaes
  ([f opts] (cmaes (cmaes-data f opts)))
  ([{:keys [^MaxIter max-iter ^double stop-fitness ^boolean active-cma? ^long diagonal-only ^long check-feasible-count ^RandomGenerator rng
            checker stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys-cmaes)
         ^CMAESOptimizer bo (CMAESOptimizer. (.getMaxIter max-iter) stop-fitness active-cma? diagonal-only check-feasible-count rng false checker)]
     (multivariate-optimize bo opt-data stats?))))

;;

(def ^:private multivariate-optimization-data-keys-simplex (-> multivariate-optimization-data-keys (conj :simplex) (disj :bounds)))

(defn- simplex-data
  [f {:keys [bounds initial length ^double rel ^double abs]
      :or {rel 1.0e-10 abs 1.0e-10}
      :as opts}]
  (let [steps (if bounds
                (bounds->steps bounds (double (or length 0.4)))
                (do (when-not initial (throw (ex-info "Bounds or initial should be provided." {:bounds bounds :initial initial})))
                    (double-array (repeat (count initial) (double (or length 1.0))))))]
    (assoc (multivariate-common-data f opts)
           :steps steps
           :rel rel
           :abs abs)))

(defn nelder-mead-data
  [f {:keys [^double rho ^double khi ^double gamma ^double sigma]
      :or {rho 1.0 khi 2.0 gamma 0.5 sigma 0.5}
      :as opts}]
  (let [{:keys [^doubles steps] :as opts} (simplex-data f opts)]
    (assoc opts :simplex (NelderMeadSimplex. steps rho khi gamma sigma))))

(defn nelder-mead
  ([f opts] (nelder-mead (nelder-mead-data f opts)))
  ([{:keys [^double rel ^double abs stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys-simplex)
         ^SimplexOptimizer so (SimplexOptimizer. rel abs)]
     (multivariate-optimize so opt-data stats?))))

(defn multidirectional-simplex-data
  [f {:keys [^double khi ^double gamma]
      :or {khi 2.0 gamma 0.5}
      :as opts}]
  (let [{:keys [^doubles steps] :as opts} (simplex-data f opts)]
    (assoc opts :simplex (MultiDirectionalSimplex. steps khi gamma))))

(defn multidirectional-simplex
  ([f opts] (multidirectional-simplex (multidirectional-simplex-data f opts)))
  ([{:keys [^double rel ^double abs stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys-simplex)
         ^SimplexOptimizer so (SimplexOptimizer. rel abs)]
     (multivariate-optimize so opt-data stats?))))

;;

(defn powell-data
  [f {:keys [^double rel ^double abs line-rel line-abs]
      :or {rel 1.0e-10 abs 1.0e-10}
      :as opts}]
  (assoc (multivariate-common-data f opts)
         :rel rel
         :abs abs
         :line-rel (or line-rel (m/sqrt rel))
         :line-abs (or line-abs (m/sqrt abs))))

(def ^:private multivariate-optimization-data-keys-powell (-> multivariate-optimization-data-keys (disj :bounds)))

(defn powell
  ([f opts] (powell (powell-data f opts)))
  ([{:keys [^double rel ^double abs ^double line-rel ^double line-abs stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys-powell)
         ^PowellOptimizer po (PowellOptimizer. rel abs line-rel line-abs)]
     (multivariate-optimize po opt-data stats?))))

;;

(defn- multivariate-function-gradient
  ^MultivariateVectorFunction [f]
  (reify MultivariateVectorFunction (value [_ x] (m/seq->double-array (f x)))))

(defn- objective-function-gradient ^ObjectiveFunctionGradient [^MultivariateVectorFunction mfg] (ObjectiveFunctionGradient. mfg))

(defn- hessian-preconditioner
  [f ^double h]
  (let [hd (finite/hessian-diagonal f {:h h})]
    (reify Preconditioner
      (precondition [_ point r]
        (let [^doubles diagonal (hd point)
              ^doubles nr (double-array r)]
          (doseq [^long i (range (alength nr))
                  :let [v (Array/aget diagonal i)]
                  :when (m/> v 1.0e-6)]
            (Array/aset nr i (m// (Array/aget r i) v)))
          nr)))))

(defn non-linear-gradient-data
  ([f {:keys [^double gradient-h ^long gradient-acc]
       :or {gradient-h 1.0e-6 gradient-acc 2}
       :as opts}]
   (non-linear-gradient-data f (finite/gradient f {:h gradient-h :acc gradient-acc}) opts))
  ([f gradient {:keys [^double rel ^double abs line-rel line-abs formula ^double bracketing-range preconditioner ^double hessian-h]
                :or {rel 1.0e-10 abs 1.0e-10 formula :polak-ribiere bracketing-range 1.0e-10 preconditioner :identity hessian-h 5.0e-3}
                :as opts}]
   (assoc (multivariate-common-data f opts)
          :rel rel
          :abs abs
          :line-rel (or line-rel (m/sqrt rel))
          :line-abs (or line-abs (m/sqrt abs))
          :bracketing-range bracketing-range
          :checker (SimpleValueChecker. rel abs)
          :formula (case formula
                     :polak-ribiere NonLinearConjugateGradientOptimizer$Formula/POLAK_RIBIERE
                     :fletcher-reeves NonLinearConjugateGradientOptimizer$Formula/FLETCHER_REEVES)
          :preconditioner (case preconditioner
                            :identity (NonLinearConjugateGradientOptimizer$IdentityPreconditioner.)
                            :hessian (hessian-preconditioner f hessian-h))
          :objective-function-gradient (objective-function-gradient (multivariate-function-gradient gradient)))))

(def ^:private multivariate-optimization-data-keys-non-linear-gradient (-> multivariate-optimization-data-keys
                                                                           (conj :objective-function-gradient)
                                                                           (disj :bounds)))

(defn non-linear-gradient
  ([f gradient opts] (non-linear-gradient (non-linear-gradient-data f gradient opts)))
  ([f opts] (non-linear-gradient (non-linear-gradient-data f opts)))
  ([{:keys [^SimpleValueChecker checker ^double line-rel ^double line-abs ^double bracketing-range formula ^Preconditioner preconditioner stats?]
     :as data}]
   (let [opt-data (map data multivariate-optimization-data-keys-non-linear-gradient)
         
         ^NonLinearConjugateGradientOptimizer nlcgo (NonLinearConjugateGradientOptimizer. formula checker line-rel line-abs bracketing-range preconditioner)]
     (multivariate-optimize nlcgo opt-data stats?))))


