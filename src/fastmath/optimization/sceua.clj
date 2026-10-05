(ns fastmath.optimization.sceua
  "Shuffled Complex Evolution (SCE-UA), a global optimizer for box constrained problems.

  The method keeps a population of points sampled from the bounds. In every loop the population is sorted, divided into complexes and every complex is evolved independently by the competitive complex evolution: a few points are drawn from the complex, with a probability higher for better points, and the worst of them is replaced by a reflection or a contraction of the centroid of the others, or by a random point when both fail. Then the complexes are merged and shuffled again. The shuffling spreads the information found by every complex over the whole population. Optionally the number of complexes shrinks during the run, and the population is repaired when it collapses into a subspace (see the option `:pca-recovery?` of [[sceua]]).

  The optimizer is a function of the objective and an options map, `(sceua f opts)`. A call keeps its state to itself, so the function can be called from many threads.

  The algorithm is described in Duan, Sorooshian and Gupta (1992), Effective and efficient global optimization for conceptual rainfall-runoff models, Water Resources Research 28(4), 1015-1031, and in Duan, Sorooshian and Gupta (1994), A shuffled complex evolution approach for effective and efficient global minimization, Journal of Optimization Theory and Applications 76(3), 501-521. The repair of a collapsed population follows the idea of Chu, Gao and Sorooshian (2010), Water Resources Research 46(9).

  Functions:

  - [[sceua]] - the optimizer."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.matrix :as mat]
            [fastmath.vector :as v]
            [fastmath.optimization.common :as common])
  (:import [java.util.concurrent ExecutionException]
           [org.apache.commons.math3.random RandomDataGenerator RandomGenerator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; A point of the population is a map `{:x coordinates :value value}`: the coordinates are a double array which is
;; never changed, the value is the value of the minimized objective. A population and a complex are vectors of points
;; sorted by the value. The settings are the validated options with the key `:evaluate` added (see `evaluator`).

(def ^:private ^:const max-dimensions 1000)
(def ^:private ^:const pca-ratio 1.0e-3)
(def ^:private ^:const pca-fraction 0.1)
(def ^:private ^:const pca-step 0.05)

(defn- parse-initial
  "Validates the initial point and returns it as an array or `nil`."
  [initial lo hi]
  (when (some? initial)
    (let [xs (if (number? initial) [initial] (seq initial))
          inside? (fn [x l h]
                    (let [x (double x)]
                      (and (m/valid-double? x) (m/<= (double l) x (double h)))))]
      (when-not (every? number? xs)
        (common/throw-invalid-option :initial initial "must be a number or a sequence of numbers"))
      (when-not (every? true? (map inside? xs lo hi))
        (common/throw-invalid-option :initial initial "must be finite and lie within the bounds"))
      (double-array xs))))

(defn- parse-options
  "Validates the options of [[sceua]] and returns them as a map with the defaults applied."
  [opts]
  (let [opts (or opts {})
        {:keys [bounds initial goal rng parallel? pca-recovery? stats?]} opts
        bounds (common/normalize-bounds :sceua bounds initial)
        n (count bounds)
        _ (when (m/> n max-dimensions)
            (common/throw-invalid-option :bounds bounds (str "at most " max-dimensions " dimensions are supported")))
        lo (double-array (map first bounds))
        hi (double-array (map second bounds))
        complexes (common/positive-integer-option opts :complexes 5)
        complex-size (common/positive-integer-option opts :complex-size (m/inc (m/* 2 n)))
        subcomplex-size (common/positive-integer-option opts :subcomplex-size (m/min (m/inc n) complex-size))
        min-complexes (common/positive-integer-option opts :min-complexes complexes)
        jitter (common/nonnegative-number-option opts :jitter 0.25)]
    (when (m/< complex-size 2)
      (common/throw-invalid-option :complex-size complex-size "must be at least 2"))
    (when-not (m/<= 2 subcomplex-size complex-size)
      (common/throw-invalid-option :subcomplex-size subcomplex-size "must be at least 2 and not greater than :complex-size"))
    (when (m/> min-complexes complexes)
      (common/throw-invalid-option :min-complexes min-complexes "must not be greater than :complexes"))
    (when (m/> jitter 1.0)
      (common/throw-invalid-option :jitter jitter "must not be greater than 1"))
    {:lo lo
     :hi hi
     :initial (parse-initial initial lo hi)
     :sign (if (= :maximize (common/parse-goal goal)) -1.0 1.0)
     :complexes complexes
     :complex-size complex-size
     :subcomplex-size subcomplex-size
     :evolution-steps (common/positive-integer-option opts :evolution-steps complex-size)
     :min-complexes min-complexes
     :stop-loops (common/positive-integer-option opts :stop-loops 40)
     :stop-improvement (common/nonnegative-number-option opts :stop-improvement 1.0e-5)
     :stop-range (common/nonnegative-number-option opts :stop-range 1.0e-3)
     :max-evals (common/positive-integer-option opts :max-evals 10000)
     :max-iters (common/positive-integer-option opts :max-iters 10000)
     :jitter jitter
     :rng (r/ensure-rng rng)
     :parallel? (boolean parallel?)
     :pca-recovery? (boolean pca-recovery?)
     :stats? (boolean stats?)}))

;; points

(defn- value-of
  "Value of the objective stored in a point."
  ^double [point]
  (:value point))

(defn- sorted-by-value
  "Returns the points as a vector sorted by the value, the best first. The order of equal values is kept."
  [points]
  (vec (sort-by :value points)))

(defn- clip
  "Moves the coordinates `x` into the box `[lo, hi]`."
  [x lo hi]
  (v/emn (v/emx x lo) hi))

(defn- evaluator
  "Creates the function which evaluates the objective at a point.

  Called with coordinates (a double array) it returns a new point. The value of the point is `sign` times the value of `f`, `NaN` becomes `+Infinity`, the worst value. The call number `max-evals` plus one throws `ex-info`. Called without arguments it returns the number of the calls made so far.

  The counter is atomic, so the function can be called from many threads. `f` receives the array itself and must not change it."
  [f ^double sign ^long max-evals]
  (let [calls (atom 0)]
    (fn
      ([] @calls)
      ([x]
       (when (m/> (long (swap! calls inc)) max-evals)
         (throw (ex-info "Maximum number of evaluations exceeded" {:max-evals max-evals :evaluations max-evals})))
       (let [value (m/* sign (double (f x)))]
         {:x x :value (if (m/nan? value) ##Inf value)})))))

(defn- initial-population
  "Creates the sorted population: the initial point, if given, and points of a jittered low discrepancy sequence scaled to the bounds."
  [{:keys [lo hi initial evaluate jitter rng ^long complexes ^long complex-size]}]
  (let [n (count lo)
        generated (long (m/- (m/* complexes complex-size) (if initial 1 0)))
        scale (fn [u] (v/einterpolate lo hi (v/vec->array u)))
        points (map (comp evaluate scale)
                    (take generated (r/jittered-sequence-generator (if (m/< n 15) :r2 :sobol) n jitter rng)))]
    (sorted-by-value (if initial (conj points (evaluate initial)) points))))

;; competitive complex evolution

(defn- triangular-index
  "Zero-based index drawn with the probability proportional to `m - index` from a uniform `u` in [0,1)."
  ^long [^long m ^double u]
  (let [mh (m/+ m 0.5)
        root (m/safe-sqrt (m/- (m/* mh mh) (m/* m (m/inc m) u)))
        index (long (m/floor (m/- mh root)))]
    (m/constrain index 0 (m/dec m))))

(defn- select-parents
  "Returns a sorted vector of `q` different indices drawn from `m` ones with the triangular probability (the lower the index, the more probable).

  An index which is drawn again is skipped, so every next one is drawn from the remaining ones."
  [^long m ^long q rng]
  (loop [chosen (sorted-set)]
    (if (m/== (count chosen) q)
      (vec chosen)
      (recur (conj chosen (triangular-index m (r/drandom rng)))))))

(defn- random-point
  "Returns a random point of the box `[lo, hi]`."
  [lo hi rng]
  (v/einterpolate lo hi (double-array (repeatedly (count lo) #(r/drandom rng)))))

(defn- evolve-complex
  "Evolves a complex (a vector of points sorted by value) and returns a new sorted vector of the same size.

  Every step draws a sub-complex of `subcomplex-size` points, and replaces the worst of them with the reflection of that point through the centroid of the others, when it is better. Otherwise with the contraction, when it is better. Otherwise with a random point of the box of the complex. Reflections are kept within the bounds."
  [complex {:keys [^long subcomplex-size evolution-steps evaluate lo hi]} rng]
  (let [size (count complex)
        xs (map :x complex)
        box-lo (reduce v/emn xs)
        box-hi (reduce v/emx xs)
        step (fn [complex]
               (let [parents (select-parents size subcomplex-size rng)
                     worst (complex (peek parents))
                     value-of-worst (value-of worst)
                     centroid (v/average-vectors (map #(:x (complex %)) (pop parents)))
                     better-than-worst (fn [x]
                                         (let [point (evaluate x)]
                                           (when (m/< (value-of point) value-of-worst)
                                             point)))
                     replacement (or (better-than-worst (clip (v/interpolate (:x worst) centroid 2.0) lo hi))
                                     (better-than-worst (v/interpolate (:x worst) centroid 0.5))
                                     (evaluate (random-point box-lo box-hi rng)))]
                 (sorted-by-value (assoc complex (peek parents) replacement))))]
    (nth (iterate step complex) evolution-steps)))

(defn- complex-of
  "The `k`-th complex of `ngs`: every `ngs`-th point of the sorted population starting at `k`."
  [population ^long k ^long ngs]
  (vec (take-nth ngs (drop k population))))

(defn- unwrap-execution
  "Runs `thunk` and rethrows the exception of a worker thread in place of its `ExecutionException` wrapper."
  [thunk]
  (try
    (thunk)
    (catch ExecutionException e
      (throw (or (ex-cause e) e)))))

(defn- evolve-population
  "One shuffling loop: divides the sorted population into complexes, evolves every complex and merges them into a new sorted population.

  The complexes are evolved in parallel when `:parallel?`, each with its own generator created by `r/child-rngs`. Otherwise `rng` is used directly."
  [population {:keys [^long complex-size parallel?] :as settings} rng]
  (let [ngs (m/long-quot (count population) complex-size)
        complexes (mapv #(complex-of population % ngs) (range ngs))
        evolve (fn [complex complex-rng] (evolve-complex complex settings complex-rng))
        evolved (if parallel?
                  (let [rngs (r/child-rngs rng ngs)]
                    (unwrap-execution #(doall (pmap evolve complexes rngs))))
                  (map #(evolve % rng) complexes))]
    (sorted-by-value (mapcat identity evolved))))

;; termination

(defn- population-range
  "Largest range of the population over the dimensions, as a fraction of the range of the bounds."
  ^double [population {:keys [lo hi]}]
  (let [xs (map :x population)
        spread (v/sub (reduce v/emx xs) (reduce v/emn xs))]
    (double (reduce m/max 0.0 (map m// spread (v/sub hi lo))))))

(defn- improvement-converged?
  "True when the best value, with `stop-loops` loops back, differs from the current one by not more than the tolerance relative to the current one. Non-finite values never converge."
  [bests ^long stop-loops ^double tolerance]
  (let [loops (m/dec (count bests))]
    (and (m/>= loops stop-loops)
         (let [previous (double (bests (m/- loops stop-loops)))
               current (double (peek bests))]
           (and (m/valid-double? previous)
                (m/valid-double? current)
                (m/<= (m/abs (m/- previous current)) (m/* tolerance (m/+ (m/abs current) tolerance))))))))

(defn- stalled?
  "True when the last loop did not improve the best value."
  [bests]
  (m/>= (double (peek bests)) (double (peek (pop bests)))))

;; recovery of a collapsed population

(defn- lost-directions
  "Returns the unit directions, in the population normalized to the unit cube, with almost no variance.

  A direction is lost when its eigenvalue of the covariance matrix is below `pca-ratio` of the largest one. Returns an empty vector when the population has no variance."
  [population {:keys [lo hi]}]
  (let [extent (v/sub hi lo)
        rows (map #(mapv m// (v/sub (:x %) lo) extent) population)
        eigen (mat/eigen-decomposition (mat/rows->mat (stats/covariance-matrix (apply mapv vector rows)))
                                       {:eigenvectors-scaling :raw})
        values (vec (mat/decomposition-component eigen :real-eigenvalues))
        largest (stats/maximum values)]
    (if (m/pos? largest)
      (let [limit (m/* pca-ratio largest)]
        (into []
              (comp (keep-indexed (fn [i direction] (when (m/< (double (values i)) limit) direction)))
                    (map v/vec->array))
              (mat/decomposition-component eigen :eigenvectors)))
      [])))

(defn- recover-dimensions
  "Moves some points of a collapsed population along the lost directions, evaluates them and returns the sorted population.

  For every lost direction, a fraction of the points except the best one is moved by a random step in both senses along the direction, in the population normalized to the unit cube, and then kept within the bounds. The best point is never moved."
  [population {:keys [lo hi evaluate] :as settings} rng]
  (let [directions (lost-directions population settings)]
    (if (empty? directions)
      population
      (let [size (count population)
            moved (m/max 1 (long (m/floor (m/* pca-fraction size))))
            extent (v/sub hi lo)
            sampler (RandomDataGenerator. ^RandomGenerator rng)
            moves (vec (for [direction directions
                             i (.nextPermutation sampler (int (m/dec size)) (int moved))]
                         [(v/emult direction extent) (m/inc (long i))]))
            move (fn [population [direction i]]
                   (let [step (m/* pca-step (m/dec (m/* 2.0 (r/drandom rng))))]
                     (assoc population i (evaluate (clip (v/add (:x (population i)) (v/mult direction step)) lo hi)))))]
        (sorted-by-value (reduce move population moves))))))

;; optimizer

(defn- reduce-complexes
  "Removes the worst complex of a sorted population when there are more than `:min-complexes` complexes."
  [population {:keys [^long complex-size ^long min-complexes]}]
  (if (m/> (m/quot (count population) complex-size) min-complexes)
    (subvec population 0 (m/long-sub (count population) complex-size))
    population))

(defn- result
  "Builds the result of the run from the final sorted population."
  [population {:keys [sign complex-size stats? evaluate]} status iteration]
  (let [best (first population)
        point (vec (:x best))
        value (m/* (double sign) (value-of best))]
    (if stats?
      {:point point
       :value value
       :evaluations (evaluate)
       :iterations iteration
       :status status
       :complexes (quot (count population) (long complex-size))}
      [point value])))

(defn sceua
  "Finds a minimum or a maximum of a function in a box with the Shuffled Complex Evolution method.

  The method is stochastic: use `:rng` with a seeded generator for reproducible results. It needs many evaluations of the objective, so it fits cheap functions with many local extrema, or a problem where a derivative is not available.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence of numbers (a `double` array which it must not change), or as separate arguments when `:vector-arg?` is `false`, and returns a number. `NaN` is treated as the worst value.
  - `opts` (map):
    - `:bounds` (required) - sequence of `[lo hi]` pairs, one for each dimension (`[lo hi]` for one dimension). The bounds are finite and `lo < hi`, up to 1000 dimensions.
    - `:initial` - a point within the bounds which joins the initial population, default: none.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - how `f` receives the point, default: `true`.
    - `:complexes` - number of complexes, default: `5`.
    - `:complex-size` - number of points in a complex, default: `2n+1` for `n` dimensions. The population has `:complexes` times `:complex-size` points.
    - `:subcomplex-size` - number of points drawn from a complex in an evolution step, from `2` to `:complex-size`, default: `n+1`, but not more than `:complex-size`.
    - `:evolution-steps` - number of evolution steps of a complex in one loop, default: `:complex-size`.
    - `:min-complexes` - the number of complexes drops by one after every loop, with the worst points removed, down to this number. Default: `:complexes`, so the number of complexes stays constant.
    - `:stop-loops`, `:stop-improvement` - the run ends when the best value, `:stop-loops` loops back, differs from the current one by not more than `:stop-improvement` relative to the current value (absolute near zero), defaults: `40` and `1.0e-5`. A shorter patience stops runs on flat parts of the objective before the minimum is reached.
    - `:stop-range` - the run ends when the range of the population is not greater than this fraction of the range of the bounds in every dimension, default: `1.0e-3`. Every tenfold decrease makes the result about a hundred times more accurate and needs about 20% more evaluations.
    - `:max-evals` - the maximum number of evaluations of `f`, default: `10000`. Exceeding it throws `ex-info`.
    - `:max-iters` - the maximum number of loops, default: `10000`. The run ends and returns the best point found.
    - `:rng` - random number generator (see [[fastmath.random/rng]]), default: a new `JDKRandomGenerator`. It draws the jitter of the initial population and the random choices of the evolution.
    - `:jitter` - jitter of the sequence which samples the initial population, from `0.0` to `1.0`, default: `0.25`. See [[fastmath.random/jittered-sequence-generator]].
    - `:parallel?` - evolve the complexes in parallel, default: `false`. The function `f` has to be thread safe. Every complex gets its own generator created from `:rng` (see [[fastmath.random/child-rngs]]), so the result can differ from the sequential one, a generator made by [[fastmath.random/synced-rng]] gives children of the default type.
    - `:pca-recovery?` - repair a population which collapsed into a subspace, default: `false`. When the best value did not improve in a loop, the covariance matrix of the population in the unit cube is decomposed. Some points are moved along the directions with almost no variance (eigenvalue below `1.0e-3` of the largest one) and evaluated. With the default sizes the population rarely collapses and the option changes little. It helped, with the cost of 4% to 25% more evaluations, when the population is small for the problem (for example `:complexes` 2 or 3 for ten and more dimensions).
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`, where `point` is a vector and `value` the value of `f` there, also when maximizing. With `:stats?` returns a map with:

  - `:point`, `:value` - the same,
  - `:evaluations` - number of evaluations of `f`,
  - `:iterations` - number of loops,
  - `:status` - why the run ended: `:converged-range`, `:converged-improvement` or `:max-iterations`,
  - `:complexes` - number of complexes in the last loop.

  The population is evaluated and the first loop is always done, so at least `:complexes` times `:complex-size` evaluations are made, and `:max-evals` below that number throws.

  Throws `ex-info` for invalid bounds, options or `:initial`, and when `:max-evals` is exceeded.

  See also [[fastmath.optimization/minimize]], [[fastmath.optimization/scan-and-minimize]], [[fastmath.optimization.acm/cmaes]]."
  [f opts]
  (let [options (parse-options opts)
        {:keys [sign max-evals stop-loops stop-improvement ^double stop-range ^long max-iters rng pca-recovery? lo]} options
        n (count lo)
        evaluate (evaluator (common/->vector-fn f (common/resolve-vector-arg? :sceua (:vector-arg? opts))) sign max-evals)
        settings (assoc options :evaluate evaluate)]
    (loop [population (initial-population settings)
           iteration 0
           bests [(value-of (first population))]]
      (let [evolved (evolve-population population settings rng)
            iteration (m/inc iteration)
            bests (conj bests (value-of (first evolved)))
            status (cond
                     (m/<= (population-range evolved settings) stop-range) :converged-range
                     (improvement-converged? bests stop-loops stop-improvement) :converged-improvement
                     (m/>= iteration max-iters) :max-iterations)]
        (if status
          (result evolved settings status iteration)
          (recur (-> (if (and pca-recovery? (m/>= n 2) (m/> (count evolved) n) (stalled? bests))
                       (recover-dimensions evolved settings rng)
                       evolved)
                     (reduce-complexes settings))
                 iteration
                 bests))))))
