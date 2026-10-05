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
            [fastmath.optimization.common :as common])
  (:import [java.util Arrays Comparator]
           [java.util.concurrent ExecutionException]
           [java.util.concurrent.atomic AtomicLong]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; A point of the population is a double[] of length n+1: the n coordinates followed by the value of the
;; (minimized) objective. A population and a complex are arrays of such points sorted by the value.

(def ^:private ^:const max-dimensions 1000)
(def ^:private ^:const pca-ratio 1.0e-3)
(def ^:private ^:const pca-fraction 0.1)
(def ^:private ^:const pca-step 0.05)

(defn- option-error
  "Throws `ex-info` about an invalid option."
  [option value reason]
  (throw (ex-info (str "Invalid option " option ": " reason) {:option option :value value :reason reason})))

(defn- integer-option
  "Reads a positive integer option, `nil` gives `default`."
  ^long [options option default]
  (let [v (get options option)]
    (cond
      (nil? v) (long default)
      (and (int? v) (m/pos? (long v))) (long v)
      :else (option-error option v "must be an integer not less than 1"))))

(defn- nonnegative-option
  "Reads a finite option not less than zero, `nil` gives `default`."
  ^double [options option default]
  (let [v (get options option)]
    (cond
      (nil? v) (double default)
      (and (number? v) (m/valid-double? (double v)) (m/not-neg? (double v))) (double v)
      :else (option-error option v "must be a finite number not less than 0"))))

(defn- parse-initial
  "Validates the initial point and returns it as an array or `nil`."
  [initial ^doubles lo ^doubles hi]
  (when (some? initial)
    (let [xs (if (number? initial) [initial] (seq initial))]
      (when-not (every? number? xs) (option-error :initial initial "must be a number or a sequence of numbers"))
      (let [x (double-array xs)]
        (dotimes [d (alength x)]
          (let [v (aget x d)]
            (when-not (and (m/valid-double? v) (m/<= (aget lo d) v (aget hi d)))
              (option-error :initial initial "must be finite and lie within the bounds"))))
        x))))

(defn- parse-options
  "Validates the options of [[sceua]] and returns them as a map with the defaults applied."
  [f opts]
  (let [opts (or opts {})
        {:keys [bounds initial goal vector-arg? rng parallel? pca-recovery? stats?]} opts
        bounds (common/normalize-bounds :sceua bounds initial)
        n (count bounds)
        _ (when (m/> n max-dimensions)
            (option-error :bounds bounds (str "at most " max-dimensions " dimensions are supported")))
        lo (double-array (map first bounds))
        hi (double-array (map second bounds))
        sign (if (= :maximize (common/parse-goal goal)) -1.0 1.0)
        complexes (integer-option opts :complexes 5)
        complex-size (integer-option opts :complex-size (m/inc (m/* 2 n)))
        subcomplex-size (integer-option opts :subcomplex-size (m/min (m/inc n) complex-size))
        evolution-steps (integer-option opts :evolution-steps complex-size)
        min-complexes (integer-option opts :min-complexes complexes)
        jitter (nonnegative-option opts :jitter 0.25)]
    (when (m/< complex-size 2)
      (option-error :complex-size complex-size "must be at least 2"))
    (when-not (m/<= 2 subcomplex-size complex-size)
      (option-error :subcomplex-size subcomplex-size "must be at least 2 and not greater than :complex-size"))
    (when (m/> min-complexes complexes)
      (option-error :min-complexes min-complexes "must not be greater than :complexes"))
    (when (m/> jitter 1.0)
      (option-error :jitter jitter "must not be greater than 1"))
    {:n n
     :lo lo
     :hi hi
     :initial (parse-initial initial lo hi)
     :sign sign
     :f (common/->vector-fn f (common/resolve-vector-arg? :sceua vector-arg?))
     :complexes complexes
     :complex-size complex-size
     :subcomplex-size subcomplex-size
     :evolution-steps evolution-steps
     :min-complexes min-complexes
     :stop-loops (integer-option opts :stop-loops 7)
     :stop-improvement (nonnegative-option opts :stop-improvement 1.0e-5)
     :stop-range (nonnegative-option opts :stop-range 1.0e-3)
     :max-evals (integer-option opts :max-evals 10000)
     :max-iters (integer-option opts :max-iters 10000)
     :jitter jitter
     :rng (r/ensure-rng rng)
     :parallel? (boolean parallel?)
     :pca-recovery? (boolean pca-recovery?)
     :stats? (boolean stats?)}))

;; points

(def ^:private ^Comparator by-value
  (reify Comparator
    (compare [_ a b]
      (let [^doubles a a
            ^doubles b b]
        (Double/compare (aget a (m/dec (alength a))) (aget b (m/dec (alength b))))))))

(defn- value-of
  "Value of the objective stored in a point."
  ^double [^doubles point]
  (aget point (m/dec (alength point))))

(defn- coordinates
  "Coordinates of a point as a vector."
  [^doubles point]
  (vec (Arrays/copyOf point (int (m/dec (alength point))))))

(defn- evaluator
  "Creates the function of a coordinates array which returns the value of the minimized objective.

  Every call is counted. The first call beyond `max-evals` throws `ex-info`. `NaN` is replaced by `+Infinity`, the worst value. The objective receives a new vector."
  [f ^double sign ^AtomicLong counter ^long max-evals]
  (fn ^double [^doubles x]
    (when (m/> (.incrementAndGet counter) max-evals)
      (throw (ex-info "Maximum number of evaluations exceeded" {:max-evals max-evals :evaluations max-evals})))
    (let [v (m/* sign (double (f (vec x))))]
      (if (m/nan? v) ##Inf v))))

(defn- make-point
  "Evaluates the objective at the coordinates and returns a new point."
  ^doubles [evaluate ^doubles x]
  (let [n (alength x)
        point (Arrays/copyOf x (int (m/inc n)))]
    (aset point n (double (evaluate x)))
    point))

(defn- unit-point->array
  ^doubles [p ^long n]
  (if (m/one? n) (double-array [p]) (double-array (seq p))))

(defn- initial-population
  "Creates the sorted population: the initial point, if given, and points of a jittered low discrepancy sequence scaled to the bounds."
  ^objects [evaluate ^doubles lo ^doubles hi initial size jitter rng]
  (let [n (alength lo)
        size (long size)
        jitter (double jitter)
        population (object-array size)
        start (if (some? initial)
                (do (aset population 0 (make-point evaluate initial)) 1)
                0)]
    (loop [i start
           points (r/jittered-sequence-generator (if (m/< n 15) :r2 :sobol) n jitter rng)]
      (when (m/< i size)
        (let [^doubles u (unit-point->array (first points) n)
              x (double-array n)]
          (dotimes [d n]
            (let [l (aget lo d)]
              (aset x d (m/constrain (m/+ l (m/* (aget u d) (m/- (aget hi d) l))) l (aget hi d)))))
          (aset population i (make-point evaluate x))
          (recur (m/inc i) (rest points)))))
    (Arrays/sort population by-value)
    population))

;; competitive complex evolution

(defn- triangular-index
  "Zero-based index drawn with the probability proportional to `m - index` from a uniform `u` in [0,1)."
  ^long [^long m ^double u]
  (let [mh (m/+ m 0.5)
        root (m/sqrt (m/max 0.0 (m/- (m/* mh mh) (m/* m (m/inc m) u))))
        index (long (m/floor (m/- mh root)))]
    (m/min (m/max index 0) (m/dec m))))

(defn- select-parents
  "Returns `q` different indices in ascending order, drawn from `m` ones with the triangular probability (the lower the index, the more probable)."
  ^longs [^long m ^long q rng]
  (let [chosen (boolean-array m)]
    (loop [k 0]
      (when (m/< k q)
        (let [i (triangular-index m (r/drandom rng))]
          (if (aget chosen i)
            (recur k)
            (do (aset chosen i true)
                (recur (m/inc k)))))))
    (let [result (long-array q)]
      (loop [i 0
             j 0]
        (if (m/< j q)
          (if (aget chosen i)
            (do (aset result j i)
                (recur (m/inc i) (m/inc j)))
            (recur (m/inc i) j))
          result)))))

(defn- bounding-box
  "Returns `[lower upper]` arrays of the smallest box which contains the points of the complex."
  [^objects complex ^long n]
  (let [lower (double-array n)
        upper (double-array n)
        size (alength complex)]
    (dotimes [d n]
      (loop [i 0
             mn ##Inf
             mx ##-Inf]
        (if (m/< i size)
          (let [v (aget ^doubles (aget complex i) d)]
            (recur (m/inc i) (m/min mn v) (m/max mx v)))
          (do (aset lower d mn)
              (aset upper d mx)))))
    [lower upper]))

(defn- evolve-complex
  "Evolves a complex (an array of points sorted by value) and returns a new sorted array of the same size.

  Every step draws a sub-complex of `subcomplex-size` points, and replaces the worst of them with the reflection of that point through the centroid of the others, when it is better. Otherwise with the contraction, when it is better. Otherwise with a random point of the box of the complex. Points are kept within the bounds. The input array is not changed."
  ^objects [^objects complex {:keys [subcomplex-size evolution-steps evaluate lo hi]} rng]
  (let [^doubles lo lo
        ^doubles hi hi
        q (long subcomplex-size)
        n (alength lo)
        ^objects evolved (aclone complex)
        size (alength evolved)
        [^doubles box-lo ^doubles box-hi] (bounding-box complex n)]
    (dotimes [_ (long evolution-steps)]
      (let [^longs parents (select-parents size q rng)
            worst-index (aget parents (m/dec q))
            ^doubles worst (aget evolved worst-index)
            worst-value (value-of worst)
            centroid (double-array n)]
        (dotimes [j (m/dec q)]
          (let [^doubles p (aget evolved (aget parents j))]
            (dotimes [d n]
              (aset centroid d (m/+ (aget centroid d) (aget p d))))))
        (dotimes [d n]
          (aset centroid d (m// (aget centroid d) (m/dec q))))
        (let [x (double-array n)]
          (dotimes [d n]
            (aset x d (m/constrain (m/- (m/* 2.0 (aget centroid d)) (aget worst d)) (aget lo d) (aget hi d))))
          (let [reflection (make-point evaluate x)]
            (if (m/< (value-of reflection) worst-value)
              (aset evolved worst-index reflection)
              (let [y (double-array n)]
                (dotimes [d n]
                  (aset y d (m/* 0.5 (m/+ (aget centroid d) (aget worst d)))))
                (let [contraction (make-point evaluate y)]
                  (if (m/< (value-of contraction) worst-value)
                    (aset evolved worst-index contraction)
                    (let [z (double-array n)]
                      (dotimes [d n]
                        (aset z d (m/+ (aget box-lo d) (m/* (r/drandom rng) (m/- (aget box-hi d) (aget box-lo d))))))
                      (aset evolved worst-index (make-point evaluate z)))))))))
        (Arrays/sort evolved by-value)))
    evolved))

(defn- complex-of
  "The `k`-th complex of `ngs`: every `ngs`-th point of the sorted population starting at `k`."
  ^objects [^objects population ^long k ^long ngs ^long size]
  (let [complex (object-array size)]
    (dotimes [i size]
      (aset complex i (aget population (m/+ k (m/* i ngs)))))
    complex))

(defn- unwrap-execution
  "Runs `thunk` and rethrows the exception of a worker thread in place of its `ExecutionException` wrapper."
  [thunk]
  (try
    (thunk)
    (catch ExecutionException e
      (throw (or (ex-cause e) e)))))

(defn- evolve-population
  "One shuffling loop: divides the sorted population into complexes, evolves every complex and merges them into a new sorted population.

  The complexes are evolved in parallel when `parallel?`, each with its own generator created by `r/child-rngs`. Otherwise `rng` is used directly."
  ^objects [^objects population ngs settings rng parallel?]
  (let [ngs (long ngs)
        size (long (:complex-size settings))
        complexes (mapv #(complex-of population % ngs size) (range ngs))
        evolve (fn [complex complex-rng] (evolve-complex complex settings complex-rng))
        evolved (if parallel?
                  (let [rngs (r/child-rngs rng ngs)]
                    (unwrap-execution #(vec (doall (pmap evolve complexes rngs)))))
                  (mapv #(evolve % rng) complexes))
        merged (object-array (m/* ngs size))]
    (dotimes [k ngs]
      (System/arraycopy ^objects (evolved k) 0 merged (int (m/* k size)) (int size)))
    (Arrays/sort merged by-value)
    merged))

;; termination

(defn- population-range
  "Largest range of the population over the dimensions, as a fraction of the range of the bounds."
  ^double [^objects population ^doubles lo ^doubles hi]
  (let [size (alength population)
        n (alength lo)]
    (loop [d 0
           largest 0.0]
      (if (m/< d n)
        (let [spread (double (loop [i 0
                                    mn ##Inf
                                    mx ##-Inf]
                               (if (m/< i size)
                                 (let [v (aget ^doubles (aget population i) d)]
                                   (recur (m/inc i) (m/min mn v) (m/max mx v)))
                                 (m/- mx mn))))]
          (recur (m/inc d) (m/max largest (m// spread (m/- (aget hi d) (aget lo d))))))
        largest))))

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
  [^objects population ^doubles lo ^doubles hi]
  (let [size (alength population)
        n (alength lo)
        columns (mapv (fn [d]
                        (let [d (long d)
                              l (aget lo d)
                              extent (m/- (aget hi d) l)]
                          (mapv (fn [i] (m// (m/- (aget ^doubles (aget population (long i)) d) l) extent))
                                (range size))))
                      (range n))
        eigen (mat/eigen-decomposition (mat/rows->mat (stats/covariance-matrix columns))
                                       {:eigenvectors-scaling :raw})
        values (double-array (seq (mat/decomposition-component eigen :real-eigenvalues)))
        vectors (mat/decomposition-component eigen :eigenvectors)
        largest (double (loop [i 0
                               mx ##-Inf]
                          (if (m/< i (alength values))
                            (recur (m/inc i) (m/max mx (aget values i)))
                            mx)))]
    (if (m/pos? largest)
      (let [limit (m/* pca-ratio largest)]
        (into []
              (comp (keep-indexed (fn [i v] (when (m/< (aget values (long i)) limit) v)))
                    (map #(double-array (seq %))))
              vectors))
      [])))

(defn- sample-indices
  "Returns `k` different indices from `from` (inclusive) to `to` (exclusive)."
  ^longs [rng ^long from ^long to ^long k]
  (let [pool (long-array (range from to))
        size (alength pool)]
    (dotimes [j k]
      (let [i (m/+ j (r/irandom rng (m/- size j)))
            tmp (aget pool j)]
        (aset pool j (aget pool i))
        (aset pool i tmp)))
    (Arrays/copyOf pool (int k))))

(defn- recover-dimensions!
  "Moves some points of a collapsed population along the lost directions, evaluates them and sorts the population again.

  For every lost direction, a fraction of the points except the best one is moved by a random step in both senses along the direction, in the population normalized to the unit cube, and then kept within the bounds. The best point is never moved. Changes `population` and returns it."
  ^objects [^objects population rng evaluate ^doubles lo ^doubles hi]
  (let [directions (lost-directions population lo hi)
        size (alength population)
        n (alength lo)
        moved (long (m/max 1 (long (m/floor (m/* pca-fraction size)))))]
    (when (seq directions)
      (doseq [^doubles direction directions
              i (sample-indices rng 1 size moved)]
        (let [^doubles point (aget population i)
              step (m/* pca-step (m/dec (m/* 2.0 (r/drandom rng))))
              x (double-array n)]
          (dotimes [d n]
            (let [l (aget lo d)
                  extent (m/- (aget hi d) l)
                  normalized (m/+ (m// (m/- (aget point d) l) extent) (m/* step (aget direction d)))]
              (aset x d (m/constrain (m/+ l (m/* normalized extent)) l (aget hi d)))))
          (aset population (int i) (make-point evaluate x))))
      (Arrays/sort population by-value))
    population))

;; optimizer

(defn- result
  "Builds the result of the run from the final sorted population."
  [^objects population settings status ngs bests ^AtomicLong counter]
  (let [{:keys [^double sign stats?]} settings
        report (fn [^doubles point] [(coordinates point) (m/* sign (value-of point))])
        [point value] (report (aget population 0))]
    (if stats?
      {:point point
       :value value
       :evaluations (.get counter)
       :iterations (m/dec (count bests))
       :status status
       :complexes ngs
       :history (mapv #(m/* sign (double %)) (rest bests))
       :population (mapv report population)}
      [point value])))

(defn sceua
  "Finds a minimum or a maximum of a function in a box with the Shuffled Complex Evolution method.

  The method is stochastic: use `:rng` with a seeded generator for reproducible results. It needs many evaluations of the objective, so it fits cheap functions with many local extrema, or a problem where a derivative is not available.

  Parameters:

  - `f` (function): the objective. It receives the point as one vector, or as separate arguments when `:vector-arg?` is `false`, and returns a number. `NaN` is treated as the worst value.
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
    - `:stop-loops`, `:stop-improvement` - the run ends when the best value, `:stop-loops` loops back, differs from the current one by not more than `:stop-improvement` relative to the current value (absolute near zero), defaults: `7` and `1.0e-5`.
    - `:stop-range` - the run ends when the range of the population is not greater than this fraction of the range of the bounds in every dimension, default: `1.0e-3`.
    - `:max-evals` - the maximum number of evaluations of `f`, default: `10000`. Exceeding it throws `ex-info`.
    - `:max-iters` - the maximum number of loops, default: `10000`. The run ends and returns the best point found.
    - `:rng` - random number generator (see [[fastmath.random/rng]]), default: a new `JDKRandomGenerator`. It draws the jitter of the initial population and the random choices of the evolution.
    - `:jitter` - jitter of the sequence which samples the initial population, from `0.0` to `1.0`, default: `0.25`. See [[fastmath.random/jittered-sequence-generator]].
    - `:parallel?` - evolve the complexes in parallel, default: `false`. The function `f` has to be thread safe. Every complex gets its own generator created from `:rng` (see [[fastmath.random/child-rngs]]), so the result can differ from the sequential one, a generator made by [[fastmath.random/synced-rng]] gives children of the default type.
    - `:pca-recovery?` - repair a population which collapsed into a subspace, default: `false`. When the best value did not improve in a loop, the covariance matrix of the population in the unit cube is decomposed. Some points are moved along the directions with almost no variance (eigenvalue below `1.0e-3` of the largest one) and evaluated.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`, where `point` is a vector and `value` the value of `f` there, also when maximizing. With `:stats?` returns a map with:

  - `:point`, `:value` - the same,
  - `:evaluations` - number of evaluations of `f`,
  - `:iterations` - number of loops,
  - `:status` - why the run ended: `:converged-range`, `:converged-improvement` or `:max-iterations`,
  - `:complexes` - number of complexes in the last loop,
  - `:history` - the best value after every loop, a vector of `:iterations` values,
  - `:population` - the final population, `[point value]` pairs sorted from the best one.

  The population is evaluated and the first loop is always done, so at least `:complexes` times `:complex-size` evaluations are made, and `:max-evals` below that number throws.

  Throws `ex-info` for invalid bounds, options or `:initial`, and when `:max-evals` is exceeded.

  See also [[fastmath.optimization/minimize]], [[fastmath.optimization/scan-and-minimize]], [[fastmath.optimization.acm/cmaes]]."
  [f opts]
  (let [settings (parse-options f opts)
        {:keys [^long n ^doubles lo ^doubles hi initial ^double sign ^long complexes ^long complex-size
                ^long min-complexes ^long stop-loops ^double stop-improvement ^double stop-range
                ^long max-evals ^long max-iters ^double jitter rng parallel? pca-recovery?]} settings
        counter (AtomicLong. 0)
        evaluate (evaluator (:f settings) sign counter max-evals)
        settings (assoc settings :evaluate evaluate)
        population (initial-population evaluate lo hi initial (m/* complexes complex-size) jitter rng)]
    (loop [population population
           ngs complexes
           bests [(value-of (aget population 0))]]
      (let [evolved (evolve-population population ngs settings rng parallel?)
            bests (conj bests (value-of (aget evolved 0)))
            loops (m/dec (count bests))
            status (cond
                     (m/<= (population-range evolved lo hi) stop-range) :converged-range
                     (improvement-converged? bests stop-loops stop-improvement) :converged-improvement
                     (m/>= loops max-iters) :max-iterations)]
        (if status
          (result evolved settings status ngs bests counter)
          (let [recovered (if (and pca-recovery?
                                   (m/>= n 2)
                                   (m/> (alength evolved) n)
                                   (stalled? bests))
                            (recover-dimensions! evolved rng evaluate lo hi)
                            evolved)]
            (if (m/> ngs min-complexes)
              (recur (Arrays/copyOf recovered (int (m/- (alength recovered) complex-size))) (m/dec ngs) bests)
              (recur recovered ngs bests))))))))
