(ns fastmath.ml.clustering
  "Clustering of points: k-means++, fuzzy k-means and DBSCAN, with cluster validation helpers.

  Algorithms are backed by Apache Commons Math. Points are sequences of numbers (or arrays); a plain number is a 1-dimensional point. All points have to have the same, non-zero dimension and finite coordinates.

  Every clustering function returns a vector of cluster maps. A cluster map contains:

  - `:representative` - the centre of the cluster (centroid for k-means++ and fuzzy k-means, the mean of cluster points for DBSCAN),
  - `:size` - the number of points,
  - `:indices` - positions of the cluster points in the input sequence (the order within a cluster is not specified),
  - `:data` - the cluster points (vectors of doubles), unless `:add-data?` is `false`.

  The order of clusters and the choice of representatives depend on the random number generator. The default `:rng` is a new generator for every call, pass a seeded generator (`fastmath.random/rng`) for reproducible results.

  Options are validated: an invalid option throws `ex-info` with the keys `:option`, `:value` and `:expected`. `:value` is the invalid value; it is `nil` when the value is the whole input or a whole result (data can be large), and a single point when the points differ in dimension. Unknown option keys are ignored.

  `:distance` is a keyword (`:euclidean`, `:chebyshev`, `:manhattan`, `:canberra`, `:earth-movers`), a function of two points (ie. a function from `fastmath.distance`) or an Apache Commons Math `DistanceMeasure`. A distance has to be non-negative and equal to 0 for equal points (`d(x, x) = 0`): the radius selection of [[dbscan]] and [[silhouette]] rely on it and return wrong values otherwise.

  Functions: [[kmeans++]], [[fuzzy-kmeans]], [[dbscan]], [[infer-dbscan-radius]] (clustering and the radius selection), [[predict]], [[wcss]], [[silhouette]] (working on a clustering result)."
  (:require [fastmath.distance :as dists]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [org.apache.commons.math3.ml.clustering
            KMeansPlusPlusClusterer
            KMeansPlusPlusClusterer$EmptyClusterStrategy
            MultiKMeansPlusPlusClusterer
            FuzzyKMeansClusterer
            DBSCANClusterer
            Clusterable Cluster Clusterer CentroidCluster]
           [org.apache.commons.math3.ml.distance
            DistanceMeasure
            EuclideanDistance ManhattanDistance ChebyshevDistance CanberraDistance EarthMoversDistance]
           [org.apache.commons.math3.linear RealMatrix]
           [org.apache.commons.math3.random RandomGenerator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Validation

(defn- invalid-option
  [option value expected]
  (throw (ex-info (str "Invalid " (name option) ", expected " expected ".")
                  {:option option :value value :expected expected})))

(defn- check-integer
  "Return `value` as a long when it is an integer not smaller than `min-value`."
  ^long [option value ^long min-value expected]
  (if (and (int? value) (m/>= (long value) min-value))
    (long value)
    (invalid-option option value expected)))

(defn- check-double
  "Return `value` as a double when it is a finite number greater than `lower` (or not smaller when `inclusive?`)."
  [option value lower inclusive? expected]
  (if (and (number? value)
           (m/valid-double? (double value))
           (if inclusive? (m/>= (double value) (double lower)) (m/> (double value) (double lower))))
    (double value)
    (invalid-option option value expected)))

(defn- check-rng
  [rng]
  (if (instance? RandomGenerator rng) rng (invalid-option :rng rng "a random number generator")))

(defn- check-member
  [option value allowed]
  (if (contains? allowed value)
    value
    (invalid-option option value (str "one of " (sort allowed)))))

;; ## Distances

(def ^:private acm-distances
  {:euclidean #(EuclideanDistance.)
   :chebyshev #(ChebyshevDistance.)
   :manhattan #(ManhattanDistance.)
   :canberra #(CanberraDistance.)
   :earth-movers #(EarthMoversDistance.)})

(defn- ->distance
  "Resolve a distance keyword, function or `DistanceMeasure` to a `DistanceMeasure`."
  ^DistanceMeasure [distance]
  (cond
    (keyword? distance) (if-let [constructor (acm-distances distance)]
                          (dists/->acm-distance (constructor))
                          (invalid-option :distance distance
                                          (str "one of " (sort (keys acm-distances))
                                               ", a function of two points or a DistanceMeasure")))
    (or (fn? distance) (instance? DistanceMeasure distance)) (dists/->acm-distance distance)
    :else (invalid-option :distance distance
                          (str "one of " (sort (keys acm-distances))
                               ", a function of two points or a DistanceMeasure"))))

;; ## Points

(deftype IndexedPoint [^doubles point ^long index]
  Clusterable
  (getPoint [_] point))

(defn- point-array
  "Convert a number, a sequence or an array to a double array of finite numbers."
  ^doubles [option x]
  (let [expected "a number, a sequence or an array of finite numbers"
        arr (try (v/vec->array x)
                 (catch Exception _ (invalid-option option x expected)))
        ^doubles arr (if (instance? (Class/forName "[D") arr) arr (invalid-option option x expected))]
    (when (or (m/zero? (alength arr))
              (not (every? m/valid-double? arr)))
      (invalid-option option x expected))
    arr))

(defn- ->points
  "Validate the data and return a vector of `IndexedPoint`."
  [xss]
  (when-not (and (seqable? xss) (not (string? xss)) (not (map? xss)) (seq xss))
    (invalid-option :xss (when-not (coll? xss) xss) "a non-empty sequence of points"))
  (let [points (into [] (map-indexed (fn [^long i x] (IndexedPoint. (point-array :xss x) i))) xss)
        dim (alength ^doubles (.-point ^IndexedPoint (first points)))]
    (when-let [^IndexedPoint different (first (remove #(m/== dim (alength ^doubles (.-point ^IndexedPoint %))) points))]
      (invalid-option :xss (vec (.-point different)) (str "points of the same dimension, the first point has " dim)))
    points))

(defn- point-dim
  ^long [points]
  (alength ^doubles (.-point ^IndexedPoint (first points))))

(defn- points->arrays
  "Vector of double arrays of `points`."
  [points]
  (mapv #(.-point ^IndexedPoint %) points))

;; ## Results

(defn- point->vec [^Clusterable p] (vec (.getPoint p)))

(defn- cluster->map
  [add-data? ^Cluster cluster]
  (let [points (.getPoints cluster)
        data (mapv point->vec points)
        res {:representative (if (instance? CentroidCluster cluster)
                               (point->vec (.getCenter ^CentroidCluster cluster))
                               (v/average-vectors data))
             :size (count data)
             :indices (mapv #(.-index ^IndexedPoint %) points)}]
    (if add-data? (assoc res :data data) res)))

(defn- ->clusters
  [points ^Clusterer clusterer add-data?]
  (mapv (partial cluster->map add-data?) (.cluster clusterer points)))

;; ## Radius selection

(defn- nth-smallest
  "The element at the position `k` of the sorted numbers `xs`.

  One pass: the `k+1` smallest numbers are kept in a sorted multiset (number -> count), the largest of them is the answer."
  [xs ^long k]
  (let [size (m/inc k)
        drop-largest (fn [kept]
                       (let [[largest cnt] (first (rseq kept))]
                         (if (m/== (long cnt) 1) (dissoc kept largest) (assoc kept largest (m/long-dec (long cnt))))))
        add (fn [kept x] (update kept x (fnil m/inc 0)))
        [kept _] (reduce (fn [[kept ^long n] ^double x]
                           (cond
                             (m/< n size) [(add kept x) (m/inc n)]
                             (m/>= x (double (key (first (rseq kept))))) [kept n]
                             :else [(add (drop-largest kept) x) n]))
                         [(sorted-map) 0] xs)]
    (key (first (rseq kept)))))

(defn- k-distances
  "Sorted (ascending) distances from every point to its `k`-th nearest neighbour.

  A point is at the distance 0 from itself, which is the first of its distances, so the `k`-th neighbour is the element at the position `k`."
  [arrays ^DistanceMeasure dist ^long k]
  (->> arrays
       (map (fn [^doubles a] (nth-smallest (map (fn [^doubles b] (.compute dist a b)) arrays) k)))
       (sort)
       (vec)))

(defn- difference-curve
  "Kneedle difference curve of an increasing convex curve `ys` (a vector): the curve normalised to the unit square, flipped to be concave, minus the diagonal."
  [ys]
  (let [n- (double (m/long-dec (count ys)))
        ymin (double (first ys))
        yrange (m/- (double (peek ys)) ymin)]
    (vec (map-indexed (fn [^long i ^double y] (m/- (m/- 1.0 (m// (m/- y ymin) yrange)) (m// i n-))) (rseq ys)))))

(defn- difference-rows
  "Rows `[i x next maximum? minimum?]` of the difference curve `d`; the first and the last element are their own neighbours."
  [d]
  (map (fn [i ^double p ^double x ^double nx]
         [i x nx (and (m/>= x p) (m/>= x nx)) (and (m/<= x p) (m/<= x nx))])
       (range) (cons (first d) d) d (concat (rest d) [(peek d)])))

(defn- kneedle
  "Index of the knee of an increasing convex curve `ys` (a sorted vector of doubles) or `nil`.

  Kneedle, Satopaa et al. 2011, offline version as in kneed 0.9.0 (after a local minimum of the difference curve the detection waits for the next local maximum), `x` is the index of `y`."
  [ys ^double sensitivity]
  (let [n (count ys)]
    (when (and (m/>= n 3) (m/pos? (m/- (double (peek ys)) (double (first ys)))))
      (let [n- (m/dec n)
            step (m// 1.0 (double n-))
            rows (->> (difference-rows (difference-curve ys))
                      (drop-while (comp not #(nth % 3)))
                      (take-while #(m/< (long (first %)) n-)))
            ;; the state is the threshold, its index and the detection flag (a local maximum
            ;; switches the detection on, a local minimum off), `reduced` returns the knee index
            knee (reduce (fn [[^double threshold ^long threshold-index active?] [^long i ^double x ^double nx
                                                                                maximum? minimum?]]
                           (let [[^double threshold ^long threshold-index active?] (if maximum?
                                                                                     [(m/- x (m/* sensitivity step)) i true]
                                                                                     [threshold threshold-index active?])
                                 [^double threshold active?] (if minimum? [0.0 false] [threshold active?])]
                             (if (and active? (m/< nx threshold))
                               (reduced (m/- n- threshold-index))
                               [threshold threshold-index active?])))
                         [0.0 0 true] rows)]
        (when (number? knee) knee)))))

(defn- knee-radius
  [points dist k ^double sensitivity]
  (let [curve (k-distances (points->arrays points) dist k)
        idx (kneedle curve sensitivity)
        eps (when idx (curve idx))]
    (if (and eps (m/pos? (double eps)))
      eps
      (invalid-option :eps nil "an explicit :eps, no knee was found in the curve of distances to the neighbours"))))

;; ## Clustering

(def ^:private kmeans-strategies
  {:largest-variance KMeansPlusPlusClusterer$EmptyClusterStrategy/LARGEST_VARIANCE
   :largest-points-number KMeansPlusPlusClusterer$EmptyClusterStrategy/LARGEST_POINTS_NUMBER
   :farthest-point KMeansPlusPlusClusterer$EmptyClusterStrategy/FARTHEST_POINT
   :error KMeansPlusPlusClusterer$EmptyClusterStrategy/ERROR})

(defn- ->int
  "Clamp to the range of int, for the Apache Commons Math constructors."
  ^long [^long x]
  (m/max -1 (m/min x Integer/MAX_VALUE)))

(defn- check-clusters
  ^long [clusters points]
  (let [k (check-integer :clusters clusters 1 "an integer, at least 1")]
    (if (m/> k (count points))
      (invalid-option :clusters clusters "an integer, at most the number of points")
      k)))

(defn- check-distinct-points
  "Throws when there are fewer than `k` distinct points: k-means++ with the unlimited number of iterations never returns then (an empty cluster cannot be filled)."
  [clusters points ^long k]
  (when (m/< (count (take k (distinct (map #(mapv (fn [^double x] (m/+ 0.0 x)) (.-point ^IndexedPoint %)) points)))) k)
    (invalid-option :clusters clusters "an integer, at most the number of distinct points")))

(defn kmeans++
  "K-means++ clustering: partitions points into `:clusters` groups minimising the within-cluster variance.

  Parameters:

  - `xss` - a non-empty sequence of points (sequences, arrays or numbers), see the namespace documentation
  - options (a map):
    - `:clusters` - the number of clusters, an integer between 1 and the number of distinct points, default: `1`
    - `:max-iters` - the maximum number of iterations, a negative value means no limit, default: `-1`. A distance that is 0 for different points (eg. a constant function) can keep the algorithm from converging, set `:max-iters` then
    - `:distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`
    - `:rng` - random number generator, default: a new default generator for every call (seed it for reproducible results)
    - `:empty-cluster-strategy` - what to do with an empty cluster: `:largest-variance` (default), `:largest-points-number`, `:farthest-point` or `:error`
    - `:trials` - the number of runs, the best (the lowest sum of cluster variances) is returned, default: `1`
    - `:add-data?` - add `:data` to the clusters, default: `true`

  Returns a vector of cluster maps (`:representative`, `:size`, `:indices` and `:data`), see the namespace documentation. The order of clusters depends on the generator.

  Throws `ex-info` for an invalid option or data, also when `:clusters` is greater than the number of distinct points. With `:empty-cluster-strategy` set to `:error` an empty cluster throws an exception.

  See also [[fuzzy-kmeans]], [[dbscan]], [[predict]], [[wcss]], [[silhouette]]."
  ([xss] (kmeans++ xss nil))
  ([xss {:keys [clusters max-iters distance rng empty-cluster-strategy trials add-data?]
         :or {clusters 1 max-iters -1 distance :euclidean empty-cluster-strategy :largest-variance
              trials 1 rng (r/rng :default) add-data? true}}]
   (let [points (->points xss)
         k (check-clusters clusters points)
         _ (check-distinct-points clusters points k)
         max-iters (->int (check-integer :max-iters max-iters Long/MIN_VALUE "an integer"))
         trials (->int (check-integer :trials trials 1 "an integer, at least 1"))
         strategy (kmeans-strategies (check-member :empty-cluster-strategy empty-cluster-strategy
                                                   (set (keys kmeans-strategies))))
         dist (->distance distance)
         rng (check-rng rng)
         clusterer (KMeansPlusPlusClusterer. (int k) (int max-iters) dist rng strategy)
         clusterer (if (m/== trials 1)
                     clusterer
                     (MultiKMeansPlusPlusClusterer. clusterer (int trials)))]
     (->clusters points clusterer add-data?))))

(defn fuzzy-kmeans
  "Fuzzy k-means clustering: every point has a degree of membership in every cluster, points are grouped by their highest membership.

  Parameters:

  - `xss` - a non-empty sequence of points (sequences, arrays or numbers), see the namespace documentation
  - options (a map):
    - `:clusters` - the number of clusters, an integer between 1 and the number of points, default: `1`
    - `:fuzziness` - the fuzziness coefficient, a finite number greater than 1, default: `2.0`
    - `:max-iters` - the maximum number of iterations, a negative value means no limit, default: `-1`
    - `:epsilon` - convergence threshold (a change of the membership matrix), a finite number greater than 0, default: `1.0e-3`
    - `:distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`
    - `:rng` - random number generator, default: a new default generator for every call (seed it for reproducible results)
    - `:add-data?` - add `:data` to the clusters, default: `true`
    - `:add-memberships?` - add `:memberships` to the clusters, default: `false`

  Returns a vector of cluster maps (`:representative`, `:size`, `:indices` and `:data`), see the namespace documentation. `:size` and `:data` describe the hard assignment (the highest membership). With `:add-memberships?`, the cluster map contains `:memberships`, the vector with the membership of every input point in the cluster (in the order of input). The vector has the metadata `:objective`, the value of the objective function.

  A cluster that gets no point (eg. all points are equal and there is more than one cluster) has `:size` `0`, empty `:indices` and `:data`, memberships `0.0` and the representative `[##NaN ##NaN ...]`: its centre is a weighted mean with zero weights. Check `:size` before using the representatives of such a result. [[predict]], [[wcss]] and [[silhouette]] throw `ex-info` for such a result.

  Throws `ex-info` for an invalid option or data.

  See also [[kmeans++]], [[dbscan]], [[predict]], [[wcss]], [[silhouette]]."
  ([xss] (fuzzy-kmeans xss nil))
  ([xss {:keys [clusters fuzziness max-iters distance epsilon rng add-data? add-memberships?]
         :or {clusters 1 fuzziness 2.0 max-iters -1 distance :euclidean epsilon 1.0e-3
              rng (r/rng :default) add-data? true add-memberships? false}}]
   (let [points (->points xss)
         k (check-clusters clusters points)
         fuzziness (check-double :fuzziness fuzziness 1.0 false "a finite number greater than 1")
         max-iters (->int (check-integer :max-iters max-iters Long/MIN_VALUE "an integer"))
         epsilon (check-double :epsilon epsilon 0.0 false "a finite number greater than 0")
         dist (->distance distance)
         rng (check-rng rng)
         clusterer (FuzzyKMeansClusterer. (int k) fuzziness (int max-iters) dist epsilon rng)
         clusters (->clusters points clusterer add-data?)]
     (if add-memberships?
       (let [^RealMatrix memberships (.getMembershipMatrix clusterer)]
         (with-meta (into [] (map-indexed (fn [^long j cluster]
                                            (assoc cluster :memberships (vec (.getColumn memberships j)))))
                          clusters)
           {:objective (.getObjectiveFunctionValue clusterer)}))
       (with-meta clusters {:objective (.getObjectiveFunctionValue clusterer)})))))

;; ## DBSCAN

(defn- select-radius
  "DBSCAN radius for `points` and `k` neighbours with `method` (`:knee` or `:pairs`, validated by the caller)."
  [points ^DistanceMeasure dist k method sensitivity rng]
  (let [n (count points)
        k (long k)]
    (when (m/< n 3)
      (invalid-option :xss n "at least 3 points to infer the radius"))
    (case method
      :knee (do (when (m/> k (m/dec n))
                  (invalid-option :points k "an integer, at most the number of points minus 1"))
                (knee-radius points dist k (double sensitivity)))
      :pairs (let [arrays (points->arrays points)
                   distances (vec (for [i (range n)
                                        j (range (m/long-inc (long i)) n)]
                                    (.compute dist ^doubles (arrays i) ^doubles (arrays j))))]
               (when (m/< (count (take 2 (distinct (map #(m/+ 0.0 (double %)) distances)))) 2)
                 (invalid-option :eps nil "an explicit :eps, all distances between points are equal"))
               (->> (kmeans++ distances {:clusters 2 :add-data? false :rng rng})
                    (map (comp first :representative))
                    (reduce m/min))))))

(def ^:private radius-methods #{:knee :pairs})

(defn- default-min-points
  ^long [points points-option]
  (if (nil? points-option)
    (m/* 2 (point-dim points))
    (check-integer :points points-option 1 "an integer, at least 1, or nil")))

(defn infer-dbscan-radius
  "Infers the `:eps` radius for [[dbscan]] from data.

  Parameters:

  - `xss` - a sequence of at least 3 points, see the namespace documentation
  - `dist` - distance (keyword, function or `DistanceMeasure`)
  - options (a map):
    - `:method` - `:knee` (default) or `:pairs`
    - `:points` - the number of neighbours (as `:points` in [[dbscan]]), an integer at least 1, default: twice the dimension
    - `:sensitivity` - sensitivity of the knee detection, a finite number not smaller than 0, default: `1.0` (`:knee` only)
    - `:rng` - random number generator (`:pairs` only), default: a new default generator for every call

  The `:knee` method sorts the distances from every point to its `:points`-th nearest neighbour and returns the distance at the knee of this curve (the Kneedle algorithm, Satopaa et al. 2011, for an increasing convex curve). It needs O(n^2) time and O(n) memory. It throws `ex-info` when no knee is found (eg. for uniform data), use an explicit `:eps` then. The distance has to be 0 for equal points.

  The `:pairs` method is a heuristic: kmeans++ with two clusters on all pairwise distances, the lower centre is the radius. It needs O(n^2) memory and depends on the generator. It throws `ex-info` when all distances between points are equal (eg. identical points).

  Returns the radius as a double. Throws `ex-info` for an invalid option or data.

  See also [[dbscan]]."
  ([xss dist] (infer-dbscan-radius xss dist nil))
  ([xss dist {:keys [method points sensitivity rng]
              :or {method :knee sensitivity 1.0 rng (r/rng :default)}}]
   (let [pts (->points xss)
         dist (->distance dist)
         k (default-min-points pts points)
         method (check-member :method method radius-methods)
         sensitivity (check-double :sensitivity sensitivity 0.0 true "a finite number, at least 0")]
     (select-radius pts dist k method sensitivity (check-rng rng)))))

(defn dbscan
  "DBSCAN density-based clustering: points with at least `:points` neighbours within `:eps` start (and extend) clusters, other points are noise.

  Parameters:

  - `xss` - a non-empty sequence of points (sequences, arrays or numbers), see the namespace documentation
  - options (a map):
    - `:eps` - the radius of a neighbourhood (inclusive), a finite number greater than 0 or `nil` to select it, default: `nil`
    - `:eps-method` - how to select `:eps` when it is `nil`: `:knee` (default) or `:pairs`, see [[infer-dbscan-radius]]. `:pairs` uses a new default generator, so it is not reproducible: call [[infer-dbscan-radius]] with `:rng` and pass the result as `:eps`
    - `:points` - the number of neighbours (the point itself is not counted) required in the neighbourhood of a core point, an integer at least 1, default: twice the dimension
    - `:distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`
    - `:add-data?` - add `:data` to the clusters, default: `true`
    - `:noise?` - add the noise points as the last map, default: `false`

  Returns a vector of cluster maps (`:representative`, `:size`, `:indices` and `:data`), see the namespace documentation. The `:representative` is the mean of the points of a cluster. With `:noise?` the last element is `{:noise? true :size :indices}` (and `:data`), it does not have `:representative` and it is present even when there is no noise. Cluster and noise indices are disjoint and cover all input points.

  Without `:eps` the radius is inferred, which needs at least 3 points and throws when no knee is found. The number of neighbours of the Scikit-learn `min_samples` is `:points` plus one.

  Throws `ex-info` for an invalid option or data.

  See also [[kmeans++]], [[infer-dbscan-radius]], [[predict]], [[wcss]], [[silhouette]]."
  ([xss] (dbscan xss nil))
  ([xss {:keys [eps eps-method points distance add-data? noise?]
         :or {eps-method :knee distance :euclidean add-data? true noise? false}}]
   (let [pts (->points xss)
         k (default-min-points pts points)
         dist (->distance distance)
         eps-method (check-member :eps-method eps-method radius-methods)
         radius (if (nil? eps)
                  (select-radius pts dist k eps-method 1.0 (r/rng :default))
                  (check-double :eps eps 0.0 false "a finite number greater than 0, or nil"))
         clusters (->clusters pts (DBSCANClusterer. radius (int (m/min k Integer/MAX_VALUE)) dist) add-data?)]
     (if noise?
       (let [clustered (into #{} (mapcat :indices) clusters)
             indices (filterv (complement clustered) (range (count pts)))
             noise {:noise? true :size (count indices) :indices indices}]
         (conj clusters (if add-data?
                          (assoc noise :data (mapv #(point->vec (nth pts %)) indices))
                          noise)))
       clusters))))

;; ## Results

(defn- result-clusters
  "Cluster maps of a result without the noise map."
  [result]
  (let [clusters (when (sequential? result) (into [] (remove :noise?) result))]
    (if (and (seq clusters) (every? #(and (map? %) (:representative %)) clusters))
      clusters
      (invalid-option :result nil "a result of kmeans++, fuzzy-kmeans or dbscan"))))

(defn predict
  "Finds the cluster of a point: the index of the nearest representative.

  Parameters:

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]] or [[dbscan]] (the noise map is ignored)
  - `x` - a point of the dimension of the representatives
  - `distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`; use the one used for the clustering

  Returns the position of the cluster in `result` (ties go to the first cluster). It does not return noise: for [[dbscan]] the representatives are means of the clusters, so the result is meaningful for compact (convex) clusters only.

  Throws `ex-info` for an invalid result or a point of another dimension.

  See also [[wcss]], [[silhouette]]."
  ([result x] (predict result x :euclidean))
  ([result x distance]
   (let [clusters (result-clusters result)
         dist (->distance distance)
         p (point-array :x x)
         reps (mapv #(point-array :result (:representative %)) clusters)]
     (when-not (every? #(m/== (alength p) (alength ^doubles %)) reps)
       (invalid-option :x x "a point of the dimension of the representatives"))
     (loop [i 0 best 0 best-d Double/POSITIVE_INFINITY]
       (if (m/< i (count reps))
         (let [d (.compute dist p ^doubles (nth reps i))]
           (if (m/< d best-d)
             (recur (m/inc i) i d)
             (recur (m/inc i) best best-d)))
         best)))))

(defn wcss
  "Within-cluster sum of squares: the sum over all points of the squared distance to the representative of the cluster.

  Parameters:

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]] or [[dbscan]] with `:data` (the noise map is ignored)
  - `distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`

  Returns a double, lower is better for the same number of clusters.

  Throws `ex-info` for an invalid result or a result without `:data` (`:add-data?` set to `false`).

  See also [[silhouette]], [[predict]]."
  ([result] (wcss result :euclidean))
  ([result distance]
   (let [clusters (result-clusters result)
         dist (->distance distance)]
     (when-not (every? :data clusters)
       (invalid-option :result nil "a result with :data (:add-data? true)"))
     (reduce (fn ^double [^double acc cluster]
               (let [rep (point-array :result (:representative cluster))]
                 (reduce (fn ^double [^double acc x]
                           (m/+ acc (m/sq (.compute dist (point-array :result x) rep))))
                         acc (:data cluster))))
             0.0 clusters))))

(defn- distance-sum
  "Sum of distances from the point `p` to `arrays`."
  ^double [^DistanceMeasure dist ^doubles p arrays]
  (reduce (fn ^double [^double acc a] (m/+ acc (.compute dist p ^doubles a))) 0.0 arrays))

(defn- silhouette-score
  "Silhouette of a point of the cluster `ci`.

  `sums` are the sums of distances from the point to the points of every cluster (the point itself is at the distance 0), `sizes` are the sizes of clusters."
  ^double [sums sizes ^long ci]
  (let [size (long (nth sizes ci))]
    (if (m/== size 1)
      0.0
      (let [a (m// (double (nth sums ci)) (double (m/long-dec size)))
            b (double (reduce m/min (keep-indexed (fn [c s] (when-not (m/== (long c) ci)
                                                             (m// (double s) (double (nth sizes c)))))
                                                  sums)))
            den (m/max a b)]
        (if (m/zero? den) 0.0 (m// (m/- b a) den))))))

(defn silhouette
  "Silhouette of a clustering: how well points fit their cluster compared to the nearest other cluster.

  For a point, `a` is the mean distance to other points of its cluster and `b` is the lowest mean distance to the points of another cluster, the silhouette is `(b - a) / max(a, b)`. It is in the range from -1 to 1 (higher is better), `0` for a point in a one-point cluster.

  Parameters:

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]] or [[dbscan]] with `:data` and at least 2 clusters (the noise map is ignored)
  - `distance` - distance (keyword, function or `DistanceMeasure`), default: `:euclidean`

  Returns a map:

  - `:mean` - the mean silhouette of all clustered points
  - `:samples` - silhouettes of clustered points, ordered by position in the input
  - `:indices` - positions of the points in `:samples` in the input sequence (noise points are not included)

  It needs O(n^2) time and O(n) memory.

  Throws `ex-info` for an invalid result, a result without `:data` (`:add-data?` set to `false`), with less than 2 clusters or with an empty cluster (see [[fuzzy-kmeans]]). The distance has to be 0 for equal points.

  See also [[wcss]], [[predict]]."
  ([result] (silhouette result :euclidean))
  ([result distance]
   (let [clusters (result-clusters result)]
     (when (m/< (count clusters) 2)
       (invalid-option :result nil "a result with at least 2 clusters"))
     (when-not (every? #(and (:data %) (:indices %)) clusters)
       (invalid-option :result nil "a result with :data (:add-data? true)"))
     (when (some (comp empty? :data) clusters)
       (invalid-option :result nil "a result without empty clusters"))
     (let [dist (->distance distance)
           arrays (mapv (fn [cluster] (mapv #(point-array :result %) (:data cluster))) clusters)
           sizes (mapv count arrays)
           ;; [index-in-input cluster-position point] for every clustered point
           rows (mapcat (fn [ci cluster cluster-arrays] (map (fn [idx p] [idx ci p]) (:indices cluster) cluster-arrays))
                        (range) clusters arrays)
           scored (sort-by first (map (fn [[idx ci p]]
                                        [idx (silhouette-score (mapv #(distance-sum dist p %) arrays) sizes ci)])
                                      rows))
           samples (mapv second scored)]
       {:mean (m// (v/sum samples) (count samples))
        :samples samples
        :indices (mapv first scored)}))))
