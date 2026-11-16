(ns fastmath.ml.clustering
  "Clustering of points: k-means++, fuzzy k-means, DBSCAN and hierarchical clustering, with cluster validation helpers.

  K-means++, fuzzy k-means and DBSCAN are backed by Apache Commons Math. Points are sequences of numbers (or arrays); a plain number is a 1-dimensional point. All points have to have the same, non-zero dimension and finite coordinates.

  Every clustering function returns a vector of cluster maps ([[hclust]] returns it as `:clusters`). A cluster map contains:

  - `:representative` - the centre of the cluster (centroid for k-means++ and fuzzy k-means, the mean of cluster points for DBSCAN and hierarchical clustering),
  - `:size` - the number of points,
  - `:indices` - positions of the cluster points in the input sequence (the order within a cluster is not specified),
  - `:data` - the cluster points (vectors of doubles), unless `:add-data?` is `false`.

  The order of clusters and the choice of representatives of k-means++, fuzzy k-means and DBSCAN depend on the random number generator. The default `:rng` is a new generator for every call, pass a seeded generator (`fastmath.random/rng`) for reproducible results.

  Options are validated: an invalid option throws `ex-info` with the keys `:option`, `:value` and `:expected`. `:value` is the invalid value; it is `nil` when the value is the whole input or a whole result (data can be large), and a single point when the points differ in dimension. Unknown option keys are ignored.

  `:distance` is a keyword (`:euclidean`, `:chebyshev`, `:manhattan`, `:canberra`, `:earth-movers`), a function of two points (ie. a function from `fastmath.distance`) or an Apache Commons Math `DistanceMeasure` ([[hclust]] takes a keyword or a function only, or `:precomputed` for a distance matrix). A distance has to be non-negative and equal to 0 for equal points (`d(x, x) = 0`): the radius selection of [[dbscan]] and [[silhouette]] rely on it and return wrong values otherwise.

  Hierarchical clustering ([[hclust]]) returns a map with the cluster maps (`:clusters`) and the tree of merges (`:tree`, a dendrogram in the format of R `hclust`); [[cut]] cuts the tree into another number of clusters, [[cophenetic]], [[cophenetic-correlation]] and [[->linkage-matrix]] work on the tree.

  Functions: [[kmeans++]], [[fuzzy-kmeans]], [[dbscan]], [[infer-dbscan-radius]], [[hclust]] (clustering and the radius selection), [[predict]], [[wcss]], [[silhouette]] (working on a clustering result), [[cut]], [[cophenetic]], [[cophenetic-correlation]], [[->linkage-matrix]] (working on a result of [[hclust]])."
  (:require [fastmath.distance :as dists]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.ml.clustering.hierarchical :as hier])
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
           [org.apache.commons.math3.random RandomGenerator]
           [fastmath.java Array]))

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
      (invalid-option :result nil "a result of kmeans++, fuzzy-kmeans, dbscan or the :clusters of hclust (on points)"))))

(defn predict
  "Finds the cluster of a point: the index of the nearest representative.

  Parameters:

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]], [[dbscan]] or the `:clusters` of [[hclust]] on points (the noise map is ignored)
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

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]], [[dbscan]] or the `:clusters` of [[hclust]] on points with `:data` (the noise map is ignored)
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

  - `result` - a result of [[kmeans++]], [[fuzzy-kmeans]], [[dbscan]] or the `:clusters` of [[hclust]] on points with `:data` and at least 2 clusters (the noise map is ignored)
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

;; ## Hierarchical clustering

(def ^{:private true :const true :tag 'long} max-observations 65536)

(def ^:private hclust-methods #{:single :complete :average :weighted :centroid :median :ward :ward-d :diana})
(def ^:private hclust-algorithms #{:generic :fast})

(def ^:private point-distances
  {:euclidean dists/euclidean
   :chebyshev dists/chebyshev
   :manhattan dists/manhattan
   :canberra dists/canberra
   :earth-movers dists/earth-movers})

(defn- ->point-distance
  "Resolve a distance keyword or a function of two points to a function of two points."
  [distance]
  (let [expected (str "one of " (sort (keys point-distances)) ", :precomputed or a function of two points")]
    (cond
      (keyword? distance) (or (point-distances distance) (invalid-option :distance distance expected))
      (fn? distance) distance
      :else (invalid-option :distance distance expected))))

(defn- check-observation-count
  "Throws when `xss` has more than 65536 elements, before any element is converted."
  [xss]
  (when (and (seqable? xss) (not (string? xss)) (not (map? xss)))
    (let [n (long (bounded-count (m/inc max-observations) xss))]
      (when (m/> n max-observations)
        (invalid-option :xss nil "at most 65536 observations")))))

(defn- matrix-row
  ^doubles [row expected]
  (let [arr (try (v/vec->array row) (catch Exception _ nil))]
    (if (instance? (Class/forName "[D") arr) arr (invalid-option :xss nil expected))))

(defn- ->distance-rows
  "Validate a square distance matrix and return its rows as double arrays."
  [xss]
  (let [expected "a square matrix of finite, non-negative and symmetric distances with zeros on the diagonal"]
    (when-not (and (sequential? xss) (seq xss))
      (invalid-option :xss nil expected))
    (let [rows (mapv #(matrix-row % expected) xss)
          n (count rows)
          entry (fn ^double [^long i ^long j] (Array/aget ^doubles (rows i) j))]
      (when-not (and (every? #(m/== n (alength ^doubles %)) rows)
                     (every? (fn [^doubles row] (every? #(and (m/valid-double? %) (m/>= (double %) 0.0)) row)) rows)
                     (every? (fn [^long i] (m/zero? (Array/aget ^doubles (rows i) i))) (range n))
                     (every? (fn [[^long i ^long j]] (m/delta-eq (entry i j) (entry j i)))
                             (for [i (range n) j (range i)] [i j])))
        (invalid-option :xss nil expected))
      rows)))

(defn- hierarchy-input
  "Observations of `xss` for the `distance` option: their number, the condensed distances (a function, every call gives a new array), the distance between two observations and the points (`nil` for a distance matrix)."
  [xss distance]
  (check-observation-count xss)
  (if (= distance :precomputed)
    (let [rows (->distance-rows xss)]
      {:n (count rows)
       :condensed (fn [] (hier/matrix-distances rows))
       :pair-distance (fn ^double [^long i ^long j] (Array/aget ^doubles (rows i) j))})
    (let [dist (->point-distance distance)
          arrays (points->arrays (->points xss))]
      {:n (count arrays)
       :arrays arrays
       :condensed (fn [] (hier/data-distances arrays dist))
       :pair-distance (fn ^double [^long i ^long j] (hier/point-distance dist (arrays i) (arrays j)))})))

(defn- check-cut-options
  "Validates `:k` and `:h` (not both, `:k` defaults to 1)."
  [k h]
  (cond
    (and (some? k) (some? h)) (invalid-option :k k "either :k or :h, not both")
    (some? h) {:h (check-double :h h 0.0 true "a finite number, at least 0")}
    :else {:k (check-integer :k (or k 1) 1 "an integer between 1 and the number of observations")}))

(defn- cluster-count
  "The number of clusters for the validated cut options of a tree of `n` observations."
  ^long [{:keys [k h]} {:keys [height]} ^long n]
  (if (some? h)
    (if (every? true? (map <= height (rest height)))
      (hier/clusters-for-height height n h)
      (invalid-option :h h "a hierarchy with non-decreasing heights (no inversions), use :k"))
    (if (m/> (long k) n)
      (invalid-option :k k "an integer between 1 and the number of observations")
      (long k))))

(defn- tree-clusters
  "Cluster maps of the cut of `tree` into `k` clusters."
  [tree {:keys [n arrays pair-distance]} ^long k add-data?]
  (mapv (fn [indices]
          (let [res {:medoid (hier/medoid indices pair-distance)
                     :size (count indices)
                     :indices indices}]
            (if arrays
              (let [data (mapv #(vec (arrays %)) indices)
                    res (assoc res :representative (vec (v/average-vectors data)))]
                (if add-data? (assoc res :data data) res))
              res)))
        (hier/clusters-at (:merge tree) n k)))

(defn- hclust-result
  "A result of [[hclust]]: the clusters of the cut of `tree` and the tree."
  [tree input cut-options add-data?]
  (let [k (cluster-count cut-options tree (:n input))]
    {:clusters (tree-clusters tree input k add-data?)
     :tree tree}))

(defn- result-tree
  "The tree of a result of [[hclust]]."
  [result]
  (let [{:keys [merge height order n] :as tree} (when (map? result) (:tree result))]
    (if (and (vector? (:clusters result)) (vector? merge) (vector? height) (vector? order) (int? n))
      tree
      (invalid-option :result nil "a result of hclust (a map with :clusters and :tree)"))))

(defn- check-input-size
  [input ^long n]
  (when-not (m/== n (long (:n input)))
    (invalid-option :xss nil (str n " observations, as clustered in the result"))))

(defn hclust
  "Hierarchical clustering: builds a tree of nested clusters (a dendrogram) and cuts it into flat clusters.

  Agglomerative methods start with every observation as a cluster and repeatedly merge the two closest clusters, the distance between clusters is defined by the linkage `:method`. The divisive method (`:diana`, DIANA) starts with one cluster and splits clusters until every observation is alone: the most distant observation starts a splinter group and observations closer to it on average move to it; the height of a split is the diameter of the split cluster. The result follows R `hclust`: the merge heights and the encoding of merges are the same.

  Parameters:

  - `xss` - a non-empty sequence of points (sequences, arrays or numbers, see the namespace documentation) or, with `:distance` `:precomputed`, a square matrix (a sequence of rows) of distances between observations: finite, non-negative, symmetric, with zeros on the diagonal. At most 65536 observations.
  - options (a map):
    - `:method` - the linkage, default: `:complete`
      - `:single` - the shortest distance between points of two clusters
      - `:complete` - the longest distance between points of two clusters
      - `:average` - the mean distance between points of two clusters (UPGMA)
      - `:weighted` - the mean of distances of the two merged clusters (WPGMA, R `mcquitty`)
      - `:centroid` - the distance between centroids (UPGMC), as SciPy `centroid` (R `centroid` on squared distances, the heights are square roots)
      - `:median` - the distance between weighted centroids (WPGMC), as SciPy `median` (R `median` on squared distances, the heights are square roots)
      - `:ward` - the increase of the within-cluster variance, Ward's method on squared distances (R `ward.D2`, SciPy `ward`)
      - `:ward-d` - Ward's update applied to distances, not squared (R `ward.D`)
      - `:diana` - divisive analysis (R `cluster::diana`)
    - `:distance` - a keyword (`:euclidean`, `:chebyshev`, `:manhattan`, `:canberra`, `:earth-movers`), a function of two points (ie. a function from `fastmath.distance`, it gets double arrays) or `:precomputed` for a distance matrix, default: `:euclidean`
    - `:algorithm` - `:generic` (default) or `:fast` (agglomerative methods only, it throws for `:diana`), see below
    - `:k` - the number of clusters, an integer between 1 and the number of observations, default: `1`
    - `:h` - the height of the cut, a finite number not smaller than 0: merges at a height not greater than `:h` are kept. Only one of `:k` and `:h` can be given.
    - `:add-data?` - add `:data` to the clusters, default: `true`

  Any distance works with any method. The `:centroid`, `:median` and `:ward` methods assume the Euclidean distance: with another distance the heights have no geometric meaning, and `:centroid` or `:median` can give `NaN` heights (a negative squared distance).

  `:algorithm` `:generic` reproduces R `hclust` exactly, also for ties. `:fast` uses the minimum spanning tree (`:single`) or the nearest-neighbour chain (`:complete`, `:average`, `:weighted`, `:ward`, `:ward-d`), as SciPy, and is faster for many observations. Without ties the tree is the same and the heights agree up to rounding; with ties another tree with the same heights can be returned. `:centroid` and `:median` always use the generic algorithm.

  Returns a map:

  - `:clusters` - a vector of cluster maps, ordered by the first observation of every cluster (as R `cutree`)
  - `:tree` - the tree of merges (a dendrogram)

  A cluster map contains:

  - `:indices` - positions of the observations in `xss`, ascending
  - `:size` - the number of observations
  - `:medoid` - the position of the observation with the lowest sum of distances to the other observations of the cluster (the lowest position on ties)
  - `:representative` - the mean of the points of the cluster (only for points)
  - `:data` - the points of the cluster (only for points, unless `:add-data?` is `false`)

  The tree is a map:

  - `:merge` - `n-1` merges `[a b]` as in R `hclust`: a negative number `-i` is the observation at the position `i-1` of `xss`, a positive number `j` is the cluster made by the merge `j` (counted from 1). A single observation goes first, two clusters are ordered.
  - `:height` - `n-1` heights of the merges, non-decreasing except for `:centroid` and `:median` (inversions)
  - `:order` - positions of the observations in the order of the leaves of the dendrogram (for plotting)
  - `:n` - the number of observations
  - `:method`, `:distance`, `:algorithm` - the options used (no `:algorithm` for `:diana`)
  - `:dc` - the divisive coefficient (`:diana` only)

  The `:clusters` of a result of points work with [[predict]], [[wcss]] and [[silhouette]]; for a distance matrix they have no `:representative` and `:data`, so these functions throw.

  It needs O(n^2) memory, the time is between O(n^2) and O(n^3).

  Throws `ex-info` for an invalid option or data, for `:h` with inversions, and for `:algorithm` with `:diana`.

  See also [[cut]], [[cophenetic]], [[cophenetic-correlation]], [[->linkage-matrix]]."
  ([xss] (hclust xss nil))
  ([xss {:keys [method distance algorithm k h add-data?]}]
   (let [method (check-member :method (or method :complete) hclust-methods)
         diana? (= method :diana)
         _ (when (and diana? (some? algorithm))
             (invalid-option :algorithm algorithm "no :algorithm for :diana"))
         algorithm (when-not diana? (check-member :algorithm (or algorithm :generic) hclust-algorithms))
         distance (or distance :euclidean)
         cut-options (check-cut-options k h)
         input (hierarchy-input xss distance)
         n (long (:n input))
         d ((:condensed input))
         tree (cond-> (assoc (if diana? (hier/diana d n) (hier/agglomerate d n method algorithm))
                             :n n :method method :distance distance)
                algorithm (assoc :algorithm algorithm))]
     (hclust-result tree input cut-options (if (nil? add-data?) true add-data?)))))

(defn cut
  "Cuts the tree of a result of [[hclust]] into another number of clusters, without clustering again.

  Parameters:

  - `result` - a result of [[hclust]]
  - `xss` - the data clustered by [[hclust]] (points or the distance matrix), needed for `:representative`, `:medoid` and `:data`
  - options (a map):
    - `:k` - the number of clusters, an integer between 1 and the number of observations, default: `1`
    - `:h` - the height of the cut, a finite number not smaller than 0 (not for a tree with inversions)
    - `:add-data?` - add `:data` to the clusters, default: `true`

  Returns a new result, as [[hclust]], with the same `:tree`.

  Throws `ex-info` for an invalid result, options, or `xss` with another number of observations.

  See also [[hclust]]."
  ([result xss] (cut result xss nil))
  ([result xss {:keys [k h add-data?]}]
   (let [tree (result-tree result)
         cut-options (check-cut-options k h)
         input (hierarchy-input xss (:distance tree))]
     (check-input-size input (:n tree))
     (hclust-result tree input cut-options (if (nil? add-data?) true add-data?)))))

(defn cophenetic
  "Cophenetic distances of a result of [[hclust]]: the height of the merge which first joins two observations.

  Parameters:

  - `result` - a result of [[hclust]]

  Returns a square matrix (a vector of `n` vectors of doubles) with zeros on the diagonal, as R `cophenetic`.

  Throws `ex-info` for an invalid result.

  See also [[cophenetic-correlation]], [[hclust]]."
  [result]
  (let [{:keys [merge height n]} (result-tree result)
        n (long n)
        d (hier/cophenetic-distances merge height n)]
    (mapv (fn [i] (mapv (fn [j] (hier/distance-at d n i j)) (range n))) (range n))))

(defn cophenetic-correlation
  "Cophenetic correlation of a result of [[hclust]]: how well the tree preserves the distances between observations.

  The Pearson correlation between the cophenetic distances (see [[cophenetic]]) and the distances of all pairs of observations. A value close to 1 means the dendrogram is a faithful summary of the distances.

  Parameters:

  - `result` - a result of [[hclust]]
  - `xss` - the data clustered by [[hclust]] (points or the distance matrix)

  Returns a double, `NaN` for less than 3 observations or when the distances or the cophenetic distances are all equal.

  Throws `ex-info` for an invalid result or `xss` with another number of observations.

  See also [[cophenetic]], [[hclust]]."
  ^double [result xss]
  (let [{:keys [merge height n distance]} (result-tree result)
        n (long n)
        input (hierarchy-input xss distance)]
    (check-input-size input n)
    (if (m/< n 3)
      ##NaN
      (let [c (hier/cophenetic-distances merge height n)
            d ((:condensed input))]
        (if (or (apply m/== (seq c)) (apply m/== (seq d)))
          ##NaN
          (stats/pearson-correlation c d))))))

(defn ->linkage-matrix
  "Converts the tree of a result of [[hclust]] to the SciPy linkage matrix.

  Parameters:

  - `result` - a result of [[hclust]]

  Returns a vector of `n-1` rows `[a b height count]`: clusters `a < b` merged at `height` into a cluster of `count` observations. Clusters `0..n-1` are observations, the cluster made by the row `i` (counted from 0) is `n+i`.

  Throws `ex-info` for an invalid result.

  See also [[hclust]]."
  [result]
  (let [{:keys [merge height n]} (result-tree result)]
    (hier/linkage-matrix merge height n)))
