(ns fastmath.clustering
  "Clustering algorithms.

  Various clustering algrorithms backed by SMILE library.

  Currently implemented: only partition clustering.
  
  ### Input data

  It's always sequence of n-sized samples as sequences.

  For example, 2d samples `[[1 2] [2 2] [3 3] ...]`

  For 1d data you can pass sequence of numbers of sequence of 1d seqs of numbers

  ```clojure
  [1 2 3]
  ;; or
  [[1] [2] [3]]
  ```

  ### Distances

  Some of the methods use distance functions, use [[fastmath.distance]] namespace to create one.
  
  ### Output

  Every function returns record which contains:

  * `:type` - name of the method used
  * `:data` - input data
  * `:clustering` - sequence of cluster ids
  * `:sizes` - sizes of clusters
  * `:clusters` - number of clusters
  * `:predict` - predicting function (see below), qualify additional sample
  * `:representatives` - list of centroids or medoids if available
  * `:info` - additional statistics for your samples (like distortion)
  * `:obj` - SMILE object

  Cluster id is a integer ranging from 0 to the number of clusters minus 1. Some methods mark outliers with [[outlier-id]].
  
  Record acts as function and can qualify additonal sample by calling `:predict` function (or just call [[predict]]), for example (`data` is sequence of 3d samples):

  ```clojure
  (let [cl (k-means data 10)] (cl [0 1 2]))
  ```

  See [[k-means]]
  
  #### Regrouping

  Clustering record can be regroupped to the list of individual clusters. Call [[regroup]] and get list of maps with following structure:
  
  * `:key` - cluster id
  * `:data` - samples which belong to the cluster
  * `:outliers?` - does it contain outliers or not
  * `:representative` - centroid/medoid or average vector if the former is not available
  * `:size` - size of cluster"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.string :as s]
            [fastmath.stats :as stat]
            [fastmath.distance :as d]
            [fastmath.protocols :as prot])
  (:import [smile.clustering PartitionClustering KMeans GMeans XMeans DeterministicAnnealing
            DENCLUE CLARANS DBSCAN MEC]
           [smile.vq NeuralGas]
           [clojure.lang IFn]))

(set! *warn-on-reflection* false)

(def ^:const ^{:doc "Id of the cluster which contain outliers."} outlier-id PartitionClustering/OUTLIER)

(defrecord ClusteringResult [type data clustering sizes clusters predict representatives info obj]
  IFn
  (invoke [_ in] (predict in)))

(defn predict
  "Predict cluster for given vector"
  [^ClusteringResult cluster in]
  ((.predict cluster) in))

(defn- structurize 
  "Pack result of the clustering function into the structure"
  [type data in repr info]
  (->ClusteringResult type
                      data
                      (seq (.getClusterLabel ^PartitionClustering in))
                      (seq (.getClusterSize ^PartitionClustering in))
                      (.getNumClusters ^PartitionClustering in)
                      (fn [x] (.predict ^PartitionClustering in (m/seq->double-array x)))
                      (when repr (m/double-double-array->seq (repr in)))
                      (into {} (map (fn [[k v]] [k (v in)]) info))
                      in))

(def ^:private clustering-classes {:k-means ['KMeans 'centroids 'distortion]
                                   :g-means ['GMeans 'centroids 'distortion]
                                   :x-means ['XMeans 'centroids 'distortion]
                                   :deterministic-annealing ['DeterministicAnnealing 'centroids 'distortion 'getAlpha]
                                   :neural-gas ['NeuralGas 'centroids 'distortion]
                                   :denclue ['DENCLUE nil 'getSigma]
                                   :clarans ['CLARANS 'medoids 'distortion 'getMaxNeighbor 'getNumLocalMinima]
                                   :dbscan ['DBSCAN nil 'getMinPts 'getRadius]
                                   :mec ['MEC nil 'entropy]})

(def ^{:doc "List of clustering methods."} clustering-methods-list (keys clustering-classes))

(defn- symbol->keyword
  "Convert java method into the keyword."
  [s]
  (let [s (name s)]
    (-> (if (= "get" (subs s 0 3)) (subs s 3) s)
        (s/lower-case)
        (keyword))))

(defmacro ^:private clustering
  "Analyze clustering method and pack into the structure."
  [clustering-method data & params]
  (let [[nm repr & rest] (clustering-classes clustering-method)
        obj (gensym "obj")
        repr (when repr `(fn [~obj] (. ~obj ~repr)))
        info (mapv #(vector (symbol->keyword %) `(fn [~obj] (. ~obj ~%))) rest)]
    `(structurize ~clustering-method ~data (new ~nm (m/seq->double-double-array ~data) ~@params) ~repr ~info)))

(defn k-means
  "K-Means++ algorithm.

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * max-iter (optional) - maximum number of iterations
  * runs (optional) - maximum number of runs

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/KMeans.html)"
  ([data clusters] (clustering :k-means data clusters))
  ([data clusters max-iter] (clustering :k-means data clusters max-iter))
  ([data clusters max-iter runs] (clustering :k-means data clusters max-iter runs)))

(defn g-means
  "G-Means algorithm.

  Input:

  * data - sequence of samples
  * max-clusters - maximum number of clusters

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/GMeans.html)"
  [data max-clusters] (clustering :g-means data max-clusters))

(defn x-means
  "X-Means algorithm.

  Input:

  * data - sequence of samples
  * max-clusters - number of clusters

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/XMeans.html)"
  [data max-clusters] (clustering :x-means data max-clusters))

(defn deterministic-annealing
  "Deterministic Annealing algorithm.

  Input:

  * data - sequence of samples
  * max-clusters - number of clusters
  * alpha (optional) - temperature decreasing factor (valued from 0 to 1)

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/DeterministicAnnealing.html)"
  ([data max-clusters] (clustering :deterministic-annealing data max-clusters))
  ([data max-clusters alpha] (clustering :deterministic-annealing data max-clusters alpha)))

(defn neural-gas
  "Neural Gas algorithm.

  Input:

  * data - sequence of samples
  * clusters - number of clusters

  Optional:

  * lambda-i - intial lambda value (soft learning radius/rate)
  * lambda-f - final lambda value
  * eps-i - initial epsilon value (learning rate)
  * eps-f - final epsilon value
  * steps - number of iterations

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/vq/NeuralGas.html)"
  ([data clusters] (clustering :neural-gas data clusters))
  ([data clusters lambda-i lambda-f eps-i eps-f steps] (clustering :neural-gas data clusters lambda-i lambda-f eps-i eps-f steps)))

(defn denclue
  "DENsity CLUstering algorithm.

  Input:

  * data - sequence of samples
  * sigma - gaussian kernel parameter
  * m - number of selected samples, much slower than number of all samples

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/DENCLUE.html)"
  [data sigma m] (clustering :denclue data sigma m))

(defn clarans
  "Clustering Large Applications based upon RANdomized Search algorithm.

  Input:

  * data - sequence of samples
  * clusters - numbe of clusters

  Optional:

  * dist - distance method, default `:euclidean`
  * max-neighbor - maximum number of neighbors checked during random search
  * num-local - the number of local minima to search for

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/CLARANS.html)"
  ([data clusters] (clarans data d/euclidean clusters))
  ([data dist clusters] (clustering :clarans data dist clusters))
  ([data dist clusters max-neighbor] (clustering :clarans data dist clusters max-neighbor))
  ([data dist clusters max-neighbor num-local] (clustering :clarans data dist clusters max-neighbor num-local)))

(defn dbscan
  "Density-Based Spatial Clustering of Applications with Noise algorithm.

  Input:

  * data - sequence of samples
  * dist (optional) - distance method, default `:euclidean`
  * min-pts - minimum number of neighbors
  * radius - the neighborhood radius

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/DBSCAN.html)"
  ([data min-pts radius] (dbscan data d/euclidean min-pts radius))
  ([data dist min-pts ^double radius] (clustering :dbscan data dist (int min-pts) radius)))

(defn mec
  "Nonparametric Minimum Conditional Entropy Clustering algorithm.

  Input:

  * data - sequence of samples
  * dist (optional) - distance method, default `:euclidean`
  * max-clusters - maximum number of clusters
  * radius - the neighborhood radius

  See more in [SMILE doc](https://haifengl.github.io/smile/api/java/smile/clustering/MEC.html)"
  ([data max-clusters radius] (mec data d/euclidean max-clusters radius))
  ([data dist max-clusters ^double radius] (clustering :mec data dist max-clusters radius)))

(defn regroup
  "Transform clustering result into list of clusters as separate maps.

  Every map contain:

  * `:key` - cluster id
  * `:data` - samples which belong to the cluster
  * `:outliers?` - does it contain outliers or not
  * `:representative` - centroid/medoid or average vector if the former is not available
  * `:size` - size of cluster

  Representative is always a n-dimensional sequence even if input is a list of numbers.

  Empty clusters are skipped."
  [clustered-data]
  (let [mvector? (satisfies? prot/VectorProto (first (:data clustered-data))) ;;required to fix missing representative
        mseqable? (sequential? (first (:data clustered-data)))]
    (for [[k lst] (group-by first (map vector (:clustering clustered-data) (:data clustered-data)))
          :let [d (map second lst)
                outliers? (== k outlier-id)
                r (:representatives clustered-data)]]
      {:key k
       :outliers? outliers?
       :data d
       :representative (if (and r (not outliers?))
                         (nth r k)
                         (cond ;; if representative is missing, calculate average
                           mvector? (v/average-vectors d)
                           mseqable? (v/average-vectors (map vec d))
                           :else (list (stat/mean d))))
       :size (if outliers?
               (count d)
               (nth (:sizes clustered-data) k))})))
