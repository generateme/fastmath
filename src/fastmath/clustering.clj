(ns fastmath.clustering
  "Clustering.

  Various clustering algrorithms backed by SMILE library.

  Only partition clustering is implemented.
  
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
  * `:representatives` - list of centroids or averages
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
  
  * `:key` - cluster id or `:outliers`
  * `:data` - samples which belong to the cluster
  * `:representative` - centroid or average vector if the former is not available
  * `:size` - size of cluster"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.distance :as d])
  (:import [smile.clustering CentroidClustering PartitionClustering KMeans GMeans XMeans DeterministicAnnealing
            DENCLUE CLARANS DBSCAN MEC SpectralClustering]
           [smile.math.distance Distance]
           [clojure.lang IFn]))

(def ^:const ^{:doc "Id of the cluster which contain outliers."} outlier-id PartitionClustering/OUTLIER)

(defrecord ClusteringResult [type data clustering sizes clusters predict representatives info obj]
  IFn
  (invoke [_ in] (predict in)))

(defn predict
  "Predict cluster for given vector"
  [^ClusteringResult cluster in]
  (when-let [pred (.predict cluster)] (pred in)))

(defn regroup
  "Transform clustering result into list of clusters as separate maps.

  Every map contain:

  * `:key` - cluster id or `:outliers`
  * `:data` - samples which belong to the cluster
  * `:representative` - centroid/medoid or average vector if the former is not available
  * `:size` - size of cluster

  Representative is always a n-dimensional sequence even if input is a list of numbers.

  Empty clusters are skipped."
  [{:keys [clustering data representatives sizes]}]
  (for [[k lst] (group-by first (map vector clustering data))
        :let [d (map second lst)
              outliers? (== k outlier-id)]]
    {:key (if outliers? :outliers k)
     :data d
     :representative (if (and representatives (not outliers?))
                       (nth representatives k)
                       (v/average-vectors d))
     :size (if outliers?
             (count d)
             (nth sizes k))}))

;;

(defn- centroid-data
  [^CentroidClustering in]
  {:representatives (m/double-double-array->seq (.centroids in))
   :distortion (.distortion in)})

(defn- mec-data [^MEC in] {:entropy (.entropy in)})
(defn- spectral-data [^SpectralClustering in] {:distortion (.distortion in)})

(defn- create-clusters
  [type data ^PartitionClustering clusters info-fn predict-fn]
  (let [additional-data (if info-fn (info-fn clusters) {})]
    (->ClusteringResult type
                        data
                        (seq (.y clusters))
                        (seq (.size clusters))
                        (.k clusters)
                        predict-fn
                        (:representatives additional-data)
                        (dissoc additional-data :representatives)
                        clusters)))

(def ^:private clustering-classes {:k-means '[KMeans centroid-data]
                                   :lloyd '[KMeans centroid-data lloyd]
                                   :g-means '[GMeans centroid-data]
                                   :x-means '[XMeans centroid-data]
                                   :deterministic-annealing '[DeterministicAnnealing centroid-data]
                                   :clarans '[CLARANS centroid-data]
                                   :denclue '[DENCLUE nil]
                                   :dbscan '[DBSCAN nil]
                                   :mec '[MEC mec-data]
                                   :spectral '[SpectralClustering spectral-data nil true]})

(def ^{:doc "List of clustering methods."} clustering-methods-list (keys clustering-classes))

(defmacro ^:private clustering
  "Analyze clustering method and pack into the structure."
  [clustering-method data & params]
  (let [[clss data-fn fit missing-predict?] (clustering-classes clustering-method)
        obj (with-meta (gensym "obj") {:tag clss})]
    `(let [~obj (. ~clss ~(or fit 'fit) (m/seq->double-double-array ~data) ~@params)]
       (create-clusters ~clustering-method ~data ~obj ~data-fn
                        ~(if missing-predict?
                           `nil
                           `(fn [in#] (.predict ~obj (m/seq->double-array in#))))))))

(defn k-means
  "K-Means++ algorithm.

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of convergence test

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/KMeans.html)"
  ([data clusters] (clustering :k-means data clusters))
  ([data clusters max-iter] (k-means data clusters max-iter 1.0e-4))
  ([data clusters max-iter tolerance] (clustering :k-means data clusters max-iter tolerance)))

(defn lloyd
  "K-Means algorithm, lloyd.

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of convergence test

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/KMeans.html)"
  ([data clusters] (clustering :lloyd data clusters))
  ([data clusters max-iter] (lloyd data clusters max-iter 1.0e-4))
  ([data clusters max-iter tolerance] (clustering :lloyd data clusters max-iter tolerance)))

(defn g-means
  "G-Means

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of convergence test

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/GMeans.html)"
  ([data clusters] (clustering :g-means data clusters))
  ([data clusters max-iter] (g-means data clusters max-iter 1.0e-4))
  ([data clusters max-iter tolerance] (clustering :g-means data clusters max-iter tolerance)))

(defn x-means
  "X-Means

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of convergence test

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/XMeans.html)"
  ([data clusters] (clustering :x-means data clusters))
  ([data clusters max-iter] (x-means data clusters max-iter 1.0e-4))
  ([data clusters max-iter tolerance] (clustering :x-means data clusters max-iter tolerance)))

(defn deterministic-annealing
  "Deterministic Annealing algorithm.

  Input:

  * data - sequence of samples
  * max-clusters - number of clusters
  * alpha (optional) - temperature decreasing factor (valued from 0 to 1)
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of convergence test
  * split-tolerance (optional) - tolerance to split a cluster
  
  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/DeterministicAnnealing.html)"
  ([data max-clusters] (clustering :deterministic-annealing data max-clusters))
  ([data max-clusters alpha] (deterministic-annealing data max-clusters alpha 100))
  ([data max-clusters alpha max-iter]
   (deterministic-annealing data max-clusters alpha max-iter 1.0e-4))
  ([data max-clusters alpha max-iter tolerance]
   (deterministic-annealing data max-clusters alpha max-iter tolerance 1.0e-2))
  ([data max-clusters alpha max-iter tolerance split-tolerance]
   (clustering :deterministic-annealing data max-clusters alpha max-iter tolerance split-tolerance)))

(defn clarans
  "Clustering Large Applications based upon RANdomized Search algorithm.

  Input:

  * data - sequence of samples
  * dist (optional) - distance method, default `euclidean`
  * clusters - number of clusters
  * max-neighbor (optional) - maximum number of neighbors checked during random search

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/CLARANS.html)"
  ([data clusters] (clarans data d/euclidean clusters))
  ([data dist ^long clusters] (clarans data dist clusters (max 1 (m/round (* 0.0125 clusters (- (count data) clusters))))))
  ([data dist clusters max-neighbor] (clustering :clarans data dist clusters max-neighbor)))

(defn denclue
  "DENsity CLUstering algorithm.

  Input:

  * data - sequence of samples
  * sigma - gaussian kernel parameter
  * m - number of selected samples, much smaller than number of all samples
  * tolerance (optional) - tolerance of hill-climbing procedure
  * min-pts (optional) - minimum number of neighbors for a core attractor
  
  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/DENCLUE.html)"
  ([data sigma m] (clustering :denclue data sigma m))
  ([data sigma m tolerance] (denclue data sigma m tolerance (max 2 (int (/ (count data) 200)))))
  ([data sigma m tolerance min-pts] (clustering :denclue data sigma m tolerance min-pts)))

(defn dbscan
  "Density-Based Spatial Clustering of Applications with Noise algorithm.

  Input:

  * data - sequence of samples
  * dist (optional) - distance method, default `euclidean`
  * min-pts - minimum number of neighbors
  * radius - the neighborhood radius

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/DBSCAN.html)"
  ([data min-pts radius] (dbscan data d/euclidean min-pts radius))
  ([data ^Distance dist min-pts ^double radius] (clustering :dbscan data dist (int min-pts) radius)))

(defn mec
  "Nonparametric Minimum Conditional Entropy Clustering algorithm.

  Input:

  * data - sequence of samples
  * dist (optional) - distance method, default `:euclidean`
  * max-clusters - maximum number of clusters
  * radius - the neighborhood radius

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/MEC.html)"
  ([data max-clusters radius] (mec data d/euclidean max-clusters radius))
  ([data dist max-clusters ^double radius] (clustering :mec data dist max-clusters radius)))

(defn spectral
  "Spectral clustering

  Input:

  * data - sequence of samples
  * clusters - number of clusters
  * samples (optional) - number of random samples for Nystrom approximation
  * sigma - width parameter for Gaussian kernel
  * max-iter (optional) - maximum number of iterations
  * tolerance (optional) - tolerance of k-means convergence test

  See more in [SMILE doc](https://haifengl.github.io/api/java/smile/clustering/SpectralClustering.html)"
  ([data clusters sigma] (clustering :spectral data clusters sigma))
  ([data clusters sigma max-iters tolerance] (clustering :spectral data clusters sigma max-iters tolerance))
  ([data ^long clusters ^long samples ^double sigma] (clustering :spectral data clusters samples sigma))
  ([data clusters samples sigma max-iters tolerance] (clustering :spectral data clusters samples sigma max-iters tolerance)))
