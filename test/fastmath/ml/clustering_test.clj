(ns fastmath.ml.clustering-test
  (:require [fastmath.ml.clustering :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.distance :as d]
            [clojure.edn :as edn]
            [clojure.java.io :as io]))

;; The Apache Commons Math clusterers are not tested here, only the wrapper:
;; options, results, the radius selection (Kneedle), `predict`, `wcss`, `silhouette`.
;; Reference values: `utils/fastmath/dev/clustering_ref_gen.py` (scikit-learn 1.9.0, kneed 0.9.0).

(def ^:private reference
  (delay (-> "clustering/clustering_reference.edn" io/resource slurp edn/read-string)))

(defn- dataset [name] (get-in @reference [:data (keyword name)]))

(defn- rng [seed] (r/rng :mersenne seed))

(defn- thrown-option
  "The `:option` of the `ex-info` thrown by `f`, `::no-exception` when `f` returns."
  [f]
  (try (f) ::no-exception
       (catch clojure.lang.ExceptionInfo e (:option (ex-data e)))))

(defn- cluster-partition
  "Set of sets of input indices of the clusters (noise map excluded)."
  [result]
  (set (map (comp set :indices) (remove :noise? result))))

(defn- label-partition
  "Set of sets of indices for `labels`, `skip?` labels are not included."
  ([labels] (label-partition labels (constantly false)))
  ([labels skip?]
   (->> (map-indexed vector labels)
        (remove (comp skip? second))
        (group-by second)
        vals
        (map #(set (map first %)))
        set)))

(defn- mean-point [points]
  (v/div (reduce v/add points) (count points)))

(defn- result-from-labels
  "A result map as the clustering functions return, built from labels (clusters ordered by label)."
  [x labels]
  (->> (map-indexed vector labels)
       (group-by second)
       (sort-by key)
       (mapv (fn [[_ ps]]
               (let [indices (mapv first ps)
                     data (mapv (comp vec x) indices)]
                 {:representative (vec (mean-point data))
                  :size (count indices)
                  :indices indices
                  :data data})))))

(def ^:private blobs (delay (dataset :blobs3-clean)))
(def ^:private blobs-3d (delay (dataset :blobs4-3d)))

(def ^:private tiny [[0 0] [0 1] [1 0] [10 10] [10 11] [11 10]])

;; ## kmeans++

(t/deftest kmeans++-options
  (let [x @blobs
        n (count x)
        run (fn [option value] (sut/kmeans++ x (assoc {:clusters 2} option value)))]
    (t/testing "invalid options throw ex-info with the option name"
      (doseq [[option value] [[:clusters 0] [:clusters -1] [:clusters (inc n)] [:clusters 1.5]
                              [:clusters nil] [:clusters "2"] [:clusters 2N] [:clusters Double/NaN]
                              [:trials 0] [:trials -1] [:trials nil] [:trials 1.0]
                              [:max-iters 1.5] [:max-iters nil] [:max-iters "1"]
                              [:empty-cluster-strategy :bogus] [:empty-cluster-strategy nil]
                              [:distance :cosine] [:distance 5] [:distance nil] [:distance "euclidean"]
                              [:rng 42] [:rng nil] [:rng "mt"]]]
        (t/is (= option (thrown-option #(run option value))) (str option " " (pr-str value)))))
    (t/testing "the ex-info carries :option, :value and :expected"
      (let [data (try (run :clusters 0) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
        (t/is (= [:clusters 0] ((juxt :option :value) data)))
        (t/is (string? (:expected data)))))
    (t/testing "boundary valid options"
      (doseq [[option value] [[:clusters 1] [:clusters n] [:clusters (int 2)] [:clusters (short 2)]
                              [:trials 1] [:trials 2] [:max-iters 0] [:max-iters -1] [:max-iters -5]
                              [:max-iters Long/MIN_VALUE] [:max-iters Long/MAX_VALUE]
                              [:empty-cluster-strategy :largest-variance]
                              [:empty-cluster-strategy :largest-points-number]
                              [:empty-cluster-strategy :farthest-point]
                              [:empty-cluster-strategy :error]
                              [:distance :euclidean] [:distance :chebyshev] [:distance :manhattan]
                              [:distance :canberra] [:distance :earth-movers]
                              [:distance d/euclidean] [:distance (fn [a b] (d/manhattan a b))]
                              [:rng (rng 1)] [:rng (r/rng :default)]]]
        (t/is (vector? (run option value)) (str option " " (pr-str value)))))
    (t/testing "unknown option keys are ignored"
      (t/is (= 2 (count (run :unknown 1)))))
    (t/testing "an empty option map equals no options"
      (t/is (= (sut/kmeans++ x {:rng (rng 1)}) (sut/kmeans++ x {:rng (rng 1)})))
      (t/is (= 1 (count (sut/kmeans++ x))))
      (t/is (= 1 (count (sut/kmeans++ x {})))))
    (t/testing "boolean options follow Clojure truthiness"
      (doseq [falsy [nil false]]
        (t/is (not (contains? (first (run :add-data? falsy)) :data))))
      (doseq [truthy [true 0 :yes "no"]]
        (t/is (contains? (first (run :add-data? truthy)) :data))))))

(t/deftest clustering-data-validation
  (doseq [[label f] {"kmeans++" #(sut/kmeans++ % {:clusters 1})
                     "fuzzy-kmeans" #(sut/fuzzy-kmeans % {:clusters 1})
                     "dbscan" #(sut/dbscan % {:eps 1.0})
                     "infer-dbscan-radius" #(sut/infer-dbscan-radius % :euclidean)}
          [what bad] {"nil" nil "empty" [] "empty list" '() "string" "abc" "map" {:a 1}
                      "ragged" [[1 2] [3]] "NaN" [[1 2] [Double/NaN 3]]
                      "Inf" [[1 2] [Double/POSITIVE_INFINITY 3]] "-Inf" [[1 2] [Double/NEGATIVE_INFINITY 3]]
                      "empty points" [[] []] "non-numeric" [["a" 1] [2 3]] "nil point" [[1 2] nil]
                      "number" 5}]
    (t/is (= :xss (thrown-option #(f bad))) (str label " " what))))

(t/deftest kmeans++-result
  (t/testing "well separated blobs are recovered as scikit-learn does"
    (doseq [{:keys [data clusters labels inertia]} (:kmeans @reference)
            :let [x (dataset data)
                  res (sut/kmeans++ x {:clusters clusters :trials 10 :rng (rng 1)})]]
      (t/is (= (label-partition labels) (cluster-partition res)) data)
      (t/testing "the representative is the mean of the cluster points"
        (doseq [{:keys [representative indices]} res]
          (t/is (v/delta-eq (mean-point (map x indices)) (vec representative) 1.0e-9))))
      (t/testing "within-cluster sum of squares equals the inertia"
        (t/is (m/delta-eq (sut/wcss res) inertia 1.0e-8)))))
  (t/testing "structure of the result"
    (let [x @blobs
          res (sut/kmeans++ x {:clusters 3 :rng (rng 1)})]
      (t/is (vector? res))
      (t/is (every? #(= #{:representative :size :indices :data} (set (keys %))) res))
      (t/is (= (range (count x)) (sort (mapcat :indices res))) "indices partition the input")
      (t/is (every? #(= (:size %) (count (:indices %)) (count (:data %))) res))
      (t/is (every? (fn [{:keys [indices data]}] (= (mapv (comp vec x) indices) data)) res)
            ":data holds the input points at :indices")
      (t/is (every? #(and (vector? (:representative %)) (vector? (:data %))) res))
      (t/is (= (count x) (reduce + (map :size res))))))
  (t/testing "one cluster is the mean of all points, k = n gives singletons"
    (let [x @blobs
          one (sut/kmeans++ x)]
      (t/is (= 1 (count one)))
      (t/is (v/delta-eq (mean-point (map vec x)) (vec (:representative (first one))) 1.0e-9))
      (let [all (sut/kmeans++ x {:clusters (count x) :rng (rng 1)})]
        (t/is (= (count x) (count all)))
        (t/is (every? #(= 1 (:size %)) all))
        (t/is (= (set (map vec x)) (set (map :representative all))))
        (t/is (m/delta-eq (sut/wcss all) 0.0 1.0e-20)))))
  (t/testing "one point"
    (t/is (= [{:representative [1.0 2.0] :size 1 :indices [0] :data [[1.0 2.0]]}]
             (sut/kmeans++ [[1 2]]))))
  (t/testing "options of the result"
    (t/is (= #{:representative :size :indices} (set (keys (first (sut/kmeans++ @blobs {:clusters 2 :add-data? false})))))))
  (t/testing "the same seed gives the same result, trials do not change a separable result"
    (let [f #(sut/kmeans++ @blobs {:clusters 3 :rng (rng 42)})]
      (t/is (= (f) (f))))
    (t/is (= (cluster-partition (sut/kmeans++ @blobs {:clusters 3 :rng (rng 1)}))
             (cluster-partition (sut/kmeans++ @blobs {:clusters 3 :trials 7 :rng (rng 2)})))))
  (t/testing "distances and strategies recover separable blobs"
    (let [expected (cluster-partition (sut/kmeans++ @blobs {:clusters 3 :trials 5 :rng (rng 1)}))]
      (doseq [distance [:euclidean :manhattan :chebyshev d/euclidean d/manhattan (fn [a b] (d/chebyshev a b))]]
        (t/is (= expected (cluster-partition (sut/kmeans++ @blobs {:clusters 3 :trials 5 :distance distance :rng (rng 1)})))
              (str distance)))
      (doseq [distance [:canberra :earth-movers]
              :let [res (sut/kmeans++ @blobs {:clusters 3 :distance distance :rng (rng 1)})]]
        (t/is (= (range 90) (sort (mapcat :indices res))) (str distance)))
      (doseq [strategy [:largest-variance :largest-points-number :farthest-point :error]]
        (t/is (= expected (cluster-partition (sut/kmeans++ @blobs {:clusters 3 :trials 5 :empty-cluster-strategy strategy :rng (rng 1)})))
              (str strategy)))))
  (t/testing "a function and its keyword give the same result"
    (t/is (= (sut/kmeans++ @blobs {:clusters 3 :distance :manhattan :rng (rng 5)})
             (sut/kmeans++ @blobs {:clusters 3 :distance d/manhattan :rng (rng 5)}))))
  (t/testing "max-iters"
    (doseq [iters [0 1 -1 -100]]
      (t/is (= (range 90) (sort (mapcat :indices (sut/kmeans++ @blobs {:clusters 3 :max-iters iters :rng (rng 1)}))))
            (str iters)))))

(t/deftest kmeans++-input-forms
  (let [expected (cluster-partition (sut/kmeans++ tiny {:clusters 2 :rng (rng 1)}))
        run #(cluster-partition (sut/kmeans++ % {:clusters 2 :rng (rng 1)}))]
    (t/is (= #{#{0 1 2} #{3 4 5}} expected))
    (t/testing "vectors, lists, lazy seqs, double arrays, vector types, longs"
      (t/is (= expected (run (map #(apply list %) tiny))))
      (t/is (= expected (run (map #(map double %) tiny))))
      (t/is (= expected (run (map double-array tiny))))
      (t/is (= expected (run (into-array (class (double-array 0)) (map double-array tiny)))))
      (t/is (= expected (run (map #(apply v/vec2 %) tiny))))
      (t/is (= expected (run (map #(apply v/vec2 (map long %)) tiny)))))
    (t/testing "numbers are 1-d points, results are doubles"
      (let [res (sut/kmeans++ [1 2 3 10 11 12] {:clusters 2 :rng (rng 1)})]
        (t/is (= #{#{0 1 2} #{3 4 5}} (cluster-partition res)))
        (t/is (= #{[2.0] [11.0]} (set (map :representative res))))
        (t/is (= #{[1.0] [2.0] [3.0] [10.0] [11.0] [12.0]} (set (mapcat :data res))))))
    (t/testing "3-d points"
      (t/is (= 2 (count (sut/kmeans++ [(v/vec3 1 2 3) (v/vec3 9 9 9)] {:clusters 2})))))
    (t/testing "input arrays are not changed"
      (let [a (double-array [1 2])
            b (double-array [9 9])]
        (sut/kmeans++ [a b] {:clusters 2})
        (t/is (= [1.0 2.0] (vec a)))
        (t/is (= [9.0 9.0] (vec b)))))
    (t/testing "a set is accepted, the indices follow its iteration order"
      (let [res (sut/kmeans++ #{[1 1] [2 2] [10 10]} {:clusters 2 :rng (rng 1)})]
        (t/is (= 3 (reduce + (map :size res))))))
    (t/testing "extreme coordinates"
      (let [res (sut/kmeans++ [[1e300 0] [1e300 1] [-1e300 0] [-1e300 1]] {:clusters 2 :rng (rng 1)})]
        (t/is (= #{#{0 1} #{2 3}} (cluster-partition res)))))))

;; ## fuzzy-kmeans

(t/deftest fuzzy-kmeans-options
  (let [x @blobs
        run (fn [option value] (sut/fuzzy-kmeans x (assoc {:clusters 2 :rng (rng 1)} option value)))]
    (t/testing "invalid options"
      (doseq [[option value] [[:clusters 0] [:clusters -1] [:clusters 91] [:clusters nil] [:clusters 1.5]
                              [:fuzziness 1] [:fuzziness 1.0] [:fuzziness 0.5] [:fuzziness 0] [:fuzziness -2]
                              [:fuzziness nil] [:fuzziness "2"] [:fuzziness Double/NaN]
                              [:fuzziness Double/POSITIVE_INFINITY]
                              [:epsilon 0] [:epsilon 0.0] [:epsilon -1] [:epsilon -1.0e-300] [:epsilon nil]
                              [:epsilon Double/NaN] [:epsilon Double/POSITIVE_INFINITY] [:epsilon "1"]
                              [:max-iters 1.5] [:max-iters nil]
                              [:distance :cosine] [:distance nil]
                              [:rng 7] [:rng nil]]]
        (t/is (= option (thrown-option #(run option value))) (str option " " (pr-str value)))))
    (t/testing "boundary valid options"
      (doseq [[option value] [[:fuzziness 1.0000001] [:fuzziness 100.0] [:fuzziness 3]
                              [:epsilon Double/MIN_VALUE] [:epsilon 1] [:epsilon 1.0e-3]
                              [:max-iters 0] [:max-iters -1] [:max-iters Long/MAX_VALUE]
                              [:clusters 1] [:clusters 90]
                              [:distance :manhattan] [:distance d/euclidean]]]
        (t/is (vector? (run option value)) (str option " " (pr-str value)))))
    (t/testing "boolean options follow Clojure truthiness"
      (t/is (not (contains? (first (run :add-data? nil)) :data)))
      (t/is (contains? (first (run :add-memberships? 1)) :memberships))
      (t/is (not (contains? (first (run :add-memberships? false)) :memberships))))))

(t/deftest fuzzy-kmeans-result
  (let [x @blobs
        n (count x)
        eps 1.0e-10]
    (t/testing "separable blobs and the structure of the result"
      (let [res (sut/fuzzy-kmeans x {:clusters 3 :rng (rng 3)})]
        (t/is (= (cluster-partition (sut/kmeans++ x {:clusters 3 :trials 5 :rng (rng 1)}))
                 (cluster-partition res)))
        (t/is (= (range n) (sort (mapcat :indices res))))
        (t/is (every? #(= #{:representative :size :indices :data} (set (keys %))) res))
        (t/is (map? (meta res)))
        (t/is (number? (:objective (meta res))))))
    (doseq [fuzziness [1.5 2.0 3.0]]
      (let [res (sut/fuzzy-kmeans x {:clusters 3 :fuzziness fuzziness :epsilon eps :add-memberships? true :rng (rng 3)})
            centres (mapv (comp vec :representative) res)
            u (mapv :memberships res)
            power (/ 2.0 (dec fuzziness))]
        (t/testing (str "memberships, fuzziness " fuzziness)
          (t/is (every? #(= n (count %)) u) "one membership per input point")
          (t/is (every? #(m/delta-eq (reduce + %) 1.0 1.0e-12) (apply map vector u)) "memberships of a point sum to 1")
          (t/is (every? (fn [{:keys [indices memberships]}]
                          (every? #(<= 0.0 (nth memberships %) 1.0) indices))
                        res))
          (t/testing "memberships follow from the centres"
            (doseq [i (range n) j (range 3)
                    :let [dij (d/euclidean (vec (nth x i)) (nth centres j))
                          expected (if (zero? dij)
                                     1.0
                                     (/ 1.0 (reduce + (map #(Math/pow (/ dij (d/euclidean (vec (nth x i)) (nth centres %))) power)
                                                           (range 3)))))]]
              (t/is (m/delta-eq expected (get-in u [j i]) 1.0e-12))))
          (t/testing "the hard assignment is the highest membership"
            (doseq [[j {:keys [indices]}] (map-indexed vector res)
                    i indices]
              (t/is (= j (apply max-key #(get-in u [% i]) (range 3))))))
          (t/testing "centres follow from the memberships (the fixed point, error of the order of epsilon)"
            (doseq [j (range 3)
                    :let [w (map #(Math/pow (get-in u [j %]) fuzziness) (range n))
                          c (v/div (reduce v/add (map #(v/mult (vec (nth x %1)) %2) (range n) w))
                                   (reduce + w))]]
              ;; observed error 3e-11 for epsilon 1e-10
              (t/is (v/delta-eq c (nth centres j) 1.0e-9))))
          (t/testing "the objective is the weighted sum of squared distances"
            (t/is (m/delta-eq (reduce + (for [i (range n) j (range 3)]
                                          (* (Math/pow (get-in u [j i]) fuzziness)
                                             (m/sq (d/euclidean (vec (nth x i)) (nth centres j))))))
                              (:objective (meta res))
                              1.0e-9))))))
    (t/testing "one cluster: membership 1, centre is the mean"
      (let [res (sut/fuzzy-kmeans x {:add-memberships? true})]
        (t/is (= (repeat n 1.0) (:memberships (first res))))
        (t/is (v/delta-eq (mean-point (map vec x)) (vec (:representative (first res))) 1.0e-9))))
    (t/testing "the objective is also in the metadata without memberships"
      (t/is (number? (:objective (meta (sut/fuzzy-kmeans x {:clusters 2}))))))
    (t/testing "the same seed gives the same result"
      (let [f #(sut/fuzzy-kmeans x {:clusters 3 :add-memberships? true :rng (rng 7)})]
        (t/is (= (f) (f)))))
    (t/testing "a distance function and its keyword"
      (t/is (= (sut/fuzzy-kmeans x {:clusters 2 :distance :manhattan :rng (rng 5)})
               (sut/fuzzy-kmeans x {:clusters 2 :distance d/manhattan :rng (rng 5)}))))))

;; ## dbscan

(t/deftest dbscan-options
  (let [x @blobs
        run (fn [option value] (sut/dbscan x (merge {:eps 1.0 :points 4} {option value})))]
    (t/testing "invalid options"
      (doseq [[option value] [[:eps 0] [:eps 0.0] [:eps -1.0] [:eps Double/NaN] [:eps Double/POSITIVE_INFINITY]
                              [:eps "1"] [:eps :a]
                              [:points 0] [:points -1] [:points 1.5] [:points "4"] [:points 2N]
                              [:eps-method :bogus] [:eps-method nil] [:eps-method "knee"]
                              [:distance :cosine] [:distance 1]]]
        (t/is (= option (thrown-option #(run option value))) (str option " " (pr-str value)))))
    (t/testing "boundary valid options"
      (doseq [[option value] [[:eps Double/MIN_VALUE] [:eps 1] [:eps 3/2] [:eps 1.5M] [:eps 1.0e300]
                              [:points 1] [:points (int 4)] [:points Long/MAX_VALUE]
                              [:eps-method :knee] [:eps-method :pairs]
                              [:distance :manhattan] [:distance d/euclidean]]]
        (t/is (vector? (run option value)) (str option " " (pr-str value)))))
    (t/testing "nil :eps and :points mean the default"
      (t/is (= (sut/dbscan x {:eps 1.0}) (sut/dbscan x {:eps 1.0 :points nil}))))
    (t/testing "boolean options follow Clojure truthiness"
      (t/is (not (contains? (first (run :add-data? nil)) :data)))
      (t/is (:noise? (last (run :noise? :yes))))
      (t/is (not-any? :noise? (run :noise? nil)))
      (t/is (not-any? :noise? (run :noise? false))))))

(t/deftest dbscan-reference
  (doseq [{:keys [data eps points distance labels core]} (:dbscan @reference)
          :let [x (dataset data)
                res (sut/dbscan x {:eps eps :points points :distance (keyword distance) :noise? true})
                clusters (butlast res)
                noise (last res)
                core? (set core)
                noise-labels (set (keep-indexed (fn [i l] (when (= l -1) i)) labels))
                label (str data " " eps " " points " " distance)]]
    ;; Border points can be assigned to a different cluster than in scikit-learn (it depends on
    ;; the order of expansion), so clusters are compared on core points, noise on all points.
    (t/is (= noise-labels (set (:indices noise))) (str "noise " label))
    (t/is (= (set (map #(set (filter core? %)) (label-partition labels #{-1})))
             (set (map #(set (filter core? (:indices %))) clusters)))
          (str "clusters " label))
    (t/is (= (count (disj (set labels) -1)) (count clusters)) (str "number of clusters " label))
    (t/testing "clusters and noise partition the input"
      (t/is (= (range (count x)) (sort (mapcat :indices res)))))
    (t/testing "the representative is the mean of the cluster"
      (doseq [{:keys [representative data]} clusters]
        (t/is (v/delta-eq (mean-point (map vec data)) (vec representative) 1.0e-9))))))

(t/deftest dbscan-semantics
  (t/testing "the neighbourhood is inclusive"
    (let [line [[0 0] [1 0] [2 0]]]
      (t/is (= [3] (map :size (sut/dbscan line {:eps 1.0 :points 2}))))
      (t/is (= [] (sut/dbscan line {:eps 0.999999 :points 2})))))
  (t/testing ":points does not count the point itself"
    (let [star [[0 0] [1 0] [-1 0] [0 1] [0 -1]]]
      (t/is (= [5] (map :size (sut/dbscan star {:eps 1.0 :points 4}))))
      (t/is (= [] (sut/dbscan star {:eps 1.0 :points 5})))))
  (t/testing "noise map"
    (let [x (conj (vec @blobs) [50 50] [-40 30])
          res (sut/dbscan x {:eps 1.0 :points 4 :noise? true})
          plain (sut/dbscan x {:eps 1.0 :points 4})]
      (t/is (= [30 30 30 2] (map :size res)))
      (t/is (= {:noise? true :size 2 :indices [90 91] :data [[50.0 50.0] [-40.0 30.0]]} (last res)))
      (t/is (= plain (butlast res)) "without :noise? the clusters are the same")
      (t/is (= 3 (count plain)))
      (t/is (= #{:noise? :size :indices} (set (keys (last (sut/dbscan x {:eps 1.0 :points 4 :noise? true :add-data? false}))))))
      (t/is (= (range 92) (sort (mapcat :indices res))))))
  (t/testing "no noise: the noise map is present and empty"
    (t/is (= {:noise? true :size 0 :indices []}
             (last (sut/dbscan @blobs {:eps 1.0 :points 4 :noise? true :add-data? false})))))
  (t/testing "all points are noise"
    (let [res (sut/dbscan @blobs {:eps 0.001 :points 4 :noise? true :add-data? false})]
      (t/is (= 1 (count res)))
      (t/is (= {:noise? true :size 90 :indices (vec (range 90))} (first res))))
    (t/is (= [] (sut/dbscan @blobs {:eps 0.001 :points 4})))
    (t/is (= [] (sut/dbscan @blobs {:eps 1.0 :points Long/MAX_VALUE}))))
  (t/testing "duplicates and one point"
    (t/is (= [10] (map :size (sut/dbscan (repeat 10 [1 1]) {:eps 0.5 :points 3}))))
    (t/is (= [1] (map :size (sut/dbscan [[1 1]] {:eps 1.0 :points 1 :noise? true}))))
    (t/is (= [] (sut/dbscan [[1 1]] {:eps 1.0 :points 1}))))
  (t/testing "1-d numbers"
    (t/is (= [3 3] (map :size (sut/dbscan [1 1.1 1.2 5 5.1 5.2] {:eps 0.5 :points 2})))))
  (t/testing "a distance function and its keyword"
    (t/is (= (sut/dbscan @blobs {:eps 1.2 :points 4 :distance :manhattan})
             (sut/dbscan @blobs {:eps 1.2 :points 4 :distance d/manhattan})))))

(t/deftest dbscan-default-radius
  (let [x @blobs]
    (t/testing "the default is the knee radius"
      (t/is (= (sut/dbscan x {:eps (sut/infer-dbscan-radius x :euclidean)})
               (sut/dbscan x))))
    (t/testing ":eps-method"
      (t/is (= (sut/dbscan x {:eps (sut/infer-dbscan-radius x :euclidean)})
               (sut/dbscan x {:eps-method :knee})))
      (t/is (vector? (sut/dbscan x {:eps-method :pairs}))))
    (t/testing "an explicit :eps wins"
      (t/is (= (sut/dbscan x {:eps 1.0 :points 4}) (sut/dbscan x {:eps 1.0 :points 4 :eps-method :pairs}))))
    (t/testing "the default :points is twice the dimension"
      (t/is (= (sut/dbscan x {:eps 1.0 :points 4}) (sut/dbscan x {:eps 1.0})))
      (t/is (= (sut/dbscan @blobs-3d {:eps 1.0 :points 6}) (sut/dbscan @blobs-3d {:eps 1.0}))))
    (t/testing "the radius needs data: fewer than 3 points, or more neighbours than points"
      (t/is (= :xss (thrown-option #(sut/dbscan [[1 1] [2 2]]))))
      (t/is (= :xss (thrown-option #(sut/dbscan [[1 1]]))))
      (t/is (= :points (thrown-option #(sut/dbscan [[0 0] [1 1] [2 2]])))))
    (t/testing "no knee throws, an explicit :eps works"
      (t/is (= :eps (thrown-option #(sut/dbscan (repeat 10 [1 1])))))
      (t/is (= [10] (map :size (sut/dbscan (repeat 10 [1 1]) {:eps 0.5 :points 3})))))))

;; ## infer-dbscan-radius

(t/deftest infer-dbscan-radius-knee
  (t/testing "the Kneedle radius equals kneed"
    ;; Cases whose knee index changes with S perturbed by 1e-9 are excluded (`:robust false`): the
    ;; k-distances of mutual nearest neighbours are exactly equal, which makes the threshold test a
    ;; tie decided by the last digit of the distances (scikit-learn and Apache Commons Math differ).
    (let [cases (filter :robust (:knee @reference))]
      (t/is (>= (count cases) 60))
      (doseq [{:keys [data points sensitivity distance eps]} cases
              :let [got (try (sut/infer-dbscan-radius (dataset data) (keyword distance)
                                                      {:points points :sensitivity sensitivity})
                             (catch clojure.lang.ExceptionInfo e (ex-data e)))]]
        (t/is (and (number? got) (m/delta-eq got eps 1.0e-12))
              (str data " k=" points " S=" sensitivity " " distance " kneed " eps " got " got)))))
  (t/testing "a keyword, a function and a DistanceMeasure give the same radius"
    (let [x @blobs
          expected (sut/infer-dbscan-radius x :euclidean)]
      (t/is (= expected (sut/infer-dbscan-radius x d/euclidean)))
      (t/is (= expected (sut/infer-dbscan-radius x (fn [a b] (d/euclidean a b)))))
      (t/is (= expected (sut/infer-dbscan-radius x (org.apache.commons.math3.ml.distance.EuclideanDistance.))))
      (t/is (= (sut/infer-dbscan-radius x :manhattan) (sut/infer-dbscan-radius x d/manhattan)))))
  (t/testing "the 2-arity uses :points twice the dimension, the knee method and sensitivity 1"
    (let [x @blobs]
      (t/is (= (sut/infer-dbscan-radius x :euclidean)
               (sut/infer-dbscan-radius x :euclidean {:method :knee :points 4 :sensitivity 1.0})))
      (t/is (= (sut/infer-dbscan-radius @blobs-3d :euclidean)
               (sut/infer-dbscan-radius @blobs-3d :euclidean {:points 6})))))
  (t/testing "the radius is a distance between points"
    (let [x (mapv vec @blobs)
          ds (set (for [a x b x] (d/euclidean a b)))]
      (t/is (some #(m/delta-eq % (sut/infer-dbscan-radius @blobs :euclidean) 1.0e-12) ds))))
  (t/testing "a knee on a clear elbow"
    ;; k-distances of 1-d points 0..20 step 1 and 5 outliers 100, 200, ...: the elbow is before the outliers
    (let [x (concat (range 21) [100 200 300 400 500])
          eps (sut/infer-dbscan-radius x :euclidean {:points 1})]
      (t/is (<= 1.0 eps 100.0)))))

(t/deftest infer-dbscan-radius-options
  (let [x @blobs]
    (t/testing "invalid options"
      (doseq [[option value] [[:method :bogus] [:method nil] [:method "knee"]
                              [:points 0] [:points -1] [:points 1.5] [:points 90] [:points 1000]
                              [:sensitivity -1] [:sensitivity -1.0e-12] [:sensitivity Double/NaN]
                              [:sensitivity Double/POSITIVE_INFINITY] [:sensitivity "1"] [:sensitivity nil]
                              [:rng 1] [:rng nil]]]
        (t/is (= option (thrown-option #(sut/infer-dbscan-radius x :euclidean {option value}))) (str option " " (pr-str value)))))
    (t/testing "boundary valid options"
      (doseq [[option value] [[:points 1] [:points 89] [:sensitivity 0] [:sensitivity 0.0] [:sensitivity 3]
                              [:method :pairs] [:rng (rng 1)]]]
        (t/is (number? (sut/infer-dbscan-radius x :euclidean {option value})) (str option " " (pr-str value)))))
    (t/testing "the distance"
      (doseq [distance [:cosine nil 5 "euclidean"]]
        (t/is (= :distance (thrown-option #(sut/infer-dbscan-radius x distance))) (pr-str distance))))
    (t/testing "the number of points: fewer than 3 throws, 3 is the minimum"
      (doseq [n [1 2]]
        (t/is (= :xss (thrown-option #(sut/infer-dbscan-radius (take n x) :euclidean {:points 1}))) (str n)))
      (t/is (= :eps (thrown-option #(sut/infer-dbscan-radius [[0 0] [1 0] [5 5]] :euclidean {:points 1})))
            "no knee on 3 points"))
    (t/testing "no knee throws ex-info asking for :eps"
      (doseq [[label x] {"identical points" (repeat 10 [1 1])
                         "evenly spaced points" (range 20)}]
        (t/is (= :eps (thrown-option #(sut/infer-dbscan-radius x :euclidean {:points 1}))) label)))
    (t/testing "the sensitivity changes the radius"
      (t/is (not= (sut/infer-dbscan-radius x :euclidean {:sensitivity 0})
                  (sut/infer-dbscan-radius x :euclidean {:sensitivity 5}))))))

(t/deftest infer-dbscan-radius-pairs
  (let [x @blobs
        ;; the definition: kmeans++ with two clusters on all pairwise distances, the lower centre
        by-definition (fn [seed]
                        (let [pts (mapv vec x)
                              ds (for [i (range (count pts)) j (range (inc i) (count pts))]
                                   (d/euclidean (pts i) (pts j)))]
                          (->> (sut/kmeans++ ds {:clusters 2 :add-data? false :rng (rng seed)})
                               (map (comp first :representative))
                               (reduce min))))]
    (doseq [seed [1 2 3]]
      ;; the distances differ in the last digit (Apache Commons Math vs fastmath.distance)
      (t/is (m/delta-eq (by-definition seed) (sut/infer-dbscan-radius x :euclidean {:method :pairs :rng (rng seed)}) 1.0e-9) (str seed)))
    (t/testing "a keyword and a function"
      (t/is (= (sut/infer-dbscan-radius x :euclidean {:method :pairs :rng (rng 1)})
               (sut/infer-dbscan-radius x d/euclidean {:method :pairs :rng (rng 1)}))))
    (t/testing "the radius is between the smallest and the largest distance"
      (let [r (sut/infer-dbscan-radius x :euclidean {:method :pairs :rng (rng 1)})]
        (t/is (< 0.0 r 20.0))))
    (t/testing ":points is not used"
      (t/is (= (sut/infer-dbscan-radius x :euclidean {:method :pairs :rng (rng 1)})
               (sut/infer-dbscan-radius x :euclidean {:method :pairs :points 1000 :rng (rng 1)}))))))

;; ## predict

(t/deftest predict-test
  (let [x @blobs
        res (sut/kmeans++ x {:clusters 3 :trials 5 :rng (rng 1)})
        reps (map (comp vec :representative) res)]
    (t/testing "the index of the nearest representative"
      (doseq [p (take 30 x)]
        (t/is (= (apply min-key #(d/euclidean (vec p) (nth reps %)) (range 3))
                 (sut/predict res p)))))
    (t/testing "every input point is predicted into its own cluster"
      (doseq [[j {:keys [indices]}] (map-indexed vector res)
              i indices]
        (t/is (= j (sut/predict res (nth x i))))))
    (t/testing "other distances"
      (doseq [[distance f] {:manhattan d/manhattan :chebyshev d/chebyshev}
              p [[0 0] [3 3] [6 0] [100 100]]]
        (t/is (= (apply min-key #(f (vec p) (nth reps %)) (range 3)) (sut/predict res p distance)))
        (t/is (= (sut/predict res p distance) (sut/predict res p f)))))
    (t/testing "ties go to the first cluster"
      (t/is (= 0 (sut/predict [{:representative [0 0]} {:representative [2 0]}] [1 0])))
      (t/is (= 0 (sut/predict [{:representative [2 0]} {:representative [0 0]}] [1 0]))))
    (t/testing "forms of the point"
      (t/is (= (sut/predict res [0 0]) (sut/predict res (double-array [0 0])) (sut/predict res (v/vec2 0 0)) (sut/predict res '(0 0))))
      (t/is (= 1 (sut/predict [{:representative [0]} {:representative [10]}] 9))))
    (t/testing "the noise map is ignored"
      (let [db (sut/dbscan x {:eps 1.0 :points 4 :noise? true})]
        (t/is (= (sut/predict (butlast db) [6 0]) (sut/predict db [6 0])))))
    (t/testing "fuzzy results"
      (let [f (sut/fuzzy-kmeans x {:clusters 3 :rng (rng 1)})]
        (t/is (every? #(<= 0 % 2) (map #(sut/predict f %) x)))))
    (t/testing "invalid arguments"
      (t/is (= :x (thrown-option #(sut/predict res [1]))))
      (t/is (= :x (thrown-option #(sut/predict res [1 2 3]))))
      (t/is (= :x (thrown-option #(sut/predict res [Double/NaN 1]))))
      (t/is (= :x (thrown-option #(sut/predict res nil))))
      (t/is (= :distance (thrown-option #(sut/predict res [1 1] :cosine))))
      (doseq [bad [nil [] [1 2] "a" 5 {:a 1}
                   [{:noise? true :size 2 :indices [0 1]}]]]
        (t/is (= :result (thrown-option #(sut/predict bad [1 1]))) (pr-str bad))))
    (t/testing "a map without :noise? is a cluster, :noise? false is a cluster"
      (t/is (= 0 (sut/predict [{:noise? false :representative [0 0]}] [1 1]))))))

;; ## wcss

(t/deftest wcss-test
  (t/testing "hand computed"
    (t/is (= 2.0 (sut/wcss [{:representative [1 0] :size 2 :data [[0 0] [2 0]]}
                            {:representative [5 5] :size 1 :data [[5 5]]}])))
    (t/is (= 8.0 (sut/wcss [{:representative [1 1] :size 2 :data [[0 0] [2 2]]}] :manhattan)))
    (t/is (= 2.0 (sut/wcss [{:representative [1 1] :size 2 :data [[0 0] [2 2]]}] :chebyshev)))
    (t/is (m/delta-eq 4.0 (sut/wcss [{:representative [1 1] :size 2 :data [[0 0] [2 2]]}]) 1.0e-14)))
  (t/testing "equals the sum of squared distances to the mean, scikit-learn labels"
    (doseq [{:keys [data labels distance wcss]} (:silhouette @reference)
            :let [res (result-from-labels (dataset data) labels)]]
      (t/is (m/delta-eq (sut/wcss res (keyword distance)) wcss 1.0e-9 1.0e-12) (str data " " distance))))
  (t/testing "k = 1 is the total sum of squares, it decreases with the number of clusters"
    (let [x @blobs
          total (let [mu (mean-point (map vec x))] (reduce + (map #(m/sq (d/euclidean (vec %) mu)) x)))
          values (map #(sut/wcss (sut/kmeans++ x {:clusters % :trials 10 :rng (rng 1)})) [1 2 3 4 5])]
      (t/is (m/delta-eq (first values) total 1.0e-9 1.0e-12))
      (t/is (apply > values))))
  (t/testing "the noise map and :noise? clusters"
    (let [db (sut/dbscan (conj (vec @blobs) [50 50]) {:eps 1.0 :points 4 :noise? true})]
      (t/is (= (sut/wcss (butlast db)) (sut/wcss db)))))
  (t/testing "invalid arguments"
    (let [res (sut/kmeans++ @blobs {:clusters 2 :add-data? false})]
      (t/is (= :result (thrown-option #(sut/wcss res)))))
    (doseq [bad [nil [] [1 2] [{:noise? true :size 1 :indices [0] :data [[1 1]]}]]]
      (t/is (= :result (thrown-option #(sut/wcss bad))) (pr-str bad)))
    (t/is (= :distance (thrown-option #(sut/wcss [{:representative [0] :data [[1]]}] :cosine))))))

;; ## silhouette

(t/deftest silhouette-reference
  (doseq [{:keys [data labels distance mean samples]} (:silhouette @reference)
          :let [res (result-from-labels (dataset data) labels)
                got (sut/silhouette res (keyword distance))]]
    (t/is (= #{:mean :samples :indices} (set (keys got))))
    (t/is (= (range (count labels)) (:indices got)))
    (t/is (m/delta-eq (:mean got) mean 1.0e-12) (str data " " distance " mean"))
    (t/is (= (count samples) (count (:samples got))))
    ;; observed maximal error 2e-15
    (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-12) (:samples got) samples)) (str data " " distance " samples"))
    (t/is (m/delta-eq (:mean got) (/ (reduce + (:samples got)) (count samples)) 1.0e-12))
    (t/is (every? #(<= -1.0 % 1.0) (:samples got)))))

(t/deftest silhouette-test
  (t/testing "hand computed: points 0, 1 and 5, clusters {0 1} and {5}"
    (let [res [{:representative [0.5] :data [[0] [1]] :indices [0 1]}
               {:representative [5] :data [[5]] :indices [2]}]
          got (sut/silhouette res)]
      ;; point 0: a = 1, b = 5 -> 0.8; point 1: a = 1, b = 4 -> 0.75; a singleton scores 0
      (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) [0.8 0.75 0.0] (:samples got))))
      (t/is (m/delta-eq (/ 1.55 3.0) (:mean got) 1.0e-15))
      (t/is (= [0 1 2] (:indices got)))))
  (t/testing "singletons and identical points score 0"
    (let [all (sut/kmeans++ (map vector (range 10)) {:clusters 10 :rng (rng 1)})]
      (t/is (= [0.0] (distinct (:samples (sut/silhouette all)))))
      (t/is (= 0.0 (:mean (sut/silhouette all)))))
    (let [same [{:representative [1] :data [[1] [1]] :indices [0 1]}
                {:representative [1] :data [[1] [1]] :indices [2 3]}]]
      (t/is (= {:mean 0.0 :samples [0.0 0.0 0.0 0.0] :indices [0 1 2 3]} (sut/silhouette same)))))
  (t/testing "samples are ordered by the position in the input, not by cluster"
    (let [res [{:representative [10] :data [[10] [11]] :indices [2 3]}
               {:representative [0.5] :data [[0] [1]] :indices [0 1]}]
          got (sut/silhouette res)]
      (t/is (= [0 1 2 3] (:indices got)))
      ;; point 0: a = 1, b = 10.5; point 1: a = 1, b = 9.5; symmetric for 10 and 11
      (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) (:samples got) [(/ 19.0 21.0) (/ 17.0 19.0) (/ 17.0 19.0) (/ 19.0 21.0)])))))
  (t/testing "the order of clusters does not change the result"
    (let [res (sut/kmeans++ @blobs {:clusters 3 :trials 5 :rng (rng 1)})]
      (t/is (= (sut/silhouette res) (sut/silhouette (reverse res))))
      (t/is (= (sut/silhouette res) (sut/silhouette (vec (shuffle res)))))))
  (t/testing "a distance function and its keyword"
    (let [res (sut/kmeans++ @blobs {:clusters 3 :rng (rng 1)})]
      (t/is (= (sut/silhouette res :manhattan) (sut/silhouette res d/manhattan)))
      (t/is (= (sut/silhouette res) (sut/silhouette res :euclidean)))))
  (t/testing "dbscan: noise points are not scored"
    (let [x (conj (vec @blobs) [50 50] [-40 30])
          db (sut/dbscan x {:eps 1.0 :points 4 :noise? true})
          got (sut/silhouette db)]
      (t/is (= (range 90) (:indices got)))
      (t/is (= got (sut/silhouette (butlast db))))
      (t/is (= (:samples got)
               (:samples (sut/silhouette (sut/dbscan @blobs {:eps 1.0 :points 4})))))))
  (t/testing "separable blobs have a high silhouette"
    (t/is (< 0.7 (:mean (sut/silhouette (sut/kmeans++ @blobs {:clusters 3 :trials 5 :rng (rng 1)}))))))
  (t/testing "invalid arguments"
    (let [one (sut/kmeans++ @blobs {:clusters 1})]
      (t/is (= :result (thrown-option #(sut/silhouette one)))))
    (t/is (= :result (thrown-option #(sut/silhouette (sut/kmeans++ @blobs {:clusters 2 :add-data? false})))))
    (doseq [bad [nil [] [1 2]
                 [{:noise? true :size 1 :indices [0] :data [[1]]}]
                 [{:representative [0] :data [[0]] :indices [0]}
                  {:noise? true :size 1 :indices [1] :data [[1]]}]]]
      (t/is (= :result (thrown-option #(sut/silhouette bad))) (pr-str bad)))
    (t/is (= :distance (thrown-option #(sut/silhouette [{:representative [0] :data [[0]] :indices [0]}
                                                        {:representative [1] :data [[1]] :indices [1]}]
                                                       :cosine))))))

;; ## Degenerate data

(t/deftest kmeans++-distinct-points
  ;; With the unlimited number of iterations Apache Commons Math never returns when there are fewer
  ;; distinct points than clusters (an empty cluster cannot be filled), so it is rejected up front.
  (doseq [[label x] {"identical points" (repeat 10 [1 1])
                     "two distinct points" [[1 1] [1 1] [2 2]]
                     "1-d" [1 1 1 2 2]
                     "negative zero equals zero" [[0.0] [-0.0] [1.0]]
                     "identical longs and doubles" [[1] [1.0] [2]]}
          :let [distinct-count (count (distinct (map #(mapv (fn [v] (+ 0.0 (double v))) (if (number? %) [%] %)) x)))]]
    (t/is (= :clusters (thrown-option #(sut/kmeans++ x {:clusters (inc distinct-count)}))) label)
    (t/is (= :clusters (thrown-option #(sut/kmeans++ x {:clusters (count x)}))) label)
    (t/is (= distinct-count (count (sut/kmeans++ x {:clusters distinct-count :rng (rng 1)}))) label)
    (t/is (= (count x) (reduce + (map :size (sut/kmeans++ x {:clusters distinct-count :rng (rng 1)})))) label))
  (t/testing "the same limit for every empty-cluster strategy and number of trials"
    (doseq [strategy [:largest-variance :largest-points-number :farthest-point :error]]
      (t/is (= :clusters (thrown-option #(sut/kmeans++ [[1 1] [1 1] [2 2]] {:clusters 3 :empty-cluster-strategy strategy :trials 3}))))))
  (t/testing "the error carries the value and the expected range"
    (let [data (try (sut/kmeans++ [[1 1] [1 1] [2 2]] {:clusters 3}) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
      (t/is (= 3 (:value data)))
      (t/is (re-find #"distinct" (:expected data)))))
  (t/testing "one cluster of identical points"
    (t/is (= [{:representative [1.0 1.0] :size 4 :indices [0 1 2 3] :data (vec (repeat 4 [1.0 1.0]))}]
             (sut/kmeans++ (repeat 4 [1 1]))))))

(t/deftest fuzzy-kmeans-identical-points
  ;; Documented: a cluster without points has the size 0 and the NaN representative.
  (let [res (sut/fuzzy-kmeans (repeat 5 [1 1]) {:clusters 2 :add-memberships? true :rng (rng 1)})
        [full empty-cluster] (sort-by (comp - :size) res)]
    (t/is (= 5 (:size full)))
    (t/is (= [1.0 1.0] (:representative full)))
    (t/is (= [1.0 1.0 1.0 1.0 1.0] (:memberships full)))
    (t/is (= {:size 0 :indices [] :data [] :memberships [0.0 0.0 0.0 0.0 0.0]}
             (dissoc empty-cluster :representative)))
    (t/is (every? #(Double/isNaN (double %)) (:representative empty-cluster)))
    (t/is (= 2 (count (:representative empty-cluster)))))
  (t/testing "empty clusters are not valid for the functions working on representatives"
    (let [res (sut/fuzzy-kmeans (repeat 5 [1 1]) {:clusters 2 :rng (rng 1)})]
      (t/is (= :result (thrown-option #(sut/predict res [1 1]))))
      (t/is (= :result (thrown-option #(sut/wcss res))))
      (t/is (= :result (thrown-option #(sut/silhouette res)))))))

;; ## Review fixes

(t/deftest ex-data-has-no-data
  (let [thrown (fn [f] (try (f) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))
        big (conj (vec (repeat 200000 [1 2])) [1])]
    (t/testing "a point of another dimension: the first different point, not the data"
      (let [data (thrown #(sut/kmeans++ big))]
        (t/is (= :xss (:option data)))
        (t/is (= [1.0] (:value data)))
        (t/is (< (count (pr-str data)) 300))))
    (t/testing "a collection that is not a sequence of points: no value"
      (t/is (nil? (:value (thrown #(sut/kmeans++ {:a 1}))))))
    (t/testing "a small invalid value is returned"
      (t/is (= "abc" (:value (thrown #(sut/kmeans++ "abc")))))
      (t/is (= ["a" 1] (:value (thrown #(sut/kmeans++ [[1 2] ["a" 1]]))))))
    (t/testing "a whole result is not returned"
      (let [res (sut/kmeans++ @blobs {:clusters 2 :add-data? false :rng (rng 1)})]
        (doseq [f [#(sut/silhouette res) #(sut/wcss res) #(sut/silhouette (take 1 res))
                   #(sut/predict [1 2] [1 1])]]
          (t/is (nil? (:value (thrown f)))))))))

(t/deftest pairs-radius-equal-distances
  ;; all distances between points are equal: the radius is not defined, the error names :eps
  (doseq [x [(repeat 6 [1 1]) [[0 0] [1 1]] [[0] [0] [0]]]
          :when (>= (count x) 3)]
    (t/is (= :eps (thrown-option #(sut/infer-dbscan-radius x :euclidean {:method :pairs}))))
    (t/is (= :eps (thrown-option #(sut/dbscan x {:eps-method :pairs})))))
  (t/testing "two different distances are enough"
    (t/is (number? (sut/infer-dbscan-radius [[0 0] [1 0] [3 0]] :euclidean {:method :pairs :rng (rng 1)})))))

(t/deftest silhouette-empty-cluster
  (let [res (conj (sut/kmeans++ @blobs {:clusters 3 :rng (rng 1)})
                  {:representative [0 0] :size 0 :indices [] :data []})]
    (t/is (= :result (thrown-option #(sut/silhouette res))))
    (t/is (map? (sut/silhouette (butlast res))) "without the empty cluster it works")))

(t/deftest k-distances-selection
  ;; the k-th neighbour distance equals the sorted row element, also with duplicate points and ties
  (let [naive (fn [x k dist]
                (vec (sort (for [a x] (nth (sort (for [b x] (dist a b))) k)))))
        k-distances #'sut/k-distances
        ->arrays #'sut/points->arrays
        ->points #'sut/->points
        ->distance #'sut/->distance
        r (rng 11)
        random-data (fn [n dim] (vec (repeatedly n #(vec (repeatedly dim (fn [] (r/drandom r)))))))
        with-duplicates (vec (concat (repeat 5 [1.0 1.0]) (repeat 5 [2.0 2.0]) [[3.0 3.0] [3.0 4.0]]))
        grid (vec (for [i (range 6) j (range 6)] [(double i) (double j)]))]
    (doseq [[label x] {"random 2d" (random-data 60 2) "random 5d" (random-data 40 5) "duplicates" with-duplicates
                       "grid (ties)" grid "1-d" (random-data 30 1)}
            k [1 2 4 8 (dec (count x))]
            :when (<= k (dec (count x)))
            [dname dist] {:euclidean d/euclidean :manhattan d/manhattan}]
      (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-12)
                               (naive x k dist)
                               (k-distances (->arrays (->points x)) (->distance dname) k)))
            (str label " k=" k " " dname)))))
