(ns fastmath.ml.clustering.hierarchical
  "Hierarchical clustering algorithms behind [[fastmath.ml.clustering/hclust]].

  Agglomerative clustering (a port of R `stats::hclust`, and the minimum spanning tree and nearest-neighbour chain algorithms as in SciPy), divisive clustering (a port of DIANA from R `cluster`), tree cutting, cophenetic distances and the conversion to the SciPy linkage matrix.

  Distances are stored in a condensed array: the pairs `i < j` (0-based) of `n` observations in row order. A tree is a map with `:merge` (rows of the R encoding: `-i` is the observation `i` counted from 1, `j` is the cluster made in the step `j`), `:height` and `:order` (0-based)."
  {:no-doc true}
  (:require [fastmath.core :as m])
  (:import [fastmath.java Array]
           [clojure.lang IFn$LLD IFn$OOD]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Condensed distances

(defn condensed-size
  "Number of pairs of `n` observations."
  ^long [^long n]
  (m/long-quot (m/long-mult n (m/long-dec n)) 2))

(defn condensed-index
  "Position of the pair `i < j` in the condensed array of `n` observations."
  ^long [^long n ^long i ^long j]
  (m/long-add (m/long-sub (m/long-mult i n) (m/long-quot (m/long-mult i (m/long-add i 3)) 2))
              (m/long-dec j)))

(defn distance-at
  "Distance between observations `i` and `j` stored in the condensed array `d`, `0.0` for `i = j`."
  ^double [^doubles d ^long n ^long i ^long j]
  (cond
    (m/< i j) (Array/aget d (condensed-index n i j))
    (m/> i j) (Array/aget d (condensed-index n j i))
    :else 0.0))

(defn point-distance
  "Distance between points `a` and `b` for `dist`, a function of two points (a primitive call for a function returning a double). `NaN` when `dist` does not return a number."
  ^double [dist a b]
  (if (instance? IFn$OOD dist)
    (.invokePrim ^IFn$OOD dist a b)
    (let [x (dist a b)]
      (if (number? x) (double x) ##NaN))))

(defn data-distances
  "Condensed distances between `points` (a vector of double arrays) for `dist`, a function of two points."
  ^doubles [points dist]
  (let [n (count points)
        d (double-array (condensed-size n))]
    (loop [i (long 0) idx (long 0)]
      (if (m/< i n)
        (let [p (nth points i)
              next-idx (long (loop [j (m/long-inc i) idx idx]
                               (if (m/< j n)
                                 (do (Array/aset d idx (point-distance dist p (nth points j)))
                                     (recur (m/long-inc j) (m/long-inc idx)))
                                 idx)))]
          (recur (m/long-inc i) next-idx))
        d))))

(defn valid-distances?
  "`true` when all condensed distances are finite and non-negative."
  [^doubles d]
  (let [size (alength d)]
    (loop [i (long 0)]
      (cond
        (m/>= i size) true
        (let [x (Array/aget d i)] (and (m/valid-double? x) (m/>= x 0.0))) (recur (m/long-inc i))
        :else false))))

(defn constant?
  "`true` when all values of `xs` are equal (`false` for `NaN`)."
  [^doubles xs]
  (let [size (alength xs)
        x0 (if (m/pos? size) (Array/aget xs 0) 0.0)]
    (loop [i (long 1)]
      (cond
        (m/>= i size) true
        (m/== x0 (Array/aget xs i)) (recur (m/long-inc i))
        :else false))))

(defn matrix-distances
  "Condensed distances from the lower triangle of a square matrix (a vector of double arrays)."
  ^doubles [rows]
  (let [n (count rows)
        d (double-array (condensed-size n))]
    (dotimes [j n]
      (let [^doubles row (nth rows j)]
        (dotimes [i j]
          (Array/aset d (condensed-index n i j) (Array/aget row i)))))
    d))

;; ## Tree encoding

(defn- canonical-row
  "A merge row in the R order: a singleton before a cluster, two singletons by observation, two clusters by step."
  [^long a ^long b]
  (cond
    (and (m/neg? a) (m/neg? b)) (if (m/> a b) [a b] [b a])
    (m/neg? a) [a b]
    (m/neg? b) [b a]
    :else [(m/long-min a b) (m/long-max a b)]))

(defn subtree-leaves
  "Observations (0-based) under `node` of `merge` (`-i` an observation, `j` a step), left branch first."
  [merge ^long node]
  (loop [stack [node]
         leaves (transient [])]
    (if-let [x (peek stack)]
      (let [x (long x)
            stack (pop stack)]
        (if (m/neg? x)
          (recur stack (conj! leaves (m/long-sub -1 x)))
          (let [[a b] (merge (m/long-dec x))]
            (recur (conj stack b a) leaves))))
      (persistent! leaves))))

(defn- root-node
  ^long [merge]
  (if (seq merge) (count merge) -1))

(defn leaf-order
  "Order of the leaves of a dendrogram (0-based), as R `hclust` `order`."
  [merge]
  (subtree-leaves merge (root-node merge)))

(defn- find-root
  "Root of `x` in the union-find `parent` (with path halving)."
  ^long [^longs parent ^long x]
  (loop [x x]
    (let [p (Array/get parent x)]
      (if (m/== p x)
        x
        (let [gp (Array/get parent p)]
          (Array/set parent x gp)
          (recur gp))))))

(deftype ^:private Merges [^longs as ^longs bs ^doubles heights])

(defn- merges->tree
  "Tree from merges of clusters given by any of their observations (`as`, `bs`, 0-based) at `heights`, in the order of `steps` (indices of the merges)."
  [^long n ^Merges merges steps]
  (let [^longs as (.as merges)
        ^longs bs (.bs merges)
        ^doubles heights (.heights merges)
        parent (long-array (range n))
        label (long-array (map (fn [^long i] (m/long-sub -1 i)) (range n)))
        merge (loop [steps steps
                     step (long 1)
                     merge (transient [])]
                (if-let [t (first steps)]
                  (let [t (long t)
                        ra (find-root parent (Array/get as t))
                        rb (find-root parent (Array/get bs t))
                        row (canonical-row (Array/get label ra) (Array/get label rb))]
                    (Array/set parent rb ra)
                    (Array/set label ra step)
                    (recur (rest steps) (m/long-inc step) (conj! merge row)))
                  (persistent! merge)))]
    {:merge merge
     :height (mapv (fn [^long t] (Array/aget heights t)) steps)
     :order (leaf-order merge)}))

;; ## Agglomerative clustering

(def ^:private squared-methods #{:ward :centroid :median})

(defmacro ^:private lance-williams
  "Lance-Williams update of the distance from the cluster `k` to the union of clusters `i` and `j` (R `hclust` formulas)."
  [method dik djk dij mi mj mk]
  `(let [dik# ~dik djk# ~djk dij# ~dij mi# ~mi mj# ~mj mk# ~mk]
     (case ~method
       :single (m/min dik# djk#)
       :complete (m/max dik# djk#)
       :average (m// (m/+ (m/* mi# dik#) (m/* mj# djk#)) (m/+ mi# mj#))
       :weighted (m/* 0.5 (m/+ dik# djk#))
       (:ward :ward-d) (m// (m/- (m/+ (m/* (m/+ mi# mk#) dik#) (m/* (m/+ mj# mk#) djk#)) (m/* mk# dij#))
                           (m/+ mi# mj# mk#))
       :centroid (m// (m/- (m/+ (m/* mi# dik#) (m/* mj# djk#)) (m// (m/* mi# mj# dij#) (m/+ mi# mj#)))
                      (m/+ mi# mj#))
       :median (m/* 0.5 (m/- (m/+ dik# djk#) (m/* 0.5 dij#))))))

(defn- square!
  [^doubles d]
  (dotimes [i (alength d)]
    (let [x (Array/aget d i)] (Array/aset d i (m/* x x)))))

(defn- height
  "Merge height from a working distance: the root for methods working on squared distances."
  ^double [method ^double x]
  (if (squared-methods method) (m/sqrt x) x))

(defn- generic-merges
  "R `hclust` (Murtagh's nearest-neighbour list, ties to the lowest pair). Consumes `d`."
  ^Merges [^doubles d ^long n method]
  (let [steps (m/long-dec n)
        active (boolean-array n true)
        size (double-array n 1.0)
        nn (long-array n -1)
        nn-dist (double-array n ##Inf)
        as (long-array steps)
        bs (long-array steps)
        heights (double-array steps)
        nearest! (fn [^long i]
                   (loop [j (m/long-inc i) best ##Inf bj (long -1)]
                     (if (m/< j n)
                       (if (aget active j)
                         (let [x (Array/aget d (condensed-index n i j))]
                           (if (m/< x best) (recur (m/long-inc j) x j) (recur (m/long-inc j) best bj)))
                         (recur (m/long-inc j) best bj))
                       (do (Array/set nn i bj)
                           (Array/aset nn-dist i best)))))]
    (when (squared-methods method) (square! d))
    (dotimes [i steps] (nearest! i))
    (dotimes [step steps]
      (let [i2 (long (loop [i (long 0) best ##Inf bi (long 0)]
                       (if (m/< i steps)
                         (if (and (aget active i) (m/< (Array/aget nn-dist i) best))
                           (recur (m/long-inc i) (Array/aget nn-dist i) i)
                           (recur (m/long-inc i) best bi))
                         bi)))
            j2 (Array/get nn i2)
            dij (Array/aget d (condensed-index n i2 j2))
            mi (Array/aget size i2)
            mj (Array/aget size j2)]
        (Array/set as step i2)
        (Array/set bs step j2)
        (Array/aset heights step (height method dij))
        (aset active j2 false)
        (loop [k (long 0) best ##Inf bk (long -1)]
          (if (m/< k n)
            (if (and (aget active k) (m/not== k i2))
              (let [ik (if (m/< i2 k) (condensed-index n i2 k) (condensed-index n k i2))
                    jk (if (m/< j2 k) (condensed-index n j2 k) (condensed-index n k j2))
                    v (lance-williams method (Array/aget d ik) (Array/aget d jk) dij mi mj (Array/aget size k))]
                (Array/aset d ik v)
                (if (m/< i2 k)
                  (if (m/< v best) (recur (m/long-inc k) v k) (recur (m/long-inc k) best bk))
                  (do (when (m/< v (Array/aget nn-dist k))
                        (Array/aset nn-dist k v)
                        (Array/set nn k i2))
                      (recur (m/long-inc k) best bk))))
              (recur (m/long-inc k) best bk))
            (do (Array/aset size i2 (m/+ mi mj))
                (Array/aset nn-dist i2 best)
                (Array/set nn i2 bk))))
        (dotimes [i steps]
          (when (and (aget active i)
                     (let [x (Array/get nn i)] (or (m/== x i2) (m/== x j2))))
            (nearest! i)))))
    (Merges. as bs heights)))

(defn- mst-merges
  "Single linkage by Prim's minimum spanning tree (SciPy `mst_single_linkage`)."
  ^Merges [^doubles d ^long n]
  (let [steps (m/long-dec n)
        merged (boolean-array n)
        dmin (double-array n ##Inf)
        as (long-array steps)
        bs (long-array steps)
        heights (double-array steps)]
    (loop [k (long 0) x (long 0)]
      (when (m/< k steps)
        (aset merged x true)
        (let [y (long (loop [i (long 0) best ##Inf y (long 0)]
                        (if (m/< i n)
                          (if (aget merged i)
                            (recur (m/long-inc i) best y)
                            (let [di (m/min (Array/aget dmin i) (distance-at d n x i))]
                              (Array/aset dmin i di)
                              (if (m/< di best) (recur (m/long-inc i) di i) (recur (m/long-inc i) best y))))
                          y)))]
          (Array/set as k x)
          (Array/set bs k y)
          (Array/aset heights k (Array/aget dmin y))
          (recur (m/long-inc k) y))))
    (Merges. as bs heights)))

(defn- nearest-in-chain
  "The nearest active cluster to `x` and its distance (`[y distance]`), preferring `prev` (the previous element of the chain, `-1` for none) on ties."
  [^doubles d ^doubles size ^long x ^long prev]
  (let [n (alength size)]
    (loop [i (long 0)
           y prev
           cur (if (m/neg? prev) ##Inf (distance-at d n x prev))]
      (if (m/< i n)
        (if (or (m/zero? (Array/aget size i)) (m/== i x))
          (recur (m/long-inc i) y cur)
          (let [di (distance-at d n x i)]
            (if (m/< di cur) (recur (m/long-inc i) i di) (recur (m/long-inc i) y cur))))
        [y cur]))))

(defn- nn-chain-merges
  "Nearest-neighbour chain for reducible linkages (SciPy `nn_chain`). Consumes `d`."
  ^Merges [^doubles d ^long n method]
  (let [steps (m/long-dec n)
        size (double-array n 1.0)
        chain (long-array n)
        as (long-array steps)
        bs (long-array steps)
        heights (double-array steps)]
    (when (squared-methods method) (square! d))
    (loop [k (long 0) len (long 0)]
      (when (m/< k steps)
        (when (m/zero? len)
          (Array/set chain 0 (long (first (filter #(m/pos? (Array/aget size (long %))) (range n))))))
        (let [[x y dxy len] (loop [len (m/long-max len 1)]
                              (let [x (Array/get chain (m/long-dec len))
                                    prev (if (m/> len 1) (Array/get chain (m/long-sub len 2)) -1)
                                    [y cur] (nearest-in-chain d size x prev)
                                    y (long y)]
                                (if (and (m/> len 1) (m/== y prev))
                                  [x y cur (m/long-sub len 2)]
                                  (do (Array/set chain len y)
                                      (recur (m/long-inc len))))))
              x (long x)
              y (long y)
              dxy (double dxy)
              lo (m/long-min x y)
              hi (m/long-max x y)
              nx (Array/aget size lo)
              ny (Array/aget size hi)]
          (Array/set as k lo)
          (Array/set bs k hi)
          (Array/aset heights k (height method dxy))
          (Array/aset size lo 0.0)
          (Array/aset size hi (m/+ nx ny))
          (dotimes [i n]
            (let [ni (Array/aget size i)]
              (when (and (m/pos? ni) (m/not== i hi))
                (let [ix (if (m/< i lo) (condensed-index n i lo) (condensed-index n lo i))
                      iy (if (m/< i hi) (condensed-index n i hi) (condensed-index n hi i))]
                  (Array/aset d iy (lance-williams method (Array/aget d ix) (Array/aget d iy) dxy nx ny ni))))))
          (recur (m/long-inc k) (long len)))))
    (Merges. as bs heights)))

(defn agglomerate
  "Agglomerative tree of `n` observations with condensed distances `d` (consumed) for a linkage `method` and an `algorithm` (`:generic` or `:fast`)."
  [^doubles d ^long n method algorithm]
  (if (m/< n 2)
    {:merge [] :height [] :order (vec (range n))}
    (let [fast? (and (= algorithm :fast) (not (#{:centroid :median} method)))
          ^Merges merges (cond
                           (not fast?) (generic-merges d n method)
                           (= method :single) (mst-merges d n)
                           :else (nn-chain-merges d n method))
          ^doubles heights (.heights merges)
          steps (range (m/long-dec n))
          steps (if fast? (sort-by (fn [^long t] (Array/aget heights t)) steps) steps)]
      (merges->tree n merges steps))))

;; ## Divisive clustering (DIANA)

;; The state of R `cluster` `splyt` (1-based arrays of length `n + 1`): `kwan[j]` the size of the
;; cluster starting at the position `j` (0 inside a cluster), `ner` the observations (1-based) by
;; position, `ban[j]` the banner height between the positions `j - 1` and `j`, `cs` the diameter.
(deftype ^:private DianaState [^doubles d ^long n ^longs kwan ^longs ner ^doubles ban ^double cs])

(defn- position-distance
  "Distance between the observations at the positions `a` and `b`."
  ^double [^DianaState s ^long a ^long b]
  (let [^longs ner (.ner s)]
    (distance-at (.d s) (.n s) (m/long-dec (Array/get ner a)) (m/long-dec (Array/get ner b)))))

(defn- diameter
  "The largest distance between the observations at the positions `ja..jb`."
  ^double [^DianaState s ^long ja ^long jb]
  (loop [k ja dm 0.0]
    (if (m/< k jb)
      (recur (m/long-inc k)
             (double (loop [j (m/long-inc k) dm dm]
                       (if (m/<= j jb)
                         (recur (m/long-inc j) (m/max dm (position-distance s k j)))
                         dm))))
      dm)))

(defn- distance-sum
  "Sum of distances from the observation at the position `l` to the positions `from..to`."
  ^double [^DianaState s ^long l ^long from ^long to]
  (loop [j from acc 0.0]
    (if (m/<= j to)
      (recur (m/long-inc j) (m/+ acc (position-distance s l j)))
      acc)))

(defn- move-to!
  "Moves `ner[from]` to the position `to` (`from <= to`), shifting the elements between to the left."
  [^longs ner ^long from ^long to]
  (when (m/not== from to)
    (let [x (Array/get ner from)]
      (loop [l from]
        (when (m/< l to)
          (Array/set ner l (Array/get ner (m/long-inc l)))
          (recur (m/long-inc l))))
      (Array/set ner to x))))

(defn- most-distant
  "Position in `ja..jb` with the largest sum of distances to the cluster (the first on ties)."
  ^long [^DianaState s ^long ja ^long jb]
  (loop [l ja best -1.0 bl (long -1)]
    (if (m/<= l jb)
      (let [sd (distance-sum s l ja jb)]
        (if (m/< best sd) (recur (m/long-inc l) sd l) (recur (m/long-inc l) best bl)))
      bl)))

(defn- move-splinters!
  "Moves observations to the splinter group at the end of `ja..jb` while they are closer to it; returns the last position of the remaining group."
  ^long [^DianaState s ^long ja ^long jb]
  (let [^longs kwan (.kwan s)
        ^longs ner (.ner s)]
    (loop [splyn (long 1) jma (m/long-dec jb)]
      (let [rest (m/long-sub jma ja)
            [bdyff jaway] (loop [l ja bdyff -1.0 jaway (long -1)]
                            (if (m/<= l jma)
                              (let [da (m// (distance-sum s l ja jma) (double rest))
                                    db (m// (distance-sum s l (m/long-inc jma) jb) (double splyn))
                                    dyff (m/- da db)]
                                (if (m/< bdyff dyff) (recur (m/long-inc l) dyff l) (recur (m/long-inc l) bdyff jaway)))
                              [bdyff jaway]))]
        (if (m/<= (double bdyff) 0.0)
          jma
          (let [jmb (m/long-inc jma)]
            (move-to! ner (long jaway) jma)
            (loop [lxx jmb]
              (when (m/<= lxx jb)
                (let [l1 (m/long-dec lxx)
                      a (Array/get ner l1)
                      b (Array/get ner lxx)]
                  (when-not (m/< a b)
                    (Array/set ner l1 b)
                    (Array/set ner lxx a)
                    (recur (m/long-inc lxx))))))
            (Array/add kwan ja -1)
            (Array/set kwan jma (m/long-inc (Array/get kwan jmb)))
            (Array/set kwan jmb 0)
            (let [jma (m/long-dec jma)]
              (if (m/== jma ja) jma (recur (m/long-inc splyn) jma)))))))))

(defn- lowest-first!
  "Swaps the two parts `ja..jma` and `jma+1..jb` when the second starts with a lower observation; returns the start of the second part."
  ^long [^DianaState s ^long ja ^long jma ^long jb]
  (let [^longs kwan (.kwan s)
        ^longs ner (.ner s)
        jmb (m/long-inc jma)]
    (if (m/>= (Array/get ner ja) (Array/get ner jmb))
      (do (loop [lgrb jmb lxxa ja]
            (when (m/<= lgrb jb)
              (let [lxxa (m/long-inc lxxa)
                    lchan (Array/get ner lgrb)]
                (loop [lxf lgrb]
                  (when (m/>= lxf lxxa)
                    (Array/set ner lxf (Array/get ner (m/long-dec lxf)))
                    (recur (m/long-dec lxf))))
                (Array/set ner (m/long-dec lxxa) lchan)
                (recur (m/long-inc lgrb) lxxa))))
          (let [llq (Array/get kwan jmb)
                _ (Array/set kwan jmb 0)
                jmb (m/long-sub (m/long-add ja jb) jma)]
            (Array/set kwan jmb (Array/get kwan ja))
            (Array/set kwan ja llq)
            jmb))
      jmb)))

(defn- split!
  "Splits the cluster starting at the position `ja` (one step of R `cluster` `splyt`), `first?` for the whole data."
  [^DianaState s ^long ja first?]
  (let [^longs kwan (.kwan s)
        ^longs ner (.ner s)
        ^doubles ban (.ban s)
        jb (m/long-dec (m/long-add ja (Array/get kwan ja)))]
    (if (m/== (Array/get kwan ja) 2)
      (do (Array/set kwan ja 1)
          (Array/set kwan jb 1)
          (Array/aset ban jb (position-distance s ja jb)))
      (let [lndsd (most-distant s ja jb)]
        (Array/add kwan ja -1)
        (Array/set kwan jb 1)
        (move-to! ner lndsd jb)
        (let [jmb (lowest-first! s ja (move-splinters! s ja jb) jb)]
          (Array/aset ban jmb (if first? (.cs s) (diameter s ja jb))))))))

(defn- next-split
  "Start of the next cluster to split after splitting the one at `ja` ending at `jb` (R `splyt` labels 420)."
  ^long [^longs kwan ^long ja ^long jb]
  (let [n (m/long-dec (alength kwan))
        scan (fn ^long [^long ja]
               (loop [ja (m/long-add ja (Array/get kwan ja))]
                 (cond
                   (m/> ja n) -1
                   (m/<= (Array/get kwan ja) 1) (recur (m/long-add ja (Array/get kwan ja)))
                   :else ja)))
        found (if (m/not== jb n) (.invokePrim ^clojure.lang.IFn$LL scan ja) -1)]
    (cond
      (m/pos? found) found
      (m/== (Array/get kwan 1) 1) (.invokePrim ^clojure.lang.IFn$LL scan 1)
      :else 1)))

(defn- banner-coefficient
  "Divisive coefficient from the banner `ban` (1-based, R `cluster` `bncoef`)."
  ^double [^doubles ban ^long n]
  (let [b (fn ^double [^long k] (Array/aget ban (m/long-inc k)))
        sup (reduce m/max 0.0 (map b (range 1 n)))]
    (m// (double (reduce m/+ 0.0 (map (fn [^long k]
                                (m/- 1.0 (m// (m/min (Array/aget ban (m/long-inc (m/long-max k 1)))
                                                     (Array/aget ban (m/long-inc (m/long-min (m/long-inc k) (m/long-dec n)))))
                                              (double sup))))
                              (range n))))
         (double n))))

(defn- banner-merges
  "Merges `[a b height]` of the split tree, by the ascending banner height (R `splyt` label 500)."
  [^DianaState s]
  (let [n (.n s)
        ^doubles ban (.ban s)
        ^longs ner (.ner s)
        step-of (fn [merge ^long l]
                  (reduce (fn [^long l [^long j [a b]]] (if (or (m/== l (long a)) (m/== l (long b))) (m/long-inc j) l))
                          l (map-indexed vector merge)))]
    (loop [used #{}
           merge []]
      (if (m/< (count merge) (m/long-dec n))
        (let [nj (long (reduce (fn [^long best ^long j]
                                 (if (and (not (used j))
                                          (m/>= (if (m/neg? best) (.cs s) (Array/aget ban best)) (Array/aget ban j)))
                                   j best))
                               -1 (range 2 (m/long-inc n))))
              l1 (step-of merge (m/long-sub 0 (Array/get ner (m/long-dec nj))))
              l2 (step-of merge (m/long-sub 0 (Array/get ner nj)))]
          (recur (conj used nj) (conj merge [l1 l2 (Array/aget ban nj)])))
        merge))))

(defn diana
  "Divisive tree of `n` observations with condensed distances `d` (R `cluster` `diana`, mapped as `as.hclust`). Adds `:dc`, the divisive coefficient (`NaN` for less than 2 observations)."
  [^doubles d ^long n]
  (if (m/< n 2)
    {:merge [] :height [] :order (vec (range n)) :dc ##NaN}
    (let [kwan (long-array (m/long-inc n))
          ner (long-array (range (m/long-inc n)))
          ban (double-array (m/long-inc n))
          s (DianaState. d n kwan ner ban (double (reduce m/max 0.0 d)))]
      (Array/set kwan 1 n)
      (loop [ja (long 1) nclu (long 1)]
        (let [jb (m/long-dec (m/long-add ja (Array/get kwan ja)))]
          (split! s ja (m/== nclu 1))
          (when (m/< (m/long-inc nclu) n)
            (recur (next-split kwan ja jb) (m/long-inc nclu)))))
      (let [merges (banner-merges s)]
        {:merge (mapv (fn [[a b]] (canonical-row a b)) merges)
         :height (vec (sort (map #(nth % 2) merges)))
         :order (mapv (fn [^long i] (m/long-dec (Array/get ner i))) (range 1 (m/long-inc n)))
         :dc (banner-coefficient ban n)}))))

;; ## Cutting

(defn clusters-at
  "Clusters (vectors of observations, ascending, ordered by the first observation) of the tree `merge` of `n` observations cut into `k` clusters."
  [merge ^long n ^long k]
  (let [last-kept (m/long-sub n k)
        roots (if (m/== k 1)
                [(root-node merge)]
                (->> (subvec merge last-kept)
                     (mapcat identity)
                     (filter (fn [^long x] (m/<= x last-kept)))))]
    (->> roots
         (map #(vec (sort (subtree-leaves merge %))))
         (sort-by first)
         (vec))))

(defn clusters-for-height
  "Number of clusters of the cut at the height `h` (R `cutree` convention: merges at `h` stay merged)."
  ^long [heights ^long n ^double h]
  (let [idx (or (first (keep-indexed (fn [i ^double x] (when (m/> x h) i)) heights))
                (m/long-dec n))]
    (m/long-sub n (long idx))))

;; ## Medoid

(defn medoid
  "The observation of `indices` with the lowest sum of distances (`pair-distance` of two indices) to the others, the lowest one on ties."
  ^long [indices ^IFn$LLD pair-distance]
  (let [ids (long-array indices)
        size (alength ids)
        sums (double-array size)]
    (dotimes [a size]
      (loop [b (m/long-inc a)]
        (when (m/< b size)
          (let [x (.invokePrim pair-distance (Array/get ids a) (Array/get ids b))]
            (Array/add sums a x)
            (Array/add sums b x))
          (recur (m/long-inc b)))))
    (loop [a (long 1) best (Array/aget sums 0) ba (long 0)]
      (if (m/< a size)
        (let [x (Array/aget sums a)]
          (if (m/< x best) (recur (m/long-inc a) x a) (recur (m/long-inc a) best ba)))
        (Array/get ids ba)))))

;; ## Cophenetic distances and the linkage matrix

(defn cophenetic-distances
  "Condensed cophenetic distances: the height of the first merge joining two observations."
  ^doubles [merge heights ^long n]
  (let [d (double-array (condensed-size n))
        members (fn [step-members ^long x] (if (m/neg? x) [(m/long-sub -1 x)] (step-members (m/long-dec x))))]
    (reduce (fn [step-members [[a b] ^double h]]
              (let [ma (members step-members a)
                    mb (members step-members b)]
                (doseq [^long i ma ^long j mb]
                  (Array/aset d (if (m/< i j) (condensed-index n i j) (condensed-index n j i)) h))
                (conj step-members (into ma mb))))
            [] (map vector merge heights))
    d))

(defn linkage-matrix
  "SciPy linkage matrix rows `[a b height count]` (observations `0..n-1`, the step `j` as `n + j - 1`)."
  [merge heights ^long n]
  (let [id (fn ^long [^long x] (if (m/neg? x) (m/long-sub -1 x) (m/long-add n (m/long-dec x))))]
    (first (reduce (fn [[rows sizes] [[a b] h]]
                     (let [size (fn ^long [^long x] (if (m/neg? x) 1 (sizes (m/long-dec x))))
                           cnt (m/long-add (.invokePrim ^clojure.lang.IFn$LL size a) (.invokePrim ^clojure.lang.IFn$LL size b))
                           ia (.invokePrim ^clojure.lang.IFn$LL id a)
                           ib (.invokePrim ^clojure.lang.IFn$LL id b)]
                       [(conj rows [(m/long-min ia ib) (m/long-max ia ib) h cnt]) (conj sizes cnt)]))
                   [[] []] (map vector merge heights)))))
