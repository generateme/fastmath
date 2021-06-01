;; ## n-dimentional vector utilities
;;
;; Main goal for this namespace is to provide various utility functions which operate on
;; mathematical vectors: 2d, 3d and 4d 
;;
;; Concept for API is taken from [Proceessing](https://github.com/processing/processing/blob/master/core/src/processing/core/PVector.java) and [openFrameworks](https://github.com/openframeworks/openFrameworks/tree/master/libs/openFrameworks/math)
;;
;; All vectors are equipped with Counted (`count`), Sequential, Sequable (`seq`) and IFn protocols. Additionally Clojure vector is equipped with defined here `VectorProto`.

(ns fastmath.vector
  "Mathematical vector operations.

  ### Types

  * Fixed size (custom types):
      * Number - 1d vector
      * Vec2 - 2d vector, creator [[vec2]]
      * Vec3 - 3d vector, creator [[vec3]]
      * Vec4 - 4d vector, creator [[vec4]]
      * ArrayVec - fixed size double array wrapper, n-dimensional, creator [[array-vec]]
  * Fixed size
      * doubles - double array itself
  * Variable size:
      * Clojure's IPersistentVector, creator `[]`
      * Clojure's ISeq

  [[VectorProto]] defines most of the functions.

  Vectors implements also:

  * `Sequable`
  * `Sequencial`
  * `IFn`
  * `Counted`
  * `Reversible`
  * `Indexed`
  * `ILookup`
  * `equals` and `toString` from `Object`
  * `IPersistentVector`
  * `Associative`
  * `clojure.core.matrix.protocols`

  That means that vectors can be destructured, treated as sequence or called as a function. See [[vec2]] for examples."
  {:metadoc/categories {:gen "Creators"
                        :geom "Geometric"
                        :dist "Distance / length"
                        :op "Operations"
                        :mop "Math operations"}}
  (:require [fastmath.core :as m]
            [clojure.string :as s]
            [fastmath.protocols :as prot])
  (:import [clojure.lang Counted IFn ISeq IPersistentVector IPersistentCollection Seqable Sequential Reversible Indexed ILookup
            Associative MapEntry]
           [clojure.core Vec]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; ## Vector definitions

(declare angle-between)
(declare normalize)
(declare div)

(defn- find-idx-reducer-fn 
  "Helper function for reduce to find index for maximum/minimum value in vector."
  [f]
  #(let [[^long midx ^long curr v] %1]
     (if (f %2 v)
       [curr (inc curr) %2]
       [midx (inc curr) v])))

(defn- near-zero?
  "Is your value less than `tol` or `1.0e-6`"
  ([^double tol ^double v]
   (< (m/abs v) tol))
  ([^double v]
   (< (m/abs v) 1.0e-6)))

;; Add `VectorProto` to Clojure vector using map/reduce terms.
(extend ISeq
  prot/VectorProto
  {:to-vec #(apply conj (vector-of :double) %1)
   :as-vec (fn
             ([v xs] (take (count v) xs))
             ([v] (prot/as-vec v (repeat 0.0))))
   :fmap #(map %2 %1)
   :approx (fn
             ([v] (map m/approx v))
             ([v d] (map #(m/approx ^double % d) v)))
   :magsq (fn ^double [v] (reduce (fn ^double [^double b ^double x] (+ b (* x x))) 0.0 v))
   :mag (fn ^double [v] (m/sqrt (prot/magsq v)))
   :dot (fn ^double [v1 v2] (reduce m/fast+ (map m/fast* v1 v2)))
   :add (fn [v1 v2] (map m/fast+ v1 v2))
   :sub (fn [v1 v2] (map m/fast- v1 v2))
   :mult (fn [v1 v] (map #(m/fast* ^double % ^double v) v1))
   :emult #(map m/fast* %1 %2)
   :abs #(map m/abs %)
   :mx #(reduce m/fast-max %)
   :mn #(reduce m/fast-min %)
   :emx #(mapv m/fast-max %1 %2)
   :emn #(mapv m/fast-min %1 %2)
   :maxdim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/>) [0 0 (first v)] v)))
   :mindim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/<) [0 0 (first v)] v)))
   :sum (fn ^double [v] (reduce m/fast+ v))
   :permute #(map (fn [idx] (%1 idx)) %2)
   :reciprocal #(map (fn [^double v] (/ v)) %)
   :heading (fn ^double [v] (angle-between v (reduce conj [1.0] (repeatedly (dec (count v)) (constantly 0.0)))))
   :interpolate (fn [v1 v2 t f] (map #(f %1 %2 t) v1 v2))
   :einterpolate (fn [v1 v2 v f] (map #(f %1 %2 %3) v1 v2 v))
   :econstrain (fn [v val1 val2] (map #(m/constrain ^double %1 ^double val1 ^double val2) v))
   :is-zero? #(every? clojure.core/zero? %)
   :is-near-zero? (fn
                    ([v] (every? near-zero? v))
                    ([v tol] (every? (partial near-zero? tol) v)))})

;; Add `VectorProto` to Clojure vector using mapv/reduce terms.
(extend IPersistentVector
  prot/VectorProto
  {:to-vec #(apply conj (vector-of :double) %1)
   :as-vec (fn
             ([v xs] (vec (take (count v) xs)))
             ([v] (prot/as-vec v (repeat 0.0))))
   :fmap #(mapv %2 %1)
   :approx (fn
             ([v] (mapv m/approx v))
             ([v d] (mapv #(m/approx ^double % d) v)))
   :magsq (fn ^double [v] (reduce (fn ^double [^double b ^double x] (+ b (* x x))) 0.0 v))
   :mag (fn ^double [v] (m/sqrt (prot/magsq v)))
   :dot (fn ^double [v1 v2] (reduce m/fast+ (map m/fast* v1 v2)))
   :add (fn [v1 v2] (mapv m/fast+ v1 v2))
   :sub (fn [v1 v2] (mapv m/fast- v1 v2))
   :mult (fn [v1 v] (mapv #(m/fast* ^double % ^double v) v1))
   :emult #(mapv m/fast* %1 %2)
   :abs #(mapv m/abs %)
   :mx #(reduce m/fast-max %)
   :mn #(reduce m/fast-min %)
   :emx #(mapv m/fast-max %1 %2)
   :emn #(mapv m/fast-min %1 %2)
   :maxdim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/>) [0 0 (first v)] v)))
   :mindim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/<) [0 0 (first v)] v)))
   :sum (fn ^double [v] (reduce m/fast+ v))
   :permute #(mapv (fn [idx] (%1 idx)) %2)
   :reciprocal #(mapv (fn [^double v] (/ v)) %)
   :heading #(angle-between % (reduce conj [1.0] (repeatedly (dec (count %)) (constantly 0.0))))
   :interpolate (fn [v1 v2 t f] (mapv #(f %1 %2 t) v1 v2))
   :einterpolate (fn [v1 v2 v f] (mapv #(f %1 %2 %3) v1 v2 v))
   :econstrain (fn [v val1 val2] (mapv #(m/constrain ^double %1 ^double val1 ^double val2) v))
   :is-zero? #(every? clojure.core/zero? %)
   :is-near-zero? (fn
                    ([v] (every? near-zero? v))
                    ([v tol] (every? (partial near-zero? tol) v)))})

(defn- aevery
  "Array version of every"
  [^doubles arr pred]
  (let [s (alength arr)]
    (loop [idx (unchecked-long 0)]
      (if (< idx s)
        (if (pred (aget arr idx))
          (recur (inc idx))
          false)
        true))))

(extend-type (Class/forName "[D")
  prot/VectorProto
  (to-vec [arr] (let [^Vec v (vector-of :double)]
                  (Vec. (.am v) (alength ^doubles arr) (.shift v) (.root v) arr (.meta v))))
  (as-vec
    ([v xs] (double-array (take (alength ^doubles v) xs)))
    ([v] (prot/as-vec v (repeat 0.0))))
  (fmap [arr f] (amap ^doubles arr idx ret ^double (f (aget ^doubles arr idx))))
  (approx
    ([arr] (amap ^doubles arr idx ret ^double (m/approx (aget ^doubles arr idx))))
    ([arr d] (amap ^doubles arr idx ret ^double (m/approx (aget ^doubles arr idx) d))))
  (magsq [arr] (smile.math.MathEx/dot ^doubles arr ^doubles arr))
  (mag [v1] (m/sqrt (prot/magsq v1)))
  (dot [arr v2] (smile.math.MathEx/dot ^doubles arr ^doubles v2))
  (add [arr v2] (let [b (aclone ^doubles arr)]
                  (smile.math.MathEx/add b ^doubles v2)
                  b))
  (sub [arr v2] (let [b (aclone ^doubles arr)]
                  (smile.math.MathEx/sub b ^doubles v2)
                  b))
  (mult [arr v] (let [b (aclone ^doubles arr)]
                  (smile.math.MathEx/scale ^double v ^doubles b)
                  b))
  (emult [arr v2] (amap ^doubles arr idx ret (* (aget ^doubles arr idx) ^double (v2 idx))))
  (abs [arr] (amap ^doubles arr idx ret (m/abs (aget ^doubles arr idx))))
  (mx [arr] (smile.math.MathEx/max ^doubles arr))
  (mn [arr] (smile.math.MathEx/min ^doubles arr))
  (maxdim [arr] (smile.math.MathEx/whichMax ^doubles arr))
  (mindim [arr] (smile.math.MathEx/whichMin ^doubles arr))
  (emx [arr v2] (amap ^doubles arr idx ret (max (aget ^doubles arr idx) ^double (v2 idx))))
  (emn [arr v2] (amap ^doubles arr idx ret (min (aget ^doubles arr idx) ^double (v2 idx))))
  (sum [arr] (smile.math.MathEx/sum ^doubles arr))
  (heading [arr] (let [v (double-array (alength ^doubles arr) 0.0)]
                   (aset v 0 1.0)
                   (angle-between arr v)))
  (reciprocal [arr] (amap ^doubles arr idx ret (/ (aget ^doubles arr idx))))
  (interpolate [arr v2 t f] (amap ^doubles arr idx ret ^double (f (aget ^doubles arr idx) (v2 idx) t)))
  (einterpolate [arr v2 v f] (amap ^doubles arr idx ret ^double (f (aget ^doubles arr idx) (v2 idx) (v idx))))
  (econstrain [arr val1 val2] (amap ^doubles arr idx ret ^double (m/constrain ^double (aget ^doubles arr idx) ^double val1 ^double val2)))
  (is-zero? [arr] (aevery arr #(zero? ^double %)))
  (is-near-zero?
    ([arr] (aevery arr near-zero?))
    ([arr tol] (aevery arr (partial near-zero? tol)))))

;; Array Vector
(deftype ArrayVec [^doubles array]
  Object
  (toString [_] (str "#arrayvec " (if (> (alength array) 10)
                                    (str "[" (s/join " " (take 10 array)) "...]")
                                    (vec array))))
  (equals [_ v]
    (bool-and (instance? ArrayVec v)
              (smile.math.MathEx/equals array ^doubles (.array ^ArrayVec v) m/MACHINE-EPSILON)))
  (hashCode [_]
    (mix-collection-hash (java.util.Arrays/hashCode array) (alength array)))
  clojure.lang.IHashEq 
  (hasheq [_]
    (mix-collection-hash (java.util.Arrays/hashCode array) (alength array)))
  Sequential
  Seqable
  (seq [_] (seq array))
  Indexed
  (nth [v id] (v id))
  (nth [v id not-found] (or (v id) not-found))
  ILookup
  (valAt [v id] (v id))
  (valAt [v id not-found] (or (v id) not-found))
  Associative
  (containsKey [v id] (< (alength array) (unchecked-int id)))
  (assoc [v k vl] (let [^ArrayVec v v
                        arr (aclone ^doubles (.array v))]
                    (aset arr (unchecked-int k) ^double vl)
                    (ArrayVec. arr)))
  (entryAt [v k] (MapEntry. k (v k)))
  IFn
  (invoke [_ n]
    (aget array ^long n))
  Counted
  (count [_] (alength array))
  IPersistentVector
  (length [_] (alength array))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/VectorProto
  (to-vec [_] (let [^Vec v (vector-of :double)]
                (Vec. (.am v) (alength array) (.shift v) (.root v) array (.meta v))))
  (as-vec [v xs] (ArrayVec. (prot/as-vec array xs)))
  (as-vec [v] (prot/as-vec v (repeat 0.0)))
  (fmap [_ f] (ArrayVec. (prot/fmap array f)))
  (approx [_] (ArrayVec. (prot/approx array)))
  (approx [_ d] (ArrayVec. (prot/approx array d)))
  (magsq [_] (smile.math.MathEx/dot array array))
  (mag [v1] (m/sqrt (prot/magsq v1)))
  (dot [_ v2] (smile.math.MathEx/dot array ^doubles (.array ^ArrayVec v2)))
  (add [_ v2] (ArrayVec. (prot/add array (.array ^ArrayVec v2))))
  (sub [_ v2] (ArrayVec. (prot/sub array (.array ^ArrayVec v2))))
  (mult [_ v] (ArrayVec. (prot/mult array v)))
  (emult [_ v2] (ArrayVec. (prot/emult array v2)))
  (abs [_] (ArrayVec. (prot/abs array)))
  (mx [_] (smile.math.MathEx/max array))
  (mn [_] (smile.math.MathEx/min array))
  (maxdim [_] (smile.math.MathEx/whichMax array))
  (mindim [_] (smile.math.MathEx/whichMin array))
  (emx [_ v2] (ArrayVec. (prot/emx array v2)))
  (emn [_ v2] (ArrayVec. (prot/emn array v2)))
  (sum [_] (smile.math.MathEx/sum array))
  (heading [v1] (prot/heading array))
  (reciprocal [_] (ArrayVec. (prot/reciprocal array)))
  (interpolate [_ v2 t f] (ArrayVec. (prot/interpolate array v2 t f))) 
  (einterpolate [_ v2 v f] (ArrayVec. (prot/einterpolate array v2 v f)))
  (econstrain [_ val1 val2] (ArrayVec. (prot/econstrain array val1 val2)))
  (is-zero? [_] (aevery array #(zero? ^double %)))
  (is-near-zero? [_] (aevery array near-zero?))  
  (is-near-zero? [_ tol] (aevery array (partial near-zero? tol))))

(extend-type Number
  prot/VectorProto
  (to-vec [v] (vector-of :double (double v)))
  (as-vec
    ([_] 0.0)
    ([_ xs] (double (first xs))))
  (fmap [v f] (f v))
  (approx
    ([v] (m/approx v))
    ([v d] (m/approx v d)))
  (magsq [v] (m/abs v))
  (mag [v] (m/sq v))
  (dot [v1 v2] (* (double v1) (double v2)))
  (add [v1 v2] (+ (double v1) (double v2)))
  (sub [v1 v2] (- (double v1) (double v2)))
  (mult [v1 ^double v] (* (double v1) v))
  (emult [v1 v2] (* (double v1) (double v2)))
  (abs [v] (m/abs v))
  (mx [v] v)
  (mn [v] v)
  (emx [v1 v2] (max (double v1) (double v2)))
  (emn [v1 v2] (max (double v1) (double v2)))
  (maxdim [_] 0)
  (mindim [_] 0)
  (sum [v] v)
  (reciprocal [v] (/ (double v)))
  (interpolate [v1 v2 t f] (f v1 v2 t))
  (einterpolate [v1 v2 t f] (f v1 v2 t))
  (econstrain [v val1 val2] (m/constrain (double v) (double val1) (double val2)))
  (is-zero? [v] (zero? (double v)))
  (is-near-zero?
    ([v] (near-zero? v))
    ([v tol] (near-zero? tol v))))

(defn- dhash-code
  "double hashcode"
  (^long [^long state ^double a]
   (let [abits (Double/doubleToLongBits a)
         elt (bit-xor abits (>>> abits 32))]
     (+ elt (* 31 state))))
  (^long [^double a]
   (let [abits (Double/doubleToLongBits a)
         elt (bit-xor abits (>>> abits 32))]
     (+ elt 31))))

;; Create Vec2 and add all necessary protocols
(deftype Vec2 [^double x ^double y]
  Object
  (toString [_] (str "#vec2 [" x ", " y "]"))
  (equals [_ v]
    (bool-and (instance? Vec2 v)
              (let [^Vec2 v v]
                (bool-and (== x (.x v))
                          (== y (.y v))))))
  (hashCode [_] (mix-collection-hash (unchecked-int (dhash-code (dhash-code x) y)) 2))
  clojure.lang.IHashEq 
  (hasheq [_] (mix-collection-hash (unchecked-int (dhash-code (dhash-code x) y)) 2))
  Sequential
  Seqable
  (seq [_] (list x y))
  Reversible
  (rseq [_] (list y x))
  Indexed
  (nth [v id] (v id))
  (nth [v id not-found] (or (v id) not-found))
  ILookup
  (valAt [v id] (v id))
  (valAt [v id not-found] (or (v id) not-found))
  Associative
  (containsKey [v id] (#{0 1} id))
  (assoc [v k vl] (case (unchecked-int k)
                    0 (Vec2. vl y)
                    1 (Vec2. x vl)
                    vl))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 2)
  IFn
  (invoke [_ id]
    (case (unchecked-int id)
      0 x
      1 y
      nil))
  IPersistentVector
  (length [_] 2)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/VectorProto
  (to-vec [_] (vector-of :double x y))
  (as-vec [_ [x y]] (Vec2. x y))
  (as-vec [_] (Vec2. 0.0 0.0))
  (fmap [_ f] (Vec2. (f x) (f y)))
  (approx [_] (Vec2. (m/approx x) (m/approx y)))
  (approx [_ d] (Vec2. (m/approx x d) (m/approx y d)))
  (magsq [_] (+ (* x x) (* y y)))
  (mag [_] (m/hypot-sqrt x y))
  (dot [_ v2] 
    (let [^Vec2 v2 v2] (+ (* x (.x v2)) (* y (.y v2)))))
  (add [_ v2] 
    (let [^Vec2 v2 v2] (Vec2. (+ x (.x v2)) (+ y (.y v2)))))
  (sub [_ v2]
    (let [^Vec2 v2 v2] (Vec2. (- x (.x v2)) (- y (.y v2)))))
  (mult [_ v] (Vec2. (* x ^double v) (* y ^double v)))
  (emult [_ v] 
    (let [^Vec2 v v] (Vec2. (* x (.x v)) (* y (.y v)))))
  (abs [_] (Vec2. (m/abs x) (m/abs y)))
  (mx [_] (max x y))
  (mn [_] (min x y))
  (emx [_ v]
    (let [^Vec2 v v] (Vec2. (max (.x v) x) (max (.y v) y))))
  (emn [_ v]
    (let [^Vec2 v v] (Vec2. (min (.x v) x) (min (.y v) y))))
  (maxdim [_]
    (if (> x y) 0 1))
  (mindim [_]
    (if (< x y) 0 1))
  (base-from [v]
    [v (prot/perpendicular v)])
  (sum [_] (+ x y))
  (permute [p [^long i1 ^long i2]]
    (Vec2. (p i1) (p i2)))
  (reciprocal [_] (Vec2. (/ x) (/ y)))
  (interpolate [_ v2 t f]
    (let [^Vec2 v2 v2] (Vec2. (f x (.x v2) t)
                              (f y (.y v2) t))))
  (einterpolate [_ v2 v f]
    (let [^Vec2 v2 v2
          ^Vec2 v v]
      (Vec2. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v)))))
  (econstrain [_ val1 val2] (Vec2. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y)))
  (is-near-zero? [_] (m/bool-and (near-zero? x) (near-zero? y)))
  (is-near-zero? [_ tol] (m/bool-and (near-zero? tol x) (near-zero? tol y)))
  (heading [_] (m/atan2 y x))
  (cross [_ v]
    (let [^Vec2 v v]
      (- (* x (.y v)) (* y (.x v)))))
  (rotate [_ angle]
    (let [sa (m/sin angle)
          ca (m/cos angle)
          nx (- (* x ca) (* y sa))
          ny (+ (* x sa) (* y ca))]
      (Vec2. nx ny)))
  (perpendicular [_]
    (normalize (Vec2. (- y) x)))
  (transform [_ o vx vy]
    (let [^Vec2 o o
          ^Vec2 vx vx
          ^Vec2 vy vy]
      (Vec2. (+ (.x o) (* x (.x vx)) (* y (.x vy))) (+ (.y o) (* x (.y vx)) (* y (.y vy))))))
  (to-polar [v]
    (Vec2. (prot/mag v) (prot/heading v)))
  (from-polar [_]
    (Vec2. (* x (m/cos y))
           (* x (m/sin y)))))

;; Create Vec3 and add all necessary protocols
(deftype Vec3 [^double x ^double y ^double z]
  Object
  (toString [_] (str "#vec3 [" x ", " y ", " z "]"))
  (equals [_ v]
    (bool-and (instance? Vec3 v)
              (let [^Vec3 v v]
                (bool-and (== x (.x v))
                          (== y (.y v))
                          (== z (.z v))))))
  (hashCode [_] (mix-collection-hash (unchecked-int (dhash-code (dhash-code (dhash-code x) y) z)) 3))
  clojure.lang.IHashEq 
  (hasheq [_] (mix-collection-hash (unchecked-int (dhash-code (dhash-code (dhash-code x) y) z)) 3))
  Sequential
  Seqable
  (seq [_] (list x y z))
  Reversible
  (rseq [_] (list z y x))
  Indexed
  (nth [v id] (v id))
  (nth [v id not-found] (or (v id) not-found))
  ILookup
  (valAt [v id] (v id))
  (valAt [v id not-found] (or (v id) not-found))
  Associative
  (containsKey [v id] (#{0 1 2} id))
  (assoc [v k vl] (case (unchecked-int k)
                    0 (Vec3. vl y z)
                    1 (Vec3. x vl z)
                    2 (Vec3. x y vl)
                    vl))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 3)
  IFn
  (invoke [_ id]
    (case (unchecked-int id)
      0 x
      1 y
      2 z
      nil))
  IPersistentVector
  (length [_] 3)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/VectorProto
  (to-vec [_] (vector-of :double x y z))
  (as-vec [_ [x y z]] (Vec3. x y z))
  (as-vec [_] (Vec3. 0.0 0.0 0.0))
  (fmap [_ f] (Vec3. (f x) (f y) (f z)))
  (approx [_] (Vec3. (m/approx x) (m/approx y) (m/approx z)))
  (approx [_ d] (Vec3. (m/approx x d) (m/approx y d) (m/approx z d)))
  (magsq [_] (+ (* x x) (* y y) (* z z)))
  (mag [_] (m/hypot-sqrt x y z))
  (dot [_ v2]
    (let [^Vec3 v2 v2] (+ (* x (.x v2)) (* y (.y v2)) (* z (.z v2)))))
  (add [_ v2] 
    (let [^Vec3 v2 v2] (Vec3. (+ x (.x v2)) (+ y (.y v2)) (+ z (.z v2)))))
  (sub [_ v2]
    (let [^Vec3 v2 v2] (Vec3. (- x (.x v2)) (- y (.y v2)) (- z (.z v2)))))
  (mult [_ v] (Vec3. (* x ^double v) (* y ^double v) (* z ^double v)))
  (emult [_ v] 
    (let [^Vec3 v v] (Vec3. (* x (.x v)) (* y (.y v)) (* z (.z v)))))
  (abs [_] (Vec3. (m/abs x) (m/abs y) (m/abs z)))
  (mx [_] (max x y z))
  (mn [_] (min x y z))
  (emx [_ v]
    (let [^Vec3 v v] (Vec3. (max (.x v) x) (max (.y v) y) (max (.z v) z))))
  (emn [_ v]
    (let [^Vec3 v v] (Vec3. (min (.x v) x) (min (.y v) y) (min (.z v) z))))
  (maxdim [_]
    (if (> x y)
      (if (> x z) 0 2)
      (if (> y z) 1 2)))
  (mindim [_]
    (if (< x y)
      (if (< x z) 0 2)
      (if (< y z) 1 2)))
  (base-from [v]
    (let [v2 (if (> (m/abs x) (m/abs y))
               (div (Vec3. (- z) 0.0 x) (m/hypot-sqrt x z))
               (div (Vec3. 0.0 z (- y)) (m/hypot-sqrt y z)))]
      [v v2 (prot/cross v v2)]))
  (sum [_] (+ x y z))
  (permute [p [^long i1 ^long i2 ^long i3]]
    (Vec3. (p i1) (p i2) (p i3)))
  (reciprocal [_] (Vec3. (/ x) (/ y) (/ z)))
  (interpolate [_ v2 t f]
    (let [^Vec3 v2 v2] (Vec3. (f x (.x v2) t)
                              (f y (.y v2) t)
                              (f z (.z v2) t))))
  (einterpolate [_ v2 v f]
    (let [^Vec3 v2 v2
          ^Vec3 v v]
      (Vec3. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v))
             (f z (.z v2) (.z v)))))
  (econstrain [_ val1 val2] (Vec3. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)
                                   (m/constrain z ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y) (zero? z)))
  (is-near-zero? [_] (bool-and (near-zero? x) (near-zero? y) (near-zero? z)))
  (is-near-zero? [_ tol] (bool-and (near-zero? tol x) (near-zero? tol y) (near-zero? tol z)))
  (heading [v1] (angle-between v1 (Vec3. 1 0 0)))
  (cross [_ v2]
    (let [^Vec3 v2 v2
          cx (- (* y (.z v2)) (* (.y v2) z))
          cy (- (* z (.x v2)) (* (.z v2) x))
          cz (- (* x (.y v2)) (* (.x v2) y))]
      (Vec3. cx cy cz)))
  (perpendicular [v1 v2] (normalize (prot/cross v1 v2)))
  (transform [_ o vx vy vz]
    (let [^Vec3 o o
          ^Vec3 vx vx
          ^Vec3 vy vy
          ^Vec3 vz vz]
      (Vec3. (+ (.x o) (* x (.x vx)) (* y (.x vy)) (* z (.x vz)))
             (+ (.y o) (* x (.y vx)) (* y (.y vy)) (* z (.y vz)))
             (+ (.z o) (* x (.z vx)) (* y (.z vy)) (* z (.z vz))))))
  (axis-rotate [_ angle axis]
    (let [^Vec3 axis axis
          ^Vec3 ax (normalize axis)
          axx (.x ax)
          axy (.y ax)
          axz (.z ax)
          cosa (m/cos angle)
          ^Vec3 sa (prot/mult ax (m/sin angle))
          sax (.x sa)
          say (.y sa)
          saz (.z sa)
          ^Vec3 cb (prot/mult ax (- 1.0 cosa))
          cbx (.x cb)
          cby (.y cb)
          cbz (.z cb)
          nx (+ (* x (+ (* axx cbx) cosa))
                (* y (- (* axx cby) saz))
                (* z (+ (* axx cbz) say)))
          ny (+ (* x (+ (* axy cbx) saz))
                (* y (+ (* axy cby) cosa))
                (* z (- (* axy cbz) sax)))
          nz (+ (* x (- (* axz cbx) say))
                (* y (+ (* axz cby) sax))
                (* z (+ (* axz cbz) cosa)))]
      (Vec3. nx ny nz)))
  (axis-rotate [v1 angle axis pivot]
    (prot/add (prot/axis-rotate (prot/sub v1 pivot) angle axis) pivot))
  (rotate [_ anglex angley anglez]
    (let [a (m/cos anglex)
          b (m/sin anglex)
          c (m/cos angley)
          d (m/sin angley)
          e (m/cos anglez)
          f (m/sin anglez)
          cex (* c x e)
          cf (* c f)
          dz (* d z)
          nx (+ (- cex cf) dz)
          af (* a f)
          de (* d e)
          bde (* b de)
          ae (* a e)
          bdf (* b d f)
          bcz (* b c z)
          ny (- (+ (* (+ af bde) x) (* (- ae bdf) y)) bcz)
          bf (* b f)
          ade (* a de)
          adf (* a d f)
          be (* b e)
          acz (* a c z)
          nz (+ (* (- bf ade) x) (* (+ adf be) y) acz)]
      (Vec3. nx ny nz)))
  (to-polar [v1]
    (let [^double r (prot/mag v1)
          zr (/ z r)
          theta (cond
                  (<= zr -1.0) m/PI
                  (>= zr 1.0) 0
                  :else (m/acos zr))
          phi (m/atan2 y x)]
      (Vec3. r theta phi)))
  (from-polar [_]
    (let [st (m/sin y)
          ct (m/cos y)
          sp (m/sin z)
          cp (m/cos z)]
      (Vec3. (* x st cp)
             (* x st sp)
             (* x ct)))))

;; Create Vec4 and add all necessary protocols
(deftype Vec4 [^double x ^double y ^double z ^double w]
  Object
  (toString [_] (str "#vec4 [" x ", " y ", " z ", " w "]"))
  (equals [_ v]
    (and (instance? Vec4 v)
         (let [^Vec4 v v]
           (bool-and (== x (.x v))
                     (== y (.y v))
                     (== z (.z v))
                     (== w (.w v))))))
  (hashCode [_]
    (mix-collection-hash (unchecked-int (dhash-code (dhash-code (dhash-code (dhash-code x) y) z) w)) 4))
  clojure.lang.IHashEq 
  (hasheq [_]
    (mix-collection-hash (unchecked-int (dhash-code (dhash-code (dhash-code (dhash-code x) y) z) w)) 4))
  Sequential
  Seqable
  (seq [_] (list x y z w))
  Reversible
  (rseq [_] (list w z y x))
  Indexed
  (nth [v id] (v id))
  (nth [v id not-found] (or (v id) not-found))
  ILookup
  (valAt [v id] (v id))
  (valAt [v id not-found] (or (v id) not-found))
  Associative
  (containsKey [v id] (#{0 1 2 3} id))
  (assoc [v k vl] (case (unchecked-int k)
                    0 (Vec4. vl y z w)
                    1 (Vec4. x vl z w)
                    2 (Vec4. x y vl w)
                    3 (Vec4. x y z vl)
                    vl))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 4)
  IFn
  (invoke [_ id]
    (case (unchecked-int id)
      0 x
      1 y
      2 z
      3 w
      nil))
  IPersistentVector
  (length [_] 4)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/VectorProto
  (to-vec [_] (vector-of :double x y z w))
  (as-vec [_ [x y z w]] (Vec4. x y z w))
  (as-vec [_] (Vec4. 0.0 0.0 0.0 0.0))
  (fmap [_ f] (Vec4. (f x) (f y) (f z) (f w)))
  (approx [_] (Vec4. (m/approx x) (m/approx y) (m/approx z) (m/approx w)))
  (approx [_ d] (Vec4. (m/approx x d) (m/approx y d) (m/approx z d) (m/approx w d)))
  (magsq [_] (+ (* x x) (* y y) (* z z) (* w w)))
  (mag [v1] (m/sqrt (prot/magsq v1)))
  (dot [_ v2]
    (let [^Vec4 v2 v2] (+ (* x (.x v2)) (* y (.y v2)) (* z (.z v2)) (* w (.w v2)))))
  (add [_ v2]
    (let [^Vec4 v2 v2] (Vec4. (+ x (.x v2)) (+ y (.y v2)) (+ z (.z v2)) (+ w (.w v2)))))
  (sub [_ v2] 
    (let [^Vec4 v2 v2] (Vec4. (- x (.x v2)) (- y (.y v2)) (- z (.z v2)) (- w (.w v2)))))
  (mult [_ v] (Vec4. (* x ^double v) (* y ^double v) (* z ^double v) (* w ^double v)))
  (emult [_ v]
    (let [^Vec4 v v] (Vec4. (* x (.x v)) (* y (.y v)) (* z (.z v)) (* w (.w v)))))
  (abs [_] (Vec4. (m/abs x) (m/abs y) (m/abs z) (m/abs w)))
  (mx [_] (max x y z w))
  (mn [_] (min x y z w))
  (emx [_ v]
    (let [^Vec4 v v] (Vec4. (max (.x v) x) (max (.y v) y) (max (.z v) z) (max (.w v) w))))
  (emn [_ v]
    (let [^Vec4 v v] (Vec4. (min (.x v) x) (min (.y v) y) (min (.z v) z) (min (.w v) w))))
  (maxdim [_]
    (max-key [x y z w] 0 1 2 3))
  (mindim [_]
    (min-key [x y z w] 0 1 2 3))
  (sum [_] (+ x y z w))
  (permute [p [^long i1 ^long i2 ^long i3 ^long i4]]
    (Vec4. (p i1) (p i2) (p i3) (p i4)))
  (reciprocal [_] (Vec4. (/ x) (/ y) (/ z) (/ w)))
  (interpolate [_ v2 t f]
    (let [^Vec4 v2 v2] (Vec4. (f x (.x v2) t)
                              (f y (.y v2) t)
                              (f z (.z v2) t)
                              (f w (.w v2) t))))
  (einterpolate [_ v2 v f]
    (let [^Vec4 v2 v2
          ^Vec4 v v]
      (Vec4. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v))
             (f z (.z v2) (.z v))
             (f w (.w v2) (.w v)))))
  (econstrain [_ val1 val2] (Vec4. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)
                                   (m/constrain z ^double val1 ^double val2)
                                   (m/constrain w ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y) (zero? z) (zero? w)))
  (is-near-zero? [_] (bool-and (near-zero? x) (near-zero? y) (near-zero? z) (near-zero? w)))
  (is-near-zero? [_ tol] (bool-and (near-zero? tol x) (near-zero? tol y) (near-zero? tol z) (near-zero? tol w)))
  (heading [v1] (angle-between v1 (Vec4. 1 0 0 0))))

;;

(def ^{:deprecated "v1.3.0" :metadoc/categories #{:op} :doc "Same as [[fmap]]. Deprecated."} applyf prot/fmap)
(def ^{:deprecated "v1.5.0" :metadoc/categories #{:gen} :doc "Same as [[vec->Vec]]. Deprecated."} to-vec prot/to-vec)

;; protocol methods mapped

(defn vec->Vec
  "Convert to Clojure primitive vector `Vec`."
  {:metadoc/categories #{:gen}}
  [v] 
  (prot/to-vec v))

(defn as-vec
  "Create vector from sequence as given type. If there is no sequence fill with `0.0`."
  {:metadoc/categories #{:gen}}
  ([v] (prot/as-vec v))
  ([v xs] (prot/as-vec v xs)))

(defn fmap
  "Apply function to all vector values (like map but returns the same type)."
  {:metadoc/categories #{:op}}
  [v f]
  (prot/fmap v f))

(defn approx
  "Round to 2 (or `d`) decimal places"
  {:metadoc/categories #{:op}}
  ([v] (prot/approx v))
  ([v d] (prot/approx v d)))

(defn magsq
  "Length of the vector squared."
  {:metadoc/categories #{:dist :geom}}
  ^double [v] (prot/magsq v))

(defn mag
  "Length of the vector."
  {:metadoc/categories #{:dist :geom}}
  ^double [v] (prot/mag v))

(defn dot
  "Dot product of two vectors."
  {:metadoc/categories #{:geom}}
  ^double [v1 v2] (prot/dot v1 v2))

(defn add
  "Sum of two vectors."
  {:metadoc/categories #{:op}}
  ([v] v)
  ([v1 v2] (prot/add v1 v2)))

(defn sub
  "Subtraction of two vectors."
  {:metadoc/categories #{:op}}
  ([v] (prot/mult v -1.0))
  ([v1 v2] (prot/sub v1 v2)))

(defn mult
  "Multiply vector by number `x`."
  {:metadoc/categories #{:op}}
  [v x] (prot/mult v x))

(defn emult
  "Element-wise vector multiplication (Hadamard product)."
  {:metadoc/categories #{:op}}
  [v1 v2] (prot/emult v1 v2))

(defn abs
  "Absolute value of vector elements"
  {:metadoc/categories #{:op}}
  [v] (prot/abs v))

(defn mx
  "Maximum value of vector elements"
  {:metadoc/categories #{:op}}
  ^double [v] (prot/mx v))

(defn mn
  "Minimum value of vector elements"
  {:metadoc/categories #{:op}}
  ^double [v] (prot/mn v))

(defn emx
  "Element-wise max from two vectors."
  {:metadoc/categories #{:op}}
  [v1 v2] (prot/emx v1 v2))

(defn emn
  "Element-wise min from two vectors."
  {:metadoc/categories #{:op}}
  [v1 v2] (prot/emn v1 v2))

(defn maxdim
  "Index of maximum value."
  {:metadoc/categories #{:op}}
  ^long [v] (prot/maxdim v))

(defn mindim
  "Index of minimum value."
  {:metadoc/categories #{:op}}
  ^long [v] (prot/mindim v))

(defn base-from
  "List of perpendicular vectors (basis). Works only for `Vec2` and `Vec3` types."
  {:metadoc/categories #{:geom}}
  [v] (prot/base-from v))

(defn sum
  "Sum of elements"
  {:metadoc/categories #{:op}}
  ^double [v] (prot/sum v))

(defn permute
  "Permute vector elements with given indices."
  {:metadoc/categories #{:op}}
  [v idxs] (prot/permute v idxs))

(defn reciprocal
  "Reciprocal of elements."
  {:metadoc/categories #{:op}}
  [v] (prot/reciprocal v))

(defn interpolate 
  "Interpolate vectors, optionally set interpolation fn (default: lerp)"
  {:metadoc/categories #{:op}}
  ([v1 v2 t] (prot/interpolate v1 v2 t m/lerp))
  ([v1 v2 t f] (prot/interpolate v1 v2 t f)))

(defn einterpolate 
  "Interpolate vector selement-wise, optionally set interpolation fn (default: lerp)"
  {:metadoc/categories #{:op}}
  ([v1 v2 v] (prot/einterpolate v1 v2 v m/lerp))
  ([v1 v2 v f] (prot/einterpolate v1 v2 v f)))

(defn econstrain
  "Element-wise constrain"
  {:metadoc/categories #{:op}}
  [v mn mx] (prot/econstrain v mn mx))

(defn is-zero?
  "Is vector zero?"
  {:metadoc/categories #{:op}}
  [v] (prot/is-zero? v))

(defn is-near-zero?
  "Is vector almost zero? (all absolute values of elements are less than `tol` tolerance or `1.0e-6`)"
  {:metadoc/categories #{:op}}
  ([v] (prot/is-near-zero? v))
  ([v tol] (prot/is-near-zero? v tol)))

(defn heading
  "Angle between vector and unit vector `[1,0,...]`"
  {:metadoc/categories #{:geom}}
  ^double [v] (prot/heading v))

(defn cross
  "Cross product"
  {:metadoc/categories #{:geom}}
  [v1 v2] (prot/cross v1 v2))

(defn rotate
  "Rotate vector. Only for `Vec2` and `Vec3` types."
  {:metadoc/categories #{:geom}}
  ([v angle] (prot/rotate v angle))
  ([v angle-x angle-y angle-z] (prot/rotate v angle-x angle-y angle-z)))

(defn axis-rotate
  "Rotate vector. Only for `Vec3` types"
  {:metadoc/categories #{:geom}}
  ([v angle axis] (prot/axis-rotate v angle axis))
  ([v angle axis pivot] (prot/axis-rotate v angle axis pivot)))

(defn perpendicular
  "Perpendicular vector. Only for `Vec2` and `Vec3` types."
  {:metadoc/categories #{:geom}}
  ([v] (prot/perpendicular v))
  ([v1 v2] (prot/perpendicular v1 v2)))

(defn transform
  "Transform vector; map point to coordinate system defined by origin, vx and vy (as bases), Only for `Vec2` and `Vec3` types."
  {:metadoc/categories #{:geom}}
  ([v o vx vy] (prot/transform v o vx vy))
  ([v o vx vy vz] (prot/transform v o vx vy vz)))

(defn to-polar
  "To polar coordinates (2d, 3d only), first element is length, the rest angle."
  {:metadoc/categories #{:geom}}
  [v] (prot/to-polar v))

(defn from-polar
  "From polar coordinates (2d, 3d only)"
  {:metadoc/categories #{:geom}}
  [v] (prot/from-polar v))


;; creators

(defn vec2
  "Make 2d vector."
  {:metadoc/categories #{:gen}} 
  ([x y] (Vec2. x y))
  ([] (Vec2. 0.0 0.0)))

(defn vec3
  "Make Vec2 vector"
  {:metadoc/categories #{:gen}} 
  ([x y z] (Vec3. x y z))
  ([^Vec2 v z] (Vec3. (.x v) (.y v) z))
  ([] (Vec3. 0.0 0.0 0.0)))

(defn vec4
  "Make Vec4 vector"
  {:metadoc/categories #{:gen}} 
  ([x y z w] (Vec4. x y z w))
  ([^Vec3 v w] (Vec4. (.x v) (.y v) (.z v) w))
  ([^Vec2 v z w] (Vec4. (.x v) (.y v) z w))
  ([] (Vec4. 0.0 0.0 0.0 0.0)))

(defn array-vec
  "Make ArrayVec type based on provided sequence `xs`."
  {:metadoc/categories #{:gen}} 
  [xs]
  (ArrayVec. (double-array xs)))

(defn make-vector
  "Returns fixed size vector for given number of dimensions.

  Proper type is used."
  {:metadoc/categories #{:gen}} 
  ([dims xs] (prot/as-vec (make-vector dims) xs))
  ([^long dims]
   (when (pos? dims)
     (case dims
       2 (vec2)
       3 (vec3)
       4 (vec4)
       (array-vec dims)))))

;; ## Common vector functions

(defn div
  "Vector division or reciprocal."
  {:metadoc/categories #{:op}}
  ([v1 ^double v] (prot/mult v1 (/ v)))
  ([v1] (prot/reciprocal v1)))

(defn ediv
  "Element-wise division of two vectors."
  {:metadoc/categories #{:op}} 
  [v1 v2]
  (prot/emult v1 (prot/reciprocal v2)))

(defn zero-count
  "Count zeros in vector"
  {:metadoc/categories #{:op}}
  [v]
  (count (filter #(zero? ^double %) v)))

(defn clamp
  "Clamp elements."
  {:metadoc/categories #{:op}}
  ([v mn mx] (prot/econstrain v mn mx))
  ([v] (prot/econstrain v 0 Double/MAX_VALUE)))

(defn nonzero-count
  "Count non zero velues in vector"
  {:metadoc/categories #{:op}}
  [v]
  (count (remove (fn [^double v] (zero? v)) v)))

(defn average-vectors
  "Average / centroid of vectors. Input: initial vector (optional), list of vectors"
  {:metadoc/categories #{:op}} 
  ([init vs]
   (div (reduce prot/add init vs) (inc (count vs))))
  ([vs] (average-vectors (first vs) (rest vs))))

(defn dist
  "Euclidean distance between vectors"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (mag (prot/sub v1 v2)))

(defn dist-sq
  "Squared Euclidean distance between vectors"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (magsq (prot/sub v1 v2)))

(defn dist-abs
  "Manhattan distance between vectors"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (sum (prot/abs (prot/sub v1 v2))))

(defn dist-cheb
  "Chebyshev distance between 2d vectors"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (mx (prot/abs (prot/sub v1 v2))))

(defn dist-discrete
  "Discrete distance between 2d vectors"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (sum (prot/fmap (prot/sub v1 v2) (fn [^double v] (if (zero? v) 0.0 1.0)))))

(defn dist-canberra
  "Canberra distance"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (let [num (prot/abs (prot/sub v1 v2))
        denom (prot/fmap (prot/add (prot/abs v1) (prot/abs v2)) (fn [^double v] (if (zero? v) 0.0 (/ v))))]
    (sum (prot/emult num denom))))

(defn dist-emd
  "Earth Mover's Distance"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (first (reduce #(let [[^double s ^double l] %1
                        [^double a ^double b] %2
                        n (- (+ a l) b)]
                    [(+ s (m/abs l)) n]) [0.0 0.0] (map vector v1 v2))))

(defn dist-ang
  "Angular distance"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (* (m/acos (/ (dot v1 v2) (* (mag v1) (mag v2)))) m/M_1_PI))

(defn sim-cos
  "Cosine similarity"
  {:metadoc/categories #{:dist}} 
  ^double [v1 v2]
  (/ (dot v1 v2) (* (mag v1) (mag v2))))


;; List of distance fn
(def ^{:metadoc/categories #{:dist}}
  distances {:euclid dist
             :euclid-sq dist-sq
             :abs dist-abs
             :cheb dist-cheb
             :canberra dist-canberra
             :emd dist-emd
             :angular dist-ang
             :discrete dist-discrete})

(defn normalize
  "Normalize vector (set length = 1.0)"
  {:metadoc/categories #{:dist}} 
  [v]
  (let [m (mag v)]
    (if (zero? m)
      (Vec2. 0.0 0.0)
      (div v m))))

(defn set-mag
  "Set length of the vector"
  {:metadoc/categories #{:dist}} 
  [v len]
  (prot/mult (normalize v) len))

(defn limit
  "Limit length of the vector by given value"
  {:metadoc/categories #{:dist}} 
  [v ^double len]
  (if (> (magsq v) (* len len))
    (set-mag v len)
    v))

(defn angle-between
  "Angle between two vectors

  See also [[relative-angle-between]]."
  {:metadoc/categories #{:geom}} 
  ^double [v1 v2]
  (if (or (is-zero? v1) (is-zero? v2))
    0
    (let [d (dot v1 v2)
          amt (/ d (* (mag v1) (mag v2)))]
      (cond
        (<= amt -1.0) m/PI
        (>= amt 1.0) 0
        :else (m/acos amt)))))

(defn relative-angle-between
  "Angle between two vectors relative to each other.

  See also [[angle-between]]."
  {:metadoc/categories #{:geom}} 
  ^double [v1 v2]
  (- (heading v2) (heading v1)))

(defn aligned?
  "Are vectors aligned (have the same direction)?"
  {:metadoc/categories #{:geom}}
  ([v1 v2 ^double tol]
   (< (angle-between v1 v2) tol))
  ([v1 v2]
   (< (angle-between v1 v2) 1.0e-6)))

(defn faceforward
  "Flip normal `n` to match the same direction as `v`."
  {:metadoc/categories #{:geom}} 
  [n v]
  (if (neg? (dot n v)) 
    (sub n)
    n))

(defn generate-vec2
  "Generate Vec2 with fn(s)"
  {:metadoc/categories #{:gen}} 
  ([f1 f2]
   (Vec2. (f1) (f2)))
  ([f]
   (Vec2. (f) (f))))

(defn generate-vec3
  "Generate Vec3 with fn(s)"
  {:metadoc/categories #{:gen}} 
  ([f1 f2 f3]
   (Vec3. (f1) (f2) (f3)))
  ([f]
   (Vec3. (f) (f) (f))))

(defn generate-vec4
  "Generate Vec4 with fn(s)"
  {:metadoc/categories #{:gen}} 
  ([f1 f2 f3 f4]
   (Vec4. (f1) (f2) (f3) (f4)))
  ([f]
   (Vec4. (f) (f) (f) (f))))

(defn array->vec2
  "Doubles array to Vec2"
  {:metadoc/categories #{:gen}} 
  [^doubles arr]
  (Vec2. (aget arr 0) (aget arr 1)))

(defn array->vec3
  "Doubles array to Vec3"
  {:metadoc/categories #{:gen}} 
  [^doubles arr]
  (Vec3. (aget arr 0) (aget arr 1) (aget arr 2)))

(defn array->vec4
  "Doubles array to Vec4"
  {:metadoc/categories #{:gen}} 
  [^doubles arr]
  (Vec4. (aget arr 0) (aget arr 1) (aget arr 2) (aget arr 3)))

(defn seq->vec2
  "Any seq to Vec2"
  {:metadoc/categories #{:gen}} 
  [xs]
  (Vec2. (nth xs 0 0.0) (nth xs 1 0.0)))

(defn seq->vec3
  "Any seq to Vec3"
  {:metadoc/categories #{:gen}} 
  [xs]
  (Vec3. (nth xs 0 0.0) (nth xs 1 0.0) (nth xs 2 0.0)))

(defn seq->vec4
  "Any seq to Vec4"
  {:metadoc/categories #{:gen}} 
  [xs]
  (Vec4. (nth xs 0 0.0) (nth xs 1 0.0) (nth xs 2 0.0) (nth xs 3 0.0)))

;; primitive functions

(defmacro ^:private primitive-ops
  "Generate primitive functions operating on vectors"
  [fns]
  (let [v (symbol "vector")]
    `(do ~@(for [f fns
                 :let [nm (with-meta (symbol (name f)) {:metadoc/categories #{:mop}})
                       doc (str "Apply " nm " to vector elements.")
                       wfn (if (:macro (meta (resolve f))) `(fn [v#] (~f v#)) f)]] ;; wrap macro into function
             `(defn ~nm ~doc
                [~v]
                (prot/fmap ~v ~wfn))))))

(primitive-ops [m/sin m/cos m/tan m/asin m/acos m/atan m/sinh m/cosh m/tanh m/asinh m/acosh m/atanh
                m/cot m/sec m/csc m/acot m/asec m/acsc m/coth m/sech m/csch m/acoth m/asech m/csch
                m/sq m/safe-sqrt m/sqrt m/cbrt m/exp m/log m/log10 m/log2 m/ln m/log1p m/expm1
                m/radians m/degrees m/sinc m/sigmoid
                m/floor m/ceil m/round m/rint m/trunc m/frac m/sfrac m/signum m/sgn])
;;

(defmethod print-method Vec2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec2 [v w] (print-method v w))

(defmethod print-method Vec3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec3 [v w] (print-method v w))

(defmethod print-method Vec4 [^Vec4 v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec4 [v w] (print-method v w))
