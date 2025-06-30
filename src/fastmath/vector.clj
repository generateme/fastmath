(ns fastmath.vector
  "Provides mathematical vector operations for various vector types, inspired by libraries like Processing and openFrameworks.

  Supported vector representations include:

  *   **Fixed-size custom types**: `Vec2` (2D), `Vec3` (3D), `Vec4` (4D), and `ArrayVec` (N-dimensional wrapper around `double[]`).
  *   **Built-in array type**: `double` arrays (`[D`).
  *   **Variable-size Clojure collections**: persistent vectors (`[]`) and sequences (`ISeq`).
  *   **Numbers**: Treated as 1D vectors.

  Most operations are defined via the [[VectorProto]] protocol, which is extended for all supported types.

  These vector types also implement standard Clojure/Java protocols like `Seqable`, `Indexed`, `IFn`, `Counted`, `Associative`, `IReduce`, `ILookup`, etc., enabling them to be used flexibly as sequences, arrays, or functions."
  (:refer-clojure :exclude [abs zero?])
  (:require [fastmath.core :as m]
            [clojure.string :as s]
            [fastmath.protocols :as prot])
  (:import [clojure.lang Counted IFn IPersistentVector IPersistentCollection Seqable Sequential Reversible Indexed ILookup Associative MapEntry IReduce IReduceInit IPersistentStack]
           [clojure.core Vec]
           [org.apache.commons.math3.linear ArrayRealVector RealVector]
           [org.apache.commons.math3.analysis UnivariateFunction]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators #{'abs 'zero?})

;; ## Vector definitions

(declare angle-between)
(declare normalize)
(declare div)

(defn- find-idx-reducer-fn 
  "Helper function for reduce to find index for maximum/minimum value in vector."
  [f]
  (fn [[^long midx ^long curr ^double v] ^double x]
    (if (f x v)
      [curr (inc curr) x]
      [midx (inc curr) v])))

;; Add `VectorProto` to Clojure vector using map/reduce terms.
(extend Seqable
  prot/VectorProto
  {:to-double-array m/seq->double-array 
   :to-acm-vec (fn [v] (ArrayRealVector. (m/seq->double-array v)))
   :to-vec (fn [v] (apply vector-of :double v))
   :as-vec (fn
             ([v xs] (take (count v) (concat xs (repeat 0.0))))
             ([v] (prot/as-vec v nil)))
   :fmap (fn [v f] (map f v))
   :approx (fn
             ([v] (map m/approx v))
             ([v ^long d] (map (fn [^double x] (m/approx x d)) v)))
   :magsq (fn ^double [v] (reduce (fn ^double [^double b ^double x] (+ b (* x x))) 0.0 v))
   :mag (fn ^double [v] (m/sqrt (prot/magsq v)))
   :dot (fn ^double [v1 v2] (reduce m/+ (map m/* v1 v2)))
   :add (fn [v1 v2] (map m/+ v1 v2))
   :sub (fn [v1 v2] (map m/- v1 v2))
   :shift (fn [v1 ^double v] (map (fn [^double x] (m/+ x v)) v1))
   :mult (fn [v1 ^double v] (map (fn [^double x] (m/* x v)) v1))
   :emult (fn [v1 v2] (map m/* v1 v2))
   :abs (fn [v] (map m/abs v))
   :mx (fn ^double [v] (reduce m/max v))
   :mn (fn ^double [v] (reduce m/min v))
   :emx (fn [v1 v2] (map m/max v1 v2))
   :emn (fn [v1 v2] (map m/min v1 v2))
   :maxdim (fn ^long [v] (first (reduce (find-idx-reducer-fn m/>) [0 0 (first v)] v)))
   :mindim (fn ^long [v] (first (reduce (find-idx-reducer-fn m/<) [0 0 (first v)] v)))
   :sum (fn ^double [v] (reduce m/+ 0.0 v))
   :prod (fn ^double [v] (reduce m/* v))
   :size count
   :permute (fn [v1 v2] (map (fn [^long idx] (nth v1 idx)) v2))
   :reciprocal (fn [v] (map (fn [^double v] (m// v)) v))
   :heading (fn ^double [v] (angle-between v (conj (repeat (dec (count v)) 0.0) 1.0)))
   :interpolate (fn [v1 v2 ^double t f] (map (fn [^double x1 ^double x2] (f x1 x2 t)) v1 v2))
   :einterpolate (fn [v1 v2 v f] (map (fn [^double x1 ^double x2 ^double t] (f x1 x2 t)) v1 v2 v))
   :econstrain (fn [v ^double val1 ^double val2] (map (fn [^double x] (m/constrain x val1 val2)) v))
   :is-zero? (fn [v] (every? m/zero? v))})

;; Add `VectorProto` to Clojure vector using mapv/reduce terms.
(extend IPersistentVector
  prot/VectorProto
  {:to-double-array m/seq->double-array
   :to-acm-vec (fn [v] (ArrayRealVector. (m/seq->double-array v)))
   :to-vec (fn [v] (apply vector-of :double v))
   :as-vec (fn
             ([v xs] (vec (take (count v) (concat xs (repeat 0.0)))))
             ([v] (prot/as-vec v nil)))
   :fmap (fn [v f] (mapv f v))
   :approx (fn
             ([v] (mapv m/approx v))
             ([v ^long d] (mapv (fn [^double x] (m/approx x d)) v)))
   :magsq (fn ^double [v] (reduce (fn ^double [^double b ^double x] (+ b (* x x))) 0.0 v))
   :mag (fn ^double [v] (m/sqrt (prot/magsq v)))
   :dot (fn ^double [v1 v2] (reduce m/+ (map m/* v1 v2)))
   :add (fn [v1 v2] (mapv m/+ v1 v2))
   :sub (fn [v1 v2] (mapv m/- v1 v2))
   :shift (fn [v1 ^double v] (mapv (fn [^double x] (m/+ x v)) v1))
   :mult (fn [v1 ^double v] (mapv (fn [^double x] (m/* x v)) v1))
   :emult (fn [v1 v2] (mapv m/* v1 v2))
   :abs (fn [v] (mapv m/abs v))
   :mx (fn ^double [v] (reduce m/max v))
   :mn (fn ^double [v] (reduce m/min v))
   :emx (fn [v1 v2] (mapv m/max v1 v2))
   :emn (fn [v1 v2] (mapv m/min v1 v2))
   :maxdim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/>) [0 0 (first v)] v)))
   :mindim (fn ^long [v] (first (reduce (find-idx-reducer-fn clojure.core/<) [0 0 (first v)] v)))
   :sum (fn ^double [v] (reduce m/+ 0.0 v))
   :prod (fn ^double [v] (reduce m/* v))
   :size count
   :permute (fn [v1 v2] (mapv (fn [^long idx] (v1 idx)) v2))
   :reciprocal (fn [v] (mapv (fn [^double v] (m// v)) v))
   :heading (fn ^double [v] (angle-between v (conj (repeat (dec (count v)) 0.0) 1.0)))
   :interpolate (fn [v1 v2 ^double t f] (mapv (fn [^double x1 ^double x2] (f x1 x2 t)) v1 v2))
   :einterpolate (fn [v1 v2 v f] (mapv (fn [^double x1 ^double x2 ^double t] (f x1 x2 t)) v1 v2 v))
   :econstrain (fn [v ^double val1 ^double val2] (mapv (fn [^double x] (m/constrain x val1 val2)) v))
   :is-zero? (fn [v] (every? m/zero? v))})

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

(extend (Class/forName "[D")
  prot/VectorProto
  {:to-double-array identity
   :to-acm-vec (fn [^doubles arr] (ArrayRealVector. arr))
   :to-vec (fn [^doubles arr] (let [^Vec v (vector-of :double)]
                               (Vec. (.am v) (alength arr) (.shift v) (.root v) arr (.meta v))))
   :as-vec (fn ([^doubles v xs] (double-array (take (alength v) (concat xs (repeat 0.0)))))
             ([v] (prot/as-vec v nil)))
   :fmap (fn [^doubles arr f] (amap arr idx _ret (double (f (Array/aget arr idx)))))
   :approx (fn ([^doubles arr] (amap arr idx _ret (m/approx (Array/aget arr idx))))
             ([^doubles arr ^long d] (amap arr idx _ret (m/approx (Array/aget arr idx) d))))
   :magsq (fn [^doubles arr] (Array/dot arr arr))
   :mag (fn [^doubles arr] (m/sqrt (Array/dot arr arr)))
   :dot (fn [^doubles arr ^doubles v2] (Array/dot arr v2))
   :add (fn [^doubles arr ^doubles v2] (Array/add arr v2))
   :sub (fn [^doubles arr ^doubles v2] (Array/sub arr v2))
   :shift (fn [^doubles arr ^double v] (Array/shift arr v))
   :mult (fn [^doubles arr ^double v] (Array/scale arr v))
   :emult (fn [^doubles arr ^doubles v2] (Array/emult arr v2))
   :abs (fn [^doubles arr] (Array/abs arr))
   :mx (fn [^doubles arr] (Array/max arr))
   :mn (fn [^doubles arr] (Array/min arr))
   :maxdim (fn [^doubles arr] (Array/which_max arr))
   :mindim (fn [^doubles arr] (Array/which_min arr))
   :emx (fn [^doubles arr ^doubles v2] (Array/emax arr v2))
   :emn (fn [^doubles arr ^doubles v2] (Array/emin arr v2))
   :sum (fn [^doubles arr] (Array/sum arr))
   :prod (fn [^doubles arr] (Array/product arr))
   :size alength
   :heading (fn [^doubles arr] (let [v (double-array (alength arr) 0.0)]
                                (Array/aset v 0 1.0)
                                (angle-between arr v)))
   :reciprocal (fn [^doubles arr] (amap arr idx _ret (m// (Array/aget arr idx))))
   :interpolate (fn [^doubles arr ^doubles v2 ^double t f] (amap arr idx _ret
                                                                (double (f (Array/aget arr idx)
                                                                           (Array/aget v2 idx) t))))
   :einterpolate (fn [^doubles arr ^doubles v2 ^doubles v f] (amap arr idx _ret
                                                                  (double (f (Array/aget arr idx)
                                                                             (Array/aget v2 idx)
                                                                             (Array/aget v idx)))))
   :econstrain (fn [^doubles arr ^double val1 ^double val2] (amap arr idx _ret
                                                                 (double (m/constrain (Array/aget arr idx) val1 val2))))
   :is-zero? (fn [^doubles arr] (aevery arr m/zero?))})

(extend ArrayRealVector
  prot/VectorProto
  {:to-double-array (fn [this] (.getDataRef ^ArrayRealVector this))
   :to-acm-vec identity
   :to-vec (fn [v] (apply vector-of :double (.getDataRef ^ArrayRealVector v)))
   :as-vec (fn
             ([^ArrayRealVector v xs] (ArrayRealVector. (double-array (take (.getDimension v)
                                                                            (concat xs (repeat 0.0))))))
             ([v] (prot/as-vec v nil)))
   :fmap (fn [^ArrayRealVector v f] (.map v (reify UnivariateFunction
                                             (value [_ v] (f v)))))
   :approx (fn
             ([v] (prot/fmap v m/approx))
             ([^ArrayRealVector v ^long d] (.map v (reify UnivariateFunction
                                                     (value [_ v] (m/approx v d))))))
   :magsq (fn [^ArrayRealVector v] (.dotProduct v v))
   :mag (fn [^ArrayRealVector v] (.getNorm v))
   :dot (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (.dotProduct v1 v2))
   :add (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (.add v1 v2))
   :sub (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (.subtract v1 v2))
   :shift (fn [^ArrayRealVector v1 ^double v2] (.mapAddToSelf (.copy v1) v2))
   :mult (fn [^ArrayRealVector v1 ^double v2] (.mapMultiplyToSelf (.copy v1) v2))
   :emult (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (.ebeMultiply v1 v2))
   :abs (fn [v] (prot/fmap v m/abs))
   :mx (fn [^ArrayRealVector v1] (.getMaxValue v1))
   :mn (fn [^ArrayRealVector v1] (.getMinValue v1))
   :emx (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (ArrayRealVector. ^doubles (prot/emx (.getDataRef v1)
                                                                                          (.getDataRef v2))))
   :emn (fn [^ArrayRealVector v1 ^ArrayRealVector v2] (ArrayRealVector. ^doubles (prot/emn (.getDataRef v1)
                                                                                          (.getDataRef v2))))
   :maxdim (fn [^ArrayRealVector v1] (.getMaxIndex v1))
   :mindim (fn [^ArrayRealVector v1] (.getMinIndex v1))
   :sum (fn [^ArrayRealVector v1] (prot/sum (.getDataRef v1)))
   :prod (fn [^ArrayRealVector v1] (prot/prod (.getDataRef v1)))
   :size (fn ^long [^ArrayRealVector v1] (.getDimension v1))
   :heading (fn [^ArrayRealVector v1] (prot/heading (.getDataRef v1)))
   :reciprocal (fn [v1] (prot/fmap v1 /))
   :interpolate (fn [^ArrayRealVector v1 ^ArrayRealVector v2 ^double t f]
                  (ArrayRealVector. ^doubles (prot/interpolate (.getDataRef v1) (.getDataRef v2) t f)))
   :einterpolate (fn [^ArrayRealVector v1 ^ArrayRealVector v2 ^ArrayRealVector t f]
                   (ArrayRealVector. ^doubles (prot/einterpolate (.getDataRef v1)
                                                                 (.getDataRef v2)
                                                                 (.getDataRef t) f)))
   :econstrain (fn [^ArrayRealVector v ^double v1 ^double v2]
                 (ArrayRealVector. v ^doubles (prot/econstrain (.getDataRef v) v1 v2)))
   :is-zero? (fn [^ArrayRealVector v] (prot/is-zero? (.getDataRef v)))})

(defn- vec-id-check
  [^long len id]
  (and (number? id) (< (unchecked-int id) len)))

(defn- assert-number
  [n]
  (when-not (number? n) (throw (IllegalArgumentException. "Key must be a number"))))

;; Array Vector
(deftype ArrayVec [^doubles array]
  Object
  (toString [_] (str "#arrayvec " (if (> (alength array) 10)
                                    (str "[" (s/join " " (take 10 array)) "...]")
                                    (vec array))))
  (equals [_ v]
    (and (instance? ArrayVec v)
         (java.util.Arrays/equals array ^doubles (.array ^ArrayVec v))))
  (hashCode [_]
    (mix-collection-hash (java.util.Arrays/hashCode array) (alength array)))
  clojure.lang.IHashEq 
  (hasheq [_]
    (mix-collection-hash (java.util.Arrays/hashCode array) (alength array)))
  Sequential
  Seqable
  (seq [_] (seq array))
  Reversible
  (rseq [_] (reverse array))
  Indexed
  (nth [_ id] (Array/aget array (unchecked-long id)))
  (nth [_ id not-found]
    (let [id (unchecked-long id)]
      (if (< id (alength array)) (Array/aget array id) not-found)))
  ILookup
  (valAt [_ id] (when (vec-id-check (alength array) id) (Array/aget array (unchecked-long id))))
  (valAt [_ id not-found] (if (vec-id-check (alength array) id) (Array/aget array (unchecked-long id)) not-found))
  IReduce
  (reduce [_ f] (reduce f array))
  IReduceInit
  (reduce [_ f start] (reduce f start array))
  Associative
  (containsKey [_ id] (vec-id-check (alength array) id))
  (assoc [v k vl]
    (assert-number k)
    (let [^ArrayVec v v
          arr (aclone ^doubles (.array v))]
      (Array/aset arr (unchecked-long k) (double vl))
      (ArrayVec. arr)))
  (entryAt [v k] (MapEntry. k (v k)))
  IFn
  (invoke [_ n]
    (assert-number n)
    (Array/aget array (unchecked-long n)))
  Counted
  (count [_] (alength array))
  IPersistentVector
  (length [_] (alength array))
  (cons [_ x] (cons x array))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  (empty [_] (double-array 0))
  prot/VectorProto
  (to-double-array [_] array)
  (to-acm-vec [_] (ArrayRealVector. array))
  (to-vec [_] (let [^Vec v (vector-of :double)]
                (Vec. (.am v) (alength array) (.shift v) (.root v) array (.meta v))))
  (as-vec [_ xs] (ArrayVec. (prot/as-vec array xs)))
  (as-vec [v] (prot/as-vec v (repeat 0.0)))
  (fmap [_ f] (ArrayVec. (prot/fmap array f)))
  (approx [_] (ArrayVec. (prot/approx array)))
  (approx [_ d] (ArrayVec. (prot/approx array d)))
  (magsq [_] (Array/dot array array))
  (mag [_] (m/sqrt (Array/dot array array)))
  (dot [_ v2] (Array/dot array ^doubles (.array ^ArrayVec v2)))
  (add [_ v2] (ArrayVec. (prot/add array (.array ^ArrayVec v2))))
  (sub [_ v2] (ArrayVec. (prot/sub array (.array ^ArrayVec v2))))
  (shift [_ v] (ArrayVec. (prot/shift array v)))
  (mult [_ v] (ArrayVec. (prot/mult array v)))
  (emult [_ v2] (ArrayVec. (prot/emult array (.array ^ArrayVec v2))))
  (abs [_] (ArrayVec. (prot/abs array)))
  (mx [_] (Array/max array))
  (mn [_] (Array/min array))
  (maxdim [_] (Array/which_max array))
  (mindim [_] (Array/which_min array))
  (emx [_ v2] (ArrayVec. (prot/emx array (.array ^ArrayVec v2))))
  (emn [_ v2] (ArrayVec. (prot/emn array (.array ^ArrayVec v2))))
  (sum [_] (Array/sum array))
  (prod [_] (Array/product array))
  (size [_] (alength array))
  (heading [_] (prot/heading array))
  (reciprocal [_] (ArrayVec. (prot/reciprocal array)))
  (interpolate [_ v2 t f] (ArrayVec. (prot/interpolate array (.array ^ArrayVec v2) t f))) 
  (einterpolate [_ v2 v f] (ArrayVec. (prot/einterpolate array (.array ^ArrayVec v2) (.array ^ArrayVec v) f)))
  (econstrain [_ val1 val2] (ArrayVec. (prot/econstrain array val1 val2)))
  (is-zero? [_] (aevery array m/zero?)))

(extend-type Number
  prot/VectorProto
  (to-double-array [v] (double-array [v]))
  (to-acm-vec [v] (ArrayRealVector. 1 (double v)))
  (to-vec [v] (vector-of :double (double v)))
  (as-vec
    ([_] 0.0)
    ([_ xs] (double (or (first xs) 0.0))))
  (fmap [v f] (f v))
  (approx
    ([v] (m/approx v))
    ([v d] (m/approx v d)))
  (magsq [v] (m/sq v))
  (mag [v] (m/abs (double v)))
  (dot [v1 v2] (* (double v1) (double v2)))
  (add [v1 v2] (+ (double v1) (double v2)))
  (sub [v1 v2] (- (double v1) (double v2)))
  (shift [v1 v2] (+ (double v1) (double v2)))
  (mult [v1 v] (* (double v1) (double v)))
  (emult [v1 v2] (* (double v1) (double v2)))
  (abs [v] (m/abs (double v)))
  (mx [v] v)
  (mn [v] v)
  (emx [v1 v2] (max (double v1) (double v2)))
  (emn [v1 v2] (min (double v1) (double v2)))
  (maxdim [_] 0)
  (mindim [_] 0)
  (sum [v] v)
  (prod [v] v)
  (size [_] 1)
  (reciprocal [v] (/ (double v)))
  (interpolate [v1 v2 t f] (f v1 v2 t))
  (einterpolate [v1 v2 t f] (f v1 v2 t))
  (econstrain [v val1 val2] (m/constrain (double v) (double val1) (double val2)))
  (is-zero? [v] (m/zero? (double v))))

(defn dhash-code
  "double hashcode"
  (^long [^long state ^double a]
   (let [abits (Double/doubleToLongBits a)
         elt (bit-xor abits (m/>>> abits 32))]
     (+ elt (* 31 state))))
  (^long [^double a]
   (let [abits (Double/doubleToLongBits a)
         elt (bit-xor abits (m/>>> abits 32))]
     (+ elt 31))))

(defn- vec-throw-ioobe
  [^long id len]
  (throw (IndexOutOfBoundsException. (str "Index " id " out of bounds for length " len))))

(declare vec3)

;; Create Vec4 and add all necessary protocols
(deftype Vec4 [^double x ^double y ^double z ^double w]
  Object
  (toString [_] (str "#vec4 [" x ", " y ", " z ", " w "]"))
  (equals [_ v]
    (and (instance? Vec4 v)
         (let [^Vec4 v v]
           (and (== x (.x v))
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
  (nth [_ id] (case (unchecked-int id) 0 x 1 y 2 z 3 w (vec-throw-ioobe id 4)))
  (nth [_ id not-found] (case (unchecked-int id) 0 x 1 y 2 z 3 w not-found))
  ILookup
  (valAt [_ id] (when (vec-id-check 4 id) (case (unchecked-int id) 0 x 1 y 2 z 3 w)))
  (valAt [_ id not-found] (if (number? id) (case (unchecked-int id) 0 x 1 y 2 z 3 w not-found) not-found))
  Associative
  (containsKey [_ id] (boolean (m/<= 0 (unchecked-long id) 3)))
  (assoc [_ k vl]
    (assert-number k)
    (case (unchecked-int k)
      0 (Vec4. vl y z w)
      1 (Vec4. x vl z w)
      2 (Vec4. x y vl w)
      3 (Vec4. x y z vl)
      (vec-throw-ioobe k 2)))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 4)
  IFn
  (invoke [_ id]
    (assert-number id)
    (case (unchecked-int id)
      0 x
      1 y
      2 z
      3 w
      (vec-throw-ioobe id 2)))
  IReduce
  (reduce [_ f] (f (f (f x y) z) w))
  IReduceInit
  (reduce [_ f start] (f (f (f (f start x) y) z) w))
  IPersistentVector
  (length [_] 4)
  (assocN [v i x] (assoc v i x))
  (cons [_ v] [x y z w v])
  IPersistentStack
  (pop [_] (vec3 x y z))
  (peek [_] w)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  (empty [_] [])
  prot/VectorProto
  (to-double-array [v] (double-array v))
  (to-acm-vec [v] (ArrayRealVector. (double-array v)))
  (to-vec [_] (vector-of :double x y z w))
  (as-vec [_ [x y z w]] (Vec4. (or x 0.0) (or y 0.0) (or z 0.0) (or w 0.0)))
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
  (shift [_ v] (let [dv (double v)] (Vec4. (+ x dv) (+ y dv) (+ z dv) (+ w dv))))
  (mult [_ v] (let [dv (double v)] (Vec4. (* x dv) (* y dv) (* z dv) (* w dv))))
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
  (prod [_] (* x y z w))
  (size [_] 4)
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
  (econstrain [_ val1 val2] (let [dval1 (double val1)
                                  dval2 (double val2)]
                              (Vec4. (m/constrain x dval1 dval2)
                                     (m/constrain y dval1 dval2)
                                     (m/constrain z dval1 dval2)
                                     (m/constrain w dval1 dval2))))
  (is-zero? [_] (and (m/zero? x) (m/zero? y) (m/zero? z) (m/zero? w)))
  (heading [v1] (angle-between v1 (Vec4. 1 0 0 0))))

(declare vec2)

;; Create Vec3 and add all necessary protocols
(deftype Vec3 [^double x ^double y ^double z]
  Object
  (toString [_] (str "#vec3 [" x ", " y ", " z "]"))
  (equals [_ v]
    (and (instance? Vec3 v)
         (let [^Vec3 v v]
           (and (== x (.x v))
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
  (nth [_ id] (case (unchecked-int id) 0 x 1 y 2 z (vec-throw-ioobe id 3)))
  (nth [_ id not-found] (case (unchecked-int id) 0 x 1 y 2 z not-found))
  ILookup
  (valAt [_ id] (when (vec-id-check 3 id) (case (unchecked-int id) 0 x 1 y 2 z)))
  (valAt [_ id not-found] (if (number? id) (case (unchecked-int id) 0 x 1 y 2 z not-found) not-found))
  Associative
  (containsKey [_ id] (boolean (m/<= 0 (unchecked-long id) 2)))
  (assoc [_ k vl]
    (assert-number k)
    (case (unchecked-int k)
      0 (Vec3. vl y z)
      1 (Vec3. x vl z)
      2 (Vec3. x y vl)
      (vec-throw-ioobe k 2)))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 3)
  IFn
  (invoke [_ id]
    (assert-number id)
    (case (unchecked-int id)
      0 x
      1 y
      2 z
      (vec-throw-ioobe id 2)))
  IReduce
  (reduce [_ f] (f (f x y) z))
  IReduceInit
  (reduce [_ f start] (f (f (f start x) y) z))
  IPersistentVector
  (length [_] 3)
  (assocN [v i x] (assoc v i x))
  (cons [_ v] (Vec4. x y z (double v)))
  IPersistentStack
  (pop [_] (vec2 x y))
  (peek [_] z)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  (empty [_] [])
  prot/VectorProto
  (to-double-array [v] (double-array v))
  (to-acm-vec [v] (ArrayRealVector. (double-array v)))
  (to-vec [_] (vector-of :double x y z))
  (as-vec [_ [x y z]] (Vec3. (or x 0.0) (or y 0.0) (or z 0.0)))
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
  (shift [_ v] (let [dv (double v)] (Vec3. (+ x dv) (+ y dv) (+ z dv))))
  (mult [_ v] (let [dv (double v)] (Vec3. (* x dv) (* y dv) (* z dv))))
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
  (prod [_] (* x y z))
  (size [_] 3)
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
  (econstrain [_ val1 val2] (let [dval1 (double val1)
                                  dval2 (double val2)]
                              (Vec3. (m/constrain x dval1 dval2)
                                     (m/constrain y dval1 dval2)
                                     (m/constrain z dval1 dval2))))
  (is-zero? [_] (and (m/zero? x) (m/zero? y) (m/zero? z)))
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
          cex (* c e x)
          cfy (* c f y)
          dz (* d z)
          nx (+ (- cex cfy) dz)
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

;; Create Vec2 and add all necessary protocols
(deftype Vec2 [^double x ^double y]
  Object
  (toString [_] (str "#vec2 [" x ", " y "]"))
  (equals [_ v]
    (and (instance? Vec2 v)
         (let [^Vec2 v v]
           (and (== x (.x v))
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
  (nth [_ id] (case (unchecked-int id) 0 x 1 y (vec-throw-ioobe id 2)))
  (nth [_ id not-found] (case (unchecked-int id) 0 x 1 y not-found))
  ILookup
  (valAt [_ id] (when (vec-id-check 2 id) (case (unchecked-int id) 0 x 1 y)))
  (valAt [_ id not-found] (if (number? id) (case (unchecked-int id) 0 x 1 y not-found) not-found))
  Associative
  (containsKey [_ id] (boolean (m/<= 0 (unchecked-long id) 1)))
  (assoc [_ k vl]
    (assert-number k)
    (case (unchecked-int k)
      0 (Vec2. vl y)
      1 (Vec2. x vl)
      (vec-throw-ioobe k 2)))
  (entryAt [v k] (MapEntry. k (v k)))
  Counted
  (count [_] 2)
  IFn
  (invoke [_ id]
    (assert-number id)
    (case (unchecked-int id)
      0 x
      1 y
      (vec-throw-ioobe id 2)))
  IReduce
  (reduce [_ f] (f x y))
  IReduceInit
  (reduce [_ f start] (f (f start x) y))
  IPersistentVector
  (length [_] 2)
  (assocN [v i x] (assoc v i x))
  (cons [_ v] (Vec3. x y (double v)))
  IPersistentStack
  (pop [_] [x])
  (peek [_] y)
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  (empty [_] [])
  prot/VectorProto
  (to-double-array [v] (double-array v))
  (to-acm-vec [v] (ArrayRealVector. (double-array v)))
  (to-vec [_] (vector-of :double x y))
  (as-vec [_ [x y]] (Vec2. (or x 0.0) (or y 0.0)))
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
  (shift [_ v] (let [dv (double v)] (Vec2. (+ x dv) (+ y dv))))
  (mult [_ v] (let [dv (double v)] (Vec2. (* x dv) (* y dv))))
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
  (prod [_] (* x y))
  (size [_] 2)
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
  (econstrain [_ val1 val2] (let [dval1 (double val1)
                                  dval2 (double val2)]
                              (Vec2. (m/constrain x dval1 dval2)
                                     (m/constrain y dval1 dval2))))
  (is-zero? [_] (and (m/zero? x) (m/zero? y)))
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

;;

(def ^{:deprecated "v1.3.0" :doc "Same as [[fmap]]. Deprecated."} applyf prot/fmap)
(def ^{:deprecated "v1.5.0" :doc "Same as [[vec->Vec]]. Deprecated."} to-vec prot/to-vec)

;; protocol methods mapped

(defn vec->array
  "Converts to double array"
  ^doubles [v] (prot/to-double-array v))

(defn vec->seq
  "Converts to sequence (same as seq)"
  [v]
  (cond
    (instance? RealVector v) (seq (prot/to-double-array v))
    (number? v) (list v)
    :else (seq v)))

(defn vec->RealVector
  "Converts to Apache Commons Math RealVector"
  ^RealVector [v]
  (prot/to-acm-vec v))

(defn real-vector
  "Converts to Apache Commons Math RealVector"
  ^RealVector [v]
  (prot/to-acm-vec v))

(defn vec->Vec
  "Converts to Clojure primitive vector `Vec`."
  [v] 
  (prot/to-vec v))

(defn as-vec
  "Creates vector from sequence as given type. If there is no sequence fill with `0.0`."
  ([v] (prot/as-vec v))
  ([v xs] (prot/as-vec v xs)))

(defn size
  "Returns elements count of the vector."
  ^long [v] (prot/size v))

(defn fmap
  "Applies function to all vector values (like map but returns the same type)."
  [v f]
  (prot/fmap v f))

(defn magsq
  "Returns length of the vector squared."
  ^double [v] (prot/magsq v))

(defn mag
  "Returns length of the vector."
  ^double [v] (prot/mag v))

(defn approx
  "Rounds to `d` (default: 2) decimal places"
  ([v] (prot/approx v 2))
  ([v d] (prot/approx v d)))

(defn delta-eq
  "Equality with given absolute (and/or relative) toleance."
  ([v1 v2] (delta-eq v1 v2 1.0e-6))
  ([v1 v2 ^double abs-tol] (m/near-zero? (mag (prot/sub v1 v2)) abs-tol))
  ([v1 v2 ^double abs-tol ^double rel-tol] (m/near-zero? (mag (prot/sub v1 v2)) abs-tol rel-tol)))

(defn edelta-eq
  "Element-wise equality with given absolute (and/or relative) toleance."
  ([v1 v2] (edelta-eq v1 v2 1.0e-6))
  ([v1 v2 ^double abs-tol] (every? identity (map (fn [^double a ^double b]
                                                   (m/delta-eq a b abs-tol)) v1 v2)))
  ([v1 v2 ^double abs-tol ^double rel-tol]
   (every? identity (map (fn [^double a ^double b]
                           (m/delta-eq a b abs-tol rel-tol)) v1 v2))))

(defn dot
  "Dot product of two vectors."
  ^double [v1 v2] (prot/dot v1 v2))

(defn add
  "Sum of two vectors."
  ([v] v)
  ([v1 v2] (prot/add v1 v2)))

(defn sub
  "Subtraction of two vectors."
  ([v] (prot/mult v -1.0))
  ([v1 v2] (prot/sub v1 v2)))

(defn shift
  "Adds a value to every vector element."
  ([v] v)
  ([v x] (prot/shift v x)))

(defn mult
  "Multiplies vector by number `x`."
  ([v] v)
  ([v x] (prot/mult v x)))

(defn emult
  "Element-wise vector multiplication (Hadamard product)."
  ([v1] v1)
  ([v1 v2] (prot/emult v1 v2)))

(defn abs
  "Absolute value of vector elements"
  [v] (prot/abs v))

(defn mx
  "Maximum value of vector elements"
  ^double [v] (prot/mx v))

(defn mn
  "Minimum value of vector elements"
  ^double [v] (prot/mn v))

(defn emx
  "Element-wise max from two vectors."
  [v1 v2] (prot/emx v1 v2))

(defn emn
  "Element-wise min from two vectors."
  [v1 v2] (prot/emn v1 v2))

(defn maxdim
  "Index of maximum value."
  ^long [v] (prot/maxdim v))

(defn mindim
  "Index of minimum value."
  ^long [v] (prot/mindim v))

(defn base-from
  "List of perpendicular vectors (basis). Works only for `Vec2` and `Vec3` types."
  [v] (prot/base-from v))

(defn sum
  "Sum of elements"
  (^double [] 0.0)
  (^double [v] (prot/sum v)))

(defn prod
  "Product of elements"
  (^double [] 1.0)
  (^double [v] (prot/prod v)))

(defn permute
  "Permutes vector elements with given indices."
  [v idxs] (prot/permute v idxs))

(defn reciprocal
  "Reciprocal of elements."
  [v] (prot/reciprocal v))

(defn interpolate 
  "Interpolates vectors, optionally set interpolation fn (default: lerp)"
  ([v1 v2 t] (prot/interpolate v1 v2 t m/lerp))
  ([v1 v2 t f] (prot/interpolate v1 v2 t f)))

(defn lerp
  "Linear interpolation of vectors"
  [v1 v2 t] (prot/interpolate v1 v2 t m/lerp))

(defn einterpolate 
  "Interpolates vector selement-wise, optionally set interpolation fn (default: lerp)"
  ([v1 v2 v] (prot/einterpolate v1 v2 v m/lerp))
  ([v1 v2 v f] (prot/einterpolate v1 v2 v f)))

(defn econstrain
  "Element-wise constrain"
  [v mn mx] (prot/econstrain v mn mx))

(defn is-zero?
  "Is vector zero?"
  [v] (prot/is-zero? v))

(defn zero?
  "Is vector zero?"
  [v] (prot/is-zero? v))

(defn is-near-zero?
  "Equality to zero `0` with given absolute (and/or relative) toleance."
  ([v] (is-near-zero? v 1.0e-6))
  ([v ^double abs-tol] (m/near-zero? (mag v) abs-tol))
  ([v ^double abs-tol ^double rel-tol] (m/near-zero? (mag v) abs-tol rel-tol)))

(defn near-zero?
  "Equality to zero `0` with given absolute (and/or relative) toleance."
  ([v] (near-zero? v 1.0e-6))
  ([v ^double abs-tol] (m/near-zero? (mag v) abs-tol))
  ([v ^double abs-tol ^double rel-tol] (m/near-zero? (mag v) abs-tol rel-tol)))

(defn heading
  "Angle between vector and unit vector `[1,0,...]`"
  ^double [v] (prot/heading v))

(defn cross
  "Cross product"
  [v1 v2] (prot/cross v1 v2))

(defn rotate
  "Rotates vector. Only for `Vec2` and `Vec3` types."
  ([v angle] (prot/rotate v angle))
  ([v angle-x angle-y angle-z] (prot/rotate v angle-x angle-y angle-z)))

(defn axis-rotate
  "Rotates vector. Only for `Vec3` types"
  ([v angle axis] (prot/axis-rotate v angle axis))
  ([v angle axis pivot] (prot/axis-rotate v angle axis pivot)))

(defn perpendicular
  "Perpendicular vector. Only for `Vec2` and `Vec3` types."
  ([v] (prot/perpendicular v))
  ([v1 v2] (prot/perpendicular v1 v2)))

(defn transform
  "Transforms vector.

  Maps point to a coordinate system defined by origin, vx, vy and vz (as bases).

  Only for `Vec2` and `Vec3` types."
  ([v o vx vy] (prot/transform v o vx vy))
  ([v o vx vy vz] (prot/transform v o vx vy vz)))

(defn to-polar
  "To polar coordinates (2d, 3d only), first element is length, the rest angle."
  [v] (prot/to-polar v))

(defn from-polar
  "From polar coordinates (2d, 3d only)"
  [v] (prot/from-polar v))

(defn triple-product
  "a o (b x c)"
  ^double [a b c]
  (dot a (cross b c)))

;; creators

(defn vec2
  "Creates 2d vector."
  ([x y] (Vec2. x y))
  ([[x y]] (Vec2. x y))
  ([] (Vec2. 0.0 0.0)))

(defn vec3
  "Creates Vec2 vector"
  ([x y z] (Vec3. x y z))
  ([^Vec2 v z] (Vec3. (.x v) (.y v) z))
  ([[x y z]] (Vec3. x y z))
  ([] (Vec3. 0.0 0.0 0.0)))

(defn vec4
  "Creates Vec4 vector"
  ([x y z w] (Vec4. x y z w))
  ([^Vec3 v w] (Vec4. (.x v) (.y v) (.z v) w))
  ([^Vec2 v z w] (Vec4. (.x v) (.y v) z w))
  ([[x y z w]] (Vec4. x y z w))
  ([] (Vec4. 0.0 0.0 0.0 0.0)))

(defn array-vec
  "Creates ArrayVec type based on provided sequence `xs`."
  [xs]
  (ArrayVec. (double-array xs)))

(defn make-vector
  "Returns fixed size vector for given number of dimensions.

  Proper type is used."
  ([dims xs] (prot/as-vec (make-vector dims) xs))
  ([^long dims]
   (when (pos? dims)
     (case dims
       1 0.0
       2 (vec2)
       3 (vec3)
       4 (vec4)
       (array-vec dims)))))

;; ## Common vector functions

(defn div
  "Vector division or reciprocal."
  ([v1 ^double v] (prot/mult v1 (/ v)))
  ([v1] (prot/reciprocal v1)))

(defn ediv
  "Element-wise division of two vectors."
  [v1 v2]
  (prot/emult v1 (prot/reciprocal v2)))

(defn zero-count
  "Counts zeros in vector"
  [v]
  (reduce (fn [^long cnt ^double x]
            (if (m/zero? x) (m/inc cnt) cnt)) 0 v))

(defn clamp
  "Clamps the elements of vector `v` element-wise between a minimum (`mn`) and maximum (`mx`) value.

  In the single-arity version `[v]`, elements are clamped between 0.0 and `Double/MAX_VALUE`."
  ([v mn mx] (prot/econstrain v mn mx))
  ([v] (prot/econstrain v 0.0 Double/MAX_VALUE)))

(defn nonzero-count
  "Counts non zero velues in vector"
  [v]
  (reduce (fn [^long cnt ^double x]
            (if (m/zero? x) cnt (m/inc cnt))) 0 v))

(defn average-vectors
  "Average / centroid of vectors. Input: initial vector (optional), list of vectors"
  ([init vs]
   (div (reduce prot/add init vs) (inc (size vs))))
  ([vs] (average-vectors (first vs) (rest vs))))

(defn average
  "Mean or weighted average of the vector"
  (^double [v] (/ (sum v) (size v)))
  (^double [v weights] (/ (dot v weights) (sum weights))))

(defn dist
  "Euclidean distance between vectors"
  ^double [v1 v2]
  (mag (prot/sub v1 v2)))

(defn dist-sq
  "Squared Euclidean distance between vectors"
  ^double [v1 v2]
  (magsq (prot/sub v1 v2)))

(defn dist-abs
  "Manhattan distance between vectors"
  ^double [v1 v2]
  (sum (prot/abs (prot/sub v1 v2))))

(defn dist-cheb
  "Chebyshev distance between 2d vectors"
  ^double [v1 v2]
  (mx (prot/abs (prot/sub v1 v2))))

(defn dist-discrete
  "Computes the discrete distance between two vectors `v1` and `v2`.

  This is the number of positions where the corresponding elements are considered different.
  With an optional absolute tolerance `eps`, elements at index `i` are considered different if `|v1_i - v2_i| > eps`.

  Returns the count of differing elements."
  (^double [v1 v2]
   (sum (prot/fmap (prot/sub v1 v2) (fn [^double v] (if (m/zero? v) 0.0 1.0)))))
  (^double [v1 v2 ^double eps]
   (sum (prot/fmap (prot/abs (prot/sub v1 v2)) (fn [^double v] (if (m/<= v eps) 0.0 1.0))))))

(defn dist-canberra
  "Canberra distance"
  ^double [v1 v2]
  (let [num (prot/abs (prot/sub v1 v2))
        denom (prot/fmap (prot/add (prot/abs v1) (prot/abs v2)) (fn [^double v] (if (m/zero? v) 0.0 (/ v))))]
    (sum (prot/emult num denom))))

(defn dist-emd
  "Earth Mover's Distance"
  ^double [v1 v2]
  (first (reduce (fn [[^double s ^double l] [^double a ^double b]]
                   [(+ s (m/abs l)) (- (+ a l) b)]) [0.0 0.0] (map vector v1 v2))))

(defn dist-ang
  "Angular distance"
  ^double [v1 v2]
  (* (m/acos (/ (dot v1 v2) (* (mag v1) (mag v2)))) m/M_1_PI))

(defn sim-cos
  "Cosine similarity"
  ^double [v1 v2]
  (/ (dot v1 v2) (* (mag v1) (mag v2))))

;; List of distance fn
(def distances {:euclid dist
              :euclid-sq dist-sq
              :abs dist-abs
              :cheb dist-cheb
              :canberra dist-canberra
              :emd dist-emd
              :angular dist-ang
              :discrete dist-discrete})

(defn normalize
  "Returns a new vector with the same direction as `v` but with a magnitude of 1.

  If `v` is a zero vector (magnitude is zero), returns a zero vector of the same type."
  [v]
  (let [m (mag v)]
    (if (m/zero? m)
      (as-vec v)
      (div v m))))

(defn set-mag
  "Sets length of the vector."
  [v ^double len]
  (prot/mult (normalize v) len))

(defn limit
  "Limits length of the vector by given value"
  [v ^double len]
  (if (> (magsq v) (* len len))
    (set-mag v len)
    v))

(defn angle-between
  "Returns the angle between two vectors in radians.

  The angle is always positive and lies between 0 and π.
  Returns 0 if either vector is zero."
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
  "Returns the difference between the heading of `v2` and the heading of `v1`.
  The heading is the angle relative to the positive primary axis ([1,0,...]).
  For 2D vectors, this is the difference between their polar angles.

  Returns value from $-2\\pi$ to $2\\pi$.

  See also [[angle-between]] (absolute angle between vectors) and [[heading]]."
  ^double [v1 v2]
  (- (heading v2) (heading v1)))

(defn aligned?
  "Checks if two vectors are aligned, meaning they point in approximately the same direction.
  This is determined by checking if the angle between them (calculated using [[angle-between]])
  is less than the specified `tolerance`.

  Defaults to an absolute tolerance of `1.0e-6`."
  ([v1 v2 ^double tol]
   (< (angle-between v1 v2) tol))
  ([v1 v2]
   (< (angle-between v1 v2) 1.0e-6)))

(defn faceforward
  "Flips vector `n` if `dot(n, v)` is negative. Returns `n` if `dot(n, v)` is non-negative. Useful for ensuring consistent vector orientation, such as making a normal vector face towards a reference vector `v`."
  [n v]
  (if (neg? (dot n v)) 
    (sub n)
    n))

(defn project
  "Calculates the vector projection of `v1` onto `v2`.
  The projection is a vector along the direction of `v2` that represents the component of `v1` in that direction."
  [v1 v2]
  (mult v2 (/ (dot v1 v2) (magsq v2))))

(defn generate-vec2
  "Generates Vec2 with fn(s)"
  ([f1 f2]
   (Vec2. (f1) (f2)))
  ([f]
   (Vec2. (f) (f))))

(defn generate-vec3
  "Generates Vec3 with fn(s)"
  ([f1 f2 f3]
   (Vec3. (f1) (f2) (f3)))
  ([f]
   (Vec3. (f) (f) (f))))

(defn generate-vec4
  "Generates Vec4 with fn(s)"
  ([f1 f2 f3 f4]
   (Vec4. (f1) (f2) (f3) (f4)))
  ([f]
   (Vec4. (f) (f) (f) (f))))

(defn array->vec2
  "Converts doubles array to Vec2"
  [^doubles arr]
  (Vec2. (Array/aget arr 0) (Array/aget arr 1)))

(defn array->vec3
  "Converts doubles array to Vec3"
  [^doubles arr]
  (Vec3. (Array/aget arr 0) (Array/aget arr 1) (Array/aget arr 2)))

(defn array->vec4
  "Converts doubles array to Vec4"
  [^doubles arr]
  (Vec4. (Array/aget arr 0) (Array/aget arr 1) (Array/aget arr 2) (Array/aget arr 3)))

(defn seq->vec2
  "Converts any seq to Vec2"
  [xs]
  (Vec2. (nth xs 0 0.0) (nth xs 1 0.0)))

(defn seq->vec3
  "Converts any seq to Vec3"
  [xs]
  (Vec3. (nth xs 0 0.0) (nth xs 1 0.0) (nth xs 2 0.0)))

(defn seq->vec4
  "Converts any seq to Vec4"
  [xs]
  (Vec4. (nth xs 0 0.0) (nth xs 1 0.0) (nth xs 2 0.0) (nth xs 3 0.0)))

;;

(defn differences
  "Computes the lagged differences of a vector or sequence.

  This function calculates the difference between elements separated by a fixed distance (`lag`)
  and can apply this operation `diffs` times.

  It supports multiple arities:
  - `[v]`: Computes the 1st difference with lag 1 (default: `diffs=1`, `lag=1`).
  - `[v diffs]`: Computes the `diffs`-th difference with lag 1 (default: `lag=1`).
  - `[v diffs lag]`: Computes the `diffs`-th difference with the specified `lag`.

  Returns a sequence of the computed differences."
  ([v] (differences v 1))
  ([v ^long diffs] (differences v diffs 1))
  ([v ^long diffs ^long lag]
   (loop [v (vec->Vec v)
          d diffs]
     (if (m/zero? d)
       v
       (let [l (count v)]
         (recur (sub (subvec v lag)
                     (subvec v 0 (- l lag)))
                (m/dec d)))))))

;; primitive functions

(defmacro ^:private primitive-ops
  "Generate primitive functions operating on vectors"
  [fns]
  (let [v (symbol "vector")]
    `(do ~@(for [f fns
                 :let [nm (symbol (name f))
                       doc (str "Applies " nm " to vector elements.")]]
             `(defn ~nm ~doc
                [~v]
                (prot/fmap ~v ~f))))))

(primitive-ops [m/sin m/cos m/tan m/asin m/acos m/atan m/sinh m/cosh m/tanh m/asinh m/acosh m/atanh
                m/cot m/sec m/csc m/acot m/asec m/acsc m/coth m/sech m/csch m/acoth m/asech m/acsch
                m/sq m/cb m/safe-sqrt m/sqrt m/cbrt m/exp m/log m/log10 m/log2 m/ln m/log1p m/expm1
                m/log1pexp m/log1mexp m/log1psq m/log1pmx m/logmxp1 m/logexpm1
                m/pow10
                m/radians m/degrees m/sinc m/sigmoid m/logit m/xlogx
                m/floor m/ceil m/round m/rint m/trunc m/frac m/sfrac m/signum m/sgn])

(defn pow
  "Applies power to a vector elements."
  [v ^double exponent]
  (fmap v (fn [^double x] (m/pow x exponent))))

(defn softmax
  "Calculates the softmax function for the given vector `v`.

  The softmax function transforms a vector of arbitrary real values into a probability distribution, where each element is a value between 0 and 1, and all elements sum to 1.

  It has two arities:
  - `[v]`: Computes the standard softmax: `softmax(v)_i = exp(v_i) / sum_j(exp(v_j))`.
  - `[v t]`: Computes the temperature-scaled softmax: `softmax(v)_i = exp(v_i / t) / sum_j(exp(v_j / t))`. The `t` parameter (temperature, defaults to 1) controls the shape of the output distribution; higher temperatures produce softer, more uniform distributions, while lower temperatures produce sharper distributions closer to a one-hot encoding.

  The implementation uses a numerically stable approach (by shifting values based on the vector's maximum) to prevent overflow when exponentiating large numbers.

  Returns a new vector of the same type and dimension as the input."
  ([v]
   (let [nv (exp (shift v (- (mx v))))
         sm (sum nv)]
     (div nv sm)))
  ([v ^double t]
   (let [nv (exp (div (shift v (- (mx v))) t))
         sm (sum nv)]
     (div nv sm))))

(defn logsoftmax
  "Calculates the element-wise natural logarithm of the softmax function for the given vector `v`.

  The standard log-softmax of an element $v_i$ is $\\log(\\operatorname{softmax}(v)_i)$, which is mathematically equivalent to $v_i - \\log(\\sum_j \\exp(v_j))$.
  This function provides a numerically stable implementation of this calculation, particularly useful for avoiding overflow and underflow issues when dealing with large or small exponential values.

  It supports two arities:
  - `[v]`: Computes the standard log-softmax using the numerically stable form.
  - `[v t]`: Computes the temperature-scaled log-softmax. Elements are divided by the temperature `t` *before* applying the log-softmax formula, i.e., $\\frac{v_i}{t} - \\log(\\sum_j \\wxp(\\fac{v_j}{t}))$. The temperature influences the shape of the corresponding softmax output distribution: $t > 1$ softens it, $t < 1$ sharpens it.

  Log-softmax is frequently used in numerical computations, especially in machine learning (e.g., as part of the cross-entropy loss function) for its superior numerical properties.

  Returns a new vector of the same type and dimension as the input."
  ([v]
   (let [shifted (shift v (- (mx v)))
         lsm (- (m/log (sum (exp shifted))))]
     (shift shifted lsm)))
  ([v ^double t]
   (let [shifted (div (shift v (- (mx v))) t)
         lsm (- (m/log (sum (exp shifted))))]
     (shift shifted lsm))))

(defn logsumexp
  "Calculates the LogSumExp of the vector `v`.

  This is the numerically stable computation of `log(sum(exp(x_i)))` for all elements `x_i` in `v`.
  It is often used to prevent overflow when exponentiating large numbers, especially in machine learning and statistics.

  Returns a double value."
  ^double [v]
  (let [m (mx v)]
    (+ m (m/log (sum (exp (shift v (- m))))))))

(defn logmeanexp
  "Calculates the numerically stable log of the mean of the exponential of each element in vector `v`.

  This is equivalent to `log(mean(exp(v_i)))` for all elements `v_i` in `v`.
  It provides a numerically stable way to compute `log(sum(exp(v_i))) - log(count(v))` by shifting values to prevent overflow during exponentiation.

  Often used in machine learning and statistical contexts where dealing with large values inside exponentials is common.

  Returns a double value."
  ^double [v]
  (let [m (mx v)]
    (+ m (m/log (average (exp (shift v (- m))))))))

;;

;; orthogonal and orthonormal polynomial vectors

(defn orthogonal-polynomials
  "Generates a sequence of orthogonal vectors by evaluating orthogonal polynomials at the points specified in the input sequence `xs`.

  The output vectors represent the values of orthogonal polynomials evaluated at the elements of `xs`. Orthogonality is defined with respect to the discrete inner product based on summation over the points in `xs`.

  The sequence starts with the vector for the polynomial of degree 1 and includes vectors up to degree `(count xs) - 1`.

  This function is useful for constructing orthogonal bases for polynomial approximation or regression on discrete data."
  [xs]
  (let [cnt (count xs)
        m (m// (sum xs) cnt)
        p0 (repeat cnt 1.0)
        p1 (shift xs (m/- m))]
    (->> (iterate (fn [[prev curr ^double nprev ^double ncurr]]
                    (let [sigma (m// (sum (emult xs (sq curr))) ncurr)
                          gamma (m// ncurr nprev)
                          np (sub (emult curr (shift xs (- sigma)))
                                  (mult prev gamma))]
                      [curr np ncurr (sum (sq np))]))
                  [p0 p1 cnt (sum (sq p1))])
         (rest) ;; drop zero degree
         (map first)
         (take (m/dec (count xs))))))

(defn orthonormal-polynomials
  "Generates a sequence of orthonormal vectors by evaluating orthogonal polynomials at the points specified in the input sequence `xs` and normalizing the resulting vectors.

  The sequence starts with the vector for the polynomial of degree 1 and includes vectors up to degree `(count xs) - 1`. Orthogonality (and thus orthonormality after normalization) is with respect to the discrete inner product based on summation over the points in `xs`.

  This function produces an orthonormal basis suitable for polynomial approximation or regression on discrete data points defined by `xs`, derived from the orthogonal basis computed by [[orthogonal-polynomials]]."
  [xs]
  (map normalize (orthogonal-polynomials xs)))

;;

(defmethod print-method Vec2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method Vec3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method Vec4 [v ^java.io.Writer w] (.write w (str v)))

(m/unuse-primitive-operators #{'abs 'zero?})

