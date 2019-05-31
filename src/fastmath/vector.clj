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
            [clojure.string :as s])
  (:import [clojure.lang Counted IFn ISeq IPersistentVector IPersistentCollection Seqable Sequential Reversible Indexed ILookup
            Associative MapEntry]
           [clojure.core Vec]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Tolerance (epsilon), used in `is-near-zero?` fn
(def ^:const ^{:doc "Tolerance used in [[is-near-zero?]]. Values less than this value are treated as zero."} ^double TOLERANCE 1.0e-6)

;; ## Vector definitions

(defprotocol VectorProto
  "Vector operations"
  (^{:metadoc/categories #{:gen}} to-vec [v] "Convert to Clojure primitive vector `Vec`.")
  (^{:metadoc/categories #{:gen}} as-vec [v] [v xs] "Create vector from sequence as given type.")
  (^{:metadoc/categories #{:op}} fmap [v f] "Apply function to all vector values (like map but returns the same type).")
  (^{:metadoc/categories #{:op}} approx [v] [v d] "Round to 2 (or `d`) decimal places")
  (^{:metadoc/categories #{:dist :geom}} magsq [v1] "Length of the vector squared.")
  (^{:metadoc/categories #{:dist :geom}} mag [v1] "length of the vector.")
  (^{:metadoc/categories #{:geom}} dot [v1 v2] "Dot product of two vectors.")
  (^{:metadoc/categories #{:op}} add [v1] [v1 v2] "Sum of two vectors.")
  (^{:metadoc/categories #{:op}} sub [v1] [v1 v2] "Subtraction of two vectors.")
  (^{:metadoc/categories #{:op}} mult [v1 v] "Multiply vector by number `v`.")
  (^{:metadoc/categories #{:op}} emult [v1 v2] "Element-wise vector multiplication (Hadamard product).")
  (^{:metadoc/categories #{:mop}} abs [v1] "Absolute value of vector elements")
  (^{:metadoc/categories #{:op}} mx [v1] "Maximum value of vector elements")
  (^{:metadoc/categories #{:op}} mn [v1] "Minimum value of vector elements")
  (^{:metadoc/categories #{:op}} emx [v1 v2] "Element-wise max from two vectors.")
  (^{:metadoc/categories #{:op}} emn [v1 v2] "Element-wise min from two vectors.")
  (^{:metadoc/categories #{:op}} maxdim [v] "Index of maximum value.")
  (^{:metadoc/categories #{:op}} mindim [v] "Index of minimum value.")
  (^{:metadoc/categories #{:geom}} base-from [v] "List of perpendicular vectors (basis)")
  (^{:metadoc/categories #{:op}} sum [v1] "Sum of elements")
  (^{:metadoc/categories #{:op}} permute [v idxs] "Permute vector elements with given indices.")
  (^{:metadoc/categories #{:op}} reciprocal [v] "Reciprocal of elements.")
  (^{:metadoc/categories #{:op}} interpolate [v1 v2 t] [v1 v2 t f] "Interpolate vectors, optionally set interpolation fn")
  (^{:metadoc/categories #{:op}} einterpolate [v1 v2 v] [v1 v2 v f] "Interpolate vectors element-wise, optionally set interpolation fn")
  (^{:metadoc/categories #{:op}} econstrain [v val1 val2] "Element-wise constrain")
  (^{:metadoc/categories #{:op}} is-zero? [v1] "Is vector zero?")
  (^{:metadoc/categories #{:op}} is-near-zero? [v1] "Is vector almost zero? (all absolute values of elements are less than `TOLERANCE`)")
  (^{:metadoc/categories #{:geom}} heading [v1] "Angle between vector and unit vector `[1,0,...]`")
  (^{:metadoc/categories #{:geom}} cross [v1 v2] "Cross product")
  (^{:metadoc/categories #{:geom}} rotate [v1 angle] [v1 anglex angley anglez] "Rotate vector")
  (^{:metadoc/categories #{:geom}} perpendicular [v1] [v1 v2] "Perpendicular vector (only 2d).")
  (^{:metadoc/categories #{:geom}} axis-rotate [v1 angle axis] [v1 angle axis pivot] "Rotate around axis, 3d only")
  (^{:metadoc/categories #{:geom}} transform [v1 o vx vy] [v1 o vx vy vz] "Transform vector; map point to coordinate system defined by origin, vx and vy (as bases), 2d and 3d only.")
  (^{:metadoc/categories #{:geom}} to-polar [v1] "To polar coordinates (2d, 3d only), first element is length, the rest angle.")
  (^{:metadoc/categories #{:geom}} from-polar [v1] "From polar coordinates (2d, 3d only)"))

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
  "Is your value less than TOLERANCE?"
  [^double v]
  (< (m/abs v) TOLERANCE))

;; Add `VectorProto` to Clojure vector using map/reduce terms.
(extend ISeq
  VectorProto
  {:to-vec #(apply conj (vector-of :double) %1)
   :as-vec (fn
             ([v xs] (take (count v) xs))
             ([v] (as-vec v (repeat 0.0))))
   :fmap #(map %2 %1)
   :approx (fn
             ([v] (map m/approx v))
             ([v d] (map #(m/approx ^double % d) v)))
   :magsq (fn [v] (reduce #(+ ^double %1 (* ^double %2 ^double %2)) 0.0 v))
   :mag #(m/sqrt (magsq %))
   :dot #(reduce clojure.core/+ (map clojure.core/* %1 %2))
   :add (fn
          ([v] v)
          ([v1 v2] (map clojure.core/+ v1 v2)))
   :sub (fn
          ([v] (map clojure.core/- v))
          ([v1 v2] (map clojure.core/- v1 v2)))
   :mult (fn [v1 v] (map #(clojure.core/* ^double % ^double v) v1))
   :emult #(map clojure.core/* %1 %2)
   :abs #(map m/abs %)
   :mx #(reduce clojure.core/max %)
   :mn #(reduce clojure.core/min %)
   :emx #(mapv clojure.core/max %1 %2)
   :emn #(mapv clojure.core/min %1 %2)
   :maxdim #(first (reduce (find-idx-reducer-fn clojure.core/>) [0 0 (first %)] %))
   :mindim #(first (reduce (find-idx-reducer-fn clojure.core/<) [0 0 (first %)] %))
   :sum #(reduce clojure.core/+ %)
   :permute #(map (fn [idx] (%1 idx)) %2)
   :reciprocal #(map (fn [^double v] (/ v)) %)
   :heading #(angle-between % (reduce conj [1.0] (repeatedly (dec (count %)) (constantly 0.0))))
   :interpolate (fn
                  ([v1 v2 t f]
                   (map #(f %1 %2 t) v1 v2))
                  ([v1 v2 t] (interpolate v1 v2 t m/lerp)))
   :einterpolate (fn
                   ([v1 v2 v f]
                    (map #(f %1 %2 %3) v1 v2 v))
                   ([v1 v2 v] (einterpolate v1 v2 v m/lerp)))
   :econstrain (fn [v val1 val2] (map #(m/constrain ^double %1 ^double val1 ^double val2) v))
   :is-zero? #(every? clojure.core/zero? %)
   :is-near-zero? #(every? near-zero? %)})

;; Add `VectorProto` to Clojure vector using mapv/reduce terms.
(extend IPersistentVector
  VectorProto
  {:to-vec #(apply conj (vector-of :double) %1)
   :as-vec (fn
             ([v xs] (vec (take (count v) xs)))
             ([v] (as-vec v (repeat 0.0))))
   :fmap #(mapv %2 %1)
   :approx (fn
             ([v] (mapv m/approx v))
             ([v d] (mapv #(m/approx ^double % d) v)))
   :magsq (fn [v] (reduce #(+ ^double %1 (* ^double %2 ^double %2)) 0.0 v))
   :mag #(m/sqrt (magsq %))
   :dot #(reduce clojure.core/+ (map clojure.core/* %1 %2))
   :add (fn
          ([v] v)
          ([v1 v2] (mapv clojure.core/+ v1 v2)))
   :sub (fn
          ([v] (mapv clojure.core/- v))
          ([v1 v2] (mapv clojure.core/- v1 v2)))
   :mult (fn [v1 v] (mapv #(clojure.core/* ^double % ^double v) v1))
   :emult #(mapv clojure.core/* %1 %2)
   :abs #(mapv m/abs %)
   :mx #(reduce clojure.core/max %)
   :mn #(reduce clojure.core/min %)
   :emx #(mapv clojure.core/max %1 %2)
   :emn #(mapv clojure.core/min %1 %2)
   :maxdim #(first (reduce (find-idx-reducer-fn clojure.core/>) [0 0 (first %)] %))
   :mindim #(first (reduce (find-idx-reducer-fn clojure.core/<) [0 0 (first %)] %))
   :sum #(reduce clojure.core/+ %)
   :permute #(mapv (fn [idx] (%1 idx)) %2)
   :reciprocal #(mapv (fn [^double v] (/ v)) %)
   :heading #(angle-between % (reduce conj [1.0] (repeatedly (dec (count %)) (constantly 0.0))))
   :interpolate (fn
                  ([v1 v2 t f]
                   (mapv #(f %1 %2 t) v1 v2))
                  ([v1 v2 t] (interpolate v1 v2 t m/lerp)))
   :einterpolate (fn
                   ([v1 v2 v f]
                    (mapv #(f %1 %2 %3) v1 v2 v))
                   ([v1 v2 v] (einterpolate v1 v2 v m/lerp)))
   :econstrain (fn [v val1 val2] (mapv #(m/constrain ^double %1 ^double val1 ^double val2) v))
   :is-zero? #(every? clojure.core/zero? %)
   :is-near-zero? #(every? near-zero? %)})

(defn- aevery
  "Array version of every"
  [^doubles array pred]
  (let [s (alength array)]
    (loop [idx (unchecked-long 0)]
      (if (< idx s)
        (if (pred (aget array idx))
          (recur (inc idx))
          false)
        true))))

(extend-type (Class/forName "[D")
  VectorProto
  (to-vec [array] (let [^Vec v (vector-of :double)]
                    (Vec. (.am v) (alength ^doubles array) (.shift v) (.root v) array (.meta v))))
  (as-vec
    ([v xs] (double-array (take (alength ^doubles v) xs)))
    ([v] (as-vec v (repeat 0.0))))
  (fmap [array f] (amap ^doubles array idx ret ^double (f (aget ^doubles array idx))))
  (approx
    ([array] (amap ^doubles array idx ret ^double (m/approx (aget ^doubles array idx))))
    ([array d] (amap ^doubles array idx ret ^double (m/approx (aget ^doubles array idx) d))))
  (magsq [array] (smile.math.Math/dot ^doubles array ^doubles array))
  (mag [v1] (m/sqrt (magsq v1)))
  (dot [array v2] (smile.math.Math/dot ^doubles array ^doubles v2))
  (add
    ([v] v)
    ([array v2] (let [b (aclone ^doubles array)]
                  (smile.math.Math/plus b ^doubles v2)
                  b)))
  (sub
    ([array] (amap ^doubles array idx ret (- (aget ^doubles array idx))))
    ([array v2] (let [b (aclone ^doubles array)]
                  (smile.math.Math/minus b ^doubles v2)
                  b)))
  (mult [array v] (let [b (aclone ^doubles array)]
                    (smile.math.Math/scale ^double v ^doubles b)
                    b))
  (emult [array v2] (amap ^doubles array idx ret (* (aget ^doubles array idx) ^double (v2 idx))))
  (abs [array] (amap ^doubles array idx ret (m/abs (aget ^doubles array idx))))
  (mx [array] (smile.math.Math/max ^doubles array))
  (mn [array] (smile.math.Math/min ^doubles array))
  (maxdim [array] (smile.math.Math/whichMax ^doubles array))
  (mindim [array] (smile.math.Math/whichMin ^doubles array))
  (emx [array v2] (amap ^doubles array idx ret (max (aget ^doubles array idx) ^double (v2 idx))))
  (emn [array v2] (amap ^doubles array idx ret (min (aget ^doubles array idx) ^double (v2 idx))))
  (sum [array] (smile.math.Math/sum ^doubles array))
  (heading [array] (let [v (double-array (alength ^doubles array) 0.0)]
                     (aset v 0 1.0)
                     (angle-between array v)))
  (reciprocal [array] (amap ^doubles array idx ret (/ (aget ^doubles array idx))))
  (interpolate
    ([v1 v2 t] (interpolate v1 v2 t m/lerp))
    ([array v2 t f] (amap ^doubles array idx ret ^double (f (aget ^doubles array idx) (v2 idx) t))))
  (einterpolate
    ([v1 v2 v] (einterpolate v1 v2 v m/lerp))
    ([array v2 v f] (amap ^doubles array idx ret ^double (f (aget ^doubles array idx) (v2 idx) (v idx)))))
  (econstrain [array val1 val2] (amap ^doubles array idx ret ^double (m/constrain ^double (aget ^doubles array idx) ^double val1 ^double val2)))
  (is-zero? [array] (aevery array #(zero? ^double %)))
  (is-near-zero? [array] (aevery array near-zero?)))

;; Array Vector
(deftype ArrayVec [^doubles array]
  Object
  (toString [_] (str "#arrayvec " (if (> (alength array) 10)
                                    (str "[" (s/join " " (take 10 array)) "...]")
                                    (vec array))))
  (equals [_ v]
    (bool-and (instance? ArrayVec v)
              (smile.math.Math/equals array ^doubles (.array ^ArrayVec v) m/MACHINE-EPSILON)))
  (hashCode [_]
    (java.util.Arrays/hashCode array))
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
  VectorProto
  (to-vec [_] (let [^Vec v (vector-of :double)]
                (Vec. (.am v) (alength array) (.shift v) (.root v) array (.meta v))))
  (as-vec [v xs] (ArrayVec. (as-vec array xs)))
  (as-vec [v] (as-vec v (repeat 0.0)))
  (fmap [_ f] (ArrayVec. (fmap array f)))
  (approx [_] (ArrayVec. (approx array)))
  (approx [_ d] (ArrayVec. (approx array d)))
  (magsq [_] (smile.math.Math/dot array array))
  (mag [v1] (m/sqrt (magsq v1)))
  (dot [_ v2] (smile.math.Math/dot array ^doubles (.array ^ArrayVec v2)))
  (add [v] v)
  (add [_ v2] (ArrayVec. (add array (.array ^ArrayVec v2))))
  (sub [_] (ArrayVec. (sub array)))
  (sub [_ v2] (ArrayVec. (sub array (.array ^ArrayVec v2))))
  (mult [_ v] (ArrayVec. (mult array v)))
  (emult [_ v2] (ArrayVec. (emult array v2)))
  (abs [_] (ArrayVec. (abs array)))
  (mx [_] (smile.math.Math/max array))
  (mn [_] (smile.math.Math/min array))
  (maxdim [_] (smile.math.Math/whichMax array))
  (mindim [_] (smile.math.Math/whichMin array))
  (emx [_ v2] (ArrayVec. (emx array v2)))
  (emn [_ v2] (ArrayVec. (emn array v2)))
  (sum [_] (smile.math.Math/sum array))
  (heading [v1] (heading array))
  (reciprocal [_] (ArrayVec. (reciprocal array)))
  (interpolate [v1 v2 t] (interpolate v1 v2 t m/lerp))
  (interpolate [_ v2 t f] (ArrayVec. (interpolate array v2 t f))) 
  (einterpolate [v1 v2 v] (einterpolate v1 v2 v m/lerp))
  (einterpolate [_ v2 v f] (ArrayVec. (einterpolate array v2 v f)))
  (econstrain [_ val1 val2] (ArrayVec. (econstrain array val1 val2)))
  (is-zero? [_] (aevery array #(zero? ^double %)))
  (is-near-zero? [_] (aevery array near-zero?)))

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
  (hashCode [_]
    (unchecked-int (dhash-code (dhash-code x) y)))
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
  VectorProto
  (to-vec [_] (vector-of :double x y))
  (as-vec [_ [x y]] (Vec2. x y))
  (as-vec [_] (Vec2. 0.0 0.0))
  (fmap [_ f] (Vec2. (f x) (f y)))
  (approx [_] (Vec2. (m/approx x) (m/approx y)))
  (approx [_ d] (Vec2. (m/approx x d) (m/approx y d)))
  (magsq [_] (+ (* x x) (* y y)))
  (mag [_] (m/hypot x y))
  (dot [_ v2] 
    (let [^Vec2 v2 v2] (+ (* x (.x v2)) (* y (.y v2)))))
  (add [v] v)
  (add [_ v2] 
    (let [^Vec2 v2 v2] (Vec2. (+ x (.x v2)) (+ y (.y v2)))))
  (sub [_] (Vec2. (- x) (- y)))
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
    [v (perpendicular v)])
  (sum [_] (+ x y))
  (permute [p [^long i1 ^long i2]]
    (Vec2. (p i1) (p i2)))
  (reciprocal [_] (Vec2. (/ x) (/ y)))
  (interpolate [_ v2 t f]
    (let [^Vec2 v2 v2] (Vec2. (f x (.x v2) t)
                              (f y (.y v2) t))))
  (interpolate [v1 v2 t] (interpolate v1 v2 t m/lerp))
  (einterpolate [_ v2 v f]
    (let [^Vec2 v2 v2
          ^Vec2 v v]
      (Vec2. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v)))))
  (einterpolate [v1 v2 v] (einterpolate v1 v2 v m/lerp))
  (econstrain [_ val1 val2] (Vec2. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y)))
  (is-near-zero? [_] (m/bool-and (near-zero? x) (near-zero? y)))
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
    (Vec2. (mag v) (heading v)))
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
  (hashCode [_]
    (unchecked-int (dhash-code (dhash-code (dhash-code x) y) z)))
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
  VectorProto
  (to-vec [_] (vector-of :double x y z))
  (as-vec [_ [x y z]] (Vec3. x y z))
  (as-vec [_] (Vec3. 0.0 0.0 0.0))
  (fmap [_ f] (Vec3. (f x) (f y) (f z)))
  (approx [_] (Vec3. (m/approx x) (m/approx y) (m/approx z)))
  (approx [_ d] (Vec3. (m/approx x d) (m/approx y d) (m/approx z d)))
  (magsq [_] (+ (* x x) (* y y) (* z z)))
  (mag [_] (m/hypot x y z))
  (dot [_ v2]
    (let [^Vec3 v2 v2] (+ (* x (.x v2)) (* y (.y v2)) (* z (.z v2)))))
  (add [v] v)
  (add [_ v2] 
    (let [^Vec3 v2 v2] (Vec3. (+ x (.x v2)) (+ y (.y v2)) (+ z (.z v2)))))
  (sub [_] (Vec3. (- x) (- y) (- z)))
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
               (div (Vec3. (- z) 0.0 x) (m/hypot x z))
               (div (Vec3. 0.0 z (- y)) (m/hypot y z)))]
      [v v2 (cross v v2)]))
  (sum [_] (+ x y z))
  (permute [p [^long i1 ^long i2 ^long i3]]
    (Vec3. (p i1) (p i2) (p i3)))
  (reciprocal [_] (Vec3. (/ x) (/ y) (/ z)))
  (interpolate [_ v2 t f]
    (let [^Vec3 v2 v2] (Vec3. (f x (.x v2) t)
                              (f y (.y v2) t)
                              (f z (.z v2) t))))
  (interpolate [v1 v2 t] (interpolate v1 v2 t m/lerp))
  (einterpolate [_ v2 v f]
    (let [^Vec3 v2 v2
          ^Vec3 v v]
      (Vec3. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v))
             (f z (.z v2) (.z v)))))
  (einterpolate [v1 v2 v] (einterpolate v1 v2 v m/lerp))
  (econstrain [_ val1 val2] (Vec3. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)
                                   (m/constrain z ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y) (zero? z)))
  (is-near-zero? [_] (bool-and (near-zero? x) (near-zero? y) (near-zero? z)))
  (heading [v1] (angle-between v1 (Vec3. 1 0 0)))
  (cross [_ v2]
    (let [^Vec3 v2 v2
          cx (- (* y (.z v2)) (* (.y v2) z))
          cy (- (* z (.x v2)) (* (.z v2) x))
          cz (- (* x (.y v2)) (* (.x v2) y))]
      (Vec3. cx cy cz)))
  (perpendicular [v1 v2]
    (normalize (cross v1 v2)))
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
          ^Vec3 sa (mult ax (m/sin angle))
          sax (.x sa)
          say (.y sa)
          saz (.z sa)
          ^Vec3 cb (mult ax (- 1.0 cosa))
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
    (add (axis-rotate (sub v1 pivot) angle axis) pivot))
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
    (let [^double r (mag v1)
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
    (unchecked-int (dhash-code (dhash-code (dhash-code (dhash-code x) y) z) w)))
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
  VectorProto
  (to-vec [_] (vector-of :double x y z w))
  (as-vec [_ [x y z w]] (Vec4. x y z w))
  (as-vec [_] (Vec4. 0.0 0.0 0.0 0.0))
  (fmap [_ f] (Vec4. (f x) (f y) (f z) (f w)))
  (approx [_] (Vec4. (m/approx x) (m/approx y) (m/approx z) (m/approx w)))
  (approx [_ d] (Vec4. (m/approx x d) (m/approx y d) (m/approx z d) (m/approx w d)))
  (magsq [_] (+ (* x x) (* y y) (* z z) (* w w)))
  (mag [v1] (m/sqrt (magsq v1)))
  (dot [_ v2]
    (let [^Vec4 v2 v2] (+ (* x (.x v2)) (* y (.y v2)) (* z (.z v2)) (* w (.w v2)))))
  (add [v] v)
  (add [_ v2]
    (let [^Vec4 v2 v2] (Vec4. (+ x (.x v2)) (+ y (.y v2)) (+ z (.z v2)) (+ w (.w v2)))))
  (sub [_] (Vec4. (- x) (- y) (- z) (- w)))
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
  (interpolate [v1 v2 t] (interpolate v1 v2 t m/lerp))
  (einterpolate [_ v2 v f]
    (let [^Vec4 v2 v2
          ^Vec4 v v]
      (Vec4. (f x (.x v2) (.x v))
             (f y (.y v2) (.y v))
             (f z (.z v2) (.z v))
             (f w (.w v2) (.w v)))))
  (einterpolate [v1 v2 v] (einterpolate v1 v2 v m/lerp))
  (econstrain [_ val1 val2] (Vec4. (m/constrain x ^double val1 ^double val2)
                                   (m/constrain y ^double val1 ^double val2)
                                   (m/constrain z ^double val1 ^double val2)
                                   (m/constrain w ^double val1 ^double val2)))
  (is-zero? [_] (bool-and (zero? x) (zero? y) (zero? z) (zero? w)))
  (is-near-zero? [_] (bool-and (near-zero? x) (near-zero? y) (near-zero? z) (near-zero? w)))
  (heading [v1] (angle-between v1 (Vec4. 1 0 0 0))))

;;

(def ^{:deprecated "v1.3.0" :metadoc/categories #{:op} :doc "Same as [[fmap]]. Deprecated."} applyf fmap)

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
  ([dims xs] (as-vec (make-vector dims) xs))
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
  ([v1 ^double v] (mult v1 (/ v)))
  ([v1] (reciprocal v1)))

(defn ediv
  "Element-wise division of two vectors."
  {:metadoc/categories #{:op}} 
  [v1 v2]
  (emult v1 (reciprocal v2)))

(defn zero-count
  "Count zeros in vector"
  {:metadoc/categories #{:op}}
  [v]
  (count (filter #(zero? ^double %) v)))

(defn clamp
  "Clamp elements."
  {:metadoc/categories #{:op}}
  ([v ^double mn ^double mx] (econstrain v mn mx))
  ([v] (econstrain v 0 Double/MAX_VALUE)))

(defn nonzero-count
  "Count non zero velues in vector"
  {:metadoc/categories #{:op}}
  [v]
  (count (remove #(zero? ^double %) v)))

(defn average-vectors
  "Average / centroid of vectors. Input: initial vector (optional), list of vectors"
  {:metadoc/categories #{:op}} 
  ([init vs]
   (div (reduce add init vs) (inc (count vs))))
  ([vs] (average-vectors (first vs) (rest vs))))

(defn dist
  "Euclidean distance between vectors"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (mag (sub v1 v2)))

(defn dist-sq
  "Squared Euclidean distance between vectors"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (magsq (sub v1 v2)))

(defn dist-abs
  "Manhattan distance between vectors"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (sum (abs (sub v1 v2))))

(defn dist-cheb
  "Chebyshev distance between 2d vectors"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (mx (abs (sub v1 v2))))

(defn dist-discrete
  "Discrete distance between 2d vectors"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (sum (fmap (sub v1 v2) #(if (zero? ^double %) 0.0 1.0))))

(defn dist-canberra
  "Canberra distance"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (let [num (abs (sub v1 v2))
        denom (fmap (add (abs v1) (abs v2)) #(if (zero? ^double %) 0.0 (/ ^double %)))]
    (sum (emult num denom))))

(defn dist-emd
  "Earth Mover's Distance"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (first (reduce #(let [[^double s ^double l] %1
                        [^double a ^double b] %2
                        n (- (+ a l) b)]
                    [(+ s (m/abs l)) n]) [0.0 0.0] (map vector v1 v2))))

(defn dist-cos
  "Cosine distance"
  {:metadoc/categories #{:dist}} 
  [v1 v2]
  (- 1.0 (/ ^double (dot v1 v2) (* ^double (mag v1) ^double (mag v2)))))

;; List of distance fn
(def ^{:metadoc/categories #{:dist}}
  distances {:euclid dist
             :euclid-sq dist-sq
             :abs dist-abs
             :cheb dist-cheb
             :canberra dist-canberra
             :emd dist-emd
             :cosine dist-cos
             :discrete dist-discrete})

(defn normalize
  "Normalize vector (set length = 1.0)"
  {:metadoc/categories #{:dist}} 
  [v]
  (let [^double m (mag v)]
    (if (zero? m)
      (Vec2. 0.0 0.0)
      (div v m))))

(defn set-mag
  "Set length of the vector"
  {:metadoc/categories #{:dist}} 
  [v len]
  (mult (normalize v) len))

(defn limit
  "Limit length of the vector by given value"
  {:metadoc/categories #{:dist}} 
  [v ^double len]
  (if (> ^double (magsq v) (* len len))
    (set-mag v len)
    v))

(defn angle-between
  "Angle between two vectors

  See also [[relative-angle-between]]."
  {:metadoc/categories #{:geom}} 
  ^double [v1 v2]
  (if (or (is-zero? v1) (is-zero? v2))
    0
    (let [^double d (dot v1 v2)
          amt (/ d (* ^double (mag v1) ^double (mag v2)))]
      (cond
        (<= amt -1.0) m/PI
        (>= amt 1.0) 0
        :else (m/acos amt)))))

(defn relative-angle-between
  "Angle between two vectors relative to each other.

  See also [[angle-between]]."
  {:metadoc/categories #{:geom}} 
  ^double [v1 v2]
  (- ^double (heading v2) ^double (heading v1)))

(defn aligned?
  "Are vectors aligned (have the same direction)?"
  {:metadoc/categories #{:geom}} 
  [v1 v2]
  (< (angle-between v1 v2) TOLERANCE))

(defn faceforward
  "Flip normal `n` to match the same direction as `v`."
  {:metadoc/categories #{:geom}} 
  [n v]
  (if (neg? ^double (dot n v)) 
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
  [^doubles a]
  (Vec2. (aget a 0) (aget a 1)))

(defn array->vec3
  "Doubles array to Vec3"
  {:metadoc/categories #{:gen}} 
  [^doubles a]
  (Vec3. (aget a 0) (aget a 1) (aget a 2)))

(defn array->vec4
  "Doubles array to Vec4"
  {:metadoc/categories #{:gen}} 
  [^doubles a]
  (Vec4. (aget a 0) (aget a 1) (aget a 2) (aget a 3)))

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
                (fmap ~v ~wfn))))))

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
