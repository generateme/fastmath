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
      * ArrayVec - fixed size vector, n-dimensional, creator [[arrayvec]]
  * Variable size:
      * Clojure's PersistentVector, creator `[]`.

  [[VectorProto]] defines most of the functions.

  Vectors implements also:

  * `Sequable`
  * `Sequencial`
  * `IFn`
  * `Counted`
  * `equals` and `toString` from `Object`"
  {:metadoc/categories {:gen "Creators"
                        :geom "Geometric"
                        :dist "Distance / length"
                        :op "Operations"}}
  (:require [fastmath.core :as m]
            [metadoc.examples :refer :all] 
            [clojure.string :as s])
  (:import [clojure.lang Counted IFn PersistentVector Seqable Sequential]
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
  (^{:metadoc/categories #{:op}} applyf [v f] "Apply function to all vector values (like map).")
  (^{:metadoc/categories #{:op}} approx [v] [v d] "Round to 2 (or `d`) decimal places")
  (^{:metadoc/categories #{:dist :geom}} magsq [v1] "Length of the vector squared.")
  (^{:metadoc/categories #{:dist :geom}} mag [v1] "length of the vector.")
  (^{:metadoc/categories #{:geom}} dot [v1 v2] "Dot product of two vectors.")
  (^{:metadoc/categories #{:op}} add [v1] [v1 v2] "Sum of two vectors.")
  (^{:metadoc/categories #{:op}} sub [v1] [v1 v2] "Subtraction of two vectors.")
  (^{:metadoc/categories #{:op}} mult [v1 v] "Multiply vector by number `v`.")
  (^{:metadoc/categories #{:op}} emult [v1 v] "Element-wise vector multiplication (Hadamard product).")
  (^{:metadoc/categories #{:op}} div [v1 v] "Divide vector by number `v`")
  (^{:metadoc/categories #{:op}} abs [v1] "Absolute value of vector elements")
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
(extend PersistentVector
  VectorProto
  {:to-vec #(apply conj (vector-of :double) %1)
   :applyf #(mapv %2 %1)
   :approx (fn
             ([v] (map m/approx v))
             ([v d] (map #(m/approx ^double % d) v)))
   :magsq (fn [v] (reduce #(+ ^double %1 (* ^double %2 ^double %2)) (double 0) v))
   :mag #(m/sqrt (magsq %))
   :dot #(reduce clojure.core/+ (map clojure.core/* %1 %2))
   :add (fn
          ([v] v)
          ([v1 v2] (mapv clojure.core/+ v1 v2)))
   :sub (fn
          ([v] (mapv clojure.core/- v))
          ([v1 v2] (mapv clojure.core/- v1 v2)))
   :mult (fn [v1 v] (mapv #(clojure.core/* (double %) ^double v) v1))
   :emult #(mapv clojure.core/* %1 %2)
   :div #(mult %1 (/ (double %2)))
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

;; Array Vector
(deftype ArrayVec [^doubles array]
  Object
  (toString [_] (str "#arrayvec " (if (> (alength array) 10)
                                    (str "[" (s/join " " (take 10 array)) "...]")
                                    (vec array))))
  (equals [_ v]
    (smile.math.Math/equals array ^doubles (.array ^ArrayVec v) m/MACHINE-EPSILON))
  Sequential
  Seqable
  (seq [_] (seq array))
  IFn
  (invoke [_ n]
    (aget array ^long n))
  Counted
  (count [_] (alength array))
  VectorProto
  (to-vec [_] (let [^Vec v (vector-of :double)]
                (Vec. (.am v) (alength array) (.shift v) (.root v) array (.meta v))))
  (applyf [_ f] (ArrayVec. (amap array idx ret ^double (f (aget array idx)))))
  (approx [_] (ArrayVec. (amap array idx ret ^double (m/approx (aget array idx)))))
  (approx [_ d] (ArrayVec. (amap array idx ret ^double (m/approx (aget array idx) d))))
  (magsq [_] (smile.math.Math/dot array array))
  (mag [v1] (m/sqrt (magsq v1)))
  (dot [_ v2] (smile.math.Math/dot array ^doubles (.array ^ArrayVec v2)))
  (add [v] v)
  (add [_ v2] (let [b (double-array array)]
                (smile.math.Math/plus b ^doubles (.array ^ArrayVec v2))
                (ArrayVec. b)))
  (sub [_] (ArrayVec. (amap array idx ret (- (aget array idx)))))
  (sub [_ v2] (let [b (double-array array)]
                (smile.math.Math/minus b ^doubles (.array ^ArrayVec v2))
                (ArrayVec. b)))
  (mult [_ v] (let [b (double-array array)]
                (smile.math.Math/scale ^double v b)
                (ArrayVec. b)))
  (emult [_ v2] (ArrayVec. (amap array idx ret (* (aget array idx) ^double (v2 idx)))))
  (div [av v] (mult av (/ ^double v)))
  (abs [_] (ArrayVec. (amap array idx ret (m/abs (aget array idx)))))
  (mx [_] (smile.math.Math/max array))
  (mn [_] (smile.math.Math/min array))
  (maxdim [_] (smile.math.Math/whichMax array))
  (mindim [_] (smile.math.Math/whichMin array))
  (emx [_ v2] (ArrayVec. (amap array idx ret (max (aget array idx) ^double (v2 idx)))))
  (emn [_ v2] (ArrayVec. (amap array idx ret (min (aget array idx) ^double (v2 idx)))))
  (sum [_] (smile.math.Math/sum array))
  (heading [v1] (let [v (double-array (alength array) 0.0)]
                  (aset v 0 1.0)
                  (angle-between v1 (ArrayVec. v))))
  (reciprocal [_] (ArrayVec. (amap array idx ret (/ (aget array idx)))))
  (interpolate [v1 v2 t]
    (interpolate v1 v2 t m/lerp))
  (interpolate [_ v2 t f]
    (ArrayVec. (amap array idx ret ^double (f (aget array idx) (v2 idx) t)))) 
  (einterpolate [v1 v2 v]
    (einterpolate v1 v2 v m/lerp))
  (einterpolate [_ v2 v f]
    (ArrayVec. (amap array idx ret ^double (f (aget array idx) (v2 idx) (v idx)))))
  (econstrain [_ val1 val2] (ArrayVec. (amap array idx ret ^double (m/constrain ^double (aget array idx) ^double val1 ^double val2))))
  (is-zero? [_] (aevery array #(zero? ^double %)))
  (is-near-zero? [_] (aevery array near-zero?)))

;; Create Vec2 and add all necessary protocols
(deftype Vec2 [^double x ^double y]
  Object
  (toString [_] (str "#vec2 [" x ", " y "]"))
  (equals [_ v]
    (and (instance? Vec2 v)
         (let [^Vec2 v v]
           (bool-and (== x (.x v))
                     (== y (.y v))))))
  Sequential
  Seqable
  (seq [_] (list x y))
  Counted
  (count [_] 2)
  IFn
  (invoke [_ id]
    (case (unchecked-int id)
      0 x
      1 y
      nil))
  VectorProto
  (to-vec [_] (vector-of :double x y))
  (applyf [_ f] (Vec2. (f x) (f y)))
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
  (div [_ v] 
    (let [v1 (/ 1.0 ^double v)] (Vec2. (* x v1) (* y v1))))
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
    (and (instance? Vec3 v)
         (let [^Vec3 v v]
           (bool-and (== x (.x v))
                     (== y (.y v))
                     (== z (.z v))))))
  Sequential
  Seqable
  (seq [_] (list x y z))
  Counted
  (count [_] 3)
  IFn
  (invoke [_ id]
    (case (unchecked-int id)
      0 x
      1 y
      2 z
      nil))
  VectorProto
  (to-vec [_] (vector-of :double x y z))
  (applyf [_ f] (Vec3. (f x) (f y) (f z)))
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
  (div [_ v] 
    (let [v1 (/ 1.0 ^double v)] (Vec3. (* x v1) (*  y v1) (* z v1))))
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
  Sequential
  Seqable
  (seq [_] (list x y z w))
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
  VectorProto
  (to-vec [_] (vector-of :double x y z w))
  (applyf [_ f] (Vec4. (f x) (f y) (f z) (f w)))
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
  (div [_ v]
    (let [v1 (/ 1.0 ^double v)] (Vec4. (* x v1) (* y v1) (* z v1) (* w v1))))
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

;; creators

(defn vec2
  "Make 2d vector"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example "Usage" (vec2 0.5 -0.5))]} 
  [x y] (Vec2. x y))

(defn vec3
  "Make Vec2 vector"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage"
                        (vec3 0.5 -0.5 1.0)
                        (let [v (vec2 1 2)]
                          (vec3 v -1.0)))]} 
  ([x y z] (Vec3. x y z))
  ([^Vec2 v z] (Vec3. (.x v) (.y v) z)))

(defn vec4
  "Make Vec4 vector"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage"
                        (vec4 0.5 -0.5 1.0 -1.0)
                        (let [v (vec2 1 2)]
                          (vec4 v -1.0 0.1))
                        (let [v (vec3 0 1 2)]
                          (vec4 v 0.1)))]} 
  ([x y z w] (Vec4. x y z w))
  ([^Vec3 v w] (Vec4. (.x v) (.y v) (.z v) w))
  ([^Vec2 v z w] (Vec4. (.x v) (.y v) z w)))

(defn array-vec
  "Make ArrayVec type based on provided sequence `xs`."
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage"
                        (array-vec [1 2 3 4 5 6 7])
                        (array-vec (range 0.0 1.0 0.25)))
                      (example-session "Operations"
                        (nth (array-vec [9 8 7 6]) 2)
                        (count (array-vec (range 0.1 1.0 0.05)))
                        (seq (array-vec [1 2])))]} 
  [xs]
  (ArrayVec. (double-array xs)))

;; ## Common vector functions

(defn ediv
  "Element-wise division of two vectors."
  {:metadoc/categories #{:op}
   :metadoc/examples [(example "Usage" (ediv (vec2 1 3) (vec2 2 4)))]} 
  [v1 v2]
  (emult v1 (reciprocal v2)))

(defn average-vectors
  "Average / centroid of vectors. Input: initial vector (optional), list of vectors"
  {:metadoc/categories #{:op}
   :metadoc/examples [(example-session "Usage"
                        (average-vectors [[1 2] [0 1] [3 4] [1 2] [4 -1]])
                        (average-vectors (vec2 0 0) [(vec2 1 1) (vec2 1 1) (vec2 1 1)]))]} 
  ([init vs]
   (div (reduce add init vs) (inc (count vs))))
  ([vs] (average-vectors (first vs) (rest vs))))

(defn dist
  "Euclidean distance between vectors"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (mag (sub v1 v2)))

(defn dist-sq
  "Squared Euclidean distance between vectors"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-sq (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-sq [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (magsq (sub v1 v2)))

(defn dist-abs
  "Manhattan distance between vectors"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-abs (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-abs [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (sum (abs (sub v1 v2))))

(defn dist-cheb
  "Chebyshev distance between 2d vectors"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-cheb (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-cheb [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (mx (abs (sub v1 v2))))

(defn dist-discrete
  "Discrete distance between 2d vectors"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-discrete (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-discrete [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (sum (applyf (sub v1 v2) #(if (zero? ^double %) 0.0 1.0))))

(defn dist-canberra
  "Canberra distance"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-canberra (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-canberra [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (let [num (abs (sub v1 v2))
        denom (applyf (add (abs v1) (abs v2)) #(if (zero? ^double %) 0.0 (/ ^double %)))]
    (sum (emult num denom))))

(defn dist-emd
  "Earth Mover's Distance"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-emd (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-emd [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (first (reduce #(let [[^double s ^double l] %1
                        [^double a ^double b] %2
                        n (- (+ a l) b)]
                    [(+ s (m/abs l)) n]) [0.0 0.0] (map vector v1 v2))))

(defn dist-cos
  "Cosine distance"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (dist-cos (vec4 1.0 2.0 -1.0 -2.0) (vec4 1.0 -2.0 3.0 4.0))
                        (dist-cos [9 8 7 6 5 4 3 2 1] [9 8 7 6 5 5 5 5 5]))]} 
  [v1 v2]
  (- 1.0 (/ ^double (dot v1 v2) (* ^double (mag v1) ^double (mag v2)))))

;; List of distance fn
(def ^{:metadoc/categories #{:dist}
       :metadoc/examples [(example "List of distances" (sort (keys distances)))]}
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
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example "Usage" (normalize (vec2 1.0 -1.0)))]} 
  [v]
  (let [^double m (mag v)]
    (if (zero? m)
      (Vec2. 0.0 0.0)
      (div v m))))

(defn set-mag
  "Set length of the vector"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (set-mag (vec2 0.22 0.22) (m/sqrt 2.0))
                        (set-mag (vec2 1.0 1.0) 0.0))]} 
  [v len]
  (mult (normalize v) len))

(defn limit
  "Limit length of the vector by given value"
  {:metadoc/categories #{:dist}
   :metadoc/examples [(example-session "Usage"
                        (limit (vec3 1.0 1.0 1.0) 1.0)
                        (limit (vec3 1.0 1.0 1.0) 2.0))]} 
  [v ^double len]
  (if (> ^double (magsq v) (* len len))
    (set-mag v len)
    v))

(defn angle-between
  "Angle between two vectors

  See also [[relative-angle-between]]."
  {:metadoc/categories #{:geom}
   :metadoc/examples [(example-session "Usage"
                        (m/degrees (angle-between (vec3 1.0 0.0 0.0) (vec3 0.0 1.0 0.0)))
                        (m/degrees (angle-between (vec (repeatedly 50 rand)) (vec (repeatedly 50 rand)))))]} 
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
  {:metadoc/categories #{:geom}
   :metadoc/examples [(example-session "Usage"
                        (m/degrees (relative-angle-between (vec3 1.0 0.0 0.0) (vec3 0.0 1.0 0.0)))
                        (m/degrees (relative-angle-between (vec (repeatedly 50 rand)) (vec (repeatedly 50 rand)))))]} 
  ^double [v1 v2]
  (- ^double (heading v2) ^double (heading v1)))

(defn aligned?
  "Are vectors aligned (have the same direction)?"
  {:metadoc/categories #{:geom}
   :metadoc/examples [(example-session "Usage"
                        (aligned? (vec2 1.0 1.0) (vec2 2.0 2.000001))
                        (aligned? (vec2 1.0 1.0) (vec2 2.0 2.00001)))]} 
  [v1 v2]
  (< (angle-between v1 v2) TOLERANCE))

(defn faceforward
  "Flip normal `n` to match the same direction as `v`."
  {:metadoc/categories #{:geom}
   :metadoc/examples [(example-session "Usage"
                        (faceforward (vec2 1.0 1.0) (vec2 -1.0 -3.0))
                        (faceforward (vec2 1.0 1.0) (vec2 1.0 0.0)))]} 
  [n v]
  (if (neg? ^double (dot n v)) 
    (sub n)
    n))

(defn generate-vec2
  "Generate Vec2 with fn(s)"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage" (generate-vec2 (constantly 2)) (generate-vec2 rand (constantly 1)))]} 
  ([f1 f2]
   (Vec2. (f1) (f2)))
  ([f]
   (Vec2. (f) (f))))

(defn generate-vec3
  "Generate Vec3 with fn(s)"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage" (generate-vec3 rand) (generate-vec3 rand (constantly 1) (constantly 2)))]} 
  ([f1 f2 f3]
   (Vec3. (f1) (f2) (f3)))
  ([f]
   (Vec3. (f) (f) (f))))

(defn generate-vec4
  "Generate Vec4 with fn(s)"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example-session "Usage" (generate-vec4 rand) (generate-vec4 rand rand (constantly 1) (constantly 2)))]} 
  ([f1 f2 f3 f4]
   (Vec4. (f1) (f2) (f3) (f4)))
  ([f]
   (Vec4. (f) (f) (f) (f))))

(defn array->vec2
  "Doubles array to Vec2"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example "Usage" (array->vec2 (double-array [11 22 33 44 55])))]} 
  [^doubles a]
  (Vec2. (aget a 0) (aget a 1)))

(defn array->vec3
  "Doubles array to Vec3"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example "Usage" (array->vec3 (double-array [11 22 33 44 55])))]} 
  [^doubles a]
  (Vec3. (aget a 0) (aget a 1) (aget a 2)))

(defn array->vec4
  "Doubles array to Vec4"
  {:metadoc/categories #{:gen}
   :metadoc/examples [(example "Usage" (array->vec4 (double-array [11 22 33 44 55])))]} 
  [^doubles a]
  (Vec4. (aget a 0) (aget a 1) (aget a 2) (aget a 3)))

;;

(defmethod print-method Vec2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec2 [v w] (print-method v w))

(defmethod print-method Vec3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec3 [v w] (print-method v w))

(defmethod print-method Vec4 [^Vec4 v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Vec4 [v w] (print-method v w))

;;

(add-examples mag (example "Length of the vector" {:test-value (m/sqrt 2.0)} (mag (vec2 1 1))))
(add-examples magsq (example "Length of the vector, squared" {:test-value 10.0} (magsq [1 2 1 2])))
(add-examples abs (example "Usage" (abs (array-vec [-1 2 -2 1]))))
(add-examples add (example "Usage" (add (vec2 1 2) (vec2 -1 -2))))
(add-examples applyf (example "Usage" (applyf [5 3] (fn [v] (m/sin v)))))
(add-examples div (example "Usage" (div [5 4 3 5] 4.0)))
(add-examples dot (example "Usage" (dot (vec4 1 0 0 0) (vec4 -1 0 -1 0))))
(add-examples econstrain (example "Usage" (econstrain (vec3 2 0 -2) -1 1)))
(add-examples emn (example "Usage" (emn [-1 2] [1 -2])))
(add-examples emx (example "Usage" (emx [-1 2] [1 -2])))
(add-examples emult (example "Usage" (emult (vec3 1 2 3) (vec3 9 9 9))))
(add-examples is-near-zero? (example-session "Usage" (is-near-zero? [0 0.0000001]) (is-near-zero? [0 0.000001])))
(add-examples is-zero? (example-session "Usage" (is-zero? [0 0.0000001]) (is-zero? [0 0.0])))
(add-examples mn (example "Usage" (mn (vec4 -1 -2 3 4))))
(add-examples mx (example "Usage" (mx (vec4 -1 -2 3 4))))
(add-examples permute (example "Usage" (permute (vec4 1 2 3 4) (vec4 0 3 2 1))))
(add-examples reciprocal (example "Usage" (reciprocal (vec3 1 2 5))))
(add-examples sub (example "Usage" (sub (vec2 1 2) (vec2 -1 -2))))
(add-examples sum (example "Usage" (sum (vec (range 1000)))))

(add-examples transform
  (example-session "Usage"
    (transform (vec2 1 1) (vec2 -1 -1) (vec2 1.0 0.0) (vec2 0.0 1.0))
    (transform (vec3 1 1 1) (vec3 -1 -1 0) (vec3 1.0 0.0 1.0) (vec3 0.0 1.0 0.0) (vec3 0.0 1.0 1.0))))

(add-examples rotate
  (example-session "Usage"
    (rotate (vec2 1 2) (m/radians 90))
    (rotate (vec3 1 2 3) (m/radians 90) 0 0)))

(add-examples perpendicular
  (example-session "Usage"
    (perpendicular (vec2 -4 0))
    (perpendicular (vec3 1.0 0.0 0.0) (vec3 0.0 -1.0 0.0))))

(add-examples maxdim
  (example-session "Usage"
    (let [v (vec (repeatedly 100 (fn [] (- (int (rand-int 200)) 100))))
          mdim (maxdim v)]
      [mdim (v mdim)])
    (maxdim (vec3 1 2 3))))

(add-examples mindim
  (example-session "Usage"
    (let [v (vec (repeatedly 100 (fn [] (- (int (rand-int 200)) 100))))
          mdim (mindim v)]
      [mdim (v mdim)])
    (mindim (vec3 1 2 3))))

(add-examples heading
  (example-session "Usage"
    (m/degrees (heading (vec2 1.0 1.0)))
    (m/degrees (heading (vec3 1.0 0.0 1.0)))
    (m/degrees (heading (vec4 1.0 -1.0 1.0 -1.0)))
    (heading [1 2 3 4 5 6 7 8 9])
    (heading (array-vec [1 2 3 4 5 6 7 8 9]))))

(add-examples from-polar
  (example-session "Usage"
    (from-polar (vec2 1.0 (m/radians 90)))
    (from-polar (vec3 1.0 (m/radians 90) (m/radians 45)))))

(add-examples from-polar
  (example-session "Usage"
    (to-polar (vec2 1.0 1.0))
    (to-polar (vec3 1.0 0.0 1.0))))

(add-examples interpolate
  (example-session "Usage"
    (interpolate (vec2 -1 -1) (vec2 1 1) 0.25)
    (interpolate (vec2 -1 -1) (vec2 1 1) 0.25 m/smooth-interpolation)))

(add-examples einterpolate
  (example-session "Usage"
    (einterpolate (vec2 -1 -1) (vec2 1 1) (vec2 0.25 0.75))
    (einterpolate (vec2 -1 -1) (vec2 1 1) (vec2 0.25 0.75) m/smooth-interpolation)))

(add-examples axis-rotate
  (example-session "Usage"
    (axis-rotate (vec3 1.0 0.0 0.0) m/HALF_PI (vec3 0.0 1.0 0.0))
    (axis-rotate (vec3 1.0 0.0 0.0) m/HALF_PI (vec3 0.0 1.0 0.0) (vec3 1.0 1.0 1.0))))

(add-examples base-from
  (example-session "Usage"
    (base-from (vec3 0.1 1.0 -1.0))
    (base-from (vec2 1.0 0.0))))

(add-examples cross
  (example-session "Usage"
    (cross (vec3 1 2 3)
           (vec3 4 3 2))
    (cross (vec2 1 2)
           (vec2 -1 2))))

(add-examples to-vec
  (example-session "Check types"
    (type (to-vec [1 2 3]))
    (type (to-vec (vec2 1 2)))
    (type (to-vec (vec3 1 2 3)))
    (type (to-vec (vec4 1 2 3 4)))
    (type (to-vec (array-vec 1)))))
