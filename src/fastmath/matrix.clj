(ns fastmath.matrix
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.protocols.matrix :as prot])
  (:import [clojure.lang Counted IFn IPersistentCollection Seqable Reversible ILookup]
           [org.apache.commons.math3.linear Array2DRowRealMatrix ArrayRealVector]
           [fastmath.vector Vec2 Vec3 Vec4]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators #{'abs})

(defn- mat-throw-ioobe
  [dim id]
  (throw (IndexOutOfBoundsException. (str "Index " id " out of bounds for mat" dim))))

(def ^:private mat2x2-throw-ioobe (partial mat-throw-ioobe "2x2"))
(def ^:private mat3x3-throw-ioobe (partial mat-throw-ioobe "3x3"))
(def ^:private mat4x4-throw-ioobe (partial mat-throw-ioobe "4x4"))

(defn- gen-sym [i j] (symbol (str "a" i j)))

(defmacro ^:private gen-condition
  [size v else]
  (let [r (range size)]
    `(condp = ~v ~@(mapcat identity (for [x r y r
                                          :let [s (gen-sym x y)]]
                                      [`[~x ~y] `~s])) ~else)))

(defmacro ^:private gen-mulm
  [size clss t1? t2?]
  (let [r (range size)]
    `(new ~clss ~@(for [m1 r m2 r]
                    `(+ ~@(for [i r]
                            `(* ~(if t1? (gen-sym i m1) (gen-sym m1 i))
                                (. ~'m2 ~(if t2? (gen-sym m2 i) (gen-sym i m2))))))))))

(defmacro ^:private gen-hashc
  [^long size]
  (let [r (range size)]
    `(mix-collection-hash (unchecked-int (-> ~@(for [x r y r]
                                                 `(v/dhash-code ~(gen-sym x y))))) ~(* size size))))

(defmacro ^:private gen-det2
  [a b c d]
  `(- (* ~a ~d) (* ~b ~c)))

(defmacro ^:private gen-det3
  [a b c d e f g h i]
  `(+ (* ~a (gen-det2 ~e ~f ~h ~i))
      (* (- ~b) (gen-det2 ~d ~f ~g ~i))
      (* ~c (gen-det2 ~d ~e ~g ~h))))

(defn- in-range?
  [x y ^long mx]
  (and (number? x) (pos? ^long x) (< ^long x mx)
       (number? y) (pos? ^long y) (< ^long y mx)))

(deftype Mat2x2 [^double a00 ^double a01
                 ^double a10 ^double a11]
  Object
  (toString [_] (str "#mat2x2 [[" a00 ", " a01 "] [" a10 ", " a11 "]]"))
  (equals [_ m]
    (and (instance? Mat2x2 m)
         (let [^Mat2x2 m m]
           (and (== a00 (.a00 m)) (== a01 (.a01 m))
                (== a10 (.a10 m)) (== a11 (.a11 m))))))
  (hashCode [_] (gen-hashc 2))
  clojure.lang.IHashEq 
  (hasheq [_] (gen-hashc 2))
  Seqable
  (seq [_] (list a00 a01 a10 a11))
  Reversible
  (rseq [_] (list a11 a10 a01 a00))
  Counted
  (count [_] 4)
  ILookup
  (valAt [_ [x y :as v]] (when (in-range? x y 2)
                           (gen-condition 2 v nil)))
  (valAt [_ [x y :as v] not-found] (if (in-range? x y 2)
                                     (gen-condition 2 v not-found)
                                     not-found))
  IFn
  (invoke [_ id]
    (gen-condition 2 id (mat2x2-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition 2 [x y] (mat2x2-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (fmap [_ f] (Mat2x2. (f a00) (f a01) (f a10) (f a11)))
  (cols [_] [(Vec2. a00 a10)
             (Vec2. a01 a11)])
  (rows [_] [(Vec2. a00 a01)
             (Vec2. a10 a11)])
  (to-array2d [_] (m/seq->double-double-array [[a00 a01] [a10 a11]]))
  (to-acm-mat [_] (Array2DRowRealMatrix. (m/seq->double-double-array [[a00 a01] [a10 a11]])))
  (nrow [_] 2)
  (ncol [_] 2)
  (column [_ id]
    (case (unchecked-int id)
      0 (Vec2. a00 a10)
      1 (Vec2. a01 a11)
      (mat2x2-throw-ioobe id)))
  (row [_ id]
    (case (unchecked-int id)
      0 (Vec2. a00 a01)
      1 (Vec2. a10 a11)
      (mat2x2-throw-ioobe id)))
  (symmetric? [_] (== a01 a10))
  (transpose [_] (Mat2x2. a00 a10 a01 a11))
  (inverse [_] (let [d (gen-det2 a00 a01 a10 a11)]
                 (when-not (zero? d)
                   (let [d (/ d)]
                     (Mat2x2. (* d a11) (* d (- a01)) (* d (- a10)) (* d a00))))))
  (diag [_] (Vec2. a00 a11))
  (trace [_] (+ a00 a11))
  (det [_] (gen-det2 a00 a01 a10 a11))
  (add [_ m] (let [^Mat2x2 m m]
               (Mat2x2. (+ a00 (.a00 m))
                        (+ a01 (.a01 m))
                        (+ a10 (.a10 m))
                        (+ a11 (.a11 m)))))
  (adds [_ v] (let [v (double v)]
                (Mat2x2. (+ v a00) (+ v a01) (+ v a10) (+ v a11))))
  (sub [_] (Mat2x2. (- a00) (- a01) (- a10) (- a11)))
  (sub [_ m] (let [^Mat2x2 m m]
               (Mat2x2. (- a00 (.a00 m))
                        (- a01 (.a01 m))
                        (- a10 (.a10 m))
                        (- a11 (.a11 m)))))
  (emul [_ m] (let [^Mat2x2 m m]
                (Mat2x2. (* a00 (.a00 m))
                         (* a01 (.a01 m))
                         (* a10 (.a10 m))
                         (* a11 (.a11 m)))))
  (mulm [_ m2] (let [^Mat2x2 m2 m2] (gen-mulm 2 Mat2x2 false false)))
  (mulm [_ t1? m2 t2?] (let [^Mat2x2 m2 m2]
                         (cond
                           (not (or t1? t2?)) (gen-mulm 2 Mat2x2 false false)
                           (and t1? (not t2?)) (gen-mulm 2 Mat2x2 true false)
                           (and (not t1?) t2?) (gen-mulm 2 Mat2x2 false true)
                           :else (gen-mulm 2 Mat2x2 true true))))
  (muls [_ s] (let [s (double s)]
                (Mat2x2. (* s a00) (* s a01) (* s a10) (* s a11))))
  (mulv [_ v] (let [^Vec2 v v]
                (Vec2. (+ (* a00 (.x v)) (* a01 (.y v)))
                       (+ (* a10 (.x v)) (* a11 (.y v))))))
  (vtmul [_ v] (let [^Vec2 v v]
                 (Vec2. (+ (* a00 (.x v)) (* a10 (.y v)))
                        (+ (* a01 (.x v)) (* a11 (.y v)))))))


(deftype Mat3x3 [^double a00 ^double a01 ^double a02
                 ^double a10 ^double a11 ^double a12
                 ^double a20 ^double a21 ^double a22]
  Object
  (toString [_] (str "#mat2x2 [[" a00 ", " a01 ", " a02  "] [" a10 ", " a11 ", " a12 "] [" a20 ", " a21 ", " a22 "]]"))
  (equals [_ m]
    (and (instance? Mat3x3 m)
         (let [^Mat3x3 m m]
           (and (== a00 (.a00 m)) (== a01 (.a01 m)) (== a02 (.a02 m))
                (== a10 (.a10 m)) (== a11 (.a11 m)) (== a12 (.a12 m))
                (== a20 (.a20 m)) (== a21 (.a21 m)) (== a22 (.a22 m))))))
  (hashCode [_] (gen-hashc 3))
  clojure.lang.IHashEq 
  (hasheq [_] (gen-hashc 3))
  Seqable
  (seq [_] (list a00 a01 a02 a10 a11 a12 a20 a21 a22))
  Reversible
  (rseq [_] (list a22 a21 a20 a12 a11 a10 a02 a01 a00))
  Counted
  (count [_] 9)
  ILookup
  (valAt [_ [x y :as v]] (when (in-range? x y 3)
                           (gen-condition 3 v nil)))
  (valAt [_ [x y :as v] not-found] (if (in-range? x y 3)
                                     (gen-condition 3 v not-found)
                                     not-found))
  IFn
  (invoke [_ id]
    (gen-condition 3 id (mat3x3-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition 3 [x y] (mat3x3-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (fmap [_ f] (Mat3x3. (f a00) (f a01) (f a02) (f a10) (f a11) (f a12) (f a20) (f a21) (f a22)))
  (cols [_] [(Vec3. a00 a10 a20)
             (Vec3. a01 a11 a21)
             (Vec3. a02 a12 a22)])
  (rows [_] [(Vec3. a00 a01 a02)
             (Vec3. a10 a11 a12)
             (Vec3. a20 a21 a22)])
  (to-array2d [_] (m/seq->double-double-array [[a00 a01 a02] [a10 a11 a12] [a20 a21 a22]]))
  (to-acm-mat [_] (Array2DRowRealMatrix.
                   (m/seq->double-double-array [[a00 a01 a02] [a10 a11 a12] [a20 a21 a22]])))
  (nrow [_] 3)
  (ncol [_] 3)
  (column [_ id]
    (case (unchecked-int id)
      0 (Vec3. a00 a10 a20)
      1 (Vec3. a01 a11 a21)
      2 (Vec3. a02 a12 a22)
      (mat3x3-throw-ioobe id)))
  (row [_ id]
    (case (unchecked-int id)
      0 (Vec3. a00 a01 a02)
      1 (Vec3. a10 a11 a12)
      2 (Vec3. a20 a21 a22)
      (mat3x3-throw-ioobe id)))
  (symmetric? [_] (and (== a01 a10) (== a12 a21) (== a20 a02)))
  (transpose [_] (Mat3x3. a00 a10 a20 a01 a11 a21 a02 a12 a22))
  (inverse [_] (let [d (gen-det3 a00 a01 a02 a10 a11 a12 a20 a21 a22)]
                 (when-not (zero? d)
                   (let [adj (Mat3x3. (gen-det2 a11 a12 a21 a22)
                                      (gen-det2 a02 a01 a22 a21)
                                      (gen-det2 a01 a02 a11 a12)
                                      (gen-det2 a12 a10 a22 a20)
                                      (gen-det2 a00 a02 a20 a22)
                                      (gen-det2 a02 a00 a12 a10)
                                      (gen-det2 a10 a11 a20 a21)
                                      (gen-det2 a01 a00 a21 a20)
                                      (gen-det2 a00 a01 a10 a11))]
                     (prot/muls adj (/ d))))))
  (diag [_] (Vec3. a00 a11 a22))
  (trace [_] (+ a00 a11 a22))
  (det [_] (gen-det3 a00 a01 a02 a10 a11 a12 a20 a21 a22))
  
  (add [_ m] (let [^Mat3x3 m m]
               (Mat3x3. (+ a00 (.a00 m)) (+ a01 (.a01 m)) (+ a02 (.a02 m))
                        (+ a10 (.a10 m)) (+ a11 (.a11 m)) (+ a12 (.a12 m))
                        (+ a20 (.a20 m)) (+ a21 (.a21 m)) (+ a22 (.a22 m)))))
  (adds [_ v] (let [v (double v)]
                (Mat3x3. (+ v a00) (+ v a01) (+ v a02)
                         (+ v a10) (+ v a11) (+ v a12)
                         (+ v a20) (+ v a21) (+ v a22))))
  (sub [_] (Mat3x3. (- a00) (- a01) (- a02)
                    (- a10) (- a11) (- a12)
                    (- a20) (- a21) (- a22)))
  (sub [_ m] (let [^Mat3x3 m m]
               (Mat3x3. (- a00 (.a00 m)) (- a01 (.a01 m)) (- a02 (.a02 m))
                        (- a10 (.a10 m)) (- a11 (.a11 m)) (- a12 (.a12 m))
                        (- a20 (.a20 m)) (- a21 (.a21 m)) (- a22 (.a22 m)))))
  (emul [_ m] (let [^Mat3x3 m m]
                (Mat3x3. (* a00 (.a00 m)) (* a01 (.a01 m)) (* a02 (.a02 m))
                         (* a10 (.a10 m)) (* a11 (.a11 m)) (* a12 (.a12 m))
                         (* a20 (.a20 m)) (* a21 (.a21 m)) (* a22 (.a22 m)))))
  (mulm [_ m2] (let [^Mat3x3 m2 m2] (gen-mulm 3 Mat3x3 false false)))
  (mulm [_ t1? m2 t2?] (let [^Mat3x3 m2 m2]
                         (cond
                           (not (or t1? t2?)) (gen-mulm 3 Mat3x3 false false)
                           (and t1? (not t2?)) (gen-mulm 3 Mat3x3 true false)
                           (and (not t1?) t2?) (gen-mulm 3 Mat3x3 false true)
                           :else (gen-mulm 3 Mat3x3 true true))))
  (muls [_ s] (let [s (double s)]
                (Mat3x3. (* s a00) (* s a01) (* s a02)
                         (* s a10) (* s a11) (* s a12)
                         (* s a20) (* s a21) (* s a22))))
  (mulv [_ v] (let [^Vec3 v v]
                (Vec3. (+ (* a00 (.x v)) (* a01 (.y v)) (* a02 (.z v)))
                       (+ (* a10 (.x v)) (* a11 (.y v)) (* a12 (.z v)))
                       (+ (* a20 (.x v)) (* a21 (.y v)) (* a22 (.z v))))))
  (vtmul [_ v] (let [^Vec3 v v]
                 (Vec3. (+ (* a00 (.x v)) (* a10 (.y v)) (* a20 (.z v)))
                        (+ (* a01 (.x v)) (* a11 (.y v)) (* a21 (.z v)))
                        (+ (* a02 (.x v)) (* a12 (.y v)) (* a22 (.z v)))))))

(deftype Mat4x4 [^double a00 ^double a01 ^double a02 ^double a03
                 ^double a10 ^double a11 ^double a12 ^double a13
                 ^double a20 ^double a21 ^double a22 ^double a23
                 ^double a30 ^double a31 ^double a32 ^double a33]
  Object
  (toString [_] (str "#mat2x2 [[" a00 ", " a01 ", " a02 ", " a03 "] [" a10 ", " a11 ", " a12 ", " a13 "] [" a20 ", " a21 ", " a22 ", " a23 "] [" a30 ", " a31 ", " a32 ", " a33 "]]"))
  (equals [_ m]
    (and (instance? Mat4x4 m)
         (let [^Mat4x4 m m]
           (and (== a00 (.a00 m)) (== a01 (.a01 m)) (== a02 (.a02 m)) (== a03 (.a03 m))
                (== a10 (.a10 m)) (== a11 (.a11 m)) (== a12 (.a12 m)) (== a13 (.a13 m))
                (== a20 (.a20 m)) (== a21 (.a21 m)) (== a22 (.a22 m)) (== a23 (.a23 m))
                (== a30 (.a30 m)) (== a31 (.a31 m)) (== a32 (.a32 m)) (== a33 (.a33 m))))))
  (hashCode [_] (gen-hashc 4))
  clojure.lang.IHashEq 
  (hasheq [_] (gen-hashc 4))
  Seqable
  (seq [_] (list a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))
  Reversible
  (rseq [_] (list a33 a32 a31 a30 a23 a22 a21 a20 a13 a12 a11 a10 a03 a02 a01 a00))
  Counted
  (count [_] 16)
  ILookup
  (valAt [_ [x y :as v]] (when (in-range? x y 4)
                           (gen-condition 4 v nil)))
  (valAt [_ [x y :as v] not-found] (if (in-range? x y 4)
                                     (gen-condition 4 v not-found)
                                     not-found))
  IFn
  (invoke [_ id]
    (gen-condition 4 id (mat4x4-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition 4 [x y] (mat4x4-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (fmap [_ f] (Mat4x4. (f a00) (f a01) (f a02) (f a03)
                       (f a10) (f a11) (f a12) (f a13)
                       (f a20) (f a21) (f a22) (f a23)
                       (f a30) (f a31) (f a32) (f a33)))
  (cols [_] [(Vec4. a00 a10 a20 a30)
             (Vec4. a01 a11 a21 a31)
             (Vec4. a02 a12 a22 a32)
             (Vec4. a03 a13 a23 a33)])
  (rows [_] [(Vec4. a00 a01 a02 a03)
             (Vec4. a10 a11 a12 a13)
             (Vec4. a20 a21 a22 a23)
             (Vec4. a30 a31 a32 a33)])
  (to-array2d [_] (m/seq->double-double-array [[a00 a01 a02 a03]
                                               [a10 a11 a12 a13]
                                               [a20 a21 a22 a23]
                                               [a30 a31 a32 a33]]))
  (to-acm-mat [_] (Array2DRowRealMatrix. (m/seq->double-double-array [[a00 a01 a02 a03]
                                                                      [a10 a11 a12 a13]
                                                                      [a20 a21 a22 a23]
                                                                      [a30 a31 a32 a33]])))
  (nrow [_] 4)
  (ncol [_] 4)
  (column [_ id]
    (case (unchecked-int id)
      0 (Vec4. a00 a10 a20 a30)
      1 (Vec4. a01 a11 a21 a31)
      2 (Vec4. a02 a12 a22 a32)
      3 (Vec4. a03 a13 a23 a33)
      (mat4x4-throw-ioobe id)))
  (row [_ id]
    (case (unchecked-int id)
      0 (Vec4. a00 a01 a02 a03)
      1 (Vec4. a10 a11 a12 a13)
      2 (Vec4. a20 a21 a22 a23)
      3 (Vec4. a30 a31 a32 a33)
      (mat4x4-throw-ioobe id)))
  (symmetric? [_] (and (== a01 a10) (== a12 a21) (== a20 a02) (== a30 a03) (== a31 a13) (== a32 a23)))
  (transpose [_] (Mat4x4. a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33))
  (inverse [_] (let [d (+ (* a00 (gen-det3 a11 a12 a13 a21 a22 a23 a31 a32 a33))
                          (* (- a10) (gen-det3 a01 a02 a03 a21 a22 a23 a31 a32 a33 ))
                          (* a20 (gen-det3 a01 a02 a03 a11 a12 a13 a31 a32 a33))
                          (* (- a30) (gen-det3 a01 a02 a03 a11 a12 a13 a21 a22 a23)))]
                 (when-not (zero? d)
                   (let [adj (Mat4x4. (gen-det3 a11 a12 a13 a21 a22 a23 a31 a32 a33)
                                      (- (gen-det3 a01 a02 a03 a21 a22 a23 a31 a32 a33))
                                      (gen-det3 a01 a02 a03 a11 a12 a13 a31 a32 a33)
                                      (- (gen-det3 a01 a02 a03 a11 a12 a13 a21 a22 a23))

                                      (- (gen-det3 a10 a12 a13 a20 a22 a23 a30 a32 a33))
                                      (gen-det3 a00 a02 a03 a20 a22 a23 a30 a32 a33)
                                      (- (gen-det3 a00 a02 a03 a10 a12 a13 a30 a32 a33))
                                      (gen-det3 a00 a02 a03 a10 a12 a13 a20 a22 a23)

                                      (gen-det3 a10 a11 a13 a20 a21 a23 a30 a31 a33)
                                      (- (gen-det3 a00 a01 a03 a20 a21 a23 a30 a31 a33))
                                      (gen-det3 a00 a01 a03 a10 a11 a13 a30 a31 a33)
                                      (- (gen-det3 a00 a01 a03 a10 a11 a13 a20 a21 a23))

                                      (- (gen-det3 a10 a11 a12 a20 a21 a22 a30 a31 a32))
                                      (gen-det3 a00 a01 a02 a20 a21 a22 a30 a31 a32)
                                      (- (gen-det3 a00 a01 a02 a10 a11 a12 a30 a31 a32))
                                      (gen-det3 a00 a01 a02 a10 a11 a12 a20 a21 a22))]
                     (prot/muls adj (/ d))))))
  (diag [_] (Vec4. a00 a11 a22 a33))
  (trace [_] (+ a00 a11 a22 a33))
  (det [_] (+ (* a00 (gen-det3 a11 a12 a13 a21 a22 a23 a31 a32 a33))
              (* (- a10) (gen-det3 a01 a02 a03 a21 a22 a23 a31 a32 a33 ))
              (* a20 (gen-det3 a01 a02 a03 a11 a12 a13 a31 a32 a33))
              (* (- a30) (gen-det3 a01 a02 a03 a11 a12 a13 a21 a22 a23))))
  
  (add [_ m] (let [^Mat4x4 m m]
               (Mat4x4. (+ a00 (.a00 m)) (+ a01 (.a01 m)) (+ a02 (.a02 m)) (+ a03 (.a03 m))
                        (+ a10 (.a10 m)) (+ a11 (.a11 m)) (+ a12 (.a12 m)) (+ a13 (.a13 m))
                        (+ a20 (.a20 m)) (+ a21 (.a21 m)) (+ a22 (.a22 m)) (+ a23 (.a23 m))
                        (+ a30 (.a30 m)) (+ a31 (.a31 m)) (+ a32 (.a32 m)) (+ a33 (.a33 m)))))
  (adds [_ v] (let [v (double v)]
                (Mat4x4. (+ v a00) (+ v a01) (+ v a02) (+ v a03)
                         (+ v a10) (+ v a11) (+ v a12) (+ v a13)
                         (+ v a20) (+ v a21) (+ v a22) (+ v a23)
                         (+ v a30) (+ v a31) (+ v a32) (+ v a33))))
  (sub [_] (Mat4x4. (- a00) (- a01) (- a02) (- a03)
                    (- a10) (- a11) (- a12) (- a13)
                    (- a20) (- a21) (- a22) (- a23)
                    (- a30) (- a31) (- a32) (- a33)))
  (sub [_ m] (let [^Mat4x4 m m]
               (Mat4x4. (- a00 (.a00 m)) (- a01 (.a01 m)) (- a02 (.a02 m)) (- a03 (.a03 m))
                        (- a10 (.a10 m)) (- a11 (.a11 m)) (- a12 (.a12 m)) (- a13 (.a13 m))
                        (- a20 (.a20 m)) (- a21 (.a21 m)) (- a22 (.a22 m)) (- a23 (.a23 m))
                        (- a30 (.a30 m)) (- a31 (.a31 m)) (- a32 (.a32 m)) (- a33 (.a33 m)))))
  (emul [_ m] (let [^Mat4x4 m m]
                (Mat4x4. (* a00 (.a00 m)) (* a01 (.a01 m)) (* a02 (.a02 m)) (* a03 (.a03 m))
                         (* a10 (.a10 m)) (* a11 (.a11 m)) (* a12 (.a12 m)) (* a13 (.a13 m))
                         (* a20 (.a20 m)) (* a21 (.a21 m)) (* a22 (.a22 m)) (* a23 (.a23 m))
                         (* a30 (.a30 m)) (* a31 (.a31 m)) (* a32 (.a32 m)) (* a33 (.a33 m)))))
  (mulm [_ m2] (let [^Mat4x4 m2 m2] (gen-mulm 4 Mat4x4 false false)))
  (mulm [_ t1? m2 t2?] (let [^Mat4x4 m2 m2]
                         (cond
                           (not (or t1? t2?)) (gen-mulm 4 Mat4x4 false false)
                           (and t1? (not t2?)) (gen-mulm 4 Mat4x4 true false)
                           (and (not t1?) t2?) (gen-mulm 4 Mat4x4 false true)
                           :else (gen-mulm 4 Mat4x4 true true))))
  (muls [_ s] (let [s (double s)]
                (Mat4x4. (* s a00) (* s a01) (* s a02) (* s a03)
                         (* s a10) (* s a11) (* s a12) (* s a13)
                         (* s a20) (* s a21) (* s a22) (* s a23)
                         (* s a30) (* s a31) (* s a32) (* s a33))))
  (mulv [_ v] (let [^Vec4 v v]
                (Vec4. (+ (* a00 (.x v)) (* a01 (.y v)) (* a02 (.z v)) (* a03 (.w v)))
                       (+ (* a10 (.x v)) (* a11 (.y v)) (* a12 (.z v)) (* a13 (.w v)))
                       (+ (* a20 (.x v)) (* a21 (.y v)) (* a22 (.z v)) (* a23 (.w v)))
                       (+ (* a30 (.x v)) (* a31 (.y v)) (* a32 (.z v)) (* a33 (.w v))))))
  (vtmul [_ v] (let [^Vec4 v v]
                 (Vec4. (+ (* a00 (.x v)) (* a10 (.y v)) (* a20 (.z v)) (* a30 (.w v)))
                        (+ (* a01 (.x v)) (* a11 (.y v)) (* a21 (.z v)) (* a31 (.w v)))
                        (+ (* a02 (.x v)) (* a12 (.y v)) (* a22 (.z v)) (* a32 (.w v)))
                        (+ (* a03 (.x v)) (* a13 (.y v)) (* a23 (.z v)) (* a33 (.w v)))))))

(defn mat2x2
  "Create 2x2 matrix.

  Arity:

  * 1 - fills matrix with given value
  * 2 - creates diagonal matrix
  * 4 - creates row ordered matrix"
  (^Mat2x2 [^double v] (Mat2x2. v v v v))
  (^Mat2x2 [^double d1 ^double d2] (Mat2x2. d1 0.0 0.0 d2))
  (^Mat2x2 [^double a00 ^double a01 ^double a10 ^double a11] (Mat2x2. a00 a01 a10 a11)))

(defn mat3x3
  "Create 3x3 matrix.

  Arity:

  * 1 - fills matrix with given value
  * 3 - creates diagonal matrix
  * 9 - creates row ordered matrix"
  (^Mat3x3 [^double v] (Mat3x3. v v v v v v v v v))
  (^Mat3x3 [^double d1 ^double d2 ^double d3] (Mat3x3. d1 0.0 0.0
                                                       0.0 d2 0.0
                                                       0.0 0.0 d3))
  (^Mat3x3 [a00 a01 a02 a10 a11 a12 a20 a21 a22] (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22)))

(defn mat4x4
  "Create 4x4 matrix.

  Arity:

  * 1 - fills matrix with given value
  * 4 - creates diagonal matrix
  * 16 - creates row ordered matrix"
  (^Mat4x4 [^double v] (Mat4x4. v v v v v v v v v v v v v v v v ))
  (^Mat4x4 [^double d1 ^double d2 ^double d3 ^double d4] (Mat4x4. d1 0.0 0.0 0.0
                                                                  0.0 d2 0.0 0.0
                                                                  0.0 0.0 d3 0.0
                                                                  0.0 0.0 0.0 d4))
  (^Mat4x4 [a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33]
   (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33)))

(def eye
  [nil 1.0
   (Mat2x2. 1.0 0.0 0.0 1.0)
   (Mat3x3. 1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0)
   (Mat4x4. 1.0 0.0 0.0 0.0
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0)])

(defn diagonal
  "Create diagonal matrix"
  [diag]
  (case (count diag)
    1 (first diag)
    2 (apply mat2x2 diag)
    3 (apply mat3x3 diag)
    4 (apply mat4x4 diag)))

(defn solve
  "Solve linear equation Ax=b"
  [A b]
  (when-let [A-1 (prot/inverse A)]
    (prot/mulv A-1 b)))

(defn outer
  "Outer project for two vectors."
  [v1 v2]
  (condp instance? v1
    Vec2 (let [^Vec2 v1 v1
               ^Vec2 v2 v2]
           (Mat2x2. (* (.x v1) (.x v2)) (* (.x v1) (.y v2))
                    (* (.y v1) (.x v2)) (* (.y v1) (.y v2))))
    Vec3 (let [^Vec3 v1 v1
               ^Vec3 v2 v2]
           (Mat3x3. (* (.x v1) (.x v2)) (* (.x v1) (.y v2)) (* (.x v1) (.z v2))
                    (* (.y v1) (.x v2)) (* (.y v1) (.y v2)) (* (.y v1) (.z v2))
                    (* (.z v1) (.x v2)) (* (.z v1) (.y v2)) (* (.z v1) (.z v2))))
    Vec4 (let [^Vec4 v1 v1
               ^Vec4 v2 v2]
           (Mat4x4. (* (.x v1) (.x v2)) (* (.x v1) (.y v2)) (* (.x v1) (.z v2)) (* (.x v1) (.w v2))
                    (* (.y v1) (.x v2)) (* (.y v1) (.y v2)) (* (.y v1) (.z v2)) (* (.y v1) (.w v2))
                    (* (.z v1) (.x v2)) (* (.z v1) (.y v2)) (* (.z v1) (.z v2)) (* (.z v1) (.w v2))
                    (* (.w v1) (.x v2)) (* (.w v1) (.y v2)) (* (.w v1) (.z v2)) (* (.w v1) (.w v2))))
    ArrayRealVector (let [^ArrayRealVector v1 v1
                          ^ArrayRealVector v2 v2]
                      (.outerProduct v1 v2))))

;;

(defn fmap
  "Apply a function `f` to each matrix element."
  [m f] (prot/fmap m f))

(defn cols
  "Return matrix columns"
  [m] (prot/cols m))

(defn rows
  "Return matrix rows"
  [m] (prot/rows m))

(defn mat->array2d
  "Return doubles of doubles"
  [m] (prot/to-array2d m))

(defn mat->RealMatrix
  "Return  Apache Commons Math Array2DRowMatrix"
  [m] (prot/to-acm-mat m))

(defn nrow
  "Return number of rows"
  [m] (prot/nrow m))

(defn ncol
  "Return number of rows"
  [m] (prot/ncol m))

(defn row
  "Return row as a vector"
  [m ^long r] (prot/row m r))

(defn col
  "Return column as a vector"
  [m ^long c] (prot/column m c))

(defn symmetric?
  "Check if matrix is symmetric"
  [m] (prot/symmetric? m))

(defn transpose
  "Transpose matrix"
  [m] (prot/transpose m))

(defn inverse
  "Matrix inversion.

  Returns `nil` if inversion doesn't exist."
  [m]  (prot/inverse m))

(defn add
  "Add matrices, C=A+B."
  ([m] m)
  ([m1 m2] (prot/add m1 m2)))

(defn mulm
  "Multiply two matrices, C=AxB."
  [m1 m2] (prot/mulm m1 m2))

(defmacro ^:private primitive-ops
  "Generate primitive functions operating on vectors"
  [fns]
  (let [v (symbol "vector")]
    `(do ~@(for [f fns
                 :let [nm (with-meta (symbol (name f)) {:metadoc/categories #{:mop}})
                       doc (str "Apply " nm " to matrix elements.")
                       wfn (if (:macro (meta (resolve f))) `(fn [v#] (~f v#)) f)]] ;; wrap macro into function
             `(defn ~nm ~doc
                [~v]
                (prot/fmap ~v ~wfn))))))

(primitive-ops [m/sin m/cos m/tan m/asin m/acos m/atan m/sinh m/cosh m/tanh m/asinh m/acosh m/atanh
                m/cot m/sec m/csc m/acot m/asec m/acsc m/coth m/sech m/csch m/acoth m/asech m/csch
                m/sq m/safe-sqrt m/sqrt m/cbrt m/exp m/log m/log10 m/log2 m/ln m/log1p m/expm1
                m/log1pexp m/log1mexp m/log1psq m/log1pmx m/logmxp1 m/logexpm1
                m/radians m/degrees m/sinc m/jinc m/sigmoid m/logit m/xlogx
                m/floor m/ceil m/round m/rint m/trunc m/frac m/sfrac m/signum m/sgn])

(defmethod print-method Mat2x2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat2x2 [v w] (print-method v w))

(defmethod print-method Mat3x3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat3x3 [v w] (print-method v w))

(defmethod print-method Mat4x4 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat4x4 [v w] (print-method v w))

