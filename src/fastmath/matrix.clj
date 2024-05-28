(ns fastmath.matrix
  "Fixed size (2x2, 3x3, 4x4) matrix types."
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.protocols.matrix :as prot])
  (:import [clojure.lang Counted IFn IPersistentCollection Seqable Reversible ILookup]
           [org.apache.commons.math3.linear Array2DRowRealMatrix RealVector ArrayRealVector EigenDecomposition
            RealMatrix MatrixUtils LUDecomposition CholeskyDecomposition]
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

(defmacro ^:private gen-condition-2
  [size x y else]
  (let [r (range size)]
    `(case (int ~x)
       ~@(mapcat identity (for [xx r]
                            [xx `(case (int ~y)
                                   ~@(mapcat identity (for [yy r :let [s (gen-sym xx yy)]]
                                                        [yy `~s]))
                                   ~else)]))
       ~else)))

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
  (toString [_] (str "#mat2x2 [[" a00 ", " a01 "]\n         [" a10 ", " a11 "]]"))
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
  (valAt [_ [^long x ^long y]] (when (in-range? x y 2)
                                 (gen-condition-2 2 x y nil)))
  (valAt [_ [^long x ^long y] not-found] (if (in-range? x y 2)
                                           (gen-condition-2 2 x y not-found)
                                           not-found))
  IFn
  (invoke [_ id]
    (gen-condition 2 id (mat2x2-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition-2 2 ^long x ^long y (mat2x2-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (->seq [_] (list a00 a01 a10 a11))
  (entry [_ x y] (gen-condition-2 2 ^long x ^long y (mat2x2-throw-ioobe [x y])))
  (fmap [_ f] (Mat2x2. (f a00) (f a01) (f a10) (f a11)))
  (cols [_] [(Vec2. a00 a10)
             (Vec2. a01 a11)])
  (rows [_] [(Vec2. a00 a01)
             (Vec2. a10 a11)])
  (to-double-array2d [_] (fastmath.java.Array/mat2array2d a00 a01
                                                          a10 a11))
  (to-float-array2d [_] (fastmath.java.Array/mat2array2d (float a00) (float a01)
                                                         (float a10) (float a11)))

  (to-double-array [_] (fastmath.java.Array/mat2array a00 a01
                                                      a10 a11))
  (to-float-array [_] (fastmath.java.Array/mat2array (float a00) (float a01)
                                                     (float a10) (float a11)))

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
  (symmetric? [_ tol] (m/delta-eq a01 a10 tol))
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
  (emulm [_ m] (let [^Mat2x2 m m]
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
                        (+ (* a01 (.x v)) (* a11 (.y v))))))
  (cholesky [m] (do
                  (assert (prot/symmetric? m) "Matrix is not symmetric.")
                  (let [a (m/sqrt a00)
                        b (/ a10 a)
                        c (m/sqrt (- a11 (* b b)))]
                    (Mat2x2. a 0.0 b c))))
  (norm [_ t] (if (sequential? t)
                (let [[^double p ^double q] t
                      qp (/ q p)]
                  (m/pow (+ (m/pow (+ (m/pow (m/abs a00) p) (m/pow (m/abs a10) p)) qp)
                            (m/pow (+ (m/pow (m/abs a01) p) (m/pow (m/abs a11) p)) qp)) (/ q)))
                (condp = t
                  :inf (m/max (+ (m/abs a00) (m/abs a01))
                              (+ (m/abs a10) (m/abs a11)))
                  :max (m/max (m/abs a00) (m/abs a01)
                              (m/abs a10) (m/abs a11))
                  (m/max (+ (m/abs a00) (m/abs a10))
                         (+ (m/abs a01) (m/abs a11)))))))

(deftype Mat3x3 [^double a00 ^double a01 ^double a02
                 ^double a10 ^double a11 ^double a12
                 ^double a20 ^double a21 ^double a22]
  Object
  (toString [_] (str "#mat3x3 [[" a00 ", " a01 ", " a02  "]\n         [" a10 ", " a11 ", " a12 "]\n         [" a20 ", " a21 ", " a22 "]]"))
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
  (valAt [_ [^long x ^long y]] (when (in-range? x y 3)
                                 (gen-condition-2 3 x y nil)))
  (valAt [_ [^long x ^long y] not-found] (if (in-range? x y 3)
                                           (gen-condition-2 3 x y not-found)
                                           not-found))
  IFn
  (invoke [_ id]
    (gen-condition 3 id (mat3x3-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition-2 3 ^long x ^long y (mat3x3-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (->seq [_] (list a00 a01 a02 a10 a11 a12 a20 a21 a22))
  (entry [_ x y] (gen-condition-2 3 ^long x ^long y (mat3x3-throw-ioobe [x y])))
  (fmap [_ f] (Mat3x3. (f a00) (f a01) (f a02) (f a10) (f a11) (f a12) (f a20) (f a21) (f a22)))
  (cols [_] [(Vec3. a00 a10 a20)
             (Vec3. a01 a11 a21)
             (Vec3. a02 a12 a22)])
  (rows [_] [(Vec3. a00 a01 a02)
             (Vec3. a10 a11 a12)
             (Vec3. a20 a21 a22)])
  (to-double-array2d [_] (fastmath.java.Array/mat2array2d a00 a01 a02
                                                          a10 a11 a12
                                                          a20 a21 a22))
  (to-float-array2d [_] (fastmath.java.Array/mat2array2d (float a00) (float a01) (float a02)
                                                         (float a10) (float a11) (float a12)
                                                         (float a20) (float a21) (float a22)))
  (to-double-array [_] (fastmath.java.Array/mat2array a00 a01 a02
                                                      a10 a11 a12
                                                      a20 a21 a22))
  (to-float-array [_] (fastmath.java.Array/mat2array (float a00) (float a01) (float a02)
                                                     (float a10) (float a11) (float a12)
                                                     (float a20) (float a21) (float a22)))

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
  (symmetric? [_ tol] (and (m/delta-eq a01 a10 tol) (m/delta-eq a12 a21 tol) (m/delta-eq a20 a02 tol)))
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
  (emulm [_ m] (let [^Mat3x3 m m]
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
                        (+ (* a02 (.x v)) (* a12 (.y v)) (* a22 (.z v))))))
  (cholesky [m] (do
                  (assert (prot/symmetric? m) "Matrix is not symmetric.")
                  (let [a (m/sqrt a00)
                        b (/ a10 a)
                        c (m/sqrt (- a11 (* b b)))
                        d (/ a20 a)
                        e (/ (- a21 (* b d)) c)
                        f (m/sqrt (- a22 (* d d) (* e e)))]
                    (Mat3x3. a 0.0 0.0 b c 0.0 d e f))))
  (norm [_ t] (if (sequential? t)
                (let [[^double p ^double q] t
                      qp (/ q p)]
                  (m/pow (+ (m/pow (+ (m/pow (m/abs a00) p) (m/pow (m/abs a10) p) (m/pow (m/abs a20) p)) qp)
                            (m/pow (+ (m/pow (m/abs a01) p) (m/pow (m/abs a11) p) (m/pow (m/abs a21) p)) qp)
                            (m/pow (+ (m/pow (m/abs a02) p) (m/pow (m/abs a12) p) (m/pow (m/abs a22) p)) qp))
                         (/ q)))
                (condp = t
                  :inf (m/max (+ (m/abs a00) (m/abs a01) (m/abs a02))
                              (+ (m/abs a10) (m/abs a11) (m/abs a12))
                              (+ (m/abs a20) (m/abs a21) (m/abs a22)))
                  :max (m/max (m/abs a00) (m/abs a01) (m/abs a02)
                              (m/abs a10) (m/abs a11) (m/abs a12)
                              (m/abs a20) (m/abs a21) (m/abs a22))
                  (m/max (+ (m/abs a00) (m/abs a10) (m/abs a20))
                         (+ (m/abs a01) (m/abs a11) (m/abs a21))
                         (+ (m/abs a02) (m/abs a12) (m/abs a22)))))))

(deftype Mat4x4 [^double a00 ^double a01 ^double a02 ^double a03
                 ^double a10 ^double a11 ^double a12 ^double a13
                 ^double a20 ^double a21 ^double a22 ^double a23
                 ^double a30 ^double a31 ^double a32 ^double a33]
  Object
  (toString [_] (str "#mat4x4 [[" a00 ", " a01 ", " a02 ", " a03 "]\n         [" a10 ", " a11 ", " a12 ", " a13 "]\n         [" a20 ", " a21 ", " a22 ", " a23 "]\n         [" a30 ", " a31 ", " a32 ", " a33 "]]"))
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
  (valAt [_ [^long x ^long y]] (when (in-range? x y 4)
                                 (gen-condition-2 4 x y nil)))
  (valAt [_ [^long x ^long y] not-found] (if (in-range? x y 4)
                                           (gen-condition-2 4 x y not-found)
                                           not-found))
  IFn
  (invoke [_ id]
    (gen-condition 4 id (mat4x4-throw-ioobe id)))
  (invoke [_ x y]
    (gen-condition-2 4 ^long x ^long y (mat4x4-throw-ioobe [x y])))
  IPersistentCollection
  (equiv [v1 v2] (.equals v1 v2))
  prot/MatrixProto
  (->seq [_] (list a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))
  (entry [_ x y] (gen-condition-2 4 ^long x ^long y (mat4x4-throw-ioobe [x y])))
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
  (to-double-array2d [_] (fastmath.java.Array/mat2array2d a00 a01 a02 a03
                                                          a10 a11 a12 a13
                                                          a20 a21 a22 a23
                                                          a30 a31 a32 a33))
  (to-float-array2d [_] (fastmath.java.Array/mat2array2d (float a00) (float a01) (float a02) (float a03)
                                                         (float a10) (float a11) (float a12) (float a13)
                                                         (float a20) (float a21) (float a22) (float a23)
                                                         (float a30) (float a31) (float a32) (float a33)))
  (to-double-array [_] (fastmath.java.Array/mat2array a00 a01 a02 a03
                                                      a10 a11 a12 a13
                                                      a20 a21 a22 a23
                                                      a30 a31 a32 a33))
  (to-float-array [_] (fastmath.java.Array/mat2array (float a00) (float a01) (float a02) (float a03)
                                                     (float a10) (float a11) (float a12) (float a13)
                                                     (float a20) (float a21) (float a22) (float a23)
                                                     (float a30) (float a31) (float a32) (float a33)))
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
  (symmetric? [_ tol] (and (m/delta-eq a01 a10 tol) (m/delta-eq a12 a21 tol)
                           (m/delta-eq a20 a02 tol) (m/delta-eq a30 a03 tol)
                           (m/delta-eq a31 a13 tol) (m/delta-eq a32 a23 tol)))
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
  (emulm [_ m] (let [^Mat4x4 m m]
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
                        (+ (* a03 (.x v)) (* a13 (.y v)) (* a23 (.z v)) (* a33 (.w v))))))
  (cholesky [m] (do
                  (assert (prot/symmetric? m) "Matrix is not symmetric.")
                  (let [a (m/sqrt a00)
                        b (/ a10 a)
                        c (m/sqrt (- a11 (* b b)))
                        d (/ a20 a)
                        e (/ (- a21 (* b d)) c)
                        f (m/sqrt (- a22 (* d d) (* e e)))
                        g (/ a30 a)
                        h (/ (- a31 (* b g)) c)
                        i (/ (- a32 (* d g) (* e h)) f)
                        j (m/sqrt (- a33 (* g g) (* h h) (* i i)))]
                    (Mat4x4. a 0.0 0.0 0.0 b c 0.0 0.0 d e f 0.0 g h i j))))
  (norm [_ t] (if (sequential? t)
                (let [[^double p ^double q] t
                      qp (/ q p)]
                  (m/pow (+ (m/pow (+ (m/pow (m/abs a00) p) (m/pow (m/abs a10) p)
                                      (m/pow (m/abs a20) p) (m/pow (m/abs a30) p)) qp)
                            (m/pow (+ (m/pow (m/abs a01) p) (m/pow (m/abs a11) p)
                                      (m/pow (m/abs a21) p) (m/pow (m/abs a31) p)) qp)
                            (m/pow (+ (m/pow (m/abs a02) p) (m/pow (m/abs a12) p)
                                      (m/pow (m/abs a22) p) (m/pow (m/abs a32) p)) qp)
                            (m/pow (+ (m/pow (m/abs a03) p) (m/pow (m/abs a13) p)
                                      (m/pow (m/abs a23) p) (m/pow (m/abs a33) p)) qp))
                         (/ q)))
                (condp = t
                  :inf (m/max (+ (m/abs a00) (m/abs a01) (m/abs a02) (m/abs a03))
                              (+ (m/abs a10) (m/abs a11) (m/abs a12) (m/abs a13))
                              (+ (m/abs a20) (m/abs a21) (m/abs a22) (m/abs a23))
                              (+ (m/abs a30) (m/abs a31) (m/abs a32) (m/abs a33)))
                  :max (m/max (m/abs a00) (m/abs a01) (m/abs a02) (m/abs a03)
                              (m/abs a10) (m/abs a11) (m/abs a12) (m/abs a13)
                              (m/abs a20) (m/abs a21) (m/abs a22) (m/abs a23)
                              (m/abs a30) (m/abs a31) (m/abs a32) (m/abs a33))
                  (m/max (+ (m/abs a00) (m/abs a10) (m/abs a20) (m/abs a30))
                         (+ (m/abs a01) (m/abs a11) (m/abs a21) (m/abs a31))
                         (+ (m/abs a02) (m/abs a12) (m/abs a22) (m/abs a32))
                         (+ (m/abs a03) (m/abs a13) (m/abs a23) (m/abs a33)))))))

(extend (Class/forName "[[D")
  prot/MatrixProto
  {:->seq (fn [arrs] (m/double-double-array->seq arrs))
   :entry (fn [^"[[D" arrs ^long x ^long y] (fastmath.java.Array/aget2d arrs x y))
   :fmap (fn [arrs f] (into-array (map #(v/fmap % f) arrs)))
   :rows seq
   :cols (fn [arrs] (fastmath.java.Array/mat2cols arrs))
   :to-double-array2d identity
   :to-float-array2d (fn [arrs] (into-array (map float-array arrs)))
   :to-double-array (fn [arrs] (double-array (mapcat seq arrs)))
   :to-float-array (fn [arrs] (float-array (mapcat seq arrs)))
   :nrow alength 
   :ncol (comp alength first)
   :row (fn [^"[[D" arrs ^long id] (aget arrs id))
   :column (fn [arrs ^long id] (fastmath.java.Array/mat2column arrs id))
   :symmetric? (fn ([^"[[D" arrs] (let [nr (prot/nrow arrs)
                                       nc (prot/ncol arrs)]
                                   (every? identity (for [^long r (range nr)
                                                          c (range (inc r) nc)]
                                                      (== ^double (aget arrs r c)
                                                          ^double (aget arrs c r))))))
                 ([^"[[D" arrs ^double tol] (let [nr (prot/nrow arrs)
                                                  nc (prot/ncol arrs)]
                                              (every? identity (for [^long r (range nr)
                                                                     c (range (inc r) nc)]
                                                                 (m/delta-eq (aget arrs r c)
                                                                             (aget arrs c r) tol))))))
   :transpose (fn [arrs] (into-array (fastmath.java.Array/mat2cols arrs)))
   :inverse (fn [^"[[D" arrs] (->  arrs
                                  (Array2DRowRealMatrix.)
                                  (MatrixUtils/inverse)
                                  (.getData)))
   :diag (fn [arrs] (fastmath.java.Array/mat2diag arrs))
   :det (fn [^"[[D" arrs] (-> arrs
                             (Array2DRowRealMatrix.)
                             (LUDecomposition.)
                             (.getDeterminant)))
   :add (fn [arrs1 arrs2] (fastmath.java.Array/matadd arrs1 arrs2))
   :adds (fn [arrs1 ^double v] (fastmath.java.Array/matadds arrs1 v))
   :sub (fn ([arrs1 arrs2] (fastmath.java.Array/matsub arrs1 arrs2))
          ([arrs1] (fastmath.java.Array/matsub arrs1)))
   :emulm (fn [arrs1 arrs2] (fastmath.java.Array/matemulm arrs1 arrs2))
   :mulm (fn ([^"[[D" arrs1 t1? ^"[[D" arrs2 t2?]
             (let [m1 (Array2DRowRealMatrix. arrs1)
                   m2 (Array2DRowRealMatrix. arrs2)
                   ^Array2DRowRealMatrix m1 (if t1? (.transpose m1) m1)
                   ^Array2DRowRealMatrix m2 (if t2? (.transpose m2) m2)]
               (.getDataRef (.multiply m1 m2))))
           ([^"[[D" arrs1 ^"[[D" arrs2]
            (let [m1 (Array2DRowRealMatrix. arrs1)
                  m2 (Array2DRowRealMatrix. arrs2)]
              (.getDataRef (.multiply m1 m2)))))
   :mulv (fn [^"[[D" arrs ^doubles v]
           (-> arrs
               (Array2DRowRealMatrix.)
               (.operate v)))
   :vtmul (fn [^"[[D" arrs ^doubles v]
            (-> arrs
                (Array2DRowRealMatrix.)
                (.preMultiply v)))
   :muls (fn [arrs1 ^double v] (fastmath.java.Array/matmuls arrs1 v))
   :trace (fn [arrs] (v/sum (fastmath.java.Array/mat2diag arrs)))
   :cholesky (fn [^"[[D" arrs] (-> arrs
                                  (Array2DRowRealMatrix.)
                                  (CholeskyDecomposition.)
                                  (.getL)
                                  (.getData)))
   :norm (fn [^"[[D" arrs t] (prot/norm (Array2DRowRealMatrix. arrs) t))})

(extend RealMatrix
  prot/MatrixProto
  {:->seq (fn [^RealMatrix m] (m/double-double-array->seq (.getData m)))
   :entry (fn ^double [^RealMatrix m ^long x ^long y] (.getEntry m x y))
   :fmap (fn [^RealMatrix m f] (Array2DRowRealMatrix. ^"[[D" (prot/fmap (.getData m) f)))
   :rows (fn [^RealMatrix m] (map (fn [idx] (.getRowVector m (unchecked-int idx)))
                                 (range (.getRowDimension m))))
   :cols (fn [^RealMatrix m] (map (fn [idx] (.getColumnVector m (unchecked-int idx)))
                                 (range (.getColumnDimension m))))
   :to-double-array2d (fn [^RealMatrix m] (.getData m))
   :to-float-array2d (fn [^RealMatrix m] (prot/to-float-array2d (.getData m)))
   :to-double-array (fn [^RealMatrix m] (prot/to-double-array (.getData m)))
   :to-float-array (fn [^RealMatrix m] (prot/to-float-array (.getData m)))
   :nrow (fn [^RealMatrix m] (.getRowDimension m))
   :ncol (fn [^RealMatrix m] (.getColumnDimension m))
   :row (fn [^RealMatrix m id] (.getRowVector m id))
   :column (fn [^RealMatrix m id] (.getColumnVector m id))
   :symmetric? (fn ([^RealMatrix m] (MatrixUtils/isSymmetric m 0.0))
                 ([^RealMatrix m ^double tol] (MatrixUtils/isSymmetric m tol)))
   :transpose (fn [^RealMatrix m] (.transpose m))
   :inverse (fn [^RealMatrix m] (MatrixUtils/inverse m))
   :diag (fn [^RealMatrix m] (let [size (.getRowDimension m)
                                  v (ArrayRealVector. size)]
                              (doseq [^int idx (range size)]
                                (.setEntry v idx (.getEntry m idx idx)))
                              v))
   :det (fn [^RealMatrix m] (.getDeterminant (LUDecomposition. m)))
   :add (fn [^RealMatrix m1 ^RealMatrix m2] (.add m1 m2))
   :adds (fn [^RealMatrix m ^double v] (.scalarAdd m v))
   :sub (fn ([^RealMatrix m1 ^RealMatrix m2] (.subtract m1 m2))
          ([^RealMatrix m] (.subtract (Array2DRowRealMatrix. (.getRowDimension m)
                                                             (.getColumnDimension m)) m)))
   :emulm (fn [^RealMatrix m1 ^RealMatrix m2] (let [m (.copy m1)]
                                               (doseq [^int r (range (.getRowDimension m))
                                                       ^int c (range (.getColumnDimension m))]
                                                 (.setEntry m r c (* (.getEntry m1 r c)
                                                                     (.getEntry m2 r c))))
                                               m))
   :mulm (fn ([^RealMatrix m1 t1? ^RealMatrix m2 t2?]
             (let [m1 (if t1? (.transpose m1) m1)
                   m2 (if t2? (.transpose m2) m2)]
               (.multiply m1 m2)))
           ([^RealMatrix m1 ^RealMatrix m2] (.multiply m1 m2)))
   :mulv (fn [^RealMatrix m1 v] (if (instance? RealVector v)
                                 (.operate m1 ^RealVector v)
                                 (.operate m1 ^doubles v)))
   :vtmul (fn [^RealMatrix m1 v] (if (instance? RealVector v)
                                  (.preMultiply m1 ^RealVector v)
                                  (.preMultiply m1 ^doubles v)))
   :muls (fn [^RealMatrix m1 ^double v] (.scalarMultiply m1 v))
   :trace (fn [^RealMatrix m] (.getTrace m))
   :cholesky (fn [^RealMatrix m] (.getL (CholeskyDecomposition. m)))
   :norm (fn [^RealMatrix m t] (if (sequential? t)
                                (let [[^double p ^double q] t]
                                  (if (and (== p 2.0) (== q 2.0))
                                    (.getFrobeniusNorm m)
                                    (let [qp (/ q p)]
                                      (-> (->> (prot/cols m)
                                               (map (fn [c] (-> (v/abs c)
                                                               (v/fmap (fn [v] (m/pow v p)))
                                                               (v/sum)
                                                               (m/pow qp)))))
                                          (v/sum)
                                          (m/pow (/ q))))))
                                (condp = t
                                  :inf (.getNorm (.transpose m))
                                  :max (reduce m/max (for [^int r (range (.getRowDimension m))
                                                           ^int c (range (.getColumnDimension m))]
                                                       (m/abs (.getEntry m r c))))
                                  (.getNorm m))))})


(defn mat2x2
  "Create 2x2 matrix.

  Arity:

  * 1 - fills matrix with given value
  * 2 - creates diagonal matrix
  * 4 - creates row ordered matrix"
  (^Mat2x2 [^double v] (Mat2x2. v v v v))
  (^Mat2x2 [^double d1 ^double d2] (Mat2x2. d1 0.0 0.0 d2))
  (^Mat2x2 [^double a00 ^double a01 ^double a10 ^double a11] (Mat2x2. a00 a01 a10 a11)))

(defn rows->mat2x2
  "Create 2x2 matrix from 2d vectors (rows)."
  ^Mat2x2 [[^double a00 ^double a01]
           [^double a10 ^double a11]]
  (Mat2x2. a00 a01 a10 a11))

(defn cols->mat2x2
  "Create 2x2 matrix from 2d vectors (columns)."
  ^Mat2x2 [[^double a00 ^double a10]
           [^double a01 ^double a11]]
  (Mat2x2. a00 a01 a10 a11))

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

(defn rows->mat3x3
  "Create 3x3 matrix from 3d vectors (rows)."
  ^Mat3x3 [[^double a00 ^double a01 ^double a02]
           [^double a10 ^double a11 ^double a12]
           [^double a20 ^double a21 ^double a22]]
  (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))

(defn cols->mat3x3
  "Create 3x3 matrix from 3d vectors (columns)."
  ^Mat3x3 [[^double a00 ^double a10 ^double a20]
           [^double a01 ^double a11 ^double a21]
           [^double a02 ^double a12 ^double a22]]
  (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))

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

(defn rows->mat4x4
  "Create 4x4 matrix from 4d vectors (rows)."
  ^Mat4x4 [[^double a00 ^double a01 ^double a02 ^double a03]
           [^double a10 ^double a11 ^double a12 ^double a13]
           [^double a20 ^double a21 ^double a22 ^double a23]
           [^double a30 ^double a31 ^double a32 ^double a33]]
  (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))

(defn cols->mat4x4
  "Create 4x4 matrix from 4d vectors (columns)."
  ^Mat4x4 [[^double a00 ^double a10 ^double a20 ^double a30]
           [^double a01 ^double a11 ^double a21 ^double a31]
           [^double a02 ^double a12 ^double a22 ^double a32]
           [^double a03 ^double a13 ^double a23 ^double a33]]
  (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))

(defn rows->RealMatrix
  "Return Apache Commons Math Array2DRowMatrix from sequence of rows"
  [rows]
  (if (= (type rows) m/double-double-array-type)
    (Array2DRowRealMatrix. ^"[[D" rows)
    (Array2DRowRealMatrix. ^"[[D" (m/seq->double-double-array (map v/vec->array rows)))))

(defn cols->RealMatrix
  "Return Apache Commons Math Array2DRowMatrix from sequence of columns"
  [cols]
  (prot/transpose (rows->RealMatrix cols)))

(defn mat
  "Create mat2x2, mat3x3 or mat4x4 or RealMatrix from rows"
  ([real-matrix-rows] (rows->RealMatrix real-matrix-rows))
  ([^double a00 ^double a01 ^double a10 ^double a11] (Mat2x2. a00 a01 a10 a11))
  ([a00 a01 a02 a10 a11 a12 a20 a21 a22] (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))
  ([a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33]
   (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33)))

(defn rows->mat
  "Create nxn matrix from nd vectors (rows)."
  ([real-matrix-rows] (rows->RealMatrix real-matrix-rows))
  ([[^double a00 ^double a01]
    [^double a10 ^double a11]]
   (Mat2x2. a00 a01 a10 a11))
  ([[^double a00 ^double a01 ^double a02]
    [^double a10 ^double a11 ^double a12]
    [^double a20 ^double a21 ^double a22]]
   (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))
  ([[^double a00 ^double a01 ^double a02 ^double a03]
    [^double a10 ^double a11 ^double a12 ^double a13]
    [^double a20 ^double a21 ^double a22 ^double a23]
    [^double a30 ^double a31 ^double a32 ^double a33]]
   (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33)))

(defn cols->mat
  "Create nxn matrix from nd vectors (columns)."
  ([real-matrix-cols] (cols->RealMatrix real-matrix-cols))
  ([[^double a00 ^double a10]
    [^double a01 ^double a11]]
   (Mat2x2. a00 a01 a10 a11))
  ([[^double a00 ^double a10 ^double a20]
    [^double a01 ^double a11 ^double a21]
    [^double a02 ^double a12 ^double a22]]
   (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))
  ([[^double a00 ^double a10 ^double a20 ^double a30]
    [^double a01 ^double a11 ^double a21 ^double a31]
    [^double a02 ^double a12 ^double a22 ^double a32]
    [^double a03 ^double a13 ^double a23 ^double a33]]
   (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33)))

(def ^{:docs "Identity matrix for given size"} eye
  [nil 1.0
   (Mat2x2. 1.0 0.0 0.0 1.0)
   (Mat3x3. 1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0)
   (Mat4x4. 1.0 0.0 0.0 0.0
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0)])

(def ^{:docs "Zero matrix for given size"} zero
  [nil 0.0
   (mat2x2 0.0)
   (mat3x3 0.0)
   (mat4x4 0.0)])

(defn diagonal
  "Create diagonal matrix"
  ([v] (case (count v)
         1 (first v)
         2 (apply mat2x2 v)
         3 (apply mat3x3 v)
         4 (apply mat4x4 v)))
  ([^double a11 ^double a22] (mat2x2 a11 a22))
  ([^double a11 ^double a22 ^double a33] (mat3x3 a11 a22 a33))
  ([^double a11 ^double a22 ^double a33 ^double a44] (mat4x4 a11 a22 a33 a44)))

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
                      (.outerProduct v1 v2))
    (let [v1 (ArrayRealVector. (double-array v1))
          v2 (ArrayRealVector. (double-array v2))]
      (.getData (.outerProduct v1 v2)))))

;;

(defn entry
  "Get entry at given row and column"
  ^double [A row col]
  (prot/entry A row col))

(defn fmap
  "Apply a function `f` to each matrix element."
  [A f] (prot/fmap A f))

(defn cols
  "Return matrix columns"
  [A] (prot/cols A))

(defn rows
  "Return matrix rows"
  [A] (prot/rows A))

(defn mat->seq
  "Return flat sequence of entries (row order)"
  [A] (prot/->seq A))

(defn mat->array2d
  "Return doubles of doubles"
  [A] (prot/to-double-array2d A))

(defn mat->float-array2d
  "Return doubles of doubles"
  [A] (prot/to-float-array2d A))

(defn mat->array
  "Return flat double array of entries (row order)"
  [A] (prot/to-double-array A))

(defn mat->float-array
  "Return flat float array of entries (row order)"
  [A] (prot/to-float-array A))

(defn mat->RealMatrix
  "Return Apache Commons Math Array2DRowMatrix from a 2x2, 3x3 or 4x4 matrix"
  [A]
  (if (instance? RealMatrix A)
    A
    (Array2DRowRealMatrix. ^"[[D" (prot/to-double-array2d A))))

(defn nrow
  "Return number of rows"
  ^long [A] (prot/nrow A))

(defn ncol
  "Return number of rows"
  ^long [A] (prot/ncol A))

(defn row
  "Return row as a vector"
  [A ^long r] (prot/row A r))

(defn col
  "Return column as a vector"
  [A ^long c] (prot/column A c))

(defn symmetric?
  "Check if matrix is symmetric"
  [A] (prot/symmetric? A))

(defn transpose
  "Transpose matrix, C=A^T"
  [A] (prot/transpose A))

(defn inverse
  "Matrix inversion.

  Returns `nil` if inversion doesn't exist."
  [m]  (prot/inverse m))

(defn diag
  "Return diagonal of the matrix as a vector."
  [A] (prot/diag A))

(defn det
  "Return determinant of the matrix."
  ^double [A] (prot/det A))

(defn add
  "Add matrices, C=A+B."
  ([A] A)
  ([A B] (prot/add A B)))

(defn adds
  "Add scalar to all matrix elements"
  [A s] (prot/adds A s))

(defn sub
  "Subract matrices, C=A-B."
  ([A] (prot/sub A))
  ([A B] (prot/sub A B)))

(defn negate
  "Negate all matrix elements, C=-A"
  [A] (prot/sub A))

(defn mulm
  "Multiply two matrices, C=AxB.

  Optionally you can request transposition of matrices."
  ([A B] (prot/mulm A B))
  ([A transposeA? B transposeB?]
   (prot/mulm A transposeA? B transposeB?)))

(defn mulmt
  "Multiply with transposed matrix, C=AxB^T"
  [A B] (prot/mulm A false B true))

(defn tmulm
  "Transpose and multiply, C=A^TxB"
  [A B] (prot/mulm A true B false))

(defn tmulmt
  "Transpose both and multiply, C=A^TxB^T"
  [A B] (prot/mulm A true B true))

(defn emulm
  "Multiply two matrices element-wise, Hadamard product, C=AoB"
  [A B] (prot/emulm A B))

(defn mulv
  "Multply matrix by vector, x=Av"
  [A v] (prot/mulv A v))

(defn muls
  "Multply matrix by a scalar, C=sA"
  [A s] (prot/muls A s))

(defn vtmul
  "Multiply transposed vector by matrix, C=v^T A"
  [A v] (prot/vtmul A v))

(defn cholesky
  "Calculate L (lower by default) triangular for where L * L^T = A.

  Checks only for symmetry, can return NaNs when A is not positive-definite."
  ([A] (prot/cholesky A))
  ([A upper?]
   (if upper?
     (prot/transpose (cholesky A))
     (cholesky A))))

(defn trace
  "Return trace of the matrix (sum of diagonal elements)"
  ^double [A] (prot/trace A))

(defn eigenvalues
  "Return complex eigenvalues for given matrix as a sequence"
  [A]
  (let [^EigenDecomposition m (-> A mat->RealMatrix (EigenDecomposition.))
        re (.getRealEigenvalues m)
        im (.getImagEigenvalues m)]
    (mapv v/vec2 re im)))

(defn singular-values
  "Returun singular values of the matrix as sqrt of eigenvalues of A^T * A matrix."
  [A]
  (->> (mulm A true A false)
       (eigenvalues)
       (map first)
       (map (fn [^double x] (m/sqrt x)))))

(defn eigenvalues-matrix
  "Return eigenvalues for given matrix as a diagonal or block diagonal matrix"
  [A]
  (->> (mat->RealMatrix A)
       (EigenDecomposition.)
       ^RealMatrix (.getD)
       (.getData)
       (m/double-double-array->seq)
       (apply rows->mat)))

(defn normalize
  "Normalize columns (or rows)"
  ([A] (normalize A false))
  ([A rows?]
   (if rows?
     (let [nA (map v/normalize (rows A))]
       (if (instance? RealMatrix A)
         (rows->RealMatrix nA)
         (apply rows->mat nA)))
     (let [nA (map v/normalize (cols A))]
       (if (instance? RealMatrix A)
         (cols->RealMatrix nA)
         (apply cols->mat nA))))))

(defn shift-rows
  "Shift rows by a value or a result of the function (mean by default)"
  ([A] (shift-rows A v/average))
  ([A shift]
   (let [sf (comp - (if (fn? shift) shift (constantly (double shift))))
         nA (map #(v/shift % (sf (v/vec->seq %))) (rows A))]
     (if (instance? RealMatrix A)
       (rows->RealMatrix nA)
       (apply rows->mat nA)))))

(defn shift-cols
  "Shift columns by a value or a result of the function (mean by default)"
  ([A] (shift-cols A v/average))
  ([A shift]
   (let [sf (comp - (if (fn? shift) shift (constantly (double shift))))
         nA (map #(v/shift % (sf (v/vec->seq %))) (cols A))]
     (if (instance? RealMatrix A)
       (cols->RealMatrix nA)
       (apply cols->mat nA)))))

(defn- default-scaler [v] (m/sqrt (/ (v/dot v v) (dec (v/size v)))))

(defn scale-rows
  "Shift rows by a value (default: sqrt(sum(x^2)/(n-1))) or a result of the function"
  ([A] (scale-rows A default-scaler))
  ([A scale]
   (let [sf (if (fn? scale) scale (constantly (double scale)))
         nA (map #(v/div % (sf %)) (rows A))]
     (if (instance? RealMatrix A)
       (rows->RealMatrix nA)
       (apply rows->mat nA)))))

(defn scale-cols
  "Shift columns by a value (default: sqrt(sum(x^2)/(n-1))) or a result of the function"
  ([A] (scale-cols A default-scaler))
  ([A scale]
   (let [sf (if (fn? scale) scale (constantly (double scale)))
         nA (map #(v/div % (sf %)) (cols A))]
     (if (instance? RealMatrix A)
       (cols->RealMatrix nA)
       (apply cols->mat nA)))))

(defn eigenvectors
  "Return eigenvectors as a matrix (columns). Vectors can be normalized."
  ([A] (eigenvectors A false))
  ([A normalize?]
   (let [evs (->> (mat->RealMatrix A)
                  (EigenDecomposition.)
                  ^RealMatrix (.getV)
                  (.getData)
                  (m/double-double-array->seq)
                  (apply rows->mat))]
     (if normalize? (normalize evs) evs))))

;;

(defn norm
  "Calculate norm of the matrix for given type, default: 1 (maximum absolute column sum).

  All norm types are:
  * 1 - maximum absolute column sum
  * :inf -  maximum absolute row sum
  * 2 - spectral norm, maximum singular value
  * :max - maximum absolute value
  * :frobenius - Frobenius norm
  * [p,q] - generalized L_pq norm, [2,2] - Frobenius norm, [p,p] - entrywise p-norm
  * [p] - Shatten p-norm, [1] - nuclear/trace norm"
  (^double [A] (norm A 1))
  (^double [A norm-type]
   (cond
     (= norm-type :frobenius) (prot/norm A [2 2])
     (= norm-type 2) (reduce m/max (singular-values A))
     (and (sequential? norm-type)
          (= 1 (count norm-type))) (let [p (double (first norm-type))]
                                     (m/pow (->> (singular-values A)
                                                 (map (fn [^double s] (m/pow s p)))
                                                 (reduce m/+)) (/ p)))
     :else (prot/norm A norm-type))))

(defn condition
  "Condition number calculated for L2 norm by default (see [[norm]] for other norm types).

   Cond(A) = norm(A) * norm(inv(A))"
  (^double [A] (condition A 2))
  (^double [A norm-type] (* (norm A norm-type) (norm (inverse A) norm-type))))

;;

(defn rotation-matrix-2d
  "Create rotation matrix for a plane"
  ^Mat2x2 [^double theta]
  (let [st (m/sin theta)
        ct (m/cos theta)]
    (Mat2x2. ct (- st) st ct)))

(defn rotation-matrix-3d
  "Create rotation matrix for a 3d space. Tait–Bryan angles z-y′-x″"
  (^Mat3x3 [[^double x ^double y ^double z]] (rotation-matrix-3d x y z))
  (^Mat3x3 [^double x ^double y ^double z]
   (let [sx (m/sin x)
         cx (m/cos x)
         sy (m/sin y)
         cy (m/cos y)
         sz (m/sin z)
         cz (m/cos z)
         sycz (* sy cz)
         sysz (* sy sz)]
     (Mat3x3.  (* cy cz) (- (* cy sz)) sy
               (+ (* cx sz) (* sx sycz)) (- (* cx cz) (* sx sysz)) (- (* sx cy))
               (- (* sx sz) (* cx sycz)) (+ (* sx cz) (* cx sysz)) (* cx cy)))))

(defn rotation-matrix-3d-x
  "Create rotation matrix for a 3d space, x-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. 1.0 0.0 0.0
             0.0 ca (- sa)
             0.0 sa ca)))

(defn rotation-matrix-3d-y
  "Create rotation matrix for a 3d space, y-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. ca 0.0 sa
             0.0 1.0 0.0
             (- sa) 0.0 ca)))

(defn rotation-matrix-3d-z
  "Create rotation matrix for a 3d space, z-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. ca (- sa) 0.0
             sa ca 0.0
             0.0 0.0 1.0)))

(defn rotation-matrix-axis-3d
  "Create 3d rotation matrix for axis ratation."
  ^Mat3x3 [^double angle ^Vec3 axis]
  (let [^Vec3 axis (v/normalize axis)
        e1 (.x axis)
        e2 (.y axis)
        e3 (.z axis)
        sa (m/sin angle)
        ca (m/cos angle)]
    (add (add (mat3x3 ca ca ca)
              (muls (outer axis axis) (- 1.0 ca)))
         (muls (mat3x3 0.0 (- e3) e2
                       e3 0.0 (- e1)
                       (- e2) e1 0.0) sa))))

(defmacro ^:private primitive-ops
  "Generate primitive functions operating on vectors"
  [fns]
  (let [v (symbol "vector")]
    `(do ~@(for [f fns
                 :let [nm (symbol (name f))
                       doc (str "Apply " nm " to matrix elements.")]]
             `(defn ~nm ~doc [~v]
                (prot/fmap ~v ~f))))))

(primitive-ops [m/sin m/cos m/tan m/asin m/acos m/atan m/sinh m/cosh m/tanh m/asinh m/acosh m/atanh
                m/cot m/sec m/csc m/acot m/asec m/acsc m/coth m/sech m/csch m/acoth m/asech m/acsch
                m/sq m/cb m/safe-sqrt m/sqrt m/cbrt m/exp m/log m/log10 m/log2 m/ln m/log1p m/expm1
                m/log1pexp m/log1mexp m/log1psq m/log1pmx m/logmxp1 m/logexpm1
                m/radians m/degrees m/sinc m/jinc m/sigmoid m/logit m/xlogx
                m/floor m/ceil m/round m/rint m/trunc m/frac m/sfrac m/signum m/sgn])

(defmethod print-method Mat2x2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat2x2 [v w] (print-method v w))

(defmethod print-method Mat3x3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat3x3 [v w] (print-method v w))

(defmethod print-method Mat4x4 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-dup Mat4x4 [v w] (print-method v w))

(m/unuse-primitive-operators #{'abs})
