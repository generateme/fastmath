(ns fastmath.matrix
  "Provides tools for working with various matrix types, including fixed-size (2x2, 3x3, 4x4), Java `double[][]` arrays, and Apache Commons Math `RealMatrix`.

  It offers efficient mathematical operations for linear algebra, geometric transformations, and data manipulation, unifying different representations under a common protocol approach where appropriate."
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.protocols.matrix :as prot])
  (:import [clojure.lang Counted IFn Seqable Reversible ILookup]
           [org.apache.commons.math3.linear Array2DRowRealMatrix RealVector ArrayRealVector RealMatrix MatrixUtils
            EigenDecomposition QRDecomposition RRQRDecomposition LUDecomposition CholeskyDecomposition SingularValueDecomposition
            DecompositionSolver]
           [fastmath.vector Vec2 Vec3 Vec4]
           [org.apache.commons.math3.stat StatUtils]
           [fastmath.java Array]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators #{'abs})

(defn- mat-throw-ioobe
  [dim id]
  (throw (IndexOutOfBoundsException. (str "Index " id " out of bounds for mat" dim))))

(def ^:private mat1x1-throw-ioobe (partial mat-throw-ioobe "1x1"))
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
  (and (number? x) (m/not-neg? ^long x) (< ^long x mx)
       (number? y) (m/not-neg? ^long y) (< ^long y mx)))

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
  prot/MatrixProto
  (->seq [_] (list a00 a01 a10 a11))
  (entry [_ x y] (gen-condition-2 2 ^long x ^long y (mat2x2-throw-ioobe [x y])))
  (fmap [_ f] (Mat2x2. (f a00) (f a01) (f a10) (f a11)))
  (cols [_] [(Vec2. a00 a10)
             (Vec2. a01 a11)])
  (rows [_] [(Vec2. a00 a01)
             (Vec2. a10 a11)])
  (to-double-array2d [_] (Array/mat2array2d a00 a01
                                            a10 a11))
  (to-float-array2d [_] (Array/mat2array2d (float a00) (float a01)
                                           (float a10) (float a11)))

  (to-double-array [_] (Array/mat2array a00 a01
                                        a10 a11))
  (to-float-array [_] (Array/mat2array (float a00) (float a01)
                                       (float a10) (float a11)))
  (to-real-matrix [_] (Array2DRowRealMatrix. (Array/mat2array2d a00 a01
                                                                a10 a11)))

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
  (singular? [m] (m/zero? (double (prot/det m))))
  (solve [m b] (prot/mulv (prot/inverse m) b))
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
  (to-double-array2d [_] (Array/mat2array2d a00 a01 a02
                                            a10 a11 a12
                                            a20 a21 a22))
  (to-float-array2d [_] (Array/mat2array2d (float a00) (float a01) (float a02)
                                           (float a10) (float a11) (float a12)
                                           (float a20) (float a21) (float a22)))
  (to-double-array [_] (Array/mat2array a00 a01 a02
                                        a10 a11 a12
                                        a20 a21 a22))
  (to-float-array [_] (Array/mat2array (float a00) (float a01) (float a02)
                                       (float a10) (float a11) (float a12)
                                       (float a20) (float a21) (float a22)))
  (to-real-matrix [_] (Array2DRowRealMatrix. (Array/mat2array2d a00 a01 a02
                                                                a10 a11 a12
                                                                a20 a21 a22)))

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
  (singular? [m] (m/zero? (double (prot/det m))))
  (solve [m b] (prot/mulv (prot/inverse m) b))
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
  (to-double-array2d [_] (Array/mat2array2d a00 a01 a02 a03
                                            a10 a11 a12 a13
                                            a20 a21 a22 a23
                                            a30 a31 a32 a33))
  (to-float-array2d [_] (Array/mat2array2d (float a00) (float a01) (float a02) (float a03)
                                           (float a10) (float a11) (float a12) (float a13)
                                           (float a20) (float a21) (float a22) (float a23)
                                           (float a30) (float a31) (float a32) (float a33)))
  (to-double-array [_] (Array/mat2array a00 a01 a02 a03
                                        a10 a11 a12 a13
                                        a20 a21 a22 a23
                                        a30 a31 a32 a33))
  (to-float-array [_] (Array/mat2array (float a00) (float a01) (float a02) (float a03)
                                       (float a10) (float a11) (float a12) (float a13)
                                       (float a20) (float a21) (float a22) (float a23)
                                       (float a30) (float a31) (float a32) (float a33)))
  (to-real-matrix [_] (Array2DRowRealMatrix. (Array/mat2array2d a00 a01 a02 a03
                                                                a10 a11 a12 a13
                                                                a20 a21 a22 a23
                                                                a30 a31 a32 a33)))
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
  (singular? [m] (m/zero? (double (prot/det m))))
  (solve [m b] (prot/mulv (prot/inverse m) b))
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
  {:->seq (fn [arrs] (mapcat identity arrs))
   :entry (fn [^"[[D" arrs ^long x ^long y] (Array/aget2d arrs x y))
   :fmap (fn [arrs f] (into-array (map #(v/fmap % f) arrs)))
   :rows seq
   :cols (fn [arrs] (Array/mat2cols arrs))
   :to-double-array2d identity
   :to-float-array2d (fn [arrs] (into-array (map float-array arrs)))
   :to-double-array (fn [arrs] (double-array (mapcat seq arrs)))
   :to-float-array (fn [arrs] (float-array (mapcat seq arrs)))
   :to-real-matrix (fn [arrs] (Array2DRowRealMatrix. ^"[[D" arrs))
   :nrow alength 
   :ncol (comp alength first)
   :row (fn [^"[[D" arrs ^long id] (aget arrs id))
   :column (fn [arrs ^long id] (Array/mat2column arrs id))
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
   :transpose (fn [arrs] (into-array (Array/mat2cols arrs)))
   :inverse (fn [^"[[D" arrs] (->  arrs
                                  (Array2DRowRealMatrix.)
                                  (MatrixUtils/inverse)
                                  (.getData)))
   :diag (fn [arrs] (Array/mat2diag arrs))
   :det (fn [^"[[D" arrs] (-> arrs
                             (Array2DRowRealMatrix.)
                             (LUDecomposition.)
                             (.getDeterminant)))
   :singular? (fn [^"[[D" arrs] (not (-> arrs
                                        (Array2DRowRealMatrix.)
                                        (LUDecomposition.)
                                        (.getSolver)
                                        (.isNonSingular))))
   :solve (fn [^"[[D" arrs ^doubles v] (-> arrs
                                          (Array2DRowRealMatrix.)
                                          (MatrixUtils/inverse)
                                          (.operate v)))
   :add (fn [arrs1 arrs2] (Array/matadd arrs1 arrs2))
   :adds (fn [arrs1 ^double v] (Array/matadds arrs1 v))
   :sub (fn ([arrs1 arrs2] (Array/matsub arrs1 arrs2))
          ([arrs1] (Array/matsub arrs1)))
   :emulm (fn [arrs1 arrs2] (Array/matemulm arrs1 arrs2))
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
   :muls (fn [arrs1 ^double v] (Array/matmuls arrs1 v))
   :trace (fn [arrs] (v/sum (Array/mat2diag arrs)))
   :cholesky (fn [^"[[D" arrs] (-> arrs
                                  (Array2DRowRealMatrix.)
                                  (CholeskyDecomposition.)
                                  (.getL)
                                  (.getData)))
   :norm (fn [^"[[D" arrs t] (prot/norm (Array2DRowRealMatrix. arrs) t))})

(extend RealMatrix
  prot/MatrixProto
  {:->seq (fn [^RealMatrix m] (mapcat identity (.getData m)))
   :entry (fn ^double [^RealMatrix m ^long x ^long y] (.getEntry m x y))
   :fmap (fn [^RealMatrix m f] (Array2DRowRealMatrix. ^"[[D" (prot/fmap (.getData m) f)))
   :rows (fn [^RealMatrix m] (mapv (fn [idx] (.getRowVector m (unchecked-int idx)))
                                  (range (.getRowDimension m))))
   :cols (fn [^RealMatrix m] (mapv (fn [idx] (.getColumnVector m (unchecked-int idx)))
                                  (range (.getColumnDimension m))))
   :to-double-array2d (fn [^RealMatrix m] (.getData m))
   :to-float-array2d (fn [^RealMatrix m] (prot/to-float-array2d (.getData m)))
   :to-double-array (fn [^RealMatrix m] (prot/to-double-array (.getData m)))
   :to-float-array (fn [^RealMatrix m] (prot/to-float-array (.getData m)))
   :to-real-matrix identity
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
   :singular? (fn [^RealMatrix m] (not (.isNonSingular (.getSolver (LUDecomposition. m)))))
   :solve (fn [^RealMatrix m ^RealVector v] (.solve (.getSolver (LUDecomposition. m)) v))
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

(extend Number
  prot/MatrixProto
  {:->seq list
   :entry (fn ^double [^double n ^long x ^long y] (if (and (m/zero? x) (m/zero? y)) n (mat1x1-throw-ioobe [x y])))
   :fmap (fn ^double [^double n f] (f n))
   :rows (fn [n] [[n]])
   :cols (fn [n] [[n]])
   :to-double-array2d (fn [n] (m/seq->double-double-array [[n]]))
   :to-float-array2d (fn [n] (into-array [(float-array [n])]))
   :to-double-array (fn [n] (double-array [n]))
   :to-float-array (fn [n] (float-array [n]))
   :to-real-matrix (fn [n] (Array2DRowRealMatrix. (m/seq->double-double-array [[n]])))
   :nrow (constantly 1)
   :ncol (constantly 1)
   :row (fn [n ^long id] (if (m/zero? id) [n] (mat1x1-throw-ioobe id)))
   :column (fn [n ^long id] (if (m/zero? id) [n] (mat1x1-throw-ioobe id)))
   :symmetric? (constantly true)
   :transpose identity
   :inverse m//
   :diag vector
   :det identity
   :singular? m/zero?
   :solve (fn [^double m ^double v] (m// v m))
   :add m/+
   :adds m/+
   :sub m/-
   :emulm m/*
   :mulm m/*
   :mulv m/*
   :vtmul m/*
   :muls m/*
   :trace identity
   :cholesky m/sqrt
   :norm (fn [n _] n)})

(defn mat2x2
  "Creates 2x2 matrix.

  Arity:

  * 1 - fills matrix with given value
  * 2 - creates diagonal matrix
  * 4 - creates row ordered matrix"
  (^Mat2x2 [^double v] (Mat2x2. v v v v))
  (^Mat2x2 [^double d1 ^double d2] (Mat2x2. d1 0.0 0.0 d2))
  (^Mat2x2 [^double a00 ^double a01 ^double a10 ^double a11] (Mat2x2. a00 a01 a10 a11)))

(defn rows->mat2x2
  "Creates 2x2 matrix from 2d vectors (rows)."
  ^Mat2x2 [[^double a00 ^double a01]
           [^double a10 ^double a11]]
  (Mat2x2. a00 a01 a10 a11))

(defn cols->mat2x2
  "Create 2x2 matrix from 2d vectors (columns)."
  ^Mat2x2 [[^double a00 ^double a10]
           [^double a01 ^double a11]]
  (Mat2x2. a00 a01 a10 a11))

(defn diag->mat2x2
  "Creates 2x2 diagonal matrix."
  (^Mat2x2 [^double d] (Mat2x2. d 0.0 0.0 d))
  (^Mat2x2 [^double d1 ^double d2] (Mat2x2. d1 0.0 0.0 d2)))

(defn mat3x3
  "Creates 3x3 matrix.

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
  "Creates 3x3 matrix from 3d vectors (rows)."
  ^Mat3x3 [[^double a00 ^double a01 ^double a02]
           [^double a10 ^double a11 ^double a12]
           [^double a20 ^double a21 ^double a22]]
  (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))

(defn cols->mat3x3
  "Creates 3x3 matrix from 3d vectors (columns)."
  ^Mat3x3 [[^double a00 ^double a10 ^double a20]
           [^double a01 ^double a11 ^double a21]
           [^double a02 ^double a12 ^double a22]]
  (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))

(defn diag->mat3x3
  "Creates 3x3 diagonal matrix."
  (^Mat3x3 [^double d] (Mat3x3. d 0.0 0.0 0.0 d 0.0 0.0 0.0 d))
  (^Mat3x3 [^double d1 ^double d2 ^double d3] (Mat3x3. d1 0.0 0.0 0.0 d2 0.0 0.0 0.0 d3)))

(defn mat4x4
  "Creates 4x4 matrix.

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
  "Creates 4x4 matrix from 4d vectors (rows)."
  ^Mat4x4 [[^double a00 ^double a01 ^double a02 ^double a03]
           [^double a10 ^double a11 ^double a12 ^double a13]
           [^double a20 ^double a21 ^double a22 ^double a23]
           [^double a30 ^double a31 ^double a32 ^double a33]]
  (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))

(defn cols->mat4x4
  "Creates 4x4 matrix from 4d vectors (columns)."
  ^Mat4x4 [[^double a00 ^double a10 ^double a20 ^double a30]
           [^double a01 ^double a11 ^double a21 ^double a31]
           [^double a02 ^double a12 ^double a22 ^double a32]
           [^double a03 ^double a13 ^double a23 ^double a33]]
  (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33))

(defn diag->mat4x4
  "Creates 4x4 diagonal matrix."
  (^Mat4x4 [^double d] (Mat4x4. d 0.0 0.0 0.0
                                0.0 d 0.0 0.0
                                0.0 0.0 d 0.0
                                0.0 0.0 0.0 d))
  (^Mat4x4 [^double d1 ^double d2 ^double d3 ^double d4] (Mat4x4. d1 0.0 0.0 0.0
                                                                  0.0 d2 0.0 0.0
                                                                  0.0 0.0 d3 0.0
                                                                  0.0 0.0 0.0 d4)))
(defn real-matrix
  "Creates Apache Commons Math Array2DRowMatrix from sequence of rows"
  [rows]
  (if (= (type rows) m/double-double-array-type)
    (Array2DRowRealMatrix. ^"[[D" rows)
    (Array2DRowRealMatrix. ^"[[D" (m/seq->double-double-array (map v/vec->array rows)))))

(defn rows->RealMatrix
  "Returns Apache Commons Math Array2DRowMatrix from sequence of rows"
  [rows]
  (if (= (type rows) m/double-double-array-type)
    (Array2DRowRealMatrix. ^"[[D" rows)
    (Array2DRowRealMatrix. ^"[[D" (m/seq->double-double-array (map v/vec->array rows)))))

(defn cols->RealMatrix
  "Returns Apache Commons Math Array2DRowMatrix from sequence of columns"
  [cols]
  (prot/transpose (rows->RealMatrix cols)))

(defn mat
  "Creates mat2x2, mat3x3 or mat4x4 or RealMatrix from rows"
  ([real-matrix-rows] (rows->RealMatrix real-matrix-rows))
  ([^double a00 ^double a01 ^double a10 ^double a11] (Mat2x2. a00 a01 a10 a11))
  ([a00 a01 a02 a10 a11 a12 a20 a21 a22] (Mat3x3. a00 a01 a02 a10 a11 a12 a20 a21 a22))
  ([a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33]
   (Mat4x4. a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33)))

(defn rows->mat
  "Creates nxn matrix from nd vectors (rows)."
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
  "Creates nxn matrix from nd vectors (columns)."
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

(defn array2d->mat2x2
  "Creates 2x2 matrix from 2d double array."
  [^"[[D" arrs]
  (Mat2x2. (Array/aget2d arrs 0 0)
           (Array/aget2d arrs 0 1)
           (Array/aget2d arrs 1 0)
           (Array/aget2d arrs 1 1)))

(defn array2d->mat3x3
  "Creates 3x3 matrix from 2d double array."
  [^"[[D" arrs]
  (Mat3x3. (Array/aget2d arrs 0 0)
           (Array/aget2d arrs 0 1)
           (Array/aget2d arrs 0 2)
           (Array/aget2d arrs 1 0)
           (Array/aget2d arrs 1 1)
           (Array/aget2d arrs 1 2)
           (Array/aget2d arrs 2 0)
           (Array/aget2d arrs 2 1)
           (Array/aget2d arrs 2 2)))

(defn array2d->mat4x4
  "Creates 4x4 matrix from 2d double array."
  [^"[[D" arrs]
  (Mat4x4. (Array/aget2d arrs 0 0)
           (Array/aget2d arrs 0 1)
           (Array/aget2d arrs 0 2)
           (Array/aget2d arrs 0 3)
           (Array/aget2d arrs 1 0)
           (Array/aget2d arrs 1 1)
           (Array/aget2d arrs 1 2)
           (Array/aget2d arrs 1 3)
           (Array/aget2d arrs 2 0)
           (Array/aget2d arrs 2 1)
           (Array/aget2d arrs 2 2)
           (Array/aget2d arrs 2 3)
           (Array/aget2d arrs 3 0)
           (Array/aget2d arrs 3 1)
           (Array/aget2d arrs 3 2)
           (Array/aget2d arrs 3 3)))

(defn array2d->RealMatrix
  "Creates RealMatrix matrix from 2d double array."
  [^"[[D" arrs]
  (Array2DRowRealMatrix. arrs))

(defn eye
  "Creates identity matrix for given size."
  ([^long size real-matrix?]
   (if real-matrix? (MatrixUtils/createRealIdentityMatrix size) (eye size)))
  ([^long size]
   (case (int size)
     2 (Mat2x2. 1.0 0.0 0.0 1.0)
     3 (Mat3x3. 1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0)
     4 (Mat4x4. 1.0 0.0 0.0 0.0
                0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0)
     (MatrixUtils/createRealIdentityMatrix size))))

(defn zero
  "Creates zero matrix for given size."
  ([^long rows ^long cols real-matrix?]
   (if (and (not real-matrix?)
            (m/== rows cols)
            (m/<= 2 rows 4)
            (m/<= 2 cols 4))
     (zero rows)
     (Array2DRowRealMatrix. rows cols)))
  ([^long size real-matrix?]
   (if real-matrix? (Array2DRowRealMatrix. size size) (zero size)))
  ([^long size]
   (case (int size)
     2 (mat2x2 0.0)
     3 (mat3x3 0.0)
     4 (mat4x4 0.0)
     (Array2DRowRealMatrix. size size))))

(defn diagonal
  "Creates diagonal matrix."
  ([v] (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array v)))
  ([^double a11 ^double a22] (mat2x2 a11 a22))
  ([^double a11 ^double a22 ^double a33] (mat3x3 a11 a22 a33))
  ([^double a11 ^double a22 ^double a33 ^double a44] (mat4x4 a11 a22 a33 a44)))

(defn outer
  "Outer project for two vectors."
  [v1 v2]
  (let [i1 (class v1)
        i2 (class v2)]
    (condp = [i1 i2]
      [Vec2 Vec2] (let [^Vec2 v1 v1
                        ^Vec2 v2 v2]
                    (Mat2x2. (* (.x v1) (.x v2)) (* (.x v1) (.y v2))
                             (* (.y v1) (.x v2)) (* (.y v1) (.y v2))))
      [Vec3 Vec3] (let [^Vec3 v1 v1
                        ^Vec3 v2 v2]
                    (Mat3x3. (* (.x v1) (.x v2)) (* (.x v1) (.y v2)) (* (.x v1) (.z v2))
                             (* (.y v1) (.x v2)) (* (.y v1) (.y v2)) (* (.y v1) (.z v2))
                             (* (.z v1) (.x v2)) (* (.z v1) (.y v2)) (* (.z v1) (.z v2))))
      [Vec4 Vec4] (let [^Vec4 v1 v1
                        ^Vec4 v2 v2]
                    (Mat4x4. (* (.x v1) (.x v2)) (* (.x v1) (.y v2)) (* (.x v1) (.z v2)) (* (.x v1) (.w v2))
                             (* (.y v1) (.x v2)) (* (.y v1) (.y v2)) (* (.y v1) (.z v2)) (* (.y v1) (.w v2))
                             (* (.z v1) (.x v2)) (* (.z v1) (.y v2)) (* (.z v1) (.z v2)) (* (.z v1) (.w v2))
                             (* (.w v1) (.x v2)) (* (.w v1) (.y v2)) (* (.w v1) (.z v2)) (* (.w v1) (.w v2))))

      (let [v1 (v/vec->RealVector v1)
            v2 (v/vec->RealVector v2)]
        (.outerProduct v1 v2)))))

(defn kronecker
  "Returns Kronecker product of two matrices."
  ^RealMatrix [mat1 mat2]
  (let [r1 (long (prot/nrow mat1))
        r2 (long (prot/nrow mat2))
        c1 (long (prot/ncol mat1))
        c2 (long (prot/ncol mat2))
        target (Array2DRowRealMatrix. (m/* r1 r2) (m/* c1 c2))]
    (doseq [^long r (range r1)
            ^long c (range c1)]
      (.setSubMatrix target (->> (prot/entry mat1 r c)
                                 (prot/muls mat2)
                                 (prot/to-double-array2d)) (m/* r r2) (m/* c c2)))
    target))

;;

(defn entry
  "Returns entry at given row and column"
  ^double [A row col]
  (prot/entry A row col))

(defn fmap
  "Applies a function `f` to each matrix element."
  [A f] (prot/fmap A f))

(defn cols
  "Returns matrix columns"
  [A] (prot/cols A))

(defn rows
  "Returns matrix rows"
  [A] (prot/rows A))

(defn mat->seq
  "Returns flat sequence of entries (row order)"
  [A] (prot/->seq A))

(defn mat->array2d
  "Returns doubles of doubles"
  [A] (prot/to-double-array2d A))

(defn mat->float-array2d
  "Returns doubles of doubles"
  [A] (prot/to-float-array2d A))

(defn mat->array
  "Returns flat double array of entries (row order)"
  [A] (prot/to-double-array A))

(defn mat->float-array
  "Returns flat float array of entries (row order)"
  [A] (prot/to-float-array A))

(defn mat->RealMatrix
  "Returns Apache Commons Math Array2DRowMatrix from a 2x2, 3x3 or 4x4 matrix"
  [A]
  (prot/to-real-matrix A))

(defn nrow
  "Returns number of rows"
  ^long [A] (prot/nrow A))

(defn ncol
  "Returns number of rows"
  ^long [A] (prot/ncol A))

(defn shape
  "Returns [nrow, ncol] pair."
  [A] [(nrow A) (ncol A)])

(defn row
  "Returns row as a vector"
  [A ^long r] (prot/row A r))

(defn col
  "Returns column as a vector"
  [A ^long c] (prot/column A c))

(defn symmetric?
  "Checks if matrix is symmetric"
  ([A] (prot/symmetric? A))
  ([A ^double tolerance] (prot/symmetric? A tolerance)))

(defn transpose
  "Transposes matrix, C=A^T"
  [A] (prot/transpose A))

(defn inverse
  "Matrix inversion.

  Returns `nil` if inversion doesn't exist."
  [m]  (prot/inverse m))

(defn diag
  "Returns diagonal of the matrix as a vector."
  [A] (prot/diag A))

(defn det
  "Returns determinant of the matrix."
  ^double [A] (prot/det A))

(defn add
  "Adds matrices, C=A+B."
  ([A] A)
  ([A B] (prot/add A B)))

(defn adds
  "Adds scalar to all matrix elements"
  [A s] (prot/adds A s))

(defn sub
  "Subracts matrices, C=A-B."
  ([A] (prot/sub A))
  ([A B] (prot/sub A B)))

(defn negate
  "Negates all matrix elements, C=-A"
  [A] (prot/sub A))

(defn mulm
  "Multiplies two matrices, C=AxB.

  Optionally you can apply transposition of matrices."
  ([A B] (prot/mulm A B))
  ([A transposeA? B transposeB?]
   (prot/mulm A transposeA? B transposeB?)))

(defn mulmt
  "Multiplies with transposed matrix, C=AxB^T"
  [A B] (prot/mulm A false B true))

(defn tmulm
  "Transposes and multiplies, C=A^TxB"
  [A B] (prot/mulm A true B false))

(defn tmulmt
  "Transposes both and multiplies, C=A^TxB^T"
  [A B] (prot/mulm A true B true))

(defn emulm
  "Multiplies two matrices element-wise, Hadamard product, C=AoB"
  [A B] (prot/emulm A B))

(defn mulv
  "Multplies matrix by a vector, x=Av"
  [A v] (prot/mulv A v))

(defn muls
  "Multplies matrix by a scalar, C=sA"
  [A s] (prot/muls A s))

(defn vtmul
  "Multiplies transposed vector by matrix, C=v^T A"
  [A v] (prot/vtmul A v))

(defn cholesky
  "Calculates L (lower by default) triangular for where L * L^T = A.

  Checks only for symmetry, can return NaNs when A is not positive-definite."
  {:deprecated "Use `cholesky-decomposition` instead."}
  ([A] (prot/cholesky A))
  ([A upper?]
   (if upper?
     (prot/transpose (cholesky A))
     (cholesky A))))

(defn trace
  "Returns trace of the matrix (sum of diagonal elements)"
  ^double [A] (prot/trace A))

(defn eigenvalues
  "Returns complex eigenvalues for given matrix as a sequence"
  [A]
  (let [^EigenDecomposition m (-> A mat->RealMatrix (EigenDecomposition.))
        re (.getRealEigenvalues m)
        im (.getImagEigenvalues m)]
    (mapv v/vec2 re im)))

(defn singular-values
  "Returuns singular values of the matrix as sqrt of eigenvalues of A^T * A matrix."
  [A]
  (->> (mulm A true A false)
       (eigenvalues)
       (map first)
       (map (fn [^double x] (m/sqrt x)))))

(defn eigenvalues-matrix
  "Returns eigenvalues for given matrix as a diagonal or block diagonal matrix"
  [A]
  (->> (mat->RealMatrix A)
       (EigenDecomposition.)
       ^RealMatrix (.getD)
       (.getData)
       (m/double-double-array->seq)
       (apply rows->mat)))

(defn square?
  "Is matrix square?"
  [m]
  (m/== (nrow m) (ncol m)))

(defn block-diagonal
  "Creates block diagonal matrix (RealMatrix) from a sequence of square matrices."
  ([m & r] (block-diagonal (conj r m)))
  ([mats]
   (if (every? square? mats)
     (let [trs (reductions m/+ 0 (map prot/nrow mats))
           tr (int (last trs))
           target (Array2DRowRealMatrix. tr tr)]
       (doseq [[^long pos m] (map vector trs mats)]
         (.setSubMatrix target (prot/to-double-array2d m) pos pos))
       target)
     (throw (ex-info "Only squared matrices can be used" {:shapes (map shape mats)})))))

(defn bind-cols
  "Creates matrix from columns of given matrices."
  ([m & r] (bind-cols (conj r m)))
  ([mats]
   (let [max-rows (int (reduce m/max (map prot/nrow mats)))
         tcs (reductions m/+ 0 (map prot/ncol mats))
         tc (int (last tcs))
         target (Array2DRowRealMatrix. max-rows tc)]
     (doseq [[^long pos m] (map vector tcs mats)]
       (.setSubMatrix target (prot/to-double-array2d m) 0 pos))
     target)))

(defn bind-rows
  "Creates matrix from rows of given matrices."
  ([m & r] (bind-rows (conj r m)))
  ([mats]
   (let [max-cols (int (reduce m/max (map prot/ncol mats)))
         trs (reductions m/+ 0 (map prot/nrow mats))
         tr (int (last trs))
         target (Array2DRowRealMatrix. tr max-cols)]
     (doseq [[^long pos m] (map vector trs mats)]
       (.setSubMatrix target (prot/to-double-array2d m) pos 0))
     target)))

(defn map-cols
  "Operate on columns, f should return a column"
  [f A]
  (let [nA (map f (cols A))]
    (if (instance? RealMatrix A)
      (cols->RealMatrix nA)
      (apply cols->mat nA))))

(defn map-rows
  "Operate on rows, f should return a row"
  [f A]
  (let [nA (map f (rows A))]
    (if (instance? RealMatrix A)
      (rows->RealMatrix nA)
      (apply rows->mat nA))))

(defn differences
  "Apply lagged differences on matrix columns, always returns RealMatrix."
  ([m] (differences m 1))
  ([m ^long diffs] (differences m diffs 1))
  ([m ^long diffs ^long lag]
   (map-cols (fn [col] (v/differences col diffs lag)) (mat->RealMatrix m))))

(defn normalize
  "Normalizes columns (or rows)"
  ([A] (normalize A false))
  ([A rows?]
   ((if rows? map-rows map-cols) v/normalize A)))

(defn demean
  "Subracts mean from columns (or rows)"
  ([A] (demean A false))
  ([A rows?]
   ((if rows? map-rows map-cols)
    (fn [v] (v/shift v (m/- (v/average v)))) A)))

(defn standardize
  "Normalizes columns (or rows) to have mean = 0 and stddev = 1"
  ([A] (standardize A false))
  ([A rows?]
   ((if rows? map-rows map-cols)
    (fn [v] (StatUtils/normalize (m/seq->double-array v))) A)))

(defn shift-rows
  "Shifts rows by a value or a result of the function (nagetive of mean by default)"
  ([A] (shift-rows A (comp m/- v/average)))
  ([A shift]
   (let [sf (if (fn? shift) shift (constantly (double shift)))]
     (map-rows (fn [v] (v/shift v (sf (v/vec->seq v)))) A))))

(defn shift-cols
  "Shifts columns by a value or a result of the function (negative of  mean by default)"
  ([A] (shift-cols A (comp m/- v/average)))
  ([A shift]
   (let [sf (comp - (if (fn? shift) shift (constantly (double shift))))]
     (map-cols (fn [v] (v/shift v (sf (v/vec->seq v)))) A))))

(defn- default-scaler [v] (m// (m/sqrt (/ (v/dot v v) (dec (v/size v))))))

(defn scale-rows
  "Multiplies rows by a value (default: 1/(sqrt(sum(x^2)/(n-1)))) or a result of the function"
  ([A] (scale-rows A default-scaler))
  ([A scale]
   (let [sf (if (fn? scale) scale (constantly (double scale)))]
     (map-rows (fn [v] (v/mult v (sf v))) A))))

(defn scale-cols
  "Multiplies cols by a value (default: 1/(sqrt(sum(x^2)/(n-1)))) or a result of the function"
  ([A] (scale-cols A default-scaler))
  ([A scale]
   (let [sf (if (fn? scale) scale (constantly (double scale)))]
     (map-cols (fn [v] (v/mult v (sf v))) A))))

;;

(defn eigenvectors
  "Returns eigenvectors as a matrix (columns). Vectors can be normalized."
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
  "Calculates norm of the matrix for given type, default: 1 (maximum absolute column sum).

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
  "Creates rotation matrix for a plane"
  ^Mat2x2 [^double theta]
  (let [st (m/sin theta)
        ct (m/cos theta)]
    (Mat2x2. ct (- st) st ct)))

(defn rotation-matrix-3d
  "Creates rotation matrix for a 3d space. Tait–Bryan angles z-y′-x″"
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
  "Creates rotation matrix for a 3d space, x-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. 1.0 0.0 0.0
             0.0 ca (- sa)
             0.0 sa ca)))

(defn rotation-matrix-3d-y
  "Creates rotation matrix for a 3d space, y-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. ca 0.0 sa
             0.0 1.0 0.0
             (- sa) 0.0 ca)))

(defn rotation-matrix-3d-z
  "Creates rotation matrix for a 3d space, z-axis, right hand rule."
  ^Mat3x3 [^double a]
  (let [sa (m/sin a)
        ca (m/cos a)]
    (Mat3x3. ca (- sa) 0.0
             sa ca 0.0
             0.0 0.0 1.0)))

(defn rotation-matrix-axis-3d
  "Creates 3d rotation matrix for axis ratation."
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

;; Matrix decomposition

(defn- ->mat
  [^long s m]
  (case s
    2 (array2d->mat2x2 (prot/to-double-array2d m))
    3 (array2d->mat3x3 (prot/to-double-array2d m))
    4 (array2d->mat4x4 (prot/to-double-array2d m))
    m))

(defn- ->vec
  [^long s v]
  (case s
    2 (v/array->vec2 (v/vec->array v))
    3 (v/array->vec3 (v/vec->array v))
    4 (v/array->vec4 (v/vec->array v))
    v))

(defn- ->mat-size
  ^long [m]
  (condp instance? m
    Mat2x2 2
    Mat3x3 3
    Mat4x4 4
    0))

(defrecord MatrixDecomposition [source components ^DecompositionSolver solver singular? ^int s]
  prot/MatrixDecompositionProto
  (component [_ c]
    (let [component (components c)]
      (if (delay? component) (deref component) component)))
  prot/MatrixProto
  (solve [_ v]
    (let [rv (v/vec->RealVector v)]
      (->vec s (.solve ^DecompositionSolver solver rv))))
  (singular? [_] singular?)
  (inverse [_]
    (->mat s (.getInverse ^DecompositionSolver solver))))

(defn qr-decomposition
  "Performs QR decomposition.

  A = Q x R, Q is orthogonal (QT x Q = I), R is upper triangular.
  
  Solver minimizes using least squares method.

  Components, access with `decomposition-component` function:

  * `:H` (Hauseholder reflecors), `:Q`, `:QT`, `:R`  - matrices.

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  ([mat] (qr-decomposition mat 0.0))
  ([mat ^double threshold]
   (let [s (->mat-size mat)
         ^QRDecomposition qr (QRDecomposition. (mat->RealMatrix mat) threshold)
         ^DecompositionSolver solver (.getSolver qr)]
     (->MatrixDecomposition qr
                            {:H (delay (->mat s (.getH qr)))
                             :Q (delay (->mat s (.getQ qr)))
                             :QT (delay (->mat s (.getQT qr)))
                             :R (delay (->mat s (.getR qr)))}
                            solver
                            (not (.isNonSingular solver)) s))))

(defn rrqr-decomposition
  "Performs Rank-Revealing QR decomposition.

  A = Q x R x inv(P), Q is orthogonal (QT x Q = I), R is upper triangular.
  
  Solver minimizes using least squares method.

  Components, access with `decomposition-component` function:

  * `:H` (Hauseholder reflecors), `:Q`, `:QT`, `:R`, `:P` (permutation)  - matrices
  * `:rank-fn` - calculates numerical matrix rank, accepts threshold for rank computation (default: 0.0).

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  ([mat] (rrqr-decomposition mat 0.0))
  ([mat ^double threshold]
   (let [s (->mat-size mat)
         ^RRQRDecomposition rrqr (RRQRDecomposition. (mat->RealMatrix mat) threshold)
         ^DecompositionSolver solver (.getSolver rrqr)]
     (->MatrixDecomposition rrqr
                            {:H (delay (->mat s (.getH rrqr)))
                             :Q (delay (->mat s (.getQ rrqr)))
                             :QT (delay (->mat s (.getQT rrqr)))
                             :R (delay (->mat s (.getR rrqr)))
                             :P (delay (->mat s (.getP rrqr)))
                             :rank-fn (fn ([^double drop-threshold] (.getRank rrqr drop-threshold))
                                        ([] (.getRank rrqr 0.0)))}
                            solver
                            (not (.isNonSingular solver)) s))))

(defn cholesky-decomposition
  "Performs Cholesky decomposition.

  Decomposition of real symmetric positive-definite matrix. A = L x LT
  
  Solver minimizes using least squares method.

  Components, access with `decomposition-component` function:

  * `:L`, `:LT` - matrices
  * `:det` - determinant.

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  ([mat] (cholesky-decomposition mat CholeskyDecomposition/DEFAULT_RELATIVE_SYMMETRY_THRESHOLD CholeskyDecomposition/DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD))
  ([mat ^double symmetry_threshold positivity_threshold]
   (let [s (->mat-size mat)
         ^CholeskyDecomposition chol (CholeskyDecomposition. (mat->RealMatrix mat) symmetry_threshold positivity_threshold)
         ^DecompositionSolver solver (.getSolver chol)]
     (->MatrixDecomposition chol
                            {:L (delay (->mat s (.getL chol)))
                             :LT (delay (->mat s (.getLT chol)))
                             :det (delay (.getDeterminant chol))}
                            solver
                            (not (.isNonSingular solver)) s))))

(defn sv-decomposition
  "Performs Singular Value Decomposition (SVD).

  A = U x S x VT

  Solver minimizes using least squares method.

  Components, access with `decomposition-component` function:

  * `:S`, `:U`, `:UT`, `:V` `:VT` - matrices
  * `:singular-values` - vector of singular values
  * `:info`
      - `:condition-number`
      - `:inv-condition-number`
      - `:norm` - L2 norm
      - `:rank` - effective numerical rank
  * `coviariance-fn` - convariance function, returns covariance V x J x VT, where J is inverse of squares of singular values. When `min-sv` argument is provided, ignores singular values lower than its value.

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  [mat]
  (let [s (->mat-size mat)
        ^SingularValueDecomposition svd (SingularValueDecomposition. (mat->RealMatrix mat))
        ^DecompositionSolver solver (.getSolver svd)]
    (->MatrixDecomposition svd
                           {:S (delay (->mat s (.getS svd)))
                            :U (delay (->mat s (.getU svd)))
                            :UT (delay (->mat s (.getUT svd)))
                            :V (delay (->mat s (.getV svd)))
                            :VT (delay (->mat s (.getVT svd)))
                            :singular-values (->vec s (.getSingularValues svd))
                            :info (delay {:condition-number (.getConditionNumber svd)
                                          :inv-condition-number (.getInverseConditionNumber svd)
                                          :norm (.getNorm svd)
                                          :rank (.getRank svd)})
                            :covariance-fn (fn ([^double min-sv] (->mat s (.getCovariance svd min-sv)))
                                             ([] (->mat s (.getCovariance svd 0.0))))}
                           solver
                           (not (.isNonSingular solver)) s)))

(defn lu-decomposition
  "Performs QR decomposition.

  A = inv(P) x L x U, L is lower triangular, U is upper triangular.
  
  Solver is exact.

  Components, access with `decomposition-component` function:

  * `:L`, `:U`, `:P` (permutation)  - matrices.
  * `:det` - determinant
  * `:pivot` - pivot permutation vector

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  ([mat] (lu-decomposition mat 1.0e-11))
  ([mat ^double threshold]
   (let [s (->mat-size mat)
         ^LUDecomposition lu (LUDecomposition. (mat->RealMatrix mat) threshold)
         ^DecompositionSolver solver (.getSolver lu)]
     (->MatrixDecomposition lu
                            {:L (delay (->mat s (.getL lu)))
                             :P (delay (->mat s (.getP lu)))
                             :U (delay (->mat s (.getU lu)))
                             :det (delay (.getDeterminant lu))
                             :pivot (delay (->vec s (double-array (.getPivot lu))))}
                            solver
                            (not (.isNonSingular solver)) s))))

(defn eigen-decomposition
  "Performs Eigen decomposition.

  A = V x D x VT, D contains eigenvalues (diagonal: real values, subdiagonal: imaginary), V - eigenvectors.

  Solver is exact.

  Components, access with `decomposition-component` function:

  * `:D`, `:V`, `:VT`, `:sqrt`  - matrices.
  * `:det` - determinant
  * `:real-eigenvalues`, `imag-eigenvalues` - eigenvalues
  * `:eigenvectors` - sequence of eigenvectors
  * `:complex?` - are eigenvalues complex?

  Can be used as input for `solve`, `inverse` and `singular?` functions."
  [mat]
  (let [s (->mat-size mat)
        ^EigenDecomposition eigen (EigenDecomposition. (mat->RealMatrix mat))
        complex? (.hasComplexEigenvalues eigen)
        ^DecompositionSolver solver (when-not complex? (.getSolver eigen))]
    (->MatrixDecomposition eigen
                           {:D (delay (->mat s (.getD eigen)))
                            :V (delay (->mat s (.getV eigen)))
                            :VT (delay (->mat s (.getVT eigen)))
                            :real-eigenvalues (delay (->vec s (.getRealEigenvalues eigen)))
                            :imag-eigenvalues (delay (->vec s (.getImagEigenvalues eigen)))
                            :sqrt (delay (->mat s (.getSquareRoot eigen)))
                            :det (delay (.getDeterminant eigen))
                            :complex? complex?
                            :eigenvectors (delay (mapv #(->vec s (.getEigenvector eigen %)) (range s)))}
                           solver
                           (when-not complex? (not (.isNonSingular solver))) s)))

(defn singular?
  "Returns singularity of the matrix"
  [mat]
  (prot/singular? mat))

(defn decomposition-component
  "Returns value of the component of matrix decomposition"
  [decomposition-matrix component-name]
  (prot/component decomposition-matrix component-name))

(defn solve
  "Solve linear equation Ax=b, if decomposition is provided, decomposition is used.

  Some decompositions solve using least squares method."
  [A b]
  (prot/solve A b))

(defmacro ^:private primitive-ops
  "Generate primitive functions operating on vectors"
  [fns]
  (let [v (symbol "vector")]
    `(do ~@(for [f fns
                 :let [nm (symbol (name f))
                       doc (str "Applies " nm " to matrix elements.")]]
             `(defn ~nm ~doc [~v]
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
  [m ^double exponent]
  (fmap m (fn [^double x] (m/pow x exponent))))

(defmethod print-method Mat2x2 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method Mat3x3 [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method Mat4x4 [v ^java.io.Writer w] (.write w (str v)))

(m/unuse-primitive-operators #{'abs})

