;; https://web.archive.org/web/20170705123142/http://www.lce.hut.fi/~ssarkka/pub/quat.pdf

(ns fastmath.quaternion
  (:refer-clojure :exclude [vector zero?])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat3x3]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ZERO (v/vec4))
(def ONE (Vec4. 1.0 0.0 0.0 0.0))
(def I (Vec4. 0.0 1.0 0.0 0.0))
(def J (Vec4. 0.0 0.0 1.0 0.0))
(def K (Vec4. 0.0 0.0 0.0 1.0))
(def -I (Vec4. 0.0 -1.0 0.0 0.0))
(def -J (Vec4. 0.0 0.0 -1.0 0.0))
(def -K (Vec4. 0.0 0.0 0.0 -1.0))

(defn quaternion
  "Create quaternion from individual values or scalar and vector parts, reprezented as `Vec4`."
  (^Vec4 [^double a ^double b ^double c ^double d] (Vec4. a b c d))
  (^Vec4 [^double scalar [^double i ^double j ^double k]] (Vec4. scalar i j k))
  (^Vec4 [^double a] (Vec4. a 0.0 0.0 0.0)))

(defn complex->quaternion
  "Create quaternion from complex number"
  ^Vec4 [^Vec2 z]
  (Vec4. (.x z) (.y z) 0.0 0.0))

(defn scalar
  "Returns scalar part of quaternion, double"
  ^double [^Vec4 quaternion] (.x quaternion))

(defn re
  "Returns scalar part of quaternion"
  ^double [^Vec4 quaternion] (.x quaternion))

(defn vector
  "Returns vector part of quaternion, `Vec3` type"
  ^Vec3 [^Vec4 quaternion] (Vec3. (.y quaternion)
                                  (.z quaternion)
                                  (.w quaternion)))

(defn im-i "Return i imaginary part" ^double [^Vec4 quaternion] (.y quaternion))
(defn im-j "Return j imaginary part" ^double [^Vec4 quaternion] (.z quaternion))
(defn im-k "Return k imaginary part" ^double [^Vec4 quaternion] (.w quaternion))

(defn real? "Is q is a real number?" [quaternion] (v/is-zero? (vector quaternion)))
(defn imaginary? "Is q is a pure imaginary number?"  [^Vec4 quaternion] (m/zero? (.x quaternion)))
(defn zero? "Is zero?" [^Vec4 quaternion] (v/is-zero? quaternion))
(defn inf? "Is infinitive?" [^Vec4 quaternion] (or (m/inf? (.x quaternion))
                                                (m/inf? (.y quaternion))
                                                (m/inf? (.z quaternion))
                                                (m/inf? (.w quaternion))))
(defn nan? "Is NaN?" [^Vec4 quaternion] (or (m/nan? (.x quaternion))
                                         (m/nan? (.y quaternion))
                                         (m/nan? (.z quaternion))
                                         (m/nan? (.w quaternion))))

(defn delta-eq
  "Compare quaternions with given accuracy (10e-6 by default)"
  ([q1 q2]
   (v/delta-eq q1 q2))
  ([q1 q2 ^double accuracy]
   (v/delta-eq q1 q2 accuracy)))

(defn arg
  "Argument of quaternion, atan2(|vector(q)|, re(q))"
  ^double [^Vec4 quaternion]
  (m/atan2 (v/mag (vector quaternion)) (.x quaternion)))

(defn norm
  "Norm of the quaternion, length of the vector"
  ^Vec4 [quaternion] (v/mag quaternion))

(defn normalize
  "Normalize quaternion"
  ^Vec4 [quaternion] (v/normalize quaternion))

(defn add
  "Sum of two quaternions"
  ^Vec4 [q1 q2] (v/add q1 q2))

(defn adds
  "Adds scalar to a quaternion"
  ^Vec4 [q1 ^double s] (v/add q1 (Vec4. s 0.0 0.0 0.0)))

(defn sub
  "Difference of two quaternions"
  ^Vec4 [q1 q2] (v/sub q1 q2))

(defn scale
  "Scale the quaternion"
  ^Vec4 [quaternion ^double scale] (v/mult quaternion scale))

(defn conjugate
  "Returns conjugate of quaternion"
  ^Vec4 [^Vec4 quaternion]
  (Vec4. (.x quaternion)
         (- (.y quaternion))
         (- (.z quaternion))
         (- (.w quaternion))))

(defn reciprocal
  ^Vec4 [^Vec4 quaternion]
  (v/mult (conjugate quaternion) (m// (v/magsq quaternion))))

#_(defn mult
    ^Vec4 [^Vec4 q1 ^Vec4 q2]
    (let [a1 (.x q1) b1 (.y q1) c1 (.z q1) d1 (.w q1)
          a2 (.x q2) b2 (.y q2) c2 (.z q2) d2 (.w q2)]
      (Vec4. (m/- (m/* a1 a2) (m/* b1 b2) (m/* c1 c2) (m/* d1 d2))
             (m/- (m/+ (m/* a1 b2) (m/* b1 a2) (m/* c1 d2)) (m/* d1 c2))
             (m/+ (m/- (m/* a1 c2) (m/* b1 d2)) (m/* c1 a2) (m/* d1 b2))
             (m/+ (m/- (m/+ (m/* a1 d2) (m/* b1 c2)) (m/* c1 b2)) (m/* d1 a2)))))

(defn mult
  "Multiply two quaternions."
  ^Vec4 [^Vec4 q1 ^Vec4 q2]
  (let [r1 (.x q1)
        r2 (.x q2)
        v1 (vector q1)
        v2 (vector q2)]
    (quaternion (m/- (m/* r1 r2) (v/dot v1 v2))
                (v/add (v/add (v/mult v2 r1)
                              (v/mult v1 r2))
                       (v/cross v1 v2)))))

(defn div
  "Divide two quaternions"
  ^Vec4 [q1 q2]
  (mult q1 (reciprocal q2)))

(defn neg
  "Negation of quaternion."
  ^Vec4 [quaternion] (v/sub quaternion))

(defn sq
  "Square of quaternion."
  ^Vec4 [^Vec4 quaternion]
  (let [a (.x quaternion)
        b (.y quaternion)
        c (.z quaternion)
        d (.w quaternion)
        aa (* 2.0 a)]
    (Vec4. (m/- (* a a) (* b b) (* c c) (* d d))
           (* aa b)
           (* aa c)
           (* aa d))))

(defn qsgn
  "sgn of the quaternion.

  Returns `0` for `0+0i+0j+0k` or calls `m/sgn` on real part otherwise."
  (^double [^double re ^double im-i ^double im-j ^double im-k]
   (if (and (m/zero? re)
            (m/zero? im-i)
            (m/zero? im-j)
            (m/zero? im-k)) 0.0 (m/sgn re)))
  (^double [^Vec4 q]
   (if (zero? q) 0.0 (m/sgn (.x q)))))


(defmacro ^:private gen-from-complex
  [sym]
  (let [src (symbol "fastmath.complex" (str sym))
        q (with-meta (symbol "q") {:tag 'fastmath.vector.Vec4})
        z (with-meta (symbol "z") {:tag 'fastmath.vector.Vec2})]
    `(defn ~sym [~q]
       (let [a# (.x ~q)
             b# (.y ~q)
             c# (.z ~q)
             d# (.w ~q)
             absim# (m/sqrt (m/+ (m/* b# b#)
                                 (m/* c# c#)
                                 (m/* d# d#)))
             ~z (~src (c/complex a# absim#))
             multpl# (if (m/zero? absim#)
                       (.y ~z)
                       (/ (.y ~z) absim#))]
         (quaternion (.x ~z) (m/* b# multpl#) (m/* c# multpl#) (m/* d# multpl#))))))

(gen-from-complex sqrt)
(gen-from-complex exp)
(gen-from-complex log)

(defn logb
  "log with base b"
  ^Vec2 [quaternion b]
  (div (log quaternion) (log b)))

(defn pow
  "Quaternion power"
  ^Vec4 [^Vec4 q ^Vec4 p]
  (exp (mult (log q) p)))

(defn rotation-quaternion
  "Create rotation quaternion around vector u and angle alpha"
  ^Vec4 [^double angle u]
  (let [half (* 0.5 angle)]
    (quaternion (m/cos half)
                (v/mult (v/normalize u) (m/sin half)))))

(defn rotate
  "Rotate 3d `in` vector around axis `u`, the same as `fastmath.vector/axis-rotate`."
  (^Vec3 [in ^Vec4 rotq]
   (let [in (apply v/vec3 in)
         qw (.x rotq)
         q (vector rotq)
         t (v/mult (v/cross q in) 2.0)]
     (v/add (v/add in (v/mult t qw)) (v/cross q t))))
  (^Vec3 [^Vec3 in ^double angle ^Vec3 u]
   (rotate in (rotation-quaternion angle u))))

(defn slerp
  "Interpolate quaternions"
  ^Vec4 [^Vec4 q1 ^Vec4 q2 ^double t]
  (let [t (m/constrain t 0.0 1.0)]
    (cond
      (m/zero? t) q1
      (m/one? t) q2
      :else (mult (exp (v/mult (log (div q2 q1)) t)) q1))))

;; https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Euler_angles_(in_3-2-1_sequence)_to_quaternion_conversion
(defn to-euler
  "Convert quaternion to Euler ZYX (body 3-2-1). Quaternion will be normalized before calculations.

  Output will contain roll (x), pitch (y) and yaw (z) angles."
  ^Vec3 [^Vec4 q]
  (let [^Vec4 q (v/normalize q)
        w (.x q) x (.y q) y (.z q) z (.w q)
        sinr-cosp (m/* 2.0 (m/+ (m/* w x) (m/* y z)))
        cosr-cosp (m/- 1.0 (m/* 2.0 (m/+ (m/* x x) (m/* y y))))
        p (m/* 2.0 (m/- (m/* w y) (m/* x z)))
        siny-cosp (m/* 2.0 (m/+ (m/* w z) (m/* x y)))
        cosy-cosp (m/- 1.0 (m/* 2.0 (m/+ (m/* y y) (m/* z z))))]
    (Vec3. (m/atan2 sinr-cosp cosr-cosp)
           (if (>= (m/abs p) 1.0)
             (m/copy-sign m/HALF_PI p)
             (m/asin p))
           (m/atan2 siny-cosp cosy-cosp))))

(defn from-euler
  "Convert Euler ZYX (body 3-2-1) representation to quaternion

  Input should be 3d vector contating roll (x), pitch (y) and yaw (z) angles, or individual values.

  * roll and yaw should be from `[-pi, pi]` range
  * pitch should be from `[-pi/2, pi/2]` range"
  (^Vec4 [[^double roll ^double pitch ^double yaw]] (from-euler roll pitch yaw))
  (^Vec4 [^double roll ^double pitch ^double yaw]
   (let [hr (m/* roll 0.5)
         hp (m/* pitch 0.5)
         hy (m/* yaw 0.5)
         cr (m/cos hr)
         sr (m/sin hr)
         cp (m/cos hp)
         sp (m/sin hp)
         cy (m/cos hy)
         sy (m/sin hy)]
     (Vec4. (m/+ (m/* cr cp cy) (m/* sr sp sy))
            (m/- (m/* sr cp cy) (m/* cr sp sy))
            (m/+ (m/* cr sp cy) (m/* sr cp sy))
            (m/- (m/* cr cp sy) (m/* sr sp cy))))))

;; OpenGL representation
;; https://ntrs.nasa.gov/api/citations/19770024290/downloads/19770024290.pdf
(defn to-angles
  "Convert quaternion to Tait–Bryan angles, z-y′-x\"."
  ^Vec3 [^Vec4 q]
  (let [^Vec4 q (v/normalize q)
        w (.x q) x (.y q) y (.z q) z (.w q)
        ww (m/* w w) xx (m/* x x) yy (m/* y y) zz (m/* z z)
        y1 (m/* 2.0 (m/- (m/* w x) (m/* y z)))
        x1 (m/+ (m/- ww xx yy) zz)
        a2 (m/* 2.0 (m/+ (m/* w y) (m/* x z)))
        y3 (m/* 2.0 (m/- (m/* w z) (m/* x y)))
        x3 (m/+ ww (m/- xx yy zz))]
    (Vec3. (m/atan2 y1 x1)
           (if (>= (m/abs a2) 1.0)
             (m/copy-sign m/HALF_PI a2)
             (m/asin a2))
           (m/atan2 y3 x3))))

(defn from-angles
  "Convert Tait–Bryan angles z-y′-x\" to quaternion."
  (^Vec4 [[^double x ^double y ^double z]] (from-angles x y z))
  (^Vec4 [^double x ^double y ^double z]
   (let [hx (m/* x 0.5)
         hy (m/* y 0.5)
         hz (m/* z 0.5)
         cx (m/cos hx)
         sx (m/sin hx)
         cy (m/cos hy)
         sy (m/sin hy)
         cz (m/cos hz)
         sz (m/sin hz)]
     (Vec4. (m/- (m/* cx cy cz) (m/* sx sy sz))
            (m/+ (m/* sx cy cz) (m/* cx sy sz))
            (m/- (m/* cx sy cz) (m/* sx cy sz))
            (m/+ (m/* sx sy cz) (m/* cx cy sz))))))

(defn to-rotation-matrix
  "Convert quaternion to rotation 3x3 matrix"
  ^Mat3x3 [^Vec4 q]
  (let [^Vec4 q (v/normalize q)
        a (.x q) b (.y q) c (.z q) d (.w q)
        aa (m/* a a) bb (m/* b b) cc (m/* c c) dd (m/* d d)
        bc (m/* b c) ad (m/* a d) bd (m/* b d)
        ac (m/* a c) cd (m/* c d) ab (m/* a b)]
    (Mat3x3. (m/+ aa (m/- bb cc dd)) (m/* 2.0 (m/- bc ad)) (m/* 2.0 (m/+ bd ac))
             (m/* 2.0 (m/+ bc ad)) (m/+ cc (m/- aa bb dd)) (m/* 2.0 (m/- cd ab))
             (m/* 2.0 (m/- bd ac)) (m/* 2.0 (m/+ cd ab)) (m/+ dd (m/- aa bb cc)))))

;; https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
(defn from-rotation-matrix
  "Convert rotation 3x3 matrix to a quaternion"
  ^Vec4 [^Mat3x3 m]
  (let [[^double t q]
        (if (m/neg? (.a22 m))
          (if (m/> (.a00 m) (.a11 m))
            (let [t (m/inc (m/- (.a00 m) (.a11 m) (.a22 m)))]
              [t (Vec4. (m/- (.a21 m) (.a12 m)) t (m/+ (.a10 m) (.a01 m)) (m/+ (.a20 m) (.a02 m)))])
            (let [t (m/inc (m/- (.a11 m) (.a00 m) (.a22 m)))]
              [t (Vec4. (m/- (.a02 m) (.a20 m)) (m/+ (.a10 m) (.a01 m)) t (m/+ (.a21 m) (.a12 m)))]))
          (if (m/< (.a00 m) (- (.a11 m)))
            (let [t (m/inc (m/- (.a22 m) (.a00 m) (.a11 m)))]
              [t (Vec4. (m/- (.a10 m) (.a01 m)) (m/+ (.a02 m) (.a20 m)) (m/+ (.a21 m) (.a12 m)) t)])
            (let [t (m/inc (m/+ (.a00 m) (.a11 m) (.a22 m)))]
              [t (Vec4. t (m/- (.a21 m) (.a12 m)) (m/- (.a02 m) (.a20 m)) (m/- (.a10 m) (.a01 m)))])))]
    (v/mult q (m// 0.5 (m/sqrt t)))))

;; trig

;; https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/Quaternion.cpp

(gen-from-complex sin)
(gen-from-complex cos)
(gen-from-complex tan)
(gen-from-complex sec)
(gen-from-complex csc)
(gen-from-complex cot)
(gen-from-complex sinh)
(gen-from-complex cosh)
(gen-from-complex tanh)
(gen-from-complex sech)
(gen-from-complex csch)
(gen-from-complex coth)
(gen-from-complex asin)
(gen-from-complex acos)
(gen-from-complex atan)
(gen-from-complex asec)
(gen-from-complex acsc)
(gen-from-complex acot)
(gen-from-complex asinh)
(gen-from-complex acosh)
(gen-from-complex atanh)
(gen-from-complex asech)
(gen-from-complex acsch)
(gen-from-complex acoth)

;; https://math.stackexchange.com/questions/1499095/how-to-calculate-sin-cos-tan-of-a-quaternion

