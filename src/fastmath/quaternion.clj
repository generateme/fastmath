;; https://web.archive.org/web/20170705123142/http://www.lce.hut.fi/~ssarkka/pub/quat.pdf

(ns fastmath.quaternion
  "Operations for quaternions.

  Quaternions extend complex numbers and are widely used in fields like 3D graphics and physics for representing rotations.

  In `fastmath`, quaternions are represented as 4-dimensional vectors ([[Vec4]]) where the components correspond to the scalar part and the three imaginary parts (`i`, `j`, `k`): `q = a + bi + cj + dk` is `(Vec4. a b c d)`.

  The namespace provides functions for creating quaternions, accessing scalar and vector parts, predicates (e.g., [[real?]], [[zero?]], [[inf?]], [[nan?]]), and fundamental properties (magnitude, argument, normalization).

  A comprehensive set of operations is included:
  - **Arithmetic:** Addition, subtraction, multiplication (not commutative), division, negation, square, reciprocal, scaling, conjugation, dot product. Functions accept numbers and complex numbers where a quaternion is expected, `add`, `sub` and `mult` accept any number of arguments.
  - **Transcendental Functions:** Extensions of standard complex functions like exponential, logarithm, power, trigonometric, hyperbolic functions, and their inverses.
  - **Rotations:** Functions for creating rotation quaternions (from an axis and an angle, from two vectors), rotating 3D vectors ([[rotate]]), interpolation ([[slerp]], [[nlerp]]), and conversions between quaternions, Euler angles (ZYX body 3-2-1 and z-y'-x''), axis and angle, and rotation matrices.

  The transcendental functions are computed with the complex functions of `fastmath.complex` applied to the scalar part and the magnitude of the vector part, so they share their branch cuts and accuracy. A quaternion with a zero vector part is treated as a complex number with the imaginary part on the `i` axis (for example the square root of `-4` is `2i`).

  The implementation correctly handles floating-point special values, including `##Inf` and `##NaN`."
  (:refer-clojure :exclude [vector zero?])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat3x3]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^{:doc "The zero quaternion, 0+0i+0j+0k."} ZERO (v/vec4))
(def ^{:doc "The real unit quaternion, 1+0i+0j+0k."} ONE (Vec4. 1.0 0.0 0.0 0.0))
(def ^{:doc "The imaginary unit i, 0+1i+0j+0k."} I (Vec4. 0.0 1.0 0.0 0.0))
(def ^{:doc "The imaginary unit j, 0+0i+1j+0k."} J (Vec4. 0.0 0.0 1.0 0.0))
(def ^{:doc "The imaginary unit k, 0+0i+0j+1k."} K (Vec4. 0.0 0.0 0.0 1.0))
(def ^{:doc "Minus the unit i, 0-1i+0j+0k."} -I (Vec4. 0.0 -1.0 0.0 0.0))
(def ^{:doc "Minus the unit j, 0+0i-1j+0k."} -J (Vec4. 0.0 0.0 -1.0 0.0))
(def ^{:doc "Minus the unit k, 0+0i+0j-1k."} -K (Vec4. 0.0 0.0 0.0 -1.0))

(defn quaternion
  "Creates a quaternion represented as [[Vec4]].

  Parameters:

  - `a`, `b`, `c`, `d`: the scalar part and the `i`, `j`, `k` parts.
  - `scalar` and `[i j k]`: the scalar part and the vector part (any sequence of three numbers, for example `Vec3`).
  - `a`: the scalar part only, the vector part is zero.

  See also [[complex->quaternion]], [[ensure-quaternion]]."
  (^Vec4 [^double a ^double b ^double c ^double d] (Vec4. a b c d))
  (^Vec4 [^double scalar [^double i ^double j ^double k]] (Vec4. scalar i j k))
  (^Vec4 [^double a] (Vec4. a 0.0 0.0 0.0)))

(defn seq->quaternion
  "Creates a quaternion from the sequence of numbers."
  ^Vec4 [[^double a ^double b ^double c ^double d]]
  (Vec4. a b c d))

(defn complex->quaternion
  "Creates a quaternion from the complex number `z` (a [[Vec2]]), `a+bi` becomes `a+bi+0j+0k`."
  ^Vec4 [^Vec2 z]
  (Vec4. (.x z) (.y z) 0.0 0.0))

(defn ensure-quaternion
  "Converts `q` to a quaternion.

  A number becomes a quaternion with the scalar part only, a complex number ([[Vec2]]) is converted by [[complex->quaternion]]. Any sequence value is returned converted to a [[Vec4]]. The arithmetic functions of this namespace use it for their arguments."
  ^Vec4 [q]
  (cond
    (instance? Vec4 q) q
    (number? q) (quaternion q)
    (instance? Vec2 q) (complex->quaternion q)
    (sequential? q) (seq->quaternion q)
    :else (throw (ex-info "Can't infer quaternion from input." {:q q}))))

(defn scalar
  "Returns the scalar part of the quaternion as a double. Alias of [[re]]."
  ^double [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (.x quaternion)))

(defn re
  "Returns the scalar (real) part of the quaternion as a double. Alias of [[scalar]]."
  ^double [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (.x quaternion)))

(defn vector
  "Returns the vector part (the `i`, `j`, `k` parts) of the quaternion as [[Vec3]]."
  ^Vec3 [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (Vec3. (.y quaternion)
           (.z quaternion)
           (.w quaternion))))

(defn im-i "Returns the `i` imaginary part as a double." ^double [^Vec4 quaternion] (.y quaternion))
(defn im-j "Returns the `j` imaginary part as a double." ^double [^Vec4 quaternion] (.z quaternion))
(defn im-k "Returns the `k` imaginary part as a double." ^double [^Vec4 quaternion] (.w quaternion))

(defn real?
  "Returns `true` when the vector part of the quaternion is zero, a real number."
  [quaternion] (v/is-zero? (vector quaternion)))

(defn imaginary?
  "Returns `true` when the scalar part of the quaternion is zero, a pure imaginary quaternion (the zero quaternion included)."
  [quaternion] (m/zero? (.x ^Vec4 (ensure-quaternion quaternion))))

(defn zero?
  "Returns `true` when all parts of the quaternion are zero."
  [quaternion] (v/is-zero? quaternion))

(defn inf?
  "Returns `true` when any part of the quaternion is infinite (also when another part is `NaN`)."
  [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (or (m/inf? (.x quaternion))
        (m/inf? (.y quaternion))
        (m/inf? (.z quaternion))
        (m/inf? (.w quaternion)))))

(defn nan?
  "Returns `true` when any part of the quaternion is `NaN`."
  [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (or (m/nan? (.x quaternion))
        (m/nan? (.y quaternion))
        (m/nan? (.z quaternion))
        (m/nan? (.w quaternion)))))

(defn invalid?
  "Returns `true` when the quaternion has an infinite or `NaN` part. See also [[valid?]]."
  [z] (or (inf? z) (nan? z)))

(defn valid?
  "Returns `true` when all parts of the quaternion are finite numbers. See also [[invalid?]]."
  [z] (not (invalid? z)))

(defn delta-eq
  "Checks whether quaternions `q1` and `q2` are equal within a tolerance. Returns `true` or `false`.

  The tolerance is applied to the length of the difference, `|q1 - q2|`, not to each part separately:

  * `(delta-eq q1 q2)` - absolute tolerance `1.0e-6`,
  * `(delta-eq q1 q2 abs-tol)` - the quaternions are equal when `|q1 - q2| <= abs-tol`,
  * `(delta-eq q1 q2 abs-tol rel-tol)` - the quaternions are equal when `|q1 - q2| < max(abs-tol, rel-tol * max(|q1|, |q2|))`.

  A part that is `NaN` makes the result `false`. Note that `q` and `-q` describe the same rotation but they are not equal. Delegates to [[fastmath.vector/delta-eq]]. See also [[norm]]."
  ([q1 q2]
   (v/delta-eq q1 q2))
  ([q1 q2 ^double abs-tol]
   (v/delta-eq q1 q2 abs-tol))
  ([q1 q2 ^double abs-tol ^double rel-tol]
   (v/delta-eq q1 q2 abs-tol rel-tol)))

(defn arg
  "Argument of the quaternion, `atan2(|vector(q)|, re(q))`, in the range `[0, pi]`.

  See also [[norm]], [[to-axis-angle]]."
  ^double [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)]
    (m/atan2 (v/mag (vector quaternion)) (.x quaternion))))

(defn norm
  "Norm (magnitude, length) of the quaternion, `sqrt(a^2+b^2+c^2+d^2)`, as a double.

  The value is correct over the whole range of doubles. Note that [[fastmath.complex/norm]] is the squared magnitude.

  See also [[normalize]], [[fastmath.vector/magsq]]."
  ^double [quaternion] (v/mag quaternion))

(defn normalize
  "Normalizes the quaternion to the unit length (a versor). The zero quaternion stays zero.

  See also [[norm]]."
  ^Vec4 [quaternion] (v/normalize quaternion))

(defn add
  "Sum of quaternions.

  Arguments can be quaternions, complex numbers or numbers. Without arguments returns `0`, with one argument the argument as a quaternion. More than two arguments are added from left to right.

  See also [[adds]], [[sub]], [[mult]]."
  (^Vec4 [] ZERO)
  (^Vec4 [q] (ensure-quaternion q))
  (^Vec4 [q1 q2]
   (let [^Vec4 q1 (ensure-quaternion q1)
         ^Vec4 q2 (ensure-quaternion q2)]
     (Vec4. (m/+ (.x q1) (.x q2))
            (m/+ (.y q1) (.y q2))
            (m/+ (.z q1) (.z q2))
            (m/+ (.w q1) (.w q2)))))
  (^Vec4 [q1 q2 & more]
   (reduce add (add q1 q2) more)))

(defn adds
  "Adds the real number `s` to the scalar part of the quaternion `q1`."
  ^Vec4 [q1 ^double s] (add q1 (Vec4. s 0.0 0.0 0.0)))

(defn neg
  "Negation of the quaternion, `-q`."
  ^Vec4 [quaternion] (v/sub (ensure-quaternion quaternion)))

(defn sub
  "Difference of quaternions, `q1 - q2 - ...`.

  Arguments can be quaternions, complex numbers or numbers. With one argument returns the negation [[neg]]. More than two arguments are subtracted from the first one, from left to right.

  See also [[add]], [[neg]]."
  (^Vec4 [q] (neg q))
  (^Vec4 [q1 q2]
   (let [^Vec4 q1 (ensure-quaternion q1)
         ^Vec4 q2 (ensure-quaternion q2)]
     (Vec4. (m/- (.x q1) (.x q2))
            (m/- (.y q1) (.y q2))
            (m/- (.z q1) (.z q2))
            (m/- (.w q1) (.w q2)))))
  (^Vec4 [q1 q2 & more]
   (reduce sub (sub q1 q2) more)))

(defn scale
  "Multiplies the quaternion by the real number `scale`."
  ^Vec4 [quaternion ^double scale] (v/mult (ensure-quaternion quaternion) scale))

(defn conjugate
  "Conjugate of the quaternion, `a-bi-cj-dk`."
  ^Vec4 [^Vec4 quaternion]
  (Vec4. (.x quaternion)
         (- (.y quaternion))
         (- (.z quaternion))
         (- (.w quaternion))))

(defn reciprocal
  "Reciprocal (multiplicative inverse) of the quaternion, `conjugate(q) / |q|^2`.

  Computed without overflow or underflow in the squared norm. The reciprocal of the zero quaternion has `NaN` parts, the reciprocal of an infinite quaternion is zero.

  See also [[inverse]], [[div]]."
  ^Vec4 [quaternion]
  (let [^Vec4 q (ensure-quaternion quaternion)
        mg (v/mag q)]
    (cond
      (m/inf? mg) (Vec4. (m/copy-sign 0.0 (.x q))
                         (m/- (m/copy-sign 0.0 (.y q)))
                         (m/- (m/copy-sign 0.0 (.z q)))
                         (m/- (m/copy-sign 0.0 (.w q))))
      (m/< 0.0 mg 1.0e-300)
      ;; 1/|q| overflows for a subnormal length: scale by a power of two (exact), invert, scale back
      (v/mult ^Vec4 (reciprocal (v/mult q 4.149515568880993E180)) 4.149515568880993E180)
      ;; numerical stability
      :else (let [r (m// 1.0 mg)]
              (Vec4. (m/* r (.x q) r)
                     (m/* r (m/- (.y q)) r)
                     (m/* r (m/- (.z q)) r)
                     (m/* r (m/- (.w q)) r))))))

(defn inverse
  "Multiplicative inverse of the quaternion, alias of [[reciprocal]]. For a unit quaternion (a rotation) it is the conjugate, the opposite rotation."
  ^Vec4 [quaternion] (reciprocal quaternion))

(defn mult
  "Product of quaternions (the Hamilton product).

  The multiplication is not commutative, `i*j = k` and `j*i = -k`. Arguments can be quaternions, complex numbers or numbers. Without arguments returns `1`, with one argument the argument as a quaternion. More than two arguments are multiplied from left to right: `(mult a b c)` is `(a*b)*c`.

  See also [[div]], [[sq]], [[scale]]."
  (^Vec4 [] ONE)
  (^Vec4 [q] (ensure-quaternion q))
  (^Vec4 [q1 q2]
   (let [^Vec4 q1 (ensure-quaternion q1)
         ^Vec4 q2 (ensure-quaternion q2)
         a1 (.x q1) b1 (.y q1) c1 (.z q1) d1 (.w q1)
         a2 (.x q2) b2 (.y q2) c2 (.z q2) d2 (.w q2)]
     (Vec4. (m/- (m/* a1 a2) (m/* b1 b2) (m/* c1 c2) (m/* d1 d2))
            (m/+ (m/* a1 b2) (m/* b1 a2) (m/* c1 d2) (m/- (m/* d1 c2)))
            (m/+ (m/* a1 c2) (m/- (m/* b1 d2)) (m/* c1 a2) (m/* d1 b2))
            (m/+ (m/* a1 d2) (m/* b1 c2) (m/- (m/* c1 b2)) (m/* d1 a2)))))
  (^Vec4 [q1 q2 & more]
   (reduce mult (mult q1 q2) more)))

(defn div
  "Division of quaternions, `q1 * q2^-1` (the divisor is applied from the right).

  Arguments can be quaternions, complex numbers or numbers. See also [[mult]], [[reciprocal]]."
  ^Vec4 [q1 q2]
  (mult q1 (reciprocal q2)))

(defn dot
  "Dot product of quaternions treated as 4-dimensional vectors, `a1*a2 + b1*b2 + c1*c2 + d1*d2`.

  For unit quaternions it is the cosine of half of the angle between the rotations, the sign shows the shortest arc (see [[slerp]]). Arguments can be quaternions, complex numbers or numbers."
  ^double [q1 q2]
  (v/dot (ensure-quaternion q1) (ensure-quaternion q2)))

(defn sq
  "Square of the quaternion, `q*q`."
  ^Vec4 [quaternion]
  (let [^Vec4 quaternion (ensure-quaternion quaternion)
        a (.x quaternion)
        b (.y quaternion)
        c (.z quaternion)
        d (.w quaternion)
        aa (* 2.0 a)]
    (Vec4. (m/- (* a a) (* b b) (* c c) (* d d))
           (* aa b)
           (* aa c)
           (* aa d))))

(defn qsgn
  "Computes the signum of a quaternion.

  Returns `0.0` for the zero quaternion (`0+0i+0j+0k`). For any other quaternion, returns the sign of its scalar part (`-1.0` or `1.0`), not a unit quaternion (see [[normalize]])."
  (^double [^double re ^double im-i ^double im-j ^double im-k]
   (if (and (m/zero? re)
            (m/zero? im-i)
            (m/zero? im-j)
            (m/zero? im-k)) 0.0 (m/sgn re)))
  (^double [q]
   (let [^Vec4 q (ensure-quaternion q)]
     (if (zero? q) 0.0 (m/sgn (.x q))))))

;; transcendental functions: f(a + b i + c j + d k) = f(a + |v| i) with the imaginary part of the result
;; spread over the direction of the vector part v

(defn- vector-part
  "Component of the vector part for the direction `b/absim` and the imaginary part `imz` of the complex result."
  ^double [^double b ^double absim ^double imz]
  (let [u (/ b absim)]
    (if (m/zero? u) u (* u imz))))

(defmacro ^:private gen-from-complex
  [sym desc]
  (let [nm (str sym)
        src (symbol "fastmath.complex" nm)
        q (with-meta (symbol "q") {:tag 'fastmath.vector.Vec4})
        z (with-meta (symbol "z") {:tag 'fastmath.vector.Vec2})]
    `(defn ~sym
       ~(str desc " of the quaternion `q`.\n\n"
             "  The value is computed with [[fastmath.complex/" nm "]] applied to `re(q) + |vector(q)| i`, the imaginary part of the result is spread over the direction of the vector part of `q`. "
             "A quaternion with a zero vector part is treated as a complex number with the imaginary part on the `i` axis. "
             "Branch cuts and the accuracy are the ones of the complex function.\n\n"
             "  See also [[fastmath.complex/" nm "]].")
       ^Vec4 [~'quaternion]
       (let [~q (ensure-quaternion ~'quaternion)
             a# (.x ~q)
             b# (.y ~q)
             c# (.z ~q)
             d# (.w ~q)
             absim# (m/hypot b# c# d#)
             ~z (~src (c/complex a# absim#))]
         (if (m/pos? absim#)
           (let [imz# (.y ~z)]
             (Vec4. (.x ~z)
                    (vector-part b# absim# imz#)
                    (vector-part c# absim# imz#)
                    (vector-part d# absim# imz#)))
           (Vec4. (.x ~z) (.y ~z) 0.0 0.0))))))

(gen-from-complex sqrt "Principal square root")
(gen-from-complex exp "Exponential")
(gen-from-complex log "Principal natural logarithm")

(defn logb
  "Logarithm of the quaternion `q` with the base `b`, `log(q) / log(b)`.

  Both arguments can be quaternions, complex numbers or numbers. See also [[log]]."
  ^Vec4 [quaternion b]
  (div (log (ensure-quaternion quaternion)) (log (ensure-quaternion b))))

(defn pow
  "Power of a quaternion, `exp(log(q) * p)`.

  Both arguments can be quaternions, complex numbers or numbers. Note that the exponent multiplies `log(q)` from the right and that the quaternion product is not commutative. The power of the zero quaternion has `NaN` parts.

  See also [[exp]], [[log]], [[sq]]."
  ^Vec4 [q p]
  (exp (mult (log (ensure-quaternion q)) p)))

(defn- check-direction
  "Throws `ex-info` when the vector `u` can't give a direction: the zero vector or a vector with an infinite or `NaN` part.

  A finite vector is accepted also when its length is not representable (`(1.7e308, 1.7e308, 1.7e308)`): then the length of the half of the vector is finite."
  [option u]
  (let [mg (v/mag u)
        finite? (or (m/zero? mg)
                    (not (or (m/inf? mg) (m/nan? mg)))
                    (let [mh (v/mag (v/mult u 0.5))]
                      (not (or (m/inf? mh) (m/nan? mh)))))]
    (when (or (m/zero? mg) (not finite?))
      (throw (ex-info "The vector must be finite and non zero" {:option option :value u :expected "finite non zero vector"})))))

(defn rotation-quaternion
  "Creates a unit quaternion representing a rotation.

  The rotation is defined by an `angle` (in radians) and a 3D vector `u`
  specifying the axis of rotation. The axis vector `u` is normalized
  internally to ensure a unit quaternion result.

  Parameters:

  * `angle`: The rotation angle in radians (double).
  * `u`: The axis of rotation (Vec3). It will be normalized before use.

  Returns A unit quaternion (Vec4) representing the rotation. Throws `ex-info` when `u` is the zero vector or has an infinite or `NaN` part (a finite vector is fine also when its length is not representable).

  See also [[from-axis-angle]], [[to-axis-angle]], [[rotate]]."
  ^Vec4 [^double angle u]
  (check-direction :u u)
  (let [half (m/* 0.5 angle)]
    (quaternion (m/cos half)
                (v/mult (v/normalize u) (m/sin half)))))

(defn from-axis-angle
  "Creates a unit quaternion representing a rotation by `angle` (radians) around the `axis`.

  Alias of [[rotation-quaternion]], the arguments have the same order. The axis is normalized, a zero axis throws `ex-info`. See also [[to-axis-angle]]."
  ^Vec4 [^double angle axis]
  (rotation-quaternion angle axis))

(defn to-axis-angle
  "Converts a quaternion to the axis and the angle of the rotation.

  The input quaternion is normalized. Returns a map with the keys:

  * `:axis`: the unit axis of rotation (`Vec3`),
  * `:angle`: the angle in radians in the range `[0, pi]`.

  The sign of the quaternion is chosen so that the angle is the shortest one (for `-q` the result is the same as for `q`). For a quaternion with a zero vector part (no rotation) the axis is `[1 0 0]` and the angle is `0`.

  See also [[from-axis-angle]], [[arg]]."
  [quaternion]
  (let [^Vec4 q (normalize quaternion)
        flip (if (m/neg? (.x q)) -1.0 1.0)
        a (m/* flip (.x q))
        vv (v/mult (vector q) flip)
        s (v/mag vv)]
    (if (m/zero? s)
      {:axis (Vec3. 1.0 0.0 0.0) :angle 0.0}
      {:axis (v/div vv s) :angle (m/* 2.0 (m/atan2 s a))})))

(defn from-two-vectors
  "Creates a unit quaternion representing the shortest rotation which takes the direction of the vector `u` to the direction of the vector `v`.

  The lengths of the vectors do not matter (also subnormal ones). For opposite directions the rotation is by `pi` around some axis perpendicular to `u`. Throws `ex-info` when any vector is the zero vector or has an infinite or `NaN` part (a finite vector is fine also when its length is not representable).

  See also [[rotate]], [[rotation-quaternion]]."
  ^Vec4 [u v]
  (check-direction :u u)
  (check-direction :v v)
  (let [^Vec3 nu (v/vec3 (v/normalize u))
        ^Vec3 nv (v/vec3 (v/normalize v))
        d (v/dot nu nv)]
    (if (m/< d -0.999999999999999)
      ;; opposite directions: a half turn around any axis perpendicular to u
      (let [^Vec3 ax (v/normalize (if (m/> (m/abs (.x nu)) 0.9)
                                    (v/cross nu (Vec3. 0.0 1.0 0.0))
                                    (v/cross nu (Vec3. 1.0 0.0 0.0))))]
        (quaternion 0.0 ax))
      (normalize (quaternion (m/+ 1.0 d) (v/cross nu nv))))))

(defn rotate
  "Rotates the 3d vector `in` by the rotation quaternion `rotq` or around the axis `u` by the `angle` (radians), the same as `fastmath.vector/axis-rotate`.

  The rotation quaternion is normalized before use, so a non-unit quaternion represents the same rotation as its unit version (the vector is not scaled). The zero quaternion does not rotate.

  Parameters:

  * `in`: the vector to rotate, any sequence of three numbers.
  * `rotq`: the rotation quaternion (see [[rotation-quaternion]], [[from-euler]]).
  * `angle`, `u`: the angle (radians) and the axis of rotation (`Vec3`).

  Returns the rotated vector as `Vec3`."
  (^Vec3 [in rotq]
   (let [in (apply v/vec3 in)
         ^Vec4 rotq (v/vec4 (v/normalize rotq))
         qw (.x rotq)
         q (vector rotq)
         t (v/mult (v/cross q in) 2.0)]
     (v/add (v/add in (v/mult t qw)) (v/cross q t))))
  (^Vec3 [in ^double angle u]
   (rotate in (rotation-quaternion angle u))))

(defn nlerp
  "Normalized linear interpolation between the unit quaternions `q1` and `q2`.

  It follows the shortest arc (when `dot(q1, q2)` is negative `q2` is negated), the result is normalized. It is faster than [[slerp]] and the same on the path, the angular velocity is not constant. `t` is constrained to `[0, 1]`, `NaN` gives `NaN`.

  See also [[slerp]], [[dot]]."
  ^Vec4 [q1 q2 ^double t]
  (if (m/nan? t)
    (Vec4. ##NaN ##NaN ##NaN ##NaN)
    (let [^Vec4 q1 (ensure-quaternion q1)
          ^Vec4 q2 (ensure-quaternion q2)
          t (m/constrain t 0.0 1.0)
          s (if (m/neg? (v/dot q1 q2)) -1.0 1.0)]
      (v/normalize (v/add (v/mult q1 (m/- 1.0 t))
                          (v/mult q2 (m/* s t)))))))

(defn slerp
  "Performs Spherical Linear Interpolation (SLERP) between two unit quaternions.

  SLERP interpolates along the shortest arc on the unit sphere, providing smooth
  interpolation between rotations represented by unit quaternions with a constant angular velocity. When `dot(q1, q2)` is negative `q2` is negated (the same rotation) to take the shorter way. For almost parallel quaternions [[nlerp]] is used, see below.

  Parameters:

  * `q1`: The starting quaternion (Vec4).
  * `q2`: The ending quaternion (Vec4).
  * `t`: The interpolation parameter (double). Should be in the range [0.0, 1.0].
       - If `t=0.0`, returns `q1`.
       - If `t=1.0`, returns `q2`.
       - For `0.0 < t < 1.0`, returns an interpolated unit quaternion.
       The parameter `t` is internally constrained to [0.0, 1.0]; `NaN` gives a quaternion with `NaN` parts.

  The quaternions are expected to have the unit length, they are not normalized (the arguments are returned unchanged for `t=0` and `t=1`). Both can be numbers or complex numbers too.

  For almost parallel quaternions (`dot(q1, q2) > 0.9995`, an angle below about 3.6 degrees) the normalized linear interpolation is used: the angular velocity is then not exactly constant, the result differs from the exact slerp by less than 1e-6 (4e-7 measured at `dot = 0.99955`).

  See also [[nlerp]], [[dot]]."
  ^Vec4 [q1 q2 ^double t]
  (let [^Vec4 q1 (ensure-quaternion q1)
        ^Vec4 q2 (ensure-quaternion q2)
        t (m/constrain t 0.0 1.0)]
    (cond
      (m/nan? t) (Vec4. ##NaN ##NaN ##NaN ##NaN)
      (m/zero? t) q1
      (m/one? t) q2
      (and (zero? q1) (zero? q2)) ZERO
      :else (let [d (v/dot q1 q2)
                  s (if (m/neg? d) -1.0 1.0)
                  d (m/abs d)]
              (if (m/> d 0.9995)
                (nlerp q1 q2 t)
                (let [theta (m/acos (m/min d 1.0))
                      sin-theta (m/sin theta)
                      w1 (m// (m/sin (m/* (m/- 1.0 t) theta)) sin-theta)
                      w2 (m// (m/* s (m/sin (m/* t theta))) sin-theta)]
                  (v/add (v/mult q1 w1) (v/mult q2 w2))))))))

;; the angle at the gimbal lock is twice an arctangent of the (half angle) quaternion parts: wrapped to [-pi, pi]
(defn- wrap-angle
  ^double [^double a]
  (cond
    (m/> a m/PI) (m/- a m/TWO_PI)
    (m/< a (m/- m/PI)) (m/+ a m/TWO_PI)
    :else a))

;; https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Euler_angles_(in_3-2-1_sequence)_to_quaternion_conversion
(defn to-euler
  "Converts a quaternion `q` to ZYX (body 3-2-1) Euler angles.
  The input quaternion is normalized before calculation.
  Returns a `Vec3` representing the angles `[roll pitch yaw]`.
  Roll is the angle around the x-axis, pitch around the y'-axis, and yaw around the z''-axis.
  Output angles are in [-pi, pi] for roll/yaw and [-pi/2, pi/2] for pitch.

  At the gimbal lock (pitch is `+-pi/2`, only the sum or the difference of roll and yaw is defined) the roll is set to zero and the yaw takes the whole rotation, so [[from-euler]] of the result represents the same rotation. The lock is detected when `cos(pitch)` is below `1.0e-8`.

  See also [[from-euler]], [[to-angles]]."
  ^Vec3 [q]
  (let [^Vec4 q (v/normalize (ensure-quaternion q))
        w (.x q) x (.y q) y (.z q) z (.w q)
        sinr-cosp (m/* 2.0 (m/+ (m/* w x) (m/* y z)))
        cosr-cosp (m/- 1.0 (m/* 2.0 (m/+ (m/* x x) (m/* y y))))
        p (m/* 2.0 (m/- (m/* w y) (m/* x z)))
        siny-cosp (m/* 2.0 (m/+ (m/* w z) (m/* x y)))
        cosy-cosp (m/- 1.0 (m/* 2.0 (m/+ (m/* y y) (m/* z z))))
        cosp (m/hypot sinr-cosp cosr-cosp)
        pitch (m/atan2 p cosp)]
    (if (m/< cosp 1.0e-8)
      (Vec3. 0.0 pitch (wrap-angle (m/* 2.0 (if (m/pos? p) (m/atan2 (m/- x) w) (m/atan2 x w)))))
      (Vec3. (m/atan2 sinr-cosp cosr-cosp)
             pitch
             (m/atan2 siny-cosp cosy-cosp)))))

(defn from-euler
  "Converts Euler angles (ZYX body 3-2-1 convention) to a quaternion.

  The rotation sequence is intrinsic Z-Y'-X'' (yaw, pitch, roll), applied to the body frame.
  The order of input parameters or vector components is [roll pitch yaw].

  Parameters:

  * `roll`: Rotation around the x-axis (radians), expected range `[-pi, pi]`.
  * `pitch`: Rotation around the y'-axis (radians), expected range `[-pi/2, pi/2]`.
  * `yaw`: Rotation around the z''-axis (radians), expected range `[-pi, pi]`.

  Can accept individual double values or a `Vec3` containing [roll pitch yaw].

  Returns A quaternion `Vec4` representing the rotation (a unit quaternion).

  See also [[to-euler]], [[from-angles]]."
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
  "Converts a quaternion `q` to Tait–Bryan angles using the z-y′-x'' intrinsic rotation sequence.
  The input quaternion is normalized before calculation.
  Returns a `Vec3` representing the angles `[x y z]`. These angles correspond to rotations around the local (body) x-axis, followed by the intermediate local y'-axis, and finally the local z''-axis.
  Output angles are in the range `[-pi, pi]` for x and z, and `[-pi/2, pi/2]` for y.

  At the gimbal lock (y is `+-pi/2`, only the sum or the difference of x and z is defined) x is set to zero and z takes the whole rotation, so [[from-angles]] of the result represents the same rotation. The lock is detected when `cos(y)` is below `1.0e-8`.

  See also [[from-angles]], [[to-euler]]."
  ^Vec3 [q]
  (let [^Vec4 q (v/normalize (ensure-quaternion q))
        w (.x q) x (.y q) y (.z q) z (.w q)
        ww (m/* w w) xx (m/* x x) yy (m/* y y) zz (m/* z z)
        y1 (m/* 2.0 (m/- (m/* w x) (m/* y z)))
        x1 (m/+ (m/- ww xx yy) zz)
        a2 (m/* 2.0 (m/+ (m/* w y) (m/* x z)))
        y3 (m/* 2.0 (m/- (m/* w z) (m/* x y)))
        x3 (m/+ ww (m/- xx yy zz))
        cosy (m/hypot y1 x1)
        ay (m/atan2 a2 cosy)]
    (if (m/< cosy 1.0e-8)
      (Vec3. 0.0 ay (wrap-angle (m/* 2.0 (if (m/pos? a2)
                                           (m/atan2 (m/+ x z) (m/+ y w))
                                           (m/- (m/atan2 (m/- x z) (m/- w y)))))))
      (Vec3. (m/atan2 y1 x1)
             ay
             (m/atan2 y3 x3)))))

(defn from-angles
  "Converts Tait–Bryan angles (z-y′-x'' intrinsic rotation sequence) to a quaternion.

  The angles `[x y z]` correspond to rotations around the local (body) x-axis, followed by the intermediate local y'-axis, and finally the local z''-axis.

  Parameters:

  * `x`: Rotation around the x-axis (radians).
  * `y`: Rotation around the y'-axis (radians).
  * `z`: Rotation around the z''-axis (radians).

  Can accept individual double values or a `Vec3` containing `[x y z]`.

  The rotation matrix of the result is `Rx(x) Ry(y) Rz(z)`, the same as [[fastmath.matrix/rotation-matrix-3d]]. Returns A quaternion `Vec4` representing the rotation (a unit quaternion).

  See also [[to-angles]], [[from-euler]]."
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
  "Converts a quaternion to a 3x3 rotation matrix.

  The input quaternion is normalized internally to ensure a valid rotation matrix.
  Unit quaternions are used to represent rotations. Both `q` and `-q` give the same matrix.

  Parameters:

  * `q`: The quaternion (Vec4) to convert.

  Returns A 3x3 rotation matrix (Mat3x3).

  See also [[from-rotation-matrix]]."
  ^Mat3x3 [q]
  (let [^Vec4 q (v/normalize (ensure-quaternion q))
        a (.x q) b (.y q) c (.z q) d (.w q)
        aa (m/* a a) bb (m/* b b) cc (m/* c c) dd (m/* d d)
        bc (m/* b c) ad (m/* a d) bd (m/* b d)
        ac (m/* a c) cd (m/* c d) ab (m/* a b)]
    (Mat3x3. (m/+ aa (m/- bb cc dd)) (m/* 2.0 (m/- bc ad)) (m/* 2.0 (m/+ bd ac))
             (m/* 2.0 (m/+ bc ad)) (m/+ cc (m/- aa bb dd)) (m/* 2.0 (m/- cd ab))
             (m/* 2.0 (m/- bd ac)) (m/* 2.0 (m/+ cd ab)) (m/+ dd (m/- aa bb cc)))))

;; https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
(defn from-rotation-matrix
  "Converts a 3x3 rotation matrix to a quaternion.

  Takes a [[Mat3x3]] rotation matrix as input.
  Returns a [[Vec4]] representing the quaternion that represents the same rotation.

  The resulting quaternion is a unit quaternion if the input matrix is a valid rotation matrix. The method handles numerical stability and normalization. A rotation corresponds to two quaternions, `q` and `-q`; the sign of the result is not defined (the sign of the scalar part is not fixed).

  See also [[to-rotation-matrix]]."
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

(gen-from-complex sin "Sine")
(gen-from-complex cos "Cosine")
(gen-from-complex tan "Tangent")
(gen-from-complex sec "Secant")
(gen-from-complex csc "Cosecant")
(gen-from-complex cot "Cotangent")
(gen-from-complex sinh "Hyperbolic sine")
(gen-from-complex cosh "Hyperbolic cosine")
(gen-from-complex tanh "Hyperbolic tangent")
(gen-from-complex sech "Hyperbolic secant")
(gen-from-complex csch "Hyperbolic cosecant")
(gen-from-complex coth "Hyperbolic cotangent")
(gen-from-complex asin "Inverse sine")
(gen-from-complex acos "Inverse cosine")
(gen-from-complex atan "Inverse tangent")
(gen-from-complex asec "Inverse secant")
(gen-from-complex acsc "Inverse cosecant")
(gen-from-complex acot "Inverse cotangent")
(gen-from-complex asinh "Inverse hyperbolic sine")
(gen-from-complex acosh "Inverse hyperbolic cosine")
(gen-from-complex atanh "Inverse hyperbolic tangent")
(gen-from-complex asech "Inverse hyperbolic secant")
(gen-from-complex acsch "Inverse hyperbolic cosecant")
(gen-from-complex acoth "Inverse hyperbolic cotangent")

;; https://math.stackexchange.com/questions/1499095/how-to-calculate-sin-cos-tan-of-a-quaternion
