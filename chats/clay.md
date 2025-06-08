### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace. 

### Please write a short paragraph for `clay/complex_quaternions.clj` book. Do not copy content of the `clay/complex_quaternions.clj` file, write only asked chapter. Put text as a Clojure comment block. Put actual Clojure code without comment block and markdown wrapper.  Follow the convention of other chapters.

- Write a paragraph for the "## Quaternions" chapter.
- Describe all mentioned functions, group them in points and describe them in brief
- When possible add latex formulas
- Illustrate with examples
- Rely on the structure and content of `Complex numbers` chapter

Write only asked chapter. Refer to previous chapters as examples of chapter structure.



```clojure
;; ## Quaternions
;;
;; Quaternions extend complex numbers with three imaginary units: $i$, $j$, and $k$, satisfying the relations $i^2 = j^2 = k^2 = ijk = -1$. They are particularly useful for representing 3D rotations concisely and avoiding gimbal lock.
;;
;; In `fastmath`, quaternions are represented as 4-dimensional vectors using the `fastmath.vector/Vec4` type. A quaternion $q = a + bi + cj + dk$ corresponds to the vector `(Vec4. a b c d)`, where $a$ is the scalar (real) part, and $(b, c, d)$ is the vector (imaginary) part.
;;
;; You can create quaternions using the dedicated `quaternion` function, or by converting numbers or complex numbers using `ensure-quaternion`.
;;
;; The implementation handles floating-point subtleties, including special values (`##Inf`, `##NaN`), for robust calculations across a wide range of inputs.
;;
;; This section covers constants, creation, accessors, predicates, basic operations, and advanced functions including power, logarithms, trigonometric functions, and functions specific to 3D rotations.
;;

;; ::: {.callout-tip title="Defined functions"}
;; * Constants: `ZERO`, `ONE`, `I`, `J`, `K`, `-I`, `-J`, `-K`
;; * Creation and Accessors: `quaternion`, `complex->quaternion`, `ensure-quaternion`, `scalar`, `re`, `vector`, `im-i`, `im-j`, `im-k`
;; * Predicates and Properties: `real?`, `imaginary?`, `zero?`, `inf?`, `nan?`, `invalid?`, `valid?`, `delta-eq`, `qsgn`, `arg`, `norm`, `normalize`, `conjugate`, `reciprocal`
;; * Arithmetic Operations: `add`, `adds`, `sub`, `scale`, `mult`, `div`, `neg`, `sq`
;; * Power and Logarithms: `pow`, `exp`, `log`, `logb`
;; * Trigonometric and Hyperbolic Functions: `sin`, `cos`, `tan`, `sec`, `csc`, `cot`, `sinh`, `cosh`, `tanh`, `sech`, `csch`, `coth`
;; * Inverse Trigonometric and Hyperbolic Functions: `asin`, `acos`, `atan`, `asec`, `acsc`, `acot`, `asinh`, `acosh`, `atanh`, `asech`, `acsch`, `acoth`
;; * Rotation Functions: `rotation-quaternion`, `rotate`, `slerp`, `to-euler`, `from-euler`, `to-angles`, `from-angles`, `to-rotation-matrix`, `from-rotation-matrix`
;; :::
;;
;; We will use these quaternions to illustrate operations:

(def q1 (quat/quaternion 0.5 -1.5 2.0 -0.5))
(def q2 (quat/quaternion -0.5 2.5 -1.0 1.5))
(def q3 (quat/quaternion 1.5 0.5 0.0 0.0)) ;; A real quaternion
(def q4 (quat/quaternion 0.0 1.0 -2.0 3.0)) ;; A pure imaginary quaternion

(utls/examples-note q1 q2 q3 q4)

;; ### Constants
;;
;; Predefined quaternion constants for common values.

(kind/table
 {:column-names ["constant" "q" "value"]
  :row-vectors [['ZERO (kind/md "$0+0i+0j+0k$") (str quat/ZERO)]
                ['ONE (kind/md "$1+0i+0j+0k$") (str quat/ONE)]
                ['I (kind/md "$0+i+0j+0k$") (str quat/I)]
                ['J (kind/md "$0+0i+j+0k$") (str quat/J)]
                ['K (kind/md "$0+0i+0j+k$") (str quat/K)]
                ['-I (kind/md "$0-i+0j+0k$") (str quat/-I)]
                ['-J (kind/md "$0+0i-j+0k$") (str quat/-J)]
                ['-K (kind/md "$0+0i+0j-k$") (str quat/-K)]]})

;; ### Creation and Accessors
;;
;; Functions to create quaternions and access their scalar and vector components.
;;
;; *   `quaternion`: Creates a quaternion from individual scalar and imaginary parts or from a scalar and a 3D vector.
;; *   `complex->quaternion`: Converts a complex number (real part `a`, imaginary part `b`) into a quaternion $a+bi+0j+0k$.
;; *   `ensure-quaternion`: Converts a number or a complex number into a quaternion.
;; *   `scalar`, `re`: Returns the scalar (real) part $a$ of a quaternion $a+bi+cj+dk$.
;; *   `vector`: Returns the vector (imaginary) part $(b, c, d)$ as a `Vec3`.
;; *   `im-i`, `im-j`, `im-k`: Return the coefficients of the $i$, $j$, and $k$ components, respectively.

(utls/examples-note
  (quat/quaternion 1.0 -2.0 3.0 -4.0)
  (quat/quaternion 5.0 [1.0 2.0 3.0])
  (quat/quaternion 7.0)
  (quat/complex->quaternion (cplx/complex 1.0 -2.0))
  (quat/ensure-quaternion 5.0)
  (quat/ensure-quaternion (cplx/complex 1.0 -2.0))
  (quat/scalar q1)
  (quat/re q1)
  (quat/vector q1)
  (quat/im-i q1)
  (quat/im-j q1)
  (quat/im-k q1))

;; ### Predicates and Properties
;;
;; Functions to check the type or state of a quaternion and compute fundamental properties.
;;
;; *   `real?`: Checks if the quaternion has a zero vector part ($b=c=d=0$).
;; *   `imaginary?`: Checks if the quaternion has a zero scalar part ($a=0$).
;; *   `zero?`: Checks if all components are zero ($a=b=c=d=0$).
;; *   `inf?`, `nan?`: Checks if any component is infinite or NaN, respectively.
;; *   `invalid?`, `valid?`: Checks if the quaternion is not a finite, non-NaN value.
;; *   `delta-eq`: Checks if two quaternions are approximately equal within a tolerance.
;; *   `qsgn`: Computes the complex signum, returning the sign of the scalar part or 0 for zero.
;; *   `arg`: Computes the argument $\theta$ such that $q = |q| (\cos \theta + \mathbf{u} \sin \theta)$, where $\mathbf{u}$ is a unit 3D vector. It's calculated as $\operatorname{atan2}(|\operatorname{vector}(q)|, \operatorname{scalar}(q))$.
;; *   `norm`: Calculates the magnitude (Euclidean norm) of the quaternion $|q| = \sqrt{a^2+b^2+c^2+d^2}$.
;; *   `normalize`: Returns the unit quaternion $\hat{q} = q / |q|$.
;; *   `conjugate`: Returns the conjugate $\bar{q} = a - bi - cj - dk$.
;; *   `reciprocal`: Returns the reciprocal $q^{-1} = \bar{q} / |q|^2$.

(utls/examples-note
  (quat/real? q3)
  (quat/real? q1)
  (quat/imaginary? q4)
  (quat/imaginary? q1)
  (quat/zero? quat/ZERO)
  (quat/zero? q1)
  (quat/inf? (quat/quaternion ##Inf 1 2 3))
  (quat/nan? (quat/quaternion 1 ##NaN 2 3))
  (quat/invalid? (quat/quaternion 1 2 3 4))
  (quat/valid? (quat/quaternion 1 2 3 4))
  (quat/delta-eq (quat/quaternion 1.0000001) quat/ONE)
  (quat/delta-eq (quat/quaternion 1.0001) quat/ONE)
  (quat/arg q1)
  (quat/norm q1)
  (quat/normalize q1)
  (quat/norm (quat/normalize q1))
  (quat/conjugate q1)
  (quat/reciprocal q1)
  (quat/mult q1 (quat/reciprocal q1)) ;; Should be close to ONE
  (quat/qsgn q1)
  (quat/qsgn q2)
  (quat/qsgn quat/ZERO))

;; ### Arithmetic Operations
;;
;; Basic arithmetic extended to quaternions.
;;
;; *   `add`: Sum of two quaternions, $q_1 + q_2 = (a_1+a_2) + (b_1+b_2)i + (c_1+c_2)j + (d_1+d_2)k$.
;; *   `adds`: Adds a real scalar $s$ to a quaternion $q = a+bi+cj+dk$: $q+s = (a+s)+bi+cj+dk$.
;; *   `sub`: Difference of two quaternions, $q_1 - q_2 = (a_1-a_2) + (b_1-b_2)i + (c_1-c_2)j + (d_1-d_2)k$.
;; *   `scale`: Multiplies a quaternion $q$ by a real scalar $s$: $q \cdot s = (as) + (bs)i + (cs)j + (ds)k$.
;; *   `mult`: Quaternion multiplication. Note that quaternion multiplication is non-commutative ($q_1 q_2 \neq q_2 q_1$ in general).
;;     If $q_1 = a_1 + \mathbf{v}_1$ and $q_2 = a_2 + \mathbf{v}_2$, where $\mathbf{v}_1, \mathbf{v}_2$ are vector parts:
;;     $$ q_1 q_2 = (a_1 a_2 - \mathbf{v}_1 \cdot \mathbf{v}_2) + (a_1 \mathbf{v}_2 + a_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2) $$
;; *   `div`: Quaternion division, $q_1 / q_2 = q_1 \cdot q_2^{-1}$.
;; *   `neg`: Negation, $-q = -a - bi - cj - dk$.
;; *   `sq`: Square of a quaternion, $q^2 = q \cdot q$.

(utls/examples-note
  (quat/add q1 q2)
  (quat/adds q1 5.0)
  (quat/sub q1 q2)
  (quat/scale q1 3.0)
  (quat/mult q1 q2)
  (quat/mult q2 q1) ;; Note non-commutativity
  (quat/div q1 q2)
  (quat/neg q1)
  (quat/sq q1))

;; ### Power and Logarithms
;;
;; Extensions of exponential, logarithm, and power functions to quaternions.
;;
;; *   `exp`: Computes the quaternion exponential $e^q$. If $q=a+\mathbf{v}$: $$e^q = e^{a+\mathbf{v}} = e^a (\cos |\mathbf{v}| + \mathbf{\hat{v}} \sin |\mathbf{v}|)$$ where $\mathbf{\hat{v}} = \mathbf{v}/|\mathbf{v}|$ is the unit vector part.
;; *   `log`: Computes the principal value of the quaternion natural logarithm $\log q$. If $q=a+\mathbf{v}$: $$\log q = \ln|q| + \mathbf{\hat{v}} \operatorname{arg}(q)$$ where $\operatorname{arg}(q) = \operatorname{atan2}(|\mathbf{v}|, a)$.
;; *   `logb`: Computes the logarithm of $q$ with a quaternion base $b$: $\log_b q = (\log q) (\log b)^{-1}$.
;; *   `pow`: Computes the quaternion power $q_1^{q_2}$ defined as $e^{(\log q_1) q_2}$.

(utls/examples-note
  (quat/exp q1)
  (quat/log q1)
  (quat/logb q1 q2)
  (quat/pow q1 q2)
  (quat/pow q3 q1))

;; ### Trigonometric and Hyperbolic Functions
;;
;; Extensions of standard and hyperbolic trigonometric functions to quaternions.
;; If $q=a+\mathbf{v}$:
;;
;; *   `sin`: $\sin q = \sin a \cosh |\mathbf{v}| + \mathbf{\hat{v}} \cos a \sinh |\mathbf{v}|$.
;; *   `cos`: $\cos q = \cos a \cosh |\mathbf{v}| - \mathbf{\hat{v}} \sin a \sinh |\mathbf{v}|$.
;; *   `tan`: $\tan q = (\sin q) (\cos q)^{-1}$.
;; *   `sec`, `csc`, `cot`: Reciprocals of `cos`, `sin`, `tan`.
;; *   `sinh`: $\sinh q = \sinh a \cos |\mathbf{v}| + \mathbf{\hat{v}} \cosh a \sin |\mathbf{v}|$.
;; *   `cosh`: $\cosh q = \cosh a \cos |\mathbf{v}| + \mathbf{\hat{v}} \sinh a \sin |\mathbf{v}|$.
;; *   `tanh`: $\tanh q = (\sinh q) (\cosh q)^{-1}$.
;; *   `sech`, `csch`, `coth`: Reciprocals of `cosh`, `sinh`, `tanh`.

(utls/examples-note
  (quat/sin q1)
  (quat/cos q1)
  (quat/tan q1)
  (quat/sinh q1)
  (quat/cosh q1)
  (quat/tanh q1))

;; ### Inverse Trigonometric and Hyperbolic Functions
;;
;; Implementations of the multi-valued inverse trigonometric and hyperbolic functions, returning the principal values. These are defined based on the quaternion logarithm.
;;
;; *   `asin`, `acos`, `atan`, `asec`, `acsc`, `acot`
;; *   `asinh`, `acosh`, `atanh`, `asech`, `acsch`, `acoth`

(utls/examples-note
  (quat/asin q1)
  (quat/acos q1)
  (quat/atan q1)
  (quat/asinh q1)
  (quat/acosh q1)
  (quat/atanh q1))

;; ### Rotation Functions
;;
;; Functions specifically for using quaternions to represent and manipulate 3D rotations.
;;
;; *   `rotation-quaternion`: Creates a unit quaternion representing a rotation by a given `angle` around a `u` (Vec3) axis.
;; *   `rotate`: Rotates a 3D vector `in` (Vec3) using a quaternion `rotq`. Can also create the rotation quaternion directly from angle and axis.
;; *   `slerp`: Performs Spherical Linear Interpolation between two unit quaternions `q1` and `q2` based on a parameter `t` in [0, 1]. Useful for animating rotations.
;; *   `to-euler`, `from-euler`: Converts between a quaternion and ZYX (body 3-2-1) Euler angles [roll pitch yaw].
;; *   `to-angles`, `from-angles`: Converts between a quaternion and Tait-Bryan angles (z-y'-x'' intrinsic) [x y z].
;; *   `to-rotation-matrix`, `from-rotation-matrix`: Converts between a quaternion and a 3x3 rotation matrix.

(utls/examples-note
  (def axis (v/vec3 1 1 0))
  (def angle m/HALF_PI)
  (def rot-q (quat/rotation-quaternion angle axis))
  rot-q
  (def initial-v (v/vec3 1 0 0))
  (quat/rotate initial-v rot-q)
  (quat/rotate initial-v angle axis)
  (def q-start quat/ONE)
  (def q-end (quat/rotation-quaternion m/PI (v/vec3 0 0 1)))
  (quat/slerp q-start q-end 0.0)
  (quat/slerp q-start q-end 0.5)
  (quat/slerp q-start q-end 1.0)
  (def euler-angles (v/vec3 m/HALF_PI 0 m/HALF_PI)) ;; roll=90, pitch=0, yaw=90
  (def q-from-euler (quat/from-euler euler-angles))
  q-from-euler
  (quat/to-euler q-from-euler) ;; Should be close to euler-angles
  (def angles (v/vec3 m/HALF_PI 0 m/HALF_PI)) ;; x=90, y=0, z=90
  (def q-from-angles (quat/from-angles angles))
  q-from-angles
  (quat/to-angles q-from-angles) ;; Should be close to angles
  (def rot-mat (quat/to-rotation-matrix rot-q))
  rot-mat
  (quat/from-rotation-matrix rot-mat) ;; Should be close to rot-q
  )
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-flash-preview-04-17 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (783 12279))) -->
<!-- End: -->
