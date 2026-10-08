(ns fastmath.polynomials
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.complex :as cplx]
            [clojure.string :as str]
            [fastmath.protocols.polynomials :as prot])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]
           [java.math BigDecimal BigInteger MathContext]
           [java.text DecimalFormat]
           [clojure.lang IFn]
           [umontreal.ssj.functionfit PolInterp]
           [org.apache.commons.math3.linear MatrixUtils RealMatrix EigenDecomposition]
           [org.apache.commons.math3.special Gamma]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmacro mevalpoly
  "Evaluates a real polynomial at `x` for coefficients given explicitly in the code; the macro version of [[evalpoly]].

  The polynomial is `c0 + c1*x + c2*x^2 + ...`, with coefficients in ascending order of power.

  Parameters:

  - `x` (number): the point of evaluation. It can be evaluated several times, so pass a symbol or a literal rather than an expensive expression.
  - `coeffs` (numbers): the coefficients `c0`, `c1`, ... as separate arguments (literals or symbols), not as a collection. Any numeric type is accepted and converted to a double.

  Returns a double. Without coefficients the result is `0.0`; with a single coefficient it is that coefficient as a double, whatever `x` is. An exactly zero result is always `+0.0`. For an infinite or NaN `x` the result follows floating point arithmetic, so `##NaN` gives `##NaN`.

  The result is identical, bit for bit, to the one from [[evalpoly]] and [[makepoly]] for the same coefficients.

  See also [[evalpoly]], [[makepoly]], [[mevalpoly-complex]], [[mevalpoly-scalar-complex]]."
  [x & coeffs]
  (let [cnt (count coeffs)]
    (condp m/== cnt
      0 0.0
      1 `(double ~(first coeffs))
      2 (let [[z y] coeffs] `(m/muladd ~x ~y ~z))
      `(m/muladd ~x (mevalpoly ~x ~@(rest coeffs)) ~(first coeffs)))))

(defn evalpoly
  "Evaluates a real polynomial at `x` for the given coefficients in ascending order of power.

  The polynomial is `c0 + c1*x + c2*x^2 + ...`. The roundoff error is bounded by about `2*n*eps` times the sum of `|ci*x^i|` for `n` coefficients, so the result is accurate except close to a root.

  Parameters:

  - `x` (number): the point of evaluation, converted to a double.
  - `coeffs` (numbers): the coefficients `c0`, `c1`, ... as separate arguments; use `apply` for a collection. Any numeric type is accepted.

  Returns a double. Without coefficients the result is `0.0`; with a single coefficient it is that coefficient as a double, whatever `x` is. An exactly zero result is always `+0.0`. For an infinite or NaN `x` the result follows floating point arithmetic, so `##NaN` gives `##NaN` and a product of an infinity and a zero gives `##NaN`.

  Calls with explicit coefficients give the same result as [[mevalpoly]] and as the function from [[makepoly]], bit for bit.

  See also [[mevalpoly]], [[makepoly]], [[evalpoly-complex]], [[evalpoly-scalar-complex]], [[polynomial]]."
  {:inline (fn [x & coeffs] `(let [x# ~x] (mevalpoly x# ~@coeffs)))
   :inline-arities (fn [^long a] (m/>= a 1))}
  [x & coeffs]
  (if-not (seq coeffs)
    0.0
    (let [rc (reverse coeffs)
          xx (double x)]
      (loop [rcoeffs (rest rc)
             ex (double (first rc))]
        (if-not (seq rcoeffs)
          ex
          (recur (rest rcoeffs)
                 (m/muladd xx ex (double (first rcoeffs)))))))))

(defn makepoly
  "Creates a function evaluating the real polynomial with the given coefficients.

  The polynomial is `c0 + c1*x + c2*x^2 + ...`.

  Parameters:

  - `coeffs` (sequence of numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Any numeric type is accepted.

  Returns a function of one number `x` that returns a double equal to `(apply evalpoly x coeffs)`. An empty or `nil` sequence gives a function that always returns `0.0`; a single coefficient gives a function that always returns that coefficient as a double. The coefficients are fixed when the function is created.

  See also [[evalpoly]], [[mevalpoly]], [[makepoly-complex]], [[makepoly-scalar-complex]], [[polynomial]]."
  [coeffs]
  (cond
    (not (seq coeffs)) (constantly 0.0)
    (m/== 1 (count coeffs)) (constantly (double (first coeffs)))
    :else (let [rc (reverse coeffs)]
            (fn [^double x]
              (loop [rcoeffs (rest rc)
                     ex (double (first rc))]
                (if-not (seq rcoeffs)
                  ex
                  (recur (rest rcoeffs)
                         (m/muladd x ex (double (first rcoeffs))))))))))

;; complex

;; polynomials

(defmacro mevalpoly-complex
  "Evaluates a complex polynomial at `z` for coefficients given explicitly in the code; the macro version of [[evalpoly-complex]].

  The polynomial is `c0 + c1*z + c2*z^2 + ...`, with coefficients in ascending order of power.

  Parameters:

  - `z` (`Vec2`): the point of evaluation. It can be evaluated several times, so pass a symbol or a literal rather than an expensive expression.
  - `coeffs` (`Vec2`s): the complex coefficients `c0`, `c1`, ... as separate arguments, not as a collection.

  The arguments must already be complex numbers (`Vec2`, for example from `fastmath.complex/complex`); real numbers are not converted and fail at runtime with a `ClassCastException`. Use [[evalpoly-complex]] when they may be real numbers.

  Returns a `Vec2`. Without coefficients the result is `fastmath.complex/ZERO`; with a single coefficient it is that coefficient.

  See also [[evalpoly-complex]], [[makepoly-complex]], [[mevalpoly-scalar-complex]], [[mevalpoly]]."
  [z & coeffs]
  (let [cnt (count coeffs)]
    (case (unchecked-int cnt)
      0 `cplx/ZERO
      1 (first coeffs)
      2 (let [[c0 c1] coeffs]
          `(cplx/muladd ~z ~c1 ~c0))
      `(cplx/muladd ~z (mevalpoly-complex ~z ~@(rest coeffs)) ~(first coeffs)))))

(defmacro mevalpoly-scalar-complex
  "Evaluates a complex polynomial with real coefficients at `z` for coefficients given explicitly in the code; the macro version of [[evalpoly-scalar-complex]].

  The polynomial is `c0 + c1*z + c2*z^2 + ...`, with coefficients in ascending order of power.

  Parameters:

  - `z` (`Vec2`): the point of evaluation. It can be evaluated several times, so pass a symbol or a literal rather than an expensive expression.
  - `coeffs` (numbers): the real coefficients `c0`, `c1`, ... as separate arguments (literals or symbols), not as a collection. Any numeric type is accepted and converted to a double.

  The point `z` must already be a complex number (`Vec2`); a real number is not converted and fails at runtime with a `ClassCastException`. Use [[evalpoly-scalar-complex]] when it may be a real number.

  Returns a `Vec2`. Without coefficients the result is `fastmath.complex/ZERO`; with a single coefficient it is that coefficient with a zero imaginary part. Unlike [[evalpoly-scalar-complex]] it works with plain complex arithmetic only, so huge values of `z` are handled the same way as in [[mevalpoly-complex]].

  See also [[evalpoly-scalar-complex]], [[makepoly-scalar-complex]], [[mevalpoly-complex]], [[mevalpoly]]."
  [z & coeffs]
  (let [cnt (count coeffs)]
    (condp clojure.core/= cnt
      0 `cplx/ZERO
      1 `(Vec2. (double ~(first coeffs)) 0.0)
      `(cplx/muladd ~z (mevalpoly-scalar-complex ~z ~@(rest coeffs)) (Vec2. (double ~(first coeffs)) 0.0)))))

(defn evalpoly-complex
  "Evaluates a complex polynomial at `z` for the given coefficients in ascending order of power.

  The polynomial is `c0 + c1*z + c2*z^2 + ...`. The roundoff error is bounded by about `4*n` times the machine epsilon times the sum of `|ci*z^i|` for `n` coefficients.

  Parameters:

  - `z` (complex number): the point of evaluation; a real number is treated as a complex number with a zero imaginary part.
  - `coeffs` (complex numbers): the coefficients `c0`, `c1`, ... as separate arguments; use `apply` for a collection. Each can be a `Vec2` or a real number.

  Returns a `Vec2`. Without coefficients the result is `fastmath.complex/ZERO`; with a single coefficient it is that coefficient as a `Vec2`, whatever `z` is. For infinite or NaN parts the result follows complex floating point arithmetic, so for example an infinite real `z` can produce `##NaN` in the imaginary part.

  See also [[evalpoly-scalar-complex]] (real coefficients), [[mevalpoly-complex]], [[makepoly-complex]], [[evalpoly]]."
  [z & coeffs]
  (if-not (seq coeffs)
    cplx/ZERO
    (let [z (cplx/ensure-complex z)
          rc (reverse (map cplx/ensure-complex coeffs))]
      (loop [rcoeffs (rest rc)
             ex (first rc)]
        (if-not (seq rcoeffs)
          ex
          (recur (rest rcoeffs)
                 (cplx/muladd z ex (first rcoeffs))))))))

(defn- scalar-complex-horner
  "Value at complex `z` of the polynomial with real coefficients `rc`, a non-empty vector in descending order (highest degree first).

  Uses the squared modulus of `z`; when that is not a finite number (`z` is huge, infinite or NaN) falls back to plain complex Horner evaluation, which does not need it."
  ^Vec2 [rc ^Vec2 z]
  (let [cnt (dec (count rc))
        q (cplx/norm z)]
    (cond
      (m/zero? cnt) (Vec2. (double (rc 0)) 0.0)
      (m/invalid-double? q) (apply evalpoly-complex z (rseq rc))
      :else (let [p (m/* -2.0 (.x z))]
              (loop [i (long 1)
                     r (double (rc 0))
                     s 0.0]
                (if (m/== i cnt)
                  (Vec2. (m/- (m/+ (double (rc i))
                                   (m/* (.x z) r))
                              (m/* q s))
                         (m/* (.y z) r))
                  (recur (m/inc i)
                         (m/- (double (rc i)) (m/* p r) (m/* q s))
                         r)))))))

(defn evalpoly-scalar-complex
  "Evaluates a complex polynomial with real coefficients at `z` for the given coefficients in ascending order of power.

  The polynomial is `c0 + c1*z + c2*z^2 + ...`. It gives the same result as [[evalpoly-complex]], within roundoff, with less work when all coefficients are real. The roundoff error is bounded by about `8*n` times the machine epsilon times the sum of `|ci*z^i|` for `n` coefficients.

  Parameters:

  - `z` (complex number): the point of evaluation; a real number is treated as a complex number with a zero imaginary part.
  - `coeffs` (numbers): the real coefficients `c0`, `c1`, ... as separate arguments; use `apply` for a collection. Any numeric type is accepted.

  Returns a `Vec2`. Without coefficients the result is `fastmath.complex/ZERO`; with a single coefficient it is that coefficient with a zero imaginary part, whatever `z` is. It is also valid for huge `z` (modulus above about 1e154), where the polynomial value is still representable, and for infinite or NaN parts, where it returns what [[evalpoly-complex]] returns.

  See also [[evalpoly-complex]], [[mevalpoly-scalar-complex]], [[makepoly-scalar-complex]], [[evalpoly]]."
  [z & coeffs]
  (if-not (seq coeffs)
    cplx/ZERO
    (scalar-complex-horner (vec (reverse coeffs)) (cplx/ensure-complex z))))

(defn makepoly-complex
  "Creates a function evaluating the complex polynomial with the given coefficients.

  The polynomial is `c0 + c1*z + c2*z^2 + ...`.

  Parameters:

  - `coeffs` (sequence of complex numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Each can be a `Vec2` or a real number.

  Returns a function of one complex number `z` (a real number is accepted) that returns a `Vec2` equal to `(apply evalpoly-complex z coeffs)`. An empty or `nil` sequence gives a function that always returns `fastmath.complex/ZERO`; a single coefficient gives a function that always returns that coefficient as a `Vec2`. The coefficients are fixed when the function is created.

  See also [[evalpoly-complex]], [[mevalpoly-complex]], [[makepoly-scalar-complex]], [[makepoly]]."
  [coeffs]
  (let [coeffs (map cplx/ensure-complex coeffs)]
    (cond
      (not (seq coeffs)) (constantly cplx/ZERO)
      (= 1 (count coeffs)) (constantly (first coeffs))
      :else (let [rc (reverse coeffs)]
              (fn [x]
                (let [x (cplx/ensure-complex x)]
                  (loop [rcoeffs (rest rc)
                         ex (first rc)]
                    (if-not (seq rcoeffs)
                      ex
                      (recur (rest rcoeffs)
                             (cplx/muladd x ex (first rcoeffs)))))))))))

(defn makepoly-scalar-complex
  "Creates a function evaluating the complex polynomial with the given real coefficients.

  The polynomial is `c0 + c1*z + c2*z^2 + ...`.

  Parameters:

  - `coeffs` (sequence of numbers): the real coefficients `c0`, `c1`, ... in ascending order of power. Any numeric type is accepted.

  Returns a function of one complex number `z` (a real number is accepted) that returns a `Vec2` equal to `(apply evalpoly-scalar-complex z coeffs)`. An empty or `nil` sequence gives a function that always returns `fastmath.complex/ZERO`; a single coefficient gives a function that always returns that coefficient with a zero imaginary part. The coefficients are fixed when the function is created.

  See also [[evalpoly-scalar-complex]], [[mevalpoly-scalar-complex]], [[makepoly-complex]], [[makepoly]]."
  [coeffs]
  (cond
    (not (seq coeffs)) (constantly cplx/ZERO)
    (= 1 (count coeffs)) (constantly (cplx/complex (first coeffs)))
    :else (let [rc (vec (reverse coeffs))]
            (fn [z] (scalar-complex-horner rc (cplx/ensure-complex z))))))

;;

(defn- degrees->vars
  ([^long v] (degrees->vars v '()))
  ([^long v buff]
   (if (m/neg? v)
     buff
     (recur (m/dec v) (conj buff (case (int v)
                                   0 ""
                                   1 "x"
                                   (str "x^" v)))))))

(defn- polynomial->str
  [coeffs ^long degree]
  (let [^DecimalFormat f (doto (DecimalFormat. "0.####" (java.text.DecimalFormatSymbols. java.util.Locale/ROOT))
                           (.setPositivePrefix "+")
                           (.setNegativePrefix "-"))
        ;; a negative coefficient that rounds to zero is formatted as -0
        format-coefficient (fn [^double v] (let [s (.format f v)] (if (= s "-0") "+0" s)))
        numbers+var (apply str (take 20 (interleave (map format-coefficient coeffs)
                                                    (degrees->vars degree))))
        res (str "#polynomial{" degree "}(x) = " (if (str/starts-with?  numbers+var "+")
                                                   (subs numbers+var 1) numbers+var))]
    (if (m/> degree 10) (str res "+...") res))  )


(defn- throw-incompatible-polynomials
  [p1 p2]
  (let [type-name (fn [p] (if (nil? p) "nil" (.getSimpleName (class p))))]
    (throw (IllegalArgumentException.
            (str (type-name p1) " and " (type-name p2)
                 " cannot be combined; both must be polynomials of the same kind"
                 " (from polynomial or from ratio-polynomial)")))))

(deftype Polynomial [^doubles cfs ^long d]
  Object
  (toString [_] (polynomial->str cfs d))
  (equals [_ poly]
    (and (instance? Polynomial poly)
         (java.util.Arrays/equals cfs ^doubles (.cfs ^Polynomial poly))))
  (hashCode [_]
    (mix-collection-hash (java.util.Arrays/hashCode cfs) d))
  IFn
  (invoke [p v] (prot/evaluate p v))
  (applyTo [p args]
    (if (and args (nil? (next args)))
      (prot/evaluate p (first args))
      (throw (clojure.lang.ArityException. (int (count args)) "Polynomial"))))
  prot/PolynomialProto
  (degree [_] d)
  (coeffs [_] (seq cfs))
  (add [p1 p2]
    (when-not (instance? Polynomial p2) (throw-incompatible-polynomials p1 p2))
    (let [[^Polynomial poly-min ^Polynomial poly-max] (if (m/> (.d ^Polynomial p1)
                                                               (.d ^Polynomial p2))
                                                        [p2 p1] [p1 p2])
          ^doubles target (double-array (.cfs poly-max))
          ^doubles source (.cfs poly-min)]
      (dotimes [i (m/inc (.d poly-min))]
        (Array/add target i (Array/aget source i)))
      (Polynomial. target (.d poly-max))))
  (negate [_] (Polynomial. (v/sub cfs) d))
  (scale [_ v] (Polynomial. (v/mult cfs v) d))
  (mult [p1 p2]
    (when-not (instance? Polynomial p2) (throw-incompatible-polynomials p1 p2))
    (let [^Polynomial p2 p2]
      (cond
        (m/zero? d) (prot/scale p2 (Array/aget cfs 0))
        (m/zero? (.d p2)) (prot/scale p1 (Array/aget ^doubles (.cfs p2) 0))
        :else (let [nd (m/+ d (.d p2))
                    target (double-array (m/inc nd))]
                (doseq [^long i1 (range (m/inc d))
                        ^long i2 (range (m/inc (.d p2)))]
                  (Array/add target (m/+ i1 i2) (m/* (Array/aget cfs i1)
                                                     (Array/aget ^doubles (.cfs p2) i2))))
                (Polynomial. target nd)))))
  (derivative [p order]
    (let [order (long order)]
      (cond
        (m/neg? order) (throw (IllegalArgumentException.
                               (str "Derivative order must not be negative, got " order)))
        (m/zero? order) p
        (m/> order d) (Polynomial. (double-array 1) 0)
        :else (let [size (m/inc (m/- d order))
                    ^doubles target (double-array size)]
                (loop [i (long 0)
                       pos order
                       fact (m/factorial order)]
                  (if (m/== i size)
                    (Polynomial. target (m/dec size))
                    (let [i+ (m/inc i)
                          pos+ (m/inc pos)]
                      (Array/aset target i (m/* fact (Array/aget cfs pos)))
                      (recur i+ pos+ (m/* (m// fact i+) pos+)))))))))
  (evaluate [_ x]
    (loop [i d
           ex (Array/aget cfs i)]
      (if (m/zero? i)
        ex
        (let [i- (m/dec i)]
          (recur i- (m/muladd x ex (Array/aget cfs i-))))))))

(set! *unchecked-math* false)

(deftype PolynomialR [cfs ^long d]
  Object
  (toString [_] (polynomial->str cfs d))
  (equals [_ poly]
    (and (instance? PolynomialR poly)
         (= cfs (.cfs ^PolynomialR poly))))
  (hashCode [_]
    (mix-collection-hash (hash cfs) d))
  IFn
  (invoke [p v] (prot/evaluate p v))
  (applyTo [p args]
    (if (and args (nil? (next args)))
      (prot/evaluate p (first args))
      (throw (clojure.lang.ArityException. (int (count args)) "PolynomialR"))))
  prot/PolynomialProto
  (degree [_] d)
  (coeffs [_] cfs)
  (add [p1 p2]
    (when-not (instance? PolynomialR p2) (throw-incompatible-polynomials p1 p2))
    (let [[^PolynomialR poly-min ^PolynomialR poly-max] (if (m/> (.d ^PolynomialR p1)
                                                                 (.d ^PolynomialR p2))
                                                          [p2 p1] [p1 p2])
          target (transient (vec (.cfs poly-max)))
          res (->> (.coeffs poly-min)
                   (reduce (fn [[^long id t] s]
                             [(m/inc id) (assoc! t id (+' (t id) s))]) [0 target])
                   (second)
                   (persistent!))]
      (PolynomialR. res (.d poly-max))))
  (negate [_] (PolynomialR. (mapv -' cfs) d))
  (scale [_ v] (let [rv (rationalize v)]
                 (PolynomialR. (mapv (fn [v] (*' v rv)) cfs) d)))
  (mult [p1 p2]
    (when-not (instance? PolynomialR p2) (throw-incompatible-polynomials p1 p2))
    (let [^PolynomialR p2 p2]
      (cond
        (m/zero? d) (prot/scale p2 (cfs 0))
        (m/zero? (.d p2)) (prot/scale p1 ((.cfs p2) 0))
        :else (let [nd (m/long-add d (.d p2))
                    target (transient (vec (repeat (m/long-inc nd) 0)))
                    res (->> (for [i1 (range (m/inc d))
                                   i2 (range (m/inc (.d p2)))]
                               [i1 i2])
                             (reduce (fn [t [^long i1 ^long i2]]                                       
                                       (let [pos (m/+ i1 i2)]
                                         (assoc! t pos (+' (t pos) (*' (cfs i1) ((.cfs p2) i2)))))) target)
                             (persistent!))]
                (PolynomialR. res nd)))))
  (derivative [p order]
    (let [order (long order)]
      (cond
        (m/neg? order) (throw (IllegalArgumentException.
                               (str "Derivative order must not be negative, got " order)))
        (m/zero? order) p
        (m/> order d) (PolynomialR. [0] 0)
        :else (let [size (m/long-inc (m/long-sub d order))]
                (loop [i (long 0)
                       pos order
                       fact (reduce *' 1 (range 1 (m/long-inc order)))
                       target (transient (vec (repeat size 0)))]
                  (if (m/== i size)
                    (PolynomialR. (persistent! target) (m/dec size))
                    (let [i+ (m/inc i)
                          pos+ (m/inc pos)]
                      (recur i+ pos+ (*' (/ fact i+) pos+)
                             (assoc! target i (*' fact (cfs pos)))))))))))
  (evaluate [_ x]
    (let [dx (double x)]
      (if (m/invalid-double? dx)
        ;; no rational number for NaN and infinities: plain double evaluation
        (loop [i d
               ex (double (cfs i))]
          (if (m/zero? i)
            ex
            (let [i- (m/dec i)]
              (recur i- (m/muladd dx ex (double (cfs i-)))))))
        (let [rx (rationalize x)]
          (loop [i d
                 ex (cfs i)]
            (if (m/zero? i)
              ex
              (let [i- (m/dec i)]
                (recur i- (+ (* rx ex) (cfs i-)))))))))))

(set! *unchecked-math* :warn-on-boxed)

(alter-meta! #'->Polynomial assoc :doc
             "Positional constructor of `Polynomial`, the polynomial object with double coefficients. Use [[polynomial]] instead.

  Parameters: `cfs` (array of doubles): the coefficients in ascending order of power, at least one; `d` (long): the nominal degree, which must equal the number of coefficients minus one. Neither is validated.")

(alter-meta! #'->PolynomialR assoc :doc
             "Positional constructor of `PolynomialR`, the polynomial object with exact rational coefficients. Use [[ratio-polynomial]] instead.

  Parameters: `cfs` (vector of rational numbers): the coefficients in ascending order of power, at least one; `d` (long): the nominal degree, which must equal the number of coefficients minus one. Neither is validated.")

(def ^:private RONE (PolynomialR. [1] 0))

(defmethod print-method Polynomial [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method PolynomialR [v ^java.io.Writer w] (.write w (str v)))

(defn polynomial
  "Creates a polynomial object with double coefficients, from the coefficients or from points to interpolate.

  The object is a `Polynomial`: `c0 + c1*x + c2*x^2 + ...`, with coefficients stored as doubles in ascending order of power. It can be called as a function of one number, and used with [[add]], [[sub]], [[scale]], [[mult]], [[derivative]], [[evaluate]], [[coeffs]] and [[degree]].

  Parameters:

  - `coeffs` (sequence of numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Any numeric type is accepted and converted to a double. An empty or `nil` sequence gives the zero polynomial (degree 0, the single coefficient `0.0`).
  - `xs`, `ys` (sequences of numbers): the abscissae and ordinates of the points. The result is the interpolating polynomial of degree `n-1` for `n` points.

  Returns a `Polynomial`. Its degree is nominal: the number of coefficients minus one, so trailing zero coefficients are kept and no operation removes them. Two polynomials are equal when their coefficient arrays are equal as doubles, compared bit by bit (so `0.0` and `-0.0` differ), and they have equal hashes. A `Polynomial` is never equal to a `PolynomialR`.

  Interpolation needs at least two points, `xs` and `ys` of the same length and distinct `xs`; otherwise an `IllegalArgumentException` is thrown. The coefficients are computed in double precision, so the result is ill-conditioned for many points or for points close to each other.

  See also [[ratio-polynomial]] (exact rational coefficients), [[coeffs->polynomial]], [[evalpoly]], [[makepoly]]."
  (^Polynomial [xs ys]
   (polynomial (PolInterp/getCoefficients (m/seq->double-array xs) (m/seq->double-array ys))))
  (^Polynomial [coeffs]
   (if (seq coeffs)
     (Polynomial. (double-array coeffs) (m/dec (count coeffs)))
     (Polynomial. (double-array 1) 0))))

(defn ratio-polynomial
  "Creates a polynomial object with exact rational coefficients, from the coefficients or from points to interpolate.

  The object is a `PolynomialR`: `c0 + c1*x + c2*x^2 + ...`, with coefficients stored as ratios (arbitrary precision integers and fractions) in ascending order of power. [[add]], [[sub]], [[scale]], [[mult]] and [[derivative]] are exact for any degree and order. It can be called as a function of one number and used with [[evaluate]], [[coeffs]] and [[degree]].

  Parameters:

  - `coeffs` (sequence of numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Each is converted with `rationalize`, so a double becomes the exact decimal number it prints as (`0.1` gives `1/10`), not its binary value. An empty or `nil` sequence gives the zero polynomial (degree 0, the single coefficient `0`). A NaN or infinite coefficient throws an `IllegalArgumentException`.
  - `xs`, `ys` (sequences of numbers): the abscissae and ordinates of the points. The interpolating polynomial is computed in double precision, as for [[polynomial]], and then converted to ratios, so only the operations on the result are exact.

  Returns a `PolynomialR`. Its degree is nominal: the number of coefficients minus one, and no operation removes trailing zero coefficients. Two polynomials are equal when their coefficients are equal as numbers (`1` and `1N` are equal), and they have equal hashes. A `PolynomialR` is never equal to a `Polynomial`.

  Calling the object with a finite number returns the exact ratio; with NaN or an infinity it returns a double, as [[polynomial]] does.

  See also [[polynomial]], [[coeffs->ratio-polynomial]]."
  (^PolynomialR [xs ys]
   (ratio-polynomial (PolInterp/getCoefficients (m/seq->double-array xs) (m/seq->double-array ys))))
  (^PolynomialR [coeffs]
   (if (seq coeffs)
     (PolynomialR. (mapv rationalize coeffs) (m/dec (count coeffs)))
     (PolynomialR. [0] 0))))

(defn coeffs->polynomial
  "Creates a polynomial object with double coefficients from coefficients given as separate arguments.

  Parameters: `coeffs` (numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Without arguments the result is the zero polynomial.

  Returns a `Polynomial`, the same as `(polynomial coeffs)`.

  See also [[polynomial]], [[coeffs->ratio-polynomial]]."
  [& coeffs] (polynomial coeffs))

(defn coeffs->ratio-polynomial
  "Creates a polynomial object with exact rational coefficients from coefficients given as separate arguments.

  Parameters: `coeffs` (numbers): the coefficients `c0`, `c1`, ... in ascending order of power. Without arguments the result is the zero polynomial. A NaN or infinite coefficient throws an `IllegalArgumentException`.

  Returns a `PolynomialR`, the same as `(ratio-polynomial coeffs)`.

  See also [[ratio-polynomial]], [[coeffs->polynomial]]."
  [& coeffs] (ratio-polynomial coeffs))

(defn add
  "Adds polynomials.

  Parameters: `poly`, `poly1`, `poly2` (polynomial objects): from [[polynomial]] or from [[ratio-polynomial]]. Both arguments of the two-argument form must be of the same kind; combining a `Polynomial` with a `PolynomialR`, or with anything that is not a polynomial object, throws an `IllegalArgumentException`.

  Returns the polynomial itself for one argument, otherwise a new polynomial whose nominal degree is the larger of the two degrees (cancellation of leading coefficients does not lower it). The sum of `PolynomialR` objects is exact; for `Polynomial` objects each coefficient has one rounding error.

  See also [[sub]], [[mult]], [[scale]]."
  ([poly] poly)
  ([poly1 poly2] (prot/add poly1 poly2)))

(defn sub
  "Subtracts polynomials, or negates a polynomial.

  Parameters: `poly`, `poly1`, `poly2` (polynomial objects): from [[polynomial]] or from [[ratio-polynomial]]. Both arguments of the two-argument form must be of the same kind; combining a `Polynomial` with a `PolynomialR`, or with anything that is not a polynomial object, throws an `IllegalArgumentException`.

  Returns the negated polynomial for one argument, otherwise `poly1 - poly2` as a new polynomial whose nominal degree is the larger of the two degrees (cancellation of leading coefficients does not lower it, so subtracting a polynomial from itself gives zero coefficients with the same degree). The difference of `PolynomialR` objects is exact.

  See also [[add]], [[mult]], [[scale]]."
  ([poly] (prot/negate poly))
  ([poly1 poly2]
   (prot/add poly1 (prot/negate poly2))))

(defn scale
  "Multiplies a polynomial by a number.

  Parameters:

  - `poly` (polynomial object): from [[polynomial]] or from [[ratio-polynomial]].
  - `v` (number): the factor. For a `PolynomialR` it is converted with `rationalize` (a double becomes the exact decimal number it prints as), so a NaN or infinite `v` throws an `IllegalArgumentException`; for a `Polynomial` any double is accepted and NaN or infinity propagate into the coefficients.

  Returns a new polynomial of the same kind and the same nominal degree; scaling by zero gives zero coefficients, not a lower degree. The product with a `PolynomialR` is exact.

  See also [[mult]], [[add]]."
  [poly v] (prot/scale poly v))

(defn mult
  "Multiplies polynomials.

  Parameters: `poly`, `poly1`, `poly2` (polynomial objects): from [[polynomial]] or from [[ratio-polynomial]]. Both arguments of the two-argument form must be of the same kind; combining a `Polynomial` with a `PolynomialR`, or with anything that is not a polynomial object, throws an `IllegalArgumentException`.

  Returns the polynomial itself for one argument, otherwise the product as a new polynomial whose nominal degree is the sum of the two degrees (zero coefficients keep the degree). The product of `PolynomialR` objects is exact; each coefficient of a `Polynomial` product is a sum of products of double coefficients.

  See also [[scale]], [[add]], [[sub]]."
  ([poly] poly)
  ([poly1 poly2]
   (prot/mult poly1 poly2)))

(defn coeffs
  "Returns the coefficients of a polynomial object in ascending order of power: `c0`, `c1`, ..., `cn`.

  Parameters: `poly` (polynomial object): from [[polynomial]] or from [[ratio-polynomial]].

  Returns a sequence of doubles for a `Polynomial` and a vector of ratios (integers and fractions) for a `PolynomialR`. The number of coefficients is [[degree]] plus one, and it is at least one: the zero polynomial has the single coefficient zero. Trailing zero coefficients are included.

  See also [[degree]], [[coeffs->polynomial]]."
  [poly] (prot/coeffs poly))

(defn degree
  "Returns the nominal degree of a polynomial object: the number of its coefficients minus one.

  Parameters: `poly` (polynomial object): from [[polynomial]] or from [[ratio-polynomial]].

  The degree is nominal because leading zero coefficients are counted and no operation removes them: `[1 2 0 0]` has degree 3, the difference of a polynomial and itself keeps its degree, and the zero polynomial has degree 0. [[derivative]] lowers the degree by the order, and [[mult]] adds the degrees.

  See also [[coeffs]]."
  ^long [poly] (prot/degree poly))

(defn derivative
  "Returns the derivative of a polynomial object, optionally of a higher order.

  Parameters:

  - `poly` (polynomial object): from [[polynomial]] or from [[ratio-polynomial]].
  - `order` (non-negative integer, default 1): the order of the derivative. A non-integer is truncated; a negative order throws an `IllegalArgumentException`.

  Returns a polynomial of the same kind. Order 0 returns `poly` itself. For an order not above the degree the result has the nominal degree `degree - order`. For a larger order, and for any derivative of a constant or of the zero polynomial, the result is the zero polynomial (degree 0). The derivative of a `PolynomialR` is exact for any order; the coefficients of a `Polynomial` have a relative error of a few units of roundoff.

  See also [[degree]], [[evaluate]]."
  ([poly] (derivative poly 1))
  ([poly order] (prot/derivative poly order)))

(defn evaluate
  "Evaluates a polynomial object at a point and returns a double.

  Parameters:

  - `poly` (polynomial object): from [[polynomial]] or from [[ratio-polynomial]].
  - `x` (number): the point of evaluation, converted to a double.

  For a `PolynomialR` the value is computed exactly at the decimal value of `x` and rounded to a double once. For NaN or infinite `x` both kinds of polynomial follow floating point arithmetic and give the same result. A polynomial object can also be called directly as a function of exactly one number (also through `apply` with one argument); then a `Polynomial` returns a double and a `PolynomialR` returns the exact ratio for a finite argument. No other call form is supported.

  The monomial form loses accuracy for high degrees and close to a root. Use the `eval-*` functions of the orthogonal families (for example [[eval-legendre-P]]) for accurate values of those polynomials.

  See also [[evalpoly]], [[derivative]]."
  ^double [poly ^double x]
  (prot/evaluate poly x))

;; Orthogonal polynomials

(defn- check-degree!
  "Throws `IllegalArgumentException` when `degree` is negative (no such polynomial exists) or not below
  `Integer/MAX_VALUE` (the evaluators dispatch on the degree as an int, which would wrap around)."
  [^long degree]
  (cond
    (m/neg? degree)
    (throw (IllegalArgumentException. (str "Degree must not be negative, got " degree)))
    (m/>= degree Integer/MAX_VALUE)
    (throw (IllegalArgumentException. (str "Degree must be below " Integer/MAX_VALUE ", got " degree)))))

(defn eval-bernstein
  ^double [^long degree ^long order ^double x]
  (case (int degree)
    0 1.0
    1 (if (m/zero? order) (m/- 1.0 x) x)
    (m/* (m/combinations degree order) (m/fpow x order) (m/fpow (m/- 1.0 x) (m/long-sub degree order)))))


(defn bernstein
  [^long degree ^long order]
  (->> (range (m/inc degree))
       (map (fn [^long l]
              (if (m/< l order)
                0.0
                (m/* (if (m/even? (m/long-sub l order)) 1.0 -1.0)
                     (m/combinations degree l)
                     (m/combinations l order)))))
       (polynomial)))

;;

(defn eval-laguerre-L
  "Evaluate generalized Laguerre polynomial"
  (^double [^long degree ^double x] (eval-laguerre-L degree 0.0 x))
  (^double [^long degree ^double order ^double x]
   (case (int degree)
     0 1.0
     1 (m/- (m/inc order) x)
     (loop [i (long 2)
            pprev 1.0
            prev (m/- (m/inc order) x)]
       (if (m/> i degree)
         prev
         (recur (m/inc i) prev
                (m// (m/- (m/* (m/+ order (m/- (m/* 2.0 i) 1.0 x)) prev)
                          (m/* (m/+ order (m/dec i)) pprev)) i)))))))

(defn laguerre-L-ratio
  [^long degree ^double order]
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [(m/inc order) -1])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [(m/inc order) -1])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (scale (sub (mult prev (ratio-polynomial [(m/+ order (m/* 2 i) -1) -1]))
                           (scale pprev (m/dec (m/+ order i)))) (/ 1 i)))))))

(defn laguerre-L
  "Generalized Laguerre polynomials"
  ([^long degree] (laguerre-L degree 0.0))
  ([^long degree ^double order]
   (polynomial (coeffs (laguerre-L-ratio degree order)))))

;;

(defn eval-chebyshev-T
  "Evaluates the Chebyshev polynomial of the first kind `T_n` at `x`.

  `T_n(cos(a)) = cos(n*a)`. The polynomials are orthogonal on `[-1, 1]` with the weight `1/sqrt(1 - x^2)` and satisfy `T_0 = 1`, `T_1 = x` and `T_(n+1) = 2*x*T_n - T_(n-1)`. They are polynomials in `x` for every real `x`: outside `[-1, 1]` the value is the continuation of the polynomial, `T_n(x) = cosh(n*acosh(x))` for `x > 1`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree and an infinite `x` gives an infinity (negative for `##-Inf` and an odd degree). The error is a few units of roundoff times `n` inside `[-1, 1]` (absolute, for values of at most 1) and a few units times `n*log(2|x|)` (relative) outside. This form stays accurate for any degree, unlike evaluating the coefficients from [[chebyshev-T]].

  See also [[chebyshev-T]], [[chebyshev-T-ratio]], [[eval-chebyshev-U]], [[eval-chebyshev-V]], [[eval-chebyshev-W]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (case (int degree)
    0 1.0
    1 x
    2 (m/dec (m/* 2.0 x x))
    3 (m/* x (m/- (m/* 4.0 x x) 3.0))
    4 (let [x2 (m/* x x)] (m/inc (m/* 8.0 x2 (m/dec x2))))
    (cond
      (m/> x 1.0) (m/cosh (m/* degree (m/acosh x)))
      (m/< x -1.0) (m/* (m/fpow -1.0 degree) (m/cosh (m/* degree (m/acosh (m/- x)))))
      :else (m/cos (* degree (m/acos x))))))

(defn chebyshev-T-ratio
  "Creates the Chebyshev polynomial of the first kind `T_n` with exact integer coefficients.

  `T_0 = 1`, `T_1 = x` and `T_(n+1) = 2*x*T_n - T_(n-1)`; see [[eval-chebyshev-T]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients grow like `2^n` and alternate in sign.

  See also [[chebyshev-T]] (double coefficients), [[eval-chebyshev-T]] (direct evaluation)."
  [^long degree]
  (check-degree! degree)
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [0 1])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [0 1])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-T
  "Creates the Chebyshev polynomial of the first kind `T_n` as a polynomial object with double coefficients.

  See [[eval-chebyshev-T]] for the definition and [[chebyshev-T-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-T]] for values at a high degree.

  See also [[chebyshev-T-ratio]], [[eval-chebyshev-T]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-T-ratio degree))))

(defn- chebyshev-U-above-one
  "Value at `x` > 1 (or at infinity) of the second kind Chebyshev polynomial of a degree of at least 1, as
  exp(n t) (1 - exp(-2(n+1) t)) / (1 - exp(-2t)) with t = acosh(x): no intermediate value exceeds the result."
  ^double [^long degree ^double x]
  (let [t (m/acosh x)]
    (m/* (m/exp (m/* degree t))
         (m// (m/expm1 (m/* -2.0 (m/inc degree) t))
              (m/expm1 (m/* -2.0 t))))))

(defn- chebyshev-U-within-one
  "Value at 0 <= `x` <= 1 (or at NaN) of the second kind Chebyshev polynomial of a degree of at least 1, as
  sin((n+1) acos(x)) / sqrt(1 - x^2), and as its two term expansion next to x = 1 where the quotient is 0/0."
  ^double [^long degree ^double x]
  (let [degree+ (m/inc degree)
        ;; (1 - x)(1 + x), not 1 - x^2: the difference 1 - x is exact next to 1
        near-one (m/* (m/- 1.0 x) (m/+ 1.0 x))]
    (if (m/< near-one (m// 1.0E-7 (m/* degree+ degree+)))
      (m/* degree+ (m/- 1.0 (m/* m/SIXTH degree (m/+ degree 2) near-one)))
      ;; the angle from its sine and cosine: acos(x) has a small absolute, but a large relative error next to x = 1;
      ;; java.lang.Math/sin, since m/sin has a small absolute, but not a small relative error for a small argument
      (let [sine (m/sqrt near-one)]
        (m// (Math/sin (m/* degree+ (m/atan2 sine x))) sine)))))

(defn eval-chebyshev-U
  "Evaluates the Chebyshev polynomial of the second kind `U_n` at `x`.

  `U_n(cos(a)) = sin((n+1)*a)/sin(a)`. The polynomials are orthogonal on `[-1, 1]` with the weight `sqrt(1 - x^2)` and satisfy `U_0 = 1`, `U_1 = 2*x` and `U_(n+1) = 2*x*U_n - U_(n-1)`. They are polynomials in `x` for every real `x`: outside `[-1, 1]` the value is the continuation of the polynomial, `U_n(x) = sinh((n+1)*acosh(x))/sinh(acosh(x))` for `x > 1`. `U_n(1) = n+1` and `U_n(-x) = (-1)^n*U_n(x)`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree and an infinite `x` gives an infinity (negative for `##-Inf` and an odd degree). The value is finite whenever the polynomial value is representable. The error is a few units of roundoff times `n` inside `[-1, 1]` (absolute, for values of at most `n+1`) and a few units times `n*log(2|x|)` (relative) outside. This form stays accurate for any degree, unlike evaluating the coefficients from [[chebyshev-U]].

  See also [[chebyshev-U]], [[chebyshev-U-ratio]], [[eval-chebyshev-T]], [[eval-chebyshev-V]], [[eval-chebyshev-W]], [[eval-gegenbauer-C]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (case (int degree)
    0 1.0
    1 (m/* 2.0 x)
    2 (m/dec (m/* 4.0 x x))
    3 (m/* 4.0 x (m/dec (m/* 2.0 x x)))
    4 (let [x2 (m/* x x)] (m/inc (m/* x2 (m/- (m/* 16.0 x2) 12.0))))
    ;; U(n, -x) = (-1)^n U(n, x): evaluate at |x|, where the trigonometric and hyperbolic forms are accurate
    (let [ax (m/abs x)
          u (if (m/> ax 1.0)
              (chebyshev-U-above-one degree ax)
              (chebyshev-U-within-one degree ax))]
      (if (and (m/neg? x) (m/odd? degree)) (m/- u) u))))

(defn chebyshev-U-ratio
  "Creates the Chebyshev polynomial of the second kind `U_n` with exact integer coefficients.

  `U_0 = 1`, `U_1 = 2*x` and `U_(n+1) = 2*x*U_n - U_(n-1)`; see [[eval-chebyshev-U]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients grow like `2^n` and alternate in sign.

  See also [[chebyshev-U]] (double coefficients), [[eval-chebyshev-U]] (direct evaluation)."
  [^long degree]
  (check-degree! degree)
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [0 2])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [0 2])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-U
  "Creates the Chebyshev polynomial of the second kind `U_n` as a polynomial object with double coefficients.

  See [[eval-chebyshev-U]] for the definition and [[chebyshev-U-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-U]] for values at a high degree.

  See also [[chebyshev-U-ratio]], [[eval-chebyshev-U]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-U-ratio degree))))

(defn eval-chebyshev-V
  "Evaluates the Chebyshev polynomial of the third kind `V_n` at `x`.

  `V_n(cos(a)) = cos((n+1/2)*a)/cos(a/2)`. The polynomials are orthogonal on `[-1, 1]` with the weight `sqrt((1+x)/(1-x))` and satisfy `V_0 = 1`, `V_1 = 2*x - 1` and `V_(n+1) = 2*x*V_n - V_(n-1)`. Also `V_n = U_n - U_(n-1)`, `V_n(1) = 1`, `V_n(-1) = (-1)^n*(2*n+1)` and `V_n(-x) = (-1)^n*W_n(x)`. They are polynomials in `x` for every real `x`, and the value is the polynomial for any argument outside `[-1, 1]` too.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree and an infinite `x` gives an infinity (negative for `##-Inf` and an odd degree). The error is a few units of roundoff times `n` inside `[-1, 1]` (absolute, for values of at most `2*n+1`), up to several times that next to `x = 1`, where the result is a difference of two nearly equal values, and a few units times `n*log(2|x|)` (relative) outside. For arguments just above 1 and degrees in the thousands the result can overflow to infinity slightly before the true value does.

  See also [[chebyshev-V]], [[chebyshev-V-ratio]], [[eval-chebyshev-W]], [[eval-chebyshev-U]], [[eval-chebyshev-T]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (if (m/zero? degree)
    1.0
    ;; V(n) = U(n) - U(n-1). V(n) overflows together with U(n), except just above x = 1, where V(n) is
    ;; smaller than U(n) by the factor 1 - exp(-acosh(x)); the result is then infinite a little too early
    (let [u (eval-chebyshev-U degree x)]
      (if (m/inf? u) u (m/- u (eval-chebyshev-U (m/long-dec degree) x))))))

(defn chebyshev-V-ratio
  "Creates the Chebyshev polynomial of the third kind `V_n` with exact integer coefficients.

  `V_0 = 1`, `V_1 = 2*x - 1` and `V_(n+1) = 2*x*V_n - V_(n-1)`; see [[eval-chebyshev-V]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument.

  See also [[chebyshev-V]] (double coefficients), [[eval-chebyshev-V]] (direct evaluation)."
  [^long degree]
  (check-degree! degree)
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [-1 2])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [-1 2])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-V
  "Creates the Chebyshev polynomial of the third kind `V_n` as a polynomial object with double coefficients.

  See [[eval-chebyshev-V]] for the definition and [[chebyshev-V-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-V]] for values at a high degree.

  See also [[chebyshev-V-ratio]], [[eval-chebyshev-V]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-V-ratio degree))))

(defn eval-chebyshev-W
  "Evaluates the Chebyshev polynomial of the fourth kind `W_n` at `x`.

  `W_n(cos(a)) = sin((n+1/2)*a)/sin(a/2)`. The polynomials are orthogonal on `[-1, 1]` with the weight `sqrt((1-x)/(1+x))` and satisfy `W_0 = 1`, `W_1 = 2*x + 1` and `W_(n+1) = 2*x*W_n - W_(n-1)`. Also `W_n = U_n + U_(n-1)`, `W_n(1) = 2*n+1`, `W_n(-1) = (-1)^n` and `W_n(-x) = (-1)^n*V_n(x)`. They are polynomials in `x` for every real `x`, and the value is the polynomial for any argument outside `[-1, 1]` too.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree and an infinite `x` gives an infinity (negative for `##-Inf` and an odd degree). The error is a few units of roundoff times `n` inside `[-1, 1]` (absolute, for values of at most `2*n+1`), up to several times that next to `x = -1`, where the result is a difference of two nearly equal values, and a few units times `n*log(2|x|)` (relative) outside. For arguments just below -1 and degrees in the thousands the result can overflow to infinity slightly before the true value does.

  See also [[chebyshev-W]], [[chebyshev-W-ratio]], [[eval-chebyshev-V]], [[eval-chebyshev-U]], [[eval-chebyshev-T]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (if (m/zero? degree)
    1.0
    ;; W(n) = U(n) + U(n-1); as for V, the result is infinite a little too early just below x = -1
    (let [u (eval-chebyshev-U degree x)]
      (if (m/inf? u) u (m/+ u (eval-chebyshev-U (m/long-dec degree) x))))))

(defn chebyshev-W-ratio
  "Creates the Chebyshev polynomial of the fourth kind `W_n` with exact integer coefficients.

  `W_0 = 1`, `W_1 = 2*x + 1` and `W_(n+1) = 2*x*W_n - W_(n-1)`; see [[eval-chebyshev-W]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument.

  See also [[chebyshev-W]] (double coefficients), [[eval-chebyshev-W]] (direct evaluation)."
  [^long degree]
  (check-degree! degree)
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [1 2])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [1 2])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-W
  "Creates the Chebyshev polynomial of the fourth kind `W_n` as a polynomial object with double coefficients.

  See [[eval-chebyshev-W]] for the definition and [[chebyshev-W-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-W]] for values at a high degree.

  See also [[chebyshev-W-ratio]], [[eval-chebyshev-W]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-W-ratio degree))))

;;

(declare gegenbauer-C-ratio jacobi-P-ratio)

;; Exact arithmetic helpers. Ratios and big integers are boxed numbers, so `clojure.core` arithmetic is used
;; on them (the `fastmath.core` functions would convert them to doubles) and the boxed math warning is off;
;; longs and doubles still use `fastmath.core`.

(set! *unchecked-math* true)

(defn- exact-rational
  "The exact rational value of a finite double: its binary value, not the decimal number it prints as."
  [^double x]
  (rationalize (BigDecimal. x)))

(defn- rational->double
  "Converts a rational number to the nearest double (the quotient is formed with 40 significant digits)."
  ^double [r]
  (if (ratio? r)
    (.doubleValue (.divide (BigDecimal. ^BigInteger (biginteger (numerator r)))
                           (BigDecimal. ^BigInteger (biginteger (denominator r)))
                           (MathContext. 40)))
    (double r)))

(defn- value-at-infinity
  "The value at an infinite `x` of the polynomial with the exact ascending coefficients `cfs`: the limit
  of its leading non-zero term, which is infinite with the sign of that term. A constant polynomial gives
  its constant, the zero polynomial gives 0."
  ^double [cfs ^double x]
  (loop [k (m/long-dec (count cfs))]
    (cond
      (m/neg? k) 0.0
      (zero? (cfs k)) (recur (m/long-dec k))
      (m/zero? k) (rational->double (cfs 0))
      :else (if (= (pos? (cfs k)) (or (m/pos? x) (m/even? k))) ##Inf ##-Inf))))

(set! *unchecked-math* :warn-on-boxed)

(defn eval-legendre-P
  "Evaluates the Legendre polynomial `P_n` at `x`.

  The polynomials are orthogonal on `[-1, 1]` with the weight 1 and satisfy `P_0 = 1`, `P_1 = x` and `n*P_n = (2n-1)*x*P_(n-1) - (n-1)*P_(n-2)`. They are polynomials in `x` for every real `x`, so outside `[-1, 1]` the value is the continuation of the polynomial. `P_n(1) = 1` and `P_n(-x) = (-1)^n*P_n(x)`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree. An infinite `x` gives `x` for degree 1 and for a higher degree the infinity of the leading term: `##Inf` for an even degree, `x` for an odd one. The error is a few units of roundoff times `n` (absolute, for values of at most `n+1`) inside `[-1, 1]`, and relative outside, so this form stays accurate for any degree, unlike evaluating the coefficients from [[legendre-P]].

  See also [[legendre-P]], [[legendre-P-ratio]], [[eval-gegenbauer-C]] (`P_n` is the Gegenbauer polynomial of order 0.5), [[eval-jacobi-P]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (case (int degree)
    0 1.0
    1 x
    (if (m/inf? x)
      (if (m/even? degree) ##Inf x)
      (loop [i (long 2)
             pprev 1.0
             prev x]
        (if (m/> i degree)
          prev
          (recur (m/inc i) prev
                 (m// (m/- (m/* (m/dec (m/* 2.0 i)) x prev)
                           (m/* (m/dec i) pprev)) i)))))))

(defn legendre-P-ratio
  "Creates the Legendre polynomial `P_n` with exact rational coefficients.

  `P_0 = 1`, `P_1 = x` and `n*P_n = (2n-1)*x*P_(n-1) - (n-1)*P_(n-2)`; see [[eval-legendre-P]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients are fractions with a power of two as the denominator, vanish for the powers of the wrong parity and alternate in sign.

  See also [[legendre-P]] (double coefficients), [[eval-legendre-P]] (direct evaluation), [[gegenbauer-C-ratio]]."
  [^long degree]
  (check-degree! degree)
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [0 1])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [0 1])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (scale (sub (mult prev (ratio-polynomial [0 (m/dec (m/* 2.0 i))]))
                           (scale pprev (m/dec i))) (/ 1 i)))))))

(defn legendre-P
  "Creates the Legendre polynomial `P_n` as a polynomial object with double coefficients.

  See [[eval-legendre-P]] for the definition and [[legendre-P-ratio]] for the exact rational coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-legendre-P]] for values at a high degree.

  See also [[legendre-P-ratio]], [[eval-legendre-P]]."
  [^long degree]
  (polynomial (coeffs (legendre-P-ratio degree))))

;;

(defn eval-gegenbauer-C
  "Evaluates the Gegenbauer (ultraspherical) polynomial `C_n^(a)` of order `a` at `x`.

  The polynomials have the generating function `(1 - 2*x*t + t^2)^(-a)` and satisfy `C_0 = 1`, `C_1 = 2*a*x` and `n*C_n = 2*(n+a-1)*x*C_(n-1) - (n+2*a-2)*C_(n-2)`. For `a > -1/2` and `a` not 0 they are orthogonal on `[-1, 1]` with the weight `(1 - x^2)^(a - 1/2)`. They are polynomials in `x` for every real `x` and every real order, and `C_n(-x) = (-1)^n*C_n(x)`. The order 1 gives the Chebyshev polynomials of the second kind and the order 0.5 the Legendre polynomials; these two orders return exactly the values of [[eval-chebyshev-U]] and [[eval-legendre-P]].

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `order` (double): the order `a`, any real number, default 1.0 (the two-argument form). A negative order is allowed. For the order 0 the recurrence gives `C_0 = 1` and `C_n = 0` for every higher degree (the same values as `scipy` and `mpmath` give); the Chebyshev polynomials of the first kind are the limit of `C_n/a`, not the value.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `order` and `x` are. `##NaN` as `x` or as `order` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading non-zero term of the polynomial (so 0 for the order 0, and `##NaN` for an infinite order). That case uses the exact coefficients and is slower for a high degree (about a second for degree 400). The error is a few units of roundoff times `n` (absolute, for values of at most `n+1`) inside `[-1, 1]`, and relative outside; unlike evaluating the coefficients from [[gegenbauer-C]] it stays accurate for any degree. The value is continuous in the order, also at 1 and 0.5.

  See also [[gegenbauer-C]], [[gegenbauer-C-ratio]], [[eval-chebyshev-U]], [[eval-legendre-P]], [[eval-jacobi-P]]."
  (^double [^long degree ^double x] (eval-gegenbauer-C degree 1.0 x))
  (^double [^long degree ^double order ^double x]
   (check-degree! degree)
   (cond
     (m/== order 1.0) (eval-chebyshev-U degree x)
     (m/== order 0.5) (eval-legendre-P degree x)
     (m/zero? degree) 1.0
     (m/inf? x) (if (m/invalid-double? order)
                  ##NaN
                  (value-at-infinity (coeffs (gegenbauer-C-ratio degree order)) x))
     (m/== degree 1) (m/* 2.0 order x)
     :else (let [o2 (m/* 2.0 order)]
             (loop [i (long 2)
                    pprev 1.0
                    prev (m/* 2.0 order x)]
               (if (m/> i degree)
                 prev
                 (recur (m/inc i) prev
                        (m// (m/- (m/* 2.0 (m/dec (m/+ order i)) x prev)
                                  (m/* (m/+ i o2 -2.0) pprev)) i))))))))

(set! *unchecked-math* true)

(defn gegenbauer-C-ratio
  "Creates the Gegenbauer (ultraspherical) polynomial `C_n^(a)` with exact rational coefficients.

  `C_0 = 1`, `C_1 = 2*a*x` and `n*C_n = 2*(n+a-1)*x*C_(n-1) - (n+2*a-2)*C_(n-2)`; see [[eval-gegenbauer-C]] for the definition, the orthogonality and the special orders.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `order` (double): the order `a`, any real number. It is converted with `rationalize`, so a double becomes the exact decimal number it prints as (`0.3` gives `3/10`) and the coefficients are exact for that number. The orders 1 and 0.5 return [[chebyshev-U-ratio]] and [[legendre-P-ratio]]. For the order 0 every polynomial of a positive degree is the zero polynomial (zero coefficients).

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. A NaN or infinite `order` throws an `IllegalArgumentException` for a positive degree; degree 0 returns the constant 1 for any `order`.

  See also [[gegenbauer-C]] (double coefficients), [[eval-gegenbauer-C]] (direct evaluation)."
  [^long degree ^double order]
  (check-degree! degree)
  (cond
    (m/== order 1.0) (chebyshev-U-ratio degree)
    (m/== order 0.5) (legendre-P-ratio degree)
    (m/zero? degree) RONE
    :else (let [alpha (rationalize order)
                alpha2 (*' 2 alpha)]
            (loop [i (long 2)
                   pprev RONE
                   prev (ratio-polynomial [0 alpha2])]
              (if (m/> i degree)
                prev
                (recur (m/inc i) prev
                       (scale (sub (mult prev (ratio-polynomial [0 (*' 2 (+' alpha (m/long-dec i)))]))
                                   (scale pprev (+' i alpha2 -2)))
                              (/ 1 i))))))))

(set! *unchecked-math* :warn-on-boxed)

(defn gegenbauer-C
  "Creates the Gegenbauer (ultraspherical) polynomial `C_n^(a)` as a polynomial object with double coefficients.

  See [[eval-gegenbauer-C]] for the definition and [[gegenbauer-C-ratio]] for the exact rational coefficients, which are converted to doubles.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `order` (double): the order `a`, any real number, default 1.0 (the Chebyshev polynomials of the second kind). A NaN or infinite order throws an `IllegalArgumentException` for a positive degree; see [[gegenbauer-C-ratio]] for the conversion of the order.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-gegenbauer-C]] for values at a high degree.

  See also [[gegenbauer-C-ratio]], [[eval-gegenbauer-C]], [[chebyshev-U]], [[legendre-P]]."
  ([^long degree] (gegenbauer-C degree 1.0))
  ([^long degree ^double order] (polynomial (coeffs (gegenbauer-C-ratio degree order)))))

;;

(defn eval-hermite-H
  "Hermite polynomials"
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 (m/* 2.0 x)
    (loop [i (long 2)
           pprev 1.0
           prev (m/* 2.0 x)]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (m/* 2.0 (m/- (m/* x prev)
                             (m/* (m/dec i) pprev))))))))

(defn hermite-H-ratio
  "Hermite polynomials"
  [^long degree]
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [0 2])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [0 2])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (scale (sub (mult prev (ratio-polynomial [0 1]))
                           (scale pprev (m/dec i))) 2))))))

(defn hermite-H
  [^long degree]
  (polynomial (coeffs (hermite-H-ratio degree))))

(defn eval-hermite-He
  "Hermite polynomials"
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 x
    (loop [i (long 2)
           pprev 1.0
           prev x]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (m/- (m/* x prev)
                    (m/* (m/dec i) pprev)))))))

(defn hermite-He-ratio
  "Hermite polynomials"
  [^long degree]
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [0 1])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [0 1])]
      (if (m/> i degree)
        prev
        (recur (m/inc i) prev
               (sub (mult prev (ratio-polynomial [0 1]))
                    (scale pprev (m/dec i))))))))

(defn hermite-He
  [^long degree]
  (polynomial (coeffs (hermite-He-ratio degree))))


;;

;; Jacobi polynomials with the three term recurrence, which divides by
;; 2 i (i + s) (2 i + s - 2) at step i, s = alpha + beta. For s a negative integer the divisor is zero for
;; some i <= degree and the recurrence is undefined; close to such an s it loses digits. In both cases the
;; explicit sum  sum_k C(n+alpha, n-k) C(n+beta, k) ((x-1)/2)^k ((x+1)/2)^(n-k)  in exact arithmetic is used.

(def ^:private ^:const jacobi-recurrence-min-gap 0.1)

(defn- jacobi-recurrence-gap
  "The smallest distance to zero of the two factors `i + s` and `2 i + s - 2` of the recurrence divisor,
  over the steps `i` from 2 to `degree`; infinite when there are no steps."
  ^double [^long degree ^double s]
  (loop [i (long 2)
         gap ##Inf]
    (if (m/> i degree)
      gap
      (recur (m/inc i)
             (m/min gap (m/abs (m/+ i s)) (m/abs (m/+ (m/* 2.0 i) s -2.0)))))))

(set! *unchecked-math* true)

(defn- generalized-binomials
  "The exact binomial coefficients `C(top, k)` for `k` from 0 to `n`, for a rational `top`."
  [top ^long n]
  (loop [k (long 0)
         c 1N
         cs (transient [])]
    (if (m/> k n)
      (persistent! cs)
      (recur (m/inc k) (/ (*' c (-' top k)) (m/inc k)) (conj! cs c)))))

(defn- powers-up-to
  "The vector of `base^k` for `k` from 0 to `n`, with the multiplication `multiply` and the unit `one`."
  [base ^long n multiply one]
  (vec (take (m/long-inc n) (iterate #(multiply % base) one))))

(defn- jacobi-explicit-value
  "The value of the Jacobi polynomial at the exact binary values of the doubles, as the nearest double."
  ^double [^long degree ^double alpha ^double beta ^double x]
  (let [rx (exact-rational x)
        lowers (powers-up-to (/ (-' rx 1) 2) degree *' 1)
        uppers (powers-up-to (/ (+' rx 1) 2) degree *' 1)
        from-alpha (generalized-binomials (+' degree (exact-rational alpha)) degree)
        from-beta (generalized-binomials (+' degree (exact-rational beta)) degree)]
    (rational->double
     (reduce +' (map (fn [^long k]
                       (let [rest-degree (m/long-sub degree k)]
                         (*' (from-alpha rest-degree) (from-beta k) (lowers k) (uppers rest-degree))))
                     (range (m/inc degree)))))))

(defn- jacobi-explicit-ratio
  "The Jacobi polynomial of the rational parameters `alpha`, `beta` as an exact `PolynomialR`."
  [^long degree alpha beta]
  (let [lowers (powers-up-to (ratio-polynomial [-1/2 1/2]) degree mult RONE)
        uppers (powers-up-to (ratio-polynomial [1/2 1/2]) degree mult RONE)
        from-alpha (generalized-binomials (+' degree alpha) degree)
        from-beta (generalized-binomials (+' degree beta) degree)]
    (reduce add (map (fn [^long k]
                       (let [rest-degree (m/long-sub degree k)]
                         (scale (mult (lowers k) (uppers rest-degree))
                                (*' (from-alpha rest-degree) (from-beta k)))))
                     (range (m/inc degree))))))

(defn- jacobi-degenerate?
  "True when the recurrence divisor is zero for a step up to `degree`, for rational `alpha` and `beta`:
  `alpha + beta` is an integer from `-degree` to `-2`, or an even integer from `2 - 2 degree` to `-2`."
  [^long degree alpha beta]
  (let [s (+' alpha beta)]
    (and (integer? s)
         (<= s -2)
         (or (<= (m/long-sub 0 degree) s)
             (and (even? s) (<= (m/long-sub 2 (m/long-mult 2 degree)) s))))))

(set! *unchecked-math* :warn-on-boxed)

(defn eval-jacobi-P
  "Evaluates the Jacobi polynomial `P_n^(a,b)` at `x`.

  The polynomials satisfy `P_0 = 1`, `P_1 = (a+1) + (a+b+2)*(x-1)/2` and a three term recurrence in `n`. They are equal to the sum over `k` of `C(n+a, n-k)*C(n+b, k)*((x-1)/2)^k*((x+1)/2)^(n-k)`, with generalized binomial coefficients, which defines them for every real `a` and `b`. For `a > -1` and `b > -1` they are orthogonal on `[-1, 1]` with the weight `(1-x)^a*(1+x)^b`. `P_n(1) = C(n+a, n)` and `P_n^(a,b)(-x) = (-1)^n*P_n^(b,a)(x)`. For `a = b = 0` they are the Legendre polynomials, and for `a = b = c - 1/2` they are proportional to the Gegenbauer polynomials of order `c`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `alpha` (double): the parameter `a`, any real number.
  - `beta` (double): the parameter `b`, any real number.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever the other arguments are. `##NaN` as `x`, `alpha` or `beta` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading non-zero term of the polynomial (a finite constant when the polynomial is constant, as for degree 1 and `a + b = -2`); that case uses the exact coefficients and is slower for a high degree.

  The three term recurrence divides by zero when `a + b` is a negative integer from `-degree` to `-2` (or an even one down to `2 - 2*degree`) and loses digits when `a + b` is within 0.1 of such a value. In both cases the value is computed in exact arithmetic from the sum above, at the exact binary values of the arguments, and rounded once; this takes time growing quickly with the degree (about a second for degree 200). Otherwise the error is a few units of roundoff times `n` (absolute, for values of at most `n+1`) inside `[-1, 1]`, measured at most 1.1 on the reference grid. At a multiple root at 1 or -1 (a negative integer `alpha` or `beta` of at most the degree in absolute value) the error is not relative to the tiny value next to the root but to the size of the polynomial elsewhere.

  See also [[jacobi-P]], [[jacobi-P-ratio]], [[eval-gegenbauer-C]], [[eval-legendre-P]]."
  ^double [^long degree ^double alpha ^double beta ^double x]
  (check-degree! degree)
  (cond
    (m/zero? degree) 1.0
    (m/inf? x) (if (or (m/invalid-double? alpha) (m/invalid-double? beta))
                 ##NaN
                 (value-at-infinity (coeffs (jacobi-P-ratio degree alpha beta)) x))
    (m/== degree 1) (m/+ (m/inc alpha) (m/* 0.5 (m/+ alpha beta 2.0) (m/dec x)))
    (and (m/valid-double? x)
         (m/< (jacobi-recurrence-gap degree (m/+ alpha beta)) jacobi-recurrence-min-gap))
    (jacobi-explicit-value degree alpha beta x)
    :else (loop [i (long 2)
                 pprev 1.0
                 prev (m/+ (m/inc alpha) (m/* 0.5 (m/+ alpha beta 2.0) (m/dec x)))]
            (if (m/> i degree)
              prev
              (let [a (m/+ i alpha)
                    b (m/+ i beta)
                    c (m/+ a b)]
                (recur (m/inc i) prev
                       (m// (m/- (m/* (m/dec c) (m/+ (m/* c (m/- c 2.0) x)
                                                     (m/* (m/- a b) (m/- c (m/* 2.0 i)))) prev)
                                 (m/* 2.0 (m/dec a) (m/dec b) c pprev))
                            (m/* 2.0 i (m/- c i) (m/- c 2.0)))))))))

(set! *unchecked-math* true)

(defn jacobi-P-ratio
  "Creates the Jacobi polynomial `P_n^(a,b)` with exact rational coefficients.

  See [[eval-jacobi-P]] for the definition, the orthogonality and the parameters.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `alpha`, `beta` (doubles): the parameters `a` and `b`, any real numbers. They are converted with `rationalize`, so a double becomes the exact decimal number it prints as (`0.1` gives `1/10`) and the coefficients are exact for those numbers.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. When `a + b` is a negative integer from `-n` to `-2` (or an even one down to `2 - 2*n`), where the three term recurrence is undefined, the explicit sum is used; the time then grows with the cube of the degree. A NaN or infinite parameter throws an `IllegalArgumentException` for a positive degree; degree 0 returns the constant 1 for any parameters.

  See also [[jacobi-P]] (double coefficients), [[eval-jacobi-P]] (direct evaluation), [[gegenbauer-C-ratio]]."
  [^long degree ^double alpha ^double beta]
  (check-degree! degree)
  (if (m/zero? degree)
    RONE
    (let [alpha (rationalize alpha)
          beta (rationalize beta)]
      (if (jacobi-degenerate? degree alpha beta)
        (jacobi-explicit-ratio degree alpha beta)
        (loop [i (long 2)
               pprev RONE
               prev (let [ab22 (/ (+' alpha beta 2) 2)]
                      (ratio-polynomial [(-' (inc' alpha) ab22) ab22]))]
          (if (m/> i degree)
            prev
            (let [a (+' i alpha)
                  b (+' i beta)
                  c (+' a b)]
              (recur (m/inc i) prev
                     (scale (sub (scale (mult prev (ratio-polynomial [(*' (-' a b) (-' c (m/long-mult 2 i)))
                                                                      (*' c (-' c 2))])) (dec' c))
                                 (scale pprev (*' 2 (dec' a) (dec' b) c)))
                            (/ 1 (*' 2 i (-' c i) (-' c 2))))))))))))

(set! *unchecked-math* :warn-on-boxed)

(defn jacobi-P
  "Creates the Jacobi polynomial `P_n^(a,b)` as a polynomial object with double coefficients.

  See [[eval-jacobi-P]] for the definition and [[jacobi-P-ratio]] for the exact rational coefficients, which are converted to doubles.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated.
  - `alpha`, `beta` (doubles): the parameters `a` and `b`, any real numbers. A NaN or infinite parameter throws an `IllegalArgumentException` for a positive degree; see [[jacobi-P-ratio]] for the conversion.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-jacobi-P]] for values at a high degree.

  See also [[jacobi-P-ratio]], [[eval-jacobi-P]], [[gegenbauer-C]], [[legendre-P]]."
  [^long degree ^double alpha ^double beta]
  (polynomial (coeffs (jacobi-P-ratio degree alpha beta))))

;;

(defn eval-bessel-y
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 (m/inc x)
    (loop [i (long 2)
           pprev 1.0
           prev (m/inc x)]
      (if (> i degree)
        prev
        (recur (inc i) pprev
               (m/+ (m/* (m/dec (m/* 2 i)) x prev)
                    pprev))))))


(defn bessel-y-ratio
  [^long degree]
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [1 1])
    (loop [i (long 2)
           pprev RONE
           prev (ratio-polynomial [1 1])]
      (if (> i degree)
        prev
        (recur (inc i) prev
               (add (mult prev (ratio-polynomial [0 (m/dec (m/* 2 i))])) pprev))))))

(defn bessel-y
  [^long degree]
  (polynomial (coeffs (bessel-y-ratio degree))))

(defn eval-bessel-t
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 (m/inc x)
    (loop [i (long 2)
           pprev 1.0
           prev (m/inc x)]
      (if (> i degree)
        prev
        (recur (inc i) prev
               (m/+ (m/* (m/dec (m/* 2 i)) prev)
                    (m/* x x pprev)))))))

(defn bessel-t-ratio
  [^long degree]
  (case (int degree)
    0 RONE
    1 (ratio-polynomial [1 1])
    (let [p001 (ratio-polynomial [0 0 1])]
      (loop [i (long 2)
             pprev RONE
             prev (ratio-polynomial [1 1])]
        (if (> i degree)
          prev
          (recur (inc i) prev
                 (add (scale prev (m/dec (m/* 2 i)))
                      (mult pprev p001))))))))

(defn bessel-t
  [^long degree]
  (polynomial (coeffs (bessel-t-ratio degree))))

;;

(defn eval-meixner-pollaczek-P
  ^double [^long degree ^double lambda ^double phi ^double x]
  (case (int degree)
    0 1.0
    1 (m/* 2.0 (m/+ (m/* lambda (m/cos phi))
                    (m/* x (m/sin phi))))
    (let [cp (m/cos phi)
          sp (m/sin phi)
          l2 (m/* 2.0 lambda)]
      (loop [i (long 2)
             pprev 1.0
             prev (m/* 2.0 (m/+ (m/* lambda cp)
                                (m/* x sp)))]
        (if (> i degree)
          prev
          (recur (inc i) prev
                 (m// (m/- (m/* 2.0 (m/+ (m/* x sp)
                                         (m/* (m/dec (m/+ i lambda)) cp)) prev)
                           (m/* (m/+ i l2 -2) pprev)) i)))))))

(defn meixner-pollaczek-P-ratio
  [^long degree ^double lambda ^double phi]
  (case (int degree)
    0 RONE
    1 (scale (ratio-polynomial [(m/* lambda (m/cos phi))
                                (m/sin phi)]) 2)
    (let [cp (m/cos phi)
          sp (m/sin phi)
          l2 (m/* 2.0 lambda)]
      (loop [i (long 2)
             pprev RONE
             prev (scale (ratio-polynomial [(m/* lambda (m/cos phi))
                                            (m/sin phi)]) 2)]
        (if (> i degree)
          prev
          (recur (inc i) prev
                 (scale (sub (scale (mult prev (ratio-polynomial [(m/* (m/dec (m/+ i lambda)) cp)
                                                                  sp])) 2)
                             (scale pprev (m/+ i l2 -2))) (/ 1 i))))))))

(defn meixner-pollaczek-P
  [^long degree ^double lambda ^double phi]
  (polynomial (coeffs (meixner-pollaczek-P-ratio degree lambda phi))))

;; Ince polynomials

;; https://www.mathworks.com/matlabcentral/fileexchange/44932-ince-polynomials
;; https://dlmf.nist.gov/28.31
;; Miguel A. Bandres and Julio C. Gutierrez-Vega Ince–Gaussian modes of the paraxial wave equation and stable resonators

(defn- ince-ev
  ^doubles [^RealMatrix m ^long order]
  (let [ed (EigenDecomposition. m)
        ;; be sure the order of eigenvalues is increasing
        ro (long (nth (m/order (.getRealEigenvalues ed)) order))]
    (-> (.getEigenvector ed ro)
        (v/vec->array))))

(defn- ince-pmv-gamma
  [^long start ^long p]
  (map (fn [^long mv]
         (m/exp (m/* 0.5 (m/+ (Gamma/logGamma (m/inc (m/* 0.5 (m/+ p mv))))
                              (Gamma/logGamma (m/inc (m/* 0.5 (m/- p mv)))))))) (range start (m/inc p) 2)))

(defn- ince-c-coeffs-even
  ([^long p ^long m ^double e normalization]
   (let [n (m// p 2)
         N (m/inc n)
         order (m/long-div m 2)
         ^RealMatrix mat (MatrixUtils/createRealMatrix N N)]
     (doseq [^long i (range 1 N)]
       (.setEntry mat i i (m/* 4.0 i i)))
     (doseq [^long i (range 0 n)
             :let [i+ (m/inc i)]]
       (.setEntry mat i i+ (m/* e (m/+ n i+)))
       (.setEntry mat i+ i (m/* e (m/- n i))))
     (.setEntry mat 1 0 (m/* 2.0 (.getEntry mat 1 0)))
     (let [^doubles a (ince-ev mat order)
           sgn (m/sgn (v/sum a))
           norm (case normalization
                  :trigonometric (let [^doubles a2 (v/sq a)]
                                   (m/sqrt (m/+ (Array/aget a2 0) (v/sum a2))))
                  :millers (m/sqrt (m/+ (m/* 2.0 (m/sq (Array/aget a 0)) (m/sq (Gamma/gamma N)))
                                        (v/sum (map-indexed (fn [^long id ^double gs]
                                                              (m/sq (m/* gs (Array/aget a (m/inc id)))))
                                                            (ince-pmv-gamma 2 p)))))
                  1.0)]
       (v/mult (v/div a norm) sgn)))))

(defn- ince-s-coeffs-even
  [^long p ^long m ^double e normalization]
  (let [n (m// p 2)
        order (m/long-dec (m/long-div m 2))
        ^RealMatrix mat (MatrixUtils/createRealMatrix n n)]
    (doseq [^long i (range 0 n)]
      (.setEntry mat i i (m/* 4.0 (m/sq (m/inc i)))))
    (doseq [^long i (range 0 (m/dec n))
            :let [i+ (m/inc i)]]
      (.setEntry mat i i+ (m/* e (m/+ n i 2.0)))
      (.setEntry mat i+ i (m/* e (m/- n i+))))
    (let [scaler (m/seq->double-array (range 1 (m/inc n)))
          ^doubles a (ince-ev mat order)
          sgn (m/sgn (v/sum (v/emult a scaler)))
          norm (case normalization
                 :trigonometric (m/sqrt (v/sum (v/sq a)))
                 :millers (m/sqrt (v/sum (map (fn [^double v ^double gs]
                                                (m/sq (m/* v gs))) a (ince-pmv-gamma 2 p))))
                 1.0)]
      (v/mult (v/div a norm) sgn))))

(defn- ince-c-coeffs-odd
  [^long p ^long m ^double e normalization]
  (let [n (m// (m/dec p) 2)
        N (m/inc n)
        order (m/long-div (m/long-dec m) 2)
        ^RealMatrix mat (MatrixUtils/createRealMatrix N N)
        he (m/* 0.5 e)]
    (doseq [^long i (range 1 N)]
      (.setEntry mat i i (m/sq (m/inc (m/* 2.0 i)))))
    (.setEntry mat 0 0 (m/+ he (m/* he p) 1.0))
    (doseq [^long i (range 0 n)
            :let [i+ (m/inc i)]]
      (.setEntry mat i i+ (m/* he (m/+ p (m/* 2.0 i) 3.0)))
      (.setEntry mat i+ i (m/* he (m/- p (m/* 2.0 i) 1.0))))
    (let [^doubles a (ince-ev mat order)
          sgn (m/sgn (v/sum a))
          norm (case normalization
                 :trigonometric (m/sqrt (v/sum (v/sq a)))
                 :millers (m/sqrt (v/sum (map (fn [^double v ^double gs]
                                                (m/sq (m/* v gs))) a (ince-pmv-gamma 1 p))))
                 1.0)]
      (v/mult (v/div a norm) sgn))))

(defn- ince-s-coeffs-odd
  [^long p ^long m ^double e normalization]
  (let [n (m// (m/dec p) 2)
        N (m/inc n)
        order (m/long-div (m/long-dec m) 2)
        ^RealMatrix mat (MatrixUtils/createRealMatrix N N)
        he (m/* 0.5 e)]
    (doseq [^long i (range 1 N)]
      (.setEntry mat i i (m/sq (m/inc (m/* 2.0 i)))))
    (.setEntry mat 0 0 (m/- 1.0 he (m/* he p)))
    (doseq [^long i (range 0 n)
            :let [i+ (m/inc i)]]
      (.setEntry mat i i+ (m/* he (m/+ p (m/* 2.0 i) 3.0)))
      (.setEntry mat i+ i (m/* he (m/- p (m/* 2.0 i) 1.0))))
    (let [scaler (m/seq->double-array (map (fn [^long v] (m/inc (m/* 2.0 v))) (range N)))
          ^doubles a (ince-ev mat order)
          sgn (m/sgn (v/sum (v/emult a scaler)))
          norm (case normalization
                 :trigonometric (m/sqrt (v/sum (v/sq a)))
                 :millers (m/sqrt (v/sum (map (fn [^double v ^double gs]
                                                (m/sq (m/* v gs))) a (ince-pmv-gamma 1 p))))
                 1.0)]
      (v/mult (v/div a norm) sgn))))

(defn ince-C-coeffs
  [^long p ^long m ^double e normalization]
  (if (m/even? p)
    (ince-c-coeffs-even p m e normalization)
    (ince-c-coeffs-odd p m e normalization)))

(defn ince-S-coeffs
  [^long p ^long m ^double e normalization]
  (if (m/even? p)
    (ince-s-coeffs-even p m e normalization)
    (ince-s-coeffs-odd p m e normalization)))

(defmacro ^:private ince-loop
  [coeffs s x f]
  `(loop [r# (long 0)
          sum# (double 0.0)]
     (if (m/== r# ~s)
       sum#
       (recur (m/inc r#) (m/+ sum# (m/* (Array/aget ~coeffs r#) (~f r# ~x)))))))

(defn- ince-cos-even ^double [^long r ^double x] (m/cos (m/* 2.0 r x)))
(defn- ince-cos-odd ^double [^long r ^double x] (m/cos (m/* (m/inc (m/* 2.0 r)) x)))

(defn ince-C
  "Ince C polynomial of order p and degree m.

  `normalization` parameter can be `:none` (default), `:trigonometric` or `millers`."
  ([^long p ^long m ^double e] (ince-C p m e :none))
  ([^long p ^long m ^double e normalization]
   (assert (m/even? (m/long-sub p m)) "p and m must be the same parity!")
   (let [^doubles coeffs (ince-C-coeffs p m e normalization)
         s (alength coeffs)]
     (if (m/even? p)
       (fn ^double [^double x] (ince-loop coeffs s x ince-cos-even))
       (fn ^double [^double x] (ince-loop coeffs s x ince-cos-odd))))))

(defn- ince-sin-even ^double [^long r ^double x] (m/sin (m/* 2.0 (m/inc r) x)))
(defn- ince-sin-odd ^double [^long r ^double x] (m/sin (m/* (m/inc (m/* 2.0 r)) x)))

(defn ince-S
  "Ince S polynomial of order p and degree m.

  `normalization` parameter can be `:none` (default), `:trigonometric` or `millers`."
  ([^long p ^long m ^double e] (ince-S p m e :none))
  ([^long p ^long m ^double e normalization]
   (assert (m/even? (m/long-sub p m)) "p and m must be the same parity!")
   (let [^doubles coeffs (ince-S-coeffs p m e normalization)
         s (alength coeffs)]
     (if (m/even? p)
       (fn ^double [^double x] (ince-loop coeffs s x ince-sin-even))
       (fn ^double [^double x] (ince-loop coeffs s x ince-sin-odd))))))

(defn- ince-cosh-even ^double [^long r ^double x] (m/cosh (m/* 2.0 r x)))
(defn- ince-cosh-odd ^double [^long r ^double x] (m/cosh (m/* (m/inc (m/* 2.0 r)) x)))

(defn ince-C-radial
  "Ince C polynomial of order p and degree m.

  `normalization` parameter can be `:none` (default), `:trigonometric` or `millers`."
  ([^long p ^long m ^double e] (ince-C-radial p m e :none))
  ([^long p ^long m ^double e normalization]
   (assert (m/even? (m/long-sub p m)) "p and m must be the same parity!")
   (let [^doubles coeffs (ince-C-coeffs p m e normalization)
         s (alength coeffs)]
     (if (m/even? p)
       (fn ^double [^double x] (ince-loop coeffs s x ince-cosh-even))
       (fn ^double [^double x] (ince-loop coeffs s x ince-cosh-odd))))))

(defn- ince-sinh-even ^double [^long r ^double x] (m/sinh (m/* 2.0 (m/inc r) x)))
(defn- ince-sinh-odd ^double [^long r ^double x] (m/sinh (m/* (m/inc (m/* 2.0 r)) x)))

(defn ince-S-radial
  "Ince S polynomial of order p and degree m.

  `normalization` parameter can be `:none` (default), `:trigonometric` or `millers`."
  ([^long p ^long m ^double e] (ince-S-radial p m e :none))
  ([^long p ^long m ^double e normalization]
   (assert (m/even? (m/long-sub p m)) "p and m must be the same parity!")
   (let [^doubles coeffs (ince-S-coeffs p m e normalization)
         s (alength coeffs)]
     (if (m/even? p)
       (fn ^double [^double x] (ince-loop coeffs s x ince-sinh-even))
       (fn ^double [^double x] (ince-loop coeffs s x ince-sinh-odd))))))
