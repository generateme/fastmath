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
  "Evaluates the generalized Laguerre polynomial `L_n^(a)` of order `a` at `x`.

  The polynomials are `L_n^(a)(x) = sum over k of (-1)^k*C(n+a, n-k)*x^k/k!` with generalized binomial coefficients, and satisfy `L_0 = 1`, `L_1 = 1 + a - x` and `n*L_n = (2n-1+a-x)*L_(n-1) - (n-1+a)*L_(n-2)`. For `a > -1` they are orthogonal on `[0, Inf)` with the weight `x^a*exp(-x)`. `L_n(0) = C(n+a, n)`. The order 0 gives the Laguerre polynomials (the default of the two-argument form).

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `order` (double): the order `a`, any real number, default 0.0 (the two-argument form). A negative order is allowed: the polynomials are then defined by the same sum and recurrence, but they are not orthogonal.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `order` and `x` are. `##NaN` as `x` or as `order` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading term `(-1)^n*x^n/n!`: `##-Inf` for `##Inf` and an odd degree, `##Inf` otherwise, and `##NaN` for a NaN or infinite `order`. The error of the recurrence grows like the square of the degree: relative to the largest value of `L_0 ... L_n` at `x` it is about `(n+1)^2` times the roundoff (roughly 1e-13 at degree 100 close to `x = 0`, much smaller at a low degree), which is larger than for [[eval-chebyshev-T]]. For a negative order of at most -3 and a high degree the intermediate values are far larger than the result and the relative error of the result can be much larger. Unlike evaluating the coefficients from [[laguerre-L]] it does not suffer from the cancellation of the large alternating coefficients.

  See also [[laguerre-L]], [[laguerre-L-ratio]], [[eval-hermite-H]]."
  (^double [^long degree ^double x] (eval-laguerre-L degree 0.0 x))
  (^double [^long degree ^double order ^double x]
   (check-degree! degree)
   (cond
     (m/zero? degree) 1.0
     ;; the leading term is (-1)^n x^n / n!
     (m/inf? x) (if (m/invalid-double? order)
                  ##NaN
                  (if (and (m/pos? x) (m/odd? degree)) ##-Inf ##Inf))
     (m/== degree 1) (m/- (m/inc order) x)
     :else (loop [i (long 2)
                  pprev 1.0
                  prev (m/- (m/inc order) x)]
             (if (m/> i degree)
               prev
               (recur (m/inc i) prev
                      (m// (m/- (m/* (m/+ order (m/- (m/* 2.0 i) 1.0 x)) prev)
                                (m/* (m/+ order (m/dec i)) pprev)) i)))))))

(set! *unchecked-math* true)

(defn laguerre-L-ratio
  "Creates the generalized Laguerre polynomial `L_n^(a)` with exact rational coefficients.

  `L_0 = 1`, `L_1 = 1 + a - x` and `n*L_n = (2n-1+a-x)*L_(n-1) - (n-1+a)*L_(n-2)`; see [[eval-laguerre-L]] for the definition and the orthogonality.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `order` (double): the order `a`, any real number. It is converted with `rationalize`, so a double becomes the exact decimal number it prints as (`0.3` gives `3/10`) and the coefficients are exact for that number.

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients alternate in sign for an order above -1. A NaN or infinite `order` throws an `IllegalArgumentException` for a positive degree; degree 0 returns the constant 1 for any `order`.

  See also [[laguerre-L]] (double coefficients), [[eval-laguerre-L]] (direct evaluation), [[hermite-H-ratio]]."
  [^long degree ^double order]
  (check-degree! degree)
  (if (m/zero? degree)
    RONE
    (let [alpha (rationalize order)]
      (loop [i (long 2)
             pprev RONE
             prev (ratio-polynomial [(inc' alpha) -1])]
        (if (m/> i degree)
          prev
          (recur (m/inc i) prev
                 (scale (sub (mult prev (ratio-polynomial [(-' (+' alpha (m/long-mult 2 i)) 1) -1]))
                             (scale pprev (+' alpha (m/long-dec i))))
                        (/ 1 i))))))))

(set! *unchecked-math* :warn-on-boxed)

(defn laguerre-L
  "Creates the generalized Laguerre polynomial `L_n^(a)` as a polynomial object with double coefficients.

  See [[eval-laguerre-L]] for the definition and [[laguerre-L-ratio]] for the exact rational coefficients, which are converted to doubles.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `order` (double): the order `a`, any real number, default 0.0. A NaN or infinite order throws an `IllegalArgumentException` for a positive degree; see [[laguerre-L-ratio]] for the conversion of the order.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-laguerre-L]] for values at a high degree.

  See also [[laguerre-L-ratio]], [[eval-laguerre-L]], [[hermite-H]]."
  ([^long degree] (laguerre-L degree 0.0))
  ([^long degree ^double order]
   (polynomial (coeffs (laguerre-L-ratio degree order)))))

;;

(defn eval-chebyshev-T
  "Evaluates the Chebyshev polynomial of the first kind `T_n` at `x`.

  `T_n(cos(a)) = cos(n*a)`. The polynomials are orthogonal on `[-1, 1]` with the weight `1/sqrt(1 - x^2)` and satisfy `T_0 = 1`, `T_1 = x` and `T_(n+1) = 2*x*T_n - T_(n-1)`. They are polynomials in `x` for every real `x`: outside `[-1, 1]` the value is the continuation of the polynomial, `T_n(x) = cosh(n*acosh(x))` for `x > 1`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-U]] for values at a high degree.

  See also [[chebyshev-U-ratio]], [[eval-chebyshev-U]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-U-ratio degree))))

(defn eval-chebyshev-V
  "Evaluates the Chebyshev polynomial of the third kind `V_n` at `x`.

  `V_n(cos(a)) = cos((n+1/2)*a)/cos(a/2)`. The polynomials are orthogonal on `[-1, 1]` with the weight `sqrt((1+x)/(1-x))` and satisfy `V_0 = 1`, `V_1 = 2*x - 1` and `V_(n+1) = 2*x*V_n - V_(n-1)`. Also `V_n = U_n - U_(n-1)`, `V_n(1) = 1`, `V_n(-1) = (-1)^n*(2*n+1)` and `V_n(-x) = (-1)^n*W_n(x)`. They are polynomials in `x` for every real `x`, and the value is the polynomial for any argument outside `[-1, 1]` too.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-V]] for values at a high degree.

  See also [[chebyshev-V-ratio]], [[eval-chebyshev-V]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-V-ratio degree))))

(defn eval-chebyshev-W
  "Evaluates the Chebyshev polynomial of the fourth kind `W_n` at `x`.

  `W_n(cos(a)) = sin((n+1/2)*a)/sin(a/2)`. The polynomials are orthogonal on `[-1, 1]` with the weight `sqrt((1-x)/(1+x))` and satisfy `W_0 = 1`, `W_1 = 2*x + 1` and `W_(n+1) = 2*x*W_n - W_(n-1)`. Also `W_n = U_n + U_(n-1)`, `W_n(1) = 2*n+1`, `W_n(-1) = (-1)^n` and `W_n(-x) = (-1)^n*V_n(x)`. They are polynomials in `x` for every real `x`, and the value is the polynomial for any argument outside `[-1, 1]` too.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-chebyshev-W]] for values at a high degree.

  See also [[chebyshev-W-ratio]], [[eval-chebyshev-W]]."
  [^long degree]
  (polynomial (coeffs (chebyshev-W-ratio degree))))

;;

;; Exact and high precision arithmetic helpers. Ratios, big integers and big decimals are boxed numbers, so
;; `clojure.core` arithmetic is used on them (the `fastmath.core` functions would convert them to doubles) and
;; the boxed math warning is off; longs and doubles still use `fastmath.core`.

(set! *unchecked-math* true)

(defn- exact-rational
  "The exact rational value of a finite double: its binary value, not the decimal number it prints as."
  [^double x]
  (rationalize (BigDecimal. x)))

(defn- rational->double
  "Converts a rational number to the nearest double, a tie going to the even neighbour (also for subnormal
  results); an out of range magnitude gives an infinity or zero."
  ^double [r]
  (if (ratio? r)
    (let [^BigInteger n (biginteger (numerator r))
          ^BigInteger d (biginteger (denominator r))
          ^BigInteger a (.abs n)
          ;; scale the dividend so that the integer quotient `q` has 55 or 56 bits
          shift (- 55 (- (.bitLength a) (.bitLength d)))
          ^BigInteger scaled-a (if (neg? shift) a (.shiftLeft a (int shift)))
          ^BigInteger scaled-d (if (neg? shift) (.shiftLeft d (int (- shift))) d)
          ^"[Ljava.math.BigInteger;" qr (.divideAndRemainder scaled-a scaled-d)
          ^BigInteger q (aget qr 0)
          ;; one more low bit, set when the division left a remainder: the value is q+ * 2^(-shift-1), and a
          ;; set bit keeps the dropped part away from an exact tie
          ^BigInteger q+ (.or (.shiftLeft q 1) (if (zero? (.signum ^BigInteger (aget qr 1))) BigInteger/ZERO BigInteger/ONE))
          ;; the result is an integer of at most 53 bits times 2^quantum-exponent (subnormals: 2^-1074)
          quantum-exponent (max (- (.bitLength q+) shift 2 52) -1074)
          drop-bits (+ shift 1 quantum-exponent)
          ^BigInteger kept (.shiftRight q+ (int drop-bits))
          ^BigInteger dropped (.subtract q+ (.shiftLeft kept (int drop-bits)))
          order-to-half (.compareTo dropped (.shiftLeft BigInteger/ONE (int (dec drop-bits))))
          ^BigInteger rounded (if (or (pos? order-to-half) (and (zero? order-to-half) (.testBit kept 0)))
                                (.add kept BigInteger/ONE)
                                kept)
          magnitude (Math/scalb (.doubleValue rounded) (int quantum-exponent))]
      (if (neg? (.signum n)) (- magnitude) magnitude))
    (double r)))

(def ^:private ^:const adaptive-start-digits 64)
(def ^:private ^:const adaptive-max-digits 65536)

(defn- adaptive-decimal-value
  "Evaluates a sum in decimal arithmetic whose precision is raised until the result is settled.

  `sum-fn` takes a `MathContext` and returns the sum as a `BigDecimal`. The number of significant digits
  starts at 64 and is doubled until two successive results differ by at most 2^-60 of the later one (or by
  at most 1e-300), so the digits lost to cancellation are recovered; the result is the nearest double of the
  last sum. The precision stops growing at 65536 digits."
  ^double [sum-fn]
  (let [tolerance (BigDecimal. "8.0E-19")
        tiny (BigDecimal. "1E-300")]
    (loop [digits adaptive-start-digits
           ^BigDecimal previous (sum-fn (MathContext. (int adaptive-start-digits)))]
      (let [next-digits (* 2 digits)
            ^BigDecimal current (sum-fn (MathContext. (int next-digits)))
            ^BigDecimal difference (.abs (.subtract current previous))]
        (if (or (> next-digits adaptive-max-digits)
                (<= (.compareTo difference (.multiply (.abs current) tolerance)) 0)
                (<= (.compareTo difference tiny) 0))
          (.doubleValue current)
          (recur next-digits current))))))

(defn- decimal-powers
  "The vector of `base^k` for `k` from 0 to `n`, each rounded to `mc`."
  [^BigDecimal base ^long n ^MathContext mc]
  (loop [k (long 0)
         power BigDecimal/ONE
         powers (transient [])]
    (if (> k n)
      (persistent! powers)
      (recur (inc k) (.multiply power base mc) (conj! powers power)))))

(defn- decimal-binomials
  "The vector of the generalized binomial coefficients `C(top, k)` for `k` from 0 to `n`, rounded to `mc`."
  [^BigDecimal top ^long n ^MathContext mc]
  (loop [k (long 0)
         c BigDecimal/ONE
         cs (transient [])]
    (if (> k n)
      (persistent! cs)
      (recur (inc k)
             (.divide (.multiply c (.subtract top (BigDecimal. k) mc) mc) (BigDecimal. (inc k)) mc)
             (conj! cs c)))))

(set! *unchecked-math* :warn-on-boxed)

(defn eval-legendre-P
  "Evaluates the Legendre polynomial `P_n` at `x`.

  The polynomials are orthogonal on `[-1, 1]` with the weight 1 and satisfy `P_0 = 1`, `P_1 = x` and `n*P_n = (2n-1)*x*P_(n-1) - (n-1)*P_(n-2)`. They are polynomials in `x` for every real `x`, so outside `[-1, 1]` the value is the continuation of the polynomial. `P_n(1) = 1` and `P_n(-x) = (-1)^n*P_n(x)`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

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

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-legendre-P]] for values at a high degree.

  See also [[legendre-P-ratio]], [[eval-legendre-P]]."
  [^long degree]
  (polynomial (coeffs (legendre-P-ratio degree))))

;;

;; Gegenbauer polynomials. The three term recurrence is accurate for an order of -0.9 and above (at most 4
;; units on random orders, down to -0.9375) and loses digits below it (a negative integer order -k even gives
;; noise for the polynomials of degree above 2k, which are 0). There the explicit sum
;;   C_n^(a)(x) = sum_j (-1)^j (a)_(n-j) (2x)^(n-2j) / (j! (n-2j)!),  (a)_m = a (a+1) ... (a+m-1)
;; is evaluated in decimal arithmetic of adaptive precision.

(def ^:private ^:const gegenbauer-recurrence-min-order -0.9)

(defn- gegenbauer-decimal-sum
  "The explicit sum of the Gegenbauer polynomial for a finite `order` and `x`, evaluated with the precision `mc`."
  ^BigDecimal [^long degree ^double order ^double x ^MathContext mc]
  (let [alpha (BigDecimal. order)
        two-x (.multiply (BigDecimal. x) (BigDecimal. 2) mc)
        powers (decimal-powers two-x degree mc)
        rising (loop [m (long 0)
                      p BigDecimal/ONE
                      ps (transient [])]
                 (if (m/> m degree)
                   (persistent! ps)
                   (recur (m/inc m) (.multiply ^BigDecimal p (.add alpha (BigDecimal. m) mc) mc) (conj! ps p))))
        factorials (loop [k (long 0)
                          f BigDecimal/ONE
                          fs (transient [])]
                     (if (m/> k degree)
                       (persistent! fs)
                       (recur (m/inc k) (.multiply ^BigDecimal f (BigDecimal. (m/inc k)) mc) (conj! fs f))))]
    (loop [j (long 0)
           total BigDecimal/ZERO]
      (if (m/> (m/long-mult 2 j) degree)
        total
        (let [m (m/long-sub degree (m/long-mult 2 j))
              term (.divide (.multiply ^BigDecimal (rising (m/long-sub degree j)) ^BigDecimal (powers m) mc)
                            (.multiply ^BigDecimal (factorials j) ^BigDecimal (factorials m) mc)
                            mc)]
          (recur (m/inc j)
                 (if (m/even? j) (.add ^BigDecimal total ^BigDecimal term mc) (.subtract ^BigDecimal total ^BigDecimal term mc))))))))

(defn- gegenbauer-value-at-infinity
  "The limit of the Gegenbauer polynomial of a positive degree and a finite `order` at an infinite `x`.

  It is the infinity of the leading non-zero term, found from the signs alone. The term of `x^(n-2j)` has the
  coefficient `(-1)^j (a)_(n-j) 2^(n-2j) / (j! (n-2j)!)`, which vanishes when the order is a negative integer
  `-k` and `n-j > k`. Then the polynomial is 0 for `n > 2k` and the constant 1 for `n = 2k`."
  ^double [^long degree ^double order ^double x]
  (let [terminating? (and (m/<= order 0.0) (m/== order (m/floor order)))
        ;; the index j of the first non-vanishing term
        j (if (and terminating? (m/< (m/- order) degree)) (m/long-sub degree (long (m/- order))) 0)
        power (m/long-sub degree (m/long-mult 2 j))]
    (cond
      (m/neg? power) 0.0
      (m/zero? power) 1.0
      :else (let [factors (m/long-sub degree j)
                  ;; factors of (a)_(n-j) that are negative: those a+i, i < -a
                  negative-factors (long (m/max 0.0 (m/min (m/ceil (m/- order)) (double factors))))
                  positive-coefficient? (m/even? (m/long-add j negative-factors))
                  flipped-by-x? (and (m/neg? x) (m/odd? power))]
              (if (not= positive-coefficient? flipped-by-x?) ##Inf ##-Inf)))))

(defn eval-gegenbauer-C
  "Evaluates the Gegenbauer (ultraspherical) polynomial `C_n^(a)` of order `a` at `x`.

  The polynomials have the generating function `(1 - 2*x*t + t^2)^(-a)` and satisfy `C_0 = 1`, `C_1 = 2*a*x` and `n*C_n = 2*(n+a-1)*x*C_(n-1) - (n+2*a-2)*C_(n-2)`. For `a > -1/2` and `a` not 0 they are orthogonal on `[-1, 1]` with the weight `(1 - x^2)^(a - 1/2)`. They are polynomials in `x` for every real `x` and every real order, and `C_n(-x) = (-1)^n*C_n(x)`. The order 1 gives the Chebyshev polynomials of the second kind and the order 0.5 the Legendre polynomials; these two orders return exactly the values of [[eval-chebyshev-U]] and [[eval-legendre-P]].

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `order` (double): the order `a`, any real number, default 1.0 (the two-argument form). A negative order is allowed. For the order 0 the recurrence gives `C_0 = 1` and `C_n = 0` for every higher degree (the same values as `scipy` and `mpmath` give); the Chebyshev polynomials of the first kind are the limit of `C_n/a`, not the value.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `order` and `x` are. `##NaN` as `x` or as `order` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading non-zero term of the polynomial (so 0 for the order 0, and `##NaN` for an infinite order). The leading term is found from the signs alone, so that case is fast for any degree; the limit is 0 for a negative integer order `-k` and a degree above `2k`, and the constant 1 for the degree `2k`.

  For an order of -0.9 and above the three term recurrence is used and the error is a few units of roundoff times `n` (absolute, for values of at most `n+1`) inside `[-1, 1]`, and relative outside; unlike evaluating the coefficients from [[gegenbauer-C]] it stays accurate for any degree. Below -0.9 the recurrence loses digits (the polynomial is small against its intermediate values and a negative integer order gives noise for a polynomial that is 0), so the explicit sum is evaluated in decimal arithmetic whose precision is raised until the result settles: the result is then the nearest double to about 2^-60, and exactly 0 for a negative integer order `-k` and a degree above `2k`. That takes more time for a high degree (about 0.1 s for degree 1000, a second or more for degree 5000). The value is continuous in the order, also at 1 and 0.5.

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
                  (gegenbauer-value-at-infinity degree order x))
     (m/== degree 1) (m/* 2.0 order x)
     ;; the recurrence is not accurate for an order below -0.9: C_n is O(a + 1) for n >= 3 next to a = -1,
     ;; small against the intermediate values, and all the more below -1
     (and (m/< order gegenbauer-recurrence-min-order) (m/valid-double? order) (m/valid-double? x))
     (adaptive-decimal-value (fn [mc] (gegenbauer-decimal-sum degree order x mc)))
     :else (let [o2 (m/* 2.0 order)]
             (loop [i (long 2)
                    pprev 1.0
                    prev (m/* 2.0 order x)]
               (if (m/> i degree)
                 prev
                 (recur (m/inc i) prev
                        ;; the factors (a + i - 1) and (i - 2 + 2a) are formed from the small part first: adding
                        ;; i and subtracting afterwards loses the relative accuracy of a small order or of an
                        ;; order next to -1 (error about eps / |a| or eps / (a + 1))
                        (m// (m/- (m/* 2.0 (m/+ order (m/long-dec i)) x prev)
                                  (m/* (m/+ (m/long-sub i 2) o2) pprev)) i))))))))

(set! *unchecked-math* true)

(defn gegenbauer-C-ratio
  "Creates the Gegenbauer (ultraspherical) polynomial `C_n^(a)` with exact rational coefficients.

  `C_0 = 1`, `C_1 = 2*a*x` and `n*C_n = 2*(n+a-1)*x*C_(n-1) - (n+2*a-2)*C_(n-2)`; see [[eval-gegenbauer-C]] for the definition, the orthogonality and the special orders.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `order` (double): the order `a`, any real number, default 1.0 (the Chebyshev polynomials of the second kind). A NaN or infinite order throws an `IllegalArgumentException` for a positive degree; see [[gegenbauer-C-ratio]] for the conversion of the order.

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-gegenbauer-C]] for values at a high degree.

  See also [[gegenbauer-C-ratio]], [[eval-gegenbauer-C]], [[chebyshev-U]], [[legendre-P]]."
  ([^long degree] (gegenbauer-C degree 1.0))
  ([^long degree ^double order] (polynomial (coeffs (gegenbauer-C-ratio degree order)))))

;;

(defn eval-hermite-H
  "Evaluates the Hermite polynomial `H_n` (physicists' convention) at `x`.

  The polynomials are orthogonal on the whole real line with the weight `exp(-x^2)` and satisfy `H_0 = 1`, `H_1 = 2*x` and `H_n = 2*x*H_(n-1) - 2*(n-1)*H_(n-2)`; the leading coefficient is `2^n`. `H_n(-x) = (-1)^n*H_n(x)`. The probabilists' polynomials `He_n(x) = 2^(-n/2)*H_n(x/sqrt(2))` are given by [[eval-hermite-He]].

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading term: `##-Inf` for `##-Inf` and an odd degree, `##Inf` otherwise. The error is a few units of roundoff times `n` relative to the value (absolute next to a zero), so this form stays accurate for any degree, unlike evaluating the coefficients from [[hermite-H]].

  See also [[hermite-H]], [[hermite-H-ratio]], [[eval-hermite-He]], [[eval-laguerre-L]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (case (int degree)
    0 1.0
    1 (m/* 2.0 x)
    (if (m/inf? x)
      ;; the leading term is 2^n x^n
      (if (and (m/neg? x) (m/odd? degree)) ##-Inf ##Inf)
      (loop [i (long 2)
             pprev 1.0
             prev (m/* 2.0 x)]
        (if (m/> i degree)
          prev
          (recur (m/inc i) prev
                 (m/* 2.0 (m/- (m/* x prev)
                               (m/* (m/dec i) pprev)))))))))

(defn hermite-H-ratio
  "Creates the Hermite polynomial `H_n` (physicists' convention) with exact integer coefficients.

  `H_0 = 1`, `H_1 = 2*x` and `H_n = 2*x*H_(n-1) - 2*(n-1)*H_(n-2)`; see [[eval-hermite-H]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients of the powers of the wrong parity are zero, the signs alternate and the leading coefficient is `2^n`.

  See also [[hermite-H]] (double coefficients), [[eval-hermite-H]] (direct evaluation), [[hermite-He-ratio]]."
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
               (scale (sub (mult prev (ratio-polynomial [0 1]))
                           (scale pprev (m/dec i))) 2))))))

(defn hermite-H
  "Creates the Hermite polynomial `H_n` (physicists' convention) as a polynomial object with double coefficients.

  See [[eval-hermite-H]] for the definition and [[hermite-H-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-hermite-H]] for values at a high degree.

  See also [[hermite-H-ratio]], [[eval-hermite-H]], [[hermite-He]]."
  [^long degree]
  (polynomial (coeffs (hermite-H-ratio degree))))

(defn eval-hermite-He
  "Evaluates the Hermite polynomial `He_n` (probabilists' convention) at `x`.

  The polynomials are orthogonal on the whole real line with the weight `exp(-x^2/2)` and satisfy `He_0 = 1`, `He_1 = x` and `He_n = x*He_(n-1) - (n-1)*He_(n-2)`; the leading coefficient is 1. `He_n(-x) = (-1)^n*He_n(x)`. The physicists' polynomials `H_n(x) = 2^(n/2)*He_n(sqrt(2)*x)` are given by [[eval-hermite-H]].

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever `x` is. `##NaN` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading term: `##-Inf` for `##-Inf` and an odd degree, `##Inf` otherwise. The error is a few units of roundoff times `n` relative to the value (absolute next to a zero), so this form stays accurate for any degree, unlike evaluating the coefficients from [[hermite-He]].

  See also [[hermite-He]], [[hermite-He-ratio]], [[eval-hermite-H]]."
  ^double [^long degree ^double x]
  (check-degree! degree)
  (case (int degree)
    0 1.0
    1 x
    (if (m/inf? x)
      ;; the leading term is x^n
      (if (and (m/neg? x) (m/odd? degree)) ##-Inf ##Inf)
      (loop [i (long 2)
             pprev 1.0
             prev x]
        (if (m/> i degree)
          prev
          (recur (m/inc i) prev
                 (m/- (m/* x prev)
                      (m/* (m/dec i) pprev))))))))

(defn hermite-He-ratio
  "Creates the Hermite polynomial `He_n` (probabilists' convention) with exact integer coefficients.

  `He_0 = 1`, `He_1 = x` and `He_n = x*He_(n-1) - (n-1)*He_(n-2)`; see [[eval-hermite-He]] for the definition and the orthogonality.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `PolynomialR` (see [[ratio-polynomial]]) of degree `n`: operations on it are exact, and it can be evaluated exactly at a rational argument. The coefficients of the powers of the wrong parity are zero, the signs alternate and the leading coefficient is 1.

  See also [[hermite-He]] (double coefficients), [[eval-hermite-He]] (direct evaluation), [[hermite-H-ratio]]."
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
               (sub (mult prev (ratio-polynomial [0 1]))
                    (scale pprev (m/dec i))))))))

(defn hermite-He
  "Creates the Hermite polynomial `He_n` (probabilists' convention) as a polynomial object with double coefficients.

  See [[eval-hermite-He]] for the definition and [[hermite-He-ratio]] for the exact integer coefficients, which are converted to doubles.

  Parameters: `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).

  Returns a `Polynomial` (see [[polynomial]]) of degree `n`, which can be differentiated, multiplied and added. Evaluating the monomial form loses accuracy as the degree grows (the coefficients are large and alternate in sign), so use [[eval-hermite-He]] for values at a high degree.

  See also [[hermite-He-ratio]], [[eval-hermite-He]], [[hermite-H]]."
  [^long degree]
  (polynomial (coeffs (hermite-He-ratio degree))))


;;

;; Jacobi polynomials. The three term recurrence divides by 2 i (i + s) (2 i + s - 2) at step i, s = alpha + beta,
;; so it is undefined where s is a negative integer from -degree to -2, and it loses digits next to such an s
;; and, much more, for a parameter below -1 (thousands of units of roundoff). It is accurate for
;; alpha > -1, beta > -1 and s > -1.9 (at most 2.5 units on random parameters). Elsewhere the explicit sum
;;   P_n^(a,b)(x) = sum_k C(n+a, n-k) C(n+b, k) ((x-1)/2)^k ((x+1)/2)^(n-k)
;; is evaluated in decimal arithmetic of adaptive precision.

(def ^:private ^:const jacobi-recurrence-min-sum -1.9)

(defn- jacobi-recurrence-accurate?
  "True when the three term recurrence is accurate for the parameters: both above -1 and their sum above -1.9."
  [^double alpha ^double beta]
  (and (m/> alpha -1.0) (m/> beta -1.0) (m/> (m/+ alpha beta) jacobi-recurrence-min-sum)))

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

(defn- jacobi-decimal-sum
  "The explicit sum of the Jacobi polynomial for finite parameters and `x`, evaluated with the precision `mc`."
  ^BigDecimal [degree alpha beta x ^MathContext mc]
  (let [n (long degree)
        bx (BigDecimal. (double x))
        half (BigDecimal. "0.5")
        lowers (decimal-powers (.multiply (.subtract bx BigDecimal/ONE mc) half mc) n mc)
        uppers (decimal-powers (.multiply (.add bx BigDecimal/ONE mc) half mc) n mc)
        from-alpha (decimal-binomials (.add (BigDecimal. n) (BigDecimal. (double alpha)) mc) n mc)
        from-beta (decimal-binomials (.add (BigDecimal. n) (BigDecimal. (double beta)) mc) n mc)]
    (loop [k (long 0)
           total BigDecimal/ZERO]
      (if (m/> k n)
        total
        (let [rest-degree (m/long-sub n k)]
          (recur (m/inc k)
                 (.add ^BigDecimal total
                       (.multiply (.multiply ^BigDecimal (from-alpha rest-degree) ^BigDecimal (from-beta k) mc)
                                  (.multiply ^BigDecimal (lowers k) ^BigDecimal (uppers rest-degree) mc)
                                  mc)
                       mc)))))))

(defn- jacobi-explicit-value
  "The value of the Jacobi polynomial for finite parameters and `x` from the explicit sum, to double precision."
  ^double [^long degree ^double alpha ^double beta ^double x]
  (adaptive-decimal-value (fn [mc] (jacobi-decimal-sum degree alpha beta x mc))))

(defn- jacobi-constant-term
  "The binomial coefficient `C(n+a, n)`, the value of `P_n^(a,b)` at `x = 1`, as the nearest double."
  ^double [^long degree ^double alpha]
  (rational->double (peek (generalized-binomials (+' degree (exact-rational alpha)) degree))))

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

(defn- minus-integer-in?
  "True when `v` is an integer and `-v` is in `[lo, hi]`: the product of the factors `v + t` for `t` from `lo` to
  `hi` then has a zero factor."
  [^double v ^long lo ^long hi]
  (and (m/== v (m/floor v)) (m/>= (m/- v) lo) (m/<= (m/- v) hi)))

(defn- count-negative-factors
  "The number of integers `t` in `[lo, hi]` for which `v + t` is negative."
  ^long [^double v ^long lo ^long hi]
  (long (m/max 0.0 (m/min (m/- (m/ceil (m/- v)) lo) (double (m/long-inc (m/long-sub hi lo)))))))

(defn- jacobi-value-at-infinity
  "The limit of the Jacobi polynomial of a positive degree and finite parameters at an infinite `x`.

  It is the infinity of the leading non-zero term, found from the signs alone. In powers of `(x-1)/2` the
  coefficient of the power `j` is `(a+j+1)_(n-j) (n+a+b+1)_j / ((n-j)! j!)`, with the rising factorial
  `(c)_m = c (c+1) ... (c+m-1)`. The highest power with a non-zero coefficient decides; when it is the power
  0 the polynomial is the constant `C(n+a, n)`, and when no coefficient is non-zero it is 0."
  ^double [^long degree ^double alpha ^double beta ^double x]
  (let [s (m/+ alpha beta)]
    (loop [j degree]
      (let [after-j (m/long-inc j)
            end-of-s (m/long-add degree j)]
        (cond
          (m/neg? j) 0.0
          (or (minus-integer-in? alpha after-j degree)
              (minus-integer-in? s (m/long-inc degree) end-of-s)) (recur (m/long-dec j))
          (m/zero? j) (jacobi-constant-term degree alpha)
          :else (let [negatives (m/long-add (count-negative-factors alpha after-j degree)
                                            (count-negative-factors s (m/long-inc degree) end-of-s))
                      positive-coefficient? (m/even? negatives)
                      flipped-by-x? (and (m/neg? x) (m/odd? j))]
                  (if (not= positive-coefficient? flipped-by-x?) ##Inf ##-Inf)))))))

(defn eval-jacobi-P
  "Evaluates the Jacobi polynomial `P_n^(a,b)` at `x`.

  The polynomials satisfy `P_0 = 1`, `P_1 = (a+1) + (a+b+2)*(x-1)/2` and a three term recurrence in `n`. They are equal to the sum over `k` of `C(n+a, n-k)*C(n+b, k)*((x-1)/2)^k*((x+1)/2)^(n-k)`, with generalized binomial coefficients, which defines them for every real `a` and `b`. For `a > -1` and `b > -1` they are orthogonal on `[-1, 1]` with the weight `(1-x)^a*(1+x)^b`. `P_n(1) = C(n+a, n)` and `P_n^(a,b)(-x) = (-1)^n*P_n^(b,a)(x)`. For `a = b = 0` they are the Legendre polynomials, and for `a = b = c - 1/2` they are proportional to the Gegenbauer polynomials of order `c`.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
  - `alpha` (double): the parameter `a`, any real number.
  - `beta` (double): the parameter `b`, any real number.
  - `x` (double): the argument, any real number.

  Returns a double. The result for degree 0 is `1.0` whatever the other arguments are. `##NaN` as `x`, `alpha` or `beta` gives `##NaN` for a positive degree. An infinite `x` gives the infinity of the leading non-zero term of the polynomial (a finite constant when the polynomial is constant, as for degree 1 and `a + b = -2`); the leading term is found from the signs alone, so that case is fast for any degree.

  For `a > -1`, `b > -1` and `a + b > -1.9` the three term recurrence is used, and the error is a few units of roundoff times `n` (absolute, for values of at most `n+1`) inside `[-1, 1]`; it was at most 2.5 on random parameters of this range. Otherwise the recurrence is not used: it divides by zero when `a + b` is a negative integer from `-degree` to `-2` (or an even one down to `2 - 2*degree`) and loses digits, up to thousands of units, for a parameter below -1 or for `a + b` next to such a value. Then the sum above is evaluated in decimal arithmetic whose precision is raised until the result settles, so the result is the nearest double to about 2^-60, also next to a multiple root at 1 or -1 (a negative integer `a` or `b`). That takes more time for a high degree (about 0.2 s for degree 1000, a few seconds for degree 5000).

  See also [[jacobi-P]], [[jacobi-P-ratio]], [[eval-gegenbauer-C]], [[eval-legendre-P]]."
  ^double [^long degree ^double alpha ^double beta ^double x]
  (check-degree! degree)
  (cond
    (m/zero? degree) 1.0
    (m/inf? x) (if (or (m/invalid-double? alpha) (m/invalid-double? beta))
                 ##NaN
                 (jacobi-value-at-infinity degree alpha beta x))
    (m/== degree 1) (m/+ (m/inc alpha) (m/* 0.5 (m/+ alpha beta 2.0) (m/dec x)))
    (and (m/valid-double? x) (m/valid-double? alpha) (m/valid-double? beta)
         (not (jacobi-recurrence-accurate? alpha beta)))
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
                                 ;; (a - 1) = alpha + i - 1 formed from alpha: small next to alpha = -1
                                 (m/* 2.0 (m/+ alpha (m/long-dec i)) (m/+ beta (m/long-dec i)) c pprev))
                            (m/* 2.0 i (m/- c i) (m/- c 2.0)))))))))

(set! *unchecked-math* true)

(defn jacobi-P-ratio
  "Creates the Jacobi polynomial `P_n^(a,b)` with exact rational coefficients.

  See [[eval-jacobi-P]] for the definition, the orthogonality and the parameters.

  Parameters:

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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

  - `degree` (non-negative integer): the degree `n`. A degree that is negative or not below `Integer/MAX_VALUE` (2147483647) throws an `IllegalArgumentException`; a non-integer is truncated toward zero (so `##NaN` and `-0.5` give degree 0).
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
