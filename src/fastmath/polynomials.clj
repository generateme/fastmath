(ns fastmath.polynomials
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.complex :as cplx]
            [clojure.string :as str]
            [fastmath.protocols.polynomials :as prot])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]
           [java.text DecimalFormat]
           [clojure.lang IFn]
           [umontreal.ssj.functionfit PolInterp]
           [org.apache.commons.math3.linear MatrixUtils RealMatrix EigenDecomposition]
           [org.apache.commons.math3.special Gamma]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defmacro mevalpoly
  "Evaluate polynomial macro version in the form coeffs[0]+coeffs[1]*x+coeffs[2]*x^2+...."
  [x & coeffs]
  (let [cnt (count coeffs)]
    (condp m/== cnt
      0 0.0
      1 `~(first coeffs)
      2 (let [[z y] coeffs]
          `(m/muladd ~x ~y ~z))
      `(m/muladd ~x (mevalpoly ~x ~@(rest coeffs)) ~(first coeffs)))))

(defn evalpoly
  "Evaluate polynomial for given coefficients"
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
  "Create polynomial function for given coefficients"
  [coeffs]
  (cond
    (not (seq coeffs)) (constantly 0.0)
    (m/== 1 (count coeffs)) (constantly (first coeffs))
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
  "Evaluate complex polynomial macro version in the form coeffs[0]+coeffs[1]*x+coeffs[2]*x^2+...."
  [x & coeffs]
  (let [cnt (count coeffs)]
    (case (unchecked-int cnt)
      0 `cplx/ZERO
      1 (first coeffs)
      2 (let [[z y] coeffs]
          `(cplx/muladd ~x ~y ~z))
      `(cplx/muladd ~x (mevalpoly-complex ~x ~@(rest coeffs)) ~(first coeffs)))))

(defmacro mevalpoly-scalar-complex
  "Evaluate complex polynomial macro version in the form coeffs[0]+coeffs[1]*x+coeffs[2]*x^2+....

  Coefficients are scalars"
  [x & coeffs]
  (let [cnt (count coeffs)]
    (condp clojure.core/= cnt
      0 `cplx/ZERO
      1 `(Vec2. ~(first coeffs) 0.0)
      `(cplx/muladd ~x (mevalpoly-scalar-complex~x ~@(rest coeffs)) (Vec2. ~(first coeffs) 0.0)))))

(defn evalpoly-complex
  "Evaluate complex polynomial"
  [x & coeffs]
  (if-not (seq coeffs)
    cplx/ZERO
    (let [x (cplx/ensure-complex x)
          rc (reverse (map cplx/ensure-complex coeffs))]
      (loop [rcoeffs (rest rc)
             ex (first rc)]
        (if-not (seq rcoeffs)
          ex
          (recur (rest rcoeffs)
                 (cplx/muladd x ex (first rcoeffs))))))))

(defn makepoly-complex
  "Create complex polynomial function for given coefficients"
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
  (let [^DecimalFormat f (doto ^DecimalFormat (DecimalFormat/getInstance)
                           (.setMaximumFractionDigits 4)
                           (.setPositivePrefix "+")
                           (.setNegativePrefix "-"))
        numbers+var (apply str (take 20 (interleave (map (fn [^double v] (.format f v)) coeffs)
                                                    (degrees->vars degree))))
        res (str "#polynomial{" degree "}(x) = " (if (str/starts-with?  numbers+var "+")
                                                   (subs numbers+var 1) numbers+var))]
    (if (m/> degree 10) (str res "+...") res))  )


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
  prot/PolynomialProto
  (degree [_] d)
  (coeffs [_] (seq cfs))
  (add [p1 p2]
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
    (let [^Polynomial p2 p2]
      (cond
        (m/zero? d) (prot/scale p2 (Array/aget cfs 0))
        (m/zero? (.d p2)) (prot/scale p1 (Array/aget ^doubles (.cfs p2) 0))
        :else (let [nd (m/+ d (.d p2))
                    target (double-array (m/inc d))]
                (doseq [^long i1 (range (m/inc d))
                        ^long i2 (range (m/inc (.d p2)))]
                  (Array/add target (m/+ i1 i2) (m/* (Array/aget cfs i1)
                                                     (Array/aget ^doubles (.cfs p2) i2))))
                (Polynomial. target nd)))))
  (derivative [p order]
    (let [order (long order)]
      (if (m/zero? order) p
          (let [size (m/inc (m/- d order))
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
         (.equals cfs (.cfs ^PolynomialR poly))))
  (hashCode [_]
    (mix-collection-hash cfs d))
  IFn
  (invoke [p v] (prot/evaluate p v))
  prot/PolynomialProto
  (degree [_] d)
  (coeffs [_] cfs)
  (add [p1 p2]
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
    (let [^PolynomialR p2 p2]
      (cond
        (m/zero? d) (prot/scale p2 (cfs 0))
        (m/zero? (.d p2)) (prot/scale p1 ((.cfs p2) 0))
        :else (let [nd (m/+ d (.d p2))
                    target (transient (vec (repeat (m/inc nd) 0)))
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
      (if (m/zero? order) p
          (let [size (m/inc (m/- d order))]
            (loop [i (long 0)
                   pos order
                   fact (rationalize (m/factorial order))
                   target (transient (vec (repeat size 0)))]
              (if (m/== i size)
                (PolynomialR. (persistent! target) (m/dec size))
                (let [i+ (m/inc i)
                      pos+ (m/inc pos)]
                  (recur i+ pos+ (* (/ fact i+) pos+)
                         (assoc! target i (*' fact (cfs pos)))))))))))
  (evaluate [_ x]
    (let [rx (rationalize x)]
      (loop [i d
             ex (cfs i)]
        (if (m/zero? i)
          ex
          (let [i- (m/dec i)]
            (recur i- (+ (* rx ex) (cfs i-)))))))))

(set! *unchecked-math* :warn-on-boxed)

(def ^:private RONE (PolynomialR. [1] 0))

(defmethod print-method Polynomial [v ^java.io.Writer w] (.write w (str v)))
(defmethod print-method PolynomialR [v ^java.io.Writer w] (.write w (str v)))

(defn polynomial
  "Create polynomial object from coeffs or fitted points."
  (^Polynomial [xs ys]
   (polynomial (PolInterp/getCoefficients (m/seq->double-array xs) (m/seq->double-array ys))))
  (^Polynomial [coeffs]
   (Polynomial. (double-array coeffs) (m/dec (count coeffs)))))

(defn ratio-polynomial
  "Create polynomial operating on ratios from coeffs or fitted points."
  (^Polynomial [xs ys]
   (ratio-polynomial (PolInterp/getCoefficients (m/seq->double-array xs) (m/seq->double-array ys))))
  (^PolynomialR [coeffs]
   (PolynomialR. (mapv rationalize coeffs) (m/dec (count coeffs)))))

(defn coeffs->polynomial
  "Create polynomial object for unrolled coefficients."
  [& coeffs] (polynomial coeffs))

(defn coeffs->ratio-polynomial
  "Create ratio based polynomial object for unrolled coefficients."
  [& coeffs] (ratio-polynomial coeffs))

(defn add
  "Add two polynomials."
  ([poly] poly)
  ([poly1 poly2] (prot/add poly1 poly2)))

(defn sub
  "Subtract two polynomials"
  ([poly] (prot/negate poly))
  ([poly1 poly2]
   (prot/add poly1 (prot/negate poly2))))

(defn scale
  "Multiply polynomial by scalar"
  [poly v] (prot/scale poly v))

(defn mult
  "Multiply two polynomials."
  ([poly] poly)
  ([poly1 poly2]
   (prot/mult poly1 poly2)))

(defn coeffs
  "Coefficients of polynomial"
  [poly] (prot/coeffs poly))

(defn degree
  ^long [poly] (prot/degree poly))

(defn derivative
  "Derivative of the polynomial."
  ([poly] (derivative poly 1))
  ([poly order] (prot/derivative poly order)))

(defn evaluate
  "Evaluate polynomial"
  ^double [poly ^double x]
  (prot/evaluate poly x))

;; Orthogonal polynomials

(defn eval-bernstein
  ^double [^long degree ^long order ^double x]
  (case (int degree)
    0 1.0
    1 (if (m/zero? order) (m/- 1.0 x) x)
    (m/* (m/combinations degree order) (m/fpow x order) (m/fpow (m/- 1.0 x) (m/- degree order)))))


(defn bernstein
  [^long degree ^long order]
  (->> (range (m/inc degree))
       (map (fn [^long l]
              (if (m/< l order)
                0.0
                (m/* (if (m/even? (m/- l order)) 1.0 -1.0)
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

(defn- laguerre-L-ratio
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
  "Chebyshev polynomial of the first kind"
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 x
    2 (dec (* 2.0 x x))
    3 (* x (- (* 4.0 x x) 3.0))
    4 (let [x2 (* x x)] (inc (* 8.0 x2 (dec x2))))
    (m/cos (* degree (m/acos x)))))

(defn- chebyshev-T-ratio
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
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-T
  [^long degree]
  (polynomial (coeffs (chebyshev-T-ratio degree))))

(defn eval-chebyshev-U
  "Chebyshev polynomials of the second kind"
  ^double [^long degree ^double x]
  (case (int degree)
    0 1.0
    1 (m/* 2.0 x)
    2 (m/dec (m/* 4.0 x x))
    3 (m/* 4.0 x (m/dec (m/* 2.0 x x)))
    4 (let [x2 (m/* x x)] (m/inc (m/- (m/* 16.0 x2 x2) (m/* 12.0 x2))))
    (let [near-one (m/- 1.0 (m/* x x))
          degree+ (inc degree)]
      (if (m/< near-one (m// (m/sqrt 1.2E-8) (m/* degree+ degree+)))
        (let [v (m/* degree+ (m/- 1.0 (m/* m/SIXTH degree (m/+ degree 2) near-one)))]
          (if (and (m/odd? degree) (m/neg? x)) (m/- v) v))
        (let [t (m/acos x)]
          (m// (m/sin (m/* t degree+))
               (m/sin t)))))))

(defn- chebyshev-U-ratio
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
               (sub (mult prev (ratio-polynomial [0 2])) pprev))))))

(defn chebyshev-U
  [^long degree]
  (polynomial (coeffs (chebyshev-U-ratio degree))))

(defn eval-chebyshev-V
  "Chebyshev polynomials of the third kind"
  ^double [^long degree ^double x]
  (m/* (m/sqrt (m// 2.0 (m/inc x)))
       (eval-chebyshev-T (m/inc (m/* 2 degree)) (m/sqrt (m/* 0.5 (m/inc x))))))

(defn- chebyshev-V-ratio
  [^long degree]
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
  [^long degree]
  (polynomial (coeffs (chebyshev-V-ratio degree))))

(defn eval-chebyshev-W
  "Chebyshev polynomials of the fourth kind"
  ^double [^long degree ^double x]
  (m/* (eval-chebyshev-U (m/* 2 degree) (m/sqrt (m/* 0.5 (m/inc x))))))

(defn- chebyshev-W-ratio
  [^long degree]
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
  [^long degree]
  (polynomial (coeffs (chebyshev-W-ratio degree))))

;;

(defn eval-legendre-P
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
               (m// (m/- (m/* (m/dec (m/* 2.0 i)) x prev)
                         (m/* (m/dec i) pprev)) i))))))

(defn- legendre-P-ratio
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
               (scale (sub (mult prev (ratio-polynomial [0 (m/dec (m/* 2.0 i))]))
                           (scale pprev (m/dec i))) (/ 1 i)))))))

(defn legendre-P
  [^long degree]
  (polynomial (coeffs (legendre-P-ratio degree))))

;;

(defn eval-gegenbauer-C
  "Gegenbauer (ultraspherical) polynomials"
  (^double [^long degree ^double x] (eval-gegenbauer-C degree 1.0 x))
  (^double [^long degree ^double order ^double x]
   (condp == order
     1.0 (eval-chebyshev-U degree x)
     0.5 (eval-legendre-P degree x)
     (case (int degree)
       0 1.0
       1 (m/* 2.0 order x)
       (let [o2 (m/* 2.0 order)]
         (loop [i (long 2)
                pprev 1.0
                prev (m/* 2.0 order x)]
           (if (m/> i degree)
             prev
             (recur (m/inc i) prev
                    (m// (m/- (m/* 2.0 (m/dec (m/+ order i)) x prev)
                              (m/* (m/+ i o2 -2.0) pprev)) i)))))))))

(defn- gegenbauer-C-ratio
  [^long degree ^double order]
  (condp == order
    1.0 (chebyshev-U-ratio degree)
    0.5 (legendre-P-ratio degree)
    (case (int degree)
      0 RONE
      1 (ratio-polynomial [0 (m/* 2.0 order)])
      (let [o2 (m/* 2.0 order)]
        (loop [i (long 2)
               pprev RONE
               prev (ratio-polynomial [0 (m/* 2.0 order)])]
          (if (m/> i degree)
            prev
            (recur (m/inc i) prev
                   (scale (sub (mult prev (ratio-polynomial [0 (m/* 2.0 (m/dec (m/+ order i)))]))
                               (scale pprev (m/+ i o2 -2.0)))
                          (/ 1 i)))))))))

(defn gegenbauer-C
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

(defn- hermite-H-ratio
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

(defn- hermite-He-ratio
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

(defn eval-jacobi-P
  "Jacobi polynomials"
  ^double [^long degree ^double alpha ^double beta ^double x]
  (case (int degree)
    0 1.0
    1 (m/+ (m/inc alpha) (m/* 0.5 (m/+ alpha beta 2.0) (m/dec x)))
    (loop [i (long 2)
           pprev 1.0
           prev (m/+ (m/inc alpha) (m/* 0.5 (m/+ alpha beta 2.0) (m/dec x)))]
      (if (m/> i degree)
        prev
        (let [a (m/+ i alpha)
              b (m/+ i beta)
              c (m/+ a b)]
          (recur (m/inc i) prev
                 (m// (m/- (m/* (dec c) (m/+ (m/* c (m/- c 2.0) x)
                                             (m/* (m/- a b) (m/- c (m/* 2.0 i)))) prev)
                           (m/* 2.0 (m/dec a) (m/dec b) c pprev))
                      (m/* 2.0 i (m/- c i) (m/- c 2.0)))))))))

(set! *unchecked-math* true)

(defn- jacobi-P-ratio
  "Jacobi polynomials"
  [^long degree ^double alpha ^double beta]
  (case (int degree)
    0 RONE
    1 (let [ab22 (m/* 0.5 (m/+ alpha beta 2.0))]
        (ratio-polynomial [(m/- (m/inc alpha) ab22) ab22]))
    (let [alpha (rationalize alpha)
          beta (rationalize beta)]
      (loop [i (long 2)
             pprev RONE
             prev (let [ab22 (/ (+ alpha beta 2) 2)]
                    (ratio-polynomial [(- (inc alpha) ab22) ab22]))]
        (if (m/> i degree)
          prev
          (let [a (+ i alpha)
                b (+ i beta)
                c (+ a b)]
            (recur (m/inc i) prev
                   (scale (sub (scale (mult prev (ratio-polynomial [(* (- a b) (- c (* 2 i)))
                                                                    (* c (- c 2))])) (dec c))
                               (scale pprev (* 2 (dec a) (dec b) c))) (/ 1 (* 2 i (- c i) (- c 2)))))))))))

(set! *unchecked-math* :warn-on-boxed)

(defn jacobi-P
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


(defn- bessel-y-ratio
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

(defn- bessel-t-ratio
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

(defn- meixner-pollaczek-P-ratio
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
         order (m// m 2)
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
        order (m/dec (m// m 2))
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
        order (m// (m/dec m) 2)
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
        order (m// (m/dec m) 2)
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
   (assert (m/even? (m/- p m)) "p and m must be the same parity!")
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
   (assert (m/even? (m/- p m)) "p and m must be the same parity!")
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
   (assert (m/even? (m/- p m)) "p and m must be the same parity!")
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
   (assert (m/even? (m/- p m)) "p and m must be the same parity!")
   (let [^doubles coeffs (ince-S-coeffs p m e normalization)
         s (alength coeffs)]
     (if (m/even? p)
       (fn ^double [^double x] (ince-loop coeffs s x ince-sinh-even))
       (fn ^double [^double x] (ince-loop coeffs s x ince-sinh-odd))))))
