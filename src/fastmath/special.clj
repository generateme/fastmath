(ns fastmath.special
  "Special functions for real arguments and value.

  * Bessel J, Y, jinc
  * Modified Bessel I, K 
  * Spherical Bessel j, y
  * Modified spherical Bessel i1, i2, k
  * Elliptic K, E, Pi, Rf, Rd, Rj, Rg, Rc
  * Jacobi am, sn, cn, dn, sc, sd, cs, cd, ds, dc, ns, nc, nd
  * Gamma, log, digamma, trigamma, polygamma, regularized, lower/upper incomplete
  * Beta, log, regularized, incomplete
  * Erf, inverse
  * Airy A, B with derivatives
  * Zeta (Riemann, Hurwitz), Eta (Dirichlet), Xi (Landau), Beta (Dirichlet)
  * Integrals: Si, Ci, li/Li, Ei, En, Ein
  * Hypergeometric 0F0, 0F1, 1F0, 1F1, 2F1, 2F0, 0F2, pFq, Kummers M, Tricomis U, Whittaker M and W
  * Lambert W (0 and -1)
  * Minkowski
  * Harmonic H
  * Owen's T
  * Complex: (log)Gamma, Hypergeometric pFq, Tricomis U, (scaled) Bessel K of half-odd order."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special.poly :as spoly]
            [fastmath.special.airy :as airy]
            [fastmath.special.hypergeometric :as hg]
            [fastmath.special.ellip :as ellip]
            [fastmath.polynomials :as poly]
            [fastmath.complex :as cplx])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]
           [fastmath.special.hypergeometric PfQWenigerResultCplx]
           [org.apache.commons.math3.special Gamma Erf Beta]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; Erf

(defn erf
  "Error function.

  For two arguments returns a difference between `(erf x2)` and `(erf x1)`."
  {:inline (fn ([x] `(. Erf (erf (double ~x))))
             ([x1 x2] `(. Erf (erf (double ~x1) (double ~x2)))))
   :inline-arities #{1 2}}
  (^double [^double x] (. Erf (erf x)))
  (^double [^double x1 ^double x2] (. Erf (erf x1 x2))))

(defn erfc
  "Complementary error function."
  {:inline (fn [x] `(. Erf (erfc (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Erf (erfc x)))

(defn inv-erf
  "Inverse of [[erf]] function."
  {:inline (fn [x] `(. Erf (erfInv (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Erf (erfInv x)))

(defn inv-erfc
  "Inverse of [[erfc]] function."
  {:inline (fn [x] `(. Erf (erfcInv (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Erf (erfcInv x)))

;; Gamma

(defn gamma
  "Gamma function $\\Gamma(x)$. Extension of the factorial."
  {:inline (fn [x] `(. Gamma (gamma (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Gamma (gamma x)))

(defn log-gamma
  "Log of Gamma function $\\log\\Gamma(x)$."
  {:inline (fn [x] `(. Gamma (logGamma (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Gamma (logGamma x)))

(defn log-gamma-1p
  "$\\ln\\Gamma(1+x)$ for $-0.5≤x≤1.5$."
  {:inline (fn [x] `(. Gamma (logGamma1p (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Gamma (logGamma1p x)))

(defn inv-gamma-1pm1
  "$\\frac{1}{\\Gamma(1+x)}-1$ for $-0.5≤x≤1.5$."
  {:inline (fn [x] `(. Gamma (invGamma1pm1 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Gamma (invGamma1pm1 x)))

(defn digamma
  "First derivative of log of Gamma function."
  ^double [^double x]
  (let [^Vec2 psix (if (m/not-pos? x)
                     (Vec2. (m/* m/-PI (m/cot (m/* m/PI x))) (m/- 1.0 x))
                     (Vec2. 0.0 x))
        x (.y psix)
        ^Vec2 psix (if (m/< x 8.0)
                     (let [n (unchecked-long (m/- 8.0 (m/floor x)))]
                       (loop [v (long 1)
                              psi (.x psix)]
                         (if (m/== v n)
                           (Vec2. (m/- psi (m// x)) (m/+ x n))
                           (recur (m/inc v) (m/- psi (m// (m/+ x v)))))))
                     psix)
        t (m// (.y psix))
        psi (m/+ (.x psix) (m/- (m/log (.y psix)) (m/* 0.5 t)))
        t (m/* t t)]
    (m/- psi (m/* t (poly/mevalpoly t 0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686)))))

(defn trigamma
  "Second derivative of log of Gamma function."
  ^double [^double x]
  (if (m/not-pos? x)
    (m/- (m/sq (m// m/PI (m/sin (m/* m/PI x)))) (trigamma (m/- 1.0 x)))
    (let [^Vec2 psix (if (m/< x 10.0)
                       (let [n (unchecked-long (m/- 10.0 (m/floor x)))]
                         (loop [v (long 1)
                                psi (m/sq (m// x))]
                           (if (m/== v n)
                             (Vec2. psi (m/+ x n))
                             (recur (m/inc v) (m/+ psi (m/sq (m// 1.0 (m/+ x v))))))))
                       (Vec2. 0.0 x))
          t (m// (.y psix))
          w (m/* t t)
          psi (m/+ (.x psix) t (m/* 0.5 w))]
      (m/+ psi (m/* t w (poly/mevalpoly w 0.16666666666666666,-0.03333333333333333,0.023809523809523808,-0.03333333333333333,0.07575757575757576,-0.2531135531135531,1.1666666666666667,-7.092156862745098))))))

;; Beta

(defn- nonpos-int?
  [^double x] (and (m/not-pos? x) (m/integer? x)))

(defn beta
  "Computes the Beta function of `p` and `q`.

  The Beta function is defined as `Gamma(p)*Gamma(q)/Gamma(p+q)`, and is
  closely related to the Gamma function and to binomial coefficients.

  Parameters:

  - `p`, `q` (double): real arguments.

  Returns a double. Defined via analytic continuation for negative `p`
  and/or `q` as well. Returns `##NaN` where `p` or `q` is a nonpositive
  integer (a true pole), except that if `p+q` is also a nonpositive integer
  while neither `p` nor `q` individually is, the two singularities cancel
  and `0.0` is returned instead.

  See also [[log-beta]], [[regularized-beta]], [[incomplete-beta]],
  [[gamma]]."
  ^double [^double p ^double q]
  (let [s (m/+ p q)]
    (cond
      (and (m/pos? p) (m/pos? q)) (m/exp (. Beta (logBeta p q)))
      (and (nonpos-int? s) (not (nonpos-int? p)) (not (nonpos-int? q))) 0.0
      :else (m// (m/* (gamma p) (gamma q)) (gamma s)))))

(defn log-beta
  "Computes the natural logarithm of the absolute value of the Beta
  function of `p` and `q`.

  Parameters:

  - `p`, `q` (double): real arguments, same domain as [[beta]].

  Returns a double. Equal to the logarithm of [[beta]] wherever [[beta]] is
  positive (always the case for `p,q>0`). Elsewhere, since [[beta]] can be
  negative, zero, or undefined there, the result is `##NaN` where [[beta]]
  is negative or at its poles, and `##-Inf` at its removable zero.

  See also [[beta]]."
  ^double [^double p ^double q]
  (if (and (m/pos? p) (m/pos? q))
    (. Beta (logBeta p q))
    (m/log (beta p q))))

(defn regularized-beta
  "Regularized Beta I_x(a,b)"
  ^double [^double x ^double a ^double b]
  (if (and (m/pos? a) (m/pos? b))
    (. Beta (regularizedBeta x a b))
    (m// (m/* (m// (m/pow x a) a)
              (hg/hypergeometric-2F1 a (m/- 1.0 b) (m/inc a) x))
         (m// (m/* (gamma a) (gamma b)) (gamma (m/+ a b))))))

(defn incomplete-beta
  "Incomplete Beta B(x,a,b)"
  ^double [^double x ^double a ^double b]
  (if (and (m/pos? a) (m/pos? b))
    (m/exp (m/+ (m/log (. Beta (regularizedBeta x a b)))
                (. Beta (logBeta a b))))
    (m/* (m// (m/pow x a) a)
         (hg/hypergeometric-2F1 a (m/- 1.0 b) (m/inc a) x))))

;; zeta
;; https://github.com/JuliaMath/SpecialFunctions.jl/blob/master/src/gamma.jl

(defmacro ^:private horner
  [x m & p]
  (let [^doubles p (double-array p)
        cnt (count p)
        ka (m/dec (m/* 2 cnt))
        kb (m/dec ka)]
    (loop [k (m/dec cnt)
           ex `(m/* (m/+ ~m ~ka)
                    (m/+ ~m ~kb)
                    ~(m// (Array/aget p (m/dec cnt))
                          (m/* ka kb)))]
      (if (m/>= k 2)
        (let [ka (m/dec (m/* 2 k))
              kb (m/dec ka)
              cdiv (m// 1.0 (m/* ka kb))]
          (recur (m/dec k) `(m/* ~cdiv (m/+ ~m ~ka) (m/+ ~m ~kb)
                                 (m/+ ~(Array/aget p (m/dec k)) (m/* ~x ~ex)))))
        `(m/* (m/inc ~m) (m/+ ~(Array/aget p 0) (m/* ~x ~ex)))))))

(defn- zeta-sz-inner
  ^Vec2 [^double s ^double z ^double cutoff]
  (let [zf (m/floor z)
        nz (unchecked-int zf)
        n (m/ceil (m/- cutoff nz))
        -s (m/- s)
        zt (double (if (m/neg? nz)
                     (let [-z (m/- z)
                           -nz (m/- nz)
                           zt (m/pow -z -s)
                           zt (if (m/not== zf z) (m/+ zt (m/pow (m/- z nz) -s)) zt)]
                       (if (m/pos? s)
                         (loop [v (m/dec -nz) zt zt]
                           (let [nzt (m/+ zt (m/pow (m/- -z v) -s))]
                             (if (or (m/== nzt zt) (m/zero? v))
                               zt (recur (m/dec v) nzt))))
                         (loop [v (long 1) zt zt]
                           (let [nzt (m/+ zt (m/pow (m/- -z v) -s))]
                             (if (or (m/== nzt zt) (m/== v -nz))
                               zt (recur (m/inc v) nzt))))))
                     (m/pow z -s)))
        mnz (m/max 1 (m/- 1 nz))
        smnz (m/dec mnz)]
    (Vec2. (if (m/pos? s)
             (loop [v mnz zt zt]
               (let [nzt (m/+ zt (m/pow (m/+ z v) -s))]
                 (if (or (m/== nzt zt) (m/== v n))
                   zt (recur (m/inc v) nzt))))
             (loop [v (m/dec n) zt zt]
               (let [nzt (m/+ zt (m/pow (m/+ z v) -s))]
                 (if (or (m/== nzt zt) (m/== v smnz))
                   zt (recur (m/dec v) nzt)))))
           (m/+ z n))))

(defn- zeta-neg-s
  ^double [^double s ^double z]
  (let [s (m/- 1.0 s)
        z (if (m/neg? z) (m/- 1.0 (m/frac z)) (m/frac z))
        f (m/exp (m/- (m/+ m/M_LN2 (log-gamma s)) (m/* s m/LOG_TWO_PI)))
        f1 (m/* m/HALF_PI s)
        f2 (m/* m/TWO_PI z)
        -s (m/- s)]
    (loop [n (long 1)
           sum 0.0
           np 1.0]
      (if (m/< np m/MACHINE-EPSILON)
        (m/* f sum)
        (let [nn (m/inc n)]
          (recur nn (m/+ sum (m/* np (m/cos (m/- f1 (m/* f2 n))))) (m/pow nn -s)))))))

(defn- zeta-shift-z
  ^double [^double s ^double z]
  (cond
    (m/> z 1.0) (let [m (long z)
                      a (m/frac z)
                      -s (m/- s)]
                  (m/- (v/sum (map (fn [^long n] (m/pow (m/+ a n) -s)) (range m)))))
    (m/neg? z) (let [m (m/inc (m/abs (long z)))
                     -s (m/- s)]
                 (v/sum (map (fn [^long n] (m/pow (m/+ z n) -s)) (range m))))
    :else 0.0))

(defn zeta
  "Computes the Riemann zeta function of `s`, or the Hurwitz (generalized)
  zeta function of `s` and `z`.

  The Riemann zeta function is `sum(n^-s)` for `n` from `1` to infinity,
  analytically continued to the whole real line except for its pole at
  `s=1`. The Hurwitz zeta function generalizes it with an additional real
  offset `z`, as `sum((z+n)^-s)` for `n` from `0` to infinity; it reduces
  to the Riemann zeta function when `z` is `0` or `1`.

  Parameters:

  - `s` (double): the order.
  - `z` (double, two-argument arity only): the offset.

  Returns a double. `##NaN` at the pole `s=1` (single-argument arity), or
  where `z<0` together with a non-integer `s` (two-argument arity, a domain
  where no real result exists). Values that legitimately exceed double
  precision saturate to `##Inf` or `##-Inf`.

  Accuracy of the two-argument (Hurwitz) form degrades for very negative
  `s` (roughly below `-5`), particularly combined with `z` close to `0`;
  the single-argument (Riemann) form is accurate across its whole domain.

  See also [[eta]], [[dirichlet-beta]], [[xi]], [[polygamma]]."
  (^double [^double s]
   (cond
     (m/zero? s) -0.5
     (or (m/one? s) (m/nan? s) (m/neg-inf? s)) ##NaN
     (m/pos-inf? s) 1.0
     ;; trivial zeros at negative even integers, returned exactly: detecting
     ;; them via `(sin (* HALF_PI s))` in the reflection branch below loses
     ;; all precision for large |s| (HALF_PI*s is an imprecise multiple of
     ;; pi long before s gets large), corrupting the result instead of
     ;; giving (approximately) zero
     (and (nonpos-int? s) (m/even? (long s))) 0.0
     (m/< (m/abs s) 1.0e-3) (poly/mevalpoly s -0.5,
                                            -0.918938533204672741780329736405617639861,
                                            -1.0031782279542924256050500133649802190,
                                            -1.00078519447704240796017680222772921424,
                                            -0.9998792995005711649578008136558752359121)
     (m/< s 0.5) (let [oms (m/- 1.0 s)
                       zoms (zeta oms)
                       sinv (m/sin (m/* m/HALF_PI s))]
                   ;; computed in log-space (using log-gamma instead of gamma) to
                   ;; avoid premature intermediate overflow of `(gamma oms)` for
                   ;; very negative `s` (i.e. very large `oms`), which would
                   ;; otherwise poison an otherwise representable finite result
                   (if (or (m/zero? zoms) (m/zero? sinv))
                     0.0
                     (m/* (m/sgn zoms) (m/sgn sinv)
                          (m/exp (m/+ (m/log (m/abs zoms))
                                      (log-gamma oms)
                                      (m/log (m/abs sinv))
                                      (m/* s (m/log m/TWO_PI))
                                      (m/log m/INV_PI))))))
     :else (let [m (m/dec s)
                 zt (m/inc (m/+ (m/pow 0.5 s)
                                (m/pow m/THIRD s)
                                (m/pow 0.25 s)
                                (m/pow 0.2 s)
                                (m/pow m/SIXTH s)))
                 w (m/pow 0.14285714285714285 m) ;; (1/7)^m
                 zt (m/+ zt (m/* w (m/+ (m// 1.0 m) 0.07142857142857142)))] ;; 0.5*(1/7)
             ;; 1/49 = 0.02040816326530612
             (m/+ zt (m/* w 0.02040816326530612 (horner 0.02040816326530612 m 0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686,3.0539543302701198))))))
  (^double [^double s ^double z]
   (cond
     (or (m/zero? z) (m/one? z)) (zeta s)
     ;; zeta(s,-1) = zeta(s,0) + (-1)^-s; for s<>0 that's zeta(s) + (-1)^-s
     ;; (only real when s is a nonpositive integer, otherwise (-1)^-s is
     ;; complex, or, for s>0, the recursion's 0^-s term diverges/is a pole);
     ;; s=0 is a special case since 0^-s=0^0=1 by convention rather than 0,
     ;; handled directly via the closed form zeta(0,a)=0.5-a
     (m/== z -1.0) (cond
                     (m/zero? s) 1.5
                     (nonpos-int? s) (m/+ (zeta s) (m/pow -1.0 (m/- s)))
                     :else ##NaN)
     (m/== s 2.0) (trigamma z)
     (or (m/nan? s) (m/nan? z) (m/neg-inf? s)
         (and (m/inf? s) (m/neg? z))) ##NaN
     (and (m/inf? s) (m/> z 1.0)) 0.0
     (m/inf? s) ##Inf
     ;; for s<=-2, z<0 and non-integer s: zeta(s,z) is genuinely complex
     ;; (analytic continuation of z^-s for negative real z and non-integer
     ;; exponent leaves the reals), so no real value exists
     (and (m/<= s -2.0) (m/neg? z) (not (m/integer? s))) ##NaN
     (m/<= s -2.0) (let [shift (zeta-shift-z s z)]
                     (m/+ shift (zeta-neg-s s z)))
     :else (let [m (m/dec s)
                 cutoff (m/+ 7.0 m)
                 ^Vec2 ztz (if (m/< z cutoff) (zeta-sz-inner s z cutoff) (Vec2. 0.0 z))
                 t (m// 1.0 (.y ztz))
                 w (m/pow t m)
                 zt (m/+ (.x ztz) (m/* w (m/+ (m// 1.0 m) (m/* 0.5 t))))
                 t (m/* t t)]
             (m/+ zt (m/* w t (horner t m 0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686,3.0539543302701198)))))))

(defn eta
  "Dirichlet Eta function"
  ^double [^double x]
  (if (m/zero? x)
    0.5
    (let [dx (m/- 1.0 x)]
      (if (m/< (m/abs dx) 7.0e-3)
        (m/* 0.6931471805599453094172321214581765
             (poly/mevalpoly dx 1.0,
                             -0.23064207462156020589789602935331414700440,
                             -0.047156357547388879740146103148112380421254,
                             -0.002263576552598880778433550956278702759143568,
                             0.001081837223249910136105931217561387128141157))
        (m/* (m/- (zeta x)) (m/expm1 (m/* 0.6931471805599453094 dx)))))))

(defn dirichlet-beta
  "Dirichlet Beta function"
  ^double [^double x]
  (cond
    ;; removable singularity: both zeta(x,0.25) and zeta(x,0.75) individually
    ;; have a pole at x=1 (with equal residues, so they cancel); the closed
    ;; form pi/4 is used directly instead of relying on the (Inf - Inf)
    ;; cancellation
    (m/one? x) (m/* 0.25 m/PI)
    ;; the alternating L-series converges to 1 extremely quickly as x grows
    ;; (already indistinguishable from 1.0 at double precision by x~40); the
    ;; general formula below is a (0 * (Inf - Inf)) indeterminate form both
    ;; at x=Inf and, for large enough finite x, numerically too: `4^-x`
    ;; underflows to exactly 0.0 around x>745, while `zeta(x,0.25)` (whose
    ;; leading term is `4^x`) overflows around x>=512, i.e. for x in
    ;; [512, 745) the product of an exact-0.0 and an Inf-Inf NaN yields NaN
    ;; instead of the correct ~1.0
    (m/>= x 100.0) 1.0
    :else (m/* (m/exp (m/* -1.38629436111989061883 x))
               (m/- (zeta x 0.25) (zeta x 0.75)))))

(defn xi
  "Riemann (Landau's) Xi function"
  ^double [^double s]
  (cond
    (m/neg? s) (xi (m/- 1.0 s))
    (m/one? s) 0.5
    (m/zero? s) 0.5
    ;; xi grows without bound (monotonically, staying positive) as s -> +Inf;
    ;; -Inf is handled by the reflection above, recursing into this branch
    (m/pos-inf? s) ##Inf
    :else (let [hs (m/* 0.5 s)
                sm1 (m/dec s)
                z (zeta s)]
            ;; computed in log-space (using log-gamma instead of gamma) to
            ;; avoid premature intermediate overflow of `(gamma hs)` for
            ;; large `s`, which would otherwise poison an otherwise
            ;; representable finite result (`hs`, `Gamma(hs)` and
            ;; `INV_PI^hs` are always positive for `s>0`; only `(s-1)` and
            ;; `zeta(s)` can be negative)
            (if (or (m/zero? sm1) (m/zero? z))
              0.0
              (m/* (m/sgn sm1) (m/sgn z)
                   (m/exp (m/+ (m/log hs)
                               (m/log (m/abs sm1))
                               (m/* hs (m/log m/INV_PI))
                               (log-gamma hs)
                               (m/log (m/abs z)))))))))

(def ^:private cotderiv-q-memo
  (memoize
   (fn ^doubles [^long m]
     (case (int m)
       0 (double-array [1.0])
       1 (double-array [1.0 1.0])
       (let [^doubles q- (cotderiv-q-memo (m/dec m))
             d (m/dec (alength q-))]
         (if (m/odd? (m/long-dec m))
           (let [rm (m// 2.0 m)
                 ^doubles q (double-array (alength q-))]
             (Array/aset q d (m/* d rm (Array/aget q- d)))
             (dotimes [i d]
               (let [i+ (m/inc i)]
                 (Array/aset q i (m/* rm (m/+ (m/* i (Array/aget q- i))
                                              (m/* i+ (Array/aget q- i+)))))))
             q)
           (let [rm (m// 1.0 m)
                 ^doubles q (double-array (m/inc (alength q-)))
                 end (alength q-)]
             (Array/aset q 0 (m/* rm (Array/aget q- 0)))
             (Array/aset q end (m/* (m/inc (m/* 2.0 d)) rm (Array/aget q- d)))
             (dotimes [i d]
               (let [i+ (m/inc i)]
                 (Array/aset q i+ (m/* rm (m/+ (m/* (m/inc (m/* 2.0 i+)) (Array/aget q- i+))
                                               (m/* (m/inc (m/* 2.0 i)) (Array/aget q- i)))))))
             q)))))))

(def ^{:private true :tag "[[D"} cotderiv-q (into-array (map cotderiv-q-memo (range 100))))

(defn- cotderiv
  ^double [^long m ^double z]
  (cond
    (m/neg? m) ##NaN
    (m/zero? m) (m/* m/PI (m/cot (m/* m/PI z)))
    (m/< m 100) (let [^doubles q (Array/arrayget2d cotderiv-q m)
                      lq (alength q)
                      x (m/cot (m/* m/PI z))
                      y (m/* x x)]
                  (loop [i (long 2)
                         s (m/+ (Array/aget q 0)
                                (m/* (Array/aget q 1) y))
                         t y]
                    (if (m/== i lq)
                      (m/* (m/pow m/PI (m/inc m))
                           (if (m/odd? m) s (m/* x s)))
                      (let [newt (m/* t y)]
                        (recur (m/inc i) (m/+ s (m/* (Array/aget q i) newt)) newt)))))
    :else (let [p (m/inc m)
                z (m/- z (m/round z))]
            (loop [n (long 1)
                   s (m// 1.0 (m/pow z p))]
              (let [a (m/pow (m/+ z n) p)
                    b (m/pow (m/- z n) p)
                    news (m/+ s (m// (m/+ a b) (m/* a b)))]
                (if (m/== s news) s (recur (m/inc n) news)))))))

(defn polygamma
  "Computes the polygamma function of order `m` at a real argument `x`.

  The polygamma function is the `m`-th derivative of the digamma function,
  generalizing [[digamma]] (order `0`) and [[trigamma]] (order `1`) to
  arbitrary nonnegative integer order.

  Parameters:

  - `m` (long): the order, must be nonnegative.
  - `x` (double): real argument.

  Returns a double. `##NaN` for negative `m`. Has poles at `x` equal to `0`
  or any negative integer: at `x=0` this is `##Inf` or `##-Inf` (sign
  depending on the parity of `m`), while at negative-integer poles a very
  large finite value approximating the pole is returned instead. For
  `x>0`, results are positive for odd `m` and negative for even `m`, and
  decay towards `0` as `x` grows.

  Accuracy degrades for large even `m` combined with nonpositive `x`,
  particularly near half-integer `x`.

  See also [[digamma]], [[trigamma]], [[zeta]]."
  ^double [^long m ^double x]
  (if (m/neg? m)
    ##NaN
    (case (int m)
      0 (digamma x)
      1 (trigamma x)
      (let [s (m/inc m)
            lgs (log-gamma s)]
        (if (m/not-pos? x)
          ;; parity already folded into `inner` (added to zeta before
          ;; multiplying by -Gamma(s)); no further sign flip needed
          (let [v (cotderiv m x)
                inner (m/+ (zeta s (m/- 1.0 x)) (if (m/even? m) v (m/- v)))]
            (if (m/zero? inner)
              0.0
              (m/* (m/- (m/sgn inner)) (m/exp (m/+ (m/log (m/abs inner)) lgs)))))
          ;; zeta(s,x) is always positive here (s>1, x>0); parity applied
          ;; as a separate sign flip on the zeta*Gamma(s) product
          (let [zs (zeta s x)]
            (if (m/zero? zs)
              0.0
              (m/* (if (m/even? m) -1.0 1.0) (m/exp (m/+ (m/log zs) lgs))))))))))

;; https://github.com/JuliaMath/Bessels.jl/blob/master/src/BesselFunctions/besselk.jl

;; Bessel J

(defn bessel-J0
  "Bessel function of the first kind of order 0, J_0(x)"
  ^double [^double x]
  (let [x (m/abs x)]
    (cond
      (m/zero? x) 1.0
      (m/< x m/HALF_PI) (let [x2 (m/* x x)]
                          (poly/mevalpoly x2
                                          1.0, -0.25, 0.01562499999999994, -0.00043402777777725544, 6.781684026082576e-6,
                                          -6.781683757550061e-8, 4.709479394601058e-10, -2.4016837144506874e-12,
                                          9.104258208703104e-15))
      (m/< x 26.0) (let [n (m/long-dec (unchecked-long (m/* m/M_2_PI x)))
                         ^Vec2 root (spoly/j0-roots n)]
                     (spoly/j0-polys n (m/- x (.x root) (.y root))))
      (m/pos-inf? x) 0.0
      :else (let [xinv (m// x)
                  x2 (m/* xinv xinv)
                  ^Vec2 pq (if (m/< x 125)
                             (Vec2. (poly/mevalpoly x2 1.0 -0.0625, 0.103515625, -0.5428466796875,
                                                    5.848699569702148, -106.8867939710617, 2968.142937842757,
                                                    -116538.4796968361)
                                    (poly/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518
                                                    -23.47512774997287 535.640519510616 -17837.27968894748))
                             (Vec2. (poly/mevalpoly x2 1.0 -0.0625, 0.103515625, -0.5428466796875)
                                    (poly/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518)))
                  a (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq))
                  xn (m/* xinv (.y pq))
                  b (m/sin (m/+ x m/QUARTER_PI xn))]
              (m/* a b)))))

(defn bessel-J1
  "Bessel function of the first kind of order 1, J_1(x)"
  ^double [^double x]
  (let [s (m/sgn x)
        x (m/abs x)]
    (cond
      (m/zero? x) 0.0
      (m/<= x m/HALF_PI) (let [x2 (m/* x x)]
                           (m/* s x
                                (poly/mevalpoly x2
                                                0.5, -0.0624999999999989, 0.002604166666657291, -5.42534721917933e-5,
                                                6.781683542660179e-7, -5.651361336587487e-9, 3.36191211106159e-11,
                                                -1.4511302591871352e-13)))
      (m/< x 26.0) (let [n (m/long-dec (unchecked-int (m/* m/M_2_PI x)))
                         ^Vec2 root (spoly/j1-roots n)]
                     (m/* s (spoly/j1-polys n (m/- x (.x root) (.y root)))))
      (m/pos-inf? x) 0.0
      :else (let [xinv (m// x)
                  x2 (m/* xinv xinv)
                  ^Vec2 pq (if (m/< x 125)
                             (Vec2. (poly/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625 -7.739953994750977
                                                    132.7618242502213 -3543.303665366024 135394.2285691809)
                                    (poly/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054
                                                    30.6240119934082 -659.185221823779 21156.31404552781))
                             (Vec2. (poly/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625)
                                    (poly/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054)))
                  a (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq))
                  xn (m/* xinv (.y pq))
                  b (m/sin (m/+ x m/-QUARTER_PI xn))]
              (m/* s a b)))))

;; jinc-c4 (/ (* PI PI PI PI) 192.0)
;; jinc-c2 (/ (* PI PI) -8.0)

(defn jinc
  "Besselj1 devided by `x`"
  ^double [^double x]
  (if (m/< (m/abs x) 0.002)
    (let [x2 (m/* x x)]
      (poly/mevalpoly x2 1.0 -1.2337005501361697 0.5073390158020964))
    (let [pix (m/* m/PI x)]
      (m/* 2.0 (m// (bessel-J1 pix) pix)))))


;;;;;;


(defn- a-ap-asymptotic
  ^Vec2 [^double v ^double x]
  (cond
    (m/> x (m/* 5.0 v)) (spoly/a-ap-poly-10 v x)
    (m/> x (m/* 2.0 v)) (spoly/a-ap-poly-20 v x)
    :else (spoly/a-ap-poly-30 v x)))

(defn- bessel-jy-debye
  ^Vec2 [^double v ^double x]
  (let [vmx (m/* (m/+ v x) (m/- v x))
        vs (m/sqrt vmx)
        sqvs (m// (m/sqrt vs))
        n (m/muladd v (m/- (m/log (m// x (m/+ v vs)))) (m/- vs))
        coeff (Vec2. (m/* m/INV_SQRT2PI (m/exp (m/- n)) sqvs)
                     (m/* -1.0 m/SQRT_2_PI (m/exp n) sqvs))
        p (m// v vs)
        p2 (m// (m/* v v) vmx)
        ^Vec2 res (spoly/split-poly (m/- (m// p v)) (spoly/uk-poly-jn p2))]
    (v/emult coeff res)))

(defn- bessel-jy-large-argument
  ^Vec2 [^double v ^double x]
  (let [^Vec2 aap (a-ap-asymptotic v x)
        vp2 (m/* m/HALF_PI v)
        b (m// m/SQRT_2_PI (m/sqrt (m/* (.y aap) x)))
        s (m/sin vp2)
        c (m/cos vp2)
        sa (m/sin (.x aap))
        ca (m/cos (.x aap))
        cms (m/- c s)
        cps (m/+ c s)
        s1 (m/* cms ca)
        s2 (m/* cps sa)
        s3 (m/* cms sa)
        s4 (m/* cps ca)]
    (v/mult (Vec2. (m/+ s1 s2) (m/- s3 s4)) (m/* m/SQRT2_2 b))))

(defn- hankel-dabye
  ^Vec2 [^double v ^double x]
  (let [vmx (m/* (m/+ x v) (m/- x v))
        vs (m/sqrt vmx)
        sqvs (m// (m/sqrt vs))
        n (m/- vs (m/* v (m/acos (m// v x))) m/QUARTER_PI)
        coef (cplx/scale (cplx/exp (Vec2. 0.0 n)) (m/* m/SQRT_2_PI sqvs))
        p (m// v vs)
        p2- (m/- (m// (m/* v v) vmx))
        poly (if (m/< v (m/+ 5.0 (m/* 0.998 x) (m/* 10.542 (m/cbrt (m/- x)))))
               (spoly/uk-poly-10 p2-)
               (spoly/uk-poly-20 p2-))
        uk-yn (second (spoly/split-poly-c (v/div (Vec2. 0.0 (- p)) v) poly))]
    (cplx/mult coef uk-yn)))

(defn- bessel-j-power-series
  ^double [^double v ^double x]
  (let [hx (m/* 0.5 x)
        t2 (m/* hx hx)]
    (loop [i (long 0)
           out 0.0
           a (m// (Math/pow hx v) (gamma (m/inc v)))]
      (let [nout (m/+ a out)]
        (if (or (m/> i 3000)
                (m/< (m/abs a) (m/* (m/abs nout) m/MACHINE-EPSILON)))
          nout
          (let [i+ (m/inc i)]
            (recur i+ nout (m/* a t2 (m// -1.0 (m/* i+ (m/+ v i+)))))))))))

(defn- jy-debye-fit
  ^double [^double x]
  (m/max 15.0 (m/+ 2.0 (m/* 1.00035 x) (m/* 6.714 (m/cbrt x)))))

(defn- bessel-j-up-recurrence
  ^Vec2 [^double x ^Vec2 jn ^Vec2 nu]
  (let [x2 (m// 2.0 x)
        end (m/+ (.y nu) 0.5)]
    (loop [start (.x nu)
           jnu (.x jn)
           jnup1 (.y jn)]
      (if (m/< start end)
        (recur (m/inc start) (m/muladd (m/* start x2) jnu (m/- jnup1)) jnu)
        (Vec2. jnup1 jnu)))))

(defn- bessel-j-down-recurrence
  ^Vec2 [^double x ^Vec2 jn ^Vec2 nu]
  (let [x2 (m// 2.0 x)
        end (m/- (.y nu) 0.5)]
    (loop [start (.x nu)
           jnu (.x jn)
           jnup1 (.y jn)]
      (if (m/> start end)
        (recur (m/dec start) (m/muladd (m/* start x2) jnu (m/- jnup1)) jnu)
        (Vec2. jnup1 jnu)))))

(defn- bessel-j-recurrence
  ^double [^double nu ^double x]
  (let [debye-cutoff (m/ceil (jy-debye-fit x))
        nu-shift (unchecked-int (m/ceil (m/- debye-cutoff nu)))
        v (m/+ nu nu-shift)
        jnu (.x ^Vec2 (bessel-jy-debye v x))
        jnup1 (.x ^Vec2 (bessel-jy-debye (m/inc v) x))]
    (.x ^Vec2 (bessel-j-down-recurrence x (Vec2. jnu jnup1) (Vec2. v nu)))))

(defn- hankel-debye-fit
  ^double [^double x]
  (m/+ 0.2 x (m/* 7.435 (m/cbrt (m/- x)))))

(defn- bessel-j-positive-args
  ^double [^double v ^double x]
  (cond
    (m/> v (jy-debye-fit x)) (.x ^Vec2 (bessel-jy-debye v x))
    (m/> x (m/max 20.0 (m/* 1.65 v))) (.x ^Vec2 (bessel-jy-large-argument v x))
    (m/< v (hankel-debye-fit x)) (cplx/re (hankel-dabye v x))
    (or (m/< x 7.0)
        (m/> v (poly/mevalpoly x 2.0 0.109 0.062))) (bessel-j-power-series v x) 
    :else (bessel-j-recurrence v x)))

(declare bessel-y-positive-args)

(defn- bessel-j-integer-order
  ^double [^long order ^double x]
  (let [abs-v (m/abs order)
        abs-x (m/abs x)
        sgn (if (m/even? (long abs-v)) 1.0 -1.0)
        bessel-j-val (bessel-j-positive-args abs-v abs-x)]
    (if (m/not-neg? order)
      (if (m/not-neg? x) bessel-j-val (m/* sgn bessel-j-val))
      (if (m/not-neg? x)
        (m/* sgn bessel-j-val)
        (let [bessel-y-val (double (bessel-y-positive-args abs-v abs-x))
              piao (m/* m/PI abs-v)
              so (m/round (m/sin piao))
              co (m/round (m/cos piao))]
          (m/* sgn (m/- (m/* bessel-j-val co) (m/* bessel-y-val so))))))))

(defn bessel-J
  "Bessel function of the first kind of order v, J_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-J0 x)
    (m/one? order) (bessel-J1 x)
    (m/invalid-double? x) x
    (m/integer? order) (bessel-j-integer-order (long order) x)
    (m/neg? x) ##NaN
    (m/not-neg? order) (bessel-j-positive-args (m/abs order) (m/abs x))
    :else (let [ao (m/abs order)
                ax (m/abs x)
                j (bessel-j-positive-args ao ax)
                y (double (bessel-y-positive-args ao ax))
                piao (m/* m/PI ao)
                so (m/sin piao)
                co (m/cos piao)]
            (m/- (m/* j co) (m/* y so)))))

;; Bessel Y

(defn bessel-Y0
  "Bessel function of the second kind of order 0, Y_0(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##-Inf
    (m/pos-inf? x) 0.0
    (or (m/neg? x) (m/invalid-double? x)) ##NaN
    (m/< x 5.0) (let [z (m/* x x)
                      w (m// (poly/mevalpoly z -1.84950800436986690637E16, 4.42733268572569800351E16,
                                             -3.46628303384729719441E15, 8.75906394395366999549E13,
                                             -9.82136065717911466409E11, 5.43526477051876500413E9,
                                             -1.46639295903971606143E7, 1.55924367855235737965E4)
                             (poly/mevalpoly z 2.50596256172653059228E17, 3.17157752842975028269E15,
                                             2.02979612750105546709E13, 8.64002487103935000337E10,
                                             2.68919633393814121987E8, 6.26107330137134956842E5,
                                             1.04128353664259848412E3, 1.00000000000000000000E0))]
                  (m/+ w (m/* m/M_2_PI (m/log x) (bessel-J0 x))))
    (m/< x 25.0) (let [w (m// 5.0 x)
                       z (m/* w w)
                       p (m// (poly/mevalpoly z 9.99999999999999997821E-1, 5.30324038235394892183E0,
                                              8.74716500199817011941E0, 5.44725003058768775090E0,
                                              1.23953371646414299388E0, 8.28352392107440799803E-2,
                                              7.96936729297347051624E-4)
                              (poly/mevalpoly z 1.00000000000000000218E0, 5.30605288235394617618E0,
                                              8.76190883237069594232E0, 5.47097740330417105182E0,
                                              1.25352743901058953537E0, 8.56288474354474431428E-2,
                                              9.24408810558863637013E-4))
                       q (m// (poly/mevalpoly z -6.05014350600728481186E0, -5.14105326766599330220E1,
                                              -1.47077505154951170175E2, -1.77681167980488050595E2,
                                              -9.32060152123768231369E1, -1.95539544257735972385E1,
                                              -1.28252718670509318512E0, -1.13663838898469149931E-2)
                              (poly/mevalpoly z 2.42005740240291393179E2, 2.06209331660327847417E3,
                                              5.93072701187316984827E3, 7.24046774195652478189E3,
                                              3.88240183605401609683E3, 8.56430025976980587198E2,
                                              6.43178256118178023184E1, 1.00000000000000000000E0))
                       xn (m/- x m/QUARTER_PI)
                       s (m/sin xn)
                       c (m/cos xn)]
                   (m// (m/* m/SQRT_2_PI (m/+ (m/* p s) (m/* w q c)))
                        (m/sqrt x)))
    :else (let [xinv (m// x)
                x2 (m/* xinv xinv)
                ^Vec2 pq (if (m/< x 125.0)
                           (Vec2. (poly/mevalpoly x2 1.0 -0.0625 0.103515625 -0.5428466796875 5.848699569702148
                                                  -106.8867939710617 2968.142937842757 -116538.4796968361)
                                  (poly/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518
                                                  -23.47512774997287 535.640519510616 -17837.27968894748))
                           (Vec2. (poly/mevalpoly x2 1.0 -0.0625 0.103515625 -0.5428466796875)
                                  (poly/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518)))]
            (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq)
                 (m/sin (m/+ x m/-QUARTER_PI (m/* xinv (.y pq))))))))

(defn bessel-Y1
  "Bessel function of the second kind of order 1, Y_1(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##-Inf
    (m/pos-inf? x) 0.0
    (or (m/neg? x) (m/invalid-double? x)) ##NaN
    (m/< x 5.0) (let [z (m/* x x)
                      w (* x (m// (poly/mevalpoly z -7.78877196265950026825E17, 2.02439475713594898196E17,
                                                  -8.12770255501325109621E15, 1.14509511541823727583E14,
                                                  -6.47355876379160291031E11, 1.26320474790178026440E9)
                                  (poly/mevalpoly z 3.97270608116560655612E18, 6.87141087355300489866E16,
                                                  6.20557727146953693363E14, 3.88231277496238566008E12,
                                                  1.87601316108706159478E10, 7.34811944459721705660E7,
                                                  2.35564092943068577943E5, 5.94301592346128195359E2,
                                                  1.00000000000000000000E0)))]
                  (m/+ w (m/* m/M_2_PI (m/- (m/* (bessel-J1 x) (m/log x)) (m// x)))))
    (m/< x 25.0) (let [w (m// 5.0 x)
                       z (m/* w w)
                       p (m// (poly/mevalpoly z 1.00000000000000000254E0, 5.21451598682361504063E0,
                                              8.42404590141772420927E0, 5.11207951146807644818E0,
                                              1.12719608129684925192E0, 7.31397056940917570436E-2,
                                              7.62125616208173112003E-4)
                              (poly/mevalpoly z 9.99999999999999997461E-1, 5.20982848682361821619E0,
                                              8.39985554327604159757E0, 5.07386386128601488557E0,
                                              1.10514232634061696926E0, 6.88455908754495404082E-2,
                                              5.71323128072548699714E-4))
                       q (m// (poly/mevalpoly z 2.52070205858023719784E1, 2.11688757100572135698E2,
                                              5.97489612400613639965E2, 7.10856304998926107277E2,
                                              3.66779609360150777800E2, 7.58238284132545283818E1,
                                              4.98213872951233449420E0, 5.10862594750176621635E-2)
                              (poly/mevalpoly z 3.36093607810698293419E2, 2.82619278517639096600E3,
                                              7.99704160447350683650E3, 9.56231892404756170795E3,
                                              4.98641058337653607651E3, 1.05644886038262816351E3,
                                              7.42373277035675149943E1, 1.00000000000000000000E0))
                       xn (m/- x m/M_3PI_4)
                       s (m/sin xn)
                       c (m/cos xn)]
                   (m// (m/* m/SQRT_2_PI (m/+ (m/* p s) (m/* w q c)))
                        (m/sqrt x)))
    :else (let [xinv (m// x)
                x2 (m/* xinv xinv)
                ^Vec2 pq (if (m/< x 135.0)
                           (Vec2. (poly/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625 -7.739953994750977
                                                  132.7618242502213 -3543.303665366024 135394.2285691809)
                                  (poly/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054
                                                  30.6240119934082 -659.185221823779 21156.31404552781))
                           (Vec2. (poly/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625)
                                  (poly/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054)))]
            (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq)
                 (m/sin (m/- m/-QUARTER_PI x (m/* xinv (.y pq))))))))

(defn- bessel-y-power-series
  ^Vec2 [^double v ^double x]
  (let [hx (m/* 0.5 x)
        a (m/pow hx v)]
    (if (m/zero? a)
      (Vec2. ##-Inf a)
      (let [b (m// a)
            t2 (m/* hx hx)
            vpi (m/* m/PI v)
            s (m/sin vpi)
            c (m/cos vpi)]
        (loop [i (long 0)
               out 0.0
               out2 0.0
               a (m// a (gamma (m/inc v)))
               b (m// b (gamma (m/- 1.0 v)))]
          (let [nout (m/+ out a)
                nout2 (m/+ out2 b)]
            (if (or (m/> i 3000)
                    (m/< (m/abs b) (m/* (m/abs nout2) m/MACHINE-EPSILON)))
              (Vec2. (m// (m/- (m/* out c) out2) s) out)
              (let [i+ (m/inc i)]
                (recur i+ nout nout2
                       (m/* a t2 -1.0 (m// (m/* (m/+ v i+) i+)))
                       (m/* b t2 -1.0 (m// (m/* (m/- i+ v) i+))))))))))))

(defn- bessel-y-chebyshev-low-orders
  ^Vec2 [^double v ^double x]
  (let [x1 (m/dec (m// (m/* 2.0 (m/- x 6.0)) 13.0))
        v1 (m/dec v)
        v2 v
        a (double-array (map (fn [ws] (spoly/clenshaw-chebyshev x1 ws)) spoly/bessel-y-chebyshev-weights))]
    (Vec2. (spoly/clenshaw-chebyshev v1 a)
           (spoly/clenshaw-chebyshev v2 a))))

(defn- bessel-y-chebyshev
  ^Vec2 [^double v ^double x]
  (let [v-floor (m/frac v)
        ^Vec2 y (bessel-y-chebyshev-low-orders v-floor x)]
    (bessel-j-up-recurrence x (Vec2. (.y y) (.x y)) (Vec2. (m/inc v-floor) v))))

(defn- bessel-y-fallback
  ^Vec2 [^double v ^double x]
  (if (m/<= 6.0 x 19.0)
    (bessel-y-chebyshev v x)
    (let [shift (unchecked-int (m/- (m/ceil v) (m/floor (hankel-debye-fit x)) -4.0))
          v2 (m/max (m/- v shift) (m/inc (m/frac v)))]
      (bessel-j-up-recurrence x (Vec2. (cplx/im (hankel-dabye v2 x))
                                       (cplx/im (hankel-dabye (m/dec v2) x)))
                              (Vec2. v2 v)))))

(defn- bessel-y-positive-args
  ^double [^double v ^double x]
  (cond
    (and (m/integer? v)
         (m/< v 250)) (.x ^Vec2 (bessel-j-up-recurrence x (Vec2. (bessel-Y1 x) (bessel-Y0 x))
                                                        (Vec2. 1.0 v)))
    (m/> v (jy-debye-fit x)) (.y ^Vec2 (bessel-jy-debye v x))
    (m/> x (m/max 20.0 (m/* 1.65 v))) (.y ^Vec2 (bessel-jy-large-argument v x))
    (m/< v (hankel-debye-fit x)) (cplx/im (hankel-dabye v x))
    (or (m/< x 7.0)
        (m/> v (m/- (m/* 1.35 x) 4.5))) (.x ^Vec2 (bessel-y-power-series v x)) 
    :else (.x ^Vec2 (bessel-y-fallback v x))))

(defn- bessel-y-integer-order
  ^double [^long order ^double x]
  (let [ao (m/long-abs order)
        y (bessel-y-positive-args ao x)]
    (if (and (m/neg? order) (m/odd? ao)) (m/- y) y)))

(defn bessel-Y
  "Bessel function of the second kind of order v, Y_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-Y0 x)
    (m/one? order) (bessel-Y1 x)
    (or (m/nan? order) (m/nan? x) (m/neg? x)) ##NaN
    (m/integer? order) (bessel-y-integer-order (long order) x)
    (m/not-neg? order) (bessel-y-positive-args (m/abs order) x)
    :else (let [ao (m/abs order)
                aopi (m/* ao m/PI)
                y (bessel-y-positive-args (m/abs order) x)
                j (bessel-j-positive-args (m/abs order) x)]
            (m/+ (m/* y (m/cos aopi))
                 (m/* j (m/sin aopi))))))

;; Bessel K

;; N.M.Temme, On the numerical evaluation of the modified bessel function of the third kind
;; (formulas 1.6 and 1.9)
;; https://www.researchgate.net/publication/242441899_On_the_numerical_evaluation_of_the_modified_bessel_function_of_the_third_kind

(defn bessel-K-half-odd
  "Bessel K_a function for a = order/2

  Function accepts only odd integers for order"
  ^double [^long odd-numerator ^double x]
  (case (int odd-numerator)
    1 (m/* (m/sqrt (m// m/HALF_PI x)) (m/exp (m/- x)))
    3 (m/* (m/sqrt (m// m/HALF_PI x)) (m/exp (m/- x)) (m/inc (m// x)))
    (loop [i (long 5)
           ^Vec2 pair (let [b1 (m/* (m/sqrt (m// m/HALF_PI x)) (m/exp (m/- x)))
                            b3 (m/* b1 (m/inc (m// x)))]
                        (Vec2. b1 b3))]
      (if (m/> i odd-numerator)
        (.y pair)
        (recur (m/+ i 2) (Vec2. (.y pair) (m/+ (m/* (.y pair) (m// (m/- i 2.0) x))
                                               (.x pair))))))))

(defn bessel-K-half-odd-scaled
  "Bessel K_a function scaled by e^x for a = order/2

  Function accepts only odd integers for order"
  ^double [^long odd-numerator ^double x]
  (case (int odd-numerator)
    1 (m/* (m/sqrt (m// m/HALF_PI x)))
    3 (m/* (m/sqrt (m// m/HALF_PI x)) (m/inc (m// x)))
    (loop [i (long 5)
           ^Vec2 pair (let [b1 (m/* (m/sqrt (m// m/HALF_PI x)))
                            b3 (m/* b1 (m/inc (m// x)))]
                        (Vec2. b1 b3))]
      (if (m/> i odd-numerator)
        (.y pair)
        (recur (m/+ i 2) (Vec2. (.y pair) (m/+ (m/* (.y pair) (m// (m/- i 2.0) x))
                                               (.x pair))))))))


(defn bessel-K0
  "Modified Bessel function of the second kind of order 0, K_0(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##Inf
    (m/pos-inf? x) 0.0
    (or (m/nan? x) (m/neg? x)) ##NaN
    (m/<= x 1.0) (let [x2 (m/* x x)
                       a (m/* 0.25 x2)
                       s (m/muladd
                          (poly/mevalpoly a -1.372509002685546267e-1, 2.574916117833312855e-1,
                                          1.395474602146869316e-2, 5.445476986653926759e-4,
                                          7.125159422136622118e-6)
                          (m// (poly/mevalpoly a 1.000000000000000000e+00, -5.458333438017788530e-02,
                                               1.291052816975251298e-03, -1.367653946978586591e-05))
                          1.137250900268554688)
                       a (m/muladd s a 1.0)]
                   (m/muladd (m/- a)
                             (m/log x)
                             (poly/mevalpoly x2 1.159315156584124484e-01, 2.789828789146031732e-01,
                                             2.524892993216121934e-02, 8.460350907213637784e-04,
                                             1.491471924309617534e-05, 1.627106892422088488e-07,
                                             1.208266102392756055e-09, 6.611686391749704310e-12)))
    :else (let [rx (m// x)
                a (m/muladd
                   (poly/mevalpoly rx 2.533141373155002416e-1, 3.628342133984595192e0,
                                   1.868441889406606057e1, 4.306243981063412784e1,
                                   4.424116209627428189e1, 1.562095339356220468e1,
                                   -1.810138978229410898e0, -1.414237994269995877e0,
                                   -9.369168119754924625e-2)
                   (m// (poly/mevalpoly rx 1.000000000000000000e0, 1.494194694879908328e1,
                                        8.265296455388554217e1, 2.162779506621866970e2,
                                        2.845145155184222157e2, 1.851714491916334995e2,
                                        5.486540717439723515e1, 6.118075837628957015e0,
                                        1.586261269326235053e-1))
                   1.0)]
            (m/* (m// a (m/sqrt x)) (m/exp (m/- x))))))

(defn bessel-K1
  "Modified Bessel function of the second kind of order 1, K_1(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##Inf
    (m/pos-inf? x) 0.0
    (or (m/nan? x) (m/neg? x)) ##NaN
    (m/<= x 1.0) (let [x2 (m/* x x)
                       a (m/* 0.25 x2)
                       pq (m/muladd
                           (poly/mevalpoly a -3.62137953440350228e-3, 7.11842087490330300e-3,
                                           1.00302560256614306e-5, 1.77231085381040811e-6)
                           (m// (poly/mevalpoly a 1.00000000000000000e0, -4.80414794429043831e-2,
                                                9.85972641934416525e-4, -8.91196859397070326e-6))
                           8.69547128677368164e-2)
                       pq (m/muladd (m/* pq a) a (m/inc (m/* a 0.5)))
                       a (m/* pq x 0.5)
                       pq (m/muladd
                           (m// (poly/mevalpoly x2 -3.07965757829206184e-1, -7.80929703673074907e-02,
                                                -2.70619343754051620e-3, -2.49549522229072008e-5)
                                (poly/mevalpoly x2 1.00000000000000000e0, -2.36316836412163098e-2,
                                                2.64524577525962719e-4, -1.49749618004162787e-6))
                           x (m// x))]
                   (m/muladd a (m/log x) pq))
    :else (let [rx (m// x)
                a (m/muladd
                   (poly/mevalpoly rx -1.97028041029226295e-1, -2.32408961548087617e0,
                                   -7.98269784507699938e0, -2.39968410774221632e0,
                                   3.28314043780858713e1, 5.67713761158496058e1,
                                   3.30907788466509823e1, 6.62582288933739787e0,
                                   3.08851840645286691e-1)
                   (m// (poly/mevalpoly rx 1.00000000000000000e0, 1.41811409298826118e1,
                                        7.35979466317556420e1, 1.77821793937080859e2,
                                        2.11014501598705982e2, 1.19425262951064454e2,
                                        2.88448064302447607e1, 2.27912927104139732e0,
                                        2.50358186953478678e-2))
                   1.45034217834472656)]
            (m/* (m// a (m/sqrt x)) (m/exp (m/- x))))))

(defn- bessel-k-large-args
  ^double [^double v ^double x]
  (let [v2 (m/* 4.0 v v)
        invx (m// (m/* 8.0 x))]
    (loop [i (long 1)
           t 1.0
           s 1.0]
      (if (or (m/< (m/abs t) m/MACHINE-EPSILON)
              (m/> i 75))
        (m/* s (m/sqrt (m// m/HALF_PI x)))
        (let [newt (m/* t invx (m// (m/- v2 (m/sq (m/dec (m/* 2.0 i)))) i))]
          (recur (m/inc i) newt (m/+ s newt)))))))

(defn- levin-scale
  ^double [^long n ^long k]
  (let [n+k (m/+ n k)
        n+2k (m/+ n+k k)]
    (m// (m/* -1.0 (m/inc n+k) n+k)
         (m/* (m/inc n+2k) n+2k))))

(defn- levin-transform
  ^double [^long N ^doubles s ^doubles w]
  (dotimes [i N]
    (let [si (Array/aget s i)
          wi (Array/aget w i)]
      (Array/aset s i (m// si wi))
      (Array/aset w i (m// wi))))
  (let [len (m/dec N)]
    (dotimes [k len]
      (dotimes [i (m/- len k)]
        (let [i+ (m/long-inc i)
              ls (levin-scale i+ k)]
          (Array/aset s i (m/muladd (Array/aget s i) ls (Array/aget s i+)))
          (Array/aset w i (m/muladd (Array/aget w i) ls (Array/aget w i+)))))))
  (m// (Array/aget s 0)
       (Array/aget w 0)))


(defn- bessel-k-levin
  ^double [^long N ^double v ^double x]
  (let [v2 (m/* 4.0 v v)
        invx (m// (m/* 8.0 x))
        buff-s (double-array N)
        buff-w (double-array N)
        val (m/sqrt (m// m/HALF_PI x))]
    (loop [i (long 0)
           s 0.0
           t 1.0]
      (cond
        (m/zero? t) (m/* s val)
        (m/== i 16) (m/* (levin-transform N buff-s buff-w) val)
        :else (let [i+ (m/inc i)
                    new-s (m/+ s t)
                    b (m// (m/- v2 (m/sq (m/dec (m/* 2.0 i+)))) i+)              
                    new-t (m/* t invx b)]
                (Array/aset buff-s i new-s)
                (Array/aset buff-w i new-t)
                (recur (m/inc i) new-s new-t))))))

(defn- bessel-k-large-orders
  ^double [^double v ^double x]
  (let [z (m// x v)
        zs (m/hypot 1.0 z)
        n (m/+ zs (m/- (m/log z) (m/log1p zs)))
        coeff (m// (m/* m/SQRT_HALFPI (m/sqrt (m// v)) (m/exp (m/* -1.0 v n)))
                   (m/sqrt zs))
        p (m// zs)
        max-vx (m/max v x)
        p2 (m// (m/* v v)
                (m/muladd max-vx max-vx (m/sq (m/min v x))))
        ^Vec2 res (spoly/split-poly (m/- (m// p v)) (spoly/uk-poly-10 p2))]
    (m/* coeff (.y res))))

(defn- bessel-k-up-recurrence
  ^Vec2 [^double x ^Vec2 kv ^Vec2 se]
  (let [x2 (m// 2.0 x)
        end (m/+ (.y se) 0.5)]
    (loop [start (.x se)
           jnum1 (.x kv)
           jnu (.y kv)]
      (if (m/>= start end)
        (Vec2. jnum1 jnu)
        (recur (m/inc start)
               jnu
               (m/muladd (m/* start x2) jnu jnum1))))))

(defn- f0-local-expansion-v0
  ^double [^double v ^double x]
  (let [l2dx (m/- m/M_LN2 (m/log x))
        mu (m/* v l2dx)
        vv (m/* v v)
        mu2 (m/* mu mu)
        sp (poly/mevalpoly vv 1.0, 1.6449340668482264, 1.8940656589944918, 1.9711021825948702)
        g1 (poly/mevalpoly vv -0.5772156649015329, 0.04200263503409518, 0.042197734555544306)
        g2 (poly/mevalpoly vv 1.0, -0.6558780715202539, 0.16653861138229145)
        sh (poly/mevalpoly mu2 1.0, 0.16666666666666666, 0.008333333333333333, 0.0001984126984126984, 2.7557319223985893e-6)]
    (m/* sp (m/+ (m/* g1 (m/cosh mu))
                 (m/* g2 sh l2dx)))))

(defn- bessel-k-temme-series
  ^Vec2 [^double v ^double x]
  (let [z (m/* x 0.5)
        zz (m/* z z)
        zv (Math/pow z v)
        negv (m/- v)]
    (loop [k (long 1)
           fk (f0-local-expansion-v0 v x)
           pk (m// (poly/mevalpoly v 1.0, -0.5772156649015329, 0.9890559953279725, -0.23263776388631713)
                   (m/* 2.0 zv))
           qk (m/* (poly/mevalpoly negv 1.0, -0.5772156649015329, 0.9890559953279725, -0.23263776388631713)
                   0.5 zv)
           ck 1.0
           out-v 0.0
           out-vp1 0.0]
      (let [term-v (m/* ck fk)
            term-vp1 (m/* ck (m/- pk (m/* (m/dec k) fk)))
            new-out-v (m/+ out-v term-v)
            new-out-vp1 (m/+ out-vp1 term-vp1)]
        (if (or (m/> k 500)
                (and (m/< (m/abs term-v) m/MACHINE-EPSILON)
                     (m/< (m/abs term-vp1) m/MACHINE-EPSILON)))
          (Vec2. out-v (m// out-vp1 z))
          (recur (m/inc k)
                 (m// (m/+ (m/* k fk) pk qk)
                      (m/- (m/* k k) (m/* v v)))
                 (m// pk (m/- k v))
                 (m// qk (m/+ k v))
                 (m/* ck (m// zz k))
                 new-out-v
                 new-out-vp1))))))

(defn- bessel-k-power-series
  ^double [^double v ^double x]
  (let [gam (gamma v)
        ngam (m// m/PI (m/* (m/sin (m/* m/-PI (m/abs v))) gam v))
        x2 (m/* x x)]
    (loop [k (long 1)
           s1 0.0
           s2 0.0
           t1 1.0
           t2 1.0]
      (if (or (m/> k 80)
              (m/< (m/abs t1) m/MACHINE-EPSILON))
        (let [xpv (Math/pow (m/* 0.5 x) v)
              s (m/+ (m/* gam s1)
                     (m/* xpv xpv ngam s2))]
          (m// s (m/* 2.0 xpv)))
        (let [ns1 (m/+ s1 t1)
              ns2 (m/+ s2 t2)
              nt1 (m/* t1 (m// x2 (m/* 4.0 k (m/- k v))))
              nt2 (m/* t2 (m// x2 (m/* 4.0 k (m/+ k v))))]
          (recur (m/inc k) ns1 ns2 nt1 nt2))))))

(defn bessel-K
  "Modified Bessel function of the second kind and real order v, K_v(x)"
  ^double [^double order ^double x]
  (let [v (m/abs order)]
    (cond
      (m/zero? v) (bessel-K0 x)
      (m/one? v) (bessel-K1 x)
      (m/nan? x) ##NaN
      (m/zero? x) ##Inf
      (m/neg? x) ##NaN
      (m/> x (m/+ 18.0 (m// (m/* v v) 36.0))) (m/* (m/exp (m/- x)) (bessel-k-large-args v x))
      (m/> x (m/+ 1.5 (m// (m/sq (m/sq v)) 2401.0))) (m/* (m/exp (m/- x)) (bessel-k-levin 16 v x))
      (or (m/> v 25.0) (m/> x 35)) (bessel-k-large-orders v x)
      :else (let [v-floor (m/frac v)]
              (cond
                (m/> x 1.5) (let [v-floor+1 (m/inc v-floor)
                                  kv (bessel-k-levin 16 v-floor x)
                                  kvp1 (bessel-k-levin 16 v-floor+1 x)
                                  ^Vec2 res (bessel-k-up-recurrence x (Vec2. kv kvp1) (Vec2. v-floor+1 v))]
                              (m/* (m/exp (m/- x)) (.x res)))
                
                (m/< (m/abs (m/- v (m/rint v))) 1.0e-5)
                (let [v-floor (if (m/> v-floor 0.5) (m/dec v-floor) v-floor)
                      kv (bessel-k-temme-series v-floor x)
                      ^Vec2 res (bessel-k-up-recurrence x kv (Vec2. (m/inc v-floor) v))]
                  (.x res))
                
                :else (bessel-k-power-series v x))))))

;; Bessel I

(defn bessel-I0
  "Modified Bessel function of the first kind of order 0, I_0(x)"
  ^double [^double x]
  (cond
    (m/invalid-double? x) ##NaN
    (m/zero? x) 1.0
    :else (let [x (m/abs x)]
            (if (m/< x 7.75)
              (let [a (m/* 0.25 x x)]
                (m/muladd a (poly/mevalpoly a  0.9999999999999998, 0.2500000000000052, 0.027777777777755364,
                                            0.001736111111149161, 6.94444444107536e-5, 1.9290123635366806e-6,
                                            3.9367592765038015e-8, 6.151201574092085e-10, 7.593827956729909e-12,
                                            7.596677643342155e-14, 6.255282299620455e-16, 4.470993793303175e-18,
                                            2.1859737023077178e-20, 2.0941557335286373e-22) 1.0))
              (let [invx (m// x)]
                (m/* (m/exp x)
                     (m// (poly/mevalpoly invx 0.3989422804014326, 0.04986778505064754, 0.028050628512954097,
                                          0.02921968830978531, 0.04466889626137549, 0.10220642174207666,
                                          -0.9937439085650689, 91.25330271974727, -4901.408890977662,
                                          199209.2752981982, -6.181516298413396e6, 1.4830278710991925e8,
                                          -2.7695254643719645e9, 4.0351394830842026e10, -4.5768930327229974e11,
                                          4.0134844243070063e12, -2.6862476523182016e13, 1.3437999451218112e14,
                                          -4.856333741437621e14, 1.1962791200680235e15, -1.796269414464399e15,
                                          1.239942074380968e15)
                          (m/sqrt x))))))))

(defn bessel-I1
  "Modified Bessel function of the first kind of order 1, I_0(x)"
  ^double [^double x]
  (cond
    (m/invalid-double? x) ##NaN
    (m/zero? x) 0.0
    :else (let [z (m/abs x)
                z (double
                   (if (m/< z 7.75)
                     (let [a (m/* 0.25 z z)
                           inner (poly/mevalpoly a 0.08333333333333334, 0.006944444444444374,
                                                 0.00034722222222248526, 1.1574074073690356e-5,
                                                 2.7557319253050506e-7, 4.920949730519126e-9,
                                                 6.834656365321179e-11, 7.593985414952446e-13,
                                                 6.904652315442046e-15, 5.2213850252454655e-17,
                                                 3.405120412140281e-19, 1.6398527256182257e-21,
                                                 1.3161876924566675e-23)]
                       (m/* 0.5 z (poly/mevalpoly a 1.0 0.5 inner)))
                     (let [invz (m// z)]
                       (m/* (m/exp z)
                            (m// (poly/mevalpoly invz 0.39894228040143276, -0.149603355151029,
                                                 -0.04675104787903509, -0.04090746353279043,
                                                 -0.05744911840910781,-0.12283724006390319,
                                                 1.0023632527650936, -94.90954045770921,
                                                 5084.06899084327, -206253.5613716743,
                                                 6.387439538535799e6, -1.529244018923123e8,
                                                 2.849523551208316e9, -4.141884344471782e10,
                                                 4.6860149658304974e11, -4.097852944580042e12,
                                                 2.7345051110005453e13, -1.3634833112030364e14,
                                                 4.909983186948892e14, -1.2048200837913132e15,
                                                 1.8014682382937435e15, -1.2377987428989558e15) 
                                 (m/sqrt z))))))]
            (if (m/neg? x) (m/- z) z))))

(defn- bessel-i-large-args
  ^double [^double v ^double x]
  (let [-invx (m// (m/* -8.0 x))
        v2 (m/* 4.0 v v)]
    (loop [i (long 1)
           t 1.0
           s 1.0]
      (let [new-t (m/* t -invx (m// (m/- v2 (m/sq (m/dec (m/* 2.0 i)))) i))
            new-s (m/+ s new-t)]
        (if (or (m/> i 1000)
                (m/< (m/abs new-t) m/MACHINE-EPSILON))
          (let [exh (m/exp (m/* 0.5 x))]
            (m/* exh (m// (m/* new-s) (m/sqrt (m/* m/TWO_PI x))) exh))
          (recur (m/inc i) new-t new-s))))))

(defn- bessel-i-large-orders
  ^double [^double v ^double x]
  (let [z (m// x v)
        zs (m/hypot 1.0 z)
        n (m/+ zs (m/- (m/log z) (m/log1p zs)))
        coeff (m// (m/* m/INV_SQRT2PI (m/sqrt (m// v)) (m/exp (m/* v n)))
                   (m/sqrt zs))
        p (m// zs)
        max-vx (m/max v x)
        p2 (m// (m/* v v)
                (m/muladd max-vx max-vx (m/sq (m/min v x))))
        ^Vec2 res (spoly/split-poly (m/- (m// p v)) (spoly/uk-poly-10 p2))]
    (m/* coeff (.x res))))

(defn- bessel-i-power-series
  ^double [^double v ^double x]
  (let [xx (m/* 0.25 x x)]
    (loop [i (long 0)
           s 0.0
           t 1.0]
      (let [new-s (m/+ s t)]
        (if (or (m/> i 3000)
                (m/<= (m/abs t) (m/* s m/MACHINE-EPSILON)))
          (m/* new-s (m// (Math/pow (m/* 0.5 x) v)
                          (gamma (m/inc v))))
          (let [i+ (m/inc i)]
            (recur i+
                   new-s
                   (m/* t (m// xx (m/* (m/+ v i+) i+))))))))))

(defn- bessel-i-positive-args
  ^double [^double v ^double x]
  (cond
    (m/> x (m/+ 19.0 (m/* 0.5 v v))) (bessel-i-large-args v x)
    (or (m/> v 25.0) (m/> x 35)) (bessel-i-large-orders v x)
    :else (bessel-i-power-series v x)))

(defn- bessel-i-integer-order
  ^double [^long v ^double x]
  (let [bessel-i-val (bessel-i-positive-args v (m/abs x))]
    (if (m/not-neg? x)
      bessel-i-val
      (if (m/even? v)
        bessel-i-val
        (m/- bessel-i-val)))))

(defn bessel-I
  "Modified Bessel function of the first kind of order v, I_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-I0 x)
    (m/one? order) (bessel-I1 x)
    (m/invalid-double? x) x
    (m/integer? order) (bessel-i-integer-order (unchecked-long (m/abs order)) x)
    (m/neg? x) ##NaN
    (m/zero? x) 0.0
    (m/not-neg? order) (bessel-i-positive-args (m/abs order) x)
    :else (let [a (m/abs order)]
            (m/+ (bessel-i-positive-args a x)
                 (m/* m/M_2_PI (m/sin (m/* m/PI a)) (bessel-K a x))))))

;; spherical

(defn spherical-bessel-j0
  "Spherical Bessel function of the first kind and order 0."
  ^double [^double x]
  (if (m/zero? x) 1.0
      (m// (m/sin x) x)))

(defn spherical-bessel-j1
  "Spherical Bessel function of the first kind and order 1."
  ^double [^double x]
  (if (m/zero? x) 0.0
      (m/- (m// (m/sin x) (m/* x x))
           (m// (m/cos x) x))))

(defn spherical-bessel-j2
  "Spherical Bessel function of the first kind and order 2."
  ^double [^double x]
  (if (m/zero? x) 0.0
      (let [x32 (m// 3.0 (m/* x x))]
        (m/- (m/* (m/- (m// x32 x) (m// x)) (m/sin x))
             (m/* x32 (m/cos x))))))

(defn spherical-bessel-j
  "Spherical Bessel function of the first kind."
  ^double [^double order ^double x]
  (condp m/== order
    0.0 (spherical-bessel-j0 x)
    1.0 (spherical-bessel-j1 x)
    2.0 (spherical-bessel-j2 x)
    (m/* (m/sqrt (m// m/HALF_PI x)) (bessel-J (m/+ order 0.5) x))))

(defn spherical-bessel-y0
  "Spherical Bessel function of the second kind and order 0."
  ^double [^double x]
  (if (m/zero? x) ##-Inf
      (m/- (m// (m/cos x) x))))

(defn spherical-bessel-y1
  "Spherical Bessel function of the second kind and order 1."
  ^double [^double x]
  (if (m/zero? x) ##-Inf
      (m/- (m/- (m// (m/cos x) (m/* x x)))
           (m// (m/sin x) x))))

(defn spherical-bessel-y2
  "Spherical Bessel function of the second kind and order 2."
  ^double [^double x]
  (if (m/zero? x) ##-Inf
      (let [x32 (m// 3.0 (m/* x x))]
        (m/- (m/* (m/- (m// x) (m// x32 x)) (m/cos x))
             (m/* x32 (m/sin x))))))

(defn spherical-bessel-y
  "Spherical Bessel function of the second kind."
  ^double [^double order ^double x]
  (condp m/== order
    0.0 (spherical-bessel-y0 x)
    1.0 (spherical-bessel-y1 x)
    2.0 (spherical-bessel-y2 x)
    (m/* (m/sqrt (m// m/HALF_PI x)) (bessel-Y (m/+ order 0.5) x))))

(defn spherical-bessel-1-i0
  "First modified spherical Bessel function of the first kind and order 0."
  ^double [^double x]
  (if (m/zero? x) 1.0
      (m// (m/sinh x) x)))

(defn spherical-bessel-1-i1
  "First modified spherical Bessel function of the first kind and order 1."
  ^double [^double x]
  (if (m/zero? x) 0.0
      (m/- (m// (m/cosh x) x) (m// (m/sinh x) (m/* x x)))))

(defn spherical-bessel-1-i2
  "First modified spherical Bessel function of the first kind and order 2."
  ^double [^double x]
  (if (m/zero? x) 0.0
      (let [x32 (m// 3.0 (m/* x x))]
        (m/- (m/* (m/+ (m// x) (m// x32 x)) (m/sinh x))
             (m/* x32 (m/cosh x))))))

(defn spherical-bessel-1-i
  "First modified spherical Bessel function of the first kind."
  ^double [^double order ^double x]
  (condp m/== order
    0.0 (spherical-bessel-1-i0 x)
    1.0 (spherical-bessel-1-i1 x)
    2.0 (spherical-bessel-1-i2 x)
    (m/* (m/sqrt (m// m/HALF_PI x)) (bessel-I (m/+ order 0.5) x))))

(defn spherical-bessel-2-i0
  "Second modified spherical Bessel function of the first kind and order 0."
  ^double [^double x] (m// (m/cosh x) x))

(defn spherical-bessel-2-i1
  "Second modified spherical Bessel function of the first kind and order 1."
  ^double [^double x] (m/- (m// (m/sinh x) x) (m// (m/cosh x) (m/* x x))))

(defn spherical-bessel-2-i2
  "Second modified spherical Bessel function of the first kind and order 2."
  ^double [^double x] (let [x32 (m// 3.0 (m/* x x))]
                        (m/- (m/* (m/+ (m// x) (m// x32 x)) (m/cosh x))
                             (m/* x32 (m/sinh x)))))

(defn spherical-bessel-2-i
  "Second modified spherical Bessel function of the first kind."
  ^double [^double order ^double x]
  (condp m/== order
    0.0 (spherical-bessel-2-i0 x)
    1.0 (spherical-bessel-2-i1 x)
    2.0 (spherical-bessel-2-i2 x)
    (m/* (m/sqrt (m// m/HALF_PI x)) (bessel-I (m/- (m/+ order 0.5)) x))))

(defn spherical-bessel-k0
  "Modified spherical Bessel function of the second kind and order 0."
  ^double [^double x] (m/* m/HALF_PI (m// (m/exp (m/- x)) x)))

(defn spherical-bessel-k1
  "Modified spherical Bessel function of the second kind and order 1."
  ^double [^double x] (m/* m/HALF_PI (m/exp (m/- x)) (m/+ (m// x) (m// (m/* x x)))))

(defn spherical-bessel-k2
  "Modified spherical Bessel function of the second kind and order 2."
  ^double [^double x] (let [x32 (m// 3.0 (m/* x x))]
                        (m/* m/HALF_PI (m/exp (m/- x)) (m/+ (m// x) x32 (m// x32 x)))))

(defn spherical-bessel-k
  "Modified spherical Bessel function of the second kind."
  ^double [^double order ^double x]
  (condp m/== order
    0.0 (spherical-bessel-k0 x)
    1.0 (spherical-bessel-k1 x)
    2.0 (spherical-bessel-k2 x)
    (m/* (m/sqrt (m// m/HALF_PI x)) (bessel-K (m/+ order 0.5) x))))

;;

(defn hankel-1
  "Hankel function of the first kind, returns complex number."
  ^Vec2 [^double order ^double x]
  (Vec2. (bessel-J order x) (bessel-Y order x)))

(defn hankel-2
  "Hankel function of the second kind, returns complex number."
  ^Vec2 [^double order ^double x]
  (Vec2. (bessel-J order x) (m/- (bessel-Y order x))))

(defn spherical-hankel-1
  "Spherical Hankel function of the first kind, returns complex number."
  ^Vec2 [^double order ^double x]
  (Vec2. (spherical-bessel-j order x) (spherical-bessel-y order x)))

(defn spherical-hankel-2
  "Spherical Hankel function of the second kind, returns complex number."
  ^Vec2 [^double order ^double x]
  (Vec2. (spherical-bessel-j order x) (m/- (spherical-bessel-y order x))))


;;

(defn minkowski
  "Minkowski's question mark function ?(x)"
  (^double [^double x]
   (loop [it (long 0) p 0.0 q 1.0 r 1.0 s 1.0 d 1.0 y 0.0]
     (if (m/< d (m/ulp y))
       (m/+ y d)
       (let [d (m/* d 0.5)
             m (m/+ p r)
             n (m/+ q s)]
         (if (m/< x (m// m n))
           (recur (m/inc it) p q m n d y)
           (recur (m/inc it) m n r s d (m/+ y d))))))))

;; sinint / cosint

(defn Si
  "Sine integral, `Si(x) = integral_0^x sin(t)/t dt`.

  `Si` is an entire, odd function that oscillates around and converges to
  `pi/2` as `x` grows without bound (and to `-pi/2` as `x` decreases without
  bound), with its largest overshoot (the Gibbs phenomenon peak, about
  1.18*pi/2) at its first local maximum `x = pi`. Near zero, `Si(x)` behaves
  like `x`, since `sin(t)/t -> 1` as `t -> 0`.

  Parameters:

  - `x` (double): evaluation point, any real number.

  Returns `Si(x)` as a double. `Si(0.0)` is `0.0` and the function saturates
  to the exact limits `+-HALF_PI` for very large `|x|`. Returns `##NaN` for a
  `##NaN` input.

  See also [[si]] (the same function shifted by `-pi/2`), [[Ci]], [[Cin]]
  (the related cosine integrals)."
  ^double [^double x]
  (if (m/nan? x)
    ##NaN
    (let [t (m/* x x)]
      (cond
        (m/<= t 36.0) (m/* x (m// (poly/mevalpoly t 1.00000000000000000000E0 -0.44663998931312457298E-1 0.11209146443112369449E-2
                                                  -0.13276124407928422367E-4 0.85118014179823463879E-7 -0.29989314303147656479E-9
                                                  0.55401971660186204711E-12 -0.42406353433133212926E-15)
                                  (poly/mevalpoly t 1.00000000000000000000E0 0.10891556624243098264E-1 0.59334456769186835896E-4
                                                  0.21231112954641805908E-6 0.54747121846510390750E-9 0.10378561511331814674E-11
                                                  0.13754880327250272679E-14 0.10223981202236205703E-17)))
        (m/<= t 144.0) (let [invt (m// t)
                             p (if (m/neg? x) m/-HALF_PI m/HALF_PI)]
                         (m/- p
                              (m/* (m/cos x)
                                   (m// (poly/mevalpoly invt 0.99999999962173909991E0 0.36451060338631902917E3
                                                        0.44218548041288440874E5 0.22467569405961151887E7
                                                        0.49315316723035561922E8 0.43186795279670283193E9
                                                        0.11847992519956804350E10 0.45573267593795103181E9)
                                        (m/* x (poly/mevalpoly invt 1.00000000000000000000E0 0.36651060273229347594E3
                                                               0.44927569814970692777E5 0.23285354882204041700E7
                                                               0.53117852017228262911E8 0.50335310667241870372E9
                                                               0.16575285015623175410E10 0.11746532837038341076E10))))
                              (m/* (m/sin x) invt
                                   (m// (poly/mevalpoly invt 0.99999999920484901956E0 0.51385504875307321394E3
                                                        0.92293483452013810811E5 0.74071341863359841727E7
                                                        0.28142356162841356551E9 0.49280890357734623984E10
                                                        0.35524762685554302472E11 0.79194271662085049376E11
                                                        0.17942522624413898907E11)
                                        (poly/mevalpoly invt 1.00000000000000000000E0 0.51985504708814870209E3
                                                        0.95292615508125947321E5 0.79215459679762667578E7
                                                        0.31977567790733781460E9 0.62273134702439012114E10
                                                        0.54570971054996441467E11 0.18241750166645704670E12
                                                        0.15407148148861454434E12)))))
        (m/< t ##Inf) (let [invt (m// t)
                            p (if (m/neg? x) m/-HALF_PI m/HALF_PI)]
                        (m/- p
                             (m/* (m// (m/cos x) x)
                                  (m/- 1.0 (m/* invt
                                                (m// (poly/mevalpoly invt 0.19999999999999978257E1 0.22206119380434958727E4
                                                                     0.84749007623988236808E6 0.13959267954823943232E9
                                                                     0.10197205463267975592E11 0.30229865264524075951E12
                                                                     0.27504053804288471142E13 0.21818989704686874983E13)
                                                     (poly/mevalpoly invt 1.00000000000000000000E0 0.11223059690217167788E4
                                                                     0.43685270974851313242E6 0.74654702140658116258E8
                                                                     0.58580034751805687471E10 0.20157980379272098841E12
                                                                     0.26229141857684496445E13 0.87852907334918467516E13)))))
                             (m/* (m/sin x) invt
                                  (m/- 1.0 (m/* invt
                                                (m// (poly/mevalpoly invt 0.59999999999999993089E1 0.96527746044997139158E4
                                                                     0.56077626996568834185E7 0.15022667718927317198E10
                                                                     0.19644271064733088465E12 0.12191368281163225043E14
                                                                     0.31924389898645609533E15 0.25876053010027485934E16
                                                                     0.12754978896268878403E16)
                                                     (poly/mevalpoly invt 1.00000000000000000000E0 0.16287957674166143196E4
                                                                     0.96636303195787870963E6 0.26839734750950667021E9
                                                                     0.37388510548029219241E11 0.26028585666152144496E13
                                                                     0.85134283716950697226E14 0.11304079361627952930E16
                                                                     0.42519841479489798424E16)))))))
        :else (if (neg? x) m/-HALF_PI m/HALF_PI)))))

(defn si
  "Sine integral [[Si]] shifted down by `pi/2`, i.e. `si(x) = Si(x) - pi/2`.

  This shifted form is convenient when working with the large-`x` behavior of
  the sine integral, since `si(x) -> 0` as `x` grows without bound (rather
  than `Si`'s `pi/2`). Unlike [[Si]], `si` is not an odd function (`si(-x) =
  -Si(x) - pi/2`, not `-si(x)`), because the shift by the constant `pi/2`
  breaks the antisymmetry.

  Parameters:

  - `x` (double): evaluation point, any real number.

  Returns `si(x)` as a double. `si(0.0)` is `-HALF_PI`. Converges to `0.0`
  as `x` grows without bound and to `-pi` as `x` decreases without bound.
  Returns `##NaN` for a `##NaN` input.

  See also [[Si]]."
  ^double [^double x] (m/- (Si x) m/HALF_PI))

(def ^:private ^:const ^{:tag 'double} ci-r0 0.616505485620716233797110404100)
(def ^:private ^:const ^{:tag 'double} ci-r1 3.384180422551186426397851146402)
(def ^:private ^:const ^{:tag 'double} ci-r01 0.6162109375)
(def ^:private ^:const ^{:tag 'double} ci-r02 0.29454812071623379711E-3)
(def ^:private ^:const ^{:tag 'double} ci-r11 3.3837890625)
(def ^:private ^:const ^{:tag 'double} ci-r12 0.39136005118642639785E-3)

(defn Ci
  "Cosine integral, `Ci(x) = -integral_x^Inf cos(t)/t dt`, for `x >= 0`.

  Equivalently, `Ci(x) = gamma + ln(x) + integral_0^x (cos(t)-1)/t dt`, where
  `gamma` is the Euler-Mascheroni constant. `Ci` has a logarithmic
  singularity at `0` (`Ci(x) -> -Inf` as `x -> 0+`), oscillates with
  decreasing amplitude for increasing `x`, and converges to `0` as `x` grows
  without bound.

  `Ci` is only defined here for non-negative `x`: for negative `x` the
  natural continuation is genuinely complex (real part `Ci(|x|)`, imaginary
  part `pi`), which a real-valued double-returning function cannot represent,
  so a negative `x` throws rather than silently returning just the real
  part.

  Parameters:

  - `x` (double): evaluation point, must be non-negative.

  Returns `Ci(x)` as a double. Throws an assertion error if `x` is negative.
  Returns `##-Inf` at `x = 0.0`, `0.0` at `x = ##Inf`, and `##NaN` for a
  `##NaN` input.

  See also [[Cin]] (the closely related entire cosine integral, defined for
  all real `x`), [[Si]], [[si]] (the sine integral)."
  ^double [^double x]
  (assert (not (neg? x)) "x must be non-negative")
  (if (m/nan? x)
    ##NaN
    (let [t (m/* x x)]
      (cond
        (m/<= x 3.0) (m/+ (m/log (m// x ci-r0))
                          (m/* (m/- (m/- x ci-r01) ci-r02)
                               (m/+ x ci-r0)
                               (m// (poly/mevalpoly t -0.24607411378767540707E0 0.72113492241301534559E-2
                                                    -0.11867127836204767056E-3 0.90542655466969866243E-6
                                                    -0.34322242412444409037E-8 0.51950683460656886834E-11)
                                    (poly/mevalpoly t 1.00000000000000000000E0 0.12670095552700637845E-1
                                                    0.78168450570724148921E-4 0.29959200177005821677E-6
                                                    0.73191677761328838216E-9 0.94351174530907529061E-12))))
        (m/<= x 6.0) (m/+ (m/log (m// x ci-r1))
                          (m/* (m/- (m/- x ci-r11) ci-r12)
                               (m/+ x ci-r1)
                               (m// (poly/mevalpoly t -0.15684781827145408780E0 0.66253165609605468916E-2
                                                    -0.12822297297864512864E-3 0.12360964097729408891E-5
                                                    -0.66450975112876224532E-8 0.20326936466803159446E-10
                                                    -0.33590883135343844613E-13 0.23686934961435015119E-16)
                                    (poly/mevalpoly t 1.00000000000000000000E0 0.96166044388828741188E-2
                                                    0.45257514591257035006E-4 0.13544922659627723233E-6
                                                    0.27715365686570002081E-9 0.37718676301688932926E-12
                                                    0.27706844497155995398E-15))))
        (m/<= x 12.0) (let [invt (m// t)]
                        (m/- (m/* (m/sin x)
                                  (m// (poly/mevalpoly invt 0.99999999962173909991E0 0.36451060338631902917E3
                                                       0.44218548041288440874E5 0.22467569405961151887E7
                                                       0.49315316723035561922E8 0.43186795279670283193E9
                                                       0.11847992519956804350E10 0.45573267593795103181E9)
                                       (m/* x (poly/mevalpoly invt 1.00000000000000000000E0 0.36651060273229347594E3
                                                              0.44927569814970692777E5 0.23285354882204041700E7
                                                              0.53117852017228262911E8 0.50335310667241870372E9
                                                              0.16575285015623175410E10 0.11746532837038341076E10))))
                             (m/* (m/cos x) invt
                                  (m// (poly/mevalpoly invt 0.99999999920484901956E0 0.51385504875307321394E3
                                                       0.92293483452013810811E5 0.74071341863359841727E7
                                                       0.28142356162841356551E9 0.49280890357734623984E10
                                                       0.35524762685554302472E11 0.79194271662085049376E11
                                                       0.17942522624413898907E11)
                                       (poly/mevalpoly invt 1.00000000000000000000E0 0.51985504708814870209E3
                                                       0.95292615508125947321E5 0.79215459679762667578E7
                                                       0.31977567790733781460E9 0.62273134702439012114E10
                                                       0.54570971054996441467E11 0.18241750166645704670E12
                                                       0.15407148148861454434E12)))))
        (m/< x ##Inf) (let [invt (m// t)]
                        (m/- (m/* (m// (m/sin x) x)
                                  (m/- 1.0 (m/* invt
                                                (m// (poly/mevalpoly invt 0.19999999999999978257E1 0.22206119380434958727E4
                                                                     0.84749007623988236808E6 0.13959267954823943232E9
                                                                     0.10197205463267975592E11 0.30229865264524075951E12
                                                                     0.27504053804288471142E13 0.21818989704686874983E13)
                                                     (poly/mevalpoly invt 1.00000000000000000000E0 0.11223059690217167788E4
                                                                     0.43685270974851313242E6 0.74654702140658116258E8
                                                                     0.58580034751805687471E10 0.20157980379272098841E12
                                                                     0.26229141857684496445E13 0.87852907334918467516E13)))))
                             (m/* (m/cos x) invt
                                  (m/- 1.0 (m/* invt
                                                (m// (poly/mevalpoly invt 0.59999999999999993089E1 0.96527746044997139158E4
                                                                     0.56077626996568834185E7 0.15022667718927317198E10
                                                                     0.19644271064733088465E12 0.12191368281163225043E14
                                                                     0.31924389898645609533E15 0.25876053010027485934E16
                                                                     0.12754978896268878403E16)
                                                     (poly/mevalpoly invt 1.00000000000000000000E0 0.16287957674166143196E4
                                                                     0.96636303195787870963E6 0.26839734750950667021E9
                                                                     0.37388510548029219241E11 0.26028585666152144496E13
                                                                     0.85134283716950697226E14 0.11304079361627952930E16
                                                                     0.42519841479489798424E16)))))))        
        :else 0.0))))

(defn Cin
  "Entire cosine integral, `Cin(x) = integral_0^x (1-cos(t))/t dt`.

  Related to [[Ci]] by `Cin(x) = gamma + ln(|x|) - Ci(|x|)` for `x != 0`
  (where `gamma` is the Euler-Mascheroni constant), but unlike `Ci`, `Cin`
  has no singularity at `0` (the integrand's `0/0` there is removable) and
  is defined for all real `x`, including negative values. `Cin` is an even
  function (`Cin(-x) = Cin(x)`), non-negative, and grows without bound
  (like `ln(|x|)`) as `|x|` grows without bound.

  Parameters:

  - `x` (double): evaluation point, any real number.

  Returns `Cin(x)` as a double. `Cin(0.0)` is `0.0`, and `Cin(x)` diverges to
  `##Inf` as `x` approaches either `##Inf` or `##-Inf`. Returns `##NaN` for a
  `##NaN` input.

  See also [[Ci]], [[Si]], [[si]]."
  ^double [^double x]
  (cond
    (m/nan? x) ##NaN
    (m/zero? x) 0.0
    :else (let [ax (m/abs x)]
            (m/- (m/+ m/GAMMA (m/log ax)) (Ci ax)))))

;; ei

(defn- e1-cf-poly-approx
  [^long n]
  (let [x (poly/ratio-polynomial [0 1])]
    (loop [i n
           p x
           q (poly/ratio-polynomial [1])]
      (if (m/zero? i)
        [p (poly/add (poly/mult x p) q)]
        (let [newp (poly/add (poly/mult x p)
                             (poly/scale q (m/inc i)))
              newq p]
          (recur (m/dec i) (poly/add newp (poly/scale newq i)) newp))))))

(defmacro ^:private e1-cf64
  [x n]
  (let [[p q] (e1-cf-poly-approx n)]
    `(let [num# (poly/mevalpoly ~x ~@(map double (poly/coeffs p)))
           den# (poly/mevalpoly ~x ~@(map double (poly/coeffs q)))]
       (m// num# den#))))

(defn- e1-taylor-coefficients-step
  ^double [^double term ^long k]
  (m// (m/* (m/- term) (m/dec k)) (m/* k k)))

(defn- e1-taylor-coefficients
  [^long n]
  (case (int n)
    0 '()
    1 (list (m/- m/GAMMA))
    (conj (reductions e1-taylor-coefficients-step 1.0 (range 2 (m/inc n))) (m/- m/GAMMA))))

(defmacro ^:private e1-taylor64
  [x n]
  `(m/- (poly/mevalpoly ~x ~@(e1-taylor-coefficients n))
        (m/log ~x)))

(defmacro ^:private ein-taylor64
  [x n]
  `(m/+ (poly/mevalpoly ~x ~@(e1-taylor-coefficients n)) m/GAMMA))

(defn E0
  "Exponential integral E0, `E0(x) = exp(-x)/x`, the `n=0` case of [[En]].

  The simplest member of the generalized exponential integral family
  `En(x) = integral_1^Inf exp(-x*t)/t^n dt`; for `n=0` this integral has a
  closed elementary form. `E0` has a pole at `0` and decays to `0` for large
  positive `x`; for large negative `x` it diverges to `-Inf`.

  Parameters:

  - `x` (double): evaluation point, any real number.

  Returns `E0(x)` as a double. `E0(0.0)` is `##Inf`, `E0(##-Inf)` is `##-Inf`,
  `E0(##Inf)` is `0.0`. Returns `##NaN` for a `##NaN` input.

  See also [[E1]], [[En]] (generalizations to other orders), [[Ei]]."
  ^double [^double x]
  (cond
    (m/zero? x) ##Inf
    (m/neg-inf? x) ##-Inf
    :else (let [e (m/exp (m/- x))]
            (if (m/inf? e)
              ;; `exp(-x)` alone overflows double range for very negative `x`,
              ;; even where the final ratio `exp(-x)/x` would not: recompute
              ;; in log-space to avoid the premature overflow.
              (let [mag (m/exp (m/- (m/- x) (m/log (m/abs x))))]
                (if (m/neg? x) (m/- mag) mag))
              (m// e x)))))

(defn E1
  "Exponential integral E1, `E1(x) = integral_1^Inf exp(-x*t)/t dt`, for
  `x >= 0` (the `n=1` case of [[En]]).

  `E1` has a logarithmic singularity at `0` (`E1(x) -> Inf` as `x -> 0+`)
  and decreases monotonically to `0` as `x` grows without bound. It is only
  defined here for non-negative `x`: for negative `x` the natural
  continuation is complex, which a real-valued double-returning function
  cannot represent, so a negative `x` returns `##NaN` rather than silently
  returning only part of the true (complex) result.

  Parameters:

  - `x` (double): evaluation point, must be non-negative.

  Returns `E1(x)` as a double. `E1(0.0)` is `##Inf`, `E1(##Inf)` is `0.0`.
  Returns `##NaN` for a `##NaN` input or a negative `x`.

  See also [[Ein]] (the closely related entire exponential integral),
  [[E0]], [[En]] (generalizations to other orders), [[Ei]]."
  ^double [^double x]
  (cond
    (m/neg? x) ##NaN
    (m/zero? x) ##Inf
    (m/pos-inf? x) 0.0
    :else (if (m/> x 2.15)
            (let [mult (m/exp (m/- x))]
              (cond
                (m/< x 4.0) (m/* mult
                                 (m// (poly/mevalpoly x 3.600530862438501481559423277418128014798, 28.73031134165011013783185685393062481126, 46.04314409968065653003548224846863877481, 21.47189493062368074985000918414086604187, 2.719957622921613321844755385973197500235, 1.508750885580864752293599048121678438568e-6)
                                      (poly/mevalpoly x 1.0, 18.06743589038646654075831055159865459831, 61.19456872238615922515377354442679999566, 64.81772518730116701207299231777089576859, 24.19034591054828198408354214931112970741, 2.720026796991567940556850921390829046015)))
                (m/< x 10.0) (m/* mult
                                  (m// (poly/mevalpoly x 3.149019890512432117647119992448352099575, 14.77395058090815888166486507990452627136, 14.78214309058953358717796744960600201013, 4.559401130686434886620102186841739864936, 0.4027858394909585103775445204576054721422, 2.302781920509468929446800686773538387432e-9)
                                       (poly/mevalpoly x 1.0, 11.65960373479520446458792926669115987821, 26.20023773503894444535165299118172674849, 18.93899893550582921168134979000319186841, 4.962178168140565906794561070524079194193, 0.4027860481050182948737116109665534358805)))
                (m/< x 20.0) (m/* mult
                                  (m// (poly/mevalpoly x 2.564801308922428705318296668924369454617, 5.482252510134574167659359513298970499768, 2.379528224853089764405551768869103662657, 0.2523431276282591480166881146593592650031, 1.444719769329975045925053905197199934930e-9, -8.977332061106791470409502623202677045468e-12)
                                       (poly/mevalpoly x 1.0, 6.421214405318108272004472721910023284626, 7.609584052290707052877352911548283916247, 2.631866613851763364839413823973711355890, 0.2523432345539146248733539983749722854603)))
                (m/< x 200.0) (m/* mult (e1-cf64 x 8))
                :else (m/* mult (e1-cf64 x 4))))
            (cond
              (m/> x 0.6) (e1-taylor64 x 37)
              (m/> x 0.053) (e1-taylor64 x 15)
              (m/> x 4.4e-3) (e1-taylor64 x 8)
              :else (e1-taylor64 x 4)))))

(defn Ein
  "Entire exponential integral, `Ein(x) = integral_0^x (1-exp(-t))/t dt`.

  Related to [[E1]] by `Ein(x) = E1(x) + ln(x) + gamma` for `x > 0` (where
  `gamma` is the Euler-Mascheroni constant), but unlike `E1`, `Ein` has no
  singularity at `0` (the integrand's `0/0` there is removable) and is, in
  principle, defined for every real (and complex) `x`. `Ein` grows without
  bound as `x` grows without bound, and approaches `-##Inf` as `x` decreases
  without bound (mirroring the corresponding limits of `E1 + ln(x)`).

  Parameters:

  - `x` (double): evaluation point. Positive `x` is fully supported; negative
    `x` is only supported down to `-2.15`, a known limitation of the current
    implementation (below that, the result is `##NaN`).

  Returns `Ein(x)` as a double. `Ein(0.0)` is `0.0`. Returns `##NaN` for a
  `##NaN` input or for `x < -2.15`.

  See also [[E1]], [[Cin]] (the analogous entire form of the cosine
  integral), [[Ei]]."
  ^double [^double x]
  (cond
    (m/nan? x) ##NaN
    (m/zero? x) 0.0
    (m/neg? x) (let [ax (m/abs x)]
                 (if (m/<= ax 2.15)
                   (cond
                     (m/> ax 0.6) (ein-taylor64 x 37)
                     (m/> ax 0.053) (ein-taylor64 x 15)
                     (m/> ax 4.4e-3) (ein-taylor64 x 8)
                     :else (ein-taylor64 x 4))
                   ##NaN))
    :else (m/+ (E1 x) (m/log x) m/GAMMA)))

(defn- en-safe-expfact
  ^double [^long v ^double x]
  (if (m/< v 100)
    (let [-x (m/- x)]
      (loop [i (long 1)
             powerterm 1.0]
        (if (m/> i v)
          powerterm
          (recur (m/inc i) (m/* powerterm (m// -x i))))))
    (let [sgn (if (m/not-pos? x) 1.0 (if (m/odd? v) -1.0 1.0))]
      (m/* sgn (m/exp (m/- (m/* v (m/log (m/abs x))) (log-gamma (m/inc v))))))))

(defn- en-expand-origin-posint
  ^double [^long v ^double x]
  (let [gamma-term (m/* (en-safe-expfact (m/long-dec v) x)
                        (m/- (digamma v) (m/log x)))
        sum-term (if (m/one? v) 0.0 (m// 1.0 (m/- 1.0 v)))
        eps (m/* 10.0 (m/ulp sum-term))
        v- (m/dec v)
        -x (m/- x)]
    (loop [k (long 1)
           frac 1.0
           sum-term sum-term]
      (let [nfrac (m/* frac (m// -x k))]
        (if (m/not== k v-)
          (let [nsum-term (m/+ sum-term (m// nfrac (m/- k v-)))]
            (if (or (m/> k 1000) (m/delta-eq sum-term nsum-term eps))
              (m/- gamma-term nsum-term)
              (recur (m/inc k) nfrac nsum-term)))
          (recur (m/inc k) nfrac sum-term))))))

(def ^{:private true :const true :tag 'double} SQPI 9.869604401089358)
(def ^{:private true :const true :tag 'double} SQPI2 19.739208802178716)
(def ^{:private true :const true :tag 'double} SQPI10 98.69604401089359)
(def ^{:private true :const true :tag 'double} SQSQPI7 681.863637238017)

(defn- en-expand-origin-general
  ^double [^double v ^double x]
  (let [omv (m/- 1.0 v)
        invomv (m// omv)
        aomv (m/abs omv)
        -x (m/- x)        
        gamma-term (m/* (gamma omv) (m/pow x (m/dec v)))
        ^Vec2 bs (loop [k (long 1)
                        frac 1.0
                        blowup (if (m/< aomv 0.5) invomv 0.0)
                        sum-term (if (m/< aomv 0.5) 0.0 invomv)]                   
                   (let [nfrac (m/* frac (m// -x k))
                         den (m/+ k omv)]
                     (if (m/< (m/abs den) 0.5)
                       (recur (m/inc k) nfrac (m/+ blowup (m// nfrac den)) sum-term)
                       (let [nsum-term (m/+ sum-term (m// nfrac den))]
                         (if (or (m/< (m/abs (m/- nsum-term sum-term))
                                      (m/* m/MACHINE-EPSILON10 (m/abs sum-term)))
                                 (m/== k 1000))
                           (Vec2. blowup sum-term)
                           (recur (m/inc k) nfrac blowup nsum-term))))))]
    (if (m/< (m/abs (m/- gamma-term (.x bs))) (m/* 1.0e-3 (m/abs (.x bs))))
      (let [delta (m/- (m/round v) v)
            delta2 (m/* delta delta)
            delta3 (m/* delta2 delta)
            delta4 (m/* delta3 delta)
            n (m/long-dec (m/round v))
            n+ (m/inc n)
            logx (m/log x)
            logx2 (m/* logx logx)
            logx3 (m/* logx2 logx)
            logx4 (m/* logx3 logx)
            logx5 (m/* logx4 logx)
            series1 (m/- (m/- logx)
                         (m// (m/* logx2 delta) 2.0)
                         (m// (m/* logx3 delta2) 6.0)
                         (m// (m/* logx4 delta3) 24.0)
                         (m// (m/* logx5 delta4) 120.0))
            psi0 (polygamma 0 n+)
            psi02 (m/* psi0 psi0)
            psi03 (m/* psi02 psi0)
            psi04 (m/* psi03 psi0)
            psi05 (m/* psi04 psi0)
            psi1 (polygamma 1 n+)
            psi13 (m/* 3.0 psi1)
            psi2 (polygamma 2 n+)
            psi3 (polygamma 3 n+)
            psi4 (polygamma 4 n+)
            series2 (m/+ psi0
                         (m// (m/* delta (m/+ (m/* 3.0 psi02) SQPI (m/* -3.0 psi1))) 6.0)
                         (m// (m/* delta2 (m/+ psi03 psi2 (m/* psi0 (m/- SQPI psi13)))) 6.0)
                         (m// (m/* delta3 (m/+ SQSQPI7 (m/* -15.0 psi3)
                                               (m/* 15.0 (m/+ psi04 (m/* 4.0 psi0 psi2)
                                                              (m/* 2.0 psi02 (m/- SQPI psi13))
                                                              (m/* psi1 (m/- psi13 SQPI2)))))) 360.0)
                         (m// (m/* delta4 (m/+ (m/* 3.0 psi05) (m/* -30.0 psi1 psi2)
                                               (m/* SQPI10 psi2) (m/* 3.0 psi4)
                                               (m/* psi03 (m/- SQPI10 (m/* 30.0 psi1)))
                                               (m/* 30.0 psi02 psi2)
                                               (m/* psi0 (m/+ (m/* 45.0 psi1 psi1)
                                                              (m/* -3.0 SQPI10 psi1)
                                                              (m/* -15.0 psi3)
                                                              SQSQPI7)))) 360.0))]
        (m/- (m/* (m/+ series1 series2)
                  (en-safe-expfact n x)
                  (m/pow x (m/- v n 1.0)))
             (.y bs)))
      (m/- gamma-term (.x bs) (.y bs)))))

(defn- en-safe-gamma-term
  ^double [^double v ^double x]
  (let [v1 (m/- 1.0 v)
        g (gamma v1)]
    (m/* (m/sgn g) (m/exp (m/+ (m/* (m/dec v) (m/log x)) (m/log (m/abs g)))))))

(def ^{:private true :const true :tag 'double} SQRTMAXDOUBLE 1.3407807929942596E154)

(defn- en-cf-gamma
  ^Vec2 [^double v ^double x]
  (loop [i (long 1)
         A (m/- 1.0 v)
         B 1.0
         Ap 1.0
         Bp 0.0]
    (let [i+ (m/inc i)
          a (if (m/even? i)
              (m/* x (m// i 2))
              (m/* -1.0 x (m/- (m// i+ 2) v)))
          b (m/- i+ v)
          nA (m/+ (m/* b A) (m/* a Ap))
          nB (m/+ (m/* b B) (m/* a Bp))
          q (m/* A nB)]
      (if (or (m/< (m/abs (m/- q (m/* nA B))) (m/* m/MACHINE-EPSILON10 (m/abs q)))
              (m/== i 1000))
        (Vec2. (en-safe-gamma-term v x)
               (m// (m/- B) A))
        (if (m/> (m/abs nA) SQRTMAXDOUBLE)
          (recur (m/inc i)
                 (m// nA SQRTMAXDOUBLE) (m// nB SQRTMAXDOUBLE)
                 (m// A SQRTMAXDOUBLE) (m// B SQRTMAXDOUBLE))
          (recur (m/inc i) nA nB A B))))))

(defn- en-cf-no-gamma
  ^double [^double v ^double x]
  (let [B (m/+ v x)
        eps (m/* 10.0 (m/ulp B))]
    (loop [i (long 2)
           A 1.0
           B B
           Ap 1.0
           Bp x]
      (let [i- (m/dec i)
            nA (m/+ (m/* x A) (m/* i- Ap))
            nB (m/+ (m/* x B) (m/* i- Bp))]
        (if (or (m/inf? nA) (m/inf? nB))
          (m// nA nB)
          (let [v+ (m/+ v i-)
                nAp nA
                nA (m/+ nA (m/* v+ A))
                nBp nB
                nB (m/+ nB (m/* v+ B))]
            (if (or (and (m/> i 4)
                         (m/< (m/abs (m/- (m/* nAp nB) (m/* nA nBp))) (m/* eps (m/abs (m/* nB nBp)))))
                    (m/== i 1000))
              (m// nA nB)
              (if (m/> (m/abs nA) SQRTMAXDOUBLE)
                (recur (m/inc i)
                       (m// nA SQRTMAXDOUBLE) (m// nB SQRTMAXDOUBLE)
                       (m// nAp SQRTMAXDOUBLE) (m// nBp SQRTMAXDOUBLE))
                (recur (m/inc i) nA nB nAp nBp)))))))))

(defn- en-cf
  ^Vec2 [^double v ^double x]
  (if (m/pos? (m/- 1.0 v))
    (let [^Vec2 gcf (en-cf-gamma v x)
          ag (m/abs (.x gcf))
          acf (m/abs (.y gcf))]
      (if (and (m/valid-double? ag) (m/> ag 1.0) (m/> ag acf))
        gcf
        (Vec2. 0.0 (en-cf-no-gamma v x))))
    (Vec2. 0.0 (en-cf-no-gamma v x))))

(defn En
  "Generalized exponential integral, `En(x) = integral_1^Inf exp(-x*t)/t^n
  dt`, for any real order `n`.

  Includes [[E0]] and [[E1]] as the `n=0` and `n=1` special cases. For
  positive `x`, `En` decreases monotonically to `0` as `x` grows without
  bound. As `x -> 0+`, `En` diverges to `##Inf` when `n <= 1` (the defining
  integral itself diverges there), but converges to the finite value
  `1/(n-1)` when `n > 1`.

  Negative `x` is supported only when `n` is a non-positive integer (`n <=
  0`): in that case `En` reduces to an elementary function of `x` and
  `exp(-x)` with no branch cut, so it stays real. For any other order
  (fractional, or a positive integer) with negative `x`, the natural
  continuation is complex, which a real-valued double-returning function
  cannot represent, so the result is `##NaN`.

  Parameters:

  - `n` (double): the order, any real number (integer or fractional).
  - `x` (double): evaluation point.

  Returns `En(x)` as a double. At `x = 0.0`: `##Inf` if `n <= 1`, `1/(n-1)`
  otherwise. Returns `##NaN` for a `##NaN` `n` or `x`, or for a negative `x`
  outside the non-positive-integer-order case described above.

  See also [[E0]], [[E1]], [[Ei]]."
  ^double [^double n ^double x]
  (cond
    (m/zero? n) (E0 x)
    (m/one? n) (E1 x)
    (and (m/zero? x) (m/< n 1.0)) ##Inf
    (m/zero? x) (m// 1.0 (m/dec n))
    (or (m/nan? n) (m/nan? x)
        (if (m/integer? n)
          (and (m/neg? x) (m/pos? n))
          (m/neg? x))) ##NaN
    (m/> x 745.0) 0.0
    (m/< (m/sq x) 9.0) (if (and (m/integer? n) (m/pos? n))
                         (en-expand-origin-posint (unchecked-long n) x)
                         (en-expand-origin-general n x))
    :else (let [^Vec2 gcf (if (m/pos? x)
                            (en-cf n x)
                            (Vec2. 0.0 (en-cf-no-gamma n x)))
                cf (.y gcf)
                e (m/exp (m/- x))
                em (if (or (m/inf? e) (m/zero? e))
                     (m/* (m/sgn cf) (m/exp (m/- (m/log (m/abs cf)) x)))
                     (m/* e cf))]
            (m/+ em (.x gcf)))))

(defmacro ^:private ei-taylor64
  [x n]
  (let [coeffs (e1-taylor-coefficients n)]
    `(m/+ (poly/mevalpoly ~x ~@(map-indexed (fn [^long i ^double c]
                                              (m/* (m/- c) (m/fpow -1.0 i))) coeffs))
          (m/log ~x))))


(defn Ei
  "Exponential integral, `Ei(x) = -integral_(-x)^Inf exp(-t)/t dt` (Cauchy
  principal value), for `x != 0`.

  Related to [[E1]] by `Ei(x) = -E1(-x)` for `x < 0`. `Ei` has a logarithmic
  singularity at `0` (`-Inf` from the left, `+Inf` from the right), a single
  real zero near `x = 0.3725`, and grows without bound (like `exp(x)/x`) as
  `x` grows without bound; as `x` decreases without bound, `Ei(x)`
  approaches `0` from below.

  Parameters:

  - `x` (double): evaluation point, any real number except `0`.

  Returns `Ei(x)` as a double. `Ei(0.0)` is `##-Inf`, `Ei(##Inf)` is `##Inf`,
  `Ei(##-Inf)` is `-0.0`. Returns `##NaN` for a `##NaN` input.

  See also [[E1]], [[E0]], [[En]], [[li]] (related by `li(x) = Ei(ln x)`)."
  ^double [^double x]
  (cond
    (m/neg? x) (m/- (E1 (m/- x)))
    (m/zero? x) ##-Inf
    (m/pos-inf? x) ##Inf
    (m/neg-inf? x) -0.0
    :else (if (m/> x 2.15)
            (cond
              (m/< x 4.0) (m// (poly/mevalpoly x -2.43791466332154621,3.09402100064798205,9.35202477109609955,0.152659977028953397,0.0157273683896079142,0.0345566671015011426,-0.000421531433157416203)
                               (poly/mevalpoly x 1.0,4.28055563991564399,0.537599625698465573,-0.511064414527643313,0.0867748262262250088,-0.00623913330836521800,0.000172066498182538260))
              (m/< x 10.0) (m/* (m/exp x)
                                (m// (poly/mevalpoly x -1.58447018083420958,4.71806833998906997,-0.587691572500210206,0.125012472861504555,-0.00178055441724967428,0.000633648975971195928,0.0000147213934578379204,2.12391754244415544e-6)
                                     (poly/mevalpoly x 1.0,1.93297600031287800,0.660790440069106542,0.198322636197663277,0.0272447293513279631,0.00399501571688512611,0.000362510989191199243,0.0000182930089855534336,2.06800780072894204e-6)))
              (m/< x 20.0) (m/* (m/exp x)
                                (m// (poly/mevalpoly x -1.77183291754640123,0.795659966861260409,-0.221223333413388642,0.0328877243243796815,-0.00331846947191676458,0.000180945604349930285,-5.97641401680304362e-6,2.42151808626299747e-11)
                                     (poly/mevalpoly x 1.0,-2.10926998628216150,0.933357955421497965,-0.245433884954174394,0.0356954809772243699,-0.00348034743685382360,0.000186615220328647350,-5.97232581353392099e-6)))
              :else (let [xinv (m// x)]
                      (if (m/< x 200.0)
                        (m/* (m/exp x)
                             (m// (poly/mevalpoly xinv -5.29842699621003563e-14, +1.00000000004732488, -60.4361334939888359, +1327.83891720487710, -6810.63668974273961, -177755.383525765400,+3.00773484037048848e6, -1.53642380695372707e7, +2.08174653368702692e7)
                                  (poly/mevalpoly xinv  1.0, -61.4361334756161381, +1387.27504658395142, -8081.03888544858393, -172104.333927401741, +3.18903576285551101e6, -1.81873890267574206e7, +3.37312131843327704e7, -1.22198734384213631e7)))
                        (let [e (m/exp x)
                              polyval (poly/mevalpoly xinv 1,1,2,6,24,120,720,5040)]
                          (if (m/inf? e)
                            ;; `exp(x)` alone overflows double range before the
                            ;; true `exp(x)/x * polyval` product does: recompute
                            ;; in log-space to avoid the premature overflow.
                            (m/exp (m/+ (m/- x (m/log x)) (m/log polyval)))
                            (m/* e xinv polyval))))))
            (cond
              (m/> x 0.6) (ei-taylor64 x 37)
              (m/> x 0.053) (ei-taylor64 x 15)
              (m/> x 4.4e-3) (ei-taylor64 x 8)
              :else (ei-taylor64 x 4)))))

(defn li
  "Logarithmic integral, `li(x) = integral_0^x dt/ln(t)`, for `x > 0`.

  Related to [[Ei]] by `li(x) = Ei(ln(x))`. `li` has a singularity at `x =
  1` (`-Inf` from the left, `+Inf` from the right, since `ln(1) = 0`), a
  single real zero away from `0` near `x = 1.4514` (the Ramanujan-Soldner
  constant), and grows without bound (roughly like `x/ln(x)`) as `x` grows
  without bound. By the prime number theorem, `li` is asymptotic to the
  prime-counting function `pi(x)` (the count of primes up to `x`).

  Parameters:

  - `x` (double): evaluation point, must be positive.

  Returns `li(x)` as a double. `li(0.0)` is `-0.0`, `li(1.0)` is `##-Inf`.
  Returns `##NaN` for a `##NaN` input or a negative `x` (`ln(x)` is not
  real there).

  See also [[Ei]], [[Li]] (the offset variant used in number theory)."
  ^double [^double x]
  (Ei (m/ln x)))

(defn Li
  "Offset logarithmic integral, `Li(x) = integral_2^x dt/ln(t)`, for `x > 0`.

  Related to [[li]] by `Li(x) = li(x) - li(2)`, so that `Li(2.0)` is `0.0`
  (unlike `li`, whose corresponding reference point is `0`). Shares `li`'s
  singularity at `x = 1` and its asymptotic growth for large `x`; `Li` is
  the form more commonly used in number theory as an estimate of the
  prime-counting function `pi(x)`.

  Parameters:

  - `x` (double): evaluation point, must be positive.

  Returns `Li(x)` as a double. `Li(2.0)` is `0.0`, `Li(1.0)` is `##-Inf`.
  Returns `##NaN` for a `##NaN` input or a negative `x`.

  See also [[li]], [[Ei]]."
  ^double [^double x]
  (m/- (Ei (m/ln x)) 1.04516378011749278484))

;;

(defn upper-incomplete-gamma
  "Upper incomplete gamma function"
  ^double [^double s ^double x]
  (if (pos? s)
    (m/exp (m/+ (m/log (. Gamma (regularizedGammaQ s x))) (. Gamma (logGamma s))))
    (m/* (m/pow x s) (En (m/- 1.0 s) x))))

(defn lower-incomplete-gamma
  "Lower incomplete gamma function"
  ^double [^double s ^double x] (m/- (gamma s) (upper-incomplete-gamma s x)))

(defn regularized-gamma-p
  "Regularized gamma P(a,x)"
  ^double [^double a ^double x]
  (if (m/pos? a)
    (. Gamma (regularizedGammaP a x))
    (m// (lower-incomplete-gamma a x) (gamma a))))

(defn regularized-gamma-q
  "Regularized gamma Q(a,x)"
  ^double [^double a ^double x]
  (if (m/pos? a)
    (. Gamma (regularizedGammaQ a x))
    (m// (upper-incomplete-gamma a x) (gamma a))))

;; Airy

(defn airy-Ai
  "Airy Ai function"
  ^double [^double x]
  (cond
    (m/pos-inf? x) 0.0
    (m/not-neg? x) (airy/ai-pos-args x)
    (m/> x -1.0e8) (airy/ai-neg-args x)
    :else ##NaN))

(defn airy-Ai'
  "First derivative of the Airy Ai function"
  ^double [^double x]
  (cond
    (m/pos-inf? x) 0.0
    (m/not-neg? x) (airy/aip-pos-args x)
    (m/> x -1.0e8) (airy/aip-neg-args x)
    :else ##NaN))

(defn airy-Bi
  "Airy Bi function"
  ^double [^double x]
  (cond
    (m/pos-inf? x) ##Inf
    (m/not-neg? x) (airy/bi-pos-args x)
    (m/> x -1.0e8) (airy/bi-neg-args x)
    :else ##NaN))

(defn airy-Bi'
  "First derivative of the Airy Bi function"
  ^double [^double x]
  (cond
    (m/pos-inf? x) ##Inf
    (m/not-neg? x) (airy/bip-pos-args x)
    (m/> x -1.0e8) (airy/bip-neg-args x)
    :else ##NaN))

;;

(defn harmonic-number
  "Harmonic number H_n or generalized harmonic number H_n,m"
  (^double [^double n]
   (if (m/zero? n)
     0.0
     (m/+ (digamma (m/inc n)) m/GAMMA)))
  (^double [^double n ^double m]
   (cond
     (m/zero? m) n 
     (m/one? m) (harmonic-number n)
     :else (m/- (zeta m) (zeta m (m/inc n))))))

;;

(def ^{:private true :const true :tag 'double} -INVE -0.36787944117144232159552)

;;  approximated by the quadratic-rate recursive formula of R. Iacono and J.P. Boyd

(defn- lambert-W-recursive
  ^double [^double w0 ^double x]
  (loop [i (long 0)
         w w0]
    (let [nw (m/* (m// w (m/inc w))
                  (m/inc (m/log (m// x w))))]
      (if (or (m/delta-eq w nw m/MACHINE-EPSILON m/MACHINE-EPSILON)
              (m/== i 1000))
        nw
        (recur (m/inc i) nw)))))

(defn lambert-W
  "Lambert W_0 function. W(xe^x)=x for x>=-1.0."
  ^double [^double x]
  (cond
    (m/< x -INVE) ##NaN
    (m/== x -INVE) -1.0
    (m/zero? x) 0.0
    (m/one? x) 0.567143290409783873
    (m/== m/E x ) 1.0
    (m/neg? x) (let [ex (m/* m/E x)
                     f (m/inc (m/sqrt (m/inc ex)))]
                 (lambert-W-recursive (m// (m/* ex (m/log f)) (m/+ ex f)) x))
    (m/< x m/E) (lambert-W-recursive (m// x m/E) x)
    :else (let [lx (m/log x)]
            (lambert-W-recursive (m/- lx (m/log lx)) x))))

(defn lambert-W-1
  "Lambert W_1 function. W_1(xe^x)=x for x<=-1.0."
  ^double [^double x]
  (cond
    (or (m/< x -INVE) (m/pos? x)) ##NaN
    (m/== -INVE x ) -1.0
    (m/zero? x) ##-Inf
    (m/< x -0.25) (lambert-W-recursive (m/- -1.0 (m/* m/SQRT2 (m/sqrt (m/inc (m/* m/E x))))) x)
    :else (let [lx (m/log (m/- x))]
            (lambert-W-recursive (m/- lx (m/log (m/- lx))) x))))


;;

(defn kummers-M
  "Kummer's (confluent hypergeometric, 1F1) function for real arguments."
  ^double [^double a ^double b ^double x]
  (cond
    ;; x=0 and a=0 always win, even over an otherwise-genuine pole, matching
    ;; this function's own documented behavior
    (m/near-zero? x m/MACHINE-EPSILON) 1.0
    (m/near-zero? a (m/ulp a)) 1.0
    ;; a genuine, unavoidable pole: b a negative integer, and a is not a
    ;; non-positive integer capable of terminating the series at or before
    ;; reaching it (a need not itself be an integer for this to apply --
    ;; only a non-positive INTEGER a can ever terminate the series)
    (and (m/neg? b)
         (m/integer? b)
         (not (and (m/not-pos? a) (m/integer? a) (m/>= a b)))) ##NaN
    (m/zero? b) (m/copy-sign ##Inf (m/* a x))
    (and (m/neg? a) (m/integer? a))
    ;; a is a non-positive integer: the series always terminates to an exact
    ;; finite polynomial in x (numerator vanishes identically after |a|+1
    ;; terms), for any b and any x -- evaluated directly here instead of via
    ;; any convergence-based approximation, both for precision and to
    ;; correctly resolve the case a = b on a non-positive integer, a
    ;; removable 0/0 coincidence in the series (not exp(x)).
    (let [nmax (unchecked-long (m/- a))]
      (loop [n (long 1) term 1.0 sum 1.0]
        (if (m/> n nmax)
          sum
          (let [nterm (m// (m/* term (m/+ a (m/dec n)) x) (m/* (m/+ b (m/dec n)) n))]
            (recur (m/inc n) nterm (m/+ sum nterm))))))
    (m/== a b) (m/exp x)
    (and (m/one? a) (m/== b 2.0)) (let [hx (m/* 0.5 x)]
                                    (m/* (m// (m/exp hx) hx) (m/sinh hx)))
    (m/pos? x) (loop [i (long 1)
                      s0 1.0
                      s1 (m/inc (m// (m/* a x) b))]
                 (cond
                   (m/inf? s1) s1
                   (or (and (m/valid-double? s0) (m/valid-double? s1)
                            (m/delta-eq s0 s1 m/MACHINE-EPSILON m/MACHINE-EPSILON))
                       (m/== i 1000000)) s1
                   :else (let [rj (m// (m/* (m/+ a i) x) (m/* (m/+ b i) (m/inc i)))]
                           (recur (m/inc i) s1 (m/+ s1 (m/* (m/- s1 s0) rj))))))
    :else (hg/weniger-1F1 a b x)))

(defn whittaker-M
  "Whittaker's M function.

  A standard solution of Whittaker's differential equation, expressed via Kummer's confluent hypergeometric function [[kummers-M]]: `M(kappa,mu,x) = exp(-x/2) x^(mu+1/2) kummers-M(mu-kappa+1/2, 2mu+1, x)`.

  Parameters:

  - `kappa` (double): the first parameter.
  - `mu` (double): the second parameter.
  - `x` (double): the argument, restricted to `x >= 0.0` (for `x < 0.0`, `x^(mu+1/2)` generally leaves the real line, so `##NaN` is returned).

  As `x` approaches `0.0` from above, the result approaches `0.0` when `mu > -0.5`, diverges to `##Inf` when `mu < -0.5`, and approaches `1.0` when `mu = -0.5` exactly (`x^(mu+1/2)` reduces to the constant `1` there). Elsewhere, inherits [[kummers-M]]'s domain and pole structure through the transformed parameters `mu-kappa+1/2` and `2mu+1`.

  See also [[kummers-M]], [[whittaker-W]], [[tricomis-U]]."
  ^double [^double kappa ^double mu ^double x]
  (let [mu+05 (m/+ 0.5 mu)]
    (cond
      (m/neg? x) ##NaN
      (and (m/zero? x) (m/zero? mu+05)) (kummers-M (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x)
      :else (let [z (m/exp (m/* 0.5 (m/+ (m/* -0.5 x) (m/* mu+05 (m/log x)))))]
              (m/* z (kummers-M (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x) z)))))

(defn hypergeometric-0F0
  "Hypergeometric ₀F₀ function.

  The degenerate generalized hypergeometric series with no numerator and no denominator parameters, `sum_n x^n / n!`, identically equal to `exp(x)` for every real `x`.

  Parameters:

  - `x` (double): the argument.

  Returns `exp(x)`. Well-defined and finite for every finite `x`; returns `##Inf` for `x = ##Inf`, `0.0` for `x = ##-Inf`, and `##NaN` for `##NaN`.

  See also [[hypergeometric-1F0]], [[hypergeometric-0F1]], [[hypergeometric-pFq]]."
  ^double [^double x] (m/exp x))

(defn hypergeometric-1F0
  "Hypergeometric ₁F₀ function.

  The generalized hypergeometric series with one numerator parameter and no denominator parameters, equal on the real line to the binomial series `(1-x)^(-a)`. For `a` a non-positive integer the series terminates, giving a polynomial in `x` that is real and finite for every `x`.

  Parameters:

  - `a` (double): the numerator parameter.
  - `x` (double): the argument.

  Returns `(1-x)^(-a)`. Diverges to `##Inf` (or `##-Inf`, depending on sign) as `x` approaches `1.0` from below when `a` is positive and not a non-positive integer. For `x >= 1.0` the result stays on the real line only when `a` is an integer (of any sign); otherwise the true value is complex and `##NaN` is returned. Returns `1.0` for `a = 0.0` regardless of `x`.

  See also [[hypergeometric-0F0]], [[hypergeometric-0F1]], [[hypergeometric-2F1]], [[hypergeometric-pFq]]."
  ^double [^double a ^double x]
  (m/pow (m/- 1.0 x) (m/- a)))

(defn hypergeometric-0F1
  "Confluent hypergeometric ₀F₁ limit function.

  The generalized hypergeometric series with no numerator parameters and one denominator parameter, `sum_n x^n / ((a)_n n!)` (using the Pochhammer symbol `(a)_n`), closely related to the Bessel functions: `0F1(;a;x) = Gamma(a) x^((1-a)/2) I_(a-1)(2 sqrt(x))` for `x > 0`, with `I` replaced by the ordinary Bessel `J` for `x < 0`.

  Parameters:

  - `a` (double): the denominator parameter.
  - `x` (double): the argument.

  Returns `1.0` for `x = 0.0`, for every `a`. Has poles at `a` equal to a non-positive integer (0, -1, -2, and so on) for any `x != 0.0`, where `##NaN` is returned; elsewhere the result is real and finite.

  See also [[hypergeometric-0F0]], [[hypergeometric-1F0]], [[hypergeometric-1F1]], [[bessel-I]], [[bessel-J]], [[hypergeometric-pFq]]."
  ^double [^double a ^double x]
  (cond
    (m/zero? x) 1.0
    (m/neg? x) (let [xx (m/* 2.0 (m/sqrt (m/abs x)))
                     a- (m/dec a)]
                 (m// (m/* (bessel-J a- xx) (gamma a))
                      (m/pow (m/* 0.5 xx ) a-)))
    :else (let [xx (m/* 2.0 (m/sqrt x))
                a- (m/dec a)]
            (m// (m/* (bessel-I a- xx) (gamma a))
                 (m/pow (m/* 0.5 xx ) a-)))))

(defn hypergeometric-1F1
  "Confluent hypergeometric ₁F₁ function of the first kind, also known as Kummer's function M.

  The generalized hypergeometric series with one numerator and one denominator parameter, `sum_n (a)_n x^n / ((b)_n n!)` (using the Pochhammer symbol `(a)_n`). Alias for [[kummers-M]].

  Parameters:

  - `a` (double): the numerator parameter.
  - `b` (double): the denominator parameter.
  - `x` (double): the argument.

  Returns `1.0` for `x = 0.0` regardless of `a` and `b`, and also `1.0` whenever `a = 0.0` regardless of `x` and `b` (the series then has only its constant term).

  Has a genuine pole (returns `##NaN`) when `b` is a non-positive integer, unless `a` is also a non-positive integer no more negative than `b`, in which case the series terminates to a finite polynomial in `x` before ever reaching the pole. When `a` equals `b` on a negative integer, numerator and denominator vanish together at the same term, a removable coincidence whose value is the corresponding truncated exponential series, not `exp(x)`; `exp(x)` remains the correct result for every other case of `a = b` (zero, positive, or non-integer). For `b = 0.0` with `a != 0.0`, diverges to a signed `##Inf` or `##-Inf` depending on the sign of `a * x`. Returns a signed `##Inf` (rather than `##NaN`) whenever the true magnitude exceeds the double-precision range.

  See also [[kummers-M]], [[whittaker-M]], [[tricomis-U]], [[hypergeometric-0F1]], [[hypergeometric-2F1]], [[hypergeometric-pFq]]."
  ^double [^double a ^double b ^double x]
  (kummers-M a b x))

(defn hypergeometric-0F2
  "Generalized hypergeometric ₀F₂ function.

  The generalized hypergeometric series with no numerator parameters and two denominator parameters, `sum_n x^n / ((a)_n (b)_n n!)` (using the Pochhammer symbol `(a)_n`).

  Parameters:

  - `a`, `b` (double): the two denominator parameters.
  - `x` (double): the argument.

  Returns `1.0` for `x = 0.0`, for every `a` and `b`. Has a genuine pole (returns `##NaN`) whenever `a` or `b` is a non-positive integer and `x != 0.0`, since both parameters appear only in the denominator and so cannot terminate the series the way a non-positive-integer numerator parameter would; elsewhere the result is real and finite.

  See also [[hypergeometric-0F1]], [[hypergeometric-1F1]], [[hypergeometric-2F0]], [[hypergeometric-pFq]]."
  ^double [^double a ^double b ^double x]
  (cond
    (m/zero? x) 1.0
    (or (and (m/not-pos? a) (m/integer? a))
        (and (m/not-pos? b) (m/integer? b))) ##NaN
    (m/pos? x) (hg/maclaurin-0F2 a b x)
    :else (hg/weniger-0F2 a b x)))

(defn hypergeometric-2F0
  "Generalized hypergeometric ₂F₀ function.

  The generalized hypergeometric series with two numerator parameters and no denominator parameters, `sum_n (a)_n (b)_n x^n / n!` (using the Pochhammer symbol `(a)_n`). Unlike [[hypergeometric-0F1]] or [[hypergeometric-0F2]], this series diverges for every `x != 0.0` (its radius of convergence is `0`), except when `a` or `b` is a non-positive integer, in which case it terminates to a finite polynomial; otherwise it is understood here via resummation of the divergent series, the same technique underlying [[tricomis-U]].

  Parameters:

  - `a`, `b` (double): the two numerator parameters.
  - `x` (double): the argument.

  Returns `1.0` for `x = 0.0`. When `a` or `b` is a non-positive integer, the result is a finite polynomial in `x`, real for every `x`. Otherwise, the resummed value matches the standard real result for `x <= 0.0`; for `x > 0.0` the resummation of a divergent series is inherently branch-dependent (other conventions, such as Borel summation, can give a genuinely complex result there instead), so the real value returned should not be assumed to match every other convention.

  See also [[hypergeometric-1F1]], [[hypergeometric-2F1]], [[tricomis-U]], [[hypergeometric-pFq]]."
  ^double [^double a ^double b ^double x]
  (let [a-term? (and (m/not-pos? a) (m/integer? a))
        b-term? (and (m/not-pos? b) (m/integer? b))]
    (if (or a-term? b-term?)
      ;; a or b is a non-positive integer: the series always terminates to
      ;; an exact finite polynomial in x (both parameters are numerator-only
      ;; here, so no pole is possible), evaluated directly for both
      ;; precision and robustness (the general resummation below is prone
      ;; to isolated NaNs at specific coincidental (a, b, x) here).
      (let [nmax (unchecked-long (cond
                                    (and a-term? b-term?) (m/min (m/- a) (m/- b))
                                    a-term? (m/- a)
                                    :else (m/- b)))]
        (loop [n (long 1) term 1.0 sum 1.0]
          (if (m/> n nmax)
            sum
            (let [nterm (m// (m/* term (m/+ a (m/dec n)) (m/+ b (m/dec n)) x) n)]
              (recur (m/inc n) nterm (m/+ sum nterm))))))
      (hg/weniger-2F0 a b x))))

(defn tricomis-U
  "Confluent hypergeometric function U of the second kind, also known as Tricomi's function.

  The second, generally unbounded-as-x-approaches-0 solution of Kummer's differential equation (the first being [[kummers-M]]), related to it for non-integer `b` via `U(a,b,x) = Gamma(1-b)/Gamma(a-b+1) kummers-M(a,b,x) + Gamma(b-1)/Gamma(a) x^(1-b) kummers-M(a-b+1,2-b,x)`.

  Parameters:

  - `a`, `b` (double): the two parameters.
  - `x` (double): the argument. The validated domain is `x >= 0.0`; `x < 0.0` is not a domain this function is intended for (typically `##NaN`, though a handful of specific parameter coincidences may return an unvalidated real number instead).

  At `x = 0.0`: returns `Gamma(1-b)/Gamma(a-b+1)` for `b < 1.0`. For `b >= 1.0`, `U` genuinely diverges there, returning a signed `##Inf` following the sign of `Gamma(a)`, unless `a` is itself a non-positive integer, in which case the divergence is removable and the result is the finite value `(-1)^n (b)_n` (writing `a = -n`, and `(b)_n` the rising Pochhammer symbol). Elsewhere, real and finite.

  See also [[kummers-M]], [[whittaker-W]], [[hypergeometric-2F0]]."
  ^double [^double a ^double b ^double x]
  (cond
    ;; wolfram alpha
    (m/zero? x) (cond
                  ;; a = -n a non-positive integer: U(-n,b,x) is finite even
                  ;; at x=0 for any b (a degenerate case of the general
                  ;; asymptotics below, where the divergent term's
                  ;; coefficient 1/gamma(a) vanishes at gamma's pole);
                  ;; confirmed against independent limit evaluations at
                  ;; small x>0 to equal (-1)^n (b)_n (rising Pochhammer
                  ;; symbol), for both b<1 and b>=1
                  (and (m/not-pos? a) (m/integer? a))
                  (let [n (unchecked-long (m/- a))]
                    (m/* (if (m/odd? n) -1.0 1.0) (m/rising-factorial-int n b)))
                  (m/< b 1.0) (m// (gamma (m/- 1.0 b)) (gamma (m/inc (m/- a b))))
                  ;; b >= 1.0, a not a non-positive integer: genuinely
                  ;; diverges as x -> 0+ (a pole at x=0), for b>1 with
                  ;; sign(gamma(b-1)/gamma(a)) = sign of gamma(a) alone
                  ;; (gamma(b-1) is always positive there, b-1>0), and,
                  ;; confirmed separately, the same sign(gamma(a)) rule also
                  ;; holds for the logarithmic divergence at b=1 exactly
                  :else (m/copy-sign ##Inf (gamma a)))
    (m/== a b) (let [ex (m/exp x)]
                 (if (m/inf? ex)
                   ;; exp(x) alone overflows to Infinity even though the
                   ;; true value (exp(x) * a rapidly-decaying incomplete
                   ;; gamma) is a tiny, perfectly finite double there;
                   ;; falls back to the general asymptotic-series formula
                   ;; below, which never computes exp(x) directly and was
                   ;; confirmed to already agree with this branch to ~13
                   ;; significant digits well before the overflow boundary
                   (m/* (m/pow x (m/- a)) (hypergeometric-2F0 a 1.0 (m/- (m// x))))
                   (m/* ex (upper-incomplete-gamma (m/- 1.0 a) x))))
    (m/== a (m/dec b)) (m/pow x (m/- a))    
    :else (m/* (m/pow x (m/- a)) (hypergeometric-2F0 a (m/inc (m/- a b)) (m/- (m// x))))))

(defn whittaker-W
  "Whittaker's W function.

  A standard solution of Whittaker's differential equation, expressed via the confluent hypergeometric function of the second kind [[tricomis-U]]: `W(kappa,mu,x) = exp(-x/2) x^(mu+1/2) tricomis-U(mu-kappa+1/2, 2mu+1, x)`.

  Parameters:

  - `kappa` (double): the first parameter.
  - `mu` (double): the second parameter.
  - `x` (double): the argument, restricted to `x >= 0.0` (for `x < 0.0`, `x^(mu+1/2)` generally leaves the real line, so `##NaN` is returned).

  As `x` approaches `0.0` from above, the result approaches `0.0` when `mu > -0.5`, diverges (with a sign that follows [[tricomis-U]]'s own limiting value there, unlike [[whittaker-M]] this is not always positive) when `mu < -0.5`, and approaches [[tricomis-U]]'s own value at that point when `mu = -0.5` exactly (`x^(mu+1/2)` reduces to the constant `1` there). Elsewhere, inherits [[tricomis-U]]'s domain and pole structure through the transformed parameters `mu-kappa+1/2` and `2mu+1`.

  See also [[tricomis-U]], [[kummers-M]], [[whittaker-M]]."
  ^double [^double kappa ^double mu ^double x]
  (let [mu+05 (m/+ 0.5 mu)]
    (cond
      (m/neg? x) ##NaN
      (and (m/zero? x) (m/zero? mu+05)) (tricomis-U (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x)
      :else (let [z (m/exp (m/* 0.5 (m/+ (m/* -0.5 x) (m/* mu+05 (m/log x)))))]
              (m/* z (tricomis-U (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x) z)))))

(defn hypergeometric-2F1
  "Gauss's hypergeometric ₂F₁ function.

  The classical generalized hypergeometric series with two numerator and one denominator parameter, `sum_n (a)_n (b)_n / (c)_n * x^n / n!` (using the Pochhammer symbol `(a)_n`), convergent for `|x| < 1` and analytically continued elsewhere on the real line.

  Parameters:

  - `a`, `b` (double): the two numerator parameters.
  - `c` (double): the denominator parameter.
  - `x` (double): the argument.

  Returns `1.0` for `x = 0.0`. When `a` or `b` is a non-positive integer, the series terminates to a finite polynomial in `x`, real for every `x` (subject to the usual pole in `c`, described below). Otherwise: for `x < 1.0`, real and finite; at `x = 1.0`, follows Gauss's summation theorem, finite when `c-a-b > 0` and diverging to `##Inf` when `c-a-b <= 0`; for `x > 1.0`, the analytic continuation leaves the real line for generic `a`, `b`, so `##NaN` is returned there (except in the terminating case above, which remains real for any `x`).

  Has a genuine pole (returns a plain `##Inf`, matching mpmath's own `hyp2f1` convention there) whenever `c` is a non-positive integer, unless the series terminates (`a` or `b` a non-positive integer, reaching that same magnitude no later than `c` does) before reaching it.

  See also [[hypergeometric-1F1]], [[hypergeometric-2F0]], [[hypergeometric-pFq]]."
  ^double [^double a ^double b ^double c ^double x]
  (hg/hypergeometric-2F1 a b c x))

;;

(defn- pfq-cancel-equal-pairs
  "Cancels any numerator parameter that exactly equals a denominator
  parameter, PROVIDED the shared value is not a non-positive integer:
  `(a_i)_k / (b_j)_k = 1` identically for every k when `a_i = b_j` and
  neither Pochhammer symbol ever hits zero, so such a pair cancels in the
  pFq definition, reducing `pFq` to `(p-1)F(q-1)`. When the shared value IS
  a non-positive integer, this simple cancellation does NOT hold (both
  Pochhammer symbols vanish together at the same term, a distinct
  removable coincidence handled separately by
  `pfq-integer-coincidence-pair`/`pfq-truncated-sum`); such pairs are left
  untouched here. Returns `[ps' qs']` with every safe pair removed
  (repeated until none remain)."
  [ps qs]
  (loop [ps (vec ps) qs (vec qs)]
    (let [match (first (for [i (range (count ps))
                              j (range (count qs))
                              :let [a (double (ps i))]
                              :when (and (== a (double (qs j)))
                                         (not (nonpos-int? a)))]
                          [i j]))]
      (if match
        (let [[^long i ^long j] match]
          (recur (into (subvec ps 0 i) (subvec ps (inc i)))
                 (into (subvec qs 0 j) (subvec qs (inc j)))))
        [ps qs]))))

(defn- pfq-terminating-n
  "Minimum `n` such that some numerator parameter equals `-n` (a
  non-positive integer): the series terminates to a finite polynomial at
  term index `n`. Returns `nil` if no numerator parameter terminates it."
  [ps]
  (let [ns (keep (fn [^double a] (when (nonpos-int? a) (long (m/- a)))) ps)]
    (when (seq ns) (long (apply min ns)))))

(defn- pfq-pole-m
  "Minimum `m` such that some denominator parameter equals `-m` (a
  non-positive integer): a genuine pole at term index `m`, unless
  preempted by numerator termination at or before it. Returns `nil` if no
  denominator parameter is a non-positive integer."
  [qs]
  (let [ms (keep (fn [^double b] (when (nonpos-int? b) (long (m/- b)))) qs)]
    (when (seq ms) (long (apply min ms)))))

(defn- pfq-integer-coincidence-pair
  "Finds a numerator/denominator index pair `[i j]` where `ps[i] = qs[j]`
  and that shared value is a non-positive integer (the case
  `pfq-cancel-equal-pairs` deliberately leaves untouched). Returns `nil`
  if there is none."
  [ps qs]
  (first (for [i (range (count ps))
               j (range (count qs))
               :let [a (double (nth ps i))]
               :when (and (== a (double (nth qs j))) (nonpos-int? a))]
           [i j])))

(defn- pfq-truncated-sum
  "Evaluates `sum_{k=0}^{n} prod(ps)_k / prod(qs)_k * x^k / k!` directly,
  term by term, bounded to `k <= n` (never computing a term beyond that
  bound). Used to resolve the removable `a_i = b_j` non-positive-integer
  coincidence: with that pair excluded from `ps`/`qs` here, this is safe
  even though `ps`/`qs` may still contain other non-positive integers
  themselves, since none of their own Pochhammer symbols can reach zero
  within `k <= n`, as confirmed against mpmath."
  ^double [ps qs ^double x ^long n]
  (let [^doubles ps (m/seq->double-array ps)
        ^doubles qs (m/seq->double-array qs)
        lp (alength ps)
        lq (alength qs)]
    (loop [k (long 0) term 1.0 s 1.0]
      (if (m/== k n)
        s
        (let [num (double (loop [i (long 0) v x]
                            (if (m/== i lp) v (recur (m/inc i) (m/* v (m/+ (Array/aget ps i) k))))))
              den (double (loop [i (long 0) v (double (m/inc k))]
                            (if (m/== i lq) v (recur (m/inc i) (m/* v (m/+ (Array/aget qs i) k))))))
              term (m/* term (m// num den))]
          (recur (m/inc k) term (m/+ s term)))))))

(defn hypergeometric-pFq
  "Generalized hypergeometric function pFq with p numerator and q denominator parameters.

  Parameters:

  - `ps` (sequence of double): the `p` numerator parameters.
  - `qs` (sequence of double): the `q` denominator parameters.
  - `x` (double): the argument.
  - `max-iters` (long, optional): maximum series/acceleration iterations, `1048576` by default.

  Whenever a numerator parameter exactly equals a denominator parameter, that pair cancels identically, reducing to a lower `(p-1)F(q-1)`; this is applied before anything else below. If that shared value is also a non-positive integer, the pair does not simply cancel (both Pochhammer symbols vanish together at the same term instead); the series is then evaluated as a truncated sum of the further-reduced coefficients.

  Convergence and evaluation strategy depend on `p` compared to `q`: for `p <= q` the series converges for every `x`; for `p = q + 1` it converges for `x < 1` (and follows the usual boundary/continuation behavior at and beyond `x = 1`); for `p > q + 1` the series formally diverges for any nonzero `x` unless it terminates (see below), in which case only an asymptotic resummation is used.

  Whenever a numerator parameter is a non-positive integer, the series terminates to an exact finite polynomial in `x`, valid for any `x`, regardless of the general convergence class above.

  Has a genuine pole (returns `##NaN`) whenever a denominator parameter is a non-positive integer that is reached before the series would otherwise terminate.

  Returns `1.0` for `x = 0.0`.

  For the generic (non-terminating), formally divergent `p > q + 1` case, the result relies on Weniger-acceleration resummation. A missing finite-value fallback at that acceleration's own loop exit (present in this project's other single-arity Weniger ports, and in the Julia `HypergeometricFunctions.jl` reference implementation this was ported from, but dropped from this general port) used to turn many cases of an already-converged intermediate value overflowing into `##NaN`; fixed. A minority of `p > q + 1` points remain genuinely unreliable regardless -- confirmed, by comparison against Julia's own raw Weniger kernel and its production `pFq` dispatcher (neither of which offers an alternative method for this regime either), to be an inherent limitation of the resummation itself, not specific to this port; not corrected here.

  See also [[hypergeometric-0F0]], [[hypergeometric-1F0]], [[hypergeometric-0F1]], [[hypergeometric-1F1]], [[hypergeometric-0F2]], [[hypergeometric-2F0]], [[hypergeometric-2F1]], [[kummers-M]], [[tricomis-U]]."
  (^double [ps qs ^double x] (hypergeometric-pFq ps qs x 1048576))
  (^double [ps qs ^double x ^long max-iters]
   (let [[ps qs] (pfq-cancel-equal-pairs ps qs)
         p (count ps) q (count qs)
         n (pfq-terminating-n ps)
         m (pfq-pole-m qs)]
     (cond
       ;; x=0 always wins, even over an otherwise-genuine pole or a NaN
       ;; parameter, matching the convention used throughout this namespace
       ;; (e.g. hypergeometric-2F1)
       (m/< (m/abs x) m/MACHINE-EPSILON10) 1.0

       ;; the degenerate p=0, q=0 case is identically exp(x); the general
       ;; Weniger acceleration below isn't scaled for it and loses all
       ;; precision once exp(x) gets very small for negative x
       (and (m/zero? p) (m/zero? q)) (m/exp x)

       ;; a numerator parameter still equals a denominator parameter and
       ;; that shared value is a non-positive integer -n: both Pochhammer
       ;; symbols vanish together at term n+1, a removable coincidence
       ;; whose value is the truncated (that pair removed) series up to
       ;; k=n, confirmed against mpmath; NOT a plain cancellation (which
       ;; would wrongly give the full, non-truncated reduced series), and
       ;; NOT safe to evaluate via hypergeometric-pFq-maclaurin either
       ;; (its ratio recursion hits the same 0/0 exactly at k=n)
       (and m n (m/== (long m) (long n)))
       (let [[^long i ^long j] (pfq-integer-coincidence-pair ps qs)]
         (pfq-truncated-sum (into (subvec ps 0 i) (subvec ps (inc i)))
                            (into (subvec qs 0 j) (subvec qs (inc j)))
                            x n))

       ;; a denominator parameter is a genuine, unavoidable pole only if
       ;; reached before any numerator termination
       (and m (or (nil? n) (m/< (long m) (long n)))) ##NaN

       ;; series terminates to an exact finite polynomial for any x; the
       ;; MacLaurin loop evaluates this exactly, unlike the Weniger
       ;; resummation path, which is built for infinite series and loses
       ;; precision as |x| grows in this case
       n (hg/hypergeometric-pFq-maclaurin ps qs x max-iters)

       ;; entire function (p<=q, converges for any x): MacLaurin is accurate
       ;; for x>=0, but suffers catastrophic cancellation for x<0, where
       ;; Weniger acceleration is needed instead
       (m/<= p q) (if (m/pos? x)
                    (hg/hypergeometric-pFq-maclaurin ps qs x max-iters)
                    (hg/hypergeometric-pFq-weniger ps qs x max-iters))
       (m/== p (m/inc q)) (if (m/< (m/abs x) 0.72)
                            (hg/hypergeometric-pFq-maclaurin ps qs x max-iters)
                            (hg/hypergeometric-pFq-weniger ps qs x max-iters))
       :else (hg/hypergeometric-pFq-weniger ps qs x max-iters)))))

(set! *unchecked-math* true)

(defn- local-abs [a] (if (pos? a) a (- a)))

(def ^:private RMEPS (rationalize m/MACHINE-EPSILON10))
(def ^:pribate RATIO-1 (rationalize 1))

(defn hypergeometric-pFq-ratio
  "Hypergeometric-pFq using MacLaurin series on ratios. Can be very slow.

  `max-iters` is set to 10000 by default."
  (^double [ps qs z] (hypergeometric-pFq-ratio ps qs z 10000))
  (^double [ps qs z ^long max-iters]
   (let [a (map rationalize ps)
         b (map rationalize qs)
         z (rationalize z)
         eps RMEPS]
     (loop [k (long 1)
            s0 RATIO-1
            s1 (+ 1 (/ (* z (reduce * RATIO-1 a)) (reduce * RATIO-1 b)))]
       (if (and (m/< k max-iters)
                (> (local-abs (- s1 s0)) (max eps (* eps (max s1 s0)))))
         (let [rk (/ z (+ k 1))
               rk (reduce (fn [r va] (* r (+ va k))) rk a)
               rk (reduce (fn [r vb] (/ r (+ vb k))) rk b)]
           (recur (m/inc k) s1 (+ s1 (* (- s1 s0) rk))))
         s1)))))

(set! *unchecked-math* :warn-on-boxed)

;; Owen's T
;; https://github.com/JuliaStats/StatsFuns.jl/issues/99

(defn owens-t
  "Owens' T function"
  ^double [^double h ^double a]
  (cond
    (m/zero? a) 0.0
    (m/zero? h) (m/* m/INV_TWO_PI (m/atan a))
    (m/neg? h) (recur (m/- h) a)
    (m/neg? a) (m/- (owens-t h (m/- a)))
    (m/one? a) (let [h' (m/* h m/INV_SQRT_2)]
                 (m/* 0.125 (erfc (m/- h')) (erfc h')))
    (m/inf? a) (m/* 0.25 (erfc (m/* (m/abs h) m/INV_SQRT_2) ))
    (m/> a 1.0) (let [h' (m/* -1.0 h m/INV_SQRT_2)
                      e1 (erfc h')
                      e2 (erfc (m/* a h'))]
                  (m/- (m/* 0.25 (m/+ e1 e2))
                       (m/* 0.25 e1 e2)
                       (owens-t (m/* a h) (m// a))))
    (m/< a 0.999999) (let [hh2 (m/* -0.5 h h)
                           t2 (map (fn [^double x]
                                     (let [ax2+1 (m/inc (m/* a a x x))]
                                       (m/* m/INV_FOUR_PI a
                                            (m// (m/exp (m/* hh2 ax2+1))
                                                 ax2+1))))
                                   [-0.9987710072524261, -0.9935301722663508, -0.9841245837228269, -0.9705915925462473, -0.9529877031604309, -0.9313866907065543, -0.9058791367155696
                                    , -0.8765720202742479, -0.8435882616243935, -0.8070662040294426, -0.7671590325157404, -0.7240341309238146, -0.6778723796326639, -0.6288673967765136
                                    , -0.5772247260839727, -0.5231609747222331, -0.4669029047509584, -0.4086864819907167, -0.34875588629216075, -0.28736248735545555, -0.22476379039468905
                                    , -0.1612223560688917, -0.0970046992094627, -0.03238017096286937, 0.03238017096286937, 0.0970046992094627, 0.1612223560688917, 0.22476379039468905
                                    , 0.28736248735545555, 0.34875588629216075, 0.4086864819907167, 0.4669029047509584, 0.5231609747222331, 0.5772247260839727, 0.6288673967765136
                                    , 0.6778723796326639, 0.7240341309238146, 0.7671590325157404, 0.8070662040294426, 0.8435882616243935, 0.8765720202742479, 0.9058791367155696
                                    , 0.9313866907065543, 0.9529877031604309, 0.9705915925462473, 0.9841245837228269, 0.9935301722663508, 0.9987710072524261])]
                       (v/dot t2 [0.0031533460523059122, 0.0073275539012762885, 0.011477234579234613, 0.015579315722943824, 0.01961616045735561, 0.023570760839324363
                                  , 0.027426509708356944, 0.031167227832798003, 0.03477722256477052, 0.038241351065830737, 0.04154508294346467, 0.0446745608566943, 0.04761665849249045
                                  , 0.05035903555385445, 0.05289018948519363, 0.055199503699984116, 0.05727729210040322, 0.05911483969839564, 0.06070443916589387, 0.06203942315989268
                                  , 0.06311419228625402, 0.06392423858464813, 0.06446616443594998, 0.06473769681268386, 0.06473769681268386, 0.06446616443594998, 0.06392423858464813
                                  , 0.06311419228625402, 0.06203942315989268, 0.06070443916589387, 0.05911483969839564, 0.05727729210040322, 0.055199503699984116
                                  , 0.05289018948519363, 0.05035903555385445, 0.04761665849249045, 0.0446745608566943, 0.04154508294346467, 0.038241351065830737, 0.03477722256477052
                                  , 0.031167227832798003, 0.027426509708356944, 0.023570760839324363, 0.01961616045735561, 0.015579315722943824, 0.011477234579234613, 0.0073275539012762885
                                  , 0.0031533460523059122]))
    :else (let [j (m/* 0.5 (erfc (m/* -1.0 h m/INV_SQRT_2)))
                a- (m/- 1.0 a)
                k (m/atan (m// a- (m/inc a)))]
            (m/- (m/* 0.5 j (m/- 1.0 j))
                 (m/* m/INV_TWO_PI k (m/exp (m// (m/* -0.5 a- h h) k)))))))


;;;;;

;; complex
;;;;;

(defn- cplx-nonpos-int?
  "Is `z` a non-positive real integer (zero imaginary part)?"
  [^Vec2 z]
  (and (m/zero? (cplx/im z)) (nonpos-int? (cplx/re z))))

(defn- pfq-cancel-equal-pairs-complex
  "Complex analogue of [[pfq-cancel-equal-pairs]]: cancels any numerator
  parameter that exactly equals a denominator parameter, provided the
  shared value is not a non-positive real integer."
  [ps qs]
  (loop [ps ps qs qs]
    (let [match (first (for [i (range (count ps))
                             j (range (count qs))
                             :when (and (= (ps i) (qs j))
                                        (not (cplx-nonpos-int? (ps i))))]
                         [i j]))]
      (if match
        (let [[^long i ^long j] match]
          (recur (into (subvec ps 0 i) (subvec ps (inc i)))
                 (into (subvec qs 0 j) (subvec qs (inc j)))))
        [ps qs]))))

(defn- pfq-terminating-n-complex
  "Complex analogue of [[pfq-terminating-n]]."
  [ps]
  (let [ns (keep (fn [z] (when (cplx-nonpos-int? z) (long (m/- (cplx/re z))))) ps)]
    (when (seq ns) (long (apply m/min ns)))))

(defn- pfq-pole-m-complex
  "Complex analogue of [[pfq-pole-m]]."
  [qs]
  (let [ms (keep (fn [z] (when (cplx-nonpos-int? z) (long (m/- (cplx/re z))))) qs)]
    (when (seq ms) (long (apply m/min ms)))))

(defn- pfq-integer-coincidence-pair-complex
  "Complex analogue of [[pfq-integer-coincidence-pair]]."
  [ps qs]
  (first (for [i (range (count ps))
               j (range (count qs))
               :when (and (= (ps i) (qs j)) (cplx-nonpos-int? (ps i)))]
           [i j])))

(defn- pfq-truncated-sum-complex
  "Complex analogue of [[pfq-truncated-sum]]."
  ^Vec2 [ps qs z ^long n]
  (let [ps (mapv cplx/ensure-complex ps)
        qs (mapv cplx/ensure-complex qs)
        lp (count ps)
        lq (count qs)]
    (loop [k (long 0) term cplx/ONE s cplx/ONE]
      (if (m/== k n)
        s
        (let [num (loop [i (long 0) v z]
                    (if (m/== i lp) v (recur (m/inc i) (cplx/mult v (cplx/adds (ps i) k)))))
              den (loop [i (long 0) v (cplx/complex (double (m/inc k)) 0.0)]
                    (if (m/== i lq) v (recur (m/inc i) (cplx/mult v (cplx/adds (qs i) k)))))
              term (cplx/mult term (cplx/div num den))]
          (recur (m/inc k) term (cplx/add s term)))))))

(defrecord PfQComplexData [done? value ps qs method])

(defn- pfq-complex-data-true [value] (PfQComplexData. true value nil nil nil))
(defn- pfq-complex-data-false [ps qs method] (PfQComplexData. false nil ps qs method))

(defn- pfq-complex-route
  "Classifies a pFq-complex call before any numeric acceleration method is
  attempted -- the exact pre-dispatch checks [[hypergeometric-pFq-complex]]
  itself performs (equal-pair cancellation, exact termination, poles, the
  degenerate `p=0,q=0` case, `z=0`), factored out so a second caller can
  reuse them without going through the public dispatcher.

  Returns either a directly-known final value (`{:done? true :value v}`),
  or a residual call still needing numeric evaluation (`{:done? false
  :method (:maclaurin|:weniger) :ps :qs}`, with `:ps`/`:qs` already
  reduced by any equal-pair cancellation).

  Precondition: `ps`/`qs` are sequences of real or [[Vec2]] numbers, `z`
  is already a [[Vec2]] (via [[fastmath.complex/ensure-complex]]).
  Postcondition: exactly one of `:done?`'s two shapes above, never both.

  Shared by [[hypergeometric-pFq-complex]] (the public dispatcher) and
  `tricomis-U-complex-asymptotic` (which calls Weniger acceleration
  directly, bypassing the dispatcher, but still needs these same checks
  -- see [[tricomis-U-complex-asymptotic]]'s own docstring)."
  [ps qs ^Vec2 z ^long max-iters]
  (let [ps (mapv cplx/ensure-complex ps)
        qs (mapv cplx/ensure-complex qs)
        [ps qs] (pfq-cancel-equal-pairs-complex ps qs)
        p (count ps) q (count qs)
        n (pfq-terminating-n-complex ps)
        m (pfq-pole-m-complex qs)]
    (cond
      ;; z=0 always wins, even over an otherwise-genuine pole
      (m/< (cplx/abs z) m/MACHINE-EPSILON10) (pfq-complex-data-true cplx/ONE)

      ;; the degenerate p=0, q=0 case is identically exp(z); the general
      ;; Weniger acceleration below isn't scaled for it and loses all
      ;; precision once exp(z) gets very small
      (and (m/zero? p) (m/zero? q)) (pfq-complex-data-true (cplx/exp z))

      ;; a numerator parameter still equals a denominator parameter and
      ;; that shared value is a non-positive real integer -n: a removable
      ;; coincidence resolved as a truncated sum, see
      ;; [[pfq-truncated-sum-complex]] and [[hypergeometric-pFq]]'s own
      ;; docstring for the real-valued derivation this mirrors
      (and m n (m/== (long m) (long n)))
      (let [[^long i ^long j] (pfq-integer-coincidence-pair-complex ps qs)]
        (pfq-complex-data-true (pfq-truncated-sum-complex (into (subvec ps 0 i) (subvec ps (inc i)))
                                                          (into (subvec qs 0 j) (subvec qs (inc j)))
                                                          z n)))

      ;; a denominator parameter is a genuine, unavoidable pole only if
      ;; reached before any numerator termination
      (and m (or (nil? n) (m/< (long m) (long n)))) (pfq-complex-data-true (Vec2. ##NaN ##NaN))

      ;; series terminates to an exact finite polynomial for any z
      n (pfq-complex-data-true (hg/hypergeometric-pFq-maclaurin-complex ps qs z max-iters))

      ;; entire function (p<=q, converges for any z): MacLaurin is
      ;; accurate for a non-negative real part, Weniger acceleration is
      ;; needed otherwise (catastrophic cancellation)
      (m/<= p q) (pfq-complex-data-false ps qs (if (m/pos? (cplx/re z)) :maclaurin :weniger))

      (m/== p (m/inc q)) (pfq-complex-data-false ps qs (if (m/< (cplx/abs z) 0.72) :maclaurin :weniger))

      :else (pfq-complex-data-false ps qs :weniger))))

(defn hypergeometric-pFq-complex
  "Generalized hypergeometric function pFq with p numerator and q denominator complex parameters.

  Complex counterpart of [[hypergeometric-pFq]]; see that docstring for the general shape and dispatch logic (entire for `p <= q`, radius-1 for `p = q+1`, formally divergent unless terminating for `p > q+1`), which carries over here with `x` replaced by the complex argument `z` and the `p <= q`/`p = q+1` MacLaurin-vs-Weniger split now made on the sign of `z`'s real part / `|z|` respectively.

  Parameters:

  - `ps` (sequence of double or [[Vec2]]): the `p` numerator parameters, real or complex (promoted via `ensure-complex`).
  - `qs` (sequence of double or [[Vec2]]): the `q` denominator parameters.
  - `z` (double or [[Vec2]]): the argument, real or complex.
  - `max-iters` (long, optional): maximum series/acceleration iterations, `1048576` by default.

  Whenever a numerator parameter exactly equals a denominator parameter, that pair cancels identically, reducing to a lower `(p-1)F(q-1)`; this is applied before anything else below. If that shared value is also a non-positive real integer (zero imaginary part), the pair does not simply cancel (both Pochhammer symbols vanish together at the same term instead); the series is then evaluated as a truncated sum of the further-reduced coefficients.

  Whenever a numerator parameter is a non-positive real integer, the series terminates to an exact finite polynomial in `z`, valid for any `z`.

  Has a genuine pole (returns `(Vec2. ##NaN ##NaN)`) whenever a denominator parameter is a non-positive real integer that is reached before the series would otherwise terminate.

  Returns `1.0+0.0i` for `z = 0.0+0.0i`.

  For the generic (non-terminating), formally divergent `p > q + 1` case, the result relies entirely on Weniger-acceleration resummation (no MacLaurin fallback exists there, since the underlying series is genuinely divergent). A missing finite-value fallback at that acceleration's own loop exit used to turn many cases of an already-converged intermediate value overflowing into `##NaN` regardless of the sign of `z`; fixed (see [[fastmath.special.hypergeometric/hypergeometric-pFq-weniger-complex-with-reason]] for the underlying mechanism). A minority of `p > q + 1` points remain genuinely unreliable regardless -- confirmed to be an inherent limitation of the resummation itself, not specific to this port (see [[hypergeometric-pFq]]'s own docstring). For `p = q + 1` beyond the MacLaurin radius (`|z| >= 0.72`), Weniger acceleration is markedly more reliable on both sides of the real axis, though (as for the real-valued [[hypergeometric-pFq]]) precision still degrades gradually as `|z|` grows very large; not corrected here.

  See also [[hypergeometric-pFq]]."
  (^Vec2 [ps qs z] (hypergeometric-pFq-complex ps qs z 1048576))
  (^Vec2 [ps qs z ^long max-iters]
   (let [z (cplx/ensure-complex z)
         ^PfQComplexData result (pfq-complex-route ps qs z max-iters)]
     (if (.done? result)
       (.value result)
       (case (.method result)
         :maclaurin (hg/hypergeometric-pFq-maclaurin-complex (.ps result) (.qs result) z max-iters)
         :weniger (hg/hypergeometric-pFq-weniger-complex (.ps result) (.qs result) z max-iters))))))

(defn- complex-log-gamma-asymptotic
  ^Vec2 [^Vec2 z]
  (let [zinv (cplx/reciprocal z)
        t (cplx/sq zinv)]
    (cplx/add (cplx/adds (cplx/sub (cplx/mult (cplx/adds z -0.5)
                                              (cplx/log z))
                                   z)
                         9.1893853320467274178032927e-01)
              (cplx/mult zinv (poly/mevalpoly-scalar-complex t
                                8.3333333333333333333333368e-02,-2.7777777777777777777777776e-03,
                                7.9365079365079365079365075e-04,-5.9523809523809523809523806e-04,
                                8.4175084175084175084175104e-04,-1.9175269175269175269175262e-03,
                                6.4102564102564102564102561e-03,-2.9550653594771241830065352e-02)))))

(defn log-gamma-complex
  "Logarithm of the complex gamma function (principal branch).

  Parameters:

  - `z` ([[Vec2]]): the argument.

  Returns `log(Gamma(z))`. At the poles of Gamma (`z` a non-positive real integer, zero included), Gamma genuinely diverges there without a single well-defined direction-independent limit in the complex plane; this returns `(Vec2. ##Inf ##Inf)` for every such `z` except `z = 0.0+0.0i` exactly, which instead returns a finite imaginary part (`##Inf` real part only) following the standard convention for the principal branch approaching along the positive real axis. `##NaN`/`##Inf` real or imaginary input parts propagate following the usual IEEE conventions.

  See also [[gamma-complex]], [[gamma]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)
        yabs (m/abs y)]
    (cond
      (and (m/inf? x) (m/valid-double? y)) (Vec2. x (if (m/pos? x)
                                                      (if (m/zero? y) y (m/copy-sign ##Inf y))
                                                      (m/copy-sign ##Inf (m/- y))))
      (and (m/valid-double? x) (m/inf? y)) (Vec2. ##-Inf y)
      (or (m/invalid-double? x) (m/invalid-double? y)) (Vec2. ##NaN ##NaN)
      (and (m/zero? x) (m/zero? y)) (Vec2. ##Inf (if (m/negative-zero? x) (m/copy-sign m/PI (m/- y)) (m/- y)))
      (and (m/integer? x) (m/neg? x) (m/zero? y)) (Vec2. ##Inf ##Inf)

      (or (m/> x 7.0) (m/> yabs 7.0)) (complex-log-gamma-asymptotic z)
      
      (m/< x 0.1) (cplx/sub
                   (cplx/sub (Vec2. m/LOG_PI (m/* (m/copy-sign m/TWO_PI y) (m/floor (m/+ 0.25 (m/* 0.5 x)))))
                             (cplx/log (cplx/sin (cplx/scale z m/PI))))
                   (log-gamma-complex (cplx/sub cplx/ONE z)))

      (m/< (m/+ (m/abs (m/dec x)) yabs) 0.1)
      (let [w (cplx/adds z -1.0)]
        (cplx/mult w (poly/mevalpoly-scalar-complex w
                       -5.7721566490153286060651188e-01,8.2246703342411321823620794e-01,
                       -4.0068563438653142846657956e-01,2.705808084277845478790009e-01,
                       -2.0738555102867398526627303e-01,1.6955717699740818995241986e-01,
                       -1.4404989676884611811997107e-01,1.2550966952474304242233559e-01,
                       -1.1133426586956469049087244e-01,1.000994575127818085337147e-01,
                       -9.0954017145829042232609344e-02,8.3353840546109004024886499e-02,
                       -7.6932516411352191472827157e-02,7.1432946295361336059232779e-02,
                       -6.6668705882420468032903454e-02)))

      (m/< (m/+ (m/abs (m/- x 2.0)) yabs) 0.1)
      (let [w (cplx/adds z -2.0)]
        (cplx/mult w (poly/mevalpoly-scalar-complex w
                       4.2278433509846713939348812e-01,3.2246703342411321823620794e-01,
                       -6.7352301053198095133246196e-02,2.0580808427784547879000897e-02,
                       -7.3855510286739852662729527e-03,2.8905103307415232857531201e-03,
                       -1.1927539117032609771139825e-03,5.0966952474304242233558822e-04,
                       -2.2315475845357937976132853e-04,9.9457512781808533714662972e-05,
                       -4.4926236738133141700224489e-05,2.0507212775670691553131246e-05)))

      :else (loop [^Vec2 shiftprod (Vec2. x yabs)
                   x (m/inc x)
                   sb false
                   signflips (long 0)]
              (if (m/<= x 7.0)
                (let [^Vec2 nsp (cplx/mult shiftprod (Vec2. x yabs))
                      nsb (or (m/neg? (.y nsp)) (m/negative-zero? (.y nsp)))]
                  (recur nsp (m/inc x) nsb (if (and nsb (not= sb nsb)) (m/inc signflips) signflips)))
                (let [^Vec2 s (cplx/log shiftprod)
                      shift (Vec2. (.x s)
                                   (if (or (m/neg? y) (m/negative-zero? y))
                                     (m/- (m/* m/-TWO_PI signflips) (.y s))
                                     (m/+ (m/* m/TWO_PI signflips) (.y s))))]
                  (cplx/sub (complex-log-gamma-asymptotic (Vec2. x y)) shift)))))))

(defn gamma-complex
  "Complex version of the gamma function (principal branch), `exp(log-gamma-complex z)`.

  Parameters:

  - `z` ([[Vec2]]): the argument.

  Returns `Gamma(z)`. At most poles of Gamma (`z` a negative real integer), returns `##NaN`, matching the real-valued [[gamma]]'s own convention at its poles (unlike [[log-gamma-complex]], which encodes those same poles as an infinite log value with an equally infinite imaginary part; exponentiating that back does not recover any single well-defined complex infinity, since Gamma has no direction-independent limit there). At `z = 0.0+0.0i` specifically, [[log-gamma-complex]] instead encodes the pole with a finite imaginary part (approaching along the positive real axis), so this returns a genuine signed infinity there (`##Inf` real part) instead of `##NaN`, unlike [[gamma]]'s own `x = 0.0` convention.

  See also [[log-gamma-complex]], [[gamma]]."
  ^Vec2 [^Vec2 z] (cplx/exp (log-gamma-complex z)))

(defn- reciprocal-gamma-complex
  "1/Gamma(z), entire (no poles): identically 0 whenever z is a
  non-positive real integer (Gamma's own poles), computed elsewhere as
  `exp(-log-gamma-complex z)` rather than `1/gamma-complex z`, since the
  latter propagates ##NaN even away from a pole's exact real-integer
  coordinates (an intermediate 0*Inf artifact whenever exp is applied to
  log-gamma-complex's own pole encoding)."
  ^Vec2 [^Vec2 z]
  (if (cplx-nonpos-int? z)
    cplx/ZERO
    (cplx/exp (cplx/neg (log-gamma-complex z)))))

(defn- cplx-int?
  "Is `z` a real integer of any sign, including zero (zero imaginary part)?"
  ^Boolean [^Vec2 z]
  (and (m/zero? (cplx/im z)) (m/integer? (cplx/re z))))

(defn- tricomis-U-complex-raw
  "Tricomi's `U(a,b,z)` reflection formula in terms of two Kummer `M`
  functions, valid only for `b` not an integer: its `pi/sin(pi b)`
  prefactor, and both bracketed terms individually via
  [[reciprocal-gamma-complex]], vanish identically at integer `b` (a
  removable singularity of this particular formula, not of `U` itself),
  resolved by [[tricomis-U-complex]] instead of here."
  ^Vec2 [^Vec2 a ^Vec2 b ^Vec2 z]
  (let [p1 (cplx/sub (cplx/add cplx/ONE a) b)
        p2 (cplx/sub cplx/TWO b)]
    (-> (cplx/sub (-> (hypergeometric-pFq-complex [a] [b] z)
                      (cplx/mult (reciprocal-gamma-complex p1))
                      (cplx/mult (reciprocal-gamma-complex b)))
                  (-> (hypergeometric-pFq-complex [p1] [p2] z)
                      (cplx/mult (reciprocal-gamma-complex a))
                      (cplx/mult (reciprocal-gamma-complex p2))
                      (cplx/mult (cplx/pow z (cplx/sub cplx/ONE b)))))
        (cplx/mult (cplx/div cplx/PI (cplx/sin (cplx/scale b m/PI)))))))

(defn- tricomis-U-complex-raw-limit
  "[[tricomis-U-complex-raw]], but safe for integer `b` (its removable
  singularity there, see that function's docstring): evaluated via a
  Richardson-extrapolated numerical limit instead -- evaluating the
  (elsewhere correct) formula at `b +/- ` a small offset and extrapolating
  to the `b = integer` limit, since `U` is itself continuous there (just
  this one formula for it isn't) -- verified against mpmath to about
  1e-7..1e-12 relative accuracy across several integer `b`, safer than
  transcribing the multi-term classical closed form (the \"logarithmic
  case\", involving `log(z)` and digamma terms) from the literature."
  ^Vec2 [^Vec2 a ^Vec2 b ^Vec2 z]
  (if (cplx-int? b)
    (let [h 1.0e-5
          bh (cplx/adds b h)
          bh2 (cplx/adds b (m/* 0.5 h))
          f-h (tricomis-U-complex-raw a bh z)
          f-h2 (tricomis-U-complex-raw a bh2 z)]
      (cplx/sub (cplx/scale f-h2 2.0) f-h))
    (tricomis-U-complex-raw a b z)))

(defn- tricomis-U-complex-asymptotic
  "Tricomi's `U(a,b,z)` via the standard asymptotic-series formula
  `z^(-a) pFq([a, 1+a-b], [], -1/z)`, mirroring the real-valued
  [[tricomis-U]]'s own general-case formula. Confirmed against mpmath to
  be reliable across a very wide magnitude range of `z` (0.1 to 300+,
  pure-imaginary `z`, integer `b`), UNLIKE [[tricomis-U-complex-raw]],
  which instead computes a difference of two individually exp(z)-scaled
  terms and loses essentially all precision for `|z|` beyond about 20-30.

  Attempted for ANY `z`, no `Re(z)` restriction -- mirroring mpmath's own
  `hyperu`, which tries this same asymptotic formula unconditionally and
  only falls back to a reflection formula on failure (confirmed by
  reading mpmath's source: two of three previously-documented `Re(z) < 0`
  failures of [[tricomis-U-complex]] turned out to be exactly this --
  the asymptotic formula itself works fine there, it was just never
  tried). Calls the underlying `p > q + 1` pFq's Weniger acceleration
  DIRECTLY (`hypergeometric.hypergeometric-pFq-weniger-complex-with-
  reason`), bypassing the public [[hypergeometric-pFq-complex]]
  dispatcher, so the caller can see whether it actually converged --
  still runs the exact same [[pfq-complex-route]] pre-checks that
  dispatcher would (terminating/pole/degenerate/`z=0`), so no coverage is
  lost by bypassing it. `[a, 1+a-b]`/`[]` always has `p=2 > q+1=1`, so
  `pfq-complex-route` never selects `:maclaurin` here.

  Returns `{:value :reason}`; `:reason` is `:converged` whenever a
  `pfq-complex-route` short-circuit fired (always trustworthy) or the
  Weniger acceleration itself converged, and something else otherwise.
  [[tricomis-U-complex]] falls back to [[tricomis-U-complex-raw-limit]]
  whenever `:reason` isn't `:converged`."
  ^PfQWenigerResultCplx [^Vec2 a ^Vec2 b ^Vec2 z]
  (let [p1 (cplx/sub (cplx/add cplx/ONE a) b)
        zpow (cplx/pow z (cplx/neg a))
        arg (cplx/neg (cplx/reciprocal z))
        ^PfQComplexData res (pfq-complex-route [a p1] [] arg 1048576)]
    (if (.done? res)
      (PfQWenigerResultCplx. (cplx/mult zpow (.value res)) :converged)
      (let [^PfQWenigerResultCplx res2 (hg/hypergeometric-pFq-weniger-complex-with-reason (.ps res) (.qs res) arg 1048576)]
        (PfQWenigerResultCplx. (cplx/mult zpow (.value res2)) (.reason res2))))))

(defn tricomis-U-complex
  "Complex version of Tricomi's confluent hypergeometric function U(a,b,z) of the second kind.

  Arguments `a`, `b` and `z` can be real or complex numbers; plain numbers are promoted to
  complex automatically via `ensure-complex`.

  - Input: `a`, `b`, `z` — real or complex numbers (scalars or [[Vec2]] complex pairs)
  - Returns: [[Vec2]] complex number

  At `z = 0.0+0.0i` with `a` and `b` both real: mirrors the real-valued [[tricomis-U]]'s own `x = 0` limit exactly (finite for `a` a non-positive integer or `b < 1.0`, a signed `##Inf`-valued complex number otherwise). At `z = 0.0+0.0i` with `a` or `b` genuinely complex, the limit is path-dependent (branch-cut sensitive) and not resolved here; `(Vec2. ##NaN ##NaN)` is returned instead.

  For `z != 0.0+0.0i`: tries the asymptotic-series formula first, for ANY `z` (see [[tricomis-U-complex-asymptotic]] -- this mirrors mpmath's own `hyperu`, which does the same). This alone resolves every case checked so far except `Re(z) < 0.0` at moderate-to-large `|z|` with `b` also an integer, where it falls back to a Kummer-`M`-function reflection formula instead. That reflection formula remains accurate for small to moderate `|z|`, but its accuracy becomes unpredictable for larger `|z|` there -- confirmed to occasionally lose several or more digits, or return a non-finite value, for specific parameter combinations even at moderate `|z|` (around 5-10), with no clean magnitude threshold separating safe from unsafe cases; this is a known, uncorrected limitation, inherited from the same underlying resummation instability documented for [[hypergeometric-pFq-complex]], and confirmed present in mpmath's own `hyperu` at at least one such point too (it recovers there only by raising its working precision arbitrarily, an option not available at fixed double precision). The reflection formula also has its own removable singularity whenever `b` is an integer, resolved via a small Richardson-extrapolated numerical limit.

  See also the real-valued [[tricomis-U]], [[hypergeometric-pFq-complex]]."
  [a b z]
  (let [a (cplx/ensure-complex a)
        b (cplx/ensure-complex b)
        z (cplx/ensure-complex z)]
    (cond
      (and (m/zero? (cplx/re z)) (m/zero? (cplx/im z))
           (m/zero? (cplx/im a)) (m/zero? (cplx/im b)))
      (cplx/complex (tricomis-U (cplx/re a) (cplx/re b) 0.0) 0.0)

      (and (m/zero? (cplx/re z)) (m/zero? (cplx/im z))) (Vec2. ##NaN ##NaN)

      :else
      (let [^PfQWenigerResultCplx res (tricomis-U-complex-asymptotic a b z)]
        (if (= (.reason res) :converged) (.value res) (tricomis-U-complex-raw-limit a b z))))))

(def ^:private CPLX_HALF_PI (cplx/complex m/HALF_PI 0.0))

(defn bessel-K-half-odd-complex
  "Bessel K_a function for a = order/2 for complex numbers

  Function accepts only odd integers for order"
  ^Vec2 [^long odd-numerator ^Vec2 x]
  (case (int odd-numerator)
    1 (cplx/mult (cplx/sqrt (cplx/div CPLX_HALF_PI x)) (cplx/exp (cplx/neg x)))
    3 (cplx/mult (cplx/mult (cplx/sqrt (cplx/div CPLX_HALF_PI x)) (cplx/exp (cplx/neg x)))
                 (cplx/add (cplx/reciprocal x) cplx/ONE))
    (loop [i (long 5)
           [b1 b3] (let [b1 (cplx/mult (cplx/sqrt (cplx/div CPLX_HALF_PI x)) (cplx/exp (cplx/neg x)))
                         b3 (cplx/mult b1 (cplx/add (cplx/reciprocal x) cplx/ONE))]
                     [b1 b3])]
      (if (m/> i odd-numerator)
        b3
        (recur (m/+ i 2) [b3 (cplx/add b1 (cplx/mult b3 (cplx/div (Vec2. (m/- i 2.0) 0.0) x)))])))))

(defn bessel-K-half-odd-scaled-complex
  "Bessel K_a function scaled by e^x for a = order/2 for complex numbers

  Function accepts only odd integers for order"
  [^long odd-numerator x]
  (case (int odd-numerator)
    1 (cplx/sqrt (cplx/div CPLX_HALF_PI x))
    3 (cplx/mult (cplx/sqrt (cplx/div CPLX_HALF_PI x))
                 (cplx/add (cplx/reciprocal x) cplx/ONE))
    (loop [i (long 5)
           [b1 b3] (let [b1 (cplx/sqrt (cplx/div CPLX_HALF_PI x))
                         b3 (cplx/mult b1 (cplx/add (cplx/reciprocal x) cplx/ONE))]
                     [b1 b3])]
      (if (m/> i odd-numerator)
        b3
        (recur (m/+ i 2) [b3 (cplx/add b1 (cplx/mult b3 (cplx/scale (cplx/reciprocal x) (m/- i 2.0))))])))))

;; elliptic

(defn elliptic-K
  "Elliptic K - complete (K) and incomplete (F) elliptic integral of the first kind."
  (^double [^double m] (ellip/K m))
  (^double [^double phi ^double m] (ellip/K phi m)))

(defn elliptic-F
  "Elliptic F - incomplete elliptic integral of the first kind."
  ^double [^double phi ^double m] (ellip/K phi m))

(defn elliptic-E
  "Elliptic E - complete and incomplete elliptic integral of the second kind"
  (^double [^double m] (ellip/E m))
  (^double [^double phi ^double m] (ellip/E phi m)))

(defn elliptic-PI
  "Elliptic PI - complete and incomplete elliptic integral of the third kind"
  (^double [^double n ^double m] (ellip/PI n m))
  (^double [^double n ^double phi ^double m] (ellip/PI n phi m)))

(defn elliptic-D
  "Elliptic D - complete and incomplete elliptic integral of Legendre’s type"
  (^double [^double x]
   (m// (m/- (ellip/K x) (ellip/E x))
        (m/* x x)))
  (^double [^double phi ^double x]
   (m// (m/- (ellip/K phi x) (ellip/E phi x))
        (m/* x x))))

(defn elliptic-Rf
  "Symmetric Rf elliptic intergral of the first kind."
  ^double [^double x ^double y ^double z]
  (ellip/Rf x y z))

(defn elliptic-Rd
  "Symmetric Rd elliptic intergral, symmetry on two variables."
  ^double [^double x ^double y ^double z]
  (ellip/Rd x y z))

(defn elliptic-Rg
  "Symmetric Rg elliptic intergral of the second kind"
  ^double [^double x ^double y ^double z]
  (ellip/Rg x y z))

(defn elliptic-Rj
  "Symmetric Rj elliptic intergral of the third kind"
  ^double [^double x ^double y ^double z ^double p]
  (ellip/Rj x y z p))

(defn elliptic-Rc
  "Rc elliptic intergral"
  ^double [^double x ^double y]
  (ellip/Rc x y))

;; Jacobi

(defn jacobi-am
  "Amplitude phi=am(u,m) such that u=F(phi,m)

  Inverse of the elliptic incomplete intergral of the first kind."
  ^double [^double u ^double m]
  (ellip/am u m))

;; Jacobi am, sn, cn, dn, sc, sd, cs, cd, ds, dc, ns, nc, nd

(defn jacobi-sn "Jacobi sn(u,m)" ^double [^double u ^double k] (ellip/jsn u k))
(defn jacobi-cn "Jacobi cn(u,m)" ^double [^double u ^double k] (ellip/jcn u k))
(defn jacobi-dn "Jacobi dn(u,m)" ^double [^double u ^double k] (ellip/jdn u k))
(defn jacobi-sc "Jacobi sc(u,m)" ^double [^double u ^double k] (ellip/jsc u k))
(defn jacobi-sd "Jacobi sd(u,m)" ^double [^double u ^double k] (ellip/jsd u k))
(defn jacobi-cs "Jacobi cs(u,m)" ^double [^double u ^double k] (ellip/jcs u k))
(defn jacobi-cd "Jacobi cd(u,m)" ^double [^double u ^double k] (ellip/jcd u k))
(defn jacobi-ds "Jacobi ds(u,m)" ^double [^double u ^double k] (ellip/jds u k))
(defn jacobi-dc "Jacobi dc(u,m)" ^double [^double u ^double k] (ellip/jdc u k))
(defn jacobi-ns "Jacobi ns(u,m)" ^double [^double u ^double k] (ellip/jns u k))
(defn jacobi-nc "Jacobi nc(u,m)" ^double [^double u ^double k] (ellip/jnc u k))
(defn jacobi-nd "Jacobi nd(u,m)" ^double [^double u ^double k] (ellip/jnd u k))

;; Inverse of Jacobi am, sn, cn, dn, sc, sd, cs, cd, ds, dc, ns, nc, nd

(defn jacobi-asn "Jacobi arcsn(u,m)" ^double [^double x ^double k] (ellip/jasn x k))
(defn jacobi-acn "Jacobi arccn(u,m)" ^double [^double x ^double k] (ellip/jacn x k))
(defn jacobi-adn "Jacobi arcdn(u,m)" ^double [^double x ^double k] (ellip/jadn x k))
(defn jacobi-asc "Jacobi arcsc(u,m)" ^double [^double x ^double k] (ellip/jasc x k))
(defn jacobi-asd "Jacobi arcsd(u,m)" ^double [^double x ^double k] (ellip/jasd x k))
(defn jacobi-acs "Jacobi arccs(u,m)" ^double [^double x ^double k] (ellip/jacs x k))
(defn jacobi-acd "Jacobi arccd(u,m)" ^double [^double x ^double k] (ellip/jacd x k))
(defn jacobi-ads "Jacobi arcds(u,m)" ^double [^double x ^double k] (ellip/jads x k))
(defn jacobi-adc "Jacobi arcdc(u,m)" ^double [^double x ^double k] (ellip/jadc x k))
(defn jacobi-ans "Jacobi arcns(u,m)" ^double [^double x ^double k] (ellip/jans x k))
(defn jacobi-anc "Jacobi arcnc(u,m)" ^double [^double x ^double k] (ellip/janc x k))
(defn jacobi-and "Jacobi arcnd(u,m)" ^double [^double x ^double k] (ellip/jand x k))
