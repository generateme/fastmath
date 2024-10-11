(ns fastmath.special
  "Special functions for real arguments and value.

  * Bessel J, Y, jinc
  * Modified Bessel I, K 
  * Spherical Bessel j, y
  * Modified spherical Bessel i1, i2, k
  * Gamma, log, digamma, trigamma, polygamma, regularized, lower/upper incomplete
  * Beta, log, regularized, incomplete
  * Erf, inverse
  * Airy A, B with derivatives
  * Zeta (Riemann, Hurwitz), Eta (Dirichlet), Xi (Landau), Beta (Dirichlet)
  * Integrals: Si, Ci, li/Li, Ei, En, Ein
  * Hypergeometric 0F0, 0F1, 1F0, 1F1, 2F1, 2F0, 0F2, pFq, Kummers M, Tricomis U, Whittaker M and W
  * Lambert W (0 and -1)
  * Minkowski
  * Harmonic H"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special.poly :as spoly]
            [fastmath.special.airy :as airy]
            [fastmath.special.hypergeometric :as hg]
            [fastmath.polynomials :as poly]
            [fastmath.complex :as cplx])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]
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

(defn log-beta
  "Logarithm of Beta function."
  {:inline (fn ([p q] `(. Beta (logBeta (double ~p) (double ~q)))))
   :inline-arities #{2}}
  ^double [^double p ^double q] (. Beta (logBeta p q)))

(defn beta
  "Beta function"
  ^double [^double p ^double q]
  (if (and (m/pos? p) (m/pos? q))
    (m/exp (. Beta (logBeta p q)))
    (m// (m/* (gamma p) (gamma q)) (gamma (m/+ p q)))))

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

(defn zeta
  "Riemann and Hurwitz zeta functions for real arguments"
  (^double [^double s]
   (cond
     (m/zero? s) -0.5
     (or (m/one? s) (m/nan? s) (m/neg-inf? s)) ##NaN
     (m/pos-inf? s) 1.0
     (m/< (m/abs s) 1.0e-3) (poly/mevalpoly s -0.5,
                                            -0.918938533204672741780329736405617639861,
                                            -1.0031782279542924256050500133649802190,
                                            -1.00078519447704240796017680222772921424,
                                            -0.9998792995005711649578008136558752359121)
     (m/< s 0.5) (let [oms (m/- 1.0 s)]
                   (m/* (zeta oms) (gamma oms) (m/sin (m/* m/HALF_PI s)) (m/pow m/TWO_PI s) m/INV_PI))
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
     (m/== s 2.0) (trigamma z)
     (or (m/nan? s) (m/nan? z) (m/neg-inf? s)
         (and (m/inf? s) (m/neg? z))) ##NaN
     (and (m/inf? s) (m/> z 1.0)) 0.0
     (m/inf? s) ##Inf
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
  (m/* (m/exp (m/* -1.38629436111989061883 x))
       (m/- (zeta x 0.25) (zeta x 0.75))))

(defn xi
  "Riemann (Landau's) Xi function"
  ^double [^double s]
  (cond
    (m/neg? s) (xi (m/- 1.0 s))
    (m/one? s) 0.5
    (m/zero? s) 0.5
    :else (let [hs (m/* 0.5 s)]
            (m/* hs (m/dec s) (m/pow m/INV_PI hs) (gamma hs) (zeta s)))))

(def ^:private cotderiv-q-memo
  (memoize
   (fn ^doubles [^long m]
     (case (int m)
       0 (double-array [1.0])
       1 (double-array [1.0 1.0])
       (let [^doubles q- (cotderiv-q-memo (m/dec m))
             d (m/dec (alength q-))]
         (if (m/odd? (m/dec m))
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
  "Polygamma function of order `m` and real argument."
  ^double [^long m ^double x]
  (if (m/neg? m)
    ##NaN
    (case (int m)
      0 (digamma x)
      1 (trigamma x)
      (let [s (m/inc m)]
        (if (m/not-pos? x)
          (let [v (cotderiv m x)]
            (m/* (m/+ (zeta s (m/- 1.0 x)) (if (m/even? m) v (m/- v)))
                 (m/- (gamma s))))
          (let [v (m/* (zeta s x) (m/- (gamma s)))]
            (if (m/even? m) v (m/- v))))))))

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
      (m/< x 26.0) (let [n (m/dec (unchecked-int (m/* m/M_2_PI x)))
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
      (m/< x 26.0) (let [n (m/dec (unchecked-int (m/* m/M_2_PI x)))
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
        sgn (if (m/even? abs-v) 1.0 -1.0)
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
    (m/integer? order) (bessel-j-integer-order order x)
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
  ^double [^double order ^double x]
  (let [ao (m/abs order)
        y (bessel-y-positive-args ao x)]
    (if (and (m/neg? order) (m/odd? ao)) (m/- y) y)))

(defn bessel-Y
  "Bessel function of the second kind of order v, Y_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-Y0 x)
    (m/one? order) (bessel-Y1 x)
    (or (m/nan? order) (m/nan? x) (m/neg? x)) ##NaN
    (m/integer? order) (bessel-y-integer-order order x)
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
        (let [i+ (m/inc i)
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
    (m/integer? order) (bessel-i-integer-order (m/abs order) x)
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
  "Sine integral"
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
  "Sine integral, Si shifted by -pi/2"
  ^double [^double x] (m/- (Si x) m/HALF_PI))

(def ^:private ^:const ^{:tag 'double} ci-r0 0.616505485620716233797110404100)
(def ^:private ^:const ^{:tag 'double} ci-r1 3.384180422551186426397851146402)
(def ^:private ^:const ^{:tag 'double} ci-r01 0.6162109375)
(def ^:private ^:const ^{:tag 'double} ci-r02 0.29454812071623379711E-3)
(def ^:private ^:const ^{:tag 'double} ci-r11 3.3837890625)
(def ^:private ^:const ^{:tag 'double} ci-r12 0.39136005118642639785E-3)

(defn Ci
  "Cosine integral"
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
  "Cosine integral, alternative definition"
  ^double [^double x]
  (m/- (m/+ m/GAMMA (m/log x)) (Ci x)))

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

(defn E0
  "Exponential integral E0"
  ^double [^double x]
  (if (m/zero? x) ##Inf (m// (m/exp (m/- x)) x)))

(defn E1
  "Exponential integral E1 for positive real numbers"
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
  "Exponential integral, alternative definition"
  ^double [^double x]
  (if (m/zero? x)
    0.0
    (m/+ (E1 x) (m/log x) m/GAMMA)))

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
  (let [gamma-term (m/* (en-safe-expfact (m/dec v) x)
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
            n (m/dec (m/round v))
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
  "Generalized exponential integral En"
  ^double [^double n ^double x]
  (cond
    (m/zero? n) (E0 x)
    (m/one? n) (E1 x)
    (and (m/zero? x) (m/neg? n)) ##Inf
    (m/zero? x) (m// 1.0 (m/dec n))
    (or (m/nan? n) (m/nan? x)
        (if (m/integer? n)
          (and (m/neg? x) (m/pos? n))
          (m/neg? x))) ##NaN
    (m/> x 745.0) 0.0
    (m/< (m/sq x) 9.0) (if (and (m/integer? n) (m/pos? n))
                         (en-expand-origin-posint n x)
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
  "Exponential integral"
  ^double [^double x]
  (cond
    (m/neg? x) (m/- (E1 (m/- x)))
    (m/zero? x) ##-Inf
    (m/> x 710.0) ##Inf
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
                        (m/* (m/exp x) xinv (poly/mevalpoly xinv 1,1,2,6,24,120,720,5040)))))
            (let [dx (m/- x 0.37250741078136663446)]
              (if (m/< (m/abs dx) 0.3)
                (m/* dx (poly/mevalpoly dx 3.896215733907167310, -3.281607866398561671, 6.52237614543892570, -12.96969738353651704, 27.88629796294204998, -62.3788015289154187, 143.5349488096750988, -337.155827178746892, 804.531839982138251, -1943.79664572349884, 4743.76565040243084, -11673.46399116716364, 28926.9553054354509))
                (cond
                  (m/> x 0.6) (ei-taylor64 x 37)
                  (m/> x 0.053) (ei-taylor64 x 15)
                  (m/> x 4.4e-3) (ei-taylor64 x 8)
                  :else (ei-taylor64 x 4)))))))

(defn li
  "Logarythmic integral"
  ^double [^double x]
  (Ei (m/ln x)))

(defn Li
  "Offset logarythmic integral"
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
    (or (m/== a b -1.0)
        (and (m/neg? b)
             (m/integer? a)
             (m/integer? b)
             (or (m/pos? a)
                 (and (m/neg? a) (m/< a b))))) ##NaN
    (m/near-zero? a (m/ulp a)) 1.0
    (m/zero? b) (m/copy-sign ##Inf (m/* a x))
    (m/near-zero? x m/MACHINE-EPSILON) 1.0
    (m/== a b) (m/exp x)
    (m/== a -1.0) (m/- 1.0 (m// x b))
    (and (m/one? a) (m/== b 2.0)) (let [hx (m/* 0.5 x)]
                                    (m/* (m// (m/exp hx) hx) (m/sinh hx)))
    (m/pos? x) (loop [i (long 1)
                      s0 1.0
                      s1 (m/inc (m// (m/* a x) b))]
                 (if (or (and (m/valid-double? s0) (m/valid-double? s1)
                              (m/delta-eq s0 s1 m/MACHINE-EPSILON m/MACHINE-EPSILON))
                         (m/== i 1000000))
                   s1
                   (let [rj (m// (m/* (m/+ a i) x) (m/* (m/+ b i) (m/inc i)))]
                     (recur (m/inc i) s1 (m/+ s1 (m/* (m/- s1 s0) rj))))))
    :else (hg/weniger-1F1 a b x)))

(defn whittaker-M
  "Whittaker's M"
  ^double [^double kappa ^double mu ^double x]
  (let [mu+05 (m/+ 0.5 mu)
        z (m/exp (m/* 0.5 (m/+ (m/* -0.5 x) (m/* mu+05 (m/log x)))))]
    (m/* z (kummers-M (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x) z)))

(defn hypergeometric-0F0
  "Hypergeometric ₀F₀ function, exp(x)"
  ^double [^double x] (m/exp x))

(defn hypergeometric-1F0
  "Hypergeometric ₁F₀ function."
  ^double [^double a ^double x]
  (m/pow (m/- 1.0 x) (m/- a)))

(defn hypergeometric-0F1
  "Confluent hypergeometric ₀F₁ limit function."
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
  "Confluent hypergeometric ₁F₁ function of the first kind, Kummer's M."
  ^double [^double a ^double b ^double x]
  (kummers-M a b x))

(defn hypergeometric-0F2
  "Generalized hypergeometric ₀F₂ function."
  ^double [^double a ^double b ^double x]
  (if (m/pos? x)
    (hg/maclaurin-0F2 a b x)
    (hg/weniger-0F2 a b x)))

(defn hypergeometric-2F0
  "Generalized hypergeometric ₂F₀ function."
  ^double [^double a ^double b ^double x]
  (hg/weniger-2F0 a b x))

(defn tricomis-U
  "Confluent hypergeometric function U of the second kind."
  ^double [^double a ^double b ^double x]
  (cond
    ;; wolfram alpha
    (m/zero? x) (if (m/< b 1.0)
                  (m// (gamma (m/- 1.0 b)) (gamma (m/inc (m/- a b))))
                  ##NaN)
    (m/== a b) (m/* (m/exp x) (upper-incomplete-gamma (m/- 1.0 a) x))
    (m/== a (m/dec b)) (m/pow x (m/- a))    
    :else (m/* (m/pow x (m/- a)) (hg/weniger-2F0 a (m/inc (m/- a b)) (m/- (m// x))))))

(defn whittaker-W
  "Whittaker's W"
  ^double [^double kappa ^double mu ^double x]
  (let [mu+05 (m/+ 0.5 mu)
        z (m/exp (m/* 0.5 (m/+ (m/* -0.5 x) (m/* mu+05 (m/log x)))))]
    (m/* z (tricomis-U (m/- mu+05 kappa) (m/inc (m/* 2.0 mu)) x) z)))

(defn hypergeometric-2F1
  "Gauss's hypergeometric ₂F₁ function."
  ^double [^double a ^double b ^double c ^double x]
  (hg/hypergeometric-2F1 a b c x))

;;

(defn- reduce-m-pFq
  ^double [^doubles a ^long l ^long k ^double r]
  (loop [i (long 0) r r]
    (if (m/== i l) r
        (recur (m/inc i) (m/* r (m/+ (Array/aget a i) k))))))

(defn- reduce-d-pFq
  ^double [^doubles a ^long l ^long k ^double r]
  (loop [i (long 0) r r]
    (if (m/== i l) r
        (recur (m/inc i) (m// r (m/+ (Array/aget a i) k))))))

(defn hypergeometric-pFq
  "hypergeometric-pFq using MacLaurin series.

  `max-iters` is set to 1000000 by default."
  (^double [ps qs ^double z] (hypergeometric-pFq ps qs z 1000000))
  (^double [ps qs ^double z ^long max-iters]
   (let [a (double-array ps)
         la (alength a)
         b (double-array qs)
         lb (alength b)]
     (loop [k (long 1)
            s0 1.0
            s1 (m/inc (m// (m/* z (v/prod a)) (v/prod b)))]
       (if (and (m/< k max-iters)
                (not (m/delta-eq s0 s1 m/MACHINE-EPSILON10 m/MACHINE-EPSILON10)))
         (let [rk (m// z (m/inc k)) 
               rk (reduce-m-pFq a la k rk)
               rk (reduce-d-pFq b lb k rk)]
           (recur (m/inc k) s1 (m/+ s1 (m/* (m/- s1 s0) rk))))
         s1)))))

(set! *unchecked-math* true)

(defn hypergeometric-pFq-ratio
  "Hypergeometric-pFq using MacLaurin series on ratios. Can be very slow.

  `max-iters` is set to 10000 by default."
  (^double [ps qs z] (hypergeometric-pFq-ratio ps qs z 10000))
  (^double [ps qs z ^long max-iters]
   (let [a (map rationalize ps)
         b (map rationalize qs)
         z (rationalize z)
         eps (rationalize m/MACHINE-EPSILON10)]
     (loop [k (long 1)
            s0 (rationalize 1)
            s1 (+ 1 (/ (* z (reduce * a)) (reduce * b)))]
       (if (and (m/< k max-iters)
                (> (abs (- s1 s0)) (max eps (* eps (max s1 s0)))))
         (let [rk (/ z (+ k 1))
               rk (reduce (fn [r va] (* r (+ va k))) rk a)
               rk (reduce (fn [r vb] (/ r (+ vb k))) rk b)]
           (recur (m/inc k) s1 (+ s1 (* (- s1 s0) rk))))
         s1)))))
