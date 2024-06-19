(ns fastmath.special
  "Bessel functions implemented directly in Clojure.

  Relative errors comparing to R are not worse than:
  
  K: 1.0e-12
  I: 1.0e-14"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special.poly :as spoly]
            [fastmath.complex :as cplx])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; https://github.com/JuliaMath/Bessels.jl/blob/master/src/BesselFunctions/besselk.jl

;; Bessel J

(defn bessel-j0
  "Bessel function of the first kind of order 0, J_0(x)"
  ^double [^double x]
  (let [x (m/abs x)]
    (cond
      (m/zero? x) 1.0
      (m/< x m/HALF_PI) (m/mevalpoly (m/* x x)
                          1.0, -0.25, 0.01562499999999994, -0.00043402777777725544, 6.781684026082576e-6,
                          -6.781683757550061e-8, 4.709479394601058e-10, -2.4016837144506874e-12,
                          9.104258208703104e-15)
      (m/< x 26.0) (let [n (m/dec (unchecked-int (m/* m/M_2_PI x)))
                         ^Vec2 root (spoly/j0-roots n)]
                     (spoly/j0-polys n (m/- x (.x root) (.y root))))
      (m/pos-inf? x) 0.0
      :else (let [xinv (m// x)
                  x2 (m/* xinv xinv)
                  ^Vec2 pq (if (m/< x 125)
                             (Vec2. (m/mevalpoly x2 1.0 -0.0625, 0.103515625, -0.5428466796875,
                                                 5.848699569702148, -106.8867939710617, 2968.142937842757,
                                                 -116538.4796968361)
                                    (m/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518
                                                 -23.47512774997287 535.640519510616 -17837.27968894748))
                             (Vec2. (m/mevalpoly x2 1.0 -0.0625, 0.103515625, -0.5428466796875)
                                    (m/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518)))
                  a (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq))
                  xn (m/* xinv (.y pq))
                  b (m/sin (m/+ x m/QUARTER_PI xn))]
              (m/* a b)))))

(defn bessel-j1
  "Bessel function of the first kind of order 1, J_1(x)"
  ^double [^double x]
  (let [s (m/sgn x)
        x (m/abs x)]
    (cond
      (m/zero? x) 0.0
      (m/<= x m/HALF_PI) (m/* s x
                              (m/mevalpoly (m/* x x)
                                0.5, -0.0624999999999989, 0.002604166666657291, -5.42534721917933e-5,
                                6.781683542660179e-7, -5.651361336587487e-9, 3.36191211106159e-11,
                                -1.4511302591871352e-13))
      (m/< x 26.0) (let [n (m/dec (unchecked-int (m/* m/M_2_PI x)))
                         ^Vec2 root (spoly/j1-roots n)]
                     (m/* s (spoly/j1-polys n (m/- x (.x root) (.y root)))))
      (m/pos-inf? x) 0.0
      :else (let [xinv (m// x)
                  x2 (m/* xinv xinv)
                  ^Vec2 pq (if (m/< x 125)
                             (Vec2. (m/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625 -7.739953994750977
                                                 132.7618242502213 -3543.303665366024 135394.2285691809)
                                    (m/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054
                                                 30.6240119934082 -659.185221823779 21156.31404552781))
                             (Vec2. (m/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625)
                                    (m/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054)))
                  a (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq))
                  xn (m/* xinv (.y pq))
                  b (m/sin (m/+ x m/-QUARTER_PI xn))]
              (m/* s a b)))))

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
           a (m// (Math/pow hx v) (m/gamma (m/inc v)))]
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
        (m/> v (m/mevalpoly x 2.0 0.109 0.062))) (bessel-j-power-series v x) 
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

(defn bessel-j
  "Bessel function of the first kind of order v, J_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-j0 x)
    (m/one? order) (bessel-j1 x)
    (m/invalid-double? x) x
    (m/== order (m/round order)) (bessel-j-integer-order order x)
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

(defn bessel-y0
  "Bessel function of the second kind of order 0, Y_0(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##-Inf
    (m/pos-inf? x) 0.0
    (or (m/neg? x) (m/invalid-double? x)) ##NaN
    (m/< x 5.0) (let [z (m/* x x)
                      w (m// (m/mevalpoly z -1.84950800436986690637E16, 4.42733268572569800351E16,
                                          -3.46628303384729719441E15, 8.75906394395366999549E13,
                                          -9.82136065717911466409E11, 5.43526477051876500413E9,
                                          -1.46639295903971606143E7, 1.55924367855235737965E4)
                             (m/mevalpoly z 2.50596256172653059228E17, 3.17157752842975028269E15,
                                          2.02979612750105546709E13, 8.64002487103935000337E10,
                                          2.68919633393814121987E8, 6.26107330137134956842E5,
                                          1.04128353664259848412E3, 1.00000000000000000000E0))]
                  (m/+ w (m/* m/M_2_PI (m/log x) (bessel-j0 x))))
    (m/< x 25.0) (let [w (m// 5.0 x)
                       z (m/* w w)
                       p (m// (m/mevalpoly z 9.99999999999999997821E-1, 5.30324038235394892183E0,
                                           8.74716500199817011941E0, 5.44725003058768775090E0,
                                           1.23953371646414299388E0, 8.28352392107440799803E-2,
                                           7.96936729297347051624E-4)
                              (m/mevalpoly z 1.00000000000000000218E0, 5.30605288235394617618E0,
                                           8.76190883237069594232E0, 5.47097740330417105182E0,
                                           1.25352743901058953537E0, 8.56288474354474431428E-2,
                                           9.24408810558863637013E-4))
                       q (m// (m/mevalpoly z -6.05014350600728481186E0, -5.14105326766599330220E1,
                                           -1.47077505154951170175E2, -1.77681167980488050595E2,
                                           -9.32060152123768231369E1, -1.95539544257735972385E1,
                                           -1.28252718670509318512E0, -1.13663838898469149931E-2)
                              (m/mevalpoly z 2.42005740240291393179E2, 2.06209331660327847417E3,
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
                           (Vec2. (m/mevalpoly x2 1.0 -0.0625 0.103515625 -0.5428466796875 5.848699569702148
                                               -106.8867939710617 2968.142937842757 -116538.4796968361)
                                  (m/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518
                                               -23.47512774997287 535.640519510616 -17837.27968894748))
                           (Vec2. (m/mevalpoly x2 1.0 -0.0625 0.103515625 -0.5428466796875)
                                  (m/mevalpoly x2 -0.125 0.06510416666666667 -0.2095703125 1.638065883091518)))]
            (m/* m/SQRT_2_PI (m/sqrt xinv) (.x pq)
                 (m/sin (m/+ x m/-QUARTER_PI (m/* xinv (.y pq))))))))

(defn bessel-y1
  "Bessel function of the second kind of order 1, Y_1(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##-Inf
    (m/pos-inf? x) 0.0
    (or (m/neg? x) (m/invalid-double? x)) ##NaN
    (m/< x 5.0) (let [z (m/* x x)
                      w (* x (m// (m/mevalpoly z -7.78877196265950026825E17, 2.02439475713594898196E17,
                                               -8.12770255501325109621E15, 1.14509511541823727583E14,
                                               -6.47355876379160291031E11, 1.26320474790178026440E9)
                                  (m/mevalpoly z 3.97270608116560655612E18, 6.87141087355300489866E16,
                                               6.20557727146953693363E14, 3.88231277496238566008E12,
                                               1.87601316108706159478E10, 7.34811944459721705660E7,
                                               2.35564092943068577943E5, 5.94301592346128195359E2,
                                               1.00000000000000000000E0)))]
                  (m/+ w (m/* m/M_2_PI (m/- (m/* (bessel-j1 x) (m/log x)) (m// x)))))
    (m/< x 25.0) (let [w (m// 5.0 x)
                       z (m/* w w)
                       p (m// (m/mevalpoly z 1.00000000000000000254E0, 5.21451598682361504063E0,
                                           8.42404590141772420927E0, 5.11207951146807644818E0,
                                           1.12719608129684925192E0, 7.31397056940917570436E-2,
                                           7.62125616208173112003E-4)
                              (m/mevalpoly z 9.99999999999999997461E-1, 5.20982848682361821619E0,
                                           8.39985554327604159757E0, 5.07386386128601488557E0,
                                           1.10514232634061696926E0, 6.88455908754495404082E-2,
                                           5.71323128072548699714E-4))
                       q (m// (m/mevalpoly z 2.52070205858023719784E1, 2.11688757100572135698E2,
                                           5.97489612400613639965E2, 7.10856304998926107277E2,
                                           3.66779609360150777800E2, 7.58238284132545283818E1,
                                           4.98213872951233449420E0, 5.10862594750176621635E-2)
                              (m/mevalpoly z 3.36093607810698293419E2, 2.82619278517639096600E3,
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
                           (Vec2. (m/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625 -7.739953994750977
                                               132.7618242502213 -3543.303665366024 135394.2285691809)
                                  (m/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054
                                               30.6240119934082 -659.185221823779 21156.31404552781))
                           (Vec2. (m/mevalpoly x2 1.0 0.1875 -0.193359375 0.8052978515625)
                                  (m/mevalpoly x2 0.375 -0.1640625 0.3708984375 -2.369397844587054)))]
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
               a (m// a (m/gamma (m/inc v)))
               b (m// b (m/gamma (m/- 1.0 v)))]
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
    (and (m/== v (m/round v))
         (m/< v 250)) (.x ^Vec2 (bessel-j-up-recurrence x (Vec2. (bessel-y1 x) (bessel-y0 x))
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

(defn bessel-y
  "Bessel function of the second kind of order v, Y_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-y0 x)
    (m/one? order) (bessel-y1 x)
    (or (m/nan? order) (m/nan? x) (m/neg? x)) ##NaN
    (m/== order (m/round order)) (bessel-y-integer-order order x)
    (m/not-neg? order) (bessel-y-positive-args (m/abs order) x)
    :else (let [ao (m/abs order)
                aopi (m/* ao m/PI)
                y (bessel-y-positive-args (m/abs order) x)
                j (bessel-j-positive-args (m/abs order) x)]
            (m/+ (m/* y (m/cos aopi))
                 (m/* j (m/sin aopi))))))

;; Bessel K

(defn bessel-k0
  "Modified Bessel function of the second kind of order 0, K_0(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##Inf
    (m/pos-inf? x) 0.0
    (or (m/nan? x) (m/neg? x)) ##NaN
    (m/<= x 1.0) (let [x2 (m/* x x)
                       a (m/* 0.25 x2)
                       s (m/muladd
                          (m/mevalpoly a -1.372509002685546267e-1, 2.574916117833312855e-1,
                                       1.395474602146869316e-2, 5.445476986653926759e-4,
                                       7.125159422136622118e-6)
                          (m// (m/mevalpoly a 1.000000000000000000e+00, -5.458333438017788530e-02,
                                            1.291052816975251298e-03, -1.367653946978586591e-05))
                          1.137250900268554688)
                       a (m/muladd s a 1.0)]
                   (m/muladd (m/- a)
                             (m/log x)
                             (m/mevalpoly x2 1.159315156584124484e-01, 2.789828789146031732e-01,
                                          2.524892993216121934e-02, 8.460350907213637784e-04,
                                          1.491471924309617534e-05, 1.627106892422088488e-07,
                                          1.208266102392756055e-09, 6.611686391749704310e-12)))
    :else (let [rx (m// x)
                a (m/muladd
                   (m/mevalpoly rx 2.533141373155002416e-1, 3.628342133984595192e0,
                                1.868441889406606057e1, 4.306243981063412784e1,
                                4.424116209627428189e1, 1.562095339356220468e1,
                                -1.810138978229410898e0, -1.414237994269995877e0,
                                -9.369168119754924625e-2)
                   (m// (m/mevalpoly rx 1.000000000000000000e0, 1.494194694879908328e1,
                                     8.265296455388554217e1, 2.162779506621866970e2,
                                     2.845145155184222157e2, 1.851714491916334995e2,
                                     5.486540717439723515e1, 6.118075837628957015e0,
                                     1.586261269326235053e-1))
                   1.0)]
            (m/* (m// a (m/sqrt x)) (m/exp (m/- x))))))

(defn bessel-k1
  "Modified Bessel function of the second kind of order 1, K_1(x)"
  ^double [^double x]
  (cond
    (m/zero? x) ##Inf
    (m/pos-inf? x) 0.0
    (or (m/nan? x) (m/neg? x)) ##NaN
    (m/<= x 1.0) (let [x2 (m/* x x)
                       a (m/* 0.25 x2)
                       pq (m/muladd
                           (m/mevalpoly a -3.62137953440350228e-3, 7.11842087490330300e-3,
                                        1.00302560256614306e-5, 1.77231085381040811e-6)
                           (m// (m/mevalpoly a 1.00000000000000000e0, -4.80414794429043831e-2,
                                             9.85972641934416525e-4, -8.91196859397070326e-6))
                           8.69547128677368164e-2)
                       pq (m/muladd (m/* pq a) a (m/inc (m/* a 0.5)))
                       a (m/* pq x 0.5)
                       pq (m/muladd
                           (m// (m/mevalpoly x2 -3.07965757829206184e-1, -7.80929703673074907e-02,
                                             -2.70619343754051620e-3, -2.49549522229072008e-5)
                                (m/mevalpoly x2 1.00000000000000000e0, -2.36316836412163098e-2,
                                             2.64524577525962719e-4, -1.49749618004162787e-6))
                           x (m// x))]
                   (m/muladd a (m/log x) pq))
    :else (let [rx (m// x)
                a (m/muladd
                   (m/mevalpoly rx -1.97028041029226295e-1, -2.32408961548087617e0,
                                -7.98269784507699938e0, -2.39968410774221632e0,
                                3.28314043780858713e1, 5.67713761158496058e1,
                                3.30907788466509823e1, 6.62582288933739787e0,
                                3.08851840645286691e-1)
                   (m// (m/mevalpoly rx 1.00000000000000000e0, 1.41811409298826118e1,
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
        sp (m/mevalpoly vv 1.0, 1.6449340668482264, 1.8940656589944918, 1.9711021825948702)
        g1 (m/mevalpoly vv -0.5772156649015329, 0.04200263503409518, 0.042197734555544306)
        g2 (m/mevalpoly vv 1.0, -0.6558780715202539, 0.16653861138229145)
        sh (m/mevalpoly (m/* mu mu) 1.0, 0.16666666666666666, 0.008333333333333333, 0.0001984126984126984, 2.7557319223985893e-6)]
    (m/* sp (m/+ (m/* g1 (m/cosh mu))
                 (m/* g2 sh l2dx)))))

(defn- bessel-k-temme-series
  ^Vec2 [^double v ^double x]
  (let [z (m/* x 0.5)
        zz (m/* z z)
        zv (Math/pow z v)]
    (loop [k (long 1)
           fk (f0-local-expansion-v0 v x)
           pk (m// (m/mevalpoly v 1.0, -0.5772156649015329, 0.9890559953279725, -0.23263776388631713)
                   (m/* 2.0 zv))
           qk (m/* (m/mevalpoly (m/- v) 1.0, -0.5772156649015329, 0.9890559953279725, -0.23263776388631713)
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
  (let [gam (m/gamma v)
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

(defn bessel-k
  "Modified Bessel function of the second kind and real order v, K_v(x)"
  ^double [^double order ^double x]
  (let [v (m/abs order)]
    (cond
      (m/zero? v) (bessel-k0 x)
      (m/one? v) (bessel-k1 x)
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
                
                (m/< (m/abs (m/- v (m/round v))) 1.0e-5)
                (let [v-floor (if (m/> v-floor 0.5) (m/dec v-floor) v-floor)
                      kv (bessel-k-temme-series v-floor x)
                      ^Vec2 res (bessel-k-up-recurrence x kv (Vec2. (m/inc v-floor) v))]
                  (.x res))
                
                :else (bessel-k-power-series v x))))))

;; Bessel I

(defn bessel-i0
  "Modified Bessel function of the first kind of order 0, I_0(x)"
  ^double [^double x]
  (cond
    (m/invalid-double? x) ##NaN
    (m/zero? x) 1.0
    :else (let [x (m/abs x)]
            (if (m/< x 7.75)
              (let [a (m/* 0.25 x x)]
                (m/muladd a (m/mevalpoly a  0.9999999999999998, 0.2500000000000052, 0.027777777777755364,
                                         0.001736111111149161, 6.94444444107536e-5, 1.9290123635366806e-6,
                                         3.9367592765038015e-8, 6.151201574092085e-10, 7.593827956729909e-12,
                                         7.596677643342155e-14, 6.255282299620455e-16, 4.470993793303175e-18,
                                         2.1859737023077178e-20, 2.0941557335286373e-22) 1.0))
              (m/* (m/exp x)
                   (m// (m/mevalpoly (m// x) 0.3989422804014326, 0.04986778505064754, 0.028050628512954097,
                                     0.02921968830978531, 0.04466889626137549, 0.10220642174207666,
                                     -0.9937439085650689, 91.25330271974727, -4901.408890977662,
                                     199209.2752981982, -6.181516298413396e6, 1.4830278710991925e8,
                                     -2.7695254643719645e9, 4.0351394830842026e10, -4.5768930327229974e11,
                                     4.0134844243070063e12, -2.6862476523182016e13, 1.3437999451218112e14,
                                     -4.856333741437621e14, 1.1962791200680235e15, -1.796269414464399e15,
                                     1.239942074380968e15)
                        (m/sqrt x)))))))

(defn bessel-i1
  "Modified Bessel function of the first kind of order 1, I_0(x)"
  ^double [^double x]
  (cond
    (m/invalid-double? x) ##NaN
    (m/zero? x) 0.0
    :else (let [z (m/abs x)
                z (double
                   (if (m/< z 7.75)
                     (let [a (m/* 0.25 z z)
                           inner (m/mevalpoly a 0.08333333333333334, 0.006944444444444374,
                                              0.00034722222222248526, 1.1574074073690356e-5,
                                              2.7557319253050506e-7, 4.920949730519126e-9,
                                              6.834656365321179e-11, 7.593985414952446e-13,
                                              6.904652315442046e-15, 5.2213850252454655e-17,
                                              3.405120412140281e-19, 1.6398527256182257e-21,
                                              1.3161876924566675e-23)]
                       (m/* 0.5 z (m/mevalpoly a 1.0 0.5 inner)))
                     (m/* (m/exp z)
                          (m// (m/mevalpoly (m// z) 0.39894228040143276, -0.149603355151029,
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
                               (m/sqrt z)))))]
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
                          (m/gamma (m/inc v))))
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

(defn bessel-i
  "Modified Bessel function of the first kind of order v, I_v(x)"
  ^double [^double order ^double x]
  (cond
    (m/zero? order) (bessel-i0 x)
    (m/one? order) (bessel-i1 x)
    (m/invalid-double? x) x
    (m/== order (m/round order)) (bessel-i-integer-order (m/abs order) x)
    (m/neg? x) ##NaN
    (m/zero? x) 0.0
    (m/not-neg? order) (bessel-i-positive-args (m/abs order) x)
    :else (let [a (m/abs order)]
            (m/+ (bessel-i-positive-args a x)
                 (m/* m/M_2_PI (m/sin (m/* m/PI a)) (bessel-k a x))))))
