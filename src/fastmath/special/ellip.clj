(ns fastmath.special.ellip
  (:require [fastmath.core :as m]
            [fastmath.polynomials :as poly])
  (:import [fastmath.java Array]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; 

;; https://arxiv.org/pdf/math/9409227v1
;; https://dlmf.nist.gov/19.20

(defn Rc-closed
  "Rc - Carlson symmetric form, closed formula."
  ^double [^double x ^double y]
  (if (m/neg? y)
    (let [y- (m/- y)]
      (m/* (m/sqrt (m// x (m/+ x y-)))
           (Rc-closed (m/+ x y-) y-)))
    (let [absdiff (m/> (m/abs (m/- x y)) (m/* 2.0 (m/max (m/ulp x) (m/ulp y))))]
      (cond
        (and absdiff (m/< x y)) (m// (m/acos (m/sqrt (m// x y)))
                                     (m/sqrt (m/- y x)))
        (and absdiff (m/> x y)) (m// (m/acosh (m/sqrt (m// x y)))
                                     (m/sqrt (m/- x y)))
        :else (m// (m/sqrt (m/* 0.5 (m/+ x y))))))))

(def ^{:private true :const true :tag (quote double)} Rc-error (m/pow (m// m/MACHINE-EPSILON 16.0) m/SIXTH))

(defn Rc
  "Rc - Carlson symmetric form, iterative.

  `x>=0`; for `y<0` the Cauchy principal value. Returns `Inf` for `y=0` (the integral diverges), `0.0` for an infinite argument and `NaN` for `x<0` or a `NaN` argument."
  ^double [^double x ^double y]
  (cond
    (or (m/nan? x) (m/nan? y) (m/neg? x)) ##NaN
    (m/== x y) (m// (m/sqrt x))
    (and (m/zero? x) (m/zero? y)) ##Inf
    (and (m/zero? x) (m/pos? y)) (m// m/HALF_PI (m/sqrt y))
    (m/zero? y) ##Inf
    (or (m/inf? x) (m/inf? y)) 0.0
    (m/neg? y) (let [y- (m/- y)]
                 (m/* (m/sqrt (m// x (m/+ x y-)))
                      (Rc (m/+ x y-) y-)))
    :else (let [mu0 (m/* m/THIRD (m/+ x y y))]
            (loop [i (long 0)
                   xn x
                   yn y
                   mu mu0
                   sn (m/- (m// (m/+ y mu0) mu0) 2.0)]
              (cond
                (m/< (m/max (m/abs sn)) Rc-error)
                (m// (m/inc (m/* sn sn (poly/mevalpoly sn 0.3 0.14285714285714285 0.375 0.4090909090909091)))
                     (m/sqrt mu))
                
                (m/> i 1000) (throw (ex-info "Rc didn't converge" {:sn sn}))

                :else (let [lambda (m/+ (m/* 2.0 (m/sqrt xn) (m/sqrt yn)) yn)
                            nxn (m/* 0.25 (m/+ xn lambda))
                            nyn (m/* 0.25 (m/+ yn lambda))
                            nmu (m/* m/THIRD (m/+ nxn nyn nyn))]
                        (recur (m/inc i)
                               nxn nyn nmu
                               (m/- (m// (m/+ nyn nmu) nmu) 2.0))))))))

(def ^{:private true :const true :tag (quote double)} Rd-error (m/pow (m// m/MACHINE-EPSILON 4.0) m/SIXTH))

(defn Rd
  "Rd - Carlson symmetric form of the second kind, iterative. See [[fastmath.special/elliptic-Rd]]."
  ^double [^double x ^double y ^double z]
  (cond
    (or (m/nan? x) (m/nan? y) (m/nan? z)
        (m/neg? x) (m/neg? y) (m/neg? z)) ##NaN
    (or (m/zero? z) (and (m/zero? x) (m/zero? y))) ##Inf
    (or (m/inf? x) (m/inf? y) (m/inf? z)) 0.0
    (and (m/zero? x) (m/== y z)) (m// (m/* 3.0 m/QUARTER_PI)
                                      (m/pow y 1.5))
    (m/== x y z) (m/pow x -1.5)
    :else (let [mu0 (m/* 0.2 (m/+ x y (m/* 3.0 z)))]
            (loop [i (long 0)
                   xn x
                   yn y
                   zn z
                   mu mu0
                   xndev (m// (m/- mu0 x) mu0)
                   yndev (m// (m/- mu0 y) mu0)
                   zndev (m// (m/- mu0 z) mu0)
                   sigma 0.0
                   power4 1.0]
              (cond
                (m/< (m/max (m/abs xndev) (m/abs yndev) (m/abs zndev)) Rd-error)
                (let [ea (m/* xndev yndev)
                      eb (m/* zndev zndev)
                      ec (m/- ea eb)
                      ed (m/- ea (m/* 6.0 eb))
                      ef (m/+ ed ec ec)
                      s1 (m/* ed (m/- (m/* 0.10227272727272728 ed)
                                      (m/* 0.1730769230769231 zndev ef)
                                      0.2142857142857143))
                      s2 (m/* zndev (m/+ (m/* 0.1666666666666667 ef)
                                         (m/* zndev (m/- (m/* 0.1153846153846154 zndev ea)
                                                         (m/* 0.4090909090909091 ec)))))]
                  (m/+ (m/* 3.0 sigma)
                       (m// (m/* power4 (m/+ 1.0 s1 s2))
                            (m/* mu (m/sqrt mu)))))

                (m/> i 1000) (throw (ex-info "Rf didn't converge" {:xndev xndev :yndev yndev :zndev zndev}))

                :else (let [xnroot (m/sqrt xn)
                            ynroot (m/sqrt yn)
                            znroot (m/sqrt zn)
                            lambda (m/+ (m/* xnroot (m/+ ynroot znroot))
                                        (m/* ynroot znroot))
                            nxn (m/* 0.25 (m/+ xn lambda))
                            nyn (m/* 0.25 (m/+ yn lambda))
                            nzn (m/* 0.25 (m/+ zn lambda))
                            nmu (m/* 0.2 (m/+ nxn nyn (m/* 3.0 nzn)))]
                        (recur (m/inc i)
                               nxn nyn nzn nmu
                               (m// (m/- nmu nxn) nmu)
                               (m// (m/- nmu nyn) nmu)
                               (m// (m/- nmu nzn) nmu)
                               (m/+ sigma (m// power4 (m/* znroot (m/+ zn lambda))))
                               (m/* 0.25 power4))))))))

(def ^{:private true :const true :tag (quote double)} Rf-error (m/pow (m/* 3 m/MACHINE-EPSILON) m/SIXTH))

(defn Rf
  "Rf - Carlson symmetric form of the first kind, iterative. See [[fastmath.special/elliptic-Rf]]."
  ^double [^double x ^double y ^double z]
  (cond
    (or (m/nan? x) (m/nan? y) (m/nan? z)
        (m/neg? x) (m/neg? y) (m/neg? z)) ##NaN
    ;; two zero arguments: the integral diverges at t = 0
    (or (and (m/zero? x) (m/zero? y))
        (and (m/zero? z) (or (m/zero? x) (m/zero? y)))) ##Inf
    (or (m/inf? x) (m/inf? y) (m/inf? z)) 0.0
    (m/== y z) (Rc x y)
    :else (let [mu0 (m/* m/THIRD (m/+ x y z))]
            (loop [i (long 0)
                   xn x
                   yn y
                   zn z
                   mu mu0
                   xndev (m/- 2.0 (m// (m/+ mu0 x) mu0))
                   yndev (m/- 2.0 (m// (m/+ mu0 y) mu0))
                   zndev (m/- 2.0 (m// (m/+ mu0 z) mu0))]
              (cond
                (m/< (m/max (m/abs xndev) (m/abs yndev) (m/abs zndev)) Rf-error)
                (let [e2 (m/- (m/* xndev yndev) (m/* zndev zndev))
                      e3 (m/* xndev yndev zndev)]
                  (m// (m/+ 1.0
                            (m/* e2 (m/- (m/* 0.041666666666666664 e2)
                                         (m/* 0.06818181818181818 e3)
                                         0.1))
                            (m/* 0.07142857142857142 e3)) (m/sqrt mu)))

                (m/> i 1000) (throw (ex-info "Rf didn't converge" {:xndev xndev :yndev yndev :zndev zndev}))

                :else (let [xnroot (m/sqrt xn)
                            ynroot (m/sqrt yn)
                            znroot (m/sqrt zn)
                            lambda (m/+ (m/* xnroot (m/+ ynroot znroot))
                                        (m/* ynroot znroot))
                            nxn (m/* 0.25 (m/+ xn lambda))
                            nyn (m/* 0.25 (m/+ yn lambda))
                            nzn (m/* 0.25 (m/+ zn lambda))
                            nmu (m/* m/THIRD (m/+ nxn nyn nzn))]
                        (recur (m/inc i)
                               nxn nyn nzn nmu
                               (m/- 2.0 (m// (m/+ nmu nxn) nmu))
                               (m/- 2.0 (m// (m/+ nmu nyn) nmu))
                               (m/- 2.0 (m// (m/+ nmu nzn) nmu)))))))))

;; largest argument below which Rj and Rg are rescaled to 1 (homogeneity): products of the arguments underflow
(def ^{:private true :const true :tag (quote double)} CARLSON-RESCALE-BELOW 1.0e-50)
(def ^{:private true :const true :tag (quote double)} Rj-error (m/pow (m// m/MACHINE-EPSILON 3.0) m/SIXTH))

(defn- rj-principal-value-overflow?
  "True when a term of the principal value of Rj (`p<0`) overflows: `q-`, `x*z` or `p*q`."
  [^double x ^double y ^double z ^double p]
  (let [q- (m// (m/* (m/- z y) (m/- y x)) (m/- y p))]
    (or (m/inf? q-) (m/inf? (m/* x z)) (m/inf? (m/* p (m/+ y q-))))))

(defn Rj
  "Rj - Carlson symmetric form of the third kind, iterative. See [[fastmath.special/elliptic-Rj]]."
  ^double [^double x ^double y ^double z ^double p]
  (cond
    (or (m/nan? x) (m/nan? y) (m/nan? z) (m/nan? p)
        (m/neg? x) (m/neg? y) (m/neg? z)) ##NaN
    (m/zero? p) ##Inf
    (m/== z p) (Rd x y z)
    (m/== x y z) (Rd p p x)
    (and (m/zero? x) (m/zero? y)) ##Inf
    (or (m/inf? x) (m/inf? y) (m/inf? z) (m/inf? p)) 0.0
    ;; homogeneity Rj(lx,ly,lz,lp) = Rj(x,y,z,p)/l^(3/2), the arguments divided by the largest one
    (or (m/< (m/max x y z (m/abs p)) CARLSON-RESCALE-BELOW)
        (and (m/neg? p) (rj-principal-value-overflow? x y z p)))
    (let [s (m/max x y z (m/abs p))]
      (m// (m// (Rj (m// x s) (m// y s) (m// z s) (m// p s)) s) (m/sqrt s)))
    (m/neg? p) (let [p+ (m/- p)
                     y+p+ (m/+ y p+)
                     q- (m// (m/* (m/- z y) (m/- y x)) y+p+)
                     q (m/+ y q-)]
                 (if (m/pos? q)
                   (m// (m/+ (m/- (m/* q- (Rj x y z q))
                                  (m/* 3.0 (Rf x y z)))
                             (m/* 3.0 (m/sqrt y) (Rc (m/* x z) (m/* p q))))
                        y+p+)
                   ;; y not the median: q <= 0 recursed forever; symmetric in x, y, z, the median as y gives q > 0
                   (Rj (m/min x y z) (m/max (m/min x y) (m/min (m/max x y) z)) (m/max x y z) p)))
    :else (let [mu0 (m/* 0.2 (m/+ x y z p p))]
            (loop [i (long 0)
                   xn x
                   yn y
                   zn z
                   pn p
                   mu mu0
                   xndev (m// (m/- mu0 x) mu0)
                   yndev (m// (m/- mu0 y) mu0)
                   zndev (m// (m/- mu0 z) mu0)
                   pndev (m// (m/- mu0 p) mu0)
                   sigma 0.0
                   power4 1.0]
              (cond
                (m/< (m/max (m/abs xndev) (m/abs yndev)
                            (m/abs zndev) (m/abs pndev)) Rj-error)
                (let [ea (m/+ (m/* xndev (m/+ yndev zndev))
                              (m/* yndev zndev))
                      eb (m/* xndev yndev zndev)
                      ec (m/sq pndev)
                      e2 (m/- ea (m/* 3.0 ec))
                      e3 (m/+ eb (m/* 2.0 pndev (m/- ea ec)))
                      s1 (m/inc (m/* e2 (m/- (m/* 0.1022727272727273 e2)
                                             (m/* 0.1730769230769231 e3)
                                             0.2142857142857143)))
                      s2 (m/* eb (m/+ 0.16666666666666666
                                      (m/* pndev (m/- (m/* 0.1153846153846154 pndev)
                                                      0.2727272727272728))))
                      s3 (m/- (m/* pndev ea (m/- 0.3333333333333333 (m/* 0.1363636363636364 pndev)))
                              (m/* 0.3333333333333333 pndev ec))]
                  (m/+ (m/* 3.0 sigma)
                       (m// (m/* power4 (m/+ s1 s2 s3))
                            (m/* mu (m/sqrt mu)))))
                
                (m/> i 1000) (throw (ex-info "Rj didn't converge" {:xndev xndev :yndev yndev :zndev zndev}))

                :else (let [xnroot (m/sqrt xn)
                            ynroot (m/sqrt yn)
                            znroot (m/sqrt zn)
                            lambda (m/+ (m/* xnroot (m/+ ynroot znroot))
                                        (m/* ynroot znroot))
                            alpha (m/sq (m/+ (m/* pn (m/+ xnroot ynroot znroot))
                                             (m/* xnroot ynroot znroot)))
                            beta (m/* pn (m/sq (m/+ pn lambda)))
                            
                            nxn (m/* 0.25 (m/+ xn lambda))
                            nyn (m/* 0.25 (m/+ yn lambda))
                            nzn (m/* 0.25 (m/+ zn lambda))
                            npn (m/* 0.25 (m/+ pn lambda))
                            nmu (m/* 0.2 (m/+ nxn nyn nzn npn npn))]
                        (if (and (m/zero? i) (or (m/inf? alpha) (m/inf? beta)))
                          ;; the first (largest) alpha, beta overflow: homogeneity as above
                          (let [s (m/max x y z p)]
                            (m// (m// (Rj (m// x s) (m// y s) (m// z s) (m// p s)) s) (m/sqrt s)))
                          (recur (m/inc i)
                                 nxn nyn nzn npn nmu
                                 (m// (m/- nmu nxn) nmu)
                                 (m// (m/- nmu nyn) nmu)
                                 (m// (m/- nmu nzn) nmu)
                                 (m// (m/- nmu npn) nmu)
                                 (m/+ sigma (m/* power4 (Rc alpha beta)))
                                 (m/* 0.25 power4)))))))))

(defn Rg
  "Rg - Carlson symmetric form of the second kind, iterative. See [[fastmath.special/elliptic-Rg]]."
  ^double [^double x ^double y ^double z]
  (cond
    (or (m/nan? x) (m/nan? y) (m/nan? z)
        (m/neg? x) (m/neg? y) (m/neg? z)) ##NaN
    (or (m/inf? x) (m/inf? y) (m/inf? z)) ##Inf
    (m/== x y z) (m/sqrt x)
    (and (m/zero? x) (m/== y z)) (m/* m/QUARTER_PI (m/sqrt y))
    (and (m/zero? x) (m/zero? y) (m/pos? z)) (m/* 0.5 (m/sqrt z))
    (m/zero? z) (Rg z x y) ;; symmetric: the branches above take the zero first (Rc(x, 0) is infinite)
    (m/== y z) (m/* 0.5 (m/+ (m/* y (Rc x y)) (m/sqrt x)))
    ;; z not the median: (x-z)(y-z) > 0 cancels in the formula below (DLMF 19.21.10); symmetric, the median as z
    (m/pos? (m/* (m/- x z) (m/- y z))) (Rg (m/min x y z) (m/max x y z) (m/max (m/min x y) (m/min (m/max x y) z)))
    ;; x*y/z or (x-z)(y-z) overflows, or all arguments tiny (the products underflow): homogeneity
    ;; Rg(lx,ly,lz) = sqrt(l) Rg(x,y,z) with l = 1/max
    (or (m/inf? (m// (m/* x y) z)) (m/inf? (m/* (m/- x z) (m/- y z))) (m/< (m/max x y z) CARLSON-RESCALE-BELOW))
    (let [hi (m/max x y z)]
      (m/* (m/sqrt hi) (Rg (m// x hi) (m// y hi) (m// z hi))))
    :else (m/* 0.5 (m/+ (m/- (m/* z (Rf x y z))
                             (m/* m/THIRD (m/- x z) (m/- y z) (Rd x y z)))
                        (m/sqrt (m// (m/* x y) z))))))

;; from Julia Special.jl

(defmacro ^:private ellip-K-from-poly
  [x m negx? & coeffs]
  `(let [t# ~x
         t# (poly/mevalpoly t# ~@coeffs)]
     (if ~negx? (m// t# (m/sqrt (m/- 1.0 ~m))) t#)))

(defn K
  "Complete elliptic integral K(m), and the incomplete one F(phi|m) with two arguments (quasi-periodic in phi); `m` is the parameter, `m<=1`. See [[fastmath.special/elliptic-K]]."
  (^double [^double phi ^double m]
   (cond
     (m/inf? phi) (if (and (m/<= m 1.0) (m/> m ##-Inf)) phi ##NaN)
     (m/neg-inf? m) (m/* 0.0 phi) ;; F(phi|-Inf) = 0
     (m/> (m/abs phi) m/HALF_PI)
     (let [phi2 (m/+ phi m/HALF_PI)]
       (m/- (m/* 2.0 (m/floor (m// phi2 m/PI)) (K m))
            (K (m/- m/HALF_PI (m/mod phi2 m/PI)) m)))
     :else (let [sinphi (m/sin phi)]
             (if (and (m/one? (m/abs sinphi)) (m/one? m))
               (m/copy-sign ##Inf sinphi)
               (let [sinphi2 (m/* sinphi sinphi)]
                 (m/* sinphi (Rf (m/- 1.0 sinphi2) (m/- 1.0 (m/* m sinphi2)) 1.0)))))))
  
  (^double [^double m]
   (let [negx? (m/neg? m)
         x (if negx? (m// m (m/dec m)) m)
         ;; 1 - x, exact for m < 0 (x rounds to 1.0 below about -2^53)
         td (if negx? (m// (m/- 1.0 m)) (m/- 1.0 x))]
     (cond
       (m/neg-inf? m) 0.0
       (or (m/> x 1.0) (m/nan? x)) ##NaN
       (m/zero? td) ##Inf
       (m/zero? x) m/HALF_PI
       (m/< x 0.1) (ellip-K-from-poly (m/- x 0.05 ) m negx?
                     1.591003453790792180 0.416000743991786912 0.245791514264103415
                     0.179481482914906162 0.144556057087555150 0.123200993312427711
                     0.108938811574293531 0.098853409871592910 0.091439629201749751
                     0.085842591595413900 0.081541118718303215)
       (m/< x 0.2) (ellip-K-from-poly (m/- x 0.15) m negx?
                     1.635256732264579992 0.471190626148732291 0.309728410831499587
                     0.252208311773135699 0.226725623219684650 0.215774446729585976
                     0.213108771877348910 0.216029124605188282 0.223255831633057896
                     0.234180501294209925 0.248557682972264071 0.266363809892617521)
       (m/< x 0.3) (ellip-K-from-poly (m/- x 0.25) m negx?
                     1.685750354812596043 0.541731848613280329 0.401524438390690257
                     0.369642473420889090 0.376060715354583645 0.405235887085125919
                     0.453294381753999079 0.520518947651184205 0.609426039204995055
                     0.724263522282908870 0.871013847709812357 1.057652872753547036)
       (m/< x 0.4) (ellip-K-from-poly (m/- x 0.35) m negx?
                     1.744350597225613243 0.634864275371935304 0.539842564164445538
                     0.571892705193787391 0.670295136265406100 0.832586590010977199
                     1.073857448247933265 1.422091460675497751 1.920387183402304829
                     2.632552548331654201 3.652109747319039160 5.115867135558865806
                     7.224080007363877411)
       (m/< x 0.5) (ellip-K-from-poly (m/- x 0.45) m negx?
                     1.813883936816982644 0.763163245700557246 0.761928605321595831
                     0.951074653668427927 1.315180671703161215 1.928560693477410941
                     2.937509342531378755 4.594894405442878062 7.330071221881720772
                     11.87151259742530180 19.45851374822937738 32.20638657246426863
                     53.73749198700554656 90.27388602940998849)
       (m/< x 0.6) (ellip-K-from-poly (m/- x 0.55) m negx?
                     1.898924910271553526 0.950521794618244435 1.151077589959015808
                     1.750239106986300540 2.952676812636875180 5.285800396121450889
                     9.832485716659979747 18.78714868327559562 36.61468615273698145
                     72.45292395127771801 145.1079577347069102 293.4786396308497026
                     598.3851815055010179 1228.420013075863451 2536.529755382764488)
       (m/< x 0.7) (ellip-K-from-poly (m/- x 0.65) m negx?
                     2.007598398424376302 1.248457231212347337 1.926234657076479729
                     3.751289640087587680 8.119944554932045802 18.66572130873555361
                     44.60392484291437063 109.5092054309498377 274.2779548232413480
                     697.5598008606326163 1795.716014500247129 4668.381716790389910
                     12235.76246813664335 32290.17809718320818 85713.07608195964685
                     228672.1890493117096 612757.2711915852774)
       (m/< x 0.8) (ellip-K-from-poly (m/- x 0.75) m negx?
                     2.156515647499643235 1.791805641849463243 3.826751287465713147
                     10.38672468363797208 31.40331405468070290 100.9237039498695416
                     337.3268282632272897 1158.707930567827917 4060.990742193632092
                     14454.00184034344795 52076.66107599404803 189493.6591462156887
                     695184.5762413896145 2567994.048255284686 9541921.966748386322
                     35634927.44218076174 133669298.4612040871 503352186.6866284541
                     1901975729.538660119 7208915015.330103756)
       (m/< x 0.85) (ellip-K-from-poly (m/- x 0.825) m negx?
                      2.318122621712510589 2.616920150291232841 7.897935075731355823
                      30.50239715446672327 131.4869365523528456 602.9847637356491617
                      2877.024617809972641 14110.51991915180325 70621.44088156540229
                      358977.2665825309926 1847238.263723971684 9600515.416049214109
                      50307677.08502366879 265444188.6527127967 1408862325.028702687
                      7515687935.373774627)
       (m/< x 0.9) (ellip-K-from-poly (m/- x 0.875) m negx?
                     2.473596173751343912 3.727624244118099310 15.60739303554930496
                     84.12850842805887747 506.9818197040613935 3252.277058145123644
                     21713.24241957434256 149037.0451890932766 1043999.331089990839
                     7427974.817042038995 53503839.67558661151 389249886.9948708474
                     2855288351.100810619 21090077038.76684053 156699833947.7902014
                     1170222242422.439893 8777948323668.937971 66101242752484.95041
                     499488053713388.7989 37859743397240299.20)
       :else (let [td1 (m/- td 0.05)
                   qd (poly/mevalpoly td
                        0.0 0.0625 0.03125 0.0205078125 0.01513671875 0.01193428039550781,
                        0.009816169738769531 0.008315593004226685 0.007199153304100037,
                        0.00633745662344154 0.00565311038371874
                        0.005097046040418718 0.004636680381850056
                        0.004249547423822886 0.003919665602267974)
                   kmd (poly/mevalpoly td1
                         1.591003453790792180 0.416000743991786912 0.245791514264103415
                         0.179481482914906162 0.144556057087555150 0.123200993312427711
                         0.108938811574293531 0.098853409871592910 0.091439629201749751
                         0.085842591595413900 0.081541118718303215)
                   t (m/* (m/- (m/log qd)) kmd m/INV_PI)]
               (if negx? (m// t (m/sqrt (m/- 1.0 m))) t))))))

;;

(defmacro ^:private ellip-E-from-poly
  [x m negx? & coeffs]
  `(let [t# ~x
         t# (poly/mevalpoly t# ~@coeffs)]
     (if ~negx? (m/* t# (m/sqrt (m/- 1.0 ~m))) t#)))

(defn E
  "Complete elliptic integral E(m), and the incomplete one E(phi|m) with two arguments (quasi-periodic in phi); `m` is the parameter, `m<=1`. See [[fastmath.special/elliptic-E]]."
  (^double [^double phi ^double m]
   (cond
     (m/inf? phi) (if (m/<= m 1.0) phi ##NaN)
     (and (m/neg-inf? m) (not (m/nan? phi))) (if (m/zero? phi) phi (m/copy-sign ##Inf phi))
     (m/> (m/abs phi) m/HALF_PI)
     (let [phi2 (m/+ phi m/HALF_PI)]
       (m/- (m/* 2.0 (m/floor (m// phi2 m/PI)) (E m))
            (E (m/- m/HALF_PI (m/mod phi2 m/PI)) m)))
     :else (let [sinphi (m/sin phi)]
             (if (and (m/one? (m/abs sinphi)) (m/one? m))
               sinphi ;; E(+-pi/2|1) = +-1 (Rf and Rd are both infinite there)
               (let [sinphi2 (m/* sinphi sinphi)
                     cosphi2 (m/- 1.0 sinphi2)
                     y (m/- 1.0 (m/* m sinphi2))
                     drf (Rf cosphi2 y 1.0)
                     drd (Rd cosphi2 y 1.0)]
                 (m/* sinphi (m/- drf (m/* m/THIRD m sinphi2 drd))))))))
  (^double [^double m]
   (let [negx? (m/neg? m)
         x (if negx? (m// m (m/dec m)) m)
         ;; 1 - x, exact for m < 0 (x rounds to 1.0 below about -2^53)
         td (if negx? (m// (m/- 1.0 m)) (m/- 1.0 x))]
     (cond
       (m/neg-inf? m) ##Inf
       (or (m/> x 1.0) (m/nan? x)) ##NaN
       (m/zero? td) 1.0
       (m/zero? x) m/HALF_PI
       (m/< x 0.1) (ellip-E-from-poly (m/- x 0.05 ) m negx?
                     +1.550973351780472328 -0.400301020103198524 -0.078498619442941939
                     -0.034318853117591992 -0.019718043317365499 -0.013059507731993309
                     -0.009442372874146547 -0.007246728512402157 -0.005807424012956090
                     -0.004809187786009338)
       (m/< x 0.2) (ellip-E-from-poly (m/- x 0.15) m negx?
                     +1.510121832092819728 -0.417116333905867549 -0.090123820404774569
                     -0.043729944019084312 -0.027965493064761785 -0.020644781177568105
                     -0.016650786739707238 -0.014261960828842520 -0.012759847429264803
                     -0.011799303775587354 -0.011197445703074968)
       (m/< x 0.3) (ellip-E-from-poly (m/- x 0.25) m negx?
                     +1.467462209339427155 -0.436576290946337775 -0.105155557666942554
                     -0.057371843593241730 -0.041391627727340220 -0.034527728505280841
                     -0.031495443512532783 -0.030527000890325277 -0.030916984019238900
                     -0.032371395314758122 -0.034789960386404158)
       (m/< x 0.4) (ellip-E-from-poly (m/- x 0.35) m negx?
                     +1.422691133490879171 -0.459513519621048674 -0.125250539822061878,
                     -0.078138545094409477 -0.064714278472050002 -0.062084339131730311,
                     -0.065197032815572477 -0.072793895362578779 -0.084959075171781003,
                     -0.102539850131045997 -0.127053585157696036 -0.160791120691274606)
       (m/< x 0.5) (ellip-E-from-poly (m/- x 0.45) m negx?
                     +1.375401971871116291 -0.487202183273184837 -0.153311701348540228
                     -0.111849444917027833 -0.108840952523135768 -0.122954223120269076
                     -0.152217163962035047 -0.200495323642697339 -0.276174333067751758
                     -0.393513114304375851 -0.575754406027879147 -0.860523235727239756
                     -1.308833205758540162)
       (m/< x 0.6) (ellip-E-from-poly (m/- x 0.55) m negx?
                     +1.325024497958230082 -0.521727647557566767 -0.194906430482126213
                     -0.171623726822011264 -0.202754652926419141 -0.278798953118534762
                     -0.420698457281005762 -0.675948400853106021 -1.136343121839229244
                     -1.976721143954398261 -3.531696773095722506 -6.446753640156048150
                     -11.97703130208884026)
       (m/< x 0.7) (ellip-E-from-poly (m/- x 0.65) m negx?
                     +1.270707479650149744 -0.566839168287866583 -0.262160793432492598
                     -0.292244173533077419 -0.440397840850423189 -0.774947641381397458
                     -1.498870837987561088 -3.089708310445186667 -6.667595903381001064
                     -14.89436036517319078 -34.18120574251449024 -80.15895841905397306
                     -191.3489480762984920 -463.5938853480342030 -1137.380822169360061)
       (m/< x 0.8) (ellip-E-from-poly (m/- x 0.75) m negx?
                     +1.211056027568459525 -0.630306413287455807 -0.387166409520669145
                     -0.592278235311934603 -1.237555584513049844 -3.032056661745247199
                     -8.181688221573590762 -23.55507217389693250 -71.04099935893064956
                     -221.8796853192349888 -712.1364793277635425 -2336.125331440396407
                     -7801.945954775964673 -26448.19586059191933 -90799.48341621365251
                     -315126.0406449163424 -1104011.344311591159)
       (m/< x 0.85) (ellip-E-from-poly (m/- x 0.825) m negx?
                      +1.161307152196282836 -0.701100284555289548 -0.580551474465437362
                      -1.243693061077786614 -3.679383613496634879 -12.81590924337895775
                      -49.25672530759985272 -202.1818735434090269 -869.8602699308701437
                      -3877.005847313289571 -17761.70710170939814 -83182.69029154232061
                      -396650.4505013548170 -1920033.413682634405)
       (m/< x 0.9) (ellip-E-from-poly (m/- x 0.875) m negx?
                     +1.124617325119752213 -0.770845056360909542 -0.844794053644911362
                     -2.490097309450394453 -10.23971741154384360 -49.74900546551479866
                     -267.0986675195705196 -1532.665883825229947 -9222.313478526091951
                     -57502.51612140314030 -368596.1167416106063 -2415611.088701091428
                     -16120097.81581656797 -109209938.5203089915 -749380758.1942496220
                     -5198725846.725541393 -36409256888.12139973)
       :else (let [td1 (m/- td 0.05)
                   kdm (poly/mevalpoly td1
                         1.591003453790792180 0.416000743991786912 0.245791514264103415
                         0.179481482914906162 0.144556057087555150 0.123200993312427711
                         0.108938811574293531 0.098853409871592910 0.091439629201749751
                         0.085842591595413900 0.081541118718303215)
                   edm (poly/mevalpoly td1
                         +1.550973351780472328 -0.400301020103198524 -0.078498619442941939
                         -0.034318853117591992 -0.019718043317365499 -0.013059507731993309
                         -0.009442372874146547 -0.007246728512402157 -0.005807424012956090
                         -0.004809187786009338)
                   hdm (m/- kdm edm)
                   ;; K(x) = sqrt(1-m) K(m) for m < 0, through the exact complement
                   km (if negx? (m/* (K m) (m/sqrt (m/- 1.0 m))) (K x))
                   t (m// (m/+ m/HALF_PI (m/* hdm km)) kdm)]
               (if negx? (m/* t (m/sqrt (m/- 1.0 m))) t))))))

;;

(defn PI
  "Complete and incomplete elliptic integral of the third kind.

  Parameters: `n` (characteristic, `n<=1`), `phi` (amplitude, any real, optional) and `m` (parameter, `m<=1`).

  Returns `NaN` for `n>1` or `m>1` (the value is complex or a principal value), `NaN` at the pole `n=1`, `Inf` for `m=1` where the integral diverges.
  The incomplete form is extended to `|phi|>π/2` by the quasi-periodicity `Π(n;φ+kπ|m)=Π(n;φ|m)+2kΠ(n|m)`."
  (^double [^double n ^double phi ^double m]
   (cond
     (m/neg-inf? n) (if (or (m/nan? phi) (m/nan? m) (m/> m 1.0)) ##NaN 0.0)
     (m/one? n) ##NaN
     (or (m/> m 1.0) (m/> n 1.0)) ##NaN
     ;; quasi-periodic with Pi(n|m) > 0 for n < 1
     (m/inf? phi) (if (and (m/< n 1.0) (m/<= m 1.0) (m/> m ##-Inf)) phi ##NaN)
     (m/> (m/abs phi) m/HALF_PI) (let [k (m/floor (m// (m/+ phi m/HALF_PI) m/PI))
                                       phi' (m/- phi (m/* k m/PI))]
                                   (m/+ (PI n phi' m) (m/* 2.0 k (PI n m))))
     :else (let [sinp (m/sin phi)]
             (if (and (m/one? (m/abs sinp)) (m/one? m))
               (m/copy-sign ##Inf sinp)
               (let [sinp2 (m/* sinp sinp)
                     cosp2 (m/- 1.0 sinp2)
                     y (m/- 1.0 (m/* m sinp2))]
                 (m/* sinp (m/+ (Rf cosp2 y 1.0)
                                (m/* sinp2 m/THIRD n
                                     (Rj cosp2 y 1.0 (m/- 1.0 (m/* n sinp2)))))))))))
  (^double [^double n ^double m]
   (cond
     (m/neg-inf? n) (if (or (m/nan? m) (m/> m 1.0)) ##NaN 0.0)
     (m/one? n) ##NaN
     (or (m/> m 1.0) (m/> n 1.0)) ##NaN
     (m/one? m) ##Inf
     (m/zero? n) (K m)
     (m/zero? m) (m// m/PI (m/* 2.0 (m/sqrt (m/- 1.0 n))))
     :else (m/+ (m/* m/THIRD n (Rj 0.0 (m/- 1.0 m) 1.0 (m/- 1.0 n)))
                (K m)))))

;; Jacobi functions

(defn- am-agm
  "am(u|k) for `0<=k<=1`, `k<1` or `k1>0`, with the complement `k1=1-k` given separately (exact when the caller knows it): the first-order series in `k` below `sqrt(tol)`, the AGM (A&S 16.4) otherwise."
  ^double [^double u ^double k ^double k1 ^double tol]
  (cond
    (or (m/nan? k) (m/nan? k1)) ##NaN ;; the AGM would not converge (e.g. m = -Inf in the callers)
    (m/< k (m/sqrt tol)) (m/- u (m/* 0.25 k (m/- u (m/* 0.5 (m/sin (m/* 2.0 u))))))
    :else
    (let [buf (double-array 20)]
      (loop [a 1.0
             b (m/sqrt k1)
             c (m/sqrt k)
             n (long 0)]
        (if (m/<= (m/abs c) tol)
          (loop [i n
                 phi (m/* a u (m/exp2 n))]
            (if (m/zero? i) phi
                (recur (m/dec i) (m/* 0.5 (m/+ phi (m/asin (m/* (Array/aget buf i)
                                                                (m/sin phi))))))))
          (let [na (m/* 0.5 (m/+ a b))
                nb (m/sqrt (m/* a b))
                nc (m/* 0.5 (m/- a b))
                nn (m/inc n)]
            (Array/aset buf nn (m// nc na))
            (recur na nb nc nn)))))))

;;  Phase u*sqrt(1-m) of am for m < 0 above which am is its linear part pi u/(2K(m)).
(def ^{:private true :const true :tag 'double} AM-LINEAR-PHASE 1.0e290)

(defn am
  "Jacobi amplitude am(u|m), the inverse of the incomplete elliptic integral F; `k` is the parameter `m` (not the modulus), `tol` the optional tolerance of the AGM iteration. For `m<0` the imaginary-modulus transformation, for `m>1` the reciprocal-parameter one (`asin(sn)`, bounded). See [[fastmath.special/jacobi-am]]."
  (^double [^double u ^double k] (am u k m/MACHINE-EPSILON))
  (^double [^double u ^double k ^double tol]
   (cond
     (or (m/nan? u) (m/nan? k)) ##NaN
     (m/zero? u) u ;; odd: keeps the sign of zero
     (m/inf? u) (cond
                  (m/< k 1.0) u
                  (m/one? k) (m/copy-sign m/HALF_PI u)
                  :else ##NaN)
     (m/one? k) (m/* 2.0 (m/atan (m/tanh (m/* 0.5 u)))) ;; gd(u)
     (m/pos-inf? k) 0.0 ;; the limit: |am| <= asin(1/sqrt(m))
     ;; tan am(u|m) = c tan am(u/c|mu), c^2 = 1/(1-m) = 1-mu, mu = -m/(1-m); the angle form keeps the quadrant
     (m/neg? k) (let [ku1 (m// (m/- 1.0 k))
                      c (m/sqrt ku1)
                      uc (m// u c)]
                  (if (m/> (m/abs uc) AM-LINEAR-PHASE)
                    ;; am = pi u/(2K) + a periodic part below pi/2, far below the ulp here (am > 1e287); the AGM
                    ;; scales the phase by 2^n and overflows next to the double range
                    (m/* m/HALF_PI (m// u (K k)))
                    (let [phi (am-agm uc (m/* (m/- k) ku1) ku1 tol)
                          s (m/sin phi)
                          co (m/cos phi)]
                      (m/+ phi (m/atan2 (m/* (m/- c 1.0) s co)
                                        (m/+ (m/* co co) (m/* c s s)))))))
     ;; the complement of 1/m is (m-1)/m, exact (1 - 1/m loses digits next to m = 1)
     ;; am = atan2(sn, cn) with cn > 0 (asin(sn) loses digits at the peaks, where sn is next to 1)
     (m/> k 1.0) (let [sk (m/sqrt k)
                       ku1 (m// (m/dec k) k)
                       phi (am-agm (m/* u sk) (m// k) ku1 tol)
                       s (m/sin phi)
                       c (m/cos phi)]
                   (m/atan2 (m// s sk) (m/sqrt (m/+ (m/* c c) (m/* ku1 s s)))))
     :else (am-agm u k (m/- 1.0 k) tol))))

(defn jsn
  "Jacobi elliptic function sn(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-sn]]."
  ^double [^double u ^double k]
  (cond
    (or (m/nan? u) (m/nan? k)) ##NaN
    (m/zero? u) u ;; odd: keeps the sign of zero
    (m/pos-inf? k) (if (m/inf? u) ##NaN 0.0) ;; the limit for a finite u: |sn| <= 1/sqrt(m)
    ;; 1 - ku sin^2 = ku1 + ku cos^2 (no cancellation for ku near 1)
    (m/neg? k) (let [ku1 (m// (m/- 1.0 k))
                     ku (m/* -1.0 k ku1)
                     sqrtku1 (m/sqrt ku1)
                     phi (am-agm (m// u sqrtku1) ku ku1 m/MACHINE-EPSILON)
                     s (m/sin phi)
                     c (m/cos phi)]
                 (m// (m/* s sqrtku1)
                      (m/sqrt (m/+ ku1 (m/* ku c c)))))
    (m/> k 1.0) (let [ku (m// k)]
                  (m/* (m/sqrt ku)
                       (m/sin (am-agm (m/* u (m/sqrt k)) ku (m// (m/dec k) k) m/MACHINE-EPSILON))))
    :else (m/sin (am u k))))

(defn jcn
  "Jacobi elliptic function cn(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-cn]]."
  ^double [^double u ^double k]
  (cond
    (or (m/nan? u) (m/nan? k)) ##NaN
    (m/zero? u) 1.0
    (m/pos-inf? k) (if (m/inf? u) ##NaN 1.0) ;; the limit for a finite u: sqrt(1-1/m) <= cn <= 1
    (m/neg? k) (let [ku1 (m// (m/- 1.0 k))
                     ku (m/* -1.0 k ku1)
                     sqrtku1 (m/sqrt ku1)
                     phi (am-agm (m// u sqrtku1) ku ku1 m/MACHINE-EPSILON)
                     c (m/cos phi)]
                 (m// c (m/sqrt (m/+ ku1 (m/* ku c c)))))
    ;; 1 - sin^2/m = cos^2 + (m-1)/m sin^2 (no cancellation next to m = 1)
    (m/> k 1.0) (let [ku (m// k)
                      ku1 (m// (m/dec k) k)
                      phi (am-agm (m/* u (m/sqrt k)) ku ku1 m/MACHINE-EPSILON)
                      c (m/cos phi)
                      s (m/sin phi)]
                  (m/sqrt (m/+ (m/* c c) (m/* ku1 s s))))
    (m/one? k) (m// (m/cosh u))
    :else (m/cos (am u k))))

(defn jdn
  "Jacobi elliptic function dn(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-dn]]."
  ^double [^double u ^double k]
  (cond
    (or (m/nan? u) (m/nan? k)) ##NaN
    (m/zero? u) 1.0
    (m/pos-inf? k) ##NaN ;; no limit: dn oscillates between -1 and 1
    (m/neg? k) (let [ku1 (m// (m/- 1.0 k))
                     ku (m/* -1.0 k ku1)
                     sqrtku1 (m/sqrt ku1)
                     c (m/cos (am-agm (m// u sqrtku1) ku ku1 m/MACHINE-EPSILON))]
                 (m// (m/sqrt (m/+ ku1 (m/* ku c c)))))
    (m/> k 1.0) (m/cos (am-agm (m/* u (m/sqrt k)) (m// k) (m// (m/dec k) k) m/MACHINE-EPSILON))
    (m/one? k) (m// (m/cosh u))
    ;; 1 - m sn^2 = cn^2 + (1-m) sn^2: no cancellation next to the minimum sqrt(1-m) of dn (eps/dn^2 otherwise)
    (m/> k 0.5) (let [phi (am u k)
                                                 c (m/cos phi)
                                                 s (m/sin phi)]
                                             (m/sqrt (m/+ (m/* c c) (m/* (m/- 1.0 k) s s))))
    :else (m/sqrt (m/- 1.0 (m/* k (m/sq (m/sin (am u k))))))))

(defn jsc
  "Jacobi elliptic function sc(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-sc]]."
  ^double [^double u ^double k] (m// (jsn u k) (jcn u k)))

(defn jsd
  "Jacobi elliptic function sd(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-sd]]."
  ^double [^double u ^double k] (m// (jsn u k) (jdn u k)))

(defn jcs
  "Jacobi elliptic function cs(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-cs]]."
  ^double [^double u ^double k]
  (if (and (m/pos-inf? k) (not (m/zero? u))) ##NaN (m// (jcn u k) (jsn u k)))) ;; sn -> 0 through both signs at m = Inf

(defn jcd
  "Jacobi elliptic function cd(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-cd]]."
  ^double [^double u ^double k] (m// (jcn u k) (jdn u k)))

(defn jds
  "Jacobi elliptic function ds(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-ds]]."
  ^double [^double u ^double k] (m// (jdn u k) (jsn u k)))

(defn jdc
  "Jacobi elliptic function dc(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-dc]]."
  ^double [^double u ^double k] (m// (jdn u k) (jcn u k)))

(defn jns
  "Jacobi elliptic function ns(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-ns]]."
  ^double [^double u ^double k]
  (if (and (m/pos-inf? k) (not (m/zero? u))) ##NaN (m// (jsn u k)))) ;; sn -> 0 through both signs at m = Inf

(defn jnc
  "Jacobi elliptic function nc(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-nc]]."
  ^double [^double u ^double k] (m// (jcn u k)))

(defn jnd
  "Jacobi elliptic function nd(u|m): the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-nd]]."
  ^double [^double u ^double k] (m// (jdn u k)))

;; inverses, https://dlmf.nist.gov/19.25#v
;;
;; Branches at negative x (DLMF 22.15): where f takes both signs on [0, 2K] the inverse is on [0, 2K] (af(-x) = 2K - af(x));
;; an odd f with positive values only on (0, 2K) is extended by oddness; an even f with positive values only has no negative x.
;; For x^2 overflowing, the Carlson homogeneity Rf(lx,ly,lz) = Rf(x,y,z)/sqrt(l) with l = 1/x^2 keeps the argument finite.

(defn- real-quarter-period
  "Real quarter period of the Jacobi functions of the parameter `k`: K(k) for `k<=1`, K(1/k)/sqrt(k) for `k>1`."
  ^double [^double k]
  (if (m/> k 1.0) (m// (Rf 0.0 (m// (m/dec k) k) 1.0) (m/sqrt k)) (K k))) ;; K(1/k) = Rf(0, (k-1)/k, 1), exact complement

(defn- reflect-half-period
  "The inverse on the second quarter period, `2K - u`, for the value at `-x` from `u` at `x`."
  ^double [^double u ^double k]
  (m/- (m/* 2.0 (real-quarter-period k)) u))

;; sq

;; dc,nc
(defn jasc
  "Inverse of the Jacobi function sc: `u` with sc(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-asc]]."
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [w (m/sq (m// x))] (m/copy-sign (Rf w (m/+ w (m/- 1.0 k)) (m/+ w 1.0)) x))
      (m/* x (Rf 1.0 (m/inc (m/* x2 (m/- 1.0 k))) (m/inc x2))))))

;; cd,nd
(defn jasd
  "Inverse of the Jacobi function sd: `u` with sd(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-asd]]."
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [w (m/sq (m// x))] (m/copy-sign (Rf w (m/+ w (m/dec k)) (m/+ w k)) x))
      (m/* x (Rf 1.0 (m/inc (m/* x2 (m/dec k))) (m/inc (m/* x2 k)))))))

;; cn,dn
(defn jasn
  "Inverse of the Jacobi function sn: `u` with sn(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). See [[fastmath.special/jacobi-asn]]."
  ^double [^double x ^double k]
  (let [x2 (m/* x x)] (m/* x (Rf 1.0 (m/inc (m/- x2)) (m/inc (m/* x2 (m/- k)))))))

;; ps, the value for x >= 0

;; dc,nc
(defn- jacs-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [w (m/sq (m// x))] (m// (Rf 1.0 (m/inc (m/* w (m/- 1.0 k))) (m/inc w)) (m/abs x)))
      (Rf x2 (m/+ x2 (m/- 1.0 k)) (m/inc x2)))))

;; cd,nd
(defn- jads-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [w (m/sq (m// x))] (m// (Rf 1.0 (m/inc (m/* w (m/dec k))) (m/inc (m/* w k))) (m/abs x)))
      (Rf x2 (m/+ x2 (m/dec k)) (m/+ x2 k)))))

;; cn,dn
(defn- jans-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [w (m/sq (m// x))] (m// (Rf 1.0 (m/- 1.0 w) (m/- 1.0 (m/* w k))) (m/abs x)))
      (Rf x2 (m/dec x2) (m/+ x2 (m/- k))))))

(defn jacs
  "Inverse of the Jacobi function cs: `u` with cs(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). For `m<1` on (0, 2K) (DLMF 22.15.23), odd for `m>=1`. See [[fastmath.special/jacobi-acs]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (jacs-nonneg x k)
    (m/< k 1.0) (reflect-half-period (jacs-nonneg (m/- x) k) k)
    :else (m/- (jacs-nonneg (m/- x) k))))

(defn jads
  "Inverse of the Jacobi function ds: `u` with ds(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). Odd for `m<=1`, on (0, 2K) for `m>1`. See [[fastmath.special/jacobi-ads]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (jads-nonneg x k)
    (m/> k 1.0) (reflect-half-period (jads-nonneg (m/- x) k) k)
    :else (m/- (jads-nonneg (m/- x) k))))

(defn jans
  "Inverse of the Jacobi function ns: `u` with ns(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). Odd. See [[fastmath.special/jacobi-ans]]."
  ^double [^double x ^double k]
  (if (m/neg? x) (m/- (jans-nonneg (m/- x) k)) (jans-nonneg x k)))

;; pq, the value for x >= 0

(defn- jacd-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [t (m/sq (m// x))
            w (m// (m/dec t) (m/- 1.0 k))]
        (m/* (m/sqrt w) (Rf 1.0 t (m/+ t (m/* w k)))))
      (let [w (m// (m/- 1.0 x2) (m/- 1.0 k))]
        (m/* (m/sqrt w) (Rf x2 1.0 (m/inc (m/* w k))))))))

(defn- jadc-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [t (m/sq (m// x))
            w (m// (m/dec t) (m/dec k))]
        (m/* (m/sqrt w) (Rf 1.0 t (m/+ t w))))
      (let [w (m// (m/- 1.0 x2) (m/dec k))]
        (m/* (m/sqrt w) (Rf x2 1.0 (m/inc w)))))))

(defn- jacn-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)
        w (m/- 1.0 x2)]
    ;; 1 - k (1 - x^2) = (1 - k) + k x^2 (no cancellation next to k = 1)
    (m/* (m/sqrt w) (Rf x2 1.0 (m/+ (m/- 1.0 k) (m/* k x2))))))

(defn- janc-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [t (m/sq (m// x))
            w (m/- 1.0 t)]
        (m/* (m/sqrt w) (Rf 1.0 t (m/+ t (m/* w (m/- 1.0 k))))))
      (let [w (m/dec x2)]
        (m/* (m/sqrt w) (Rf x2 1.0 (m/inc (m/* w (m/- 1.0 k)))))))))

(defn- jadn-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)
        w (m// (m/- 1.0 x2) k)]
    ;; 1 - (1 - x^2)/k = ((k - 1) + x^2)/k (no cancellation next to k = 1)
    (m/* (m/sqrt w) (Rf x2 1.0 (m// (m/+ (m/dec k) x2) k)))))

(defn- jand-nonneg
  ^double [^double x ^double k]
  (let [x2 (m/* x x)]
    (if (m/inf? x2)
      (let [t (m/sq (m// x))
            w (m// (m/- 1.0 t) k)]
        (m/* (m/sqrt w) (Rf 1.0 t (m/+ t (m/* w (m/dec k))))))
      (let [w (m// (m/- 1.0 x2) (m/- k))]
        (m/* (m/sqrt w) (Rf x2 1.0 (m/inc (m/* w (m/dec k)))))))))

(defn jacd
  "Inverse of the Jacobi function cd: `u` with cd(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K]. See [[fastmath.special/jacobi-acd]]."
  ^double [^double x ^double k]
  (if (m/neg? x) (reflect-half-period (jacd-nonneg (m/- x) k) k) (jacd-nonneg x k)))

(defn jadc
  "Inverse of the Jacobi function dc: `u` with dc(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K]. See [[fastmath.special/jacobi-adc]]."
  ^double [^double x ^double k]
  (if (m/neg? x) (reflect-half-period (jadc-nonneg (m/- x) k) k) (jadc-nonneg x k)))

(defn jacn
  "Inverse of the Jacobi function cn: `u` with cn(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K] for `m<1` (DLMF 22.15.6); `cn>0` for `m>=1`. See [[fastmath.special/jacobi-acn]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (jacn-nonneg x k)
    (m/< k 1.0) (reflect-half-period (jacn-nonneg (m/- x) k) k)
    :else ##NaN))

(defn janc
  "Inverse of the Jacobi function nc: `u` with nc(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K] for `m<1`; `nc>0` for `m>=1`. See [[fastmath.special/jacobi-anc]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (janc-nonneg x k)
    (m/< k 1.0) (reflect-half-period (janc-nonneg (m/- x) k) k)
    :else ##NaN))

(defn jadn
  "Inverse of the Jacobi function dn: `u` with dn(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K] for `m>1`; `dn>0` for `m<=1`. See [[fastmath.special/jacobi-adn]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (jadn-nonneg x k)
    (m/> k 1.0) (reflect-half-period (jadn-nonneg (m/- x) k) k)
    :else ##NaN))

(defn jand
  "Inverse of the Jacobi function nd: `u` with nd(u|m)=x; the second argument `k` is the parameter `m=k^2` (not the modulus). On [0, 2K] for `m>1`; `nd>0` for `m<=1`. See [[fastmath.special/jacobi-and]]."
  ^double [^double x ^double k]
  (cond
    (not (m/neg? x)) (jand-nonneg x k)
    (m/> k 1.0) (reflect-half-period (jand-nonneg (m/- x) k) k)
    :else ##NaN))
