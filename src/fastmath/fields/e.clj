(ns fastmath.fields.e
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.complex :as c]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn ecollide
  "Ecollide"
  ([] {:type :regular
       :config (fn [] {:num (r/randval (u/sirand 1 8) (u/sdrand 0.2 5.0))
                      :a (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double num ^double a]}]
   (let [ecn-pi (m/* num m/M_1_PI)
         pi-ecn (m// m/PI num)
         eca (m/* a m/PI)
         eca-ecn (m// eca num)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             tmp2 (m/* 2.0 (.x v))
             xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                     (m/safe-sqrt (m/- tmp tmp2)))))
             t (m/constrain (m// (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             alt (int (m/* nu ecn-pi))
             altp (m/* alt pi-ecn)
             nu (if (m/even? alt)
                  (m/+ altp (m/mod (m/+ nu eca-ecn) pi-ecn))
                  (m/+ altp (m/mod (m/- nu eca-ecn) pi-ecn)))
             nu (if-not (m/pos? (.y v)) (m/- nu) nu)]
         (Vec2. (m/* amount xmax (m/cos nu))
                (m/* amount (m/sqrt (m/dec xmax)) (m/sqrt (m/inc xmax)) (m/sin nu))))))))


(defn edisc
  "edisc"
  ([] {:type :regular})
  ([^double amount _]
   (let [w (m// amount 11.57034632)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             tmp2 (m/* 2.0 (.x v))
             xmax (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                            (m/safe-sqrt (m/- tmp tmp2))))
             a1 (m/log (m/+ xmax (m/sqrt (m/dec xmax))))
             a2 (m/- (m/acos (m// (.x v) xmax)))
             snv (m/sin a1)
             snv (if (m/pos? (.y v)) (m/- snv) snv)]
         (Vec2. (m/* w (m/cosh a2) (m/cos a1))
                (m/* w (m/sinh a2) snv)))))))

(defn ejulia
  "EJulia"
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 4) (u/sdrand 0.2 4.1))})})
  ([^double amount {:keys [^double power]}]
   (let [pp (m// m/TWO_PI power)
         sign (m/sgn power)]
     (fn [^Vec2 v]
       (let [r2 (if (m/pos? sign) (v/magsq v) (m// (v/magsq v)))
             x (if (m/pos? sign) (.x v) (m/* r2 (.x v)))
             tmp (m/inc r2)
             tmp2 (m/* 2.0 x)
             xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                     (m/safe-sqrt (m/- tmp tmp2)))))
             mu (m/acosh xmax)
             t (m/constrain (m// x xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (m/pos? (.y v)) (m/- nu) nu)
             nu (m/+ (m// nu power) (m/* pp (m/floor (r/drand power))))
             mu (m// mu power)]
         (Vec2. (m/* amount (m/cosh mu) (m/cos nu))
                (m/* amount (m/sinh mu) (m/sin nu))))))))


(defn emod
  "eMod by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:radius (r/drand 0.1 4.0)
                      :distance (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double radius ^double distance]}]
   (let [radius2 (m/* 2.0 radius)
         rdr (m/+ radius (m/* distance radius))]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             tmp2 (m/* 2.0 (.x v))
             xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                     (m/safe-sqrt (m/- tmp tmp2)))))
             mu (m/acosh xmax)          
             t (m/constrain (m// (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (m/pos? (.y v)) (m/- nu) nu)
             mu (if (and (m/< mu radius) (m/< (m/- mu) radius))
                  (if (m/pos? nu)
                    (m/- (m/mod (m/+ mu rdr) radius2) radius)
                    (m/+ (m/mod (m/- mu rdr) radius2) radius))
                  mu)
             xx (m/* amount (m/cosh mu) (m/cos nu))
             yy (m/* amount (m/sinh mu) (m/sin nu))]
         (Vec2. xx yy))))))

(defn emotion
  "eMotion by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:move (r/drand -2.0 2.0)
                      :rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double move ^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (m/inc (v/magsq v))
           tmp2 (m/* 2.0 (.x v))
           xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                   (m/safe-sqrt (m/- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (m// (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (m/pos? (.y v)) (m/- nu) nu)
           mu (if (m/neg? nu) (m/+ mu move) (m/- mu move))
           mu (if-not (m/pos? mu) (m/- mu) mu)
           nu (if-not (m/pos? mu) (m/- nu) nu)
           nu (m/+ nu rotate)
           xx (m/* amount (m/cosh mu) (m/cos nu))
           yy (m/* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn epush
  "ePush by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:push (r/drand -1.5 1.5)
                      :dist (u/sdrand 0.5 1.5)
                      :rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double push ^double dist ^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (m/inc (v/magsq v))
           tmp2 (m/* 2.0 (.x v))
           xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                   (m/safe-sqrt (m/- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (m// (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (m/pos? (.y v)) (m/- nu) nu)
           nu (m/+ nu rotate)
           mu (m/+ push (m/* mu dist))
           xx (m/* amount (m/cosh mu) (m/cos nu))
           yy (m/* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn erotate
  "eRotate by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (m/inc (v/magsq v))
           tmp2 (m/* 2.0 (.x v))
           xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                   (m/safe-sqrt (m/- tmp tmp2)))))
           t (m/constrain (m// (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (m/pos? (.y v)) (m/- nu) nu)
           nu (m/- (m/mod (m/+ nu rotate m/PI) m/TWO_PI) m/PI)
           
           xx (m/* amount xmax (m/cos nu))
           yy (m/* amount (m/sqrt (m/dec xmax)) (m/sqrt (m/inc xmax)) (m/sin nu))]
       (Vec2. xx yy)))))

(defn escale
  "eScale by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:scale (u/sdrand 0.5 1.5)
                      :angle (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double scale ^double angle]}]
   (let [angle+pi (m/+ m/PI angle)
         scale2pi (m/* m/TWO_PI scale)
         scalepi (m/* m/PI scale)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             tmp2 (m/* 2.0 (.x v))
             xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                     (m/safe-sqrt (m/- tmp tmp2)))))
             mu (m/acosh xmax)          
             t (m/constrain (m// (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (m/pos? (.y v)) (m/- nu) nu)
             mu (m/* mu scale)
             nu (m/mod (m/- (m/mod (m/* scale (m/+ nu angle+pi)) scale2pi) angle scalepi) m/TWO_PI)
             nu (if (m/> nu m/PI)
                  (m/- nu m/TWO_PI)
                  (if (m/< nu m/-PI)
                    (m/+ nu m/TWO_PI)
                    nu))
             xx (m/* amount (m/cosh mu) (m/cos nu))
             yy (m/* amount (m/sinh mu) (m/sin nu))]
         (Vec2. xx yy))))))

(defn eswirl
  "eSwirl by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:in (r/drand m/-PI m/PI)
                      :out (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double in ^double out]}]
   (fn [^Vec2 v]
     (let [tmp (m/inc (v/magsq v))
           tmp2 (m/* 2.0 (.x v))
           xmax (m/max 1.0 (m/* 0.5 (m/+ (m/safe-sqrt (m/+ tmp tmp2))
                                   (m/safe-sqrt (m/- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (m// (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (m/pos? (.y v)) (m/- nu) nu)
           nu (m/+ nu (m/* mu out) (m// in mu))
           
           xx (m/* amount (m/cosh mu) (m/cos nu))
           yy (m/* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn eclipse
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double shift]}]
   (fn [^Vec2 v]
     (if (m/< (m/abs (.y v)) amount)
       (let [c2 (m/sqrt (m/- (m/sq amount) (m/sq (.y v))))]
         (if (m/< (m/abs (.x v)) c2)
           (let [x (m/+ (.x v) (m/* shift amount))]
             (if (m/>= (m/abs x) c2)
               (Vec2. (m/* (m/- amount) (.x v)) (m/* amount (.y v)))
               (Vec2. (m/* amount x) (m/* amount (.y v)))))
           (v/mult v amount)))
       (v/mult v amount)))))

(defn elliptic2
  "Elliptic2"
  ([] {:type :random
       :config (fn [] {:a1 (r/drand -2.0 2.0)
                      :a2 (u/sdrand 0.5 2.0)
                      :a3 (r/drand -1.0 1.0)
                      :b1 (u/sdrand 0.5 3.0)
                      :b2 (u/sdrand 0.5 2.0)
                      :c (r/drand 0.2 1.2)
                      :d (r/drand 0.5 2.0)
                      :e (r/drand 0.1 0.9)
                      :f (r/drand -2.0 2.0)
                      :g (r/drand -2.0 2.0)
                      :h (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double a1 ^double a2 ^double a3
                           ^double b1 ^double b2
                           ^double c ^double d ^double e ^double f ^double g ^double h]}]
   (let [v- (m// amount (m// m/HALF_PI h))
         ps (m/* a3 (m/- m/M_PI_2))]
     (fn [^Vec2 v]
       (let [tmp (m/+ (v/magsq v) a1)
             x2 (m/* b1 (.x v))
             xmax (m/* c (m/+ (m/sqrt (m/+ tmp x2)) (m/sqrt (m/- tmp x2))))
             a (m/* a2 (m// (.x v) xmax))
             b (m/* b2 (m/safe-sqrt (m/- d (m/* a a))))
             x (m/+ ps (m/* v- (m/atan2 a b))) 
             y (r/randval e (m/* v- (m/log (m/+ xmax (m/safe-sqrt (m/- xmax f)))))
                          (m/- (m/* v- (m/log (m/+ xmax (m/safe-sqrt (m/- xmax g)))))))]
         (Vec2. x y))))))

(defn- sqrt1pm1
  ^double [^double x]
  (if (m/< -0.0625 x 0.0625)
    (let [num 0.03125
          den 0.00390625
          num (m/* num x)
          den (m/* den x)
          num (m/+ num 0.3125)
          den (m/+ den 0.15625)
          num (m/* num x)
          den (m/* den x)
          num (m/+ num 0.75)
          den (m/+ den 0.9375)
          num (m/* num x)
          den (m/* den x)
          num (m/+ num 0.5)
          den (m/+ den 1.75)
          num (m/* num x)
          den (m/* den x)]
      (m// num (m/inc den)))
    (m/dec (m/sqrt (m/inc x)))))

(defn ellipticprecision
  "Elliptic high precision version"
  ([] {:type :regular})
  ([^double amount _]
   (let [-a (m// amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [sq (v/magsq v)
             x2 (m/+ (.x v) (.x v))
             xmaxm1 (m/* 0.5 (m/+ (sqrt1pm1 (m/+ sq x2)) (sqrt1pm1 (m/- sq x2))))
             ssx (m/safe-sqrt xmaxm1)
             a (m// (.x v) (m/inc xmaxm1))
             l (m/log1p (m/+ xmaxm1 ssx))
             x (m/* -a (m/asin (m/constrain a -1.0 1.0))) 
             y (if (m/pos? (.y v)) (m/* -a l) (m/- (m/* -a l)))]
         (Vec2. x y))))))

(defn ellipticapo
  "Elliptic APO version"
  ([] {:type :regular})
  ([^double amount _]
   (let [-a (m// amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             x2 (m/+ (.x v) (.x v))
             xmax (m/* 0.5 (m/+ (m/sqrt (m/+ tmp x2)) (m/sqrt (m/- tmp x2))))
             a (m// (.x v) xmax)
             b (m/safe-sqrt (m/- 1.0 (m/* a a)))
             l (m/log (m/+ xmax (m/safe-sqrt (m/dec xmax))))
             x (m/* -a (m/atan2 a b)) 
             y (if (m/pos? (.y v)) (m/* -a l) (m/- (m/* -a l)))]
         (Vec2. x y))))))

(defn elliptic
  "Elliptic"
  ([] {:type :random})
  ([^double amount _]
   (let [-a (m// amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             x2 (m/+ (.x v) (.x v))
             xmax (m/* 0.5 (m/+ (m/sqrt (m/+ tmp x2)) (m/sqrt (m/- tmp x2))))
             a (m// (.x v) xmax)
             b (m/safe-sqrt (m/- 1.0 (m/* a a)))
             l (m/log (m/+ xmax (m/safe-sqrt (m/dec xmax))))
             x (m/* -a (m/atan2 a b)) 
             y (r/randval (m/* -a l) (m/- (m/* -a l)))]
         (Vec2. x y))))))

(defn ennepers
  "Ennepers"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sx (m/* (.x v) (.x v))
           sy (m/* (.y v) (.y v))
           x (m/+ (m/- (.x v) (m/* m/THIRD sx (.x v))) (m/* (.x v) sy))
           y (m/+ (m/- (.y v) (m/* m/THIRD sy (.y v))) (m/* (.y v) sx))]
       (Vec2. (m/* amount x) (m/* amount y))))))

(defn epispiral
  "Epispiral"
  ([] {:type :random
       :config (fn [] {:n (r/randval (r/irand -20 21) (r/drand -20.0 21.0))
                      :thickness (r/drand -1.5 1.5)
                      :holes (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double n ^double thickness ^double holes]}]
   (let [t (m/- holes)]
     (fn [^Vec2 v]
       (let [theta (v/heading v)
             d (m// (m/cos (m/* n theta)))
             t (m/* amount (if (m/> (m/abs thickness) m/EPSILON)
                           (m/+ t (m/* d (r/drand thickness)))
                           (m/+ t d)))]
         (Vec2. (m/* t (m/cos theta))
                (m/* t (m/sin theta))))))))

(defn epispiralwf
  ([] {:type :regular
       :config (fn [] {:waves (r/randval (r/irand -20 21) (r/drand -20.0 21.0))})})
  ([^double amount {:keys [^double waves]}]
   (fn [^Vec2 v]
     (let [a (v/heading v)
           r (m/* amount (m// 0.5 (m/cos (m/* waves a))))]
       (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a)))))))

(defn erf
  "Erf"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (special/erf (.x v)))
            (m/* amount (special/erf (.y v)))))))

(defn erf2
  "Erf"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (special/erf (.x v) (.y v)))
            (m/* amount (v/heading v))))))


(defn escher
  "Escher"
  ([] {:type :regular
       :config (fn [] {:beta (r/drand m/TWO_PI)})})
  ([^double amount {:keys [^double beta]}]
   (let [seb (m/sin beta)
         ceb (m/cos beta)
         vc (m/* 0.5 (m/inc ceb))
         vd (m/* 0.5 seb)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             lnr (m/* 0.5 (m/log (v/magsq v)))
             m (m/* amount (m/exp (m/- (m/* vc lnr)
                                   (m/* vd a))))
             n (m/+ (m/* vc a)
                  (m/* vd lnr))]
         (Vec2. (m/* m (m/cos n))
                (m/* m (m/sin n))))))))

(defn ex
  "Ex"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           h (v/heading v)
           n0 (m/sin (m/+ h r))
           n1 (m/cos (m/- h r))
           m0 (m/* n0 n0 n0)
           m1 (m/* n1 n1 n1)
           ar (m/* amount r)]
       (Vec2. (m/* ar (m/+ m0 m1))
              (m/* ar (m/- m0 m1)))))))

(defn exp2bs
  "Exp"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.5 2.0)
                      :y1 (u/sdrand 0.1 4.0)
                      :y2 (u/sdrand 0.1 4.0)})})
  ([^double amount {:keys [^double x1 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [e (m/* amount (m/exp (m/* x1 (.x v))))]
       (Vec2. (m/* e (m/cos (m/* y1 (.y v))))
              (m/* e (m/sin (m/* y2 (.y v)))))))))

(defn exp
  "Exp"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [e (m/* amount (m/exp (.x v)))]
       (Vec2. (m/* e (m/cos (.y v)))
              (m/* e (m/sin (.y v))))))))

(defn expmulti
  "Exp Multi"
  ([] {:type :random
       :config (fn [] {:re (u/sdrand 0.1 2.0)
                      :im (r/drand -2.0 2.0)
                      :re-add (r/drand -2.0 2.0)
                      :im-add (r/drand -2.0 2.0)
                      :sqr (r/randval 0.0 (u/sdrand 0.8 3.0))
                      :asech (r/randval 0.0 (u/sdrand 0.8 3.0))
                      :acosech (r/randval 0.0 (u/sdrand 0.8 3.0))
                      :acoth (r/randval 0.0 (u/sdrand 0.8 3.0))})})
  ([^double amount {:keys [^double re ^double im ^double re-add ^double im-add
                           ^double sqr ^double asech ^double acosech ^double acoth]}]
   (let [scale (m/* amount m/M_2_PI)
         z2 (c/complex re im)
         z3 (c/complex re-add im-add)
         cacoth (c/complex acoth 0.0)
         cacosech (c/complex acosech 0.0)
         casech (c/complex asech 0.0)
         csqr (c/complex sqr 0.0)]
     (fn [^Vec2 v]
       (let [z (-> (c/exp v)
                   (c/div z2)
                   (as-> t (if-not (m/zero? sqr) (-> (c/sq t) (c/mult csqr)) t))
                   (c/add z3)
                   (c/sqrt)
                   (as-> t (if-not (m/zero? asech) (-> (c/asech t) (c/mult casech)) t))
                   (as-> t (if-not (m/zero? acosech) (-> (c/acsch t) (c/mult cacosech)) t))
                   (as-> t (if-not (m/zero? acoth) (-> (c/acoth t) (c/mult cacoth)) t))
                   (c/scale scale))]
         (r/randval z (v/sub z)))))))

(defn exponential
  "Exponential"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [e (m/* amount (m/exp (m/dec (.x v))))
           r (m/* m/PI (.y v))]
       (Vec2. (m/* e (m/cos r))
              (m/* e (m/sin r)))))))

(defn eyefish
  "Eyefish"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m// (m/* amount 4.0) (m/inc (v/mag v)))]
       (Vec2. (m/* r (.x v)) (m/* r (.y v)))))))

