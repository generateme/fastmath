(ns fastmath.fields.e
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.complex :as c]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn ecollide
  "Ecollide"
  ([] {:type :regular
       :config (fn [] {:num (r/randval (u/sirand 1 8) (u/sdrand 0.2 5.0))
                      :a (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double num ^double a]}]
   (let [ecn-pi (* num m/M_1_PI)
         pi-ecn (/ m/PI num)
         eca (* a m/PI)
         eca-ecn (/ eca num)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             tmp2 (* 2.0 (.x v))
             xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                     (m/safe-sqrt (- tmp tmp2)))))
             t (m/constrain (/ (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             alt (int (* nu ecn-pi))
             altp (* alt pi-ecn)
             nu (if (even? alt)
                  (+ altp (mod (+ nu eca-ecn) pi-ecn))
                  (+ altp (mod (- nu eca-ecn) pi-ecn)))
             nu (if-not (pos? (.y v)) (- nu) nu)]
         (Vec2. (* amount xmax (m/cos nu))
                (* amount (m/sqrt (dec xmax)) (m/sqrt (inc xmax)) (m/sin nu))))))))


(defn edisc
  "edisc"
  ([] {:type :regular})
  ([^double amount _]
   (let [w (/ amount 11.57034632)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             tmp2 (* 2.0 (.x v))
             xmax (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                            (m/safe-sqrt (- tmp tmp2))))
             a1 (m/log (+ xmax (m/sqrt (dec xmax))))
             a2 (- (m/acos (/ (.x v) xmax)))
             snv (m/sin a1)
             snv (if (pos? (.y v)) (- snv) snv)]
         (Vec2. (* w (m/cosh a2) (m/cos a1))
                (* w (m/sinh a2) snv)))))))

(defn ejulia
  "EJulia"
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 4) (u/sdrand 0.2 4.1))})})
  ([^double amount {:keys [^double power]}]
   (let [pp (/ m/TWO_PI power)
         sign (m/sgn power)]
     (fn [^Vec2 v]
       (let [r2 (if (pos? sign) (v/magsq v) (/ (v/magsq v)))
             x (if (pos? sign) (.x v) (* r2 (.x v)))
             tmp (inc r2)
             tmp2 (* 2.0 x)
             xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                     (m/safe-sqrt (- tmp tmp2)))))
             mu (m/acosh xmax)
             t (m/constrain (/ x xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (pos? (.y v)) (- nu) nu)
             nu (+ (/ nu power) (* pp (m/floor (r/drand power))))
             mu (/ mu power)]
         (Vec2. (* amount (m/cosh mu) (m/cos nu))
                (* amount (m/sinh mu) (m/sin nu))))))))


(defn emod
  "eMod by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:radius (r/drand 0.1 4.0)
                      :distance (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double radius ^double distance]}]
   (let [radius2 (* 2.0 radius)
         rdr (+ radius (* distance radius))]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             tmp2 (* 2.0 (.x v))
             xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                     (m/safe-sqrt (- tmp tmp2)))))
             mu (m/acosh xmax)          
             t (m/constrain (/ (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (pos? (.y v)) (- nu) nu)
             mu (if (and (< mu radius) (< (- mu) radius))
                  (if (pos? nu)
                    (- (mod (+ mu rdr) radius2) radius)
                    (+ (mod (- mu rdr) radius2) radius))
                  mu)
             xx (* amount (m/cosh mu) (m/cos nu))
             yy (* amount (m/sinh mu) (m/sin nu))]
         (Vec2. xx yy))))))

(defn emotion
  "eMotion by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:move (r/drand -2.0 2.0)
                      :rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double move ^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (inc (v/magsq v))
           tmp2 (* 2.0 (.x v))
           xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                   (m/safe-sqrt (- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (/ (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (pos? (.y v)) (- nu) nu)
           mu (if (neg? nu) (+ mu move) (- mu move))
           mu (if-not (pos? mu) (- mu) mu)
           nu (if-not (pos? mu) (- nu) nu)
           nu (+ nu rotate)
           xx (* amount (m/cosh mu) (m/cos nu))
           yy (* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn epush
  "ePush by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:push (r/drand -1.5 1.5)
                      :dist (u/sdrand 0.5 1.5)
                      :rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double push ^double dist ^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (inc (v/magsq v))
           tmp2 (* 2.0 (.x v))
           xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                   (m/safe-sqrt (- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (/ (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (pos? (.y v)) (- nu) nu)
           nu (+ nu rotate)
           mu (+ push (* mu dist))
           xx (* amount (m/cosh mu) (m/cos nu))
           yy (* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn erotate
  "eRotate by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:rotate (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double rotate]}]
   (fn [^Vec2 v]
     (let [tmp (inc (v/magsq v))
           tmp2 (* 2.0 (.x v))
           xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                   (m/safe-sqrt (- tmp tmp2)))))
           t (m/constrain (/ (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (pos? (.y v)) (- nu) nu)
           nu (- (mod (+ nu rotate m/PI) m/TWO_PI) m/PI)
           
           xx (* amount xmax (m/cos nu))
           yy (* amount (m/sqrt (dec xmax)) (m/sqrt (inc xmax)) (m/sin nu))]
       (Vec2. xx yy)))))

(defn escale
  "eScale by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:scale (u/sdrand 0.5 1.5)
                      :angle (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double scale ^double angle]}]
   (let [angle+pi (+ m/PI angle)
         scale2pi (* m/TWO_PI scale)
         scalepi (* m/PI scale)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             tmp2 (* 2.0 (.x v))
             xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                     (m/safe-sqrt (- tmp tmp2)))))
             mu (m/acosh xmax)          
             t (m/constrain (/ (.x v) xmax) -1.0 1.0)
             nu (m/acos t)
             nu (if-not (pos? (.y v)) (- nu) nu)
             mu (* mu scale)
             nu (mod (- (mod (* scale (+ nu angle+pi)) scale2pi) angle scalepi) m/TWO_PI)
             nu (if (> nu m/PI)
                  (- nu m/TWO_PI)
                  (if (< nu m/-PI)
                    (+ nu m/TWO_PI)
                    nu))
             xx (* amount (m/cosh mu) (m/cos nu))
             yy (* amount (m/sinh mu) (m/sin nu))]
         (Vec2. xx yy))))))

(defn eswirl
  "eSwirl by Michael Faber, http://michaelfaber.deviantart.com/art/eSeries-306044892"
  ([] {:type :regular
       :config (fn [] {:in (r/drand m/-PI m/PI)
                      :out (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double in ^double out]}]
   (fn [^Vec2 v]
     (let [tmp (inc (v/magsq v))
           tmp2 (* 2.0 (.x v))
           xmax (max 1.0 (* 0.5 (+ (m/safe-sqrt (+ tmp tmp2))
                                   (m/safe-sqrt (- tmp tmp2)))))
           mu (m/acosh xmax)          
           t (m/constrain (/ (.x v) xmax) -1.0 1.0)
           nu (m/acos t)
           nu (if-not (pos? (.y v)) (- nu) nu)
           nu (+ nu (* mu out) (/ in mu))
           
           xx (* amount (m/cosh mu) (m/cos nu))
           yy (* amount (m/sinh mu) (m/sin nu))]
       (Vec2. xx yy)))))

(defn eclipse
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double shift]}]
   (fn [^Vec2 v]
     (if (< (m/abs (.y v)) amount)
       (let [c2 (m/sqrt (- (m/sq amount) (m/sq (.y v))))]
         (if (< (m/abs (.x v)) c2)
           (let [x (+ (.x v) (* shift amount))]
             (if (>= (m/abs x) c2)
               (Vec2. (* (- amount) (.x v)) (* amount (.y v)))
               (Vec2. (* amount x) (* amount (.y v)))))
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
   (let [v- (/ amount (/ m/HALF_PI h))
         ps (* a3 (- m/M_PI_2))]
     (fn [^Vec2 v]
       (let [tmp (+ (v/magsq v) a1)
             x2 (* b1 (.x v))
             xmax (* c (+ (m/sqrt (+ tmp x2)) (m/sqrt (- tmp x2))))
             a (* a2 (/ (.x v) xmax))
             b (* b2 (m/safe-sqrt (- d (* a a))))
             x (+ ps (* v- (m/atan2 a b))) 
             y (r/randval e (* v- (m/log (+ xmax (m/safe-sqrt (- xmax f)))))
                          (- (* v- (m/log (+ xmax (m/safe-sqrt (- xmax g)))))))]
         (Vec2. x y))))))

(defn- sqrt1pm1
  ^double [^double x]
  (if (< -0.0625 x 0.0625)
    (let [num 0.03125
          den 0.00390625
          num (* num x)
          den (* den x)
          num (+ num 0.3125)
          den (+ den 0.15625)
          num (* num x)
          den (* den x)
          num (+ num 0.75)
          den (+ den 0.9375)
          num (* num x)
          den (* den x)
          num (+ num 0.5)
          den (+ den 1.75)
          num (* num x)
          den (* den x)]
      (/ num (inc den)))
    (dec (m/sqrt (inc x)))))

(defn ellipticprecision
  "Elliptic high precision version"
  ([] {:type :regular})
  ([^double amount _]
   (let [-a (/ amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [sq (v/magsq v)
             x2 (+ (.x v) (.x v))
             xmaxm1 (* 0.5 (+ (sqrt1pm1 (+ sq x2)) (sqrt1pm1 (- sq x2))))
             ssx (m/safe-sqrt xmaxm1)
             a (/ (.x v) (inc xmaxm1))
             l (m/log1p (+ xmaxm1 ssx))
             x (* -a (m/asin (m/constrain a -1.0 1.0))) 
             y (if (pos? (.y v)) (* -a l) (- (* -a l)))]
         (Vec2. x y))))))

(defn ellipticapo
  "Elliptic APO version"
  ([] {:type :regular})
  ([^double amount _]
   (let [-a (/ amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             x2 (+ (.x v) (.x v))
             xmax (* 0.5 (+ (m/sqrt (+ tmp x2)) (m/sqrt (- tmp x2))))
             a (/ (.x v) xmax)
             b (m/safe-sqrt (- 1.0 (* a a)))
             l (m/log (+ xmax (m/safe-sqrt (dec xmax))))
             x (* -a (m/atan2 a b)) 
             y (if (pos? (.y v)) (* -a l) (- (* -a l)))]
         (Vec2. x y))))))

(defn elliptic
  "Elliptic"
  ([] {:type :random})
  ([^double amount _]
   (let [-a (/ amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             x2 (+ (.x v) (.x v))
             xmax (* 0.5 (+ (m/sqrt (+ tmp x2)) (m/sqrt (- tmp x2))))
             a (/ (.x v) xmax)
             b (m/safe-sqrt (- 1.0 (* a a)))
             l (m/log (+ xmax (m/safe-sqrt (dec xmax))))
             x (* -a (m/atan2 a b)) 
             y (r/randval (* -a l) (- (* -a l)))]
         (Vec2. x y))))))

(defn ennepers
  "Ennepers"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sx (* (.x v) (.x v))
           sy (* (.y v) (.y v))
           x (+ (- (.x v) (* m/THIRD sx (.x v))) (* (.x v) sy))
           y (+ (- (.y v) (* m/THIRD sy (.y v))) (* (.y v) sx))]
       (Vec2. (* amount x) (* amount y))))))

(defn epispiral
  "Epispiral"
  ([] {:type :random
       :config (fn [] {:n (r/randval (r/irand -20 21) (r/drand -20.0 21.0))
                      :thickness (r/drand -1.5 1.5)
                      :holes (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double n ^double thickness ^double holes]}]
   (let [t (- holes)]
     (fn [^Vec2 v]
       (let [theta (v/heading v)
             d (/ (m/cos (* n theta)))
             t (* amount (if (> (m/abs thickness) m/EPSILON)
                           (+ t (* d (r/drand thickness)))
                           (+ t d)))]
         (Vec2. (* t (m/cos theta))
                (* t (m/sin theta))))))))

(defn epispiralwf
  ([] {:type :regular
       :config (fn [] {:waves (r/randval (r/irand -20 21) (r/drand -20.0 21.0))})})
  ([^double amount {:keys [^double waves]}]
   (fn [^Vec2 v]
     (let [a (v/heading v)
           r (* amount (/ 0.5 (m/cos (* waves a))))]
       (Vec2. (* r (m/cos a)) (* r (m/sin a)))))))

(defn erf
  "Erf"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (special/erf (.x v)))
            (* amount (special/erf (.y v)))))))

(defn erf2
  "Erf"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (special/erf (.x v) (.y v)))
            (* amount (v/heading v))))))


(defn escher
  "Escher"
  ([] {:type :regular
       :config (fn [] {:beta (r/drand m/TWO_PI)})})
  ([^double amount {:keys [^double beta]}]
   (let [seb (m/sin beta)
         ceb (m/cos beta)
         vc (* 0.5 (inc ceb))
         vd (* 0.5 seb)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             lnr (* 0.5 (m/log (v/magsq v)))
             m (* amount (m/exp (- (* vc lnr)
                                   (* vd a))))
             n (+ (* vc a)
                  (* vd lnr))]
         (Vec2. (* m (m/cos n))
                (* m (m/sin n))))))))

(defn ex
  "Ex"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           h (v/heading v)
           n0 (m/sin (+ h r))
           n1 (m/cos (- h r))
           m0 (* n0 n0 n0)
           m1 (* n1 n1 n1)
           ar (* amount r)]
       (Vec2. (* ar (+ m0 m1))
              (* ar (- m0 m1)))))))

(defn exp2bs
  "Exp"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.5 2.0)
                      :y1 (u/sdrand 0.1 4.0)
                      :y2 (u/sdrand 0.1 4.0)})})
  ([^double amount {:keys [^double x1 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [e (* amount (m/exp (* x1 (.x v))))]
       (Vec2. (* e (m/cos (* y1 (.y v))))
              (* e (m/sin (* y2 (.y v)))))))))

(defn exp
  "Exp"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [e (* amount (m/exp (.x v)))]
       (Vec2. (* e (m/cos (.y v)))
              (* e (m/sin (.y v))))))))

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
   (let [scale (* amount m/M_2_PI)
         z2 (c/complex re im)
         z3 (c/complex re-add im-add)
         cacoth (c/complex acoth 0.0)
         cacosech (c/complex acosech 0.0)
         casech (c/complex asech 0.0)
         csqr (c/complex sqr 0.0)]
     (fn [^Vec2 v]
       (let [z (-> (c/exp v)
                   (c/div z2)
                   (as-> t (if-not (zero? sqr) (-> (c/sq t) (c/mult csqr)) t))
                   (c/add z3)
                   (c/sqrt)
                   (as-> t (if-not (zero? asech) (-> (c/asech t) (c/mult casech)) t))
                   (as-> t (if-not (zero? acosech) (-> (c/acsch t) (c/mult cacosech)) t))
                   (as-> t (if-not (zero? acoth) (-> (c/acoth t) (c/mult cacoth)) t))
                   (c/scale scale))]
         (r/randval z (v/sub z)))))))

(defn exponential
  "Exponential"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [e (* amount (m/exp (dec (.x v))))
           r (* m/PI (.y v))]
       (Vec2. (* e (m/cos r))
              (* e (m/sin r)))))))

(defn eyefish
  "Eyefish"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (/ (* amount 4.0) (inc (v/mag v)))]
       (Vec2. (* r (.x v)) (* r (.y v)))))))

(m/unuse-primitive-operators)
