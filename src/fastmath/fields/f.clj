(ns fastmath.fields.f
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn fdisc
  ([] {:type :regular
       :config (fn [] {:ashift (r/drand m/-PI m/PI)
                      :rshift (r/drand m/-PI m/PI)
                      :xshift (r/drand -2.0 2.0)
                      :yshift (r/drand -2.0 2.0)
                      :term1 (r/randval 0.0 (r/drand -2.0 2.0))
                      :term2 (r/randval 0.0 (r/drand -2.0 2.0))
                      :term3 (r/randval 0.0 (r/drand -2.0 2.0))
                      :term4 (u/sdrand 0.5 1.5)})})
  ([^double amount {:keys [^double ashift ^double rshift ^double xshift ^double yshift
                           ^double term1 ^double term2 ^double term3 ^double term4]}]
   (fn [^Vec2 v]
     (let [afactor (m// m/TWO_PI (m/+ (v/mag v) ashift))
           r (m/* 0.5 (m/+ (m/* (v/heading v) m/M_1_PI) rshift))
           xfactor (m/cos (m/+ afactor xshift))
           yfactor (m/sin (m/+ afactor yshift))
           pr (m/* amount r)
           t3pr (m/* term3 pr)
           prx (m/* pr xfactor)
           pry (m/* pr yfactor)]
       (Vec2. (m/+ (m/* term1 prx) (m/* term2 prx (.x v)) (m/* t3pr (.x v)) (m/* term4 (.x v)))
              (m/+ (m/* term1 pry) (m/* term2 pry (.y v)) (m/* t3pr (.y v)) (m/* term4 (.y v))))))))

(defn fan2
  "Fan2"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -2.0 2.0)
                      :y (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           angle (v/heading v)
           ac (m/+ angle y)
           dx (m/+ m/EPSILON (m/* m/PI x x))
           dx2 (m/* 0.5 dx)
           t (m/- ac (m/* dx (long (m// ac dx))))
           a (if (m/> t dx2)
               (m/- angle dx2)
               (m/+ angle dx2))]
       (Vec2. (m/* amount r (m/sin a))
              (m/* amount r (m/cos a)))))))

(defn fan
  "Fan"
  ([] {:type :regular
       :config (fn [] {:coeff20 (r/drand -2.0 2.0)
                      :coeff21 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double coeff20 ^double coeff21]}]
   (let [dx (m/+ m/EPSILON (m/* m/PI (m/sq coeff20)))
         dx2 (m/* 0.5 dx)]
     (fn [^Vec2 v]
       (let [angle (v/heading v)
             r (m/* amount (v/mag v))
             ac (m/+ angle coeff21)
             a (if (m/> ^double (m/mod ac dx) dx2)
                 (m/- angle dx2)
                 (m/+ angle dx2))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(def ^{:const true :private true :tag 'double} fib-fnatlog (m/log m/PHI))

(defn fibonacci2
  ([] {:type :regular
       :config (fn [] {:sc (u/sdrand 0.2 2.0)
                      :sc2 (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double sc ^double sc2]}]
   (let [affive (m// amount m/SQRT5)]
     (fn [^Vec2 v]
       (let [a (m/* (.y v) fib-fnatlog)
             snum1 (m/sin a)
             cnum1 (m/cos a)
             b (m/- (m/+ (m/* (.x v) m/PI) a))
             snum2 (m/sin b)
             cnum2 (m/cos b)
             aa (m/* (.x v) fib-fnatlog)
             eradius1 (m/* sc (m/exp (m/* sc2 aa)))
             eradius2 (m/* sc (m/exp (m/* sc2 -1.0 (m/- aa (m/* (.y v) m/PI)))))]
         (Vec2. (m/* affive (m/- (m/* eradius1 cnum1) (m/* eradius2 cnum2)))
                (m/* affive (m/- (m/* eradius1 snum1) (m/* eradius2 snum2)))))))))

(defn fisheye
  "Fisheye"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m// (m/* amount 4.0) (m/inc (v/mag v)))]
       (Vec2. (m/* r (.y v)) (m/* r (.x v)))))))

(defn flipcircle
  ([] {:type :regular})
  ([^double amount _]
   (let [samount (m/* amount amount)]
     (fn [^Vec2 v]
       (if (m/> (v/magsq v) samount)
         (v/mult v amount)
         (Vec2. (m/* amount (.x v)) (m/* amount -1.0 (.y v))))))))

(defn flipy
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (if-not (m/pos? (.x v))
       (v/mult v amount)
       (Vec2. (m/* amount (.x v)) (m/* amount -1.0 (.y v)))))))

(defn flower
  "Flower"
  ([] {:type :random
       :config (fn [] {:petals (r/randval (u/sirand 1 11)(u/sdrand 0.1 10.0))
                      :holes (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double petals ^double holes]}]
   (fn [v]
     (let [theta (v/heading v)
           d (m// 1.0 (m/+ m/EPSILON (v/mag v)))
           r (m/* amount (m/- (r/drand) holes) (m/cos (m/* petals theta)) d)]
       (v/mult v r)))))

(defn flux
  "Flux"
  ([] {:type :regular
       :config (fn [] {:spread (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double spread]}]
   (let [aspread2 (m/* amount (m/+ 2.0 spread))]
     (fn [^Vec2 v]
       (let [xpw (m/+ (.x v) amount)
             xmw (m/- (.x v) amount)
             y2 (m/* (.y v) (.y v))
             avgr (m/* aspread2 (m/sqrt (m// (m/sqrt (m/+ y2 (m/* xpw xpw)))
                                         (m/sqrt (m/+ y2 (m/* xmw xmw))))))
             avga (m/* 0.5 (m/- (m/atan2 (.y v) xmw) (m/atan2 (.y v) xpw)))]
         (Vec2. (m/* avgr (m/cos avga))
                (m/* avgr (m/sin avga))))))))

(defn foci
  "Foci"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [expx (m/* 0.5 (m/exp (.x v)))
           expnx (m// 0.25 expx)
           sy (m/sin (.y v))
           cy (m/cos (.y v))
           tmp (m/- (m/+ expx expnx) cy)
           tmp (m// amount (if (m/zero? tmp) m/EPSILON tmp))]
       (Vec2. (m/* tmp (m/- expx expnx))
              (m/* tmp sy))))))

(defn fourth
  ([] {:type :regular
       :config (fn [] {:spin (r/drand m/-TWO_PI m/TWO_PI)
                      :space (u/sdrand 0.1 2.0)
                      :twist (r/drand m/-TWO_PI m/TWO_PI)
                      :x (r/randval 0.0 (r/drand -0.5 0.5))
                      :y (r/randval 0.0 (r/drand -0.5 0.5))})})
  ([^double amount {:keys [^double spin ^double space ^double twist
                           ^double x ^double y]}]
   (let [sqrvvar (m/* amount amount)
         ^Vec2 x-y (Vec2. (m/- x) y)
         ^Vec2 xy- (Vec2. x (m/- y))]
     (fn [^Vec2 v]
       (cond
         (and (m/pos? (.x v))
              (m/pos? (.y v))) (let [a (v/heading v)
                                   r (m// (v/mag v))
                                   s (m/sin a)
                                   c (m/cos a)]
                               (Vec2. (m/* amount r c) (m/* amount r s)))
         (and (m/pos? (.x v))
              (m/neg? (.y v))) (let [r2 (v/magsq v)]
                               (if (m/< r2 sqrvvar)
                                 (let [r (m/* amount (m/sqrt (m/dec (m// sqrvvar r2))))]
                                   (v/mult v r))
                                 (v/mult v amount)))
         (and (m/neg? (.x v))
              (m/pos? (.y v))) (let [xy (v/add v x-y)
                                   r (v/mag xy)]
                               (if (m/< r amount)
                                 (let [a (m/+ (v/heading xy) spin (m/* twist (m/- amount r)))
                                       r (m/* amount r)]
                                   (v/add (Vec2. (m/* r (m/cos a))
                                                 (m/* r (m/sin a))) xy-))
                                 (let [r (m/* amount (m/inc (m// space r)))]
                                   (v/add (v/mult xy r) xy-))))
         :else (v/mult v amount))))))

(defn funnel
  ([] {:type :regular
       :config (fn [] {:effect (r/drand -1.6 1.6)})})
  ([^double amount {:keys [^double effect]}]
   (let [effect (m/* effect m/PI)]
     (fn [v]
       (v/emult (v/mult (v/tanh v) amount)
                (v/shift (v/reciprocal (v/cos v)) effect))))))

;;


(defn foucaut
  "Foucaut"
  ([] {:type :regular})
  ([^double amount _]
   (let [as (m/* amount m/SQRTPI)]
     (fn [^Vec2 v]
       (let [k (m/* 0.5 (.y v))
             cosk (m/cos k)
             xx (->> cosk
                     (m/* cosk)
                     (m/* (m/cos (.y v)))
                     (m/* (m// (.x v) m/SQRTPI))
                     (m/* 2.0)
                     (m/* amount))
             yy (m/* as (m/tan k))]
         (Vec2. xx yy))))))

