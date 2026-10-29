(ns fastmath.fields.m
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn mcarpet
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.2 1.5)
                      :y (u/sdrand 0.2 1.5)
                      :twist (u/sdrand 0.2 1.5)
                      :tilt (u/sdrand 0.2 1.5)})})
  ([^double amount {:keys [^double x ^double y ^double twist ^double tilt]}]
   (fn [^Vec2 v]
     (let [T (m/inc (m/* 0.25 (v/magsq v)))
           r (m// amount T)]
       (Vec2. (m/+ (m/* (.x v) r x) (m/* amount (m/+ (m/- 1.0 (m/* twist (.x v) (.x v))) (.y v))))
              (m/+ (m/* (.y v) r y) (m/* amount tilt (.x v))))))))

(defn mask
  ([] {:type :regular
       :config (fn [] {:xshift (r/randval 0.3 0.0 (r/drand -2.0 2.0))
                      :yshift (r/randval 0.3 0.0 (r/drand -2.0 2.0))
                      :ushift (r/drand -3.0 3.0)
                      :xscale (r/drand -5.0 5.0)
                      :yscale (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double xshift ^double yshift ^double ushift
                           ^double xscale ^double yscale]}]
   (fn [^Vec2 v]
     (let [xfactor (m/+ xshift (m/* xscale (.x v)))
           yfactor (m/+ yshift (m/* yscale (.y v)))
           r (m/* (m// amount (m/+ m/EPSILON (v/magsq v)))
                (m/+ ushift (m/cosh yfactor))
                (m/sq (m/sin xfactor)))]
       (Vec2. (m/* r (m/sin xfactor))
              (m/* r (m/cos xfactor)))))))

(defn minkqm
  ([] {:type :regular
       :config (fn [] {:a (r/randval (r/drand 1.0) (r/drand -1.01 3.0))
                      :b (r/randval (r/drand 1.0) (r/drand -1.01 3.0))
                      :c (r/randval (r/drand 1.0) (r/drand -1.01 3.0))
                      :dd (r/drand -5.0 5.0)
                      :e (r/drand -1.0 1.0)
                      :f (r/irand 3 20)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double dd ^double e ^long f]}]
   (let [minkowski (fn ^double [^double x]
                     (loop [it (long 0) p 0.0 q a r b s c d dd y 0.0]
                       (if (m/< it f)
                         (let [d (m/* d e)
                               m (m/+ p r)
                               n (m/+ q s)
                               p? (m/< x (m// m n))]
                           (recur (m/inc it)
                                  (if p? p m)
                                  (if p? q n)
                                  (if p? m r)
                                  (if p? n s)
                                  d
                                  (if p? y (m/+ y d))))
                         (m/+ y d))))]
     (fn [^Vec2 v]
       (let [isnx? (m/neg? (.x v))
             isny? (m/neg? (.y v))
             mnkx (if isnx? (m/- (.x v)) (.x v))
             mnkx (if (and (m/pos? mnkx) (m/< mnkx 1.0)) ^double (minkowski mnkx) mnkx)
             mnkx (if isnx? (m/- mnkx) mnkx)
             mnky (if isny? (m/- (.y v)) (.y v))
             mnky (if (and (m/pos? mnky) (m/< mnky 1.0)) ^double (minkowski mnky) mnky)
             mnky (if isny? (m/- mnky) mnky)]
         (Vec2. (m/* amount mnkx)
                (m/* amount mnky)))))))

(defn- minkosine
  ^double [altwave? ^double x]
  (let [lp (m/mod (m/abs x) 4.0)
        p (m/mod (m/abs x) 2.0)
        p (if (m/> p 1.0) (m/- 2.0 p) p)
        mink (if altwave? (m/- (special/minkowski p) p) (special/minkowski p))]
    (if (m/bool-xor (m/< lp 2.0) (m/pos? x))
      mink (m/- mink))))

(defn- minkocosine
  ^double [altwave? ^double x]
  (minkosine altwave? (m/dec x)))

(defn minkowskope
  ([] {:type :regular
       :config (fn [] {:separation (r/drand 0.1 2.0)
                      :frequencyx (r/drand -6.0 6.0)
                      :frequencyy (r/drand -6.0 6.0)
                      :amplitude (u/sdrand 0.3 3.0)
                      :perturbation (r/drand -4.0 4.0)
                      :damping (r/randval 0.0 (r/drand))})})
  ([^double amount {:keys [^double separation ^double frequencyx ^double frequencyy
                           ^double amplitude ^double perturbation ^double damping]}]
   (let [tpf (m/* 0.5 frequencyx)
         tpf2 (m/* 0.5 frequencyy)
         nodamping? (m/< (m/abs damping) m/EPSILON)
         altwave? (not (m/pos? frequencyx))]
     (fn [^Vec2 v]
       (let [pt (m/* perturbation (minkosine altwave? (m/* tpf2 (.y v))))
             t (m/+ separation
                  (m/* amplitude (minkocosine altwave? (m/+ pt (m/* tpf (.x v))))
                     (if nodamping? 1.0 (m/exp (m/* damping (m/- (m/abs (.x v))))))))]
         (if (m/< (m/abs (.y v)) t)
           (v/mult v (m/- amount))
           (v/mult v amount)))))))

(defn mobius
  ([] {:type :regular
       :config (fn [] {:re-a (r/drand -1.5 1.5)
                      :re-b (r/drand -1.5 1.5)
                      :re-c (r/drand -1.5 1.5)
                      :re-d (r/drand -1.5 1.5)
                      :im-a (r/drand -1.5 1.5)
                      :im-b (r/drand -1.5 1.5)
                      :im-c (r/drand -1.5 1.5)
                      :im-d (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double re-a ^double re-b ^double re-c ^double re-d
                           ^double im-a ^double im-b ^double im-c ^double im-d]}]
   (fn [^Vec2 v]
     (let [re-u (m/+ (m/- (m/* re-a (.x v)) (m/* im-a (.y v))) re-b)
           im-u (m/+ (m/+ (m/* re-a (.y v)) (m/* im-a (.x v))) im-b)
           re-v (m/+ (m/- (m/* re-c (.x v)) (m/* im-c (.y v))) re-d)
           im-v (m/+ (m/+ (m/* re-c (.y v)) (m/* im-c (.x v))) im-d)
           radv (m// amount (m/+ (m/* re-v re-v) (m/* im-v im-v)))]
       (Vec2. (m/* radv (m/+ (m/* re-u re-v) (m/* im-u im-v)))
              (m/* radv (m/- (m/* im-u re-v) (m/* re-u im-v))))))))

(defn mobiusn
  ([] {:type :random
       :config (fn [] {:re-a (r/drand -1.5 1.5)
                      :re-b (r/drand -1.5 1.5)
                      :re-c (r/drand -1.5 1.5)
                      :re-d (r/drand -1.5 1.5)
                      :im-a (r/drand -1.5 1.5)
                      :im-b (r/drand -1.5 1.5)
                      :im-c (r/drand -1.5 1.5)
                      :im-d (r/drand -1.5 1.5)
                      :power (r/randval (u/sirand 1 10) (u/sdrand 1.0 10.0))
                      :dist (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double re-a ^double re-b ^double re-c ^double re-d
                           ^double im-a ^double im-b ^double im-c ^double im-d
                           ^double power ^double dist]}]
   (let [z (m// (m/* 4.0 dist) power)
         zr (m// z)
         fpower (m/floor power)]
     (fn [^Vec2 v]
       (let [r (m/pow (v/magsq v) z)
             alpha (m/* power (v/heading v))
             sina (m/sin alpha)
             cosa (m/cos alpha)
             x (m/* r cosa)
             y (m/* r sina)
             realu (m/+ (m/- (m/* re-a x) (m/* im-a y)) re-b)
             imagu (m/+ (m/+ (m/* re-a y) (m/* im-a x)) im-b)
             realv (m/+ (m/- (m/* re-c x) (m/* im-c y)) re-d)
             imagv (m/+ (m/+ (m/* re-c y) (m/* im-c x)) im-d)
             radv (m// (m/+ (m/* realv realv) (m/* imagv imagv)))
             x (m/* radv (m/+ (m/* realu realv) (m/* imagu imagv)))
             y (m/* radv (m/- (m/* imagu realv) (m/* realu imagv)))
             r (m/* amount (m/pow (m/+ (m/* x x) (m/* y y)) zr))
             n (m/floor (r/drand power))
             alpha (m// (m/+ (m/atan2 y x) (m/* n m/TWO_PI)) fpower)]
         (Vec2. (m/* r (m/cos alpha))
                (m/* r (m/sin alpha))))))))

(defn modulus
  "Modulus"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.01 2.0)
                      :y (u/sdrand 0.01 2.0)})})
  ([amount {:keys [^double x ^double y]}]
   (let [xr (m/+ x x)
         yr (m/+ y y)]
     (fn [^Vec2 v]
       (v/mult (Vec2. (cond
                        (m/> (.x v) x) (m/+ (m/- x) (m/mod (m/+ (.x v) x) xr))
                        (m/< (.x v) (m/- x)) (m/- x (m/mod (m/- x (.x v)) xr))
                        :else (.x v))
                      (cond
                        (m/> (.y v) y) (m/+ (m/- y) (m/mod (m/+ (.y v) y) yr))
                        (m/< (.y v) (m/- y)) (m/- y (m/mod (m/- y (.y v)) yr))
                        :else (.y v))) amount)))))

(defn murl2
  ([] {:type :regular
       :config (fn [] {:c (r/randval (r/drand -0.5 5.0) (r/drand -5.0 5.0))
                      :power (r/randval (u/sirand 1 5) (u/sdrand 0.01 5.0))})})
  ([^double amount {:keys [^double c ^double power]}]
   (let [p2 (m/* 0.5 power)
         invp (if (m/zero? power) 1.0e10 (m// power))
         c (if (m/== c -1.0) m/EPSILON (m/inc c))
         vp (if-not (m/zero? power)
              (m/* amount (m/pow (m/abs c) (m// 2.0 power)))
              (m/* amount (m/sq (m/sq c))))]
     (fn [^Vec2 v]
       (let [a (m/* (v/heading v) power)
             sina (m/sin a)
             cosa (m/cos a)
             r (m/* c (m/pow (v/magsq v) p2))
             re (m/inc (m/* r cosa))
             im (m/* r sina)
             r (m/pow (m/+ (m/* re re) (m/* im im)) invp)
             a (m/* (m/atan2 im re) 2.0 invp)
             re (m/* r (m/cos a))
             im (m/* r (m/sin a))
             r1 (m// vp (m/* r r))]
         (Vec2. (m/* r1 (m/+ (m/* (.x v) re) (m/* (.y v) im)))
                (m/* r1 (m/- (m/* (.y v) re) (m/* (.x v) im)))))))))

(defn murl
  ([] {:type :regular
       :config (fn [] {:c (r/randval (r/drand -0.5 5.0) (r/drand -5.0 5.0))
                      :power (r/randval (u/sirand 1 5) (u/sdrand 0.01 5.0))})})
  ([^double amount {:keys [^double c ^double power]}]
   (let [p2 (m/* 0.5 power)
         c (if (not= power 1.0) (m// c (m/dec power)) c)
         vp (if (m/== c -1.0) (m/* m/EPSILON amount) (m/* amount (m/inc c)))]
     (fn [^Vec2 v]
       (let [a (m/* (v/heading v) power)
             sina (m/sin a)
             cosa (m/cos a)
             r (m/* c (m/pow (v/magsq v) p2))
             re (m/inc (m/* r cosa))
             im (m/* r sina)
             r1 (m// vp (m/* r r))]
         (Vec2. (m/* r1 (m/+ (m/* (.x v) re) (m/* (.y v) im)))
                (m/* r1 (m/- (m/* (.y v) re) (m/* (.x v) im)))))))))

;; https://github.com/d3/d3-geo-projection/tree/master/src

(defn miller
  "Miller"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (Vec2. (.x v)
                    (->> (m/constrain (.y v) -1.9634 1.9634)
                         (m/* 0.4)
                         (m/+ m/QUARTER_PI)
                         (m/tan)
                         (m/log)
                         (m/* 1.25))) amount))))

(defn millerrev
  "Millerrev"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (Vec2. (.x v)
                    (-> (.y v)
                        (m/* 0.8)
                        (m/exp)
                        (m/atan)
                        (m/* 2.5)
                        (m/- (m/* 0.625 m/PI)))) amount))))

