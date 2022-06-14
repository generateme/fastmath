(ns fastmath.fields.m
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn mcarpet
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.2 1.5)
                      :y (u/sdrand 0.2 1.5)
                      :twist (u/sdrand 0.2 1.5)
                      :tilt (u/sdrand 0.2 1.5)})})
  ([^double amount {:keys [^double x ^double y ^double twist ^double tilt]}]
   (fn [^Vec2 v]
     (let [T (inc (* 0.25 (v/magsq v)))
           r (/ amount T)]
       (Vec2. (+ (* (.x v) r x) (* amount (+ (- 1.0 (* twist (.x v) (.x v))) (.y v))))
              (+ (* (.y v) r y) (* amount tilt (.x v))))))))

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
     (let [xfactor (+ xshift (* xscale (.x v)))
           yfactor (+ yshift (* yscale (.y v)))
           r (* (/ amount (+ m/EPSILON (v/magsq v)))
                (+ ushift (m/cosh yfactor))
                (m/sq (m/sin xfactor)))]
       (Vec2. (* r (m/sin xfactor))
              (* r (m/cos xfactor)))))))

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
                       (if (< it f)
                         (let [d (* d e)
                               m (+ p r)
                               n (+ q s)
                               p? (< x (/ m n))]
                           (recur (inc it)
                                  (if p? p m)
                                  (if p? q n)
                                  (if p? m r)
                                  (if p? n s)
                                  d
                                  (if p? y (+ y d))))
                         (+ y d))))]
     (fn [^Vec2 v]
       (let [isnx? (neg? (.x v))
             isny? (neg? (.y v))
             mnkx (if isnx? (- (.x v)) (.x v))
             mnkx (if (and (pos? mnkx) (< mnkx 1.0)) ^double (minkowski mnkx) mnkx)
             mnkx (if isnx? (- mnkx) mnkx)
             mnky (if isny? (- (.y v)) (.y v))
             mnky (if (and (pos? mnky) (< mnky 1.0)) ^double (minkowski mnky) mnky)
             mnky (if isny? (- mnky) mnky)]
         (Vec2. (* amount mnkx)
                (* amount mnky)))))))

(defn- minkosine
  ^double [altwave? ^double x]
  (let [lp (mod (m/abs x) 4.0)
        p (mod (m/abs x) 2.0)
        p (if (> p 1.0) (- 2.0 p) p)
        mink (if altwave? (- (m/minkowski p) p) (m/minkowski p))]
    (if (m/bool-xor (< lp 2.0) (pos? x))
      mink (- mink))))

(defn- minkocosine
  ^double [altwave? ^double x]
  (minkosine altwave? (dec x)))

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
   (let [tpf (* 0.5 frequencyx)
         tpf2 (* 0.5 frequencyy)
         nodamping? (< (m/abs damping) m/EPSILON)
         altwave? (not (pos? frequencyx))]
     (fn [^Vec2 v]
       (let [pt (* perturbation (minkosine altwave? (* tpf2 (.y v))))
             t (+ separation
                  (* amplitude (minkocosine altwave? (+ pt (* tpf (.x v))))
                     (if nodamping? 1.0 (m/exp (* damping (- (m/abs (.x v))))))))]
         (if (< (m/abs (.y v)) t)
           (v/mult v (- amount))
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
     (let [re-u (+ (- (* re-a (.x v)) (* im-a (.y v))) re-b)
           im-u (+ (+ (* re-a (.y v)) (* im-a (.x v))) im-b)
           re-v (+ (- (* re-c (.x v)) (* im-c (.y v))) re-d)
           im-v (+ (+ (* re-c (.y v)) (* im-c (.x v))) im-d)
           radv (/ amount (+ (* re-v re-v) (* im-v im-v)))]
       (Vec2. (* radv (+ (* re-u re-v) (* im-u im-v)))
              (* radv (- (* im-u re-v) (* re-u im-v))))))))

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
   (let [z (/ (* 4.0 dist) power)
         zr (/ z)
         fpower (m/floor power)]
     (fn [^Vec2 v]
       (let [r (m/pow (v/magsq v) z)
             alpha (* power (v/heading v))
             sina (m/sin alpha)
             cosa (m/cos alpha)
             x (* r cosa)
             y (* r sina)
             realu (+ (- (* re-a x) (* im-a y)) re-b)
             imagu (+ (+ (* re-a y) (* im-a x)) im-b)
             realv (+ (- (* re-c x) (* im-c y)) re-d)
             imagv (+ (+ (* re-c y) (* im-c x)) im-d)
             radv (/ (+ (* realv realv) (* imagv imagv)))
             x (* radv (+ (* realu realv) (* imagu imagv)))
             y (* radv (- (* imagu realv) (* realu imagv)))
             r (* amount (m/pow (+ (* x x) (* y y)) zr))
             n (m/floor (r/drand power))
             alpha (/ (+ (m/atan2 y x) (* n m/TWO_PI)) fpower)]
         (Vec2. (* r (m/cos alpha))
                (* r (m/sin alpha))))))))

(defn modulus
  "Modulus"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.01 2.0)
                      :y (u/sdrand 0.01 2.0)})})
  ([amount {:keys [^double x ^double y]}]
   (let [xr (+ x x)
         yr (+ y y)]
     (fn [^Vec2 v]
       (v/mult (Vec2. (cond
                        (> (.x v) x) (+ (- x) (mod (+ (.x v) x) xr))
                        (< (.x v) (- x)) (- x (mod (- x (.x v)) xr))
                        :else (.x v))
                      (cond
                        (> (.y v) y) (+ (- y) (mod (+ (.y v) y) yr))
                        (< (.y v) (- y)) (- y (mod (- y (.y v)) yr))
                        :else (.y v))) amount)))))

(defn murl2
  ([] {:type :regular
       :config (fn [] {:c (r/randval (r/drand -0.5 5.0) (r/drand -5.0 5.0))
                      :power (r/randval (u/sirand 1 5) (u/sdrand 0.01 5.0))})})
  ([^double amount {:keys [^double c ^double power]}]
   (let [p2 (* 0.5 power)
         invp (if (zero? power) 1.0e10 (/ power))
         c (if (== c -1.0) m/EPSILON (inc c))
         vp (if-not (zero? power)
              (* amount (m/pow (m/abs c) (/ 2.0 power)))
              (* amount (m/sq (m/sq c))))]
     (fn [^Vec2 v]
       (let [a (* (v/heading v) power)
             sina (m/sin a)
             cosa (m/cos a)
             r (* c (m/pow (v/magsq v) p2))
             re (inc (* r cosa))
             im (* r sina)
             r (m/pow (+ (* re re) (* im im)) invp)
             a (* (m/atan2 im re) 2.0 invp)
             re (* r (m/cos a))
             im (* r (m/sin a))
             r1 (/ vp (* r r))]
         (Vec2. (* r1 (+ (* (.x v) re) (* (.y v) im)))
                (* r1 (- (* (.y v) re) (* (.x v) im)))))))))

(defn murl
  ([] {:type :regular
       :config (fn [] {:c (r/randval (r/drand -0.5 5.0) (r/drand -5.0 5.0))
                      :power (r/randval (u/sirand 1 5) (u/sdrand 0.01 5.0))})})
  ([^double amount {:keys [^double c ^double power]}]
   (let [p2 (* 0.5 power)
         c (if (not= power 1.0) (/ c (dec power)) c)
         vp (if (== c -1.0) (* m/EPSILON amount) (* amount (inc c)))]
     (fn [^Vec2 v]
       (let [a (* (v/heading v) power)
             sina (m/sin a)
             cosa (m/cos a)
             r (* c (m/pow (v/magsq v) p2))
             re (inc (* r cosa))
             im (* r sina)
             r1 (/ vp (* r r))]
         (Vec2. (* r1 (+ (* (.x v) re) (* (.y v) im)))
                (* r1 (- (* (.y v) re) (* (.x v) im)))))))))
