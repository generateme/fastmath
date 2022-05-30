(ns fastmath.fields.f
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *warn-on-reflection* true)

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
     (let [afactor (/ m/TWO_PI (+ (v/mag v) ashift))
           r (* 0.5 (+ (* (v/heading v) m/M_1_PI) rshift))
           xfactor (m/cos (+ afactor xshift))
           yfactor (m/sin (+ afactor yshift))
           pr (* amount r)
           t3pr (* term3 pr)
           prx (* pr xfactor)
           pry (* pr yfactor)]
       (Vec2. (+ (* term1 prx) (* term2 prx (.x v)) (* t3pr (.x v)) (* term4 (.x v)))
              (+ (* term1 pry) (* term2 pry (.y v)) (* t3pr (.y v)) (* term4 (.y v))))))))

(defn fan2
  "Fan2"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -2.0 2.0)
                      :y (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           angle (v/heading v)
           ac (+ angle y)
           dx (+ m/EPSILON (* m/PI x x))
           dx2 (* 0.5 dx)
           t (- ac (* dx (long (/ ac dx))))
           a (if (> t dx2)
               (- angle dx2)
               (+ angle dx2))]
       (Vec2. (* amount r (m/sin a))
              (* amount r (m/cos a)))))))

(defn fan
  "Fan"
  ([] {:type :regular
       :config (fn [] {:coeff20 (r/drand -2.0 2.0)
                      :coeff21 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double coeff20 ^double coeff21]}]
   (let [dx (+ m/EPSILON (* m/PI (m/sq coeff20)))
         dx2 (* 0.5 dx)]
     (fn [^Vec2 v]
       (let [angle (v/heading v)
             r (* amount (v/mag v))
             ac (+ angle coeff21)
             a (if (> ^double (mod ac dx) dx2)
                 (- angle dx2)
                 (+ angle dx2))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(def ^:const ^:private ^double fib-fnatlog (m/log m/PHI))

(defn fibonacci2
  ([] {:type :regular
       :config (fn [] {:sc (u/sdrand 0.2 2.0)
                      :sc2 (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double sc ^double sc2]}]
   (let [affive (/ amount m/SQRT5)]
     (fn [^Vec2 v]
       (let [a (* (.y v) fib-fnatlog)
             snum1 (m/sin a)
             cnum1 (m/cos a)
             b (- (+ (* (.x v) m/PI) a))
             snum2 (m/sin b)
             cnum2 (m/cos b)
             aa (* (.x v) fib-fnatlog)
             eradius1 (* sc (m/exp (* sc2 aa)))
             eradius2 (* sc (m/exp (* sc2 -1.0 (- aa (* (.y v) m/PI)))))]
         (Vec2. (* affive (- (* eradius1 cnum1) (* eradius2 cnum2)))
                (* affive (- (* eradius1 snum1) (* eradius2 snum2)))))))))

(defn fisheye
  "Fisheye"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (/ (* amount 4.0) (inc (v/mag v)))]
       (Vec2. (* r (.y v)) (* r (.x v)))))))

(defn flower
  "Flower"
  ([] {:type :random
       :config (fn [] {:petals (r/randval (u/sirand 1 11)(u/sdrand 0.1 10.0))
                      :holes (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double petals ^double holes]}]
   (fn [v]
     (let [theta (v/heading v)
           d (/ 1.0 (+ m/EPSILON (v/mag v)))
           r (* amount (- (r/drand) holes) (m/cos (* petals theta)) d)]
       (v/mult v r)))))

(defn foci
  "Foci"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [expx (* 0.5 (m/exp (.x v)))
           expnx (/ 0.25 expx)
           sy (m/sin (.y v))
           cy (m/cos (.y v))
           tmp (- (+ expx expnx) cy)
           tmp (/ amount (if (zero? tmp) m/EPSILON tmp))]
       (Vec2. (* tmp (- expx expnx))
              (* tmp sy))))))

