(ns fastmath.fields.l
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn lace
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r0 (m/- (m/* 0.5 amount (v/mag v)))
           weight (r/drand)]
       (cond
         (m/> weight 0.75) (let [w (m/atan2 (.y v) (m/dec (.x v)))]
                           (Vec2. (m/* r0 (m/sin w))
                                  (m/inc (m/* r0 (m/cos w)))))
         (m/> weight 0.5) (let [w (m/atan2 (m/- (.y v) m/SQRT3_2) (m/+ (.x v) 0.5))]
                          (Vec2. (m/- (m/* r0 (m/sin w)) 0.5)
                                 (m/+ (m/* r0 (m/cos w)) m/SQRT3_2)))
         (m/> weight 0.25) (let [w (m/atan2 (m/+ (.y v) m/SQRT3_2) (m/+ (.x v) 0.5))]
                           (Vec2. (m/- (m/* r0 (m/sin w)) 0.5)
                                  (m/- (m/* r0 (m/cos w)) m/SQRT3_2)))
         :else (let [w (v/heading v)]
                 (Vec2. (m/* r0 (m/sin w)) (m/* r0 (m/cos w)))))))))

(defn layeredspiral
  ([] {:type :regular
       :config (fn [] {:radius (u/sdrand 0.4 1.5)})})
  ([^double amount {:keys [^double radius]}]
   (fn [^Vec2 v]
     (let [a (m/* (.x v) radius)
           t (m/+ (v/magsq v) m/EPSILON)]
       (Vec2. (m/* amount a (m/sin t))
              (m/* amount a (m/cos t)))))))

(defn lazyjess
  ([] {:type :regular
       :config (fn [] {:n (let [u (r/drand 8.0 25.0)]
                           (r/randval (r/irand 2 (int u)) (r/drand 2.0 u)))
                      :spin (r/drand m/TWO_PI)
                      :space (r/drand -1.0 1.0)
                      :corner (r/randval (r/irand -10 10) (r/drand -10 10))})})
  ([^double amount {:keys [^double n ^double spin ^double space ^double corner]}]
   (let [vertex (m// (m/* m/PI (m/- n 2.0)) (m/* 2.0 n))
         sin-vertex (m/sin vertex)
         pie-slice (m// m/TWO_PI n)
         half-slice (m// pie-slice 2.0)
         corner-rotation (m/* (m/dec corner) pie-slice)
         rnomin (m/* amount m/SQRT2 sin-vertex)
         rdenom (m/- m/PI vertex)
         spin2pi (m/+ spin m/TWO_PI)]
     (fn [^Vec2 v]
       (let [modulus (v/mag v)
             amodulus (m/* amount modulus)]
         (if (m/== n 2.0)
           (if (m/< (m/abs (.x v)) amount)
             (let [theta (m/+ (v/heading v) spin)
                   sina (m/sin theta)
                   cosa (m/cos theta)
                   x (m/* amodulus cosa)
                   y (m/* amodulus sina)]
               (if (m/< (m/abs x) amount)
                 (Vec2. x y)
                 (let [theta (m/+ (m/- (m/atan2 y x) spin) corner-rotation)
                       sina (m/sin theta)
                       cosa (m/cos theta)]
                   (Vec2. (m/* amodulus cosa) (m/- (m/* amodulus sina))))))
             (v/mult v (m/* amount (m/inc (m// space modulus)))))
           (let [theta (m/+ (v/heading v) m/TWO_PI)
                 theta-diff (m/mod (m/+ theta half-slice) pie-slice)
                 r (m// rnomin (m/sin (m/- rdenom theta-diff)))]
             (if (m/< modulus r)
               (let [theta (m/+ (v/heading v) spin2pi)
                     sina (m/sin theta)
                     cosa (m/cos theta)
                     x (m/* amodulus cosa)
                     y (m/* amodulus sina)
                     theta-diff (m/mod (m/+ theta half-slice) pie-slice)
                     r (m// rnomin (m/sin (m/- rdenom theta-diff)))
                     modulus (m/hypot-sqrt x y)]
                 (if (m/< modulus r)
                   (Vec2. x y)
                   (let [theta (m/+ (m/- (m/atan2 y x) spin) corner-rotation m/TWO_PI)
                         sina (m/sin theta)
                         cosa (m/cos theta)
                         amodulus (m/* amount modulus)]
                     (Vec2. (m/* amodulus cosa) (m/- (m/* amodulus sina))))))
               (v/mult v (m/* amount (m/inc (m// space modulus))))))))))))

(defn lazysusan
  "Lazysusan"
  ([] {:type :regular
       :config (fn [] {:twist (r/drand -6.0 6.0)
                      :spin (r/drand -4.0 4.0)
                      :space (r/drand -2.0 2.0)
                      :x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double twist ^double spin ^double space ^double x ^double y]}]
   (fn [^Vec2 v]
     (let [xx (m/- (.x v) x)
           yy (m/+ (.y v) y)
           rr (m/hypot-sqrt xx yy)]
       (if (m/< rr amount)
         (let [a (m/+ (m/atan2 yy xx) spin (m/* twist (m/- amount rr)))
               nr (m/* amount rr)]
           (Vec2. (m/+ (m/* nr (m/cos a)) x)
                  (m/- (m/* nr (m/sin a)) y)))
         (let [nr (m/* amount (m/inc (m// space rr)))]
           (Vec2. (m/+ (m/* nr xx) x)
                  (m/- (m/* nr yy) y))))))))

(defn lazytravis
  ([] {:type :regular
       :config (fn [] {:spin-in (r/drand -3.0 3.0)
                      :spin-out (r/drand -3.0 3.0)
                      :space (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double spin-in ^double spin-out ^double space]}]
   (let [-spin-in (m/* 4.0 spin-in)
         -spin-out (m/* 4.0 spin-out)]
     (fn [^Vec2 v]
       (let [x (m/abs (.x v))
             y (m/abs (.y v))]
         (if (or (m/> x amount) (m/> y amount))
           (let [^Vec2 sp (if (m/> x y)
                            (if (m/pos? (.x v))
                              (Vec2. x (m/+ x (.y v) (m/* x -spin-out)))
                              (Vec2. x (m/+ (m/- (m/* 5.0 x) (.y v)) (m/* x -spin-out))))
                            (if (m/pos? (.y v))
                              (Vec2. y (m/+ (m/- (m/* 3.0 y) (.x v)) (m/* y -spin-out)))
                              (Vec2. y (m/+ (m/* 7.0 y) (.x v) (m/* y -spin-out)))))
                 s (.x sp)
                 p (m/mod (.y sp) (m/* 8.0 s))]
             (-> (cond
                   (m/<= p (m/* 2.0 s)) (let [y2 (m/- p s)]
                                      (Vec2. (m/+ s space)
                                             (m/+ y2 (m/* (m// y2 s) space))))
                   (m/<= p (m/* 4.0 s)) (let [x2 (m/- (m/* 3.0 s) p)]
                                      (Vec2. (m/+ x2 (m/* (m// x2 s) space))
                                             (m/+ s space)))
                   (m/<= p (m/* 6.0 s)) (let [y2 (m/- (m/* 5.0 s) p)]
                                      (Vec2. (m/- (m/+ s space))
                                             (m/+ y2 (m/* (m// y2 s) space))))
                   :else            (let [x2 (m/- (m/* 7.0 s) p)]
                                      (Vec2. (m/+ x2 (m/* (m// x2 s) space))
                                             (m/- (m/+ s space)))))
                 (v/mult amount)))
           (let [^Vec2 sp (if (m/> x y)
                            (if (m/pos? (.x v))
                              (Vec2. x (m/+ x (.y v) (m/* x -spin-in)))
                              (Vec2. x (m/+ (m/- (m/* 5.0 x) (.y v)) (m/* x -spin-in))))
                            (if (m/pos? (.y v))
                              (Vec2. y (m/+ (m/- (m/* 3.0 y) (.x v)) (m/* y -spin-in)))
                              (Vec2. y (m/+ (m/* 7.0 y) (.x v) (m/* y -spin-in)))))
                 s (.x sp)
                 p (m/mod (.y sp) (m/* 8.0 s))]
             (-> (cond
                   (m/<= p (m/* 2.0 s)) (Vec2. (m/* amount s) (m/* amount (m/- p s)))
                   (m/<= p (m/* 4.0 s)) (Vec2. (m/* amount (m/- (m/* 3.0 s) p)) (m/* amount s))
                   (m/<= p (m/* 6.0 s)) (Vec2. (m/* amount -1.0 s) (m/* amount (m/- (m/* 5.0 s) p)))
                   :else            (Vec2. (m/* amount (m/- p (m/* 7.0 s))) (m/* amount -1.0 s)))
                 (v/mult amount)))))))))

(defn lineart
  ([] {:type :regular
       :config (fn [] {:powx (u/sdrand 0.2 2.0)
                      :powy (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double powx ^double powy]}]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (m/sgn (.x v)) (m/pow (m/abs (.x v)) powx))
            (m/* amount (m/sgn (.y v)) (m/pow (m/abs (.y v)) powy))))))

(defn lissajous
  ([] {:type :pattern
       :config (fn [] {:tmin (r/randval m/-PI (r/drand m/-TWO_PI m/-HALF_PI))
                      :tmax (r/randval m/PI (r/drand m/HALF_PI m/TWO_PI))
                      :a (r/randval (r/irand -6 7) (r/drand -6.0 7.0))
                      :b (r/randval (r/irand -6 7) (r/drand -6.0 7.0))
                      :c (r/drand -0.5 0.5)
                      :d (r/drand m/TWO_PI)
                      :e (m/sq (r/drand 0.0 1.0))})})
  ([^double amount {:keys [^double tmin ^double tmax
                           ^double a ^double b ^double c ^double d ^double e]}]
   (let [diff (m/- tmax tmin)]
     (fn [_]
       (let [t (m/+ tmin (r/drand diff))
             y (r/drand -0.5 0.5)
             x1 (m/sin (m/+ d (m/* a t)))
             y1 (m/sin (m/* b t))
             z (m/+ (m/* c t) (m/* e y))]
         (Vec2. (m/* amount (m/+ x1 z))
                (m/* amount (m/+ y1 z))))))))

(defn logapo
  "LogApo"
  ([] {:type :regular
       :config (fn [] {:base (r/drand 0.01 10)})})
  ([^double amount {:keys [^double base]}]
   (let [adenom (m/* amount (m// 0.5 (m/log base)))]
     (fn [v]
       (Vec2. (m/* adenom (m/log (v/magsq v)))
              (m/* amount (v/heading v)))))))

(defn logdb
  ([] {:type :random
       :config (fn [] {:base (r/drand 0.01 10)
                      :fix-period (m/sq (r/drand 0.001 0.5))})})
  ([^double amount {:keys [^double base ^double fix-period]}]
   (let [denom (m/* amount (if (m/> base 1.0e-20) (m// 0.5 (m/log (m/* m/E base))) 0.5))
         fixpe (if (m/> base 1.0e-20) (m/* m/PI fix-period) m/PI)]
     (fn [^Vec2 v]
       (let [^double fix-atan-period (loop [fap 0.0
                                            i (long 0)]
                                       (if (m/< i 7)
                                         (let [adp (m/rint (r/drand -5.0 5.0))
                                               adp (if (m/> (m/abs adp) 3.0) 0.0 adp)]
                                           (recur (m/+ fap adp) (m/inc i)))
                                         (m/* fixpe fap)))]
         (Vec2. (m/* denom (v/magsq v))
                (m/* amount (m/+ (v/heading v) fix-atan-period))))))))

(defn log
  "Log"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount 0.5 (m/log (v/magsq v)))
            (m/* amount (v/heading v))))))

(defn loonie2
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :star (u/sdrand 0.01 1.2)
                      :circle (u/sdrand 0.01 1.2)})})
  ([^double amount {:keys [^double sides ^double star ^double circle]}]
   (let [sqrvvar (m/* amount amount)
         a (m// m/TWO_PI sides)
         sina (m/sin a)
         cosa (m/cos a)
         a (m/* m/M_PI_2 -1.0 star)
         sins (m/sin a)
         coss (m/cos a)
         a (m/* m/M_PI_2 circle)
         sinc (m/sin a)
         cosc (m/cos a)
         sides- (m/abs (m/dec sides))]
     (fn [^Vec2 v]
       (let [xrt (.x v)
             yrt (.y v)
             r2 (m/+ (m/* xrt coss) (m/* (m/abs yrt) sins))
             circle (v/mag v)
             ^double r2 (loop [r2 r2 xrt xrt yrt yrt i (long 0)]
                          (if (m/< i sides-)
                            (let [nxrt (m/- (m/* xrt cosa) (m/* yrt sina))
                                  nyrt (m/+ (m/* xrt sina) (m/* yrt cosa))]
                              (recur (m/max r2 (m/+ (m/* nxrt coss) (m/* (m/abs nyrt) sins))) nxrt nyrt (m/inc i)))
                            (m/+ (m/* r2 cosc) (m/* circle sinc))))
             r2 (if (m/> sides- 1) (m/* r2 r2) (m/* (m/abs r2) r2))]
         (cond
           (and (m/pos? r2) (m/< r2 sqrvvar)) (let [r (m/* amount (m/sqrt (m/abs (m/dec (m// sqrvvar r2)))))]
                                            (v/mult v r))
           (m/neg? r2) (let [r (m// amount (m/sqrt (m/dec (m/abs (m// sqrvvar r2)))))]
                       (v/mult v r))
           :else (v/mult v amount)))))))

(defn loonie3
  "Loonie"
  ([] {:type :regular})
  ([^double amount _]
   (let [sqrvvar (m/sq amount)
         r2 (m/* 2.0 sqrvvar)]
     (fn [^Vec2 v]
       (let [r2 (if (m/> (.x v) m/EPSILON) (m/sq (m// (v/magsq v) (.x v))) r2)]
         (if (m/< r2 sqrvvar)
           (let [r (m/* amount (m/sqrt (m/dec (m// sqrvvar r2))))]
             (v/mult v r))
           (v/mult v amount)))))))

(defn loonie
  "Loonie"
  ([] {:type :regular})
  ([^double amount _]
   (let [w2 (m/sq amount)]
     (fn [v]
       (let [r2 (v/magsq v)]
         (if (and (m/< r2 w2) (not (m/zero? r2)))
           (let [r (m/* amount (m/sqrt (m/dec (m// w2 r2))))]
             (v/mult v r))
           (v/mult v amount)))))))

(defn lozi
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.3 2.0)
                      :b (u/sdrand 0.3 2.0)
                      :c (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c]}]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (m/+ (m/- c (m/* a (m/abs (.x v)))) (.y v)))
            (m/* amount b (.x v))))))

