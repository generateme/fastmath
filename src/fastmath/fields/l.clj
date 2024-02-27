(ns fastmath.fields.l
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn lace
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r0 (- (* 0.5 amount (v/mag v)))
           weight (r/drand)]
       (cond
         (> weight 0.75) (let [w (m/atan2 (.y v) (dec (.x v)))]
                           (Vec2. (* r0 (m/sin w))
                                  (inc (* r0 (m/cos w)))))
         (> weight 0.5) (let [w (m/atan2 (- (.y v) m/SQRT3_2) (+ (.x v) 0.5))]
                          (Vec2. (- (* r0 (m/sin w)) 0.5)
                                 (+ (* r0 (m/cos w)) m/SQRT3_2)))
         (> weight 0.25) (let [w (m/atan2 (+ (.y v) m/SQRT3_2) (+ (.x v) 0.5))]
                           (Vec2. (- (* r0 (m/sin w)) 0.5)
                                  (- (* r0 (m/cos w)) m/SQRT3_2)))
         :else (let [w (v/heading v)]
                 (Vec2. (* r0 (m/sin w)) (* r0 (m/cos w)))))))))

(defn layeredspiral
  ([] {:type :regular
       :config (fn [] {:radius (u/sdrand 0.4 1.5)})})
  ([^double amount {:keys [^double radius]}]
   (fn [^Vec2 v]
     (let [a (* (.x v) radius)
           t (+ (v/magsq v) m/EPSILON)]
       (Vec2. (* amount a (m/sin t))
              (* amount a (m/cos t)))))))

(defn lazyjess
  ([] {:type :regular
       :config (fn [] {:n (let [u (r/drand 8.0 25.0)]
                           (r/randval (r/irand 2 (int u)) (r/drand 2.0 u)))
                      :spin (r/drand m/TWO_PI)
                      :space (r/drand -1.0 1.0)
                      :corner (r/randval (r/irand -10 10) (r/drand -10 10))})})
  ([^double amount {:keys [^double n ^double spin ^double space ^double corner]}]
   (let [vertex (/ (* m/PI (- n 2.0)) (* 2.0 n))
         sin-vertex (m/sin vertex)
         pie-slice (/ m/TWO_PI n)
         half-slice (/ pie-slice 2.0)
         corner-rotation (* (dec corner) pie-slice)
         rnomin (* amount m/SQRT2 sin-vertex)
         rdenom (- m/PI vertex)
         spin2pi (+ spin m/TWO_PI)]
     (fn [^Vec2 v]
       (let [modulus (v/mag v)
             amodulus (* amount modulus)]
         (if (== n 2.0)
           (if (< (m/abs (.x v)) amount)
             (let [theta (+ (v/heading v) spin)
                   sina (m/sin theta)
                   cosa (m/cos theta)
                   x (* amodulus cosa)
                   y (* amodulus sina)]
               (if (< (m/abs x) amount)
                 (Vec2. x y)
                 (let [theta (+ (- (m/atan2 y x) spin) corner-rotation)
                       sina (m/sin theta)
                       cosa (m/cos theta)]
                   (Vec2. (* amodulus cosa) (- (* amodulus sina))))))
             (v/mult v (* amount (inc (/ space modulus)))))
           (let [theta (+ (v/heading v) m/TWO_PI)
                 theta-diff (mod (+ theta half-slice) pie-slice)
                 r (/ rnomin (m/sin (- rdenom theta-diff)))]
             (if (< modulus r)
               (let [theta (+ (v/heading v) spin2pi)
                     sina (m/sin theta)
                     cosa (m/cos theta)
                     x (* amodulus cosa)
                     y (* amodulus sina)
                     theta-diff (mod (+ theta half-slice) pie-slice)
                     r (/ rnomin (m/sin (- rdenom theta-diff)))
                     modulus (m/hypot-sqrt x y)]
                 (if (< modulus r)
                   (Vec2. x y)
                   (let [theta (+ (- (m/atan2 y x) spin) corner-rotation m/TWO_PI)
                         sina (m/sin theta)
                         cosa (m/cos theta)
                         amodulus (* amount modulus)]
                     (Vec2. (* amodulus cosa) (- (* amodulus sina))))))
               (v/mult v (* amount (inc (/ space modulus))))))))))))

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
     (let [xx (- (.x v) x)
           yy (+ (.y v) y)
           rr (m/hypot-sqrt xx yy)]
       (if (< rr amount)
         (let [a (+ (m/atan2 yy xx) spin (* twist (- amount rr)))
               nr (* amount rr)]
           (Vec2. (+ (* nr (m/cos a)) x)
                  (- (* nr (m/sin a)) y)))
         (let [nr (* amount (inc (/ space rr)))]
           (Vec2. (+ (* nr xx) x)
                  (- (* nr yy) y))))))))

(defn lazytravis
  ([] {:type :regular
       :config (fn [] {:spin-in (r/drand -3.0 3.0)
                      :spin-out (r/drand -3.0 3.0)
                      :space (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double spin-in ^double spin-out ^double space]}]
   (let [-spin-in (* 4.0 spin-in)
         -spin-out (* 4.0 spin-out)]
     (fn [^Vec2 v]
       (let [x (m/abs (.x v))
             y (m/abs (.y v))]
         (if (or (> x amount) (> y amount))
           (let [^Vec2 sp (if (> x y)
                            (if (pos? (.x v))
                              (Vec2. x (+ x (.y v) (* x -spin-out)))
                              (Vec2. x (+ (- (* 5.0 x) (.y v)) (* x -spin-out))))
                            (if (pos? (.y v))
                              (Vec2. y (+ (- (* 3.0 y) (.x v)) (* y -spin-out)))
                              (Vec2. y (+ (* 7.0 y) (.x v) (* y -spin-out)))))
                 s (.x sp)
                 p (mod (.y sp) (* 8.0 s))]
             (-> (cond
                   (<= p (* 2.0 s)) (let [y2 (- p s)]
                                      (Vec2. (+ s space)
                                             (+ y2 (* (/ y2 s) space))))
                   (<= p (* 4.0 s)) (let [x2 (- (* 3.0 s) p)]
                                      (Vec2. (+ x2 (* (/ x2 s) space))
                                             (+ s space)))
                   (<= p (* 6.0 s)) (let [y2 (- (* 5.0 s) p)]
                                      (Vec2. (- (+ s space))
                                             (+ y2 (* (/ y2 s) space))))
                   :else            (let [x2 (- (* 7.0 s) p)]
                                      (Vec2. (+ x2 (* (/ x2 s) space))
                                             (- (+ s space)))))
                 (v/mult amount)))
           (let [^Vec2 sp (if (> x y)
                            (if (pos? (.x v))
                              (Vec2. x (+ x (.y v) (* x -spin-in)))
                              (Vec2. x (+ (- (* 5.0 x) (.y v)) (* x -spin-in))))
                            (if (pos? (.y v))
                              (Vec2. y (+ (- (* 3.0 y) (.x v)) (* y -spin-in)))
                              (Vec2. y (+ (* 7.0 y) (.x v) (* y -spin-in)))))
                 s (.x sp)
                 p (mod (.y sp) (* 8.0 s))]
             (-> (cond
                   (<= p (* 2.0 s)) (Vec2. (* amount s) (* amount (- p s)))
                   (<= p (* 4.0 s)) (Vec2. (* amount (- (* 3.0 s) p)) (* amount s))
                   (<= p (* 6.0 s)) (Vec2. (* amount -1.0 s) (* amount (- (* 5.0 s) p)))
                   :else            (Vec2. (* amount (- p (* 7.0 s))) (* amount -1.0 s)))
                 (v/mult amount)))))))))

(defn lineart
  ([] {:type :regular
       :config (fn [] {:powx (u/sdrand 0.2 2.0)
                      :powy (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double powx ^double powy]}]
   (fn [^Vec2 v]
     (Vec2. (* amount (m/sgn (.x v)) (m/pow (m/abs (.x v)) powx))
            (* amount (m/sgn (.y v)) (m/pow (m/abs (.y v)) powy))))))

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
   (let [diff (- tmax tmin)]
     (fn [_]
       (let [t (+ tmin (r/drand diff))
             y (r/drand -0.5 0.5)
             x1 (m/sin (+ d (* a t)))
             y1 (m/sin (* b t))
             z (+ (* c t) (* e y))]
         (Vec2. (* amount (+ x1 z))
                (* amount (+ y1 z))))))))

(defn logapo
  "LogApo"
  ([] {:type :regular
       :config (fn [] {:base (r/drand 0.01 10)})})
  ([^double amount {:keys [^double base]}]
   (let [adenom (* amount (/ 0.5 (m/log base)))]
     (fn [v]
       (Vec2. (* adenom (m/log (v/magsq v)))
              (* amount (v/heading v)))))))

(defn logdb
  ([] {:type :random
       :config (fn [] {:base (r/drand 0.01 10)
                      :fix-period (m/sq (r/drand 0.001 0.5))})})
  ([^double amount {:keys [^double base ^double fix-period]}]
   (let [denom (* amount (if (> base 1.0e-20) (/ 0.5 (m/log (* m/E base))) 0.5))
         fixpe (if (> base 1.0e-20) (* m/PI fix-period) m/PI)]
     (fn [^Vec2 v]
       (let [^double fix-atan-period (loop [fap 0.0
                                            i (long 0)]
                                       (if (< i 7)
                                         (let [adp (m/rint (r/drand -5.0 5.0))
                                               adp (if (> (m/abs adp) 3.0) 0.0 adp)]
                                           (recur (+ fap adp) (inc i)))
                                         (* fixpe fap)))]
         (Vec2. (* denom (v/magsq v))
                (* amount (+ (v/heading v) fix-atan-period))))))))

(defn log
  "Log"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount 0.5 (m/log (v/magsq v)))
            (* amount (v/heading v))))))

(defn loonie2
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :star (u/sdrand 0.01 1.2)
                      :circle (u/sdrand 0.01 1.2)})})
  ([^double amount {:keys [^double sides ^double star ^double circle]}]
   (let [sqrvvar (* amount amount)
         a (/ m/TWO_PI sides)
         sina (m/sin a)
         cosa (m/cos a)
         a (* m/M_PI_2 -1.0 star)
         sins (m/sin a)
         coss (m/cos a)
         a (* m/M_PI_2 circle)
         sinc (m/sin a)
         cosc (m/cos a)
         sides- (m/abs (dec sides))]
     (fn [^Vec2 v]
       (let [xrt (.x v)
             yrt (.y v)
             r2 (+ (* xrt coss) (* (m/abs yrt) sins))
             circle (v/mag v)
             ^double r2 (loop [r2 r2 xrt xrt yrt yrt i (long 0)]
                          (if (< i sides-)
                            (let [nxrt (- (* xrt cosa) (* yrt sina))
                                  nyrt (+ (* xrt sina) (* yrt cosa))]
                              (recur (max r2 (+ (* nxrt coss) (* (m/abs nyrt) sins))) nxrt nyrt (inc i)))
                            (+ (* r2 cosc) (* circle sinc))))
             r2 (if (> sides- 1) (* r2 r2) (* (m/abs r2) r2))]
         (cond
           (and (pos? r2) (< r2 sqrvvar)) (let [r (* amount (m/sqrt (m/abs (dec (/ sqrvvar r2)))))]
                                            (v/mult v r))
           (neg? r2) (let [r (/ amount (m/sqrt (dec (m/abs (/ sqrvvar r2)))))]
                       (v/mult v r))
           :else (v/mult v amount)))))))

(defn loonie3
  "Loonie"
  ([] {:type :regular})
  ([^double amount _]
   (let [sqrvvar (m/sq amount)
         r2 (* 2.0 sqrvvar)]
     (fn [^Vec2 v]
       (let [r2 (if (> (.x v) m/EPSILON) (m/sq (/ (v/magsq v) (.x v))) r2)]
         (if (< r2 sqrvvar)
           (let [r (* amount (m/sqrt (dec (/ sqrvvar r2))))]
             (v/mult v r))
           (v/mult v amount)))))))

(defn loonie
  "Loonie"
  ([] {:type :regular})
  ([^double amount _]
   (let [w2 (m/sq amount)]
     (fn [v]
       (let [r2 (v/magsq v)]
         (if (and (< r2 w2) (not (zero? r2)))
           (let [r (* amount (m/sqrt (dec (/ w2 r2))))]
             (v/mult v r))
           (v/mult v amount)))))))

(defn lozi
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.3 2.0)
                      :b (u/sdrand 0.3 2.0)
                      :c (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c]}]
   (fn [^Vec2 v]
     (Vec2. (* amount (+ (- c (* a (m/abs (.x v)))) (.y v)))
            (* amount b (.x v))))))

(m/unuse-primitive-operators)
