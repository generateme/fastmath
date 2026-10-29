(ns fastmath.fields.o
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn ortho
  ([] {:type :regular
       :config (fn [] {:in (r/drand -4.0 4.0)
                      :out (r/drand -4.0 4.0)})})
  ([^double amount {:keys [^double in ^double out]}]
   (fn [^Vec2 v]
     (let [r (v/magsq v)]
       (if (m/< r 1.0)
         (if-not (m/neg? (.x v))
           (let [xo (m// (m/inc r) (m/* 2.0 (.x v)))
                 ro (m/sqrt (m/+ (m/sq (m/- (.x v) xo)) (m/* (.y v) (.y v))))
                 theta (m/atan2 1.0 ro)
                 a (m/- (m/mod (m/+ (m/* in theta)
                              (m/atan2 (.y v) (m/- xo (.x v)))
                              theta) (m/* 2.0 theta)) theta)
                 aro (m/* amount ro)]
             (Vec2. (m/* (m/- xo (m/* (m/cos a) aro)))
                    (m/* (m/sin a) aro)))
           (let [xo (m// (m/inc r) (m/* -2.0 (.x v)))
                 ro (m/sqrt (m/+ (m/sq (m/- (m/- (.x v)) xo)) (m/* (.y v) (.y v))))
                 theta (m/atan2 1.0 ro)
                 a (m/- (m/mod (m/+ (m/* in theta)
                              (m/atan2 (.y v) (m/+ xo (.x v)))
                              theta) (m/* 2.0 theta)) theta)
                 aro (m/* amount ro)]
             (Vec2. (m/- (m/* (m/- xo (m/* (m/cos a) aro))))
                    (m/* (m/sin a) aro))))
         (let [r (m// (m/sqrt r))
               ta (v/heading v)
               ts (m/sin ta)
               tc (m/cos ta)
               x (m/* r tc)
               y (m/* r ts)]
           (if-not (m/neg? x)
             (let [xo (m// (m/inc (m/+ (m/* x x) (m/* y y))) (m/* 2.0 x))
                   ro (m/sqrt (m/+ (m/sq (m/- x xo)) (m/* y y)))
                   theta (m/atan2 1.0 ro)
                   a (m/- (m/mod (m/+ (m/* out theta)
                                (m/atan2 y (m/- xo x))
                                theta) (m/* 2.0 theta)) theta)
                   s (m/sin a)
                   c (m/cos a)
                   x (m/- xo (m/* c ro))
                   y (m/* s ro)
                   ta (m/atan2 y x)
                   ts (m/sin ta)
                   tc (m/cos ta)
                   r (m// amount (m/hypot-sqrt x y))]
               (Vec2. (m/* r tc) (m/* r ts)))
             (let [xo (m// (m/inc (m/+ (m/* x x) (m/* y y))) (m/* -2.0 x))
                   ro (m/sqrt (m/+ (m/sq (m/- (m/- x) xo)) (m/* y y)))
                   theta (m/atan2 1.0 ro)
                   a (m/- (m/mod (m/+ (m/* out theta)
                                (m/atan2 y (m/+ xo x))
                                theta) (m/* 2.0 theta)) theta)
                   s (m/sin a)
                   c (m/cos a)
                   x (m/- xo (m/* c ro))
                   y (m/* s ro)
                   ta (m/atan2 y x)
                   ts (m/sin ta)
                   tc (m/cos ta)
                   r (m// amount (m/hypot-sqrt x y))]
               (Vec2. (m/- (m/* r tc)) (m/* r ts))))))))))

(defn- octapol-hits-circle-around-origin
  ^double [^double radius ^Vec2 p]
  (if (m/zero? radius) 0.0 (v/mag p)))

(defn- octapol-hits-square-around-origin
  [^double a ^Vec2 XY]
  (and (m/<= (m/abs (.x XY)) a)
       (m/<= (m/abs (.y XY)) a)))

(defn- octapol-hits-rect
  [^Vec2 t1 ^Vec2 br ^Vec2 p]
  (and (m/>= (.x p) (.x t1))
       (m/>= (.y p) (.y t1))
       (m/<= (.x p) (.x br))
       (m/<= (.y p) (.y br))))

(defn- octapol-hits-triangle
  [^Vec2 a ^Vec2 b ^Vec2 c ^Vec2 p]
  (let [v0 (v/sub c a)
        v1 (v/sub b a)
        v2 (v/sub p a)
        d00 (v/dot v0 v0)
        d01 (v/dot v0 v1)
        d02 (v/dot v0 v2)
        d11 (v/dot v1 v1)
        d12 (v/dot v1 v2)
        denom (m/- (m/* d00 d11) (m/* d01 d01))
        ^Vec2 uv (if (m/zero? denom)
                   (Vec2. 0.0 0.0)
                   (v/div (Vec2. (m/- (m/* d11 d02) (m/* d01 d12))
                                 (m/- (m/* d00 d12) (m/* d01 d02))) denom))]
    (and (m/< (v/sum uv) 1.0) (m/pos? (.x uv)) (m/pos? (.y uv)))))

(defn octapol
  ([] {:type :regular
       :config (fn [] {:polarweight (r/randval 0.2 0.0 (r/drand 0.01 2.0))
                      :radius (r/randval 0.2 0.0 (r/drand 0.2 2.0))
                      :s (r/drand 0.01 1.0)
                      :t (r/drand 0.01 1.0)
                      :scale (u/sdrand 0.05 0.5)})})
  ([^double amount {:keys [^double polarweight ^double radius ^double s ^double t ^double scale]}]
   (let [-hs (m/* -0.5 s)
         hs (m/* 0.5 s)
         a (m/+ hs t)
         b (m/- -hs t)
         rad (m/* (m/abs radius) (m// s m/SQRT2))
         A (Vec2. -hs a)
         B (Vec2. hs a)
         C (Vec2. t hs)
         D (Vec2. t -hs)
         E (Vec2. hs b)
         F (Vec2. -hs b)
         G (Vec2. (m/- t) -hs)
         H (Vec2. (m/- t) hs)
         I (Vec2. -hs hs)
         J (Vec2. hs hs)
         K (Vec2. -hs -hs)
         L (Vec2. hs -hs)]
     (fn [^Vec2 v]
       (let [^Vec2 XY (v/mult v scale)
             r (octapol-hits-circle-around-origin rad XY)]
         (cond
           (and (m/pos? rad) (m/<= r rad)) (let [rd (m/log (m/sq (m// r rad)))
                                             phi (v/heading XY)
                                             t (m/* rd polarweight)]
                                         (v/mult (v/add XY (Vec2. (m/lerp (.x XY) phi t)
                                                                  (m/lerp (.y XY) r t))) amount))
           (and (octapol-hits-square-around-origin a XY)
                (or (octapol-hits-rect H K XY)
                    (octapol-hits-rect J D XY)
                    (octapol-hits-rect A J XY)
                    (octapol-hits-rect K E XY)
                    (octapol-hits-triangle I A H XY)
                    (octapol-hits-triangle J B C XY)
                    (octapol-hits-triangle L D E XY)
                    (octapol-hits-triangle K F G XY))) (v/mult (v/add XY XY) amount)
           :else (v/mult v amount)))))))

(defn oscilloscope
  ([] {:type :regular
       :config (fn [] {:separation (r/drand 0.1 2.5)
                      :frequency (u/sdrand 0.01 10.0)
                      :amplitude (u/sdrand 0.1 2.0)
                      :damping (r/randval 0.5 0.0 (m/sqrt (r/drand 0.01 0.99)))})})
  ([^double amount {:keys [^double separation ^double frequency ^double amplitude ^double damping]}]
   (let [tpf (m/* m/TWO_PI frequency)]
     (fn [^Vec2 v]
       (let [t (if (m/zero? damping)
                 (m/+ separation (m/* amplitude (m/cos (m/* tpf (.x v)))))
                 (m/+ separation (m/* amplitude (m/* (m/exp (m/* damping (m/- (m/abs (.x v)))))
                                               (m/cos (m/* tpf (.x v)))))))]
         (if (m/<= (m/abs (.y v)) t)
           (Vec2. (m/* amount (.x v))
                  (m/* -1.0 amount (.y v)))
           (v/mult v amount)))))))

(defn oscilloscope2
  ([] {:type :regular
       :config (fn [] {:separation (r/drand 0.1 2.5)
                      :frequency-x (u/sdrand 0.01 10.0)
                      :frequency-y (u/sdrand 0.01 10.0)
                      :amplitude (u/sdrand 0.1 2.0)
                      :perturbation (u/sdrand 0.1 2.0)
                      :damping (r/randval 0.5 0.0 (m/sqrt (r/drand 0.01 0.99)))})})
  ([^double amount {:keys [^double separation ^double frequency-x ^double frequency-y ^double perturbation ^double amplitude ^double damping]}]
   (let [tpf (m/* m/TWO_PI frequency-x)
         tpf2 (m/* m/TWO_PI frequency-y)]
     (fn [^Vec2 v]
       (let [pt (m/* perturbation (m/sin (m/* tpf2 (.y v))))
             t (if (m/zero? damping)
                 (m/+ separation (m/* amplitude (m/cos (m/+ pt (m/* tpf (.x v))))))
                 (m/+ separation (m/* amplitude (m/* (m/exp (m/* damping (m/- (m/abs (.x v)))))
                                               (m/cos (m/+ pt (m/* tpf (.x v))))))))]
         (if (m/<= (m/abs (.y v)) t)
           (v/mult v (m/- amount))
           (v/mult v amount)))))))

(defn ovoid
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.01 m/PI)
                      :y (u/sdrand 0.01 m/PI)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [r (m// amount (m/+ m/EPSILON (v/magsq v)))]
       (v/emult (v/mult v r) (Vec2. x y))))))

