(ns fastmath.fields.i
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn idisc
  ([] {:type :regular})
  ([^double amount _]
   (let [v- (m/* amount m/M_1_PI)]
     (fn [^Vec2 v]
       (let [a (m// m/PI (m/inc (v/mag v)))
             r (m/* v- (v/heading v))]
         (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a))))))))

(defn intersection
  ([] {:type :random
       :config (fn [] {:xwidth (u/sdrand 0.2 4.0)
                      :xtilesize (u/sdrand 0.2 2.0)
                      :xmod1 (u/sdrand 0.1 2.0)
                      :xmod2 (u/sdrand 0.1 2.0)
                      :xheight (u/sdrand 0.2 4.0)
                      :ywidth (u/sdrand 0.2 4.0)
                      :ytilesize (u/sdrand 0.2 2.0)
                      :ymod1 (u/sdrand 0.1 2.0)
                      :ymod2 (u/sdrand 0.1 2.0)
                      :yheight (u/sdrand 0.2 4.0)})})
  ([^double amount {:keys [^double xwidth ^double xtilesize ^double xmod1 ^double xmod2 ^double xheight
                           ^double ywidth ^double ytilesize ^double ymod1 ^double ymod2 ^double yheight]}]
   (let [xr1 (m/* xmod2 xmod1)
         yr1 (m/* ymod2 ymod1)]
     (fn [^Vec2 v]
       (-> (r/randval
            (let [x (r/randval xwidth (m/- xwidth))]
              (Vec2. (m/* xtilesize (m/+ (.x v) (m/round (m/* x (m/log (r/drand))))))
                     (cond
                       (m/> (.y v) xmod1) (m/* xheight (m/- (m/mod (m/+ (.y v) xmod1) xr1) xmod1))
                       (m/< (.y v) (m/- xmod1)) (m/* xheight (m/- xmod1 (m/mod (m/- xmod1 (.y v)) xr1)))
                       :else (m/* xheight (.y v)))))
            (let [y (r/randval yheight (m/- yheight))]
              (Vec2. (cond
                       (m/> (.x v) ymod1) (m/* ywidth (m/- (m/mod (m/+ (.x v) ymod1) yr1) ymod1))
                       (m/< (.x v) (m/- ymod1)) (m/* ywidth (m/- ymod1 (m/mod (m/- ymod1 (.x v)) yr1)))
                       :else (m/* ywidth (.x v)))
                     (m/* ytilesize (m/+ (.y v) (m/round (m/* y (m/log (r/drand)))))))))
           (v/mult amount))))))

(defn invsquircular
  ([] {:type :regular})
  ([^double amount _]
   (let [amount2 (m/* amount amount)]
     (fn [^Vec2 v]
       (let [r (v/magsq v)
             r2 (m/sqrt (m// (m/* r (m/- (m/* amount2 r) (m/* 4.0 (.x v) (.x v) (.y v) (.y v)))) amount))
             r (m// (m/sqrt (m/- r r2)) m/M_SQRT2)]
         (Vec2. (m// r (.x v)) (m// r (.y v))))))))

(defn invtree
  "InvTree"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (cond
       (r/brand 0.333) (v/mult v (m/* 0.5 amount))
       (r/brand 0.666) (v/mult (Vec2. (m// (m/inc (.x v)))
                                      (m// (.y v) (m/inc (.y v)))) amount)
       :else (v/mult (Vec2. (m// (.x v) (m/inc (.x v)))
                            (m// (m/inc (.y v)))) amount)))))

(defn inverted-julia
  ([] {:type :random
       :config (fn [] {:power (u/sdrand 0.1 3.0)
                      :y2-mult (r/drand 0.1 2.0)
                      :a2x-mult (u/sdrand 0.5 2.0)
                      :a2y-mult (u/sdrand 0.5 2.0)
                      :a2y-add (r/drand -1.0 1.0)
                      :cos-mult (r/drand -4.0 4.0)
                      :y-mult (u/sdrand 0.5 2.0)
                      :center (u/sdrand 0.5 4.0)
                      :x2y2-add (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double power ^double y2-mult ^double a2x-mult ^double a2y-mult ^double a2y-add
                           ^double cos-mult ^double y-mult ^double center ^double x2y2-add]}]
   (fn [^Vec2 v]
     (let [xs (m/* (.x v) (.x v))
           ys (m/* (.y v) (.y v))
           z (m/+ x2y2-add (m/pow (m/+ xs (m/* ys y2-mult)) power))
           q (m/+ (m/* 0.5 (m/atan2 (m/* (.x v) a2x-mult) (m/+ (m/* (.y v) a2y-mult) a2y-add)))
                (m/* m/PI (r/lrand 2)))
           r (m/* amount (m/cos (m/* z cos-mult)) (m// (m/* z center)))]
       (Vec2. (m/* r (m/sin q)) (m/* r (m/cos q) y-mult))))))

(defn invpolar
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [ny (m/* amount (m/inc (.y v)))
           nx (m/* (.x v) m/PI)]
       (Vec2. (m/* ny (m/sin nx))
              (m/* ny (m/cos ny)))))))

