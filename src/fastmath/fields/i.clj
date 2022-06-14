(ns fastmath.fields.i
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn idisc
  ([] {:type :regular})
  ([^double amount _]
   (let [v- (* amount m/M_1_PI)]
     (fn [^Vec2 v]
       (let [a (/ m/PI (inc (v/mag v)))
             r (* v- (v/heading v))]
         (Vec2. (* r (m/cos a)) (* r (m/sin a))))))))

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
   (let [xr1 (* xmod2 xmod1)
         yr1 (* ymod2 ymod1)]
     (fn [^Vec2 v]
       (-> (r/randval
            (let [x (r/randval xwidth (- xwidth))]
              (Vec2. (* xtilesize (+ (.x v) (m/round (* x (m/log (r/drand))))))
                     (cond
                       (> (.y v) xmod1) (* xheight (- (mod (+ (.y v) xmod1) xr1) xmod1))
                       (< (.y v) (- xmod1)) (* xheight (- xmod1 (mod (- xmod1 (.y v)) xr1)))
                       :else (* xheight (.y v)))))
            (let [y (r/randval yheight (- yheight))]
              (Vec2. (cond
                       (> (.x v) ymod1) (* ywidth (- (mod (+ (.x v) ymod1) yr1) ymod1))
                       (< (.x v) (- ymod1)) (* ywidth (- ymod1 (mod (- ymod1 (.x v)) yr1)))
                       :else (* ywidth (.x v)))
                     (* ytilesize (+ (.y v) (m/round (* y (m/log (r/drand)))))))))
           (v/mult amount))))))

(defn invsquircular
  ([] {:type :regular})
  ([^double amount _]
   (let [amount2 (* amount amount)]
     (fn [^Vec2 v]
       (let [r (v/magsq v)
             r2 (m/sqrt (/ (* r (- (* amount2 r) (* 4.0 (.x v) (.x v) (.y v) (.y v)))) amount))
             r (/ (m/sqrt (- r r2)) m/M_SQRT2)]
         (Vec2. (/ r (.x v)) (/ r (.y v))))))))

(defn invtree
  "InvTree"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (cond
       (r/brand 0.333) (v/mult v (* 0.5 amount))
       (r/brand 0.666) (v/mult (Vec2. (/ (inc (.x v)))
                                      (/ (.y v) (inc (.y v)))) amount)
       :else (v/mult (Vec2. (/ (.x v) (inc (.x v)))
                            (/ (inc (.y v)))) amount)))))

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
     (let [xs (* (.x v) (.x v))
           ys (* (.y v) (.y v))
           z (+ x2y2-add (m/pow (+ xs (* ys y2-mult)) power))
           q (+ (* 0.5 (m/atan2 (* (.x v) a2x-mult) (+ (* (.y v) a2y-mult) a2y-add)))
                (* m/PI (r/irand 2)))
           r (* amount (m/cos (* z cos-mult)) (/ (* z center)))]
       (Vec2. (* r (m/sin q)) (* r (m/cos q) y-mult))))))

(defn invpolar
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [ny (* amount (inc (.y v)))
           nx (* (.x v) m/PI)]
       (Vec2. (* ny (m/sin nx))
              (* ny (m/cos ny)))))))
