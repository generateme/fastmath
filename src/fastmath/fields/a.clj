(ns fastmath.fields.a
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn acosech
  "Acosech"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (-> v c/acsch c/flip (c/scale (m/* amount m/M_2_PI)))]
       (r/randval z (c/neg z))))))

(defn acosh
  "Acosh"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (-> v c/acosh (c/scale (m/* amount m/M_2_PI)))]
       (r/randval z (c/neg z))))))

(defn acoth
  "Acoth"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> v c/acoth c/flip (c/scale (m/* amount m/M_2_PI))))))

(defn anamorphcyl
  "AnamorphCyl"
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.01 2.0)
                      :b (r/drand -2.0 2.0)
                      :k (u/sdrand 0.01 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double k]}]
   (fn [^Vec2 v]
     (let [by (m/* a amount (m/+ b (.y v)))
           kx (m/* k (.x v))]
       (Vec2. (m/* by (m/cos kx))
              (m/* by (m/sin kx)))))))

(def ^{:private true :const true :tag 'double} apocarpet_r (m// (m/inc m/SQRT2)))
(def ^:private apocarpet+r-r (Vec2. apocarpet_r (m/- apocarpet_r)))
(def ^:private apocarpet-r+r (Vec2. (m/- apocarpet_r) apocarpet_r))

(defn apocarpet
  "Apolonian carpet"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [denom (v/magsq v)]
       (-> (case (unchecked-int (r/irand 6))
             0 (v/div (Vec2. (m/* 2.0 (.x v) (.y v))
                             (m/- (m/sq (.x v)) (m/sq (.y v)))) denom)
             1 (v/shift (v/mult v apocarpet_r) (m/- apocarpet_r))
             2 (v/shift (v/mult v apocarpet_r) apocarpet_r)
             3 (v/add (v/mult v apocarpet_r) apocarpet+r-r)
             4 (v/add (v/mult v apocarpet_r) apocarpet-r+r)
             5 (v/div (Vec2. (m/- (m/sq (.x v)) (m/sq (.y v)))
                             (m/* 2.0 (.x v) (.y v))) denom))
           (v/mult amount))))))

(def ^{:private true :const true :tag 'double} apollony-inc-sqrt3 (m/inc m/SQRT3))
(def ^{:private true :const true :tag 'double} apollony-inc-sqrt3-div (m// apollony-inc-sqrt3 (m/inc apollony-inc-sqrt3)))

(defn apollony
  "Apollony"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sqy (m/sq (.y v))
           sx (m/- apollony-inc-sqrt3 (.x v))
           p (m/+ (m/sq sx) sqy)
           a (m/- (m// (m/* 3.0 sx) p)
                apollony-inc-sqrt3-div)
           b (m// (m/* 3.0 (.y v)) p)
           r (m/mod (long (r/drand 4.0)) 3)]
       (-> (if (m/zero? r)
             (Vec2. a b)
             (let [denom (m/+ (m/* a a) (m/* b b))
                   f1x (m// a denom)
                   f1y (m// (m/- b) denom)]
               (-> (if (m/one? r)
                     (Vec2. (m/- (m/- f1x ) (m/* m/SQRT3 f1y))
                            (m/- (m/* m/SQRT3 f1x) f1y))
                     (Vec2. (m/+ (m/- f1x ) (m/* m/SQRT3 f1y))
                            (m/- (m/- (m/* m/SQRT3 f1x)) f1y)))
                   (v/mult 0.5))))
           (v/mult amount))))))

(defn arctruchet
  "Arc Truchet"
  ([] {:type :pattern
       :config (fn [] (let [radius (m/sq (r/drand 0.25 0.99))]
                       {:seed (r/irand)
                        :thickness (r/drand 0.01 0.99)
                        :radius radius
                        :legacy? (r/brand)
                        :tiles-per-row (r/irand 4 (m// 3.0 radius))
                        :tiles-per-col (r/irand 4 (m// 3.0 radius))}))})
  ([^double amount {:keys [seed ^double thickness ^double radius
                           ^long tiles-per-row ^long tiles-per-col legacy?]}]
   (let [rng (r/rng :well19937c seed)
         tilesize (m/* 2.0 radius)
         number-of-tiles (m/long-mult tiles-per-row tiles-per-col)
         tilt-array (vec (repeatedly number-of-tiles #(m/radians (m/* 90 (r/lrandom rng 4)))))
         r+t (m/+ radius thickness)
         hthickness (m/* 0.5 thickness)
         gamma (m// (m/* thickness (m/+ tilesize thickness)) r+t)
         shift (Vec2. (m/* 0.5 tilesize tiles-per-row)
                      (m/* 0.5 tilesize tiles-per-col))
         rfun (if legacy? #(m/- r+t (r/drand gamma)) #(m/- radius (r/drand (m/- hthickness) hthickness)))]
     (fn [_]
       (let [side (r/lrand 2)
             phi1 (if (m/zero? side) 0.0 m/PI)
             i (r/lrandom rng tiles-per-row)
             j (r/lrandom rng tiles-per-col)
             ^double tilt (tilt-array (m/+ i (m/* j tiles-per-row)))
             ^double r (rfun)
             phi (m/+ phi1 (r/drand m/HALF_PI))
             p (v/rotate (Vec2. (m/* r (m/cos phi)) (m/* r (m/sin phi))) tilt)
             radio (if (m/zero? side) radius (m/- radius))
             vradio (v/rotate (Vec2. radio radio) tilt)]
         (-> p
             (v/sub vradio)
             (v/add (Vec2. (m/+ radius (m/* i tilesize)) (m/+ radius (m/* j tilesize))))
             (v/sub shift)
             (v/mult amount)))))))

(defn arch
  "Arch"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [ang (m/* amount (r/drand m/PI))
           sinr (m/sin ang)
           cosr (m/cos ang)]
       (if (m/zero? cosr)
         u/zerov
         (Vec2. (m/* amount sinr)
                (m/* amount (m// (m/sq sinr) cosr))))))))

(defn arcsech2
  "Arcsech2"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (c/reciprocal v)
           z2 (-> z (c/sub c/ONE) c/sqrt)
           z3 (-> z (c/add c/ONE) c/sqrt (c/mult z2))
           z (-> z (c/add z3) c/log (c/scale (m/* amount m/M_2_PI)))]
       (if (m/neg? (c/im z))
         (Vec2. (c/re z) (m/inc (c/im z)))
         (Vec2. (m/- (c/re z)) (m/dec (c/im z))))))))

(defn arcsinh
  "Arcsinh"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> v c/asinh (c/scale (m/* amount m/M_2_PI))))))

(def ^:deprecated asinh arcsinh)

(defn arctanh
  "Arctanh"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z2 (-> v (c/scale -1.0) (c/add c/ONE))]
       (-> v (c/add c/ONE) (c/div z2) c/log (c/scale (m/* amount m/M_2_PI)))))))

(def ^:deprecated atanh arctanh)

(defn asteria
  "asteria by dark-beam, http://jwildfire.org/forum/viewtopic.php?f=23&t=1464"
  ([] {:type :random
       :config (fn [] {:alpha (r/drand -5 5)})})
  ([^double amount {:keys [^double alpha]}]
   (let [sina (m/sin (m/* m/PI alpha))
         cosa (m/cos (m/* m/PI alpha))]
     (fn [^Vec2 v]
       (let [^Vec2 v0 (v/mult v amount)
             r (v/magsq v0)
             xx (m/sq (m/dec (m/abs (.x v0))))
             yy (m/sq (m/dec (m/abs (.y v0))))
             r2 (m/sqrt (m/+ xx yy))
             in (m/< r 1.0)
             in1 (if (and (m/< r 1.0) (m/< r2 1.0))
                   (r/brand 0.65)
                   (not in))]
         (if in1
           v0
           (let [xx (m/- (m/* cosa (.x v0))
                       (m/* sina (.y v0)))
                 yy (m/+ (m/* sina (.x v0))
                       (m/* cosa (.y v0)))
                 nx (m/* (m// xx (m/sqrt (m/- 1.0 (m/* yy yy))))
                       (m/- 1.0 (m/sqrt (m/- 1.0 (m/sq (m/inc (m/- (m/abs yy))))))))]
             (Vec2. (m/+ (m/* cosa nx)
                       (m/* sina yy))
                    (m/+ (m/- (m/* sina nx))
                       (m/* cosa yy))))))))))

(defn atan2spirals
  "Atan2 Spirals"
  ([] {:type :regular
       :config (fn [] {:r-mult (u/sdrand 0.01 2.0)
                      :r-add (r/drand 0.9 1.1)
                      :xy2-mult (u/sdrand 0.1 1.5)
                      :xy2-add (u/sdrand 0.1 3.0)
                      :x-mult (r/drand 0.7 3.0)
                      :x-add (r/drand -0.5 3.0)
                      :yx-div (u/sdrand 0.1 2.0)
                      :yx-add (r/drand -0.5 0.5)
                      :yy-div (u/sdrand 0.1 2.0)
                      :yy-add (r/drand -0.5 0.5)
                      :sin-add (r/drand m/-TWO_PI m/TWO_PI)
                      :y-mult (u/sdrand 0.9 2.0)
                      :r-power (r/drand 0.05 0.6)
                      :x2y2-power (r/drand 0.9 4.0)})})
  ([^double amount {:keys [^double r-mult ^double r-add ^double xy2-mult ^double xy2-add
                           ^double x-mult ^double x-add ^double yx-div ^double yx-add
                           ^double yy-div ^double yy-add ^double sin-add ^double y-mult
                           ^double r-power ^double x2y2-power]}]
   (fn [^Vec2 v]
     (let [xs+ys (v/magsq v)
           xy2 (m/pow xs+ys x2y2-power)
           r (m/pow xs+ys r-power)
           fx (m/- (m/+ x-add (m/* x-mult (m/atan2 (m/+ r-add (m/* r r-mult))
                                             (m/+ xy2-add (m/* xy2 xy2-mult))))) m/PI)
           fy (m/* y-mult (m/sin (m/+ sin-add (m/atan2 (m/+ yy-add (m// (.y v) yy-div))
                                                   (m/+ yx-add (m// (.x v) yx-div))))))]
       (v/mult (Vec2. (if (m/neg? (.x v)) fx (m/- fx)) fy) amount)))))

(defn atan
  "Atan"
  ([] {:type :regular
       :config (fn [] {:mode (r/lrand 3)
                      :stretch (u/sdrand 0.01 2.0)})})
  ([^double amount {:keys [^long mode ^double stretch]}]
   (let [m (m/constrain mode 0 2)
         norm (m// 1.0 (m/* m/M_PI_2 amount))]
     (fn [^Vec2 v]
       (case m
         0 (Vec2. (.x v)
                  (m/* norm (m/atan (m/* stretch (.y v)))))
         1 (Vec2. (m/* norm (m/atan (m/* stretch (.x v))))
                  (.y v))
         (Vec2. (m/* norm (m/atan (m/* stretch (.x v))))
                (m/* norm (m/atan (m/* stretch (.y v))))))))))

(defn auger
  "Auger by Xyrus02"
  ([] {:type :regular
       :config (fn [] {:freq (r/drand -5.0 5.0)
                      :weight (r/drand -1.0 1.0)
                      :sym (r/drand -2.0 2.0)
                      :scale (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double freq ^double weight ^double sym ^double scale]}]
   (fn [^Vec2 v]
     (let [x (.x v)
           y (.y v)
           s (m/sin (m/* freq x))
           t (m/sin (m/* freq y))
           dy (m/+ y (m/* weight (m/+ (m/abs y) (m/* 0.5 s scale)) s))
           dx (m/+ x (m/* weight (m/+ (m/abs x) (m/* 0.5 t scale)) t))
           xx (m/* amount (m/+ x (m/* sym (m/- dx x))))
           yy (m/* amount dy)]
       (Vec2. xx yy)))))

