(ns fastmath.fields.a
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn acosech
  "Acosech"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (-> v c/acsch c/flip (c/scale (* amount m/M_2_PI)))]
       (r/randval z (c/neg z))))))

(defn acosh
  "Acosh"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (-> v c/acosh (c/scale (* amount m/M_2_PI)))]
       (r/randval z (c/neg z))))))

(defn acoth
  "Acoth"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> v c/acoth c/flip (c/scale (* amount m/M_2_PI))))))

(defn anamorphcyl
  "AnamorphCyl"
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.01 2.0)
                      :b (r/drand -2.0 2.0)
                      :k (u/sdrand 0.01 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double k]}]
   (fn [^Vec2 v]
     (let [by (* a amount (+ b (.y v)))
           kx (* k (.x v))]
       (Vec2. (* by (m/cos kx))
              (* by (m/sin kx)))))))

(def ^:private ^:const ^double apocarpet_r (/ (inc m/SQRT2)))
(def ^:private apocarpet+r-r (Vec2. apocarpet_r (- apocarpet_r)))
(def ^:private apocarpet-r+r (Vec2. (- apocarpet_r) apocarpet_r))

(defn apocarpet
  "Apolonian carpet"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [denom (v/magsq v)]
       (-> (case (unchecked-int (r/irand 6))
             0 (v/div (Vec2. (* 2.0 (.x v) (.y v))
                             (- (m/sq (.x v)) (m/sq (.y v)))) denom)
             1 (v/shift (v/mult v apocarpet_r) (- apocarpet_r))
             2 (v/shift (v/mult v apocarpet_r) apocarpet_r)
             3 (v/add (v/mult v apocarpet_r) apocarpet+r-r)
             4 (v/add (v/mult v apocarpet_r) apocarpet-r+r)
             5 (v/div (Vec2. (- (m/sq (.x v)) (m/sq (.y v)))
                             (* 2.0 (.x v) (.y v))) denom))
           (v/mult amount))))))

(def ^:private ^:const ^double apollony-inc-sqrt3 (inc m/SQRT3))
(def ^:private ^:const ^double apollony-inc-sqrt3-div (/ apollony-inc-sqrt3 (inc apollony-inc-sqrt3)))

(defn apollony
  "Apollony"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sqy (m/sq (.y v))
           sx (- apollony-inc-sqrt3 (.x v))
           p (+ (m/sq sx) sqy)
           a (- (/ (* 3.0 sx) p)
                apollony-inc-sqrt3-div)
           b (/ (* 3.0 (.y v)) p)
           r (mod (long (r/drand 4.0)) 3)]
       (-> (if (zero? r)
             (Vec2. a b)
             (let [denom (+ (* a a) (* b b))
                   f1x (/ a denom)
                   f1y (/ (- b) denom)]
               (-> (if (m/one? r)
                     (Vec2. (- (- f1x ) (* m/SQRT3 f1y))
                            (- (* m/SQRT3 f1x) f1y))
                     (Vec2. (+ (- f1x ) (* m/SQRT3 f1y))
                            (- (- (* m/SQRT3 f1x)) f1y)))
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
                        :tiles-per-row (r/irand 4 (/ 3.0 radius))
                        :tiles-per-col (r/irand 4 (/ 3.0 radius))}))})
  ([^double amount {:keys [seed ^double thickness ^double radius
                           ^long tiles-per-row ^long tiles-per-col legacy?]}]
   (let [rng (r/rng :well19937c seed)
         tilesize (* 2.0 radius)
         number-of-tiles (* tiles-per-row tiles-per-col)
         tilt-array (vec (repeatedly number-of-tiles #(m/radians (* 90 (r/lrandom rng 4)))))
         r+t (+ radius thickness)
         hthickness (* 0.5 thickness)
         gamma (/ (* thickness (+ tilesize thickness)) r+t)
         shift (Vec2. (* 0.5 tilesize tiles-per-row)
                      (* 0.5 tilesize tiles-per-col))
         rfun (if legacy? #(- r+t (r/drand gamma)) #(- radius (r/drand (- hthickness) hthickness)))]
     (fn [_]
       (let [side (r/lrand 2)
             phi1 (if (zero? side) 0.0 m/PI)
             i (r/lrandom rng tiles-per-row)
             j (r/lrandom rng tiles-per-col)
             ^double tilt (tilt-array (+ i (* j tiles-per-row)))
             ^double r (rfun)
             phi (+ phi1 (r/drand m/HALF_PI))
             p (v/rotate (Vec2. (* r (m/cos phi)) (* r (m/sin phi))) tilt)
             radio (if (zero? side) radius (- radius))
             vradio (v/rotate (Vec2. radio radio) tilt)]
         (-> p
             (v/sub vradio)
             (v/add (Vec2. (+ radius (* i tilesize)) (+ radius (* j tilesize))))
             (v/sub shift)
             (v/mult amount)))))))

(defn arch
  "Arch"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [ang (* amount (r/drand m/PI))
           sinr (m/sin ang)
           cosr (m/cos ang)]
       (if (zero? cosr)
         u/zerov
         (Vec2. (* amount sinr)
                (* amount (/ (m/sq sinr) cosr))))))))

(defn arcsech2
  "Arcsech2"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z (c/reciprocal v)
           z2 (-> z (c/sub c/ONE) c/sqrt)
           z3 (-> z (c/add c/ONE) c/sqrt (c/mult z2))
           z (-> z (c/add z3) c/log (c/scale (* amount m/M_2_PI)))]
       (if (neg? (c/im z))
         (Vec2. (c/re z) (inc (c/im z)))
         (Vec2. (- (c/re z)) (dec (c/im z))))))))

(defn arcsinh
  "Arcsinh"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> v c/asinh (c/scale (* amount m/M_2_PI))))))

(def ^:deprecated asinh arcsinh)

(defn arctanh
  "Arctanh"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [z2 (-> v (c/scale -1.0) (c/add c/ONE))]
       (-> v (c/add c/ONE) (c/div z2) c/log (c/scale (* amount m/M_2_PI)))))))

(def ^:deprecated atanh arctanh)

(defn asteria
  "asteria by dark-beam, http://jwildfire.org/forum/viewtopic.php?f=23&t=1464"
  ([] {:type :random
       :config (fn [] {:alpha (r/drand -5 5)})})
  ([^double amount {:keys [^double alpha]}]
   (let [sina (m/sin (* m/PI alpha))
         cosa (m/cos (* m/PI alpha))]
     (fn [^Vec2 v]
       (let [^Vec2 v0 (v/mult v amount)
             r (v/magsq v0)
             xx (m/sq (dec (m/abs (.x v0))))
             yy (m/sq (dec (m/abs (.y v0))))
             r2 (m/sqrt (+ xx yy))
             in (< r 1.0)
             in1 (if (and (< r 1.0) (< r2 1.0))
                   (r/brand 0.65)
                   (not in))]
         (if in1
           v0
           (let [xx (- (* cosa (.x v0))
                       (* sina (.y v0)))
                 yy (+ (* sina (.x v0))
                       (* cosa (.y v0)))
                 nx (* (/ xx (m/sqrt (- 1.0 (* yy yy))))
                       (- 1.0 (m/sqrt (- 1.0 (m/sq (inc (- (m/abs yy))))))))]
             (Vec2. (+ (* cosa nx)
                       (* sina yy))
                    (+ (- (* sina nx))
                       (* cosa yy))))))))))

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
           fx (- (+ x-add (* x-mult (m/atan2 (+ r-add (* r r-mult))
                                             (+ xy2-add (* xy2 xy2-mult))))) m/PI)
           fy (* y-mult (m/sin (+ sin-add (m/atan2 (+ yy-add (/ (.y v) yy-div))
                                                   (+ yx-add (/ (.x v) yx-div))))))]
       (v/mult (Vec2. (if (neg? (.x v)) fx (- fx)) fy) amount)))))

(defn atan
  "Atan"
  ([] {:type :regular
       :config (fn [] {:mode (r/lrand 3)
                      :stretch (u/sdrand 0.01 2.0)})})
  ([^double amount {:keys [^long mode ^double stretch]}]
   (let [m (m/constrain mode 0 2)
         norm (/ 1.0 (* m/M_PI_2 amount))]
     (fn [^Vec2 v]
       (case m
         0 (Vec2. (.x v)
                  (* norm (m/atan (* stretch (.y v)))))
         1 (Vec2. (* norm (m/atan (* stretch (.x v))))
                  (.y v))
         (Vec2. (* norm (m/atan (* stretch (.x v))))
                (* norm (m/atan (* stretch (.y v))))))))))

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
           s (m/sin (* freq x))
           t (m/sin (* freq y))
           dy (+ y (* weight (+ (m/abs y) (* 0.5 s scale)) s))
           dx (+ x (* weight (+ (m/abs x) (* 0.5 t scale)) t))
           xx (* amount (+ x (* sym (- dx x))))
           yy (* amount dy)]
       (Vec2. xx yy)))))
