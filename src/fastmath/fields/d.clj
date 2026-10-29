(ns fastmath.fields.d
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]
           [fastmath.java Array]))


(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn dlawf
  ([] {:type :pattern
       :config (fn [] {:buffer-size (r/irand 200 1000)
                      :max-iter (r/irand 1000 10000)
                      :seed (r/irand)
                      :scale (r/drand 5.0 15.0)
                      :jitter (m/sq (r/drand 0.15))})})
  ([^double amount {:keys [^long buffer-size ^long max-iter ^long seed
                           ^double scale ^double jitter]}]
   (let [jitter-radius (m/constrain jitter 0.0 1.0)
         rng (r/rng :mersenne seed)
         centre (m// buffer-size 2)
         size2 (m/- buffer-size 2)
         bss (m// scale buffer-size)
         q (int-array (m/* buffer-size buffer-size))]
     (Array/set2d q buffer-size centre centre 1)
     (loop [i (long 0)
            ^Vec2 rr (Vec2. 3.0 9.0)]
       (when (m/< i max-iter)
         (let [phi (r/drandom rng m/TWO_PI)
               ri (m/* (.x rr) (m/cos phi))
               rj (m/* (.x rr) (m/sin phi))
               ci (m/+ centre (long (m/+ ri 0.5)))
               cj (m/+ centre (long (m/+ rj 0.5)))
               ^double r (loop [nci ci ncj cj]
                           (if (or (m/< nci 1) (m/> nci size2)
                                   (m/< ncj 1) (m/> ncj size2))
                             -1.0
                             (let [idx (m/+ nci (m/* buffer-size ncj))
                                   sum (m/+ (Array/get q (m/dec idx))
                                          (Array/get q (m/inc idx))
                                          (Array/get q (m/- idx buffer-size))
                                          (Array/get q (m/+ idx buffer-size)))
                                   r3 (m/dist centre centre nci ncj)]
                               (if-not (m/zero? sum)
                                 (do
                                   (Array/set q idx 1)
                                   (m/max r3 (.x rr)))
                                 (if (m/> r3 (.y rr))
                                   -1.0
                                   (case (r/irandom rng 4)
                                     0 (recur (m/inc nci) ncj)
                                     1 (recur (m/dec nci) ncj)
                                     2 (recur nci (m/inc ncj))
                                     3 (recur nci (m/dec ncj))))))))]
           (if (m/neg? r)
             (recur i rr)
             (recur (m/inc i) (Vec2. r (m/* 2.1 r)))))))
     (let [r (range buffer-size)
           points (vec (for [^long i r ^long j r
                             :when (not (m/zero? (Array/get2d q buffer-size i j)))
                             :let [arnd (if (m/zero? jitter-radius) 0.0 (r/drandom rng))
                                   v (Vec2. (m/* (m/- i centre) bss)
                                            (m/* (m/- j centre) bss))]]
                         (v/mult (if (m/zero? jitter-radius)
                                   v
                                   (let [alpha (m/* arnd m/TWO_PI)]
                                     (v/add v (Vec2. (m/* jitter-radius (m/cos alpha))
                                                     (m/* jitter-radius (m/sin alpha)))))) amount)))
           cnt (count points)]
       (fn [_]
         (points (r/irand cnt)))))))

(defn dspherical
  ([] {:type :random
       :config (fn [] {:d-spher-weight (r/drand)})})
  ([^double amount {:keys [^double d-spher-weight]}]
   (fn [^Vec2 v]
     (if (m/< (r/drand) d-spher-weight)
       (let [r (m// amount (v/magsq v))]
         (v/mult v r))
       (v/mult v amount)))))

(def ^:private deltaa-unit-x (Vec2. 1.0 0.0))

(defn deltaa
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [v- (v/sub v deltaa-unit-x)
           v+ (v/add v deltaa-unit-x)
           avgr (m/* amount (m// (v/mag v+) (v/mag v-)))
           avga (m/* 0.5 (m/- (v/heading v-) (v/heading v+)))]
       (Vec2. (m/* avgr (m/cos avga))
              (m/* avgr (m/sin avga)))))))

(defn devilwarp
  ([] {:type :regular
       :config (fn [] (let [rmin (r/drand -20.0 0.0)]
                       {:a (r/drand 0.1 3.0)
                        :b (r/drand 0.1 3.0)
                        :effect (u/sdrand 0.2 2.0)
                        :warp (u/sdrand 0.1 3.0)
                        :rmin rmin
                        :rmax (r/drand (m/inc rmin) (m/+ rmin 40.0))}))})
  ([^double amount {:keys [^double a ^double b ^double effect
                           ^double warp ^double rmin ^double rmax]}]
   (let [rmn (m/min rmin rmax)
         rmx (m/max rmin rmax)]
     (fn [^Vec2 v]
       (let [r2 (m// (v/magsq v))
             r (-> (m/- (m/pow (m/+ (m/sq (.x v)) (m/* (m/sq (.y v)) r2 b)) warp)
                      (m/pow (m/+ (m/sq (.y v)) (m/* (m/sq (.x v)) r2 a)) warp))
                   (m/constrain rmn rmx)
                   (m/* effect)
                   (m/inc)
                   (m/* amount))]
         (v/mult v r))))))

(defn diamond
  "Diamond"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [length (v/mag v)
           sina (m// (.x v) length)
           cosa (m// (.y v) length)
           sinr (m/sin length)
           cosr (m/cos length)]
       (Vec2. (m/* amount sina cosr)
              (m/* amount cosa sinr))))))

(defn disc2
  "Disc2"
  ([] {:type :regular
       :config (fn [] {:rot (r/drand -3.0 3.0)
                      :twist (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double rot ^double twist]}]
   (let [timespi (m/* rot m/PI) 
         k1 (if (m/> twist m/TWO_PI) (m/- (m/inc twist) m/TWO_PI) 1.0)
         k2 (if (m/< twist (m/- m/TWO_PI)) (m/+ (m/inc twist) m/TWO_PI) 1.0) 
         sinadd (m/* (m/sin twist) k1 k2)
         cosadd (m/* (m/dec (m/cos twist)) k1 k2)]
     (fn [^Vec2 v]
       (let [t (m/* timespi (m/+ (.x v) (.y v)))
             sinr (m/sin t)
             cosr (m/cos t)
             r (m// (m/* amount (v/heading v)) m/PI)]
         (Vec2. (m/* r (m/+ cosadd sinr))
                (m/* r (m/+ sinadd cosr))))))))

(defn disc3
  "Disc"
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.5 1.5)
                      :b (u/sdrand 0.5 1.5)
                      :c (u/sdrand 0.5 2.0)
                      :d (r/drand -2.0 2.0)
                      :e (r/drand -2.0 2.0)
                      :f (r/drand -2.0 2.0)
                      :g (r/drand -2.0 2.0)
                      :h (u/sdrand 0.5 1.5)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d
                           ^double e ^double f ^double g ^double h]}]
   (let [api (m/* c (m// amount m/PI))
         de (m/abs (m/* d e))
         fg (m/abs (m/* f g))]
     (fn [^Vec2 v]
       (let [rpi (m/* m/PI (m/sqrt (m/+ (m/* (.x v) (.x v) de)
                                    (m/* (.y v) (.y v) fg))))
             sinr (m/* a (m/sin rpi))
             cosr (m/* b (m/cos rpi))
             r (m/* api (v/heading v))]
         (Vec2. (m/* r sinr h) (m/* r cosr h)))))))

(defn disc
  "Disc"
  ([] {:type :regular})
  ([^double amount _]
   (let [api (m// amount m/PI)]
     (fn [^Vec2 v]
       (let [rpi (m/* m/PI (v/mag v))
             sinr (m/sin rpi)
             cosr (m/cos rpi)
             r (m/* api (v/heading v))]
         (Vec2. (m/* r sinr) (m/* r cosr)))))))

(defn dustpoint
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [p (r/randval 1.0 -1.0)
           r (v/mag v)
           w (r/drand)]
       (v/mult (cond
                 (m/< w 0.5) (Vec2. (m/dec (m// (.x v) r))
                                  (m/* p (m// (.y v) r)))
                 (m/< w 0.75) (v/div v 3.0)
                 :else (Vec2. (m/+ m/TWO_THIRD (m// (.x v) 3.0))
                              (m// (.y v) 3.0))) amount)))))

