(ns fastmath.fields.d
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]
           [fastmath.java Array]))


(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
         centre (/ buffer-size 2)
         size2 (- buffer-size 2)
         bss (/ scale buffer-size)
         q (int-array (* buffer-size buffer-size))]
     (Array/set2d q buffer-size centre centre 1)
     (loop [i (long 0)
            ^Vec2 rr (Vec2. 3.0 9.0)]
       (when (< i max-iter)
         (let [phi (r/drandom rng m/TWO_PI)
               ri (* (.x rr) (m/cos phi))
               rj (* (.x rr) (m/sin phi))
               ci (+ centre (long (+ ri 0.5)))
               cj (+ centre (long (+ rj 0.5)))
               ^double r (loop [nci ci ncj cj]
                           (if (or (< nci 1) (> nci size2)
                                   (< ncj 1) (> ncj size2))
                             -1.0
                             (let [idx (+ nci (* buffer-size ncj))
                                   sum (+ (Array/get q (dec idx))
                                          (Array/get q (inc idx))
                                          (Array/get q (- idx buffer-size))
                                          (Array/get q (+ idx buffer-size)))
                                   r3 (m/dist centre centre nci ncj)]
                               (if-not (zero? sum)
                                 (do
                                   (Array/set q idx 1)
                                   (max r3 (.x rr)))
                                 (if (> r3 (.y rr))
                                   -1.0
                                   (case (r/irandom rng 4)
                                     0 (recur (inc nci) ncj)
                                     1 (recur (dec nci) ncj)
                                     2 (recur nci (inc ncj))
                                     3 (recur nci (dec ncj))))))))]
           (if (neg? r)
             (recur i rr)
             (recur (inc i) (Vec2. r (* 2.1 r)))))))
     (let [r (range buffer-size)
           points (vec (for [^long i r ^long j r
                             :when (not (zero? (Array/get2d q buffer-size i j)))
                             :let [arnd (if (zero? jitter-radius) 0.0 (r/drandom rng))
                                   v (Vec2. (* (- i centre) bss)
                                            (* (- j centre) bss))]]
                         (v/mult (if (zero? jitter-radius)
                                   v
                                   (let [alpha (* arnd m/TWO_PI)]
                                     (v/add v (Vec2. (* jitter-radius (m/cos alpha))
                                                     (* jitter-radius (m/sin alpha)))))) amount)))
           cnt (count points)]
       (fn [_]
         (points (r/irand cnt)))))))

(defn dspherical
  ([] {:type :random
       :config (fn [] {:d-spher-weight (r/drand)})})
  ([^double amount {:keys [^double d-spher-weight]}]
   (fn [^Vec2 v]
     (if (< (r/drand) d-spher-weight)
       (let [r (/ amount (v/magsq v))]
         (v/mult v r))
       (v/mult v amount)))))

(def ^:private deltaa-unit-x (Vec2. 1.0 0.0))

(defn deltaa
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [v- (v/sub v deltaa-unit-x)
           v+ (v/add v deltaa-unit-x)
           avgr (* amount (/ (v/mag v+) (v/mag v-)))
           avga (* 0.5 (- (v/heading v-) (v/heading v+)))]
       (Vec2. (* avgr (m/cos avga))
              (* avgr (m/sin avga)))))))

(defn devilwarp
  ([] {:type :regular
       :config (fn [] (let [rmin (r/drand -20.0 0.0)]
                       {:a (r/drand 0.1 3.0)
                        :b (r/drand 0.1 3.0)
                        :effect (u/sdrand 0.2 2.0)
                        :warp (u/sdrand 0.1 3.0)
                        :rmin rmin
                        :rmax (r/drand (inc rmin) (+ rmin 40.0))}))})
  ([^double amount {:keys [^double a ^double b ^double effect
                           ^double warp ^double rmin ^double rmax]}]
   (let [rmn (min rmin rmax)
         rmx (max rmin rmax)]
     (fn [^Vec2 v]
       (let [r2 (/ (v/magsq v))
             r (-> (- (m/pow (+ (m/sq (.x v)) (* (m/sq (.y v)) r2 b)) warp)
                      (m/pow (+ (m/sq (.y v)) (* (m/sq (.x v)) r2 a)) warp))
                   (m/constrain rmn rmx)
                   (* effect)
                   (inc)
                   (* amount))]
         (v/mult v r))))))

(defn diamond
  "Diamond"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [length (v/mag v)
           sina (/ (.x v) length)
           cosa (/ (.y v) length)
           sinr (m/sin length)
           cosr (m/cos length)]
       (Vec2. (* amount sina cosr)
              (* amount cosa sinr))))))

(defn disc2
  "Disc2"
  ([] {:type :regular
       :config (fn [] {:rot (r/drand -3.0 3.0)
                      :twist (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double rot ^double twist]}]
   (let [timespi (* rot m/PI) 
         k1 (if (> twist m/TWO_PI) (- (inc twist) m/TWO_PI) 1.0)
         k2 (if (< twist (- m/TWO_PI)) (+ (inc twist) m/TWO_PI) 1.0) 
         sinadd (* (m/sin twist) k1 k2)
         cosadd (* (dec (m/cos twist)) k1 k2)]
     (fn [^Vec2 v]
       (let [t (* timespi (+ (.x v) (.y v)))
             sinr (m/sin t)
             cosr (m/cos t)
             r (/ (* amount (v/heading v)) m/PI)]
         (Vec2. (* r (+ cosadd sinr))
                (* r (+ sinadd cosr))))))))

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
   (let [api (* c (/ amount m/PI))
         de (m/abs (* d e))
         fg (m/abs (* f g))]
     (fn [^Vec2 v]
       (let [rpi (* m/PI (m/sqrt (+ (* (.x v) (.x v) de)
                                    (* (.y v) (.y v) fg))))
             sinr (* a (m/sin rpi))
             cosr (* b (m/cos rpi))
             r (* api (v/heading v))]
         (Vec2. (* r sinr h) (* r cosr h)))))))

(defn disc
  "Disc"
  ([] {:type :regular})
  ([^double amount _]
   (let [api (/ amount m/PI)]
     (fn [^Vec2 v]
       (let [rpi (* m/PI (v/mag v))
             sinr (m/sin rpi)
             cosr (m/cos rpi)
             r (* api (v/heading v))]
         (Vec2. (* r sinr) (* r cosr)))))))

(defn dustpoint
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [p (r/randval 1.0 -1.0)
           r (v/mag v)
           w (r/drand)]
       (v/mult (cond
                 (< w 0.5) (Vec2. (dec (/ (.x v) r))
                                  (* p (/ (.y v) r)))
                 (< w 0.75) (v/div v 3.0)
                 :else (Vec2. (+ m/TWO_THIRD (/ (.x v) 3.0))
                              (/ (.y v) 3.0))) amount)))))

(m/unuse-primitive-operators)
