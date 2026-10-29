(ns fastmath.fields.b
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- bseries-calc
  "Common calculations for bSeries "
  [^double amount ^double tau ^double sigma]
  (let [sinht (m/sinh tau)
        cosht (m/cosh tau)
        sins (m/sin sigma)
        coss (m/cos sigma)
        temp (m/- cosht coss)]
    (if (m/zero? temp)
      (Vec2. 0.0 0.0)
      (Vec2. (m/* amount (m// sinht temp))
             (m/* amount (m// sins temp))))))

(defn bcollide
  "bCollide by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:num (u/sdrand 1.0 30.0)
                      :a (r/drand 2.0)})})
  ([^double amount {:keys [^double num ^double a]}]
   (let [bcn-pi (m/* num m/M_1_PI)
         pi-bcn (m// m/PI num)
         bca-bcn (m// (m/* m/PI a) num)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (m/- 1.0 (.x v)) (.y v))
             tau (m/* 0.5 (m/- (m/log (v/magsq v+))
                           (m/log (v/magsq v-))))
             pre-sigma (m/- m/PI (v/heading v+) (v/heading v-))
             alt (double (int (m/* pre-sigma bcn-pi)))
             sigma (if (m/even? (int alt))
                     (m/+ (m/* alt pi-bcn) (m/mod (m/+ pre-sigma bca-bcn) pi-bcn))
                     (m/+ (m/* alt pi-bcn) (m/mod (m/- pre-sigma bca-bcn) pi-bcn)))]
         (bseries-calc amount tau sigma))))))

(defn bmod
  "bMod by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:radius (r/drand 0.5 2.0)
                      :distance (r/drand 2.0)})})
  ([^double amount {:keys [^double radius ^double distance]}]
   (let [rd (m/* radius distance)
         r2 (m/+ radius radius)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (m/- 1.0 (.x v)) (.y v))
             pre-tau (m/* 0.5 (m/- (m/log (v/magsq v+))
                               (m/log (v/magsq v-))))
             sigma (m/- m/PI (v/heading v+) (v/heading v-))
             tau (if (and (m/< pre-tau radius) (m/< (m/- pre-tau) radius))
                   (m/- (m/mod (m/+ pre-tau radius rd) r2) radius)
                   pre-tau)]
         (bseries-calc amount tau sigma))))))

(defn bsplit
  "Raykoid666, transcribed and modded by Nic Anderson, chronologicaldot"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -2.0 2.0)
                      :y (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [xx (m/+ x (.x v))]      
       (Vec2. (m/* (m// amount (m/tan xx)) (m/cos (m/+ y (.y v))))
              (m/* (m// amount (m/sin xx)) (m/- y (.y v))))))))

(defn bswirl
  "bSwirl by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:in (r/drand -2.0 2.0)
                      :out (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double in ^double out]}]
   (fn [^Vec2 v]
     (let [v+ (v/add v u/unitx)
           v- (Vec2. (m/- 1.0 (.x v)) (.y v))
           tau (m/* 0.5 (m/- (m/log (v/magsq v+))
                         (m/log (v/magsq v-))))
           pre-sigma (m/- m/PI (v/heading v+) (v/heading v-))
           sigma (m/+ pre-sigma (m/* tau out) (m// in tau))]
       (bseries-calc amount tau sigma)))))

(defn btransform
  "bTransform by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :random
       :config (fn [] {:rotate (r/drand m/TWO_PI)
                      :power (r/randval (r/drand 10.0) (r/irand 10.0))
                      :move (r/drand -2.0 2.0)
                      :split (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double rotate ^double power ^double move ^double split]}]
   (let [mp (m// m/TWO_PI power)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (m/- 1.0 (.x v)) (.y v))
             pre-tau (m/+ (m// (m/* 0.5 (m/- (m/log (v/magsq v+))
                                     (m/log (v/magsq v-)))) power) move)
             pre-sigma (m/+ (m/- m/PI (v/heading v+) (v/heading v-)) rotate)
             tau (if (m/neg? (.x v))
                   (m/- pre-tau split)
                   (m/+ pre-tau split))
             sigma (m/+ (m// pre-sigma power)
                      (m/* mp (m/floor (m/* power (r/drand)))))]
         (bseries-calc amount tau sigma))))))

(defn- bwrands-bytemix ^long [^long a ^long b] (m/bit-xor (m/bit-and a 0x5a5a) (m/bit-and b 0xa5a5)))
(defn- bwrands-bytexim ^long [^long a ^long b] (m/bit-xor (m/bit-and a 0xaaaa) (m/bit-and b 0x5555)))
(defn- bwrands-byteshf ^long [^long a ^long b] (m/bit-and (m/bit-xor (m/<< a 8) (m/>> b 8)) 0xffff))
(defn- bwrands-byteprimes ^long [^long a ^long b]
  (m/bit-and 0xffff (m/bit-xor (m/long-sub (m/long-mult a 857) 4) (m/long-add (m/long-mult b -977) 8))))

(defn bwrands
  "bwrands"
  ([] {:type :random
       :config (fn [] (let [min-petals (r/lrand 1 5)]
                       {:cellsize (u/sdrand 0.1 2.0)
                        :space (r/drand -1.0 1.0)
                        :gain (r/drand -2.0 2.0)
                        :inner-twist (r/drand -2.0 2.0)
                        :outer-twist (r/drand -2.0 2.0)
                        :seed (r/irand)
                        :rrot (r/drand -1.0 1.0)
                        :rmin (r/drand)
                        :loonie-chance (r/drand 0.5 0.9)
                        :petal-chance (r/drand 0.5 0.9)
                        :min-petals min-petals
                        :max-petals (m/+ min-petals (r/lrand 15))}))})
  ([^double amount {:keys [^double cellsize ^double space ^double gain ^double inner-twist ^double outer-twist
                           ^long seed ^double rrot ^double rmin ^double loonie-chance ^double petal-chance
                           ^long min-petals ^long max-petals]}]
   (let [radius (m/* 0.5 (m// cellsize (m/inc (m/sq space))))
         g2 (m/+ m/EPSILON (m// (m/sq gain) radius))
         max-bubble (as-> (m/* g2 radius) a
                      (if (m/> a 2.0) 1.0 (m/* a (m// (m/inc (m/* 0.25 a a))))))
         r2 (m/* radius radius)
         rfactor (m// radius max-bubble)
         petx min-petals
         pety (m/- max-petals min-petals)]
     (fn [^Vec2 v]
       (if (m/< (m/abs cellsize) m/EPSILON)
         (v/mult v amount)
         (let [^Vec2 I (v/floor (v/div v cellsize))
               C (v/mult (v/shift I 0.5) cellsize)
               xx (m/bit-xor (int (.x I)) 0xb641)
               yy (m/bit-xor (int (.y I)) 0x9d81)
               xy (m/bit-and 0xffff (m/long-add seed (m/long-mult xx yy)))
               xx (m/bit-and 0xffff xx)
               yy (m/bit-and 0xffff yy)
               tt (bwrands-bytemix xx yy)
               yy (bwrands-bytemix yy xx)
               xx tt
               tt (bwrands-byteshf xx yy)
               yy (bwrands-byteshf xy yy)
               xx tt
               tt (bwrands-bytexim xx yy)
               yy (bwrands-bytexim yy xx)
               xx tt
               ssz (-> (m// xx 65536.0)
                       (m/* (m/- 1.0 rmin))
                       (m/+ rmin))
               aan (m// (m/* rrot m/TWO_PI yy) 65536.0)
               tt (bwrands-byteprimes xx yy)
               yy (bwrands-byteprimes yy xx)
               xx tt
               LC (m/+ (m// xx -65536.0) loonie-chance)
               PC (if (m/neg? LC) (m/+ LC petal-chance) 0.0)
               LC (if (m/<= LC PC) -1.0 LC)
               L (v/sub v C)
               vv2 (m/* ssz r2)
               mL (v/magsq L)]
           (if (m/> mL vv2)
             (v/mult v amount)
             (let [nL (if (m/pos? LC)
                        (v/mult L (m/sqrt (m/dec (m// vv2 mL))))
                        (if (m/pos? PC)
                          (let [NPetals (if (m/zero? pety)
                                          petx
                                          (m/+ petx (m/mod (m/>> yy 3) (m/inc pety))))
                                flrw (as-> (m/* (m/+ m/PI (v/heading L))
                                              (m// NPetals m/TWO_PI)) a
                                       (m/- a (int a))
                                       (m/* 2.0 (m/abs (m/- a 0.5))))
                                r (m/sqrt mL)]
                            (r/randval (m/* 0.5 (m/+ flrw 0.5))
                                       (v/mult L (m/* (m/- 1.0 r) (m/* 1.1 flrw)))
                                       (v/mult L (m// (m/- (m/sqrt vv2)
                                                       (m/* r (m/- 1.0 flrw))) (m/+ r m/EPSILON)))))
                          (let [nL (v/mult L g2)
                                r (m// rfactor (m/inc (m// (v/magsq nL) (m/* 4.0 ssz))))]
                            (v/mult nL r))))
                   nvvl (if-not (or (m/pos? LC) (m/pos? PC)) (m/sqrt ssz) 1.0)
                   r (m// (m/* nvvl (v/mag nL)) r2)
                   theta (m/lerp inner-twist outer-twist r)]
               (v/mult (v/add C (v/rotate nL (m/+ aan theta))) amount)))))))))

(defn bwraps7
  "http://slobo777.deviantart.com/art/Bubble-Wrap-WIP-Plugin-112370125"
  ([] {:type :regular
       :config (fn [] {:cellsize (u/sdrand 0.5 2.0)
                      :space (r/drand -1.0 1.0)
                      :gain (r/drand -2.0 2.0)
                      :inner-twist (r/drand -2.0 2.0)
                      :outer-twist (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double cellsize ^double space ^double gain
                           ^double inner-twist ^double outer-twist]}]
   (let [radius (m/* 0.5 (m// cellsize (m/inc (m/sq space))))
         g2 (m/+ m/EPSILON (m/* gain gain))
         max-bubble- (m/* g2 radius)
         max-bubble (if (m/> max-bubble- 2.0)
                      1.0
                      (m/* max-bubble- (m// 1.0 (m/inc (m/* 0.25 (m/sq max-bubble-))))))
         r2 (m/sq radius)
         rfactor (m// radius max-bubble)]
     (fn [^Vec2 v]
       (if (m/< (m/abs cellsize) m/EPSILON)
         (v/mult v amount)
         (let [^Vec2 C (-> v
                           (v/div cellsize)
                           (v/fmap m/floor)
                           (v/add (Vec2. 0.5 0.5))
                           (v/mult cellsize))
               L (v/sub v C)]
           (if (m/> (v/magsq L) r2)
             (v/mult v amount)
             (let [L (v/mult L g2)
                   r (m// rfactor (m/inc (m/* 0.25 (v/magsq L))))
                   ^Vec2 L (v/mult L r)
                   r (m// (v/magsq L) r2)
                   theta (m/+ (m/* inner-twist (m/- 1.0 r))
                            (m/* outer-twist r))
                   s (m/sin theta)
                   c (m/cos theta)
                   vx (m/+ (.x C) (m/* c (.x L)) (m/* s (.y L)))
                   vy (m/+ (.y C) (m/* -1.0 s (.x L)) (m/* c (.y L)))]
               (v/mult (Vec2. vx vy) amount)))))))))

(defn barycentroid
  "barycentroid from Xyrus02, http://xyrusworx.deviantart.com/art/Barycentroid-Plugin-144832371?q=sort%3Atime+favby%3Amistywisp&qo=0&offset=10"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (let [dot00 (m/+ (m/* a a) (m/* b b))
         dot01 (m/+ (m/* a c) (m/* b d))
         dot11 (m/+ (m/* c c) (m/* d d))
         inv-denom (m// 1.0 (m/- (m/* dot00 dot11) (m/* dot01 dot01)))]
     (fn [^Vec2 v]
       (let [dot02 (m/+ (m/* a (.x v)) (m/* b (.y v)))
             dot12 (m/+ (m/* c (.x v)) (m/* d (.y v)))
             u (m/* inv-denom (m/- (m/* dot11 dot02) (m/* dot01 dot12)))
             vv (m/* inv-denom (m/- (m/* dot00 dot12) (m/* dot01 dot02)))
             um (m/* (m/signum u) (m/sqrt (m/+ (m/* u u) (m/sq (.x v)))))
             vm (m/* (m/signum vv) (m/sqrt (m/+ (m/* vv vv) (m/sq (.y v)))))]
         (Vec2. (m/* amount um)
                (m/* amount vm)))))))

(defn bent
  "Bent"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [nx (if (m/neg? (.x v)) (m/+ (.x v) (.x v)) (.x v))
           ny (if (m/neg? (.y v)) (m/* (.y v) 0.5) (.y v))]
       (Vec2. (m/* amount nx)
              (m/* amount ny))))))

(defn bent2
  "Bent2"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.5 2.0)
                      :y (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [nx (if (m/neg? (.x v)) (m/* (.x v) x) (.x v))
           ny (if (m/neg? (.y v)) (m/* (.y v) y) (.y v))]
       (Vec2. (m/* amount nx)
              (m/* amount ny))))))

(defn bilinear
  "Bilinear"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (.y v))
            (m/* amount (.x v))))))

(defn bipolar
  "Bipolar"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shift]}]
   (let [ps (m/* (m/- m/HALF_PI) shift)]
     (fn [^Vec2 v]
       (let [x2y2 (v/magsq v)
             t (m/inc x2y2)
             x2 (m/+ (.x v) (.x v))
             pre-y (m/+ ps (m/* 0.5 (m/atan2 (m/+ (.y v) (.y v))
                                         (m/dec x2y2))))
             y (if (m/> pre-y m/HALF_PI)
                 (m/- (m/mod (m/+ pre-y m/HALF_PI) m/PI) m/HALF_PI)
                 (if (m/< pre-y (m/- m/HALF_PI))
                   (m/- m/HALF_PI (m/mod (m/- m/HALF_PI pre-y) m/PI))
                   pre-y))
             f (m/+ t x2)
             g (m/- t x2)]
         (if (or (m/zero? g)
                 (not (m/pos? (m// f g))))
           (Vec2. 0.0 0.0)
           (Vec2. (m/* amount m/M_2_PI 0.25 (m/log (m// f g)))
                  (m/* amount m/M_2_PI y))))))))

(defn bipolar2
  "Bipolar2"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -2.0 2.0)
                      :a (r/drand -2.0 2.0)
                      :b (u/sdrand 1.5 3.0)
                      :c (u/sdrand 0.8 2.0)
                      :d (r/drand -1.0 2.0)
                      :e (u/sdrand 0.5 3.0)
                      :f1 (u/sdrand 0.8 3.0)
                      :g1 (r/drand 0.5 2.0)
                      :h (u/sdrand 0.8 2.0)})})
  ([^double amount {:keys [^double shift ^double a ^double b ^double c ^double d
                           ^double e ^double f1 ^double g1 ^double h]}]
   (let [ps (m/* (m/- m/HALF_PI) shift)]
     (fn [^Vec2 v]
       (let [x2y2 (m/* g1 (v/magsq v))
             t (m/+ a x2y2)
             x2 (m/* b (.x v))
             pre-y (m/+ ps (m/* c (m/atan2 (m/* e (.y v))
                                       (m/- x2y2 d))))
             y (if (m/> pre-y m/HALF_PI)
                 (m/- (m/mod (m/+ pre-y m/HALF_PI) m/PI) m/HALF_PI)
                 (if (m/< pre-y (m/- m/HALF_PI))
                   (m/- m/HALF_PI (m/mod (m/- m/HALF_PI pre-y) m/PI))
                   pre-y))
             f (m/+ t x2)
             g (m/- t x2)]
         (if (or (m/zero? g)
                 (not (m/pos? (m// f g))))
           (Vec2. 0.0 0.0)
           (Vec2. (m/* amount f1 m/M_2_PI 0.25 (m/log (m// f g)))
                  (m/* amount m/M_2_PI y h))))))))

(defn blade
  "Blade"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* (r/drand amount) (v/mag v))
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (m/* amount (.x v) (m/+ cosr sinr))
              (m/* amount (.x v) (m/- cosr sinr)))))))

(defn blade2
  "Blade2"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* (r/drand amount) (v/mag v))
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (m/* amount (.x v) (m/+ cosr sinr))
              (m/* amount (.y v) (m/- cosr sinr)))))))

(defn blob
  "Blob"
  ([] {:type :regular
       :config (fn [] {:low (u/sdrand 0.1 2.0)
                      :high (u/sdrand 0.1 2.0)
                      :waves (r/randval (r/irand -6 6) (r/drand -6.0 6.0))})})
  ([^double amount {:keys [^double low ^double high ^double waves]}]
   (let [hl (m/- high low)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             rr (->> (m/* a waves)
                     (m/sin)
                     (m/* 0.5)
                     (m/+ 0.5)
                     (m/* hl)
                     (m/+ low)
                     (m/* r))]
         (Vec2. (m/* amount rr (m/sin a))
                (m/* amount rr (m/cos a))))))))

(defn blocky
  "blocky from FracFx, http://fracfx.deviantart.com/art/FracFx-Plugin-Pack-171806681"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.5 1.5)
                      :y (u/sdrand 0.5 1.5)
                      :mp (u/sdrand 0.001 6.0)})})
  ([^double amount {:keys [^double x ^double y ^double mp]}]
   (let [vv (m// amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [T (m/inc (m// (m/+ (m/cos (.x v)) (m/cos (.y v))) mp))
             r (m// amount T)
             tmp (m/inc (v/magsq v))
             x2 (m/+ (.x v) (.x v))
             y2 (m/+ (.y v) (.y v))
             xmax (m/* 0.5 (m/+ (m/sqrt (m/+ tmp x2)) (m/sqrt (m/- tmp x2))))
             ymax (m/* 0.5 (m/+ (m/sqrt (m/+ tmp y2)) (m/sqrt (m/- tmp y2))))
             ax (m// (.x v) xmax)
             bx (m/safe-sqrt (m/- 1.0 (m/sq ax)))
             ay (m// (.y v) ymax)
             by (m/safe-sqrt (m/- 1.0 (m/sq ay)))]
         (Vec2. (m/* vv (m/atan2 ax bx) r x)
                (m/* vv (m/atan2 ay by) r y)))))))

(defn blurcircle
  "Blur circle"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [x (r/drand -1.0 1.0)
           y (r/drand -1.0 1.0)
           absx (m/abs x)
           absy (m/abs y)
           ^Vec2 ps (if (m/>= absx absy)
                      (Vec2. (if (m/>= x absy)
                               (m/+ absx y)
                               (m/- (m/* 5.0 absx) y)) absx)
                      (Vec2. (if (m/>= y absx)
                               (m/- (m/* 3.0 absy) x)
                               (m/+ (m/* 7.0 absy) x)) absy))
           r (m/* amount (.y ps))
           a (-> m/M_PI_4
                 (m/* (.x ps))
                 (m// (.y ps))
                 (m/- m/M_PI_4))
           sa (m/sin a)
           ca (m/cos a)]
       (Vec2. (m/* r ca) (m/* r sa))))))

(defn blur
  "Blur"
  ([] {:type :random})
  ([^double amount _]
   (fn [_]
     (let [r (r/drand m/TWO_PI)
           sr (m/sin r)
           cr (m/cos r)
           r2 (r/drand amount)]
       (Vec2. (m/* r2 cr) (m/* r2 sr))))))

(defn blurlinear
  "Blur Linear"
  ([] {:type :random
       :config (fn [] {:length (u/sdrand 0.1 2.0)
                      :angle (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double length ^double angle]}]
   (let [s (m/sin angle)
         c (m/cos angle)]
     (fn [^Vec2 v]
       (let [r (m/* length (r/drand))]
         (Vec2. (m/* amount (m/+ (.x v) (m/* c r)))
                (m/* amount (m/+ (.y v) (m/* s r)))))))))

(defn blurpixelize
  "Blur Pixelize from Apo7X15C"
  ([] {:type :random
       :config (fn [] {:size (u/sdrand 0.01 1.2)
                      :scale (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double size ^double scale]}]
   (let [inv-size (m// 1.0 size)
         av (m/* amount size)
         half (Vec2. 0.5 0.5)]
     (fn [v]
       (-> v
           (v/mult inv-size)
           (v/fmap m/floor)
           (v/add (-> (v/generate-vec2 r/drand)
                      (v/sub half)
                      (v/mult scale)))
           (v/add half)
           (v/mult av))))))

(defn blurzoom
  "Blur Zoom from Apo7X15C"
  ([] {:type :random
       :config (fn [] {:length (r/drand -1.2 1.2)
                      :x (r/drand -1.2 1.2)
                      :y (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double length ^double x ^double y]}]
   (let [xy (Vec2. x y)
         xy- (Vec2. x (m/- y))]
     (fn [v]
       (-> v
           (v/sub xy)
           (v/mult (m/inc (r/drand length)))
           (v/add xy-)
           (v/mult amount))))))

(defn boarders
  "Boarders"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [roundx (m/rint (.x v))
           roundy (m/rint (.y v))
           offsetx (m/- (.x v) roundx)
           offsety (m/- (.y v) roundy)
           hoffsetx (m/* 0.5 offsetx)
           hoffsety (m/* 0.5 offsety)]
       (r/randval 0.75
                  (Vec2. (m/* amount (m/+ roundx hoffsetx))
                         (m/* amount (m/+ roundy hoffsety)))
                  (if (m/>= (m/abs offsetx) (m/abs offsety))
                    
                    (if (m/>= offsetx 0.0)
                      (Vec2. (m/* amount (m/+ hoffsetx roundx 0.25))
                             (m/* amount (m/+ hoffsety roundy (m// (m/* 0.25 offsety) offsetx))))
                      (Vec2. (m/* amount (m/- (m/+ hoffsetx roundx) 0.25))
                             (m/* amount (m/- (m/+ hoffsety roundy) (m// (m/* 0.25 offsety) offsetx)))))
                    
                    (if (m/>= offsety 0.0)
                      (Vec2. (m/* amount (m/+ hoffsetx roundx (m// (m/* 0.25 offsetx) offsety)))
                             (m/* amount (m/+ hoffsety roundy 0.25)))
                      (Vec2. (m/* amount (m/- (m/+ hoffsetx roundx) (m// (m/* 0.25 offsetx) offsety)))
                             (m/* amount (m/- (m/+ hoffsety roundy) 0.25))))))))))

(defn boarders2
  "Boarders"
  ([] {:type :random
       :config (fn [] {:c (r/drand -1.2 1.2)
                      :left (r/drand -1.2 1.2)
                      :right (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double c ^double left ^double right]}]
   (let [cc (m/abs c)
         cl (m/abs left)
         cr (m/abs right)
         cc (if (m/zero? cc) m/EPSILON cc)
         cl (if (m/zero? cl) m/EPSILON cl)
         cr (if (m/zero? cr) m/EPSILON cr)
         cl (m/* cc cl)
         cr (m/+ cc (m/* cc cr))]
     (fn [^Vec2 v]
       (let [roundx (m/rint (.x v))
             roundy (m/rint (.y v))
             offsetx (m/- (.x v) roundx)
             offsety (m/- (.y v) roundy)
             coffsetx (m/* c offsetx)
             coffsety (m/* c offsety)]
         (r/randval cr
                    (Vec2. (m/* amount (m/+ roundx coffsetx))
                           (m/* amount (m/+ roundy coffsety)))
                    (if (m/>= (m/abs offsetx) (m/abs offsety))
                      
                      (if (m/>= offsetx 0.0)
                        (Vec2. (m/* amount (m/+ coffsetx roundx cl))
                               (m/* amount (m/+ coffsety roundy (m// (m/* cl offsety) offsetx))))
                        (Vec2. (m/* amount (m/- (m/+ coffsetx roundx) cl))
                               (m/* amount (m/- (m/+ coffsety roundy) (m// (m/* cl offsety) offsetx)))))
                      
                      (if (m/>= offsety 0.0)
                        (Vec2. (m/* amount (m/+ coffsetx roundx (m// (m/* cl offsetx) offsety)))
                               (m/* amount (m/+ coffsety roundy cl)))
                        (Vec2. (m/* amount (m/- (m/+ coffsetx roundx) (m// (m/* cl offsetx) offsety)))
                               (m/* amount (m/- (m/+ coffsety roundy) cl)))))))))))

(defn- xypoints
  ([rng ^double var ^long n] (xypoints rng [u/zerov u/zerov] (m// var m/SQRT2) n))
  ([rng [v0 v1 :as v] ^double var ^long n]
   (if (m/zero? n)
     v
     (let [^Vec2 mid (v/add (v/mult (v/add v0 v1) 0.5)
                            (v/generate-vec2 #(r/grandom rng var)))
           nvar (m// var 2.7)
           nn (m/long-dec n)]
       (mapcat identity [(xypoints rng [v0 mid] nvar nn)
                         (xypoints rng [mid v1] nvar nn)])))))

(defn brownian
  "Brownian"
  ([] {:type :pattern
       :config (fn [] {:level (r/irand 2 20)
                      :variation (r/drand 3.0 10.0)
                      :seed (r/irand)
                      :line-thickness (r/drand 1.0 100.0)
                      :point-thickness (r/drand 1.0 100.0)
                      :show-lines (r/drand)
                      :show-points (r/drand)})})
  ([^double amount {:keys [^long level ^double variation ^long seed
                           ^double line-thickness ^double show-lines
                           ^double point-thickness ^double show-points]}]
   (let [line-thickness (m// line-thickness 100.0)
         point-thickness (m// point-thickness 100.0)
         show-sum (m/+ show-lines show-points)
         line-fraction (m// show-lines show-sum)
         line-threshold line-fraction
         rng (r/rng :mersenne seed)
         level (m/long-min 15 level)
         xys (atom (partition 2 1 (cycle (xypoints rng variation level))))]
     (fn [_]
       (let [[^Vec2 p1 ^Vec2 p2] (first @xys)
             out (if (m/< (r/drand) line-threshold)
                   (let [^Vec2 diff (v/sub p2 p1)                         
                         m (if (m/zero? (.x diff)) 10000.0 (m// (.y diff) (.x diff)))
                         line-length (v/mag diff)
                         d (r/drand line-length)
                         xoffset (m// d (m/sqrt (m/inc (m/* m m))))
                         xoffset (if (m/< (.x p2) (.x p1)) (m/- xoffset) xoffset)
                         yoffset (m/abs (m/* m xoffset))
                         yoffset (if (m/< (.y p2) (.y p1)) (m/- yoffset) yoffset)]
                     (-> (v/add p1 (Vec2. xoffset yoffset))
                         (v/add (v/generate-vec2 #(m/* line-thickness (r/drand -0.5 0.5))))))
                   (let [roffset (r/drand point-thickness)
                         rangle (r/drand m/TWO_PI)]
                     (v/add p1 (Vec2. (m/* roffset (m/cos rangle))
                                      (m/* roffset (m/sin rangle))))))]
         (swap! xys rest)
         (v/mult out amount))))))

(defn bubble
  "Bubble"
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (v/mult v (m// amount (m/inc (m/* 0.25 (v/mag v))))))))

(defn bulge
  "Bulge"
  ([] {:type :regular
       :config (fn [] {:N (r/drand -4.0 4.0)})})
  ([^double amount {:keys [^double N]}]
   (fn [v]
     (let [r (v/mag v)
           rn (m/pow r N)]
       (v/mult v (m// (m/* amount rn) r))))))

(defn butterflyfay
  "Butterfly Fay"
  ([] {:type :random
       :config (fn [] {:offset (r/randval 0.0 (r/drand -5.0 5.0))
                      :unified-inner-outer (r/brand)
                      :outer-mode (r/irand 6)
                      :inner-mode (r/irand 6)
                      :outer-spread (r/drand -1.1 1.1)
                      :inner-spread (r/drand -1.1 1.1)
                      :outer-spread-ratio (u/sdrand 0.01 1.1)
                      :inner-spread-ratio (u/sdrand 0.01 1.1)
                      :spread-split (u/sdrand 0.5 1.5)
                      :fill (r/randval 0.0 (r/drand -2.0 2.0))})})
  ([^double amount {:keys [^double spread-split ^double offset ^double fill
                           unified-inner-outer
                           ^long outer-mode ^double outer-spread ^double outer-spread-ratio
                           ^long inner-mode ^double inner-spread ^double inner-spread-ratio]}]
   (let [oo (m/* outer-spread outer-spread-ratio)
         ii (m/* inner-spread inner-spread-ratio)
         outer-mode (int outer-mode)
         inner-mode (int inner-mode)]
     (fn [^Vec2 v]
       (let [theta (v/heading v)
             t (m/* m/PI m/PI theta)
             rin (m/* spread-split (v/mag v))
             r (m/* 0.5 (m/+ offset
                         (m/- (m/exp (m/cos t))
                            (m/* 2.0 (m/cos (m/* 4.0 t)))
                            (m/pow (m/abs (m/sin (m// t 12.0))) 5.0))))
             r (if (m/zero? fill) r (m/+ r (m/* fill (r/drand -0.5) 0.5)))
             x (m/* r (m/sin t))
             y (m/- (m/* r (m/cos t)))
             res (if (or unified-inner-outer
                         (m/> (m/abs rin) (m/abs r)))
                   (case outer-mode
                     1 (let [rinx (m/inc (m/- (m/* rin oo) oo))
                             riny (m/inc (m/- (m/* rin outer-spread) outer-spread))]
                         (Vec2. (m/* x rinx) (m/* y riny)))
                     2 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (m/+ x (m/* oo (m/- xin x)))
                                (m/+ y (m/* outer-spread (m/- yin y)))))
                     3 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (m/+ x (m/* oo xin))
                                (m/+ y (m/* outer-spread yin))))
                     4 (let [hrin (m/* 0.5 rin)
                             rinx (m/+ hrin oo)
                             riny (m/+ hrin outer-spread)]
                         (Vec2. (m/* x rinx) (m/* y riny)))
                     5 (Vec2. (m/+ x (m/* oo (.x v)))
                              (m/+ y (m/* outer-spread (.y v))))
                     (Vec2. x y))
                   (case inner-mode
                     1 (let [rinx (m/inc (m/- (m/* rin ii) ii))
                             riny (m/inc (m/- (m/* rin inner-spread) inner-spread))]
                         (Vec2. (m/* x rinx) (m/* y riny)))
                     2 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (m/+ x (m/* ii (m/- xin x)))
                                (m/+ y (m/* inner-spread (m/- yin y)))))
                     3 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (m/+ x (m/* ii xin))
                                (m/+ y (m/* inner-spread yin))))
                     4 (let [hrin (m/* 0.5 rin)
                             rinx (m/+ hrin ii)
                             riny (m/+ hrin inner-spread)]
                         (Vec2. (m/* x rinx) (m/* y riny)))
                     5 (Vec2. (m/+ x (m/* ii (.x v)))
                              (m/+ y (m/* inner-spread (.y v))))
                     (Vec2. x y)))]
         (v/mult res amount))))))

(defn butterfly
  "Butterfly"
  ([] {:type :regular})
  ([^double amount _]
   (let [wx (m/* amount 1.3029400317411197908970256609023)]
     (fn [^Vec2 v]
       (let [y2 (m/* 2.0 (.y v))
             r (m/* wx (m/sqrt (m// (m/abs (m/* (.y v) (.x v)))
                                (m/+ m/EPSILON (m/sq (.x v)) (m/sq y2)))))]
         (Vec2. (m/* r (.x v))
                (m/* r y2)))))))

;; 

(defn besselj
  "Bessel"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (v/mag v) (special/bessel-J (m/abs (.x v)) (m/min 10000.0 (m/abs (.y v)))))
            (m/* amount (v/heading v))))))

(defn beta
  "Beta"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (special/log-beta (m/+ m/EPSILON (m/abs (.x v))) (m/+ m/EPSILON (m/abs (.y v)))))
            (m/* amount (v/heading v))))))

