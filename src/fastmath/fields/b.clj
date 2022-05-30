(ns fastmath.fields.b
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- bseries-calc
  "Common calculations for bSeries "
  [^double amount ^double tau ^double sigma]
  (let [sinht (m/sinh tau)
        cosht (m/cosh tau)
        sins (m/sin sigma)
        coss (m/cos sigma)
        temp (- cosht coss)]
    (if (zero? temp)
      (Vec2. 0.0 0.0)
      (Vec2. (* amount (/ sinht temp))
             (* amount (/ sins temp))))))

(defn bcollide
  "bCollide by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:num (u/sdrand 1.0 30.0)
                      :a (r/drand 2.0)})})
  ([^double amount {:keys [^double num ^double a]}]
   (let [bcn-pi (* num m/M_1_PI)
         pi-bcn (/ m/PI num)
         bca-bcn (/ (* m/PI a) num)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (- 1.0 (.x v)) (.y v))
             tau (* 0.5 (- (m/log (v/magsq v+))
                           (m/log (v/magsq v-))))
             pre-sigma (- m/PI (v/heading v+) (v/heading v-))
             alt (double (int (* pre-sigma bcn-pi)))
             sigma (if (even? (int alt))
                     (+ (* alt pi-bcn) (rem (+ pre-sigma bca-bcn) pi-bcn))
                     (+ (* alt pi-bcn) (rem (- pre-sigma bca-bcn) pi-bcn)))]
         (bseries-calc amount tau sigma))))))

(defn bmod
  "bMod by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:radius (r/drand 0.5 2.0)
                      :distance (r/drand 2.0)})})
  ([^double amount {:keys [^double radius ^double distance]}]
   (let [rd (* radius distance)
         r2 (+ radius radius)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (- 1.0 (.x v)) (.y v))
             pre-tau (* 0.5 (- (m/log (v/magsq v+))
                               (m/log (v/magsq v-))))
             sigma (- m/PI (v/heading v+) (v/heading v-))
             tau (if (and (< pre-tau radius) (< (- pre-tau) radius))
                   (- (rem (+ pre-tau radius rd) r2) radius)
                   pre-tau)]
         (bseries-calc amount tau sigma))))))

(defn bsplit
  "Raykoid666, transcribed and modded by Nic Anderson, chronologicaldot"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -2.0 2.0)
                      :y (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [xx (+ x (.x v))]      
       (Vec2. (* (/ amount (m/tan xx)) (m/cos (+ y (.y v))))
              (* (/ amount (m/sin xx)) (- y (.y v))))))))

(defn bswirl
  "bSwirl by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :regular
       :config (fn [] {:in (r/drand -2.0 2.0)
                      :out (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double in ^double out]}]
   (fn [^Vec2 v]
     (let [v+ (v/add v u/unitx)
           v- (Vec2. (- 1.0 (.x v)) (.y v))
           tau (* 0.5 (- (m/log (v/magsq v+))
                         (m/log (v/magsq v-))))
           pre-sigma (- m/PI (v/heading v+) (v/heading v-))
           sigma (+ pre-sigma (* tau out) (/ in tau))]
       (bseries-calc amount tau sigma)))))

(defn btransform
  "bTransform by Michael Faber, http://michaelfaber.deviantart.com/art/bSeries-320574477"
  ([] {:type :random
       :config (fn [] {:rotate (r/drand m/TWO_PI)
                      :power (r/randval (r/drand 10.0) (r/irand 10.0))
                      :move (r/drand -2.0 2.0)
                      :split (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double rotate ^double power ^double move ^double split]}]
   (let [mp (/ m/TWO_PI power)]
     (fn [^Vec2 v]
       (let [v+ (v/add v u/unitx)
             v- (Vec2. (- 1.0 (.x v)) (.y v))
             pre-tau (+ (/ (* 0.5 (- (m/log (v/magsq v+))
                                     (m/log (v/magsq v-)))) power) move)
             pre-sigma (+ (- m/PI (v/heading v+) (v/heading v-)) rotate)
             tau (if (neg? (.x v))
                   (- pre-tau split)
                   (+ pre-tau split))
             sigma (+ (/ pre-sigma power)
                      (* mp (m/floor (* power (r/drand)))))]
         (bseries-calc amount tau sigma))))))

(defn- bwrands-bytemix ^long [^long a ^long b] (bit-xor (bit-and a 0x5a5a) (bit-and b 0xa5a5)))
(defn- bwrands-bytexim ^long [^long a ^long b] (bit-xor (bit-and a 0xaaaa) (bit-and b 0x5555)))
(defn- bwrands-byteshf ^long [^long a ^long b] (bit-and (bit-xor (<< a 8) (>> b 8)) 0xffff))
(defn- bwrands-byteprimes ^long [^long a ^long b]
  (bit-and 0xffff (bit-xor (- (* a 857) 4) (+ (* b -977) 8))))

(defn bwrands
  "bwrands"
  ([] {:type :random
       :config (fn [] (let [min-petals (r/irand 1 5)]
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
                        :max-petals (+ min-petals (r/irand 15))}))})
  ([^double amount {:keys [^double cellsize ^double space ^double gain ^double inner-twist ^double outer-twist
                           ^long seed ^double rrot ^double rmin ^double loonie-chance ^double petal-chance
                           ^long min-petals ^long max-petals]}]
   (let [radius (* 0.5 (/ cellsize (inc (m/sq space))))
         g2 (+ m/EPSILON (/ (m/sq gain) radius))
         max-bubble (as-> (* g2 radius) a
                      (if (> a 2.0) 1.0 (* a (/ (inc (* 0.25 a a))))))
         r2 (* radius radius)
         rfactor (/ radius max-bubble)
         petx min-petals
         pety (- max-petals min-petals)]
     (fn [^Vec2 v]
       (if (< (m/abs cellsize) m/EPSILON)
         (v/mult v amount)
         (let [^Vec2 I (v/floor (v/div v cellsize))
               C (v/mult (v/shift I 0.5) cellsize)
               xx (bit-xor (int (.x I)) 0xb641)
               yy (bit-xor (int (.y I)) 0x9d81)
               xy (bit-and 0xffff (+ seed (* xx yy)))
               xx (bit-and 0xffff xx)
               yy (bit-and 0xffff yy)
               tt (bwrands-bytemix xx yy)
               yy (bwrands-bytemix yy xx)
               xx tt
               tt (bwrands-byteshf xx yy)
               yy (bwrands-byteshf xy yy)
               xx tt
               tt (bwrands-bytexim xx yy)
               yy (bwrands-bytexim yy xx)
               xx tt
               ssz (-> (/ xx 65536.0)
                       (* (- 1.0 rmin))
                       (+ rmin))
               aan (/ (* rrot m/TWO_PI yy) 65536.0)
               tt (bwrands-byteprimes xx yy)
               yy (bwrands-byteprimes yy xx)
               xx tt
               LC (+ (/ xx -65536.0) loonie-chance)
               PC (if (neg? LC) (+ LC petal-chance) 0.0)
               LC (if (<= LC PC) -1.0 LC)
               L (v/sub v C)
               vv2 (* ssz r2)
               mL (v/magsq L)]
           (if (> mL vv2)
             (v/mult v amount)
             (let [nL (if (pos? LC)
                        (v/mult L (m/sqrt (dec (/ vv2 mL))))
                        (if (pos? PC)
                          (let [NPetals (if (zero? pety)
                                          petx
                                          (+ petx (mod (>> yy 3) (inc pety))))
                                flrw (as-> (* (+ m/PI (v/heading L))
                                              (/ NPetals m/TWO_PI)) a
                                       (- a (int a))
                                       (* 2.0 (m/abs (- a 0.5))))
                                r (m/sqrt mL)]
                            (r/randval (* 0.5 (+ flrw 0.5))
                                       (v/mult L (* (- 1.0 r) (* 1.1 flrw)))
                                       (v/mult L (/ (- (m/sqrt vv2)
                                                       (* r (- 1.0 flrw))) (+ r m/EPSILON)))))
                          (let [nL (v/mult L g2)
                                r (/ rfactor (inc (/ (v/magsq nL) (* 4.0 ssz))))]
                            (v/mult nL r))))
                   nvvl (if-not (or (pos? LC) (pos? PC)) (m/sqrt ssz) 1.0)
                   r (/ (* nvvl (v/mag nL)) r2)
                   theta (m/lerp inner-twist outer-twist r)]
               (v/mult (v/add C (v/rotate nL (+ aan theta))) amount)))))))))

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
   (let [radius (* 0.5 (/ cellsize (inc (m/sq space))))
         g2 (+ m/EPSILON (* gain gain))
         max-bubble- (* g2 radius)
         max-bubble (if (> max-bubble- 2.0)
                      1.0
                      (* max-bubble- (/ 1.0 (inc (* 0.25 (m/sq max-bubble-))))))
         r2 (m/sq radius)
         rfactor (/ radius max-bubble)]
     (fn [^Vec2 v]
       (if (< (m/abs cellsize) m/EPSILON)
         (v/mult v amount)
         (let [^Vec2 C (-> v
                           (v/div cellsize)
                           (v/fmap m/floor)
                           (v/add (Vec2. 0.5 0.5))
                           (v/mult cellsize))
               L (v/sub v C)]
           (if (> (v/magsq L) r2)
             (v/mult v amount)
             (let [L (v/mult L g2)
                   r (/ rfactor (inc (* 0.25 (v/magsq L))))
                   ^Vec2 L (v/mult L r)
                   r (/ (v/magsq L) r2)
                   theta (+ (* inner-twist (- 1.0 r))
                            (* outer-twist r))
                   s (m/sin theta)
                   c (m/cos theta)
                   vx (+ (.x C) (* c (.x L)) (* s (.y L)))
                   vy (+ (.y C) (* -1.0 s (.x L)) (* c (.y L)))]
               (v/mult (Vec2. vx vy) amount)))))))))

(defn barycentroid
  "barycentroid from Xyrus02, http://xyrusworx.deviantart.com/art/Barycentroid-Plugin-144832371?q=sort%3Atime+favby%3Amistywisp&qo=0&offset=10"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (let [dot00 (+ (* a a) (* b b))
         dot01 (+ (* a c) (* b d))
         dot11 (+ (* c c) (* d d))
         inv-denom (/ 1.0 (- (* dot00 dot11) (* dot01 dot01)))]
     (fn [^Vec2 v]
       (let [dot02 (+ (* a (.x v)) (* b (.y v)))
             dot12 (+ (* c (.x v)) (* d (.y v)))
             u (* inv-denom (- (* dot11 dot02) (* dot01 dot12)))
             vv (* inv-denom (- (* dot00 dot12) (* dot01 dot02)))
             um (* (m/signum u) (m/sqrt (+ (* u u) (m/sq (.x v)))))
             vm (* (m/signum vv) (m/sqrt (+ (* vv vv) (m/sq (.y v)))))]
         (Vec2. (* amount um)
                (* amount vm)))))))

(defn bent
  "Bent"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [nx (if (neg? (.x v)) (+ (.x v) (.x v)) (.x v))
           ny (if (neg? (.y v)) (* (.y v) 0.5) (.y v))]
       (Vec2. (* amount nx)
              (* amount ny))))))

(defn bent2
  "Bent2"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.5 2.0)
                      :y (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (let [nx (if (neg? (.x v)) (* (.x v) x) (.x v))
           ny (if (neg? (.y v)) (* (.y v) y) (.y v))]
       (Vec2. (* amount nx)
              (* amount ny))))))

(defn bilinear
  "Bilinear"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (.y v))
            (* amount (.x v))))))

(defn bipolar
  "Bipolar"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shift]}]
   (let [ps (* (- m/HALF_PI) shift)]
     (fn [^Vec2 v]
       (let [x2y2 (v/magsq v)
             t (inc x2y2)
             x2 (+ (.x v) (.x v))
             pre-y (+ ps (* 0.5 (m/atan2 (+ (.y v) (.y v))
                                         (dec x2y2))))
             y (if (> pre-y m/HALF_PI)
                 (- (rem (+ pre-y m/HALF_PI) m/PI) m/HALF_PI)
                 (if (< pre-y (- m/HALF_PI))
                   (- m/HALF_PI (rem (- m/HALF_PI pre-y) m/PI))
                   pre-y))
             f (+ t x2)
             g (- t x2)]
         (if (or (zero? g)
                 (not (pos? (/ f g))))
           (Vec2. 0.0 0.0)
           (Vec2. (* amount m/M_2_PI 0.25 (m/log (/ f g)))
                  (* amount m/M_2_PI y))))))))

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
   (let [ps (* (- m/HALF_PI) shift)]
     (fn [^Vec2 v]
       (let [x2y2 (* g1 (v/magsq v))
             t (+ a x2y2)
             x2 (* b (.x v))
             pre-y (+ ps (* c (m/atan2 (* e (.y v))
                                       (- x2y2 d))))
             y (if (> pre-y m/HALF_PI)
                 (- (rem (+ pre-y m/HALF_PI) m/PI) m/HALF_PI)
                 (if (< pre-y (- m/HALF_PI))
                   (- m/HALF_PI (rem (- m/HALF_PI pre-y) m/PI))
                   pre-y))
             f (+ t x2)
             g (- t x2)]
         (if (or (zero? g)
                 (not (pos? (/ f g))))
           (Vec2. 0.0 0.0)
           (Vec2. (* amount f1 m/M_2_PI 0.25 (m/log (/ f g)))
                  (* amount m/M_2_PI y h))))))))

(defn blade
  "Blade"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (* (r/drand amount) (v/mag v))
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (* amount (.x v) (+ cosr sinr))
              (* amount (.x v) (- cosr sinr)))))))

(defn blade2
  "Blade2"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (* (r/drand amount) (v/mag v))
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (* amount (.x v) (+ cosr sinr))
              (* amount (.y v) (- cosr sinr)))))))

(defn blob
  "Blob"
  ([] {:type :regular
       :config (fn [] {:low (u/sdrand 0.1 2.0)
                      :high (u/sdrand 0.1 2.0)
                      :waves (r/randval (r/irand -6 6) (r/drand -6.0 6.0))})})
  ([^double amount {:keys [^double low ^double high ^double waves]}]
   (let [hl (- high low)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             rr (->> (* a waves)
                     (m/sin)
                     (* 0.5)
                     (+ 0.5)
                     (* hl)
                     (+ low)
                     (* r))]
         (Vec2. (* amount rr (m/sin a))
                (* amount rr (m/cos a))))))))

(defn blocky
  "blocky from FracFx, http://fracfx.deviantart.com/art/FracFx-Plugin-Pack-171806681"
  ([] {:type :regular
       :config (fn [] {:x (u/sdrand 0.5 1.5)
                      :y (u/sdrand 0.5 1.5)
                      :mp (u/sdrand 0.001 6.0)})})
  ([^double amount {:keys [^double x ^double y ^double mp]}]
   (let [vv (/ amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [T (inc (/ (+ (m/cos (.x v)) (m/cos (.y v))) mp))
             r (/ amount T)
             tmp (inc (v/magsq v))
             x2 (+ (.x v) (.x v))
             y2 (+ (.y v) (.y v))
             xmax (* 0.5 (+ (m/sqrt (+ tmp x2)) (m/sqrt (- tmp x2))))
             ymax (* 0.5 (+ (m/sqrt (+ tmp y2)) (m/sqrt (- tmp y2))))
             ax (/ (.x v) xmax)
             bx (m/safe-sqrt (- 1.0 (m/sq ax)))
             ay (/ (.y v) ymax)
             by (m/safe-sqrt (- 1.0 (m/sq ay)))]
         (Vec2. (* vv (m/atan2 ax bx) r x)
                (* vv (m/atan2 ay by) r y)))))))

(defn blurcircle
  "Blur circle"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [x (r/drand -1.0 1.0)
           y (r/drand -1.0 1.0)
           absx (m/abs x)
           absy (m/abs y)
           ^Vec2 ps (if (>= absx absy)
                      (Vec2. (if (>= x absy)
                               (+ absx y)
                               (- (* 5.0 absx) y)) absx)
                      (Vec2. (if (>= y absx)
                               (- (* 3.0 absy) x)
                               (+ (* 7.0 absy) x)) absy))
           r (* amount (.y ps))
           a (-> m/M_PI_4
                 (* (.x ps))
                 (/ (.y ps))
                 (- m/M_PI_4))
           sa (m/sin a)
           ca (m/cos a)]
       (Vec2. (* r ca) (* r sa))))))

(defn blur
  "Blur"
  ([] {:type :random})
  ([^double amount _]
   (fn [_]
     (let [r (r/drand m/TWO_PI)
           sr (m/sin r)
           cr (m/cos r)
           r2 (r/drand amount)]
       (Vec2. (* r2 cr) (* r2 sr))))))

(defn blurlinear
  "Blur Linear"
  ([] {:type :random
       :config (fn [] {:length (u/sdrand 0.1 2.0)
                      :angle (r/drand m/-PI m/PI)})})
  ([^double amount {:keys [^double length ^double angle]}]
   (let [s (m/sin angle)
         c (m/cos angle)]
     (fn [^Vec2 v]
       (let [r (* length (r/drand))]
         (Vec2. (* amount (+ (.x v) (* c r)))
                (* amount (+ (.y v) (* s r)))))))))

(defn blurpixelize
  "Blur Pixelize from Apo7X15C"
  ([] {:type :random
       :config (fn [] {:size (u/sdrand 0.01 1.2)
                      :scale (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double size ^double scale]}]
   (let [inv-size (/ 1.0 size)
         av (* amount size)
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
         xy- (Vec2. x (- y))]
     (fn [v]
       (-> v
           (v/sub xy)
           (v/mult (inc (r/drand length)))
           (v/add xy-)
           (v/mult amount))))))

(defn boarders
  "Boarders"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [roundx (m/rint (.x v))
           roundy (m/rint (.y v))
           offsetx (- (.x v) roundx)
           offsety (- (.y v) roundy)
           hoffsetx (* 0.5 offsetx)
           hoffsety (* 0.5 offsety)]
       (r/randval 0.75
                  (Vec2. (* amount (+ roundx hoffsetx))
                         (* amount (+ roundy hoffsety)))
                  (if (>= (m/abs offsetx) (m/abs offsety))
                    
                    (if (>= offsetx 0.0)
                      (Vec2. (* amount (+ hoffsetx roundx 0.25))
                             (* amount (+ hoffsety roundy (/ (* 0.25 offsety) offsetx))))
                      (Vec2. (* amount (- (+ hoffsetx roundx) 0.25))
                             (* amount (- (+ hoffsety roundy) (/ (* 0.25 offsety) offsetx)))))
                    
                    (if (>= offsety 0.0)
                      (Vec2. (* amount (+ hoffsetx roundx (/ (* 0.25 offsetx) offsety)))
                             (* amount (+ hoffsety roundy 0.25)))
                      (Vec2. (* amount (- (+ hoffsetx roundx) (/ (* 0.25 offsetx) offsety)))
                             (* amount (- (+ hoffsety roundy) 0.25))))))))))

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
         cc (if (zero? cc) m/EPSILON cc)
         cl (if (zero? cl) m/EPSILON cl)
         cr (if (zero? cr) m/EPSILON cr)
         cl (* cc cl)
         cr (+ cc (* cc cr))]
     (fn [^Vec2 v]
       (let [roundx (m/rint (.x v))
             roundy (m/rint (.y v))
             offsetx (- (.x v) roundx)
             offsety (- (.y v) roundy)
             coffsetx (* c offsetx)
             coffsety (* c offsety)]
         (r/randval cr
                    (Vec2. (* amount (+ roundx coffsetx))
                           (* amount (+ roundy coffsety)))
                    (if (>= (m/abs offsetx) (m/abs offsety))
                      
                      (if (>= offsetx 0.0)
                        (Vec2. (* amount (+ coffsetx roundx cl))
                               (* amount (+ coffsety roundy (/ (* cl offsety) offsetx))))
                        (Vec2. (* amount (- (+ coffsetx roundx) cl))
                               (* amount (- (+ coffsety roundy) (/ (* cl offsety) offsetx)))))
                      
                      (if (>= offsety 0.0)
                        (Vec2. (* amount (+ coffsetx roundx (/ (* cl offsetx) offsety)))
                               (* amount (+ coffsety roundy cl)))
                        (Vec2. (* amount (- (+ coffsetx roundx) (/ (* cl offsetx) offsety)))
                               (* amount (- (+ coffsety roundy) cl)))))))))))

(defn- xypoints
  ([rng ^double var ^long n] (xypoints rng [u/zerov u/zerov] (/ var m/SQRT2) n))
  ([rng [v0 v1 :as v] ^double var ^long n]
   (if (zero? n)
     v
     (let [^Vec2 mid (v/add (v/mult (v/add v0 v1) 0.5)
                            (v/generate-vec2 #(r/grandom rng var)))
           nvar (/ var 2.7)
           nn (dec n)]
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
   (let [line-thickness (/ line-thickness 100.0)
         point-thickness (/ point-thickness 100.0)
         show-sum (+ show-lines show-points)
         line-fraction (/ show-lines show-sum)
         line-threshold line-fraction
         rng (r/rng :mersenne seed)
         level (min 15 level)
         xys (atom (cycle (xypoints rng variation level)))]
     (fn [_]
       (let [[^Vec2 p1 ^Vec2 p2] (first @xys)
             out (if (< (r/drand) line-threshold)
                   (let [^Vec2 diff (v/sub p2 p1)
                         m (if (zero? (.x diff)) 10000.0 (/ (.y diff) (.x diff)))
                         line-length (v/mag diff)
                         d (r/drand line-length)
                         xoffset (/ d (m/sqrt (inc (* m m))))
                         xoffset (if (< (.x p2) (.x p1)) (- xoffset) xoffset)
                         yoffset (m/abs (* m xoffset))
                         yoffset (if (< (.y p2) (.y p1)) (- yoffset) yoffset)]
                     (-> (v/add p1 (Vec2. xoffset yoffset))
                         (v/add (v/generate-vec2 #(* line-thickness (r/drand -0.5 0.5))))))
                   (let [roffset (r/drand point-thickness)
                         rangle (r/drand m/TWO_PI)]
                     (v/add p1 (Vec2. (* roffset (m/cos rangle))
                                      (* roffset (m/sin rangle))))))]
         (swap! xys rest)
         (v/mult out amount))))))

(defn bubble
  "Bubble"
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (v/mult v (/ amount (inc (* 0.25 (v/mag v))))))))

(defn bulge
  "Bulge"
  ([] {:type :regular
       :config (fn [] {:N (r/drand -4.0 4.0)})})
  ([^double amount {:keys [^double N]}]
   (fn [v]
     (let [r (v/mag v)
           rn (m/pow r N)]
       (v/mult v (/ (* amount rn) r))))))

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
   (let [oo (* outer-spread outer-spread-ratio)
         ii (* inner-spread inner-spread-ratio)
         outer-mode (int outer-mode)
         inner-mode (int inner-mode)]
     (fn [^Vec2 v]
       (let [theta (v/heading v)
             t (* m/PI m/PI theta)
             rin (* spread-split (v/mag v))
             r (* 0.5 (+ offset
                         (- (m/exp (m/cos t))
                            (* 2.0 (m/cos (* 4.0 t)))
                            (m/pow (m/abs (m/sin (/ t 12.0))) 5.0))))
             r (if (zero? fill) r (+ r (* fill (r/drand -0.5) 0.5)))
             x (* r (m/sin t))
             y (- (* r (m/cos t)))
             res (if (or unified-inner-outer
                         (> (m/abs rin) (m/abs r)))
                   (case outer-mode
                     1 (let [rinx (inc (- (* rin oo) oo))
                             riny (inc (- (* rin outer-spread) outer-spread))]
                         (Vec2. (* x rinx) (* y riny)))
                     2 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (+ x (* oo (- xin x)))
                                (+ y (* outer-spread (- yin y)))))
                     3 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (+ x (* oo xin))
                                (+ y (* outer-spread yin))))
                     4 (let [hrin (* 0.5 rin)
                             rinx (+ hrin oo)
                             riny (+ hrin outer-spread)]
                         (Vec2. (* x rinx) (* y riny)))
                     5 (Vec2. (+ x (* oo (.x v)))
                              (+ y (* outer-spread (.y v))))
                     (Vec2. x y))
                   (case inner-mode
                     1 (let [rinx (inc (- (* rin ii) ii))
                             riny (inc (- (* rin inner-spread) inner-spread))]
                         (Vec2. (* x rinx) (* y riny)))
                     2 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (+ x (* ii (- xin x)))
                                (+ y (* inner-spread (- yin y)))))
                     3 (let [xin (m/copy-sign (.x v) x)
                             yin (m/copy-sign (.y v) y)]
                         (Vec2. (+ x (* ii xin))
                                (+ y (* inner-spread yin))))
                     4 (let [hrin (* 0.5 rin)
                             rinx (+ hrin ii)
                             riny (+ hrin inner-spread)]
                         (Vec2. (* x rinx) (* y riny)))
                     5 (Vec2. (+ x (* ii (.x v)))
                              (+ y (* inner-spread (.y v))))
                     (Vec2. x y)))]
         (v/mult res amount))))))

(defn butterfly
  "Butterfly"
  ([] {:type :regular})
  ([^double amount _]
   (let [wx (* amount 1.3029400317411197908970256609023)]
     (fn [^Vec2 v]
       (let [y2 (* 2.0 (.y v))
             r (* wx (m/sqrt (/ (m/abs (* (.y v) (.x v)))
                                (+ m/EPSILON (m/sq (.x v)) (m/sq y2)))))]
         (Vec2. (* r (.x v))
                (* r y2)))))))

;; 

(defn besselj
  "Bessel"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (v/mag v) (m/bessel-j (m/abs (.x v)) (m/abs (.y v))))
            (* amount (v/heading v))))))

(defn beta
  "Beta"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (m/log-beta (+ m/EPSILON (m/abs (.x v))) (+ m/EPSILON (m/abs (.y v)))))
            (* amount (v/heading v))))))
