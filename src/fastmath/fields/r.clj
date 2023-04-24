(ns fastmath.fields.r
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn r-circleblur
  "R Circle Blur"
  ([] {:type :random
       :config (fn [] {:n (u/sdrand 0.9 7.0)
                      :seed (r/drand Float/MAX_VALUE)
                      :dist (r/drand -1.0 1.0)
                      :mn (r/drand -3.0 3.0)
                      :mx (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double n ^double seed ^double dist ^double mn ^double mx]}]
   (let [dm (- mx mn)]
     (fn [v]
       (let [angle (v/heading v)
             rad (v/mag v)
             rad (mod rad n)
             by (m/sin (+ angle rad))
             bx (m/cos (+ angle rad))
             by (m/round (* by rad))
             bx (m/round (* bx rad))
             rad2 (* 0.5 (m/sqrt (r/drand)))
             angle2 (r/drand m/TWO_PI)
             a1 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 127.1) (* by 311.7) seed))))
             a2 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 269.5) (* by 183.3) seed))))
             a3 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 78.233) (* by 12.9898) seed))))
             a3 (+ mn (* a3 dm))
             rad2 (* rad2 a3)]
         (Vec2. (* amount (+ bx (* rad2 (m/cos angle2)) (* dist a1)))
                (* amount (+ by (* rad2 (m/sin angle2)) (* dist a2)))))))))

(defn radialblur
  "Radial blur"
  ([] {:type :random
       :config (fn [] {:angle (r/drand (- m/TWO_PI) m/TWO_PI)})})
  ([^double amount {:keys [^double angle]}]
   (let [spin (* amount (m/sin (* angle m/HALF_PI)))
         zoom (* amount (m/cos (* angle m/HALF_PI)))]
     (fn [^Vec2 v]
       (let [rnd-g (+ (r/drand) (r/drand) (r/drand) (r/drand) -2.0)
             ra (v/mag v)
             alpha (+ (* spin rnd-g) (v/heading v))
             rz (dec (* zoom rnd-g))]
         (Vec2. (+ (* rz (.x v)) (* ra (m/cos alpha)))
                (+ (* rz (.y v)) (* ra (m/sin alpha)))))))))

(defn rational3
  "Rational3"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -3.0 3.0)
                      :b (r/drand -3.0 3.0)
                      :c (r/drand -3.0 3.0)
                      :d (r/drand -3.0 3.0)
                      :e (r/drand -3.0 3.0)
                      :f (r/drand -3.0 3.0)
                      :g (r/drand -3.0 3.0)
                      :h (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d
                           ^double e ^double f ^double g ^double h]}]
   (fn [^Vec2 v]
     (let [^Vec2 sqr (v/emult v v)
           ^Vec2 cb (v/emult sqr v)
           zt3 (- (.x cb) (* 3.0 (.x v) (.y sqr)))
           zt2 (- (.x sqr) (.y sqr))
           zb3 (- (* 3.0 (.x sqr) (.y v)) (.y cb))
           zb2 (* 2.0 (.x v) (.y v))
           tr (+ (* a zt3) (* b zt2) (* c (.x v)) d)
           ti (+ (* a zb3) (* b zb2) (* c (.y v )))
           br (+ (* e zt3) (* f zt2) (* g (.x v)) h)
           bi (+ (* e zb3) (* f zb2) (* g (.y v )))
           r3den (/ amount (+ (* br br) (* bi bi)))]
       (Vec2. (* r3den (+ (* tr br) (* ti bi)))
              (* r3den (- (* ti br) (* tr bi))))))))

(defn rays1
  "Rays1"
  ([] {:type :regular})
  ([^double amount _]
   (let [pa (* amount (m/sq m/M_2_PI))]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             u (+ pa (/ (m/tan (m/sqrt t))))
             r (* amount u t)]
         (Vec2. (/ r (.x v))
                (/ r (.y v))))))))

(defn rays2
  "Rays2"
  ([] {:type :regular})
  ([^double amount _]
   (let [a10 (/ amount 10.0)]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             u (/ (m/cos (* (+ t m/EPSILON) (m/tan (/ (+ m/EPSILON t))))))
             r (* a10 t u)]
         (Vec2. (/ r (.x v))
                (/ r (.y v))))))))

(defn rays3
  "Rays3"
  ([] {:type :regular})
  ([^double amount _]
   (let [a10 (/ amount 10.0)]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             t2 (* t t)
             u (/ (m/sqrt (m/cos (m/sin (* (+ m/EPSILON t2) (m/sin (/ (+ t2 m/EPSILON))))))))
             r (* a10 t u)]
         (Vec2. (/ (* r (m/cos t)) (.x v))
                (/ (* r (m/tan t)) (.y v))))))))

(defn rays
  "Rays"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [ang (* amount (r/drand m/PI))
           r (/ amount (+ m/EPSILON (v/magsq v)))
           tanr (* amount r (m/tan ang))]
       (Vec2. (* tanr (m/cos (.x v)))
              (* tanr (m/sin (.y v))))))))

(defn rectangles
  "Rectangles"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (Vec2. (if (< (m/abs (.x v)) m/EPSILON)
              (* amount (.x v))
              (* amount (-> (.x v)
                            (/ x)
                            (m/floor)
                            (* 2.0)
                            inc
                            (* x)
                            (- (.x v)))))
            (if (< (m/abs (.y v)) m/EPSILON)
              (* amount (.y v))
              (* amount (-> (.y v)
                            (/ y)
                            (m/floor)
                            (* 2.0)
                            inc
                            (* y)
                            (- (.y v)))))))))

(defn rhodonea
  "Rhodonea"
  ([] {:type :random
       :config (fn [] {:knumer (r/randval (int (u/sdrand 1.0 30.0)) (u/sdrand 1.0 30.0))
                      :kdenom (r/randval (int (u/sdrand 1.0 30.0)) (u/sdrand 1.0 30.0))
                      :radial-offset (r/drand -1.0 1.0)
                      :inner-mode (r/irand 7)
                      :outer-mode (r/irand 7)
                      :inner-spread (r/drand -1.0 1.0)
                      :outer-spread (r/drand -1.0 1.0)
                      :inner-spread-ratio (r/drand -2.0 2.0)
                      :outer-spread-ratio (r/drand -2.0 2.0)
                      :spread-split (r/drand -1.5 1.5)
                      :cycle-offset (r/drand m/TWO_PI)
                      :cycles-param (r/randval 0.0 (r/drand 100.0))
                      :metacycle-expansion (r/drand -1.0 1.0)
                      :metacycles (r/drand 10.0)
                      :fill (r/randval 0.0 (r/drand))})})
  ([^double amount {:keys [^double knumer ^double kdenom ^double radial-offset
                           ^long inner-mode ^long outer-mode
                           ^double inner-spread ^double outer-spread
                           ^double inner-spread-ratio ^double outer-spread-ratio
                           ^double spread-split ^double cycles-param
                           ^double cycle-offset ^double metacycle-expansion
                           ^double metacycles ^double fill]}]
   (let [kn knumer
         kd kdenom
         k (/ kn kd)
         cycles-to-close (double (if (zero? (mod k 1.0))
                                   (if (zero? (mod k 2.0))
                                     1.0
                                     (if (or (not (zero? radial-offset))
                                             (not (zero? inner-spread))
                                             (not (zero? outer-spread))
                                             (not (zero? fill)))
                                       1.0
                                       0.5))
                                   (if (and (zero? (mod kn 1.0))
                                            (zero? (mod kd 1.0)))
                                     (let [lkn (long kn)
                                           lkd (long kd)
                                           gcd (m/gcd lkn lkd)
                                           [^long kn ^long kd] (if (m/not== gcd 1)
                                                                 [(/ lkn gcd) (/ lkd gcd)]
                                                                 [lkn lkd])]
                                       (if (or (zero? (mod kn 2.0))
                                               (zero? (mod kd 2.0)))
                                         kd
                                         (/ kd 2)))
                                     (if (< cycles-param 16)
                                       16
                                       (* 2 kd kn)))))
         cycles (if (zero? cycles-param)
                  (* cycles-to-close metacycles)
                  cycles-param)]
     (fn [^Vec2 v]
       (let [rin (* spread-split (v/mag v))
             tin (v/heading v)
             t (* cycles (+ tin (* cycle-offset m/TWO_PI)))
             r (+ radial-offset (m/cos (* t k)))
             r (if-not (zero? fill)
                 (+ r (* fill (- (r/drand) 0.5)))
                 r)
             x (* r (m/cos t))
             y (* r (m/sin t))
             expansion (m/floor (/ (* cycles (+ tin m/PI))
                                   (* cycles-to-close m/TWO_PI)))
             adjusted-amount (+ amount (* expansion metacycle-expansion))]
         (if (> (m/abs rin) (m/abs r))
           (case (long outer-mode)
             0 (Vec2. (* adjusted-amount x)
                      (* adjusted-amount y))
             1 (let [rinx (inc (* (dec rin) outer-spread outer-spread-ratio))
                     riny (inc (* (dec rin) outer-spread))]
                 (Vec2. (* adjusted-amount rinx x)
                        (* adjusted-amount riny y)))
             2 (let [xin (* (m/sgn x) (m/abs (.x v)))
                     yin (* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (* adjusted-amount (+ x (* outer-spread outer-spread-ratio (- xin x))))
                        (* adjusted-amount (+ y (* outer-spread (- yin y))))))
             3 (let [xin (* (m/sgn x) (m/abs (.x v)))
                     yin (* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (* adjusted-amount (+ x (* outer-spread outer-spread-ratio xin)))
                        (* adjusted-amount (+ y (* outer-spread yin)))))
             4 (let [rinx (+ (* 0.5 rin) (* outer-spread outer-spread-ratio))
                     riny (+ (* 0.5 rin) outer-spread)]
                 (Vec2. (* adjusted-amount rinx x)
                        (* adjusted-amount riny y)))
             5 v
             6 u/zerov)
           (case (long inner-mode)
             0 (Vec2. (* adjusted-amount x)
                      (* adjusted-amount y))
             1 (let [rinx (inc (* (dec rin) inner-spread inner-spread-ratio))
                     riny (inc (* (dec rin) inner-spread))]
                 (Vec2. (* adjusted-amount rinx x)
                        (* adjusted-amount riny y)))
             2 (let [xin (* (m/sgn x) (m/abs (.x v)))
                     yin (* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (* adjusted-amount (+ x (* inner-spread inner-spread-ratio (- xin x))))
                        (* adjusted-amount (+ y (* inner-spread (- yin y))))))
             3 (let [xin (* (m/sgn x) (m/abs (.x v)))
                     yin (* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (* adjusted-amount (+ x (* inner-spread inner-spread-ratio xin)))
                        (* adjusted-amount (+ y (* inner-spread yin)))))
             4 (let [rinx (+ (* 0.5 rin) (* inner-spread inner-spread-ratio))
                     riny (+ (* 0.5 rin) inner-spread)]
                 (Vec2. (* adjusted-amount rinx x)
                        (* adjusted-amount riny y)))
             5 v
             6 u/zerov)))))))

(defn rings2
  "Rings2"
  ([] {:type :regular
       :config (fn [] {:val (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double val]}]
   (let [dx (+ m/EPSILON (m/sq val))]
     (fn [^Vec2 v]
       (let [l (v/mag v)
             r (* amount (- 2.0 (* dx (inc (/ (* 2.0 (double (int (* 0.5 (inc (/ l dx)))))) l)))))]
         (v/mult v r))))))

(defn rings
  "Rings"
  ([] {:type :regular
       :config (fn [] {:coeff20 (r/drand 1.3)})})
  ([^double amount {:keys [^double coeff20]}]
   (let [dx (+ m/EPSILON (m/sq coeff20))
         dx2 (+ dx dx)
         rdx (/ 1.0 dx2)
         dx- (- 1.0 dx)]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             rr (* amount (+ (- r (* dx2 (double (int (* (+ r dx) rdx))))) (* r dx-)))]
         (Vec2. (* rr (/ (.x v) r))
                (* rr (/ (.y v) r))))))))

(defn ripple
  "Ripple"
  ([] {:type :regular
       :config (fn [] {:frequency (r/drand -3.0 3.0)
                      :velocity (r/drand -3.0 3.0)
                      :amplitude (r/drand -2.0 2.0)
                      :centerx (r/drand -0.5 0.5)
                      :centery (r/drand -0.5 0.5)
                      :phase (r/drand -1.0 1.0)
                      :scale (u/sdrand 0.5 4.0)
                      :fixed-dist-calc (r/irand 4)})})
  ([^double amount {:keys [^double frequency ^double velocity ^double amplitude ^double centerx
                           ^double centery ^double phase ^double scale ^int fixed-dist-calc]}]
   (let [f (* frequency 5.0)
         a (* amplitude 0.01)
         p (- (* phase m/TWO_PI) m/PI)
         s (if (zero? scale) m/EPSILON scale)
         is (/ s)
         vxp (* velocity p)
         pxa (* p a)
         pixa (* (- m/PI p) a)]
     (fn [^Vec2 v]
       (let [x (- (* s (.x v)) centerx)
             y (- (* s (.y v)) centery)
             d (case (unchecked-int fixed-dist-calc)
                 0 (v/mag v)
                 1 (m/sqrt (* (m/sq (.x v)) (m/sq (.y v))))
                 2 (max (m/abs (.x v)) (m/abs (.y v)))
                 3 (+ (m/abs (.x v)) (m/abs (.y v))))
             d (if (< d m/EPSILON) m/EPSILON d)
             nx (/ x d)
             ny (/ y d)
             wave (m/cos (- (* f d) vxp))
             d1 (+ d (* wave pxa))
             d2 (+ d (* wave pixa))
             u1 (+ centerx (* nx d1))
             u2 (+ centerx (* nx d2))
             v1 (- (* ny d1) centery)
             v2 (- (* ny d2) centery)]
         (Vec2. (* amount is (m/lerp u1 u2 p))
                (* amount is (m/lerp v1 v2 p))))))))

(defn rippled
  "Rippled"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d (+ m/EPSILON (v/magsq v))]
       (Vec2. (* (* (m/tanh d) (* 2.0 (.x v))) (/ amount 2.0))
              (* (* (m/cos d) (* 2.0 (.y v))) (/ amount 2.0)))))))


(defn rosewf
  ([] {:type :random
       :config (fn [] {:amp (u/sdrand 0.8 3.0)
                      :waves (r/randval (u/sirand 1 10) (u/sdrand 0.1 10))
                      :filled (r/randval 0.2 0.0 (r/drand))})})
  ([^double amount {:keys [^double amp ^double waves ^double filled]}]
   (fn [v]
     (let [a (v/heading v)
           r (* amount amp (m/cos (* waves a)))
           r (if (and (pos? filled)
                      (> filled (r/drand)))
               (* r (r/drand)) r)]
       (Vec2. (* r (m/cos a))
              (* r (m/sin a)))))))

(defn rosoni
  ([] {:type :regular
       :config (fn [] (let [mi (r/irand 1 30)]
                       {:maxiter mi
                        :sweetiter (r/irand mi)
                        :altshapes (r/brand)
                        :cutoff (u/sdrand 0.5 3.0)
                        :radius (u/sdrand 0.5 3.0)
                        :dx (r/drand -2.0 2.0)
                        :dy (r/randval 0.2 0.0 (r/drand -2.0 2.0))}))})
  ([^double amount {:keys [^long maxiter ^long sweetiter altshapes
                           ^double cutoff ^double radius ^double dx ^double dy]}]
   (let [phi (/ m/TWO_PI maxiter)
         sina (m/sin phi)
         cosa (m/cos phi)
         sweetiter (m/constrain sweetiter 0 (dec maxiter))
         radius2 (* radius radius)]
     (fn [^Vec2 v]
       (let [r (if (neg? cutoff)
                 (+ (max (m/abs (.x v)) (m/abs (.y v))) cutoff)
                 (- (v/mag v) cutoff))]
         (if (pos? r)
           (v/mult v amount)
           (let [y-dy (- (.y v) dy)
                 ^Vec3 sweet+cerc (loop [xrt (.x v)
                                         yrt (.y v)
                                         cerc (pos? dx)
                                         i (long 0)
                                         sweetx (.x v)
                                         sweety (.y v)]
                                    (if (< i maxiter)
                                      (let [r2 (if-not altshapes
                                                 (if (neg? radius)
                                                   (+ (max (m/abs (- xrt dx))
                                                           (m/abs (- yrt dy)))
                                                      radius)
                                                   (- (+ (m/sq (- xrt dx))
                                                         (m/sq (- yrt dy)))
                                                      radius2))
                                                 (if (neg? radius)
                                                   (+ (* (m/abs (m/atan2 y-dy (- xrt dx))) m/M_1_PI)
                                                      radius)
                                                   (let [diff (- xrt dx)]
                                                     (if (neg? diff)
                                                       (- diff)
                                                       (let [diff2 (* diff diff)]
                                                         (- (m/sq (- yrt dy))
                                                            (* diff2 (- radius2 diff2))))))))]
                                        (recur (- (* xrt cosa) (* yrt sina))
                                               (+ (* xrt sina) (* yrt cosa))
                                               (m/bool-xor cerc (<= r2 0.0))
                                               (inc i)
                                               (if (== i sweetiter) xrt sweetx)
                                               (if (== i sweetiter) yrt sweety)))
                                      (Vec3. sweetx sweety (if cerc 1.0 0.0))))]
             (if (pos? (.z sweet+cerc))
               (if (zero? sweetiter)
                 (if-not (zero? dy)
                   (Vec2. (* amount (- (.x sweet+cerc))) (* amount (- (.y sweet+cerc))))
                   (Vec2. (* amount (.x sweet+cerc)) (* amount (- (.y sweet+cerc)))))
                 (Vec2. (* amount (.x sweet+cerc)) (* amount (.y sweet+cerc))))
               (v/mult v amount)))))))))

(defn roundspher
  ([] {:type :regular})
  ([^double amount _]
   (let [s (m/sq m/M_2_PI)]
     (fn [v]
       (let [d (v/magsq v)
             re (/ (+ s (/ d)))
             ad (/ amount d)]
         (v/mult v (* amount ad re)))))))

