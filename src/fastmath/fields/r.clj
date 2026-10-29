(ns fastmath.fields.r
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn r-circleblur
  "R Circle Blur"
  ([] {:type :random
       :config (fn [] {:n (u/sdrand 0.9 7.0)
                      :seed (r/drand Float/MAX_VALUE)
                      :dist (r/drand -1.0 1.0)
                      :mn (r/drand -3.0 3.0)
                      :mx (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double n ^double seed ^double dist ^double mn ^double mx]}]
   (let [dm (m/- mx mn)]
     (fn [v]
       (let [angle (v/heading v)
             rad (v/mag v)
             rad (m/mod rad n)
             by (m/sin (m/+ angle rad))
             bx (m/cos (m/+ angle rad))
             by (m/round (m/* by rad))
             bx (m/round (m/* bx rad))
             rad2 (m/* 0.5 (m/sqrt (r/drand)))
             angle2 (r/drand m/TWO_PI)
             a1 (m/sfrac (m/* 43758.5453 (m/sin (m/+ (m/* bx 127.1) (m/* by 311.7) seed))))
             a2 (m/sfrac (m/* 43758.5453 (m/sin (m/+ (m/* bx 269.5) (m/* by 183.3) seed))))
             a3 (m/sfrac (m/* 43758.5453 (m/sin (m/+ (m/* bx 78.233) (m/* by 12.9898) seed))))
             a3 (m/+ mn (m/* a3 dm))
             rad2 (m/* rad2 a3)]
         (Vec2. (m/* amount (m/+ bx (m/* rad2 (m/cos angle2)) (m/* dist a1)))
                (m/* amount (m/+ by (m/* rad2 (m/sin angle2)) (m/* dist a2)))))))))

(defn radialblur
  "Radial blur"
  ([] {:type :random
       :config (fn [] {:angle (r/drand (m/- m/TWO_PI) m/TWO_PI)})})
  ([^double amount {:keys [^double angle]}]
   (let [spin (m/* amount (m/sin (m/* angle m/HALF_PI)))
         zoom (m/* amount (m/cos (m/* angle m/HALF_PI)))]
     (fn [^Vec2 v]
       (let [rnd-g (m/+ (r/drand) (r/drand) (r/drand) (r/drand) -2.0)
             ra (v/mag v)
             alpha (m/+ (m/* spin rnd-g) (v/heading v))
             rz (m/dec (m/* zoom rnd-g))]
         (Vec2. (m/+ (m/* rz (.x v)) (m/* ra (m/cos alpha)))
                (m/+ (m/* rz (.y v)) (m/* ra (m/sin alpha)))))))))

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
           zt3 (m/- (.x cb) (m/* 3.0 (.x v) (.y sqr)))
           zt2 (m/- (.x sqr) (.y sqr))
           zb3 (m/- (m/* 3.0 (.x sqr) (.y v)) (.y cb))
           zb2 (m/* 2.0 (.x v) (.y v))
           tr (m/+ (m/* a zt3) (m/* b zt2) (m/* c (.x v)) d)
           ti (m/+ (m/* a zb3) (m/* b zb2) (m/* c (.y v )))
           br (m/+ (m/* e zt3) (m/* f zt2) (m/* g (.x v)) h)
           bi (m/+ (m/* e zb3) (m/* f zb2) (m/* g (.y v )))
           r3den (m// amount (m/+ (m/* br br) (m/* bi bi)))]
       (Vec2. (m/* r3den (m/+ (m/* tr br) (m/* ti bi)))
              (m/* r3den (m/- (m/* ti br) (m/* tr bi))))))))

(defn rays1
  "Rays1"
  ([] {:type :regular})
  ([^double amount _]
   (let [pa (m/* amount (m/sq m/M_2_PI))]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             u (m/+ pa (m// (m/tan (m/sqrt t))))
             r (m/* amount u t)]
         (Vec2. (m// r (.x v))
                (m// r (.y v))))))))

(defn rays2
  "Rays2"
  ([] {:type :regular})
  ([^double amount _]
   (let [a10 (m// amount 10.0)]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             u (m// (m/cos (m/* (m/+ t m/EPSILON) (m/tan (m// (m/+ m/EPSILON t))))))
             r (m/* a10 t u)]
         (Vec2. (m// r (.x v))
                (m// r (.y v))))))))

(defn rays3
  "Rays3"
  ([] {:type :regular})
  ([^double amount _]
   (let [a10 (m// amount 10.0)]
     (fn [^Vec2 v]
       (let [t (v/magsq v)
             t2 (m/* t t)
             u (m// (m/sqrt (m/cos (m/sin (m/* (m/+ m/EPSILON t2) (m/sin (m// (m/+ t2 m/EPSILON))))))))
             r (m/* a10 t u)]
         (Vec2. (m// (m/* r (m/cos t)) (.x v))
                (m// (m/* r (m/tan t)) (.y v))))))))

(defn rays
  "Rays"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [ang (m/* amount (r/drand m/PI))
           r (m// amount (m/+ m/EPSILON (v/magsq v)))
           tanr (m/* amount r (m/tan ang))]
       (Vec2. (m/* tanr (m/cos (.x v)))
              (m/* tanr (m/sin (.y v))))))))

(defn rectangles
  "Rectangles"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (Vec2. (if (m/< (m/abs (.x v)) m/EPSILON)
              (m/* amount (.x v))
              (m/* amount (-> (.x v)
                            (m// x)
                            (m/floor)
                            (m/* 2.0)
                            m/inc
                            (m/* x)
                            (m/- (.x v)))))
            (if (m/< (m/abs (.y v)) m/EPSILON)
              (m/* amount (.y v))
              (m/* amount (-> (.y v)
                            (m// y)
                            (m/floor)
                            (m/* 2.0)
                            m/inc
                            (m/* y)
                            (m/- (.y v)))))))))

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
         k (m// kn kd)
         cycles-to-close (double (if (m/zero? (m/mod k 1.0))
                                   (if (m/zero? (m/mod k 2.0))
                                     1.0
                                     (if (or (not (m/zero? radial-offset))
                                             (not (m/zero? inner-spread))
                                             (not (m/zero? outer-spread))
                                             (not (m/zero? fill)))
                                       1.0
                                       0.5))
                                   (if (and (m/zero? (m/mod kn 1.0))
                                            (m/zero? (m/mod kd 1.0)))
                                     (let [lkn (long kn)
                                           lkd (long kd)
                                           gcd (m/gcd lkn lkd)
                                           [^long kn ^long kd] (if (m/not== gcd 1)
                                                                 [(m// lkn gcd) (m// lkd gcd)]
                                                                 [lkn lkd])]
                                       (if (or (m/zero? (m/mod kn 2.0))
                                               (m/zero? (m/mod kd 2.0)))
                                         kd
                                         (m// kd 2)))
                                     (if (m/< cycles-param 16)
                                       16
                                       (m/* 2 kd kn)))))
         cycles (if (m/zero? cycles-param)
                  (m/* cycles-to-close metacycles)
                  cycles-param)]
     (fn [^Vec2 v]
       (let [rin (m/* spread-split (v/mag v))
             tin (v/heading v)
             t (m/* cycles (m/+ tin (m/* cycle-offset m/TWO_PI)))
             r (m/+ radial-offset (m/cos (m/* t k)))
             r (if-not (m/zero? fill)
                 (m/+ r (m/* fill (m/- (r/drand) 0.5)))
                 r)
             x (m/* r (m/cos t))
             y (m/* r (m/sin t))
             expansion (m/floor (m// (m/* cycles (m/+ tin m/PI))
                                   (m/* cycles-to-close m/TWO_PI)))
             adjusted-amount (m/+ amount (m/* expansion metacycle-expansion))]
         (if (m/> (m/abs rin) (m/abs r))
           (case (long outer-mode)
             0 (Vec2. (m/* adjusted-amount x)
                      (m/* adjusted-amount y))
             1 (let [rinx (m/inc (m/* (m/dec rin) outer-spread outer-spread-ratio))
                     riny (m/inc (m/* (m/dec rin) outer-spread))]
                 (Vec2. (m/* adjusted-amount rinx x)
                        (m/* adjusted-amount riny y)))
             2 (let [xin (m/* (m/sgn x) (m/abs (.x v)))
                     yin (m/* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (m/* adjusted-amount (m/+ x (m/* outer-spread outer-spread-ratio (m/- xin x))))
                        (m/* adjusted-amount (m/+ y (m/* outer-spread (m/- yin y))))))
             3 (let [xin (m/* (m/sgn x) (m/abs (.x v)))
                     yin (m/* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (m/* adjusted-amount (m/+ x (m/* outer-spread outer-spread-ratio xin)))
                        (m/* adjusted-amount (m/+ y (m/* outer-spread yin)))))
             4 (let [rinx (m/+ (m/* 0.5 rin) (m/* outer-spread outer-spread-ratio))
                     riny (m/+ (m/* 0.5 rin) outer-spread)]
                 (Vec2. (m/* adjusted-amount rinx x)
                        (m/* adjusted-amount riny y)))
             5 v
             6 u/zerov)
           (case (long inner-mode)
             0 (Vec2. (m/* adjusted-amount x)
                      (m/* adjusted-amount y))
             1 (let [rinx (m/inc (m/* (m/dec rin) inner-spread inner-spread-ratio))
                     riny (m/inc (m/* (m/dec rin) inner-spread))]
                 (Vec2. (m/* adjusted-amount rinx x)
                        (m/* adjusted-amount riny y)))
             2 (let [xin (m/* (m/sgn x) (m/abs (.x v)))
                     yin (m/* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (m/* adjusted-amount (m/+ x (m/* inner-spread inner-spread-ratio (m/- xin x))))
                        (m/* adjusted-amount (m/+ y (m/* inner-spread (m/- yin y))))))
             3 (let [xin (m/* (m/sgn x) (m/abs (.x v)))
                     yin (m/* (m/sgn y) (m/abs (.y v)))]
                 (Vec2. (m/* adjusted-amount (m/+ x (m/* inner-spread inner-spread-ratio xin)))
                        (m/* adjusted-amount (m/+ y (m/* inner-spread yin)))))
             4 (let [rinx (m/+ (m/* 0.5 rin) (m/* inner-spread inner-spread-ratio))
                     riny (m/+ (m/* 0.5 rin) inner-spread)]
                 (Vec2. (m/* adjusted-amount rinx x)
                        (m/* adjusted-amount riny y)))
             5 v
             6 u/zerov)))))))

(defn rings2
  "Rings2"
  ([] {:type :regular
       :config (fn [] {:val (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double val]}]
   (let [dx (m/+ m/EPSILON (m/sq val))]
     (fn [^Vec2 v]
       (let [l (v/mag v)
             r (m/* amount (m/- 2.0 (m/* dx (m/inc (m// (m/* 2.0 (double (int (m/* 0.5 (m/inc (m// l dx)))))) l)))))]
         (v/mult v r))))))

(defn rings
  "Rings"
  ([] {:type :regular
       :config (fn [] {:coeff20 (r/drand 1.3)})})
  ([^double amount {:keys [^double coeff20]}]
   (let [dx (m/+ m/EPSILON (m/sq coeff20))
         dx2 (m/+ dx dx)
         rdx (m// 1.0 dx2)
         dx- (m/- 1.0 dx)]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             rr (m/* amount (m/+ (m/- r (m/* dx2 (double (int (m/* (m/+ r dx) rdx))))) (m/* r dx-)))]
         (Vec2. (m/* rr (m// (.x v) r))
                (m/* rr (m// (.y v) r))))))))

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
   (let [f (m/* frequency 5.0)
         a (m/* amplitude 0.01)
         p (m/- (m/* phase m/TWO_PI) m/PI)
         s (if (m/zero? scale) m/EPSILON scale)
         is (m// s)
         vxp (m/* velocity p)
         pxa (m/* p a)
         pixa (m/* (m/- m/PI p) a)]
     (fn [^Vec2 v]
       (let [x (m/- (m/* s (.x v)) centerx)
             y (m/- (m/* s (.y v)) centery)
             d (case (unchecked-int fixed-dist-calc)
                 0 (v/mag v)
                 1 (m/sqrt (m/* (m/sq (.x v)) (m/sq (.y v))))
                 2 (m/max (m/abs (.x v)) (m/abs (.y v)))
                 3 (m/+ (m/abs (.x v)) (m/abs (.y v))))
             d (if (m/< d m/EPSILON) m/EPSILON d)
             nx (m// x d)
             ny (m// y d)
             wave (m/cos (m/- (m/* f d) vxp))
             d1 (m/+ d (m/* wave pxa))
             d2 (m/+ d (m/* wave pixa))
             u1 (m/+ centerx (m/* nx d1))
             u2 (m/+ centerx (m/* nx d2))
             v1 (m/- (m/* ny d1) centery)
             v2 (m/- (m/* ny d2) centery)]
         (Vec2. (m/* amount is (m/lerp u1 u2 p))
                (m/* amount is (m/lerp v1 v2 p))))))))

(defn rippled
  "Rippled"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d (m/+ m/EPSILON (v/magsq v))]
       (Vec2. (m/* (m/* (m/tanh d) (m/* 2.0 (.x v))) (m// amount 2.0))
              (m/* (m/* (m/cos d) (m/* 2.0 (.y v))) (m// amount 2.0)))))))


(defn rosewf
  ([] {:type :random
       :config (fn [] {:amp (u/sdrand 0.8 3.0)
                      :waves (r/randval (u/sirand 1 10) (u/sdrand 0.1 10))
                      :filled (r/randval 0.2 0.0 (r/drand))})})
  ([^double amount {:keys [^double amp ^double waves ^double filled]}]
   (fn [v]
     (let [a (v/heading v)
           r (m/* amount amp (m/cos (m/* waves a)))
           r (if (and (m/pos? filled)
                      (m/> filled (r/drand)))
               (m/* r (r/drand)) r)]
       (Vec2. (m/* r (m/cos a))
              (m/* r (m/sin a)))))))

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
   (let [phi (m// m/TWO_PI maxiter)
         sina (m/sin phi)
         cosa (m/cos phi)
         sweetiter (m/constrain sweetiter 0 (m/dec maxiter))
         radius2 (m/* radius radius)]
     (fn [^Vec2 v]
       (let [r (if (m/neg? cutoff)
                 (m/+ (m/max (m/abs (.x v)) (m/abs (.y v))) cutoff)
                 (m/- (v/mag v) cutoff))]
         (if (m/pos? r)
           (v/mult v amount)
           (let [y-dy (m/- (.y v) dy)
                 ^Vec3 sweet+cerc (loop [xrt (.x v)
                                         yrt (.y v)
                                         cerc (m/pos? dx)
                                         i (long 0)
                                         sweetx (.x v)
                                         sweety (.y v)]
                                    (if (m/< i maxiter)
                                      (let [r2 (if-not altshapes
                                                 (if (m/neg? radius)
                                                   (m/+ (m/max (m/abs (m/- xrt dx))
                                                           (m/abs (m/- yrt dy)))
                                                      radius)
                                                   (m/- (m/+ (m/sq (m/- xrt dx))
                                                         (m/sq (m/- yrt dy)))
                                                      radius2))
                                                 (if (m/neg? radius)
                                                   (m/+ (m/* (m/abs (m/atan2 y-dy (m/- xrt dx))) m/M_1_PI)
                                                      radius)
                                                   (let [diff (m/- xrt dx)]
                                                     (if (m/neg? diff)
                                                       (m/- diff)
                                                       (let [diff2 (m/* diff diff)]
                                                         (m/- (m/sq (m/- yrt dy))
                                                            (m/* diff2 (m/- radius2 diff2))))))))]
                                        (recur (m/- (m/* xrt cosa) (m/* yrt sina))
                                               (m/+ (m/* xrt sina) (m/* yrt cosa))
                                               (m/bool-xor cerc (m/<= r2 0.0))
                                               (m/inc i)
                                               (if (m/== i sweetiter) xrt sweetx)
                                               (if (m/== i sweetiter) yrt sweety)))
                                      (Vec3. sweetx sweety (if cerc 1.0 0.0))))]
             (if (m/pos? (.z sweet+cerc))
               (if (m/zero? sweetiter)
                 (if-not (m/zero? dy)
                   (Vec2. (m/* amount (m/- (.x sweet+cerc))) (m/* amount (m/- (.y sweet+cerc))))
                   (Vec2. (m/* amount (.x sweet+cerc)) (m/* amount (m/- (.y sweet+cerc)))))
                 (Vec2. (m/* amount (.x sweet+cerc)) (m/* amount (.y sweet+cerc))))
               (v/mult v amount)))))))))

(defn roundspher
  ([] {:type :regular})
  ([^double amount _]
   (let [s (m/sq m/M_2_PI)]
     (fn [v]
       (let [d (v/magsq v)
             re (m// (m/+ s (m// d)))
             ad (m// amount d)]
         (v/mult v (m/* amount ad re)))))))

