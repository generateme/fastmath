(ns fastmath.fields.c
  (:refer-clojure :exclude [chunk])
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn cpow2
  "CPow2"
  ([] {:type :random
       :config (fn [] {:r (u/sdrand 0.1 1.5)
                      :a (u/sdrand 0.1 3.1)
                      :divisor (u/sdrand 0.2 2.0)
                      :range (r/randval (u/sdrand 0.1 10.0)
                                        (let [v (double (r/irand 1 10))
                                              v (r/randval v (m// v))]
                                          (r/randval (m/* -1.0 v) v)))})})
  ([^double amount {:keys [^double r ^double a ^double divisor ^double range]}]
   (let [ang (m// m/TWO_PI divisor)
         c (m// (m/* r (m/cos (m/* m/HALF_PI a))) divisor)
         d (m// (m/* r (m/sin (m/* m/HALF_PI a))) divisor)
         half-c (m/* 0.5 c)
         half-d (m/* 0.5 d)
         inv-range (m// 0.5 range)
         full-range (m/* m/TWO_PI range)]
     (fn [^Vec2 v]
       (let [ai (v/heading v)
             n (long (r/drand range))
             n (if (m/neg? ai) (m/inc n) n)
             ai (m/+ ai (m/* m/TWO_PI n))
             ai (if (m/< (m/cos (m/* ai inv-range)) (r/drand -1.0 1.0)) (m/- ai full-range) ai)
             lnr2 (m/log (v/magsq v))
             r (m/* amount (m/exp (m/- (m/* half-c lnr2) (m/* d ai))))
             th (m/+ (m/* c ai) (m/* half-d lnr2) (m/* ang (m/floor (r/drand divisor))))]
         (Vec2. (m/* r (m/cos th))
                (m/* r (m/sin th))))))))


(defn cpow3
  "CPow3"
  ([] {:type :random
       :config (fn [] {:r (u/sdrand 0.1 1.5)
                      :d (u/sdrand 0.1 24.0)
                      :divisor (u/sdrand 0.5 1.5)
                      :spread (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double r ^double d ^double divisor ^double spread]}]
   (let [ang (m// m/TWO_PI divisor)
         pa (m/atan2 (if (m/neg? d) (m/- (m/log (m/- d))) (m/log d)) m/TWO_PI)
         cpa (m/cos pa)
         spa (m/sin pa)
         tc (m// (m/* r cpa cpa) divisor)
         td (m// (m/* r cpa spa) divisor)
         half-c (m/* 0.5 tc)
         half-d (m/* 0.5 td)
         coeff (if (m/zero? td) 0.0 (m// (m/* -0.095 spread) td))
         fac (m/* half-d ang)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             a (if (m/neg? a) (m/+ a m/TWO_PI) a)
             a (if (m/< (m/cos (m// a 2.0)) (r/drand -1.0 1.0)) (m/- a m/TWO_PI) a)
             a (m/+ a (m/* (r/randval m/TWO_PI m/-TWO_PI) coeff (m/round (m/log (r/drand)))))
             lnr2 (m/log (v/magsq v))
             r (m/* amount (m/exp (m/- (m/* half-c lnr2) (m/* td a))))
             th (m/+ (m/* tc a) (m/* fac lnr2) (m/* ang (m/floor (r/drand divisor))))]
         (Vec2. (m/* r (m/cos th))
                (m/* r (m/sin th))))))))

(defn cpow3wf
  "CPow3WF"
  ([] {:type :random
       :config (fn [] {:r (u/sdrand 0.1 1.5)
                      :a (u/sdrand 0.1 3.1)
                      :divisor (u/sdrand 0.5 1.5)
                      :spread (r/drand -3.0 3.0)
                      :discrete-spread (r/brand 0.7)
                      :spread2 (r/drand -0.5 0.5)
                      :offset2 (u/sdrand 0.15 3.5)})})
  ([^double amount {:keys [^double r ^double a ^double divisor ^double spread
                           discrete-spread ^double spread2 ^double offset2]}]
   (let [ang (m// m/TWO_PI divisor)
         c (m// (m/* r (m/cos (m/* m/HALF_PI a))) divisor)
         d (m// (m/* r (m/sin (m/* m/HALF_PI a))) divisor)
         half-c (m/* 0.5 c)
         half-d (m/* 0.5 d)
         inv-spread (m// 0.5 spread)
         full-spread (m/* m/TWO_PI spread)
         fac (m/* c half-d ang)]
     (fn [^Vec2 v]
       (let [ai (v/heading v)
             n (r/drand spread)
             n (double (if discrete-spread (int n) n))
             n (if (m/neg? ai) (m/inc n) n)
             ai (m/+ ai (m/* m/TWO_PI n))
             ai (if (m/< (m/cos (m/* ai inv-spread)) (r/drand -1.0 1.0)) (m/- ai full-spread) ai)
             lnr2 (m/log (v/magsq v))
             ri (m/* amount (m/exp (m/- (m/* half-c lnr2) (m/* d ai))))
             ang2 (m/* fac ai lnr2 (m/+ (r/drand spread2) offset2))]
         (Vec2. (m/* ri (m/cos ang2))
                (m/* ri (m/sin ang2))))))))

(defn cpow
  "CPow"
  ([] {:type :random
       :config (fn [] {:r (r/drand -2.0 2.0)
                      :i (r/drand -2.0 2.0)
                      :power (u/sdrand 0.1 12)})})
  ([^double amount {:keys [^double r ^double i ^double power]}]
   (let [va (m// m/TWO_PI power)
         vc (m// r power)
         vd (m// i power)]
     (fn [v]
       (let [a (v/heading v)
             lnr (m/* 0.5 (m/log (v/magsq v)))
             ang (m/+ (m/* a vc)
                    (m/* vd lnr)
                    (m/* va (m/floor (r/drand power))))
             m (m/* amount (m/exp (m/- (m/* vc lnr)
                                   (m/* vd a))))]
         (Vec2. (m/* m (m/cos ang)) (m/* m (m/sin ang))))))))

(defn csin
  "CSin by zephyrtronium, http://fractal-resources.deviantart.com/art/CSin-Apophysis-Plugin-158332287"
  ([] {:type :regular
       :config (fn [] {:stretch (u/sdrand 0.3 3.0)})})
  ([^double amount {:keys [^double stretch]}]
   (let [s-cx (Vec2. stretch 0.0)]
     (fn [^Vec2 v]
       (v/mult (->> v
                    (c/mult s-cx)
                    (c/sin)) amount)))))


(defn cannabiswf
  "CannabisWF"
  ([] {:type :pattern
       :config (fn [] {:filled (r/randval 0.3 0.0 (r/randval 0.3 1.0 (r/drand)))})})
  ([^double amount {:keys [^double filled]}]
   (fn [v]
     (let [a (v/heading v)
           r (m/* (m/inc (m/* 0.9 (m/cos (m/* 8.0 a))))
                (m/inc (m/* 0.1 (m/cos (m/* 24.0 a))))
                (m/+ 0.9 (m/* 0.1 (m/cos (m/* 200.0 a))))
                (m/inc (m/sin a)))
           r (m/* amount (if (m/> filled (r/drand)) (m/* r (r/drand)) r))
           a (m/+ a m/HALF_PI)]
       (Vec2. (m/* r (m/sin a))
              (m/* r (m/cos a)))))))

(defn cardioid
  "Cardioid"
  ([] {:type :regular
       :config (fn [] {:a (r/randval 0.7 (u/sirand 1 9) (u/sdrand 0.5 9.0))})})
  ([^double amount {:keys [^double a]}]
   (fn [^Vec2 v]
     (let [ai (v/heading v)
           r (m/* amount (m/sqrt (m/+ (v/magsq v) (m/sin (m/* a ai)) 1.0)))]
       (Vec2. (m/* r (m/cos ai))
              (m/* r (m/sin ai)))))))

(defn cell
  "Cell"
  ([] {:type :regular
       :config (fn [] {:size (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double size]}]
   (let [inv-cell-size (m// 1.0 size)
         av (Vec2. amount (m/- amount))]
     (fn [v]
       (let [^Vec2 xy (-> (v/mult v inv-cell-size)
                          (v/floor))
             dxy (v/sub v (v/mult xy size))
             newxy (if-not (m/neg? (.y xy))
                     (if-not (m/neg? (.x xy))
                       (v/mult xy 2.0)
                       (Vec2. (m/- (m/inc (m/* 2.0 (.x xy))))
                              (m/* 2.0 (.y xy))))
                     (if-not (m/neg? (.x xy))
                       (Vec2. (m/* 2.0 (.x xy))
                              (m/- (m/inc (m/* 2.0 (.y xy)))))
                       (Vec2. (m/- (m/inc (m/* 2.0 (.x xy))))
                              (m/- (m/inc (m/* 2.0 (.y xy)))))))]
         (-> newxy
             (v/mult size)
             (v/add dxy)
             (v/emult av)))))))

(defn checks
  "Checks"
  ([] {:type :random
       :config (fn [] {:x (r/drand -3.0 3.0)
                      :y (r/drand -3.0 3.0)
                      :size (r/drand -2.0 2.0)
                      :rnd (r/randval 0.0 (u/sdrand 0.1 1.0))})})
  ([^double amount {:keys [^double x ^double y ^double size ^double rnd]}]
   (let [cs (m// 1.0 (m/+ size m/EPSILON))
         ncx (m/* -1.0 x)
         ncy (m/* -1.0 y)]
     (fn [^Vec2 v]
       (let [isxy (m/long-add (m/round (m/* cs (.x v)))
                          (m/round (m/* cs (.y v))))
             dxy (if (m/even? isxy)
                   (Vec2. (m/+ ncx (r/drand rnd)) ncy)
                   (Vec2. x (m/+ y (r/drand rnd))))]
         (-> v
             (v/add dxy)
             (v/mult amount)))))))

(defn chrysantemum
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [u (r/drand (m/* 21.0 m/PI))
           p4 (m/sin (m/* u 5.666666666666667))
           p4 (m/* p4 p4)
           p4 (m/* p4 p4)
           p8 (m/sin (m/- (m/* 2.0 (m/cos (m/* 3.0 u))) (m/* 28.0 u)))
           p8 (m/* p8 p8)
           p8 (m/* p8 p8)
           p8 (m/* p8 p8)
           r (m/* m/THIRD amount (m/- (m/* 5.0 (m/inc (m/sin (m/* u 2.2))))
                                  (m/* 4.0 p4 p8)))]
       (Vec2. (m/* r (m/cos u))
              (m/* r (m/sin u)))))))

(defn chunk
  "Chunk, by zephyrtronium https://zephyrtronium.deviantart.com/art/Chunk-Apophysis-Plugin-Pack-182375397"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -1.2 1.2)
                      :b (r/drand -1.2 1.2)
                      :c (r/drand -1.2 1.2)
                      :d (r/drand -1.2 1.2)
                      :e (r/drand -1.2 1.2)
                      :f (r/drand -1.2 1.2)
                      :mode (r/brand)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d ^double e ^double f mode]}]
   (fn [^Vec2 v]
     (let [r (m/* amount (m/+ (m/* a (m/sq (.x v)))
                          (m/* b (.x v) (.y v))
                          (m/* c (m/sq (.y v)))
                          (m/* d (.x v))
                          (m/* e (.y v))
                          f))]
       (if mode
         (if-not (m/pos? r) v u/zerov)
         (if (m/pos? r) v u/zerov))))))

(defn circleblur
  "Circle blur"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [rad (m/sqrt (r/drand))
           a (r/drand m/TWO_PI)]
       (Vec2. (m/* amount (m/cos a) rad)
              (m/* amount (m/sin a) rad))))))

(defn circlecrop
  "CircleCrop"
  ([] {:type :random
       :config (fn [] {:radius (r/drand 0.1 2.0)
                      :x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)
                      :scatter-area (r/drand -1.0 1.0)
                      :zero (r/brand)})})
  ([^double amount {:keys [^double radius ^double x ^double y ^double scatter-area zero]}]
   (let [ca (m/constrain scatter-area -1.0 1.0)
         v0 (Vec2. x y)]
     (fn [^Vec2 v]
       (let [^Vec2 v (v/sub v v0)
             rad (v/mag v)
             ang (v/heading v)
             rdc (m/+ radius (r/drand (m/* 0.5 ca)))
             s (m/sin ang)
             c (m/cos ang)]
         (cond
           (m/<= rad radius) (Vec2. (m/+ x (m/* amount (.x v)))
                                  (m/+ y (m/* amount (.y v))))
           (and (not zero)
                (m/> rad radius)) (Vec2. (m/+ x (m/* amount rdc c))
                                       (m/+ y (m/* amount rdc s)))
           :else v))))))

(defn circlelinear
  "CircleLinear by eralex, http://eralex61.deviantart.com/art/Circles-Plugins-126273412"
  ([] {:type :regular
       :config (fn [] {:Sc (r/drand)
                      :K (r/drand -2.0 2.0)
                      :Dens1 (r/drand)
                      :Dens2 (r/drand)
                      :Reverse (r/drand -1.0 1.0)
                      :Seed (r/irand)})})
  ([^double amount {:keys [^double Sc ^double K ^double Dens1 ^double Dens2 ^double Reverse ^long Seed]}]
   (let [dd (m/* Dens1 Dens2)]
     (fn [^Vec2 v]
       (let [M (->> Sc
                    (m// (.x v))
                    (m/* 0.5)
                    (m/floor))
             N (->> Sc
                    (m// (.y v))
                    (m/* 0.5)
                    (m/floor))
             X (m/- (.x v) (->> M
                              (m/* 2.0)
                              (m/inc)
                              (m/* Sc)))
             Y (m/- (.y v) (->> N
                              (m/* 2.0)
                              (m/inc)
                              (m/* Sc)))
             U (m/hypot X Y)
             V (->> (r/discrete-noise (unchecked-long (m/+ M 10.0)) (unchecked-long (m/+ N 3.0)))
                    (m/* 0.7)
                    (m/+ 0.3)
                    (m/* Sc))
             Z1 (r/discrete-noise (unchecked-long (m/+ M Seed)) N)
             ^Vec2 XY (if (and (m/< Z1 Dens1) (m/< U V)) 
                        (if (m/pos? Reverse)
                          (if (m/< Z1 dd)
                            (Vec2. (m/* K X) (m/* K Y))
                            (let [Z (->> K
                                         (m/- 1.0)
                                         (m/* U)
                                         (m// V)
                                         (m/+ K))]
                              (Vec2. (m/* Z X) (m/* Z Y))))
                          (if (m/> Z1 dd)
                            (Vec2. (m/* K X) (m/* K Y))
                            (let [Z (->> K
                                         (m/- 1.0)
                                         (m/* U)
                                         (m// V)
                                         (m/+ K))]
                              (Vec2. (m/* Z X) (m/* Z Y)))))
                        (Vec2. X Y))]
         (Vec2. (->> 2.0
                     (m/* M)
                     (m/inc)
                     (m/* Sc)
                     (m/+ (.x XY))
                     (m/* amount))
                (->> 2.0
                     (m/* N)
                     (m/inc)
                     (m/* Sc)
                     (m/+ (.y XY))
                     (m/* amount))))))))

(defn circlerand
  "Circle Rand http://eralex61.deviantart.com/art/Circles-Plugins-126273412"
  ([] {:type :random
       :config (fn [] {:Sc (r/drand 2.0)
                      :Dens (r/drand)
                      :X (u/sdrand 0.5 12.0)
                      :Y (u/sdrand 0.5 12.0)
                      :Seed (r/irand)})})
  ([^double amount {:keys [^double Sc ^double Dens ^double X ^double Y ^double Seed]}]
   (let [xy (Vec2. X Y)] 
     (fn [_]
       (loop [iter (int 0)]
         (let [XY (-> (v/generate-vec2 #(r/drand -1.0 1.0))
                      (v/emult xy))
               ^Vec2 MN (-> XY
                            (v/mult 0.5)
                            (v/div Sc)
                            (v/fmap m/floor))
               XY (v/sub XY (-> MN
                                (v/mult 2.0)
                                (v/fmap #(m/inc ^double %))
                                (v/mult Sc)))]
           (if (and (m/< iter 60)
                    (or (m/> (r/discrete-noise (unchecked-long (m/+ Seed (.x MN))) (.y MN)) Dens)
                        (m/> (v/mag XY) (-> (r/discrete-noise (unchecked-long (m/+ 10.0 (.x MN))) (unchecked-long (m/+ 3.0 (.y MN))))
                                          (m/* 0.7)
                                          (m/+ 0.3)
                                          (m/* Sc)))))
             (recur (m/inc iter))
             (-> MN
                 (v/mult 2.0)
                 (v/fmap #(m/inc ^double %))
                 (v/mult Sc)
                 (v/add XY)
                 (v/mult amount)))))))))

(defn circletrans1
  ([] {:type :random
       :config (fn [] {:Sc (r/drand)
                      :X (r/drand -5.0 5.0)
                      :Y (r/drand -5.0 5.0)
                      :Dens (r/drand)
                      :Seed (r/irand)})})
  ([^double amount {:keys [^double Sc ^double X ^double Y ^double Dens ^long Seed]}]
   (let [AB (Vec2. X Y)
         aAB (v/abs AB)]
     (fn [^Vec2 v]
       (let [^Vec2 Uxy (v/add (v/mult (v/sub v AB) 0.5) AB)
             ^Vec2 MN (v/floor (v/div (v/mult Uxy 0.5) Sc))
             M (long (.x MN))
             N (long (.y MN))
             XY (v/sub Uxy (v/mult (v/shift (v/mult MN 2.0) 1.0) Sc))
             U (v/mag XY)]
         (v/mult (if-not (or (m/> (r/discrete-noise (m/long-add M Seed) N) Dens)
                             (m/> U (m/* Sc (m/+ 0.3 (m/* 0.7 (r/discrete-noise (m/long-add M 10) (m/long-add N 3)))))))
                   (loop [i (long 0)]
                     (let [XY (v/emult aAB (Vec2. (r/drand -1.0 1.0) (r/drand -1.0 1.0)))
                           ^Vec2 MN (v/floor (v/div (v/mult XY 0.5) Sc))
                           M (long (.x MN))
                           N (long (.y MN))]
                       (if (and (m/< i 100) (m/> (r/discrete-noise (m/long-add M Seed) N) Dens))
                         (recur (m/inc i))
                         (let [alpha (r/drand m/TWO_PI)
                               s (m/sin alpha)
                               c (m/cos alpha)
                               U (m/+ 0.3 (m/* 0.7 (r/discrete-noise (m/long-add M 10) (m/long-add N 3))))]
                           (v/add (Vec2. (m/* U c) (m/* U s))
                                  (v/mult (v/shift (v/mult MN 2.0) 1.0) Sc))))))
                   Uxy) amount))))))

(defn circlesplit
  ([] {:type :regular
       :config (fn [] (let [r (r/drand 0.1 2.0)]
                       {:radius r
                        :split (r/drand -2.0 r)}))})
  ([^double amount {:keys [^double radius ^double split]}]
   (let [diff (m/- radius split)]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (v/mult (if (m/< r diff)
                   v
                   (let [a (v/heading v)
                         len (m/+ r split)]
                     (Vec2. (m/* len (m/cos a))
                            (m/* len (m/sin a))))) amount))))))

(defn circlize2
  ([] {:type :regular
       :config (fn [] {:hole (r/drand -2.2 2.2)})})
  ([^double amount {:keys [^double hole]}]
   (fn [^Vec2 v]
     (let [absx (m/abs (.x v))
           absy (m/abs (.y v))
           side (if (m/>= absx absy) absx absy)
           perimeter (if (m/>= absx absy)
                       (if (m/> (.x v) absy)
                         (m/+ absx (.y v))
                         (m/- (m/* 5.0 absx) (.y v)))
                       (if (m/> (.y v) absx)
                         (m/- (m/* 3.0 absy) (.x v))
                         (m/+ (m/* 7.0 absy) (.x v))))
           r (m/* amount (m/+ side hole))
           a (m/- (m// (m/* m/M_PI_4 perimeter) side) m/M_PI_4)
           sina (m/sin a)
           cosa (m/cos a)]
       (Vec2. (m/* r cosa) (m/* r sina))))))

(defn circlize
  ([] {:type :regular
       :config (fn [] {:hole (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double hole]}]
   (let [var4-pi (m// amount m/M_PI_4)]
     (fn [^Vec2 v]
       (let [absx (m/abs (.x v))
             absy (m/abs (.y v))
             side (if (m/>= absx absy) absx absy)
             perimeter (if (m/>= absx absy)
                         (if (m/> (.x v) absy)
                           (m/+ absx (.y v))
                           (m/- (m/* 5.0 absx) (.y v)))
                         (if (m/> (.y v) absx)
                           (m/- (m/* 3.0 absy) (.x v))
                           (m/+ (m/* 7.0 absy) (.x v))))
             r (m/+ (m/* var4-pi side) hole)
             a (m/- (m// (m/* m/M_PI_4 perimeter) side) m/M_PI_4)
             sina (m/sin a)
             cosa (m/cos a)]
         (Vec2. (m/* r cosa) (m/* r sina)))))))

(defn circular2
  ([] {:type :random
       :config (fn [] {:angle (r/drand -45.0 45.0)
                      :seed (r/drand m/TWO_PI)
                      :xx (r/drand -100.0 100.0)
                      :yy (r/drand -100.0 100.0)})})
  ([^double amount {:keys [^double angle ^double seed ^double xx ^double yy]}]
   (let [ca (m/radians angle)]
     (fn [^Vec2 v]
       (let [aux (m/frac (m/* 43758.5453 (m/sin (m/+ (m/* xx (.x v))
                                                 (m/* yy (.y v))
                                                 seed))))
             rnd (m/* (m/- (m/* 2.0 (m/+ (r/drand) aux)) 2.0) ca)
             rad (m/* amount (v/mag v))
             ang (m/+ rnd (v/heading v))
             by (m/sin ang)
             bx (m/cos ang)]
         (Vec2. (m/* rad bx)
                (m/* rad by)))))))

(defn circular
  ([] {:type :random
       :config (fn [] {:angle (r/drand -45.0 45.0)
                      :seed (r/drand m/TWO_PI)})})
  ([^double amount {:keys [^double angle ^double seed]}]
   (let [ca (m/radians angle)]
     (fn [^Vec2 v]
       (let [aux (m/frac (m/* 43758.5453 (m/sin (m/+ (m/* 12.9898 (.x v))
                                                 (m/* 78.233 (.y v))
                                                 seed))))
             rnd (m/* (m/- (m/* 2.0 (m/+ (r/drand) aux)) 2.0) ca)
             rad (m/* amount (v/mag v))
             ang (m/+ rnd (v/heading v))
             by (m/sin ang)
             bx (m/cos ang)]
         (Vec2. (m/* rad bx)
                (m/* rad by)))))))

(defn circus
  ([] {:type :regular
       :config (fn [] {:scale (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double scale]}]
   (let [scale-1 (m// scale)]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             vn (v/div v r)
             r (if (m/<= r 1.0) (m/* r scale) (m/* r scale-1))]
         (v/mult vn (m/* amount r)))))))

(defn clifford
  "Clifford Pickover attractor"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (v/mult (Vec2. (m/+ (m/sin (m/* a (.y v))) (m/* c (m/cos (m/* a (.x v)))))
                    (m/+ (m/sin (m/* b (.x v))) (m/* d (m/cos (m/* b (.y v)))))) amount))))

(defn cloverleafwf
  "CloverLeafWF"
  ([] {:type :pattern
       :config (fn [] {:filled (r/randval 0.3 0.0 (r/randval 0.3 1.0 (r/drand)))})})
  ([^double amount {:keys [^double filled]}]
   (fn [v]
     (let [a (v/heading v)
           r (m/+ (m/sin (m/* 2.0 a))
                (m/* 0.25 (m/sin (m/* 6.0 a))))
           r (m/* amount (if (m/> filled (r/drand)) (m/* r (r/drand)) r))]
       (Vec2. (m/* r (m/sin a))
              (m/* r (m/cos a)))))))

(defn collideoscope
  ([] {:type :regular
       :config (fn [] {:a (r/drand)
                      :num (r/randval (u/sirand 1 10) (u/sdrand 0.2 10.0))})})
  ([^double amount {:keys [^double a ^long num]}]
   (let [kn-pi (m/* num m/M_1_PI)
         pi-kn (m// m/PI num)
         ka (m/* m/PI a)
         ka-kn (m// ka num)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (m/* amount (v/mag v))
             a (if (m/>= a 0.0)
                 (let [alt (long (m/* a kn-pi))]
                   (if (m/zero? (m/mod alt 2))
                     (m/+ (m/* alt pi-kn) (m/mod (m/+ ka-kn a) pi-kn))
                     (m/+ (m/* alt pi-kn) (m/mod (m/- a ka-kn) pi-kn))))
                 (let [alt (long (m/* -1.0 a kn-pi))]
                   (if-not (m/zero? (m/mod alt 2))
                     (m/- (m/+ (m/* alt pi-kn) (m/mod (m/- (m/- ka-kn) a) pi-kn)))
                     (m/- (m/+ (m/* alt pi-kn) (m/mod (m/- ka-kn a) pi-kn))))))]
         (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a))))))))

(defn conic
  "Conic"
  ([] {:type :random
       :config (fn [] {:eccentricity (r/drand -3.0 3.0)
                      :holes (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double eccentricity ^double holes]}]
   (fn [^Vec2 v]
     (let [magr (m// 1.0 (m/+ (v/mag v) m/EPSILON))
           ct (m/* (.x v) magr)
           r (m/* (m// (m/* (m/* amount (m/- (r/drand) holes)) eccentricity) (m/inc (m/* eccentricity ct))) magr)]
       (v/mult v r)))))

(defn corners
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)
                      :mult-x (u/sdrand 0.5 1.5)
                      :mult-y (u/sdrand 0.5 1.5)
                      :x-power (u/sdrand 0.25 2.0)
                      :y-power (u/sdrand 0.25 2.0)
                      :xy-power-add (r/randval 0.0 (r/drand -2.0 2.0))
                      :log-mode (r/brand)
                      :log-base (r/drand 1.1 20.0)})})
  ([^double amount {:keys [^double x ^double y ^double mult-x ^double mult-y
                           ^double x-power ^double y-power ^double xy-power-add
                           log-mode ^double log-base]}]
   (fn [^Vec2 v]
     (let [xs (m/sq (.x v))
           ys (m/sq (.y v))
           ex (if log-mode
                (m/- (m/pow (m/logb log-base (m/+ 3.0 (m/* xs mult-x)))
                          (m/+ x-power 2.25 xy-power-add)) 1.33)
                (m/* mult-x (m/pow xs (m/+ x-power xy-power-add))))
           ey (if log-mode
                (m/- (m/pow (m/logb log-base (m/+ 3.0 (m/* ys mult-y)))
                          (m/+ y-power 2.25 xy-power-add)) 1.33)
                (m/* mult-y (m/pow ys (m/+ y-power xy-power-add))))]
       (Vec2. (if (m/pos? (.x v)) (m/+ (m/* amount ex) x) (m/- (m/* amount (m/- ex)) x))
              (if (m/pos? (.y v)) (m/+ (m/* amount ey) y) (m/- (m/* amount (m/- ey)) y)))))))

(defn cos2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cossin (m/sin (m/* (.x v) x1))
           coscos (m/cos (m/* (.x v) x2))
           cossinh (m/sinh (m/* (.y v) y1))
           coscosh (m/cosh (m/* (.y v) y2))]
       (Vec2. (m/* amount coscos coscosh)
              (m/* -1.0 amount cossin cossinh))))))

(defn cos
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cossin (m/sin (.x v))
           coscos (m/cos (.x v))
           cossinh (m/sinh (.y v))
           coscosh (m/cosh (.y v))]
       (Vec2. (m/* amount coscos coscosh)
              (m/* -1.0 amount cossin cossinh))))))

(defn cosh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [coshsin (m/sin (m/* (.y v) y1))
           coshcos (m/cos (m/* (.y v) y2))
           coshsinh (m/sinh (m/* (.x v) x1))
           coshcosh (m/cosh (m/* (.x v) x2))]
       (Vec2. (m/* amount coshcos coshcosh)
              (m/* amount coshsin coshsinh))))))

(defn cosh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [coshsin (m/sin (.y v))
           coshcos (m/cos (.y v))
           coshsinh (m/sinh (.x v))
           coshcosh (m/cosh (.x v))]
       (Vec2. (m/* amount coshcos coshcosh)
              (m/* amount coshsin coshsinh))))))

(defn cosine
  "Cosine"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* m/PI (.x v))]
       (Vec2. (m/* amount (m/cos r) (m/cosh (.y v)))
              (m/- (m/* amount (m/sin r) (m/sinh (.y v)))))))))


(defn cot2bs
  "Cot2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cotsin (m/sin (m/* x1 (.x v)))
           cotcos (m/cos (m/* x2 (.x v)))
           cotsinh (m/sinh (m/* y1 (.y v)))
           cotcosh (m/cosh (m/* y2 (.y v)))
           cotden (m// (m/- cotcosh cotcos))]
       (Vec2. (m/* amount cotden cotsin)
              (m/* amount cotden -1.0 cotsinh))))))

(defn cot
  "Cot"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cotsin (m/sin (m/* 2.0 (.x v)))
           cotcos (m/cos (m/* 2.0 (.x v)))
           cotsinh (m/sinh (m/* 2.0 (.y v)))
           cotcosh (m/cosh (m/* 2.0 (.y v)))
           cotden (m// (m/- cotcosh cotcos))]
       (Vec2. (m/* amount cotden cotsin)
              (m/* amount cotden -1.0 cotsinh))))))

(defn coth2bs
  "Coth2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cothsin (m/sin (m/* y1 (.y v)))
           cothcos (m/cos (m/* y2 (.y v)))
           cothsinh (m/sinh (m/* x1 (.x v)))
           cothcosh (m/cosh (m/* x2 (.x v)))
           cothden (m// (m/- cothcosh cothcos))]
       (Vec2. (m/* amount cothden cothsinh)
              (m/* amount cothden cothsin))))))

(defn coth
  "Coth"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cothsin (m/sin (m/* 2.0 (.y v)))
           cothcos (m/cos (m/* 2.0 (.y v)))
           cothsinh (m/sinh (m/* 2.0 (.x v)))
           cothcosh (m/cosh (m/* 2.0 (.x v)))
           cothden (m// (m/- cothcosh cothcos))]
       (Vec2. (m/* amount cothden cothsinh)
              (m/* amount cothden cothsin))))))

(defn- crackle-position
  ^Vec2 [^Vec2 xy ^double z ^double s ^double d]
  (let [^Vec2 E (v/mult xy 2.5)
        z (m/+ (m/* z 2.5) 19.8)
        sv (if (m/zero? d)
             u/zerov
             (Vec2. (m/* d (r/simplex (.x E) (.y E) z))
                    (m/* d (r/simplex (m/+ (.x E) 30.2) (m/- (.y E) 12.1) z))))]
    (v/mult (v/add xy sv) s)))

(defn crackle
  "Crackle"
  ([] {:type :pattern
       :config (fn [] {:cellsize (u/sdrand 0.1 2.0)
                      :power (u/sdrand 0.01 2.0)
                      :distort (r/drand -3.0 3.0)
                      :scale (u/sdrand 0.1 1.5)
                      :z (r/drand -4.0 4.0)})})
  ([^double amount {:keys [^double cellsize ^double power ^double distort ^double scale ^double z]}]
   (let [s (m// cellsize 2.0)
         cache (memoize (fn ^Vec2 [xy] (crackle-position xy z s distort)))]
     (doseq [x (range -10 11)
             y (range -10 11)]
       (cache (Vec2. x y)))
     (fn [_]
       (let [blurr (m/* 2.0 (m/+ (m// (m/+ (r/drand) (r/drand)) 2.0)
                             (r/drand -0.125 0.125)))
             theta (r/drand m/TWO_PI)
             U (Vec2. (m/* blurr (m/sin theta))
                      (m/* blurr (m/cos theta)))
             Cv (v/floor (v/div U s))
             P (mapv #(cache (v/add Cv %)) u/offsets)
             q (u/closest P 9 U)
             Cv (v/add Cv (u/offsets q))
             P (mapv #(cache (v/add Cv %)) u/offsets)
             L (m/+ (u/voronoi P 9 4 U) 1.0e-100)
             P4 (P 4)
             Do (v/sub U P4)
             trgL (m/* (m/pow L power) scale)
             R (m// trgL L)]
         (v/mult (v/add (v/mult Do R) P4) amount))))))

(defn crob
  "Crob"
  ([] {:type :random
       :config (fn [] {:left (r/drand -0.1 -3.5)
                      :right (r/drand 0.1 3.5)
                      :top (r/drand -0.1 -3.5)
                      :bottom (r/drand 0.1 3.5)
                      :blur (r/brand 0.75)
                      :ratioblur (u/sdrand 0.01 0.3)
                      :directblur (r/drand 0.1 3.0)})})
  ([^double amount {:keys [^double left ^double right ^double top ^double bottom
                           blur ^double ratioblur ^double directblur]}]
   (let [[^double top ^double bottom] (if (m/> top bottom) [bottom top] [top bottom])
         [^double left ^double right] (if (m/> left right) [right left] [left right])
         [^double top ^double bottom] (if (m/== top bottom) [-1.0 1.0] [top bottom])
         [^double left ^double right] (if (m/== left right) [-1.0 1.0] [left right])
         directblur (m/max 0.0 directblur)
         xinterval (m/abs (m/- right left))
         yinterval (m/abs (m/- bottom top))
         xint2 (m// xinterval 2.0)
         yint2 (m// yinterval 2.0)
         minint2 (m/min xint2 yint2)
         x0 (m/- right xint2)
         y0 (m/+ top yint2)
         [^double x0c ^double y0c] (cond
                                     (m/> xint2 yint2) [(m/- right minint2) y0]
                                     (m/< xint2 yint2) [x0 (m/+ top minint2)]
                                     :else [x0 y0])
         setprob (m// yinterval (m/+ xinterval yinterval))
         setprobq (m/* 0.25 setprob)
         setprobh (m/* 0.5 setprob)
         setprobtq (m/* 0.75 setprob)
         setcompprob (m/- 1.0 setprob)
         setcompprobq (m/+ setprob (m/* 0.25 setcompprob))
         setcompprobh (m/+ setprob (m/* 0.5 setcompprob))
         setcompprobtq (m/+ setprob (m/* 0.75 setcompprob))
         mr (m/* minint2 ratioblur)
         [^double top-border ^double bottom-border
          ^double left-border ^double right-border] (if-not blur
                                                      [top bottom left right]
                                                      [(m/+ top mr) (m/- bottom mr)
                                                       (m/+ left mr) (m/- right mr)])]
     (fn [^Vec2 v]
       (if (or (m/< (.x v) left-border)
               (m/> (.x v) right-border)
               (m/< (.y v) top-border)
               (m/> (.y v) bottom-border))
         (if-not blur
           u/zerov
           (let [sectmp (r/drand)
                 [^double xtmp ^double ytmp] (if (m/< sectmp setprob)
                                               (loop []
                                                 (let [ytmp (m/+ top (r/drand yint2))
                                                       xtmp (m/- right (m/* (m/pow (r/drand) directblur)
                                                                        ratioblur minint2))]
                                                   (if (m/< (m// (m/- ytmp y0c)
                                                             (m/- xtmp x0c)) -1.0)
                                                     (recur)
                                                     [(if (m/< sectmp setprobh)
                                                        (m/- (m/+ left right) xtmp)
                                                        xtmp)
                                                      (if (and (m/> sectmp setprobq)
                                                               (m/< sectmp setprobtq))
                                                        (m/- (m/+ bottom top) ytmp)
                                                        ytmp)])))
                                               (loop []
                                                 (let [xtmp (m/- right (r/drand xint2))
                                                       ytmp (m/+ top (m/* (m/pow (r/drand) directblur)
                                                                      ratioblur minint2))
                                                       gradtmp (m// (m/- ytmp y0c)
                                                                  (m/- xtmp x0c))]
                                                   (if (and (m/<= gradtmp 0.0)
                                                            (m/> gradtmp -1.0))
                                                     (recur)
                                                     [(if (and (m/> sectmp setcompprobq)
                                                               (m/< sectmp setcompprobtq))
                                                        (m/- (m/+ left right) xtmp)
                                                        xtmp)
                                                      (if (m/> sectmp setcompprobh)
                                                        (m/- (m/+ bottom top) ytmp)
                                                        ytmp)]))))]
             (Vec2. (m/* amount xtmp) (m/* amount ytmp))))
         (v/mult v amount))))))

(defn cross
  "Cross"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [s (m/- (m/sq (.x v)) (m/sq (.y v)))
           r (m/* amount (m/sqrt (m// 1.0 (m/+ m/EPSILON (m/* s s)))))]
       (Vec2. (m/* (.x v) r) (m/* (.y v) r))))))

(defn csc2bs
  "Csc2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cscsin (m/sin (m/* x1 (.x v)))
           csccos (m/cos (m/* x2 (.x v)))
           cscsinh (m/sinh (m/* y1 (.y v)))
           csccosh (m/cosh (m/* y2 (.y v)))
           d (m/- (m/cosh (m/* 2.0 (.y v)))
                (m/cos (m/* 2.0 (.x v))))
           cscden (m// 2.0 d)]
       (Vec2. (m/* amount cscden cscsin csccosh)
              (m/* amount -1.0 cscden csccos cscsinh))))))

(defn csc
  "Csc"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cscsin (m/sin (.x v))
           csccos (m/cos (.x v))
           cscsinh (m/sinh (.y v))
           csccosh (m/cosh (.y v))
           d (m/- (m/cosh (m/* 2.0 (.y v)))
                (m/cos (m/* 2.0 (.x v))))
           cscden (m// 2.0 d)]
       (Vec2. (m/* amount cscden cscsin csccosh)
              (m/* amount -1.0 cscden csccos cscsinh))))))

(defn cscsquared
  "Csc squared"
  ([] {:type :regular
       :config (fn [] {:csc-div (u/sdrand 0.2 2.2)
                      :cos-div (u/sdrand 0.2 2.2)
                      :tan-div (u/sdrand 0.2 2.2)
                      :csc-pow (r/drand -3.0 1.0)
                      :pi-mult (r/drand 3.0)
                      :csc-add (r/drand)
                      :scale-y (u/sdrand 0.5 1.5)})})
  ([^double amount {:keys [^double csc-div ^double cos-div ^double tan-div
                           ^double csc-pow ^double pi-mult ^double csc-add
                           ^double scale-y]}]
   (fn [^Vec2 v]
     (let [csc (m// (m// csc-div (m/cos (m// (.x v) cos-div)))
                  (m/tan (m// (.x v) tan-div)))
           fx (m/* amount (m/+ (m/pow (m/+ (m/* csc csc) (m/* m/PI pi-mult)) csc-pow)
                           csc-add))]
       (Vec2. (m/* fx (.x v))
              (m/* fx (.y v) scale-y))))))

(defn csch2bs
  "Csch2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cschsin (m/sin (m/* y1 (.y v)))
           cschcos (m/cos (m/* y2 (.y v)))
           cschsinh (m/sinh (m/* x1 (.x v)))
           cschcosh (m/cosh (m/* x2 (.x v)))
           d (m/- (m/cosh (m/* 2.0 (.x v)))
                (m/cos (m/* 2.0 (.y v))))
           cschden (m// 2.0 d)]
       (Vec2. (m/* amount cschden cschsinh cschcos)
              (m/* amount -1.0 cschden cschcosh cschsin))))))

(defn csch
  "Csch"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cschsin (m/sin (.y v))
           cschcos (m/cos (.y v))
           cschsinh (m/sinh (.x v))
           cschcosh (m/cosh (.x v))
           d (m/- (m/cosh (m/* 2.0 (.x v)))
                (m/cos (m/* 2.0 (.y v))))
           cschden (m// 2.0 d)]
       (Vec2. (m/* amount cschden cschsinh cschcos)
              (m/* amount -1.0 cschden cschcosh cschsin))))))

(defn curl
  "Curl"
  ([] {:type :regular
       :config (fn [] {:c1 (r/drand -1.0 1.0)
                      :c2 (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double c1 ^double c2]}]
   (fn [^Vec2 v]
     (let [re (m/inc (m/+ (m/* c1 (.x v)) (m/* c2 (m/- (m/sq (.x v)) (m/sq (.y v))))))
           im (m/+ (m/* c1 (.y v)) (m/* c2 2.0 (.x v) (.y v)))
           r (m// amount (m/+ (m/sq re) (m/sq im)))]
       (Vec2. (m/* r (m/+ (m/* (.x v) re) (m/* (.y v) im)))
              (m/* r (m/- (m/* (.y v) re) (m/* (.x v) im))))))))

(defn curve
  ([] {:type :regular
       :config (fn [] {:xamp (u/sdrand 0.1 3.0)
                      :yamp (u/sdrand 0.1 3.0)
                      :xlength (r/drand 0.1 2.0)
                      :ylength (r/drand 0.1 2.0)})})
  ([^double amount {:keys [^double xamp ^double yamp ^double xlength ^double ylength]}]
   (let [pc-xlen (m// (m/- (m/max m/EPSILON (m/* xlength xlength))))
         pc-ylen (m// (m/- (m/max m/EPSILON (m/* ylength ylength))))]
     (fn [^Vec2 v]
       (Vec2. (m/* amount (m/+ (.x v) (m/* xamp (m/exp (m/* (.y v) (.y v) pc-xlen)))))
              (m/* amount (m/+ (.y v) (m/* yamp (m/exp (m/* (.x v) (.x v) pc-ylen))))))))))

(defn cylinder2
  "Cylinder"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (Vec2. (m// (.x v) (m/sqrt (m/inc (m/* (.x v) (.x v))))) (.y v)) amount))))

(defn cylinder
  "Cylinder"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (Vec2. (m/sin (.x v)) (.y v)) amount))))


;; 

(defn cayley
  "Cayley transform"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (if (m/== (.y v) -1.0)
       u/zerov
       (v/mult (c/div (c/add v c/I-)
                      (c/add v c/I)) amount)))))

(defn circle-inverse
  "Circle inverse"
  ([] {:type :regular
       :config (fn [] {:x0 (r/randval 0.3 0.0 (r/randval (r/irand -3 3) (r/drand -3.0 3.0)))
                      :y0 (r/randval 0.3 0.0 (r/randval (r/irand -3 3) (r/drand -3.0 3.0)))
                      :r (r/drand 0.01 3.0)})})
  ([^double amount {:keys [^double x0 ^double y0 ^double r]}]
   (let [o (Vec2. x0 y0)
         r2 (m/* r r)]
     (fn [^Vec2 v]
       (let [diff (v/sub v o)]
         (v/mult (v/add o (v/div (v/mult diff r2) (v/magsq diff))) amount))))))

