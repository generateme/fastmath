(ns fastmath.fields.s
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn stwin-jw
  "STwin JWildfire version"
  ([] {:type :regular
       :config (fn [] {:distort (r/drand -6.0 6.0)
                      :multiplier (u/sdrand 0.0001 3.0)
                      :offset-xy (r/drand -1.0 1.0)
                      :offset-x2 (r/drand -1.0 1.0)
                      :offset-y2 (r/drand -1.0 1.0)
                      :multiplier2 (u/sdrand 0.0001 2.0)
                      :multiplier3 (u/sdrand 0.0001 2.0)})})
  ([^double amount {:keys [^double distort ^double multiplier
                           ^double offset-xy ^double offset-x2 ^double offset-y2
                           ^double multiplier2 ^double multiplier3]}]
   (let [amultiplier (m/* amount multiplier)
         om-x2 (m/* offset-x2 multiplier2)
         om-y2 (m/* offset-y2 multiplier2)
         om-xy (m/* offset-xy multiplier3)]
     (fn [^Vec2 v]
       (let [x (m/* (.x v) amultiplier)
             y (m/* (.y v) amultiplier)
             x2 (m/+ (m/* x x) om-x2)
             y2 (m/+ (m/* y y) om-y2)
             x2+y2 (m/+ x2 y2)
             x2-y2 (m/- x2 y2)
             div (if (m/zero? x2+y2) 1.0 x2+y2)
             result (m// (m/* x2-y2 (m/sin (m/* m/TWO_PI distort (m/+ x y om-xy)))) div)]
         (Vec2. (m/+ (m/* amount (.x v)) result)
                (m/+ (m/* amount (.y v)) result)))))))


(defn stwin
  "STwin by Xyrus-02, http://timothy-vincent.deviantart.com/art/STwin-Plugin-136504836"
  ([] {:type :regular
       :config (fn [] {:distort (r/drand -6.0 6.0)
                      :multiplier (u/sdrand 0.001 3.0)})})
  ([^double amount {:keys [^double distort ^double multiplier]}]
   (let [amultiplier (m/* amount multiplier)]
     (fn [^Vec2 v]
       (let [x (m/* (.x v) amultiplier)
             y (m/* (.y v) amultiplier)
             x2 (m/* x x)
             y2 (m/* y y)
             x2+y2 (m/+ x2 y2)
             x2-y2 (m/- x2 y2)
             div (if (m/zero? x2+y2) 1.0 x2+y2)
             result (m// (m/* x2-y2 (m/sin (m/* m/TWO_PI distort (m/+ x y)))) div)]
         (Vec2. (m/+ (m/* amount (.x v)) result)
                (m/+ (m/* amount (.y v)) result)))))))


(defn sattractor
  ([] {:type :random
       :config (fn [] {:m (r/randval (r/irand 2 15) (r/drand 2.0 15.0))})})
  ([^double amount {:keys [^double m]}]
   (let [m (m/max 2.0 m)
         im (int m)
         a (mapv #(m/cos (m/* m/TWO_PI (m// ^long % m))) (range im))
         b (mapv #(m/sin (m/* m/TWO_PI (m// ^long % m))) (range im))
         hamount (m/* 0.5 amount)]
     (fn [^Vec2 v]
       (let [l (r/irand im)
             ^double al (a l)
             ^double bl (b l)]
         (-> (r/randval
              (Vec2. (m/+ (m/* 0.5 (.x v)) al)
                     (m/+ (m/* 0.5 (.y v)) bl))
              (let [xx (m/* (.x v) (.x v))]
                (Vec2. (m/+ (m/* (.x v) al)
                          (m/* (.y v) bl)
                          (m/* xx bl))
                       (m/+ (m/* (.y v) al)
                          (m/* (.x v) bl -1.0)
                          (m/* xx al)))))
             (v/mult hamount)))))))

(defn- scrambly-mx-randflip
  [^long idxmin ^long idxmax ^long seed v]
  (reduce (fn [v [^long j ^long prn]]
            (let [ridx (m/inc j)]
              (if (m/> idxmax ridx)
                (let [ridx (m/+ ridx (m/mod prn (m/- idxmax ridx)))]
                  (-> v
                      (assoc ridx (v j))
                      (assoc j (v ridx))))
                (reduced v)))) v (map vector
                                      (range idxmin Integer/MAX_VALUE)
                                      (next (iterate (fn [^long prn]
                                                       (let [prn (m/long-add (m/long-mult prn 1103515245) 12345)
                                                             prn (m/bit-or (m/bit-and prn 0xffff0000)
                                                                         (m/bit-and 0xff00
                                                                                  (m/bit-shift-left prn 8))
                                                                         (m/bit-and 0xff
                                                                                  (m/bit-shift-right prn 8)))
                                                             prn (if-not (m/zero? (m/bit-and prn 4))
                                                                   (m/- prn seed)
                                                                   (m/bit-xor prn seed))]
                                                         (if (m/neg? prn) (m/- prn) prn))) 1)))))

(defn scrambly
  ([] {:type :regular
       :config (fn [] {:l (r/irand 3 26)
                      :seed (r/randval (r/irand 51) (r/irand))
                      :byrows (r/brand)
                      :cellsize (u/sdrand 1.0 5.0)})})
  ([^double amount {:keys [^long l ^long seed byrows ^double cellsize]}]
   (let [LL (m/constrain (long (m/abs l)) 3 25)
         LL2 (m/long-mult LL LL)
         mx-rd (int-array (if (m/< seed 50)
                            (map (fn [^long j] (m/mod (m/+ seed j 1) LL2)) (range LL2))
                            (let [tmp (vec (range LL2))]
                              (if byrows
                                (scrambly-mx-randflip 0 LL2 seed tmp)
                                (reduce (fn [t ^long j]
                                          (scrambly-mx-randflip (m/long-mult j LL)
                                                                (m/long-mult (m/long-inc j) LL)
                                                                (m/long-add seed j)
                                                                t)) tmp (range LL))))))
         cella (m// cellsize LL)
         mzcella (m/* 0.5 cellsize)
         cellainv (m// cella)]
     (fn [^Vec2 v]
       (let [^Vec2 V (-> (v/shift v mzcella)
                         (v/mult cellainv))
             ^Vec2 I (v/floor V)]
         (if (or (m/neg? (.x I))
                 (m/neg? (.y I))
                 (m/>= (.x I) LL)
                 (m/>= (.y I) LL))
           (v/mult v amount)
           (let [V (v/sub V I)
                 swp (aget mx-rd (m/+ (.x I) (m/* LL (.y I))))
                 I (Vec2. (m/quot swp LL) (m/mod swp LL))]
             (-> (v/add V I)
                 (v/mult cella)
                 (v/shift (m/- mzcella))
                 (v/mult amount)))))))))


(defn scry2
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :star (u/sdrand 0.01 1.2)
                      :circle (u/sdrand 0.01 1.2)})})
  ([^double amount {:keys [^double sides ^double star ^double circle]}]
   (let [a (m// m/TWO_PI sides)
         sina (m/sin a)
         cosa (m/cos a)
         a (m/* m/M_PI_2 -1.0 star)
         sins (m/sin a)
         coss (m/cos a)
         a (m/* m/M_PI_2 circle)
         sinc (m/sin a)
         cosc (m/cos a)
         sides- (m/abs (m/dec sides))
         ramount (m// amount)]
     (fn [^Vec2 v]
       (let [xrt (.x v)
             yrt (.y v)
             r2 (m/+ (m/* xrt coss) (m/* (m/abs yrt) sins))
             circle (v/mag v)
             ^double r2 (loop [r2 r2 xrt xrt yrt yrt i (long 0)]
                          (if (m/< i sides-)
                            (let [nxrt (m/- (m/* xrt cosa) (m/* yrt sina))
                                  nyrt (m/+ (m/* xrt sina) (m/* yrt cosa))]
                              (recur (m/max r2 (m/+ (m/* nxrt coss) (m/* (m/abs nyrt) sins))) nxrt nyrt (m/inc i)))
                            (m/+ (m/* r2 cosc) (m/* circle sinc))))
             r1 r2
             r2 (if (m/> sides- 1) (m/* r2 r2) (m/* (m/abs r2) r2))]
         (v/mult v (m// (m/* r1 (m/+ r2 ramount)))))))))

(defn scry
  "Scry"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [t (v/magsq v)
           d (-> 1.0
                 (m// amount)
                 (m/+ t)
                 (m/* (m/sqrt t))
                 (m/+ m/EPSILON))
           r (m// 1.0 d)]
       (v/mult v r)))))

(defn sec2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [secsin (m/sin (m/* (.x v) x1))
           seccos (m/cos (m/* (.x v) x2))
           secsinh (m/sinh (m/* (.y v) y1))
           seccosh (m/cosh (m/* (.y v) y2))
           secden (m// (m/* amount 2.0) (m/+ (m/cos (m/* 2.0 (.x v)))
                                       (m/cosh (m/* 2.0 (.y v)))))]
       (Vec2. (m/* secden seccos seccosh)
              (m/* secden secsin secsinh))))))

(defn sec
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [secsin (m/sin (.x v))
           seccos (m/cos (.x v))
           secsinh (m/sinh (.y v))
           seccosh (m/cosh (.y v))
           secden (m// (m/* amount 2.0) (m/+ (m/cos (m/* 2.0 (.x v)))
                                       (m/cosh (m/* 2.0 (.y v)))))]
       (Vec2. (m/* secden seccos seccosh)
              (m/* secden secsin secsinh))))))

(defn secant2
  "Secant2"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* amount (v/mag v))
           cr (m/cos r)
           icr (m// 1.0 (if (m/zero? cr) m/EPSILON cr))
           ny (if (m/neg? cr)
                (m/* amount (m/inc icr))
                (m/* amount (m/dec icr)))]
       (Vec2. (m/* amount (.x v)) ny)))))

(defn sech2bs
  "Sech"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [sn (m/sin (m/* (.y v) y1))
           cn (m/cos (m/* (.y v) y2))
           snh (m/sinh (m/* (.x v) x1))
           cnh (m/cosh (m/* (.x v) x2))
           d (m/+ (m/cos (m/* 2.0 (.y v)))
                (m/cosh (m/* 2.0 (.x v))))
           d (if (m/zero? d) m/EPSILON d)
           den (m// 2.0 d)]
       (Vec2. (m/* amount den cn cnh)
              (m/* (m/- amount) den sn snh))))))

(defn sech
  "Sech"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sn (m/sin (.y v))
           cn (m/cos (.y v))
           snh (m/sinh (.x v))
           cnh (m/cosh (.x v))
           d (m/+ (m/cos (m/* 2.0 (.y v)))
                (m/cosh (m/* 2.0 (.x v))))
           d (if (m/zero? d) m/EPSILON d)
           den (m// 2.0 d)]
       (Vec2. (m/* amount den cn cnh)
              (m/* (m/- amount) den sn snh))))))

(defn separation
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.2 1.2)
                      :y (r/drand -1.2 1.2)
                      :xinside (r/drand -1.2 1.2)
                      :yinside (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double x ^double y ^double xinside ^double yinside]}]
   (let [sx2 (m/* x x)
         sy2 (m/* y y)]
     (fn [^Vec2 v]
       (let [xx (m/sqrt (m/+ (m/* (.x v) (.x v)) sx2))
             xi (m/* (.x v) xinside)
             yy (m/sqrt (m/+ (m/* (.y v) (.y v)) sy2))
             yi (m/* (.y v) yinside)]
         (v/mult (Vec2. (if (m/pos? (.x v)) (m/- xx xi) (m/- (m/+ xx xi)))
                        (if (m/pos? (.y v)) (m/- yy yi) (m/- (m/+ yy yi)))) amount))))))

(defn shift
  ([] {:type :regular
       :config (fn [] {:angle (r/drand 360.0)
                      :shift-x (r/drand -1.0 1.0)
                      :shift-y (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double angle ^double shift-x ^double shift-y]}]
   (let [ang (m/radians angle)
         sn (m/sin ang)
         cs (m/cos ang)]
     (fn [^Vec2 v]
       (v/mult (Vec2. (m/- (m/+ (.x v) (m/* cs shift-x)) (m/* sn shift-y))
                      (m/- (.y v) (m/* cs shift-y) (m/* sn shift-x))) amount)))))

(defn shredlin
  ([] {:type :regular
       :config (fn [] {:xdistance (u/sdrand 0.2 2.0)
                      :ydistance (u/sdrand 0.2 2.0)
                      :xwidth (u/sdrand 0.1 2.0)
                      :ywidth (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double xdistance ^double ydistance ^double xwidth ^double ywidth]}]
   (let [asxd (m/* amount xdistance)
         asyd (m/* amount ydistance)
         sxw- (m/- 1.0 xwidth)
         syw- (m/- 1.0 ywidth)]
     (fn [^Vec2 v]
       (let [xpos (if (m/neg? (.x v)) -0.5 0.5)
             ypos (if (m/neg? (.y v)) -0.5 0.5)
             xrng (m// (.x v) xdistance)
             ixrng (int xrng)
             yrng (m// (.y v) ydistance)
             iyrng (int yrng)]
         (Vec2. (m/* asxd (m/+ (m/* (m/- xrng ixrng) xwidth) ixrng (m/* xpos sxw-)))
                (m/* asyd (m/+ (m/* (m/- yrng iyrng) ywidth) iyrng (m/* ypos syw-)))))))))

(defn shredrad
  "ShreadRad"
  ([] {:type :regular
       :config (fn [] {:n (r/randval (u/sirand 1 9) (u/sdrand 0.0001 8.0))
                      :width (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double n ^double width]}]
   (let [sa (m// m/TWO_PI n)
         sa2 (m/* 0.5 sa)
         sa2sw (m/* sa2 width)
         pi3 (m/* 3.0 m/PI)]
     (fn [v]
       (let [ang (v/heading v)
             rad (v/mag v)
             xang (m// (m/+ ang pi3 sa2) sa)
             ixang (unchecked-int xang)
             zang (m/- (m/* sa (m/+ ixang (m/* width (m/- xang ixang)))) m/PI sa2sw)]
         (Vec2. (m/* amount rad (m/cos zang))
                (m/* amount rad (m/sin zang))))))))


(defn sigmoid
  ([] {:type :regular
       :config (fn [] {:shiftx (r/drand -2.0 2.0)
                      :shifty (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shiftx ^double shifty]}]
   (let [[ax sx] (if (m/< -1.0 shiftx 1.0)
                   (if (m/zero? shiftx)
                     [1.0 m/EPSILON]
                     [(if (m/neg? shiftx) -1.0 1.0) (m// shiftx)])
                   [1.0 shiftx])
         [ay sy] (if (m/< -1.0 shifty 1.0)
                   (if (m/zero? shifty)
                     [1.0 m/EPSILON]
                     [(if (m/neg? shifty) -1.0 1.0) (m// shifty)])
                   [1.0 shifty])
         s (v/mult (Vec2. sx sy) -5.0)
         a (Vec2. ax ay)
         vv (m/* 2.0 (m/abs amount))]
     (fn [^Vec2 v]
       (-> (v/ediv a (v/shift (v/exp (v/emult v s)) 1.0))
           (v/shift -0.5)
           (v/mult vv))))))

(defn sin2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [sinsin (m/sin (m/* (.x v) x1))
           sincos (m/cos (m/* (.x v) x2))
           sinsinh (m/sinh (m/* (.y v) y1))
           sincosh (m/cosh (m/* (.y v) y2))]
       (Vec2. (m/* amount sinsin sincosh)
              (m/* amount sincos sinsinh))))))

(defn sin
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sinsin (m/sin (.x v) )
           sincos (m/cos (.x v) )
           sinsinh (m/sinh (.y v))
           sincosh (m/cosh (.y v))]
       (Vec2. (m/* amount sinsin sincosh)
              (m/* amount sincos sinsinh))))))

(defn sineblur
  ([] {:type :pattern
       :config (fn [] {:power (r/randval 0.2 1.0 (r/drand 0.1 2.2))})})
  ([^double amount {:keys [^double power]}]
   (let [power (m/max 0.0 power)
         am (m// amount m/PI)]
     (fn [_]
       (let [ang (r/drand m/TWO_PI)
             r (m/* am (if (m/one? power)
                       (m/acos (r/drand -1.0 1.0))
                       (m/acos (m/dec (m/* 2.0 (m/exp (m/* (m/log (r/drand)) power)))))))]
         (Vec2. (m/* r (m/cos ang))
                (m/* r (m/sin ang))))))))

(defn sinh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [sinhsin (m/sin (m/* (.y v) y1))
           sinhcos (m/cos (m/* (.y v) y2))
           sinhsinh (m/sinh (m/* (.x v) x1))
           sinhcosh (m/cosh (m/* (.x v) x2))]
       (Vec2. (m/* amount sinhsinh sinhcos)
              (m/* amount sinhcosh sinhsin))))))

(defn sinh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sinhsin (m/sin (.y v))
           sinhcos (m/cos (.y v))
           sinhsinh (m/sinh (.x v))
           sinhcosh (m/cosh (.x v))]
       (Vec2. (m/* amount sinhsinh sinhcos)
              (m/* amount sinhcosh sinhsin))))))

(defn sintrange
  ([] {:type :regular
       :config (fn [] {:w (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double w]}]
   (fn [^Vec2 v]
     (let [xx (m/* (.x v) (.x v))
           yy (m/* (.y v) (.y v))
           wv (m/- w (m/* w (m/+ xx yy)))]
       (Vec2. (m/* amount (m/sin (.x v)) (m/+ xx wv))
              (m/* amount (m/sin (.y v)) (m/+ yy wv)))))))

(defn sinusgrid
  ([] {:type :regular
       :config (fn [] {:ampx (u/sdrand 0.1 0.9)
                      :ampy (u/sdrand 0.1 0.9)
                      :freqx (u/sdrand 0.1 5.0)
                      :freqy (u/sdrand 0.1 5.0)})})
  ([^double amount {:keys [^double ampx ^double ampy ^double freqx ^double freqy]}]
   (let [fx (m/* freqx m/TWO_PI)
         fy (m/* freqy m/TWO_PI)]
     (fn [^Vec2 v]
       (let [sx (m/- (m/cos (m/* (.x v) fx)))
             sy (m/- (m/cos (m/* (.y v) fy)))
             tx (m/lerp (.x v) sx ampx)
             ty (m/lerp (.y v) sy ampy)]
         (Vec2. (m/* amount tx) (m/* amount ty)))))))

(defn sinusoidal
  "Sinusoidal"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (m/sin (.x v))) (m/* amount (m/sin (.y v)))))))

(defn spherical
  "Spherical"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult v (m// amount (m/+ m/EPSILON (v/magsq v)))))))

(defn sphericaln
  ([] {:type :random
       :config (fn [] {:power (u/sdrand 1.0 8.0)
                      :dist (u/sdrand 0.1 3.0)})})
  ([^double amount {:keys [^double power ^double dist]}]
   (let [fpower (m// m/TWO_PI (m/floor power))]
     (fn [^Vec2 v]
       (let [R (m// amount (m/pow (v/mag v) dist))
             N (int (m/floor (r/drand power)))
             alpha (m/+ (v/heading v) (m/* N fpower))]
         (Vec2. (m/* R (m/cos alpha))
                (m/* R (m/sin alpha))))))))

(defn spiral
  "Spiral"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/+ m/EPSILON (v/mag v))
           revr (m// 1.0 r)
           sina (m/* (.x v) revr)
           cosa (m/* (.y v) revr)
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (m/* amount revr (m/+ cosa sinr))
              (m/* amount revr (m/- sina cosr)))))))

(defn spiralwing
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [c1 (m/* (.x v) (.x v))
           c2 (m/* (.y v) (.y v))
           d (m// amount (m/+ c1 c2 m/EPSILON))
           c2 (m/sin c2)]
       (Vec2. (m/* d (m/cos c1) c2)
              (m/* d (m/sin c1) c2))))))

(defn spligon
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :r (u/sdrand 0.5 1.5)
                      :i (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double sides ^double r ^double i]}]
   (let [th (m// sides m/TWO_PI)
         thi (m// th)
         j (m// (m/* m/PI i) (m/* -2.0 sides))]
     (fn [^Vec2 v]
       (let [t (m/+ j (m/* thi (m/floor (m/* (v/heading v) th))))
             dx (m/* r (m/sin t))
             dy (m/* r (m/cos t))]
         (Vec2.  (m/* amount (m/+ (.x v) dx))
                 (m/* amount (m/+ (.y v) dy))))))))

(defn splipticbs
  ([] {:type :random
       :config (fn [] {:x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (let [v- (m// amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (m/inc (v/magsq v))
             x2 (m/* 2.0 (.x v))
             xmax (m/* 0.5 (m/+ (m/sqrt (m/+ tmp x2))
                            (m/sqrt (m/- tmp x2))))
             a (m// (.x v) xmax)
             b (m/safe-sqrt (m/- 1.0 (m/* a a)))
             xx (m/+ (m/* v- (m/atan2 a b)) (if-not (m/neg? (.x v)) x (m/- x)))
             yy (m/+ (m/* v- (m/log (m/+ xmax (m/safe-sqrt (m/dec xmax))))) y)]
         (Vec2. xx (r/randval yy (m/- yy))))))))

(defn splitbrdr
  ([] {:type :random
       :config (fn [] {:x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)
                      :px (r/drand -0.5 0.5)
                      :py (r/drand -0.5 0.5)})})
  ([^double amount {:keys [^double x ^double y ^double px ^double py]}]
   (let [p (Vec2. px py)]
     (fn [^Vec2 v]
       (let [B (m/inc (m/* 0.25 (v/magsq v)))
             b (m// amount B)
             V (v/add (v/mult v b) (v/emult v p))
             ^Vec2 round (Vec2. (m/rint (.x v)) (m/rint (.y v)))
             ^Vec2 offset (v/mult (Vec2. (m/- (.x v) (.x round))
                                         (m/- (.y v) (.y round))) 0.5)
             roffset (v/add offset round)]
         (r/randval 0.25
                    (-> roffset
                        (v/mult amount)
                        (v/add V))
                    (if (m/>= (m/abs (.x offset))
                            (m/abs (.y offset)))
                      (if-not (m/neg? (.x offset))
                        (-> roffset
                            (v/add (Vec2. x (m// (m/* y (.y offset)) (.x offset))))
                            (v/mult amount)
                            (v/add V))
                        (-> roffset
                            (v/sub (Vec2. y (m// (m/* y (.y offset)) (.x offset))))
                            (v/mult amount)
                            (v/add V)))
                      (if-not (m/neg? (.y offset))
                        (-> roffset
                            (v/add (Vec2. (m// (m/* y (.x offset)) (.y offset)) y))
                            (v/mult amount)
                            (v/add V))
                        (-> roffset
                            (v/sub (Vec2. (m// (m/* x (.x offset)) (.y offset)) y))
                            (v/mult amount)
                            (v/add V))))))))))

(defn split
  "Split"
  ([] {:type :regular
       :config (fn [] {:xsplit (r/drand m/-TWO_PI m/TWO_PI)
                      :ysplit (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double xsplit ^double ysplit]}]
   (fn [^Vec2 v]
     (Vec2. (if (m/pos? (m/cos (m/* (.x v) xsplit)))
              (m/* amount (.y v))
              (m/- (m/* amount (.y v))))
            (if (m/pos? (m/cos (m/* (.y v) ysplit)))
              (m/* amount (.x v))
              (m/- (m/* amount (.x v))))))))

(defn splits
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (Vec2. (if (m/pos? (.x v))
              (m/* amount (m/+ (.x v) x))
              (m/* amount (m/- (.x v) x)))
            (if (m/pos? (.y v))
              (m/* amount (m/+ (.y v) y))
              (m/* amount (m/- (.y v) y)))))))


(defn sqrt-acosech
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acsch)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-acosh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acosh)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-acoth
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acoth)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-asech
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/asech)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-asinh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/asinh)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-atanh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (m/* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/atanh)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn square
  "Square"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (Vec2. (m/* amount (r/drand -0.5 0.5))
            (m/* amount (r/drand -0.5 0.5))))))

(defn squarize
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [s (v/mag v)
           a (v/heading v)
           a (if (m/neg? a) (m/+ a m/TWO_PI) a)
           p (m/* 4.0 s a m/M_1_PI)]
       (-> (cond
             (m/<= p s) (Vec2. s p)
             (m/<= p (m/* 3.0 s)) (Vec2. (m/- (m/* 2.0 s) p) s)
             (m/<= p (m/* 5.0 s)) (Vec2. (m/- s) (m/- (m/* 4.0 s) p))
             (m/<= p (m/* 7.0 s)) (Vec2. (m/- p (m/* 6.0 s)) (m/- s))
             :else (Vec2. s (m/- p (m/* 8.0 s))))
           (v/mult amount))))))

(defn squircular
  ([] {:type :regular})
  ([^double amount _]
   (let [a2 (m/* amount amount)]
     (fn [^Vec2 v]
       (let [u (.x v)
             u2 (m/* u u)
             v (.y v)
             v2 (m/* v v)
             r (m/+ u2 v2)
             rs (m/sqrt r)
             xs (m/sgn u)
             r (m/sqrt (m/- (m/* a2 r) (m/* 4.0 u2 v2)))
             r (m/sqrt (m/inc (m/- (m// u2 v2) (m// (m/* rs r) (m/* amount v2)))))
             r (m// r m/SQRT2)]
         (Vec2. (m/* xs r)
                (m/* (m// v u) r)))))))

(defn squirrel
  "Squirrel"
  ([] {:type :regular
       :config (fn [] {:a (r/drand m/EPSILON 4.0)
                      :b (r/drand m/EPSILON 4.0)})})
  ([^double amount {:keys [^double a ^double b]}]
   (fn [^Vec2 v]
     (let [u (m/sqrt (m/+ (m/* a (m/sq (.x v)))
                        (m/* b (m/sq (.y v)))))]
       (Vec2. (m/* amount (m/cos u) (m/tan (.x v)))
              (m/* amount (m/sin u) (m/tan (.y v))))))))

(defn squish
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double power]}]
   (let [inv-power (m// power)]
     (fn [^Vec2 v]
       (let [x (m/abs (.x v))
             y (m/abs (.y v))
             ^Vec2 sp (if (m/> x y)
                        (Vec2. x (if (m/pos? (.x v))
                                   (.y v)
                                   (m/- (m/* 4.0 x) (.y v))))
                        (Vec2. y (if (m/pos? (.y v))
                                   (m/- (m/* 2.0 y) (.x v))
                                   (m/+ (m/* 6.0 y) (.x v)))))
             s (.x sp)
             p (m/* inv-power (m/+ (.y sp) (m/* 8.0 s (m/floor (r/drand power)))))]
         (-> (cond
               (m/<= p s) (Vec2. s p)
               (m/<= p (m/* 3.0 s)) (Vec2. (m/- (m/* 2.0 s) p) s)
               (m/<= p (m/* 5.0 s)) (Vec2. (m/- s) (m/- (m/* 4.0 s) p))
               (m/<= p (m/* 7.0 s)) (Vec2. (m/- p (m/* 6.0 s)) (m/- s))
               :else (Vec2. s (m/- p (m/* 8.0 s))))
             (v/mult amount)))))))

(defn starblur
  ([] {:type :pattern
       :config (fn [] {:power (r/randval (u/sirand 2 10) (u/sdrand 0.1 10.0))
                      :range (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double power ^double range]}]
   (let [alpha (m// m/PI power)
         length (m/sqrt (m/inc (m/- (m/* range range) (m/* 2.0 range (m/cos alpha)))))
         alpha (m/asin (m// (m/* (m/sin alpha) range) length))
         calpha (m/cos alpha)
         salpha (m/sin alpha)
         power2 (m/* 2.0 power)
         ppower (m// m/TWO_PI power)]
     (fn [_]
       (let [f (r/drand power2)
             angle (m/trunc f)
             iangle (int angle)
             f (m/- f angle)
             x (m/* f length)
             z (m/sqrt (m/inc (m/- (m/* x x) (m/* 2.0 x calpha))))
             angle (m/- (if (m/even? iangle)
                        (m/+ (m/* ppower (m/quot iangle 2)) (m/asin (m// (m/* x salpha) z)))
                        (m/- (m/* ppower (m/quot iangle 2)) (m/asin (m// (m/* x salpha) z))))
                      m/M_PI_2)
             z (m/* z (m/sqrt (r/drand)))]
         (Vec2. (m/* amount z (m/cos angle))
                (m/* amount z (m/sin angle))))))))

(defn stripes
  ([] {:type :regular
       :config (fn [] {:space (r/drand -2.0 2.0)
                      :warp (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double space ^double warp]}]
   (let [space- (m/- 1.0 space)]
     (fn [^Vec2 v]
       (let [roundx (m/floor (m/+ (.x v) 0.5))
             offsetx (m/- (.x v) roundx)]
         (-> (Vec2. (m/+ (m/* offsetx space-) roundx)
                    (m/+ (.y v) (m/* offsetx offsetx warp)))
             (v/mult amount)))))))

(defn stripfit
  ([] {:type :regular
       :config (fn [] {:dx (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double dx]}]
   (let [dxp (m/* -0.5 dx)]
     (fn [^Vec2 v]
       (cond
         (m/> (.y v) 1.0) (let [fity (m/mod (m/inc (.y v)) 2.0)]
                          (Vec2. (m/+ (m/* amount (.x v)) (m/* dxp (m/inc (m/- (.y v) fity))))
                                 (m/* amount (m/dec fity))))
         (m/< (.y v) -1.0) (let [fity (m/mod (m/- 1.0 (.y v)) 2.0)]
                           (Vec2. (m/+ (m/* amount (.x v)) (m/* dxp (m/dec (m/+ (.y v) fity))))
                                  (m/* amount (m/- 1.0 fity))))
         :else (v/mult v amount))))))

(defn supershape
  "Supershape"
  ([] {:type :random
       :config (fn [] {:rnd (r/drand -1.0 1.0)
                      :m (r/drand m/TWO_PI)
                      :n1 (r/drand -5.0 5.0)
                      :n2 (r/drand -5.0 5.0)
                      :n3 (r/drand -5.0 5.0)
                      :holes (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double rnd ^double m ^double n1 ^double n2 ^double n3 ^double holes]}]
   (let [pm-4 (m// m 4.0)
         pneg1-n1 (m// -1.0 n1)]
     (fn [^Vec2 v]
       (let [theta (m/+ (m/* pm-4 (v/heading v)) m/M_PI_4)
             st (m/sin theta)
             ct (m/cos theta)
             t1 (m/pow (m/abs ct) n2)
             t2 (m/pow (m/abs st) n3)
             mag (v/mag v)
             r (m// (m/* (m/* amount (m/- (m/+ (r/drand rnd) (m/* (m/- 1.0 rnd) mag)) holes))
                     (m/pow (m/+ t1 t2) pneg1-n1)) mag)]
         (v/mult v r))))))

(defn svensson
  ([] {:type :regular
       :config (fn [] {:a (r/drand -3.0 3.0)
                      :b (r/drand -3.0 3.0)
                      :c (r/drand -3.0 3.0)
                      :d (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (let [ax (m/* a (.x v))
           by (m/* b (.y v))]
       (Vec2. (m/* amount (m/- (m/* d (m/sin ax)) (m/sin by)))
              (m/* amount (m/+ (m/* c (m/cos ax)) (m/cos by))))))))

(defn swirl3
  "Swirl3"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double shift]}]
   (fn [^Vec2 v]
     (let [rad (v/magsq v)
           ang (m/+ (v/heading v) (m/* shift (m/log rad)))
           s (m/sin ang)
           c (m/cos ang)
           arad (m/* amount rad)]
       (Vec2. (m/* arad c) (m/* arad s))))))

(defn swirl
  "Swirl"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/magsq v)
           c1 (m/sin r)
           c2 (m/cos r)]
       (Vec2. (m/* amount (m/- (m/* c1 (.x v)) (m/* c2 (.y v))))
              (m/* amount (m/+ (m/* c2 (.x v)) (m/* c1 (.y v)))))))))

(defmacro ^:private sym-transform
  [a b c d e f]
  (let [v (vary-meta (gensym "v") assoc :tag 'Vec2)]
    `(fn [~v] (Vec2. (m/+ (m/* ~a (.x ~v)) (m/* ~b (.y ~v)) ~c)
                    (m/+ (m/* ~d (.x ~v)) (m/* ~e (.y ~v)) ~f)))))

(defn- get-symband
  [{:keys [id ^double stepx ^double stepy]}]
  (let [sx  (m// stepx 2.0)
        sy  (m// stepy 2.0)
        -sx (m/- sx)
        -sy (m/- sy)]
    (case (int id)
      0 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 -sy)
         (sym-transform 1.0 0.0 sx 0.0 1.0 sy)]
      1 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 -sy)
         (sym-transform 1.0 0.0 sx 0.0 -1.0 (m/inc sy))]
      2 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 (m/- -sy 0.5))
         (sym-transform -1.0 0.0 (m/inc sx) 0.0 -1.0 (m/+ sy 0.5))]
      3 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 (m/- -sy 0.5))
         (sym-transform -1.0 0.0 (m/inc sx) 0.0 1.0 (m/+ sy 0.5))]
      4 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 (m/- -sy 0.5))
         (sym-transform 1.0 0.0 (m/dec sx) 0.0 -1.0 (m/+ sy 0.5))]
      5 [(sym-transform 1.0 0.0 (m/dec -sx) 0.0 1.0 (m/- -sy 0.5))
         (sym-transform -1.0 0.0 (m/inc sx) 0.0 1.0 (m/- sy 0.5))
         (sym-transform 1.0 0.0 (m/dec -sx) 0.0 -1.0 (m/+ -sy 0.5))
         (sym-transform -1.0 0.0 (m/inc sx) 0.0 -1.0 (m/+ sy 0.5))]
      6 [(sym-transform 1.0 0.0 (m/- -sx 2.0) 0.0 1.0 (m/- -sy 0.5))
         (sym-transform -1.0 0.0 (m/+ sx 2.0) 0.0 1.0 (m/- sy 0.5))
         (sym-transform 1.0 0.0 sx 0.0 -1.0 (m/+ sy 0.5))
         (sym-transform -1.0 0.0 -sx 0.0 -1.0 (m/+ -sy 0.5))])))

(defn symband
  ([] {:type :random
       :config (fn [] {:stepx (r/drand -2.0 2.0)
                      :stepy (r/drand -2.0 2.0)
                      :id (r/irand 7)})})
  ([^double amount opts]
   (let [transform (get-symband opts)]
     (fn [v]
       (let [f (rand-nth transform)]
         (v/mult (f v) amount))))))

(defn- get-symnet
  [{:keys [id ^double spacex ^double spacey ^double stepx ^double stepy]}]
  (let [spx (m// spacex 2.0)
        spy (m// spacey 2.0)
        -spx (m/- spx)
        -spy (m/- spy)
        stx (m// stepx 2.0)
        sty (m// stepy 2.0)
        -stx (m/- stx)
        -sty (m/- sty)]
    (case (int id)
      0[(sym-transform -1.0 0.0 -spx 0.0 -1.0 -spy)
        (sym-transform 0.0 1.0 -spx -1.0 0.0 -spy)
        (sym-transform 0.0 -1.0 spx 1.0 0.0 spy)
        (sym-transform 1.0 0.0 spx 0.0 1.0 spy)]
      1 [(sym-transform 1.0 0.0 (m/+ spacex stx) 0.0 1.0 (m/+ spacey sty))
         (sym-transform -1.0 0.0 (m/- stx spacex) 0.0 -1.0 (m/- sty spacey))
         (sym-transform 0.0 1.0 (m/+ spacex stx) -1.0 0.0 (m/- sty spacey))
         (sym-transform 0.0 -1.0 (m/- stx spacex) 1.0 0.0 (m/+ spacey sty))
         (sym-transform -1.0 0.0 (m/- -stx spacex) 0.0 1.0 (m/- spacey sty))
         (sym-transform 0.0 -1.0 (m/- -stx spacex) -1.0 0.0 (m/- -sty spacey))
         (sym-transform 0.0 1.0 (m/- spacex stx) 1.0 0.0 (m/- spacey sty))
         (sym-transform 1.0 0.0 (m/- spacex stx) 0.0 -1.0 (m/- spacey -sty))])))

;; TODO: add more
(defn symnet
  ([] {:type :random
       :config (fn [] {:space (r/drand -3.0 3.0) 
                      :spacex (r/drand -3.0 3.0)
                      :spacey (r/drand -3.0 3.0)
                      :stepx (r/drand -2.0 2.0)
                      :stepy (r/drand -2.0 2.0)
                      :id (r/irand 2)})})
  ([^double amount {:keys [^double space] :as opts}]
   (let [transform (get-symnet opts)
         vspace (Vec2. space space)]
     (fn [v]
       (let [f (rand-nth transform)]
         (v/mult (f (v/add vspace v)) amount))))))

;;

(defn secant
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* amount (v/mag v))
           cr (m/* amount (m/cos r))
           icr (m// 1.0 (if (m/zero? cr) m/EPSILON cr))]
       (Vec2. (m/* amount (.x v)) icr)))))

