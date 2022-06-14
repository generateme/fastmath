(ns fastmath.fields.s
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
   (let [amultiplier (* amount multiplier)
         om-x2 (* offset-x2 multiplier2)
         om-y2 (* offset-y2 multiplier2)
         om-xy (* offset-xy multiplier3)]
     (fn [^Vec2 v]
       (let [x (* (.x v) amultiplier)
             y (* (.y v) amultiplier)
             x2 (+ (* x x) om-x2)
             y2 (+ (* y y) om-y2)
             x2+y2 (+ x2 y2)
             x2-y2 (- x2 y2)
             div (if (zero? x2+y2) 1.0 x2+y2)
             result (/ (* x2-y2 (m/sin (* m/TWO_PI distort (+ x y om-xy)))) div)]
         (Vec2. (+ (* amount (.x v)) result)
                (+ (* amount (.y v)) result)))))))


(defn stwin
  "STwin by Xyrus-02, http://timothy-vincent.deviantart.com/art/STwin-Plugin-136504836"
  ([] {:type :regular
       :config (fn [] {:distort (r/drand -6.0 6.0)
                      :multiplier (u/sdrand 0.001 3.0)})})
  ([^double amount {:keys [^double distort ^double multiplier]}]
   (let [amultiplier (* amount multiplier)]
     (fn [^Vec2 v]
       (let [x (* (.x v) amultiplier)
             y (* (.y v) amultiplier)
             x2 (* x x)
             y2 (* y y)
             x2+y2 (+ x2 y2)
             x2-y2 (- x2 y2)
             div (if (zero? x2+y2) 1.0 x2+y2)
             result (/ (* x2-y2 (m/sin (* m/TWO_PI distort (+ x y)))) div)]
         (Vec2. (+ (* amount (.x v)) result)
                (+ (* amount (.y v)) result)))))))


(defn sattractor
  ([] {:type :random
       :config (fn [] {:m (r/randval (r/irand 2 15) (r/drand 2.0 15.0))})})
  ([^double amount {:keys [^double m]}]
   (let [m (max 2.0 m)
         im (int m)
         a (mapv #(m/cos (* m/TWO_PI (/ ^long % m))) (range im))
         b (mapv #(m/sin (* m/TWO_PI (/ ^long % m))) (range im))
         hamount (* 0.5 amount)]
     (fn [^Vec2 v]
       (let [l (r/irand im)
             ^double al (a l)
             ^double bl (b l)]
         (-> (r/randval
              (Vec2. (+ (* 0.5 (.x v)) al)
                     (+ (* 0.5 (.y v)) bl))
              (let [xx (* (.x v) (.x v))]
                (Vec2. (+ (* (.x v) al)
                          (* (.y v) bl)
                          (* xx bl))
                       (+ (* (.y v) al)
                          (* (.x v) bl -1.0)
                          (* xx al)))))
             (v/mult hamount)))))))

(defn- scrambly-mx-randflip
  [^long idxmin ^long idxmax ^long seed v]
  (reduce (fn [v [^long j ^long prn]]
            (let [ridx (inc j)]
              (if (> idxmax ridx)
                (let [ridx (+ ridx (mod prn (- idxmax ridx)))]
                  (-> v
                      (assoc ridx (v j))
                      (assoc j (v ridx))))
                (reduced v)))) v (map vector
                                      (range idxmin Integer/MAX_VALUE)
                                      (next (iterate (fn [^long prn]
                                                       (let [prn (+ (* prn 1103515245) 12345)
                                                             prn (bit-or (bit-and prn 0xffff0000)
                                                                         (bit-and 0xff00
                                                                                  (bit-shift-left prn 8))
                                                                         (bit-and 0xff
                                                                                  (bit-shift-right prn 8)))
                                                             prn (if-not (zero? (bit-and prn 4))
                                                                   (- prn seed)
                                                                   (bit-xor prn seed))]
                                                         (if (neg? prn) (- prn) prn))) 1)))))

(defn scrambly
  ([] {:type :regular
       :config (fn [] {:l (r/irand 3 26)
                      :seed (r/randval (r/irand 51) (r/irand))
                      :byrows (r/brand)
                      :cellsize (u/sdrand 1.0 5.0)})})
  ([^double amount {:keys [^long l ^long seed byrows ^double cellsize]}]
   (let [LL (m/constrain (long (m/abs l)) 3 25)
         LL2 (* LL LL)
         mx-rd (int-array (if (< seed 50)
                            (map (fn [^long j] (mod (+ seed j 1) LL2)) (range LL2))
                            (let [tmp (vec (range LL2))]
                              (if byrows
                                (scrambly-mx-randflip 0 LL2 seed tmp)
                                (reduce (fn [t ^long j]
                                          (scrambly-mx-randflip (* j LL)
                                                                (* (inc j) LL)
                                                                (+ seed j)
                                                                t)) tmp (range LL))))))
         cella (/ cellsize LL)
         mzcella (* 0.5 cellsize)
         cellainv (/ cella)]
     (fn [^Vec2 v]
       (let [^Vec2 V (-> (v/shift v mzcella)
                         (v/mult cellainv))
             ^Vec2 I (v/floor V)]
         (if (or (neg? (.x I))
                 (neg? (.y I))
                 (>= (.x I) LL)
                 (>= (.y I) LL))
           (v/mult v amount)
           (let [V (v/sub V I)
                 swp (aget mx-rd (+ (.x I) (* LL (.y I))))
                 I (Vec2. (quot swp LL) (mod swp LL))]
             (-> (v/add V I)
                 (v/mult cella)
                 (v/shift (- mzcella))
                 (v/mult amount)))))))))


(defn scry2
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :star (u/sdrand 0.01 1.2)
                      :circle (u/sdrand 0.01 1.2)})})
  ([^double amount {:keys [^double sides ^double star ^double circle]}]
   (let [a (/ m/TWO_PI sides)
         sina (m/sin a)
         cosa (m/cos a)
         a (* m/M_PI_2 -1.0 star)
         sins (m/sin a)
         coss (m/cos a)
         a (* m/M_PI_2 circle)
         sinc (m/sin a)
         cosc (m/cos a)
         sides- (m/abs (dec sides))
         ramount (/ amount)]
     (fn [^Vec2 v]
       (let [xrt (.x v)
             yrt (.y v)
             r2 (+ (* xrt coss) (* (m/abs yrt) sins))
             circle (v/mag v)
             ^double r2 (loop [r2 r2 xrt xrt yrt yrt i (long 0)]
                          (if (< i sides-)
                            (let [nxrt (- (* xrt cosa) (* yrt sina))
                                  nyrt (+ (* xrt sina) (* yrt cosa))]
                              (recur (max r2 (+ (* nxrt coss) (* (m/abs nyrt) sins))) nxrt nyrt (inc i)))
                            (+ (* r2 cosc) (* circle sinc))))
             r1 r2
             r2 (if (> sides- 1) (* r2 r2) (* (m/abs r2) r2))]
         (v/mult v (/ (* r1 (+ r2 ramount)))))))))

(defn scry
  "Scry"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [t (v/magsq v)
           d (-> 1.0
                 (/ amount)
                 (+ t)
                 (* (m/sqrt t))
                 (+ m/EPSILON))
           r (/ 1.0 d)]
       (v/mult v r)))))

(defn sec2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [secsin (m/sin (* (.x v) x1))
           seccos (m/cos (* (.x v) x2))
           secsinh (m/sinh (* (.y v) y1))
           seccosh (m/cosh (* (.y v) y2))
           secden (/ (* amount 2.0) (+ (m/cos (* 2.0 (.x v)))
                                       (m/cosh (* 2.0 (.y v)))))]
       (Vec2. (* secden seccos seccosh)
              (* secden secsin secsinh))))))

(defn sec
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [secsin (m/sin (.x v))
           seccos (m/cos (.x v))
           secsinh (m/sinh (.y v))
           seccosh (m/cosh (.y v))
           secden (/ (* amount 2.0) (+ (m/cos (* 2.0 (.x v)))
                                       (m/cosh (* 2.0 (.y v)))))]
       (Vec2. (* secden seccos seccosh)
              (* secden secsin secsinh))))))

(defn secant2
  "Secant2"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (* amount (v/mag v))
           cr (m/cos r)
           icr (/ 1.0 (if (zero? cr) m/EPSILON cr))
           ny (if (neg? cr)
                (* amount (inc icr))
                (* amount (dec icr)))]
       (Vec2. (* amount (.x v)) ny)))))

(defn sech2bs
  "Sech"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [sn (m/sin (* (.y v) y1))
           cn (m/cos (* (.y v) y2))
           snh (m/sinh (* (.x v) x1))
           cnh (m/cosh (* (.x v) x2))
           d (+ (m/cos (* 2.0 (.y v)))
                (m/cosh (* 2.0 (.x v))))
           d (if (zero? d) m/EPSILON d)
           den (/ 2.0 d)]
       (Vec2. (* amount den cn cnh)
              (* (- amount) den sn snh))))))

(defn sech
  "Sech"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sn (m/sin (.y v))
           cn (m/cos (.y v))
           snh (m/sinh (.x v))
           cnh (m/cosh (.x v))
           d (+ (m/cos (* 2.0 (.y v)))
                (m/cosh (* 2.0 (.x v))))
           d (if (zero? d) m/EPSILON d)
           den (/ 2.0 d)]
       (Vec2. (* amount den cn cnh)
              (* (- amount) den sn snh))))))

(defn separation
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.2 1.2)
                      :y (r/drand -1.2 1.2)
                      :xinside (r/drand -1.2 1.2)
                      :yinside (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double x ^double y ^double xinside ^double yinside]}]
   (let [sx2 (* x x)
         sy2 (* y y)]
     (fn [^Vec2 v]
       (let [xx (m/sqrt (+ (* (.x v) (.x v)) sx2))
             xi (* (.x v) xinside)
             yy (m/sqrt (+ (* (.y v) (.y v)) sy2))
             yi (* (.y v) yinside)]
         (v/mult (Vec2. (if (pos? (.x v)) (- xx xi) (- (+ xx xi)))
                        (if (pos? (.y v)) (- yy yi) (- (+ yy yi)))) amount))))))

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
       (v/mult (Vec2. (- (+ (.x v) (* cs shift-x)) (* sn shift-y))
                      (- (.y v) (* cs shift-y) (* sn shift-x))) amount)))))

(defn shredlin
  ([] {:type :regular
       :config (fn [] {:xdistance (u/sdrand 0.2 2.0)
                      :ydistance (u/sdrand 0.2 2.0)
                      :xwidth (u/sdrand 0.1 2.0)
                      :ywidth (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double xdistance ^double ydistance ^double xwidth ^double ywidth]}]
   (let [asxd (* amount xdistance)
         asyd (* amount ydistance)
         sxw- (- 1.0 xwidth)
         syw- (- 1.0 ywidth)]
     (fn [^Vec2 v]
       (let [xpos (if (neg? (.x v)) -0.5 0.5)
             ypos (if (neg? (.y v)) -0.5 0.5)
             xrng (/ (.x v) xdistance)
             ixrng (int xrng)
             yrng (/ (.y v) ydistance)
             iyrng (int yrng)]
         (Vec2. (* asxd (+ (* (- xrng ixrng) xwidth) ixrng (* xpos sxw-)))
                (* asyd (+ (* (- yrng iyrng) ywidth) iyrng (* ypos syw-)))))))))

(defn shredrad
  "ShreadRad"
  ([] {:type :regular
       :config (fn [] {:n (r/randval (u/sirand 1 9) (u/sdrand 0.0001 8.0))
                      :width (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double n ^double width]}]
   (let [sa (/ m/TWO_PI n)
         sa2 (* 0.5 sa)
         sa2sw (* sa2 width)
         pi3 (* 3.0 m/PI)]
     (fn [v]
       (let [ang (v/heading v)
             rad (v/mag v)
             xang (/ (+ ang pi3 sa2) sa)
             ixang (unchecked-int xang)
             zang (- (* sa (+ ixang (* width (- xang ixang)))) m/PI sa2sw)]
         (Vec2. (* amount rad (m/cos zang))
                (* amount rad (m/sin zang))))))))


(defn sigmoid
  ([] {:type :regular
       :config (fn [] {:shiftx (r/drand -2.0 2.0)
                      :shifty (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shiftx ^double shifty]}]
   (let [[ax sx] (if (< -1.0 shiftx 1.0)
                   (if (zero? shiftx)
                     [1.0 m/EPSILON]
                     [(if (neg? shiftx) -1.0 1.0) (/ shiftx)])
                   [1.0 shiftx])
         [ay sy] (if (< -1.0 shifty 1.0)
                   (if (zero? shifty)
                     [1.0 m/EPSILON]
                     [(if (neg? shifty) -1.0 1.0) (/ shifty)])
                   [1.0 shifty])
         s (v/mult (Vec2. sx sy) -5.0)
         a (Vec2. ax ay)
         vv (* 2.0 (m/abs amount))]
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
     (let [sinsin (m/sin (* (.x v) x1))
           sincos (m/cos (* (.x v) x2))
           sinsinh (m/sinh (* (.y v) y1))
           sincosh (m/cosh (* (.y v) y2))]
       (Vec2. (* amount sinsin sincosh)
              (* amount sincos sinsinh))))))

(defn sin
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sinsin (m/sin (.x v) )
           sincos (m/cos (.x v) )
           sinsinh (m/sinh (.y v))
           sincosh (m/cosh (.y v))]
       (Vec2. (* amount sinsin sincosh)
              (* amount sincos sinsinh))))))

(defn sineblur
  ([] {:type :pattern
       :config (fn [] {:power (r/randval 0.2 1.0 (r/drand 0.1 2.2))})})
  ([^double amount {:keys [^double power]}]
   (let [power (max 0.0 power)
         am (/ amount m/PI)]
     (fn [_]
       (let [ang (r/drand m/TWO_PI)
             r (* am (if (m/one? power)
                       (m/acos (r/drand -1.0 1.0))
                       (m/acos (dec (* 2.0 (m/exp (* (m/log (r/drand)) power)))))))]
         (Vec2. (* r (m/cos ang))
                (* r (m/sin ang))))))))

(defn sinh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [sinhsin (m/sin (* (.y v) y1))
           sinhcos (m/cos (* (.y v) y2))
           sinhsinh (m/sinh (* (.x v) x1))
           sinhcosh (m/cosh (* (.x v) x2))]
       (Vec2. (* amount sinhsinh sinhcos)
              (* amount sinhcosh sinhsin))))))

(defn sinh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [sinhsin (m/sin (.y v))
           sinhcos (m/cos (.y v))
           sinhsinh (m/sinh (.x v))
           sinhcosh (m/cosh (.x v))]
       (Vec2. (* amount sinhsinh sinhcos)
              (* amount sinhcosh sinhsin))))))

(defn sintrange
  ([] {:type :regular
       :config (fn [] {:w (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double w]}]
   (fn [^Vec2 v]
     (let [xx (* (.x v) (.x v))
           yy (* (.y v) (.y v))
           wv (- w (* w (+ xx yy)))]
       (Vec2. (* amount (m/sin (.x v)) (+ xx wv))
              (* amount (m/sin (.y v)) (+ yy wv)))))))

(defn sinusgrid
  ([] {:type :regular
       :config (fn [] {:ampx (u/sdrand 0.1 0.9)
                      :ampy (u/sdrand 0.1 0.9)
                      :freqx (u/sdrand 0.1 5.0)
                      :freqy (u/sdrand 0.1 5.0)})})
  ([^double amount {:keys [^double ampx ^double ampy ^double freqx ^double freqy]}]
   (let [fx (* freqx m/TWO_PI)
         fy (* freqy m/TWO_PI)]
     (fn [^Vec2 v]
       (let [sx (- (m/cos (* (.x v) fx)))
             sy (- (m/cos (* (.y v) fy)))
             tx (m/lerp (.x v) sx ampx)
             ty (m/lerp (.y v) sy ampy)]
         (Vec2. (* amount tx) (* amount ty)))))))

(defn sinusoidal
  "Sinusoidal"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (m/sin (.x v))) (* amount (m/sin (.y v)))))))

(defn spherical
  "Spherical"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult v (/ amount (+ m/EPSILON (v/magsq v)))))))

(defn sphericaln
  ([] {:type :random
       :config (fn [] {:power (u/sdrand 1.0 8.0)
                      :dist (u/sdrand 0.1 3.0)})})
  ([^double amount {:keys [^double power ^double dist]}]
   (let [fpower (/ m/TWO_PI (m/floor power))]
     (fn [^Vec2 v]
       (let [R (/ amount (m/pow (v/mag v) dist))
             N (int (m/floor (r/drand power)))
             alpha (+ (v/heading v) (* N fpower))]
         (Vec2. (* R (m/cos alpha))
                (* R (m/sin alpha))))))))

(defn spiral
  "Spiral"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (+ m/EPSILON (v/mag v))
           revr (/ 1.0 r)
           sina (* (.x v) revr)
           cosa (* (.y v) revr)
           sinr (m/sin r)
           cosr (m/cos r)]
       (Vec2. (* amount revr (+ cosa sinr))
              (* amount revr (- sina cosr)))))))

(defn spiralwing
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [c1 (* (.x v) (.x v))
           c2 (* (.y v) (.y v))
           d (/ amount (+ c1 c2 m/EPSILON))
           c2 (m/sin c2)]
       (Vec2. (* d (m/cos c1) c2)
              (* d (m/sin c1) c2))))))

(defn spligon
  ([] {:type :regular
       :config (fn [] {:sides (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :r (u/sdrand 0.5 1.5)
                      :i (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double sides ^double r ^double i]}]
   (let [th (/ sides m/TWO_PI)
         thi (/ th)
         j (/ (* m/PI i) (* -2.0 sides))]
     (fn [^Vec2 v]
       (let [t (+ j (* thi (m/floor (* (v/heading v) th))))
             dx (* r (m/sin t))
             dy (* r (m/cos t))]
         (Vec2.  (* amount (+ (.x v) dx))
                 (* amount (+ (.y v) dy))))))))

(defn splipticbs
  ([] {:type :random
       :config (fn [] {:x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double x ^double y]}]
   (let [v- (/ amount m/HALF_PI)]
     (fn [^Vec2 v]
       (let [tmp (inc (v/magsq v))
             x2 (* 2.0 (.x v))
             xmax (* 0.5 (+ (m/sqrt (+ tmp x2))
                            (m/sqrt (- tmp x2))))
             a (/ (.x v) xmax)
             b (m/safe-sqrt (- 1.0 (* a a)))
             xx (+ (* v- (m/atan2 a b)) (if-not (neg? (.x v)) x (- x)))
             yy (+ (* v- (m/log (+ xmax (m/safe-sqrt (dec xmax))))) y)]
         (Vec2. xx (r/randval yy (- yy))))))))

(defn splitbrdr
  ([] {:type :random
       :config (fn [] {:x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)
                      :px (r/drand -0.5 0.5)
                      :py (r/drand -0.5 0.5)})})
  ([^double amount {:keys [^double x ^double y ^double px ^double py]}]
   (let [p (Vec2. px py)]
     (fn [^Vec2 v]
       (let [B (inc (* 0.25 (v/magsq v)))
             b (/ amount B)
             V (v/add (v/mult v b) (v/emult v p))
             ^Vec2 round (Vec2. (m/rint (.x v)) (m/rint (.y v)))
             ^Vec2 offset (v/mult (Vec2. (- (.x v) (.x round))
                                         (- (.y v) (.y round))) 0.5)
             roffset (v/add offset round)]
         (r/randval 0.25
                    (-> roffset
                        (v/mult amount)
                        (v/add V))
                    (if (>= (m/abs (.x offset))
                            (m/abs (.y offset)))
                      (if-not (neg? (.x offset))
                        (-> roffset
                            (v/add (Vec2. x (/ (* y (.y offset)) (.x offset))))
                            (v/mult amount)
                            (v/add V))
                        (-> roffset
                            (v/sub (Vec2. y (/ (* y (.y offset)) (.x offset))))
                            (v/mult amount)
                            (v/add V)))
                      (if-not (neg? (.y offset))
                        (-> roffset
                            (v/add (Vec2. (/ (* y (.x offset)) (.y offset)) y))
                            (v/mult amount)
                            (v/add V))
                        (-> roffset
                            (v/sub (Vec2. (/ (* x (.x offset)) (.y offset)) y))
                            (v/mult amount)
                            (v/add V))))))))))

(defn split
  "Split"
  ([] {:type :regular
       :config (fn [] {:xsplit (r/drand m/-TWO_PI m/TWO_PI)
                      :ysplit (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double xsplit ^double ysplit]}]
   (fn [^Vec2 v]
     (Vec2. (if (pos? (m/cos (* (.x v) xsplit)))
              (* amount (.y v))
              (- (* amount (.y v))))
            (if (pos? (m/cos (* (.y v) ysplit)))
              (* amount (.x v))
              (- (* amount (.x v))))))))

(defn splits
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double x ^double y]}]
   (fn [^Vec2 v]
     (Vec2. (if (pos? (.x v))
              (* amount (+ (.x v) x))
              (* amount (- (.x v) x)))
            (if (pos? (.y v))
              (* amount (+ (.y v) y))
              (* amount (- (.y v) y)))))))


(defn sqrt-acosech
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acosech)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-acosh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acosh)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-acoth
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/acoth)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-asech
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/asech)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-asinh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
     (fn [z]
       (let [z (-> (c/sqrt z)
                   (c/asinh)
                   (c/scale atwopi))]
         (r/randval z (c/neg z)))))))

(defn sqrt-atanh
  ([] {:type :random})
  ([^double amount _]
   (let [atwopi (* amount m/M_2_PI)]
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
     (Vec2. (* amount (r/drand -0.5 0.5))
            (* amount (r/drand -0.5 0.5))))))

(defn squarize
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [s (v/mag v)
           a (v/heading v)
           a (if (neg? a) (+ a m/TWO_PI) a)
           p (* 4.0 s a m/M_1_PI)]
       (-> (cond
             (<= p s) (Vec2. s p)
             (<= p (* 3.0 s)) (Vec2. (- (* 2.0 s) p) s)
             (<= p (* 5.0 s)) (Vec2. (- s) (- (* 4.0 s) p))
             (<= p (* 7.0 s)) (Vec2. (- p (* 6.0 s)) (- s))
             :else (Vec2. s (- p (* 8.0 s))))
           (v/mult amount))))))

(defn squircular
  ([] {:type :regular})
  ([^double amount _]
   (let [a2 (* amount amount)]
     (fn [^Vec2 v]
       (let [u (.x v)
             u2 (* u u)
             v (.y v)
             v2 (* v v)
             r (+ u2 v2)
             rs (m/sqrt r)
             xs (m/sgn u)
             r (m/sqrt (- (* a2 r) (* 4.0 u2 v2)))
             r (m/sqrt (inc (- (/ u2 v2) (/ (* rs r) (* amount v2)))))
             r (/ r m/SQRT2)]
         (Vec2. (* xs r)
                (* (/ v u) r)))))))

(defn squirrel
  "Squirrel"
  ([] {:type :regular
       :config (fn [] {:a (r/drand m/EPSILON 4.0)
                      :b (r/drand m/EPSILON 4.0)})})
  ([^double amount {:keys [^double a ^double b]}]
   (fn [^Vec2 v]
     (let [u (m/sqrt (+ (* a (m/sq (.x v)))
                        (* b (m/sq (.y v)))))]
       (Vec2. (* amount (m/cos u) (m/tan (.x v)))
              (* amount (m/sin u) (m/tan (.y v))))))))

(defn squish
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double power]}]
   (let [inv-power (/ power)]
     (fn [^Vec2 v]
       (let [x (m/abs (.x v))
             y (m/abs (.y v))
             ^Vec2 sp (if (> x y)
                        (Vec2. x (if (pos? (.x v))
                                   (.y v)
                                   (- (* 4.0 x) (.y v))))
                        (Vec2. y (if (pos? (.y v))
                                   (- (* 2.0 y) (.x v))
                                   (+ (* 6.0 y) (.x v)))))
             s (.x sp)
             p (* inv-power (+ (.y sp) (* 8.0 s (m/floor (r/drand power)))))]
         (-> (cond
               (<= p s) (Vec2. s p)
               (<= p (* 3.0 s)) (Vec2. (- (* 2.0 s) p) s)
               (<= p (* 5.0 s)) (Vec2. (- s) (- (* 4.0 s) p))
               (<= p (* 7.0 s)) (Vec2. (- p (* 6.0 s)) (- s))
               :else (Vec2. s (- p (* 8.0 s))))
             (v/mult amount)))))))

(defn starblur
  ([] {:type :pattern
       :config (fn [] {:power (r/randval (u/sirand 2 10) (u/sdrand 0.1 10.0))
                      :range (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double power ^double range]}]
   (let [alpha (/ m/PI power)
         length (m/sqrt (inc (- (* range range) (* 2.0 range (m/cos alpha)))))
         alpha (m/asin (/ (* (m/sin alpha) range) length))
         calpha (m/cos alpha)
         salpha (m/sin alpha)
         power2 (* 2.0 power)
         ppower (/ m/TWO_PI power)]
     (fn [_]
       (let [f (r/drand power2)
             angle (m/trunc f)
             iangle (int angle)
             f (- f angle)
             x (* f length)
             z (m/sqrt (inc (- (* x x) (* 2.0 x calpha))))
             angle (- (if (even? iangle)
                        (+ (* ppower (quot iangle 2)) (m/asin (/ (* x salpha) z)))
                        (- (* ppower (quot iangle 2)) (m/asin (/ (* x salpha) z))))
                      m/M_PI_2)
             z (* z (m/sqrt (r/drand)))]
         (Vec2. (* amount z (m/cos angle))
                (* amount z (m/sin angle))))))))

(defn stripes
  ([] {:type :regular
       :config (fn [] {:space (r/drand -2.0 2.0)
                      :warp (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double space ^double warp]}]
   (let [space- (- 1.0 space)]
     (fn [^Vec2 v]
       (let [roundx (m/floor (+ (.x v) 0.5))
             offsetx (- (.x v) roundx)]
         (-> (Vec2. (+ (* offsetx space-) roundx)
                    (+ (.y v) (* offsetx offsetx warp)))
             (v/mult amount)))))))

(defn stripfit
  ([] {:type :regular
       :config (fn [] {:dx (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double dx]}]
   (let [dxp (* -0.5 dx)]
     (fn [^Vec2 v]
       (cond
         (> (.y v) 1.0) (let [fity (mod (inc (.y v)) 2.0)]
                          (Vec2. (+ (* amount (.x v)) (* dxp (inc (- (.y v) fity))))
                                 (* amount (dec fity))))
         (< (.y v) -1.0) (let [fity (mod (- 1.0 (.y v)) 2.0)]
                           (Vec2. (+ (* amount (.x v)) (* dxp (dec (+ (.y v) fity))))
                                  (* amount (- 1.0 fity))))
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
   (let [pm-4 (/ m 4.0)
         pneg1-n1 (/ -1.0 n1)]
     (fn [^Vec2 v]
       (let [theta (+ (* pm-4 (v/heading v)) m/M_PI_4)
             st (m/sin theta)
             ct (m/cos theta)
             t1 (m/pow (m/abs ct) n2)
             t2 (m/pow (m/abs st) n3)
             mag (v/mag v)
             r (/ (* (* amount (- (+ (r/drand rnd) (* (- 1.0 rnd) mag)) holes))
                     (m/pow (+ t1 t2) pneg1-n1)) mag)]
         (v/mult v r))))))

(defn svensson
  ([] {:type :regular
       :config (fn [] {:a (r/drand -3.0 3.0)
                      :b (r/drand -3.0 3.0)
                      :c (r/drand -3.0 3.0)
                      :d (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (let [ax (* a (.x v))
           by (* b (.y v))]
       (Vec2. (* amount (- (* d (m/sin ax)) (m/sin by)))
              (* amount (+ (* c (m/cos ax)) (m/cos by))))))))

(defn swirl3
  "Swirl3"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double shift]}]
   (fn [^Vec2 v]
     (let [rad (v/magsq v)
           ang (+ (v/heading v) (* shift (m/log rad)))
           s (m/sin ang)
           c (m/cos ang)
           arad (* amount rad)]
       (Vec2. (* arad c) (* arad s))))))

(defn swirl
  "Swirl"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/magsq v)
           c1 (m/sin r)
           c2 (m/cos r)]
       (Vec2. (* amount (- (* c1 (.x v)) (* c2 (.y v))))
              (* amount (+ (* c2 (.x v)) (* c1 (.y v)))))))))

(defmacro ^:private sym-transform
  [a b c d e f]
  (let [v (vary-meta (gensym "v") assoc :tag 'Vec2)]
    `(fn [~v] (Vec2. (+ (* ~a (.x ~v)) (* ~b (.y ~v)) ~c)
                    (+ (* ~d (.x ~v)) (* ~e (.y ~v)) ~f)))))

(defn- get-symband
  [{:keys [id ^double stepx ^double stepy]}]
  (let [sx  (/ stepx 2.0)
        sy  (/ stepy 2.0)
        -sx (- sx)
        -sy (- sy)]
    (case (int id)
      0 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 -sy)
         (sym-transform 1.0 0.0 sx 0.0 1.0 sy)]
      1 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 -sy)
         (sym-transform 1.0 0.0 sx 0.0 -1.0 (inc sy))]
      2 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 (- -sy 0.5))
         (sym-transform -1.0 0.0 (inc sx) 0.0 -1.0 (+ sy 0.5))]
      3 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 (- -sy 0.5))
         (sym-transform -1.0 0.0 (inc sx) 0.0 1.0 (+ sy 0.5))]
      4 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 (- -sy 0.5))
         (sym-transform 1.0 0.0 (dec sx) 0.0 -1.0 (+ sy 0.5))]
      5 [(sym-transform 1.0 0.0 (dec -sx) 0.0 1.0 (- -sy 0.5))
         (sym-transform -1.0 0.0 (inc sx) 0.0 1.0 (- sy 0.5))
         (sym-transform 1.0 0.0 (dec -sx) 0.0 -1.0 (+ -sy 0.5))
         (sym-transform -1.0 0.0 (inc sx) 0.0 -1.0 (+ sy 0.5))]
      6 [(sym-transform 1.0 0.0 (- -sx 2.0) 0.0 1.0 (- -sy 0.5))
         (sym-transform -1.0 0.0 (+ sx 2.0) 0.0 1.0 (- sy 0.5))
         (sym-transform 1.0 0.0 sx 0.0 -1.0 (+ sy 0.5))
         (sym-transform -1.0 0.0 -sx 0.0 -1.0 (+ -sy 0.5))])))

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
  (let [spx (/ spacex 2.0)
        spy (/ spacey 2.0)
        -spx (- spx)
        -spy (- spy)
        stx (/ stepx 2.0)
        sty (/ stepy 2.0)
        -stx (- stx)
        -sty (- sty)]
    (case (int id)
      0[(sym-transform -1.0 0.0 -spx 0.0 -1.0 -spy)
        (sym-transform 0.0 1.0 -spx -1.0 0.0 -spy)
        (sym-transform 0.0 -1.0 spx 1.0 0.0 spy)
        (sym-transform 1.0 0.0 spx 0.0 1.0 spy)]
      1 [(sym-transform 1.0 0.0 (+ spacex stx) 0.0 1.0 (+ spacey sty))
         (sym-transform -1.0 0.0 (- stx spacex) 0.0 -1.0 (- sty spacey))
         (sym-transform 0.0 1.0 (+ spacex stx) -1.0 0.0 (- sty spacey))
         (sym-transform 0.0 -1.0 (- stx spacex) 1.0 0.0 (+ spacey sty))
         (sym-transform -1.0 0.0 (- -stx spacex) 0.0 1.0 (- spacey sty))
         (sym-transform 0.0 -1.0 (- -stx spacex) -1.0 0.0 (- -sty spacey))
         (sym-transform 0.0 1.0 (- spacex stx) 1.0 0.0 (- spacey sty))
         (sym-transform 1.0 0.0 (- spacex stx) 0.0 -1.0 (- spacey -sty))])))

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
     (let [r (* amount (v/mag v))
           cr (* amount (m/cos r))
           icr (/ 1.0 (if (zero? cr) m/EPSILON cr))]
       (Vec2. (* amount (.x v)) icr)))))
