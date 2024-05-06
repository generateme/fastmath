(ns fastmath.fields.c
  (:refer-clojure :exclude [chunk])
  (:require [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators #{'chunk})

(defn cpow2
  "CPow2"
  ([] {:type :random
       :config (fn [] {:r (u/sdrand 0.1 1.5)
                      :a (u/sdrand 0.1 3.1)
                      :divisor (u/sdrand 0.2 2.0)
                      :range (r/randval (u/sdrand 0.1 10.0)
                                        (let [v (double (r/irand 1 10))
                                              v (r/randval v (/ v))]
                                          (r/randval (* -1.0 v) v)))})})
  ([^double amount {:keys [^double r ^double a ^double divisor ^double range]}]
   (let [ang (/ m/TWO_PI divisor)
         c (/ (* r (m/cos (* m/HALF_PI a))) divisor)
         d (/ (* r (m/sin (* m/HALF_PI a))) divisor)
         half-c (* 0.5 c)
         half-d (* 0.5 d)
         inv-range (/ 0.5 range)
         full-range (* m/TWO_PI range)]
     (fn [^Vec2 v]
       (let [ai (v/heading v)
             n (long (r/drand range))
             n (if (neg? ai) (inc n) n)
             ai (+ ai (* m/TWO_PI n))
             ai (if (< (m/cos (* ai inv-range)) (r/drand -1.0 1.0)) (- ai full-range) ai)
             lnr2 (m/log (v/magsq v))
             r (* amount (m/exp (- (* half-c lnr2) (* d ai))))
             th (+ (* c ai) (* half-d lnr2) (* ang (m/floor (r/drand divisor))))]
         (Vec2. (* r (m/cos th))
                (* r (m/sin th))))))))


(defn cpow3
  "CPow3"
  ([] {:type :random
       :config (fn [] {:r (u/sdrand 0.1 1.5)
                      :d (u/sdrand 0.1 24.0)
                      :divisor (u/sdrand 0.5 1.5)
                      :spread (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double r ^double d ^double divisor ^double spread]}]
   (let [ang (/ m/TWO_PI divisor)
         pa (m/atan2 (if (neg? d) (- (m/log (- d))) (m/log d)) m/TWO_PI)
         cpa (m/cos pa)
         spa (m/sin pa)
         tc (/ (* r cpa cpa) divisor)
         td (/ (* r cpa spa) divisor)
         half-c (* 0.5 tc)
         half-d (* 0.5 td)
         coeff (if (zero? td) 0.0 (/ (* -0.095 spread) td))
         fac (* half-d ang)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             a (if (neg? a) (+ a m/TWO_PI) a)
             a (if (< (m/cos (/ a 2.0)) (r/drand -1.0 1.0)) (- a m/TWO_PI) a)
             a (+ a (* (r/randval m/TWO_PI m/-TWO_PI) coeff (m/round (m/log (r/drand)))))
             lnr2 (m/log (v/magsq v))
             r (* amount (m/exp (- (* half-c lnr2) (* td a))))
             th (+ (* tc a) (* fac lnr2) (* ang (m/floor (r/drand divisor))))]
         (Vec2. (* r (m/cos th))
                (* r (m/sin th))))))))

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
   (let [ang (/ m/TWO_PI divisor)
         c (/ (* r (m/cos (* m/HALF_PI a))) divisor)
         d (/ (* r (m/sin (* m/HALF_PI a))) divisor)
         half-c (* 0.5 c)
         half-d (* 0.5 d)
         inv-spread (/ 0.5 spread)
         full-spread (* m/TWO_PI spread)
         fac (* c half-d ang)]
     (fn [^Vec2 v]
       (let [ai (v/heading v)
             n (r/drand spread)
             n (double (if discrete-spread (int n) n))
             n (if (neg? ai) (inc n) n)
             ai (+ ai (* m/TWO_PI n))
             ai (if (< (m/cos (* ai inv-spread)) (r/drand -1.0 1.0)) (- ai full-spread) ai)
             lnr2 (m/log (v/magsq v))
             ri (* amount (m/exp (- (* half-c lnr2) (* d ai))))
             ang2 (* fac ai lnr2 (+ (r/drand spread2) offset2))]
         (Vec2. (* ri (m/cos ang2))
                (* ri (m/sin ang2))))))))

(defn cpow
  "CPow"
  ([] {:type :random
       :config (fn [] {:r (r/drand -2.0 2.0)
                      :i (r/drand -2.0 2.0)
                      :power (u/sdrand 0.1 12)})})
  ([^double amount {:keys [^double r ^double i ^double power]}]
   (let [va (/ m/TWO_PI power)
         vc (/ r power)
         vd (/ i power)]
     (fn [v]
       (let [a (v/heading v)
             lnr (* 0.5 (m/log (v/magsq v)))
             ang (+ (* a vc)
                    (* vd lnr)
                    (* va (m/floor (r/drand power))))
             m (* amount (m/exp (- (* vc lnr)
                                   (* vd a))))]
         (Vec2. (* m (m/cos ang)) (* m (m/sin ang))))))))

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
           r (* (inc (* 0.9 (m/cos (* 8.0 a))))
                (inc (* 0.1 (m/cos (* 24.0 a))))
                (+ 0.9 (* 0.1 (m/cos (* 200.0 a))))
                (inc (m/sin a)))
           r (* amount (if (> filled (r/drand)) (* r (r/drand)) r))
           a (+ a m/HALF_PI)]
       (Vec2. (* r (m/sin a))
              (* r (m/cos a)))))))

(defn cardioid
  "Cardioid"
  ([] {:type :regular
       :config (fn [] {:a (r/randval 0.7 (u/sirand 1 9) (u/sdrand 0.5 9.0))})})
  ([^double amount {:keys [^double a]}]
   (fn [^Vec2 v]
     (let [ai (v/heading v)
           r (* amount (m/sqrt (+ (v/magsq v) (m/sin (* a ai)) 1.0)))]
       (Vec2. (* r (m/cos ai))
              (* r (m/sin ai)))))))

(defn cell
  "Cell"
  ([] {:type :regular
       :config (fn [] {:size (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double size]}]
   (let [inv-cell-size (/ 1.0 size)
         av (Vec2. amount (- amount))]
     (fn [v]
       (let [^Vec2 xy (-> (v/mult v inv-cell-size)
                          (v/floor))
             dxy (v/sub v (v/mult xy size))
             newxy (if-not (neg? (.y xy))
                     (if-not (neg? (.x xy))
                       (v/mult xy 2.0)
                       (Vec2. (- (inc (* 2.0 (.x xy))))
                              (* 2.0 (.y xy))))
                     (if-not (neg? (.x xy))
                       (Vec2. (* 2.0 (.x xy))
                              (- (inc (* 2.0 (.y xy)))))
                       (Vec2. (- (inc (* 2.0 (.x xy))))
                              (- (inc (* 2.0 (.y xy)))))))]
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
   (let [cs (/ 1.0 (+ size m/EPSILON))
         ncx (* -1.0 x)
         ncy (* -1.0 y)]
     (fn [^Vec2 v]
       (let [isxy (+ (m/round (* cs (.x v)))
                     (m/round (* cs (.y v))))
             dxy (if (even? isxy)
                   (Vec2. (+ ncx (r/drand rnd)) ncy)
                   (Vec2. x (+ y (r/drand rnd))))]
         (-> v
             (v/add dxy)
             (v/mult amount)))))))

(defn chrysantemum
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [u (r/drand (* 21.0 m/PI))
           p4 (m/sin (* u 5.666666666666667))
           p4 (* p4 p4)
           p4 (* p4 p4)
           p8 (m/sin (- (* 2.0 (m/cos (* 3.0 u))) (* 28.0 u)))
           p8 (* p8 p8)
           p8 (* p8 p8)
           p8 (* p8 p8)
           r (* m/THIRD amount (- (* 5.0 (inc (m/sin (* u 2.2))))
                                  (* 4.0 p4 p8)))]
       (Vec2. (* r (m/cos u))
              (* r (m/sin u)))))))

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
     (let [r (* amount (+ (* a (m/sq (.x v)))
                          (* b (.x v) (.y v))
                          (* c (m/sq (.y v)))
                          (* d (.x v))
                          (* e (.y v))
                          f))]
       (if mode
         (if-not (pos? r) v u/zerov)
         (if (pos? r) v u/zerov))))))

(defn circleblur
  "Circle blur"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [rad (m/sqrt (r/drand))
           a (r/drand m/TWO_PI)]
       (Vec2. (* amount (m/cos a) rad)
              (* amount (m/sin a) rad))))))

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
             rdc (+ radius (r/drand (* 0.5 ca)))
             s (m/sin ang)
             c (m/cos ang)]
         (cond
           (<= rad radius) (Vec2. (+ x (* amount (.x v)))
                                  (+ y (* amount (.y v))))
           (and (not zero)
                (> rad radius)) (Vec2. (+ x (* amount rdc c))
                                       (+ y (* amount rdc s)))
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
   (let [dd (* Dens1 Dens2)]
     (fn [^Vec2 v]
       (let [M (->> Sc
                    (/ (.x v))
                    (* 0.5)
                    (m/floor))
             N (->> Sc
                    (/ (.y v))
                    (* 0.5)
                    (m/floor))
             X (- (.x v) (->> M
                              (* 2.0)
                              (inc)
                              (* Sc)))
             Y (- (.y v) (->> N
                              (* 2.0)
                              (inc)
                              (* Sc)))
             U (m/hypot X Y)
             V (->> (r/discrete-noise (+ M 10.0) (+ N 3.0))
                    (* 0.7)
                    (+ 0.3)
                    (* Sc))
             Z1 (r/discrete-noise (+ M Seed) N)
             ^Vec2 XY (if (and (< Z1 Dens1) (< U V)) 
                        (if (pos? Reverse)
                          (if (< Z1 dd)
                            (Vec2. (* K X) (* K Y))
                            (let [Z (->> K
                                         (- 1.0)
                                         (* U)
                                         (/ V)
                                         (+ K))]
                              (Vec2. (* Z X) (* Z Y))))
                          (if (> Z1 dd)
                            (Vec2. (* K X) (* K Y))
                            (let [Z (->> K
                                         (- 1.0)
                                         (* U)
                                         (/ V)
                                         (+ K))]
                              (Vec2. (* Z X) (* Z Y)))))
                        (Vec2. X Y))]
         (Vec2. (->> 2.0
                     (* M)
                     (inc)
                     (* Sc)
                     (+ (.x XY))
                     (* amount))
                (->> 2.0
                     (* N)
                     (inc)
                     (* Sc)
                     (+ (.y XY))
                     (* amount))))))))

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
                                (v/fmap #(inc ^double %))
                                (v/mult Sc)))]
           (if (and (< iter 60)
                    (or (> (r/discrete-noise (+ Seed (.x MN)) (.y MN)) Dens)
                        (> (v/mag XY) (-> (r/discrete-noise (+ 10.0 (.x MN)) (+ 3.0 (.y MN)))
                                          (* 0.7)
                                          (+ 0.3)
                                          (* Sc)))))
             (recur (inc iter))
             (-> MN
                 (v/mult 2.0)
                 (v/fmap #(inc ^double %))
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
         (v/mult (if-not (or (> (r/discrete-noise (+ M Seed) N) Dens)
                             (> U (* Sc (+ 0.3 (* 0.7 (r/discrete-noise (+ M 10) (+ N 3)))))))
                   (loop [i (long 0)]
                     (let [XY (v/emult aAB (Vec2. (r/drand -1.0 1.0) (r/drand -1.0 1.0)))
                           ^Vec2 MN (v/floor (v/div (v/mult XY 0.5) Sc))
                           M (long (.x MN))
                           N (long (.y MN))]
                       (if (and (< i 100) (> (r/discrete-noise (+ M Seed) N) Dens))
                         (recur (inc i))
                         (let [alpha (r/drand m/TWO_PI)
                               s (m/sin alpha)
                               c (m/cos alpha)
                               U (+ 0.3 (* 0.7 (r/discrete-noise (+ M 10) (+ N 3))))]
                           (v/add (Vec2. (* U c) (* U s))
                                  (v/mult (v/shift (v/mult MN 2.0) 1.0) Sc))))))
                   Uxy) amount))))))

(defn circlesplit
  ([] {:type :regular
       :config (fn [] (let [r (r/drand 0.1 2.0)]
                       {:radius r
                        :split (r/drand -2.0 r)}))})
  ([^double amount {:keys [^double radius ^double split]}]
   (let [diff (- radius split)]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (v/mult (if (< r diff)
                   v
                   (let [a (v/heading v)
                         len (+ r split)]
                     (Vec2. (* len (m/cos a))
                            (* len (m/sin a))))) amount))))))

(defn circlize2
  ([] {:type :regular
       :config (fn [] {:hole (r/drand -2.2 2.2)})})
  ([^double amount {:keys [^double hole]}]
   (fn [^Vec2 v]
     (let [absx (m/abs (.x v))
           absy (m/abs (.y v))
           side (if (>= absx absy) absx absy)
           perimeter (if (>= absx absy)
                       (if (> (.x v) absy)
                         (+ absx (.y v))
                         (- (* 5.0 absx) (.y v)))
                       (if (> (.y v) absx)
                         (- (* 3.0 absy) (.x v))
                         (+ (* 7.0 absy) (.x v))))
           r (* amount (+ side hole))
           a (- (/ (* m/M_PI_4 perimeter) side) m/M_PI_4)
           sina (m/sin a)
           cosa (m/cos a)]
       (Vec2. (* r cosa) (* r sina))))))

(defn circlize
  ([] {:type :regular
       :config (fn [] {:hole (r/drand -1.2 1.2)})})
  ([^double amount {:keys [^double hole]}]
   (let [var4-pi (/ amount m/M_PI_4)]
     (fn [^Vec2 v]
       (let [absx (m/abs (.x v))
             absy (m/abs (.y v))
             side (if (>= absx absy) absx absy)
             perimeter (if (>= absx absy)
                         (if (> (.x v) absy)
                           (+ absx (.y v))
                           (- (* 5.0 absx) (.y v)))
                         (if (> (.y v) absx)
                           (- (* 3.0 absy) (.x v))
                           (+ (* 7.0 absy) (.x v))))
             r (+ (* var4-pi side) hole)
             a (- (/ (* m/M_PI_4 perimeter) side) m/M_PI_4)
             sina (m/sin a)
             cosa (m/cos a)]
         (Vec2. (* r cosa) (* r sina)))))))

(defn circular2
  ([] {:type :random
       :config (fn [] {:angle (r/drand -45.0 45.0)
                      :seed (r/drand m/TWO_PI)
                      :xx (r/drand -100.0 100.0)
                      :yy (r/drand -100.0 100.0)})})
  ([^double amount {:keys [^double angle ^double seed ^double xx ^double yy]}]
   (let [ca (m/radians angle)]
     (fn [^Vec2 v]
       (let [aux (m/frac (* 43758.5453 (m/sin (+ (* xx (.x v))
                                                 (* yy (.y v))
                                                 seed))))
             rnd (* (- (* 2.0 (+ (r/drand) aux)) 2.0) ca)
             rad (* amount (v/mag v))
             ang (+ rnd (v/heading v))
             by (m/sin ang)
             bx (m/cos ang)]
         (Vec2. (* rad bx)
                (* rad by)))))))

(defn circular
  ([] {:type :random
       :config (fn [] {:angle (r/drand -45.0 45.0)
                      :seed (r/drand m/TWO_PI)})})
  ([^double amount {:keys [^double angle ^double seed]}]
   (let [ca (m/radians angle)]
     (fn [^Vec2 v]
       (let [aux (m/frac (* 43758.5453 (m/sin (+ (* 12.9898 (.x v))
                                                 (* 78.233 (.y v))
                                                 seed))))
             rnd (* (- (* 2.0 (+ (r/drand) aux)) 2.0) ca)
             rad (* amount (v/mag v))
             ang (+ rnd (v/heading v))
             by (m/sin ang)
             bx (m/cos ang)]
         (Vec2. (* rad bx)
                (* rad by)))))))

(defn circus
  ([] {:type :regular
       :config (fn [] {:scale (u/sdrand 0.2 2.0)})})
  ([^double amount {:keys [^double scale]}]
   (let [scale-1 (/ scale)]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             vn (v/div v r)
             r (if (<= r 1.0) (* r scale) (* r scale-1))]
         (v/mult vn (* amount r)))))))

(defn clifford
  "Clifford Pickover attractor"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (v/mult (Vec2. (+ (m/sin (* a (.y v))) (* c (m/cos (* a (.x v)))))
                    (+ (m/sin (* b (.x v))) (* d (m/cos (* b (.y v)))))) amount))))

(defn cloverleafwf
  "CloverLeafWF"
  ([] {:type :pattern
       :config (fn [] {:filled (r/randval 0.3 0.0 (r/randval 0.3 1.0 (r/drand)))})})
  ([^double amount {:keys [^double filled]}]
   (fn [v]
     (let [a (v/heading v)
           r (+ (m/sin (* 2.0 a))
                (* 0.25 (m/sin (* 6.0 a))))
           r (* amount (if (> filled (r/drand)) (* r (r/drand)) r))]
       (Vec2. (* r (m/sin a))
              (* r (m/cos a)))))))

(defn collideoscope
  ([] {:type :regular
       :config (fn [] {:a (r/drand)
                      :num (r/randval (u/sirand 1 10) (u/sdrand 0.2 10.0))})})
  ([^double amount {:keys [^double a ^long num]}]
   (let [kn-pi (* num m/M_1_PI)
         pi-kn (/ m/PI num)
         ka (* m/PI a)
         ka-kn (/ ka num)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (* amount (v/mag v))
             a (if (>= a 0.0)
                 (let [alt (long (* a kn-pi))]
                   (if (zero? (mod alt 2))
                     (+ (* alt pi-kn) (mod (+ ka-kn a) pi-kn))
                     (+ (* alt pi-kn) (mod (- a ka-kn) pi-kn))))
                 (let [alt (long (* -1.0 a kn-pi))]
                   (if-not (zero? (mod alt 2))
                     (- (+ (* alt pi-kn) (mod (- (- ka-kn) a) pi-kn)))
                     (- (+ (* alt pi-kn) (mod (- ka-kn a) pi-kn))))))]
         (Vec2. (* r (m/cos a)) (* r (m/sin a))))))))

(defn conic
  "Conic"
  ([] {:type :random
       :config (fn [] {:eccentricity (r/drand -3.0 3.0)
                      :holes (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double eccentricity ^double holes]}]
   (fn [^Vec2 v]
     (let [magr (/ 1.0 (+ (v/mag v) m/EPSILON))
           ct (* (.x v) magr)
           r (* (/ (* (* amount (- (r/drand) holes)) eccentricity) (inc (* eccentricity ct))) magr)]
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
                (- (m/pow (m/logb log-base (+ 3.0 (* xs mult-x)))
                          (+ x-power 2.25 xy-power-add)) 1.33)
                (* mult-x (m/pow xs (+ x-power xy-power-add))))
           ey (if log-mode
                (- (m/pow (m/logb log-base (+ 3.0 (* ys mult-y)))
                          (+ y-power 2.25 xy-power-add)) 1.33)
                (* mult-y (m/pow ys (+ y-power xy-power-add))))]
       (Vec2. (if (pos? (.x v)) (+ (* amount ex) x) (- (* amount (- ex)) x))
              (if (pos? (.y v)) (+ (* amount ey) y) (- (* amount (- ey)) y)))))))

(defn cos2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cossin (m/sin (* (.x v) x1))
           coscos (m/cos (* (.x v) x2))
           cossinh (m/sinh (* (.y v) y1))
           coscosh (m/cosh (* (.y v) y2))]
       (Vec2. (* amount coscos coscosh)
              (* -1.0 amount cossin cossinh))))))

(defn cos
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cossin (m/sin (.x v))
           coscos (m/cos (.x v))
           cossinh (m/sinh (.y v))
           coscosh (m/cosh (.y v))]
       (Vec2. (* amount coscos coscosh)
              (* -1.0 amount cossin cossinh))))))

(defn cosh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [coshsin (m/sin (* (.y v) y1))
           coshcos (m/cos (* (.y v) y2))
           coshsinh (m/sinh (* (.x v) x1))
           coshcosh (m/cosh (* (.x v) x2))]
       (Vec2. (* amount coshcos coshcosh)
              (* amount coshsin coshsinh))))))

(defn cosh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [coshsin (m/sin (.y v))
           coshcos (m/cos (.y v))
           coshsinh (m/sinh (.x v))
           coshcosh (m/cosh (.x v))]
       (Vec2. (* amount coshcos coshcosh)
              (* amount coshsin coshsinh))))))

(defn cosine
  "Cosine"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (* m/PI (.x v))]
       (Vec2. (* amount (m/cos r) (m/cosh (.y v)))
              (- (* amount (m/sin r) (m/sinh (.y v)))))))))


(defn cot2bs
  "Cot2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cotsin (m/sin (* x1 (.x v)))
           cotcos (m/cos (* x2 (.x v)))
           cotsinh (m/sinh (* y1 (.y v)))
           cotcosh (m/cosh (* y2 (.y v)))
           cotden (/ (- cotcosh cotcos))]
       (Vec2. (* amount cotden cotsin)
              (* amount cotden -1.0 cotsinh))))))

(defn cot
  "Cot"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cotsin (m/sin (* 2.0 (.x v)))
           cotcos (m/cos (* 2.0 (.x v)))
           cotsinh (m/sinh (* 2.0 (.y v)))
           cotcosh (m/cosh (* 2.0 (.y v)))
           cotden (/ (- cotcosh cotcos))]
       (Vec2. (* amount cotden cotsin)
              (* amount cotden -1.0 cotsinh))))))

(defn coth2bs
  "Coth2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cothsin (m/sin (* y1 (.y v)))
           cothcos (m/cos (* y2 (.y v)))
           cothsinh (m/sinh (* x1 (.x v)))
           cothcosh (m/cosh (* x2 (.x v)))
           cothden (/ (- cothcosh cothcos))]
       (Vec2. (* amount cothden cothsinh)
              (* amount cothden cothsin))))))

(defn coth
  "Coth"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cothsin (m/sin (* 2.0 (.y v)))
           cothcos (m/cos (* 2.0 (.y v)))
           cothsinh (m/sinh (* 2.0 (.x v)))
           cothcosh (m/cosh (* 2.0 (.x v)))
           cothden (/ (- cothcosh cothcos))]
       (Vec2. (* amount cothden cothsinh)
              (* amount cothden cothsin))))))

(defn- crackle-position
  ^Vec2 [^Vec2 xy ^double z ^double s ^double d]
  (let [^Vec2 E (v/mult xy 2.5)
        z (+ (* z 2.5) 19.8)
        sv (if (zero? d)
             u/zerov
             (Vec2. (* d (r/simplex (.x E) (.y E) z))
                    (* d (r/simplex (+ (.x E) 30.2) (- (.y E) 12.1) z))))]
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
   (let [s (/ cellsize 2.0)
         cache (memoize (fn ^Vec2 [xy] (crackle-position xy z s distort)))]
     (doseq [x (range -10 11)
             y (range -10 11)]
       (cache (Vec2. x y)))
     (fn [_]
       (let [blurr (* 2.0 (+ (/ (+ (r/drand) (r/drand)) 2.0)
                             (r/drand -0.125 0.125)))
             theta (r/drand m/TWO_PI)
             U (Vec2. (* blurr (m/sin theta))
                      (* blurr (m/cos theta)))
             Cv (v/floor (v/div U s))
             P (mapv #(cache (v/add Cv %)) u/offsets)
             q (u/closest P 9 U)
             Cv (v/add Cv (u/offsets q))
             P (mapv #(cache (v/add Cv %)) u/offsets)
             L (+ (u/voronoi P 9 4 U) 1.0e-100)
             P4 (P 4)
             Do (v/sub U P4)
             trgL (* (m/pow L power) scale)
             R (/ trgL L)]
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
   (let [[^double top ^double bottom] (if (> top bottom) [bottom top] [top bottom])
         [^double left ^double right] (if (> left right) [right left] [left right])
         [^double top ^double bottom] (if (== top bottom) [-1.0 1.0] [top bottom])
         [^double left ^double right] (if (== left right) [-1.0 1.0] [left right])
         directblur (max 0.0 directblur)
         xinterval (m/abs (- right left))
         yinterval (m/abs (- bottom top))
         xint2 (/ xinterval 2.0)
         yint2 (/ yinterval 2.0)
         minint2 (min xint2 yint2)
         x0 (- right xint2)
         y0 (+ top yint2)
         [^double x0c ^double y0c] (cond
                                     (> xint2 yint2) [(- right minint2) y0]
                                     (< xint2 yint2) [x0 (+ top minint2)]
                                     :else [x0 y0])
         setprob (/ yinterval (+ xinterval yinterval))
         setprobq (* 0.25 setprob)
         setprobh (* 0.5 setprob)
         setprobtq (* 0.75 setprob)
         setcompprob (- 1.0 setprob)
         setcompprobq (+ setprob (* 0.25 setcompprob))
         setcompprobh (+ setprob (* 0.5 setcompprob))
         setcompprobtq (+ setprob (* 0.75 setcompprob))
         mr (* minint2 ratioblur)
         [^double top-border ^double bottom-border
          ^double left-border ^double right-border] (if-not blur
                                                      [top bottom left right]
                                                      [(+ top mr) (- bottom mr)
                                                       (+ left mr) (- right mr)])]
     (fn [^Vec2 v]
       (if (or (< (.x v) left-border)
               (> (.x v) right-border)
               (< (.y v) top-border)
               (> (.y v) bottom-border))
         (if-not blur
           u/zerov
           (let [sectmp (r/drand)
                 [^double xtmp ^double ytmp] (if (< sectmp setprob)
                                               (loop []
                                                 (let [ytmp (+ top (r/drand yint2))
                                                       xtmp (- right (* (m/pow (r/drand) directblur)
                                                                        ratioblur minint2))]
                                                   (if (< (/ (- ytmp y0c)
                                                             (- xtmp x0c)) -1.0)
                                                     (recur)
                                                     [(if (< sectmp setprobh)
                                                        (- (+ left right) xtmp)
                                                        xtmp)
                                                      (if (and (> sectmp setprobq)
                                                               (< sectmp setprobtq))
                                                        (- (+ bottom top) ytmp)
                                                        ytmp)])))
                                               (loop []
                                                 (let [xtmp (- right (r/drand xint2))
                                                       ytmp (+ top (* (m/pow (r/drand) directblur)
                                                                      ratioblur minint2))
                                                       gradtmp (/ (- ytmp y0c)
                                                                  (- xtmp x0c))]
                                                   (if (and (<= gradtmp 0.0)
                                                            (> gradtmp -1.0))
                                                     (recur)
                                                     [(if (and (> sectmp setcompprobq)
                                                               (< sectmp setcompprobtq))
                                                        (- (+ left right) xtmp)
                                                        xtmp)
                                                      (if (> sectmp setcompprobh)
                                                        (- (+ bottom top) ytmp)
                                                        ytmp)]))))]
             (Vec2. (* amount xtmp) (* amount ytmp))))
         (v/mult v amount))))))

(defn cross
  "Cross"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [s (- (m/sq (.x v)) (m/sq (.y v)))
           r (* amount (m/sqrt (/ 1.0 (+ m/EPSILON (* s s)))))]
       (Vec2. (* (.x v) r) (* (.y v) r))))))

(defn csc2bs
  "Csc2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cscsin (m/sin (* x1 (.x v)))
           csccos (m/cos (* x2 (.x v)))
           cscsinh (m/sinh (* y1 (.y v)))
           csccosh (m/cosh (* y2 (.y v)))
           d (- (m/cosh (* 2.0 (.y v)))
                (m/cos (* 2.0 (.x v))))
           cscden (/ 2.0 d)]
       (Vec2. (* amount cscden cscsin csccosh)
              (* amount -1.0 cscden csccos cscsinh))))))

(defn csc
  "Csc"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cscsin (m/sin (.x v))
           csccos (m/cos (.x v))
           cscsinh (m/sinh (.y v))
           csccosh (m/cosh (.y v))
           d (- (m/cosh (* 2.0 (.y v)))
                (m/cos (* 2.0 (.x v))))
           cscden (/ 2.0 d)]
       (Vec2. (* amount cscden cscsin csccosh)
              (* amount -1.0 cscden csccos cscsinh))))))

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
     (let [csc (/ (/ csc-div (m/cos (/ (.x v) cos-div)))
                  (m/tan (/ (.x v) tan-div)))
           fx (* amount (+ (m/pow (+ (* csc csc) (* m/PI pi-mult)) csc-pow)
                           csc-add))]
       (Vec2. (* fx (.x v))
              (* fx (.y v) scale-y))))))

(defn csch2bs
  "Csch2 bs"
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [cschsin (m/sin (* y1 (.y v)))
           cschcos (m/cos (* y2 (.y v)))
           cschsinh (m/sinh (* x1 (.x v)))
           cschcosh (m/cosh (* x2 (.x v)))
           d (- (m/cosh (* 2.0 (.x v)))
                (m/cos (* 2.0 (.y v))))
           cschden (/ 2.0 d)]
       (Vec2. (* amount cschden cschsinh cschcos)
              (* amount -1.0 cschden cschcosh cschsin))))))

(defn csch
  "Csch"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [cschsin (m/sin (.y v))
           cschcos (m/cos (.y v))
           cschsinh (m/sinh (.x v))
           cschcosh (m/cosh (.x v))
           d (- (m/cosh (* 2.0 (.x v)))
                (m/cos (* 2.0 (.y v))))
           cschden (/ 2.0 d)]
       (Vec2. (* amount cschden cschsinh cschcos)
              (* amount -1.0 cschden cschcosh cschsin))))))

(defn curl
  "Curl"
  ([] {:type :regular
       :config (fn [] {:c1 (r/drand -1.0 1.0)
                      :c2 (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double c1 ^double c2]}]
   (fn [^Vec2 v]
     (let [re (inc (+ (* c1 (.x v)) (* c2 (- (m/sq (.x v)) (m/sq (.y v))))))
           im (+ (* c1 (.y v)) (* c2 2.0 (.x v) (.y v)))
           r (/ amount (+ (m/sq re) (m/sq im)))]
       (Vec2. (* r (+ (* (.x v) re) (* (.y v) im)))
              (* r (- (* (.y v) re) (* (.x v) im))))))))

(defn curve
  ([] {:type :regular
       :config (fn [] {:xamp (u/sdrand 0.1 3.0)
                      :yamp (u/sdrand 0.1 3.0)
                      :xlength (r/drand 0.1 2.0)
                      :ylength (r/drand 0.1 2.0)})})
  ([^double amount {:keys [^double xamp ^double yamp ^double xlength ^double ylength]}]
   (let [pc-xlen (/ (- (max m/EPSILON (* xlength xlength))))
         pc-ylen (/ (- (max m/EPSILON (* ylength ylength))))]
     (fn [^Vec2 v]
       (Vec2. (* amount (+ (.x v) (* xamp (m/exp (* (.y v) (.y v) pc-xlen)))))
              (* amount (+ (.y v) (* yamp (m/exp (* (.x v) (.x v) pc-ylen))))))))))

(defn cylinder2
  "Cylinder"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (Vec2. (/ (.x v) (m/sqrt (inc (* (.x v) (.x v))))) (.y v)) amount))))

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
     (if (== (.y v) -1.0)
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
         r2 (* r r)]
     (fn [^Vec2 v]
       (let [diff (v/sub v o)]
         (v/mult (v/add o (v/div (v/mult diff r2) (v/magsq diff))) amount))))))

(m/unuse-primitive-operators #{'chunk})
