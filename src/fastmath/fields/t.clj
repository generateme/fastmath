(ns fastmath.fields.t
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn tqmirror
  ([] {:type :regular
       :config (fn [] {:a (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :b (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :c (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :d (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :e (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :f (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :g (r/randval (u/sdrand 0.2 1.5) (r/irand -1 2))
                      :h (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :i (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :j (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :k (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :l (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :m (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :n (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :o (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :p (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :q (r/randval 0.1 (r/drand -1.0 1.0) (r/irand -1 2))
                      :r (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :s (r/randval (r/drand -1.2 1.2) (u/sirand 1 2))
                      :mode (r/irand 3)
                      :type (r/brand)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d ^double e
                           ^double f ^double g ^double h ^double i ^double j
                           ^double k ^double l ^double m ^double n ^double o
                           ^double p ^double q ^double r ^double s
                           ^long mode type]}]
   (let [aa (m/* amount a)
         ab (m/* amount b)
         ac (m/* amount c)
         ad (m/* amount d)
         ae (m/* amount e)
         af (m/* amount f)
         ag (m/* amount g)
         mode (int mode)]
     (fn [^Vec2 v]
       (let [x (.x v)
             y (.y v)]
         (if (or (m/< (m/+ ad x) l)
                 (m/< (m/+ ae y) m))
           (if type
             (Vec2. (m/* x r) (m/* y s))
             (Vec2. (m/* y r) (m/* x s)))
           (if (and (m/< x n) (m/< y o))
             (Vec2. (m/+ x af) (m/+ y ag))
             (case mode
               0 (if (and (m/< (m/+ x q) amount) (m/< y aa) (m/> (m/+ x ab) l) (m/> (m/+ y ac) p))
                   (Vec2. (m/- (m/* y h)) (m/- (m/* x i)))
                   (Vec2. (m/* x j) (m/* y k)))
               1 (if (and (m/< (m/+ x q) amount) (m/< y aa) (m/> (m/+ x ab) l) (m/> (m/+ y ac) p))
                   (Vec2. (m/- (m/* y h)) (m/- (m/* x i)))
                   (Vec2. (m/- (m/* x j)) (m/- (m/* y k))))
               (if (and (m/< (m/+ x q) amount) (m/< y aa) (m/> (m/+ x ab) l) (m/> (m/+ y ac) p))
                 (Vec2. (m/* y h) (m/- (m/* x i)))
                 (Vec2. (m/* x j) (m/- (m/* y k))))))))))))

(defn tan2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [tansin (m/sin (m/* (.x v) x1))
           tancos (m/cos (m/* (.x v) x2))
           tansinh (m/sinh (m/* (.y v) y1))
           tancosh (m/cosh (m/* (.y v) y2))
           tanden (m//  amount (m/+ tancos tancosh))]
       (Vec2. (m/* tanden tansin)
              (m/* tanden tansinh))))))

(defn tancos
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d1 (m/+ m/EPSILON (v/magsq v))
           d2 (m// amount d1)]
       (Vec2. (m/* d2 (m/tanh d1) 2.0 (.x v))
              (m/* d2 (m/cos d1) 2.0 (.y v)))))))

(defn tan
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [tansin (m/sin (.x v))
           tancos (m/cos (.x v))
           tansinh (m/sinh (.y v))
           tancosh (m/cosh (.y v))
           tanden (m// amount (m/+ tancos tancosh))]
       (Vec2. (m/* tanden tansin)
              (m/* tanden tansinh))))))

(defn tangent
  "Tangent"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d (m/cos (.y v))
           id (m// 1.0 (if (m/zero? d) m/EPSILON d))]
       (Vec2. (m/* amount (m/sin (.x v)) id)
              (m/* amount (m/tan (.y v))))))))

(defn tanh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [tanhsin (m/sin (m/* (.y v) y1))
           tanhcos (m/cos (m/* (.y v) y2))
           tanhsinh (m/sinh (m/* (.x v) x1))
           tanhcosh (m/cosh (m/* (.x v) x2))
           tanhden (m//  amount (m/+ tanhcos tanhcosh))]
       (Vec2. (m/* tanhden tanhsin)
              (m/* tanhden tanhsinh))))))

(defn tanh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [tanhsin (m/sin (m/* (.y v) 2.0))
           tanhcos (m/cos (m/* (.y v) 2.0))
           tanhsinh (m/sinh (m/* (.x v) 2.0))
           tanhcosh (m/cosh (m/* (.x v) 2.0))
           tanhden (m//  amount (m/+ tanhcos tanhcosh))]
       (Vec2. (m/* tanhden tanhsin)
              (m/* tanhden tanhsinh))))))

(defn target
  ([] {:type :regular
       :config (fn [] {:even (r/drand m/TWO_PI)
                      :odd (r/drand m/TWO_PI)
                      :size (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double even ^double odd ^double size]}]
   (let [tsize (m/* 0.5 size)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             t (as-> (m/log r) t
                 (if (m/neg? t) (m/- t tsize) t)
                 (m/mod (m/abs t) size))
             a (if (m/< t tsize) (m/+ a even) (m/+ a odd))
             ar (m/* amount r)]
         (Vec2. (m/* ar (m/cos a)) (m/* ar (m/sin a))))))))

(defn targetsp
  ([] {:type :regular
       :config (fn [] {:twist (r/drand -1.0 1.0)
                      :n-of-sp (u/sirand 1 8)
                      :size (u/sdrand 0.1 1.5)
                      :tightness (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double twist ^double n-of-sp ^double size ^double tightness]}]
   (let [tsize (m/* 0.5 size)
         nsp (m/* n-of-sp m/M_1_PI)
         rota (m/* m/PI twist)
         rotb (m/- rota m/PI)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             t (as-> (-> (m/log r)
                         (m/* tightness)
                         (m/+ (m/* nsp (m/+ a m/PI)))) t
                 (if (m/neg? t) (m/- t tsize) t)
                 (m/mod (m/abs t) size))
             a (if (m/< t tsize) (m/+ a rota) (m/+ a rotb))
             ar (m/* amount r)]
         (Vec2. (m/* ar (m/cos a)) (m/* ar (m/sin a))))))))

(defn taurus
  "Taurus"
  ([] {:type :regular
       :config (fn [] {:r (r/drand -5.0 5.0)
                      :n (r/drand -5.0 5.0)
                      :inv (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double r ^double n ^double inv]}]
   (let [rinv (m/* r inv)
         revinv (m/- 1.0 inv)]
     (fn [^Vec2 v]
       (let [sx (m/sin (.x v))
             cx (m/cos (.x v))
             sy (m/sin (.y v))
             ir (m/+ rinv (m/* revinv r (m/cos (m/* n (.x v)))))
             irsy (m/+ ir sy)]
         (Vec2. (m/* amount cx irsy)
                (m/* amount sx irsy)))))))

(defn threepointifs
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (cond
               (m/< (r/drand) m/THIRD) (Vec2. (m/+ 0.5 (m/- (m/* 0.5 (.x v))
                                                      (m/* 0.5 (.y v))))
                                            (m/+ 0.5 (m/- (m/* -0.5 (.x v))
                                                      (m/* 0.5 (.y v)))))
               (m/< (r/drand) m/TWO_THIRD) (Vec2. (.y v) (.x v))
               :else (Vec2. (m/+ 0.5 (m/* -0.5 (.y v)))
                            (m/+ 0.5 (m/* -0.5 (.x v)))))
             amount))))

(defn tilehlp
  ([] {:type :random
       :config (fn [] {:width (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double width]}]
   (let [width2 (m/* amount width)
         -width2 (m/- width2)]
     (fn [^Vec2 v]
       (let [x (m// (.x v) width)
             aux (-> (if (m/pos? x)
                       (m/- x (int x))
                       (m/+ x (int x)))
                     (m/* m/PI)
                     (m/cos))
             aux2 (if (m/< aux (r/drand -1.0 1.0)) -width2 width2)]
         (Vec2. (m/+ aux2 (m/* amount (.x v)))
                (m/* amount (.y v))))))))

(defn tile-log
  ([] {:type :random
       :config (fn [] {:spread (u/sdrand 0.5 4.0)})})
  ([^double amount {:keys [^double spread]}]
   (let [-spread (m/- spread)]
     (fn [^Vec2 v]
       (let [x (r/randval spread -spread)]
         (Vec2. (m/* amount (m/+ (.x v) (m/round (m/* x (m/log (r/drand))))))
                (m/* amount (.y v))))))))

(defn tile-reverse
  ([] {:type :random
       :config (fn [] {:space (u/sdrand 0.2 1.5)
                      :reversal (r/brand)
                      :vertical (r/brand)})})
  ([^double amount {:keys [^double space reversal vertical]}]
   (let [rev (if reversal -1.0 1.0)]
     (fn [^Vec2 v]
       (r/randval
        (if vertical
          (Vec2. (m/* amount (.x v)) (m/+ space (m/* amount rev (.y v))))
          (Vec2. (m/+ space (m/* amount rev (.x v))) (m/* amount (.y v))))
        (if vertical
          (Vec2. (m/* amount (.x v)) (m/- (m/* amount rev (.y v)) space))
          (Vec2. (m/- (m/* amount rev (.x v)) space) (m/* amount (.y v)))))))))

(defn trade
  "trade by Michael Faber,  http://michaelfaber.deviantart.com/art/The-Lost-Variations-258913970"
  ([] {:type :regular
       :config (fn [] {:r1 (r/drand 0.1 3.0)
                      :r2 (r/drand 0.1 3.0)
                      :d1 (r/drand -2.0 2.0)
                      :d2 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double r1 ^double r2 ^double d1 ^double d2]}]
   (let [c1 (m/+ r1 d1)
         c2 (m/+ r2 d2)]
     (fn [^Vec2 v]
       (let [^Vec4 cr (if (m/pos? (.x v))
                        (Vec4. c1 (m/- c2) (m// r2 r1) r1)
                        (Vec4. (m/- c2) c1 (m// r1 r2) r2))
             nv (Vec2. (m/- (.x cr) (.x v)) (.y v))
             rm (v/mag nv)
             r (m/* rm (.z cr))
             a (v/heading nv)
             res (Vec2. (m/+ (.y cr) (m/* r (m/cos a)))
                        (m/* r (m/sin a)))]
         (if (m/<= rm (.w cr))
           (v/mult res amount)
           (v/mult v amount)))))))

(defn- triantruchet-get-n
  [rng]
  (case (int (r/irandom rng 4))
    0 [0.0 1.0]
    1 [1.0 0.0]
    2 [0.0 -1.0]
    [-1.0 0.0]))

(defn triantruchet
  ([] {:type :pattern
       :config (fn [] {:seed (r/irand)
                      :size (u/sdrand 1.0 4.0)
                      :numbertilesperrow (r/irand 3 20)
                      :numbertilespercolumn (r/irand 3 20)})})
  ([^double amount {:keys [^long seed ^double size ^long numbertilesperrow ^long numbertilespercolumn]}]
   (let [tilesize (m// size numbertilesperrow)
         ts2 (m/* 0.5 tilesize)
         -ts2 (m/- ts2)
         rng (r/rng :mersenne seed)
         [^double s1 ^double c1] (triantruchet-get-n rng)
         [^double s2 ^double c2] (triantruchet-get-n rng) 
         [^double s3 ^double c3] (triantruchet-get-n rng) 
         [^double s4 ^double c4] (triantruchet-get-n rng) 
         x1 -ts2 y1 -ts2 x2 -ts2
         y2 ts2 x3 ts2 y3 ts2
         shift (Vec2. (m/* ts2 numbertilesperrow)
                      (m/* ts2 numbertilespercolumn))]
     (fn [_]
       (let [i (r/irandom rng numbertilesperrow)
             j (r/irandom rng numbertilespercolumn)
             x (m/+ ts2 (m/* i tilesize))
             y (m/+ ts2 (m/* j tilesize))
             stilt (if (m/even? i)
                     (if (m/even? j) s1 s2)
                     (if (m/even? j) s3 s4))
             ctilt (if (m/even? i)
                     (if (m/even? j) c1 c2)
                     (if (m/even? j) c3 c4))
             p1 (Vec2. (m/+ (m/* x1 ctilt) (m/* y1 stilt) x)
                       (m/+ (m/- (m/* y1 ctilt) (m/* x1 stilt)) y))
             p2 (Vec2. (m/+ (m/* x2 ctilt) (m/* y2 stilt) x)
                       (m/+ (m/- (m/* y2 ctilt) (m/* x2 stilt)) y))
             p3 (Vec2. (m/+ (m/* x3 ctilt) (m/* y3 stilt) x)
                       (m/+ (m/- (m/* y3 ctilt) (m/* x3 stilt)) y))
             sqrt-r1 (m/sqrt (r/drand))
             r2 (r/drand)
             a (m/- 1.0 sqrt-r1)
             b (m/* sqrt-r1 (m/- 1.0 r2))
             c (m/* r2 sqrt-r1)]
         (-> (Vec2. (m/+ (m/* a (.x p1)) (m/* b (.x p2)) (m/* c (.x p3)))
                    (m/+ (m/* a (.y p1)) (m/* b (.y p2)) (m/* c (.y p3))))
             (v/sub shift)
             (v/mult amount)))))))

(defn truchet2
  ([] {:type :regular
       :config (fn [] {:exponent1 (r/drand -1.0 3.0)
                      :exponent2 (r/drand -1.0 3.0)
                      :width1 (r/drand -1.0 2.0)
                      :width2 (r/drand -1.0 2.0)
                      :scale (u/sdrand 0.1 10.0)
                      :seed (r/drand 1000.0)
                      :inverse (r/brand)})})
  ([^double amount {:keys [^double exponent1 ^double exponent2 ^double width1 ^double width2
                           ^double scale ^double seed inverse]}]
   (let [seed (m/abs seed)
         seed2 (m/* 0.25 (m// (m/sqrt (m/+ m/EPSILON (m/* 1.5 seed)))
                          (m/+ m/EPSILON (m/* 0.5 seed))))
         hseed2 (m/* 0.5 seed2)]
     (fn [^Vec2 v]
       (let [xs (m// (.x v) scale)
             xp (m/* 2.0 (m/abs (m/- xs (m/floor xs) 0.5)))
             xp- (m/- 1.0 xp)
             width (m/min 1.0 (m/+ (m/* width1 xp-) (m/* width2 xp)))]
         (if-not (m/pos? width)
           (v/mult v amount)
           (let [xp2 (m/+ (m/* exponent1 xp-) (m/* exponent2 xp))
                 n (m/min 2.0 xp2)]
             (if-not (m/pos? n)
               (v/mult v amount)
               (let [onen (m// n)
                     intx (m/round (.x v))
                     inty (m/round (.y v))
                     r (m/- (.x v) intx)
                     x (if (m/neg? r) (m/inc r) r)
                     r (m/- (.y v) inty)
                     y (if (m/neg? r) (m/inc r) r)
                     tiletype (let [xrand (m/* intx seed2)
                                    yrand (m/* inty seed2)
                                    niter (m/+ xrand yrand (m/* xrand yrand))
                                    randint (m/* (m/+ niter seed) hseed2)
                                    randint (m/mod (m/+ (m/* randint 32747.0) 12345.0) 65535.0)]
                                (m/mod randint 2.0))
                     r0 (if (m/< tiletype 1.0)
                          (m/pow (m/+ (m/pow (m/abs x) n)
                                    (m/pow (m/abs y) n)) onen)
                          (m/pow (m/+ (m/pow (m/abs (m/dec x)) n)
                                    (m/pow (m/abs y) n)) onen))
                     r1 (if (m/< tiletype 1.0)
                          (m/pow (m/+ (m/pow (m/abs (m/dec x)) n)
                                    (m/pow (m/abs (m/dec y)) n)) onen)
                          (m/pow (m/+ (m/pow (m/abs x) n)
                                    (m/pow (m/abs (m/dec y)) n)) onen))
                     rmax (m/* 0.5 (m/dec (m/pow 2.0 onen)) width)
                     r00 (m// (m/abs (m/- r0 0.5)) rmax)
                     r11 (m// (m/abs (m/- r1 0.5)) rmax)
                     xy (Vec2. x y)]
                 (if (or (and (not inverse) (or (m/< r00 1.0) (m/< r11 1.0)))
                         (and inverse (m/> r00 1.0) (m/> r11 1.0)))
                   (v/mult (v/add (v/floor v) xy) amount)
                   u/zerov))))))))))

(defn truchetfill
  ([] {:type :regular
       :config (fn [] {:pexponent (r/drand 2.0)
                      :arc-width (r/drand)
                      :seed (r/drand 1000.0)})})
  ([^double amount {:keys [^double pexponent ^double arc-width ^double seed]}]
   (let [exponent (m/constrain pexponent 0.001 2.0)
         onen (m// exponent)
         width (m/constrain arc-width 0.001 1.0)
         seed2 (m/* 0.25 (m// (m/sqrt (m/* 1.5 seed)) (m/* 0.5 seed)))
         hseed2 (m/* 0.5 seed2)
         rmax (m/* 0.5 (m/pow 2.0 onen) width)
         scale (m// amount)]
     (fn [^Vec2 v]
       (let [xy (v/mult v scale)
             intxy (v/round xy)
             ^Vec2 r (v/sub xy intxy)             
             x (if (m/neg? (.x r)) (m/inc (.x r)) (.x r))
             y (if (m/neg? (.y r)) (m/inc (.y r)) (.y r))
             tiletype (let [^Vec2 xyrand (v/mult (v/round (v/abs v)) seed2)
                            niter (m/+ (.x xyrand) (.y xyrand) (m/* (.x xyrand) (.y xyrand)))
                            randint (m/* (m/+ niter seed) hseed2)
                            randint (m/mod (m/+ (m/* randint 32747.0) 12345.0) 65535.0)]
                        (m/mod randint 2.0))
             r0 (if (m/< tiletype 1.0)
                  (m/pow (m/+ (m/pow (m/abs x) exponent)
                            (m/pow (m/abs y) exponent)) onen)
                  (m/pow (m/+ (m/pow (m/abs (m/dec x)) exponent)
                            (m/pow (m/abs y) exponent)) onen))
             r1 (if (m/< tiletype 1.0)
                  (m/pow (m/+ (m/pow (m/abs (m/dec x)) exponent)
                            (m/pow (m/abs (m/dec y)) exponent)) onen)
                  (m/pow (m/+ (m/pow (m/abs x) exponent)
                            (m/pow (m/abs (m/dec y)) exponent)) onen))
             xy (Vec2. x y)
             r00 (m// (m/abs (m/- r0 0.5)) rmax)
             xy1 (if (m/< r00 1.0)
                   (v/mult (v/add xy (v/floor v)) 2.0)
                   u/zerov)
             r11 (m// (m/abs (m/- r1 0.5)) rmax)]
         (if (m/< r11 1.0)
           (v/sub (v/add xy1 (v/mult (v/add xy (v/floor v)) 2.0)) v)
           (v/sub xy1 v)))))))

(def ^:private tf-m (Vec2. 234.234 823.923))

(defn truchetflow
  ([] {:type :pattern
       :config (fn [] {:seed (r/irand)
                      :zoom (u/sdrand 0.1 20.0)
                      :invert (r/brand)})})
  ([^double amount {:keys [^double zoom ^long seed invert]}]
   (let [seed (m/* seed (m/sin (m/radians seed)))
         amount2 (m/* 2.0 amount)]
     (fn [_]
       (let [^Vec2 UV (v/generate-vec2 r/drand)
             uv (v/shift UV -0.5)
             UV2 (v/abs uv)
             uv (v/mult (v/shift uv seed) zoom)
             ^Vec2 gv (v/shift (v/frac uv) -0.5)
             id (v/floor uv)
             ;; N21
             p (v/frac (v/emult id tf-m))
             ^Vec2 p (v/shift p (v/dot p (v/shift p 42.34)))
             n (m/frac (m/* (.x p) (.y p)))
             band (m/- (m/pow (m/* 0.5 (.y UV)) 1.5) 0.03)
             width (m/+ 0.01 (m/* 0.1 (.x UV)))
             gv (if (m/> n 0.5) (Vec2. (m/- (.x gv)) (.y gv)) gv)
             s (m/* -0.5 (m/signum (m/+ (v/sum gv) 0.001)))
             d (-> (m/- (v/mag (v/shift gv s)) 0.5)
                   (m/abs)
                   (m/- band) (m/abs)
                   (m/- band) (m/abs)
                   (m/- band) (m/abs)
                   (m/- width))
             fuzzy (m/+ 0.01 (m/* 0.2 (m/sq (v/sum UV2))))
             col (m/smoothstep fuzzy 0.0 d)]
         (if (or (and invert (m/<= col 0.0))
                 (and (not invert) (m/> col 0.0)))
           u/zerov
           (v/mult (v/shift UV -0.5) amount2)))))))

(def ^{:private true :const true :tag 'double} thc-coeff 0.477464829275686)

(defn truchethexcrop
  ([] {:type :random
       :config (fn [] {:seed (r/irand)
                      :inv (r/brand)
                      :mode (r/irand 3)
                      :wd (r/drand 0.05 0.45)})})
  ([^double amount {:keys [^long seed inv ^long mode ^double wd]}]
   (let [d1 (m/+ 0.5 wd)
         d2 (m/- 0.5 wd)
         mode (int mode)]
     (fn [^Vec2 v]
       (let [x (m/- (m/* m/SQRT3_3 (.x v))
                  (m/* m/THIRD (.y v)))
             z (m/* m/TWO_THIRD (.y v))
             y (m/- (m/- x) z)
             xyz (Vec3. x y z)
             ^Vec3 r (v/round xyz)
             ^Vec3 diff (v/abs (v/sub r xyz))
             ^Vec2 r (cond
                       (and (m/> (.x diff) (.y diff))
                            (m/> (.x diff) (.z diff))) (Vec2. (m/- (m/- (.y r)) (.z r)) (.z r))
                       (m/> (.y diff) (.z diff)) (Vec2. (.x r) (.z r))
                       :else (Vec2. (.x r) (m/- (m/- (.x r)) (.y r))))
             ^Vec2 F-h (Vec2. (m/+ (m/* m/SQRT3 (.x r)) (m/* m/SQRT3_2 (.y r)))
                              (m/* 1.5 (.y r)))
             F (v/sub v F-h)
             add (if (m/one? seed)
                   (if (and (m/zero? (m/mod (.x r) 2.0))
                            (m/zero? (m/mod (.y r) 2.0)))
                     1.0471975511965976 0.0)
                   (let [hashf (m/* (m/sin (m/+ seed (m/* (.x F-h) 12.9898) (m/* (.y F-h) 78.233))) 43758.5453)
                         hashf (m/- hashf (m/floor hashf))]
                     (if (m/< hashf 0.5) 1.0471975511965976 0.0)))
             angle (m/- (m/+ (v/heading F) 0.5235987755982988) add)
             angle2 (m/+ (m// (m/floor (m/* angle thc-coeff)) thc-coeff) 0.5235987755982988 add)
             ^Vec2 xy0 (Vec2. (m/cos angle2) (m/sin angle2))
             diff (v/sub F xy0)
             dist (v/mag diff)
             F (if (or (and inv (or (m/> dist d1) (m/< dist d2)))
                       (and (not inv) (m/< dist d1) (m/> dist d2)))
                 (case mode
                   0 u/zerov
                   1 xy0
                   2 (let [rangle (v/heading diff)
                           D (r/randval d1 d2)]
                       (Vec2. (m/+ (.x xy0) (m/* D (m/cos rangle)))
                              (m/+ (.y xy0) (m/* D (m/sin rangle))))))
                 F)]
         (v/mult (v/add F F-h) amount))))))

(defn truchethexfill
  ([] {:type :random
       :config (fn [] {:seed (r/irand)
                      :flipx (r/brand)
                      :flipy (r/brand)
                      :n (u/sirand 1 20)
                      :spreadx (u/sdrand 0.1 3.0)
                      :spready (u/sdrand 0.1 3.0)})})
  ([^double amount {:keys [^double spreadx ^double spready ^long n ^long seed flipx flipy]}]
   (let [-spreadx (m/- spreadx)
         -spready (m/- spready)
         nnn (m/* 3.0 n)
         nnnn (m// nnn)
         nnnn2pi (m/* m/TWO_PI nnnn)
         s1 (m/exp (m/* m/PI nnnn))
         k-factor (m// 2.0 (m/+ s1 (m// s1)))
         a-factor 0.477464829275686]
     (fn [^Vec2 v]
       (let [rx (m/floor (m/* (r/randval spreadx -spreadx) (m/log (r/drand))))
             rz (m/floor (m/* (r/randval spready -spready) (m/log (r/drand))))
             FX-h (m/+ (m/* m/SQRT3 rx) (m/* m/SQRT3_2 (r/randval rz (m/- rz))))
             FY-h (m/* 1.5 rz)
             add (not (or (and (not (m/one? seed))
                               (let [hashf (m/* (m/sin (m/+ (m/* FX-h 12.9898)
                                                        (m/* FY-h 78.233)
                                                        seed)) 43758.5453)
                                     hashf (m/- hashf (m/floor hashf))]
                                 (m/< hashf 0.5)))
                          (and (m/one? seed)
                               (m/zero? (m/mod rx 2.0))
                               (m/zero? (m/mod rz 2.0)))))
             rangle (m/* (m/floor (r/drand nnn)) nnnn2pi)
             y-aux (if-not flipy (if add (.y v) (m/- (.y v))) (.y v))
             x-aux (if-not flipx (if add (.x v) (m/- (.x v))) (.x v))
             FX (m/* y-aux nnnn)
             FY (m/* x-aux nnnn)
             aFYpi (m/+ rangle (m/* FY m/PI))
             sn (m/sin aFYpi)
             cs (m/cos aFYpi)
             a (m/* k-factor (m/exp (m/* FX m/PI)))
             FX (m/* a cs)
             FY (m/* a sn)
             A (m/atan2 FY FX)
             A (if (m/neg? A) (m/+ m/TWO_PI A) A)
             A (m/floor (m/* A a-factor))
             ang (m/* m/THIRD (m/+ m/PI (m/* A m/TWO_PI)))
             sn2 (m/sin ang)
             cs2 (m/cos ang)
             FX-new (m/- FX (m/* 2.0 cs2))
             FY-new (m/- FY (m/* 2.0 sn2))
             FX (if add
                  (m/- (m/* m/SQRT3_2 FX-new) (m/* 0.5 FY-new))
                  (m/+ (m/* m/SQRT3_2 FX-new) (m/* 0.5 FY-new)))
             FY (if add
                  (m/+ (m/* m/SQRT3_2 FY-new) (m/* 0.5 FX-new))
                  (m/- (m/* m/SQRT3_2 FY-new) (m/* 0.5 FX-new)))]
         (Vec2. (m/* amount (m/+ (m/* FX 0.5) FX-h))
                (m/* amount (m/+ (m/* FY 0.5) FY-h))))))))

(defn twintrian
  "Twintrian"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/* amount (r/drand) (v/mag v))
           sinr (m/sin r)
           diff (m/+ (m/cos r) (m/log10 (m/sq sinr)))]
       (Vec2. (m/* amount diff (.x v))
              (m/* amount (.x v) (m/- diff (m/* m/PI sinr))))))))

(defn twoface
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult v (if (m/pos? (.x v))
                 (m// amount (v/magsq v))
                 amount)))))

