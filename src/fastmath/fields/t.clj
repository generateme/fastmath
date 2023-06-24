(ns fastmath.fields.t
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
   (let [aa (* amount a)
         ab (* amount b)
         ac (* amount c)
         ad (* amount d)
         ae (* amount e)
         af (* amount f)
         ag (* amount g)
         mode (int mode)]
     (fn [^Vec2 v]
       (let [x (.x v)
             y (.y v)]
         (if (or (< (+ ad x) l)
                 (< (+ ae y) m))
           (if type
             (Vec2. (* x r) (* y s))
             (Vec2. (* y r) (* x s)))
           (if (and (< x n) (< y o))
             (Vec2. (+ x af) (+ y ag))
             (case mode
               0 (if (and (< (+ x q) amount) (< y aa) (> (+ x ab) l) (> (+ y ac) p))
                   (Vec2. (- (* y h)) (- (* x i)))
                   (Vec2. (* x j) (* y k)))
               1 (if (and (< (+ x q) amount) (< y aa) (> (+ x ab) l) (> (+ y ac) p))
                   (Vec2. (- (* y h)) (- (* x i)))
                   (Vec2. (- (* x j)) (- (* y k))))
               (if (and (< (+ x q) amount) (< y aa) (> (+ x ab) l) (> (+ y ac) p))
                 (Vec2. (* y h) (- (* x i)))
                 (Vec2. (* x j) (- (* y k))))))))))))

(defn tan2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [tansin (m/sin (* (.x v) x1))
           tancos (m/cos (* (.x v) x2))
           tansinh (m/sinh (* (.y v) y1))
           tancosh (m/cosh (* (.y v) y2))
           tanden (/  amount (+ tancos tancosh))]
       (Vec2. (* tanden tansin)
              (* tanden tansinh))))))

(defn tancos
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d1 (+ m/EPSILON (v/magsq v))
           d2 (/ amount d1)]
       (Vec2. (* d2 (m/tanh d1) 2.0 (.x v))
              (* d2 (m/cos d1) 2.0 (.y v)))))))

(defn tan
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [tansin (m/sin (.x v))
           tancos (m/cos (.x v))
           tansinh (m/sinh (.y v))
           tancosh (m/cosh (.y v))
           tanden (/ amount (+ tancos tancosh))]
       (Vec2. (* tanden tansin)
              (* tanden tansinh))))))

(defn tangent
  "Tangent"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d (m/cos (.y v))
           id (/ 1.0 (if (zero? d) m/EPSILON d))]
       (Vec2. (* amount (m/sin (.x v)) id)
              (* amount (m/tan (.y v))))))))

(defn tanh2bs
  ([] {:type :regular
       :config (fn [] {:x1 (u/sdrand 0.1 2.0)
                      :x2 (u/sdrand 0.1 2.0)
                      :y1 (u/sdrand 0.1 2.0)
                      :y2 (u/sdrand 0.1 2.0)})})
  ([^double amount {:keys [^double x1 ^double x2 ^double y1 ^double y2]}]
   (fn [^Vec2 v]
     (let [tanhsin (m/sin (* (.y v) y1))
           tanhcos (m/cos (* (.y v) y2))
           tanhsinh (m/sinh (* (.x v) x1))
           tanhcosh (m/cosh (* (.x v) x2))
           tanhden (/  amount (+ tanhcos tanhcosh))]
       (Vec2. (* tanhden tanhsin)
              (* tanhden tanhsinh))))))

(defn tanh
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [tanhsin (m/sin (* (.y v) 2.0))
           tanhcos (m/cos (* (.y v) 2.0))
           tanhsinh (m/sinh (* (.x v) 2.0))
           tanhcosh (m/cosh (* (.x v) 2.0))
           tanhden (/  amount (+ tanhcos tanhcosh))]
       (Vec2. (* tanhden tanhsin)
              (* tanhden tanhsinh))))))

(defn target
  ([] {:type :regular
       :config (fn [] {:even (r/drand m/TWO_PI)
                      :odd (r/drand m/TWO_PI)
                      :size (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double even ^double odd ^double size]}]
   (let [tsize (* 0.5 size)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             t (as-> (m/log r) t
                 (if (neg? t) (- t tsize) t)
                 (mod (m/abs t) size))
             a (if (< t tsize) (+ a even) (+ a odd))
             ar (* amount r)]
         (Vec2. (* ar (m/cos a)) (* ar (m/sin a))))))))

(defn targetsp
  ([] {:type :regular
       :config (fn [] {:twist (r/drand -1.0 1.0)
                      :n-of-sp (u/sirand 1 8)
                      :size (u/sdrand 0.1 1.5)
                      :tightness (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double twist ^double n-of-sp ^double size ^double tightness]}]
   (let [tsize (* 0.5 size)
         nsp (* n-of-sp m/M_1_PI)
         rota (* m/PI twist)
         rotb (- rota m/PI)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             r (v/mag v)
             t (as-> (-> (m/log r)
                         (* tightness)
                         (+ (* nsp (+ a m/PI)))) t
                 (if (neg? t) (- t tsize) t)
                 (mod (m/abs t) size))
             a (if (< t tsize) (+ a rota) (+ a rotb))
             ar (* amount r)]
         (Vec2. (* ar (m/cos a)) (* ar (m/sin a))))))))

(defn taurus
  "Taurus"
  ([] {:type :regular
       :config (fn [] {:r (r/drand -5.0 5.0)
                      :n (r/drand -5.0 5.0)
                      :inv (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double r ^double n ^double inv]}]
   (let [rinv (* r inv)
         revinv (- 1.0 inv)]
     (fn [^Vec2 v]
       (let [sx (m/sin (.x v))
             cx (m/cos (.x v))
             sy (m/sin (.y v))
             ir (+ rinv (* revinv r (m/cos (* n (.x v)))))
             irsy (+ ir sy)]
         (Vec2. (* amount cx irsy)
                (* amount sx irsy)))))))

(defn threepointifs
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult (cond
               (< (r/drand) m/THIRD) (Vec2. (+ 0.5 (- (* 0.5 (.x v))
                                                      (* 0.5 (.y v))))
                                            (+ 0.5 (- (* -0.5 (.x v))
                                                      (* 0.5 (.y v)))))
               (< (r/drand) m/TWO_THIRD) (Vec2. (.y v) (.x v))
               :else (Vec2. (+ 0.5 (* -0.5 (.y v)))
                            (+ 0.5 (* -0.5 (.x v)))))
             amount))))

(defn tilehlp
  ([] {:type :random
       :config (fn [] {:width (u/sdrand 0.1 1.5)})})
  ([^double amount {:keys [^double width]}]
   (let [width2 (* amount width)
         -width2 (- width2)]
     (fn [^Vec2 v]
       (let [x (/ (.x v) width)
             aux (-> (if (pos? x)
                       (- x (int x))
                       (+ x (int x)))
                     (* m/PI)
                     (m/cos))
             aux2 (if (< aux (r/drand -1.0 1.0)) -width2 width2)]
         (Vec2. (+ aux2 (* amount (.x v)))
                (* amount (.y v))))))))

(defn tile-log
  ([] {:type :random
       :config (fn [] {:spread (u/sdrand 0.5 4.0)})})
  ([^double amount {:keys [^double spread]}]
   (let [-spread (- spread)]
     (fn [^Vec2 v]
       (let [x (r/randval spread -spread)]
         (Vec2. (* amount (+ (.x v) (m/round (* x (m/log (r/drand))))))
                (* amount (.y v))))))))

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
          (Vec2. (* amount (.x v)) (+ space (* amount rev (.y v))))
          (Vec2. (+ space (* amount rev (.x v))) (* amount (.y v))))
        (if vertical
          (Vec2. (* amount (.x v)) (- (* amount rev (.y v)) space))
          (Vec2. (- (* amount rev (.x v)) space) (* amount (.y v)))))))))

(defn trade
  "trade by Michael Faber,  http://michaelfaber.deviantart.com/art/The-Lost-Variations-258913970"
  ([] {:type :regular
       :config (fn [] {:r1 (r/drand 0.1 3.0)
                      :r2 (r/drand 0.1 3.0)
                      :d1 (r/drand -2.0 2.0)
                      :d2 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double r1 ^double r2 ^double d1 ^double d2]}]
   (let [c1 (+ r1 d1)
         c2 (+ r2 d2)]
     (fn [^Vec2 v]
       (let [^Vec4 cr (if (pos? (.x v))
                        (Vec4. c1 (- c2) (/ r2 r1) r1)
                        (Vec4. (- c2) c1 (/ r1 r2) r2))
             nv (Vec2. (- (.x cr) (.x v)) (.y v))
             rm (v/mag nv)
             r (* rm (.z cr))
             a (v/heading nv)
             res (Vec2. (+ (.y cr) (* r (m/cos a)))
                        (* r (m/sin a)))]
         (if (<= rm (.w cr))
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
   (let [tilesize (/ size numbertilesperrow)
         ts2 (* 0.5 tilesize)
         -ts2 (- ts2)
         rng (r/rng :mersenne seed)
         [^double s1 ^double c1] (triantruchet-get-n rng)
         [^double s2 ^double c2] (triantruchet-get-n rng) 
         [^double s3 ^double c3] (triantruchet-get-n rng) 
         [^double s4 ^double c4] (triantruchet-get-n rng) 
         x1 -ts2 y1 -ts2 x2 -ts2
         y2 ts2 x3 ts2 y3 ts2
         shift (Vec2. (* ts2 numbertilesperrow)
                      (* ts2 numbertilespercolumn))]
     (fn [_]
       (let [i (r/irandom rng numbertilesperrow)
             j (r/irandom rng numbertilespercolumn)
             x (+ ts2 (* i tilesize))
             y (+ ts2 (* j tilesize))
             stilt (if (even? i)
                     (if (even? j) s1 s2)
                     (if (even? j) s3 s4))
             ctilt (if (even? i)
                     (if (even? j) c1 c2)
                     (if (even? j) c3 c4))
             p1 (Vec2. (+ (* x1 ctilt) (* y1 stilt) x)
                       (+ (- (* y1 ctilt) (* x1 stilt)) y))
             p2 (Vec2. (+ (* x2 ctilt) (* y2 stilt) x)
                       (+ (- (* y2 ctilt) (* x2 stilt)) y))
             p3 (Vec2. (+ (* x3 ctilt) (* y3 stilt) x)
                       (+ (- (* y3 ctilt) (* x3 stilt)) y))
             sqrt-r1 (m/sqrt (r/drand))
             r2 (r/drand)
             a (- 1.0 sqrt-r1)
             b (* sqrt-r1 (- 1.0 r2))
             c (* r2 sqrt-r1)]
         (-> (Vec2. (+ (* a (.x p1)) (* b (.x p2)) (* c (.x p3)))
                    (+ (* a (.y p1)) (* b (.y p2)) (* c (.y p3))))
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
         seed2 (* 0.25 (/ (m/sqrt (+ m/EPSILON (* 1.5 seed)))
                          (+ m/EPSILON (* 0.5 seed))))
         hseed2 (* 0.5 seed2)]
     (fn [^Vec2 v]
       (let [xs (/ (.x v) scale)
             xp (* 2.0 (m/abs (- xs (m/floor xs) 0.5)))
             xp- (- 1.0 xp)
             width (min 1.0 (+ (* width1 xp-) (* width2 xp)))]
         (if-not (pos? width)
           (v/mult v amount)
           (let [xp2 (+ (* exponent1 xp-) (* exponent2 xp))
                 n (min 2.0 xp2)]
             (if-not (pos? n)
               (v/mult v amount)
               (let [onen (/ n)
                     intx (m/round (.x v))
                     inty (m/round (.y v))
                     r (- (.x v) intx)
                     x (if (neg? r) (inc r) r)
                     r (- (.y v) inty)
                     y (if (neg? r) (inc r) r)
                     tiletype (let [xrand (* intx seed2)
                                    yrand (* inty seed2)
                                    niter (+ xrand yrand (* xrand yrand))
                                    randint (* (+ niter seed) hseed2)
                                    randint (mod (+ (* randint 32747.0) 12345.0) 65535.0)]
                                (mod randint 2.0))
                     r0 (if (< tiletype 1.0)
                          (m/pow (+ (m/pow (m/abs x) n)
                                    (m/pow (m/abs y) n)) onen)
                          (m/pow (+ (m/pow (m/abs (dec x)) n)
                                    (m/pow (m/abs y) n)) onen))
                     r1 (if (< tiletype 1.0)
                          (m/pow (+ (m/pow (m/abs (dec x)) n)
                                    (m/pow (m/abs (dec y)) n)) onen)
                          (m/pow (+ (m/pow (m/abs x) n)
                                    (m/pow (m/abs (dec y)) n)) onen))
                     rmax (* 0.5 (dec (m/pow 2.0 onen)) width)
                     r00 (/ (m/abs (- r0 0.5)) rmax)
                     r11 (/ (m/abs (- r1 0.5)) rmax)
                     xy (Vec2. x y)]
                 (if (or (and (not inverse) (or (< r00 1.0) (< r11 1.0)))
                         (and inverse (> r00 1.0) (> r11 1.0)))
                   (v/mult (v/add (v/floor v) xy) amount)
                   u/zerov))))))))))

(defn truchetfill
  ([] {:type :regular
       :config (fn [] {:pexponent (r/drand 2.0)
                      :arc-width (r/drand)
                      :seed (r/drand 1000.0)})})
  ([^double amount {:keys [^double pexponent ^double arc-width ^double seed]}]
   (let [exponent (m/constrain pexponent 0.001 2.0)
         onen (/ exponent)
         width (m/constrain arc-width 0.001 1.0)
         seed2 (* 0.25 (/ (m/sqrt (* 1.5 seed)) (* 0.5 seed)))
         hseed2 (* 0.5 seed2)
         rmax (* 0.5 (m/pow 2.0 onen) width)
         scale (/ amount)]
     (fn [^Vec2 v]
       (let [xy (v/mult v scale)
             intxy (v/round xy)
             ^Vec2 r (v/sub xy intxy)             
             x (if (neg? (.x r)) (inc (.x r)) (.x r))
             y (if (neg? (.y r)) (inc (.y r)) (.y r))
             tiletype (let [^Vec2 xyrand (v/mult (v/round (v/abs v)) seed2)
                            niter (+ (.x xyrand) (.y xyrand) (* (.x xyrand) (.y xyrand)))
                            randint (* (+ niter seed) hseed2)
                            randint (mod (+ (* randint 32747.0) 12345.0) 65535.0)]
                        (mod randint 2.0))
             r0 (if (< tiletype 1.0)
                  (m/pow (+ (m/pow (m/abs x) exponent)
                            (m/pow (m/abs y) exponent)) onen)
                  (m/pow (+ (m/pow (m/abs (dec x)) exponent)
                            (m/pow (m/abs y) exponent)) onen))
             r1 (if (< tiletype 1.0)
                  (m/pow (+ (m/pow (m/abs (dec x)) exponent)
                            (m/pow (m/abs (dec y)) exponent)) onen)
                  (m/pow (+ (m/pow (m/abs x) exponent)
                            (m/pow (m/abs (dec y)) exponent)) onen))
             xy (Vec2. x y)
             r00 (/ (m/abs (- r0 0.5)) rmax)
             xy1 (if (< r00 1.0)
                   (v/mult (v/add xy (v/floor v)) 2.0)
                   u/zerov)
             r11 (/ (m/abs (- r1 0.5)) rmax)]
         (if (< r11 1.0)
           (v/sub (v/add xy1 (v/mult (v/add xy (v/floor v)) 2.0)) v)
           (v/sub xy1 v)))))))

(def ^:private tf-m (Vec2. 234.234 823.923))

(defn truchetflow
  ([] {:type :pattern
       :config (fn [] {:seed (r/irand)
                      :zoom (u/sdrand 0.1 20.0)
                      :invert (r/brand)})})
  ([^double amount {:keys [^double zoom ^long seed invert]}]
   (let [seed (* seed (m/sin (m/radians seed)))
         amount2 (* 2.0 amount)]
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
             n (m/frac (* (.x p) (.y p)))
             band (- (m/pow (* 0.5 (.y UV)) 1.5) 0.03)
             width (+ 0.01 (* 0.1 (.x UV)))
             gv (if (> n 0.5) (Vec2. (- (.x gv)) (.y gv)) gv)
             s (* -0.5 (m/signum (+ (v/sum gv) 0.001)))
             d (-> (- (v/mag (v/shift gv s)) 0.5)
                   (m/abs)
                   (- band) (m/abs)
                   (- band) (m/abs)
                   (- band) (m/abs)
                   (- width))
             fuzzy (+ 0.01 (* 0.2 (m/sq (v/sum UV2))))
             col (m/smoothstep fuzzy 0.0 d)]
         (if (or (and invert (<= col 0.0))
                 (and (not invert) (> col 0.0)))
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
   (let [d1 (+ 0.5 wd)
         d2 (- 0.5 wd)
         mode (int mode)]
     (fn [^Vec2 v]
       (let [x (- (* m/SQRT3_3 (.x v))
                  (* m/THIRD (.y v)))
             z (* m/TWO_THIRD (.y v))
             y (- (- x) z)
             xyz (Vec3. x y z)
             ^Vec3 r (v/round xyz)
             ^Vec3 diff (v/abs (v/sub r xyz))
             ^Vec2 r (cond
                       (and (> (.x diff) (.y diff))
                            (> (.x diff) (.z diff))) (Vec2. (- (- (.y r)) (.z r)) (.z r))
                       (> (.y diff) (.z diff)) (Vec2. (.x r) (.z r))
                       :else (Vec2. (.x r) (- (- (.x r)) (.y r))))
             ^Vec2 F-h (Vec2. (+ (* m/SQRT3 (.x r)) (* m/SQRT3_2 (.y r)))
                              (* 1.5 (.y r)))
             F (v/sub v F-h)
             add (if (m/one? seed)
                   (if (and (zero? (mod (.x r) 2.0))
                            (zero? (mod (.y r) 2.0)))
                     1.0471975511965976 0.0)
                   (let [hashf (* (m/sin (+ seed (* (.x F-h) 12.9898) (* (.y F-h) 78.233))) 43758.5453)
                         hashf (- hashf (m/floor hashf))]
                     (if (< hashf 0.5) 1.0471975511965976 0.0)))
             angle (- (+ (v/heading F) 0.5235987755982988) add)
             angle2 (+ (/ (m/floor (* angle thc-coeff)) thc-coeff) 0.5235987755982988 add)
             ^Vec2 xy0 (Vec2. (m/cos angle2) (m/sin angle2))
             diff (v/sub F xy0)
             dist (v/mag diff)
             F (if (or (and inv (or (> dist d1) (< dist d2)))
                       (and (not inv) (< dist d1) (> dist d2)))
                 (case mode
                   0 u/zerov
                   1 xy0
                   2 (let [rangle (v/heading diff)
                           D (r/randval d1 d2)]
                       (Vec2. (+ (.x xy0) (* D (m/cos rangle)))
                              (+ (.y xy0) (* D (m/sin rangle))))))
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
   (let [-spreadx (- spreadx)
         -spready (- spready)
         nnn (* 3.0 n)
         nnnn (/ nnn)
         nnnn2pi (* m/TWO_PI nnnn)
         s1 (m/exp (* m/PI nnnn))
         k-factor (/ 2.0 (+ s1 (/ s1)))
         a-factor 0.477464829275686]
     (fn [^Vec2 v]
       (let [rx (m/floor (* (r/randval spreadx -spreadx) (m/log (r/drand))))
             rz (m/floor (* (r/randval spready -spready) (m/log (r/drand))))
             FX-h (+ (* m/SQRT3 rx) (* m/SQRT3_2 (r/randval rz (- rz))))
             FY-h (* 1.5 rz)
             add (not (or (and (not (m/one? seed))
                               (let [hashf (* (m/sin (+ (* FX-h 12.9898)
                                                        (* FY-h 78.233)
                                                        seed)) 43758.5453)
                                     hashf (- hashf (m/floor hashf))]
                                 (< hashf 0.5)))
                          (and (m/one? seed)
                               (zero? (mod rx 2.0))
                               (zero? (mod rz 2.0)))))
             rangle (* (m/floor (r/drand nnn)) nnnn2pi)
             y-aux (if-not flipy (if add (.y v) (- (.y v))) (.y v))
             x-aux (if-not flipx (if add (.x v) (- (.x v))) (.x v))
             FX (* y-aux nnnn)
             FY (* x-aux nnnn)
             aFYpi (+ rangle (* FY m/PI))
             sn (m/sin aFYpi)
             cs (m/cos aFYpi)
             a (* k-factor (m/exp (* FX m/PI)))
             FX (* a cs)
             FY (* a sn)
             A (m/atan2 FY FX)
             A (if (neg? A) (+ m/TWO_PI A) A)
             A (m/floor (* A a-factor))
             ang (* m/THIRD (+ m/PI (* A m/TWO_PI)))
             sn2 (m/sin ang)
             cs2 (m/cos ang)
             FX-new (- FX (* 2.0 cs2))
             FY-new (- FY (* 2.0 sn2))
             FX (if add
                  (- (* m/SQRT3_2 FX-new) (* 0.5 FY-new))
                  (+ (* m/SQRT3_2 FX-new) (* 0.5 FY-new)))
             FY (if add
                  (+ (* m/SQRT3_2 FY-new) (* 0.5 FX-new))
                  (- (* m/SQRT3_2 FY-new) (* 0.5 FX-new)))]
         (Vec2. (* amount (+ (* FX 0.5) FX-h))
                (* amount (+ (* FY 0.5) FY-h))))))))

(defn twintrian
  "Twintrian"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (* amount (r/drand) (v/mag v))
           sinr (m/sin r)
           diff (+ (m/cos r) (m/log10 (m/sq sinr)))]
       (Vec2. (* amount diff (.x v))
              (* amount (.x v) (- diff (* m/PI sinr))))))))

(defn twoface
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (v/mult v (if (pos? (.x v))
                 (/ amount (v/magsq v))
                 amount)))))
