(ns fastmath.fields.w
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]
           [fastmath.fields.utils JacobiData]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn wdisk
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [a (/ m/PI (inc (v/mag v)))
           r (* amount (* (v/heading v) m/M_1_PI))
           a (if (pos? r) (- m/PI a) a)]
       (Vec2. (* r (m/cos a))
              (* r (m/sin a)))))))


(defn w
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :hypergon (r/randval 0.0 (r/drand 1.0 8.0))
                      :hypergon-n (r/randval (u/sirand 1 9) (u/sdrand 0.1 9.0))
                      :hypergon-r (u/sdrand 0.2 1.5)
                      :star (r/randval 0.0 (r/drand 0.5 4.0))
                      :star-n (r/randval (u/sirand 1 9) (u/sdrand 0.1 9.0))                      
                      :star-slope (r/drand -2.0 2.0)
                      :lituus (r/randval 0.0 (r/drand 0.5 8.0))
                      :lituus-a (r/drand -2.0 2.0)
                      :super (r/randval 0.0 (r/drand 0.5 4.0))
                      :super-m (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n1 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n2 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n3 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))})})
  ([^double amount {:keys [^double angle ^double hypergon ^double hypergon-n ^double hypergon-r
                           ^double star ^double star-n ^double star-slope
                           ^double lituus ^double lituus-a
                           ^double super ^double super-m ^double super-n1 ^double super-n2 ^double super-n3]}]
   (let [-hypergon-d (m/sqrt (inc (m/sq hypergon-r)))
         -lituus-a (- lituus-a)
         -star-slope (m/tan star-slope)
         -super-m (* 0.25 super-m)
         -super-n1 (/ -1.0 super-n1)
         twopi_hypergon-n (/ m/TWO_PI hypergon-n)
         pi_hypergon-n (/ m/PI hypergon-n)
         sq-hypergon-d (m/sq -hypergon-d)
         twopi_star-n (/ m/TWO_PI star-n)
         pi_star-n (/ m/PI star-n)
         sq-star-slope (m/sq -star-slope)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             absa (m/abs a)
             r (v/mag v)
             a2 (+ a angle)
             ^double a2 (cond
                          (< a2 m/-PI) (+ a2 m/TWO_PI)
                          (> a2 m/PI) (- a2 m/TWO_PI)
                          :else a2)
             absa2 (m/abs a2)
             total (as-> 0.0 total
                     (if (zero? hypergon) total
                         (let [temp1 (- (mod absa twopi_hypergon-n) pi_hypergon-n)
                               temp2 (inc (m/sq (m/tan temp1)))]
                           (if (>= temp2 sq-hypergon-d)
                             hypergon
                             (/ (* hypergon
                                   (- -hypergon-d (m/sqrt (- sq-hypergon-d temp2))))
                                (m/sqrt temp2)))))
                     (if (zero? star) total
                         (let [temp1 (m/tan (m/abs
                                             (- (mod absa twopi_star-n)
                                                pi_star-n)))]
                           (+ total (* star (m/sqrt (/ (* sq-star-slope (inc (m/sq temp1)))
                                                       (m/sq (+ temp1 -star-slope))))))))
                     (if (zero? lituus) total
                         (+ total (* lituus (m/pow (inc (/ a m/PI)) -lituus-a))))
                     (if (zero? super) total
                         (let [ang (* a -super-m)
                               as (m/abs (m/sin ang))
                               ac (m/abs (m/cos ang))]
                           (+ total (* super (m/pow (+ (m/pow ac super-n2)
                                                       (m/pow as super-n3)) -super-n1))))))]
         (if (< r total)
           (let [total2 (as-> 0.0 total
                          (if (zero? hypergon) total
                              (let [temp1 (- (mod absa2 twopi_hypergon-n) pi_hypergon-n)
                                    temp2 (inc (m/sq (m/tan temp1)))]
                                (if (>= temp2 sq-hypergon-d)
                                  hypergon
                                  (/ (* hypergon
                                        (- -hypergon-d (m/sqrt (- sq-hypergon-d temp2))))
                                     (m/sqrt temp2)))))
                          (if (zero? star) total
                              (let [temp1 (m/tan (m/abs
                                                  (- (mod absa2 twopi_star-n)
                                                     pi_star-n)))]
                                (+ total (* star (m/sqrt (/ (* sq-star-slope (inc (m/sq temp1)))
                                                            (m/sq (+ temp1 -star-slope))))))))
                          (if (zero? lituus) total
                              (+ total (* lituus (m/pow (inc (/ a2 m/PI)) -lituus-a))))
                          (if (zero? super) total
                              (let [ang (* a2 -super-m)
                                    as (m/abs (m/sin ang))
                                    ac (m/abs (m/cos ang))]
                                (+ total (* super (m/pow (+ (m/pow ac super-n2)
                                                            (m/pow as super-n3)) -super-n1))))))
                 r (/ (* amount total2 r) total)]
             (Vec2. (* r (m/cos a2)) (* r (m/sin a2))))
           (v/mult v amount)))))))


(defn wallpaper
  ([] {:type :random
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double a ^double b ^double c]}]
   (fn [^Vec2 v]
     (-> (r/randval v
                    (Vec2. (- (.y v) (* (m/signum (.x v)) (m/sqrt (m/abs (- (* b (.x v)) c)))))
                           (- a (.x v))))
         (v/mult amount)))))


(defn waves22
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.01 1.0)
                      :scaley (u/sdrand 0.01 1.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :modex (r/brand)
                      :modey (r/brand)
                      :powerx (r/randval (u/sirand 1 4) (u/sdrand 0.05 3.0))
                      :powery (r/randval (u/sirand 1 4) (u/sdrand 0.05 3.0))})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double modex ^double modey ^double powerx ^double powery]}]
   (fn [^Vec2 v]
     (let [sinx (if modex
                  (m/sin (* (.y v) freqx))
                  (m/sin (* 0.5 (inc (m/sin (* (.y v) freqx))))))
           offsetx (* scalex (m/pow (m/abs sinx) powerx))
           siny (if modey
                  (m/sin (* (.x v) freqy))
                  (m/sin (* 0.5 (inc (m/sin (* (.x v) freqy))))))
           offsety (* scaley (m/pow (m/abs siny) powery))]
       (v/mult (v/add v (Vec2. offsetx offsety)) amount)))))

(defn- waves23-shift
  ^double [^double x]
  (if (> x 0.5) (- 0.5 x) x))

(defn waves23
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.01 1.0)
                      :scaley (u/sdrand 0.01 1.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy]}]
   (let [freq (Vec2. (* freqx m/M_1_PI)
                     (* freqy m/M_1_PI))
         scale (Vec2. scalex scaley)]
     (fn [^Vec2 v]
       (let [m (v/emult v freq)]
         (-> (v/sub m (v/floor m))
             (v/fmap waves23-shift)
             (v/emult scale)
             (v/add v)
             (v/mult amount)))))))

(defn- waves2b-safediv
  ^double [^double q ^double r]
  (if (< r 1.0e-10) (/ r) (/ q r)))

(defn waves2b
  ([] {:type :regular
       :config (fn [] {:freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :pwx (r/randval m/THIRD 1.0e-6
                                      (r/randval -1.0e-6 (u/sdrand 1.0e-4 10.0)))
                      :pwy (r/randval m/THIRD 1.0e-6
                                      (r/randval -1.0e-6 (u/sdrand 1.0e-4 10.0)))
                      :scalex (u/sdrand 0.1 1.5)
                      :scaleinfx (u/sdrand 0.1 1.5)
                      :scaley (u/sdrand 0.1 1.5)
                      :scaleinfy (u/sdrand 0.1 1.5)
                      :unity (u/sdrand 0.1 2.0)
                      :jacok (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double freqx ^double freqy ^double pwx ^double pwy
                           ^double scalex ^double scaleinfx ^double scaley ^double scaleinfy
                           ^double unity ^double jacok]}]
   (let [six (- scalex scaleinfx)
         siy (- scaley scaleinfy)
         modex (int (cond
                      (and (>= pwx 0.0) (< pwx 1.0e-4)) 0
                      (and (neg? pwx) (> pwx -1.0e-4)) 1
                      :else 2))
         modey (int (cond
                      (and (>= pwy 0.0) (< pwy 1.0e-4)) 0
                      (and (neg? pwy) (> pwy -1.0e-4)) 1
                      :else 2))]
     (fn [^Vec2 v]
       (let [csx (+ scaleinfx (* (waves2b-safediv unity (+ unity (m/sq (.x v)))) six))
             csy (+ scaleinfy (* (waves2b-safediv unity (+ unity (m/sq (.y v)))) siy))]
         (Vec2. (case modex
                  0 (let [^JacobiData jacobi (u/jacobi-elliptic (* (.y v) freqx) jacok)]
                      (* amount (+ (.x v) (* csx (.sn jacobi)))))
                  1 (* amount (+ (.x v) (* csx (m/bessel-j 1 (m/abs (* (.y v) freqx))))))
                  2 (* amount (+ (.x v) (* csx (m/sin (* (m/sgn (.y v)) (m/pow (+ (m/abs (.y v)) m/EPSILON) pwx) freqx))))))
                (case modey
                  0 (let [^JacobiData jacobi (u/jacobi-elliptic (* (.x v) freqy) jacok)]
                      (* amount (+ (.y v) (* csy (.sn jacobi)))))
                  1 (* amount (+ (.y v) (* csy (m/bessel-j 1 (m/abs (* (.x v) freqy))))))
                  2 (* amount (+ (.y v) (* csy (m/sin (* (m/sgn (.x v)) (m/pow (+ (m/abs (.x v)) m/EPSILON) pwy) freqy))))))))))))

(defn waves2
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy]}]
   (fn [^Vec2 v]
     (Vec2. (* amount (+ (.x v) (* scalex (m/sin (* (.y v) freqx)))))
            (* amount (+ (.y v) (* scaley (m/sin (* (.x v) freqy)))))))))

(defn waves2radial
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :null (r/drand 0.5 2.0)
                      :distance (r/drand 0.5 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double null ^double distance]}]
   (let [diff (- distance null)]
     (fn [^Vec2 v]
       (let [dist (v/mag v)
             factor (if (< dist distance) (/ (- dist null) diff) 1.0)
             factor (if (< dist null) 0.0 factor)]
         (Vec2. (* amount (+ (.x v) (* scalex factor (m/sin (* (.y v) freqx)))))
                (* amount (+ (.y v) (* scaley factor (m/sin (* (.x v) freqy)))))))))))

(defn waves2wf
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :dampx (r/drand -2.0 1.5)
                      :dampy (r/drand -2.0 1.5)
                      :use-cos-x (r/brand)
                      :use-cos-y (r/brand)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double dampx ^double dampy use-cos-x use-cos-y]}]
   (let [dampingx (if (< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (Vec2. (if use-cos-x
                (* amount dampingx (+ (.x v) (* scalex dampingx (m/cos (* (.y v) freqx)))))
                (* amount dampingx (+ (.x v) (* scalex dampingx (m/sin (* (.y v) freqx))))))
              (if use-cos-y
                (* amount dampingy (+ (.y v) (* scaley dampingy (m/cos (* (.x v) freqy)))))
                (* amount dampingy (+ (.y v) (* scaley dampingy (m/sin (* (.x v) freqy)))))))))))

(defn waves3
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :sx-freq (u/sdrand 0.01 10.0)
                      :sy-freq (u/sdrand 0.01 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double sx-freq ^double sy-freq]}]
   (let [hscalex (* 0.5 scalex)
         hscaley (* 0.5 scaley)]
     (fn [^Vec2 v]
       (let [scalexx (* hscalex (inc (m/sin (* (.y v) sx-freq))))
             scaleyy (* hscaley (inc (m/sin (* (.x v) sy-freq))))]
         (Vec2. (* amount (+ (.x v) (* scalexx (m/sin (* (.y v) freqx)))))
                (* amount (+ (.y v) (* scaleyy (m/sin (* (.x v) freqy)))))))))))

(defn waves3wf
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :dampx (r/drand -2.0 1.5)
                      :dampy (r/drand -2.0 1.5)
                      :use-cos-x (r/brand)
                      :use-cos-y (r/brand)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double dampx ^double dampy use-cos-x use-cos-y]}]
   (let [dampingx (if (< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (Vec2. (if use-cos-x
                (* amount dampingx (+ (.x v) (* scalex dampingx (m/sq (m/cos (* (.y v) freqx))))))
                (* amount dampingx (+ (.x v) (* scalex dampingx (m/sq (m/sin (* (.y v) freqx)))))))
              (if use-cos-y
                (* amount dampingy (+ (.y v) (* scaley dampingy (m/sq (m/cos (* (.x v) freqy))))))
                (* amount dampingy (+ (.y v) (* scaley dampingy (m/sq (m/sin (* (.x v) freqy))))))))))))


(defn waves42
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :cont (r/brand 0.2)
                      :yfact (r/drand -2.0 2.0)
                      :freqx2 (u/sdrand 0.01 5.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           cont ^double yfact ^double freqx2]}]
   (let [lowyfact (* 0.001 yfact)]
     (fn [^Vec2 v]
       (let [ax (as-> (m/floor (* (.y v) freqx2)) ax
                  (* 43758.5453 (m/sin (+ (* ax 91.2228) 1.0 (* lowyfact (.y v)))))
                  (- ax (int ax))
                  (if cont (if (> ax 0.5) 1.0 0.0) ax))]
         (Vec2. (* amount (+ (.x v) (* scalex ax ax (m/sin (* (.y v) freqx)))))
                (* amount (+ (.y v) (* scaley (m/sin (* (.x v) freqy)))))))))))

(defn waves4
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :cont (r/brand 0.2)
                      :yfact (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           cont ^double yfact]}]
   (let [lowyfact (* 0.001 yfact)
         freqx2pi (/ freqx m/TWO_PI)]
     (fn [^Vec2 v]
       (let [ax (as-> (m/floor (* (.y v) freqx2pi)) ax
                  (* 43758.5453 (m/sin (+ (* ax 91.2228) 1.0 (* lowyfact (.y v)))))
                  (- ax (int ax))
                  (if cont (if (> ax 0.5) 1.0 0.0) ax))]
         (Vec2. (* amount (+ (.x v) (* scalex ax ax (m/sin (* (.y v) freqx)))))
                (* amount (+ (.y v) (* scaley (m/sin (* (.x v) freqy)))))))))))

(defn waves4wf
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)
                      :dampx (r/drand -2.0 1.5)
                      :dampy (r/drand -2.0 1.5)
                      :use-cos-x (r/brand)
                      :use-cos-y (r/brand)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy
                           ^double dampx ^double dampy use-cos-x use-cos-y]}]
   (let [dampingx (if (< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (let [fx (* (.y v) freqx)
             fy (* (.y v) freqy)]
         (Vec2. (if use-cos-x
                  (* amount dampingx (+ (.x v) (* scalex dampingx (m/sq (m/cos fy)) (m/sin fy))))
                  (* amount dampingx (+ (.x v) (* scalex dampingx (m/sq (m/sin fy)) (m/cos fy)))))
                (if use-cos-y
                  (* amount dampingy (+ (.y v) (* scaley dampingy (m/sq (m/cos fx)) (m/sin fx))))
                  (* amount dampingy (+ (.y v) (* scaley dampingy (m/sq (m/sin fx)) (m/cos fx)))))))))))

(defn waves
  "Waves"
  ([] {:type :regular
       :config (fn [] {:coeff10 (r/drand -2.0 2.0)
                      :coeff11 (r/drand -2.0 2.0)
                      :coeff20 (r/drand -2.0 2.0)
                      :coeff21 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double coeff10 ^double coeff11 ^double coeff20 ^double coeff21]}]
   (let [c202 (+ m/EPSILON (m/sq coeff20))
         c212 (+ m/EPSILON (m/sq coeff21))]
     (fn [^Vec2 v]
       (Vec2. (->> c202
                   (/ (.y v))
                   (m/sin)
                   (* coeff10)
                   (+ (.x v))
                   (* amount))
              (->> c212
                   (/ (.x v))
                   (m/sin)
                   (* coeff11)
                   (+ (.y v))
                   (* amount)))))))

(defn wedge
  "Wedge"
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :hole (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :swirl (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double angle ^double hole ^double count ^double swirl]}]
   (let [hm1p (* m/M_1_PI 0.5)
         comp-fac (- 1.0 (* angle count hm1p))]
     (fn [v]
       (let [r (v/mag v)
             a (+ (v/heading v) (* r swirl))
             c (m/floor (* (+ (* count a) m/PI) hm1p))
             a (+ (* a comp-fac) (* c angle))
             r (* amount (+ r hole))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(defn wedgejulia
  ([] {:type :random
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :dist (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :power (r/randval (u/sirand 1 11) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double angle ^double dist ^double count ^double power]}]
   (let [cf (- 1.0 (* angle count m/M_1_PI 0.5))
         rn (m/abs power)
         cn (* 0.5 (/ dist power))]
     (fn [^Vec2 v]
       (let [r (* amount (m/pow (v/magsq v) cn))
             t-rnd (int (r/drand rn))
             a (/ (+ (v/heading v) (* m/TWO_PI t-rnd)) power)
             c (m/floor (/ (+ (* count a) m/PI) m/TWO_PI))
             a (+ (* a cf) (* c angle))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(defn wedgesph
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :hole (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :swirl (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double angle ^double hole ^double count ^double swirl]}]
   (let [hm1p (* m/M_1_PI 0.5)
         comp-fac (- 1.0 (* angle count hm1p))]
     (fn [v]
       (let [r (/ (+ (v/mag v) m/EPSILON))
             a (+ (v/heading v) (* r swirl))
             c (m/floor (* (+ (* count a) m/PI) hm1p))
             a (+ (* a comp-fac) (* c angle))
             r (* amount (+ r hole))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(defn whorl
  ([] {:type :regular
       :config (fn [] {:inside (r/drand m/-TWO_PI m/TWO_PI)
                      :outside (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double inside ^double outside]}]
   (fn [v]
     (let [r (v/mag v)
           a (if (< r amount)
               (+ (v/heading v) (/ inside (- amount r)))
               (+ (v/heading v) (/ outside (- amount r))))]
       (Vec2. (* amount r (m/cos a))
              (* amount r (m/sin a)))))))

(defn woggle
  ([] {:type :random
       :config (fn [] {:m (r/irand 2 15)})})
  ([^double amount {:keys [^long m]}]
   (let [m (max 2 m)
         a (mapv (fn [^long i] (m/cos (/ (* m/TWO_PI i) m)))(range m))
         b (mapv (fn [^long i] (m/sin (/ (* m/TWO_PI i) m)))(range m))
         r (* (m/sqrt 1.25) (m/sqrt m))
         rr (/ r)
         -rr (- rr)]
     (fn [^Vec2 v]
       (let [c (int (r/drand m))
             ra (/ (* m/SQRT3 (v/mag v)))
             x (if (even? c)
                 (+ (* (.x v) -rr) (* ra (.y v) rr) ^double (a c))
                 (+ (* (.x v) rr)  (* ra (.y v) rr) ^double (a c)))
             y (if (even? c)
                 (+ (* ra (.x v) -rr) (* (.y v) -rr) ^double (b c))
                 (+ (* ra (.x v) -rr) (* (.y v) rr)  ^double (b c)))]
         (Vec2. (* amount x) (* amount y)))))))
