(ns fastmath.fields.w
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]
           [fastmath.fields.utils JacobiData]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn wdisk
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [a (m// m/PI (m/inc (v/mag v)))
           r (m/* amount (m/* (v/heading v) m/M_1_PI))
           a (if (m/pos? r) (m/- m/PI a) a)]
       (Vec2. (m/* r (m/cos a))
              (m/* r (m/sin a)))))))


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
   (let [hypergon (if (and (m/zero? hypergon) (m/zero? star) (m/zero? lituus) (m/zero? super)) 1.0 hypergon)
         -hypergon-d (m/sqrt (m/inc (m/sq hypergon-r)))
         -lituus-a (m/- lituus-a)
         -star-slope (m/tan star-slope)
         -super-m (m/* 0.25 super-m)
         -super-n1 (m// -1.0 super-n1)
         twopi_hypergon-n (m// m/TWO_PI hypergon-n)
         pi_hypergon-n (m// m/PI hypergon-n)
         sq-hypergon-d (m/sq -hypergon-d)
         twopi_star-n (m// m/TWO_PI star-n)
         pi_star-n (m// m/PI star-n)
         sq-star-slope (m/sq -star-slope)]
     (fn [^Vec2 v]
       (let [a (v/heading v)
             absa (m/abs a)
             r (v/mag v)
             a2 (m/+ a angle)
             ^double a2 (cond
                          (m/< a2 m/-PI) (m/+ a2 m/TWO_PI)
                          (m/> a2 m/PI) (m/- a2 m/TWO_PI)
                          :else a2)
             absa2 (m/abs a2)
             total (as-> 0.0 total
                     (if (m/zero? hypergon) total
                         (let [temp1 (m/- (m/mod absa twopi_hypergon-n) pi_hypergon-n)
                               temp2 (m/inc (m/sq (m/tan temp1)))]
                           (if (m/>= temp2 sq-hypergon-d)
                             hypergon
                             (m// (m/* hypergon
                                   (m/- -hypergon-d (m/sqrt (m/- sq-hypergon-d temp2))))
                                (m/sqrt temp2)))))
                     (if (m/zero? star) total
                         (let [temp1 (m/tan (m/abs
                                             (m/- (m/mod absa twopi_star-n)
                                                pi_star-n)))]
                           (m/+ total (m/* star (m/sqrt (m// (m/* sq-star-slope (m/inc (m/sq temp1)))
                                                       (m/sq (m/+ temp1 -star-slope))))))))
                     (if (m/zero? lituus) total
                         (m/+ total (m/* lituus (m/pow (m/inc (m// a m/PI)) -lituus-a))))
                     (if (m/zero? super) total
                         (let [ang (m/* a -super-m)
                               as (m/abs (m/sin ang))
                               ac (m/abs (m/cos ang))]
                           (m/+ total (m/* super (m/pow (m/+ (m/pow ac super-n2)
                                                       (m/pow as super-n3)) -super-n1))))))]
         (if (m/< r total)
           (let [total2 (as-> 0.0 total
                          (if (m/zero? hypergon) total
                              (let [temp1 (m/- (m/mod absa2 twopi_hypergon-n) pi_hypergon-n)
                                    temp2 (m/inc (m/sq (m/tan temp1)))]
                                (if (m/>= temp2 sq-hypergon-d)
                                  hypergon
                                  (m// (m/* hypergon
                                        (m/- -hypergon-d (m/sqrt (m/- sq-hypergon-d temp2))))
                                     (m/sqrt temp2)))))
                          (if (m/zero? star) total
                              (let [temp1 (m/tan (m/abs
                                                  (m/- (m/mod absa2 twopi_star-n)
                                                     pi_star-n)))]
                                (m/+ total (m/* star (m/sqrt (m// (m/* sq-star-slope (m/inc (m/sq temp1)))
                                                            (m/sq (m/+ temp1 -star-slope))))))))
                          (if (m/zero? lituus) total
                              (m/+ total (m/* lituus (m/pow (m/inc (m// a2 m/PI)) -lituus-a))))
                          (if (m/zero? super) total
                              (let [ang (m/* a2 -super-m)
                                    as (m/abs (m/sin ang))
                                    ac (m/abs (m/cos ang))]
                                (m/+ total (m/* super (m/pow (m/+ (m/pow ac super-n2)
                                                            (m/pow as super-n3)) -super-n1))))))
                 r (m// (m/* amount total2 r) total)]
             (Vec2. (m/* r (m/cos a2)) (m/* r (m/sin a2))))
           (v/mult v amount)))))))


(defn wallpaper
  ([] {:type :random
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -3.0 3.0)})})
  ([^double amount {:keys [^double a ^double b ^double c]}]
   (fn [^Vec2 v]
     (-> (r/randval v
                    (Vec2. (m/- (.y v) (m/* (m/signum (.x v)) (m/sqrt (m/abs (m/- (m/* b (.x v)) c)))))
                           (m/- a (.x v))))
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
                  (m/sin (m/* (.y v) freqx))
                  (m/sin (m/* 0.5 (m/inc (m/sin (m/* (.y v) freqx))))))
           offsetx (m/* scalex (m/pow (m/abs sinx) powerx))
           siny (if modey
                  (m/sin (m/* (.x v) freqy))
                  (m/sin (m/* 0.5 (m/inc (m/sin (m/* (.x v) freqy))))))
           offsety (m/* scaley (m/pow (m/abs siny) powery))]
       (v/mult (v/add v (Vec2. offsetx offsety)) amount)))))

(defn- waves23-shift
  ^double [^double x]
  (if (m/> x 0.5) (m/- 0.5 x) x))

(defn waves23
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.01 1.0)
                      :scaley (u/sdrand 0.01 1.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy]}]
   (let [freq (Vec2. (m/* freqx m/M_1_PI)
                     (m/* freqy m/M_1_PI))
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
  (if (m/< r 1.0e-10) (m// r) (m// q r)))

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
   (let [six (m/- scalex scaleinfx)
         siy (m/- scaley scaleinfy)
         modex (int (cond
                      (and (m/>= pwx 0.0) (m/< pwx 1.0e-4)) 0
                      (and (m/neg? pwx) (m/> pwx -1.0e-4)) 1
                      :else 2))
         modey (int (cond
                      (and (m/>= pwy 0.0) (m/< pwy 1.0e-4)) 0
                      (and (m/neg? pwy) (m/> pwy -1.0e-4)) 1
                      :else 2))]
     (fn [^Vec2 v]
       (let [csx (m/+ scaleinfx (m/* (waves2b-safediv unity (m/+ unity (m/sq (.x v)))) six))
             csy (m/+ scaleinfy (m/* (waves2b-safediv unity (m/+ unity (m/sq (.y v)))) siy))]
         (Vec2. (case modex
                  0 (let [^JacobiData jacobi (u/jacobi-elliptic (m/* (.y v) freqx) jacok)]
                      (m/* amount (m/+ (.x v) (m/* csx (.sn jacobi)))))
                  1 (m/* amount (m/+ (.x v) (m/* csx (special/bessel-J1 (m/abs (m/* (.y v) freqx))))))
                  2 (m/* amount (m/+ (.x v) (m/* csx (m/sin (m/* (m/sgn (.y v)) (m/pow (m/+ (m/abs (.y v)) m/EPSILON) pwx) freqx))))))
                (case modey
                  0 (let [^JacobiData jacobi (u/jacobi-elliptic (m/* (.x v) freqy) jacok)]
                      (m/* amount (m/+ (.y v) (m/* csy (.sn jacobi)))))
                  1 (m/* amount (m/+ (.y v) (m/* csy (special/bessel-J1 (m/abs (m/* (.x v) freqy))))))
                  2 (m/* amount (m/+ (.y v) (m/* csy (m/sin (m/* (m/sgn (.x v)) (m/pow (m/+ (m/abs (.x v)) m/EPSILON) pwy) freqy))))))))))))

(defn waves2
  ([] {:type :regular
       :config (fn [] {:scalex (u/sdrand 0.1 2.0)
                      :scaley (u/sdrand 0.1 2.0)
                      :freqx (u/sdrand 0.01 10.0)
                      :freqy (u/sdrand 0.01 10.0)})})
  ([^double amount {:keys [^double scalex ^double scaley ^double freqx ^double freqy]}]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (m/+ (.x v) (m/* scalex (m/sin (m/* (.y v) freqx)))))
            (m/* amount (m/+ (.y v) (m/* scaley (m/sin (m/* (.x v) freqy)))))))))

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
   (let [diff (m/- distance null)]
     (fn [^Vec2 v]
       (let [dist (v/mag v)
             factor (if (m/< dist distance) (m// (m/- dist null) diff) 1.0)
             factor (if (m/< dist null) 0.0 factor)]
         (Vec2. (m/* amount (m/+ (.x v) (m/* scalex factor (m/sin (m/* (.y v) freqx)))))
                (m/* amount (m/+ (.y v) (m/* scaley factor (m/sin (m/* (.x v) freqy)))))))))))

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
   (let [dampingx (if (m/< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (m/< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (Vec2. (if use-cos-x
                (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/cos (m/* (.y v) freqx)))))
                (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/sin (m/* (.y v) freqx))))))
              (if use-cos-y
                (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/cos (m/* (.x v) freqy)))))
                (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/sin (m/* (.x v) freqy)))))))))))

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
   (let [hscalex (m/* 0.5 scalex)
         hscaley (m/* 0.5 scaley)]
     (fn [^Vec2 v]
       (let [scalexx (m/* hscalex (m/inc (m/sin (m/* (.y v) sx-freq))))
             scaleyy (m/* hscaley (m/inc (m/sin (m/* (.x v) sy-freq))))]
         (Vec2. (m/* amount (m/+ (.x v) (m/* scalexx (m/sin (m/* (.y v) freqx)))))
                (m/* amount (m/+ (.y v) (m/* scaleyy (m/sin (m/* (.x v) freqy)))))))))))

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
   (let [dampingx (if (m/< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (m/< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (Vec2. (if use-cos-x
                (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/sq (m/cos (m/* (.y v) freqx))))))
                (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/sq (m/sin (m/* (.y v) freqx)))))))
              (if use-cos-y
                (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/sq (m/cos (m/* (.x v) freqy))))))
                (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/sq (m/sin (m/* (.x v) freqy))))))))))))


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
   (let [lowyfact (m/* 0.001 yfact)]
     (fn [^Vec2 v]
       (let [ax (as-> (m/floor (m/* (.y v) freqx2)) ax
                  (m/* 43758.5453 (m/sin (m/+ (m/* ax 91.2228) 1.0 (m/* lowyfact (.y v)))))
                  (m/- ax (int ax))
                  (if cont (if (m/> ax 0.5) 1.0 0.0) ax))]
         (Vec2. (m/* amount (m/+ (.x v) (m/* scalex ax ax (m/sin (m/* (.y v) freqx)))))
                (m/* amount (m/+ (.y v) (m/* scaley (m/sin (m/* (.x v) freqy)))))))))))

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
   (let [lowyfact (m/* 0.001 yfact)
         freqx2pi (m// freqx m/TWO_PI)]
     (fn [^Vec2 v]
       (let [ax (as-> (m/floor (m/* (.y v) freqx2pi)) ax
                  (m/* 43758.5453 (m/sin (m/+ (m/* ax 91.2228) 1.0 (m/* lowyfact (.y v)))))
                  (m/- ax (int ax))
                  (if cont (if (m/> ax 0.5) 1.0 0.0) ax))]
         (Vec2. (m/* amount (m/+ (.x v) (m/* scalex ax ax (m/sin (m/* (.y v) freqx)))))
                (m/* amount (m/+ (.y v) (m/* scaley (m/sin (m/* (.x v) freqy)))))))))))

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
   (let [dampingx (if (m/< (m/abs dampx) m/EPSILON) 1.0 (m/exp dampx))
         dampingy (if (m/< (m/abs dampy) m/EPSILON) 1.0 (m/exp dampy))]
     (fn [^Vec2 v]
       (let [fx (m/* (.y v) freqx)
             fy (m/* (.y v) freqy)]
         (Vec2. (if use-cos-x
                  (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/sq (m/cos fy)) (m/sin fy))))
                  (m/* amount dampingx (m/+ (.x v) (m/* scalex dampingx (m/sq (m/sin fy)) (m/cos fy)))))
                (if use-cos-y
                  (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/sq (m/cos fx)) (m/sin fx))))
                  (m/* amount dampingy (m/+ (.y v) (m/* scaley dampingy (m/sq (m/sin fx)) (m/cos fx)))))))))))

(defn waves
  "Waves"
  ([] {:type :regular
       :config (fn [] {:coeff10 (r/drand -2.0 2.0)
                      :coeff11 (r/drand -2.0 2.0)
                      :coeff20 (r/drand -2.0 2.0)
                      :coeff21 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double coeff10 ^double coeff11 ^double coeff20 ^double coeff21]}]
   (let [c202 (m/+ m/EPSILON (m/sq coeff20))
         c212 (m/+ m/EPSILON (m/sq coeff21))]
     (fn [^Vec2 v]
       (Vec2. (->> c202
                   (m// (.y v))
                   (m/sin)
                   (m/* coeff10)
                   (m/+ (.x v))
                   (m/* amount))
              (->> c212
                   (m// (.x v))
                   (m/sin)
                   (m/* coeff11)
                   (m/+ (.y v))
                   (m/* amount)))))))

(defn wedge
  "Wedge"
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :hole (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :swirl (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double angle ^double hole ^double count ^double swirl]}]
   (let [hm1p (m/* m/M_1_PI 0.5)
         comp-fac (m/- 1.0 (m/* angle count hm1p))]
     (fn [v]
       (let [r (v/mag v)
             a (m/+ (v/heading v) (m/* r swirl))
             c (m/floor (m/* (m/+ (m/* count a) m/PI) hm1p))
             a (m/+ (m/* a comp-fac) (m/* c angle))
             r (m/* amount (m/+ r hole))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(defn wedgejulia
  ([] {:type :random
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :dist (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :power (r/randval (u/sirand 1 11) (u/sdrand 0.1 10.0))})})
  ([^double amount {:keys [^double angle ^double dist ^double count ^double power]}]
   (let [cf (m/- 1.0 (m/* angle count m/M_1_PI 0.5))
         rn (m/abs power)
         cn (m/* 0.5 (m// dist power))]
     (fn [^Vec2 v]
       (let [r (m/* amount (m/pow (v/magsq v) cn))
             t-rnd (int (r/drand rn))
             a (m// (m/+ (v/heading v) (m/* m/TWO_PI t-rnd)) power)
             c (m/floor (m// (m/+ (m/* count a) m/PI) m/TWO_PI))
             a (m/+ (m/* a cf) (m/* c angle))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(defn wedgesph
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :hole (r/drand -2.0 2.0)
                      :count (r/drand -5.0 5.0)
                      :swirl (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double angle ^double hole ^double count ^double swirl]}]
   (let [hm1p (m/* m/M_1_PI 0.5)
         comp-fac (m/- 1.0 (m/* angle count hm1p))]
     (fn [v]
       (let [r (m// (m/+ (v/mag v) m/EPSILON))
             a (m/+ (v/heading v) (m/* r swirl))
             c (m/floor (m/* (m/+ (m/* count a) m/PI) hm1p))
             a (m/+ (m/* a comp-fac) (m/* c angle))
             r (m/* amount (m/+ r hole))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(defn whorl
  ([] {:type :regular
       :config (fn [] {:inside (r/drand m/-TWO_PI m/TWO_PI)
                      :outside (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double inside ^double outside]}]
   (fn [v]
     (let [r (v/mag v)
           a (if (m/< r amount)
               (m/+ (v/heading v) (m// inside (m/- amount r)))
               (m/+ (v/heading v) (m// outside (m/- amount r))))]
       (Vec2. (m/* amount r (m/cos a))
              (m/* amount r (m/sin a)))))))

(defn woggle
  ([] {:type :random
       :config (fn [] {:m (r/irand 2 15)})})
  ([^double amount {:keys [^long m]}]
   (let [m (m/max 2 m)
         a (mapv (fn [^long i] (m/cos (m// (m/* m/TWO_PI i) m)))(range m))
         b (mapv (fn [^long i] (m/sin (m// (m/* m/TWO_PI i) m)))(range m))
         r (m/* (m/sqrt 1.25) (m/sqrt m))
         rr (m// r)
         -rr (m/- rr)]
     (fn [^Vec2 v]
       (let [c (int (r/drand m))
             ra (m// (m/* m/SQRT3 (v/mag v)))
             x (if (m/even? c)
                 (m/+ (m/* (.x v) -rr) (m/* ra (.y v) rr) ^double (a c))
                 (m/+ (m/* (.x v) rr)  (m/* ra (.y v) rr) ^double (a c)))
             y (if (m/even? c)
                 (m/+ (m/* ra (.x v) -rr) (m/* (.y v) -rr) ^double (b c))
                 (m/+ (m/* ra (.x v) -rr) (m/* (.y v) rr)  ^double (b c)))]
         (Vec2. (m/* amount x) (m/* amount y)))))))

