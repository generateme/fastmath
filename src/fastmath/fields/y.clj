(ns fastmath.fields.y
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn y
  ([] {:type :regular
       :config (fn [] {:hypergon (r/randval 0.0 (r/drand -3.0 3.0))
                      :hypergon-n (r/randval (u/sirand 1 9) (u/sdrand 0.1 9.0))
                      :hypergon-r (u/sdrand 0.2 1.5)
                      :star (r/randval 0.0 (r/drand -3.0 3.0))
                      :star-n (r/randval (u/sirand 1 9) (u/sdrand 0.1 9.0))                      
                      :star-slope (r/drand -2.0 2.0)
                      :lituus (r/randval 0.0 (r/drand -3.0 3.0))
                      :lituus-a (r/drand -2.0 2.0)
                      :super (r/randval 0.0 (r/drand -3.0 3.0))
                      :super-m (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n1 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n2 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :super-n3 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))})})
  ([^double amount {:keys [^double hypergon ^double hypergon-n ^double hypergon-r
                           ^double star ^double star-n ^double star-slope
                           ^double lituus ^double lituus-a
                           ^double super ^double super-m ^double super-n1 ^double super-n2 ^double super-n3]}]
   (let [hypergon (if (and (zero? hypergon) (zero? star) (zero? lituus) (zero? super)) 1.0 hypergon)
         -hypergon-d (m/sqrt (inc (m/sq hypergon-r)))
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
                                                       (m/pow as super-n3)) -super-n1))))))
             r (* amount (/ (m/sq total) (v/mag v)))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(defn yinyang
  ([] {:type :random
       :config (fn [] {:radius (u/sdrand 0.1 2.0)
                      :dual-t (r/brand)
                      :outside (r/brand)
                      :ang1 (r/drand -2.0 2.0)
                      :ang2 (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double radius ^double ang1 ^double ang2 dual-t outside]}]
   (let [sina (m/sin (* m/PI ang1))
         cosa (m/cos (* m/PI ang1))
         sinb (m/sin (* m/PI ang2))
         cosb (m/sin (* m/PI ang2))]
     (fn [^Vec2 v]
       (let [R2 (v/magsq v)]
         (cond
           (< R2 1.0) (let [d? (and dual-t (r/brand))
                            xx (if d?
                                 (- (* (.x v) cosb) (* (.y v) sinb))
                                 (- (* (.x v) cosa) (* (.y v) sina)))
                            yy (if d?
                                 (+ (* (.x v) sinb) (* (.y v) cosb))
                                 (+ (* (.x v) sina) (* (.y v) cosa)))
                            RR (if d? (- 1.0 radius) radius)
                            iamount (if d? (- amount) amount)]
                        (if (pos? yy)
                          (let [t (m/sqrt (- 1.0 (* yy yy)))
                                k (/ xx t)
                                t1 (* 2.0 (- t 0.5))
                                alpha (* 0.5 (- 1.0 k))
                                beta (- 1.0 alpha)
                                dx (* alpha (dec RR))
                                k1 (+ (* alpha RR) beta)]
                            (Vec2. (* iamount (+ (* t1 k1) dx))
                                   (* iamount (m/sqrt (- 1.0 (* t1 t1))) k1)))
                          (let [RR- (- 1.0 RR)]
                            (Vec2. (* iamount (+ (* xx RR-) RR))
                                   (* iamount yy RR-)))))
           outside (v/mult v amount)
           :else u/zerov))))))
