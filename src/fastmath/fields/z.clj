(ns fastmath.fields.z
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *warn-on-reflection* true)

(set! *unchecked-math* :warn-on-boxed)
(defn z
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
                                                       (m/pow as super-n3)) -super-n1))))))
             r (m/* amount (m/+ total (v/mag v)))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(defn zsymmetry
  ([] {:type :regular
       :config (fn [] {:p1 (u/sdrand 0.1 5.0)
                      :p2 (u/sdrand 0.1 5.0)})})
  ([^double amount {:keys [^double p1 ^double p2]}]
   (fn [^Vec2 v]
     (let [t (m/* p1 (v/heading v))]
       (-> (Vec2. (m/cos t) (m/sin t))
           (v/mult (m/pow (v/mag v) p2))
           (v/mult amount))))))

