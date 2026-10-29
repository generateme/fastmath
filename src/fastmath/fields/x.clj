(ns fastmath.fields.x
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn x
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
             r (m/* amount (m/sqrt (m/+ (v/magsq v) (m/sq total))))]
         (Vec2. (m/* r (m/cos a))
                (m/* r (m/sin a))))))))

(defn xheart
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/-TWO_PI m/TWO_PI)
                      :ratio (r/drand -8.0 8.0)})})
  ([^double amount {:keys [^double angle ^double ratio]}]
   (let [ang (m/+ m/M_PI_4 (m/* 0.5 m/M_PI_4 angle))
         sina (m/sin ang)
         cosa (m/cos ang)
         rat (m/+ 6.0 ratio ratio)]
     (fn [^Vec2 v]
       (let [r2-4 (m/+ 4.0 (v/magsq v))
             r2-4 (if (m/zero? r2-4) 1.0 r2-4)
             bx (m// 4.0 r2-4)
             by (m// rat r2-4)
             x (m/- (m/* cosa bx (.x v)) (m/* sina by (.y v)))
             y (m/+ (m/* sina bx (.x v)) (m/* cosa by (.y v)))]
         (if (m/pos? x)
           (Vec2. (m/* amount x) (m/* amount y))
           (Vec2. (m/* amount x) (m/* -1.0 amount y))))))))

(defn xtrb
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 6) (u/sdrand 1.0 5.0))
                      :radius (u/sdrand 0.3 2.0)
                      :width (r/drand -2.0 2.0)
                      :dist (u/sdrand 0.1 1.5)
                      :a (r/drand m/-TWO_PI m/TWO_PI)
                      :b (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double power ^double radius ^double dist ^double width ^double a ^double b]}]
   (let [angle-br (m/+ 0.047 a)
         angle-cr (m/+ 0.047 b)
         angle-ar (m/- m/PI angle-br angle-cr)
         sina2 (m/sin (m/* 0.5 angle-ar))
         cosa2 (m/cos (m/* 0.5 angle-ar))
         sinb2 (m/sin (m/* 0.5 angle-br))
         cosb2 (m/cos (m/* 0.5 angle-br))
         sinc2 (m/sin (m/* 0.5 angle-cr))
         cosc2 (m/cos (m/* 0.5 angle-cr))
         sinc (m/sin angle-cr)
         cosc (m/cos angle-cr)
         a (m/* radius (m/+ (m// sinc2 cosc2) (m// sinb2 cosb2)))
         b (m/* radius (m/+ (m// sinc2 cosc2) (m// sina2 cosa2)))
         c (m/* radius (m/+ (m// sina2 cosa2) (m// sina2 cosa2)))
         width1 (m/- 1.0 width)
         width2 (m/* 2.0 width)
         width3 (m/- 1.0 (m/* width width))
         s2 (m/* radius (m/+ a b c))
         ha (m// s2 a 6.0)
         hb (m// s2 b 6.0)
         hc (m// s2 c 6.0)
         ab (m// a b)
         ac (m// a c)
         ba (m// b a)
         bc (m// b c)
         ca (m// c a)
         cb (m// c b)
         s2a (m/* 6.0 ha)
         s2b (m/* 6.0 hb)
         s2c (m/* 6.0 hc)
         s2bc (m// s2 (m/+ b c) 6.0)
         s2ab (m// s2 (m/+ a b) 6.0)
         s2ac (m// s2 (m/+ a c) 6.0)
         absn (long (m/abs power))
         cn (m// dist power 2.0)
         twopi_power (m// m/TWO_PI power)
         direct-trilinear (fn [^double x ^double y]
                            (-> (Vec2. y (m/- (m/* x sinc) (m/* y cosc)))
                                (v/shift radius)))
         inverse-trilinear (fn [^double al ^double be]
                             (let [iny (m/- al radius)
                                   inx (m// (m/+ (m/- be radius)
                                             (m/* iny cosc)) sinc)
                                   angle (m/+ (m/atan2 iny inx)
                                            (m/* twopi_power (r/lrand absn)))
                                   r (m/* amount (m/pow (m/+ (m/* inx inx) (m/* iny iny)) cn))]
                               (Vec2. (m/* r (m/cos angle))
                                      (m/* r (m/sin angle)))))
         hex (fn [^double al ^double be ^double ga]
               (let [R (r/drand)]
                 (if (m/< be al)
                   (cond
                     (m/< ga be) (let [de1 (if (m/>= R width3)
                                           (m/* width be)
                                           (m/+ (m/* width be) (m/* width2 s2ab (m/- 3.0 (m// ga be)))))
                                     ga1 (if (m/>= R width3)
                                           (m/* width ga)
                                           (m/+ (m/* width ga) (m/* width2 hc (m// ga be))))]
                                 (Vec2. (m/- s2a (m/* ba de1) (m/* ca ga1)) de1))
                     (m/< ga al) (let [de1 (if (m/>= R width3)
                                           (m/* width be)
                                           (m/+ (m/* width be) (m/* width2 hb (m// be ga))))
                                     ga1 (if (m/>= R width3)
                                           (m/* width ga)
                                           (m/+ (m/* width ga) (m/* width2 s2ac (m/- 3.0 (m// be ga)))))]
                                 (Vec2. (m/- s2a (m/* ba de1) (m/* ca ga1)) de1))
                     (m/>= R width3) (Vec2. (m/* width al) (m/* width be))
                     :else (Vec2. (m/+ (m/* width1 al) (m/* width2 s2ac (m/- 3.0 (m// be al))))
                                  (m/+ (m/* width1 be) (m/* width2 hb (m// be al)))))
                   (cond
                     (m/< ga al) (let [de1 (if (m/>= R width3)
                                           (m/* width al)
                                           (m/+ (m/* width al) (m/* width2 s2ab (m/- 3.0 (m// ga al)))))
                                     ga1 (if (m/>= R width3)
                                           (m/* width ga)
                                           (m/+ (m/* width ga) (m/* width2 hc (m// ga al))))]
                                 (Vec2. (m/- s2b (m/* ab de1) (m/* cb ga1)) de1))
                     (m/< ga be) (let [de1 (if (m/>= R width3)
                                           (m/* width al)
                                           (m/+ (m/* width al) (m/* width2 ha (m// al ga))))
                                     ga1 (if (m/>= R width3)
                                           (m/* width ga)
                                           (m/+ (m/* width ga) (m/* width2 s2bc (m/- 3.0 (m// al ga)))))]
                                 (Vec2. (m/- s2b (m/* ab de1) (m/* cb ga1)) de1))
                     (m/>= R width3) (Vec2. (m/* width al) (m/* width be))
                     :else (Vec2. (m/+ (m/* width1 al) (m/* width2 ha (m// al be)))
                                  (m/+ (m/* width1 be) (m/* width2 s2bc (m/- 3.0 (m// al be)))))))))]
     (fn [^Vec2 v]
       (let [^Vec2 to (direct-trilinear (.y v) (.x v))
             Alpha (.x to)
             Beta (.y to)
             M (int (m/floor (m// Alpha s2a)))
             ms2a (m/* M s2a)
             OffsetAl (m/- Alpha ms2a)
             N (int (m/floor (m// Beta s2b)))
             ns2b (m/* N s2b)
             OffsetBe (m/- Beta ns2b)
             OffsetGa (m/- s2c (m/* ac OffsetAl) (m/* bc OffsetBe))
             ^Vec2 alphabeta (if (m/pos? OffsetGa)
                               (hex OffsetAl OffsetBe OffsetGa)
                               (let [^Vec2 alphabeta (hex (m/- s2a OffsetAl)
                                                          (m/- s2b OffsetBe)
                                                          (m/- OffsetGa))]
                                 (Vec2. (m/- s2a (.x alphabeta))
                                        (m/- s2b (.y alphabeta)))))
             Alpha (m/+ (.x alphabeta) ms2a)
             Beta (m/+ (.y alphabeta) ns2b)]
         (v/mult (inverse-trilinear Alpha Beta) amount))))))

