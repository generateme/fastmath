(ns fastmath.fields.x
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
             r (* amount (m/sqrt (+ (v/magsq v) (m/sq total))))]
         (Vec2. (* r (m/cos a))
                (* r (m/sin a))))))))

(defn xheart
  ([] {:type :regular
       :config (fn [] {:angle (r/drand m/-TWO_PI m/TWO_PI)
                      :ratio (r/drand -8.0 8.0)})})
  ([^double amount {:keys [^double angle ^double ratio]}]
   (let [ang (+ m/M_PI_4 (* 0.5 m/M_PI_4 angle))
         sina (m/sin ang)
         cosa (m/cos ang)
         rat (+ 6.0 ratio ratio)]
     (fn [^Vec2 v]
       (let [r2-4 (+ 4.0 (v/magsq v))
             r2-4 (if (zero? r2-4) 1.0 r2-4)
             bx (/ 4.0 r2-4)
             by (/ rat r2-4)
             x (- (* cosa bx (.x v)) (* sina by (.y v)))
             y (+ (* sina bx (.x v)) (* cosa by (.y v)))]
         (if (pos? x)
           (Vec2. (* amount x) (* amount y))
           (Vec2. (* amount x) (* -1.0 amount y))))))))

(defn xtrb
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 6) (u/sdrand 1.0 5.0))
                      :radius (u/sdrand 0.3 2.0)
                      :width (r/drand -2.0 2.0)
                      :dist (u/sdrand 0.1 1.5)
                      :a (r/drand m/-TWO_PI m/TWO_PI)
                      :b (r/drand m/-TWO_PI m/TWO_PI)})})
  ([^double amount {:keys [^double power ^double radius ^double dist ^double width ^double a ^double b]}]
   (let [angle-br (+ 0.047 a)
         angle-cr (+ 0.047 b)
         angle-ar (- m/PI angle-br angle-cr)
         sina2 (m/sin (* 0.5 angle-ar))
         cosa2 (m/cos (* 0.5 angle-ar))
         sinb2 (m/sin (* 0.5 angle-br))
         cosb2 (m/cos (* 0.5 angle-br))
         sinc2 (m/sin (* 0.5 angle-cr))
         cosc2 (m/cos (* 0.5 angle-cr))
         sinc (m/sin angle-cr)
         cosc (m/cos angle-cr)
         a (* radius (+ (/ sinc2 cosc2) (/ sinb2 cosb2)))
         b (* radius (+ (/ sinc2 cosc2) (/ sina2 cosa2)))
         c (* radius (+ (/ sina2 cosa2) (/ sina2 cosa2)))
         width1 (- 1.0 width)
         width2 (* 2.0 width)
         width3 (- 1.0 (* width width))
         s2 (* radius (+ a b c))
         ha (/ s2 a 6.0)
         hb (/ s2 b 6.0)
         hc (/ s2 c 6.0)
         ab (/ a b)
         ac (/ a c)
         ba (/ b a)
         bc (/ b c)
         ca (/ c a)
         cb (/ c b)
         s2a (* 6.0 ha)
         s2b (* 6.0 hb)
         s2c (* 6.0 hc)
         s2bc (/ s2 (+ b c) 6.0)
         s2ab (/ s2 (+ a b) 6.0)
         s2ac (/ s2 (+ a c) 6.0)
         absn (long (m/abs power))
         cn (/ dist power 2.0)
         twopi_power (/ m/TWO_PI power)
         direct-trilinear (fn [^double x ^double y]
                            (-> (Vec2. y (- (* x sinc) (* y cosc)))
                                (v/shift radius)))
         inverse-trilinear (fn [^double al ^double be]
                             (let [iny (- al radius)
                                   inx (/ (+ (- be radius)
                                             (* iny cosc)) sinc)
                                   angle (+ (m/atan2 iny inx)
                                            (* twopi_power (r/lrand absn)))
                                   r (* amount (m/pow (+ (* inx inx) (* iny iny)) cn))]
                               (Vec2. (* r (m/cos angle))
                                      (* r (m/sin angle)))))
         hex (fn [^double al ^double be ^double ga]
               (let [R (r/drand)]
                 (if (< be al)
                   (cond
                     (< ga be) (let [de1 (if (>= R width3)
                                           (* width be)
                                           (+ (* width be) (* width2 s2ab (- 3.0 (/ ga be)))))
                                     ga1 (if (>= R width3)
                                           (* width ga)
                                           (+ (* width ga) (* width2 hc (/ ga be))))]
                                 (Vec2. (- s2a (* ba de1) (* ca ga1)) de1))
                     (< ga al) (let [de1 (if (>= R width3)
                                           (* width be)
                                           (+ (* width be) (* width2 hb (/ be ga))))
                                     ga1 (if (>= R width3)
                                           (* width ga)
                                           (+ (* width ga) (* width2 s2ac (- 3.0 (/ be ga)))))]
                                 (Vec2. (- s2a (* ba de1) (* ca ga1)) de1))
                     (>= R width3) (Vec2. (* width al) (* width be))
                     :else (Vec2. (+ (* width1 al) (* width2 s2ac (- 3.0 (/ be al))))
                                  (+ (* width1 be) (* width2 hb (/ be al)))))
                   (cond
                     (< ga al) (let [de1 (if (>= R width3)
                                           (* width al)
                                           (+ (* width al) (* width2 s2ab (- 3.0 (/ ga al)))))
                                     ga1 (if (>= R width3)
                                           (* width ga)
                                           (+ (* width ga) (* width2 hc (/ ga al))))]
                                 (Vec2. (- s2b (* ab de1) (* cb ga1)) de1))
                     (< ga be) (let [de1 (if (>= R width3)
                                           (* width al)
                                           (+ (* width al) (* width2 ha (/ al ga))))
                                     ga1 (if (>= R width3)
                                           (* width ga)
                                           (+ (* width ga) (* width2 s2bc (- 3.0 (/ al ga)))))]
                                 (Vec2. (- s2b (* ab de1) (* cb ga1)) de1))
                     (>= R width3) (Vec2. (* width al) (* width be))
                     :else (Vec2. (+ (* width1 al) (* width2 ha (/ al be)))
                                  (+ (* width1 be) (* width2 s2bc (- 3.0 (/ al be)))))))))]
     (fn [^Vec2 v]
       (let [^Vec2 to (direct-trilinear (.y v) (.x v))
             Alpha (.x to)
             Beta (.y to)
             M (int (m/floor (/ Alpha s2a)))
             ms2a (* M s2a)
             OffsetAl (- Alpha ms2a)
             N (int (m/floor (/ Beta s2b)))
             ns2b (* N s2b)
             OffsetBe (- Beta ns2b)
             OffsetGa (- s2c (* ac OffsetAl) (* bc OffsetBe))
             ^Vec2 alphabeta (if (pos? OffsetGa)
                               (hex OffsetAl OffsetBe OffsetGa)
                               (let [^Vec2 alphabeta (hex (- s2a OffsetAl)
                                                          (- s2b OffsetBe)
                                                          (- OffsetGa))]
                                 (Vec2. (- s2a (.x alphabeta))
                                        (- s2b (.y alphabeta)))))
             Alpha (+ (.x alphabeta) ms2a)
             Beta (+ (.y alphabeta) ns2b)]
         (v/mult (inverse-trilinear Alpha Beta) amount))))))

(m/unuse-primitive-operators)
