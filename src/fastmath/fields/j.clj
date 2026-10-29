(ns fastmath.fields.j
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2 Vec4]
           [fastmath.fields.utils JacobiData]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private jacasn-mpi4 (c/complex m/M_PI_4 0.0))

(defn jacasn
  ([] {:type :regular
       :config (fn [] {:kr (r/drand -1.5 1.5)
                      :ki (r/drand -1.5 1.5)
                      :type (r/irand 8)})})
  ([^double amount {:keys [^double kr ^double ki ^long type]}]
   (let [type (m/max 0 (m/rem type 8))
         modul (c/complex kr ki)
         geo (c/complex amount 0.0)]
     (fn [^Vec2 v]
       (let [prephi (case type
                      (0 4) (-> (c/sin v)
                                (c/mult modul)
                                (c/asin)
                                (c/add v)
                                (c/scale 0.5))
                      (1 5) (-> (c/mult modul v)
                                (c/asin)
                                (c/add (c/asin v))
                                (c/scale 0.5))
                      (2 6) (-> (c/mult v v)
                                (->> (c/sub c/ONE))
                                (c/sqrt)
                                (c/mult modul)
                                (c/asin)
                                (c/add (c/asin v))
                                (c/scale 0.5))
                      (3 7) (let [phi (c/mult c/I (c/asinh v))]
                              (-> (c/sin phi)
                                  (c/mult modul)
                                  (c/asin)
                                  (c/add phi)
                                  (c/scale 0.5))))
             phi (if (m/> type 3) modul prephi)
             modul (if (m/> type 3) prephi modul)
             ^fastmath.fields.utils.Pair pg (loop [jj (long 0)
                                                   phi phi
                                                   geo geo
                                                   modul modul]
                                              (if (or (m/> jj 9)
                                                      (m/< (v/magsq (c/sub c/ONE modul)) m/EPSILON))
                                                (fastmath.fields.utils.Pair. phi geo)
                                                (let [geo1 (c/div c/TWO (c/add c/ONE modul))]
                                                  (recur (m/inc jj)
                                                         (if (m/zero? jj) phi (-> (c/sin phi)
                                                                                (c/mult modul)
                                                                                (c/asin)
                                                                                (c/add phi)
                                                                                (c/scale 0.5)))
                                                         (c/mult geo geo1)
                                                         (c/mult geo1 (c/sqrt modul))))))]
         (-> (c/scale (.a pg) 0.5)
             (c/tan)
             (c/add jacasn-mpi4)
             (c/log)
             (c/mult (.b pg))))))))

(defn jaccn
  ([] {:type :regular
       :config (fn [] {:k (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double k]}]
   (fn [^Vec2 v]
     (let [^JacobiData jac-x (u/jacobi-elliptic (.x v) k)
           ^JacobiData jac-y (u/jacobi-elliptic (.y v) (m/- 1.0 k))
           numx (m/* (.cn jac-x) (.cn jac-y))
           numy (m/* (m/- (.dn jac-x)) (.sn jac-x) (.dn jac-y) (.sn jac-y))
           denom (m// amount (m/+ (m/* (.sn jac-x) (.sn jac-x) (.sn jac-y) (.sn jac-y) k)
                              (m/* (.cn jac-y) (.cn jac-y))
                              m/EPSILON))]
       (Vec2. (m/* denom numx) (m/* denom numy))))))

(defn jacdn
  ([] {:type :regular
       :config (fn [] {:k (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double k]}]
   (fn [^Vec2 v]
     (let [^JacobiData jac-x (u/jacobi-elliptic (.x v) k)
           ^JacobiData jac-y (u/jacobi-elliptic (.y v) (m/- 1.0 k))
           numx (m/* (.dn jac-x) (.cn jac-y) (.dn jac-y))
           numy (m/* (m/- (.cn jac-x)) (.sn jac-x) (.sn jac-y) k)
           denom (m// amount (m/+ (m/* (.sn jac-x) (.sn jac-x) (.sn jac-y) (.sn jac-y) k)
                              (m/* (.cn jac-y) (.cn jac-y))
                              m/EPSILON))]
       (Vec2. (m/* denom numx) (m/* denom numy))))))

(defn- jacelk-rf
  ^double [^double X ^double Y ^double Z]
  (let [^Vec4 d (loop [X X Y Y Z Z i (long 0)]
                  (let [lambda (m/+ (m/sqrt (m/* X Y))
                                  (m/sqrt (m/* Y Z))
                                  (m/sqrt (m/* Z X)))
                        X (m/* 0.25 (m/+ X lambda))
                        Y (m/* 0.25 (m/+ Y lambda))
                        Z (m/* 0.25 (m/+ Z lambda))
                        A (m/* m/THIRD (m/+ X Y Z))
                        dx (m/- 1.0 (m// X A))
                        dy (m/- 1.0 (m// Y A))
                        dz (m/- 1.0 (m// Z A))]
                    (if (or (m/> i 25)
                            (m/<= (m/max (m/abs dx) (m/abs dy) (m/abs dz)) 1.0e-5))
                      (Vec4. dx dy dz A)
                      (recur X Y Z (m/inc i)))))
        E2 (m/+ (m/* (.x d) (.y d))
              (m/* (.y d) (.z d))
              (m/* (.z d) (.x d)))
        E22 (m/* E2 E2)
        E3 (m/* (.x d) (.y d) (.z d))]
    (-> (m/- 1.0 (m/* 0.1 E2))
        (m/+ (m/* 0.07142857142857142 E3))
        (m/+ (m/* 0.04166666666666666 E22))
        (m/- (m/* 0.06818181818181818 E2 E3))
        (m/- (m/* 0.02403846153846154 E22 E2))
        (m/+ (m/* 0.028846153846153848 E3 E3))
        (m/+ (m/* 0.0625 E22 E3))
        (m// (m/sqrt (.w d))))))

(defn- jacelk-nonz
  ^double [^double v]
  (if (m/< (m/abs v) m/EPSILON) (m/* m/EPSILON (m/sgn v)) v))

(defn jacelk
  ([] {:type :regular
       :config (fn [] {:k (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double k]}]
   (fn [^Vec2 v]
     (let [phi (.x v)
           psi (.y v)
           cotphi (m/cot phi)
           cotphi2 (m/* cotphi cotphi)
           b (m/* 0.5 (m/- (m/+ (m/* k (m/sq (m// (m/sinh psi) (m/+ m/EPSILON (m/sin phi)))))
                          cotphi2 k -1.0)))
           -b (m/- b)
           c (m/sqrt (m/abs (m/- (m/* b b) (m/* -1.0 (m/- 1.0 k) cotphi2))))
           X1 (m/max (m/+ -b c) (m/- -b c))
           mu (m/* (m/sgn psi)
                 (m/sqrt (m/abs (m// (m/dec (m// X1 (jacelk-nonz cotphi2))) (jacelk-nonz k)))))
           lambda (m/* (m/sgn phi) (m/sqrt (m/abs X1)))
           sina (m// (m/sgn lambda) (m/sqrt (m/inc (m/* lambda lambda))))
           cosa (m/* lambda sina)
           e-phi (m/* sina (jacelk-rf (m/* cosa cosa) (m/- 1.0 (m/* k sina sina)) 1.0))
           cosa (m// (m/sqrt (m/inc (m/* mu mu))))
           sina (m/* mu cosa)
           e-psi (m/* sina (jacelk-rf (m/* cosa cosa) (m/- 1.0 (m/* (m/- 1.0 k) sina sina)) 1.0))]
       (Vec2. (m/* amount e-phi) (m/* amount e-psi))))))


(defn jacsn
  ([] {:type :regular
       :config (fn [] {:k (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double k]}]
   (fn [^Vec2 v]
     (let [^JacobiData jac-x (u/jacobi-elliptic (.x v) k)
           ^JacobiData jac-y (u/jacobi-elliptic (.y v) (m/- 1.0 k))
           numx (m/* (.sn jac-x) (.dn jac-y))
           numy (m/* (.cn jac-x) (.dn jac-x) (.cn jac-y) (.sn jac-y))
           denom (m// amount (m/+ (m/* (.sn jac-x) (.sn jac-x) (.sn jac-y) (.sn jac-y) k)
                              (m/* (.cn jac-y) (.cn jac-y))
                              m/EPSILON))]
       (Vec2. (m/* denom numx) (m/* denom numy))))))

(defn joukowski
  ([] {:type :regular
       :config (fn [] {:p1 (u/sdrand 0.1 3.0)
                      :p2 (r/drand -3.0 3.0)
                      :inverse (r/brand)})})
  ([^double amount {:keys [^double p1 ^double p2 inverse]}]
   (let [z0 (Vec2. (m/- p2) 0.2)
         z1 (Vec2. (m/* p1 p1) 0.0)]
     (fn [^Vec2 z]
       (v/mult (if-not inverse
                 (let [zp (c/add z z0)]
                   (c/add zp (c/div z1 zp)))
                 (let [c2 (v/div z 2.0)
                       r (-> (c/mult z z)
                             (v/div 4.0)
                             (v/sub z1)
                             (c/sqrt))
                       r1 (-> (c/add c2 r)
                              (c/sub z0))
                       r2 (-> (c/sub c2 r)
                              (c/sub z0))]
                   (if (m/> (v/magsq r1) (v/magsq r2)) r1 r2))) amount)))))


(defn juliac
  "JuliaC"
  ([] {:type :random
       :config (fn [] {:re (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :im (r/drand -2.0 2.0)
                      :dist (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double re ^double im ^double dist]}]
   (let [im (m/* 0.01 im)
         rre (m// (m/+ m/EPSILON re))]
     (fn [^Vec2 v]
       (let [arg (m/+ (v/heading v)
                    (m/* m/TWO_PI (m/mod (r/lrand) re)))
             lnmod (m/* dist (m/log (v/magsq v)))
             a (m/+ (m/* arg rre)
                  (m/* lnmod im))
             mod2 (m/* amount (m/exp (m/- (m/* lnmod rre)
                                      (m/* arg im))))]
         (Vec2. (m/* mod2 (m/cos a))
                (m/* mod2 (m/sin a))))))))

(defn julia
  "Julia"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [a (m/+ (m/* 0.5 (v/heading v)) (m/* m/PI (r/lrand 2)))
           r (m/* amount (m/sqrt (v/mag v)))]
       (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a)))))))

(defn julian2
  ([] {:type :random
       :config (fn [] {:a (r/drand -1.5 1.5)
                      :b (r/drand -1.5 1.5)
                      :c (r/drand -1.5 1.5)
                      :d (r/drand -1.5 1.5)
                      :e (r/randval 0.0 (r/drand -1.5 1.5))
                      :f (r/randval 0.0 (r/drand -1.5 1.5))
                      :power (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :dist (u/sdrand 0.1 2.5)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d ^double e ^double f
                           ^double power ^double dist]}]
   (let [absn (m/long-max 1 (long (m/abs power)))
         cn (m/* 0.5 (m// dist power))]
     (fn [^Vec2 v]
       (let [x (m/+ (m/* a (.x v)) (m/* b (.y v)) e)
             y (m/+ (m/* c (.x v)) (m/* d (.y v)) f)
             angle (m/+ (v/heading v) (m/* m/TWO_PI (r/lrand absn)))
             r (m/* amount (m/pow (m/+ (m/* x x) (m/* y y)) cn))]
         (Vec2. (m/* r (m/cos angle))
                (m/* r (m/sin angle))))))))

(defn julian
  "JuliaN"
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :dist (u/sdrand 0.1 2.5)})})
  ([^double amount {:keys [^double power ^double dist]}]
   (let [abspower (m/max 1 (int (m/abs power)))
         cpower (m/* 0.5 (m// dist power))]
     (fn [v]
       (let [a (m// (m/+ (v/heading v) (m/* m/TWO_PI (double (r/irand abspower)))) power)
             r (m/* amount (m/pow (v/magsq v) cpower))]
         (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a))))))))

(defn juliaoutside
  ([] {:type :random
       :config (fn [] {:mode (r/irand 3)
                      :re-div (u/sdrand 0.1 2.0)
                      :im-div (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^long mode ^double re-div ^double im-div]}]
   (let [mode (m/long-mod mode 3)
         z3 (c/complex re-div im-div)
         me? (m/even? mode)
         ml2? (m/< mode 2)]
     (fn [^Vec2 v]
       (let [z (as-> v z
                 (if me? (c/sqrt z) z)
                 (c/add z c/ONE)
                 (if me? (c/sq z) z)
                 (c/div z (as-> v z2
                            (if me? (c/sqrt z2) z2)
                            (c/sub z2 c/ONE)
                            (if me? (c/sq z2) z2)))
                 (if ml2? (c/sqrt z) z)
                 (c/div z z3)
                 (c/scale z amount))]
         (if ml2? (r/randval z (c/neg z)) z))))))

(defn juliaq
  "juliaq by Zueuk, http://zueuk.deviantart.com/art/juliaq-Apophysis-plugins-340813357"
  ([] {:type :random
       :config (fn [] {:power (u/sirand 1 10)
                      :divisor (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double divisor ^double power]}]
   (let [inv-power (m// divisor power)
         half-inv-power (m/* 0.5 inv-power)
         inv-power-2pi (m// m/TWO_PI power)]
     (fn [^Vec2 v]
       (let [a (m/+ (m/* inv-power (v/heading v))
                  (m/* inv-power-2pi (double (r/irand))))
             r (m/* amount (m/pow (v/magsq v) half-inv-power))]
         (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a))))))))

(defn juliascope
  "JuliaScope"
  ([] {:type :random
       :config (fn [] {:power (r/randval (u/sirand 1 10) (u/sdrand 0.1 10.0))
                      :dist (u/sdrand 0.1 2.5)})})
  ([^double amount {:keys [^double power ^double dist]}]
   (let [abspower (m/max 1 (int (m/abs power)))
         cpower (m/* 0.5 (m// dist power))]
     (fn [v]
       (let [rnd (double (r/irand abspower))
             a (if (m/zero? (m/bit-and (long rnd) 1))
                 (m// (m/+ (m/* m/TWO_PI rnd) (v/heading v)) power)
                 (m// (m/- (m/* m/TWO_PI rnd) (v/heading v)) power))
             r (m/* amount (m/pow (v/magsq v) cpower))]
         (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a))))))))

;;

(defn julia2
  "Julia with different angle calc"
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [a (m/+ (m/* 0.5 (m/atan2 (.x v) (.y v))) (m/* m/PI (r/lrand 2)))
           r (m/* amount (m/sqrt (v/mag v)))]
       (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a)))))))

