(ns fastmath.fields.utils
  (:require [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec3]))

#_(set! *warn-on-reflection* true)

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def unitx (Vec2. 1.0 0.0))
(def zerov (Vec2. 0.0 0.0))

(deftype Pair [a b])
(deftype Triplet [a b c])
(deftype Tetrad [a b c d])

(defn sdrand
  "Symetric random from [-mx -mn] and [mn mx]"
  ^double  [^double mn ^double mx]
  (let [rand (r/drand mn mx)]
    (r/randval rand (m/* -1.0 rand))))

(defn sirand
  "Symmetric irand"
  (^long [^long mx] (sirand 1 mx))
  (^long [^long mn ^long mx]
   (let [rand (r/lrand mn mx)]
     (r/randval rand (m/* -1 rand)))))

(defn spow
  ^double [^double x ^double p]
  (if (m/neg? x)
    (m/- (m/pow (m/- x) p))
    (m/pow x p)))

;;

(def offsets (vec (for [x (range -1 2) y (range -1 2)] (Vec2. x y))))

(defn closest
  ^long [P ^long n U]
  (loop [i (long 0)
         d2min Double/MAX_VALUE
         j (long 0)]
    (if (m/< i n)
      (let [d2 (v/dist-sq (P i) U)
            low? (m/< d2 d2min)]
        (recur (m/inc i) (if low? d2 d2min) (if low? i j)))
      j)))

(defn vratio
  ^double [^Vec2 P ^Vec2 Q ^Vec2 U]
  (let [PmQ (v/sub P Q)]
    (if (v/is-zero? PmQ)
      1.0
      (-> (v/sub U Q)
          (v/emult PmQ)
          (v/sum)
          (m// (v/magsq PmQ))
          (m/* 2.0)))))

(defn voronoi
  ^double [P ^long n ^long q U]
  (let [Pq (P q)]
    (loop [i (long 0)
           ratiomax Double/MIN_VALUE]
      (if (m/< i n)
        (if (m/== i q)
          (recur (m/inc i) ratiomax)
          (recur (m/inc i) (m/max ratiomax (vratio (P i) Pq U))))
        ratiomax))))

;;

(defrecord JacobiData [^double sn ^double cn ^double dn])

(defn- jacobi-data [^Vec3 v] (JacobiData. (.x v) (.y v) (.z v)))

;; order: sn,cn,dn
;; emmc = 1-k^2 (m=k^2, m=1-emmc)
(defn jacobi-elliptic
  (^JacobiData [^double uu ^double emmc] (jacobi-elliptic uu emmc 8 0.003))
  (^JacobiData [^double uu ^double emmc ^long iters ^double accuracy]
   (jacobi-data (if-not (m/zero? emmc)
                  (let [bo? (m/neg? emmc)
                        ^Vec3 emcud (if bo?
                                      (let [d (m/- 1.0 emmc)
                                            emc (m// (m/- emmc) d)
                                            d (m/sqrt d)]
                                        (Vec3. emc (m/* d uu) d))
                                      (Vec3. emmc uu 0.0))
                        em (double-array iters)
                        en (double-array iters)
                        ^Vec3 emccl (loop [a 1.0
                                           l (long 0)
                                           emc (.x emcud)]
                                      (if (m/< l iters)
                                        (do
                                          (aset ^doubles em l a)
                                          (let [emc (m/sqrt emc)]
                                            (aset ^doubles en l emc)
                                            (let [c (m/* 0.5 (m/+ a emc))]
                                              (if (m/< (m/abs (m/- a emc)) (m/* accuracy a))
                                                (Vec3. emc c l)
                                                (recur c (m/inc l) (m/* a emc))))))
                                        (Vec3. emc a (m/dec iters))))
                        u (m/* (.y emccl) (.y emcud))
                        sn (m/sin u)
                        cn (m/cos u)
                        ^Vec3 scd (if-not (m/zero? sn)
                                    (loop [dn 1.0
                                           i (long (.z emccl))
                                           a (m// cn sn)
                                           c (m/* (.y emccl) a)]
                                      (if-not (m/neg? i)
                                        (let [b (aget ^doubles em i)
                                              a (m/* c a)
                                              c (m/* dn c)]
                                          (recur (m// (m/+ (aget ^doubles en i) a)
                                                    (m/+ a b)) (m/dec i) (m// c b) c))
                                        (let [a (m// (m/sqrt (m/inc (m/* c c))))
                                              sn (if (m/neg? sn) (m/- a) a)]
                                          (Vec3. sn (m/* c sn) dn))))
                                    (Vec3. sn cn 1.0))]
                    (if bo?
                      (Vec3. (m// (.x scd) (.z emcud)) (.z scd) (.y scd))
                      scd))
                  (let [cn (m// (m/cosh uu))]
                    (Vec3. (m/tanh uu) cn cn))))))

;; Local noise generating functions (for various noise implementations)

(defn make-noise-variation
  "Calculate shift by angle taken from noise"
  ([^double amount ^double scale noise-fn] (make-noise-variation amount scale 1.0 noise-fn))
  ([^double amount ^double scale ^double m noise-fn]
   (let [m (m/* m/TWO_PI m)]
     (fn [v]
       (let [^Vec2 vv (v/mult v scale)
             angle (m/* m ^double (noise-fn (.x vv) (.y vv)))]
         (v/add v (Vec2. (m/* amount (m/cos angle))
                         (m/* amount (m/sin angle)))))))))

(defn make-noise-variation2
  "Calculate shift by vector taken from noise"
  [^double amount scale noise-fn]
  (fn [v]
    (let [^Vec2 vv (v/mult v scale)
          x1 (m/- ^double (noise-fn (.x vv) (.y vv)) 0.5)
          y1 (m/- ^double (noise-fn (.y vv) m/E (.x vv)) 0.5)]
      (v/add v (Vec2. (m/* amount x1)
                      (m/* amount y1))))))

