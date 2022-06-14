(ns fastmath.fields.utils
  (:require [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec3]))

#_(set! *warn-on-reflection* true)

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def unitx (Vec2. 1.0 0.0))
(def zerov (Vec2. 0.0 0.0))

(deftype Pair [a b])
(deftype Triplet [a b c])
(deftype Tetrad [a b c d])

(defn sdrand
  "Symetric random from [-mx -mn] and [mn mx]"
  ^double  [^double mn ^double mx]
  (let [rand (r/drand mn mx)]
    (r/randval rand (* -1.0 rand))))

(defn sirand
  "Symmetric irand"
  (^long [^long mx] (sirand 1 mx))
  (^long [^long mn ^long mx]
   (let [rand (r/irand mn mx)]
     (r/randval rand (* -1 rand)))))

(defn spow
  ^double [^double x ^double p]
  (if (neg? x)
    (- (m/pow (- x) p))
    (m/pow x p)))

;;

(def offsets (vec (for [x (range -1 2) y (range -1 2)] (Vec2. x y))))

(defn closest
  ^long [P ^long n U]
  (loop [i (long 0)
         d2min Double/MAX_VALUE
         j (long 0)]
    (if (< i n)
      (let [d2 (v/dist-sq (P i) U)
            low? (< d2 d2min)]
        (recur (inc i) (if low? d2 d2min) (if low? i j)))
      j)))

(defn vratio
  ^double [^Vec2 P ^Vec2 Q ^Vec2 U]
  (let [PmQ (v/sub P Q)]
    (if (v/is-zero? PmQ)
      1.0
      (-> (v/sub U Q)
          (v/emult PmQ)
          (v/sum)
          (/ (v/magsq PmQ))
          (* 2.0)))))

(defn voronoi
  ^double [P ^long n ^long q U]
  (let [Pq (P q)]
    (loop [i (long 0)
           ratiomax Double/MIN_VALUE]
      (if (< i n)
        (if (== i q)
          (recur (inc i) ratiomax)
          (recur (inc i) (max ratiomax (vratio (P i) Pq U))))
        ratiomax))))

;;

(defrecord JacobiData [^double sn ^double cn ^double dn])

(defn- jacobi-data [^Vec3 v] (JacobiData. (.x v) (.y v) (.z v)))

;; order: sn,cn,dn
;; emmc = 1-k^2 (m=k^2, m=1-emmc)
(defn jacobi-elliptic
  (^JacobiData [^double uu ^double emmc] (jacobi-elliptic uu emmc 8 0.003))
  (^JacobiData [^double uu ^double emmc ^long iters ^double accuracy]
   (jacobi-data (if-not (zero? emmc)
                  (let [bo? (neg? emmc)
                        ^Vec3 emcud (if bo?
                                      (let [d (- 1.0 emmc)
                                            emc (/ (- emmc) d)
                                            d (m/sqrt d)]
                                        (Vec3. emc (* d uu) d))
                                      (Vec3. emmc uu 0.0))
                        em (double-array iters)
                        en (double-array iters)
                        ^Vec3 emccl (loop [a 1.0
                                           l (long 0)
                                           emc (.x emcud)]
                                      (if (< l iters)
                                        (do
                                          (aset ^doubles em l a)
                                          (let [emc (m/sqrt emc)]
                                            (aset ^doubles en l emc)
                                            (let [c (* 0.5 (+ a emc))]
                                              (if (< (m/abs (- a emc)) (* accuracy a))
                                                (Vec3. emc c l)
                                                (recur c (inc l) (* a emc))))))
                                        (Vec3. emc a (dec iters))))
                        u (* (.y emccl) (.y emcud))
                        sn (m/sin u)
                        cn (m/cos u)
                        ^Vec3 scd (if-not (zero? sn)
                                    (loop [dn 1.0
                                           i (long (.z emccl))
                                           a (/ cn sn)
                                           c (* (.y emccl) a)]
                                      (if-not (neg? i)
                                        (let [b (aget ^doubles em i)
                                              a (* c a)
                                              c (* dn c)]
                                          (recur (/ (+ (aget ^doubles en i) a)
                                                    (+ a b)) (dec i) (/ c b) c))
                                        (let [a (/ (m/sqrt (inc (* c c))))
                                              sn (if (neg? sn) (- a) a)]
                                          (Vec3. sn (* c sn) dn))))
                                    (Vec3. sn cn 1.0))]
                    (if bo?
                      (Vec3. (/ (.x scd) (.z emcud)) (.z scd) (.y scd))
                      scd))
                  (let [cn (/ (m/cosh uu))]
                    (Vec3. (m/tanh uu) cn cn))))))

;; Local noise generating functions (for various noise implementations)

(defn make-noise-variation
  "Calculate shift by angle taken from noise"
  ([^double amount ^double scale noise-fn] (make-noise-variation amount scale 1.0 noise-fn))
  ([^double amount ^double scale ^double m noise-fn]
   (let [m (* m/TWO_PI m)]
     (fn [v]
       (let [^Vec2 vv (v/mult v scale)
             angle (* m ^double (noise-fn (.x vv) (.y vv)))]
         (v/add v (Vec2. (* amount (m/cos angle))
                         (* amount (m/sin angle)))))))))

(defn make-noise-variation2
  "Calculate shift by vector taken from noise"
  [^double amount scale noise-fn]
  (fn [v]
    (let [^Vec2 vv (v/mult v scale)
          x1 (- ^double (noise-fn (.x vv) (.y vv)) 0.5)
          y1 (- ^double (noise-fn (.y vv) m/E (.x vv)) 0.5)]
      (v/add v (Vec2. (* amount x1)
                      (* amount y1))))))
