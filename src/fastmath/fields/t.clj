(ns fastmath.fields.t
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c]
            [fastmath.interpolation :as i])
  (:import [fastmath.vector Vec2]))

(set! *warn-on-reflection* true)

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

(defn tangent
  "Tangent"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [d (m/cos (.y v))
           id (/ 1.0 (if (zero? d) m/EPSILON d))]
       (Vec2. (* amount (m/sin (.x v)) id)
              (* amount (m/tan (.y v))))))))

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
       (let [[^double cc1 ^double cc2 ^double fr ^double rr] (if (pos? (.x v))
                                                               [c1 (- c2) (/ r2 r1) r1]
                                                               [(- c2) c1 (/ r1 r2) r2])
             nv (Vec2. (- cc1 (.x v)) (.y v))
             rm (v/mag nv)
             r (* rm fr)
             a (v/heading nv)
             res (Vec2. (+ cc2 (* r (m/cos a)))
                        (* r (m/sin a)))]
         (if (<= rm rr)
           (v/mult res amount)
           (v/mult v amount)))))))

