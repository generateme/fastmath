(ns fastmath.fields.n
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; TODO nblur

(defn npolar
  ([] {:type :random
       :config (fn [] {:parity (r/irand)
                      :n (u/sirand 1 10)})})
  ([^double amount {:keys [^long parity ^long n]}]
   (let [nnz (double (if (m/zero? n) 1.0 n))
         vvar (m// amount m/PI)
         vvar2 (m/* 0.5 vvar)
         absn (m/abs nnz)
         iabsn (int absn)
         cn (m// (m/* 2.0 nnz))
         isodd? (m/odd? parity)
         p (if isodd? (m/mod parity 6.0) 1.0)]
     (fn [^Vec2 v]
       (let [x (if isodd? (.x v) (m/* vvar (m/atan2 (.x v) (.y v))))
             y (if isodd? (.y v) (m/* vvar2 (m/log (v/magsq v))))
             angle (m// (m/+ (m/atan2 y x)
                         (m/* m/TWO_PI (r/lrand iabsn))) nnz)
             r (m/* amount (m/pow (m/+ (m/* x x) (m/* y y)) cn) p)
             cosa (m/* r (m/cos angle))
             sina (m/* r (m/sin angle))]
         (Vec2. (if isodd? cosa (m/* vvar2 (m/log (m/+ (m/* sina sina) (m/* cosa cosa)))))
                (if isodd? sina (m/* vvar (m/atan2 cosa sina)))))))))

(defn ngon
  "Ngon"
  ([] {:type :regular
       :config (fn [] {:circle (r/drand -2.0 2.0)
                      :corners (r/drand -2.0 2.0)
                      :power (r/drand -10.0 10.0)
                      :sides (r/drand -10.0 10.0)})})
  ([^double amount {:keys [^double circle ^double corners ^double power ^double sides]}]
   (let [b (m// m/TWO_PI sides)
         hb (m// b 2.0)
         hpower (m// power 2.0)]
     (fn [v]
       (let [r-factor (m/pow (v/magsq v) hpower)
             theta (v/heading v)
             phi (m/- theta (m/* b (m/floor (m// theta b))))
             phi (if (m/> phi hb) (m/- phi b) phi)
             amp (m// (m/+ circle (m/* corners (m/dec (m// 1.0 (m/+ (m/cos phi) m/EPSILON))))) (m/+ r-factor m/EPSILON))]
         (v/mult v (m/* amount amp)))))))

(defn noise
  "Noise"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [a (r/drand m/TWO_PI)
           r (r/drand amount)]
       (Vec2. (m/* r (m/cos a))
              (m/* r (m/sin a)))))))

;;

