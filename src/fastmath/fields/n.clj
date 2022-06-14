(ns fastmath.fields.n
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; TODO nblur

(defn npolar
  ([] {:type :random
       :config (fn [] {:parity (r/irand)
                      :n (u/sirand 1 10)})})
  ([^double amount {:keys [^long parity ^long n]}]
   (let [nnz (double (if (zero? n) 1.0 n))
         vvar (/ amount m/PI)
         vvar2 (* 0.5 vvar)
         absn (m/abs nnz)
         iabsn (int absn)
         cn (/ (* 2.0 nnz))
         isodd? (odd? parity)
         p (if isodd? (mod parity 6.0) 1.0)]
     (fn [^Vec2 v]
       (let [x (if isodd? (.x v) (* vvar (m/atan2 (.x v) (.y v))))
             y (if isodd? (.y v) (* vvar2 (m/log (v/magsq v))))
             angle (/ (+ (m/atan2 y x)
                         (* m/TWO_PI (mod (r/irand) iabsn))) nnz)
             r (* amount (m/pow (+ (* x x) (* y y)) cn) p)
             cosa (* r (m/cos angle))
             sina (* r (m/sin angle))]
         (Vec2. (if isodd? cosa (* vvar2 (m/log (+ (* sina sina) (* cosa cosa)))))
                (if isodd? sina (* vvar (m/atan2 cosa sina)))))))))

(defn ngon
  "Ngon"
  ([] {:type :regular
       :config (fn [] {:circle (r/drand -2.0 2.0)
                      :corners (r/drand -2.0 2.0)
                      :power (r/drand -10.0 10.0)
                      :sides (r/drand -10.0 10.0)})})
  ([^double amount {:keys [^double circle ^double corners ^double power ^double sides]}]
   (let [b (/ m/TWO_PI sides)
         hb (/ b 2.0)
         hpower (/ power 2.0)]
     (fn [v]
       (let [r-factor (m/pow (v/magsq v) hpower)
             theta (v/heading v)
             phi (- theta (* b (m/floor (/ theta b))))
             phi (if (> phi hb) (- phi b) phi)
             amp (/ (+ circle (* corners (dec (/ 1.0 (+ (m/cos phi) m/EPSILON))))) (+ r-factor m/EPSILON))]
         (v/mult v (* amount amp)))))))

(defn noise
  "Noise"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [a (r/drand m/TWO_PI)
           r (r/drand amount)]
       (Vec2. (* r (m/cos a))
              (* r (m/sin a)))))))

;;

