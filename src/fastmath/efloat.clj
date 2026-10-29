(ns fastmath.efloat
  "(re)Implementation of EFloat/Interval from pbrt-v3/pbrt-v4.

  A floating point number structure which keeps a track of error caused by operations."
  (:refer-clojure :exclude [min max abs])
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defrecord EFloat [^double v ^double low ^double high])

;; helpers

(defn- add+ ^double [^double a ^double b] (m/next-double (m/+ a b)))
(defn- add- ^double [^double a ^double b] (m/prev-double (m/+ a b)))
(defn- sub+ ^double [^double a ^double b] (m/next-double (m/- a b)))
(defn- sub- ^double [^double a ^double b] (m/prev-double (m/- a b)))
(defn- mul+ ^double [^double a ^double b] (m/next-double (m/* a b)))
(defn- mul- ^double [^double a ^double b] (m/prev-double (m/* a b)))
(defn- div+ ^double [^double a ^double b] (m/next-double (m// a b)))
(defn- div- ^double [^double a ^double b] (m/prev-double (m// a b)))
(defn- sqrt+ ^double [^double a] (m/next-double (m/sqrt a)))
(defn- sqrt- ^double [^double a] (m/prev-double (m/sqrt a)))
(defn- fma+ ^double [^double a ^double b ^double c] (m/next-double (m/fma a b c)))
(defn- fma- ^double [^double a ^double b ^double c] (m/prev-double (m/fma a b c)))

(defn efloat
  "Create EFloat object from a single value or low and high values."
  (^EFloat [^double v] (EFloat. v v v))
  (^EFloat [^double v ^double err]
   (if (m/zero? err)
     (EFloat. v v v)
     (EFloat. v (sub- v err) (add+ v err))))
  (^EFloat [^double v ^double low ^double high] (EFloat. v low high)))

(defn upper-bound ^double [^EFloat ev] (.high ev))
(defn lower-bound ^double [^EFloat ev] (.low ev))
(defn mid-point ^double [^EFloat ev] (m/* 0.5 (m/+ (.low ev) (.high ev))))
(defn width ^double [^EFloat ev] (m/- (.high ev) (.low ev)))
(defn ->double ^double [^EFloat ev] (.v ev))
(defn absolute-error ^double [^EFloat ev] (m/next-double (m/max (m/abs (m/- (.high ev) (.v ev)))
                                                             (m/abs (m/- (.v ev) (.low ev))))))
(defn relative-error ^double [^EFloat ev] (m/next-double (m// (m/max (m/abs (m/- (.high ev) (.v ev)))
                                                                (m/abs (m/- (.v ev) (.low ev))))
                                                         (m/abs (.v ev)))))

(defn equals?
  ^EFloat [^EFloat ev ^double v]
  (and (m/== v (.low ev))
       (m/== v (.high ev))))

(defn in-range?
  [^EFloat ev ^double v]
  (and (m/>= v (.low ev))
       (m/<= v (.high ev))))

(defn neg ^EFloat [^EFloat ev] (EFloat. (m/- (.v ev)) (m/- (.high ev)) (m/- (.low ev))))

(defn add ^EFloat [^EFloat ev1 ^EFloat ev2]
  (EFloat. (m/+ (.v ev1) (.v ev2)) (add- (.low ev1) (.low ev2)) (add+ (.high ev1) (.high ev2))))

(defn sub ^EFloat [^EFloat ev1 ^EFloat ev2]
  (EFloat. (m/- (.v ev1) (.v ev2)) (sub- (.low ev1) (.low ev2)) (sub+ (.high ev1) (.high ev2))))

(defn mul ^EFloat [^EFloat ev1 ^EFloat ev2]
  (let [m1 (m/* (.low ev1) (.low ev2))
        m2 (m/* (.high ev1) (.low ev2))
        m3 (m/* (.low ev1) (.high ev2))
        m4 (m/* (.high ev1) (.high ev2))]
    (EFloat. (m/* (.v ev1) (.v ev2))
             (m/min (m/prev-double m1) (m/prev-double m2) (m/prev-double m3) (m/prev-double m4))
             (m/max (m/next-double m1) (m/next-double m2) (m/next-double m3) (m/next-double m4)))))

(defn div ^EFloat [^EFloat ev1 ^EFloat ev2]
  (if (in-range? ev2 0.0)
    (EFloat. (m// (.v ev1) 0.0) ##-Inf ##Inf)
    (let [m1 (m// (.low ev1) (.low ev2))
          m2 (m// (.high ev1) (.low ev2))
          m3 (m// (.low ev1) (.high ev2))
          m4 (m// (.high ev1) (.high ev2))]
      (EFloat. (m// (.v ev1) (.v ev2))
               (m/min (m/prev-double m1) (m/prev-double m2) (m/prev-double m3) (m/prev-double m4))
               (m/max (m/next-double m1) (m/next-double m2) (m/next-double m3) (m/next-double m4))))))

(defn sq ^EFloat [^EFloat ev]
  (let [alow (m/abs (.low ev))
        ahigh (m/abs (.high ev))
        alow' (if (m/> alow ahigh) ahigh alow)
        ahigh' (if (m/> alow ahigh) alow ahigh)]
    (if (in-range? ev 0.0)
      (EFloat. (m/* (.v ev) (.v ev)) 0.0 (mul+ ahigh' ahigh'))
      (EFloat. (m/* (.v ev) (.v ev)) (mul- alow' alow') (mul+ ahigh' ahigh')))))
(defn sqrt ^EFloat [^EFloat ev] (EFloat. (m/sqrt (.v ev)) (sqrt- (.low ev)) (sqrt+ (.high ev))))

(defn fma ^EFloat [^EFloat ev1 ^EFloat ev2 ^EFloat ev3]
  (let [low (m/min (fma- (.low ev1) (.low ev2) (.low ev3))
                   (fma- (.high ev1) (.low ev2) (.low ev3))
                   (fma- (.low ev1) (.high ev2) (.low ev3))
                   (fma- (.high ev1) (.high ev2) (.low ev3)))
        high (m/max (fma+ (.low ev1) (.low ev2) (.high ev3))
                    (fma+ (.high ev1) (.low ev2) (.high ev3))
                    (fma+ (.low ev1) (.high ev2) (.high ev3))
                    (fma+ (.high ev1) (.high ev2) (.high ev3)))]
    (EFloat. (m/fma (.v ev1) (.v ev2) (.v ev3)) low high)))

(defn difference-of-products
  ^EFloat [^EFloat a ^EFloat b ^EFloat c ^EFloat d]
  (let [ab0 (m/* (.low a) (.low b))
        ab1 (m/* (.high a) (.low b))
        ab2 (m/* (.low a) (.high b))
        ab3 (m/* (.high a) (.high b))
        ab-low (m/min ab0 ab1 ab2 ab3)
        ab-high (m/max ab0 ab2 ab2 ab3)
        a-low (if (or (m/== ab-low ab0) (m/== ab-low ab2)) (.low a) (.high a))
        b-low (if (or (m/== ab-low ab0) (m/== ab-low ab1)) (.low b) (.high b))
        a-high (if (or (m/== ab-high ab0) (m/== ab-high ab2)) (.low a) (.high a))
        b-high (if (or (m/== ab-high ab0) (m/== ab-high ab1)) (.low b) (.high b))
        
        cd0 (m/* (.low c) (.low d))
        cd1 (m/* (.high c) (.low d))
        cd2 (m/* (.low c) (.high d))
        cd3 (m/* (.high c) (.high d))
        cd-low (m/min cd0 cd1 cd2 cd3)
        cd-high (m/max cd0 cd2 cd2 cd3)
        c-low (if (or (m/== cd-low cd0) (m/== cd-low cd2)) (.low c) (.high c))
        d-low (if (or (m/== cd-low cd0) (m/== cd-low cd1)) (.low d) (.high d))
        c-high (if (or (m/== cd-high cd0) (m/== cd-high cd2)) (.low c) (.high c))
        d-high (if (or (m/== cd-high cd0) (m/== cd-high cd1)) (.low d) (.high d))

        low (m/difference-of-products a-low b-low c-high d-high)
        high (m/difference-of-products a-high b-high c-low d-low)]
    (EFloat. (m/difference-of-products (.v a) (.v b) (.v c) (.v d))
             (m/prev-double (m/prev-double low))
             (m/next-double (m/next-double high)))))

(defn abs ^EFloat [^EFloat ev]
  (cond
    (m/>= (.low ev) 0.0) ev
    (m/<= (.high ev) 0.0) (neg ev)
    :else (EFloat. (m/abs (.v ev)) 0.0 (m/max (m/- (.low ev)) (.high ev)))))

(defn acos ^EFloat [^EFloat ev]
  (let [low (m/acos (m/min 1.0 (.high ev)))
        high (m/acos (m/max -1.0 (.low ev)))]
    (EFloat. (m/acos (.v ev)) (m/max 0.0 (m/prev-double low)) (m/next-double high))))

(def ^{:const true :private true :tag 'double} PI32 (m/* 1.5 m/PI))

(defn sin ^EFloat [^EFloat ev]
  (let [low (m/sin (m/max 0.0 (.low ev)))
        high (m/sin (.high ev))
        low (if (m/> low high) high low)
        high (if (m/> low high) low high)
        low (m/max -1.0 (m/prev-double low))
        high (m/min 1.0 (m/next-double high))
        s (m/sin (.v ev))]
    (cond
      (in-range? ev m/HALF_PI) (EFloat. s low 1.0)
      (in-range? ev PI32) (EFloat. s -1.0 high)
      :else (EFloat. s low high))))

(defn cos ^EFloat [^EFloat ev]
  (let [low (m/cos (m/max 0.0 (.low ev)))
        high (m/cos (.high ev))
        low (if (m/> low high) high low)
        high (if (m/> low high) low high)
        low (m/max -1.0 (m/prev-double low))
        high (m/min 1.0 (m/next-double high))]
    (EFloat. (m/cos (.v ev)) (if (in-range? ev m/PI) -1.0 low) high)))

(defn sum-of-products
  ^EFloat [^EFloat a ^EFloat b ^EFloat c ^EFloat d]
  (difference-of-products a b (neg c) d))

(defn addf ^EFloat [^EFloat ev ^double v] (add ev (EFloat. v v v)))
(defn subf ^EFloat [^EFloat ev ^double v] (sub ev (EFloat. v v v)))
(defn mulf ^EFloat [^EFloat ev ^double v]
  (if (m/pos? v)
    (EFloat. (m/* (.v ev) v) (mul- (.low ev) v) (mul+ (.high ev) v))
    (EFloat. (m/* (.v ev) v) (mul- (.high ev) v) (mul+ (.low ev) v))))
(defn divf ^EFloat [^EFloat ev ^double v]
  (if (m/zero? v)
    (EFloat. (m// (.v ev) v) ##-Inf ##Inf)
    (if (m/pos? v)
      (EFloat. (m// (.v ev) v) (div- (.low ev) v) (div+ (.high ev) v))
      (EFloat. (m// (.v ev) v) (div- (.high ev) v) (div+ (.low ev) v)))))

(defn floor ^double [^EFloat ev] (m/floor (.low ev)))
(defn ceil ^double [^EFloat ev] (m/ceil (.high ev)))

(defn min ^double [^EFloat ev1 ^EFloat ev2] (m/min (.low ev1) (.low ev2)))
(defn max ^double [^EFloat ev1 ^EFloat ev2] (m/max (.low ev1) (.low ev2)))

;;

(defn- mulpow2
  ^EFloat [^EFloat ev ^double s]
  (let [low (m/* s (.low ev))
        high (m/* s (.high ev))]
    (EFloat. (m/* s (.v ev)) (m/min low high) (m/max low high))))

(defrecord Pair [a b])

(defn quadratic
  ^Pair [^EFloat a ^EFloat b ^EFloat c]
  (let [discrim (difference-of-products b b (mulpow2 a 4.0) c)]
    (when-not (m/neg? (.low discrim))
      (let [float-root-discrim (sqrt discrim)
            q (if (m/neg? (->double b))
                (mulpow2 (sub b float-root-discrim) -0.5)
                (mulpow2 (add b float-root-discrim) -0.5))
            t0 (div q a)
            t1 (div c q)]
        (if (m/> (.low t0) (.low t1))
          (Pair. t1 t0) (Pair. t0 t1))))))

;; https://github.com/mmp/pbrt-v3/blob/master/src/core/efloat.h

;; pbrt source code is Copyright(c) 1998-2016
;; Matt Pharr, Greg Humphreys, and Wenzel Jakob.
;; This file is part of pbrt.
;; Redistribution and use in source and binary forms, with or without
;; modification, are permitted provided that the following conditions are
;; met:
;; - Redistributions of source code must retain the above copyright
;; notice, this list of conditions and the following disclaimer.
;; - Redistributions in binary form must reproduce the above copyright
;; notice, this list of conditions and the following disclaimer in the
;; documentation and/or other materials provided with the distribution.
;; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
;;     IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
;; TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
;; PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
;; HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

;; https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/math.h

;; pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
;; The pbrt source code is licensed under the Apache License, Version 2.0.
;; SPDX: Apache-2.0
