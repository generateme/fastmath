(ns fastmath.interpolation.sprague
  (:require [fastmath.core :as m]
            [fastmath.polynomials :as poly]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- fix-search
  ^long [^long id]
  (m/dec (if (m/neg? id) (m/dec (m/abs id)) id)))

(def ^:private coeffs0 (double-array [4.229665071770334 -9.37799043062201 14.511961722488037
                                    -12.669856459330143 5.167464114832535 -0.8612440191387559]))
(def ^:private coeffs1 (double-array [2.430622009569378 -2.5837320574162677 2.3349282296650715
                                    -1.755980861244019 0.6889952153110047 -0.11483253588516745]))
(def ^:private coeffs2 (double-array [-0.11483253588516745 0.6889952153110047 -1.755980861244019
                                    2.3349282296650715 -2.5837320574162677 2.430622009569378]))
(def ^:private coeffs3 (double-array [-0.8612440191387559 5.167464114832535 -12.669856459330143
                                    14.511961722488037 -9.37799043062201 4.229665071770334]))

(def ^:private w0 (double-array [0.08333333333333333 -0.6666666666666666 0.0
                               0.6666666666666666 -0.08333333333333333 0.0]))
(def ^:private w1 (double-array [-0.041666666666666664 0.6666666666666666 -1.25
                               0.6666666666666666 -0.041666666666666664 0.0]))
(def ^:private w2 (double-array [-0.375 1.625 -2.9166666666666665
                               2.75 -1.375 0.29166666666666663]))
(def ^:private w3 (double-array [0.5416666666666666 -2.6666666666666665 5.25
                               -5.166666666666666 2.5416666666666665 -0.5]))
(def ^:private w4 (double-array [-0.20833333333333331 1.0416666666666665 -2.083333333333333
                               2.083333333333333 -1.0416666666666665 0.20833333333333331]))

(defn sprague
  [xs ys]
  (let [ax (m/seq->double-array xs)
        len (alength ax)
        nlen (m/+ len 4)
        last (m/dec len)
        x0 (Array/aget ax 0)
        diff0 (m/- (Array/aget ax 1) x0)
        x1 (Array/aget ax last)
        diff1 (m/- x1 (Array/aget ax (m/dec last)))
        ^doubles nx (double-array nlen)
        
        ay (m/seq->double-array ys)
        ytemp (double-array 6)
        ny (double-array nlen)]
    
    (Array/aset nx 0 (m/- x0 (m/* 2.0 diff0)))
    (Array/aset nx 1 (m/- x0 diff0))
    (Array/aset nx (m/+ 2 len) (m/+ x1 diff1))
    (Array/aset nx (m/+ 3 len) (m/+ x1 (m/* 2.0 diff1)))
    (System/arraycopy ax 0 nx 2 len)

    (System/arraycopy ay 0 ytemp 0 6)
    (Array/aset ny 0 (v/dot ytemp coeffs0))
    (Array/aset ny 1 (v/dot ytemp coeffs1))
    (System/arraycopy ay (m/- len 6 ) ytemp 0 6)
    (Array/aset ny (m/+ 2 len) (v/dot ytemp coeffs2))
    (Array/aset ny (m/+ 3 len) (v/dot ytemp coeffs3))
    (System/arraycopy ay 0 ny 2 len)
    
    (fn ^double [^double x]
      (let [i (m/constrain (fix-search (java.util.Arrays/binarySearch nx x)) 2 len)
            xx (m// (m/- x (Array/aget nx i))
                    (m/- (Array/aget nx (m/inc i))
                         (Array/aget nx i)))
            ytemp (double-array 6)]
        (System/arraycopy ny (m/- i 2) ytemp 0 6)
        (poly/mevalpoly xx (Array/aget ny i)
                        (v/dot w0 ytemp)
                        (v/dot w1 ytemp)
                        (v/dot w2 ytemp)
                        (v/dot w3 ytemp)
                        (v/dot w4 ytemp))))))

