(ns fastmath.rbf
  "Radial Basis Function

  Create with multifunction [[rbf]].

  All of them accept scaling factor `scale`.
  Only polyharmonic is defined with integer exponent `k`.
  See [[rbfs-list]] for all names.

  [[rbf-obj]] returns SMILE library object for defined function."
  {:metadoc/categories {:rbf "Radial Basis Function"}}
  (:require [fastmath.core :as m])
  (:import [smile.math.rbf RadialBasisFunction]))

;; RBF
(defmulti rbf
  "Create Radial Basis Function

  Optional parameter `scale`.

  `:polyharmonic` has also exponent `k` parameter."
  {:metadoc/categories #{:rbf}}
  (fn [n & _] n))

(defmethod rbf :linear
  ([_] (fn ^double [^double x] x))
  ([_ ^double scale] (fn ^double [^double x] (* scale x))))

(defmethod rbf :gaussian
  ([_] (fn ^double [^double x] (m/exp (- (* x x)))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (m/exp (- (* s2 x x)))))))

(defmethod rbf :multiquadratic
  ([_] (fn ^double [^double x] (m/sqrt (inc (* x x)))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (m/sqrt (inc (* s2 x x)))))))

(defmethod rbf :inverse-multiquadratic
  ([_] (fn ^double [^double x] (/ (m/sqrt (inc (* x x))))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (/ (m/sqrt (inc (* s2 x x))))))))

(defmethod rbf :inverse-quadratic
  ([_] (fn ^double [^double x] (/ (inc (* x x)))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (/ (inc (* s2 x x)))))))

(defmethod rbf :polyharmonic
  ([_ ^long k] (rbf :polyharmonic k 1.0))
  ([_ ^long k ^double scale] (if (even? k)
                               (fn ^double [^double x]
                                 (if (pos? x)
                                   (let [sx (* scale x)] 
                                     (* (m/log sx) (m/pow sx k)))
                                   0.0))
                               (fn ^double [^double x] (m/pow (* scale x) k)))))


(defmethod rbf :thinplate
  ([_] (fn ^double [^double x] (if (pos? x)
                                 (* x x (m/log x))
                                 0.0)))
  ([_ ^double scale]
   (let [s2 (* scale scale)]
     (fn ^double [^double x] (if (pos? x)
                               (* s2 x x (m/log (* x scale)))
                               0.0)))))

(defmethod rbf :wendland
  ([_] (rbf :wendland 4.0))
  ([_ ^double scale] (let [rr (/ scale)]
                       (fn ^double [^double x]
                         (let [xrr (* x rr)]
                           (if (< xrr 1.0)
                             (* (m/pow (- 1.0 xrr) 4.0)
                                (inc (* 4.0 xrr)))
                             0.0))))))

(defmethod rbf :wu
  ([_] (rbf :wu 2.0))
  ([_ ^double scale] (let [rr (/ scale)]
                       (fn ^double [^double x]
                         (let [xrr (* x rr)]
                           (if (< xrr 1.0)
                             (let [xrr2 (* xrr xrr)
                                   xrr3 (* xrr xrr2)]
                               (* (m/pow (- 1.0 xrr) 4.0)
                                  (+ 4.0 (* 16.0 xrr) (* 12.0 xrr2) (* 3.0 xrr3))))
                             0.0))))))


(defn rbf-obj
  "Create RBF Smile object.

  Used to pass to Smile constructors/functions."
  {:metadoc/categories #{:rbf}}
  [rbf-fn]
  (reify RadialBasisFunction
    (^double f [_ ^double x] (rbf-fn x))))

(def ^{:doc "Radial Basis function names"
       :metadoc/categories #{:rbf}}
  rbfs-list (keys (methods rbf)))

