(ns fastmath.interpolation
  (:require [fastmath.interpolation.acm :as acm]
            [fastmath.interpolation.ssj :as ssj]
            [fastmath.interpolation.linear :as linear]
            [fastmath.interpolation.cubic :as cubic]
            [fastmath.interpolation.barycentric :as bc]
            [fastmath.interpolation.shepard :as shepard]
            [fastmath.interpolation.rbf :as rbf]
            [fastmath.interpolation.kriging :as kriging]
            [fastmath.interpolation.gp :as gp]
            [fastmath.interpolation.step :as step]
            [fastmath.interpolation.monotone :as monotone]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmulti interpolation (fn [interpolation-name & _] interpolation-name))

;; 1d

(defmethod interpolation :linear [_ xs ys] (linear/linear xs ys))
(defmethod interpolation :cubic [_ xs ys] (cubic/cubic xs ys))
(defmethod interpolation :monotone [_ xs ys] (monotone/monotone xs ys))
(defmethod interpolation :akima [_ xs ys] (acm/akima xs ys))
(defmethod interpolation :neville [_ xs ys] (acm/neville xs ys))
(defmethod interpolation :divided-difference [_ xs ys] (acm/divided-difference xs ys))
(defmethod interpolation :polynomial [_ xs ys] (ssj/polynomial xs ys))

(defmethod interpolation :step-before [_ xs ys] (step/step-before xs ys))
(defmethod interpolation :step-after [_ xs ys] (step/step-after xs ys))
(defmethod interpolation :step
  ([_ xs ys] (step/step xs ys))
  ([_ xs ys params] (step/step xs ys params)))

(defmethod interpolation :b-spline
  ([_ xs ys] (ssj/b-spline xs ys))
  ([_ xs ys params] (ssj/b-spline xs ys params)))

(defmethod interpolation :barycentric
  ([_ xs ys] (bc/barycentric xs ys))
  ([_ xs ys params] (bc/barycentric xs ys params)))

(defmethod interpolation :loess
  ([_ xs ys] (acm/loess xs ys))
  ([_ xs ys params] (acm/loess xs ys params)))

(defmethod interpolation :cubic-smoothing
  ([_ xs ys] (ssj/cubic-smoothing xs ys))
  ([_ xs ys params] (ssj/cubic-smoothing xs ys params)))

;; 2d grid

(defmethod interpolation :bilinear [_ xs ys vss] (linear/bilinear xs ys vss))
(defmethod interpolation :bicubic [_ xs ys vss] (acm/bicubic xs ys vss))
(defmethod interpolation :cubic-2d [_ xs ys vss] (cubic/cubic-2d xs ys vss))

;; multidim and kernel

(defmethod interpolation :microsphere-projection
  ([_ xss ys] (acm/microsphere-projection xss ys))
  ([_ xss ys params] (acm/microsphere-projection xss ys params)))

(defmethod interpolation :shepard
  ([_ xss ys] (shepard/shepard xss ys))
  ([_ xss ys params] (shepard/shepard xss ys params)))

(defmethod interpolation :rbf
  ([_ xss ys] (rbf/rbf xss ys))
  ([_ xss ys kernel] (rbf/rbf xss ys kernel))
  ([_ xss ys kernel params] (rbf/rbf xss ys kernel params)))

(defmethod interpolation :kriging
  ([_ xss ys] (kriging/kriging xss ys))
  ([_ xss ys variogram] (kriging/kriging xss ys variogram))
  ([_ xss ys variogram params] (kriging/kriging xss ys variogram params)))

(defmethod interpolation :gp
  ([_ xss ys] (gp/gp xss ys))
  ([_ xss ys kernel] (gp/gp xss ys kernel))
  ([_ xss ys kernel params] (gp/gp xss ys kernel params)))
