(ns fastmath.kernel.mercer
  "Kernel functions."
  (:require [fastmath.core :as m])
  (:import [clojure.lang IFn]
           [smile.math.kernel MercerKernel GaussianKernel HellingerKernel HyperbolicTangentKernel LaplacianKernel LinearKernel PearsonKernel, PolynomialKernel ThinPlateSplineKernel]))

(set! *warn-on-reflection* true)

(defn- gen-smile
  [^MercerKernel dst]
  (reify
    IFn (invoke [_ x y] (.k dst (double-array x) (double-array y)))
    MercerKernel (k [_ x y] (.k dst x y))))

(defmulti kernel "Create kernel object."
  (fn [k & r] k))

(defmethod kernel :gaussian [_ & [sigma]] (gen-smile (GaussianKernel. sigma)))
(defmethod kernel :hellinger [_ & _] (gen-smile (HellingerKernel.)))
(defmethod kernel :hyperbolic-tangent [_ & [scale offset]]
  (if-not (or scale offset)
    (gen-smile (HyperbolicTangentKernel.))
    (gen-smile (HyperbolicTangentKernel. scale offset))))
(defmethod kernel :laplacian [_ & [sigma]] (gen-smile (LaplacianKernel. sigma)))
(defmethod kernel :linear [_ & _] (gen-smile (LinearKernel.)))
(defmethod kernel :pearson [_ & [omega sigma]]
  (if-not (or omega sigma)
    (gen-smile (PearsonKernel.))
    (gen-smile (PearsonKernel. omega sigma))))
(defmethod kernel :polynomial [_ & [degree scale offset]]
  (if-not (or scale offset)
    (gen-smile (PolynomialKernel. degree))
    (gen-smile (PolynomialKernel. degree scale offset))))
(defmethod kernel :thinplate [_ & [sigma]] (gen-smile (ThinPlateSplineKernel. sigma)))

(def ^{:doc "List of kernel names."} kernel-names (sort (keys (methods kernel))))
