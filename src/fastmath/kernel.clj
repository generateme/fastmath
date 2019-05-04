(ns fastmath.kernel
  (:require [fastmath.core :as m]
            [fastmath.distance :as d]
            [fastmath.vector :as v])
  (:import [smile.math.rbf RadialBasisFunction]
           [clojure.lang IFn]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;;
;; RBF kernels

;; https://www.math.unipd.it/~demarchi/RBF/LectureNotes.pdf
;; http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol4/Choosing_Basis_Functions_and_Shape_Parameters.pdf

(defmacro ^:private make-scaled-x
  [scale formula]
  (let [s (symbol "xs")]
    `(fn [^double x#] (let [~s (/ (m/abs x#) ~scale)] ~formula))))

(defmacro ^:private make-scaled-x+
  [scale formula]
  (let [s (symbol "xs")]
    `(fn [^double x#] (let [~s (/ (m/abs x#) ~scale)]
                       (if (< ~s 1.0) ~formula 0.0)))))

(defmacro ^:private emit-simple-rbf
  [nm form]
  `(defmethod rbf ~nm
     ([a#] (rbf ~nm 1.0))
     ([a# ^double scale#] (make-scaled-x scale# ~form))))

(defmacro ^:private emit-simple-rbf+
  [nm form]
  `(defmethod rbf ~nm
     ([a#] (rbf ~nm 1.0))
     ([a# ^double scale#] (make-scaled-x+ scale# ~form))))

(defmacro ^:private emit-beta-rbf
  [nm form]
  (let [beta (symbol "beta")]
    `(defmethod rbf ~nm
       ([a#] (rbf ~nm 1.0))
       ([a# ^double scale#] (rbf ~nm 1.0 scale#))
       ([a# ^double ~beta ^double scale#] (make-scaled-x scale# ~form)))))

(defmulti rbf (fn [k & _] k))

(emit-simple-rbf :linear xs)
(emit-simple-rbf :gaussian (m/exp (- (* xs xs))))

(defmethod rbf :multiquadratic
  ([_] (rbf :multiquadratic 1.0))
  ([_ scale] (rbf :multiquadratic 0.5 scale))
  ([_ ^double beta ^double scale] (let [s2 (* scale scale)]
                                    (if (== beta 0.5)
                                      (fn [^double x] (m/sqrt (+ s2 (* x x))))
                                      (fn [^double x] (m/pow (+ s2 (* x x)) beta))))))

(defmethod rbf :inverse-multiquadratic
  ([_] (rbf :inverse-multiquadratic 1.0))
  ([_ scale] (rbf :inverse-multiquadratic 0.5 scale))
  ([_ ^double beta ^double scale] (let [s2 (* scale scale)
                                        beta- (- beta)]
                                    (if (== beta 0.5)
                                      (fn [^double x] (/ (m/sqrt (+ s2 (* x x)))))
                                      (fn [^double x] (m/pow (+ s2 (* x x)) beta-))))))

(defmethod rbf :truncated-power
  ([_] (rbf :truncated-power 1.0))
  ([_ scale] (rbf :truncated-power 1.0 scale))
  ([_ ^double k ^double scale] (cond
                                 (== k 1.0) (make-scaled-x+ scale (- 1.0 xs))
                                 (== k 2.0) (make-scaled-x+ scale (m/sq (- 1.0 xs)))
                                 (== k 3.0) (make-scaled-x+ scale (m/pow3 (- 1.0 xs)))
                                 :else (make-scaled-x+ scale (m/pow (- 1.0 xs) k)))))


;; https://www.math.unipd.it/~demarchi/RBF/LectureNotes.pdf
;; page 33

(emit-simple-rbf :gaussians-laguerre-11 (let [x2 (* xs xs)]
                                          (* (- 1.5 x2) (m/exp (- x2)))))


(emit-simple-rbf :gaussians-laguerre-12 (let [x2 (* xs xs)]
                                          (* (+ 1.875 (* -2.5 x2) (* 0.5 x2 x2)) (m/exp (- x2)))))

(emit-simple-rbf :gaussians-laguerre-21 (let [x2 (* xs xs)]
                                          (* (- 2.0 x2) (m/exp (- x2)))))

(emit-simple-rbf :gaussians-laguerre-22 (let [x2 (* xs xs)]
                                          (* (+ 3.0 (* -3.0 x2) (* 0.5 x2 x2)) (m/exp (- x2)))))

;; page 34

(emit-simple-rbf :poisson-2 (m/bessel-j 0 xs))
(emit-simple-rbf :poisson-3 (* 0.7978845608028654 (m/sinc xs)))
(emit-simple-rbf :poisson-4 (/ (m/bessel-j 1 xs) xs))

;; page 35
;; also http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol4/Choosing_Basis_Functions_and_Shape_Parameters.pdf
;; page 193

(emit-simple-rbf :mattern-c0 (m/exp (- xs)))
(emit-simple-rbf :mattern-c2 (* (inc xs) (m/exp (- xs))))
(emit-simple-rbf :mattern-c4 (* (+ 3.0 (* 3.0 xs) (* xs xs)) (m/exp (- xs))))

;; page 37

(emit-simple-rbf :whittaker-02 (+ 1.0 (- xs) (* xs (m/exp (/ -1.0 xs)))))
(emit-simple-rbf :whittaker-12 (+ 1.0 (* -2.0 xs) (* (inc (* 2.0 xs)) (m/exp (/ -1.0 xs)))))
(emit-simple-rbf :whittaker-03 (let [x2 (* 2.0 xs xs)]
                                 (+ 1.0 (* -2.0 xs) x2 (* (- x2) (m/exp (/ -1.0 xs))))))
(emit-simple-rbf :whittaker-13 (let [x2 (* 6.0 xs xs)]
                                 (+ 1.0 (* -4.0 xs) x2 (* (- (+ xs xs x2)) (m/exp (/ -1.0 xs))))))

;; page 43

(emit-beta-rbf :radial-powers (m/pow xs beta))

(defmethod rbf :thin-plate
  ([_] (rbf :thin-plate 1.0))
  ([_ scale] (rbf :thin-plate 1.0 scale))
  ([_ ^double beta ^double scale] (if (== beta 1.0)
                                    (make-scaled-x scale (* xs xs (m/log xs)))
                                    (make-scaled-x scale (* (m/pow (* xs xs) beta) (m/log xs))))))

;; https://www.researchgate.net/profile/Zongmin_Wu/publication/246909840_Multivariate_compactly_supported_positive_definite_radial_functions/links/542247e30cf26120b7a0209b/Multivariate-compactly-supported-positive-definite-radial-functions.pdf
;; page 9

(emit-simple-rbf+ :wu-00 (- 1.0 xs))
(emit-simple-rbf+ :wu-10 (let [r- (- 1.0 xs)]
                           (* r- r- r- (inc (* xs (+ 3.0 xs))))))
(emit-simple-rbf+ :wu-11 (let [r- (- 1.0 xs)]
                           (* r- r- (+ 2.0 xs))))
(emit-simple-rbf+ :wu-20 (let [r- (- 1.0 xs)
                               r2- (* r- r-)
                               xs2 (* xs xs)]
                           (* r- r2- r2- (inc (+ (* 5 xs)
                                                 (* 9 xs2)
                                                 (* 5 xs xs2)
                                                 (* xs2 xs2))))))
(emit-simple-rbf+ :wu-21 (let [r- (- 1.0 xs)
                               r2- (* r- r-)
                               xs2 (* xs xs)]
                           (* r2- r2- (+ 4.0
                                         (* 16.0 xs)
                                         (* 12.0 xs2)
                                         (* 3.0 xs xs2)))))

(emit-simple-rbf+ :wu-22 (let [r- (- 1.0 xs)
                               r2- (* r- r-)]
                           (* r- r2- (+ 8.0
                                        (* 9.0 xs)
                                        (* 3.0 xs xs)))))

(emit-simple-rbf+ :wu-30 (let [r- (- 1.0 xs)
                               r3- (* r- r- r-)
                               xs2 (* xs xs)
                               xs3 (* xs2 xs)]
                           (* r- r3- r3- (+ 5.0
                                            (* 35.0 xs)
                                            (* 101.0 xs2)
                                            (* 147.0 xs xs2)
                                            (* 101.0 xs2 xs2)
                                            (* 35.0 xs3 xs2)
                                            (* 5.0 xs3 xs3)))))

(emit-simple-rbf+ :wu-31 (let [r- (- 1.0 xs)
                               r3- (* r- r- r-)
                               xs2 (* xs xs)
                               xs4 (* xs2 xs2)]
                           (* r3- r3- (+ 6.0
                                         (* 36.0 xs)
                                         (* 82.0 xs2)
                                         (* 72.0 xs xs2)
                                         (* 30.0 xs4)
                                         (* 5.0 xs xs4)))))

(emit-simple-rbf+ :wu-32 (let [r- (- 1.0 xs)
                               r2- (* r- r-)
                               xs2 (* xs xs)]
                           (* r- r2- r2- (+ 8.0
                                            (* 40.0 xs)
                                            (* 48.0 xs2)
                                            (* 25.0 xs xs2)
                                            (* 5.0 xs2 xs2)))))

(emit-simple-rbf+ :wu-33 (let [r- (- 1.0 xs)
                               r2- (* r- r-)
                               xs2 (* xs xs)]
                           (* r2- r2- (+ 16.0
                                         (* 29.0 xs)
                                         (* 20.0 xs2)
                                         (* 5.0 xs xs2)))))

;; wendland

(emit-simple-rbf+ :wendland-10 (- 1.0 xs))
(emit-simple-rbf+ :wendland-21 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)]
                                 (* r- r2- (inc (* 3.0 xs)))))
(emit-simple-rbf+ :wendland-32 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)]
                                 (* r- r2- r2- (inc (+ (* 8.0 xs xs)
                                                       (* 5.0 xs))))))
(emit-simple-rbf+ :wendland-20 (m/sq (- 1.0 xs)))
(emit-simple-rbf+ :wendland-31 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)]
                                 (* r2- r2- (inc (* 4.0 xs)))))
(emit-simple-rbf+ :wendland-42 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)]
                                 (* r2- r2- r2- (+ (* 35.0 xs xs)
                                                   (* 18.0 xs)
                                                   3.0))))
(emit-simple-rbf+ :wendland-53 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)
                                     r4- (* r2- r2-)
                                     xs2 (* xs xs)]
                                 (* r4- r4- (inc (+ (* 32.0 xs2 xs)
                                                    (* 25.0 xs2)
                                                    (* 8.0 xs))))))

(emit-simple-rbf+ :wendland-30 (m/pow3 (- 1.0 xs)))
(emit-simple-rbf+ :wendland-41 (let [r- (- 1.0 xs)
                                     r2- (* r- r-)]
                                 (* r- r2- r2- (inc (* 5.0 xs)))))
(emit-simple-rbf+ :wendland-52 (let [r- (- 1.0 xs)
                                     r3- (* r- r- r-)]
                                 (* r- r3- r3- (inc (+ (* 16.0 xs xs)
                                                       (* 7.0 xs))))))


(defn smile-rbf
  "Create RBF Smile object.

  Used to pass to Smile constructors/functions."
  {:metadoc/categories #{:rbf}}
  [rbf-fn]
  (reify
    RadialBasisFunction (^double f [_ ^double x] (rbf-fn x))
    IFn (invoke ^double [_ x] (rbf-fn x))))


;;;;;;;;;;;;;;;;;;;;;
;; Mattern kernels / Covariance kernels

(defn rbf->mercer
  "Treat RBF kernel as Mercer kernel for given distance."
  ([rbf-kernel] (rbf->mercer rbf-kernel d/euclidean))
  ([rbf-kernel distance]
   (fn [x y] (rbf-kernel (distance x y)))))

;; http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/

(defmulti mercer (fn [k & _] k))

(defmethod mercer :linear
  ([_] (fn [x y] (v/dot x y)))
  ([_ ^double alpha] (fn [x y] (* alpha ^double (v/dot x y)))))

(defmethod mercer :polynomial
  ([_] (fn [x y] (m/sq (v/dot x y))))
  ([_ ^double alpha ^double c ^double d] (fn [x y] (m/pow (+ c (* alpha ^double (v/dot x y))) d))))

(defmethod mercer :gaussian
  ([_] (mercer :gaussian 1.0))
  ([_ sigma] (mercer :gaussian sigma d/euclidean))
  ([_ ^double sigma distance]
   (let [s2r  (/ (* sigma sigma))]
     (fn [x y] (m/exp (* -0.5 s2r (m/sq (distance x y))))))))

(defmethod mercer :exponential
  ([_] (mercer :exponential 1.0))
  ([_ sigma] (mercer :exponential sigma d/euclidean))
  ([_ ^double sigma distance]
   (let [s2r (/ (* sigma sigma))]
     (fn [x y] (m/exp (* -0.5 s2r ^double (distance x y)))))))

(defmethod mercer :laplacian
  ([_] (mercer :laplacian 1.0))
  ([_ sigma] (mercer :laplacian sigma d/euclidean))
  ([_ ^double sigma distance] (fn [x y] (m/exp (- (/ ^double (distance x y) sigma))))))

(defmethod mercer :anova
  ([_] (mercer :anova 1.0))
  ([_ ^double sigma] (mercer :anova sigma 1.0 1.0))
  ([_ ^double sigma ^double k ^double d]
   (let [powk #(m/pow ^double % k)
         powd #(m/pow ^double % d)]
     (fn [x y]
       (let [xk (v/fmap x powk)
             yk (v/fmap y powk)]
         (-> (v/sub xk yk)
             (v/sq)
             (v/mult (- sigma))
             (v/exp)
             (v/fmap powd)
             (v/sum)))))))


(defmethod mercer :hyperbolic-tangent
  ([_] (mercer :hyperbolic-tangent 1.0))
  ([_ ^double alpha] (mercer :hyperbolic-tangent alpha 0.0))
  ([_ ^double alpha ^double c]
   (fn [x y] (m/tanh (+ c (* alpha ^double (v/dot x y)))))))

(defmethod mercer :rational-quadratic
  ([_] (mercer :rational-quadratic 1.0))
  ([_ ^double c] (mercer :rational-quadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (- 1.0 (/ d (+ d c)))))))

(defmethod mercer :multiquadratic
  ([_] (mercer :multiquadratic 1.0))
  ([_ ^double c] (mercer :multiquadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (m/sqrt (+ d (* c c)))))))

(defmethod mercer :inverse-multiquadratic
  ([_] (mercer :inverse-multiquadratic 1.0))
  ([_ ^double c] (mercer :inverse-multiquadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (/ (m/sqrt (+ d (* c c))))))))

(defmethod mercer :circular
  ([_] (mercer :circular 1.0))
  ([_ ^double sigma] (mercer :circular sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (if (>= d sigma) 0.0
                  (let [ds (/ d sigma)]
                    (* 0.6366197723675814 (- (m/acos ds)
                                             (* ds (m/sqrt (- 1.0 (* ds ds))))))))))))

(defmethod mercer :spherical
  ([_] (mercer :spherical 1.0))
  ([_ ^double sigma] (mercer :spherical sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (if (>= d sigma) 0.0
                  (let [ds (/ d sigma)]
                    (+ (- 1.0 (* 1.5 ds)) (* 0.5 ds ds ds))))))))

(defmethod mercer :wave
  ([_] (mercer :wave 1.0))
  ([_ ^double sigma] (mercer :wave sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (* (/ sigma d) (m/sin (/ d sigma)))))))



