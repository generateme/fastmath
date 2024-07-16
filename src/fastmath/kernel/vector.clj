(ns fastmath.kernel.vector
  "Various vector based kernels"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special :as special]))

;; http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/
;; Marc G. Genton, Classes of Kernels for Machine Learning: A Statistics Perspective
;; http://www.jmlr.org/papers/volume2/genton01a/genton01a.pdf


;; https://www.cnblogs.com/chenying99/p/5226126.html

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn linear
  "Dot product kernel"
  ([] (fn ^double [x y] (v/dot x y)))
  ([_] (fn ^double [x y] (v/dot x y))))

(defn polynomial
  "Polynomial (power of dot product).

  Parameters:

  * `:alpha` - dot product multiplier (default: 1.0)
  * `:c` - shift (default: 0.0)
  * `:p` - exponent (default: 2.0)"
  ([] (fn ^double [x y] (m/sq (v/dot x y))))
  ([{:keys [^double alpha ^double c ^double p]
     :or {alpha 1.0 c 0.0 p 2.0}}]
   (fn ^double [x y] (m/pow (m/+ c (m/* alpha (v/dot x y))) p))))

(defn gaussian
  "Gaussian kernel.

  Parameters:

  * `:sigma` - shape of the kernel (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (gaussian nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv (m/sq (double (distance x y)))))))))

(defn exponential
  "Exponential kernel.

  Parameters:

  * `:sigma` - shape of the kernel (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (exponential nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv ^double (distance x y)))))))

(defn laplacian
  "Laplacian kernel.

  Parameters:

  * `:sigma` - shape of the kernel (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (laplacian nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-sinv (m// -1.0 sigma)]
     (fn ^double [x y] (m/exp (m/* -sinv ^double (distance x y)))))))

(defn anova
  "Anova kernel.

  Parameters:

  * `:sigma` - multiplier (default: 1.0)
  * `:k` and `:d` - exponents (default: 1.0)"
  ([] (anova nil))
  ([{:keys [^double sigma ^double k ^double d]
     :or {sigma 1.0 k 1.0 d 1.0}}]
   (let [powk (fn ^double [^double v] (m/pow v k))
         powd (fn ^double [^double v] (m/pow v d))
         -s (m/- sigma)]
     (fn ^double [x y]
       (let [xk (v/fmap x powk)
             yk (v/fmap y powk)]
         (-> (v/sub xk yk)
             (v/sq)
             (v/mult -s)
             (v/exp)
             (v/fmap powd)
             (v/sum)))))))

(defn hyperbolic-tangent
  "Hyperbolic tangent of the dot product.

  Parameters:

  * `:alpha` - dot product multiplier (default: 1.0)
  * `:c` - shift (default: 0.0)"
  ([] (hyperbolic-tangent nil))
  ([{:keys [^double alpha ^double c]
     :or {alpha 1.0 c 0.0}}]
   (fn ^double [x y] (m/tanh (m/+ c (m/* alpha (v/dot x y)))))))

(defn rational-quadratic
  "Rational quadratic kernel.

  Parameters:

  `:c` - shift (default: 1.0)
  `:distance` - distance function (default: euclidean)"
  ([] (rational-quadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                      (m/- 1.0 (m// dist (m/+ dist c)))))))

(defn multiquadratic
  "Multiquadratic kernel.

  Parameters:

  `:c` - shift (default: 1.0)
  `:distance` - distance function (default: euclidean)"
  ([] (multiquadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (let [c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                        (m/sqrt (m/+ dist c2)))))))

(defn inverse-multiquadratic
  "Inverse multiquadratic kernel.

  Parameters:

  `:c` - shift (default: 1.0)
  `:distance` - distance function (default: euclidean)"
  ([] (inverse-multiquadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (let [c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                        (m// (m/sqrt (m/+ dist c2))))))))

;; http://perso.lcpc.fr/tarel.jean-philippe/publis/jpt-icann05b.pdf

(defn- geometric-internal
  ^double [^long n ^double ds]
  (loop [k (long 3)
         p1 (m/- (m/acos ds) (m/* ds (m/sqrt (m/- 1.0 (m/* ds ds)))))
         p2 (m/- 1.0 ds)]
    (if (m/> k n)
      p1
      (let [k- (m/- k 1.0)
            curr (m/- (m/* (m// k- k) p2)
                      (m/* (m// 1.0 k) ds
                           (m/pow (m/- 1.0 (m/* ds ds)) (m/* 0.5 k-))))]
        (recur (m/inc k) curr p1)))))

(defn geometric
  "Geometric Compactly Supported kernel

  Parameters:

  * `:n` - dimension
  * `:r` - shape (default: 1.0)
  * `:distance` - distance function (default: euclidean)

  Specific kernel names for `:n`:
  
  * 1 - triangular
  * 2 - circular
  * 3 - spherical"
  ([] (geometric nil))
  ([{:keys [^long n ^double r distance]
     :or {r 1.0 distance v/dist}}]
   (let [r2 (m/* 2.0 r)]
     (case (int n)
       1 (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (>= dist r2) 0.0
                                (m/- 1.0 (m// dist r2)))))
       2 (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (>= dist r2) 0.0
                                (let [ds (m// dist r2)]
                                  (m/- (m/acos ds)
                                       (m/* ds (m/sqrt (m/- 1.0 (m/* ds ds)))))))))
       3 (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (>= dist r2) 0.0
                                (let [ds (m// dist r2)]
                                  (m/+ (m/- 1.0 (m/* 1.5 ds)) (m/* 0.5 ds ds ds))))))
       (let [r2 (m/* 2.0 r)
             z (geometric-internal n 0.0)]
         (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (m/>= dist r2) 0.0
                                (m// (geometric-internal n (m// dist r2)) z)))))))))

(defn wave
  "Wave (sinc) kernel.

  Parameters:

  * `:sigma` - scale (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (wave nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (fn ^double [x y] (let [^double dist (distance x y)]
                      (if (m/zero? dist) 1.0
                          (m/* (m// sigma dist) (m/sin (m// dist sigma))))))))

(defn periodic
  "Periodic kernel.

  Parameters:

  * `:sigma` - scale (default: 1.0)
  * `:periodicity` - periodicity (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (periodic nil))
  ([{:keys [^double sigma ^double periodicity distance]
     :or {sigma 1.0 periodicity 1.0 distance v/dist}}]
   (let [p (m// m/PI periodicity)
         s2 (m/* sigma sigma)]
     (fn ^double [x y] (let [^double dist (distance x y)]
                        (m/exp (m// (m/* -2.0 (m/sq (m/sin (m/* p dist)))) s2)))))))

(defn power
  "Power (negative) of distance.

  Parameters:

  * `:p` - exponent (default: 2.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (power nil))
  ([{:keys [^double p distance]
     :or {p 2.0 distance v/dist}}]
   (fn [x y] (m/- (m/pow (distance x y) p)))))

(defn log
  "Logarithmic.

  Parameters:

  * `:p` - exponent (default: 2.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (log nil))
  ([{:keys [^double p distance]}]
   (fn [x y] (m/- (m/log (m/inc (m/pow (distance x y) p)))))))

(defn spline
  "Spline kernel."
  ([] (spline nil))
  ([_] (fn [x y]
         (reduce m/* 1.0 (map (fn [^double xi ^double yi]
                                (let [xiyi (m/* xi yi)
                                      m (min xi yi)
                                      m2 (m/* m m)]
                                  (inc (m/+ xiyi
                                            (m/* xiyi m)
                                            (m/* -0.5 (m/+ xi yi) m2)
                                            (m/* m/THIRD m2 m)))))
                              (if (number? x) [x] x)
                              (if (number? y) [y] y))))))


;; https://www.researchgate.net/publication/2538959_Inverse_B-spline_interpolation

(defn- b-spline-beta
  [^long n ^double x]
  (let [n!inv (m/inv-factorial n)
        n+ (m/inc n)
        nx (m/+ x (m/* 0.5 n+))]
    (loop [k (long 0)
           neg 1.0
           curr 0.0]
      (if (or (m/> k n+)
              (m/<= nx k))
        (m/* n!inv curr)
        (recur (m/inc k)
               (m/* -1.0 neg)
               (m/+ curr (m/* neg (m/combinations n+ k)
                              (m/fpow (m/- nx k) n))))))))

(defn b-spline
  "B-spline kernel with degree `:n` (default: 2.0)."
  ([] (b-spline nil))
  ([{:keys [^long n]
     :or {n 2.0}}]
   (let [nn (m/inc (m/* 2 n))]
     (fn ^double [x y]
       (-> (v/sub x y)
           (v/abs)
           (v/fmap (fn [^double x] (b-spline-beta nn x)))
           (v/prod))))))

(defn bessel
  "Bessel (of the first kind) kernel

  Parameters:

  * `:sigma` - shape (default: 1.0)
  * `:n` - exponent factor (default: 2.0)
  * `:v` - Bessel J order - 1 (default: -1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (bessel nil))
  ([{:keys [^double sigma ^double n ^double v distance]
     :or {sigma 1.0 n 2.0 v -1.0 distance v/dist}}]
   (fn ^double [x y] (let [^double dist (distance x y)
                          v+ (m/inc v)]
                      (/ (special/bessel-J v+ (m/* sigma dist))
                         (m/pow dist (m/- (m/* n v+))))))))

;; R kernlab
(defn bessel2
  "Bessel (of the first kind) kernel, R kernlab implementation.

  Parameters:

  * `:sigma` - shape (default: 1.0)
  * `:degree` - exponent (default: 1.0)
  * `:order` - Bessel J order (default: 0.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (bessel2 nil))
  ([{:keys [^double sigma ^double order ^double degree distance]
     :or {sigma 1.0 order 0.0 degree 1.0 distance v/dist}}]
   (fn ^double [x y]
     (let [lim (m// (m/* (special/gamma (m/inc order)) (m/pow 2.0 order)))
           bkt (m/* sigma ^double (distance x y))]
       (if (m/< bkt 1.0e-5)
         1.0
         (m/pow (m// (m/* (special/bessel-J order bkt)
                          (m/pow bkt (m/- order))) lim) degree))))))

(defn cauchy
  "Cauchy kernel.

  Parameters:

  * `sigma` - scale (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (cauchy nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (fn ^double [x y] (m// (m/inc (m/sq (m// ^double (distance x y) sigma)))))))

(defn chi-square
  "Chi-square kernel."
  ([] (chi-square nil))
  ([_] (fn ^double [x y] (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                               (m// (m/* 2.0 xi yi)
                                                    (m/+ xi yi)))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y))))))

(defn chi-square2
  "Chi-square kernel, second version"
  ([] (chi-square2 nil))
  ([_] (fn ^double [x y]
         (- 1.0 ^double (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                               (m// (m/sq (m/- xi yi))
                                                    (m/* 0.5 (m/+ xi yi))))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y)))))))

(defn histogram
  "Histogram kernel."
  ([] (histogram nil))
  ([_] (fn ^double [x y] (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                               (min xi yi))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y))))))

(defn generalized-histogram
  "Generalized histogram with `:p` exponent (default: 2.0)."
  ([] (generalized-histogram nil))
  ([{:keys [^double p]
     :or {p 2.0}}]
   (fn ^double [x y] (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                           (min (m/pow (m/abs xi) p)
                                                (m/pow (m/abs yi) p)))
                                         (if (number? x) [x] x)
                                         (if (number? y) [y] y))))))

(defn generalized-t-student
  "Generalized t-student.

  Parameters:

  * `:p` - exponent
  * `:distance` - distance function (default: euclidean)"
  ([] (generalized-t-student nil))
  ([{:keys [^double p distance]}]
   (fn ^double [x y] (m// (m/inc (m/pow (distance x y) p))))))

(defn- dirichlet-dim
  ^double [^double nhalf ^double delta]
  (if (m/zero? delta)
    nhalf
    (let [den (m/* 2.0 (m/sin (m/* 0.5 delta)))]
      (if (m/zero? den)
        nhalf
        (m// (m/sin (m/* nhalf delta)) den)))))

(defn dirichlet
  "Dirichlet kernel with `:n` dimensionality (default: 1.0)."
  ([] (dirichlet nil))
  ([{:keys [^double n]
     :or {n 1.0}}]
   (let [nhalf (m/+ 0.5 n)]
     (fn ^double [x y]
       (let [diff (v/sub x y)]
         (if (number? diff)
           (dirichlet-dim nhalf diff)
           (reduce m/* 1.0 (map (partial dirichlet-dim nhalf) diff))))))))

(defn hellinger
  "Hellinger kernel."
  ([] (hellinger nil))
  ([_] (fn ^double [x y] (-> (v/emult x y) v/safe-sqrt v/sum))))

(defn pearson
  "Pearson VII kernel

  Parameters:

  * `:sigma` - scale (default: 1.0)
  * `:omega` - exponent (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (pearson nil))
  ([{:keys [^double sigma ^double omega distance]
     :or {sigma 1.0 omega 1.0 distance v/dist}}]
   (let [c (m// (m/* 4.0 (m/sqrt (m/dec (m/pow 2.0 (m// omega)))))
                (m/sq sigma))]
     (fn ^double [x y]
       (m// (m/pow (m/inc (m/* c (m/sq (double (distance x y))))) omega))))))

(defn hyperbolic-secant
  "Hyperbolic secant kernel.

  Parameters:

  * `:a` scaling factor (default: 1.0)
  * `:distance` - distance function (default: euclidean)"
  ([] (hyperbolic-secant nil))
  ([{:keys [^double a distance]
     :or {a 1.0 distance v/dist}}]
   (fn ^double [x y]
     (m/sech (m/* a ^double (distance x y))))))

(defn matern
  "Matern kernel.

  Parameters:

  * `:order` - order of the kernel, should be odd (default: 1).
  * `:theta` - shape (default: 1.0)
  * `:distance` - distance function (default: euclidean)

  Order of the Bessel K function is a half of `:order` parameter. For example to get Matern 5/2 kernel, call `(matern 5)`."
  ([] (matern nil))
  ([{:keys [^long order ^double theta distance]
     :or {order 1 theta 1.0 distance v/dist}}]
   (let [mu (m// order 2.0)
         gf (m// (m/* (m/pow 2.0 (m/- 1.0 mu))) (special/gamma mu))
         s (m// (m/* 2.0 (m/sqrt mu)) theta)]
     (fn ^double [x y]
       (let [v (m/* s ^double (distance x y))]
         (if (m/< v 1.0e-16)
           1.0
           (m/* gf (m/pow v mu) (special/bessel-K-half-odd order v))))))))

(defn rbf->kernel
  "Convert RBF kernel as vector kernel using a `distance` function (default: euclidean)."
  ([rbf-kernel] (rbf->kernel rbf-kernel v/dist))
  ([rbf-kernel distance]
   (fn [x y] (rbf-kernel (distance x y)))))
