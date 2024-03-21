(ns fastmath.kernel.vector
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]))

;; https://www.cnblogs.com/chenying99/p/5226126.html

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn linear [] (fn ^double [x y] (v/dot x y)))

(defn polynomial
  ([] (fn ^double [x y] (m/sq (v/dot x y))))
  ([^double alpha ^double c ^double p]
   (fn ^double [x y] (m/pow (m/+ c (m/* alpha (v/dot x y))) p))))

(defn gaussian
  ([] (gaussian 1.0))
  ([^double sigma] (gaussian sigma nil))
  ([^double sigma distance]
   (let [d (or distance v/dist)
         -s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv (m/sq (d x y))))))))

(defn exponential
  ([] (exponential 1.0))
  ([^double sigma] (exponential sigma nil))
  ([^double sigma distance]
   (let [d (or distance v/dist)
         -s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv ^double (d x y)))))))

(defn laplacian
  ([] (laplacian 1.0))
  ([^double sigma] (laplacian sigma nil))
  ([^double sigma distance]
   (let [d (or distance v/dist)
         -sinv (m// -1.0 sigma)]
     (fn ^double [x y] (m/exp (m/* -sinv ^double (d x y)))))))

(defn anova
  ([] (anova 1.0))
  ([^double sigma] (anova sigma 1.0 1.0))
  ([^double sigma ^double k ^double d]
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
  ([] (hyperbolic-tangent 1.0))
  ([^double alpha] (hyperbolic-tangent alpha 0.0))
  ([^double alpha ^double c]
   (fn ^double [x y] (m/tanh (m/+ c (m/* alpha (v/dot x y)))))))

(defn rational-quadratic
  ([] (rational-quadratic 1.0))
  ([^double c] (rational-quadratic c nil))
  ([^double c distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y] (let [dist (m/sq (d x y))]
                        (m/- 1.0 (m// dist (m/+ dist c))))))))

(defn multiquadratic
  ([] (multiquadratic 1.0))
  ([^double c] (multiquadratic c nil))
  ([^double c distance]
   (let [d (or distance v/dist)
         c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (d x y))]
                        (m/sqrt (m/+ dist c2)))))))

(defn inverse-multiquadratic
  ([] (inverse-multiquadratic 1.0))
  ([^double c] (inverse-multiquadratic c nil))
  ([^double c distance]
   (let [d (or distance v/dist)
         c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (d x y))]
                        (m// (m/sqrt (m/+ dist c2))))))))

;; http://perso.lcpc.fr/tarel.jean-philippe/publis/jpt-icann05b.pdf

(defn triangular
  ([] (triangular 1.0))
  ([^double r] (triangular r nil))
  ([^double r distance]
   (let [d (or distance v/dist)
         r2 (m/* 2.0 r)]
     (fn ^double [x y] (let [^double dist (d x y)]
                        (if (>= dist r2) 0.0
                            (m/- 1.0 (m// dist r2))))))))

(defn circular
  ([] (circular 1.0))
  ([^double r] (circular r nil))
  ([^double r distance]
   (let [d (or distance v/dist)
         r2 (m/* 2.0 r)]
     (fn ^double [x y] (let [^double dist (d x y)]
                        (if (>= dist r2) 0.0
                            (let [ds (m// dist r2)]
                              (m/- (m/acos ds)
                                   (m/* ds (m/sqrt (m/- 1.0 (m/* ds ds))))))))))))

(defn spherical
  ([] (spherical 1.0))
  ([^double r] (spherical r nil))
  ([^double r distance]
   (let [d (or distance v/dist)
         r2 (m/* 2.0 r)]
     (fn ^double [x y] (let [^double dist (d x y)]
                        (if (>= dist r2) 0.0
                            (let [ds (m// dist r2)]
                              (m/+ (m/- 1.0 (m/* 1.5 ds)) (m/* 0.5 ds ds ds)))))))))

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
  "Geometric Compactly Supported kernel"
  ([^long n] (geometric n 1.0))
  ([^long n ^double r] (geometric n r nil))
  ([^long n ^double r distance]
   (case (int n)
     1 (triangular r distance)
     2 (circular r distance)
     3 (spherical r distance)
     (let [d (or distance v/dist)
           r2 (m/* 2.0 r)
           z (geometric-internal n 0.0)]
       (fn ^double [x y] (let [^double dist (d x y)]
                          (if (m/>= dist r2) 0.0
                              (m// (geometric-internal n (m// dist r2)) z))))))))

(defn wave
  ([] (wave 1.0))
  ([^double sigma] (wave sigma nil))
  ([^double sigma distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y] (let [^double dist (d x y)]
                        (if (m/zero? dist) 1.0
                            (m/* (m// sigma dist) (m/sin (m// dist sigma)))))))))

(defn periodic
  ([] (periodic 1.0))
  ([^double sigma] (periodic sigma 1.0))
  ([^double sigma periodicity] (periodic sigma periodicity nil))
  ([^double sigma ^double periodicity distance]
   (let [d (or distance v/dist)
         p (m// m/PI periodicity)
         s2 (m/* sigma sigma)]
     (fn ^double [x y] (let [^double dist (d x y)]
                        (m/exp (m// (m/* -2.0 (m/sq (m/sin (m/* p dist)))) s2)))))))

(defn power
  ([] (power 2.0))
  ([^double p] (power p nil))
  ([^double p distance]
   (let [d (or distance v/dist)]
     (fn [x y] (m/- (m/pow (d x y) p))))))

(defn log
  ([] (log 2.0))
  ([^double p] (log p nil))
  ([^double p distance]
   (let [d (or distance v/dist)]
     (fn [x y] (m/- (m/log (m/inc (m/pow (d x y) p))))))))

(defn spline
  [] (fn [x y]
       (reduce m/fast* 1.0 (map (fn [^double xi ^double yi]
                                  (let [xiyi (m/* xi yi)
                                        m (min xi yi)
                                        m2 (m/* m m)]
                                    (inc (m/+ xiyi
                                              (m/* xiyi m)
                                              (m/* -0.5 (m/+ xi yi) m2)
                                              (m/* m/THIRD m2 m)))))
                                (if (number? x) [x] x)
                                (if (number? y) [y] y)))))


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
  ([] (b-spline 2))
  ([^long n]
   (let [nn (m/inc (m/* 2 n))]
     (fn ^double [x y]
       (-> (v/sub x y)
           (v/abs)
           (v/fmap (fn [^double x] (b-spline-beta nn x)))
           (v/prod))))))

(defn bessel
  ([] (bessel 1.0))
  ([^double sigma] (bessel sigma 2.0))
  ([^double sigma ^double n] (bessel sigma n -1.0))
  ([^double sigma ^double n ^double v] (bessel sigma n v nil))
  ([^double sigma ^double n ^double v distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y] (let [^double dist (d x y)
                            v+ (m/inc v)]
                        (/ (m/bessel-j v+ (m/* sigma dist))
                           (m/pow dist (m/- (m/* n v+)))))))))

;; R kernlab
(defn bessel2
  ([] (bessel2 1.0))
  ([^double sigma] (bessel2 sigma 0.0))
  ([^double sigma ^double order] (bessel2 sigma order 1.0))
  ([^double sigma ^double order ^double degree] (bessel2 sigma order degree nil))
  ([^double sigma ^double order ^double degree distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y]
       (let [lim (m// (m/* (m/gamma (m/inc order)) (m/pow 2.0 order)))
             bkt (m/* sigma ^double (d x y))]
         (if (m/< bkt 1.0e-5)
           1.0
           (m/pow (m// (m/* (m/bessel-j order bkt)
                            (m/pow bkt (m/- order))) lim) degree)))))))

(defn cauchy
  ([] (cauchy 1.0))
  ([^double sigma] (cauchy sigma nil))
  ([^double sigma distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y] (m// (m/inc (m/sq (m// ^double (d x y) sigma))))))))

(defn chi-square []
  (fn ^double [x y] (reduce m/fast+ 0.0 (map (fn [^double xi ^double yi]
                                              (m// (m/* 2.0 xi yi)
                                                   (m/+ xi yi)))
                                            (if (number? x) [x] x)
                                            (if (number? y) [y] y)))))

(defn chi-square2 []
  (fn ^double [x y]
    (- 1.0 ^double (reduce m/fast+ 0.0 (map (fn [^double xi ^double yi]
                                              (m// (m/sq (m/- xi yi))
                                                   (m/* 0.5 (m/+ xi yi))))
                                            (if (number? x) [x] x)
                                            (if (number? y) [y] y))))))

(defn histogram []
  (fn ^double [x y] (reduce m/fast+ 0.0 (map (fn [^double xi ^double yi]
                                              (min xi yi))
                                            (if (number? x) [x] x)
                                            (if (number? y) [y] y)))))

(defn generalized-histogram
  ([] (generalized-histogram 2.0))
  ([^double p]
   (fn ^double [x y] (reduce m/fast+ 0.0 (map (fn [^double xi ^double yi]
                                               (min (m/pow (m/abs xi) p)
                                                    (m/pow (m/abs yi) p)))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y))))))

(defn generalized-t-student
  ([] (generalized-t-student 1.0))
  ([^double p] (generalized-t-student p nil))
  ([^double p distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y] (m// (m/inc (m/pow (d x y) p)))))))

(defn- dirichlet-dim
  ^double [^double nhalf ^double delta]
  (if (m/zero? delta)
    nhalf
    (let [den (m/* 2.0 (m/sin (m/* 0.5 delta)))]
      (if (m/zero? den)
        nhalf
        (m// (m/sin (m/* nhalf delta)) den)))))

(defn dirichlet
  ([] (dirichlet 1.0))
  ([^double n]
   (let [nhalf (m/+ 0.5 n)]
     (fn ^double [x y]
       (let [diff (v/sub x y)]
         (if (number? diff)
           (dirichlet-dim nhalf diff)
           (reduce m/fast* 1.0 (map (partial dirichlet-dim nhalf) diff))))))))

(defn hellinger []
  (fn ^double [x y] (-> (v/emult x y) v/safe-sqrt v/sum)))

(defn pearson
  ([] (pearson 1.0 1.0))
  ([^double sigma ^double omega] (pearson sigma omega nil))
  ([^double sigma ^double omega distance]
   (let [d (or distance v/dist)
         c (m// (m/* 4.0 (m/sqrt (m/dec (m/pow 2.0 (m// omega)))))
                (m/sq sigma))]
     (fn ^double [x y]
       (m// (m/pow (m/inc (m/* c (m/sq (d x y)))) omega))))))

(defn hyperbolic-secant
  ([] (hyperbolic-secant 1.0))
  ([^double a] (hyperbolic-secant a nil))
  ([^double a distance]
   (let [d (or distance v/dist)]
     (fn ^double [x y]
       (m/sech (m/* a ^double (d x y)))))))

(defn matern
  ([] (matern 1))
  ([^long order] (matern order nil))
  ([^long order distance] (matern order 1.0 distance))
  ([^long order ^double theta distance]
   (let [d (or distance v/dist)
         mu (/ order 2.0)
         gf (m// (m/* (m/pow 2.0 (m/- 1.0 mu))) (m/gamma mu))
         s (m// (m/* 2.0 (m/sqrt mu)) theta)]
     (fn ^double [x y]
       (let [v (m/* s ^double (d x y))]
         (if (m/< v 1.0e-16)
           1.0
           (m/* gf (m/pow v mu) (m/bessel-k-half order v))))))))

