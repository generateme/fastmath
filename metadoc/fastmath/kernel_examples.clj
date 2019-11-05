(ns fastmath.kernel-examples
  (:require [fastmath.kernel :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.distance :as d]
            [fastmath.fields :as fields]
            [fastmath.vector :as v]))

(add-examples approx
  (example-session "Usage"
    (let [ak (approx (kernel :gaussian))]
      (ak [1 2] [3 4]))
    (let [ak (approx (kernel :gaussian) 6)]
      (ak [1 2] [3 4])))
  (example-image "Plot of `:gaussian` kernel pproximated to one decimal place."
    "images/k/approx.jpg"))

(add-examples cpd->pd
  (example (let [k (kernel :circular)
                 pd (cpd->pd k)]
             (pd [0.0 -0.1] [0.2 0.4])))
  (example-image "Plot of `:periodic` kernel converted with `cpd->pd` function."
    "images/k/cpdpd.jpg"))

(add-examples kernel-density-list (example "List of density kernels" kernel-density-list))
(add-examples rbf-list (example "List of RBF kernels" rbf-list))
(add-examples kernels-list (example "List of vector kernels" kernels-list))

(add-examples exp
  (example (let [k (exp (kernel :laplacian))]
             (k [1 2] [3 4])))
  (example "Exp with scale=0.5" (let [k (exp (kernel :laplacian) 0.5)]
                                  (k [1 2] [3 4])))
  (example-image "Plot of exp of `:dirichlet` kernel with scaling=5.0"
    "images/k/exp.jpg"))

(add-examples fields
  (example (let [k (kernel :laplacian)
                 fld (fields/field :horseshoe)
                 kf (fields k fld)]
             (kf (v/vec2 2 3) (v/vec2 1 2)))))

(add-examples mult
  (example (let [k1 (kernel :laplacian)
                 k2 (kernel :periodic)
                 res (mult k1 k2)]
             (res [1 2 3 4] [-1 2 3 5])))
  (example-image "Product of two kernels" "images/k/mult.jpg"))

(add-examples wadd
  (example (let [k1 (kernel :laplacian)
                 k2 (kernel :periodic)
                 res (wadd [k1 k2])]
             (res [1 2 3 4] [-1 2 3 5])))
  (example "Weighted sum" (let [k1 (kernel :laplacian)
                                k2 (kernel :periodic)
                                res (wadd [0.2 0.8] [k1 k2])]
                            (res [1 2 3 4] [-1 2 3 5])))
  (example-image "Weighted sum of two kernels" "images/k/wadd.jpg"))

(add-examples scale
  (example (let [k (kernel :laplacian)
                 sk (scale k 2.0)]
             {:kernel (k [1 2 3 4] [-1 2 3 5])
              :scaled-kernel (sk [1 2 3 4] [-1 2 3 5])})))

(add-examples kernel->rbf
  (example (let [k (kernel :mattern-52)
                 r (kernel->rbf k)]
             (r 0.234))))

(add-examples rbf->kernel
  (example (let [r (rbf :mattern-c4)
                 k (rbf->kernel r)]
             (k [1 2 3 4] [-1 2 3 5]))))

(add-examples rbf
  (example-session "Usage"
    (let [k (rbf :gaussian)] (k 1))
    (let [k (rbf :gaussian 0.5)] (k 1)))
  (example "Thin-plate RBF with `beta=2` and `scale=1`"
    (let [k (rbf :thin-plate 2.0 1.0)] (k 0.5)))
  (example-image "Plot of above RBF" "images/k/thin-plate.jpg"))

(add-examples kernel
  (example-session "Usage"
    (let [k (kernel :gaussian)] (k [1 2 3 4] [-1 2 3 5]))
    (let [k (kernel :gaussian 1.0 d/chebyshev)] (k [1 2 3 4] [-1 2 3 5]))
    (let [k (kernel :thin-plate)] (k [1 2 3 4] [-1 2 3 5]))
    (let [k (kernel :mattern-52 0.5)] (k [1 2 3 4] [-1 2 3 5]))))

(add-examples kernel-density
  (example-session "Usage"
    (let [k (kernel-density :epanechnikov (repeatedly 1000 rand))] (k 0.5))
    (let [k (kernel-density :gaussian (repeatedly 1000 rand) 2)] (k 0.5))))

(add-examples kernel-density-ci
  (example-session "Usage"
    (let [k (kernel-density-ci :epanechnikov (repeatedly 1000 rand))] (k 0.5))
    (let [k (kernel-density-ci :gaussian (repeatedly 1000 rand) 2)] (k 0.5)))
  (example-image "Kernel density with confidence intervals" "images/k/ci.jpg"))

(add-examples smile-rbf
  (example (smile-rbf (rbf :mattern-c2))))

(add-examples smile-mercer
  (example (smile-mercer (kernel :mattern-52))))


(defmacro ^:private add-image-examples
  [f pref post lst]
  `(do
     ~@(for [x lst]
         `(add-examples ~f
            (example-image ~(str "Plot of " (name x)) ~(str "images/k/" pref "_" (name x) post))))))

(add-image-examples kernel "k" ".jpg" [:anova :bessel :cauchy :chi-square-cpd :chi-square-pd :circular :dirichlet :exponential :gaussian :generalized-histogram :generalized-t-student :hellinger :histogram :hyperbolic-secant :hyperbolic-tangent :inverse-multiquadratic :laplacian :linear :log :mattern-12 :mattern-32 :mattern-52 :multiquadratic :pearson :periodic :polynomial :power :rational-quadratic :scalar-functions :spherical :spline :thin-plate :variance-function :wave])

(add-image-examples kernel-density "d" ".jpg" [:cauchy :cosine :default :epanechnikov :gaussian :laplace :logistic :quartic :sigmoid :silverman :smile :triangular :tricube :triweight :uniform :wigner])

(add-image-examples rbf "rbf" ".png" [:gaussian :gaussians-laguerre-11 :gaussians-laguerre-12 :gaussians-laguerre-21 :gaussians-laguerre-22 :inverse-multiquadratic :linear :mattern-c0 :mattern-c2 :mattern-c4 :multiquadratic :poisson-2 :poisson-3 :poisson-4 :radial-powers :thin-plate :truncated-power :wendland-10 :wendland-20 :wendland-21 :wendland-30 :wendland-31 :wendland-32 :wendland-41 :wendland-42 :wendland-52 :wendland-53 :whittaker-02 :whittaker-03 :whittaker-12 :whittaker-13 :wu-00 :wu-10 :wu-11 :wu-20 :wu-21 :wu-22 :wu-30 :wu-31 :wu-32 :wu-33])

