(ns fastmath.solver
  (:require [fastmath.vector :as v]
            [fastmath.core :as m])
  (:import [fastmath.vector Vec2 Vec3]
           [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.solvers UnivariateSolver BrentSolver BisectionSolver IllinoisSolver MullerSolver MullerSolver2 PegasusSolver RegulaFalsiSolver RiddersSolver SecantSolver]))

(defn- wrap-univariate-function ^UnivariateFunction [f] (reify UnivariateFunction (value [_ x] (f x))))

(defmacro ^:private solver-constructor
  [clazz]
  `(cond
     ~'_relative-accuracy (new ~clazz (or ~'_absolute-accuracy 1.0e-8) ~'_relative-accuracy)
     ~'_absolute-accuracy (new ~clazz (or ~'_absolute-accuracy 1.0e-8))
     :else (new ~clazz)))

(defn- get-solver
  ^UnivariateSolver [solver-type _absolute-accuracy _relative-accuracy]
  (case solver-type
    :brent (solver-constructor BrentSolver)
    :bisection (solver-constructor BisectionSolver)
    :illinois (solver-constructor IllinoisSolver)
    :muller (solver-constructor MullerSolver)
    :muller2 (solver-constructor MullerSolver2)
    :pegasus (solver-constructor PegasusSolver)
    :regula-falsi (solver-constructor RegulaFalsiSolver)
    :ridders (solver-constructor RiddersSolver)
    :secant (solver-constructor SecantSolver)
    (solver-constructor BrentSolver)))

(defn find-root
  "Find zero (root) of a function `f` in given range [`lower-bound`, `upper-bound`].

  Optional parameters:

  * `:absolute-accuracy` - default 1.0e-8
  * `:relative-accuracy`
  * `:max-iters` - maximum iterations (default: 100)
  * `:initial-value` - algorithm starting value
  * `:solver` - one of: `:brent` (default), `:bisection`, `:illinois`, `:muller`, `:muller2`, `:pegasus`, `:regula-falsi`, `:ridders` and `:secant`."
  (^double [f lower-bound upper-bound] (find-root f lower-bound upper-bound {}))
  (^double [f lower-bound upper-bound {:keys [absolute-accuracy relative-accuracy
                                              ^int max-iters initial-value solver]
                                       :or {max-iters 100 solver :brent}}]
   (let [^UnivariateSolver solver (get-solver solver absolute-accuracy relative-accuracy)]
     (if initial-value
       (.solve solver max-iters (wrap-univariate-function f) (double lower-bound) (double upper-bound) (double initial-value))
       (.solve solver max-iters (wrap-univariate-function f) (double lower-bound) (double upper-bound))))))

(defn quadratic
  "Real solutions of quadratic formula ax^2+bx+c=0. Returns `NaNs` where no real solutions are found."
  ^Vec2 [^double a ^double b ^double c]
  (if (zero? a)
    (if-not (zero? b)
      (let [r (m/- (m// c b))]
        (Vec2. r ##NaN))
      (Vec2. ##NaN ##NaN))
    (let [discrim (m/difference-of-products b b (m/* 4.0 a) c)]
      (when-not (neg? discrim)
        (let [discrim-root (m/sqrt discrim)
              q (m/* -0.5 (m/+ b (m/copy-sign discrim-root b)))
              r0 (m// q a)
              r1 (m// c q)]
          (if (m/< r0 r1)
            (Vec2. r0 r1)
            (Vec2. r1 r0)))))))

(defn cubic
  "Real solution of cubic formula ax^3+bx^2+cx+d=0. Returns `NaNs` where no real solutions are found."
  ^Vec3 [^double a ^double b ^double c ^double d]
  (if (zero? a)
    (v/vec3 (quadratic b c d) ##NaN)
    (let [a2 (* a a)
          b2 (* b b)
          p (/ (- (* 3.0 a c) b2) (* 3.0 a2))
          q (/ (+ (- (* 2.0 b2 b) (* 9.0 a b c)) (* 27.0 a2 d)) (* 27.0 a2 a))
          l (* 4.0 p p p)
          r (* -27.0 q q)]
      (v/shift (cond
                 (zero? p) (let [root (m/cbrt (- q))]
                             (Vec3. root root root))
                 (< l r) (let [sqp3 (m/sqrt (/ p -3.0))
                               theta (m/acos (/ (* 3.0 q) (* 2.0 p sqp3)))]
                           (-> (Vec3. (m/cos (/ theta 3.0))
                                      (m/cos (/ (+ theta m/TWO_PI) 3.0))
                                      (m/cos (/ (- theta m/TWO_PI) 3.0)))
                               (v/mult (* 2.0 sqp3))))
                 (== l r) (let [qp (/ (* 3.0 q) p)
                                qp2 (* -0.5 qp)]
                            (Vec3. qp qp2 qp2))
                 :else (let [hq (/ q -2.0)
                             v (m/sqrt (/ (- l r) 108.0))
                             root (+ (m/cbrt (+ hq v))
                                     (m/cbrt (- hq v)))]
                         (Vec3. root ##NaN ##NaN)))
               (/ b (* -3.0 a))))))

