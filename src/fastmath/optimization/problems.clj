(ns fastmath.optimization.problems
  "A collection of test functions for optimization with known global minima.

  Univariate problems follow the collection of one dimensional global optimization problems by Gavana (https://infinity77.net/global_optimization/test_functions_1d.html). For every problem there are three functions: `problemNN` takes a number and returns the value, `problemNN-bounds` returns the search domain as a vector with one `[lo hi]` pair, and `dproblemNN` is the derivative, which takes a sequence with one number (the form used by the optimizers) and returns a vector with the derivative. The numbers follow the collection, so some of them are missing.

  Multivariate problems follow the virtual library of simulation experiments by Surjanovic and Bingham (https://www.sfu.ca/~ssurjano/optimization.html): [[->ackley]], [[bukin-no-6]], [[cross-in-tray]], [[drop-wave]], [[egg-holder]], [[rosenbrock]], [[himmelblau]], [[beale]] and [[sphere]]. They take the point as one sequence of numbers and return a double. Every problem has a `-bounds` function with its usual search domain (a function of the number of dimensions for the problems with any number of dimensions). The gradients of [[rosenbrock]], [[himmelblau]], [[beale]] and [[sphere]] are given by [[rosenbrock-gradient]], [[himmelblau-gradient]], [[beale-gradient]] and [[sphere-gradient]].

  The docstrings give the global minimum of every function."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; univariate
;; https://infinity77.net/global_optimization/test_functions_1d.html

(defn problem02-bounds
  "Returns the search domain of [[problem02]], the interval `[2.7 7.5]`, as a vector with one `[lo hi]` pair.

  See also [[problem02]], [[dproblem02]]."
  [] [[2.7 7.5]])

(defn problem02
  "Univariate test function number 02 of the collection of one dimensional global optimization problems: `sin(x) + sin(10x/3)`.

  Domain: `[2.7, 7.5]`, see [[problem02-bounds]]. Global minimum: `f(x*) = -1.899599 at x* = 5.145735`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem02]] (the derivative)."
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* 3.333333333333333 x))))

(defn dproblem02
  "Derivative of the univariate test function [[problem02]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem02]]."
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* 3.333333333333333 (m/cos (m/* 3.333333333333333 x))))])

;;

(defn problem03-bounds
  "Returns the search domain of [[problem03]], the interval `[-10.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem03]], [[dproblem03]]."
  [] [[-10.0 10.0]])

(defn problem03
  "Univariate test function number 03 of the collection of one dimensional global optimization problems: `-sum of k sin((k+1)x + k) for k = 1..5`.

  Domain: `[-10, 10]`, see [[problem03-bounds]]. Global minimum: `f(x*) = -12.031249 at x* = -6.774576, -0.491391 and 5.791794` (three minima of the same depth).

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem03]] (the derivative)."
  ^double [^double x]
  (m/- (m/+ (m/sin (m/+ (m/* 2.0 x) 1.0))
            (m/* 2.0 (m/sin (m/+ (m/* 3.0 x) 2.0)))
            (m/* 3.0 (m/sin (m/+ (m/* 4.0 x) 3.0)))
            (m/* 4.0 (m/sin (m/+ (m/* 5.0 x) 4.0)))
            (m/* 5.0 (m/sin (m/+ (m/* 6.0 x) 5.0))))))

(defn dproblem03
  "Derivative of the univariate test function [[problem03]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem03]]."
  [[^double x]]
  [(m/- (m/+ (m/* 2.0 (m/cos (m/+ (m/* 2.0 x) 1.0)))
             (m/* 6.0 (m/cos (m/+ (m/* 3.0 x) 2.0)))
             (m/* 12.0 (m/cos (m/+ (m/* 4.0 x) 3.0)))
             (m/* 20.0 (m/cos (m/+ (m/* 5.0 x) 4.0)))
             (m/* 30.0 (m/cos (m/+ (m/* 6.0 x) 5.0)))))])

;;

(defn problem04-bounds
  "Returns the search domain of [[problem04]], the interval `[1.9 3.9]`, as a vector with one `[lo hi]` pair.

  See also [[problem04]], [[dproblem04]]."
  [] [[1.9 3.9]])

(defn problem04
  "Univariate test function number 04 of the collection of one dimensional global optimization problems: `-(16x^2 - 24x + 5) exp(-x)`.

  Domain: `[1.9, 3.9]`, see [[problem04-bounds]]. Global minimum: `f(x*) = -3.850450 at x* = 2.868034`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem04]] (the derivative)."
  ^double [^double x]
  (m/- (m/* (m/+ (m/* 16.0 x x)
                 (m/* -24.0 x) 5.0)
            (m/exp (m/- x)))))

(defn dproblem04
  "Derivative of the univariate test function [[problem04]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem04]]."
  [[^double x]]
  [(m/* (m/+ (m/* 16.0 x x)
             (m/* -56.0 x) 29.0)
        (m/exp (m/- x)))])

;;

(defn problem05-bounds
  "Returns the search domain of [[problem05]], the interval `[0.0 1.2]`, as a vector with one `[lo hi]` pair.

  See also [[problem05]], [[dproblem05]]."
  [] [[0.0 1.2]])

(defn problem05
  "Univariate test function number 05 of the collection of one dimensional global optimization problems: `-(1.4 - 3x) sin(18x)`.

  Domain: `[0, 1.2]`, see [[problem05-bounds]]. Global minimum: `f(x*) = -1.489073 at x* = 0.966090`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem05]] (the derivative)."
  ^double [^double x]
  (m/- (m/* (m/- 1.4 (m/* 3.0 x))
            (m/sin (m/* 18.0 x)))))

(defn dproblem05
  "Derivative of the univariate test function [[problem05]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem05]]."
  [[^double x]]
  [(let [x18 (m/* 18.0 x)]
     (m/- (m/* 3.0 (m/sin x18))
          (m/* 18.0 (m/- 1.4 (m/* 3.0 x)) (m/cos x18))))])

;;

(defn problem06-bounds
  "Returns the search domain of [[problem06]], the interval `[-10.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem06]], [[dproblem06]]."
  [] [[-10.0 10.0]])

(defn problem06
  "Univariate test function number 06 of the collection of one dimensional global optimization problems: `-(x + sin(x)) exp(-x^2)`.

  Domain: `[-10, 10]`, see [[problem06-bounds]]. Global minimum: `f(x*) = -0.824239 at x* = 0.679560`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem06]] (the derivative)."
  ^double [^double x]
  (m/- (m/* (m/+ x (m/sin x))
            (m/exp (m/- (m/* x x))))))

(defn dproblem06
  "Derivative of the univariate test function [[problem06]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem06]]."
  [[^double x]]
  [(m/* (m/- (m/+ (m/* 2.0 x x) (m/* 2.0 x (m/sin x))) 1.0 (m/cos x))
        (m/exp (m/- (m/* x x))))])

;;

(defn problem07-bounds
  "Returns the search domain of [[problem07]], the interval `[2.7 7.5]`, as a vector with one `[lo hi]` pair.

  See also [[problem07]], [[dproblem07]]."
  [] [[2.7 7.5]])

(defn problem07
  "Univariate test function number 07 of the collection of one dimensional global optimization problems: `sin(x) + sin(10x/3) + ln(x) - 0.84x + 3`.

  Domain: `[2.7, 7.5]`, see [[problem07-bounds]]. Global minimum: `f(x*) = -1.601308 at x* = 5.199780`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem07]] (the derivative)."
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* 3.333333333333333 x))
       (m/log x)
       (m/* -0.84 x)
       3.0))

(defn dproblem07
  "Derivative of the univariate test function [[problem07]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem07]]."
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* 3.333333333333333 (m/cos (m/* 3.333333333333333 x)))
        (m// x)
        -0.84)])

;;

(defn problem08-bounds
  "Returns the search domain of [[problem08]], the interval `[-10.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem08]], [[dproblem08]]."
  [] [[-10.0 10.0]])

(defn problem08
  "Univariate test function number 08 of the collection of one dimensional global optimization problems: `-sum of k cos((k+1)x + k) for k = 1..5`.

  Domain: `[-10, 10]`, see [[problem08-bounds]]. Global minimum: `f(x*) = -14.508008 at x* = -7.083506, -0.800321 and 5.482864` (three minima of the same depth).

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem08]] (the derivative)."
  ^double [^double x]
  (m/- (m/+ (m/cos (m/+ (m/* 2.0 x) 1.0))
            (m/* 2.0 (m/cos (m/+ (m/* 3.0 x) 2.0)))
            (m/* 3.0 (m/cos (m/+ (m/* 4.0 x) 3.0)))
            (m/* 4.0 (m/cos (m/+ (m/* 5.0 x) 4.0)))
            (m/* 5.0 (m/cos (m/+ (m/* 6.0 x) 5.0))))))

(defn dproblem08
  "Derivative of the univariate test function [[problem08]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem08]]."
  [[^double x]]
  [(m/+ (m/* 2.0 (m/sin (m/+ (m/* 2.0 x) 1.0)))
        (m/* 6.0 (m/sin (m/+ (m/* 3.0 x) 2.0)))
        (m/* 12.0 (m/sin (m/+ (m/* 4.0 x) 3.0)))
        (m/* 20.0 (m/sin (m/+ (m/* 5.0 x) 4.0)))
        (m/* 30.0 (m/sin (m/+ (m/* 6.0 x) 5.0))))])

;; 

(defn problem09-bounds
  "Returns the search domain of [[problem09]], the interval `[3.1 20.4]`, as a vector with one `[lo hi]` pair.

  See also [[problem09]], [[dproblem09]]."
  [] [[3.1 20.4]])

(defn problem09
  "Univariate test function number 09 of the collection of one dimensional global optimization problems: `sin(x) + sin(2x/3)`.

  Domain: `[3.1, 20.4]`, see [[problem09-bounds]]. Global minimum: `f(x*) = -1.905961 at x* = 17.039199`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem09]] (the derivative)."
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* m/TWO_THIRDS x))))

(defn dproblem09
  "Derivative of the univariate test function [[problem09]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem09]]."
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* m/TWO_THIRDS (m/cos (m/* m/TWO_THIRDS x))))])

;;

(defn problem10-bounds
  "Returns the search domain of [[problem10]], the interval `[0.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem10]], [[dproblem10]]."
  [] [[0.0 10.0]])

(defn problem10
  "Univariate test function number 10 of the collection of one dimensional global optimization problems: `-x sin(x)`.

  Domain: `[0, 10]`, see [[problem10-bounds]]. Global minimum: `f(x*) = -7.916727 at x* = 7.978666`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem10]] (the derivative)."
  ^double [^double x]
  (m/- (m/* x (m/sin x))))

(defn dproblem10
  "Derivative of the univariate test function [[problem10]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem10]]."
  [[^double x]]
  [(m/- (m/+ (m/sin x) (m/* x (m/cos x))))])

;;

(defn problem11-bounds
  "Returns the search domain of [[problem11]], the interval `[-pi/2 2 pi]`, as a vector with one `[lo hi]` pair.

  See also [[problem11]], [[dproblem11]]."
  [] [[m/-HALF_PI m/TWO_PI]])

(defn problem11
  "Univariate test function number 11 of the collection of one dimensional global optimization problems: `2 cos(x) + cos(2x)`.

  Domain: `[-pi/2, 2 pi]`, see [[problem11-bounds]]. Global minimum: `f(x*) = -1.5 at x* = 2 pi / 3 and x* = 4 pi / 3`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem11]] (the derivative)."
  ^double [^double x]
  (m/+ (m/* 2.0 (m/cos x))
       (m/cos (m/* 2.0 x))))

(defn dproblem11
  "Derivative of the univariate test function [[problem11]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem11]]."
  [[^double x]]
  [(m/- (m/+ (m/* 2.0 (m/sin x))
             (m/* 2.0 (m/sin (m/* 2.0 x)))))])

;;

(defn problem12-bounds
  "Returns the search domain of [[problem12]], the interval `[0 2 pi]`, as a vector with one `[lo hi]` pair.

  See also [[problem12]], [[dproblem12]]."
  [] [[0.0 m/TWO_PI]])

(defn problem12
  "Univariate test function number 12 of the collection of one dimensional global optimization problems: `sin(x)^3 + cos(x)^3`.

  Domain: `[0, 2 pi]`, see [[problem12-bounds]]. Global minimum: `f(x*) = -1 at x* = pi and x* = 3 pi / 2`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem12]] (the derivative)."
  ^double [^double x]
  (let [sx (m/sin x)
        cx (m/cos x)]
    (m/+ (m/* sx sx sx) (m/* cx cx cx))))

(defn dproblem12
  "Derivative of the univariate test function [[problem12]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem12]]."
  [[^double x]]
  [(let [sx (m/sin x)
         cx (m/cos x)]
     (m/- (m/* 3.0 sx sx cx) (m/* 3.0 cx cx sx)))])

;;

(defn problem13-bounds
  "Returns the search domain of [[problem13]], the interval `[0.001 0.99]`, as a vector with one `[lo hi]` pair.

  See also [[problem13]], [[dproblem13]]."
  [] [[0.001 0.99]])

(defn problem13
  "Univariate test function number 13 of the collection of one dimensional global optimization problems: `-x^(2/3) - (1 - x^2)^(1/3)`.

  Domain: `[0.001, 0.99]`, see [[problem13-bounds]]. Global minimum: `f(x*) = -1.587401 at x* = 1 / sqrt(2) = 0.707107`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem13]] (the derivative)."
  ^double [^double x]
  (m/- (m/+ (m/pow x m/TWO_THIRDS)
            (m/pow (m/- 1.0 (m/* x x)) m/THIRD))))

(defn dproblem13
  "Derivative of the univariate test function [[problem13]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem13]]."
  [[^double x]]
  [(m/+ (m/* (m/- m/TWO_THIRDS) (m/pow x (m/- m/THIRD)))
        (m/* m/TWO_THIRDS x (m/pow (m/- 1.0 (m/* x x)) (m/- m/TWO_THIRDS))))])

;;

(defn problem14-bounds
  "Returns the search domain of [[problem14]], the interval `[0.0 4.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem14]], [[dproblem14]]."
  [] [[0.0 4.0]])

(defn problem14
  "Univariate test function number 14 of the collection of one dimensional global optimization problems: `-exp(-x) sin(2 pi x)`.

  Domain: `[0, 4]`, see [[problem14-bounds]]. Global minimum: `f(x*) = -0.788685 at x* = 0.224885`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem14]] (the derivative)."
  ^double [^double x]
  (m/- (m/* (m/exp (m/- x)) (m/sin (m/* m/TWO_PI x)))))

(defn dproblem14
  "Derivative of the univariate test function [[problem14]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem14]]."
  [[^double x]]
  [(let [px (m/* m/TWO_PI x)]
     (m/* (m/exp (m/- x))
          (m/- (m/sin px) (m/* m/TWO_PI (m/cos px)))))])

;;

(defn problem15-bounds
  "Returns the search domain of [[problem15]], the interval `[-5.0 5.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem15]], [[dproblem15]]."
  [] [[-5.0 5.0]])

(defn problem15
  "Univariate test function number 15 of the collection of one dimensional global optimization problems: `(x^2 - 5x + 6) / (x^2 + 1)`.

  Domain: `[-5, 5]`, see [[problem15-bounds]]. Global minimum: `f(x*) = -0.035534 at x* = 1 + sqrt(2) = 2.414214`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem15]] (the derivative)."
  ^double [^double x]
  (let [x2 (m/* x x)]
    (m// (m/+ x2 (m/* -5.0 x) 6.0)
         (m/inc x2))))

(defn dproblem15
  "Derivative of the univariate test function [[problem15]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem15]]."
  [[^double x]]
  [(let [x2 (m/* x x)]
     (m// (m/+ (m/* 5.0 x2) (m/* -10.0 x) -5.0)
          (m/sq (m/inc x2))))])

;;

(defn problem18-bounds
  "Returns the search domain of [[problem18]], the interval `[0.0 6.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem18]], [[dproblem18]]."
  [] [[0.0 6.0]])

(defn problem18
  "Univariate test function number 18 of the collection of one dimensional global optimization problems: `(x - 2)^2 for x <= 3, otherwise 2 ln(x - 2) + 1`.

  Domain: `[0, 6]`, see [[problem18-bounds]]. Global minimum: `f(x*) = 0 at x* = 2`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem18]] (the derivative)."
  ^double [^double x]
  (if (m/<= x 3.0)
    (m/sq (m/- x 2.0))
    (m/+ (m/* 2.0 (m/log (m/- x 2.0))) 1.0)))

(defn dproblem18
  "Derivative of the univariate test function [[problem18]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem18]]."
  [[^double x]]
  [(if (m/<= x 3.0)
     (m/* 2.0 (m/- x 2.0))
     (m// 2.0 (m/- x 2.0)))])

;;

(defn problem20-bounds
  "Returns the search domain of [[problem20]], the interval `[-10.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem20]], [[dproblem20]]."
  [] [[-10.0 10.0]])

(defn problem20
  "Univariate test function number 20 of the collection of one dimensional global optimization problems: `-(x - sin(x)) exp(-x^2)`.

  Domain: `[-10, 10]`, see [[problem20-bounds]]. Global minimum: `f(x*) = -0.063491 at x* = 1.195137`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem20]] (the derivative)."
  ^double [^double x]
  (m/- (m/* (m/- x (m/sin x))
            (m/exp (m/- (m/* x x))))))

(defn dproblem20
  "Derivative of the univariate test function [[problem20]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem20]]."
  [[^double x]]
  (let [x2 (m/* x x)]
    [(m/* (m/+ (m/- (m/* 2.0 x2) (m/* 2.0 x (m/sin x))) -1.0 (m/cos x))
          (m/exp (m/- x2)))]))

;;

(defn problem21-bounds
  "Returns the search domain of [[problem21]], the interval `[0.0 10.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem21]], [[dproblem21]]."
  [] [[0.0 10.0]])

(defn problem21
  "Univariate test function number 21 of the collection of one dimensional global optimization problems: `x sin(x) + x cos(2x)`.

  Domain: `[0, 10]`, see [[problem21-bounds]]. Global minimum: `f(x*) = -9.508350 at x* = 4.795400`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem21]] (the derivative)."
  ^double [^double x]
  (m/+ (m/* x (m/sin x))
       (m/* x (m/cos (m/* 2.0 x)))))

(defn dproblem21
  "Derivative of the univariate test function [[problem21]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem21]]."
  [[^double x]]
  (let [x2 (m/* 2.0 x)]
    [(m/+ (m/sin x)
          (m/* x (m/cos x))
          (m/cos x2)
          (m/* -2.0 x (m/sin x2)))]))

;;

(defn problem22-bounds
  "Returns the search domain of [[problem22]], the interval `[0.0 20.0]`, as a vector with one `[lo hi]` pair.

  See also [[problem22]], [[dproblem22]]."
  [] [[0.0 20.0]])

(defn problem22
  "Univariate test function number 22 of the collection of one dimensional global optimization problems: `exp(-3x) - sin(x)^3`.

  Domain: `[0, 20]`, see [[problem22-bounds]]. Global minimum: `f(x*) = -1 (to about 1e-10) at x* = 5 pi / 2 = 7.853982 and further minima of the same depth`.

  Parameters:

  - `x` (number): the point.

  Returns the value of the function as a double.

  See also [[dproblem22]] (the derivative)."
  ^double [^double x]
  (m/- (m/exp (m/* -3.0 x))
       (m/fpow (m/sin x) 3)))

(defn dproblem22
  "Derivative of the univariate test function [[problem22]].

  Parameters:

  - `[x]` (sequence with one number): the point, in the form used by the optimizers.

  Returns a vector with the derivative at `x` as the only element.

  See also [[problem22]]."
  [[^double x]]
  (let [sx (m/sin x)]
    [(m/- (m/* -3.0 (m/exp (m/* -3.0 x)))
          (m/* 3.0 sx sx (m/cos x)))]))

;; https://www.sfu.ca/~ssurjano/optimization.html

(defn ackley-bounds
  "Returns the usual search domain of the Ackley function [[->ackley]] for `N` dimensions: `[-32.768 32.768]` in every dimension.

  Parameters:

  - `N` (long): number of dimensions.

  Returns a sequence of `N` pairs `[lo hi]`."
  [^long N]
  (repeat N [-32.768 32.768]))

(defn ->ackley
  "Creates the Ackley function for any number of dimensions: `-a exp(-b sqrt(mean(x^2))) - exp(mean(cos(c x))) + a + e`.

  The function has many local minima and one global minimum `f(x*) = 0` at `x* = (0, ..., 0)`. The usual domain is given by [[ackley-bounds]].

  Parameters:

  - `opts` (optional map): `:a` (default: `20.0`), `:b` (default: `0.2`) and `:c` (default: `2 pi`).

  Returns a function of one sequence of numbers (the point) which returns the value as a double."
  ([] (->ackley {:a 20.0 :b 0.2 :c m/TWO_PI}))
  ([{:keys [^double a ^double b ^double c]}]
   (let [ae (m/+ a m/E)
         -b (m/- b)]
     (fn ^double [vs]
       (m/- ae
            (m/* a (m/exp (m/* -b (m/sqrt (m// (v/magsq vs) (double (count vs)))))))
            (m/exp (v/average (v/cos (v/mult vs c)))))))))

;;

(defn bukin-no-6-bounds
  "Returns the usual search domain of [[bukin-no-6]]: `[-15.0 -5.0]` for x1 and `[-3.0 3.0]` for x2."
  [] [[-15.0 -5.0] [-3.0 3.0]])

(defn bukin-no-6
  "Bukin function number 6: `100 sqrt(|x2 - 0.01 x1^2|) + 0.01 |x1 + 10|`.

  The function has a long narrow valley of local minima and is not differentiable on it. Domain: `[-15, -5] x [-3, 3]`, see [[bukin-no-6-bounds]]. Global minimum: `f(x*) = 0` at `x* = (-10, 1)`.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns the value as a double."
  ^double [[^double x1 ^double x2]]
  (m/+ (m/* 100.0 (m/sqrt (m/abs (m/- x2 (m/* 0.01 x1 x1)))))
       (m/* 0.01 (m/abs (m/+ x1 10.0)))))

;;

(defn cross-in-tray-bounds
  "Returns the usual search domain of [[cross-in-tray]]: `[-10.0 10.0]` for x1 and `[-10.0 10.0]` for x2."
  [] [[-10.0 10.0] [-10.0 10.0]])

(defn cross-in-tray
  "Cross-in-tray function: `-0.0001 (|sin(x1) sin(x2) exp(|100 - sqrt(x1^2 + x2^2) / pi|)| + 1)^0.1`.

  The function has four global minima. Domain: `[-10, 10]^2`, see [[cross-in-tray-bounds]]. Global minimum: `f(x*) = -2.06261` at `x* = (1.34941, 1.34941)`, `(-1.34941, 1.34941)`, `(1.34941, -1.34941)` and `(-1.34941, -1.34941)`.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns the value as a double."
  ^double [[^double x1 ^double x2 :as v]]
  (m/* -0.0001 (m/pow (m/inc (m/abs (m/* (m/sin x1)
                                         (m/sin x2)
                                         (m/exp (m/abs (m/- 100.0 (m// (v/mag v) m/PI))))))) 0.1)))

;;

(defn drop-wave-bounds
  "Returns the usual search domain of [[drop-wave]]: `[-5.12 5.12]` for x1 and `[-5.12 5.12]` for x2."
  [] [[-5.12 5.12] [-5.12 5.12]])

(defn drop-wave
  "Drop-wave function: `-(1 + cos(12 sqrt(x1^2 + x2^2))) / (0.5 (x1^2 + x2^2) + 2)`.

  The function is multimodal. Domain: `[-5.12, 5.12]^2`, see [[drop-wave-bounds]]. Global minimum: `f(x*) = -1` at `x* = (0, 0)`.

  Parameters:

  - `v` (sequence of two numbers): the point.

  Returns the value as a double."
  ^double [v]
  (let [ms (v/magsq v)]
    (m/- (m// (m/inc (m/cos (m/* 12.0 (m/sqrt ms))))
              (m/+ (m/* 0.5 ms) 2.0)))))

;;

(defn egg-holder-bounds
  "Returns the usual search domain of [[egg-holder]]: `[-512.0 512.0]` for x1 and `[-512.0 512.0]` for x2."
  [] [[-512.0 512.0] [-512.0 512.0]])

(defn egg-holder
  "Egg-holder function: `-(x2 + 47) sin(sqrt(|x2 + x1 / 2 + 47|)) - x1 sin(sqrt(|x1 - (x2 + 47)|))`.

  The function is very rugged with many local minima. Domain: `[-512, 512]^2`, see [[egg-holder-bounds]]. Global minimum: `f(x*) = -959.6407` at `x* = (512, 404.2319)`.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns the value as a double."
  ^double [[^double x1 ^double x2]]
  (let [x47 (m/+ x2 47.0)]
    (m/- (m/* -1.0 x47 (m/sin (m/sqrt (m/abs (m/+ x47 (m/* 0.5 x1))))))
         (m/* x1 (m/sin (m/sqrt (m/abs (m/- x1 x47))))))))

;;

(defn gramacy-lee-2012-bounds [] [[0.5 2.5]])

(defn gramacy-lee-2012
  ^double [^double x]
  (m/+ (m// (m/sinpi (m/* 10.0 x)) (m/* 2.0 x))
       (m/fpow (m/dec x) 4)))

;;

(defn griewank-bounds [^long N] (repeat N [-600.0 600.0]))

(defn griewank
  ^double [v]
  (m/inc (m/- (m// (v/magsq v) 4000.0)
              (v/prod (map-indexed (fn [^long i ^double x] (m/cos (m// x (m/sqrt (m/inc i))))) v)))))

;;

(defn holder-table-bounds [] [[-10.0 10.0] [-10.0 10.0]])

(defn holder-table
  ^double [[^double x1 ^double x2 :as v]]
  (m/- (m/abs (m/* (m/sin x1) (m/sin x2) (m/exp (m/abs (m/- 1.0 (m// (v/mag v) m/PI))))))))

;;

(defn ->langerman-general
  [c A-rows]
  (when-not (m/== (count c) (count A-rows)) (throw (ex-info "c and A dimensionality do not agree."
                                                            {:m (count c) :shape-A [(count A-rows)
                                                                                    (count (first A-rows))]})))
  (fn ^double [v]
    (v/dot c (map (fn [row]
                    (let [d (v/magsq (v/sub row v))]
                      (m/* (m/exp (m/- (m// d m/PI)))
                           (m/cospi d)))) A-rows))))

(defn langerman-bounds [^long N] (repeat N [0.0 10.0]))

;;  As defined in Molga/Smutnicki
(def langerman (->langerman-general [1.0 2.0 5.0 2.0 3.0]
                                    [[3.0 5.0] [5.0 2.0] [2.0 1.0] [1.0 4.0] [7.0 9.0]]))

;; (sceua/sceua langerman {:bounds (langerman-bounds 2)})
;; => [[2.7933531577960338 1.5972900206456435] -4.155809123348931]

;;

(def ^:private langerman-A
  [[9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020]
   [9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374]
   [8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982]
   [2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426]
   [8.074 8.777 3.467 1.863 6.708 6.349 4.534 0.276 7.633 1.567]
   [7.650 5.658 0.720 2.764 3.278 5.283 7.474 6.274 1.409 8.208]
   [1.256 3.605 8.623 6.905 0.584 8.133 6.071 6.888 4.187 5.448]
   [8.314 2.261 4.224 1.781 4.124 0.932 8.129 8.658 1.208 5.762]
   [0.226 8.858 1.420 0.945 1.622 4.698 6.228 9.096 0.972 7.637]
   [0.305 2.228 1.242 5.928 9.133 1.826 4.060 5.204 8.713 8.247]
   [0.652 7.027 0.508 4.876 8.807 4.632 5.808 6.937 3.291 7.016]
   [2.699 3.516 5.874 4.119 4.461 7.496 8.817 0.690 6.593 9.789]
   [8.327 3.897 2.017 9.570 9.825 1.150 1.395 3.885 6.354 0.109]
   [2.132 7.006 7.136 2.641 1.882 5.943 7.273 7.691 2.880 0.564]
   [4.707 5.579 4.080 0.581 9.698 8.542 8.077 8.515 9.231 4.670]
   [8.304 7.559 8.567 0.322 7.128 8.392 1.472 8.524 2.277 7.826]
   [8.632 4.409 4.832 5.768 7.050 6.715 1.711 4.323 4.405 4.591]
   [4.887 9.112 0.170 8.967 9.693 9.867 7.508 7.770 8.382 6.740]
   [2.440 6.686 4.299 1.007 7.008 1.427 9.398 8.480 9.950 1.675]
   [6.306 8.583 6.084 1.138 4.350 3.134 7.853 6.061 7.457 2.258]
   [0.652 2.343 1.370 0.821 1.310 1.063 0.689 8.819 8.833 9.070]
   [5.558 1.272 5.756 9.857 2.279 2.764 1.284 1.677 1.244 1.234]
   [3.352 7.549 9.817 9.437 8.687 4.167 2.570 6.540 0.228 0.027]
   [8.798 0.880 2.370 0.168 1.701 3.680 1.231 2.390 2.499 0.064]
   [1.460 8.057 1.336 7.217 7.914 3.615 9.981 9.198 5.292 1.224]
   [0.432 8.645 8.774 0.249 8.081 7.461 4.416 0.652 4.002 4.644]
   [0.679 2.800 5.523 3.049 2.968 7.225 6.730 4.199 9.614 9.229]
   [4.263 1.074 7.286 5.599 8.291 5.200 9.214 8.272 4.398 4.506]
   [9.496 4.830 3.150 8.270 5.079 1.231 5.731 9.494 1.883 9.732]
   [4.138 2.562 2.532 9.661 5.611 5.500 6.886 2.341 9.699 6.500]])

(def ^:private langerman-c [0.806,0.517,1.5,0.908,0.965,0.669,0.524,0.902,0.531,0.876,0.462,
                            0.491,0.463,0.714,0.352,0.869,0.813,0.811,0.828,0.964,0.789,
                            0.360,0.369,0.992,0.332,0.817,0.632,0.883,0.608,0.326])

(defn ->langerman
  "As defined in Adorio, 1<=N<=10; 1<=m<=30"
  [^long N ^long m]
  (let [c (subvec langerman-c 0 m)
        A-rows (take m (map (fn [row] (subvec row 0 N)) langerman-A))]
    (->langerman-general c A-rows)))

;;

(defn levy-bounds [^long N] (repeat N [-10.0 10.0]))

(defn levy
  ^double [v]
  (let [[w1 & r] (v/mult (v/shift v 3.0) 0.25)
        t1 (m/sq (m/sinpi w1))]
    (if-not (seq r)
      t1
      (let [wd (double (last r))
            td (m/* (m/sq (m/dec wd)) (m/inc (m/sq (m/sinpi (m/* 2.0 wd)))))]
        (m/+ t1 td (v/sum (map (fn ^double [^double wi]
                                 (m/* (m/sq (m/dec wi))
                                      (m/inc (m/* 10.0 (m/sq (m/sin (m/inc (m/* wi m/PI)))))))) (butlast r))))))))

;;

(defn levy-13-bounds [] [[-10.0 10.0] [-10.0 10.0]])

(defn levy-13
  ^double [[^double x1 ^double x2]]
  (m/+ (m/sq (m/sin (m/* 3.0 m/PI x1)))
       (m/* (m/sq (m/dec x1)) (m/inc (m/sq (m/sin (m/* 3.0 m/PI x2)))))
       (m/* (m/sq (m/dec x2)) (m/inc (m/sq (m/sin (m/* 2.0 m/PI x2)))))))

;;

(defn rastrigin-bounds [^long N] (repeat N [-5.12 5.12]))

(defn rastrigin
  ^double [v]
  (m/+ (m/* 10.0 (count v))
       (v/sum (map (fn [^double xi]
                     (m/- (m/* xi xi) (m/* 10.0 (m/cos (m/* m/TWO_PI xi))))) v))))


;;

(defn schaffer-2-bounds [] [[-100.0 100.0] [-100.0 100.0]])

(defn schaffer-2
  ^double [[^double x1 ^double x2]]
  (m/+ 0.5 (m// (m/- (m/sq (m/sin (m/- (m/* x1 x1) (m/* x2 x2)))) 0.5)
                (m/sq (m/inc (m/* 0.001 (m/+ (m/* x1 x1) (m/* x2 x2))))))))

;;

(defn schaffer-4-bounds [] [[-100.0 100.0] [-100.0 100.0]])

(defn schaffer-4
  ^double [[^double x1 ^double x2]]
  (m/+ 0.5 (m// (m/- (m/sq (m/cos (m/sin (m/abs (m/- (m/* x1 x1) (m/* x2 x2)))))) 0.5)
                (m/sq (m/inc (m/* 0.001 (m/+ (m/* x1 x1) (m/* x2 x2))))))))

;;

(defn schwefel-bounds [^long N] (repeat N [-500.0 500.0]))

(defn schwefel
  ^double [v]
  (m/- (m/* 418.9829 (count v))
       (v/sum (map (fn [^double xi] (m/* xi (m/sin (m/sqrt (m/abs xi))))) v))))

;;

(defn shubert-bounds [] [[-5.12 5.12] [-5.12 5.12]])

(defn shubert
  ^double [[^double x1 ^double x2]]
  (m/* (m/+ (m/cos (m/inc (m/* 2.0 x1)))
            (m/* 2.0 (m/cos (m/inc (m/* 3.0 x1))))
            (m/* 3.0 (m/cos (m/inc (m/* 4.0 x1))))
            (m/* 4.0 (m/cos (m/inc (m/* 5.0 x1))))
            (m/* 5.0 (m/cos (m/inc (m/* 6.0 x1)))))
       (m/+ (m/cos (m/inc (m/* 2.0 x2)))
            (m/* 2.0 (m/cos (m/inc (m/* 3.0 x2))))
            (m/* 3.0 (m/cos (m/inc (m/* 4.0 x2))))
            (m/* 4.0 (m/cos (m/inc (m/* 5.0 x2))))
            (m/* 5.0 (m/cos (m/inc (m/* 6.0 x2)))))))

;;

(defn bohachevsky-bounds [] [[-100.0 100.0] [-100.0 100.0]])

(defn bohachevsky-1
  ^double [[^double x1 ^double x2]]
  (m/- (m/+ (m/* x1 x1) (m/* 2.0 x2 x2) 0.7)
       (m/* 0.3 (m/cos (m/* 3.0 m/PI x1)))
       (m/* 0.4 (m/cos (m/* 4.0 m/PI x2)))))

(defn bohachevsky-2
  ^double [[^double x1 ^double x2]]
  (m/- (m/+ (m/* x1 x1) (m/* 2.0 x2 x2) 0.3)
       (m/* 0.3 (m/cos (m/* 3.0 m/PI x1)) (m/cos (m/* 4.0 m/PI x2)))))

(defn bohachevsky-3
  ^double [[^double x1 ^double x2]]
  (m/- (m/+ (m/* x1 x1) (m/* 2.0 x2 x2) 0.3)
       (m/* 0.3 (m/cos (m/+ (m/* 3.0 m/PI x1)
                            (m/* 4.0 m/PI x2))))))

;;

(defn perm-0db-bounds [^long N] (repeat N [(m/- N) N]))

(defn ->perm-0db
  [^double beta]
  (fn ^double [v]
    (v/sum (for [^long i (range (count v))
                 :let [i+ (m/inc i)]]
             (m/sq (v/sum (map-indexed (fn [^long j ^double xj]
                                         (let [j+ (m/inc j)]
                                           (m/* (m/+ j+ beta)
                                                (m/- (m/fpow xj i+)
                                                     (m// 1.0 (m/fpow j+ i+)))))) v)))))))

;;

(defn rotated-hyper-ellipsoid-bounds [^long N] (repeat N [-65.536 65.536]))

(defn rotated-hyper-ellipsoid
  ^double [v]
  (let [v (vec v)]
    (v/sum (for [^long i (range 1 (m/inc (count v)))
                 ^long j (range i)]
             (m/sq (v j))))))

;;

(defn sphere-bounds
  "Returns the usual search domain of [[sphere]] for `N` dimensions: `[-5.12 5.12]` in every dimension.

  Parameters:

  - `N` (long): number of dimensions.

  Returns a sequence of `N` pairs `[lo hi]`."
  [^long N]
  (repeat N [-5.12 5.12]))

(defn sphere
  "Sphere function for any number of dimensions: the sum of squares of the coordinates.

  The function is convex with one global minimum. Domain: `[-5.12, 5.12]^N`, see [[sphere-bounds]]. Global minimum: `f(x*) = 0` at `x* = (0, ..., 0)`.

  Parameters:

  - `v` (sequence of numbers): the point.

  Returns the value as a double.

  See also [[sphere-gradient]]."
  ^double [v]
  (v/magsq v))

(defn sphere-gradient
  "Gradient of the [[sphere]] function: twice the point.

  Parameters:

  - `v` (sequence of numbers): the point.

  Returns the partial derivatives in the form of the argument (a vector for a vector, an array for an array)."
  [v]
  (v/mult v 2.0))

;;

(defn sum-of-different-powers-bounds [^long N] (repeat N [-1.0 1.0]))

(defn sum-of-different-powers
  ^double [v]
  (v/sum (map-indexed (fn [^long i ^double xi]
                        (m/fpow (m/abs xi) (m/long-add i 2))) v)))

;;

(defn sum-squares-bounds [^long N] (repeat N [-10.0 10.0]))

(defn sum-squares
  ^double [v]
  (v/sum (map-indexed (fn [^long i ^double xi]
                        (m/* (m/inc i) xi xi)) v)))


;;

(defn trid-bounds [^long N] (let [n2 (m/sq N)] (repeat N [(m/- n2) n2])))

(defn trid
  ^double [v]
  (m/- (v/sum (map (fn [^double x] (m/sq (m/dec x))) v))
       (v/sum (map (fn [[^double xp ^double xn]] (m/* xp xn)) (partition 2 1 v)))))

;;

(defn booth-bounds [] [[-10.0 10.0] [-10.0 10.0]])

(defn booth
  ^double [[^double x1 ^double x2]]
  (m/+ (m/sq (m/+ x1 (m/* 2.0 x2) -7.0))
       (m/sq (m/+ (m/* 2.0 x1) x2 -5.0))))

;;

(defn matyas-bounds [] [[-10.0 10.0] [-10.0 10.0]])

(defn matyas
  ^double [[^double x1 ^double x2]]
  (m/- (m/* 0.26 (m/+ (m/* x1 x1) (m/* x2 x2)))
       (m/* 0.48 x1 x2)))

;;

(defn mccormick-bounds [] [[-1.5 4.0] [-3.0 4.0]])

(defn mccormick 
  ^double [[^double x1 ^double x2]]
  (m/inc (m/+ (m/sin (m/+ x1 x2))
              (m/sq (m/- x1 x2))
              (m/* -1.5 x1)
              (m/* 2.5 x2))))


;;

(defn power-sum-bounds [^long N] (repeat N [0.0 N]))

(defn ->power-sum
  [b]
  (fn ^double [v]
    (v/sum (map-indexed (fn [^long i ^double bi]
                          (m/sq (m/- (v/sum (v/pow v (m/inc i))) bi))) b))))


(defn power-sum-4-bounds [] (power-sum-bounds 4))

(def power-sum-4 (->power-sum [8.0 18.0 44.0 114.0]))

;;

(defn zakharov-bounds [^long N] (repeat N [-5.0 10.0]))

(defn zakharov
  ^double [v]
  (let [s (v/sum (map-indexed (fn [^long i ^double xi]
                                (m/* 0.5 (m/inc i) xi)) v))
        ss (m/sq s)]
    (m/+ (v/magsq s) ss (m/sq ss))))

;;

(defn three-hump-camel-bounds [] [[-5.0 5.0] [-5.0 5.0]])

(defn three-hump-camel
  ^double [[^double x1 ^double x2]]
  (let [x12 (m/* x1 x1)
        x14 (m/* x12 x12)
        x16 (m/* x12 x14)]
    (m/+ (m/* 2.0 x12)
         (m/* -1.05 x14)
         (m/* m/SIXTH x16)
         (m/* x1 x2)
         (m/* x2 x2))))

;;

(defn six-hump-camel-bounds [] [[-3.0 3.0] [-2.0 2.0]])

(defn six-hump-camel
  ^double [[^double x1 ^double x2]]
  (let [x12 (m/* x1 x1)
        x14 (m/* x12 x12)
        x22 (m/* x2 x2)]
    (m/+ (m/* x12 (m/+ 4.0
                       (m/* -2.1 x12)
                       (m/* m/THIRD x14)))
         (m/* x1 x2)
         (m/* x22 (m/+ -4.0
                       (m/* 4.0 x22))))))

;;

(defn dixon-price-bounds [^long N] (repeat N [-10.0 10.0]))

(defn dixon-price
  ^double [[^double x1 :as v]]
  (m/+ (m/sq (m/dec x1))
       (->> (partition 2 1 v)
            (map-indexed (fn [^long i [^double xp ^double xn]]
                           (m/* (m/+ i 2)
                                (m/sq (m/- (m/* 2.0 xn xn) xp)))))
            (v/sum))))

;;

(defn rosenbrock-bounds
  "Returns the usual search domain of [[rosenbrock]] for `N` dimensions: `[-5 10]` in every dimension.

  Parameters:

  - `N` (long): number of dimensions.

  Returns a sequence of `N` pairs `[lo hi]`. A smaller domain `[-2.048, 2.048]` is also common."
  [^long N]
  (repeat N [-5.0 10.0]))

(defn rosenbrock
  "Rosenbrock function for any number of dimensions: the sum of `100 (x[i+1] - x[i]^2)^2 + (x[i] - 1)^2` for consecutive pairs of coordinates.

  The function has a long curved valley and one global minimum (for more than three dimensions also a local one close to `(-1, 1, ..., 1)`). Domain: `[-5, 10]^N`, see [[rosenbrock-bounds]]. Global minimum: `f(x*) = 0` at `x* = (1, ..., 1)`.

  Parameters:

  - `v` (sequence of numbers): the point.

  Returns the value as a double. A point with fewer than two coordinates gives `0.0`.

  See also [[rosenbrock-gradient]]."
  ^double [v]
  (->> (partition 2 1 v)
       (map (fn [[^double x1 ^double x2]]
              (m/+ (m/* 100.0 (m/sq (m/- x2 (m/* x1 x1))))
                   (m/sq (m/dec x1)))))
       (v/sum)))

(defn rosenbrock-gradient
  "Gradient of the [[rosenbrock]] function.

  Parameters:

  - `v` (sequence of numbers): the point.

  Returns a vector of doubles with the partial derivatives, one for every coordinate. A point with fewer than two coordinates gives zeros."
  [v]
  (let [xs (m/seq->double-array v)
        n (alength xs)
        g (double-array n)]
    (dotimes [i (m/dec n)]
      (let [x (aget xs i)
            x-next (Array/aget xs (m/inc i))
            d (m/- x-next (m/* x x))]
        (Array/aset g i (m/+ (Array/aget g i) (m/* -400.0 x d) (m/* 2.0 (m/dec x))))
        (Array/aset g (m/inc i) (m/+ (Array/aget g (m/inc i)) (m/* 200.0 d)))))
    g))

;;

(defn de-jong-5-bounds [] [[-65.536 65.536] [-65.536 65.536]])

(def ^:private de-jong-a [[-32.0 -32.0] [-16.0 -32.0] [0.0 -32.0] [16.0 -32.0] [32.0 -32.0] [-32.0 -16.0] [-16.0 -16.0] [0.0 -16.0] [16.0 -16.0] [32.0 -16.0] [-32.0 0.0] [-16.0 0.0] [0.0 0.0] [16.0 0.0] [32.0 0.0] [-32.0 16.0] [-16.0 16.0] [0.0 16.0] [16.0 16.0] [32.0 16.0] [-32.0 32.0] [-16.0 32.0] [0.0 32.0] [16.0 32.0] [32.0 32.0]])

(defn de-jong-5
  ^double [[^double x1 ^double x2]]
  (m// (m/+ 0.002 (->> de-jong-a
                       (map-indexed (fn [^long i [^double a1 ^double a2]]
                                      (m// (m/+ 1 i (m/fpow (m/- x1 a1) 6) (m/fpow (m/- x2 a2) 6)))))
                       (v/sum)))))

;;

(defn easom-bounds [] [[-100.0 100.0] [-100.0 100.0]])

(defn easom
  ^double [[^double x1 ^double x2]]
  (m/- (m/* (m/cos x1)
            (m/cos x2)
            (m/exp (m/- (m/- (m/sq (m/- x1 m/PI)))
                        (m/sq (m/- x2 m/PI)))))))

;;

(defn michalewicz-bounds [^long N] (repeat N [0.0 m/PI]))

(defn ->michalewicz
  [^double m]
  (let [m2 (m/* 2.0 m)]
    (fn ^double [v]
      (m/- (v/sum (map-indexed (fn [^long i ^double xi]
                                 (m/* (m/sin xi)
                                      (m/pow (m/sin (m// (m/* (m/inc i) xi xi) m/PI)) m2))) v))))))

(def michalewicz (->michalewicz 10.0))

;;

(defn beale-bounds
  "Returns the usual search domain of [[beale]]: `[-4.5 4.5]` for x1 and `[-4.5 4.5]` for x2."
  [] [[-4.5 4.5] [-4.5 4.5]])

(defn beale
  "Beale function: `(1.5 - x1 + x1 x2)^2 + (2.25 - x1 + x1 x2^2)^2 + (2.625 - x1 + x1 x2^3)^2`.

  The function has sharp peaks at the corners of the domain and flat valleys. Domain: `[-4.5, 4.5]^2`, see [[beale-bounds]]. Global minimum: `f(x*) = 0` at `x* = (3, 0.5)`.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns the value as a double.

  See also [[beale-gradient]]."
  ^double [[^double x1 ^double x2]]
  (let [x1x2 (m/* x1 x2)
        x1x22 (m/* x1x2 x2)]
    (m/+ (m/sq (m/+ (m/- 1.5 x1) x1x2))
         (m/sq (m/+ (m/- 2.25 x1) x1x22))
         (m/sq (m/+ (m/- 2.625 x1) (m/* x1x22 x2))))))

(defn beale-gradient
  "Gradient of the [[beale]] function.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns a vector of two doubles with the partial derivatives."
  [[^double x1 ^double x2]]
  (let [x22 (m/* x2 x2)
        x23 (m/* x22 x2)
        a (m/+ (m/- 1.5 x1) (m/* x1 x2))
        b (m/+ (m/- 2.25 x1) (m/* x1 x22))
        c (m/+ (m/- 2.625 x1) (m/* x1 x23))]
    [(m/* 2.0 (m/+ (m/* a (m/dec x2)) (m/* b (m/dec x22)) (m/* c (m/dec x23))))
     (m/* 2.0 x1 (m/+ a (m/* 2.0 x2 b) (m/* 3.0 x22 c)))]))

;;

(defn branin-bounds [] [[-5.0 10.0] [0.0 15.0]])

(defn ->branin
  [[^double a ^double b ^double c ^double r ^double s ^double t]]
  (let [st (m/* s (m/- 1.0 t))]
    (fn ^double [[^double x1 ^double x2]]
      (m/+ (m/* a (m/sq (m/- (m/+ x2 (m/* c x1))
                             (m/* b x1 x1)
                             r)))
           (m/* st (m/cos x1))
           s))))

(def branin (->branin [1.0 (m// 5.1 (m/* 4.0 m/PI2)) (m// 5.0 m/PI)
                       6.0 10.0 (m// (m/* 8.0 m/PI))]))

;;

(defn colville-bounds [] (repeat 4 [-10.0 10.0]))

(defn colville
  ^double [[^double x1 ^double x2 ^double x3 ^double x4]]
  (let [x1- (m/dec x1)
        x2- (m/dec x2)
        x3- (m/dec x3)
        x4- (m/dec x4)]
    (m/+ (m/* 100.0 (m/sq (m/- (m/* x1 x1) x2)))
         (m/* x1- x1-)
         (m/* x3- x3-)
         (m/* 90.0 (m/sq (m/- (m/* x3 x3) x4)))
         (m/* 10.1 (m/+ (m/* x2- x2-)
                        (m/* x4- x4-)))
         (m/* 19.8 x2- x4-))))

;;

(defn forrester-2018-bounds [] [[0.0 1.0]])

(defn forrester-2018
  ^double [[^double x]]
  (m/* (m/sq (m/- (m/* 6.0 x) 2.0))
       (m/sin (m/- (m/* 12.0 x) 4.0))))

;;

(defn goldstein-price-bounds [] [[-2.0 2.0] [-2.0 2.0]])

(defn goldstein-price
  ^double [[^double x1 ^double x2]]
  (let [x12 (m/* x1 x1)
        x22 (m/* x2 x2)]
    (m/* (m/inc (m/* (m/sq (m/+ x1 x2 1.0))
                     (m/+ 19.0
                          (m/* -14.0 x1)
                          (m/* 3.0 x12)
                          (m/* -14.0 x2)
                          (m/* 6.0 x1 x2)
                          (m/* 3.0 x22))))
         (m/+ 30.0 (m/* (m/sq (m/- (m/* 2.0 x1)
                                   (m/* 3.0 x2)))
                        (m/+ 18.0
                             (m/* -32.0 x1)
                             (m/* 12.0 x12)
                             (m/* 48.0 x2)
                             (m/* -36.0 x1 x2)
                             (m/* 27.0 x22)))))))

;;

(defn perm-db-bounds [^long N] (repeat N [(m/- N) N]))

(defn ->perm-db
  [^double beta]
  (fn ^double [v]
    (v/sum (for [^long i (range (count v))
                 :let [i+ (m/inc i)]]
             (m/sq (v/sum (map-indexed (fn [^long j ^double xj]
                                         (let [j+ (m/inc j)]
                                           (m/* (m/+ (m/fpow j+ i+) beta)
                                                (m/dec (m/fpow (m// xj j+) i+))))) v)))))))

;;

(defn styblinski-tang-bounds [^long N] (repeat N [-5.0 5.0]))

(defn styblinski-tang
  ^double [v]
  (m/* 0.5 (v/sum (map (fn [^double xi]
                         (let [xi2 (m/* xi xi)]
                           (m/+ (m/* xi2 xi2)
                                (m/* -16.0 xi2)
                                (m/* 5.0 xi)))) v))))

;;

(defn himmelblau-bounds
  "Returns the usual search domain of [[himmelblau]]: `[-5.0 5.0]` for x1 and `[-5.0 5.0]` for x2."
  [] [[-5.0 5.0] [-5.0 5.0]])

(defn himmelblau
  "Himmelblau function: `(x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2`.

  The function has four global minima. Domain: `[-5, 5]^2`, see [[himmelblau-bounds]]. Global minimum: `f(x*) = 0` at `x* = (3, 2)`, `(-2.805118, 3.131312)`, `(-3.779310, -3.283186)` and `(3.584428, -1.848126)`.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns the value as a double.

  See also [[himmelblau-gradient]]."
  ^double [[^double x1 ^double x2]]
  (m/+ (m/sq (m/+ (m/* x1 x1) x2 -11.0))
       (m/sq (m/+ x1 (m/* x2 x2) -7.0))))

(defn himmelblau-gradient
  "Gradient of the [[himmelblau]] function.

  Parameters:

  - `[x1 x2]` (sequence of two numbers): the point.

  Returns a vector of two doubles with the partial derivatives."
  [[^double x1 ^double x2]]
  (let [a (m/+ (m/* x1 x1) x2 -11.0)
        b (m/+ x1 (m/* x2 x2) -7.0)]
    [(m/+ (m/* 4.0 x1 a) (m/* 2.0 b))
     (m/+ (m/* 2.0 a) (m/* 4.0 x2 b))]))
