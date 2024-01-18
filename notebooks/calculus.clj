^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true}
(ns calculus
  {:clj-kondo/config '{:config-in-call {utils/table2 {:ignore [:unresolved-symbol]}}}}
  (:require [fastmath.calculus :as calc]
            [fastmath.solver :as solver]
            [nextjournal.clerk :as clerk]
            [utils :as u]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.n :as n]))

;; # Calculus

;; Integration, differentiation (finite difference), solvers and Richardson extrapolation.

#_:clj-kondo/ignore
^{::clerk/visibility {:result :hide}}
(require '[fastmath.calculus :as calc])

;; ## Richardson extrapolation

;; [Richardson extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation) is a way to increase accuracy and convergence of limits of the form:

;; $$f(x) = \lim_{h\to 0}g(x,h)$$

;; If $x = \infty$ the following limit is calculated (similar for negative infinity):

;; $$f(\infty) = \lim_{h\to 0}g(0,\frac{1}{h})$$

;; There are two helpers which create $g(x,h)$ from any single arity function: $g_+(x,h) = f(x+h)$ and $g_-(x,h) = f(x-h)$.

;; List of the functions:

^{::clerk/visibility :hide}
(u/table2
 [[extrapolate "Create extrapolated function f(x) from a function g(x,h)"]
  [fx->gx+h "Create g+"]
  [fx->gx-h "Create g-"]])

;; `extrapolate` accepts following options:

;; * `:contract` - `h` shrinkage factor, default=`1/2`
;; * `:power` - set to `2.0` for even functions around `x0`, default `1.0`
;; * `:init-h` - initial step `h`, default=`1/2`
;; * `:abs` - absolute error, default: machine epsilon
;; * `:rel` - relative error, default: ulp for init-h
;; * `:tol` - tolerance for current error, default: `2.0` 
;; * `:max-evals` - maximum evaluations, default: maximum integer

;; ### Example - finite difference

;; As the first example let's see how extrapolation increases accuracy for finite differences method. Let's define.

;; $$f(x)=x\sin{x}$$

^{::clerk/visibility {:result :hide}}
(defn xsinx [x] (* x (m/sin x)))

;; And actual derivative

;; $$f'(x)=\sin{x}+x\cos{x}$$

^{::clerk/visibility {:result :hide}}
(defn xsinx' [x] (+ (m/sin x) (* x (m/cos x))))

;; Value of derivative for $x=0.5$ equals:

(xsinx' 0.5)

;; Now let's define finite difference version of derivative

;; $$f'(x,h) \approx f'_{fd}(x,h) = \frac{f(x+h)-f(x)}{h}$$

^{::clerk/visibility {:result :hide}}
(defn xsinx-finite-diff [x h] (/ (- (xsinx (+ x h)) (xsinx x)) h))

;; For low `h` the result is close to the actual value:

(xsinx-finite-diff 0.5 1.0e-6)

;; We list what is the error when `h` goes to `0.0`. Observe, that the best results is for `h` equal `1.0e-6` with error around `1.0e-8`.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["h" "finite difference" "error"]
  :rows (map (fn [p]
               (let [h (double (.divide (bigdec 1.0) (.pow (bigdec 10.0) p)))
                     v (xsinx-finite-diff 0.5 h)]
                 [h v (m/abs (- v (xsinx' 0.5)))])) (range 18))})

;; Now, let's construct extrapolated version of our finite difference derivative

^{::clerk/visibility {:result :hide}}
(def xsinx-finite-diff+ (calc/extrapolate xsinx-finite-diff))

;; Extrapolated value is:

(xsinx-finite-diff+ 0.5)

;; with the error close to `6.0e-15`, which is much better than the default approach.

(- (xsinx-finite-diff+ 0.5) (xsinx' 0.5))

;; ### Example - infinite sum

;; In the second example we will try to calculate the following infinite sum:

;; $$\sum_{n=1}^\infty\frac{1}{n^2} = \frac{\pi^2}{6}$$

^{::clerk/visibility {:result :hide}}
(defn infinite-sum
  [n]
  (->> (range 1 (inc n))
       (map (fn [^double x] (/ 1.0 (* x x))))
       (reduce + )))

;; We list value and error for various `n`. Worth to mention that the last position for `n=1.0e9` takes around 20 seconds on my PC.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["n" "infinite sum" "error"]
  :rows (map (fn [p]
               (let [h (unchecked-long (m/pow 10 p))
                     v (infinite-sum h)]
                 [h v (m/abs (- v (m// (m/sq m/PI) 6.0)))])) (range 1 9))})

;; Now let's evaluate extrapolated version at infinity. Relative tolerance is set to $0.0$ to stop only on absolute tolerance.

^{::clerk/visibility {:result :hide}}
(def infinite-sum+ (-> infinite-sum
                     (calc/fx->gx+h) ;; convert to g(x,h)=f(x+h)
                     (calc/extrapolate {:rel 0})))

;; From the debugging we know that evaluation stops for `n=2048`

(infinite-sum+ ##Inf)

;; and the error is around `5.0e-15`

(- (infinite-sum+ ##Inf) (/ (* m/PI m/PI) 6.0))

;; ### Example - limit

;; In another example we will calculate the following limit:

;; $$\lim_{x\to \infty}\frac{\ln(x)}{x} = 0$$

^{::clerk/visibility {:result :hide}}
(defn limit-example [x] (/ (m/ln x) x))

(limit-example Double/MAX_VALUE)

;; First approach shows that this limit behaves badly when extrapolated. The result is far from beeing accurate.

^{::clerk/visibility {:result :hide}}
(def limit-example+-bad (-> limit-example
                          (calc/fx->gx+h)
                          (calc/extrapolate)))

(limit-example+-bad ##Inf)

;; To fix it we need to raise the initial `h` and turn off the absolute tolerance checking (and rely only on the relative tolerance)

^{::clerk/visibility {:result :hide}}
(def limit-example+-good (-> limit-example
                           (calc/fx->gx+h)
                           (calc/extrapolate {:init-h 10 :abs 0})))

(limit-example+-good ##Inf)

;; ## Differentiation

;; Differentiation is based on [finite difference](https://en.wikipedia.org/wiki/Finite_difference) method. Forward, backward and central difference methods are possible. Any order with any [accuracy coefficients](https://en.wikipedia.org/wiki/Finite_difference_coefficient) are supported.

;; Please note that higher orders and higher accuracy can be not stable.

^{::clerk/visibility :hide}
(u/table2
 [[derivative "Generate nth order derivative"]
  [f' "First derivative"]
  [f'' "Second derivative"]
  [f''' "Third derivative"]
  [gradient "Grandient of the scalar valued function of several variables"]
  [hessian "Matrix of second derivatives of scalar valued function of several variables"]])

;; ### Derivative

;; `derivative` creates a `n`th order (default: first) derivative of the input function as another function. Possible options are:
;;
;; * `:h` - step, default: `0.0` - automatic
;; * `:acc` - order of accuracy, default: `2`
;; * `:method` - one of: `:forward`, `:backward` or `:central` (default)
;; * `:extrapolate?` - if true creates an extrapolated version, optionally a map of `extrapolate` options can be passed

^{::clerk/visibility {:result :hide}}
(def xsinx-diff (calc/derivative xsinx))

(xsinx-diff 0.5)

;; Let's list all the possible values for first derivative and different accuracy, method and extrapolation

^{::clerk/visibility :hide}
(clerk/table
 {:head ["accuracy" "method" "extrapolation?" "f'(0.5)" "error"]
  :rows (let [res (xsinx' 0.5)]
          (for [acc [2 4 6]
                m [:forward :backward :central]
                extr? [false true]
                :let [d (calc/derivative xsinx 1 {:acc acc :extrapolate? extr? :method m})
                      dres (d 0.5)]]
            [acc m extr? dres (m/abs (- res dres))]))})

;; Similar for third derivative

^{::clerk/visibility {:result :hide}}
(defn xsinx''' [x] (- (* -3.0 (m/sin x))
                   (* x (m/cos x))))

(xsinx''' 0.5)

^{::clerk/visibility {:result :hide}}
(def xsinx-diff-3 (calc/derivative xsinx 3))

(xsinx-diff-3 0.5)

^{::clerk/visibility :hide}
(clerk/table
 {:head ["accuracy" "method" "extrapolation?" "f'''(0.5)" "error"]
  :rows (let [res (xsinx''' 0.5)]
          (for [acc [2 4 6]
                m [:forward :backward :central]
                extr? [false true]
                :let [d (calc/derivative xsinx 3 {:acc acc :extrapolate? extr? :method m})
                      dres (d 0.5)]]
            [acc m extr? dres (m/abs (- res dres))]))})

;; ### Gradient and Hessian

;; For multivariate real-valued function there are two functions: `gradient` and `hessian`. Both functions accept `:h` option, `gradient` also `:acc` (2 or 4).

^{::clerk/visibility :hide}
(u/table2
 [[gradient "first partial derivatives, returns vector"]
  [hessian "second partial derivatives, returns vector of vectors"]])

;; Let's define multivariate function x2siny

;; $$f(x,y) = x^2\sin{y}$$

;; $$\nabla f = \left(2x\sin{y}, x^2\cos{y} \right)$$

;; $$Hf = \begin{bmatrix}
;;   2\sin{y} & 2x\cos{y} \\
;;   2x\cos{y} & -x^2\sin{y}
;; \end{bmatrix}$$

^{::clerk/visibility {:result :hide}}
(defn x2siny [[x y]] (* x x (m/sin y)))
^{::clerk/visibility {:result :hide}}
(defn x2siny' [[x y]] [(* 2 x (m/sin y))
                    (* x x (m/cos y))])
^{::clerk/visibility {:result :hide}}
(defn x2siny'' [[x y]] [[(* 2 (m/sin y)) (* 2 x (m/cos y))]
                     [(* 2 x (m/cos y)) (* -1 x x (m/sin y))]])

^{::clerk/visibility {:result :hide}}
(def x2siny-fd' (calc/gradient x2siny))
^{::clerk/visibility {:result :hide}}
(def x2siny-fd'' (calc/hessian x2siny {:h 1.0e-4}))

;; Values of gradient and hessian evaluated at `[1,2]`

^{::clerk/visibility :hide ::clerk/auto-expand-results? true}
(clerk/table
 {:head ["method" "gradient" "hessian"]
  :rows [["actual" (x2siny' [1 2]) (x2siny'' [1 2])]
         ["finite difference" (x2siny-fd' [1 2]) (x2siny-fd'' [1 2])]]})



;; ## Integration

;; Univariate $R\to R$ and multivariate $R^n\to R$ functions can be integrated. 

;; Improper integrals are automatically handled and substituted in the following way. 

;; $$\int\limits_{-\infty}^a f(x)\,dx = \int\limits_{-1}^0 \frac{1}{(1+t)^2} f\left(a+\frac{t}{1+t}\right)\,dt$$
;; $$\int\limits_a^\infty f(x)\,dx = \int\limits_0^1 \frac{1}{(1-t)^2} f\left(a+\frac{t}{1-t}\right)\,dt$$
;; $$\int\limits_{-\infty}^\infty f(x)\,dx = \int\limits_{-1}^1 \frac{1+t^2}{(1-t^2)^2} f\left(\frac{t}{1-t^2}\right)\,dt$$

;; The same applies for multivariate case.
;;
;; ### Univariate

;; Integration of $f: R\to R$ functions.

^{::clerk/visibility :hide}
(u/table2
 [[integate "Integrate univariate function"]])

;; `integrate` accepts function, lower and upper bound and the following options:

;;   * `:integrator` - integration algorithm, one of: `:romberg`, `:trapezoid`, `:midpoint`, `:simpson`, `:gauss-legendre` and `:gauss-kronrod` (default).
;; * `:min-iters` - minimum number of iterations (default: 3), not used in `:gauss-kronrod`
;; * `:max-iters` - maximum number of iterations (default: 32 or 64)
;; * `:max-evals` - maximum number of evaluations, (default: maximum integer)
;; * `:rel` - relative error
;; * `:abs` - absolute error
;; * `:integration-points` - number of integration (quadrature) points for `:gauss-legendre` and `:gauss-kronrod`, default 7
;; * `:initdiv` - initial number of subdivisions for `:gauss-kronrod`, default: 1
;; * `:info?` - return full information about integration, default: false

;; `:gauss-kronrod` is a h-adaptive implementation based on [QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/)

;; #### Examples

;; First function will is the following:

;; $$f_1(x) = \frac{\ln(1-x)}{x}$$

(defn f1 ^double [^double x] (/ (m/ln (- 1.0 x)) x))

;; $$\int\limits_0^1 f_1(x)\,dx = -\frac{\pi^2}{6} \approx -1.6449341$$

(/ (m/sq m/PI) -6.0)

;; Default, Gauss-Kronrod quadrature is used, the other integrators fail.

^{::clerk/auto-expand-results? true}
(calc/integrate f1 0.0 1.0 {:max-iters 100000
                            :integration-points 64
                            :info? true
                            :abs 1.0e-8
                            :rel 0})

;; Another case is the following:

;; $$f_2(x) =  \frac{e^{-x^2} - e^{-x}}{x}$$

(defn f2 ^double [^double x] (/ (- (m/exp (- (m/sq x)))
                                (m/exp (- x))) x))

;; $$\int\limits_0^\infty f_2(x)\,dx = \frac{\gamma}{2} \approx 0.2886078$$

(/ m/GAMMA 2)

;; Most of the algorithms fail, only `:midpoint` and quadrature based algorithms work.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["method" "result" "error"]
  :rows (map (fn [method]
               (let [res (calc/integrate f2 0.0 ##Inf {:integrator method})]
                 [method res (m/abs (m/- res (m// m/GAMMA 2.0)))])) [:midpoint
                                                                     :gauss-kronrod :gauss-legendre])})

;; To overcome singularity we can slightly modify integral limits, setting next possible double value of `0.0`. 
(calc/integrate f2 (m/next-double 0.0) ##Inf {:integrator :romberg})

;; The next function is a cosine with high frequency oscillations.

;; $$f_3(x) = \cos(200x)$$

(defn f3 ^double [^double x] (m/cos (* 200.0 x)))

;; $$\int\limits_0^1 f_3(x)\,dx = \frac{\sin(200)}{200} \approx -0.0043665$$

(def f3-integral01 (/ (m/sin 200.0) 200.0))

;; As you can see some of the algorithms fail giving really wrong result.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["method" "result" "error"]
  :rows (map (fn [method]
               (let [res (calc/integrate f3 0.0 1.0 {:integrator method})]
                 [method res (m/abs (m/- res f3-integral01))])) [:trapezoid :midpoint :romberg :simpson
                                                                 :gauss-kronrod :gauss-legendre])})

;; The last example is to show the ability to use Richardson extrapolation instead of substitution for improper integrals, relying on the following equality:

;; $$\int\limits_0^\infty f(x)\,dx = \lim_{t\to\infty} \int\limits_0^t f(x)\,dx$$

;; $$f_4 = \frac{e^{-x}}{x+1}$$

(defn f4 [x] (/ (m/exp (- x)) (inc x)))

(calc/integrate f4 0 ##Inf)

;; Let's build extrapolation helper

(defn f4-extrapolated [x] (calc/integrate f4 0 x))

((calc/extrapolate (calc/fx->gx+h f4-extrapolated)) ##Inf)

;; ### Multivariate

;; There are two algorithms to calculate an integral of multivariate function. Both methods accept function (with one, vector, argument), lower and upper bounds and options. Integration can be indefinite.

^{::clerk/visibility :hide}
(u/table2
 [[cubature "h-adaptive multi-dimensional integration, based on Genz-Malik algorithm"]
  [vegas "VEGAS/VEGAS+ Monte-Carlo integration"]])

;; #### Cubature

;; `cubature` is based on algorithm described by [A.C.Genz and A.A.Malik](https://www.sciencedirect.com/science/article/pii/0771050X8090039X) and is a reimplementation of [HCubature.jl package](https://github.com/JuliaMath/HCubature.jl)
;;
;; `cubature` is an iterative and adaptative method which subdivides volume into smaller volumes to increase accuracy.

;; Options are:

;; * `:initvid` - initial subdivision per dimension, default: 2.
;; * `:max-evals` - maximum number of evaluations, default: max integer value.
;; * `:max-iters` - maximum number of iterations, default: 64.
;; * `:rel` - relative error, 1.0e-7
;; * `:abs` - absolute error, 1.0e-7
;; * `:info?` - return full information about integration, default: false

;; If `:info?` is set to `true`, more additional info is returned:

;; * `:result` - integration value
;; * `:error` - integration error
;; * `:iterations` - number of iterations
;; * `:evaluations` - number of evaluations
;; * `:subdivisions` - final number of volumes used in integration
;; * `:fail?` - set to `:max-evals` or `:max-iters` when one of the limits has been reached without desired the convergence.

;; As the first example, lets start with a 2d multi-normal distribution. If we integrate a PDF the result should be `1.0`. We need to highly raise number of iterations to achieve convergence.

(def multi-normal (partial r/pdf (r/distribution :multi-normal {:means [-1 1]})))

^{::clerk/auto-expand-results? true}
(calc/cubature multi-normal [##-Inf ##-Inf] [##Inf ##Inf] {:info? true :max-iters 10000})

;; Another example will be a volume of 3d ball, which is:

;; $$\int\limits_0^R\int\limits_0^{2\pi}\int\limits_0^\pi r^2\sin\theta\,d\theta d\varphi dr = \frac{4}{3}\pi R^3$$

(defn volume-3d-ball-integrant
  ^double [[^double theta ^double _phi ^double r]]
  (* r r (m/sin theta)))

(defn volume-3d-ball [R]
  (calc/cubature volume-3d-ball-integrant [0 0 0] [m/PI m/TWO_PI R] {:rel 0 :abs 1.0e-10
                                                                     :max-iters 1000}))

(volume-3d-ball 1)

;; Here are more result for various R compared to the real volume value, ie. $\frac{4}{3}\pi R^3$

^{::clerk/visibility :hide}
(clerk/table
 {:head ["R" "integral" "expected" "error"]
  :rows (for [R [0.1 0.5 1.0 2.0 5.0 20.0]
              :let [i (volume-3d-ball R)
                    r (* (/ 4.0 3.0) m/PI R R R)]]
          [R i r (m/abs (- i r))])})

;; #### VEGAS

;; [VEGAS and VEGAS+](https://arxiv.org/abs/2009.05112) are adaptive multidimensional (Quasi) Monte Carlo integration based on adaptive importance and stratified (VEGAS+) sampling.

;; Options are:

;; * `:max-iters` - maximum number of iterations, default: 10
;; * `:nevals` - number of evaluations per iteration, default: 10000
;; * `:nintervals` - number of grid intervals per dimension (default: 1000)
;; * `:nstrats` - number of stratifications per dimension (calculated if set to `-1`)
;; * `:warmup` - number of warmup iterations (results are used to train stratification and grid spacings, default: `0`
;; * `:alpha` - grid refinement parameter, 0.5 slow (default for vegas+), 1.5 moderate/fast (defatult for vegas)
;; * `:beta` - stratification damping parameter for startification adaptation, default: 0.75
;; * `:rel` - relative accuracy, default: 5.0e-4
;; * `:abs` - absolute accuracy, default: 5.0e-4
;; * `:random-sequence` - random sequence used for generating samples: `:uniform` (default), low-discrepancy sequences: `:r2`, `:sobol` and `:halton`.
;; * `:jitter` - jittering factor for low-discrepancy random sequence, default: 0.75
;; * `:info?` - return full information about integration, default: false
;; * `:record-data?` - stores samples, number of strata, x and dx for each iteration, default: false (requires, `:info?` to be set to `true`)

;; Function returns a map with following keys (if info? is true, returns result otherwise):

;; * `:result` - value of integral
;; * `:iterations` - number of iterations (excluding warmup)
;; * `:sd` - standard deviation of results
;; * `:nintervals` - actual grid size
;; * `:nstrats` - number of stratitfications per dimension
;; * `:nhcubes` - number of hypercubes
;; * `:evaluations` - number of function calls
;; * `:chi2-avg` - average of chi2
;; * `:dof` - degrees of freedom
;; * `:Q` - goodness of fit indicator, 1 - very good, <0.25 very poor
;; * `:data` - recorded data (if available)

;; Importance and stratified sampling allows algorithm to sample function in regions where the function vary the most. In the following example we'll see how this adaptation is performed.

;; Original VEGAS used only importance sampling. The options which are used to adjust this part of the algorithm are: `nintervals` and `alpha`. The former sets the size of the grid which is later adjusted to condense sampling in the areas of the highest variation. The latter is responsible for speed of the adjustment (`0.5` is quite low, while `1.5` is normal). The grid subdivision can be high, by default it's `1000`.

;; VEGAS+ add stratified sampling on the top of VEGAS. It splits volume to the number of hypercuebes or hyperrectangles. Each subvolume gets the number of points proportional to the variance of the function values. `nstrats` and `beta` are options used to control adaptation of the algorithm. `nstrats` can be a number or a list of number of strats per dimension (`-1` for automatic selection). Set to `0` or `1` to turn off stratification. `beta` controls pace of the adaptation.

;; Let's illustrate it in the following example. Let's define a function which creates nice shape as a graph (white - near `0`, dark - the highest value). The function has its highest values near `(1.5,0)` and is around `2.16`:

;; $$f(x,y)=e^{-2(\cos(y)-\frac{x^2}{2})^2-\cos(x+1)}$$

(defn func ^double [[^double x ^double y]]
  (m/exp (m/- (m/* -2.0 (m/sq (m/- (m/cos y) (m// (m/sq x) 2))))
              (m/cos (m/+ 1.0 x)))))

^{::clerk/visibility :hide ::clerk/no-cache true}
(u/graph2d (fn [x y] (func [x y])) [m/-PI m/PI] [m/-PI m/PI])

;; The value of the integral over $[-\pi,\pi]\times[-\pi,\pi]$ region is:

;; $$\int\limits_{-\pi}^\pi \int\limits_{-\pi}^\pi e^{-2(\cos(y)-\frac{x^2}{2})^2-\cos(x+1)} \,dx\,dy \approx 8.52775$$

;; Cubature algorithm resturns the expected value:

^{::clerk/auto-expand-results? true}
(calc/cubature func [m/-PI m/-PI] [m/PI m/PI] {:info? true :max-iters 10000})

;; Below `vegas` is set to the original algorithm (without stratification) with number of points is set to `5000`, small number of importance sampling intervals. We will use Quasi Monte Carlo by using jittered Halton low-discrepancy sequence to generate samples.

^{::clerk/auto-expand-results? true}
(calc/vegas func [m/-PI m/-PI] [m/PI m/PI] {:info? true :nstrats 0 :nevals 5000
                                            :nintervals 20 :max-iters 51
                                            :random-sequence :halton :jitter 0.5})


;; See how adaptation is performed at selected iterations. The 20x20 grid is squeezed and stretched to adopt to the function shape.

^{::clerk/visibility :hide ::clerk/no-cache true}
(clerk/row
 (clerk/caption "initial" (clerk/image "notebooks/images/vegas/0.jpg"))
 (clerk/caption "1st iteration" (clerk/image "notebooks/images/vegas/1.jpg"))
 (clerk/caption "5th" (clerk/image "notebooks/images/vegas/5.jpg"))
 (clerk/caption "10th" (clerk/image "notebooks/images/vegas/10.jpg"))
 (clerk/caption "50th" (clerk/image "notebooks/images/vegas/50.jpg")))

;; When we add stratification the coverage is even better. Number of stratification is `10*10=100`. `beta` is set to a little bit lower value than default (`0.5` vs `0.75`). The rest is as above.

^{::clerk/auto-expand-results? true}
(calc/vegas func [m/-PI m/-PI] [m/PI m/PI] {:info? true :nstrats 10
                                            :beta 0.5 :nevals 5000
                                            :nintervals 20 :max-iters 51
                                            :random-sequence :halton :jitter 0.5})

;; You can observe that flat surfaces are less sampled than the areas where function varies most.

^{::clerk/visibility :hide ::clerk/no-cache true}
(clerk/row
 (clerk/caption "initial" (clerk/image "notebooks/images/vegasp/0.jpg"))
 (clerk/caption "1st iteration" (clerk/image "notebooks/images/vegasp/1.jpg"))
 (clerk/caption "5th" (clerk/image "notebooks/images/vegasp/5.jpg"))
 (clerk/caption "10th" (clerk/image "notebooks/images/vegasp/10.jpg"))
 (clerk/caption "50th" (clerk/image "notebooks/images/vegasp/50.jpg")))

;; Now let's see the result of examples explored in `cubature` section: multi-normal and 3d ball volume.

^{::clerk/auto-expand-results? true}
(repeatedly 5 #(calc/vegas multi-normal
                           [##-Inf ##-Inf]
                           [##Inf ##Inf]
                           {:max-iters 20 :warmup 5
                            :nstrats 5
                            :abs 1.0e-6 :rel 1.0e-6
                            :max-evals 50000}))

(defn volume-3d-ball-vegas [R]
  (calc/vegas volume-3d-ball-integrant [0 0 0] [m/PI m/TWO_PI R] {:rel 1.0e-6 :abs 1.0e-6
                                                                  :warmup 5
                                                                  :max-iters 20}))

(volume-3d-ball-vegas 1)

^{::clerk/visibility :hide}
(clerk/table
 {:head ["R" "integral" "expected" "error"]
  :rows (for [R [0.1 0.5 1.0 2.0 5.0 20.0]
              :let [i (volume-3d-ball-vegas R)
                    r (* (/ 4.0 3.0) m/PI R R R)]]
          [R i r (m/abs (- i r))])})

;; The last example show the ability to calculate integral when the function is really high dimensional. Let's consider the function:

;; $$f(\vec{v}) = e^{-\|\vec{v}\|_2^2}$$

(defn exp-magsq-v
  ^double [v] (m/exp (m/- (v/magsq v))))

;; And the intergal over volume `V`.

;; $$\mathop{\int\dots\int}\limits_V e^{-\|\vec{v}\|_2^2}\,dV = \left[\frac{\sqrt\pi}{2}\left(erf(1)+erf(3)\right)\right]^n$$
;; $$\textrm{ where } V=[-1,3]\times\dots\times[-1,3]$$

(defn exp-magsq-v-result
  [n] (m/fpow (* m/SQRTPI 0.5 (+ (m/erf 1) (m/erf 3))) n))

;; Now we compare two results, one from `cubature` and second from `vegas` in 2d.

(exp-magsq-v-result 2)

^{::clerk/auto-expand-results? true}
(def cubature-result (calc/cubature exp-magsq-v [-1 -1] [3 3] {:max-iters 1000
                                                             :abs 1.0e-8
                                                             :rel 1.0e-8
                                                             :info? true}))

^{::clerk/auto-expand-results? true ::clerk/no-cache true}
(def vegas-result (calc/vegas exp-magsq-v [-1 -1] [3 3] {:warmup 5
                                                       :max-iters 50 :nevals 100000
                                                       :random-sequence :r2
                                                       :info? true}))

^{::clerk/auto-expand-results? true ::clerk/no-cache true}
{:vegas-error (m/abs (- (exp-magsq-v-result 2) (:result vegas-result)))
 :cubature-error (m/abs (- (exp-magsq-v-result 2) (:result cubature-result)))}

;; And now in 16-dimensional volume:

(exp-magsq-v-result 16)

;; Let's limit hypercubes in stratification and subdivide every second dimension. Also we lower adaptation ratio for stratification.

^{::clerk/auto-expand-results? true}
(calc/vegas exp-magsq-v (repeat 16 -1) (repeat 16 3) {:nstrats [4 1 2 1]
                                                      :beta 0.2
                                                      :warmup 1 :nevals 150000
                                                      :max-iters 10 :info? true})



;; ## Solvers

;; Solvers are used to find zeros of given $R\to R$ function. 

#_:clj-kondo/ignore
(require '[fastmath.solver :as solver])

^{::clerk/visibility :hide}
(u/table2
 [[find-root "Find zero of a function in given interval"]
  [quadratic "Solve quadratic equation"]])

;; ### Solver

;; `find-root` function accepts a function and bounds. Additional options are:

;; * `:absolute-accuracy` - absolute tolerance
;; * `:relative-accuracy` - relative tolerance
;; * `:max-iters` - maximum number of iterations
;; * `initial-value` - starting point
;; * `solver` - algorithm used to find root of the function, one of: `:brent` (default), `:bisection`, `:illinois`, `:muller`, `:muller2`, `:pegasus`, `:regula-falsi`, `:ridders` and `:secant`

;; For given function:

;; $$f(x) = \cos{x}\left( x-1 \right)^3$$

(defn solver-target ^double [^double x]
  (* (m/cos x) (m/fpow (dec x) 3)))

;; we have a root for $x=1$ on a `[-1.1, 1.1]` range. As you can see, `:regula-falsi` fails to find a root for given absolute accuracy (`1.0e-10` in this case)

^{::clerk/visibility :hide}
(clerk/table
 {:head ["solver" "result"]
  :rows (for [solver [:brent :bisection :illinois :muller :muller2 :pegasus :regula-falsi :ridders :secant]]
          (do (println solver)
              [solver (try (solver/find-root solver-target -1.1 1.1 {:solver solver :absolute-accuracy 1.0e-10})
                           (catch Exception e (.getMessage e)))]))})

;; `:regula-falsi` solver doesn't work well for above example but gives the result for below polynomial.

(solver/find-root (fn [x] (- (* x x) 0.1)) 0 1 {:solver :regula-falsi})

;; ### Quadratic

;; This is a special case which analytically find roots of any quadratic equations with given coefficients

;; $$q(x) = ax^2+bx+c$$

;; Returns a pair of doubles or `nil` when no solution is found

^{::clerk/visibility :hide}
(clerk/example
  (solver/quadratic 1 0 -0.1)
  (solver/quadratic 1 -2 1)
  (solver/quadratic 1 1 1))

;; # List of symbols

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.calculus {:exclude-vars [#"^\->.*"]})

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.solver)


#_(do
    (require '[clojure2d.core]
             '[clojure.java.shell])

    (clojure.java.shell/sh "rm" "-fr" "notebooks/images/vegas/" "notebooks/images/vegasp/")

    (let [res1 (:data (calc/vegas func [m/-PI m/-PI] [m/PI m/PI] {:info? true :nstrats 0 :warmup 0 :nevals 5000 :nintervals 20 :record-data? true :max-iters 51 :random-sequence :halton :jitter 0.5}))
          res2 (:data (calc/vegas func [m/-PI m/-PI] [m/PI m/PI] {:info? true :beta 0.5 :nstrats 10 :warmup 0 :nevals 5000 :nintervals 20 :record-data? true :max-iters 51 :random-sequence :halton :jitter 0.5}))]
      (doseq [pos [0 1 5 10 50]]
        (clojure2d.core/save (u/scatter (:samples (nth res1 pos)) [m/-PI m/PI] [m/-PI m/PI])
                             (str "notebooks/images/vegas/" pos ".jpg"))
        (clojure2d.core/save (u/scatter (:samples (nth res2 pos)) [m/-PI m/PI] [m/-PI m/PI])
                             (str "notebooks/images/vegasp/" pos ".jpg")))
      ))
