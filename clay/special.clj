^:kindly/hide-code
(ns special
  (:require [fastmath.core :as m]
            [fastmath.special :as special]
            [ggplot :as gg]
            [codox]))

;; # Special functions {.unnumbered}

;; Collection of special functions for real arguments and real returned value. Most of the functions are implemented natively in Clojure, some are based on Apache Commons Math.

;; Native implementation is based on Julia packages ([SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions), [Bessel](https://github.com/JuliaMath/Bessels.jl), [HypergeometricFunctions](https://github.com/JuliaMath/HypergeometricFunctions.jl)) or scientific papers and books ([NIST](https://dlmf.nist.gov/), Meshfree Approximation Methods with Matlab by G. E. Fasshauer]).

(require '[fastmath.special :as special]
         '[fastmath.core :as m])

;; ## Gamma

;; ### Gamma function

;; `gamma` $\Gamma(x)$ function is an extension of the factorial.

;; $$\Gamma(x) = \int_0^\infty t^{x-1}e^{-t}\,dt$$

;; For positive integer `n`

;; $$\Gamma(n) = (n-1)!$$

(gg/->image (gg/function special/gamma {:x [-4.5 4.5] :ylim [-5 5] :steps 1000 :title "Gamma function"}))

(special/gamma 1.5)
(special/gamma -1.5)

;; Gamma for negative integers is not defined.

(special/gamma -2.0)
(special/gamma 15)
(m/factorial 14)

;; ---

;; Additionally reciprocal gamma function `inv-gamma-1pm1` is defined as:

;; $$\frac{1}{\Gamma(x+1)}-1\text{ for } -0.5\leq x\leq 1.5$$

(gg/->image (gg/function special/inv-gamma-1pm1 {:x [-0.5 1.5] :title "Reciprocal of gamma"}))

(special/inv-gamma-1pm1 -0.5)
(special/inv-gamma-1pm1 0.5)

;; ### Log gamma 

;; Logartihm of gamma `log-gamma` $\log\Gamma(x)$ with derivatives: `digamma` $\psi$, `trigamma` $\psi_1$  and `polygamma` $\psi^{(m)}$.

(gg/->image (gg/function special/log-gamma {:x [0 5] :title "Log-gamma (log-gamma)"}))

(special/log-gamma 1.01)

;; ---

;; `log-gamma-1p` is more accurate function defined as $\log\Gamma(x+1)$ for $-0.5\leq x 1.5$

(gg/->image (gg/function special/log-gamma-1p {:x [-0.5 1.5] :title "Log-gamma (log-gamma-1p)"}))

(special/log-gamma-1p 0.01)

;; ---

;; Derivatives of `log-gamma` function. First derivative `digamma`.

;; $$\operatorname{digamma}(x)=\psi(x)=\psi^{(0)}(x)=\frac{d}{dx}\log\Gamma(x)=\frac{\Gamma'(x)}{\Gamma(x)}$$

(gg/->image (gg/function special/digamma {:x [-3 3] :ylim [-10 10] :title "digamma"}))

(special/digamma 0.5)

;; ---

;; Second derivative `trigamma`

;; $$\operatorname{trigamma}(x)=\psi_1(x)=\psi^{(1)}(x)=\frac{d}{dx}\psi(x)=\frac{d^2}{dx^2}\log\Gamma(x)$$

(gg/->image (gg/function special/trigamma {:x [-3 3] :ylim [0 10] :title "trigamma"}))

(special/trigamma 0.5)

;; ---

;; `polygamma` as m^th^ derivative of `digamma`

;; $$\operatorname{polygamma}(m,x)=\psi^{(m)}=\frac{d^m}{dx^m}\psi(x)=\frac{d^{m+1}}{dx^{m+1}}\log\Gamma(x)$$

(gg/->image (gg/functions [["0" (partial special/polygamma 0)]
                           ["1" (partial special/polygamma 1)]
                           ["2" (partial special/polygamma 2)]
                           ["3" (partial special/polygamma 3)]]
                          {:x [-2 8]
                           :ylim [-10 13]
                           :legend-name "m"
                           :title "polygamma(m,x)"
                           :palette gg/palette-blue-1}))

(special/trigamma 0.5)

;; ### Incomplete and regularized

;; `upper-incomplete-gamma` $\Gamma(s,x)$ and `lower-incomplete-gamma` $\gamma(s,x)$ are versions of gamma function with parametrized integral limits.

;; Upper incomplete gamma is defined as:

;; $$\Gamma(s,x) = \int_x^\infty t^{s-1}e^{-t}\,dt$$

^:kind/table
[[(gg/->image (gg/functions [["1" (partial special/upper-incomplete-gamma 1)]
                             ["2" (partial special/upper-incomplete-gamma 2)]
                             ["2.5" (partial special/upper-incomplete-gamma 2.5)]
                             ["3" (partial special/upper-incomplete-gamma 3)]]
                            {:x [0 8]
                             :ylim [0 2]
                             :legend-name "s"
                             :title "upper incomplete gamma, positive s"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-1" (partial special/upper-incomplete-gamma -1)]
                             ["-2" (partial special/upper-incomplete-gamma -2)]
                             ["-2.5" (partial special/upper-incomplete-gamma -2.5)]
                             ["-3" (partial special/upper-incomplete-gamma -3)]]
                            {:x [0 5]
                             :ylim [0 1]
                             :legend-name "s"
                             :title "upper incomplete gamma, negative s"
                             :palette gg/palette-blue-1}))]]

(special/upper-incomplete-gamma 2 0.5)
(special/upper-incomplete-gamma -2 0.5)

;; ---

;; Lower incomplete gamma is defined as:

;; $$\gamma(s,x) = \int_0^x t^{s-1}e^{-t}\,dt$$

^:kind/table
[[(gg/->image (gg/functions [["1" (partial special/lower-incomplete-gamma 1)]
                             ["2" (partial special/lower-incomplete-gamma 2)]
                             ["2.5" (partial special/lower-incomplete-gamma 2.5)]
                             ["3" (partial special/lower-incomplete-gamma 3)]]
                            {:x [0 8]
                             :legend-name "s"
                             :title "Lower incomplete gamma, positive s"
                             :palette gg/palette-blue-1}))

  (gg/->image (gg/functions [["-0.5" (partial special/lower-incomplete-gamma -0.5)]
                             ["-1.5" (partial special/lower-incomplete-gamma -1.5)]
                             ["-2.5" (partial special/lower-incomplete-gamma -2.5)]]
                            {:x [0 8]
                             :ylim [-10 4]
                             :steps 1000
                             :legend-name "s"
                             :title "Lower incomplete gamma, negative s"
                             :palette gg/palette-blue-1}))]]

(special/lower-incomplete-gamma 0.5 3)
(special/lower-incomplete-gamma -0.5 3)

;; ---

;; `regularized-gamma-p` $P(s,x)$ and `regularized-gamma-q` $Q(s,x)$ are normalized incomplete gamma functions by gamma of `s`. `s` can be negative non-integer.

;; $$P(s,x)=\frac{\gamma(s,x)}{\Gamma(x)}$$

(gg/->image (gg/functions [["1" (partial special/regularized-gamma-p 1)]
                           ["2" (partial special/regularized-gamma-p 2)]
                           ["2.5" (partial special/regularized-gamma-p 2.5)]
                           ["3" (partial special/regularized-gamma-p 3)]]
                          {:x [0 8]
                           :legend-name "s"
                           :title "Regularized gamma P(s,x)"
                           :palette gg/palette-blue-1}))

(special/regularized-gamma-p 2.5 0.5)
(special/regularized-gamma-p -2.5 0.5)

;; ---

;; $$Q(s,x)=\frac{\Gamma(s,x)}{\Gamma(x)}=1-P(s,x)$$

(gg/->image (gg/functions [["1" (partial special/regularized-gamma-q 1)]
                           ["2" (partial special/regularized-gamma-q 2)]
                           ["2.5" (partial special/regularized-gamma-q 2.5)]
                           ["3" (partial special/regularized-gamma-q 3)]]
                          {:x [0 8]
                           :legend-name "s"
                           :title "Regularized gamma Q(s,x)"
                           :palette gg/palette-blue-1}))

(special/regularized-gamma-q 2.5 0.5)
(special/regularized-gamma-q -2.5 0.5)

;; ## Beta

;; `beta` $B(p,q)$ function, defined also for negative non-integer `p` and `q`.

;; $$B(p,q) = \int_0^1 t^{p-1}(1-t)^{q-1}\,dt = \frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}$$

(gg/->image (gg/functions [["-0.5" (partial special/beta -0.5)]
                           ["0.5" (partial special/beta 0.5)]
                           ["1.5" (partial special/beta 1.5)]]
                          {:x [0.01 2]
                           :ylim [-5 50]
                           :legend-name "p"
                           :title "Beta(p,x)"
                           :palette gg/palette-blue-1}))

(special/beta 2 3)
(special/beta -1.2 0.1)
(special/beta -1.2 -0.1)

;; ---

;; `log-beta` is log of $B(p,q)$ for positive `p` and `q`

(gg/->image (gg/functions [["0.5" (partial special/log-beta 0.5)]
                           ["1.5" (partial special/log-beta 1.5)]
                           ["2.5" (partial special/log-beta 2.5)]]
                          {:x [0.01 2]
                           :legend-name "p"
                           :title "log of Beta(p,x)"
                           :palette gg/palette-blue-1}))

[(special/log-beta 2 3)
 (m/log (special/beta 2 3))]

;; ---

;; `incomplete-beta` $B(x,a,b)$ and `regularized-beta` $I_x(a,b)$. Both are defined also for negative non-integer `a` and `b`.

;; $$B(x,a,b)=\int_0^x t^{a-1}(1-t)^{b-1}\,dt$$

(gg/->image (gg/functions [["1,1" #(special/incomplete-beta % 1 1)]
                           ["1,2" #(special/incomplete-beta % 1 2)]
                           ["1.5,0.5" #(special/incomplete-beta % 1.5 0.5)]]
                          {:x [0.01 0.99]
                           :legend-name "a, b"
                           :title "Beta(x,a,b)"
                           :palette gg/palette-blue-1}))

(special/incomplete-beta 0.5 0.1 0.2)
(special/incomplete-beta 0.5 -0.1 -0.2)

;; $$I_x(a,b)=\frac{B(x,a,b)}{B(a,b)}$$

(gg/->image (gg/functions [["1,1" #(special/regularized-beta % 1 1)]
                           ["1,2" #(special/regularized-beta % 1 2)]
                           ["1.5,0.5" #(special/regularized-beta % 1.5 0.5)]]
                          {:x [0.01 0.99]
                           :legend-name "a, b"
                           :title "Beta(x,a,b)"
                           :palette gg/palette-blue-1}))


(special/regularized-beta 0.99 -0.5 -0.7)
(special/regularized-beta 0.01 0.5 0.5)

;; ## Bessel

;; ### J

;; Bessel functions of the first kind, `bessel-J` $J_\alpha(x)$. `bessel-J0` and `bessel-J1` are bessel of orders `0` and `1`. An order should be integer for negative arguments.

(gg/->image (gg/functions [["0" special/bessel-J0]
                           ["1" special/bessel-J1]
                           ["2" (partial special/bessel-J 2)]
                           ["3" (partial special/bessel-J 3)]]
                          {:x [-5 20]
                           :legend-name "order"
                           :title "J(a,x)"
                           :palette gg/palette-blue-1}))

(special/bessel-J0 2.3)
(special/bessel-J1 2.3)
(special/bessel-J 2.1 3)

(special/bessel-J -3 -3.2)
(special/bessel-J 3 -3.2)

(special/bessel-J -3.1 -3.2)
(special/bessel-J 3.1 -3.2)

;; ---

;; Additional `jinc` (sombrero) function is defined as:

;; $$\operatorname{jinc}(x)=\frac{2J_1(\pi x)}{\pi x}$$

(gg/->image (gg/function special/jinc {:x [-5 5] :title "jinc(x)"}))

(special/jinc 2.3)

;; ### Y

;; Bessel functions of the second kind, `bessel-Y`, $Y_\alpha(x)$. `bessel-Y0` and `bessel-Y1` are functions of orders `0` and `1`. They are defined for positive argument only and any order.

(gg/->image (gg/functions [["0" special/bessel-Y0]
                           ["1" special/bessel-Y1]
                           ["2" (partial special/bessel-Y 2)]
                           ["3" (partial special/bessel-Y 3)]]
                          {:x [0.01 20]
                           :ylim [-1 nil]
                           :legend-name "order"
                           :title "Y(a,x)"
                           :palette gg/palette-blue-1}))

(special/bessel-Y0 2.3)
(special/bessel-Y1 2.3)
(special/bessel-Y 2 2.3)
(special/bessel-Y -2.1 2.3)
(special/bessel-Y 3 -1)

;; ### I

;; Modified Bessel functions of the first kind, `bessel-I`, $I_\alpha(x)$, `bessel-I0` and `bessel-I1` are functions of orders `0` and `1`. An order should be integer for negative arguments.

(gg/->image (gg/functions [["0" special/bessel-I0]
                           ["1" special/bessel-I1]
                           ["2" (partial special/bessel-I 2)]
                           ["3" (partial special/bessel-I 3)]]
                          {:x [-5 5]
                           :legend-name "order"
                           :title "I(a,x)"
                           :palette gg/palette-blue-1}))


(special/bessel-I0 2.3)
(special/bessel-I1 2.3)
(special/bessel-I 2 2.3)
(special/bessel-I -2.1 2.3)
(special/bessel-I -3 -1)
(special/bessel-I -3.1 -1)

;; ### K

;; Modified Bessel functions of the second kind, `bessel-K`, $K_\alpha(x)$, `bessel-K0` and `bessel-K1` are functions of orders `0` and `1`. They are defined for positive argument only and any order.

(gg/->image (gg/functions [["0" special/bessel-K0]
                           ["1" special/bessel-K1]
                           ["2" (partial special/bessel-K 2)]
                           ["3" (partial special/bessel-K 3)]]
                          {:x [0.01 5]
                           :ylim [nil 5]
                           :legend-name "order"
                           :title "K(a,x)"
                           :palette gg/palette-blue-1}))

(special/bessel-K0 2.3)
(special/bessel-K1 2.3)
(special/bessel-K 2 2.3)
(special/bessel-K -2.1 2.3)
(special/bessel-K 3 -1)

;; ---

;; Additionally `bessel-K-half-odd` function is optimized version for order of the half of odd integer, ie `1/2`, `3/2`, `5/2` and so on. First argument is an odd numerator. 

(gg/->image (gg/functions [["1/2" (partial special/bessel-K-half-odd 1)]
                           ["3/2" (partial special/bessel-K-half-odd 3)]
                           ["5/2" (partial special/bessel-K-half-odd 5)]
                           ["7/2" (partial special/bessel-K-half-odd 7)]]
                          {:x [0.01 5]
                           :ylim [nil 5]
                           :legend-name "order"
                           :title "K(a/2,x)"
                           :palette gg/palette-blue-1}))


(special/bessel-K-half-odd 1 2.3)
(special/bessel-K-half-odd 3 2.3)
(special/bessel-K-half-odd 5 2.3)

;; ## Erf

;; ## Airy

;; ## Integrals

;; ## Hypergeometric

;; ## Zeta

;; ## Other

;; ### Lambert W

;; ### Harmonic H

;; ### Minowski

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.special)
