^:kindly/hide-code
(ns special
  (:require [fastmath.core :as m]
            [fastmath.special :as special]
            [fastmath.complex :as complex]
            
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]))


;; # Special functions {.unnumbered}

;; Collection of special functions for real arguments and real returned value. Most of the functions are implemented natively in Clojure, some are based on Apache Commons Math.

;; Native implementation is based on Julia packages ([SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions), [Bessel](https://github.com/JuliaMath/Bessels.jl), [HypergeometricFunctions](https://github.com/JuliaMath/HypergeometricFunctions.jl)) or scientific papers and books ([NIST](https://dlmf.nist.gov/), Meshfree Approximation Methods with Matlab by G. E. Fasshauer]).

(require '[fastmath.special :as special]
         '[fastmath.core :as m]
         '[fastmath.complex :as complex])

;; ## Gamma

;; Gamma and related functions

;; ::: {.callout-tip title="Defined functions"}
;; * `gamma`, `log-gamma`
;; * `inv-gamma-1pm1`, `log-gamma-1p`
;; * `upper-incomplete-gamma`, `lower-incomplete-gamma`
;; * `regularized-gamma-p`, `regularized-gamma-q`
;; * `digamma`, `trigamma`, `polygamma`
;; * `gamma-complex`, `log-gamma-complex`
;; :::

;; ### Gamma function

;; `gamma` $\Gamma(x)$ function is an extension of the factorial.

;; $$\Gamma(x) = \int_0^\infty t^{x-1}e^{-t}\,dt$$

;; For positive integer `n`

;; $$\Gamma(n) = (n-1)!$$

;; Gamma for negative integers is not defined.

(gg/->image (gg/function special/gamma {:x [-4.5 4.5] :ylim [-5 5] :steps 1000 :title "Gamma function"}))

(utls/examples-note
  (special/gamma 1.5)
  (special/gamma -1.5)
  (special/gamma -2.0)
  (special/gamma 15)
  (m/factorial 14))

;; ---

;; Additionally reciprocal gamma function `inv-gamma-1pm1` is defined as:

;; $$\frac{1}{\Gamma(x+1)}-1\text{ for } -0.5\leq x\leq 1.5$$

(gg/->image (gg/function special/inv-gamma-1pm1 {:x [-0.5 1.5] :title "Reciprocal of gamma"}))

(utls/examples-note
  (special/inv-gamma-1pm1 -0.5)
  (special/inv-gamma-1pm1 0.5))

;; ---

;; For complex argument call `gamma-complex`.

(gg/->image (gg/complex-function special/gamma-complex {:x [-5.1 5.1] :y [-5.1 5.1]
                                                        :steps 500 :title "Complex gamma"}))

(utls/examples-note
  (special/gamma-complex (complex/complex 1.5 0.0))
  (special/gamma-complex (complex/complex 0.0 1.0))
  (special/gamma-complex (complex/complex -1.5 -1.0)))

;; ### Log gamma 

;; Logartihm of gamma `log-gamma` $\log\Gamma(x)$ with derivatives: `digamma` $\psi$, `trigamma` $\psi_1$  and `polygamma` $\psi^{(m)}$.

(gg/->image (gg/function special/log-gamma {:x [0 5] :title "Log-gamma (log-gamma)"}))

(utls/examples-note
  (special/log-gamma 1.01)
  (special/log-gamma 0.5)
  (m/exp (special/log-gamma 5)))

;; ---

;; For complex argument call `log-gamma-complex`.

(gg/->image (gg/complex-function special/log-gamma-complex {:x [-5.1 5.1] :y [-5.1 5.1]
                                                            :steps 500 :title "Complex log-gamma"}))

(utls/examples-note
  (special/log-gamma-complex (complex/complex 1.01 0.0))
  (special/log-gamma-complex (complex/complex 0.0 1.0))
  (special/log-gamma-complex (complex/complex -1.5 -1.0)))

;; ---

;; `log-gamma-1p` is more accurate function defined as $\log\Gamma(x+1)$ for $-0.5\leq x 1.5$

(gg/->image (gg/function special/log-gamma-1p {:x [-0.5 1.5] :title "Log-gamma (log-gamma-1p)"}))

(utls/examples-note
  (special/log-gamma-1p -0.1)
  (special/log-gamma-1p 0.01)
  (special/log-gamma-1p 1.01))

;; ---

;; Derivatives of `log-gamma` function. First derivative `digamma`.

;; $$\operatorname{digamma}(x)=\psi(x)=\psi^{(0)}(x)=\frac{d}{dx}\log\Gamma(x)=\frac{\Gamma'(x)}{\Gamma(x)}$$

(gg/->image (gg/function special/digamma {:x [-3 3] :ylim [-10 10] :title "digamma"}))

;; ---

;; Second derivative `trigamma`

;; $$\operatorname{trigamma}(x)=\psi_1(x)=\psi^{(1)}(x)=\frac{d}{dx}\psi(x)=\frac{d^2}{dx^2}\log\Gamma(x)$$

(gg/->image (gg/function special/trigamma {:x [-3 3] :ylim [0 10] :title "trigamma"}))

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

(utls/examples-note
  (special/digamma 0.5)
  (special/trigamma 0.5)
  (special/polygamma 0 0.5)
  (special/polygamma 1 0.5)
  (special/polygamma 2 0.5)
  (special/polygamma 3 0.5))

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

(utls/examples-note
  (special/upper-incomplete-gamma 2 0.5)
  (special/upper-incomplete-gamma -2 0.5))

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

(utls/examples-note
  (special/lower-incomplete-gamma 0.5 3)
  (special/lower-incomplete-gamma -0.5 3))

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

(utls/examples-note
  (special/regularized-gamma-p 2.5 0.5)
  (special/regularized-gamma-p -2.5 0.5))

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

(utls/examples-note
  (special/regularized-gamma-q 2.5 0.5)
  (special/regularized-gamma-q -2.5 0.5))

;; ## Beta

;; Beta and related functions

;; ::: {.callout-tip title="Defined functions"}
;; * `beta`, `log-beta`
;; * `incomplete-beta`, `regularized-beta`
;; :::

;; ### Beta function

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

(utls/examples-note
  (special/beta 2 3)
  (special/beta -1.2 0.1)
  (special/beta -1.2 -0.1))

;; ---

;; `log-beta` is log of $B(p,q)$ for positive `p` and `q`

(gg/->image (gg/functions [["0.5" (partial special/log-beta 0.5)]
                           ["1.5" (partial special/log-beta 1.5)]
                           ["2.5" (partial special/log-beta 2.5)]]
                          {:x [0.01 2]
                           :legend-name "p"
                           :title "log of Beta(p,x)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/log-beta 2 3)
  (special/log-beta 2 3)
  (m/log (special/beta 2 3)))

;; ### Incomplete and regularized

;; `incomplete-beta` $B(x,a,b)$ and `regularized-beta` $I_x(a,b)$. Both are defined also for negative non-integer `a` and `b`.

;; $$B(x,a,b)=\int_0^x t^{a-1}(1-t)^{b-1}\,dt$$

(gg/->image (gg/functions [["1,1" #(special/incomplete-beta % 1 1)]
                           ["1,2" #(special/incomplete-beta % 1 2)]
                           ["1.5,0.5" #(special/incomplete-beta % 1.5 0.5)]]
                          {:x [0.01 0.99]
                           :legend-name "a, b"
                           :title "Beta(x,a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/incomplete-beta 0.5 0.1 0.2)
  (special/incomplete-beta 0.5 -0.1 -0.2))

;; $$I_x(a,b)=\frac{B(x,a,b)}{B(a,b)}$$

(gg/->image (gg/functions [["1,1" #(special/regularized-beta % 1 1)]
                           ["1,2" #(special/regularized-beta % 1 2)]
                           ["1.5,0.5" #(special/regularized-beta % 1.5 0.5)]]
                          {:x [0.01 0.99]
                           :legend-name "a, b"
                           :title "Beta(x,a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/regularized-beta 0.99 -0.5 -0.7)
  (special/regularized-beta 0.01 0.5 0.5))

;; ## Bessel

;; * Bessel functions of the first ($J_\alpha$) and the second ($Y_\alpha$) kind
;; * Modified Bessel functions of the first ($I_\alpha$) and the second ($K_\alpha$) kind
;; * Spherical Bessel functions of the first ($j_\alpha$) and the second ($y_\alpha$) kind
;; * Modified spherical Bessel functions of the first ($i_\alpha^{(1)}$, $i_\alpha^{(2)}$) and the second ($k_\alpha$) kind
;; * Sombrero function `jinc`

;; ::: {.callout-tip title="Defined functions"}
;; * `bessel-J0`, `bessel-J1`, `bessel-J`, `jinc`
;; * `bessel-Y0`, `bessel-Y1`, `bessel-Y`
;; * `bessel-I0`, `bessel-I1`, `bessel-I`
;; * `bessel-K0`, `bessel-K1`, `bessel-K`, `bessel-K-half-odd`
;; * `spherical-bessel-j0`, `spherical-bessel-j1`, `spherical-bessel-j2`,`spherical-bessel-j`
;; * `spherical-bessel-y0`, `spherical-bessel-y1`, `spherical-bessel-y2`,`spherical-bessel-y`
;; * `spherical-bessel-1-i0`, `spherical-bessel-1-i1`, `spherical-bessel-1-i2`,`spherical-bessel-1-i`
;; * `spherical-bessel-2-i0`, `spherical-bessel-2-i1`, `spherical-bessel-2-i2`,`spherical-bessel-2-1`
;; * `spherical-bessel-k0`, `spherical-bessel-k1`, `spherical-bessel-k2`,`spherical-bessel-k`
;; :::

;; ### Bessel J, j

;; Bessel functions of the first kind, `bessel-J` $J_\alpha(x)$. `bessel-J0` and `bessel-J1` are functions of orders `0` and `1`. An order should be integer for negative arguments.

^:kind/table
[[(gg/->image (gg/functions [["0" special/bessel-J0]
                             ["1" special/bessel-J1]
                             ["1.5" (partial special/bessel-J 1.5)]
                             ["2" (partial special/bessel-J 2)]
                             ["3" (partial special/bessel-J 3)]]
                            {:x [-5 20]
                             :legend-name "order"
                             :title "J(a,x)"
                             :palette gg/palette-blue-1}))

  (gg/->image (gg/functions [["-0.5" (partial special/bessel-J -0.5)]
                             ["-1" (partial special/bessel-J -1)]
                             ["-1.5" (partial special/bessel-J -1.5)]
                             ["-2" (partial special/bessel-J -2)]
                             ["-2.5" (partial special/bessel-J -2.5)]]
                            {:x [-5 20]
                             :ylim [-1 1]
                             :legend-name "order"
                             :title "J(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/bessel-J0 2.3)
  (special/bessel-J1 2.3)
  (special/bessel-J 2.1 3)

  (special/bessel-J -3 -3.2)
  (special/bessel-J 3 -3.2)

  (special/bessel-J -3.1 -3.2)
  (special/bessel-J 3.1 -3.2))

;; ---

;; Spherical Bessel functions of the first kind `spherical-bessel-j` $j_\alpha(x)$, `spherical-bessel-j0`, `spherical-bessel-j1` and `spherical-bessel-j2` are functions of orders `0`,`1` and `2`. Functions are defined for positive argument (only functions with orders `0`, `1` and `2` accept non positive argument).

;; $$j_\alpha(x)=\sqrt{\frac{\pi}{2x}}J_{\alpha+\frac{1}{2}}(x)$$

^:kind/table
[[(gg/->image (gg/functions [["0" special/spherical-bessel-j0]
                             ["1" special/spherical-bessel-j1]
                             ["1.5" (partial special/spherical-bessel-j 1.5)]
                             ["2" special/spherical-bessel-j2]
                             ["3" (partial special/spherical-bessel-j 3)]]
                            {:x [0 20]
                             :legend-name "order"
                             :title "j(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial   special/spherical-bessel-j -0.5)]
                             ["-1" (partial   special/spherical-bessel-j -1)]
                             ["-1.5" (partial special/spherical-bessel-j -1.5)]
                             ["-2" (partial   special/spherical-bessel-j -2)]
                             ["-2.5" (partial special/spherical-bessel-j -2.5)]]
                            {:x [0 20]
                             :ylim [-1 1]
                             :legend-name "order"
                             :title "j(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/spherical-bessel-j0 2.3)
  (special/spherical-bessel-j1 2.3)
  (special/spherical-bessel-j2 2.3)
  (special/spherical-bessel-j 3.1 3.2)
  (special/spherical-bessel-j -3.1 3.2))

;; ---

;; additional `jinc` (sombrero) function is defined as:

;; $$\operatorname{jinc}(x)=\frac{2J_1(\pi x)}{\pi x}$$

(gg/->image (gg/function special/jinc {:x [-5 5] :title "jinc(x)"}))

(utls/examples-note
  (special/jinc 0.0)
  (special/jinc -2.3)
  (special/jinc 2.3))

;; ### Bessel Y, y

;; Bessel functions of the second kind, `bessel-Y`, $Y_\alpha(x)$. `bessel-Y0` and `bessel-Y1` are functions of orders `0` and `1`. They are defined for positive argument only and any order.

^:kind/table
[[(gg/->image (gg/functions [["0" special/bessel-Y0]
                             ["1" special/bessel-Y1]
                             ["1.5" (partial special/bessel-Y 1.5)]
                             ["2" (partial special/bessel-Y 2)]
                             ["3" (partial special/bessel-Y 3)]]
                            {:x [m/MACHINE-EPSILON 20]
                             :ylim [-1 nil]
                             :legend-name "order"
                             :title "Y(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/bessel-Y -0.5)]
                             ["-1" (partial special/bessel-Y -1)]
                             ["-1.5" (partial special/bessel-Y -1.5)]
                             ["-2" (partial special/bessel-Y -2)]
                             ["-2.5" (partial special/bessel-Y -2.5)]]
                            {:x [m/MACHINE-EPSILON 20]
                             :ylim [-1 1]
                             :legend-name "order"
                             :title "Y(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/bessel-Y0 2.3)
  (special/bessel-Y1 2.3)
  (special/bessel-Y 2 2.3)
  (special/bessel-Y -2.1 2.3)
  (special/bessel-Y 3 -1))

;; ---

;; Spherical Bessel functions of the second kind `spherical-bessel-y` $y_\alpha(x)$, `spherical-bessel-y0`, `spherical-bessel-y1` and `spherical-bessel-y2` are functions of orders `0`,`1` and `2`. Functions are defined for positive argument.

;; $$y_\alpha(x)=\sqrt{\frac{\pi}{2x}}Y_{\alpha+\frac{1}{2}}(x)$$

^:kind/table
[[(gg/->image (gg/functions [["0" special/spherical-bessel-y0]
                             ["1" special/spherical-bessel-y1]
                             ["1.5" (partial special/spherical-bessel-y 1.5)]
                             ["2" special/spherical-bessel-y2]
                             ["3" (partial special/spherical-bessel-y 3)]]
                            {:x [m/MACHINE-EPSILON 20]
                             :ylim [-1 nil]
                             :legend-name "order"
                             :title "y(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/spherical-bessel-y -0.5)]
                             ["-1" (partial   special/spherical-bessel-y -1)]
                             ["-1.5" (partial special/spherical-bessel-y -1.5)]
                             ["-2" (partial   special/spherical-bessel-y -2)]
                             ["-2.5" (partial special/spherical-bessel-y -2.5)]]
                            {:x [0.000001 20]
                             :ylim [-1 1]
                             :legend-name "order"
                             :title "y(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/spherical-bessel-j0 2.3)
  (special/spherical-bessel-j1 2.3)
  (special/spherical-bessel-j2 2.3)
  (special/spherical-bessel-j 3.1 3.2)
  (special/spherical-bessel-j -3.1 3.2))

;; ### Bessel I, i

;; Modified Bessel functions of the first kind, `bessel-I`, $I_\alpha(x)$, `bessel-I0` and `bessel-I1` are functions of orders `0` and `1`. An order should be integer for negative arguments.

^:kind/table
[[(gg/->image (gg/functions [["0" special/bessel-I0]
                             ["1" special/bessel-I1]
                             ["1.5" (partial special/bessel-I 1.5)]
                             ["2" (partial special/bessel-I 2)]
                             ["3" (partial special/bessel-I 3)]]
                            {:x [-5.0 5.0]
                             :legend-name "order"
                             :title "I(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/bessel-I -0.5)]
                             ["-1" (partial special/bessel-I -1)]
                             ["-1.5" (partial special/bessel-I -1.5)]
                             ["-2" (partial special/bessel-I -2)]
                             ["-2.5" (partial special/bessel-I -2.5)]]
                            {:x [-5.0 5.0]
                             :ylim [-5.0 5.0]
                             :legend-name "order"
                             :title "I(a,x), negative order"
                             :palette gg/palette-blue-1}))]]



(utls/examples-note
  (special/bessel-I0 2.3)
  (special/bessel-I1 2.3)
  (special/bessel-I 2 2.3)
  (special/bessel-I -2.1 2.3)
  (special/bessel-I -3 -1)
  (special/bessel-I -3.1 -1))

;; ---

;; Two modfified spherical Bessel functions of the first kind `spherical-bessel-1-i` $i_\alpha^{(1)}(x)$ and `spherical-bessel-2-i` $i_\alpha^{(2)}(x)$. `spherical-bessel-1-i0`, `spherical-bessel-1-i1`, `spherical-bessel-1-i2`, `spherical-bessel-2-i0`, `spherical-bessel-2-i1` and `spherical-bessel-2-i2` are functions of orders `0`,`1` and `2`. Functions are defined for positive argument.

;; $$i_\alpha^{(1)}(x)=\sqrt{\frac{\pi}{2x}}I_{\alpha+\frac{1}{2}}(x)$$
;; $$i_\alpha^{(2)}(x)=\sqrt{\frac{\pi}{2x}}I_{-\alpha-\frac{1}{2}}(x)$$

^:kind/table
[[(gg/->image (gg/functions [["0" special/spherical-bessel-1-i0]
                             ["1" special/spherical-bessel-1-i1]
                             ["1.5" (partial special/spherical-bessel-1-i 1.5)]
                             ["2" special/spherical-bessel-1-i2]
                             ["3" (partial special/spherical-bessel-1-i 3)]]
                            {:x [0 5]
                             :ylim [nil 5]
                             :legend-name "order"
                             :title "i(1)(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/spherical-bessel-1-i -0.5)]
                             ["-1" (partial   special/spherical-bessel-1-i -1)]
                             ["-1.5" (partial special/spherical-bessel-1-i -1.5)]
                             ["-2" (partial   special/spherical-bessel-1-i -2)]
                             ["-2.5" (partial special/spherical-bessel-1-i -2.5)]]
                            {:x [0.000001 5]
                             :ylim [-3 7]
                             :legend-name "order"
                             :title "i(1)(a,x), negative order"
                             :palette gg/palette-blue-1}))]
 [(gg/->image (gg/functions [["0" special/spherical-bessel-2-i0]
                             ["1" special/spherical-bessel-2-i1]
                             ["1.5" (partial special/spherical-bessel-2-i 1.5)]
                             ["2" special/spherical-bessel-2-i2]
                             ["3" (partial special/spherical-bessel-2-i 3)]]
                            {:x [0 5]
                             :ylim [-5 5]
                             :legend-name "order"
                             :title "i(2)(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/spherical-bessel-2-i -0.5)]
                             ["-1" (partial   special/spherical-bessel-2-i -1)]
                             ["-1.5" (partial special/spherical-bessel-2-i -1.5)]
                             ["-2" (partial   special/spherical-bessel-2-i -2)]
                             ["-2.5" (partial special/spherical-bessel-2-i -2.5)]]
                            {:x [0.000001 5]
                             :ylim [nil 7]
                             :legend-name "order"
                             :title "i(2)(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/callout "note" "Examples"
  (kind/md "$i_\\alpha^{(1)}$")
  (utls/examples
    (special/spherical-bessel-1-i0 2.3)
    (special/spherical-bessel-1-i1 2.3)
    (special/spherical-bessel-1-i2 2.3)
    (special/spherical-bessel-1-i 3.1 3.2)
    (special/spherical-bessel-1-i -3.1 3.2))
  (kind/md "$i_\\alpha^{(2)}$")
  (utls/examples
    (special/spherical-bessel-2-i0 2.3)
    (special/spherical-bessel-2-i1 2.3)
    (special/spherical-bessel-2-i2 2.3)
    (special/spherical-bessel-2-i 3.1 3.2)
    (special/spherical-bessel-2-i -3.1 3.2)))

;; ### Bessel K, k

;; Modified Bessel functions of the second kind, `bessel-K`, $K_\alpha(x)$, `bessel-K0` and `bessel-K1` are functions of orders `0` and `1`. They are defined for positive argument only and any order.

^:kind/table
[[(gg/->image (gg/functions [["0" special/bessel-K0]
                             ["1" special/bessel-K1]
                             ["1.5" (partial special/bessel-K 1.5)]
                             ["2" (partial special/bessel-K 2)]
                             ["3" (partial special/bessel-K 3)]]
                            {:x [m/MACHINE-EPSILON 5.0]
                             :ylim [nil 5]
                             :legend-name "order"
                             :title "K(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/bessel-K -0.5)]
                             ["-1" (partial special/bessel-K -1)]
                             ["-1.5" (partial special/bessel-K -1.5)]
                             ["-2" (partial special/bessel-K -2)]
                             ["-2.5" (partial special/bessel-K -2.5)]]
                            {:x [m/MACHINE-EPSILON 5.0]
                             :ylim [nil 5.0]
                             :legend-name "order"
                             :title "K(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/bessel-K0 2.3)
  (special/bessel-K1 2.3)
  (special/bessel-K 2 2.3)
  (special/bessel-K -2.1 2.3)
  (special/bessel-K 3 -1))

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

(utls/examples-note
  [(special/bessel-K-half-odd 1 2.3)
   (special/bessel-K 0.5 2.3)]
  [(special/bessel-K-half-odd 3 2.3)
   (special/bessel-K 1.5 2.3)]
  [(special/bessel-K-half-odd 5 2.3)
   (special/bessel-K 2.5 2.3)])

;; ---

;; Modified spherical Bessel functions of the second kind `spherical-bessel-k` $k_\alpha(x)$, `spherical-bessel-k0`, `spherical-bessel-k1` and `spherical-bessel-k2` are functions of orders `0`,`1` and `2`. Functions are defined for positive argument.

;; $$k_\alpha(x)=\sqrt{\frac{\pi}{2x}}K_{\alpha+\frac{1}{2}}(x)$$

^:kind/table
[[(gg/->image (gg/functions [["0" special/spherical-bessel-k0]
                             ["1" special/spherical-bessel-k1]
                             ["1.5" (partial special/spherical-bessel-k 1.5)]
                             ["2" special/spherical-bessel-k2]
                             ["3" (partial special/spherical-bessel-k 3)]]
                            {:x [m/MACHINE-EPSILON 5]
                             :ylim [nil 1]
                             :legend-name "order"
                             :title "k(a,x)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-0.5" (partial special/spherical-bessel-k -0.5)]
                             ["-1" (partial   special/spherical-bessel-k -1)]
                             ["-1.5" (partial special/spherical-bessel-k -1.5)]
                             ["-2" (partial   special/spherical-bessel-k -2)]
                             ["-2.5" (partial special/spherical-bessel-k -2.5)]]
                            {:x [m/MACHINE-EPSILON 5]
                             :ylim [nil 1]
                             :legend-name "order"
                             :title "k(a,x), negative order"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
 (special/spherical-bessel-k0 2.3)
 (special/spherical-bessel-k1 2.3)
 (special/spherical-bessel-k2 2.3)
 (special/spherical-bessel-k 3.1 3.2)
 (special/spherical-bessel-k -3.1 3.2))

;; ## Hankel

;; ::: {.callout-tip title="Defined functions"}
;; * `hankel-1`, `hankel-2`
;; * `spherical-hankel-1`, `spherical-hankel-2`
;; :::

;; Hankel functions of the first and second kind.

;; $$H^{(1)}_\alpha(x)=J_\alpha(x)+ i Y_\alpha(x)$$
;; $$H^{(2)}_\alpha(x)=J_\alpha(x)- i Y_\alpha(x)$$

(utls/examples-note
 (special/hankel-1 1 2.3)
 (special/hankel-2 1 2.3))

;; Spherical functions of the first and second kind.

;; $$h^{(1)}_\alpha(x)=j_\alpha(x)+ i j_\alpha(x)$$
;; $$h^{(2)}_\alpha(x)=j_\alpha(x)- i j_\alpha(x)$$

(utls/examples-note
 (special/spherical-hankel-1 1 2.3)
 (special/spherical-hankel-2 1 2.3))

;; ## Erf

;; ::: {.callout-tip title="Defined functions"}
;; * `erf`, `erfc`
;; * `inv-erf`, `inv-erfc`
;; :::

;; Error functions

;; $$\operatorname{erf}(x)=\frac{2}{\sqrt\pi}\int_0^x e^{-t^2}\,dt$$
;; $$\operatorname{erfc}(x)=1-\operatorname{erf}(x)$$

^:kind/table
[[(gg/->image (gg/function special/erf {:x [-3 3] :title "erf(x)"}))
  (gg/->image (gg/function special/erfc {:x [-3 3] :title "erfc(x)"}))]]

;; When two arguments are passed, difference between erf of two values is calculated $\operatorname{erf}(x_2)-\operatorname{erf}(x_1)$

(utls/examples-note
  (special/erf 0.4)
  (special/erfc 0.4)

  [(special/erf 0.5 0.4)
   (- (special/erf 0.4) (special/erf 0.5))])

;; ---

;; Inverse of error functions.

;; * `inv-erf` - inverse of `erf` defined on $(-1,1)$
;; * `inv-erfc`- inverse of `erfc` defined on $(0,2)$

^:kind/table
[[(gg/->image (gg/function special/inv-erf {:x [-1.1 1.1] :ylim [-2 2] :title "inv-erf(x)"}))
  (gg/->image (gg/function special/inv-erfc {:x [-0.1 2.1] :ylim [-2 2] :title "inv-erfc(x)"}))]]

(utls/examples-note
  (special/inv-erf 0.42839235504666855)
  (special/inv-erfc (- 1 0.42839235504666855)))

;; ## Airy

;; Airy functions and derivatives

;; ::: {.callout-tip title="Defined functions"}
;; * `airy-Ai`, `airy-Bi`
;; * `airy-Ai'` `airy-Bi'`
;; :::

^:kind/table
[[(gg/->image (gg/functions [["Ai" special/airy-Ai]
                             ["Bi" special/airy-Bi]]
                            {:x [-10 5]
                             :ylim [nil 2]
                             :title "Airy functions"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["Ai'" special/airy-Ai']
                             ["Bi'" special/airy-Bi']]
                            {:x [-10 5]
                             :ylim [nil 2]
                             :title "Derivatives of Airy functions"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/airy-Ai 2.3)
  (special/airy-Bi 2.3)
  (special/airy-Ai' 2.3)
  (special/airy-Bi' 2.3))

;; ## Integrals

;; Trigonometric, exponential and logarithmic integrals

;; ::: {.callout-tip title="Defined functions"}
;; * `Si`, `si`, `Ci`, `Cin`
;; * `E0`, `E1`, `Ei`, `Ein`, `En`
;; * `li`, `Li` (offset)
;; :::

;; ### Trigonometric

;; Sine and cosine integrals

;; $$\operatorname{Si}(x)=\int_0^x\frac{\sin t}{t}\, dt$$
;; $$\operatorname{si}(x)=-\int_x^\infty\frac{\sin t}{t}\, dt = \operatorname{Si}(x)-\frac{\pi}{2}$$
;; $$\operatorname{Ci}(x)=-\int_x^\infty\frac{\cos t}{t}\, dt$$
;; $$\operatorname{Cin}(x)=\int_0^x\frac{1-\cos t}{t}\, dt = \gamma + \ln x- \operatorname{Ci}(x)$$

^:kind/table
[[(gg/->image (gg/functions [["Si" special/Si]
                             ["si" special/si]]
                            {:x [-15 15]
                             :title "Sine integrals"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["Ci" special/Ci]
                             ["Cin" special/Cin]]
                            {:x [0 15]
                             :title "Cosine integrals"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/Si 0.5)
  (special/si 0.5)
  (special/Ci 0.5)
  (special/Cin 0.5))

;; ### Exponential

;; Exponential integrals

;; $$E_0(x)=\frac{e^{-x}}{x}$$
;; $$E_1(x)=\int_x^\infty\frac{e^{-t}}{t}\,dt$$
;; $$E_i(x)=-\int_{-x}^\infty\frac{e^{-t}}{t}\,dt$$
;; $$E_{in}(x)=\int_0^x\frac{1-e^{-t}}{t}\,dt$$
;; $$E_n(x)=\int_1^\infty\frac{e^{-xt}}{t^n}\,dt$$


^:kind/table
[[(gg/->image (gg/functions [["E0" special/E0]
                             ["E1" special/E1]]
                            {:x [0.0001 5]
                             :ylim [nil 2]
                             :title "E0 and E1 exponential integrals"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["0" (partial special/En 0)]
                             ["1" (partial special/En 1)]
                             ["1.5" (partial special/En 1.5)]
                             ["2" (partial special/En 2)]
                             ["-1" (partial special/En -1)]
                             ["-2" (partial special/En -2)]]
                            {:x [0.0001 5]
                             :ylim [nil 2]
                             :title "En exponential integrals"
                             :legend-name "n"
                             :palette gg/palette-blue-1}))]
 [(gg/->image (gg/function special/Ei
                           {:x [-5 2.5]
                            :ylim [-5 5]
                            :title "Ei exponential integral"}))
  (gg/->image (gg/function special/Ein
                           {:x [0 10]
                            :title "Ein exponential integral"}))]]

(utls/examples-note
  (special/E0 0.6)
  (special/E1 0.6)
  (special/Ei 0.6)
  (special/Ein 0.6)
  (special/En 2 0.6)
  (special/En -2 0.6))

;; ### Logarithmic

;; Logarithmic integrals

;; $$\operatorname{li}(x)=\int_0^x\frac{dt}{\ln t}$$
;; $$\operatorname{Li}(x)=\int_2^x\frac{dt}{\ln t}=\operatorname{li}(x)-\operatorname{li}(2)$$

(gg/->image (gg/functions [["li" special/li]
                           ["Li" special/Li]]
                          {:x [0 5]
                           :title "li and Li logarithmic intergrals"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/li 0.5)
  (special/Li 0.5))

;; ## Zeta

;; Zeta function and related

;; ::: {.callout-tip title="Defined functions"}
;; * `zeta` - Riemann and Hurwitz zeta
;; * `xi` - Riemann (Landau) xi
;; * `eta` - Dirichlet eta
;; * `dirichlet-beta` - Dirichlet beta
;; :::

;; ### Riemann zeta

;; $$\zeta(s)=\sum_{n=1}^\infty\frac{1}{n^s}$$

(gg/->image (gg/function special/zeta
                         {:x [-5 7]
                          :ylim [-5 5]
                          :title "Riemann zeta"}))

(utls/examples-note
  (special/zeta 0.0)
  (special/zeta 2.2)
  (special/zeta -2.2))

;; ### Hurwitz zeta

;; $$\zeta(s,z)=\sum_{n=1}^\infty\frac{1}{(n+z)^s}$$

^:kind/table
[[(gg/->image (gg/functions [["1" #(special/zeta % 1)]
                             ["2.5" #(special/zeta % 2.5)]
                             ["5" #(special/zeta % 5)]]
                            {:x [-5 7]
                             :ylim [-7 7]
                             :title "Hurwitz zeta"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["-1.5" #(special/zeta % -1.5)]
                             ["-2.5" #(special/zeta % -2.5)]
                             ["-5" #(special/zeta % -5)]]
                            {:x [-5 5]
                             :ylim [-10 20]
                             :title "Hurwitz zeta"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/zeta 0.0 3)
  (special/zeta 2.2 3)
  (special/zeta -2.2 3)
  (special/zeta 2.2 -3)
  (special/zeta -2.2 -3)
  (special/zeta 0.0 -3))

;; ### xi

;; Landau's Xi function, symmetrical along $x=0.5$

;; $$\xi(s)=\frac{1}{2}s(s-1)\pi^{-\frac{s}{2}}\Gamma\left(\frac{s}{2}\right)\zeta(s)$$
;; $$\xi(s)=\xi(1-s)$$

(gg/->image (gg/function special/xi
                         {:x [-5 5]
                          :title "Landau's xi"}))

(utls/examples-note
  (special/xi 0.0)
  (special/xi 3.5)
  (special/xi (- 1.0 3.5)))

;; ### eta

;; Dirichlet eta function

;; $$\eta(s)=\sum_{n=1}^\infty\frac{(-1)^{n-1}}{n^s}=(1-2^{1-s})\zeta(s)$$

(gg/->image (gg/function special/eta
                         {:x [-8.5 5]
                          :title "Dirichlet eta"}))

(utls/examples-note
  (special/eta 0.0)
  (special/eta 3.5)
  (special/eta -3.5))

;; ### beta

;; Dirichlet (Catalan) beta function

;; $$\beta(s)=\sum_{n=0}^\infty\frac{(-1)^n}{(2n+1)^s}$$

(gg/->image (gg/function special/dirichlet-beta
                         {:x [-5 5]
                          :title "Dirichlet beta"}))

(utls/examples-note
  (special/dirichlet-beta 0.0)
  (special/dirichlet-beta 3.5)
  (special/dirichlet-beta -3.5))

;; ## Hypergeometric

;; Selection of hypergeometric functions ${}_pF_q$

;; $${{}_{p}F_{q}}\left({a_{1},\dots,a_{p}\atop b_{1},\dots,b_{q}};x\right)=\sum_{k=0}^{\infty}\frac{{\left(a_{1}\right)_{k}}\cdots{\left(a_{p}\right)_{k}}}{{\left(b_{1}\right)_{k}}\cdots{\left(b_{q}\right)_{k}}}\frac{x^{k}}{k!}.$$

;; where $(a_p)_k$ and $(b_q)_k$ are k^th^ rising factorials

;; ::: {.callout-tip title="Defined functions"}
;; * `hypergeometric-pFq`, `hypergeometric-pFq-ratio`, `hypergeometric-pFq-complex`
;; * `hypergeometric-0F0`, `hypergeometric-0F1`, `hypergeometric-0F2`
;; * `hypergeometric-1F0`, `hypergeometric-1F1`
;; * `hypergeometric-2F0`, `hypergeometric-2F1`
;; * `kummers-M`, `tricomis-U`, `tricomis-U-complex`
;; * `whittaker-M`, `whittaker-W`
;; :::

;; Functions are implemented using various recursive formulas, Maclaurin series and Weniger acceleration.

;; ### ~p~F~q~, generalized

;; Two implementations of general ${}_pF_q$ hypergeometric functions using Maclaurin series. One implementation operates on doubles (`hypergeometric-pFq`), second on Clojure `ratio` type which is accurate but slow (`hypergeometric-pFq-ratio`), third on complex numbers..

(utls/examples-note
  (special/hypergeometric-2F0 0.1 0.1 0.01)
  (special/hypergeometric-pFq [0.1 0.1] [] 0.01)
  (special/hypergeometric-pFq-ratio [0.1 0.1] [] 0.01)
  (special/hypergeometric-pFq-complex [(complex/complex 0.2 0.2)] [1 (complex/complex -1 1)] 0.1))

;; Both functions accept optional `max-iters` parameter to control number of iterations.

;; Every implementation but ratio is unstable. Take a look at the example of ${}_1F_1(-50;3;19.5)$, only ratio version gives a valid result.

(utls/examples-note
  (special/hypergeometric-1F1 -50 3 19.5)
  (special/hypergeometric-pFq [-50] [3] 19.5)
  (special/hypergeometric-pFq-ratio [-50] [3] 19.5))

;; Following plot shows stable (but slow) implementation `hypergeometric-pFq-ratio` vs unstable (but fast) `hypergeometric-1F1` and Maclaurin seriers `hypergeometric-pFq`.

(gg/->image (gg/functions [["stable" (partial special/hypergeometric-pFq-ratio [-50] [3])]
                           ["1F1" (partial special/hypergeometric-1F1 -50 3)]
                           ["pFq" (partial special/hypergeometric-pFq [-50] [3])]]
                          {:x [0.5 15]
                           :ylim [-1 1]
                           :title "Numerical stability of hypergeometric implementations"}))

;; ### ~0~F~0~, exp

;; ${}_0F_0$ is simply exponential function.

;; $${}_0F_0(;;x)=\sum_{k=0}^\infty\frac{x^k}{k!}=e^x$$

(utls/examples-note
  (special/hypergeometric-0F0 2.3)
  (m/exp 2.3))

;; ### ~0~F~1~

;; ${}_0F_1$ is called confluent hypergeometric limit function.

;; $${}_0F_1(;a;x)=\sum_{k=0}^\infty\frac{x^k}{(a)_k k!}=
;; \begin{cases}
;; 1.0 & x=0 \\
;; \frac{J_{a-1}\left(2\sqrt{|x|}\right)\Gamma(a)}{|x|^\frac{a-1}{2}} & x<0 \\
;; \frac{I_{a-1}\left(2\sqrt{|x|}\right)\Gamma(a)}{|x|^\frac{a-1}{2}} & x>0
;; \end{cases}
;; $$

(gg/->image (gg/functions [["-1.5" (partial special/hypergeometric-0F1 -1.5)]
                           ["-0.5" (partial special/hypergeometric-0F1 -0.5)]
                           ["0.5" (partial special/hypergeometric-0F1 0.5)]
                           ["1" (partial special/hypergeometric-0F1 1)]
                           ["2" (partial special/hypergeometric-0F1 2)]]
                          {:x [-6 3]
                           :ylim [-6 6]
                           :title "0F1"
                           :legend-name "a"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [] [-0.5] -2)
  (special/hypergeometric-0F1 -0.5 -2))

;; ### ~0~F~2~

;; $${}_0F_2(;a,b;x)=\sum_{k=0}^\infty\frac{x^k}{(a)_k(b)_k k!}$$

(gg/->image (gg/functions [["(-1.5,-1.5)" (partial special/hypergeometric-0F2 -1.5 -1.5)]
                           ["(-0.5,0.5)" (partial special/hypergeometric-0F2 -0.5 0.5)]
                           ["(0.5,0.5)" (partial special/hypergeometric-0F2 0.5 0.5)]
                           ["(1.5,1.5)" (partial special/hypergeometric-0F2 1.5 1.5)]]
                          {:x [-6 3]
                           :ylim [-6 6]
                           :title "0F2"
                           :legend-name "(a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [] [-0.5, 0.5] -2)
  (special/hypergeometric-0F2 -0.5 0.5 -2))

;; ### ~1~F~0~

;; $${}_1F_0(a;;x)=\sum_{k=0}^\infty\frac{(a)_k x^k}{k!}=(1-x)^{-a}$$

(gg/->image (gg/functions [["-1.5" (partial special/hypergeometric-1F0 -1.5)]
                           ["-0.5" (partial special/hypergeometric-1F0 -0.5)]
                           ["0.5" (partial special/hypergeometric-1F0 0.5)]
                           ["1" (partial special/hypergeometric-1F0 1)]
                           ["2" (partial special/hypergeometric-1F0 2)]]
                          {:x [-6 3]
                           :ylim [-6 6]
                           :title "1F0"
                           :legend-name "a"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [0.5] [] -0.5)
  (special/hypergeometric-1F0 0.5 -0.5))

;; ### ~1~F~1~, M

;; Confluent hypergeometric function of the first kind, Kummer's M.

;; $${}_1F_1(a;b;x)=M(a,b,x)=\sum_{k=0}^\infty\frac{(a)_k x^k}{(b)_k k!}$$

(gg/->image (gg/functions [["(-1.5,1)" (partial special/hypergeometric-1F1 -1.5 1)]
                           ["(-1.5,-1.5)" (partial special/hypergeometric-1F1 -1.5 -1.5)]
                           ["(0.2, 0.3)" (partial special/hypergeometric-1F1 0.2 0.3)]
                           ["(0.2, -0.3)" (partial special/hypergeometric-1F1 0.2 -0.3)]]
                          {:x [-3 3]
                           :ylim [-6 6]
                           :title "1F1, Kummer's M"
                           :legend-name "(a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [0.5] [1] -0.5)
  (special/hypergeometric-1F1 0.5 1 -0.5)
  (special/kummers-M 0.5 1 -0.5))

;; ### ~2~F~0~, U

;; ${}_2F_0$ is related to the confluent hypergeometric function of the second kind, Tricomi's U

;; $${}_2F_0(a,b;;x)=\sum_{k=0}^\infty\frac{(a)_k (b)_k x^k}{k!}$$

(gg/->image (gg/functions [["(0.1, -0.1)" (partial special/hypergeometric-1F1 0.1 -0.1)]
                           ["(-0.1, 0.1)" (partial special/hypergeometric-1F1 -0.1 0.1)]]
                          {:x [-3 2]
                           :ylim [-6 nil]
                           :title "2F0"
                           :legend-name "(a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [0.1 -0.1] [] 0.01)
  (special/hypergeometric-2F0 0.1 -0.1 0.01)
  (special/hypergeometric-2F0 0.1 -0.1 1.2))


;; ---

;; $$U(a,b,x) \sim x^{-a}{}_2F_0(a,b;;-\frac{1}{x})$$

(gg/->image (gg/functions [["(0.1, -0.1)" (partial special/tricomis-U 0.1 -0.1)]
                           ["(-0.1, 0.1)" (partial special/tricomis-U -0.1 0.1)]
                           ["(0.5, 0.5)" (partial special/tricomis-U 0.5 0.5)]
                           ["(-0.5, -0.5)" (partial special/tricomis-U -0.5 -0.5)]]
                          {:x [0 4]
                           #_#_:ylim [-6 nil]
                           :title "Tricomi's U"
                           :legend-name "(a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/tricomis-U 0.5 0.2 0.1)
  (special/tricomis-U -0.5 0.2 0.1))

;; Complex variant accepts both complex numbers and real numbers as arguments.

(gg/->image (gg/complex-function (partial special/tricomis-U-complex (complex/complex -1 1)
                                          (complex/complex 1 -1))
                                 {:x [-2.5 3.5] :y [-3 3]
                                  :steps 500 :title "Complex Tricomis U (i-1, 1-i)"}))

(utls/examples-note
  (special/tricomis-U-complex 0.5 0.2 0.1)
  (special/tricomis-U-complex (complex/complex -1 1) (complex/complex 1 -1) 1))

;; ### ~2~F~1~, Gauss

;; Gauss' hypergeometric function ${}_2F_1$.

;; $${}_2F_1(a,b;c;x)=\sum_{k=0}^\infty\frac{(a)_k (b)_k x^k}{(c)_k k!}$$

(gg/->image (gg/functions [["(1.1, -0.1, 0.2)" (partial special/hypergeometric-2F1 1.1 -0.1 0.2)]
                           ["(-0.1, 0.1, -1.2)" (partial special/hypergeometric-2F1 -0.1 0.1 -1.2)]
                           ["(0.5, 2, 0.5)" (partial special/hypergeometric-2F1 0.5 2 0.5)]
                           ["(-0.5,-3, -0.5)" (partial special/hypergeometric-2F1 -0.5 -3 -0.5)]]
                          {:x [-2 0.9999]
                           :ylim [-2 2]
                           :title "2F1, Gauss hypergeometric function"
                           :legend-name "(a,b)"
                           :palette gg/palette-blue-1}))

(utls/examples-note
  (special/hypergeometric-pFq-ratio [1.1 -0.1] [0.2] 0.5)
  (special/hypergeometric-2F1 1.1 -0.1 0.2 0.5))

;; ### Whittaker M and W

;; Modified hypergeometric functions by Whittaker

;; $$M_{\kappa,\mu}\left(x\right)=e^{-\frac{1}{2}x}x^{\frac{1}{2}+\mu}M\left(\tfrac{1}{2}+\mu-\kappa,1+2\mu,x\right)$$
;; $$W_{\kappa,\mu}\left(x\right)=e^{-\frac{1}{2}x}x^{\frac{1}{2}+\mu}U\left(\tfrac{1}{2}+\mu-\kappa,1+2\mu,x\right)$$

^:kind/table
[[(gg/->image (gg/functions [["(0.5, -0.4)" (partial special/whittaker-M 0.5 -0.25)]
                             ["(0, 1)" (partial special/whittaker-M 0 1)]
                             ["(-0.5, 1)" (partial special/whittaker-M -0.5 1)]
                             ["(-0.2, -1.2)" (partial special/whittaker-M -0.2 -1.2)]]
                            {:x [0 5]
                             :ylim [-6 10]
                             :title "Whittaker M"
                             :legend-name "(kappa, mu)"
                             :palette gg/palette-blue-1}))
  (gg/->image (gg/functions [["(0.5, -0.4)" (partial special/whittaker-W 0.5 -0.25)]
                             ["(0, 1)" (partial special/whittaker-W 0 1)]
                             ["(-0.5, 1)" (partial special/whittaker-W -0.5 1)]
                             ["(-0.2, -1.2)" (partial special/whittaker-W -0.2 -1.2)]]
                            {:x [0 5]
                             :ylim [nil 3]
                             :title "Whittaker W"
                             :legend-name "(kappa, mu)"
                             :palette gg/palette-blue-1}))]]

(utls/examples-note
  (special/whittaker-M 0.3 0.4 1.2)
  (special/whittaker-W 0.3 0.4 1.2))

;; ## Other

;; ::: {.callout-tip title="Defined functions"}
;; * `lambert-W` ($W_0$), `lambert-W-1` ($W_{-1}$)
;; * `harmonic-number`
;; * `minkowski` - $?(x)$
;; * `owens-t`
;; :::

;; ### Lambert W

;; Lambert W is a function for which $W(xe^x)=x$. There are two branches $W_0$ (`lambert-W`) and $W_{-1}$ (`lambert-W-1`).

;; $$
;; \begin{align}
;; W_0(xe^x)=x & \text{ for } x\ge -1 \\
;; W_{-1}(xe^x)=x & \text{ for } x\le -1
;; \end{align}
;; $$

;; domain of functions

;; $$
;; \begin{align}
;; W_0(t) & \text{ for } t\in(-1/e,\infty) \\
;; W_{-1}(t) & \text{ for } t\in(-1/e,0)
;; \end{align}
;; $$


(gg/->image (gg/functions [["W0" special/lambert-W]
                           ["W-1" special/lambert-W-1]]
                          {:x [-1 5]
                           :ylim [-5 nil]
                           :steps 2000}))

(utls/examples-note
  (special/lambert-W m/E)
  (special/lambert-W (* 2.3 (m/exp 2.3)))
  (special/lambert-W-1 (* -2 (m/exp -2))))

;; ### Harmonic H

;; Harmonic numbers

;; $$H_n=\int_0^1\frac{1-x^n}{1-x}\,dx=\operatorname{digamma}(x+1)+\gamma$$

;; For non-negative integers

;; $$H_n=\sum_{k=1}^n\frac{1}{k}$$

(gg/->image (gg/function special/harmonic-number
                         {:x [-1 5]
                          :ylim [-2 nil]
                          :title "Harmonic number"}))

(utls/examples-note
  (special/harmonic-number 2)
  (special/harmonic-number 2.5)
  (special/harmonic-number 3)
  (special/harmonic-number -0.5))

;; ---

;; Generalized harmonic numbers for $m\neq0$ or $m\neq1$

;; $$H_{n,m}=\zeta(m)-\zeta(m,n+1)$$
;; $$H_{n,0}=n\text{, }H_{n,1}=H_n$$

;; For non-negative integer `n`

;; $$H_{n,m}=\sum_{k=1}^n\frac{1}{k^m}$$

^:kind/table
[[(gg/->image (gg/functions [["0.5" #(special/harmonic-number % 0.5)]
                             ["1.5" #(special/harmonic-number % 1.5)]
                             ["2.5" #(special/harmonic-number % 2.5)]]
                            {:x [-1 5]
                             :ylim [-2 nil]
                             :legend-name "m"
                             :palette gg/palette-blue-1
                             :title "Generalized armonic number"}))
  (gg/->image (gg/functions [["-0.5" #(special/harmonic-number % -0.5)]
                             ["-1.5" #(special/harmonic-number % -1.5)]
                             ["-2.5" #(special/harmonic-number % -2.5)]]
                            {:x [-1 3]
                             :ylim [nil 5]
                             :legend-name "m"
                             :palette gg/palette-blue-1
                             :title "Generalized armonic number, negative m"}))]]

(utls/examples-note
  (special/harmonic-number 2.2 -0.5)
  (special/harmonic-number 2.2 0.5))

;; ### Minowski

;; Minkowski's question mark $?(x)$ function.

(gg/->image (gg/function special/minkowski
                         {:title "Minkowski function"
                          :steps 1000}))

(utls/examples-note
  (special/minkowski 0.5)
  [(special/minkowski 0.2)
   (- 1.0 (special/minkowski 0.8))]

  [(special/minkowski (/ 0.5 1.5))
   (/ (special/minkowski 0.5) 2)])

;; ### Owen's T

;; Owen's T function

;; $$T(h,a) = \frac{1}{2\pi } \int_{0}^{a} \frac{e^{-\frac{1}{2}h^2(1+x^2)}}{1+x^2}dx\quad(-\infty < h,a < +\infty)$$

(utls/examples-note
  (special/owens-t 1 0)
  (special/owens-t 0 1)
  (special/owens-t 1 1)
  (special/owens-t -1 1)
  (special/owens-t 1 -1)
  (special/owens-t -1 -1))

(gg/->image (gg/function2d special/owens-t {:varg? false :x [-2 2] :y [-2 2]
                                            :xlab "h" :ylab "a"
                                            :title "Owen's T"}))

^:kind/table
[[(gg/->image (gg/functions [["0.0" (partial special/owens-t 0.0)]
                             ["0.5" (partial special/owens-t 0.5)]
                             ["1.0" (partial special/owens-t 1.0)]
                             ["2.0" (partial special/owens-t 2.0)]]
                            {:x [-6 6]
                             :legend-name "h"
                             :xlab "a"
                             :palette gg/palette-blue-1
                             :title "Owen's T (x=a)"}))
  (gg/->image (gg/functions [["1.0" #(special/owens-t % 1.0)]
                             ["2.0" #(special/owens-t % 2.0)]
                             ["5.0" #(special/owens-t % 5.0)]]
                            {:x [-3 3]
                             :legend-name "a"
                             :xlab "h"
                             :palette gg/palette-blue-1
                             :title "Owen's T (x=h)"}))]]

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.special)
