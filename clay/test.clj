^{:kindly/hide-code true :clay {:quarto {:title "My title for this namespace"}}}
(ns clay.test
  (:require [scicloj.clay.v2.api :as clay]
            [scicloj.kindly.v4.kind :as kind]))

;; # ABC

(require '[fastmath.calculus :as calc])

;; ## Richardson extrapolation

;; [Richardson extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation) is a way to increase accuracy and convergence of limits of the form:

;; $$f(x) = \lim_{h\to 0}g(x,h)$$

;; If $x = \infty$ the following limit is calculated (similar for negative infinity):

;; $$f(\infty) = \lim_{h\to 0}g(0,\frac{1}{h})$$

;; There are two helpers which create $g(x,h)$ from any single arity function: $g_+(x,h) = f(x+h)$ and $g_-(x,h) = f(x-h)$.

^:kind/table
{:column-names [:zzz/symbol-> "description" "blah"]
 :row-vectors [['extrapolate "Create extrapolated function f(x) from a function g(x,h)"]
               ['fx->gx+h "Create g+" (kind/html "<span data-qmd=\"$\\infty$\"></span>")]
               ['fx->gx+h "Create g+" (kind/md "$\\infty$")]
               ['fx->gx-h "Create g-"]]}

(kind/hiccup [:span {:data-qmd "$\\infty$"}]  {:element/max-height nil})

;; | fruit  | price  |
;; |--------|--------|
;; | apple  | 2.05   |
;; | pear   | $\infty+1$   |
;; | orange | 3.09   |

(kind/md "$x = \\infty$")


;; List of the functions:

;; `extrapolate` accepts following options:
;;

;; * `:contract`: `h` shrinkage factor, default=`1/2`
;; * `:power`: set to `2.0` for even functions around `x0`, default `1.0`
;; * `:init-h` - initial step `h`, default=`1/2`
;; * `:abs` - absolute error, default: machine epsilon
;; * `:rel` - relative error, default: ulp for init-h
;; * `:tol` - tolerance for current error, default: `2.0` 
;; * `:max-evals` - maximum evaluations, default: maximum integer

;;

;; sadfasdfs^[Inlines notes are easier to write,
;; since you don't have to pick an identifier and move down to
;; type the note.]

;; $$f(\infty) = \lim_{h\to 0}g\left(0,\frac{1}{h}\right)$$

['some->fancy-symbol']


:some->keyword

:namespaced/keyword

;; ### dsfs

