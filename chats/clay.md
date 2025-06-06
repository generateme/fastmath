### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace. 

### Please write a short paragraph for `clay/complex_quaternions.clj` book. Do not copy content of the `clay/complex_quaternions.clj` file, write only asked chapter. Put text as a Clojure comment block. Put actual Clojure code without comment block and markdown wrapper.  Follow the convention of other chapters.

- Write a paragraph for the "### Trigonometric" chapter.
- Describe all mentioned functions, group them in points and describe them in brief
- When possible add latex formulas
- Compare these functions
- Illustrate with examples

Write only asked chapter. Refer to previous chapters as examples of chapter structure.



```clojure
;; ### Trigonometric
;;
;; Complex trigonometric and hyperbolic functions extend their real-valued counterparts to the complex plane. They are defined in terms of the complex exponential function:
;;
;; $$ \sin z = \frac{e^{iz} - e^{-iz}}{2i} \quad \cos z = \frac{e^{iz} + e^{-iz}}{2} $$
;; $$ \sinh z = \frac{e^{z} - e^{-z}}{2} \quad \cosh z = \frac{e^{z} + e^{-z}}{2} $$
;;
;; This leads to relationships between complex trigonometric and hyperbolic functions, such as $\sin(z) = -i \sinh(iz)$ and $\cos(z) = \cosh(iz)$. These functions are implemented using formulas that provide numerical stability, especially for large arguments.
;;
;; ::: {.callout-tip title="Defined functions"}
;; * `sin`, `cos`, `tan`
;; * `sec`, `csc`, `cot`
;; * `sinh`, `cosh`, `tanh`
;; * `sech`, `csch`, `coth`
;; * `asin`, `acos`, `atan`
;; * `asec`, `acsc`, `acot`
;; * `asinh`, `acosh`, `atanh`
;; * `asech`, `acsch`, `acoth`
;; :::
;;
;; *   **Standard Trigonometric Functions**: `sin`, `cos`, `tan`, `sec`, `csc`, `cot`
;;     *   These are the direct extensions of real trigonometric functions. For $z = x + iy$:
;;         $$ \sin z = \sin x \cosh y + i \cos x \sinh y $$
;;         $$ \cos z = \cos x \cosh y - i \sin x \sinh y $$
;;         $$ \tan z = \frac{\sin x \cos x + i \sinh y \cosh y}{\cos^2 x + \sinh^2 y} $$
;;     *   `sec(z)`, `csc(z)`, and `cot(z)` are computed as the reciprocals of `cos(z)`, `sin(z)`, and `tan(z)`, respectively.
;;
;; *   **Hyperbolic Trigonometric Functions**: `sinh`, `cosh`, `tanh`, `sech`, `csch`, `coth`
;;     *   These are also extensions of real hyperbolic functions. For $z = x + iy$:
;;         $$ \sinh z = \sinh x \cos y + i \cosh x \sin y $$
;;         $$ \cosh z = \cosh x \cos y + i \sinh x \sin y $$
;;         $$ \tanh z = \frac{\sinh x \cosh x + i \sin y \cos y}{\sinh^2 x + \cos^2 y} $$
;;     *   `sech(z)`, `csch(z)`, and `coth(z)` are computed as the reciprocals of `cosh(z)`, `sinh(z)`, and `tanh(z)`, respectively.
;;
;; *   **Inverse Trigonometric Functions**: `asin`, `acos`, `atan`, `asec`, `acsc`, `acot`
;;     *   These are the multi-valued inverse functions. The implementations return the principal values, which are related to complex logarithms. For example:
;;         $$ \arcsin z = -i \log(iz + \sqrt{1-z^2}) $$
;;         $$ \arccos z = -i \log(z + i\sqrt{1-z^2}) $$
;;         $$ \arctan z = \frac{i}{2} \log\left(\frac{1-iz}{1+iz}\right) $$
;;
;; *   **Inverse Hyperbolic Trigonometric Functions**: `asinh`, `acosh`, `atanh`, `asech`, `acsch`, `acoth`
;;     *   These are the multi-valued inverse hyperbolic functions. The implementations return the principal values, which are also related to complex logarithms. For example:
;;         $$ \operatorname{arsinh} z = \log(z + \sqrt{z^2+1}) $$
;;         $$ \operatorname{arcosh} z = \log(z + \sqrt{z^2-1}) $$
;;         $$ \operatorname{artanh} z = \frac{1}{2} \log\left(\frac{1+z}{1-z}\right) $$
;;
;; These functions are continuous and analytic everywhere in the complex plane, except at certain points (poles) where the denominator is zero for tangent, cotangent, secant, cosecant, and at branch cuts for the inverse functions.
;;

(utls/examples-note
  (cplx/sin z1)
  (cplx/cos z1)
  (cplx/tan z1)
  (cplx/sinh z2)
  (cplx/cosh z2)
  (cplx/tanh z2)
  (cplx/asin z3)
  (cplx/acos z3)
  (cplx/atan z3))

(kind/table
 [[(gg/->image (gg/complex-function cplx/sin {:title "f(z)=sin(z)"}))
   (gg/->image (gg/complex-function cplx/cos {:title "f(z)=cos(z)"}))]
  [(gg/->image (gg/complex-function cplx/tan {:title "f(z)=tan(z)"}))
   (gg/->image (gg/complex-function cplx/sinh {:title "f(z)=sinh(z)"}))]
  [(gg/->image (gg/complex-function cplx/cosh {:title "f(z)=cosh(z)"}))
   (gg/->image (gg/complex-function cplx/tanh {:title "f(z)=tanh(z)"}))]
  [(gg/->image (gg/complex-function cplx/asin {:title "f(z)=asin(z)"}))
   (gg/->image (gg/complex-function cplx/acos {:title "f(z)=acos(z)"}))]
  [(gg/->image (gg/complex-function cplx/atan {:title "f(z)=atan(z)"}))]])
```

### What about formulas for inverse functions: `asec`, `acsc`, `acot`, `asech`, `acsch`, `acoth`?



```clojure
;; ### Trigonometric
;;
;; Complex trigonometric and hyperbolic functions extend their real-valued counterparts to the complex plane. They are defined in terms of the complex exponential function:
;;
;; $$ \sin z = \frac{e^{iz} - e^{-iz}}{2i} \quad \cos z = \frac{e^{iz} + e^{-iz}}{2} $$
;; $$ \sinh z = \frac{e^{z} - e^{-z}}{2} \quad \cosh z = \frac{e^{z} + e^{-z}}{2} $$
;;
;; This leads to relationships between complex trigonometric and hyperbolic functions, such as $\sin(z) = -i \sinh(iz)$ and $\cos(z) = \cosh(iz)$. These functions are implemented using formulas that provide numerical stability, especially for large arguments.
;;
;; ::: {.callout-tip title="Defined functions"}
;; * `sin`, `cos`, `tan`
;; * `sec`, `csc`, `cot`
;; * `sinh`, `cosh`, `tanh`
;; * `sech`, `csch`, `coth`
;; * `asin`, `acos`, `atan`
;; * `asec`, `acsc`, `acot`
;; * `asinh`, `acosh`, `atanh`
;; * `asech`, `acsch`, `acoth`
;; :::
;;
;; *   **Standard Trigonometric Functions**: `sin`, `cos`, `tan`, `sec`, `csc`, `cot`
;;     *   These are the direct extensions of real trigonometric functions. For $z = x + iy$:
;;         $$ \sin z = \sin x \cosh y + i \cos x \sinh y $$
;;         $$ \cos z = \cos x \cosh y - i \sin x \sinh y $$
;;         $$ \tan z = \frac{\sin x \cos x + i \sinh y \cosh y}{\cos^2 x + \sinh^2 y} $$
;;     *   `sec(z)`, `csc(z)`, and `cot(z)` are computed as the reciprocals of `cos(z)`, `sin(z)`, and `tan(z)`, respectively.
;;
;; *   **Hyperbolic Trigonometric Functions**: `sinh`, `cosh`, `tanh`, `sech`, `csch`, `coth`
;;     *   These are also extensions of real hyperbolic functions. For $z = x + iy$:
;;         $$ \sinh z = \sinh x \cos y + i \cosh x \sin y $$
;;         $$ \cosh z = \cosh x \cos y + i \sinh x \sin y $$
;;         $$ \tanh z = \frac{\sinh x \cosh x + i \sin y \cos y}{\sinh^2 x + \cos^2 y} $$
;;     *   `sech(z)`, `csch(z)`, and `coth(z)` are computed as the reciprocals of `cosh(z)`, `sinh(z)`, and `tanh(z)`, respectively.
;;
;; *   **Inverse Trigonometric Functions**: `asin`, `acos`, `atan`, `asec`, `acsc`, `acot`
;;     *   These are the multi-valued inverse functions. The implementations return the principal values, which are related to complex logarithms. For example:
;;         $$ \arcsin z = -i \log(iz + \sqrt{1-z^2}) $$
;;         $$ \arccos z = -i \log(z + i\sqrt{1-z^2}) $$
;;         $$ \arctan z = \frac{i}{2} \log\left(\frac{1-iz}{1+iz}\right) $$
;;         $$ \operatorname{arcsec} z = \arccos(1/z) $$
;;         $$ \operatorname{arccsc} z = \arcsin(1/z) $$
;;         $$ \operatorname{arccot} z = \operatorname{atan}(1/z) $$
;;
;; *   **Inverse Hyperbolic Trigonometric Functions**: `asinh`, `acosh`, `atanh`, `asech`, `acsch`, `acoth`
;;     *   These are the multi-valued inverse hyperbolic functions. The implementations return the principal values, which are also related to complex logarithms. For example:
;;         $$ \operatorname{arsinh} z = \log(z + \sqrt{z^2+1}) $$
;;         $$ \operatorname{arcosh} z = \log(z + \sqrt{z^2-1}) $$
;;         $$ \operatorname{artanh} z = \frac{1}{2} \log\left(\frac{1+z}{1-z}\right) $$
;;         $$ \operatorname{arsech} z = \operatorname{arcosh}(1/z) $$
;;         $$ \operatorname{arcsch} z = \operatorname{arsinh}(1/z) $$
;;         $$ \operatorname{arcoth} z = \operatorname{artanh}(1/z) $$
;;
;; These functions are continuous and analytic everywhere in the complex plane, except at certain points (poles) where the denominator is zero for tangent, cotangent, secant, cosecant, and at branch cuts for the inverse functions.
;;

(utls/examples-note
  (cplx/sin z1)
  (cplx/cos z1)
  (cplx/tan z1)
  (cplx/sinh z2)
  (cplx/cosh z2)
  (cplx/tanh z2)
  (cplx/asin z3)
  (cplx/acos z3)
  (cplx/atan z3))

(kind/table
 [[(gg/->image (gg/complex-function cplx/sin {:title "f(z)=sin(z)"}))
   (gg/->image (gg/complex-function cplx/cos {:title "f(z)=cos(z)"}))]
  [(gg/->image (gg/complex-function cplx/tan {:title "f(z)=tan(z)"}))
   (gg/->image (gg/complex-function cplx/sinh {:title "f(z)=sinh(z)"}))]
  [(gg/->image (gg/complex-function cplx/cosh {:title "f(z)=cosh(z)"}))
   (gg/->image (gg/complex-function cplx/tanh {:title "f(z)=tanh(z)"}))]
  [(gg/->image (gg/complex-function cplx/asin {:title "f(z)=asin(z)"}))
   (gg/->image (gg/complex-function cplx/acos {:title "f(z)=acos(z)"}))]
  [(gg/->image (gg/complex-function cplx/atan {:title "f(z)=atan(z)"}))]])
```

### Good, now write a latex formulas for `sqrt` and `sqrt1z` functions.



<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-flash-preview-04-17 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (747 4745) (4794 4816) (4818 4843) (4848 9236))) -->
<!-- End: -->
