^:kindly/hide-code
(ns transform
  (:require [fastmath.transform :as t]
            [fastmath.signal :as signal]
            
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]

            [fastmath.core :as m]
            [fastmath.random :as r]

            [fastmath.kernel.window :as window]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.polynomials :as poly]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.transform.pad :as pad]
            [fastmath.transform.wavelets :as wv]))

;; # Transforms {.unnumbered}

;; ::: {.callout-tip title="Defined functions"}
;; * `transformer`
;; * `forward-1d`, `forward-2d`
;; * `reverse-1d`, `reverse-2d`
;; :::

;; Collection of 1d and 2d transforms for signal processing, including:
;;
;; * DFT - Discrete Fourier Transform, real and complex
;; * DST/DCT - Discrete Sine and Cosine Transforms
;; * DHT - Discrete Hadamard Transform
;; * DWT/WPT - Discrete Wavelet Transform and Wavelet Packet Transform

;; ## Fast Fourier

(def fft-real (t/transformer :real :fft))

(utls/examples-note
  (seq (t/forward-1d fft-real [1 2 -10 1]))
  (seq (t/reverse-1d fft-real [-6 -12 11 -1])))

(defn magnitudes
  [xs]
  (-> (->> (t/forward-1d fft-real xs)
           (partition 2)
           (map v/mag))
      (v/div (count xs))
      (v/mult 2)))

(def sig (map (fn [id]
              (let [p (/ id 256.0)]
                (+ (m/* 12 (m/sin (* m/TWO_PI p 28)))
                   (m/* 6 (m/sin (* m/TWO_PI p 29)))
                   (m/* 8 (m/sin (* m/TWO_PI p 117)))
                   (m/* 2 (if (> (m/frac (m/* m/PI 50 p)) 0.5) -1.0 1.0))))) (range 256)))

(kind/table
 [[(gg/->image (gg/line (range) sig {:title "Signal"}))
   (gg/->image (gg/line (range) (:spectrum (signal/spectrum sig {:method :amplitude})) {:title "FFT of the signal"}))]])


;; ### Windows

;; A collection of window functions used in spectral analysis. A `fastmath.kernel/window` function allows to generate window coefficients for selected window type and optional parameters. Namespace `fastmath.kernel.window` contains implementations of all window functions.

;; Most of windows are based on continuous functions which are later sampled to get discrete coefficients.

;; `k/window` common options are:

;; * `:symmetric?` (default: `true`)
;;     * `true` (default) - window coefficients are symmetric
;;     * `false` - periodic window (last coefficient is skipped)
;; * `:normalize?` (default: `true`), for continuous window functions
;;     * `true` - maximum value is `1`
;;     * `false` - intergral of function is `1`

(utls/examples
  (k/window :welch 5)
  (k/window :welch 4 {:symmetric? false})
  (k/window :welch 5 {:normalize? false}))

^:kindly/hide-code
(defn window-plots
  ([window] (window-plots window nil))
  ([window opts]
   (partition 2 (map #(gg/->image % {:height 300}) (gg/fft-window-plots window sig opts)))))

;; Each window is illustrated by 4 charts:

;; * a window function
;; * FFT of a signal after applying a window
;; * FFT of a window (padded by zeros for interpolation)
;; * main- and side-lobes 

;; #### Rectangular

;; There are two rectangular windows: `:rectangular` and `:rectangular05`. The latter has values `0.5` at both ends.

(k/window :rectangular 8)
(kind/table (window-plots :rectangular))

(k/window :rectangular05 8)
(kind/table (window-plots :rectangular05))

;; #### Triangular

;; Triangular window has a `:shift` parameter, which shifts both ends of the triangle to a non-zero value. By default `:shift` is set to `2` (which results in the same coefficients as in Scipy).

;; $$w(n) = 1 - \left|\frac{n - \frac{N}{2}}{\frac{N+shift}{2}}\right|,\quad 0\le n \le N$$

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :triangular 5 {:shift 0})
;; :::
;; ::: {.g-col-6}
(k/window :triangular 8 {:shift 0})
;; :::
;; ::::

(kind/table (window-plots :triangular {:shift 0}))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :triangular 5)
;; :::
;; ::: {.g-col-6}
(k/window :triangular 8)
;; :::
;; ::::

(kind/table (window-plots :triangular))

;; #### Parzen

;; There are two Parzen window definitions, based on continuous function and discrete formulation ([wikipedia](https://en.wikipedia.org/wiki/Window_function#Parzen_window)) when `:discrete?` option is set to `true` (default).

;; Again, the difference is in shifted ends of the window. Continuous version has ends set to `0.0`.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen 5)
;; :::
;; ::: {.g-col-6}
(k/window :parzen 8)
;; :::
;; ::::

(kind/table (window-plots :parzen))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen 5 {:discrete? false})
;; :::
;; ::: {.g-col-6}
(k/window :parzen 8 {:discrete? false})
;; :::
;; ::::

(kind/table (window-plots :parzen {:discrete? false}))

;; 

;; #### B-spline

;; B-spline family is a generalization of above windows. `:order` parameter controls number of convolutions of rectangle window.

;; Order:

;; * `1` - rectangular window (with the first coefficient set to `0.0`)
;; * `2` - triangular window
;; * `4` - Parzen window

;; By default, `:order` is set to a value of `3`.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :b-spline 5)
;; :::
;; ::: {.g-col-6}
(k/window :b-spline 8)
;; :::
;; ::::

(kind/table (window-plots :b-spline))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :b-spline 5 {:order 5})
;; :::
;; ::: {.g-col-6}
(k/window :b-spline 8 {:order 5})
;; :::
;; ::::

(kind/table (window-plots :b-spline {:order 5}))

;; #### Welch

;; Polynomial window, inverted parabola.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :welch 5)
;; :::
;; ::: {.g-col-6}
(k/window :welch 8)
;; :::
;; ::::

(kind/table (window-plots :welch))

;; #### Connes

;; Polynomial family, square of inverted parabola. 

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :connes 5)
;; :::
;; ::: {.g-col-6}
(k/window :connes 8)
;; :::
;; ::::

(kind/table (window-plots :connes))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :connes 5 {:alpha 1.1})
;; :::
;; ::: {.g-col-6}
(k/window :connes 8 {:alpha 1.1})
;; :::
;; ::::

(kind/table (window-plots :connes {:alpha 1.1}))

;; #### Parzen algebraic

;; Algebraic family

;; $$w(x)=1-\gamma\left|2x\right|^u$$

;; Parameters:

;; * `:gamma` - $0<\gamma\le 1.0$, default: `1.0`
;; * `:u` - $u>0$, default: `3.0`

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen-algebraic 5)
;; :::
;; ::: {.g-col-6}
(k/window :parzen-algebraic 8)
;; :::
;; ::::

(kind/table (window-plots :parzen-algebraic))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen-algebraic 5 {:gamma 0.9 :u 1.5})
;; :::
;; ::: {.g-col-6}
(k/window :parzen-algebraic 8 {:gamma 0.9 :u 1.5})
;; :::
;; ::::

(kind/table (window-plots :parzen-algebraic {:gamma 0.9 :u 1.5}))

;; #### Singla and Singh

;; Family of desired order continuous polynomial time window functions, see [paper](https://www.eng.buffalo.edu/Research/code/jrnl/Window.pdf).

;; Parameter `:order`, default `1`.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :singla-singh 5)
;; :::
;; ::: {.g-col-6}
(k/window :singla-singh 8)
;; :::
;; ::::

(kind/table (window-plots :singla-singh))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :singla-singh 5 {:order 4})
;; :::
;; ::: {.g-col-6}
(k/window :singla-singh 8 {:order 4})
;; :::
;; ::::

(kind/table (window-plots :singla-singh {:order 4}))

;; #### Sinc Lobe

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :sinc 5)
;; :::
;; ::: {.g-col-6}
(k/window :sinc 8)
;; :::
;; ::::

(kind/table (window-plots :sinc))

;; #### Fejer

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :fejer 5)
;; :::
;; ::: {.g-col-6}
(k/window :fejer 8)
;; :::
;; ::::

(kind/table (window-plots :fejer))

;; #### de la Vallee Poussin

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :de-la-vallee-poussin 5)
;; :::
;; ::: {.g-col-6}
(k/window :de-la-vallee-poussin 8)
;; :::
;; ::::

(kind/table (window-plots :de-la-vallee-poussin))

;; #### Lanczos

;; Lanczos family, power of sinc function with parameter `:L` (power, default: `3.0`).

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :lanczos 5)
;; :::
;; ::: {.g-col-6}
(k/window :lanczos 8)
;; :::
;; ::::

(kind/table (window-plots :lanczos))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :lanczos 5 {:L 5})
;; :::
;; ::: {.g-col-6}
(k/window :lanczos 8 {:L 5})
;; :::
;; ::::

(kind/table (window-plots :lanczos {:L 5}))

;; #### Hamming

;; Two Hamming windows:

;; * `:hamming` - raised cosine window with rounded alpha, $\alpha=0.54$
;; * `:hamming-exact` - raised cosine window with $\alpha=\frac{25}{46}$

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :hamming 5)
;; :::
;; ::: {.g-col-6}
(k/window :hamming 8)
;; :::
;; ::::

(kind/table (window-plots :hamming))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :hamming-exact 5)
;; :::
;; ::: {.g-col-6}
(k/window :hamming-exact 8)
;; :::
;; ::::

(kind/table (window-plots :hamming-exact))

;; #### Hann

;; Raised cosine window with $\alpha=0.5$

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :hann 5)
;; :::
;; ::: {.g-col-6}
(k/window :hann 8)
;; :::
;; ::::

(kind/table (window-plots :hann))

;; #### Raised Cosine

;; Raised Cosine family with parameter `:alpha`, $1/2\leq\alpha\leq 1$

;; For `:alpha`

;; * `0.5` - Hann window (default)
;; * `0.54` - Hamming
;; * `25/46` - Hammin exact
;; * `1` - Rectangular

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :raised-cosine 5 {:alpha 0.75})
;; :::
;; ::: {.g-col-6}
(k/window :raised-cosine 8 {:alpha 0.75})
;; :::
;; ::::

(kind/table (window-plots :raised-cosine {:alpha 0.75}))

;; #### Webseter-Hamming

;; Generalized Hamming window, parameter `:v` (default: `1.0`). $v\geq 0.0$.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :webster-hamming 5)
;; :::
;; ::: {.g-col-6}
(k/window :webster-hamming 8)
;; :::
;; ::::

(kind/table (window-plots :webster-hamming))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :webster-hamming 5 {:v 0.5})
;; :::
;; ::: {.g-col-6}
(k/window :webster-hamming 8 {:v 0.5})
;; :::
;; ::::

(kind/table (window-plots :webster-hamming {:v 0.5}))

;; #### Power of Cosine

;; Cosine lobe raised to a power of parameter `:m`, default `1.0`.

;; For `:m`:

;; * `0.0` - Rectangular window
;; * `1.0` - Cosine lobe (default)
;; * `2.0` - Hann window

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :power-of-cosine 5)
;; :::
;; ::: {.g-col-6}
(k/window :power-of-cosine 8)
;; :::
;; ::::

(kind/table (window-plots :power-of-cosine))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :power-of-cosine 5 {:m 4.0})
;; :::
;; ::: {.g-col-6}
(k/window :power-of-cosine 8 {:m 4.0})
;; :::
;; ::::

(kind/table (window-plots :power-of-cosine {:m 4.0}))

;; #### Raised Power of Cosine

;; Combinatin of raised cosine and power of cosine windows. Two parameters `:alpha` (default: `0.05`) and `:m` (power, default `1.0`)

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :raised-power-of-cosine 5)
;; :::
;; ::: {.g-col-6}
(k/window :raised-power-of-cosine 8)
;; :::
;; ::::

(kind/table (window-plots :raised-power-of-cosine))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :raised-power-of-cosine 5 {:m 4.0 :alpha 0.1})
;; :::
;; ::: {.g-col-6}
(k/window :raised-power-of-cosine 8 {:m 4.0 :alpha 0.1})
;; :::
;; ::::

(kind/table (window-plots :raised-power-of-cosine {:m 4.0 :alpha 0.1}))

;; #### Parzen Cosine

;; Parzen Cosine family. Parameters `:gamma` (default: `1.0`) and `:m` (default `2.0`)

;; $$w(x) = 1+\cos(\pi\gamma\left|2x\right|^m)$$

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen-cosine 5)
;; :::
;; ::: {.g-col-6}
(k/window :parzen-cosine 8)
;; :::
;; ::::

(kind/table (window-plots :parzen-cosine))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :parzen-cosine 5 {:m 4.0 :gamma 0.9})
;; :::
;; ::: {.g-col-6}
(k/window :parzen-cosine 8 {:m 4.0 :gamma 0.9})
;; :::
;; ::::

(kind/table (window-plots :parzen-cosine {:m 4.0 :gamma 0.9}))

;; #### Bohman

;; Cosine lobe convolved with itself.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :bohman 5)
;; :::
;; ::: {.g-col-6}
(k/window :bohman 8)
;; :::
;; ::::

(kind/table (window-plots :bohman))

;; #### Trapezoid

;; Trapezoid window. Combination of Triangular and Rectangular windows.

;; Parameter `:alpha` controls flatness of trapezoid (default: `0.25`).

;; For `:alpha`:

;; * `0.0` - Triangular window
;; * `0.5` - Rectangular window

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :trapezoid 5)
;; :::
;; ::: {.g-col-6}
(k/window :trapezoid 8)
;; :::
;; ::::

(kind/table (window-plots :trapezoid))

;; #### Tukey

;; Tukey window family, combination of Rectangular and Hann windows.

;; For `:alpha`:

;; * `0.0` - Rectangular window
;; * `1.0` - Hann window

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :tukey 5)
;; :::
;; ::: {.g-col-6}
(k/window :tukey 8)
;; :::
;; ::::

(kind/table (window-plots :tukey))

;; #### Bartlett-Hann

;; Combination of Triangular and Hann windows.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :bartlett-hann 5)
;; :::
;; ::: {.g-col-6}
(k/window :bartlett-hann 8)
;; :::
;; ::::

(kind/table (window-plots :bartlett-hann))

;; #### Blackman-Harris, Nutall

;; A collection of Blackman-Harris family, a weighted sum of cosines.

;; $$w(x)=\sum_{l=0}^{L-1}\alpha_l\cos(2\pi l x)$$

;; A `:blackman-harris-family` accepts `:coeffs` parameter as a list of coefficients.

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris-family 5 {:coeffs [0.4 0.5 0.1]})
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris-family 8 {:coeffs [0.4 0.5 0.1]})
;; :::
;; ::::

(kind/table (window-plots :blackman-harris-family {:coeffs [0.4 0.5 0.1]}))

;; There is a predefined collection of the following windows from the family:

;; * `:blackman`
;; * `:blackman-exact`
;; * `:blackman-harris` - three-term, minimum side-lobe
;; * `:blackman-harris-61db` - -61dB, three-term
;; * `:blackman-harris-67db` - -67dB, three-term
;; * `:blackman-harris-74db` - -74dB, four-term
;; * `:blackman-harris-92db` - -92dB, four-term
;; * `:nutall-3-1st` - three-term Nutall, continuous 1st derivative
;; * `:nutall-3-3rd` - three-term Nutall, continuous 3rd derivative
;; * `:blackman-nutall` - four-term, minimum side-lobe
;; * `:nutall-1st` - four-term Nutall, continuous 1st derivative
;; * `:nutall-3rd` - four-term Nutall, continuous 3rd derivative
;; * `:nutall-5th` - four-term Nutall, continuous 5th derivative
;; * `:mottaghi-kashtiban-shayesteh`  - four-term

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman 8)
;; :::
;; ::::

(kind/table (window-plots :blackman))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-exact 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-exact 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-exact))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-harris))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris-61db 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris-61db 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-harris-61db))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris-67db 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris-67db 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-harris-67db))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris-74db 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris-74db 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-harris-74db))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-harris-92db 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-harris-92db 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-harris-92db))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :nutall-3-1st 5)
;; :::
;; ::: {.g-col-6}
(k/window :nutall-3-1st 8)
;; :::
;; ::::

(kind/table (window-plots :nutall-3-1st))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :nutall-3-3rd 5)
;; :::
;; ::: {.g-col-6}
(k/window :nutall-3-3rd 8)
;; :::
;; ::::

(kind/table (window-plots :nutall-3-3rd))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :blackman-nutall 5)
;; :::
;; ::: {.g-col-6}
(k/window :blackman-nutall 8)
;; :::
;; ::::

(kind/table (window-plots :blackman-nutall))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :nutall-1st 5)
;; :::
;; ::: {.g-col-6}
(k/window :nutall-1st 8)
;; :::
;; ::::

(kind/table (window-plots :nutall-1st))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :nutall-3rd 5)
;; :::
;; ::: {.g-col-6}
(k/window :nutall-3rd 8)
;; :::
;; ::::

(kind/table (window-plots :nutall-3rd))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :nutall-5th 5)
;; :::
;; ::: {.g-col-6}
(k/window :nutall-5th 8)
;; :::
;; ::::

(kind/table (window-plots :nutall-5th))

;; :::: {.grid}
;; ::: {.g-col-6}
(k/window :mottaghi-kashtiban-shayesteh 5)
;; :::
;; ::: {.g-col-6}
(k/window :mottaghi-kashtiban-shayesteh 8)
;; :::
;; ::::

(kind/table (window-plots :mottaghi-kashtiban-shayesteh))

;;

;; ## Sine / Cosine

;; ## Hadamard

;; ## Wavelets

;; ::: {.callout-tip title="Defined functions"}
;; * `wavelet`
;; * `wavelet-name`
;; * `coeff-size`
;; * `phi`, `psi`
;; * `scaling-function`, `wavelet-function`
;; :::

;; ### DWT

;; The Discrete Wavelet Transform (DWT) provides a time–frequency (or time–scale) representation of a signal, analogous to the Short-Time Fourier Transform (STFT). Unlike the STFT, however, the DWT offers multi-resolution analysis, meaning that time and frequency resolutions vary across the spectrum:

;; * At low frequencies, the DWT provides higher frequency resolution and lower time resolution.
;; * At high frequencies, it provides higher time resolution and lower frequency resolution.

;; In the DWT, the signal is successively passed through a pair of filters: a low-pass filter, producing the approximation coefficients, and a high-pass filter, producing the detail coefficients. After filtering, each output is downsampled by a factor of two. The approximation coefficients are then recursively decomposed in the same manner, yielding multiple levels of detail.

;; The resulting representation consists of detail coefficients at different decomposition levels, each corresponding to a progressively narrower frequency band. At each level, the number of coefficients is halved, while the associated frequency range becomes smaller, collectively capturing both coarse and fine signal structures across scales.

;; DWT works only on a singal which length is the power of 2.

;; DWT can be applied to desired level with a `{:level n}` option. Maximum number of levels can be calculated by applying $log_2()$ to the length of the signal. For example signal of the length of 8192 will contain 13 levels. 

(m/round (m/log2 8192))

;; Layout of the coefficients after full DWT goes as follows:

;; * Id 0 - zero frequency term (approximation)
;; * Id 1 - last level N (detail), lowest time resolution, lowest frequencies, narrowest band
;; * Ids 2 and 3 - level N-1 (detail)
;; * ids from 4 to 7 - level N-2 (detail)
;; * ...
;; * ids from $2^{N-1}$ to $2^N-1$ - Level 1 (detail), highest time resolution, highest frequencies, widest band (half of the spectrum)

;; Let's see how signal is decomposed into lower and higher frequencies using different wavelets. First define three wavelet transforms:

(def haar (t/transformer :dwt "haar"))
(def daubechies-2 (t/transformer :dwt "db2"))
(def coiflet-1 (t/transformer :dwt :coiflet-1))
(def biorthogonal-39 (t/transformer :dwt "bior3.9"))

(utls/examples-note
  (seq (t/forward-1d haar (range 8)))
  (seq (t/forward-1d daubechies-2 (range 8)))
  (seq (t/forward-1d coiflet-1 (range 8)))
  (seq (t/forward-1d biorthogonal-39 (range 8))))

;; First we construct a linear chrip signal with added low frequency square wave to the second half of a chirp.

(defn chirp+sin
  [^long levels]
  (let [size (m/exp2 levels)
        s4 (m// size 4.0)
        s2 (m// size 2.0)
        square (m/- levels 4)]
    (map (fn [^long x] (m/+ (m/sin (m/* m/TWO_PI (m/cb (m// x s4))))
                           (if (m/> x s2)
                             (if (m/even? (m/>> x square)) -0.75 0.75)
                             0.0))) (range size))))

(def chirp+sin11 (chirp+sin 11))

;; Length of the signal

(count chirp+sin11)

(gg/->image (gg/line (range) chirp+sin11 {:title "Chirp signal"}))

;; Scaleogram shows time-frequency coefficients (L1, normalized to be comparable) for each level. You can see that first level is most detailed but covers half of the band, while last level contains only one coefficient and represents the lowest frequency. Level 0 (DC) is not shown.

;; You can see the pattern of orange stripes representing raising frequencies over the time. Also, low frequency is recognized at level 8 (especially when HAAR and BiOrthogonal 3.9 wavelets are used). Dark region at top-left means thare are no higher frequencies at the beginning of the signal.

(kind/table
 [[(gg/->image (gg/dwt-scaleogram (t/forward-1d haar chirp+sin11) {:wavelet "HAAR" :normalize? true}))
   (gg/->image (gg/dwt-scaleogram (t/forward-1d daubechies-2 chirp+sin11) {:wavelet "Daubechies 2" :normalize? true}))]
  [(gg/->image (gg/dwt-scaleogram (t/forward-1d coiflet-1 chirp+sin11) {:wavelet "Coiflet-1" :normalize? true}))
   (gg/->image (gg/dwt-scaleogram (t/forward-1d biorthogonal-39 chirp+sin11) {:wavelet "BiOrthogonal-39" :normalize? true}))]])

;; Each level (band) can be reconstructed and analyzed separately. To do this it's enough to set all other bands to zero value and reverse DWT. Lowest frequencies reveal the shape of the wavelet.

(defn dwt-band
  "Return detail band at given level"
  [coeffs level]
  (let [N (count coeffs)
        target (double-array N)
        rlevel (int (m/exp2 (m/- (m/round (m/log2 N)) level)))]
    (System/arraycopy coeffs rlevel target rlevel rlevel) ;; copy given level
    target))

;; Daubechies 2

(kind/table
 (let [c (t/forward-1d daubechies-2 chirp+sin11)]
   (partition 3 (map #(gg/->image (gg/line (range) (t/reverse-1d daubechies-2 (dwt-band c %))
                                           {:title (str "Level: " %)})) (range 1 10)))))

;; HAAR

(kind/table
 (let [c (t/forward-1d haar chirp+sin11)]
   (partition 3 (map #(gg/->image (gg/line (range) (t/reverse-1d haar (dwt-band c %))
                                           {:title (str "Level: " %)})) (range 1 10)))))


;; ### WPT

;; Wavelet Packet Transform (or Decomposition) applies filters also on high frequencies (on details). This way the whole band is splitted evenly into 2^N bands at level N. There is no reason to perform full transform. For the best result it's better to stop early. Let's see how it works.

(def haar-wpt (t/transformer :wpt :haar))
(def coiflet-1-wpt (t/transformer :wpt :coiflet-1))

(defn wpt-bands
  "Return all bands at given level"
  [coeffs level]
  (let [N (count coeffs)
        cnt (m/round (m/exp2 level))
        len (m// N cnt)]
    (map (fn [band]
           (let [target (double-array N)
                 pos (m/* band len)]
             (System/arraycopy coeffs pos target pos len)
             target)) (range cnt))))

;; HAAR WPT at level 2 (4 bands)

(kind/table
 [(map-indexed #(gg/->image (gg/line (range) (t/reverse-1d haar-wpt %2) {:title (str "Band " (inc %1))})) (wpt-bands (t/forward-1d haar-wpt chirp+sin11) 2))])

;; Coiflet-1 WPT at level 4 (16 bands)

(kind/table
 (partition 4 (map-indexed #(gg/->image (gg/line (range) (t/reverse-1d coiflet-1-wpt %2) {:title (str "Band " (inc %1))})) (wpt-bands (t/forward-1d coiflet-1-wpt chirp+sin11) 4))))

;; ### Decomposed

;; ### List of wavelets

;; There are two sources of wavelets, one based on [JWave](https://github.com/graetz23/JWave) and second based on Matlab.

;; ::: {.callout-caution}
;; JWave wavelets have some issues with coefficients: some of them are not reversible or have reversed coefficients.
;; :::

;; `keyword` based names are for `JWave` wavelets, `string` based are for Matlab (built-in).

;; #### Haar

;; Haar (step, or Daubechies 1, db1) wavelet. `:haar-orthogonal` is not normalized version of a wavelet.

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["haar" "db1"]
                 [:haar :haar-orthogonal]]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "haar"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair :haar-orthogonal))])

;; #### Daubechies

;; Daubechies, extremal phase,  wavelets.

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[(mapv #(str "db" %)(range 2 46))
                 (mapv #(keyword (str "daubechies-" %)) (range 2 21))]]})
(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "db2"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "db3"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "db10"))])

;; #### Symlet

;; Symlet, least asymmetric, wavelets.

;; Symlet 2 and 3 are the same as Daubechies. Wavelets starting with `la` are defined in `wavelets` R package and are the same as symlets ("sym4"="la8", "sym5"="la10", etc.).

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[(conj (mapv #(str "sym" %)(range 2 46))
                       "la8" "la10" "la12" "la14" "la16" "la18" "la20")
                 (mapv #(keyword (str "symlet-" %)) (range 2 21))]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "sym2"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "sym4"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "sym10"))])

;; #### Coiflet

;; Coiflet wavelets.

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[(mapv #(str "coif" %)(range 1 6))
                 (mapv #(keyword (str "coiflet-" %)) (range 1 6))]]})
(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "coif1"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "coif2"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "coif4"))])

;; #### Best-localized Daubechies

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[(mapv #(str "bl" %) [7 9 10]) '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "bl7"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "bl9"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "bl10"))])

;; #### Fejér-Korovkin

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[(mapv #(str "fk" %) [4 6 8 14 18 22]) '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "fk4"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "fk8"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "fk18"))])

;; #### Morris minimum-bandwidth

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["mb4.2" "mb8.2" "mb8.3" "mb8.4" "mb10.3" "mb12.3" "mb14.3" "mb16.3" "mb18.3" "mb24.3" "mb32.3"] '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "mb8.2"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "mb8.3"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "mb24.3"))])

;; #### Han linear-phase moments

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["han2.3" "han3.3" "han4.5" "han5.5"] '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "han2.3"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "han3.3"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "han4.5"))])

;; #### Beylkin

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["beyl"] '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "beyl"))])

;; #### Discrete Meyer

;; There is a difference in number of coefficients.

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["dmey"] [:discrete-mayer]]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "dmey"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair :discrete-mayer))])

;; #### Vaidyanathan

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["vaid"] '-]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "vaid"))])

;; #### Bi-orthogonal

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [[["bior1.1" "bior1.3" "bior1.5"
                  "bior2.2" "bior2.4" "bior2.6" "bior2.8"
                  "bior3.1" "bior3.3" "bior3.5" "bior3.7" "bior3.9"
                  "bior4.4" "bior5.5" "bior6.8"
                  "rbio1.1" "rbio1.3" "rbio1.5"
                  "rbio2.2" "rbio2.4" "rbio2.6" "rbio2.8"
                  "rbio3.1" "rbio3.3" "rbio3.5" "rbio3.7" "rbio3.9"
                  "rbio4.4" "rbio5.5" "rbio6.8"]
                 [:biorthogonal-11 :biorthogonal-13 :biorthogonal-15 :biorthogonal-22
                  :biorthogonal-24 :biorthogonal-26 :biorthogonal-28 :biorthogonal-31
                  :biorthogonal-33 :biorthogonal-35 :biorthogonal-37 :biorthogonal-39
                  :biorthogonal-44 :biorthogonal-55 :biorthogonal-68]]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair "bior1.3"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "bior1.3" {:reconstruction? true}))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "bior4.4"))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair "bior4.4" {:reconstruction? true}))])


;; #### Other

;; All below (but `:legendre-1`) are not orthogonal (nor bi-orthogonal) and probably are defined with wrong coefficients. Listing them for the reference only.

(kind/table
 {:column-names ["Built-in" "JWave"]
  :row-vectors [['- [:legendre-1 :legendre-2 :legendre-3
                     :cdf-53 :cdf-97 :battle-23]]]})

(kind/table
 [(map #(gg/->image % {:height 200}) (gg/wavelet-pair :legendre-2))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair :cdf-53))
  (map #(gg/->image % {:height 200}) (gg/wavelet-pair :battle-23))])

;; ## Padding

(def partial-signal (m/sample m/sin 0.25 2.5 20))

(defn padding-chart
  [signal pad-method]
  (let [len (count signal)
        start (pad/pad-position :both 90 len)]
    (gg/scatters [[(name pad-method) (range) (t/pad signal 90 pad-method)]
                  ["orignal" (range start (+ start len)) signal]]
                 {:legend-name "Padding" :size 3 :palette [(last gg/palette-blue-0) gg/color-main]
                  :title (str "Padding: " (name pad-method))})))

(kind/table
 (partition 2 (map #(gg/->image (padding-chart partial-signal %)) [:zero :edge
                                                                   :linear :periodic 
                                                                   :symmetric :antisymmetric
                                                                   :reflect :antireflect])))




;; DWT HAAR decomposition for zero and perioding padding (4x original size).

(kind/table
 [[(gg/->image (gg/dwt-scaleogram (t/forward-1d haar (t/pad chirp+sin11 4096 :zero)) {:normalize? true :wavelet "HAAR, zero padding"}))
   (gg/->image (gg/dwt-scaleogram (t/forward-1d haar (t/pad chirp+sin11 4096 :antireflect)) {:normalize? true :wavelet "HAAR, anti-reflect padding"}))]])



;; ----

(gg/->image (gg/functions [["basic" m/sin]
                           ["noisy" (fn [x] (+ (m/sin x) (* 0.2 (- (rand) 0.5))))]]
                          {:x [m/-TWO_PI m/TWO_PI]
                           :ylim [-2 2]
                           :steps 500
                           :palette gg/palette-blue-0}))

(gg/->image (gg/function2d (fn [[x y]] (m/sin (m/* x (m/cos y)))) {:x [m/-TWO_PI m/TWO_PI]
                                                                  :y [m/-TWO_PI m/TWO_PI]
                                                                  :title "sin(x*cos(y))"
                                                                  :legend-name "value"}))

(gg/function2d (fn [[x y]] (m/sin (m/* x (m/cos y)))) {:x [m/-TWO_PI m/TWO_PI]
                                                      :y [m/-TWO_PI m/TWO_PI]
                                                      :title "sin(x*cos(y))"
                                                      :legend-name "value"})

(gg/function m/tan {:x [m/-TWO_PI m/TWO_PI]
                    :title "tan(x)"
                    :ylab "y=tan(x)"
                    :ylim [-2 2] ;; we need to limit y axis
                    :steps 500})

(let [xs (repeatedly 2000 r/grand)
      ys (map (fn [x] (+ (r/grand (+ 0.1 (* x 0.5))) (m/sin (* 2 x)))) xs)]
  (gg/scatter xs ys {:title "Scatter"}))

(let [xy (take 1000 (r/sequence-generator :r2 2))]
  (gg/scatter xy nil {:title "R2 low-discrepancy sequence generator"}))

(gg/functions [["tan" m/tan]
               ["cot" m/cot]
               ["sin" m/sin]
               ["cos" m/cos]]
              {:x [m/-TWO_PI m/TWO_PI]
               :title "Basic trigonometric functions"
               :ylim [-2 2]
               :steps 500
               :palette gg/palette-blue-1})

(gg/->image (gg/function m/sec
                         {:x [m/-TWO_PI m/TWO_PI]
                          :ylim [-2 2]
                          :steps 500}))


;; ## FFT

;; Details about FFT and use-cases

;; Some examples:



;; ## Wavelets

;; ## Compression and denoising

;; An use case with charts

(def domain (m/slice-range 0 10 512))
(def signal (map (fn [x] (+ (Math/sin x)
                         (* 0.1 (- (rand) 0.5)))) ;; add some noise
               domain))
(def denoised-signal (t/denoise fft-real signal {:method :hard}))

^:kind/table
[[(gg/->image (gg/line domain signal {:title "Original signal"}))
  (gg/->image (gg/line domain denoised-signal {:title "Denoised signal"}))]]

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.transform)
