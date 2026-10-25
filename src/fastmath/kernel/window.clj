(ns fastmath.kernel.window
  "Tapering window functions for spectral analysis, filter design and signal processing.

  Every window is a function `(window N)` or `(window N options)` which returns a sequence of `N` doubles; `(window)` returns `256` samples. `N` must be an integer valued number `>= 0`; anything else throws `ex-info` (`:type :invalid-input`, `:function`, `:property`, `:value`, `:expected`). `N = 0` returns `()` and `N = 1` returns `(1.0)` for every window (with the default `:normalize?`). Most windows are registered by name in `fastmath.kernel/window`, which also handles the periodic variant.

  Conventions:

  - Continuous model windows (most of them) are a continuous function `f` on `[-0.5, 0.5]` with the integral `1` (the `-continuous` functions and the constructors `->name`), sampled at `N` evenly spaced points including both ends. The window functions have support `[-0.5, 0.5]`; the kernels of `fastmath.kernel.density` are the standard forms with support `[-1, 1]`, so a window is the kernel shape at `2x` (for example `welch` is the Epanechnikov shape).
  - Discrete windows (`triangular` and `parzen` default paths, `dolph-chebyshev`, `ultraspherical`, `saramaki`, `legendre`, `shayesteh-kashtiban`, `kaiser-bessel-derived`, `vorbis`, `cosine`) are defined directly on the index.
  - Raw scale (`:normalize? false`): continuous windows have the integral `1`, the discrete windows of Doerry (`dolph-chebyshev`, `ultraspherical`, `saramaki`, `legendre`, `shayesteh-kashtiban`) have the mean `1`; `kaiser-bessel-derived` (power complementary), `vorbis` and `cosine` (formula values) keep their own scale and have the default `:normalize? false`, as `rectangular05`.
  - Option `:normalize?`: `true` divides by the value of the window function at the middle for the continuous model windows (for even `N` the maximum is below `1`, as in scipy and Octave) and by the maximum (`:LInf`) for the discrete windows; `false` is the raw scale; `:LInf` (alias `:Linf`), `:L1`, `:L2`, `:N` and a finite non-zero number are accepted, see [[normalize-coefficients]]. `nil` or a missing key selects the default of the window, any other value throws.
  - The parameters of the windows (`:alpha`, `:order`, ...) are validated against their documented domains; `nil` is a missing parameter and a value outside of the domain throws `ex-info`. No window returns `NaN` for an accepted parameter.
  - The option `:symmetric?` is ignored by the functions of this namespace, they always return the symmetric window. Only `fastmath.kernel/window` handles it: a periodic window of `N` samples is the symmetric window of `N + 1` samples without the last one.
  - `:blackman-harris` is the Nuttall three-term minimum sidelobe window (`0.4243801 0.4973406 0.0782793`, peak sidelobe -71.5 dB), not the window called blackmanharris by scipy, Matlab and Octave, which is `blackman-harris-92db`. `blackman-nutall` is scipy `nuttall` and `nutall-1st` is Octave `nuttallwin`. The names are spelled `nutall-*`.
  - `lanczos` has the default `:L 3.0`; the usual Lanczos window of scipy is `:L 1.0`. `taylor` is sampled at `-0.5, ..., 0.5` including both ends, which differs from the grid of scipy `taylor` and Octave `taylorwin`.
  - The `:level` of `dolph-chebyshev` is the achieved peak sidelobe level. For `ultraspherical`, `saramaki`, `legendre` and `taylor` it is Doerry's design parameter and the achieved peak sidelobe level differs from it by up to about 17 dB (higher for `ultraspherical`, `saramaki` and `taylor` with a small `:n`, about 3 dB lower for `legendre`); the measured levels are in the docstrings of these windows.

  Functions:

  - helpers: [[normalize-coefficients]], [[sample-window]], [[continuous->window-fn]], [[discrete->window-fn]], [[idft]]
  - rectangular and triangular: [[rectangular]], [[rectangular05]], [[triangular]], [[trapezoid]], [[tukey]], [[welch]], [[parzen]], [[b-spline]], [[parzen-algebraic]], [[singla-singh]], [[connes]]
  - sinc family: [[sinc]], [[fejer]], [[de-la-vallee-poussin]], [[lanczos]]
  - cosine family: [[cosine]], [[hann]], [[hamming]], [[hamming-exact]], [[raised-cosine]], [[webster-hamming]], [[power-of-cosine]], [[raised-power-of-cosine]], [[parzen-cosine]], [[bohman]], [[bartlett-hann]], [[vorbis]]
  - sum of cosines: [[blackman]], [[blackman-exact]], [[blackman-harris]], [[blackman-harris-61db]], [[blackman-harris-67db]], [[blackman-harris-74db]], [[blackman-harris-92db]], [[blackman-nutall]], [[nutall-3-1st]], [[nutall-3-3rd]], [[nutall-1st]], [[nutall-3rd]], [[nutall-5th]], [[mottaghi-kashtiban-shayesteh]], [[blackman-harris-family]] (user coefficients, scipy general_cosine)
  - flat-top (amplitude accuracy): [[flat-top]] (scipy), [[flat-top-3]], [[flat-top-hp]], [[flat-top-srs]], [[flat-top-bk]] (Octave flattopwin), [[hft70]], [[hft90d]], [[hft95]], [[hft116d]], [[hft144d]], [[hft169d]], [[hft196d]], [[hft223d]], [[hft248d]]
  - Kaiser-Bessel family and DPSS approximations: [[kaiser]], [[kaiser-bessel]], [[kaiser-bessel-derived]], [[cosh]], [[avci-nacaroglu]], [[knab]], [[bessel-I1]]
  - exponential family: [[gaussian]], [[general-gaussian]], [[exponential]], [[hanning-poisson]], [[cauchy]], [[parzen-exponential]], [[parzen-geometric]]
  - sidelobe level design: [[dolph-chebyshev]], [[ultraspherical]], [[saramaki]], [[legendre]], [[taylor]], [[shayesteh-kashtiban]]
  - every continuous window has `name-continuous` (fixed windows) or `->name` (parametrized windows) returning the continuous function, which can be used with [[sample-window]].

  scipy.signal.windows and Octave signal equivalents (default options): boxcar and rectwin are `rectangular`; triang is `triangular`; bartlett is `triangular` with `:shift 0`; parzen and parzenwin are `parzen`; bohman and bohmanwin, barthann and barthannwin, blackman, hamming, hann, tukey and tukeywin have the same names; cosine is `cosine`; lanczos is `lanczos` with `:L 1.0`; general_cosine is `blackman-harris-family`; general_hamming is `raised-cosine`; blackmanharris is `blackman-harris-92db`; nuttall is `blackman-nutall`; nuttallwin is `nutall-1st`; flattop is `flat-top`; flattopwin of Octave is `flat-top-bk`; gaussian and gausswin are `gaussian` (`:alpha` is `(N - 1) / (2 std)` for scipy and `a` for Octave); general_gaussian is `general-gaussian`; exponential is `exponential` (`:alpha` is `(N - 1) / (2 tau)`); kaiser is `kaiser`; kaiser_bessel_derived is `kaiser-bessel-derived` (`:alpha` is `beta / pi`); chebwin is `dolph-chebyshev` (`:level` is `-at`); taylor and taylorwin are `taylor` (different sampling grid).

  Literature: A. W. Doerry, Catalog of Window Taper Functions for Sidelobe Control, SAND2017-4042; F. J. Harris, On the use of windows for harmonic analysis with the discrete Fourier transform, Proc. IEEE 66(1), 1978; A. H. Nuttall, Some windows with very good sidelobe behavior, IEEE Trans. ASSP 29(1), 1981; G. Heinzel, A. Ruediger and R. Schilling, Spectrum and spectral density estimation by the Discrete Fourier Transform (DFT), including a comprehensive list of window functions and some new flat-top windows, 2002."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.quadrature :as qint]
            [fastmath.polynomials :as poly]
            [fastmath.special :as special]
            [fastmath.special.gamma :as gamma])
  (:import [fastmath.java Array]))

;; https://www.researchgate.net/profile/Armin_Doerry/publication/316281181_Catalog_of_Window_Taper_Functions_for_Sidelobe_Control/links/58f92cb2a6fdccb121c9d54d/Catalog-of-Window-Taper-Functions-for-Sidelobe-Control.pdf

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmacro ^:private clamp [x body] `(if (m/<= -0.5 ~x 0.5) (m/max 0.0 (double ~body)) 0.0))
(defmacro ^:private clamp- [x body] `(if (m/<= -0.5 ~x 0.5) (double ~body) 0.0))

;; input validation helpers

(defn- invalid-input!
  "Throws `ex-info` for an invalid `property` of the window `function`."
  [function property message value expected]
  (throw (ex-info (str "window " (name function) ", " (name property) ": " message)
                  {:type :invalid-input
                   :function function
                   :property property
                   :value value
                   :expected expected})))

(defn- check-n!
  "Returns `N` as a long or throws `ex-info` unless `N` is an integer valued number `>= 0`."
  ^long [function N]
  (if (number? N)
    (let [d (double N)]
      (if (and (m/valid-double? d)
               (m/== d (m/floor d))
               (m/<= 0.0 d 2.147483647E9))
        (long d)
        (invalid-input! function :N "must be an integer >= 0" N "integer >= 0")))
    (invalid-input! function :N "must be an integer >= 0" N "integer >= 0")))

(def ^:private open-positive [0.0 false ##Inf false])
(def ^:private non-negative [0.0 true ##Inf false])

(defn- interval-str
  [[lo lo-closed? hi hi-closed?]]
  (str (if lo-closed? "[" "(") lo ", " hi (if hi-closed? "]" ")")))

(defn- in-interval?
  [^double x interval]
  (let [[lo lo-closed? hi hi-closed?] interval
        lo (double lo)
        hi (double hi)]
    (boolean (and (if lo-closed? (m/>= x lo) (m/> x lo))
                  (if hi-closed? (m/<= x hi) (m/< x hi))))))

(defn- param
  "Returns the number under `property` in `opts` as a double (`default` when missing or `nil`); throws unless it is in `interval` `[lo lo-closed? hi hi-closed?]`."
  [function opts property default interval]
  (let [v (get opts property)
        v (if (nil? v) default v)]
    (if (and (number? v) (in-interval? (double v) interval))
      (double v)
      (invalid-input! function property
                      (str "must be a number in " (interval-str interval))
                      v (str "number in " (interval-str interval))))))

(defn- int-param
  "Like [[param]] but the number must have an integral value; returns a long."
  [function opts property default interval]
  (let [v (get opts property)
        v (if (nil? v) default v)]
    (if (and (number? v)
             (let [d (double v)]
               (and (m/valid-double? d) (m/== d (m/floor d)) (in-interval? d interval))))
      (long (double v))
      (invalid-input! function property
                      (str "must be an integer in " (interval-str interval))
                      v (str "integer in " (interval-str interval))))))

(defn- boolean-param
  "Returns the boolean under `property` in `opts` (`default` when missing or `nil`), throws for other values."
  [function opts property default]
  (let [v (get opts property)]
    (cond
      (nil? v) default
      (boolean? v) v
      :else (invalid-input! function property "must be true or false" v "true or false"))))

;; normalization

(defn- canonical-normalize
  "Validates `normalize?` and returns `true`, `false`, `:LInf`, `:L1`, `:L2`, `:N` or a double; throws otherwise."
  [function normalize?]
  (cond
    (true? normalize?) true
    (false? normalize?) false
    (#{:LInf :Linf} normalize?) :LInf
    (#{:L1 :L2 :N} normalize?) normalize?
    (and (number? normalize?)
         (let [d (double normalize?)] (and (m/valid-double? d) (not (m/zero? d))))) (double normalize?)
    :else (invalid-input! function :normalize?
                          "must be true, false, :LInf, :L1, :L2, :N or a finite non-zero number"
                          normalize?
                          "true, false, :LInf, :Linf, :L1, :L2, :N or a finite non-zero number")))

(defn- window-mode
  "Canonical `:normalize?` option of `opts`; `nil` or missing means `default`."
  [function opts default]
  (let [n (get opts :normalize?)]
    (canonical-normalize function (if (nil? n) default n))))

(defn- apply-normalize
  "Applies the canonical (not `true`) `mode` to `coeffs`; empty input is returned unchanged."
  [coeffs mode]
  (cond
    (false? mode) coeffs
    (empty? coeffs) coeffs
    (number? mode) (v/div coeffs mode)
    :else (case mode
            :LInf (v/normalize-LInf coeffs)
            :L1 (v/normalize-L1 coeffs)
            :L2 (v/normalize coeffs)
            :N (v/mult (v/normalize-L1 coeffs) (count coeffs)))))

(defn normalize-coefficients
  "Normalizes window coefficients.

  Parameters:

  - `coeffs` - a sequence (or vector) of numbers
  - `normalize?` - one of:
      - `true` or `:LInf` (alias `:Linf`) - divides by the maximum absolute value, the maximum becomes `1.0`
      - `:L1` - divides by the sum of the absolute values
      - `:L2` - divides by the euclidean length
      - `:N` - the `:L1` result multiplied by the number of coefficients
      - a finite non-zero number - divides all coefficients by this number
      - `false` - returns `coeffs` unchanged

  Returns the normalized coefficients (the same container type as `coeffs`). Empty input is returned unchanged. All-zero input stays all-zero for the norm based modes.

  Throws `ex-info` for `nil`, a zero, `NaN` or infinite number and for any other value of `normalize?`.

  The windows of this namespace resolve `true` before they call this function: for continuous model windows `true` means division by the window function value at `0` (see [[sample-window]]), the default `nil` is replaced by the default of the window.

  See also [[sample-window]]."
  [coeffs normalize?]
  (let [mode (canonical-normalize :normalize-coefficients normalize?)]
    (apply-normalize coeffs (if (true? mode) :LInf mode))))

(defn- sample*
  "Samples continuous `f` at `N` (long) points of [-0.5, 0.5] and normalizes with canonical `mode`."
  [function f ^long N mode]
  (if (m/zero? N)
    ()
    (let [coeffs (m/sample f -0.5 0.5 N)]
      (if (true? mode)
        (let [f0 (double (f 0.0))]
          (if (and (m/valid-double? f0) (not (m/zero? f0)))
            (if (m/== N 1) (list 1.0) (v/div coeffs f0))
            (invalid-input! function :normalize? "the value of the window function at 0 is zero or not finite, can't divide by it" f0 "window function with a finite non-zero value at 0")))
        (apply-normalize coeffs mode)))))

(defn- continuous-window
  "Samples the continuous window function `f` of the window `function`, `N` is validated."
  [function f N mode]
  (sample* function f (check-n! function N) mode))

(defn- discrete-window
  "Window defined on the index: `make` returns the coefficients for `N >= 2`, `true-mode` replaces `true` (`:LInf` or `false`)."
  [function make N mode true-mode]
  (let [N (check-n! function N)
        mode (if (true? mode) true-mode mode)]
    (cond
      (m/zero? N) ()
      (m/== N 1) (apply-normalize (list 1.0) mode)
      :else (apply-normalize (make N) mode))))

(defn sample-window
  "Samples a continuous window function and normalizes the samples.

  Parameters:

  - `f` - a function of one double, defined on `[-0.5, 0.5]`
  - `N` - number of samples, an integer `>= 0`
  - `normalize?` (optional, default `true`) - `true` divides the samples by `(f 0.0)` (the value at the middle, for even `N` the maximum of the samples can be lower than `1.0`), `false` returns the samples as they are; the other accepted values are `:LInf` (alias `:Linf`), `:L1`, `:L2`, `:N` and a finite non-zero number, see [[normalize-coefficients]]. `nil` means `true`.

  Returns a sequence of `N` doubles at the points `-0.5, ..., 0.5` (both ends included; for `N = 1` the point `0.0`), `()` for `N = 0`.

  Throws `ex-info` for an invalid `N`, an invalid `normalize?` or when `normalize?` is `true` and `(f 0.0)` is zero or not finite.

  See also [[normalize-coefficients]], [[continuous->window-fn]]."
  ([f N] (sample-window f N true))
  ([f N normalize?]
   (sample* :sample-window f (check-n! :sample-window N)
            (canonical-normalize :sample-window (if (nil? normalize?) true normalize?)))))

(defn idft
  "Calculates the real part of the inverse discrete Fourier transform of window spectrum samples (direct O(N^2) sum).

  Used to build the windows defined by their spectrum (`dolph-chebyshev`, `ultraspherical`, `saramaki`, `legendre`). The sample `n` of the result is `(1/N) sum_k Wn[k] cos(2 pi k (n - (N-1)/2) / N)`, so a symmetric spectrum gives a symmetric window centered at `(N-1)/2`.

  Parameters:

  - `Wn` - sequence of `N` numbers, the spectrum samples `k = 0, ..., N-1`

  Returns a vector of `N` doubles. For `N > 200` the samples are computed in parallel (`pmap`). The cost is quadratic: about 10 ms for `N = 1000`, 0.5 s for `N = 10000`. An empty `Wn` gives an empty vector.

  See also [[dolph-chebyshev]], [[ultraspherical]], [[legendre]]."
  [Wn]
  (let [N (count Wn)
        hN (m/* 0.5 (m/dec N))]
    (-> ((if (m/> N 200) pmap map) (fn [^long n]
                                     (let [f (m/* m/TWO_PI (m// (m/- n hN) N))]
                                       (->> Wn
                                            (map-indexed (fn [^long k ^double wn]
                                                           (m/* wn (m/cos (m/* k f)))))
                                            (v/sum)))) (range N))
        (v/div N)
        (vec))))

;; numeric helpers

(def ^:private quad-options {:rel 1.0e-14 :abs 1.0e-15 :max-iters 400 :initdiv 4})

(defn- integrate-half
  "Integral over [-0.5, 0.5] of the even function `f` (twice the integral over [0, 0.5]; the cusp at 0 is an end point)."
  ^double [f]
  (m/* 2.0 (double (qint/gk-quadrature f 0.0 0.5 quad-options))))

(defn- integrate-phi
  "Integral over [-0.5, 0.5] of `F(s)`, `s = sqrt(1 - 4x^2)`, where `h(s) = s F(s)` is smooth: with `2x = sin(phi)` the integral is `integral from 0 to pi/2 of h(cos(phi)) dphi` (no square root singularity)."
  ^double [h]
  (double (qint/gk-quadrature (fn ^double [^double phi] (double (h (m/cos phi)))) 0.0 m/HALF_PI quad-options)))

(defn- gamma-ratio
  "`Gamma(a) / Gamma(b)` for `a, b > 0`, through log-gamma when the gamma function overflows."
  ^double [^double a ^double b]
  (gamma/gamma-ratio [a] [b]))

(defn- bessel-scaled-asymptotic
  "`exp(-z) I_nu(z)` from the asymptotic series for large `z`, `mu = 4 nu^2`."
  ^double [^double mu ^double z]
  (let [z8 (m/* 8.0 z)]
    (loop [k (long 1)
           term 1.0
           sum 1.0]
      (if (m/> k 14)
        (m// sum (m/sqrt (m/* m/TWO_PI z)))
        (let [kk (m/dec (m/* 2 k))
              term (m// (m/* term (m/- (m/- mu (m/* kk kk)))) (m/* k z8))]
          (recur (m/inc k) term (m/+ sum term)))))))

(defn- bessel-i0e
  "`exp(-z) I_0(z)` for `z >= 0`."
  ^double [^double z]
  (if (m/<= z 500.0)
    (m/* (m/exp (m/- z)) (special/bessel-I0 z))
    (bessel-scaled-asymptotic 0.0 z)))

(defn- bessel-i1e
  "`exp(-z) I_1(z)` for `z >= 0`."
  ^double [^double z]
  (if (m/<= z 500.0)
    (m/* (m/exp (m/- z)) (special/bessel-I1 z))
    (bessel-scaled-asymptotic 4.0 z)))

;;

(defn rectangular05-continuous
  "Rectangular window with the ends set to `0.5`, continuous function: `1` for `abs(x) < 0.5`, `0.5` for `abs(x) = 0.5`, `0` outside (integral `1`)."
  ^double [^double x]
  (let [ax (m/abs x)]
    (cond
      (m/< ax 0.5) 1.0
      (m/== ax 0.5) 0.5
      :else 0.0)))

(defn rectangular05
  "Rectangular window with both end samples set to `0.5`.

  The window is `1` inside and `0.5` at the two ends (the function value at `abs(x) = 0.5`). It is the rectangular window seen as the limit of a symmetric smooth window, useful when the ends should carry half weight. Unlike the other continuous model windows the default `:normalize?` is `false`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `false`: the raw samples (`0.5 1 ... 1 0.5`, the continuous function has the integral `1`); `true` divides by the value at the middle (the maximum becomes `1`); other values: `:LInf` (alias `:Linf`), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[rectangular]], [[triangular]], [[tukey]]."
  ([] (rectangular05 256))
  ([N] (rectangular05 N nil))
  ([N opts]
   (continuous-window :rectangular05 rectangular05-continuous N (window-mode :rectangular05 opts false))))

;;

(defn rectangular-continuous
  "Rectangular window, continuous function: `1` on `[-0.5, 0.5]`, `0` outside (integral `1`)."
  ^double [^double x]
  (clamp x 1.0))

(defn rectangular
  "Rectangular (boxcar) window, all samples are equal to `1`.

  Equal to `scipy.signal.windows.boxcar` and Octave `rectwin`. It has the narrowest main lobe and the highest sidelobes (-13.3 dB) of all windows.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Equal to `tukey` with `:alpha 0`, `kaiser` with `:beta 0` and `power-of-cosine` with `:m 0`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[rectangular05]], [[tukey]], [[triangular]], [[hann]]."
  ([] (rectangular 256))
  ([N] (rectangular N nil))
  ([N opts]
   (continuous-window :rectangular rectangular-continuous N (window-mode :rectangular opts true))))

;;

(defn triangular-continuous
  "Triangular window, continuous function `2 (1 - 2 abs(x))` on `[-0.5, 0.5]` (integral `1`, zero at the ends)."
  ^double [^double x]
  (clamp x (m/* 2.0 (m/- 1.0 (m/* 2.0 (m/abs x))))))

(defn ->triangular-discrete
  "Returns a function of `N` (a long `>= 2`) which returns the discrete triangular window `1 - abs(n - (N-1)/2) / floor((N - 1 + shift) / 2)`, `n = 0, ..., N-1`, for the integer `shift >= 1`; `shift = 2` is equal to `scipy.signal.windows.triang`."
  [^long shift]
  (fn [^long N]
    (let [N- (dec N)
          L (m/+ N- shift)
          hN- (m// N- 2.0)
          hL (m// L 2)]
      (map (fn [^long n] (m/- 1.0 (m/abs (m// (m/- n hN-) hL)))) (range N)))))

(defn triangular
  "Triangular window with a selectable end behaviour.

  The continuous triangle `1 - abs(2x)` (zero ends, Bartlett) is used for `:shift 0`; for a positive `:shift` the window is the discrete formula `1 - abs(n - (N-1)/2) / floor((N - 1 + shift) / 2)`, `n = 0, ..., N-1`, whose ends are non-zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:shift` - integer `>= 0`, default `2`
          - `0` - the Bartlett window, zero at both ends; it is the continuous triangle sampled at `-0.5, ..., 0.5` (equal to `scipy.signal.windows.bartlett` and to `b-spline` of order `2`)
          - `1` - discrete formula with the divisor `floor(N / 2)`
          - `2` (default) - the discrete formula equal to `scipy.signal.windows.triang` and Octave `triang` (non-zero ends, `1/N` for even `N` and `2/(N+1)` for odd `N`)
          - larger values widen the base of the triangle (the ends are larger)
      - `:normalize?` - default `true`. The discrete paths (`:shift > 0`) have the peak `1` already, so `true` and `false` return the formula as it is (the maximum is `1` for odd `N` and below `1` for even `N`); the modes `:LInf`, `:L1`, `:L2`, `:N` and a finite non-zero number are honoured (see [[normalize-coefficients]]). For `:shift 0` `true` divides by the middle value, `false` returns the raw function with the integral `1` (the peak is `2`).

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`. Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[b-spline]], [[parzen]], [[welch]], [[trapezoid]]."
  ([] (triangular 256))
  ([N] (triangular N nil))
  ([N opts]
   (let [shift (int-param :triangular opts :shift 2 non-negative)
         mode (window-mode :triangular opts true)]
     (if (m/pos? (long shift))
       (discrete-window :triangular (->triangular-discrete shift) N mode false)
       (continuous-window :triangular triangular-continuous N mode)))))

;;

(defn parzen-continuous
  "Parzen window, continuous function (the cubic B-spline of order 4 scaled to `[-0.5, 0.5]`, integral `1`, zero ends)."
  ^double [^double x]
  (clamp x
    (let [ax (m/abs x)
          ax2 (m/* ax ax)]
      (if (m/<= ax 0.25)
        (m/* 8.0 (m/+ (m/- m/THIRD (m/* 8.0 ax2))
                      (m/* 16.0 ax2 ax)))
        (m/* 8.0 (m/+ (m/- m/TWO_THIRDS (m/* 4.0 ax) (m/* 5.333333333333333 ax2 ax))
                      (m/* 8.0 ax2)))))))

(defn parzen-discrete
  "Returns the coefficients of the discrete Parzen window of `N` samples (a sequence of doubles, peak `1`, equal to `scipy.signal.windows.parzen`); `N` is a long `>= 2`."
  [^long N]
  (let [N- (m/dec N)
        hN- (m// N- 2.0)
        hL (m// N 2.0)
        qL (m// N 4.0)]
    (map (fn [^long n] (let [an (m/abs (m/- n hN-))
                            nn (m// an hL)]
                        (if (m/<= an qL)
                          (m/- 1.0 (m/* 6.0 (m/sq nn) (m/- 1.0 nn)))
                          (m/* 2.0 (m/cb (m/- 1.0 nn)))))) (range N))))

(defn parzen
  "Parzen window (de la Vallee-Poussin), a piecewise cubic bell.

  The default is the discrete definition, `1 - 6 u^2 (1 - u)` for `abs(n - (N-1)/2) <= N/4` and `2 (1 - u)^3` otherwise, `u = abs(n - (N-1)/2) / (N/2)`. With `:discrete? false` the continuous function (the cubic B-spline of order 4, equal to `b-spline` with `:order 4`) is sampled instead.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:discrete?` - boolean, default `true`: the discrete definition (equal to `scipy.signal.windows.parzen` and Octave `parzenwin`, peak `1`, non-zero small ends); `false`: the continuous function sampled at `-0.5, ..., 0.5` (zero ends)
      - `:normalize?` - default `true`. The discrete path has the peak `1` already, so `true` and `false` return the formula as it is; `:LInf`, `:L1`, `:L2`, `:N` and a finite non-zero number are honoured (see [[normalize-coefficients]]). For `:discrete? false`, `true` divides by the middle value and `false` returns the raw function with the integral `1`.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`. Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[b-spline]], [[triangular]], [[welch]], [[parzen-algebraic]], [[parzen-cosine]]."
  ([] (parzen 256))
  ([N] (parzen N nil))
  ([N opts]
   (let [discrete? (boolean-param :parzen opts :discrete? true)
         mode (window-mode :parzen opts true)]
     (if discrete?
       (discrete-window :parzen parzen-discrete N mode false)
       (continuous-window :parzen parzen-continuous N mode)))))

;; B-spline

(defn ->b-spline
  "Returns the b-spline window function of the integer `order >= 1`: `order * B(order * x)` where `B` is the centered cardinal b-spline of `order` (degree `order - 1`), the integral over `[-0.5, 0.5]` is `1`; computed with the stable Cox-de Boor recurrence in `O(order^2)` per call."
  [^long order]
  (if (m/<= order 1)
    rectangular-continuous
    (let [dn (double order)]
      (fn ^double [^double x]
        (let [ax (m/abs x)]
          (if (m/>= ax 0.5)
            0.0
            (let [u (m/* dn (m/- 0.5 ax))
                  j (long (m/floor u))
                  t (m/- u (double j))
                  ^doubles b (double-array order)]
              (Array/aset b 0 1.0)
              (loop [k (long 2)]
                (when (m/<= k order)
                  (let [km1 (m// 1.0 (double (m/dec k)))
                        dk (double k)]
                    (loop [i (m/dec k)]
                      (when (m/>= i 0)
                        (let [di (double i)
                              a (m/* (m/+ t di) (Array/aget b i))
                              c (if (m/pos? i) (m/* (m/- dk t di) (Array/aget b (m/dec i))) 0.0)]
                          (Array/aset b i (m/* (m/+ a c) km1))
                          (recur (m/dec i))))))
                  (recur (m/inc k))))
              (m/* dn (Array/aget b j)))))))))

(defn b-spline
  "B-spline window: the centered cardinal b-spline of the given order, scaled to the support `[-0.5, 0.5]`.

  The window is `order * B(order * x)` with `B` the cardinal b-spline of `order` (degree `order - 1`), the convolution of `order` rectangles; it is non-negative and has the integral `1`. The orders `1`, `2` and `4` are the rectangular, the triangular (Bartlett) and the continuous Parzen windows, a high order approaches a gaussian. The recurrence is stable, the order is limited only to keep the cost `O(order^2)` per sample reasonable.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:order` - integer in `[1, 1000]`, default `3`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  The cost is `O(order^2)` per sample.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[rectangular]], [[triangular]], [[parzen]], [[gaussian]]."
  ([] (b-spline 256))
  ([N] (b-spline N nil))
  ([N opts]
   (let [order (int-param :b-spline opts :order 3 [1.0 true 1000.0 true])]
     (continuous-window :b-spline (->b-spline order) N (window-mode :b-spline opts true)))))

;;

(defn welch-continuous
  "Welch window, continuous function `1.5 (1 - 4 x^2)` on `[-0.5, 0.5]` (integral `1`, zero ends)."
  ^double [^double x]
  (clamp x (m/* 1.5 (m/- 1.0 (m/* 4.0 x x)))))

(defn welch
  "Welch (parabolic) window `1.5 (1 - 4 x^2)`.

  The parabola with zero ends, the window shape of the Epanechnikov kernel (fastmath kernels are the standard forms with support `[-1, 1]`, the windows have support `[-0.5, 0.5]`: `welch(x)` is the Epanechnikov kernel at `2x`, up to scale). It is equal to `parzen-algebraic` with `:gamma 1.0 :u 2.0`. The window is also called Riesz or Parzen window of the second order.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[parzen-algebraic]], [[power-of-cosine]], [[hann]], [[triangular]]."
  ([] (welch 256))
  ([N] (welch N nil))
  ([N opts]
   (continuous-window :welch welch-continuous N (window-mode :welch opts true))))

;;

(defn ->connes
  "Returns the Connes window function `A (1 - (2x/alpha)^2)^2` (integral `1`) for `alpha > 0`."
  [^double alpha]
  (let [a2 (m/* alpha alpha)
        small? (m/<= alpha 1.0)
        b (if small? 0.0 (m// 1.0 a2))
        den (if small?
              (m/+ 3.0 (m/* -10.0 a2) (m/* 15.0 a2 a2))
              (m/+ (m/* 3.0 b b) (m/* -10.0 b) 15.0))]
    (if small?
      (fn ^double [^double x]
        (clamp x (m// (m/* 15.0 (m/sq (m/- a2 (m/* 4.0 x x)))) den)))
      (fn ^double [^double x]
        (clamp x (m// (m/* 15.0 (m/sq (m/- 1.0 (m/* 4.0 x x b)))) den))))))

(defn connes
  "Connes window `A (1 - (2x/alpha)^2)^2`, the square of the Welch parabola.

  `alpha = 1` is the Connes window with zero ends; for `alpha > 1` the ends are non-zero (a pedestal). For `alpha < 1` the function is `A (alpha^2 - 4 x^2)^2`, which is zero at `x = +-alpha/2` and rises again toward the ends: it is not a taper.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[welch]], [[parzen-algebraic]], [[power-of-cosine]]."
  ([] (connes 256))
  ([N] (connes N nil))
  ([N opts]
   (continuous-window :connes (->connes (param :connes opts :alpha 1.0 open-positive)) N (window-mode :connes opts true))))

;;

(defn ->parzen-algebraic
  "Returns the Parzen algebraic window function `A (1 - gamma abs(2x)^u)` (integral `1`) for `0 < gamma <= 1`, `u > 0`."
  [^double gamma ^double u]
  (let [A (m// (m/inc u) (m/+ u (m/- 1.0 gamma)))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/- 1.0 (m/* gamma (m/pow (m/abs (m/* 2.0 x)) u))))))))

(defn parzen-algebraic
  "Parzen algebraic family window `A (1 - gamma abs(2x)^u)`, `A = (u + 1) / (u + 1 - gamma)`.

  `gamma = 1` makes the ends zero; `u = 1` is the triangular window and `u = 2` the Welch window (for `gamma = 1`); `gamma < 1` gives a pedestal (non-zero ends).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:gamma` - number in `(0, 1]`, default `1.0`
      - `:u` - power, number `> 0`, default `3.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[welch]], [[triangular]], [[power-of-cosine]], [[trapezoid]]."
  ([] (parzen-algebraic 256))
  ([N] (parzen-algebraic N nil))
  ([N opts]
   (let [gamma (param :parzen-algebraic opts :gamma 1.0 [0.0 false 1.0 true])
         u (param :parzen-algebraic opts :u 3.0 open-positive)]
     (continuous-window :parzen-algebraic (->parzen-algebraic gamma u) N (window-mode :parzen-algebraic opts true)))))

;;

;; https://www.eng.buffalo.edu/Research/code/jrnl/Window.pdf

(defn ->singla-singh
  "Returns the Singla-Singh window function of the integer `order >= 0`: `2 (1 - I_t(order + 1, order + 1))` with `t = abs(2x)` and `I` the regularized incomplete beta function (integral `1`)."
  [^long order]
  (let [a (double (m/inc order))]
    (fn ^double [^double x]
      (clamp x
        (let [t (m/abs (m/* 2.0 x))]
          (if (m/<= t 0.5)
            (m/* 2.0 (m/- 1.0 (special/regularized-beta t a a)))
            (m/* 2.0 (special/regularized-beta (m/- 1.0 t) a a))))))))

(defn singla-singh
  "Singla-Singh family of windows, the survival function of the symmetric beta distribution.

  The window is `2 (1 - I_t(m + 1, m + 1))`, `t = abs(2x)`, `m` the order, `I` the regularized incomplete beta function (Singla and Singh, Doerry catalog). Order `0` is the triangular window; the window becomes narrower and smoother for higher orders (the main lobe gets wider, the sidelobes lower).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:order` - integer in `[0, 1000]` (a double with an integral value is accepted), default `1`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[triangular]], [[b-spline]], [[gaussian]]."
  ([] (singla-singh 256))
  ([N] (singla-singh N nil))
  ([N opts]
   (let [order (int-param :singla-singh opts :order 1 [0.0 true 1000.0 true])]
     (continuous-window :singla-singh (->singla-singh order) N (window-mode :singla-singh opts true)))))

;;

(defn- sinc-factor
  ^double [^double L]
  (m// 1.0 (integrate-half (fn [^double x] (m/pow (m/sinc (m/* 2.0 x)) L)))))

(defn sinc-continuous
  "Sinc lobe window, continuous function `A sinc(2x)` on `[-0.5, 0.5]` (integral `1`, zero ends), `sinc(z) = sin(pi z) / (pi z)`."
  ^double [^double x] (clamp x (m/* 1.6963819856764424 (m/sinc (m/* 2.0 x)))))
(defn fejer-continuous
  "Fejer window, continuous function `A sinc(2x)^2` on `[-0.5, 0.5]` (integral `1`, zero ends)."
  ^double [^double x] (clamp x (m/* 2.2152728287035983 (m/sq (m/sinc (m/* 2.0 x))))))
(defn de-la-vallee-poussin-continuous
  "De la Vallee-Poussin window, continuous function `A sinc(2x)^4` on `[-0.5, 0.5]` (integral `1`, zero ends)."
  ^double [^double x] (clamp x (m/* 3.008860051924727 (m/fpow (m/sinc (m/* 2.0 x)) 4))))

(defn sinc
  "Sinc lobe window `A sinc(2x)`, `sinc(z) = sin(pi z) / (pi z)`: the main lobe of the sinc function.

  Equal to `lanczos` with `:L 1.0` and to `scipy.signal.windows.lanczos`. The ends are zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[lanczos]], [[fejer]], [[de-la-vallee-poussin]]."
  ([] (sinc 256))
  ([N] (sinc N nil))
  ([N opts]
   (continuous-window :sinc sinc-continuous N (window-mode :sinc opts true))))

(defn fejer
  "Fejer window `A sinc(2x)^2`, the square of the sinc lobe.

  Equal to `lanczos` with `:L 2.0`. The ends are zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[lanczos]], [[sinc]], [[de-la-vallee-poussin]]."
  ([] (fejer 256))
  ([N] (fejer N nil))
  ([N opts]
   (continuous-window :fejer fejer-continuous N (window-mode :fejer opts true))))

(defn de-la-vallee-poussin
  "De la Vallee-Poussin window `A sinc(2x)^4`, the fourth power of the sinc lobe.

  Equal to `lanczos` with `:L 4.0`. The ends are zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[lanczos]], [[sinc]], [[fejer]]."
  ([] (de-la-vallee-poussin 256))
  ([N] (de-la-vallee-poussin N nil))
  ([N opts]
   (continuous-window :de-la-vallee-poussin de-la-vallee-poussin-continuous N (window-mode :de-la-vallee-poussin opts true))))

(defn ->lanczos
  "Returns the Lanczos window function `A sinc(2x)^L` (integral `1`) for `L >= 0`."
  [^double L]
  (let [A (sinc-factor L)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/pow (m/sinc (m/* 2.0 x)) L))))))

(defn lanczos
  "Lanczos family window `A sinc(2x)^L`, the power `L` of the sinc lobe, `sinc(z) = sin(pi z) / (pi z)`.

  `L = 1` is the usual Lanczos window (equal to `scipy.signal.windows.lanczos` and to `sinc`), `L = 2` is `fejer` and `L = 4` is `de-la-vallee-poussin`; `L = 0` is the rectangular window. The default is `L = 3` (the cube of the sinc lobe), not the `L = 1` of scipy.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:L` - power, number `>= 0`, default `3.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[sinc]], [[fejer]], [[de-la-vallee-poussin]]."
  ([] (lanczos 256))
  ([N] (lanczos N nil))
  ([N opts]
   (continuous-window :lanczos (->lanczos (param :lanczos opts :L 3.0 non-negative)) N (window-mode :lanczos opts true))))

;; 

(defn hamming-continuous
  "Hamming window, continuous function `A (1 + (23/27) cos(2 pi x))`, that is `0.54 + 0.46 cos(2 pi x)` scaled to the integral `1`."
  ^double [^double x]
  (clamp x (m/inc (m/* 0.8518518518518518 (m/cos (m/* m/TWO_PI x))))))

(defn hamming
  "Hamming window `0.54 + 0.46 cos(2 pi x)`, the raised cosine with `alpha = 0.54`.

  Equal to `scipy.signal.windows.hamming` and Octave `hamming`. The ends are non-zero (`0.08`), the peak sidelobe level is about -42.7 dB. The coefficient `0.54` is a rounded value of `25/46`, see `hamming-exact`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hamming-exact]], [[hann]], [[raised-cosine]], [[blackman]]."
  ([] (hamming 256))
  ([N] (hamming N nil))
  ([N opts]
   (continuous-window :hamming hamming-continuous N (window-mode :hamming opts true))))

(defn hamming-exact-continuous
  "Hamming window with the exact coefficients, continuous function `A (1 + 0.84 cos(2 pi x))`, that is `25/46 + 21/46 cos(2 pi x)` scaled to the integral `1`."
  ^double [^double x]
  (clamp x (m/inc (m/* 0.84 (m/cos (m/* m/TWO_PI x))))))

(defn hamming-exact
  "Hamming window with the exact coefficients `25/46 + 21/46 cos(2 pi x)`, the raised cosine with `alpha = 25/46`.

  The coefficients cancel the first sidelobe of the Hann window (a zero of the spectrum at 2.5 bins); the peak sidelobe level is about -41.7 dB. Equal to `scipy.signal.windows.general_hamming` with `alpha = 25/46` and to `webster-hamming` with `:v 0`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hamming]], [[hann]], [[raised-cosine]], [[webster-hamming]]."
  ([] (hamming-exact 256))
  ([N] (hamming-exact N nil))
  ([N opts]
   (continuous-window :hamming-exact hamming-exact-continuous N (window-mode :hamming-exact opts true))))

;;

(defn hann-continuous
  "Hann window, continuous function `1 + cos(2 pi x)` on `[-0.5, 0.5]` (integral `1`, zero ends)."
  ^double [^double x]
  (clamp x (m/inc (m/cos (m/* m/TWO_PI x)))))

(defn hann
  "Hann (Hanning) window `0.5 + 0.5 cos(2 pi x)`, the raised cosine with `alpha = 0.5`.

  Equal to `scipy.signal.windows.hann` and Octave `hann` (symmetric; the Octave function `hanning` is a different window without the zero ends). The ends are zero, the peak sidelobe level is -31.5 dB, the sidelobes drop as `1/f^3`. Equal to `power-of-cosine` with `:m 2`, to `tukey` with `:alpha 1` and to `blackman-harris-family` with `:coeffs [0.5 0.5]`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hamming]], [[raised-cosine]], [[tukey]], [[power-of-cosine]], [[blackman]]."
  ([] (hann 256))
  ([N] (hann N nil))
  ([N opts]
   (continuous-window :hann hann-continuous N (window-mode :hann opts true))))

;;

(defn ->raised-cosine
  "Returns the raised cosine window function `1 + ((1 - alpha) / alpha) cos(2 pi x)` (integral `1`) for `0.5 <= alpha <= 1`."
  [^double alpha]
  (let [f (m// (m/- 1.0 alpha) alpha)]
    (fn ^double [^double x]
      (clamp x (m/inc (m/* f (m/cos (m/* m/TWO_PI x))))))))

(defn raised-cosine
  "Raised cosine family `alpha + (1 - alpha) cos(2 pi x)` (up to scale).

  `alpha = 0.5` is `hann`, `alpha = 0.54` is `hamming` and `alpha = 25/46` is `hamming-exact`. Equal to `scipy.signal.windows.general_hamming` with the same `alpha`. For `alpha < 0.5` the window has negative lobes and for `alpha > 1` the maximum is not at the middle, so the range is limited.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number in `[0.5, 1]`, default `0.5`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hann]], [[hamming]], [[hamming-exact]], [[blackman-harris-family]]."
  ([] (raised-cosine 256))
  ([N] (raised-cosine N nil))
  ([N opts]
   (continuous-window :raised-cosine (->raised-cosine (param :raised-cosine opts :alpha 0.5 [0.5 true 1.0 true])) N (window-mode :raised-cosine opts true))))

;; 

(defn ->webster-hamming
  "Returns the Webster (generalized Hamming) window function (integral `1`) for `v > -0.5`."
  [^double v]
  (let [v+2 (m/+ v 2.0)
        alpha (m/- 1.0 (m// (m/+ 21.0 (m/* 6.0 v))
                            (poly/mevalpoly v 23.0 9.0 1.0)))
        alpha' (m/- 1.0 alpha)
        A (m// (m/* 2.0 m/SQRTPI (gamma-ratio (m/* 0.5 (m/+ 4.0 v)) (m/* 0.5 (m/inc v))))
               (m/inc (m/+ alpha v)))]
    (fn ^double [^double x]
      (let [cpx (m/cos (m/* m/PI x))]
        (clamp x (m/* A (m/+ (m/* alpha (m/pow cpx v))
                             (m/* alpha' (m/pow cpx v+2)))))))))

(defn webster-hamming
  "Webster generalized Hamming window `A (alpha(v) cos(pi x)^v + (1 - alpha(v)) cos(pi x)^(v + 2))`.

  `alpha(v) = (v^2 + 3 v + 2) / (v^2 + 9 v + 23)` (Doerry, eq. 147). `v = 0` is `hamming-exact`; a larger `v` gives lower sidelobes (peak sidelobe level -41.7, -53.9, -71.6 and -123.2 dB for `v` 0, 1, 3 and 10, `N` = 2049) and a wider main lobe. For `v < 0` the window is unbounded at the ends (the end samples are large finite values of the sampled `cos(pi/2)`).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:v` - number in `(-0.5, 1000]`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hamming-exact]], [[power-of-cosine]], [[raised-cosine]]."
  ([] (webster-hamming 256))
  ([N] (webster-hamming N nil))
  ([N opts]
   (continuous-window :webster-hamming (->webster-hamming (param :webster-hamming opts :v 1.0 [-0.5 false 1000.0 true])) N (window-mode :webster-hamming opts true))))

;;

(defn ->power-of-cosine
  "Returns the power of cosine window function `A cos(pi x)^m` (integral `1`) for `m >= 0`."
  [^double m]
  (let [A (m/* m/SQRTPI (gamma-ratio (m/* 0.5 (m/+ m 2.0)) (m/* 0.5 (m/inc m))))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/pow (m/cos (m/* m/PI x)) m))))))

(defn power-of-cosine
  "Power of cosine window `A cos(pi x)^m`, `A = sqrt(pi) Gamma((m + 2) / 2) / Gamma((m + 1) / 2)`.

  `m = 0` is rectangular, `m = 1` is the cosine lobe (sampled at `-0.5, ..., 0.5` with zero ends; this is not `scipy.signal.windows.cosine`, which samples the half-sample grid, see [[cosine]]), `m = 2` is `hann`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:m` - power, number in `[0, 1000]`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[cosine]], [[hann]], [[raised-power-of-cosine]], [[welch]]."
  ([] (power-of-cosine 256))
  ([N] (power-of-cosine N nil))
  ([N opts]
   (continuous-window :power-of-cosine (->power-of-cosine (param :power-of-cosine opts :m 1.0 [0.0 true 1000.0 true])) N (window-mode :power-of-cosine opts true))))

;;

(defn ->raised-power-of-cosine
  "Returns the raised power of cosine window function `A (alpha + (1 - alpha) cos(pi x)^m)` (integral `1`) for `0 <= alpha <= 1`, `m >= 0`."
  [^double m ^double alpha]
  (let [alpha' (m/- 1.0 alpha)
        A (m// (m/+ alpha (m// (m/* alpha' (gamma-ratio (m/* 0.5 (m/inc m)) (m/* 0.5 (m/+ m 2.0))))
                               m/SQRTPI)))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/+ alpha (m/* alpha' (m/pow (m/cos (m/* m/PI x)) m))))))))

(defn raised-power-of-cosine
  "Raised power of cosine window `A (alpha + (1 - alpha) cos(pi x)^m)`, a pedestal `alpha` plus a power of cosine.

  `alpha = 0` is `power-of-cosine`, `alpha = 1` is rectangular; `m = 2` gives the raised cosine (`hann` for `alpha = 0.5`).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - pedestal, number in `[0, 1]`, default `0.05`
      - `:m` - power, number in `[0, 1000]`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[power-of-cosine]], [[raised-cosine]], [[hamming]]."
  ([] (raised-power-of-cosine 256))
  ([N] (raised-power-of-cosine N nil))
  ([N opts]
   (let [alpha (param :raised-power-of-cosine opts :alpha 0.05 [0.0 true 1.0 true])
         m (param :raised-power-of-cosine opts :m 1.0 [0.0 true 1000.0 true])]
     (continuous-window :raised-power-of-cosine (->raised-power-of-cosine m alpha) N (window-mode :raised-power-of-cosine opts true)))))

;;

(defn- parzen-cosine-helper
  ^double [^double m ^double gamma ^double x]
  (m/inc (m/cos (m/* m/PI gamma (m/pow (m/abs (m/* 2.0 x)) m)))))

(defn ->parzen-cosine
  "Returns the Parzen cosine window function `A (1 + cos(pi gamma abs(2x)^m))` (integral `1`, the constant by quadrature) for `0 <= gamma <= 1`, `m > 0`."
  [^double m ^double gamma]
  (let [A (m// (integrate-half (partial parzen-cosine-helper m gamma)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-cosine-helper m gamma x))))))

(defn parzen-cosine
  "Parzen cosine family window `A (1 + cos(pi gamma abs(2x)^m))`.

  `gamma = 1` makes the ends zero, `m = 2` gives a flatter top than the Hann window. `gamma < 1` leaves a pedestal at the ends.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:gamma` - number in `[0, 1]`, default `1.0`
      - `:m` - power, number `> 0`, default `2.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hann]], [[raised-cosine]], [[parzen-algebraic]]."
  ([] (parzen-cosine 256))
  ([N] (parzen-cosine N nil))
  ([N opts]
   (let [gamma (param :parzen-cosine opts :gamma 1.0 [0.0 true 1.0 true])
         m (param :parzen-cosine opts :m 2.0 open-positive)]
     (continuous-window :parzen-cosine (->parzen-cosine m gamma) N (window-mode :parzen-cosine opts true)))))

;;

(defn bohman-continuous
  "Bohman window, continuous function `A ((1 - u) cos(pi u) + sin(pi u) / pi)`, `u = abs(2x)` (integral `1`, zero ends)."
  ^double [^double x]
  (let [x2 (m/* 2.0 (m/abs x))
        px2 (m/* m/PI x2)]
    (clamp x (m/+ (m/* 2.4674011002723395 (m/- 1.0 x2) (m/cos px2))
                  (m/* m/QUARTER_PI (m/sin px2))))))

(defn bohman
  "Bohman window `(1 - u) cos(pi u) + sin(pi u) / pi`, `u = abs(2x)`: the convolution of two half cosine lobes.

  Equal to `scipy.signal.windows.bohman` and Octave `bohmanwin`. The ends and the first two derivatives there are zero, the sidelobes drop as `1/f^4`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hann]], [[parzen]], [[cosine]]."
  ([] (bohman 256))
  ([N] (bohman N nil))
  ([N opts]
   (continuous-window :bohman bohman-continuous N (window-mode :bohman opts true))))

;;

(defn ->trapezoid
  "Returns the trapezoid window function (flat for `abs(x) <= alpha`, linear to zero at the ends, integral `1`) for `0 <= alpha <= 0.5`."
  [^double alpha]
  (let [alpha2 (m/* 2.0 alpha)
        salpha2 (m/- 1.0 (m/sq alpha2))]
    (fn ^double [^double x]
      (let [ax (m/abs x)]
        (cond
          (m/<= ax alpha) (m// 2.0 (m/inc alpha2))
          (m/<= ax 0.5) (m// (m/* 2.0 (m/- 1.0 (m/* 2.0 ax)))
                             salpha2)
          :else 0.0)))))

(defn trapezoid
  "Trapezoid window, flat for `abs(x) <= alpha` and linear to zero at `abs(x) = 0.5`.

  `alpha = 0` is the triangular (Bartlett) window and `alpha = 0.5` the rectangular window. It is the convolution of a rectangle and a triangle.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - half width of the flat part, number in `[0, 0.5]`, default `0.25`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[triangular]], [[rectangular]], [[tukey]]."
  ([] (trapezoid 256))
  ([N] (trapezoid N nil))
  ([N opts]
   (continuous-window :trapezoid (->trapezoid (param :trapezoid opts :alpha 0.25 [0.0 true 0.5 true])) N (window-mode :trapezoid opts true))))

;;

(defn ->tukey
  "Returns the Tukey (tapered cosine) window function (integral `1`); `alpha >= 0` is the tapered fraction of the window, `alpha >= 1` is the Hann window."
  [^double alpha]
  (if (m/zero? alpha)
    rectangular-continuous
    (let [a (m/min alpha 1.0)
          w (m/* 0.5 a)
          a' (m/- 0.5 w)
          h (m// 2.0 (m/- 2.0 a))
          hh (m/* 0.5 h)
          cw (m// m/PI w)]
      (fn ^double [^double x]
        (let [ax (m/abs x)]
          (cond
            (m/<= ax a') h
            (m/<= ax 0.5) (m/* hh (m/inc (m/cos (m/* cw (m/- ax a')))))
            :else 0.0))))))

(defn tukey
  "Tukey (tapered cosine) window: a flat middle part with raised cosine tapers at both ends.

  `alpha` is the fraction of the window length which is tapered (`alpha / 2` at each end). `alpha = 0` is rectangular and `alpha = 1` is `hann`; like `scipy.signal.windows.tukey` every `alpha > 1` gives the Hann window. Equal to scipy `tukey` and Octave `tukeywin` for `0 <= alpha <= 1`. The raw window has a plateau `2 / (2 - alpha)` and the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `>= 0`, default `0.5`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hann]], [[rectangular]], [[trapezoid]], [[raised-cosine]]."
  ([] (tukey 256))
  ([N] (tukey N nil))
  ([N opts]
   (continuous-window :tukey (->tukey (param :tukey opts :alpha 0.5 non-negative)) N (window-mode :tukey opts true))))

;;

(defn bartlett-hann-continuous
  "Bartlett-Hann window, continuous function `2 (0.62 - 0.48 abs(x) + 0.38 cos(2 pi x))` on `[-0.5, 0.5]` (integral `1`)."
  ^double [^double x]
  (clamp x (m/* 2.0 (m/+ 0.62 (m/* -0.48 (m/abs x))
                         (m/* 0.38 (m/cos (m/* m/TWO_PI x)))))))

(defn bartlett-hann
  "Bartlett-Hann window `0.62 - 0.48 abs(x) + 0.38 cos(2 pi x)`, a linear combination of the triangular and the Hann window.

  Equal to `scipy.signal.windows.barthann` and Octave `barthannwin`. The ends are almost zero (`0.0`, up to rounding).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[triangular]], [[hann]], [[blackman]]."
  ([] (bartlett-hann 256))
  ([N] (bartlett-hann N nil))
  ([N opts]
   (continuous-window :bartlett-hann bartlett-hann-continuous N (window-mode :bartlett-hann opts true))))

;;

(defn blackman-continuous
  "Blackman window, continuous function `A (0.42 + 0.5 cos(2 pi x) + 0.08 cos(4 pi x))` (integral `1`, zero ends)."
  ^double [^double x]
  (let [px (m/* m/TWO_PI x)]
    (clamp x (m/inc (m/+ (m/* 1.1904761904761905 (m/cos px))
                         (m/* 0.1904761904761905 (m/cos (m/* 2.0 px))))))))


(defn blackman
  "Blackman window `0.42 + 0.5 cos(2 pi x) + 0.08 cos(4 pi x)`, three cosine terms (centered form).

  Equal to `scipy.signal.windows.blackman` and Octave `blackman`. The ends are zero, the peak sidelobe level is about -58.1 dB, the sidelobes drop as `1/f^3`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-exact]], [[hann]], [[hamming]], [[blackman-harris-family]]."
  ([] (blackman 256))
  ([N] (blackman N nil))
  ([N opts]
   (continuous-window :blackman blackman-continuous N (window-mode :blackman opts true))))


(defn blackman-exact-continuous
  "Exact Blackman window, continuous function `A (7938 + 9240 cos(2 pi x) + 1430 cos(4 pi x)) / 18608` (integral `1`)."
  ^double [^double x]
  (let [px (m/* m/TWO_PI x)]
    (clamp x (m/inc (m/+ (m/* 1.164021164021164 (m/cos px))
                         (m/* 0.18014613252708492 (m/cos (m/* 2.0 px))))))))

(defn blackman-exact
  "Exact Blackman window `(7938 + 9240 cos(2 pi x) + 1430 cos(4 pi x)) / 18608`, the three-term cosine sum which cancels the first two sidelobes.

  Equal to `scipy.signal.windows.general_cosine` with the coefficients `[7938 9240 1430] / 18608`. The ends are small (`0.0069`), the peak sidelobe level is about -68.2 dB.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman]], [[hamming-exact]], [[blackman-harris-family]]."
  ([] (blackman-exact 256))
  ([N] (blackman-exact N nil))
  ([N opts]
   (continuous-window :blackman-exact blackman-exact-continuous N (window-mode :blackman-exact opts true))))


(defn ->blackman-harris
  "Returns the sum of cosines window function `A (a0 + a1 cos(2 pi x) + a2 cos(4 pi x) + ...)` with `A = 1 / a0` (integral `1`); `coeffs` are the positive magnitudes `[a0 a1 ...]` (the centered form of the classical alternating signs formula, as `scipy.signal.windows.general_cosine`). With `clip?` `true` (default) the negative values of the sum are replaced by `0.0`, with `false` the sum is kept as it is."
  ([coeffs] (->blackman-harris coeffs true))
  ([coeffs clip?]
   (let [^doubles cs (m/seq->double-array coeffs)
         n (alength cs)
         a0 (Array/aget cs 0)
         A (m// a0)]
     (fn ^double [^double x]
       (if (m/<= -0.5 x 0.5)
         (let [v (m/* m/TWO_PI x)
               r (m/* A (double (loop [id (long 1)
                                       sum a0]
                                  (if (m/< id n)
                                    (recur (m/inc id)
                                           (m/+ sum (m/* (Array/aget cs id) (m/cos (m/* id v)))))
                                    sum))))]
           (if clip? (m/max 0.0 r) r))
         0.0)))))

(defn blackman-harris-family
  "Sum of cosines window `A (a0 + a1 cos(2 pi x) + a2 cos(4 pi x) + ...)` for `x` in `[-0.5, 0.5]` with user coefficients; `A = 1 / a0`, so the integral is `1`.

  The terms are added with the positive magnitudes `a0 a1 ...` (the centered form). The classical formula on the sample index, `a0 - a1 cos(2 pi n / (N-1)) + a2 cos(4 pi n / (N-1)) - ...`, has alternating signs; both are the same window, and this one is `scipy.signal.windows.general_cosine`. Hann is `[0.5 0.5]`, Hamming `[0.54 0.46]`, Blackman `[0.42 0.5 0.08]`, the 4-term Blackman-Harris `[0.35875 0.48829 0.14128 0.01168]`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:coeffs` - non-empty sequence of finite numbers `[a0 a1 ...]` with `a0 > 0`, default `[0.5 0.5]` (Hann)
      - `:normalize?` - default `true`: the samples are divided by the value at the middle, `sum(a)` (the maximum is `1` for odd `N` when the window is positive in the middle); `false` returns the raw samples of the function with the integral `1`. Other values: `:LInf` (alias `:Linf`), `:L1` (divides by the sum of the absolute values), `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  The negative values of the sum are kept: the family is not clipped (the flat-top coefficient sets have negative lobes), unlike the named positive windows such as `blackman-harris`. Throws `ex-info` (`:type :invalid-input`) for empty or non-numeric `:coeffs`, `a0 <= 0` and an invalid `N` or `:normalize?`.

  See also [[->blackman-harris]], [[hann]], [[blackman]], [[flat-top]], [[flat-top-bk]], [[hft248d]]."
  ([] (blackman-harris-family 256))
  ([N] (blackman-harris-family N nil))
  ([N opts]
   (let [coeffs (let [c (get opts :coeffs)]
                  (cond
                    (nil? c) [0.5 0.5]
                    (and (or (sequential? c) (and (some? c) (.isArray (class c))))
                         (seq c)
                         (every? (fn [a] (and (number? a) (m/valid-double? (double a)))) c)
                         (m/pos? (double (first c)))) c
                    :else (invalid-input! :blackman-harris-family :coeffs
                                          "must be a non-empty sequence of finite numbers with a positive first element"
                                          c "non-empty sequence of finite numbers, first element > 0")))]
     (continuous-window :blackman-harris-family (->blackman-harris coeffs false) N (window-mode :blackman-harris-family opts true)))))

(defmacro make-blackman-harris
  "Macro which expands to the clipped or unclipped sum of cosines `A (a0 + a1 cos(2 pi x) + a2 cos(4 pi x) + ...)` with literal coefficients, `A = 1 / a0`.

  Parameters:

  - `clamp-fn` - one of the private clamping macros of this namespace: `clamp` replaces negative values by `0.0`, `clamp-` keeps them (both give `0.0` outside of `[-0.5, 0.5]`)
  - `v` - the argument symbol
  - `coeffs` - the positive magnitudes `a0 a1 ...` as number literals

  Returns the code of the continuous window function body, used by the `-continuous` functions of the cosine sum windows. See also [[->blackman-harris]]."
  [clamp-fn v & coeffs]
  (let [a0 (double (first coeffs))
        A (m// a0)]
    `(~clamp-fn ~v
      (let [~'xx (m/* m/TWO_PI ~v)]
        (m/* ~A (m/+ ~a0 ~@(map-indexed (fn [^long id ^double c]
                                          (if (m/zero? id)
                                            `(m/* ~c (m/cos ~'xx))
                                            `(m/* ~c (m/cos (m/* ~(m/inc id) ~'xx))))) (rest coeffs))))))))

(defn blackman-harris-continuous
  "Blackman-Harris (Nuttall three-term) window, continuous function `A (0.4243801 + 0.4973406 cos(2 pi x) + 0.0782793 cos(4 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.4243801 0.4973406 0.0782793))
(defn blackman-harris-61db-continuous
  "Blackman-Harris -61 dB window, continuous function `A (0.44959 + 0.49364 cos(2 pi x) + 0.05677 cos(4 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.44959 0.49364 0.05677))
(defn blackman-harris-67db-continuous
  "Blackman-Harris -67 dB window, continuous function `A (0.42323 + 0.49755 cos(2 pi x) + 0.07922 cos(4 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.42323 0.49755 0.07922))
(defn blackman-harris-74db-continuous
  "Blackman-Harris -74 dB window, continuous function `A (0.40217 + 0.49703 cos(2 pi x) + 0.09892 cos(4 pi x) + 0.00188 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.40217 0.49703 0.09892 0.00188))
(defn blackman-harris-92db-continuous
  "Blackman-Harris -92 dB window, continuous function `A (0.35875 + 0.48829 cos(2 pi x) + 0.14128 cos(4 pi x) + 0.01168 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.35875 0.48829 0.14128 0.01168))
(defn nutall-3-1st-continuous
  "Nuttall three-term (first derivative continuous) window, continuous function `A (0.40897 + 0.5 cos(2 pi x) + 0.09103 cos(4 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.40897 0.5 0.09103))
(defn nutall-3-3rd-continuous
  "Nuttall three-term (third derivative continuous) window, continuous function `A (0.375 + 0.5 cos(2 pi x) + 0.125 cos(4 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.375 0.5 0.125))
(defn blackman-nutall-continuous
  "Blackman-Nuttall window, continuous function `A (0.3635819 + 0.4891775 cos(2 pi x) + 0.1365995 cos(4 pi x) + 0.0106411 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.3635819 0.4891775 0.1365995 0.0106411))
(defn nutall-1st-continuous
  "Nuttall four-term (first derivative continuous) window, continuous function `A (0.355768 + 0.487396 cos(2 pi x) + 0.144232 cos(4 pi x) + 0.012604 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.355768 0.487396 0.144232 0.012604))
(defn nutall-3rd-continuous
  "Nuttall four-term (third derivative continuous) window, continuous function `A (0.338946 + 0.481973 cos(2 pi x) + 0.161054 cos(4 pi x) + 0.018027 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.338946 0.481973 0.161054 0.018027))
(defn nutall-5th-continuous
  "Nuttall four-term (fifth derivative continuous) window, continuous function `A (0.3125 + 0.46875 cos(2 pi x) + 0.1875 cos(4 pi x) + 0.03125 cos(6 pi x))` on `[-0.5, 0.5]` (integral `1`), `0` outside; negative values are clipped to `0`."
  ^double [^double x] (make-blackman-harris clamp x 0.3125 0.46875 0.1875 0.03125))

(defn blackman-harris
  "Blackman-Harris window in the Nuttall three-term minimum sidelobe form (Heinzel et al. 2002: Nuttall3b).

  The window is `A (0.4243801 + 0.4973406 cos(2 pi x) + 0.0782793 cos(4 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -71.5 dB. This is not `scipy.signal.windows.blackmanharris` or Octave `blackmanharris` (those are the four-term -92 dB window, [[blackman-harris-92db]]), and not scipy `nuttall` ([[blackman-nutall]]). Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-harris-92db]], [[blackman-nutall]], [[blackman-harris-67db]], [[blackman-harris-family]]."
  ([] (blackman-harris 256))
  ([N] (blackman-harris N nil))
  ([N opts]
   (continuous-window :blackman-harris blackman-harris-continuous N (window-mode :blackman-harris opts true))))

(defn blackman-harris-61db
  "Blackman-Harris three-term window (Harris 1978), named -61 dB after the rounded sidelobe level of the paper.

  The window is `A (0.44959 + 0.49364 cos(2 pi x) + 0.05677 cos(4 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -62.0 dB. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-harris]], [[blackman-harris-67db]], [[blackman-harris-74db]]."
  ([] (blackman-harris-61db 256))
  ([N] (blackman-harris-61db N nil))
  ([N opts]
   (continuous-window :blackman-harris-61db blackman-harris-61db-continuous N (window-mode :blackman-harris-61db opts true))))

(defn blackman-harris-67db
  "Blackman-Harris three-term window (Harris 1978), named -67 dB after the rounded sidelobe level of the paper.

  The window is `A (0.42323 + 0.49755 cos(2 pi x) + 0.07922 cos(4 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -70.8 dB. The name is Harris' label; the sidelobe level of the coefficients is lower than the label (the measured level is the value above). Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-harris]], [[blackman-harris-61db]], [[blackman-harris-74db]]."
  ([] (blackman-harris-67db 256))
  ([N] (blackman-harris-67db N nil))
  ([N opts]
   (continuous-window :blackman-harris-67db blackman-harris-67db-continuous N (window-mode :blackman-harris-67db opts true))))

(defn blackman-harris-74db
  "Blackman-Harris four-term window (Harris 1978), named -74 dB after the rounded sidelobe level of the paper.

  The window is `A (0.40217 + 0.49703 cos(2 pi x) + 0.09892 cos(4 pi x) + 0.00188 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -74.4 dB. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-harris]], [[blackman-harris-92db]], [[blackman-harris-67db]]."
  ([] (blackman-harris-74db 256))
  ([N] (blackman-harris-74db N nil))
  ([N opts]
   (continuous-window :blackman-harris-74db blackman-harris-74db-continuous N (window-mode :blackman-harris-74db opts true))))

(defn blackman-harris-92db
  "Blackman-Harris four-term window (Harris 1978) with the sidelobe level of -92 dB.

  The window is `A (0.35875 + 0.48829 cos(2 pi x) + 0.14128 cos(4 pi x) + 0.01168 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -92.0 dB. Equal to `scipy.signal.windows.blackmanharris`, Octave `blackmanharris` and Matlab `blackmanharris`. The ends are `6e-5` of the middle value (almost zero), the sidelobes drop as `1/f`. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-harris]], [[blackman-nutall]], [[nutall-1st]], [[blackman-harris-family]]."
  ([] (blackman-harris-92db 256))
  ([N] (blackman-harris-92db N nil))
  ([N opts]
   (continuous-window :blackman-harris-92db blackman-harris-92db-continuous N (window-mode :blackman-harris-92db opts true))))

(defn nutall-3-1st
  "Nuttall three-term window with continuous first derivative (Nuttall 1981; Heinzel et al. 2002: Nuttall3a).

  The window is `A (0.40897 + 0.5 cos(2 pi x) + 0.09103 cos(4 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -64.2 dB. The name keeps the spelling `nutall` of the library. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[nutall-3-3rd]], [[blackman-harris]], [[nutall-1st]]."
  ([] (nutall-3-1st 256))
  ([N] (nutall-3-1st N nil))
  ([N opts]
   (continuous-window :nutall-3-1st nutall-3-1st-continuous N (window-mode :nutall-3-1st opts true))))

(defn nutall-3-3rd
  "Nuttall three-term window with continuous third derivative (Nuttall 1981; Heinzel et al. 2002: Nuttall3).

  The window is `A (0.375 + 0.5 cos(2 pi x) + 0.125 cos(4 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -46.7 dB. It is `cos(pi x)^4` up to scale, equal to `power-of-cosine` with `:m 4`; the sidelobes drop as `1/f^5`. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[nutall-3-1st]], [[power-of-cosine]], [[hann]]."
  ([] (nutall-3-3rd 256))
  ([N] (nutall-3-3rd N nil))
  ([N opts]
   (continuous-window :nutall-3-3rd nutall-3-3rd-continuous N (window-mode :nutall-3-3rd opts true))))

(defn blackman-nutall
  "Blackman-Nuttall four-term window (Nuttall 1981; Heinzel et al. 2002: Nuttall4c).

  The window is `A (0.3635819 + 0.4891775 cos(2 pi x) + 0.1365995 cos(4 pi x) + 0.0106411 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -98.2 dB. Equal to `scipy.signal.windows.nuttall`. Octave `nuttallwin` is a different window, [[nutall-1st]]. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[nutall-1st]], [[blackman-harris-92db]], [[nutall-3rd]], [[blackman-harris-family]]."
  ([] (blackman-nutall 256))
  ([N] (blackman-nutall N nil))
  ([N opts]
   (continuous-window :blackman-nutall blackman-nutall-continuous N (window-mode :blackman-nutall opts true))))

(defn nutall-1st
  "Nuttall four-term window with continuous first derivative (Nuttall 1981; Heinzel et al. 2002: Nuttall4b).

  The window is `A (0.355768 + 0.487396 cos(2 pi x) + 0.144232 cos(4 pi x) + 0.012604 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -93.3 dB. Equal to Octave `nuttallwin`. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[blackman-nutall]], [[nutall-3rd]], [[nutall-5th]], [[blackman-harris-92db]]."
  ([] (nutall-1st 256))
  ([N] (nutall-1st N nil))
  ([N opts]
   (continuous-window :nutall-1st nutall-1st-continuous N (window-mode :nutall-1st opts true))))

(defn nutall-3rd
  "Nuttall four-term window with continuous third derivative (Nuttall 1981; Heinzel et al. 2002: Nuttall4a).

  The window is `A (0.338946 + 0.481973 cos(2 pi x) + 0.161054 cos(4 pi x) + 0.018027 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -82.6 dB. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[nutall-1st]], [[nutall-5th]], [[blackman-nutall]]."
  ([] (nutall-3rd 256))
  ([N] (nutall-3rd N nil))
  ([N opts]
   (continuous-window :nutall-3rd nutall-3rd-continuous N (window-mode :nutall-3rd opts true))))

(defn nutall-5th
  "Nuttall four-term window with continuous fifth derivative (Nuttall 1981; Heinzel et al. 2002: Nuttall4).

  The window is `A (0.3125 + 0.46875 cos(2 pi x) + 0.1875 cos(4 pi x) + 0.03125 cos(6 pi x))` (the centered form with the positive magnitudes, `A = 1 / a0`). The peak sidelobe level is about -61.0 dB. It is `cos(pi x)^6` up to scale, equal to `power-of-cosine` with `:m 6`; the sidelobes drop as `1/f^7`. Negative values would be replaced by zero (there are none for these coefficients). The raw window has the integral `1`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[nutall-3rd]], [[power-of-cosine]], [[nutall-3-3rd]]."
  ([] (nutall-5th 256))
  ([N] (nutall-5th N nil))
  ([N opts]
   (continuous-window :nutall-5th nutall-5th-continuous N (window-mode :nutall-5th opts true))))

(defn mottaghi-kashtiban-shayesteh
  "Mottaghi-Kashtiban-Shayesteh window: a four-term cosine sum whose coefficients depend on the window length.

  The window is `a0 + a1 cos(2 pi x) + a2 cos(4 pi x) + a3 cos(6 pi x)` (centered form, scaled by `1 / a0`) with `a0 = 0.5363 - 0.14 / (N - 1)`, `a1 = 0.996 - a0`, `a2 = 0`, `a3 = 0.004` (Doerry, catalog, table 1). It improves on the Hamming window (peak sidelobe level about -44.6 dB for `N` = 256 against -42.7 dB).

  Parameters:

  - `N` - number of samples, an integer `>= 0`; the coefficients are computed for the given `N`
  - `options` (optional map):
      - `:normalize?` - default `true`: division by the value at the middle; `false` returns the raw samples, which have the integral `1` of the function; other values: `:LInf`, `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`. The formula is not clipped: the window is non-negative for `N >= 5`, for `N` from `2` to `4` the end values are negative. Because the coefficients depend on `N`, a periodic window of `fastmath.kernel/window` (`:symmetric? false`) is the first `N` samples of the window designed for `N + 1`, not a periodic design for `N`. Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hamming]], [[blackman-harris-family]], [[shayesteh-kashtiban]]."
  ([] (mottaghi-kashtiban-shayesteh 256))
  ([N] (mottaghi-kashtiban-shayesteh N nil))
  ([N opts]
   (let [mode (window-mode :mottaghi-kashtiban-shayesteh opts true)
         n (check-n! :mottaghi-kashtiban-shayesteh N)]
     (if (m/< n 2)
       (continuous-window :mottaghi-kashtiban-shayesteh rectangular-continuous n mode)
       (let [a0 (m/- 0.5363 (m// 0.14 (m/dec n)))]
         (continuous-window :mottaghi-kashtiban-shayesteh (->blackman-harris [a0 (m/- 0.996 a0) 0.0 0.004] false) n mode))))))

;; 

(defn ->exponential
  "Returns the exponential window function `A exp(-2 alpha abs(x))` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [A (m// alpha (m/- (m/expm1 (m/- alpha))))
        f (m/* -2.0 alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f (m/abs x)))))))) 

(defn exponential
  "Exponential (Poisson) window `A exp(-2 alpha abs(x))`, `A = alpha / (1 - exp(-alpha))`.

  The value at the ends is `exp(-alpha)` of the middle value. `alpha` is the decay over the half length. Equal to `scipy.signal.windows.exponential` with the center at the middle and `tau = (N - 1) / (2 alpha)`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - decay, number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hanning-poisson]], [[gaussian]], [[cauchy]], [[parzen-exponential]]."
  ([] (exponential 256))
  ([N] (exponential N nil))
  ([N opts]
   (continuous-window :exponential (->exponential (param :exponential opts :alpha 1.0 open-positive)) N (window-mode :exponential opts true))))

;;

(defn ->hanning-poisson
  "Returns the Hanning-Poisson window function (the Hann window times `exp(-2 alpha abs(x))`, integral `1`), `alpha > 0`."
  [^double alpha]
  (let [a2 (m/* alpha alpha)
        A (m// (m/* alpha (m/+ a2 m/PI2))
               (m/- (m/* 2.0 a2) (m/* m/PI2 (m/expm1 (m/- alpha)))))
        f (m/* -2.0 alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f (m/abs x))) (m/inc (m/cos (m/* m/TWO_PI x))))))))

(defn hanning-poisson
  "Hanning-Poisson window, the product of the Hann and the exponential window `A (1 + cos(2 pi x)) exp(-2 alpha abs(x))`.

  The ends are zero; `alpha` is the decay of the exponential factor.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - decay, number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[hann]], [[exponential]], [[gaussian]]."
  ([] (hanning-poisson 256))
  ([N] (hanning-poisson N nil))
  ([N opts]
   (continuous-window :hanning-poisson (->hanning-poisson (param :hanning-poisson opts :alpha 1.0 open-positive)) N (window-mode :hanning-poisson opts true))))

;;

(defn ->gaussian
  "Returns the truncated gaussian window function `A exp(-2 alpha^2 x^2)` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [A (m// (m/* alpha m/SQRT_2_PI)
               (special/erf (m// alpha m/SQRT2)))
        f (m/* -2.0 alpha alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f x x)))))))

(defn gaussian
  "Truncated gaussian window `A exp(-2 alpha^2 x^2)`, `A = alpha sqrt(2/pi) / erf(alpha / sqrt(2))`.

  The standard deviation is `sigma = 1 / (2 alpha)` of the window length, so `alpha` is the number of standard deviations from the middle to each end. Equal to `scipy.signal.windows.gaussian` with `std = (N - 1) / (2 alpha)` and to Octave `gausswin(N, alpha)`. The ends are non-zero (`exp(-alpha^2 / 2)` of the middle).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - truncation in standard deviations, number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[general-gaussian]], [[exponential]], [[kaiser]], [[parzen-exponential]]."
  ([] (gaussian 256))
  ([N] (gaussian N nil))
  ([N opts]
   (continuous-window :gaussian (->gaussian (param :gaussian opts :alpha 1.0 open-positive)) N (window-mode :gaussian opts true))))

;; generalized gaussian

(defn- general-gaussian-integral
  "Integral over [-0.5, 0.5] of `exp(-0.5 abs(x / sig)^(2 p))`: with `a = 1 / (2 p)` and `z = 0.5 (0.5 / sig)^(2 p)` it is `exp(-z) sum_k z^k / ((a + 1) ... (a + k))` for moderate `z` (positive terms only) and `2 sig 2^a Gamma(1 + a) P(a, z)` otherwise (`P` the regularized lower incomplete gamma function)."
  ^double [^double p ^double sig]
  (let [a (m// 0.5 p)
        z (m/exp (m/- (m/* 2.0 p (m/log (m// 0.5 sig))) m/LN2))]
    (if (m/<= z 50.0)
      (m/* (m/exp (m/- z))
           (double (loop [k (long 1)
                          term 1.0
                          sum 1.0]
                     (let [term (m// (m/* term z) (m/+ a k))
                           sum (m/+ sum term)]
                       (if (or (m/<= term (m/* 1.0e-17 sum)) (m/> k 100000))
                         sum
                         (recur (m/inc k) term sum))))))
      (m/* 2.0 sig
           (m/exp (m/+ (m/* a m/LN2) (special/log-gamma (m/inc a))))
           (if (m/valid-double? z) (special/regularized-gamma-p a z) 1.0)))))

(defn ->general-gaussian
  "Returns the generalized gaussian window function `A exp(-0.5 abs(x / sig)^(2 p))` (integral `1`), `p > 0`, `sig > 0` (a fraction of the window length). Throws `ex-info` for a `sig` so small that the window has no non-zero sample."
  [^double p ^double sig]
  (let [A (m// (general-gaussian-integral p sig))
        e (m/* 2.0 p)]
    (when-not (and (m/valid-double? A) (m/pos? A))
      (invalid-input! :general-gaussian :sig "too small, the window has no non-zero sample" sig "sig for which the window has a non-zero sample"))
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* -0.5 (m/pow (m/abs (m// x sig)) e))))))))

(defn general-gaussian
  "Generalized gaussian window `A exp(-0.5 abs(x / sig)^(2 p))` for `x` in `[-0.5, 0.5]`; `p = 1` is the gaussian window, a larger `p` gives a flatter top and steeper sides, a smaller `p` a more peaked shape.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:p` - shape, number in `(0, 1000]`, default `1.5`
      - `:sig` - width, number `> 0`, default `0.25`. It is the width in units of the window length, so the window is `N` independent: `scipy.signal.windows.general_gaussian(N, p, sig * (N - 1))` is the same window (the `sig` of scipy is in samples). `p = 1` and `sig = 1 / (2 alpha)` is `gaussian` with `alpha`.
      - `:normalize?` - default `true` (division by the value at the middle); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` for a `p` or `sig` outside of its domain and for a `sig` so small that the window has no non-zero sample.

  See also [[gaussian]], [[parzen-exponential]] (`r = 2 p`, the same family with a different parametrization)."
  ([] (general-gaussian 256))
  ([N] (general-gaussian N nil))
  ([N opts]
   (let [p (param :general-gaussian opts :p 1.5 [0.0 false 1000.0 true])
         sig (param :general-gaussian opts :sig 0.25 open-positive)]
     (continuous-window :general-gaussian (->general-gaussian p sig) N (window-mode :general-gaussian opts true)))))

;;

(defn- parzen-exponential-fn
  ^double [^double alpha ^double r ^double x]
  (m/exp (m/- (m/pow (m/abs (m/* 2.0 alpha x)) r))))

(defn ->parzen-exponential
  "Returns the Parzen exponential window function `A exp(-abs(2 alpha x)^r)` (integral `1`, the constant by quadrature)."
  [^double alpha ^double r]
  (let [A (m// (integrate-half (partial parzen-exponential-fn alpha r)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-exponential-fn alpha r x))))))

(defn parzen-exponential
  "Parzen exponential family window `A exp(-abs(2 alpha x)^r)`.

  `r = 1` is `exponential` (same `alpha`), `r = 2` is `gaussian` with `alpha * sqrt(2)`. It is `general-gaussian` with `p = r / 2`.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`
      - `:r` - power, number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[exponential]], [[gaussian]], [[general-gaussian]], [[parzen-geometric]]."
  ([] (parzen-exponential 256))
  ([N] (parzen-exponential N nil))
  ([N opts]
   (let [alpha (param :parzen-exponential opts :alpha 1.0 open-positive)
         r (param :parzen-exponential opts :r 1.0 open-positive)]
     (continuous-window :parzen-exponential (->parzen-exponential alpha r) N (window-mode :parzen-exponential opts true)))))

;;

(defn- eta
  "dBc to linear"
  ^double [^double level]
  (m/exp10 (m// (m/abs level) 20.0)))

(def ^:private level-interval [-5000.0 true 5000.0 true])

(defn ->dolph-chebyshev
  "Returns a function of `N >= 2`, which returns the Dolph-Chebyshev window of `N` coefficients (mean `1`) with the sidelobe `level` in dB."
  [^double level]
  (let [e (eta level)
        ace (m/acosh e)]
    (fn [^long N]
      (let [N- (m/long-dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (poly/eval-chebyshev-T N- x0))
            dN (double N)]
        (->> (range N)
             (map (fn [^long n] (m/* f (poly/eval-chebyshev-T N- (m/* x0 (m/cospi (m// n dN)))))))
             (idft))))))

(defn dolph-chebyshev
  "Dolph-Chebyshev window: all sidelobes of the spectrum have the same level (equiripple).

  The window is the inverse DFT of the Chebyshev polynomial `T_(N-1)(x0 cos(pi k / N))` spectrum, `x0 = cosh(acosh(10^(abs(level)/20)) / (N - 1))`. Equal to `scipy.signal.windows.chebwin(N, at)` and Octave `chebwin(N, at)` with `level = -at`. The peak sidelobe level is the requested level (measured -30.00, -50.00, -80.00 and -100.00 dB for these levels, `N` = 2049). The cost is quadratic in `N` (about 10 ms for `N` = 1000).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:level` - sidelobe level in dB, finite number, `abs` is taken (so `50` and `-50` are the same), `abs(level) <= 5000`, default `-50.0`
      - `:normalize?` - default `true` (`:LInf`: the coefficients are divided by the largest absolute value, the maximum becomes `1`); `false` returns the raw window (mean `1`, sum `N`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  The accuracy is limited to about `1e-16 N^2` by the conditioning of the Chebyshev polynomial near its argument `1`: 2.5e-14 for `N` = 64, 4.7e-13 for `N` = 256 and 5e-11 for `N` = 2000 against a 60 digit reference. scipy and Octave have the same error, so windows of several hundred samples differ between the libraries by about `1e-13`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[ultraspherical]], [[saramaki]], [[legendre]], [[taylor]], [[kaiser]]."
  ([] (dolph-chebyshev 256))
  ([N] (dolph-chebyshev N nil))
  ([N opts]
   (discrete-window :dolph-chebyshev (->dolph-chebyshev (param :dolph-chebyshev opts :level -50.0 level-interval))
                    N (window-mode :dolph-chebyshev opts true) :LInf)))

;;

(defn- taylor-f-denominator
  ^double [^long n ^long m]
  (let [m2 (m/sq m)]
    (reduce (fn [^double p ^double n']
              (if (m/== n' m)
                p
                (m/* p (m/- 1.0 (m// m2 (m/* n' n')))))) 1.0 (range 1 n))))

(defn- taylor-f-numerator
  ^double [^double A2 ^double sigma2 ^long n ^long m]
  (let [msigma2 (m// (m/* m m) sigma2)        ]
    (reduce (fn [^double p ^double n']
              (m/* p (m/- 1.0 (m// msigma2 (m/+ A2 (m/sq (m/- n' 0.5))))))) 1.0 (range 1 n))))

(defn- taylor-f
  ^double [^double A2 ^double sigma2 ^long n ^long m]
  (let [s (if (m/even? m) -0.5 0.5)]
    (m// (m/* 2.0 s (taylor-f-numerator A2 sigma2 n m))
         (taylor-f-denominator n m))))

(defn ->taylor
  "Returns the Taylor window function (not clipped) with the design sidelobe `level` in dB and `n >= 2` nearly constant level sidelobes."
  [^double level ^long n]
  (let [e (eta level)
        A (m// (m/acosh e) m/PI)
        A2 (m/* A A)
        sigma2 (m// (m/* n n)
                    (m/+ A2 (m/sq (m/- n 0.5))))
        coeffs (vec (conj (map (partial taylor-f A2 sigma2 n) (range 1 n)) 1.0))]
    (->blackman-harris coeffs false)))

(defn taylor
  "Taylor window: a cosine sum with `n - 1` terms which approximates the Dolph-Chebyshev window, with the first `n - 1` sidelobes at nearly the design level and a decay beyond them (Doerry, eq. 208-212).

  Parameters:

  - `N` - number of samples, an integer `>= 0`. The window is the continuous function sampled at `-0.5, ..., 0.5` including both ends, which differs from the grid of `scipy.signal.windows.taylor` and Octave `taylorwin` (they sample `(n - N/2 + 1/2) / N`): the values differ from scipy by `O(1/N)` (0.45 for `N` = 2, 0.013 for `N` = 64, 0.003 for `N` = 256 for `n` = 4 and level -50).
  - `options` (optional map):
      - `:level` - design sidelobe level in dB (finite number, `abs` is taken, `abs(level) <= 5000`, default `-50.0`); it is `-sll` of scipy. It is the design parameter, the achieved peak sidelobe level depends on `n` and is higher for a small `n`. Measured (dB, `N` = 2049) for the levels -30, -50, -80, -100: `n` 2: -25.7, -30.3, -33.3, -34.0; `n` 4 (default): -30.4, -44.7, -56.0, -60.6; `n` 6: -30.3, -49.3, -66.9, -75.3; `n` 8: -30.2, -50.2, -72.7, -83.9. A level so small that the middle of the window is not positive (for example `0`, below about 20 dB the ends are higher than the middle) throws.
      - `:n` - integer in `[2, 400]`, the number of nearly constant level sidelobes (`nbar` of scipy), default `4`
      - `:normalize?` - default `true` (division by the value at the middle); `false` returns the raw samples of the function with the integral `1`; other values: `:LInf`, `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`. The cost of the coefficients is `O(n^2)`. Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?`, `:level` or `:n`.

  See also [[dolph-chebyshev]], [[ultraspherical]], [[blackman-harris-family]]."
  ([] (taylor 256))
  ([N] (taylor N nil))
  ([N opts]
   (let [level (param :taylor opts :level -50.0 level-interval)
         n (int-param :taylor opts :n 4 [2.0 true 400.0 true])
         f (->taylor level n)
         f0 (double (f 0.0))]
     (when-not (m/> f0 0.0)
       (invalid-input! :taylor :level
                       "too small, the window has a non-positive value in the middle"
                       (get opts :level) "level large enough for a positive window middle"))
     (continuous-window :taylor f N (window-mode :taylor opts true)))))

;;

(defn ->cauchy
  "Returns the Cauchy window function `A / (1 + (2 alpha x)^2)` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [A (m// alpha (m/atan alpha))]
    (fn ^double [^double x]
      (clamp x (m/* A (m// (m/inc (m/sq (m/* 2.0 alpha x)))))))))

(defn cauchy
  "Cauchy (Lorentzian) window `A / (1 + (2 alpha x)^2)`, `A = alpha / atan(alpha)`.

  Equal to `parzen-geometric` with `:r 2`. The ends are non-zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[parzen-geometric]], [[gaussian]], [[exponential]]."
  ([] (cauchy 256))
  ([N] (cauchy N nil))
  ([N opts]
   (continuous-window :cauchy (->cauchy (param :cauchy opts :alpha 1.0 open-positive)) N (window-mode :cauchy opts true))))

;;

(defn- parzen-geometric-fn
  ^double [^double alpha ^double r ^double x]
  (m// (m/inc (m/pow (m/abs (m/* 2.0 alpha x)) r))))

(defn ->parzen-geometric
  "Returns the Parzen geometric window function `A / (1 + abs(2 alpha x)^r)` (integral `1`, the constant by quadrature)."
  [^double alpha ^double r]
  (let [A (m// (integrate-half (partial parzen-geometric-fn alpha r)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-geometric-fn alpha r x))))))

(defn parzen-geometric
  "Parzen geometric family window `A / (1 + abs(2 alpha x)^r)`.

  `r = 2` is `cauchy` (same `alpha`).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`
      - `:r` - power, number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[cauchy]], [[parzen-exponential]], [[exponential]]."
  ([] (parzen-geometric 256))
  ([N] (parzen-geometric N nil))
  ([N opts]
   (let [alpha (param :parzen-geometric opts :alpha 1.0 open-positive)
         r (param :parzen-geometric opts :r 1.0 open-positive)]
     (continuous-window :parzen-geometric (->parzen-geometric alpha r) N (window-mode :parzen-geometric opts true)))))

;;

(defn ->kaiser-bessel
  "Returns the Kaiser-Bessel window function `A I0(pi alpha sqrt(1 - 4x^2))` (integral `1`), `alpha > 0`; `alpha = beta / pi` of the Kaiser window of scipy."
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        A (m// (m/* -2.0 pa) (m/expm1 (m/* -2.0 pa)))]
    (fn ^double [^double x]
      (clamp x (let [z (m/* pa (m/sqrt (m/- 1.0 (m/sq (m/* 2.0 x)))))]
                 (m/* A (m/exp (m/- z pa)) (bessel-i0e z)))))))

(defn kaiser-bessel
  "Kaiser-Bessel window `A I0(pi alpha sqrt(1 - 4x^2))`, an approximation of the DPSS (Slepian) window, with `I0` the modified Bessel function of the first kind and order 0.

  The window is [[kaiser]] with `beta = pi alpha`, so `alpha = beta / pi` (equal to `scipy.signal.windows.kaiser` and Octave `kaiser(N, pi alpha)`). The ends have the value `1 / I0(pi alpha)` of the middle. The peak sidelobe level is about -24.6 dB for the default `alpha = 1` and -45.9 dB for `alpha = 2`; a larger `alpha` gives lower sidelobes and a wider main lobe.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser]], [[kaiser-bessel-derived]], [[cosh]], [[knab]], [[avci-nacaroglu]], [[bessel-I1]]."
  ([] (kaiser-bessel 256))
  ([N] (kaiser-bessel N nil))
  ([N opts]
   (continuous-window :kaiser-bessel (->kaiser-bessel (param :kaiser-bessel opts :alpha 1.0 open-positive)) N (window-mode :kaiser-bessel opts true))))

(defn ->kaiser
  "Returns the Kaiser window function `A I0(beta sqrt(1 - 4x^2))` (integral `1`), `beta >= 0`; the same function as [[->kaiser-bessel]] with `alpha = beta / pi`, `beta = 0` is the rectangular window function."
  [^double beta]
  (if (m/zero? beta)
    rectangular-continuous
    (->kaiser-bessel (m// beta m/PI))))

(defn kaiser
  "Kaiser window `A I0(beta sqrt(1 - 4x^2))` for `x` in `[-0.5, 0.5]`, `I0` the modified Bessel function of the first kind and order 0; the window is a computable approximation of the DPSS (Slepian) window.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:beta` - shape, number in `[0, 1e300]`, default `pi`. A larger `beta` gives lower sidelobes and a wider main lobe. `beta = 0` is the rectangular window. It is the `beta` of `scipy.signal.windows.kaiser` and of Octave `kaiser`.
      - `:normalize?` - default `true` (division by the value at the middle, as scipy: the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`. The ends have the non-zero value `1 / I0(beta)` (relative to the middle).

  Throws `ex-info` for a `beta` which is negative, larger than `1e300`, `NaN` or infinite.

  The window is the same as [[kaiser-bessel]] with `:alpha` equal to `beta / pi`.

  See also [[kaiser-bessel]], [[kaiser-bessel-derived]], [[cosh]], [[knab]]."
  ([] (kaiser 256))
  ([N] (kaiser N nil))
  ([N opts]
   (continuous-window :kaiser (->kaiser (param :kaiser opts :beta m/PI [0.0 true 1.0e300 true])) N (window-mode :kaiser opts true))))

;;

(defn- cosh-ratio
  "`cosh(pa s) / cosh(pa)` without overflow."
  ^double [^double pa ^double s]
  (m// (m/+ (m/exp (m/* pa (m/dec s))) (m/exp (m/* (m/- pa) (m/inc s))))
       (m/inc (m/exp (m/* -2.0 pa)))))

(defn ->cosh
  "Returns the cosh window function `A cosh(pi alpha sqrt(1 - 4x^2)) / cosh(pi alpha)` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        A (m// (integrate-phi (fn [^double s] (m/* s (cosh-ratio pa s)))))]
    (fn ^double [^double x]
      (clamp x (m/* A (cosh-ratio pa (m/sqrt (m/- 1.0 (m/* 4.0 x x)))))))))

(defn cosh
  "Cosh window `A cosh(pi alpha s) / cosh(pi alpha)`, `s = sqrt(1 - 4x^2)`, an approximation of the DPSS (Slepian) window.

  It is the Kaiser-Bessel window with the Bessel function replaced by the hyperbolic cosine (Doerry, catalog).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `2.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel]], [[knab]], [[avci-nacaroglu]], [[bessel-I1]]."
  ([] (cosh 256))
  ([N] (cosh N nil))
  ([N opts]
   (continuous-window :cosh (->cosh (param :cosh opts :alpha 2.0 open-positive)) N (window-mode :cosh opts true))))

;;

(defn- avci-nacaroglu-fn
  ^double [^double alpha ^double x]
  (let [pa (m/* m/PI alpha)]
    (m/exp (m/* pa (m/dec (m/sqrt (m/- 1.0 (m/* 4.0 x x))))))))

(defn ->avci-nacaroglu
  "Returns the Avci-Nacaroglu window function `A exp(pi alpha (sqrt(1 - 4x^2) - 1))` (integral `1`, the constant by quadrature), `alpha > 0`."
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        A (m// (integrate-phi (fn [^double s] (m/* s (m/exp (m/* pa (m/dec s)))))))]
    (fn ^double [^double x]
      (clamp x (m/* A (avci-nacaroglu-fn alpha x))))))

(defn avci-nacaroglu
  "Avci-Nacaroglu window `A exp(pi alpha (s - 1))`, `s = sqrt(1 - 4x^2)`, an exponential approximation of the DPSS (Slepian) window.

  The ends are `exp(-pi alpha)` of the middle value.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `2.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel]], [[cosh]], [[knab]], [[bessel-I1]]."
  ([] (avci-nacaroglu 256))
  ([N] (avci-nacaroglu N nil))
  ([N opts]
   (continuous-window :avci-nacaroglu (->avci-nacaroglu (param :avci-nacaroglu opts :alpha 2.0 open-positive)) N (window-mode :avci-nacaroglu opts true))))

;;

(defn- knab-ratio
  "`sinh(pa s) / sinh(pa)` without overflow, `s >= 0`."
  ^double [^double pa ^double s]
  (m// (m/* (m/exp (m/* pa (m/dec s))) (m/expm1 (m/* -2.0 pa s)))
       (m/expm1 (m/* -2.0 pa))))

(defn- knab-fn
  "`sinh(pa s) / (sinh(pa) s)` with `s = sqrt(1 - 4x^2)`, equal to `pa / sinh(pa)` at `s = 0`."
  ^double [^double alpha ^double x]
  (let [pa (m/* m/PI alpha)
        s (m/sqrt (m/- 1.0 (m/* 4.0 x x)))]
    (if (m/zero? s)
      (m// (m/* -2.0 pa (m/exp (m/- pa))) (m/expm1 (m/* -2.0 pa)))
      (m// (knab-ratio pa s) s))))

(defn ->knab
  "Returns the Knab window function `A sinh(pi alpha s) / (sinh(pi alpha) s)`, `s = sqrt(1 - 4x^2)` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        A (m// (integrate-phi (fn [^double s] (knab-ratio pa s))))]
    (fn ^double [^double x]
      (clamp x (m/* A (knab-fn alpha x))))))

(defn knab
  "Knab window `A sinh(pi alpha s) / (sinh(pi alpha) s)`, `s = sqrt(1 - 4x^2)`, an approximation of the DPSS (Slepian) window.

  At the ends (`s = 0`) the function has the limit value `pi alpha / sinh(pi alpha)` of the middle value, not zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `2.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel]], [[cosh]], [[avci-nacaroglu]], [[bessel-I1]]."
  ([] (knab 256))
  ([N] (knab N nil))
  ([N opts]
   (continuous-window :knab (->knab (param :knab opts :alpha 2.0 open-positive)) N (window-mode :knab opts true))))

;;

(defn ->ultraspherical
  "Returns a function of `N >= 2`, which returns the ultraspherical window of `N` coefficients (mean `1`) for the parameter `level` (Doerry's `mu`) and `alpha > -1`; `alpha = 0` is the Dolph-Chebyshev window."
  [^double level ^double alpha]
  (let [e (eta level)
        ace (m/acosh e)
        cheb? (m/zero? alpha)
        poly-fn (fn ^double [^long n ^double x]
                  (if cheb?
                    (poly/eval-chebyshev-T n x)
                    (poly/eval-gegenbauer-C n alpha x)))]
    (fn [^long N]
      (let [N- (m/long-dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (double (poly-fn N- x0)))
            dN (double N)]
        (->> (range N)
             (mapv (fn [^long n] (m/* f (double (poly-fn N- (m/* x0 (m/cospi (m// n dN))))))))
             (idft))))))

(def ^:private ultraspherical-alpha-interval [-1.0 false 150.0 true])

(defn ultraspherical
  "Ultraspherical (Gegenbauer) window, a family between the Dolph-Chebyshev (`alpha = 0`) and the Saramaki (`alpha = 1`) windows.

  The window is the inverse DFT of the Gegenbauer polynomial spectrum `C_(N-1)^(alpha)(x0 cos(pi k / N))` (Doerry, catalog). `alpha = 0` is `dolph-chebyshev`, `alpha = 1` is `saramaki`, `alpha = 0.5` is the Legendre window (with a different level mapping, see [[legendre]]); for `alpha > 0` the sidelobes decay away from the main lobe (the Dolph-Chebyshev sidelobes stay at one level), faster for a larger `alpha`: for the level -50 and `N` = 513 the sidelobes between 40 and 60 bins are below -100 dB for `alpha` 2 and below -145 dB for `alpha` 4.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:level` - Doerry's sidelobe parameter `mu` in dB (finite number, `abs` is taken, `abs(level) <= 5000`, default `-50.0`). It is a design parameter, not the achieved peak sidelobe level. Measured achieved levels (dB, `N` = 2049, odd and even `N` agree within 0.1 dB) for the levels -30, -50, -80, -100: `alpha` 0.5: -30.1, -48.0, -76.1, -95.2; `alpha` 2 (default): -33.8, -47.0, -70.2, -86.7; `alpha` 4: -41.1, -50.5, -68.7, -82.6
      - `:alpha` - number in `(-1, 150]`, default `2.0`. Values above `150` overflow for `N` about `1000` (the limit gets lower for larger `N`); values close to `-1` are slow.
      - `:normalize?` - default `true` (`:LInf`: the coefficients are divided by the largest absolute value, the maximum becomes `1`); `false` returns the raw window (mean `1`, sum `N`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[dolph-chebyshev]], [[saramaki]], [[legendre]], [[taylor]]."
  ([] (ultraspherical 256))
  ([N] (ultraspherical N nil))
  ([N opts]
   (let [level (param :ultraspherical opts :level -50.0 level-interval)
         alpha (param :ultraspherical opts :alpha 2.0 ultraspherical-alpha-interval)]
     (discrete-window :ultraspherical (->ultraspherical level alpha) N (window-mode :ultraspherical opts true) :LInf))))

(defn saramaki
  "Saramaki window, the ultraspherical window with `alpha = 1`.

  Equal to `ultraspherical` with `:alpha 1.0`. As for all the ultraspherical windows the level is a design parameter and the achieved peak sidelobe level is higher than requested.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:level` - Doerry's sidelobe parameter in dB (finite number, `abs` is taken, `abs(level) <= 5000`, default `-50.0`). Measured achieved peak sidelobe levels (dB, `N` = 2049) for the levels -30, -50, -80, -100: -30.9, -47.1, -73.4, -91.6
      - `:normalize?` - default `true` (`:LInf`: the coefficients are divided by the largest absolute value, the maximum becomes `1`); `false` returns the raw window (mean `1`, sum `N`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[ultraspherical]], [[dolph-chebyshev]], [[legendre]]."
  ([] (saramaki 256))
  ([N] (saramaki N nil))
  ([N opts]
   (let [level (param :saramaki opts :level -50.0 level-interval)]
     (discrete-window :saramaki (->ultraspherical level 1.0) N (window-mode :saramaki opts true) :LInf))))

;;

(defn ->legendre
  "Returns a function of `N >= 2`, which returns the Legendre window of `N` coefficients (mean `1`) with the design `level` (Doerry, eq. 266: `1.0754 abs(level) + 1.7388` dB)."
  [^double level]
  (let [ace (m/acosh (m/exp10 (m// (m/+ (m/* 1.0754 (m/abs level)) 1.7388) 20.0)))]
    (fn [^long N]
      (let [N- (m/long-dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (poly/eval-legendre-P N- x0))
            dN (double N)]
        (->> (range N)
             (map (fn [^long n] (m/* f (poly/eval-legendre-P N- (m/* x0 (m/cospi (m// n dN)))))))
             (idft))))))

(defn legendre
  "Legendre window, the ultraspherical window with `alpha = 0.5` and Doerry's level mapping.

  The spectrum is the Legendre polynomial `P_(N-1)(x0 cos(pi k / N))` with the sidelobe parameter `1.0754 abs(level) + 1.7388` dB (Doerry, eq. 266), so `legendre` with `level` L is `ultraspherical` with `:alpha 0.5` and `:level` that value.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:level` - design sidelobe level in dB (finite number, `abs` is taken, `abs(level) <= 5000`, default `-50.0`). The achieved peak sidelobe level is about 3.4 dB lower than requested; measured (dB, `N` = 2049) for the levels -30, -50, -80, -100: -33.6, -53.1, -83.5, -104.1
      - `:normalize?` - default `true` (`:LInf`: the coefficients are divided by the largest absolute value, the maximum becomes `1`); `false` returns the raw window (mean `1`, sum `N`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[ultraspherical]], [[dolph-chebyshev]], [[saramaki]]."
  ([] (legendre 256))
  ([N] (legendre N nil))
  ([N opts]
   (discrete-window :legendre (->legendre (param :legendre opts :level -50.0 level-interval))
                    N (window-mode :legendre opts true) :LInf)))

;;

(defn ->bessel-I1
  "Returns the Bessel I1 window function `A I1(pi alpha s) / (I1(pi alpha) s)`, `s = sqrt(1 - 4x^2)` (integral `1`), `alpha > 0`."
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        c (m// (m/* 2.0 pa pa) (m/sq (m/- (m/expm1 (m/- pa)))))]
    (fn ^double [^double x]
      (clamp x (let [z (m/* pa (m/sqrt (m/- 1.0 (m/* 4.0 x x))))
                     q (if (m/zero? z) 0.5 (m// (bessel-i1e z) z))]
                 (m/* c q (m/exp (m/- z pa))))))))

(defn bessel-I1
  "Modified Bessel window of the first order `A I1(pi alpha s) / (I1(pi alpha) s)`, `s = sqrt(1 - 4x^2)`, an approximation of the DPSS (Slepian) window.

  `I1` is the modified Bessel function of the first kind and order 1; at the ends (`s = 0`) the function has the limit value `pi alpha / (2 I1(pi alpha))` of the middle value.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `2.0`
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel]], [[cosh]], [[knab]], [[avci-nacaroglu]]."
  ([] (bessel-I1 256))
  ([N] (bessel-I1 N nil))
  ([N opts]
   (continuous-window :bessel-I1 (->bessel-I1 (param :bessel-I1 opts :alpha 2.0 open-positive)) N (window-mode :bessel-I1 opts true))))

;;

(defn- shayesteh-kashtiban-n
  ^double [^long n ^long N-]
  (m/pow (m/sinc (m// (m/- n (m/* 0.5 N-))
                      (m/* 0.654 N-))) 2.5))


(defn shayesteh-kashtiban-discrete
  "Returns the coefficients of the Shayesteh-Kashtiban window of `N >= 2` samples as a vector, scaled to the mean `1` (Doerry, eq. 276, 277)."
  [^long N]
  (let [N- (m/long-dec N)
        endpoint (m/+ 0.02 (m/* 0.001 N-) (m// (m/+ 50.0 N- N-)))
        vs (mapv (fn [^long n]
                   (if (or (m/zero? n) (m/== n N-))
                     endpoint
                     (shayesteh-kashtiban-n n N-))) (range N))
        A (m// (double N) (v/sum vs))]
    (v/mult vs A)))

(defn shayesteh-kashtiban
  "Shayesteh-Kashtiban window: a sinc power window with enlarged end samples (Doerry, eq. 276, 277).

  The interior samples are `sinc((n - (N-1)/2) / (0.654 (N-1)))^2.5`; the two end samples are `0.02 + 0.001 (N-1) + 1 / (50 + 2 (N-1))`, before the common scale. The ends grow relative to the interior with `N`. The peak sidelobe level is about -47.9 dB (no sidelobe control).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true` (`:LInf`: the coefficients are divided by the largest absolute value, the maximum becomes `1`); `false` returns the raw window (mean `1`, sum `N`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[mottaghi-kashtiban-shayesteh]], [[dolph-chebyshev]], [[sinc]]."
  ([] (shayesteh-kashtiban 256))
  ([N] (shayesteh-kashtiban N nil))
  ([N opts]
   (discrete-window :shayesteh-kashtiban shayesteh-kashtiban-discrete N (window-mode :shayesteh-kashtiban opts true) :LInf)))

;;

(defn ->kaiser-bessel-derived
  "Returns a function of `N >= 2`, which returns the Kaiser-Bessel-derived window of `N` coefficients; `alpha = beta / pi`."
  [^double alpha]
  (fn [^long N]
    (let [kb (kaiser-bessel (long (m/inc (m/ceil (m// N 2.0)))) {:alpha alpha :normalize? false})
          total (double (v/sum kb))
          _ (when-not (m/pos? total)
              (invalid-input! :kaiser-bessel-derived :alpha
                              "too large, the kernel window has no non-zero sample"
                              alpha "alpha for which the Kaiser-Bessel kernel has a non-zero sample"))
          v (v/sqrt (v/div (butlast (reductions m/+ kb)) total))]
      (vec (concat v (if (m/odd? N) (rest (reverse v)) (reverse v)))))))

(defn kaiser-bessel-derived
  "Kaiser-Bessel-derived (KBD) window: the square root of the normalized cumulative sum of a Kaiser-Bessel window, power complementary.

  For even `N` the window satisfies `w[n]^2 + w[n + N/2]^2 = 1` (Princen-Bradley condition for the MDCT), equal to `scipy.signal.windows.kaiser_bessel_derived(N, pi alpha)`. For odd `N` the window is a symmetric extension, which is not power complementary (scipy rejects odd `N`).

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:alpha` - number `> 0`, default `1.0`; `alpha = beta / pi` of scipy. A value so large that the kernel window has no non-zero sample throws.
      - `:normalize?` - default `false`: the raw power complementary window; `true` is `:LInf` (the maximum becomes `1`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel]], [[kaiser]], [[vorbis]]."
  ([] (kaiser-bessel-derived 256))
  ([N] (kaiser-bessel-derived N nil))
  ([N opts]
   (discrete-window :kaiser-bessel-derived (->kaiser-bessel-derived (param :kaiser-bessel-derived opts :alpha 1.0 open-positive))
                    N (window-mode :kaiser-bessel-derived opts false) :LInf)))

;;

(defn vorbis-continuous
  "Vorbis window, continuous function `A sin(pi/2 cos^2(pi x))` on `[-0.5, 0.5]` (integral `1`), the sampled continuous model, not the discrete Vorbis window."
  ^double [^double x]
  (clamp x (m/* 1.660592492619163 (m/sin (m/* m/HALF_PI (m/sq (m/cospi x)))))))

(defn vorbis-discrete
  "Returns the coefficients of the Vorbis window of `N >= 2` samples as a sequence: `sin(pi/2 sin^2(pi (n + 1/2) / N))`, `n = 0, ..., N-1`."
  [^long N]
  (map (fn [^long n]
         (m/sin (m/* m/HALF_PI (m/sq (m/sin (m/* m/PI (m// (m/+ n 0.5) N)))))))
       (range N)))

(defn vorbis
  "Vorbis window (Xiph.org) `sin(pi/2 sin^2(pi (n + 1/2) / N))`, `n = 0, ..., N-1`: the power complementary window of the Vorbis audio codec.

  For even `N` the window satisfies `w[n]^2 + w[n + N/2]^2 = 1` (Princen-Bradley). It is defined on the half-sample grid (it is not a sample of the continuous function [[vorbis-continuous]] at `-0.5, ..., 0.5`); the ends are non-zero.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `false` (the values of the formula); `true` is `:LInf` (the maximum becomes `1`); other values: `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[kaiser-bessel-derived]], [[cosine]], [[hann]]."
  ([] (vorbis 256))
  ([N] (vorbis N nil))
  ([N opts]
   (discrete-window :vorbis vorbis-discrete N (window-mode :vorbis opts false) :LInf)))

;;

(defn cosine-continuous
  "Cosine window, continuous function `pi/2 cos(pi x)` on `[-0.5, 0.5]` (integral `1`). Sampled at the points `-0.5, ..., 0.5` it is [[power-of-cosine]] with `m = 1`, which is not the discrete [[cosine]] window."
  ^double [^double x]
  (clamp x (m/* m/HALF_PI (m/cospi x))))

(defn- cosine-discrete
  "Cosine window coefficients for `N >= 2`: `sin(pi (n + 1/2) / N)`, `n = 0, ..., N-1`."
  [^long N]
  (map (fn [^long n] (m/sin (m/* m/PI (m// (m/+ n 0.5) N)))) (range N)))

(defn cosine
  "Cosine (sine) window `sin(pi (n + 1/2) / N)`, `n = 0, ..., N-1`: the half-sample shifted cosine lobe, all values are positive and the ends are non-zero (`sin(pi / (2 N))`).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `false`, the values of the formula (the maximum is `1` for odd `N`, as in `scipy.signal.windows.cosine`, and `cos(pi / (2 N))` for even `N`); `true` is `:LInf` (the maximum becomes `1`). The other accepted values are `:L1`, `:L2`, `:N` and a finite non-zero number, see [[normalize-coefficients]].

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  This is a discrete window defined on the index (it is not a sample of a continuous function at the points `-0.5, ..., 0.5`), the raw scale is the formula, not the mean `1` scale. It is equal to `scipy.signal.windows.cosine` (symmetric and, through `fastmath.kernel/window` with `:symmetric? false`, periodic) and it is not `power-of-cosine` with `m = 1` (that one has zero ends, see [[cosine-continuous]]).

  See also [[cosine-continuous]], [[power-of-cosine]], [[vorbis]], [[hann]]."
  ([] (cosine 256))
  ([N] (cosine N nil))
  ([N opts]
   (discrete-window :cosine cosine-discrete N (window-mode :cosine opts false) :LInf)))

;;

(defn flat-top-continuous
  "Flat-top window of scipy, continuous function `A (0.21557895 + 0.41663158 cos(2 pi x) + 0.277263158 cos(4 pi x) + 0.083578947 cos(6 pi x) + 0.006947368 cos(8 pi x))` (integral `1`); it has negative values."
  ^double [^double x] (make-blackman-harris clamp- x 0.21557895 0.41663158 0.277263158 0.083578947 0.006947368))

(defn flat-top
  "Flat-top window of scipy and Matlab, a sum of five cosine terms for amplitude measurements.

  The window is `A (0.21557895 + 0.41663158 cos(2 pi x) + 0.277263158 cos(4 pi x) + 0.083578947 cos(6 pi x) + 0.006947368 cos(8 pi x))` (the centered form of the alternating sign formula). Equal to `scipy.signal.windows.flattop` and Matlab `flattopwin`: the coefficients sum to `1.000000003`, so the default `:normalize?` differs from scipy by `3e-9`. The window has negative lobes (not clipped) and a small jump at the ends; the peak sidelobe level is about -93.0 dB, the sidelobes drop as `1/f`. Octave `flattopwin` is another set, [[flat-top-bk]].

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[flat-top-3]], [[flat-top-bk]], [[flat-top-hp]], [[hft95]], [[blackman-harris-family]]."
  ([] (flat-top 256))
  ([N] (flat-top N nil))
  ([N opts]
   (continuous-window :flat-top flat-top-continuous N (window-mode :flat-top opts true))))

(defn flat-top-3-continuous
  "Flat-top window FTNI, continuous function `A (0.2810639 + 0.5208972 cos(2 pi x) + 0.1980399 cos(4 pi x))` (integral `1`); it has negative values."
  ^double [^double x] (make-blackman-harris clamp- x 0.2810639 0.5208972 0.1980399))

(defn flat-top-3
  "Flat-top window FTNI of the National Instruments software, a sum of three cosine terms.

  The window is `A (0.2810639 + 0.5208972 cos(2 pi x) + 0.1980399 cos(4 pi x))` (Heinzel et al. 2002, section D.2.1). The peak sidelobe level is about -44.4 dB, the sidelobes drop as `1/f`. The window has negative lobes (not clipped) and a jump at the ends.

  Parameters:

  - `N` - number of samples, an integer `>= 0`
  - `options` (optional map):
      - `:normalize?` - default `true`: the samples are divided by the value of the window function at the middle (the maximum is `1` for odd `N`, slightly lower for even `N`); `false` returns the raw samples of the continuous function, which has the integral `1`. Other accepted values: `:LInf` (alias `:Linf`, divides by the maximum), `:L1`, `:L2`, `:N` or a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  Returns a sequence of `N` doubles: the continuous function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included (`N = 1` gives the middle point). `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  Throws `ex-info` (`:type :invalid-input`) for an invalid `N`, `:normalize?` or parameter.

  See also [[flat-top]], [[flat-top-hp]], [[flat-top-bk]], [[blackman-harris-family]]."
  ([] (flat-top-3 256))
  ([N] (flat-top-3 N nil))
  ([N opts]
   (continuous-window :flat-top-3 flat-top-3-continuous N (window-mode :flat-top-3 opts true))))

;; flat-top windows of Heinzel, Ruediger and Schilling (2002), sections D.2 and D.3, and of Gade and Herlufsen (1987).
;; The sets are sums of cosines in the centered form with positive magnitudes (the report has alternating signs
;; and the window peak at j = N / 2). All of them have negative lobes and no clipping.

(defn flat-top-hp-continuous
  "Flat-top HP window (FTHP), continuous function `A (1 + 1.912510941 cos(2 pi x) + 1.079173272 cos(4 pi x) + 0.1832630879 cos(6 pi x))` with the integral `1`; it has negative values."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.912510941 1.079173272 0.1832630879))

(defn flat-top-hp
  "Flat-top window of the old Hewlett-Packard spectrum analyzers (FTHP of Heinzel et al. 2002, section D.2.2), a sum of four cosine terms with the positive magnitudes `1, 1.912510941, 1.079173272, 0.1832630879` (centered form).

  The peak sidelobe level is -70.4 dB (first zero at 4 bins, sidelobes drop as `1/f`), the normalized equivalent noise bandwidth is 3.4279 bins and the amplitude flatness error 0.0096 dB. The window has negative lobes and a jump at the ends (about `-0.004` of the middle value for the default `:normalize?`). The window has the same sidelobe level as [[hft70]], which has a smaller flatness error.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft70]], [[flat-top-srs]], [[flat-top-bk]], [[flat-top]], [[flat-top-3]], [[blackman-harris-family]]."
  ([] (flat-top-hp 256))
  ([N] (flat-top-hp N nil))
  ([N opts]
   (continuous-window :flat-top-hp flat-top-hp-continuous N (window-mode :flat-top-hp opts true))))

(defn flat-top-srs-continuous
  "Flat-top SRS window (FTSRS), continuous function `A (1 + 1.93 cos(2 pi x) + 1.29 cos(4 pi x) + 0.388 cos(6 pi x) + 0.028 cos(8 pi x))` with the integral `1`; it has negative values."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.93 1.29 0.388 0.028))

(defn flat-top-srs
  "Flat-top window of the Stanford Research SR785 spectrum analyzer (FTSRS of Heinzel et al. 2002, section D.2.4), a sum of five cosine terms with the positive magnitudes `1, 1.93, 1.29, 0.388, 0.028` (centered form).

  The peak sidelobe level is -76.6 dB (first zero at 4.72 bins, sidelobes drop as `1/f^3`), the normalized equivalent noise bandwidth is 3.7702 bins and the amplitude flatness error -0.0156 dB. The window is zero at the ends and has negative lobes. It differs from [[flat-top-bk]] only in the last coefficient.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[flat-top-bk]], [[flat-top-hp]], [[hft90d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (flat-top-srs 256))
  ([N] (flat-top-srs N nil))
  ([N opts]
   (continuous-window :flat-top-srs flat-top-srs-continuous N (window-mode :flat-top-srs opts true))))

(defn flat-top-bk-continuous
  "Flat-top BK window, continuous function `A (1 + 1.93 cos(2 pi x) + 1.29 cos(4 pi x) + 0.388 cos(6 pi x) + 0.0322 cos(8 pi x))` with the integral `1`; it has negative values."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.93 1.29 0.388 0.0322))

(defn flat-top-bk
  "Flat-top window of Gade and Herlufsen (1987, Bruel and Kjaer), a sum of five cosine terms with the positive magnitudes `1, 1.93, 1.29, 0.388, 0.0322` (centered form); with the default `:normalize?` it is equal to Octave `flattopwin` (the coefficients `1, 1.93, 1.29, 0.388, 0.0322` divided by their sum, the value of the middle).

  The peak sidelobe level is -68.3 dB (first zero at 5 bins), the normalized equivalent noise bandwidth is 3.7703 bins. The window has negative lobes and a small jump at the ends (about `0.0009` of the middle value). It is not the window of `scipy.signal.windows.flattop` (that one is [[flat-top]]).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[flat-top-srs]], [[flat-top]], [[flat-top-hp]], [[blackman-harris-family]]."
  ([] (flat-top-bk 256))
  ([N] (flat-top-bk N nil))
  ([N opts]
   (continuous-window :flat-top-bk flat-top-bk-continuous N (window-mode :flat-top-bk opts true))))

(defn hft70-continuous
  "HFT70 flat-top window, continuous function `A (1 + 1.90796 cos(2 pi x) + 1.07349 cos(4 pi x) + 0.18199 cos(6 pi x))` with the integral `1`; it has negative values."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.90796 1.07349 0.18199))

(defn hft70
  "HFT70 flat-top window (Heinzel et al. 2002, section D.3.1), a sum of four cosine terms with the positive magnitudes `1, 1.90796, 1.07349, 0.18199` (centered form), optimized for the lowest sidelobe level of three cosines.

  The peak sidelobe level is -70.4 dB (first zero at 4 bins, sidelobes drop as `1/f`), the normalized equivalent noise bandwidth is 3.4129 bins and the amplitude flatness error -0.0065 dB. The window has negative lobes and a jump at the ends (about `-0.004` of the middle value for the default `:normalize?`).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[flat-top-hp]], [[hft95]], [[hft90d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft70 256))
  ([N] (hft70 N nil))
  ([N opts]
   (continuous-window :hft70 hft70-continuous N (window-mode :hft70 opts true))))

(defn hft95-continuous
  "HFT95 flat-top window, continuous function `A (1 + 1.9383379 cos(2 pi x) + 1.3045202 cos(4 pi x) + 0.402827 cos(6 pi x) + 0.0350665 cos(8 pi x))` with the integral `1`; it has negative values."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.9383379 1.3045202 0.4028270 0.0350665))

(defn hft95
  "HFT95 flat-top window (Heinzel et al. 2002, section D.3.2), a sum of five cosine terms with the positive magnitudes `1, 1.9383379, 1.3045202, 0.402827, 0.0350665` (centered form), optimized for the lowest sidelobe level of four cosines. It is very close to the flat-top window of the HP/Agilent 35670A analyzers.

  The peak sidelobe level is -95.0 dB (first zero at 5 bins, sidelobes drop as `1/f`), the normalized equivalent noise bandwidth is 3.8112 bins and the amplitude flatness error 0.0044 dB. The window has negative lobes and a jump at the ends (about `-0.0003` of the middle value for the default `:normalize?`).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft90d]], [[hft70]], [[hft116d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft95 256))
  ([N] (hft95 N nil))
  ([N opts]
   (continuous-window :hft95 hft95-continuous N (window-mode :hft95 opts true))))

(defn hft90d-continuous
  "HFT90D flat-top window, continuous function `A (1 + 1.942604 cos(2 pi x) + 1.340318 cos(4 pi x) + 0.440811 cos(6 pi x) + 0.043097 cos(8 pi x))` with the integral `1`; it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.942604 1.340318 0.440811 0.043097))

(defn hft90d
  "HFT90D flat-top window (Heinzel et al. 2002, section D.3.3), a sum of five cosine terms with the positive magnitudes `1, 1.942604, 1.340318, 0.440811, 0.043097` (centered form). The coefficients satisfy `1 - a1 + a2 - a3 + a4 = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -90.2 dB (first zero at 5 bins), the normalized equivalent noise bandwidth is 3.8832 bins and the amplitude flatness error -0.0039 dB. The window has negative lobes.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft95]], [[hft116d]], [[flat-top-srs]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft90d 256))
  ([N] (hft90d N nil))
  ([N opts]
   (continuous-window :hft90d hft90d-continuous N (window-mode :hft90d opts true))))

(defn hft116d-continuous
  "HFT116D flat-top window, continuous function `A (1 + 1.9575375 cos(2 pi x) + 1.4780705 cos(4 pi x) + 0.6367431 cos(6 pi x) + 0.1228389 cos(8 pi x) + 0.0066288 cos(10 pi x))` with the integral `1`; it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.9575375 1.4780705 0.6367431 0.1228389 0.0066288))

(defn hft116d
  "HFT116D flat-top window (Heinzel et al. 2002, section D.3.4), a sum of six cosine terms with the positive magnitudes `1, 1.9575375, 1.4780705, 0.6367431, 0.1228389, 0.0066288` (centered form). The coefficients satisfy `1 - a1 + a2 - a3 + a4 - a5 = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -116.8 dB (first zero at 6 bins), the normalized equivalent noise bandwidth is 4.2186 bins and the amplitude flatness error -0.0028 dB. The window has negative lobes.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft90d]], [[hft144d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft116d 256))
  ([N] (hft116d N nil))
  ([N opts]
   (continuous-window :hft116d hft116d-continuous N (window-mode :hft116d opts true))))

(defn hft144d-continuous
  "HFT144D flat-top window, continuous function `A (1 + a1 cos(2 pi x) + ... + a6 cos(12 pi x))` with the integral `1` (coefficients see [[hft144d]]); it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.96760033 1.57983607 0.81123644 0.22583558 0.02773848 0.00090360))

(defn hft144d
  "HFT144D flat-top window (Heinzel et al. 2002, section D.3.5), a sum of seven cosine terms with the positive magnitudes `1, 1.96760033, 1.57983607, 0.81123644, 0.22583558, 0.02773848, 0.0009036` (centered form). The coefficients satisfy `1 - a1 + a2 - ... = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -144.1 dB (first zero at 7 bins), the normalized equivalent noise bandwidth is 4.5386 bins and the amplitude flatness error 0.0021 dB. The window has negative lobes.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft116d]], [[hft169d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft144d 256))
  ([N] (hft144d N nil))
  ([N opts]
   (continuous-window :hft144d hft144d-continuous N (window-mode :hft144d opts true))))

(defn hft169d-continuous
  "HFT169D flat-top window, continuous function `A (1 + a1 cos(2 pi x) + ... + a7 cos(14 pi x))` with the integral `1` (coefficients see [[hft169d]]); it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.97441842 1.65409888 0.95788186 0.33673418 0.06364621 0.00521942 0.00010599))

(defn hft169d
  "HFT169D flat-top window (Heinzel et al. 2002, section D.3.6), a sum of eight cosine terms with the positive magnitudes `1, 1.97441842, 1.65409888, 0.95788186, 0.33673418, 0.06364621, 0.00521942, 0.00010599` (centered form). The coefficients satisfy `1 - a1 + a2 - ... = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The fourth coefficient is `0.33673418`; the report prints `0.33673420`, for which the alternating sum is `2e-8` instead of 0 and the peak sidelobe level is only -167.9 dB. With `0.33673418` the sum is zero and the peak sidelobe level is -169.5 dB at 10.41 bins, as the report states.

  The first zero is at 8 bins, the normalized equivalent noise bandwidth is 4.8347 bins and the amplitude flatness error 0.0017 dB. The window has negative lobes.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft144d]], [[hft196d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft169d 256))
  ([N] (hft169d N nil))
  ([N opts]
   (continuous-window :hft169d hft169d-continuous N (window-mode :hft169d opts true))))

(defn hft196d-continuous
  "HFT196D flat-top window, continuous function `A (1 + a1 cos(2 pi x) + ... + a8 cos(16 pi x))` with the integral `1` (coefficients see [[hft196d]]); it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.979280420 1.710288951 1.081629853 0.448734314 0.112376628 0.015122992 0.000871252 0.000011896))

(defn hft196d
  "HFT196D flat-top window (Heinzel et al. 2002, section D.3.7), a sum of nine cosine terms with the positive magnitudes `1, 1.97928042, 1.710288951, 1.081629853, 0.448734314, 0.112376628, 0.015122992, 0.000871252, 0.000011896` (centered form). The coefficients satisfy `1 - a1 + a2 - ... = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -196.2 dB at 9.06 bins (first zero at 9 bins), the normalized equivalent noise bandwidth is 5.1134 bins and the amplitude flatness error 0.0013 dB. The window has negative lobes. Double precision samples reproduce the sidelobe level (measured -196.2 dB with an FFT of 8192 samples).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft169d]], [[hft223d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft196d 256))
  ([N] (hft196d N nil))
  ([N opts]
   (continuous-window :hft196d hft196d-continuous N (window-mode :hft196d opts true))))

(defn hft223d-continuous
  "HFT223D flat-top window, continuous function `A (1 + a1 cos(2 pi x) + ... + a9 cos(18 pi x))` with the integral `1` (coefficients see [[hft223d]]); it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.98298997309 1.75556083063 1.19037717712 0.56155440797 0.17296769663 0.03233247087 0.00324954578 0.00013801040 0.00000132725))

(defn hft223d
  "HFT223D flat-top window (Heinzel et al. 2002, section D.3.8), a sum of ten cosine terms with the positive magnitudes `1, 1.98298997309, 1.75556083063, 1.19037717712, 0.56155440797, 0.17296769663, 0.03233247087, 0.00324954578, 0.0001380104, 0.00000132725` (centered form). The coefficients satisfy `1 - a1 + a2 - ... = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -223.0 dB at 11.38 bins (first zero at 10 bins), the normalized equivalent noise bandwidth is 5.3888 bins and the amplitude flatness error -0.0011 dB. The window has negative lobes. Double precision samples reproduce the sidelobe level (measured -223.0 dB with an FFT of 8192 samples).

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft196d]], [[hft248d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft223d 256))
  ([N] (hft223d N nil))
  ([N opts]
   (continuous-window :hft223d hft223d-continuous N (window-mode :hft223d opts true))))

(defn hft248d-continuous
  "HFT248D flat-top window, continuous function `A (1 + a1 cos(2 pi x) + ... + a10 cos(20 pi x))` with the integral `1` (coefficients see [[hft248d]]); it has negative values and is zero at both ends."
  ^double [^double x]
  (make-blackman-harris clamp- x 1.0 1.985844164102 1.791176438506 1.282075284005 0.667777530266 0.240160796576 0.056656381764 0.008134974479 0.000624544650 0.000019808998 0.000000132974))

(defn hft248d
  "HFT248D flat-top window (Heinzel et al. 2002, section D.3.9), a sum of eleven cosine terms with the positive magnitudes `1, 1.985844164102, 1.791176438506, 1.282075284005, 0.667777530266, 0.240160796576, 0.056656381764, 0.008134974479, 0.00062454465, 0.000019808998, 0.000000132974` (centered form). The coefficients satisfy `1 - a1 + a2 - ... = 0`, so the window is zero at both ends and differentiable (the D in the name), the sidelobes drop as `1/f^3`.

  The peak sidelobe level is -248.4 dB at 13.37 bins (first zero at 11 bins), the normalized equivalent noise bandwidth is 5.6512 bins and the amplitude flatness error 0.0009 dB. The window has negative lobes. Double precision samples reproduce the sidelobe level (measured -248.4 dB with an FFT of 8192 samples), the sidelobes are below the precision of most signals.

  Parameters:

  - `N` - number of samples
  - `options`:
      - `:normalize?` - default `true` (division by the value at the middle, the maximum is `1` for odd `N`); `false` returns the continuous function sampled at `-0.5, ..., 0.5`, which has the integral `1`

  Returns a sequence of `N` doubles. `N = 0` returns `()` and `N = 1` returns `(1.0)`.

  See also [[hft223d]], [[hft196d]], [[flat-top]], [[blackman-harris-family]]."
  ([] (hft248d 256))
  ([N] (hft248d N nil))
  ([N opts]
   (continuous-window :hft248d hft248d-continuous N (window-mode :hft248d opts true))))

;;

(defn continuous->window-fn
  "Converts a continuous function of one variable on `[-0.5, 0.5]` to a window function with the same options as the windows of this namespace.

  Parameters:

  - `continuous-window-fn` - a function of one double, defined on `[-0.5, 0.5]` (it is called at the sample points and, for the default `:normalize?`, at `0.0`)

  Returns a function `(fn [] [N] [N options])` (`N` defaults to `256`). The returned window function takes `N`, an integer `>= 0`, and the map `options` with the key `:normalize?` (default `true`). It returns a sequence of `N` doubles, the function sampled at `N` evenly spaced points from `-0.5` to `0.5`, both ends included, as [[sample-window]] does: `true` divides the samples by the function value at `0.0` (a negative value flips the sign of the window), `false` returns the raw samples, the other accepted values are `:LInf` (alias `:Linf`), `:L1`, `:L2`, `:N` and a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default. The window function ignores `:symmetric?`; `fastmath.kernel/window` handles it.

  `N = 0` returns `()` and `N = 1` returns `(1.0)` (default normalization). Throws `ex-info` (`:type :invalid-input`) for an invalid `N` or `:normalize?` and when the value at `0.0` is zero or not finite and `:normalize?` is `true`.

  See also [[discrete->window-fn]], [[sample-window]]."
  [continuous-window-fn]
  (fn custom-continuous-window
    ([] (custom-continuous-window 256))
    ([N] (custom-continuous-window N nil))
    ([N opts]
     (continuous-window :custom-continuous-window continuous-window-fn N (window-mode :custom-continuous-window opts true)))))

(defn discrete->window-fn
  "Converts a function of the window length to a window function with the same options as the windows of this namespace.

  Parameters:

  - `discrete-window-fn` - a function of `N` (a long `>= 2`) returning a sequence of `N` numbers, the window coefficients

  Returns a function `(fn [] [N] [N options])` (`N` defaults to `256`). The returned window function takes `N`, an integer `>= 0`, and the map `options` with the key `:normalize?` (default `true`, which is `:LInf`: the coefficients are divided by the largest absolute value). `false` returns the coefficients of `discrete-window-fn` as they are (the type of the returned sequence is passed through), the other accepted values are `:L1`, `:L2`, `:N` and a finite non-zero number, see [[normalize-coefficients]]. `nil` is the default.

  `N = 0` returns `()` and `N = 1` returns `(1.0)` without calling `discrete-window-fn`. Throws `ex-info` (`:type :invalid-input`) for an invalid `N` or `:normalize?`. The window function ignores `:symmetric?`; `fastmath.kernel/window` handles it.

  See also [[continuous->window-fn]], [[normalize-coefficients]]."
  [discrete-window-fn]
  (fn custom-discrete-window
    ([] (custom-discrete-window 256))
    ([N] (custom-discrete-window N nil))
    ([N opts]
     (discrete-window :custom-discrete-window discrete-window-fn N (window-mode :custom-discrete-window opts true) :LInf))))
