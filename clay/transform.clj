^:kindly/hide-code
(ns transform
  (:require [fastmath.transform :as t]
            [fastmath.signal :as signal]
            [fastmath.signal.pad :as pad]
            [fastmath.signal.denoise :as denoise]
            [fastmath.signal.iir :as iir]
            [fastmath.signal.fir :as fir]
            [fastmath.signal.waveform :as wave]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]

            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.complex :as complex]
            [fastmath.kernel.window :as window]
            [fastmath.interpolation :as i]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.polynomials :as poly]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.transform.wavelets :as wv]))

;; # Transforms {.unnumbered}

;; ## Transformer factory

;; ::: {.callout-tip title="Defined functions"}
;; * `transformer`
;; :::

;; All transforms are created via a single multimethod `transformer`, dispatching on the transform family keyword and an optional subtype.

;; ### Forward and inverse transforms

;; ::: {.callout-tip title="Defined functions"}
;; * `forward-1d`, `reverse-1d`
;; * `forward-2d`, `reverse-2d`
;; :::

;; Collection of 1d and 2d transforms:

;; * Discrete Fourier Transform, real and complex
;; * Discrete Sine and Cosine Transforms
;; * Discrete Hadamard Transform
;; * Discrete Hartley Transform
;; * Discrete Wavelet Transform and Wavelet Packet Transform

;; All transforms can be created with `transformer` function and then used in `forward-` and `reverse-` functions.

(def fft-trans (t/transformer :real :fft))

(utls/examples-note
  (seq (t/forward-1d fft-trans [1 2 -10 1]))
  (seq (t/reverse-1d fft-trans [-6 -12 11 -1]))
  (m/double-double-array->seq (t/forward-2d fft-trans [[2 3] [-10 1]]))
  (m/double-double-array->seq (t/reverse-2d fft-trans [[-4 -12] [14 10]])))

;; Some of the transforms has its own dedicated function. For example `fft` and `ifft` for 1d DFT transform.

(utls/examples-note
  (t/fft [1 2 -10 1])
  (t/ifft [[-6 0] [11 -1] [-12 0]]))

;; See the table below with a list of all possible transforms.

(kind/table
 {:column-names ["transform" "type" "subtype" "comment"]
  :row-vectors [["Real DFT, FFT algorithm" (kind/md "`:real`") (kind/md "`:fft`") "Real input, complex output, half of the spectrum with Nyquist for even signal sizes."]
                ["Real DFT, naive implementation" (kind/md "`:real`") (kind/md "`:dft`") "Only 1d transform for signals of power of 2 sizes, half of the spectrum, no Nyquist component."]
                ["Complex DFT" (kind/md "`:complex`") (kind/md "`:fft`") "Complex input and output"]
                ["Complex DFT" (kind/md "`:complex`") (kind/md "`:rfft` or `:fftr`") "Real input, complex output (whole spectrum)"]
                ["Cosine" (kind/md "`:real`") (kind/md "`:dct`") "DCT-II/III"]
                ["Cosine" (kind/md "`:real`") (kind/md "`:cosine`") "DCT-I, signal size should be power of 2 plus 1."]
                ["Sine" (kind/md "`:real`") (kind/md "`:dst`") "DST-II/III"]
                ["Sine" (kind/md "`:real`") (kind/md "`:sine`") "DST-I, signal size should be power of 2, first element must be 0.0."]
                ["Hartley" (kind/md "`:real`") (kind/md "`:dht`") ""]
                ["Hadamard" (kind/md "`:real`") (kind/md "`:hadamard`") ""]
                ["Wavelet" (kind/md "`:dwt`") "name as a string" "Native implementation"]
                ["Wavelet" (kind/md "`:dwt`") "name as a keyword" "JWave implementation"]
                ["Wavelet Packet" (kind/md "`:wpt` or `:wpd`") "name as a string" "Native implementation"]
                ["Wavelet Packet" (kind/md "`:wpt` or `:wpd`") "name as a keyword" "JWave implementation"]]})

;; The following singal will be used as an illustration of coefficients.

(def signal (let [s1 (signal/waveform :sine {:f 28 :amplitude 12})
                s2 (signal/waveform :sine {:f 29 :amplitude 6})
                s3 (signal/waveform :sine {:f 117 :amplitude 8})
                s4 (signal/waveform :square {:f 50 :amplitude 2})]
            (-> (signal/add-waveforms s1 s2 s3 s4)
                (signal/sample-waveform 256 1.0))))

(gg/->image (gg/line (m/slice-range 0 1 256) signal {:xlab "Time (s)"
                                                     :ylab "Amplitude"
                                                     :title "Signal"}))

;; ## Fourier Transform (FFT)

;; Real and complex Fast Fourier Transform.

;; ::: {.callout-tip title="Defined functions"}
;; * `fft`, `ifft`
;; * `transformer` cases: `:real` `:fft`, `:complex` `:fftr`, `:complex` `:fft`
;; :::

;; ### `transformer`

;; When `transformer` is used directly we have following methods:

;; * `:real` `:fft` - real input, complex output, half of the spectrum, packed (see below)
;; * `:complex` `:fftr` - real input, complex output, full spectrum
;; * `:complex` `:fft` - complex input and output

;; Reverse transform also accepts `scale?` option (default `true`) to indicate if the result should be divided by number of elements or not.

;; ---

;; `:real` `:fft`

(def real-fft (t/transformer :real :fft))

;; This transofmer accepts real signal and returns half of the spectrum with the following layout:

;; **Even signal length** Nyquist component exists

(utls/examples-note
  (seq (t/forward-1d real-fft [2 1 4 -10 -2 0 2 0])))

;; * first element - real part of DC, here `-3.0`
;; * second element - real part of Nyquist, here `15.0`
;; * the rest of the result - interleaved real and imaginary parts from the first half of the spectrum, here: `11.78+4.36i`, `-6-11i` and `-3.78+8.36i`

;; **Odd signal length**  - Nyquist component doesn't exist

(utls/examples-note
  (seq (t/forward-1d real-fft [2 1 4 -10 -2 0 2])))

;; * first element - real part of DC, here `-3.0`
;; * second element - imaginary part of the last spectrum component, here `11.36i`
;; * the rest of the result - interleaved real and imaginary parts from the first half of the spectrum, last element is real part only and should be combined with second element, here: `13.79+0.35i`, `-9.75-3.54i` and real part `4.46`

;; For 2d signals, only sizes of power of two are accepted. Layout is complicated, please refer [JTransform](https://javadoc.io/static/com.github.wendykierp/JTransforms/3.2/org/jtransforms/fft/DoubleFFT_2D.html#realForward(double%5B%5D%5B%5D)) documentation.

(m/double-double-array->seq (t/forward-2d real-fft [[1 2 3 4]
                                                    [-1 0 -1 2]
                                                    [-1 1 1 -1]
                                                    [1 1 1 1]]))


;; ---

;; `:complex` `:fftr`

;; To get full spectrum without encoding DC and Nyquist components for real inputs we can use complex transformer

(def complex-fftr (t/transformer :complex :fftr))

(utls/examples-note
  (seq (t/forward-1d complex-fftr [2 1 4 -10 -2 0 2 0]))
  (seq (t/forward-1d complex-fftr [2 1 4 -10 -2 0 2])))

(m/double-double-array->seq (t/forward-2d complex-fftr [[1 2 3 4]
                                                        [-1 0 -1 2]
                                                        [-1 1 1 -1]
                                                        [1 1 1 1]]))

;; ---

;; `:complex` `:fft`

;; Complex transform for complex input. Layout for input should contain interleaved real and imaginary parts.

(def complex-fft (t/transformer :complex :fftr))

(utls/examples-note
  (seq (t/forward-1d complex-fft (mapcat identity [[1 2] [3 4] [10 -1]]))))

;; ### `fft`, `ifft`

;; `fft` and `ifft` are helper functions for 1d case. Result of `fft` is used in `fastmath.signal` functions.

;; For `fft`

;; * input can be any sequence of numbers (for `real` case) or pairs (for `complex` case)
;; * output is always a sequence of complex numbers, output sequence contains metadata which allows to perform proper reverse transform.
;; * there is one option: `spectrum` which indicates if returned sequences should contain whole spectrum (`:double-sided`, default for `complex` input) or half of the spectrum (`:single-sided`, default for `real` input)
;; * **NOTE:** real input with even data size always contains Nyquist component

;; See the following cases:

;; Real input, even data size, last element is Nyquist component

(let [result (t/fft [2 1 4 -10 -2 0 2 0])]
  {:meta (meta result)
   :result result})

;; Real input, odd data size, no Nyquist component

(let [result (t/fft [2 1 4 -10 -2 0 2])]
  {:meta (meta result)
   :result result})

;; Real input, full spectrum

(let [result (t/fft [2 1 4 -10 -2 0 2 0] {:spectrum :double-sided})]
  {:meta (meta result)
   :result result})

;; Complex input, full spectrum

(let [result (t/fft [[2 1] [4 -10] [-2 0] [2 0]])]
  {:meta (meta result)
   :result result})

;; For `ifft`:

;; * input is always a sequence of complex numbers (or pairs), preferably the result of `fft` function (which contains metadata about the type of the input).
;; * if input doesn't contain metadata, they can be provided as options:
;;   * `:kind` - the type of transform to perform: `:real` (default) or `:complex`.
;;   * `:even?` - for `:real` transforms, indicates if the original time-domain signal had an even length (crucial for correctly placing the Nyquist frequency).
;;   * `:real?` - for `:complex` transforms, indicates if the output should be narrowed to real numbers (doubles).
;; * additionaly, `:scale?` option indicates if result should be returned scaled (default: `true`) or not.

;; Reverse transform relies on meta data, however they can be overwriten when provieded in options. Let's construct result for odd input size.

(def fft-result-odd-signal (t/fft [2 1 4 -10 0]))

;; Reverse transform properly reads metadata and return original signal. However when metadata are stripped off, result of `ifft` is can be unexpected. To fix it, we can add proper information in options.

;; Reverse transform can be returned unscaled.

(utls/examples-note
  (t/ifft fft-result-odd-signal)
  (t/ifft (with-meta fft-result-odd-signal nil))
  (t/ifft (with-meta fft-result-odd-signal nil) {:even? false})
  (t/ifft fft-result-odd-signal {:scale? false}))

;; ---

;; Let's see coefficients of the example signal. Later in Signal Processing section we will show how to calculate amplitude, magnitude and other power spectrum values.

(take 10 (t/fft signal))

(let [signal-fft (t/fft signal)
      cnt (count signal-fft)
      r (i/interpolation :linear (range cnt) (map first signal-fft))
      c (i/interpolation :linear (range cnt) (map second signal-fft))]
  (gg/->image (gg/functions [["real" r]
                             ["imag" c]] {:x [0 cnt]
                                          :xlab "coeff number"
                                          :ylab "value"
                                          :title "FFT of the signal"})))

;; ## Discrete Cosine Transform (DCT)

;; DCT-I, DCT-II and DCT-III are implemented.

;; ::: {.callout-tip title="Defined functions"}
;; * `dct`, `idct`
;; * `transformer` cases: `:real` `:dct` (DCT-II/III), `:real` `:cosine` (DCT-I)
;; :::

;; ### `transformer`

;; When `transformer` is used directly we have following methods:

;; * `:real` `:dct` — DCT-II (forward) and DCT-III (inverse)
;; * `:real` `:cosine` — DCT-I; input should have size of power of 2 plus one ($N=2^p+1$)

;; Options:

;; * `:scale?`, for `:dct` — forward and reverse results can be scaled (default: `true`)
;; * `:normalization`, for `:cosine` — `:standard` or `:orthogonal` (default)

(def real-dct (t/transformer :real :dct))
(def real-cosine (t/transformer :real :cosine))

(utls/examples-note
  (seq (t/forward-1d real-dct [1 2 -1 -2 0]))
  (seq (t/forward-1d real-dct [1 2 -1 -2 0] {:scale? false}))
  (seq (t/forward-1d real-cosine [1 2 -1 -2 0]))
  (seq (t/forward-1d real-cosine [1 2 -1 -2 0] {:normalization :standard})))

;; DCT-II and DCT-III can be used as their own inverses mutually. DCT-I is it's own inverse (forward and reverse are the same functions)

(utls/examples-note
  (seq (t/forward-1d real-dct (t/reverse-1d real-dct [1 -2 3])))
  (seq (t/reverse-1d real-dct (t/forward-1d real-dct [1 -2 3])))
  (seq (t/forward-1d real-cosine [1 -2 3]))
  (seq (t/reverse-1d real-cosine [1 -2 3]))
  (seq (t/forward-1d real-cosine (t/forward-1d real-cosine [1 -2 3]))))

;; ### `dct`, `idct`

;; `dct` and `idct` are helper functions that work without constructing a transformer explicitly.

;; For `dct` options are:

;; * `:method` — `:DCT-I`, `:DCT-II` (default) or `:DCT-III`
;; * `:pad?` — for `:DCT-I`, left-pad with zeros so the size satisfies $N=2^p+1$ (default: `true`)
;; * `:scale?`

;; `idct` reads the `::dct` metadata attached by `dct` to determine the correct inverse type automatically.

;; ---

;; Let's see coefficients for all transforms for example signal.

(utls/examples-note
  (take 10 (t/dct signal {:method :DCT-I}))
  (take 10 (t/dct signal {:method :DCT-II}))
  (take 10 (t/dct signal {:method :DCT-III})))

^:kindly/hide-code
(defn- sine-cosine-coeffs-plot
  [f method]
  (let [coeffs (f signal {:method method})
        cnt (count coeffs)]
    (gg/line (range cnt) coeffs {:xlab "coeff number"
                                 :ylab "value"
                                 :title (str (name method) " of the signal.")})))

(kind/table
 [(->> [:DCT-I :DCT-II :DCT-III]
       (map (partial sine-cosine-coeffs-plot t/dct))
       (map gg/->image))])


;; ## Discrete Sine Transform (DST)

;; DST-I, DST-II and DST-III are implemented.

;; ::: {.callout-tip title="Defined functions"}
;; * `dst`, `idst`
;; * `transformer` cases: `:real` `:dst` (DST-II/III), `:real` `:sine` (DST-I)
;; :::

;; ### `transformer`

;; * `:real` `:dst` — DST-II (forward) and DST-III (inverse)
;; * `:real` `:sine` — DST-I; input must have even size ($2^p$) and first element must be `0.0`

;; Options: `:scale?` and `:normalization` (`:standard` or `:orthogonal`)

;; ### `dst`, `idst`

;; `dst` and `idst` are convenience helpers. Options for `dst`:

;; * `:method` — `:DST-I`, `:DST-II` (default) or `:DST-III`
;; * `:pad?` — for `:DST-I`, left-pad with a zero to satisfy the size requirement (default: `true`)
;; * `:scale?`

;; `idst` reads the `::dst` metadata attached by `dst` to choose the correct inverse automatically.

;; ---

;; Coefficients for all variants. Note that the DST-I result is longer due to left-padding with a zero.

(utls/examples-note
  (take 10 (t/dst signal {:method :DST-I}))
  (take 10 (t/dst signal {:method :DST-II}))
  (take 10 (t/dst signal {:method :DST-III})))

(kind/table
 [(->> [:DST-I :DST-II :DST-III]
       (map (partial sine-cosine-coeffs-plot t/dst))
       (map gg/->image))])

;; ## Discrete Hartley Transform (DHT)

;; ::: {.callout-tip title="Defined functions"}
;; * `dht`, `idht`
;; * `transformer` cases: `:real` `:dht`
;; :::

;; ### `transformer`

;; * `:real` `:dht` — Hartley transform (self-inverse up to scaling)

;; Options: `:scale?` (default `true`)

;; ### `dht`, `idht`

;; `dht` and `idht` are convenience helpers. `idht` accepts `:scale?` option.

(let [coeffs (t/dht signal)]
  (gg/->image (gg/line (range (count coeffs)) coeffs {:xlab "coeff number"
                                                      :ylab "value"
                                                      :title "DHT of the signal."})))

;; ## Fast Hadamard Transform (FHT)

;; ::: {.callout-tip title="Defined functions"}
;; * `hadamard`, `ihadamard`
;; * `transformer` cases: `:real` `:hadamard`
;; :::

;; ### `transformer`

;; * `:real` `:hadamard` — input length must be a power of 2; output is a `double-array`

;; ### `hadamard`, `ihadamard`

;; `hadamard` and `ihadamard` are convenience helpers.

(let [coeffs (t/hadamard signal)]
  (gg/->image (gg/line (range (count coeffs)) coeffs {:xlab "coeff number"
                                                      :ylab "value"
                                                      :title "Hadamard transform of the signal."})))


;; ## Wavelet Transforms

;; ::: {.callout-tip title="Defined functions"}
;; * `dwt`, `idwt`
;; * `wpt`, `iwpt`
;; * `wpd`
;; * `transformer` cases: `:dwt`, `:wpt`, `:wpd` + wavelet name (string = built-in, keyword = JWave)
;; :::

;; ### Wavelet families catalogue

;; There are two sources of wavelets: [JWave](https://github.com/graetz23/JWave) (keyword names) and a built-in Matlab-compatible library (string names).

;; `keyword` based names are for `JWave` wavelets, `string` based are for Matlab (built-in).

;; ::: {.callout-caution}
;; JWave wavelets have some issues with coefficients: some of them are not reversible or have reversed coefficients.
;; :::

;; ##### Haar

;; ::: {.callout-tip title="Wavelets"}
;; * `"haar"`, `"db1"`
;; * `:haar`, `:haar-orthogonal`
;; :::


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




(def sig (map (fn [id]
              (let [p (/ id 256.0)]
                (+ (m/* 12 (m/sin (* m/TWO_PI p 28)))
                   (m/* 6 (m/sin (* m/TWO_PI p 29)))
                   (m/* 8 (m/sin (* m/TWO_PI p 117)))
                   (m/* 2 (if (> (m/frac (m/* m/PI 50 p)) 0.5) -1.0 1.0))))) (range 256)))

(kind/table
 [[(gg/->image (gg/line (range) sig {:title "Signal"}))
   (let [s (signal/spectrum sig {:method :amplitude})]
     (gg/->image (gg/line (:freqs s) (:spectrum s)  {:title "FFT of the signal"})))]])


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

(utls/examples-note
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

;; ### Discrete Wavelet Transform (DWT)

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


;; ### Wavelet Packet Transform (WPT)

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

;; ### Wavelet Packet Decomposition (WPD)

;; `wpd` performs a full multi-level wavelet packet decomposition, returning coefficients for every level as a sequence of arrays.

(utls/examples-note
  (count (t/wpd chirp+sin11 "db2"))
  (count (first (t/wpd chirp+sin11 "db2"))))

;; ## Signal Padding

;; ::: {.callout-tip title="Defined functions"}
;; * `pad`
;; :::

;; Padding adds samples around a signal before transforms or analysis. `pad` supports eight methods (`:zero`, `:edge`, `:linear`, `:periodic`, `:symmetric`, `:antisymmetric`, `:reflect`, `:antireflect`) and three side options (`:left`, `:right`, `:both`).

(def partial-signal (m/sample m/sin 0.25 2.5 20))

(defn padding-chart
  [signal pad-method]
  (let [len (count signal)
        start (pad/pad-position :both 90 len)]
    (gg/scatters [[(name pad-method) (range) (signal/pad signal 90 pad-method)]
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
 [[(gg/->image (gg/dwt-scaleogram (t/forward-1d haar (signal/pad chirp+sin11 4096 :zero)) {:normalize? true :wavelet "HAAR, zero padding"}))
   (gg/->image (gg/dwt-scaleogram (t/forward-1d haar (signal/pad chirp+sin11 4096 :antireflect)) {:normalize? true :wavelet "HAAR, anti-reflect padding"}))]])



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


;; # Signal Processing {.unnumbered}

;; ## Utilities

;; ::: {.callout-tip title="Defined functions"}
;; * `db->linear`, `linear->db`
;; :::

;; Helpers for converting between linear amplitude and decibels.

(utls/examples-note
  (signal/db->linear 0.0)
  (signal/db->linear -6.0)
  (signal/linear->db 1.0)
  (signal/linear->db 0.5))

;; ## Waveforms and Oscillators

;; ::: {.callout-tip title="Defined functions"}
;; * `waveform`
;; * `add-waveforms`, `gain-waveform`
;; * `sample-waveform`
;; :::

^:kindly/hide-code
(defn wv->images
  ([nm] (wv->images nm nil))
  ([nm opts] (wv->images nm opts 64))
  ([nm opts fs]
   (let [opts (merge {:f 5} opts)
         wv (signal/waveform nm opts)
         s (signal/spectrum (signal/sample-waveform wv fs 4) {:method :amplitude :fs fs})]
     (kind/table
      [[(gg/->image (gg/function wv {:x [0 1] :title (str nm " signal, opts=" opts)}))
        (gg/->image (gg/lollipop (:freqs s) (:spectrum s)  {:title (str "FFT of the signal, fs=" fs)}))]]))))

^:kindly/hide-code
(defn ch->images
  ([nm] (ch->images nm nil))
  ([nm opts] (ch->images nm opts 64))
  ([nm opts fs]
   (let [opts (merge {:time 3 :f0 0 :f1 20} opts)
         wv (signal/chirp nm opts)
         s (signal/spectrum (signal/sample-waveform wv fs 3) {:method :amplitude :fs fs})]
     (kind/table
      [[(gg/->image (gg/function wv {:x [0 3] :title (str nm " chirp signal, opts=" opts)}))
        (gg/->image (gg/lollipop (:freqs s) (:spectrum s)  {:title (str "FFT of the signal, fs=" fs)}))]]))))

^:kindly/hide-code
(defn ->player
  ([f what] (->player f what nil))
  ([f what opts]
   (kind/audio {:samples (-> (f what (assoc opts :amplitude 0.5))
                             (signal/sample-waveform 10000.0 5))
                :sample-rate 10000.0})))

^:kindly/hide-code

;; ### Basic

(wv->images :sine)
(->player signal/waveform :sine {:f 440})

(wv->images :square)
(->player signal/waveform :square {:f 440})

(wv->images :square {:duty 0.8})
(->player signal/waveform :square {:f 440 :duty 0.8})

(wv->images :saw)
(->player signal/waveform :saw {:f 440})

(wv->images :saw {:up? false})
(->player signal/waveform :saw {:f 440 :up? false})

(wv->images :triangle)
(->player signal/waveform :triangle {:f 440})


;; ### Analog

(wv->images :analog-sine)
(->player signal/waveform :analog-sine {:f 440})

(wv->images :analog-saw)
(->player signal/waveform :analog-saw {:f 440})

(wv->images :analog-saw {:up? false})
(->player signal/waveform :analog-saw {:f 440 :up? false})

(wv->images :analog-triangle)
(->player signal/waveform :analog-triangle {:f 440})

;; ### Band-limited

(wv->images :band-limited-square)
(->player signal/waveform :band-limited-square {:f 440})

(wv->images :band-limited-square {:bands 50})
(->player signal/waveform :band-limited-square {:f 440 :bands 50})

(wv->images :band-limited-saw)
(->player signal/waveform :band-limited-saw {:f 440})

(wv->images :band-limited-saw {:bands 50})
(->player signal/waveform :band-limited-saw {:f 440 :bands 50})

(wv->images :band-limited-triangle)
(->player signal/waveform :band-limited-triangle {:f 440})

(wv->images :band-limited-triangle {:bands 50})
(->player signal/waveform :band-limited-triangle {:f 440 :bands 50})

;; ## Chirp Signals

;; ::: {.callout-tip title="Defined functions"}
;; * `chirp`
;; :::

;; A chirp is a sinusoid whose instantaneous frequency changes over time. Five sweep types are supported.

;; ### Linear

(ch->images :linear)

;; ### Quadratic up

(ch->images :quadratic-up)

;; ### Quadratic down

(ch->images :quadratic-down)

;; ### Logarithmic

(ch->images :logarithmic {:f0 1})

;; ### Hyperbolic

(ch->images :hyperbolic {:f0 1})

;; ## Spectral Analysis

;; ::: {.callout-tip title="Defined functions"}
;; * `spectrum`
;; * `fft-energy`, `fft-magnitude`, `fft-amplitude`, `fft-power`
;; * `fft-frequencies`
;; :::

;; `spectrum` is the main entry point. It calls `fft` internally and returns a map with `:freqs` and `:spectrum` keys. Available `:method` values:

;; * `:magnitude` — raw absolute values
;; * `:amplitude` (default) — normalized by N
;; * `:energy` — magnitude squared
;; * `:power` — normalized by N²
;; * `:psd` — Power Spectral Density
;; * `:asd` — Amplitude Spectral Density
;; * `:phase` — angle in radians
;; * `:real`, `:imag`, `:complex` — raw FFT components

;; `fft-frequencies` returns the frequency axis for N FFT bins at a given sampling rate `fs`.

;; ## Periodogram

;; ::: {.callout-tip title="Defined functions"}
;; * `periodogram`
;; :::

;; Welch's method: divide the signal into overlapping windowed frames, compute the FFT of each, and average the resulting power spectra. Options: `:method`, `:average` (`:mean`/`:median`/`:umedian`), `:window`, `:overlap`, `:fs`.

(def chirp-log (-> (signal/chirp :logarithmic {:f0 180 :f1 20 :time 30})
                 (signal/sample-waveform 400 30)))

;; first 100 samples
(gg/->image (gg/line (take 100 chirp-log)))

(def spectrum (-> chirp-log
                (t/fft)
                (signal/spectrum {:method :amplitude})
                (:spectrum)))

;; See the Spectral Analysis section above for all available `:method` options.

(def freqs (signal/fft-frequencies 400 (count chirp-log)))

;; first 100 coefficients
(gg/->image (gg/lollipop freqs (take 1000 spectrum)))

(def periodogram (-> chirp-log
                   (signal/periodogram {:fs 400
                                        :window (k/window :blackman-harris-92db 512)
                                        :overlap 0.75})))

(gg/->image (gg/line (:freqs periodogram) (:spectrum periodogram)))

;; ## Short-Time Fourier Transform (STFT)

;; ::: {.callout-tip title="Defined functions"}
;; * `stft`
;; :::

;; STFT slides a window over the signal and computes the FFT of each frame, producing a time–frequency map. Options: `:window`, `:overlap`, `:fs`, `:method`, `:db?`. The result map contains `:N`, `:spectrum`, `:freqs`, `:times`.

(def stft (-> chirp-log
            (signal/stft {:fs 400
                          :method :power
                          :window (k/window :blackman-harris-92db 512)
                          :overlap 0.75
                          :db? true})))

(-> stft
    (gg/spectrogram)
    (gg/->image))

(-> (signal/waveform :analog-triangle {:f 20 :phase m/sin})
    (signal/sample-waveform 400 1)
    (gg/line)
    (gg/->image))

(-> (signal/waveform :analog-triangle {:f 20 :phase m/sin})
    (signal/sample-waveform 400 30)
    #_(repeatedly (* 400 30) r/grand)
    (signal/stft {:fs 400
                  :method :power
                  :window (k/window :blackman-harris-92db 512)
                  :overlap 0.0
                  :db? true})
    (gg/spectrogram)
    (gg/->image))

(kind/audio {:samples (-> (signal/waveform :analog-triangle {:f 440 :phase #(m/sin (* 2000 %)) :amplitude 0.5})
                          (signal/sample-waveform 10000.0 5))
             :sample-rate 10000.0})

(v/mn (-> (signal/chirp :linear {:f1 50 :f0 440 :time 5 :amplitude 0.5})
          (signal/sample-waveform 10000.0 5)))

(v/mx (-> (signal/waveform :band-limited-square {:f 440 :amplitude 0.5})
          (signal/sample-waveform 10000.0 5)))

;; ## Convolution and Correlation

;; ::: {.callout-tip title="Defined functions"}
;; * `convolve`, `fft-convolve`
;; * `correlate`, `fft-correlate`
;; :::

;; Both direct (O(N·M)) and FFT-based (O(N log N)) implementations are provided. Output size is controlled by `:mode` — `:full` (default), `:same`, `:first`, `:last`.

(utls/examples-note
  (signal/convolve [1 2 3] [0 1 0.5])
  (signal/fft-convolve [1 2 3] [0 1 0.5])
  (signal/correlate [1 2 3] [1 0])
  (signal/fft-correlate [1 2 3] [1 0]))

;; ## Analytic Signal and Instantaneous Quantities

;; ::: {.callout-tip title="Defined functions"}
;; * `hilbert`
;; * `amplitude-envelope`
;; * `instantaneous-phase`, `instantaneous-frequency`
;; * `zero-phase`
;; * `resample`
;; :::

;; The Hilbert transform converts a real signal into its analytic (complex) form. The magnitude gives the amplitude envelope; the angle gives the instantaneous phase; its time derivative gives the instantaneous frequency.

;; `zero-phase` zeroes out the phase component in the FFT domain (keeps magnitude only). `resample` changes the signal length using FFT-based interpolation.

^:kindly/hide-code
(defn ->envelope
  "Wrap a waveform function with a slow amplitude modulation envelope."
  [sig]
  (fn [^double t]
    (m/* (double (sig t)) (m/inc (m/* 0.5 (m/sin (m/* m/TWO_PI 3.0 t)))))))

;; Instantaneous envelope, phase and frequency for a hyperbolic chirp with amplitude modulation:

(let [s (-> (signal/chirp :hyperbolic {:f0 20 :f1 100})
            (->envelope)
            (signal/sample-waveform 400))]
  [(gg/->image (gg/lines [["envelope" (range) (signal/amplitude-envelope s)]
                          ["signal"   (range) s]]) {:width 800})
   (gg/->image (gg/line (signal/instantaneous-phase s) {:title "Instantaneous phase"}))
   (gg/->image (gg/line (signal/instantaneous-frequency s 400) {:title "Instantaneous frequency"}))])

;; Same for the chirp+sin signal used in the DWT examples:

(let [s chirp+sin11]
  [(gg/->image (gg/lines [["envelope" (range) (signal/amplitude-envelope s)]
                          ["signal"   (range) s]]) {:width 800})
   (gg/->image (gg/line (signal/instantaneous-phase s) {:title "Instantaneous phase"}))
   (gg/->image (gg/line (signal/instantaneous-frequency s 400) {:title "Instantaneous frequency"}))])

;; ## Denoising

(require '[fastmath.signal.test-signals :as tsignal])

(def noisy-signal
  (v/add (signal/sample-waveform tsignal/blocks 8192)
         (r/->seq (r/distribution :normal {:sd 0.2}) 8192)))

(gg/->image (gg/line noisy-signal))

(def noisy-coeffs (t/dwt noisy-signal "db2"))

(gg/->image (gg/dwt-scaleogram noisy-coeffs {:wavelet "db2" :normalize? true}))


;; methods: `:hard`, `:soft`, `:garrote` and `:hyperbole`

;; * `:visu` - based on median absolute deviation estimate (default)
;; * `:universal` - based on standard deviation estimate
;; * `:sure` or `:rigrsure` - based on SURE estimator
;; * `:hybrid` or `:heursure` - hybrid SURE estimator
;; * `:avg` - abs coefficients average
;; * `:peaksavg` - mid point between min and max absolute value of coefficents
;; * `:topn` - keep top n/N largest coefficients
;; * `:minfdr`, `minfdr-robust` - keep top n/N coefficients based on p-values (robust - uses MAD for sigma estimation)

(denoise/threshold noisy-coeffs)

(let [denoised-coeffs (denoise/denoise noisy-coeffs {:method :garrote :thr :topn :top-n-ratio 0.02})]
  (kind/table
   [[(gg/->image (gg/dwt-scaleogram denoised-coeffs {:wavelet "db2" :normalize? true}))
     (gg/->image (gg/line (t/idwt denoised-coeffs "db2")))]]))

;; ## Filtering

;; ::: {.callout-tip title="Defined functions"}
;; * `filter-signal`, `filter-signal-1`
;; * `filter-filter-signal`
;; * `reset-IIR!`
;; :::

;; `filter-signal` applies either FIR coefficients (a `double-array`) or a stateful IIR cascade (`Cascade` from `fastmath.signal.iir`) to a full signal. `filter-signal-1` processes a single sample (or returns a curried function). `filter-filter-signal` applies the filter forward then backward for zero phase. `reset-IIR!` resets internal IIR state.

;; ### IIR filter frequency response

(let [size 1024]
  (gg/->image (gg/line (:spectrum (signal/spectrum
                                   (iir/response
                                    (iir/sos [[0.00734357  0.00193118  0.00734357  1.         -1.74382376
                                               0.82107003]
                                              [1.         -1.15125701  1.          1.         -1.62052328
                                               0.94522601]] 1) size)
                                   {:method :energy
                                    :kind :complex
                                    :db? true
                                    :domain :frequency})))))

;; ### FIR gammatone filter

(let [{:keys [freqs spectrum]} (signal/spectrum (signal/pad (fir/gammatone {:fs 16000 :cutoff 1000})
                                                            1024 :zero :right)
                                                {:method :energy
                                                 :db? true
                                                 :fs 16000})]
  (gg/->image (gg/line freqs spectrum {:xlog "log10"})))

;; ### Arbitrary-gain FIR filter

(let [{:keys [freqs spectrum]} (signal/spectrum (signal/pad
                                                 (fir/firgain
                                                  {:fs 20000
                                                   :taps 71
                                                   :cutoff [2500 4500 6000 7500]
                                                   :kind :bandstop
                                                   :window :kaiser-bessel-derived
                                                   :freq-gain-pairs [[0 0] [2500 25] [5000 0.0] [7000 1] [10000 0.0]]
                                                   :antisymmetric? false
                                                   :interpolator :step-before})
                                                 1024 :zero :right)
                                                {:method :energy
                                                 :db? true
                                                 :fs 20000})]
  (gg/->image (gg/line freqs spectrum)))

;; ## Smoothing

;; ::: {.callout-tip title="Defined functions"}
;; * `savgol-filter`
;; * `moving-average-filter`
;; * `kernel-smoothing-filter`
;; :::

;; Each function returns a filtering function `(fn [signal] ...)` that can be applied to any sequence.

;; ### Savitzky-Golay

;; `savgol-filter` fits a local polynomial of given `:order` inside a sliding window of `:window-size` samples. Optional `:derivative` returns the n-th derivative instead of the smoothed value.

;; ### Moving average

;; `moving-average-filter` computes the unweighted mean over a sliding window. Options: `:window-size`, `:side` (`:left`, `:right`, `:center`).

;; ### Kernel smoother

;; `kernel-smoothing-filter` applies Nadaraya-Watson kernel-weighted averaging. Options: `:kernel`, `:bandwidth`.

;; ## Effects and DSP Effects Chain

;; ::: {.callout-tip title="Defined functions"}
;; * `effect`, `effects-list`
;; * `compose-effects`, `reset-effects`
;; * `single-pass`, `apply-effects`
;; :::

;; DSP effects are stateful processors. Build a chain with `effect` and `compose-effects`, then apply it over a sequence with `apply-effects`.

;; ### Building an effects chain

(let [chain (signal/compose-effects
             #_(signal/effect :clipping {:method :soft})
             (signal/effect :svf3 {:kind :lowshelf :Q 0.5 :cutoff 500 :tan :exact :gaindb -20})
             #_(signal/effect :gain 5)
             #_(signal/effect :iir (iir/chebyshev {:ripple 15 :kind :lowpass :direct-form 2}))
             #_(signal/effect :simple-lowpass {:rate 22050 :cutoff 1000})
             #_(signal/effect :vcf303 {:rate 22050 :resonance 2 :gain 5}))
      signal (signal/sample-waveform (signal/waveform :saw {:f 440}) 22050 1)
      filtered (signal/apply-effects signal chain)]
  (kind/table
   [[(gg/->image (gg/line (take 1000 signal)))
     (gg/->image (gg/line (take 1000 filtered)))]
    [(kind/audio {:samples signal :sample-rate 22050.0})
     (kind/audio {:samples filtered :sample-rate 22050.0})]]))

;; ### All available effects

;; ::: {.callout-tip title="Effect keywords"}
;; * `:simple-lowpass`, `:simple-highpass`
;; * `:biquad-eq`, `:biquad-hs`, `:biquad-ls`
;; * `:biquad-lp`, `:biquad-hp`, `:biquad-bp`
;; * `:dj-eq`, `:phaser-allpass`
;; * `:divider`, `:fm`, `:bandwidth-limit`
;; * `:distort`, `:foverdrive`, `:decimator`
;; * `:basstreble`, `:echo`, `:vcf303`
;; * `:slew-limit`, `:mda-thru-zero`
;; :::

(utls/examples-note signal/effects-list)

;; ## File I/O

;; ::: {.callout-tip title="Defined functions"}
;; * `save-signal`, `load-signal`
;; :::

;; Signals are stored as 16-bit signed big-endian raw binary files.

(let [s (signal/sample-waveform (signal/waveform :sine {:f 440}) 8000 1)
      f (java.io.File/createTempFile "signal" ".raw")]
  (signal/save-signal s (.getPath f))
  (utls/examples-note
    (= (mapv int s) (mapv int (signal/load-signal (.getPath f))))))

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.transform)
(codox/make-public-fns-table-clay 'fastmath.signal)


;;

(gg/->image (gg/function #(wave/blamp % 0.2)))


(->> (wave/polyblep-hyptri2 {:fs 800 :f 1})
     (take 1600)
     (gg/line)
     (gg/->image))

(reduce m/min
        (->> (wave/polyblep-hyptri {:fs 2000 :f 1 :phase 0.4})
             (take 200000)))
