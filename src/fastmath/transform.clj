(ns fastmath.transform
  "Transforms.

  See [[transformer]] and [[TransformProto]] for details.
  
  ### Wavelet
  
  Based on [JWave](https://github.com/cscheiblich/JWave/) library.

  Be aware that some of the wavelet types doesn't work properly. `:battle-23`, `:cdf-53`, `:cdf-97`.

  ### Cos/Sin/Hadamard

  Orthogonal or standard fast sine/cosine/hadamard 1d transforms.

  ### Fourier

  DFT, FFT, DHT."
  (:require [fastmath.core :as m]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.vector :as v]
            [fastmath.signal.pad :as pad]
            [fastmath.transform.wavelets :as wv]

            [fastmath.signal.denoise :as denoise])
  (:import [jwave.transforms FastWaveletTransform WaveletPacketTransform AncientEgyptianDecomposition DiscreteFourierTransform]
           [jwave.compressions CompressorPeaksAverage CompressorMagnitude]
           [org.apache.commons.math3.transform FastSineTransformer FastCosineTransformer FastHadamardTransformer RealTransformer DstNormalization DctNormalization TransformType]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]
           [org.jtransforms.dht DoubleDHT_1D DoubleDHT_2D]
           [org.jtransforms.dct DoubleDCT_1D DoubleDCT_2D]
           [org.jtransforms.dst DoubleDST_1D DoubleDST_2D]
           [fastmath.java Array]
           [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(defmulti
  ^{:doc "Create transform object for given wavelet.

  #### Wavelets

  * `:dwt` or `:fast` for 1d or 2d Fast Wavelet Transform. Size of data should be power of `2`.
  * `:wpt` or `:packet` for 1d or 2d Wavelet Packet Transform. Size of data should be power of `2`.
  * `:decomposed-fast` for 1d Fast Wavelet Transform. Data can have any size (Ancient Egyptian Decomposition is used).
  * `:decomposed-packet` for 1d Wavelet Packet Transform. Data can have any size (Ancient Egyptian Decomposition is used).

  Second argument is wavelet name as key. See [[wavelets-list]] for all supported names.

  #### Sine/Cosine/Hadamard

  * `:standard` for 1d `:sine`, `:cosine`, `:hadamard`.
  * `:orthogonal` for 1d `:sine`, `:cosine`.

  Note that `:sine` and `:cosine` require first element to be equal `0`. Size of data should be power of 2.

  #### Fourier

  * `:standard` `:dft` - 1d Discrete Fourier Transform - returns double-array where even elements are real part, odd elements are imaginary part."}
  transformer (fn [t _] t))

;; deprecated, use :dwt
(defmethod transformer :fast [_ w] (transformer :dwt w))

(defmethod transformer :dwt [_ w] (if (keyword? w)
                                    (FastWaveletTransform. (wv/wavelet w))
                                    (wv/wavelet-reify w :dwt)))

(defmethod transformer :wpt [_ w] (if (keyword? w)
                                    (WaveletPacketTransform. (wv/wavelet w))
                                    (wv/wavelet-reify w :wpt)))

;; deprecaed use :wpt
(defmethod transformer :packet [_ w] (transformer :wpt w))


;;deprecated -> decomposed-dwt
(defmethod transformer :decomposed-fast [_ w] (AncientEgyptianDecomposition. (transformer :fast w)))
;;deprecated -> decomposed-wpt
(defmethod transformer :decomposed-packet [_ w] (AncientEgyptianDecomposition. (transformer :packet w)))

(defmethod transformer :decomposed-dwt [_ w] (AncientEgyptianDecomposition. (transformer :fast w)))
(defmethod transformer :decomposed-wpt [_ w] (AncientEgyptianDecomposition. (transformer :packet w)))

;; depracated -> :real
(defmethod transformer :standard [_ t] (case t
                                         :sine (FastSineTransformer. DstNormalization/STANDARD_DST_I)
                                         :cosine (FastCosineTransformer. DctNormalization/STANDARD_DCT_I)
                                         :hadamard (FastHadamardTransformer.)
                                         :dft (DiscreteFourierTransform.)))

;; deprecated -> :real
(defmethod transformer :orthogonal [_ t] (case t
                                           :sine (FastSineTransformer. DstNormalization/ORTHOGONAL_DST_I)
                                           :cosine (FastCosineTransformer. DctNormalization/ORTHOGONAL_DCT_I)))
(defmethod transformer :sine [_ t] (case t
                                     :orthogonal (FastSineTransformer. DstNormalization/ORTHOGONAL_DST_I)
                                     :standard (FastSineTransformer. DstNormalization/STANDARD_DST_I)))
(defmethod transformer :cosine [_ t] (case t
                                       :orthogonal (FastCosineTransformer. DctNormalization/ORTHOGONAL_DCT_I)
                                       :standard (FastCosineTransformer. DctNormalization/STANDARD_DCT_I)))

;; ACM

(defn- perform-sc-acm
  [kind xs normalization forward?]
  (let [^RealTransformer t (case kind
                             :sine (FastSineTransformer. (case normalization
                                                           :standard DstNormalization/STANDARD_DST_I
                                                           :orthogonal DstNormalization/ORTHOGONAL_DST_I))
                             :cosine (FastCosineTransformer. (case normalization
                                                               :standard DctNormalization/STANDARD_DCT_I
                                                               :orthogonal DctNormalization/ORTHOGONAL_DCT_I)))]
    (.transform t (m/seq->double-array xs)
                (if forward? TransformType/FORWARD TransformType/INVERSE))))

;; sine and cosine ACM
(defn sc-acm-reify
  [kind]
  (reify prot/TransformProto
    (forward-1d [t xs] (prot/forward-1d t xs nil))
    (forward-1d [_ xs {:keys [normalization] :or {normalization :orthogonal}}]
      (perform-sc-acm kind xs normalization true))
    (reverse-1d [t xs] (prot/reverse-1d t xs nil))
    (reverse-1d [_ xs {:keys [normalization] :or {normalization :orthogonal}}]
      (perform-sc-acm kind xs normalization false))))

;; for Hadamard
(extend RealTransformer
  prot/TransformProto
  {:forward-1d (fn [^RealTransformer t xs] (.transform t (m/seq->double-array xs) TransformType/FORWARD))
   :reverse-1d (fn [^RealTransformer t xs] (.transform t (m/seq->double-array xs) TransformType/INVERSE))})

;; JTransform

(defn- jt-forward-fft [xs _]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (count xs))
        out (double-array xs)]
    (.realForward t out)
    out))

(defn- jt-reverse-fft [xs scale?]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (count xs))
        out (double-array xs)]
    (.realInverse t out (boolean scale?))
    out))

(defn- jt-forward2-fft [xss _]
  (let [^DoubleFFT_2D t (DoubleFFT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.realForward t ^"[[D" out)
    out))

(defn- jt-reverse2-fft [xss scale?]
  (let [^DoubleFFT_2D t (DoubleFFT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.realInverse t ^"[[D" out (boolean scale?))
    out))

(defn- jt-forward-cfft [xs _]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (m// (count xs) 2))
        out (double-array xs)]
    (.complexForward t out)
    out))

(defn- jt-reverse-cfft [xs scale?]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (m// (count xs) 2))
        out (double-array xs)]
    (.complexInverse t out (boolean scale?))
    out))

(defn- jt-forward2-cfft [xss _]
  (let [s (m// (count (first xss)) 2)
        s2 (m/* s 2)
        ^DoubleFFT_2D t (DoubleFFT_2D. (count xss) s)
        out (into-array (map (fn [xs] (double-array (take s2 xs))) xss))]
    (.complexForward t ^"[[D" out)
    out))

(defn- jt-reverse2-cfft [xss scale?]
  (let [s (m// (count (first xss)) 2)
        s2 (m/* s 2)
        ^DoubleFFT_2D t (DoubleFFT_2D. (count xss) s)
        out (into-array (map (fn [xs] (double-array (take s2 xs))) xss))]
    (.complexInverse t ^"[[D" out (boolean scale?))
    out))

(defn- jt-forward-cfftr [xs _]
  (let [s (count xs)
        ^DoubleFFT_1D t (DoubleFFT_1D. s)
        in (double-array xs)
        out (double-array (* 2 s))]
    (System/arraycopy in 0 out 0 s)
    (.realForwardFull t out)
    out))

(defn- jt-forward2-cfftr [xss _]
  (let [r (count xss)
        s (count (first xss))
        s2 (m/* s 2)
        ^DoubleFFT_2D t (DoubleFFT_2D. r s)
        in (map double-array xss)
        out (repeatedly r #(double-array s2))]
    (doseq [[i o] (map vector in out)]
      (System/arraycopy i 0 o 0 s))
    (.realForwardFull t ^"[[D" (into-array out))
    out))

;; Hartley

(defn- jt-forward-dht [xs _]
  (let [^DoubleDHT_1D t (DoubleDHT_1D. (count xs))
        out (double-array xs)]
    (.forward t out)
    out))

(defn- jt-reverse-dht [xs scale?]
  (let [^DoubleDHT_1D t (DoubleDHT_1D. (count xs))
        out (double-array xs)]
    (.inverse t out (boolean scale?))
    out))

(defn- jt-forward2-dht [xss _]
  (let [^DoubleDHT_2D t (DoubleDHT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out)
    out))

(defn- jt-reverse2-dht [xss scale?]
  (let [^DoubleDHT_2D t (DoubleDHT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out (boolean scale?))
    out))

;; cosine II/III

(defn- jt-forward-dct [xs scale?]
  (let [^DoubleDCT_1D t (DoubleDCT_1D. (count xs))
        out (double-array xs)]
    (.forward t out (boolean scale?))
    out))

(defn- jt-reverse-dct [xs scale?]
  (let [^DoubleDCT_1D t (DoubleDCT_1D. (count xs))
        out (double-array xs)]
    (.inverse t out (boolean scale?))
    out))

(defn- jt-forward2-dct [xss scale?]
  (let [^DoubleDCT_2D t (DoubleDCT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out (boolean scale?))
    out))

(defn- jt-reverse2-dct [xss scale?]
  (let [^DoubleDCT_2D t (DoubleDCT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out (boolean scale?))
    out))

;; sine

(defn- jt-forward-dst [xs scale?]
  (let [^DoubleDST_1D t (DoubleDST_1D. (count xs))
        out (double-array xs)]
    (.forward t out (boolean scale?))
    out))

(defn- jt-reverse-dst [xs scale?]
  (let [^DoubleDST_1D t (DoubleDST_1D. (count xs))
        out (double-array xs)]
    (.inverse t out (boolean scale?))
    out))

(defn- jt-forward2-dst [xss scale?]
  (let [^DoubleDST_2D t (DoubleDST_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out (boolean scale?))
    out))

(defn- jt-reverse2-dst [xss scale?]
  (let [^DoubleDST_2D t (DoubleDST_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out (boolean scale?))
    out))

(defn- jt-reify
  [f r f2 r2]
  (reify prot/TransformProto
    (forward-1d [_ xs] (f xs true))
    (forward-1d [_ xs {:keys [scale?] :or {scale? true}}] (f xs scale?))
    (reverse-1d [_ xs] (r xs true))
    (reverse-1d [_ xs {:keys [scale?] :or {scale? true}}] (r xs scale?))
    (forward-2d [_ xss] (f2 xss true))
    (forward-2d [_ xss {:keys [scale?] :or {scale? true}}] (f2 xss scale?))
    (reverse-2d [_ xss] (r2 xss true))
    (reverse-2d [_ xss {:keys [scale?] :or {scale? true}}] (r2 xss scale?))))

(defmethod transformer :real [_ t]
  (case t
    :fft (jt-reify jt-forward-fft jt-reverse-fft jt-forward2-fft jt-reverse2-fft)
    :dht (jt-reify jt-forward-dht jt-reverse-dht jt-forward2-dht jt-reverse2-dht)
    :dct (jt-reify jt-forward-dct jt-reverse-dct jt-forward2-dct jt-reverse2-dct)
    :dst (jt-reify jt-forward-dst jt-reverse-dst jt-forward2-dst jt-reverse2-dst)
    :sine (sc-acm-reify :sine)
    :cosine (sc-acm-reify :cosine)
    :hadamard (FastHadamardTransformer.)
    :dft (DiscreteFourierTransform.)))

(defmethod transformer :complex [_ t]
  (case t
    :fft (jt-reify jt-forward-cfft jt-reverse-cfft jt-forward2-cfft jt-reverse2-cfft)
    :fftr (jt-reify jt-forward-cfftr jt-reverse-cfft jt-forward2-cfftr jt-reverse2-cfft)
    :rfft (transformer :complex :fftr)))

;;

(defn forward-1d
  "Forward transform of sequence or array."
  ([t xs] (prot/forward-1d t xs))
  ([t xs options] (prot/forward-1d t xs options)))

(defn reverse-1d
  "Forward transform of sequence or array."
  ([t xs] (prot/reverse-1d t xs))
  ([t xs options] (prot/reverse-1d t xs options)))

(defn forward-2d
  "Forward transform of sequence or array."
  [t xss] (prot/forward-2d t xss))

(defn reverse-2d
  "Forward transform of sequence or array."
  [t xss] (prot/reverse-2d t xss))


(set! *warn-on-reflection* false)
;; 1d or 2d unknown in the compilation time

(defn compress
  "Compress transformed signal `xs` with given magnitude `mag`."
  {:deprecated "Use `fastmath.signal.denoise/denoise` with `:avg` threshold."}
  ([trans xs ^double mag]
   (let [[fwd rev] (if (seqable? (first xs))
                     [prot/forward-2d prot/reverse-2d]
                     [prot/forward-1d prot/reverse-1d])]
     (->> xs
          (fwd trans)
          (.compress (CompressorMagnitude. mag))
          (rev trans))))
  ([xs ^double mag]
   (.compress (CompressorMagnitude. mag) xs)))

(defn compress-peaks-average
  "Compress transformed signal `xs` with peaks average as a magnitude"
  {:deprecated "Use `fastmath.signal.denoise/denoise` with `:peakavg` threshold."}
  ([trans xs]
   (let [[fwd rev] (if (seqable? (first xs))
                     [prot/forward-2d prot/reverse-2d]
                     [prot/forward-1d prot/reverse-1d])]
     (->> xs
          (fwd trans)
          (.compress (CompressorPeaksAverage.))
          (rev trans))))
  ([xs]
   (.compress (CompressorPeaksAverage.) xs)))

(set! *warn-on-reflection* true)

(def ^{:deprecated "Use `fastmath.signal.denoise/threshold`"} denoise-threshold denoise/threshold)
(def ^{:deprecated "Use `fastmath.signal.denoise/denoise`"} denoise denoise/denoise)

(m/unuse-primitive-operators)

(defn- rfft
  [xs]
  (let [t (transformer :real :fft)
        ^doubles txs (forward-1d t xs)
        len (alength txs)
        e? (m/even? len)
        len- (m/dec len)
        res (reduce (fn [buff ^long id]
                      (conj buff (Vec2. (Array/aget txs id)
                                        (Array/aget txs (m/inc id))))) [(Vec2. (Array/aget txs 0) 0.0)] (range 2 (if e? len len-) 2))]
    (with-meta (if (m/one? len) ;; single value
                 res
                 (conj res (if e?
                             (Vec2. (Array/aget txs 1) 0.0)
                             (Vec2. (Array/aget txs len-) (Array/aget txs 1)))))
      {::fft {:kind :real :nyquist? e?}})))

(defn- crfft
  ([xs] (crfft xs false))
  ([xs real?]
   (let [t (transformer :complex (if real? :rfft :fft))
         ^doubles txs (forward-1d t (if real? xs (mapcat identity xs)))
         len (alength txs)]
     (with-meta (reduce (fn [buff ^long id]
                          (conj buff (Vec2. (Array/aget txs id)
                                            (Array/aget txs (m/inc id))))) [] (range 0 len 2))
       {::fft {:kind :complex :real? real?}}))))

(defn fft
  "Compute the Fast Fourier Transform of a 1D signal.

  This function converts a time-domain signal into its frequency-domain representation. It automatically detects the input type: real-valued signals (sequence of numbers) are processed using an efficient Real FFT, while complex-valued signals (sequence of Complex numbers, see `fastmath.complex`) use a Complex FFT.

  For real signals, it defaults to a single-sided spectrum, exploiting Hermitian symmetry to save memory and computation.

  Input parameters:
  * `xs` - Input signal as a sequence of doubles (real) or a sequence of complex objects (`Vec2` type).
  * `options` - A map of configuration options:
      * `:spectrum` - Determines the output format for real-valued inputs. Options are `:single-sided` (default), which returns `(N/2)+1` coefficients, or `:double-sided`, which returns the full `N` length complex spectrum.

  Returns a sequence of `fastmath.vector.Vec2` representing complex coefficients in the frequency domain. The result includes metadata (e.g., `:kind`, `:nyquist?`) required by [[ifft]] to correctly perform the inverse transform.

  These are:

  * `:kind` - can be `:real` or `:complex`
  * `:nyquist?` - `true`, when real signal length was even and contained nyquist frequency
  * `:real?` - `true` when complex fft was performed

  Returned sequence contains Nyquist frequency coefficient for real and even signals."
  ([xs] (fft xs nil))
  ([xs {:keys [spectrum]
        :or {spectrum :single-sided}}]
   (if (number? (first xs))
     (if (= :double-sided spectrum)
       (crfft xs true)
       (rfft xs))
     (crfft xs))))

(defn ifft
  "Compute the Inverse Fast Fourier Transform (IFFT).

  Converts a frequency-domain signal (spectrum) back into its original time-domain representation. This function is the inverse operation of [[fft]] and automatically handles the reconstruction logic based on the type of the input spectrum. It utilizes metadata (such as signal symmetry and original length parity) attached to the input sequence by [[fft]] to ensure the resulting signal is restored with the correct dimensions and data type.

  Input parameters:
  * `xs` - A sequence of complex coefficients, typically as `fastmath.vector.Vec2` objects.
  * `options` - An optional map of configuration keys (usually inferred from `xs` metadata):
    * `:kind` - The type of transform to perform: `:real` (default) or `:complex`.
    * `:nyquist?` - For `:real` transforms, indicates if the original time-domain signal had an even length (crucial for correctly placing the Nyquist frequency).
    * `:real?` - For `:complex` transforms, indicates if the output should be narrowed to real numbers (doubles).
    * `:scale?` - For scaling the output, default: `true`.

  Output:
  Returns a sequence representing the time-domain signal. For `:real` kind, it returns a sequence of doubles. For `:complex` kind, it returns a sequence of `Vec2` (complex numbers) unless `:real?` is set to true."
  ([xs] (ifft xs nil))
  ([xs opts]
   (let [{:keys [kind real? nyquist?]
          :or {kind :real real? true nyquist? true}} (merge (::fft (meta xs)) opts)]
     (if (= :real kind)
       (let [t (transformer :real :fft)
             xs (if nyquist?
                  (let [[^double nyquist] (last xs)
                        ^doubles xs (double-array (mapcat identity (butlast xs)))]
                    (Array/aset xs 1 nyquist)
                    xs)
                  (let [xs (mapcat identity xs)
                        im (double (last xs))
                        ^doubles xs (double-array (butlast xs))]                    
                    (when-not (m/one? (alength xs)) (Array/aset xs 1 im)) ;; single value
                    xs))]
         (seq (reverse-1d t xs opts)))
       (let [t (transformer :complex (if real? :rfft :fft))
             ^doubles res (reverse-1d t (mapcat identity xs) opts)]
         (if real?
           (take-nth 2 res)
           (for [^long id (range 0 (alength res) 2)]
             (Vec2. (Array/aget res id)
                    (Array/aget res (m/inc id))))))))))

;; DCT/DST

(defn- maybe-pad-dct-i
  [xs]
  (let [cnt (count xs)]
    (if (m/power-of-two? (m/dec cnt))
      [xs 0]
      (let [nct (m/inc (m/round-up-pow2 cnt))]
        [(pad/zero (m/seq->double-array xs) nct :left) (m/- nct cnt)]))))

(defn dct
  "Compute the Discrete Cosine Transform (DCT) of a 1D signal.

  The DCT transforms a time-domain signal into the frequency domain using a sum of cosine functions. 

  Input parameters:
  * `xs` - Sequence of real numbers (doubles) representing the input signal.
  * `options` - A map of configuration keys:
      * `:method` - Specifies the DCT variant: `:DCT-I`, `:DCT-II` (standard forward DCT, default), or `:DCT-III` (often used as the inverse for DCT-II).
      * `:pad?` - Boolean (default `true`). When using `:DCT-I`, it automatically pads the signal with zeros to the nearest $2^k+1$ size if necessary.
      * `:scale?` - Boolean (default `true`). Indicates whether the transform should be normalized (for `:DCT-I` sets `:orthogonal` normalization when `true`)

  Returns a vector of doubles representing the frequency coefficients. The result includes metadata (under the `::dct` key) storing the `:method` and the number of `:padded` elements required to correctly reverse the transform using [[idct]]."
  ([xs] (dct xs nil))
  ([xs {:keys [method pad? scale?] :or {method :DCT-II pad? true scale? true} :as options}]
   (let [[xs padded] (if (and pad? (= :DCT-I method)) (maybe-pad-dct-i xs) [xs 0])
         t (case method
             :DCT-I #(forward-1d (transformer :real :cosine) % (assoc options :normalization
                                                                      (if scale? :orthogonal :standard)))
             :DCT-II #(forward-1d (transformer :real :dct) % options)
             :DCT-III #(reverse-1d (transformer :real :dct) % options))]
     (with-meta (vec (t xs)) {::dct {:method method
                                     :padded padded}}))))

(defn idct
  "Compute the Inverse Discrete Cosine Transform (IDCT) of a 1D signal.

  Converts frequency-domain coefficients back into the original time domain. This function is the inverse operation of [[dct]] and utilizes metadata (such as the specific transform method and any padding applied) attached to the input sequence to ensure the signal is restored with its original dimensions and scaling.

  Input parameters:
  * `xs` - Sequence of real numbers (frequency coefficients), typically the result of the [[dct]] function.
  * `options` - An optional map of configuration keys (usually inferred from `xs` metadata):
    * `:method` - The DCT variant was used for a transform: `:DCT-I`, `:DCT-II`, or `:DCT-III`.
    * `:padded` - The number of elements to drop from the beginning of the result to reverse padding added during the forward transform.
    * `:scale?` - Boolean (default `true`). Indicates whether the transform should be normalized. For `:DCT-I` sets `:orthogonal` normalization when `true`.

  Returns a sequence of doubles representing the reconstructed signal in the time/spatial domain."
  ([xs] (idct xs nil))
  ([xs options]
   (let [{:keys [method padded scale?] :or {method :DCT-II padded 0 scale? true} :as options} (merge (::dct (meta xs)) options)
         t (case method
             :DCT-I #(reverse-1d (transformer :real :cosine) % (assoc options :normalization
                                                                      (if scale? :orthogonal :standard)))
             :DCT-II #(reverse-1d (transformer :real :dct) % options)
             :DCT-III #(forward-1d (transformer :real :dct) % options))]
     (drop padded (t xs)))))

;;

(defn- maybe-pad-dst-i
  ([xs]
   (if-not (m/zero? (double (first xs)))
     (maybe-pad-dst-i (cons 0.0 xs) 1)
     (maybe-pad-dst-i xs 0)))
  ([xs ^long padded]
   (let [cnt (count xs)]
     (if (m/power-of-two? cnt)
       [xs padded]
       (let [nct (m/round-up-pow2 cnt)]
         [(pad/zero (m/seq->double-array xs) nct :left) (m/+ padded (m/- nct cnt))])))))

(defn dst
  "Compute the Discrete Sine Transform (DST) of a 1D signal.

  The DST transforms a time-domain signal into the frequency domain using a sum of sine functions. 

  Input parameters:
  * `xs` - Sequence of real numbers (doubles) representing the input signal.
  * `options` - A map of configuration keys:
      * `:method` - Specifies the DST variant: `:DST-I`, `:DST-II` (standard forward DST, default), or `:DST-III` (often used as the inverse for DST-II).
      * `:pad?` - Boolean (default `true`). When using `:DST-I`, it automatically pads the signal with zeros to the nearest $2^k$ size if necessary, also prepends with 0.0 if necessary.
      * `:scale?` - Boolean (default `true`). Indicates whether the transform should be normalized (for `:DST-I` sets `:orthogonal` normalization when `true`)

  Returns a vector of doubles representing the frequency coefficients. The result includes metadata (under the `::dst` key) storing the `:method` and the number of `:padded` elements required to correctly reverse the transform using [[idst]]."
  ([xs] (dst xs nil))
  ([xs {:keys [method pad? scale?] :or {method :DST-II pad? true scale? true} :as options}]
   (let [[xs padded] (if (and pad? (= :DST-I method)) (maybe-pad-dst-i xs) [xs 0])
         t (case method
             :DST-I #(forward-1d (transformer :real :sine) % (assoc options :normalization
                                                                    (if scale? :orthogonal :standard)))
             :DST-II #(forward-1d (transformer :real :dst) % options)
             :DST-III #(reverse-1d (transformer :real :dst) % options))]
     (with-meta (vec (t xs)) {::dst {:method method
                                     :padded padded}}))))

(defn idst
  "Compute the Inverse Discrete Sine Transform (IDST) of a 1D signal.

  Converts frequency-domain coefficients back into the original time domain. This function is the inverse operation of [[dst]] and utilizes metadata (such as the specific transform method and any padding applied) attached to the input sequence to ensure the signal is restored with its original dimensions and scaling.

  Input parameters:
  * `xs` - Sequence of real numbers (frequency coefficients), typically the result of the [[dst]] function.
  * `options` - An optional map of configuration keys (usually inferred from `xs` metadata):
    * `:method` - The DST variant was used for a transform: `:DST-I`, `:DST-II`, or `:DST-III`.
    * `:padded` - The number of elements to drop from the beginning of the result to reverse padding added during the forward transform.
    * `:scale?` - Boolean (default `true`). Indicates whether the transform should be normalized. For `:DST-I` sets `:orthogonal` normalization when `true`.

  Returns a sequence of doubles representing the reconstructed signal in the time/spatial domain."
  ([xs] (idst xs nil))
  ([xs options]
   (let [{:keys [method padded scale?] :or {method :DST-II padded 0 scale? true} :as options} (merge (::dst (meta xs)) options)
         t (case method
             :DST-I #(reverse-1d (transformer :real :sine) % (assoc options :normalization
                                                                    (if scale? :orthogonal :standard)))
             :DST-II #(reverse-1d (transformer :real :dst) % options)
             :DST-III #(forward-1d (transformer :real :dst) % options))]
     (drop padded (t xs)))))

;; Hartley

(defn dht
  "Compute the Discrete Hartley Transform (DHT) of a 1D signal.

  Input parameters:
  * `xs` - Input signal as a sequence or array of real numbers (doubles).

  Returns a sequence of doubles representing the Hartley coefficients in the frequency domain."
  [xs]
  (seq (forward-1d (transformer :real :dht) xs)))

(defn idht
  "Compute the Inverse Discrete Hartley Transform (IDHT) of a 1D signal.

  Input parameters:
  * `xs` - Sequence or array of real numbers representing the Hartley coefficients.
  * `options` - An optional map of configuration keys:
      * `:scale?` - Boolean (default `true`). When true, the resulting signal is divided by the signal length $N$ to ensure the inverse mapping returns the data to its original scale.

  Returns a sequence of doubles representing the reconstructed signal in the time/spatial domain."
  ([xs] (idht xs nil))
  ([xs options]
   (seq (reverse-1d (transformer :real :dht) xs options))))

;; Hadamard

(defn hadamard
  "Compute the Fast Hadamard Transform (FHT) of a 1D signal.

  Input parameters:
  * `xs` - Input signal as a sequence or array of real numbers (doubles). The length of the input must be a power of 2.

  Returns a double array of Hadamard coefficients."
  [xs]
  (forward-1d (transformer :real :hadamard) xs))

(defn ihadamard
  "Compute the Inverse Fast Hadamard Transform (IFHT) of a 1D signal.

  Input parameters:
  * `xs` - Sequence or array of real numbers (doubles) representing the Hadamard coefficients. The length of the input must be a power of 2.

  Returns a sequence of doubles representing the reconstructed signal."
  [xs]
  (seq (reverse-1d (transformer :real :hadamard) xs)))


;; Wavelets

(defn- maybe-wavelet-pad
  [xs pad?]
  (let [len (count xs)
        po2? (m/power-of-two? len)]
    (if-not (or pad? po2?)
      (throw (ex-info "Length of the signal should be power of 2." {:length len}))
      (if po2?
        [xs 0]
        (let [nlen (m/round-up-pow2 len)]
          [(pad/zero (m/seq->double-array xs) nlen :left) (m/- nlen len)])))))

(defn dwt
  "Compute the Discrete Wavelet Transform (DWT) of a 1D signal.

  Performs a multi-resolution analysis by decomposing the input signal into approximation and detail coefficients. This allows for the analysis of signal features at different scales and positions, providing a time-frequency representation where the frequency resolution is high for low-frequency components and time resolution is high for high-frequency components.

  Input parameters:
  * `xs` - Input signal as a sequence. The length must typically be a power of 2.
  * `wavelet` - The wavelet to use for decomposition. Can be a string for built-in wavelets (e.g., \"haar\", \"db4\", \"sym5\") or a keyword for JWave-based wavelets (e.g., :haar, :daubechies-4). See [[fastmath.transform.wavelets/wavelet-names]] for all supported names.
  * `options` - A map of configuration keys:
      * `:level` - The number of decomposition levels to perform.
      * `:decompose?` - Boolean (default `false`). If `true`, the resulting flat coefficient array is restructured into a sequence of bands (approximation followed by details from coarsest to finest) using [[fastmath.transform.wavelets/decompose-dwt]]
      * `:pad?` - Boolean (default `true`). When true, automatically pads the signal with zeros to the nearest power of 2 size..

  Returns a double array of wavelet coefficients. If `:decompose?` is set to `true`, returns a sequence of sequences representing the structured decomposition levels. The result includes metadata (under the `::dwt` key) storing the number of `:padded` elements required to correctly reverse the transform using [[idwt]]."
  ([xs wavelet] (dwt xs wavelet nil))
  ([xs wavelet {:keys [level decompose? pad?]
                :or {pad? true}
                :as opts}]
   (let [[xs padded] (maybe-wavelet-pad xs pad?)
         t (transformer :dwt wavelet)
         coeffs (forward-1d t xs opts)]
     (with-meta (if-not decompose?
                  (vec coeffs)
                  (wv/decompose-dwt coeffs level))
       {::dwt {:padded padded}}))))

(defn idwt
  "Compute the Inverse Discrete Wavelet Transform (IDWT) to reconstruct a signal from its wavelet coefficients.

  Reverses the multi-resolution decomposition performed by [[dwt]].

  Input parameters:
  * `coeffs` - A sequence or array of wavelet coefficients. It supports both flat representations (concatenated coefficients) and structured sequences of bands (nested sequences as returned by [[dwt]] with the `:decompose?` option).
  * `wavelet` - The wavelet used for reconstruction. This must match the wavelet used in the forward transform. Can be a string for built-in wavelets (e.g., \"haar\", \"db4\") or a keyword for JWave-based wavelets (e.g., :haar, :daubechies-4).
  * `options` - An optional map of configuration keys:
      * `:level` - The number of reconstruction levels to perform.
      * `:padded` - The number of leading elements to drop from the result (usually inferred from `::dwt` metadata attached to `coeffs` by the [[dwt]] function).

  Returns a sequence of doubles representing the reconstructed signal in the time domain."
  ([coeffs wavelet] (idwt coeffs wavelet nil))
  ([coeffs wavelet options]
   (let [{:keys [^long padded] :or {padded 0} :as options} (merge (::dwt (meta coeffs)) options)
         t (transformer :dwt wavelet)
         coeffs (if (sequential? (first coeffs)) (flatten coeffs) coeffs)]
     (drop padded (reverse-1d t coeffs options)))))

(defn wpt
  "Compute the Wavelet Packet Transform (WPT) of a 1D signal.

  Performs a full decomposition of the signal by recursively applying low-pass and high-pass filters to both approximation and detail coefficients at each level.

  Input parameters:
  * `xs` - Input signal as a sequence or array. The length must be a power of 2.
  * `wavelet` - The wavelet to use for decomposition. Can be a string for built-in wavelets (e.g., \"haar\", \"db4\") or a keyword for JWave-based wavelets (e.g., :haar, :daubechies-4). See [[fastmath.transform.wavelets/wavelet-names]] for all supported names.
  * `options` - A map of configuration keys:
      * `:level` - The number of decomposition levels to perform.
      * `:decompose?` - Boolean (default `false`). If `true`, the resulting flat coefficient array is restructured into a sequence of $2^{level}$ equal-width frequency bands using [[fastmath.transform.wavelets/decompose-wpt]].
      * `:pad?` - Boolean (default `true`). When `true`, automatically pads the signal with zeros to the nearest power of 2.

  Returns a double array of wavelet packet coefficients. If `:decompose?` is `true`, it returns a sequence of sequences representing the structured frequency sub-bands. The result includes metadata (under the `::wpt` key) containing the number of `:padded` elements to facilitate correct reconstruction via [[iwpt]]."
  ([xs wavelet] (wpt xs wavelet nil))
  ([xs wavelet {:keys [level decompose? pad?]
                :or {pad? true}
                :as options}]
   (let [[xs padded] (maybe-wavelet-pad xs pad?)
         t (transformer :wpt wavelet)
         coeffs (forward-1d t xs options)]
     (with-meta (if-not decompose?
                  (vec coeffs)
                  (wv/decompose-wpt coeffs level))
       {::wpt {:padded padded}}))))

(defn iwpt
  "Compute the Inverse Wavelet Packet Transform (IWPT) to reconstruct a signal.

  Reverses the full decomposition tree performed by [[wpt]], transforming wavelet packet coefficients back into the time domain.

  Input parameters:
  * `coeffs` - A sequence or array of wavelet packet coefficients. It supports both flat representations and structured sequences of frequency bands (as returned by [[wpt]] with the `:decompose?` option).
  * `wavelet` - The wavelet used for reconstruction. This must match the wavelet used in the forward transform. Can be a string (e.g., \"db4\") or a keyword (e.g., :daubechies-4). See [[fastmath.transform.wavelets/wavelet-names]] for supported names.
  * `options` - An optional map of configuration keys:
    * `:level` - The number of decomposition levels to perform during reconstruction.
    * `:padded` - The number of leading elements to drop from the result to reverse padding added during the forward transform (usually inferred from metadata attached to `coeffs`).

  Returns a sequence of doubles representing the reconstructed signal in the time domain."
  ([coeffs wavelet] (iwpt coeffs wavelet nil))
  ([coeffs wavelet options]
   (let [{:keys [^long padded] :or {padded 0} :as options} (merge (::dwt (meta coeffs)) options)
         t (transformer :wpt wavelet)
         coeffs (if (sequential? (first coeffs)) (flatten coeffs) coeffs)]
     (drop padded (reverse-1d t coeffs options)))))

(defn wpd
  "Compute the full Wavelet Packet Decomposition (WPD) for all levels of a 1D signal.

  This function performs a comprehensive decomposition by recursively applying high-pass and low-pass filters to both approximation and detail coefficients at every possible scale. While [[wpt]] typically targets a specific level, `wpd` returns the coefficients for every level from 0 (the original or padded signal) to the maximum depth, providing a complete hierarchical view of the signal's frequency structure across all resolutions.

  Input parameters:
  * `xs` - Input signal as a sequence or array of real numbers.
  * `wavelet` - The wavelet used for decomposition. Can be a string for built-in wavelets (e.g., \"db4\") or a keyword for JWave-based wavelets (e.g., :daubechies-4).
  * `options` - A map of configuration keys:
      * `:decompose?` - Boolean (default `true`). If `true`, restructures the flat coefficient array of each level into a sequence of equal-width frequency bands using [[fastmath.transform.wavelets/decompose-wpt]].
      * `:pad?` - Boolean (default `true`). When `true`, automatically pads the signal with zeros to the nearest power of 2.

  Returns a sequence of levels, where each level contains the corresponding wavelet packet coefficients. If `:decompose?` is `true`, each level is returned as a sequence of sequences representing frequency sub-bands. The result includes metadata (under the `::wpd` key) storing the number of `:padded` elements used during the transform."
  ([xs wavelet] (wpd xs wavelet nil))
  ([xs wavelet {:keys [decompose? pad?]
                :or {decompose? true pad? true}}]
   (let [[xs padded] (maybe-wavelet-pad xs pad?)
         arr (m/seq->double-array xs)
         levels (range 1 (m/inc (m/log2int (alength arr))))
         wv (wv/wavelet wavelet)
         res (->> levels
                  (reduce (fn [[^long prev-level buff] ^long level]
                            [level (conj buff (wv/wpt-forward-1d wv (last buff) prev-level level))]) [0 [arr]])
                  (second))]
     (with-meta (if-not decompose?
                  res
                  (mapv (fn [coeffs ^long level]
                          (wv/decompose-wpt coeffs level)) res (conj levels 0)))
       {::wpd {:padded padded}}))))
