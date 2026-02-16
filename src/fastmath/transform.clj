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
            [fastmath.stats :as stat]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.optimization :as optim]
            [fastmath.vector :as v]

            [fastmath.transform.pad :as pad]
            [fastmath.transform.wavelets :as wv])
  (:import [jwave.transforms FastWaveletTransform WaveletPacketTransform AncientEgyptianDecomposition
            BasicTransform DiscreteFourierTransform]
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

(defmethod transformer :fast [_ w] (transformer :dwt w))
(defmethod transformer :dwt [_ w] (if (keyword? w)
                                    (FastWaveletTransform. (wv/wavelet w))
                                    (wv/wavelet-reify w :dwt)))
(defmethod transformer :packet [_ w] (transformer :wpt w))
(defmethod transformer :wpd [_ w] (transformer :wpt w))
(defmethod transformer :wpt [_ w] (if (keyword? w)
                                    (WaveletPacketTransform. (wv/wavelet w))
                                    (wv/wavelet-reify w :wpt)))

(defmethod transformer :decomposed-fast [_ w] (AncientEgyptianDecomposition. (transformer :fast w)))
(defmethod transformer :decomposed-packet [_ w] (AncientEgyptianDecomposition. (transformer :packet w)))

(defmethod transformer :standard [_ t] (case t
                                         :sine (FastSineTransformer. DstNormalization/STANDARD_DST_I)
                                         :cosine (FastCosineTransformer. DctNormalization/STANDARD_DCT_I)
                                         :hadamard (FastHadamardTransformer.)
                                         :dft (DiscreteFourierTransform.)))


(defmethod transformer :orthogonal [_ t] (case t
                                           :sine (FastSineTransformer. DstNormalization/ORTHOGONAL_DST_I)
                                           :cosine (FastCosineTransformer. DctNormalization/ORTHOGONAL_DCT_I)))

(extend BasicTransform
  prot/TransformProto
  {:forward-1d (fn ([^BasicTransform t xs] (.forward t (m/seq->double-array xs)))
                 ([^BasicTransform t xs {:keys [^long level]}] (.forward t (m/seq->double-array xs) level)))
   :reverse-1d (fn ([^BasicTransform t xs] (.reverse t (m/seq->double-array xs)))
                 ([^BasicTransform t xs {:keys [^long level]}] (.reverse t (m/seq->double-array xs) level)))
   :forward-2d (fn [^BasicTransform t xss] (.forward t (m/seq->double-double-array xss)))
   :reverse-2d (fn [^BasicTransform t xss] (.reverse t (m/seq->double-double-array xss)))})

(extend RealTransformer
  prot/TransformProto
  {:forward-1d (fn [^RealTransformer t xs] (.transform t (m/seq->double-array xs) TransformType/FORWARD))
   :reverse-1d (fn [^RealTransformer t xs] (.transform t (m/seq->double-array xs) TransformType/INVERSE))})

;; jtransform

(defn- jt-forward-fft [xs]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (count xs))
        out (double-array xs)]
    (.realForward t out)
    out))

(defn- jt-reverse-fft [xs]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (count xs))
        out (double-array xs)]
    (.realInverse t out true)
    out))

(defn- jt-forward2-fft [xss]
  (let [^DoubleFFT_2D t (DoubleFFT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.realForward t ^"[[D" out)
    out))

(defn- jt-reverse2-fft [xss]
  (let [^DoubleFFT_2D t (DoubleFFT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.realInverse t ^"[[D" out true)
    out))

(defn- jt-forward-cfft [xs]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (m// (count xs) 2))
        out (double-array xs)]
    (.complexForward t out)
    out))

(defn- jt-reverse-cfft [xs]
  (let [^DoubleFFT_1D t (DoubleFFT_1D. (m// (count xs) 2))
        out (double-array xs)]
    (.complexInverse t out true)
    out))

(defn- jt-forward2-cfft [xss]
  (let [s (m// (count (first xss)) 2)
        s2 (m/* s 2)
        ^DoubleFFT_2D t (DoubleFFT_2D. (count xss) s)
        out (into-array (map (fn [xs] (double-array (take s2 xs))) xss))]
    (.complexForward t ^"[[D" out)
    out))

(defn- jt-reverse2-cfft [xss]
  (let [s (m// (count (first xss)) 2)
        s2 (m/* s 2)
        ^DoubleFFT_2D t (DoubleFFT_2D. (count xss) s)
        out (into-array (map (fn [xs] (double-array (take s2 xs))) xss))]
    (.complexInverse t ^"[[D" out true)
    out))

(defn- jt-forward-cfftr [xs]
  (let [s (count xs)
        ^DoubleFFT_1D t (DoubleFFT_1D. s)
        in (double-array xs)
        out (double-array (* 2 s))]
    (System/arraycopy in 0 out 0 s)
    (.realForwardFull t out)
    out))

(defn- jt-forward2-cfftr [xss]
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

(defn- jt-forward-dht [xs]
  (let [^DoubleDHT_1D t (DoubleDHT_1D. (count xs))
        out (double-array xs)]
    (.forward t out)
    out))

(defn- jt-reverse-dht [xs]
  (let [^DoubleDHT_1D t (DoubleDHT_1D. (count xs))
        out (double-array xs)]
    (.inverse t out true)
    out))

(defn- jt-forward2-dht [xss]
  (let [^DoubleDHT_2D t (DoubleDHT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out)
    out))

(defn- jt-reverse2-dht [xss]
  (let [^DoubleDHT_2D t (DoubleDHT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out true)
    out))

(defn- jt-forward-dct [xs]
  (let [^DoubleDCT_1D t (DoubleDCT_1D. (count xs))
        out (double-array xs)]
    (.forward t out true)
    out))

(defn- jt-reverse-dct [xs]
  (let [^DoubleDCT_1D t (DoubleDCT_1D. (count xs))
        out (double-array xs)]
    (.inverse t out true)
    out))

(defn- jt-forward2-dct [xss]
  (let [^DoubleDCT_2D t (DoubleDCT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out true)
    out))

(defn- jt-reverse2-dct [xss]
  (let [^DoubleDCT_2D t (DoubleDCT_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out true)
    out))

(defn- jt-forward-dst [xs]
  (let [^DoubleDST_1D t (DoubleDST_1D. (count xs))
        out (double-array xs)]
    (.forward t out true)
    out))

(defn- jt-reverse-dst [xs]
  (let [^DoubleDST_1D t (DoubleDST_1D. (count xs))
        out (double-array xs)]
    (.inverse t out true)
    out))

(defn- jt-forward2-dst [xss]
  (let [^DoubleDST_2D t (DoubleDST_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.forward t ^"[[D" out true)
    out))

(defn- jt-reverse2-dst [xss]
  (let [^DoubleDST_2D t (DoubleDST_2D. (count xss) (count (first xss)))
        out (into-array (map double-array xss))]
    (.inverse t ^"[[D" out true)
    out))

(defn- jt-reify
  [f r f2 r2]
  (reify prot/TransformProto
    (forward-1d [_ xs] (f xs))
    (reverse-1d [_ xs] (r xs))
    (forward-2d [_ xss] (f2 xss))
    (reverse-2d [_ xss] (r2 xss))))

(defmethod transformer :real [_ t]
  (case t
    :fft (jt-reify jt-forward-fft jt-reverse-fft jt-forward2-fft jt-reverse2-fft)
    :dht (jt-reify jt-forward-dht jt-reverse-dht jt-forward2-dht jt-reverse2-dht)
    :dct (jt-reify jt-forward-dct jt-reverse-dct jt-forward2-dct jt-reverse2-dct)
    :dst (jt-reify jt-forward-dst jt-reverse-dst jt-forward2-dst jt-reverse2-dst)
    :sine (FastSineTransformer. DstNormalization/STANDARD_DST_I)
    :cosine (FastCosineTransformer. DctNormalization/STANDARD_DCT_I)
    :hadamard (FastHadamardTransformer.)
    :dft (DiscreteFourierTransform.)))

(defmethod transformer :complex [_ t]
  (case t
    :fft (jt-reify jt-forward-cfft jt-reverse-cfft jt-forward2-cfft jt-reverse2-cfft)
    :fftr (jt-reify jt-forward-cfftr jt-reverse-cfft jt-forward2-cfftr jt-reverse2-cfft)
    :rfft (transformer :complex :fftr)))

(defn ->complex
  "Convert transformed signal to complex numbers."
  [complex-signal]
  (let [fd (m/seq->double-array complex-signal)]
    (map (fn [^long id]
           (v/vec2 (Array/aget fd id)
                   (Array/aget fd (m/inc id)))) (range 0 (m// (alength fd) 2) 2))))

(defn fft-magnitudes
  [freq-domain]
  (map v/mag (->complex freq-domain)))

(defn fft-phases
  [freq-domain]
  (map v/heading (->complex freq-domain)))

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

;; padding

(defn pad
  "Pad signal."
  ([signal] (pad (m/<< 1 (m/high-2-exp (count signal)))))
  ([signal ^long N] (pad signal N :periodic))
  ([signal ^long N pad-method] (pad signal N pad-method :both))
  ([signal ^long N pad-method side]
   (if (m/> (count signal) N)
     (throw (ex-info "New length of the signal is lower than signal size."
                     {:N N :signal-length (count signal)}))
     (let [asignal (m/seq->double-array signal)]
       (case pad-method
         :zero (pad/zero asignal N side)
         :edge (pad/edge asignal N side)
         :linear (pad/linear asignal N side)
         :periodic (pad/periodic asignal N side)
         :symmetric (pad/symmetric asignal N side)
         :antisymmetric (pad/antisymmetric asignal N side)
         :reflect (pad/reflect asignal N side)
         :antireflect (pad/antireflect asignal N side)
         (throw (ex-info "Unknown padding method" {:pad-method pad-method})))))))

(set! *warn-on-reflection* false)
;; 1d or 2d unknown in the compilation time

(defn compress
  "Compress transformed signal `xs` with given magnitude `mag`."
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


;; https://www.diva-portal.org/smash/get/diva2:1003644/FULLTEXT01.pdf

(defn- sure
  ^double [xs ^long n]
  (let [sxs (double-array (sort (map m/abs xs)))]
    (loop [k (long 1)
           curr (Array/get sxs 0)
           s (* curr curr)
           minrisk ##Inf]
      (if (> k n)
        curr
        (let [v (Array/get sxs (dec k))
              v2 (* v v)
              risk (+ (- n (* 2 k)) s (* (- n k) v2))]
          (if (< risk minrisk)
            (recur (inc k) v (+ s v2) risk)
            (recur (inc k) curr (+ s v2) minrisk)))))))

(defn denoise-threshold
  "Calculate optimal denoise threshold.

  `threshold` is one of the following
  
  * `:visu` - based on median absolute deviation estimate (default)
  * `:universal` - based on standard deviation estimate
  * `:sure` or `:rigrsure` - based on SURE estimator
  * `:hybrid` or `:heursure` - hybrid SURE estimator"
  ^double [xs threshold]
  (let [n (count xs)]
    (if (number? threshold)
      threshold
      (case threshold
        :visu (-> (drop (/ n 2) xs)
                  (stat/median-absolute-deviation)
                  (/ 0.6745)
                  (* (m/sqrt (* 2.0 (m/log n)))))
        :universal (-> (drop (/ n 2) xs)
                       (stat/stddev)
                       (* (m/sqrt (* 2.0 (m/log n)))))
        (:sure :rigrsure) (sure xs n)
        (:hybrid :heursure) (let [eta (/ (- (v/dot xs xs) n) n)
                                  crit (/ (m/pow (m/log2 n) 1.5) (m/sqrt n))]
                              (if (< eta crit)
                                (m/sqrt (* 2.0 (m/log n)))
                                (min (sure xs n) (m/sqrt (* 2.0 (m/log n))))))))))

(defn denoise
  "Wavelet shrinkage with some threshold.

  Methods can be:
  * `:hard` (default)  
  * `:soft`
  * `:garrote`
  * `:hyperbole`

  `:threshold` can be a number of one of the [[denoise-threshold]] methods (default: `:visu`)

  `:skip` can be used to leave `:skip` number of coefficients unaffected (default: 0)

  Use on transformed sequences or call with transformer object."
  ([xs {:keys [method threshold ^long skip]
        :or {method :hard threshold :universal skip 0}}]
   (let [t (double-array xs)
         n (alength t)
         lambda (denoise-threshold xs threshold )
         ids (range skip n)]
     (case method
       :soft (doseq [^long i ids]
               (let [v (Array/aget t i)]
                 (Array/aset t i (* (m/signum v) (max (- (m/abs v) lambda) 0.0)))))
       :hard (doseq [^long i ids]
               (let [v (Array/aget t i)]
                 (when (< (m/abs v) lambda) (Array/aset t i 0.0))))
       :garrote (let [l2 (m/sq lambda)]
                  (doseq [^long i ids]
                    (let [v (Array/aget t i)]
                      (Array/aset t i (if (> (m/abs v) lambda)
                                        (- v (/ l2 v))
                                        0.0)))))
       :hyperbole (let [l2 (m/sq lambda)]
                    (doseq [^long i ids]
                      (let [v (Array/aget t i)]
                        (Array/aset t i (if (> (m/abs v) lambda)
                                          (* (m/signum v) (m/sqrt (- (* v v) l2)))
                                          0.0))))))
     t))
  ([trans xs method]
   (let [v (prot/forward-1d trans xs)]
     (prot/reverse-1d trans (denoise v method))))
  ([xs] (denoise xs nil)))

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
    (with-meta (conj res (if e?
                           (Vec2. (Array/aget txs 1) 0.0)
                           (Vec2. (Array/aget txs len-) (Array/aget txs 1))))
      {::fft {:kind :real :even? e?}})))

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

  Returns a sequence of `fastmath.vector.Vec2` representing complex coefficients in the frequency domain. The result includes metadata (e.g., `:kind`, `:even?`) required by [[ifft]] to correctly perform the inverse transform.

  These are:

  * `:kind` - can be `:real` or `:complex`
  * `:even?` - `true`, when signal length was even
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
    * `:even?` - For `:real` transforms, indicates if the original time-domain signal had an even length (crucial for correctly placing the Nyquist frequency).
    * `:real?` - For `:complex` transforms, indicates if the output should be narrowed to real numbers (doubles).

  Output:
  Returns a sequence representing the time-domain signal. For `:real` kind, it returns a sequence of doubles. For `:complex` kind, it returns a sequence of `Vec2` (complex numbers) unless `:real?` is set to true."
  ([xs] (ifft xs (::fft (meta xs))))
  ([xs {:keys [kind real? even?]
        :or {kind :real real? true even? true}}]
   (if (= :real kind)
     (let [t (transformer :real :fft)
           xs (if even?
                (let [[^double nyquist] (last xs)
                      ^doubles xs (double-array (mapcat identity (butlast xs)))]
                  (Array/aset xs 1 nyquist)
                  xs)
                (let [xs (mapcat identity xs)
                      im (double (last xs))
                      ^doubles xs (double-array (butlast xs))]
                  (Array/aset xs 1 im)
                  xs))]
       (reverse-1d t xs))
     (let [t (transformer :complex (if real? :rfft :fft))
           ^doubles res (reverse-1d t (mapcat identity xs))]
       (if real?
         (take-nth 2 res)
         (reduce (fn [buff ^long id]
                   (conj buff (Vec2. (Array/aget res id)
                                     (Array/aget res (m/inc id))))) [] (range 0 (alength res) 2)))))))

;;





















(comment
  (require '[ggplot])

  
  (defn- error [x1 x2]
    (reduce m/+ (map #(m/sq (- %1 %2)) x1 x2)))

  (def xx (m/slice-range 0 10 512))
  (def s (map #(m/sin %) xx))
  (def d (map #(+ % (* 0.2 (- (rand) 0.5))) s))
  (def t (transformer :fast :daubechies-10))


  (optim/minimize :lbfgsb (fn [x]
                            (error s (denoise t d {:method :garrote :threshold x :skip 0})))
                  {:bounds [[0 (m/sqrt (* 2 (m/log 512)))]]
                   :initial [0.1]})

  (ggplot/->file (ggplot/function
                  (fn [x]
                    (error s (denoise t d {:method :hyperbole :threshold x :skip 0})))
                  {:x [0 (m/sqrt (* 2 (m/log 512)))]}))

  (let [res (denoise t d {:method :hyperbole :threshold :heursure :skip 0})
        data (map (fn [x y] {:x x :y y}) xx res)
        data2 (map (fn [x s y] {:x x :y (- s y)}) xx s res)]
    (ggplot/->file (ggplot/ggaes+ (ggplot/aes :x :x :y :y) (ggplot/line data) (ggplot/line data2 :color "red")))
    (error s res))
  ;; => [0.0054888943941000514 0.0022685134312533306 0.0034123146526252967 0.0027147518947475057]
  ;; => [0.005699029380586008 0.002268315211977693 0.0034120129587574197 0.002714501585972715]


  (seq (denoise (transformer :packet :daubechies-4) [2 3 1 2 3 1 -1 3] :hyperbole))

  (def cfft (transformer :complex :fftr))
  (def res-1 (forward-1d cfft [1 2 3]))
  res-1
  ;; => [6.0, 0.0, -1.5, 0.8660254037844387, -1.5, -0.8660254037844387]
  (seq (reverse-1d cfft res-1))
  ;; => (1.0 0.0 2.0 0.0 3.0 0.0)

  (def rfft (transformer :real :fft))
  (def rres-1 (forward-1d rfft [1 2 3]))
  rres-1
  ;; => [6.0, 0.8660254037844387, -1.5]
  (seq (reverse-1d rfft rres-1))
  ;; => (1.0 2.0 3.0)

  (def packet-symlet-5 (transformer :packet :symlet-5))

  (def res-symlet (forward-1d packet-symlet-5 [1 2 -1 -2]))

  res-symlet
  ;; => [-6.9087374710008476E-12, 2.9275731889992005, 2.82142711971467E-12,
  ;;     -1.1955397203974003]

  (reverse-1d packet-symlet-5 res-symlet)
  ;; => [0.9999999999994849, 1.9999999999989695, -0.9999999999994843,
  ;;     -1.9999999999989686]

  (def res-symlet-2d (forward-2d packet-symlet-5 [[1 2 3 1] [4 0 -1 2] [5 5 9 1] [9 8 -2 -3]]))

  res-symlet-2d
  ;; => [[10.999999999985063, 5.921050975270259, 3.0000000000306537,
  ;;      -1.3932535120463374],
  ;;     [-5.087986539484961, -4.487244758825454, -2.0945087533205813,
  ;;      3.2633503021589165],
  ;;     [2.5000000000411826, -6.782439224284976, 1.500000000025788,
  ;;      -4.061836797399661],
  ;;     [-1.1672159070148709, 4.7633503021476145, -1.3649296986674757,
  ;;      -1.0127552411825524]]

  (reverse-2d packet-symlet-5 res-symlet-2d)
  ;; => [[1.0000000000015503, 2.000000000000519, 3.0000000000010374,
  ;;      1.0000000000002596],
  ;;     [4.000000000000007, 2.5801027980776325E-12, -0.9999999999989707,
  ;;      1.9999999999981963],
  ;;     [5.000000000000007, 4.9999999999982006, 8.999999999997437,
  ;;      1.0000000000010316],
  ;;     [8.999999999995888, 7.9999999999951115, -1.999999999996907,
  ;;      -2.9999999999958775]]

  (def haar (transformer :fast :haar))
  (seq (forward-1d haar [1 2 -1 -3 0 -3 -1 2] {:level 1}))


  ;; => (2.1213203435596424
  ;;     -2.82842712474619
  ;;     -2.1213203435596424
  ;;     0.7071067811865475
  ;;     -0.7071067811865475
  ;;     1.414213562373095
  ;;     2.1213203435596424
  ;;    -2.1213203435596424)
  ;; => (2.1213203435596424
  ;;     -2.82842712474619
  ;;     0.7071067811865475
  ;;     0.7071067811865475
  ;;     -0.7071067811865475
  ;;     1.414213562373095
  ;;     -0.7071067811865475
  ;;    2.1213203435596424)
  ;; => (2.1213203435596424
  ;;     -2.82842712474619
  ;;     1.414213562373095
  ;;     1.414213562373095
  ;;     -0.7071067811865475
  ;;     1.414213562373095
  ;;     -1.414213562373095
  ;;    1.414213562373095)

  (seq (forward-1d haar [1 2 -1 -3 0 2 2 1] {:level 3}))
  ;; => (1.0606601717798212
  ;;     -1.7677669529663682
  ;;     3.499999999999999
  ;;     0.0
  ;;     -0.7071067811865475
  ;;     1.414213562373095
  ;;     -1.414213562373095
  ;;    1.414213562373095)
  ;; => Execution error (JWaveFailure) at jwave.transforms.FastWaveletTransform/forward (FastWaveletTransform.java:82).
  ;;    JWave: Failure: FastWaveletTransform#forward - given level is out of range for given array

  (seq (forward-1d (transformer :dwt "db2") (double-array [1,2,3,10,-3,-2,3,0])))
  ;; => (7.071067811865473
  ;;     3.708909791235272
  ;;     0.006569860407205863
  ;;     -2.274519052838328
  ;;     -4.440892098500626E-16
  ;;     7.563321910700776
  ;;     -4.191872704379808
  ;;    -0.5430220815747799)
  ;; => (7.071067811865473
  ;;     3.708909791235272
  ;;     0.006569860407205863
  ;;     -2.274519052838328
  ;;     -4.440892098500626E-16
  ;;     7.563321910700776
  ;;     -4.191872704379808
  ;;    -0.5430220815747799)
  ;; => (1.5343318992335864
  ;;     12.020815280171306
  ;;     7.14041816598649
  ;;     4.760278777324325
  ;;     -2.8977774788672046
  ;;     5.606086266752902
  ;;     -1.294095225512604
  ;;    -1.414213562373095)


  
  )
;; => nil
;; => nil
