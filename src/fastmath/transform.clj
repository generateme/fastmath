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
            [fastmath.protocols :as prot]
            [fastmath.optimization :as optim]
            [fastmath.vector :as v])
  (:import [jwave.transforms FastWaveletTransform WaveletPacketTransform AncientEgyptianDecomposition
            BasicTransform DiscreteFourierTransform]
           [jwave.exceptions JWaveFailure]
           [jwave.compressions CompressorPeaksAverage CompressorMagnitude]
           [org.apache.commons.math3.transform FastSineTransformer FastCosineTransformer FastHadamardTransformer RealTransformer DstNormalization DctNormalization TransformType]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]
           [org.jtransforms.dht DoubleDHT_1D DoubleDHT_2D]
           [org.jtransforms.dct DoubleDCT_1D DoubleDCT_2D]
           [org.jtransforms.dst DoubleDST_1D DoubleDST_2D]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

;;

(defmulti
  ^{:doc "Create wavelet object."
    :private true} wavelet identity)

(defmethod wavelet :default [n] (throw (JWaveFailure. (str "Unknown wavelet: " n))))

(defmethod wavelet :haar [_] (jwave.transforms.wavelets.haar.Haar1.))
(defmethod wavelet :haar-orthogonal [_] (jwave.transforms.wavelets.haar.Haar1Orthogonal.))

(defmethod wavelet :biorthogonal-11 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal11.))
(defmethod wavelet :biorthogonal-13 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal13.))
(defmethod wavelet :biorthogonal-15 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal15.))
(defmethod wavelet :biorthogonal-22 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal22.))
(defmethod wavelet :biorthogonal-24 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal24.))
(defmethod wavelet :biorthogonal-26 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal26.))
(defmethod wavelet :biorthogonal-28 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal28.))
(defmethod wavelet :biorthogonal-31 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal31.))
(defmethod wavelet :biorthogonal-33 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal33.))
(defmethod wavelet :biorthogonal-35 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal35.))
(defmethod wavelet :biorthogonal-37 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal37.))
(defmethod wavelet :biorthogonal-39 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal39.))
(defmethod wavelet :biorthogonal-44 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal44.))
(defmethod wavelet :biorthogonal-55 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal55.))
(defmethod wavelet :biorthogonal-68 [_] (jwave.transforms.wavelets.biorthogonal.BiOrthogonal68.))

(defmethod wavelet :coiflet-1 [_] (jwave.transforms.wavelets.coiflet.Coiflet1.))
(defmethod wavelet :coiflet-2 [_] (jwave.transforms.wavelets.coiflet.Coiflet2.))
(defmethod wavelet :coiflet-3 [_] (jwave.transforms.wavelets.coiflet.Coiflet3.))
(defmethod wavelet :coiflet-4 [_] (jwave.transforms.wavelets.coiflet.Coiflet4.))
(defmethod wavelet :coiflet-5 [_] (jwave.transforms.wavelets.coiflet.Coiflet5.))

(defmethod wavelet :daubechies-2 [_] (jwave.transforms.wavelets.daubechies.Daubechies2.))
(defmethod wavelet :daubechies-3 [_] (jwave.transforms.wavelets.daubechies.Daubechies3.))
(defmethod wavelet :daubechies-4 [_] (jwave.transforms.wavelets.daubechies.Daubechies4.))
(defmethod wavelet :daubechies-5 [_] (jwave.transforms.wavelets.daubechies.Daubechies5.))
(defmethod wavelet :daubechies-6 [_] (jwave.transforms.wavelets.daubechies.Daubechies6.))
(defmethod wavelet :daubechies-7 [_] (jwave.transforms.wavelets.daubechies.Daubechies7.))
(defmethod wavelet :daubechies-8 [_] (jwave.transforms.wavelets.daubechies.Daubechies8.))
(defmethod wavelet :daubechies-9 [_] (jwave.transforms.wavelets.daubechies.Daubechies9.))
(defmethod wavelet :daubechies-10 [_] (jwave.transforms.wavelets.daubechies.Daubechies10.))
(defmethod wavelet :daubechies-11 [_] (jwave.transforms.wavelets.daubechies.Daubechies11.))
(defmethod wavelet :daubechies-12 [_] (jwave.transforms.wavelets.daubechies.Daubechies12.))
(defmethod wavelet :daubechies-13 [_] (jwave.transforms.wavelets.daubechies.Daubechies13.))
(defmethod wavelet :daubechies-14 [_] (jwave.transforms.wavelets.daubechies.Daubechies14.))
(defmethod wavelet :daubechies-15 [_] (jwave.transforms.wavelets.daubechies.Daubechies15.))
(defmethod wavelet :daubechies-16 [_] (jwave.transforms.wavelets.daubechies.Daubechies16.))
(defmethod wavelet :daubechies-17 [_] (jwave.transforms.wavelets.daubechies.Daubechies17.))
(defmethod wavelet :daubechies-18 [_] (jwave.transforms.wavelets.daubechies.Daubechies18.))
(defmethod wavelet :daubechies-19 [_] (jwave.transforms.wavelets.daubechies.Daubechies19.))
(defmethod wavelet :daubechies-20 [_] (jwave.transforms.wavelets.daubechies.Daubechies20.))

(defmethod wavelet :legendre-1 [_] (jwave.transforms.wavelets.legendre.Legendre1.))
(defmethod wavelet :legendre-2 [_] (jwave.transforms.wavelets.legendre.Legendre2.))
(defmethod wavelet :legendre-3 [_] (jwave.transforms.wavelets.legendre.Legendre3.))

(defmethod wavelet :symlet-2 [_] (jwave.transforms.wavelets.symlets.Symlet2.))
(defmethod wavelet :symlet-3 [_] (jwave.transforms.wavelets.symlets.Symlet3.))
(defmethod wavelet :symlet-4 [_] (jwave.transforms.wavelets.symlets.Symlet4.))
(defmethod wavelet :symlet-5 [_] (jwave.transforms.wavelets.symlets.Symlet5.))
(defmethod wavelet :symlet-6 [_] (jwave.transforms.wavelets.symlets.Symlet6.))
(defmethod wavelet :symlet-7 [_] (jwave.transforms.wavelets.symlets.Symlet7.))
(defmethod wavelet :symlet-8 [_] (jwave.transforms.wavelets.symlets.Symlet8.))
(defmethod wavelet :symlet-9 [_] (jwave.transforms.wavelets.symlets.Symlet9.))
(defmethod wavelet :symlet-10 [_] (jwave.transforms.wavelets.symlets.Symlet10.))
(defmethod wavelet :symlet-11 [_] (jwave.transforms.wavelets.symlets.Symlet11.))
(defmethod wavelet :symlet-12 [_] (jwave.transforms.wavelets.symlets.Symlet12.))
(defmethod wavelet :symlet-13 [_] (jwave.transforms.wavelets.symlets.Symlet13.))
(defmethod wavelet :symlet-14 [_] (jwave.transforms.wavelets.symlets.Symlet14.))
(defmethod wavelet :symlet-15 [_] (jwave.transforms.wavelets.symlets.Symlet15.))
(defmethod wavelet :symlet-16 [_] (jwave.transforms.wavelets.symlets.Symlet16.))
(defmethod wavelet :symlet-17 [_] (jwave.transforms.wavelets.symlets.Symlet17.))
(defmethod wavelet :symlet-18 [_] (jwave.transforms.wavelets.symlets.Symlet18.))
(defmethod wavelet :symlet-19 [_] (jwave.transforms.wavelets.symlets.Symlet19.))
(defmethod wavelet :symlet-20 [_] (jwave.transforms.wavelets.symlets.Symlet20.))

(defmethod wavelet :battle-23 [_] (jwave.transforms.wavelets.other.Battle23.))
(defmethod wavelet :cdf-53 [_] (jwave.transforms.wavelets.other.CDF53.))
(defmethod wavelet :cdf-97 [_] (jwave.transforms.wavelets.other.CDF97.))
(defmethod wavelet :discrete-mayer [_] (jwave.transforms.wavelets.other.DiscreteMayer.))

(def ^{:doc "List of all possible wavelets."}
  wavelets-list (remove #{:default} (keys (methods wavelet))))

(defmulti
  ^{:doc "Create transform object for given wavelet.

  #### Wavelets

  * `:fast` for 1d or 2d Fast Wavelet Transform. Size of data should be power of `2`.
  * `:packet` for 1d or 2d Wavelet Packet Transform. Size of data should be power of `2`.
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

(defmethod transformer :fast [_ w] (FastWaveletTransform. (wavelet w)))
(defmethod transformer :packet [_ w] (WaveletPacketTransform. (wavelet w)))
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
  {:forward-1d (fn [^BasicTransform t xs] (.forward t (m/seq->double-array xs)))
   :reverse-1d (fn [^BasicTransform t xs] (.reverse t (m/seq->double-array xs)))
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
    :fftr (jt-reify jt-forward-cfftr jt-reverse-cfft jt-forward2-cfftr jt-reverse2-cfft)))

;;

(defn forward-1d
  "Forward transform of sequence or array."
  [t xs] (prot/forward-1d t xs))

(defn reverse-1d
  "Forward transform of sequence or array."
  [t xs] (prot/reverse-1d t xs))

(defn forward-2d
  "Forward transform of sequence or array."
  [t xss] (prot/forward-2d t xss))

(defn reverse-2d
  "Forward transform of sequence or array."
  [t xss] (prot/reverse-2d t xss))

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
     (println "lambda:" lambda)
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


(comment
  (require '[ggplot])


  (defn- error [x1 x2]
    (reduce m/fast+ (map #(m/sq (- %1 %2)) x1 x2)))

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


  (seq (denoise (transformer :packet :daubechies-4) [2 3 1 2 3 1 -1 3] :hyperbole)))

