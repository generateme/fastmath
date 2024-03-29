(ns fastmath.transform
  "Transforms.

  See [[transformer]] and [[TransformProto]] for details.
  
  ### Wavelet
  
  Based on [JWave](https://github.com/cscheiblich/JWave/) library.

  Be aware that some of the wavelet types doesn't work properly. `:battle-23`, `:cdf-53`, `:cdf-97`.

  ### Cos/Sin/Hadamard

  Orthogonal or standard fast sine/cosine/hadamard 1d transforms.

  ### Fourier

  DFT."
  (:require [fastmath.core :as m]
            [fastmath.stats :as stat]
            [fastmath.protocols :as prot])
  (:import [jwave.transforms FastWaveletTransform WaveletPacketTransform AncientEgyptianDecomposition
            BasicTransform DiscreteFourierTransform]
           [jwave.exceptions JWaveFailure]
           [org.apache.commons.math3.transform FastSineTransformer FastCosineTransformer FastHadamardTransformer RealTransformer
            DstNormalization DctNormalization TransformType]))

(set! *unchecked-math* :warn-on-boxed)
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

(defmethod transformer :standard [_ t] (cond
                                         (= t :sine) (FastSineTransformer. DstNormalization/STANDARD_DST_I)
                                         (= t :cosine) (FastCosineTransformer. DctNormalization/STANDARD_DCT_I)
                                         (= t :hadamard) (FastHadamardTransformer.)
                                         (= t :dft) (DiscreteFourierTransform.)))

(defmethod transformer :orthogonal [_ t] (cond
                                           (= t :sine) (FastSineTransformer. DstNormalization/ORTHOGONAL_DST_I)
                                           (= t :cosine) (FastCosineTransformer. DctNormalization/ORTHOGONAL_DCT_I)))

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
          (.compress (jwave.compressions.CompressorMagnitude. mag))
          (rev trans))))
  ([xs ^double mag]
   (.compress (jwave.compressions.CompressorMagnitude. mag) xs)))

(defn denoise
  "Adaptive denoising of time series (1d).

  Use on transformed sequences or call with transformer object.

  SMILE implementation of WaveletShrinkage denoise function."
  ([xs soft?]
   (let [t (double-array xs)
         n (alength t)
         nh (m/>> n 1)
         wc (double-array nh)]
     (System/arraycopy t nh wc 0 nh)
     (let [error (/ (stat/median-absolute-deviation wc) 0.6745)
           lambda (* error (m/sqrt (* 2.0 (m/log n))))]
       (if soft?
         (dotimes [i (- n 2)]
           (let [i2 (+ i 2)
                 v (aget t i2)]
             (aset-double t i2 (* (m/signum v) (max (- (m/abs v) lambda) 0.0)))))
         (dotimes [i (- n 2)]
           (let [i2 (+ i 2)
                 v (aget t i2)]
             (when (< (m/abs v) lambda) (aset-double t i2 0.0)))))
       t)))
  ([trans xs soft?]
   (let [v (prot/forward-1d trans xs)]
     (prot/reverse-1d trans (denoise v soft?))))
  ([xs] (denoise xs false)))
