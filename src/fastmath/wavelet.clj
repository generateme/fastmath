(ns fastmath.wavelet
  "Wavelet transformations.

  Based on [JWave](https://github.com/cscheiblich/JWave/) library.

  Be aware that some of the wavelet types doesn't work properly. `:battle-23`, `:cdf-53`, `:cdf-97`."
  {:metadoc/categories {:w "Wavelet"}}
  (:require [fastmath.core :as m]
            [metadoc.examples :refer :all])
  (:import [jwave.transforms FastWaveletTransform WaveletPacketTransform AncientEgyptianDecomposition BasicTransform]
           [jwave.exceptions JWaveFailure]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;;

(defmulti
  ^{:doc "Create wavelet object.

  Shouldn't be used directly"
    :metadoc/categories #{:w}
    :metadoc/examples [(example "Usage" (wavelet :haar))]} wavelet identity)

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

(def ^{:metadoc/categories #{:w}
       :doc "List of all possible wavelets."
       :metadoc/examples [(example "List of wavelets" (sort wavelets-list))]}
  wavelets-list (remove #{:default} (keys (methods wavelet))))

(defmulti
  ^{:doc "Create transform object for given wavelet.

  You can use:

  * `:fast` for 1d or 2d Fast Wavelet Transform. Size of data should be power of `2`.
  * `:packet` for 1d or 2d Wavelet Packet Transform. Size of data should be power of `2`.
  * `:decomposed-fast` for 1d Fast Wavelet Transform. Data can have any size (Ancient Egyptian Decomposition is used).
  * `:decomposed-packet` for 1d Wavelet Packet Transform. Data can have any size (Ancient Egyptian Decomposition is used).

  Second argument is wavelet name as key. See [[wavelets-list]] for all supported names."
    :metadoc/categories #{:w}
    :metadoc/examples [(example "Usage" (transformer :packet :discrete-mayer))]}
  transformer (fn [t w] t))

(defmethod transformer :fast [_ w] (FastWaveletTransform. (wavelet w)))
(defmethod transformer :packet [_ w] (WaveletPacketTransform. (wavelet w)))
(defmethod transformer :decomposed-fast [_ w] (AncientEgyptianDecomposition. (transformer :fast w)))
(defmethod transformer :decomposed-packet [_ w] (AncientEgyptianDecomposition. (transformer :packet w)))

(defprotocol TransformProto
  "Transformer functions."
  (^{:metadoc/categories #{:w}} forward-1d [t xs] "Forward transform of sequence or array. Returns double array.")
  (^{:metadoc/categories #{:w}} reverse-1d [t xs] "Reverse transform of sequence or array. Returns double array.")
  (^{:metadoc/categories #{:w}} forward-2d [t xss] "Forward transform of sequence of sequences or 2d double array. Returns 2d double array.")
  (^{:metadoc/categories #{:w}} reverse-2d [t xss] "Reverse transform of sequence of sequences or 2d double array. Returns 2d double array."))

(extend BasicTransform
  TransformProto
  {:forward-1d (fn [^BasicTransform t xs] (.forward t (m/seq->double-array xs)))
   :reverse-1d (fn [^BasicTransform t xs] (.reverse t (m/seq->double-array xs)))
   :forward-2d (fn [^BasicTransform t xss] (.forward t (m/seq->double-double-array xss)))
   :reverse-2d (fn [^BasicTransform t xss] (.reverse t (m/seq->double-double-array xss)))})

(add-examples forward-1d
  (example-session "Usage"
    (seq (forward-1d (transformer :packet :haar-orthogonal) [-1 8 7 6]))
    (seq (forward-1d (transformer :fast :haar-orthogonal) [-1 8 7 6]))
    (seq (forward-1d (transformer :decomposed-fast :haar-orthogonal) [10 -1 8 7 6]))))

(add-examples reverse-1d
  (example-session "Usage"
    (let [t (transformer :packet :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))
    (let [t (transformer :fast :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))
    (let [t (transformer :decomposed-fast :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [10 -1 8 7 6]))))))

(add-examples forward-2d
  (example-session "Usage"
    (m/double-double-array->seq (forward-2d (transformer :packet :daubechies-6) [[-1 8] [7 6]]))
    (m/double-double-array->seq (forward-2d (transformer :fast :daubechies-6) [[-1 8] [7 6]]))))

(add-examples reverse-2d
  (example-session "Usage"
    (let [t (transformer :packet :daubechies-20)]
      (m/double-double-array->seq (reverse-2d t (forward-2d t [[-1 8] [7 6]]))))
    (let [t (transformer :fast :daubechies-20)]
      (m/double-double-array->seq (reverse-2d t (forward-2d t [[-1 8] [7 6]]))))))

(set! *warn-on-reflection* false)

(defn compress
  "Compress transformed signal `xs` with given magnitude `mag`."
  {:metadoc/categories #{:w}
   :metadoc/examples [(example "Compress 1d"
                        (let [t (transformer :fast :symlet-5)]
                          (->> [1 2 3 4]
                               (forward-1d t)
                               (compress 0.3)
                               (reverse-1d t)
                               (seq))))
                      (example "Compress 2d"
                        (let [t (transformer :fast :symlet-5)]
                          (->> [[1 2] [3 4]]
                               (forward-2d t)
                               (compress 0.5)
                               (reverse-2d t)
                               (m/double-double-array->seq))))]}
  [^double mag xs]
  (.compress (jwave.compressions.CompressorMagnitude. mag) xs))
