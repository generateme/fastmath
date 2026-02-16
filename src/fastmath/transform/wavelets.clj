(ns fastmath.transform.wavelets
  (:require [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.vector :as v]
            [fastmath.interpolation.linear :as il])
  (:import [fastmath.java Array]
           [java.util Arrays]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defrecord Wavelet [name ^long length ^doubles de-low ^doubles de-high ^doubles re-low ^doubles re-high]
  prot/WaveletProto
  (phi [_] de-low)
  (phi [_ kind] (case kind
                  (:deconstruction :de :dec :d) de-low
                  (:reconstruction :re :rec :r) re-low
                  (throw (ex-info "Unknown phi kind. Can be :rec or :dec." {:kind kind}))))
  (psi [_] de-high)
  (psi [_ kind] (case kind
                  (:deconstruction :de :dec :d) de-high
                  (:reconstruction :re :rec :r) re-high
                  (throw (ex-info "Unknown phi kind. Can be :rec or :dec." {:kind kind}))))

  (coeffs-size [_] length)
  (wavelet-name [_] name)
  (wavelet-forward [_ signal len] (let [len (long len)
                                        target (double-array len)
                                        h (m/>> len 1)]
                                    (dotimes [i h]
                                      (loop [j (long 0)
                                             ai (double 0.0)
                                             aih (double 0.0)]
                                        (if (m/== j length)
                                          (do
                                            (Array/aset target i ai)
                                            (Array/aset target (+ i h) aih))
                                          (let [k (m/mod (m/+ (m/<< i 1) j) len)
                                                v (Array/aget ^doubles signal k)]
                                            (recur (m/inc j)
                                                   (m/+ ai (m/* v (Array/aget de-low j)))
                                                   (m/+ aih (m/* v (Array/aget de-high j))))))))
                                    target))
  
  (wavelet-reverse [_ coeffs len] (let [len (long len)
                                        target (double-array len)
                                        h (m/>> len 1)]                                    
                                    (dotimes [i h]
                                      (let [ai (Array/aget ^doubles coeffs i)
                                            aih (Array/aget ^doubles coeffs (m/+ i h))]
                                        (dotimes [j length]
                                          (let [k (m/mod (m/+ (m/<< i 1) j) len)]
                                            (Array/aset target k (m/+ (Array/aget target k)
                                                                      (m/* ai (Array/aget re-low j))
                                                                      (m/* aih (Array/aget re-high j))))))))
                                    target)))

(defn qmf
  "Quadrature mirror filter used to create scaling coefficients from wavelet coefficients.

  If `inverse?` is true, creates wavelet coefficients from scaling coefficients."
  ([coeffs] (qmf coeffs false))
  ([coeffs inverse?]
   (map-indexed (if inverse?
                  (fn [^long id ^double x] (if (m/odd? id) (m/- x) x))
                  (fn [^long id ^double x] (if (m/even? id) (m/- x) x))) (reverse coeffs))))

(defn ->wavelet
  "Creates wavelet object from coeffients."
  ([nm de-low] (->wavelet nm de-low (qmf de-low true)))
  ([nm de-low de-high] (->wavelet nm de-low de-high de-low de-high))
  ([nm de-low de-high re-low re-high] (->wavelet nm (count de-low) de-low de-high re-low re-high))
  ([nm length de-low de-high re-low re-high]
   (->Wavelet nm length
              (m/seq->double-array de-low)
              (m/seq->double-array de-high)
              (m/seq->double-array re-low)
              (m/seq->double-array re-high))))

(extend-type jwave.transforms.wavelets.Wavelet
  prot/WaveletProto
  (coeffs-size [w] (.getMotherWavelength ^jwave.transforms.wavelets.Wavelet w))
  (wavelet-name [w] (.getName ^jwave.transforms.wavelets.Wavelet w))
  (phi
    ([w] (.getScalingDeComposition ^jwave.transforms.wavelets.Wavelet w))
    ([w kind] (case kind
                (:deconstruction :de :dec :d) (.getScalingDeComposition ^jwave.transforms.wavelets.Wavelet w)
                (:reconstruction :re :rec :r) (.getScalingReConstruction ^jwave.transforms.wavelets.Wavelet w)
                (throw (ex-info "Unknown phi kind. Can be :rec or :dec." {:kind kind})))))
  (psi
    ([w] (.getWaveletDeComposition ^jwave.transforms.wavelets.Wavelet w))
    ([w kind] (case kind
                (:deconstruction :de :dec :d) (.getWaveletDeComposition ^jwave.transforms.wavelets.Wavelet w)
                (:reconstruction :re :rec :r) (.getWaveletReConstruction ^jwave.transforms.wavelets.Wavelet w)
                (throw (ex-info "Unknown phi kind. Can be :rec or :dec." {:kind kind})))))
  (wavelet-forward [w signal len] (.forward ^jwave.transforms.wavelets.Wavelet w signal (int len)))
  (wavelet-reverse [w coeffs len] (.reverse ^jwave.transforms.wavelets.Wavelet w coeffs (int len))))

(defn coeffs-size
  "Returns number of wavelet coefficients."
  ^long [wv] (prot/coeffs-size wv))

(defn wavelet-name
  "Returns wavelet name."
  [wv] (prot/wavelet-name wv))

(defn phi
  "Returns scaling (low-pass) wavelet coefficients, default deconstruction.

  Kind can be:
  
  * `:deconstruction`, `:dec`, `:de` or `:d` - for deconstrunction (default)
  * `:reconstruction`, `:rec`, `:re` or `:r` - for reconstrunction"
  ([wv] (prot/phi wv))
  ([wv kind] (prot/phi wv kind)))

(defn psi
  "Returns wavelet (high-pass) wavelet coefficients, default deconstruction.

  Kind can be:
  
  * `:deconstruction`, `:dec`, `:de` or `:d` - for deconstrunction (default)
  * `:reconstruction`, `:rec`, `:re` or `:r` - for reconstrunction"
  ([wv] (prot/psi wv))
  ([wv kind] (prot/psi wv kind)))

;; matlab and R wavelets coeffs 

(defn- read-wavelets []
  (->> (for [[nm [dl dh rl rh]] (read-string (slurp (io/resource "wavelets/coeffs.edn")))
             :let [s (double-array dl)
                   w (double-array dh)]]
         [nm (if rl
               (->Wavelet nm (alength s) s w (double-array rl) (double-array rh))
               (->Wavelet nm (alength s) s w s w))])
       (into {})))

;; jwave

(defn- jwave-wavelet
  [wavelet]
  (case wavelet
    :haar (jwave.transforms.wavelets.haar.Haar1.)
    :haar-orthogonal (jwave.transforms.wavelets.haar.Haar1Orthogonal.)
    :biorthogonal-11 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal11.)
    :biorthogonal-13 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal13.)
    :biorthogonal-15 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal15.)
    :biorthogonal-22 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal22.)
    :biorthogonal-24 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal24.)
    :biorthogonal-26 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal26.)
    :biorthogonal-28 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal28.)
    :biorthogonal-31 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal31.)
    :biorthogonal-33 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal33.)
    :biorthogonal-35 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal35.)
    :biorthogonal-37 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal37.)
    :biorthogonal-39 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal39.)
    :biorthogonal-44 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal44.)
    :biorthogonal-55 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal55.)
    :biorthogonal-68 (jwave.transforms.wavelets.biorthogonal.BiOrthogonal68.)
    :coiflet-1 (jwave.transforms.wavelets.coiflet.Coiflet1.)
    :coiflet-2 (jwave.transforms.wavelets.coiflet.Coiflet2.)
    :coiflet-3 (jwave.transforms.wavelets.coiflet.Coiflet3.)
    :coiflet-4 (jwave.transforms.wavelets.coiflet.Coiflet4.)
    :coiflet-5 (jwave.transforms.wavelets.coiflet.Coiflet5.)
    :daubechies-2 (jwave.transforms.wavelets.daubechies.Daubechies2.)
    :daubechies-3 (jwave.transforms.wavelets.daubechies.Daubechies3.)
    :daubechies-4 (jwave.transforms.wavelets.daubechies.Daubechies4.)
    :daubechies-5 (jwave.transforms.wavelets.daubechies.Daubechies5.)
    :daubechies-6 (jwave.transforms.wavelets.daubechies.Daubechies6.)
    :daubechies-7 (jwave.transforms.wavelets.daubechies.Daubechies7.)
    :daubechies-8 (jwave.transforms.wavelets.daubechies.Daubechies8.)
    :daubechies-9 (jwave.transforms.wavelets.daubechies.Daubechies9.)
    :daubechies-10 (jwave.transforms.wavelets.daubechies.Daubechies10.)
    :daubechies-11 (jwave.transforms.wavelets.daubechies.Daubechies11.)
    :daubechies-12 (jwave.transforms.wavelets.daubechies.Daubechies12.)
    :daubechies-13 (jwave.transforms.wavelets.daubechies.Daubechies13.)
    :daubechies-14 (jwave.transforms.wavelets.daubechies.Daubechies14.)
    :daubechies-15 (jwave.transforms.wavelets.daubechies.Daubechies15.)
    :daubechies-16 (jwave.transforms.wavelets.daubechies.Daubechies16.)
    :daubechies-17 (jwave.transforms.wavelets.daubechies.Daubechies17.)
    :daubechies-18 (jwave.transforms.wavelets.daubechies.Daubechies18.)
    :daubechies-19 (jwave.transforms.wavelets.daubechies.Daubechies19.)
    :daubechies-20 (jwave.transforms.wavelets.daubechies.Daubechies20.)
    :legendre-1 (jwave.transforms.wavelets.legendre.Legendre1.)
    :legendre-2 (jwave.transforms.wavelets.legendre.Legendre2.)
    :legendre-3 (jwave.transforms.wavelets.legendre.Legendre3.)
    :symlet-2 (jwave.transforms.wavelets.symlets.Symlet2.)
    :symlet-3 (jwave.transforms.wavelets.symlets.Symlet3.)
    :symlet-4 (jwave.transforms.wavelets.symlets.Symlet4.)
    :symlet-5 (jwave.transforms.wavelets.symlets.Symlet5.)
    :symlet-6 (jwave.transforms.wavelets.symlets.Symlet6.)
    :symlet-7 (jwave.transforms.wavelets.symlets.Symlet7.)
    :symlet-8 (jwave.transforms.wavelets.symlets.Symlet8.)
    :symlet-9 (jwave.transforms.wavelets.symlets.Symlet9.)
    :symlet-10 (jwave.transforms.wavelets.symlets.Symlet10.)
    :symlet-11 (jwave.transforms.wavelets.symlets.Symlet11.)
    :symlet-12 (jwave.transforms.wavelets.symlets.Symlet12.)
    :symlet-13 (jwave.transforms.wavelets.symlets.Symlet13.)
    :symlet-14 (jwave.transforms.wavelets.symlets.Symlet14.)
    :symlet-15 (jwave.transforms.wavelets.symlets.Symlet15.)
    :symlet-16 (jwave.transforms.wavelets.symlets.Symlet16.)
    :symlet-17 (jwave.transforms.wavelets.symlets.Symlet17.)
    :symlet-18 (jwave.transforms.wavelets.symlets.Symlet18.)
    :symlet-19 (jwave.transforms.wavelets.symlets.Symlet19.)
    :symlet-20 (jwave.transforms.wavelets.symlets.Symlet20.)
    :battle-23 (jwave.transforms.wavelets.other.Battle23.)
    :cdf-53 (jwave.transforms.wavelets.other.CDF53.)
    :cdf-97 (jwave.transforms.wavelets.other.CDF97.)
    :discrete-mayer (jwave.transforms.wavelets.other.DiscreteMayer.)
    nil))

(def ^:private matlab-wavelets
  (let [wv (read-wavelets)]
    (-> (assoc wv "haar" (wv "db1")))))

(defn wavelet
  "Returns wavelet object containing coefficients"
  [wavelet-name]
  (or (matlab-wavelets wavelet-name) (jwave-wavelet wavelet-name)
      (throw (ex-info "Unknown wavelet." {:wavelet wavelet-name}))))

(defn dwt-forward-1d
  ([wv ^doubles signal]
   (dwt-forward-1d wv signal (m/round (m/log2 (alength signal)))))
  ([wv ^doubles signal ^long level]
   (let [len (alength signal)
         target (Arrays/copyOf signal len)]
     (loop [h len
            l (long 0)]
       (if (or (m/one? h) (m/== l level))
         target
         (let [^doubles step (prot/wavelet-forward wv target h)]
           (System/arraycopy step 0 target 0 h)
           (recur (m/>> h 1) (m/inc l))))))))

(defn dwt-reverse-1d
  ([wv ^doubles signal]
   (dwt-reverse-1d wv signal (m/round (m/log2 (alength signal)))))
  ([wv ^doubles signal ^long level]
   (let [len (alength signal)
         max-level (m/round (m/log2 (alength signal)))
         target (Arrays/copyOf signal len)]
     (loop [h (m/<< 2 (m/max 0 (m/- max-level level)))]
       (if (m/> h len)
         target
         (let [^doubles step (prot/wavelet-reverse wv target h)]
           (System/arraycopy step 0 target 0 h)
           (recur (m/<< h 1))))))))

(defn wpt-forward-1d
  ([wv ^doubles signal]
   (wpt-forward-1d wv signal (m/round (m/log2 (alength signal)))))
  ([wv ^doubles signal ^long level]
   (let [len (alength signal)
         target (Arrays/copyOf signal len)
         tmp (double-array len)]
     (loop [h len
            l (long 0)]
       (if (or (m/one? h) (m/== l level))
         target
         (do (dotimes [p (m// len h)]
               (let [pos (m/* p h)]
                 (System/arraycopy target pos tmp 0 h)
                 (let [^doubles step (prot/wavelet-forward wv tmp h)]
                   (System/arraycopy step 0 target pos h))))
             (recur (m/>> h 1) (m/inc l))))))))

(defn wpt-reverse-1d
  ([wv ^doubles signal]
   (wpt-reverse-1d wv signal (m/round (m/log2 (alength signal)))))
  ([wv ^doubles signal ^long level]
   (let [len (alength signal)
         max-level (m/round (m/log2 (alength signal)))
         target (Arrays/copyOf signal len)
         tmp (double-array len)]
     (loop [h (m/<< 2 (m/max 0 (m/- max-level level)))]
       (if (m/> h len)
         target
         (do (dotimes [p (m// len h)]
               (let [pos (m/* p h)]
                 (System/arraycopy target pos tmp 0 h)
                 (let [^doubles step (prot/wavelet-reverse wv tmp h)]
                   (System/arraycopy step 0 target pos h))))
             (recur (m/<< h 1))))))))

(defn- pow2? [^long v] (and (m/pos? v) (m/zero? (m/bit-and v (m/dec v)))))

(defn- call-transform
  [f wv signal level]
  (let [len (count signal)]
    (if-not (pow2? len)
      (throw (ex-info "Length of the signal should be power of 2." {:length len}))
      (let [s (double-array signal)]
        (if level (f wv s level) (f wv s))))))

(defn wavelet-reify
  [wavelet-name type]
  (let [wv (wavelet wavelet-name)]
    (case type
      :dwt (reify prot/TransformProto
             (forward-1d [o xs] (prot/forward-1d o xs nil))
             (forward-1d [_ xs {:keys [level]}] (call-transform dwt-forward-1d wv xs level))
             (reverse-1d [o xs] (prot/reverse-1d o xs nil))
             (reverse-1d [_ xs {:keys [level]}] (call-transform dwt-reverse-1d wv xs level)))
      :wpt (reify prot/TransformProto
             (forward-1d [o xs] (prot/forward-1d o xs nil))
             (forward-1d [_ xs {:keys [level]}] (call-transform wpt-forward-1d wv xs level))
             (reverse-1d [o xs] (prot/reverse-1d o xs nil))
             (reverse-1d [_ xs {:keys [level]}] (call-transform wpt-reverse-1d wv xs level))))))

;;

;; https://github.com/vincentherrmann/vincentherrmann.github.io/blob/master/scripts/wavelets.js

(defn- upscale
  [coeffs buff]
  (let [buff-cnt (count buff)
        coeffs-cnt (count coeffs)
        nsize (m/+ (m/* 2 buff-cnt) coeffs-cnt -2)]
    (mapv (fn [^long j]
            (reduce (fn [^double s ^long k]
                      (if (m/even? (m/+ k j))
                        (let [i (m// (m/- j k) 2)]
                          (if (and (m/not-neg? i) (m/< i buff-cnt))
                            (m/+ s (m/* (double (coeffs k)) (double (buff i))))
                            s))
                        s)) 0.0 (range coeffs-cnt))) (range nsize))))

(defn- upscale-final
  [coeffs buff scale]
  (let [coeffs-cnt (count coeffs)
        offset (m/round (m/exp2 scale))]
    (->> coeffs
         (map-indexed (fn [^long i ^double coeff]
                        (concat (repeat (m/* i offset) 0.0)
                                (v/mult buff coeff)
                                (repeat (m/* (m/- coeffs-cnt i 1) offset) 0.0))))
         (reduce v/add))))

(defn- function-coeffs
  [w f reconstruction?]
  (-> (if reconstruction? (f w :rec) (f w))
      (reverse)
      (v/mult m/SQRT2)
      (vec)))

(defn- function-step1
  [wv ^long level reconstruction?]
  (let [w (if (or (string? wv) (keyword? wv)) (wavelet wv) wv)
        low (function-coeffs w phi reconstruction?)]
    [w (nth (iterate (partial upscale low) [1.0]) (m/dec level)) low]))

(defn- function-interpolator
  [w y]
  (let [s (m/dec (coeffs-size w))
        x (m/slice-range 0 s (count y))
        i (il/linear x y)]
    (fn [^double x] (if (m/<= 0 x s) (i x) 0.0))))

(defn scaling-function
  "Returns approximation of scaling function, phi."
  ([wv] (scaling-function wv 8))
  ([wv ^long level] (scaling-function wv level false))
  ([wv ^long level reconstruction?]
   (let [[w step-1 low] (function-step1 wv level reconstruction?)]
     (function-interpolator w (upscale-final low step-1 (m/dec level))))))

(defn wavelet-function
  "Returns approximation of wavelet function, psi."
  ([wv] (wavelet-function wv 8))
  ([wv ^long level] (wavelet-function wv level false))
  ([wv ^long level reconstruction?]
   (let [[w step-1] (function-step1 wv level reconstruction?)
         ;; negate?
         high (function-coeffs w psi reconstruction?)]
     (function-interpolator w (upscale-final high step-1 (m/dec level))))))
