(ns fastmath.transform
  "Discrete transforms of 1d (and some 2d) signals.

  High-level functions work on sequences and return vectors with metadata needed by their inverses:

  * Fourier: [[fft]], [[ifft]] - real input (single- or double-sided spectrum) and complex input (sequence of `Vec2`, see `fastmath.complex`)
  * Cosine: [[dct]], [[idct]] - types I, II and III
  * Sine: [[dst]], [[idst]] - types I, II and III
  * Hartley: [[dht]], [[idht]]
  * Walsh-Hadamard: [[hadamard]], [[ihadamard]]
  * Wavelets: [[dwt]], [[idwt]], [[wpt]], [[iwpt]], [[wpd]] (see also `fastmath.transform.wavelets`)

  Fourier, cosine, sine and Hartley functions accept a `:norm` option with the meaning used by `scipy.fft`:

  * `:backward` - forward transform not scaled, inverse scaled by `1/N` (default for [[fft]] and [[dht]])
  * `:ortho` - both directions scaled to an orthonormal transform (default for [[dct]] and [[dst]])
  * `:forward` - forward transform scaled by `1/N`, inverse not scaled

  For cosine and sine transforms `N` is the scipy normalization length (`2N`, `2(N-1)` or `2(N+1)` depending on the type), so the results equal `scipy.fft.dct` / `dst` with the same `type` and `norm`.

  Low-level transform objects are created by [[transformer]] and used with [[forward-1d]], [[reverse-1d]], [[forward-2d]] and [[reverse-2d]]. They wrap JTransforms (FFT, DCT-II/III, DST-II/III, DHT, also 2d), Apache Commons Math (DCT-I/DST-I of power-of-two sizes, Hadamard) and the wavelet transforms.

  Fourier transforms are computed by JTransforms."
  (:require [fastmath.core :as m]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.transform.wavelets :as wv])
  (:import [org.apache.commons.math3.transform FastSineTransformer FastCosineTransformer FastHadamardTransformer RealTransformer DstNormalization DctNormalization TransformType]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]
           [org.jtransforms.dht DoubleDHT_1D DoubleDHT_2D]
           [org.jtransforms.dct DoubleDCT_1D DoubleDCT_2D]
           [org.jtransforms.dst DoubleDST_1D DoubleDST_2D]
           [fastmath.java Array]
           [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

;; guards

(defn- unknown
  [what value allowed]
  (throw (ex-info (str "Unknown " what ".") {what value :allowed allowed})))

(defn- non-empty!
  "Throws when `xs` is empty, returns `xs`."
  [xs]
  (when (or (nil? xs) (zero? (count xs)))
    (throw (ex-info "Signal is empty." {:length 0})))
  xs)

(defn- even-length!
  "Throws when interleaved complex data has odd length, returns `xs`."
  [xs]
  (when (odd? (count xs))
    (throw (ex-info "Interleaved complex data needs an even number of values." {:length (count xs)})))
  xs)

(defn- non-empty-2d!
  [xss]
  (when (or (nil? xss) (zero? (count xss)) (zero? (count (first xss))))
    (throw (ex-info "2d signal is empty." {:rows (count xss)})))
  xss)

(defmulti
  ^{:doc "Creates a low-level transform object, used with [[forward-1d]], [[reverse-1d]], [[forward-2d]] and [[reverse-2d]].

  Prefer the high-level functions ([[fft]], [[dct]], [[dst]], [[dht]], [[hadamard]], [[dwt]], [[wpt]]) unless you need the raw layout or a 2d transform.

  First argument selects the family, second the variant:

  * `:real`
      * `:fft` - real FFT in the JTransforms packed layout (1d and 2d)
      * `:dct` - DCT-II forward, DCT-III reverse (1d and 2d)
      * `:dst` - DST-II forward, DST-III reverse (1d and 2d)
      * `:dht` - discrete Hartley transform (1d and 2d)
      * `:sine`, `:cosine` - DST-I / DCT-I from Apache Commons Math, options `{:normalization :orthogonal}` (default) or `:standard`; size must be a power of two (`:sine`, first value must be 0) or a power of two plus one (`:cosine`)
      * `:hadamard` - Walsh-Hadamard transform, size must be a power of two
  * `:complex`
      * `:fft` - complex FFT of interleaved `[re0 im0 re1 im1 ...]` data (1d and 2d)
      * `:fftr` (alias `:rfft`) - forward: full complex spectrum of real data; reverse: complex inverse of interleaved data
  * `:sine`, `:cosine` - with second argument `:standard` or `:orthogonal`, Apache Commons Math DST-I / DCT-I with fixed normalization
  * `:dwt`, `:wpt` - wavelet transforms ([[dwt]], [[wpt]]) with flat `double-array` coefficients, second argument is a wavelet name, see `fastmath.transform.wavelets/wavelet-names`; options `:level` and `:mode` as in [[dwt]], the reverse transform also needs the signal length `:n` (and the same `:level` and `:mode`). 1d only.

  JTransforms based objects take the option `:scale?` (default `true`): when `true` the reverse FFT/DHT is divided by `N` and DCT/DST are orthonormal; with `false` the raw JTransforms scaling is used (for DCT/DST it differs between power-of-two and other sizes). The high-level functions use `:norm` instead.

  2d transforms need more than one row and column, the real 2d FFT (`:real :fft`) also power-of-two sizes; JTransforms throws `IllegalArgumentException` otherwise.

  Throws `ex-info` for an unknown family or variant."}
  transformer (fn [t _] t))

(defmethod transformer :default [t _]
  (unknown :transformer t [:real :complex :sine :cosine :dwt :wpt]))

;; Apache Commons Math sine/cosine (DST-I / DCT-I)

(defn- no-2d
  [t & _]
  (throw (ex-info "2d transform is not supported by this transformer." {:transformer (type t)})))

(defn- acm-dst-normalization
  [normalization]
  (case normalization
    :standard DstNormalization/STANDARD_DST_I
    :orthogonal DstNormalization/ORTHOGONAL_DST_I
    (unknown :normalization normalization [:standard :orthogonal])))

(defn- acm-dct-normalization
  [normalization]
  (case normalization
    :standard DctNormalization/STANDARD_DCT_I
    :orthogonal DctNormalization/ORTHOGONAL_DCT_I
    (unknown :normalization normalization [:standard :orthogonal])))

(defn- acm-transformer
  ^RealTransformer [kind normalization]
  (case kind
    :sine (FastSineTransformer. (acm-dst-normalization normalization))
    :cosine (FastCosineTransformer. (acm-dct-normalization normalization))))

(defn- acm-sine-cosine-transformer
  "Transform object for ACM DST-I / DCT-I with normalization given in options."
  [kind]
  (letfn [(run [xs normalization forward?]
            (.transform (acm-transformer kind (or normalization :orthogonal))
                        (m/seq->double-array (non-empty! xs))
                        (if forward? TransformType/FORWARD TransformType/INVERSE)))]
    (reify prot/TransformProto
      (forward-1d [_ xs] (run xs nil true))
      (forward-1d [_ xs {:keys [normalization]}] (run xs normalization true))
      (reverse-1d [_ xs] (run xs nil false))
      (reverse-1d [_ xs {:keys [normalization]}] (run xs normalization false))
      (forward-2d [t _] (no-2d t))
      (forward-2d [t _ _] (no-2d t))
      (reverse-2d [t _] (no-2d t))
      (reverse-2d [t _ _] (no-2d t)))))

(defmethod transformer :sine [_ normalization] (acm-transformer :sine normalization))
(defmethod transformer :cosine [_ normalization] (acm-transformer :cosine normalization))

;; Hadamard and fixed-normalization ACM transforms

(extend RealTransformer
  prot/TransformProto
  {:forward-1d (fn
                 ([^RealTransformer t xs] (.transform t (m/seq->double-array (non-empty! xs)) TransformType/FORWARD))
                 ([^RealTransformer t xs _] (.transform t (m/seq->double-array (non-empty! xs)) TransformType/FORWARD)))
   :reverse-1d (fn
                 ([^RealTransformer t xs] (.transform t (m/seq->double-array (non-empty! xs)) TransformType/INVERSE))
                 ([^RealTransformer t xs _] (.transform t (m/seq->double-array (non-empty! xs)) TransformType/INVERSE)))
   :forward-2d no-2d
   :reverse-2d no-2d})

;; JTransforms

(defn- jt-forward-fft [xs _]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.realForward (DoubleFFT_1D. (alength out)) out)
    out))

(defn- jt-reverse-fft [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.realInverse (DoubleFFT_1D. (alength out)) out (boolean scale?))
    out))

(defn- ->double-double-array
  ^"[[D" [xss]
  (into-array (map m/seq->double-array (non-empty-2d! xss))))

(defn- jt-forward2-fft [xss _]
  (let [out (->double-double-array xss)]
    (.realForward (DoubleFFT_2D. (alength out) (alength ^doubles (aget out 0))) out)
    out))

(defn- jt-reverse2-fft [xss scale?]
  (let [out (->double-double-array xss)]
    (.realInverse (DoubleFFT_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
    out))

(defn- jt-forward-cfft [xs _]
  (let [out (m/seq->double-array (even-length! (non-empty! xs)))]
    (.complexForward (DoubleFFT_1D. (m/long-quot (alength out) 2)) out)
    out))

(defn- jt-reverse-cfft [xs scale?]
  (let [out (m/seq->double-array (even-length! (non-empty! xs)))]
    (.complexInverse (DoubleFFT_1D. (m/long-quot (alength out) 2)) out (boolean scale?))
    out))

(defn- ->interleaved-rows
  "Rows of interleaved complex data, truncated to an even length."
  ^"[[D" [xss]
  (let [s2 (m/long-mult 2 (m/long-quot (count (first (non-empty-2d! xss))) 2))]
    (into-array (map (fn [xs] (m/seq->double-array (take s2 xs))) xss))))

(defn- jt-forward2-cfft [xss _]
  (let [out (->interleaved-rows xss)]
    (.complexForward (DoubleFFT_2D. (alength out) (m/long-quot (alength ^doubles (aget out 0)) 2)) out)
    out))

(defn- jt-reverse2-cfft [xss scale?]
  (let [out (->interleaved-rows xss)]
    (.complexInverse (DoubleFFT_2D. (alength out) (m/long-quot (alength ^doubles (aget out 0)) 2)) out (boolean scale?))
    out))

(defn- jt-forward-cfftr [xs _]
  (let [in (m/seq->double-array (non-empty! xs))
        s (alength in)
        out (double-array (m/long-mult 2 s))]
    (System/arraycopy in 0 out 0 s)
    (.realForwardFull (DoubleFFT_1D. s) out)
    out))

(defn- jt-forward2-cfftr [xss _]
  (let [rows (count (non-empty-2d! xss))
        s (count (first xss))
        ^"[[D" out (into-array (map (fn [xs]
                                      (let [o (double-array (m/long-mult 2 s))]
                                        (System/arraycopy (m/seq->double-array xs) 0 o 0 s)
                                        o)) xss))]
    (.realForwardFull (DoubleFFT_2D. rows s) out)
    out))

(defn- jt-forward-dht [xs _]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.forward (DoubleDHT_1D. (alength out)) out)
    out))

(defn- jt-reverse-dht [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.inverse (DoubleDHT_1D. (alength out)) out (boolean scale?))
    out))

(defn- jt-forward2-dht [xss _]
  (let [out (->double-double-array xss)]
    (.forward (DoubleDHT_2D. (alength out) (alength ^doubles (aget out 0))) out)
    out))

(defn- jt-reverse2-dht [xss scale?]
  (let [out (->double-double-array xss)]
    (.inverse (DoubleDHT_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
    out))

(defn- jt-forward-dct [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.forward (DoubleDCT_1D. (alength out)) out (boolean scale?))
    out))

(defn- jt-reverse-dct [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.inverse (DoubleDCT_1D. (alength out)) out (boolean scale?))
    out))

(defn- jt-forward2-dct [xss scale?]
  (let [out (->double-double-array xss)]
    (.forward (DoubleDCT_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
    out))

(defn- jt-reverse2-dct [xss scale?]
  (let [out (->double-double-array xss)]
    (.inverse (DoubleDCT_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
    out))

(defn- jt-forward-dst [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.forward (DoubleDST_1D. (alength out)) out (boolean scale?))
    out))

(defn- jt-reverse-dst [xs scale?]
  (let [out (m/seq->double-array (non-empty! xs))]
    (.inverse (DoubleDST_1D. (alength out)) out (boolean scale?))
    out))

(defn- jt-forward2-dst [xss scale?]
  (let [out (->double-double-array xss)]
    (.forward (DoubleDST_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
    out))

(defn- jt-reverse2-dst [xss scale?]
  (let [out (->double-double-array xss)]
    (.inverse (DoubleDST_2D. (alength out) (alength ^doubles (aget out 0))) out (boolean scale?))
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
    :sine (acm-sine-cosine-transformer :sine)
    :cosine (acm-sine-cosine-transformer :cosine)
    :hadamard (FastHadamardTransformer.)
    (unknown :variant t [:fft :dht :dct :dst :sine :cosine :hadamard])))

(defmethod transformer :complex [_ t]
  (case t
    :fft (jt-reify jt-forward-cfft jt-reverse-cfft jt-forward2-cfft jt-reverse2-cfft)
    (:fftr :rfft) (jt-reify jt-forward-cfftr jt-reverse-cfft jt-forward2-cfftr jt-reverse2-cfft)
    (unknown :variant t [:fft :fftr :rfft])))

;;

(defn forward-1d
  "Forward transform of a sequence or array with a transform object created by [[transformer]].

  `options` are passed to the transform object, see [[transformer]]. Returns a `double-array` (raw layout of the underlying library)."
  ([t xs] (prot/forward-1d t xs))
  ([t xs options] (prot/forward-1d t xs options)))

(defn reverse-1d
  "Reverse (inverse) transform of a sequence or array with a transform object created by [[transformer]].

  `options` are passed to the transform object, see [[transformer]]. Returns a `double-array`."
  ([t xs] (prot/reverse-1d t xs))
  ([t xs options] (prot/reverse-1d t xs options)))

(defn forward-2d
  "Forward 2d transform of a sequence of rows with a transform object created by [[transformer]].

  `options` are passed to the transform object, see [[transformer]]. Returns an array of `double-array` rows (raw layout of the underlying library)."
  ([t xss] (prot/forward-2d t xss))
  ([t xss options] (prot/forward-2d t xss options)))

(defn reverse-2d
  "Reverse (inverse) 2d transform of a sequence of rows with a transform object created by [[transformer]].

  `options` are passed to the transform object, see [[transformer]]. Returns an array of `double-array` rows."
  ([t xss] (prot/reverse-2d t xss))
  ([t xss options] (prot/reverse-2d t xss options)))

;; normalization

(def ^:private norms [:backward :ortho :forward])

(defn- norm-scale
  "Scale applied by a Fourier-type transform of length `n` (normalization length) for `norm` and direction."
  ^double [norm ^long n inverse?]
  (case norm
    :backward (if inverse? (/ 1.0 n) 1.0)
    :ortho (/ 1.0 (m/sqrt n))
    :forward (if inverse? 1.0 (/ 1.0 n))
    (unknown :norm norm norms)))

(defn- scale!
  "Multiplies `arr` in place by `s`, returns `arr`."
  ^doubles [^doubles arr ^double s]
  (when-not (== s 1.0)
    (dotimes [i (alength arr)]
      (Array/aset arr i (* s (Array/aget arr i)))))
  arr)

;; Fourier

(defn- complex-seq?
  [xs]
  (instance? Vec2 (first xs)))

(defn- ->complex-vec
  "Interleaved `[re im ...]` array to a vector of `Vec2`."
  [^doubles arr]
  (mapv (fn [^long id] (Vec2. (Array/aget arr id) (Array/aget arr (inc id)))) (range 0 (alength arr) 2)))

(defn- interleave-complex
  ^doubles [xs]
  (let [zs (vec xs)
        n (count zs)
        arr (double-array (* 2 n))]
    (dotimes [id n]
      (let [^Vec2 z (zs id)]
        (Array/aset arr (* 2 id) (.x z))
        (Array/aset arr (inc (* 2 id)) (.y z))))
    arr))

(defn- real-fft
  "Single-sided spectrum of real `xs`: `N/2+1` coefficients."
  [xs norm]
  (let [^doubles packed (jt-forward-fft xs nil)
        n (alength packed)
        s (norm-scale norm n false)
        last-full (if (even? n) (dec (m/long-quot n 2)) (m/long-quot n 2))
        res (transient [(Vec2. (* s (Array/aget packed 0)) 0.0)])]
    (loop [k 1 res res]
      (if (> k last-full)
        (persistent! (if (even? n)
                       (conj! res (Vec2. (* s (Array/aget packed 1)) 0.0))
                       res))
        (recur (inc k)
               (conj! res (if (and (odd? n) (== k last-full))
                            (Vec2. (* s (Array/aget packed (dec n))) (* s (Array/aget packed 1)))
                            (Vec2. (* s (Array/aget packed (* 2 k))) (* s (Array/aget packed (inc (* 2 k))))))))))))

(defn fft
  "Computes the discrete Fourier transform of a 1d signal.

  Input `xs` is a sequence of numbers (real signal) or a sequence of complex numbers (`Vec2`, see `fastmath.complex`).

  Options:

  * `:spectrum` - for a real signal: `:single-sided` (default) returns the `N/2+1` non-negative frequency coefficients (as `numpy.fft.rfft`), `:double-sided` returns all `N` coefficients (as `numpy.fft.fft`). Ignored for complex input.
  * `:norm` - `:backward` (default, not scaled), `:ortho` (scaled by `1/sqrt(N)`), `:forward` (scaled by `1/N`).

  Returns a vector of `Vec2`. Metadata under `:fastmath.transform/fft` holds `{:kind :n :norm :real?}` used by [[ifft]]: `:kind` is `:real` for a single-sided spectrum and `:complex` otherwise, `:n` is the signal length, `:real?` is `true` for the double-sided spectrum of a real signal.

  Throws `ex-info` for an empty signal or an unknown option value."
  ([xs] (fft xs nil))
  ([xs {:keys [spectrum norm]
        :or {spectrum :single-sided norm :backward}}]
   (non-empty! xs)
   (norm-scale norm 1 false)
   (cond
     (complex-seq? xs) (let [arr (jt-forward-cfft (interleave-complex xs) nil)
                             n (count xs)]
                         (with-meta (->complex-vec (scale! arr (norm-scale norm n false)))
                           {::fft {:kind :complex :n n :norm norm :real? false}}))
     (number? (first xs)) (let [n (count xs)]
                            (case spectrum
                              :single-sided (with-meta (real-fft xs norm)
                                              {::fft {:kind :real :n n :norm norm}})
                              :double-sided (with-meta (->complex-vec (scale! (jt-forward-cfftr xs nil)
                                                                              (norm-scale norm n false)))
                                              {::fft {:kind :complex :n n :norm norm :real? true}})
                              (unknown :spectrum spectrum [:single-sided :double-sided])))
     :else (throw (ex-info "Signal should contain numbers or complex numbers (Vec2)." {:first (first xs)})))))

(defn- real-ifft
  [xs ^long n norm]
  (let [m (inc (m/long-quot n 2))]
    (when-not (== m (count xs))
      (throw (ex-info "Single-sided spectrum length should be N/2+1." {:n n :length (count xs) :expected m})))
    (let [packed (double-array n)
          zs (vec xs)
          ^Vec2 z0 (zs 0)]
      (Array/aset packed 0 (.x z0))
      (doseq [k (range 1 m)
              :let [k (long k)
                    ^Vec2 z (zs k)]]
        (cond
          (and (even? n) (== k (m/long-quot n 2))) (Array/aset packed 1 (.x z))
          (and (odd? n) (== k (m/long-quot n 2))) (do (Array/aset packed (dec n) (.x z))
                                               (when (> n 1) (Array/aset packed 1 (.y z))))
          :else (do (Array/aset packed (* 2 k) (.x z))
                    (Array/aset packed (inc (* 2 k)) (.y z)))))
      ;; JTransforms' unscaled real inverse is N/2 times the signal for even N, the scaled one is exact
      (scale! (jt-reverse-fft packed true) (* n (norm-scale norm n true))))))

(defn ifft
  "Computes the inverse discrete Fourier transform, the inverse of [[fft]].

  `xs` is a sequence of complex coefficients (`Vec2`). The way it is inverted is read from the metadata attached by [[fft]]; `options` override the metadata:

  * `:kind` - `:real` for a single-sided spectrum (`N/2+1` coefficients) of a real signal, `:complex` for a full spectrum.
  * `:n` - length of the original signal, required for `:real` (`N` cannot be inferred from a single-sided spectrum). For `:complex` it is the number of coefficients.
  * `:real?` - for `:complex`, return only real parts (default `false`).
  * `:norm` - `:backward` (default, scaled by `1/N`), `:ortho` (scaled by `1/sqrt(N)`), `:forward` (not scaled).

  Returns a vector of doubles for `:real` (or `:real? true`), otherwise a vector of `Vec2`.

  Throws `ex-info` when `xs` is empty, when there is no metadata and `:kind` is not given, when `:n` is missing for `:kind :real` or does not match the spectrum length."
  ([xs] (ifft xs nil))
  ([xs options]
   (non-empty! xs)
   (let [{:keys [kind n real? norm] :or {norm :backward}} (merge (::fft (meta xs)) options)]
     (case kind
       :real (do (when-not n
                   (throw (ex-info "Signal length :n is required to invert a single-sided spectrum." {:options options})))
                 (vec (real-ifft xs (long n) norm)))
       :complex (let [cnt (count xs)
                      res (scale! (jt-reverse-cfft (interleave-complex xs) true) (* cnt (norm-scale norm cnt true)))]
                  (if real?
                    (vec (take-nth 2 res))
                    (->complex-vec res)))
       nil (throw (ex-info "Unknown spectrum kind: provide metadata from fft or the :kind option." {:options options}))
       (unknown :kind kind [:real :complex])))))

;; DCT / DST

(defn- type-23-weights
  "Per-element multipliers which turn the orthonormal DCT/DST-II/III into the scipy `norm` variant.

  For type II: forward `y = ortho(x) / w`, inverse `x = ortho-inverse(y * w)`.
  For type III: forward `y = ortho-inverse(x * w)`, inverse `x = ortho(y) / w`.
  `special` is the index with a different weight (0 for cosine, `n-1` for sine)."
  ^doubles [transform-type norm ^long n ^long special]
  (let [w (double-array n)]
    (dotimes [i n]
      (let [base (if (== 2 (long transform-type))
                   (if (== i special) (m/sqrt (/ 0.25 n)) (m/sqrt (/ 0.5 n)))
                   (if (== i special) (m/sqrt n) (m/sqrt (* 2.0 n))))]
        (Array/aset w i (double (case norm
                                  :ortho 1.0
                                  :backward base
                                  :forward (if (== 2 (long transform-type)) (* base 2.0 n) (/ base (* 2.0 n)))
                                  (unknown :norm norm norms))))))
    w))

(defn- emul!
  ^doubles [^doubles a ^doubles w]
  (dotimes [i (alength a)] (Array/aset a i (* (Array/aget a i) (Array/aget w i)))) a)

(defn- ediv!
  ^doubles [^doubles a ^doubles w]
  (dotimes [i (alength a)] (Array/aset a i (/ (Array/aget a i) (Array/aget w i)))) a)

(defn- type-23
  "Type II/III transform (`transform-type` 2 or 3, `inverse?`) on top of the orthonormal JTransforms one."
  [ortho-forward ortho-inverse special-fn transform-type norm inverse? xs]
  (let [arr (m/seq->double-array (non-empty! xs))
        n (alength arr)
        w (type-23-weights transform-type norm n (special-fn n))]
    (if (== 2 (long transform-type))
      (if inverse?
        (ortho-inverse (emul! arr w) true)
        (ediv! (ortho-forward arr true) w))
      (if inverse?
        (ediv! (ortho-forward arr true) w)
        (ortho-inverse (emul! arr w) true)))))

(defn- dct-1
  "Unnormalized DCT-I (scipy `backward`) via FFT of the even extension, length `n >= 2`."
  ^doubles [^doubles x]
  (let [n (alength x)
        len (* 2 (dec n))
        ext (double-array len)]
    (System/arraycopy x 0 ext 0 n)
    (loop [i 1]
      (when (< i (dec n))
        (Array/aset ext (- len i) (Array/aget x i))
        (recur (inc i))))
    (let [^doubles spec (jt-forward-cfftr ext nil)
          out (double-array n)]
      (dotimes [k n] (Array/aset out k (Array/aget spec (* 2 k))))
      out)))

(defn- dst-1
  "Unnormalized DST-I (scipy `backward`) via FFT of the odd extension, length `n >= 1`."
  ^doubles [^doubles x]
  (let [n (alength x)
        len (* 2 (inc n))
        ext (double-array len)]
    (dotimes [i n]
      (Array/aset ext (inc i) (Array/aget x i))
      (Array/aset ext (- len i 1) (- (Array/aget x i))))
    (let [^doubles spec (jt-forward-cfftr ext nil)
          out (double-array n)]
      (dotimes [k n] (Array/aset out k (- (Array/aget spec (inc (* 2 (inc k)))))))
      out)))

(defn- type-1
  "DCT-I or DST-I with scipy `norm`; the inverse equals the forward transform with `:backward` and `:forward` swapped."
  [kind norm inverse? xs]
  (let [norm (if inverse? (case norm :backward :forward :forward :backward norm) norm)
        ^doubles x (m/seq->double-array (non-empty! xs))
        n (alength x)]
    (case kind
      :cosine (do
                (when (< n 2) (throw (ex-info "DCT-I needs at least 2 values." {:length n})))
                (let [len (* 2 (dec n))]
                  (case norm
                    :backward (dct-1 x)
                    :forward (scale! (dct-1 x) (/ 1.0 len))
                    :ortho (let [r2 (m/sqrt 2.0)
                                 xs (doto ^doubles (aclone x)
                                      (Array/aset 0 (* r2 (Array/aget x 0)))
                                      (Array/aset (dec n) (* r2 (Array/aget x (dec n)))))
                                 ^doubles y (scale! (dct-1 xs) (/ 1.0 (m/sqrt len)))]
                             (Array/aset y 0 (/ (Array/aget y 0) r2))
                             (Array/aset y (dec n) (/ (Array/aget y (dec n)) r2))
                             y)
                    (unknown :norm norm norms))))
      :sine (let [len (* 2 (inc n))]
              (case norm
                :backward (dst-1 x)
                :forward (scale! (dst-1 x) (/ 1.0 len))
                :ortho (scale! (dst-1 x) (/ 1.0 (m/sqrt len)))
                (unknown :norm norm norms))))))

(def ^:private dct-methods [:DCT-I :DCT-II :DCT-III])
(def ^:private dst-methods [:DST-I :DST-II :DST-III])

(defn- cosine-sine
  [kind method norm inverse? xs]
  (case kind
    :cosine (case method
              :DCT-I (type-1 :cosine norm inverse? xs)
              :DCT-II (type-23 jt-forward-dct jt-reverse-dct (constantly 0) 2 norm inverse? xs)
              :DCT-III (type-23 jt-forward-dct jt-reverse-dct (constantly 0) 3 norm inverse? xs)
              (unknown :method method dct-methods))
    :sine (case method
            :DST-I (type-1 :sine norm inverse? xs)
            :DST-II (type-23 jt-forward-dst jt-reverse-dst dec 2 norm inverse? xs)
            :DST-III (type-23 jt-forward-dst jt-reverse-dst dec 3 norm inverse? xs)
            (unknown :method method dst-methods))))

(defn dct
  "Computes the discrete cosine transform (DCT) of a 1d signal.

  Options:

  * `:method` - `:DCT-I`, `:DCT-II` (default) or `:DCT-III`
  * `:norm` - `:ortho` (default, orthonormal), `:backward` (not scaled) or `:forward`, as in `scipy.fft.dct`

  The result equals `scipy.fft.dct(xs, type, norm)`. `:DCT-I` needs at least 2 values, other types at least 1; any length is accepted.

  Returns a vector of doubles. Metadata under `:fastmath.transform/dct` holds `{:method :norm}` used by [[idct]].

  Throws `ex-info` for an empty signal or an unknown option value."
  ([xs] (dct xs nil))
  ([xs {:keys [method norm] :or {method :DCT-II norm :ortho}}]
   (with-meta (vec (cosine-sine :cosine method norm false xs))
     {::dct {:method method :norm norm}})))

(defn idct
  "Computes the inverse discrete cosine transform, the inverse of [[dct]].

  `:method` and `:norm` are read from the metadata attached by [[dct]]; `options` override them. Without metadata the defaults are `:DCT-II` and `:ortho`. The result equals `scipy.fft.idct(xs, type, norm)`.

  Returns a vector of doubles. Throws `ex-info` for an empty input or an unknown option value."
  ([xs] (idct xs nil))
  ([xs options]
   (let [{:keys [method norm] :or {method :DCT-II norm :ortho}} (merge (::dct (meta xs)) options)]
     (vec (cosine-sine :cosine method norm true xs)))))

(defn dst
  "Computes the discrete sine transform (DST) of a 1d signal.

  Options:

  * `:method` - `:DST-I`, `:DST-II` (default) or `:DST-III`
  * `:norm` - `:ortho` (default, orthonormal), `:backward` (not scaled) or `:forward`, as in `scipy.fft.dst`

  The result equals `scipy.fft.dst(xs, type, norm)`. Any length of at least 1 is accepted.

  Returns a vector of doubles. Metadata under `:fastmath.transform/dst` holds `{:method :norm}` used by [[idst]].

  Throws `ex-info` for an empty signal or an unknown option value."
  ([xs] (dst xs nil))
  ([xs {:keys [method norm] :or {method :DST-II norm :ortho}}]
   (with-meta (vec (cosine-sine :sine method norm false xs))
     {::dst {:method method :norm norm}})))

(defn idst
  "Computes the inverse discrete sine transform, the inverse of [[dst]].

  `:method` and `:norm` are read from the metadata attached by [[dst]]; `options` override them. Without metadata the defaults are `:DST-II` and `:ortho`. The result equals `scipy.fft.idst(xs, type, norm)`.

  Returns a vector of doubles. Throws `ex-info` for an empty input or an unknown option value."
  ([xs] (idst xs nil))
  ([xs options]
   (let [{:keys [method norm] :or {method :DST-II norm :ortho}} (merge (::dst (meta xs)) options)]
     (vec (cosine-sine :sine method norm true xs)))))

;; Hartley

(defn dht
  "Computes the discrete Hartley transform `H[k] = sum x[n] (cos(2 pi k n / N) + sin(2 pi k n / N))` of a 1d signal.

  Option `:norm`: `:backward` (default, not scaled), `:ortho` (scaled by `1/sqrt(N)`), `:forward` (scaled by `1/N`).

  Returns a vector of doubles with metadata `{:norm}` under `:fastmath.transform/dht` used by [[idht]]. Throws `ex-info` for an empty signal or an unknown `:norm`."
  ([xs] (dht xs nil))
  ([xs {:keys [norm] :or {norm :backward}}]
   (let [n (count (non-empty! xs))]
     (with-meta (vec (scale! (jt-forward-dht xs nil) (norm-scale norm n false)))
       {::dht {:norm norm}}))))

(defn idht
  "Computes the inverse discrete Hartley transform, the inverse of [[dht]].

  `:norm` is read from the metadata attached by [[dht]]; `options` override it (default `:backward`, the result is divided by `N`).

  Returns a vector of doubles. Throws `ex-info` for an empty input or an unknown `:norm`."
  ([xs] (idht xs nil))
  ([xs options]
   (let [{:keys [norm] :or {norm :backward}} (merge (::dht (meta xs)) options)
         n (count (non-empty! xs))]
     (vec (scale! (jt-reverse-dht xs true) (* n (norm-scale norm n true)))))))

;; Hadamard

(defn- power-of-two!
  [xs]
  (let [n (count (non-empty! xs))]
    (when-not (m/power-of-two? n)
      (throw (ex-info "Length of the signal should be a power of 2." {:length n})))
    xs))

(defn hadamard
  "Computes the fast Walsh-Hadamard transform (natural / Sylvester ordering, not normalized) of a 1d signal.

  The length must be a power of 2. Returns a vector of doubles. Throws `ex-info` for an empty signal or another length."
  [xs]
  (vec (forward-1d (FastHadamardTransformer.) (power-of-two! xs))))

(defn ihadamard
  "Computes the inverse fast Walsh-Hadamard transform (the forward transform divided by `N`).

  The length must be a power of 2. Returns a vector of doubles. Throws `ex-info` for an empty signal or another length."
  [xs]
  (vec (reverse-1d (FastHadamardTransformer.) (power-of-two! xs))))

(m/unuse-primitive-operators)

;; Wavelets

(defn- wavelet-options
  "Validated `[wavelet mode]` and the filter length."
  [wavelet mode]
  (let [w (wv/wavelet wavelet)]
    (when-not (wv/modes mode)
      (throw (ex-info "Unknown wavelet mode." {:mode mode :allowed (sort wv/modes)})))
    [w (wv/coeffs-size w)]))

(defn- resolve-level
  ^long [level ^long n ^long filter-length]
  (cond
    (nil? level) (wv/max-level n filter-length)
    (and (integer? level) (m/not-neg? (long level))) (long level)
    :else (throw (ex-info "Level should be a non-negative integer." {:level level}))))

(defn- signal->array
  ^doubles [xs]
  (let [arr (m/seq->double-array xs)]
    (when (m/zero? (alength arr)) (throw (ex-info "Signal is empty." {:length 0})))
    arr))

(defn- dwt-bands
  "Bands `[cA_level cD_level ... cD_1]` as `double-array`s."
  [w ^doubles x mode ^long level]
  (loop [a x l 0 details ()]
    (if (m/== l level)
      (vec (cons a details))
      (let [[ca cd] (wv/dwt-level w a mode)]
        (recur ca (m/long-inc l) (cons cd details))))))

(defn- band-lengths
  "Lengths of the [[dwt]] bands for a signal of length `n`."
  [^long n ^long filter-length mode ^long level]
  (let [cs (rest (take (m/long-inc level) (iterate #(wv/coeffs-length % filter-length mode) n)))]
    (if (empty? cs)
      [n]
      (vec (cons (last cs) (reverse cs))))))

(defn- split-flat
  [coeffs lengths]
  (when-not (== (count coeffs) (long (reduce + lengths)))
    (throw (ex-info "Flat coefficients do not match band lengths." {:length (count coeffs) :lengths lengths})))
  (loop [cs (vec coeffs) ls lengths acc []]
    (if-let [l (first ls)]
      (recur (subvec cs l) (rest ls) (conj acc (subvec cs 0 l)))
      acc)))

(defn- trim-to
  "First `n` values of `arr` when `n` is given, otherwise `arr`."
  ^doubles [^doubles arr n]
  (cond
    (nil? n) arr
    (== (alength arr) (long n)) arr
    (> (alength arr) (long n)) (java.util.Arrays/copyOf arr (int n))
    :else (throw (ex-info "Reconstructed signal is shorter than the requested length." {:length (alength arr) :n n}))))

(defn dwt
  "Computes the multilevel discrete wavelet transform of a 1d signal, as `pywt.wavedec`.

  * `xs` - signal, any length
  * `wavelet` - name from `fastmath.transform.wavelets/wavelet-names` (e.g. `\"db4\"`, `\"sym8\"`, `:daubechies-4`) or a wavelet object
  * `options`:
      * `:level` - number of decomposition levels, default: the maximum useful level `floor(log2(N / (L - 1)))` (`pywt.dwt_max_level`)
      * `:mode` - boundary extension, default `:symmetric`; see `fastmath.transform.wavelets` for all modes, `:periodization` keeps the total length
      * `:flat?` - return one concatenated vector instead of bands (default `false`)

  Returns a vector of bands `[cA_n cD_n cD_n-1 ... cD_1]` (approximation at the coarsest level followed by details from coarsest to finest), each a vector of doubles; with `:level 0` the only band is the signal. Metadata under `:fastmath.transform/dwt` holds `{:wavelet :mode :level :n}` (and `:lengths` of the bands for flat output), used by [[idwt]].

  Throws `ex-info` for an empty signal, an unknown wavelet or mode, a negative level, or a one-sample signal with `:reflect` / `:antireflect`."
  ([xs wavelet] (dwt xs wavelet nil))
  ([xs wavelet {:keys [level mode flat?] :or {mode :symmetric}}]
   (let [[w l] (wavelet-options wavelet mode)
         x (signal->array xs)
         n (alength x)
         level (resolve-level level n l)
         bands (mapv vec (dwt-bands w x mode level))
         info {:wavelet wavelet :mode mode :level level :n n}]
     (if flat?
       (with-meta (vec (apply concat bands)) {::dwt (assoc info :lengths (mapv count bands))})
       (with-meta bands {::dwt info})))))

(defn idwt
  "Reconstructs a signal from the [[dwt]] bands, as `pywt.waverec`.

  * `coeffs` - vector of bands `[cA_n cD_n ... cD_1]`, or a flat vector together with `:lengths`
  * `wavelet` - the wavelet used by [[dwt]]
  * `options` (merged over the metadata attached by [[dwt]]):
      * `:mode` - boundary mode used by [[dwt]] (default `:symmetric`)
      * `:n` - length of the original signal; without it the result can be one sample longer, as in PyWavelets
      * `:lengths` - band lengths, required for flat coefficients

  At every level an approximation one sample longer than the next detail band is trimmed (as in PyWavelets). Returns a vector of doubles.

  Throws `ex-info` for empty input, an unknown wavelet or mode, flat coefficients without `:lengths`, or band lengths that do not fit together."
  ([coeffs wavelet] (idwt coeffs wavelet nil))
  ([coeffs wavelet options]
   (let [{:keys [mode n lengths] :or {mode :symmetric}} (merge (::dwt (meta coeffs)) options)
         [w _] (wavelet-options wavelet mode)
         _ (when (empty? coeffs) (throw (ex-info "Coefficients are empty." {})))
         bands (if (number? (first coeffs))
                 (if lengths
                   (split-flat coeffs lengths)
                   (throw (ex-info "Flat coefficients need :lengths of the bands." {})))
                 coeffs)
         res (reduce (fn [^doubles a d]
                       (let [^doubles d (m/seq->double-array d)
                             la (alength a)
                             ld (alength d)
                             ^doubles a (cond
                                          (== la ld) a
                                          (== la (inc ld)) (java.util.Arrays/copyOf a ld)
                                          :else (throw (ex-info "Band lengths do not fit together." {:approximation la :detail ld})))]
                         (wv/idwt-level w a d mode)))
                     (m/seq->double-array (first bands)) (rest bands))]
     (vec (trim-to res n)))))

(defn- packet-nodes
  "Nodes of the wavelet packet tree at each level `0..level`, natural order."
  [w ^doubles x mode ^long level]
  (take (m/long-inc level)
        (iterate (fn [nodes] (vec (mapcat (fn [node] (wv/dwt-level w node mode)) nodes))) [x])))

(defn wpt
  "Computes the wavelet packet transform of a 1d signal at one level, as the level nodes of `pywt.WaveletPacket`.

  Every node is split into approximation and detail by one level of the discrete wavelet transform, down to `:level`.

  * `xs` - signal, any length
  * `wavelet` - name from `fastmath.transform.wavelets/wavelet-names` or a wavelet object
  * `options`:
      * `:level` - decomposition level, default the maximum useful level (as for [[dwt]])
      * `:mode` - boundary extension, default `:symmetric`
      * `:flat?` - return one concatenated vector instead of nodes (default `false`)

  Returns the `2^level` nodes of the level in natural order (paths `a...a`, `a...ad`, ..., `d...d`, approximation before detail at every split), each a vector of doubles of the same length. Metadata under `:fastmath.transform/wpt` holds `{:wavelet :mode :level :n}`, used by [[iwpt]].

  Throws `ex-info` for an empty signal, an unknown wavelet or mode, or a negative level."
  ([xs wavelet] (wpt xs wavelet nil))
  ([xs wavelet {:keys [level mode flat?] :or {mode :symmetric}}]
   (let [[w l] (wavelet-options wavelet mode)
         x (signal->array xs)
         n (alength x)
         level (resolve-level level n l)
         nodes (mapv vec (last (packet-nodes w x mode level)))]
     (with-meta (if flat? (vec (apply concat nodes)) nodes)
       {::wpt {:wavelet wavelet :mode mode :level level :n n}}))))

(defn iwpt
  "Reconstructs a signal from the [[wpt]] nodes of one level.

  * `coeffs` - vector of `2^level` nodes in natural order, or a flat vector together with `:level`
  * `wavelet` - the wavelet used by [[wpt]]
  * `options` (merged over the metadata attached by [[wpt]]):
      * `:mode` - boundary mode used by [[wpt]] (default `:symmetric`)
      * `:level` - required for flat coefficients
      * `:n` - length of the original signal; every reconstructed node is trimmed to its length in the tree (as PyWavelets does with a known data shape), without it nodes can be one sample longer

  Returns a vector of doubles. Throws `ex-info` for empty input, an unknown wavelet or mode, a node count which is not a power of two, nodes of different lengths or lengths which do not match `:n`."
  ([coeffs wavelet] (iwpt coeffs wavelet nil))
  ([coeffs wavelet options]
   (let [{:keys [mode n level] :or {mode :symmetric}} (merge (::wpt (meta coeffs)) options)
         [w l] (wavelet-options wavelet mode)
         _ (when (empty? coeffs) (throw (ex-info "Coefficients are empty." {})))
         nodes (if (number? (first coeffs))
                 (if level
                   (let [parts (long (m/exp2 level))]
                     (when-not (zero? (m/long-mod (count coeffs) parts))
                       (throw (ex-info "Flat coefficients can't be split into 2^level equal nodes." {:length (count coeffs) :level level})))
                     (mapv vec (partition (quot (count coeffs) parts) coeffs)))
                   (throw (ex-info "Flat coefficients need :level." {})))
                 coeffs)
         cnt (count nodes)
         lvl (m/log2int cnt)]
     (when-not (m/power-of-two? cnt)
       (throw (ex-info "Number of nodes should be a power of two." {:nodes cnt})))
     (when-not (apply == (map count nodes))
       (throw (ex-info "Nodes should have the same length." {:lengths (map count nodes)})))
     (let [tree-lengths (when n (vec (take (inc lvl) (iterate #(wv/coeffs-length % l mode) (long n)))))]
       (when (and tree-lengths (not= (count (first nodes)) (tree-lengths lvl)))
         (throw (ex-info "Node length does not match the signal length :n." {:node-length (count (first nodes)) :expected (tree-lengths lvl)})))
       (loop [nodes (mapv m/seq->double-array nodes) k lvl]
         (if (zero? k)
           (vec ^doubles (first nodes))
           (recur (mapv (fn [[a d]] (trim-to (wv/idwt-level w a d mode) (when tree-lengths (tree-lengths (dec k)))))
                        (partition 2 nodes))
                  (dec k))))))))

(defn wpd
  "Computes the full wavelet packet decomposition: the [[wpt]] nodes of every level from 0 (the signal) to `:level`.

  * `xs` - signal, any length
  * `wavelet` - name from `fastmath.transform.wavelets/wavelet-names` or a wavelet object
  * `options`:
      * `:level` - deepest level, default the maximum useful level (as for [[dwt]])
      * `:mode` - boundary extension, default `:symmetric`
      * `:flat?` - every level as one concatenated vector instead of nodes (default `false`)

  Returns a vector of levels; level `k` holds `2^k` nodes in natural order, as [[wpt]] with `:level k`. Metadata under `:fastmath.transform/wpd` holds `{:wavelet :mode :level :n}`.

  Throws `ex-info` for an empty signal, an unknown wavelet or mode, or a negative level."
  ([xs wavelet] (wpd xs wavelet nil))
  ([xs wavelet {:keys [level mode flat?] :or {mode :symmetric}}]
   (let [[w l] (wavelet-options wavelet mode)
         x (signal->array xs)
         n (alength x)
         level (resolve-level level n l)]
     (with-meta (mapv (fn [nodes] (let [nodes (mapv vec nodes)]
                                    (if flat? (vec (apply concat nodes)) nodes)))
                      (packet-nodes w x mode level))
       {::wpd {:wavelet wavelet :mode mode :level level :n n}}))))

(defn- wavelet-transformer
  "Transform object over [[dwt]] / [[wpt]] with flat `double-array` coefficients; the reverse transform needs `:n`."
  [kind wavelet]
  (let [l (wv/coeffs-size (wv/wavelet wavelet))]
    (letfn [(fwd [xs opts]
              (double-array ((if (= kind :dwt) dwt wpt) xs wavelet (assoc opts :flat? true))))
            (rev [xs {:keys [n mode level] :or {mode :symmetric} :as opts}]
              (when-not n
                (throw (ex-info "Reverse wavelet transform needs the signal length :n." {:options opts})))
              (let [level (resolve-level level n l)]
                (double-array (if (= kind :dwt)
                                (idwt (vec xs) wavelet (assoc opts :lengths (band-lengths n l mode level) :level level))
                                (iwpt (vec xs) wavelet (assoc opts :level level))))))]
      (reify prot/TransformProto
        (forward-1d [_ xs] (fwd xs nil))
        (forward-1d [_ xs opts] (fwd xs opts))
        (reverse-1d [_ xs] (rev xs nil))
        (reverse-1d [_ xs opts] (rev xs opts))
        (forward-2d [t _] (no-2d t))
        (forward-2d [t _ _] (no-2d t))
        (reverse-2d [t _] (no-2d t))
        (reverse-2d [t _ _] (no-2d t))))))

(defmethod transformer :dwt [_ w] (wavelet-transformer :dwt w))
(defmethod transformer :wpt [_ w] (wavelet-transformer :wpt w))
