(ns fastmath.transform.wavelets
  "Wavelet filter banks and the discrete wavelet transform engine used by `fastmath.transform`.

  A wavelet is a filter bank `[dec_lo dec_hi rec_lo rec_hi]` in PyWavelets orientation, see [[wavelet]] and [[wavelet-names]]. The transforms follow PyWavelets: [[dwt-level]] / [[idwt-level]] is one level of `pywt.dwt` / `pywt.idwt` (convolution with the decomposition filters and downsampling), with the signal extended at both ends according to a boundary `mode`:

  * `:symmetric` (default in `fastmath.transform`), `:reflect`, `:periodic`, `:zero`, `:constant`, `:smooth`, `:antisymmetric`, `:antireflect` - output length `floor((N + L - 1) / 2)` for a signal of length `N` and filters of length `L`; the extension is done by `fastmath.signal.pad`
  * `:periodization` - output length `ceil(N / 2)`, an odd-length signal is extended by its last value

  Multi-level transforms, wavelet packets and their inverses are in `fastmath.transform` ([[fastmath.transform/dwt]], [[fastmath.transform/wpt]]).

  [[scaling-function]] and [[wavelet-function]] approximate the functions phi and psi by the cascade algorithm."
  (:require [clojure.java.io :as io]
            [clojure.edn :as edn]
            [fastmath.core :as m]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.signal.pad :as pad]
            [fastmath.interpolation.linear :as il])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defrecord Wavelet [name ^long length ^doubles dec-lo ^doubles dec-hi ^doubles rec-lo ^doubles rec-hi]
  prot/WaveletProto
  (phi [_] dec-lo)
  (phi [_ kind] (case kind
                  (:deconstruction :de :dec :d) dec-lo
                  (:reconstruction :re :rec :r) rec-lo
                  (throw (ex-info "Unknown phi kind. Can be :rec or :dec." {:kind kind}))))
  (psi [_] dec-hi)
  (psi [_ kind] (case kind
                  (:deconstruction :de :dec :d) dec-hi
                  (:reconstruction :re :rec :r) rec-hi
                  (throw (ex-info "Unknown psi kind. Can be :rec or :dec." {:kind kind}))))
  (coeffs-size [_] length)
  (wavelet-name [_] name))

(defn qmf
  "Quadrature mirror filter of `coeffs`, as `pywt.qmf`: the filter reversed, with every odd element negated.

  For an orthogonal wavelet the reconstruction high-pass filter is `(qmf rec-lo)`."
  [coeffs]
  (vec (map-indexed (fn [^long id ^double x] (if (m/odd? id) (m/- x) x)) (reverse coeffs))))

(defn ->wavelet
  "Creates a wavelet (filter bank) object, PyWavelets orientation.

  * `[name dec-lo]` - orthogonal wavelet from its decomposition low-pass filter: `rec-lo` is `dec-lo` reversed, `rec-hi` is `(qmf rec-lo)`, `dec-hi` is `rec-hi` reversed
  * `[name dec-lo dec-hi]` - orthogonal wavelet, reconstruction filters are the reversed decomposition filters
  * `[name dec-lo dec-hi rec-lo rec-hi]` - general (biorthogonal) filter bank

  All filters must have the same length, at least 2. Throws `ex-info` otherwise."
  ([nm dec-lo] (let [rec-lo (reverse dec-lo)
                     rec-hi (qmf rec-lo)]
                 (->wavelet nm dec-lo (reverse rec-hi) rec-lo rec-hi)))
  ([nm dec-lo dec-hi] (->wavelet nm dec-lo dec-hi (reverse dec-lo) (reverse dec-hi)))
  ([nm dec-lo dec-hi rec-lo rec-hi]
   (let [lengths (map count [dec-lo dec-hi rec-lo rec-hi])]
     (when-not (apply = lengths)
       (throw (ex-info "All filters should have the same length." {:lengths lengths})))
     (when (m/< (long (first lengths)) 2)
       (throw (ex-info "Filters should have at least 2 coefficients." {:lengths lengths})))
     (->Wavelet nm (count dec-lo)
                (m/seq->double-array dec-lo)
                (m/seq->double-array dec-hi)
                (m/seq->double-array rec-lo)
                (m/seq->double-array rec-hi)))))

(defn coeffs-size
  "Returns the filter length of a wavelet."
  ^long [wv] (prot/coeffs-size wv))

(defn wavelet-name
  "Returns the wavelet name."
  [wv] (prot/wavelet-name wv))

(defn phi
  "Returns the scaling (low-pass) filter, decomposition by default.

  Kind can be:
  
  * `:deconstruction`, `:dec`, `:de` or `:d` - decomposition filter (default)
  * `:reconstruction`, `:rec`, `:re` or `:r` - reconstruction filter"
  ([wv] (prot/phi wv))
  ([wv kind] (prot/phi wv kind)))

(defn psi
  "Returns the wavelet (high-pass) filter, decomposition by default.

  Kind can be:
  
  * `:deconstruction`, `:dec`, `:de` or `:d` - decomposition filter (default)
  * `:reconstruction`, `:rec`, `:re` or `:r` - reconstruction filter"
  ([wv] (prot/psi wv))
  ([wv kind] (prot/psi wv kind)))

;; filter table

(def ^:private filter-table
  (delay (into {} (for [[nm [dl dh rl rh]] (edn/read-string (slurp (io/resource "wavelets/coeffs.edn")))]
                    [nm (->wavelet nm dl dh rl rh)]))))

(def ^:private aliases
  "Names resolved to an entry of the filter table."
  (merge {"la8" "sym4" "la10" "sym5" "la12" "sym6" "la14" "sym7" "la16" "sym8" "la18" "sym9" "la20" "sym10"
          :haar "haar" :haar-orthogonal "haar" :legendre-1 "haar"
          :cdf-53 "bior2.2" :cdf-97 "bior4.4" :battle-23 "battle3" :discrete-mayer "dmey"}
         (into {} (for [n (range 2 21)] [(keyword (str "daubechies-" n)) (str "db" n)]))
         (into {} (for [n (range 2 21)] [(keyword (str "symlet-" n)) (str "sym" n)]))
         (into {} (for [n (range 1 6)] [(keyword (str "coiflet-" n)) (str "coif" n)]))
         (into {} (for [nm ["11" "13" "15" "22" "24" "26" "28" "31" "33" "35" "37" "39" "44" "55" "68"]]
                    [(keyword (str "biorthogonal-" nm)) (str "bior" (subs nm 0 1) "." (subs nm 1))]))))

(def ^{:doc "Names of all available wavelets.

  Strings are the filter table, in PyWavelets / MATLAB naming:

  * from PyWavelets (exact copy): `haar`, `db1`-`db38`, `sym2`-`sym20`, `coif1`-`coif17`, `bior1.1`-`bior6.8`, `rbio1.1`-`rbio6.8`, `dmey` (62 taps)
  * computed by spectral factorisation (double precision): `db39`-`db45`, `sym21`-`sym45` (root choice as in MATLAB)
  * `la8`-`la20` - least asymmetric, equal to `sym4`-`sym10`
  * `battle3` - cubic Battle-Lemarie scaling filter truncated to 161 symmetric taps and one zero (162 taps)
  * MATLAB values copied by hand, lower precision: `beyl`, `bl7`, `bl9`, `bl10`, `fk4`, `fk6`, `fk8`, `fk14`, `fk18`, `fk22`, `han2.3`, `han3.3`, `han4.5`, `han5.5`, `mb4.2`, `mb8.2`, `mb8.3`, `mb8.4`, `mb10.3`, `mb12.3`, `mb14.3`, `mb16.3`, `mb18.3`, `mb24.3`, `mb32.3`, `vaid`

  Perfect reconstruction (`idwt` of `dwt`, any mode, default level) holds to about `1e-14` for most filters. Measured maximum errors (signals of length 301 and 512) above `1e-10`: `dmey` 1.3e-1 (truncated filter, about 1e-2 except in `:smooth` mode; the same in PyWavelets), `mb14.3` 3.9e-4, `mb8.2` 7.5e-6, `mb8.4` 5.9e-6, `mb4.2` 4.1e-6, `mb16.3` 2.4e-6, `mb10.3` 2.2e-6, `mb12.3` 1.4e-6, `mb24.3` 9.1e-7, `mb8.3` 5.4e-7, `mb32.3` 5.2e-7, `mb18.3` 3.6e-7, `vaid` 3.6e-8, `fk8` 2.4e-8, `fk22` 1.4e-8, `bl10` 1.4e-8, `bl9` 5.1e-9, `bl7` 4.9e-10, `sym20` 1.6e-10 (PyWavelets table precision). Above the maximum useful level (filters long compared to the signal) reconstruction is not exact, as in PyWavelets.

  Keywords are former JWave names, aliases of the strings: `:haar`, `:daubechies-2`...`:daubechies-20` (`dbN`), `:symlet-2`...`:symlet-20` (`symN`), `:coiflet-1`...`:coiflet-5`, `:biorthogonal-11`...`:biorthogonal-68` (`biorX.Y`), `:discrete-mayer` (`dmey`), `:cdf-53` (`bior2.2`), `:cdf-97` (`bior4.4`), `:battle-23` (`battle3`), `:haar-orthogonal` and `:legendre-1` (`haar`)."}
  wavelet-names (sort-by str (concat (edn/read-string (slurp (io/resource "wavelets/names.edn"))) (keys aliases))))

(defn wavelet
  "Returns a wavelet (filter bank) object for a name from [[wavelet-names]] (string or keyword); a wavelet object is returned as is.

  Throws `ex-info` for an unknown name."
  [wavelet-name]
  (if (instance? Wavelet wavelet-name)
    wavelet-name
    (let [nm (get aliases wavelet-name wavelet-name)]
      (or (@filter-table nm)
          (throw (ex-info "Unknown wavelet." {:wavelet wavelet-name}))))))

;; engine

(def ^:private extension-modes
  {:zero pad/zero :constant pad/edge :symmetric pad/symmetric :reflect pad/reflect
   :periodic pad/periodic :smooth pad/linear :antisymmetric pad/antisymmetric :antireflect pad/antireflect})

(def ^{:doc "Boundary modes accepted by the transforms."} modes (conj (set (keys extension-modes)) :periodization))

(defn- check-mode!
  [mode]
  (when-not (modes mode)
    (throw (ex-info "Unknown wavelet mode." {:mode mode :allowed (sort modes)}))))

(defn coeffs-length
  "Length of the coefficient bands after one level of [[dwt-level]] for a signal of length `n`, filter length `filter-length` and `mode`."
  ^long [^long n ^long filter-length mode]
  (check-mode! mode)
  (if (= mode :periodization)
    (m/long-quot (m/long-inc n) 2)
    (m/long-quot (m/long-add n (m/long-dec filter-length)) 2)))

(defn max-level
  "Maximum useful decomposition level, as `pywt.dwt_max_level`: the largest `k` with `(filter-length - 1) * 2^k <= n`, 0 when `n < filter-length - 1` or `filter-length < 2`."
  ^long [^long n ^long filter-length]
  (let [l- (m/long-dec filter-length)]
    (if (or (m/< l- 1) (m/< n l-))
      0
      (loop [k 0 v l-]
        (if (m/<= (m/long-mult v 2) n) (recur (m/long-inc k) (m/long-mult v 2)) k)))))

(defn- check-signal!
  [^doubles x mode]
  (let [n (alength x)]
    (when (m/zero? n) (throw (ex-info "Signal is empty." {:length 0})))
    (when (and (m/one? n) (#{:reflect :antireflect} mode))
      (throw (ex-info "Signal length must be greater than 1 for :reflect and :antireflect modes." {:mode mode})))))

(defn dwt-level
  "One level of the discrete wavelet transform of `signal` (a `double-array`), as `pywt.dwt`.

  Returns `[cA cD]`, two `double-array`s of length [[coeffs-length]]."
  [wv ^doubles signal mode]
  (check-mode! mode)
  (check-signal! signal mode)
  (let [w (wavelet wv)
        ^doubles lo (phi w)
        ^doubles hi (psi w)
        l (alength lo)
        n (alength signal)
        out-n (coeffs-length n l mode)
        ca (double-array out-n)
        cd (double-array out-n)]
    (if (= mode :periodization)
      (let [np (m/long-mult 2 out-n)
            x (double-array np)
            half (m/long-quot l 2)]
        (System/arraycopy signal 0 x 0 n)
        (when (m/< n np) (Array/aset x n (Array/aget signal (m/long-dec n))))
        (dotimes [i out-n]
          (loop [j 0 a 0.0 d 0.0]
            (if (m/== j l)
              (do (Array/aset ca i a) (Array/aset cd i d))
              (let [v (Array/aget x (m/long-mod (m/long-sub (m/long-add (m/long-mult 2 i) half) j) np))]
                (recur (m/long-inc j) (m/+ a (m/* v (Array/aget lo j))) (m/+ d (m/* v (Array/aget hi j)))))))))
      (let [^doubles x ((extension-modes mode) signal (m/long-add n (m/long-mult 2 (m/long-dec l))) :both)]
        (dotimes [i out-n]
          (let [k (m/long-add l (m/long-mult 2 i))]
            (loop [j 0 a 0.0 d 0.0]
              (if (m/== j l)
                (do (Array/aset ca i a) (Array/aset cd i d))
                (let [v (Array/aget x (m/long-sub k j))]
                  (recur (m/long-inc j) (m/+ a (m/* v (Array/aget lo j))) (m/+ d (m/* v (Array/aget hi j)))))))))))
    [ca cd]))

(defn idwt-level
  "One level of the inverse discrete wavelet transform, as `pywt.idwt`: reconstructs a signal from bands `ca` and `cd` (`double-array`s of the same length `n`).

  Returns a `double-array` of length `2n - L + 2` (`2n` for `:periodization`), which can be one longer than the original signal; the caller trims it."
  [wv ^doubles ca ^doubles cd mode]
  (check-mode! mode)
  (let [w (wavelet wv)
        ^doubles lo (phi w :rec)
        ^doubles hi (psi w :rec)
        l (alength lo)
        n (alength ca)]
    (when-not (m/== n (alength cd))
      (throw (ex-info "Approximation and detail bands should have the same length." {:approximation n :detail (alength cd)})))
    (if (= mode :periodization)
      (let [m (m/long-mult 2 n)
            out (double-array m)
            shift (m/long-dec (m/long-quot l 2))]
        (dotimes [i n]
          (let [a (Array/aget ca i)
                d (Array/aget cd i)]
            (dotimes [j l]
              (let [k (m/long-mod (m/long-sub (m/long-add (m/long-mult 2 i) j) shift) m)]
                (Array/aset out k (m/+ (Array/aget out k) (m/* a (Array/aget lo j)) (m/* d (Array/aget hi j))))))))
        out)
      (let [m (m/long-add (m/long-sub (m/long-mult 2 n) l) 2)]
        (when-not (m/pos? m)
          (throw (ex-info "Bands are too short for this wavelet." {:length n :filter-length l})))
        (let [out (double-array m)]
          (dotimes [k m]
            (loop [i 0 s 0.0]
              (if (m/== i n)
                (Array/aset out k s)
                (let [idx (m/long-sub (m/long-add k (m/long-sub l 2)) (m/long-mult 2 i))]
                  (recur (m/long-inc i)
                         (if (and (m/not-neg? idx) (m/< idx l))
                           (m/+ s (m/* (Array/aget ca i) (Array/aget lo idx)) (m/* (Array/aget cd i) (Array/aget hi idx)))
                           s))))))
          out)))))

;; scaling and wavelet functions, as pywt.Wavelet.wavefun

(defn- upsample-convolve
  "One full reconstruction step: `c` upsampled by 2 and convolved with `filt`, length `2n + L - 2`."
  ^doubles [^doubles c ^doubles filt]
  (let [n (alength c)
        l (alength filt)
        out (double-array (m/long-sub (m/long-add (m/long-mult 2 n) l) 2))]
    (dotimes [i n]
      (let [v (Array/aget c i)
            base (m/long-mult 2 i)]
        (dotimes [j l]
          (let [k (m/long-add base j)]
            (when (m/< k (alength out))
              (Array/aset out k (m/+ (Array/aget out k) (m/* v (Array/aget filt j)))))))))
    out))

(defn- upcoef
  "`pywt.upcoef` of a single coefficient `c` at `level`: the first step uses the reconstruction low-pass (`approximation?`) or high-pass filter, the others the low-pass one."
  ^doubles [[^doubles rec-lo ^doubles rec-hi] approximation? ^double c ^long level]
  (loop [i 0 coeffs (double-array [c])]
    (if (m/== i level)
      coeffs
      (recur (m/long-inc i) (upsample-convolve coeffs (if (or approximation? (m/pos? i)) rec-lo rec-hi))))))

(defn- biorthogonal-name?
  [nm]
  (and (string? nm) (or (.startsWith ^String nm "bior") (.startsWith ^String nm "rbio"))))

(defn- psi-sign
  "PyWavelets negates psi of a biorthogonal wavelet when the number of vanishing moments of psi (the first digit of `biorX.Y` / `rbioX.Y`) modulo 4 is not 1."
  ^double [nm]
  (if (and (biorthogonal-name? nm)
           (not (m/== 1 (m/long-mod (Long/parseLong (subs nm 4 5)) 4))))
    -1.0 1.0))

(defn- orthogonal?
  [^Wavelet w]
  (and (not (biorthogonal-name? (.name w)))
       (let [^doubles dl (.dec-lo w) ^doubles rl (.rec-lo w)
             ^doubles dh (.dec-hi w) ^doubles rh (.rec-hi w)
             l (alength dl)]
         (every? (fn [^long i] (and (m/== (Array/aget dl i) (Array/aget rl (m/long-sub (m/long-dec l) i)))
                                    (m/== (Array/aget dh i) (Array/aget rh (m/long-sub (m/long-dec l) i)))))
                 (range l)))))

(defn- reversed
  ^doubles [^doubles a]
  (let [l (alength a) out (double-array l)]
    (dotimes [i l] (Array/aset out i (Array/aget a (m/long-sub (m/long-dec l) i))))
    out))

(defn- wavefun-samples
  "Samples `[xs ys]` of phi (`scaling?`) or psi as computed by `pywt.Wavelet.wavefun(level)`."
  [wv ^long level scaling? reconstruction?]
  (when-not (m/pos? level)
    (throw (ex-info "Level must be positive." {:level level})))
  (let [^Wavelet w (wavelet wv)
        orthogonal (orthogonal? w)
        [^doubles lo ^doubles hi] (if (or orthogonal reconstruction?)
                                    [(.rec-lo w) (.rec-hi w)]
                                    [(reversed (.dec-lo w)) (reversed (.dec-hi w))])
        l (alength lo)
        p (m/<< 1 level)
        c (m/* (m/pow m/SQRT2 level) (if (or scaling? orthogonal) 1.0 (psi-sign (.name w))))
        keep-length (loop [i 0 k 1] (if (m/== i level) k (recur (m/long-inc i) (m/long-add (m/long-mult 2 k) (m/long-sub l 2)))))
        output-length (m/long-max (m/long-add keep-length 2)
                             (if orthogonal
                               (m/long-inc (m/long-mult (m/long-dec l) p))
                               (m/long-mult (m/long-dec l) p)))
        ^doubles rec (upcoef [lo hi] scaling? c level)
        left (m/long-quot (m/long-sub (alength rec) keep-length) 2)
        ys (double-array output-length)]
    (System/arraycopy rec (m/long-max 0 left) ys 1 (m/long-min keep-length (alength rec)))
    [(m/slice-range 0.0 (m// (double (m/long-dec output-length)) p) output-length) (vec ys)]))

(defn- samples->function
  [[xs ys]]
  (let [x0 (double (first xs))
        x1 (double (last xs))
        f (il/linear xs ys)]
    (fn ^double [^double x] (if (m/<= x0 x x1) (f x) 0.0))))

(defn scaling-function
  "Returns the scaling function phi of a wavelet, as a function of `x`, computed as `pywt.Wavelet.wavefun(level)`.

  * `wv` - wavelet name or object
  * `level` - refinement level (default 8): phi is sampled at `2^level` points per unit by `level` reconstruction steps of a single coefficient `2^(level/2)` and linearly interpolated
  * `reconstruction?` - for biorthogonal wavelets the reconstruction (`true`) or decomposition (`false`, default) function; ignored for orthogonal wavelets

  The function is 0 outside the PyWavelets grid `[0, (M-1)/2^level]`, `M = max((L-1) 2^level, K + 2)` for biorthogonal and `M = max((L-1) 2^level + 1, K + 2)` for orthogonal wavelets (`L` the filter length, `K` the length of the reconstructed part, `K_0 = 1`, `K_i = 2 K_(i-1) + L - 2`). For `L >= 4` (biorthogonal) or `L >= 3` (orthogonal) the first term wins: `[0, L-1-1/2^level]`, `[0, L-1]`; for shorter filters the grid is `[0, (K+1)/2^level]` (haar: `[0, 1 + 1/2^level]`). Wavelets named `biorX.Y` / `rbioX.Y` are biorthogonal; another wavelet is orthogonal when its reconstruction filters are its reversed decomposition filters.

  Throws `ex-info` for a level smaller than 1."
  ([wv] (scaling-function wv 8))
  ([wv ^long level] (scaling-function wv level false))
  ([wv ^long level reconstruction?]
   (samples->function (wavefun-samples wv level true reconstruction?))))

(defn wavelet-function
  "Returns the wavelet function psi of a wavelet, as a function of `x`, computed as `pywt.Wavelet.wavefun(level)`.

  Arguments as in [[scaling-function]]. Like PyWavelets, psi of a biorthogonal wavelet is negated when the number of vanishing moments of psi (`X` in `biorX.Y`) modulo 4 is not 1."
  ([wv] (wavelet-function wv 8))
  ([wv ^long level] (wavelet-function wv level false))
  ([wv ^long level reconstruction?]
   (samples->function (wavefun-samples wv level false reconstruction?))))
