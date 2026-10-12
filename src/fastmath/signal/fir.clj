(ns fastmath.signal.fir
  "FIR (finite impulse response) filter design.

  * [[firwin]] - window method for low pass, high pass, band pass and band stop filters, as `scipy.signal.firwin`, except that an even number of taps for a filter passing the Nyquist frequency is increased by one (scipy throws `ValueError`)
  * [[firgain]] - frequency sampling method for an arbitrary piecewise gain, as `scipy.signal.firwin2`
  * [[gammatone]] - gammatone filter of auditory models, as `scipy.signal.gammatone` with `ftype='fir'`

  Filters are vectors of coefficients (impulse responses); apply them with `fastmath.signal/filter-signal` or `fastmath.signal/convolve`. Windows come from `fastmath.kernel/window` (symmetric by default); the default window is `:hann` (scipy uses `hamming` in `firwin`)."
  (:require [fastmath.core :as m]
            [fastmath.complex :as cplx]
            [fastmath.vector :as v]
            [fastmath.kernel :as kernel]
            [fastmath.transform :as t]
            [fastmath.interpolation :as interp]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^:private kinds [:lowpass :highpass :bandpass :bandstop])

(defn firwin-cutoff
  "Converts cutoff frequencies to the band edges used by [[firwin]], normalized to the Nyquist frequency `fs/2`.

  `cutoff` is one frequency or a sequence of frequencies (sorted, duplicates removed), each in `(0, fs/2)`. Returns a flat sequence of pass band edges `[low1 high1 low2 high2 ...]` in `[0, 1]`:

  * `:lowpass` - `[0 c]`, `c` the largest cutoff
  * `:highpass` - `[c 1]`, `c` the smallest cutoff
  * `:bandpass` - pass bands between consecutive cutoffs, starting with a pass band at the first cutoff; an odd count passes up to the Nyquist frequency
  * `:bandstop` - stop bands between consecutive cutoffs, starting with a pass band at 0; an odd count stops up to the Nyquist frequency

  Throws `ex-info` for an empty cutoff, a cutoff outside `(0, fs/2)` or an unknown `kind`."
  [kind cutoff ^double fs]
  (when (and (sequential? cutoff) (empty? cutoff))
    (throw (ex-info "Cutoff can't be empty." {:cutoff cutoff})))
  (let [nyquist (m/* 0.5 fs)
        ncutoff (-> (if (number? cutoff) [cutoff] (sort cutoff))
                    (distinct)
                    (v/div nyquist))
        elen? (m/even? (count ncutoff))]
    (when-not (every? (fn [^double c] (m/< 0 c 1.0)) ncutoff)
      (throw (ex-info "Cutoff values should be 0 < v < fs/2" {:cutoff cutoff})))
    (case kind
      :lowpass [0.0 (reduce max ncutoff)] ;; ------______
      :highpass [(reduce min ncutoff) 1.0] ;; ______------
      :bandstop (if elen?
                  (concat [0.0] ncutoff [1.0]) ;; --____------
                  (concat [0.0] ncutoff))      ;; --____--____
      :bandpass (if elen?
                  ncutoff                    ;; __----______
                  (concat ncutoff [1.0]))    ;; __----__----
      (throw (ex-info "Unknown filter kind." {:kind kind :allowed kinds})))))

(defn firwin-scale
  "Scales filter coefficients `h` (at positions `m`, centered at 0) to unit gain at one frequency, as `scipy.signal.firwin` with `scale=True`.

  * `cutoff` - normalized band edges from [[firwin-cutoff]]; the frequency is 0 when the first band starts at 0, 1 (Nyquist) when it ends at 1, otherwise the center of the first band
  * `scale?` - `false` returns `h` unchanged, `true` uses the frequency above, a number is the normalized frequency (0 to 1) to use"
  [h m cutoff scale?]
  (if (not scale?)
    h
    (let [sf (double (if (number? scale?)
                       scale?
                       (let [[^double low ^double high] cutoff]
                         (cond
                           (m/zero? low) 0.0
                           (m/one? high) 1.0
                           :else (m/* 0.5 (m/+ low high))))))
          s (if (m/zero? sf)
              (v/sum h)
              (->> (m/* m/PI sf)
                   (v/mult m)
                   (v/cos)
                   (v/emult h)
                   (v/sum)))]
      (v/div h s))))

(defn- get-window
  [window ^long taps opts]
  (let [w (cond
            (not window) (repeat taps 1.0)
            (sequential? window) window
            (fn? window) (window taps)
            :else (kernel/window window taps opts))]
    (when-not (m/== taps (count w))
      (throw (ex-info "Window length should equal the number of taps." {:taps taps :window-length (count w)})))
    w))

(defn firwin
  "Designs a Finite Impulse Response (FIR) filter using the window method.

  This function constructs a linear-phase FIR filter by calculating the ideal impulse response (sinc function) for specified frequency bands and tapering it with a window function to reduce spectral leakage (Gibbs phenomenon). It provides a robust way to create standard filter types—low-pass, high-pass, band-pass, and band-stop—with controlled transition regions.

  Input parameters:

  * `taps` - The number of filter coefficients (length of the filter). If `nil`, it is automatically estimated based on the `:p` parameter and the highest cutoff frequency. If the filter includes the Nyquist frequency, an even tap count is automatically incremented to preserve symmetry.
  * `kind` - A keyword specifying the filter type: `:lowpass`, `:highpass`, `:bandstop`, or `:bandpass`.
  * `opts` - A map of configuration options:
    * `:fs` - Sampling frequency. Default is `1.0`.
    * `cutoff` - A single frequency or a sequence of frequencies defining the filter transitions. Values must be between 0 and `fs/2` (not inclusive). Default is `0.1`.
    * `:window` - The windowing strategy. Can be a keyword (e.g., `:hann`, `:hamming`, `:blackman`), a window function, or a sequence of pre-calculated coefficients. Default is `:hann`.
    * `:p` - Number of cycles of the cutoff frequency used for automatic tap estimation `taps = floor(p * fs / max_cutoff)` (with cutoffs below `fs/2`). Default is `4`.
    * `:scale?` - Whether to normalize the filter coefficients. Can be `true` (auto-scale to unity gain at passband center), `false`, or a specific frequency ratio (0.0 to 1.0) to scale against. Default is `true`.
    * Other keys are passed directly to the window creation function.

  For `:lowpass` with several cutoffs the largest is used, for `:highpass` the smallest.

  Output:

  Returns a vector of doubles representing the FIR filter coefficients (impulse response). Equals `scipy.signal.firwin(taps, cutoff, window, pass_zero, scale, fs)` for the same window, except for an even `taps` with a filter passing the Nyquist frequency: it is increased by one here, scipy throws `ValueError`.

  Throws `ex-info` for an empty or out of range cutoff, an unknown `kind`, or a window sequence of another length than `taps`."
  ([taps kind opts] (firwin (assoc opts :taps taps :kind kind)))
  ([{:keys [taps kind cutoff ^double fs window ^long p scale?]
     :or {kind :lowpass cutoff 0.1 fs 1.0 window :hann p 4 scale? true}
     :as opts}]
   (let [cutoff (firwin-cutoff kind cutoff fs)
         taps (long (or taps
                        (when (sequential? window) (count window))
                        (m/floor (m// (m/* 2.0 p) (->> (filter (fn [^double v] (m/< 0.0 v 1.0)) cutoff)
                                                       (apply m/max)
                                                       (double))))))
         taps (if (and (m/even? taps) (m/one? (double (last cutoff)))) (m/inc taps) taps)
         m (v/shift (range taps) (m/* -0.5 (m/dec taps)))
         h (reduce (fn [buff [^double low ^double high]]
                     (let [curr (v/sub (v/mult (v/sinc (v/mult m high)) high)
                                       (v/mult (v/sinc (v/mult m low)) low))]
                       (if buff (v/add buff curr) curr))) nil (partition 2 cutoff))
         window (get-window window taps opts)]
     (vec (-> (v/emult h window)
              (firwin-scale m cutoff scale?))))))

;;

(defn- validate-firgain
  [freqs gains ^double nyq ^long ftype]
  (when-not (and (m/zero? (double (first freqs)))
                 (m/== nyq (double (last freqs))))
    (throw (ex-info "First frequency should be 0.0 and last should be nyquist=fs/2"
                    {:nyquist nyq
                     :first (first freqs) :last (last freqs)})))
  (when-not (= freqs (distinct freqs))
    (throw (ex-info "Frequencies can't contain duplicates"
                    {:frequencies freqs})))
  (let [g0 (double (first gains))
        gn (double (last gains))]
    (condp m/== ftype
      2 (when-not (m/zero? gn)
          (throw (ex-info "Type II filter must have 0 at nyquist frequency"
                          {:gain-nyquist gn})))
      3 (when-not (and (m/zero? g0) (m/zero? gn))
          (throw (ex-info "Type III filter must have 0 at DC and nyquist frequency"
                          {:gain-DC g0 :gain-nyquist gn})))
      4 (when-not (m/zero? g0)
          (throw (ex-info "Type IV filter must have 0 at DC"
                          {:gain-DC g0})))
      nil)))

(defn firgain-shift
  "Linear phase terms of a [[firgain]] filter at frequencies `xs`: `exp(-i pi/2 (taps-1) x / nyquist)`, multiplied by `i` for antisymmetric filters (types III and IV, `ftype` 3 or 4)."
  [xs ^double nyquist ^long taps ^long ftype]
  (let [fac (m// (m/* m/-HALF_PI (m/dec taps)) nyquist)
        res (map (fn [^double x]
                   (cplx/exp (cplx/complex 0.0 (m/* fac x)))) xs)]
    (if (m/> ftype 2)
      (map (fn [z] (cplx/mult-I z)) res)
      res)))

(defn firgain
  "Design an FIR filter from a specified frequency-gain response using the frequency sampling method.

  This function constructs an FIR filter by interpolating a set of desired frequency-gain points to create a continuous frequency response. It automatically handles the design of the four types of linear-phase FIR filters based on the number of taps and the symmetry requirements.

  Parameters:

  * `taps`: The number of filter coefficients to generate.
  * `freq-gain-pairs`: A sequence of `[frequency, gain]` pairs defining the desired filter shape. The first frequency must be 0.0 (DC) and the last must be the Nyquist frequency (`fs/2`).
  * `opts`: A map of optional configuration:
    * `:fs`: Sampling frequency. Default is `1.0`.
    * `:window`: Window function name (keyword), function, or sequence of coefficients. Default is `:hann`.
    * `:antisymmetric?`: Boolean. If `true`, designs an antisymmetric (Type III or IV) filter; otherwise, designs a symmetric (Type I or II) filter. Default is `false`.
    * `:interpolator`: The interpolation method used to bridge the provided frequency points (e.g., `:linear`, `:cubic`).
    * `:nfreqs`: Number of points of the frequency grid from 0 to Nyquist used for the inverse FFT (at least `taps + 1`). Default `1 + 2^ceil(log2(1.25 taps))`; `scipy.signal.firwin2` uses `1 + 2^ceil(log2(taps))`, so results differ slightly unless `:nfreqs` is given.

  Output:

  Returns a vector of doubles representing the FIR filter coefficients (impulse response). With the same `:nfreqs` and window it equals `scipy.signal.firwin2(taps, freqs, gains, nfreqs, window, antisymmetric, fs)`.

  Throws `ex-info` when the frequencies don't start at 0 and end at Nyquist, contain duplicates, or the gains violate the type constraints (type II: 0 at Nyquist, type III: 0 at DC and Nyquist, type IV: 0 at DC)."
  ([^long taps freq-gain-pairs opts] (firgain (assoc opts :taps taps :freq-gain-pairs freq-gain-pairs)))
  ([{:keys [^long taps freq-gain-pairs ^double fs window antisymmetric? interpolator ^long nfreqs]
     :or {taps 15 fs 1.0 window :hann antisymmetric? false
          interpolator :linear}
     :as opts}]
   (let [ftype (if (m/even? taps)
                 (if antisymmetric? 4 2)
                 (if antisymmetric? 3 1))
         fg (sort-by first freq-gain-pairs)
         nyquist (m/* 0.5 fs)
         freqs (mapv first fg)
         gains (mapv second fg)]
     (validate-firgain freqs gains nyquist ftype)
     (let [interpolator (if (fn? interpolator)
                          (interpolator freqs gains)
                          (interp/interpolation interpolator freqs gains))
           nfreqs (max (m/inc taps)
                       (double (or nfreqs (m/inc (m/exp2 (m/ceil (m/log2 (m/* 1.25 taps))))))))
           xs (m/slice-range 0.0 nyquist nfreqs)
           coeffs (t/ifft (map (fn [^double v z]
                                 (cplx/scale z v))
                               (map interpolator xs)
                               (firgain-shift xs nyquist taps ftype))
                          {:kind :real :n (m/long-mult 2 (m/long-dec (count xs)))})
           w (get-window window taps opts)
           res (mapv m/* coeffs w)]
       (if (m/== ftype 3)
         (assoc res (m// (count res) 2) 0.0)
         res)))))

;;

(defn gammatone
  "FIR gammatone filter (auditory filter model), as `scipy.signal.gammatone(cutoff, 'fir', order, taps, fs)`.

  The impulse response `t^(order-1) exp(-2 pi b t) cos(2 pi cutoff t)`, with the equivalent rectangular bandwidth `b = 1.019 (24.7 + cutoff / 9.26449)` in Hz, is sampled at `t = n / fs` and normalized.

  Options:

  * `:cutoff` - center frequency in Hz, `0 < cutoff < fs/2` (default 1000.0)
  * `:fs` - sample rate in Hz (default 16000.0)
  * `:order` - filter order, at least 1 (default 4)
  * `:taps` - number of coefficients (default `max(15, 0.015 fs)`)

  Frequencies are in Hz because the bandwidth formula is. Returns a vector of doubles. Throws `ex-info` for an invalid `:cutoff` or `:order`."
  [{:keys [^double cutoff ^double fs ^long taps ^long order]
    :or {cutoff 1000.0 fs 16000.0 order 4}}]
  (when-not (m/< 0.0 cutoff (m/* 0.5 fs))
    (throw (ex-info "Cutoff should be in (0, fs/2)." {:cutoff cutoff :fs fs})))
  (when-not (m/pos? order)
    (throw (ex-info "Order should be at least 1." {:order order})))
  (let [taps (double (or taps (m/max 15 (long (m/* 0.015 fs)))))
        order- (m/dec order)
        
        ts (v/div (range taps) fs)
        bw (m/* 1.019 (m/+ 24.7 (m// cutoff 9.26449)))
        b (map (fn [^double t]
                 (if (m/zero? t)
                   (if (m/zero? order-) 1.0 0.0) ;; t^(order-1) at t = 0, 0^0 = 1
                   (let [tt (m/* m/TWO_PI t)]
                     (m/* (m/exp (m/- (m/* order- (m/log t))
                                      (m/* tt bw)))
                          (m/cos (m/* tt cutoff)))))) ts)
        scale (m/exp (m/- (m/+ m/LN2 (m/* order (m/log (m/* m/TWO_PI bw))))
                          (m/log-factorial order-)
                          (m/log fs)))]
    (vec (v/mult b scale))))
