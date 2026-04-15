(ns fastmath.signal.fir
  (:require [fastmath.core :as m]
            [fastmath.complex :as cplx]
            [fastmath.vector :as v]
            [fastmath.kernel :as kernel]
            [fastmath.transform :as t]
            [fastmath.interpolation :as interp]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn firwin-cutoff
  "Prepare bands from cutoff information."
  [kind cutoff ^double fs]
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
                  (concat ncutoff [1.0])))))       ;; __----__----

(defn firwin-scale
  "Calculate the scale factor for firwin."
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
  [window taps opts]
  (cond
    (not window) (repeat taps 1.0)
    (sequential? window) window
    (fn? window) (window taps)
    :else (kernel/window window taps opts)))

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
    * `:p` - Number of cycles of the cutoff frequency used for automatic tap estimation (taps = floor(2 * p * fs / max_cutoff)). Default is `4`.
    * `:scale?` - Whether to normalize the filter coefficients. Can be `true` (auto-scale to unity gain at passband center), `false`, or a specific frequency ratio (0.0 to 1.0) to scale against. Default is `true`.
    * Other keys are passed directly to the window creation function.

  Output:

  Returns a sequence of doubles representing the FIR filter coefficients (impulse response)."
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
     (-> (v/emult h window)
         (firwin-scale m cutoff scale?)))))

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
  "Calculate phases for firgain filter"
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
    * `:nfreqs`: The size of the FFT grid used for the IFFT calculation. Defaults to a power of two slightly larger than the tap count.

  Output:

  Returns a sequence of doubles representing the FIR filter coefficients (impulse response)."
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
                               (firgain-shift xs nyquist taps ftype)))
           w (get-window window taps opts)
           res (mapv m/* coeffs w)]
       (if (m/== ftype 3)
         (assoc res (m// (count res) 2) 0.0)
         res)))))

;;

(defn gammatone
  ""
  [{:keys [^double cutoff ^double fs ^long taps ^long order]
    :or {cutoff 0.1 fs 1.0 order 4}}]
  (let [taps (double (or taps (m/max 15 (long (m/* 0.015 fs)))))
        order- (m/dec order)
        
        ts (v/div (range taps) fs)
        bw (m/* 1.019 (m/+ 24.7 (m// cutoff 9.26449)))
        b (map (fn [^double t]
                 (if (m/zero? t)
                   0.0
                   (let [tt (m/* m/TWO_PI t)]
                     (m/* (m/exp (m/- (m/* order- (m/log t))
                                      (m/* tt bw)))
                          (m/cos (m/* tt cutoff)))))) ts)
        scale (m/exp (m/- (m/+ m/LN2 (m/* order (m/log (m/* m/TWO_PI bw))))
                          (m/log-factorial order-)
                          (m/log fs)))]
    (v/mult b scale)))
