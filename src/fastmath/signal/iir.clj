(ns fastmath.signal.iir
  "IIR filter design, based on the iirj library.

  Filters are iirj `Cascade` objects (cascades of second order sections) with mutable state; process samples with `fastmath.signal/filter-signal`, `fastmath.signal/filter-signal-1` or the `Cascade.filter` method, reset the state with `fastmath.signal/reset-IIR!`.

  * [[butterworth]], [[chebyshev]] (types I and II), [[bessel]] - classic designs, as `scipy.signal.butter`, `cheby1`, `cheby2`, `bessel` (with `norm='phase'`) in second order sections
  * [[biquad]] - one Audio EQ Cookbook section (see `fastmath.signal.biquad`)
  * [[sos]] - a cascade from given section coefficients

  Common options:

  * `:kind` - `:lowpass` (default), `:highpass`, `:bandpass`, `:bandstop`
  * `:order` - filter order (default 5)
  * `:fs` - sample rate (default 1.0)
  * `:cutoff` - corner frequency for low and high pass, center frequency for band filters (default 0.1)
  * `:bandwidth` - width of band filters in frequency units: the pass or stop band is `cutoff ± bandwidth/2` (default 0.1)
  * `:direct-form` - realization of the sections, 1 or 2 (default 2)

  Frequencies are in the units of `:fs` and must be below `fs/2`.

  Throws `ex-info` for an unknown `:kind`, `:type` or `:direct-form`."
  (:require [fastmath.core :as m]
            [fastmath.complex :as cplx]
            [fastmath.signal.biquad :as biquad])
  (:import [uk.me.berndporr.iirj Bessel Butterworth ChebyshevI ChebyshevII SOSCascade Cascade Biquad 
            DirectFormAbstract]
           [fastmath.signal.biquad BiquadConf]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^{:doc "Direct form realizations of the sections: `1` - direct form I, `2` - direct form II."}
  direct-forms {1 DirectFormAbstract/DIRECT_FORM_I
                2 DirectFormAbstract/DIRECT_FORM_II})

(defn- direct-form-constant
  ^long [df]
  (or (direct-forms df)
      (throw (ex-info "Unknown direct form, use 1 or 2." {:direct-form df}))))

(def ^:private kinds [:lowpass :highpass :bandpass :bandstop])

(defn- unknown-kind
  [kind]
  (throw (ex-info "Unknown filter kind." {:kind kind :allowed kinds})))

(defn- positive-order!
  [^long order]
  (when-not (m/pos? order)
    (throw (ex-info "Filter order should be positive." {:order order}))))

(defn bessel
  "Bessel low pass, high pass, band pass or band stop filter, as `scipy.signal.bessel` with `norm='phase'`.

  Options: `:kind`, `:order`, `:fs`, `:cutoff`, `:bandwidth`, `:direct-form`, see the namespace docstring. Returns a `Cascade`."
  (^Cascade [kind ^long order opts] (bessel (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ^double bandwidth ^double cutoff ^double fs direct-form]
              :or {kind :lowpass order 5 bandwidth 0.1 cutoff 0.1 fs 1.0 direct-form 2}}]
   (positive-order! order)
   (let [df (direct-form-constant direct-form)]
     (case kind
       :lowpass  (doto (Bessel.) (.lowPass order fs cutoff (int df)))
       :highpass (doto (Bessel.) (.highPass order fs cutoff (int df)))
       :bandpass (doto (Bessel.) (.bandPass order fs cutoff bandwidth (int df)))
       :bandstop (doto (Bessel.) (.bandStop order fs cutoff bandwidth (int df)))
       (unknown-kind kind)))))

(defn butterworth
  "Butterworth low pass, high pass, band pass or band stop filter, as `scipy.signal.butter`.

  Options: `:kind`, `:order`, `:fs`, `:cutoff`, `:bandwidth`, `:direct-form`, see the namespace docstring. Returns a `Cascade`."
  (^Cascade [kind ^long order opts] (butterworth (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ^double bandwidth ^double cutoff ^double fs direct-form]
              :or {kind :lowpass order 5 bandwidth 0.1 cutoff 0.1 fs 1.0 direct-form 2}}]
   (positive-order! order)
   (let [df (direct-form-constant direct-form)]
     (case kind
       :lowpass  (doto (Butterworth.) (.lowPass order fs cutoff (int df)))
       :highpass (doto (Butterworth.) (.highPass order fs cutoff (int df)))
       :bandpass (doto (Butterworth.) (.bandPass order fs cutoff bandwidth (int df)))
       :bandstop (doto (Butterworth.) (.bandStop order fs cutoff bandwidth (int df)))
       (unknown-kind kind)))))

(defn chebyshev
  "Chebyshev type I or II low pass, high pass, band pass or band stop filter, as `scipy.signal.cheby1` / `cheby2`.

  Options as in the namespace docstring, plus:

  * `:type` - `:type-I` (default, ripple in the pass band) or `:type-II` (ripple in the stop band)
  * `:ripple` - in dB: pass band ripple for type I (`rp`, default 0.5), minimum stop band attenuation for type II (`rs`, default 40)

  Returns a `Cascade`."
  (^Cascade [kind ^long order opts] (chebyshev (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ripple type ^double bandwidth ^double cutoff ^double fs direct-form]
              :or {kind :lowpass order 5 type :type-I bandwidth 0.1 cutoff 0.1 fs 1.0 direct-form 2}}]
   (positive-order! order)
   (let [df (direct-form-constant direct-form)]
     (case type
       :type-I (let [ripple (double (or ripple 0.5))]
                 (case kind
                   :lowpass  (doto (ChebyshevI.) (.lowPass order fs cutoff ripple (int df)))
                   :highpass (doto (ChebyshevI.) (.highPass order fs cutoff ripple (int df)))
                   :bandpass (doto (ChebyshevI.) (.bandPass order fs cutoff bandwidth ripple (int df)))
                   :bandstop (doto (ChebyshevI.) (.bandStop order fs cutoff bandwidth ripple (int df)))
                   (unknown-kind kind)))
       :type-II (let [ripple (double (or ripple 40.0))]
                  (case kind
                    :lowpass  (doto (ChebyshevII.) (.lowPass order fs cutoff ripple (int df)))
                    :highpass (doto (ChebyshevII.) (.highPass order fs cutoff ripple (int df)))
                    :bandpass (doto (ChebyshevII.) (.bandPass order fs cutoff bandwidth ripple (int df)))
                    :bandstop (doto (ChebyshevII.) (.bandStop order fs cutoff bandwidth ripple (int df)))
                    (unknown-kind kind)))
       (throw (ex-info "Unknown Chebyshev type." {:type type :allowed [:type-I :type-II]}))))))

(defn sos
  "Creates an IIR filter (`Cascade`) from a sequence of second order sections `([b0 b1 b2 a0 a1 a2] ...)`, as rows of `scipy.signal` `sos` arrays.

  `direct-form` is 1 or 2 (default)."
  (^Cascade [seq-of-coeffs] (sos seq-of-coeffs 2))
  (^Cascade [seq-of-coeffs df]
   (doto (SOSCascade.) (.setup (m/seq->double-double-array seq-of-coeffs) (int (direct-form-constant df))))))

;; RBJ https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt

(def ^:private biquad-kinds [:equalizer :highshelf :lowshelf :lowpass :highpass :bandpass :notch :allpass])

(defn biquad
  "One Audio EQ Cookbook biquad section as a `Cascade`, see `fastmath.signal.biquad`.

  Options:

  * `:kind` - `:lowpass` (default), `:highpass`, `:bandpass`, `:notch`, `:allpass`, `:equalizer` (peaking), `:lowshelf`, `:highshelf`
  * `:cutoff` - corner or center frequency (default 0.1), `:fs` - sample rate (default 1.0)
  * `:bandwidth` - in octaves (default 0.5), or `:Q` - quality factor (converted with the digital cookbook formula)
  * `:gain` - in dB for `:equalizer` and shelves (default 1.0)
  * `:direct-form` - 1 or 2 (default)"
  (^Cascade [kind opts] (biquad (assoc opts :kind kind)))
  (^Cascade [{:keys [kind ^double cutoff ^double fs direct-form ^double gain ^double bandwidth Q]
              :or {kind :lowpass cutoff 0.1 fs 1.0 direct-form 2 gain 1.0 bandwidth 0.5}}]
   (let [bandwidth (if Q (biquad/Q->bandwidth Q cutoff fs) bandwidth)
         ^BiquadConf bq (case kind
                          :equalizer (biquad/equalizer cutoff gain bandwidth fs)
                          :highshelf (biquad/highshelf cutoff gain (biquad/bandwidth->Q bandwidth cutoff fs) fs)
                          :lowshelf (biquad/lowshelf cutoff gain (biquad/bandwidth->Q bandwidth cutoff fs) fs)
                          :lowpass (biquad/lowpass cutoff bandwidth fs)
                          :highpass (biquad/highpass cutoff bandwidth fs)
                          :bandpass (biquad/bandpass cutoff bandwidth fs)
                          :notch (biquad/notch cutoff bandwidth fs)
                          :allpass (biquad/allpass cutoff bandwidth fs)
                          (throw (ex-info "Unknown biquad kind." {:kind kind :allowed biquad-kinds})))
         a0 (.a0 bq)]
     (sos [[(m/* (.b0 bq) a0)
            (m/* (.b1 bq) a0)
            (m/* (.b2 bq) a0)
            a0
            (m/* (.a1 bq) a0)
            (m/* (.a2 bq) a0)]] direct-form))))

(defn response
  "Complex frequency response of the filter at `N/2` frequencies `k fs / N`, `k = 0 ... N/2 - 1` (from 0 to below the Nyquist frequency), as `scipy.signal.sosfreqz` at these frequencies.

  Returns a vector of complex numbers (`Vec2`)."
  [^Cascade iir ^long N]
  (let [dN (double N)]
    (mapv (fn [^long n]
            (cplx/ensure-complex (.response iir (m// n dN)))) (range (m/long-quot N 2)))))

(defn sos-coeffs
  "Returns the second order sections of a filter: `([b0 b1 b2 a0 a1 a2] ...)`, as rows of a `scipy.signal` `sos` array."
  [^Cascade irr]
  (let [N (.getNumBiquads irr)]
    (vec (for [^int i (range N)
               :let [^Biquad b (.getBiquad irr i)]]
           [(.getB0 b) (.getB1 b) (.getB2 b)
            (.getA0 b) (.getA1 b) (.getA2 b)]))))
