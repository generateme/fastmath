(ns fastmath.signal.iir
  (:require [fastmath.core :as m]
            [fastmath.complex :as cplx]
            [fastmath.signal.biquad :as biquad])
  (:import [uk.me.berndporr.iirj Bessel Butterworth ChebyshevI ChebyshevII SOSCascade Cascade Biquad 
            DirectFormAbstract]
           [fastmath.signal.biquad BiquadConf]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def direct-forms  {1 DirectFormAbstract/DIRECT_FORM_I
                  2 DirectFormAbstract/DIRECT_FORM_II})

(defn bessel
  (^Cascade [kind ^long order opts] (bessel (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ^double bandwidth ^double cutoff ^double fs ^long direct-form]
              :or {kind :lowpass order 5 bandwidth 0.1 cutoff 0.1 fs 1.0 direct-form 2}}]
   (let [df (direct-forms direct-form)]
     (case kind
       :lowpass  (doto (Bessel.) (.lowPass order fs cutoff (int df)))
       :highpass (doto (Bessel.) (.highPass order fs cutoff (int df)))
       :bandpass (doto (Bessel.) (.bandPass order fs cutoff bandwidth (int df)))
       :bandstop (doto (Bessel.) (.bandStop order fs cutoff bandwidth (int df)))))))

(defn butterworth
  (^Cascade [kind ^long order opts] (butterworth (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ^double bandwidth ^double cutoff ^double fs ^long direct-form]
              :or {kind :lowpass order 5 bandwidth 0.1 cutoff 0.1 fs 1.0 direct-form 2}}]
   (let [df (direct-forms direct-form)]
     (case kind
       :lowpass  (doto (Butterworth.) (.lowPass order fs cutoff (int df)))
       :highpass (doto (Butterworth.) (.highPass order fs cutoff (int df)))
       :bandpass (doto (Butterworth.) (.bandPass order fs cutoff bandwidth (int df)))
       :bandstop (doto (Butterworth.) (.bandStop order fs cutoff bandwidth (int df)))))))

(defn chebyshev
  (^Cascade [kind ^long order opts] (chebyshev (assoc opts :kind kind :order order)))
  (^Cascade [{:keys [kind ^long order ^double ripple type ^double bandwidth ^double cutoff ^double fs
                     ^long direct-form]
              :or {kind :lowpass order 5 ripple 0.5 type :type-I bandwidth 0.1 cutoff 0.1 fs 1.0
                   direct-form 2}}]
   (let [df (direct-forms direct-form)]
     (case type
       :type-I (case kind
                 :lowpass  (doto (ChebyshevI.) (.lowPass order fs cutoff ripple (int df)))
                 :highpass (doto (ChebyshevI.) (.highPass order fs cutoff ripple (int df)))
                 :bandpass (doto (ChebyshevI.) (.bandPass order fs cutoff bandwidth ripple (int df)))
                 :bandstop (doto (ChebyshevI.) (.bandStop order fs cutoff bandwidth ripple (int df))))
       :type-II (case kind
                  :lowpass  (doto (ChebyshevII.) (.lowPass order fs cutoff ripple (int df)))
                  :highpass (doto (ChebyshevII.) (.highPass order fs cutoff ripple (int df)))
                  :bandpass (doto (ChebyshevII.) (.bandPass order fs cutoff bandwidth ripple (int df)))
                  :bandstop (doto (ChebyshevII.) (.bandStop order fs cutoff bandwidth ripple (int df))))))))

(defn sos
  "Create IIR filter from a sequence of SOS coefficients: ([b0 b1 b2 a0 a1 a2],...)."
  (^Cascade [seq-of-coeffs] (sos seq-of-coeffs 2))
  (^Cascade [seq-of-coeffs direct-form]
   (doto (SOSCascade.) (.setup (m/seq->double-double-array seq-of-coeffs) (direct-forms direct-form)))))

;; RBJ https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt

(defn biquad
  (^Cascade [kind opts] (biquad (assoc opts :kind kind)))
  (^Cascade [{:keys [kind ^double cutoff ^double fs direct-form ^double gain ^double bandwidth Q]
              :or {kind :lowpass cutoff 0.1 fs 1.0 direct-form 2 gain 1.0 bandwidth 0.5}}]
   (let [bandwidth (if Q (biquad/Q->bandwidth Q) bandwidth)
         ^BiquadConf bq (case kind
                          :equalizer (biquad/equalizer cutoff gain bandwidth fs)
                          :highshelf (biquad/highshelf cutoff gain (biquad/bandwidth->Q bandwidth) fs)
                          :lowshelf (biquad/lowshelf cutoff gain (biquad/bandwidth->Q bandwidth) fs)
                          :lowpass (biquad/lowpass cutoff bandwidth fs)
                          :highpass (biquad/highpass cutoff bandwidth fs)
                          :bandpass (biquad/bandpass cutoff bandwidth fs)
                          :notch (biquad/notch cutoff bandwidth fs)
                          :allpass (biquad/allpass cutoff bandwidth fs))
         a0 (.a0 bq)]
     (sos [[(m/* (.b0 bq) a0)
            (m/* (.b1 bq) a0)
            (m/* (.b2 bq) a0)
            a0
            (m/* (.a1 bq) a0)
            (m/* (.a2 bq) a0)]] direct-form))))

(defn response
  "Returns complex response spectrum (only first half).

  To convert to DB amplitude use `signal/spectrum` with the following parameters:

  {:method :energy
    :kind :complex
    :db? true
    :domain :frequency}  "
  [^Cascade iir ^long N]
  (let [dN (double N)]
    (map (fn [^long n]
           (cplx/ensure-complex (.response iir (m// n dN)))) (range (m// N 2)))))

(defn sos-coeffs
  "Returns sequence of SOS coefficients: ([b0 b1 b2 a0 a1 a2],...)."
  [^Cascade irr]
  (let [N (.getNumBiquads irr)]
    (for [^int i (range N)
          :let [^Biquad b (.getBiquad irr i)]]
      [(.getB0 b) (.getB1 b) (.getB2 b)
       (.getA0 b) (.getA1 b) (.getA2 b)])))
