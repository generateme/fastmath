(ns fastmath.signal-examples
  (:require [metadoc.examples :refer :all]
            [fastmath.signal :refer :all]
            [fastmath.core :as m]))

(add-examples compose-effects
  (example "Usage" (let [lpf (effect :simple-lowpass)
                         hpf(effect :simple-highpass {:cutoff 2000.0})
                         composed (compose-effects lpf hpf)
                         result (composed 0.5)]
                     {"Composed effects" composed
                      "After call on sample" result
                      "Extracted sample" (result)})))

(add-examples apply-effects
  (example "Usage" (let [lpf (effect :simple-lowpass {:cutoff 10000})
                         sgnal [-1.0 1.0 -0.5 0.5 -0.1 0.1 0 0]]
                     (apply-effects sgnal lpf)))
  (example "Reset state every two samples" (let [lpf (effect :simple-lowpass {:cutoff 10000})
                                                 sgnal [-1.0 1.0 -0.5 0.5 -0.1 0.1 0 0]]
                                             (apply-effects sgnal lpf 2))))

(add-examples apply-effects-raw
  (example "Usage" (let [lpf (effect :simple-lowpass {:cutoff 10000})
                         sgnal (m/seq->double-array [-1.0 1.0 -0.5 0.5 -0.1 0.1 0 0])]
                     (apply-effects-raw sgnal lpf)))
  (example "Reset state every two samples" (let [lpf (effect :simple-lowpass {:cutoff 10000})
                                                 sgnal (m/seq->double-array [-1.0 1.0 -0.5 0.5 -0.1 0.1 0 0])]
                                             (seq (apply-effects-raw sgnal lpf 2)))))

(add-examples db->linear (example (db->linear 0.5)))
(add-examples linear->db (example (linear->db 0.5)))

(add-examples effect
  (example-session "Basic usage" (effect :fm) (effect :fm {:quant 5}))
  (example "Use as a function" (let [fm (effect :fm)]
                                 (fm 0.5)))
  (example "Extract last result." (let [fm (effect :fm)
                                        result (fm 0.5)]
                                    (result))))

(add-examples effects-list
  (example "List of all effects" effects-list))

(add-examples save-signal
  (example (let [s (oscillator->signal (oscillator :sin 1 1 0) 11050 5)]
             (save-signal s "signal.raw"))))

(add-examples load-signal
  (example (let [s (oscillator->signal (oscillator :sin 1 1 0) 11050 5)]
             (save-signal s "signal.raw")
             {:out (nth (load-signal "signal.raw") 20)
              :in (nth s 20)})))

(add-examples oscillator
  (example-session "Usage"
    (let [wave-sin (oscillator :sin 1.0 1.0 0.5)]
      (wave-sin 0.1))
    (let [wave-noise (oscillator :triangle 0.1 1.0 0.123)]
      (wave-noise 0.1))))

(add-examples oscillator-gain
  (example
    (let [wave-sin (oscillator :sin 1.0 1.0 0.5)
          gained (oscillator-gain wave-sin 0.22)]
      {:original (wave-sin 0.75)
       :gained (gained 0.75)})))


(add-examples oscillators
  (example "List of oscillators" oscillators))

(add-examples oscillators-sum
  (example (let [w1 (oscillator :triangle 1.5 0.5 0.5)
                 w2 (oscillator :sin 1 0.5 0)
                 sum (oscillators-sum w1 w2)]
             (sum 0.5)))
  (example-image "Signal plot" "images/s/sum.jpg"))

(add-examples oscillator->signal
  (example (let [w (oscillator :sin 0.2 0.9 0.2)]
             (oscillator->signal w 44100 5))))

(add-examples signal->oscillator
  (example (let [w (oscillator :sin 0.2 0.9 0.2)
                 w2 (signal->oscillator (oscillator->signal w 44100 5) 5)]
             {:original (w 0.5)
              :converted (w2 0.5)})))

(defmacro oscillator-plots
  []
  `(add-examples oscillator
     ~@(for [o (disj (set oscillators) :constant)]
         (example-image (str "Plot of " (name o)) (str "images/s/" (name o) ".jpg")))))

(oscillator-plots)

(defmacro effect-plots
  []
  `(add-examples effect
     ~@(for [o effects-list]
         (example-image (str "Plot of " (name o)) (str "images/s/" (name o) ".jpg")))))

(effect-plots)

(add-examples savgol-filter
  (example (let [signal [-1 2 4 99 4 2 -1]
                 savgol (savgol-filter)
                 savgol-1st-deriv (savgol-filter 5 2 1)]
             {:smoothed (savgol signal)
              :deriv (savgol-1st-deriv signal)}))
  (example-image "Savitzky-Golay smoothing" "images/s/savgol.jpg"))
