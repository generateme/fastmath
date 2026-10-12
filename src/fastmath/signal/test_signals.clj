(ns fastmath.signal.test-signals
  "Donoho-Johnstone test signals for denoising and wavelet analysis (WaveLab `MakeSignal`, also in Julia `Wavelets.jl`).

  Every function takes time `t` and is periodic with period 1 (the value at `t` is the value at the fractional part `t - floor(t)`, so the signal at `t = 1` equals the one at `t = 0`). WaveLab samples `t = (1..n)/n`."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; https://github.com/JuliaDSP/Wavelets.jl/blob/master/src/mod/Util.jl#L463

(defn- cycle-frac
  "Fractional part in [0, 1), periodic also for negative values."
  ^double [^double t]
  (m/- t (m/floor t)))

(def ^{:doc "Jump positions (`x`) and half heights (`y`) of [[blocks]]."}
  blocks-coeffs (mapv v/vec2
                      [0.1 0.13 0.15 0.23 0.25 0.4 0.44 0.65 0.76 0.78 0.81]
                      [2.0 -2.5 1.5 -2.0 2.5 -2.1 1.05 2.15 -1.55 1.05 -2.1]))

(defn blocks
  "Piecewise constant signal with 11 jumps: `sum h_j (1 + sign(t - p_j))`, see [[blocks-coeffs]]."
  ^double [^double t]
  (let [t (cycle-frac t)]
    (reduce (fn [^double v ^Vec2 c]
              (m/+ v (m/* (.y c) (m/inc (m/signum (m/- t (.x c))))))) 0.0 blocks-coeffs)))

(def ^{:doc "Positions (`x`), heights (`y`) and widths (`z`) of [[bumps]]."}
  bumps-coeffs (mapv v/vec3
                     [0.1 0.13 0.15 0.23 0.25 0.4 0.44 0.65 0.76 0.78 0.81]
                     [4.0 5.0 3.0 4.0 5.0 4.2 2.1 4.3 3.1 5.1 4.2]
                     [0.005 0.005 0.006 0.01 0.01 0.03 0.01 0.01 0.005 0.008 0.005]))

(defn bumps
  "Sum of 11 sharp bumps: `sum h_j / (1 + |t - p_j| / w_j)^4`, see [[bumps-coeffs]]."
  ^double [^double t]
  (let [t (cycle-frac t)]
    (reduce (fn [^double v ^Vec3 c]
              (m/+ v (m// (.y c)
                          (m/fpow (m/inc (m/abs (m// (m/- t (.x c)) (.z c)))) 4)))) 0.0 bumps-coeffs)))

(defn heavisine
  "Sine with two jumps: `4 sin(4 pi t) - sign(t - 0.3) - sign(0.72 - t)`."
  ^double [^double t]
  (let [t (cycle-frac t)]
    (m/- (m/* 4.0 (m/sinpi (m/* 4.0 t)))
         (m/signum (m/- t 0.3))
         (m/signum (m/- 0.72 t)))))

(defn doppler
  "Chirp-like signal with decreasing frequency: `sqrt(t (1 - t)) sin(2 pi 1.05 / (t + 0.05))`."
  ^double [^double t]
  (let [t (cycle-frac t)]
    (m/* (m/sqrt (m/* t (m/- 1.0 t)))
         (m/sinpi (m// 2.1 (m/+ t 0.05))))))
