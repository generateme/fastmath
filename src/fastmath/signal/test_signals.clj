(ns fastmath.signal.test-signals
  "Collection of test signals"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec3]))

;; https://github.com/JuliaDSP/Wavelets.jl/blob/master/src/mod/Util.jl#L463

(def blocks-coeffs (mapv v/vec2
                       [0.1 0.13 0.15 0.23 0.25 0.4 0.44 0.65 0.76 0.78 0.81]
                       [2.0 -2.5 1.5 -2.0 2.5 -2.1 1.05 2.15 -1.55 1.05 -2.1]))

(defn blocks
  ^double [^double t]
  (let [t (m/frac t)]
    (reduce (fn [^double v ^Vec2 c]
              (m/+ v (m/* (.y c) (m/inc (m/signum (m/- t (.x c))))))) 0.0 blocks-coeffs)))

(def bumps-coeffs (mapv v/vec3
                      [0.1 0.13 0.15 0.23 0.25 0.4 0.44 0.65 0.76 0.78 0.81]
                      [4.0 5.0 3.0 4.0 5.0 4.2 2.1 4.3 3.1 5.1 4.2]
                      [0.005 0.005 0.006 0.01 0.01 0.03 0.01 0.01 0.005 0.008 0.005]))

(defn bumps
  ^double [^double t]
  (let [t (m/frac t)]
    (reduce (fn [^double v ^Vec3 c]
              (m/+ v (m// (.y c)
                          (m/fpow (m/inc (m/abs (m// (m/- t (.x c)) (.z c)))) 4)))) 0.0 bumps-coeffs)))


(defn heavisine
  ^double [^double t]
  (let [t (m/frac t)]
    (m/- (m/* 4.0 (m/sinpi (m/* 4.0 t)))
         (m/signum (m/- t 0.3))
         (m/signum (m/- 0.72 t)))))

(defn doppler
  ^double [^double t]
  (let [t (m/frac t)]
    (m/* (m/sqrt (m/* t (m/- 1.0 t)))
         (m/sinpi (m// 2.1 (m/+ t 0.05))))))


(require '[fastmath.dev.ggplot :as gg])

(gg/->file (gg/function doppler {:x [0 1] :steps 1000}))

