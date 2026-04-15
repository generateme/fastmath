(ns fastmath.signal.biquad
  (:require [fastmath.core :as m]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; https://atparisi.com/html/digitalBiquadFilter.html
;; https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt

;; Store biquad effect configuration in `BiquadConf` type
(deftype BiquadConf [^double b0 ^double b1 ^double b2 ^double a0 ^double a1 ^double a2])

(defn bandwidth->Q ^double [^double bw] (m// (m/* 2.0 (m/sinh (m/* bw m/LN2_2)))))
(defn Q->bandwidth ^double [^double Q] (m// (m/asinh (m// (m// Q) 2.0)) m/LN2_2))

(defn equalizer
  "Calculate configuration for biquad equalizer
   fc - center frequency
   gain
   bw - bandwidth
   fs - sample rate"
  [^double fc ^double gain ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        cw (m/cos omega)
        sw (m/sin omega)
        A (m/exp10 (m/* gain 0.025))
        alpha (m// sw (* 2.0 (bandwidth->Q bw)))
        adivA (m// alpha A)
        amulA (m/* alpha A)

        a0 (m/inc adivA)
        a0r (m// a0)

        b0 (m/* a0r (m/inc amulA))
        b1 (m/* a0r -2.0 cw)
        b2 (m/* a0r (m/- 1.0 amulA))
        a1 b1
        a2 (m/* a0r (m/- 1.0 adivA))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn highshelf
  "Calculate configuration for biquad high shelf
   fc - center frequency
   gain
   Q - quality factor
   fs - sample rate"
  [^double fc ^double gain ^double Q ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        cw (m/cos omega)
        sw (m/sin omega)
        A (m/exp10 (m/* gain 0.025))
        A+ (m/inc A)
        A- (m/dec A)
        alpha (m/* sw (m// (m/sqrt A) Q))
        cA+ (m/* cw A+)
        cA- (m/* cw A-)
        
        a0 (m/- (m/+ A+ alpha) cA-)
        a0r (m// a0)
        
        b0 (m/* a0r A (m/+ A+ cA- alpha))
        b1 (m/* a0r -2.0 A (m/+ A- cA+))
        b2 (m/* a0r A (m/- (m/+ A+ cA-) alpha))
        a1 (m/* a0r 2.0 (m/- A- cA+))
        a2 (m/* a0r (m/- A+ cA- alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn lowshelf
  "Calculate configuration for biquad low shelf
   fc - center frequency
   gain
   Q - quality factor  
   fs - sample rate"
  [^double fc ^double gain ^double Q ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        cw (m/cos omega)
        sw (m/sin omega)
        A (m/exp10 (m/* gain 0.025))
        A+ (m/inc A)
        A- (m/dec A)
        alpha (m/* sw (m// (m/sqrt A) Q))
        cA+ (m/* cw A+)
        cA- (m/* cw A-)
        
        a0 (m/+ A+ cA- alpha)
        a0r (m// a0)
        
        b0 (m/* a0r A (m/- (m/+ A+ alpha) cA-))
        b1 (m/* a0r 2.0 A (m/- A- cA+))
        b2 (m/* a0r A (m/- A+ cA- alpha))
        a1 (m/* a0r -2.0 (m/+ A- cA+))
        a2 (m/* a0r (m/- (m/+ A+ cA-) alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn lowpass
  "Calculate configuration for biquad low pass"
  [^double fc ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        sw (m/sin omega)
        cw (m/cos omega)

        alpha (m// sw (m/* 2.0 (bandwidth->Q bw)))

        a0 (m/inc alpha)
        a0r (m// a0)
        
        cw- (m/- 1.0 cw)
        
        b0 (m/* a0r 0.5 cw-)
        b1 (m/* a0r cw-)
        b2 b0
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn highpass
  "Calculate configuration for biquad high pass"
  [^double fc ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        sw (m/sin omega)
        cw (m/cos omega)
        
        alpha (m// sw (m/* 2.0 (bandwidth->Q bw)))

        a0 (m/inc alpha)
        a0r (m// a0)
        cw+ (m/inc cw)
        
        b0 (m/* a0r 0.5 cw+)
        b1 (m/* a0r (m/- cw+))
        b2 b0
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn bandpass
  "Calculate configuration for biquad band pass"
  [^double fc ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        sw (m/sin omega)
        cw (m/cos omega)
        
        alpha (m// sw (m/* 2.0 (bandwidth->Q bw)))
        
        a0 (m/inc alpha)
        a0r (m// a0)
        
        b0 (m/* a0r alpha)
        b1 0.0
        b2 (m/* a0r (m/- alpha))
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn notch
  "Calculate configuration for biquad notch"
  [^double fc ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        sw (m/sin omega)
        cw (m/cos omega)
        
        alpha (m// sw (m/* 2.0 (bandwidth->Q bw)))
        
        a0 (m/inc alpha)
        a0r (m// a0)
        cw2 (m/* -2.0 cw)
        
        b0 a0r
        b1 (m/* a0r cw2)
        b2 a0r
        a1 b1
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn allpass
  "Calculate configuration for biquad notch"
  [^double fc ^double bw ^double fs]
  (let [omega (m/* m/TWO_PI (m// fc fs))
        sw (m/sin omega)
        cw (m/cos omega)
        
        alpha (m// sw (m/* 2.0 (bandwidth->Q bw)))
        
        a0 (m/inc alpha)
        a0r (m// a0)
        cw2 (m/* -2.0 cw)
        
        b0 (m/* a0r (m/- 1.0 alpha))
        b1 (m/* a0r cw2)
        b2 (m/* a0r (m/inc alpha))
        a1 b1
        a2 b0]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))
