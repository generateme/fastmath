^:kindly/hide-code
(ns transform
  (:require [fastmath.transform :as t]

            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]

            [fastmath.core :as m]
            [fastmath.random :as r]))

(gg/->file (gg/functions [["basic" m/sin]
                          ["noisy" (fn [x] (+ (m/sin x) (* 0.2 (- (rand) 0.5))))]]
                         {:x [m/-TWO_PI m/TWO_PI]
                          :ylim [-2 2]
                          :steps 500
                          :palette gg/palette-blue-0}))

(gg/->file (gg/function2d (fn [[x y]] (m/sin (m/* x (m/cos y)))) {:x [m/-TWO_PI m/TWO_PI]
                                                                 :y [m/-TWO_PI m/TWO_PI]
                                                                 :title "sin(x*cos(y))"
                                                                 :legend-name "value"}))

(gg/function2d (fn [[x y]] (m/sin (m/* x (m/cos y)))) {:x [m/-TWO_PI m/TWO_PI]
                                                      :y [m/-TWO_PI m/TWO_PI]
                                                      :title "sin(x*cos(y))"
                                                      :legend-name "value"})

(gg/function m/tan {:x [m/-TWO_PI m/TWO_PI]
                    :title "tan(x)"
                    :ylab "y=tan(x)"
                    :ylim [-2 2] ;; we need to limit y axis
                    :steps 500})

(let [xs (repeatedly 2000 r/grand)
      ys (map (fn [x] (+ (r/grand (+ 0.1 (* x 0.5))) (m/sin (* 2 x)))) xs)]
  (gg/scatter xs ys {:title "Scatter"}))

(let [xy (take 1000 (r/sequence-generator :r2 2))]
  (gg/scatter xy nil {:title "R2 low-discrepancy sequence generator"}))

(gg/functions [["tan" m/tan]
               ["cot" m/cot]
               ["sin" m/sin]
               ["cos" m/cos]]
              {:x [m/-TWO_PI m/TWO_PI]
               :title "Basic trigonometric functions"
               :ylim [-2 2]
               :steps 500
               :palette gg/palette-blue-1})

(gg/->file (gg/function m/sec
                        {:x [m/-TWO_PI m/TWO_PI]
                         :ylim [-2 2]
                         :steps 500}))


;; # Transforms {.unnumbered}

;; General description of the topic

;; ::: {.callout-tip title="Defined functions"}
;; * `transformer`
;; * `forward-1d`, `forward-2d`
;; * `reverse-1d`, `reverse-2d`
;; :::

;; ## FFT

;; Details about FFT and use-cases

;; Some examples:

(def fft-real (t/transformer :real :fft ))

(utls/examples-note
  (seq (t/forward-1d fft-real [1 2 -10 1]))
  (seq (t/reverse-1d fft-real [-6 -12 11 -1])))

;; ## Wavelets

;; ## Compression and denoising

;; An use case with charts

(def domain (m/slice-range 0 10 512))
(def signal (map (fn [x] (+ (Math/sin x)
                         (* 0.1 (- (rand) 0.5)))) ;; add some noise
               domain))
(def denoised-signal (t/denoise fft-real signal {:method :hard}))

^:kind/table
[[(gg/->image (gg/line domain signal {:title "Original signal"}))
  (gg/->image (gg/line domain denoised-signal {:title "Denoised signal"}))]]

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.transform)
