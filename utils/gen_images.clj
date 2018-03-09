(ns fastmath.utils.gen-images
  "Generate images from examples attached to metadata."
  (:require [clojure2d.core :refer :all]
            [clojure2d.color :refer :all]
            [clojure2d.pixels :as p]
            [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]))

;; core

(defn symbol->fn
  "Convert symbol to function"
  [s] (eval `(fn [x#] (~s x#))))

(defn symbol->fn2
  "Convert symbol to 2 parameter function"
  [s] (eval `(fn [x# y#] (~s x# y#))))


(defn generate-math-graph
  "Generate math graphs for function"
  [canvas f a b]
  (let [mid-x (m/norm 0 a b 0 200)
        mid-y (- 200 mid-x)]
    (-> canvas
        (set-color 60 100 120)
        (set-stroke 0.8)
        (line 0 mid-y 200 mid-y)
        (line mid-x 0 mid-x 200)
        (set-color :white 130) 
        (set-stroke 1.5))
    (dotimes [x 400]
      (let [xx (m/norm x 0 400 a b) 
            y (m/norm (f xx) b a 0 200)]
        (point canvas (* 0.5 x) y))))
  canvas)

(doseq [s m/single-list]
  (let [n (name s)
        f (symbol->fn s)
        c (with-canvas-> (make-canvas 200 200)
            (set-background 0x30426a)
            (generate-math-graph f (- m/PI) m/PI))]
    (save c (str "docs/images/m/" n ".png"))))

(doseq [s m/interp-list]
  (let [n (name s)
        f (eval `(fn [x#] (~s 0.0 1.0 x#)))
        c (with-canvas-> (make-canvas 200 200)
            (set-background 0x30426a)
            (generate-math-graph f -0.1 1.1))]
    (save c (str "docs/images/m/" n ".png"))))

;; complex

(defn generate-complex-graph
  "Generate graph for complex fn."
  [canvas f]
  (let [w (width canvas)
        h (height canvas)]
    (set-stroke canvas 1.5)
    (set-color canvas :white 30)
    (dotimes [x w]
      (dotimes [y h]
        (let [xx (m/norm x 0 w (- m/PI) m/PI)
              yy (m/norm y 0 h (- m/PI) m/PI)
              res (f (c/complex xx yy))
              resx (m/norm (res 0) (- m/PI) m/PI 0 w)
              resy (m/norm (res 1) (- m/PI) m/PI 0 h)]
          (point canvas resx resy)))))
  canvas)

(binding [*jpeg-image-quality* 0.85]
  (doseq [s c/fn-list]
    (let [n (name s)
          f (symbol->fn s)
          c (with-canvas-> (make-canvas 200 200)
              (set-background 0x30426a)
              (generate-complex-graph f))]
      (save c (str "docs/images/c/" n ".jpg")))))

;; random/noise

;; discrete noise
(let [canvas (make-canvas 200 200)]
  (with-canvas [c canvas]
    (set-background c 0x30426a)
    (dotimes [x 180]
      (dotimes [y 180]
        (let [n (r/discrete-noise x y)
              col (* n 255.0)]
          (set-color c col col col)
          (point c (+ x 10) (+ y 10))))))
  (save canvas (str "docs/images/n/discrete_noise.jpg")))

(defn draw-and-save-noise
  "Draw noise and save to the file"
  [f nm]
  (let [canvas (make-canvas 200 200)]
    (with-canvas [c canvas]
      (set-background c 0x30426a)
      (dotimes [x 180]
        (dotimes [y 180]
          (let [xx (/ x 50.0)
                yy (/ y 50.0)
                n (f xx yy)
                col (* n 255.0)]
            (set-color c col col col)
            (point c (+ x 10) (+ y 10))))))
    (save canvas (str "docs/images/n/" nm ".jpg"))))

;; noise

(doseq [s r/noise-fn-list]
  (let [nm (name s)
        f (symbol->fn2 s)]
    (draw-and-save-noise f nm)))

;; combined noise (fbm...)
