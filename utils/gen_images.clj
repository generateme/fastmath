(ns fastmath.utils.gen-images
  "Generate images from examples attached to metadata."
  (:require [clojure2d.core :refer :all]
            [clojure2d.color :refer :all]
            [clojure2d.pixels :as p]
            [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]
            [fastmath.vector :as v]))

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
        c (with-canvas-> (canvas 200 200)
            (set-background 0x30426a)
            (generate-math-graph f (- m/PI) m/PI))]
    (save c (str "docs/images/m/" n ".png"))))

(doseq [s m/interp-list]
  (let [n (name s)
        f (eval `(fn [x#] (~s 0.0 1.0 x#)))
        c (with-canvas-> (canvas 200 200)
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
          c (with-canvas-> (canvas 200 200)
              (set-background 0x30426a)
              (generate-complex-graph f))]
      (save c (str "docs/images/c/" n ".jpg")))))

;; random/noise

(def background-vec (to-color 0x30426a))
(def white (to-color :white))

(defn icolor 
  [t]
  (v/interpolate background-vec white t))

;; discrete noise
(let [canvas (canvas 200 200)]
  (with-canvas [c canvas]
    (set-background c 0x30426a)
    (dotimes [x 180]
      (dotimes [y 180]
        (let [n (r/discrete-noise x y)]
          (set-color c (icolor n))
          (point c (+ x 10) (+ y 10))))))
  (save canvas (str "docs/images/n/discrete_noise.jpg")))

(defn draw-and-save-noise
  "Draw noise and save to the file"
  [f nm]
  (let [canvas (canvas 200 200)]
    (with-canvas [c canvas] 
      (set-background c 0x30426a)
      (dotimes [x 180]
        (dotimes [y 180]
          (let [xx (/ x 50.0)
                yy (/ y 50.0)
                n (f xx yy)]
            (set-color c (icolor n))
            (point c (+ x 10) (+ y 10))))))
    (save canvas (str "docs/images/n/" nm ".jpg"))))

;; noise

(doseq [s `(r/noise r/vnoise r/simplex)]
  (let [nm (name s)
        f (symbol->fn2 s)]
    (draw-and-save-noise f nm)))

;; combined noise (fbm...)

(draw-and-save-noise (r/make-fbm-noise {:interpolation :linear
                                        :noise-type :value}) "fbm")

(draw-and-save-noise (r/make-ridgedmulti-noise {:octaves 3
                                                :lacunarity 2.1
                                                :gain 0.7
                                                :noise-type :simplex}) "ridgedmulti")

(draw-and-save-noise (r/make-billow-noise {:seed 12345
                                           :interpolation :linear}) "billow")

(draw-and-save-noise (r/make-single-noise {:interpolation :linear}) "single")

(draw-and-save-noise (r/make-random-noise-fn) "random1")
(draw-and-save-noise (r/make-random-noise-fn) "random2")
(draw-and-save-noise (r/make-random-noise-fn) "random3")


;; random distributions
;; random series

(defn draw-and-save-random-2d
  "Draw noise and save to the file"
  [f nm]
  (let [canvas (canvas 200 200)]
    (with-canvas [c canvas]
      (set-background c 0x30426a)
      (set-color c :white 200)
      (set-stroke c 2.0)
      (doseq [^fastmath.vector.Vec2 res (take 1000 f)]
        (let [x (* (.x res) 200)
              y (* (.y res) 200)]
          (point c x y))))
    (save canvas (str "docs/images/r/" nm ".jpg"))))

(draw-and-save-random-2d ((r/make-sequence-generator :sobol 2)) "sobol")
(draw-and-save-random-2d ((r/make-sequence-generator :halton 2)) "halton")
(draw-and-save-random-2d ((r/make-sequence-generator :default 2)) "default")
(draw-and-save-random-2d (map #(v/add (v/mult % 0.1) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :gaussian 2))) "gaussian")
(draw-and-save-random-2d (map #(v/add (v/mult % 0.4) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :sphere 2))) "sphere")

