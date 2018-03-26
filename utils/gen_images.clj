(ns fastmath.utils.gen-images
  "Generate images from examples attached to metadata."
  (:require [clojure2d.core :refer :all]
            [clojure2d.color :refer :all]
            [clojure2d.pixels :as p]
            [fastmath.core :as m]
            [fastmath.complex :as c]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.interpolation :as i]))

;; core

(def ^:const bg-color 0x30426a)

(defn save-canvas
  ""
  ([c prefix n suff]
   (save c (str "docs/images/" prefix "/" n suff)))
  ([c prefix n]
   (save-canvas c prefix n ".jpg")))


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
            (set-background bg-color)
            (generate-math-graph f (- m/PI) m/PI))]
    (save-canvas c "m" n ".png")))

(doseq [s m/interp-list]
  (let [n (name s)
        f (eval `(fn [x#] (~s 0.0 1.0 x#)))
        c (with-canvas-> (canvas 200 200)
            (set-background bg-color)
            (generate-math-graph f -0.1 1.1))]
    (save-canvas c "m" n ".png")))

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
              (set-background bg-color)
              (generate-complex-graph f))]
      (save-canvas c "c" n))))

;; random/noise

(def background-vec (to-color bg-color))
(def white (to-color :white))

(defn icolor 
  [t]
  (v/interpolate background-vec white t))

;; discrete noise
(let [canvas (canvas 200 200)]
  (with-canvas [c canvas]
    (set-background c bg-color)
    (dotimes [x 180]
      (dotimes [y 180]
        (let [n (r/discrete-noise x y)]
          (set-color c (icolor n))
          (point c (+ x 10) (+ y 10))))))
  (save-canvas canvas "n" "discrete_noise"))

(defn draw-and-save-noise
  "Draw noise and save to the file"
  [f nm]
  (let [canvas (canvas 200 200)]
    (with-canvas [c canvas] 
      (set-background c bg-color)
      (dotimes [x 180]
        (dotimes [y 180]
          (let [xx (/ x 50.0)
                yy (/ y 50.0)
                n (f xx yy)]
            (set-color c (icolor n))
            (point c (+ x 10) (+ y 10))))))
    (save-canvas canvas "n" nm)))

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
      (set-background c bg-color)
      (set-color c :white 200)
      (set-stroke c 2.0)
      (doseq [^fastmath.vector.Vec2 res (take 1000 f)]
        (let [x (* (.x res) 200)
              y (* (.y res) 200)]
          (point c x y))))
    (save-canvas canvas "r" nm)))

(draw-and-save-random-2d ((r/make-sequence-generator :sobol 2)) "sobol")
(draw-and-save-random-2d ((r/make-sequence-generator :halton 2)) "halton")
(draw-and-save-random-2d ((r/make-sequence-generator :default 2)) "default")
(draw-and-save-random-2d (map #(v/add (v/mult % 0.1) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :gaussian 2))) "gaussian")
(draw-and-save-random-2d (map #(v/add (v/mult % 0.4) (v/vec2 0.5 0.5)) ((r/make-sequence-generator :sphere 2))) "sphere")

;;;; interpolations

(defn real-fn
  ""
  [x]
  (let [a (m/norm x 0 200 0 (* 2.2 m/PI))
        f (m/sin (* a (* 0.5 (m/cos (inc a)))))]
    (m/norm f -1.0 1.0 20 180)))

(def xs [20 50 58 66 100 121 140 150 160 170 180])
(def ys (map real-fn xs))

(def real-fn-path (for [x (range 200)]
                    (v/vec2 x (real-fn x))))

(defn draw-interpolation
  "Draw interpolation line"
  [interp-name interp]
  (let [canvas (canvas 200 200)
        p (remove nil? (for [x (range 200)]
                         (try
                           (v/vec2 x (interp x))
                           (catch Exception e nil))))]
    (with-canvas [c canvas]
      (set-background c bg-color)
      (set-color c :cyan)
      (doseq [[x y] (map vector xs ys)]
        (ellipse c x y 8 8 true))
      (set-color c 255 100 100 200)
      (path c real-fn-path)
      (set-color c :white 200)
      (set-stroke c 2.0)
      (path c p))
    (save-canvas canvas "i" (name interp-name))))

(draw-interpolation :shepard (i/shepard-interpolator xs ys))
(draw-interpolation :shepard-5 (i/shepard-interpolator xs ys 5))
(draw-interpolation :shepard-09 (i/shepard-interpolator xs ys 0.9))

(draw-interpolation :linear-smile (i/linear-smile-interpolator xs ys))
(draw-interpolation :kriging (i/kriging-spline-interpolator xs ys))
(draw-interpolation :cubic-spline (i/cubic-spline-interpolator xs ys))
(draw-interpolation :spline (i/spline-interpolator xs ys))
(draw-interpolation :neville (i/neville-interpolator xs ys))
(draw-interpolation :divided-diff (i/divided-difference-interpolator xs ys))
(draw-interpolation :akima (i/akima-spline-interpolator xs ys))
(draw-interpolation :linear (i/linear-interpolator xs ys))

(draw-interpolation :loess (i/loess-interpolator xs ys))
(draw-interpolation :loess-07-4 (i/loess-interpolator xs ys 0.7 4))
(draw-interpolation :loess-02-1 (i/loess-interpolator xs ys 0.2 1))

(draw-interpolation :rbf-gaussian (i/rbf-interpolator xs ys (i/rbf :gaussian 120)))
(draw-interpolation :rbf-imultiquadratic (i/rbf-interpolator xs ys (i/rbf :inverse-multiquadratic 80)))
(draw-interpolation :rbf-multiquadratic (i/rbf-interpolator xs ys (i/rbf :multiquadratic 120)))
(draw-interpolation :rbf-thinplate (i/rbf-interpolator xs ys (i/rbf :thinplate 80)))

(draw-interpolation :microsphere (i/microsphere-projection-interpolator xs ys 8 0.9 0.0000001 1 1.5 false 1))

;; draw 2d interpolation

(def xs [20 50 58 66 100 121 140 150 160 170 180])
(def ys [20 30 58 66 90  121 140 152 170 172 180])

(def vs (map (fn [_] (repeatedly (count xs) #(r/irand 256))) ys))

(defn draw-interpolation-2d
  "Draw interpolation map"
  [interp-name interp]
  (let [canvas (canvas 200 200)]
    (with-canvas [c canvas]
      (set-background c bg-color)
      (dotimes [x 200]
        (dotimes [y 200]
          (try
            (set-color c :white (interp x y))
            (rect c x y 1 1)
            (catch Exception e nil))))
      (set-color c 255 100 100 200)
      (doseq [x xs
              y ys]
        (ellipse c x y 5 5 true)))
    (save-canvas canvas "i" (name interp-name))))

(draw-interpolation-2d :bicubic (i/bicubic-interpolator xs ys vs))
(draw-interpolation-2d :piecewise-bicubic (i/piecewise-bicubic-interpolator xs ys vs))
(draw-interpolation-2d :microsphere-2d (i/microsphere-2d-projection-interpolator xs ys vs 32 0.99 0.00001 1 2 true 1))
(draw-interpolation-2d :bilinear (i/bilinear-interpolator xs ys vs))
(draw-interpolation-2d :bicubic-smile (i/bicubic-smile-interpolator xs ys vs))
(draw-interpolation-2d :cubic-2d (i/cubic-2d-interpolator xs ys vs))

;;

(defn draw-distribution
  "Draw real distribution."
  [a b s distr-name & distr]
  (let [canvas (canvas 200 200)
        p (for [d distr]
            (for [x (range 1 200)]
              (let [xx (m/norm x 0 200 a b)
                    y (- 200 (* s (r/pdf d xx)))]
                (v/vec2 x y))))]
    (with-canvas [c canvas]
      (set-background c bg-color)
      (doseq [[p col] (map vector p [:white
                                     (color 102 190 141)
                                     (color 258 164 88)
                                     (color 244 77 81)])] 
        (set-color c col)
        (path c p)))
    (save-canvas canvas "d" (name distr-name))))

(draw-distribution 0 3 200 :levy
                   (r/real-distribution :levy {:c 0.5})
                   (r/real-distribution :levy {:c 1})
                   (r/real-distribution :levy {:c 2}))

(draw-distribution 0.0001 1 50 :beta
                   (r/real-distribution :beta {:alpha 2 :beta 5})
                   (r/real-distribution :beta {:alpha 0.5 :beta 0.5})
                   (r/real-distribution :beta {:alpha 2.5 :beta 2}))

(draw-distribution -4 4 250 :cauchy
                   (r/real-distribution :cauchy {})
                   (r/real-distribution :cauchy {:mean 0.5 :scale 3.0})
                   (r/real-distribution :cauchy {:scale 0.5}))

(draw-distribution 0 8 250 :chi-squared
                   (r/real-distribution :chi-squared {:degrees-of-freedom 1})
                   (r/real-distribution :chi-squared {:degrees-of-freedom 5})
                   (r/real-distribution :chi-squared {:degrees-of-freedom 9}))

(draw-distribution -1.1 1.1 60 :empirical
                   (r/real-distribution :empirical {:bin-count 10
                                                    :data (repeatedly 10000 #(m/sin (r/drand 0 m/TWO_PI)))})
                   (r/real-distribution :empirical {:bin-count 15
                                                    :data (repeatedly 10000 #(* (r/drand -1 1) (r/drand -1 1) (r/drand -1 1)))}))

(draw-distribution 0.0001 2.5 50 :weibull
                   (r/real-distribution :weibull {:alpha 0.5 :beta 1.0})
                   (r/real-distribution :weibull {:alpha 1.5 :beta 1.0})
                   (r/real-distribution :weibull {:alpha 5.0 :beta 1.0}))
