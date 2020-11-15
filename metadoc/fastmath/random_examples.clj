(ns fastmath.random-examples
  (:require [fastmath.random :refer :all]
            [metadoc.examples :refer :all]))

(add-examples rngs-list
  (example "Contains" (sort rngs-list)))

(add-examples rng
  (example-session "Creation" (rng :mersenne) (rng :isaac 1234))
  (example "Usage" (irandom (rng :mersenne 999) 15 25)))

(defsnippet rngproto-snippet
  "Show [[RNGProto]] methods."
  (let [rng (rng :well44497b)]
    (f rng)))

(add-examples irandom (example-snippet "integer" rngproto-snippet irandom))
(add-examples lrandom (example-snippet "long" rngproto-snippet lrandom))
(add-examples drandom (example-snippet "double" rngproto-snippet drandom))
(add-examples frandom (example-snippet "float" rngproto-snippet frandom))
(add-examples grandom (example-snippet "gaussian double" rngproto-snippet grandom))
(add-examples brandom (example-snippet "boolean" rngproto-snippet brandom))

(add-examples set-seed!
              (example "Set seed for the RNG object" {:test-value 10} (let [rng (rng :isaac)]
                                                                        (set-seed! rng 1234)
                                                                        (irandom rng 10 15)))
              (example "Set seed for the distribution object" {:test-value 2} (let [d (distribution :enumerated-int {:data [1 1 1 2 3]})] 
                                                                                (set-seed! d 1234)
                                                                                (irandom d)))
              (example "Set seed to SMILE (global) RNG" (set-seed! :smile 1234))
              (example "Set seed to both SMILE and fastmath default RNGs" (set-seed! 1234))
              (example "Set to random seed both default RNGs" (set-seed!)))


(add-examples default-rng
  (example-session "Usage"
    (set-seed! default-rng 111)
    (irandom default-rng)
    (set-seed! default-rng 999)
    (irandom default-rng)
    (set-seed! default-rng 111)
    (irandom default-rng)))

(add-examples frand (example-session "Usage" (frand) (frand 10) (frand 10 20)))

(add-examples brand
  (example-session "Usage" (brand) (brand 0.1))
  (example "Count number of `true` values with probability 0.15" (count (filter true? (repeatedly 100000 #(brand 0.15))))))

(add-examples drand (example-session "Usage" (drand) (drand 10) (drand 10 20)))
(add-examples grand (example-session "Usage" (grand) (grand 10) (grand 10 20)))
(add-examples irand (example-session "Usage" (irand) (irand 10) (irand 10 20)))
(add-examples lrand (example-session "Usage" (lrand) (lrand 10) (lrand 10 20)))

(add-examples randval
  (example-session "Usage" (randval :val-one :val-two) (randval 0.001 :low-probability :high-probability))
  (example "Check probability of nil (should return value around 1000)." (count (filter nil? (repeatedly 1000000 #(randval 0.001 nil 101))))))

(add-examples flip
  (example-session "Usage" (flip) (flip 0.2) (repeatedly 10 #(flip 0.1))))

(add-examples flipb
  (example-session "Usage" (flipb) (flipb 0.2) (repeatedly 10 #(flipb 0.1))))

;; generators

(add-examples sequence-generators-list
  (example "Generator names." (sort sequence-generators-list)))

(add-examples sequence-generator
  (example "Usage (2d)" (let [gen (sequence-generator :halton 2)]
                          (take 5 gen)))
  (example "Usage (1d)" (let [gen (sequence-generator :sobol 1)]
                          (take 5 gen)))
  (example "Usage (10d)" (second (sequence-generator :halton 10)))
  (example "Usage, R2 sequence" (take 5 (sequence-generator :r2 3)))
  (example-image "R2 plot (500 samples)" "images/r/r2.jpg")
  (example-image "Halton plot (500 samples)" "images/r/halton.jpg")
  (example-image "Sobol plot (500 samples)" "images/r/sobol.jpg")
  (example-image "Sphere plot (500 samples)" "images/r/sphere.jpg")
  (example-image "Gaussian plot (500 samples)" "images/r/gaussian.jpg")
  (example-image "Default plot (500 samples)" "images/r/default.jpg"))

(add-examples jittered-sequence-generator
  (example (let [gen1 (jittered-sequence-generator :r2 2 0.5)
                 gen2 (jittered-sequence-generator :r2 2 0.5)]
             [(first gen1) (first gen2)]))
  (example-image "Jittered (0.5) R2 plot (500 samples)" "images/r/jr2.jpg")
  (example-image "Jittered (0.5) Halton plot (500 samples)" "images/r/jhalton.jpg")
  (example-image "Jittered (0.5) Sobol plot (500 samples)" "images/r/jsobol.jpg")
  (example-image "Jittered (0.5) Sphere plot (500 samples)" "images/r/jsphere.jpg")
  (example-image "Jittered (0.5) Gaussian plot (500 samples)" "images/r/jgaussian.jpg")
  (example-image "Jittered (0.5) Default plot (500 samples)" "images/r/jdefault.jpg"))

;; noise

(add-examples noise-interpolations (example "List of names (keys)" (keys noise-interpolations)))
(add-examples noise-generators (example "List of names (keys)" (keys noise-generators)))
(add-examples noise-types (example "List of names (keys)" (keys noise-types)))

(add-examples vnoise
  (example-session "Usage"
    (vnoise 3.3)
    (vnoise 3.3 1.1)
    (vnoise 3.3 0.0 -0.1))
  (example-image "2d noise" "images/n/vnoise.jpg"))

(add-examples noise
  (example-session "Usage"
    (noise 3.3)
    (noise 3.3 1.1)
    (noise 3.3 0.0 -0.1))
  (example-image "2d noise" "images/n/noise.jpg"))

(add-examples simplex
  (example-session "Usage"
    (simplex 3.3)
    (simplex 3.3 1.1)
    (simplex 3.3 0.0 -0.1))
  (example-image "2d noise" "images/n/simplex.jpg"))

(add-examples random-noise-cfg
  (example "Random configuration" (random-noise-cfg)))

(add-examples random-noise-fn
  (example-session "Create function"
    (random-noise-fn)
    (random-noise-fn (random-noise-cfg)))
  (example-image "One" "images/n/random1.jpg")
  (example-image "Two" "images/n/random2.jpg")
  (example-image "Three" "images/n/random3.jpg"))

(add-examples single-noise
  (example "Usage"
    (let [n (single-noise {:interpolation :linear})]
      (n 0.5 1.1 -1.3)))
  (example-image "2d noise" "images/n/single.jpg"))

(add-examples fbm-noise
  (example "Usage"
    (let [n (fbm-noise {:interpolation :linear
                        :noise-type :value})]
      (n 0.5 1.1 -1.3)))
  (example-image "2d noise" "images/n/fbm.jpg"))

(add-examples billow-noise
  (example "Usage"
    (let [n (billow-noise {:seed 12345
                           :interpolation :none})]
      (n 0.5 1.1 -1.3)))
  (example-image "2d noise" "images/n/billow.jpg"))

(add-examples ridgedmulti-noise
  (example "Usage"
    (let [n (ridgedmulti-noise {:octaves 3
                                :lacunarity 2.1
                                :gain 0.7
                                :noise-type :simplex})]
      (n 0.5 1.1 -1.3)))
  (example-image "2d noise" "images/n/ridgedmulti.jpg"))

(add-examples discrete-noise
  (example-session "Example calls"
    (discrete-noise 123 444)
    (discrete-noise 123 444)
    (discrete-noise 123 445)
    (discrete-noise 123))
  (example-image "Draw noise for [0-180] range." "images/n/discrete_noise.jpg"))

(add-examples warp-noise-fn
  (example (let [n (warp-noise-fn simplex 2.0 2.0)]
             [(n 0.0)
              (n 1.0 0.5)
              (n 2 2 2)]))
  (example-image "Default warp (noise=vnoise, scale=4.0, depth=1.0)." "images/n/warp.jpg"))

(add-examples distribution
  (example-session "Usage"
    (distribution :beta)
    (distribution :beta {:alpha 1.0 :beta 1.0}))
  (example "All parameters"
    (into (sorted-map) (map #(vector % (sort (distribution-parameters (distribution %)))) (keys (methods distribution))))))

(add-examples distributions-list
  (example-session "Number and list of distributions" distributions-list (count distributions-list)))

(doseq [n (disj distributions-list
                :categorical-distribution :cramer-von-mises :logarithmic :normal-inverse-gaussian
                :multi-normal :dirichlet :kde)]
  (add-examples distribution (example-image (str "PDFs of " (name n)) (str "images/d/pdf-" (name n) ".jpg")))
  (add-examples distribution (example-image (str "CDFs of " (name n)) (str "images/d/cdf-" (name n) ".jpg")))
  (add-examples distribution (example-image (str "ICDFs of " (name n)) (str "images/d/icdf-" (name n) ".jpg"))))

(add-examples distribution (example-image "2d multidimensional normal (mean=[0,0], covariances=I)"
                             (str "images/d/multi-normal.jpg")))
(add-examples distribution (example-image "2d multidimensional normal (mean=[0,0], covariances=[[1 -1] [-1 2]])"
                             (str "images/d/multi-normal2.jpg")))
(add-examples distribution (example-image "2d dirichlet (alpha=[2,0.8])"
                             (str "images/d/dirichlet.jpg")))

;;

(add-examples cdf (example-session "Usage" (cdf (distribution :gamma) 1) (cdf (distribution :gamma) 1 4)))
(add-examples pdf (example-session "Usage"
                    (pdf (distribution :gamma) 1)
                    (pdf (distribution :pascal) 1)))
(add-examples probability (example-session "Usage"
                            (probability (distribution :gamma) 1)
                            (probability (distribution :pascal) 1)))
(add-examples lpdf (example "Usage" (lpdf (distribution :gamma) 1)))
(add-examples icdf (example "Usage" (icdf (distribution :gamma) 0.5)))
(add-examples mean (example "Usage" {:test-value 4.0} (mean (distribution :gamma))))
(add-examples variance (example "Usage" {:test-value 8.0} (variance (distribution :gamma))))
(add-examples lower-bound (example "Usage"  {:test-value 0.0} (lower-bound (distribution :gamma))))
(add-examples upper-bound (example "Usage" {:test-value ##Inf} (upper-bound (distribution :gamma))))
(add-examples sample (example "Random value from distribution" (sample (distribution :gamma))))
(add-examples ->seq (example "Sequence of random values from distribution" (->seq (distribution :gamma) 5)))
(add-examples log-likelihood (example "Usage" (log-likelihood (distribution :gamma) [10 0.5 0.5 1 2])))
(add-examples likelihood (example "Usage" (likelihood (distribution :gamma) [10 0.5 0.5 1 2])))

(add-examples drandom (example "Double random value from distribution" (drandom (distribution :gamma))))
(add-examples irandom (example "Integer random value from distribution (sample cast to `int`)" (irandom (distribution :gamma))))
(add-examples frandom (example "Float random value from distribution (sample cast to `float`)" (frandom (distribution :gamma))))
(add-examples lrandom (example "Long random value from distribution (sample cast to `long`)" (lrandom (distribution :gamma))))

(add-examples continuous?
  (example-session "Usage" {:test-values [true false]}
    (continuous? (distribution :gamma))
    (continuous? (distribution :pascal))))

(add-examples covariance
  (example-session "Usage" {:test-values [[[1.0 0.0] [0.0 1.0]] [[0.05 -0.05] [-0.05 0.05]]]}
    (covariance (distribution :multi-normal))
    (covariance (distribution :dirichlet {:alpha [2 2]}))))

(add-examples means
  (example-session "Usage" {:test-values [[0.0 0.0] [0.5 0.5]]}
    (means (distribution :multi-normal))
    (means (distribution :dirichlet {:alpha [2 2]}))))

(add-examples default-normal
  (example-session "Usage"
    (sample default-normal)
    (set-seed! default-normal 1234)
    (sample default-normal)
    (irandom default-normal)
    (mean default-normal)
    (variance default-normal)))

(add-examples dimensions
  (example-session "Usage" {:test-values [1 30]}
    (dimensions (distribution :gamma))
    (dimensions (distribution :dirichlet {:alpha (repeat 30 2.0)}))))

(add-examples distribution-id
  (example-session "Usage" {:test-values [:gamma :normal]}
    (distribution-id (distribution :gamma))
    (distribution-id default-normal)))

(add-examples distribution-parameters
  (example-session "Usage" {:test-values [[:scale :shape] [:rng :shape :scale :inverse-cumm-accuracy] [:sd :mu]]}
    (distribution-parameters (distribution :gamma))
    (distribution-parameters (distribution :gamma) true)
    (distribution-parameters default-normal)))

(add-examples observe (example "Usage" (observe (distribution :gamma) [10 0.5 0.5 1 2])))
(add-examples observe1 (example "Usage" (observe1 (distribution :gamma) 10)))

(add-examples source-object (example (source-object default-normal)))

(add-examples synced-rng (example (drandom (synced-rng :mersenne 1234))))

