(ns fastmath.easings-examples
  (:require [fastmath.easings :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.core :as m]
            [incanter.charts :as c]
            [incanter.core :as i]))

(defsnippet fastmath.easings save-graph
  "Save incanter graph"
  (let [fname (str "images/e/" (first opts) ".png")]
    (i/save (c/function-plot f 0.0 1.0 :y-label "easing value") (str "docs/" fname) :width 500 :height 250)
    fname))

(add-examples linear
  (example "Usage" {:test-value 0.5} (linear 0.5))
  (example-snippet "Plot" save-graph :image linear))

(add-examples back-in
  (example "Usage" {:test-value -0.087698} (m/approx (back-in 0.5) 6))
  (example "Different overshot" (back-in 2.0 0.5))
  (example-snippet "Plot" save-graph :image back-in)
  (example-snippet "Plot, different `s` value." save-graph :image (partial back-in 5.0)))

(add-examples back-out
  (example "Usage" {:test-value (- 1.0 -0.087698)} (m/approx (back-out 0.5) 6))
  (example "Different overshot" (back-out 2.0 0.5))
  (example-snippet "Plot" save-graph :image back-out)
  (example-snippet "Plot, different `s` value." save-graph :image (partial back-out 5.0)))

(add-examples back-in-out
  (example "Usage" {:test-value -0.043849} (m/approx (back-in-out 0.25) 6))
  (example "Different overshot" (back-in-out 2.0 0.75))
  (example-snippet "Plot" save-graph :image back-in-out)
  (example-snippet "Plot, different `s` value." save-graph :image (partial back-in-out 5.0)))

(add-examples bounce-out
  (example "Usage" {:test-value 0.3025} (m/approx (bounce-out 0.2) 6))
  (example-snippet "Plot" save-graph :image bounce-out))

(add-examples bounce-in
  (example "Usage" {:test-value 0.06} (m/approx (bounce-in 0.2) 6))
  (example-snippet "Plot" save-graph :image bounce-in))

(add-examples bounce-in-out
  (example "Usage" (m/approx (bounce-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image bounce-in-out))

(add-examples circle-in
  (example "Usage" {:test-value 0.020204} (m/approx (circle-in 0.2) 6))
  (example-snippet "Plot" save-graph :image circle-in))

(add-examples circle-out
  (example "Usage" {:test-value 0.6} (m/approx (circle-out 0.2) 6))
  (example-snippet "Plot" save-graph :image circle-out))

(add-examples circle-in-out
  (example "Usage" (m/approx (circle-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image circle-in-out))

(add-examples cubic-in
  (example "Usage" {:test-value 0.008} (m/approx (cubic-in 0.2) 6))
  (example-snippet "Plot" save-graph :image cubic-in))

(add-examples cubic-out
  (example "Usage" {:test-value 0.488} (m/approx (cubic-out 0.2) 6))
  (example-snippet "Plot" save-graph :image cubic-out))

(add-examples cubic-in-out
  (example "Usage" (m/approx (cubic-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image cubic-in-out))

(add-examples elastic-in
  (example "Usage" {:test-value -0.001953} (m/approx (elastic-in 0.2) 6))
  (example "Create custom elastic easing" (m/approx ((elastic-in 1.0 0.1) 0.2) 6))
  (example-snippet "Plot" save-graph :image elastic-in)
  (example-snippet "Plot, different settings" save-graph :image (elastic-in 1.1 0.1)))

(add-examples elastic-out
  (example "Usage" {:test-value 1.125} (m/approx (elastic-out 0.2) 6))
  (example "Create custom elastic easing" (m/approx ((elastic-out 1.0 0.1) 0.2) 6))
  (example-snippet "Plot" save-graph :image elastic-out)
  (example-snippet "Plot, different settings" save-graph :image (elastic-out 1.1 0.1)))

(add-examples elastic-in-out
  (example "Usage" (m/approx (elastic-in-out 0.2) 6))
  (example "Create custom elastic easing" (m/approx ((elastic-in-out 1.0 0.1) 0.2) 6))
  (example-snippet "Plot" save-graph :image elastic-in-out)
  (example-snippet "Plot, different settings" save-graph :image (elastic-in-out 1.1 0.1)))

(add-examples exp-in
  (example "Usage" {:test-value 0.003906} (m/approx (exp-in 0.2) 6))
  (example-snippet "Plot" save-graph :image exp-in))

(add-examples exp-out
  (example "Usage" {:test-value 0.75} (m/approx (exp-out 0.2) 6))
  (example-snippet "Plot" save-graph :image exp-out))

(add-examples exp-in-out
  (example "Usage" (m/approx (exp-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image exp-in-out))

(add-examples poly-in
  (example "Usage" {:test-value 0.008} (m/approx (poly-in 0.2) 6))
  (example "Other exponent" (poly-in 8.0 0.5))
  (example-snippet "Plot" save-graph :image poly-in)
  (example-snippet "Plot" save-graph :image (partial poly-in 8.0)))

(add-examples poly-out
  (example "Usage" {:test-value 0.488} (m/approx (poly-out 0.2) 6))
  (example "Other exponent" (poly-out 8.0 0.5))
  (example-snippet "Plot" save-graph :image poly-out)
  (example-snippet "Plot" save-graph :image (partial poly-out 8.0)))

(add-examples poly-in-out
  (example "Usage" (m/approx (poly-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image poly-in-out)
  (example-snippet "Plot" save-graph :image (partial poly-in-out 8.0)))

(add-examples quad-in
  (example "Usage" {:test-value 0.04} (m/approx (quad-in 0.2) 6))
  (example-snippet "Plot" save-graph :image quad-in))

(add-examples quad-out
  (example "Usage" {:test-value 0.36} (m/approx (quad-out 0.2) 6))
  (example-snippet "Plot" save-graph :image quad-out))

(add-examples quad-in-out
  (example "Usage" (m/approx (quad-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image quad-in-out))

(add-examples sin-in
  (example "Usage" {:test-value 0.048943} (m/approx (sin-in 0.2) 6))
  (example-snippet "Plot" save-graph :image sin-in))

(add-examples sin-out
  (example "Usage" {:test-value 0.309017} (m/approx (sin-out 0.2) 6))
  (example-snippet "Plot" save-graph :image sin-out))

(add-examples sin-in-out
  (example "Usage" (m/approx (sin-in-out 0.2) 6))
  (example-snippet "Plot" save-graph :image sin-in-out))

(add-examples out
  (example "Usage" (let [outeasing (out sin-in)]
                     (== ^double (sin-out 0.75) ^double (outeasing 0.75))))
  (example-snippet "Create out easing" save-graph :image (out sin-in-out)))

(add-examples in-out
  (example "Usage" (let [inouteasing (in-out quad-in)]
                     (== ^double (quad-in-out 0.75) ^double (inouteasing 0.75))))
  (example-snippet "Create in-out easing" save-graph :image (in-out sin-in-out)))

(add-examples reflect
  (example "Usage" (let [neasing (reflect (partial back-in 2.0) 0.2)]
                     [(neasing 0.1) (neasing 0.5) (neasing 0.9)]))
  (example "For `center=0.5` function returns regular in-out easing"
    {:test-value true} (== (back-in-out 0.4) ^double ((reflect back-in 0.5) 0.4)))
  (example-snippet "Reflect one easing in center=0.2 " save-graph :image (reflect elastic-in-out 0.2)))

(add-examples easings-list
  (example "List of easings" (sort (keys easings-list)))
  (example "Access to function" ((easings-list :back-in) 3.0 0.25)))
