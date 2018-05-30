(ns fastmath.rbf-examples
  (:require [fastmath.rbf :refer :all]
            [metadoc.examples :refer :all]
            [incanter.charts :as c]))

(defsnippet fastmath.rbf save-graph
  "Save graph"
  (let [fname (str "images/rbf/" (first opts) ".png")
        [x1 x2] params]
    (incanter.core/save (c/function-plot f x1 x2 :y-label "value") (str "docs/" fname) :width 250 :height 250)
    fname))

(add-examples rbf
  (example-snippet "Linear" save-graph :image (rbf :linear) -2.0 2.0)
  (example-snippet "Gaussian" save-graph :image (rbf :gaussian) -2.0 2.0)
  (example-snippet "Multiquadratic" save-graph :image (rbf :multiquadratic) -3.0 3.0)
  (example-snippet "Inverse multiquadratic" save-graph :image (rbf :inverse-multiquadratic) -3.0 3.0)
  (example-snippet "Inverse quadratic" save-graph :image (rbf :inverse-quadratic) -3.0 3.0)
  (example-snippet "Thinplate" save-graph :image (rbf :thinplate) -1.0 3.0)
  (example-snippet "Polyharmonic 3" save-graph :image (rbf :polyharmonic 3) -3.0 3.0)
  (example-snippet "Polyharmonic 4" save-graph :image (rbf :polyharmonic 4) -3.0 3.0)
  (example-snippet "Wendland" save-graph :image (rbf :wendland) -3.0 3.0)
  (example-snippet "Wu" save-graph :image (rbf :wu) -3.0 3.0))

(add-examples rbfs-list
  (example "List of names" (sort rbfs-list)))

(comment incanter.core/view (c/function-plot (rbf :polyharmonic 2) -3.0 3.0 :y-label "value") :width 250 :height 250)
