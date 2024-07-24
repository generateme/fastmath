^:kindly/hide-code
(ns index
  (:require [scicloj.clay.v2.api :as clay]))

;; # Preface {.unnumbered}

^:kindly/hide-code
(comment
  (clay/make! {:source-path ["index.clj"
                             "core.clj"
                             "vector_matrix.clj"
                             "random.clj"
                             "stats.clj"
                             "complex_quaternions.clj"
                             "special.clj"]
               :book {:favicon "clay/resources/favicon.png"
                      :title "Fastmath documentation"}}))


