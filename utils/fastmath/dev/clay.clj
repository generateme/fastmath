(ns fastmath.dev.clay
  (:require [clojure.string :as str]
            [scicloj.clay.v2.api :as clay]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.kindly.v4.api :as kindly]))

(defmacro callout
  [what title & forms]
  `(kind/fragment
    [(kind/md (str "::: {.callout-" ~what " title=\"" ~title "\"}"))
     ~@forms
     (kind/md (str ":::"))]))

(defmacro examples [& forms]
  `(->> (map (fn [form# value#] (str form# " ;; => " (pr-str value#)))
             ~(mapv (fn [form#] `'~form#) forms)
             ~(vec forms))
        (str/join "\n")
        kind/code
        kindly/hide-code))

(defmacro examples-note
  [& forms]
  `(callout "note" "Examples" (examples ~@forms)))

;; build whole book
(comment
  (clay/make! {:source-path ["index.clj"
                             "core.clj"
                             "vector_matrix.clj"
                             "random.clj"
                             "stats.clj"
                             "polynomials.clj"
                             "special.clj"
                             "calculus.clj"
                             "interpolation.clj"
                             "optimization.clj"
                             "transform.clj"
                             "ml.clj"
                             "complex_quaternions.clj"
                             "easings.clj"
                             "fields.clj"
                             "efloat.clj"]
               :format [:quarto :html]
               :book {:favicon "clay/resources/favicon.png"
                      :title "Fastmath documentation"}}))
