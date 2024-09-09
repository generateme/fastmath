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


(defmacro symbol-info-table
  [rows]
  `(kind/table {:column-names ["Symbol" "Info"]
                :row-vectors [~@(for [[s i] rows]
                                  `[(kind/code (pr-str (quote ~s))) (or ~i "-")])]}))


;; build whole book
(comment
  (clay/make! {:source-path ["index.clj"
                             "core.clj"
                             "vector_matrix.clj"
                             "random.clj"
                             "stats.clj"
                             "calculus.clj"
                             "complex_quaternions.clj"
                             "special.clj"]
               :format [:quarto :html]
               :book {:favicon "clay/resources/favicon.png"
                      :title "Fastmath documentation"}}))
