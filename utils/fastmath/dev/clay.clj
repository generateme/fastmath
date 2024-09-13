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


(defmacro dgraph-table
  "Helper macro to render `rows` to a table containing PDF, CDF, iCDF graphs images of distributions.
  
  The `rows` has even length. Starting from the first row, a pair of adjacent rows perform a group.
  The first row in a group has columns containing the parameter of a distribution.
  The next row in the same group has columns containing a vector of PDF, CDF, iCDF charts of a distribution
  whose parameters are in the same column in the previous row.

  Generates a table like this:
  | Param-1      |
  | Charts-1     |
  | Param-2      |
  | Charts-2     |
  "
  [rows]
  `(kind/table ~(->> (map (fn [p g]
                            [(kind/hiccup [:dl [:dt nil `(kind/code (str ~p))] [:dd nil `(kind/table [~g])]])])
                          (mapcat identity (take-nth 2 rows)) (mapcat identity (take-nth 2 (rest rows))))
                     (into []))))


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
