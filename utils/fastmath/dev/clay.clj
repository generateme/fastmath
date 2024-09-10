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

(defmacro dgraph-table-v1
  "Helper macro to render rows to a table containing PDF, CDF, iCDF graphs images of distributions.
  
  The rows has even length. Starting from the first row, a pair of adjacent rows perform a group.
  The first row in a group has columns containing the parameter of a distribution.
  The next row in the same group has columns containing PDF, CDF, iCDF graphs images of a distribution
  whose parameters are in the same column in the previous row. 

  the images are displayed side-by-side in the same column."
  [rows]
  `(kind/table ~(->> (map-indexed (fn [i r]
                                    (if (even? i)
                                      r
                                      (into []
                                            (map (fn [v]
                                                   `(kind/hiccup (into [:div nil] ~v))))
                                            r)))
                                  rows)
                     (into []))))


(defmacro dgraph-table-v2
  "Helper macro to render rows to a table containing PDF, CDF, iCDF graphs images of distributions.
  
  The rows has even length. Starting from the first row, a pair of adjacent rows perform a group.
  The first row in a group has columns containing the parameter of a distribution.
  The next row in the same group has columns containing PDF, CDF, iCDF graphs images of a distribution
  whose parameters are in the same column in the previous row.
  
  The images are put into columns inside a table."
  [rows]
  `(kind/table ~(->> (map-indexed (fn [i r]
                                    (if (even? i)
                                      r
                                      (into []
                                            (map (fn [v]
                                                   `(kind/table [~v])))
                                            r)))
                                  rows)
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
