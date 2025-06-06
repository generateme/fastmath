(ns fastmath.dev.clay
  (:require [clojure.string :as str]
            [scicloj.clay.v2.api :as clay]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.kindly.v4.api :as kindly]
            [zprint.core :as zp]
            [clojure.walk :as walk]))

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

(defmacro zp
  [c]
  `(kind/fragment
    [(kind/code ~(str c))
     (kind/code (zp/zprint-str ~c))]))

(defn wrap-into-markdown
  [coll]
  (walk/postwalk (fn [form] (if (string? form) (kind/md form) form)) coll))

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
  `(kind/table ~(mapv (fn [p g]
                        [(kind/hiccup [:dl [:dt nil `(kind/code (str ~p))] [:dd nil `(kind/table [~g])]])])
                      (mapcat identity (take-nth 2 rows)) (mapcat identity (take-nth 2 (rest rows))))))


(defn gen-constants
  ([ns] (gen-constants ns #{}))
  ([ns skip]
   (let [lst (->> (ns-publics ns)
                  (filter (comp :const meta second))
                  (remove (comp skip first))
                  (sort-by first)
                  (map (fn [[s v]] [(name s) (var-get v) (kind/md (:doc (meta v)))])))]
     (-> {:column-names ["Constant symbol" "Value" "Description"]
          :row-vectors lst}
         (kind/table)
         (kindly/hide-code)))))

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
                      :title "Fastmath documentation"}})
  )
