(ns fastmath.dev.clay
  (:require [clojure.string :as str]
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
