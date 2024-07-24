^:kindly/hide-code
(ns index
  (:require [scicloj.kindly.v4.kind :as kind]
            [scicloj.kindly.v4.api :as kindly]
            [scicloj.clay.v2.api :as clay]
            
            [fastmath.dev.clay :as utls]
            [fastmath.core :as m]))

^:kindly/hide-code
(def md (comp kindly/hide-code
            kind/md))

;; zzz

(md "# Preface {.unnumbered}")

^{:kindly/kind :kind/table}
[[1 2 3] [:a :b :c]]

^{:kindly/kind :kind/table}
[{:a 1 :b 2} {:c 3 :d 4} {:a 0 :b 0 :c 0 :g 0}]

^{:kindly/kind :kind/table}
{:a [1 2 3 4]
 :b [4 5 6 7]}

^{:kindly/kind :kind/table}
{:row-maps [{:a 1 :b 2} {:c 3 :d 4} {:a 0 :b 0 :c 0 :g 0}]}

^{:kindly/kind :kind/table}
{:row-vectors [[1 2 3] [:a :b :c]]}

^{:kindly/kind :kind/table}
{:column-names [:a :g]
 :row-maps [{:a 1 :b 2} {:c 3 :d 4} {:a 0 :b 0 :c 0 :g 0}]}

^{:kindly/kind :kind/table}
{:column-names ["A" "B" "C"]
 :row-vectors [[1 2 3] [:a :b :c]]}

;; ::: {layout="[50,50]"}

;; ::: {#first-column}

^:kindly/hide-code
(kind/code "(+ 1 2 3)")

^:kindly/hide-code
(+ 1 2 3)

^:kindly/hide-code
(kind/code "(* 3 4 5 6 7)")

^:kindly/hide-code
(* 3 4 5 6 7)

^:kindly/hide-code
(kind/code "(range 25)")

;; :::

;; ::: {#second-column}

^:kindly/hide-code
(range 25)

^:kindly/hide-code
(kind/code
 "(filter (fn [v] (odd? v))
(map (fn [v] (inc v) 
(range 5)))")

^:kindly/hide-code
(filter odd? (map inc (range 5)))

^:kindly/hide-code
(kind/code
 "(filter (fn [v] (odd? v))
(map (fn [v] (inc v) 
(range 5)))")

^:kindly/hide-code
(filter odd? (map inc (range 5)))

;; :::

;; :::



(utls/examples-note (+ 1 2) [1] m/PI)


;; some

^:kindly/hide-code
(comment
  (clay/make! {:source-path ["index.clj"
                             "core.clj"
                             "vector_matrix.clj"
                             "random.clj"
                             "stats.clj"
                             "integration.clj"
                             "complex_quaternions.clj"
                             "special.clj"]
               :book {:favicon "clay/resources/favicon.png"
                      :title "Fastmath documentation"}}))


