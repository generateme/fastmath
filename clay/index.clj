^:kindly/hide-code
(ns index
  (:require [scicloj.kindly.v4.kind :as kind]
            [scicloj.kindly.v4.api :as kindly]
            [scicloj.clay.v2.api :as clay]))

^:kindly/hide-code
(def md (comp kindly/hide-code
            kind/md))

(md "# Preface ")

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


;; some

^:kindly/hide-code
(comment
  (clay/make! {:source-path ["index.clj"
                             "kernel.clj"]
               :book {:title "Fastmath documentation"}}))
