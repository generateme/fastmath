(ns fastmath.dev.llm
  (:require [clojure.string :as str]
            [fastmath.core]
            [fastmath.complex]
            [fastmath.quaternion]
            [fastmath.stats]
            [fastmath.vector]
            [fastmath.matrix]
            [fastmath.random]
            [clojure.data.json :as json]))

(defn generate-toc
  [namespace]
  (for [sv (->> (ns-publics namespace)
                (map (comp meta second))
                (remove :deprecated))
        :let [m (-> (dissoc sv :line :column :file :ns :tag :inline-arities :inline)
                    (update :doc (comp first (fnil str/split-lines ""))))]]
    m))

(defn info-prompt
  [project]
  (str "Here is the list of all namespaces and public symbols (global vars, functions and macros) in a " project
       " project. Format is a JSON with following structure. Main object contains:

- key: a namespace name
- value: an array containing objects:
     - name: symbol name
     - doc: short documention of a symbol
     - macro: for macros only, marks if symbol represents a macro
     - arglists: for functions/macros list of accepted arguments
     - const: for constant symbols, marks if symbol is a constant

"))

(defn gen-llm-symbols
  ([output project nss]
   (->> (into {} (for [namespace nss]
                   [(name namespace) (generate-toc namespace)]))
        (json/write-str)
        (str (info-prompt project))
        (spit output))))

(gen-llm-symbols "llm_functions_list.txt" "fastmath" '[fastmath.core
                                                       fastmath.stats
                                                       fastmath.matrix
                                                       fastmath.vector
                                                       fastmath.complex
                                                       fastmath.quaternion
                                                       fastmath.random])


(set (mapcat keys (generate-toc 'fastmath.stats)))
;; => #{:name :const :macro :arglists :doc}
