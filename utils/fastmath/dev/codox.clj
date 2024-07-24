(ns fastmath.dev.codox
  (:require [clojure.string :as str]
            [scicloj.kindly.v4.kind :as kind]))

(def source-files "https://github.com/generateme/fastmath/blob/master/src/")

(defn args->call [f args] (conj (seq args) f))
(defn fix-tex [s]  (str/replace s #"\\\\\(|\\\\\)" "\\$"))
(defn fix-anchor [s] (str/replace s #"\[\[(.+?)\]\]" "[$1](#LOS-$1)"))

(defn remover
  [exclude-vars]
  (->> exclude-vars
       (map (fn [cnd]
              (if (symbol? cnd)
                #{cnd}
                #(re-matches cnd (str %)))))
       (apply juxt (constantly nil))))

(defn make-public-fns-table-clay
  ([ns] (make-public-fns-table-clay ns nil))
  ([ns {:keys [exclude-vars level]
        :or {level "###"}}]
   (let [l1 (str level " ")
         l2 (str "#" l1)]
     (kind/fragment
      [(kind/md (str l1 (name ns)))
       (kind/md (->> (or (:doc (meta (the-ns ns))) "\n") fix-tex fix-anchor))
       (kind/fragment (for [v (->> (ns-publics ns)
                                   (sort-by first)
                                   (remove (comp (partial some identity)
                                                 (remover exclude-vars)
                                                 first))
                                   (map second))
                            :let [{:keys [name file macro arglists doc line const deprecated]} (meta v)]]
                        (kind/fragment
                         (remove nil? [(kind/hiccup [:span {:id (str "#LOS-" (clojure.core/name name))}])
                                       (kind/md (str l2 name
                                                     (when deprecated " ^~DEPRECATED~^")
                                                     (when macro " ^~MACRO~^")
                                                     (when const " ^~CONST~^")))
                                       (when (and deprecated (string? deprecated))
                                         (kind/md
                                          (str "*" (str "Deprecated: "
                                                        (->> deprecated fix-tex fix-anchor)) "*")))
                                       (when const (kind/md (str "`" name " = " (var-get v) "`")))
                                       (when arglists (kind/md
                                                       (str/join "\n"
                                                                 (for [args arglists]
                                                                   (str "+ `" (args->call name args) "`")))))
                                       (when doc (kind/md (->> doc fix-tex fix-anchor)))
                                       (kind/hiccup [:div {:style "text-align: right"}
                                                     [:small
                                                      [:a {:href (str source-files file "#L" line)} "source"]]
                                                     [:hr {:style "margin: 0"}]])]))))]))))

