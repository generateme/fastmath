;; Dumps the coefficient tables of `fastmath.special.poly` (and the Chebyshev tables used by
;; `fastmath.special.airy` / `fastmath.special.hypergeometric`) to JSON, as the input of
;; `verify_special_poly_tables.py`, which regenerates every table from its closed form
;; (high-precision `mpmath` / exact rational arithmetic) and compares coefficient by coefficient.
;;
;; Run from the repo root:
;;   clojure -M utils/fastmath/dev/dump_special_poly_tables.clj /tmp/special_poly_tables.json
;; (or evaluate it in a REPL with `*command-line-args*` bound).
;;
;; The tables are read from the source text (not from the loaded vars): the `uk-poly-*`, `a-ap-poly-*`
;; and `j0-polys`/`j1-polys` coefficients live in function bodies as literal arguments of `mevalpoly`.

(require '[clojure.java.io :as io]
         '[clojure.string :as str])

(def out-file (or (first *command-line-args*) "/tmp/special_poly_tables.json"))

(defn read-forms [path]
  (with-open [r (java.io.PushbackReader. (io/reader path))]
    (binding [*read-eval* false]
      (loop [forms []]
        (let [f (read {:eof ::eof :read-cond :allow} r)]
          (if (= f ::eof) forms (recur (conj forms f))))))))

(defn num-str [x]
  (cond (instance? Double x) (pr-str x)
        (integer? x) (str x)
        :else (str x)))

(defn mevalpoly-coeffs
  "Coefficient list of a `(poly/mevalpoly var c0 c1 ...)` or `(mevalpoly-m var [ints] denom)` form."
  [form]
  (when (and (seq? form) (symbol? (first form)))
    (case (name (first form))
      "mevalpoly" {:coeffs (mapv num-str (drop 2 form))}
      "mevalpoly-m" {:coeffs (mapv num-str (nth form 2)) :denominator (num-str (nth form 3))}
      nil)))

(defn find-defn [forms nm]
  (first (filter #(and (seq? %) (#{'defn 'def} (first %)) (= nm (name (second %)))) forms)))

(defn let-tables
  "{binding-name coeffs} for every let binding of a function body whose value is a mevalpoly form."
  [defn-form]
  (let [res (atom {})]
    (clojure.walk/postwalk
     (fn [x]
       (when (and (seq? x) (= 'let (first x)) (vector? (second x)))
         (doseq [[k v] (partition 2 (second x))]
           (when-let [c (mevalpoly-coeffs v)]
             (swap! res assoc (name k) c))))
       x)
     defn-form)
    @res))

(require 'clojure.walk)

(defn case-tables
  "{index coeffs} of a `(case (int n) 0 (poly/mevalpoly r ...) 1 ...)` body."
  [defn-form]
  (let [res (atom {})]
    (clojure.walk/postwalk
     (fn [x]
       (when (and (seq? x) (= 'case (first x)))
         (doseq [[k v] (partition 2 (drop 2 x))]
           (when (and (integer? k) (mevalpoly-coeffs v))
             (swap! res assoc (str k) (mevalpoly-coeffs v)))))
       x)
     defn-form)
    @res))

(defn double-array-literal [forms nm]
  (let [f (find-defn forms nm)
        lst (last f)
        items (if (and (seq? lst) (= 'double-array (first lst)))
                (second lst)
                lst)
        items (if (and (seq? items) (= 'quote (first items))) (second items) items)]
    (mapv num-str items)))

(defn roots [forms nm]
  (->> (last (find-defn forms nm))
       (mapv (fn [v] [(num-str (nth v 1)) (num-str (nth v 2))]))))

(let [poly (read-forms "src/fastmath/special/poly.clj")
      table (fn [nm] (let-tables (find-defn poly nm)))
      data {"uk-poly-10" (table "uk-poly-10")
            "uk-poly-20" (table "uk-poly-20")
            "uk-poly-jn" (table "uk-poly-jn")
            "a-ap-poly-10" (table "a-ap-poly-10")
            "a-ap-poly-20" (table "a-ap-poly-20")
            "a-ap-poly-30" (table "a-ap-poly-30")
            "j0-roots" (roots poly "j0-roots")
            "j1-roots" (roots poly "j1-roots")
            "j0-polys" (case-tables (find-defn poly "j0-polys"))
            "j1-polys" (case-tables (find-defn poly "j1-polys"))
            "abi-1" (double-array-literal poly "abi-1")
            "abi-2" (double-array-literal poly "abi-2")
            "abip-1" (double-array-literal poly "abip-1")
            "abip-2" (double-array-literal poly "abip-2")
            "hg-2F1-poly" (double-array-literal poly "hg-2F1-poly")
            "bessel-y-chebyshev-weights" (let [f (find-defn poly "bessel-y-chebyshev-weights")]
                                           (mapv (fn [da] (mapv num-str (second (second da)))) (last f)))}
      json-str (fn json [x]
                 (cond (map? x) (str "{" (str/join "," (map (fn [[k v]] (str (pr-str (if (keyword? k) (name k) (str k))) ":" (json v))) x)) "}")
                       (sequential? x) (str "[" (str/join "," (map json x)) "]")
                       :else (pr-str (str x))))]
  (spit out-file (json-str data))
  (println "written" out-file (.length (io/file out-file)) "bytes"))
