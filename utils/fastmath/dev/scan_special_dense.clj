;; Compares `fastmath.special` with the dense scan grids of `generate_special_dense_scan.py`
;; (audit WP 9.2) and prints, per function and interval, the worst error.
;;
;; The error is relative to the value, except for the oscillating functions where it is relative
;; to the local envelope (so that a zero crossing does not dominate): J0, J1, Y0, Y1, Airy (x < 0),
;; Si, Ci (x > 1), digamma/trigamma/zeta (near their zeros and poles) use `max(|ref|, envelope * 1e-3)`.
;;
;; Run in a REPL started in the repo root:
;;   (load-file "utils/fastmath/dev/scan_special_dense.clj")

(require '[clojure.edn :as edn]
         '[fastmath.special :as sp])

(def scan (edn/read-string (slurp "/tmp/special_dense_scan.edn")))

(def fns
  {"bessel-J0" sp/bessel-J0, "bessel-J1" sp/bessel-J1, "bessel-Y0" sp/bessel-Y0, "bessel-Y1" sp/bessel-Y1
   "bessel-I0" sp/bessel-I0, "bessel-I1" sp/bessel-I1, "bessel-K0" sp/bessel-K0, "bessel-K1" sp/bessel-K1
   "digamma" sp/digamma, "trigamma" sp/trigamma, "zeta" (fn [s] (sp/zeta s))
   "airy-Ai" sp/airy-Ai, "airy-Ai'" sp/airy-Ai', "airy-Bi" sp/airy-Bi, "airy-Bi'" sp/airy-Bi'
   "Si" sp/Si, "Ci" sp/Ci, "E1" sp/E1, "Ei" sp/Ei, "si" sp/si, "Cin" sp/Cin, "E0" sp/E0, "Ein" sp/Ein, "li" sp/li, "Li" sp/Li})

(defn envelope ^double [name ^double x]
  (let [ax (Math/abs x)]
    (case name
      ("bessel-J0" "bessel-J1" "bessel-Y0" "bessel-Y1") (if (< ax 1.0) 1.0 (Math/sqrt (/ 2.0 (* Math/PI ax))))
      ("airy-Ai" "airy-Bi") (if (< x -1.0) (/ (Math/pow ax -0.25) (Math/sqrt Math/PI)) 0.0)
      ("airy-Ai'" "airy-Bi'") (if (< x -1.0) (/ (Math/pow ax 0.25) (Math/sqrt Math/PI)) 0.0)
      ("Si" "Ci" "si") (if (> ax 1.0) (/ 1.0 ax) 0.0)
      0.0)))

(defn err ^double [name ^double x ^double v ^double r]
  (let [scale (max (Math/abs r) (* 1e-3 (envelope name x)))]
    (cond (and (Double/isNaN v) (Double/isNaN r)) 0.0
          (Double/isNaN v) ##Inf
          (== v r) 0.0
          (zero? scale) (Math/abs (- v r))
          :else (/ (Math/abs (- v r)) scale))))

(def decade-edges [0.0 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 700.0 1e4 1e6 1e8])

(doseq [[name pts] (sort-by key scan)
        :let [f (fns name)
              res (for [[x r] pts
                        :let [v (try (double (f x)) (catch Throwable e ##NaN))]]
                    [(err name x v r) x v r])
              worst (apply max-key first res)
              by-interval (group-by (fn [[_ x]] (last (take-while #(<= % (Math/abs ^double x)) decade-edges))) res)]]
  (println (format "%-10s n=%-6d worst=%.2e at x=%s (got %s ref %s)" name (count res) (first worst) (second worst) (nth worst 2) (nth worst 3)))
  (doseq [[lo rs] (sort-by (comp #(or % 0.0) key) by-interval)
          :let [w (apply max-key first rs)]
          :when (> (first w) 3e-15)]
    (println (format "    |x| >= %-8s worst=%.2e at x=%s" lo (first w) (second w)))))

;; ---- real orders: [order x reference]
(def scan-orders (edn/read-string (slurp "/tmp/special_dense_scan_orders.edn")))

(doseq [[name pts] (sort-by key scan-orders)
        :let [f (case name "bessel-J" sp/bessel-J "bessel-Y" sp/bessel-Y "bessel-I" sp/bessel-I "bessel-K" sp/bessel-K)
              res (for [[v x r] pts
                        :let [got (try (double (f v x)) (catch Throwable e ##NaN))
                              scale (max (Math/abs r) (if (#{"bessel-J" "bessel-Y"} name) (* 1e-3 (envelope "bessel-J0" x)) 0.0))
                              e (cond (Double/isNaN got) ##Inf (== got r) 0.0 :else (/ (Math/abs (- got r)) scale))]]
                    [e v x got r])
              by-order (group-by second res)]]
  (println (format "%-10s n=%-5d worst=%.2e" name (count res) (first (apply max-key first res))))
  (doseq [[v rs] (sort-by key by-order)
          :let [w (apply max-key first rs)]
          :when (> (first w) 1e-12)]
    (println (format "    order %-6s worst=%.2e at x=%s (got %s ref %s)" v (first w) (nth w 2) (nth w 3) (nth w 4)))))
