(ns fastmath.dev.special-ref-gen
  "One-off generator of R reference values for `fastmath.special-test`.

  Produces `test/resources/special/reference.edn`, a map keyed by test name
  containing, per comparison block, the exact input grids (`:order`, `:arg`)
  and the R reference outputs (`:ref`). The test simply maps the `sut` function
  over the stored inputs and compares to `:ref`, so no live R (clojisr) session
  is needed to run the suite.

  Reference source per input:
   - x >= 0 : R stats `besselJ/Y/I/K` (vectorised, element-wise recycling)
   - x <  0 : `Bessel::Bessel{J,Y,I,K}` (per element, handles negative argument)
   - Airy   : `Bessel::Airy{A,B}` with `:deriv`.

  Run from a dev REPL:
    (require '[fastmath.dev.special-ref-gen :as g] :reload)
    (g/-main)"
  (:require [clojisr.v1.r :as rr]
            [fastmath.core :as m]
            [clojure.java.io :as io]
            [clojure.pprint :as pp]))

(rr/require-r '[Bessel] '[base])

(defn- rvec [rres] (mapv double (rr/r->clj rres)))

(defn- bessel-scalar
  "Single Bessel:: evaluation returning a double (used for negative arguments)."
  [rbess a b]
  (double (first (rr/r->clj (rbess (double a) (double b))))))

(defn- ref-vals
  "Element-wise reference for a bessel family.
   `args` are the x-values, `orders` are the nu-values (same length).
   Uses `rbase` (R stats, vectorised) when all args are non-negative,
   otherwise `rbess` (Bessel::) per element."
  [rbase rbess args orders]
  (let [args (mapv double args)
        orders (mapv double orders)]
    (if (every? #(>= (double %) 0.0) args)
      (rvec (rbase args orders))
      (mapv (partial bessel-scalar rbess) args orders))))

;; family -> [rbase rbess]
(def ^:private families
  {:J [base/besselJ Bessel/BesselJ]
   :Y [base/besselY Bessel/BesselY]
   :I [base/besselI Bessel/BesselI]
   :K [base/besselK Bessel/BesselK]})

(defn- block
  "Build a {:order :arg :ref} block from x-values `xs` and nu-values `nus`.
   `argf` maps [x nu] -> argument passed to the special fn (e.g. x*nu or x).
   Iteration order matches `(doseq [x xs nu nus] ...)` (x outer)."
  [family xs nus argf]
  (let [[rbase rbess] (families family)
        pairs (for [x xs nu nus] [(double nu) (double (argf x nu))])
        orders (mapv first pairs)
        args (mapv second pairs)]
    {:order orders :arg args :ref (ref-vals rbase rbess args orders)}))

(defn- vblock
  "Fixed x, varying nu (order sweep)."
  [family order-seq xval]
  (let [[rbase rbess] (families family)
        orders (mapv double order-seq)
        args (mapv (constantly (double xval)) orders)]
    {:order orders :arg args :ref (ref-vals rbase rbess args orders)}))

(defn- sblock
  "Single-order family (J0/J1/... ): fixed integer order, varying x."
  [family order xs]
  (let [[rbase rbess] (families family)
        args (mapv double xs)
        orders (mapv (constantly (double order)) args)]
    {:arg args :ref (ref-vals rbase rbess args orders)}))

(defn- pblock
  "Explicit (order,arg) pairs (from the `are` tables)."
  [family pairs]
  (let [[rbase rbess] (families family)
        orders (mapv (comp double first) pairs)
        args (mapv (comp double second) pairs)]
    {:order orders :arg args :ref (ref-vals rbase rbess args orders)}))

(defn- airy-block [rf xs d]
  (let [args (mapv double xs)]
    {:arg args :ref (rvec (rf args :deriv d))}))

;; ---- input grids (reduced but spanning the original ranges) ----

(def ^:private xpos (vec (range 1.0e-6 100.0 0.5)))     ; 200 pts (was step 0.01)
(def ^:private xneg (mapv - xpos))

(defn build []
  {;; -------- single-order families --------
   :bessel-J0 {:pos (sblock :J 0 xpos)
               :neg {:arg xneg :ref (rvec (Bessel/BesselJ xneg 0.0))}}
   :bessel-J1 {:pos (sblock :J 1 xpos)
               :neg {:arg xneg :ref (rvec (Bessel/BesselJ xneg 1.0))}}
   :bessel-Y0 {:pos (sblock :Y 0 xpos)}
   :bessel-Y1 {:pos (sblock :Y 1 xpos)}
   :bessel-I0 {:pos (sblock :I 0 xpos)
               :neg {:arg xneg :ref (rvec (Bessel/BesselI xneg 0.0))}}
   :bessel-I1 {:pos (sblock :I 1 xpos)
               :neg {:arg xneg :ref (rvec (Bessel/BesselI xneg 1.0))}}
   :bessel-K0 {:pos (sblock :K 0 xpos)}
   :bessel-K1 {:pos (sblock :K 1 xpos)}

   ;; -------- bessel-J --------
   :bessel-J
   (let [xs-i [0.1 0.5 0.9 0.95 0.99 1.0 1.01 1.05]
         nus-i [2 6 15 30 60 100]
         xs-f [0.1 0.5 0.9 0.99 1.0 1.05 1.2 1.5 2.0 3.0]
         nus-f [0.1 0.8123 1.5 4.1234 12.3 28.2345 51.23 80.5 104.2]
         xs-l [0.5 0.8 0.9 0.95 0.99 1.0 1.01 1.1 1.2]
         nus-l [150 200 500 1000 10000 50000]
         xs-n [0.05 0.3 0.6 0.9 1.0 1.2 2.0 3.0 5.1]
         nus-n [-100 -80 -60 -40 -20]]
     {:int-xx (block :J xs-i nus-i (fn [x nu] (* x nu)))
      :int-x  (block :J xs-i nus-i (fn [x _] x))
      :frac-xx (block :J xs-f nus-f (fn [x nu] (* x nu)))
      :frac-x  (block :J xs-f nus-f (fn [x _] x))
      :large-xx (block :J xs-l nus-l (fn [x nu] (* x nu)))
      :large-x  (block :J xs-l nus-l (fn [x _] x))
      :vs-pos-015 (vblock :J (range 0.0 250.0 2.5) 0.15)
      :vs-pos-21  (vblock :J (range 0.0 250.0 2.5) 2.1)
      :vs-pos-421 (vblock :J (range 0.0 250.0 2.5) 42.1)
      :vs-pos-1421 (vblock :J (range 0.0 250.0 2.5) 142.1)
      :vs-neg-015 (vblock :J (range -100.0 0.25) 0.15)
      :vs-neg-21  (vblock :J (range -100.0 0.25) 2.1)
      :vs-neg-421 (vblock :J (range -100.0 0.25) 42.1)
      :vs-neg-1421 (vblock :J (range -100.0 0.25) 142.1)
      :neg-ord (block :J xs-n nus-n (fn [x nu] (* x nu)))
      :are (pblock :J [[-5.0 -5.1] [-7.3 19.1] [-14.0 21.3] [-13.0 21.3]
                       [-14.0 -21.3] [-13.0 -21.3] [7.3 19.1] [14.0 21.3]
                       [13.0 21.3] [14.0 -21.3] [13.0 -21.3]])})

   ;; -------- bessel-Y --------
   :bessel-Y
   (let [xs-i [0.05 0.4 0.7 0.9 0.99 1.0 1.05 1.5 3.0 10.0]
         nus-i [0 1 2 6 20 50 100 200]
         xs-f [0.05 0.3 0.6 0.9 1.0 1.1 1.5 2.0 3.0 5.1]
         nus-f [0.1 0.8123 1.5 4.1234 12.3 28.2345 51.23 80.5]
         xs-l [0.5 0.8 0.9 0.95 0.99 1.0 1.01 1.1 1.2]
         nus-l [150 200 500 1000 20000]]
     {:int-xx (block :Y xs-i nus-i (fn [x nu] (* x nu)))
      :frac-xx (block :Y xs-f nus-f (fn [x nu] (* x nu)))
      :frac-x  (block :Y xs-f nus-f (fn [x _] x))
      :large-xx (block :Y xs-l nus-l (fn [x nu] (* x nu)))
      :vs-pos-015 (vblock :Y (range 0.0 100.0 2.0) 0.15)
      :vs-pos-21  (vblock :Y (range 0.0 100.0 2.0) 2.1)
      :vs-pos-421 (vblock :Y (range 0.0 100.0 2.0) 42.1)
      :vs-pos-1421 (vblock :Y (range 0.0 100.0 2.0) 142.1)
      :vs-neg-015 (vblock :Y (range -100.01 0.0 2.0) 0.15)
      :vs-neg-21  (vblock :Y (range -100.01 0.0 2.0) 2.1)
      :vs-neg-421 (vblock :Y (range -100.01 0.0 2.0) 42.1)
      :vs-neg-1421 (vblock :Y (range -100.01 0.0 2.0) 142.1)
      :are (pblock :Y [[-6.2 18.6] [-8.0 23.2] [-7.0 23.2] [-6.0 23.2] [-0.1 2.2]])})

   ;; -------- bessel-I --------
   :bessel-I
   (let [xs-g [0.05 0.4 0.7 0.9 1.0 1.1 1.5 2.0 3.0 4.0]
         nus-g [0.01 0.5 1 2 5.23 10 20 50 100 160.789]
         xs-f [0.05 0.3 0.6 0.9 1.0 1.1 1.5 2.0 3.0 5.1]
         nus-f [0.1 0.8123 1.5 4.1234 12.3 28.2345 51.23 80.5]
         xs-n [0.05 0.3 0.6 0.9 1.0 1.2 2.0 3.0 5.1]
         nus-n [-100 -80 -60 -40 -20]]
     {:gen-xx (block :I xs-g nus-g (fn [x nu] (* x nu)))
      :gen-x  (block :I xs-g nus-g (fn [x _] x))
      :frac-xx (block :I xs-f nus-f (fn [x nu] (* x nu)))
      :frac-x  (block :I xs-f nus-f (fn [x _] x))
      :vs-pos-015 (vblock :I (range 0.0 250.0 2.5) 0.15)
      :vs-pos-21  (vblock :I (range 0.0 250.0 2.5) 2.1)
      :vs-pos-421 (vblock :I (range 0.0 250.0 2.5) 42.1)
      :vs-pos-1421 (vblock :I (range 0.0 250.0 2.5) 142.1)
      :vs-neg-015 (vblock :I (range -100.0 0.0 2.0) 0.15)
      :vs-neg-21  (vblock :I (range -100.0 0.0 2.0) 2.1)
      :vs-neg-421 (vblock :I (range -100.0 0.0 2.0) 42.1)
      :vs-neg-1421 (vblock :I (range -100.0 0.0 2.0) 142.1)
      :negord-xx (block :I xs-n nus-n (fn [x nu] (* x nu)))
      :negord-axx (block :I xs-n nus-n (fn [x nu] (m/abs (* x nu))))
      :are (pblock :I [[12.0 3.2] [13.0 -1.0] [-8.0 4.2] [12.3 8.2] [-12.3 8.2] [-14.0 -9.9]])})

   ;; -------- bessel-K --------
   :bessel-K
   (let [xs-f [0.02 0.1 0.3 0.6 0.9 1.0 1.1 1.5 2.0 3.0 5.1]
         nus-f [0.1 0.8123 1.5 4.1234 12.3 28.2345 51.23 72.23435]
         xs-n [0.05 0.3 0.6 0.9 1.0 1.2 2.0 3.0 5.1]]
     {:nu-sweep (let [nus (range -36.0 82.0 6.0)
                      xs (rest (m/slice-range 0.0 30.0 16))]
                  ;; doseq [nu nus x xs] -> nu outer
                  (let [pairs (for [nu nus x xs] [(double nu) (double x)])]
                    {:order (mapv first pairs) :arg (mapv second pairs)
                     :ref (ref-vals base/besselK Bessel/BesselK
                                    (mapv second pairs) (mapv first pairs))}))
      :frac-xx (block :K xs-f nus-f (fn [x nu] (* x nu)))
      :frac-x  (block :K xs-f nus-f (fn [x _] x))
      :vs-pos-015 (vblock :K (range 0.0 100.0 2.0) 0.15)
      :vs-pos-21  (vblock :K (range 0.0 100.0 2.0) 2.1)
      :vs-pos-421 (vblock :K (range 0.0 100.0 2.0) 42.1)
      :vs-pos-1421 (vblock :K (range 0.0 100.0 2.0) 142.1)
      :vs-neg-015 (vblock :K (range -100.0 0.0 2.0) 0.15)
      :vs-neg-21  (vblock :K (range -100.0 0.0 2.0) 2.1)
      :vs-neg-421 (vblock :K (range -100.0 0.0 2.0) 42.1)
      :vs-neg-1421 (vblock :K (range -100.0 0.0 2.0) 142.1)
      :negord-axx (block :K xs-n (range -50.0 0.0 2.0) (fn [x nu] (m/abs (* x nu))))
      :are (pblock :K [[12.0 3.2] [-8.0 4.2] [12.3 8.2] [-12.3 8.2]])})

   ;; -------- airy --------
   :airy
   (let [rsmall (range -20.00001 21.00001 0.25)
         rlargea (range -12345.5 5678.5 91.23)
         rlargeb (range -12345.5 101.5 91.23)
         rlargea' (range -3234.5 5678.5 91.23)]
     {:ai-small (airy-block Bessel/AiryA rsmall 0)
      :bi-small (airy-block Bessel/AiryB rsmall 0)
      :ai'-small (airy-block Bessel/AiryA rsmall 1)
      :bi'-small (airy-block Bessel/AiryB rsmall 1)
      :ai-large (airy-block Bessel/AiryA rlargea 0)
      :bi-large (airy-block Bessel/AiryB rlargeb 0)
      :ai'-large (airy-block Bessel/AiryA rlargea' 1)
      :bi'-large (airy-block Bessel/AiryB rlargeb 1)})})

(defn -main [& _]
  (let [data (build)
        f (io/file "test/resources/special/reference.edn")]
    (io/make-parents f)
    (binding [*print-length* nil]
      (spit f (with-out-str (pp/pprint data))))
    (println "Wrote" (str f))
    (doseq [[k v] data]
      (println (format "  %-12s blocks=%d vals=%d"
                       (name k) (count v)
                       (reduce + (map (comp count :ref) (vals v))))))
    (println "done")))
