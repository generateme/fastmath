;; Step 1 of 2 of the noise reference pipeline (Random Audit, work package 9.4).
;;
;; The pipeline
;;   1. this script  : dumps the seed-dependent lookup tables of the library (the permutation `perm`, `perm12` = perm mod 12
;;                     and the value-noise table `valueLUT`, produced by `java.util.Random` inside `fastmath.java.noise.NoiseConfig`)
;;                     through the public fields of `NoiseConfig` into a small JSON file. It evaluates NO noise function.
;;   2. generate_random_noise_reference.py : an independent Python re-implementation of the noise algorithms (Gustavson simplex,
;;                     Perlin-style fade/lerp gradient and value noise, libnoise-style fbm/billow/ridged multifractal) which reads the
;;                     tables and writes test/resources/random/noise_reference.edn.
;;
;; The tables are the only thing taken from the library: the hash, the gradient sets and the lattice walk are re-implemented
;; from the published formulas. Both steps are deterministic; the EDN file is byte-identical on rerun.
;;
;; Run (from the repository root):
;;   java -cp "$(cat /tmp/fm_cp.txt)" clojure.main utils/fastmath/dev/generate_random_noise_reference.clj [tables.json]
;;   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_noise_reference.py [tables.json]
;; (`/tmp/fm_cp.txt` is the project classpath, e.g. `lein classpath`, including target/classes with the compiled Java.)
;; Default tables file: /tmp/fastmath_noise_tables.json.

(ns generate-random-noise-reference
  (:require [clojure.string :as str])
  (:import [fastmath.java.noise NoiseConfig]))

(set! *warn-on-reflection* true)

(def seed 12345)

(def out-file (or (first *command-line-args*) "/tmp/fastmath_noise_tables.json"))

(defn- js
  "JSON text of numbers and sequences (doubles are printed round-trip exact)."
  [x]
  (cond
    (number? x) (str x)
    (string? x) (str "\"" x "\"")
    :else (str "[" (str/join "," (map js x)) "]")))

(let [c (NoiseConfig. (int seed) (int NoiseConfig/NOISE_VALUE) (int NoiseConfig/INTERPOLATE_HERMITE)
                      (int 5) 2.1 0.45 true)
      other (NoiseConfig. (int seed) (int NoiseConfig/NOISE_SIMPLEX) (int NoiseConfig/INTERPOLATE_NONE)
                          (int 6) 2.0 0.5 false)
      perm (vec (take 256 (.perm c)))
      perm12 (vec (take 256 (.perm12 c)))
      lut (vec (.valueLUT c))]
  ;; the tables depend on the seed only
  (assert (= perm (vec (take 256 (.perm other)))))
  (assert (= lut (vec (.valueLUT other))))
  (assert (= perm (vec (drop 256 (.perm c)))) "perm is stored twice")
  (assert (= perm12 (vec (drop 256 (.perm12 c)))))
  (spit out-file (str "{\"seed\":" seed ",\"perm\":" (js perm) ",\"perm12\":" (js perm12) ",\"valueLUT\":" (js lut) "}\n"))
  (println "written" out-file "perm" (count perm) "perm12" (count perm12) "valueLUT" (count lut)))

(shutdown-agents)
