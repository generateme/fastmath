(ns fastmath.kernel.bandwidth-test
  "Tests of the bandwidth selectors of `fastmath.kernel.density`.
  Oracles: R `bw.nrd0`/`bw.SJ` (utils/fastmath/dev/generate_bandwidth_reference.R, .py: `bandwidth_reference.edn`) and an independent
  numpy implementation of the cross validation criteria with a fine grid optimum (`bandwidth_reference.edn`, `rlcv_reference.edn`)."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.kernel.density :as d]))

(def ^:private bandwidth-reference
  (delay (edn/read-string (slurp (io/resource "kernel/bandwidth_reference.edn")))))

(def ^:private rlcv-reference
  (delay (edn/read-string (slurp (io/resource "kernel/rlcv_reference.edn")))))

(defn- rel-error
  "`abs(a - b) / abs(b)`"
  ^double [a b]
  (/ (m/abs (- (double a) (double b))) (m/abs (double b))))

(defn- invalid-input?
  [property thunk]
  (try (thunk) false
       (catch clojure.lang.ExceptionInfo e
         (let [data (ex-data e)]
           (and (= :invalid-input (:type data))
                (= property (:property data))
                (contains? data :function)
                (contains? data :value)
                (contains? data :expected))))))

(defn- value-of
  "The number of an R outcome `{:value v}` or `nil` for `{:error true}`."
  [outcome]
  (when-not (:error outcome) (:value outcome)))

;; ---------------------------------------------------------------------------------------------
;; rules of thumb

(deftest nrd-reference-test
  (testing ":nrd and :nrd0 equal R bw.nrd0 (the type 7 IQR, the fallback chain of R) on 32 datasets (KB-07, KB-08)"
    (doseq [{:keys [name data n py nrd]} (:cases @bandwidth-reference)
            :when (>= n 2)]
      (is (< (rel-error (d/bandwidth :gaussian data :nrd) (:value (:nrd-chain py))) 1.0e-12) (str name))
      (is (< (rel-error (d/bandwidth :gaussian data :nrd0) (:value (:nrd0-chain py))) 1.0e-12) (str name))
      (when-let [r (value-of nrd)]
        (when (pos? r)
          (is (< (rel-error (d/bandwidth :gaussian data :nrd) r) 1.0e-12) (str name " vs R bw.nrd"))))))
  (testing "values of the type 7 quantile: [-1 0 1] gives 0.635 (the type 6 quantile gave 0.851)"
    (is (m/delta-eq 0.6350045190043616 (d/bandwidth :gaussian [-1 0 1] :nrd) 1.0e-12))
    (is (m/delta-eq (* 0.9 (/ 1.5 1.34) (Math/pow 4 -0.2)) (d/bandwidth :gaussian [0 1 2 3] :nrd0) 1.0e-12)
        "min(sd = 1.29, IQR/1.34 = 1.5/1.34)")))

(deftest nrd-fallback-chain-test
  (testing "zero spread: sd, then the first value, then 1"
    (let [h (fn [data] (d/bandwidth :gaussian data :nrd))]
      (is (m/delta-eq (* 1.06 (Math/sqrt (/ 18.0 99)) (Math/pow 100 -0.2)) ;; 98 zeros, one 3 and one -3: IQR is 0, sd > 0
                      (h (concat (repeat 98 0.0) [3.0 -3.0])) 1.0e-12))
      (is (m/delta-eq (* 1.06 2.5 (Math/pow 3 -0.2)) (h [2.5 2.5 2.5]) 1.0e-12))
      (is (m/delta-eq (* 1.06 1.0 (Math/pow 3 -0.2)) (h [0.0 0.0 0.0]) 1.0e-12))))
  (testing "n < 2 throws"
    (doseq [selector [:nrd :nrd0 :nrd-adjust :sj :sj-dpi :lscv :lcv :rlcv]]
      (is (invalid-input? :data #(d/bandwidth :gaussian [1.0] selector)) (str selector)))
    (is (invalid-input? :data #(d/bandwidth :gaussian [] :nrd)))))

(deftest nrd-adjust-test
  (testing "the 12 kernels with a canonical bandwidth: nrd * nrd-scale"
    (let [data [1.0 2.0 2.5 3.0 4.0 7.0 7.5 3.3]]
      (doseq [kernel (keys d/kde-data)
              :when (not (#{:silverman :cauchy} kernel))]
        (is (m/delta-eq (* (d/bandwidth kernel data :nrd) (:nrd-scale (d/kde-data kernel)))
                        (d/bandwidth kernel data :nrd-adjust) 1.0e-14) (str kernel)))))
  (testing "KB-10: silverman, cauchy, custom kernel and an unknown kernel throw (no silent :nrd)"
    (let [data [1.0 2.0 2.5 3.0 4.0 7.0]]
      (is (invalid-input? :bandwidth #(d/bandwidth :silverman data :nrd-adjust)))
      (is (invalid-input? :bandwidth #(d/bandwidth :cauchy data :nrd-adjust)))
      (is (invalid-input? :bandwidth #(d/bandwidth (fn [x] (d/gaussian x)) data :nrd-adjust)))
      (is (invalid-input? :kernel #(d/bandwidth :unknown data :nrd-adjust))))))

(deftest bandwidth-arguments-test
  (let [data [1.0 2.0 2.5 3.0 4.0 7.0]]
    (testing "a number is returned as a double, nil is :nrd"
      (is (= 0.25 (d/bandwidth :gaussian data 0.25)))
      (is (= 2.0 (d/bandwidth :gaussian data 2)))
      (is (= (d/bandwidth :gaussian data :nrd) (d/bandwidth :gaussian data nil))))
    (testing "invalid selectors and numbers"
      (is (invalid-input? :bandwidth #(d/bandwidth :gaussian data :foo)))
      (doseq [bad [0 -1.0 ##NaN ##Inf "a"]]
        (is (invalid-input? :bandwidth #(d/bandwidth :gaussian data bad)) (pr-str bad))))
    (testing "the input is neither sorted nor modified"
      (let [a (double-array [3.0 1.0 2.0 9.0 4.0])]
        (d/bandwidth :gaussian a :lscv)
        (is (= [3.0 1.0 2.0 9.0 4.0] (vec a)))))
    (testing "the wrapper in fastmath.kernel"
      (is (= (d/bandwidth :gaussian data :nrd) (k/bandwidth :gaussian data :nrd)))
      (is (= (d/bandwidth :gaussian data :lscv {:info? true}) (k/bandwidth :gaussian data :lscv {:info? true}))))))

;; ---------------------------------------------------------------------------------------------
;; Sheather-Jones

(deftest sj-reference-test
  (testing ":sj and :sj-dpi equal R bw.SJ(method = 'ste' / 'dpi') on 32 datasets, errors of R are errors here (D7)"
    (doseq [{:keys [name data sj-ste sj-dpi]} (:cases @bandwidth-reference)]
      (if-let [r (value-of sj-ste)]
        (is (< (rel-error (d/bandwidth :gaussian data :sj) r) 1.0e-10) (str name " ste"))
        (is (invalid-input? :data #(d/bandwidth :gaussian data :sj)) (str name " ste: R fails")))
      (if-let [r (value-of sj-dpi)]
        (is (< (rel-error (d/bandwidth :gaussian data :sj-dpi) r) 1.0e-10) (str name " dpi"))
        (is (invalid-input? :data #(d/bandwidth :gaussian data :sj-dpi)) (str name " dpi: R fails")))))
  (testing "other kernels: the gaussian result times the canonical bandwidth factor"
    (let [{:keys [data]} (first (filter #(= "lognorm100" (:name %)) (:cases @bandwidth-reference)))]
      (doseq [kernel [:epanechnikov :triangular :cosine :logistic :laplace]]
        (is (m/delta-eq (* (:nrd-scale (d/kde-data kernel)) (d/bandwidth :gaussian data :sj))
                        (d/bandwidth kernel data :sj) 1.0e-14))
        (is (m/delta-eq (* (:nrd-scale (d/kde-data kernel)) (d/bandwidth :gaussian data :sj-dpi))
                        (d/bandwidth kernel data :sj-dpi) 1.0e-14)))
      (is (invalid-input? :bandwidth #(d/bandwidth :silverman data :sj)))
      (is (invalid-input? :bandwidth #(d/bandwidth :cauchy data :sj-dpi)))))
  (testing "the equation is solved with the tolerance of R (0.1 * lower bound): the tight solution differs by up to 1 percent"
    (doseq [{:keys [name sj-ste sj-ste-tight]} (:cases @bandwidth-reference)
            :when (and (value-of sj-ste) (value-of sj-ste-tight))]
      (is (< (rel-error (value-of sj-ste) (value-of sj-ste-tight)) 0.02) (str name)))))

;; ---------------------------------------------------------------------------------------------
;; cross validation

(deftest cv-gaussian-battery-test
  (testing "gaussian kernel, 32 datasets (sizes 20 .. 500, scales 0.01 .. 1000, offset 1e6, skewed, heavy tails, ties):
            the optimum of the fine grid of the independent implementation.
            lscv, lcv and rlcv :continuous agree within 1e-3 (measured 1e-4), rlcv :printed within 3% (measured 2.2%:
            its criterion jumps at f = a_n, so the optimum is not smooth)"
    (doseq [{:keys [name data] :as c} (:cv-gaussian-battery @bandwidth-reference)
            [criterion selector opts tol] [[:lscv :lscv nil 1.0e-3]
                                           [:lcv :lcv nil 1.0e-3]
                                           [:rlcv-continuous :rlcv {:rlcv-form :continuous} 1.0e-3]
                                           [:rlcv-printed :rlcv {:rlcv-form :printed} 0.03]]
            :let [info (d/bandwidth :gaussian data selector (assoc opts :info? true))
                  reference (get c criterion)]]
      (is (<= (rel-error (:h info) (:h reference)) tol) (str name " " criterion " " (:h info) " " (:h reference)))
      (is (= (:at-bound reference) (:at-bound? info)) (str name " " criterion " at-bound?"))
      (is (= :ok (:status info)) (str name " " criterion))))
  (testing "KB-01: sd 1000, the search was in the wrong units (an h of the order 6e4 for a similar sample)"
    (let [{:keys [data lscv sd]} (first (filter #(= :normal-sd1000-n100 (:name %)) (:cv-gaussian-battery @bandwidth-reference)))
          h (d/bandwidth :gaussian data :lscv)]
      (is (> sd 500.0))
      (is (< (rel-error h (:h lscv)) 1.0e-3))
      (is (< h (* 2.0 sd))))))

(deftest cv-other-kernels-test
  (testing "the 14 kernels on 3 datasets against the optimum of the independent fine grid (3%).
            Not compared: uniform and wigner (a jump and a square root edge: the criteria are not smooth, several
            local optima of the same value), the likelihood criterion of silverman (negative densities, undefined
            where one is negative), triangular on the exponential sample (the least squares criterion changes by 4e-4
            over h in [0.2, 0.6], a flat sawtooth), the printed rlcv which is compared with 8% (measured 5.7%)"
    (let [datasets (into {} (map (juxt :name identity) (:datasets @rlcv-reference)))]
      (doseq [{:keys [dataset kernel] :as c} (:cases @rlcv-reference)
              :when (and (#{:normal-n100 :lognormal-n100 :bimodal-n200} dataset)
                         (not (#{:uniform :wigner} kernel)))
              :let [{:keys [data sd]} (datasets dataset)
                    flat? (and (= kernel :triangular) (= dataset :exponential-n150))]
              [criterion selector opts tol] [[:lscv :lscv nil 0.03]
                                             [:lcv :lcv nil 0.03]
                                             [:rlcv-continuous :rlcv {:rlcv-form :continuous} 0.03]
                                             [:rlcv-printed :rlcv {:rlcv-form :printed} 0.08]]
              :when (not (and (= kernel :silverman) (= criterion :lcv)))
              :let [info (d/bandwidth kernel data selector (assoc opts :info? true))
                    reference (get c criterion)]]
        (is (or flat? (<= (rel-error (:h info) (* sd (:h reference))) tol))
            (str dataset " " kernel " " criterion " " (:h info) " " (* sd (:h reference))))))))

(deftest cv-criteria-values-test
  (testing "the criteria at the fixed bandwidths (standardised data) equal the independent implementation:
            the criteria with a smooth integrand agree to 1e-9 relative (measured 2e-12), the continuous rlcv with a kink in the
            integrand (the bias term at f = a_n) to 1e-7 (measured 1e-8), the printed rlcv with a jump to 5e-4 (measured 1.2e-4,
            the trapezoid error of a jump is of the first order in the step h / 40)"
    (let [datasets (into {} (map (juxt :name identity) (:datasets @rlcv-reference)))]
      (doseq [{:keys [dataset kernel fixed-h fixed]} (:cases @rlcv-reference)
              :when (and (#{:normal-n100 :bimodal-n200} dataset)
                         (#{:gaussian :logistic :sigmoid} kernel))
              :let [{:keys [data mean sd]} (datasets dataset)
                    z (mapv #(/ (- % mean) sd) data)
                    kdata (#'d/->kdata :test z kernel false nil)]
              [i h] (map-indexed vector fixed-h)]
        (doseq [[criterion form key] [[:lscv :printed :lscv]
                                      [:lcv :printed :lcv]
                                      [:rlcv :printed :rlcv-printed]
                                      [:rlcv :continuous :rlcv-continuous]]
                :let [value (#'d/cv-criterion kdata criterion form h)
                      expected (nth (get fixed key) i)]]
          (is (< (rel-error value expected) (case key :rlcv-printed 5.0e-4 :rlcv-continuous 1.0e-7 1.0e-9))
              (str dataset " " kernel " " key " " h " " value " " expected))))))
  (testing "the two forms of the rlcv bias term differ by a/2 times the measure of {f >= a}"
    (let [{:keys [data mean sd]} (first (:datasets @rlcv-reference))
          z (mapv #(/ (- % mean) sd) data)
          kdata (#'d/->kdata :test z :gaussian false nil)
          c (first (filter #(and (= :gaussian (:kernel %)) (= (:name (first (:datasets @rlcv-reference))) (:dataset %)))
                           (:cases @rlcv-reference)))
          a (#'d/rlcv-a (count z))]
      ;; the jump of the printed form: 2.7e-4 measured
      (doseq [[i h] (map-indexed vector (:fixed-h c))
              :let [printed (#'d/cv-criterion kdata :rlcv :printed h)
                    continuous (#'d/cv-criterion kdata :rlcv :continuous h)
                    mass (nth (:measure-above-a (:fixed c)) i)]]
        (is (< (m/abs (- (- continuous printed) (* 0.5 a mass))) 5.0e-4) (str h))))))

(deftest cv-invariance-test
  (testing "h(a x + b) = a h(x): the data are standardised (KB-01, KB-02)"
    (let [x (:data (first (filter #(= :normal-n100 (:name %)) (:cv-gaussian-battery @bandwidth-reference))))]
      (doseq [selector [:lscv :lcv :rlcv]
              kernel [:gaussian :epanechnikov]
              :let [h0 (d/bandwidth kernel x selector)]
              ;; the optimum of a smooth criterion is known to sqrt(eps) = 1e-8: measured 1.2e-7 at most,
              ;; 1e-5 where the offset 1e6 leaves 8 digits of data of the scale 1e-2
              [a b tol] [[1.0e-2 0.0 1.0e-6] [1.0e3 0.0 1.0e-6] [1.0 1.0e6 1.0e-6] [1.0e3 1.0e6 1.0e-6] [1.0e-2 1.0e6 1.0e-5]]
              :let [y (map #(+ (* a %) b) x)
                    h (d/bandwidth kernel y selector)]]
        (is (< (rel-error h (* a h0)) tol) (str kernel " " selector " " a " " b " " (rel-error h (* a h0)))))))
  (testing "the rules of thumb and Sheather-Jones are equivariant"
    (let [x (:data (first (filter #(= "lognorm100" (:name %)) (:cases @bandwidth-reference))))]
      (doseq [selector [:nrd :nrd0 :sj :sj-dpi]]
        (is (< (rel-error (d/bandwidth :gaussian (map #(* 1000.0 %) x) selector)
                          (* 1000.0 (d/bandwidth :gaussian x selector)))
               1.0e-9) (str selector))))))

(deftest cv-determinism-test
  (testing "no random search: the same input gives bitwise the same bandwidth"
    (let [x (:data (first (filter #(= :lognormal-n100 (:name %)) (:cv-gaussian-battery @bandwidth-reference))))]
      (doseq [selector [:lscv :lcv :rlcv]]
        (is (= (d/bandwidth :gaussian x selector) (d/bandwidth :gaussian x selector) (d/bandwidth :gaussian x selector))
            (str selector))))))

(deftest cv-options-test
  (let [x (:data (first (filter #(= :normal-n100 (:name %)) (:cv-gaussian-battery @bandwidth-reference))))]
    (testing ":info? map"
      (let [info (d/bandwidth :gaussian x :lscv {:info? true})]
        (is (= #{:h :value :at-bound? :status} (set (keys info))))
        (is (= (:h info) (d/bandwidth :gaussian x :lscv)))
        (is (= :ok (:status info)))
        (is (false? (:at-bound? info))))
      (let [info (d/bandwidth :gaussian x :nrd {:info? true})]
        (is (= (:h info) (d/bandwidth :gaussian x :nrd)))
        (is (= :closed-form (:status info)))
        (is (false? (:at-bound? info)))))
    (testing "KB-15: the search interval is reported when it clips the optimum"
      (let [free (d/bandwidth :gaussian x :lscv {:info? true})
            scale (d/bandwidth :gaussian x :nrd0 {:info? true})
            clipped (d/bandwidth :gaussian x :lscv {:info? true :bounds [0.01 0.2]})
            clipped-low (d/bandwidth :gaussian x :lscv {:info? true :bounds [2.0 5.0]})]
        (is (true? (:at-bound? clipped)))
        (is (true? (:at-bound? clipped-low)))
        (is (< (:h clipped) (:h free)))
        (is (> (:h clipped-low) (:h free)))
        (is (< 0.0 (:h scale)))))
    (testing "the grid option, a coarse and a fine grid find the same optimum of a smooth criterion"
      (is (< (rel-error (d/bandwidth :gaussian x :lscv {:grid 20}) (d/bandwidth :gaussian x :lscv {:grid 200})) 1.0e-5)))
    (testing "invalid options"
      (doseq [grid [2 0 -1]]
        (is (invalid-input? :grid #(d/bandwidth :gaussian x :lscv {:grid grid})) (str grid)))
      (doseq [bounds [[1.0 0.5] [0.0 1.0] [-1.0 1.0] [1.0] [1.0 ##Inf] 1.0]]
        (is (invalid-input? :bounds #(d/bandwidth :gaussian x :lscv {:bounds bounds})) (pr-str bounds)))
      (is (invalid-input? :rlcv-form #(d/bandwidth :gaussian x :rlcv {:rlcv-form :foo}))))
    (testing "nil option = missing option"
      (is (= (d/bandwidth :gaussian x :rlcv)
             (d/bandwidth :gaussian x :rlcv {:bounds nil :grid nil :bins nil :rlcv-form nil :info? nil :radius nil :extent nil}))))
    (testing "constant data cannot be cross validated"
      (is (invalid-input? :data #(d/bandwidth :gaussian [2.0 2.0 2.0 2.0] :lscv))))
    ;; the sum over every point (radius Inf) against the radius 6.65: 8e-8 measured
    (testing "a custom kernel function works with the cross validation"
      (is (< (rel-error (d/bandwidth (fn [t] (d/gaussian t)) x :lscv) (d/bandwidth :gaussian x :lscv)) 1.0e-6)))))

(deftest cv-degenerate-test
  (testing "lcv of a compact kernel when a point is isolated: the criterion is -Inf for every bandwidth, reported"
    (let [x (concat (repeatedly 99 #(* 0.01 (rand))) [100.0])
          info (d/bandwidth :epanechnikov x :lcv {:info? true :bounds [0.01 0.05]})]
      (is (= :degenerate (:status info)))
      (is (true? (:at-bound? info)))
      (is (pos? (:h info)))))
  (testing "heavy ties: lscv drives the bandwidth to the lower end of the search interval (KB-09)"
    (let [{:keys [data lscv]} (first (filter #(= :ties-integers-n100 (:name %)) (:cv-gaussian-battery @bandwidth-reference)))
          info (d/bandwidth :gaussian data :lscv {:info? true})]
      (is (true? (:at-bound? info)))
      (is (true? (:at-bound lscv))))))

(deftest ns-publics-test
  (testing "the targets of the old optimizer are not public any more (D4)"
    (is (not-any? (set (keys (ns-publics 'fastmath.kernel.density))) '[lscv-target lcv-target rlcv-target]))))
