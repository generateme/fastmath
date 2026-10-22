;; Pre-change baseline of fastmath.random (Random Audit, commit 34b05ef).
;;
;; Records outputs of the code BEFORE the audit changes, so that tests can state which seeded
;; outputs must stay bitwise equal and which are allowed to change (Spec 1.1, 1.3, 3.2, 3.3, 4.1, 5.1, 5.2, 7.x).
;; Run on the unchanged build only:
;;   java -cp "$(cat /tmp/fm_cp.txt)" clojure.main utils/fastmath/dev/generate_random_baseline.clj
;; Output: test/resources/random/baseline_reference.edn (deterministic, byte-identical on rerun).

(ns generate-random-baseline
  (:require [clojure.pprint :as pp]
            [fastmath.random :as r]))

(def algorithms [:jdk :mersenne :isaac :well512a :well1024a :well19937a :well19937c :well44497a :well44497b])

(defn- n [x] (if (number? x) x (vec x)))

(defn- draws [k f] (vec (repeatedly k f)))

(def rng-streams
  (into {}
        (for [a algorithms]
          [a (let [mk #(r/rng a 42)]
               {:irandom-mx (into {} (for [mx [1 2 10 1000000 2147483647]]
                                       [mx (let [g (mk)] (draws 20 #(r/irandom g mx)))]))
                :irandom-range (into {} (for [[lo hi] [[-3 4] [0 100]]]
                                          [[lo hi] (let [g (mk)] (draws 20 #(r/irandom g lo hi)))]))
                :lrandom-mx (into {} (for [mx [1 2 10 4294967303 1000000000000]]
                                       [mx (let [g (mk)] (draws 20 #(r/lrandom g mx)))]))
                :lrandom-range (into {} (for [[lo hi] [[-5 5] [-1000000000000 1000000000000]]]
                                          [[lo hi] (let [g (mk)] (draws 20 #(r/lrandom g lo hi)))]))
                :drandom (let [g (mk)] (draws 20 #(r/drandom g)))
                :drandom-range (let [g (mk)] (draws 20 #(r/drandom g -2.0 6.0)))
                :frandom (let [g (mk)] (draws 10 #(double (r/frandom g))))
                :grandom (let [g (mk)] (draws 20 #(r/grandom g)))
                :grandom-std (let [g (mk)] (draws 20 #(r/grandom g 2.5)))
                :grandom-mean-0 (let [g (mk)] (draws 20 #(r/grandom g 0.0 2.5)))
                :grandom-mean-nonzero (let [g (mk)] (draws 20 #(r/grandom g 3.7 1.3)))
                :brandom (let [g (mk)] (draws 20 #(r/brandom g)))
                :brandom-p (let [g (mk)] (draws 20 #(r/brandom g 0.3)))})])))

(def signal (let [g (r/rng :mersenne 7)] (vec (repeatedly 200 #(r/grandom g)))))

(def series
  {:signal signal
   :ma-1 (vec (r/ma 0.5 signal))
   :ma-2 (vec (r/ma [0.5 -0.3] signal))
   :ma-4 (vec (r/ma [0.4 0.3 -0.2 0.1] signal))
   :ar-1 (vec (r/ar 0.7 signal))
   :ar-1-high (vec (r/ar 0.99 signal))
   :ar-2 (vec (r/ar [0.5 -0.3] signal))
   :ar-3 (vec (r/ar [0.2 0.3 0.1] signal))
   :ar-4 (vec (r/ar [0.3 0.2 0.1 0.05] signal))
   :ar-zero (vec (r/ar 0.0 signal))
   :ar-zeros (vec (r/ar [0.0 0.0] signal))
   :ar-lead-zero (vec (r/ar [0.5 0.0] signal))
   :arma-11 (vec (r/arma 0.7 0.4 signal))
   :arma-21 (vec (r/arma [0.5 0.2] 0.3 signal))
   :arma-22 (vec (r/arma [0.5 -0.2] [0.3 0.1] signal))
   :fi-exact (into {} (for [d [0.1 0.3 0.45 -0.3 0.7]] [d (vec (r/fi d (take 120 signal)))]))
   :fi-trunc (into {} (for [d [0.2 0.4]] [d (vec (r/fi d 20 (take 120 signal)))]))
   :fi-int (into {} (for [d [1 2 3 -1 -2]] [d (vec (r/fi d (take 60 signal)))]))
   :fi-zero (vec (take 20 (r/fi 0.0 signal)))
   :arfima-exact (vec (r/arfima 0.5 0.2 0.3 (take 120 signal)))
   :arfima-trunc (vec (r/arfima [0.5 0.1] 0.3 20 [0.2] signal))})

(def noise-defs
  (for [t [:value :gradient :simplex]
        i [:none :linear :hermite :quintic]
        g [:single :fbm :billow :ridgemulti]
        c [{:octaves 5 :lacunarity 2.1 :gain 0.45 :normalize? true}
           {:octaves 6 :lacunarity 2.0 :gain 0.5 :normalize? false}]]
    [t i g c]))

(def noise-points
  {1 [[0.3] [-4.7] [123.456] [-0.0001]]
   2 [[0.3 0.4] [-4.7 2.2] [123.456 -78.9] [10.0 10.0]]
   3 [[0.3 0.4 0.5] [-4.7 2.2 9.9] [123.4 -78.9 1.1] [0.5 0.5 0.5]]})

(def noise-values
  (vec (for [[t i g c] noise-defs]
         (let [cfg (merge {:seed 12345 :noise-type t :interpolation i} c)
               f ((case g :single r/single-noise :fbm r/fbm-noise :billow r/billow-noise :ridgemulti r/ridgedmulti-noise) cfg)]
           {:cfg cfg :generator g
            :values (into {} (for [[d pts] noise-points]
                               [d (mapv #(double (apply f %)) pts)]))}))))

(def sequences
  (into {}
        (concat
         (for [g [:sobol :halton :r2 :default :gaussian :sphere :ball] d [1 2 3 6]]
           [[g d] (let [rg (r/rng :mersenne 3)] (mapv n (take 10 (r/sequence-generator g d rg))))])
         (for [g [:r2 :halton :sobol :default :ball] d [1 2 6] j [0.0 0.25 1.0]]
           [[:jittered g d j] (let [rg (r/rng :mersenne 3)] (mapv n (take 20 (r/jittered-sequence-generator g d j rg))))]))))

(def seq-methods
  (into {}
        (for [mth [:uniform :systematic :stratified :antithetic :r2 :sobol :halton]]
          [mth (vec (r/->seq (r/rng :mersenne 5) 10 mth))])))

(def dist-samples
  (let [s (fn [d k] (vec (r/->seq d k)))]
    {:inverse-gaussian-1-50 (s (r/inverse-gaussian {:mu 1.0 :lambda 50.0 :rng (r/rng :mersenne 9)}) 10)
     :inverse-gaussian-2-100 (s (r/inverse-gaussian {:mu 2.0 :lambda 100.0 :rng (r/rng :mersenne 9)}) 10)
     :dirichlet-234 (s (r/dirichlet {:alpha [2.0 3.0 4.0] :rng (r/rng :mersenne 9)}) 10)
     :truncated-n-3-4 (s (r/truncated {:distr (r/normal {:rng (r/rng :mersenne 3)}) :left 3.0 :right 4.0}) 10)
     :truncated-exp (s (r/truncated {:distr (r/exponential {:mean 2.0 :rng (r/rng :mersenne 3)}) :left 0.5 :right 3.0}) 10)
     :truncated-n-moderate (s (r/truncated {:distr (r/normal {:rng (r/rng :mersenne 3)}) :left -1.0 :right 2.0}) 10)
     :mixture-2 (s (r/mixture {:distrs [(r/normal {:rng (r/rng :mersenne 1)}) (r/normal {:mu 5.0 :rng (r/rng :mersenne 2)})]
                               :weights [0.3 0.7] :rng (r/rng :mersenne 3)}) 10)
     :logarithmic-5 (s (r/logarithmic {:p 0.5 :rng (r/rng :mersenne 9)}) 20)
     :zip-default (s (r/zero-inflated-poisson {:rng (r/rng :mersenne 9)}) 20)
     :nbi (s (r/nbi {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20)
     :nbii (s (r/nbii {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20)
     :enumerated-real-sorted (s (r/enumerated-real {:data [1.0 2.0 3.0 5.0] :probabilities [0.1 0.4 0.3 0.2] :rng (r/rng :mersenne 9)}) 20)
     :enumerated-real-unsorted (s (r/enumerated-real {:data [3.0 1.0 2.0 5.0] :probabilities [0.1 0.4 0.3 0.2] :rng (r/rng :mersenne 9)}) 20)}))

(def interior-xs [-3.5 -1.2 -0.3 0.4 0.9 1.7 2.6 4.1 7.3 12.5])

(def dist-values
  (let [f (fn [d xs] {:xs xs :pdf (mapv #(r/pdf d %) xs) :cdf (mapv #(r/cdf d %) xs)})
        pos [0.2 0.5 0.9 1.7 2.6 4.1 7.3 12.5]]
    {:logistic (assoc (f (r/logistic {:mu 0.5 :s 2.0}) interior-xs) :mean (r/mean (r/logistic {:mu 0.5 :s 2.0})))
     :logistic-default (assoc (f (r/logistic) interior-xs) :variance (r/variance (r/logistic)))
     :gumbel (assoc (f (r/gumbel {:mu 1.0 :beta 2.0}) interior-xs) :variance (r/variance (r/gumbel {:mu 1.0 :beta 2.0})))
     :nakagami (f (r/nakagami {:mu 1.5 :omega 2.0}) pos)
     :f (f (r/f {:numerator-degrees-of-freedom 5.0 :denominator-degrees-of-freedom 7.0}) pos)
     :levy (f (r/levy {:mu 0.1 :c 1.5}) [0.2 0.5 0.9 1.7 2.6 4.1 7.3 12.5])
     :inverse-gaussian-1-50 (f (r/inverse-gaussian {:mu 1.0 :lambda 50.0}) [0.5 0.8 1.0 1.2 1.5 2.0])
     :truncated-n-moderate (let [d (r/truncated {:distr (r/normal) :left -1.0 :right 2.0})]
                             (assoc (f d [-0.9 -0.3 0.0 0.5 1.5 1.99]) :mean (r/mean d) :variance (r/variance d)
                                    :icdf (mapv #(r/icdf d %) [0.01 0.25 0.5 0.75 0.99])))
     :mixture (let [d (r/mixture {:distrs [(r/normal) (r/normal {:mu 5.0})] :weights [0.3 0.7]})]
                (assoc (f d [-1.0 0.0 1.0 3.0 5.0 6.0]) :mean (r/mean d) :variance (r/variance d)
                       :icdf (mapv #(r/icdf d %) [0.1 0.25 0.9])))
     :multi-normal (let [d (r/multi-normal {:means [0.0 1.0] :covariances [[1.0 0.5] [0.5 2.0]]})]
                     {:points [[0.0 0.0] [1.0 -1.0] [2.5 3.0]] :pdf (mapv #(r/pdf d %) [[0.0 0.0] [1.0 -1.0] [2.5 3.0]])
                      :lpdf (mapv #(r/lpdf d %) [[0.0 0.0] [1.0 -1.0] [2.5 3.0]])})
     :logarithmic (let [d (r/logarithmic {:p 0.5})] {:ks [1 2 3 5 10 40] :pdf (mapv #(r/pdf d %) [1 2 3 5 10 40])
                                                      :cdf (mapv #(r/cdf d %) [1 2 3 5 10 40])
                                                      :icdf (mapv #(r/icdf d %) [0.1 0.5 0.9 0.99])})
     :zip (let [d (r/zero-inflated-poisson)] {:ks [0 1 2 5 9] :pdf (mapv #(r/pdf d %) [0 1 2 5 9])
                                              :cdf (mapv #(r/cdf d %) [0 1 2 5 9])
                                              :icdf (mapv #(r/icdf d %) (map #(/ % 20.0) (range 1 20)))})
     :dirichlet (let [d (r/dirichlet {:alpha [2.0 3.0 4.0]})]
                  {:pdf (mapv #(r/pdf d %) [[0.2 0.3 0.5] [0.1 0.1 0.8]]) :lpdf (mapv #(r/lpdf d %) [[0.2 0.3 0.5] [0.1 0.1 0.8]])})
     :multinomial (let [d (r/multinomial {:n 10 :ps [1 2 3]})]
                    {:pdf (mapv #(r/pdf d %) [[1 3 6] [2 2 6] [0 5 5]]) :cdf (mapv #(r/cdf d %) [[1 3 6] [5 5 5]])})
     :neg-binomial (let [d (r/negative-binomial {:r 5.0 :p 0.4})] {:pdf (mapv #(r/pdf d %) [0 1 5 20]) :cdf (mapv #(r/cdf d %) [0 1 5 20])})
     :nbi (let [d (r/nbi {:mu 4.0 :sigma 0.5})] {:pdf (mapv #(r/pdf d %) [0 1 5 20]) :cdf (mapv #(r/cdf d %) [0 1 5 20])})
     :beta-binomial (let [d (r/beta-binomial {:alpha 2.0 :beta 3.0 :n 10})] {:pdf (mapv #(r/pdf d %) [0 1 5 10]) :cdf (mapv #(r/cdf d %) [0 1 5 10])})
     :integer-discrete (let [d (r/integer-discrete {:data [1 2 4] :probabilities [0.2 0.3 0.5]})] {:pdf (mapv #(r/pdf d %) [0 1 2 3 4 5]) :cdf (mapv #(r/cdf d %) [0 1 2 3 4 5])})
     :ccdf-normal (let [d (r/normal)] {:xs [0.0 1.0 3.0 5.0] :ccdf (mapv #(r/ccdf d %) [0.0 1.0 3.0 5.0])})}))

(def baseline
  {:commit "34b05ef"
   :rng-streams rng-streams
   :series series
   :noise noise-values
   :sequences sequences
   :seq-methods seq-methods
   :dist-samples dist-samples
   :dist-values dist-values})

(binding [*print-length* nil *print-level* nil]
  (spit "test/resources/random/baseline_reference.edn"
        (with-out-str (pp/write baseline :right-margin 160 :pretty true))))
(println "written")
(shutdown-agents)
