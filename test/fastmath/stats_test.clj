(ns fastmath.stats-test
  (:require [fastmath.stats :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]))

(defn transform [spec data]
  (if (sequential? spec)
    (let [[nm f] spec] [nm (f data)])
    [spec (read-string data)]))
(defn zip [header row] (into {} (map transform header row)))
(defn parse-rows  [header rows]  (map (partial zip header) rows))
(defn read-csv
  [fname header]
  (->> (io/resource fname)
       (slurp)
       (csv/read-csv)
       (rest)
       (parse-rows header)))

(defn data->fn
  [data]
  (memoize (fn
             ([] data)
             ([selector] (map selector data))
             ([selector filter-pred] (map selector (filter filter-pred data))))))

(defn by
  ([data f]
   (mapv second (sort-by first (group-by f (data)))))
  ([data f selector]
   (vec (for [group (by data f)]
          (map selector group)))))

(def iris (data->fn (read-csv "iris.csv"
                            [:sepal-length :sepal-width :petal-length :petal-width [:species keyword]])))
(def mtcars (data->fn (read-csv "mtcars.csv"
                              [[:name identity] :mpg :cyl :disp :hp :drat :wt :qsec :vs :am :gear :carb])))


;; basic

(t/deftest basic-test
  (t/are [f res] (m/delta= res (f (mtcars :mpg)))
    sut/minimum 10.4
    sut/maximum 33.9
    sut/sum 642.9))

;; means

(t/deftest means-test
  (t/are [f res] (m/delta= res (f (mtcars :carb)))
    sut/mean 2.8125
    sut/geomean 2.395987
    sut/harmean 2.026385)
  (t/are [power res] (m/delta= res (sut/powmean (mtcars :carb) power))
    -10.0 1.163977
    -1.0 2.026385
    0.0 2.395987
    0.1 2.43644
    0.5 2.60179
    1.0 2.8125
    2.0 3.230712
    10.0 5.69326))

;; deviations

(t/deftest deviations-test
  (t/are [f res] (m/delta= res (f (mtcars :drat)))
    sut/population-variance 0.2769476
    sut/variance 0.2858814
    sut/population-stddev 0.5262581
    sut/stddev 0.5346787
    sut/variation 0.1486638
    sut/mean-absolute-deviation 0.4532422
    sut/median-absolute-deviation 0.475
    sut/sem 0.09451874)
  (let [d (mtcars :drat)]
    (t/is (m/delta= 0.4459375 (sut/mean-absolute-deviation d (sut/median d))))
    (t/is (m/delta= 3.5965625 (sut/mean-absolute-deviation d 0.0)))
    (t/is (m/delta= 0.465 (sut/median-absolute-deviation d (sut/mean d))))
    (t/is (m/delta= 0.46 (sut/median-absolute-deviation d nil :r1)))
    (t/is (m/delta= 3.695 (sut/median-absolute-deviation d 0)))))

;; effect size

(def mpgs (by mtcars :am :mpg))

(t/deftest effect-size-test
  (t/are [res f attr] (m/delta= res (apply f (mpgs 0) (mpgs 1) attr))
    -1.477947 sut/cohens-d nil
    -1.477947 sut/cohens-d [:unbiased]
    -1.526417 sut/cohens-d [:biased]
    -1.411046 sut/cohens-d [:avg]
    -1.440688 sut/cohens-d-corrected nil
    -1.440688 sut/cohens-d-corrected [:unbiased]
    -1.490360 sut/cohens-d-corrected [:biased]
    -1.4779470958015888 sut/hedges-g nil
    -1.4406879253191958 sut/hedges-g-corrected nil
    -1.4406354024018648 sut/hedges-g* nil
    -1.174886 sut/glass-delta nil
    0.7029826 sut/means-ratio nil
    0.7021799 sut/means-ratio-corrected nil
    -0.659919 sut/cliffs-delta nil
    0.1700405 sut/ameasure nil
    0.204878 sut/wmw-odds nil
    0.8025362 sut/cohens-u2 [:r7] ;; same as R
    0.1538462 sut/cohens-u3 nil))

;; effect size - correlations

(t/deftest effect-size-correlation-test
  (t/are [res f] (m/delta= res (f (mtcars :cyl) (mtcars :mpg)))
    -0.852162 sut/pearson-r
    0.72618 sut/r2-determination
    0.72618 sut/eta-sq
    0.7105671 sut/omega-sq
    0.7170527 sut/epsilon-sq)
  (t/are [method res] (m/delta= res (sut/cohens-f2 (mtcars :cyl) (mtcars :mpg) method))
    :eta 2.652034
    :omega 2.455032
    :epsilon 2.534227)
  (t/are [method res] (m/delta= res (sut/cohens-f (mtcars :cyl) (mtcars :mpg) method))
    :eta 1.628507
    :omega 1.566854
    :epsilon 1.5919255)
  (t/is (m/delta= 2.0921305 (sut/cohens-q (sut/pearson-correlation
                                           (iris :petal-width) (iris :petal-length))
                                          (sut/pearson-correlation
                                           (iris :sepal-width) (iris :sepal-length)))))
  (t/is (m/delta= 2.0921305 (sut/cohens-q (iris :petal-width) (iris :petal-length)
                                          (iris :sepal-width) (iris :sepal-length))))
  (t/is (m/delta= -1.9568811 (sut/cohens-q (mtcars :mpg)
                                           (mtcars :cyl) (mtcars :am)))))

;; kruskal effect size

(t/deftest effect-size-kruskal
  (t/is (m/delta= 0.8305211 (sut/rank-epsilon-sq (by mtcars :cyl :mpg))))
  (t/is (m/delta= 0.818833 (sut/rank-eta-sq (by mtcars :cyl :mpg)))))

;; 2x2 contingency

(def c2x2 (sut/contingency-2x2-measures-all 70 2 4 40))
(def c2x2p (sut/contingency-2x2-measures-all 70 40 40 40))
(def c2x2rr (sut/contingency-2x2-measures-all 15 135 100 150))
(def c2x2ri (sut/contingency-2x2-measures-all 75 75 100 150))

(defn seq-delta-eq
  ([a b] (seq-delta-eq a b 1.0e-6))
  ([a b ^double acc] (every? identity (map #(m/delta-eq %1 %2 acc) a b))))

(defn map-delta-eq
  ([keys a b] (map-delta-eq keys a b 1.0e-6))
  ([keys a b acc]
   (let [f (apply juxt keys)]
     (seq-delta-eq (f a) (f b) acc))))

;; source: wiki, https://statpages.info/ctab2x2.html, https://arxiv.org/pdf/2203.09628.pdf 
(t/deftest contingency-2x2-measures-test
  (t/is (map-delta-eq [:chi2 :yates :cochran-mantel-haenszel]
                      (:p-values c2x2p)
                      {:chi2 0.06 :yates 0.083 :cochran-mantel-haenszel 0.061} 1.0e-3))
  (t/are [ks vs] (seq-delta-eq ((juxt :a :b :c :d) (get-in c2x2 ks)) vs 1.0e-3)
    [:expected] [45.931 26.069 28.069 15.931]
    [:proportions :table] [0.603 0.017 0.034 0.345]
    [:proportions :rows]  [0.972 0.028 0.091 0.909]
    [:proportions :cols]  [0.946 0.048 0.054 0.952])
  (t/is (map-delta-eq [:row1 :row2 :col1 :col2]
                      (:marginals (:proportions c2x2))
                      {:row1 0.621 :row2 0.379 :col1 0.638 :col2 0.362} 1.0e-3))
  (t/is (= (:table c2x2) {:a 70 :b 2 :c 4 :d 40}))
  (t/is (= (:marginals c2x2) {:row1 72 :row2 44 :col1 74 :col2 42 :total 116}))
  (t/are [v k] (m/delta= v (c2x2 k))
    350.0 :OR
    (m/log 350.0) :lOR
    116 :n
    10.6944444 :RR
    (m/sqrt (+ (/ 70.0) (/ 2.0) (/ 4.0) (/ 40.0))) :SE)
  (t/are [v k] (m/delta= v (get-in c2x2 [:measures k]))
    0.889172 :cohens-kappa
    0.994302 :yules-q
    0.8897794 :phi
    0.8891367 :scotts-pi
    0.9589041 :F1
    0.9030371 :gwets-ac1
    0.8965517 :holley-guilfords-g
    0.6666667 :mcnemars-chi2
    0.987322 :TCC
    0.664735 :PCC
    0.940078 :PCC-adjusted
    0.9589041 :F1
    0.8985198 :yules-y
    0.9030371 :gwets-ac1
    0.8813131 :youdens-j
    0.803805 :huberts-gamma
    0.9057971 :bangdiwalas-b
    91.838046 :chi2
    88.062057 :yates
    91.046338 :cochran-mantel-haenszel
    (m/sqrt (/ 91.838046 116.0)) :cramers-v)
  (t/are [v k] (m/delta= v (get-in c2x2rr [:risk k]))
    150 :ES
    250 :CS
    0.1 :EER
    0.4 :CER
    0.3 :ARR
    3.333333 :NNT
    0.25 :RR
    0.75 :RRR
    0.75 :PFu
    (- (/ 15.0 150.0)
       (/ 100.0 250.0)) :RD)
  (t/are [v k] (m/delta= v (get-in c2x2ri [:risk k]))
    150 :ES
    250 :CS
    0.5 :EER
    0.4 :CER
    0.1 :ARI
    10.0 :NNH
    1.25 :RR
    0.25 :RRI
    0.2 :AFe
    (- (/ 75.0 150.0)
       (/ 100.0 250.0)) :RD))

;; histogram

(def one {:size 1, :step 0.0, :samples 1, :min 1.0, :max 1.0, :bins '([1.0 1])
        :bins-maps '({:min 1.0, :max 1.0, :count 1, :step 0.0, :avg 1.0}), :frequencies {1.0 1}})
(def one2 {:size 1, :step 1.0, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 2])
         :bins-maps '({:min 1.0, :max 2.0, :count 2, :step 1.0, :avg 1.5}), :frequencies {1.5 2}})
(def two {:size 2, :step 0.5, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 1] [1.5 1])
        :bins-maps '({:min 1.0, :max 1.5, :count 1, :step 0.5, :avg 1.0}
                     {:min 1.5, :max 2.0, :count 1, :step 0.5, :avg 2.0}), :frequencies {1.0 1, 2.0 1}})

(t/deftest histogram-tests
  (t/are [in method res] (= res (sut/histogram in method))
    [1 1 1 1] :sqrt  {:size 1, :step 0.0, :samples 4, :min 1.0, :max 1.0, :bins '([1.0 4])
                      :bins-maps '({:min 1.0, :max 1.0, :count 4, :step 0.0, :avg 1.0}),
                      :frequencies {1.0 4}}
    [1 1] :sturges {:size 1, :step 0.0, :samples 2, :min 1.0, :max 1.0, :bins '([1.0 2])
                    :bins-maps '({:min 1.0, :max 1.0, :count 2, :step 0.0, :avg 1.0}),
                    :frequencies {1.0 2}}
    [1] :rice one
    [1] :doane one
    [1 2] :sqrt one2
    [1 2] :sturges two
    [1 2] :rice two
    [1 2] :doane one2
    [1 2] :scott one2
    [1 2] :freedman-diaconis one2))
