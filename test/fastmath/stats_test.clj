(ns fastmath.stats-test
  (:require [fastmath.stats :as sut]
            [fastmath.stats.binary :as bm]
            [clojure.test :as t]
            [fastmath.core :as m]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [fastmath.vector :as v]
            [fastmath.random :as r]))

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
    #_#_#_0.8025362 sut/cohens-u2 nil ;; same as R - changed optimization method, current value differs
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
        :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 1, :step 0.0, :avg 1.0 :probability 1.0})
        :frequencies {1.0 1} :intervals '(1.0 1.0)})
(def one2 {:size 1, :step 1.0, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 2])
         :bins-maps '({:min 1.0, :mid 1.5, :max 2.0, :count 2, :step 1.0, :avg 1.5 :probability 1.0})
         :frequencies {1.5 2} :intervals '(1.0 2.0)})
(def two {:size 2, :step 0.5, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 1] [1.5 1])
        :bins-maps '({:min 1.0, :mid 1.25, :max 1.5, :count 1, :step 0.5, :avg 1.0 :probability 0.5}
                     {:min 1.5, :mid 1.75, :max 2.0, :count 1, :step 0.5, :avg 2.0 :probability 0.5})
        :frequencies {1.0 1, 2.0 1} :intervals '(1.0 1.5 2.0)})

(t/deftest histogram-tests
  (t/are [in method res] (= res (sut/histogram in method))
    [1 1 1 1] :sqrt  {:size 1, :step 0.0, :samples 4, :min 1.0, :max 1.0, :bins '([1.0 4])
                      :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 4, :step 0.0, :avg 1.0
                                    :probability 1.0}),
                      :frequencies {1.0 4} :intervals '(1.0 1.0)}
    [1 1] :sturges {:size 1, :step 0.0, :samples 2, :min 1.0, :max 1.0, :bins '([1.0 2])
                    :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 2, :step 0.0, :avg 1.0 :probability 1.0}),
                    :frequencies {1.0 2} :intervals '(1.0 1.0)}
    [1] :rice one
    [1] :doane one
    [1 2] :sqrt one2
    [1 2] :sturges two
    [1 2] :rice two
    [1 2] :doane one2
    [1 2] :scott one2
    [1 2] :freedman-diaconis one2))

;; chatgpt :)
(t/deftest kruskal-wallis-test
  (t/testing "Basic Kruskal-Wallis test with distinct groups"
    (let [result (sut/kruskal-test [[1 2 3 4] [2 3 4 5] [6 7 8 9]])]
      (t/is (map? result))
      (t/is (contains? result :stat))
      (t/is (contains? result :p-value))
      (t/is (contains? result :df))
      (t/is (= (:df result) 2))
      (t/is (= (:k result) 3))
      (t/is (= (:sides result) :right))
      (t/is (< (:p-value result) 0.05) "Should detect significant difference")))

  (t/testing "Kruskal-Wallis test with identical groups (no difference expected)"
    (let [result (sut/kruskal-test [[5 5 5 5] [5 5 5 5] [5 5 5 1]])]
      (t/is (map? result))
      (t/is (contains? result :stat))
      (t/is (contains? result :p-value))
      (t/is (contains? result :df))
      (t/is (= (:df result) 2))
      (t/is (= (:k result) 3))
      (t/is (> (:p-value result) 0.05) "Should not detect significant difference"))))

(t/deftest transformations
  (let [ys [0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]
        ys- [-1 0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]]
    (t/is (v/delta-eq (sut/box-cox-transformation ys 0.65)
                      [-0.9300129 -0.4227523  0.0000000  0.7804770  6.1394250 12.3319046 25.8838035
                       40.9374869 53.8112999]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:alpha 2.0})
                      [0 1.06018885 1.331677747 1.603605588 2.175413313 6.95789528 12.93691302
                       26.30609562 41.27180377 54.10140563]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys 0.65 {:scaled? true})
                      [-1.9941886  -0.9064904   0.0000000   1.6735450  13.1645181  26.4427989
                       55.5015817  87.7805795 115.3853705]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:scaled? true :alpha 2.0})
                      [0 2.578237269 3.238461899 3.899753984 5.290313778 16.9206693 31.46083955
                       63.97290082 100.3674984 131.5673715]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:negative? true})
                      [-3.0769230769230766 -0.9300128623708235 -0.4227522822458623 0.0 0.7804770384569027
                       6.139425022491503 12.33190461964338 25.883803510261142 40.937486879512214
                       53.8112999422251]))
    (t/is (v/delta-eq (sut/yeo-johnson-transformation ys- 0.65)
                      [-1.147497226 0.2308761917 0.5581765097 0.8756433781 1.521329444 6.554235326
                       12.63614713 26.09538393 41.10482174 53.95645486]))
    ;; same as Julia BoxCoxTrans
    (t/is (m/delta-eq (sut/box-cox-infer-lambda ys) 0.01235497))
    ;; same as Python Scipy
    (t/is (m/delta-eq (sut/yeo-johnson-infer-lambda ys-) -0.004272723))))

(t/deftest l-moments
  (let [ys (sort [-1.7728323 , -1.34577139, -0.06550838, -3.06360177,  0.28895007,
                  1.36727242, -2.39761596, -0.63039948,  0.86878218,  0.68354586,
                  2.82116707,  1.04778304, -0.17578999, -2.32261317,  2.91884183,
                  0.59641717, -0.53668897, -1.09227543,  0.51320672, -4.122064  ,
                  0.20198867,  0.70337326, -0.09336329,  0.21434579,  0.06893888,
                  1.69565577,  2.27596985,  0.60637383,  0.27447569,  0.18262615,
                  -0.24239664, -0.47757852, -4.22660875,  0.92345263,  1.58138485,
                  -1.68156114,  0.0093263 , -0.5041279 ,  3.81543149, -0.82702232,
                  -0.43378152, -2.02045542, -1.07154893, -2.3172839 ,  0.96072838,
                  -0.18502813,  0.526105  ,  0.404828  , -0.18582972, -2.65091393,
                  -0.31282541, -0.19424187, -0.19954772, -0.16113933,  0.42566295,
                  1.44705338, -0.09691397,  0.28710832,  1.28102335,  0.09519895,
                  -1.61865923,  1.12958884,  0.06267352,  0.07283967, -0.66745461,
                  -0.48761821,  2.66665185, -0.27096777,  0.28232669,  0.72006836,
                  -1.39375443, -4.38246225,  7.93752809,  2.78756317, -2.2537608 ,
                  -4.89626616,  0.64561325, -2.4040938 ,  0.02983087, -1.28748567,
                  -0.74376858, -0.16892732, -0.03936628, -0.80297845,  0.69381382,
                  -0.22758131,  1.96957425,  0.29550671,  0.30425838,  3.40986819,
                  0.43560846,  1.83062494,  0.47820415, -2.35097606, -0.06590768,
                  -1.17337123,  0.39161468,  0.48678179,  3.5270126 ])]
    (t/testing "First 8 L-moments"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true}))
        0 1.0
        1 -0.01412282
        2 0.94063132
        3 -0.00167452
        4 0.27196273
        5 0.04125875
        6 0.0800777
        7 0.02874381
        8 0.04678928))
    (t/testing "L-moment ratios"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :ratio? true}))
        3 -0.0017802
        4 0.28912787
        5 0.04386283
        6 0.08513187
        7 0.030558
        8 0.04974242))
    (t/testing "First 8 L-moments, trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 1 :t 1}) 1.0e-8)
        0 1.00000000e+00
        1 -1.24483014e-02
        2 4.01201153e-01
        3 -2.04444131e-02,
        4 7.99520964e-02
        5 4.77843216e-03
        6 1.19496896e-02
        7 4.17554371e-03,
        8 -2.80236176e-04))
    (t/testing "First 8 L-moments, left trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 1 :t 0}))
        0 1.0
        1 0.9265085
        2 0.7042176
        3 0.18019214
        4 0.19576343,
        5 0.07280187
        6 0.06347921
        7 0.04316176
        8 0.03563686))
    (t/testing "First 8 L-moments, right trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 0 :t 1}))
        0 1.0
        1 -0.95475414
        2 0.70672938
        3 -0.18242483
        4 0.14418999,
        5 -0.02329137
        6 0.02994477
        7 -0.0103117
        8 0.01700108))
    (t/testing "L-variation"
      (t/is (m/delta-eq -66.603658638321 (sut/l-variation ys) 1.0e-5)))))

(t/deftest expectile
  (let [ys [-1 0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]]
    (t/is (m/delta-eq (sut/expectile ys 0.5) (sut/mean ys)))
    (t/is (m/delta-eq 25.9 (sut/expectile ys 0.25)))
    (t/is (m/delta-eq 97.54428571428572 (sut/expectile ys 0.75)))))

(t/deftest wmode
  (let [ys [4 1 4 2 4 2]
        ys2 [:a :b :a :c :a :c]
        ws [1, 3, 0.5, 1.5, 1, 2]]
    (t/is (== 4 (sut/wmode ys)))
    (t/is (== 4 (sut/mode ys)))
    (t/is (== 2 (sut/wmode ys ws)))
    (t/is (= :a (sut/wmode ys2)))
    (t/is (= :c (sut/wmode ys2 ws)))))

;; comparison with R test
(t/deftest ks-two-samples
  (t/testing "distinct values"
    (let [x (range -100 101)
          y [21.125 22.925 21.525 18.825 18.225 14.425 24.525 19.325 17.925 16.525 17.425 15.325 10.525 14.825 32.525 30.525 34.025 21.625 15.625 13.425 27.425 26.125 15.925 19.825 15.125]]
      (t/testing "exact method"
        (let [res (sut/ks-test-two-samples x y)
              res+ (sut/ks-test-two-samples x y {:sides :right})
              res- (sut/ks-test-two-samples x y {:sides :left})]
          (t/is (m/delta-eq (:p-value res) 6.675949e-07 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:stat res) 0.5522388))
          (t/is (m/delta-eq (:p-value res+) 3.337974e-07 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:stat res+) 0.5522388))
          (t/is (m/delta-eq (:p-value res-) 0.005897501))
          (t/is (m/delta-eq (:stat res-) 0.3283582))))
      (t/testing "approximate method"
        (let [res (sut/ks-test-two-samples x y {:method :approximate})
              res+ (sut/ks-test-two-samples x y {:sides :right :method :approximate})
              res- (sut/ks-test-two-samples x y {:sides :left :method :approximate})]
          (t/is (m/delta-eq (:p-value res) 2.57807e-06 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:d res) 0.5522388))
          (t/is (m/delta-eq (:p-value res+)  1.289035e-06 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:dp res+) 0.5522388))
          (t/is (m/delta-eq (:p-value res-) 0.008274217))
          (t/is (m/delta-eq (:dn res-) 0.3283582))))))
  (t/testing "ties"
    (let [{:keys [p-value stat]} (sut/ks-test-two-samples [1 1 1 1 1 20] [1 1 1 -1])]
      (t/is (m/delta-eq 0.6666667 p-value))
      (t/is (m/delta-eq 0.25 stat)))
    (let [{:keys [p-value stat]} (sut/ks-test-two-samples [1 1 1 -1] [1 1 1 1 1 20])]
      (t/is (m/delta-eq 0.6666667 p-value))
      (t/is (m/delta-eq 0.25 stat)))))

(t/deftest ks-one-sample
  (t/testing "normality"
    (t/are [side distr s pv]
        (let [{:keys [stat p-value]} (sut/ks-test-one-sample (iris :petal-width)
                                                             (r/distribution distr)
                                                             {:sides side})]
          (and (m/delta-eq stat s)
               (m/delta-eq p-value pv 1.0e-6 1.0e-14)))
      :both :normal 0.5686175 2.981394e-07
      :left :normal 0.5686175 1.490697e-07
      :right :normal 0.006209665 0.9929283))
  (t/testing "cauchy"
    (t/are [side distr s pv]
        (let [{:keys [stat p-value]} (sut/ks-test-one-sample (iris :petal-width)
                                                             (r/distribution distr)
                                                             {:sides side})]
          (and (m/delta-eq stat s)
               (m/delta-eq p-value pv 1.0e-6 1.0e-14)))
      :both :cauchy 0.5317255 2.469265e-06
      :left :cauchy 0.5317255  1.234632e-06
      :right :cauchy 0.1211189 0.4858579)))

(t/deftest auc
  (let [labels [0 1 0 1 1]
        scores [0 1 1 2 3]
        m (sut/binary-measures-thr labels scores)]
    (t/is (m/delta-eq 0.9166666667 (sut/auc m) 1.0e-10))
    (t/is (m/delta-eq 0.9513888889 (sut/auc m :recall :precision) 1.0e-10))
    (t/is (m/delta-eq 0.3585624955 (sut/auc m (bm/->f-beta 2) :recall) 1.0e-10))
    (t/is (m/delta-eq 2.25 (sut/auc m :tp [1 1 0.5 0.5 0.5]) 1.0e-10))
    (t/testing "U based auc-roc"
      (t/is (m/delta-eq 0.9166666667 (sut/auc-roc labels scores) 1.0e-10)))))

(t/deftest multiclass-measures
  (t/are [res m] (= res (sut/multiclass-measure [0, 1, 2, 1, 1, 2] [0, 1, 1, 0, 0, 2] m))
    {0 0.5, 1 0.4, 2 0.6666666666666666} {:average nil}
    {0 0.6666666666666666, 1 0.5, 2 0.8333333333333334} {:metric :accuracy :average nil}
    {0 2.6632016632016633, 1 2.1683991683991684, 2 1.1683991683991684} {:metric :f-inv-beta :beta 0.45 :average nil})
  (t/are [res m] (m/delta-eq res (sut/multiclass-measure [0, 1, 2, 1, 1, 2] [0, 1, 1, 0, 0, 2] m) 1.0e-10)
    0.5222222222222223 {}
    0.5 {:average sut/harmean}
    0.5108729549290354 {:average sut/geomean}
    0.5055555555555555 {:weighted? true}
    0.6111111111111112 {:metric :recall}
    0.3777777777777777 {:metric :mk}
    0.25 {:metric :mk :average :micro}
    0.5641764963265006 {:metric (bm/->f-beta 0.45)}
    0.5641764963265006 {:metric :f-beta :beta 0.45}
    2.0 {:metric :f-inv-beta :beta 0.45}
    0.5 {:average :micro}
    1.0 {:metric :tp}))

;; python

(t/deftest multiclass-auc
  (let [cls [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        probs (zipmap [0 1 2] (apply map vector [[9.81571915e-01, 1.84280700e-02, 1.45178280e-08],
                                                 [9.71303935e-01, 2.86960349e-02, 3.02480322e-08],
                                                 [9.85256375e-01, 1.47436122e-02, 1.23637595e-08],
                                                 [9.76035540e-01, 2.39644199e-02, 3.97833127e-08],
                                                 [9.85224691e-01, 1.47752968e-02, 1.20202360e-08],
                                                 [9.70237188e-01, 2.97627377e-02, 7.39922282e-08],
                                                 [9.86761074e-01, 1.32389059e-02, 2.00173574e-08],
                                                 [9.76132059e-01, 2.38679136e-02, 2.77705876e-08],
                                                 [9.79594414e-01, 2.04055553e-02, 3.06806396e-08],
                                                 [9.68730674e-01, 3.12692937e-02, 3.17879211e-08],
                                                 [9.76228194e-01, 2.37717864e-02, 1.93851435e-08],
                                                 [9.75192202e-01, 2.48077541e-02, 4.40211987e-08],
                                                 [9.74194503e-01, 2.58054754e-02, 2.15639157e-08],
                                                 [9.91855741e-01, 8.14425510e-03, 3.90224543e-09],
                                                 [9.88016494e-01, 1.19835030e-02, 2.84388421e-09],
                                                 [9.86664404e-01, 1.33355827e-02, 1.29289141e-08],
                                                 [9.87967825e-01, 1.20321655e-02, 9.26622416e-09],
                                                 [9.81326796e-01, 1.86731838e-02, 1.97939406e-08],
                                                 [9.56133678e-01, 4.38662530e-02, 6.89813101e-08],
                                                 [9.83979353e-01, 1.60206265e-02, 2.06537553e-08],
                                                 [9.46162653e-01, 5.38372604e-02, 8.69361054e-08],
                                                 [9.81572481e-01, 1.84274864e-02, 3.30183499e-08],
                                                 [9.95955798e-01, 4.04420107e-03, 1.31942828e-09],
                                                 [9.51852640e-01, 4.81471223e-02, 2.37432063e-07],
                                                 [9.51599615e-01, 4.84001781e-02, 2.06841197e-07],
                                                 [9.50935312e-01, 4.90646012e-02, 8.71140079e-08],
                                                 [9.69326068e-01, 3.06738459e-02, 8.66087146e-08],
                                                 [9.74634012e-01, 2.53659631e-02, 2.51118716e-08],
                                                 [9.77037142e-01, 2.29628404e-02, 1.75183184e-08],
                                                 [9.70949888e-01, 2.90500534e-02, 5.86690765e-08],
                                                 [9.63897392e-01, 3.61025375e-02, 7.06065611e-08],
                                                 [9.64477230e-01, 3.55227119e-02, 5.78003442e-08],
                                                 [9.88293432e-01, 1.17065609e-02, 7.09086352e-09],
                                                 [9.88961024e-01, 1.10389705e-02, 5.35397576e-09],
                                                 [9.68320259e-01, 3.16796972e-02, 4.33328420e-08],
                                                 [9.84418938e-01, 1.55810539e-02, 8.03169211e-09],
                                                 [9.78616140e-01, 2.13838501e-02, 9.70819839e-09],
                                                 [9.86732169e-01, 1.32678227e-02, 8.55691000e-09],
                                                 [9.85673700e-01, 1.43262842e-02, 1.55572690e-08],
                                                 [9.73798331e-01, 2.62016402e-02, 2.85802792e-08],
                                                 [9.86459575e-01, 1.35404140e-02, 1.14219407e-08],
                                                 [9.61558525e-01, 3.84414086e-02, 6.60696542e-08],
                                                 [9.88894343e-01, 1.11056459e-02, 1.13028599e-08],
                                                 [9.72218536e-01, 2.77813260e-02, 1.37484231e-07],
                                                 [9.60092902e-01, 3.99068746e-02, 2.23467087e-07],
                                                 [9.73508856e-01, 2.64911043e-02, 4.00774441e-08],
                                                 [9.80164438e-01, 1.98355369e-02, 2.54747511e-08],
                                                 [9.83156866e-01, 1.68431138e-02, 2.01948304e-08],
                                                 [9.78350311e-01, 2.16496704e-02, 1.88317688e-08],
                                                 [9.78396413e-01, 2.16035675e-02, 1.93702474e-08],
                                                 [2.13987432e-03, 8.73951959e-01, 1.23908167e-01],
                                                 [5.85216912e-03, 8.59695656e-01, 1.34452175e-01],
                                                 [1.06758790e-03, 7.25097759e-01, 2.73834653e-01],
                                                 [1.54552953e-02, 9.39649064e-01, 4.48956405e-02],
                                                 [2.38989912e-03, 8.15215488e-01, 1.82394613e-01],
                                                 [7.01568413e-03, 8.60085940e-01, 1.32898376e-01],
                                                 [3.79562357e-03, 7.16470249e-01, 2.79734128e-01],
                                                 [1.47908548e-01, 8.49036974e-01, 3.05447806e-03],
                                                 [2.79169950e-03, 8.96639557e-01, 1.00568743e-01],
                                                 [4.16030884e-02, 9.11777354e-01, 4.66195574e-02],
                                                 [5.58421379e-02, 9.37689054e-01, 6.46880852e-03],
                                                 [1.52500349e-02, 8.98756062e-01, 8.59939028e-02],
                                                 [9.11131071e-03, 9.76561002e-01, 1.43276872e-02],
                                                 [3.06586048e-03, 7.79254339e-01, 2.17679800e-01],
                                                 [7.46682507e-02, 9.14880929e-01, 1.04508203e-02],
                                                 [5.31148888e-03, 9.26323353e-01, 6.83651586e-02],
                                                 [8.75944810e-03, 7.74608924e-01, 2.16631628e-01],
                                                 [1.64929247e-02, 9.65142381e-01, 1.83646947e-02],
                                                 [1.81849807e-03, 8.01061537e-01, 1.97119965e-01],
                                                 [2.40187127e-02, 9.59384631e-01, 1.65966560e-02],
                                                 [2.32188990e-03, 4.40086504e-01, 5.57591606e-01],
                                                 [1.68885300e-02, 9.56721430e-01, 2.63900402e-02],
                                                 [7.19118786e-04, 5.96336663e-01, 4.02944218e-01],
                                                 [3.05498792e-03, 8.59895356e-01, 1.37049656e-01],
                                                 [7.10773506e-03, 9.42926499e-01, 4.99657660e-02],
                                                 [5.10304023e-03, 9.20078724e-01, 7.48182358e-02],
                                                 [1.12624118e-03, 8.01558033e-01, 1.97315726e-01],
                                                 [5.82649343e-04, 4.81159124e-01, 5.18258227e-01],
                                                 [5.51186421e-03, 8.13058474e-01, 1.81429662e-01],
                                                 [6.17716388e-02, 9.34849170e-01, 3.37919145e-03],
                                                 [2.92013412e-02, 9.57187455e-01, 1.36112039e-02],
                                                 [3.72110271e-02, 9.55257338e-01, 7.53163510e-03],
                                                 [2.52496079e-02, 9.56412697e-01, 1.83376948e-02],
                                                 [4.52127608e-04, 3.49726620e-01, 6.49821253e-01],
                                                 [1.02808960e-02, 7.50809181e-01, 2.38909923e-01],
                                                 [1.00508595e-02, 7.88600958e-01, 2.01348183e-01],
                                                 [2.27821616e-03, 8.05162838e-01, 1.92558946e-01],
                                                 [2.77362712e-03, 9.13185469e-01, 8.40409035e-02],
                                                 [2.71657101e-02, 9.28408022e-01, 4.44262680e-02],
                                                 [1.99679531e-02, 9.38027271e-01, 4.20047761e-02],
                                                 [8.74954409e-03, 8.97865098e-01, 9.33853579e-02],
                                                 [4.67010088e-03, 8.28230808e-01, 1.67099091e-01],
                                                 [1.76280668e-02, 9.56967313e-01, 2.54046202e-02],
                                                 [1.21786612e-01, 8.75164475e-01, 3.04891349e-03],
                                                 [1.45196867e-02, 9.20414079e-01, 6.50662344e-02],
                                                 [2.00671541e-02, 9.38012634e-01, 4.19202121e-02],
                                                 [1.71686071e-02, 9.25354578e-01, 5.74768153e-02],
                                                 [8.53517335e-03, 9.35086554e-01, 5.63782730e-02],
                                                 [2.43745722e-01, 7.54963149e-01, 1.29112920e-03],
                                                 [1.92301930e-02, 9.35969076e-01, 4.48007311e-02],
                                                 [9.12432594e-07, 3.91398404e-03, 9.96085104e-01],
                                                 [2.44409256e-04, 1.62562032e-01, 8.37193559e-01],
                                                 [2.50456633e-06, 2.55801238e-02, 9.74417372e-01],
                                                 [3.14684752e-05, 8.16928482e-02, 9.18275683e-01],
                                                 [3.76935037e-06, 1.74447068e-02, 9.82551524e-01],
                                                 [5.59246256e-08, 4.64120216e-03, 9.95358742e-01],
                                                 [5.79752413e-03, 5.13669097e-01, 4.80533378e-01],
                                                 [6.27420734e-07, 2.13546452e-02, 9.78644727e-01],
                                                 [5.26779981e-06, 5.33042449e-02, 9.46690487e-01],
                                                 [6.61332577e-07, 5.74260730e-03, 9.94256731e-01],
                                                 [3.04457419e-04, 2.10442697e-01, 7.89252845e-01],
                                                 [7.32666069e-05, 1.37317523e-01, 8.62609210e-01],
                                                 [2.14474813e-05, 6.52667101e-02, 9.34711842e-01],
                                                 [2.30768730e-04, 1.45237851e-01, 8.54531380e-01],
                                                 [6.95555419e-05, 4.34975598e-02, 9.56432885e-01],
                                                 [5.20470570e-05, 5.40241368e-02, 9.45923816e-01],
                                                 [5.60392603e-05, 1.22915359e-01, 8.77028602e-01],
                                                 [8.56887541e-08, 3.56143008e-03, 9.96438484e-01],
                                                 [3.18606600e-09, 1.00100928e-03, 9.98998988e-01],
                                                 [3.91566006e-04, 4.52005093e-01, 5.47603341e-01],
                                                 [5.65329877e-06, 2.38481898e-02, 9.76146157e-01],
                                                 [6.17556478e-04, 1.90422288e-01, 8.08960156e-01],
                                                 [3.17040922e-08, 4.65751504e-03, 9.95342453e-01],
                                                 [5.89152668e-04, 3.93054898e-01, 6.06355950e-01],
                                                 [1.29504371e-05, 3.86041120e-02, 9.61382938e-01],
                                                 [4.90165820e-06, 5.15012509e-02, 9.48493847e-01],
                                                 [1.07899404e-03, 4.56445696e-01, 5.42475310e-01],
                                                 [1.02969477e-03, 3.85315989e-01, 6.13654316e-01],
                                                 [1.07399358e-05, 3.63419125e-02, 9.63647348e-01],
                                                 [1.70281075e-05, 1.42022336e-01, 8.57960636e-01],
                                                 [1.07347467e-06, 2.92100364e-02, 9.70788890e-01],
                                                 [7.11029023e-07, 1.74223056e-02, 9.82576983e-01],
                                                 [7.94942894e-06, 2.72639644e-02, 9.72728086e-01],
                                                 [5.32052441e-04, 4.75589879e-01, 5.23878068e-01],
                                                 [6.30714650e-05, 1.88649594e-01, 8.11287335e-01],
                                                 [3.97353603e-07, 1.17478837e-02, 9.88251719e-01],
                                                 [1.17294484e-05, 1.73241394e-02, 9.82664131e-01],
                                                 [6.81993651e-05, 1.19491575e-01, 8.80440225e-01],
                                                 [1.63122392e-03, 4.40321653e-01, 5.58047123e-01],
                                                 [4.00086159e-05, 9.34804514e-02, 9.06479540e-01],
                                                 [6.36270339e-06, 2.02880236e-02, 9.79705614e-01],
                                                 [1.00474356e-04, 1.20620465e-01, 8.79279060e-01],
                                                 [2.44409256e-04, 1.62562032e-01, 8.37193559e-01],
                                                 [2.06886328e-06, 1.25871870e-02, 9.87410744e-01],
                                                 [3.84717412e-06, 1.21026066e-02, 9.87893546e-01],
                                                 [5.63724539e-05, 8.00990354e-02, 9.19844592e-01],
                                                 [2.28496543e-04, 2.51916686e-01, 7.47854818e-01],
                                                 [1.39414414e-04, 1.57115825e-01, 8.42744760e-01],
                                                 [4.60438469e-05, 3.84191847e-02, 9.61534771e-01],
                                                 [4.78871489e-04, 2.34864482e-01, 7.64656646e-01]]))]
    (t/are [res m] (m/delta-eq res (sut/multiclass-auc cls probs m) 1.0e-10)
      0.9991777777777777 {:average :micro}
      0.9983333333333334 {}
      0.0016666666666666 {:metric :det}
      0.9967842720318594 {:metric :pr}
      0.1778127938399307 {:metric [(bm/->f-beta 2) :recall]})))
