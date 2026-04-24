(ns fastmath.stats.binary-test
  (:require [fastmath.stats.binary :as sut]
            [fastmath.stats :as stats]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t]))

(t/deftest binary-confusion
  (t/testing "Confusion keywords"
    (t/is (= :tp (sut/binary-confusion true true)))
    (t/is (= :fn (sut/binary-confusion true false)))
    (t/is (= :tn (sut/binary-confusion false false)))
    (t/is (= :fp (sut/binary-confusion false true))))
  (t/testing "Processing list to true/false values."
    (t/testing "Just numbers"
      (t/is (= [true false true true] (sut/binary-process-list [1 0 1 1])))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 2 3])))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 1 1] nil))))
    (t/testing "Indicate true value"
      (t/is (= [true false true true] (sut/binary-process-list [1 0 1 1] #{1})))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 1 1] [1])))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 1 1] 1)))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 2 3] pos?)))
      (t/is (= [true false true true] (sut/binary-process-list [1 0 2 3] [1 2 3])))
      (t/is (= [false true false false] (sut/binary-process-list [1 0 1 1] #{0})))
      (t/is (= [false true false false] (sut/binary-process-list [1 0 1 1] [0])))
      (t/is (= [false true false false] (sut/binary-process-list [1 0 1 1] 0)))
      (t/is (= [false true false false] (sut/binary-process-list [1 0 1 1] zero?)))))
  (t/testing "Confusion matrix inference"
    (let [res {:tp 0 :fn 1 :fp 2 :tn 3}]
      (t/are [in] (= res (sut/infer-confusion-matrix in))
        {[:f :n] 1 [:f :p] 2 [:t :n] 3}
        {[:t :p] 0 [:f :n] 1 [:f :p] 2 [:t :n] 3}
        {:fn 1 :fp 2 :tn 3}
        res
        [[0 1] [2 3]]
        [0 1 2 3]))
    (t/is (= :a (try (or (sut/infer-confusion-matrix :a) false)
                     (catch clojure.lang.ExceptionInfo e (:input (ex-data e))))))))

;; R metrica

(t/deftest binary-measures
  (t/testing "wikipedia example"
    (let [mat {:tp 20 :fn 10 :fp 180 :tn 1820}
          res (sut/binary-measures-all-calc mat)]
      (t/are [k v] (and (m/delta-eq v (res k) 1.0e-2)
                        (m/delta-eq v ((sut/measures k) mat) 1.0e-2))
        :p             30
        :n             2000
        :pp            200
        :pn            1830
        :total         2030
        :fdr           0.9
        :accuracy      0.9064
        :fnr           0.3333
        :miss-rate     0.3333
        :recall        0.6667
        :ppv           0.1
        :tnr           0.91
        :f1-score      0.174
        :for           0.0055
        :precision     0.1
        :lr+           7.41
        :lr-           0.366
        :prevalence    0.0148
        :sensitivity   0.6667
        :npv           0.9945
        :specificity   0.91
        :dor           20.22
        :fpr           0.09
        :fall-out      0.09
        :selectivity   0.91
        :tpr           0.6667)))
  (let [mat (stats/confusion-matrix
             [1 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 0]
             [0 1 1 0 1 1 0 1 1 0 0 0 1 1 1 1 1])
        res (sut/binary-measures-all-calc mat)]
    (t/are [k v] (and (m/delta-eq v (res k) 1.0e-10)
                      (m/delta-eq v ((sut/measures k) mat) 1.0e-10))
      :p             9.0
      :n             8.0
      :pp            11.0
      :pn            6.0
      :total         17.0
      :kappa        -0.1971830986
      :csi           0.3333333333
      :error         0.5882352941
      :gmean         0.37267799625
      :fdr           0.5454545455
      :p4            0.3636363636
      :accuracy      0.4117647059
      :fnr           0.4444444444
      :miss-rate     0.4444444444
      :mcc          -0.2030905986
      :ba            0.4027777778
      :mk           -0.21212121212
      :markedness   -0.21212121212
      :recall        0.55555555556
      :deltaP       -0.21212121212
      :informedness -0.19444444444
      :tnr           0.25
      :f1-score      0.5
      :for           0.66666666667
      :precision     0.45454545455
      :ppv           0.45454545455
      :lr+           0.74074074074
      :lr-           1.7777777778
      :prevalence    0.52941176471
      :sensitivity   0.55555555556
      :phi          -0.20309059861
      :npv           0.33333333333
      :jaccard       0.33333333333
      :ts            0.33333333333
      :specificity   0.25
      :dor           0.41666666667
      :fpr           0.75
      :fall-out      0.75
      :selectivity   0.25
      :agf           0.40770504566
      :hit-rate      0.55555555556
      :bm           -0.19444444444
      :tpr           0.55555555556
      :fm            0.50251890763
      :pt            0.537442846107)
    (t/testing "fbeta and it's inverse"
      (t/are [beta v] (and (m/delta-eq v ((sut/->f-beta beta) mat) 1.0e-10)
                           (m/delta-eq v ((:f-beta res) beta) 1.0e-10))
        1.0   0.5
        0.5   0.47169811321
        2.0   0.53191489362
        -1.0  0.5
        -0.5  0.47169811321
        -2.0  0.53191489362
        0.0   0.45454545455
        100   0.55554321139
        0.01  0.45455371833)
      (t/are [beta v] (m/delta-eq (m// v) ((sut/->f-inv-beta beta) mat) 1.0e-10)
        1.0   0.5
        0.5   0.47169811321
        2.0   0.53191489362
        -1.0  0.5
        -0.5  0.47169811321
        -2.0  0.53191489362
        0.0   0.45454545455
        100   0.55554321139
        0.01  0.45455371833))))

;; R precrec

(t/deftest binary-measures-thr
  (let [labels [0 1 0 1 1]
        scores [0 1 1 2 3]
        res (sut/binary-measures-thr labels scores)]
    (t/are [mes vs] (v/delta-eq vs (remove m/invalid-double? (res mes)) 1.0e-9)
      :error       [0.6 0.4 0.2 0.2 0.2 0.4]
      :accuracy    [0.4 0.6 0.8 0.8 0.8 0.6]
      :specificity [1.0 1.0 1.0 0.75 0.5 0.0]
      :sensitivity [0 0.3333333333 0.6666666667 0.8333333333 1 1]
      :precision   [1 1 1 0.8333333333 0.75 0.6]
      :mcc         [0.4082482905 0.6666666667 0.5833333333 0.6123724357] ;; filtered NaNs at both ends
      :f1-score    [0 0.5 0.8 0.8333333333 0.8571428571 0.75]
      :thr         [3 2 1 1 0] ;; filterd ##Inf at beginning
      :tn          [2.0 2.0 2.0 1.5 1.0 0.0]
      :tp          [0.0 1.0 2.0 2.5 3.0 3.0]
      :fn          [3.0 2.0 1.0 0.5 0.0 0.0]
      :fp          [0.0 0.0 0.0 0.5 1.0 2.0]
      :tnr         [1.0 1.0 1.0 0.75 0.5 0.0]
      :tpr         [0.0 0.3333333333333333 0.6666666666666666 0.8333333333333333 1.0 1.0]
      :fnr         [1.0 0.6666666666666666 0.3333333333333333 0.16666666666666666 0.0 0.0]
      :fpr         [0.0 0.0 0.0 0.25 0.5 1.0])
    (t/are [mes v] (m/delta-eq v (res mes) 1.0e-9)
      :n 2 :p 3 :total 5 :prevalence 0.6)))

