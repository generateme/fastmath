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
      (t/is (= [false true false false] (sut/binary-process-list [1 0 1 1] zero?))))
    (t/testing "Non-numeric labels, no true-value: boolean-coerced, not returned unchanged"
      (t/is (= [true false false true] (sut/binary-process-list [:a nil false "x"])))
      (t/is (= [true true true] (sut/binary-process-list ["yes" "no" "yes"])))))
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

;; Group 1: counts & aggregates -- manual identities, both map and 4-arg arities

(t/deftest counts-and-aggregates
  (let [cm {:tp 5.0 :fp 3.0 :fn 2.0 :tn 10.0}]
    (t/are [f v] (and (m/delta-eq v (f cm)) (m/delta-eq v (apply f ((juxt :tp :fp :fn :tn) cm))))
      sut/tp    5.0
      sut/fp    3.0
      sut/fn    2.0
      sut/tn    10.0
      sut/p     7.0   ;; tp+fn
      sut/n     13.0  ;; fp+tn
      sut/pp    8.0   ;; tp+fp
      sut/pn    12.0  ;; fn+tn
      sut/total 20.0)))

;; Group 2: predictive-value rates -- R metrica::precision/FDR/FOR/npv
;; #+begin_src R
;; obs/pred built from tp/fp/fn/tn counts, pos_level=2 ("Pos")
;; precision(obs=obs, pred=pred, pos_level=2); FDR(...); FOR(...); npv(...)
;; #+end_src
;; Matrix A (wikipedia tp=20 fp=180 fn=10 tn=1820): precision=0.1, FDR=0.9, FOR=0.005464481, npv=0.9945355
;; Matrix B (tp=5 fp=6 fn=4 tn=2): precision=0.4545455, FDR=0.5454545, FOR=0.6666667, npv=0.3333333

(t/deftest predictive-value-rates
  (let [mA {:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
        mB {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}]
    (t/are [f vA vB] (and (m/delta-eq vA (f mA) 1.0e-6) (m/delta-eq vB (f mB) 1.0e-6))
      sut/precision 0.1        0.4545455
      sut/ppv       0.1        0.4545455
      sut/fdr       0.9        0.5454545
      sut/for       0.005464481 0.6666667
      sut/npv       0.9945355  0.3333333)))

;; Group 3: sensitivity/specificity family -- R metrica::recall/TPR/FNR/FPR/TNR
;; Matrix A: recall=TPR=0.6666667, FNR=0.3333333, FPR=0.09, TNR=0.91
;; Matrix B: recall=TPR=0.5555556, FNR=0.4444444, FPR=0.75, TNR=0.25
;; sensitivity/hit-rate/miss-rate/fall-out/specificity/selectivity are literal aliases
;; (same value returned) of recall/tpr/fnr/fpr/tnr -- checked equal to the canonical fn too.

(t/deftest sensitivity-specificity-family
  (let [mA {:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
        mB {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}]
    (t/are [f vA vB] (and (m/delta-eq vA (f mA) 1.0e-6) (m/delta-eq vB (f mB) 1.0e-6))
      sut/recall      0.6666667 0.5555556
      sut/sensitivity 0.6666667 0.5555556
      sut/tpr         0.6666667 0.5555556
      sut/hit-rate    0.6666667 0.5555556
      sut/fnr         0.3333333 0.4444444
      sut/miss-rate   0.3333333 0.4444444
      sut/fall-out    0.09      0.75
      sut/fpr         0.09      0.75
      sut/specificity 0.91      0.25
      sut/tnr         0.91      0.25
      sut/selectivity 0.91      0.25)))

;; Group 4: aggregate scores -- R metrica::accuracy/balacc/mk/deltap; prevalence/error/
;; informedness/bm have no direct metrica counterpart (metrica's `preval` returns NULL
;; under this call shape, and `TSS` is a regression metric -- Total Sum of Squares, not
;; the True Skill Statistic/informedness despite the name), verified by manual formula
;; (prevalence=p/total, error=1-accuracy, informedness=tpr+tnr-1) instead.
;; Matrix A: accuracy=0.9064039, balacc=0.7883333, mk=deltap=0.09453552
;; Matrix B: accuracy=0.4117647, balacc=0.4027778, mk=deltap=-0.2121212

(t/deftest aggregate-scores
  (let [mA {:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
        mB {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}]
    (t/are [f vA vB] (and (m/delta-eq vA (f mA) 1.0e-6) (m/delta-eq vB (f mB) 1.0e-6))
      sut/accuracy      0.9064039  0.4117647
      sut/ba             0.7883333  0.4027778
      sut/error          0.0935961  0.5882353
      sut/prevalence     0.014778325123152709 0.5294117647058824
      sut/mk             0.09453552 -0.2121212
      sut/markedness      0.09453552 -0.2121212
      sut/deltaP          0.09453552 -0.2121212
      sut/informedness    0.5766667  -0.1944444
      sut/bm              0.5766667  -0.1944444)))

;; Group 5: odds-ratio family -- R metrica::posLr/negLr/dor/p4/khat.
;; mcc/phi cross-checked vs metrica::mcc on Matrix B, and vs Python
;; sklearn.matthews_corrcoef on Matrix A (metrica's own mcc overflows R's 32-bit
;; integers on Matrix A's larger counts -- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) exceeds
;; .Machine$integer.max -- and returns NA there; not usable for that matrix).
;; pt (prevalence threshold) verified via manual re-derivation of the published
;; Balayla formula, NOT metrica::preval_t: reading metrica's own R source showed
;; preval_t computes TN/(TN+FP) (specificity) and mislabels it "FPR", when the
;; formula requires the actual FP/(FP+TN) -- a bug in the metrica package itself,
;; confirmed by source inspection, not a fastmath issue.
;; Matrix A: posLr=7.407407, negLr=0.3663004, dor=20.22222, mcc(sklearn)=0.2334855,
;;           p4=0.2940226, khat=0.1521213, pt(manual)=0.2686976
;; Matrix B: posLr=0.7407407, negLr=1.777778, dor=0.4166667, mcc(metrica)=-0.2030906,
;;           p4=0.3636364, khat=-0.1971831, pt(manual, matches existing fixture)=0.5374428
;; degenerate (tp=4,fp=fn=tn=0): fastmath mcc/phi -> ##NaN, matching metrica::mcc
;; (also NaN there) -- NOT scikit-learn's 0.0-with-warning convention (documented).

(t/deftest odds-ratio-family
  (let [mA {:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
        mB {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}]
    (t/are [f vA vB] (and (m/delta-eq vA (f mA) 1.0e-6) (m/delta-eq vB (f mB) 1.0e-6))
      sut/lr+    7.407407  0.7407407
      sut/lr-    0.3663004 1.777778
      sut/dor    20.2222222 0.4166667
      sut/mcc    0.2334855 -0.2030906
      sut/phi    0.2334855 -0.2030906
      sut/p4     0.2940226 0.3636364
      sut/kappa  0.1521213 -0.1971831
      sut/pt     0.2686976 0.5374428))
  (t/testing "degenerate zero-denominator confusion matrix -> ##NaN (documented)"
    (t/is (m/invalid-double? (sut/mcc 4.0 0.0 0.0 0.0)))
    (t/is (m/invalid-double? (sut/phi 4.0 0.0 0.0 0.0)))))

;; Group 6: F-scores & composite indices.
;; f1-score/->f-beta vs both R metrica::fscore and Python sklearn.fbeta_score (agree).
;; agf vs metrica::agf, fm vs metrica::fmi, gmean vs metrica::gmean.
;; jaccard/ts/csi vs Python sklearn.jaccard_score, NOT metrica::jaccardindex: reading
;; metrica's own R source showed its binary-case formula is TP/(TP+TN+FP), not the
;; standard TP/(TP+FP+FN) -- a bug in the metrica package, confirmed by source
;; inspection and by sklearn independently matching fastmath's existing value instead.
;; fm verified against metrica::fmi, NOT sklearn.fowlkes_mallows_score: sklearn's
;; version is a *different* metric (pair-counting cluster-agreement FM, for comparing
;; two clusterings) that happens to share a name with the classification-FM-index
;; (sqrt(precision*recall)) fastmath/metrica both compute -- confirmed by sklearn
;; giving 0.9085 on Matrix A vs fastmath/metrica's 0.2582, an unrelated number.
;; Matrix A: f1=0.173913, f0.5=0.1204819, f2=0.3125, agf=0.5523798, fmi=0.2581989,
;;           gmean=0.7788881, jaccard(sklearn)=0.09523810
;; Matrix B: f1=0.5, f0.5=0.4716981, f2=0.5319149, agf=0.407705, fmi=0.5025189,
;;           gmean=0.372678, jaccard(sklearn)=0.3333333

(t/deftest f-scores-and-composite-indices
  (let [mA {:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
        mB {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}]
    (t/are [f vA vB] (and (m/delta-eq vA (f mA) 1.0e-6) (m/delta-eq vB (f mB) 1.0e-6))
      sut/f1-score         0.173913  0.5
      (sut/->f-beta 1.0)   0.173913  0.5
      (sut/->f-beta 0.5)   0.1204819 0.4716981
      (sut/->f-beta 2.0)   0.3125    0.5319149
      sut/adj-f-score      0.5523798 0.407705
      sut/agf              0.5523798 0.407705
      sut/fm               0.2581989 0.5025189
      sut/gmean            0.7788881 0.372678
      sut/jaccard          0.09523810 0.3333333
      sut/ts               0.09523810 0.3333333
      sut/csi              0.09523810 0.3333333)
    (t/testing "->f-inv-beta is the algebraic reciprocal of ->f-beta"
      (doseq [beta [1.0 0.5 2.0 -1.0 0.0]]
        (t/is (m/delta-eq (/ 1.0 ((sut/->f-beta beta) mA)) ((sut/->f-inv-beta beta) mA) 1.0e-9))
        (t/is (m/delta-eq (/ 1.0 ((sut/->f-beta beta) mB)) ((sut/->f-inv-beta beta) mB) 1.0e-9))))
    (t/testing "->f-beta: beta sign doesn't matter (only beta^2 used)"
      (t/is (= ((sut/->f-beta 2.0) mA) ((sut/->f-beta -2.0) mA))))))

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
      :adj-f-score   0.40770504566 ;; alias of :agf, previously untested by this key name (Group 7 gap)
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

;; Group 7: bulk orchestration.
;; R precrec (re-confirmed this audit: evalmod(scores, labels, mode="basic") on the
;; same labels/scores reproduces every value below exactly, closing Unknown C --
;; the NaN/Inf boundary-guard idiom (ppv/fdr guarded at the start via rest+prepend,
;; for/npv guarded at the end via butlast+append) matches precrec's own tie-handling
;; convention point-for-point).

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
      :n 2 :p 3 :total 5 :prevalence 0.6)
    (t/testing "previously-untested sequences: :kappa/:ts/:jaccard need no boundary
                guard on this dataset (no 0/0 hit), verified exact against the
                already-verified scalar fns applied per-threshold to res's own tp/fp/fn/tn"
      (t/is (v/delta-eq (map sut/kappa (:tp res) (:fp res) (:fn res) (:tn res)) (:kappa res) 1.0e-9))
      (t/is (v/delta-eq (map sut/jaccard (:tp res) (:fp res) (:fn res) (:tn res)) (:ts res) 1.0e-9))
      (t/is (= (:ts res) (:jaccard res))))
    (t/testing ":fm/:fdr/:for/:npv match the per-threshold scalar fns everywhere except
                the one boundary point each already-established guard (ppv's start-guard,
                for's end-guard) affects -- documented, not a bug (see Group 5/6 notes)"
      (t/is (= (rest (map sut/fm (:tp res) (:fp res) (:fn res) (:tn res))) (rest (:fm res)))
            "fm: all but index 0 (ppv's guarded start cascades into fm)")
      (t/is (m/delta= 0.0 (first (:fm res))) "fm's guarded start value is 0.0, not NaN")
      (t/is (= (rest (map sut/fdr (:tp res) (:fp res) (:fn res) (:tn res))) (rest (:fdr res)))
            "fdr: all but index 0 (ppv's own guarded start)")
      (t/is (= (butlast (map sut/for (:tp res) (:fp res) (:fn res) (:tn res))) (butlast (:for res)))
            "for: all but the last index (for's own guarded end)")
      (t/is (= (butlast (map sut/npv (:tp res) (:fp res) (:fn res) (:tn res))) (butlast (:npv res)))
            "npv: all but the last index (for's guard, npv=1-for)"))))

;; measures map / binary-measures-all-calc stress test across 5 diverse confusion
;; matrices (incl. a degenerate all-equal one), closing Unknown D: binary-measures-all-calc
;; independently reimplements most formulas rather than delegating to the top-level fns
;; (only mcc/agf/kappa/p4 are delegated) -- confirms no divergence beyond the two
;; pre-existing fixtures. Also closes the :adj-f-score gap (previously only :agf tested
;; by name in the fixture table above) via a direct key-identity check here, in addition
;; to the numeric :adj-f-score row added to the fixture table.

(def ^:private stress-matrices
  [{:tp 20.0 :fp 180.0 :fn 10.0 :tn 1820.0}
   {:tp 5.0 :fp 6.0 :fn 4.0 :tn 2.0}
   {:tp 100.0 :fp 1.0 :fn 1.0 :tn 100.0}
   {:tp 1.0 :fp 1.0 :fn 1.0 :tn 1.0}
   {:tp 7.0 :fp 13.0 :fn 23.0 :tn 3.0}])

(t/deftest orchestration-stress-test
  (t/testing "measures map's key set matches binary-measures-all-calc's (minus the
              latter's extra cp/cn/pcp/pcn/f-measure/f-beta convenience keys)"
    (t/is (= (set (keys sut/measures))
             (-> (sut/binary-measures-all-calc (first stress-matrices))
                 keys set (disj :cp :cn :pcp :pcn :f-measure :f-beta)))))
  (t/testing "every measures-map key agrees with binary-measures-all-calc's same key,
              across 5 diverse matrices (incl. degenerate 1/1/1/1)"
    (doseq [mat stress-matrices]
      (let [all (sut/binary-measures-all-calc mat)]
        (doseq [[k f] sut/measures]
          (let [a (f mat), b (get all k)]
            (t/is (or (and (m/invalid-double? a) (m/invalid-double? b))
                      (m/delta-eq a b 1.0e-9))
                  (str k " on " mat)))))))
  (t/testing ":adj-f-score key present and dispatches identically to :agf (Group 7 gap closed)"
    (doseq [mat stress-matrices]
      (t/is (= ((sut/measures :adj-f-score) mat) ((sut/measures :agf) mat)))
      (t/is (= (:adj-f-score (sut/binary-measures-all-calc mat))
               (:agf (sut/binary-measures-all-calc mat)))))))

