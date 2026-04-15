(ns fastmath.stats.binary
  "Binary measures functions"
  (:refer-clojure :exclude [for])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn p
  "Real positive"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (p tp fp fn tn))
  (^double [^double tp ^double _fp ^double fn ^double _tn] (m/+ tp fn)))

(defn n
  "Real negative"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (n tp fp fn tn))
  (^double [^double _tp ^double fp ^double _fn ^double tn] (m/+ fp tn)))

(defn pp
  "Predicted positive"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (pp tp fp fn tn))
  (^double [^double tp ^double fp ^double _fn ^double _tn] (m/+ tp fp)))

(defn pn
  "Predicted negative"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (pn tp fp fn tn))
  (^double [^double _tp ^double _fp ^double fn ^double tn] (m/+ fn tn)))

(defn total
  "Total"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (total tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/+ tp fp fn tn)))

;;

(defn precision
  "Precision, positive predictive value, PPV"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (precision tp fp fn tn))
  (^double [^double tp ^double fp ^double _fn ^double _tn] (m// tp (m/+ tp fp))))

(defn ppv
  "Precision, positive predictive value, PPV"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (precision tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (precision tp fp fn tn)))

(defn fdr
  "False discovery rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fdr tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (ppv tp fp fn tn))))

(defn for
  "False omission rate, FOR"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (for tp fp fn tn))
  (^double [^double _tp ^double _fp ^double fn ^double tn] (m// fn (m/+ tn fn))))

(defn npv
  "Negative predictive value"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (npv tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (for tp fp fn tn))))

(defn recall
  "True positive rate, TPR, recall, sensitivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double _fp ^double fn ^double _tn] (m// tp (m/+ tp fn))))

(defn sensitivity
  "True positive rate, TPR, recall, sensitivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (recall tp fp fn tn)))

(defn tpr
  "True positive rate, TPR, recall, sensitivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (recall tp fp fn tn)))

(defn fnr
  "False negative rate, FNR"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fnr tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (recall tp fp fn tn))))

(defn fall-out
  "False positive rate, FPR, fall-out"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fall-out tp fp fn tn))
  (^double [^double _tp ^double fp ^double _fn ^double tn] (m// fp (m/+ fp tn))))

(defn fpr
  "False positive rate, FPR, fall-out"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fall-out tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (fall-out tp fp fn tn)))

(defn specificity
  "True negative rate, TNR, specificity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (specificity tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (fall-out tp fp fn tn))))

(defn tnr
  "True negative rate, TNR, specificity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (specificity tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (specificity tp fp fn tn)))

;;

(defn prevalance
  "Prevalence, p/(p+n)"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (prevalance tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [-p (p tp fp fn tn)
         -n (n tp fp fn tn)]
     (m// -p (m/+ -p -n)))))

(defn accuracy
  "Accuracy"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (accuracy tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [num (m/+ tp tn)]
     (m// num (m/+ num fp fn)))))

(defn error
  "Error"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (error tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [num (m/+ fp fn)]
     (m// num (m/+ num tp tn)))))

(defn ba
  "Balanced accuracy"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (ba tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/* 0.5 (m/+ (m// tp (m/+ tp fn))
                 (m// tn (m/+ tn fp))))))

;;

(defn lr+
  "Positive likelihood ratio"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (lr+ tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m// (tpr tp fp fn tn)
        (fpr tp fp fn tn))))

(defn lr-
  "Negative likelihood ratio"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (lr- tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m// (fnr tp fp fn tn)
        (tnr tp fp fn tn))))

(defn dor
  "Diagnostic odds ratio"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (dor tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m// (m/* tp tn)
        (m/* fp fn))))

;;

(defn informedness
  "Bookmarker informedness, BM"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (informedness tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/dec (m/+ (tpr tp fp fn tn)
               (tnr tp fp fn tn)))))

(defn bm
  "Bookmarker informedness, BM"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (informedness tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (informedness tp fp fn tn)))

(defn markedness
  "Markedness, MK, deltaP"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (markedness tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/dec (m/+ (ppv tp fp fn tn)
               (npv tp fp fn tn)))))

(defn mk
  "Markedness, MK, deltaP"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (markedness tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (markedness tp fp fn tn)))

(defn deltaP
  "Markedness, MK, deltaP"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (markedness tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (markedness tp fp fn tn)))

;;

(defn mcc
  "Matthews correlcation coefficient, MCC, phi "
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (mcc tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] 
   (m// (m/- (m/* tp tn)
             (m/* fp fn))
        (m/sqrt (m/* (m/+ tp fp)
                     (m/+ tp fn)
                     (m/+ tn fp)
                     (m/+ tn fn))))))

(defn phi
  "Matthews correlcation coefficient, MCC, phi "
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (mcc tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (mcc tp fp fn tn)))

;;

(defn f1-score
  "Matthews correlcation coefficient, MCC, phi "
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (f1-score tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double _tn]
   (let [tp2 (m/* 2.0 tp)]
     (m// tp2 (m/+ tp2 fp fn)))))

(defn ->f-beta
  "f-beta score creator, returns f-beta measure function for given `beta`"
  [^double beta]
  (let [beta2 (m/* beta beta)]
    (fn f-beta-score
      (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (f-beta-score tp fp fn tn))
      (^double [^double tp ^double fp ^double fn ^double _tn]
       (m// (m/* (m/inc beta2) tp)
            (m/+ (m/* beta2 (m/+ tp fn)) tp fp))))))

;;

(defn fm
  "Fowlkes-Mallows index, FM"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fm tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/sqrt (m/* (ppv tp fp fn tn)
                (tpr tp fp fn tn)))))

(defn jaccard
  "Jaccard index, threat score, TS, critical success index CSI"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (jaccard tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double _tn]
   (m// tp (m/+ tp fn fp))))

(defn ts
  "Jaccard index, threat score, TS, critical success index CSI"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (jaccard tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (jaccard tp fp fn tn)))

(defn csi
  "Jaccard index, threat score, TS, critical success index CSI"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (jaccard tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (jaccard tp fp fn tn)))

;;

(defn pt
  "Prevalence threshold"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (pt tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [-tpr (tpr tp fp fn tn)
         -fpr (fpr tp fp fn tn)]
     (m// (m/- (m/sqrt (m/* -tpr -fpr)) -fpr)
          (m/- -tpr -fpr)))))

(defn kappa
  "Cohen's kappa"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (kappa tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m// (m/* 2.0 (m/- (m/* tp tn) (m/* fp fn)))
        (m/+ (m/* (m/+ tp fp)
                  (m/+ fp tn))
             (m/* (m/+ tp fn)
                  (m/+ fn tn))))))

;;;;;;

(defn binary-confusion
  "Assign category to a pair of true and predicted values.

  Returned category is one of the following:

  * `:tp` - true positive
  * `:fn` - false negative
  * `:fp` - false positive
  * `:tn` - true negative"
  [t p]
  (cond
    (and t p) :tp
    (and t (not p)) :fn
    (and (not t) p) :fp
    :else :tn))

(defn binary-process-list
  "Convert binary labels to a true/false values.

  Optional `true-value` when is:

  * `nil` - if labels are numbers, all non-zero values are treated as true. Otherwise returns labels unchanged.
  * a sequence - all values in the sequence are treated as true.
  * a function - function is used to map true values, function should return true value (any value) and `false`/`nil`"
  [xs true-value]
  (if-not true-value
    (if (every? number? xs) (map m/not-zero? xs) xs)
    (let [f (if (sequential? true-value) (set true-value) true-value)]
      (map f xs))))

(defn infer-confusion-matrix
  [confusion-matrix]
  (cond

    (and (map? confusion-matrix)
         (every? #{[:t :p] [:t :n] [:f :p] [:f :n]} (keys confusion-matrix)))
    (into {} (map (fn [[[a b] v]] [(keyword (str (name a) (name b))) v]) confusion-matrix))

    (and (sequential? confusion-matrix)
         (m/== 2 (count confusion-matrix))
         (every? sequential? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] (flatten confusion-matrix))

    (and (sequential? confusion-matrix)
         (every? number? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] confusion-matrix)
    
    :else confusion-matrix))

(defn binary-measures-all-calc
  [{:keys [^double tp ^double fp ^double fn ^double tn]
    :or {tp 0.0 fp 0.0 fn 0.0 tn 0.0}
    :as details}]
  (let [cp (m/+ tp fn)
        cn (m/+ fp tn)
        total (m/+ cp cn)
        pcp (m/+ tp fp)
        pcn (m/+ fn tn)
        ppv (m// tp pcp)
        npv (m// tn pcn)
        tpr (m// tp cp)
        fpr (m// fp cn)
        tnr (m/- 1.0 fpr)
        fnr (m/- 1.0 tpr)
        lr+ (m// tpr fpr)
        lr- (m// fnr tnr)
        ts (m// tp (m/+ tp fn fp))
        f-beta (clojure.core/fn ^double [^double beta] (let [b2 (m/* beta beta)]
                                                         (m/* (inc b2) (m// (m/* ppv tpr)
                                                                            (m/+ (m/* b2 ppv) tpr)))))
        f1-score (f-beta 1.0)
        mcc (mcc tp fp fn tn)
        bm (m/dec (m/+ tpr tnr))
        mk (m/dec (m/+ ppv npv))]
    (merge details {:cp cp :p cp
                    :cn cn :n cn
                    :pcp pcp :pp pcp
                    :pcn pcn :pn pcn
                    :total total
                    :tpr tpr
                    :recall tpr
                    :sensitivity tpr
                    :hit-rate tpr
                    :fnr fnr
                    :miss-rate fnr
                    :fpr fpr
                    :fall-out fpr
                    :tnr tnr
                    :specificity tnr
                    :selectivity tnr
                    :prevalence (m// cp total)
                    :accuracy (m// (m/+ tp tn) total)
                    :error (m// (m/+ fp fn) total)
                    :ba (m// (m/+ tpr tnr) 2.0)
                    :ppv ppv
                    :precision ppv
                    :fdr (m/- 1.0 ppv)
                    :npv npv
                    :for (m/- 1.0 npv)
                    :lr+ lr+
                    :lr- lr-
                    :dor (m// lr+ lr-)
                    :fm (m/sqrt (m/* ppv tpr))
                    :pt (m// (m/- (m/sqrt (m/* tpr fpr)) fpr)
                             (m/- tpr fpr))
                    :ts ts :jaccard ts :csi ts
                    :f-measure f1-score
                    :f1-score f1-score
                    :f-beta f-beta
                    :mcc mcc :phi mcc
                    :bm bm :informedness bm
                    :kappa (kappa tp fp fn tn)
                    :mk mk :markedness mk :deltaP mk})))

(defn binary-measures-thr
  "Calculate binary metrics at various score thresholds for given labels (true/false)."
  ([labels scores true-value]
   (let [labels (binary-process-list labels true-value)
         [tp fp thr] (->> (map vector labels scores)
                          (group-by second)
                          (sort-by first m/>)
                          (reduce (fn [curr [score lst]]
                                    (let [[^long tpn ^long fpn] (->> (map first lst)
                                                                     (reduce (fn [[^long t ^long f] l]
                                                                               (if l
                                                                                 [(m/inc t) f]
                                                                                 [t (m/inc f)])) [0 0]))
                                          cnt (double (m/+ tpn fpn))
                                          tp-step (m// tpn cnt)
                                          fp-step (m// fpn cnt)]
                                      (->> (range cnt)
                                           (reduce (fn [[[^double tp-last :as tp]
                                                        [^double fp-last :as fp]
                                                        thr] _]
                                                     [(conj tp (m/+ tp-last tp-step))
                                                      (conj fp (m/+ fp-last fp-step))
                                                      (conj thr score)]) curr)))) ['(0.0) '(0.0) '(##Inf)]))
         p (double (first tp))
         n (double (first fp))
         total (m/+ p n)
         -tp (reverse tp)
         -fp (reverse fp)
         -fn (map (fn [^double v] (m/- p v)) -tp)
         -tn (map (fn [^double v] (m/- n v)) -fp)
         tpr (when (m/pos? p) (v/div -tp p))
         fnr (when (m/pos? p) (v/div -fn p))
         fpr (when (m/pos? n) (v/div -fp n))
         tnr (when (m/pos? n) (v/div -tn n))
         ppv (let [tmp (rest (map (fn [^double tp ^double fp] (m// tp (m/+ tp fp))) -tp -fp))]
               (conj tmp (first tmp)))
         fdr (map (fn [^double v] (m/- 1.0 v)) ppv)
         for (let [tmp (butlast (map (fn [^double tn ^double fn] (m// fn (m/+ tn fn))) -tn -fn))]
               (concat tmp [(last tmp)]))
         npv (map (fn [^double v] (m/- 1.0 v)) for)
         ts (map (fn [^double tp ^double fn ^double fp]
                   (m// tp (m/+ tp fn fp))) -tp -fn -fp)] 
     {:p p :n n :total total :tp -tp :fp -fp :fn -fn :tn -tn
      :accuracy (map accuracy -tp -fp -fn -tn)
      :error (map error -tp -fp -fn -tn)
      :tpr tpr :recall tpr :sensitivity tpr
      :fnr fnr :miss-rate fnr
      :fpr fpr :fall-out fpr
      :tnr tnr :specificity tnr
      :ppv ppv :precision ppv
      :fdr fdr :for for :npv npv
      :mcc (map mcc -tp -tn -fp -fn)
      :f1-score (map (fn [^double tp ^double fp ^double fn]
                       (let [tp2 (m/* 2.0 tp)]
                         (m// tp2 (m/+ tp2 fp fn)))) -tp -fp -fn)
      :kappa (map kappa -tp -tn -fp -fn)
      :fm (v/sqrt (v/emult ppv tpr))
      :ts ts :jaccard ts
      :thr (reverse thr)})))

;;

#_(let [labels [1 1 0 0 0 1 1 1 1 0 1 0 1 0 0 0 1 1 1 0 0 0 0 1 0 1 0 0 1 1 0 1 1 1 0 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 0 0 0 1 1 1 1 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 0 1 0 1 0 1 1 0 1 0 0 0 1 0 0 1 0 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 0 0 1 1 1 1 1 0 1 1 0 0 0 0 1 1 0 1 0 1 0 1 1 1 1 1 0 0 0 1 1 0 1 0 0 0 0 1 0 0 1 0 0 0 0 1 1 0 1 1 1 0 1 1 0 1 1 0 1 0 0 0 1 0 0 0 1 0 1 1 0 1 0 1 0]
        scores [0.6125478 0.364271 0.4321361 0.1402911 0.3848959 0.2444155 0.9706413 0.8901728 0.7817814 0.8687518 0.7166806 0.3601688 0.5479834 0.3852405 0.4237394 0.1017 0.6280956 0.74477 0.6577326 0.4901199 0.07236992 0.1727417 0.1057221 0.8900782 0.9455489 0.9846673 0.3601804 0.4486873 0.0148236 0.5435338 0.2923684 0.7015615 0.7154593 0.7149859 0.1206047 0.3196722 0.9117236 0.7573256 0.09098828 0.5294022 0.257403 0.5899093 0.7084121 0.3266729 0.08654628 0.8794599 0.3626936 0.2301572 0.779772 0.8760862 0.353281 0.2120146 0.7032935 0.6890757 0.6270125 0.2409111 0.402802 0.1347941 0.1204734 0.6654447 0.5363395 0.6234946 0.8851797 0.3537774 0.4089399 0.2656861 0.9321598 0.2485005 0.8588767 0.4917356 0.151351 0.6944575 0.4965132 0.1235049 0.4997881 0.3107186 0.9076511 0.3400782 0.195098 0.371937 0.5173086 0.4195601 0.865639 0.0185276 0.539086 0.005422562 0.7727288 0.7038851 0.3482135 0.2776569 0.4586742 0.05904587 0.1332578 0.08368588 0.5319582 0.4296504 0.7178455 0.5370913 0.2124049 0.9308469 0.08304838 0.4686102 0.3933781 0.6633676 0.3495409 0.1943984 0.8444154 0.9594178 0.2113788 0.9434322 0.5981629 0.834804 0.5768362 0.3803965 0.1618743 0.9123258 0.6429336 0.392174 0.122284 0.5868578 0.1806317 0.08599322 0.7005014 0.06041363 0.531464 0.08425479 0.4484847 0.938583 0.5310065 0.7852131 0.905121 0.7484381 0.6052354 0.8429743 0.8359819 0.3642886 0.4925969 0.4881797 0.259279 0.9910964 0.757364 0.2882583 0.7733362 0.040907 0.110241 0.7607261 0.9845992 0.2532711 0.6972353 0.6205011 0.814586 0.3009731 0.3780921 0.01669441 0.6988265 0.6586926 0.470206 0.5014893 0.2391433 0.05099914 0.08845098 0.1070318 0.7465881 0.4801002 0.3365921 0.5795111 0.1185553 0.2331608 0.4611508 0.3705493 0.7701785 0.537336 0.4632275 0.7902402 0.8834314 0.7451107 0.007746305 0.01265352 0.8683312 0.4394 0.5402213 0.5670432 0.0358154 0.8065439 0.2487075 0.6967022 0.08143913 0.3363153 0.1264804 0.6367285 0.03023506 0.2681383 0.9834944 0.7285364 0.7395543 0.5223845 0.8589705 0.383808 0.6069602 0.1383871]
        ]
    )
