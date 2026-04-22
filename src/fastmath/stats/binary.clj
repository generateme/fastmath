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
  "True positive rate, TPR, recall, sensitivity, hit-rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double _fp ^double fn ^double _tn] (m// tp (m/+ tp fn))))

(defn sensitivity
  "True positive rate, TPR, recall, sensitivity, hit-rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (recall tp fp fn tn)))

(defn tpr
  "True positive rate, TPR, recall, sensitivity, hit-rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (recall tp fp fn tn)))

(defn hit-rate
  "True positive rate, TPR, recall, sensitivity, hit-rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (recall tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (recall tp fp fn tn)))

(defn fnr
  "False negative rate, FNR, miss-rate"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fnr tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (recall tp fp fn tn))))

(defn miss-rate
  "False negative rate, FNR, miss-rate"
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
  "True negative rate, TNR, specificity, selectivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (specificity tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (m/- 1.0 (fall-out tp fp fn tn))))

(defn tnr
  "True negative rate, TNR, specificity, selectivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (specificity tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (specificity tp fp fn tn)))

(defn selectivity
  "True negative rate, TNR, specificity, selectivity"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (specificity tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (specificity tp fp fn tn)))

;;

(defn prevalence
  "Prevalence, p/(p+n)"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (prevalence tp fp fn tn))
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

(defn p4
  "P4 metric, "
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (p4 tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [v4 (m/* 4.0 tn tp)]
     (m// v4 (m/+ v4 (m/* (m/+ tp tn) (m/+ fp fn)))))))

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

(defn adj-f-score
  "Adjusted f-score, agf"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (adj-f-score tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (let [rec (m// tp (m/+ tp fn))
         prec (m// tp (m/+ tp fp))
         npv (m// tn (m/+ fn tn))
         spec (m// tn (m/+ tn fp))]
     (m/sqrt (m/* (m// (m/* 5.0 rec prec)
                       (m/+ (m/* 4.0 prec) rec))
                  (m// (m/* 1.25 npv spec)
                       (m/+ (m/* 0.25 npv) spec)))))))

(defn agf
  "Adjusted f-score, agf"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (adj-f-score tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn] (adj-f-score tp fp fn tn)))

;;

(defn fm
  "Fowlkes-Mallows index, FM"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (fm tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/sqrt (m/* (ppv tp fp fn tn)
                (tpr tp fp fn tn)))))

(defn gmean
  "Geometric mean score, sqrt(specificity * recall)"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (gmean tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double tn]
   (m/sqrt (m/* (specificity tp fp fn tn)
                (recall tp fp fn tn)))))


(defn jaccard
  "Jaccard index, threat score, TS, critical success index CSI"
  (^double [{:keys [^double tp ^double fp ^double fn ^double tn]}] (jaccard tp fp fn tn))
  (^double [^double tp ^double fp ^double fn ^double _tn]
   (m// tp (m/+ tp fp fn))))

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

(def measures
  {:p p :n n :pp pp :pn pn :total total
   :precision precision :ppv ppv
   :fdr fdr :for for :npv npv
   :recall recall :sensitivity sensitivity :hit-rate hit-rate
   :tpr tpr :miss-rate miss-rate :fnr fnr :fall-out fall-out :fpr fpr
   :specificity specificity :tnr tnr :selectivity selectivity
   :prevalence prevalence
   :accuracy accuracy :error error :ba ba
   :lr+ lr+ :lr- lr- :dor dor
   :informedness informedness :bm bm
   :markedness markedness :mk mk :deltaP deltaP
   :mcc mcc :phi phi :f1-score f1-score
   :fm fm :jaccard jaccard :ts ts :csi csi
   :pt pt :kappa kappa
   :gmean gmean :p4 p4
   :adj-f-score adj-f-score :agf agf})

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
  * a function - function is used to map true values, function should return true value (any value) and `false`/`nil` for false.
  * a value - just a single value."
  ([xs] (binary-process-list xs nil))
  ([xs true-value]
   (let [f (if-not true-value
             (if (every? number? xs) m/not-zero? boolean)
             (cond
               (sequential? true-value) (comp boolean (set true-value))
               (ifn? true-value) (comp boolean true-value)
               :else #(= % true-value)))]
     (map f xs))))

(defn infer-confusion-matrix
  "Construct a confusion matrix from the input.

  Returns a map with `:tp`, `:fn`, `:fp` and `:tn` keys.  

  Input can be one of:

  * a map with keys [:t :p], [:t :n], [:f :p] and [:f :n]
  * pair of pairs [[tp fn] [fp tn]]
  * a sequence of values [tp fn fp tn]
  * a confusion matrix"
  [confusion-matrix]
  (cond

    (and (map? confusion-matrix)
         (every? #{[:t :p] [:t :n] [:f :p] [:f :n]} (keys confusion-matrix))
         (every? number? (vals confusion-matrix)))
    (merge {:tp 0 :fn 0 :fp 0 :tn 0}
           (into {} (map (fn [[[a b] v]] [(keyword (str (name a) (name b))) v]) confusion-matrix)))

    (and (map? confusion-matrix)
         (every? #{:tp :tn :fp :fn} (keys confusion-matrix))
         (every? number? (vals confusion-matrix)))
    (merge {:tp 0 :fn 0 :fp 0 :tn 0} confusion-matrix )
    
    (and (sequential? confusion-matrix)
         (m/== 2 (count confusion-matrix))
         (every? sequential? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] (flatten confusion-matrix))

    (and (sequential? confusion-matrix)
         (every? number? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] confusion-matrix)
    
    :else (throw (ex-info "Can't infer confusion matrix from the input." {:input confusion-matrix}))))

(defn binary-measures-all-calc
  [{:keys [^double tp ^double fp ^double fn ^double tn]
    :or {tp 0.0 fp 0.0 fn 0.0 tn 0.0}}]
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
        mk (m/dec (m/+ ppv npv))
        agf (agf tp fp fn tn)]
    {:tp tp :fp fp :fn fn :tn tn
     :cp cp :p cp
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
     :mk mk :markedness mk :deltaP mk
     :gmean (m/sqrt (m/* tnr tpr))
     :p4 (p4 tp fp fn tn)
     :adj-f-score agf :agf agf}))

(defn binary-measures-thr
  "Calculate binary metrics at various score thresholds for given labels (true/false)."
  ([labels scores] (binary-measures-thr labels scores nil))
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
     {:p p :n n :total total :prevalence (m// p total)
      :tp -tp :fp -fp :fn -fn :tn -tn
      :accuracy (map accuracy -tp -fp -fn -tn)
      :error (map error -tp -fp -fn -tn)
      :tpr tpr :recall tpr :sensitivity tpr
      :fnr fnr :miss-rate fnr
      :fpr fpr :fall-out fpr
      :tnr tnr :specificity tnr
      :ppv ppv :precision ppv
      :fdr fdr :for for :npv npv
      :mcc (map mcc -tp -fp -fn -tn)
      :f1-score (map (fn [^double tp ^double fp ^double fn]
                       (let [tp2 (m/* 2.0 tp)]
                         (m// tp2 (m/+ tp2 fp fn)))) -tp -fp -fn)
      :kappa (map kappa -tp -fp -fn -tn)
      :fm (v/sqrt (v/emult ppv tpr))
      :ts ts :jaccard ts
      :thr (reverse thr)})))

