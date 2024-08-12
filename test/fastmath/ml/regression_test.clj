(ns fastmath.ml.regression-test
  (:require [fastmath.ml.regression :as sut]
            [clojure.test :as t]
            [clojisr.v1.r :as rr]
            [tablecloth.api :as tc]
            [tablecloth.column.api :as tcc]
            [fastmath.core :as m]
            [clojure.pprint]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [clojisr.v1.impl.protocols]))

;; to move to clojisr
(extend nil
  clojisr.v1.impl.protocols/Clojable
  {:->clj (fn [_] nil)})

(defn one-hot
  [dataset column values]
  (reduce (fn [ds v]
            (tc/map-columns ds (keyword (str (name column) (name v))) column
                            (fn [k] (if (= v k) 1 0)))) dataset values))

(defn mult-columns
  ([dataset column1 column2] (mult-columns dataset column1 false column2 false))
  ([dataset column1 sort1? column2 sort2?]
   (let [pre-cs1 (distinct (dataset column1))
         pre-cs2 (distinct (dataset column2))
         cs1 (rest (if sort1? (sort pre-cs1) pre-cs1))
         cs2 (rest (if sort2? (sort pre-cs2) pre-cs2))
         ds (-> dataset
                (one-hot column1 cs1)
                (one-hot column2 cs2))]
     (->> (for [c2 cs2 c1 cs1
                :let [n1 (keyword (str (name column1) (name c1)))
                      n2 (keyword (str (name column2) (name c2)))
                      n (keyword (str (name n1) (name n2)))]]
            [n (tcc/* (ds n1) (ds n2))])
          (reduce (fn [d [n c]]
                    (tc/add-column d n c)) ds)))))

(defn add-poly
  [ds column n]
  (let [poly (take n (v/orthonormal-polynomials (seq (ds column))))]
    (reduce (fn [d [i p]]
              (tc/add-column d (keyword (str (name column) "-poly-" i)) p))
            ds (map-indexed vector poly))))

(defn load-r []
  (rr/require-r '[GLMsData] '[base] '[utils] '[moments] '[stats] '[car] '[MASS] '[statmod]))

(load-r)

(defn init-r []
  (rr/discard-all-sessions)
  (load-r)
  (base/options :contrasts ["contr.treatment" "contr.treatment"])
  (utils/data 'gestation)
  (utils/data 'lungcap)
  (rr/r "lungcap$Smokef <- factor(lungcap$Smoke,levels=c(0, 1), labels=c(\"Non-smoker\",\"Smoker\"))")
  (utils/data 'dental)
  (utils/data 'cheese)
  (utils/data 'turbines)
  (utils/data 'germ)
  (utils/data 'mammary)
  (utils/data 'nminer)
  (utils/data 'deposit)
  (utils/data 'danishlc)
  (rr/r '(<- ($ danishlc Rate) (* 1000 (/ ($ danishlc Cases) ($ danishlc Pop)))))
  (rr/r '(<- ($ danishlc Age) (ordered ($ danishlc Age) :levels ["40-54", "55-59", "60-64", "65-69", "70-74", ">74"])))
  (rr/r '(<- ($ danishlc City) (abbreviate ($ danishlc City) 1)))
  (rr/r '(<- ($ danishlc AgeNum) (rep [40, 55, 60, 65, 70, 75] 4)))
  (utils/data 'kstones)
  (utils/data 'pock)
  (utils/data 'hcrabs)
  (utils/data 'lime)
  (utils/data 'perm)
  (rr/r '(<- ($ perm Day) (factor ($ perm Day))))
  (utils/data 'yieldden)
  (rr/r '(<- ($ yieldden Var) (factor ($ yieldden Var))))
  (rr/r '(<- ($ yieldden YD) (with yieldden (* Yield Dens)))))

(defn with-r-session [f]
  (init-r)
  (f)
  (rr/discard-all-sessions))

(t/use-fixtures :once with-r-session)

(def ys [1 2.1 3.2 4 4.9 5.5 10 7 8.1 8.7])
(def nys (v/normalize ys))
(def xss-1 '(-0.5 -1 -2 -3 -4 -5 -6 -7 -8 -9))
(def xss-2 '(1 2 0 2 1 0 1 2 0 2))
(def xss-3 (v/emult xss-1 xss-2))
(def xss (map vector xss-1 xss-2 xss-3))
(def weights [2 1 2 1 2 1 2 1 2 1])
(def offset (range 10))

(defn lm-tests
  ([ys xss lm-r] (lm-tests ys xss nil lm-r))
  ([ys xss options lm-r]
   (let [lm (sut/lm ys xss options)
         rlm (rr/r lm-r)
         rlmdata (rr/r->clj rlm)
         analysis (sut/analysis lm)
         summary (rr/r->clj `(summary ~rlm :correlation true))
         influence (rr/r->clj `(influence ~rlm))
         lobs (m/log (:observations lm))]
     (t/is (m/delta-eq (:intercept lm) (if (:intercept? lm) (first (:coefficients rlmdata)) 0.0)))
     (t/is (v/delta-eq (:beta lm) (if (:intercept? lm)
                                    (rest (:coefficients rlmdata))
                                    (:coefficients rlmdata))))
     (t/is (v/delta-eq (get-in lm [:residuals :raw] [##Inf]) (:residuals rlmdata)))
     (t/is (v/delta-eq (get-in lm [:residuals :weighted] [##Inf]) (:residuals summary)))
     (t/is (v/delta-eq (get-in lm [:fitted] [##Inf]) (:fitted.values rlmdata)))
     (t/is (== (get-in lm [:df :residual] [##Inf]) (first (:df.residual rlmdata))))
     (t/is (== (get-in lm [:df :model] [##Inf]) (second (:fstatistic summary))))
     (t/is (m/delta-eq (:sigma lm) (first (:sigma summary))))
     (t/is (m/delta-eq (:sigma2 lm) (m/sq (first (:sigma summary)))))
     (t/is (m/delta-eq (:r-squared lm) (first (:r.squared summary))))
     (t/is (m/delta-eq (:adjusted-r-squared lm) (first (:adj.r.squared summary))))
     (t/is (m/delta-eq (:f-statistic lm) (first (:fstatistic summary))))
     (t/is (v/delta-eq (mat/mat->seq (:xtxinv lm)) (:cov.unscalesd summary)))
     (t/is (v/delta-eq (concat (map :estimate (:coefficients lm))
                               (map :stderr (:coefficients lm))
                               (map :t-value (:coefficients lm))
                               (map :p-value (:coefficients lm)))
                       (:coefficients summary)))
     (t/is (m/delta-eq (get-in lm [:ll :log-likelihood] [##Inf]) (-> rlm stats/logLik rr/r->clj first)))
     (t/is (m/delta-eq (get-in lm [:ll :aic-rss] [##Inf]) (-> rlm stats/extractAIC rr/r->clj second)))
     (t/is (m/delta-eq (get-in lm [:ll :aic] [##Inf]) (-> rlm stats/AIC rr/r->clj first)))
     (t/is (m/delta-eq (get-in lm [:ll :bic-rss] [##Inf])
                       (-> (stats/extractAIC rlm :k lobs) rr/r->clj second)))
     (t/is (m/delta-eq (get-in lm [:ll :bic] [##Inf]) (-> (stats/AIC rlm :k lobs) rr/r->clj first)))
     (t/is (v/delta-eq (flatten (:correlation analysis)) (:correlation summary)))
     (t/is (m/delta-eq (get-in analysis [:normality :skewness] [##Inf])
                       (-> rlm stats/weighted-residuals moments/skewness rr/r->clj first)))
     (t/is (m/delta-eq (get-in analysis [:normality :kurtosis] [##Inf])
                       (-> rlm stats/weighted-residuals moments/kurtosis rr/r->clj first)))
     (t/is (m/delta-eq (get-in analysis [:normality :durbin-watson :raw] [##Inf])
                       (-> rlm car/durbinWatsonTest rr/r->clj :dw first)))
     (t/is (v/delta-eq (get-in analysis [:residuals :standardized] [##Inf])
                       (-> rlm stats/rstandard rr/r->clj)))
     (t/is (v/delta-eq (get-in analysis [:residuals :studentized] [##Inf])
                       (-> rlm stats/rstudent rr/r->clj)))
     (t/is (v/delta-eq (get-in analysis [:laverage :hat] [##Inf]) (:hat influence)))
     (t/is (v/delta-eq (get-in analysis [:laverage :sigmas] [##Inf]) (:sigma influence)))
     (t/is (v/delta-eq (flatten (get-in analysis [:laverage :coefficients] [##Inf]))
                       (:coefficients influence)))
     (t/is (v/delta-eq (get-in analysis [:influence :cooks-distance] [##Inf])
                       (-> rlm stats/cooks-distance rr/r->clj)))
     (t/is (v/delta-eq (get-in analysis [:influence :dffits] [##Inf])
                       (-> rlm stats/dffits rr/r->clj)))
     (t/is (v/delta-eq (get-in analysis [:influence :covratio] [##Inf])
                       (-> rlm stats/covratio rr/r->clj)))
     (t/is (v/delta-eq (flatten (get-in analysis [:influence :dfbetas] [##Inf]))
                       (-> rlm stats/dfbetas rr/r->clj vals rest flatten))))))

(t/deftest basic-weights-intercept
  (lm-tests ys xss {:weights weights}
            `(lm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights)))

(t/deftest basic-weights-intercept-transformer
  (lm-tests ys xss {:weights weights :transformer (fn [[a b c]] [(m/exp a) a b c])}
            `(lm (formula ~ys (+ (exp ~xss-1) (* ~xss-1 ~xss-2))) :weights ~weights)))

(t/deftest basic-weights-intercept-offset
  (lm-tests ys xss {:weights weights :offset offset}
            `(lm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights :offset ~offset)))

(t/deftest basic-weights-no-intercept
  (lm-tests ys xss {:weights weights :intercept? false}
            `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights)))

(t/deftest basic-weights-no-intercept-offset
  (lm-tests ys xss {:weights weights :intercept? false :offset offset}
            `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :offset ~offset)))

(t/deftest basic-no-weights-intercept
  (lm-tests ys xss `(lm (formula ~ys (* ~xss-1 ~xss-2)))))

(t/deftest basic-no-weights-intercept-offset
  (lm-tests ys xss {:offset offset} `(lm (formula ~ys (* ~xss-1 ~xss-2)) :offset ~offset)))

(t/deftest basic-no-weights-no-intercept
  (lm-tests ys xss {:intercept? false} `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))))))

(t/deftest basic-no-weights-no-intercept-offset
  (lm-tests ys xss {:intercept? false :offset offset}
            `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :offset ~offset)))

(t/deftest gestation-data
  (let [ds (rr/r->clj 'gestation)]
    (lm-tests (:Weight ds) (:Age ds) {:weights (:Births ds)}
              '(lm (formula Weight Age) :weights Births :data gestation))
    (lm-tests (:Weight ds) (:Age ds)
              '(lm (formula Weight Age) :data gestation))))

(t/deftest lungcap-data
  (let [ds (-> (rr/r->clj 'lungcap)
               (tc/map-columns :GenderM [:Gender] (fn [v] (if (= v :F) 0 1)))
               (tc/map-columns :GenderF [:Gender] (fn [v] (if (= v :F) 1 0)))
               (tc/map-columns :Gender (fn [v] (if (= v :F) 0 1)))
               (tc/log :logFEV :FEV)
               (tc/* :HtSmoke [:Ht :Smoke]))]
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Age :Ht :Gender :Smoke]))
              '(lm (formula (log FEV) (+ Age Ht Gender Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Age :Ht :GenderF :GenderM :Smoke]))
              {:intercept? false}
              '(lm (formula (log FEV) (+ 0 Age Ht Gender Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Ht :Gender :Smoke]))
              '(lm (formula (log FEV) (+ Ht Gender Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Gender :Smoke]))
              '(lm (formula (log FEV) (+ Gender Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Age :Smoke]))
              '(lm (formula (log FEV) (+ Age Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Ht :Smoke]))
              '(lm (formula (log FEV) (+ Ht Smokef)) :data lungcap))
    (lm-tests (:logFEV ds) (tc/rows (tc/select-columns ds [:Ht :Smoke :HtSmoke]))
              '(lm (formula (log FEV) (* Ht Smokef)) :data lungcap))))

(t/deftest dental-data
  (let [ds (-> (rr/r->clj 'dental)
               (tc/map-columns :Indus (fn [v] (if (= v :Ind) 0 1)))
               (tc/* :SugarIndus [:Sugar :Indus]))]
    (lm-tests (:DMFT ds) (tc/rows (tc/select-columns ds [:Sugar :Indus :SugarIndus]))
              `(lm (formula DMFT (* Sugar Indus)) :data dental))))

(t/deftest cheese-data
  (let [ds (-> (rr/r->clj 'cheese)
               (tc/log :lH2S :H2S)
               (tc/* :Acetic-lH2S [:Acetic :lH2S])
               (tc/* :Acetic-Lactic [:Acetic :Lactic])
               (tc/* :lH2S-Lactic [:lH2S :Lactic])
               (tc/* :Acetic-lH2S-Lactic [:Acetic :lH2S :Lactic]))]
    (lm-tests (:Taste ds) (tc/rows (tc/select-columns ds [:Acetic :lH2S :Lactic
                                                          :Acetic-lH2S :Acetic-Lactic
                                                          :lH2S-Lactic :Acetic-lH2S-Lactic]))
              `(lm (formula Taste (* Acetic (log H2S) Lactic)) :data cheese))))

(defn glm-tests
  ([ys xss glm-r] (glm-tests ys xss nil glm-r))
  ([ys xss options glm-r]
   (let [options (merge {:epsilon 1.0e-16 :max-iters 100} options)
         glm (sut/glm ys xss options)
         rglm (rr/r (concat glm-r (list :epsilon (:epsilon options) :maxit (:max-iters options))))
         rglmdata (rr/r->clj rglm)
         analysis (sut/analysis glm)
         summary (rr/r->clj `(summary ~rglm :correlation true))
         influence (rr/r->clj `(influence ~rglm))]
     (when (:print? options) (clojure.pprint/pprint glm))
     (t/is (m/delta-eq (:intercept glm) (if (:intercept? glm) (first (:coefficients rglmdata)) 0.0)))
     (t/is (v/delta-eq (:beta glm) (if (:intercept? glm)
                                     (rest (:coefficients rglmdata))
                                     (:coefficients rglmdata))))
     (t/is (v/delta-eq (get-in glm [:weights :weights] [##Inf]) (:weights rglmdata)))
     (t/is (v/delta-eq (get-in glm [:weights :initial] [##Inf]) (:prior.weights rglmdata)))
     (t/is (v/delta-eq (get-in glm [:residuals :working] [##Inf]) (:residuals rglmdata)))
     (t/is (v/delta-eq (get-in glm [:residuals :deviance] [##Inf]) (:deviance.resid summary)))
     (t/is (v/delta-eq (get-in glm [:residuals :pearson] [##Inf])
                       (rr/r->clj `(residuals ~rglm :type "pearson"))))
     (t/is (v/delta-eq (get-in glm [:residuals :raw] [##Inf])
                       (rr/r->clj `(residuals ~rglm :type "response")) 1.0e-4))
     (t/is (v/delta-eq (get-in glm [:fitted] [##Inf]) (:fitted.values rglmdata) 1.0e-4))
     (t/is (m/delta-eq (get-in glm [:deviance :residual] [##Inf]) (first (:deviance rglmdata))))
     (when (:compare-deviance-null? options)
       (t/is (m/delta-eq (get-in glm [:deviance :null] [##Inf]) (first (:null.deviance rglmdata)))))
     (t/is (== (get-in glm [:df :residual] [##Inf]) (first (:df.residual rglmdata))))
     (t/is (== (get-in glm [:df :null] [##Inf]) (first (:df.null rglmdata))))
     (t/is (m/delta-eq (:dispersion glm) (first (:dispersion summary))))
     (t/is (v/delta-eq (mat/mat->seq (:xtxinv glm)) (:cov.unscalesd summary)))
     (t/is (v/delta-eq (concat (map :estimate (:coefficients glm))
                               (map :stderr (:coefficients glm))
                               (map (if (:estimated-dispersion? glm) :t-value :z-value)
                                    (:coefficients glm))
                               (map :p-value (:coefficients glm)))
                       (:coefficients summary)))
     (when-not (:no-ll? options)
       (t/is (m/delta-eq (get-in glm [:ll :log-likelihood] [##Inf]) (-> rglm stats/logLik rr/r->clj first)))
       (t/is (m/delta-eq (get-in glm [:ll :aic] [##Inf]) (-> rglm stats/AIC rr/r->clj first)))
       (t/is (m/delta-eq (get-in glm [:ll :bic] [##Inf]) (-> rglm stats/BIC rr/r->clj first))))
     (when-not (:no-analysis? options)
       (t/is (v/delta-eq (flatten (:correlation analysis)) (:correlation summary)))
       (t/is (v/delta-eq (get-in analysis [:residuals :studentized] [##Inf])
                         (-> rglm stats/rstudent rr/r->clj)))
       (t/is (v/delta-eq (get-in analysis [:residuals :standardized :deviance] [##Inf])
                         (-> (stats/rstandard rglm :type "deviance") rr/r->clj)))
       (t/is (v/delta-eq (get-in analysis [:residuals :standardized :pearson] [##Inf])
                         (-> (stats/rstandard rglm :type "pearson") rr/r->clj)))
       (t/is (v/delta-eq (get-in analysis [:laverage :hat] [##Inf]) (:hat influence)))
       (t/is (v/delta-eq (get-in analysis [:laverage :sigmas] [##Inf]) (:sigma influence)))
       (t/is (v/delta-eq (flatten (get-in analysis [:laverage :coefficients] [##Inf]))
                         (:coefficients influence)))
       (t/is (v/delta-eq (get-in analysis [:influence :cooks-distance] [##Inf])
                         (-> rglm stats/cooks-distance rr/r->clj)))
       (t/is (v/delta-eq (get-in analysis [:influence :dffits] [##Inf])
                         (-> rglm stats/dffits rr/r->clj)))
       (t/is (v/delta-eq (get-in analysis [:influence :covratio] [##Inf])
                         (-> rglm stats/covratio rr/r->clj)))
       (t/is (v/delta-eq (flatten (get-in analysis [:influence :dfbetas] [##Inf]))
                         (-> rglm stats/dfbetas rr/r->clj vals rest flatten))))
     (when (> (count (:coefficeints glm)) 1)
       (let [dglm (sut/dose glm)
             rdglm (MASS/dose-p rglm)
             dattr (rr/r->clj (base/attributes rdglm))]
         (t/is (m/delta-eq (:dose dglm) (first (rr/r->clj rdglm))))
         (t/is (m/delta-eq (:stderr dglm) (first (:SE dattr))))))
     (when (:qres? options)
       (t/is (v/delta-eq (sut/quantile-residuals glm)
                         (rr/r->clj (statmod/qresid rglm))))))))

(t/deftest dummy-data
  (glm-tests nys xss {:weights weights :family :binomial}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family binomial))
  (glm-tests nys xss {:weights weights :family :binomial :link :probit}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link probit)))
  (glm-tests nys xss {:weights weights :family :binomial :link :cloglog}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link cloglog)))
  (glm-tests nys xss {:weights weights :family :binomial :link :log}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link log)))
  (glm-tests nys xss {:family :binomial}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family binomial))
  (glm-tests nys xss {:family :binomial :link :probit}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link probit)))
  (glm-tests nys xss {:family :binomial :link :cloglog}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link cloglog)))
  (glm-tests nys xss {:family :binomial :link :log}
             `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link log)))
  (glm-tests nys xss {:weights weights :family :binomial :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family binomial))
  (glm-tests nys xss {:weights weights :family :binomial :link :probit :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family (binomial :link probit)))
  (glm-tests nys xss {:weights weights :family :binomial :link :cloglog :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family (binomial :link cloglog)))
  (glm-tests nys xss {:family :binomial :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family binomial))
  (glm-tests nys xss {:family :binomial :link :probit :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family (binomial :link probit)))
  (glm-tests nys xss {:family :binomial :link :cloglog :intercept? false}
             `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family (binomial :link cloglog)))
  (glm-tests nys xss {:family :binomial
                      :transformer (fn [[a b c]] [(m/exp a) a b c])}
             `(glm (formula ~nys (+ (exp ~xss-1) (* ~xss-1 ~xss-2))) :family binomial)))

(t/deftest turbunes-data
  (let [turbines (-> (rr/r->clj 'turbines)
                     (tc/convert-types :Fissures :double)
                     (tc// :ratio [:Fissures :Turbines])
                     (tc/- :diff [:Turbines :Fissures]))]
    (glm-tests (:ratio turbines) (:Hours turbines) {:family :binomial :weights (:Turbines turbines)}
               '(glm (formula (/ Fissures Turbines) Hours)
                     :family binomial :data turbines :weights Turbines))
    (glm-tests (tc/rows (tc/select-columns turbines [:Fissures :diff])) (:Hours turbines) {:family :binomial}
               '(glm (formula (cbind Fissures (- Turbines Fissures)) Hours)
                     :family binomial :data turbines))
    (glm-tests (:ratio turbines) (:Hours turbines)
               {:family :binomial :weights (:Turbines turbines) :link :probit}
               '(glm (formula (/ Fissures Turbines) Hours)
                     :family (binomial :link "probit") :data turbines :weights Turbines))
    (glm-tests (:ratio turbines) (:Hours turbines)
               {:family :binomial :weights (:Turbines turbines) :link :cloglog}
               '(glm (formula (/ Fissures Turbines) Hours)
                     :family (binomial :link "cloglog") :data turbines :weights Turbines))))

(t/deftest germ-data
  (let [germ (-> (rr/r->clj 'germ)
                 (tc/map-columns :SeedsOA75 :Seeds #(if (= :OA75 %) 1 0))
                 (tc/map-columns :ExtractCucumber :Extract #(if (= :Cucumber %) 1 0))
                 (tc/map-columns :ratio [:Germ :Total] m//)
                 (tc/* :ES [:ExtractCucumber :SeedsOA75]))]
    (glm-tests (:ratio germ) (tc/rows (tc/select-columns germ [:SeedsOA75 :ExtractCucumber]))
               {:family :binomial :weights (:Total germ)}
               '(glm (formula (/ Germ Total) (+ Seeds Extract))
                     :family binomial :data germ :weights Total))

    (glm-tests (:ratio germ) (tc/rows (tc/select-columns germ [:SeedsOA75 :ExtractCucumber :ES]))
               {:family :binomial :weights (:Total germ)}
               '(glm (formula (/ Germ Total) (* Seeds Extract))
                     :family binomial :data germ :weights Total))

    (glm-tests (:ratio germ) (tc/rows (tc/select-columns germ [:SeedsOA75 :ExtractCucumber :ES]))
               {:family :quasi-binomial :weights (:Total germ) :no-ll? true}
               '(glm (formula (/ Germ Total) (* Seeds Extract))
                     :family quasibinomial :data germ :weights Total))))

(t/deftest mammary-data
  (let [mammary (-> (rr/r->clj 'mammary)
                    (tc// :ratio [:N.Outgrowths :N.Assays])
                    (tc/log :lncells :N.Cells))]

    ;; only intercept and offset
    (glm-tests (:ratio mammary) :intercept {:family :binomial :link :cloglog
                                            :weights (:N.Assays mammary)
                                            :offset (:lncells mammary)}
               '(glm (formula (/ N.Outgrowths N.Assays) (offset (log N.Cells)))
                     :family (binomial :link cloglog) :weights N.Assays
                     :data mammary))

    ;; intercept, offset, data point
    (glm-tests (:ratio mammary) (:N.Outgrowths mammary) {:family :binomial :link :cloglog
                                                         :weights (:N.Assays mammary)
                                                         :offset (:lncells mammary)}
               '(glm (formula (/ N.Outgrowths N.Assays) N.Outgrowths)
                     :family (binomial :link cloglog) :weights N.Assays
                     :data mammary :offset (log N.Cells)))

    ;; no intercept, offset, data point
    (glm-tests (:ratio mammary) (:N.Outgrowths mammary) {:family :binomial :link :cloglog
                                                         :weights (:N.Assays mammary)
                                                         :offset (:lncells mammary)
                                                         :intercept? false
                                                         :compare-deviance-null? false}
               '(glm (formula (/ N.Outgrowths N.Assays) (+ 0 N.Outgrowths))
                     :family (binomial :link cloglog) :weights N.Assays
                     :data mammary :offset (log N.Cells)))

    (glm-tests (:ratio mammary) (:lncells mammary) {:family :binomial :link :cloglog
                                                    :weights (:N.Assays mammary)}
               '(glm (formula (/ N.Outgrowths N.Assays) (log N.Cells))
                     :family (binomial :link cloglog) :weights N.Assays
                     :data mammary))))

(t/deftest deposit-data
  (let [deposit (-> (rr/r->clj 'deposit)
                    (tc/map-columns :ratio [:Killed :Number] m//)
                    (tc/log :logDep :Deposit)
                    (one-hot :Insecticide [:A :B :C])
                    (add-poly :logDep 2))
        model (sut/glm (:ratio deposit) (tc/rows (tc/select-columns
                                                  deposit [:Deposit :InsecticideB :InsecticideC]))
                       {:family :binomial :weights (:Number deposit)})
        model2 (sut/glm (:ratio deposit) (tc/rows (tc/select-columns
                                                   deposit [:Deposit :InsecticideA :InsecticideB :InsecticideC]))
                        {:family :binomial :weights (:Number deposit) :intercept? false})
        rmodel (rr/r '(glm (formula (/ Killed Number) (+ Deposit Insecticide))
                           :family binomial :weights Number :data deposit))
        newdata (m/slice-range 2.0 8.0 100)]
    (glm-tests (:ratio deposit) (tc/rows (tc/select-columns deposit [:Deposit :InsecticideB :InsecticideC]))
               {:family :binomial :weights (:Number deposit)}
               '(glm (formula (/ Killed Number) (+ Deposit Insecticide))
                     :family binomial :weights Number :data deposit))
    (glm-tests (:ratio deposit) (tc/rows (tc/select-columns
                                          deposit [:Deposit :InsecticideA :InsecticideB :InsecticideC]))
               {:family :binomial :weights (:Number deposit) :intercept? false}
               '(glm (formula (/ Killed Number) (+ 0 Deposit Insecticide))
                     :family binomial :weights Number :data deposit))
    (glm-tests (:ratio deposit) (tc/rows (tc/select-columns
                                          deposit [:logDep :InsecticideA :InsecticideB :InsecticideC]))
               {:family :binomial :weights (:Number deposit) :intercept? false}
               '(glm (formula (/ Killed Number) (+ 0 (log Deposit) Insecticide))
                     :family binomial :weights Number :data deposit))
    (t/is (v/delta-eq (map model (map vector newdata (repeat 0) (repeat 0)))
                      (rr/r->clj `(predict ~rmodel :type "response"
                                           :newdata (data.frame :Deposit ~newdata :Insecticide "A")))))
    (t/is (v/delta-eq (map model (map vector newdata (repeat 1) (repeat 0)))
                      (rr/r->clj `(predict ~rmodel :type "response"
                                           :newdata (data.frame :Deposit ~newdata :Insecticide "B")))))
    (t/is (v/delta-eq (map model (map vector newdata (repeat 0) (repeat 1)))
                      (rr/r->clj `(predict ~rmodel :type "response"
                                           :newdata (data.frame :Deposit ~newdata :Insecticide "C")))))
    (t/are [id v] (v/delta-eq (:dose (sut/dose model2 0.5 id 0)) v)
      1 5.099708
      2 4.514714
      3 0.8443372)
    (glm-tests (:ratio deposit)
               (tc/rows (tc/select-columns deposit
                                           [:logDep-poly-0 :logDep-poly-1 :InsecticideB :InsecticideC]))
               {:family :binomial :weights (:Number deposit)}

               '(glm (formula (/ Killed Number) (+ (poly (log Deposit) 2) Insecticide))
                     :family binomial :weights Number :data deposit))))

(t/deftest gaussian-glm
  (glm-tests ys xss {:weights weights :qres? true}
             `(glm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights))
  (glm-tests ys xss {:weights weights :offset offset :qres? true}
             `(glm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights :offset ~offset))
  (glm-tests ys xss {:weights weights :intercept? false :qres? true}
             `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights))
  (glm-tests ys xss {:weights weights :intercept? false :offset offset :qres? true}
             `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :offset ~offset))
  (glm-tests ys xss {:qres? true} `(glm (formula ~ys (* ~xss-1 ~xss-2))))
  (glm-tests ys xss {:offset offset :qres? true} `(glm (formula ~ys (* ~xss-1 ~xss-2)) :offset ~offset))
  (glm-tests ys xss {:intercept? false :qres? true} `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2)))))
  (glm-tests ys xss {:intercept? false :offset offset :qres? true}
             `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :offset ~offset)))

;; poisson

(t/deftest danishlc-data
  (let [danishlc (-> (rr/r->clj 'danishlc)
                     (mult-columns :City :Age)
                     (add-poly :AgeNum 2))]
    ;; df = 0
    (glm-tests (:Cases danishlc) (tc/rows (tc/drop-columns danishlc [:Cases :Pop :Age :City :Rate :AgeNum
                                                                     :AgeNum-poly-0 :AgeNum-poly-1]))
               {:family :poisson :offset (v/log (seq (:Pop danishlc))) :no-analysis? true}
               '(glm (formula Cases (+ (offset (log Pop)) (* City Age))) :family poisson :data danishlc))
    (glm-tests (:Cases danishlc) (tc/rows (tc/select-columns danishlc [:Age55-59 :Age60-64 :Age65-69 :Age70-74 :Age>74]))
               {:family :poisson :offset (v/log (seq (:Pop danishlc)))}
               '(glm (formula Cases (+ (offset (log Pop)) Age)) :family poisson :data danishlc))
    (glm-tests (:Cases danishlc) (:AgeNum danishlc)
               {:family :poisson :offset (v/log (seq (:Pop danishlc)))}
               '(glm (formula Cases (+ (offset (log Pop)) AgeNum)) :family poisson :data danishlc))
    (glm-tests (:Cases danishlc) (tc/rows (tc/select-columns danishlc [:AgeNum-poly-0 :AgeNum-poly-1]))
               {:family :poisson :offset (v/log (seq (:Pop danishlc)))}
               '(glm (formula Cases (+ (offset (log Pop)) (poly AgeNum 2))) :family poisson :data danishlc))))


(t/deftest counts2x2-data
  (let [rcounts2x2 (rr/r '(data.frame :Counts [263 258 151 222]
                                      :Att (gl 2 2 4 :labels ["For" "Against"])
                                      :Inc (gl 2 1 4 :labels ["High" "Low"])))
        counts2x2 (-> (rr/r->clj rcounts2x2)
                      (mult-columns :Att :Inc))]
    (glm-tests (:Counts counts2x2) (tc/rows (tc/select-columns counts2x2 [:AttAgainst :IncLow]))
               {:family :poisson :no-analysis? true}
               `(glm (formula Counts (+ Att Inc)) :family poisson :data ~rcounts2x2))
    (glm-tests (:Counts counts2x2) (tc/rows (tc/select-columns counts2x2
                                                               [:AttAgainst :IncLow :AttAgainstIncLow]))
               {:family :poisson :no-analysis? true}
               `(glm (formula Counts (* Att Inc)) :family poisson :data ~rcounts2x2))
    (glm-tests (:AttAgainst counts2x2) (:IncLow counts2x2) {:family :binomial :weights (:Counts counts2x2)}
               `(glm (formula (ifelse (== Att "Against") 1 0) Inc)
                     :family binomial :weight Counts :data ~rcounts2x2))))

(t/deftest kstones-data
  (let [kstones (-> (rr/r->clj 'kstones)
                    (one-hot :Size [:Small])
                    (one-hot :Method [:B])
                    (one-hot :Outcome [:Success])
                    (tc/* :SizeSmallMethodB [:SizeSmall :MethodB])
                    (tc/* :SizeSmallOutcomeSuccess [:SizeSmall :OutcomeSuccess])
                    (tc/* :MethodBOutcomeSuccess [:MethodB :OutcomeSuccess])
                    (tc/* :SizeSmallMethodBOutcomeSuccess [:SizeSmall :MethodB :OutcomeSuccess]))]
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones [:SizeSmall :MethodB :OutcomeSuccess]))
               {:family :poisson}
               '(glm (formula Counts (+ Size Method Outcome)) :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:SizeSmall :MethodB
                                                              :OutcomeSuccess :SizeSmallMethodB]))
               {:family :poisson}
               '(glm (formula Counts (+ (* Size Method) Outcome)) :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:SizeSmall :OutcomeSuccess
                                                              :MethodB :SizeSmallOutcomeSuccess]))
               {:family :poisson}
               '(glm (formula Counts (+ (* Size Outcome) Method)) :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:OutcomeSuccess :MethodB
                                                              :SizeSmall :MethodBOutcomeSuccess]))
               {:family :poisson}
               '(glm (formula Counts (+ (* Outcome Method) Size)) :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:SizeSmall :MethodB :OutcomeSuccess
                                                              :SizeSmallMethodB :SizeSmallOutcomeSuccess]))
               {:family :poisson}
               '(glm "Counts ~ Size * (Method + Outcome)" :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:MethodB :OutcomeSuccess :SizeSmall
                                                              :MethodBOutcomeSuccess :SizeSmallMethodB]))
               {:family :poisson}
               '(glm "Counts ~ Method * (Outcome + Size)" :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:OutcomeSuccess :MethodB :SizeSmall
                                                              :MethodBOutcomeSuccess :SizeSmallOutcomeSuccess]))
               {:family :poisson}
               '(glm "Counts ~ Outcome * (Method + Size)" :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/select-columns kstones
                                                             [:SizeSmall :MethodB :OutcomeSuccess
                                                              :SizeSmallMethodB :SizeSmallOutcomeSuccess
                                                              :MethodBOutcomeSuccess]))
               {:family :poisson :no-analysis? true}
               '(glm "Counts ~ Size * Method * Outcome - Size:Method:Outcome" :family poisson :data kstones))
    (glm-tests (:Counts kstones) (tc/rows (tc/drop-columns kstones [:Counts :Size :Method :Outcome]))
               {:family :poisson :no-analysis? true}
               '(glm (formula Counts (* Size Method Outcome)) :family poisson :data kstones))
    (glm-tests (:OutcomeSuccess kstones) (tc/rows (tc/select-columns kstones
                                                                     [:SizeSmall :MethodB
                                                                      :SizeSmallMethodB]))
               {:family :binomial :weights (:Counts kstones)}
               '(glm (formula (ifelse (== Outcome "Success") 1 0) (* Size Method))
                     :family binomial :data kstones :weights Counts))))

(t/deftest pock-data
  (let [pock (-> (rr/r->clj 'pock)
                 (tc/map-columns :l2Dilution :Dilution m/log2))]
    (glm-tests (:Count pock) (:l2Dilution pock) {:family :poisson}
               '(glm (formula Count (log2 Dilution)) :family poisson :data pock))
    (glm-tests (:Count pock) (:l2Dilution pock) {:family :quasi-poisson :no-ll? true}
               '(glm (formula Count (log2 Dilution)) :family quasipoisson :data pock))
    (glm-tests (:Count pock) (:l2Dilution pock) {:family :nbinomial :nbinomial-theta 9.892894299757403}
               '(glm.nb (formula Count (log2 Dilution))  :data pock))
    (t/is (m/delta-eq (first (:theta (rr/r->clj '(glm.nb (formula Count (log2 Dilution)) :data pock))))
                      (:nbinomial-theta (sut/glm-nbinomial (:Count pock) (:l2Dilution pock)))))))

(t/deftest hcrabs-data
  (let [hcrabs (-> (rr/r->clj 'hcrabs)
                   (tc/log :logWt :Wt)
                   (tc/log :logWidth :Width)
                   (one-hot :Spine [:NoneOK :OneOK])
                   (one-hot :Col [:DM :LM :M]))]
    (glm-tests (:Sat hcrabs) (tc/rows (tc/drop-columns hcrabs [:Wt :Width :Spine :Col :Sat]))
               {:family :quasi-poisson :no-ll? true}
               '(glm (formula Sat (+ (log Wt) (log Width) Spine Col))
                     :family quasipoisson :data hcrabs))
    (glm-tests (:Sat hcrabs) (:logWt hcrabs)
               {:family :quasi-poisson :no-ll? true}
               '(glm (formula Sat (+ (log Wt)))
                     :family quasipoisson :data hcrabs))
    (glm-tests (:Sat hcrabs) (:logWt hcrabs)
               {:family :nbinomial :nbinomial-theta 0.9580286019527684}
               '(glm.nb (formula Sat (+ (log Wt))) :data hcrabs))
    (t/is (m/delta-eq (first (:theta (rr/r->clj '(glm.nb (formula Sat (+ (log Wt))) :data hcrabs))))
                      (:nbinomial-theta (sut/glm-nbinomial (:Sat hcrabs) (:logWt hcrabs)))))))

(t/deftest lime-data
  (let [lime (-> (rr/r->clj 'lime)
                 (one-hot :Origin [:Natural :Planted])
                 (tc/log :logDBH :DBH)
                 (tc/* :OriginNaturallogDBH [:OriginNatural :logDBH])
                 (tc/* :OriginPlantedlogDBH [:OriginPlanted :logDBH]))]
    (glm-tests (:Foliage lime) (tc/rows (tc/drop-columns lime [:Foliage :DBH :Origin :Age]))
               {:family :gamma :link :log}
               '(glm (formula Foliage (* Origin (log DBH))) :family (Gamma :link log) :data lime))
    (t/is (v/delta-eq (-> (sut/glm (:Foliage lime)
                                   (tc/rows (tc/drop-columns lime [:Foliage :DBH :Origin :Age]))
                                   {:family :gamma :link :log})
                          (sut/quantile-residuals))
                      (rr/r->clj '(qresid (glm (formula Foliage (* Origin (log DBH))) :family (Gamma :link log) :data lime))) 5.0e-3))
    (glm-tests (:Foliage lime) (tc/rows (tc/drop-columns lime [:Foliage :DBH :Origin :Age]))
               {:family :inverse-gaussian :link :log}
               '(glm (formula Foliage (* Origin (log DBH))) :family (inverse.gaussian :link log) :data lime))))

(t/deftest perm-data
  (let [lime (-> (rr/r->clj 'perm)
                 (mult-columns :Mach :Day))]
    (glm-tests (:Perm lime) (tc/rows (tc/drop-columns lime [:Mach :Day :Perm]))
               {:family :inverse-gaussian :link :log}
               '(glm (formula Perm (* Mach Day)) :family (inverse.gaussian :link log) :data perm))))

(t/deftest yieldden-data
  (let [yieldden (-> (rr/r->clj 'yieldden)
                     (tc/map-columns :rDens :Dens m//)
                     (one-hot :Var [:2 :3])
                     (tc/* :DensVar2 [:Dens :Var2])
                     (tc/* :DensVar3 [:Dens :Var3])
                     (tc/* :rDensVar2 [:rDens :Var2])
                     (tc/* :rDensVar3 [:rDens :Var3]))]
    (glm-tests (:YD yieldden) (tc/rows (tc/drop-columns yieldden [:Var :YD :Yield]))
               {:family :gamma :link :inverse}
               '(glm "YD ~ (Dens + I(1/Dens)) * Var" :family (Gamma :link inverse) :data yieldden))
    (glm-tests (:YD yieldden) (tc/rows (tc/select-columns yieldden [:Dens :rDens :Var2 :Var3]))
               {:family :gamma :link :inverse}
               '(glm "YD ~ Dens + I(1/Dens) + Var" :family (Gamma :link inverse) :data yieldden))
    (t/is (v/delta-eq (-> (sut/glm (:YD yieldden)
                                   (tc/rows (tc/select-columns yieldden [:Dens :rDens :Var2 :Var3]))
                                   {:family :inverse-gaussian})
                          (sut/quantile-residuals))
                      (rr/r->clj '(qresid (glm "YD ~ Dens + I(1/Dens) + Var"
                                               :family inverse.gaussian :data yieldden)))))))

;; from Daniel, simplest case

(t/deftest lm-simplest-case
  (let [model (sut/lm [3 3 5 6] [[1 1] [2 1] [3 2] [4 2]])]
    (t/is (m/delta-eq 4.25 (model [1 2])))))

(comment (init-r))

