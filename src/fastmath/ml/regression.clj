(ns fastmath.ml.regression
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot]
            [fastmath.random :as r])
  (:import [org.apache.commons.math3.stat.regression OLSMultipleLinearRegression]
           [fastmath.java Array]
           [clojure.lang IFn]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn predict [model xs] (prot/predict model xs))

(defrecord LMData [^OLSMultipleLinearRegression model
                   ^double intercept beta coefficients
                   residual fitted df ^long observations
                   ^double r-squared ^double adjusted-r-squared
                   ^double stderr ^double tss ^double rss ^double ess ^double mse
                   ^double f-statistic ^double p-value
                   ^double log-likelihood ^double aic ^double bic]
  IFn
  (invoke [_ xs] (m/+ (v/dot beta xs) intercept))
  prot/PredictProto
  (predict [_ xs] (m/+ (v/dot beta xs) intercept)))

(defn- residual-analysis
  [residuals ^double rss ^long df]
  (let [skew (stats/skewness residuals :g1)
        kurt (stats/kurtosis residuals :kurt)
        normality (stats/normality-test residuals skew kurt nil)]
    {:ss rss :df df
     :residuals residuals
     :ms (m// rss df)
     :skewness skew
     :kurtosis kurt
     :durbin-watson (stats/durbin-watson residuals)
     :jarque-berra (dissoc (stats/jarque-bera-test residuals skew kurt nil) :skewness :kurtosis)
     :normality (dissoc normality :skewness :kurtosis)}))

(defn- lm-analysis
  [beta stderrs ^long df]
  (let [tdistr (r/distribution :t {:degrees-of-freedom df})]
    (map (fn [^double b ^double err]
           (let [t-value (m// b err)]
             {:estimate b :stderr err :t-value t-value
              :p-value (stats/p-value tdistr t-value :both)})) beta stderrs)))

(defn lm
  ([ys xss] (lm ys xss nil))
  ([ys xss {:keys [^double qr-threshold]
            :or {qr-threshold 0.0}}]
   (let [model (OLSMultipleLinearRegression. qr-threshold)
         xss (if (number? (first xss)) (map vector xss) xss)]
     (.newSampleData model (m/seq->double-array ys) (m/seq->double-double-array xss))
     (let [^doubles beta (.estimateRegressionParameters model)
           ^doubles stderrs (.estimateRegressionParametersStandardErrors model)
           n (count xss)           
           intercept (Array/aget beta 0)
           coefficients (vec (rest beta))
           residuals (.estimateResiduals model)
           fitted (map (fn [xs] (m/+ (v/dot xs coefficients) intercept)) xss)
           df (m/- (count xss) (alength beta))
           model-df (count coefficients)
           stderr (.estimateRegressionStandardError model)
           tss (.calculateTotalSumOfSquares model)
           rss (.calculateResidualSumOfSquares model)
           ess (m/- tss rss)
           r2 (.calculateRSquared model)
           f-statistic (m// (m// ess model-df) (m/sq stderr))
           log-likelihood (m/* (m// n -2) (m/+ m/LOG_TWO_PI (m/log (m// rss n)) 1.0))]
       (->LMData model intercept coefficients (lm-analysis beta stderrs df)
                 (residual-analysis (seq residuals) rss df)
                 fitted {:residuals df :model model-df} n
                 r2 (.calculateAdjustedRSquared model)
                 stderr tss rss ess (m// ess model-df)
                 f-statistic (stats/p-value (r/distribution :f {:numerator-degrees-of-freedom model-df
                                                                :denominator-degrees-of-freedom df})
                                            f-statistic :one-sided-greater)
                 log-likelihood (m/* -2.0 (m/- log-likelihood (m/inc model-df)))
                 (m/+ (m/* -2.0 log-likelihood)
                      (m/* (m/log n) (m/inc model-df))))))))


#_(defn glm
    )


(comment
  (def dist [   2  10   4  22  16  10  18  26  34  17  28  14  20  24  28  26  34  34  46
           26  36  60  80  20  26  54  32  40  32  40  50  42  56  76  84  36  46  68
           32  48  52  56  64  66  54  70  92  93 120  85])

  (def speed [  4  4  7  7  8  9 10 10 10 11 11 12 12 12 12 13 13 13 13 14 14 14 14 15 15
            15 16 16 17 17 17 18 18 18 18 19 19 19 20 20 20 20 20 22 23 24 24 24 24 25])

  (def ctl [4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14])
  (def trt [4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69])
  (def weight (concat ctl trt))
  (def group (concat (repeat 10 0) (repeat 10 1)))

  (lm group weight)

  (dissoc (:residual (lm group weight)) :residuals)
  ;; => {:variance 0.45943421052631583,
  ;;     :skewness 0.49076129003476227,
  ;;     :kurtosis 2.545241591026403,
  ;;     :durbin-watson 2.7630301572300033,
  ;;     :jarque-berra {:p-value 0.6141108002291311, :Z 0.9751598214321114},
  ;;     :normality {:p-value 0.5591574636998055, :Z 1.162648314601708}}

  (dissoc (lm group weight) :residual :fitted)
  ;; => {:intercept 5.031999999999999,
  ;;     :beta [-0.37099999999999966],
  ;;     :adjusted-r-squared 0.021581910045406882,
  ;;     :r-squared 0.07307759899038546,
  ;;     :f-statistic 1.4191012973623132,
  ;;     :rss 8.729249999999999,
  ;;     :p-value 0.24902316597300622,
  ;;     :df {:residuals 18, :model 1},
  ;;     :tss 9.417454999999999,
  ;;     :ess 0.6882049999999983,
  ;;     :coefficients
  ;;     ({:estimate 5.031999999999999,
  ;;       :stderr 0.2202176953229084,
  ;;       :t-value 22.850116529561824,
  ;;       :p-value 9.547918011776346E-15}
  ;;      {:estimate -0.37099999999999966,
  ;;       :stderr 0.3114348514002032,
  ;;       :t-value -1.1912603818487013,
  ;;       :p-value 0.2490231659730061}),
  ;;     :stderr 0.6963894982933999,
  ;;     :model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x4a1f3279 "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@4a1f3279"]}

  (prot/predict (lm group weight) [1])

  (dissoc (lm speed dist) :residual :fitted)
  ;; => {:intercept -17.57909489051092,
  ;;     :beta [3.932408759124086],
  ;;     :adjusted-r-squared 0.6438102011907144,
  ;;     :r-squared 0.6510793807582509,
  ;;     :df 48,
  ;;     :coefficients
  ;;     ({:estimate -17.57909489051092,
  ;;       :std-err 6.7584401693792335,
  ;;       :t-value -2.601058003022252,
  ;;       :p-value 0.012318816153808889}
  ;;      {:estimate 3.932408759124086,
  ;;       :std-err 0.4155127766571222,
  ;;       :t-value 9.463989990298368,
  ;;       :p-value 1.48991929904696E-12}),
  ;;     :model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x3c3845ab "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@3c3845ab"]}


  (def xss [[ 0.00000000e+00,  0.00000000e+00,  2.50000000e+01],
          [ 4.08163265e-01,  3.96924149e-01,  2.10849646e+01],
          [ 8.16326531e-01,  7.28634783e-01,  1.75031237e+01],
          [ 1.22448980e+00,  9.40632785e-01,  1.42544773e+01],
          [ 1.63265306e+00,  9.98087482e-01,  1.13390254e+01],
          [ 2.04081633e+00,  8.91559230e-01,  8.75676801e+00],
          [ 2.44897959e+00,  6.38550320e-01,  6.50770512e+00],
          [ 2.85714286e+00,  2.80629400e-01,  4.59183673e+00],
          [ 3.26530612e+00, -1.23398137e-01,  3.00916285e+00],
          [ 3.67346939e+00, -5.07151709e-01,  1.75968347e+00],
          [ 4.08163265e+00, -8.07581691e-01,  8.43398584e-01],
          [ 4.48979592e+00, -9.75328286e-01,  2.60308205e-01],
          [ 4.89795918e+00, -9.82831204e-01,  1.04123282e-02],
          [ 5.30612245e+00, -8.28857736e-01,  9.37109538e-02],
          [ 5.71428571e+00, -5.38705288e-01,  5.10204082e-01],
          [ 6.12244898e+00, -1.60045086e-01,  1.25989171e+00],
          [ 6.53061224e+00,  2.44910071e-01,  2.34277384e+00],
          [ 6.93877551e+00,  6.09627196e-01,  3.75885048e+00],
          [ 7.34693878e+00,  8.74184299e-01,  5.50812162e+00],
          [ 7.75510204e+00,  9.95115395e-01,  7.59058726e+00],
          [ 8.16326531e+00,  9.52551848e-01,  1.00062474e+01],
          [ 8.57142857e+00,  7.53486727e-01,  1.27551020e+01],
          [ 8.97959184e+00,  4.30625870e-01,  1.58371512e+01],
          [ 9.38775510e+00,  3.70144015e-02,  1.92523948e+01],
          [ 9.79591837e+00, -3.62678429e-01,  2.30008330e+01],
          [ 1.02040816e+01, -7.02784220e-01,  2.70824656e+01],
          [ 1.06122449e+01, -9.27424552e-01,  3.14972928e+01],
          [ 1.10204082e+01, -9.99691655e-01,  3.62453145e+01],
          [ 1.14285714e+01, -9.07712248e-01,  4.13265306e+01],
          [ 1.18367347e+01, -6.66598288e-01,  4.67409413e+01],
          [ 1.22448980e+01, -3.15964115e-01,  5.24885464e+01],
          [ 1.26530612e+01,  8.65820672e-02,  5.85693461e+01],
          [ 1.30612245e+01,  4.74903061e-01,  6.49833403e+01],
          [ 1.34693878e+01,  7.85198826e-01,  7.17305289e+01],
          [ 1.38775510e+01,  9.66488646e-01,  7.88109121e+01],
          [ 1.42857143e+01,  9.88987117e-01,  8.62244898e+01],
          [ 1.46938776e+01,  8.48997803e-01,  9.39712620e+01],
          [ 1.51020408e+01,  5.69520553e-01,  1.02051229e+02],
          [ 1.55102041e+01,  1.96472687e-01,  1.10464390e+02],
          [ 1.59183673e+01, -2.08855085e-01,  1.19210746e+02],
          [ 1.63265306e+01, -5.79868557e-01,  1.28290296e+02],
          [ 1.67346939e+01, -8.55611267e-01,  1.37703040e+02],
          [ 1.71428571e+01, -9.90779466e-01,  1.47448980e+02],
          [ 1.75510204e+01, -9.63165404e-01,  1.57528113e+02],
          [ 1.79591837e+01, -7.77305991e-01,  1.67940441e+02],
          [ 1.83673469e+01, -4.63737404e-01,  1.78685964e+02],
          [ 1.87755102e+01, -7.39780734e-02,  1.89764681e+02],
          [ 1.91836735e+01,  3.27935645e-01,  2.01176593e+02],
          [ 1.95918367e+01,  6.75970465e-01,  2.12921699e+02],
          [ 2.00000000e+01,  9.12945251e-01,  2.25000000e+02]])

  (def ys [ 4.48933467,  5.34231784,  5.81681264,  6.69809894,  5.94032168,
         6.65120265,  6.59592036,  6.07958661,  6.36186999,  6.7775208 ,
         6.6836138 ,  6.56551043,  7.17869566,  7.25158876,  8.04903067,
         7.84930564,  8.39547421,  8.25715349,  9.42432217,  8.43407786,
         9.70661341,  8.17868934, 10.44401028,  8.72852257,  9.32853046,
         9.70458442,  8.58017153,  9.78906458, 10.20670793,  9.20752194,
         9.26055989, 10.66574522,  9.89388631, 11.87907875, 10.65688892,
         11.09346875, 10.77674283, 10.11844363, 10.75104318, 10.19883781,
         10.68264587, 10.73214468, 10.13324477,  9.6745733 , 10.05157595,
         10.74993283, 10.4821192 , 10.38354809, 11.52270194, 11.19546067])

  (:residual (lm ys xss))
  ;; => {:model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x5a0f9d15 "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@5a0f9d15"],
  ;;     :intercept 5.205844957505818,
  ;;     :beta
  ;;     [0.4687244824101058 0.48360119133943785 -0.017404792025998428],
  ;;     :coefficients
  ;;     ({:estimate 5.205844957505818,
  ;;       :stderr 0.17121765098613817,
  ;;       :t-value 30.404838096553984,
  ;;       :p-value 0.0}
  ;;      {:estimate 0.4687244824101058,
  ;;       :stderr 0.02640602380582178,
  ;;       :t-value 17.75066499435501,
  ;;       :p-value 0.0}
  ;;      {:estimate 0.48360119133943785,
  ;;       :stderr 0.10380518021321035,
  ;;       :t-value 4.65873851715441,
  ;;       :p-value 2.7357988015230816E-5}
  ;;      {:estimate -0.017404792025998428,
  ;;       :stderr 0.0023184680741556492,
  ;;       :t-value -7.507022511982179,
  ;;       :p-value 1.5906164352043508E-9}),
  ;;     :residual
  ;;     {:skewness 0.20739432426143098,
  ;;      :durbin-watson 2.895610842471915,
  ;;      :df 46,
  ;;      :residuals
  ;;      (-0.281390486855857
  ;;       0.1201831997789844
  ;;       0.18060503046724197
  ;;       0.7115107121796109
  ;;       -0.3161106545012595
  ;;       0.21002773477713443
  ;;       0.04664027043646435
  ;;       -0.5212639045929288
  ;;       -0.2624545488744907
  ;;       0.12571689950623455
  ;;       -0.03016566394419673
  ;;       -0.26861126500748433
  ;;       0.15253688652965813
  ;;       -0.03889808824475338
  ;;       0.43415861572664927
  ;;       -0.12695489929314618
  ;;       0.050908099821604935
  ;;       -0.43045985438778267
  ;;       0.4478982783433221
  ;;       -0.7558996851690925
  ;;       0.38794759174036564
  ;;       -1.1871812187495188
  ;;       1.0966019242028828
  ;;       -0.5604093189118498
  ;;       0.10681516991152051
  ;;       0.5270685634295127
  ;;       -0.6031849766769817
  ;;       0.5319787291205937
  ;;       0.8022621531175034
  ;;       -0.41060627454152687
  ;;       -0.6184156884468752
  ;;       0.5066167977434599
  ;;       -0.5327165039046751
  ;;       1.22853381728563
  ;;       -0.14941147324256043
  ;;       0.2140037158843029
  ;;       -0.09150837760802233
  ;;       -0.6653379860026725
  ;;       0.10278114237650193
  ;;       -0.29249481169484604
  ;;       0.337447353275687
  ;;       0.49280638651751474
  ;;       -0.0624160472758728
  ;;       -0.5503326284207137
  ;;       -0.27330354015398584
  ;;       0.2691187059729234
  ;;       -0.1856763668723076
  ;;       -0.47130760645851844
  ;;       0.5126412325563283
  ;;       0.08970285913042986),
  ;;      :ss 11.607740226533197,
  ;;      :kurtosis 3.0259454940918746,
  ;;      :ms 0.2523421788376782,
  ;;      :normality {:p-value 0.7208635552500905, :Z 0.6546108067291813},
  ;;      :jarque-berra
  ;;      {:p-value 0.8353373910613211, :Z 0.35983914918144544}},
  ;;     :fitted
  ;;     (4.770725156855857
  ;;      5.2221346402210145
  ;;      5.636207609532758
  ;;      5.9865882278203895
  ;;      6.256432334501259
  ;;      6.441174915222865
  ;;      6.549280089563535
  ;;      6.600850514592928
  ;;      6.6243245388744905
  ;;      6.651803900493766
  ;;      6.7137794639441974
  ;;      6.834121695007484
  ;;      7.026158773470343
  ;;      7.290486848244754
  ;;      7.614872054273352
  ;;      7.976260539293145
  ;;      8.344566110178395
  ;;      8.687613344387781
  ;;      8.976423891656678
  ;;      9.189977545169095
  ;;      9.318665818259635
  ;;      9.365870558749517
  ;;      9.347408355797118
  ;;      9.28893188891185
  ;;      9.22171529008848
  ;;      9.177515856570487
  ;;      9.18335650667698
  ;;      9.257085850879404
  ;;      9.404445776882497
  ;;      9.618128214541526
  ;;      9.878975578446877
  ;;      10.15912842225654
  ;;      10.426602813904674
  ;;      10.650544932714368
  ;;      10.80630039324256
  ;;      10.879465034115698
  ;;      10.868251207608022
  ;;      10.783781616002674
  ;;      10.648262037623496
  ;;      10.491332621694847
  ;;      10.345198516724313
  ;;      10.239338293482485
  ;;      10.19566081727587
  ;;      10.224905928420714
  ;;      10.324879490153986
  ;;      10.480814124027077
  ;;      10.667795566872307
  ;;      10.854855696458518
  ;;      11.010060707443671
  ;;      11.10575781086957),
  ;;     :df {:residuals 46, :model 3},
  ;;     :observations 50,
  ;;     :r-squared 0.9325049619539212,
  ;;     :adjusted-r-squared 0.9281031116465681,
  ;;     :stderr 0.5023367185839377,
  ;;     :tss 171.97916413660798,
  ;;     :rss 11.607740226533197,
  ;;     :ess 160.37142391007478,
  ;;     :mse 53.45714130335826,
  ;;     :f-statistic 211.84386038667412,
  ;;     :p-value 0.0,
  ;;     :log-likelihood -34.438154937049326,
  ;;     :aic 76.87630987409865,
  ;;     :bic 84.52440189581124}

  (defn convert-x [x] [x (m/sin x) (m/sq (- x 5))])

  (predict (lm xss ys) (convert-x 2.123))
  ;; => 6.468609171503275
  )
