(ns fastmath.kernel.density
  (:require [fastmath.core :as m]
            [fastmath.optimization.lbfgsb :as optim]
            [fastmath.calculus.quadrature :as calc])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.distribution NormalDistribution]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:const true :private true :tag 'double} gaussian-factor (/ (m/sqrt m/TWO_PI)))

(defn uniform
  "Uniform kernel"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) 0.5 0.0))

(defn gaussian
  "Gaussian kernel"
  ^double [^double x]
  (m/* gaussian-factor (m/exp (m/* -0.5 x x))))

(defn triangular
  "Triangular kernel"
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0) (m/- 1.0 absx) 0.0)))

(defn epanechnikov
  "Epanechnikov kernel"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* 0.75 (m/- 1.0 (m/* x x))) 0.0))

(defn quartic
  "Quartic kernel"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* 0.9375 (m/sq (m/- 1.0 (m/* x x)))) 0.0))

(defn triweight
  "Triweight kernel"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0)
    (let [v (m/- 1.0 (m/* x x))]
      (m/* 1.09375 v v v)) 0.0))

(defn tricube
  "Tricube kernel"
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0)
      (let [v (m/- 1.0 (m/* absx absx absx))]
        (m/* 0.8641975308641975 v v v)) 0.0)))

(defn cosine
  "Cosine kernel"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* m/QUARTER_PI (m/cos (m/* m/HALF_PI x))) 0.0))

(defn logistic
  "Logistic kernel"
  ^double [^double x]
  (m// (m/+ 2.0 (m/exp x) (m/exp (m/- x)))))

(defn sigmoid
  "Sigmoid kernel"
  ^double [^double x]
  (m/* m/TWO_INV_PI (m// (m/+ (m/exp x) (m/exp (m/- x))))))

(defn silverman
  "Silverman kernel"
  ^double [^double x]
  (let [xx (m// (m/abs x) m/SQRT2)]
    (m/* 0.5 (m/exp (m/- xx)) (m/sin (m/+ xx m/QUARTER_PI)))))

;; experimental
(defn laplace
  "Laplace kernel, experimental"
  ^double [^double x]
  (m/* 0.5 (m/exp (m/- (m/abs x)))))

(defn wigner
  "Wigner kernel, experimental"
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* m/TWO_INV_PI (m/sqrt (m/- 1.0 (m/* x x)))) 0.0))

(defn cauchy
  "Cauchy kernel, experimental"
  ^double [^double x]
  (m// (m/* m/HALF_PI (m/inc (m/sq (m/* 2.0 x))))))

;; Canonical kernels for density estimation
;; Marron, Nolan

;; factor used to scale nrd kernel size
(def ^{:private true :const true :tag 'double} GAUSS-FACT (m// (m/pow (m/* 0.5 (m// m/SQRTPI)) 0.2)))

(defn- calc-additional-consts
  [{:keys [^double k2 ^double x2k] :as d}]
  (let [delta0 (m/* (m/pow k2 0.2)
                    (m/pow x2k -0.4))]
    (assoc d :delta0 delta0
           :nrd-scale (m/* delta0 GAUSS-FACT)
           :efficiency (m/* k2 (m/sqrt x2k)))))

(def  kde-data
  "Collection of KDE data.

  Keys - kernel keyword, vals - kernel data

  * `:kernel` - kernel fn
  * `:k2` - integral of K(x)^2
  * `:x2k` - integral of x^2K(x)
  * `:delta0` - canonical bandwidth
  * `:efficiency` - absolute efficiency
  * `:radius` - radius of the function where function is positive or greater than 1.0e-10"
  (let [in {:uniform {:k2 0.5 :x2k m/THIRD :kernel uniform :radius 1.0}
            :triangular {:k2 m/TWO_THIRD :x2k m/SIXTH :kernel triangular :radius 1.0}
            :epanechnikov {:k2 0.6 :x2k 0.2 :kernel epanechnikov :radius 1.0}
            :quartic {:k2 (m// 5.0 7.0) :x2k (m// 7.0) :kernel quartic :radius 1.0}
            :triweight {:k2 (m// 350.0 429.0) :x2k (m// 9.0) :kernel triweight :radius 1.0}
            :tricube {:k2 (m// 175.0 247.0) :x2k (m// 35.0 243.0) :kernel tricube :radius 1.0}
            :gaussian {:k2 (m/* 0.5 (m// m/SQRTPI)) :x2k 1.0 :kernel gaussian :radius 6.65}
            :cosine {:k2 (m/* 0.0625 m/PI m/PI) :x2k (m/- 1.0 (m// 8.0 (m/* m/PI m/PI)))
                     :kernel cosine :radius 1.0}
            :logistic {:k2 m/SIXTH :x2k (m// (m/* m/PI m/PI) 3.0) :kernel logistic :radius 23.1}
            :sigmoid {:k2 (m// 2.0 (* m/PI m/PI)) :x2k (m// (* m/PI m/PI) 4.0) :kernel sigmoid :radius 22.6}
            :silverman {:k2 (m/* 0.0625 3.0 m/SQRT2) :x2k 0.0 :kernel silverman :radius 29.7}
            :wigner {:k2 (m// 16.0 (m/* 3.0 m/PI m/PI)) :x2k 0.25 :kernel wigner :radius 1.0}
            :cauchy {:k2 m/M_1_PI :x2k ##Inf :kernel cauchy :radius 30}
            :laplace {:k2 0.25 :x2k 2.0 :kernel laplace :radius 22.4}}]
    
    (reduce (fn [m k] (update m k calc-additional-consts)) in (keys in))))

;;

(defn- preprocess-data->map
  [^doubles adata kernel binned?]
  (let [len (alength adata)
        len- (m/dec len)]
    {:data adata
     :len len
     :mn (Array/aget adata 0)
     :mx (Array/aget adata len-)
     :binned? binned?
     :kernel-name kernel
     :kernel (if (fn? kernel) kernel
                 (get-in kde-data [kernel :kernel]))}))

(defn- preprocess-data
  "Returns preprocessed kde data."
  [data kernel binned?]
  (preprocess-data->map (let [a (m/seq->double-array data)]
                          (java.util.Arrays/sort a)
                          a) kernel binned?))


;; 

(defn- maybe-binned
  "Bin data and calculate averages and frequencies"
  [^doubles data ^double h ^long len binned?]
  (if-not binned?
    [data]
    (let [width (m// h (if (number? binned?) (double binned?) 5.0))]
      (loop [id (int 0)
             mx (m/+ width (Array/get data 0))
             arr []
             ws []
             curr-sum (double 0.0)
             curr-cnt (int 0)]
        (if (m/== id len)
          (if (zero? curr-cnt)
            [(double-array arr) (int-array ws)]
            [(double-array (conj arr (m// curr-sum curr-cnt))) (int-array (conj ws curr-cnt))])
          (let [curr-v (Array/get data id)]
            (if (m/< curr-v mx)
              (recur (inc id) mx arr ws (m/+ curr-sum curr-v) (m/inc curr-cnt))
              (if (m/zero? curr-cnt)
                (recur id (m/+ mx width) arr ws 0.0 0)
                (recur id (m/+ mx width)
                       (conj arr (m// curr-sum curr-cnt))
                       (conj ws curr-cnt)
                       0.0 0)))))))))

(defn- kde-
  "Expects sorted data as double-array"
  ([{:keys [^doubles data ^long len ^double mn ^double mx kernel kernel-name binned?]} ^double h]
   (let [[^doubles data ^ints weights] (maybe-binned data h len binned?)
         last-idx (m/dec (alength data))
         hrev (m// h)
         span (m/* 1.01 h (double (get-in kde-data [kernel-name :radius] 5.0)))
         -span (m/- span)
         factor (m// (m/* len h))
         kf (if binned?
              (fn ^double [^long i ^double diff]
                (m/* (Array/aget weights i)
                     ^double (kernel (m/* hrev diff))))
              (fn ^double [_ ^double diff]
                (kernel (m/* hrev diff))))]
     {:kde (fn ^double [^double x]
             (let [id (as-> (java.util.Arrays/binarySearch data x) id
                        (if (m/neg? id) (m/- (m/- id) 2) id)
                        (m/constrain (long id) 0 last-idx))
                   ^double k+ (loop [i id
                                     s (double 0.0)]
                                (if (m/== i len)
                                  s
                                  (let [diff (m/- x (Array/get data i))
                                        ^double kv (kf i diff)]
                                    (if (m/<= -span diff span)
                                      (recur (m/inc i) (m/+ s kv))
                                      (m/+ s kv)))))
                   ^double k-all (loop [i (m/dec id)
                                        s k+]
                                   (if (m/neg? i)
                                     s
                                     (let [diff (m/- x (Array/get data i))
                                           ^double kv (kf i diff)]
                                       (if (m/<= -span diff span)
                                         (recur (m/dec i) (m/+ s kv))
                                         (m/+ s kv)))))]
               (m/* factor k-all))
             #_(let [start (java.util.Arrays/binarySearch data (m/- x span))
                     start (m/max 0 (long (if (m/neg? start) (m/- (m/- start) 2) start)))
                     end (java.util.Arrays/binarySearch data (m/+ x span))
                     end (m/min last-idx (long (if (m/neg? end) (m/dec (m/- end)) end)))]
                 (loop [i start
                        sum (double 0.0)]
                   (if (m/<= i end)
                     (recur (m/inc i) (m/+ sum ^double (kf i x)))
                     (m/* factor sum)))))
      :factor factor
      :h h
      :mn (m/- mn span)
      :mx (m/+ mx span)})))

;; bandwidth estimation

(defn- select-one
  ^double [^double a ^double b]
  (let [mn (m/min a b)]
    (cond
      (m/pos? mn) mn
      (m/pos? a) a
      (m/pos? b) b
      :else 1.0)))

;; nrd0 0.9
;; nrd 1.06
(defn- nrd
  (^double [{:keys [^doubles data ^long len]} ^double sd ^double scale]
   (let [iqr- (m// (m/- (StatUtils/percentile data 75.0)
                        (StatUtils/percentile data 25.0)) 1.34)]
     (m/* scale (select-one sd iqr-) (m/pow len -0.2))))
  (^double [kdata ^double scale]
   (let [sd (m/sqrt (StatUtils/variance ^doubles (:data kdata)))]
     (nrd kdata sd scale))))

(defn- nrd-adjust
  ^double [{:keys [kernel-name] :as kdata}]
  (let [nrd-raw (nrd kdata 1.06)
        {:keys [^double nrd-scale]} (kde-data kernel-name)]
    (if (and nrd-scale (m/pos? nrd-scale) (m/valid-double? nrd-scale))
      (m/* nrd-raw nrd-scale)
      nrd-raw)))

;;

(defn- kde-loo
  "KDE with removed x0 and corrected scaling"
  [{:keys [^doubles data kernel ^long len] :as kdata} ^double h]
  (let [{kf :kde ^double factor :factor} (kde- kdata h)
        ^double k0 (kernel 0.0)
        fact (m// len (m/- len 1.0))]
    (map (fn [^double x]
           (let [^double v (kf x)]
             (m/* fact (m/- v (m/* factor k0))))) data)))


;; https://www.researchgate.net/publication/322358754_Robust_Likelihood_Cross_Validation_for_Kernel_Density_Estimation

(defn- rlcv-a
  ^double [{:keys [^long len]}]
  (m// (m/* m/INV_SQRT_2 (m/sqrt (m/ln len))) len))

(defn- rlcv-b
  ^double [{:keys [kde ^double mn ^double mx]} ^double a]
  (let [a- (m// (m/* 2.0 a))
        bhx (fn ^double [^double x]
              (let [^double v (kde x)]
                (if (m/>= v a) v (m/* v v a-))))]
    (calc/gk-quadrature bhx mn mx)))

(defn- ->rlcv-log
  [^double a ^double la-]
  (fn ^double [^double x]
    (if (m/>= x a)
      (m/ln x)
      (m/+ la- (m// x a)))))

(defn rlcv-target
  "Create target function to estimate bandwidth using Robust Likelihood Cross Validation."
  [kdata]
  (let [a (rlcv-a kdata)
        la- (m/dec (m/ln a))]
    (fn ^double [^double h]
      (let [k (kde- kdata h)
            b (rlcv-b k a)
            lf (map (->rlcv-log a la-) (kde-loo kdata h))]
        (m/- (StatUtils/mean (double-array lf)) b)))))

(defn lcv-target
  "Create target function to estimate bandwidth using Likelihood Cross Validation."
  [kdata]
  (fn ^double [^double h]
    (->> (kde-loo kdata h)
         #_ (filter (fn [^double v] (m/pos? v)))
         (map (fn [^double v] (m/log (m/max v 1.0e-3))))
         (double-array)
         (StatUtils/mean))))

(defn lscv-target
  "Create target function to estimate bandwidth using Least Squares Cross Validation."
  [kdata]
  (fn ^double [^double h]
    (let [{:keys [kde ^double mn ^double mx]} (kde- kdata h)
          in (double (calc/gk-quadrature (fn [^double v] (m/sq (kde v))) mn mx))
          lf (kde-loo kdata h)]
      (m/- in (m/* 2.0 (StatUtils/mean (double-array lf)))))))

(defn- cv
  "Cross validation method optimizer."
  [kdata opt target]
  (let [sd (m/sqrt (StatUtils/variance (:data kdata)))
        nkdata (update kdata :data (fn [^doubles d] (StatUtils/normalize d)))
        f (target nkdata)
        init-h (nrd kdata sd 1.0) 
        mn (m/* 0.2 init-h)
        mx (m/* 3.0 init-h)]
    (-> (opt f {:initial [init-h]
                :bounds [[mn mx]]})
        ^double (ffirst)
        (m/* sd))))

(defn- infer-h
  ^double [kdata h]
  (let [kdata (assoc kdata :binned? false)]
    (if (keyword? h)
      (case h
        :nrd (nrd kdata 1.06)
        :nrd-adjust (nrd-adjust kdata)
        :nrd0 (nrd kdata 0.9)
        :rlcv (cv kdata optim/maximize rlcv-target)
        :lcv (cv kdata optim/maximize lcv-target)
        :lscv (cv kdata optim/minimize lscv-target))
      (or h (nrd kdata 1.06)))))

(defn bandwidth
  "Returns infered bandwidth (h).

  h can be one of:

  * `:nrd` - rule-of-thumb (scale=1.06)
  * `:nrd0` - rule-of-thumb (scake=0.9)
  * `:nrd-adjust` - kernel specific adjustment of `:nrd`, doesn't work for `silverman` and `cauchy` 
  * `:rlcv` - robust likelihood cross-validation
  * `:lcv` - likelihood cross-validation
  * `:lscv` - least squares cross-validation"
  [kernel data h]
  (-> data
      (preprocess-data kernel nil)
      (infer-h h)))

(defn kernel-density+
  "Returns kernel density estimation function with additional information, 1d.

  Returns a map:

  * `:kde` - density function
  * `:factor` - 1/nh
  * `:h` - provided or infered bandwidth
  * `:mn` and `:mx` - infered extent of the kde`

  For arguments see [[kernel-density]]"
  ([kernel data] (kernel-density+ kernel data nil))
  ([kernel data {:keys [^double bandwidth binned?]}] (let [kdata (preprocess-data data kernel binned?)
                                                           h (infer-h kdata bandwidth)]
                                                       (kde- kdata h))))

(defn kernel-density
  "Returns kernel density estimation function, 1d
  
  Arguments:
  * `kernel` - kernel name or kernel function  
  * `data` - data
  * `params` - a map containing:
      * `:bandwidth` - bandwidth h
      * `:binned?` - if data should be binned, if `true` the width of the bin is `bandwidth` divided by 5, if is a number then it will be used as denominator. Default: `false`.

  `:bandwidth` can be a number or one of:

  * `:nrd` - rule-of-thumb (scale=1.06)
  * `:nrd0` - rule-of-thumb (scake=0.9)
  * `:nrd-adjust` - kernel specific adjustment of `:nrd`, doesn't work for `silverman` and `cauchy` 
  * `:rlcv` - robust likelihood cross-validation
  * `:lcv` - likelihood cross-validation
  * `:lscv` - least squares cross-validation"
  ([kernel data] (kernel-density kernel data nil))
  ([kernel data params] (:kde (kernel-density+ kernel data params))))

(defn kernel-density-ci
  "Create function which returns confidence intervals for given kde method.

  Check 6.1.5 http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html

  Arguments:

  * `data` - sequence of data values
  * `kernel` - kernel name  
  * `bandwidth` - as in `kde`
  * `alpha` - confidence level parameter

  Returns three values: density, lower confidence value and upper confidence value"
  ([kernel data] (kernel-density-ci kernel data nil))
  ([kernel data {:keys [^double alpha]
                 :or {alpha 0.05}
                 :as params}]
   (assert (contains? kde-data kernel) "No information about given kernel for CI calculation. ")
   (let [^double kded (:k2 (kde-data kernel))
         ^NormalDistribution local-normal (NormalDistribution.)
         za (.inverseCumulativeProbability local-normal (- 1.0 (* 0.5 (or alpha 0.05))))
         {:keys [kde ^double factor]} (kernel-density+ kernel data params)]
     (fn [^double x]
       (let [^double fx (kde x)
             band (* za (m/sqrt (* factor ^double kded fx)))]
         [fx (- fx band) (+ fx band)])))))

(comment

  (require '[ggplot])
  (require '[fastmath.random :as r])

  (def mixture (let [d1 (r/distribution :laplace {:mu 0})
                   d2 (r/distribution :normal {:mu -3})
                   d3 (r/distribution :laplace {:mu 3})
                   d4 (r/distribution :gamma)]
               (r/distribution :mixture {:distrs [d1 d2 d3 d4]})))

  (def mdata (r/->seq mixture 1000))
  (def mdata2 (repeatedly 500 r/grand))
  
  (defn error
    [{:keys [kde ^double mn ^double mx]} fb]
    (m/sqrt (calc/gk-quadrature (fn [x] (m/sq (m/- ^double (kde x) ^double (fb x)))) mn mx)))

  (optim/minimize :lbfgsb
                  (let [z (partial r/pdf mixture)]
                    (fn [^double h]
                      (error (kernel-density+ :triangular mdata h) z))) {:bounds [[0.001 10]]})
  ;; => 0.031145928740074845
  ;; => 0.032285504370544026

  (defn diff-b []
    (let [f1 (kernel-density :epanechnikov mdata {:bandwidth :nrd-adjust})
          f2 (kernel-density :epanechnikov mdata {:bandwidth :nrd-adjust
                                                  :binned? 5})]
      (fn [^double x] (m/abs (- (f1 x) (f2 x))))))
  
  (time (-> (ggplot/function (kernel-density :triangular mdata {:bandwidth :lscv
                                                                :binned? false})
                             {:x [-10 15]
                              :palette (ggplot/paletter-d :wesanderson/Royal1)})
            (ggplot/->file)))

  (require '[criterium.core :as crit])
  
  (crit/quick-bench (let [kdata (preprocess-data mdata :gaussian false)]
                      (reduce m/fast+ (kde-loo kdata 0.2))))

  (take 5 (drop 150 (kde-loo2 (preprocess-data mdata :gaussian false) 0.5)))
  ;; => (0.047971176833162335
  ;;     0.048028494833802386
  ;;     0.048073927630125524
  ;;     0.04809236379383182
  ;;     0.048167195335382065)
  ;; => (0.04797117683316232
  ;;     0.048028494833802386
  ;;     0.048073927630125524
  ;;     0.04809236379383182
  ;;     0.04816719533538206)
  
  ;; => 458.3142188267523
  ;; => 458.31421882335604
  
  (let [data (concat mdata mdata2)
        k :epanechnikov
        kdata (preprocess-data data k false)
        init-h (nrd kdata 1.0)]
    (time (println init-h (bandwidth k data :rlcv)))
    (-> (ggplot/function (rlcv-target kdata) {:x [(m/* 0.2 init-h)
                                                  (m/* 3.0 init-h)]
                                              :steps 200})
        (ggplot/->file)))


  (time ((lcv-target (preprocess-data mdata :epanechnikov false)) 0.5))
  
  (prof/profile (time (bandwidth :gaussian mdata :lcv)))

  (time (bandwidth :epanechnikov mdata :lscv))
  ;; => 0.4199722051981471
  ;; => 0.4199722051981471
  ;; => 0.5077160785575244
  ;; => 0.507716078554091

  
  (prof/profile (reduce m/fast+ (mapv epanechnikov (repeatedly 100000 r/grand))))

  (ggplot/gg-> (ggplot/histogram {:x mdata2} :mapping (ggplot/aes :x :x))
               (ggplot/->file))


  (require '[clj-async-profiler.core :as prof])

  (prof/profile (dotimes [i 10000] (reduce + (range i))))

  (prof/serve-ui 8080)
  
  (-> (ggplot/function (partial r/pdf mixture) {:x [-7 15]})
      (ggplot/->file)))
;; => nil
;; => nil
;; => nil
;; => nil
