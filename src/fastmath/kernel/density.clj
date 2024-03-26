(ns fastmath.kernel.density
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.stats :as stats]
            [fastmath.optimization :as optim]
            [fastmath.calculus :as calc])
  (:import [org.apache.commons.math3.stat StatUtils]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:const true :private true :tag 'double} gaussian-factor (/ (m/sqrt m/TWO_PI)))

(defn uniform
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) 0.5 0.0))

(defn gaussian
  ^double [^double x]
  (m/* gaussian-factor (m/exp (m/* -0.5 x x))))

(defn triangular
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0) (m/- 1.0 absx) 0.0)))

(defn epanechnikov
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* 0.75 (m/- 1.0 (m/* x x))) 0.0))

(defn quartic
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* 0.9375 (m/sq (m/- 1.0 (m/* x x)))) 0.0))

(defn triweight
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0)
    (let [v (m/- 1.0 (m/* x x))]
      (m/* 1.09375 v v v)) 0.0))

(defn tricube
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0)
      (let [v (m/- 1.0 (m/* absx absx absx))]
        (m/* 0.8641975308641975 v v v)) 0.0)))

(defn cosine
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* m/QUARTER_PI (m/cos (m/* m/HALF_PI x))) 0.0))

(defn logistic
  ^double [^double x]
  (m// (m/+ 2.0 (m/exp x) (m/exp (m/- x)))))

(defn sigmoid
  ^double [^double x]
  (m/* m/TWO_INV_PI (m// (m/+ (m/exp x) (m/exp (m/- x))))))

(defn silverman
  ^double [^double x]
  (let [xx (m// (m/abs x) m/SQRT2)]
    (m/* 0.5 (m/exp (m/- xx)) (m/sin (m/+ xx m/QUARTER_PI)))))

;; experimental
(defn laplace
  ^double [^double x]
  (m/* 0.5 (m/exp (m/- (m/abs x)))))

(defn wigner
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* m/TWO_INV_PI (m/sqrt (m/- 1.0 (m/* x x)))) 0.0))

(defn cauchy ^double [^double x]
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
  * `:efficiency` - absolute efficiency"
  (let [in {:uniform {:k2 0.5 :x2k m/THIRD :kernel uniform}
            :triangular {:k2 m/TWO_THIRD :x2k m/SIXTH :kernel triangular}
            :epanechnikov {:k2 0.6 :x2k 0.2 :kernel epanechnikov}
            :quartic {:k2 (m// 5.0 7.0) :x2k (m// 7.0) :kernel quartic}
            :triweight {:k2 (m// 350.0 429.0) :x2k (m// 9.0) :kernel triweight}
            :tricube {:k2 (m// 175.0 247.0) :x2k (m// 35.0 243.0) :kernel tricube}
            :gaussian {:k2 (m/* 0.5 (m// m/SQRTPI)) :x2k 1.0 :kernel gaussian}
            :cosine {:k2 (m/* 0.0625 m/PI m/PI) :x2k (m/- 1.0 (m// 8.0 (m/* m/PI m/PI))) :kernel cosine}
            :logistic {:k2 m/SIXTH :x2k (m// (m/* m/PI m/PI) 3.0) :kernel logistic}
            :sigmoid {:k2 (m// 2.0 (* m/PI m/PI)) :x2k (m// (* m/PI m/PI) 4.0) :kernel sigmoid}
            :silverman {:k2 (m/* 0.0625 3.0 m/SQRT2) :x2k 0.0 :kernel silverman}
            :wigner {:k2 (m// 16.0 (m/* 3.0 m/PI m/PI)) :x2k 0.25 :kernel wigner}
            :cauchy {:k2 m/M_1_PI :x2k ##Inf :kernel cauchy}
            :laplace {:k2 0.25 :x2k 2.0 :kernel laplace}}]
    
    (reduce (fn [m k] (update m k calc-additional-consts)) in (keys in))))

;; 

(defn- kde-
  "Expects sorted data as double-array"
  ([{:keys [^doubles data ^long len ^double mn ^double mx ^long last-idx kernel]} ^double h]
   (let [hrev (m// h)
         span (m/* 6.0 h)
         factor (m// (* len h))]
     {:kde (fn ^double [^double x]
             (let [start (java.util.Arrays/binarySearch data (- x span))
                   start (long (if (neg? start) (dec (- start)) start))
                   end (java.util.Arrays/binarySearch data (+ x span))
                   end (min last-idx (long (if (neg? end) (dec (- end)) end)))]
               (loop [i start
                      sum 0.0]
                 (if (<= i end)
                   (recur (inc i) (+ sum ^double (kernel (* hrev (- x (aget data i))))))
                   (* factor sum)))))
      :factor factor
      :h h
      :mn (- mn span)
      :mx (+ mx span)})))

(defn- preprocess-data->map
  [^doubles adata kernel]
  (let [len (alength adata)
        len- (m/dec len)]
    {:data adata
     :len len
     :last-idx len-
     :mn (aget adata 0)
     :mx (aget adata len-)
     :kernel-name kernel
     :kernel (if (fn? kernel) kernel
                 (get-in kde-data [kernel :kernel]))}))

(defn preprocess-data
  "Returns preprocessed kde data."
  [data kernel]
  (preprocess-data->map (let [a (m/seq->double-array data)]
                          (java.util.Arrays/sort a)
                          a) kernel))

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
   (let [iqr- (/ (m/- (StatUtils/percentile data 75.0)
                      (StatUtils/percentile data 25.0)) 1.34)]
     (* scale (select-one sd iqr-) (m/pow len -0.2))))
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

(defn- array-loo
  [^doubles arr ^long cnt ^long idx]
  (let [narr (double-array (m/dec cnt))]
    (System/arraycopy arr 0 narr 0 idx)
    (System/arraycopy arr (m/inc idx) narr idx (m/- cnt idx 1))
    narr))

(defn- kde-i
  "Loo kde values"
  [{:keys [data kernel ^long len]} ^double h]
  (let [kdes (->> (range len)
                  (map (partial array-loo data len))
                  (map (fn [arr]
                         (:kde (kde- (preprocess-data->map arr kernel) h)))))]
    (map #(%1 %2) kdes data)))

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
    (calc/integrate bhx mn mx)))

(defn- ->rlcv-log
  [^double a ^double la-]
  (fn ^double [^double x]
    (if (m/>= x a)
      (m/ln x)
      (m/+ la- (m// x a)))))

(defn- rlcv-target
  [kdata]
  (let [a (rlcv-a kdata)
        la- (m/dec (m/ln a))]
    (fn ^double [^double h]
      (let [k (kde- kdata h)
            b (rlcv-b k a)
            lf (map (->rlcv-log a la-) (kde-i kdata h))]
        (m/- (stats/mean lf) b)))))

(defn- rlcv
  [kdata]
  (let [sd (m/sqrt (StatUtils/variance (:data kdata)))
        init-h (nrd kdata sd 1.06)
        nkdata (update kdata :data (fn [^doubles d] (StatUtils/normalize d)))]
    (-> (optim/maximize :lbfgsb
                        (rlcv-target nkdata)
                        {:bounds [[(m/* 0.2 init-h)
                                   (m/* 5.0 init-h)]]
                         :init [init-h]})
        ^double (ffirst)
        (m/* sd))))

(defn- lcv-target
  [kdata]
  (fn ^double [^double h]
    (->> (kde-i kdata h)
         (filter pos?)
         (map (fn [^double v] (m/log v)))
         (stats/mean))))

(defn- lcv
  [kdata]
  (let [sd (m/sqrt (StatUtils/variance (:data kdata)))
        init-h (nrd kdata sd 1.06)
        nkdata (update kdata :data (fn [^doubles d] (StatUtils/normalize d)))]
    (-> (optim/maximize :lbfgsb
                        (lcv-target nkdata)
                        {:bounds [[(m/* 0.2 init-h)
                                   (m/* 5.0 init-h)]]
                         :init [init-h]})
        ^double (ffirst)
        (m/* sd))))


(defn- lscv-target
  [kdata]
  (fn ^double [^double h]
    (let [{:keys [kde ^double mn ^double mx]} (kde- kdata h)
          in (double (calc/integrate (fn [^double v] (m/sq (kde v))) mn mx))
          lf (stats/sum (kde-i kdata h))]
      (m/- in (m/* 2.0 (stats/mean lf))))))

(defn- lscv
  [kdata]
  (let [sd (m/sqrt (StatUtils/variance (:data kdata)))
        init-h (nrd kdata sd 1.06)
        nkdata (update kdata :data (fn [^doubles d] (StatUtils/normalize d)))]
    (-> (optim/minimize :lbfgsb
                        (lscv-target nkdata)
                        {:bounds [[(m/* 0.2 init-h)
                                   (m/* 5.0 init-h)]]
                         :init [init-h]})
        ^double (ffirst)
        (m/* sd))))

(defn- infer-h
  ^double [kdata h]
  (if (keyword? h)
    (case h
      :nrd (nrd kdata 1.06)
      :nrd-adjust (nrd-adjust kdata)
      :nrd0 (nrd kdata 0.9)
      :rlcv (rlcv kdata)
      :lcv (lcv kdata)
      :lscv (lscv kdata))
    (or h (nrd kdata 1.06))))

(defn bandwidth
  "Returns infered bandwidth (h).

  h can be one of:

  * `:nrd` - rule-of-thumb (scale=1.06)
  * `:nrd0` - rule-of-thumb (scake=0.9)
  * `:nrd-adjust` - kernel specific adjustment of `:nrd`
  * `:rlcv` - robust "
  [data k h]
  (-> data
      (preprocess-data k)
      (infer-h h)))

(defn kde
  "Returns kernel density estimation function"
  ([data k] (kde data k nil))
  ([data k h] (let [kdata (preprocess-data data k)
                    h (infer-h kdata h)]
                (println "h=" h)
                (kde- kdata h))))








(m/* 1.06 (m// (m// (:delta0 (:gaussian kde-data)) (:delta0 (:quartic kde-data)))))






(require '[fastmath.calculus])
(require '[ggplot])
(require '[fastmath.stats.bootstrap :as boot])

(def mixture (let [d1 (r/distribution :laplace {:mu 0})
                 d2 (r/distribution :normal {:mu -3})
                 d3 (r/distribution :laplace {:mu 3})
                 d4 (r/distribution :gamma)]
             (r/distribution :mixture {:distrs [d1 d2 d3 d4]})))

(def mdata (r/->seq mixture 500))
(def mdata2 (repeatedly 5000 r/grand))

(defn error
  [{:keys [kde ^double mn ^double mx]} fb]
  (m/sqrt (calc/integrate (fn [x] (m/sq (m/- ^double (kde x) ^double (fb x)))) mn mx)))

(optim/minimize :lbfgsb
                (let [z (partial r/pdf mixture)]
                  (fn [^double h]
                    (error (kde mdata triangular h) z))) {:bounds [[0.001 10]]})
;; => 0.031145928740074845
;; => 0.032285504370544026

(-> (ggplot/function (:kde (kde mdata :laplace :nrd-adjust)) {:x [-5 5]})
    (ggplot/->file))



(m/* 1.03 (get-in kde-data [:triangular :nrd-scale]))

(ggplot/gg-> (ggplot/histogram {:x mdata2} :mapping (ggplot/aes :x :x))
             (ggplot/->file))


(-> (ggplot/function (partial r/pdf mixture) {:x [-7 15]})
    (ggplot/->file))
