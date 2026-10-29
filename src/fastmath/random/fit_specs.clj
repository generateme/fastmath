(ns fastmath.random.fit-specs
  "Fitting specifications of the distributions, the values of `fastmath.random/fit-spec`.

  A specification is a map: `:params` (the fitted parameters with their domains and roles), `:valid?` (the validation of the data) and, optionally, `:estimators` (closed form estimators by method), `:initial` (starting values), `:moments` and `:lmoments`. See `fastmath.random/fit-spec` for the meaning of the keys.

  The functions of `:valid?` take the sample as given by the user. The functions of `:estimators` and `:initial` take the sample as a sorted `double[]`, they return a map of parameters or `nil` when the sample has no estimate (for example when the variance is too large for the method of moments).

  The namespace has the data checks used by the specifications and one specification per distribution kind, defined with the distribution (`fastmath.random.distributions` for the Apache Commons Math classes and the custom distributions, `fastmath.random` for the constructors which wrap SSJ distributions)."
  (:require [fastmath.core :as m]
            [fastmath.special :as special])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; data checks

(defn real-data-error
  "Message when `data` is not a sequence of at least `min-count` finite numbers, otherwise `nil`."
  [data ^long min-count]
  (cond
    (not (every? (fn [x] (and (number? x) (m/valid-double? x))) data)) "data must be finite numbers"
    (m/< (count data) min-count) (str "data must have at least " min-count " observations")))

(defn positive-data-error
  "Message when `data` is not a sequence of at least `min-count` finite numbers above `0`, otherwise `nil`."
  [data ^long min-count]
  (or (real-data-error data min-count)
      (when-not (every? m/pos? data) "data must be positive")))

(defn non-negative-data-error
  "Message when `data` is not a sequence of at least `min-count` finite numbers not below `0`, otherwise `nil`."
  [data ^long min-count]
  (or (real-data-error data min-count)
      (when-not (every? m/not-neg? data) "data must be non-negative")))

(defn unit-data-error
  "Message when `data` is not a sequence of at least `min-count` numbers inside the open interval from `0` to `1`, otherwise `nil`."
  [data ^long min-count]
  (or (real-data-error data min-count)
      (when-not (every? (fn [^double x] (m/< 0.0 x 1.0)) data) "data must be inside the interval (0, 1)")))

;; sample statistics (population versions: the maximum likelihood and moment estimates divide by n)

(defn- mean-of
  ^double [data]
  (m// (double (reduce m/+ 0.0 data)) (count data)))

(defn- central-moment-of
  "Population central moment of the order `k`."
  ^double [data ^long k]
  (let [mean (mean-of data)]
    (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/fpow (m/- x mean) k))) 0.0 data))
         (count data))))

(defn- l2-of
  "Second L-moment of the sorted `double[]`."
  ^double [^doubles xs]
  (let [n (alength xs)]
    (m// (double (loop [i (long 0)
                        sum 0.0]
                   (if (m/< i n)
                     (recur (m/inc i) (m/+ sum (m/* (m/- (m/* 2.0 (m/inc i)) n 1.0) (Array/aget xs i))))
                     sum)))
         (m/* n (m/dec n)))))

;; normal

(def normal
  "Specification of the normal distribution (`:mu`, `:sd`): the maximum likelihood estimate is the sample mean and the population standard deviation."
  {:params [{:name :mu :domain :real :role :location}
            {:name :sd :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [data] {:mu (mean-of data) :sd (m/sqrt (central-moment-of data 2))})}}})

;; exponential

(def exponential
  "Specification of the exponential distribution (`:mean`): the maximum likelihood estimate is the sample mean."
  {:params [{:name :mean :domain :positive :role :scale}]
   :valid? #(non-negative-data-error % 1)
   :estimators {:mle {:exact? true
                      :estimate (fn [data] {:mean (mean-of data)})}}})

;; gamma

(defn- gamma-moments
  [data]
  (let [mean (mean-of data)
        variance (central-moment-of data 2)]
    (when (m/pos? variance)
      {:shape (m// (m/sq mean) variance) :scale (m// variance mean)})))

(def gamma
  "Specification of the gamma distribution (`:shape`, `:scale`): the method of moments has a closed form."
  {:params [{:name :shape :domain :positive :role :shape}
            {:name :scale :domain :positive :role :scale}]
   :valid? #(positive-data-error % 2)
   :initial gamma-moments
   :estimators {:mme {:exact? true :estimate gamma-moments}}})

;; weibull

(def weibull
  "Specification of the Weibull distribution (`:alpha` the shape, `:beta` the scale): the L-moment estimate has a closed form, the starting values come from the logarithms of the sample."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(positive-data-error % 2)
   ;; the logarithm of a Weibull variable has the mean ln(beta) - gamma / alpha and the variance pi^2 / (6 alpha^2)
   :initial (fn [data]
              (let [logs (map m/log data)
                    shape (m// m/PI (m/* (m/sqrt 6.0) (m/sqrt (central-moment-of logs 2))))]
                (when (m/valid-double? shape)
                  {:alpha shape :beta (m/exp (m/+ (mean-of logs) (m// m/GAMMA shape)))})))
   ;; lambda_2 / lambda_1 = 1 - 2^(-1 / alpha), lambda_1 = beta * Gamma(1 + 1 / alpha)
   :estimators {:lmme {:exact? true
                       :estimate (fn [data]
                                   (let [l1 (mean-of data)
                                         ratio (m// (l2-of data) l1)
                                         shape (m// (m/- (m/log 2.0)) (m/log (m/- 1.0 ratio)))]
                                     (when (m/pos? shape)
                                       {:alpha shape :beta (m// l1 (special/gamma (m/inc (m// shape))))})))}}})

;; beta

(defn- beta-moments
  [data]
  (let [mean (mean-of data)
        variance (central-moment-of data 2)
        common (m/dec (m// (m/* mean (m/- 1.0 mean)) variance))]
    (when (m/pos? common)
      {:alpha (m/* mean common) :beta (m/* (m/- 1.0 mean) common)})))

(def beta
  "Specification of the beta distribution (`:alpha`, `:beta`): the method of moments has a closed form (when the variance is below `mean * (1 - mean)`)."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :shape}]
   :valid? #(unit-data-error % 2)
   :initial beta-moments
   :estimators {:mme {:exact? true :estimate beta-moments}}})

;; log-normal

(def log-normal
  "Specification of the log-normal distribution (`:scale` the mean and `:shape` the standard deviation of the logarithm): the maximum likelihood estimate is the mean and the population standard deviation of the logarithms."
  {:params [{:name :scale :domain :real :role :other}
            {:name :shape :domain :positive :role :other}]
   :valid? #(positive-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [data]
                                  (let [logs (map m/log data)]
                                    {:scale (mean-of logs) :shape (m/sqrt (central-moment-of logs 2))}))}}})

;; cauchy

(def cauchy
  "Specification of the Cauchy distribution (`:median`, `:scale`): no closed form, the starting values are the sample median and half of the interquartile range."
  {:params [{:name :median :domain :real :role :location}
            {:name :scale :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [^doubles data]
              (let [n (alength data)
                    q1 (Array/aget data (m/quot n 4))
                    q3 (Array/aget data (m/quot (m/* 3 n) 4))]
                {:median (Array/aget data (m/quot n 2)) :scale (m/max 1.0e-8 (m/* 0.5 (m/- q3 q1)))}))})

;; half-normal

(def half-normal
  "Specification of the half-normal distribution (`:mu` the lower bound, `:sigma`): the maximum likelihood estimate is the sample minimum and the root mean square distance from it."
  {:params [{:name :mu :domain :real :role :location}
            {:name :sigma :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  (let [mu (Array/aget data 0)
                                        n (alength data)
                                        sigma (m/sqrt (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/sq (m/- x mu)))) 0.0 data)) n))]
                                    (when (m/pos? sigma) {:mu mu :sigma sigma})))}}})

;; uniform

(def uniform-real
  "Specification of the continuous uniform distribution (`:lower`, `:upper`): the maximum likelihood estimate is the sample minimum and maximum."
  {:params [{:name :lower :domain :real :role :location}
            {:name :upper :domain :real :role :location}]
   :valid? #(real-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  (let [lower (Array/aget data 0)
                                        upper (Array/aget data (m/dec (alength data)))]
                                    (when (m/< lower upper) {:lower lower :upper upper})))}}})

;; chi-squared

(def chi-squared
  "Specification of the chi-squared distribution (`:degrees-of-freedom`): the starting value is the sample mean, which is the mean of the distribution."
  {:params [{:name :degrees-of-freedom :domain :positive :role :shape}]
   :valid? #(positive-data-error % 1)
   :initial (fn [data] {:degrees-of-freedom (mean-of data)})})

;; erlang

(def erlang
  "Specification of the Erlang distribution (`:k` the integer shape, which has to be fixed, `:lambda` the rate)."
  {:params [{:name :k :domain :positive :role :integer}
            {:name :lambda :domain :positive :role :other}]
   :valid? #(positive-data-error % 1)})

;; rayleigh

(def rayleigh
  "Specification of the Rayleigh distribution (`:a` the location, `:beta` the scale): no closed form with a free location, the starting values put the location slightly below the sample minimum."
  {:params [{:name :a :domain :real :role :location}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [^doubles data]
              (let [n (alength data)
                    lowest (Array/aget data 0)
                    a (m/- lowest (m// (m/- (Array/aget data (m/dec n)) lowest) n))
                    beta (m/sqrt (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/sq (m/- x a)))) 0.0 data)) (m/* 2.0 n)))]
                (when (m/pos? beta) {:a a :beta beta})))})

;; inverse gamma

(def inverse-gamma
  "Specification of the inverse gamma distribution (`:alpha` the shape, `:beta` the scale): the starting values come from the method of moments."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(positive-data-error % 2)
   :initial (fn [data]
              (let [mean (mean-of data)
                    variance (central-moment-of data 2)
                    alpha (m/+ 2.0 (m// (m/sq mean) variance))]
                (when (m/valid-double? alpha) {:alpha alpha :beta (m/* mean (m/dec alpha))})))})

;; inverse gaussian

(def inverse-gaussian
  "Specification of the inverse Gaussian distribution (`:mu` the mean, `:lambda` the shape): the maximum likelihood estimate is the sample mean and the reciprocal of the mean difference of `1 / x` and `1 / mean`."
  {:params [{:name :mu :domain :positive :role :scale}
            {:name :lambda :domain :positive :role :shape}]
   :valid? #(positive-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [data]
                                  (let [mu (mean-of data)
                                        s (m/- (mean-of (map m// data)) (m// mu))]
                                    (when (m/pos? s) {:mu mu :lambda (m// s)})))}}})

;; nakagami

(def nakagami
  "Specification of the Nakagami distribution (`:mu` the shape, not below `0.5`, `:omega` the spread): the starting values come from the method of moments."
  {:params [{:name :mu :domain [0.5 ##Inf] :role :shape}
            {:name :omega :domain :positive :role :scale}]
   :valid? #(positive-data-error % 2)
   :initial (fn [data]
              (let [squares (map m/sq data)
                    omega (mean-of squares)
                    variance (central-moment-of squares 2)]
                (when (m/pos? variance) {:mu (m/max 0.5 (m// (m/sq omega) variance)) :omega omega})))})

;; f

(def f
  "Specification of the F distribution (`:numerator-degrees-of-freedom`, `:denominator-degrees-of-freedom`): the starting values come from the method of moments, or are 5 and 10 when the sample has no such estimate."
  {:params [{:name :numerator-degrees-of-freedom :domain :positive :role :shape}
            {:name :denominator-degrees-of-freedom :domain :positive :role :shape}]
   :valid? #(positive-data-error % 2)
   :initial (fn [data]
              (let [mean (mean-of data)
                    variance (central-moment-of data 2)
                    d2 (if (m/> mean 1.0) (m// (m/* 2.0 mean) (m/dec mean)) ##NaN)
                    d1 (if (m/> d2 4.0)
                         (m// (m/* 2.0 (m/sq d2) (m/- d2 2.0))
                              (m/- (m/* variance (m/sq (m/- d2 2.0)) (m/- d2 4.0)) (m/* 2.0 (m/sq d2))))
                         ##NaN)]
                (if (m/pos? d1)
                  {:numerator-degrees-of-freedom d1 :denominator-degrees-of-freedom d2}
                  {:numerator-degrees-of-freedom 5.0 :denominator-degrees-of-freedom 10.0})))})

;; t

(def t
  "Specification of the Student t distribution (`:degrees-of-freedom`, the standard t distribution has no location and scale): the starting value comes from the excess kurtosis of the sample (`6 / (df - 4)`), or is 30 for a sample which is not heavy-tailed."
  {:params [{:name :degrees-of-freedom :domain :positive :role :shape}]
   :valid? #(real-data-error % 3)
   :initial (fn [data]
              (let [excess (m/- (m// (central-moment-of data 4) (m/sq (central-moment-of data 2))) 3.0)]
                {:degrees-of-freedom (if (m/pos? excess) (m/max 4.5 (m/+ 4.0 (m// 6.0 excess))) 30.0)}))})
