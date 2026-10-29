(ns fastmath.random.fit-specs
  "Fitting specifications of the distributions, the values of `fastmath.random/fit-spec`.

  A specification is a map: `:params` (the fitted parameters with their domains and roles), `:valid?` (the validation of the data) and, optionally, `:estimators` (closed form estimators by method), `:initial` (starting values), `:moments` and `:lmoments`. See `fastmath.random/fit-spec` for the meaning of the keys.

  The functions of `:valid?` take the sample as given by the user. The functions of `:estimators` and `:initial` take the sample as a sorted `double[]`, they return a map of parameters or `nil` when the sample has no estimate (for example when the variance is too large for the method of moments). An estimator which is declared for fixed parameters (the `:fixed` set of its entry) takes also the map of the fixed parameters and returns the other parameters. A parameter with the `:integer` role may have a `:range`, a function of the sorted sample returning `[lo hi start]`, the integers between which the parameter is searched from `start`.

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

(defn count-data-error
  "Message when `data` is not a sequence of at least `min-count` integers not below `lower`, otherwise `nil`."
  [data ^long min-count ^long lower]
  (cond
    (not (every? (fn [x] (and (number? x) (m/valid-double? x) (m/== (double x) (m/floor (double x))))) data)) "data must be integers"
    (not (every? (fn [^double x] (m/>= x lower)) data)) (str "data must not be below " lower)
    (m/< (count data) min-count) (str "data must have at least " min-count " observations")))

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

;; chi

(def chi
  "Specification of the chi distribution (`:nu` the integer degrees of freedom): `nu` is searched from the value given by the sample mean, `mean^2 + 1/2` (the mean of the distribution is about `sqrt(nu - 1/2)`)."
  {:params [{:name :nu :domain :positive :role :integer
             :range (fn [data] [1 100000 (m/max 1 (long (m/round (m/+ (m/sq (mean-of data)) 0.5))))])}]
   :valid? #(positive-data-error % 1)})

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

;; discrete distributions

(defn- variance-of
  ^double [data]
  (central-moment-of data 2))

(def poisson
  "Specification of the Poisson distribution (`:p` the rate): the maximum likelihood estimate is the sample mean."
  {:params [{:name :p :domain :positive :role :scale}]
   :valid? #(count-data-error % 1 0)
   :estimators {:mle {:exact? true
                      :estimate (fn [data] (let [mean (mean-of data)] (when (m/pos? mean) {:p mean})))}}})

(def geometric
  "Specification of the geometric distribution (`:p`, the number of failures before the first success): the maximum likelihood estimate is `1 / (1 + mean)`."
  {:params [{:name :p :domain :unit :role :probability}]
   :valid? #(count-data-error % 1 0)
   :estimators {:mle {:exact? true
                      :estimate (fn [data] {:p (m// (m/inc (mean-of data)))})}}})

(def binomial
  "Specification of the binomial distribution (`:trials` an integer, `:p`): for given trials the maximum likelihood estimate is `mean / trials`. Without fixed trials they are searched from the largest observation, starting from the method of moments estimate (a sample which is not under-dispersed has no finite maximum, the search stops at the limit of `1000000` trials)."
  {:params [{:name :trials :domain :positive :role :integer
             :range (fn [^doubles data]
                      (let [lo (m/max 1 (long (Array/aget data (m/dec (alength data)))))
                            mean (mean-of data)
                            variance (variance-of data)
                            start (if (m/< variance mean) (long (m/round (m// (m/sq mean) (m/- mean variance)))) lo)]
                        [lo 1000000 (m/max lo start)]))}
            {:name :p :domain :unit :role :probability}]
   :valid? #(count-data-error % 1 0)
   :estimators {:mle {:exact? true
                      :fixed #{:trials}
                      :estimate (fn [data fixed]
                                  (let [p (m// (mean-of data) (double (:trials fixed)))]
                                    (when (m/<= 0.0 p 1.0) {:p p})))}}})

(def pascal
  "Specification of the Pascal distribution (`:r` the integer number of successes, `:p`; the number of failures before `r` successes): for given `r` the maximum likelihood estimate is `r / (r + mean)`. Without fixed `r` it is searched, starting from the method of moments estimate (up to `1000000`)."
  {:params [{:name :r :domain :positive :role :integer
             :range (fn [data]
                      (let [mean (mean-of data)
                            variance (variance-of data)
                            start (if (m/> variance mean) (long (m/round (m// (m/sq mean) (m/- variance mean)))) 1)]
                        [1 1000000 (m/max 1 start)]))}
            {:name :p :domain :unit :role :probability}]
   :valid? #(count-data-error % 1 0)
   :estimators {:mle {:exact? true
                      :fixed #{:r}
                      :estimate (fn [data fixed]
                                  (let [r (double (:r fixed))] {:p (m// r (m/+ r (mean-of data)))}))}}})

(def negative-binomial
  "Specification of the negative binomial distribution (`:r` real, `:p`): the starting values come from the method of moments."
  {:params [{:name :r :domain :positive :role :shape}
            {:name :p :domain :unit :role :probability}]
   :valid? #(count-data-error % 2 0)
   :initial (fn [data]
              (let [mean (mean-of data)
                    variance (variance-of data)]
                (when (m/> variance mean)
                  {:r (m// (m/sq mean) (m/- variance mean)) :p (m// mean variance)})))})

(def nbi
  "Specification of the negative binomial distribution type I (`:mu` the mean, `:sigma`, the variance is `mu + sigma mu^2`): the starting values come from the method of moments."
  {:params [{:name :mu :domain :positive :role :scale}
            {:name :sigma :domain :positive :role :shape}]
   :valid? #(count-data-error % 2 0)
   :initial (fn [data]
              (let [mean (mean-of data)
                    variance (variance-of data)]
                (when (m/pos? mean)
                  {:mu mean :sigma (if (m/> variance mean) (m// (m/- variance mean) (m/sq mean)) 0.1)})))})

(def nbii
  "Specification of the negative binomial distribution type II (`:mu` the mean, `:sigma`, the variance is `mu (1 + sigma)`): the starting values come from the method of moments."
  {:params [{:name :mu :domain :positive :role :scale}
            {:name :sigma :domain :positive :role :shape}]
   :valid? #(count-data-error % 2 0)
   :initial (fn [data]
              (let [mean (mean-of data)
                    variance (variance-of data)]
                (when (m/pos? mean)
                  {:mu mean :sigma (if (m/> variance mean) (m// (m/- variance mean) mean) 0.1)})))})

(def logarithmic
  "Specification of the logarithmic distribution (`:p`, support from `1`): the starting value is chosen from the sample mean, `1 - 1 / mean` kept away from the ends."
  {:params [{:name :p :domain :unit :role :probability}]
   :valid? #(count-data-error % 1 1)
   :initial (fn [data] {:p (m/constrain (m/- 1.0 (m// (mean-of data))) 0.05 0.95)})})

(def uniform-int
  "Specification of the discrete uniform distribution (`:lower`, `:upper`, integers): the maximum likelihood estimate is the sample minimum and maximum."
  {:params [{:name :lower :domain :real :role :integer
             :range (fn [^doubles data] (let [mn (long (Array/aget data 0))] [mn mn mn]))}
            {:name :upper :domain :real :role :integer
             :range (fn [^doubles data] (let [mx (long (Array/aget data (m/dec (alength data))))] [mx mx mx]))}]
   :valid? #(count-data-error % 1 Long/MIN_VALUE)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  {:lower (Array/aget data 0) :upper (Array/aget data (m/dec (alength data)))})}}})

(def zipf
  "Specification of the Zipf distribution (`:number-of-elements` the integer highest rank, `:exponent`): the number of elements is searched at the largest observation, where the likelihood is the largest."
  {:params [{:name :number-of-elements :domain :positive :role :integer
             :range (fn [^doubles data] (let [mx (long (Array/aget data (m/dec (alength data))))] [mx mx mx]))}
            {:name :exponent :domain :non-negative :role :shape}]
   :valid? #(count-data-error % 1 1)})

(def beta-binomial
  "Specification of the beta-binomial distribution (`:alpha`, `:beta`, `:n` the integer number of trials, which has to be fixed)."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :shape}
            {:name :n :domain :positive :role :integer}]
   :valid? #(count-data-error % 2 0)})

(def bb
  "Specification of the beta-binomial distribution with the gamlss parameters (`:mu` the mean probability, `:sigma`, `:bd` the integer binomial denominator, which has to be fixed)."
  {:params [{:name :mu :domain :unit :role :probability}
            {:name :sigma :domain :positive :role :shape}
            {:name :bd :domain :positive :role :integer}]
   :valid? #(count-data-error % 2 0)})

;; heavy-tailed, location-scale and bounded distributions

(defn- sd-of
  ^double [data]
  (m/sqrt (central-moment-of data 2)))

(defn- median-of
  "Median of the sorted `double[]`."
  ^double [^doubles xs]
  (let [n (alength xs)
        h (m/quot n 2)]
    (if (m/even? n)
      (m/* 0.5 (m/+ (Array/aget xs (m/dec h)) (Array/aget xs h)))
      (Array/aget xs h))))

(defn- quantile-of
  "Value of the sorted `double[]` at the probability `p` (the nearest rank)."
  ^double [^doubles xs ^double p]
  (Array/aget xs (m/constrain (long (m/round (m/* p (m/dec (alength xs))))) 0 (m/dec (alength xs)))))

(defn- below-minimum
  "A point a little below the minimum of the sorted `double[]`: a tenth of the distance of the minimum from the median (a hundredth of the range when they are equal)."
  ^double [^doubles xs]
  (let [lowest (Array/aget xs 0)
        gap (m/- (median-of xs) lowest)
        gap (if (m/pos? gap) gap (m/* 0.01 (m/- (Array/aget xs (m/dec (alength xs))) lowest)))]
    (m/- lowest (m/* 0.1 (if (m/pos? gap) gap 1.0)))))

(defn- spacing-of
  "A tenth of the average distance of neighbouring observations of the sorted `double[]` (`1.0` for a constant sample)."
  ^double [^doubles xs]
  (let [n (alength xs)
        range (m/- (Array/aget xs (m/dec n)) (Array/aget xs 0))]
    (if (m/pos? range) (m// range (m/* 10.0 n)) 1.0)))

(def laplace
  "Specification of the Laplace distribution (`:mu`, `:beta`): the maximum likelihood estimate is the sample median and the mean absolute deviation from it."
  {:params [{:name :mu :domain :real :role :location}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  (let [mu (median-of data)
                                        beta (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/abs (m/- x mu)))) 0.0 data))
                                                  (alength data))]
                                    (when (m/pos? beta) {:mu mu :beta beta})))}}})

(def logistic
  "Specification of the logistic distribution (`:mu`, `:s`): the starting values are the sample mean and the scale from the standard deviation."
  {:params [{:name :mu :domain :real :role :location}
            {:name :s :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [data] {:mu (mean-of data) :s (m// (m/* (sd-of data) (m/sqrt 3.0)) m/PI)})})

(def gumbel
  "Specification of the Gumbel distribution of maxima (`:mu`, `:beta`): the starting values come from the method of moments."
  {:params [{:name :mu :domain :real :role :location}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [data]
              (let [beta (m// (m/* (sd-of data) (m/sqrt 6.0)) m/PI)]
                (when (m/pos? beta) {:mu (m/- (mean-of data) (m/* m/GAMMA beta)) :beta beta})))})

(def levy
  "Specification of the Levy distribution (`:mu` the lower end of the support, `:c`): the starting values put `mu` a little below the minimum and `c` is the maximum likelihood estimate for it."
  {:params [{:name :mu :domain :real :role :lower-bound}
            {:name :c :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [^doubles data]
              (let [mu (below-minimum data)]
                {:mu mu :c (m// (alength data) (double (reduce (fn [^double s ^double x] (m/+ s (m// (m/- x mu)))) 0.0 data)))}))})

(def pareto
  "Specification of the Pareto distribution (`:scale` the lowest value, `:shape`): the maximum likelihood estimate is the sample minimum and the reciprocal of the mean logarithm of the observations relative to it."
  {:params [{:name :scale :domain :positive :role :lower-bound}
            {:name :shape :domain :positive :role :shape}]
   :valid? #(positive-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  (let [lowest (Array/aget data 0)
                                        s (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/log (m// x lowest)))) 0.0 data))
                                               (alength data))]
                                    (when (m/pos? s) {:scale lowest :shape (m// s)})))}}})

(defn- log-moments
  "Mean and standard deviation of `ln(x - location)`."
  [data ^double location]
  (let [logs (map (fn [^double x] (m/log (m/- x location))) data)]
    [(mean-of logs) (sd-of logs)]))

(def frechet
  "Specification of the Frechet distribution (`:alpha` the shape, `:beta` the scale, `:delta` the lower end of the support): the starting values put `delta` a little below the minimum and use the Gumbel distribution of the logarithms."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :scale}
            {:name :delta :domain :real :role :lower-bound}]
   :valid? #(real-data-error % 3)
   :initial (fn [^doubles data]
              (let [delta (below-minimum data)
                    [^double mean ^double sd] (log-moments data delta)
                    alpha (m// m/PI (m/* (m/sqrt 6.0) sd))]
                (when (m/valid-double? alpha)
                  {:alpha alpha :beta (m/exp (m/- mean (m// m/GAMMA alpha))) :delta delta})))})

(def log-logistic
  "Specification of the log-logistic distribution (`:alpha` the shape, `:beta` the scale, the median): the starting values come from the median and the quartiles."
  {:params [{:name :alpha :domain :positive :role :shape}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(positive-data-error % 3)
   :initial (fn [^doubles data]
              (let [q1 (quantile-of data 0.25)
                    q3 (quantile-of data 0.75)
                    alpha (m// (m/* 2.0 (m/log 3.0)) (m/log (m// q3 q1)))]
                (when (m/pos? alpha) {:alpha alpha :beta (median-of data)})))})

(def triangular
  "Specification of the triangular distribution (`:a` the lower end, `:c` the mode, `:b` the upper end): the starting values put the ends slightly outside of the sample and the mode from the mean."
  {:params [{:name :a :domain :real :role :lower-bound}
            {:name :c :domain :real :role :location}
            {:name :b :domain :real :role :upper-bound}]
   :valid? #(real-data-error % 3)
   :initial (fn [^doubles data]
              (let [lowest (Array/aget data 0)
                    highest (Array/aget data (m/dec (alength data)))
                    h (spacing-of data)
                    a (m/- lowest h)
                    b (m/+ highest h)]
                {:a a :b b :c (m/constrain (m/- (m/* 3.0 (mean-of data)) a b) lowest highest)}))})

(def power
  "Specification of the power-function distribution (`:a` the lower end, `:b` the upper end, `:c` the shape): the starting values put the ends slightly outside of the sample and `c` is the maximum likelihood estimate for them."
  {:params [{:name :a :domain :real :role :lower-bound}
            {:name :b :domain :real :role :upper-bound}
            {:name :c :domain :positive :role :shape}]
   :valid? #(real-data-error % 3)
   :initial (fn [^doubles data]
              (let [h (spacing-of data)
                    a (m/- (Array/aget data 0) h)
                    b (m/+ (Array/aget data (m/dec (alength data))) h)
                    s (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/log (m// (m/- x a) (m/- b a))))) 0.0 data))
                           (alength data))]
                (when (m/neg? s) {:a a :b b :c (m/- (m// s))})))})

(def reciprocal
  "Specification of the reciprocal distribution (`:a` the lower end, `:b` the upper end): the maximum likelihood estimate is the sample minimum and maximum."
  {:params [{:name :a :domain :positive :role :lower-bound}
            {:name :b :domain :positive :role :upper-bound}]
   :valid? #(positive-data-error % 2)
   :estimators {:mle {:exact? true
                      :estimate (fn [^doubles data]
                                  (let [a (Array/aget data 0)
                                        b (Array/aget data (m/dec (alength data)))]
                                    (when (m/< a b) {:a a :b b})))}}})

(def johnson-su
  "Specification of the Johnson SU distribution (`:gamma` and `:delta` shapes, `:xi` location, `:lambda` scale): no starting values, the parameters are found from the scan."
  {:params [{:name :gamma :domain :real :role :shape}
            {:name :delta :domain :positive :role :shape}
            {:name :xi :domain :real :role :location}
            {:name :lambda :domain :positive :role :scale}]
   :valid? #(real-data-error % 4)})

(def johnson-sb
  "Specification of the Johnson SB distribution (`:gamma` and `:delta` shapes, `:xi` the lower end of the support, `:lambda` the width of the support): no starting values, the parameters are found from the scan."
  {:params [{:name :gamma :domain :real :role :shape}
            {:name :delta :domain :positive :role :shape}
            {:name :xi :domain :real :role :lower-bound}
            {:name :lambda :domain :positive :role :scale}]
   :valid? #(real-data-error % 4)})

(def johnson-sl
  "Specification of the Johnson SL distribution (`:gamma` and `:delta` shapes, `:xi` the lower end of the support, `:lambda` the scale). `gamma` and `lambda` can't be told apart (only `ln(lambda) - gamma / delta` matters): fix one of them, for example `{:fixed {:lambda 1.0}}`."
  {:params [{:name :gamma :domain :real :role :shape}
            {:name :delta :domain :positive :role :shape}
            {:name :xi :domain :real :role :lower-bound}
            {:name :lambda :domain :positive :role :scale}]
   :valid? #(real-data-error % 4)})

(def fatigue-life
  "Specification of the fatigue life distribution (`:alpha` the lower end of the support, `:beta` the scale, `:gamma` the shape): the starting values put `alpha` a little below the minimum and use the modified moment estimates of the Birnbaum-Saunders distribution."
  {:params [{:name :alpha :domain :real :role :lower-bound}
            {:name :beta :domain :positive :role :scale}
            {:name :gamma :domain :positive :role :shape}]
   :valid? #(real-data-error % 3)
   :initial (fn [^doubles data]
              (let [alpha (below-minimum data)
                    shifted (map (fn [^double x] (m/- x alpha)) data)
                    s (mean-of shifted)
                    r (m// (mean-of (map m// shifted)))
                    gamma (m/sqrt (m/* 2.0 (m/dec (m/sqrt (m// s r)))))]
                (when (m/valid-double? gamma) {:alpha alpha :beta (m/sqrt (m/* s r)) :gamma gamma})))})

(def hyperbolic-secant
  "Specification of the hyperbolic secant distribution (`:mu`, `:sigma`): the starting values are the sample mean and standard deviation."
  {:params [{:name :mu :domain :real :role :location}
            {:name :sigma :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [data] {:mu (mean-of data) :sigma (sd-of data)})})

(def half-cauchy
  "Specification of the half-Cauchy distribution (`:mu` the lower end of the support, `:scale`): the starting values put `mu` a little below the minimum and the scale is the median distance from it."
  {:params [{:name :mu :domain :real :role :lower-bound}
            {:name :scale :domain :positive :role :scale}]
   :valid? #(real-data-error % 2)
   :initial (fn [^doubles data]
              (let [mu (below-minimum data)]
                {:mu mu :scale (m/- (median-of data) mu)}))})

(def half-logistic
  "Specification of the half-logistic distribution (`:scale`): the starting value comes from the mean, which is `2 ln(2) scale`."
  {:params [{:name :scale :domain :positive :role :scale}]
   :valid? #(non-negative-data-error % 1)
   :initial (fn [data] {:scale (m// (mean-of data) (m/* 2.0 (m/log 2.0)))})})

(def folded-normal
  "Specification of the folded normal distribution (`:mu`, not negative, `:sigma` of the normal distribution before the folding): the starting values are the mean and the standard deviation of the sample."
  {:params [{:name :mu :domain :non-negative :role :location}
            {:name :sigma :domain :positive :role :scale}]
   :valid? #(non-negative-data-error % 2)
   :initial (fn [data] {:mu (mean-of data) :sigma (sd-of data)})})

(def beta-symmetrical
  "Specification of the symmetrical beta distribution (`:alpha`, `:d` is the precision of the computation and is not fitted): the starting value comes from the variance, `1 / (4 (2 alpha + 1))`."
  {:params [{:name :alpha :domain :positive :role :shape}]
   :valid? #(unit-data-error % 2)
   :initial (fn [data]
              (let [alpha (m/* 0.5 (m/dec (m// (m/* 4.0 (central-moment-of data 2)))))]
                (when (m/pos? alpha) {:alpha alpha})))})

(defn- von-mises-concentration
  "Concentration of the von Mises distribution from the mean resultant length `r` (the approximation of Best and Fisher)."
  ^double [^double r]
  (cond
    (m/< r 0.53) (m/+ (m/* 2.0 r) (m/* r r r) (m// (m/* 5.0 (m/fpow r 5)) 6.0))
    (m/< r 0.85) (m/+ -0.4 (m/* 1.39 r) (m// 0.43 (m/- 1.0 r)))
    :else (m// (m/- (m/+ (m/* r r r) (m/* 3.0 r)) (m/* 4.0 r r)))))

(def von-mises
  "Specification of the von Mises distribution (`:mu` the mean direction, `:kappa` the concentration): the maximum likelihood estimate of `mu` is the circular mean, the starting value of `kappa` comes from the mean resultant length. The distribution is not wrapped, its support is `[mu - pi, mu + pi]`: the angles have to lie in an interval of the length `2 pi` which contains this window, move them to the window around the circular mean first."
  {:params [{:name :mu :domain :real :role :location}
            {:name :kappa :domain :non-negative :role :shape}]
   :valid? #(or (real-data-error % 2)
                (when (m/> (m/- (double (reduce m/max %)) (double (reduce m/min %))) (m/* 2.0 m/PI))
                  "angles have to lie in an interval of the length 2 pi"))
   :initial (fn [data]
              (let [n (count data)
                    c (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/cos x))) 0.0 data)) n)
                    s (m// (double (reduce (fn [^double s ^double x] (m/+ s (m/sin x))) 0.0 data)) n)
                    r (m/hypot c s)]
                (when (m/pos? r) {:mu (m/atan2 s c) :kappa (von-mises-concentration (m/min r 0.999))})))})

(defn- ex-gaussian-moments
  "Method of moments of the ex-Gaussian distribution: `tau` from the skewness."
  [data]
  (let [sd (sd-of data)
        skew (m// (central-moment-of data 3) (m/fpow sd 3))
        tau (m/* sd (m/cbrt (m/* 0.5 (m/max skew 0.05))))
        variance (m/- (m/sq sd) (m/sq tau))]
    {:mu (m/- (mean-of data) tau) :sigma (m/sqrt (m/max variance (m/* 0.01 (m/sq sd)))) :tau tau}))

(def ex-gaussian
  "Specification of the ex-Gaussian distribution (`:mu`, `:sigma`, `:tau`): the starting values come from the method of moments."
  {:params [{:name :mu :domain :real :role :location}
            {:name :sigma :domain :positive :role :scale}
            {:name :tau :domain :positive :role :scale}]
   :valid? #(real-data-error % 3)
   :initial ex-gaussian-moments})

(def exgaus
  "Specification of the ex-Gaussian distribution with the gamlss parameters (`:mu`, `:sigma`, `:nu` the mean of the exponential part): the starting values come from the method of moments."
  {:params [{:name :mu :domain :real :role :location}
            {:name :sigma :domain :positive :role :scale}
            {:name :nu :domain :positive :role :scale}]
   :valid? #(real-data-error % 3)
   :initial (fn [data] (let [{:keys [mu sigma tau]} (ex-gaussian-moments data)] {:mu mu :sigma sigma :nu tau}))})

(def pearson-6
  "Specification of the Pearson type VI distribution (`:alpha1` and `:alpha2` shapes, `:beta` the scale): no starting values, the parameters are found from the scan."
  {:params [{:name :alpha1 :domain :positive :role :shape}
            {:name :alpha2 :domain :positive :role :shape}
            {:name :beta :domain :positive :role :scale}]
   :valid? #(positive-data-error % 3)})
