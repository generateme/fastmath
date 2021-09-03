(ns fastmath.optimization.bo
  (:require [fastmath.gp :as gp]
            [fastmath.random :as r]
            [fastmath.kernel :as k]
            [fastmath.core :as m]
            [fastmath.optimization :as opt]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defmulti utility-function :method)

(defmethod utility-function :default [m]
  (utility-function (assoc m :method :ucb)))

(defmethod utility-function :ucb
  [{:keys [^double kappa gp]
    :or {kappa 2.576}}]
  (fn [& x]
    (let [[^double mean ^double stddev] (gp/predict gp x true)]
      (+ mean (* kappa stddev)))))

(defmethod utility-function :ei
  [{:keys [^double y-max gp ^double xi]
    :or {xi 0.001}}]
  (let [d (+ y-max xi)]
    (fn [& x]
      (let [[^double mean ^double stddev] (gp/predict gp x true)
            diff (- mean d)
            z (/ diff stddev)]
        (+ (* diff ^double (r/cdf r/default-normal z))
           (* stddev ^double (r/pdf r/default-normal z)))))))

(defmethod utility-function :poi
  [{:keys [^double y-max gp ^double xi]
    :or {xi 0.001}}]
  (let [d (+ y-max xi)]
    (fn [& x]
      (let [[^double mean ^double stddev] (gp/predict gp x true)]
        (r/cdf r/default-normal (/ (- mean d) stddev))))))

(defmethod utility-function :pi [m]
  (utility-function (assoc m :method :poi)))

(defmethod utility-function :ei-pi
  [{:keys [^double kappa] :as m
    :or {kappa 1.0}}]
  (let [ei (utility-function (assoc m :method :ei))
        pi (utility-function (assoc m :method :poi))]
    (fn [& x]
      (let [^double vei (apply ei x)
            ^double vpi (apply pi x)]
        (+ vpi (* kappa vei))))))

;; https://www.ijcai.org/Proceedings/2020/0316.pdf
(defmethod utility-function :rgp-ucb
  [{:keys [^double scale ^long t] :as m
    :or {scale 1.0 t 0}}]
  (let [t (+ t 2)
        shape (max 0.01 (/ (m/log (/ (inc (* t t)) m/SQRT2PI))
                           (m/log (inc scale))))
        gamma (r/distribution :gamma {:shape shape scale scale})
        beta (m/sqrt (r/drandom gamma))]
    (utility-function (assoc m :kappa beta :method :ucb))))

;;

(defn- cartesian-prod
  [l1 & rst]
  (if (seq rst)
    (for [l l1
          r (apply cartesian-prod rst)]
      (conj r l))
    (map list l1)))

(defn prepare-args
  ([args] (prepare-args args true))
  ([args append?]
   (let [categorical (map first (filter (comp set? second) args))
         continuous (keys (apply dissoc args categorical))
         result (if (seq categorical)
                  (let [categorical-sel (apply juxt categorical)
                        space (map #(vector (zipmap categorical %)
                                            continuous) (apply cartesian-prod (categorical-sel args)))]
                    {::categorical categorical
                     ::continuous continuous
                     ::space space})
                  {::space continuous
                   ::continuous continuous})]
     (if append?
       (merge args result)
       result))))

(defn- ensure-prepared-args
  [args]
  (if (::space args) args (prepare-args args)))

(defn continuous-restrict
  ([args categorical-var categorical-val continuous-vars]
   (continuous-restrict args #(= (% categorical-var) categorical-val) continuous-vars))
  ([args pred continuous-vars]
   (-> args
       (ensure-prepared-args)
       (update ::space #(map (fn [[m c]]
                               (if (pred m)
                                 [m continuous-vars]
                                 [m c])) %)))))

(defn continuous-append
  ([args categorical-var categorical-val continuous-vars]
   (continuous-append args #(= (% categorical-var) categorical-val) continuous-vars))
  ([args pred continuous-vars]
   (-> args
       (ensure-prepared-args)
       (update ::space #(map (fn [[m c]]
                               (if (pred m)
                                 [m (distinct (concat c continuous-vars))]
                                 [m c])) %)))))

(defn- continuous-restrict-or-append-key
  [f args categorical-var deps]
  (reduce (fn [args [cat-val cont-vars]]
            (f args categorical-var cat-val cont-vars)) args deps))

(defn continuous-restrict-key
  [args categorical-var deps]
  (continuous-restrict-or-append-key continuous-restrict args categorical-var deps))

(defn continuous-append-key
  [args categorical-var deps]
  (continuous-restrict-or-append-key continuous-append args categorical-var deps))

(defn categorical-drop
  ([args pred]
   (let [prepared-args (ensure-prepared-args args)]
     (update prepared-args ::space #(remove (comp pred first) %))))
  ([args k pred]
   (let [vs (-> args
                (ensure-prepared-args)
                (update ::space #(->> %
                                      (map (fn [[m c]]
                                             (if (pred m)
                                               [(dissoc m k) c]
                                               [m c])))
                                      (group-by first)
                                      (map (fn [[m lst]]
                                             [m (->> (map second lst)
                                                     (mapcat identity)
                                                     (distinct))])))))]
     vs)))

(defn categorical-keep
  ([args pred]
   (categorical-drop args (complement pred)))
  ([args k pred]
   (categorical-drop args k (complement pred))))


;; (def argg {:a #{:a :b}
;;            :b #{:x :z}
;;            :c #{1 2}
;;            :z [1 2]
;;            :w [9 6]})

;; (-> (prepare-args argg)
;;     (continuous-restrict #(and (= (:a %) :b) (= (:b %) :x)) [:z])
;;     (categorical-drop :c #(not= (:b %) :z)))


(require '[fastmath.clustering :as clust])

(def example-data
  (shuffle (concat
            (r/->seq (r/distribution :weibull {:alpha 3 :beta 5}) 1000)
            (r/->seq (r/distribution :gamma {:shape 10.0}) 1000))))

(defn target
  [{:keys [clustering d1 d2] :as args}]
  (let [clusters (clust/regroup ((get {:x-means clust/x-means
                                       :clarans clust/clarans
                                       :deterministic-annealing clust/deterministic-annealing} clustering) example-data 2))
        distr1 (r/distribution d1 args)
        distr2 (r/distribution d2 args)]
    (+ (r/log-likelihood distr1 (:data (first clusters)))
       (r/log-likelihood distr2 (:data (second clusters))))))

(target {:clustering :x-means
         :d1 :gamma
         :d2 :weibull
         :alpha 3.0
         :beta 5.0
         :shape 10.0
         :scale 4.0
         :k 2.0
         :lambda 2.0})

(def args (-> {:clustering #{:x-means :clarans :deterministic-annealing}
               :d1 #{:gamma :weibull :erlang}
               :d2 #{:gamma :weibull :erlang}
               :dist #{:euclidean :manhattan :chebyshev}
               :alpha [0.01 10]
               :beta [0.01 10]
               :shape [1.0 20.0]
               :scale [1.0 10.0]
               :k [1.0 10.0]
               :lambda [0.1 10.0]}
              (continuous-restrict-key :d1 {:erlang [:k :lambda]
                                            :weibull [:alpha :beta]
                                            :gamma [:shape :scale]}) ;; replaces continuous list with provided for given pair key/val from categorical
              (continuous-append-key :d2 {:erlang [:k :lambda]
                                          :weibull [:alpha :beta]
                                          :gamma [:shape :scale]})
              (categorical-drop #(= (:d1 %) (:d2 %))) ;; removes categorical when pred is true
              (categorical-keep :dist #(= (:clustering %) :clarans)) ;; keeps given key only when pred is true
              (categorical-drop #(and (= (:d1 %) :weibull)
                                      (not= (:clustering %) :x-means)))))

(def args-cont {:alpha [0.01 10]
                :beta [0.01 10]
                :shape [1.0 20.0]
                :scale [1.0 10.0]
                :k [1.0 10.0]
                :lambda [0.1 10.0]})

(defn- make-forward-backward
  [args ks forward?]
  (into {} (map (fn [k]
                  (let [[mn mx] (args k)]
                    [k (if forward?
                         (m/make-norm mn mx 0.0 1.0)
                         (m/make-norm 0.0 1.0 mn mx))])) ks)))

(defn- make-x-normalizer
  [args]
  (let [ks (::continuous args)
        forward (make-forward-backward args ks true)
        backward (make-forward-backward args ks false)]
    (fn local-normalizer
      ([m] (local-normalizer m true))
      ([m forward?]
       (let [fm (if forward? forward backward)]
         (reduce (fn [curr k]
                   (update curr k (fm k))) m (keys m)))))))

(defn initialize-bo
  ([args] (initialize-bo args {}))
  ([args {:keys [noise normalize-y? normalize-x? utility-function-type utility-function-params
                 kernel kernel-params kernel-scale
                 optimizer optimizer-params
                 jitter]
          :or {noise 1.0e-5
               normalize-y? true
               normalize-x? true
               utility-function-type :rgp-ucb
               kernel :mattern-12
               optimizer :bfgs
               jitter 0.1}}]
   (let [args (ensure-prepared-args args)
         bandit? (::categorical args)
         x-normalizer (if-not normalize-x?
                        (fn [m & _] m)
                        (make-x-normalizer args))]
     {:args args
      :noise noise
      :normalize-y? normalize-y?
      :x-normalizer x-normalizer
      :bandit? bandit?
      :utility-function-type utility-function-type
      :utility-function-params utility-function-params
      :kernel kernel
      :kernel-params kernel-params
      :kernel-scale kernel-scale
      :optimizer optimizer
      :optimizer-params optimizer-params
      :jitter jitter})))

(comment
  (def gp (gp/gaussian-process [[1 1] [2 2] [3 3] [9 -1] [5 -2]] [-2 -1 1 -1 2] {:noise 0.001
                                                                                 :normalize? true
                                                                                 :kernel (k/kernel :mattern-52 1.48796)
                                                                                 :kscale 0.87396}))

  (def gp (gp/gaussian-process [[1 1]] [-2] {:noise 0.001
                                             :normalize? true
                                             :kernel (k/kernel :mattern-52)
                                             :kscale 0.87396}))

  (gp/predict gp '(4.62100181160906 -3.5))

  (opt/scan-and-maximize :bfgs (utility-function {:method :ucb

                                                  :gp gp}) {:bounds [[0 10] [-5 5]]})
  ;; => [(4.862100181160906 -2.349996775981009) 0.557154163177696]
  ;; => [(4.6001465276639495 -3.531503947537969) 4.333240536082748]
  ;; => [(4.981048745997669 -2.0285429857231874) 0.4899829648551593]
  ;; => [(4.750949198693047 -2.818882138909999) 0.19133431735709724]
  ;; => [(4.567607972440291 -3.6893979394097896) 4.945040161917608]
  ;; => [(4.600145774519866 -3.5315036892286065) 4.333240536082746]

  (time (opt/scan-and-maximize :bfgs (fn [^double k ^double s]
                                       ;; (println k)
                                       (let [gp (gp/gaussian-process [[1 1] [2 2] [3 3] [9 -1] [5 -2]] [-2 -1 1 -1 2]
                                                                     {:noise 0.0001
                                                                      :normalize? true
                                                                      :kernel (k/kernel :mattern-52 k)
                                                                      :kscale s})]
                                         (gp/L gp))) {:bounds [[0.001 20] [0.001 20]]
                                                      :initial [3.0]}))

  (gp/L gp)


  (gp/predict gp [2 22] true)
  (gp/predict-all gp [[1 1] [2 2]] true)

  )
