(ns fastmath.calculus.vegas
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.calculus.common :as cc]
            [fastmath.special :as special]
            [clojure.math.combinatorics :as combo])
  (:import [fastmath.vector Vec2]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; vegas / vegas+
;; https://arxiv.org/pdf/2009.05112.pdf
;; https://github.com/ranjanan/MonteCarloIntegration.jl/blob/master/src/vegas.jl

;; a variance of a single iteration is never below this value: weights of the iterations stay finite for a constant integrand
(def ^{:private true :const true :tag 'double} tiny-variance 1.0e-250)

(defn- vegas-build-initial-grid
  [lower upper ^long nintervals]
  (let [pairs (map vector lower upper)]
    [(mapv (fn [[^double lb ^double ub]]
             (m/seq->double-array (m/slice-range lb ub (m/long-inc nintervals)))) pairs)
     (mapv (fn [[^double lb ^double ub]]
             (m/seq->double-array (repeat nintervals (m// (m/- ub lb) nintervals)))) pairs)]))

(defn- vegas-random-sequence
  "Infinite sequence of points of the unit cube (vectors, also for one dimension)."
  [^long dims {:keys [random-sequence ^double jitter rng]
               :or {random-sequence :uniform jitter 0.75}}]
  (let [sequence (if (and (not (= random-sequence :uniform))
                          (m/pos? jitter))
                   (if (nil? rng)
                     (r/jittered-sequence-generator random-sequence dims jitter)
                     (r/jittered-sequence-generator random-sequence dims jitter rng))
                   (if (nil? rng)
                     (r/sequence-generator random-sequence dims)
                     (r/sequence-generator random-sequence dims rng)))]
    ;; the generators return numbers in one dimension
    (if (m/one? dims) (map vector sequence) sequence)))

(defn- Itot-sd
  "Weighted mean of the results of the iterations and its standard deviation."
  [integrals rev-sigma-squares]
  (let [sum-rev-sigma-squares (v/sum rev-sigma-squares)]
    [(m// (v/sum (v/emult integrals rev-sigma-squares)) sum-rev-sigma-squares)
     (m// (m/sqrt sum-rev-sigma-squares))]))

(defn- normalize-and-compress-d!
  [^doubles d ^double sumd ^double alpha]
  (dotimes [i (alength d)]
    (let [v (m// (Array/aget d i) sumd)
          ^double c (cond
                      (m/>= v 1.0) 1.0 ;; limit of (1-v)/ln(1/v)
                      (m/not-pos? v) 0.0
                      :else (m/pow (m// (m/- 1.0 v)
                                        (m/log (m// v))) alpha))]
      (Array/aset d i c))))

(defn- smooth-d
  "Smooths the weights of the grid intervals and compresses them; all zeros (the integrand is zero) are left as they are."
  [^doubles d ^long nintervals ^double alpha]
  (let [nintervals-1 (m/long-dec nintervals)
        buff (double-array nintervals)
        sumd (m/* 8.0 (v/sum d))]
    (if (m/== nintervals 1)
      (Array/aset buff 0 (Array/aget d 0))
      (do
        (Array/aset buff 0 (m/+ (m/* 7.0 (Array/aget d 0))
                                (Array/aget d 1)))
        (doseq [^long id (range 1 nintervals-1)]
          (Array/aset buff id (m/+ (Array/aget d (m/long-dec id))
                                   (m/* 6.0 (Array/aget d id))
                                   (Array/aget d (m/long-inc id)))))
        (Array/aset buff nintervals-1 (m/+ (Array/aget d (m/long-dec nintervals-1))
                                           (m/* 7.0 (Array/aget d nintervals-1))))))
    (when (m/pos? sumd)
      (normalize-and-compress-d! buff sumd alpha))
    buff))

(defn- calculate-d
  [dims Jsf ssamples nintervals alpha]
  (let [dims (long dims)
        nintervals (long nintervals)
        alpha (double alpha)
        nevals (count Jsf)
        ni (m// nevals (double nintervals))
        Jsf22 (mapv (fn [idxs ^double v]
                      [(int-array idxs) (m// (m/sq v) ni)]) ssamples Jsf)]
    (->> (for [dim (range dims)]
           (let [buffer (double-array nintervals)]
             (run! (fn [in]
                     (let [id (aget ^ints (in 0) dim)]
                       (Array/aset buffer id (m/+ (Array/aget buffer id)
                                                  (double (in 1)))))) Jsf22)
             buffer))
         (mapv #(smooth-d % nintervals alpha)))))

(defn- update-grid-dim
  "New grid edges of one dimension, the grid is not changed when all the weights are zero."
  [^doubles x ^doubles dx ^doubles d ^long nintervals]
  (let [cx (alength x)
        end (m/long-dec cx)
        buff (double-array cx)
        delta-d (stats/mean d)]
    (if-not (m/pos? delta-d)
      (aclone x)
      (do
        (Array/aset buff 0 (Array/aget x 0))
        (Array/aset buff end (Array/aget x end))
        (loop [i (long 1)
               j (long 0)
               sd 0.0]
          (when (m/< i nintervals)
            (let [^Vec2 nj-nsd (loop [nsd sd
                                      nj j]
                                 (if (m/< nsd delta-d)
                                   (recur (m/+ nsd (Array/aget d nj))
                                          (m/long-inc nj))
                                   (Vec2. nj nsd)))
                  nj (long (.x nj-nsd))
                  nsd (m/- (.y nj-nsd) delta-d)
                  nj- (m/long-dec nj)]
              (Array/aset buff i (m/- (Array/aget x nj)
                                      (m// (m/* nsd (Array/aget dx nj-))
                                           (Array/aget d nj-))))
              (recur (m/long-inc i) nj nsd))))
        buff))))

(defn- calc-dx
  [x]
  (m/seq->double-array (map (fn [[^double a ^double b]]
                              (m/- b a)) (partition 2 1 x))))

(defn- update-grid
  [xs dxs ds ^long nintervals]
  (let [nxs (mapv (fn [x dx d]
                    (update-grid-dim x dx d nintervals)) xs dxs ds)
        ndxs (mapv calc-dx nxs)]
    [nxs ndxs]))

;; strata

(defn- stratifications-number
  "Calculate optimal number of stratas per dimension, eq 41"
  ^long [^long nev ^long dims]
  (-> (m// nev 4.0)
      (m/pow (m// dims))
      (m/floor)
      (unchecked-long)
      (m/max 1)))

(defn- ensure-bins-count
  "Make grid size divisible by number of stratas, comment to eq 41"
  ^long [^long nintervals ^long nst]
  (if (m/zero? (m/mod nintervals nst))
    nintervals
    (m/long-mult (m/long-inc (m/long-quot nintervals nst)) nst)))

(defn- reshape-hypercube [hc] (apply mapv vector hc))

(defn- volume
  ^double [[left right]]
  (->> (map m/- right left)
       (v/prod)))

(defn- uniform-allocation
  "Number of samples of each hypercube when they are equal."
  ^long [^long nevals ^long nhcubes]
  (m/long-add (long (m/floor (m// nevals (double nhcubes)))) 2))

(defn- build-hypercubes
  [{:keys [nstrats ^long nevals ^long nintervals ^long dims]
    :or {nstrats -1 nintervals 1000}}]
  (let [nstrats (->> (cond
                       (= -1 nstrats) [(stratifications-number nevals dims)]
                       (number? nstrats) [(m/max 1 (unchecked-long nstrats))]
                       (sequential? nstrats) (map unchecked-long nstrats)
                       :else [(stratifications-number nevals dims)])
                     (cycle)
                     (take dims))

        nintervals (->> nstrats
                        (reduce m/lcm)
                        (ensure-bins-count nintervals))

        strata-1d (map (fn [^long n]
                         (->> (m/slice-range 0.0 1.0 (m/long-inc n))
                              (partition 2 1)))
                       nstrats)

        hypercubes (->> strata-1d
                        (apply combo/cartesian-product)
                        (map reshape-hypercube))

        nhcubes (count hypercubes)
        nhs (repeat nhcubes (uniform-allocation nevals nhcubes))]
    {:nstrats nstrats
     :nintervals nintervals
     :nhs nhs
     :hcubes hypercubes
     :nhcubes nhcubes
     :volume (volume (first hypercubes))}))

(defn- random-sequence-samples
  [random-sequence hcubes nhs ^long nhcubes]
  (if (m/one? nhcubes)
    (split-at (first nhs) random-sequence)
    (loop [[nh & rest-nhs] nhs
           [[left right] & rest-hcubes] hcubes
           rs random-sequence
           res []]
      (if left
        (let [[samples rest-samples] (split-at nh rs)
              interpolated (map (partial v/einterpolate left right) samples)]
          (recur rest-nhs rest-hcubes rest-samples (conj res interpolated)))
        [(mapcat identity res) rs]))))

(defn- sigmas-in-hypercubes-and-integral-mc
  "Unbiased variances of the integrand (times the hypercube volume) in every hypercube
  and the Monte Carlo estimate of the integral."
  [Jsf nhs ^double volume]
  (loop [jsf Jsf
         [^long nh & nhrest] nhs
         sigmas []
         Imc 0.0]
    (if (seq jsf)
      (let [[jsf0 rst] (split-at nh jsf)
            fdv (map (fn [^double jsf] (m/* jsf volume)) jsf0)
            fdv-sum (v/sum fdv)
            fdv2-sum (v/sum (map m/sq fdv))
            I (m// fdv-sum nh)
            ;; every hypercube has at least 2 samples; the absolute value removes the rounding noise below zero
            nsigmas (conj sigmas (m/* (m// (double nh) (double (m/long-dec nh)))
                                      (m/abs (m/- (m// fdv2-sum nh) (m/* I I)))))]
        (recur rst nhrest nsigmas (m/+ Imc I)))
      [sigmas Imc])))

(defn- recalculate-nhs
  "Allocation of samples for the next iteration, proportional to a power of the standard deviation in a hypercube;
  uniform when the integrand is constant (no variance anywhere)."
  [hcsigmas ^long nevals ^double hbeta]
  (let [dhs (map (fn [^double s] (m/pow s hbeta)) hcsigmas)
        dhssum (v/sum dhs)]
    (if (and (m/pos? dhssum) (m/valid-double? dhssum))
      (map (fn [^double dh]
             (m/long-add (long (m/floor (m/* nevals (m// dh dhssum)))) 2)) dhs)
      (repeat (count hcsigmas) (uniform-allocation nevals (count hcsigmas))))))

;; options

(defn- option-error
  [option value expected]
  (ex-info (str "Invalid option " option) {:option option :value value :expected expected}))

(defn- natural? [^long min-value] (fn [x] (and (integer? x) (m/>= (long x) min-value))))
(defn- non-negative? [x] (and (number? x) (m/not-neg? (double x))))

(def ^:private option-checks
  [[:nevals (natural? 2) "integer >= 2"]
   [:max-iters (natural? 1) "integer >= 1"]
   [:warmup (natural? 0) "integer >= 0"]
   [:nintervals (natural? 1) "integer >= 1"]
   [:rel non-negative? "number >= 0"]
   [:abs non-negative? "number >= 0"]
   [:alpha non-negative? "number >= 0"]
   [:beta non-negative? "number >= 0"]
   [:random-sequence #{:uniform :r2 :sobol :halton} "one of :uniform, :r2, :sobol, :halton"]])

(defn- check-options
  "Throws `ex-info` with `{:option :value :expected}` for the first invalid option value, unknown options are ignored."
  [options]
  (doseq [[option valid? expected] option-checks]
    (when-some [value (get options option)]
      (when-not (valid? value)
        (throw (option-error option value expected))))))

(defn- checked-integrand
  "Wraps the integrand so that a NaN or an infinite value throws `ex-info` with the point."
  [f]
  (fn ^double [x]
    (let [y (double (f x))]
      (if (m/valid-double? y)
        y
        (throw (ex-info "The integrand returned a NaN or an infinite value" {:x (vec x) :value y}))))))

(defn- zero-width-result
  [info?]
  (if info?
    {:result 0.0 :sd 0.0 :iterations 0 :evaluations 0 :nintervals 0 :nstrats [] :nhcubes 0
     :chi2-avg ##NaN :dof 0 :Q ##NaN}
    0.0))

;;

(defn- vegas-iterations
  "The iterations of the algorithm on the oriented bounds."
  [f lower upper options]
  (let [{:keys [^long max-iters ^double rel ^double abs ^long nevals alpha ^double beta ^long warmup info? record-data?]
         :or {max-iters 10 rel 5.0e-4 abs 5.0e-4 nevals 10000 beta 0.75 warmup 0}} options
        [f lower upper] (cc/subst-multi (checked-integrand f) lower upper)
        dims (count lower)
        max-iters (m/+ max-iters warmup)

        random-sequence (vegas-random-sequence dims options)

        {:keys [^long nintervals nstrats ^long nhcubes
                hcubes nhs ^double volume]} (build-hypercubes (assoc options
                                                                     :dims dims
                                                                     :nevals nevals))

        ;; lower damping alpha for stratified sampling; the original VEGAS (one stratum per dimension) refines the grid faster
        alpha (double (or alpha (if (every? m/one? nstrats) 1.5 0.5)))
        hbeta (m// beta 2.0)]

    (loop [iter (long 1)
           evaluations (long 0)
           [x dx] (vegas-build-initial-grid lower upper nintervals)
           random-sequence random-sequence
           nhs nhs
           integrals []
           rev-sigma-squares []
           data []]
      (let [[samples rst] (random-sequence-samples random-sequence hcubes nhs nhcubes)

            scaled-samples (mapv #(v/mult % nintervals) samples)

            ndata (when (and record-data? info?)
                    (conj data {:x x :dx dx :nhs nhs
                                :samples (map (fn [ys]
                                                (mapv (fn [^doubles ax ^doubles adx ^double ny]
                                                        (let [iy (unchecked-int ny)]
                                                          (m/+ (Array/get ax iy)
                                                               (m/* (Array/get adx iy)
                                                                    (m/frac ny))))) x dx ys)) scaled-samples)}))

            buff (double-array dims)
            Jsf (mapv (fn [ys]
                        (loop [dim (long 0)
                               Js 1.0]
                          (if (m/< dim dims)
                            (let [ny (double (ys dim))
                                  iy (unchecked-int ny)
                                  ^doubles adx (dx dim)]
                              ;; write mapped y->x to a buffer
                              (Array/aset buff dim (m/+ (Array/get ^doubles (x dim) iy)
                                                        (m/* (Array/get adx iy)
                                                             (m/frac ny))))
                              (recur (m/long-inc dim)
                                     ;; calculate Jacobian
                                     (m/* Js nintervals (Array/get adx iy))))
                            (m/* Js (double (f buff))))))
                      scaled-samples)

            evaluations (m/+ evaluations (long (count Jsf)))

            ;; variances of this iteration belong to the allocation which was sampled, not to the next one
            [hcsigmas ^double integral-mc] (sigmas-in-hypercubes-and-integral-mc Jsf nhs volume)
            next-nhs (recalculate-nhs hcsigmas nevals hbeta)]

        (if (m/<= iter warmup)

          (let [d (calculate-d dims Jsf scaled-samples nintervals alpha)
                x-dx (update-grid x dx d nintervals)]
            (recur (m/long-inc iter) evaluations
                   x-dx rst next-nhs integrals rev-sigma-squares ndata))

          (let [variance-mc (m/max tiny-variance
                                   (v/sum (map (fn [^double s ^double nh]
                                                 (m// s nh)) hcsigmas nhs)))
                new-integrals (conj integrals integral-mc)
                new-rev-sigma-squares (conj rev-sigma-squares (m// variance-mc))

                [^double Itot ^double sd] (Itot-sd new-integrals new-rev-sigma-squares)]

            (if (or (m/== iter max-iters)
                    (m/<= sd (m/max (m/* rel (m/abs Itot)) abs)))
              (if info?
                (let [chisq (-> new-integrals
                                (v/shift (m/- Itot))
                                (v/sq)
                                (v/emult new-rev-sigma-squares)
                                (v/sum))
                      dof (m/long-dec (count new-integrals))
                      res {:iterations (m/- iter warmup)
                           :result Itot
                           :sd sd
                           :nintervals nintervals
                           :nstrats nstrats
                           :nhcubes nhcubes
                           :evaluations evaluations
                           :chi2-avg (if (m/zero? dof) ##NaN (m// chisq dof))
                           :dof dof
                           :Q (if (m/zero? dof) ##NaN (special/regularized-gamma-q (m// dof 2.0) (m// chisq 2.0)))}]
                  (if record-data?
                    (assoc res :data ndata :hcubes hcubes)
                    res))
                Itot)

              (let [d (calculate-d dims Jsf scaled-samples nintervals alpha)
                    x-dx (update-grid x dx d nintervals)]
                (recur (m/long-inc iter) evaluations
                       x-dx rst next-nhs new-integrals new-rev-sigma-squares ndata)))))))))

(defn vegas
  "VEGAS / VEGAS+ Monte Carlo integration of `f` on the box `[lower, upper]`, see [[fastmath.calculus/vegas]] for options and result."
  ([f lower upper] (vegas f lower upper nil))
  ([f lower upper options]
   (let [options (cc/drop-nil-options options)
         [sign lower upper] (cc/normalize-bounds-multi lower upper)]
     (check-options options)
     (cc/orient sign
                (if (some true? (map m/== lower upper))
                  (zero-width-result (:info? options))
                  (vegas-iterations f lower upper options))))))
