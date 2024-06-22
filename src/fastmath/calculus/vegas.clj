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


;; vegas / vegas+
;; https://arxiv.org/pdf/2009.05112.pdf
;; https://github.com/ranjanan/MonteCarloIntegration.jl/blob/master/src/vegas.jl

(defn vegas-build-initial-grid
  [lower upper ^long nintervals]
  (let [pairs (map vector lower upper)]
    [(mapv (fn [[^double lb ^double ub]]
             (m/seq->double-array (m/slice-range lb ub (m/inc nintervals)))) pairs)
     (mapv (fn [[^double lb ^double ub]]
             (m/seq->double-array (repeat nintervals (m// (m/- ub lb) nintervals)))) pairs)]))

(defn vegas-random-sequence
  [^long dims {:keys [random-sequence ^double jitter]
               :or {random-sequence :uniform jitter 0.75}}]
  (if (and (not (= random-sequence :uniform))
           (m/pos? jitter))
    (r/jittered-sequence-generator random-sequence dims jitter)
    (r/sequence-generator random-sequence dims)))

(defn Itot-sd
  [integrals rev-sigma-squares]
  (let [sum-rev-sigma-squares (v/sum rev-sigma-squares)]
    [(m// (v/sum (v/emult integrals rev-sigma-squares)) sum-rev-sigma-squares)
     (m// (m/sqrt sum-rev-sigma-squares))]))

(defn normalize-and-compress-d!
  [^doubles d ^double sumd ^double alpha]
  (dotimes [i (alength d)]
    (let [v (m// (aget d i) sumd)]
      (aset d i (m/pow (m// (m/- 1.0 v)
                            (m/log (m// v))) alpha)))))

(defn smooth-d
  [^doubles d ^long nintervals ^double alpha]
  (let [nintervals-1 (dec nintervals)
        buff (double-array nintervals)
        sumd (m/* 8.0 (v/sum d))]
    
    (aset buff 0 (m/+ (m/* 7.0 (aget d 0))
                      (aget d 1)))
    (doseq [^long id (range 1 nintervals-1)]
      (aset buff id (m/+ (aget d (dec id))
                         (m/* 6.0 (aget d id))
                         (aget d (inc id)))))
    (aset buff nintervals-1 (m/+ (aget d (dec nintervals-1))
                                 (m/* 7.0 (aget d nintervals-1))))
    
    (normalize-and-compress-d! buff sumd alpha)

    buff))

(defn calculate-d
  [dims Jsf ssamples nintervals alpha]
  (let [nevals (count Jsf)
        ni (m// nevals (double nintervals))
        Jsf22 (mapv (fn [idxs ^double v]
                      [(int-array idxs) (m// (m/sq v) ni)]) ssamples Jsf)]
    (->> (for [dim (range dims)]
           (let [buffer (double-array nintervals)]
             (run! (fn [in]
                     (let [id (aget ^ints (in 0) dim)]
                       (aset buffer id (m/+ (aget buffer id)
                                            (double (in 1)))))) Jsf22)
             buffer))
         (mapv #(smooth-d % nintervals alpha)))))

(defn update-grid-dim [^doubles x ^doubles dx ^doubles d ^long nintervals]
  (let [cx (count x)
        end (m/dec cx)
        buff (double-array cx)
        delta-d (stats/mean d)]
    (aset buff 0 (aget x 0))
    (aset buff end (aget x end))
    (loop [i (long 1)
           j (long 0)
           sd 0.0]
      (when (m/< i nintervals)
        (let [^Vec2 nj-nsd (loop [nsd sd
                                  nj j]
                             (if (m/< nsd delta-d)
                               (recur (m/+ nsd (aget d nj))
                                      (m/inc nj))
                               (Vec2. nj nsd)))
              nj (long (.x nj-nsd))
              nsd (.y nj-nsd)
              nsd (m/- nsd delta-d)
              nj- (m/dec nj)]
          (aset buff i (m/- (aget x nj)
                            (m// (m/* nsd (aget dx nj-))
                                 (aget d nj-))))
          (recur (m/inc i) (long nj) nsd))))
    buff))

(defn calc-dx [x]
  (m/seq->double-array (map (fn [[^double a ^double b]]
                              (m/- b a)) (partition 2 1 x))))

(defn update-grid [xs dxs ds ^long nintervals]
  (let [nxs (mapv (fn [x dx d]
                    (update-grid-dim x dx d nintervals)) xs dxs ds)
        ndxs (mapv calc-dx nxs)]
    [nxs ndxs]))

;; strata

(defn stratifications-number
  "Calculate optimal number of stratas per dimension, eq 41"
  ^long [^long nev ^long dims]
  (-> (m// nev 4.0)
      (m/pow (m// dims))
      (m/floor)
      (unchecked-int)
      (m/max 1)))

(defn ensure-bins-count
  "Make grid size divisible by number of stratas, comment to eq 41"
  ^long [^long nintervals ^long nst]
  (if (m/zero? (m/mod nintervals nst))
    nintervals
    (m/* (m/inc (m// nintervals nst)) nst)))

(defn reshape-hypercube [hc] (apply mapv vector hc))

(defn volume
  ^double [[left right]]
  (->> (map m/- right left)
       (reduce m/*)))

(defn build-hypercubes
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
                         (->> (m/slice-range 0.0 1.0 (inc n))
                              (partition 2 1)))
                       nstrats)
        
        hypercubes (->> strata-1d
                        (apply combo/cartesian-product)
                        (map reshape-hypercube))
        
        nhcubes (count hypercubes)
        nhs (repeat nhcubes (-> (m// nevals (double nhcubes))
                                (m/floor)
                                (unchecked-int)
                                (m/+ 2)))]
    {:nstrats nstrats
     :nintervals nintervals
     :nhs nhs
     :hcubes hypercubes
     :nhcubes nhcubes
     :volume (volume (first hypercubes))}))

(defn random-sequence-samples
  [random-sequence hcubes nhs ^long nhcubes]
  (if (m/one? nhcubes)
    (split-at (first nhs) random-sequence)
    (loop [[nh & rest-nhs] nhs
           [[left right] & rest-hcubes] hcubes
           rs random-sequence
           res []]
      (if left
        (let [[samples rest-samples] (split-at nh rs)
              interpolated (map (partial v/einterpolate left right) samples)              ]
          (recur rest-nhs rest-hcubes rest-samples (conj res interpolated)))
        [(mapcat identity res) rs]))))

(defn sigmas-in-hypercubes-and-integral-mc
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
            nsigmas (conj sigmas (m/+ (m/abs (m/- (m// fdv2-sum nh) (m/* I I)))
                                      m/EPSILON))]
        (recur rst nhrest nsigmas (m/+ Imc I)))
      [sigmas Imc])))

(defn recalculate-nhs
  [hcsigmas ^long nevals ^double hbeta]
  (let [dhs (->> hcsigmas
                 (map (fn [^double s] (m/pow s hbeta))))
        dhssum (v/sum dhs)]
    (map (fn [^double dh]
           (unchecked-int (m/+ (m/floor (m/* nevals (m// dh dhssum))) 2))) dhs)))

;;

(defn vegas
  ([f lower upper] (vegas f lower upper nil))
  ([f lower upper {:keys [^long max-iters ^double rel ^double abs ^long nevals
                          alpha ^double beta ^long warmup info? record-data?]
                   :or {max-iters 10 rel 5.0e-4 abs 5.0e-4 nevals 10000 beta 0.75 warmup 0
                        info? false record-data? false}
                   :as options}]
   (let [[f lower upper] (cc/subst-multi f lower upper)
         dims (count lower)
         max-iters (m/+ (m/max 1 max-iters)
                        (m/max 0 warmup))
         
         random-sequence (vegas-random-sequence dims options)
         
         {:keys [^long nintervals nstrats ^long nhcubes
                 hcubes nhs ^double volume]} (build-hypercubes (assoc options
                                                                      :dims dims
                                                                      :nevals nevals))

         ;; lower damping alpha for statified sampling
         alpha (double (or alpha (if (and (= (count nstrats) 1)
                                          (= (first nstrats) 1)) 1.5 0.5)))
         hbeta (m// beta 2.0)]
     
     (loop [iter (long 1)
            evaluations (long 1)
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
                                                           (m/+ (Array/get ^doubles ax iy)
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
                               (aset ^doubles buff dim (m/+ (Array/get ^doubles (x dim) iy)
                                                            (m/* (Array/get adx iy)
                                                                 (m/frac ny))))
                               (recur (m/inc dim)
                                      ;; calculate Jacobian
                                      (m/* Js nintervals (Array/get adx iy))))
                             (m/* Js (double (f buff))))))
                       scaled-samples)
             
             real-calls (count Jsf)
             evaluations (+ evaluations real-calls)
             
             [hcsigmas integral-mc] (sigmas-in-hypercubes-and-integral-mc Jsf nhs volume)
             nhs (recalculate-nhs hcsigmas nevals hbeta)             ]
         
         (if (m/<= iter warmup)
           
           (let [d (calculate-d dims Jsf scaled-samples nintervals alpha)
                 x-dx (update-grid x dx d nintervals)]
             (recur (m/inc iter) evaluations
                    x-dx rst nhs integrals rev-sigma-squares ndata))

           (let [variance-mc (m/+ (v/sum (map (fn [^double s ^double nh]
                                                (m// s nh)) hcsigmas nhs)))
                 new-integrals (conj integrals integral-mc)
                 new-rev-sigma-squares (conj rev-sigma-squares (m// variance-mc))
                 
                 [^double Itot ^double sd] (Itot-sd new-integrals new-rev-sigma-squares)]
             
             (if (or (m/== iter max-iters)
                     (and (m/< (m/abs (m// sd Itot)) rel)
                          (m/< (m/abs sd) abs)))
               (if info?
                 (let [chisq (-> new-integrals
                                 (v/shift (m/- Itot))
                                 (v/sq)
                                 (v/emult new-rev-sigma-squares)
                                 (v/sum))
                       dof (m/dec (count new-integrals))
                       res {:iterations (m/- iter warmup)
                            :result Itot
                            :sd sd
                            :nintervals nintervals
                            :nstrats nstrats
                            :nhcubes nhcubes
                            :evaluations evaluations
                            :chi2-avg (/ chisq dof)
                            :dof dof
                            :Q (special/regularized-gamma-q (m// dof 2.0) (m// chisq 2.0))}]
                   (if record-data?
                     (assoc res :data ndata :hcubes hcubes)
                     res))
                 Itot)

               (let [d (calculate-d dims Jsf scaled-samples nintervals alpha)
                     x-dx (update-grid x dx d nintervals)]
                 (recur (m/inc iter) evaluations
                        x-dx rst nhs new-integrals new-rev-sigma-squares ndata))))))))))

