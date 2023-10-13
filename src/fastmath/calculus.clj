(ns fastmath.calculus
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; https://arxiv.org/pdf/2009.05112.pdf
;; https://github.com/ranjanan/MonteCarloIntegration.jl/blob/master/src/vegas.jl

;; vegas / vegas+

(defn- vegas-build-initial-grid
  [bounds ^long nbins]
  [(mapv (fn [[^double lb ^double ub]]
           (vec (m/slice-range lb ub (m/inc nbins)))) bounds)
   (mapv (fn [[^double lb ^double ub]]
           (vec (repeat nbins (m// (m/- ub lb) nbins)))) bounds)])

(defn- vegas-random-sequence
  [^long dims {:keys [random-sequence ^double jitter]
               :or {random-sequence :uniform jitter 0.0}}]
  (if (and (not (= random-sequence :uniform))
           (m/pos? jitter))
    (r/jittered-sequence-generator random-sequence dims jitter)
    (r/sequence-generator random-sequence dims)))

(defn- Itot-sd
  [integrals rev-sigma-squares]
  (let [sum-rev-sigma-squares (v/sum rev-sigma-squares)]
    [(m// (v/sum (v/emult integrals rev-sigma-squares)) sum-rev-sigma-squares)
     (m// (m/sqrt sum-rev-sigma-squares))]))

(defn- normalize-and-compress-d
  [d ^double sumd ^double alpha]
  (mapv (fn [^double v]
          (let [d (m// v sumd)]
            (m/pow (m// (m/- 1.0 d)
                        (m/log (m// d))) alpha))) d))

(defn- smooth-d
  [^doubles d ^long nbins ^double alpha]
  (let [nbins-1 (dec nbins)
        buff (double-array nbins)
        sumd (* 8.0 (stats/sum d))]
    
    (aset buff 0 (m/+ (m/* 7.0 (aget d 0))
                      (aget d 1)))
    (doseq [^long id (range 1 nbins-1)]
      (aset buff id (m/+ (aget d (dec id))
                         (* 6.0 (aget d id))
                         (aget d (inc id)))))
    (aset buff nbins-1 (m/+ (aget d (dec nbins-1))
                            (m/* 7.0 (aget d nbins-1))))
    
    (normalize-and-compress-d buff sumd alpha)))

(defn- calculate-d
  [dims Jsf imat nbins alpha]
  (let [ncalls (count Jsf)
        ni (m// ncalls (double nbins))
        Jsf22 (map (fn [idxs ^double v]
                     [idxs (m// (m/sq v) ni)]) imat Jsf)]
    (for [dim (range dims)]
      (let [buffer (double-array nbins)]
        (run! (fn [[idxs ^double v]]
                (let [id (idxs dim)
                      curr (aget buffer id)]
                  (aset buffer id (m/+ curr v)))) Jsf22)
        (smooth-d buffer nbins alpha)))))

(defn- update-grid-dim [x dx d ^long nbins]
  (let [cx (count x)
        end (m/dec cx)
        buff (double-array cx)
        delta-d (stats/mean d)]
    (aset buff 0 ^double (x 0))
    (aset buff end ^double (x end))
    (loop [i (long 1)
           j (long 0)
           sd 0.0]
      (when (m/< i nbins)
        (let [[^long nj ^double nsd] (loop [nsd sd
                                            nj j]
                                       (if (m/< nsd delta-d)
                                         (recur (m/+ nsd ^double (d nj))
                                                (m/inc nj))
                                         [nj nsd]))
              nsd (m/- nsd delta-d)
              nj- (m/dec nj)]
          (aset buff i (m/- ^double (x nj)
                            (m// (m/* nsd ^double (dx nj-)) ^double (d nj-))))
          (recur (m/inc i) (long nj) nsd))))
    (vec buff)))

(defn- calc-dx [x]
  (mapv (fn [[^double a ^double b]]
          (m/- b a)) (partition 2 1 x)))

(defn- update-grid [xs dxs ds ^long nbins]
  (let [nxs (mapv (fn [x dx d]
                    (update-grid-dim x dx d nbins)) xs dxs ds)
        ndxs (mapv calc-dx nxs)]
    [nxs ndxs]))

(defn vegas
  "VEGAS - Monte Carlo integration"
  ([f] (vegas f nil))
  ([f opts] (vegas f [[0.0 1.0] [0.0 1.0]] opts))
  ([f bounds {:keys [^long nbins ^long max-iters ^long ncalls ^double rel ^double abs
                     ^double alpha ^long warmup]
              :or {nbins 1000 max-iters 10 ncalls 10000 rel 5.0e-4 abs 5.0e-4
                   alpha 1.5 warmup 0}
              :as options}]
   (let [dims (count bounds)
         max-iters (m/+ max-iters warmup)
         
         ;; convert y-space to x-space
         y->x (fn ^double [x dx ^long dim ^double y]
                (let [nby (m/* nbins y)
                      i (unchecked-long (m/floor nby))
                      dy (m/- nby i)]
                  (m/+ ^double (get-in x [dim i])
                       (m/* ^double (get-in dx [dim i]) dy))))
         
         ;; Jacobian function
         J (fn ^double [dx ^long dim ^double y]
             (m/* nbins ^double (get-in dx [dim (unchecked-long (m/floor (m/* nbins y)))])))]
     
     (loop [iter (long 0)
            [x dx] (vegas-build-initial-grid bounds nbins)
            random-sequence (vegas-random-sequence dims options)
            integrals []
            rev-sigma-squares []]
       (let [[samples rst] (split-at ncalls random-sequence)
             xss (->> samples
                      (map (fn [ys]
                             (vec (map-indexed (fn [^long dim ^double y]
                                                 (y->x x dx dim y)) ys)))))
             Js (map (fn [ys]
                       (reduce m/fast* (map-indexed (partial J dx) ys))) samples)
             imat (map (fn [ys]
                         (mapv (fn [^double y]
                                 (unchecked-long (m/floor (m/* nbins y)))) ys)) samples)
             fevals (map f xss)
             Jsf (map m/fast* Js fevals)]
         
         (if (m/< iter warmup)

           (let [d (calculate-d dims Jsf imat nbins alpha)
                 x-dx (update-grid x dx d nbins)]
             (recur (m/inc iter) x-dx rst integrals rev-sigma-squares))

           (let [integral-mc (m// (v/sum Jsf) ncalls)
                 variance-mc (m/+ (m// (- (m// (v/sum (v/sq Jsf)) ncalls)
                                          (m/sq integral-mc)) (m/dec ncalls))
                                  m/MACHINE-EPSILON)
                 new-integrals (conj integrals integral-mc)
                 new-rev-sigma-squares (conj rev-sigma-squares (m// variance-mc))
                 
                 [^double Itot ^double sd] (Itot-sd new-integrals new-rev-sigma-squares)]
             
             (if (or (m/== iter max-iters)
                     (and (m/< (m/abs (m// sd Itot)) rel)
                          (m/< (m/abs sd) abs)))
               (let [chisq (-> new-integrals
                               (v/shift (m/- Itot))
                               (v/sq)
                               (v/emult new-rev-sigma-squares)
                               (v/sum))]
                 {:iters (m/- iter warmup)
                  :result Itot
                  :sd sd
                  :chi2 (/ chisq (count new-integrals))})
               (let [d (calculate-d dims Jsf imat nbins alpha)
                     x-dx (update-grid x dx d nbins)]
                 (recur (m/inc iter) x-dx rst new-integrals new-rev-sigma-squares))))))))))

;;

;; integral_0^1 integral_0^1 sin(x) cos(y) dx dy = 2 sin^2(1/2) sin(1)≈0.38682227139505565895449238867442658

(def result-sinxcosy 0.38682227139505565895449238867442658)

(defn sinxcosy
  ^double [v]
  (m/* (m/sin (v 0)) (m/cos (v 1))))

(vegas sinxcosy {:max-iters 130 :ncalls 10000 :random-sequence :r2 :jitter 0.25
                 :warmup 1})

(defn x2+y2
  ^double [v]
  (v/sum (v/sq v)))

(vegas x2+y2 [[-1 1] [-1 1]] {:random-sequence :r2 :jitter 0.2})

(defn sinu
  ^double [v]
  (v/sum (v/sin v)))

(defn cosu
  ^double [v]
  (v/sum (v/cos v)))

(vegas sinu [[1 3] [1 3] [1 3]] {:warmup 2 :max-iters 100 :random-sequence :uniform :jitter 0.25
                                 :rel 1.0e-3 :abs 1.0e-3})
