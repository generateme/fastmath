(ns fastmath.calculus.cubature
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]

            [fastmath.calculus.common :as cc]
            [clojure.math.combinatorics :as combo])
  (:import [fastmath.vector Vec2]
           [java.util TreeSet Comparator]))

;; hcubature
;; https://github.com/JuliaMath/HCubature.jl/blob/master/src/HCubature.jl
;; https://www.sciencedirect.com/science/article/pii/0771050X8090039X

(defn- combos
  [^long k ^double lambda ^long n]
  (mapv vec (combo/permutations (concat (repeat k lambda)
                                        (repeat (m/- n k) 0.0)))))

(defn- signcombos
  [^long k ^double lambda ^long n]
  (->> (apply combo/cartesian-product (repeat k [lambda (m/- lambda)]))
       (map (partial concat (repeat (m/- n k) 0.0)))
       (mapcat combo/permutations)
       (distinct)
       (mapv vec)))

(defn- genz-malik-
  [^long n]
  (let [twon (m/<< 1 n)]
    {:p [(combos 1 0.35856858280031806 n)     ;; sqrt(9/70)
         (combos 1 0.9486832980505138 n)      ;; sqrt(9/10)
         (signcombos 2 0.9486832980505138 n)  ;; sqrt(9/10)
         (signcombos n 0.6882472016116853 n)] ;; sqrt(9/19)
     :w [(m/* twon (m// (m/+ 12824.0
                             (m/* -9120.0 n)
                             (m/* 400.0 n n)) 19683.0))
         (m/* twon 0.14936747447035512) ;; 980/6561
         (m/* twon (m// (m/- 1820.0 (m/* 400.0 n)) 19683.0))
         (m/* twon 0.010161052685058172) ;; 200/19683
         0.34847330183407] ;; 6859/19683
     :w' (v/mult [(m// (m/+ 729.0 (m/* -950 n) (m/* 50.0 n n)) 729.0)
                  0.5041152263374485 ;; 245/486
                  (m// (m/- 265.0 (m/* 100.0 n)) 1458.0)
                  0.03429355281207133] ;; 25/729
                 twon)
     :gm-evals (m/inc (m/+ (m/* 4 n)
                           (m/* 2 n (m/dec n))
                           (m/<< 1 n)))}))

(def ^:private genz-malik (memoize genz-malik-))

(deftype CubatureBox [a b ^double I ^double E ^long kdiv])

(defn- integrate-gm
  ^CubatureBox [f lower upper dims gp1 gp2 gp3 gp4 w w']
  (let [c (v/mult (v/add lower upper) 0.5)
        delta (vec (v/mult (v/abs (v/sub upper lower)) 0.5))
        V (v/prod delta)
        ^double f1 (f c)
        twelfef1 (m/* 12.0 f1)
        
        f23s (map (fn [^long i]
                    (let [p2 (v/emult delta (gp1 i))
                          f2i (m/+ ^double (f (v/add c p2))
                                   ^double (f (v/sub c p2)))
                          p3 (v/emult delta (gp2 i))
                          f3i (m/+ ^double (f (v/add c p3))
                                   ^double (f (v/sub c p3)))]
                      (Vec2. f2i f3i))) (range dims))
        
        divdiff (mapv (fn [^Vec2 v]
                        (m/abs (m/+ (.y v) twelfef1 (m/* -7.0 (.x v))))) f23s)
        
        ^Vec2 f23 (reduce v/add f23s)
        
        f4 (->> gp3
                (map #(f (v/add c (v/emult delta %))))
                (v/sum))
        f5 (->> gp4
                (map #(f (v/add c (v/emult delta %))))
                (v/sum))

        I (m/* V (m/+ (m/* ^double (w 0) f1)
                      (m/* ^double (w 1) (.x f23))
                      (m/* ^double (w 2) (.y f23))
                      (m/* ^double (w 3) f4)
                      (m/* ^double (w 4) f5)))

        I' (m/* V (m/+ (m/* ^double (w' 0) f1)
                       (m/* ^double (w' 1) (.x f23))
                       (m/* ^double (w' 2) (.y f23))
                       (m/* ^double (w' 3) f4)))
        
        E (m/abs (- I I'))
        deltaf (m// E (m/* (m/pow 10.0 dims) V))
        
        kdivide (->> (range dims)
                     (reduce (fn [[^double mx ^long kdv] ^long dim]
                               (let [^double dd (divdiff dim) 
                                     d (m/- dd mx)
                                     pred (m/> d deltaf)
                                     nkdv (if pred dim kdv)
                                     nmx (if pred dd dim)
                                     nkdv (if (and (m/<= (m/abs d) deltaf)
                                                   (m/> ^double (delta dim)
                                                        ^double (delta nkdv))) dim nkdv)]
                                 [nmx nkdv])) [0.0 0])
                     (second))]
    
    (CubatureBox. lower upper I E kdivide)))

(def ^:private cb-comparator (comparator (fn [^CubatureBox x ^CubatureBox y]
                                           (compare (.E x) (.E y)))))

(defn- enumerate-boxes
  [lower upper ^long div]
  (->> (map (fn [^double a ^double b]
              (->> (m/slice-range a b (m/inc div))
                   (partition 2 1))) lower upper)
       (apply combo/cartesian-product)
       (map (partial apply map vector))))

(defn- build-initial-boxes
  [f lower upper div dims gp1 gp2 gp3 gp4 w w']
  (let [^TreeSet boxes (TreeSet. ^Comparator cb-comparator)]
    (doseq [[a b] (enumerate-boxes lower upper div)]
      (.add boxes (integrate-gm f a b dims gp1 gp2 gp3 gp4 w w')))
    boxes))

(defn- sum-boxes
  [boxes]
  (reduce (fn [[^double I ^double E] ^CubatureBox b]
            [(m/+ I (.I b))
             (m/+ E (.E b))]) [0.0 0.0] boxes))

(defn- split-box
  [^CubatureBox box]
  (let [split-id (.kdiv box)
        a (.a box)
        b (.b box)
        ^double ida (a split-id)
        ^double idb (b split-id)
        mid (* 0.5 (m/+ ida idb))]
    [a (assoc b split-id mid)
     (assoc a split-id mid) b]))

(defn cubature
  ([f lower upper] (cubature f lower upper nil))
  ([f lower upper {:keys [^double rel ^double abs ^int max-evals ^int max-iters ^int initdiv info?]
                   :or {rel 1.0e-7 abs 1.0e-7 max-evals Integer/MAX_VALUE max-iters 64 initdiv 2 info? false}}]
   (let [[f lower upper] (cc/subst-multi f lower upper)
         dims (count lower)
         {:keys [^long gm-evals p w w']} (genz-malik dims)
         [gp1 gp2 gp3 gp4] p
         ^TreeSet boxes (build-initial-boxes f lower upper initdiv dims gp1 gp2 gp3 gp4 w w')
         initial-evals (m/* gm-evals (count boxes))
         [^double I ^double E] (sum-boxes boxes)]
     (if (or (m/<= E (m/max (m/* rel (m/abs I)) abs))
             (m/>= initial-evals max-evals))
       (if info?
         {:result I :error E
          :evaluations initial-evals
          :iterations 1
          :subdivisions (count boxes)}
         I)
       (loop [iters (long 1)
              evals initial-evals
              I I
              E E]
         (let [niters (m/inc iters)
               nevals (m/+ evals gm-evals)
               ^CubatureBox max-E-box (.pollLast boxes)
               [a1 b1 a2 b2] (split-box max-E-box)
               ^CubatureBox nbox1 (integrate-gm f a1 b1 dims gp1 gp2 gp3 gp4 w w')
               ^CubatureBox nbox2 (integrate-gm f a2 b2 dims gp1 gp2 gp3 gp4 w w')
               nI (m/+ (.I nbox1) (.I nbox2) (m/- I (.I max-E-box)))
               nE (m/+ (.E nbox1) (.E nbox2) (m/- E (.E max-E-box)))
               fail? (cond
                       (m/>= nevals max-evals) :max-evals
                       (m/>= niters max-iters) :max-iters
                       :else false)]
           (.add boxes nbox1)
           (.add boxes nbox2)
           (if (or (m/<= nE (m/max (m/* rel (m/abs nI)) abs)) fail?)
             (let [[^double I ^double E] (sum-boxes boxes)]
               (if info?
                 {:result I :error E
                  :evaluations nevals
                  :iterations niters
                  :subdivisions (count boxes)
                  :fail? fail?}
                 I))
             (recur niters nevals nI nE))))))))

