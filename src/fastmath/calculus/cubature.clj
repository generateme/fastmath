(ns fastmath.calculus.cubature
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.quadrature :as quadrature]
            [clojure.math.combinatorics :as combo])
  (:import [fastmath.vector Vec2]
           [fastmath.calculus.common Estimate]
           [java.util PriorityQueue]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; hcubature
;; https://github.com/JuliaMath/HCubature.jl/blob/master/src/HCubature.jl
;; https://www.sciencedirect.com/science/article/pii/0771050X8090039X

(defn- combos
  [^long k ^double lambda ^long n]
  (mapv vec (combo/permutations (concat (repeat k lambda)
                                        (repeat (m/long-sub n k) 0.0)))))

(defn- signcombos
  [^long k ^double lambda ^long n]
  (->> (apply combo/cartesian-product (repeat k [lambda (m/- lambda)]))
       (map (partial concat (repeat (m/long-sub n k) 0.0)))
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
     :w' (v/mult [(m// (m/+ 729.0 (m/* -950.0 n) (m/* 50.0 n n)) 729.0)
                  0.5041152263374485 ;; 245/486
                  (m// (m/- 265.0 (m/* 100.0 n)) 1458.0)
                  0.03429355281207133] ;; 25/729
                 twon)
     :gm-evals (m/long-add 1 (m/long-mult 4 n) (m/long-mult 2 n (m/long-dec n)) twon)}))

(def ^:private genz-malik (memoize genz-malik-))

(deftype CubatureBox [a b ^double I ^double E ^long kdiv]
  Estimate
  (value [_] I)
  (error [_] E))

(defn- integrate-gm
  "Genz-Malik rule on the box `[lower, upper]`: integral, error estimate and the dimension to subdivide next."
  ^CubatureBox [f lower upper dims gp1 gp2 gp3 gp4 w w']
  (let [dims (long dims)
        c (v/mult (v/add lower upper) 0.5)
        delta (vec (v/mult (v/abs (v/sub upper lower)) 0.5))
        V (v/prod delta)
        ^double f1 (f c)
        twelvef1 (m/* 12.0 f1)

        f23s (mapv (fn [^long i]
                     (let [p2 (v/emult delta (gp1 i))
                           f2i (m/+ ^double (f (v/add c p2))
                                    ^double (f (v/sub c p2)))
                           p3 (v/emult delta (gp2 i))
                           f3i (m/+ ^double (f (v/add c p3))
                                    ^double (f (v/sub c p3)))]
                       (Vec2. f2i f3i))) (range dims))

        divdiff (mapv (fn [^Vec2 v]
                        (m/abs (m/+ (.y v) twelvef1 (m/* -7.0 (.x v))))) f23s)

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

        E (m/abs (m/- I I'))
        deltaf (m// E (m/* (m/pow 10.0 dims) V))

        ;; the dimension with the largest fourth difference (HCubature): `mx` is the largest difference found,
        ;; within `deltaf` of it the widest dimension wins
        kdivide (->> (range dims)
                     (reduce (fn [[^double mx ^long kdv] ^long dim]
                               (let [^double dd (divdiff dim)
                                     d (m/- dd mx)]
                                 (cond
                                   (m/> d deltaf) [dd dim]
                                   (and (m/<= (m/abs d) deltaf)
                                        (m/> ^double (delta dim)
                                             ^double (delta kdv))) [mx dim]
                                   :else [mx kdv]))) [0.0 0])
                     (second))]

    (CubatureBox. lower upper I E kdivide)))

(defn- enumerate-boxes
  [lower upper div]
  (->> (map (fn [^double a ^double b]
              (->> (m/slice-range a b (m/long-inc (long div)))
                   (partition 2 1))) lower upper)
       (apply combo/cartesian-product)
       (map (partial apply map vector))))

(defn- split-box
  [^CubatureBox box]
  (let [split-id (.kdiv box)
        a (.a box)
        b (.b box)
        ^double ida (a split-id)
        ^double idb (b split-id)
        mid (m/* 0.5 (m/+ ida idb))]
    [a (assoc b split-id mid)
     (assoc a split-id mid) b]))

(defn- check-initdiv
  ^long [n]
  (if (and (integer? n) (m/>= (long n) 1))
    (long n)
    (throw (ex-info "Option :initdiv must be an integer >= 1"
                    {:option :initdiv :value n :expected "integer >= 1"}))))

(defn- cubature-nd
  "Adaptive Genz-Malik cubature in `dims >= 2` dimensions on oriented bounds, returns the result (or the info map)."
  [f lower upper rel abs max-evals max-iters initdiv info?]
  (let [rel (double rel)
        abs (double abs)
        max-evals (long max-evals)
        max-iters (long max-iters)
        [f lower upper] (cc/subst-multi f lower upper)
        dims (count lower)
        {:keys [^long gm-evals p w w']} (genz-malik dims)
        [gp1 gp2 gp3 gp4] p
        split-evals (m/long-mult 2 gm-evals)
        ^PriorityQueue queue (cc/max-error-queue)]
    (doseq [[a b] (enumerate-boxes lower upper initdiv)]
      (.add queue (integrate-gm f a b dims gp1 gp2 gp3 gp4 w w')))
    (let [initial-evals (m/long-mult gm-evals (.size queue))
          [^double I0 ^double E0] (cc/estimates-sum queue)]
      (loop [iters (long 1)
             evals initial-evals
             I I0
             E E0]
        (let [fail? (cond
                      (m/<= E (m/max (m/* rel (m/abs I)) abs)) false
                      (m/>= iters max-iters) :max-iters
                      (m/> (m/+ evals split-evals) max-evals) :max-evals
                      :else :continue)]
          (if (identical? fail? :continue)
            (let [^CubatureBox worst (.poll queue)
                  [a1 b1 a2 b2] (split-box worst)
                  ^CubatureBox box1 (integrate-gm f a1 b1 dims gp1 gp2 gp3 gp4 w w')
                  ^CubatureBox box2 (integrate-gm f a2 b2 dims gp1 gp2 gp3 gp4 w w')]
              (.add queue box1)
              (.add queue box2)
              (recur (m/long-inc iters)
                     (m/+ evals split-evals)
                     (m/+ (.I box1) (.I box2) (m/- I (.I worst)))
                     (m/+ (.E box1) (.E box2) (m/- E (.E worst)))))
            (let [[^double I ^double E] (cc/estimates-sum queue)]
              (if info?
                {:result I :error E
                 :evaluations evals
                 :iterations iters
                 :subdivisions (.size queue)
                 :fail? fail?}
                I))))))))

(defn cubature
  "Adaptive Genz-Malik cubature of `f` on the box `[lower, upper]`, see [[fastmath.calculus/cubature]] for options and result."
  ([f lower upper] (cubature f lower upper nil))
  ([f lower upper {:keys [rel abs max-evals max-iters initdiv info?]
                   :or {rel 1.0e-7 abs 1.0e-7 max-evals Integer/MAX_VALUE max-iters 64 initdiv 2 info? false}
                   :as options}]
   (let [initdiv (check-initdiv initdiv)
         [sign lower upper] (cc/normalize-bounds-multi lower upper)]
     (cc/orient sign
                (if (m/== (count lower) 1)
                  (quadrature/gk-quadrature (fn ^double [^double x] (f [x])) (lower 0) (upper 0)
                                            (assoc options :rel rel :abs abs :max-evals max-evals :max-iters max-iters
                                                   :initdiv initdiv :info? info?))
                  (cubature-nd f lower upper rel abs max-evals max-iters initdiv info?))))))
