(ns fastmath.calculus.quadrature
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.common :as cc])
  (:import [java.util TreeSet Comparator]
           [fastmath.java Array]
           [fastmath.vector Vec2]
           [org.apache.commons.math3.linear Array2DRowRealMatrix EigenDecomposition]))



;; quadgk
;; https://github.com/JuliaMath/QuadGK.jl/blob/master/src/gausskronrod.jl
;; https://github.com/ESSS/cquadpack/blob/master/src/

(defn eigpoly
  ^double [b ^double z]
  (let [[^double p
         ^double pderiv] (reduce (fn [[^double d1 ^double d1deriv ^double d2 ^double d2deriv] ^double Hev]
                                   (let [b2 (m/sq Hev)
                                         a z
                                         d (m/- (m/* a d1)
                                                (m/* b2 d2))
                                         dderiv (m/+ d1 (m/- (m/* a d1deriv)
                                                             (m/* b2 d2deriv)))]
                                     [d dderiv d1 d1deriv])) [z 1.0 1.0 0.0] b)]
    (m// p pderiv)))

(defn newton-iterations
  [nb ^double lambda]
  (reduce (fn [^double lambda _]
            (let [dlambda (eigpoly nb lambda)]
              (if (m/valid-double? dlambda)
                (let [nlambda (m/- lambda dlambda)]
                  (if (m/< (m/abs dlambda)
                           (m/ulp nlambda))
                    (let [dlambda (eigpoly nb nlambda)]
                      (reduced (if (m/valid-double? dlambda)
                                 (m/- nlambda dlambda)
                                 nlambda)))
                    nlambda))
                (reduced lambda)))) lambda (range 1000)))

(defn get-quadrature-points
  [b ^long n]
  (let [cb (count b)
        m (m/inc cb)
        buff (double-array (m/sq m))]
    (doseq [^long i (range cb)
            :let [^double v (b i)]]
      (Array/set2d buff m i (m/inc i) v)
      (Array/set2d buff m (m/inc i) i v))
    (let [lambdas (->> buff
                       (partition m)
                       (m/seq->double-double-array)
                       (Array2DRowRealMatrix.)
                       (EigenDecomposition.)
                       (.getRealEigenvalues)
                       (reverse)
                       (take n))]
      (mapv (partial newton-iterations b) lambdas))))

(defn eigvec
  ^double [b ^double lambda]
  (let [v1 (m// lambda ^double (b 0))]
    (->> (partition 2 1 b)
         (reduce (fn [[curr ^double v-2 ^double v-1] [^double b-2 ^double b-1]]
                   (let [newv (m// (m/- (m/* lambda v-1)
                                        (m/* b-2 v-2)) b-1)]
                     [(conj curr newv) v-1 newv])) [[1.0 v1] 1.0 v1])
         (first)
         (v/mag)
         (m//))))

;; https://github.com/JuliaMath/QuadGK.jl/blob/master/src/gausskronrod.jl#L448-L497
(defn kronrodjacobi
  [b ^long n]
  (let [s (vec (repeat (m/+ (m/quot n 2) 2) 0.0))
        t (assoc (vec (repeat (count s) 0.0)) 1 (b n))
        [s t] (reduce (fn [[s t] ^long m]
                        [t (->> (range (m/quot (m/inc m) 2) -1 -1)
                                (reduce (fn [[s ^double u] ^long k]
                                          (let [nu (m/+ u (m/- (m/* ^double (b (m/+ k n))
                                                                    ^double (s k))
                                                               (if (m/> m k)
                                                                 (m/* ^double (b (m/- m k 1))
                                                                      ^double (s (m/inc k)))
                                                                 0.0)))]
                                            [(assoc s (m/inc k) nu) nu]))
                                        [s 0.0])
                                (first))])
                      [s t] (range (m/dec n)))
        s (reduce (fn [s ^long j]
                    (assoc s (m/inc j) (s j))) s (range (m/quot n 2) -1 -1))
        nb (->> (range (m/dec n) (m/+ n n -2))
                (reduce (fn [[b s t] ^long m]
                          (let [news (->> (range (m/- m n -1) (m/inc (m/quot (m/dec m) 2)))
                                          (reduce (fn [[s ^double u] ^long k]
                                                    (let [j (m/dec (m/- n (m/- m k)))
                                                          nu (m/- u (m/- (m/* ^double (b (m/+ k n))
                                                                              ^double (s (m/inc j)))
                                                                         (m/* ^double (b (m/- m k 1))
                                                                              ^double (s (m/+ j 2)))))]
                                                      [(assoc s (m/inc j) nu) nu]))
                                                  [s 0.0])
                                          (first))
                                k (m/quot (m/inc m) 2)
                                j (m/- n (m/+ (m/- m k) 2))
                                nb (if (m/not== (m/* 2 k) m)
                                     (assoc b (m/+ k n) (m// ^double (news (m/+ j 1))
                                                             ^double (news (m/+ j 2))))
                                     b)]
                            [nb t news]))
                        [b s t])
                (first)
                (v/sqrt))
        x (get-quadrature-points nb (m/inc n))
        w (mapv (fn [^double x] (m/* 2.0 (m/sq (eigvec nb x)))) x)
        nb (map-indexed (fn [^long j ^double b]
                          (let [j (m/inc j)]
                            (if (m/< j n)
                              (m// j (m/sqrt (m/dec (m/* 4.0 j j))))
                              b))) nb)
        gw (mapv (fn [^double x]
                   (m/* 2.0 (m/sq (eigvec (vec (take (m/dec n) nb)) x))))
                 (->> x rest (take-nth 2)))]
    {:x x :w w :gw gw :gk-evals (m/+ (m/* 4 n) 2)}))

(defn kronrod-
  [^long n]
  (let [last-j (m/quot (m/inc (m/* 3 n)) 2)
        b (mapv (fn [^double j]
                  (if (m/<= j last-j)
                    (let [j2 (m/sq j)]
                      (m// j2 (m/dec (m/* 4.0 j2)) ))
                    0.0)) (range 1 (inc (m/* 2 n))))]
    (kronrodjacobi b n)))

(def kronrod (memoize kronrod-))

(deftype QuadGKSegment [^Vec2 ab ^double I ^double E])

(defn integrate-gk
  ^QuadGKSegment [f ^Vec2 ab x w gw]
  (let [s (m/* 0.5 (m/- (.y ab) (.x ab)))
        lastxw (m/dec (count w))
        lastxw- (m/dec lastxw)
        lastgw (m/dec (count gw))
        n1 (m/- 1 (m/bit-and (m/inc lastxw) 1))
        a (.x ab)
        ^double f0 (f (m/+ a s))
        Igk (if (m/zero? n1)
              (Vec2. 0.0 (m/* f0 ^double (w lastxw)))
              (let [^double x-1 (x lastxw-)]
                (Vec2. (m/* f0 ^double (gw lastgw))
                       (m/+ (m/* f0 ^double (w lastxw))
                            (m/* (m/+ ^double (f (m/+ a (m/* (m/+ 1.0 x-1) s)))
                                      ^double (f (m/+ a (m/* (m/- 1.0 x-1) s))))
                                 ^double (w lastxw-))))))
        ^Vec2 Igk (-> (->> (if (m/zero? n1) gw (butlast gw))
                           (map-indexed (fn [^long i ^double gwi]
                                          (let [i (m/inc i)
                                                ii (m/dec (m/* 2 i))
                                                ii-  (m/dec ii)
                                                ^double x2i (x ii)
                                                ^double x2i-1 (x ii-)
                                                ^double w2i (w ii)
                                                ^double w2i-1 (w ii-)
                                                fg (m/+ ^double (f (m/+ a (m/* (m/+ 1.0 x2i) s)))
                                                        ^double (f (m/+ a (m/* (m/- 1.0 x2i) s))))
                                                fk (m/+ ^double (f (m/+ a (m/* (m/+ 1.0 x2i-1) s)))
                                                        ^double (f (m/+ a (m/* (m/- 1.0 x2i-1) s))))]
                                            (Vec2. (m/* fg gwi)
                                                   (m/+ (m/* fg w2i)
                                                        (m/* fk w2i-1))))))
                           (reduce v/add Igk))
                      (v/mult s))
        E (m/abs (m/- (.y Igk) (.x Igk)))]
    (QuadGKSegment. ab (.y Igk) E)))

(def gk-comparator (comparator (fn [^QuadGKSegment x ^QuadGKSegment y]
                               (compare (.E x) (.E y)))))

(defn build-initial-segments
  [f lower upper initdiv x w gw]
  (let [^TreeSet boxes (TreeSet. ^Comparator gk-comparator)
        segments (->> (m/slice-range lower upper (m/inc ^long initdiv))
                      (partition 2 1))]
    (doseq [[a b] segments]
      (.add boxes (integrate-gk f (Vec2. a b) x w gw)))
    boxes))

(defn sum-segments
  [segments]
  (reduce (fn [[^double I ^double E] ^QuadGKSegment s]
            [(m/+ I (.I s))
             (m/+ E (.E s))]) [0.0 0.0] segments))

(defn- split-segment
  [^QuadGKSegment s]
  (let [^Vec2 ab (.ab s)
        mid (m/* 0.5 (v/sum ab))]
    [(Vec2. (.x ab) mid)
     (Vec2. mid (.y ab))]))

(defn gk-quadrature
  ([f lower upper] (gk-quadrature f lower upper nil))
  ([f lower upper {:keys [^double rel ^double abs ^int max-evals ^int max-iters
                          ^int integration-points ^int initdiv info?]
                   :or {rel 1.0e-8 abs 1.0e-8 max-evals Integer/MAX_VALUE max-iters 64
                        integration-points 7 initdiv 1 info? false}}]
   (let [[f ^double lower ^double upper] (cc/subst-1d f lower upper)
         {:keys [^long gk-evals x w gw]} (kronrod integration-points)
         ^TreeSet segments (build-initial-segments f lower upper initdiv x w gw)
         initial-evals (m/* gk-evals (count segments))
         [^double I ^double E] (sum-segments segments)]
     (if (or (m/<= E (m/max (m/* rel (m/abs I)) abs))
             (m/>= initial-evals max-evals))
       (if info?
         {:result I :error E
          :evaluations initial-evals
          :iterations 1
          :subdivisions (count segments)}
         I)
       (loop [iters (long 1)
              evals initial-evals
              I I
              E E]
         (let [niters (m/inc iters)
               nevals (m/+ evals gk-evals)
               ^QuadGKSegment max-E-segment (.pollLast segments)
               [ab1 ab2] (split-segment max-E-segment)
               ^QuadGKSegment nsegment1 (integrate-gk f ab1 x w gw)
               ^QuadGKSegment nsegment2 (integrate-gk f ab2 x w gw)
               nI (m/+ (.I nsegment1) (.I nsegment2) (m/- I (.I max-E-segment)))
               nE (m/+ (.E nsegment1) (.E nsegment2) (m/- E (.E max-E-segment)))
               fail? (cond
                       (m/>= nevals max-evals) :max-evals
                       (m/>= niters max-iters) :max-iters
                       :else false)]
           (.add segments nsegment1)
           (.add segments nsegment2)
           (if (or (m/<= nE (m/max (m/* rel (m/abs nI)) abs)) fail?)
             (let [[^double I ^double E] (sum-segments segments)]
               (if info?
                 {:result I :error E
                  :evaluations nevals
                  :iterations niters
                  :subdivisions (count segments)
                  :fail? fail?}
                 I))
             (recur niters nevals nI nE))))))))
