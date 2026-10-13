(ns fastmath.calculus.quadrature
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.common :as cc]
            [fastmath.matrix :as mat])
  (:import [java.util PriorityQueue]
           [fastmath.java Array]
           [fastmath.vector Vec2]
           [fastmath.calculus.common Estimate]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; quadgk
;; https://github.com/JuliaMath/QuadGK.jl/blob/master/src/gausskronrod.jl
;; https://github.com/ESSS/cquadpack/blob/master/src/

(defn- eigpoly
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

(defn- newton-iterations
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

(defn- get-quadrature-points
  [b ^long n]
  (let [cb (count b)
        m (m/long-inc cb)
        buff (double-array (m/sq m))]
    (doseq [^long i (range cb)
            :let [^double v (b i)]]
      (Array/set2d buff m i (m/inc i) v)
      (Array/set2d buff m (m/inc i) i v))
    (let [lambdas (->> buff
                       (partition m)
                       (m/seq->double-double-array)
                       (mat/eigenvalues)
                       (map first)
                       (take n))]
      (mapv (partial newton-iterations b) lambdas))))

(defn- eigvec
  ^double [b ^double lambda]
  (if (empty? b)
    1.0 ;; first component of the normalized eigenvector of a 1x1 matrix
    (let [v1 (m// lambda ^double (b 0))]
      (->> (partition 2 1 b)
           (reduce (fn [[curr ^double v-2 ^double v-1] [^double b-2 ^double b-1]]
                     (let [newv (m// (m/- (m/* lambda v-1)
                                          (m/* b-2 v-2)) b-1)]
                       [(conj curr newv) v-1 newv])) [[1.0 v1] 1.0 v1])
           (first)
           (v/mag)
           (m//)))))

;; https://github.com/JuliaMath/QuadGK.jl/blob/master/src/gausskronrod.jl#L448-L497
(defn- kronrodjacobi
  [b ^long n]
  (let [s (vec (repeat (m/long-add (m/long-quot n 2) 2) 0.0))
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
        x (get-quadrature-points nb (m/long-inc n))
        w (mapv (fn [^double x] (m/* 2.0 (m/sq (eigvec nb x)))) x)
        nb (map-indexed (fn [^long j ^double b]
                          (let [j (m/inc j)]
                            (if (m/< j n)
                              (m// j (m/sqrt (m/dec (m/* 4.0 j j))))
                              b))) nb)
        gw (mapv (fn [^double x]
                   (m/* 2.0 (m/sq (eigvec (vec (take (m/long-dec n) nb)) x))))
                 (->> x rest (take-nth 2)))]
    {:x x :w w :gw gw :gk-evals (m/long-add (m/long-mult 4 n) 2)}))

;; precomputed from julia
(def ^:private kronrod-7-result
  {:x
   [-9.9145537112081263920685469752598e-01,
    -9.4910791234275852452618968404809e-01,
    -8.6486442335976907278971278864098e-01,
    -7.415311855993944398638647732811e-01,
    -5.8608723546769113029414483825842e-01,
    -4.0584515137739716690660641207707e-01,
    -2.0778495500789846760068940377309e-01,
    0.0],
   :w
   [2.2935322010529224963732008059913e-02,
    6.3092092629978553290700663189093e-02,
    1.0479001032225018383987632254189e-01,
    1.4065325971552591874518959051021e-01,
    1.6900472663926790282658342659795e-01,
    1.9035057806478540991325640242055e-01,
    2.0443294007529889241416199923466e-01,
    2.0948214108472782801299917489173e-01],
   :gw
   [1.2948496616886969327061143267787e-01,
    2.797053914892766679014677714229e-01,
    3.8183005050511894495036977548818e-01,
    4.1795918367346938775510204081658e-01],
   :gk-evals 30})

(defn- kronrod-
  [^long n]
  (if (m/== n 7)
    kronrod-7-result
    (let [last-j (m/quot (m/inc (m/* 3 n)) 2)
          b (mapv (fn [^double j]
                    (if (m/<= j last-j)
                      (let [j2 (m/sq j)]
                        (m// j2 (m/dec (m/* 4.0 j2)) ))
                      0.0)) (range 1 (inc (m/* 2 n))))]
      (kronrodjacobi b n))))

(def ^{:arglists '([n])
       :doc "Gauss-Kronrod rule of order `n` (`n` Gauss points inside `2n + 1` Kronrod points) on the interval `[-1, 1]`.

  Parameters:

  - `n` (integer >= 1): the order of the rule.

  Returns a map with the keys:

  - `:x` - the `n + 1` nodes of the left half of the interval, in increasing order up to `0.0`; the other half is symmetric.
  - `:w` - the weights of the Kronrod rule for the nodes `:x`.
  - `:gw` - the weights of the embedded Gauss rule for the nodes of `:x` with even indexes counted from 1 (the Gauss nodes lie between the Kronrod ones).
  - `:gk-evals` - the number of evaluations of the integrand for two intervals, `4n + 2`.

  The rules are computed with the algorithm of QuadGK.jl and match it to 1e-15. The results are memoized; the rule of order 7 is precomputed."}
  kronrod (memoize kronrod-))

(deftype QuadGKSegment [^Vec2 ab ^double I ^double E]
  Estimate
  (value [_] I)
  (error [_] E))

(defn- integrate-gk
  ^QuadGKSegment [f ^Vec2 ab x w gw]
  (let [s (m/* 0.5 (m/- (.y ab) (.x ab)))
        lastxw (m/long-dec (count w))
        lastxw- (m/long-dec lastxw)
        lastgw (m/long-dec (count gw))
        n1 (m/long-sub 1 (bit-and (m/long-inc lastxw) 1))
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

(defn check-integration-points
  "Returns the order of the Gauss-Kronrod rule as a long; throws `ex-info` unless it is an integer >= 1."
  ^long [n]
  (if (and (integer? n) (m/>= (long n) 1))
    (long n)
    (throw (ex-info "Option :integration-points must be an integer >= 1"
                    {:option :integration-points :value n :expected "integer >= 1"}))))

(defn- check-initdiv
  ^long [n]
  (if (and (integer? n) (m/>= (long n) 1))
    (long n)
    (throw (ex-info "Option :initdiv must be an integer >= 1"
                    {:option :initdiv :value n :expected "integer >= 1"}))))

(defn- breakpoints
  "Sorted, distinct boundaries of the initial pieces in the variable of the substitution:
  the substituted bounds with the user points which lie strictly inside the interval."
  [points lower upper tl tu]
  (when-not (every? number? points)
    (throw (ex-info "Option :points must be a sequence of numbers"
                    {:option :points :value points :expected "sequence of numbers"})))
  (let [lower (double lower)
        upper (double upper)
        tl (double tl)
        tu (double tu)
        inv (cc/subst-1d-inverse lower upper)
        ts (->> points
                (filter (fn [^double p] (m/< lower p upper)))
                (map (fn [^double p] (inv p)))
                (filter (fn [^double t] (m/< tl t tu))))]
    (vec (distinct (concat [tl] (sort ts) [tu])))))

(defn- initial-segments
  [f breakpoints initdiv x w gw]
  (->> (partition 2 1 breakpoints)
       (mapcat (fn [[^double a ^double b]]
                 (->> (m/slice-range a b (m/long-inc (long initdiv)))
                      (partition 2 1))))
       (map (fn [[a b]] (integrate-gk f (Vec2. a b) x w gw)))))

(defn- split-segment
  [^QuadGKSegment s]
  (let [^Vec2 ab (.ab s)
        mid (m/* 0.5 (v/sum ab))]
    [(Vec2. (.x ab) mid)
     (Vec2. mid (.y ab))]))

(defn gk-quadrature
  "Adaptive Gauss-Kronrod quadrature of `f` on `[lower, upper]`, see [[fastmath.calculus/integrate]] for options and result."
  ([f lower upper] (gk-quadrature f lower upper nil))
  ([f lower upper {:keys [^double rel ^double abs ^long max-evals ^long max-iters
                          integration-points initdiv info? points]
                   :or {rel 1.0e-8 abs 1.0e-8 max-evals Integer/MAX_VALUE max-iters 64
                        integration-points 7 initdiv 1 info? false}}]
   (let [n (check-integration-points integration-points)
         initdiv (check-initdiv initdiv)
         [^double sign ^double lower ^double upper] (cc/normalize-bounds lower upper)
         [f ^double tl ^double tu] (cc/subst-1d f lower upper)
         {:keys [x w gw]} (kronrod n)
         segment-evals (m/long-inc (m/long-mult 2 n))
         split-evals (m/long-mult 2 segment-evals)
         ^PriorityQueue queue (cc/max-error-queue)]
     (doseq [segment (initial-segments f (breakpoints points lower upper tl tu) initdiv x w gw)]
       (.add queue segment))
     (let [initial-evals (m/long-mult segment-evals (.size queue))
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
             (let [^QuadGKSegment worst (.poll queue)
                   [ab1 ab2] (split-segment worst)
                   ^QuadGKSegment s1 (integrate-gk f ab1 x w gw)
                   ^QuadGKSegment s2 (integrate-gk f ab2 x w gw)]
               (.add queue s1)
               (.add queue s2)
               (recur (m/long-inc iters)
                      (m/+ evals split-evals)
                      (m/+ (.I s1) (.I s2) (m/- I (.I worst)))
                      (m/+ (.E s1) (.E s2) (m/- E (.E worst)))))
             (let [[^double I ^double E] (cc/estimates-sum queue)
                   res (if info?
                         {:result I :error E
                          :evaluations evals
                          :iterations iters
                          :subdivisions (.size queue)
                          :fail? fail?}
                         I)]
               (cc/orient sign res)))))))))
