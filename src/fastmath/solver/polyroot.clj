(ns fastmath.solver.polyroot
  "All the complex roots of a polynomial with real or complex coefficients.

  The roots are found simultaneously by a modified Laguerre iteration (Aberth's correction for the roots already known is applied to the Laguerre step), started from points on circles whose radii come from the Newton polygon of the coefficients. The iteration of a root ends when the backward error of the root, the relative size of the polynomial value, is within the rounding error. The method follows T. R. Cameron, An effective implementation of a modified Laguerre method for the roots of a polynomial, Numerical Algorithms 82 (2019). The polynomial is evaluated in its reversed form for the roots outside of the unit circle, so large coefficients and roots do not overflow."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c]
            [fastmath.complex :as cplx]
            [fastmath.protocols.polynomials :as prot]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'double} SIGMA 0.7)
(def ^{:private true :const true :tag 'double} SNAP-EPSILON (* 8.0 m/MACHINE-EPSILON2))

;; input

(defn- invalid-coefficient
  [index reason]
  (throw (ex-info (str "Invalid coefficient " index ": " reason)
                  {:option :coefficients :index index :reason reason})))

(defn- coefficient
  "Real and imaginary part of the coefficient `x` at `index`: a number, a `Vec2` or a pair `[re im]`."
  [x index]
  (let [[re im] (cond
                  (number? x) [(double x) 0.0]
                  (instance? Vec2 x) (vec x)
                  (and (sequential? x) (== 2 (count x)) (every? number? x)) [(double (first x)) (double (second x))]
                  :else (invalid-coefficient index "must be a number, a Vec2 or a pair [re im]"))]
    (when-not (and (m/valid-double? re) (m/valid-double? im))
      (invalid-coefficient index "must be finite"))
    [re im]))

(defn- parse-coefficients
  "Arrays `[re im]` of the coefficients in ascending order of power."
  [coeffs]
  (let [xs (cond
             (satisfies? prot/PolynomialProto coeffs) (map double (prot/coeffs coeffs))
             (and (some? coeffs) (or (sequential? coeffs) (.isArray (class coeffs)))) (seq coeffs)
             :else (c/throw-invalid-option :coefficients (class coeffs) "must be a sequence of coefficients or a polynomial"))]
    (when (empty? xs) (c/throw-invalid-option :coefficients 0 "must not be empty"))
    (let [cs (vec (map-indexed (fn [i x] (coefficient x i)) xs))]
      [(double-array (map first cs)) (double-array (map second cs))])))

;; the polynomial and a root

(deftype Poly [^doubles re ^doubles im ^long n])

(deftype Root [^double re ^double im ^double berr ^double cnd ok])

;; evaluation

(defn- evaluate
  "Evaluates the polynomial `p` at `z = zr + i zi`.

  Returns a `double-array` `[ar ai br bi cr ci wsum usum reversed?]`. For `|z| <= 1`: `a = p(z)`, `b = p'(z)`, `c = p''(z) / 2`, `usum = sum |c_k| |z|^k` and `wsum` the same sum with the weights `3.8 k + 1` of the rounding error bound. For `|z| > 1` the same values of the reversed polynomial `q(u) = u^n p(1/u)` at `u = 1/z` and the sums divided by `|z|^n`; the last element is `1.0`."
  ^doubles [^Poly p ^double zr ^double zi]
  (let [^doubles cr (.re p)
        ^doubles ci (.im p)
        n (.n p)
        r (m/hypot zr zi)
        out (double-array 9)]
    (if (m/<= r 1.0)
      (loop [k (m/dec n)
             ar (Array/aget cr n) ai (Array/aget ci n) br 0.0 bi 0.0 qr 0.0 qi 0.0
             ws (m/* (m/+ (m/* 3.8 n) 1.0) (m/hypot (Array/aget cr n) (Array/aget ci n)))
             us (m/hypot (Array/aget cr n) (Array/aget ci n))]
        (if (m/< k 0)
          (do (Array/aset out 0 ar) (Array/aset out 1 ai) (Array/aset out 2 br) (Array/aset out 3 bi) (Array/aset out 4 qr) (Array/aset out 5 qi)
              (Array/aset out 6 ws) (Array/aset out 7 us) (Array/aset out 8 0.0))
          (let [nqr (m/+ (m/- (m/* zr qr) (m/* zi qi)) br)
                nqi (m/+ (m/+ (m/* zr qi) (m/* zi qr)) bi)
                nbr (m/+ (m/- (m/* zr br) (m/* zi bi)) ar)
                nbi (m/+ (m/+ (m/* zr bi) (m/* zi br)) ai)
                nar (m/+ (m/- (m/* zr ar) (m/* zi ai)) (Array/aget cr k))
                nai (m/+ (m/+ (m/* zr ai) (m/* zi ar)) (Array/aget ci k))
                ck (m/hypot (Array/aget cr k) (Array/aget ci k))]
            (recur (m/dec k) nar nai nbr nbi nqr nqi
                   (m/+ (m/* r ws) (m/* (m/+ (m/* 3.8 k) 1.0) ck))
                   (m/+ (m/* r us) ck)))))
      (let [[ur ui] (cplx/reciprocal (Vec2. zr zi))
            ur (double ur)
            ui (double ui)
            rr (m// 1.0 r)]
        (loop [k (long 1)
               ar (Array/aget cr 0) ai (Array/aget ci 0) br 0.0 bi 0.0 qr 0.0 qi 0.0
               ws (m/hypot (Array/aget cr 0) (Array/aget ci 0))
               us (m/hypot (Array/aget cr 0) (Array/aget ci 0))]
          (if (m/> k n)
            (do (Array/aset out 0 ar) (Array/aset out 1 ai) (Array/aset out 2 br) (Array/aset out 3 bi) (Array/aset out 4 qr) (Array/aset out 5 qi)
                (Array/aset out 6 ws) (Array/aset out 7 us) (Array/aset out 8 1.0))
            (let [nqr (m/+ (m/- (m/* ur qr) (m/* ui qi)) br)
                  nqi (m/+ (m/+ (m/* ur qi) (m/* ui qr)) bi)
                  nbr (m/+ (m/- (m/* ur br) (m/* ui bi)) ar)
                  nbi (m/+ (m/+ (m/* ur bi) (m/* ui br)) ai)
                  nar (m/+ (m/- (m/* ur ar) (m/* ui ai)) (Array/aget cr k))
                  nai (m/+ (m/+ (m/* ur ai) (m/* ui ar)) (Array/aget ci k))
                  ck (m/hypot (Array/aget cr k) (Array/aget ci k))]
              (recur (m/inc k) nar nai nbr nbi nqr nqi
                     (m/+ (m/* rr ws) (m/* (m/+ (m/* 3.8 k) 1.0) ck))
                     (m/+ (m/* rr us) ck)))))))
    out))

(defn- cond-number
  "Relative condition number of the root `z`: `sum |c_k| |z|^k / (|z| |p'(z)|)`."
  ^double [^doubles e ^long n ^double zr ^double zi]
  (let [us (Array/aget e 7)]
    (if (m/== (Array/aget e 8) 1.0)
      ;; z p'(z) = n p(z) - u q'(u) q-terms: the denominator is |n q(u) - u q'(u)|
      (let [[ur ui] (cplx/reciprocal (Vec2. zr zi))
            ur (double ur)
            ui (double ui)
            ;; n q - u q'
            tr (m/- (m/* n (Array/aget e 0)) (m/- (m/* ur (Array/aget e 2)) (m/* ui (Array/aget e 3))))
            ti (m/- (m/* n (Array/aget e 1)) (m/+ (m/* ur (Array/aget e 3)) (m/* ui (Array/aget e 2))))]
        (m// us (m/hypot tr ti)))
      (m// us (m/* (m/hypot zr zi) (m/hypot (Array/aget e 2) (Array/aget e 3)))))))

(defn- laguerre-terms
  "The logarithmic derivative `G = p'/p` and `H = G^2 - p''/p = sum 1/(z - z_k)^2` at `z` as a `Vec4` `[Gr Gi Hr Hi]`."
  ^Vec4 [^doubles e ^long n ^double zr ^double zi]
  (let [a (Vec2. (Array/aget e 0) (Array/aget e 1))
        b (Vec2. (Array/aget e 2) (Array/aget e 3))
        cc (Vec2. (Array/aget e 4) (Array/aget e 5))
        g (cplx/div b a)
        h2 (cplx/scale (cplx/div cc a) 2.0)]
    (if (m/== (Array/aget e 8) 1.0)
      ;; reversed polynomial, u = 1/z: G = n u - u^2 g,  H = u^2 (n - 2 u g - u^2 (h - g^2))
      (let [u (cplx/reciprocal (Vec2. zr zi))
            u2 (cplx/sq u)
            gg (cplx/sub (cplx/scale u n) (cplx/mult u2 g))
            hh (cplx/mult u2 (cplx/sub (cplx/sub (Vec2. n 0.0) (cplx/scale (cplx/mult u g) 2.0))
                                       (cplx/mult u2 (cplx/sub h2 (cplx/sq g)))))]
        (Vec4. (.x gg) (.y gg) (.x hh) (.y hh)))
      (let [hh (cplx/sub (cplx/sq g) h2)]
        (Vec4. (.x g) (.y g) (.x hh) (.y hh))))))

;; initial estimates

(defn- newton-polygon
  "Indices of the vertices of the upper convex hull of the points `(k, log |c_k|)`, in ascending order."
  [^doubles logs]
  (let [n (alength logs)
        hull (long-array n)]
    (loop [i (long 0) top (long 0)]
      (if (m/< i n)
        (if (m/> (Array/aget logs i) -1.0e300)
          (let [top (long (loop [top top]
                            (if (and (m/>= top 2)
                                     (let [k1 (Array/get hull (m/- top 2)) k2 (Array/get hull (m/dec top))
                                           e1 (Array/aget logs k1) e2 (Array/aget logs k2) e (Array/aget logs i)]
                                       (m/>= (m/- (m/* (m/- k2 k1) (m/- e e1)) (m/* (m/- e2 e1) (m/- i k1))) 0.0)))
                              (recur (m/dec top))
                              top)))]
            (Array/set hull top i)
            (recur (m/inc i) (m/inc top)))
          (recur (m/inc i) top))
        (vec (take top hull))))))

(defn- initial-estimates
  "Starting points on circles of the radii given by the Newton polygon; arrays `[re im]` of the `n` points."
  [^Poly p]
  (let [^doubles cr (.re p)
        ^doubles ci (.im p)
        n (.n p)
        logs (double-array (map (fn [^long k] (let [a (m/hypot (Array/aget cr k) (Array/aget ci k))]
                                                (if (m/> a 0.0) (m/log a) -1.0e301)))
                                (range (m/inc n))))
        hull (newton-polygon logs)
        zr (double-array n)
        zi (double-array n)
        two-pi m/TWO_PI
        offset (m// two-pi n)]
    (loop [edges (partition 2 1 hull) pos (long 0)]
      (when-let [[i j] (first edges)]
        (let [i (long i) j (long j)
              cnt (m/- j i)
              radius (m/exp (m/max -700.0 (m/min 700.0 (m// (m/- (Array/aget logs i) (Array/aget logs j)) cnt))))
              step (m// two-pi cnt)]
          (dotimes [t cnt]
            (let [angle (m/+ (m/* step (double t)) (m/* offset (double i)) SIGMA)]
              (Array/aset zr (m/+ pos t) (m/* radius (m/cos angle)))
              (Array/aset zi (m/+ pos t) (m/* radius (m/sin angle)))))
          (recur (rest edges) (m/+ pos cnt)))))
    [zr zi]))

;; the iteration

(defn- correction
  "The Laguerre step of the root `j` with the Aberth correction for the other roots; returns the step as a `Vec2` or `nil` when it is not finite."
  [^Vec4 terms ^long j ^doubles zr ^doubles zi]
  (let [n (alength zr)
        z (Vec2. (Array/aget zr j) (Array/aget zi j))
        ;; Aberth: the known root k removes 1/(z - z_k) from G and 1/(z - z_k)^2 from H
        ^Vec4 gh (loop [k (long 0) gr (.x terms) gi (.y terms) hr (.z terms) hi (.w terms)]
                   (if (m/== k n)
                     (Vec4. gr gi hr hi)
                     (if (m/== k j)
                       (recur (m/inc k) gr gi hr hi)
                       (let [t (cplx/reciprocal (cplx/sub z (Vec2. (Array/aget zr k) (Array/aget zi k))))
                             t2 (cplx/sq t)]
                         (recur (m/inc k) (m/- gr (.x t)) (m/- gi (.y t)) (m/- hr (.x t2)) (m/- hi (.y t2)))))))
        g (Vec2. (.x gh) (.y gh))
        h (Vec2. (.z gh) (.w gh))
        s (cplx/sqrt (cplx/scale (cplx/sub (cplx/scale h n) (cplx/sq g)) (m/dec n)))
        d1 (cplx/add g s)
        d2 (cplx/sub g s)
        den (if (m/> (cplx/abs d1) (cplx/abs d2)) d1 d2)]
    (when (m/> (cplx/abs den) 0.0)
      (let [^Vec2 delta (cplx/div (Vec2. n 0.0) den)]
        (when (and (m/valid-double? (.x delta)) (m/valid-double? (.y delta)))
          delta)))))

(defn- iterate-roots
  "Improves the roots in `zr`, `zi` of the polynomial `p` in place; returns the vector `[converged-flags iterations]` where a flag is `1` (converged), `0` (not) or `-1` (not finite)."
  [^Poly p ^doubles zr ^doubles zi ^long max-iters]
  (let [n (.n p)
        flags (long-array n)]
    (loop [it (long 0) done (long 0)]
      (if (or (m/== done n) (m/>= it max-iters))
        [flags it]
        (let [done (long (loop [j (long 0) done done]
                           (if (m/== j n)
                             done
                             (if (m/zero? (Array/get flags j))
                               (let [^doubles e (evaluate p (Array/aget zr j) (Array/aget zi j))
                                     small? (m/<= (m/hypot (Array/aget e 0) (Array/aget e 1)) (m/* m/MACHINE-EPSILON2 (Array/aget e 6)))]
                                 (if small?
                                   (do (Array/set flags j 1) (recur (m/inc j) (m/inc done)))
                                   (let [^Vec4 terms (laguerre-terms e n (Array/aget zr j) (Array/aget zi j))
                                         ^Vec2 delta (when (and (m/valid-double? (.x terms)) (m/valid-double? (.y terms)) (m/valid-double? (.z terms)) (m/valid-double? (.w terms)))
                                                   (correction terms j zr zi))]
                                     (if delta
                                       (do (Array/aset zr j (m/- (Array/aget zr j) (.x delta)))
                                           (Array/aset zi j (m/- (Array/aget zi j) (.y delta)))
                                           (recur (m/inc j) done))
                                       (do (Array/set flags j -1) (recur (m/inc j) (m/inc done)))))))
                               (recur (m/inc j) done)))))]
          (recur (m/inc it) done))))))

;; the result

(defn- sorted-roots
  "Sorts the `Root` records: for real coefficients the real roots ascending, then the pairs of complex conjugates by the real part, the root with the positive imaginary part first; for complex coefficients by `(re, im)`."
  [records real?]
  (if-not real?
    (vec (sort-by (fn [^Root r] [(.re r) (.im r)]) records))
    (let [real-roots (filter (fn [^Root r] (m/zero? (.im r))) records)
          others (remove (fn [^Root r] (m/zero? (.im r))) records)
          upper (vec (sort-by (fn [^Root r] [(m/- (.im r)) (.re r)]) (filter (fn [^Root r] (m/pos? (.im r))) others)))
          lower (filter (fn [^Root r] (m/neg? (.im r))) others)
          ;; the exact conjugate pair of a root of the upper half plane and its nearest root of the lower one
          conjugate-pair (fn [^Root a ^Root b]
                           (let [re (m/* 0.5 (m/+ (.re a) (.re b)))
                                 im (m/* 0.5 (m/- (.im a) (.im b)))
                                 berr (m/max (.berr a) (.berr b))
                                 cnd (m/max (.cnd a) (.cnd b))
                                 ok (and (.ok a) (.ok b))]
                             [(Root. re im berr cnd ok) (Root. re (m/- im) berr cnd ok)]))
          [pairs unmatched] (reduce (fn [[acc free] ^Root a]
                              (let [^Root best (when (seq free)
                                                 (apply min-key (fn [^Root b] (m/hypot (m/- (.re a) (.re b)) (m/+ (.im a) (.im b)))) free))]
                                (if best
                                  [(conj acc (conjugate-pair a best)) (remove #(identical? % best) free)]
                                  [(conj acc [a]) free])))
                            [[] lower] upper)]
      ;; a root of the lower half plane without a partner is kept as it is
      (vec (concat (sort-by (fn [^Root r] (.re r)) real-roots)
                   (apply concat (sort-by (fn [pair] (let [^Root r (first pair)] [(.re r) (m/- (.im r))]))
                                          (concat pairs (map vector unmatched)))))))))

(defn polyroot
  "Finds all the complex roots of a polynomial.

  Parameters:

  - `coeffs`: the coefficients `c0, c1, ..., cn` in ascending order of the power (the polynomial `c0 + c1 z + ... + cn z^n`, the order of [[fastmath.polynomials/polynomial]]), as a sequence or an array, or a polynomial object of `fastmath.polynomials`. A coefficient is a real number, a `Vec2` (the complex number `re + i im`) or a pair `[re im]`. All coefficients have to be finite.
  - `options` (map, optional):
      - `:max-iters` - limit of the iterations, the sweeps over the roots not converged yet, default `1200`.
      - `:info?` - when truthy, returns the map described below instead of the vector of roots. Default `false`.

  Returns a vector of the `n` roots as `Vec2` (the complex numbers), each root repeated according to its multiplicity (a multiple root is found with the accuracy of the order of the root of the machine epsilon of the multiplicity). Leading zero coefficients are dropped, trailing zero coefficients give exact zero roots, a constant gives no roots.

  Order: for real coefficients the real roots in ascending order, then the pairs of complex conjugates by the real part of the pair with the positive imaginary part first; the conjugates are exact and a root whose imaginary part is within the backward error of zero is real. For complex coefficients by the real part, then the imaginary part.

  With `:info?` returns the map `{:roots :backward-errors :condition-numbers :converged? :iterations :status}`. The vectors have the order of `:roots`. A backward error is `|p(z)| / sum |c_k| |z|^k` of the root, the smallest relative perturbation of the coefficients for which `z` is an exact root; a condition number `sum |c_k| |z|^k / (|z| |p'(z)|)` is the factor between the relative perturbation of the coefficients and the relative error of the root; the product bounds the relative error of the root. A zero root has the condition number `Infinity`. `:converged?` is a vector of booleans. The status is `:converged`, `:max-evals` (the iteration limit is reached) or `:diverged` (a value is not finite).

  Throws `ex-info` for an invalid argument (an empty or zero polynomial, an invalid coefficient: the data holds its `:index`, never the coefficients). A polynomial whose roots did not converge throws with the status `:max-evals` or `:diverged` and the number of the roots not converged as `:unconverged`; with `:info?` the status is returned.

  See also [[fastmath.solver/cubic]], [[fastmath.solver/quadratic]], [[fastmath.solver/find-roots]]."
  ([coeffs] (polyroot coeffs nil))
  ([coeffs options]
   (let [{:keys [^long max-evals info?]} (c/parse-options options)
         [^doubles re0 ^doubles im0] (parse-coefficients coeffs)
         real? (every? m/zero? (seq im0))
         len (alength re0)
         high (long (loop [k (m/dec len)] (if (and (m/>= k 0) (m/zero? (Array/aget re0 k)) (m/zero? (Array/aget im0 k))) (recur (m/dec k)) k)))]
     (when (m/< high 0)
       (c/throw-invalid-option :coefficients len "the polynomial is zero, all its coefficients are zero"))
     (let [low (long (loop [k (long 0)] (if (and (m/< k high) (m/zero? (Array/aget re0 k)) (m/zero? (Array/aget im0 k))) (recur (m/inc k)) k)))
           n (long (m/- high low))
           ^doubles cr (java.util.Arrays/copyOfRange re0 (int low) (int (m/inc high)))
           ^doubles ci (java.util.Arrays/copyOfRange im0 (int low) (int (m/inc high)))
           poly (Poly. cr ci n)
           zero-roots (repeat low (Root. 0.0 0.0 0.0 ##Inf true))
           [records iterations status]
           (cond
             (m/== n 0) [[] 0 :converged]
             (m/== n 1) (let [c0r (Array/aget cr 0) c0i (Array/aget ci 0) c1r (Array/aget cr 1) c1i (Array/aget ci 1)
                              ^Vec2 q (cplx/div (Vec2. c0r c0i) (Vec2. c1r c1i))
                              slope (m/* (m/hypot (.x q) (.y q)) (m/hypot c1r c1i))]
                          [[(Root. (m/- (.x q)) (m/- (.y q)) 0.0 (m// (m/+ (m/hypot c0r c0i) slope) slope) true)] 0 :converged])
             :else
             (let [[^doubles zr ^doubles zi] (initial-estimates poly)
                   [^longs flags iterations] (iterate-roots poly zr zi max-evals)
                   recs (mapv (fn [^long j]
                                (let [^doubles e (evaluate poly (Array/aget zr j) (Array/aget zi j))
                                      berr (m// (m/hypot (Array/aget e 0) (Array/aget e 1)) (Array/aget e 7))]
                                  (Root. (Array/aget zr j) (Array/aget zi j) berr (cond-number e n (Array/aget zr j) (Array/aget zi j)) (m/== (Array/get flags j) 1))))
                              (range n))
                   bad (count (remove (fn [^Root r] (.ok r)) recs))
                   status (cond (m/zero? bad) :converged
                                (some (fn [^long j] (m/== (Array/get flags j) -1)) (range n)) :diverged
                                :else :max-evals)]
               [recs iterations status]))
           ;; roots of real polynomials: an imaginary part within the rounding error of the polynomial is zero
           snapped (if real?
                     (mapv (fn [^Root r]
                             (if (and (not (m/zero? (.im r)))
                                      (let [^doubles e (evaluate poly (.re r) 0.0)]
                                        (m/<= (m// (m/hypot (Array/aget e 0) (Array/aget e 1)) (Array/aget e 7)) (m/* SNAP-EPSILON (m/+ 1.0 (m/* 3.8 n))))))
                               (Root. (.re r) 0.0 (.berr r) (.cnd r) (.ok r))
                               r))
                           records)
                     records)
           result (sorted-roots (vec (concat zero-roots snapped)) real?)
           out-roots (mapv (fn [^Root r] (Vec2. (.re r) (.im r))) result)
           unconverged (count (remove (fn [^Root r] (.ok r)) result))]
       (when (and (not= status :converged) (not info?))
         (c/throw-status status (str "polyroot did not converge: " unconverged " of " (count result) " roots")
                         {:degree n :unconverged unconverged :iterations iterations}))
       (if info?
         {:roots out-roots
          :backward-errors (mapv (fn [^Root r] (.berr r)) result)
          :condition-numbers (mapv (fn [^Root r] (.cnd r)) result)
          :converged? (mapv (fn [^Root r] (boolean (.ok r))) result)
          :iterations iterations
          :status status}
         out-roots)))))
