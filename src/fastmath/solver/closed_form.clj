(ns fastmath.solver.closed-form
  "Real roots of polynomials of degree two and three.

  [[quadratic]] and [[cubic]] return the real roots as a vector padded with `NaN`: the real roots in ascending order, then `NaN` for each missing one. The coefficients are scaled by powers of two before any product is formed, so the results are correct over the whole range of doubles whenever the roots are representable, and `NaN` or an infinity otherwise. Polynomials of a lower degree (a leading coefficient equal to zero) are solved as such."
  (:require [fastmath.core :as m]
            [fastmath.solver.brent :as brent]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2 Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- exponent
  "The binary exponent of a non zero finite `x`: `x = mantissa * 2^e` with the mantissa in `[1, 2)`, subnormals included."
  ^long [^double x]
  (let [e (m/double-exponent x)]
    (if (m/> e -1023)
      e
      (m/- (m/double-exponent (m/* x 4.503599627370496E15)) 52))))

(defn- ratio
  "`n / d` of a finite `n` and a non zero finite `d` without intermediate overflow or underflow."
  ^double [^double n ^double d]
  (if (m/zero? n)
    0.0
    (let [en (exponent n)
          ed (exponent d)]
      (m/scalb (m// (m/scalb n (m/long-sub en)) (m/scalb d (m/long-sub ed))) (m/long-sub en ed)))))

(defn- ordered2
  "The pair in ascending order, with `NaN` last."
  ^Vec2 [^double x ^double y]
  (cond
    (m/nan? y) (Vec2. x y)
    (m/nan? x) (Vec2. y x)
    (m/<= x y) (Vec2. x y)
    :else (Vec2. y x)))

(defn quadratic
  "Real roots of the quadratic equation `a x^2 + b x + c = 0`.

  Parameters: `a`, `b`, `c` (numbers): the coefficients.

  Returns a `Vec2` `[r0 r1]` with the real roots in ascending order. A double root is returned twice. `NaN` fills the places of the missing roots: `[NaN NaN]` when the roots are not real (negative discriminant), when all the coefficients but `c` are zero (`c = 0` is an identity, `c != 0` has no solution, the two cases are not distinguished), and when any coefficient is not finite. For `a = 0` the equation is linear and the result is `[-c/b NaN]`.

  The coefficients are scaled by a power of two before the discriminant is formed, so the result does not overflow or underflow when the roots are representable. The roots are computed without cancellation, each with a relative error of a few units in the last place.

  See also [[cubic]], [[fastmath.solver/polyroot]]."
  ^Vec2 [^double a ^double b ^double c]
  (cond
    (not (and (m/valid-double? a) (m/valid-double? b) (m/valid-double? c))) (Vec2. ##NaN ##NaN)
    (m/zero? a) (if (m/zero? b)
                  (Vec2. ##NaN ##NaN)
                  (Vec2. (m/+ 0.0 (m/- (ratio c b))) ##NaN))
    (m/zero? c) (ordered2 0.0 (m/+ 0.0 (m/- (ratio b a))))
    (m/zero? b) (if (m/pos? (m/* (m/signum a) (m/signum c)))
                  (Vec2. ##NaN ##NaN)
                  (let [r (m// (m/sqrt (m/abs c)) (m/sqrt (m/abs a)))]
                    (Vec2. (m/- r) r)))
    :else
    (let [ea (exponent a)
          eb (exponent b)
          ec (exponent c)
          ;; t: the exponent of the larger of b and sqrt(a c), the scaled terms b' and 4 a'' c'' are O(1)
          t (m/long-max eb (long (m/ceil (m// (m/+ ea ec) 2.0))))
          a2 (m/scalb a (m/long-sub ea))
          b2 (m/scalb b (m/long-sub t))
          c2 (m/scalb c (m/long-sub (m/long-sub (m/long-mult 2 t) ea)))
          disc (m/difference-of-products b2 b2 (m/* 4.0 a2) c2)]
      (if (m/neg? disc)
        (Vec2. ##NaN ##NaN)
        (let [q (m/* -0.5 (m/+ b2 (m/copy-sign (m/sqrt disc) b2)))
              r0 (m/scalb (m// q a2) (m/long-sub t ea))
              r1 (m/scalb (m// (m/scalb c (m/long-sub ec)) q) (m/long-sub ec t))]
          (ordered2 r0 r1))))))

;; cubic

(def ^{:private true :const true :tag 'double} NOISE (m/* 64.0 m/MACHINE-EPSILON2))
(def ^{:private true :const true :tag 'double} WINDOW 1.125899906842624E15)
(def ^{:private true :const true :tag 'double} -WINDOW -1.125899906842624E15)

(defn- value3
  "Value of `d + z (c + z (b + z a))`, the coefficients in the `Vec4` `cf` as `[d c b a]`."
  ^double [^Vec4 cf ^double z]
  (m/+ (.x cf) (m/* z (m/+ (.y cf) (m/* z (m/+ (.z cf) (m/* z (.w cf))))))))

(defn- noise3
  "The rounding error scale of the value at `z`: 64 eps times the sum of the absolute values of the terms."
  ^double [^Vec4 cf ^double z]
  (let [az (m/abs z)]
    (m/* NOISE (m/+ (m/abs (.x cf)) (m/* az (m/+ (m/abs (.y cf)) (m/* az (m/+ (m/abs (.z cf)) (m/* az (m/abs (.w cf)))))))))))

(defn- roots-in-window
  "Real roots `z` with `|z| < WINDOW` of the cubic with the coefficients `cf = [d c b a]` as a vector.

  The critical points of the cubic split the window into monotone parts and every part with a sign change holds one root, found by Brent's method. A critical value which is zero within the rounding error (a multiple root) gives the critical point as a double (or triple) root."
  [^Vec4 cf]
  (let [value (fn ^double [^double z] (value3 cf z))
        solve (fn ^double [^double lo ^double hi] (brent/brent value lo hi {:absolute-accuracy 0.0}))
        [^double c1 ^double c2] (quadratic (m/* 3.0 (.w cf)) (m/* 2.0 (.z cf)) (.y cf))
        in? (fn [^double c] (and (not (m/nan? c)) (m/< (m/abs c) WINDOW)))
        crit (cond
               (not (in? c1)) (if (in? c2) [c2] [])
               (not (in? c2)) [c1]
               (m/== c1 c2) [c1]
               :else [c1 c2])
        pts (vec (concat [-WINDOW] crit [WINDOW]))
        svals (mapv value pts)
        ;; the interior points are the critical points: a value within the rounding error means a multiple root
        last-index (m/dec (long (count pts)))
        at-zero (mapv (fn [^long i] (and (m/pos? i) (m/< i last-index)
                                        (m/<= (m/abs (double (svals i))) (noise3 cf (double (pts i))))))
                      (range (count pts)))
        multiple (cond
                   (and (m/one? (count crit)) (at-zero 1)) [(pts 1) (pts 1) (pts 1)]
                   :else (mapcat (fn [i] (when (at-zero i) [(pts i) (pts i)])) (range 1 last-index)))
        simple (keep (fn [^long i]
                       (let [lo (double (pts i)) hi (double (pts (m/inc i)))
                             vl (double (svals i)) vh (double (svals (m/inc i)))]
                         (cond
                           (or (at-zero i) (at-zero (m/inc i))) nil
                           (m/zero? vl) lo
                           (m/zero? vh) hi
                           (not (m/== (m/signum vl) (m/signum vh))) (solve lo hi)
                           :else nil)))
                     (range last-index))]
    (vec (concat multiple simple))))

(defn- hull-scales
  "The binary exponents of the magnitudes of the roots of a polynomial: the slopes of the upper convex hull of the points `(k, exponent of the coefficient k)`, rounded, sorted and distinct.

  `exps` is the vector of the exponents of the coefficients in the order of the power of `x`, `nil` for a zero coefficient."
  [exps]
  (let [n (count exps)
        hk (int-array n)
        he (int-array n)]
    (loop [i (long 0) top (long 0)]
      (if (m/< i n)
        (if-let [e (exps i)]
          (let [e (long e)
                top (long (loop [top top]
                            (if (and (m/>= top 2)
                                     (let [k1 (Array/get hk (m/- top 2)) e1 (Array/get he (m/- top 2))
                                           k2 (Array/get hk (m/dec top)) e2 (Array/get he (m/dec top))]
                                       (m/>= (m/- (m/* (m/- k2 k1) (m/- e e1)) (m/* (m/- e2 e1) (m/- i k1))) 0)))
                              (recur (m/dec top))
                              top)))]
            (Array/set hk top i)
            (Array/set he top e)
            (recur (m/inc i) (m/inc top)))
          (recur (m/inc i) top))
        (vec (sort (distinct (map (fn [^long j]
                                    (m/round (m// (m/- (Array/get he j) (Array/get he (m/inc j)))
                                                  (m/- (Array/get hk (m/inc j)) (Array/get hk j)))))
                                  (range (m/dec top))))))))))

(defn- roots-at-scale
  "The real roots of the cubic `cf = [d c b a]` whose magnitude is nearest to the scale `sigma` among `scales`, as numbers."
  [^Vec4 cf exps scales ^long sigma]
  (let [m (long (apply m/max (keep-indexed (fn [^long k e] (when e (m/+ (long e) (m/* k sigma)))) exps)))
        cz (v/vec4 (map-indexed (fn [^long k ^double c] (m/scalb c (m/long-sub (m/long-mult k sigma) m))) cf))
        nearest (fn [^long ex] (first (sort-by (fn [^long sg] (m/abs (m/- ex sg))) scales)))]
    (keep (fn [^double z]
            (let [x (m/scalb z sigma)
                  ex (cond (m/inf? x) 3000
                           (m/zero? x) -3000
                           :else (exponent x))]
              (when (m/== sigma (long (nearest ex))) x)))
          (roots-in-window cz))))

(defn cubic
  "Real roots of the cubic equation `a x^3 + b x^2 + c x + d = 0`.

  Parameters: `a`, `b`, `c`, `d` (numbers): the coefficients.

  Returns a `Vec3` `[r0 r1 r2]` with the real roots in ascending order, `NaN` after the last real root: one real root gives `[r NaN NaN]`, three real roots (a multiple root is returned as many times as its multiplicity) give all three places. `[NaN NaN NaN]` is returned when any coefficient is not finite. For `a = 0` the equation is quadratic and the result is the result of [[quadratic]] followed by `NaN`. For `d = 0` one root is exactly `0.0`.

  A polynomial which differs from one with a multiple real root by no more than the rounding error of its coefficients is reported with the multiple root (a double root is returned twice, a triple root three times), with a relative accuracy of the root of the order of the square root, respectively the cube root, of the machine epsilon.

  The magnitudes of the roots are estimated from the exponents of the coefficients (the Newton polygon) and every group of roots is solved in its own scale: the critical points of the scaled cubic split the line into monotone parts, and the root of every part with a sign change is found by Brent's method. No root is deflated from the polynomial and the coefficients are scaled by powers of two, so the roots keep a small componentwise backward error also when they differ by many orders of magnitude, and the result is correct over the whole range of doubles whenever the roots are representable. A root beyond the range is `Infinity`, a root below the range is `0.0`.

  See also [[quadratic]], [[fastmath.solver/polyroot]]."
  ^Vec3 [^double a ^double b ^double c ^double d]
  (cond
    (not (and (m/valid-double? a) (m/valid-double? b) (m/valid-double? c) (m/valid-double? d))) (Vec3. ##NaN ##NaN ##NaN)
    (m/zero? a) (let [[^double r0 ^double r1] (quadratic b c d)] (Vec3. r0 r1 ##NaN))
    (m/zero? d) (let [[^double r0 ^double r1] (quadratic a b c)
                      rs (vec (sort (remove m/nan? [r0 r1 0.0])))]
                  (Vec3. (double (rs 0)) (double (get rs 1 ##NaN)) (double (get rs 2 ##NaN))))
    :else
    (let [cf (Vec4. d c b a)
          exps (mapv (fn [^double x] (when-not (m/zero? x) (exponent x))) [d c b a])
          scales (hull-scales exps)
          rs (vec (sort (mapcat (fn [sigma] (roots-at-scale cf exps scales sigma)) scales)))]
      (Vec3. (double (get rs 0 ##NaN)) (double (get rs 1 ##NaN)) (double (get rs 2 ##NaN))))))
