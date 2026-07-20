(ns fastmath.stats.logmean
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.divided :as div])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn logmean2
  "Logarithmic mean for two values."
  (^double [[^double x ^double y]] (logmean2 x y))
  (^double [^double x ^double y]
   (let [diff (m/- x y)]
     (if (m/zero? diff)
       x
       (m// diff (m/- (m/log x) (m/log y)))))))

(defn logmean3-mean-value
  "Generalized logarithmic mean, mean value for 3 values."
  (^double [[^double x ^double y ^double z]] (logmean3-mean-value x y z))
  (^double [^double x ^double y ^double z]
   (if (m/== x y z)
     x
     (m/sqrt (m// (m/* (m/- x y) (m/- y z) (m/- z x))
                  (m/* 2.0 (m/+ (m/* (m/- y z) (m/log x))
                                (m/* (m/- z x) (m/log y))
                                (m/* (m/- x y) (m/log z)))))))))

(defn logmean3-integral
  "Generalized logarithmic mean, integral method for 3 values."
  (^double [[^double x ^double y ^double z]] (logmean3-integral x y z))
  (^double [^double x ^double y ^double z]
   (if (m/== x y z)
     x
     (let [lx (m/log x)
           ly (m/log y)
           lz (m/log z)
           dxy (m/- lx ly)
           dyz (m/- ly lz)
           dzx (m/- lz lx)]
       (m// (m/+ (m/* x dyz) (m/* y dzx) (m/* z dxy))
            (m/* -0.5 dxy dyz dzx))))))

(defn logmean-mean-value
  "Generalized logarithmic mean - mean value / divided differences method."
  (^double [vs] (logmean-mean-value vs (count vs)))
  (^double [vs ^long cnt]
   (if (m/one? cnt)(first vs)
       (let [n- (m/dec cnt)]
         (m/pow (m/* n- (-> (div/divided m/log vs) m/abs))
                (m// -1.0 n-))))))

(defn middle
  ^double [xs]
  (let [fv (double (first xs))
        ^Vec2 mm (reduce (fn [^Vec2 curr ^double v]
                           (Vec2. (min (.x curr) v) (max (.y curr) v))) (Vec2. fv fv) (rest xs))]
    (m/* 0.5 (m/+ (.x mm) (.y mm)))))

;;;;;;;;;;;;;;;;;;;;;;;

;; Based on Claude/Opus4.8 approach: https://claude.ai/share/154297ad-95dc-4f56-8c76-b912410182ce

;; https://www.survo.fi/papers/logmean.pdf


;; Generalized logarithmic mean for n positive arguments, following
;; S. Mustonen, \"Logarithmic mean for several arguments\" (2002),
;;    https://www.survo.fi/papers/logmean.pdf

;;    For two arguments the logarithmic mean is

;;        L(x1, x2) = (x1 - x2) / log(x1/x2),   L(x, x) = x.

;;    The paper generalizes it to n arguments via the series expansion (eq. 3):

;;        L(x1,...,xn) = (n-1)! * Σ_{m≥0} h_m(u) / (n+m-1)! ,   u_i = log x_i,

;;    where h_m is the complete homogeneous symmetric polynomial of degree m,
;;    h_m(u) = Σ_{i1+...+in=m} u1^i1 ... un^in.  (The paper writes this polynomial
;;    as P(n,m); its C(n+m-1, m)·m! divisor is just (n+m-1)!/(n-1)!.)

;;    There is also a closed form (eq. 4),

;;        L(x1,...,xn) = (n-1)! Σ_i x_i / ∏_{j≠i} (log x_i - log x_j),

;;    which is the (n-1)th divided difference of exp at the points log x_i times
;;    (n-1)!.  It is exact and beautiful but, as §8 of the paper stresses, blows up
;;    numerically for n ≳ 14 (alternating huge terms) and is undefined when two
;;    arguments coincide.  Hence the series is the recommended workhorse.")

;; ;; ---------------------------------------------------------------------------
;; ;; Optimal implementation: the numerically stable series expansion (eq. 3)
;; ;; ---------------------------------------------------------------------------

;; (defn logarithmic-mean
;;   "Logarithmic mean of the positive numbers in `xs` (a seq of length n ≥ 1).

;;    Based on the series expansion (eq. 3), which the paper recommends for
;;    computation. Two refinements make it fast and robust for all realistic
;;    inputs:

;;    1. Recentering. By homogeneity, L(a·x) = a·L(x) (eq. 16), so we subtract the
;;       midrange of the logs, c = (min u + max u)/2, before summing and multiply
;;       e^c back at the end. This minimizes |log x_i|, which both maximizes the
;;       convergence rate of the series and improves accuracy.

;;    2. Overflow-free term recurrence. Rather than forming each h_m and the tiny
;;       coefficient (n-1)!/(n+m-1)! separately (their magnitudes span the whole
;;       double range for large n or wide data), we propagate the *terms*
;;       T_m = (n-1)! h_m / (n+m-1)! directly. Newton's identity for h_m,
;;       m·h_m = Σ_{k=1}^m p_k h_{m-k}, becomes after scaling the logs to [-1,1]

;;           m·T_m = Σ_{k=1}^m P_k · b_{m,k} · T_{m-k},
;;           b_{m,k} = s^k / ∏_{j=1}^k (n+m-j),

;;       where P_k are the bounded power sums of the rescaled logs (|P_k| ≤ n) and
;;       b is built as a running product, so nothing huge or tiny is materialized
;;       on its own.

;;    Repeated arguments and n = 1 are handled naturally (no division by log
;;    differences). Cost is O(n·M + M²) where M, the number of terms to converge,
;;    is small (~ proportional to half the spread of the logs).


(defn logmean-integral
  "Generalized logarithmic mean - integral method."
  (^double [xs] (logmean-integral xs (count xs)))
  (^double [xs ^long n] (logmean-integral xs n 10000 1.0e-15))
  (^double [xs ^long n ^long max-iters ^double tol]
   (let [us (v/log xs)
         c (middle us)
         vs (v/shift us (m/- c)) ;; L(ax)=aL(x)
         s (v/mx (map m/abs vs))]
     (if (m/zero? s)
       (Math/exp c)
       (let [ws (v/div vs s)] ;; scale to range [-1,1]
         (loop [m (long 1)
                wpow  ws
                psums [(reduce m/+ ws)]
                terms [1.0]
                sum   1.0
                tinyterms (long 0)]
           ;; m·T_m = Σ_{k=1..m} P_k · b_{m,k} · T_{m-k}
           (let [acc (double (loop [k (long 1)
                                    b (m// s (m/+ n m -1))
                                    a 0.0]
                               (if (m/> k m)
                                 a
                                 (recur (m/inc k)
                                        (m/* b (m// s (m/- (m/+ n m) (m/inc k))))
                                        (m/+ a (m/* (double (psums (m/dec k))) b (double (terms (m/- m k)))))))))
                 t (m// acc m)
                 nsum (m/+ sum t)
                 tiny?  (m/<= (m/abs t) (m/* tol (m/abs nsum)))
                 ntinyterms (if tiny? (m/inc tinyterms) 0)]
             (if (or (m/> ntinyterms 2) (m/> m max-iters))
               (m/* (Math/exp c) nsum)
               (let [nwpow (v/emult wpow ws)]
                 (recur (inc m)
                        nwpow
                        (conj psums (v/sum nwpow))
                        (conj terms t)
                        nsum
                        ntinyterms))))))))))

