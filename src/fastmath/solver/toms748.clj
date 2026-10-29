(ns fastmath.solver.toms748
  "Algorithm 748 (TOMS748) of Alefeld, Potra and Shi for a root of a function on an interval with a sign change.

  The method encloses the root in a shrinking bracket and uses inverse cubic interpolation (a quadratic step when the values are too close), a double length secant step and, when the bracket does not shrink fast enough, a bisection step. It needs the fewest evaluations of the bracketing methods in this library on smooth functions (convergence order about 1.65 per evaluation) and keeps the guarantees of the bisection."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [fastmath.java Array]
           [fastmath.solver.common IEvaluator Problem]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'double} EPS2 (m/* 2.0 m/MACHINE-EPSILON2))
(def ^{:private true :const true :tag 'double} MIN-DIFF 1.0E-306)

;; the state is a double array: a, b, fa, fb, d, fd, e, fe; (a, b) is the bracket, d the point dropped last and e the one dropped before

(def ^{:private true :const true :tag 'long} A 0)
(def ^{:private true :const true :tag 'long} B 1)
(def ^{:private true :const true :tag 'long} FA 2)
(def ^{:private true :const true :tag 'long} FB 3)
(def ^{:private true :const true :tag 'long} D 4)
(def ^{:private true :const true :tag 'long} FD 5)
(def ^{:private true :const true :tag 'long} E 6)
(def ^{:private true :const true :tag 'long} FE 7)

(defn- safe-div
  ^double [^double num ^double denom ^double r]
  (if (m/< (m/abs denom) 1.0) 
    (if (m/>= (m/abs num) (m/* (m/abs denom) Double/MAX_VALUE)) r (m// num denom))
    (m// num denom)))

(defn- secant-interpolate
  ^double [^double a ^double b ^double fa ^double fb]
  (let [tol m/MACHINE-EPSILON10
        c (m/- a (m/* (m// fa (m/- fb fa)) (m/- b a)))]
    (if (or (m/<= c (m/+ a (m/* (m/abs a) tol))) (m/>= c (m/- b (m/* (m/abs b) tol))) (m/nan? c))
      (m/* 0.5 (m/+ a b))
      c)))

(defn- quadratic-interpolate
  "The root of the quadratic through `a`, `b`, `d` (Newton steps starting from an end of the bracket) with the state `s`."
  ^double [^doubles s ^long steps]
  (let [a (Array/aget s A) b (Array/aget s B) d (Array/aget s D)
        fa (Array/aget s FA) fb (Array/aget s FB) fd (Array/aget s FD)
        bb (safe-div (m/- fb fa) (m/- b a) Double/MAX_VALUE)
        aa0 (safe-div (m/- fd fb) (m/- d b) Double/MAX_VALUE)
        aa (safe-div (m/- aa0 bb) (m/- d a) 0.0)]
    (if (m/zero? aa)
      (secant-interpolate a b fa fb)
      (let [c0 (if (m/pos? (m/* (m/signum aa) (m/signum fa))) a b)
            c (double (loop [i (long 1) c c0]
                (if (m/> i steps)
                  c
                  (recur (m/inc i)
                         (m/- c (safe-div (m/+ fa (m/* (m/+ bb (m/* aa (m/- c b))) (m/- c a)))
                                          (m/+ bb (m/* aa (m/- (m/- (m/* 2.0 c) a) b)))
                                          (m/- (m/+ 1.0 c) a)))))))]
        (if (or (m/<= c a) (m/>= c b) (m/nan? c))
          (secant-interpolate a b fa fb)
          c)))))

(defn- cubic-interpolate
  "The root of the inverse cubic interpolation through `a`, `b`, `d`, `e` with the state `s`."
  ^double [^doubles s]
  (let [a (Array/aget s A) b (Array/aget s B) d (Array/aget s D) e (Array/aget s E)
        fa (Array/aget s FA) fb (Array/aget s FB) fd (Array/aget s FD) fe (Array/aget s FE)
        q11 (m// (m/* (m/- d e) fd) (m/- fe fd))
        q21 (m// (m/* (m/- b d) fb) (m/- fd fb))
        q31 (m// (m/* (m/- a b) fa) (m/- fb fa))
        d21 (m// (m/* (m/- b d) fd) (m/- fd fb))
        d31 (m// (m/* (m/- a b) fb) (m/- fb fa))
        q22 (m// (m/* (m/- d21 q11) fb) (m/- fe fb))
        q32 (m// (m/* (m/- d31 q21) fa) (m/- fd fa))
        d32 (m// (m/* (m/- d31 q21) fd) (m/- fd fa))
        q33 (m// (m/* (m/- d32 q22) fa) (m/- fe fa))
        cc (m/+ a q31 q32 q33)]
    (if (or (m/<= cc a) (m/>= cc b) (m/nan? cc))
      (quadratic-interpolate s 3)
      cc)))

(defn- close-values?
  [^doubles s]
  (let [fa (Array/aget s FA) fb (Array/aget s FB) fd (Array/aget s FD) fe (Array/aget s FE)]
    (or (m/< (m/abs (m/- fa fb)) MIN-DIFF) (m/< (m/abs (m/- fa fd)) MIN-DIFF) (m/< (m/abs (m/- fa fe)) MIN-DIFF)
        (m/< (m/abs (m/- fb fd)) MIN-DIFF) (m/< (m/abs (m/- fb fe)) MIN-DIFF) (m/< (m/abs (m/- fd fe)) MIN-DIFF))))

(defn- bracket!
  "Evaluates `f` at `c` (moved away from the ends) and shrinks the bracket in the state `s`; the dropped end goes to `d`."
  [^IEvaluator ev ^doubles s ^double c]
  (let [a (Array/aget s A) b (Array/aget s B) fa (Array/aget s FA) fb (Array/aget s FB)
        c (double (cond (m/< (m/- b a) (m/* 2.0 EPS2 (m/abs a))) (m/+ a (m/* 0.5 (m/- b a)))
                (m/<= c (m/+ a (m/* (m/abs a) EPS2))) (m/+ a (m/* (m/abs a) EPS2))
                (m/>= c (m/- b (m/* (m/abs b) EPS2))) (m/- b (m/* (m/abs b) EPS2))
                :else c))
        c (if (and (m/< c b) (m/> c a)) c (m/+ a (m/* 0.5 (m/- b a))))
        fc (.eval ev c)]
    (cond
      (m/zero? fc) (do (Array/aset s A c) (Array/aset s FA 0.0) (Array/aset s D 0.0) (Array/aset s FD 0.0))
      (m/== (m/signum fa) (m/signum fc)) (do (Array/aset s D a) (Array/aset s FD fa) (Array/aset s A c) (Array/aset s FA fc))
      :else (do (Array/aset s D b) (Array/aset s FD fb) (Array/aset s B c) (Array/aset s FB fc)))))

(defn- core
  ^double [^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)
        s (double-array [(.a p) (.b p) (.fa p) (.fb p) 0.0 0.0 0.0 0.0])
        done? (fn [] (or (m/zero? (Array/aget s FA))
                         (m/<= (m/* 0.5 (m/- (Array/aget s B) (Array/aget s A)))
                               (c/tolerance atol rtol (m/+ (m/* 0.5 (Array/aget s A)) (m/* 0.5 (Array/aget s B)))))))
        step! (fn [^double c] (bracket! ev s c))
        shift! (fn [] (Array/aset s E (Array/aget s D)) (Array/aset s FE (Array/aget s FD)))]
    ;; first step: secant, second: quadratic, then the cycle
    (step! (secant-interpolate (Array/aget s A) (Array/aget s B) (Array/aget s FA) (Array/aget s FB)))
    (when-not (done?)
      (shift!)
      (step! (quadratic-interpolate s 2)))
    (loop []
      (when-not (done?)
        (let [width0 (m/- (Array/aget s B) (Array/aget s A))
              interp! (fn [^long qsteps cubic?]
                        (step! (double (if (and cubic? (not (close-values? s)))
                                 (cubic-interpolate s)
                                 (quadratic-interpolate s qsteps)))))]
          (interp! 2 true)
          (shift!)
          (when-not (done?)
            (interp! 3 true)
            (when-not (done?)
              ;; double length secant step
              (let [a (Array/aget s A) b (Array/aget s B) fa (Array/aget s FA) fb (Array/aget s FB)
                    low? (m/< (m/abs fa) (m/abs fb))
                    u (if low? a b)
                    fu (if low? fa fb)
                    cc0 (m/- u (m/* 2.0 (m/* (m// fu (m/- fb fa)) (m/- b a))))
                    cc (if (or (m/nan? cc0) (m/> (m/abs (m/- cc0 u)) (m/* 0.5 (m/- b a)))) (m/+ a (m/* 0.5 (m/- b a))) cc0)]
                (shift!)
                (step! cc)
                (when-not (done?)
                  (when-not (m/< (m/- (Array/aget s B) (Array/aget s A)) (m/* 0.5 width0))
                    (shift!)
                    (step! (m/+ (Array/aget s A) (m/* 0.5 (m/- (Array/aget s B) (Array/aget s A))))))
                  (recur))))))))
    (let [fa (Array/aget s FA) fb (Array/aget s FB)]
      (if (m/<= (m/abs fa) (m/abs fb)) (Array/aget s A) (Array/aget s B)))))

(defn toms748
  "Finds a root of `f` on `[lower, upper]` with Algorithm 748 (TOMS748).

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double, the end of the final bracket with the smaller absolute value of the function. A root at an end of the interval is returned after at most two evaluations.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.brent/brent]], [[fastmath.solver.bisection/bisection]]."
  ([f lower upper] (toms748 f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core false f lower upper options)))
