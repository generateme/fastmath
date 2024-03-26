;; Radial Basis Functions
(ns fastmath.kernel.rbf
  "Radial Basis Functions.

  Collection of RBF kernels.

  Each kernel supports `scale` parameter."
  (:require [fastmath.core :as m]
            [fastmath.calculus.quadrature :as calc]
            [fastmath.calculus.finite :as finite]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Literature
;; http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol4/Choosing_Basis_Functions_and_Shape_Parameters.pdf
;; https://www.math.ucdavis.edu/~saito/data/jim/buhmann-actanumerica.pdf
;; https://www.math.utah.edu/vigre/pubs/pdfs/flw08.pdf
;; https://num.math.uni-goettingen.de/schaback/teaching/sc.pdf

(defn linear
  "Linear kernel"
  ([] (linear nil))
  ([{:keys [^double scale]
     :or {scale 1.0}}] (fn [^double x] (m/abs (/ x scale)))))

(defn gaussian
  "Gaussian kernel"
  ([] (gaussian nil))
  ([{:keys [^double scale]
     :or {scale 1.0}}] (fn [^double x] (m/exp (- (m/sq (/ x scale)))))))

;; RBF kernels
;; [1] https://www.math.unipd.it/~demarchi/RBF/LectureNotes.pdf
;; [2] Meshfree Approximation Methods with MATLAB, Gregory E. Fasshauer

;; [1] p.32

(defn truncated-power
  "Truncated power kernel with `k` exponent."
  ([] (truncated-power nil))
  ([{:keys [^double k ^double scale]
     :or {scale 1.0}}]
   (cond
     (== k m/THIRD) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0)
                                                                           (m/cbrt (- 1.0 ex)) 0.0)))
     (== k 0.5) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0) (m/sqrt (- 1.0 ex)) 0.0)))
     (== k 1.0) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0) (- 1.0 ex) 0.0)))
     (== k 2.0) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0) (m/sq (- 1.0 ex)) 0.0)))
     (== k 3.0) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0) (m/cb (- 1.0 ex)) 0.0)))
     :else (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (<= ex 1.0) (m/pow (- 1.0 ex) k) 0.0))))))

;; [1] p.33

(defn gaussians-laguerre
  "Gaussians-Laguerre kernel with `dimension` (double) and `degree` (long)"
  ([] (gaussians-laguerre nil))
  ([{:keys [^double dimension ^long degree ^double scale]
     :or {scale 1.0}}]
   (let [hdimension (* 0.5 dimension)]
     (fn ^double [^double x]
       (let [ex (m/sq (/ x scale))]
         (* (m/exp (- ex)) (m/laguerre-polynomials degree hdimension ex)))))))

;; [1] p.33

(defn poisson
  "Poisson (or Bessel) kernel for `d` (double)"
  ([] (poisson nil))
  ([{:keys [^double d ^double scale]
     :or {scale 1.0}}]
   (cond
     (== d 2) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (m/bessel-j 0 ex)))
     (== d 3) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (* 0.7978845608028654
                                                                      (if (m/near-zero? ex 1.0e-13) 1.0 (/ (m/sin ex) ex)))))
     (== d 4) (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (if (m/near-zero? ex 1.0e-13) 0.5
                                                                       (/ (m/bessel-j 1 ex) ex))))
     :else (let [p (dec (* 0.5 d))]
             (fn ^double [^double x] (let [ex (m/max 1.0e-13 (m/abs (/ x scale)))]
                                      (/ (m/bessel-j p ex) (m/pow ex p))))))))


;; [1] p. 34

;; note: matern52, from formula it's 1/3 of the example from paper

(defn matern
  "Matern kernel for odd `order` based on Bessel function of the second kind of `order/2` order.

  For example for kernel 5/2 use `order` equal 5"
  ([] (matern nil))
  ([{:keys [^long order ^double scale]
     :or {scale 1.0}}]
   (let [v (* 0.5 order)
         rg (m/exp (- (* (- 1.0 v) m/LN2) (m/log-gamma v)))]
     (fn ^double [^double x]
       (let [ex (m/abs (/ x scale))]
         (if (< ex 1.0e-20)
           1.0
           (* rg (m/pow ex v) (m/bessel-k-half order ex))))))))

;; [1] p. 43

(defn- infer-negate
  ^double [^double beta negate?]
  (cond
    (= :auto negate?) (if (even? (int (m/ceil (* 0.5 beta)))) 1.0 -1.0)
    negate? -1.0
    :else 1.0))

(defn generalized-multiquadratic
  "Generalized multiquadratic kernel with `beta` exponent.

  Can be `negated` if `:negate?` is set to true or default `:auto` (negate for cases when ceil of half beta is even)."
  ([] (generalized-multiquadratic nil))
  ([{:keys [^double beta negate? ^double scale]
     :or {scale 1.0 negate? :auto}}]
   (let [fact (infer-negate beta negate?)]
     (cond
       (== beta m/THIRD) (fn ^double [^double x] (* fact (m/cbrt (inc (m/sq (/ x scale))))))
       (== beta (- m/THIRD)) (fn ^double [^double x] (* fact (m/cbrt (/ (inc (m/sq (/ x scale)))))))
       (== beta 0.5) (fn ^double [^double x] (* fact (m/sqrt (inc (m/sq (/ x scale))))))
       (== beta -0.5) (fn ^double [^double x] (* fact (m/sqrt (/ (inc (m/sq (/ x scale)))))))
       (== beta 1.0) (fn ^double [^double x] (* fact (inc (m/sq (/ x scale)))))
       (== beta -1.0) (fn ^double [^double x] (* fact (/ (inc (m/sq (/ x scale))))))
       (== beta 2.0) (fn ^double [^double x] (* fact (m/sq (inc (m/sq (/ x scale))))))
       (== beta -2.0) (fn ^double [^double x] (* fact (m/sq (/ (inc (m/sq (/ x scale)))))))
       (== beta 3.0) (fn ^double [^double x] (* fact (m/cb (inc (m/sq (/ x scale))))))
       (== beta -3.0) (fn ^double [^double x] (* fact (m/cb (/ (inc (m/sq (/ x scale)))))))
       :else (fn ^double [^double x] (* fact (m/pow (inc (m/sq (/ x scale))) beta)))))))

;; [1] p. 43

(defn radial-powers
  "Radial powers kernel with `beta` exponent.

  Can be `negated` if `:negate?` is set to true or default `:auto` (negate for cases when ceil of half beta is even)."
  ([] (radial-powers nil))
  ([{:keys [^double beta negate? ^double scale]
     :or {scale 1.0 negate? :auto}}]
   (let [fact (infer-negate beta negate?)]
     (cond
       (== beta 1.0) (fn ^double [^double x] (* fact (m/abs (/ x scale))))
       (== beta 3.0) (fn ^double [^double x] (* fact (m/cb (m/abs (/ x scale)))))
       :else (fn ^double [^double x] (* fact (m/pow (m/abs (/ x scale)) beta)))))))

(defn thin-plate-splines
  "Thin-plate splines kernel with `beta` exponent.

  Can be `negated` if `:negate?` is set to true or default `:auto` (negate for cases when ceil of half of beta is even)."
  ([] (thin-plate-splines nil))
  ([{:keys [^double beta negate? ^double scale]
     :or {scale 1.0 negate? :auto}}]
   (let [beta2 (* 2.0 beta)
         fact (infer-negate beta negate?)]
     (fn ^double [^double x]
       (let [ex (m/abs (/ x scale))]
         (if (zero? ex) 0.0 (* fact (m/pow ex beta2) (m/log ex))))))))

;; [2] p.44

(defn- whittaker-integration [^double alpha ^double k ^double beta]
  (let [k- (dec k)]
    (if (zero? alpha)
      (fn ^double [^double r]
        (let [f (if (zero? r)
                  (fn ^double [^double t] (m/exp (- (* t beta))))
                  (fn ^double [^double t] (* (m/pow (max 0.0 (- 1.0 (* r t))) k-)
                                            (m/exp (- (* t beta))))))]
          (calc/gk-quadrature f 0.0 ##Inf {:max-iters 1000 :abs 1.0e-6 :rel 1.0e-6})))
      (fn ^double [^double r]
        (let [f (if (zero? r)
                  (fn ^double [^double t] (* (m/pow t alpha)
                                            (m/exp (- (* t beta)))))
                  (fn ^double [^double t] (* (m/pow (max 0.0 (- 1.0 (* r t))) k-)
                                            (m/pow t alpha)
                                            (m/exp (- (* t beta))))))]
          (calc/gk-quadrature f 0.0 ##Inf {:max-iters 1000 :abs 1.0e-6 :rel 1.0e-6}))))))

(defn whittaker
  "Wittaker kernel family.

  Closed form expressions are used for `alpha` equal 0 or 1 and `k` equal 2 or 3. Other cases use numerical integration"
  ([] (whittaker nil))
  ([{:keys [^double alpha ^double k ^double beta ^double scale]
     :or {scale 1.0}}]
   (condp = [alpha k]
     [0.0 2.0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                        (/ (+ (- beta ex)
                                              (* ex (m/exp (- (/ beta ex))))) beta)))
     [1.0 2.0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))
                                            ex2 (* 2.0 ex)]
                                        (/ (+ (- beta ex2)
                                              (* (+ beta ex2) (m/exp (- (/ beta ex))))) beta)))
     [0.0 3.0] (let [b (m/sq beta)]
                 (fn ^double [^double x] (let [ex (m/abs (/ x scale))
                                              ex2 (* 2.0 (m/sq ex))]
                                          (/ (+ (- b (* 2.0 beta ex)) ex2
                                                (- (* ex2 (m/exp (- (/ beta ex)))))) b))))
     [1.0 3.0] (let [b (m/sq beta)]
                 (fn ^double [^double x] (let [ex (m/abs (/ x scale))
                                              ex2 (* 6.0 (m/sq ex))
                                              bex (* 2.0 beta ex)]
                                          (/ (+ (- b (* 2.0 bex)) ex2
                                                (- (* (+ bex ex2) (m/exp (- (/ beta ex)))))) b))))
     (let [w (whittaker-integration alpha k beta)
           m (double (w 0.0))]
       (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                (/ ^double (w ex) m)))))))

;; [2] p.131

(defn shifted-surface-splines
  "Shifter surface splines kernel family with parameters `s` and `beta`"
  ([] (shifted-surface-splines nil))
  ([{:keys [^long s ^double beta ^double scale]
     :or {scale 1.0}}]
   (let [t (- beta (/ s 2))
         ct (long (m/ceil t))
         fact (if (even? ct) 1.0 -1.0)
         fact- (- fact)]
     (assert (pos? t) "beta should be greater than half of s.")
     (if (odd? s)
       (fn ^double [^double x]
         (* fact (m/pow (inc (m/sq (/ x scale))) t)))
       (fn ^double [^double x]
         (let [ex (inc (m/sq (/ x scale)))]
           (* fact- (m/pow ex t) (m/sqrt (m/log ex)))))))))

;; [2] p.88
;; normalized to have 1.0 at 0.0

(defn wendland
  "Wendland kernel family with parameters `s` and `k`. Only `k` < 5 is supported."
  ([] (wendland nil))
  ([{:keys [^double s ^long k ^double scale]
     :or {scale 1.0}}]
   (let [l (inc (+ k (* s 0.5)))
         p (+ l k)]
     (case (int k)
       0 (fn ^double [^double x] (let [r (m/abs (/ x scale))]
                                  (if (<= r 1.0) (m/pow (- 1.0 (m/abs (/ x scale))) l) 0.0)))
       1 (fn ^double [^double x] (let [r (m/abs (/ x scale))]
                                  (if (<= r 1.0)
                                    (* (m/pow (- 1.0 r) p) (inc (* p r))) 0.0)))
       2 (let [f2 (m/mevalpoly l 3.0 4.0 1.0)
               f1 (m/mevalpoly l 6.0 3.0)]
           (fn ^double [^double x] (let [r (m/abs (/ x scale))]
                                    (if (<= r 1.0)
                                      (/ (* (m/pow (- 1.0 r) p)
                                            (m/mevalpoly r 3.0 f1 f2)) 3.0) 0.0))))
       3 (let [f3 (m/mevalpoly l 15.0 23.0 9.0 1.0)
               f2 (m/mevalpoly l 45.0 36.0 6.0)
               f1 (m/mevalpoly l 45.0 15.0)]
           (fn ^double [^double x] (let [r (m/abs (/ x scale))]
                                    (if (<= r 1.0)
                                      (/ (* (m/pow (- 1.0 r) p)
                                            (m/mevalpoly r 15.0 f1 f2 f3)) 15.0) 0.0))))
       4 (let [f4 (m/mevalpoly l 105.0 176.0 86.0 16.0 1.0)
               f3 (* 5.0 (m/mevalpoly l 84.0 85.0 24.0 2.0))
               f2 (* 45.0 (m/mevalpoly l 14.0 8.0 1.0))
               f1 (* 105 (m/mevalpoly l 4.0 1.0))]
           (fn ^double [^double x] (let [r (m/abs (/ x scale))]
                                    (if (<= r 1.0)
                                      (/ (* (m/pow (- 1.0 r) p)
                                            (m/mevalpoly r 105.0 f1 f2 f3 f4)) 105.0) 0.0))))))))

;; [2] p. 91

(defn gneiting
  "Gneiting kernel family with parameters `s` and `l`."
  ([] (gneiting nil))
  ([{:keys [^double s ^double l ^double scale]
     :or {scale 1.0}}]
   (let [f (- (/ (* (inc l) (+ l s 2.0)) s))]
     (fn [^double x]
       (let [ex (m/abs (/ x scale))]
         (if (<= ex 1.0) (* (m/pow (- 1.0 ex) l)
                            (m/mevalpoly ex 1.0 l f)) 0.0))))))

;;

;; https://www.researchgate.net/publication/226871599_Compactly_supported_positive_definite_radial_function
;; https://clojurians.zulipchat.com/#narrow/stream/151924-data-science/topic/.E2.9C.94.20Calculus.20help

(defn- wu-integration [^double l]
  (fn ^double [^double r]
    (let [rr (* 2.0 r)]
      (calc/gk-quadrature (fn [^double t] (* (m/pow (- 1.0 (* t t)) l)
                                            (m/pow (- 1.0 (m/sq (- rr t))) l)))
                          (dec rr) 1.0))))

(defn- wu-differentiation [f ^long k]
  (if (zero? k)
    f
    (let [Dk (finite/derivative f) ;; derivative
          fDk (fn ^double [^double r ^double h] (let [hr (+ r h)]
                                                 (/ ^double (Dk hr) (- hr)))) ;; D operator
          eDk (finite/extrapolate fDk {:power 2.0})] ;; extrapolation
      (recur (fn ^double [^double r]
               (let [r (m/abs r)]
                 (if (< r 1.0e-3)
                   (eDk r) ;; extrapolate at 0 for small `r`
                   (fDk r 0.0))))
             (dec k)))))

(defn wu
  "Wu kernel family with parameters `l` (double) and `k (long)`.

   Integer `l` from 0 to 3 and `k` from 0 to `l` has closed from expressions. Other cases use numerical integration and differentiation and might be unstable."
  ([] (wu nil))
  ([{:keys [^double l ^long k ^double scale]
     :or {scale 1.0}}]
   (condp = [l k]
     [0.0 0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))] (max 0.0 (- 1.0 ex))))
     [1.0 0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (* (m/cb (max 0.0 (- 1.0 ex)))
                                         (m/mevalpoly ex 1.0 3.0 1.0))))
     [1.0 1] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (* 0.5 (m/sq (max 0.0 (- 1.0 ex)))
                                         (m/mevalpoly ex 2.0 1.0))))
     [2.0 0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (* (m/fpow (max 0.0 (- 1.0 ex)) 5)
                                         (m/mevalpoly ex 1.0 5.0 9.0 5.0 1.0))))
     [2.0 1] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (* 0.25 (m/fpow (max 0.0 (- 1.0 ex)) 4)
                                         (m/mevalpoly ex 4.0 16.0 12.0 3.0))))
     [2.0 2] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (* 0.125 (m/cb (max 0.0 (- 1.0 ex)))
                                         (m/mevalpoly ex 8.0 9.0 3.0))))
     [3.0 0] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (/ (* (m/fpow (max 0.0 (- 1.0 ex)) 7)
                                            (m/mevalpoly ex 5.0 35.0 101.0 147.0 101.0 35.0 5.0)) 5.0)))
     [3.0 1] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (/ (* (m/fpow (max 0.0 (- 1.0 ex)) 6)
                                            (m/mevalpoly ex 6.0 36.0 82.0 72.0 30.0 5.0)) 6.0)))
     [3.0 2] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (/ (* (m/fpow (max 0.0 (- 1.0 ex)) 5)
                                            (m/mevalpoly ex 8.0 40.0 48.0 25.0 5.0)) 8.0)))
     [3.0 3] (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                      (/ (* (m/fpow (max 0.0 (- 1.0 ex)) 4)
                                            (m/mevalpoly ex 16.0 29.0 20.0 5.0)) 16.0)))
     (let [integr (wu-integration l)
           phi (wu-differentiation integr k)
           fact (/ ^double (phi 0.0))]
       (fn ^double [^double x] (let [ex (m/abs (/ x scale))]
                                (if (<= ex 1.0) (* fact ^double (phi ex)) 0.0)))))))


(m/unuse-primitive-operators)
