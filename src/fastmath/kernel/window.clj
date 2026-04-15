(ns fastmath.kernel.window
  "Collection of window (tapering) functions."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.quadrature :as qint]            
            [fastmath.polynomials :as poly]
            [fastmath.special :as special]))

;; https://www.researchgate.net/profile/Armin_Doerry/publication/316281181_Catalog_of_Window_Taper_Functions_for_Sidelobe_Control/links/58f92cb2a6fdccb121c9d54d/Catalog-of-Window-Taper-Functions-for-Sidelobe-Control.pdf

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmacro ^:private clamp [x body] `(if (m/<= -0.5 ~x 0.5) (m/max 0.0 (double ~body)) 0.0))
(defmacro ^:private clamp- [x body] `(if (m/<= -0.5 ~x 0.5) (double ~body) 0.0))

(defn idft
  "Create window from a spectrum (frequency domain). Calculates only real (cosine) part from IDFT definition (O(n^2))."
  [Wn]
  (let [N (count Wn)
        hN (m/* 0.5 (m/dec N))]
    (-> ((if (m/> N 200) pmap map) (fn [^long n]
                                     (let [f (m/* m/TWO_PI (m// (m/- n hN) N))]
                                       (->> Wn
                                            (map-indexed (fn [^long k ^double wn]
                                                           (m/* wn (m/cos (m/* k f)))))
                                            (v/sum)))) (range N))
        (v/div N)
        (vec))))

(defn normalize-coefficients
  "Normalize window coefficients.

  `normalize?` can be:

  * `:LInf` or `true` (default) - sets maximum value to a `1.0`
  * `:L1` - sum of the coefficients is set to a `1.0`
  * `:L2` - length of the coefficients vector is set to a `1.0`
  * `:N` - sum of the coefficients is set to a `N`
  * any number - divide all coefficients by this number"
  [coeffs normalize?]
  (condp = normalize?
    :LInf (v/normalize-LInf coeffs)
    :L1 (v/normalize-L1 coeffs)
    :L2 (v/normalize coeffs)
    :N (v/mult (v/normalize-L1 coeffs) (count coeffs))
    (cond
      (number? normalize?) (v/div coeffs normalize?)
      normalize? (v/normalize-LInf coeffs)
      :else coeffs)))

(defn sample-window
  "Sample continuous window function, returns N values from -0.5 to 0.5.

  `normalize?` can be:

  * `true` - sets mid value to a `1.0`
  * `:Linf` - sets maximum value of coefficients to a `1.0`
  * `:L1` - sum of the coefficients is set to a `1.0`
  * `:L2` - length of the coefficients vector is set to a `1.0`
  * `:N` - sum of the coefficients is set to a `N`
  * any number - divide all coefficients by this number"
  ([f ^long N] (sample-window f N true))
  ([f ^long N normalize?]
   (let [coeffs (m/sample f -0.5 0.5 N)]
     (if (true? normalize?)
       (normalize-coefficients coeffs (f 0.0))
       (normalize-coefficients coeffs normalize?)))))

;;

(defn rectangular05-continuous
  "Rectangular window with ends set to 0.5, continuous function."
  ^double [^double x]
  (let [ax (m/abs x)]
    (cond
      (m/< ax 0.5) 1.0
      (m/== ax 0.5) 0.5
      :else 0.0)))

(defn rectangular05
  "Rectangular window, both ends are set to 0.5."
  ([] (rectangular05 256))
  ([^long N] (rectangular05 N nil))
  ([^long N {:keys [normalize?] :or {normalize? false}}]
   (sample-window rectangular05-continuous N normalize?)))

;;

(defn rectangular-continuous
  "Rectangular window, continuous function."
  ^double [^double x]
  (clamp x 1.0))

(defn rectangular
  "Rectangular window."
  ([] (rectangular 256))
  ([^long N] (rectangular N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window rectangular-continuous N normalize?)))

;; 

(defn triangular-continuous
  "Triangular window, continuous functions."
  ^double [^double x]
  (clamp x (m/* 2.0 (m/- 1.0 (m/* 2.0 (m/abs x))))))

(defn ->triangular-discrete
  "Returns triangular discrete window."
  [^long shift]
  (fn [^long N]
    (let [N- (dec N)
          L (m/+ N- shift)
          hN- (m// N- 2.0)
          hL (m// L 2)]
      (map (fn [^long n] (m/- 1.0 (m/abs (m// (m/- n hN-) hL)))) (range N)))))

(defn triangular
  "Triangular window.

   `shift` parameter controls scaling factor, example values:

  * `0` - triangle, Bartlett, with zeros at both ends
  * `2` (default) same as in Python Scipy."
  ([] (triangular 256))
  ([^long N] (triangular N nil))
  ([^long N {:keys [normalize? ^long shift] :or {normalize? true shift 2}}]
   (if (m/pos? shift)
     ((->triangular-discrete shift) N)
     (sample-window triangular-continuous N normalize?))))

;;

(defn parzen-continuous
  "Parzen window, continuous function."
  ^double [^double x]
  (clamp x
    (let [ax (m/abs x)
          ax2 (m/* ax ax)]
      (if (m/<= ax 0.25)
        (m/* 8.0 (m/+ (m/- m/THIRD (m/* 8.0 ax2))
                      (m/* 16.0 ax2 ax)))
        (m/* 8.0 (m/+ (m/- m/TWO_THIRDS (m/* 4.0 ax) (m/* 5.333333333333333 ax2 ax))
                      (m/* 8.0 ax2)))))))

(defn parzen-discrete
  "Parzen window, discrete definition."
  [^long N]
  (let [N- (m/dec N)
        hN- (m// N- 2.0)
        hL (m// N 2.0)
        qL (m// N 4.0)]
    (map (fn [^long n] (let [an (m/abs (m/- n hN-))
                            nn (m// an hL)]
                        (if (m/<= an qL)
                          (m/- 1.0 (m/* 6.0 (m/sq nn) (m/- 1.0 nn)))
                          (m/* 2.0 (m/cb (m/- 1.0 nn)))))) (range N))))

(defn parzen
  "Parzen window.

  There are two ways to define the window, as a continuous function (with 0s at both ends) and as a discrete function (same as in Python's Scipy).

  By default, `:discrete?` option is set to `true`."
  ([] (parzen 256))
  ([^long N] (parzen N nil))
  ([^long N {:keys [normalize? discrete?] :or {normalize? true discrete? true}}]
   (if discrete?
     (parzen-discrete N)
     (sample-window parzen-continuous N normalize?))))


;; (b-spline 1) - is not the same as rectangular (it has 0.0 at -0.5)
;; use for order > 1

(defn ->b-spline
  [^long order]
  (let [mm (m/pow order order)
        m- (m/dec order)
        hm (m// order 2.0)]
    (fn ^double [^double x]
      (clamp x
        (loop [p (long 0)
               s 1.0
               sum 0.0]
          (if (m/> p order)
            (m/* mm sum)
            (recur (m/inc p)
                   (m/* -1.0 s)
                   (m/+ sum (m/* s (m/tpow x m- (m// (m/- p hm) order))
                                 (m/inv-factorial p)
                                 (m/inv-factorial (m/- order p)))))))))))

(defn b-spline
  "General b-spline window.
  
  `:order` - b-spline order

  * `1` - rectangular, however the first value is `0.0`
  * `2` - triangular
  * `3` (default)
  * `4` - Parzen"
  ([] (b-spline 256))
  ([^long N] (b-spline N nil))
  ([^long N {:keys [normalize? ^long order] :or {normalize? true order 3}}]
   (sample-window (->b-spline order) N normalize?)))

;;

(defn welch-continuous
  "Welch window, continuous function."
  ^double [^double x]
  (clamp x (m/* 1.5 (m/- 1.0 (m/* 4.0 x x)))))

(defn welch
  "Welch window."
  ([] (welch 256))
  ([^long N] (welch N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window welch-continuous N normalize?)))

;;

(defn ->connes
  [^double alpha]
  (let [a2 (m/* alpha alpha)
        a4 (m/* a2 a2)
        sa4 (m/* 15.0 a4)
        A (m// sa4 (m/+ 3.0 (m/* -10.0 a2) sa4))]
    (fn ^double [^double x]
      (clamp x (m/* A (m// (m/sq (m/- a2 (m/* 4.0 x x))) a4))))))

(defn connes
  "Connes window."
  ([] (connes 256))
  ([^long N] (connes N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->connes alpha) N normalize?)))

;;

(defn ->parzen-algebraic
  [^double gamma ^double u]
  (let [A (m// (m/- 1.0 (m// gamma (m/inc u))))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/- 1.0 (m/* gamma (m/pow (m/abs (m/* 2.0 x)) u))))))))

(defn parzen-algebraic
  "Parzen algebraic family window.

  Two parameters: `:gamma` (default: 1.0), 0<gamma<=1.0, and `:u` (default: 3.0), u>0."
  ([] (parzen-algebraic 256))
  ([^long N] (parzen-algebraic N nil))
  ([^long N {:keys [normalize? ^double gamma ^double u] :or {normalize? true gamma 1.0 u 3.0}}]
   (sample-window (->parzen-algebraic gamma u) N normalize?)))

;;

;; https://www.eng.buffalo.edu/Research/code/jrnl/Window.pdf

(defn ->singla-singh
  [^long order]
  (let [m21 (m/inc (m/* 2.0 order))
        K (m/* (if (m/odd? order) -1.0 1.0)
               (m/factorial m21)
               (m/sq (m/inv-factorial order)))]
    (fn ^double [^double x]
      (clamp x
        (loop [n (long 0)
               s 1.0
               sum 0.0]
          (if (m/> n order)
            (m/* 2.0 (m/- 1.0 (m/* K sum)))
            (let [d (m/- m21 n)]
              (recur (m/inc n)
                     (m/* -1.0 s)
                     (m/+ sum (m// (m/* s (m/combinations order n) (m/pow (m/abs (m/* 2.0 x)) d))
                                   d))))))))))

(defn singla-singh
  "Singla and Singh family of windows.

  Parameter `:order`, default `1`."
  ([] (singla-singh 256))
  ([^long N] (singla-singh N nil))
  ([^long N {:keys [normalize? ^double order] :or {normalize? true order 1.0}}]
   (sample-window (->singla-singh order) N normalize?)))


;;

(defn- sinc-factor ^double [^double L]
  (m// 1.0 (double (qint/gk-quadrature (fn [^double x] (m/pow (m/sinc (m/* 2.0 x)) L)) -0.5 0.5))))

(defn sinc-continuous  ^double [^double x] (clamp x (m/* 1.6963819856764424 (m/sinc (m/* 2.0 x)))))
(defn fejer-continuous  ^double [^double x] (clamp x (m/* 2.2152728287035983 (m/sq (m/sinc (m/* 2.0 x))))))
(defn de-la-vallee-poussin-continuous  ^double [^double x] (clamp x (m/* 3.008860051924727 (m/fpow (m/sinc (m/* 2.0 x)) 4))))

(defn sinc
  "Sinc lobe"
  ([] (sinc 256))
  ([^long N] (sinc N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window sinc-continuous N normalize?)))

(defn fejer
  "Fejer window, square of sinc"
  ([] (fejer 256))
  ([^long N] (fejer N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window fejer-continuous N normalize?)))

(defn de-la-vallee-poussin
  "De-La-Vallee-Poussin window, fourth power of sinc"
  ([] (de-la-vallee-poussin 256))
  ([^long N] (de-la-vallee-poussin N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window de-la-vallee-poussin-continuous N normalize?)))

(defn ->lanczos
  [^double L]
  (let [A (sinc-factor L)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/pow (m/sinc (m/* 2.0 x)) L))))))

(defn lanczos
  "Lanczos family, power of sinc"
  ([] (lanczos 256))
  ([^long N] (lanczos N nil))
  ([^long N {:keys [normalize? ^double L] :or {normalize? true L 3.0}}]
   (sample-window (->lanczos L) N normalize?)))

;; 

(defn hamming-continuous
  ^double [^double x]
  (clamp x (m/inc (m/* 0.8518518518518518 (m/cos (m/* m/TWO_PI x))))))

(defn hamming
  "Hamming window, raised cosine with alpha = 0.54"
  ([] (hamming 256))
  ([^long N] (hamming N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window hamming-continuous N normalize?)))

(defn hamming-exact-continuous
  ^double [^double x]
  (clamp x (m/inc (m/* 0.84 (m/cos (m/* m/TWO_PI x))))))

(defn hamming-exact
  "Hamming window, raised cosine with alpha = 25/46"
  ([] (hamming-exact 256))
  ([^long N] (hamming-exact N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window hamming-exact-continuous N normalize?)))

;;

(defn hann-continuous
  ^double [^double x]
  (clamp x (m/inc (m/cos (m/* m/TWO_PI x)))))

(defn hann
  "Hann window, raised cosine with alpha = 0.5"
  ([] (hann 256))
  ([^long N] (hann N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window hann-continuous N normalize?)))

;;

(defn ->raised-cosine
  [^double alpha]
  (let [f (m// (m/- 1.0 alpha) alpha)]
    (fn ^double [^double x]
      (clamp x (m/inc (m/* f (m/cos (m/* m/TWO_PI x))))))))

(defn raised-cosine
  "Raised cosine family with parameter `:alpha`, 0.5<=alpha<=1, default: `0.5`"
  ([] (raised-cosine 256))
  ([^long N] (raised-cosine N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 0.5}}]
   (sample-window (->raised-cosine alpha) N normalize?)))

;; 

(defn ->webster-hamming
  [^double v]
  (let [v+2 (m/+ v 2.0)
        alpha (m// (poly/mevalpoly v 2.0 3.0 1.0)
                   (poly/mevalpoly v 23.0 9.0 1.0))
        alpha' (m/- 1.0 alpha)
        A (m// (m/* 2.0 m/SQRTPI (m/* 0.5 (m/+ 4.0 v)))
               (m/* (m/inc (m/+ alpha v)) (special/gamma (m/* 0.5 (m/inc v)))))]
    (fn ^double [^double x]
      (let [cpx (m/cos (m/* m/PI x))]
        (clamp x (m/* A (m/+ (m/* alpha (m/pow cpx v))
                             (m/* alpha' (m/pow cpx v+2)))))))))

(defn webster-hamming
  "Generalized Hamming window, parameter `:v`>=-1/2, default: 1.0"
  ([] (webster-hamming 256))
  ([^long N] (webster-hamming N nil))
  ([^long N {:keys [normalize? ^double v] :or {normalize? true v 1.0}}]
   (sample-window (->webster-hamming v) N normalize?)))

;;

(defn ->power-of-cosine
  [^double m]
  (let [A (m// (m/* m/SQRTPI (special/gamma (m/* 0.5 (m/+ m 2.0))))
               (special/gamma (m/* 0.5 (m/inc m))))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/pow (m/cos (m/* m/PI x)) m))))))

(defn power-of-cosine
  "Power of cosine, parameter `:m`, a power, default: 1.0."
  ([] (power-of-cosine 256))
  ([^long N] (power-of-cosine N nil))
  ([^long N {:keys [normalize? ^double m] :or {normalize? true m 1.0}}]
   (sample-window (->power-of-cosine m) N normalize?)))

;;

(defn ->raised-power-of-cosine
  [^double m ^double alpha]
  (let [alpha' (m/- 1.0 alpha)
        A (m// (m/+ alpha (m// (m/* alpha' (special/gamma (m/* 0.5 (m/inc m))))
                               (m/* m/SQRTPI (special/gamma (m/* 0.5 (m/+ m 2.0)))))))]
    (fn ^double [^double x]
      (clamp x (m/* A (m/+ alpha (m/* alpha' (m/pow (m/cos (m/* m/PI x)) m))))))))

(defn raised-power-of-cosine
  "Power of cosine, parameters `:alpha` (default: 0.05) and `:m` (default: 1.0)."
  ([] (raised-power-of-cosine 256))
  ([^long N] (raised-power-of-cosine N nil))
  ([^long N {:keys [normalize? ^double alpha ^double m] :or {normalize? true alpha 0.05 m 1.0}}]
   (sample-window (->raised-power-of-cosine m alpha) N normalize?)))

;;

(defn- parzen-cosine-helper
  ^double [^double m ^double gamma ^double x]
  (m/inc (m/cos (m/* m/PI gamma (m/pow (m/abs (m/* 2.0 x)) m)))))

(defn ->parzen-cosine
  [^double m ^double gamma]
  (let [A (m// (double (qint/gk-quadrature (partial parzen-cosine-helper m gamma) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-cosine-helper m gamma x))))))

(defn parzen-cosine
  "Parzen Cosine family. Parameters `:gamma` and `:m`."
  ([] (parzen-cosine 256))
  ([^long N] (parzen-cosine N nil))
  ([^long N {:keys [normalize? ^double gamma ^double m] :or {normalize? true gamma 1.0 m 2.0}}]
   (sample-window (->parzen-cosine m gamma) N normalize?)))

;;

(defn bohman-continuous
  ^double [^double x]
  (let [x2 (m/* 2.0 (m/abs x))
        px2 (m/* m/PI x2)]
    (clamp x (m/+ (m/* 2.4674011002723395 (m/- 1.0 x2) (m/cos px2))
                  (m/* m/QUARTER_PI (m/sin px2))))))

(defn bohman
  "Bohman window, cosine lobe convolved with itself."
  ([] (bohman 256))
  ([^long N] (bohman N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window bohman-continuous N normalize?)))

;;


(defn ->trapezoid
  [^double alpha]
  (let [alpha2 (m/* 2.0 alpha)
        salpha2 (m/- 1.0 (m/sq alpha2))]
    (fn ^double [^double x]
      (let [ax (m/abs x)]
        (cond
          (m/<= ax alpha) (m// 2.0 (m/inc alpha2))
          (m/<= ax 0.5) (m// (m/* 2.0 (m/- 1.0 (m/* 2.0 ax)))
                             salpha2)
          :else 0.0)))))

(defn trapezoid
  "Trapezoid window."
  ([] (trapezoid 256))
  ([^long N] (trapezoid N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 0.25}}]
   (sample-window (->trapezoid alpha) N normalize?)))

;;

(defn ->tukey
  [^double alpha]
  (if (m/zero? alpha)
    rectangular-continuous
    (let [a' (m/* 0.5 (m/- 1.0 alpha))
          ra' (m// (m/- 0.5 a'))
          ca' (m/* m/PI ra')
          ra'2 (m/* 2.0 ra')]
      (fn ^double [^double x]
        (let [ax (m/abs x)]
          (cond
            (m/<= ax a') ra'2 
            (m/<= ax 0.5) (m/* ra' (m/inc (m/cos (m/* ca' (m/- ax a')))) )
            :else 0.0))))))

(defn tukey
  "Tukey window."
  ([] (tukey 256))
  ([^long N] (tukey N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 0.5}}]
   (sample-window (->tukey alpha) N normalize?)))

;;

(defn bartlett-hann-continuous
  ^double [^double x]
  (clamp x (m/* 2.0 (m/+ 0.62 (m/* -0.48 (m/abs x))
                         (m/* 0.38 (m/cos (m/* m/TWO_PI x)))))))

(defn bartlett-hann
  "Bartlett-Hann window."
  ([] (bartlett-hann 256))
  ([^long N] (bartlett-hann N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window bartlett-hann-continuous N normalize?)))

;;

(defn blackman-continuous
  ^double [^double x]
  (let [px (m/* m/TWO_PI x)]
    (clamp x (m/inc (m/+ (m/* 1.1904761904761905 (m/cos px))
                         (m/* 0.1904761904761905 (m/cos (m/* 2.0 px))))))))


(defn blackman
  "Blackman window."
  ([] (blackman 256))
  ([^long N] (blackman N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-continuous N normalize?)))


(defn blackman-exact-continuous
  ^double [^double x]
  (let [px (m/* m/TWO_PI x)]
    (clamp x (m/inc (m/+ (m/* 1.164021164021164 (m/cos px))
                         (m/* 0.18014613252708492 (m/cos (m/* 2.0 px))))))))

(defn blackman-exact
  "Blackman window. Exact coeffients."
  ([] (blackman-exact 256))
  ([^long N] (blackman-exact N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-exact-continuous N normalize?)))


(defn ->blackman-harris
  [coeffs]
  (let [a0 (double (first coeffs))
        A (m// a0)]
    (fn ^double [^double x]
      (let [v (m/* m/TWO_PI x)]
        (loop [id (long 1)
               cs (rest coeffs)
               sum a0]
          (if-not (seq cs)
            (clamp x (m/* A sum))
            (recur (m/inc id)
                   (rest cs)
                   (m/+ sum (m/* (double (first cs)) (m/cos (m/* id v)))))))))))

(defn blackman-harris-family
  "Sum of cosines terms."
  ([] (blackman-harris-family 256))
  ([^long N] (blackman-harris-family N nil))
  ([^long N {:keys [normalize? coeffs] :or {normalize? true coeffs [0.5 0.5]}}]
   (sample-window (->blackman-harris coeffs) N normalize?)))

(defmacro make-blackman-harris
  [clamp-fn v & coeffs]
  (let [a0 (double (first coeffs))
        A (m// a0)]
    `(~clamp-fn ~v
      (let [~'xx (m/* m/TWO_PI ~v)]
        (m/* ~A (m/+ ~a0 ~@(map-indexed (fn [^long id ^double c]
                                          (if (m/zero? id)
                                            `(m/* ~c (m/cos ~'xx))
                                            `(m/* ~c (m/cos (m/* ~(m/inc id) ~'xx))))) (rest coeffs))))))))

(defn blackman-harris-continuous ^double [^double x] (make-blackman-harris clamp x 0.4243801 0.4973406 0.0782793))
(defn blackman-harris-61db-continuous ^double [^double x] (make-blackman-harris clamp x 0.44959 0.49364 0.05677))
(defn blackman-harris-67db-continuous ^double [^double x] (make-blackman-harris clamp x 0.42323 0.49755 0.07922))
(defn blackman-harris-74db-continuous ^double [^double x] (make-blackman-harris clamp x 0.40217 0.49703 0.09892 0.00188))
(defn blackman-harris-92db-continuous ^double [^double x] (make-blackman-harris clamp x 0.35875 0.48829 0.14128 0.01168))
(defn nutall-3-1st-continuous ^double [^double x] (make-blackman-harris clamp x 0.40897 0.5 0.09103))
(defn nutall-3-3rd-continuous ^double [^double x] (make-blackman-harris clamp x 0.375 0.5 0.125))
(defn blackman-nutall-continuous ^double [^double x] (make-blackman-harris clamp x 0.3635819 0.4891775 0.1365995 0.0106411))
(defn nutall-1st-continuous ^double [^double x] (make-blackman-harris clamp x 0.355768 0.487396 0.144232 0.012604))
(defn nutall-3rd-continuous ^double [^double x] (make-blackman-harris clamp x 0.338946 0.481973 0.161054 0.018027))
(defn nutall-5th-continuous ^double [^double x] (make-blackman-harris clamp x 0.3125 0.46875 0.1875 0.03125))

(defn blackman-harris
  "Blackman-Harris window."
  ([] (blackman-harris 256))
  ([^long N] (blackman-harris N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-harris-continuous N normalize?)))

(defn blackman-harris-61db
  "Blackman-Harris window, -61dB"
  ([] (blackman-harris-61db 256))
  ([^long N] (blackman-harris-61db N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-harris-61db-continuous N normalize?)))

(defn blackman-harris-67db
  "Blackman-Harris window, -67dB"
  ([] (blackman-harris-67db 256))
  ([^long N] (blackman-harris-67db N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-harris-67db-continuous N normalize?)))

(defn blackman-harris-74db
  "Blackman-Harris window, -74dB"
  ([] (blackman-harris-74db 256))
  ([^long N] (blackman-harris-74db N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-harris-74db-continuous N normalize?)))

(defn blackman-harris-92db
  "Blackman-Harris window, -92dB"
  ([] (blackman-harris-92db 256))
  ([^long N] (blackman-harris-92db N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-harris-92db-continuous N normalize?)))

(defn nutall-3-1st
  "Nutall three-term, 1st derivative continuous"
  ([] (nutall-3-1st 256))
  ([^long N] (nutall-3-1st N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window nutall-3-1st-continuous N normalize?)))

(defn nutall-3-3rd
  "Nutall three-term, 3rd derivative continuous"
  ([] (nutall-3-3rd 256))
  ([^long N] (nutall-3-3rd N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window nutall-3-3rd-continuous N normalize?)))

(defn blackman-nutall
  "Blackman-Nutall window"
  ([] (blackman-nutall 256))
  ([^long N] (blackman-nutall N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window blackman-nutall-continuous N normalize?)))

(defn nutall-1st
  "Nutall four-term, 1st derivative continuous"
  ([] (nutall-1st 256))
  ([^long N] (nutall-1st N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window nutall-1st-continuous N normalize?)))

(defn nutall-3rd
  "Nutall four-term, 3rd derivative continuous"
  ([] (nutall-3rd 256))
  ([^long N] (nutall-3rd N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window nutall-3rd-continuous N normalize?)))

(defn nutall-5th
  "Nutall four-term, 5th derivative continuous"
  ([] (nutall-5th 256))
  ([^long N] (nutall-5th N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window nutall-5th-continuous N normalize?)))

(defn mottaghi-kashtiban-shayesteh
  "Four-term Blackman-Harris window, with coefficients depending on N"
  ([] (mottaghi-kashtiban-shayesteh 256))
  ([^long N] (mottaghi-kashtiban-shayesteh N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (let [a0 (m/- 0.5363 (m// 0.14 (m/dec N)))]
     (sample-window (->blackman-harris [a0 (m/- 0.996 a0) 0.0 0.004]) N normalize?))))

;;

(defn low-sidelobe-continuous ^double [^double x] (make-blackman-harris clamp x 0.471492057 0.0 0.17553428 0.028497078 0.001261367))

(defn low-sidelobe
  "Low-sidelobe"
  ([] (low-sidelobe 256))
  ([^long N] (low-sidelobe N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window low-sidelobe N normalize?)))

;; 

(defn ->exponential
  [^double alpha]
  (let [A (m// alpha (m/- 1.0 (m/exp (m/- alpha))))
        f (m/* -2.0 alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f (m/abs x))))))))

(defn exponential
  "Exponential family, parameter `:alpha`, decay, default: 1.0"
  ([] (exponential 256))
  ([^long N] (exponential N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->exponential alpha) N normalize?)))

;;

(defn ->hanning-poisson
  [^double alpha]
  (let [a2 (m/* alpha alpha)
        A (m// (m/* alpha (m/+ a2 m/PI2))
               (m/- (m/+ (m/* 2.0 a2) m/PI2)
                    (m/* m/PI2 (m/exp (m/- alpha)))))
        f (m/* -2.0 alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f (m/abs x))) (m/inc (m/cos (m/* m/TWO_PI x))))))))

(defn hanning-poisson
  "Hanning Poisson (product of Exponential and Hann windows)."
  ([] (hanning-poisson 256))
  ([^long N] (hanning-poisson N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->hanning-poisson alpha) N normalize?)))

;;

(defn ->gaussian
  [^double alpha]
  (let [A (m// (m/* alpha m/SQRT_2_PI)
               (special/erf (m// alpha m/SQRT2)))
        f (m/* -2.0 alpha alpha)]
    (fn ^double [^double x]
      (clamp x (m/* A (m/exp (m/* f x x)))))))

(defn gaussian
  "Truncated gaussian, parameter `:alpha` controls spread."
  ([] (gaussian 256))
  ([^long N] (gaussian N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->gaussian alpha) N normalize?)))

;;

(defn- parzen-exponential-fn
  ^double [^double alpha ^double r ^double x]
  (m/exp (m/- (m/pow (m/abs (m/* 2.0 alpha x)) r))))

(defn ->parzen-exponential
  [^double alpha ^double r]
  (let [A (m// (double (qint/gk-quadrature (partial parzen-exponential-fn alpha r) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-exponential-fn alpha r x))))))

(defn parzen-exponential
  "Parzen exponential, parameters `:alpha` (default: 1.0) and `:r` (power, default: 1.0)"
  ([] (parzen-exponential 256))
  ([^long N] (parzen-exponential N nil))
  ([^long N {:keys [normalize? ^double alpha ^double r] :or {normalize? true alpha 1.0 r 1.0}}]
   (sample-window (->parzen-exponential alpha r) N normalize?)))

;;

(defn- eta
  "dBc to linear"
  ^double [^double level]
  (m/exp10 (m// (m/abs level) 20.0)))

(defn ->dolph-chebyshev
  [^double level]
  (let [e (eta level)
        ace (m/acosh e)]
    (fn [^long N]
      (let [N- (m/dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (double (poly/eval-chebyshev-T N- x0)))
            dN (double N)]
        (->> (range N)
             (map (fn [^long n] (m/* f (poly/eval-chebyshev-T N- (m/* x0 (m/cospi (m// n dN)))))))
             (idft))))))

(defn dolph-chebyshev
  "Dolph-Chebyshev, parameter `:level` controls side-lobe level (in dB, default: -50.0)"
  ([] (dolph-chebyshev 256))
  ([^long N] (dolph-chebyshev N nil))
  ([^long N {:keys [normalize? ^double level] :or {normalize? true level -50.0}}]
   (normalize-coefficients ((->dolph-chebyshev level) N) normalize?)))

;;

(defn- taylor-f-denominator
  ^double [^long n ^long m]
  (let [m2 (m/sq m)]
    (reduce (fn [^double p ^double n']
              (if (m/== n' m)
                p
                (m/* p (m/- 1.0 (m// m2 (m/* n' n')))))) 1.0 (range 1 n))))

(defn- taylor-f-numerator
  ^double [^double A2 ^double sigma2 ^long n ^long m]
  (let [msigma2 (m// (m/* m m) sigma2)        ]
    (reduce (fn [^double p ^double n']
              (m/* p (m/- 1.0 (m// msigma2 (m/+ A2 (m/sq (m/- n' 0.5))))))) 1.0 (range 1 n))))

(defn- taylor-f
  ^double [^double A2 ^double sigma2 ^long n ^long m]
  (let [s (if (m/even? m) -0.5 0.5)]
    (m// (m/* 2.0 s (taylor-f-numerator A2 sigma2 n m))
         (taylor-f-denominator n m))))

(defn ->taylor
  [^double level ^long n]
  (let [e (eta level)
        A (m// (m/acosh e) m/PI)
        A2 (m/* A A)
        sigma2 (m// (m/* n n)
                    (m/+ A2 (m/sq (m/- n 0.5))))
        coeffs (vec (conj (map (partial taylor-f A2 sigma2 n) (range 1 n)) 1.0))]
    (->blackman-harris coeffs)))

(defn taylor
  "Taylor, approximation of Dolph-Chebyshev window. Parameter `:level` controls side-lobe level (in dB, default: -50.0), `:n` number of constant side-lobes (from the main-lobe) (default: 4)."
  ([] (taylor 256))
  ([^long N] (taylor N nil))
  ([^long N {:keys [normalize? ^double level ^long n] :or {normalize? true level -50.0 n 4}}]
   (sample-window (->taylor level n) N normalize?)))

;;

(defn ->cauchy
  [^double alpha]
  (let [A (m// alpha (m/atan alpha))]
    (fn ^double [^double x]
      (clamp x (m/* A (m// (m/inc (m/sq (m/* 2.0 alpha x)))))))))

(defn cauchy
  "Cauchy, parameter `:alpha` (decay, default: 1.0)."
  ([] (cauchy 256))
  ([^long N] (cauchy N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->cauchy alpha) N normalize?)))

;;

(defn- parzen-geometric-fn
  ^double [^double alpha ^double r ^double x]
  (m// (m/inc (m/pow (m/abs (m/* 2.0 alpha x)) r))))

(defn ->parzen-geometric
  [^double alpha ^double r]
  (let [A (m// (double (qint/gk-quadrature (partial parzen-geometric-fn alpha r) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (parzen-geometric-fn alpha r x))))))

(defn parzen-geometric
  "Parzen-Geometric, parameter `:alpha` (default: 1.0) and `:r` (power, default: 1.0)"
  ([] (parzen-geometric 256))
  ([^long N] (parzen-geometric N nil))
  ([^long N {:keys [normalize? ^double alpha ^double r] :or {normalize? true alpha 1.0 r 1.0}}]
   (sample-window (->parzen-geometric alpha r) N normalize?)))

;;

(defn ->kaiser-bessel
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        A (m// pa (m/sinh pa))]
    (fn ^double [^double x]
      (clamp x (m/* A (special/bessel-I0 (m/* pa (m/sqrt (m/- 1.0 (m/sq (m/* 2.0 x)))))))))))

(defn kaiser-bessel
  "Kaiser-Bessel, approximation of DPSS/Slepian window, parameter `:alpha` (default: 1.0)"
  ([] (kaiser-bessel 256))
  ([^long N] (kaiser-bessel N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 1.0}}]
   (sample-window (->kaiser-bessel alpha) N normalize?)))

;;

(defn- cosh-fn
  ^double [^double alpha ^double x]
  (let [pa (m/* m/PI alpha)]
    (m// (m/cosh (m/* pa (m/sqrt (m/- 1.0 (m/* 4.0 x x)))))
         (m/cosh pa))))

(defn ->cosh
  [^double alpha]
  (let [A (m// (double (qint/gk-quadrature (partial cosh-fn alpha) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (cosh-fn alpha x))))))

(defn cosh
  "Cosh, approximation of DPSS/Slepian window, parameter `:alpha` (default: 2.0)"
  ([] (cosh 256))
  ([^long N] (cosh N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 2.0}}]
   (sample-window (->cosh alpha) N normalize?)))

;;

(defn avci-nacaroglu-fn
  ^double [^double alpha ^double x]
  (let [pa (m/* m/PI alpha)]
    (m/exp (m/* pa (m/dec (m/sqrt (m/- 1.0 (m/* 4.0 x x))))))))

(defn ->avci-nacaroglu
  [^double alpha]
  (let [A (m// (double (qint/gk-quadrature (partial avci-nacaroglu-fn alpha) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (avci-nacaroglu-fn alpha x))))))

(defn avci-nacaroglu
  "Avci-Nacaroglu, approximation of DPSS/Slepian window, parameter `:alpha` (default: 2.0)"
  ([] (avci-nacaroglu 256))
  ([^long N] (avci-nacaroglu N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 2.0}}]
   (sample-window (->avci-nacaroglu alpha) N normalize?)))

;;

(defn knab-fn
  ^double [^double alpha ^double x]
  (let [pa (m/* m/PI alpha)
        sx (m/sqrt (m/- 1.0 (m/* 4.0 x x)))]
    (if (m/zero? sx)
      0.0
      (m// (m/sinh (m/* pa sx))
           (m/* (m/sinh pa) sx)))))

(defn ->knab
  [^double alpha]
  (let [A (m// (double (qint/gk-quadrature (partial knab-fn alpha) -0.5 0.5)))]
    (fn ^double [^double x]
      (clamp x (m/* A (knab-fn alpha x))))))

(defn knab
  "Knab, approximation of DPSS/Slepian window, parameter `:alpha` (default: 2.0)"
  ([] (knab 256))
  ([^long N] (knab N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 2.0}}]
   (sample-window (->knab alpha) N normalize?)))

;;

(defn ->ultraspherical
  [^double level ^double alpha]
  (let [e (eta level)
        ace (m/acosh e)]
    (fn [^long N]
      (let [N- (m/dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (double (poly/eval-gegenbauer-C N- alpha x0)))
            dN (double N)]
        (->> (range N)
             (mapv (fn [^long n] (m/* f (poly/eval-gegenbauer-C N- alpha (m/* x0 (m/cospi (m// n dN)))))))
             (idft))))))

(defn ultraspherical
  "Ultraspherical, parameter `:level` controls side-lobe level (in dB, default: -50.0), `:alpha`"
  ([] (ultraspherical 256))
  ([^long N] (ultraspherical N nil))
  ([^long N {:keys [normalize? ^double level ^double alpha] :or {normalize? true level -50.0 alpha 2.0}}]
   (normalize-coefficients ((->ultraspherical level alpha) N) normalize?)))

(defn saramaki
  ([] (saramaki 256))
  ([^long N] (saramaki N nil))
  ([^long N {:keys [normalize? ^double level] :or {normalize? true level -50.0}}]
   (normalize-coefficients ((->ultraspherical level 1.0) N) normalize?)))

;;

(defn ->legendre
  [^double level]
  (let [ace (m/acosh (m/exp10 (m// (m/+ (m/* 1.0754 (m/abs level)) 1.7388) 20.0)))]
    (fn [^long N]
      (let [N- (m/dec N)
            x0 (m/cosh (m// ace N-))
            f (m// N (double (poly/eval-legendre-P N- x0)))
            dN (double N)]
        (->> (range N)
             (map (fn [^long n] (m/* f (poly/eval-legendre-P N- (m/* x0 (m/cospi (m// n dN)))))))
             (idft))))))

(defn legendre
  "Legendre, parameter `:level` controls side-lobe level (in dB, default: -50.0)"
  ([] (legendre 256))
  ([^long N] (legendre N nil))
  ([^long N {:keys [normalize? ^double level] :or {normalize? true level -50.0}}]
   (normalize-coefficients ((->legendre level) N) normalize?)))

;;

(defn ->bessel-I1
  [^double alpha]
  (let [pa (m/* m/PI alpha)
        i1pa (special/bessel-I1 pa)
        A (m// (m/* pa i1pa)
               (m/dec (m/cosh pa)))]
    (fn ^double [^double x]
      (let [x (if (m/== x 0.5) 0.49999999999999994
                  (if (m/== x -0.5) -0.49999999999999994
                      x))]
        (clamp x (let [xx (m/sqrt (m/- 1.0 (m/* 4.0 x x)))]
                   (m/* A (m// (special/bessel-I1 (m/* pa xx))
                               (m/* i1pa xx)))))))))

(defn bessel-I1
  "Modified first-order Bessel, parameter `:alpha` (default: 2.0)"
  ([] (bessel-I1 256))
  ([^long N] (bessel-I1 N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? true alpha 2.0}}]
   (sample-window (->bessel-I1 alpha) N normalize?)))

;;

(defn- shayesteh-kashtiban-n
  ^double [^long n ^long N-]
  (m/pow (m/sinc (m// (m/- n (m/* 0.5 N-))
                      (m/* 0.654 N-))) 2.5))


(defn shayesteh-kashtiban-discrete
  [^long N]
  (let [N- (m/dec N)
        endpoint (m/+ 0.02 (m/* 0.001 N-) (m// (m/+ 50.0 N- N-)))
        vs (mapv (fn [^long n]
                   (if (or (m/zero? n) (m/== n N-))
                     endpoint
                     (shayesteh-kashtiban-n n N-))) (range N))
        A (m// (v/sum vs))]
    (v/mult vs A)))

(defn shayesteh-kashtiban
  "Shayesteh-Kashtiban window."
  ([] (shayesteh-kashtiban 256))
  ([^long N] (shayesteh-kashtiban N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (normalize-coefficients (shayesteh-kashtiban-discrete N) normalize?)))

;;

(defn ->kaiser-bessel-derived
  [^double alpha]
  (fn [^long N]
    (let [kb (kaiser-bessel (m/inc (m/ceil (m// N 2.0))) {:alpha alpha :normalize? false})
          v (v/sqrt (v/div (butlast (reductions m/+ kb)) (v/sum kb)))]
      (vec (concat v (if (m/odd? N) (rest (reverse v)) (reverse v)))))))

(defn kaiser-bessel-derived
  "Kaiser-Bessel-derived, KBD, parameter `:alpha`"
  ([] (kaiser-bessel-derived 256))
  ([^long N] (kaiser-bessel-derived N nil))
  ([^long N {:keys [normalize? ^double alpha] :or {normalize? false alpha 1.0}}]
   (normalize-coefficients ((->kaiser-bessel-derived alpha) N) normalize?)))

;;

(defn vorbis-continuous
  ^double [^double x]
  (clamp x (m/* 1.660592492619163 (m/sin (m/* m/HALF_PI (m/sq (m/cospi x)))))))

(defn vorbis
  "Vorbis window"
  ([] (vorbis 256))
  ([^long N] (vorbis N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window vorbis-continuous N normalize?)))

;;

(defn flat-top-continuous ^double [^double x] (make-blackman-harris clamp- x 0.21557895 0.41663158 0.277263158 0.083578947 0.006947368))

(defn flat-top
  "Flat top window"
  ([] (flat-top 256))
  ([^long N] (flat-top N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window flat-top-continuous N normalize?)))

(defn flat-top-3-continuous ^double [^double x] (make-blackman-harris clamp- x 0.2811 0.5209 0.1980))

(defn flat-top-3
  "Flat top window, version with 3 coefficients"
  ([] (flat-top-3 256))
  ([^long N] (flat-top-3 N nil))
  ([^long N {:keys [normalize?] :or {normalize? true}}]
   (sample-window flat-top-3-continuous N normalize?)))

;;

(defn continuous->window-fn
  "Convert any y=f(x), x=[-0.5,0.5], continuous function to a window function (adds normalization)."
  [continuous-window]
  (fn custom-continuous-window
    ([] (custom-continuous-window 256))
    ([^long N] (custom-continuous-window N nil))
    ([^long N {:keys [normalize?] :or {normalize? true}}]
     (sample-window continuous-window N normalize?))))

(defn discrete->window-fn
  "Convert any y=f(n), n=0,...,N-1, discrete function to a window function (adds normalization)."
  [discrete-window]
  (fn custom-discrete-window
    ([] (custom-discrete-window 256))
    ([^long N] (custom-discrete-window N nil))
    ([^long N {:keys [normalize?] :or {normalize? true}}]
     (normalize-coefficients (discrete-window N) normalize?))))

