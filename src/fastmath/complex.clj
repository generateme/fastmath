;; # Namespace scope
;;
;; Functions to manipulate Vec2 as Complex numbers.
;; Implementation based on Apache Commons Math and https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/Complex.cpp

(ns fastmath.complex
  "Complex numbers functions.

  Complex number is represented as `Vec2` type (from [[clojure2d.math.vector]] namespace).

  To create complex number use [[complex]], [[vec2]] or [[->Vec2]].

  Implementation checks for ##Inf, ##NaN and some of the function distinguish +0.0 and -0.0" 
  (:refer-clojure :exclude [abs zero?])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators #{'abs 'zero?})

(def I (Vec2. 0.0 1.0))
(def I- (Vec2. 0.0 -1.0))
(def -I (Vec2. 0.0 -1.0))
(def ONE (Vec2. 1.0 0.0))
(def TWO (Vec2. 2.0 0.0))
(def ZERO (Vec2. 0.0 0.0))
(def PI (Vec2. m/PI 0.0))

(defn complex
  "Create complex number. Represented as `Vec2`."
  ([^double a ^double b] (Vec2. a b))
  ([^double a] (Vec2. a 0.0))
  ([] ZERO))

(defn ensure-complex
  "Convert possible number to complex or return input."
  [v] (if (number? v) (complex (double v)) v))

(defn re "Real part" ^double [^Vec2 z] (.x z))
(defn im "Imaginary part" ^double [^Vec2 z] (.y z))

(defn real? "Is z is a real number?"  [^Vec2 z] (m/zero? (.y z)))
(defn imaginary? "Is z is a pure imaginary number?"  [^Vec2 z] (m/zero? (.x z)))
(defn zero? "Is zero?" [^Vec2 z] (v/is-zero? z))
(defn inf? "Is infinite?" [^Vec2 z] (or (m/inf? (.x z))
                                     (m/inf? (.y z))))
(defn nan? "Is NaN?" [^Vec2 z] (or (m/nan? (.x z))
                                (m/nan? (.y z))))
(defn invalid? "Is NaN or Inf?" [z] (or (inf? z) (nan? z)))
(defn valid? "Is valid complex (not NaN or Inf)?" [z] (not (invalid? z)))

(defn delta-eq
  "Compare complex numbers with given accuracy (10e-6 by default)"
  ([q1 q2]
   (v/delta-eq q1 q2))
  ([q1 q2 ^double accuracy]
   (v/delta-eq q1 q2 accuracy)))

(defn csgn
  "Complex sgn.

  Returns `0` for `0+0i` or calls `m/sgn` on real part otherwise."
  (^double [^double re ^double im]
   (if (and (m/zero? re)
            (m/zero? im)) 0.0 (m/sgn re)))
  (^double [^Vec2 z]
   (if (zero? z) 0.0 (m/sgn (.x z)))))

(defn flip "Exchange imaginary and real parts" ^Vec2 [^Vec2 z] (Vec2. (.y z) (.x z)))

(defn abs
  "Absolute value, magnitude"
  ^double [^Vec2 z]
  (if (inf? z)
    ##Inf
    (m/sqrt (+ (m/sq (.x z)) (m/sq (.y z))))))

(defn norm
  "Norm (Guass) of the complex number, absolute value squared"
  ^double [^Vec2 z]
  (if (inf? z)
    ##Inf
    (+ (m/sq (.x z)) (m/sq (.y z)))))

(defn arg
  "Argument (angle) of the complex number"
  ^double [z]
  (v/heading z))

(defn add
  "Sum of two complex numbers"
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (Vec2. (+ (.x z1) (.x z2))
         (+ (.y z1) (.y z2))))

(defn adds
  "Add scalar to complex number."
  ^Vec2 [^Vec2 z ^double v]
  (Vec2. (+ (.x z) v) (.y z)))

(defn sub
  "Difference of two complex numbers"
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (Vec2. (- (.x z1) (.x z2))
         (- (.y z1) (.y z2))))

(defn scale
  "Multiply by real number"
  ^Vec2 [^Vec2 z ^double v]
  (Vec2. (* (.x z) v)
         (* (.y z) v)))

(defn mult
  "Multiply two complex numbers."
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (let [a (.x z1)
        b (.y z1)
        c (.x z2)
        d (.y z2)]
    (Vec2. (- (* a c) (* b d))
           (+ (* a d) (* b c)))))

(defn mult-I ^Vec2 [^Vec2 z] (Vec2. (- (.y z)) (.x z)))
(defn mult-I- ^Vec2 [^Vec2 z] (Vec2. (.y z) (- (.x z))))

(defn muladd
  "`(x y z)` -> `(+ z (* x y))`"
  ^Vec2 [x y z]
  (add z (mult x y)))

(defn sq
  "Square complex number. \\\\(z^2\\\\)"
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (Vec2. (- (* a a) (* b b))
           (* 2.0 a b))))

(defn conjugate
  "Complex conjugate. \\\\(\\bar{z}\\\\)"
  ^Vec2 [^Vec2 z] 
  (Vec2. (.x z) (- (.y z))))

(defn div
  "Divide two complex numbers."
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (let [a (.x z1)
        b (.y z1)
        c (.x z2)
        d (.y z2)
        den (+ (* c c) (* d d))]
    (Vec2. (/ (+ (* a c) (* b d)) den)
           (/ (- (* b c) (* a d)) den))))

(defn reciprocal
  "\\\\(\\frac{1}{z}\\\\)"
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (cond
      (zero? z) (Vec2. (m/copy-sign ##Inf a)
                       (- (m/copy-sign ##Inf b)))
      (inf? z) (Vec2. (m/copy-sign 0.0 a)
                      (- (m/copy-sign 0.0 b)))
      (nan? z) (Vec2. ##NaN ##NaN)
      :else (let [den (+ (* a a) (* b b))]
              (Vec2. (/ a den)
                     (/ (- b) den))))))
(defn neg
  "Negate complex number. \\\\(-z\\\\)"
  ^Vec2 [^Vec2 z]
  (Vec2. (- (.x z)) (- (.y z))))

(defn sqrt
  "Sqrt of complex number. \\\\(\\sqrt{z}\\\\)"
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    z
    (let [a (.x z)
          b (.y z)
          l (m/sqrt (+ (* a a) (* b b)))
          aa (m/sqrt (+ l a))
          bb (* (m/sgn b) (m/sqrt (- l a)))]
      (Vec2. (* m/SQRT2_2 aa) (* m/SQRT2_2 bb)))))

(defn sqrt1z
  "\\\\(\\sqrt{1-z^2}\\\\)"
  ^Vec2 [z]
  (->> z
       (mult z)
       (sub ONE)
       (sqrt)))

;; exp

(defn exp
  "exp"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/exp x) y)
      (m/zero? x) (Vec2. (m/cos y) (m/sin y))
      :else (let [e (m/exp x)]
              (Vec2. (* e (m/cos y))
                     (* e (m/sin y)))))))

;; log

(defn log
  "log, principal value"
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (cond
      (m/zero? b) (cond
                    (m/pos? a) (Vec2. (m/log a) b)
                    (m/zero? a) (Vec2. ##-Inf ##NaN)
                    :else (Vec2. (m/log (- a)) (m/copy-sign m/PI b)))
      (m/zero? a) (if (m/pos? b)
                    (Vec2. (m/log b) m/HALF_PI)
                    (Vec2. (m/log (- b)) m/-HALF_PI))
      :else (Vec2. (* 0.5 (m/log (+ (* a a) (* b b))))
                   (m/atan2 b a)))))

(defn logb
  "log with base b"
  ^Vec2 [z b]
  (div (log z) (log b)))

;; bunch of power branches

(defn- pow-a+0i-x+0i
  [^double r ^double zr]
  (if (pos? r) ;; real exponent
    (Vec2. (m/pow r zr) 0.0)
    (let [e (m/exp (* zr (m/log (- r))))
          p (* zr m/PI)]
      (Vec2. (* e (m/cos p))
             (* e (m/sin p))))))

(defn- pow-a+0i-0+yi
  [^double r ^double i ^double zi]
  (if (pos? r) ;; pure imaginary exponent
    (let [p (* zi (m/log r))]
      (Vec2. (m/cos p) (m/sin p)))
    (let [p (* zi (m/log (- r)))]
      (if (m/negative-zero? i)
        (let [e (m/exp (* zi m/PI))]
          (Vec2. (* e (m/cos p)) (* e (m/sin p))))
        (let [e (m/exp (* zi m/-PI))]
          (Vec2. (* e (m/cos p)) (* e (m/sin p))))))))

(defn- pow-a+0i-x+yi-neg-0
  [^double r ^double zr ^double zi]
  (let [pw (m/pow (- r) zr)
        e (m/exp (* zi m/PI))
        p (+ (* -1.0 zi (m/log (- r)))
             (* zr m/PI))]
    (Vec2. (* pw e (m/cos p))
           (* -1.0 pw e (m/sin p)))))

(defn- pow-a+0i-x+yi-pos-0
  [^double r ^double zr ^double zi]
  (let [lr (m/log (- r))
        e (m/exp (- (* zr lr)
                    (* zi m/PI)))
        p (+ (* zi lr)
             (* zr m/PI))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-a+0i-x+yi
  [^double r ^double i ^double zr ^double zi]
  (if (pos? r)
    (let [lr (m/log r)
          e (m/exp (* zr lr))
          p (* zi lr)]
      (Vec2. (* e (m/cos p)) (* e (m/sin p))))
    (if (m/negative-zero? i)
      (pow-a+0i-x+yi-neg-0 r zr zi)
      (pow-a+0i-x+yi-pos-0 r zr zi))))

(defn- pow-0+bi-x+0i
  [^double i ^double zr]
  (let [hzr (* 0.5 zr)
        e (m/exp (* hzr (m/log (* i i))))
        p (* hzr (m/sgn i) m/PI)]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-0+bi-0+yi
  [^double i ^double zi]
  (let [hzi (* 0.5 zi)
        e (m/exp (* -1.0 hzi (m/sgn i) m/PI))
        p (* hzi (m/log (* i i)))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-0+bi-x+yi
  [^double i ^double zr ^double zi]
  (let [hzr (* 0.5 zr)
        hzi (* 0.5 zi)
        lii (m/log (* i i))
        spi (* (m/sgn i) m/PI)
        e (m/exp (- (* hzr lii)
                    (* hzi spi)))
        p (+ (* hzi lii)
             (* hzr spi))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-a+bi-a+0i
  [^double r ^double i ^double zr]
  (let [e (m/exp (* 0.5 zr (m/log (+ (* r r) (* i i)))))
        p (* zr (m/atan2 i r))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-a+bi-0+yi
  [^double r ^double i ^double zi]
  (let [e (m/exp (* -1.0 zi (m/atan2 i r)))
        p (* 0.5 zi (m/log (+ (* r r) (* i i))))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn- pow-a+bi-x+yi
  [^double r ^double i ^double zr ^double zi]
  (let [lrrii (m/log (+ (* r r) (* i i)))
        at2 (m/atan2 i r)
        e (m/exp (- (* 0.5 zr lrrii)
                    (* zi at2)))
        p (+ (* 0.5 zi lrrii)
             (* zr at2))]
    (Vec2. (* e (m/cos p)) (* e (m/sin p)))))

(defn pow
  "Power. \\\\(z_1^{z_2}\\\\)"
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (if (zero? z1)
    (if (zero? z2) ;; both zero
      (Vec2. ##NaN ##NaN)
      (if (m/negative-zero? (.x z1)) ;; zero base, non-zero exponent 
        (Vec2. ##Inf ##Inf)
        ZERO))
    (if (zero? z2) ;; non-zero base and zero exponent, always one
      ONE
      (let [r (.x z1)
            i (.y z1)
            zr (.x z2)
            zi (.y z2)]
        (cond (m/zero? i) (cond
                            (m/zero? zi) (pow-a+0i-x+0i r zr)
                            (m/zero? zr) (pow-a+0i-0+yi r i zi)
                            :else (pow-a+0i-x+yi r i zr zi))
              (m/zero? r) (cond
                            (m/zero? zi) (pow-0+bi-x+0i i zr)
                            (m/zero? zr) (pow-0+bi-0+yi i zi)
                            :else (pow-0+bi-x+yi i zr zi))
              :else (cond
                      (m/zero? zi) (pow-a+bi-a+0i r i zr)
                      (m/zero? zr) (pow-a+bi-0+yi r i zi)
                      :else (pow-a+bi-x+yi r i zr zi)))))))

;;;;;; trig

(defn sin
  "sin"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sin x) (* y (m/cos x)))
      (m/zero? x) (Vec2. x (m/sinh y))
      :else (Vec2. (* (m/sin x) (m/cosh y))
                   (* (m/cos x) (m/sinh y))))))

(defn cos
  "cos"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/cos x) (* -1.0 y (m/sin x)))
      (m/zero? x) (Vec2. (m/cosh y) (* -1.0 y x))
      :else (Vec2. (* (m/cos x) (m/cosh y))
                   (* -1.0 (m/sin x) (m/sinh y))))))


(defn tan
  "tan"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/tan x) y)
      (m/zero? x) (Vec2. x (m/tanh y))
      :else (let [cosx (m/cos x)
                  sinhy (m/sinh y)
                  denom (+ (* cosx cosx) (* sinhy sinhy))]
              (Vec2. (/ (* (m/sin x) cosx) denom)
                     (/ (* sinhy (m/cosh y)) denom))))))

(defn sec
  "sec"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sec x) (- y))
      (m/zero? x) (Vec2. (m/sech y) (* x (m/sgn y)))
      :else (let [cosx (m/cos x)
                  sinhy (m/sinh y)
                  denom (+ (* cosx cosx) (* sinhy sinhy))]
              (Vec2. (/ (* cosx (m/cosh y)) denom)
                     (/ (* (m/sin x) sinhy) denom))))))

(defn csc
  "csc"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/csc x) (- y))
      (m/zero? x) (Vec2. x (- (m/csch y)))
      :else (let [sinx (m/sin x)
                  sinhy (m/sinh y)
                  denom (+ (* sinx sinx) (* sinhy sinhy))]
              (Vec2. (/ (* sinx (m/cosh y)) denom)
                     (/ (* (- (m/cos x)) sinhy) denom))))))

(defn cot
  "csc"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/cot x) y)
      (m/zero? x) (Vec2. x (- (m/coth y)))
      :else (let [sinx (m/sin x)
                  sinhy (m/sinh y)
                  denom (+ (* sinx sinx) (* sinhy sinhy))]
              (Vec2. (/ (* sinx (m/cos x)) denom)
                     (/ (* -1.0 sinhy (m/cosh y)) denom))))))

;;

(defn sinh
  "sinh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sinh x) y)
      (m/zero? x) (Vec2. x (m/sin y))
      :else (Vec2. (* (m/sinh x) (m/cos y))
                   (* (m/cosh x) (m/sin y))))))

(defn cosh
  "cosh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/cosh x) y)
      (m/zero? x) (Vec2. (m/cos y) (- x))
      :else (Vec2. (* (m/cosh x) (m/cos y))
                   (* (m/sinh x) (m/sin y))))))

(defn tanh
  "tanh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/tanh x) y)
      (m/zero? x) (Vec2. x (m/tan y))
      :else (let [sinhx (m/sinh x)
                  cosy (m/cos y)
                  denom (+ (* sinhx sinhx) (* cosy cosy))]
              (Vec2. (/ (* sinhx (m/cosh x)) denom)
                     (/ (* (m/sin y) cosy) denom))))))

(defn sech
  "sech"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sech x) (- y))
      (m/zero? x) (Vec2. (m/sec y) x)
      :else (let [sinhx (m/sinh x)
                  cosy (m/cos y)
                  denom (+ (* sinhx sinhx) (* cosy cosy))]
              (Vec2. (/ (* (m/cosh x) cosy) denom)
                     (/ (* -1.0 sinhx (m/sin y)) denom))))))

(defn csch
  "csch"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/csch x) (- y))
      (m/zero? x) (Vec2. x (- (m/csc y)))
      :else (let [sinhx (m/sinh x)
                  siny (m/sin y)
                  denom (+ (* sinhx sinhx) (* siny siny))]
              (Vec2. (/ (* sinhx (m/cos y)) denom)
                     (/ (* -1.0 (m/cosh x) siny) denom))))))

(defn coth
  "coth"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/coth x) y)
      (m/zero? x) (Vec2. x (- (m/cot y)))
      :else (let [sinhx (m/sinh x)
                  siny (m/sin y)
                  denom (+ (* sinhx sinhx) (* siny siny))]
              (Vec2. (/ (* sinhx (m/cosh x)) denom)
                     (/ (* -1.0 (m/cos y) siny) denom))))))

;;

(defn asin
  "asin"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (cond
                    (m/between? -1.0 1.0 x) (Vec2. (m/asin x) y)
                    (> x 1.0) (Vec2. m/HALF_PI (* (m/copy-sign 1.0 y)
                                                  (m/log (+ (m/sqrt (dec (* x x))) x))))
                    :else (Vec2. m/-HALF_PI (* (m/copy-sign 1.0 y)
                                               (m/log (- (m/sqrt (dec (* x x))) x)))))
      (m/zero? x) (Vec2. x (m/log (+ y (m/sqrt (inc (* y y))))))
      :else (let [ss (inc (+ (* x x) (* y y)))
                  r2 (* 2.0 x)
                  ssp2r (m/sqrt (+ ss r2))
                  ssm2r (m/sqrt (- ss r2))
                  sum (* 0.5 (+ ssp2r ssm2r))]
              (Vec2. (m/asin (* 0.5 (- ssp2r ssm2r)))
                     (* (m/sgn y) (m/log (+ sum (m/sqrt (dec (* sum sum)))))))))))

(defn acos
  "acos"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (if (m/zero? y)
      (cond
        (m/between? -1.0 1.0 x) (Vec2. (m/acos x) (- y))
        (> x 1.0) (Vec2. 0.0 (* -1.0 (m/copy-sign 1.0 y)
                                (m/log (+ (m/sqrt (dec (* x x))) x))))
        :else (Vec2. m/PI (* -1.0 (m/copy-sign 1.0 y)
                             (m/log (- (m/sqrt (dec (* x x))) x)))))
      (let [ss (inc (+ (* x x) (* y y)))
            r2 (* 2.0 x)
            ssp2r (m/sqrt (+ ss r2))
            ssm2r (m/sqrt (- ss r2))
            sum (* 0.5 (+ ssp2r ssm2r))]
        (Vec2. (m/acos (* 0.5 (- ssp2r ssm2r)))
               (* -1.0 (m/sgn y) (m/log (+ sum (m/sqrt (dec (* sum sum)))))))))))

(defn atan
  "atan"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? x) (cond
                    (or (< y -1.0)
                        (> y 1.0)) (let [ip (inc y)
                                         im (dec y)]
                                     (Vec2. (m/copy-sign m/HALF_PI x)
                                            (* 0.25 (m/log (/ (* ip ip) (* im im))))))
                    (== y -1.0) (Vec2. 0.0 ##-Inf) ;; R and WolframAlpha result, was ##NaN as real part
                    (m/zero? y) z
                    (m/one? y) (Vec2. 0.0 ##Inf) ;; R and WolframAlpha result
                    :else (let [ip (inc y)
                                im (dec y)]
                            (Vec2. x (* 0.25 (m/log (/ (* ip ip) (* im im)))))))
      (m/zero? y) (Vec2. (m/atan x) y)
      :else (let [opi (inc y)
                  omi (- 1.0 y)
                  rr (* x x)]
              (Vec2. (* 0.5 (- (m/atan2 x omi)
                               (m/atan2 (- x) opi)))
                     (* 0.25 (m/log (/ (+ rr (* opi opi))
                                       (+ rr (* omi omi))))))))))

(defn asec
  "asec"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (if (m/zero? y)
      (cond
        (or (<= x -1.0)
            (>= x 1.0)) (Vec2. (m/acos (/ x)) y)
        (m/zero? x) (Vec2. ##Inf ##Inf)
        (neg? x) (Vec2. m/PI (* (m/copy-sign 1.0 y) (m/log (+ (/ -1.0 x)
                                                              (m/sqrt (dec (/ (* x x))))))))
        :else (Vec2. 0.0 (* (m/copy-sign 1.0 y) (m/log (+ (/ x)
                                                          (m/sqrt (dec (/ (* x x)))))))))
      (let [ss (+ (* x x) (* y y))
            p1 (/ x ss)
            p1p1 (m/sq (inc p1))
            p1m1 (m/sq (dec p1))
            p2 (m/sq (/ y ss))
            p1p1p2 (m/sqrt (+ p1p1 p2))
            p1m1p2 (m/sqrt (+ p1m1 p2))
            sum (* 0.5 (+ p1p1p2 p1m1p2))]
        (Vec2. (m/acos (* 0.5 (- p1p1p2 p1m1p2)))
               (* (m/sgn y) (m/log (+ sum (m/sqrt (dec (* sum sum)))))))))))


(defn acsc
  "acsc"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y)  (cond
                     (or (<= x -1.0)
                         (>= x 1.0)) (Vec2. (m/asin (/ x)) (- y))
                     (m/zero? x) (Vec2. ##Inf ##Inf)
                     (neg? x) (Vec2. m/-HALF_PI (* -1.0 (m/copy-sign 1.0 y) (m/log (+ (/ -1.0 x)
                                                                                      (m/sqrt (dec (/ (* x x))))))))
                     :else (Vec2. m/HALF_PI (* -1.0 (m/copy-sign 1.0 y) (m/log (+ (/ x)
                                                                                  (m/sqrt (dec (/ (* x x)))))))))
      (m/zero? x) (Vec2. x (- (m/log (+ (/ y) (m/sqrt (inc (/ (* y y))))))))
      :else (let [ss (+ (* x x) (* y y))
                  p1 (/ x ss)
                  p1p1 (m/sq (inc p1))
                  p1m1 (m/sq (dec p1))
                  p2 (m/sq (/ y ss))
                  p1p1p2 (m/sqrt (+ p1p1 p2))
                  p1m1p2 (m/sqrt (+ p1m1 p2))
                  sum (* 0.5 (+ p1p1p2 p1m1p2))]
              (Vec2. (m/asin (* 0.5 (- p1p1p2 p1m1p2)))
                     (* -1.0 (m/sgn y) (m/log (+ sum (m/sqrt (dec (* sum sum)))))))))))


(defn acot
  "acot"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (- (m/copy-sign m/HALF_PI x) (m/atan x)) (- y))
      (m/zero? x) (cond
                    (or (< y -1.0)
                        (> y 1.0)) (Vec2. (if (m/negative-zero? x) m/PI 0.0)
                                          (* 0.5 (m/log (/ (- 1.0 y)
                                                           (- -1.0 y)))))
                    (== y -1.0) (Vec2. 0.0 ##Inf)
                    (== y 1.0) (Vec2. 0.0 ##-Inf)
                    :else (Vec2. (m/copy-sign m/HALF_PI x) (* -0.25 (m/log (/ (m/sq (inc y))
                                                                              (m/sq (dec y)))))))
      :else (let [opi (+ 1.0 y)
                  omi (- 1.0 y)
                  rr (* x x)]
              (Vec2. (+ (m/copy-sign m/HALF_PI x) (* 0.5 (- (m/atan2 (- x) opi) (m/atan2 x omi))))
                     (* -0.25 (m/log (/ (+ rr (* opi opi))
                                        (+ rr (* omi omi))))))))))

;;

(defn asinh
  "asinh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? x) (cond
                    (m/between? -1.0 1.0 y) (Vec2. x (m/asin y))
                    (> y 1.0) (Vec2. (* (m/copy-sign 1.0 x)
                                        (m/log (+ (m/sqrt (dec (* y y))) y)))
                                     m/HALF_PI)
                    :else (Vec2. (* (m/copy-sign 1.0 x)
                                    (m/log (- (m/sqrt (dec (* y y))) y)))
                                 m/-HALF_PI))
      (m/zero? y) (Vec2. (m/log (+ x (m/sqrt (inc (* x x))))) y)
      :else (let [ss (inc (+ (* x x) (* y y)))
                  i2 (* 2.0 y)
                  ssp2i (m/sqrt (+ ss i2))
                  ssm2i (m/sqrt (- ss i2))
                  sum (* 0.5 (+ ssp2i ssm2i))]
              (Vec2. (* (m/sgn x) (m/log (+ sum (m/sqrt (dec (* sum sum))))))
                     (m/asin (* 0.5 (- ssp2i ssm2i))))))))

(defn acosh
  "acosh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (cond
                    (> x 1.0) (Vec2. (* (m/copy-sign 1.0 y) (m/log (+ x (m/sqrt (* (dec x) (inc x)))))) y)
                    (m/one? x) (Vec2. y 0.0)
                    (m/zero? x) (Vec2. y m/HALF_PI)
                    (> x -1.0) (Vec2. y (m/acos x))
                    (== x -1.0) (Vec2. y m/PI)
                    :else (Vec2. (* (m/copy-sign 1.0 y) (m/log (- (m/sqrt (* (- -1.0 x) (- 1.0 x))) x)))
                                 m/PI))
      (m/zero? x) (Vec2. (* (m/copy-sign 1.0 y)
                            (m/log (+ (m/sqrt (inc (* y y))) (m/abs y))))
                         m/HALF_PI)
      :else (let [ss (inc (+ (* x x) (* y y)))
                  r2 (* 2.0 x)
                  ssp2r (m/sqrt (+ ss r2))
                  ssm2r (m/sqrt (- ss r2))
                  sum (* 0.5 (+ ssp2r ssm2r))]
              (Vec2. (* (m/sgn y)
                        (m/log (+ sum (m/sqrt (dec (* sum sum))))))
                     (m/acos (* 0.5 (- ssp2r ssm2r))))))))

(defn atanh
  "atanh"
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (cond
                    (or (> x 1.0)
                        (< x -1.0)) (let [rp (inc x)
                                          rm (dec x)]
                                      (Vec2. (* 0.25 (m/log (/ (* rp rp) (* rm rm))))
                                             (m/copy-sign m/HALF_PI y)))
                    (m/one? x) (Vec2. ##Inf ##NaN)
                    (m/pos? x) (Vec2. (* 0.5 (m/log (/ (inc x) (- 1.0 x)))) y)
                    (m/zero? x) z
                    (> x -1.0) (Vec2. (* -0.5 (m/log (/ (- 1.0 x) (inc x)))) y)
                    :else (Vec2. ##-Inf ##NaN))
      (m/zero? x) (Vec2. x (m/atan y))
      :else (let [opr (inc x)
                  omr (- 1.0 x)
                  ii (* y y)]
              (Vec2. (* 0.25 (m/log (/ (+ ii (* opr opr))
                                       (+ ii (* omr omr)))))
                     (* 0.5 (- (m/atan2 y opr) (m/atan2 (- y) omr))))))))


(defn asech
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (cond
                    (or (< x -1.0)
                        (> x 1.0)) (Vec2. 0.0 (m/acos (/ x)))
                    (== x -1.0) (Vec2. y m/PI)
                    (m/neg? x) (let [rr (/ x)]
                                 (Vec2. (* -1.0 (m/copy-sign 1.0 y)
                                           (m/log (- (m/sqrt (* (dec rr) (inc rr))) rr)))
                                        m/PI))
                    (m/zero? x) (Vec2. ##Inf ##NaN)
                    :else (let [rr (/ x)]
                            (Vec2. (* -1.0 (m/copy-sign 1.0 y)
                                      (m/log (+ (m/sqrt (* (dec rr) (inc rr))) rr)))
                                   0.0)))
      (m/zero? x) (let [rr (/ y)
                        srr (m/sqrt (inc (* rr rr)))]
                    (Vec2. (if (m/pos? y)
                             (- (m/log (+ srr rr)))
                             (m/log (- srr rr))) m/HALF_PI))
      :else (let [ss (+ (* x x) (* y y))
                  Rr (/ x ss)
                  Rrp1 (m/sq (inc Rr))
                  Rrm1 (m/sq (dec Rr))
                  Ri (m/sq (/ y ss))
                  Rrp1pRi (m/sqrt (+ Rrp1 Ri))
                  Rrm1pRi (m/sqrt (+ Rrm1 Ri))
                  sum (* 0.5 (+ Rrp1pRi Rrm1pRi))]
              (Vec2. (* (m/sgn (- y)) (m/log (+ sum (m/sqrt (dec (* sum sum))))))
                     (m/acos (* 0.5 (- Rrp1pRi Rrm1pRi))))))))

(defn acsch
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? x) (cond
                    (or (< y -1.0)
                        (> y 1.0)) (Vec2. x (- (m/asin (/ y))))
                    (== y -1.0) (Vec2. x m/HALF_PI)
                    (m/neg? y) (let [ri (/ y)]
                                 (Vec2. (* (m/copy-sign 1.0 x)
                                           (m/log (- (m/sqrt (dec (* ri ri))) ri)))
                                        m/HALF_PI))
                    (m/zero? y) (Vec2. ##Inf ##-Inf)
                    (m/one? y) (Vec2. x m/-HALF_PI)
                    :else (let [ri (/ y)]
                            (Vec2. (* (m/copy-sign 1.0 x)
                                      (m/log (+ (m/sqrt (dec (* ri ri))) ri)))
                                   m/-HALF_PI)))
      (m/zero? y) (let [rr (/ x)]
                    (Vec2. (m/log (+ rr (m/sqrt (inc (* rr rr))))) (- y)))
      :else (let [ss (+ (* x x) (* y y))
                  p1 (/ y ss)
                  p1p1 (m/sq (inc p1))
                  p1m1 (m/sq (dec p1))
                  p2 (m/sq (/ x ss))
                  p1p1p2 (m/sqrt (+ p1p1 p2))
                  p1m1p2 (m/sqrt (+ p1m1 p2))
                  sum (* 0.5 (+ p1p1p2 p1m1p2))]
              (Vec2. (* (m/sgn x) (m/log (+ sum (m/sqrt (dec (* sum sum))))))
                     (m/asin (* 0.5 (- p1m1p2 p1p1p2))))))))


(defn acoth
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (cond
                    (or (< x -1.0)
                        (> x 1.0)) (Vec2. (* 0.5 (m/log (/ (inc x) (dec x))))
                                          (- y))
                    (== x -1.0) (Vec2. ##-Inf ##NaN)
                    (m/neg? x) (Vec2. (* 0.5 (m/log (/ (inc x) (- 1.0 x))))
                                      (- (m/copy-sign m/HALF_PI y)))
                    (m/zero? x) (Vec2. x (- (m/copy-sign m/HALF_PI y)))
                    (< x 1.0) (Vec2. (* 0.5 (m/log (/ (inc x) (- 1.0 x))))
                                     (- (m/copy-sign m/HALF_PI y)))
                    :else (Vec2. ##Inf ##NaN))
      (m/zero? x) (if (m/pos? y)
                    (Vec2. x (- (m/atan y) m/HALF_PI))
                    (Vec2. x (+ (m/atan y) m/HALF_PI)))
      :else (let [rp1 (inc x)
                  rm1 (dec x)
                  ii (* y y)]
              (Vec2. (* 0.25 (m/log (/ (+ ii (* rp1 rp1))
                                       (+ ii (* rm1 rm1)))))
                     (* 0.5 (- (m/atan2 y rp1)
                               (m/atan2 y rm1))))))))

(m/unuse-primitive-operators #{'abs 'zero?})
