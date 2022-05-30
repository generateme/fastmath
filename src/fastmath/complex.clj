;; # Namespace scope
;;
;; Functions to manipulate Vec2 as Complex numbers.
;; Implementation based on Apache Commons Math

(ns fastmath.complex
  "Complex numbers functions.

  Complex number is represented as `Vec2` type (from [[clojure2d.math.vector]] namespace).

  To create complex number use [[complex]], [[vec2]] or [[->Vec2]].

  Simplified implementation based on Apache Commons Math. Functions don't check NaNs or INF values.
  
  Complex plane (identity) looks as follows:

  ![identity](images/c/identity.jpg)" 
  {:metadoc/categories {:trig "Trigonometry"
                        :pow "Power / logarithm"}}
  (:refer-clojure :exclude [abs])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators #{'abs})

(def I (Vec2. 0.0 1.0))
(def I- (Vec2. 0.0 -1.0))
(def -I (Vec2. 0.0 -1.0))
(def ONE (Vec2. 1.0 0.0))
(def TWO (Vec2. 2.0 0.0))
(def ZERO (Vec2. 0.0 0.0))

(defn complex
  "Create complex number. Represented as `Vec2`."
  ([^double a ^double b] (Vec2. a b))
  ([^double a] (Vec2. a 0.0))
  ([] ZERO))

(def ^{:doc "Absolute value"} abs v/mag)
(def ^{:doc "Sum of two complex numbers."} add v/add)
(def ^{:doc "Subtraction of two complex numbers."} sub v/sub)
(def ^{:doc "Argument (angle) of complex number."} arg v/heading)
(def ^{:doc "Scale number"} scale v/mult)

(defn re "Real part" ^double [^Vec2 z] (.x z))
(defn im "Imaginary part" ^double [^Vec2 z] (.y z))

(defn flip "Exchange imaginary and real parts" ^Vec2 [^Vec2 z] (Vec2. (.y z) (.x z)))

(defn conjugate
  "Complex conjugate. \\\\(\\bar{z}\\\\)"
  [^Vec2 z] 
  (Vec2. (.x z) (- (.y z))))

(defn div
  "Divide two complex numbers."
  [^Vec2 z1 ^Vec2 z2]
  (let [a (.x z1)
        b (.y z1)
        c (.x z2)
        d (.y z2)
        den (+ (* c c) (* d d))]
    (if (zero? den)
      ZERO
      (Vec2. (/ (+ (* a c) (* b d)) den)
             (/ (- (* b c) (* a d)) den)))))

(defn reciprocal
  "\\\\(\\frac{1}{z}\\\\)"
  [z]
  (div ONE z))

;; [[../../docs/images/c/reciprocal.jpg]]

(defn mult
  "Multiply two complex numbers."
  [^Vec2 z1 ^Vec2 z2]
  (let [a (.x z1)
        b (.y z1)
        c (.x z2)
        d (.y z2)]
    (Vec2. (- (* a c) (* b d))
           (+ (* a d) (* b c)))))

(defn mult-I [^Vec2 z] (Vec2. (- (.y z)) (.x z)))
(defn mult-I- [^Vec2 z] (Vec2. (.y z) (- (.x z))))

(defn neg
  "Negate complex number. \\\\(-z\\\\)"
  [z]
  (v/sub z))

(defn sq
  "Square complex number. \\\\(z^2\\\\)"
  [z]
  (mult z z))

;; [[../../docs/images/c/sq.jpg]]

(defn sqrt
  "Sqrt of complex number. \\\\(\\sqrt{z}\\\\)"
  [^Vec2 z]
  (let [x (.x z)
        y (.y z)
        ^double l (abs z)
        xx (m/sqrt (+ l x))
        yy (* (m/signum y) (m/sqrt (- l x)))]
    (Vec2. (* m/SQRT2_2 xx) (* m/SQRT2_2 yy))))

;; [[../../docs/images/c/sqrt.jpg]]

(defn sqrt1z
  "\\\\(\\sqrt{1-z^2}\\\\)"
  [z]
  (->> z
       (mult z)
       (sub ONE)
       (sqrt)))

;; [[../../docs/images/c/sqrt1z.jpg]]

(defn cos
  "cos"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (Vec2. (* (m/cos x) (m/cosh y))
           (* (- (m/sin x)) (m/sinh y)))))

;; [[../../docs/images/c/cos.jpg]]

(defn sin
  "sin"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (Vec2. (* (m/sin x) (m/cosh y))
           (* (m/cos x) (m/sinh y)))))

;; [[../../docs/images/c/sin.jpg]]

(defn cosh
  "cosh"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (Vec2. (* (m/cosh x) (m/cos y))
           (* (m/sinh x) (m/sin y)))))

;; [[../../docs/images/c/cosh.jpg]]

(defn sinh
  "sinh"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (Vec2. (* (m/sinh x) (m/cos y))
           (* (m/cosh x) (m/sin y)))))

;; [[../../docs/images/c/sinh.jpg]]

(defn tan
  "tan"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [aa (* 2.0 (.x z))
        bb (* 2.0 (.y z))
        cc (+ (m/cos aa) (m/cosh bb))]
    (Vec2. (/ (m/sin aa) cc)
           (/ (m/sinh bb) cc))))

;; [[../../docs/images/c/tan.jpg]]

(defn tanh
  "tanh"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [aa (* 2.0 (.x z))
        bb (* 2.0 (.y z))
        cc (+ (m/cosh aa) (m/cos bb))]
    (Vec2. (/ (m/sinh aa) cc)
           (/ (m/sin bb) cc))))

;; [[../../docs/images/c/tanh.jpg]]

(defn sec
  "secant"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [cc (+ (m/cos (* 2.0 (.x z)))
              (m/cosh (* 2.0 (.y z))))]
    (Vec2. (/ (* 2.0 (m/cos (.x z)) (m/cosh (.y z))) cc)
           (/ (* 2.0 (m/sin (.x z)) (m/sinh (.y z))) cc))))

;; [[../../docs/images/c/sec.jpg]]

(defn csc
  "cosecant"
  {:metadoc/categories #{:trig}}
  [^Vec2 z]
  (let [cc (- (m/cos (* 2.0 (.x z)))
              (m/cosh (* 2.0 (.y z))))]
    (Vec2. (- (/ (* 2.0 (m/cosh (.y z)) (m/sin (.x z))) cc))
           (/ (* 2.0 (m/cos (.x z)) (m/sinh (.y z))) cc))))

;; [[../../docs/images/c/csc.jpg]]

(defn exp
  "exp"
  {:metadoc/categories #{:pow}}
  [^Vec2 z]
  (let [e (m/exp (.x z))
        y (.y z)]
    (Vec2. (* e (m/cos y))
           (* e (m/sin y)))))

;; [[../../docs/images/c/exp.jpg]]

(defn log
  "log"
  {:metadoc/categories #{:pow}}
  [^Vec2 z]
  (Vec2. (m/log (abs z))
         (m/atan2 (.y z) (.x z))))

;; [[../../docs/images/c/log.jpg]]

(defn acos
  "acos"
  {:metadoc/categories #{:trig}}
  [z]
  (->> (sqrt1z z)
       (mult-I)
       (add z)
       (log)
       (mult-I-)))

;; [[../../docs/images/c/acos.jpg]]

(defn asin
  "asin"
  {:metadoc/categories #{:trig}}
  [z]
  (->> (sqrt1z z)
       (add (mult-I z))
       (log)
       (mult-I-)))

;; [[../../docs/images/c/asin.jpg]]

(def ^:private i-div-two (div I TWO))

(defn atan
  "atan"
  {:metadoc/categories #{:trig}}
  [z]
  (->> (sub I z)
       (div (add I z))
       (log)
       (mult i-div-two)))

;; [[../../docs/images/c/atan.jpg]]

(defn pow
  "Power. \\\\(z_1^{z_2}\\\\)"
  {:metadoc/categories #{:pow}}
  [z1 z2]
  (->> z1
       (log)
       (mult z2)
       (exp)))

(defn acosh
  "acosh"
  {:metadoc/categories #{:trig}}
  [z]
  (-> z
      (acos)
      (mult-I)))

(defn acosech
  "acosech"
  {:metadoc/categories #{:trig}}
  [z]
  (-> z
      (reciprocal)
      (acosh)))

(defn asinh
  "asinh"
  {:metadoc/categories #{:trig}}
  [z]
  (-> z
      (mult-I)
      (asin)
      (mult-I-)))

(defn atanh
  "atanh"
  {:metadoc/categories #{:trig}}
  [z]
  (-> z
      (mult-I)
      (atan)
      (mult-I-)))

(defn sech
  "sech"
  {:metadoc/categories #{:trig}}
  [z]
  (reciprocal (cosh z)))

(defn csch
  "cosech"
  {:metadoc/categories #{:trig}}
  [z]
  (reciprocal (sinh z)))

(defn coth
  "coth"
  {:metadoc/categories #{:trig}}
  [z]
  (reciprocal (tanh z)))

(defn asech
  "sech"
  {:metadoc/categories #{:trig}}
  [z]
  (acosh (reciprocal z)))

(defn acsch
  "sech"
  {:metadoc/categories #{:trig}}
  [z]
  (asinh (reciprocal z)))

(defn acoth
  "sech"
  {:metadoc/categories #{:trig}}
  [z]
  (atanh (reciprocal z)))

