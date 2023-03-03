(ns fastmath.easings
  "Easing functions.

  List of all are in [[easings-list]]."
  {:metadoc/categories {:lin "Linear"
                        :back "Anticipatory easings"
                        :bounce "Bounce"
                        :circle "Circular"
                        :cubic "Cubic"
                        :elastic "Elastic"
                        :exp "Exponential"
                        :poly "Polynomial"
                        :quad "Quadratic"
                        :sin "Sinusoidal"
                        :cr "Creators"}}
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;;

;; collect easings names and functions
(defonce ^:private easings-list-atom (atom {}))
(defmacro ^:private add-name [s] `(swap! easings-list-atom assoc ~(keyword s) ~s))

(reset! easings-list-atom {})

;;

(def ^{:const true :private true :tag 'double} ease-back-s 1.70158)

;;

(defn linear
  "Linear easing (identity)"
  {:metadoc/categories #{:lin}}
  ^double [^double t] t)

(defn back-in
  "BackIn easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}}
  (^double [^double t] (back-in ease-back-s t))
  (^double [^double s ^double t]
   (* t t (- (* (inc s) t) s))))
(add-name back-in)

(defn back-out
  "BackOut easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}}
  (^double [^double t] (back-out ease-back-s t))
  (^double [^double s ^double t]
   (let [t- (dec t)]
     (inc (* t- t- (+ (* (inc s) t-) s))))))
(add-name back-out)

(defn back-in-out
  "BackInOut easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}}
  (^double [^double t] (back-in-out ease-back-s t))
  (^double [^double s ^double t]
   (let [t2 (+ t t)]
     (if (< t2 1.0)
       (* 0.5 t2 t2 (- (* (inc s) t2) s))
       (let [t- (- t2 2.0)]
         (/ (+ 2.0 (* t- t- (+ (* (inc s) t-) s))) 2.0))))))
(add-name back-in-out)

;;

(def ^{:const true :private true :tag 'double} b1 (/ 4.0 11.0))
(def ^{:const true :private true :tag 'double} b2 (/ 6.0 11.0))
(def ^{:const true :private true :tag 'double} b3 (/ 8.0 11.0))
(def ^{:const true :private true :tag 'double} b4 (/ 3.0 4.0))
(def ^{:const true :private true :tag 'double} b5 (/ 9.0 11.0))
(def ^{:const true :private true :tag 'double} b6 (/ 10.0 11.0))
(def ^{:const true :private true :tag 'double} b7 (/ 15.0 16.0))
(def ^{:const true :private true :tag 'double} b8 (/ 21.0 22.0))
(def ^{:const true :private true :tag 'double} b9 (/ 63.0 64.0))
(def ^{:const true :private true :tag 'double} b0 (/ (/ b1) b1))

(defn bounce-out
  "BounceOut easing"
  {:metadoc/categories #{:bounce}}
  ^double [^double t]
  (cond
    (< t b1) (* b0 t t)
    (< t b3) (let [t- (- t b2)] (+ b4 (* b0 t- t-)))
    (< t b6) (let [t- (- t b5)] (+ b7 (* b0 t- t-)))
    :else (let [t- (- t b8)] (+ b9 (* b0 t- t-)))))
(add-name bounce-out)

(defn bounce-in
  "BounceIn easing"
  {:metadoc/categories #{:bounce}}
  ^double [^double t]
  (- 1.0 (bounce-out (- 1.0 t))))
(add-name bounce-in)

(defn bounce-in-out
  "BounceInOut easing"
  {:metadoc/categories #{:bounce}}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (- 1.0 (bounce-out (- 1.0 t2)))
             (inc (bounce-out (dec t2)))))))
(add-name bounce-in-out)

;;

(defn circle-in
  "CircleIn easing"
  {:metadoc/categories #{:circle}}
  ^double [^double t]
  (- 1.0 (m/sqrt (- 1.0 (* t t)))))
(add-name circle-in)

(defn circle-out
  "CircleIn easing"
  {:metadoc/categories #{:circle}}
  ^double [^double t]
  (let [t- (dec t)]
    (m/sqrt (- 1.0 (* t- t-)))))
(add-name circle-out)

(defn circle-in-out
  "CircleInOut easing"
  {:metadoc/categories #{:circle}}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (- 1.0 (m/sqrt (- 1.0 (* t2 t2))))
             (let [t- (- t2 2.0)]
               (inc (m/sqrt (- 1.0 (* t- t-)))))))))
(add-name circle-in-out)

;;

(defn cubic-in
  "CubicIn easing"
  {:metadoc/categories #{:cubic}}
  ^double [^double t]
  (* t t t))
(add-name cubic-in)

(defn cubic-out
  "CubicOut easing"
  {:metadoc/categories #{:cubic}}
  ^double [^double t]
  (let [t- (dec t)]
    (inc (* t- t- t-))))
(add-name cubic-out)

(defn cubic-in-out
  "CubicInOut easing"
  {:metadoc/categories #{:cubic}}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (* t2 t2 t2)
             (let [t- (- t2 2.0)]
               (+ 2.0 (* t- t- t-)))))))
(add-name cubic-in-out)

;;

(defn- calc-elastic-factor
  "Calculate elastic factor for given `amplitude` (default: 1.0) and `period` (default: 0.3)."
  (^double [] (calc-elastic-factor 1.0 0.3))
  (^double [^double amplitude ^double period]
   (* (/ period m/TWO_PI) (m/asin (/ (max 1.0 amplitude))))))

(def ^{:const true :private true :tag 'double} elastic-s (calc-elastic-factor))
(def ^{:const true :private true :tag 'double} elastic-period (/ 0.3 m/TWO_PI))

(defn elastic-in
  "ElasticIn.

  When called with `t` parameter, returns easing value (for `amplitude=1.0` and `period=0.3`).
  When called with `amplitude` and `period` returns custom easing function."
  {:metadoc/categories #{:elastic}}
  (^double [^double t]
   (let [t- (dec t)]
     (* (m/pow 2.0 (* 10.0 t-)) (m/sin (/ (- elastic-s t-) elastic-period)))))
  ([^double amplitude ^double period]
   (let [p (/ period m/TWO_PI)
         s (calc-elastic-factor amplitude period)]
     (fn ^double [^double t]
       (let [t- (dec t)] 
         (* amplitude (m/pow 2.0 (* 10.0 t-)) (m/sin (/ (- s t-) p))))))))
(add-name elastic-in)

(defn elastic-out
  "ElasticOut.

  When called with `t` parameter, returns easing value (for `amplitude=1.0` and `period=0.3`).
  When called with `amplitude` and `period` returns custom easing function."
  {:metadoc/categories #{:elastic}}
  (^double [^double t]
   (- 1.0 (* (m/pow 2.0 (* -10.0 t)) (m/sin (/ (+ elastic-s t) elastic-period)))))
  ([^double amplitude ^double period]
   (let [p (/ period m/TWO_PI)
         s (calc-elastic-factor amplitude period)]
     (fn ^double [^double t]
       (- 1.0 (* amplitude (m/pow 2.0 (* -10.0 t)) (m/sin (/ (+ s t) p))))))))
(add-name elastic-out)

(defn elastic-in-out
  "ElasticInOut.

  When called with `t` parameter, returns easing value (for `amplitude=1.0` and `period=0.3`).
  When called with `amplitude` and `period` returns custom easing function."
  {:metadoc/categories #{:elastic}}
  (^double [^double t]
   (let [t2 (dec (+ t t))]
     (* 0.5 (if (neg? t2)
              (* (m/pow 2.0 (* 10.0 t2)) (m/sin (/ (- elastic-s t2) elastic-period)))
              (- 2.0 (* (m/pow 2.0 (* -10.0 t2)) (m/sin (/ (+ elastic-s t2) elastic-period))))))))
  ([^double amplitude ^double period]
   (let [p (/ period m/TWO_PI)
         s (calc-elastic-factor amplitude period)]
     (fn ^double [^double t]
       (let [t2 (dec (+ t t))]
         (* 0.5 (if (neg? t2)
                  (* amplitude (m/pow 2.0 (* 10.0 t2)) (m/sin (/ (- s t2) p)))
                  (- 2.0 (* amplitude (m/pow 2.0 (* -10.0 t2)) (m/sin (/ (+ s t2) p)))))))))))
(add-name elastic-in-out)

;;

(defn exp-in
  "ExpIn easing"
  {:metadoc/categories #{:exp}}
  ^double [^double t]
  (m/pow 2.0 (- (* 10.0 t) 10.0)))
(add-name exp-in)

(defn exp-out
  "ExpOut easing"
  {:metadoc/categories #{:exp}}
  ^double [^double t]
  (- 1.0 (m/pow 2.0 (* -10.0 t))))
(add-name exp-out)

(defn exp-in-out
  "ExpInOut easing"
  {:metadoc/categories #{:exp}}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (m/pow 2.0 (- (* 10.0 t2) 10.0))
             (- 2.0 (m/pow 2.0 (- 10.0 (* 10.0 t2))))))))
(add-name exp-in-out)

;;

(defn poly-in
  "PolyIn easing.

  Optional exponent `e` defaults to 3.0."
  {:metadoc/categories #{:poly}}
  (^double [^double t] (m/pow t 3.0))
  (^double [^double e ^double t] (m/pow t e)))
(add-name poly-in)

(defn poly-out
  "PolyOut easing.

  Optional exponent `e` defaults to 3.0."
  {:metadoc/categories #{:poly}}
  (^double [^double t] (- 1.0 (m/pow (- 1.0 t) 3.0)))
  (^double [^double e ^double t] (- 1.0 (m/pow (- 1.0 t) e))))
(add-name poly-out)

(defn poly-in-out
  "PolyInOut easing.

  Optional exponent `e` defaults to 3.0."
  {:metadoc/categories #{:poly}}
  (^double [^double t]
   (let [t2 (+ t t)]
     (* 0.5 (if (< t 0.5)
              (m/pow t2 3.0)
              (- 2.0 (m/pow (- 2.0 t2) 3.0))))))
  (^double [^double e ^double t]
   (let [t2 (+ t t)]
     (* 0.5 (if (< t 0.5)
              (m/pow t2 e)
              (- 2.0 (m/pow (- 2.0 t2) e)))))))
(add-name poly-in-out)

;;

(defn quad-in
  "QuadIn easing"
  {:metadoc/categories #{:quad}}
  ^double [^double t]
  (* t t))
(add-name quad-in)

(defn quad-out
  "QuadOut easing"
  {:metadoc/categories #{:quad}}
  ^double [^double t]
  (* t (- 2.0 t)))
(add-name quad-out)

(defn quad-in-out
  "QuadInOut easing"
  {:metadoc/categories #{:quad}}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (* t2 t2)
             (let [t- (dec t2)]
               (inc (* t- (- 2.0 t-))))))))
(add-name quad-in-out)

;;


(defn sin-in
  "SinIn easing"
  {:metadoc/categories #{:sin}}
  ^double [^double t]
  (- 1.0 (m/cos (* t m/HALF_PI))))
(add-name sin-in)

(defn sin-out
  "SinOut easing"
  {:metadoc/categories #{:sin}}
  ^double [^double t]
  (m/sin (* t m/HALF_PI)))
(add-name sin-out)

(defn sin-in-out
  "SinInOut easing"
  {:metadoc/categories #{:sin}}
  ^double [^double t]
  (* 0.5 (- 1.0 (m/cos (* t m/PI)))))
(add-name sin-in-out)


;;

(defn out
  "Create out easing for given `easing` function."
  {:metadoc/categories #{:cr}}
  [easeing]
  (fn ^double [^double t]
    (- 1.0 ^double (easeing (- 1.0 t)))))

(defn in-out
  "Create in-out easing for given `easing` function."
  {:metadoc/categories #{:cr}}
  [easeing]
  (fn ^double [^double t]
    (* 0.5 (if (< t 0.5)
             ^double (easeing (* 2.0 t))
             (- 2.0 ^double (easeing (* 2.0 (- 1.0 t))))))))

(defn reflect
  "Create in-out easing for given `easing` function and `center`."
  {:metadoc/categories #{:cr}}
  [easing ^double center]
  (fn ^double [^double t]
    (if (< t center)
      (* center ^double (easing (/ t center)))
      (let [center- (- 1.0 center)]
        (- 1.0 (* center- ^double (easing (/ (- 1.0 t) center-))))))))

;;

(def ^{:doc "Map of easing names (as keywords) and functions."}
  easings-list @easings-list-atom)
