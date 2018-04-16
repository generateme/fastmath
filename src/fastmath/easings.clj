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
  (:require [fastmath.core :as m]
            [metadoc.examples :refer :all]
            [incanter.charts :as c]
            [incanter.core :as i]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;;

;; collect easings names and functions
(defonce ^:private easings-list-atom (atom {}))
(defmacro ^:private add-name [s] `(swap! easings-list-atom assoc ~(keyword s) ~s))

;;

(def ^:private ^:const ^double ease-back-s 1.70158)

(defsnippet save-graph
  "Save incanter graph"
  (let [fname (str "images/e/" (first opts) ".png")]
    (i/save (c/function-plot f 0.0 1.0 :y-label "easing value") (str "docs/" fname) :width 500 :height 250)
    fname))

;;

(defn linear
  "Linear easing (identity)"
  {:metadoc/categories #{:lin}
   :metadoc/examples [(example "Usage" {:test-value 0.5} (linear 0.5))
                      (example-snippet "Plot" save-graph :image linear)]}
  ^double [^double t] t)

(defn back-in
  "BackIn easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}
   :metadoc/examples [(example "Usage" {:test-value -0.087698} (m/approx (back-in 0.5) 6))
                      (example "Different overshot" (back-in 2.0 0.5))
                      (example-snippet "Plot" save-graph :image back-in)
                      (example-snippet "Plot, different `s` value." save-graph :image (partial back-in 5.0))]}
  (^double [^double t] (back-in ease-back-s t))
  (^double [^double s ^double t]
   (* t t (- (* (inc s) t) s))))
(add-name back-in)

(defn back-out
  "BackOut easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}
   :metadoc/examples [(example "Usage" {:test-value (- 1.0 -0.087698)} (m/approx (back-out 0.5) 6))
                      (example "Different overshot" (back-out 2.0 0.5))
                      (example-snippet "Plot" save-graph :image back-out)
                      (example-snippet "Plot, different `s` value." save-graph :image (partial back-out 5.0))]}
  (^double [^double t] (back-out ease-back-s t))
  (^double [^double s ^double t]
   (let [t- (dec t)]
     (inc (* t- t- (+ (* (inc s) t-) s))))))
(add-name back-out)

(defn back-in-out
  "BackInOut easing.

  Parameter `s` (default: 1.70158) defines overshoot."
  {:metadoc/categories #{:back}
   :metadoc/examples [(example "Usage" {:test-value -0.043849} (m/approx (back-in-out 0.25) 6))
                      (example "Different overshot" (back-in-out 2.0 0.75))
                      (example-snippet "Plot" save-graph :image back-in-out)
                      (example-snippet "Plot, different `s` value." save-graph :image (partial back-in-out 5.0))]}
  (^double [^double t] (back-in-out ease-back-s t))
  (^double [^double s ^double t]
   (let [t2 (+ t t)]
     (if (< t2 1.0)
       (* 0.5 t2 t2 (- (* (inc s) t2) s))
       (let [t- (- t2 2.0)]
         (/ (+ 2.0 (* t- t- (+ (* (inc s) t-) s))) 2.0))))))
(add-name back-in-out)

;;

(def ^:const ^:private ^double b1 (/ 4.0 11.0))
(def ^:const ^:private ^double b2 (/ 6.0 11.0))
(def ^:const ^:private ^double b3 (/ 8.0 11.0))
(def ^:const ^:private ^double b4 (/ 3.0 4.0))
(def ^:const ^:private ^double b5 (/ 9.0 11.0))
(def ^:const ^:private ^double b6 (/ 10.0 11.0))
(def ^:const ^:private ^double b7 (/ 15.0 16.0))
(def ^:const ^:private ^double b8 (/ 21.0 22.0))
(def ^:const ^:private ^double b9 (/ 63.0 64.0))
(def ^:const ^:private ^double b0 (/ (/ b1) b1))

(defn bounce-out
  "BounceOut easing"
  {:metadoc/categories #{:bounce}
   :metadoc/examples [(example "Usage" {:test-value 0.3025} (m/approx (bounce-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image bounce-out)]}
  ^double [^double t]
  (cond
    (< t b1) (* b0 t t)
    (< t b3) (let [t- (- t b2)] (+ b4 (* b0 t- t-)))
    (< t b6) (let [t- (- t b5)] (+ b7 (* b0 t- t-)))
    :else (let [t- (- t b8)] (+ b9 (* b0 t- t-)))))
(add-name bounce-out)

(defn bounce-in
  "BounceIn easing"
  {:metadoc/categories #{:bounce}
   :metadoc/examples [(example "Usage" {:test-value 0.06} (m/approx (bounce-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image bounce-in)]}
  ^double [^double t]
  (- 1.0 (bounce-out (- 1.0 t))))
(add-name bounce-in)

(defn bounce-in-out
  "BounceInOut easing"
  {:metadoc/categories #{:bounce}
   :metadoc/examples [(example "Usage" (m/approx (bounce-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image bounce-in-out)]}
  ^double [^double t]
  (let [t2 (+ t t)]
    (* 0.5 (if (< t 0.5)
             (- 1.0 (bounce-out (- 1.0 t2)))
             (inc (bounce-out (dec t2)))))))
(add-name bounce-in-out)

;;

(defn circle-in
  "CircleIn easing"
  {:metadoc/categories #{:circle}
   :metadoc/examples [(example "Usage" {:test-value 0.020204} (m/approx (circle-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image circle-in)]}
  ^double [^double t]
  (- 1.0 (m/sqrt (- 1.0 (* t t)))))
(add-name circle-in)

(defn circle-out
  "CircleIn easing"
  {:metadoc/categories #{:circle}
   :metadoc/examples [(example "Usage" {:test-value 0.6} (m/approx (circle-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image circle-out)]}
  ^double [^double t]
  (let [t- (dec t)]
    (m/sqrt (- 1.0 (* t- t-)))))
(add-name circle-out)

(defn circle-in-out
  "CircleInOut easing"
  {:metadoc/categories #{:circle}
   :metadoc/examples [(example "Usage" (m/approx (circle-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image circle-in-out)]}
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
  {:metadoc/categories #{:cubic}
   :metadoc/examples [(example "Usage" {:test-value 0.008} (m/approx (cubic-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image cubic-in)]}
  ^double [^double t]
  (* t t t))
(add-name cubic-in)

(defn cubic-out
  "CubicOut easing"
  {:metadoc/categories #{:cubic}
   :metadoc/examples [(example "Usage" {:test-value 0.488} (m/approx (cubic-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image cubic-out)]}
  ^double [^double t]
  (let [t- (dec t)]
    (inc (* t- t- t-))))
(add-name cubic-out)

(defn cubic-in-out
  "CubicInOut easing"
  {:metadoc/categories #{:cubic}
   :metadoc/examples [(example "Usage" (m/approx (cubic-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image cubic-in-out)]}
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

(def ^:private ^:const ^double elastic-s (calc-elastic-factor))
(def ^:private ^:const ^double elastic-period (/ 0.3 m/TWO_PI))

(defn elastic-in
  "ElasticIn.

  When called with `t` parameter, returns easing value (for `amplitude=1.0` and `period=0.3`).
  When called with `amplitude` and `period` returns custom easing function."
  {:metadoc/categories #{:elastic}
   :metadoc/examples [(example "Usage" {:test-value -0.001953} (m/approx (elastic-in 0.2) 6))
                      (example "Create custom elastic easing" (m/approx ((elastic-in 1.0 0.1) 0.2) 6))
                      (example-snippet "Plot" save-graph :image elastic-in)
                      (example-snippet "Plot, different settings" save-graph :image (elastic-in 1.1 0.1))]}
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
  {:metadoc/categories #{:elastic}
   :metadoc/examples [(example "Usage" {:test-value 1.125} (m/approx (elastic-out 0.2) 6))
                      (example "Create custom elastic easing" (m/approx ((elastic-out 1.0 0.1) 0.2) 6))
                      (example-snippet "Plot" save-graph :image elastic-out)
                      (example-snippet "Plot, different settings" save-graph :image (elastic-out 1.1 0.1))]}
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
  {:metadoc/categories #{:elastic}
   :metadoc/examples [(example "Usage" (m/approx (elastic-in-out 0.2) 6))
                      (example "Create custom elastic easing" (m/approx ((elastic-in-out 1.0 0.1) 0.2) 6))
                      (example-snippet "Plot" save-graph :image elastic-in-out)
                      (example-snippet "Plot, different settings" save-graph :image (elastic-in-out 1.1 0.1))]}
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
  {:metadoc/categories #{:exp}
   :metadoc/examples [(example "Usage" {:test-value 0.003906} (m/approx (exp-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image exp-in)]}
  ^double [^double t]
  (m/pow 2.0 (- (* 10.0 t) 10.0)))
(add-name exp-in)

(defn exp-out
  "ExpOut easing"
  {:metadoc/categories #{:exp}
   :metadoc/examples [(example "Usage" {:test-value 0.75} (m/approx (exp-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image exp-out)]}
  ^double [^double t]
  (- 1.0 (m/pow 2.0 (* -10.0 t))))
(add-name exp-out)

(defn exp-in-out
  "ExpInOut easing"
  {:metadoc/categories #{:exp}
   :metadoc/examples [(example "Usage" (m/approx (exp-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image exp-in-out)]}
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
  {:metadoc/categories #{:poly}
   :metadoc/examples [(example "Usage" {:test-value 0.008} (m/approx (poly-in 0.2) 6))
                      (example "Other exponent" (poly-in 8.0 0.5))
                      (example-snippet "Plot" save-graph :image poly-in)
                      (example-snippet "Plot" save-graph :image (partial poly-in 8.0))]}
  (^double [^double t] (m/pow t 3.0))
  (^double [^double e ^double t] (m/pow t e)))
(add-name poly-in)

(defn poly-out
  "PolyOut easing.

  Optional exponent `e` defaults to 3.0."
  {:metadoc/categories #{:poly}
   :metadoc/examples [(example "Usage" {:test-value 0.488} (m/approx (poly-out 0.2) 6))
                      (example "Other exponent" (poly-out 8.0 0.5))
                      (example-snippet "Plot" save-graph :image poly-out)
                      (example-snippet "Plot" save-graph :image (partial poly-out 8.0))]}
  (^double [^double t] (- 1.0 (m/pow (- 1.0 t) 3.0)))
  (^double [^double e ^double t] (- 1.0 (m/pow (- 1.0 t) e))))
(add-name poly-out)

(defn poly-in-out
  "PolyInOut easing.

  Optional exponent `e` defaults to 3.0."
  {:metadoc/categories #{:poly}
   :metadoc/examples [(example "Usage" (m/approx (poly-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image poly-in-out)
                      (example-snippet "Plot" save-graph :image (partial poly-in-out 8.0))]}
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
  {:metadoc/categories #{:quad}
   :metadoc/examples [(example "Usage" {:test-value 0.04} (m/approx (quad-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image quad-in)]}
  ^double [^double t]
  (* t t))
(add-name quad-in)

(defn quad-out
  "QuadOut easing"
  {:metadoc/categories #{:quad}
   :metadoc/examples [(example "Usage" {:test-value 0.36} (m/approx (quad-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image quad-out)]}
  ^double [^double t]
  (* t (- 2.0 t)))
(add-name quad-out)

(defn quad-in-out
  "QuadInOut easing"
  {:metadoc/categories #{:quad}
   :metadoc/examples [(example "Usage" (m/approx (quad-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image quad-in-out)]}
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
  {:metadoc/categories #{:sin}
   :metadoc/examples [(example "Usage" {:test-value 0.048943} (m/approx (sin-in 0.2) 6))
                      (example-snippet "Plot" save-graph :image sin-in)]}
  ^double [^double t]
  (- 1.0 (m/cos (* t m/HALF_PI))))
(add-name sin-in)

(defn sin-out
  "SinOut easing"
  {:metadoc/categories #{:sin}
   :metadoc/examples [(example "Usage" {:test-value 0.309017} (m/approx (sin-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image sin-out)]}
  ^double [^double t]
  (m/sin (* t m/HALF_PI)))
(add-name sin-out)

(defn sin-in-out
  "SinInOut easing"
  {:metadoc/categories #{:sin}
   :metadoc/examples [(example "Usage" (m/approx (sin-in-out 0.2) 6))
                      (example-snippet "Plot" save-graph :image sin-in-out)]}
  ^double [^double t]
  (* 0.5 (- 1.0 (m/cos (* t m/PI)))))
(add-name sin-in-out)


;;

(defn out
  "Create out easing for given `ease`."
  {:metadoc/categories #{:cr}
   :metadoc/examples [(example "Usage" (let [outeasing (out sin-in)]
                                         (== ^double (sin-out 0.75) ^double (outeasing 0.75))))
                      (example-snippet "Create out easing" save-graph :image (out sin-in-out))]}
  [ease]
  (fn ^double [^double t]
    (- 1.0 ^double (ease (- 1.0 t)))))

(defn in-out
  "Create in-out easing for given `ease`."
  {:metadoc/categories #{:cr}
   :metadoc/examples [(example "Usage" (let [inouteasing (in-out quad-in)]
                                         (== ^double (quad-in-out 0.75) ^double (inouteasing 0.75))))
                      (example-snippet "Create in-out easing" save-graph :image (in-out sin-in-out))]}
  [ease]
  (fn ^double [^double t]
    (* 0.5 (if (< t 0.5)
             ^double (ease (* 2.0 t))
             (- 2.0 ^double (ease (* 2.0 (- 1.0 t))))))))

(defn reflect
  "Create in-out easing for given `ease` function and `center`."
  {:metadoc/categories #{:cr}
   :metadoc/examples [(example "Usage" (let [neasing (reflect (partial back-in 2.0) 0.2)]
                                         [(neasing 0.1) (neasing 0.5) (neasing 0.9)]))
                      (example "For `center=0.5` function returns regular in-out easing"
                        {:test-value true} (== (back-in-out 0.4) ^double ((reflect back-in 0.5) 0.4)))
                      (example-snippet "Reflect one easing in center=0.2 " save-graph :image (reflect elastic-in-out 0.2))]}
  [ease ^double center]
  (fn ^double [^double t]
    (if (< t center)
      (* center ^double (ease (/ t center)))
      (let [center- (- 1.0 center)]
        (- 1.0 (* center- ^double (ease (/ (- 1.0 t) center-))))))))

;;

(def ^{:doc "Map of easing names (as keywords) and functions."
       :metadoc/examples [(example "List of easings" (sort (keys easings-list)))
                          (example "Access to function" ((easings-list :back-in) 3.0 0.25))]}
  easings-list @easings-list-atom)
