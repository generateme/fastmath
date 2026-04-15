(ns fastmath.signal.waveform
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- ->time-fn [v]  (if (number? v) (let [v (double v)] (fn ^double [^double t] (m/* t v))) v))
(defn- ->const-fn [v]  (if (number? v) (let [v (double v)] (fn ^double [_] v)) v))

(defn- exp-curve
  ^double [^double x]
  (m// (m/+ 3.0 (m/* x (m/+ -13.0 (m/* 5.0 x))))
       (m/+ 3.0 x x)))

(defn ->sine
  ([] (->sine nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (m/sin (m/* m/TWO_PI (m/+ (double (pfn t)) (double (ffn t))))))))))

(defn- quadratic-sine-approx
  ^double [^double v]
  (let [phase (m/frac v)
        hphase? (m/< phase 0.5)
        x (m/- phase (if hphase? 0.25 0.75))
        v (m/- 1.0 (m/* 16.0 x x))]
    (if hphase? v (m/- v))))

(defn ->analog-sine
  ([] (->analog-sine nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (quadratic-sine-approx (m/+ (double (pfn t)) (double (ffn t)))))))))

(defn ->square
  ([] (->square nil))
  ([{:keys [f phase amplitude duty]
     :or {f 1.0 phase 0.0 amplitude 1.0 duty 0.5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)
         wfn (->const-fn duty)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (if (m/< (m/frac (m/+ (double (pfn t)) (double (ffn t))))
                     (m/constrain (double (wfn t)) 0.0 1.0))
              1.0 -1.0))))))

(defn ->saw
  ([] (->saw nil))
  ([{:keys [f  phase amplitude up?]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (m/* (double (afn t))
              (m/dec (m/* 2.0 (m/frac (m/+ 0.5 (double (pfn t)) (double (ffn t))))))))
       (fn ^double [^double t]
         (m/* (double (afn t))
              (m/dec (m/* 2.0 (m/- 1.0 (m/frac (m/+ 0.5 (double (pfn t)) (double (ffn t)))))))))))))

(defn ->analog-saw
  ([] (->analog-saw nil))
  ([{:keys [f phase amplitude up?]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (let [x (m/frac (m/+ 0.255969349108945 (double (pfn t)) (double (ffn t))))]
           (m/* -1.0 (double (afn t)) (exp-curve x))))
       (fn ^double [^double t]
         (let [x (m/frac (m/+ 0.255969349108945 (double (pfn t)) (double (ffn t))))]
           (m/* (double (afn t)) (exp-curve x))))))))

(defn ->triangle
  ([] (->triangle nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (let [v (m/frac (m/+ (double (pfn t)) (double (ffn t))))]
         (m/* 4.0 (double (afn t))
              (double (cond
                        (m/< v 0.25) v
                        (m/< v 0.75) (m/- 0.5 v)
                        :else (m/dec v)))))))))

(defn- analog-triangle-approx
  ^double [^double phase]
  (let [x (-> (m/frac phase)
              (m/+ 0.1279846745544725)
              (m/frac))
        hx? (m/>= x 0.5)
        x (-> (m/* 2.0 x)
              (m/frac))]
    (if hx?
      (exp-curve x)
      (m/- (exp-curve x)))))

(defn ->analog-triangle
  ([] (->analog-triangle nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (analog-triangle-approx (m/+ (double (pfn t)) (double (ffn t)))))))))

;; band-limited

(defn- square-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/+ 1.0 b b)]                          
                     (m// (m/sin (m/* m/TWO_PI (double (m/* b+ t))))
                          b+)))) 0.0 (range bands)))

(defn ->band-limited-square
  ([] (->band-limited-square nil))
  ([{:keys [f phase amplitude ^long bands]
     :or {bands 5 f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m// (m/* (double (afn t))
                 (square-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
            m/QUARTER_PI)))))

(defn- saw-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/inc b)]                          
                     (-> (m/sin (m/* m/TWO_PI (double (m/* b+ t))))
                         (m// b+)
                         (m/* (if (m/even? b) 1.0 -1.0)))))) 0.0 (range bands)))

(defn ->band-limited-saw
  ([] (->band-limited-saw nil))
  ([{:keys [f phase amplitude up? ^long bands]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true bands 5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (m// (m/* (double (afn t)) (saw-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
              m/HALF_PI))
       (fn ^double [^double t]
         (m// (m/* (double (afn t)) (saw-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
              m/-HALF_PI))))))

(defn- triangle-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/+ 1 b b)
                         b+2 (m/sq b+)]
                     (-> (m/sin (m/* m/TWO_PI (double (m/* b+ t))))
                         (m// b+2)
                         (m/* (if (m/even? b) 1.0 -1.0)))))) 0.0 (range bands)))

(defn ->band-limited-triangle
  ([] (->band-limited-triangle nil))
  ([{:keys [f phase amplitude ^long bands]
     :or {f 1.0 phase 0.0 amplitude 1.0 bands 5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* 0.8105694691387022 ;; 8/pi^2
            (double (afn t))
            (triangle-bl bands (m/+ (double (pfn t)) (double (ffn t)))))))))

;; poly_blep.jsfx-inc, licence: http://www.wtfpl.net/
;; blep/blamp/bluh

(defn blep
  "Band-limited step"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/- (m/sq (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/sq (m/inc (m// (m/dec t) dt)))
    :else 0.0))

(defn blamp
  "Band-limited ramp"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/* -0.3333333333333333 (m/cb (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/* 0.3333333333333333 (m/cb (m/inc (m// (m/dec t) dt))))
    :else 0.0))

(defn bluh
  "Band-limited curve"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (let [t (m/sq (m/dec (m// t dt)))]
                 (m/* -4.0 (m/- (m/sq t) t)))
    (m/> t (m/- 1.0 dt)) (let [t (m/sq (m/inc (m// (m/dec t) dt)))]
                           (m/* 4.0 (m/- (m/sq t) t)))
    :else 0.0))

;;

(defmacro polyblep
  "Creates polyblep signal creator based on generating function f(phase, phase-step) -> amplitude."
  [n fun]
  `(defn ~(symbol (str "polyblep-" n))
     [opts#]
     (let [f# (->const-fn (or (:f opts#) 1.0))
           phase# (->const-fn (or (:phase opts#) 0.0))
           amplitude# (->const-fn (or (:amplitude opts#) 1.0))
           fs# (double (or (:fs opts#) 44100.0))]
       (->> (iterate (fn [^Vec3 stp#]
                       (let [time-in-sec# (m// (.y stp#) fs#)
                             current-phase# (double (phase# time-in-sec#))
                             phase# (m/mod (m/+ (.z stp#) current-phase#) 1.0)
                             current-frequency# (double (f# time-in-sec#))
                             phase-step# (m// current-frequency# fs#)
                             current-amplitude# (double (amplitude# time-in-sec#))]
                         (Vec3. (m/* current-amplitude# (~fun phase# phase-step#))
                                (m/inc (.y stp#))
                                (m/+ (.z stp#) phase-step#)))) (Vec3. 0.0 0.0 0.0))
            (rest)
            (map first)))))
;;

(defn- pb-hyptri
  ^double [^double phase ^double phase-step]
  (m/+ -1.3130352854993313
       (m/* phase-step (m/- (m// (blamp phase phase-step) m/QUARTER_PI)
                            (m/* (blamp (m/frac (m/+ phase 0.5)) phase-step) 9.273)))
       (m/* (m/exp (m/* 4.0 (if (m/< phase 0.5) phase (m/- 1.0 phase))))
            0.3130352854993313)))

(polyblep hyptri2 pb-hyptri)


(defn polyblep-hyptri
  [{:keys [^double fs f phase amplitude]
    :or {fs 44100.0 f 1.0 phase 0.0 amplitude 1.0}}]
  (let [f (->const-fn f)
        phase (->const-fn phase)
        amplitude (->const-fn amplitude)]
    (->> (iterate (fn [^Vec3 stp]
                    (let [time-in-sec (m// (.y stp) fs)
                          current-phase (double (phase time-in-sec))
                          phase (m/mod (m/+ (.z stp) current-phase) 1.0)
                          current-frequency(double (f time-in-sec))
                          phase-step (m// current-frequency fs)
                          current-amplitude (double (amplitude time-in-sec))]
                      (Vec3. (m/* current-amplitude (pb-hyptri phase phase-step))
                             (m/inc (.y stp))
                             (m/+ (.z stp) phase-step)))) (Vec3. 0.0 0.0 0.0))
         (rest)
         (map first))))
