(ns fastmath.signal.waveform
  (:require [fastmath.core :as m]))

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

