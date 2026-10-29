(ns fastmath.fields.v
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- v-modulate 
  "Modulate"
  ^double [^double amp ^double freq ^double x]
  (m/* amp (m/cos (m/* x freq m/TWO_PI))))

(defn vibration2
  "Vibration 2 http://fractal-resources.deviantart.com/art/Apo-Plugins-Vibration-1-and-2-252001851"
  ([] {:type :regular
       :config (fn [] {:dir (r/drand m/TWO_PI)
                      :angle (r/drand m/TWO_PI)
                      :freq (u/sdrand 0.01 2.0)
                      :amp (u/sdrand 0.1 1.0)
                      :phase (r/drand)
                      :dir2 (r/drand m/TWO_PI)
                      :angle2 (r/drand m/TWO_PI)
                      :freq2 (u/sdrand 0.01 2.0)
                      :amp2 (u/sdrand 0.1 1.0)
                      :phase2 (r/drand)
                      :dm (r/drand -0.5 0.5)
                      :dmfreq (u/sdrand 0.01 1.0)
                      :tm (r/drand -0.5 0.5)
                      :tmfreq (u/sdrand 0.01 1.0)
                      :fm (r/drand -0.5 0.5)
                      :fmfreq (u/sdrand 0.01 1.0)
                      :am (r/drand -0.5 0.5)
                      :amfreq (u/sdrand 0.01 1.0)
                      :d2m (r/drand -0.5 0.5)
                      :d2mfreq (u/sdrand 0.01 1.0)
                      :t2m (r/drand -0.5 0.5)
                      :t2mfreq (u/sdrand 0.01 1.0)
                      :f2m (r/drand -0.5 0.5)
                      :f2mfreq (u/sdrand 0.01 1.0)
                      :a2m (r/drand -0.5 0.5)
                      :a2mfreq (u/sdrand 0.01 1.0)})})
  ([^double amount {:keys [^double dir ^double angle ^double freq ^double amp ^double phase
                           ^double dir2 ^double angle2 ^double freq2 ^double amp2 ^double phase2
                           ^double dm ^double dmfreq
                           ^double tm ^double tmfreq
                           ^double fm ^double fmfreq
                           ^double am ^double amfreq
                           ^double d2m ^double d2mfreq
                           ^double t2m ^double t2mfreq
                           ^double f2m ^double f2mfreq
                           ^double a2m ^double a2mfreq]}]
   (let [cdir (m/cos dir)
         sdir (m/sin dir)
         cdir2 (m/cos dir2)
         sdir2 (m/sin dir2)]
     (fn [^Vec2 v]
       (let [d-along-dir (m/+ (m/* (.x v) cdir)
                            (m/* (.y v) sdir))
             dir-l (m/+ dir (v-modulate dm dmfreq d-along-dir))
             angle-l (m/+ angle (v-modulate tm tmfreq d-along-dir))
             freq-l (m// (v-modulate fm fmfreq d-along-dir) freq)
             amp-l (m/+ amp (m/* amp (v-modulate am amfreq d-along-dir)))
             total-angle (m/+ angle-l dir-l)
             cos-dir (m/cos dir-l)
             sin-dir (m/sin dir-l)
             cos-tot (m/cos total-angle)
             sin-tot (m/sin total-angle)
             scaled-freq (m/* m/TWO_PI freq)
             phase-shift (m// (m/* m/TWO_PI phase) freq)
             d-along-dir (m/+ (m/* (.x v) cos-dir)
                            (m/* (.y v) sin-dir))
             local-amp (m/* amp-l (m/sin (m/+ (m/* d-along-dir scaled-freq) freq-l phase-shift)))
             x (m/+ (.x v) (m/* local-amp cos-tot))
             y (m/+ (.y v) (m/* local-amp sin-tot))

             d-along-dir (m/+ (m/* (.x v) cdir2)
                            (m/* (.y v) sdir2))
             dir-l (m/+ dir2 (v-modulate d2m d2mfreq d-along-dir))
             angle-l (m/+ angle2 (v-modulate t2m t2mfreq d-along-dir))
             freq-l (m// (v-modulate f2m f2mfreq d-along-dir) freq2)
             amp-l (m/+ amp2 (m/* amp2 (v-modulate a2m a2mfreq d-along-dir)))
             total-angle (m/+ angle-l dir-l)
             cos-dir (m/cos dir-l)
             sin-dir (m/sin dir-l)
             cos-tot (m/cos total-angle)
             sin-tot (m/sin total-angle)
             scaled-freq (m/* m/TWO_PI freq2)
             phase-shift (m// (m/* m/TWO_PI phase2) freq2)
             d-along-dir (m/+ (m/* (.x v) cos-dir)
                            (m/* (.y v) sin-dir))
             local-amp (m/* amp-l (m/sin (m/+ (m/* d-along-dir scaled-freq) freq-l phase-shift)))
             x (m/+ x (m/* local-amp cos-tot))
             y (m/+ y (m/* local-amp sin-tot))]
         (Vec2. (m/* amount x)
                (m/* amount y)))))))


(defn vogel
  ([] {:type :random
       :config (fn [] {:scale (u/sdrand 0.01 0.5)
                      :n (u/sdrand 0.01 5.0)})})
  ([^double amount {:keys [^double scale ^double n]}]
   (fn [^Vec2 v]
     (let [i (int (m/inc (r/drand n)))
           a (m/* i 2.399963229728653)
           sina (m/sin a)
           cosa (m/cos a)
           r (m/* amount (m/+ (v/mag v) (m/sqrt i)))]
       (Vec2. (m/* r (m/+ cosa (m/* scale (.x v))))
              (m/* r (m/+ sina (m/* scale (.y v)))))))))


(deftype VoronResType [^double R ^double X0 ^double Y0])
(deftype VoronCalcType [^long M1 ^long N1 ^long k])

(defn voron
  "Voron by eralex61, http://eralex61.deviantart.com/art/Voronoi-Diagram-plugin-153126702"
  ([] {:type :regular
       :config (fn [] {:k (u/sdrand 0.6 1.3)
                      :step (u/sdrand 0.1 1.2)
                      :num (r/drand 0.1 25.0)
                      :xseed (r/irand)
                      :yseed (r/irand)})})
  ([^double amount {:keys [^double k ^double step ^double num ^int xseed ^int yseed]}]
   (fn [^Vec2 v]
     (let [fk (fn ^VoronCalcType [^long M1 ^long N1]
                (VoronCalcType. M1 N1
                                (long (m/inc (m/floor (m/* (r/discrete-noise (m/long-add (m/long-add (m/long-mult M1 19) (m/long-mult N1 257)) xseed) 0) num))))))
           m (long (m/floor (m// (.x v) step)))
           n (long (m/floor (m// (.y v) step)))
           m- (m/dec m)
           m+ (m/inc m)
           n- (m/dec n)
           n+ (m/inc n)
           Ks (mapv fk [m- m- m- m m m m+ m+ m+] [n- n n+ n- n n+ n- n n+])
           ^VoronResType res (reduce (fn [^VoronResType curr ^VoronCalcType calc]
                                       (loop [i (long 0)
                                              ^VoronResType currl curr]
                                         (if (m/< i (.k calc))
                                           (let [X (m/* step (m/+ (.M1 calc)
                                                              (r/discrete-noise (m/long-add
                                                                                 (m/long-add i (m/long-mult 64 (.M1 calc)))
                                                                                 (m/long-add xseed (m/long-mult 15 (.N1 calc)))) 0)))
                                                 Y (m/* step (m/+ (.N1 calc)
                                                              (r/discrete-noise (m/long-add
                                                                                 (m/long-add i (m/long-mult 21 (.M1 calc)))
                                                                                 (m/long-add yseed (m/long-mult 33 (.N1 calc)))) 0)))
                                                 R (m/hypot-sqrt (m/- (.x v) X) (m/- (.y v) Y))]
                                             (recur (unchecked-inc i)
                                                    (if (m/< R (.R currl))
                                                      (VoronResType. R X Y)
                                                      currl)))
                                           currl))) (VoronResType. 20.0 0.0 0.0) Ks)]
       (Vec2. (m/* amount (m/+ (.X0 res) (m/* k (m/- (.x v) (.X0 res)))))
              (m/* amount (m/+ (.Y0 res) (m/* k (m/- (.y v) (.Y0 res))))))))))

;;

(defn vibration
  "Vibration http://fractal-resources.deviantart.com/art/Apo-Plugins-Vibration-1-and-2-252001851"
  ([] {:type :regular
       :config (fn [] {:dir (r/drand m/TWO_PI)
                      :angle (r/drand m/TWO_PI)
                      :freq (u/sdrand 0.01 2.0)
                      :amp (u/sdrand 0.1 1.0)
                      :phase (r/drand)
                      :dir2 (r/drand m/TWO_PI)
                      :angle2 (r/drand m/TWO_PI)
                      :freq2 (u/sdrand 0.01 2.0)
                      :amp2 (u/sdrand 0.1 1.0)
                      :phase2 (r/drand)})})
  ([^double amount {:keys [^double dir ^double angle ^double freq ^double amp ^double phase
                           ^double dir2 ^double angle2 ^double freq2 ^double amp2 ^double phase2]}]
   (let [total-angle (m/+ angle dir)
         cos-dir (m/cos dir)
         sin-dir (m/sin dir)
         cos-tot (m/cos total-angle)
         sin-tot (m/sin total-angle)
         scaled-freq (m/* m/TWO_PI freq)
         phase-shift (m// (m/* m/TWO_PI phase) freq)
         total-angle2 (m/+ angle2 dir2)
         cos-dir2 (m/cos dir2)
         sin-dir2 (m/sin dir2)
         cos-tot2 (m/cos total-angle2)
         sin-tot2 (m/sin total-angle2)
         scaled-freq2 (m/* m/TWO_PI freq2)
         phase-shift2 (m// (m/* m/TWO_PI phase2) freq2)]
     (fn [^Vec2 v]
       (let [d-along-dir (m/+ (m/* (.x v) cos-dir)
                            (m/* (.y v) sin-dir))
             local-amp (m/* amp (m/sin (m/+ (m/* d-along-dir scaled-freq) phase-shift)))
             x (m/+ (.x v) (m/* local-amp cos-tot))
             y (m/+ (.y v) (m/* local-amp sin-tot))
             d-along-dir (m/+ (m/* (.x v) cos-dir2)
                            (m/* (.y v) sin-dir2))
             local-amp (m/* amp2 (m/sin (m/+ (m/* d-along-dir scaled-freq2) phase-shift2)))
             x (m/+ x (m/* local-amp cos-tot2))
             y (m/+ y (m/* local-amp sin-tot2))]
         (Vec2. (m/* amount x) (m/* amount y)))))))

