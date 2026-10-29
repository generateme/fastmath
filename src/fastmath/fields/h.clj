(ns fastmath.fields.h
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn hadamard
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> (cond
           (m/< (r/drand) m/THIRD) (v/mult v 0.5)
           (m/< (r/drand) m/TWO_THIRD) (Vec2. (m/* 0.5 (.y v)) (m/- (m/* -0.5 (.x v)) 0.5))
           :else (Vec2. (m/- (m/* -0.5 (.y v)) 0.5) (m/* 0.5 (.x v))))
         (v/mult amount)))))

(defn- hamid-pline
  ^Vec2 [^Vec4 xy] ;; x1:x y1:y x2:z y2:w
  (let [ydiff (m/- (.w xy) (.y xy))
        xdiff (m/- (.z xy) (.x xy))
        m (if (m/zero? xdiff) 10000.0 (m// ydiff xdiff))
        line-length (m/hypot-sqrt xdiff ydiff)
        d (r/drand line-length)
        xoffset (m// d (m/sqrt (m/inc (m/* m m))))
        xoffset (if (m/< (.z xy) (.x xy)) (m/- xoffset) xoffset)
        yoffset (m/abs (m/* m xoffset))
        yoffset (if (m/< (.w xy) (.y xy)) (m/- yoffset) xoffset)]
    (Vec2. (m/+ (.x xy) xoffset)
           (m/+ (.y xy) yoffset))))

(defn- hamid-pcircle
  ^Vec2 [^Vec3 xy filled]
  (let [radius (if (and filled (m/> (r/drand) 0.05)) (r/drand (.z xy)) (.z xy))
        a (r/drand m/TWO_PI)]
    (Vec2. (m/+ (.x xy) (m/* radius (m/cos a)))
           (m/+ (.y xy) (m/* radius (m/sin a))))))

(defn- hamid-A
  ^Vec2 [^long n]
  (let [n (m/rem n 400.0)]
    (cond
      (m/< n 100) (Vec2. (m// n 100.0) 0.0)
      (m/< n 200) (Vec2. 1.0 (m// (m/- n 100.0) 100.0))
      (m/< n 300) (Vec2. (m// (m/- 300.0 n) 100.0) 1.0)
      :else (Vec2. 0.0 (m// (m/- 400.0 n) 100.0)))))

(defn hamid
  ([] {:type :pattern
       :config (fn [] {:presetid (r/irand 20)
                      :npoints (r/irand 100 5000)
                      :filled (r/brand)
                      :a (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :b (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :c (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :d (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :e (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :f (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :g (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :h (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :i (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :j (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :k (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :l (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :m (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :n (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :o (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))
                      :p (r/randval (u/sirand 1 3) (u/sdrand 0.125 2.0))})})
  ([^double amount {:keys [^long presetid ^long npoints filled
                           ^double a ^double b ^double c ^double d ^double e
                           ^double f ^double g ^double h ^double i ^double j
                           ^double k ^double l ^double m ^double n ^double o ^double p]}]
   (let [npointsd (m// (double npoints))
         presetid (int (m/rem presetid 20))
         hamid-line-1 (let [fa (m/* a 4.0 m/PI) fb (m/* b 2.0 m/PI)
                            fc (m/* c 8.0 m/PI) fd (m/* d 4.0 m/PI)
                            fe (m/* e -0.5) ff (m/* f -0.5)]
                        (fn [^double in]
                          (Vec4. (m/- (m/sin (m/* fa in)))    (m/- (m/cos (m/* fb in)))
                                 (m/* fe (m/sin (m/* fc in))) (m/* ff (m/cos (m/* fd in))))))
         hamid-line-2 (let [fa (m/* a 2.0 m/PI) fb (m/* b 2.0 m/PI)
                            fc (m/* c 8.0 m/PI) fd (m/* d 12.0 m/PI)
                            fe (m/* e -0.5) ff (m/* f -0.5)]
                        (fn [^double in]
                          (Vec4. (m/- (m/sin (m/* fa in)))    (m/- (m/cos (m/* fb in)))
                                 (m/* fe (m/sin (m/* fc in))) (m/* ff (m/cos (m/* fd in))))))
         hamid-line-3 (let [fa (m/* a 8.0 m/PI) fb (m/* b 2.0 m/PI)
                            fc (m/* c 6.0 m/PI) fd (m/* d 2.0 m/PI)
                            fe (m/* e -0.5) ff (m/* f -0.5)]
                        (fn [^double in]
                          (Vec4. (m/- (m/sin (m/* fa in)))    (m/- (m/cos (m/* fb in)))
                                 (m/* fe (m/sin (m/* fc in))) (m/* ff (m/cos (m/* fd in))))))
         hamid-line-4 (let [fa (m/* a 10.0 m/PI) fb (m/* b 2.0 m/PI)
                            fc (m/* c 12.0 m/PI) fd (m/* d 2.0 m/PI)
                            fe (m/* e -0.5) ff (m/* f -0.5)]
                        (fn [^double in]
                          (Vec4. (m/- (m/sin (m/* fa in)))    (m/- (m/cos (m/* fb in)))
                                 (m/* fe (m/sin (m/* fc in))) (m/* ff (m/cos (m/* fd in))))))
         hamid-line-5 (let [fa (m/* a 10.0 m/PI) fb (m/* b 8.0 m/PI)
                            fc (m/* c 12.0 m/PI) fd (m/* d 10.0 m/PI)
                            s (m// 699.0 npoints)]
                        (fn [^double in]
                          (let [in (m/+ in s)]
                            (Vec4. (m/sin (m/* fa in)) (m/cos (m/* fb in))
                                   (m/sin (m/* fc in)) (m/cos (m/* fd in))))))
         hamid-line-6 (let [fa (m/* a -2.0) fb (m/* b 4.0 m/PI) f1 (m// npoints 1000.0)
                            fc (m/* c 0.5) fd (m/* d 6.0 m/PI) fi (m/* i 3.0)
                            fej (m// (m/* e j -2.0) 15.0) ff (m/* f 6.0 m/PI)
                            fgk (m// (m/* g k 4.0) 15.0) fh (m/* h 2.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (m/* fa (m/cos (m/* fb in f1))) (m/* fc (u/spow (m/cos (m/* fd in)) fi))
                                 (m/* fej (m/sin (m/* ff in))) (m/* fgk (m/sin (m/* fh in))))))
         hamid-line-7 (let [fa (m/* a 3.0) fb (m/* b 2.0 m/PI) fh (m/* h 3.0)
                            fc (m/* c 8.0 m/PI) fd (m/* d 1.5) fe (m/* e 2.0 m/PI)
                            fi (m/* i 3.0) ff (m/* f -0.5) fg (m/* g 6.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (m/* fa (u/spow (m/sin (m/* fb in)) fh)) (m/- (m/cos (m/* fc in)))
                                 (m/* fd (u/spow (m/sin (m/* fe in)) fi)) (m/* ff (m/sin (m/* fg in))))))
         hamid-line-8 (let [fa (m/* a 6.0 m/PI) fb (m/* 8.0 m/PI)
                            fc (m/* c 1.5) ff (m/* 3.0 f)
                            fd (m/* d -0.5) fe (m/* e 6.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (m/cos (m/* fa in)) (m/- (m/cos (m/* fb in)))
                                 (m/* fc (u/spow (m/sin (m/* m/TWO_PI in)) ff))
                                 (m/* fd (m/cos (m/* fe in))))))
         hamid-line-9 (let [fa (m/* a 20.0 m/PI) fb (m/* b 18.0 m/PI)
                            fe (m/* e 3.0) fc (m/* c 16.0 m/PI)
                            fd (m/* d 14.0 m/PI) ff (m/* f 3.0)]
                        (fn [^double in]
                          (Vec4. (m/cb (m/sin (m/* fa in))) (u/spow (m/cos (m/* fb in)) fe)
                                 (m/cb (m/sin (m/* fc in))) (u/spow (m/cos (m/* fd in)) ff))))
         hamid-line-10 (fn [^long in]
                         (let [^Vec2 out1 (hamid-A in)
                               ^Vec2 out2 (hamid-A (m/long-add in 199))]
                           (Vec4. (.x out1) (.y out1) (.x out2) (.y out2))))
         hamid-circle-10 (let [fa (m/* a m/TWO_THIRD) fb (m/* b 2.0 m/PI) fc (m/* c 36.0 m/PI)
                               fd (m/* d 30.0 m/PI) fk (m// k 3.0)
                               fl (m/* l 3.0) fm (m// m 3.0)
                               fen (m/* e n m/TWO_THIRD) ff (m/* f 2.0 m/PI) fg (m/* g 36.0 m/PI)
                               fh (m/* h 30.0 m/PI) fo (m// o 3.0) fi (m/* i 14.0 m/PI) fp (m// p 5.0)
                               fj (m/* j 36.0 m/PI)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* fk (m/sin (m/* fc in)) (m/cos (m/* fd in)))
                                       (m/* fm (u/spow (m/sin (m/* 43.982297150257104 in)) fl)))
                                    (m/+ (m/* fen (m/cos (m/* ff in)))
                                       (m/* fo (m/cos (m/* fg in)) (m/cos (m/* fh in)))
                                       (m/* fp (m/cb (m/cos (m/* fi in)))))
                                    (m/+ 0.006666666666666667 (m// (m/sq (m/sq (m/sin (m/* fj in)))) 9.0)))))
         hamid-circle-11 (let [fa (m/* a 14.0 m/PI) fb (m/* b 32.0 m/PI)
                               fc (m/* c 20.0 m/PI) fd (m// d 3.0)]
                           (fn [^double in]
                             (Vec3. (m/sin (m/* fa in)) (m/cos (m/* fb in))
                                    (m/* fd (m/sq (m/sin (m/* fc in)))))))
         hamid-circle-12 (let [fa (m/* a m/TWO_THIRD) fb (m/* b 6.0 m/PI) fc (m/* c 36.0 m/PI)
                               fd (m/* d 30.0 m/PI) fk (m/* k 3.0)
                               fe (m/* e 14.0) fl (m// l 4.0) fm (m/* m 3.0)
                               ff (m/* f m/TWO_THIRD) fg (m/* g 6.0 m/PI)
                               fh (m/* h 28.0 m/PI) fi (m/* i 30.0 m/PI) fp (m// p 20.0)
                               ffn (m// n 5.0) fj (m/* j 36.0 m/PI) fo (m/* o 4.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* 0.25 (m/sin (m/* fc in)) (m/cos (m/* fd in)))
                                       (m/* fl (u/spow (m/sin (m/* fe in)) fk)))
                                    (m/+ (m/* ff (m/cos (m/* fg in)))
                                       (m/* m/THIRD (m/cos (m/* fh in)) (m/cos (m/* fi in)))
                                       (m/* ffn (u/spow (m/cos (m/* 43.982297150257104 in)) fm)))
                                    (m/+  0.006666666666666667
                                        (m/* fp (u/spow (m/sq (m/sin (m/* fj in))) fo))))))
         hamid-circle-13 (let [fa (m/* a 10.0 m/PI) fb (m/* b 28.0 m/PI) fe (m/* e 3.0)
                               fc (m/* c 18.0 m/PI) fd (m// d 3.0)]
                           (fn [^double in]
                             (Vec3. (m/sin (m/* fa in)) (u/spow (m/cos (m/* fb in)) fe)
                                    (m/* fd (m/sq (m/sin (m/* fc in)))))))
         hamid-circle-14 (let [fa (m/* a m/TWO_THIRD) fb (m/* b 18.0 m/PI) fc (m/* c 40.0 m/PI)
                               fd (m/* d 30.0 m/PI) fe (m/* e 34.0 m/PI) fk (m/* k 3.0) fl (m// l 4.0)
                               ff (m/* f m/TWO_THIRD) fg (m/* g 18.0 m/PI) fh (m/* h 28.0 m/PI)
                               f1 125.66370614359172 fi (m/* i 34.0 m/PI) fm (m/* m 3.0) ffn (m// n 3.0)
                               fj (m/* j 56.0 m/PI) fp (m// p 8.0) fo (m/* o 4.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* m/THIRD (m/cb (m/sin (m/* fc in)))
                                          (m/* (m/cos (m/* fd in))))
                                       (m/* fl (u/spow (m/sin (m/* fe in)) fk)))
                                    (m/+ (m/* ff (m/cos (m/* fg in)))
                                       (m/* m/THIRD (m/cb (m/cos (m/* fh in)))
                                          (m/cos (m/* f1 in)))
                                       (m/* ffn (u/spow (m/cos (m/* fi in)) fm)))
                                    (m/+ 0.006666666666666667 (m/* fp (u/spow (m/sin (m/* fj in)) fo))))))
         hamid-circle-15 (let [fa (m/* a m/TWO_THIRD) fb (m/* b 14.0 m/PI) fc (m/* c 10.0 m/PI)
                               fd (m/* d 28.0 m/PI) fe (m/* e 46.0 m/PI) fm (m/* m 3.0) ffn (m// n 3.0)
                               ff (m/* f 14.0 m/PI) fg (m/* g 10.0 m/PI) fh (m/* h 28.0 m/PI)
                               fj (m/* j 34.0 m/PI) fo (m/* o 3.0)
                               fk (m/* k 72.0 m/PI) fl (m/* l 12.0) fp (m// p 20.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* m/THIRD (m/cb (m/sin (m/* fc in)))
                                          (m/* (m/cos (m/* fd in))))
                                       (m/* ffn (u/spow (m/sin (m/* fe in)) fm)))
                                    (m/+ (m/* m/TWO_THIRD (m/cos (m/* ff in)))
                                       (m/* m/THIRD (m/cb (m/cos (m/* fg in)))
                                          (m/cos (m/* fh in)))
                                       (m/* m/THIRD (u/spow (m/cos (m/* fj in)) fo)))
                                    (m/+ 0.006666666666666667
                                       (m/* 0.1111111111111111 (m/sq (m/cb (m/sin (m/* fk in)))))
                                       (m/* fp (m/sq (m/sq (m/sin (m/* fl in)))))))))
         hamid-circle-16 (let [fa (m/* a 38.0 m/PI) fb (m/* b 26.0 m/PI) fc (m/* c 18.0 m/PI)
                               fl (m/* l 3.0) fm (m// m 2.0) fd (m/* d 18.0 m/PI) fe (m/* e 38.0 m/PI)
                               ff (m/* f 36.0 m/PI) fg (m/* g 18.0 m/PI) ffn (m/* n 3.0) fo (m// o 4.0)
                               fh (m/* h 18.0 m/PI) fi (m/* i 100.0 m/PI)
                               fj (m/* j 72.0 m/PI) fk (m/* k 18.0 m/PI) fp (m// p 20.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* 0.5 (m/sin (m/* fa in)) (m/cos (m/* fb in)))
                                       (m/* fm (u/spow (m/sin (m/* fc in)) fl)))
                                    (m/+ (m/* m/TWO_THIRD (m/cos (m/* fd in)))
                                       (m/* m/THIRD (m/cos (m/* fe in)) (m/cos (m/* ff in)))
                                       (m/* fo (u/spow (m/cos (m/* fg in)) ffn)))
                                    (m/+ 0.006666666666666667
                                       (m/* 0.01 (u/spow (m/sin (m/* fh in)) 10.0)
                                          (m/cb (m/sq (m/sin (m/* fi in)))))
                                       (m/* fp (m/sq (m/sq (m/sq (m/sin (m/* fj in)))))
                                          (m/cb (m/sq (m/cos (m/* fk in)))))))))
         hamid-circle-17 (let [fa (m/* a 2.0) fb (m/* b 10.0 m/PI) fc (m/* c 6.0 m/PI)
                               fd (m/* d 28.0 m/PI) fe (m/* e 126.0 m/PI) fl (m// m 10.0)
                               ff (m/* f 6.0 m/PI) fg (m/* g 18.0 m/PI) fh (m/* h 36.0 m/PI)
                               fi (m/* i 126.0 m/PI) fm (m// m 10.0)
                               fj (m/* j 48.0 m/PI) fk (m/* k 12.0 m/PI) ffn (m/* 2.0 n)
                               fo (m// o 18.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* m/THIRD (m/cb (m/sin (m/* fc in))) (m/cos (m/* fd in)))
                                       (m/* fl (m/sin (m/* fe in))))
                                    (m/+ (m/* m/TWO_THIRD (m/cos (m/* ff in)))
                                       (m/* m/THIRD (m/cb (m/cos (m/* fg in))) (m/cos (m/* fh in)))
                                       (m/* fm (m/cos (m/* fi in))))
                                    (m/+ 0.006666666666666667
                                       (m/* fo (m/sq (m/sq (m/sin (m/* fj in))))
                                          (m/inc (u/spow (m/sin (m/* fk in)) ffn)))))))
         hamid-circle-18 (let [fa (m/* a 0.4) fb (m/* b 18.0 m/PI) fc (m/* c 34.0 m/PI)
                               fj (m/* j 3.0) fk (m// k 3.0) fd (m/* d 50.0 m/PI) fl (m/* l 3.0)
                               fm (m/* m 3.0) fe (m/* e 18.0 m/PI) ff (m/* f 34.0 m/PI) ffn (m/* n 3.0)
                               fg (m/* g 50.0 m/PI) fo (m// o 10.0) fi (m/* i 58.0 m/PI) fp (m// p 10.0)]
                           (fn [^double in]
                             (Vec3. (m/+ (m/* fa (m/sin (m/* fb in)))
                                       (m/* fk (u/spow (m/sin (m/* fc in)) fj))
                                       (m/* fm (u/spow (m/sin (m/* fd in)) fl)))
                                    (m/+ (m/* m/TWO_THIRD (m/cos (m/* fe in)))
                                       (m/* m/THIRD (m/cb (m/cos (m/* ff in))))
                                       (m/* fo (u/spow (m/cos (m/* fg in)) ffn)))
                                    (m/+ 0.006666666666666667
                                       (m/* fp (m/sq (m/sq (m/sin (m/* fi in)))))))))
         hamid-circle-19 (let [fa (m/* a 10.0 m/PI) fb (m/* b 16.0 m/PI) fi (m// i 2.0)
                               fd (m/* d 10.0 m/PI) fe (m/* e 34.0 m/PI) fk (m// k 2.0)
                               fh (m/* h 2.0) fj (m/* j 2.0) fg (m/* g 52.0 m/PI)
                               fl (m/* l 4.0) fm (m// m 10.0)]
                           (fn [^double in]
                             (Vec3. (m/* (m/cos (m/* fa in)) (m/- c (m/* fi (u/spow (m/sin (m/* fb in)) fh))))
                                    (m/* (m/sin (m/* fd in)) (m/- f (m/* fk (u/spow (m/cos (m/* fe in)) fj))))
                                    (m/+ 0.005 (m/* fm (u/spow (m/sin (m/* fg in)) fl))))))]
     (fn [_]
       (-> (let [iin (m/inc (r/lrand npoints))
                 in (m/* iin npointsd)]
             (case presetid
               0 (hamid-pline (hamid-line-1 in))
               1 (hamid-pline (hamid-line-2 in))
               2 (hamid-pline (hamid-line-3 in))
               3 (hamid-pline (hamid-line-4 in))
               4 (hamid-pline (hamid-line-5 in))
               5 (hamid-pline (hamid-line-6 in))
               6 (hamid-pline (hamid-line-7 in))
               7 (hamid-pline (hamid-line-8 in))
               8 (hamid-pline (hamid-line-9 in))
               9 (hamid-pline (hamid-line-10 iin))
               10 (hamid-pcircle (hamid-circle-10 in) filled)
               11 (hamid-pcircle (hamid-circle-11 in) filled)
               12 (hamid-pcircle (hamid-circle-12 in) filled)
               13 (hamid-pcircle (hamid-circle-13 in) filled)
               14 (hamid-pcircle (hamid-circle-14 in) filled)
               15 (hamid-pcircle (hamid-circle-15 in) filled)
               16 (hamid-pcircle (hamid-circle-16 in) filled)
               17 (hamid-pcircle (hamid-circle-17 in) filled)
               18 (hamid-pcircle (hamid-circle-18 in) filled)
               19 (hamid-pcircle (hamid-circle-19 in) filled)
               (hamid-pline (hamid-line-3 in))))
           (v/mult amount))))))

(defn handkerchief
  "Handkerchief"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [angle (v/heading v)
           r (v/mag v)]
       (Vec2. (m/* amount (m/* r (m/sin (m/+ angle r))))
              (m/* amount (m/* r (m/cos (m/- angle r)))))))))

(defn heart
  "Heart"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           theta (v/heading v)
           rt (m/* r theta)
           sr (m/sin rt)
           cr (m/cos rt)]
       (Vec2. (m/* amount r sr) (m/- (m/* amount r cr)))))))

(defn heartwf
  ([] {:type :regular
       :config (fn [] (let [scale-r-left (r/drand 0.5 2.0)]
                       {:scale-x (u/sdrand 0.5 2.0)
                        :shift-t (r/drand)
                        :scale-r-left scale-r-left
                        :scale-r-right (r/randval scale-r-left (r/drand 0.5 2.0))}))})
  ([^double amount {:keys [^double scale-x ^double shift-t
                           ^double scale-r-left ^double scale-r-right]}]
   (fn [^Vec2 v]
     (let [a (m/atan2 (.x v) (.y v))
           r (v/mag v)
           t (if (m/neg? a)
               (m/min 60.0 (m/- (m/* -19.098593171027442 a scale-r-left) shift-t))
               (m/min 60.0 (m/- (m/* 19.098593171027442 a scale-r-right) shift-t)))
           tt (m/- (m/* 40.0 t) (m/* t t))
           tr (m/radians t)
           nx (m/* 0.001 (m/+ tt 1200.0) (m/sin tr) r)
           nx (m/* amount scale-x (if (m/neg? a) (m/- nx) nx))
           ny (m/* amount -0.001 (m/+ tt 400.0) (m/cos tr) r)]
       (Vec2. nx ny)))))

(defn hemisphere
  "Hemisphere"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m// amount (m/sqrt (m/inc (v/magsq v))))]
       (Vec2. (m/* r (.x v))
              (m/* r (.y v)))))))

(defn hexmodulus
  ([] {:type :regular
       :config (fn [] {:size (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double size]}]
   (let [hsize (m// m/SQRT3_2 size)
         weight (m// amount m/SQRT3_2)]
     (fn [^Vec2 v]
       (let [^Vec2 XY (v/mult v hsize)
             x (m/- (m/* m/SQRT3_3 (.x XY))
                  (m/* m/THIRD (.y XY)))
             z (m/* m/TWO_THIRD (.y XY))
             y (m/- (m/- x) z)
             rx (m/round x)
             ry (m/round y)
             rz (m/round z)
             x-diff (m/abs (m/- rx x))
             y-diff (m/abs (m/- ry y))
             z-diff (m/abs (m/- rz z))
             ^Vec2 r (cond
                       (and (m/> x-diff y-diff)
                            (m/> x-diff z-diff)) (Vec2. (m/- (m/- ry) rz) rz)
                       (m/> y-diff z-diff) (Vec2. rx rz)
                       :else (Vec2. rx (m/- (m/- rx) ry)))
             fx-h (m/+ (m/* m/SQRT3 (.x r)) (m/* m/SQRT3_2 (.y r)))
             fy-h (m/* 1.5 (.y r))
             fx (m/- (.x XY) fx-h)
             fy (m/- (.y XY) fy-h)]
         (Vec2. (m/* fx weight) (m/* fy weight)))))))

(defn- hexes-cell-centre
  ^Vec2 [^Vec2 in ^double s]
  (Vec2. (m/* s (m/+ (m/* 1.5 (.x in)) (m/* -1.5 (.y in))))
         (m/* s m/SQRT3_2 (m/+ (.x in) (.y in)))))

(defn hexes
  ([] {:type :regular
       :config (fn [] {:cellsize (u/sdrand 0.5 2.0)
                      :power (r/drand -2.0 2.0)
                      :rotate (r/drand -1.0 1.0)
                      :scale (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double cellsize ^double power ^double rotate ^double scale]}]
   (let [a (m/* rotate m/TWO_PI)
         sr (m// cellsize)]
     (fn [^Vec2 U]
       (let [x (m/* m/THIRD (.x U))
             y (m/* m/SQRT3_3 (.y U))
             H (Vec2. (m/floor (m/* sr (m/+ x y))) (m/floor (m/* sr (m/- y x))))
             P (mapv #(hexes-cell-centre (v/add H %) cellsize) u/offsets)
             q (u/closest P 9 U)
             H (v/add H (u/offsets q))
             P (mapv #(hexes-cell-centre (v/add H (u/offsets %)) cellsize) [4 5 8 7 3 0 1])
             P0 (P 0)
             L1 (u/voronoi P 7 0 U)
             Do (v/sub U P0)
             trgL (m/* scale (m/pow (m/+ L1 m/EPSILON) power))
             V (v/rotate Do a)
             L2 (u/voronoi P 7 0 (v/add V P0))
             L (m/max L1 L2)
             R (if (m/< L 0.5)
                 (m// trgL L1)
                 (if (m/> L 0.8)
                   (m// trgL L2)
                   (m// (m/+ (m/* (m// trgL L1) (m/- 0.8 L))
                         (m/* (m// trgL L2) (m/- L 0.5))) 0.3)))]
         (v/mult (v/add (v/mult V R) P0) amount))))))

(defn hole2
  "Hole2"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand 3.0)
                      :d (r/drand -2.0 2.0)
                      :inside (r/brand)
                      :shape (r/irand 10)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d inside shape]}]
   (fn [v]
     (let [rhosq (v/magsq v)
           theta (m/* d (v/heading v))
           delta (m/* c (m/pow (m/inc (m// theta m/PI)) a))
           r (case (unchecked-int shape)
               1 (m/sqrt (m/+ rhosq delta))
               2 (m/sqrt (m/+ rhosq (m/sin (m/* b theta)) delta))
               3 (m/sqrt (m/+ rhosq (m/sin theta) delta))
               4 (m/sqrt (m/- (m/inc (m/+ rhosq (m/sin theta))) delta))
               5 (m/sqrt (m/+ rhosq (m/abs (m/tan theta)) delta))
               6 (m/sqrt (m/+ rhosq (m/inc (m/sin (m/* b theta))) delta))
               7 (m/sqrt (m/+ rhosq (m/abs (m/sin (m/* 0.5 b theta))) delta))
               8 (m/sqrt (m/+ rhosq (m/sin (m/* m/PI (m/sin (m/* b theta)))) delta))
               9 (m/sqrt (m/+ rhosq (m/* 0.5 (m/+ (m/sin (m/* b theta))
                                            (m/sin (m/+ m/M_PI_2 (m/* 2.0 b theta))))) delta))
               (m/+ delta (m/sqrt rhosq)))
           r1 (if inside
                (m// amount r)
                (m/* amount r))]
       (Vec2. (m/* r1 (m/cos theta))
              (m/* r1 (m/sin theta)))))))

(defn hole
  ([] {:type :regular
       :config (fn [] {:inside (r/brand)
                      :a (r/drand -3.0 3.0)})})
  ([^double amount {:keys [inside ^double a]}]
   (fn [^Vec2 v]
     (let [alpha (v/heading v)
           delta (m/pow (m/inc (m// alpha m/PI)) a)
           r (if inside
               (m/* amount (m// delta (m/+ (v/magsq v) delta)))
               (m/* amount (m/sqrt (m/+ (v/magsq v) delta))))]
       (Vec2. (m/* r (m/cos alpha)) (m/* r (m/sin alpha)))))))

(defn holesq
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [^Vec2 v (v/mult v amount)
           fax (m/abs (.x v))
           fay (m/abs (.y v))]
       (if (m/> (m/+ fax fay) 1.0)
         v
         (if (m/> fax fay)
           (let [t (if (m/neg? (.x v))
                     (m/* 0.5 (m/dec (m/+ (.x v) fay)))
                     (m/* 0.5 (m/inc (m/- (.x v) fay))))]
             (Vec2. (m/+ (.x v) t) (m/+ (.y v) (.y v))))
           (let [t (if (m/neg? (.y v))
                     (m/* 0.5 (m/dec (m/+ (.y v) fax)))
                     (m/* 0.5 (m/inc (m/- (.y v) fax))))]
             (Vec2. (m/+ (.x v) (.x v)) (m/+ (.y v) t)))))))))

(defn horseshoe
  "Horseshoe"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/+ m/EPSILON (v/mag v))
           sina (m// (.x v) r)
           cosa (m// (.y v) r)]
       (Vec2. (m/* amount (m/- (m/* sina (.x v)) (m/* cosa (.y v))))
              (m/* amount (m/+ (m/* cosa (.x v)) (m/* sina (.y v)))))))))

(defn hyperbolicellipse
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.05 2.0)})})
  ([^double amount {:keys [^double a]}]
   (fn [^Vec2 v]
     (let [ay (m/* a (.y v))
           ex+ (m/exp (.x v))
           ex- (m/exp (m/- (.x v)))]
       (Vec2. (m/* amount 0.5 (m/- ex+ ex-) (m/cos ay))
              (m/* amount 0.5 (m/+ ex+ ex-) (m/sin ay)))))))

(defn hyperbolic
  "Hyperbolic"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (m/+ m/EPSILON (v/mag v))
           theta (v/heading v)]
       (Vec2. (m/* amount (m// (m/sin theta) r))
              (m/* amount (m/cos theta) r))))))

(defn hypershift2
  ([] {:type :random
       :config (fn [] {:p (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pq (m// m/PI q)
         pp (m// m/PI p)
         spq (m/sin pq)
         spp (m/sin pp)
         shift (m// (m/sin (m/- m/HALF_PI pq pp))
                  (m/sqrt (m/- 1.0 (m/* spq spq) (m/* spp spp))))
         scale2 (m/* (m// (m/sqrt (m/dec (m// (m/sq (m/sin (m/+ m/HALF_PI pp))) (m/sq spq)))))
                   (m/dec (m// (m/sin (m/+ m/HALF_PI pp)) spq)))
         scale (m/- 1.0 (m/* shift shift))]
     (fn [^Vec2 v]
       (let [^Vec2 F (v/mult v scale2)
             rad (m// (v/magsq F))
             x (m/+ (m/* rad (.x v)) shift)
             y (m/* rad (.y v))
             rad (m/* amount (m// scale (m/+ (m/* x x) (m/* y y))))
             angle (m/* pp (m/inc (m/* 2.0 (m/mod (r/lrand) p))))
             X (m/+ (m/* rad x) shift)
             Y (m/* rad y)
             cosa (m/cos angle)
             sina (m/sin angle)]
         (Vec2. (m/- (m/* cosa X) (m/* sina Y))
                (m/+ (m/* sina X) (m/* cosa Y))))))))

(defn hypershift
  "Hypershift"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -2.0 2.0)
                      :stretch (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shift ^double stretch]}]
   (let [scale (m/- 1.0 (m/* shift shift))]
     (fn [^Vec2 v]
       (let [rad (m// (v/magsq v))
             x (m/+ shift (m/* rad (.x v)))
             y (m/* rad (.y v))
             r (m// (m/* amount scale) (m/+ (m/* x x) (m/* y y)))]
         (Vec2. (m/+ shift (m/* r x))
                (m/* r y stretch)))))))

(defn hypertile1
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pa (m// m/TWO_PI p)
         pq (m// m/TWO_PI q)
         r2 (m/- 1.0 (m// (m/dec (m/cos pa))
                      (m/+ (m/cos pa) (m/cos pq))))
         r (if (m/pos? r2) (m// (m/sqrt r2)) 1.0)]
     (fn [^Vec2 v]
       (let [rpa (m/* pa (r/lrand))
             sina (m/sin rpa)
             cosa (m/cos rpa)
             re (m/* r cosa)
             im (m/* r sina)
             a (m/+ (.x v) re)
             b (m/- (.y v) im)
             c (m/inc (m/- (m/* re (.x v)) (m/* im (.y v))))
             d (m/+ (m/* re (.y v)) (m/* im (.x v)))
             vr (m// amount (m/+ (m/* c c) (m/* d d)))]
         (Vec2. (m/* vr (m/+ (m/* a c) (m/* b d)))
                (m/* vr (m/- (m/* b c) (m/* a d)))))))))

(defn hypertile2
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pa (m// m/TWO_PI p)
         pq (m// m/TWO_PI q)
         r2 (m/- 1.0 (m// (m/dec (m/cos pa))
                      (m/+ (m/cos pa) (m/cos pq))))
         r (if (m/pos? r2) (m// (m/sqrt r2)) 1.0)]
     (fn [^Vec2 v]
       (let [a (m/+ (.x v) r)
             b (.y v)
             c (m/inc (m/* r (.x v)))
             d (m/* r (.y v))
             x (m/+ (m/* a c) (m/* b d))
             y (m/- (m/* b c) (m/* a d))
             vr (m// amount (m/+ (m/* c c) (m/* d d)))
             rpa (m/* pa (r/lrand))
             sina (m/sin rpa)
             cosa (m/cos rpa)]
         (Vec2. (m/* vr (m/+ (m/* x cosa) (m/* y sina)))
                (m/* vr (m/- (m/* y cosa) (m/* x sina)))))))))

(defn hypertile
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 8) (u/sdrand 4.0 8.0))
                      :q (r/randval (u/sirand 4 8) (u/sdrand 4.0 8.0))
                      :n (r/randval (r/irand -10 10) (r/drand -10.0 10.0))})})
  ([^double amount {:keys [^double p ^double q ^double n]}]
   (let [pa (m// m/TWO_PI p)
         pq (m// m/TWO_PI q)
         r2 (m/inc (m// (m/- 1.0 (m/cos pa))
                    (m/+ (m/cos pa) (m/cos pq))))
         r (if (m/pos? r2) (m// (m/sqrt r2)) 1.0)
         a (m/* n pa)
         re (m/* r (m/cos a))
         im (m/* r (m/sin a))]
     (fn [^Vec2 v]
       (let [a (m/+ (.x v) re)
             b (m/- (.y v) im)
             c (m/inc (m/- (m/* re (.x v)) (m/* im (.y v))))
             d (m/+ (m/* re (.y v)) (m/* im (.x v)))
             vr (m// amount (m/+ (m/* c c) (m/* d d)))]
         (Vec2. (m/* vr (m/+ (m/* a c) (m/* b d)))
                (m/* vr (m/- (m/* b c) (m/* a d)))))))))

