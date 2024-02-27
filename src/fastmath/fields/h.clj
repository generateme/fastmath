(ns fastmath.fields.h
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn hadamard
  ([] {:type :random})
  ([^double amount _]
   (fn [^Vec2 v]
     (-> (cond
           (< (r/drand) m/THIRD) (v/mult v 0.5)
           (< (r/drand) m/TWO_THIRD) (Vec2. (* 0.5 (.y v)) (- (* -0.5 (.x v)) 0.5))
           :else (Vec2. (- (* -0.5 (.y v)) 0.5) (* 0.5 (.x v))))
         (v/mult amount)))))

(defn- hamid-pline
  ^Vec2 [^Vec4 xy] ;; x1:x y1:y x2:z y2:w
  (let [ydiff (- (.w xy) (.y xy))
        xdiff (- (.z xy) (.x xy))
        m (if (zero? xdiff) 10000.0 (/ ydiff xdiff))
        line-length (m/hypot-sqrt xdiff ydiff)
        d (r/drand line-length)
        xoffset (/ d (m/sqrt (inc (* m m))))
        xoffset (if (< (.z xy) (.x xy)) (- xoffset) xoffset)
        yoffset (m/abs (* m xoffset))
        yoffset (if (< (.w xy) (.y xy)) (- yoffset) xoffset)]
    (Vec2. (+ (.x xy) xoffset)
           (+ (.y xy) yoffset))))

(defn- hamid-pcircle
  ^Vec2 [^Vec3 xy filled]
  (let [radius (if (and filled (> (r/drand) 0.05)) (r/drand (.z xy)) (.z xy))
        a (r/drand m/TWO_PI)]
    (Vec2. (+ (.x xy) (* radius (m/cos a)))
           (+ (.y xy) (* radius (m/sin a))))))

(defn- hamid-A
  ^Vec2 [^long n]
  (let [n (rem n 400.0)]
    (cond
      (< n 100) (Vec2. (/ n 100.0) 0.0)
      (< n 200) (Vec2. 1.0 (/ (- n 100.0) 100.0))
      (< n 300) (Vec2. (/ (- 300.0 n) 100.0) 1.0)
      :else (Vec2. 0.0 (/ (- 400.0 n) 100.0)))))

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
   (let [npointsd (/ (double npoints))
         presetid (int (rem presetid 20))
         hamid-line-1 (let [fa (* a 4.0 m/PI) fb (* b 2.0 m/PI)
                            fc (* c 8.0 m/PI) fd (* d 4.0 m/PI)
                            fe (* e -0.5) ff (* f -0.5)]
                        (fn [^double in]
                          (Vec4. (- (m/sin (* fa in)))    (- (m/cos (* fb in)))
                                 (* fe (m/sin (* fc in))) (* ff (m/cos (* fd in))))))
         hamid-line-2 (let [fa (* a 2.0 m/PI) fb (* b 2.0 m/PI)
                            fc (* c 8.0 m/PI) fd (* d 12.0 m/PI)
                            fe (* e -0.5) ff (* f -0.5)]
                        (fn [^double in]
                          (Vec4. (- (m/sin (* fa in)))    (- (m/cos (* fb in)))
                                 (* fe (m/sin (* fc in))) (* ff (m/cos (* fd in))))))
         hamid-line-3 (let [fa (* a 8.0 m/PI) fb (* b 2.0 m/PI)
                            fc (* c 6.0 m/PI) fd (* d 2.0 m/PI)
                            fe (* e -0.5) ff (* f -0.5)]
                        (fn [^double in]
                          (Vec4. (- (m/sin (* fa in)))    (- (m/cos (* fb in)))
                                 (* fe (m/sin (* fc in))) (* ff (m/cos (* fd in))))))
         hamid-line-4 (let [fa (* a 10.0 m/PI) fb (* b 2.0 m/PI)
                            fc (* c 12.0 m/PI) fd (* d 2.0 m/PI)
                            fe (* e -0.5) ff (* f -0.5)]
                        (fn [^double in]
                          (Vec4. (- (m/sin (* fa in)))    (- (m/cos (* fb in)))
                                 (* fe (m/sin (* fc in))) (* ff (m/cos (* fd in))))))
         hamid-line-5 (let [fa (* a 10.0 m/PI) fb (* b 8.0 m/PI)
                            fc (* c 12.0 m/PI) fd (* d 10.0 m/PI)
                            s (/ 699.0 npoints)]
                        (fn [^double in]
                          (let [in (+ in s)]
                            (Vec4. (m/sin (* fa in)) (m/cos (* fb in))
                                   (m/sin (* fc in)) (m/cos (* fd in))))))
         hamid-line-6 (let [fa (* a -2.0) fb (* b 4.0 m/PI) f1 (/ npoints 1000.0)
                            fc (* c 0.5) fd (* d 6.0 m/PI) fi (* i 3.0)
                            fej (/ (* e j -2.0) 15.0) ff (* f 6.0 m/PI)
                            fgk (/ (* g k 4.0) 15.0) fh (* h 2.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (* fa (m/cos (* fb in f1))) (* fc (u/spow (m/cos (* fd in)) fi))
                                 (* fej (m/sin (* ff in))) (* fgk (m/sin (* fh in))))))
         hamid-line-7 (let [fa (* a 3.0) fb (* b 2.0 m/PI) fh (* h 3.0)
                            fc (* c 8.0 m/PI) fd (* d 1.5) fe (* e 2.0 m/PI)
                            fi (* i 3.0) ff (* f -0.5) fg (* g 6.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (* fa (u/spow (m/sin (* fb in)) fh)) (- (m/cos (* fc in)))
                                 (* fd (u/spow (m/sin (* fe in)) fi)) (* ff (m/sin (* fg in))))))
         hamid-line-8 (let [fa (* a 6.0 m/PI) fb (* 8.0 m/PI)
                            fc (* c 1.5) ff (* 3.0 f)
                            fd (* d -0.5) fe (* e 6.0 m/PI)]
                        (fn [^double in]
                          (Vec4. (m/cos (* fa in)) (- (m/cos (* fb in)))
                                 (* fc (u/spow (m/sin (* m/TWO_PI in)) ff))
                                 (* fd (m/cos (* fe in))))))
         hamid-line-9 (let [fa (* a 20.0 m/PI) fb (* b 18.0 m/PI)
                            fe (* e 3.0) fc (* c 16.0 m/PI)
                            fd (* d 14.0 m/PI) ff (* f 3.0)]
                        (fn [^double in]
                          (Vec4. (m/cb (m/sin (* fa in))) (u/spow (m/cos (* fb in)) fe)
                                 (m/cb (m/sin (* fc in))) (u/spow (m/cos (* fd in)) ff))))
         hamid-line-10 (fn [^long in]
                         (let [^Vec2 out1 (hamid-A in)
                               ^Vec2 out2 (hamid-A (+ in 199))]
                           (Vec4. (.x out1) (.y out1) (.x out2) (.y out2))))
         hamid-circle-10 (let [fa (* a m/TWO_THIRD) fb (* b 2.0 m/PI) fc (* c 36.0 m/PI)
                               fd (* d 30.0 m/PI) fk (/ k 3.0)
                               fl (* l 3.0) fm (/ m 3.0)
                               fen (* e n m/TWO_THIRD) ff (* f 2.0 m/PI) fg (* g 36.0 m/PI)
                               fh (* h 30.0 m/PI) fo (/ o 3.0) fi (* i 14.0 m/PI) fp (/ p 5.0)
                               fj (* j 36.0 m/PI)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* fk (m/sin (* fc in)) (m/cos (* fd in)))
                                       (* fm (u/spow (m/sin (* 43.982297150257104 in)) fl)))
                                    (+ (* fen (m/cos (* ff in)))
                                       (* fo (m/cos (* fg in)) (m/cos (* fh in)))
                                       (* fp (m/cb (m/cos (* fi in)))))
                                    (+ 0.006666666666666667 (/ (m/sq (m/sq (m/sin (* fj in)))) 9.0)))))
         hamid-circle-11 (let [fa (* a 14.0 m/PI) fb (* b 32.0 m/PI)
                               fc (* c 20.0 m/PI) fd (/ d 3.0)]
                           (fn [^double in]
                             (Vec3. (m/sin (* fa in)) (m/cos (* fb in))
                                    (* fd (m/sq (m/sin (* fc in)))))))
         hamid-circle-12 (let [fa (* a m/TWO_THIRD) fb (* b 6.0 m/PI) fc (* c 36.0 m/PI)
                               fd (* d 30.0 m/PI) fk (* k 3.0)
                               fe (* e 14.0) fl (/ l 4.0) fm (* m 3.0)
                               ff (* f m/TWO_THIRD) fg (* g 6.0 m/PI)
                               fh (* h 28.0 m/PI) fi (* i 30.0 m/PI) fp (/ p 20.0)
                               ffn (/ n 5.0) fj (* j 36.0 m/PI) fo (* o 4.0)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* 0.25 (m/sin (* fc in)) (m/cos (* fd in)))
                                       (* fl (u/spow (m/sin (* fe in)) fk)))
                                    (+ (* ff (m/cos (* fg in)))
                                       (* m/THIRD (m/cos (* fh in)) (m/cos (* fi in)))
                                       (* ffn (u/spow (m/cos (* 43.982297150257104 in)) fm)))
                                    (+  0.006666666666666667
                                        (* fp (u/spow (m/sq (m/sin (* fj in))) fo))))))
         hamid-circle-13 (let [fa (* a 10.0 m/PI) fb (* b 28.0 m/PI) fe (* e 3.0)
                               fc (* c 18.0 m/PI) fd (/ d 3.0)]
                           (fn [^double in]
                             (Vec3. (m/sin (* fa in)) (u/spow (m/cos (* fb in)) fe)
                                    (* fd (m/sq (m/sin (* fc in)))))))
         hamid-circle-14 (let [fa (* a m/TWO_THIRD) fb (* b 18.0 m/PI) fc (* c 40.0 m/PI)
                               fd (* d 30.0 m/PI) fe (* e 34.0 m/PI) fk (* k 3.0) fl (/ l 4.0)
                               ff (* f m/TWO_THIRD) fg (* g 18.0 m/PI) fh (* h 28.0 m/PI)
                               f1 125.66370614359172 fi (* i 34.0 m/PI) fm (* m 3.0) ffn (/ n 3.0)
                               fj (* j 56.0 m/PI) fp (/ p 8.0) fo (* o 4.0)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* m/THIRD (m/cb (m/sin (* fc in)))
                                          (* (m/cos (* fd in))))
                                       (* fl (u/spow (m/sin (* fe in)) fk)))
                                    (+ (* ff (m/cos (* fg in)))
                                       (* m/THIRD (m/cb (m/cos (* fh in)))
                                          (m/cos (* f1 in)))
                                       (* ffn (u/spow (m/cos (* fi in)) fm)))
                                    (+ 0.006666666666666667 (* fp (u/spow (m/sin (* fj in)) fo))))))
         hamid-circle-15 (let [fa (* a m/TWO_THIRD) fb (* b 14.0 m/PI) fc (* c 10.0 m/PI)
                               fd (* d 28.0 m/PI) fe (* e 46.0 m/PI) fm (* m 3.0) ffn (/ n 3.0)
                               ff (* f 14.0 m/PI) fg (* g 10.0 m/PI) fh (* h 28.0 m/PI)
                               fj (* j 34.0 m/PI) fo (* o 3.0)
                               fk (* k 72.0 m/PI) fl (* l 12.0) fp (/ p 20.0)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* m/THIRD (m/cb (m/sin (* fc in)))
                                          (* (m/cos (* fd in))))
                                       (* ffn (u/spow (m/sin (* fe in)) fm)))
                                    (+ (* m/TWO_THIRD (m/cos (* ff in)))
                                       (* m/THIRD (m/cb (m/cos (* fg in)))
                                          (m/cos (* fh in)))
                                       (* m/THIRD (u/spow (m/cos (* fj in)) fo)))
                                    (+ 0.006666666666666667
                                       (* 0.1111111111111111 (m/sq (m/cb (m/sin (* fk in)))))
                                       (* fp (m/sq (m/sq (m/sin (* fl in)))))))))
         hamid-circle-16 (let [fa (* a 38.0 m/PI) fb (* b 26.0 m/PI) fc (* c 18.0 m/PI)
                               fl (* l 3.0) fm (/ m 2.0) fd (* d 18.0 m/PI) fe (* e 38.0 m/PI)
                               ff (* f 36.0 m/PI) fg (* g 18.0 m/PI) ffn (* n 3.0) fo (/ o 4.0)
                               fh (* h 18.0 m/PI) fi (* i 100.0 m/PI)
                               fj (* j 72.0 m/PI) fk (* k 18.0 m/PI) fp (/ p 20.0)]
                           (fn [^double in]
                             (Vec3. (+ (* 0.5 (m/sin (* fa in)) (m/cos (* fb in)))
                                       (* fm (u/spow (m/sin (* fc in)) fl)))
                                    (+ (* m/TWO_THIRD (m/cos (* fd in)))
                                       (* m/THIRD (m/cos (* fe in)) (m/cos (* ff in)))
                                       (* fo (u/spow (m/cos (* fg in)) ffn)))
                                    (+ 0.006666666666666667
                                       (* 0.01 (u/spow (m/sin (* fh in)) 10.0)
                                          (m/cb (m/sq (m/sin (* fi in)))))
                                       (* fp (m/sq (m/sq (m/sq (m/sin (* fj in)))))
                                          (m/cb (m/sq (m/cos (* fk in)))))))))
         hamid-circle-17 (let [fa (* a 2.0) fb (* b 10.0 m/PI) fc (* c 6.0 m/PI)
                               fd (* d 28.0 m/PI) fe (* e 126.0 m/PI) fl (/ m 10.0)
                               ff (* f 6.0 m/PI) fg (* g 18.0 m/PI) fh (* h 36.0 m/PI)
                               fi (* i 126.0 m/PI) fm (/ m 10.0)
                               fj (* j 48.0 m/PI) fk (* k 12.0 m/PI) ffn (* 2.0 n)
                               fo (/ o 18.0)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* m/THIRD (m/cb (m/sin (* fc in))) (m/cos (* fd in)))
                                       (* fl (m/sin (* fe in))))
                                    (+ (* m/TWO_THIRD (m/cos (* ff in)))
                                       (* m/THIRD (m/cb (m/cos (* fg in))) (m/cos (* fh in)))
                                       (* fm (m/cos (* fi in))))
                                    (+ 0.006666666666666667
                                       (* fo (m/sq (m/sq (m/sin (* fj in))))
                                          (inc (u/spow (m/sin (* fk in)) ffn)))))))
         hamid-circle-18 (let [fa (* a 0.4) fb (* b 18.0 m/PI) fc (* c 34.0 m/PI)
                               fj (* j 3.0) fk (/ k 3.0) fd (* d 50.0 m/PI) fl (* l 3.0)
                               fm (* m 3.0) fe (* e 18.0 m/PI) ff (* f 34.0 m/PI) ffn (* n 3.0)
                               fg (* g 50.0 m/PI) fo (/ o 10.0) fi (* i 58.0 m/PI) fp (/ p 10.0)]
                           (fn [^double in]
                             (Vec3. (+ (* fa (m/sin (* fb in)))
                                       (* fk (u/spow (m/sin (* fc in)) fj))
                                       (* fm (u/spow (m/sin (* fd in)) fl)))
                                    (+ (* m/TWO_THIRD (m/cos (* fe in)))
                                       (* m/THIRD (m/cb (m/cos (* ff in))))
                                       (* fo (u/spow (m/cos (* fg in)) ffn)))
                                    (+ 0.006666666666666667
                                       (* fp (m/sq (m/sq (m/sin (* fi in)))))))))
         hamid-circle-19 (let [fa (* a 10.0 m/PI) fb (* b 16.0 m/PI) fi (/ i 2.0)
                               fd (* d 10.0 m/PI) fe (* e 34.0 m/PI) fk (/ k 2.0)
                               fh (* h 2.0) fj (* j 2.0) fg (* g 52.0 m/PI)
                               fl (* l 4.0) fm (/ m 10.0)]
                           (fn [^double in]
                             (Vec3. (* (m/cos (* fa in)) (- c (* fi (u/spow (m/sin (* fb in)) fh))))
                                    (* (m/sin (* fd in)) (- f (* fk (u/spow (m/cos (* fe in)) fj))))
                                    (+ 0.005 (* fm (u/spow (m/sin (* fg in)) fl))))))]
     (fn [_]
       (-> (let [iin (inc (r/lrand npoints))
                 in (* iin npointsd)]
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
       (Vec2. (* amount (* r (m/sin (+ angle r))))
              (* amount (* r (m/cos (- angle r)))))))))

(defn heart
  "Heart"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (v/mag v)
           theta (v/heading v)
           rt (* r theta)
           sr (m/sin rt)
           cr (m/cos rt)]
       (Vec2. (* amount r sr) (- (* amount r cr)))))))

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
           t (if (neg? a)
               (min 60.0 (- (* -19.098593171027442 a scale-r-left) shift-t))
               (min 60.0 (- (* 19.098593171027442 a scale-r-right) shift-t)))
           tt (- (* 40.0 t) (* t t))
           tr (m/radians t)
           nx (* 0.001 (+ tt 1200.0) (m/sin tr) r)
           nx (* amount scale-x (if (neg? a) (- nx) nx))
           ny (* amount -0.001 (+ tt 400.0) (m/cos tr) r)]
       (Vec2. nx ny)))))

(defn hemisphere
  "Hemisphere"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (/ amount (m/sqrt (inc (v/magsq v))))]
       (Vec2. (* r (.x v))
              (* r (.y v)))))))

(defn hexmodulus
  ([] {:type :regular
       :config (fn [] {:size (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double size]}]
   (let [hsize (/ m/SQRT3_2 size)
         weight (/ amount m/SQRT3_2)]
     (fn [^Vec2 v]
       (let [^Vec2 XY (v/mult v hsize)
             x (- (* m/SQRT3_3 (.x XY))
                  (* m/THIRD (.y XY)))
             z (* m/TWO_THIRD (.y XY))
             y (- (- x) z)
             rx (m/round x)
             ry (m/round y)
             rz (m/round z)
             x-diff (m/abs (- rx x))
             y-diff (m/abs (- ry y))
             z-diff (m/abs (- rz z))
             ^Vec2 r (cond
                       (and (> x-diff y-diff)
                            (> x-diff z-diff)) (Vec2. (- (- ry) rz) rz)
                       (> y-diff z-diff) (Vec2. rx rz)
                       :else (Vec2. rx (- (- rx) ry)))
             fx-h (+ (* m/SQRT3 (.x r)) (* m/SQRT3_2 (.y r)))
             fy-h (* 1.5 (.y r))
             fx (- (.x XY) fx-h)
             fy (- (.y XY) fy-h)]
         (Vec2. (* fx weight) (* fy weight)))))))

(defn- hexes-cell-centre
  ^Vec2 [^Vec2 in ^double s]
  (Vec2. (* s (+ (* 1.5 (.x in)) (* -1.5 (.y in))))
         (* s m/SQRT3_2 (+ (.x in) (.y in)))))

(defn hexes
  ([] {:type :regular
       :config (fn [] {:cellsize (u/sdrand 0.5 2.0)
                      :power (r/drand -2.0 2.0)
                      :rotate (r/drand -1.0 1.0)
                      :scale (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double cellsize ^double power ^double rotate ^double scale]}]
   (let [a (* rotate m/TWO_PI)
         sr (/ cellsize)]
     (fn [^Vec2 U]
       (let [x (* m/THIRD (.x U))
             y (* m/SQRT3_3 (.y U))
             H (Vec2. (m/floor (* sr (+ x y))) (m/floor (* sr (- y x))))
             P (mapv #(hexes-cell-centre (v/add H %) cellsize) u/offsets)
             q (u/closest P 9 U)
             H (v/add H (u/offsets q))
             P (mapv #(hexes-cell-centre (v/add H (u/offsets %)) cellsize) [4 5 8 7 3 0 1])
             P0 (P 0)
             L1 (u/voronoi P 7 0 U)
             Do (v/sub U P0)
             trgL (* scale (m/pow (+ L1 m/EPSILON) power))
             V (v/rotate Do a)
             L2 (u/voronoi P 7 0 (v/add V P0))
             L (max L1 L2)
             R (if (< L 0.5)
                 (/ trgL L1)
                 (if (> L 0.8)
                   (/ trgL L2)
                   (/ (+ (* (/ trgL L1) (- 0.8 L))
                         (* (/ trgL L2) (- L 0.5))) 0.3)))]
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
           theta (* d (v/heading v))
           delta (* c (m/pow (inc (/ theta m/PI)) a))
           r (case (unchecked-int shape)
               1 (m/sqrt (+ rhosq delta))
               2 (m/sqrt (+ rhosq (m/sin (* b theta)) delta))
               3 (m/sqrt (+ rhosq (m/sin theta) delta))
               4 (m/sqrt (- (inc (+ rhosq (m/sin theta))) delta))
               5 (m/sqrt (+ rhosq (m/abs (m/tan theta)) delta))
               6 (m/sqrt (+ rhosq (inc (m/sin (* b theta))) delta))
               7 (m/sqrt (+ rhosq (m/abs (m/sin (* 0.5 b theta))) delta))
               8 (m/sqrt (+ rhosq (m/sin (* m/PI (m/sin (* b theta)))) delta))
               9 (m/sqrt (+ rhosq (* 0.5 (+ (m/sin (* b theta))
                                            (m/sin (+ m/M_PI_2 (* 2.0 b theta))))) delta))
               (+ delta (m/sqrt rhosq)))
           r1 (if inside
                (/ amount r)
                (* amount r))]
       (Vec2. (* r1 (m/cos theta))
              (* r1 (m/sin theta)))))))

(defn hole
  ([] {:type :regular
       :config (fn [] {:inside (r/brand)
                      :a (r/drand -3.0 3.0)})})
  ([^double amount {:keys [inside ^double a]}]
   (fn [^Vec2 v]
     (let [alpha (v/heading v)
           delta (m/pow (inc (/ alpha m/PI)) a)
           r (if inside
               (* amount (/ delta (+ (v/magsq v) delta)))
               (* amount (m/sqrt (+ (v/magsq v) delta))))]
       (Vec2. (* r (m/cos alpha)) (* r (m/sin alpha)))))))

(defn holesq
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [^Vec2 v (v/mult v amount)
           fax (m/abs (.x v))
           fay (m/abs (.y v))]
       (if (> (+ fax fay) 1.0)
         v
         (if (> fax fay)
           (let [t (if (neg? (.x v))
                     (* 0.5 (dec (+ (.x v) fay)))
                     (* 0.5 (inc (- (.x v) fay))))]
             (Vec2. (+ (.x v) t) (+ (.y v) (.y v))))
           (let [t (if (neg? (.y v))
                     (* 0.5 (dec (+ (.y v) fax)))
                     (* 0.5 (inc (- (.y v) fax))))]
             (Vec2. (+ (.x v) (.x v)) (+ (.y v) t)))))))))

(defn horseshoe
  "Horseshoe"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (+ m/EPSILON (v/mag v))
           sina (/ (.x v) r)
           cosa (/ (.y v) r)]
       (Vec2. (* amount (- (* sina (.x v)) (* cosa (.y v))))
              (* amount (+ (* cosa (.x v)) (* sina (.y v)))))))))

(defn hyperbolicellipse
  ([] {:type :regular
       :config (fn [] {:a (u/sdrand 0.05 2.0)})})
  ([^double amount {:keys [^double a]}]
   (fn [^Vec2 v]
     (let [ay (* a (.y v))
           ex+ (m/exp (.x v))
           ex- (m/exp (- (.x v)))]
       (Vec2. (* amount 0.5 (- ex+ ex-) (m/cos ay))
              (* amount 0.5 (+ ex+ ex-) (m/sin ay)))))))

(defn hyperbolic
  "Hyperbolic"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [r (+ m/EPSILON (v/mag v))
           theta (v/heading v)]
       (Vec2. (* amount (/ (m/sin theta) r))
              (* amount (m/cos theta) r))))))

(defn hypershift2
  ([] {:type :random
       :config (fn [] {:p (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pq (/ m/PI q)
         pp (/ m/PI p)
         spq (m/sin pq)
         spp (m/sin pp)
         shift (/ (m/sin (- m/HALF_PI pq pp))
                  (m/sqrt (- 1.0 (* spq spq) (* spp spp))))
         scale2 (* (/ (m/sqrt (dec (/ (m/sq (m/sin (+ m/HALF_PI pp))) (m/sq spq)))))
                   (dec (/ (m/sin (+ m/HALF_PI pp)) spq)))
         scale (- 1.0 (* shift shift))]
     (fn [^Vec2 v]
       (let [^Vec2 F (v/mult v scale2)
             rad (/ (v/magsq F))
             x (+ (* rad (.x v)) shift)
             y (* rad (.y v))
             rad (* amount (/ scale (+ (* x x) (* y y))))
             angle (* pp (inc (* 2.0 (mod (r/lrand) p))))
             X (+ (* rad x) shift)
             Y (* rad y)
             cosa (m/cos angle)
             sina (m/sin angle)]
         (Vec2. (- (* cosa X) (* sina Y))
                (+ (* sina X) (* cosa Y))))))))

(defn hypershift
  "Hypershift"
  ([] {:type :regular
       :config (fn [] {:shift (r/drand -2.0 2.0)
                      :stretch (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double shift ^double stretch]}]
   (let [scale (- 1.0 (* shift shift))]
     (fn [^Vec2 v]
       (let [rad (/ (v/magsq v))
             x (+ shift (* rad (.x v)))
             y (* rad (.y v))
             r (/ (* amount scale) (+ (* x x) (* y y)))]
         (Vec2. (+ shift (* r x))
                (* r y stretch)))))))

(defn hypertile1
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pa (/ m/TWO_PI p)
         pq (/ m/TWO_PI q)
         r2 (- 1.0 (/ (dec (m/cos pa))
                      (+ (m/cos pa) (m/cos pq))))
         r (if (pos? r2) (/ (m/sqrt r2)) 1.0)]
     (fn [^Vec2 v]
       (let [rpa (* pa (r/lrand))
             sina (m/sin rpa)
             cosa (m/cos rpa)
             re (* r cosa)
             im (* r sina)
             a (+ (.x v) re)
             b (- (.y v) im)
             c (inc (- (* re (.x v)) (* im (.y v))))
             d (+ (* re (.y v)) (* im (.x v)))
             vr (/ amount (+ (* c c) (* d d)))]
         (Vec2. (* vr (+ (* a c) (* b d)))
                (* vr (- (* b c) (* a d)))))))))

(defn hypertile2
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 15) (u/sdrand 4.0 15.0))
                      :q (r/randval (u/sirand 4 15) (u/sdrand 4.0 15.0))})})
  ([^double amount {:keys [^double p ^double q]}]
   (let [pa (/ m/TWO_PI p)
         pq (/ m/TWO_PI q)
         r2 (- 1.0 (/ (dec (m/cos pa))
                      (+ (m/cos pa) (m/cos pq))))
         r (if (pos? r2) (/ (m/sqrt r2)) 1.0)]
     (fn [^Vec2 v]
       (let [a (+ (.x v) r)
             b (.y v)
             c (inc (* r (.x v)))
             d (* r (.y v))
             x (+ (* a c) (* b d))
             y (- (* b c) (* a d))
             vr (/ amount (+ (* c c) (* d d)))
             rpa (* pa (r/lrand))
             sina (m/sin rpa)
             cosa (m/cos rpa)]
         (Vec2. (* vr (+ (* x cosa) (* y sina)))
                (* vr (- (* y cosa) (* x sina)))))))))

(defn hypertile
  ([] {:type :random
       :config (fn [] {:p (r/randval 0.9 (u/sirand 4 8) (u/sdrand 4.0 8.0))
                      :q (r/randval (u/sirand 4 8) (u/sdrand 4.0 8.0))
                      :n (r/randval (r/irand -10 10) (r/drand -10.0 10.0))})})
  ([^double amount {:keys [^double p ^double q ^double n]}]
   (let [pa (/ m/TWO_PI p)
         pq (/ m/TWO_PI q)
         r2 (inc (/ (- 1.0 (m/cos pa))
                    (+ (m/cos pa) (m/cos pq))))
         r (if (pos? r2) (/ (m/sqrt r2)) 1.0)
         a (* n pa)
         re (* r (m/cos a))
         im (* r (m/sin a))]
     (fn [^Vec2 v]
       (let [a (+ (.x v) re)
             b (- (.y v) im)
             c (inc (- (* re (.x v)) (* im (.y v))))
             d (+ (* re (.y v)) (* im (.x v)))
             vr (/ amount (+ (* c c) (* d d)))]
         (Vec2. (* vr (+ (* a c) (* b d)))
                (* vr (- (* b c) (* a d)))))))))

(m/unuse-primitive-operators)
