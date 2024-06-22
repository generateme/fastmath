(ns fastmath.fields.g
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- fclp ^double [^double a] (if (neg? a) (- (rem (m/abs a) 1.0)) (rem (m/abs a) 1.0)))
(defn- fscl ^double [^double a] (fclp (* 0.5 (inc a))))
(defn- fosc ^double [^double p ^double a] (fscl (- (m/cos (* p a m/TWO_PI)))))
(defn- flip ^double [^double a ^double b ^double c] (+ a (* c (- b a))))

(defn gdoffs
  "GDOffs"
  ([] {:type :regular
       :config (fn [] {:delta-x (r/drand -2.5 2.5)
                      :delta-y (r/drand -2.5 2.5)
                      :area-x (u/sdrand 0.1 2.5)
                      :area-y (u/sdrand 0.1 2.5)
                      :center-x (r/drand -1.5 1.5)
                      :center-y (r/drand -1.5 1.5)
                      :gamma (r/drand -2.5 2.5)
                      :square (r/brand)})})
  ([^double amount {:keys [^double delta-x ^double delta-y ^double area-x ^double area-y
                           ^double center-x ^double center-y ^double gamma square]}]
   (let [gdodx (* delta-x 0.1)
         gdody (* delta-y 0.1)
         gdoax (* 2.0 (if (< (m/abs area-x) 0.1) 0.1 (m/abs area-x)))
         gdoay (* 2.0 (if (< (m/abs area-y) 0.1) 0.1 (m/abs area-y)))
         gdocx center-x
         gdocy center-y
         gdog gamma
         gdos square
         gdob (/ (* gdog 2.0) (max gdoax gdoay))]
     (fn [^Vec2 v]
       (let [osc-x (fosc gdodx 1.0)
             osc-y (if gdos osc-x (fosc gdody 1.0))
             in-x (+ (.x v) gdocx)
             in-y (+ (.y v) gdocy)]
         (-> (Vec2. (flip (flip in-x (fosc in-x 4.0) osc-x) (fosc (fclp (* gdob in-x)) 4.0) osc-x)
                    (flip (flip in-y (fosc in-y 4.0) osc-y) (fosc (fclp (* gdob in-y)) 4.0) osc-y))
             (v/mult amount)))))))

(defn gamma
  "gamma by zephyrtronium, http://fractal-resources.deviantart.com/art/Gamma-Apophysis-Plugin-154808483"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (* amount (special/log-gamma (v/mag v)))
            (* amount (v/heading v))))))

(defn gaussianblur
  "Gaussian"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [a (r/drand m/TWO_PI)
           r (* amount (+ (r/drand) (r/drand) (r/drand) (r/drand) -2.0))]
       (Vec2. (* r (m/cos a)) (* r (m/sin a)))))))

(def ^:private zerov3 (Vec3. 0.0 0.0 0.0))

(defn glitchy2
  ([] {:type :regular
       :config (fn [] {:mode (r/brand)
                      :lr-spin (r/drand -3.0 3.0)
                      :lr-x (r/drand -1.0 1.0)
                      :lr-shift-x (r/drand -1.0 1.0)
                      :lr-y (r/drand -1.0 1.0)
                      :lr-shift-y (r/drand -1.0 1.0)
                      :lr-re-a (r/drand -1.0 1.0)
                      :lr-re-b (r/drand -1.0 1.0)
                      :lr-re-c (r/drand -1.0 1.0)
                      :lr-re-d (r/drand -1.0 1.0)
                      :lr-im-a (r/drand -1.0 1.0)
                      :lr-im-b (r/drand -1.0 1.0)
                      :lr-im-c (r/drand -1.0 1.0)
                      :lr-im-d (r/drand -1.0 1.0)
                      :ur-spin (r/drand -3.0 3.0)
                      :ur-x (r/drand -1.0 1.0)
                      :ur-shift-x (r/drand -1.0 1.0)
                      :ur-y (r/drand -1.0 1.0)
                      :ur-shift-y (r/drand -1.0 1.0)
                      :ur-re-a (r/drand -1.0 1.0)
                      :ur-re-b (r/drand -1.0 1.0)
                      :ur-re-c (r/drand -1.0 1.0)
                      :ur-re-d (r/drand -1.0 1.0)
                      :ur-im-a (r/drand -1.0 1.0)
                      :ur-im-b (r/drand -1.0 1.0)
                      :ur-im-c (r/drand -1.0 1.0)
                      :ur-im-d (r/drand -1.0 1.0)
                      :ll-spin (r/drand -3.0 3.0)
                      :ll-x (r/drand -1.0 1.0)
                      :ll-shift-x (r/drand -1.0 1.0)
                      :ll-y (r/drand -1.0 1.0)
                      :ll-shift-y (r/drand -1.0 1.0)
                      :ll-re-a (r/drand -1.0 1.0)
                      :ll-re-b (r/drand -1.0 1.0)
                      :ll-re-c (r/drand -1.0 1.0)
                      :ll-re-d (r/drand -1.0 1.0)
                      :ll-im-a (r/drand -1.0 1.0)
                      :ll-im-b (r/drand -1.0 1.0)
                      :ll-im-c (r/drand -1.0 1.0)
                      :ll-im-d (r/drand -1.0 1.0)
                      :ul-spin (r/drand -3.0 3.0)
                      :ul-x (r/drand -1.0 1.0)
                      :ul-shift-x (r/drand -1.0 1.0)
                      :ul-y (r/drand -1.0 1.0)
                      :ul-shift-y (r/drand -1.0 1.0)
                      :ul-re-a (r/drand -1.0 1.0)
                      :ul-re-b (r/drand -1.0 1.0)
                      :ul-re-c (r/drand -1.0 1.0)
                      :ul-re-d (r/drand -1.0 1.0)
                      :ul-im-a (r/drand -1.0 1.0)
                      :ul-im-b (r/drand -1.0 1.0)
                      :ul-im-c (r/drand -1.0 1.0)
                      :ul-im-d (r/drand -1.0 1.0)})})
  ([^double amount {:keys [mode
                           ^double lr-spin ^double lr-x ^double lr-shift-x ^double lr-y ^double lr-shift-y
                           ^double lr-re-a ^double lr-re-b ^double lr-re-c ^double lr-re-d
                           ^double lr-im-a ^double lr-im-b ^double lr-im-c ^double lr-im-d
                           ^double ur-spin ^double ur-x ^double ur-shift-x ^double ur-y ^double ur-shift-y
                           ^double ur-re-a ^double ur-re-b ^double ur-re-c ^double ur-re-d
                           ^double ur-im-a ^double ur-im-b ^double ur-im-c ^double ur-im-d
                           ^double ll-spin ^double ll-x ^double ll-shift-x ^double ll-y ^double ll-shift-y
                           ^double ll-re-a ^double ll-re-b ^double ll-re-c ^double ll-re-d
                           ^double ll-im-a ^double ll-im-b ^double ll-im-c ^double ll-im-d
                           ^double ul-spin ^double ul-x ^double ul-shift-x ^double ul-y ^double ul-shift-y
                           ^double ul-re-a ^double ul-re-b ^double ul-re-c ^double ul-re-d
                           ^double ul-im-a ^double ul-im-b ^double ul-im-c ^double ul-im-d]}]
   (let [lr-pz-sin (m/sin (* lr-spin m/M_PI_2))
         lr-pz-cos (m/sin (* lr-spin m/M_PI_2))
         ur-pz-sin (m/sin (* ur-spin m/M_PI_2))
         ur-pz-cos (m/sin (* ur-spin m/M_PI_2))
         ll-pz-sin (m/sin (* ll-spin m/M_PI_2))
         ll-pz-cos (m/sin (* ll-spin m/M_PI_2))
         ul-pz-sin (m/sin (* ul-spin m/M_PI_2))
         ul-pz-cos (m/sin (* ul-spin m/M_PI_2))]
     (fn [^Vec2 v]
       (let [lr (if (and (> (.x v) lr-shift-x) (> (.y v) lr-shift-y))
                  (let [re-u (+ (- (* lr-re-a (.x v)) (* lr-im-a (.y v))) lr-re-b)
                        im-u (+ (+ (* lr-re-a (.y v)) (* lr-im-a (.x v))) lr-im-b)
                        re-v (+ (- (* lr-re-c (.x v)) (* lr-im-c (.y v))) lr-re-d)
                        im-v (+ (+ (* lr-re-c (.y v)) (* lr-im-c (.x v))) lr-im-d)
                        d (+ (* re-v re-v) (* im-v im-v))]
                    (Vec3. (+ (* amount (+ (* lr-pz-sin (.y v)) (* lr-pz-cos (.x v))))
                              lr-x (* d (+ (* re-u re-v) (* im-u im-v))))
                           (+ (* amount (+ (* lr-pz-cos (.y v)) (* lr-pz-sin (.x v))))
                              lr-y (* d (+ (* im-u re-v) (* re-u im-v))))
                           1.0))
                  zerov3)
             ur (if (and (> (.x v) ur-shift-x) (< (.y v) ur-shift-y))
                  (let [re-u (+ (- (* ur-re-a (.x v)) (* ur-im-a (.y v))) ur-re-b)
                        im-u (+ (+ (* ur-re-a (.y v)) (* ur-im-a (.x v))) ur-im-b)
                        re-v (+ (- (* ur-re-c (.x v)) (* ur-im-c (.y v))) ur-re-d)
                        im-v (+ (+ (* ur-re-c (.y v)) (* ur-im-c (.x v))) ur-im-d)
                        d (+ (* re-v re-v) (* im-v im-v))]
                    (Vec3. (+ (* amount (+ (* ur-pz-sin (.y v)) (* ur-pz-cos (.x v))))
                              ur-x (* d (+ (* re-u re-v) (* im-u im-v))))
                           (+ (* amount (+ (* ur-pz-cos (.y v)) (* ur-pz-sin (.x v))))
                              ur-y (* d (+ (* im-u re-v) (* re-u im-v))))
                           1.0))
                  zerov3)
             ll (if (and (< (.x v) ll-shift-x) (> (.y v) ll-shift-y))
                  (let [re-u (+ (- (* ll-re-a (.x v)) (* ll-im-a (.y v))) ll-re-b)
                        im-u (+ (+ (* ll-re-a (.y v)) (* ll-im-a (.x v))) ll-im-b)
                        re-v (+ (- (* ll-re-c (.x v)) (* ll-im-c (.y v))) ll-re-d)
                        im-v (+ (+ (* ll-re-c (.y v)) (* ll-im-c (.x v))) ll-im-d)
                        d (+ (* re-v re-v) (* im-v im-v))]
                    (Vec3. (+ (* amount (+ (* ll-pz-sin (.y v)) (* ll-pz-cos (.x v))))
                              ll-x (* d (+ (* re-u re-v) (* im-u im-v))))
                           (+ (* amount (+ (* ll-pz-cos (.y v)) (* ll-pz-sin (.x v))))
                              ll-y (* d (+ (* im-u re-v) (* re-u im-v))))
                           1.0))
                  zerov3)
             ul (if (and (< (.x v) ul-shift-x) (< (.y v) ul-shift-y))
                  (let [re-u (+ (- (* ul-re-a (.x v)) (* ul-im-a (.y v))) ul-re-b)
                        im-u (+ (+ (* ul-re-a (.y v)) (* ul-im-a (.x v))) ul-im-b)
                        re-v (+ (- (* ul-re-c (.x v)) (* ul-im-c (.y v))) ul-re-d)
                        im-v (+ (+ (* ul-re-c (.y v)) (* ul-im-c (.x v))) ul-im-d)
                        d (+ (* re-v re-v) (* im-v im-v))]
                    (Vec3. (+ (* amount (+ (* ul-pz-sin (.y v)) (* ul-pz-cos (.x v))))
                              ul-x (* d (+ (* re-u re-v) (* im-u im-v))))
                           (+ (* amount (+ (* ul-pz-cos (.y v)) (* ul-pz-sin (.x v))))
                              ul-y (* d (+ (* im-u re-v) (* re-u im-v))))
                           1.0))
                  zerov3)
             ^Vec3 sum (v/add (v/add (v/add lr ur) ll) ul)]
         (if (zero? (.z sum))
           (v/mult v amount)
           (if mode
             (Vec2. (/ (/ amount (.x sum)) (.z sum))
                    (/ (/ amount (.y sum)) (.z sum)))
             (Vec2. (/ (* amount (.x sum)) (.z sum))
                    (/ (* amount (.y sum)) (.z sum))))))))))

(defn glynnlissa
  ([] {:type :random
       :config (fn [] {:radius (u/sdrand 0.2 2.0)
                      :radius1 (u/sdrand 0.2 2.0)
                      :thickness (r/drand)
                      :phi (r/drand 360.0)
                      :a (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :b (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :width (r/drand -1.0 1.0)
                      :phase (r/drand m/-TWO_PI m/TWO_PI)
                      :scale (u/sdrand 0.5 1.5)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double radius1 ^double thickness ^double phi
                           ^double a ^double b ^double width ^double phase
                           ^double scale ^double pow ^double contrast]}]
   (let [aa (m/radians phi)
         x1 (* radius (m/cos aa))
         y1 (* radius (m/sin aa))
         xy1 (Vec2. x1 y1)
         abspow (m/abs pow)
         absradius (m/abs radius)
         thickness- (- 1.0 thickness)
         radius1sq (* radius1 radius1)
         l (fn [] (let [t (r/drand m/TWO_PI)
                       y (r/drand -0.5 0.5)]
                   (-> (Vec2. (m/sin (+ phase (* a t))) (m/sin (* b t)))
                       (v/mult scale)
                       (v/shift (* width y))
                       (v/mult radius1)
                       (v/add xy1)
                       (v/mult amount))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (< r absradius)
           (if-not (neg? radius1)
             (l)
             (let [r (* radius1 (+ thickness (r/drand thickness-)))
                   phi (r/drand m/TWO_PI)]
               (v/mult (v/add (Vec2. (* r (m/cos phi)) (* r (m/sin phi))) xy1) amount)))
           (let [alpha (/ absradius r)
                 xyi (if (> (r/drand) (* contrast (m/pow alpha abspow)))
                       v
                       (v/mult v (* alpha alpha)))
                 Z (v/dist-sq xy1 xyi)]
             (if (< Z radius1sq)
               (l)
               (v/mult xyi amount)))))))))

(defn glynnsim1
  ([] {:type :random
       :config (fn [] {:radius (r/drand 0.2 2.0)
                      :radius1 (u/sdrand 0.2 2.0)
                      :phi (r/drand 360.0)
                      :thickness (r/drand)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double radius1 ^double phi
                           ^double thickness ^double pow ^double contrast]}]
   (let [a (m/radians phi)
         sinphi1 (m/sin a)
         cosphi1 (m/cos a)
         x1 (* radius sinphi1)
         y1 (* radius cosphi1)
         abspow (m/abs pow)
         circle-fn #(let [r (* radius1 (+ thickness (* (- 1.0 thickness) (r/drand))))
                          phi (r/drand m/TWO_PI)]
                      (Vec2. (+ x1 (* r (m/sin phi)))
                             (+ y1 (* r (m/cos phi)))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (/ radius (+ m/EPSILON r))]
         (if (< r radius)
           (v/mult (circle-fn) amount)
           (let [^Vec2 toolpoint (if (> (r/drand) (* contrast (m/pow alpha abspow)))
                                   v
                                   (v/mult v (* alpha alpha)))
                 z (+ (m/sq (- (.x toolpoint) x1))
                      (m/sq (- (.y toolpoint) y1)))]
             (if (< z (* radius1 radius1))
               (v/mult (circle-fn) amount)
               (v/mult toolpoint amount)))))))))

(defn glynnsim2
  ([] {:type :random
       :config (fn [] {:radius (r/drand 0.2 2.0)
                      :phi1 (r/drand 360.0)
                      :phi2 (r/drand 360.0)
                      :thickness (r/drand)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double phi1 ^double phi2 ^double thickness
                           ^double pow ^double contrast]}]
   (let [phi10 (m/radians phi1)
         phi20 (m/radians phi2)
         delta (- phi20 phi10)
         gamma (/ (* thickness (+ radius radius thickness))
                  (+ radius thickness))
         abspow (m/abs pow)
         circle-fn #(let [r (- (+ radius thickness)
                               (* gamma (r/drand)))
                          phi (+ phi10 (* delta (r/drand)))]
                      (Vec2. (* r (m/sin phi))
                             (* r (m/cos phi))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (/ radius (+ m/EPSILON r))]
         (if (< r radius)
           (v/mult (circle-fn) amount)
           (if (> (r/drand) (* contrast (m/pow alpha abspow)))
             v
             (v/mult v (* alpha alpha)))))))))

(defn glynnsim3
  ([] {:type :random
       :config (fn [] {:radius (r/drand 0.2 2.0)
                      :thickness (r/drand)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double thickness ^double pow ^double contrast]}]
   (let [radius1 (+ radius thickness)
         radius2 (/ (m/sq radius) radius1)
         gamma (/ radius1
                  (+ radius1 radius2))
         abspow (m/abs pow)
         circle-fn #(let [phi (r/drand m/TWO_PI)
                          r (if (< (r/drand) gamma) radius1 radius2)]
                      (Vec2. (* r (m/sin phi))
                             (* r (m/cos phi))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (/ radius (+ m/EPSILON r))]
         (if (< r radius1)
           (v/mult (circle-fn) amount)
           (if (> (r/drand) (* contrast (m/pow alpha abspow)))
             v
             (v/mult v (* alpha alpha)))))))))

(defn glynnsupershape
  ([] {:type :random
       :config (fn [] {:radius (u/sdrand 0.2 2.0)
                      :radius1 (u/sdrand 0.2 2.0)
                      :thickness (r/drand)
                      :phi1 (r/drand 360.0)
                      :m (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :n1 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :n2 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :n3 (r/randval (u/sirand 1 10) (u/sdrand 0.125 10.0))
                      :scale (u/sdrand 0.5 1.5)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double radius1 ^double thickness ^double phi1
                           ^double m ^double n1 ^double n2 ^double n3
                           ^double scale ^double pow ^double contrast]}]
   (let [aa (m/radians phi1)
         x1 (* radius (m/cos aa))
         y1 (* radius (m/sin aa))
         xy1 (Vec2. x1 y1)
         abspow (m/abs pow)
         absradius (m/abs radius)
         thickness- (- 1.0 thickness)
         radius1sq (* radius1 radius1)
         m4 (* m 0.25)
         n1r (/ n1)
         l (fn [] (let [phi (r/drand m/TWO_PI)
                       phim4 (* phi m4)
                       t1 (-> (m/cos phim4)
                              (m/abs)
                              (m/pow n2))
                       t2 (-> (m/sin phim4)
                              (m/abs)
                              (m/pow n3))
                       r (m/pow (+ t1 t2) n1r)]
                   (-> (if (zero? r) u/zerov (v/div (Vec2. (m/cos phi) (m/sin phi)) r))
                       (v/mult scale)
                       (v/mult radius1)
                       (v/add xy1)
                       (v/mult amount))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (< r absradius)
           (if-not (neg? radius1)
             (l)
             (let [r (* radius1 (+ thickness (r/drand thickness-)))
                   phi (r/drand m/TWO_PI)]
               (v/mult (v/add (Vec2. (* r (m/cos phi)) (* r (m/sin phi))) xy1) amount)))
           (let [alpha (/ absradius r)
                 xyi (if (> (r/drand) (* contrast (m/pow alpha abspow)))
                       v
                       (v/mult v (* alpha alpha)))
                 Z (v/dist-sq xy1 xyi)]
             (if (< Z radius1sq)
               (l)
               (v/mult xyi amount)))))))))

(defn glynnia3
  ([] {:type :random
       :config (fn [] {:rscale (r/drand 0.1 2.0)
                      :dscale (u/sdrand 0.1 2.0)
                      :rthresh (r/drand)
                      :ythresh (r/randval 0.0 (r/drand -0.5 0.5))})})
  ([^double amount {:keys [^double rscale ^double dscale ^double rthresh ^double ythresh]}]
   (let [vvar2 (* amount m/SQRT2_2)]
     (fn [^Vec2 v]
       (let [r (* rscale (v/mag v))]
         (if (and (> r rthresh) (> (.y v) ythresh))
           (if (r/brand)
             (let [d (* dscale (m/sqrt (+ r (.x v))))]
               (Vec2. (* vvar2 d)
                      (- (* (.y v) (/ vvar2 d)))))
             (let [d (* dscale (+ r (.x v)))
                   dx (m/sqrt (* r (+ (* (.y v) (.y v))
                                      (* d d))))
                   r (/ amount dx)]
               (Vec2. (* r d) (* r (.y v)))))
           (if (r/brand)
             (let [d (* dscale (m/sqrt (+ r (.x v))))]
               (Vec2. (- (* vvar2 d))
                      (- (* (.y v) (/ vvar2 d)))))
             (let [d (* dscale (+ r (.x v)))
                   dx (m/sqrt (* r (+ (* (.y v) (.y v))
                                      (* d d))))
                   r (/ amount dx)]
               (Vec2. (- (* r d)) (* r (.y v)))))))))))

(defn glynnia
  ([] {:type :random})
  ([^double amount _]
   (let [vvar2 (* amount m/SQRT2_2)]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (>= r 1.0)
           (if (r/brand)
             (let [d (m/sqrt (+ r (.x v)))]
               (Vec2. (* vvar2 d)
                      (- (* (.y v) (/ vvar2 d)))))
             (let [d (+ r (.x v))
                   dx (m/sqrt (* r (+ (* (.y v) (.y v))
                                      (* d d))))
                   r (/ amount dx)]
               (Vec2. (* r d) (* r (.y v)))))
           (if (r/brand)
             (let [d (m/sqrt (+ r (.x v)))]
               (Vec2. (- (* vvar2 d))
                      (- (* (.y v) (/ vvar2 d)))))
             (let [d (+ r (.x v))
                   dx (m/sqrt (* r (+ (* (.y v) (.y v))
                                      (* d d))))
                   r (/ amount dx)]
               (Vec2. (- (* r d)) (* r (.y v)))))))))))

(defn gridout2
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (let [x (* (m/rint (.x v)) c)
           y (* (m/rint (.y v)) d)]
       (-> (if-not (pos? y)
             (if (pos? x)
               (if (>= (- y) x)
                 (Vec2. (+ (.x v) a) (.y v))
                 (Vec2. (.x v) (+ (.y v) b)))
               (if (<= y x)
                 (Vec2. (+ (.x v) a) (.y v))
                 (Vec2. (.x v) (- (.y v) b))))
             (if (pos? x)
               (if (>= y x)
                 (Vec2. (- (.x v) a) (.y v))
                 (Vec2. (.x v) (+ (.y v) b)))
               (if (> y (- x))
                 (Vec2. (- (.x v) a) (.y v))
                 (Vec2. (.x v) (- (.y v) b)))))
           (v/mult amount))))))

(defn gridout
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [x (m/rint (.x v))
           y (m/rint (.y v))]
       (-> (if-not (pos? y)
             (if (pos? x)
               (if (>= (- y) x)
                 (Vec2. (inc (.x v)) (.y v))
                 (Vec2. (.x v) (inc (.y v))))
               (if (<= y x)
                 (Vec2. (inc (.x v)) (.y v))
                 (Vec2. (.x v) (dec (.y v)))))
             (if (pos? x)
               (if (>= y x)
                 (Vec2. (dec (.x v)) (.y v))
                 (Vec2. (.x v) (inc (.y v))))
               (if (> y (- x))
                 (Vec2. (dec (.x v)) (.y v))
                 (Vec2. (.x v) (dec (.y v))))))
           (v/mult amount))))))

;;

(defn general-noise
  "Perlin noise"
  ([] {:type :regular
       :config (fn [] (assoc (r/random-noise-cfg) :scale (u/sdrand 0.1 1.5)))})
  ([amount cfg]
   (let [n (r/random-noise-fn cfg)]
     (u/make-noise-variation amount (:scale cfg) 2.0 n))))

(defn general-noise2
  "Perlin noise"
  ([] {:type :regular
       :config (fn [] (assoc (r/random-noise-cfg) :scale (u/sdrand 0.1 1.5)))})
  ([amount cfg]
   (let [n (r/random-noise-fn cfg)]
     (u/make-noise-variation2 amount (:scale cfg) n))))

(m/unuse-primitive-operators)
