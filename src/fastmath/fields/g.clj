(ns fastmath.fields.g
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.special :as special]
            [fastmath.fields.utils :as u])
  (:import [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- fclp ^double [^double a] (if (m/neg? a) (m/- (m/rem (m/abs a) 1.0)) (m/rem (m/abs a) 1.0)))
(defn- fscl ^double [^double a] (fclp (m/* 0.5 (m/inc a))))
(defn- fosc ^double [^double p ^double a] (fscl (m/- (m/cos (m/* p a m/TWO_PI)))))
(defn- flip ^double [^double a ^double b ^double c] (m/+ a (m/* c (m/- b a))))

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
   (let [gdodx (m/* delta-x 0.1)
         gdody (m/* delta-y 0.1)
         gdoax (m/* 2.0 (if (m/< (m/abs area-x) 0.1) 0.1 (m/abs area-x)))
         gdoay (m/* 2.0 (if (m/< (m/abs area-y) 0.1) 0.1 (m/abs area-y)))
         gdocx center-x
         gdocy center-y
         gdog gamma
         gdos square
         gdob (m// (m/* gdog 2.0) (m/max gdoax gdoay))]
     (fn [^Vec2 v]
       (let [osc-x (fosc gdodx 1.0)
             osc-y (if gdos osc-x (fosc gdody 1.0))
             in-x (m/+ (.x v) gdocx)
             in-y (m/+ (.y v) gdocy)]
         (-> (Vec2. (flip (flip in-x (fosc in-x 4.0) osc-x) (fosc (fclp (m/* gdob in-x)) 4.0) osc-x)
                    (flip (flip in-y (fosc in-y 4.0) osc-y) (fosc (fclp (m/* gdob in-y)) 4.0) osc-y))
             (v/mult amount)))))))

(defn gamma
  "gamma by zephyrtronium, http://fractal-resources.deviantart.com/art/Gamma-Apophysis-Plugin-154808483"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (Vec2. (m/* amount (special/log-gamma (v/mag v)))
            (m/* amount (v/heading v))))))

(defn gaussianblur
  "Gaussian"
  ([] {:type :pattern})
  ([^double amount _]
   (fn [_]
     (let [a (r/drand m/TWO_PI)
           r (m/* amount (m/+ (r/drand) (r/drand) (r/drand) (r/drand) -2.0))]
       (Vec2. (m/* r (m/cos a)) (m/* r (m/sin a)))))))

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
   (let [lr-pz-sin (m/sin (m/* lr-spin m/M_PI_2))
         lr-pz-cos (m/sin (m/* lr-spin m/M_PI_2))
         ur-pz-sin (m/sin (m/* ur-spin m/M_PI_2))
         ur-pz-cos (m/sin (m/* ur-spin m/M_PI_2))
         ll-pz-sin (m/sin (m/* ll-spin m/M_PI_2))
         ll-pz-cos (m/sin (m/* ll-spin m/M_PI_2))
         ul-pz-sin (m/sin (m/* ul-spin m/M_PI_2))
         ul-pz-cos (m/sin (m/* ul-spin m/M_PI_2))]
     (fn [^Vec2 v]
       (let [lr (if (and (m/> (.x v) lr-shift-x) (m/> (.y v) lr-shift-y))
                  (let [re-u (m/+ (m/- (m/* lr-re-a (.x v)) (m/* lr-im-a (.y v))) lr-re-b)
                        im-u (m/+ (m/+ (m/* lr-re-a (.y v)) (m/* lr-im-a (.x v))) lr-im-b)
                        re-v (m/+ (m/- (m/* lr-re-c (.x v)) (m/* lr-im-c (.y v))) lr-re-d)
                        im-v (m/+ (m/+ (m/* lr-re-c (.y v)) (m/* lr-im-c (.x v))) lr-im-d)
                        d (m/+ (m/* re-v re-v) (m/* im-v im-v))]
                    (Vec3. (m/+ (m/* amount (m/+ (m/* lr-pz-sin (.y v)) (m/* lr-pz-cos (.x v))))
                              lr-x (m/* d (m/+ (m/* re-u re-v) (m/* im-u im-v))))
                           (m/+ (m/* amount (m/+ (m/* lr-pz-cos (.y v)) (m/* lr-pz-sin (.x v))))
                              lr-y (m/* d (m/+ (m/* im-u re-v) (m/* re-u im-v))))
                           1.0))
                  zerov3)
             ur (if (and (m/> (.x v) ur-shift-x) (m/< (.y v) ur-shift-y))
                  (let [re-u (m/+ (m/- (m/* ur-re-a (.x v)) (m/* ur-im-a (.y v))) ur-re-b)
                        im-u (m/+ (m/+ (m/* ur-re-a (.y v)) (m/* ur-im-a (.x v))) ur-im-b)
                        re-v (m/+ (m/- (m/* ur-re-c (.x v)) (m/* ur-im-c (.y v))) ur-re-d)
                        im-v (m/+ (m/+ (m/* ur-re-c (.y v)) (m/* ur-im-c (.x v))) ur-im-d)
                        d (m/+ (m/* re-v re-v) (m/* im-v im-v))]
                    (Vec3. (m/+ (m/* amount (m/+ (m/* ur-pz-sin (.y v)) (m/* ur-pz-cos (.x v))))
                              ur-x (m/* d (m/+ (m/* re-u re-v) (m/* im-u im-v))))
                           (m/+ (m/* amount (m/+ (m/* ur-pz-cos (.y v)) (m/* ur-pz-sin (.x v))))
                              ur-y (m/* d (m/+ (m/* im-u re-v) (m/* re-u im-v))))
                           1.0))
                  zerov3)
             ll (if (and (m/< (.x v) ll-shift-x) (m/> (.y v) ll-shift-y))
                  (let [re-u (m/+ (m/- (m/* ll-re-a (.x v)) (m/* ll-im-a (.y v))) ll-re-b)
                        im-u (m/+ (m/+ (m/* ll-re-a (.y v)) (m/* ll-im-a (.x v))) ll-im-b)
                        re-v (m/+ (m/- (m/* ll-re-c (.x v)) (m/* ll-im-c (.y v))) ll-re-d)
                        im-v (m/+ (m/+ (m/* ll-re-c (.y v)) (m/* ll-im-c (.x v))) ll-im-d)
                        d (m/+ (m/* re-v re-v) (m/* im-v im-v))]
                    (Vec3. (m/+ (m/* amount (m/+ (m/* ll-pz-sin (.y v)) (m/* ll-pz-cos (.x v))))
                              ll-x (m/* d (m/+ (m/* re-u re-v) (m/* im-u im-v))))
                           (m/+ (m/* amount (m/+ (m/* ll-pz-cos (.y v)) (m/* ll-pz-sin (.x v))))
                              ll-y (m/* d (m/+ (m/* im-u re-v) (m/* re-u im-v))))
                           1.0))
                  zerov3)
             ul (if (and (m/< (.x v) ul-shift-x) (m/< (.y v) ul-shift-y))
                  (let [re-u (m/+ (m/- (m/* ul-re-a (.x v)) (m/* ul-im-a (.y v))) ul-re-b)
                        im-u (m/+ (m/+ (m/* ul-re-a (.y v)) (m/* ul-im-a (.x v))) ul-im-b)
                        re-v (m/+ (m/- (m/* ul-re-c (.x v)) (m/* ul-im-c (.y v))) ul-re-d)
                        im-v (m/+ (m/+ (m/* ul-re-c (.y v)) (m/* ul-im-c (.x v))) ul-im-d)
                        d (m/+ (m/* re-v re-v) (m/* im-v im-v))]
                    (Vec3. (m/+ (m/* amount (m/+ (m/* ul-pz-sin (.y v)) (m/* ul-pz-cos (.x v))))
                              ul-x (m/* d (m/+ (m/* re-u re-v) (m/* im-u im-v))))
                           (m/+ (m/* amount (m/+ (m/* ul-pz-cos (.y v)) (m/* ul-pz-sin (.x v))))
                              ul-y (m/* d (m/+ (m/* im-u re-v) (m/* re-u im-v))))
                           1.0))
                  zerov3)
             ^Vec3 sum (v/add (v/add (v/add lr ur) ll) ul)]
         (if (m/zero? (.z sum))
           (v/mult v amount)
           (if mode
             (Vec2. (m// (m// amount (.x sum)) (.z sum))
                    (m// (m// amount (.y sum)) (.z sum)))
             (Vec2. (m// (m/* amount (.x sum)) (.z sum))
                    (m// (m/* amount (.y sum)) (.z sum))))))))))

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
         x1 (m/* radius (m/cos aa))
         y1 (m/* radius (m/sin aa))
         xy1 (Vec2. x1 y1)
         abspow (m/abs pow)
         absradius (m/abs radius)
         thickness- (m/- 1.0 thickness)
         radius1sq (m/* radius1 radius1)
         l (fn [] (let [t (r/drand m/TWO_PI)
                       y (r/drand -0.5 0.5)]
                   (-> (Vec2. (m/sin (m/+ phase (m/* a t))) (m/sin (m/* b t)))
                       (v/mult scale)
                       (v/shift (m/* width y))
                       (v/mult radius1)
                       (v/add xy1)
                       (v/mult amount))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (m/< r absradius)
           (if-not (m/neg? radius1)
             (l)
             (let [r (m/* radius1 (m/+ thickness (r/drand thickness-)))
                   phi (r/drand m/TWO_PI)]
               (v/mult (v/add (Vec2. (m/* r (m/cos phi)) (m/* r (m/sin phi))) xy1) amount)))
           (let [alpha (m// absradius r)
                 xyi (if (m/> (r/drand) (m/* contrast (m/pow alpha abspow)))
                       v
                       (v/mult v (m/* alpha alpha)))
                 Z (v/dist-sq xy1 xyi)]
             (if (m/< Z radius1sq)
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
         x1 (m/* radius sinphi1)
         y1 (m/* radius cosphi1)
         abspow (m/abs pow)
         circle-fn #(let [r (m/* radius1 (m/+ thickness (m/* (m/- 1.0 thickness) (r/drand))))
                          phi (r/drand m/TWO_PI)]
                      (Vec2. (m/+ x1 (m/* r (m/sin phi)))
                             (m/+ y1 (m/* r (m/cos phi)))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (m// radius (m/+ m/EPSILON r))]
         (if (m/< r radius)
           (v/mult (circle-fn) amount)
           (let [^Vec2 toolpoint (if (m/> (r/drand) (m/* contrast (m/pow alpha abspow)))
                                   v
                                   (v/mult v (m/* alpha alpha)))
                 z (m/+ (m/sq (m/- (.x toolpoint) x1))
                      (m/sq (m/- (.y toolpoint) y1)))]
             (if (m/< z (m/* radius1 radius1))
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
         delta (m/- phi20 phi10)
         gamma (m// (m/* thickness (m/+ radius radius thickness))
                  (m/+ radius thickness))
         abspow (m/abs pow)
         circle-fn #(let [r (m/- (m/+ radius thickness)
                               (m/* gamma (r/drand)))
                          phi (m/+ phi10 (m/* delta (r/drand)))]
                      (Vec2. (m/* r (m/sin phi))
                             (m/* r (m/cos phi))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (m// radius (m/+ m/EPSILON r))]
         (if (m/< r radius)
           (v/mult (circle-fn) amount)
           (if (m/> (r/drand) (m/* contrast (m/pow alpha abspow)))
             v
             (v/mult v (m/* alpha alpha)))))))))

(defn glynnsim3
  ([] {:type :random
       :config (fn [] {:radius (r/drand 0.2 2.0)
                      :thickness (r/drand)
                      :pow (r/drand -6.0 6.0)
                      :contrast (r/drand)})})
  ([^double amount {:keys [^double radius ^double thickness ^double pow ^double contrast]}]
   (let [radius1 (m/+ radius thickness)
         radius2 (m// (m/sq radius) radius1)
         gamma (m// radius1
                  (m/+ radius1 radius2))
         abspow (m/abs pow)
         circle-fn #(let [phi (r/drand m/TWO_PI)
                          r (if (m/< (r/drand) gamma) radius1 radius2)]
                      (Vec2. (m/* r (m/sin phi))
                             (m/* r (m/cos phi))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)
             alpha (m// radius (m/+ m/EPSILON r))]
         (if (m/< r radius1)
           (v/mult (circle-fn) amount)
           (if (m/> (r/drand) (m/* contrast (m/pow alpha abspow)))
             v
             (v/mult v (m/* alpha alpha)))))))))

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
         x1 (m/* radius (m/cos aa))
         y1 (m/* radius (m/sin aa))
         xy1 (Vec2. x1 y1)
         abspow (m/abs pow)
         absradius (m/abs radius)
         thickness- (m/- 1.0 thickness)
         radius1sq (m/* radius1 radius1)
         m4 (m/* m 0.25)
         n1r (m// n1)
         l (fn [] (let [phi (r/drand m/TWO_PI)
                       phim4 (m/* phi m4)
                       t1 (-> (m/cos phim4)
                              (m/abs)
                              (m/pow n2))
                       t2 (-> (m/sin phim4)
                              (m/abs)
                              (m/pow n3))
                       r (m/pow (m/+ t1 t2) n1r)]
                   (-> (if (m/zero? r) u/zerov (v/div (Vec2. (m/cos phi) (m/sin phi)) r))
                       (v/mult scale)
                       (v/mult radius1)
                       (v/add xy1)
                       (v/mult amount))))]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (m/< r absradius)
           (if-not (m/neg? radius1)
             (l)
             (let [r (m/* radius1 (m/+ thickness (r/drand thickness-)))
                   phi (r/drand m/TWO_PI)]
               (v/mult (v/add (Vec2. (m/* r (m/cos phi)) (m/* r (m/sin phi))) xy1) amount)))
           (let [alpha (m// absradius r)
                 xyi (if (m/> (r/drand) (m/* contrast (m/pow alpha abspow)))
                       v
                       (v/mult v (m/* alpha alpha)))
                 Z (v/dist-sq xy1 xyi)]
             (if (m/< Z radius1sq)
               (l)
               (v/mult xyi amount)))))))))

(defn glynnia3
  ([] {:type :random
       :config (fn [] {:rscale (r/drand 0.1 2.0)
                      :dscale (u/sdrand 0.1 2.0)
                      :rthresh (r/drand)
                      :ythresh (r/randval 0.0 (r/drand -0.5 0.5))})})
  ([^double amount {:keys [^double rscale ^double dscale ^double rthresh ^double ythresh]}]
   (let [vvar2 (m/* amount m/SQRT2_2)]
     (fn [^Vec2 v]
       (let [r (m/* rscale (v/mag v))]
         (if (and (m/> r rthresh) (m/> (.y v) ythresh))
           (if (r/brand)
             (let [d (m/* dscale (m/sqrt (m/+ r (.x v))))]
               (Vec2. (m/* vvar2 d)
                      (m/- (m/* (.y v) (m// vvar2 d)))))
             (let [d (m/* dscale (m/+ r (.x v)))
                   dx (m/sqrt (m/* r (m/+ (m/* (.y v) (.y v))
                                      (m/* d d))))
                   r (m// amount dx)]
               (Vec2. (m/* r d) (m/* r (.y v)))))
           (if (r/brand)
             (let [d (m/* dscale (m/sqrt (m/+ r (.x v))))]
               (Vec2. (m/- (m/* vvar2 d))
                      (m/- (m/* (.y v) (m// vvar2 d)))))
             (let [d (m/* dscale (m/+ r (.x v)))
                   dx (m/sqrt (m/* r (m/+ (m/* (.y v) (.y v))
                                      (m/* d d))))
                   r (m// amount dx)]
               (Vec2. (m/- (m/* r d)) (m/* r (.y v)))))))))))

(defn glynnia
  ([] {:type :random})
  ([^double amount _]
   (let [vvar2 (m/* amount m/SQRT2_2)]
     (fn [^Vec2 v]
       (let [r (v/mag v)]
         (if (m/>= r 1.0)
           (if (r/brand)
             (let [d (m/sqrt (m/+ r (.x v)))]
               (Vec2. (m/* vvar2 d)
                      (m/- (m/* (.y v) (m// vvar2 d)))))
             (let [d (m/+ r (.x v))
                   dx (m/sqrt (m/* r (m/+ (m/* (.y v) (.y v))
                                      (m/* d d))))
                   r (m// amount dx)]
               (Vec2. (m/* r d) (m/* r (.y v)))))
           (if (r/brand)
             (let [d (m/sqrt (m/+ r (.x v)))]
               (Vec2. (m/- (m/* vvar2 d))
                      (m/- (m/* (.y v) (m// vvar2 d)))))
             (let [d (m/+ r (.x v))
                   dx (m/sqrt (m/* r (m/+ (m/* (.y v) (.y v))
                                      (m/* d d))))
                   r (m// amount dx)]
               (Vec2. (m/- (m/* r d)) (m/* r (.y v)))))))))))

(defn gridout2
  ([] {:type :regular
       :config (fn [] {:a (r/drand -2.0 2.0)
                      :b (r/drand -2.0 2.0)
                      :c (r/drand -2.0 2.0)
                      :d (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (let [x (m/* (m/rint (.x v)) c)
           y (m/* (m/rint (.y v)) d)]
       (-> (if-not (m/pos? y)
             (if (m/pos? x)
               (if (m/>= (m/- y) x)
                 (Vec2. (m/+ (.x v) a) (.y v))
                 (Vec2. (.x v) (m/+ (.y v) b)))
               (if (m/<= y x)
                 (Vec2. (m/+ (.x v) a) (.y v))
                 (Vec2. (.x v) (m/- (.y v) b))))
             (if (m/pos? x)
               (if (m/>= y x)
                 (Vec2. (m/- (.x v) a) (.y v))
                 (Vec2. (.x v) (m/+ (.y v) b)))
               (if (m/> y (m/- x))
                 (Vec2. (m/- (.x v) a) (.y v))
                 (Vec2. (.x v) (m/- (.y v) b)))))
           (v/mult amount))))))

(defn gridout
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [x (m/rint (.x v))
           y (m/rint (.y v))]
       (-> (if-not (m/pos? y)
             (if (m/pos? x)
               (if (m/>= (m/- y) x)
                 (Vec2. (m/inc (.x v)) (.y v))
                 (Vec2. (.x v) (m/inc (.y v))))
               (if (m/<= y x)
                 (Vec2. (m/inc (.x v)) (.y v))
                 (Vec2. (.x v) (m/dec (.y v)))))
             (if (m/pos? x)
               (if (m/>= y x)
                 (Vec2. (m/dec (.x v)) (.y v))
                 (Vec2. (.x v) (m/inc (.y v))))
               (if (m/> y (m/- x))
                 (Vec2. (m/dec (.x v)) (.y v))
                 (Vec2. (.x v) (m/dec (.y v))))))
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

