(ns fastmath.curves
  "Collection of parametric curves"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- srandom
  "Symetric random from [-mx -mn] and [mn mx]"
  ^double  [^double mn ^double mx]
  (let [rand (r/drand mn mx)]
    (r/randval rand (* -1.0 rand))))

(defn- mk-random
  ([] (r/randval r/drand r/irand))
  ([prob] (r/randval prob r/drand r/irand))
  ([prob round-fn] (r/randval prob r/drand (comp round-fn r/drand))))

(defmulti parametrization
  "Return random parametrization map for given curve.

  You can pass part of the parametrization. In this case function will add remaining keys with randomly generated values.

  If field doesn't have any parametrization, empty map will be returned.
  
  See [[curve]]."
  (fn [key & _] key))

(defmethod parametrization :default
  ([_ _] {})
  ([_] {}))

(defmulti curve
  "Return vector field for given name and parametrization.

  Default parametrization is random.

  Resulting function returns [[Vec2]] type."
  (fn [key & _] key))

(defmethod curve :default
  ([_ _] identity)
  ([_] identity))

(defmacro ^:private make-config-method
  "Add new multimethod for variation parametrization"
  [sym m]
  (let [k (keyword sym)]
    `(defmethod parametrization ~k
       ([k# p#] (merge ~m p#))
       ([k#] ~m))))

(defmacro ^:private make-curve-method
  "Add new multimethod for curve factory function."
  [sym]
  (let [k (keyword sym)
        m (symbol (str "make-" sym))]
    `(defmethod curve ~k
       ([k# p#] (~m (parametrization ~k p#)))
       ([k#] (~m (parametrization ~k))))))

(make-config-method astroid {:a 1.0})
(defn- make-astroid
  [{:keys [^double a]}]
  (fn [^double t]
    (Vec2. (* a (m/cb (m/cos t)))
           (* a (m/cb (m/sin t))))))
(make-curve-method astroid)

(make-config-method astroid-pedal-curve (let [rx0 (mk-random)
                                              ry0 (mk-random)]
                                          {:a 1.0
                                           :x0 (rx0 -4.0 4.0)
                                           :y0 (ry0 -4.0 4.0)}))
(defn- make-astroid-pedal-curve
  [{:keys [^double a ^double x0 ^double y0]}]
  (let [x02 (+ x0 x0)
        y02 (+ y0 y0)]
    (fn [^double t]
      (let [sint (m/sin t)
            cost (m/cos t)]
        (Vec2. (* cost (+ (* sint (- (* a sint) y0))
                          (* x0 cost)))
               (* 0.5 sint (- (+ (* a (m/cos (+ t t)))
                                 a
                                 (* y02 sint))
                              (* x02 cost))))))))
(make-curve-method astroid-pedal-curve)

(make-config-method cardioid-pedal-curve (let [rx0 (mk-random)
                                               ry0 (mk-random)]
                                           {:a 1.0
                                            :x0 (rx0 -4.0 4.0)
                                            :y0 (ry0 -4.0 4.0)}))
(defn- make-cardioid-pedal-curve
  [{:keys [^double a ^double x0 ^double y0]}]
  (let [x02 (+ x0 x0)
        y02 (+ y0 y0)
        a3 (* 3.0 a)]
    (fn [^double t]
      (let [t2 (+ t t)
            t3 (+ t2 t)]
        (Vec2. (* 0.25 (- (+ (* (+ a x02) (m/cos t3))
                             (* a3 (m/cos t))
                             (* y02 (m/sin t3))
                             x02)
                          (* a3 (m/cos t2))
                          a))
               (* 0.25 (- (+ (* (+ a x02) (m/sin t3))
                             (* a3 (m/sin t))
                             y02)
                          (* a3 (m/sin t2))
                          (* y02 (m/cos t3)))))))))
(make-curve-method cardioid-pedal-curve)

(make-config-method lissajous (let [rkx (mk-random 0.3 m/ceil)
                                    rky (mk-random 0.3 m/ceil)]
                                {:a 1.0 :b 1.0
                                 :kx (rkx 0.1 8.0)
                                 :ky (rky 0.1 8.0)}))
(defn- make-lissajous
  [{:keys [^double a ^double b ^double kx ^double ky]}]
  (fn [^double t]
    (Vec2. (* a (m/cos (* kx t)))
           (* b (m/sin (* ky t))))))
(make-curve-method lissajous)

(make-config-method piriform {:a (srandom 0.5 1.5) :b (srandom 0.5 1.5)})
(defn- make-piriform
  [{:keys [^double a ^double b]}]
  (fn [^double t]
    (let [sint+1 (inc (m/sin t))]
      (Vec2. (* a sint+1)
             (* b sint+1 (m/cos t))))))
(make-curve-method piriform)

(make-config-method superellipse-general {:a 1.0 :b 1.0 :y (srandom 0.01 m/TWO_PI) :z (srandom 0.01 m/TWO_PI)
                                          :n1 (r/drand -5.0 5.0) :n2 (r/drand -5.0 5.0) :n3 (r/drand -5.0 5.0)})
(defn- make-superellipse-general
  [{:keys [^double a ^double b ^double y ^double z ^double n1 ^double n2 ^double n3]}]
  (let [ra (/ a)
        rb (/ b)
        rn1 (/ n1)]
    (fn [^double t]
      (v/from-polar (Vec2. (/ (m/pow (+ (m/pow (m/abs (* ra (m/cos (* 0.25 y t)))) n2)
                                        (m/pow (m/abs (* rb (m/sin (* 0.25 z t)))) n3)) rn1))
                           t)))))
(make-curve-method superellipse-general)

(make-config-method superellipse {:a 1.0 :b 1.0 :m (srandom 0.01 m/TWO_PI)
                                  :n1 (r/drand -5.0 5.0) :n2 (r/drand -5.0 5.0) :n3 (r/drand -5.0 5.0)})
(defn- make-superellipse
  [params]
  (make-superellipse-general (-> params
                                 (assoc :y (:m params))
                                 (assoc :z (:m params)))))
(make-curve-method superellipse)

(make-config-method superformula-general {:a 1.0 :b 1.0 :y (srandom 0.01 m/TWO_PI) :z (srandom 0.01 m/TWO_PI)
                                          :n1 (r/drand -5.0 5.0) :n2 (r/drand -5.0 5.0) :n3 (r/drand -5.0 5.0)})
(defn- make-superformula-general
  [{:keys [^double a ^double b ^double y ^double z ^double n1 ^double n2 ^double n3]}]
  (let [ra (/ a)
        rb (/ b)
        -rn1 (/ -1.0 n1)]
    (fn [^double t]
      (v/from-polar (Vec2. (m/pow (+ (m/pow (m/abs (* ra (m/cos (* 0.25 y t)))) n2)
                                     (m/pow (m/abs (* rb (m/sin (* 0.25 z t)))) n3)) -rn1)
                           t)))))
(make-curve-method superformula-general)

(make-config-method superformula {:a 1.0 :b 1.0 :m (srandom 0.01 m/TWO_PI)
                                  :n1 (r/drand -5.0 5.0) :n2 (r/drand -5.0 5.0) :n3 (r/drand -5.0 5.0)})
(defn- make-superformula
  [params]
  (make-superformula-general (-> params
                                 (assoc :y (:m params))
                                 (assoc :z (:m params)))))
(make-curve-method superformula)

(make-config-method hypotrochoid {:R (srandom 0.1 2.0) :r (srandom 0.1 2.0) :d (srandom 0.5 2.0)})
(defn- make-hypotrochoid
  [{:keys [^double R ^double r ^double d]}]
  (fn [^double t]
    (let [diff (- R r)
          ratio (* t (/ diff r))]
      (Vec2. (+ (* diff (m/cos t))
                (* d (m/cos ratio)))
             (- (* diff (m/sin t))
                (* d (m/sin ratio)))))))
(make-curve-method hypotrochoid)

(defn mult
  "Multiply two curves (by multiplying vectors elementwise)."
  [c1 c2]
  (fn [^double t]
    (v/emult (c1 t) (c2 t))))

(defn add
  "Add two curves (by adding resulting vectors)."
  [c1 c2]
  (fn [^double t]
    (v/add (c1 t) (c2 t))))

(def curves-list (sort (keys (methods curve))))
