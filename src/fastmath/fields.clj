(ns fastmath.fields
  "Vector field functions.

  Vector fields are functions R^2->R^2.

  Names are taken from fractal flames world where such fields are call `variations`. Most implementations are taken from [JWildfire](http://jwildfire.org/) software.

  ### Creation

  To create vector field call [[field]] multimethod with name of the field as keyword.

  Some of the vector fields require additional configuration as a map of parameters as keywords and values. Call [[parametrization]] to create random one or to merge with provided.

  Additionally you can provide `amount` parameter which is scaling factor for vector field (default: `1.0`).

  ### Derived fields

  You can use several method to derive new vector field from the other one(s). Possible options are:

  * [[derivative]], [[grad-x]], [[grad-y]] - directional derivative of the field
  * [[sum]] - sum of two fields
  * [[multiply]] - multiplication of the fields
  * [[composition]] - composition of the fields
  * [[angles]] - angles of the field vectors

  ### Scalar fields

  You can derive scalar fields from given vector field(s):

  * [[jacobian]] - determinant of jacobian matrix
  * [[divergence]] - divergence of the field
  * [[cross]] - cross product of the fields (as a determinant of the 2x2 matrix of vectors)
  * [[dot]] - dot product
  * [[angle-between]] - angle between vectors from fields.

  ### Combinations

  The other option is to create vector field using some of the above possibilities. Combination is a tree of field operations with parametrizations. Functions:

  * [[combine]] - create vector field randomly of from given parametrization.
  * [[random-configuration]] - returns random configuration as a map
  * [[randomize-configuration]] - change parametrization for given configuration."
  {:metadoc/categories {:cr "Create fields"
                        :sc "Derive scalar field from vector field"
                        :vf "Derive vector field from other vector field(s)."}}
  (:require [fastmath.core :as m]
            [fastmath.random :refer [brand discrete-noise drand fbm-noise irand lrand noise randval] :exclude [flip] :as r]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Two atoms to store variation names. One for non-random functions and second for random.

(def ^:private fields-atom (atom {:regular [:default]}))

(defn- register-field
  "Add `name` to the atom `what`"
  [type name]
  (swap! fields-atom update type conj name))

(defmulti parametrization
  "Return random parametrization map for given field.

  Optinally you can pass part of the parametrization. In this case function will add remaining keys with randomly generated values.

  If field doesn't have parametrization, empty map will be returned.
  
  See [[field]]."
  {:metadoc/categories #{:cr}}
  (fn [key & _] key))

(defmethod parametrization :default
  ([_ _] {})
  ([_] {}))

(defmulti field
  "Return vector field for given name and options: amount (scaling factor) and parametrization.

  Default scaling factor is 1.0, default parametrization is random.

  Resulting function operates on [[Vec2]] type."
  {:metadoc/categories #{:cr}}
  (fn [key & _] key))

(defmethod field :default
  ([_ s _] (fn [v] (v/mult v s)))
  ([_ s] (fn [v] (v/mult v s)))
  ([_] identity))

(defmacro ^:private make-config-method
  "Add new multimethod for variation parametrization"
  [sym m]
  (when (and sym m)
    (let [k (keyword sym)]
      `(defmethod parametrization ~k
         ([k# p#] (merge ~m p#))
         ([k#] ~m)))))

(defmacro ^:private make-field-method
  "Add new multimethod for variation factory function and store variation fn in global list."
  [sym t]
  (let [k (keyword sym)
        m (symbol (str "make-" sym))]
    `(do (defmethod field ~k
           ([k# a# p#] (~m a# (parametrization ~k p#)))
           ([k# a#] (~m a# (parametrization ~k)))
           ([k#] (~m 1.0 (parametrization ~k))))
         (register-field ~t ~k))))

;; Locally used random function for some parametrization parameters. Mostly used to avoid `0` value.
(defn- sdrand
  "Symetric random from [-mx -mn] and [mn mx]"
  ^double  [^double mn ^double mx]
  (let [rand (drand mn mx)]
    (randval rand (* -1.0 rand))))

(defn- sirand
  "Symmetric irand"
  (^long [^long mx] (sirand 1 mx))
  (^long [^long mn ^long mx]
   (let [rand (irand mn mx)]
     (randval rand (* -1 rand)))))

;; Local noise generating functions (for various noise implementations)

(defn- make-noise-variation
"Calculate shift by angle taken from noise"
[^double amount scale noise-fn]
(fn [v]
  (let [^Vec2 vv (v/mult v scale)
        angle (* m/TWO_PI ^double (noise-fn (.x vv) (.y vv)))]
    (v/add v (Vec2. (* amount (m/cos angle))
                    (* amount (m/sin angle)))))))

(defn- make-noise-variation2
  "Calculate shift by vector taken from noise"
  [^double amount scale noise-fn]
  (fn [v]
    (let [^Vec2 vv (v/mult v scale)
          x1 (- ^double (noise-fn (.x vv) (.y vv)) 0.5)
          y1 (- ^double (noise-fn (.y vv) m/E (.x vv)) 0.5)]
      (v/add v (Vec2. (* amount x1)
                      (* amount y1))))))

;;;;;

(defmacro ^:private load-fields-from-namespace
  [letter]
  (let [ns-sym (symbol (str "fastmath.fields." (name letter)))]
    (require ns-sym)
    `(do ~@(for [[sym v] (ns-publics ns-sym)
                 :let [{:keys [type config]} (v)
                       k (keyword sym)
                       f (symbol (str ns-sym "/" sym))
                       cs (gensym "config")]]
             `(let [~cs ~(when config `(:config (~f)))]
                ~(let [cfg (if config `(~cs) `{})]
                   `(do
                      (defmethod parametrization ~k
                        ([k# p#] (merge ~cfg p#))
                        ([k#] ~cfg))
                      (defmethod field ~k
                        ([k# a# p#] (~f a# (merge ~cfg p#)))
                        ([k# a#] (~f a# ~cfg))
                        ([k#] (~f 1.0 ~cfg)))
                      (register-field ~type ~k))))))))

(defmacro ^:private load-fields
  [letters]
  `(do ~@(for [l (map (comp keyword str) letters)]
           `(load-fields-from-namespace ~l))))

(load-fields "abcdef")

;; ## G

;; ### Gamma

(defn- make-gamma
  "gamma by zephyrtronium, http://fractal-resources.deviantart.com/art/Gamma-Apophysis-Plugin-154808483"
  [^double amount _]
  (fn [^Vec2 v]
    (Vec2. (* amount (m/log-gamma (v/mag v)))
           (* amount (v/heading v)))))
(make-field-method gamma :regular)

;; ### GaussianBlur

(defn- make-gaussianblur
  "Gaussian"
  [^double amount _]
  (fn [_]
    (let [a (drand m/TWO_PI)
          r (* amount (+ (drand) (drand) (drand) (drand) -2.0))]
      (Vec2. (* r (m/cos a)) (* r (m/sin a))))))
(make-field-method gaussianblur :random)

;; ### GlynnSim

(make-config-method glynnsim1 {:radius (r/drand 0.01 2.0)
                               :radius1 (sdrand 0.01 2.0)
                               :phi (r/drand 0.0 360.0)
                               :thickness (r/drand)
                               :pow (r/drand 6.0)
                               :contrast (r/drand)})

(defn- make-glynnsim1
  [^double amount {:keys [^double radius ^double radius1 ^double phi ^double thickness ^double pow ^double contrast]}]
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
              (v/mult toolpoint amount))))))))
(make-field-method glynnsim1 :random)

(make-config-method glynnsim2 {:radius (r/drand 0.01 2.0)
                               :phi1 (r/drand 0.0 360.0)
                               :phi2 (r/drand 0.0 360.0)
                               :thickness (r/drand)
                               :pow (r/drand 6.0)
                               :contrast (r/drand)})

(defn- make-glynnsim2
  [^double amount {:keys [^double radius ^double phi1 ^double phi2 ^double thickness ^double pow ^double contrast]}]
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
            (v/mult v (* alpha alpha))))))))
(make-field-method glynnsim2 :random)

(make-config-method glynnsim3 {:radius (r/drand 0.01 2.0)
                               :thickness (r/drand)
                               :pow (r/drand 6.0)
                               :contrast (r/drand)})

(defn- make-glynnsim3
  [^double amount {:keys [^double radius ^double thickness ^double pow ^double contrast]}]
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
            (v/mult v (* alpha alpha))))))))
(make-field-method glynnsim3 :random)

;; ### GDOffs

(make-config-method gdoffs {:delta-x (drand -6 6)
                            :delta-y (drand -6 6)
                            :area-x (sdrand 0.5 3)
                            :area-y (sdrand 0.5 3)
                            :center-x (drand -1 1)
                            :center-y (drand -1 1)
                            :gamma (drand -5 5)
                            :square (brand)})

(def ^:const ^:private ^double agdod-- 0.1)
(def ^:const ^:private ^double agdoa-- 2.0)
(def ^:const ^:private ^double agdoc-- 1.0)

(defn- fclp ^double [^double a] (if (neg? a) (- (rem (m/abs a) 1.0)) (rem (m/abs a) 1.0)))
(defn- fscl ^double [^double a] (fclp (* 0.5 (inc a))))
(defn- fosc ^double [^double p ^double a] (fscl (- (m/cos (* p a m/TWO_PI)))))
(defn- flip ^double [^double a ^double b ^double c] (+ a (* c (- b a))))

(defn- make-gdoffs
  "GDOffs"
  [^double amount {:keys [^double delta-x ^double delta-y ^double area-x ^double area-y ^double center-x ^double center-y ^double gamma square]}]
  (let [gdodx (* delta-x agdod--)
        gdody (* delta-y agdod--)
        gdoax (* agdoa-- (if (< (m/abs area-x) 0.1) 0.1 (m/abs area-x)))
        gdoay (* agdoa-- (if (< (m/abs area-y) 0.1) 0.1 (m/abs area-y)))
        gdocx (* center-x agdoc--)
        gdocy (* center-y agdoc--)
        gdog gamma
        gdos square
        gdob (/ (* gdog agdoa--) (max gdoax gdoay))]
    (fn [^Vec2 v]
      (let [osc-x (fosc gdodx 1.0)
            osc-y (if gdos (fosc gdody 1.0) 1.0)
            in-x (+ (.x v) gdocx)
            in-y (+ (.y v) gdocy)]
        (v/mult (if gdos
                  (Vec2. (flip (flip in-x (fosc in-x 4.0) osc-x) (fosc (fclp (* gdob in-x)) 4.0) osc-x)
                         (flip (flip in-y (fosc in-y 4.0) osc-x) (fosc (fclp (* gdob in-y)) 4.0) osc-x))
                  (Vec2. (flip (flip in-x (fosc in-x 4.0) osc-x) (fosc (fclp (* gdob in-x)) 4.0) osc-x)
                         (flip (flip in-y (fosc in-y 4.0) osc-y) (fosc (fclp (* gdob in-y)) 4.0) osc-y))) amount)))))
(make-field-method gdoffs :regular)

;; ## H
;;
;; ### Heart

(defn- make-heart
  "Heart"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (v/mag v)
          theta (v/heading v)
          rt (* r theta)
          sr (m/sin rt)
          cr (m/cos rt)]
      (Vec2. (* amount r sr) (- (* amount r cr))))))
(make-field-method heart :regular)

;; ### Handkerchief

(defn- make-handkerchief
  "Handkerchief"
  [^double amount _]
  (fn [^Vec2 v]
    (let [angle (v/heading v)
          r (v/mag v)]
      (Vec2. (* amount (* r (m/sin (+ angle r))))
             (* amount (* r (m/cos (- angle r))))))))
(make-field-method handkerchief :regular)

;; ### Hemisphere

(defn- make-hemisphere
  "Hemisphere"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (/ amount (m/sqrt (inc (v/magsq v))))]
      (Vec2. (* r (.x v))
             (* r (.y v))))))
(make-field-method hemisphere :regular)

;; ### Hole2

(make-config-method hole2 {:a (drand -2.0 2.0)
                           :b (drand -2.0 2.0)
                           :c (drand 3.0)
                           :d (drand -2.0 2.0)
                           :inside (brand)
                           :shape (irand 10)})

(defn- make-hole2
  "Hole2"
  [^double amount {:keys [^double a ^double b ^double c ^double d inside shape]}]
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
             (* r1 (m/sin theta))))))
(make-field-method hole2 :regular)

;; ### Horseshoe

(defn- make-horseshoe
  "Horseshoe"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (+ m/EPSILON (v/mag v))
          sina (/ (.x v) r)
          cosa (/ (.y v) r)]
      (Vec2. (* amount (- (* sina (.x v)) (* cosa (.y v))))
             (* amount (+ (* cosa (.x v)) (* sina (.y v))))))))
(make-field-method horseshoe :regular)

;; ### Hyperbolic

(defn- make-hyperbolic
  "Hyperbolic"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (+ m/EPSILON (v/mag v))
          theta (v/heading v)]
      (Vec2. (* amount (/ (m/sin theta) r))
             (* amount (m/cos theta) r)))))
(make-field-method hyperbolic :regular)

;;; ### Hypershift

(make-config-method hypershift {:shift (drand -2.0 2.0)
                                :stretch (drand -2.0 2.0)})

(defn- make-hypershift
  "Hypershift"
  [^double amount {:keys [^double shift ^double stretch]}]
  (let [scale (- 1.0 (* shift shift))]
    (fn [^Vec2 v]
      (let [rad (/ (v/magsq v))
            x (+ shift (* rad (.x v)))
            y (* rad (.y v))
            r (/ (* amount scale) (+ (* x x) (* y y)))]
        (Vec2. (+ shift (* r x))
               (* r y stretch))))))
(make-field-method hypershift :regular)

;; ## I

;; ### Circle inverse
;; https://mathworld.wolfram.com/Inversion.html

(make-config-method circle-inverse {:x0 (randval 0.3 0.0 (randval (irand -3 3) (drand -3.0 3.0)))
                                    :y0 (randval 0.3 0.0 (randval (irand -3 3) (drand -3.0 3.0)))
                                    :r (drand 0.01 3.0)})

(defn- make-circle-inverse
  "Circle inverse"
  [^double amount {:keys [^double x0 ^double y0 ^double r]}]
  (let [o (Vec2. x0 y0)
        r2 (* r r)]
    (fn [^Vec2 v]
      (let [diff (v/sub v o)]
        (v/mult (v/add o (v/div (v/mult diff r2) (v/magsq diff))) amount)))))
(make-field-method circle-inverse :regular)

;; ### InvTree

(defn- make-invtree
  "InvTree"
  [^double amount _]
  (fn [^Vec2 v]
    (cond
      (brand 0.333) (v/mult v (* 0.5 amount))
      (brand 0.666) (v/mult (Vec2. (/ (inc (.x v)))
                                   (/ (.y v) (inc (.y v)))) amount)
      :else (v/mult (Vec2. (/ (.x v) (inc (.x v)))
                           (/ (inc (.y v)))) amount))))
(make-field-method invtree :random)

;; ## J

;; ### Julia

(defn- make-julia
  "Julia"
  [^double amount _]
  (fn [^Vec2 v]
    (let [a (+ (* 0.5 (v/heading v)) (* m/PI (lrand 2)))
          r (* amount (m/sqrt (v/mag v)))]
      (Vec2. (* r (m/cos a)) (* r (m/sin a))))))
(make-field-method julia :random)

(defn- make-julia2
  "Julia with different angle calc"
  [^double amount _]
  (fn [^Vec2 v]
    (let [a (+ (* 0.5 (m/atan2 (.x v) (.y v))) (* m/PI (lrand 2)))
          r (* amount (m/sqrt (v/mag v)))]
      (Vec2. (* r (m/cos a)) (* r (m/sin a))))))
(make-field-method julia2 :random)

;; ### JuliaC

(make-config-method juliac {:re (int (sdrand 1.0 10.0))
                            :im (* 0.01 (drand -2.0 2.0))
                            :dist (drand -2.0 2.0)})

(defn- make-juliac
  "JuliaC"
  [^double amount {:keys [^double re ^double im ^double dist]}]
  (let [rre (/ 1.0 re)]
    (fn [^Vec2 v]
      (let [arg (+ (v/heading v)
                   (* m/TWO_PI ^double (mod (irand) re)))
            lnmod (* dist (m/log (v/magsq v)))
            a (+ (* arg rre)
                 (* lnmod im))
            mod2 (* amount (m/exp (- (* lnmod rre)
                                     (* arg im))))]
        (Vec2. (* mod2 (m/cos a))
               (* mod2 (m/sin a)))))))
(make-field-method juliac :random)

;; ### JuliaN

(make-config-method julian (let [r (sdrand 1 10)]
                             {:power (randval r (int r))
                              :dist (drand -4 4)}))

(defn- make-julian
  "JuliaN"
  [^double amount {:keys [^double power ^double dist]}]
  (let [abspower (int (m/abs power))
        cpower (* 0.5 (/ dist power))]
    (fn [v]
      (let [a (/ (+ (v/heading v) (* m/TWO_PI (double (irand abspower)))) power)
            r (* amount (m/pow (v/magsq v) cpower))]
        (Vec2. (* r (m/cos a)) (* r (m/sin a)))))))
(make-field-method julian :random)

;; ### JuliaScope

(make-config-method juliascope (let [r (sdrand 1 10)]
                                 {:power (randval r (int r))
                                  :dist (drand -4 4)}))

(defn- make-juliascope
  "JuliaScope"
  [^double amount {:keys [^double power ^double dist]}]
  (let [abspower (int (m/abs power))
        cpower (* 0.5 (/ dist power))]
    (fn [v]
      (let [rnd (double (lrand abspower))
            a (if (zero? (bit-and rnd 1))
                (/ (+ (* m/TWO_PI rnd) (v/heading v)) power)
                (/ (- (* m/TWO_PI rnd) (v/heading v)) power))
            r (* amount (m/pow (v/magsq v) cpower))]
        (Vec2. (* r (m/cos a)) (* r (m/sin a)))))))
(make-field-method juliascope :random)


;; ### JuliaQ

(make-config-method juliaq {:power (int (sdrand 1 10))
                            :divisor (sdrand 1 8)})

(defn- make-juliaq
  "juliaq by Zueuk, http://zueuk.deviantart.com/art/juliaq-Apophysis-plugins-340813357" 
  [^double amount {:keys [^double divisor ^double power]}]
  (let [inv-power (/ divisor power)
        half-inv-power (* 0.5 inv-power)
        inv-power-2pi (/ m/TWO_PI power)]
    (fn [^Vec2 v]
      (let [a (+ (* inv-power (v/heading v))
                 (* inv-power-2pi (double (irand))))
            r (* amount (m/pow (v/magsq v) half-inv-power))]
        (Vec2. (* r (m/cos a)) (* r (m/sin a)))))))
(make-field-method juliaq :random)

;; ## L

;; ### Lazy Susan

(make-config-method lazysusan {:twist (drand -6 6)
                               :spin (drand -4 4)
                               :space (drand -2 2)
                               :x (drand -1 1)
                               :y (drand -1 1)})

(defn- make-lazysusan
  "Lazysusan"
  [^double amount {:keys [^double twist ^double spin ^double space ^double x ^double y]}]
  (fn [^Vec2 v]
    (let [xx (- (.x v) x)
          yy (+ (.y v) y)
          rr (m/hypot xx yy)]
      (if (< rr amount)
        (let [a (+ (m/atan2 yy xx) spin (* twist (- amount rr)))
              nr (* amount rr)]
          (Vec2. (+ (* nr (m/cos a)) x)
                 (- (* nr (m/sin a)) y)))
        (let [nr (* amount (inc (/ space rr)))]
          (Vec2. (+ (* nr xx) x)
                 (- (* nr yy) y)))))))
(make-field-method lazysusan :regular)

;; ### LogApo

(make-config-method logapo {:base (drand 0.01 20)})

(defn- make-logapo
  "LogApo"
  [^double amount {:keys [^double base]}]
  (let [denom (/ 0.5 (m/log base))]
    (fn [v]
      (Vec2. (* amount denom (m/log (v/magsq v)))
             (* amount (v/heading v))))))
(make-field-method logapo :regular)

;; ### Log

(defn- make-log
  "Log"
  [^double amount _]
  (fn [^Vec2 v]
    (Vec2. (* amount 0.5 (m/log (v/magsq v)))
           (* amount (v/heading v)))))
(make-field-method log :regular)

;; ### Loonie

(defn- make-loonie
  "Loonie"
  [^double amount _]
  (let [w2 (m/sq amount)]
    (fn [v]
      (let [r2 (v/magsq v)]
        (if (and (< r2 w2) (not (zero? r2)))
          (let [r (* amount (m/sqrt (dec (/ w2 r2))))]
            (v/mult v r))
          (v/mult v amount))))))
(make-field-method loonie :regular)

;; ## M

;; ### Modulus

(make-config-method modulus {:x (sdrand 0.01 2.0)
                             :y (sdrand 0.01 2.0)})

(defn- make-modulus
  "Modulus"
  [amount {:keys [^double x ^double y]}]
  (let [xr (+ x x)
        yr (+ y y)]
    (fn [^Vec2 v]
      (v/mult (Vec2. (cond
                       (> (.x v) x) (+ (- x) (rem (+ (.x v) x) xr))
                       (< (.x v) (- x)) (- x (rem (- x (.x v)) xr))
                       :else (.x v))
                     (cond
                       (> (.y v) y) (+ (- y) (rem (+ (.y v) y) yr))
                       (< (.y v) (- y)) (- y (rem (- y (.y v)) yr))
                       :else (.y v))) amount))))
(make-field-method modulus :regular)


;; ## N

;; ### Ngon

(make-config-method ngon {:circle (drand -2 2)
                          :corners (drand -2 2)
                          :power (drand -10 10)
                          :sides (drand -10 10)})

(defn- make-ngon
  "Ngon"
  [^double amount {:keys [^double circle ^double corners ^double power ^double sides]}]
  (let [b (/ m/TWO_PI sides)
        hb (/ b 2.0)
        hpower (/ power 2.0)]
    (fn [v]
      (let [r-factor (m/pow (v/magsq v) hpower)
            theta (v/heading v)
            phi (- theta (* b (m/floor (/ theta b))))
            phi (if (> phi hb) (- phi b) phi)
            amp (/ (+ circle (* corners (dec (/ 1.0 (+ (m/cos phi) m/EPSILON))))) (+ r-factor m/EPSILON))]
        (v/mult v (* amount amp))))))
(make-field-method ngon :regular)

;; ### Noise

(defn- make-noise
  "Noise"
  [^double amount _]
  (fn [_]
    (let [a (drand m/TWO_PI)
          r (drand amount)]
      (Vec2. (* r (m/cos a))
             (* r (m/sin a))))))
(make-field-method noise :random)

;; ## O

;; ## Octapol

(make-config-method octapol {:polarweight (r/randval 0.2 0.0 (r/drand 0.01 2.0))
                             :radius (r/randval 0.2 0.0 (r/drand 0.2 2.0))
                             :s (r/drand 0.01 1.0)
                             :t (r/drand 0.01 1.0)
                             :scale (sdrand 0.05 0.5)})

(defn- octapol-hits-circle-around-origin
  ^double [^double radius ^Vec2 p]
  (if (zero? radius) 0.0 (v/mag p)))

(defn- octapol-hits-square-around-origin
  [^double a ^Vec2 XY]
  (and (<= (m/abs (.x XY)) a)
       (<= (m/abs (.y XY)) a)))

(defn- octapol-hits-rect
  [^Vec2 t1 ^Vec2 br ^Vec2 p]
  (and (>= (.x p) (.x t1))
       (>= (.y p) (.y t1))
       (<= (.x p) (.x br))
       (<= (.y p) (.y br))))

(defn- octapol-hits-triangle
  [^Vec2 a ^Vec2 b ^Vec2 c ^Vec2 p]
  (let [v0 (v/sub c a)
        v1 (v/sub b a)
        v2 (v/sub p a)
        d00 (v/dot v0 v0)
        d01 (v/dot v0 v1)
        d02 (v/dot v0 v2)
        d11 (v/dot v1 v1)
        d12 (v/dot v1 v2)
        denom (- (* d00 d11) (* d01 d01))
        ^Vec2 uv (if (zero? denom)
                   (Vec2. 0.0 0.0)
                   (v/div (Vec2. (- (* d11 d02) (* d01 d12))
                                 (- (* d00 d12) (* d01 d02))) denom))]
    (and (< (v/sum uv) 1.0) (pos? (.x uv)) (pos? (.y uv)))))



(defn- make-octapol
  [^double amount {:keys [^double polarweight ^double radius ^double s ^double t ^double scale]}]
  (let [-hs (* -0.5 s)
        hs (* 0.5 s)
        a (+ hs t)
        b (- -hs t)
        rad (* (m/abs radius) (/ s m/SQRT2))
        A (Vec2. -hs a)
        B (Vec2. hs a)
        C (Vec2. t hs)
        D (Vec2. t -hs)
        E (Vec2. hs b)
        F (Vec2. -hs b)
        G (Vec2. (- t) -hs)
        H (Vec2. (- t) hs)
        I (Vec2. -hs hs)
        J (Vec2. hs hs)
        K (Vec2. -hs -hs)
        L (Vec2. hs -hs)]
    (fn [^Vec2 v]
      (let [^Vec2 XY (v/mult v scale)
            r (octapol-hits-circle-around-origin rad XY)]
        (cond
          (and (pos? rad) (<= r rad)) (let [rd (m/log (m/sq (/ r rad)))
                                            phi (v/heading XY)
                                            t (* rd polarweight)]
                                        (v/mult (v/add XY (Vec2. (m/lerp (.x XY) phi t)
                                                                 (m/lerp (.y XY) r t))) amount))
          (and (octapol-hits-square-around-origin a XY)
               (or (octapol-hits-rect H K XY)
                   (octapol-hits-rect J D XY)
                   (octapol-hits-rect A J XY)
                   (octapol-hits-rect K E XY)
                   (octapol-hits-triangle I A H XY)
                   (octapol-hits-triangle J B C XY)
                   (octapol-hits-triangle L D E XY)
                   (octapol-hits-triangle K F G XY))) (v/mult (v/add XY XY) amount)
          :else (v/mult v amount))))))
(make-field-method octapol :regular)

;; ### Oscilloscope

(make-config-method oscilloscope {:separation (sdrand 0.05 2.0)
                                  :frequency (sdrand 0.01 10.0)
                                  :amplitude (sdrand 0.1 2.0)
                                  :damping (r/randval 0.5 0.0 (m/sqrt (r/drand 0.01 0.99)))})

(defn- make-oscilloscope
  [^double amount {:keys [^double separation ^double frequency ^double amplitude ^double damping]}]
  (let [tpf (* m/TWO_PI frequency)]
    (fn [^Vec2 v]
      (let [t (if (zero? damping)
                (+ separation (* amplitude (m/cos (* tpf (.x v)))))
                (+ separation (* amplitude (* (m/exp (* damping (- (m/abs (.x v)))))
                                              (m/cos (* tpf (.x v)))))))]
        (if (<= (m/abs (.y v)) t)
          (Vec2. (* amount (.x v))
                 (* -1.0 amount (.y v)))
          (v/mult v amount))))))
(make-field-method oscilloscope :regular)

(make-config-method oscilloscope2 {:separation (sdrand 0.05 2.0)
                                   :frequency-x (sdrand 0.01 10.0)
                                   :frequency-y (sdrand 0.01 10.0)
                                   :amplitude (sdrand 0.1 2.0)
                                   :perturbation (sdrand 0.1 2.0)
                                   :damping (r/randval 0.5 0.0 (m/sqrt (r/drand 0.01 0.99)))})

(defn- make-oscilloscope2
  [^double amount {:keys [^double separation ^double frequency-x ^double frequency-y ^double perturbation ^double amplitude ^double damping]}]
  (let [tpf (* m/TWO_PI frequency-x)
        tpf2 (* m/TWO_PI frequency-y)]
    (fn [^Vec2 v]
      (let [pt (* perturbation (m/sin (* tpf2 (.y v))))
            t (if (zero? damping)
                (+ separation (* amplitude (m/cos (+ pt (* tpf (.x v))))))
                (+ separation (* amplitude (* (m/exp (* damping (- (m/abs (.x v)))))
                                              (m/cos (+ pt (* tpf (.x v))))))))]
        (if (<= (m/abs (.y v)) t)
          (v/mult v (- amount))
          (v/mult v amount))))))
(make-field-method oscilloscope2 :regular)

;; ### Ovoid

(make-config-method ovoid {:x (sdrand 0.01 m/PI)
                           :y (sdrand 0.01 m/PI)})

(defn- make-ovoid
  [^double amount {:keys [^double x ^double y]}]
  (fn [^Vec2 v]
    (let [r (/ amount (+ m/EPSILON (v/magsq v)))]
      (v/emult (v/mult v r) (Vec2. x y)))))
(make-field-method ovoid :regular)

;; ## P

;; ### Panorama1

(defn- make-panorama1
  "Panorama1"
  [^double amount _]
  (fn [v]
    (let [aux (/ (m/sqrt (inc (v/magsq v))))
          nv (v/mult v aux)
          aux (v/mag nv)]
      (Vec2. (* amount m/M_1_PI (v/heading nv))
             (* amount (- aux 0.5))))))
(make-field-method panorama1 :regular)

;; ### Panorama2

(defn- make-panorama2
  "Panorama2"
  [^double amount _]
  (fn [v]
    (let [aux (/ (inc (m/sqrt (v/magsq v))))
          nv (v/mult v aux)
          aux (v/mag nv)]
      (Vec2. (* amount m/M_1_PI (v/heading nv))
             (* amount (- aux 0.5))))))
(make-field-method panorama2 :regular)

;; ### Parabola

(make-config-method parabola {:width (sdrand 0.5 2.0)
                              :height (sdrand 0.5 2.0)})

(defn- make-parabola
  "Parabola fn"
  [^double amount {:keys [^double width ^double height]}]
  (fn [v]
    (let [r (v/mag v)
          sr (m/sin r)
          cr (m/cos r)]
      (Vec2. (* amount height sr sr (drand))
             (* amount width cr (drand))))))
(make-field-method parabola :random)


;; ### Perlin

(make-config-method perlin {:seed (irand)
                            :scale (sdrand 0.1 1.5)
                            :octaves (irand 1 6)})

(defn- make-perlin
  "Perlin noise"
  [amount {:keys [seed octaves scale]}]
  (let [n (fbm-noise {:seed seed :octaves octaves})]
    (make-noise-variation amount scale n)))
(make-field-method perlin :regular)

(make-config-method perlin2 {:seed (irand)
                             :scale (sdrand 0.1 1.5)
                             :octaves (irand 1 6)})

(defn- make-perlin2
  "Perlin noise"
  [^double amount {:keys [^int seed ^int octaves ^double scale]}]
  (let [n (fbm-noise {:seed seed :octaves octaves})]
    (make-noise-variation2 amount scale n)))
(make-field-method perlin2 :regular)

;; 

(make-config-method general-noise (assoc (r/random-noise-cfg) :scale (sdrand 0.1 1.5)))

(defn- make-general-noise
  "Perlin noise"
  [amount cfg]
  (let [n (r/random-noise-fn cfg)]
    (make-noise-variation amount (:scale cfg) n)))
(make-field-method general-noise :regular)

(make-config-method general-noise2 (assoc (r/random-noise-cfg) :scale (sdrand 0.1 1.5)))

(defn- make-general-noise2
  "Perlin noise"
  [amount cfg]
  (let [n (r/random-noise-fn cfg)]
    (make-noise-variation2 amount (:scale cfg) n)))
(make-field-method general-noise2 :regular)

;; ### Petal

(defn- make-petal
  "Petal"
  [^double amount _]
  (fn [^Vec2 v]
    (let [a (m/cos (.x v))
          bx (m/pow (* (m/cos (.x v)) (m/cos (.y v))) 3.0)
          by (m/pow (* (m/sin (.x v)) (m/cos (.y v))) 3.0)]
      (Vec2. (* amount a bx)
             (* amount a by)))))
(make-field-method petal :regular)

;; ### Pie

(make-config-method pie {:slices (sdrand 0.01 7.0)
                         :rotation (drand m/TWO_PI)
                         :thickness (drand -2.0 2.0)})

(defn- make-pie
  "pie from jwildfire"
  [^double amount {:keys [^double slices ^double rotation ^double thickness]}]
  (fn [_]
    (let [sl (double (m/round (+ 0.5 (* slices (drand)))))
          a (-> thickness
                (* (drand))
                (+ sl)
                (* m/TWO_PI)
                (/ slices)
                (+ rotation))
          r (* amount (drand))]
      (Vec2. (* r (m/cos a))
             (* r (m/sin a))))))
(make-field-method pie :random)

;; ### PixelFlow

(make-config-method pixelflow {:angle (drand m/TWO_PI)
                               :len (sdrand 0.01 3.0)
                               :width (sdrand 0.1 10.0)
                               :seed (irand)})

(defn- pixelflow-hash
  ^double [^long a]
  (as-> (unchecked-int a) a
    (bit-xor (bit-xor a 61)
             (bit-shift-right a 16))
    (+ a (bit-shift-left a 3))
    (bit-xor a (bit-shift-right a 4))
    (unchecked-int (* a 0x27d4eb2d))
    (bit-xor a (bit-shift-right a 15))
    (/ (double (bit-and a 0xffffffff)) Integer/MAX_VALUE)))

(defn- make-pixelflow
  "Pixel Flow"
  [^double amount {:keys [^double angle ^double len ^double width ^long seed]}]
  (let [vin (Vec2. (m/cos angle) (m/sin angle))]
    (fn [^Vec2 v]
      (let [blockx (int (m/floor (* width (.x v))))
            blockx (+ blockx (- 2.0 (* 4.0 (pixelflow-hash (inc (* seed blockx))))))
            blocky (int (m/floor (* width (.y v))))
            blocky (+ blocky (- 2.0 (* 4.0 (pixelflow-hash (inc (* seed blocky))))))
            flen (* 0.5 (+ (pixelflow-hash (+ blocky (* blockx (- seed))))
                           (pixelflow-hash (+ blockx (* blocky (/ seed 2))))))]
        (v/add v (v/mult vin (* amount flen (m/sq (m/sq (drand))) len)))))))
(make-field-method pixelflow :random)

;; ### PDJ

(make-config-method pdj {:a (drand -6.0 6.0)
                         :b (drand -6.0 6.0)
                         :c (drand -6.0 6.0)
                         :d (drand -6.0 6.0)})

(defn- make-pdj
  "PDJ"
  [^double amount {:keys [^double a ^double b ^double c ^double d]}]
  (fn [^Vec2 v]
    (Vec2. (* amount (- (m/sin (* a (.y v))) (m/cos (* b (.x v)))))
           (* amount (- (m/sin (* c (.x v))) (m/cos (* d (.y v))))))))
(make-field-method pdj :regular)

;; ### Perspective

(make-config-method perspective {:angle (drand (- m/PI) m/PI)
                                 :dist (drand -5.0 5.0)})

(defn- make-perspective
  "Perspective"
  [^double amount {:keys [^double angle ^double dist]}]
  (let [ang (* m/HALF_PI angle)
        vsin (m/sin ang)
        vfcos (* dist (m/cos ang))]
    (fn [^Vec2 v]
      (let [t (/ amount (- dist (* (.y v) vsin)))]
        (Vec2. (* t dist (.x v))
               (* t vfcos (.y v)))))))
(make-field-method perspective :regular)

;; ### Phoenix julia

(make-config-method phoenix-julia {:power (int (sdrand 0.51 10))
                                   :dist (drand -2.0 2.0)
                                   :x-distort (drand -2.0 2.0)
                                   :y-distort (drand -2.0 2.0)})

(defn- make-phoenix-julia
  "Phoenix julia"
  [^double amount {:keys [^double power ^double dist ^double x-distort ^double y-distort]}]
  (let [inv-n (/ dist power)
        inv2pi-n (/ m/TWO_PI power)
        c-n (* 0.5 inv-n)]
    (fn [^Vec2 v]
      (let [pre-x (* (.x v) x-distort)
            pre-y (* (.y v) y-distort)
            a (+ (* (m/atan2 pre-y pre-x) inv-n)
                 (* (irand) inv2pi-n))
            sina (m/sin a)
            cosa (m/cos a)
            r (* amount (m/pow (v/magsq v) c-n))]
        (Vec2. (* r cosa)
               (* r sina))))))
(make-field-method phoenix-julia :random)

;; ### Polar

(defn- make-polar
  "Polar"
  [^double amount _]
  (fn [^Vec2 v]
    (let [ny (dec (v/mag v))]
      (Vec2. (* amount (v/heading v) m/M_1_PI)
             (* amount ny)))))
(make-field-method polar :regular)

(defn- make-polar2
  "Polar2"
  [^double amount _]
  (let [p2v (/ amount m/PI)
        p2v2 (* 0.5 p2v)]
    (fn [^Vec2 v] (Vec2. (* p2v (v/heading v)) (* p2v2 (m/log (v/magsq v)))))))
(make-field-method polar2 :regular)

;; ### PowBlock

(make-config-method powblock {:numerator (drand -20 20)
                              :denominator (drand -20 20)
                              :root (drand -6 6)
                              :correctn (drand -2 2)
                              :correctd (drand -2 2)})

(defn- make-powblock
  "PowBlock"
  [^double amount {:keys [^double numerator ^double denominator ^double root ^double correctn ^double correctd]}]
  (let [power (/ (* denominator correctn) (+ m/EPSILON (m/abs correctd)))
        power (if (< (m/abs power) m/EPSILON) m/EPSILON power)
        power (/ (* 0.5 numerator) power) 
        deneps (/ (if (< (m/abs denominator) m/EPSILON) m/EPSILON denominator))]
    (fn [^Vec2 v]
      (let [theta (v/heading v)
            r2 (* amount (m/pow (v/magsq v) power))
            ran (+ (* numerator (+ (* theta deneps) (* root m/TWO_PI (m/floor (drand denominator)) deneps))))]
        (Vec2. (* r2 (m/cos ran))
               (* r2 (m/sin ran)))))))
(make-field-method powblock :random)

;; ### Power

(defn- make-power
  "Power"
  [^double amount _]
  (fn [^Vec2 v]
    (let [theta (v/heading v)
          sa (m/sin theta)
          ca (m/cos theta)
          pow (* amount (m/pow (v/mag v) sa))]
      (Vec2. (* pow ca) (* pow sa)))))
(make-field-method power :regular)

;; ### Popcorn2

(make-config-method popcorn2 {:x (drand -1.5 1.5)
                              :y (drand -1.5 1.5)
                              :c (drand -5.0 5.0)})

(defn- make-popcorn2
  "popcorn2 from apophysis"
  [^double amount {:keys [^double x ^double y ^double c]}]
  (fn [^Vec2 v]
    (let [xx (->> (.y v)
                  (* c)
                  (m/tan)
                  (m/sin)
                  (* x)
                  (+ (.x v))
                  (* amount))
          yy (->> (.x v)
                  (* c)
                  (m/tan)
                  (m/sin)
                  (* y)
                  (+ (.y v))
                  (* amount))]
      (Vec2. xx yy))))
(make-field-method popcorn2 :regular)

;; ### Pressure Wave

(make-config-method pressure-wave {:x-freq (drand -6 6)
                                   :y-freq (drand -6 6)})

(defn- make-pressure-wave
  "Pressure Wave"
  [^double amount {:keys [^double x-freq ^double y-freq]}]
  (let [[^double pwx ^double ipwx] (if (zero? x-freq) [1.0 1.0] (let [pwx (* x-freq m/TWO_PI)] [pwx (/ pwx)]))
        [^double pwy ^double ipwy] (if (zero? y-freq) [1.0 1.0] (let [pwy (* y-freq m/TWO_PI)] [pwy (/ pwy)]))]
    (fn [^Vec2 v]
      (Vec2. (* amount (+ (.x v) (* ipwx (m/sin (* pwx (.x v))))))
             (* amount (+ (.y v) (* ipwy (m/sin (* pwy (.y v))))))))))
(make-field-method pressure-wave :regular)

;; ## R

;; ### R Circle Blur

(make-config-method r-circleblur {:n (drand -3 3)
                                  :seed (drand Float/MAX_VALUE)
                                  :dist (drand -1 1)
                                  :mn (drand -3 3)
                                  :mx (drand -3 3)})

(defn- make-r-circleblur
  "R Circle Blur"
  [^double amount {:keys [^double n ^double seed ^double dist ^double mn ^double mx]}]
  (let [dm (- mx mn)]
    (fn [v]
      (let [angle (v/heading v)
            rad (v/mag v)
            rad (mod rad n)
            by (m/sin (+ angle rad))
            bx (m/cos (+ angle rad))
            by (m/round (* by rad))
            bx (m/round (* bx rad))
            rad2 (* 0.5 (m/sqrt (drand)))
            angle2 (drand m/TWO_PI)
            a1 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 127.1) (* by 311.7) seed))))
            a2 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 269.5) (* by 183.3) seed))))
            a3 (m/sfrac (* 43758.5453 (m/sin (+ (* bx 78.233) (* by 12.9898) seed))))
            a3 (+ mn (* a3 dm))
            rad2 (* rad2 a3)]
        (Vec2. (* amount (+ bx (* rad2 (m/cos angle2)) (* dist a1)))
               (* amount (+ by (* rad2 (m/sin angle2)) (* dist a2))))))))
(make-field-method r-circleblur :random)


;; ### Radial Blur

(make-config-method radialblur {:angle (drand (- m/TWO_PI) m/TWO_PI)})

(defn- make-radialblur
  "Radial blur"
  [^double amount {:keys [^double angle]}]
  (let [spin (* amount (m/sin (* angle m/HALF_PI)))
        zoom (* amount (m/cos (* angle m/HALF_PI)))]
    (fn [^Vec2 v]
      (let [rnd-g (+ (drand) (drand) (drand) (drand) -2.0)
            ra (v/mag v)
            alpha (+ (* spin rnd-g) (v/heading v))
            rz (dec (* zoom rnd-g))]
        (Vec2. (+ (* rz (.x v)) (* ra (m/cos alpha)))
               (+ (* rz (.y v)) (* ra (m/sin alpha))))))))
(make-field-method radialblur :random)

;; ### Rational3

(make-config-method rational3 {:a (drand -3 3)
                               :b (drand -3 3)
                               :c (drand -3 3)
                               :d (drand -3 3)
                               :e (drand -3 3)
                               :f (drand -3 3)
                               :g (drand -3 3)
                               :h (drand -3 3)})

(defn- make-rational3
  "Rational3"
  [^double amount {:keys [^double a ^double b ^double c ^double d
                          ^double e ^double f ^double g ^double h]}]
  (fn [^Vec2 v]
    (let [^Vec2 sqr (v/emult v v)
          ^Vec2 cb (v/emult sqr v)
          zt3 (- (.x cb) (* 3.0 (.x v) (.y sqr)))
          zt2 (- (.x sqr) (.y sqr))
          zb3 (- (* 3.0 (.x sqr) (.y v)) (.y cb))
          zb2 (* 2.0 (.x v) (.y v))
          tr (+ (* a zt3) (* b zt2) (* c (.x v)) d)
          ti (+ (* a zb3) (* b zb2) (* c (.y v )))
          br (+ (* e zt3) (* f zt2) (* g (.x v)) h)
          bi (+ (* e zb3) (* f zb2) (* g (.y v )))
          r3den (/ amount (+ (* br br) (* bi bi)))]
      (Vec2. (* r3den (+ (* tr br) (* ti bi)))
             (* r3den (- (* ti br) (* tr bi)))))))
(make-field-method rational3 :regular)

;; ### Rays

(defn- make-rays
  "Rays"
  [^double amount _]
  (fn [^Vec2 v]
    (let [ang (* amount (drand m/PI))
          r (/ amount (+ m/EPSILON (v/magsq v)))
          tanr (* amount r (m/tan ang))]
      (Vec2. (* tanr (m/cos (.x v)))
             (* tanr (m/sin (.y v)))))))
(make-field-method rays :random)

(defn- make-rays1
  "Rays1"
  [^double amount _]
  (let [pa (* amount (m/sq m/M_2_PI))]
    (fn [^Vec2 v]
      (let [t (v/magsq v)
            u (+ pa (/ (m/tan (m/sqrt t))))
            r (* amount u t)]
        (Vec2. (/ r (.x v))
               (/ r (.y v)))))))
(make-field-method rays1 :regular)

(defn- make-rays2
  "Rays2"
  [^double amount _]
  (let [a10 (/ amount 10.0)]
    (fn [^Vec2 v]
      (let [t (v/magsq v)
            u (/ (m/cos (* (+ t m/EPSILON) (m/tan (/ (+ m/EPSILON t))))))
            r (* a10 t u)]
        (Vec2. (/ r (.x v))
               (/ r (.y v)))))))
(make-field-method rays2 :regular)

(defn- make-rays3
  "Rays3"
  [^double amount _]
  (let [a10 (/ amount 10.0)]
    (fn [^Vec2 v]
      (let [t (v/magsq v)
            t2 (* t t)
            u (/ (m/sqrt (m/cos (m/sin (* (+ m/EPSILON t2) (m/sin (/ (+ t2 m/EPSILON))))))))
            r (* a10 t u)]
        (Vec2. (/ (* r (m/cos t)) (.x v))
               (/ (* r (m/tan t)) (.y v)))))))
(make-field-method rays3 :regular)

;; ### Rectangles

(make-config-method rectangles {:x (drand -1.5 1.5)
                                :y (drand -1.5 1.5)})

(defn- make-rectangles
  "Rectangles"
  [^double amount {:keys [^double x ^double y]}]
  (fn [^Vec2 v]
    (Vec2. (if (< (m/abs (.x v)) m/EPSILON)
             (* amount (.x v))
             (* amount (-> (.x v)
                           (/ x)
                           (m/floor)
                           (* 2.0)
                           inc
                           (* x)
                           (- (.x v)))))
           (if (< (m/abs (.y v)) m/EPSILON)
             (* amount (.y v))
             (* amount (-> (.y v)
                           (/ y)
                           (m/floor)
                           (* 2.0)
                           inc
                           (* y)
                           (- (.y v))))))))
(make-field-method rectangles :regular)

;; ### Rhodonea

(make-config-method rhodonea {:knumer (randval (int (sdrand 1 30)) (sdrand 1 30))
                              :kdenom (randval (int (sdrand 1 30)) (sdrand 1 30))
                              :radial-offset (drand -1 1)
                              :inner-mode (irand 7)
                              :outer-mode (irand 7)
                              :inner-spread (drand -1 1)
                              :outer-spread (drand -1 1)
                              :inner-spread-ratio (drand -2 2)
                              :outer-spread-ratio (drand -2 2)
                              :spread-split (drand -1.5 1.5)
                              :cycle-offset (drand m/TWO_PI)
                              :cycles-param (randval 0 (drand 100))
                              :metacycle-expansion (drand -1 1)
                              :metacycles (drand 10)
                              :fill (randval 0 (drand))})


(defn- make-rhodonea
  "Rhodonea"
  [^double amount {:keys [^double knumer ^double kdenom ^double radial-offset ^long inner-mode ^long outer-mode
                          ^double inner-spread ^double outer-spread ^double inner-spread-ratio ^double outer-spread-ratio
                          ^double spread-split ^double cycles-param ^double cycle-offset ^double metacycle-expansion
                          ^double metacycles ^double fill]}]
  (let [kn knumer
        kd kdenom
        k (/ kn kd)
        cycles-to-close (double (if (zero? (mod k 1.0))
                                  (if (zero? (mod k 2.0))
                                    1.0
                                    (if (or (not (zero? radial-offset))
                                            (not (zero? inner-spread))
                                            (not (zero? outer-spread))
                                            (not (zero? fill)))
                                      1.0
                                      0.5))
                                  (if (and (zero? (mod kn 1.0))
                                           (zero? (mod kd 1.0)))
                                    (let [lkn (long kn)
                                          lkd (long kd)
                                          gcd (m/gcd lkn lkd)
                                          [^long kn ^long kd] (if (not== gcd 1)
                                                                [(/ lkn gcd) (/ lkd gcd)]
                                                                [lkn lkd])]
                                      (if (or (zero? (mod kn 2.0))
                                              (zero? (mod kd 2.0)))
                                        kd
                                        (/ kd 2)))
                                    (if (< cycles-param 16)
                                      16
                                      (* 2 kd kn)))))
        cycles (if (zero? cycles-param)
                 (* cycles-to-close metacycles)
                 cycles-param)]
    (fn [^Vec2 v]
      (let [rin (* spread-split (v/mag v))
            tin (v/heading v)
            t (* cycles (+ tin (* cycle-offset m/TWO_PI)))
            r (+ radial-offset (m/cos (* t k)))
            r (if-not (zero? fill)
                (+ r (* fill (- (drand) 0.5)))
                r)
            x (* r (m/cos t))
            y (* r (m/sin t))
            expansion (m/floor (/ (* cycles (+ tin m/PI))
                                  (* cycles-to-close m/TWO_PI)))
            adjusted-amount (+ amount (* expansion metacycle-expansion))]
        (if (> (m/abs rin) (m/abs r))
          (case (long outer-mode)
            0 (Vec2. (* adjusted-amount x)
                     (* adjusted-amount y))
            1 (let [rinx (inc (* (dec rin) outer-spread outer-spread-ratio))
                    riny (inc (* (dec rin) outer-spread))]
                (Vec2. (* adjusted-amount rinx x)
                       (* adjusted-amount riny y)))
            2 (let [xin (* (m/sgn x) (m/abs (.x v)))
                    yin (* (m/sgn y) (m/abs (.y v)))]
                (Vec2. (* adjusted-amount (+ x (* outer-spread outer-spread-ratio (- xin x))))
                       (* adjusted-amount (+ y (* outer-spread (- yin y))))))
            3 (let [xin (* (m/sgn x) (m/abs (.x v)))
                    yin (* (m/sgn y) (m/abs (.y v)))]
                (Vec2. (* adjusted-amount (+ x (* outer-spread outer-spread-ratio xin)))
                       (* adjusted-amount (+ y (* outer-spread yin)))))
            4 (let [rinx (+ (* 0.5 rin) (* outer-spread outer-spread-ratio))
                    riny (+ (* 0.5 rin) outer-spread)]
                (Vec2. (* adjusted-amount rinx x)
                       (* adjusted-amount riny y)))
            5 v
            6 (Vec2. 0.0 0.0))
          (case (long inner-mode)
            0 (Vec2. (* adjusted-amount x)
                     (* adjusted-amount y))
            1 (let [rinx (inc (* (dec rin) inner-spread inner-spread-ratio))
                    riny (inc (* (dec rin) inner-spread))]
                (Vec2. (* adjusted-amount rinx x)
                       (* adjusted-amount riny y)))
            2 (let [xin (* (m/sgn x) (m/abs (.x v)))
                    yin (* (m/sgn y) (m/abs (.y v)))]
                (Vec2. (* adjusted-amount (+ x (* inner-spread inner-spread-ratio (- xin x))))
                       (* adjusted-amount (+ y (* inner-spread (- yin y))))))
            3 (let [xin (* (m/sgn x) (m/abs (.x v)))
                    yin (* (m/sgn y) (m/abs (.y v)))]
                (Vec2. (* adjusted-amount (+ x (* inner-spread inner-spread-ratio xin)))
                       (* adjusted-amount (+ y (* inner-spread yin)))))
            4 (let [rinx (+ (* 0.5 rin) (* inner-spread inner-spread-ratio))
                    riny (+ (* 0.5 rin) inner-spread)]
                (Vec2. (* adjusted-amount rinx x)
                       (* adjusted-amount riny y)))
            5 v
            6 (Vec2. 0.0 0.0)
            ))))))
(make-field-method rhodonea :random)

;; ### Rings

(make-config-method rings {:coeff20 (drand 1.3)})

(defn- make-rings
  "Rings"
  [^double amount {:keys [^double coeff20]}]
  (let [dx (+ m/EPSILON (m/sq coeff20))
        dx2 (+ dx dx)
        rdx (/ 1.0 dx2)
        dx- (- 1.0 dx)]
    (fn [^Vec2 v]
      (let [r (v/mag v)
            rr (* amount (+ (- r (* dx2 (double (int (* (+ r dx) rdx))))) (* r dx-)))]
        (Vec2. (* rr (/ (.x v) r))
               (* rr (/ (.y v) r)))))))
(make-field-method rings :regular)

(make-config-method rings2 {:val (drand -1.0 1.0)})

(defn- make-rings2
  "Rings2"
  [^double amount {:keys [^double val]}]
  (let [dx (+ m/EPSILON (m/sq val))]
    (fn [^Vec2 v]
      (let [l (v/mag v)
            r (* amount (- 2.0 (* dx (inc (/ (* 2.0 (double (int (* 0.5 (inc (/ l dx)))))) l)))))]
        (v/mult v r)))))
(make-field-method rings2 :regular)

;; ### Ripple

(make-config-method ripple {:frequency (drand -3 3)
                            :velocity (drand -3 3)
                            :amplitude (drand -2 2)
                            :centerx (drand -0.5 0.5)
                            :centery (drand -0.5 0.5)
                            :phase (drand -1 1)
                            :scale (sdrand 0.5 4)
                            :fixed-dist-calc (irand 4)})

(defn- make-ripple
  "Ripple"
  [^double amount {:keys [^double frequency ^double velocity ^double amplitude ^double centerx
                          ^double centery ^double phase ^double scale ^int fixed-dist-calc]}]
  (let [f (* frequency 5.0)
        a (* amplitude 0.01)
        p (- (* phase m/TWO_PI) m/PI)
        s (if (zero? scale) m/EPSILON scale)
        is (/ s)
        vxp (* velocity p)
        pxa (* p a)
        pixa (* (- m/PI p) a)]
    (fn [^Vec2 v]
      (let [x (- (* s (.x v)) centerx)
            y (- (* s (.y v)) centery)
            d (case (unchecked-int fixed-dist-calc)
                0 (v/mag v)
                1 (m/sqrt (* (m/sq (.x v)) (m/sq (.y v))))
                2 (max (m/abs (.x v)) (m/abs (.y v)))
                3 (+ (m/abs (.x v)) (m/abs (.y v))))
            d (if (< d m/EPSILON) m/EPSILON d)
            nx (/ x d)
            ny (/ y d)
            wave (m/cos (- (* f d) vxp))
            d1 (+ d (* wave pxa))
            d2 (+ d (* wave pixa))
            u1 (+ centerx (* nx d1))
            u2 (+ centerx (* nx d2))
            v1 (- (* ny d1) centery)
            v2 (- (* ny d2) centery)]
        (Vec2. (* amount is (m/lerp u1 u2 p))
               (* amount is (m/lerp v1 v2 p)))))))
(make-field-method ripple :regular)

;; ### Rippled

(defn- make-rippled
  "Rippled"
  [^double amount _]
  (fn [^Vec2 v]
    (let [d (+ m/EPSILON (v/magsq v))]
      (Vec2. (* (* (m/tanh d) (* 2.0 (.x v))) (/ amount 2.0))
             (* (* (m/cos d) (* 2.0 (.y v))) (/ amount 2.0))))))
(make-field-method rippled :regular)

;; ### Round Sphere

(defn- make-roundspher
  ""
  [^double amount _]
  (let [s (m/sq m/M_2_PI)]
    (fn [v]
      (let [d (v/magsq v)
            re (/ (+ s (/ d)))
            ad (/ amount d)]
        (v/mult v (* amount ad re))))))
(make-field-method roundspher :regular)

;; ## S

;; ### Scry

(defn- make-scry
  "Scry"
  [^double amount _]
  (fn [^Vec2 v]
    (let [t (v/magsq v)
          d (-> 1.0
                (/ amount)
                (+ t)
                (* (m/sqrt t))
                (+ m/EPSILON))
          r (/ 1.0 d)]
      (v/mult v r))))
(make-field-method scry :regular)

;; ### Sech

(defn- make-sech
  "Sech"
  [^double amount _]
  (fn [^Vec2 v]
    (let [sn (m/sin (.y v))
          cn (m/cos (.y v))
          snh (m/sinh (.x v))
          cnh (m/cosh (.x v))
          d (+ (m/cos (* 2.0 (.y v)))
               (m/cosh (* 2.0 (.x v))))
          d (if (zero? d) m/EPSILON d)
          den (/ 2.0 d)]
      (Vec2. (* amount den cn cnh)
             (* (- amount) den sn snh)))))
(make-field-method sech :regular)

;; ### Shreadrad

(make-config-method shreadrad {:n (randval (int (sdrand 1 9)) (sdrand 0.0001 8.0))
                               :width (drand -2.0 2.0)})

(defn- make-shreadrad
  "ShreadRad"
  [^double amount {:keys [^double n ^double width]}]
  (let [sa (/ m/TWO_PI n)
        sa2 (* 0.5 sa)
        sa2sw (* sa2 width)
        pi3 (* 3.0 m/PI)]
    (fn [v]
      (let [ang (v/heading v)
            rad (v/mag v)
            xang (/ (+ ang pi3 sa2) sa)
            ixang (unchecked-int xang)
            zang (- (* sa (+ ixang (* width (- xang ixang)))) m/PI sa2sw)]
        (Vec2. (* amount rad (m/cos zang))
               (* amount rad (m/sin zang)))))))
(make-field-method shreadrad :regular)

;; ### Sinusoidal

(defn- make-sinusoidal
  "Sinusoidal"
  [^double amount _]
  (fn [^Vec2 v]
    (Vec2. (* amount (m/sin (.x v))) (* amount (m/sin (.y v))))))
(make-field-method sinusoidal :regular)

;; ### Secant

(defn- make-secant
  "Secant2"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (* amount (v/mag v))
          cr (* amount (m/cos r))
          icr (/ 1.0 (if (zero? cr) m/EPSILON cr))]
      (Vec2. (* amount (.x v)) icr))))
(make-field-method secant :regular)

(defn- make-secant2
  "Secant2"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (* amount (v/mag v))
          cr (m/cos r)
          icr (/ 1.0 (if (zero? cr) m/EPSILON cr))
          ny (if (neg? cr)
               (* amount (inc icr))
               (* amount (dec icr)))]
      (Vec2. (* amount (.x v)) ny))))
(make-field-method secant2 :regular)

;; ### Spherical

(defn- make-spherical
  "Spherical"
  [^double amount _]
  (fn [^Vec2 v]
    (v/mult v (/ amount (+ m/EPSILON (v/magsq v))))))
(make-field-method spherical :regular)

;; ### Spiral

(defn- make-spiral
  "Spiral"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (+ m/EPSILON (v/mag v))
          revr (/ 1.0 r)
          sina (* (.x v) revr)
          cosa (* (.y v) revr)
          sinr (m/sin r)
          cosr (m/cos r)]
      (Vec2. (* amount revr (+ cosa sinr))
             (* amount revr (- sina cosr))))))
(make-field-method spiral :regular)

;; ### Split

(make-config-method split {:xsplit (drand (- m/TWO_PI) m/TWO_PI)
                           :ysplit (drand (- m/TWO_PI) m/TWO_PI)})

(defn- make-split
  "Split"
  [^double amount {:keys [^double xsplit ^double ysplit]}]
  (fn [^Vec2 v]
    (Vec2. (if (pos? (m/cos (* (.x v) xsplit)))
             (* amount (.y v))
             (- (* amount (.y v))))
           (if (pos? (m/cos (* (.y v) ysplit)))
             (* amount (.x v))
             (- (* amount (.x v)))))))
(make-field-method split :regular)

;; ### Splits

(make-config-method splits {:x (drand -1.5 1.5)
                            :y (drand -1.5 1.5)})

(defn- make-splits
  "Splits"
  [^double amount {:keys [^double x ^double y]}]
  (fn [^Vec2 v]
    (Vec2. (if (pos? (.x v))
             (* amount (+ (.x v) x))
             (* amount (- (.x v) x)))
           (if (pos? (.y v))
             (* amount (+ (.y v) y))
             (* amount (- (.y v) y))))))
(make-field-method splits :regular)

;; ### Square

(defn- make-square
  "Square"
  [^double amount _]
  (fn [_]
    (Vec2. (* amount (drand -0.5 0.5))
           (* amount (drand -0.5 0.5)))))
(make-field-method square :random)

;; ### Squirrel

(make-config-method squirrel {:a (drand m/EPSILON 4.0)
                              :b (drand m/EPSILON 4.0)})

(defn- make-squirrel
  "Squirrel"
  [^double amount {:keys [^double a ^double b]}]
  (fn [^Vec2 v]
    (let [u (m/sqrt (+ (* a (m/sq (.x v)))
                       (* b (m/sq (.y v)))))]
      (Vec2. (* amount (m/cos u) (m/tan (.x v)))
             (* amount (m/sin u) (m/tan (.y v)))))))
(make-field-method squirrel :regular)

;; ### STwin

(make-config-method stwin {:distort (drand -6 6)
                           :multiplier (sdrand 0.001 3.0)})

(defn- make-stwin
  "STwin by Xyrus-02, http://timothy-vincent.deviantart.com/art/STwin-Plugin-136504836"
  [^double amount {:keys [^double distort ^double multiplier]}]
  (fn [^Vec2 v]
    (let [x (* (.x v) amount multiplier)
          y (* (.y v) amount multiplier)
          x2 (* x x)
          y2 (* y y)
          x2+y2 (+ x2 y2)
          x2-y2 (- x2 y2)
          div (if (zero? x2+y2) 1.0 x2+y2)
          result (/ (* x2-y2 (m/sin (* m/TWO_PI distort (+ x y)))) div)]
      (Vec2. (+ (* amount (.x v)) result)
             (+ (* amount (.y v)) result)))))
(make-field-method stwin :regular)

;; ### Supershape

(make-config-method supershape {:rnd (drand -1 1)
                                :m (drand m/TWO_PI)
                                :n1 (drand -5 5)
                                :n2 (drand -5 5)
                                :n3 (drand -5 5)
                                :holes (drand -1 1)})

(defn- make-supershape
  "Supershape"
  [^double amount {:keys [^double rnd ^double m ^double n1 ^double n2 ^double n3 ^double holes]}]
  (let [pm-4 (/ m 4.0)
        pneg1-n1 (/ -1.0 n1)]
    (fn [^Vec2 v]
      (let [theta (+ (* pm-4 (v/heading v)) m/M_PI_4)
            st (m/sin theta)
            ct (m/cos theta)
            t1 (m/pow (m/abs ct) n2)
            t2 (m/pow (m/abs st) n3)
            mag (v/mag v)
            r (/ (* (* amount (- (+ (drand rnd) (* (- 1.0 rnd) mag)) holes)) (m/pow (+ t1 t2) pneg1-n1)) mag)]
        (v/mult v r)))))
(make-field-method supershape :random)

;; ### Swirl

(defn- make-swirl
  "Swirl"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (v/magsq v)
          c1 (m/sin r)
          c2 (m/cos r)]
      (Vec2. (* amount (- (* c1 (.x v)) (* c2 (.y v))))
             (* amount (+ (* c2 (.x v)) (* c1 (.y v))))))))
(make-field-method swirl :regular)

;; ## T

;; ### Tangent

(defn- make-tangent
  "Tangent"
  [^double amount _]
  (fn [^Vec2 v]
    (let [d (m/cos (.y v))
          id (/ 1.0 (if (zero? d) m/EPSILON d))]
      (Vec2. (* amount (m/sin (.x v)) id)
             (* amount (m/tan (.y v)))))))
(make-field-method tangent :regular)

;; ### Twintrian

(defn- make-twintrian
  "Twintrian"
  [^double amount _]
  (fn [^Vec2 v]
    (let [r (* amount (drand) (v/mag v))
          sinr (m/sin r)
          diff (+ (m/cos r) (m/log10 (m/sq sinr)))]
      (Vec2. (* amount diff (.x v))
             (* amount (.x v) (- diff (* m/PI sinr)))))))
(make-field-method twintrian :random)


;; ### Taurus

(make-config-method taurus {:r (drand -5.0 5.0)
                            :n (drand -5.0 5.0)
                            :inv (drand -2.0 2.0)})

(defn- make-taurus
  "Taurus"
  [^double amount {:keys [^double r ^double n ^double inv]}]
  (let [rinv (* r inv)
        revinv (- 1.0 inv)]
    (fn [^Vec2 v]
      (let [sx (m/sin (.x v))
            cx (m/cos (.x v))
            sy (m/sin (.y v))
            ir (+ rinv (* revinv r (m/cos (* n (.x v)))))
            irsy (+ ir sy)]
        (Vec2. (* amount cx irsy)
               (* amount sx irsy))))))
(make-field-method taurus :regular)

;; ### Trade

(make-config-method trade {:r1 (drand 0.1 3.0)
                           :r2 (drand 0.1 3.0)
                           :d1 (drand -2.0 2.0)
                           :d2 (drand -2.0 2.0)})

(defn- make-trade
  "trade by Michael Faber,  http://michaelfaber.deviantart.com/art/The-Lost-Variations-258913970"
  [^double amount {:keys [^double r1 ^double r2 ^double d1 ^double d2]}]
  (let [c1 (+ r1 d1)
        c2 (+ r2 d2)]
    (fn [^Vec2 v]
      (let [[^double cc1 ^double cc2 ^double fr ^double rr] (if (pos? (.x v))
                                                              [c1 (- c2) (/ r2 r1) r1]
                                                              [(- c2) c1 (/ r1 r2) r2])
            nv (Vec2. (- cc1 (.x v)) (.y v))
            rm (v/mag nv)
            r (* rm fr)
            a (v/heading nv)
            res (Vec2. (+ cc2 (* r (m/cos a)))
                       (* r (m/sin a)))]
        (if (<= rm rr)
          (v/mult res amount)
          (v/mult v amount))))))
(make-field-method trade :regular)

;; ## V

;; ### Vibration

(make-config-method vibration {:dir (drand m/TWO_PI)
                               :angle (drand m/TWO_PI)
                               :freq (sdrand 0.01 2.0)
                               :amp (sdrand 0.1 1.0)
                               :phase (drand)
                               :dir2 (drand m/TWO_PI)
                               :angle2 (drand m/TWO_PI)
                               :freq2 (sdrand 0.01 2.0)
                               :amp2 (sdrand 0.1 1.0)
                               :phase2 (drand)})

(defn- make-vibration
  "Vibration http://fractal-resources.deviantart.com/art/Apo-Plugins-Vibration-1-and-2-252001851"
  [^double amount {:keys [^double dir ^double angle ^double freq ^double amp ^double phase
                          ^double dir2 ^double angle2 ^double freq2 ^double amp2 ^double phase2]}]
  (let [total-angle (+ angle dir)
        cos-dir (m/cos dir)
        sin-dir (m/sin dir)
        cos-tot (m/cos total-angle)
        sin-tot (m/sin total-angle)
        scaled-freq (* m/TWO_PI freq)
        phase-shift (/ (* m/TWO_PI phase) freq)
        total-angle2 (+ angle2 dir2)
        cos-dir2 (m/cos dir2)
        sin-dir2 (m/sin dir2)
        cos-tot2 (m/cos total-angle2)
        sin-tot2 (m/sin total-angle2)
        scaled-freq2 (* m/TWO_PI freq2)
        phase-shift2 (/ (* m/TWO_PI phase2) freq2)]
    (fn [^Vec2 v]
      (let [d-along-dir (+ (* (.x v) cos-dir)
                           (* (.y v) sin-dir))
            local-amp (* amp (m/sin (+ (* d-along-dir scaled-freq) phase-shift)))
            x (+ (.x v) (* local-amp cos-tot))
            y (+ (.y v) (* local-amp sin-tot))
            d-along-dir (+ (* (.x v) cos-dir2)
                           (* (.y v) sin-dir2))
            local-amp (* amp2 (m/sin (+ (* d-along-dir scaled-freq2) phase-shift2)))
            x (+ x (* local-amp cos-tot2))
            y (+ y (* local-amp sin-tot2))]
        (Vec2. (* amount x) (* amount y))))))
(make-field-method vibration :regular)

(make-config-method vibration2 {:dir (drand m/TWO_PI)
                                :angle (drand m/TWO_PI)
                                :freq (sdrand 0.01 2.0)
                                :amp (sdrand 0.1 1.0)
                                :phase (drand)
                                :dir2 (drand m/TWO_PI)
                                :angle2 (drand m/TWO_PI)
                                :freq2 (sdrand 0.01 2.0)
                                :amp2 (sdrand 0.1 1.0)
                                :phase2 (drand)
                                :dm (drand -0.5 0.5)
                                :dmfreq (sdrand 0.01 1.0)
                                :tm (drand -0.5 0.5)
                                :tmfreq (sdrand 0.01 1.0)
                                :fm (drand -0.5 0.5)
                                :fmfreq (sdrand 0.01 1.0)
                                :am (drand -0.5 0.5)
                                :amfreq (sdrand 0.01 1.0)
                                :d2m (drand -0.5 0.5)
                                :d2mfreq (sdrand 0.01 1.0)
                                :t2m (drand -0.5 0.5)
                                :t2mfreq (sdrand 0.01 1.0)
                                :f2m (drand -0.5 0.5)
                                :f2mfreq (sdrand 0.01 1.0)
                                :a2m (drand -0.5 0.5)
                                :a2mfreq (sdrand 0.01 1.0)})

(defn- v-modulate 
  "Modulate"
  ^double [^double amp ^double freq ^double x]
  (* amp (m/cos (* x freq m/TWO_PI))))

(defn- make-vibration2
  "Vibration 2 http://fractal-resources.deviantart.com/art/Apo-Plugins-Vibration-1-and-2-252001851"
  [^double amount {:keys [^double dir ^double angle ^double freq ^double amp ^double phase
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
      (let [d-along-dir (+ (* (.x v) cdir)
                           (* (.y v) sdir))
            dir-l (+ dir (v-modulate dm dmfreq d-along-dir))
            angle-l (+ angle (v-modulate tm tmfreq d-along-dir))
            freq-l (/ (v-modulate fm fmfreq d-along-dir) freq)
            amp-l (+ amp (* amp (v-modulate am amfreq d-along-dir)))
            total-angle (+ angle-l dir-l)
            cos-dir (m/cos dir-l)
            sin-dir (m/sin dir-l)
            cos-tot (m/cos total-angle)
            sin-tot (m/sin total-angle)
            scaled-freq (* m/TWO_PI freq)
            phase-shift (/ (* m/TWO_PI phase) freq)
            d-along-dir (+ (* (.x v) cos-dir)
                           (* (.y v) sin-dir))
            local-amp (* amp-l (m/sin (+ (* d-along-dir scaled-freq) freq-l phase-shift)))
            x (+ (.x v) (* local-amp cos-tot))
            y (+ (.y v) (* local-amp sin-tot))

            d-along-dir (+ (* (.x v) cdir2)
                           (* (.y v) sdir2))
            dir-l (+ dir2 (v-modulate d2m d2mfreq d-along-dir))
            angle-l (+ angle2 (v-modulate t2m t2mfreq d-along-dir))
            freq-l (/ (v-modulate f2m f2mfreq d-along-dir) freq2)
            amp-l (+ amp2 (* amp2 (v-modulate a2m a2mfreq d-along-dir)))
            total-angle (+ angle-l dir-l)
            cos-dir (m/cos dir-l)
            sin-dir (m/sin dir-l)
            cos-tot (m/cos total-angle)
            sin-tot (m/sin total-angle)
            scaled-freq (* m/TWO_PI freq2)
            phase-shift (/ (* m/TWO_PI phase2) freq2)
            d-along-dir (+ (* (.x v) cos-dir)
                           (* (.y v) sin-dir))
            local-amp (* amp-l (m/sin (+ (* d-along-dir scaled-freq) freq-l phase-shift)))
            x (+ x (* local-amp cos-tot))
            y (+ y (* local-amp sin-tot))]
        (Vec2. (* amount x)
               (* amount y))))))
(make-field-method vibration2 :regular)

;; ### Voron

(make-config-method voron {:k (sdrand 0.6 1.3)
                           :step (sdrand 0.1 1.2)
                           :num (drand 0.1 25.0)
                           :xseed (irand)
                           :yseed (irand)})

(deftype VoronResType [^double R ^double X0 ^double Y0])
(deftype VoronCalcType [^long M1 ^long N1 ^long k])

(defn- make-voron
  "Voron by eralex61, http://eralex61.deviantart.com/art/Voronoi-Diagram-plugin-153126702"
  [^double amount {:keys [^double k ^double step ^double num ^int xseed ^int yseed]}]
  (fn [^Vec2 v]
    (let [fk (fn ^VoronCalcType [^long M1 ^long N1]
               (VoronCalcType. M1 N1
                               (long (inc (m/floor (* (discrete-noise (+ (+ (* M1 19) (* N1 257)) xseed) 0) num))))))
          m (long (m/floor (/ (.x v) step)))
          n (long (m/floor (/ (.y v) step)))
          m- (dec m)
          m+ (inc m)
          n- (dec n)
          n+ (inc n)
          Ks (mapv fk [m- m- m- m m m m+ m+ m+] [n- n n+ n- n n+ n- n n+])
          ^VoronResType res (reduce (fn [^VoronResType curr ^VoronCalcType calc]
                                      (loop [i (long 0)
                                             ^VoronResType currl curr]
                                        (if (< i (.k calc))
                                          (let [X (* step (+ (double (.M1 calc)) (discrete-noise (+
                                                                                                  (+ i (* 64 (.M1 calc)))
                                                                                                  (+ xseed (* 15 (.N1 calc)))) 0)))
                                                Y (* step (+ (double (.N1 calc)) (discrete-noise (+
                                                                                                  (+ i (* 21 (.M1 calc)))
                                                                                                  (+ yseed (* 33 (.N1 calc)))) 0)))
                                                R (m/hypot (- (.x v) X) (- (.y v) Y))]
                                            (recur (unchecked-inc i)
                                                   (if (< R (.R currl))
                                                     (VoronResType. R X Y)
                                                     currl)))
                                          currl))) (VoronResType. 20.0 0.0 0.0) Ks)]
      (Vec2. (* amount (+ (.X0 res) (* k (- (.x v) (.X0 res)))))
             (* amount (+ (.Y0 res) (* k (- (.y v) (.Y0 res)))))))))
(make-field-method voron :regular)

;; ### Waves

(make-config-method waves {:coeff10 (drand -2.0 2.0)
                           :coeff11 (drand -2.0 2.0)
                           :coeff20 (drand -2.0 2.0)
                           :coeff21 (drand -2.0 2.0)})

(defn- make-waves
  "Waves"
  [^double amount {:keys [^double coeff10 ^double coeff11 ^double coeff20 ^double coeff21]}]
  (let [c202 (+ m/EPSILON (m/sq coeff20))
        c212 (+ m/EPSILON (m/sq coeff21))]
    (fn [^Vec2 v]
      (Vec2. (->> c202
                  (/ (.y v))
                  (m/sin)
                  (* coeff10)
                  (+ (.x v))
                  (* amount))
             (->> c212
                  (/ (.x v))
                  (m/sin)
                  (* coeff11)
                  (+ (.y v))
                  (* amount))))))
(make-field-method waves :regular)

;; ### Wedge

(make-config-method wedge {:angle (drand m/TWO_PI)
                           :hole (drand -2 2)
                           :count (drand -5 5)
                           :swirl (drand -2 2)})

(defn- make-wedge
  "Wedge"
  [^double amount {:keys [^double angle ^double hole ^double count ^double swirl]}]
  (let [hm1p (* m/M_1_PI 0.5)]
    (fn [v]
      (let [r (v/mag v)
            a (+ (v/heading v) (* r swirl))
            c (m/floor (* (+ (* count a) m/PI) hm1p))
            comp-fac (- 1.0 (* angle count hm1p))
            a (+ (* a comp-fac) (* c angle))
            r (* amount (+ r hole))]
        (Vec2. (* r (m/cos a))
               (* r (m/sin a)))))))
(make-field-method wedge :regular)


;; ## Additional variations
;;
;; https://github.com/d3/d3-geo-projection/tree/master/src

;; ### Miller

(defn- make-miller
  "Miller"
  [^double amount _]
  (fn [^Vec2 v]
    (v/mult (Vec2. (.x v)
                   (->> (m/constrain (.y v) -1.9634 1.9634)
                        (* 0.4)
                        (+ m/QUARTER_PI)
                        (m/tan)
                        (m/log)
                        (* 1.25))) amount)))
(make-field-method miller :regular)

(defn- make-millerrev
  "Millerrev"
  [^double amount _]
  (fn [^Vec2 v]
    (v/mult (Vec2. (.x v)
                   (-> (.y v)
                       (* 0.8)
                       (m/exp)
                       (m/atan)
                       (* 2.5)
                       (- (* 0.625 m/PI)))) amount)))
(make-field-method millerrev :regular)

;; ### Foucaut

(defn- make-foucaut
  "Foucaut"
  [^double amount _]
  (fn [^Vec2 v]
    (let [k (* 0.5 (.y v))
          cosk (m/cos k)
          xx (->> cosk
                  (* cosk)
                  (* (m/cos (.y v)))
                  (* (/ (.x v) m/SQRTPI))
                  (* 2.0)
                  (* amount))
          yy (* amount m/SQRTPI (m/tan k))]
      (Vec2. xx yy))))
(make-field-method foucaut :regular)

;; ## Lists

;; List of variations based on RNG
(def fields-list-random (concat (:random @fields-atom)
                              (:pattern @fields-atom)))

;; List of variations not random
(def fields-list-not-random (:regular @fields-atom))

;; list of all variations defined in the file
(def fields-list (concat fields-list-random fields-list-not-random))
(def fields-map @fields-atom)

(reset! fields-atom nil)

;; ## Function arithmetic
;;
(def ^{:dynamic true
     :doc "When random configuration for [[combine]] is used. Skip vector fields which are random."
     :metadoc/categories #{:vf}}
  *skip-random-fields* false)

(defn- directional-derivative
  "Compute directional derivative."
  [f dir ^double amount ^double h]
  (fn [v]
    (let [v1 (f v)
          v2 (f (v/add v dir))]
      (v/mult (v/div (v/sub v2 v1) h) amount))))

(defn derivative 
  "Calculate directional derivative of fn. Derivative is calculated along [1,1] vector with `h` as a step (default `1.0e-6`)."
  {:metadoc/categories #{:vf}}
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. h h) amount h))
  ([f ^double h]
   (derivative f 1.0 h))
  ([f]
   (derivative f 1.0 1.0e-6)))

(defn grad-x
  "Calculate gradient along x axis."
  {:metadoc/categories #{:vf}}
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. h 0.0) amount h))
  ([f ^double h]
   (grad-x f 1.0 h))
  ([f]
   (grad-x f 1.0 1.0e-6)))

(defn grad-y
  "Calculate gradient along y axis."
  {:metadoc/categories #{:vf}}
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. 0.0 h) amount h))
  ([f ^double h]
   (grad-y f 1.0 h))
  ([f]
   (grad-y f 1.0 1.0e-6)))

(defn jacobian
  "Det of Jacobian of the field"
  {:metadoc/categories #{:sc}}
  ([f] (jacobian f 1.0e-6))
  ([f ^double h]
   (let [gx (grad-x f h)
         gy (grad-y f h)]
     (fn [^Vec2 v]
       (let [^Vec2 gxv (gx v)
             ^Vec2 gyv (gy v)]
         (- (* (.x gxv) (.y gyv))
            (* (.x gyv) (.y gxv))))))))

(defn divergence
  "Divergence of the field.

  See: https://youtu.be/rB83DpBJQsE?t=855"
  {:metadoc/categories #{:sc}}
  ([f] (divergence f 1.0e-6))
  ([f ^double h]
   (let [gx (grad-x f h)
         gy (grad-y f h)]
     (fn ^double [v]
       (let [^Vec2 gxv (gx v)
             ^Vec2 gyv (gy v)]
         (+ (.x gxv) (.y gyv)))))))

(defn curl
  "Curl (2d version) of the field.

  See: https://youtu.be/rB83DpBJQsE?t=855"
  {:metadoc/categories #{:sc}}
  ([f] (curl f 1.0e-6))
  ([f ^double h]
   (let [gx (grad-x f h)
         gy (grad-y f h)]
     (fn ^double [v]
       (let [^Vec2 gxv (gx v)
             ^Vec2 gyv (gy v)]
         (- (.x gyv) (.y gxv)))))))

(defn magnitude
  "Magnitude of the vectors from field."
  {:metadoc/categories #{:sc}}
  [f]
  (fn ^double [v]
    (v/mag (f v))))

(defn heading
  "Angle of the vectors from field."
  {:metadoc/categories #{:sc}}
  [f]
  (fn ^double [v]
    (v/heading (f v))))

(defn- generate-scalar-field
  "Generate scalar field"
  [op]
  (fn
    ([f]
     (fn [v]
       (op (f v) v)))
    ([f1 f2]
     (fn [v]
       (op (f1 v) (f2 v))))))

(def ^{:doc "2d cross product (det of the 2x2 matrix) of the input vector and result of the vector field.

In case when two vector fields are given, cross product is taken from results of vector fields."
       :metadoc/categories #{:sc}}
  cross (generate-scalar-field v/cross))
(def ^{:doc "Dot product of the input vector and result of the vector field.

In case when two vector fields are given, cross product is taken from result of vector fields."
       :metadoc/categories #{:sc}}
  dot (generate-scalar-field v/dot))
(def ^{:doc "Angle between input vector and result of the vector field.

In case when two vector fields are given, cross product is taken from result of vector fields.

Resulting value is from range `[-PI,PI]`."
       :metadoc/categories #{:vf}}
  angle-between (generate-scalar-field v/angle-between))

(defn scalar->vector-field
  "Returns vector field build from scalar fields of the input vector and result of the vector field."
  {:metadoc/categories #{:vf}}
  ([scalar f] 
   (fn [v]
     (Vec2. (scalar (f v)) (scalar v))))
  ([scalar f1 f2]
   (fn [v]
     (Vec2. (scalar (f1 v)) (scalar (f2 v))))))

(defn composition
  "Compose two vector fields."
  {:metadoc/categories #{:vf}}
  ([f1 f2 ^double amount]
   (fn [v] (v/mult (f1 (f2 v)) amount)))
  ([f1 f2] (composition f1 f2 1.0)))

(defn sum
  "Add two vector fields."
  {:metadoc/categories #{:vf}}
  ([f1 f2 ^double amount]
   (fn [v] (v/mult (v/add (f1 v) (f2 v)) amount)))
  ([f1 f2] (sum f1 f2 1.0)))

(defn multiplication
  "Multiply two vector fields (as a element-wise multiplication of results)."
  {:metadoc/categories #{:vf}}
  ([f1 f2 ^double amount]
   (fn [v] (v/mult (v/emult (f1 v) (f2 v)) amount)))
  ([f1 f2] (multiplication f1 f2 1.0)))

(defn- build-random-variation-step
  "Create variation parametrization"
  []
  (let [n (rand-nth (if *skip-random-fields* fields-list-not-random fields-list))]
    {:type :variation :name n :amount 1.0 :config (parametrization n)}))

(defn- build-random-parametrization-step
  "Create parametrization tree"
  ([f1 f2]
   (let [operand (rand-nth [:comp :add :comp :add :comp :mult :comp :angles :comp])]
     {:type :operation :name operand :var1 f1 :var2 f2}))
  ([f]
   (randval 0.1 f
            (randval 0.1
                     {:type :operation :name :deriv :var f}
                     (build-random-parametrization-step f (build-random-variation-step)))))
  ([]
   (build-random-parametrization-step (build-random-variation-step) (build-random-variation-step))))

(defn randomize-configuration
  "Randomize values for given configuration. Keeps structure untouched."
  {:metadoc/categories #{:vf}}
  ([f]
   (if (= (:type f) :variation)
     (assoc f :amount 1.0 :config (parametrization (:name f) {}))
     (let [name (:name f)]
       (if (= name :deriv)
         (assoc f :amount 1.0 :step (m/sq (drand 0.01 1.0)) :var (randomize-configuration (:var f)))
         (let [amount1 (if (#{:comp :angles} name) 1.0 (drand -2.0 2.0))
               amount2 (if (#{:comp :angles} name) 1.0 (drand -2.0 2.0)) 
               amount (case name
                        :add (/ 1.0 (+ (m/abs amount1) (m/abs amount2)))
                        :mult (/ 1.0 (* amount1 amount2))
                        1.0)]
           (assoc f :amount amount
                  :var1 (assoc (randomize-configuration (:var1 f)) :amount amount1)
                  :var2 (assoc (randomize-configuration (:var2 f)) :amount amount2))))))))

(defn random-configuration
  "Create random configuration for [[combine]] function. Optionally with depth (0 = only root is created).

  See [[combine]] for structure.

  Bind `*skip-random-fields*` to true to exclude fields which are random."
  {:metadoc/categories #{:vf}}
  ([] (random-configuration (lrand 5)))
  ([depth] (random-configuration depth (build-random-variation-step)))
  ([^long depth f]
   (if (pos? depth)
     (random-configuration (dec depth) (randomize-configuration (build-random-parametrization-step f)))
     f)))

(defn combine
  "Create composite vector field function based on configuration

  Call without argument to get random vector field.

  Configuration is a tree structure where nodes are one of the following

  * `{:type :variation :name NAME :amount AMOUNT :config CONFIG}` where
      * NAME is variation name (keyword)
      * AMOUNT is scaling factor
      * CONFIG is variation parametrization
  * `{:type :operation :name OPERATION :amount AMOUNT :var1 VAR1 :var2 VAR2}` where
      * OPERATION is one of the operations (see below)
      * AMOUNT is scaling factor
      * VAR1 and VAR2 two variations to combine
  * `{:type :operation :name :derivative :amount AMOUNT :var VAR :step STEP}` where
      * AMOUNT is scaling factor
      * VAR variation, subject to calculate derivative
      * STEP dx and dy value

  Possible OPERATIONs are:

  * `:add` - sum of two variations
  * `:mult` - multiplication
  * `:comp` - composition
  * `:angles` - vector field from angles

  See [[random-configuration]] for example."
  {:metadoc/categories #{:vf}}
  ([{:keys [type name amount config var step var1 var2]}]
   (if (= type :variation)
     (field name amount config)
     (if (= name :deriv)
       (derivative (combine var) amount step)
       (let [v1 (combine var1)
             v2 (combine var2)]
         (case name
           :comp (composition v1 v2 amount)
           :add (sum v1 v2 amount)
           :mult (multiplication v1 v2 amount)
           :angles (scalar->vector-field v/heading v1 v2))))))
  ([] (combine (random-configuration))))

(defn random-field
  "Create randomized field (optional depth can be provided)."
  ([] (combine))
  ([depth] (combine (random-configuration depth))))
