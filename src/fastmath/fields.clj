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

;;;;;

(defmacro ^:private load-fields-from-namespace
  [letter-char]
  (let [ns-sym (symbol (str "fastmath.fields." letter-char))]
    (require ns-sym)
    `(do ~@(for [[sym v] (ns-publics ns-sym)
                 :when (= (first (name sym)) letter-char)
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
  `(do ~@(for [l letters]
           `(load-fields-from-namespace ~l))))

(load-fields "abcdefghijklmnopqrst")


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
