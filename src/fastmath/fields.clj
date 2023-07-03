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
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Two atoms to store variation names. One for non-random functions and second for random.

(def ^:private fields-atom (atom {}))

(defn- register-field
  "Add `name` to the atom `what`"
  [type name]
  (swap! fields-atom update type conj name))

(defmulti parametrization
  "Return random parametrization map for given field.

  Optinally you can pass part of the parametrization. In this case function will add remaining keys with randomly generated values.

  If field doesn't have parametrization, empty map will be returned.

  See [[field]]."
  (fn [key & _] key))

(defmethod parametrization :default
  ([_ _] {})
  ([_] {}))

(defmulti field
  "Return vector field for given name and options: amount (scaling factor) and parametrization.

  Default scaling factor is 1.0, default parametrization is random.

  Resulting function operates on [[Vec2]] type."
  (fn [key & _] key))

#_(defmethod field :default
    ([_ s _] (fn [v] (v/mult v s)))
    ([_ s] (fn [v] (v/mult v s)))
    ([_] identity))

;;;;;

(defmacro ^:private load-fields-from-namespace
  [letter-char]
  (let [ns-sym (symbol (str "fastmath.fields." letter-char))]
    (require ns-sym)
    `(do (require (quote ~ns-sym))
         ~@(for [[sym v] (ns-publics ns-sym)
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

(load-fields "abcdefghijklmnopqrstuvwxyz")

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
     :doc "When random configuration for [[combine]] is used. Skip vector fields which are random."}
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
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. h h) amount h))
  ([f ^double h]
   (derivative f 1.0 h))
  ([f]
   (derivative f 1.0 1.0e-6)))

(defn grad-x
  "Calculate gradient along x axis."
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. h 0.0) amount h))
  ([f ^double h]
   (grad-x f 1.0 h))
  ([f]
   (grad-x f 1.0 1.0e-6)))

(defn grad-y
  "Calculate gradient along y axis."
  ([f ^double amount ^double h]
   (directional-derivative f (Vec2. 0.0 h) amount h))
  ([f ^double h]
   (grad-y f 1.0 h))
  ([f]
   (grad-y f 1.0 1.0e-6)))

(defn jacobian
  "Det of Jacobian of the field"
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
  [f]
  (fn ^double [v]
    (v/mag (f v))))

(defn heading
  "Angle of the vectors from field."
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

In case when two vector fields are given, cross product is taken from results of vector fields."}
  cross (generate-scalar-field v/cross))
(def ^{:doc "Dot product of the input vector and result of the vector field.

In case when two vector fields are given, cross product is taken from result of vector fields."}
  dot (generate-scalar-field v/dot))
(def ^{:doc "Angle between input vector and result of the vector field.

In case when two vector fields are given, cross product is taken from result of vector fields.

Resulting value is from range `[-PI,PI]`."}
  angle-between (generate-scalar-field v/angle-between))

(defn scalar->vector-field
  "Returns vector field build from scalar fields of the input vector and result of the vector field."
  ([scalar f]
   (fn [v]
     (Vec2. (scalar (f v)) (scalar v))))
  ([scalar f1 f2]
   (fn [v]
     (Vec2. (scalar (f1 v)) (scalar (f2 v))))))

(defn composition
  "Compose two vector fields."
  ([f1 f2 ^double amount]
   (fn [v] (v/mult (f1 (f2 v)) amount)))
  ([f1 f2] (composition f1 f2 1.0)))

(defn sum
  "Add two vector fields."
  ([f1 f2 ^double amount]
   (fn [v] (v/mult (v/add (f1 v) (f2 v)) amount)))
  ([f1 f2] (sum f1 f2 1.0)))

(defn multiplication
  "Multiply two vector fields (as a element-wise multiplication of results)."
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
   (r/randval 0.1 f
              (r/randval 0.1
                         {:type :operation :name :deriv :var f}
                         (build-random-parametrization-step f (build-random-variation-step)))))
  ([]
   (build-random-parametrization-step (build-random-variation-step) (build-random-variation-step))))

(defn randomize-configuration
  "Randomize values for given configuration. Keeps structure untouched."
  ([f]
   (if (= (:type f) :variation)
     (assoc f :amount 1.0 :config (parametrization (:name f) {}))
     (let [name (:name f)]
       (if (= name :deriv)
         (assoc f :amount 1.0 :step (m/sq (r/drand 0.01 1.0)) :var (randomize-configuration (:var f)))
         (let [amount1 (if (#{:comp :angles} name) 1.0 (r/drand -2.0 2.0))
               amount2 (if (#{:comp :angles} name) 1.0 (r/drand -2.0 2.0))
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
  ([] (random-configuration (r/lrand 5)))
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
