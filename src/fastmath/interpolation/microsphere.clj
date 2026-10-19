(ns fastmath.interpolation.microsphere
  "Microsphere projection interpolation of scattered data (Apache Commons Math)."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.analysis MultivariateFunction]
           [org.apache.commons.math3.analysis.interpolation InterpolatingMicrosphere
            InterpolatingMicrosphere2D MicrosphereProjectionInterpolator]
           [org.apache.commons.math3.random RandomGenerator UnitSphereRandomVectorGenerator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- infer-elements
  ^long [x ^long cnt]
  (m/max (if (number? x) 5 (m/+ 4 (count x))) (m// cnt 5)))

(defn- check-option!
  [option value ok? expected]
  (when-not (and (number? value) (ok? (double value)))
    (ic/invalid-input! :microsphere-projection option
                       (str option " has to be " expected ", got " (pr-str value))
                       {:option option :value value :expected expected})))

(defn- ->point
  ^doubles [x]
  (m/seq->double-array x))

(defn microsphere-projection
  "Creates the microsphere projection interpolator of scattered data in any dimension.

  The method (the algorithm of Dudziak as implemented in Apache Commons Math) places a microsphere with facets in the query point and lights every facet by the nodes, with the brightness falling with the distance. The value is the sum of the illumination of the facets, weighted by the values at the nodes. It passes through the nodes (within `:no-interpolation-tolerance`), returns `:background` where the sphere is dark, and is usually less accurate than [[fastmath.interpolation.rbf/rbf]] on smooth data.

  Parameters:

  - `xss` (sequence): the nodes, at least one; numbers or sequences of one number for one dimensional data, sequences of numbers of an equal dimension otherwise
  - `ys` (sequence of numbers): the values at the nodes, as many as nodes
  - `opts` (map, optional; a `nil` value is the same as a missing key):
      - `:elements` (long, default: `max(5, 4 + dimension, n / 5)` with `n` the number of nodes): the number of facets of the microsphere; more facets give a smoother and more accurate result and a longer evaluation (the cost is proportional to the number of facets and nodes)
      - `:exponent` (double, default: `1.0`): the exponent of the brightness `1 / distance^exponent`, not negative; the original paper uses `2`, a larger value makes the result closer to the value of the nearest node
      - `:max-dark-friction` (double, default: `0.9`): the largest fraction of the facets which can be dark (not lit by any node) for the value to be computed; above it the `:background` is returned
      - `:dark-threshold` (double, default: `0.01`): the illumination below which a facet counts as dark
      - `:background` (double, default: `0.0`): the value returned where the sphere is too dark
      - `:no-interpolation-tolerance` (double, default: `1.0e-6`): the distance to a node below which the value of the node is returned
      - `:rng` (random number generator, default: a new unseeded generator for every call): the generator of the random facets for the dimensions other than 2, see [[fastmath.random/rng]]. Pass a seeded one to get the same interpolator on every construction (the generator is consumed); equal seeds give equal values and different seeds give different values. In two dimensions the facets are equally spaced angles (deterministic) and `:rng` is not used

  Returns a function of a query point (a number or a one element sequence for one dimensional data, a sequence of numbers otherwise; in two dimensions also two numbers `(f x y)`) returning a double. The function is safe to call from several threads (calls are serialised).

  The result depends on the facets: in dimensions other than 2 the sphere is random, so without a seeded `:rng` two constructions differ (by up to about 0.2 on a smooth test function), and the default number of facets is a compromise (a few hundred facets are more accurate).

  Throws `ex-info` (with `:type :invalid-input`, `:method :microsphere-projection` and the key `:property`) when `xss` and `ys` have different lengths (`:length-mismatch`), there are no nodes (`:too-few-points`), the nodes have different dimensions (`:dimension-mismatch`, also for a query of a wrong dimension), or `:elements` is below 1 or `:exponent` is negative or not finite.

  See also [[fastmath.interpolation/microsphere-projection]], [[fastmath.interpolation.rbf/rbf]], [[fastmath.interpolation.shepard/shepard]]."
  ([xss ys] (microsphere-projection xss ys nil))
  ([xss ys opts]
   (let [{:keys [elements exponent max-dark-friction dark-threshold background
                 no-interpolation-tolerance rng]
          :or {exponent 1.0
               max-dark-friction 0.9
               dark-threshold 0.01
               background 0.0
               no-interpolation-tolerance 1.0e-6}}
         (into {} (remove (fn [[_ v]] (nil? v))) opts)
         nx (count xss)
         ny (count ys)]
     (when-not (m/== nx ny)
       (ic/invalid-input! :microsphere-projection :length-mismatch
                          (str "xss and ys have different lengths: " nx " and " ny)
                          {:xs-count nx :ys-count ny}))
     (when (m/zero? nx)
       (ic/invalid-input! :microsphere-projection :too-few-points "at least 1 point is required, got 0"
                          {:axis :xss :min-points 1 :count 0}))
     (let [elements (long (or elements (infer-elements (first xss) nx)))
           _ (check-option! :elements elements (fn [^double x] (m/>= x 1.0)) "a number >= 1")
           _ (check-option! :exponent exponent (fn [^double x] (and (m/valid-double? x) (m/>= x 0.0)))
                            "a finite number >= 0")
           pts (mapv ->point xss)
           dims (alength ^doubles (first pts))
           _ (doseq [^doubles pt pts]
               (when-not (m/== dims (alength pt))
                 (ic/invalid-input! :microsphere-projection :dimension-mismatch
                                    (str "all the nodes have to have the same dimension, got " dims " and " (alength pt))
                                    {:expected dims :got (alength pt)})))
           ^InterpolatingMicrosphere sphere
           (if (m/== dims 2)
             (InterpolatingMicrosphere2D. elements max-dark-friction dark-threshold background)
             (InterpolatingMicrosphere. dims elements max-dark-friction dark-threshold background
                                        (UnitSphereRandomVectorGenerator. dims ^RandomGenerator (r/ensure-rng rng))))
           ^MultivariateFunction obj (-> (MicrosphereProjectionInterpolator. sphere exponent false no-interpolation-tolerance)
                                         (.interpolate (into-array (class (double-array 0)) pts)
                                                       (m/seq->double-array ys)))
           value (fn ^double [x]
                   (let [q (->point x)]
                     (when-not (m/== dims (alength q))
                       (ic/invalid-input! :microsphere-projection :dimension-mismatch
                                          (str "the query has the dimension " (alength q) ", the nodes have " dims)
                                          {:expected dims :got (alength q)}))
                     ;; the evaluation changes the state of the sphere
                     (locking obj (.value obj q))))]
       (if (m/== dims 2)
         (fn microsphere-interpolator
           (^double [xs] (value xs))
           (^double [^double x ^double y] (value [x y])))
         (fn ^double [x] (value x)))))))
