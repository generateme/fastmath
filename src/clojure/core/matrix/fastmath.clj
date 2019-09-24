(ns clojure.core.matrix.fastmath
  (:require [clojure.core.matrix.protocols :as mat]
            [fastmath.vector :refer [abs add as-vec clamp cross dist div dot emult fmap interpolate mag magsq make-vector mn mult mx nonzero-count normalize permute sub sum zero-count sq sigmoid log1p exp ceil floor atan cos sin log10 tan log cbrt sqrt cosh asin sinh acos round signum radians degrees tanh]]
            [fastmath.core :as m])
  (:import [fastmath.vector Vec2 Vec3 Vec4 ArrayVec]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- new-matrix-impl 
  "Return matrix as vector."
  [^long cols ^long rows]
  (cond
    (== 1 cols) (make-vector rows)
    (== 1 rows) (make-vector cols)
    :else nil))

(defmacro ^:private extend-protocol-for-types
  [p xs m]
  (let [mp (gensym "fnmap")]
    `(let [~mp ~m]
       ~@(for [t xs]
           `(extend ~t ~p ~mp)))))

(extend-protocol-for-types mat/PImplementation [Vec2 Vec3 Vec4 ArrayVec]
                           {:implementation-key (constantly :fastmath)
                            :make-vector (fn [_ d] (make-vector d))
                            :new-matrix (fn [_ r c] (new-matrix-impl r c))
                            :new-matrix-nd (fn [_ shape] (make-vector (first shape)))
                            :supports-dimensionality? (fn [_ ^long dims] (== dims 1))
                            :construct-matrix as-vec})

(extend-protocol mat/PImplementation
  Vec2 (mat/meta-info [v] {:doc (str "Fastmath 2d vector: " v)})
  Vec3 (mat/meta-info [v] {:doc (str "Fastmath 3d vector: " v)})
  Vec4 (mat/meta-info [v] {:doc (str "Fastmath 4d vector: " v)})
  ArrayVec (mat/meta-info [v] {:doc "Fastmath double array vector"}))

(extend-protocol-for-types mat/PDimensionInfo [Vec2 Vec3 Vec4 ArrayVec]
                           {:dimensionality (constantly 1)
                            :is-scalar? (constantly false)
                            :is-vector? (constantly true)
                            :dimension-count (fn [v ^long dim] (when (== 0 dim) (count v)))
                            :get-shape #(vector (count %))})

(extend-protocol-for-types mat/PIndexedAccess [Vec2 Vec3 Vec4 ArrayVec] 
                           {:get-1d (fn [m row] (m row))})

(extend-protocol-for-types mat/PIndexedSetting [Vec2 Vec3 Vec4] 
                           {:is-mutable? (constantly false)})

(extend-protocol mat/PIndexedSetting
  Vec2 (mat/set-1d [^Vec2 m ^long row ^double v] (case row
                                                   0 (Vec2. v (.y m))
                                                   1 (Vec2. (.x m) v)
                                                   m))
  Vec3 (mat/set-1d [^Vec3 m ^long row ^double v] (case row
                                                   0 (Vec3. v (.y m) (.z m))
                                                   1 (Vec3. (.x m) v (.z m))
                                                   2 (Vec3. (.x m) (.y m) v)
                                                   m))
  Vec4 (mat/set-1d [^Vec4 m ^long row ^double v] (case row
                                                   0 (Vec4. v (.y m) (.z m) (.w m))
                                                   1 (Vec4. (.x m) v (.z m) (.w m))
                                                   2 (Vec4. (.x m) (.y m) v (.w m)) 
                                                   3 (Vec4. (.x m) (.y m) (.z m) v)
                                                   m))
  ArrayVec
  (mat/set-1d [^ArrayVec m ^long row ^double v]
    (let [^doubles a (aclone ^doubles (.array m))]
      (aset a row v)
      (ArrayVec. a)))
  (mat/is-mutable? [_] true))

(extend-protocol mat/PIndexedSettingMutable
  ArrayVec
  (mat/set-1d! [^ArrayVec m ^long row ^double v]
    (aset ^doubles (.array m) row v)
    m))

(extend-protocol mat/PMatrixCloning
  ArrayVec (mat/clone [^ArrayVec m] (ArrayVec. (aclone ^doubles (.array m)))))

(extend-protocol-for-types mat/PTypeInfo [Vec2 Vec3 Vec4 ArrayVec] {:element-type (constantly Double/TYPE)})
(extend-protocol-for-types mat/PArrayMetrics [Vec2 Vec3 Vec4 ArrayVec] {:nonzero-count nonzero-count})
(extend-protocol-for-types mat/PConversion [Vec2 Vec3 Vec4 ArrayVec] {:convert-to-nested-vectors vector})
(extend-protocol-for-types mat/PZeroCount [Vec2 Vec3 Vec4 ArrayVec] {:zero-count zero-count})
(extend-protocol-for-types mat/PDoubleArrayOutput [Vec2 Vec3 Vec4]
                           {:to-double-array double-array
                            :as-double-array (constantly nil)})

(extend-protocol mat/PDoubleArrayOutput
  ArrayVec
  (mat/to-double-array [^ArrayVec m] (aclone ^doubles (.array m)))
  (mat/as-double-array [^ArrayVec m] (.array m)))

(extend-protocol-for-types mat/PValueEquality [Vec2 Vec3 Vec4 ArrayVec] {:value-equals =})
(extend-protocol-for-types mat/PMatrixEquality [Vec2 Vec3 Vec4 ArrayVec] {:matrix-equals =})
(extend-protocol-for-types mat/PMatrixEqualityEpsilon [Vec2 Vec3 Vec4 ArrayVec]
                           {:matrix-equals-epsilon (fn [m a ^double eps]
                                                     (and (= (type m) (type a))
                                                          (every? #(< (m/abs ^double %) eps) (abs (sub m a)))))})

(extend-protocol-for-types mat/PMatrixProducts [Vec2 Vec3 Vec4 ArrayVec] {:inner-product dot})
(extend-protocol-for-types mat/PMatrixMultiply [Vec2 Vec3 Vec4 ArrayVec] {:element-multiply emult})
(extend-protocol-for-types mat/PAddScaled [Vec2 Vec3 Vec4 ArrayVec] {:add-scaled #(add %1 (mult %2 %3))})
(extend-protocol-for-types mat/PMatrixDivide [Vec2 Vec3 Vec4 ArrayVec] {:element-divide div})
(extend-protocol-for-types mat/PMatrixScaling [Vec2 Vec3 Vec4 ArrayVec] {:scale mult})
(extend-protocol-for-types mat/PMatrixAdd [Vec2 Vec3 Vec4 ArrayVec] {:matrix-add add :matrix-sub sub})
(extend-protocol-for-types mat/PLerp [Vec2 Vec3 Vec4 ArrayVec] {:lerp interpolate})
(extend-protocol-for-types mat/POrder [Vec2 Vec3 Vec4 ArrayVec] {:order permute})
(extend-protocol-for-types mat/PNumerical [Vec2 Vec3 Vec4 ArrayVec] {:numerical? (constantly true)})

(extend-protocol-for-types mat/PVectorOps [Vec2 Vec3 Vec4 ArrayVec]
                           {:vector-dot dot
                            :length mag
                            :length-squared magsq
                            :normalise normalize})

(extend-protocol-for-types mat/PVectorCross [Vec2 Vec3 Vec4 ArrayVec] {:cross-product cross})
(extend-protocol-for-types mat/PVectorDistance [Vec2 Vec3 Vec4 ArrayVec] {:distance dist})
(extend-protocol-for-types mat/PVectorView [Vec2 Vec3 Vec4 ArrayVec] {:as-vector identity})
(extend-protocol-for-types mat/PVectorisable [Vec2 Vec3 Vec4 ArrayVec] {:to-vector identity})
(extend-protocol-for-types mat/PNegation [Vec2 Vec3 Vec4 ArrayVec] {:negate sub})
(extend-protocol-for-types mat/PSummable [Vec2 Vec3 Vec4 ArrayVec] {:element-sum sum})
(extend-protocol-for-types mat/PExponent [Vec2 Vec3 Vec4 ArrayVec]
                           {:element-pow (fn [m exponent] (fmap m #(m/pow % exponent)))})
(extend-protocol-for-types mat/PSquare [Vec2 Vec3 Vec4 ArrayVec] {:square sq})
(extend-protocol-for-types mat/PLogistic [Vec2 Vec3 Vec4 ArrayVec] {:logistic sigmoid})
(extend-protocol-for-types mat/PSoftplus [Vec2 Vec3 Vec4 ArrayVec] {:softplus (comp log1p exp)})
(extend-protocol-for-types mat/PReLU [Vec2 Vec3 Vec4 ArrayVec] {:relu (fn [m] (fmap m #(max 0.0 ^double %)))})
(extend-protocol-for-types mat/PSoftmax [Vec2 Vec3 Vec4 ArrayVec] {:softmax #(let [e (exp %)] (div e (sum e)))})
(extend-protocol-for-types mat/PMathsFunctions [Vec2 Vec3 Vec4 ArrayVec]
                           {:abs abs :acos acos :asin asin :atan atan
                            :cbrt cbrt :ceil ceil :cos cos :cosh cosh
                            :exp exp :floor floor :log log :log10 log10
                            :round round :signum signum :sin sin :sinh sinh
                            :sqrt sqrt :tan tan :tanh tanh :to-degrees degrees :to-radians radians})
(extend-protocol-for-types mat/PElementCount [Vec2 Vec3 Vec4 ArrayVec] {:element-count count})
(extend-protocol-for-types mat/PElementMinMax [Vec2 Vec3 Vec4 ArrayVec] {:element-min mn
                                                                         :element-max mx
                                                                         :element-clamp clamp})
(extend-protocol-for-types mat/PCompare [Vec2 Vec3 Vec4 ArrayVec]
                           {:element-compare #(signum (sub %1 %2))
                            :element-if (fn [m a b]
                                          (let [aa (if (vector? a) a (constantly a))
                                                bb (if (vector? b) b (constantly b))]
                                            (as-vec m (map-indexed (fn [idx ^double v]
                                                                     (if (pos? v) (aa idx) (bb idx))) m))))
                            :element-lt (fn [m ^double a] (fmap m #(if (< ^double % a) 1 0)))
                            :element-le (fn [m ^double a] (fmap m #(if (<= ^double % a) 1 0)))
                            :element-gt (fn [m ^double a] (fmap m #(if (> ^double % a) 1 0)))
                            :element-ge (fn [m ^double a] (fmap m #(if (>= ^double % a) 1 0)))
                            :element-ne (fn [m ^double a] (fmap m #(if (not== ^double % a) 1 0)))
                            :element-eq (fn [m ^double a] (fmap m #(if (== ^double % a) 1 0)))})

(extend-protocol-for-types mat/PFunctionalOperations [Vec2 Vec3 Vec4 ArrayVec]
                           {:element-seq identity})

;;
