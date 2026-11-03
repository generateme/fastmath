(ns ^:no-doc fastmath.calculus.autodiff.dual
  "Internal: the dual number type of [[fastmath.calculus.autodiff]], its tags, layer access and partials containers (no arithmetic)."
  (:require [fastmath.core :as m])
  (:import [clojure.lang IHashEq LazilyPersistentVector Util]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Tags

(defonce ^:private tag-counter (atom 0))

(defn fresh-tag
  "Returns a new tag, larger than all tags returned before."
  ^long []
  (long (swap! tag-counter inc)))

;; ## The dual type

;; Invariant: every tag found in `value` and in `partials` is smaller than `tag`.
;; `partials` is an `Object[]` which is never modified after the dual is created (the hash is cached and duals are
;; values); it is not exposed by the public API, results are converted to vectors.
;; `hash-cache` holds the hash once computed (0 = not yet), as Clojure collections do; it is always created by
;; `new-dual`.
(deftype Dual [^long tag value ^objects partials ^:unsynchronized-mutable ^int hash-cache]
  Object
  (toString [_] (str "#dual[" tag " " (pr-str value) " " (pr-str (vec partials)) "]"))
  ;; Structural equality (`clojure.core/=`, map keys, sets): the same tag and partials, the numeric leaves compared
  ;; as by `clojure.core/=` (`0.0` equals `-0.0`, `NaN` equals only the same object), consistent with `hasheq`.
  ;; Numeric comparison is `fastmath.calculus.autodiff/==`.
  (equals [this o]
    (or (identical? this o)
        (and (instance? Dual o)
             (m/== tag (.tag ^Dual o))
             (Util/equiv value (.value ^Dual o))
             (let [^objects q (.partials ^Dual o)
                   n (alength partials)]
               (and (m/== n (alength q))
                    (loop [i 0]
                      (cond
                        (m/== i n) true
                        (Util/equiv (aget partials i) (aget q i)) (recur (m/long-inc i))
                        :else false)))))))
  (hashCode [this] (.hasheq ^IHashEq this))
  IHashEq
  (hasheq [_]
    (when (m/zero? hash-cache)
      (set! hash-cache (Util/hashCombine (Util/hashCombine (Long/hashCode tag) (Util/hasheq value))
                                         (.hasheq ^IHashEq (LazilyPersistentVector/createOwning partials)))))
    hash-cache))

(defmacro new-dual
  "A new `Dual` with an empty hash cache."
  [tag value partials]
  `(Dual. ~tag ~value ~partials 0))

(defmethod print-method Dual [^Dual d ^java.io.Writer w]
  (.write w (.toString d)))

(defn leaf
  "A dual unchanged, a number coerced to a double."
  [x]
  (if (instance? Dual x) x (.doubleValue ^Number x)))

(defn scalar
  "Value of a number or a dual with all perturbations removed, as a double."
  ^double [x]
  (if (instance? Dual x) (recur (.value ^Dual x)) (.doubleValue ^Number x)))

;; ## Tag layers

(defn tag-of
  "Tag of a dual, `-1` for a number."
  ^long [x]
  (if (instance? Dual x) (.tag ^Dual x) -1))

(defn top-tag
  "The newer of the tags of `x` and `y`."
  ^long [x y]
  (m/long-max (tag-of x) (tag-of y)))

(defn value-at
  "Value of `x` in the layer of tag `t`: the value of a dual tagged `t`, otherwise `x` itself (a constant)."
  [x ^long t]
  (if (and (instance? Dual x) (m/== t (.tag ^Dual x))) (.value ^Dual x) x))

(defn partials-at
  "Partials of `x` in the layer of tag `t`, `nil` (zero) for a constant."
  [x ^long t]
  (when (and (instance? Dual x) (m/== t (.tag ^Dual x))) (.partials ^Dual x)))

;; ## Partials (an `Object[]`, `nil` is a zero vector)

;; Partials are plain arrays, filled by index loops: no persistent vector is allocated per operation, and no
;; `mapv`, which goes through a lazy sequence (two collections) or a transient vector (one collection). An array is
;; never modified after it is stored in a dual. The hash wraps it as a vector with
;; `LazilyPersistentVector/createOwning` (internal to Clojure but stable for a long time), without copying for up
;; to 32 partials.

(defn zero-partials
  "`k` zero partials."
  ^objects [^long k]
  (let [a (object-array k)]
    (java.util.Arrays/fill a (Double/valueOf 0.0))
    a))

(defn basis-partials
  "`k` partials, zero except `1.0` at index `i`."
  ^objects [^long k ^long i]
  (let [a (object-array k)]
    (java.util.Arrays/fill a (Double/valueOf 0.0))
    (aset a i (Double/valueOf 1.0))
    a))

(defn p-zeros
  "Zero partials of the length of `p`."
  ^objects [^objects p]
  (zero-partials (alength p)))

(defn p-map
  "`f` applied to each partial, in index order."
  ^objects [f ^objects p]
  (let [n (alength p)
        a (object-array n)]
    (loop [i 0]
      (when (m/< i n)
        (aset a i (f (aget p i)))
        (recur (m/long-inc i))))
    a))

(defn p-map2
  "`f` applied to the partials of the same index, in index order; both arrays have the same length."
  ^objects [f ^objects px ^objects py]
  (let [n (alength px)
        a (object-array n)]
    (loop [i 0]
      (when (m/< i n)
        (aset a i (f (aget px i) (aget py i)))
        (recur (m/long-inc i))))
    a))

(defn check-same-length!
  "Duals with the same tag must have the same number of partials (a precondition of the caller)."
  [^objects px ^objects py]
  (when-not (m/== (alength px) (alength py))
    (throw (ex-info "Duals with the same tag have different numbers of partials."
                    {:count-x (alength px) :count-y (alength py)}))))

(declare zero-element?)

(defn p-zero?
  "Are all partials zero (`nil` is zero)? Stops at the first nonzero element."
  [p]
  (or (nil? p)
      (let [^objects v p
            n (alength v)]
        (loop [i 0]
          (cond
            (m/== i n) true
            (zero-element? (aget v i)) (recur (m/long-inc i))
            :else false)))))

(defn zero-element?
  "Is `x` zero, including all of its partials? The value is checked first, so a nonzero scalar value is found by following the values only."
  [x]
  (if (instance? Dual x)
    (and (zero-element? (.value ^Dual x)) (p-zero? (.partials ^Dual x)))
    (m/zero? (.doubleValue ^Number x))))
