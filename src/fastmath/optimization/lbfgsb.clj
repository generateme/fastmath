(ns fastmath.optimization.lbfgsb
  (:require [fastmath.core :as m])
  (:import [org.generateme.lbfgsb Parameters Parameters$LINESEARCH LBFGSB IGradFunction]))

(m/use-primitive-operators)
(set! *unchecked-math* :warn-on-boxed)

(def ^:private line-search-methods
  {:orig Parameters$LINESEARCH/MORETHUENTE_ORIG
   :lbfgsb Parameters$LINESEARCH/MORETHUENTE_LBFGSPP})

(defn- parameters
  [{:keys [m rel abs past delta max-iters max-submin max-linesearch
           linesearch xtol min-step max-step ftol wolfe weak-wolfe? debug?]}]
  (let [^Parameters p (Parameters.)]
    (set! org.generateme.lbfgsb.Debug/DEBUG (boolean debug?))
    (when m (set! (.-m p) (int m)))
    (when abs (set! (.-epsilon p) (double abs)))
    (when rel (set! (.-epsilon_rel p) (double rel)))
    (when past (set! (.-past p) (int past)))
    (when delta (set! (.-delta p) (double delta)))
    (when max-iters (set! (.-max_iterations p) (int max-iters)))
    (when max-submin (set! (.-max_submin p) (int max-submin)))
    (when max-linesearch (set! (.-max_linesearch p) (int max-linesearch)))
    (when linesearch (set! (.-linesearch p)
                           (get line-search-methods linesearch Parameters$LINESEARCH/MORETHUENTE_ORIG)))
    (when xtol (set! (.-xtol p) (double xtol)))
    (when min-step (set! (.min_step p) (double min-step)))
    (when max-step (set! (.max_step p) (double max-step)))
    (when ftol (set! (.ftol p) (double ftol)))
    (when wolfe (set! (.wolfe p) (double wolfe)))
    (when weak-wolfe? (set! (.weak_wolfe p) (boolean weak-wolfe?)))
    (assert (m/pos? (.m p)) ":m must be positive")
    (assert (m/pos? (.epsilon p)) ":tolerance must be positive")
    (assert (m/pos? (.epsilon_rel p)) ":tolerance-rel must be positive")
    (assert (m/not-neg? (.past p)) ":past must be non-negative")
    (assert (m/not-neg? (.delta p)) ":delta must be non-negative")
    (assert (m/not-neg? (.max_iterations p)) ":max-iterations must be non-negative")
    (assert (m/not-neg? (.max_submin p)) ":max-submin must be non-negative")
    (assert (m/pos? (.max_linesearch p)) ":max-linesearch must be positive")
    (assert (m/pos? (.min_step p)) ":min-step must be positive")
    (assert (>= (.max_step p) (.min_step p)) ":max-step must be greater than :min-step ")
    (assert (< 0.0 (.ftol p) 0.5) ":ftol must satisfy 0<ftol<0.5")
    (assert (< (.ftol p) (.wolfe p) 1.0) ":wolfe must satisfy ftol<wolfe<1")
    p))

(defn grad-function
  ([f goal] (grad-function f goal nil))
  ([f goal tol]
   (if (sequential? f) ;; f and grad functions provided
     (let [[f grad] f]
       (if (= goal :minimize)
         (reify IGradFunction
           (evaluate [_ xs] (apply f xs))
           (gradient [_ xs g]
             (let [res (apply grad xs)]
               (System/arraycopy (m/seq->double-array res) 0 ^doubles g 0 (count res)))))
         (reify IGradFunction
           (evaluate [_ xs] (- ^double (apply f xs)))
           (gradient [_ xs g]
             (let [res (map (fn [^double v] (- v)) (apply grad xs))]
               (System/arraycopy (double-array res) 0 ^doubles g 0 (count res)))))))
     (if tol ;; tolerance for autograd provided
       (let [tol (double tol)]
         (if (= goal :minimize)
           (reify IGradFunction
             (evaluate [_ xs] (apply f xs))
             (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol)))
           (reify IGradFunction
             (evaluate [_ xs] (- ^double (apply f xs)))
             (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol)))))
       (if (= goal :minimize)
         (reify IGradFunction
           (evaluate [_ xs] (apply f xs)))
         (reify IGradFunction
           (evaluate [_ xs] (- ^double (apply f xs)))))))))

(defn ->lbfgsb
  "Create lbfgsb object"
  (^LBFGSB [] (->lbfgsb nil))
  (^LBFGSB [params] (LBFGSB. (parameters params))))

(m/unuse-primitive-operators)
