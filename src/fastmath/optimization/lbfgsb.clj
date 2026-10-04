(ns fastmath.optimization.lbfgsb
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [org.generateme.lbfgsb Parameters Parameters$LINESEARCH LBFGSB IGradFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn- line-search-methods
  [linesearch]
  (case linesearch
    (:orig :more-thuente) Parameters$LINESEARCH/MORETHUENTE_ORIG
    (:lbfgsb :more-thuente-lbfgspp) Parameters$LINESEARCH/MORETHUENTE_LBFGSPP
    :levis-overton Parameters$LINESEARCH/LEWISOVERTON))

(defn parameters
  [{:keys [^int m ^double rel ^double abs ^int past ^double delta ^int max-iters ^int max-submin ^int max-linesearch
           linesearch ^double xtol ^double min-step ^double max-step ^double ftol ^double wolfe ^boolean weak-wolfe? ^boolean debug?]
    :or {debug? false m 6 abs 1.0e-8 rel 1.0e-8 past 3 delta 1.0e-10 max-iters 1000 max-submin 10 max-linesearch 20
         linesearch :more-thuente xtol 1.0e-8 min-step 1.0e-20 max-step 1.0e20 ftol 1.0e-4 wolfe 0.9 weak-wolfe? true}}]
  (when-not (m/pos? m) (throw (ex-info "m must be positive" {:m m})))
  (when-not (m/pos? abs) (throw (ex-info "abs must be positive" {:abs abs})))
  (when-not (m/pos? rel) (throw (ex-info "rel must be positive" {:rel rel}))) 
  (when-not (m/pos? max-linesearch) (throw (ex-info "max-lineasearch must be positive" {:max-linesearch max-linesearch})))
  (when-not (m/pos? min-step) (throw (ex-info "min-step must be positive" {:min-step min-step})))
  (when-not (m/not-neg? past) (throw (ex-info "past must be non-negative" {:past past})))
  (when-not (m/not-neg? delta) (throw (ex-info "delta must be non-negative" {:delta delta})))
  (when-not (m/not-neg? max-iters) (throw (ex-info "max-iters must be non-negative" {:max-iters max-iters})))
  (when-not (m/not-neg? max-submin) (throw (ex-info "max-submin must be non-negative" {:max-submin max-submin})))
  (when-not (m/>= max-step min-step) (throw (ex-info "max-step must be greater than min-step" {:min-step min-step :max-step max-step})))
  (when-not (and (m/< 0.0 ftol 0.5)
                 (m/< ftol wolfe 1.0)) (throw (ex-info "ftol and wolfe must satisfy 0<ftol<0.5 and ftol<wolfe<1.0" {:ftol ftol :wolfe wolfe})))
  (let [^Parameters p (Parameters.)]
    (set! org.generateme.lbfgsb.Debug/DEBUG debug?)
    (set! (.-m p) m)
    (set! (.-epsilon p) abs)
    (set! (.-epsilon_rel p) rel)
    (set! (.-past p) past)
    (set! (.-delta p) delta)
    (set! (.-max_iterations p) max-iters)
    (set! (.-max_submin p) max-submin)
    (set! (.-max_linesearch p) max-linesearch)
    (set! (.-linesearch p) (line-search-methods linesearch))
    (set! (.-xtol p) xtol)
    (set! (.min_step p) min-step)
    (set! (.max_step p) max-step)
    (set! (.ftol p) ftol)
    (set! (.wolfe p) wolfe)
    (set! (.weak_wolfe p) weak-wolfe?)
    p))

(defn grad-function
  (^IGradFunction [f goal] (grad-function f goal true))
  (^IGradFunction [f goal vector-arg?] (grad-function f goal vector-arg? 1.0e-6))
  (^IGradFunction [f goal vector-arg? ^double tol]
   (if (sequential? f) ;; f and grad functions provided
     (let [[f grad] f]
       (if (= goal :minimize)
         (if vector-arg?
           (reify IGradFunction
             (evaluate [_ xs] (f xs))
             (gradient [_ xs g] (let [res (grad xs)] (System/arraycopy (m/seq->double-array res) 0 ^doubles g 0 (count res)))))
           (reify IGradFunction
             (evaluate [_ xs] (apply f xs))
             (gradient [_ xs g] (let [res (grad xs)] (System/arraycopy (m/seq->double-array res) 0 ^doubles g 0 (count res))))))
         (if vector-arg?
           (reify IGradFunction
             (evaluate [_ xs] (m/- (double (f xs))))
             (gradient [_ xs g]
               (let [res (v/sub (grad xs))]
                 (System/arraycopy (double-array res) 0 ^doubles g 0 (count res)))))
           (reify IGradFunction
             (evaluate [_ xs] (m/- (double (apply f xs))))
             (gradient [_ xs g]
               (let [res (v/sub (grad xs))]
                 (System/arraycopy (double-array res) 0 ^doubles g 0 (count res))))))))
     ;; no gradient provided
     (if (= goal :minimize)
       (if vector-arg?
         (reify IGradFunction
           (evaluate [_ xs] (f xs))
           (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol)))
         (reify IGradFunction
           (evaluate [_ xs] (apply f xs))
           (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol))))
       (if vector-arg?
         (reify IGradFunction
           (evaluate [_ xs] (m/- (double (f xs))))
           (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol)))
         (reify IGradFunction
           (evaluate [_ xs] (m/- (double (apply f xs))))
           (gradient [this xs grad] (.gradient ^IGradFunction this xs grad tol))))))))


(defn lbfgsb-data
  ([f opts] (lbfgsb-data f nil opts))
  ([f gradient {:keys [goal bounds gradient-h initial vector-arg? stats?]
                :or {goal :minimize gradient-h 1.0e-6 vector-arg? true stats? false}
                :as opts}]
   (let [gf (grad-function (if gradient [f gradient] f) goal vector-arg? gradient-h)
         l (double-array (map first bounds))
         u (double-array (map second bounds))
         params (parameters opts)
         initial (m/seq->double-array (or initial (v/interpolate l u 0.5)))]
     {:gf gf :l l :u u :initial initial :params params :stats? stats? :goal goal})))

(defn lbfgsb
  ([f opts] (lbfgsb (lbfgsb-data f opts)))
  ([f gradient opts] (lbfgsb (lbfgsb-data f gradient opts)))
  ([{:keys [^IGradFunction gf ^doubles l ^doubles u ^doubles initial ^Parameters params stats? goal]}]
   (let [^LBFGSB optimizer (LBFGSB. params)
         x (.minimize optimizer gf initial l u)
         res (if (= goal :minimize) (.-fx optimizer) (m/- (.-fx optimizer)))]
     (if-not stats?
       [(vec x) res]
       {:point (vec x)
        :value res
        :iterations (.-k optimizer)
        :gradient (vec (.m_grad optimizer))}))))

(defn update-initial
  [data initial]
  (assoc data :initial (m/seq->double-array initial)))
