(ns fastmath.dev.regression-ref-gen
  "Generator of R reference values for `fastmath.ml.regression-test`.

  Produces `test/resources/regression/<name>_reference.edn` files, each holding
  the R-computed oracle `:model`/`:analysis`/`:extra` fields for one deftest.
  Tier 1 (hardcoded-input deftests) only, so far -- see the topic note
  \"Fastmath Regression Test R and Tablecloth Removal\".

  Each file's variant block is preceded by a plain EDN comment recording the
  R formula/packages used to produce it (per spec: provenance as comments, not
  a structured key -- no live R session is needed to run the test suite, only
  the checked-in .edn output is read at test time.

  Run from a dev REPL:
    (require '[fastmath.dev.regression-ref-gen :as g] :reload)
    (g/-main)"
  (:require [clojisr.v1.r :as rr]
            [clojisr.v1.impl.protocols]
            [tablecloth.api :as tc]
            [fastmath.ml.regression :as reg]
            [fastmath.ml.regression.contrast :as contrast]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.java.io :as io]
            [clojure.pprint :as pp]))

;; Same workaround as fastmath.ml.regression-test ("to move to clojisr"): R
;; returns a bare NA/NULL for some quasi-family statistics (e.g. logLik on a
;; quasibinomial glm), which `rr/r->clj` can't convert without this.
(extend nil
  clojisr.v1.impl.protocols/Clojable
  {:->clj (fn [_] nil)})

(rr/require-r '[stats] '[car] '[moments] '[base] '[MASS] '[statmod] '[GLMsData] '[utils])

(defn- dataset-columns
  "Export named columns from an R data.frame `ds` (already `rr/r->clj`'d, a
  tablecloth dataset -- confirmed via `rr/r->clj`) as plain Clojure vectors,
  for the §1.3 :inputs.columns section. No forced numeric coercion: a
  categorical column comes back as a vector of keywords (R factor levels),
  a numeric column as whatever number type clojisr already produced."
  [ds col-names]
  (into {} (for [c col-names] [c (vec (seq (tc/column ds c)))])))

(defn- one-hot
  "§2.3: pure-fastmath replacement for the old tablecloth-based `one-hot`.
  `col` is a raw categorical column (a seq of keywords); `values` are the
  non-reference levels to build 0/1 columns for, matching the original
  tablecloth `one-hot`'s calling convention exactly. The reference level is
  whichever level of `col` is not in `values` (assumed exactly one -- true
  for every current call site, confirmed against each dataset's actual
  levels before use, not assumed). Delegates the actual encoding to
  `contrast/dummy`. Returns {value -> column-vector}, values in `values` order."
  [col values]
  (let [levels (distinct col)
        reference (first (remove (set values) levels))
        coding (contrast/dummy (cons reference values))
        columns (apply map vector (map (:mapping coding) col))]
    (zipmap values columns)))

(defn- one-hot-mapping
  "Same `col`/`values`/reference-level convention as `one-hot`, but returns
  the per-row encoding function directly (`raw-level -> dummy row`, in
  `values` order) instead of applying it to `col`'s own rows -- for encoding
  *new* raw values (predict newdata) with the exact training-time reference
  level, which `one-hot` itself can't do (it only maps over its own `col`)."
  [col values]
  (let [levels (distinct col)
        reference (first (remove (set values) levels))]
    (:mapping (contrast/dummy (cons reference values)))))

(defn- mult-columns
  "§2.3 extension: pure-fastmath replacement for the old tablecloth-based
  `mult-columns`. `col1`/`col2` are raw categorical columns; non-reference
  levels default to first-appearance order minus the first level (matching
  the old `mult-columns`' default `sort1?`/`sort2?` = false). Returns
  {:cs1 :cs2 :xss}, `:xss` in the same column order tablecloth produced:
  one-hot(col1) columns, then one-hot(col2) columns, then every (c2,c1)
  interaction column, c2 outer / c1 inner (matching `(for [c2 cs2 c1 cs1] ...)`)."
  [col1 col2]
  (let [cs1 (vec (rest (distinct col1)))
        cs2 (vec (rest (distinct col2)))
        oh1 (one-hot col1 cs1)
        oh2 (one-hot col2 cs2)
        onehot1-cols (map oh1 cs1)
        onehot2-cols (map oh2 cs2)
        interactions (for [c2 cs2 c1 cs1] (mapv * (oh1 c1) (oh2 c2)))]
    {:cs1 cs1 :cs2 cs2 :xss (apply map vector (concat onehot1-cols onehot2-cols interactions))}))

;; §2.4 (superseded): the original `add-poly` (a pure-fastmath, tablecloth-free
;; port wrapping `v/orthonormal-polynomials`) has been replaced by
;; `reg/polynomial-transformer` at both `danishlc-data`'s and `deposit-data`'s
;; poly call sites -- see the topic note "Fastmath Regression Test R and
;; Tablecloth Removal" (Continuation, `polynomial-transformer` bugfix). Unlike
;; `add-poly`, `polynomial-transformer` returns a genuine `:transformer`
;; closure that can be evaluated at *new* x values (it fixes its centering/
;; recurrence coefficients once, from the training `xs`), enabling a
;; newdata-predict check for both poly-encoded deftests, previously out of
;; scope for exactly this reason.

;; ---- shared conversions ----

(defn- rvec
  "Convert an R value to a flat Clojure vector of doubles."
  [rres]
  (mapv double (rr/r->clj rres)))

(defn- coef-groups
  "Split R's flattened `summary$coefficients` matrix (column-major: estimate,
  stderr, t/z-value, p-value) into the four named vectors."
  [flat]
  (let [n (/ (count flat) 4)
        [e s tv p] (partition n flat)]
    {:estimate (vec e) :stderr (vec s) :t-value (vec tv) :p-value (vec p)}))

;; ---- LM oracle extraction ----

(defn- lm-oracle
  "Build the {:model :analysis} reference map for one fitted R `lm` object `rlm`.
  `intercept?` MUST match the fastmath `:intercept?` option used for this variant:
  R's `:coefficients` omits the intercept row entirely for a no-intercept formula,
  so `first`/`rest` on it is only valid when an intercept is present."
  [rlm intercept?]
  (let [rlmdata (rr/r->clj rlm)
        summary (rr/r->clj `(summary ~rlm :correlation true))
        influence (rr/r->clj `(influence ~rlm))
        lobs (m/log (count (:residuals rlmdata)))
        rcoefs (mapv double (:coefficients rlmdata))]
    {:model
     {:intercept (if intercept? (first rcoefs) 0.0)
      :beta (vec (if intercept? (rest rcoefs) rcoefs))
      :fitted (rvec (:fitted.values rlmdata))
      :sigma (double (first (:sigma summary)))
      :sigma2 (m/sq (double (first (:sigma summary))))
      :r-squared (double (first (:r.squared summary)))
      :adjusted-r-squared (double (first (:adj.r.squared summary)))
      :f-statistic (double (first (:fstatistic summary)))
      :xtxinv (rvec (:cov.unscaled summary))
      :residuals {:raw (rvec (:residuals rlmdata))
                  :weighted (rvec (:residuals summary))}
      :df {:residual (long (first (:df.residual rlmdata)))
           :model (long (second (:fstatistic summary)))}
      :ll {:log-likelihood (double (first (rr/r->clj (stats/logLik rlm))))
           :aic-rss (double (second (rr/r->clj (stats/extractAIC rlm))))
           :aic (double (first (rr/r->clj (stats/AIC rlm))))
           :bic-rss (double (second (rr/r->clj (stats/extractAIC rlm :k lobs))))
           :bic (double (first (rr/r->clj (stats/AIC rlm :k lobs))))}
      :coefficients (coef-groups (rvec (:coefficients summary)))}
     :analysis
     {:correlation (rvec (:correlation summary))
      :normality {:skewness (double (first (rr/r->clj (moments/skewness (stats/weighted-residuals rlm)))))
                  :kurtosis (double (first (rr/r->clj (moments/kurtosis (stats/weighted-residuals rlm)))))
                  :durbin-watson (double (first (:dw (rr/r->clj (car/durbinWatsonTest rlm)))))}
      :residuals {:standardized (rvec (stats/rstandard rlm))
                  :studentized (rvec (stats/rstudent rlm))}
      :leverage {:hat (rvec (:hat influence))
                 :sigmas (rvec (:sigma influence))
                 :coefficients (rvec (:coefficients influence))}
      :influence {:cooks-distance (rvec (stats/cooks-distance rlm))
                  :dffits (rvec (stats/dffits rlm))
                  :covratio (rvec (stats/covratio rlm))
                  :dfbetas (vec (flatten (rest (vals (rr/r->clj (stats/dfbetas rlm))))))}}}))

;; ---- GLM oracle extraction ----

(defn- glm-oracle
  "Build the {:model :analysis} reference map for one fitted R `rglm` object.
  `intercept?` MUST match the fastmath `:intercept?` option (same reason as
  `lm-oracle`). `stat-key` is `:t-value` for an estimated-dispersion family
  (gaussian/gamma/inverse-gaussian/quasi-*) or `:z-value` for a fixed-dispersion
  family (binomial/poisson) -- MUST match fastmath's own `:estimated-dispersion?`
  branch, since `check-coefficients` reads whichever key is present."
  [rglm intercept? stat-key]
  (let [rglmdata (rr/r->clj rglm)
        summary (rr/r->clj `(summary ~rglm :correlation true))
        influence (rr/r->clj `(influence ~rglm))
        rcoefs (mapv double (:coefficients rglmdata))
        flat-summary-coefs (rvec (:coefficients summary))
        n (/ (count flat-summary-coefs) 4)
        [e s tv p] (partition n flat-summary-coefs)]
    {:model
     {:intercept (if intercept? (first rcoefs) 0.0)
      :beta (vec (if intercept? (rest rcoefs) rcoefs))
      :fitted (rvec (:fitted.values rglmdata))
      :dispersion (double (first (:dispersion summary)))
      :xtxinv (rvec (:cov.unscaled summary))
      :weights {:weights (rvec (:weights rglmdata))
                :initial (rvec (:prior.weights rglmdata))}
      :residuals {:working (rvec (:residuals rglmdata))
                  :deviance (rvec (:deviance.resid summary))
                  :pearson (rvec (rr/r->clj `(residuals ~rglm :type "pearson")))
                  :raw (rvec (rr/r->clj `(residuals ~rglm :type "response")))}
      :deviance {:residual (double (first (:deviance rglmdata)))
                 :null (double (first (:null.deviance rglmdata)))}
      :df {:residual (long (first (:df.residual rglmdata)))
           :null (long (first (:df.null rglmdata)))}
      ;; R's logLik/AIC/BIC are undefined for quasi-families (no true likelihood) --
      ;; return NaN rather than throwing; harmless since such variants set :no-ll?
      ;; and `check-ll` is never called against these values.
      :ll {:log-likelihood (double (or (first (rr/r->clj (stats/logLik rglm))) ##NaN))
           :aic (double (or (first (rr/r->clj (stats/AIC rglm))) ##NaN))
           :bic (double (or (first (rr/r->clj (stats/BIC rglm))) ##NaN))}
      :coefficients (merge {:estimate (vec e) :stderr (vec s) :p-value (vec p)}
                            {stat-key (vec tv)})}
     :analysis
     {:correlation (rvec (:correlation summary))
      :residuals {:studentized (rvec (stats/rstudent rglm))
                  :standardized {:deviance (rvec (rr/r->clj (stats/rstandard rglm :type "deviance")))
                                 :pearson (rvec (rr/r->clj (stats/rstandard rglm :type "pearson")))}}
      :leverage {:hat (rvec (:hat influence))
                 :sigmas (rvec (:sigma influence))
                 :coefficients (rvec (:coefficients influence))}
      :influence {:cooks-distance (rvec (stats/cooks-distance rglm))
                  :dffits (rvec (stats/dffits rglm))
                  :covratio (rvec (stats/covratio rglm))
                  :dfbetas (vec (flatten (rest (vals (rr/r->clj (stats/dfbetas rglm))))))}}}))

(defn- qresid-oracle
  "R's `statmod::qresid` output for `rglm`, for §3.3's `:qres?`-gated deftests."
  [rglm]
  (rvec (rr/r->clj (statmod/qresid rglm))))

(defn- dose-oracle
  "Attempt `MASS::dose.p` on `rglm`. Returns nil (not an error) if the family/link
  combination doesn't support it -- callers decide whether nil is expected."
  [rglm]
  (try
    (let [rdose (MASS/dose-p rglm)
          dattr (rr/r->clj (base/attributes rdose))]
      {:dose (double (first (rr/r->clj rdose)))
       :stderr (double (first (:SE dattr)))})
    (catch Exception e
      {:error (str e)})))

;; ---- New-data prediction oracle (fit/stderr/confidence-interval on genuinely
;; new points, not the training data) -- companion to `fastmath.ml.regression`'s
;; `predict`/`IFn` with `stderr?`=true. Requires `rlm`/`rglm` to have been fit
;; via a *named* data.frame formula (not literal-embedded vectors), so R's own
;; `predict(object, newdata=...)` can resolve columns by name; every Tier 2/3
;; dataset already fits this way (`:data <dataset>`), Tier 1 needs a dedicated
;; named-data.frame refit (see `predict-lm-refit`/`predict-glm-refit` at each
;; Tier 1 build site).
;;
;; `predict.lm`'s own `interval=\"confidence\"` is used directly for LM. GLM has
;; no analogous built-in interval argument, so `predict_glm_ci` reconstructs it
;; exactly the way `fastmath.ml.regression`'s own `predict` does: se.fit on the
;; *link* scale, `qt`/`qnorm` matching fastmath's own `:estimated-dispersion?`
;; branch (`stat-key-for`'s fixed- vs. estimated-dispersion split), then
;; `family$linkinv` applied to fit and to fit +/- q*se for the response-scale CI.
(def ^:private predict-r-helpers-loaded?
  (do
    (rr/r "
predict_lm_ci <- function(model, newdata, alpha=0.05) {
  p <- predict(model, newdata=newdata, se.fit=TRUE, interval='confidence', level=1-alpha)
  list(fit=as.vector(p$fit[,'fit']), lwr=as.vector(p$fit[,'lwr']), upr=as.vector(p$fit[,'upr']),
       se=as.vector(p$se.fit))
}
predict_glm_ci <- function(model, newdata, estimated_dispersion, alpha=0.05) {
  p <- predict(model, newdata=newdata, type='link', se.fit=TRUE)
  fam <- family(model)
  q <- if (estimated_dispersion) qt(1-alpha/2, df.residual(model)) else qnorm(1-alpha/2)
  list(fit=as.vector(fam$linkinv(p$fit)), se=as.vector(p$se.fit),
       lwr=as.vector(fam$linkinv(p$fit - q*p$se.fit)), upr=as.vector(fam$linkinv(p$fit + q*p$se.fit)))
}
")
    true))

(defn- predict-lm-oracle
  "`newdata-form` is an R data.frame expression (clojisr syntax-quoted list),
  columns named to match `rlm`'s own formula. Returns {:fit :stderr
  :confidence-interval}, CI as a seq of [lo hi] pairs, parallel to `newdata-form`'s
  rows."
  ([rlm newdata-form] (predict-lm-oracle rlm newdata-form 0.05))
  ([rlm newdata-form alpha]
   (let [{:keys [fit lwr upr se]} (rr/r->clj (rr/r `(predict_lm_ci ~rlm ~newdata-form ~(double alpha))))]
     {:fit (rvec fit) :stderr (rvec se) :confidence-interval (mapv vector (rvec lwr) (rvec upr))})))

(defn- predict-glm-oracle
  "As `predict-lm-oracle`, `estimated-dispersion?` MUST mirror fastmath's own
  `:estimated-dispersion?` for this family (see `stat-key-for`)."
  ([rglm newdata-form estimated-dispersion?] (predict-glm-oracle rglm newdata-form estimated-dispersion? 0.05))
  ([rglm newdata-form estimated-dispersion? alpha]
   (let [{:keys [fit lwr upr se]} (rr/r->clj (rr/r `(predict_glm_ci ~rglm ~newdata-form
                                                                     ~(boolean estimated-dispersion?)
                                                                     ~(double alpha))))]
     {:fit (rvec fit) :stderr (rvec se) :confidence-interval (mapv vector (rvec lwr) (rvec upr))})))

;; ---- Tier 1 inputs (MUST match fastmath.ml.regression-test literals exactly) ----

(def ys [1 2.1 3.2 4 4.9 5.5 10 7 8.1 8.7])
(def xss-1 '(-0.5 -1 -2 -3 -4 -5 -6 -7 -8 -9))
(def xss-2 '(1 2 0 2 1 0 1 2 0 2))
(def xss-3 (v/emult xss-1 xss-2))
(def xss (map vector xss-1 xss-2 xss-3))
(def weights [2 1 2 1 2 1 2 1 2 1])
(def offset (range 10))
(def nys (v/normalize ys))

;; ---- Tier 1, basic-* (LM) variants: [name r-form r-formula-str fastmath-options] ----
;; `r-formula-str` is a hand-written, human-readable R-syntax rendering of `r-form`,
;; used only for the §1.5 provenance comment -- `r-form` itself is a syntax-quoted
;; Clojure data structure whose `pr-str` (fully namespace-qualified) is not readable.

(def ^:private basic-lm-variants
  [["basic-weights-intercept"
    `(lm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights)
    "lm(y ~ x1*x2, weights=w)"
    {:weights weights}]
   ["basic-weights-intercept-transformer"
    `(lm (formula ~ys (+ (exp ~xss-1) (* ~xss-1 ~xss-2))) :weights ~weights)
    "lm(y ~ exp(x1) + x1*x2, weights=w)"
    {:weights weights :transformer (fn [[a b c]] [(m/exp a) a b c])}]
   ["basic-weights-intercept-offset"
    `(lm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights :offset ~offset)
    "lm(y ~ x1*x2, weights=w, offset=o)"
    {:weights weights :offset offset}]
   ["basic-weights-no-intercept"
    `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights)
    "lm(y ~ 0 + x1*x2, weights=w)"
    {:weights weights :intercept? false}]
   ["basic-weights-no-intercept-offset"
    `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :offset ~offset)
    "lm(y ~ 0 + x1*x2, weights=w, offset=o)"
    {:weights weights :intercept? false :offset offset}]
   ["basic-no-weights-intercept"
    `(lm (formula ~ys (* ~xss-1 ~xss-2)))
    "lm(y ~ x1*x2)"
    nil]
   ["basic-no-weights-intercept-offset"
    `(lm (formula ~ys (* ~xss-1 ~xss-2)) :offset ~offset)
    "lm(y ~ x1*x2, offset=o)"
    {:offset offset}]
   ["basic-no-weights-no-intercept"
    `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))))
    "lm(y ~ 0 + x1*x2)"
    {:intercept? false}]
   ["basic-no-weights-no-intercept-offset"
    `(lm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :offset ~offset)
    "lm(y ~ 0 + x1*x2, offset=o)"
    {:intercept? false :offset offset}]])

;; ---- Tier 1, dummy-data (GLM binomial) variants: [key r-form r-formula-str options] ----

(def ^:private dummy-glm-variants
  [[:call-0 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family binomial)
    "glm(y ~ x1*x2, weights=w, family=binomial)" {:weights weights :family :binomial}]
   [:call-1 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link probit))
    "glm(y ~ x1*x2, weights=w, family=binomial(link=probit))"
    {:weights weights :family :binomial :link :probit}]
   [:call-2 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link cloglog))
    "glm(y ~ x1*x2, weights=w, family=binomial(link=cloglog))"
    {:weights weights :family :binomial :link :cloglog}]
   [:call-3 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :weights ~weights :family (binomial :link log))
    "glm(y ~ x1*x2, weights=w, family=binomial(link=log))"
    {:weights weights :family :binomial :link :log}]
   [:call-4 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family binomial)
    "glm(y ~ x1*x2, family=binomial)" {:family :binomial}]
   [:call-5 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link probit))
    "glm(y ~ x1*x2, family=binomial(link=probit))" {:family :binomial :link :probit}]
   [:call-6 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link cloglog))
    "glm(y ~ x1*x2, family=binomial(link=cloglog))" {:family :binomial :link :cloglog}]
   [:call-7 `(glm (formula ~nys (* ~xss-1 ~xss-2)) :family (binomial :link log))
    "glm(y ~ x1*x2, family=binomial(link=log))" {:family :binomial :link :log}]
   [:call-8 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family binomial)
    "glm(y ~ 0 + x1*x2, weights=w, family=binomial)"
    {:weights weights :family :binomial :intercept? false}]
   [:call-9 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family (binomial :link probit))
    "glm(y ~ 0 + x1*x2, weights=w, family=binomial(link=probit))"
    {:weights weights :family :binomial :link :probit :intercept? false}]
   [:call-10 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :family (binomial :link cloglog))
    "glm(y ~ 0 + x1*x2, weights=w, family=binomial(link=cloglog))"
    {:weights weights :family :binomial :link :cloglog :intercept? false}]
   [:call-11 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family binomial)
    "glm(y ~ 0 + x1*x2, family=binomial)" {:family :binomial :intercept? false}]
   [:call-12 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family (binomial :link probit))
    "glm(y ~ 0 + x1*x2, family=binomial(link=probit))" {:family :binomial :link :probit :intercept? false}]
   [:call-13 `(glm (formula ~nys (+ 0 (* ~xss-1 ~xss-2))) :family (binomial :link cloglog))
    "glm(y ~ 0 + x1*x2, family=binomial(link=cloglog))" {:family :binomial :link :cloglog :intercept? false}]
   [:call-14 `(glm (formula ~nys (+ (exp ~xss-1) (* ~xss-1 ~xss-2))) :family binomial)
    "glm(y ~ exp(x1) + x1*x2, family=binomial)"
    {:family :binomial :transformer (fn [[a b c]] [(m/exp a) a b c])}]])

;; ---- Tier 1, gaussian-glm variants: [key r-form r-formula-str options] ----

(def ^:private gaussian-glm-variants
  [[:call-0 `(glm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights)
    "glm(y ~ x1*x2, weights=w)" {:weights weights :qres? true}]
   [:call-1 `(glm (formula ~ys (* ~xss-1 ~xss-2)) :weights ~weights :offset ~offset)
    "glm(y ~ x1*x2, weights=w, offset=o)" {:weights weights :offset offset :qres? true}]
   [:call-2 `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights)
    "glm(y ~ 0 + x1*x2, weights=w)" {:weights weights :intercept? false :qres? true}]
   [:call-3 `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :weights ~weights :offset ~offset)
    "glm(y ~ 0 + x1*x2, weights=w, offset=o)"
    {:weights weights :intercept? false :offset offset :qres? true}]
   [:call-4 `(glm (formula ~ys (* ~xss-1 ~xss-2)))
    "glm(y ~ x1*x2)" {:qres? true}]
   [:call-5 `(glm (formula ~ys (* ~xss-1 ~xss-2)) :offset ~offset)
    "glm(y ~ x1*x2, offset=o)" {:offset offset :qres? true}]
   [:call-6 `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))))
    "glm(y ~ 0 + x1*x2)" {:intercept? false :qres? true}]
   [:call-7 `(glm (formula ~ys (+ 0 (* ~xss-1 ~xss-2))) :offset ~offset)
    "glm(y ~ 0 + x1*x2, offset=o)" {:intercept? false :offset offset :qres? true}]])

;; ---- file writers (§1.5: plain EDN comment provenance, no :meta key) ----

(defn- spit-single-variant!
  [dir name {:keys [r-formula r-packages data]}]
  (let [f (io/file dir (str name "_reference.edn"))]
    (io/make-parents f)
    (spit f
          (str ";; Generated by fastmath.dev.regression-ref-gen. DO NOT hand-edit.\n"
               ";; :r-formula " (pr-str r-formula) "\n"
               ";; :r-packages " (pr-str r-packages) "\n"
               (binding [*print-length* nil *print-level* nil]
                 (with-out-str (pp/pprint data)))))
    (str f)))

(defn- spit-multi-variant!
  "Write a multi-call file: an optional shared :inputs block (§1.3 -- same raw
  columns across all variants, e.g. one R dataset) and an optional shared
  file-level :extra block (§1.4 -- for values not tied to any single
  lm-tests/glm-tests call, e.g. a standalone R-vs-fastmath comparison),
  followed by several sub-keyed variant blocks, each with its own §1.5
  provenance comment (§3.2.1: one lm-tests/glm-tests call per sub-key)."
  [dir name {:keys [inputs-comment inputs extra-comment extra variants]}]
  (let [f (io/file dir (str name "_reference.edn"))]
    (io/make-parents f)
    (spit f
          (str ";; Generated by fastmath.dev.regression-ref-gen. DO NOT hand-edit.\n"
               "{"
               (when inputs
                 (str "\n ;; " inputs-comment "\n"
                      " :inputs\n"
                      (binding [*print-length* nil *print-level* nil]
                        (with-out-str (pp/pprint inputs)))))
               (when extra
                 (str "\n ;; " extra-comment "\n"
                      " :extra\n"
                      (binding [*print-length* nil *print-level* nil]
                        (with-out-str (pp/pprint extra)))))
               (apply str
                      (for [{:keys [key r-formula r-packages data]} variants]
                        (str "\n ;; :r-formula " (pr-str r-formula) "\n"
                             " ;; :r-packages " (pr-str r-packages) "\n"
                             " " (pr-str key) "\n"
                             (binding [*print-length* nil *print-level* nil]
                               (with-out-str (pp/pprint data))))))
               "}\n"))
    (str f)))

;; ---- cross-validation (§4.5) ----

(defn- max-abs-diff
  ^double [a b]
  (->> (map (fn [^double x ^double y] (m/abs (m/- x y))) a b)
       (reduce max 0.0)))

(defn- cross-validate-lm
  "Fit fastmath's `lm` on `ys`/`xss`/`options`, diff its :model fields against
  the just-generated oracle :model map. Generalized across tiers -- `ys`/`xss`
  are explicit params, not closed-over Tier 1 globals."
  [ys xss options oracle-model]
  (let [model (reg/lm ys xss options)]
    {:intercept (m/abs (m/- (double (:intercept model)) (double (:intercept oracle-model))))
     :beta (max-abs-diff (:beta model) (:beta oracle-model))
     :fitted (max-abs-diff (:fitted model) (:fitted oracle-model))
     :sigma (m/abs (m/- (double (:sigma model)) (double (:sigma oracle-model))))
     :r-squared (m/abs (m/- (double (:r-squared model)) (double (:r-squared oracle-model))))}))

;; ---- GLM cross-validation (§4.5) ----

(defn- with-glm-control
  "Append the `:epsilon 1e-16 :maxit 100` R-side control args every glm-tests
  call uses, matching `sut/glm`'s own `{:epsilon 1.0e-16 :max-iters 100}` default."
  [r-form]
  (concat r-form (list :epsilon 1.0e-16 :maxit 100)))

(defn- fit-glm
  [ys xss options]
  (reg/glm ys xss (merge {:epsilon 1.0e-16 :max-iters 100} options)))

(defn- cross-validate-glm
  "Diff a fitted fastmath `model`'s :model fields against the just-generated
  oracle :model map."
  [model oracle-model]
  {:intercept (m/abs (m/- (double (:intercept model)) (double (:intercept oracle-model))))
   :beta (max-abs-diff (:beta model) (:beta oracle-model))
   :fitted (max-abs-diff (:fitted model) (:fitted oracle-model))
   :dispersion (m/abs (m/- (double (:dispersion model)) (double (:dispersion oracle-model))))})

(defn- stat-key-for
  "§3.3: :z-value for fixed-dispersion families (binomial/poisson/nbinomial --
  nbinomial's dispersion is fixed at 1.0, the overdispersion is absorbed into
  theta instead), :t-value for estimated-dispersion families -- MUST mirror
  `sut/glm`'s own `:estimated-dispersion?` branch, confirmed directly against
  a fitted model, not assumed."
  [family]
  (if (#{:binomial :poisson :nbinomial} family) :z-value :t-value))

;; ---- build + main ----

(def ^:private out-dir "test/resources/regression")

;; ---- Tier 1 predict (new-data fit/stderr/CI), basic-lm ----
;; `basic-lm-variants`' own R fits use literal-embedded vectors in the formula
;; (confirmed via a `model.frame` column-name probe: the column names come
;; back as the literal vector's own printed text, not `x1`/`x2`), so R's own
;; `predict(object, newdata=...)` can't resolve a newdata there. Refit each
;; variant on an equivalent named data.frame purely for prediction --
;; `:refit-diff` cross-validates each refit's own coefficients against
;; fastmath's, exactly like `:diff` already does for the literal-embedded fit.

(def ^:private basic-new-x1 [-10.0 -0.3])
(def ^:private basic-new-x2 [1.0 2.0])
(def ^:private basic-new-off [10.0 11.0])

(defn- basic-lm-predict-refits
  "[name rfit offset?] for each basic-lm variant, fit against one shared named
  data.frame `df`."
  []
  (let [df (rr/r `(data.frame :y ~ys :x1 ~xss-1 :x2 ~xss-2 :w ~weights :off ~offset))]
    [["basic-weights-intercept" (rr/r `(lm (formula y (* x1 x2)) :weights w :data ~df)) false]
     ["basic-weights-intercept-transformer"
      (rr/r `(lm (formula y (+ (exp x1) (* x1 x2))) :weights w :data ~df)) false]
     ["basic-weights-intercept-offset" (rr/r `(lm (formula y (* x1 x2)) :weights w :offset off :data ~df)) true]
     ["basic-weights-no-intercept" (rr/r `(lm (formula y (+ 0 (* x1 x2))) :weights w :data ~df)) false]
     ["basic-weights-no-intercept-offset"
      (rr/r `(lm (formula y (+ 0 (* x1 x2))) :weights w :offset off :data ~df)) true]
     ["basic-no-weights-intercept" (rr/r `(lm (formula y (* x1 x2)) :data ~df)) false]
     ["basic-no-weights-intercept-offset" (rr/r `(lm (formula y (* x1 x2)) :offset off :data ~df)) true]
     ["basic-no-weights-no-intercept" (rr/r `(lm (formula y (+ 0 (* x1 x2))) :data ~df)) false]
     ["basic-no-weights-no-intercept-offset"
      (rr/r `(lm (formula y (+ 0 (* x1 x2))) :offset off :data ~df)) true]]))

(defn build-basic-lm!
  "Fit each basic-* R form, write its reference file, cross-validate against
  fastmath's own `lm` on the same inputs/options. Also computes a new-data
  `:extra :predict` block (fit/stderr/confidence-interval) via a named-data.frame
  refit, see `basic-lm-predict-refits`. Returns a report seq."
  []
  (let [refits (into {} (map (fn [[n r o]] [n [r o]]) (basic-lm-predict-refits)))]
    (doall
     (for [[name r-form r-formula-str options] basic-lm-variants]
       (let [rlm (rr/r r-form)
             oracle (lm-oracle rlm (not (false? (:intercept? options))))
             [rfit off?] (refits name)
             newdata-form (if off?
                            `(data.frame :x1 ~basic-new-x1 :x2 ~basic-new-x2 :off ~basic-new-off)
                            `(data.frame :x1 ~basic-new-x1 :x2 ~basic-new-x2))
             predict (predict-lm-oracle rfit newdata-form)
             refit-diff (cross-validate-lm ys xss options
                                            (:model (lm-oracle rfit (not (false? (:intercept? options))))))
             oracle (assoc oracle :extra {:predict predict})
             file (spit-single-variant! out-dir name
                                         {:r-formula r-formula-str
                                          :r-packages ["stats" "car" "moments"]
                                          :data oracle})
             diff (cross-validate-lm ys xss options (:model oracle))]
         {:name name :file file :diff diff :refit-diff refit-diff})))))

;; ---- Tier 1 predict, dummy-data + gaussian-glm ----
;; Same landmine/fix as basic-lm: `dummy-glm-variants`/`gaussian-glm-variants`
;; fit on literal-embedded vectors, so a dedicated named-data.frame refit is
;; needed for `predict(..., newdata=...)` to resolve columns by name. `f`
;; below is a plain helper fn (not a macro), so its own backtick is evaluated
;; fresh per call -- no nested-syntax-quote hazard.

(defn- dummy-glm-predict-refits
  "[key rfit] per dummy-data variant, fit against one shared named data.frame
  `df` (response `nys`, matching `dummy-glm-variants`' own response). No
  variant uses offset."
  []
  (let [df (rr/r `(data.frame :y ~nys :x1 ~xss-1 :x2 ~xss-2 :w ~weights))
        f (fn [rhs family-r & {:keys [weights?]}]
            (rr/r (with-glm-control
                   `(glm (formula y ~rhs) :family ~family-r :data ~df
                         ~@(when weights? [:weights 'w])))))]
    [[:call-0 (f '(* x1 x2) 'binomial :weights? true)]
     [:call-1 (f '(* x1 x2) '(binomial :link probit) :weights? true)]
     [:call-2 (f '(* x1 x2) '(binomial :link cloglog) :weights? true)]
     [:call-3 (f '(* x1 x2) '(binomial :link log) :weights? true)]
     [:call-4 (f '(* x1 x2) 'binomial)]
     [:call-5 (f '(* x1 x2) '(binomial :link probit))]
     [:call-6 (f '(* x1 x2) '(binomial :link cloglog))]
     [:call-7 (f '(* x1 x2) '(binomial :link log))]
     [:call-8 (f '(+ 0 (* x1 x2)) 'binomial :weights? true)]
     [:call-9 (f '(+ 0 (* x1 x2)) '(binomial :link probit) :weights? true)]
     [:call-10 (f '(+ 0 (* x1 x2)) '(binomial :link cloglog) :weights? true)]
     [:call-11 (f '(+ 0 (* x1 x2)) 'binomial)]
     [:call-12 (f '(+ 0 (* x1 x2)) '(binomial :link probit))]
     [:call-13 (f '(+ 0 (* x1 x2)) '(binomial :link cloglog))]
     [:call-14 (f '(+ (exp x1) (* x1 x2)) 'binomial)]]))

(defn- gaussian-glm-predict-refits
  "[key rfit offset?] per gaussian-glm variant, fit against one shared named
  data.frame `df` (response `ys`)."
  []
  (let [df (rr/r `(data.frame :y ~ys :x1 ~xss-1 :x2 ~xss-2 :w ~weights :off ~offset))
        f (fn [rhs & {:keys [weights? offset?]}]
            (rr/r (with-glm-control
                   `(glm (formula y ~rhs) :data ~df
                         ~@(when weights? [:weights 'w])
                         ~@(when offset? [:offset 'off])))))]
    [[:call-0 (f '(* x1 x2) :weights? true) false]
     [:call-1 (f '(* x1 x2) :weights? true :offset? true) true]
     [:call-2 (f '(+ 0 (* x1 x2)) :weights? true) false]
     [:call-3 (f '(+ 0 (* x1 x2)) :weights? true :offset? true) true]
     [:call-4 (f '(* x1 x2)) false]
     [:call-5 (f '(* x1 x2) :offset? true) true]
     [:call-6 (f '(+ 0 (* x1 x2))) false]
     [:call-7 (f '(+ 0 (* x1 x2)) :offset? true) true]]))

(defn- estimated-dispersion-for
  "Inverse of `stat-key-for`'s split, for `predict-glm-oracle`'s
  `estimated-dispersion?` arg."
  [family]
  (= (stat-key-for family) :t-value))

(defn build-dummy-data!
  "Fit each dummy-data (binomial GLM) R form, cross-validate against fastmath's
  `glm`, attempt the §3.3 dose fix's oracle where fastmath's model qualifies
  (`(count (:coefficients model)) > 1`), write the one multi-variant file."
  []
  (let [refits (into {} (dummy-glm-predict-refits))
        newdata-form `(data.frame :x1 ~basic-new-x1 :x2 ~basic-new-x2)
        variants
        (doall
         (for [[key r-form r-formula-str options] dummy-glm-variants]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm nys xss options)
                 diff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle (refits key) newdata-form
                                             (estimated-dispersion-for (:family options)))
                 extra (cond-> {:predict predict} dose (assoc :dose dose))
                 data (assoc oracle :extra extra)]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff diff :dose dose})))]
    (spit-multi-variant! out-dir "dummy-data" {:variants variants})
    variants))

(defn build-gaussian-glm!
  "Fit each gaussian-glm R form, cross-validate against fastmath's `glm`,
  capture `statmod::qresid` (§3.3's `:qres?` flag) and attempt the dose oracle,
  write the one multi-variant file."
  []
  (let [refits (into {} (map (fn [[k r o]] [k [r o]]) (gaussian-glm-predict-refits)))
        variants
        (doall
         (for [[key r-form r-formula-str options] gaussian-glm-variants]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 oracle (glm-oracle rglm intercept? :t-value)
                 model (fit-glm ys xss options)
                 diff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 [rfit off?] (refits key)
                 newdata-form (if off?
                                `(data.frame :x1 ~basic-new-x1 :x2 ~basic-new-x2 :off ~basic-new-off)
                                `(data.frame :x1 ~basic-new-x1 :x2 ~basic-new-x2))
                 predict (predict-glm-oracle rfit newdata-form true)
                 extra (cond-> {:quantile-residuals (qresid-oracle rglm) :predict predict}
                        dose (assoc :dose dose))
                 data (assoc oracle :extra extra)]
             {:key key :r-formula r-formula-str :r-packages ["stats" "statmod" "MASS"] :data data
              :diff diff :dose dose})))]
    (spit-multi-variant! out-dir "gaussian-glm" {:variants variants})
    variants))

;; ---- Tier 2, gestation-data (LM) ----

(defn build-gestation-data!
  "GLMsData::gestation, two lm-tests calls (weighted, unweighted) on the same
  raw columns -- §1.3's file-level :inputs, shared across both :call-N blocks.
  Both variants are fit via `:data gestation` (named columns), so `rlm` can be
  used directly for the `:extra :predict` new-data check -- no refit needed
  (unlike Tier 1)."
  []
  (utils/data 'gestation)
  (let [ds (rr/r->clj 'gestation)
        columns (dataset-columns ds [:Weight :Age :Births])
        weight (:Weight columns) age (:Age columns) births (:Births columns)
        newdata-form `(data.frame :Age [26.0 40.0])
        variant-specs
        [[:call-0 '(lm (formula Weight Age) :weights Births :data gestation)
          "lm(Weight ~ Age, weights=Births, data=gestation)" {:weights births}]
         [:call-1 '(lm (formula Weight Age) :data gestation)
          "lm(Weight ~ Age, data=gestation)" nil]]
        variants
        (doall
         (for [[key r-form r-formula-str options] variant-specs]
           (let [rlm (rr/r r-form)
                 oracle (lm-oracle rlm true)
                 diff (cross-validate-lm weight age options (:model oracle))
                 predict (predict-lm-oracle rlm newdata-form)
                 data (assoc oracle :extra {:predict predict})]
             {:key key :r-formula r-formula-str :r-packages ["stats" "car" "moments"]
              :data data :diff diff})))]
    (spit-multi-variant! out-dir "gestation-data"
                         {:inputs-comment "raw GLMsData::gestation columns (Weight, Age, Births)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 2, dental-data (LM) ----

(defn build-dental-data!
  "GLMsData::dental. Single lm-tests call; :Indus is a 2-level factor with a
  trivial binary mapping (reference level :Ind -> 0.0) done inline -- no
  contrast/dummy needed for a single 2-level factor (§2.3 stays deferred).
  Predict newdata: Sugar=[10,50], Indus=[Ind,NonInd] (one row per level)."
  []
  (utils/data 'dental)
  (let [ds (rr/r->clj 'dental)
        columns (dataset-columns ds [:DMFT :Sugar :Indus])
        {:keys [DMFT Sugar Indus]} columns
        indus-bin (mapv #(if (= % :Ind) 0.0 1.0) Indus)
        sugar-indus (mapv * Sugar indus-bin)
        xss (map vector Sugar indus-bin sugar-indus)
        r-form '(lm (formula DMFT (* Sugar Indus)) :data dental)
        rlm (rr/r r-form)
        oracle (lm-oracle rlm true)
        diff (cross-validate-lm DMFT xss nil (:model oracle))
        predict (predict-lm-oracle rlm `(data.frame :Sugar [10.0 50.0] :Indus ["Ind" "NonInd"]))
        file (spit-single-variant! out-dir "dental-data"
                                   {:r-formula "lm(DMFT ~ Sugar*Indus, data=dental)"
                                    :r-packages ["stats" "car" "moments"]
                                    :data (assoc oracle :inputs {:columns columns} :extra {:predict predict})})]
    {:file file :diff diff}))

;; ---- Tier 2, cheese-data (LM) ----

(defn build-cheese-data!
  "GLMsData::cheese. Single lm-tests call; full 3-way factorial expansion of
  Acetic*log(H2S)*Lactic (3 main + 3 pairwise + 1 triple = 7 columns) -- plain
  elementwise products, no contrast/dummy needed (all-numeric interaction).
  Predict newdata: two new (Acetic, H2S, Lactic) rows within the training
  range."
  []
  (utils/data 'cheese)
  (let [ds (rr/r->clj 'cheese)
        columns (dataset-columns ds [:Taste :Acetic :H2S :Lactic])
        {:keys [Taste Acetic H2S Lactic]} columns
        lH2S (mapv m/log H2S)
        acetic-lH2S (mapv * Acetic lH2S)
        acetic-lactic (mapv * Acetic Lactic)
        lH2S-lactic (mapv * lH2S Lactic)
        triple (mapv * Acetic lH2S Lactic)
        xss (map vector Acetic lH2S Lactic acetic-lH2S acetic-lactic lH2S-lactic triple)
        r-form '(lm (formula Taste (* Acetic (log H2S) Lactic)) :data cheese)
        rlm (rr/r r-form)
        oracle (lm-oracle rlm true)
        diff (cross-validate-lm Taste xss nil (:model oracle))
        predict (predict-lm-oracle rlm `(data.frame :Acetic [300.0 150.0] :H2S [5000.0 100.0] :Lactic [1.3 1.0]))
        file (spit-single-variant! out-dir "cheese-data"
                                   {:r-formula "lm(Taste ~ Acetic*log(H2S)*Lactic, data=cheese)"
                                    :r-packages ["stats" "car" "moments"]
                                    :data (assoc oracle :inputs {:columns columns} :extra {:predict predict})})]
    {:file file :diff diff}))

;; ---- Tier 2, turbunes-data (GLM binomial) ----

(defn build-turbunes-data!
  "GLMsData::turbines. First Tier 2 GLM -- two response-encoding conventions
  (ratio+weights vs. cbind-pairs) across 4 variants, still no categorical
  encoding (§2.3 stays deferred). Predict newdata: Hours=[1000,3000]."
  []
  (utils/data 'turbines)
  (let [ds (rr/r->clj 'turbines)
        columns (dataset-columns ds [:Fissures :Turbines :Hours])
        {:keys [Fissures Turbines Hours]} columns
        ratio (mapv / Fissures Turbines)
        turb-fiss-diff (mapv - Turbines Fissures)
        newdata-form `(data.frame :Hours [1000.0 3000.0])
        variant-specs
        [[:call-0 '(glm (formula (/ Fissures Turbines) Hours)
                        :family binomial :data turbines :weights Turbines)
          "glm(Fissures/Turbines ~ Hours, family=binomial, weights=Turbines, data=turbines)"
          ratio Hours {:family :binomial :weights Turbines}]
         [:call-1 '(glm (formula (cbind Fissures (- Turbines Fissures)) Hours)
                        :family binomial :data turbines)
          "glm(cbind(Fissures, Turbines-Fissures) ~ Hours, family=binomial, data=turbines)"
          (mapv vector Fissures turb-fiss-diff) Hours {:family :binomial}]
         [:call-2 '(glm (formula (/ Fissures Turbines) Hours)
                        :family (binomial :link "probit") :data turbines :weights Turbines)
          "glm(Fissures/Turbines ~ Hours, family=binomial(link=probit), weights=Turbines, data=turbines)"
          ratio Hours {:family :binomial :weights Turbines :link :probit}]
         [:call-3 '(glm (formula (/ Fissures Turbines) Hours)
                        :family (binomial :link "cloglog") :data turbines :weights Turbines)
          "glm(Fissures/Turbines ~ Hours, family=binomial(link=cloglog), weights=Turbines, data=turbines)"
          ratio Hours {:family :binomial :weights Turbines :link :cloglog}]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict}) dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "turbunes-data"
                         {:inputs-comment "raw GLMsData::turbines columns (Fissures, Turbines, Hours)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 2, germ-data (GLM binomial/quasi-binomial) ----

(defn build-germ-data!
  "GLMsData::germ. First deftest exercising :no-ll? (quasi-binomial family --
  R's logLik/AIC/BIC aren't well-defined for quasi-families). Both Seeds and
  Extract are 2-level factors -- inline binary mapping suffices, same as
  dental-data's Indus (§2.3 still deferred). Predict newdata: one row per
  Seeds/Extract level combo (Seeds=[OA75,OA73], Extract=[Cucumber,Bean])."
  []
  (utils/data 'germ)
  (let [ds (rr/r->clj 'germ)
        columns (dataset-columns ds [:Germ :Total :Seeds :Extract])
        {:keys [Germ Total Seeds Extract]} columns
        seeds-oa75 (mapv #(if (= % :OA75) 1.0 0.0) Seeds)
        extract-cucumber (mapv #(if (= % :Cucumber) 1.0 0.0) Extract)
        es (mapv * extract-cucumber seeds-oa75)
        ratio (mapv / Germ Total)
        newdata-form `(data.frame :Seeds ["OA75" "OA73"] :Extract ["Cucumber" "Bean"])
        variant-specs
        [[:call-0 '(glm (formula (/ Germ Total) (+ Seeds Extract))
                        :family binomial :data germ :weights Total)
          "glm(Germ/Total ~ Seeds+Extract, family=binomial, weights=Total, data=germ)"
          ratio (mapv vector seeds-oa75 extract-cucumber) {:family :binomial :weights Total}]
         [:call-1 '(glm (formula (/ Germ Total) (* Seeds Extract))
                        :family binomial :data germ :weights Total)
          "glm(Germ/Total ~ Seeds*Extract, family=binomial, weights=Total, data=germ)"
          ratio (mapv vector seeds-oa75 extract-cucumber es) {:family :binomial :weights Total}]
         [:call-2 '(glm (formula (/ Germ Total) (* Seeds Extract))
                        :family quasibinomial :data germ :weights Total)
          "glm(Germ/Total ~ Seeds*Extract, family=quasibinomial, weights=Total, data=germ)"
          ratio (mapv vector seeds-oa75 extract-cucumber es)
          {:family :quasi-binomial :weights Total :no-ll? true}]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict}) dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "germ-data"
                         {:inputs-comment "raw GLMsData::germ columns (Germ, Total, Seeds, Extract)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 2, mammary-data (GLM binomial cloglog) ----

(defn build-mammary-data!
  "GLMsData::mammary. First deftest using the `:intercept` xss sentinel
  (intercept-only model) and (redundantly, matching the original) explicit
  :compare-deviance-null? false. All-numeric predictors, no encoding needed."
  []
  (utils/data 'mammary)
  (let [ds (rr/r->clj 'mammary)
        columns (dataset-columns ds [:N.Outgrowths :N.Assays :N.Cells])
        n-outgrowths (:N.Outgrowths columns)
        n-assays (:N.Assays columns)
        n-cells (:N.Cells columns)
        ratio (mapv / n-outgrowths n-assays)
        lncells (mapv m/log n-cells)
        ;; new N.Cells (offset source)/N.Outgrowths, within the training range
        newdata-off-only `(data.frame :N.Cells [50.0 90.0])
        newdata-with-nout `(data.frame :N.Outgrowths [5.0 10.0] :N.Cells [50.0 90.0])
        newdata-no-offset `(data.frame :N.Cells [50.0 90.0])
        variant-specs
        [[:call-0 '(glm (formula (/ N.Outgrowths N.Assays) (offset (log N.Cells)))
                        :family (binomial :link cloglog) :weights N.Assays :data mammary)
          "glm(N.Outgrowths/N.Assays ~ 1 + offset(log(N.Cells)), family=binomial(link=cloglog), weights=N.Assays, data=mammary)"
          ratio :intercept {:family :binomial :link :cloglog :weights n-assays :offset lncells}
          newdata-off-only]
         [:call-1 '(glm (formula (/ N.Outgrowths N.Assays) N.Outgrowths)
                        :family (binomial :link cloglog) :weights N.Assays :data mammary :offset (log N.Cells))
          "glm(N.Outgrowths/N.Assays ~ N.Outgrowths + offset(log(N.Cells)), family=binomial(link=cloglog), weights=N.Assays, data=mammary)"
          ratio n-outgrowths {:family :binomial :link :cloglog :weights n-assays :offset lncells}
          newdata-with-nout]
         [:call-2 '(glm (formula (/ N.Outgrowths N.Assays) (+ 0 N.Outgrowths))
                        :family (binomial :link cloglog) :weights N.Assays :data mammary :offset (log N.Cells))
          "glm(N.Outgrowths/N.Assays ~ 0 + N.Outgrowths + offset(log(N.Cells)), family=binomial(link=cloglog), weights=N.Assays, data=mammary)"
          ratio n-outgrowths {:family :binomial :link :cloglog :weights n-assays :offset lncells
                              :intercept? false :compare-deviance-null? false}
          newdata-with-nout]
         [:call-3 '(glm (formula (/ N.Outgrowths N.Assays) (log N.Cells))
                        :family (binomial :link cloglog) :weights N.Assays :data mammary)
          "glm(N.Outgrowths/N.Assays ~ log(N.Cells), family=binomial(link=cloglog), weights=N.Assays, data=mammary)"
          ratio lncells {:family :binomial :link :cloglog :weights n-assays}
          newdata-no-offset]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options newdata-form] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict}) dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "mammary-data"
                         {:inputs-comment "raw GLMsData::mammary columns (N.Outgrowths, N.Assays, N.Cells)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 2, pock-data (GLM poisson/quasi-poisson/negative-binomial) ----

(defn build-pock-data!
  "GLMsData::pock. First glm.nb/:nbinomial exercise. `:call-2`'s fastmath
  `:nbinomial-theta` option is a literal (§1's :options-stay-in-source
  decision) matching R's own glm.nb fit -- captured here only to keep the
  generator/test literal consistent, not stored in EDN. The file-level :extra
  holds R's own glm.nb theta, for a standalone comparison against fastmath's
  iterative `glm-nbinomial` (not tied to any single glm-tests-edn call)."
  []
  (utils/data 'pock)
  (let [ds (rr/r->clj 'pock)
        columns (dataset-columns ds [:Count :Dilution])
        count-col (:Count columns)
        dilution (:Dilution columns)
        l2dilution (mapv #(/ (m/log %) (m/log 2.0)) dilution)
        nbinomial-theta 9.892894299757403
        newdata-form `(data.frame :Dilution [4.0 10.0])
        variant-specs
        [[:call-0 '(glm (formula Count (log2 Dilution)) :family poisson :data pock)
          "glm(Count ~ log2(Dilution), family=poisson, data=pock)"
          count-col l2dilution {:family :poisson}]
         [:call-1 '(glm (formula Count (log2 Dilution)) :family quasipoisson :data pock)
          "glm(Count ~ log2(Dilution), family=quasipoisson, data=pock)"
          count-col l2dilution {:family :quasi-poisson :no-ll? true}]
         [:call-2 '(glm.nb (formula Count (log2 Dilution)) :data pock)
          "glm.nb(Count ~ log2(Dilution), data=pock)"
          count-col l2dilution {:family :nbinomial :nbinomial-theta nbinomial-theta}]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict}) dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))
        r-theta (double (first (:theta (rr/r->clj (rr/r '(glm.nb (formula Count (log2 Dilution)) :data pock))))))]
    (spit-multi-variant! out-dir "pock-data"
                         {:inputs-comment "raw GLMsData::pock columns (Count, Dilution)"
                          :inputs {:columns columns}
                          :extra-comment "R glm.nb theta, for a standalone comparison against fastmath's glm-nbinomial"
                          :extra {:nbinomial-theta r-theta}
                          :variants variants})
    {:variants variants :r-theta r-theta}))

;; ---- Tier 2, hcrabs-data (GLM quasi-poisson/negative-binomial) ----

(defn build-hcrabs-data!
  "GLMsData::hcrabs. First multi-level (3-/4-level) categorical exercise --
  activates §2.3's contrast/dummy-based `one-hot`. Also repeats pock-data's
  quasi-poisson (:no-ll?) / glm.nb / standalone-theta pattern."
  []
  (utils/data 'hcrabs)
  (let [ds (rr/r->clj 'hcrabs)
        columns (dataset-columns ds [:Sat :Wt :Width :Spine :Col])
        {:keys [Sat Wt Width Spine Col]} columns
        logwt (mapv m/log Wt)
        logwidth (mapv m/log Width)
        {spine-noneok :NoneOK spine-oneok :OneOK} (one-hot Spine [:NoneOK :OneOK])
        {col-dm :DM col-lm :LM col-m :M} (one-hot Col [:DM :LM :M])
        full-xss (map vector logwt logwidth spine-noneok spine-oneok col-dm col-lm col-m)
        nbinomial-theta 0.9580286019527684
        ;; new-data predict points: Wt=[2000,4000], Width=[25,30],
        ;; Spine=[NoneOK,OneOK], Col=[DM,M] -- `one-hot-mapping` encodes the
        ;; new Spine/Col values with the SAME reference level `one-hot` used above.
        new-wt [2000.0 4000.0] new-width [25.0 30.0]
        new-logwt (mapv m/log new-wt) new-logwidth (mapv m/log new-width)
        new-spine (map (one-hot-mapping Spine [:NoneOK :OneOK]) [:NoneOK :OneOK])
        new-col (map (one-hot-mapping Col [:DM :LM :M]) [:DM :M])
        new-full-xss (map (fn [lw lwd sp co] (vec (concat [lw lwd] sp co)))
                          new-logwt new-logwidth new-spine new-col)
        newdata-form-full `(data.frame :Wt ~new-wt :Width ~new-width
                                        :Spine ["NoneOK" "OneOK"] :Col ["DM" "M"])
        newdata-form-wt `(data.frame :Wt ~new-wt)
        variant-specs
        [[:call-0 '(glm (formula Sat (+ (log Wt) (log Width) Spine Col))
                        :family quasipoisson :data hcrabs)
          "glm(Sat ~ log(Wt)+log(Width)+Spine+Col, family=quasipoisson, data=hcrabs)"
          Sat full-xss {:family :quasi-poisson :no-ll? true} newdata-form-full new-full-xss]
         [:call-1 '(glm (formula Sat (+ (log Wt))) :family quasipoisson :data hcrabs)
          "glm(Sat ~ log(Wt), family=quasipoisson, data=hcrabs)"
          Sat logwt {:family :quasi-poisson :no-ll? true} newdata-form-wt new-logwt]
         [:call-2 '(glm.nb (formula Sat (+ (log Wt))) :data hcrabs)
          "glm.nb(Sat ~ log(Wt), data=hcrabs)"
          Sat logwt {:family :nbinomial :nbinomial-theta nbinomial-theta} newdata-form-wt new-logwt]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options newdata-form new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict :newdata-xss (vec new-xss)})
                        dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))
        r-theta (double (first (:theta (rr/r->clj (rr/r '(glm.nb (formula Sat (+ (log Wt))) :data hcrabs))))))]
    (spit-multi-variant! out-dir "hcrabs-data"
                         {:inputs-comment "raw GLMsData::hcrabs columns (Sat, Wt, Width, Spine, Col)"
                          :inputs {:columns columns}
                          :extra-comment "R glm.nb theta, for a standalone comparison against fastmath's glm-nbinomial"
                          :extra {:nbinomial-theta r-theta}
                          :variants variants})
    {:variants variants :r-theta r-theta}))

;; ---- Tier 2, lime-data (GLM gamma/inverse-gaussian) ----

(defn build-lime-data!
  "GLMsData::lime. `Origin` is one-hot'd over ALL its levels (no reference
  held back), intentional in the original -- `one-hot`'s reference-lookup
  degrades to nil there, giving every real level its own full dummy column,
  matching the original's construction exactly. `:call-0` (gamma) additionally
  captures `statmod::qresid` for a 5e-3-tolerance quantile-residuals check."
  []
  (utils/data 'lime)
  (let [ds (rr/r->clj 'lime)
        columns (dataset-columns ds [:Foliage :DBH :Origin])
        {:keys [Foliage DBH Origin]} columns
        logdbh (mapv m/log DBH)
        {origin-natural :Natural origin-planted :Planted} (one-hot Origin [:Natural :Planted])
        origin-natural-logdbh (mapv * origin-natural logdbh)
        origin-planted-logdbh (mapv * origin-planted logdbh)
        xss (map vector origin-natural origin-planted logdbh origin-natural-logdbh origin-planted-logdbh)
        ;; new-data predict points, shared by both variants (same formula/xss
        ;; shape, only family/link differ): DBH=[10,25], Origin=[Natural,Planted]
        new-dbh [10.0 25.0]
        new-logdbh (mapv m/log new-dbh)
        new-origin (map (one-hot-mapping Origin [:Natural :Planted]) [:Natural :Planted])
        new-xss (mapv (fn [[on op] ld] [on op ld (* on ld) (* op ld)]) new-origin new-logdbh)
        newdata-form `(data.frame :DBH ~new-dbh :Origin ["Natural" "Planted"])
        variant-specs
        [[:call-0 '(glm (formula Foliage (* Origin (log DBH))) :family (Gamma :link log) :data lime)
          "glm(Foliage ~ Origin*log(DBH), family=Gamma(link=log), data=lime)"
          Foliage xss {:family :gamma :link :log}]
         [:call-1 '(glm (formula Foliage (* Origin (log DBH))) :family (inverse.gaussian :link log) :data lime)
          "glm(Foliage ~ Origin*log(DBH), family=inverse.gaussian(link=log), data=lime)"
          Foliage xss {:family :inverse-gaussian :link :log}]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 qres (when (= key :call-0) (qresid-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 extra (cond-> {:predict predict :newdata-xss new-xss}
                        dose (assoc :dose dose) qres (assoc :quantile-residuals qres))
                 data (assoc oracle :extra extra)]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS" "statmod"] :data data
              :diff vdiff :dose dose :qres qres})))]
    (spit-multi-variant! out-dir "lime-data"
                         {:inputs-comment "raw GLMsData::lime columns (Foliage, DBH, Origin)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 2, perm-data (GLM inverse-gaussian) ----

(defn build-perm-data!
  "GLMsData::perm. First factor-x-factor interaction (Mach x Day), exercising
  the new §2.3 `mult-columns` helper alongside `one-hot`. `Day` must be
  R-side-refactored (matching the original's `factor(Day)` step) before both
  fitting and column extraction, or R would fit it as a continuous term."
  []
  (utils/data 'perm)
  (rr/r '(<- ($ perm Day) (factor ($ perm Day))))
  (let [ds (rr/r->clj 'perm)
        columns (dataset-columns ds [:Perm :Mach :Day])
        {:keys [Perm Mach Day]} columns
        {:keys [xss cs1 cs2]} (mult-columns Mach Day)
        r-form '(glm (formula Perm (* Mach Day)) :family (inverse.gaussian :link log) :data perm)
        rglm (rr/r (with-glm-control r-form))
        stat-key (stat-key-for :inverse-gaussian)
        oracle (glm-oracle rglm true stat-key)
        model (fit-glm Perm xss {:family :inverse-gaussian :link :log})
        vdiff (cross-validate-glm model (:model oracle))
        dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
        ;; new-data predict points: Mach=[B,C], Day=["3","7"] -- encoded with
        ;; the exact same cs1/cs2 (non-reference level) order `mult-columns`
        ;; used above, via `one-hot-mapping`.
        new-mach [:B :C] new-day [:3 :7]
        new-mach-oh (map (one-hot-mapping Mach cs1) new-mach)
        new-day-oh (map (one-hot-mapping Day cs2) new-day)
        new-xss (mapv (fn [moh doh] (vec (concat moh doh (for [d doh m moh] (* (double m) (double d))))))
                      new-mach-oh new-day-oh)
        newdata-form `(data.frame :Mach ["B" "C"] :Day ["3" "7"])
        predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for :inverse-gaussian))
        data (assoc oracle :extra (cond-> {:predict predict :newdata-xss new-xss} dose (assoc :dose dose)))]
    (spit-single-variant! out-dir "perm-data"
                          {:r-formula "glm(Perm ~ Mach*Day, family=inverse.gaussian(link=log), data=perm)"
                           :r-packages ["stats" "MASS"]
                           :data (assoc data :inputs {:columns columns})})
    {:diff vdiff :dose dose}))

;; ---- Tier 2, yieldden-data (GLM gamma) ----

(defn build-yieldden-data!
  "GLMsData::yieldden. `Var` needs R-side factor conversion and `YD =
  Yield*Dens` needs R-side computation, matching the original init-r steps,
  before both fitting and column extraction. The standalone quantile-
  residuals comparison refits call-1's design under :inverse-gaussian (not
  either glm-tests call's own family) -- file-level :extra, same pattern as
  pock/hcrabs-data's standalone theta comparisons."
  []
  (utils/data 'yieldden)
  (rr/r '(<- ($ yieldden Var) (factor ($ yieldden Var))))
  (rr/r '(<- ($ yieldden YD) (with yieldden (* Yield Dens))))
  (let [ds (rr/r->clj 'yieldden)
        columns (dataset-columns ds [:YD :Dens :Var])
        {:keys [YD Dens Var]} columns
        rdens (mapv / (repeat 1.0) Dens)
        {var2 :2 var3 :3} (one-hot Var [:2 :3])
        dens-var2 (mapv * Dens var2)
        dens-var3 (mapv * Dens var3)
        rdens-var2 (mapv * rdens var2)
        rdens-var3 (mapv * rdens var3)
        full-xss (map vector Dens rdens var2 var3 dens-var2 dens-var3 rdens-var2 rdens-var3)
        simple-xss (map vector Dens rdens var2 var3)
        ;; new-data predict points, shared by both variants: Dens=[5,20], Var=[2,3]
        new-dens [5.0 20.0]
        new-rdens (mapv #(/ 1.0 (double %)) new-dens)
        new-var (map (one-hot-mapping Var [:2 :3]) [:2 :3])
        new-simple-xss (mapv (fn [d rd [v2 v3]] [d rd v2 v3]) new-dens new-rdens new-var)
        new-full-xss (mapv (fn [d rd [v2 v3]] [d rd v2 v3 (* d v2) (* d v3) (* rd v2) (* rd v3)])
                           new-dens new-rdens new-var)
        newdata-form `(data.frame :Dens ~new-dens :Var ["2" "3"])
        variant-specs
        [[:call-0 '(glm "YD ~ (Dens + I(1/Dens)) * Var" :family (Gamma :link inverse) :data yieldden)
          "glm(YD ~ (Dens + I(1/Dens)) * Var, family=Gamma(link=inverse), data=yieldden)"
          YD full-xss {:family :gamma :link :inverse} new-full-xss]
         [:call-1 '(glm "YD ~ Dens + I(1/Dens) + Var" :family (Gamma :link inverse) :data yieldden)
          "glm(YD ~ Dens + I(1/Dens) + Var, family=Gamma(link=inverse), data=yieldden)"
          YD simple-xss {:family :gamma :link :inverse} new-simple-xss]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm true stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (assoc oracle :extra (cond-> {:predict predict :newdata-xss new-xss}
                                             dose (assoc :dose dose)))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))
        rig (rr/r (with-glm-control
                   '(glm "YD ~ Dens + I(1/Dens) + Var" :family inverse.gaussian :data yieldden)))
        qres (qresid-oracle rig)]
    (spit-multi-variant! out-dir "yieldden-data"
                         {:inputs-comment "raw GLMsData::yieldden columns (YD, Dens, Var)"
                          :inputs {:columns columns}
                          :extra-comment "R qresid on an inverse-gaussian refit of call-1's design, for a standalone comparison"
                          :extra {:quantile-residuals qres}
                          :variants variants})
    {:variants variants :qres-count (count qres)}))

;; ---- Tier 2, lungcap-data (LM, 7 calls) ----

(defn build-lungcap-data!
  "GLMsData::lungcap, the last Tier 2 dataset. `Smoke` is already a plain
  0/1 numeric column -- R's `Smokef` factor coincides with it numerically
  under default treatment coding, so no encoding is needed for it. `Gender`
  needs both a simple 0/1 numeric encoding (calls 1/3/4, reference F=0,
  matching R's alphabetical-default contrast) and a full all-levels one-hot
  (call 2's no-intercept model, same pattern as lime-data's Origin)."
  []
  (utils/data 'lungcap)
  (rr/r "lungcap$Smokef <- factor(lungcap$Smoke,levels=c(0, 1), labels=c(\"Non-smoker\",\"Smoker\"))")
  (let [ds (rr/r->clj 'lungcap)
        columns (dataset-columns ds [:FEV :Age :Ht :Gender :Smoke])
        {:keys [FEV Age Ht Gender Smoke]} columns
        logfev (mapv m/log FEV)
        gender-numeric (mapv #(if (= % :F) 0.0 1.0) Gender)
        {gender-f :F gender-m :M} (one-hot Gender [:F :M])
        ht-smoke (mapv * Ht Smoke)
        ;; new-data predict points, shared by all 7 calls: Age=[8,15], Ht=[55,65],
        ;; Gender=[F,M], Smoke=[0,1] (R side: Smokef=[Non-smoker,Smoker])
        new-age [8.0 15.0] new-ht [55.0 65.0]
        new-gender-numeric [0.0 1.0]
        new-gender-f [1.0 0.0] new-gender-m [0.0 1.0]
        new-smoke [0.0 1.0]
        new-ht-smoke (mapv * new-ht new-smoke)
        newdata-form `(data.frame :Age ~new-age :Ht ~new-ht :Gender ["F" "M"]
                                   :Smokef ["Non-smoker" "Smoker"])
        variant-specs
        [["call-0" '(lm (formula (log FEV) (+ Age Ht Gender Smokef)) :data lungcap)
          "lm(log(FEV) ~ Age+Ht+Gender+Smokef, data=lungcap)"
          (map vector Age Ht gender-numeric Smoke) nil
          (mapv vector new-age new-ht new-gender-numeric new-smoke)]
         ["call-1" '(lm (formula (log FEV) (+ 0 Age Ht Gender Smokef)) :data lungcap)
          "lm(log(FEV) ~ 0+Age+Ht+Gender+Smokef, data=lungcap)"
          (map vector Age Ht gender-f gender-m Smoke) {:intercept? false}
          (mapv vector new-age new-ht new-gender-f new-gender-m new-smoke)]
         ["call-2" '(lm (formula (log FEV) (+ Ht Gender Smokef)) :data lungcap)
          "lm(log(FEV) ~ Ht+Gender+Smokef, data=lungcap)"
          (map vector Ht gender-numeric Smoke) nil
          (mapv vector new-ht new-gender-numeric new-smoke)]
         ["call-3" '(lm (formula (log FEV) (+ Gender Smokef)) :data lungcap)
          "lm(log(FEV) ~ Gender+Smokef, data=lungcap)"
          (map vector gender-numeric Smoke) nil
          (mapv vector new-gender-numeric new-smoke)]
         ["call-4" '(lm (formula (log FEV) (+ Age Smokef)) :data lungcap)
          "lm(log(FEV) ~ Age+Smokef, data=lungcap)"
          (map vector Age Smoke) nil
          (mapv vector new-age new-smoke)]
         ["call-5" '(lm (formula (log FEV) (+ Ht Smokef)) :data lungcap)
          "lm(log(FEV) ~ Ht+Smokef, data=lungcap)"
          (map vector Ht Smoke) nil
          (mapv vector new-ht new-smoke)]
         ["call-6" '(lm (formula (log FEV) (* Ht Smokef)) :data lungcap)
          "lm(log(FEV) ~ Ht*Smokef, data=lungcap)"
          (map vector Ht Smoke ht-smoke) nil
          (mapv vector new-ht new-smoke new-ht-smoke)]]
        variants
        (doall
         (for [[key r-form r-formula-str vxs options new-xss] variant-specs]
           (let [rlm (rr/r r-form)
                 oracle (lm-oracle rlm (not (false? (:intercept? options))))
                 vdiff (cross-validate-lm logfev vxs options (:model oracle))
                 predict (predict-lm-oracle rlm newdata-form)
                 data (assoc oracle :extra {:predict predict :newdata-xss new-xss})]
             {:key (keyword key) :r-formula r-formula-str :r-packages ["stats" "car" "moments"]
              :data data :diff vdiff})))]
    (spit-multi-variant! out-dir "lungcap-data"
                         {:inputs-comment "raw GLMsData::lungcap columns (FEV, Age, Ht, Gender, Smoke)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 3, counts2x2-data (GLM poisson/binomial, hand-built data.frame) ----

(defn build-counts2x2-data!
  "A hand-built 2x2 contingency table (not an R package dataset) -- `Att`/`Inc`
  built via R's `gl()` exactly as the original test constructed them. Two
  poisson log-linear calls (independence, saturated) plus one binomial call
  modeling P(Att=Against) weighted by Counts -- exercises `one-hot` without
  `mult-columns` (each call needs a different sub-slice of the same two
  factors' encoding, not the full interaction design at once)."
  []
  (let [rcounts2x2 (rr/r '(data.frame :Counts [263 258 151 222]
                                      :Att (gl 2 2 4 :labels ["For" "Against"])
                                      :Inc (gl 2 1 4 :labels ["High" "Low"])))
        ds (rr/r->clj rcounts2x2)
        columns (dataset-columns ds [:Counts :Att :Inc])
        {:keys [Counts Att Inc]} columns
        att-against (:Against (one-hot Att [:Against]))
        inc-low (:Low (one-hot Inc [:Low]))
        att-inc-low (mapv * att-against inc-low)
        ;; new-data predict points: (Att,Inc)=[(For,High),(Against,Low)]
        new-att-against (map (comp first (one-hot-mapping Att [:Against])) [:For :Against])
        new-inc-low (map (comp first (one-hot-mapping Inc [:Low])) [:High :Low])
        newdata-form `(data.frame :Att ["For" "Against"] :Inc ["High" "Low"])
        variant-specs
        [[:call-0 `(glm (formula Counts (+ Att Inc)) :family poisson :data ~rcounts2x2)
          "glm(Counts ~ Att+Inc, family=poisson, data=counts2x2)"
          Counts (map vector att-against inc-low) {:family :poisson :no-analysis? true}
          (mapv vector new-att-against new-inc-low)]
         [:call-1 `(glm (formula Counts (* Att Inc)) :family poisson :data ~rcounts2x2)
          "glm(Counts ~ Att*Inc, family=poisson, data=counts2x2)"
          Counts (map vector att-against inc-low att-inc-low) {:family :poisson :no-analysis? true}
          (mapv vector new-att-against new-inc-low (map * new-att-against new-inc-low))]
         [:call-2 `(glm (formula (ifelse (== Att "Against") 1 0) Inc)
                        :family binomial :weight Counts :data ~rcounts2x2)
          "glm(ifelse(Att==\"Against\",1,0) ~ Inc, family=binomial, weight=Counts, data=counts2x2)"
          att-against inc-low {:family :binomial :weights Counts} new-inc-low]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict :newdata-xss new-xss})
                        dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "counts2x2-data"
                         {:inputs-comment "hand-built 2x2 contingency table (Counts, Att, Inc via R's gl())"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 3, kstones-data (GLM poisson/binomial, 3-way factorial) ----

(defn build-kstones-data!
  "GLMsData::kstones. Ten variants exhausting every sub-model of the full
  Size*Method*Outcome 3-way factorial (each a different column order/subset,
  matching the original's hand-picked `tc/select-columns` slices exactly),
  plus a final binomial call on the Outcome response weighted by Counts."
  []
  (utils/data 'kstones)
  (let [ds (rr/r->clj 'kstones)
        columns (dataset-columns ds [:Counts :Size :Method :Outcome])
        {:keys [Counts Size Method Outcome]} columns
        size-small (:Small (one-hot Size [:Small]))
        method-b (:B (one-hot Method [:B]))
        outcome-success (:Success (one-hot Outcome [:Success]))
        size-method (mapv * size-small method-b)
        size-outcome (mapv * size-small outcome-success)
        method-outcome (mapv * method-b outcome-success)
        size-method-outcome (mapv * size-small method-b outcome-success)
        ;; new-data predict points: (Size,Method,Outcome)=[(Small,B,Success),(Small,A,Failure)]
        new-size (map (comp first (one-hot-mapping Size [:Small])) [:Small :Small])
        new-method (map (comp first (one-hot-mapping Method [:B])) [:B :A])
        new-outcome (map (comp first (one-hot-mapping Outcome [:Success])) [:Success :Failure])
        new-size-method (mapv * new-size new-method)
        new-size-outcome (mapv * new-size new-outcome)
        new-method-outcome (mapv * new-method new-outcome)
        new-size-method-outcome (mapv * new-size new-method new-outcome)
        newdata-form `(data.frame :Size ["Small" "Small"] :Method ["B" "A"] :Outcome ["Success" "Failure"])
        variant-specs
        [[:call-0 '(glm (formula Counts (+ Size Method Outcome)) :family poisson :data kstones)
          "glm(Counts ~ Size+Method+Outcome, family=poisson, data=kstones)"
          Counts (map vector size-small method-b outcome-success) {:family :poisson}
          (mapv vector new-size new-method new-outcome)]
         [:call-1 '(glm (formula Counts (+ (* Size Method) Outcome)) :family poisson :data kstones)
          "glm(Counts ~ Size*Method+Outcome, family=poisson, data=kstones)"
          Counts (map vector size-small method-b outcome-success size-method) {:family :poisson}
          (mapv vector new-size new-method new-outcome new-size-method)]
         [:call-2 '(glm (formula Counts (+ (* Size Outcome) Method)) :family poisson :data kstones)
          "glm(Counts ~ Size*Outcome+Method, family=poisson, data=kstones)"
          Counts (map vector size-small outcome-success method-b size-outcome) {:family :poisson}
          (mapv vector new-size new-outcome new-method new-size-outcome)]
         [:call-3 '(glm (formula Counts (+ (* Outcome Method) Size)) :family poisson :data kstones)
          "glm(Counts ~ Outcome*Method+Size, family=poisson, data=kstones)"
          Counts (map vector outcome-success method-b size-small method-outcome) {:family :poisson}
          (mapv vector new-outcome new-method new-size new-method-outcome)]
         [:call-4 '(glm "Counts ~ Size * (Method + Outcome)" :family poisson :data kstones)
          "glm(Counts ~ Size*(Method+Outcome), family=poisson, data=kstones)"
          Counts (map vector size-small method-b outcome-success size-method size-outcome) {:family :poisson}
          (mapv vector new-size new-method new-outcome new-size-method new-size-outcome)]
         [:call-5 '(glm "Counts ~ Method * (Outcome + Size)" :family poisson :data kstones)
          "glm(Counts ~ Method*(Outcome+Size), family=poisson, data=kstones)"
          Counts (map vector method-b outcome-success size-small method-outcome size-method) {:family :poisson}
          (mapv vector new-method new-outcome new-size new-method-outcome new-size-method)]
         [:call-6 '(glm "Counts ~ Outcome * (Method + Size)" :family poisson :data kstones)
          "glm(Counts ~ Outcome*(Method+Size), family=poisson, data=kstones)"
          Counts (map vector outcome-success method-b size-small method-outcome size-outcome) {:family :poisson}
          (mapv vector new-outcome new-method new-size new-method-outcome new-size-outcome)]
         [:call-7 '(glm "Counts ~ Size * Method * Outcome - Size:Method:Outcome" :family poisson :data kstones)
          "glm(Counts ~ Size*Method*Outcome - Size:Method:Outcome, family=poisson, data=kstones)"
          Counts (map vector size-small method-b outcome-success size-method size-outcome method-outcome)
          {:family :poisson :no-analysis? true}
          (mapv vector new-size new-method new-outcome new-size-method new-size-outcome new-method-outcome)]
         [:call-8 '(glm (formula Counts (* Size Method Outcome)) :family poisson :data kstones)
          "glm(Counts ~ Size*Method*Outcome, family=poisson, data=kstones)"
          Counts (map vector size-small method-b outcome-success size-method size-outcome method-outcome
                       size-method-outcome)
          {:family :poisson :no-analysis? true}
          (mapv vector new-size new-method new-outcome new-size-method new-size-outcome new-method-outcome
                new-size-method-outcome)]
         [:call-9 '(glm (formula (ifelse (== Outcome "Success") 1 0) (* Size Method))
                        :family binomial :data kstones :weights Counts)
          "glm(ifelse(Outcome==\"Success\",1,0) ~ Size*Method, family=binomial, weights=Counts, data=kstones)"
          outcome-success (map vector size-small method-b size-method) {:family :binomial :weights Counts}
          (mapv vector new-size new-method new-size-method)]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (predict-glm-oracle rglm newdata-form (estimated-dispersion-for (:family options)))
                 data (cond-> (assoc oracle :extra {:predict predict :newdata-xss new-xss})
                        dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "kstones-data"
                         {:inputs-comment "raw GLMsData::kstones columns (Counts, Size, Method, Outcome)"
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 3, danishlc-data (GLM poisson, offset + City*Age interaction + poly) ----

(defn build-danishlc-data!
  "GLMsData::danishlc. R-side setup (Rate, ordered Age releveling, abbreviated
  City, AgeNum) replicated from the original test's init-r, matching exactly
  -- including `(base/options :contrasts [\"contr.treatment\"
  \"contr.treatment\"])`, WITHOUT which R's default `contr.poly` for ordered
  factors makes `Age`'s encoding disagree with this file's dummy coding
  (confirmed directly: cross-validation diverges at the 1e0 scale without it,
  drops to 1e-14 with it). `call-3` uses `reg/polynomial-transformer` (a
  2-degree fitted polynomial `:transformer`, superseding the old `add-poly`
  helper -- see the topic note's Continuation) alongside `mult-columns`'
  City*Age interaction; unlike `add-poly`, the transformer closure can be
  evaluated at new AgeNum values, so `call-3` also gets a newdata-predict check."
  []
  (base/options :contrasts ["contr.treatment" "contr.treatment"])
  (utils/data 'danishlc)
  (rr/r '(<- ($ danishlc Rate) (* 1000 (/ ($ danishlc Cases) ($ danishlc Pop)))))
  (rr/r '(<- ($ danishlc Age) (ordered ($ danishlc Age) :levels ["40-54", "55-59", "60-64", "65-69", "70-74", ">74"])))
  (rr/r '(<- ($ danishlc City) (abbreviate ($ danishlc City) 1)))
  (rr/r '(<- ($ danishlc AgeNum) (rep [40, 55, 60, 65, 70, 75] 4)))
  (let [ds (rr/r->clj 'danishlc)
        columns (dataset-columns ds [:Cases :Pop :Age :City :AgeNum])
        {:keys [Cases Pop Age City AgeNum]} columns
        loffset (v/log (seq Pop))
        {:keys [xss cs1 cs2]} (mult-columns City Age)
        age-oh (one-hot Age [:55-59 :60-64 :65-69 :70-74 :>74])
        age-xss (map vector (:55-59 age-oh) (:60-64 age-oh) (:65-69 age-oh) (:70-74 age-oh) (:>74 age-oh))
        ;; `polynomial-transformer` is fit once here (on the training AgeNum)
        ;; and applied to each row via `:transformer`, same convention as
        ;; `basic-weights-intercept-transformer` -- `call-3`'s xss stays the
        ;; RAW AgeNum column (like `call-2`'s), not a pre-expanded matrix.
        poly-transformer (reg/polynomial-transformer AgeNum 2)
        ;; new-data predict points. City=[H,V], Age=[60-64,70-74],
        ;; AgeNum=[50,68], Pop=[1000,2000] (offset)
        new-pop [1000.0 2000.0]
        new-loffset (mapv m/log new-pop)
        ;; City is a plain character column (R's `abbreviate()` returns a
        ;; character vector, not a factor) -- its Clojure values are strings,
        ;; unlike Age's keyword-valued ordered factor.
        new-city-oh (map (one-hot-mapping City cs1) ["H" "V"])
        new-age-oh (map (one-hot-mapping Age cs2) [:60-64 :70-74])
        new-call0-xss (mapv (fn [coh aoh] (vec (concat coh aoh (for [a aoh c coh] (* (double c) (double a))))))
                            new-city-oh new-age-oh)
        new-call1-xss (map (one-hot-mapping Age [:55-59 :60-64 :65-69 :70-74 :>74]) [:60-64 :70-74])
        new-agenum [50.0 68.0]
        ;; call-2 and call-3 (poly) share this exact raw-input shape --
        ;; call-3's `polynomial-transformer` is applied internally by `predict`,
        ;; after offset extraction, same as at fitting time.
        new-agenum-xss (mapv vector new-loffset new-agenum)
        newdata-form `(data.frame :City ["H" "V"] :Age ["60-64" "70-74"] :AgeNum ~new-agenum :Pop ~new-pop)
        variant-specs
        [[:call-0 '(glm (formula Cases (+ (offset (log Pop)) (* City Age))) :family poisson :data danishlc)
          "glm(Cases ~ offset(log(Pop)) + City*Age, family=poisson, data=danishlc)"
          Cases xss {:family :poisson :offset loffset :no-analysis? true}
          (mapv cons new-loffset new-call0-xss)]
         [:call-1 '(glm (formula Cases (+ (offset (log Pop)) Age)) :family poisson :data danishlc)
          "glm(Cases ~ offset(log(Pop)) + Age, family=poisson, data=danishlc)"
          Cases age-xss {:family :poisson :offset loffset}
          (mapv cons new-loffset new-call1-xss)]
         [:call-2 '(glm (formula Cases (+ (offset (log Pop)) AgeNum)) :family poisson :data danishlc)
          "glm(Cases ~ offset(log(Pop)) + AgeNum, family=poisson, data=danishlc)"
          Cases AgeNum {:family :poisson :offset loffset}
          new-agenum-xss]
         [:call-3 '(glm (formula Cases (+ (offset (log Pop)) (poly AgeNum 2))) :family poisson :data danishlc)
          "glm(Cases ~ offset(log(Pop)) + poly(AgeNum,2), family=poisson, data=danishlc)"
          Cases AgeNum {:family :poisson :offset loffset :transformer poly-transformer}
          new-agenum-xss]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (when new-xss (predict-glm-oracle rglm newdata-form
                                                            (estimated-dispersion-for (:family options))))
                 data (cond-> (assoc oracle :extra {})
                        predict (assoc-in [:extra :predict] predict)
                        new-xss (assoc-in [:extra :newdata-xss] new-xss)
                        dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))]
    (spit-multi-variant! out-dir "danishlc-data"
                         {:inputs-comment (str "raw GLMsData::danishlc columns (Cases, Pop, Age, City, AgeNum) "
                                               "after R-side Rate/Age/City/AgeNum setup")
                          :inputs {:columns columns}
                          :variants variants})
    variants))

;; ---- Tier 3, deposit-data (GLM binomial, one-hot over all levels + poly + newdata predict + dose) ----

(defn build-deposit-data!
  "GLMsData::deposit, the last Tier 3 dataset. `Insecticide` is one-hot'd over
  ALL 3 levels (same all-levels-degenerate-reference pattern as lime-data's
  Origin/lungcap-data's call-2), letting each call pick whichever k-1 (with
  intercept) or k (without) dummy columns it needs. `call-3` uses
  `reg/polynomial-transformer` wrapped in a small `deposit-poly-transformer`
  (poly-expands only the `log(Deposit)` slot, passes the two Insecticide
  dummies through unchanged) -- supersedes the old `add-poly` helper, see the
  topic note's Continuation; also gets its own newdata-predict check for
  exactly that reason (`add-poly` couldn't). The file-level :extra
  holds R's `predict(..., newdata=...)` response for each Insecticide level
  over a shared 100-point Deposit grid, for a standalone comparison against
  fastmath's own `sut/glm` model built directly in the test (not through
  glm-tests-edn, which only ever cross-validates the fitted :model/:analysis
  fields, not predictions on new data). The `sut/dose` checks against
  call-1's fitted model use literal expected values already in the original
  test (not sourced from R at all) -- unaffected by this migration, kept as
  literals at the call site."
  []
  (utils/data 'deposit)
  (let [ds (rr/r->clj 'deposit)
        columns (dataset-columns ds [:Killed :Number :Deposit :Insecticide])
        {:keys [Killed Number Deposit Insecticide]} columns
        ratio (mapv / Killed Number)
        logdep (mapv m/log Deposit)
        ins-oh (one-hot Insecticide [:A :B :C])
        {ins-a :A ins-b :B ins-c :C} ins-oh
        poly-transformer (reg/polynomial-transformer logdep 2)
        deposit-poly-transformer (fn [[ld insb insc]] (into (poly-transformer [ld]) [insb insc]))
        ;; new-data predict points for call-3 only: Deposit=[3,6], Insecticide=[B,C]
        new-deposit-poly [3.0 6.0]
        new-logdep-poly (mapv m/log new-deposit-poly)
        new-insb-poly [1.0 0.0] new-insc-poly [0.0 1.0]
        new-call3-xss (mapv vector new-logdep-poly new-insb-poly new-insc-poly)
        newdata-form-poly `(data.frame :Deposit ~new-deposit-poly :Insecticide ["B" "C"])
        variant-specs
        [[:call-0 '(glm (formula (/ Killed Number) (+ Deposit Insecticide))
                        :family binomial :weights Number :data deposit)
          "glm(Killed/Number ~ Deposit+Insecticide, family=binomial, weights=Number, data=deposit)"
          ratio (map vector Deposit ins-b ins-c) {:family :binomial :weights Number} nil]
         [:call-1 '(glm (formula (/ Killed Number) (+ 0 Deposit Insecticide))
                        :family binomial :weights Number :data deposit)
          "glm(Killed/Number ~ 0+Deposit+Insecticide, family=binomial, weights=Number, data=deposit)"
          ratio (map vector Deposit ins-a ins-b ins-c) {:family :binomial :weights Number :intercept? false} nil]
         [:call-2 '(glm (formula (/ Killed Number) (+ 0 (log Deposit) Insecticide))
                        :family binomial :weights Number :data deposit)
          "glm(Killed/Number ~ 0+log(Deposit)+Insecticide, family=binomial, weights=Number, data=deposit)"
          ratio (map vector logdep ins-a ins-b ins-c) {:family :binomial :weights Number :intercept? false} nil]
         [:call-3 '(glm (formula (/ Killed Number) (+ (poly (log Deposit) 2) Insecticide))
                        :family binomial :weights Number :data deposit)
          "glm(Killed/Number ~ poly(log(Deposit),2)+Insecticide, family=binomial, weights=Number, data=deposit)"
          ratio (map vector logdep ins-b ins-c) {:family :binomial :weights Number :transformer deposit-poly-transformer}
          new-call3-xss]]
        variants
        (doall
         (for [[key r-form r-formula-str vys vxs options new-xss] variant-specs]
           (let [rglm (rr/r (with-glm-control r-form))
                 intercept? (not (false? (:intercept? options)))
                 stat-key (stat-key-for (:family options))
                 oracle (glm-oracle rglm intercept? stat-key)
                 model (fit-glm vys vxs options)
                 vdiff (cross-validate-glm model (:model oracle))
                 dose (when (> (count (:coefficients model)) 1) (dose-oracle rglm))
                 predict (when new-xss
                           (predict-glm-oracle rglm newdata-form-poly (estimated-dispersion-for (:family options))))
                 data (cond-> (assoc oracle :extra {})
                        predict (assoc-in [:extra :predict] predict)
                        new-xss (assoc-in [:extra :newdata-xss] new-xss)
                        dose (assoc-in [:extra :dose] dose))]
             {:key key :r-formula r-formula-str :r-packages ["stats" "MASS"] :data data
              :diff vdiff :dose dose})))
        rglm0 (rr/r (with-glm-control (second (first variant-specs))))
        newdata (m/slice-range 2.0 8.0 100)
        ;; fit/stderr/confidence-interval (fixed-dispersion binomial ->
        ;; estimated-dispersion?=false), not just :fit as before -- upgrades
        ;; the original bare-fit newdata comparison to the same {:fit :stderr
        ;; :confidence-interval} shape every other dataset's :extra :predict uses.
        predict-at (fn [level]
                     (predict-glm-oracle rglm0 `(data.frame :Deposit ~newdata :Insecticide ~level) false))
        predict-a (predict-at "A") predict-b (predict-at "B") predict-c (predict-at "C")]
    (spit-multi-variant! out-dir "deposit-data"
                         {:inputs-comment "raw GLMsData::deposit columns (Killed, Number, Deposit, Insecticide)"
                          :inputs {:columns columns}
                          :extra-comment (str "R predict(..., type=\"response\"), fit/stderr/CI, on call-0's model "
                                              "over a shared 100-point Deposit grid, one block per Insecticide "
                                              "level, for a standalone comparison against fastmath's own newdata "
                                              "predictions (built directly in the test, not through glm-tests-edn)")
                          :extra {:newdata newdata :predict-a predict-a :predict-b predict-b :predict-c predict-c}
                          :variants variants})
    {:variants variants}))

(defn -main [& _]
  (let [lm-report (build-basic-lm!)]
    (doseq [{:keys [name diff]} lm-report]
      (println (format "%-40s %s" name diff)))
    (println "basic-lm:" (count lm-report) "files written"))
  (let [dummy-report (build-dummy-data!)]
    (doseq [{:keys [key diff dose]} dummy-report]
      (println (format "dummy-data %-10s %s dose=%s" key diff dose)))
    (println "dummy-data: 1 file," (count dummy-report) "variants"))
  (let [gauss-report (build-gaussian-glm!)]
    (doseq [{:keys [key diff dose]} gauss-report]
      (println (format "gaussian-glm %-10s %s dose=%s" key diff dose)))
    (println "gaussian-glm: 1 file," (count gauss-report) "variants"))
  (let [gestation-report (build-gestation-data!)]
    (doseq [{:keys [key diff]} gestation-report]
      (println (format "gestation-data %-10s %s" key diff)))
    (println "gestation-data: 1 file," (count gestation-report) "variants"))
  (println "dental-data" (build-dental-data!))
  (println "cheese-data" (build-cheese-data!))
  (println "turbunes-data" (build-turbunes-data!))
  (println "germ-data" (build-germ-data!))
  (println "mammary-data" (build-mammary-data!))
  (println "pock-data" (build-pock-data!))
  (println "hcrabs-data" (build-hcrabs-data!))
  (println "lime-data" (build-lime-data!))
  (println "perm-data" (build-perm-data!))
  (println "yieldden-data" (build-yieldden-data!))
  (println "lungcap-data" (build-lungcap-data!))
  (println "counts2x2-data" (build-counts2x2-data!))
  (println "kstones-data" (build-kstones-data!))
  (println "danishlc-data" (build-danishlc-data!))
  (println "deposit-data" (build-deposit-data!))
  (println "done"))
