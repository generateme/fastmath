### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace. 

### Please write a short paragraph for `clay/stats.clj` book. Do not copy content of the `clay/stats.clj` file, write only changed chapter. Put text as a Clojure comment block.

- Write a paragraph for a "## Moments and Shape" chapter.
- Describe all covered functions (`moment`, `skewness`, `kurtosis`, `l-moment` and `l-variation`).
- Describe all types of them
- Write formulas for each of the type
- Compare these functions
- Illustrate usage with examples.



```clojure
;;
;; ## Moments and Shape
;;
;; Moments and shape statistics describe the form of a dataset's distribution, particularly its symmetry and peakedness.
;;
;; ::: {.callout-tip title="Defined functions"}
;; * `moment`
;; * `skewness`
;; * `kurtosis`
;; * `l-moment`, `l-variation`
;; :::
;;
;; ### Conventional Moments (`moment`)
;;
;; The `moment` function calculates statistical moments of a dataset. Moments can be central (around the mean), raw (around zero), or around a specified center. They can also be absolute and/or normalized.
;;
;; *   **k-th Central Moment**: $\mu_k = E[(X - \mu)^k] \approx \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^k$. Calculated when `center` is `nil` (default) and `:mean?` is `true` (default).
;; *   **k-th Raw Moment** (about origin): $\mu'_k = E[X^k] \approx \frac{1}{n} \sum_{i=1}^{n} x_i^k$. Calculated if `center` is `0.0`.
;; *   **k-th Absolute Central Moment**: $E[|X - \mu|^k] \approx \frac{1}{n} \sum_{i=1}^{n} |x_i - \bar{x}|^k$. Calculated if `:absolute?` is `true`.
;; *   **Normalization**: If `:normalize?` is `true`, the moment is divided by $\sigma^k$ (where $\sigma$ is the standard deviation), yielding a scale-invariant measure. For example, the 3rd normalized central moment is related to skewness, and the 4th to kurtosis.
;; *   If `:mean?` is `false`, the function returns the sum $\sum (x_i - c)^k$ (or sum of absolute values) instead of the mean.
;;
;; The `order` parameter specifies $k$. For example, the 2nd central moment is the variance.
;;
(utls/examples-note
  (stats/moment mpg 2) ;; Variance (almost, as it's mean of squared devs, not sample variance)
  (stats/variance mpg)
  (stats/moment mpg 3) ;; Related to skewness
  (stats/moment mpg 4 {:normalize? true}) ;; Normalized 4th central moment, related to kurtosis
  (stats/moment mpg 1 {:absolute? true :center (stats/median mpg)}) ;; Mean absolute deviation from median
  )
;;
;; ### Skewness
;;
;; Skewness measures the asymmetry of a probability distribution about its mean. `fastmath.stats/skewness` offers several types:
;;
;; *   **Moment-based (sensitive to outliers):**
;;     *   `:G1` (Default): Sample skewness based on the 3rd standardized moment, adjusted for sample bias (via Apache Commons Math).
;;         A common formula for sample skewness is $G_1 = \frac{n}{(n-1)(n-2)} \sum_{i=1}^n \left(\frac{x_i - \bar{x}}{s}\right)^3$.
;;     *   `:g1` or `:pearson`: Pearson's moment coefficient of skewness, another bias-adjusted version of the 3rd standardized moment.
;;         $g_1 = \frac{m_3}{s^3}$, where $m_3$ is the 3rd central moment and $s$ is the sample standard deviation. Fastmath applies specific sample corrections.
;;     *   `:b1`: Sample skewness coefficient, related to $g_1$. Standard $b_1 = \frac{m_3}{m_2^{3/2}}$. Fastmath applies specific sample corrections.
;; *   **Robust (less sensitive to outliers):**
;;     *   `:median` (Median Skewness / Pearson's first skewness coefficient): $S_P = 3 \frac{\text{mean} - \text{median}}{\text{stddev}}$.
;;     *   `:mode` (Pearson's second skewness coefficient): $S_K = \frac{\text{mean} - \text{mode}}{\text{stddev}}$. Mode estimation method can be specified.
;;     *   `:bowley` or `:yule` (with $u=0.25$): Based on quartiles $Q_1, Q_2, Q_3$.
;;         $S_B = \frac{(Q_3 - Q_2) - (Q_2 - Q_1)}{Q_3 - Q_1} = \frac{Q_3 + Q_1 - 2Q_2}{Q_3 - Q_1}$.
;;     *   `:yule` (Yule's coefficient): Generalization of Bowley's, using quantiles $Q_u, Q_{0.5}, Q_{1-u}$.
;;         $S_Y(u) = \frac{(Q_{1-u} - Q_{0.5}) - (Q_{0.5} - Q_u)}{Q_{1-u} - Q_u}$.
;;     *   `:B3`: Robust measure $\frac{\text{mean} - \text{median}}{\text{mean}(|X_i - \text{median}(X)|)}$.
;;     *   `:hogg`: Based on comparing trimmed means ($U_{0.05}$: mean of top 5%, $L_{0.05}$: mean of bottom 5%, $M_{0.25}$: 25% trimmed mean).
;;         $S_H = \frac{U_{0.05} - M_{0.25}}{M_{0.25} - L_{0.05}}$.
;;     *   `:l-skewness`: $\tau_3 = \lambda_3 / \lambda_2$ (see L-moments).
;;
;; Positive skewness indicates a tail on the right side of the distribution; negative skewness indicates a tail on the left. Zero indicates symmetry.
;;
(utls/examples-note
  (stats/skewness residual-sugar)
  (stats/skewness residual-sugar :bowley)
  (stats/skewness residual-sugar :l-skewness)
  (stats/skewness residual-sugar :mode {:method :pearson})
  (stats/skewness (concat residual-sugar [-1000])) ;; Effect of an outlier
  (stats/skewness (concat residual-sugar [-1000]) :l-skewness) ;; Robustness of L-skewness
  )
;;
;; ### Kurtosis
;;
;; Kurtosis measures the "tailedness" or "peakedness" of a distribution. High kurtosis means heavy tails (more outliers) and a sharp peak (leptokurtic); low kurtosis means light tails and a flatter peak (platykurtic). `fastmath.stats/kurtosis` offers several types:
;;
;; *   **Moment-based (sensitive to outliers):**
;;     *   `:G2` (Default): Sample kurtosis (Fisher's definition, not excess), adjusted for sample bias (via Apache Commons Math). For a normal distribution, this is approximately 3.
;;     *   `:g2` or `:excess`: Sample excess kurtosis ($K_{excess} = K - 3$). For a normal distribution, this is approximately 0. Fastmath applies specific sample corrections.
;;     *   `:kurt`: Kurtosis defined as $g_2 + 3$. Similar to `:G2`.
;;     *   `:b2`: Sample kurtosis related to $m_4/m_2^2$. For a normal distribution, this is approximately 3. Fastmath applies specific sample corrections.
;; *   **Robust (less sensitive to outliers):**
;;     *   `:geary` (Geary's 'g'): $g = \frac{\text{mean}(|X_i - \bar{X}|)}{\text{population_stddev}(X)}$. Normal $\approx \sqrt{2/\pi} \approx 0.798$. Lower values indicate leptokurtosis.
;;     *   `:moors`: Based on octiles $E_i$ (quantiles $i/8$). $M_0 = \frac{(E_7-E_5) + (E_3-E_1)}{E_6-E_2}$. Centered by subtracting $1.233$ (Moors' constant for normality). Normal $\approx 0$.
;;     *   `:crow` (Crow-Siddiqui): Based on quantiles $Q_\alpha, Q_{1-\alpha}, Q_\beta, Q_{1-\beta}$. $CS(\alpha, \beta) = \frac{Q_{1-\alpha} - Q_{\alpha}}{Q_{1-\beta} - Q_{\beta}}$. Centered by subtracting $2.906$. Normal $\approx 0$.
;;     *   `:hogg`: Based on trimmed means $U_p$ (mean of top $p\%$) and $L_p$ (mean of bottom $p\%$). $K_H(\alpha, \beta) = \frac{U_{\alpha} - L_{\alpha}}{U_{\beta} - L_{\beta}}$. Centered by subtracting $2.585$. Normal $\approx 0$.
;;     *   `:l-kurtosis`: $\tau_4 = \lambda_4 / \lambda_2$ (see L-moments). Normal $\approx 0.1226$.
;;
(utls/examples-note
  (stats/kurtosis residual-sugar)
  (stats/kurtosis residual-sugar :excess)
  (stats/kurtosis residual-sugar :l-kurtosis)
  (stats/kurtosis residual-sugar :moors)
  (stats/kurtosis (concat residual-sugar [-1000 1000])) ;; Effect of outliers
  (stats/kurtosis (concat residual-sugar [-1000 1000]) :l-kurtosis) ;; Robustness of L-kurtosis
  )
;;
;; ### L-moments (`l-moment`, `l-variation`)
;;
;; L-moments are summary statistics analogous to conventional moments but are computed from linear combinations of order statistics (sorted data). They are more robust to outliers and provide better estimates for small samples compared to conventional moments.
;;
;; *   `l-moment vs order`: Calculates the L-moment of a specific `order`.
;;     *   $\lambda_1$: L-location (identical to the mean).
;;     *   $\lambda_2$: L-scale (a measure of dispersion).
;;     *   Higher orders relate to shape.
;; *   Trimmed L-moments (TL-moments) can be calculated by specifying `:s` (left trim) and `:r` (right trim) options.
;; *   L-moment Ratios: If `:ratio? true`, normalized L-moments are returned.
;;     *   $\tau_3 = \lambda_3 / \lambda_2$: Coefficient of L-skewness (same as `(stats/skewness vs :l-skewness)`).
;;     *   $\tau_4 = \lambda_4 / \lambda_2$: Coefficient of L-kurtosis (same as `(stats/kurtosis vs :l-kurtosis)`).
;; *   `l-variation vs`: Calculates the coefficient of L-variation, $\tau_2 = \lambda_2 / \lambda_1$. This is a dimensionless measure of dispersion, analogous to the coefficient of variation.
;;
;; L-moments often provide more reliable inferences about the underlying distribution shape, especially when data may contain outliers or come from heavy-tailed distributions.
;;
(utls/examples-note
  ;; L-location (mean) and L-scale for mpg
  (stats/l-moment mpg 1)
  (stats/mean mpg)
  (stats/l-moment mpg 2)

  ;; L-skewness and L-kurtosis for residual-sugar
  (stats/l-moment residual-sugar 3 {:ratio? true})
  (stats/skewness residual-sugar :l-skewness)
  (stats/l-moment residual-sugar 4 {:ratio? true})
  (stats/kurtosis residual-sugar :l-kurtosis)

  ;; Coefficient of L-variation for hp
  (stats/l-variation hp)
  )
;;
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-pro-exp-03-25 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (578 9154))) -->
<!-- End: -->
