### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace. 

### Please write a short paragraph for `clay/stats.clj` book. Do not copy content of the `clay/stats.clj` file, write only changed chapter. Put text as a Clojure comment block.

- Write a paragraph for a "## Intervals and Extents" chapter.
- Describe all covered functions.
- Compare these functions

Write only asked chapter.



```clojure
;; ## Intervals and Extents
;;
;; This section provides functions to describe the spread or define specific ranges and intervals within a dataset. Measures like `span` and `iqr` offer simple indicators of total range and the spread of the central 50%. Functions ending in `-extent` provide symmetric intervals around a central measure (like mean, median, or SEM), often representing a common spread like +/- 1 standard deviation or MAD. Quantile-based functions (`percentile-extent`, `quantile-extent`, `pi`, `pi-extent`, `hpdi-extent`) define intervals based on data ranks or probability content, useful for describing specific portions of the distribution or confidence/credible intervals. The `adjacent-values` and fence-related functions (`inner-fence-extent`, `outer-fence-extent`) define boundaries commonly used in box plots and outlier detection. Finally, `ci`, `percentile-bc-extent`, and `percentile-bca-extent` calculate confidence intervals for the mean or other statistics, using methods like the t-distribution or bootstrap resampling (bias-corrected and accelerated). These functions collectively offer various perspectives on the dispersion, typical range, and extreme values of a dataset, from simple summaries to robust or inferential measures.
```

### Can you provide better formatting. I would like to see each group as a separate point.



```clojure
;; ## Intervals and Extents
;;
;; This section provides functions to describe the spread or define specific ranges and intervals within a dataset.
;;
;; *   **Basic Range:** Functions like `span` ($max - min$) and `extent` (providing $[min, max]$ and optionally the mean) offer simple measures of the total spread of the data.
;;
;; *   **Interquartile Range:** `iqr` ($Q_3 - Q_1$) specifically measures the spread of the middle 50% of the data, providing a robust alternative to the total range.
;;
;; *   **Symmetric Spread Intervals:** Functions ending in `-extent` such as `stddev-extent`, `mad-extent`, and `sem-extent` define intervals typically centered around the mean or median. They represent a range defined by adding/subtracting a multiple (usually 1) of a measure of dispersion (Standard Deviation, Median Absolute Deviation, or Standard Error of the Mean) from the central point.
;;
;; *   **Quantile-Based Intervals:** `percentile-extent`, `quantile-extent`, `pi`, `pi-extent`, and `hpdi-extent` define intervals based on quantiles or percentiles of the data. These functions capture specific ranges containing a certain percentage of the data points (e.g., the middle 95% defined by quantiles 0.025 and 0.975). `hpdi-extent` calculates the shortest interval containing a given proportion of data, based on empirical density.
;;
;; *   **Box Plot Boundaries:** `adjacent-values` (LAV, UAV) and fence functions (`inner-fence-extent`, `outer-fence-extent`) calculate specific bounds based on quartiles and multiples of the IQR. These are primarily used in box plot visualization and as a conventional method for identifying potential outliers.
;;
;; *   **Confidence and Prediction Intervals:** `ci`, `percentile-bc-extent`, and `percentile-bca-extent` provide inferential intervals. `ci` estimates a confidence interval for the population mean using the t-distribution. `percentile-bc-extent` and `percentile-bca-extent` (Bias-Corrected and Bias-Corrected Accelerated) are advanced bootstrap methods for estimating confidence intervals for statistics, offering robustness against non-normality and bias.
```

### Write latex formula for `ci`



```latex
$\bar{x} \pm t_{\alpha/2, n-1} \frac{s}{\sqrt{n}}$
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-flash-preview-04-17 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (443 1720) (1816 3948) (3986 4049))) -->
<!-- End: -->
