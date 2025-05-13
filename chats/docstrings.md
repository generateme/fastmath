### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace.

### Write a docstring for a `variation` function.



```clojure
(defn variation
  "Calculates the coefficient of variation (CV) for a sequence `vs`.

  The CV is a standardized measure of dispersion of a probability distribution or
  frequency distribution. It is defined as the ratio of the standard deviation
  to the mean:

  `CV = stddev(vs) / mean(vs)`

  This measure is unitless and allows for comparison of variability between
  datasets with different means or different units.

  Parameters:

  - `vs`: Sequence of numbers.

  Returns the calculated coefficient of variation as a double.

  Note: The CV is undefined if the mean is zero, and may be misleading if the
  mean is close to zero or if the data can take both positive and negative values.
  All values in `vs` should ideally be positive.

  See also [[stddev]], [[mean]]."
  ^double [vs]
  (let [vs (m/seq->double-array vs)]
    (/ (stddev vs)
       (mean vs))))
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-pro-exp-03-25 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (165 1050))) -->
<!-- End: -->
