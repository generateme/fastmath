---
name: docstring
description: How to write a docstring what should be included in a Clojure function or namespace docstring. 
compatibility: Designed for a Clojure project.
---

# Instructions

First, analyze and understand the code. Use common math knowledge for reasoning.

## Docstring structure

Write a docstring with the following content:

* Write a single statement about what is a given function or namespace about.
* Optionally write a short paragraph about function or namespace context
* (namespace only) Write a summary about functions defined in the namespace
* (function only) Describe all input parameters in points. In case if input is a `map` type, describe additionally all keys and their meaning and default values (if defined).
* (function only) Describe all returned values and their interpretation if applicable
* (function only) Describe corner cases, exceptions, constrains 
* (function only) Link to other related functions using markdown wikilink syntax, ie: `[[reference]]`.
* Use inline code (`) for symbols, input, keywords, Clojure forms
* There is not column width contstrain. Do not force a newline in the middle of a paragraph.
* (var only) Put the docstring as a meta tag

## Formatting

A docstring is a double quoted text put after a function or namespace name. Example of a complete function with a desired docstring.

```clojure
(defn pearson-correlation
  "Calculates the Pearson product-moment correlation coefficient between two sequences.

  This function measures the linear relationship between two datasets. The coefficient value ranges from -1.0 (perfect negative linear correlation) to 1.0 (perfect positive linear correlation), with 0.0 indicating no linear correlation.

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both input sequences must contain only numbers and must have the same length.

  Returns the calculated Pearson correlation coefficient as a double. Returns `NaN` if either sequence has zero variance (i.e., all elements are the same).

  See also [[correlation]] (general correlation, defaults to Pearson), [[spearman-correlation]],
  [[kendall-correlation]], [[correlation-matrix]]."
  (^double [[vs1 vs2]] (pearson-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))
```

In case of a var, docstring is put in a meta map after `def`.

```clojure
(def ^{:doc "Count equal values in both seqs. Alias for [[count==]]"} L0 count=)
```

## Gotchas

* Avoid LateX formulas.
* Avoid double quotes inside the text ("") and backslashes (\)
* Avoid following markdown formatting: headings, emphasis, images, code blocks (```), footnotes, blockquotes, html, horizontal rule.
