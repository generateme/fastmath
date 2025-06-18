^:kindly/hide-code
(ns index
  (:require [scicloj.kindly.v4.kind :as kind]))

;; # Preface {.unnumbered}

;; Documentation work in progress

;; ::: {.callout-warning}
;; This notebook is written with the support of Gemini LLM models:
;; 
;; * `gemini-2.5-pro-exp-03-25`
;; * `gemini-2.5-flash-preview-04-17`
;;
;; I did my best to verify the output of LLMs however I don't guarantee absence of the model hallucinations or incorrectnesses.
;; :::

;; ## Status

;; ✓ - done

;; \+ -  partially done

;; ⇾ -  wip

;; . - awaiting

(kind/table
 {:column-names ["namespace" "clay docs" "docstrings" "included in book" "notes"]
  :row-vectors [['fastmath.core "✓" "✓" "✓" ""]
                ['fastmath.vector "✓" "✓" "✓" ""]
                ['fastmath.matrix "✓" "✓" "✓" ""]
                ['fastmath.random "+" "." "✓" "new functions added"]
                ['fastmath.stats "✓" "✓" "✓" "experimental, LLM based"]
                ['fastmath.stats.bootstrap "✓" "✓" "✓" "experimental, LLM based"]
                ['fastmath.polynomials "." "." "✓" ""]
                ['fastmath.special "✓" "." "✓" ""]
                ['fastmath.calculus "+" "." "✓" "Clerk version exists"]
                ['fastmath.solver "." "." "+" ""]
                ['fastmath.interpolation "⇾" "." "✓" ""]
                ['fastmath.kernel "." "." "+" ""]
                ['fastmath.optimization "." "." "✓" ""]
                ['fastmath.transform "⇾" "." "✓" ""]
                ['fastmath.signal "." "." "." "refactor required"]
                ['fastmath.ml.regression "." "." "✓" ""]
                ['fastmath.ml.clustering "." "." "✓" ""]
                ['fastmath.complex "✓" "✓" "✓" ""]
                ['fastmath.quaternions "✓" "✓" "✓" ""]
                ['fastmath.distance "." "." "." ""]
                ['fastmath.easings "." "." "✓" ""]
                ['fastmath.grid "." "." "." ""]
                ['fastmath.fields "." "." "✓" ""]
                ['fastmath.curves "." "." "." ""]
                ['fastmath.efloat "." "." "✓" ""]]})
