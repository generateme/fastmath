^:kindly/hide-code
(ns index
  (:require [scicloj.kindly.v4.kind :as kind]))

;; # Preface {.unnumbered}

;; Documentation work in progress

;; Statuses:

;; ✓ - done

;; \+ -  partially done

;; ⇾ -  wip

;; . - awaiting

(kind/table
 {:column-names ["namespace" "clay docs" "docstrings" "included in book" "notes"]
  :row-vectors [['fastmath.core "⇾" "⇾" "✓" ""]
                ['fastmath.vector "." "." "✓" ""]
                ['fastmath.matrix "." "." "✓" ""]
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
