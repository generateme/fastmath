^:kindly/hide-code
(ns core
  (:require [fastmath.dev.codox :as codox]
            [fastmath.dev.ggplot :as gg]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.core :as m]
            [fastmath.special :as special]
            [fastmath.transform :as t]
            [fastmath.kernel.variogram :as variogram]            
            [fastmath.vector :as v]))

;; # Core {.unnumbered}

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.core)
