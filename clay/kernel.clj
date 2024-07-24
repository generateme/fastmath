^{:kindly/hide-code true}
(ns kernel
  (:require [fastmath.kernel.density :as kde]
            [scicloj.kindly.v4.kind :as kind]))

^{:kindly/hide-code true}
(require '[ggplot]
         '[formulas])

;; # Kernels

;; ## Kernel density estimation (KDE)

^{:kindly/kind :kind/table :kindly/hide-code true}
(->> (for [[k {:keys [radius kernel]}] (sort-by first kde/kde-data)
           :let [r (* radius 1.1)]]
       (ggplot/->image (ggplot/function kernel {:x [(- r) r] :title (name k)})))
     (partition 2))

^{:kindly/hide-code true}
(kind/md (formulas/formulas :kde/uniform))


(kind/md "$$\\frac{1}{2}$$")


;; ## Radial Basis Functions (RBF)

;; ## Vector kernels

;; ## Variograms
