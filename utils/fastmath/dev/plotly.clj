(ns fastmath.dev.plotly
  (:require [fastmath.dev.ggplot :as gg]))

(def c1 "#819BB4FF")
(def c2 "#30588CFF")
(def c3 "#72A684FF")

(defn scatter3d-data
  ([vs nm] (scatter3d-data vs nm gg/color-main))
  ([vs nm color] (scatter3d-data vs nm color 5))
  ([vs nm color size]
   {:x (map first vs)
    :y (map second vs)
    :z (map last vs)
    :name nm
    :type :scatter3d
    :mode :markers
    :opacity 0.9
    :marker {:size size
             :color color}}))

(defn line3d-data
  ([vs nm] (line3d-data vs nm gg/color-main))
  ([vs nm color] (line3d-data vs nm color 5))
  ([vs nm color size]
   {:x (map first vs)
    :y (map second vs)
    :z (map last vs)
    :name nm
    :type :scatter3d
    :mode :lines
    :opacity 0.9
    :line {:width size
           :color color}}))
