(ns ggplot
  (:require [clojisr.v1.r :as r]
            [clojisr.v1.applications.plotting :as rplot]
            [tablecloth.api :as tc]
            [fastmath.core :as m]
            [fastmath.stats :as stats]))

(r/require-r '[ggplot2 :as gg]
             '[paletteer :as pal]
             '[base])

(def palette ["#000e60" "#3c5fff", "#526dff", "#637bff", "#7189ff", "#7d97ff",
            "#87a5ff", "#91b3ff", "#9ac2ff", "#a2d0ff", "#aadfff", "#b1eeff"])

(defn ->file
  ([obj] (->file obj nil))
  ([obj opts] (->file obj "img.png" opts))
  ([obj fname {:keys [width height]
               :or {width 500 height 400}}]
   (rplot/plot->file fname obj :width width :height height)))

(defn ->image
  ([obj] (->image obj nil))
  ([obj {:keys [width height]
         :or {width 500 height 400}}]
   (rplot/plot->buffered-image obj :width width :height height)))

(defmacro gg-> [& forms]  `(-> (r/r+ (gg/ggplot) (gg/theme_light)) ~@forms))
(defn aes [& aes-opts] (apply gg/aes aes-opts))
(defn line [obj data & opts] (r/r+ obj (apply gg/geom_line :data (tc/dataset data) opts)))
(defn point [obj data & opts] (r/r+ obj (apply gg/geom_point :data (tc/dataset data) opts)))
(defn histogram [obj data & opts] (r/r+ obj (apply gg/geom_histogram :data (tc/dataset data) opts)))
(defn xlim [obj [x-min x-max]] (r/r+ obj (gg/xlim (or x-min 'NA) (or x-max 'NA))))
(defn ylim [obj [y-min y-max]] (r/r+ obj (gg/ylim (or y-min 'NA) (or y-max 'NA))))
(defn title [obj title] (r/r+ obj (gg/labs :title title)))
(defn xlab [obj title] (r/r+ obj (gg/xlab title)))
(defn ylab [obj title] (r/r+ obj (gg/ylab title)))

(defn function->data
  ([f] (function->data f nil))
  ([f {:keys [x steps fname]
       :or {steps 400 fname "function"}}]
   (let [[min-x max-x] (or x [0.0 1.0])
         xs (m/slice-range min-x max-x steps)]
     (map (fn [x] {:x x :y (f x) :func fname}) xs))))

(defn functions->data
  ([fs] (functions->data fs nil))
  ([fs opts]
   (cond
     (map? fs) (functions->data (seq fs) opts)
     (fn? (first fs)) (functions->data (map-indexed (fn [id f]
                                                      [(str "function " id) f]) fs) opts)
     :else (mapcat (fn [[n f]]
                     (function->data f (assoc opts :fname n))) fs))))

;; plots

(defn function
  ([f] (function f nil))
  ([f opts]
   (gg-> (line (function->data f opts)
               :mapping (aes :x :x :y :y) :color "#000e60"))))

(defn functions
  ([fs] (functions fs nil))
  ([fs opts]
   (let [data (functions->data fs opts)
         breaks (distinct (map :func data))]
     (r/r+ (gg-> (line data :mapping (gg/aes :linetype :func :color :func :x :x :y :y)))
           (gg/scale_linetype_manual :breaks breaks :values (map inc (range (count breaks))))
           (gg/scale_color_manual :breaks breaks :values palette)))))

#_(->file (functions [["a" #(m/sin %)]
                      ["b" #(m/cos %)]
                      ["c" m/sinc]
                      ["d" #(m/sinh %)]
                      ["e" #(m/cosh %)]
                      ["f" #(m/sech %)]] {:x [-2 2]}))

(defn function+scatter
  ([f xs ys] (function+scatter f xs ys nil))
  ([f xs ys {:keys [x] :as opts}]
   (let [[x-min-p x-max-p] x
         [x-min-e x-max-e] (stats/extent xs)
         x-min (or x-min-p x-min-e)
         x-max (or x-max-p x-max-e)
         data (tc/dataset {:x xs :y ys})
         ff (if (sequential? f) functions function)]
     (r/r+ (ff f (assoc opts :x [x-min x-max]))
           (gg/geom_point :mapping (aes :x :x :y :y) :data data :color "blue" :fill "light blue"
                          :shape "circle filled" :size 3 :alpha 0.8)))))

(defn scatter
  ([xs ys]
   (-> (tc/dataset {:x xs :y ys})
       (gg/ggplot (gg/aes :x :x :y :y))
       (r/r+ (gg/geom_point :color "blue" :fill "light blue"
                            :shape "circle filled" :size 3 :alpha 0.8)))))

(defn scatters
  [data aes]
  (-> (tc/dataset data)
      (gg/ggplot (apply gg/aes aes))
      (r/r+ (gg/geom_point :shape "circle filled" :size 3))))

(defn function2d
  ([f] (function2d f nil))
  ([f {:keys [x y steps]
       :or {steps 200}}]
   (let [[min-x max-x] (or x [0.0 1.0])
         [min-y max-y] (or y [0.0 1.0])
         data (-> (for [x (m/slice-range min-x max-x steps)
                        y (m/slice-range min-y max-y steps)]
                    {:x x :y y :z (f [x y])})
                  (tc/dataset))]
     (-> data
         (gg/ggplot (gg/aes :x :x :y :y :fill :z))
         (r/r+ (gg/geom_raster :interpolate true)
               (pal/scale_fill_paletteer_c "pals::ocean.ice"))))))

(defn function2d+scatter
  ([f xs ys] (function2d+scatter f xs ys nil))
  ([f xs ys {:keys [x y] :as opts}]
   (let [[x-min-p x-max-p] x
         [x-min-e x-max-e] (stats/extent xs)
         x-min (or x-min-p x-min-e)
         x-max (or x-max-p x-max-e)
         [y-min-p y-max-p] y
         [y-min-e y-max-e] (stats/extent ys)
         y-min (or y-min-p y-min-e)
         y-max (or y-max-p y-max-e)
         data (tc/dataset {:x xs :y ys})]
     (r/r+ (function2d f (assoc opts :x [x-min x-max] :y [y-min y-max]))
           (gg/geom_point :data data :color "blue" :fill "light blue"
                          :shape "circle filled" :size 3 :alpha 0.8)))))



#_(->file (function2d+scatter (fn [[x y]] (* (m/sin (* m/TWO_PI x))
                                            (m/sin (* m/TWO_PI y)))) [0.2 0.5 0.8] [0.2 0.5 0.8]
                              {:x [0 1] :y [0 1]}))

#_(-> (function+scatter [["sinc" m/sinc]
                         ["sin" #(m/sin %)]] [-2 0 2] [-1 0 1])
      (title (r.base/expression (symbol "ABC~x^2~frac(1,2)~sqrt(2,3)")))
      (->file))

