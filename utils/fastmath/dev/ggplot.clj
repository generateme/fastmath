(ns fastmath.dev.ggplot
  "Fastmath charts"
  (:require [clojisr.v1.r :as r]
            [clojisr.v1.applications.plotting :as rplot]
            [tablecloth.api :as tc]
            [fastmath.core :as m]
            [fastmath.random :as fr]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.complex :as cplx]
            [fastmath.kernel :as kernel]
            [fastmath.signal :as signal]
            [fastmath.transform.wavelets :as wv]
            [fastmath.transform :as t]))

(r/require-r '[ggplot2 :as gg]
             '[paletteer :as pal]
             '[grDevices :as gr]
             '[sf :as sf]
             '[ozmaps :as oz]
             '[base])

(defn ->palette
  [pal]
  (if (keyword? pal)
    (if-let [n (namespace pal)]
      (str n "::" (name pal))
      (name pal))
    pal))

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn paletter-d [pal]
  (r/r->clj (pal/paletteer_d (->palette pal))))

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn paletter-c [pal n]
  (r/r->clj (pal/paletteer_c (->palette pal) n)))

(def color-main "#000e60")
(def color-light "#87a5ff")
(def palette-blue-0 [color-main "#3c5fff", "#526dff", "#637bff", "#7189ff", "#7d97ff",
                   "#87a5ff", "#91b3ff", "#9ac2ff", "#a2d0ff", "#aadfff", "#b1eeff"])
(def palette-blue-1 (repeat 12 color-main))
(def palette-blue-2 (conj (reverse '("#3c5fff", "#526dff", "#637bff", "#7189ff", "#7d97ff",
                                   "#87a5ff", "#91b3ff", "#9ac2ff", "#a2d0ff", "#aadfff", "#b1eeff"))
                        "#000e60"))


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

;; (defn ->data [data] (tc/dataset data))
;; (defmacro r+ [& forms] `(r/r+ ~@forms))
;; (defmacro gg+ [& forms]  `(r/r+ (gg/ggplot) ~@forms))
;; (defmacro ggaes+ [aes & forms]  `(r/r+ (gg/ggplot :mapping ~aes) (gg/theme_light) ~@forms))
;; (defn aes [& aes-opts] (apply gg/aes aes-opts))
;; (defn line [data & opts] (apply gg/geom_line :data (tc/dataset data) opts))
;; (defn point [data & opts] (apply gg/geom_point :data (tc/dataset data) opts))
;; (defn raster [data & opts] (apply gg/geom_raster :data (tc/dataset data) opts))
;; (defn histogram [data & opts] (apply gg/geom_histogram :data (tc/dataset data) opts))
;; (defn ribbon [data & opts] (apply gg/geom_ribbon :data (tc/dataset data) opts))
;; (defn xlim [[x-min x-max]] (gg/xlim (or x-min 'NA) (or x-max 'NA)))
;; (defn ylim [[y-min y-max]] (gg/ylim (or y-min 'NA) (or y-max 'NA)))
;; (defn title [title] (gg/labs :title title))
;; (defn xlab [title] (gg/xlab title))
;; (defn ylab [title] (gg/ylab title))

;; primitives

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn aes [& opts] (apply gg/aes opts))
(defn geom-hline [obj & opts] (r/r+ obj (apply gg/geom_hline obj opts)))

;; data processing

(defn line->data
  ([xs ys] (line->data xs ys nil))
  ([xs ys {:keys [fname]
           :or {fname "line"}}]
   (map (fn [x y] {:x x :y y :fname fname}) xs ys)))

(defn lines->data
  ([ls] (lines->data ls nil))
  ([ls opts]
   (if (map? ls)
     (lines->data (seq ls) opts)
     (mapcat (fn [[n xs ys]]
               (line->data xs ys (assoc opts :fname n))) ls))))


(defn slice
  [x steps]
  (let [[min-x max-x] (or x [0.0 1.0])]
    (m/slice-range min-x max-x (or steps 400))))

(defn function->data
  ([f] (function->data f nil))
  ([f {:keys [x steps fname]
       :or {fname "function"}}]
   (->> (slice x steps)
        (map (fn [x] {:x x :y (f x) :fname fname})))))

(defn functions->data
  ([fs] (functions->data fs nil))
  ([fs opts]
   (cond
     (map? fs) (functions->data (seq fs) opts)
     (fn? (first fs)) (functions->data (map-indexed (fn [id f]
                                                      [(str "function " id) f]) fs) opts)
     :else (mapcat (fn [[n f]]
                     (function->data f (assoc opts :fname n))) fs))))


(defn function-ci->data
  ([f] (function-ci->data f nil))
  ([f {:keys [x steps fname]
       :or {fname "function"}}]
   (->> (slice x steps)
        (map (fn [x] (let [[y ymin ymax] (f x)] {:x x :y y :ymin ymin :ymax ymax :fname fname}))))))

(defn function2d->data
  ([f] (function2d->data f nil))
  ([f {:keys [x y steps varg?]
       :or {steps 200 varg? true}}]
   (let [xs (slice x steps)
         ys (slice y steps)
         f (if varg? f (fn [[x y]] (f x y)))]
     (for [x xs y ys] {:x x :y y :z (f [x y])}))))

(defn wrap
  ^double [method ^double v]
  (case method
    :log2 (m/frac (m/log2 v))
    :log10 (m/frac (m/log10 v))
    :sin (m/norm (m/sin v) -1.0 1.0 0.0 1.0)
    :exp (- 1.0 (m/exp (- v)))
    :sigmoid (m/sigmoid v)
    (m/frac v)))

(defn ->valid [^double v] (if (m/invalid-double? v) 0.0 v))

(defn complex-function->data
  ([f] (complex-function->data f nil))
  ([f {:keys [x y steps wrap-method]
       :or {x [-4.05 4.05] y [-4.05 4.05] steps 400 wrap-method :log2}}]
   (let [xs (slice x steps)
         ys (slice y steps)]
     (for [x xs y ys
           :let [v (f (v/vec2 x y))
                 arg  (->valid (m/cnorm (cplx/arg v) m/-PI m/PI 0.0 1.0))
                 wmag (->valid (m/pow (wrap wrap-method (cplx/abs v)) 0.2))
                 sat  (->valid (m/- 1.0 (m// (m/- 1.0 wmag) 4.0)))]]
       {:x x :y y :hue arg :val wmag :sat sat}))))

(defn spectrogram->data
  ([{:keys [times freqs spectrum]}]
   (mapcat (fn [t xs]
             (map (fn [f x]
                    {:x t :y f :a x}) freqs xs)) times spectrum)))

;; plots

(def line-common {:linetype "dashed" :color color-light})

(defn add-common
  [object {:keys [title xlab ylab xlim ylim hline vline filllab ylog xlog]}]
  (r/r+ object
        (when title (gg/labs :title title))
        (when filllab (gg/labs :fill filllab))
        (when xlab (gg/xlab xlab))
        (when ylab (gg/ylab ylab))
        (when xlim (gg/xlim xlim))
        (when ylim (gg/ylim ylim))
        (when hline (apply gg/geom_hline (flatten (seq (merge line-common hline)))))
        (when vline (apply gg/geom_vline (flatten (seq (merge line-common vline)))))
        (when ylog (gg/scale_y_continuous :trans ylog))
        (when xlog (gg/scale_x_continuous :trans xlog))))

(defn function
  "Single function"
  ([f] (function f nil))
  ([f {:keys [color]
       :or {color color-main}
       :as opts}]
   (-> (r/r+ (gg/ggplot)
             (gg/theme_light)
             (gg/geom_line :data (tc/dataset (function->data f opts))
                           :mapping (gg/aes :x :x :y :y)
                           :color color))
       (add-common opts))))

(defn- narrow-range
  [mn mx p]
  (let [perc (* p (- mx mn))]
    [(+ mn perc) (- mx perc)]))

(defn density
  "Density data"
  ([xs] (density xs nil))
  ([xs {:keys [kernel bandwidth bins color fill x narrow?]
        :or {kernel :gaussian color color-main fill color-light bins :freedman-diaconis narrow? 0.05}
        :as opts}]
   (let [b (stats/estimate-bins xs bins)
         {:keys [kde mn mx]} (kernel/kernel-density kernel xs bandwidth true)
         dx (or x (if narrow? (narrow-range mn mx narrow?) [mn mx]))]
     (-> (r/r+ (gg/ggplot)
               (gg/theme_light)
               (gg/geom_histogram :mapping (gg/aes :x :x :y '(after_stat density))
                                  :data (tc/dataset {:x xs})
                                  :fill fill :bins b)
               (gg/geom_line :data (tc/dataset (function->data kde (assoc opts :x dx)))
                             :mapping (gg/aes :x :x :y :y)
                             :color color) 
               #_(gg/geom_density :color color :kernel kernel :bounds [##-Inf ##Inf]))
         (add-common opts)))))

(defn from-histogram
  "Take histogram data and convert to bars"
  ([xs] (from-histogram xs nil))
  ([xs {:keys [bins span]
        :as opts}]
   (let [h (apply stats/histogram xs bins span)]
     (-> (r/r+ (gg/ggplot)
               (gg/theme_light)
               (gg/geom_rect :mapping (gg/aes :xmin :min :xmax :max :ymin 0 :ymax :count)
                             :data (tc/dataset (:bins-maps h))
                             :fill color-light)
               (gg/xlab "x")
               (gg/ylab "counts"))
         (add-common opts)))))

(defn line
  ([ys] (line (map double (range)) ys))
  ([xs ys] (line xs ys nil))
  ([xs ys {:keys [color]
           :or {color color-main}
           :as opts}]
   (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y))
             (gg/theme_light)
             (gg/geom_line :data (tc/dataset {:x xs :y ys}) :color color))
       (add-common opts))))

(defn lines
  "Lines"
  ([ls] (lines ls nil))
  ([ls {:keys [palette legend-name linetype?]
        :or {palette palette-blue-0 legend-name "Functions" linetype? true}
        :as opts}]
   (let [data (lines->data ls opts)
         breaks (distinct (map :fname data))]
     (-> (r/r+ (gg/ggplot)
               (gg/theme_light)
               (gg/geom_line :data (tc/dataset data) :mapping (if linetype?
                                                                (gg/aes :linetype :fname :color :fname :x :x :y :y)
                                                                (gg/aes :color :fname :x :x :y :y)))
               (when linetype? (gg/scale_linetype_manual :name legend-name :breaks breaks :values (map inc (range (count breaks)))))
               (gg/scale_color_manual :name legend-name :breaks breaks :values palette))
         (add-common opts)))))

(defn lollipop
  ([xs ys] (lollipop xs ys nil))
  ([xs ys {:keys [color]
           :or {color color-main}
           :as opts}]
   (-> (r/r+ (gg/ggplot :data (tc/dataset {:x xs :y ys}) :mapping (gg/aes :x :x :y :y))
             (gg/theme_light)
             (gg/geom_segment :mapping (gg/aes :x :x :xend :x :y 0 :yend :y) :color color))
       (add-common opts))))

(defn errorbars
  ([eb] (errorbars eb nil))
  ([eb opts]
   (let [data (tc/dataset (map #(zipmap [:name :xmin :xmax :x] (flatten %)) (reverse eb)))]
     (-> (r/r+ (gg/ggplot :data data :mapping (gg/aes :x :x :xmin :xmin :xmax :xmax :y :name))
               (gg/theme_light)
               (gg/geom_point :size 2)
               (gg/geom_errorbarh :color color-main :height 0.1))
         (add-common opts)))))

(defn functions
  "Functions"
  ([fs] (functions fs nil))
  ([fs {:keys [palette legend-name linetype?]
        :or {palette palette-blue-0 legend-name "Functions" linetype? true}
        :as opts}]
   (let [data (functions->data fs opts)
         breaks (distinct (map :fname data))]
     (-> (r/r+ (gg/ggplot)
               (gg/theme_light)
               (gg/geom_line :data (tc/dataset data) :mapping (if linetype?
                                                                (gg/aes :linetype :fname :color :fname :x :x :y :y)
                                                                (gg/aes :color :fname :x :x :y :y)))
               (when linetype? (gg/scale_linetype_manual :name legend-name :breaks breaks :values (map inc (range (count breaks)))))
               (gg/scale_color_manual :name legend-name :breaks breaks :values palette))
         (add-common opts)))))

(defn function2d
  ([f] (function2d f nil))
  ([f {:keys [palette legend-name]
       :or {palette :pals/ocean.ice}
       :as opts}]
   (let [data (function2d->data f opts)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y :fill :z))
               (gg/theme_light)
               (gg/geom_raster :data (tc/dataset data) :interpolate true)
               (pal/scale_fill_paletteer_c :name legend-name (->palette palette)))
         (add-common opts)))))

(defn complex-function
  ([f] (complex-function f nil))
  ([f opts]
   (let [data (complex-function->data f opts)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y))
               (gg/theme_light)
               (gg/guides :fill false)
               (gg/scale_fill_identity)
               (gg/coord_fixed)
               (gg/geom_raster :data (tc/dataset data) :interpolate true
                               :mapping (gg/aes :fill '(hsv hue sat val)))
               (gg/xlab "Re")
               (gg/ylab "Im"))
         (add-common opts)))))

(defn spectrogram
  ([stft] (spectrogram stft nil))
  ([stft {:keys [palette legend-name]
          :or {palette :grDevices/Plasma legend-name "Power (dB)"}
          :as opts}]
   (let [data (spectrogram->data stft)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y :fill :a))
               (gg/theme_light)
               (gg/geom_raster :data (tc/dataset data) :interpolate true)
               (gg/xlab "Time")
               (gg/ylab "Frequency")
               (pal/scale_fill_paletteer_c :name legend-name (->palette palette)))
         (add-common opts)))))


(defn function-ci
  ([f] (function-ci f nil))
  ([f {:keys [color alpha fill]
       :or {color color-main alpha 0.4 fill color-light}
       :as opts}]
   (let [data (function-ci->data f opts)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y :ymin :ymin :ymax :ymax))
               (gg/theme_light)
               (gg/geom_ribbon :data (tc/dataset data) :alpha alpha :fill fill)
               (gg/geom_line data :color color))
         (add-common opts)))))


(defn function+scatter
  ([f xs ys] (function+scatter f xs ys nil))
  ([f xs ys {:keys [x dot-size dot-alpha dot-shape]
             :or {dot-size 3.0 dot-alpha 0.8 dot-shape "circle filled"}
             :as opts}]
   (let [[x-min-p x-max-p] x
         [x-min-e x-max-e] (stats/extent xs)
         x-min (or x-min-p x-min-e)
         x-max (or x-max-p x-max-e)
         data (tc/dataset {:x xs :y ys})
         ff (if (sequential? f) functions function)]
     (r/r+ (ff f (assoc opts :x [x-min x-max]))
           (gg/geom_point :color "blue" :fill "light blue" :data data :mapping (gg/aes :x :x :y :y)
                          :shape dot-shape :size dot-size :alpha dot-alpha)))))

(defn scatter
  "plot scatter graph.
  `xs` and `ys` are the data points.
  
  options `opts`:
  `aspect-ratio`: default `nil`.
                  `nil` means letting the scale of the coordinate-system follows
                  drawing-window's aspect ratio.
                  set to any positive whole number to set fixed ratio y/x
                  for the coordinate-system.
                  setting aspect-ratio to true or 1 is useful if we want
                  to draw plots depicting any perfect circle.
  `color`: default `color-main`.
  `fill-color`: default `color-light`. "
  ([pairs] (scatter pairs nil))
  ([xs ys]
   (scatter xs ys nil))
  ([xs ys {:keys [aspect-ratio color fill-color]
           :or {color color-main fill-color color-light}
           :as opts}]
   (let [ys? ys
         ys (if ys? ys (map second xs))
         xs (if ys? xs (map first xs))]
     (-> (tc/dataset {:x xs :y ys})
         (gg/ggplot (gg/aes :x :x :y :y))
         (r/r+ (gg/theme_light)
               (gg/geom_point :color color :fill fill-color
                              :shape "circle filled" :size 1 :alpha 0.8))
         (cond-> aspect-ratio (r/r+ (gg/coord_fixed :ratio aspect-ratio)))
         (add-common opts)))))

(defn scatters
  "plot scatter graph.
  `xs` and `ys` are the data points.
  
  options `opts`:
  `aspect-ratio`: default `nil`.
                  `nil` means letting the scale of the coordinate-system follows
                  drawing-window's aspect ratio.
                  set to any positive whole number to set fixed ratio y/x
                  for the coordinate-system.
                  setting aspect-ratio to true or 1 is useful if we want
                  to draw plots depicting any perfect circle.
  `color`: default `color-main`.
  `fill-color`: default `color-light`. "
  ([xsysgroups])
  ([xsysgroups {:keys [aspect-ratio color size legend-name palette]
                :or {color color-main size 1 legend-name "Groups" palette palette-blue-2}
                :as opts}]
   (let [ds (mapcat (fn [[group xs ys]]
                      (map (fn [x y] {:x x :y y :group group}) xs ys)) xsysgroups)
         breaks (distinct (map first xsysgroups))]
     (-> (tc/dataset ds)
         (gg/ggplot (gg/aes :x :x :y :y :fill :group))
         (r/r+ (gg/theme_light)
               (gg/geom_point :shape "circle filled" :color color :size size :alpha 0.8)
               (gg/scale_fill_manual :name legend-name :breaks breaks :values palette))
         (cond-> aspect-ratio (r/r+ (gg/coord_fixed :ratio aspect-ratio)))
         (add-common opts)))))

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

(defn fgraph-int
  ([f domain]
   (fgraph-int f domain {}))
  ([f domain opts]
   (function f (assoc opts :x domain))))


(defn sample-int
  [f ^long dx ^long dy]
  (map (fn [v] [(+ v 0.5) (f v)]) (range dx (inc dy))))


(defn bgraph-int
  ([f] (bgraph-int f nil))
  ([f domain] (bgraph-int f domain {}))
  ([f domain opts]
   (let [[dx dy] domain
         xsys (if (map? f)
                (seq f)
                (sample-int f dx dy))
         xs (map first xsys)
         ys (map second xsys)]
     (lollipop xs ys (assoc opts :x domain)))))


(defn dgraph-cont
  "returns a vector of [pdf-plot cdf-plot icdf-plot] of a `distr`.
  Use supplied `pdf-graph` to plot PDF.
  Use `fgraph-int` to plot CDF and iCDF.
  params:
  `opts`:
     keys:`:pdf` a vector of pdf range, default 0 to 1
          `:icdf` a vector of icdf range, default 0 0.999
          `:data` supplied data"
  ([pdf-graph distr] (dgraph-cont pdf-graph distr nil))
  ([pdf-graph distr {:keys [pdf icdf data]
                     :or {pdf [0 1] icdf [0 0.999]}}]
   (let [pdf-plot (pdf-graph (if data data (partial fr/pdf distr)) pdf {:title "PDF"})
         cdf-plot (fgraph-int (partial fr/cdf distr) pdf {:title "CDF"})
         icdf-plot (fgraph-int (partial fr/icdf distr) icdf {:title "ICDF"})]
     [pdf-plot cdf-plot icdf-plot])))


(defn dgraph
  "given a distribution `distr`, return PDF, CDF, and ICDF plot images of that distribution.
  use `fgraph-int` to plot pdf.

  params:
  `opts`: see docs for `dgraph-cont`."
  ([distr]
   (dgraph distr {}))
  ([distr opts]
   (dgraph-cont fgraph-int distr opts)))


(defn dgraphi
  "given a distribution `distr`, return PDF, CDF, and ICDF plot images of that distribution.
  use `bgraph-int` to plot pdf.

  params:
  `opts`: see docs for `dgraph-cont`"
  ([distr]
   (dgraph distr {}))
  ([distr opts]
   (dgraph-cont bgraph-int distr opts)))

;;

(defn- acf-
  ([prefix series cnt nm]
   (let [acf? (= prefix "ACF")
         r (range cnt)
         a ((if acf? stats/acf-ci stats/pacf-ci) series)]
     (-> (lollipop r (if acf? (:acf a) (:pacf a)) 
                   {:title (str prefix " for " nm)
                    :xlab "lag"
                    :ylab (if acf?
                            "autocorrelation"
                            "partial autocorrelation")
                    :ylim [nil 1]
                    :hline {:yintercept (:ci a)}})
         (add-common {:hline {:yintercept (- (:ci a))}})))))

(defn acf
  ([series nm] (acf- "ACF" series 100 nm))
  ([series cnt nm] (acf- "ACF" series cnt nm)))

(defn pacf
  ([series nm] (acf- "PACF" series 100 nm))
  ([series cnt nm] (acf- "PACF" series cnt nm)))

;;

(defn half-split
  ([coeffs] (half-split coeffs []))
  ([coeffs buff]
   (if (m/< (count coeffs) 2)
     buff
     (let [[a b] (split-at (/ (count coeffs) 2) coeffs)]
       (recur a (conj buff b))))))

(defn ->level-data
  [max-size log-coeffs? normalize? id level]
  (let [y (inc id)
        dups (* 2 (/ max-size (count level)))
        sdups (m/sqrt dups)
        tr (if normalize? (comp (fn [^double z] (m// z sdups)) m/abs) m/abs)
        tr (if log-coeffs? (comp m/log1p tr) tr)]
    (->> (mapcat (partial repeat dups) level)
         (map-indexed (fn [x ^double z] {:x x :y y :z (tr z)})))))

(defn dwt-scaleogram
  ([coeffs] (dwt-scaleogram coeffs nil))
  ([coeffs {:keys [palette log-coeffs? normalize? wavelet]
            :or {palette :pals/parula log-coeffs? false normalize? false}}]
   (let [levels (half-split coeffs)
         max-size (count (first levels))
         ds (tc/dataset (mapcat (partial ->level-data max-size log-coeffs? normalize?) (range) levels))]
     (r/r+ (gg/ggplot ds (gg/aes :x :x :y :y :fill :z))
           (gg/theme_light)
           (gg/geom_raster)
           (gg/scale_y_reverse :breaks (range 1 (m/inc (m/log2 (count coeffs)))))
           (gg/ylab "Levels")
           (gg/xlab "Time")
           (pal/scale_fill_paletteer_c (->palette palette))
           (gg/labs :title (if wavelet
                             (str "DWT Scaleogram (" wavelet ")")
                             (str "DWT Scaleogram"))
                    :fill (if-not log-coeffs? "|coeffs|" "log(1+|coeffs|)"))))))

;;

(defn wavelet-pair
  ([wv] (wavelet-pair wv nil))
  ([wv {:keys [level reconstruction?]
        :or {level 7 reconstruction? false}
        :as options}]
   (let [w (if (or (string? wv) (keyword? wv)) (wv/wavelet wv) wv)
         nm (if (or (string? wv) (keyword? wv)) wv (wv/wavelet-name w))
         max-x (dec (wv/coeffs-size w))
         x [0 max-x]
         dr (if reconstruction? "(reconstruction)." "(deconstruction).")]
     [(function (wv/scaling-function w level reconstruction?) (assoc options :x x :title (str "Scaling function of " nm " wavelet " dr)))
      (function (wv/wavelet-function w level reconstruction?) (assoc options :x x :title (str "Wavelet function of " nm " wavelet " dr)))])))

;;

(defn magnitudes
  [xs]
  (-> (->> (t/forward-1d (t/transformer :real :fft) xs)
           (partition 2)
           (map v/mag))
      (v/div (count xs))
      (v/mult 2)))

(defn calc-dB
  [window {:keys [size pad cut]
           :or {size 256 pad 4096 cut 200}
           :as opts}]
  (let [w (kernel/window window size opts)
        coeffs (if (m/pos? pad) (signal/pad w pad :zero) w)
        s (map signal/linear->db (magnitudes coeffs))
        db (if (m/pos? cut) (take cut s) s)]
    (map #(m/constrain % -150.0 ##Inf) (v/shift db (m/- (double (first db)))))))

(defn fft-window-plot
  ([window] (fft-window-plot window {}))
  ([window opts]
   (let [db (calc-dB window opts)]
     (if (:symmetry? opts)
       (let [dbs (concat (reverse (rest db)) db)
             c (count db)]
         (line (range (- c) c) dbs opts))
       (line (range) db opts)))))

(defn fft-window-plots
  ([window sig] (fft-window-plots window sig nil))
  ([window sig opts]
   (let [t (str (name window) " window")
         title (if opts (str t " " opts) t)]
     [(line (range) (kernel/window window 256 opts) {:ylim [-0.1 nil]
                                                     :xlab "samples"
                                                     :title title})
      (line (range) (magnitudes (v/emult sig (kernel/window window 256 opts)))
            {:title "FFT of a tapered signal"
             :ylab "magnitude" :xlab "bins"})
      (fft-window-plot window (assoc opts :pad 400 :cut 0 :size 128
                                     :symmetry? true
                                     :ylab "dBc" :xlab "bins"
                                     :ylim [-150 nil]
                                     :title "FFT of a padded window"))
      (fft-window-plot window (assoc opts :title "FFT of padded window (main- and side-lobes)"
                                     :ylab "dBc" :xlab "bins"
                                     :ylim [-150 nil]))])))
