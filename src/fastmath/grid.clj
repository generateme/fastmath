(ns fastmath.grid
  "Grid calculation functions.

  Convert 2d coordinates into various grid coordinates and back.

  Terms used:

  * cell type - grid cell shape: square, triangular, hexagonal, rhomboidal
  * coords - 2d euclidean coordinates (x,y)
  * cell - cell coordinates (q,r)
  * anchor - cell position in 2d euclidean space
  * corners - shape vertices
  * size - size of the cell.

  ### Grids

  Each grid is defined by cell type and size. Optionally you can provide translating vector.

  Each cell has it's own coordinates, mostly axial based (only square has offset).

  For hexagonal cell size is a radius from midpoint to corner. For the rest it is the size of the side.

  Cell types are:

  * `:square`
  * `:shifted-square`
  * `:triangle`
  * `:rhombus`
  * `:flat-hex` - flat topped
  * `:pointy-hex` - pointy topped

  ### Notes

  * Hexagonal grids are based on https://www.redblobgames.com/grids/hexagons/
  * Only hexagonal cells have anchor at the center. For the rest the anchor is at the top left vertex.
      * Anchors for triangular grids are shared between two cells: even and odd `q` coordinate. Even `q` is pointy topped, odd `q` is flat topped."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec3]
           [clojure.lang Named]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defprotocol GridProto
  "Common grid conversion functions."
  (coords->cell [g coords] "Converts 2d space coordinates to cell coordinates.")
  (cell->anchor [g cell] "Converts cell coordinates to anchor coordinates.")
  (grid-type [g] "Returns type of the cell.")
  (corners [g coords] [g coords scale] "Returns list of cell vertices for given 2d space coordinates."))

(defn coords->anchor
  "Converts 2d coordinates to cell anchor."
  [g coords]
  (cell->anchor g (coords->cell g coords)))

(defn- grid-obj
  "Create grid object."
  ([typ {:keys [from-cell to-cell vertices anchor]} ^double size ^Vec2 sv]
   (let [nm (subs (str typ) 1)
         tostr (str nm ", size=" size)
         anchor (or anchor coords->anchor)]
     (reify
       GridProto
       (coords->cell [_ [^double x ^double y]] (to-cell size (- x (.x sv)) (- y (.y sv))))
       (cell->anchor [_ [^double q ^double r]] (v/add sv (from-cell size q r)))
       (corners [g coords] (vertices size (anchor g coords)))
       (corners [g coords scale] (vertices (* ^double scale size) (anchor g coords)))
       (grid-type [_] typ)
       Named
       (getName [_] nm)
       Object
       (toString [_] tostr)))))

;; square

(defn- square->pixel
  "Square cell coords to anchor."
  [^double size ^long q ^long r]
  (Vec2. (* q size) (* r size)))

(defn- pixel->square
  "Square 2d coords to cell coords."
  [^double size ^double x ^double y]
  (Vec2. (m/floor (/ x size)) (m/floor (/ y size))))

(defn- square-corners
  "Square shape."
  ([^double size ^double x ^double y]
   (let [x+ (+ x size)
         y+ (+ y size)]
     [(v/vec2 x y)
      (v/vec2 x+ y)
      (v/vec2 x+ y+)
      (v/vec2 x y+)]))
  ([size [x y]] (square-corners size x y)))

;; hex

(defn- hex-round-coords
  "Round hex cell coordinates."
  ^Vec2 [^Vec2 v]
  (let [x (.x v)
        z (.y v)
        y (- (- x) z)
        rx (m/round x)
        ry (m/round y)
        rz (m/round z)
        diff-x (m/abs (- rx x))
        diff-y (m/abs (- ry y))
        diff-z (m/abs (- rz z))]
    (if (bool-and (> diff-x diff-y)
                  (> diff-x diff-z))
      (Vec2. (- (- ry) rz) rz)
      (if (> diff-y diff-z)
        (Vec2. rx rz)
        (Vec2. rx (- (- rx) ry))))))

(defn- pointy-hex->pixel
  "Pointy hex cell to anchor."
  [^double size ^double q ^double r]
  (v/mult (Vec2. (+ (* m/SQRT3 q) (* m/SQRT3_2 r))
                 (* 1.5 r)) size))

(defn- pixel->pointy-hex
  "2d coords to pointy topped hex cell coords."
  [^double size ^double x ^double y]
  (hex-round-coords (v/div (v/vec2 (- (* m/SQRT3_3 x) (* m/THIRD y))
                                   (* m/TWO_THIRD y)) size)))

(defn- flat-hex->pixel
  "Flat hex cell to anchor."
  [^double size ^long q ^long r]
  (v/mult (v/vec2 (* 1.5 q)
                  (+ (* m/SQRT3_2 q) (* m/SQRT3 r))) size))

(defn- pixel->flat-hex
  "2d coords to flat topped hex cell coords."
  [^double size ^double x ^double y]
  (hex-round-coords (v/div (v/vec2 (* m/TWO_THIRD x)
                                   (- (* m/SQRT3_3 y) (* m/THIRD x))) size)))

(defn ^:private angles->vectors
  "Returns list of vertices from unit circle."
  [^double angle]
  (v/from-polar (Vec2. 1.0 (m/radians angle))))

(def ^:private flat-topped-sc (map angles->vectors (range 0 360 60)))
(def ^:private pointy-topped-sc (map angles->vectors (range 30 360 60)))

(defn- hex-corners
  "Returns hex corner vertices."
  ([lst ^double size ^double x ^double y]
   (let [in (v/vec2 x y)]
     (map #(v/add in (v/mult % size)) lst)))
  ([lst ^double size center]
   (map #(v/add center (v/mult % size)) lst)))

(def ^{:doc "Function which returns vertices for pointy topped hexagon for given size and coordinates."}
  pointy-hex-corners (partial hex-corners pointy-topped-sc))
(def ^{:doc "Function which returns vertices for flat topped hexagon for given size and coordinates."}
  flat-hex-corners (partial hex-corners flat-topped-sc))

;; shifted square

(defn- shifted-square->pixel
  "Shifted square to anchor."
  [^double size ^long q ^long r]
  (v/vec2 (+ (* r (* 0.5 size)) (* q size)) (* r size)))

(defn- pixel->shifted-square
  "2d coords to shifted square cell."
  ([^double size ^double x ^double y]
   (let [yy (m/floor (/ y size))]
     (v/vec2 (m/floor (/ (- x (* yy (* 0.5 size))) size)) yy))))

;; rhombus

(def ^:private rhombus->pixel shifted-square->pixel)

(defn- pixel->rhombus
  "2d coords to rhombus cell."
  [^double size ^double x ^double y]
  (let [ys (/ y size)
        yy (m/floor ys)
        fy (if (neg? ys) (- 1.0 (m/frac ys)) (m/frac ys))
        hs (* 0.5 size)]
    (v/vec2 (m/floor (/ (+ (* fy hs)
                           (- x (* yy hs))) size))
            yy)))

(defn- rhombus-corners
  "Rhombus vertices."
  ([^double size ^double x ^double y]
   (let [hs (* 0.5 size)
         y+ (+ y size)]
     [(v/vec2 x y)
      (v/vec2 (+ x size) y)
      (v/vec2 (+ x hs) y+)
      (v/vec2 (- x hs) y+)]))
  ([size [x y]] (rhombus-corners size x y)))


;; triangle

(defn- triangle->pixel
  "Triangle cell to anchor."
  [^double size ^long q ^long r]
  (shifted-square->pixel size (>> q 1) r))

(defn- pixel->triangle
  "2d coords to triangle cell."
  [^double size ^double x ^double y]
  (let [ys (/ y size)
        yy (m/floor ys)
        fy (if (neg? ys) (- 1.0 (m/frac ys)) (m/frac ys))
        hs (* 0.5 size)
        xs (/ (+ (* fy hs)
                 (- x (* yy hs))) size)
        xx (* 2.0 (m/floor xs))
        fx (if (neg? xs) (- 1.0 (m/frac xs)) (m/frac xs))]
    (if (< fy fx)
      (Vec2. (inc xx) yy)
      (Vec2. xx yy))))

(defn- triangle-corners
  "Triangle shape."
  ([^double size ^double x ^double y ^long down?]
   (if (zero? down?)
     [(v/vec2 x y)
      (v/vec2 (+ x (/ size 2.0)) (+ y size))
      (v/vec2 (- x (/ size 2.0)) (+ y size))]
     [(v/vec2 x y)
      (v/vec2 (+ x size) y)
      (v/vec2 (+ x (/ size 2.0)) (+ y size))]))
  ([size [x y down?]] (triangle-corners size x y down?)))

(defn- coords->triangle-anchor
  "2d coordinates to triangle anchor with cell information (even or odd)."
  [g coords]
  (let [^Vec2 qr (coords->cell g coords)]
    (v/vec3 (cell->anchor g qr) (bit-and (long (.x qr)) 0x1))))

;;

(defn- grid-type-fns
  "Return functions for grid calculations."
  [type]
  (case type
    :square {:from-cell square->pixel :to-cell pixel->square :vertices square-corners}
    :pointy-hex {:from-cell pointy-hex->pixel :to-cell pixel->pointy-hex :vertices pointy-hex-corners}
    :flat-hex {:from-cell flat-hex->pixel :to-cell pixel->flat-hex :vertices flat-hex-corners}
    :shifted-square {:from-cell shifted-square->pixel :to-cell pixel->shifted-square :vertices square-corners}
    :rhombus {:from-cell rhombus->pixel :to-cell pixel->rhombus :vertices rhombus-corners}
    :triangle {:from-cell triangle->pixel :to-cell pixel->triangle :vertices triangle-corners :anchor coords->triangle-anchor}
    {:from-cell square->pixel :to-cell pixel->square :vertices square-corners}))

(defn grid
  "Create grid for given type, size and optional translating vector."
  ([type ^double size ^double sx ^double sy]
   (let [sv (Vec2. sx sy)
         fp (grid-type-fns type)]
     (grid-obj type fp size sv)))
  ([type ^double size] (grid type size 0.0 0.0))
  ([type] (grid type 10.0))
  ([] (grid :square)))

(def cell-names ^{:doc "List of cell types"} [:square :shifted-square :triangle :rhombus :flat-hex :pointy-hex])


