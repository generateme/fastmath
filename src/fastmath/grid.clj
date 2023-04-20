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
      * Anchors for triangular grids are shared between two cells: even and odd `q` coordinate. Even `q` is pointy topped, odd `q` is flat topped which is marked in the third coordinate of anchor."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot])
  (:import [fastmath.vector Vec2 Vec3]
           [clojure.lang Named]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn coords->cell
  "Converts 2d space coordinates (x,y) to cell coordinates (q,r)."
  ([g coords] (prot/coords->cell g coords))
  ([g x y] (prot/coords->cell g x y)))

(defn cell->anchor
  "Converts cell coordinates (q,r) to anchor coordinates (x,y)."
  ([g cell] (prot/cell->anchor g cell))
  ([g q r] (prot/cell->anchor g q r)))

(defn coords->mid
  "Converts 2d space coordinates (x,y) into cell midpoint (x,y)."
  ([g coords] (prot/coords->mid g coords))
  ([g x y] (prot/coords->mid g x y)))

(defn grid-type
  "Returns type of the cell."
  [g] (prot/grid-type g))

(defn corners
  "Returns list of cell vertices for given 2d space coordinates."
  ([g coords] (prot/corners g coords))
  ([g coords scale] (prot/corners g coords scale))
  ([g x y scale] (prot/corners g x y scale)))

(defn coords->anchor
  "Converts 2d coordinates (x,y) to cell anchor (x,y)."
  ([g coords] (prot/cell->anchor g (prot/coords->cell g coords)))
  ([g x y] (prot/cell->anchor g (prot/coords->cell g x y))))

(defn cell->mid
  "Converts cell coordinates (q,r) to mid point (x,y)."
  ([g cell] (prot/cell->mid g cell))
  ([g x y] (prot/cell->mid g x y)))

(defn- grid-obj
  "Create grid object."
  ([typ {:keys [from-cell to-cell to-mid vertices anchor]} ^double size ^Vec2 sv]
   (let [nm (name typ)
         tostr (str nm ", size=" size)
         anchor (or anchor coords->anchor)
         ^Vec3 sv3 (v/vec3 sv 0.0)
         triangle? (= typ :triangle)]
     (reify
       prot/GridProto
       (coords->cell [_ [^double x ^double y]] (to-cell size (- x (.x sv)) (- y (.y sv))))
       (coords->cell [_ x y] (to-cell size (- ^double x (.x sv)) (- ^double y (.y sv))))
       (cell->anchor [_ [q r]] (let [cell (from-cell size q r)]
                                 (if triangle?
                                   (v/add cell sv3)
                                   (v/add cell sv))))
       (cell->anchor [_ q r] (let [cell (from-cell size q r)]
                               (if triangle?
                                 (v/add cell sv3)
                                 (v/add cell sv))))
       (cell->mid [g cell] (to-mid size (prot/cell->anchor g cell)))
       (cell->mid [g q r] (to-mid size (prot/cell->anchor g q r)))
       (coords->mid [g coords] (to-mid size (anchor g coords)))
       (coords->mid [g x y] (to-mid size (anchor g x y)))
       (corners [g coords] (vertices size (anchor g coords)))
       (corners [g coords scale] (vertices (* ^double scale size) (anchor g coords)))
       (corners [g x y scale] (vertices (* ^double scale size) (anchor g x y)))
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

(defn- square-pixel->mid
  "Mid point of cell"
  [^double size v]
  (let [s2 (/ size 2.0)]
    (v/add v (Vec2. s2 s2))))

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
    (if (and (> diff-x diff-y)
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

(defn- hex->mid [_ v] v)

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

(defn- rhombus->pixel
  "Shifted square to anchor."
  [^double size ^long q ^long r]
  (v/vec2 (+ (* r (* 0.5 size)) (* q size)) (* r size m/SQRT3_2)))

(defn- pixel->rhombus
  "2d coords to rhombus cell."
  [^double size ^double x ^double y]
  (let [h (* size m/SQRT3_2)
        ys (/ y h)
        yy (m/floor ys)
        fy (if (neg? ys) (- 1.0 (m/frac ys)) (m/frac ys))
        hs (* 0.5 size)]
    (v/vec2 (m/floor (/ (+ (* fy hs)
                           (- x (* yy hs))) size))
            yy)))

(defn- rhombus-pixel->mid
  [^double size v]
  (let [hs (* m/SQRT3_4 size)
        qs (* 0.25 size)]
    (v/add v (Vec2. qs hs))))

(defn- rhombus-corners
  "Rhombus vertices."
  ([^double size ^double x ^double y]
   (let [hs (* 0.5 size)
         y+ (+ y (* m/SQRT3_2 size))]
     [(v/vec2 x y)
      (v/vec2 (+ x size) y)
      (v/vec2 (+ x hs) y+)
      (v/vec2 (- x hs) y+)]))
  ([size [x y]] (rhombus-corners size x y)))


;; triangle

(defn- triangle->pixel
  "Triangle cell to anchor."
  [^double size ^long q ^long r]
  (v/vec3 (rhombus->pixel size (m/>> q 1) r) (bit-and q 0x1)))

(defn- pixel->triangle
  "2d coords to triangle cell."
  [^double size ^double x ^double y]
  (let [h (* size m/SQRT3_2)
        ys (/ y h)
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
   (let [y+ (+ y (* m/SQRT3_2 size))
         hsize (* 0.5 size)]
     (if (zero? down?)
       [(v/vec2 x y)
        (v/vec2 (+ x hsize) y+)
        (v/vec2 (- x hsize) y+)]
       [(v/vec2 x y)
        (v/vec2 (+ x size) y)
        (v/vec2 (+ x hsize) y+)])))
  ([size [x y down?]] (triangle-corners size x y down?)))

(defn- triangle-pixel->mid
  [^double size [^double x ^double y ^long down?]]
  (let [h (* m/SQRT3_2 size)]
    (if (zero? down?)
      (v/vec2 x (+ y (* m/TWO_THIRD h)))
      (v/vec2 (+ x (/ size 2.0)) (+ y (* m/THIRD h))))))

(defn- coords->triangle-anchor
  "2d coordinates to triangle anchor with cell information (even or odd)."
  ([g coords]
   (let [^Vec2 qr (prot/coords->cell g coords)]
     (prot/cell->anchor g qr)))
  ([g x y]
   (let [^Vec2 qr (prot/coords->cell g x y)]
     (prot/cell->anchor g qr))))

;;

(def ^:private grid-type-fns
  {:square {:from-cell square->pixel :to-cell pixel->square :to-mid square-pixel->mid :vertices square-corners}
   :pointy-hex {:from-cell pointy-hex->pixel :to-cell pixel->pointy-hex :to-mid hex->mid :vertices pointy-hex-corners}
   :flat-hex {:from-cell flat-hex->pixel :to-cell pixel->flat-hex :to-mid hex->mid :vertices flat-hex-corners}
   :shifted-square {:from-cell shifted-square->pixel :to-cell pixel->shifted-square :to-mid square-pixel->mid :vertices square-corners}
   :rhombus {:from-cell rhombus->pixel :to-cell pixel->rhombus :to-mid rhombus-pixel->mid :vertices rhombus-corners}
   :triangle {:from-cell triangle->pixel :to-cell pixel->triangle :to-mid triangle-pixel->mid :vertices triangle-corners :anchor coords->triangle-anchor}})

(def ^:private grid-type-fns-default {:from-cell square->pixel :to-cell pixel->square :to-mid square-pixel->mid :vertices square-corners})

(defn grid
  "Create grid for given type, size and optional translating vector."
  ([type ^double size ^double sx ^double sy]
   (let [sv (Vec2. sx sy)
         size (if (or (= type :pointy-hex)
                      (= type :flat-hex)) (/ size 2.0) size)
         fp (get grid-type-fns type grid-type-fns-default)]
     (grid-obj type fp size sv)))
  ([type ^double size] (grid type size 0.0 0.0))
  ([type] (grid type 10.0))
  ([] (grid :square)))

(def cell-names ^{:doc "List of cell types"} [:square :shifted-square :triangle :rhombus :flat-hex :pointy-hex])
