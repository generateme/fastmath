(ns fastmath.grid-examples
  (:require [fastmath.grid :refer :all]
            [metadoc.examples :refer :all]))

(add-examples grid
  (example-session "Usage"
    (grid :pointy-hex 20)
    (grid :triangular 10 5.0 5.0)
    (name (grid :square))
    (str (grid :rhombus))
    (grid))
  (example-image "Pointy hex (red dot: mid point, cyan circle: anchor)" "images/g/pointy-hex.jpg")
  (example-image "Flat hex (red dot: mid point, cyan circle: anchor)" "images/g/flat-hex.jpg")
  (example-image "Rhombus (red dot: mid point, cyan circle: anchor)" "images/g/rhombus.jpg")
  (example-image "Triangle (red dot: mid point, cyan circle: anchor)" "images/g/triangle.jpg")
  (example-image "Square (red dot: mid point, cyan circle: anchor)" "images/g/square.jpg")
  (example-image "Shifted square (red dot: mid point, cyan circle: anchor)" "images/g/shifted-square.jpg"))

(add-examples coords->cell
  (example "Usage"
    (coords->cell (grid :pointy-hex 20 10 10) [50 50]))
  (let [g (grid :pointy-hex 20 10 10)]
    (example-session "More examples (grid `g` same as above)"
      (coords->cell g [0 0])
      (coords->cell g [100 100])
      (coords->cell g [-150 -150]))))

(add-examples coords->mid
  (example "Usage"
    (coords->mid (grid :pointy-hex 20 10 10) [50 50]))
  (let [g (grid :pointy-hex 20 10 10)]
    (example-session "More examples (grid `g` same as above)"
      (coords->mid g [0 0])
      (coords->mid g [100 100])
      (coords->mid g [-150 -150]))))

(add-examples cell->anchor
  (example "Usage"
    (cell->anchor (grid :triangle 20 10 10) [5 5]))
  (let [g (grid :triangle 20 10 10)]
    (example-session "More examples (grid `g` same as above). For triangle grid pair of cells share the same anchor, third coordinate indicates if it's up (`0`) or down (`1`) triangle."
      (cell->anchor g [0 0])
      (cell->anchor g [1 0])
      (cell->anchor g [2 0])
      (cell->anchor g [3 0])
      (cell->anchor (grid :shifted-square 20 10 10) [5 5]))))

(add-examples cell->mid
  (example "Usage"
    (cell->mid (grid :triangle 20 10 10) [4 5]))
  (let [g (grid :triangle 20 10 10)]
    (example-session "More examples (grid `g` same as above)."
      (cell->mid g [0 0])
      (cell->mid g [1 0])
      (cell->mid g [2 0])
      (cell->mid g [3 0])
      (cell->mid (grid :shifted-square 20 10 10) [5 5]))))

(add-examples grid-type
  (example (grid-type (grid :flat-hex 10))))

(add-examples corners
  (example-session "Usage"
    (corners (grid :flat-hex 20) [100 100])
    (corners (grid :triangle 10 10 10) [100 100])))

(add-examples coords->anchor
  (example-session "Usage"
    (coords->anchor (grid :rhombus) [105 105])
    (coords->anchor (grid :rhombus) [100 100])
    (coords->anchor (grid :rhombus) [99 99])))

(add-examples pointy-hex-corners
  (example-session "Generate vertices for pointy topped hexagon for given position and size"
    (pointy-hex-corners 100 [15 15])
    (pointy-hex-corners 1 [0 0])
    (pointy-hex-corners 100 0 0)))

(add-examples flat-hex-corners
  (example-session "Generate vertices for flat topped hexagon for given position and size"
    (flat-hex-corners 100 [15 15])
    (flat-hex-corners 1 [0 0])
    (flat-hex-corners 100 0 0)))

(add-examples cell-names
  (example "List of all grid types." cell-names))
