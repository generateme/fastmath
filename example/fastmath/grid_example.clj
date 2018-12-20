(ns fastmath.grid-example
  (:require [fastmath.grid :refer :all]
            [metadoc.examples :refer :all]))

(add-examples grid
  (example-session "Usage"
    (grid :pointy-hex 20)
    (grid :triangular 10 5.0 5.0)
    (name (grid :square))
    (str (grid :rhombus))))

(add-examples coords->cell
  (example "Usage"
    (coords->cell (grid :pointy-hex 20 10 10) [50 50]))
  (let [g (grid :pointy-hex 20 10 10)]
    (example-session "More examples (grid `g` same as above)"
      (coords->cell g [0 0])
      (coords->cell g [100 100])
      (coords->cell g [-150 -150]))))

(add-examples coords->cell
  (example "Usage"
    (cell->anchor (grid :triangle 20 10 10) [5 5]))
  (let [g (grid :triangle 20 10 10)]
    (example-session "More examples (grid `g` same as above). For triangle grid pair of cells share the same anchor."
      (cell->anchor g [0 0])
      (cell->anchor g [1 0])
      (cell->anchor g [2 0])
      (cell->anchor g [3 0])
      (cell->anchor (grid :shifted-square 20 10 10) [5 5]))))

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

(add-examples pointy-hex-corners
  (example-session "Generate vertices for flat topped hexagon for given position and size"
    (flat-hex-corners 100 [15 15])
    (flat-hex-corners 1 [0 0])
    (flat-hex-corners 100 0 0)))
