(ns fastmath.grid-test
  (:require [fastmath.grid :as sut]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t]))

;; shared test helpers (used across Groups 2-5)

(defn- point-in-polygon?
  "Ray-casting point-in-polygon test. `pts` is a seq of `[x y]` pairs."
  [[px py] pts]
  (let [pts (vec pts)
        n (count pts)]
    (loop [i 0 j (dec n) inside false]
      (if (= i n)
        inside
        (let [[xi yi] (nth pts i)
              [xj yj] (nth pts j)
              intersect (and (not= (> yi py) (> yj py))
                             (< px (+ xi (/ (* (- xj xi) (- py yi)) (- yj yi)))))]
          (recur (inc i) i (if intersect (not inside) inside)))))))

(defn- xy [c] [(nth c 0) (nth c 1)])

(defn- centroid [pts]
  (let [n (count pts)]
    [(/ (reduce + (map first pts)) n) (/ (reduce + (map second pts)) n)]))

(def ^:private sample-coords
  "Sample coordinates spanning negative/positive/zero in both axes, deliberately
  off any obvious grid alignment (irrational-ish steps) to avoid boundary bias."
  (vec (for [x (range -53.3 47 5.1) y (range -49 53.7 4.9)] [x y])))

(defn- roundtrip-ok?
  "`coords->cell` of `coords->mid`'s own output must return to the same cell."
  [gr [x y :as _coords]]
  (let [cell (sut/coords->cell gr x y)
        mid (sut/coords->mid gr x y)
        cell2 (sut/coords->cell gr (nth mid 0) (nth mid 1))]
    (= (mapv double (xy cell)) (mapv double (xy cell2)))))

(defn- contains-point?
  "The cell polygon returned by `corners` for `coords` must contain `coords`."
  [gr [x y :as coords]]
  (point-in-polygon? [x y] (map xy (sut/corners gr coords))))

;; Group 1 - Constructor & type dispatch

(t/deftest cell-names-test
  (t/is (= [:square :shifted-square :triangle :rhombus :flat-hex :pointy-hex] sut/cell-names))
  (t/is (= 6 (count sut/cell-names))))

(t/deftest grid-valid-types-construct-test
  (doseq [type sut/cell-names]
    (t/testing (str type)
      (t/is (some? (sut/grid type)))
      (t/is (some? (sut/grid type 20.0)))
      (t/is (some? (sut/grid type 20.0 5.0 7.0)))
      (t/is (= type (sut/grid-type (sut/grid type)))))))

(t/deftest grid-invalid-type-throws-test
  (doseq [call [#(sut/grid :bogus)
                #(sut/grid :bogus 20.0)
                #(sut/grid :bogus 20.0 5.0 7.0)]]
    (t/is (thrown? clojure.lang.ExceptionInfo (call)))))

(t/deftest grid-defaults-test
  (t/is (= :square (sut/grid-type (sut/grid))))
  (t/is (= "square, size=10.0" (str (sut/grid))))
  (t/is (= "square, size=10.0" (str (sut/grid :square))))
  (t/is (= (str (sut/grid :square 10.0 0.0 0.0)) (str (sut/grid :square 10.0))))
  (t/is (= (sut/coords->cell (sut/grid :square 10.0 0.0 0.0) [5 5])
           (sut/coords->cell (sut/grid :square 10.0) [5 5]))))

(t/deftest grid-hex-size-scaled-test
  ;; size passed to `grid` becomes the internal size (= circumradius) divided by sqrt(3) for
  ;; hex types, so that neighboring hex anchors end up `size` apart, same as every other type
  (t/is (= "flat-hex, size=5.773502691896258" (str (sut/grid :flat-hex 10.0))))
  (t/is (= "pointy-hex, size=5.773502691896258" (str (sut/grid :pointy-hex 10.0))))
  ;; non-hex types are unaffected
  (t/is (= "square, size=10.0" (str (sut/grid :square 10.0))))
  ;; the actual invariant this scaling exists for: anchor-to-anchor distance == size, for hex too
  (doseq [type [:flat-hex :pointy-hex]]
    (let [gr (sut/grid type 10.0)]
      (t/is (m/delta-eq 10.0 (v/dist (sut/cell->anchor gr 0 0) (sut/cell->anchor gr 1 0))))
      (t/is (m/delta-eq 10.0 (v/dist (sut/cell->anchor gr 0 0) (sut/cell->anchor gr 0 1)))))))

;; Group 2 - Square + shifted-square family

(t/deftest square-family-roundtrip-and-containment-test
  (doseq [type [:square :shifted-square]
          [sx sy] [[0.0 0.0] [7.0 -3.0]]]
    (t/testing (str type " sv=[" sx " " sy "]")
      (let [gr (sut/grid type 10.0 sx sy)]
        (doseq [c sample-coords]
          (t/is (roundtrip-ok? gr c))
          (t/is (contains-point? gr c)))))))

(t/deftest square-family-shape-test
  (doseq [type [:square :shifted-square]]
    (t/testing (str type)
      (let [gr (sut/grid type 10.0)
            pts (mapv xy (sut/corners gr [5.0 5.0]))]
        (t/is (= 4 (count pts)))
        (t/is (every? #(m/delta-eq (apply v/dist %) 10.0)
                       (map vector pts (concat (rest pts) [(first pts)]))))))))

(t/deftest square-family-mid-test
  (doseq [type [:square :shifted-square]]
    (t/testing (str type)
      (let [gr (sut/grid type 10.0)
            cell (sut/coords->cell gr 5.0 5.0)
            mid-from-coords (sut/coords->mid gr 5.0 5.0)
            mid-from-cell (sut/cell->mid gr cell)]
        (t/is (v/delta-eq mid-from-coords mid-from-cell))
        (t/is (contains-point? gr (xy mid-from-coords)))
        ;; mid == centroid of the cell's own corners
        (let [c (centroid (mapv xy (sut/corners gr [5.0 5.0])))]
          (t/is (m/delta-eq (nth mid-from-coords 0) (nth c 0)))
          (t/is (m/delta-eq (nth mid-from-coords 1) (nth c 1))))))))

(t/deftest square-family-corners-scale-test
  (doseq [type [:square :shifted-square]]
    (t/testing (str type)
      (let [gr (sut/grid type 10.0)
            full (mapv xy (sut/corners gr [5.0 5.0]))
            half (mapv xy (sut/corners gr 5.0 5.0 0.5))
            anchor (xy (sut/cell->anchor gr (sut/coords->cell gr 5.0 5.0)))]
        ;; scale halves the distance from anchor to each vertex
        (dotimes [i 4]
          (t/is (m/delta-eq (v/dist anchor (nth half i))
                             (* 0.5 (v/dist anchor (nth full i))))))))))

;; Group 3 - Hex family

(t/deftest hex-family-roundtrip-and-containment-test
  (doseq [type [:pointy-hex :flat-hex]
          [sx sy] [[0.0 0.0] [7.0 -3.0]]]
    (t/testing (str type " sv=[" sx " " sy "]")
      (let [gr (sut/grid type 10.0 sx sy)]
        (doseq [c sample-coords]
          (t/is (roundtrip-ok? gr c))
          (t/is (contains-point? gr c)))))))

(t/deftest hex-family-shape-test
  (doseq [type [:pointy-hex :flat-hex]]
    (t/testing (str type)
      ;; user size=10.0 -> internal size (= circumradius) = 10.0/sqrt(3) (anchor-to-anchor unification)
      (let [gr (sut/grid type 10.0)
            radius (/ 10.0 m/SQRT3)
            cell (sut/coords->cell gr 12.0 8.0)
            anchor (sut/cell->anchor gr cell)
            pts (mapv xy (sut/corners gr (xy anchor)))]
        (t/is (= 6 (count pts)))
        (doseq [p pts]
          (t/is (m/delta-eq (v/dist (xy anchor) p) radius)))
        ;; hex anchor is already the center: coords->mid/cell->mid equal cell->anchor
        (t/is (v/delta-eq anchor (sut/coords->mid gr (xy anchor))))
        (t/is (v/delta-eq anchor (sut/cell->mid gr cell)))))))

(t/deftest hex-family-closed-form-test
  ;; hand-derived from redblobgames axial->pixel formulas (namespace docstring reference),
  ;; internal size (circumradius) = 10.0/sqrt(3) for user size=10.0 (anchor-to-anchor unification)
  (t/testing :pointy-hex
    (let [gr (sut/grid :pointy-hex 10.0)
          anchor (sut/cell->anchor gr 1 1)]
      (t/is (m/delta-eq (nth anchor 0) 15.0))
      (t/is (m/delta-eq (nth anchor 1) (/ 15.0 m/SQRT3)))))
  (t/testing :flat-hex
    (let [gr (sut/grid :flat-hex 10.0)
          anchor (sut/cell->anchor gr 1 1)]
      (t/is (m/delta-eq (nth anchor 0) (/ 15.0 m/SQRT3)))
      (t/is (m/delta-eq (nth anchor 1) 15.0)))))

;; Group 4 - Rhombus family

(t/deftest rhombus-roundtrip-and-containment-test
  (doseq [[sx sy] [[0.0 0.0] [7.0 -3.0]]]
    (t/testing (str "sv=[" sx " " sy "]")
      (let [gr (sut/grid :rhombus 10.0 sx sy)]
        (doseq [c sample-coords]
          (t/is (roundtrip-ok? gr c))
          (t/is (contains-point? gr c)))))))

(t/deftest rhombus-shape-test
  (let [gr (sut/grid :rhombus 10.0)
        pts (mapv xy (sut/corners gr [5.0 5.0]))]
    (t/is (= 4 (count pts)))
    ;; all 4 sides equal `size` (a true rhombus, not just any parallelogram)
    (t/is (every? #(m/delta-eq (apply v/dist %) 10.0)
                   (map vector pts (concat (rest pts) [(first pts)]))))))

(t/deftest rhombus-mid-test
  (let [gr (sut/grid :rhombus 10.0)
        cell (sut/coords->cell gr 5.0 5.0)
        mid-from-coords (sut/coords->mid gr 5.0 5.0)
        mid-from-cell (sut/cell->mid gr cell)]
    (t/is (v/delta-eq mid-from-coords mid-from-cell))
    (t/is (contains-point? gr (xy mid-from-coords)))
    (let [c (centroid (mapv xy (sut/corners gr [5.0 5.0])))]
      (t/is (m/delta-eq (nth mid-from-coords 0) (nth c 0)))
      (t/is (m/delta-eq (nth mid-from-coords 1) (nth c 1))))))

;; Group 5 - Triangle family
;;
;; Methodology pitfalls (Research, see Group 5 topic note): `corners` for `:triangle`
;; routes through a coords->cell recompute internally (`coords->triangle-anchor`), so
;; using it to test `coords->cell` on a shared-anchor boundary point is circular and
;; over/under-reports failures unrelated to `pixel->triangle` itself. `roundtrip-ok?`/
;; `contains-point?` avoid this because `sample-coords` are interior points, never
;; exact grid vertices.

(t/deftest triangle-roundtrip-and-containment-test
  (doseq [[sx sy] [[0.0 0.0] [7.0 -3.0]]]
    (t/testing (str "sv=[" sx " " sy "]")
      (let [gr (sut/grid :triangle 10.0 sx sy)]
        (doseq [c sample-coords]
          (t/is (roundtrip-ok? gr c))
          (t/is (contains-point? gr c)))))))

(t/deftest triangle-interior-centroid-correctness-test
  ;; For every (q,r) in a wide range spanning both signs, the centroid of the cell's
  ;; own corners (an unambiguous interior point, never a shared vertex) must map back
  ;; to that same (q,r) via `coords->cell`. This is the fix's regression guard: it was
  ;; already true before the fix (0/289 in Research) and must remain true after it.
  ;;
  ;; Deliberately bypasses the public `corners` fn: for `:triangle` it routes through
  ;; `coords->cell` internally (`coords->triangle-anchor`) to resolve which cell's
  ;; corners to return, which is exactly the boundary-ambiguity this test must not be
  ;; confounded by (see the pitfall note above `triangle-roundtrip-and-containment-test`).
  ;; `triangle->pixel`/`triangle-corners` (private) build a specific cell's corners
  ;; directly from its own `(q,r)`, with no circular dependency on `coords->cell`.
  (let [gr (sut/grid :triangle 10.0)
        triangle->pixel @#'sut/triangle->pixel
        triangle-corners @#'sut/triangle-corners]
    (doseq [q (range -8 9) r (range -8 9)]
      (let [anchor (triangle->pixel 10.0 (long q) (long r))
            down? (long (nth anchor 2))
            pts (triangle-corners 10.0 (nth anchor 0) (nth anchor 1) down?)
            cen (centroid (mapv xy pts))
            cell-back (sut/coords->cell gr (nth cen 0) (nth cen 1))]
        (t/is (= (double q) (nth cell-back 0)))
        (t/is (= (double r) (nth cell-back 1)))))))

(t/deftest triangle-boundary-consistency-test
  ;; F3 fix: the shared anchor of an up/down pair (even `q` and odd `q+1`) cannot
  ;; round-trip to BOTH members (2 cells share 1 anchor point in 2d space; `coords->cell`
  ;; never sees the `down?` flag) - that ambiguity is structural, not a bug. What the fix
  ;; guarantees is that resolution is now deterministic and sign-symmetric: every shared
  ;; anchor consistently resolves to the even (`down?` false) member, for every (q,r) in
  ;; a range spanning both signs - not just for positive coordinates, as before the fix.
  (let [gr (sut/grid :triangle 10.0)]
    (doseq [q (range -8 9) r (range -8 9)]
      (let [anchor (sut/cell->anchor gr [q r])
            cell-back (sut/coords->cell gr (nth anchor 0) (nth anchor 1))
            expected-q (if (even? q) q (dec q))]
        (t/is (= (double expected-q) (nth cell-back 0)))
        (t/is (= (double r) (nth cell-back 1)))))))

(t/deftest triangle-up-down-pair-closed-form-test
  ;; hand-derived: cells [0 0]/[1 0] and [2 0]/[3 0] are the two up/down pairs sharing
  ;; anchors [0 0] and [10 0] respectively (rhombus->pixel with q>>1 = 0 and 1)
  (let [gr (sut/grid :triangle 10.0)]
    (t/is (v/delta-eq (sut/cell->anchor gr [0 0]) (v/vec3 0.0 0.0 0.0)))
    (t/is (v/delta-eq (sut/cell->anchor gr [1 0]) (v/vec3 0.0 0.0 1.0)))
    (t/is (v/delta-eq (sut/cell->anchor gr [2 0]) (v/vec3 10.0 0.0 0.0)))
    (t/is (v/delta-eq (sut/cell->anchor gr [3 0]) (v/vec3 10.0 0.0 1.0)))))
