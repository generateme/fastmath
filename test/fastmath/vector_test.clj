(ns fastmath.vector-test
  (:require [fastmath.vector :as sut]
            [fastmath.core :as m]
            [clojure.test :as t]))

;; test protocol

(def cv-in1 [-1.0 4.0])
(def cv-in2 [3.0 2.0])

(t/deftest clojure-vector-test
  (t/is (= [-1.0 4.0] (sut/to-vec cv-in1)))
  (t/is (= [0.0 5.0] (sut/fmap cv-in1 inc)))
  (t/is (m/approx-eq 17.0 (sut/magsq cv-in1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag cv-in1)))
  (t/is (m/approx-eq 5.0 (sut/dot cv-in1 cv-in2)))
  (t/is (= cv-in1 (sut/add cv-in1)))
  (t/is (= [2.0 6.0] (sut/add cv-in1 cv-in2)))
  (t/is (= [1.0 -4.0] (sut/sub cv-in1)))
  (t/is (= [-4.0 2.0] (sut/sub cv-in1 cv-in2)))
  (t/is (= [-2.0 8.0] (sut/mult cv-in1 2.0)))
  (t/is (= [-3.0 8.0] (sut/emult cv-in1 cv-in2)))
  (t/is (= [-0.5 2.0] (sut/div cv-in1 2.0)))
  (t/is (= [1.0 4.0] (sut/abs cv-in1)))
  (t/is (== 4.0 (sut/mx cv-in1)))
  (t/is (== -1.0 (sut/mn cv-in1)))
  (t/is (= [3.0 4.0] (sut/emx cv-in1 cv-in2)))
  (t/is (= [-1.0 2.0] (sut/emn cv-in2 cv-in1)))
  (t/is (== 1 (sut/maxdim cv-in1)))
  (t/is (== 0 (sut/mindim cv-in1)))
  (t/is (== 3.0 (sut/sum cv-in1)))
  (t/is (= [4.0 -1.0] (sut/permute cv-in1 [1 0])))
  (t/is (= [-1.0 0.25] (sut/reciprocal cv-in1)))
  (t/is (= [1.0 3.0] (sut/interpolate cv-in1 cv-in2 0.5)))
  (t/is (= [1.0 3.0] (sut/einterpolate cv-in1 cv-in2 [0.5 0.5])))
  (t/is (= [0.0 2.0] (sut/econstrain cv-in1 0.0 2.0)))
  (t/is (not (sut/is-zero? cv-in1)))
  (t/is (sut/is-zero? [0.0 0.0]))
  (t/is (not (sut/is-near-zero? cv-in1)))
  (t/is (sut/is-near-zero? [-0.0000001 0.0])))

(def av-in1 (sut/array-vec cv-in1))
(def av-in2 (sut/array-vec cv-in2))

(t/deftest array-vec-test
  (t/is (== 2 (count av-in1)))
  (t/is (== 4.0 (av-in1 1)))
  (t/is (= [-1.0 4.0] (sut/to-vec av-in1)))
  (t/is (= (sut/array-vec [0.0 5.0]) (sut/fmap av-in1 inc)))
  (t/is (m/approx-eq 17.0 (sut/magsq av-in1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag av-in1)))
  (t/is (m/approx-eq 5.0 (sut/dot av-in1 av-in2)))
  (t/is (= av-in1 (sut/add av-in1)))
  (t/is (= (sut/array-vec [2.0 6.0])  (sut/add av-in1 av-in2)))
  (t/is (= (sut/array-vec [1.0 -4.0]) (sut/sub av-in1)))
  (t/is (= (sut/array-vec [-4.0 2.0]) (sut/sub av-in1 av-in2)))
  (t/is (= (sut/array-vec [-2.0 8.0]) (sut/mult av-in1 2.0)))
  (t/is (= (sut/array-vec [-3.0 8.0]) (sut/emult av-in1 av-in2)))
  (t/is (= (sut/array-vec [-0.5 2.0]) (sut/div av-in1 2.0)))
  (t/is (= (sut/array-vec [1.0 4.0])  (sut/abs av-in1)))
  (t/is (== 4.0 (sut/mx av-in1)))
  (t/is (== -1.0 (sut/mn av-in1)))
  (t/is (= (sut/array-vec [3.0 4.0])  (sut/emx av-in1 av-in2)))
  (t/is (= (sut/array-vec [-1.0 2.0]) (sut/emn av-in2 av-in1)))
  (t/is (== 3.0 (sut/sum av-in1)))
  (t/is (= (sut/array-vec [-1.0 0.25]) (sut/reciprocal av-in1)))
  (t/is (= (sut/array-vec [1.0 3.0])   (sut/interpolate av-in1 av-in2 0.5)))
  (t/is (= (sut/array-vec [1.0 3.0])   (sut/einterpolate av-in1 av-in2 (sut/array-vec [0.5 0.5]))))
  (t/is (= (sut/array-vec [0.0 2.0])   (sut/econstrain av-in1 0.0 2.0)))
  (t/is (not (sut/is-zero? av-in1)))
  (t/is (sut/is-zero? (sut/array-vec [0.0 0.0])))
  (t/is (not (sut/is-near-zero? av-in1)))
  (t/is (sut/is-near-zero? (sut/array-vec [-0.0000001 0.0]))))

;; vec2

(def v2-in1 (sut/vec2 -1.0 4.0))
(def v2-in2 (sut/vec2 3.0 2.0))

(t/deftest vec2-test
  (t/is (== 2 (count v2-in1)))
  (t/is (== 4.0 (v2-in1 1)))
  (t/is (seqable? v2-in1))
  (t/is (== -1.0 (first v2-in1)))
  (t/is (== 4.0 (second v2-in1)))
  (t/is (= [-1.0 4.0] (sut/to-vec v2-in1)))
  (t/is (= (sut/vec2 0.0 5.0) (sut/fmap v2-in1 inc)))
  (t/is (m/approx-eq 17.0 (sut/magsq v2-in1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag v2-in1)))
  (t/is (m/approx-eq 5.0 (sut/dot v2-in1 v2-in2)))
  (t/is (= v2-in1 (sut/add v2-in1)))
  (t/is (= (sut/vec2 2.0 6.0)  (sut/add v2-in1 v2-in2)))
  (t/is (= (sut/vec2 1.0 -4.0) (sut/sub v2-in1)))
  (t/is (= (sut/vec2 -4.0 2.0) (sut/sub v2-in1 v2-in2)))
  (t/is (= (sut/vec2 -2.0 8.0) (sut/mult v2-in1 2.0)))
  (t/is (= (sut/vec2 -3.0 8.0) (sut/emult v2-in1 v2-in2)))
  (t/is (= (sut/vec2 -0.5 2.0) (sut/div v2-in1 2.0)))
  (t/is (= (sut/vec2 1.0 4.0)  (sut/abs v2-in1)))
  (t/is (== 4.0 (sut/mx v2-in1)))
  (t/is (== -1.0 (sut/mn v2-in1)))
  (t/is (= (sut/vec2 3.0 4.0)  (sut/emx v2-in1 v2-in2)))
  (t/is (= (sut/vec2 -1.0 2.0) (sut/emn v2-in2 v2-in1)))
  (t/is (== 1 (sut/maxdim v2-in1)))
  (t/is (== 0 (sut/mindim v2-in1)))
  (t/is (== 3.0 (sut/sum v2-in1)))
  (t/is (= (sut/vec2 4.0 -1.0)  (sut/permute v2-in1 [1 0])))
  (t/is (= (sut/vec2 -1.0 0.25) (sut/reciprocal v2-in1)))
  (t/is (= (sut/vec2 1.0 3.0) (sut/interpolate v2-in1 v2-in2 0.5)))
  (t/is (= (sut/vec2 1.0 3.0) (sut/einterpolate v2-in1 v2-in2 (sut/vec2 0.5 0.5))))
  (t/is (= (sut/vec2 0.0 2.0) (sut/econstrain v2-in1 0.0 2.0)))
  (t/is (not (sut/is-zero? v2-in1)))
  (t/is (sut/is-zero? (sut/vec2 0.0 0.0)))
  (t/is (not(sut/is-near-zero? v2-in1)))
  (t/is (sut/is-near-zero? (sut/vec2 -0.0000001 0.0)))
  (t/is (m/approx-eq 1.815775 (sut/heading v2-in1)))
  (t/is (m/approx-eq -14.0 (sut/cross v2-in1 v2-in2)))
  (t/is (sut/delta-eq (sut/vec2 1.0 -4.0) (sut/rotate v2-in1 m/PI)))
  (t/is (sut/delta-eq (sut/normalize (sut/vec2 -4.0 -1.0)) (sut/perpendicular v2-in1)))
  (t/is (sut/delta-eq (sut/vec2 2.0 -3.0) (sut/transform v2-in1 (sut/vec2 1 1) (sut/vec2 -1 0) (sut/vec2 0 -1))))
  (t/is (= (sut/vec2 1.0 0.0) (sut/to-polar (sut/vec2 1.0 0.0))))
  (t/is (sut/delta-eq (sut/vec2 1.0 0.0) (sut/from-polar (sut/vec2 1.0 0.0)))))

;; vec3

(def v3-in1 (sut/vec3 v2-in1 0.0))
(def v3-in2 (sut/vec3 v2-in2 0.0))

(t/deftest vec3-test
  (t/is (== 3 (count v3-in1)))
  (t/is (== 4.0 (v3-in1 1)))
  (t/is (seqable? v3-in1))
  (t/is (== -1.0 (first v3-in1)))
  (t/is (== 4.0 (second v3-in1)))
  (t/is (= [-1.0 4.0 0.0] (sut/to-vec v3-in1)))
  (t/is (= (sut/vec3 0.0 5.0 1.0) (sut/fmap v3-in1 inc)))
  (t/is (m/approx-eq 17.0 (sut/magsq v3-in1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag v3-in1)))
  (t/is (m/approx-eq 5.0 (sut/dot v3-in1 v3-in2)))
  (t/is (= v3-in1 (sut/add v3-in1)))
  (t/is (= (sut/vec3 2.0 6.0 0.0)  (sut/add v3-in1 v3-in2)))
  (t/is (= (sut/vec3 1.0 -4.0 0.0) (sut/sub v3-in1)))
  (t/is (= (sut/vec3 -4.0 2.0 0.0) (sut/sub v3-in1 v3-in2)))
  (t/is (= (sut/vec3 -2.0 8.0 0.0) (sut/mult v3-in1 2.0)))
  (t/is (= (sut/vec3 -3.0 8.0 0.0) (sut/emult v3-in1 v3-in2)))
  (t/is (= (sut/vec3 -0.5 2.0 0.0) (sut/div v3-in1 2.0)))
  (t/is (= (sut/vec3 1.0 4.0 0.0)  (sut/abs v3-in1)))
  (t/is (== 4.0 (sut/mx v3-in1)))
  (t/is (== -1.0 (sut/mn v3-in1)))
  (t/is (= (sut/vec3 3.0 4.0 0.0)  (sut/emx v3-in1 v3-in2)))
  (t/is (= (sut/vec3 -1.0 2.0 0.0) (sut/emn v3-in2 v3-in1)))
  (t/is (== 1 (sut/maxdim v3-in1)))
  (t/is (== 0 (sut/mindim v3-in1)))
  (t/is (== 3.0 (sut/sum v3-in1)))
  (t/is (= (sut/vec3 4.0 -1.0 0.0) (sut/permute v3-in1 [1 0 2])))
  (t/is (= (sut/vec3 -1.0 0.25 1.0) (sut/reciprocal (sut/vec3 v2-in1 1.0))))
  (t/is (= (sut/vec3 1.0 3.0 0.0) (sut/interpolate v3-in1 v3-in2 0.5)))
  (t/is (= (sut/vec3 1.0 3.0 0.0) (sut/einterpolate v3-in1 v3-in2 (sut/vec3 0.5 0.5 0.5))))
  (t/is (= (sut/vec3 0.0 2.0 0.0) (sut/econstrain v3-in1 0.0 2.0)))
  (t/is (not (sut/is-zero? v3-in1)))
  (t/is (sut/is-zero? (sut/vec3 0.0 0.0 0.0)))
  (t/is (not (sut/is-near-zero? v3-in1)))
  (t/is (sut/is-near-zero? (sut/vec3 -0.0000001 0.0 0.0)))
  (t/is (m/approx-eq 1.815775 (sut/heading v3-in1)))
  (t/is (sut/delta-eq (sut/vec3 0.0 0.0 -14.0) (sut/cross v3-in1 v3-in2)))
  (t/is (sut/delta-eq (sut/vec3 0.0 0.0 -1.0) (sut/perpendicular v3-in1 v3-in2)))
  (t/is (sut/delta-eq (sut/vec3 2.0 -3.0 0.0) (sut/transform v3-in1 (sut/vec3 1 1 0) (sut/vec3 -1 0 0) (sut/vec3 0 -1 0) (sut/vec3 0 0 -1.0)))))

;; rotations
;; from/to-polar


(t/deftest global-fns-test
  (t/is (sut/delta-eq (sut/vec2 -0.3333333333 2.0) (sut/ediv v2-in1 v2-in2)))
  (t/is (= (sut/vec2 1.0 3.0) (sut/average-vectors [v2-in1 v2-in2])))

  (t/is (m/approx-eq (m/sqrt 20.0) (sut/dist v2-in1 v2-in2)))
  (t/is (m/approx-eq 20.0 (sut/dist-sq v2-in1 v2-in2)))
  (t/is (m/approx-eq 6.0 (sut/dist-abs v2-in1 v2-in2)))
  (t/is (m/approx-eq 4.0 (sut/dist-cheb v2-in1 v2-in2)))
  (t/is (m/approx-eq 2.0 (sut/dist-discrete v2-in1 v2-in2)))
  (t/is (m/approx-eq 4.0 (sut/dist-emd v2-in1 v2-in2)))
  (t/is (m/approx-eq 1.3333333 (sut/dist-canberra v2-in1 v2-in2)))
  (t/is (m/approx-eq 0.390812 (sut/dist-ang v2-in1 v2-in2)))

  (t/is (sut/delta-eq (sut/vec2 -0.242535 0.9701425) (sut/normalize v2-in1)))
  (t/is (sut/edelta-eq (sut/vec2 -0.242535 0.9701425) (sut/normalize v2-in1)))
  (t/is (m/approx-eq 0.70710678 (first (sut/set-mag (sut/vec2 1 1) 1))))
  (t/is (== 1.0 (sut/mag (sut/limit v2-in1 1.0))))
  (t/is (m/approx-eq (sut/angle-between v2-in1 v2-in2) (- (sut/relative-angle-between v2-in1 v2-in2))))
  (t/is (not (sut/aligned? v2-in1 v2-in2)))
  (t/is (sut/aligned? v2-in1 (sut/add (sut/vec2 0.0000001 -0.0000001) (sut/mult v2-in1 0.555))))
  (t/is (= v2-in1 (sut/faceforward v2-in1 v2-in2))))

;; clojure contract tests

(defmacro exception?
  [& forms]
  `(try
     ~@forms
     (catch IllegalArgumentException iae# :iae)
     (catch IndexOutOfBoundsException ioobe# :ioobe)))

(defn clojure-contract-vec-tests
  [vfn in]
  (let [v (vfn in)
        cnt (count in)
        lst (dec cnt)
        lin (last in)]
    (t/is (= cnt (count v)))
    (t/is (= in v))
    (t/is (= in (seq v)))
    (t/is (= (reverse in) (rseq v)))
    (t/is (= 2.0 (nth v 1)))
    (t/is (= lin (nth v lst)))
    (t/is (= :ioobe (exception? (nth v cnt))))
    (t/is (= 2.0 (nth v 1 :not-found)))
    (t/is (= :not-found (nth v cnt :not-found)))
    (t/is (= 2.0 (get v 1)))
    (t/is (= lin (get v lst)))
    (t/is (nil? (get v cnt)))
    (t/is (= :not-found (get v cnt :not-found)))
    (t/is (= 2.0 (get v 1 :not-found)))
    (t/is (contains? v 1))
    (t/is (not (contains? v 11)))
    (t/is (= -11.0 (second (assoc v 1 -11.0))))
    (t/is (= :iae (exception? (assoc v :a -11.0))))
    (t/is (= :ioobe (exception? (assoc v cnt -11.0))))
    (t/is (= 2.0 (v 1)))
    (t/is (= lin (v lst)))
    (t/is (= :iae (exception? (v :a))))
    (t/is (= :ioobe (exception? (v cnt))))))

(t/deftest clojure-contract-vecs
  (clojure-contract-vec-tests sut/array-vec [1.0 2.0 3.0 4.0 5.0 6.0])
  (clojure-contract-vec-tests (partial apply sut/vec2) [1.0 2.0])
  (clojure-contract-vec-tests (partial apply sut/vec3) [1.0 2.0 3.0])
  (clojure-contract-vec-tests (partial apply sut/vec4) [1.0 2.0 3.0 4.0]))

;;

(t/deftest similarity
  (t/is (sut/delta-eq [1 2 3 4 5] [1 2 3 4 5]))
  (t/is (not (sut/delta-eq [1 2 3 4 5] [1 2 3 4 5.01])))
  (t/is (sut/delta-eq [1 2 3 4 5] [1 2 3 4 5.01] 1.0e-2))
  (t/is (sut/edelta-eq [1 2 3 4 5] [1 2 3 4 5]))
  (t/is (not (sut/edelta-eq [1 2 3 4 5] [1 2 3 4 5.01])))
  (t/is (sut/edelta-eq [1 2 3 4 5] [1 2 3 4 5.01] 1.0e-2)))

;;

(t/deftest unwrapping
  (t/testing "Python numpy examples"
    (t/is (sut/edelta-eq [0.0 0.78539816 1.57079633 -0.78539816 0.0]
                         (sut/unwrap [0.0 0.78539816 1.57079633 5.49778714 6.28318531] m/TWO_PI)))
    (t/is (sut/edelta-eq [0, 1, 2, 3, 4] (sut/unwrap [0, 1, 2, -1, 0] 4)))
    (t/is (sut/edelta-eq [1, 2, 3, 4, 5, 6, 7, 8, 9] (sut/unwrap [1, 2, 3, 4, 5, 6, 1, 2, 3] 6)))
    (t/is (sut/edelta-eq [2, 3, 4, 5, 6, 7, 8, 9] (sut/unwrap [2, 3, 4, 5, 2, 3, 4, 5] 4)))
    (t/is (sut/edelta-eq [-180., -140., -100.,  -60.,  -20.,   20.,   60.,  100.,  140.,
                          180.,  220.,  260.,  300.,  340.,  380.,  420.,  460.,  500.,
                          540.]
                         (sut/unwrap [-180., -140., -100.,  -60.,  -20.,   20.,   60.,  100.,  140.,
                                      -180., -140., -100.,  -60.,  -20.,   20.,   60.,  100.,  140.,
                                      -180.] 360)))))

;; ==== Fastmath Vector Audit — Group 2.1: Constructors & conversions ====
;; Reference: structural/round-trip checks (constructors/converters, not formulas);
;; behavior confirmed against fastmath.vector source via REPL exploration, 2026-09-23.

(t/deftest vec2-vec3-vec4-constructors
  (t/is (= (sut/vec2 1.0 2.0) (sut/vec2 [1.0 2.0])))
  (t/is (= (sut/vec2 0.0 0.0) (sut/vec2)))
  (t/is (= (sut/vec3 1.0 2.0 3.0) (sut/vec3 (sut/vec2 1.0 2.0) 3.0)))
  (t/is (= (sut/vec3 1.0 2.0 3.0) (sut/vec3 [1.0 2.0 3.0])))
  (t/is (= (sut/vec3 0.0 0.0 0.0) (sut/vec3)))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/vec4 (sut/vec3 1.0 2.0 3.0) 4.0)))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/vec4 (sut/vec2 1.0 2.0) 3.0 4.0)))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/vec4 [1.0 2.0 3.0 4.0])))
  (t/is (= (sut/vec4 0.0 0.0 0.0 0.0) (sut/vec4))))

(t/deftest make-vector-test
  (t/is (nil? (sut/make-vector 0)))
  (t/is (nil? (sut/make-vector -1)))
  (t/is (== 0.0 (sut/make-vector 1)))
  (t/is (= (sut/vec2 0.0 0.0) (sut/make-vector 2)))
  (t/is (= (sut/vec3 0.0 0.0 0.0) (sut/make-vector 3)))
  (t/is (= (sut/vec4 0.0 0.0 0.0 0.0) (sut/make-vector 4)))
  (t/is (= [0.0 0.0 0.0 0.0 0.0] (vec (sut/make-vector 5))))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/make-vector 4 [1 2 3 4])))
  (t/is (== 7.0 (sut/make-vector 1 [7])))
  (t/is (= [1.0 2.0 3.0 4.0 5.0] (vec (sut/make-vector 5 [1 2 3 4 5 6])))))

(t/deftest generate-vec-test
  (t/is (= (sut/vec2 5.0 5.0) (sut/generate-vec2 (constantly 5.0))))
  (let [ctr (atom 0)]
    (t/is (= (sut/vec2 1.0 2.0) (sut/generate-vec2 (fn [] (double (swap! ctr inc)))
                                                     (fn [] (double (swap! ctr inc)))))))
  (t/is (= (sut/vec3 5.0 5.0 5.0) (sut/generate-vec3 (constantly 5.0))))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/generate-vec4 (constantly 1.0) (constantly 2.0)
                                                           (constantly 3.0) (constantly 4.0)))))

(t/deftest array-seq-vec-conversion-test
  (t/is (= (sut/vec2 1.0 2.0) (sut/array->vec2 (double-array [1.0 2.0]))))
  (t/is (= (sut/vec3 1.0 2.0 3.0) (sut/array->vec3 (double-array [1.0 2.0 3.0]))))
  (t/is (= (sut/vec4 1.0 2.0 3.0 4.0) (sut/array->vec4 (double-array [1.0 2.0 3.0 4.0]))))
  (t/is (thrown? ArrayIndexOutOfBoundsException (sut/array->vec3 (double-array [1.0 2.0]))))

  (t/is (= (sut/vec2 1.0 0.0) (sut/seq->vec2 [1.0])))
  (t/is (= (sut/vec2 1.0 2.0) (sut/seq->vec2 [1.0 2.0 3.0])))
  (t/is (= (sut/vec3 0.0 0.0 0.0) (sut/seq->vec3 [])))
  (t/is (= (sut/vec4 1.0 2.0 0.0 0.0) (sut/seq->vec4 [1 2]))))

(t/deftest vec-conversion-helpers-test
  (t/is (= [1.0 2.0] (vec (sut/vec->array (sut/vec2 1.0 2.0)))))
  (t/is (= "[D" (.getName (class (sut/vec->array (sut/vec2 1.0 2.0))))))
  (t/is (= '(1.0 2.0) (sut/vec->seq (sut/vec2 1.0 2.0))))
  (t/is (= '(5.0) (sut/vec->seq 5.0)))
  (t/is (= '(1.0 2.0 3.0) (sut/vec->seq (sut/vec->RealVector [1.0 2.0 3.0]))))
  (t/is (= [1.0 2.0] (sut/vec->vector (sut/vec2 1.0 2.0))))
  (t/is (= [1.0 2.0 3.0] (sut/vec->vector (sut/array-vec [1.0 2.0 3.0]))))
  (t/is (= (sut/vec->RealVector (sut/vec2 1.0 2.0)) (sut/real-vector (sut/vec2 1.0 2.0))))
  (t/is (= org.apache.commons.math3.linear.ArrayRealVector (class (sut/vec->RealVector (sut/vec2 1.0 2.0)))))
  (t/is (= [1.0 2.0] (vec (sut/vec->Vec (sut/vec2 1.0 2.0)))))
  (t/is (= clojure.core.Vec (class (sut/vec->Vec (sut/vec2 1.0 2.0))))))

(t/deftest as-vec-test
  (t/is (= (sut/vec2 9.0 8.0) (sut/as-vec (sut/vec2 1.0 2.0) [9.0 8.0])))
  (t/is (= (sut/vec2 9.0 0.0) (sut/as-vec (sut/vec2 1.0 2.0) [9.0])))
  (t/is (= (sut/vec2 0.0 0.0) (sut/as-vec (sut/vec2 1.0 2.0))))
  (t/is (= (sut/array-vec [9.0 8.0 0.0]) (sut/as-vec (sut/array-vec [1.0 2.0 3.0]) [9.0 8.0])))
  (t/is (= [9.0 8.0 0.0] (sut/as-vec [1.0 2.0 3.0] [9.0 8.0])))
  (t/is (== 9.0 (sut/as-vec 5.0 [9.0 8.0])))
  (t/is (== 0.0 (sut/as-vec 5.0)))
  (t/is (= [9.0 0.0] (vec (sut/as-vec (double-array [1.0 2.0]) [9.0])))))

(t/deftest size-test
  (t/is (== 2 (sut/size (sut/vec2 1.0 2.0))))
  (t/is (== 3 (sut/size (sut/vec3 1.0 2.0 3.0))))
  (t/is (== 4 (sut/size (sut/vec4 1.0 2.0 3.0 4.0))))
  (t/is (== 3 (sut/size (sut/array-vec [1.0 2.0 3.0]))))
  (t/is (== 3 (sut/size [1.0 2.0 3.0])))
  (t/is (== 2 (sut/size (double-array [1.0 2.0]))))
  (t/is (== 1 (sut/size 5.0)))
  (t/is (== 3 (sut/size (sut/vec->RealVector [1.0 2.0 3.0])))))

;; ==== Fastmath Vector Audit — Group 2.2: Shared VectorProto ops ====
;; Reference: structural checks (arithmetic identities on fixture values, cross-checked
;; by hand and via REPL against fastmath.vector source, 2026-09-23); covers the
;; double[], ArrayRealVector, Vec4, Number and true-Seqable (lazy seq) representations
;; not exercised by the original 2.1-era fixtures, plus gap-fill (shift/prod/approx/
;; zero?/near-zero?) for the already-covered representations.

(def da1 (double-array [-1.0 4.0]))
(def da2 (double-array [3.0 2.0]))
(def rv1 (sut/vec->RealVector [-1.0 4.0]))
(def rv2 (sut/vec->RealVector [3.0 2.0]))
(def v4a (sut/vec4 -1.0 4.0 0.0 2.0))
(def v4b (sut/vec4 3.0 2.0 1.0 -1.0))

(t/deftest double-array-vectorproto-test
  (t/is (m/approx-eq 17.0 (sut/magsq da1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag da1)))
  (t/is (m/approx-eq 5.0 (sut/dot da1 da2)))
  (t/is (= [2.0 6.0] (seq (sut/add da1 da2))))
  (t/is (= [-4.0 2.0] (seq (sut/sub da1 da2))))
  (t/is (= [1.0 6.0] (seq (sut/shift da1 2.0))))
  (t/is (= [-2.0 8.0] (seq (sut/mult da1 2.0))))
  (t/is (= [-3.0 8.0] (seq (sut/emult da1 da2))))
  (t/is (= [1.0 4.0] (seq (sut/abs da1))))
  (t/is (= [-1.0 0.25] (seq (sut/reciprocal da1))))
  (t/is (== 4.0 (sut/mx da1)))
  (t/is (== -1.0 (sut/mn da1)))
  (t/is (= [3.0 4.0] (seq (sut/emx da1 da2))))
  (t/is (= [-1.0 2.0] (seq (sut/emn da1 da2))))
  (t/is (== 1 (sut/maxdim da1)))
  (t/is (== 0 (sut/mindim da1)))
  (t/is (== 3.0 (sut/sum da1)))
  (t/is (== -4.0 (sut/prod da1)))
  (t/is (= [-1.2 4.6] (seq (sut/approx (double-array [-1.234 4.567]) 1))))
  (t/is (= [0.0 5.0] (seq (sut/fmap da1 inc))))
  (t/is (not (sut/is-zero? da1)))
  (t/is (sut/is-zero? (double-array [0.0 0.0])))
  (t/is (not (sut/zero? da1)))
  (t/is (sut/near-zero? (double-array [0.0000001])))
  (t/is (= [1.0 3.0] (seq (sut/interpolate da1 da2 0.5))))
  (t/is (= [1.0 3.0] (seq (sut/einterpolate da1 da2 (double-array [0.5 0.5])))))
  (t/is (= [0.0 2.0] (seq (sut/econstrain da1 0.0 2.0)))))

(t/deftest array-real-vector-vectorproto-test
  (t/is (m/approx-eq 17.0 (sut/magsq rv1)))
  (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag rv1)))
  (t/is (m/approx-eq 5.0 (sut/dot rv1 rv2)))
  (t/is (= (sut/vec->RealVector [2.0 6.0]) (sut/add rv1 rv2)))
  (t/is (= (sut/vec->RealVector [-4.0 2.0]) (sut/sub rv1 rv2)))
  (t/is (= (sut/vec->RealVector [1.0 6.0]) (sut/shift rv1 2.0)))
  (t/is (= (sut/vec->RealVector [-2.0 8.0]) (sut/mult rv1 2.0)))
  (t/is (= (sut/vec->RealVector [-3.0 8.0]) (sut/emult rv1 rv2)))
  (t/is (= (sut/vec->RealVector [1.0 4.0]) (sut/abs rv1)))
  (t/is (= (sut/vec->RealVector [-1.0 0.25]) (sut/reciprocal rv1)))
  (t/is (== 4.0 (sut/mx rv1)))
  (t/is (== -1.0 (sut/mn rv1)))
  (t/is (= (sut/vec->RealVector [3.0 4.0]) (sut/emx rv1 rv2)))
  (t/is (= (sut/vec->RealVector [-1.0 2.0]) (sut/emn rv1 rv2)))
  (t/is (== 1 (sut/maxdim rv1)))
  (t/is (== 0 (sut/mindim rv1)))
  (t/is (== 3.0 (sut/sum rv1)))
  (t/is (== -4.0 (sut/prod rv1)))
  (t/is (= (sut/vec->RealVector [-1.23 4.57]) (sut/approx (sut/vec->RealVector [-1.234 4.567]) 2)))
  (t/is (= (sut/vec->RealVector [0.0 5.0]) (sut/fmap rv1 inc)))
  (t/is (not (sut/is-zero? rv1)))
  (t/is (sut/is-zero? (sut/vec->RealVector [0.0 0.0])))
  (t/is (not (sut/zero? rv1)))
  (t/is (= (sut/vec->RealVector [1.0 3.0]) (sut/interpolate rv1 rv2 0.5)))
  (t/is (= (sut/vec->RealVector [1.0 3.0]) (sut/einterpolate rv1 rv2 (sut/vec->RealVector [0.5 0.5]))))
  (t/testing "econstrain preserves dimension (regression for the ArrayRealVector concatenating-constructor bug)"
    (let [c (sut/econstrain rv1 0.0 2.0)]
      (t/is (= (sut/vec->RealVector [0.0 2.0]) c))
      (t/is (== 2 (sut/size c))))))

(t/deftest vec4-vectorproto-test
  (t/is (m/approx-eq 21.0 (sut/magsq v4a)))
  (t/is (m/approx-eq (m/sqrt 21.0) (sut/mag v4a)))
  (t/is (m/approx-eq 3.0 (sut/dot v4a v4b)))
  (t/is (= (sut/vec4 2.0 6.0 1.0 1.0) (sut/add v4a v4b)))
  (t/is (= (sut/vec4 -4.0 2.0 -1.0 3.0) (sut/sub v4a v4b)))
  (t/is (= (sut/vec4 1.0 6.0 2.0 4.0) (sut/shift v4a 2.0)))
  (t/is (= (sut/vec4 -2.0 8.0 0.0 4.0) (sut/mult v4a 2.0)))
  (t/is (= (sut/vec4 -3.0 8.0 0.0 -2.0) (sut/emult v4a v4b)))
  (t/is (= (sut/vec4 1.0 4.0 0.0 2.0) (sut/abs v4a)))
  (t/is (= (sut/vec4 -1.0 0.25 ##Inf 0.5) (sut/reciprocal v4a)))
  (t/is (== 4.0 (sut/mx v4a)))
  (t/is (== -1.0 (sut/mn v4a)))
  (t/is (= (sut/vec4 3.0 4.0 1.0 2.0) (sut/emx v4a v4b)))
  (t/is (= (sut/vec4 -1.0 2.0 0.0 -1.0) (sut/emn v4a v4b)))
  (t/is (== 1 (sut/maxdim v4a)))
  (t/is (== 0 (sut/mindim v4a)))
  (t/is (== 5.0 (sut/sum v4a)))
  (t/is (== -0.0 (sut/prod v4a)))
  (t/is (= (sut/vec4 -1.23 4.57 0.0 2.0) (sut/approx (sut/vec4 -1.234 4.567 0.0 2.0) 2)))
  (t/is (= (sut/vec4 0.0 5.0 1.0 3.0) (sut/fmap v4a inc)))
  (t/is (not (sut/is-zero? v4a)))
  (t/is (sut/is-zero? (sut/vec4)))
  (t/is (not (sut/zero? v4a)))
  (t/is (= (sut/vec4 4.0 -1.0 2.0 0.0) (sut/permute v4a [1 0 3 2])))
  (t/is (= (sut/vec4 1.0 3.0 0.5 0.5) (sut/interpolate v4a v4b 0.5)))
  (t/is (= (sut/vec4 1.0 3.0 0.5 0.5) (sut/einterpolate v4a v4b (sut/vec4 0.5 0.5 0.5 0.5))))
  (t/is (= (sut/vec4 0.0 2.0 0.0 2.0) (sut/econstrain v4a 0.0 2.0))))

(t/deftest number-vectorproto-test
  (t/is (m/approx-eq 25.0 (sut/magsq 5.0)))
  (t/is (m/approx-eq 5.0 (sut/mag -5.0)))
  (t/is (m/approx-eq 12.0 (sut/dot 3.0 4.0)))
  (t/is (== 7.0 (sut/add 3.0 4.0)))
  (t/is (== -1.0 (sut/sub 3.0 4.0)))
  (t/is (== 7.0 (sut/shift 3.0 4.0)))
  (t/is (== 12.0 (sut/mult 3.0 4.0)))
  (t/is (== 12.0 (sut/emult 3.0 4.0)))
  (t/is (== 5.0 (sut/abs -5.0)))
  (t/is (== 0.25 (sut/reciprocal 4.0)))
  (t/is (== 5.0 (sut/mx 5.0)))
  (t/is (== 5.0 (sut/mn 5.0)))
  (t/is (== 4.0 (sut/emx 3.0 4.0)))
  (t/is (== 3.0 (sut/emn 3.0 4.0)))
  (t/is (== 0 (sut/maxdim 5.0)))
  (t/is (== 0 (sut/mindim 5.0)))
  (t/is (== 5.0 (sut/sum 5.0)))
  (t/is (== 5.0 (sut/prod 5.0)))
  (t/is (== 1.23 (sut/approx 1.23456 2)))
  (t/is (== 6.0 (sut/fmap 5.0 inc)))
  (t/is (sut/is-zero? 0.0))
  (t/is (not (sut/is-zero? 5.0)))
  (t/is (sut/zero? 0.0))
  (t/is (sut/near-zero? 0.0000001))
  (t/is (== 2.0 (sut/interpolate 1.0 3.0 0.5)))
  (t/is (== 2.0 (sut/econstrain 5.0 0.0 2.0)))
  (t/testing "permute is not implemented for scalars (no notion of element order)"
    (t/is (thrown? IllegalArgumentException (sut/permute 5.0 [0])))))

(t/deftest seqable-lazy-seq-vectorproto-test
  (let [sq1 (map identity [-1.0 4.0])]
    (t/is (not (instance? clojure.lang.IPersistentVector sq1)))
    (t/is (m/approx-eq 17.0 (sut/magsq sq1)))
    (t/is (m/approx-eq (m/sqrt 17.0) (sut/mag sq1)))
    (t/is (= [1.0 6.0] (seq (sut/shift sq1 2.0))))
    (t/is (== -4.0 (sut/prod sq1)))
    (t/is (= [-1.2 4.6] (seq (sut/approx (map identity [-1.234 4.567]) 1))))
    (t/is (not (sut/zero? sq1)))
    (t/is (sut/near-zero? (map identity [0.0000001])))))

(t/deftest gap-fill-shift-prod-approx-zero-test
  (t/is (= (sut/array-vec [1.0 6.0]) (sut/shift av-in1 2.0)))
  (t/is (== -4.0 (sut/prod av-in1)))
  (t/is (= av-in1 (sut/approx av-in1 1)))
  (t/is (not (sut/zero? av-in1)))
  (t/is (sut/near-zero? (sut/array-vec [0.0000001 0.0])))

  (t/is (= (sut/vec2 1.0 6.0) (sut/shift v2-in1 2.0)))
  (t/is (== -4.0 (sut/prod v2-in1)))
  (t/is (= v2-in1 (sut/approx v2-in1 1)))
  (t/is (not (sut/zero? v2-in1)))
  (t/is (sut/near-zero? (sut/vec2 0.0000001 0.0)))

  (t/is (= (sut/vec3 1.0 6.0 2.0) (sut/shift v3-in1 2.0)))
  (t/is (== -0.0 (sut/prod v3-in1)))
  (t/is (= v3-in1 (sut/approx v3-in1 1)))
  (t/is (not (sut/zero? v3-in1)))
  (t/is (sut/near-zero? (sut/vec3 0.0000001 0.0 0.0))))

(t/deftest heading-semantics-test
  (t/testing "Vec2: signed, full-circle atan2"
    (t/is (m/approx-eq 1.8157749899217608 (sut/heading (sut/vec2 -1.0 4.0))))
    (t/is (m/approx-eq -1.8157749899217608 (sut/heading (sut/vec2 -1.0 -4.0)))))
  (t/testing "Vec3/n-dim: unsigned angle-from-primary-axis, degenerate for dim >= 3"
    (t/is (m/approx-eq 1.8157749899217608 (sut/heading (sut/vec3 -1.0 4.0 0.0))))
    (t/is (= (sut/heading (sut/vec3 -1.0 4.0 0.0)) (sut/heading (sut/vec3 -1.0 -4.0 0.0))))
    (t/is (m/approx-eq 0.0 (sut/heading [1.0 0.0 0.0 0.0])))
    (t/is (m/approx-eq m/HALF_PI (sut/heading (double-array [0.0 1.0]))))))

;; ==== Fastmath Vector Audit — Group 2.3: Type-exclusive geometric ops ====
;; Reference: structural/identity checks (orthogonality via dot products, rotate/axis-rotate
;; agreement, to-polar/from-polar round-trip and known-pole values), confirmed via REPL
;; against fastmath.vector source, 2026-09-23. Vec4 implements none of cross/rotate/
;; axis-rotate/perpendicular/transform/to-polar/from-polar/base-from (confirmed by source
;; inspection and editor diagnostics) — out of scope for this group (WONT, not applicable).

(t/deftest base-from-test
  (t/is (= [v2-in1 (sut/perpendicular v2-in1)] (sut/base-from v2-in1)))
  (let [[a b c] (sut/base-from v3-in1)]
    (t/is (= v3-in1 a))
    (t/is (sut/delta-eq (sut/vec3 0.0 0.0 -1.0) b))
    (t/is (m/approx-eq 1.0 (sut/mag b)))
    (t/is (sut/delta-eq (sut/cross v3-in1 b) c))
    (t/is (m/approx-eq 0.0 (sut/dot a b)))
    (t/is (m/approx-eq 0.0 (sut/dot a c)))
    (t/is (m/approx-eq 0.0 (sut/dot b c)))))

(t/deftest axis-rotate-test
  (t/is (sut/delta-eq (sut/rotate v3-in1 0.0 0.0 m/HALF_PI)
                      (sut/axis-rotate v3-in1 m/HALF_PI (sut/vec3 0.0 0.0 1.0))))
  (let [pivot (sut/vec3 1.0 1.0 0.0)]
    (t/is (sut/delta-eq (sut/add (sut/axis-rotate (sut/sub v3-in1 pivot) m/HALF_PI (sut/vec3 0.0 0.0 1.0)) pivot)
                        (sut/axis-rotate v3-in1 m/HALF_PI (sut/vec3 0.0 0.0 1.0) pivot)))))

(t/deftest vec3-polar-test
  (t/is (sut/delta-eq v3-in1 (sut/from-polar (sut/to-polar v3-in1))))
  (t/is (sut/delta-eq (sut/vec3 1.0 0.0 0.0) (sut/to-polar (sut/vec3 0.0 0.0 1.0))))
  (t/is (sut/delta-eq (sut/vec3 1.0 m/PI 0.0) (sut/to-polar (sut/vec3 0.0 0.0 -1.0))))
  (t/is (sut/delta-eq (sut/vec3 1.0 m/HALF_PI 0.0) (sut/to-polar (sut/vec3 1.0 0.0 0.0)))))
