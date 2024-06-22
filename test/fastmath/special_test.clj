(ns fastmath.special-test
  (:require [fastmath.special :as sut]
            [clojure.test :as t]
            [clojisr.v1.r :as rr]
            [fastmath.core :as m]
            [fastmath.vector :as v]))

(defn init-r []
  (rr/discard-all-sessions)
  (rr/require-r '[Bessel] '[base]))

(init-r)

(defn with-r-session [f]
  (init-r)
  (f)
  (rr/discard-all-sessions))

(t/use-fixtures :once with-r-session)

(def x (range 1.0e-6 100.0 0.01))
(def -x (map - x))
(def xl [1e12, 5e12, 1e13, 5e13, 1e14, 5e14, 1e15, 5e15, 1e16, 5e16, 1e17, 5e17, 1e18, 5e18, 1e19, 5e19, 1e20, 1e22, 1e25, 1e30, 1e40])

(t/deftest bessel-J0
  (t/is (m/nan? (sut/bessel-j0 ##NaN)))
  (t/is (m/one? (sut/bessel-j0 0.0)))
  (t/is (m/zero? (sut/bessel-j0 ##Inf)))
  (t/is (m/zero? (sut/bessel-j0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-j0 x) (rr/r->clj (base/besselJ x 0.0)) 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-j0 x) (rr/r->clj (Bessel/BesselJ x 0.0)) 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-j0 -x) (rr/r->clj (Bessel/BesselJ -x 0.0)) 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-j0 xl) [1.016712505004068e-7
                                             -2.1276975389854557e-7
                                             1.192648473966565e-7
                                             -7.094408384899425e-8
                                             -6.698265203680451e-8
                                             -3.4394958970536735e-8
                                             6.156638646885022e-9
                                             -1.0644174809027939e-8
                                             8.661427680921681e-10
                                             -1.8463226228625315e-9
                                             -2.4087235483673835e-9
                                             6.690193890936787e-11
                                             -4.934387036790141e-10
                                             -6.814178646566312e-11
                                             -2.3228731060101177e-10
                                             -4.7522631032150424e-11
                                             6.6980090407034224e-12
                                             -1.856105106510822e-12
                                             1.1543496219672643e-13
                                             -5.589003016686146e-16
                                             -6.538288347442135e-22] 1.0e-8 1.0e-8)))

(t/deftest bessel-J1
  (t/is (m/nan? (sut/bessel-j1 ##NaN)))
  (t/is (m/zero? (sut/bessel-j1 0.0)))
  (t/is (m/zero? (sut/bessel-j0 ##Inf)))
  (t/is (m/zero? (sut/bessel-j0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-j1 x) (rr/r->clj (base/besselJ x 1.0)) 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-j1 x) (rr/r->clj (Bessel/BesselJ x 1.0)) 1.0e-15))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-j1 -x) (rr/r->clj (Bessel/BesselJ -x 1.0)) 1.0e-15))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-j1 xl) [-7.913802683850442e-7
                                             2.8644892441665137e-7
                                             -2.2234629165382475e-7
                                             -8.774583986821614e-8
                                             4.3353454877231535e-8
                                             -9.498754768402967e-9
                                             2.4468665123771328e-8
                                             -3.745063031294846e-9
                                             7.931694266803266e-9
                                             -3.0534387532186807e-9
                                             7.511648229358563e-10
                                             -1.1263941030142736e-9
                                             -6.270071914094412e-10
                                             -3.5025797836849376e-10
                                             -9.85118397478559e-11
                                             1.0234253752537956e-10
                                             -7.95068198242545e-11
                                             -7.759951744073064e-12
                                             -2.2435852276969216e-13
                                             5.694297368089524e-16
                                             7.952011386537066e-21] 1.0e-8 1.0e-8)))

(t/deftest bessel-J
  (t/is (m/nan? (sut/bessel-j 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-j -0.5 ##NaN)))
  (t/is (every? m/zero? (map #(sut/bessel-j % 0.0) (range 0.2 100.1 0.1))))
  (t/is (every? m/zero? (map #(sut/bessel-j % 0.0) (range -100 0))))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05]
          nu [2, 4, 6, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-14)))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5, 98.5, 104.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-13)))
  (doseq [nu [150, 165.2, 200.0, 300.0, 500.0, 1000.0, 5000.2, 10000.0, 50000.0]
          x [0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92,0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,0.995, 0.999, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-j nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-12)))
  (let [vs (range 0.0 250.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-j % 0.15) vs) (rr/r->clj (base/besselJ 0.15 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 2.1) vs) (rr/r->clj (base/besselJ 2.1 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 42.1) vs) (rr/r->clj (base/besselJ 42.1 vs)) 1.0e-14))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 142.1) vs) (rr/r->clj (base/besselJ 142.1 vs)) 1.0e-14)))
  (let [vs (range -100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-j % 0.15) vs) (rr/r->clj (base/besselJ 0.15 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 2.1) vs) (rr/r->clj (base/besselJ 2.1 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 42.1) vs) (rr/r->clj (base/besselJ 42.1 vs)) 1.0e-14))
    (t/is (v/edelta-eq (map #(sut/bessel-j % 142.1) vs) (rr/r->clj (base/besselJ 142.1 vs)) 1.0e-14)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -100 0 5)
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-j nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-14 1.0e-14)))
  (t/are [v x] (m/delta-eq (sut/bessel-j v x) (first (rr/r->clj (Bessel/BesselJ x v))) 1.0e-15)
    -5.0 -5.1
    -7.3 19.1
    -14.0 21.3
    -13.0 21.3
    -14.0 -21.3
    -13.0 -21.3
    7.3 19.1
    14.0 21.3
    13.0 21.3
    14.0 -21.3
    13.0 -21.3))

;;;;

(t/deftest bessel-Y0
  (t/is (m/nan? (sut/bessel-y0 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-y0 0.0)))
  (t/is (m/zero? (sut/bessel-y0 ##Inf)))
  (t/is (m/nan? (sut/bessel-y0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-y0 x) (rr/r->clj (base/besselY x 0.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-y0 x) (rr/r->clj (Bessel/BesselY x 0.0)) 1.0e-14 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-y0 xl) [ -7.91380268385095e-7
                                             2.8644892441667265e-7
                                             -2.223462916538307e-7
                                             -8.774583986821542e-8
                                             4.335345487723187e-8
                                             -9.498754768402932e-9
                                             2.446866512377132e-8
                                             -3.745063031294845e-9
                                             7.931694266803265e-9
                                             -3.0534387532186807e-9
                                             7.511648229358563e-10
                                             -1.1263941030142736e-9
                                             -6.270071914094412e-10
                                             -3.5025797836849376e-10
                                             -9.85118397478559e-11
                                             1.0234253752537956e-10
                                             -7.95068198242545e-11
                                             -7.759951744073064e-12
                                             -2.2435852276969216e-13
                                             5.694297368089524e-16
                                             7.952011386537066e-21] 1.0e-8 1.0e-8)))

(t/deftest bessel-Y1
  (t/is (m/nan? (sut/bessel-y1 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-y1 0.0)))
  (t/is (m/zero? (sut/bessel-y1 ##Inf)))
  (t/is (m/nan? (sut/bessel-y1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-y1 x) (rr/r->clj (base/besselY x 1.0)) 1.0e-15 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-y1 x) (rr/r->clj (Bessel/BesselY x 1.0)) 1.0e-14 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-y1 xl) [-1.0167125050080249e-7
                                             2.127697538985742e-7
                                             -1.1926484739666762e-7
                                             7.094408384899338e-8
                                             6.698265203680473e-8
                                             3.439495897053673e-8
                                             -6.15663864688501e-9
                                             1.0644174809027939e-8
                                             -8.661427680921677e-10
                                             1.8463226228625315e-9
                                             2.4087235483673835e-9
                                             -6.690193890936787e-11
                                             4.934387036790141e-10
                                             6.814178646566312e-11
                                             2.3228731060101177e-10
                                             4.7522631032150424e-11
                                             -6.6980090407034224e-12
                                             1.856105106510822e-12
                                             -1.1543496219672643e-13
                                             5.589003016686146e-16
                                             6.538288347442135e-22] 1.0e-8 1.0e-8)))

(t/deftest bessel-Y
  (t/is (m/nan? (sut/bessel-y 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-y -0.5 ##NaN)))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 1.9, 2.5, 3.0, 3.5, 5.0, 10.0]
          nu [0, 1, 2, 4, 6, 10, 15, 20, 25, 30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 125, 150, 175, 200]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-y nu xx) (first (rr/r->clj (Bessel/BesselY xx nu))) 1.0e-12 1.0e-12)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-y nu xx) (first (rr/r->clj (Bessel/BesselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-y nu x) (first (rr/r->clj (base/besselY x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-y nu x) (first (rr/r->clj (Bessel/BesselY x nu))) 1.0e-13 1.0e-13)))
  (doseq [nu [150, 165.2, 200.0, 300.0, 500.0, 1000.0, 5000.2, 10000.0 20000.0]
          x [0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92,0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,0.995, 0.999, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-10 1.0e-10)))
  (let [vs (range 0.0 100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-y % 0.15) vs) (rr/r->clj (base/besselY 0.15 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-y % 2.1) vs) (rr/r->clj (base/besselY 2.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-y % 42.1) vs) (rr/r->clj (base/besselY 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-y % 142.1) vs) (rr/r->clj (base/besselY 142.1 vs)) 1.0e-14 1.0e-14)))
  (doseq [v (range -100.0 0.0 0.25)]
    (t/is (m/delta-eq (sut/bessel-y v 0.15) (first (rr/r->clj (Bessel/BesselY 0.15 v))) 1.0e-6 1.0e-6))
    (t/is (m/delta-eq (sut/bessel-y v 2.1) (first (rr/r->clj (Bessel/BesselY 2.1 v))) 1.0e-8 1.0e-8))
    (t/is (m/delta-eq (sut/bessel-y v 42.1) (first (rr/r->clj (Bessel/BesselY 42.1 v))) 1.0e-9 1.0e-9))
    (t/is (m/delta-eq (sut/bessel-y v 142.1) (first (rr/r->clj (Bessel/BesselY 142.1 v))) 1.0e-12 1.0e-12)))
  (t/are [v x] (m/delta-eq (sut/bessel-y v x) (first (rr/r->clj (Bessel/BesselY x v))) 1.0e-14 1.0e-14)
    -6.2 18.6
    -8.0 23.2
    -7.0 23.2
    -6.0 23.2
    -0.1 2.2))

;;

(t/deftest bessel-I0
  (t/is (m/nan? (sut/bessel-i0 ##NaN)))
  (t/is (m/one? (sut/bessel-i0 0.0)))
  (t/is (m/nan? (sut/bessel-i0 ##Inf)))
  (t/is (m/nan? (sut/bessel-i0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-i0 x) (rr/r->clj (base/besselI x 0.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-i0 x) (rr/r->clj (Bessel/BesselI x 0.0)) 1.0e-14 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-i0 -x) (rr/r->clj (Bessel/BesselI -x 0.0)) 1.0e-14 1.0e-14)))

(t/deftest bessel-I1
  (t/is (m/nan? (sut/bessel-i1 ##NaN)))
  (t/is (m/zero? (sut/bessel-i1 0.0)))
  (t/is (m/nan? (sut/bessel-i1 ##Inf)))
  (t/is (m/nan? (sut/bessel-i1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-i1 x) (rr/r->clj (base/besselI x 1.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-i1 x) (rr/r->clj (Bessel/BesselI x 1.0)) 1.0e-14 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-i1 -x) (rr/r->clj (Bessel/BesselI -x 1.0)) 1.0e-14 1.0e-14)))

;;

(t/deftest bessel-I
  (t/is (m/nan? (sut/bessel-i 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-i -0.5 ##NaN)))
  (doseq [x [0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 1.9, 2.5, 3.0, 3.5, 4.0]
          nu [0.01,0.1, 0.5, 0.8, 1, 1.23, 2,2.56, 4,5.23, 6,9.2, 10,12.89, 15, 19.1, 20, 25, 30, 33.123, 40, 45, 50, 51.5, 55, 60, 65, 70, 72.34, 75, 80, 82.1, 85, 88.76, 90, 92.334, 95, 99.87,100, 110, 125, 145.123, 150, 160.789]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-i nu xx) (first (rr/r->clj (base/besselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-i nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-i nu x) (first (rr/r->clj (base/besselI x nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-i nu x) (first (rr/r->clj (Bessel/BesselI x nu))) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-i nu xx) (first (rr/r->clj (base/besselI xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-i nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-i nu x) (first (rr/r->clj (base/besselI x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-i nu x) (first (rr/r->clj (Bessel/BesselI x nu))) 1.0e-13 1.0e-13)))
  (let [vs (range 0.0 250.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-i % 0.15) vs) (rr/r->clj (base/besselI 0.15 vs)) 1.0e-15 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 2.1) vs) (rr/r->clj (base/besselI 2.1 vs)) 1.0e-15 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 42.1) vs) (rr/r->clj (base/besselI 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 142.1) vs) (rr/r->clj (base/besselI 142.1 vs)) 1.0e-12 1.0e-12)))
  (let [vs (range -100.0 0.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-i % 0.15) vs) (rr/r->clj (base/besselI 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 2.1) vs) (rr/r->clj (base/besselI 2.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 42.1) vs) (rr/r->clj (base/besselI 42.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-i % 142.1) vs) (rr/r->clj (base/besselI 142.1 vs)) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -100 0 5)
          :let [xx (m/* x nu)
                axx (m/abs xx)]]
    (t/is (m/delta-eq (sut/bessel-i nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-i nu axx) (first (rr/r->clj (Bessel/BesselI axx nu))) 1.0e-12 1.0e-12)))
  (t/are [v x] (m/delta-eq (sut/bessel-i v x) (first (rr/r->clj (Bessel/BesselI x v))) 1.0e-14 1.0e-14)
    12.0 3.2
    13.0 -1.0
    -8.0 4.2
    12.3 8.2
    -12.3 8.2
    -14.0 -9.9))

;;

(t/deftest bessel-K0
  (t/is (m/nan? (sut/bessel-k0 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-k0 0.0)))
  (t/is (m/zero? (sut/bessel-k0 ##Inf)))
  (t/is (m/nan? (sut/bessel-k0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-k0 x) (rr/r->clj (base/besselK x 0.0)) 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-k0 x) (rr/r->clj (Bessel/BesselK x 0.0)) 1.0e-14)))

(t/deftest bessel-K1
  (t/is (m/nan? (sut/bessel-k1 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-k1 0.0)))
  (t/is (m/zero? (sut/bessel-k1 ##Inf)))
  (t/is (m/nan? (sut/bessel-k1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-k1 x) (rr/r->clj (base/besselK x 1.0)) 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-k1 x) (rr/r->clj (Bessel/BesselK x 1.0)) 1.0e-14)))

(t/deftest bessel-K
  (t/is (m/nan? (sut/bessel-k 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-k -0.5 ##NaN)))
  (doseq [nu (range -36.0 82.0 0.87654)
          x (rest (m/slice-range 0 30 31))]
    (t/is (m/delta-eq (sut/bessel-k nu x) (first (rr/r->clj (base/besselK x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-k nu x) (first (rr/r->clj (Bessel/BesselK x nu))) 1.0e-12 1.0e-12)))
  (doseq [x [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-k nu xx) (first (rr/r->clj (base/besselK xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-k nu xx) (first (rr/r->clj (Bessel/BesselK xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-k nu x) (first (rr/r->clj (base/besselK x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-k nu x) (first (rr/r->clj (Bessel/BesselK x nu))) 1.0e-13 1.0e-13)))
  (let [vs (range 0.0 100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-k % 0.15) vs) (rr/r->clj (base/besselK 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 2.1) vs) (rr/r->clj (base/besselK 2.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 42.1) vs) (rr/r->clj (base/besselK 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 142.1) vs) (rr/r->clj (base/besselK 142.1 vs)) 1.0e-15 1.0e-15)))
  (let [vs (range -100.0 0.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-k % 0.15) vs) (rr/r->clj (base/besselK 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 2.1) vs) (rr/r->clj (base/besselK 2.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 42.1) vs) (rr/r->clj (base/besselK 42.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-k % 142.1) vs) (rr/r->clj (base/besselK 142.1 vs)) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -50 0 0.25)
          :let [xx (m/abs (* nu x))]]
    (t/is (m/delta-eq (sut/bessel-k nu xx) (first (rr/r->clj (Bessel/BesselK xx nu))) 1.0e-13 1.0e-13)))
  (t/are [v x] (m/delta-eq (sut/bessel-k v x) (first (rr/r->clj (Bessel/BesselK x v))) 1.0e-14 1.0e-14)
    12.0 3.2
    -8.0 4.2
    12.3 8.2
    -12.3 8.2))

(t/deftest bessel-k-half
  (doseq [o (range 1 100 2)
          x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          :let [xx (* o x)
                oh (* o 0.5)]]
    (t/is (m/delta-eq (sut/bessel-k oh xx) (sut/bessel-k-half o xx) 1.0e-13 1.0e-13))))

;;

;; Abramowitz and Stegun p.511
(t/deftest kummers-m
  (t/are [a b x res] (m/delta-eq res (sut/kummers-M a b x))
    0.3 0.2 -0.1 0.8578490
    -0.1 0.2 0.1 (* 0.8578490 (m/pow m/E 0.1))
    17 16 1 2.8881744
    -1 16 -1 1.0625
    -1.3 0.2 0.1 0.3582123
    -1.3 1.2 0.1 0.8924108
    -0.3 1.2 0.1 0.9745952
    1 1 1 m/E
    2 2 2 (m/exp 2)
    0.3 0.4 0.5 (/ 1.724128 (/ 0.7 0.6))
    ;; p.533
    -1 1 9 -8
    -1 0.6 9 -14
    0 1 9 1
    1 1 9 (m/exp 9))
  (t/are [a b x res acc] (m/delta-eq res (sut/kummers-M a b x) acc)
    0.9 0.1 10 1227235 1
    -52.5 0.1 1 -16.34 0.2)
  (t/is (m/pos-inf? (sut/kummers-M 1 0 1)))
  (t/is (m/neg-inf? (sut/kummers-M -1 0 1)))
  (t/is (m/pos-inf? (sut/kummers-M -1 0 -1)))
  (t/is (m/neg-inf? (sut/kummers-M 1 0 -1)))
  (t/is (m/nan? (sut/kummers-M -1 -1 2)))) ;; example 4

(t/deftest whittaker-m
  (t/is (m/delta-eq 1.10622 (sut/whittaker-M 0 -0.4 1) 1.0e-5)))

(t/deftest besselk
  (t/are [order ress] (v/delta-eq ress (mapv (partial sut/bessel-k-half order) [0.5 1 1.33 2.5 5]))
    1 [1.075047603 0.461068504 0.287423621 0.065065943 0.003776613]
    3 [3.225142810 0.922137009 0.503531608 0.091092320 0.004531936]
    5 [20.425904466  3.227479531  1.423209202  0.174376728  0.006495775]
    7 [207.48418748  17.05953466   5.85394214   0.43984578   0.01102771]
    9 [2.925204529e+03 1.226442222e+02 3.223343101e+01 1.405944900e+00 2.193457048e-02]))
