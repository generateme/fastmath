(ns fastmath.special-test
  (:require [fastmath.special :as sut]
            [clojure.test :as t]
            [clojisr.v1.r :as rr]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]))

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
  (t/is (m/nan? (sut/bessel-J0 ##NaN)))
  (t/is (m/one? (sut/bessel-J0 0.0)))
  (t/is (m/zero? (sut/bessel-J0 ##Inf)))
  (t/is (m/zero? (sut/bessel-J0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-J0 x) (rr/r->clj (base/besselJ x 0.0)) 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-J0 x) (rr/r->clj (Bessel/BesselJ x 0.0)) 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-J0 -x) (rr/r->clj (Bessel/BesselJ -x 0.0)) 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-J0 xl) [1.016712505004068e-7
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
  (t/is (m/nan? (sut/bessel-J1 ##NaN)))
  (t/is (m/zero? (sut/bessel-J1 0.0)))
  (t/is (m/zero? (sut/bessel-J0 ##Inf)))
  (t/is (m/zero? (sut/bessel-J0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-J1 x) (rr/r->clj (base/besselJ x 1.0)) 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-J1 x) (rr/r->clj (Bessel/BesselJ x 1.0)) 1.0e-15))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-J1 -x) (rr/r->clj (Bessel/BesselJ -x 1.0)) 1.0e-15))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-J1 xl) [-7.913802683850442e-7
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
  (t/is (m/nan? (sut/bessel-J 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-J -0.5 ##NaN)))
  (t/is (every? m/zero? (map #(sut/bessel-J % 0.0) (range 0.2 100.1 0.1))))
  (t/is (every? m/zero? (map #(sut/bessel-J % 0.0) (range -100 0))))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05]
          nu [2, 4, 6, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-14)))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5, 98.5, 104.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-14))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-13)))
  (doseq [nu [150, 165.2, 200.0, 300.0, 500.0, 1000.0, 5000.2, 10000.0, 50000.0]
          x [0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92,0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,0.995, 0.999, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (base/besselJ xx nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (base/besselJ x nu))) 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-J nu x) (first (rr/r->clj (Bessel/BesselJ x nu))) 1.0e-12)))
  (let [vs (range 0.0 250.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-J % 0.15) vs) (rr/r->clj (base/besselJ 0.15 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 2.1) vs) (rr/r->clj (base/besselJ 2.1 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 42.1) vs) (rr/r->clj (base/besselJ 42.1 vs)) 1.0e-14))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 142.1) vs) (rr/r->clj (base/besselJ 142.1 vs)) 1.0e-14)))
  (let [vs (range -100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-J % 0.15) vs) (rr/r->clj (base/besselJ 0.15 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 2.1) vs) (rr/r->clj (base/besselJ 2.1 vs)) 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 42.1) vs) (rr/r->clj (base/besselJ 42.1 vs)) 1.0e-14))
    (t/is (v/edelta-eq (map #(sut/bessel-J % 142.1) vs) (rr/r->clj (base/besselJ 142.1 vs)) 1.0e-14)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -100 0 5)
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-J nu xx) (first (rr/r->clj (Bessel/BesselJ xx nu))) 1.0e-14 1.0e-14)))
  (t/are [v x] (m/delta-eq (sut/bessel-J v x) (first (rr/r->clj (Bessel/BesselJ x v))) 1.0e-15)
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
  (t/is (m/nan? (sut/bessel-Y0 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-Y0 0.0)))
  (t/is (m/zero? (sut/bessel-Y0 ##Inf)))
  (t/is (m/nan? (sut/bessel-Y0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-Y0 x) (rr/r->clj (base/besselY x 0.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-Y0 x) (rr/r->clj (Bessel/BesselY x 0.0)) 1.0e-14 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-Y0 xl) [ -7.91380268385095e-7
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
  (t/is (m/nan? (sut/bessel-Y1 ##NaN)))
  (t/is (m/neg-inf? (sut/bessel-Y1 0.0)))
  (t/is (m/zero? (sut/bessel-Y1 ##Inf)))
  (t/is (m/nan? (sut/bessel-Y1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-Y1 x) (rr/r->clj (base/besselY x 1.0)) 1.0e-15 1.0e-15))
  (t/is (v/edelta-eq (map sut/bessel-Y1 x) (rr/r->clj (Bessel/BesselY x 1.0)) 1.0e-14 1.0e-14))
  ;; large
  (t/is (v/edelta-eq (map sut/bessel-Y1 xl) [-1.0167125050080249e-7
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
  (t/is (m/nan? (sut/bessel-Y 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-Y -0.5 ##NaN)))
  (doseq [x [0.05, 0.1, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 1.9, 2.5, 3.0, 3.5, 5.0, 10.0]
          nu [0, 1, 2, 4, 6, 10, 15, 20, 25, 30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 125, 150, 175, 200]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-Y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-Y nu xx) (first (rr/r->clj (Bessel/BesselY xx nu))) 1.0e-12 1.0e-12)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-Y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-Y nu xx) (first (rr/r->clj (Bessel/BesselY xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-Y nu x) (first (rr/r->clj (base/besselY x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-Y nu x) (first (rr/r->clj (Bessel/BesselY x nu))) 1.0e-13 1.0e-13)))
  (doseq [nu [150, 165.2, 200.0, 300.0, 500.0, 1000.0, 5000.2, 10000.0 20000.0]
          x [0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92,0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,0.995, 0.999, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-Y nu xx) (first (rr/r->clj (base/besselY xx nu))) 1.0e-10 1.0e-10)))
  (let [vs (range 0.0 100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-Y % 0.15) vs) (rr/r->clj (base/besselY 0.15 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-Y % 2.1) vs) (rr/r->clj (base/besselY 2.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-Y % 42.1) vs) (rr/r->clj (base/besselY 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-Y % 142.1) vs) (rr/r->clj (base/besselY 142.1 vs)) 1.0e-14 1.0e-14)))
  (doseq [v (range -100.0 0.0 0.25)]
    (t/is (m/delta-eq (sut/bessel-Y v 0.15) (first (rr/r->clj (Bessel/BesselY 0.15 v))) 1.0e-6 1.0e-6))
    (t/is (m/delta-eq (sut/bessel-Y v 2.1) (first (rr/r->clj (Bessel/BesselY 2.1 v))) 1.0e-8 1.0e-8))
    (t/is (m/delta-eq (sut/bessel-Y v 42.1) (first (rr/r->clj (Bessel/BesselY 42.1 v))) 1.0e-9 1.0e-9))
    (t/is (m/delta-eq (sut/bessel-Y v 142.1) (first (rr/r->clj (Bessel/BesselY 142.1 v))) 1.0e-12 1.0e-12)))
  (t/are [v x] (m/delta-eq (sut/bessel-Y v x) (first (rr/r->clj (Bessel/BesselY x v))) 1.0e-14 1.0e-14)
    -6.2 18.6
    -8.0 23.2
    -7.0 23.2
    -6.0 23.2
    -0.1 2.2))

;;

(t/deftest bessel-I0
  (t/is (m/nan? (sut/bessel-I0 ##NaN)))
  (t/is (m/one? (sut/bessel-I0 0.0)))
  (t/is (m/nan? (sut/bessel-I0 ##Inf)))
  (t/is (m/nan? (sut/bessel-I0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-I0 x) (rr/r->clj (base/besselI x 0.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-I0 x) (rr/r->clj (Bessel/BesselI x 0.0)) 1.0e-14 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-I0 -x) (rr/r->clj (Bessel/BesselI -x 0.0)) 1.0e-14 1.0e-14)))

(t/deftest bessel-I1
  (t/is (m/nan? (sut/bessel-I1 ##NaN)))
  (t/is (m/zero? (sut/bessel-I1 0.0)))
  (t/is (m/nan? (sut/bessel-I1 ##Inf)))
  (t/is (m/nan? (sut/bessel-I1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-I1 x) (rr/r->clj (base/besselI x 1.0)) 1.0e-14 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-I1 x) (rr/r->clj (Bessel/BesselI x 1.0)) 1.0e-14 1.0e-14))
  ;; negative
  (t/is (v/edelta-eq (map sut/bessel-I1 -x) (rr/r->clj (Bessel/BesselI -x 1.0)) 1.0e-14 1.0e-14)))

;;

(t/deftest bessel-I
  (t/is (m/nan? (sut/bessel-I 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-I -0.5 ##NaN)))
  (doseq [x [0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 1.0, 1.001, 1.01, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 1.9, 2.5, 3.0, 3.5, 4.0]
          nu [0.01,0.1, 0.5, 0.8, 1, 1.23, 2,2.56, 4,5.23, 6,9.2, 10,12.89, 15, 19.1, 20, 25, 30, 33.123, 40, 45, 50, 51.5, 55, 60, 65, 70, 72.34, 75, 80, 82.1, 85, 88.76, 90, 92.334, 95, 99.87,100, 110, 125, 145.123, 150, 160.789]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-I nu xx) (first (rr/r->clj (base/besselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-I nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-I nu x) (first (rr/r->clj (base/besselI x nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-I nu x) (first (rr/r->clj (Bessel/BesselI x nu))) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435, 80.5]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-I nu xx) (first (rr/r->clj (base/besselI xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-I nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-I nu x) (first (rr/r->clj (base/besselI x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-I nu x) (first (rr/r->clj (Bessel/BesselI x nu))) 1.0e-13 1.0e-13)))
  (let [vs (range 0.0 250.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-I % 0.15) vs) (rr/r->clj (base/besselI 0.15 vs)) 1.0e-15 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 2.1) vs) (rr/r->clj (base/besselI 2.1 vs)) 1.0e-15 1.0e-15))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 42.1) vs) (rr/r->clj (base/besselI 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 142.1) vs) (rr/r->clj (base/besselI 142.1 vs)) 1.0e-12 1.0e-12)))
  (let [vs (range -100.0 0.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-I % 0.15) vs) (rr/r->clj (base/besselI 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 2.1) vs) (rr/r->clj (base/besselI 2.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 42.1) vs) (rr/r->clj (base/besselI 42.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-I % 142.1) vs) (rr/r->clj (base/besselI 142.1 vs)) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -100 0 5)
          :let [xx (m/* x nu)
                axx (m/abs xx)]]
    (t/is (m/delta-eq (sut/bessel-I nu xx) (first (rr/r->clj (Bessel/BesselI xx nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-I nu axx) (first (rr/r->clj (Bessel/BesselI axx nu))) 1.0e-12 1.0e-12)))
  (t/are [v x] (m/delta-eq (sut/bessel-I v x) (first (rr/r->clj (Bessel/BesselI x v))) 1.0e-14 1.0e-14)
    12.0 3.2
    13.0 -1.0
    -8.0 4.2
    12.3 8.2
    -12.3 8.2
    -14.0 -9.9))

;;

(t/deftest bessel-K0
  (t/is (m/nan? (sut/bessel-K0 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-K0 0.0)))
  (t/is (m/zero? (sut/bessel-K0 ##Inf)))
  (t/is (m/nan? (sut/bessel-K0 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-K0 x) (rr/r->clj (base/besselK x 0.0)) 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-K0 x) (rr/r->clj (Bessel/BesselK x 0.0)) 1.0e-14)))

(t/deftest bessel-K1
  (t/is (m/nan? (sut/bessel-K1 ##NaN)))
  (t/is (m/pos-inf? (sut/bessel-K1 0.0)))
  (t/is (m/zero? (sut/bessel-K1 ##Inf)))
  (t/is (m/nan? (sut/bessel-K1 ##-Inf)))
  ;; positive
  (t/is (v/edelta-eq (map sut/bessel-K1 x) (rr/r->clj (base/besselK x 1.0)) 1.0e-14))
  (t/is (v/edelta-eq (map sut/bessel-K1 x) (rr/r->clj (Bessel/BesselK x 1.0)) 1.0e-14)))

(t/deftest bessel-K
  (t/is (m/nan? (sut/bessel-K 0.5 ##NaN)))
  (t/is (m/nan? (sut/bessel-K -0.5 ##NaN)))
  (doseq [nu (range -36.0 82.0 0.87654)
          x (rest (m/slice-range 0 30 31))]
    (t/is (m/delta-eq (sut/bessel-K nu x) (first (rr/r->clj (base/besselK x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-K nu x) (first (rr/r->clj (Bessel/BesselK x nu))) 1.0e-12 1.0e-12)))
  (doseq [x [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu [0.1, 0.4567, 0.8123, 1.5, 2.5, 4.1234, 6.8, 12.3, 18.9, 28.2345, 38.1235, 51.23, 72.23435]
          :let [xx (m/* x nu)]]
    (t/is (m/delta-eq (sut/bessel-K nu xx) (first (rr/r->clj (base/besselK xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-K nu xx) (first (rr/r->clj (Bessel/BesselK xx nu))) 1.0e-13 1.0e-13))
    (t/is (m/delta-eq (sut/bessel-K nu x) (first (rr/r->clj (base/besselK x nu))) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq (sut/bessel-K nu x) (first (rr/r->clj (Bessel/BesselK x nu))) 1.0e-13 1.0e-13)))
  (let [vs (range 0.0 100.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-K % 0.15) vs) (rr/r->clj (base/besselK 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 2.1) vs) (rr/r->clj (base/besselK 2.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 42.1) vs) (rr/r->clj (base/besselK 42.1 vs)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 142.1) vs) (rr/r->clj (base/besselK 142.1 vs)) 1.0e-15 1.0e-15)))
  (let [vs (range -100.0 0.0 0.25)]
    (t/is (v/edelta-eq (map #(sut/bessel-K % 0.15) vs) (rr/r->clj (base/besselK 0.15 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 2.1) vs) (rr/r->clj (base/besselK 2.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 42.1) vs) (rr/r->clj (base/besselK 42.1 vs)) 1.0e-12 1.0e-12))
    (t/is (v/edelta-eq (map #(sut/bessel-K % 142.1) vs) (rr/r->clj (base/besselK 142.1 vs)) 1.0e-13 1.0e-13)))
  (doseq [x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          nu (range -50 0 0.25)
          :let [xx (m/abs (* nu x))]]
    (t/is (m/delta-eq (sut/bessel-K nu xx) (first (rr/r->clj (Bessel/BesselK xx nu))) 1.0e-13 1.0e-13)))
  (t/are [v x] (m/delta-eq (sut/bessel-K v x) (first (rr/r->clj (Bessel/BesselK x v))) 1.0e-14 1.0e-14)
    12.0 3.2
    -8.0 4.2
    12.3 8.2
    -12.3 8.2))

(t/deftest bessel-K-half
  (doseq [o (range 1 100 2)
          x [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5,0.55,  0.6,0.65,  0.7, 0.75, 0.8, 0.85, 0.9, 0.92, 0.95, 0.97, 0.99, 1.0, 1.01, 1.05, 1.08, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 4.5, 4.99, 5.1]
          :let [xx (* o x)
                oh (* o 0.5)]]
    (t/is (m/delta-eq (sut/bessel-K oh xx) (sut/bessel-K-half-odd o xx) 1.0e-13 1.0e-13))))

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
  (t/is (m/nan? (sut/kummers-M -1 -1 2))) ;; example 4
  (t/testing "neg b and pos a is nan" (t/is (m/nan? (sut/kummers-M 10 -10 2))))
  (t/testing "neg b and neg a and a<b is nan"(t/is (m/nan? (sut/kummers-M -20 -10 2))))) 

(t/deftest whittaker-m
  (t/is (m/delta-eq 1.10622 (sut/whittaker-M 0 -0.4 1) 1.0e-5)))

(t/deftest besselk
  (t/are [order ress] (v/delta-eq ress (mapv (partial sut/bessel-K-half-odd order) [0.5 1 1.33 2.5 5]))
    1 [1.075047603 0.461068504 0.287423621 0.065065943 0.003776613]
    3 [3.225142810 0.922137009 0.503531608 0.091092320 0.004531936]
    5 [20.425904466  3.227479531  1.423209202  0.174376728  0.006495775]
    7 [207.48418748  17.05953466   5.85394214   0.43984578   0.01102771]
    9 [2.925204529e+03 1.226442222e+02 3.223343101e+01 1.405944900e+00 2.193457048e-02]))

;;

(t/deftest airy
  (let [rsmall (range -20.00001 21.00001 0.01)
        rlargea (range -12345.5 5678.5 9.123)
        rlargeb (range -12345.5 101.5 9.123)
        rlargea' (range -3234.5 5678.5 9.123)]
    (t/is (v/edelta-eq (map sut/airy-Ai rsmall) (rr/r->clj (Bessel/AiryA rsmall)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map sut/airy-Bi rsmall) (rr/r->clj (Bessel/AiryB rsmall)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map sut/airy-Ai' rsmall) (rr/r->clj (Bessel/AiryA rsmall :deriv 1)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map sut/airy-Bi' rsmall) (rr/r->clj (Bessel/AiryB rsmall :deriv 1)) 1.0e-13 1.0e-13))
    (t/is (v/edelta-eq (map sut/airy-Ai rlargea) (rr/r->clj (Bessel/AiryA rlargea)) 1.0e-11 1.0e-11))
    (t/is (v/edelta-eq (map sut/airy-Bi rlargeb) (rr/r->clj (Bessel/AiryB rlargeb)) 1.0e-11 1.0e-11))
    (t/is (v/edelta-eq (map sut/airy-Ai' rlargea') (rr/r->clj (Bessel/AiryA rlargea' :deriv 1)) 1.0e-10 1.0e-10))
    (t/is (v/edelta-eq (map sut/airy-Bi' rlargeb) (rr/r->clj (Bessel/AiryB rlargeb :deriv 1)) 1.0e-9 1.0e-9))))

;;

(t/deftest lambert-w
  (doseq [x (repeatedly 1000 #(r/drand -1.0 500.0))]
    (t/is (v/delta-eq (sut/lambert-W (m/* x (m/exp x))) x)))
  (doseq [x (repeatedly 1000 #(r/drand -500.0 -1.0))]
    (t/is (v/delta-eq (sut/lambert-W-1(m/* x (m/exp x))) x))))

;;

;; owens_t

(def hvec [0.0625, 6.5, 7.0, 4.78125, 2.0, 1.0, 0.0625, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.125, 0.125, 0.125, 0.125, 0.0078125
         , 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0078125, 0.0625, 0.5, 0.9, 2.5, 7.33, 0.6, 1.6, 2.33, 2.33])
(def avec [0.25, 0.4375, 0.96875, 0.0625, 0.5, 0.9999975, 0.999999125, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 0.5, 1, 2, 3, 10, 100
         , 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999, 0.999999999999999
         , 0.99999])
(def cvec [0.0389119302347013668966224771378, 2.00057730485083154100907167685e-11, 6.399062719389853083219914429e-13
         , 1.06329748046874638058307112826e-7, 0.00862507798552150713113488319155, 0.0667418089782285927715589822405
         , 0.1246894855262192
         , 0.04306469112078537, 0.06674188216570097, 0.0784681869930841, 0.0792995047488726, 0.06448860284750375, 0.1066710629614485
         , 0.1415806036539784, 0.1510840430760184, 0.07134663382271778, 0.1201285306350883, 0.1666128410939293, 0.1847501847929859
         , 0.07317273327500386, 0.1237630544953746, 0.1737438887583106, 0.1951190307092811, 0.07378938035365545
         , 0.1249951430754052, 0.1761984774738108, 0.1987772386442824, 0.2340886964802671, 0.2479460829231492
         , 0.1246895548850743676554299881345328280176736760739893903915691894
         , 0.1066710629614484543187382775527753945264849005582264731161129477
         , 0.0750909978020473015080760056431386286348318447478899039422181015
         , 0.0030855526911589942124216949767707430484322201889568086810922629
         , 5.7538182971139187466647478665179637676531179007295252996453e-14, 0.0995191725772188724714794470740785702586033387786949658229016920
         , 0.0258981646643923680014142514442989928165349517076730515952020227
         , 0.0049025023268168675126146823752680242063832053551244071400100690
         , 0.0049024988349089450612896251009169062698683918433614542387524648])

(t/deftest owens_t
  (doseq [[h a res] (map vector hvec avec cvec)]
    (t/is (m/delta-eq res (sut/owens-t h a)))))
