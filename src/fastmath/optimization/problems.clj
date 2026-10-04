(ns fastmath.optimization.problems
  "A collection of test functions"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; univariate
;; https://infinity77.net/global_optimization/test_functions_1d.html

(defn problem02-bounds [] [[2.7 7.5]])

(defn problem02
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* 3.333333333333333 x))))

(defn dproblem02
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* 3.333333333333333 (m/cos (m/* 3.333333333333333 x))))])

;;

(defn problem03-bounds [] [[-10.0 10.0]])

(defn problem03
  ^double [^double x]
  (m/- (m/+ (m/sin (m/+ (m/* 2.0 x) 1.0))
            (m/* 2.0 (m/sin (m/+ (m/* 3.0 x) 2.0)))
            (m/* 3.0 (m/sin (m/+ (m/* 4.0 x) 3.0)))
            (m/* 4.0 (m/sin (m/+ (m/* 5.0 x) 4.0)))
            (m/* 5.0 (m/sin (m/+ (m/* 6.0 x) 5.0)))
            (m/* 6.0 (m/sin (m/+ (m/* 7.0 x) 6.0))))))

(defn dproblem03
  [[^double x]]
  [(m/- (m/+ (m/* 2.0 (m/cos (m/+ (m/* 2.0 x) 1.0)))
             (m/* 6.0 (m/cos (m/+ (m/* 3.0 x) 2.0)))
             (m/* 12.0 (m/cos (m/+ (m/* 4.0 x) 3.0)))
             (m/* 20.0 (m/cos (m/+ (m/* 5.0 x) 4.0)))
             (m/* 30.0 (m/cos (m/+ (m/* 6.0 x) 5.0)))
             (m/* 42.0 (m/cos (m/+ (m/* 7.0 x) 6.0)))))])

;;

(defn problem04-bounds [] [[1.9 3.9]])

(defn problem04
  ^double [^double x]
  (m/- (m/* (m/+ (m/* 16.0 x x)
                 (m/* -24.0 x) 5.0)
            (m/exp (m/- x)))))

(defn dproblem04
  [[^double x]]
  [(m/* (m/+ (m/* 16.0 x x)
             (m/* -56.0 x) 29.0)
        (m/exp (m/- x)))])

;;

(defn problem05-bounds [] [[0.0 1.2]])

(defn problem05
  ^double [^double x]
  (m/- (m/* (m/- 1.4 (m/* 3.0 x))
            (m/sin (m/* 18.0 x)))))

(defn dproblem05
  [[^double x]]
  [(let [x18 (m/* 18.0 x)]
     (m/- (m/* 3.0 (m/sin x18))
          (m/* 18.0 (m/- 1.4 (m/* 3.0 x)) (m/cos x18))))])

;;

(defn problem06-bounds [] [[-10.0 10.0]])

(defn problem06
  ^double [^double x]
  (m/- (m/* (m/+ x (m/sin x))
            (m/exp (m/- (m/* x x))))))

(defn dproblem06
  [[^double x]]
  [(m/* (m/- (m/+ (m/* 2.0 x x) (m/* 2.0 x (m/sin x))) 1.0 (m/cos x))
        (m/exp (m/- (m/* x x))))])

;;

(defn problem07-bounds [] [[2.7 7.5]])

(defn problem07
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* 3.333333333333333 x))
       (m/log x)
       (m/* -0.84 x)
       3.0))

(defn dproblem07
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* 3.333333333333333 (m/cos (m/* 3.333333333333333 x)))
        (m// x)
        -0.84)])

;;

(defn problem08-bounds [] [[-10.0 10.0]])

(defn problem08
  ^double [^double x]
  (m/- (m/+ (m/cos (m/+ (m/* 2.0 x) 1.0))
            (m/* 2.0 (m/cos (m/+ (m/* 3.0 x) 2.0)))
            (m/* 3.0 (m/cos (m/+ (m/* 4.0 x) 3.0)))
            (m/* 4.0 (m/cos (m/+ (m/* 5.0 x) 4.0)))
            (m/* 5.0 (m/cos (m/+ (m/* 6.0 x) 5.0)))
            (m/* 6.0 (m/cos (m/+ (m/* 7.0 x) 6.0))))))

(defn dproblem08
  [[^double x]]
  [(m/+ (m/* 2.0 (m/sin (m/+ (m/* 2.0 x) 1.0)))
        (m/* 6.0 (m/sin (m/+ (m/* 3.0 x) 2.0)))
        (m/* 12.0 (m/sin (m/+ (m/* 4.0 x) 3.0)))
        (m/* 20.0 (m/sin (m/+ (m/* 5.0 x) 4.0)))
        (m/* 30.0 (m/sin (m/+ (m/* 6.0 x) 5.0)))
        (m/* 42.0 (m/sin (m/+ (m/* 7.0 x) 6.0))))])

;; 

(defn problem09-bounds [] [[3.1 20.4]])

(defn problem09
  ^double [^double x]
  (m/+ (m/sin x)
       (m/sin (m/* m/TWO_THIRDS x))))

(defn dproblem09
  [[^double x]]
  [(m/+ (m/cos x)
        (m/* m/TWO_THIRDS (m/cos (m/* m/TWO_THIRDS x))))])

;;

(defn problem10-bounds [] [0.0 10.0])

(defn problem10
  ^double [^double x]
  (m/- (m/* x (m/sin x))))

(defn dproblem10
  [[^double x]]
  [(m/- (m/+ (m/sin x) (m/* x (m/cos x))))])

;;

(defn problem11-bounds [] [[m/-HALF_PI m/TWO_PI]])

(defn problem11
  ^double [^double x]
  (m/+ (m/* 2.0 (m/cos x))
       (m/cos (m/* 2.0 x))))

(defn dproblem11
  [[^double x]]
  [(m/- (m/+ (m/* 2.0 (m/sin x))
             (m/* 2.0 (m/sin (m/* 2.0 x)))))])

;;

(defn problem12-bounds [] [[0.0 m/TWO_PI]])

(defn problem12
  ^double [^double x]
  (let [sx (m/sin x)
        cx (m/cos x)]
    (m/+ (m/* sx sx sx) (m/* cx cx cx))))

(defn dproblem12
  [[^double x]]
  [(let [sx (m/sin x)
         cx (m/cos x)]
     (m/- (m/* 3.0 sx sx cx) (m/* 3.0 cx cx sx)))])

;;

(defn problem13-bounds [] [[0.001 0.99]])

(defn problem13
  ^double [^double x]
  (m/- (m/+ (m/pow x m/TWO_THIRDS)
            (m/pow (m/- 1.0 (m/* x x)) m/THIRD))))

(defn dproblem13
  [[^double x]]
  [(m/+ (m/* (m/- m/TWO_THIRDS) (m/pow x (m/- m/THIRD)))
        (m/* m/TWO_THIRDS x (m/pow (m/- 1.0 (m/* x x)) (m/- m/TWO_THIRDS))))])

;;

(defn problem14-bounds [] [[0.0 4.0]])

(defn problem14
  ^double [^double x]
  (m/- (m/* (m/exp (m/- x)) (m/sin (m/* m/TWO_PI x)))))

(defn dproblem14
  [[^double x]]
  [(let [px (m/* m/TWO_PI x)]
     (m/* (m/exp (m/- x))
          (m/- (m/sin px) (m/* m/TWO_PI (m/cos px)))))])

;;

(defn problem15-bounds [] [[-5.0 5.0]])

(defn problem15
  ^double [^double x]
  (let [x2 (m/* x x)]
    (m// (m/+ x2 (m/* -5.0 x) 6.0)
         (m/inc x2))))

(defn dproblem15
  [[^double x]]
  [(let [x2 (m/* x x)]
     (m// (m/+ (m/* 5.0 x2) (m/* -10.0 x) -5.0)
          (m/sq (m/inc x2))))])

;;

(defn problem18-bounds [] [[0.0 6.0]])

(defn problem18
  ^double [^double x]
  (if (m/<= x 3.0)
    (m/sq (m/- x 2.0))
    (m/+ (m/* 2.0 (m/log (m/- x 2.0))) 1.0)))

(defn dproblem18
  [[^double x]]
  [(if (m/<= x 3.0)
     (m/* 2.0 (m/- x 2.0))
     (m// 2.0 (m/- x 2.0)))])

;;

(defn problem20-bounds [] [[-10.0 10.0]])

(defn problem20
  ^double [^double x]
  (m/- (m/* (m/- x (m/sin x))
            (m/exp (m/- (m/* x x))))))

(defn dproblem20
  [[^double x]]
  (let [x2 (m/* x x)]
    [(m/* (m/+ (m/- (m/* 2.0 x2) (m/* 2.0 x (m/sin x))) -1.0 (m/cos x))
          (m/exp (m/- x2)))]))

;;

(defn problem21-bounds [] [[0.0 10.0]])

(defn problem21
  ^double [^double x]
  (m/+ (m/* x (m/sin x))
       (m/* x (m/cos (m/* 2.0 x)))))

(defn dproblem21
  [[^double x]]
  (let [x2 (m/* 2.0 x)]
    [(m/+ (m/sin x)
          (m/* x (m/cos x))
          (m/cos x2)
          (m/* -2.0 x (m/sin x2)))]))

;;

(defn problem22-bounds [] [[0.0 20.0]])

(defn problem22
  ^double [^double x]
  (m/- (m/exp (m/* -3.0 x))
       (m/fpow (m/sin x) 3)))

(defn dproblem22
  [[^double x]]
  (let [sx (m/sin x)]
    [(m/- (m/* -3.0 (m/exp (m/* -3.0 x)))
          (m/* 3.0 sx sx (m/cos x)))]))

;; https://www.sfu.ca/~ssurjano/optimization.html

(defn auckley-bounds [^long N] (repeat N [-32768.0 32768.0]))

(defn ->auckley
  ([] (->auckley {:a 20.0 :b 0.2 :c m/TWO_PI}))
  ([{:keys [^double a ^double b ^double c]}]
   (let [ae (m/+ a m/E)
         -b (m/- b)]
     (fn ^double [vs]
       (m/- ae
            (m/* a (m/exp (m/* -b (v/average (v/magsq vs)))))
            (m/exp (v/average (v/cos (v/mult vs c)))))))))

;;

(defn bukin-no-6-bounds [] [[-15.0 -5.0] [-3.0 3.0]])

(defn bukin-no-6
  ^double [[^double x1 ^double x2]]
  (m/+ (m/* 100.0 (m/sqrt (m/abs (m/- x2 (m/* 0.01 x1 x1)))))
       (m/* 0.01 (m/abs (m/+ x1 10.0)))))

;;

(defn cross-in-tray-bounds [] [[-10.0 10.0] [-10.0 10.0]])

(defn cross-in-tray
  ^double [[^double x1 ^double x2 :as v]]
  (m/* -0.0001 (m/pow (m/inc (m/abs (m/* (m/sin x1)
                                         (m/sin x2)
                                         (m/exp (m/abs (m/- 100.0 (m// (v/mag v) m/PI))))))) 0.1)))

;;

(defn drop-wave-bounds [] [[-5.12 5.12] [-5.12 5.12]])

(defn drop-wave
  ^double [v]
  (let [ms (v/magsq v)]
    (m/- (m// (m/inc (m/cos (m/* 12.0 (m/sqrt ms))))
              (m/+ (m/* 0.5 ms) 2.0)))))

;;

(defn egg-holder-bounds [] [[-512.0 512.0] [-512.0 512.0]])

(defn egg-holder
  ^double [[^double x1 ^double x2]]
  (let [x47 (m/+ x2 47.0)]
    (m/- (m/* -1.0 x47 (m/sin (m/sqrt (m/abs (m/+ x47 (m/* 0.5 x1))))))
         (m/* x1 (m/sin (m/sqrt (m/abs (m/- x1 x47))))))))

;;

;;

(defn rosenbrock-bounds [^long N] (repeat N [-5.0 10.0]))

(defn rosenbrock
  ^double [v]
  (->> (partition 2 1 v)
       (map (fn [[^double x1 ^double x2]]
              (m/+ (m/* 100.0 (m/sq (m/- x2 (m/* x1 x1))))
                   (m/sq (m/dec x1)))))
       (v/sum)))

;;

;;

(defn himmelblau-bounds [] [[-5.0 5.0] [-5.0 5.0]])

(defn himmelblau
  ^double [[^double x1 ^double x2]]
  (m/+ (m/sq (m/+ (m/* x1 x1) x2 -11.0))
       (m/sq (m/+ x1 (m/* x2 x2) -7.0))))

;;

(defn beale-bounds [] [[-4.5 4.5] [4.5 4.5]])

(defn beale
  ^double [[^double x1 ^double x2]]
  (let [x1x2 (m/* x1 x2)
        x1x22 (m/* x1x2 x2)]
    (m/+ (m/sq (m/+ (m/- 1.5 x1) x1x2))
         (m/sq (m/+ (m/- 2.25 x1) x1x22))
         (m/sq (m/+ (m/- 2.625 x1) (m/* x1x22 x2))))))

;;

(defn sphere-bounds [^long N] (repeat N [-5.0 6.0]))

(defn sphere ^double [v] (v/magsq v))

(defn sphere-gradient [v] (v/mult v 2.0))


