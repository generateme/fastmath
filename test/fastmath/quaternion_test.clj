;; https://github.com/MartinWeigel/Quaternion/blob/master/TestQuaternion.c
(ns fastmath.quaternion-test
  (:require [fastmath.quaternion :as sut]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.core :as m]))

(t/deftest rotation-quaternion
  (t/is (v/delta-eq (sut/rotation-quaternion m/HALF_PI (v/vec3 1 0 0))
                    (sut/quaternion 0.7071 0.7071 0 0) 1.0e-4)))

(t/deftest multiplication
  (t/are [a2 b2 c2 d2 ra rb rc rd]
      (v/delta-eq (sut/mult (sut/quaternion 1 0 1 0)
                            (sut/quaternion a2 b2 c2 d2))
                  (sut/quaternion ra rb rc rd) 1.0e-4)
    1 0.5 0.5 0.75 0.5 1.25 1.5 0.25
    1 0 1 0 0 0 2 0
    2 1 0.1 0.1 1.9 1.1 2.1 -0.9))

(t/deftest rotate
  (t/are [x y z a b c d rx ry rz]
      (v/delta-eq (v/vec3 rx ry rz)
                  (sut/rotate (v/vec3 x y z) (sut/quaternion a b c d)) 1.0e-4)
    5.1 6.8 -5.3 1 0 0 0 5.1 6.8 -5.3
    1 0 0 0.5 0.5 0.5 0.5 0 1 0
    1 0 0 0.6532815, -0.270598, 0.270598, 0.6532815 0 0.7071 -0.7071
    1 0 0 0.7071, 0, 0, 0.7071 0 1 0))

(t/deftest euler

  (t/are [r p y a b c d] (v/delta-eq (sut/from-euler (v/fmap (v/vec3 r p y) m/radians))
                                     (sut/quaternion a b c d) 1.0e-4)
    165 63 122 0.5070333, 0.3501829, 0.7724199, -0.1538071
    90 0 90 0.5 0.5 0.5 0.5
    90 0 180 0 0 0.7071 0.7071
    90 0 0 0.7071 0.7071 0 0)
  
  (t/are [r p y a b c d] (v/delta-eq (sut/to-euler (sut/quaternion a b c d))
                                     (v/fmap (v/vec3 r p y) m/radians) 1.0e-4)
    165 63 122 0.5070333, 0.3501829, 0.7724199, -0.1538071
    90 0 90 0.5 0.5 0.5 0.5
    90 0 180 0 0 0.7071 0.7071
    90 0 0 0.7071 0.7071 0 0))

(t/deftest slerp
  (let [q1 (sut/quaternion 0.6532815, -0.270598, 0.270598, 0.6532815)
        q2 (sut/quaternion 0.5 0.5 0.5 0.5)]
    (t/is (= q1 (sut/slerp q1 q2 0.0)))
    (t/is (= q2 (sut/slerp q1 q2 1.0)))
    (t/is (v/delta-eq (sut/slerp q1 q2 0.62)
                      (sut/quaternion 0.6119266025696755
                                      0.22069444274723088
                                      0.4498729015909088
                                      0.6119266025696755) 1.0e-4))))
