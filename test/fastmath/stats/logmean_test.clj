(ns fastmath.stats.logmean-test
  (:require [fastmath.stats.logmean :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]))

(t/deftest logmean
  (t/testing "single value"
    (t/is (m/delta-eq 3.22 (sut/logmean-integral [3.22])))
    (t/is (m/delta-eq 3.22 (sut/logmean-mean-value [3.22]))))
  (t/testing "two values"
    (t/is (m/delta-eq (sut/logmean2 1 10) (sut/logmean2 [1 10])))
    (t/is (m/delta-eq (sut/logmean-mean-value [1 10]) (sut/logmean2 1 10)))
    (t/is (m/delta-eq (sut/logmean-integral [1 10]) (sut/logmean2 1 10)))
    (t/is (m/delta-eq (sut/logmean-mean-value [0.1 100.5]) (sut/logmean2 0.1 100.5)))
    (t/is (m/delta-eq (sut/logmean-integral [0.1 100.5]) (sut/logmean2 0.1 100.5)))
    (t/is (m/delta-eq (sut/logmean-mean-value [0.1 0.10001]) (sut/logmean2 0.1 0.10001)))
    (t/is (m/delta-eq (sut/logmean-integral [0.1 0.10001]) (sut/logmean2 0.1 0.10001)))
    (t/testing "same values"
      (t/is (m/delta-eq 3.22 (sut/logmean2 3.22 3.22)))
      (t/is (m/delta-eq 3.22 (sut/logmean-integral [3.22 3.22])))))
  (t/testing "three values"
    (t/is (m/delta-eq (sut/logmean3-mean-value 1 10 100) (sut/logmean3-mean-value [1 10 100])))
    (t/is (m/delta-eq (sut/logmean3-integral 1 10 100) (sut/logmean3-integral [1 10 100])))
    (t/is (m/delta-eq (sut/logmean-mean-value [1 10 100]) (sut/logmean3-mean-value 1 10 100)))
    (t/is (m/delta-eq (sut/logmean-integral [1 10 100]) (sut/logmean3-integral 1 10 100)))
    (t/is (m/delta-eq (sut/logmean-mean-value [0.1 100.5 1.0e6]) (sut/logmean3-mean-value 0.1 100.5 1.0e6)))
    (t/is (m/delta-eq (sut/logmean-integral [0.1 100.5 1.0e6]) (sut/logmean3-integral 0.1 100.5 1.0e6)))
    (t/is (m/delta-eq (sut/logmean-mean-value [0.1 0.10001 0.10002]) (sut/logmean3-mean-value 0.1 0.10001 0.10002)))
    (t/is (m/delta-eq (sut/logmean-integral [0.1 0.10001 0.10002]) (sut/logmean3-integral 0.1 0.10001 0.10002)))
    (t/testing "same values"
      (t/is (m/delta-eq 3.22 (sut/logmean3-mean-value 3.22 3.22 3.22)))
      (t/is (m/delta-eq 3.22 (sut/logmean3-integral 3.22 3.22 3.22)))
      (t/is (m/delta-eq 3.22 (sut/logmean-integral [3.22 3.22 3.22])))))
  (t/testing "paper example"
    (t/is (m/delta-eq 73578.65538616560 (sut/logmean-integral (range 1 200001))))))
