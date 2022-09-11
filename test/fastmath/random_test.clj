(ns fastmath.random-test
  (:require [fastmath.random :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]))

(t/deftest zaga
  (t/are [mu sigma nu lower-tail? d p q v]
      (let [dist (sut/distribution :zaga {:mu mu
                                          :sigma sigma
                                          :nu nu
                                          :lower-tail? lower-tail?})]
        (and (m/approx-eq d (sut/pdf dist v) 6)
             (m/approx-eq p (sut/cdf dist v) 6)
             (m/approx-eq q (sut/icdf dist v) 6)))
    0.1 2.0 0.3 true 0.1169786 0.9669273 0.001805674 0.5
    0.1 2.0 0.3 true 1.063233 0.8205746 0.0 0.1
    0.1 2.0 0.3 true 0.0276922 0.9911922 0.2159042 0.9
    0.1 2.0 0.3 true 0.3 0.3 0.0 0.0))

(t/deftest zaga-mv
  (let [dist (sut/distribution :zaga)]
    (t/is (m/approx-eq 0.9 (sut/mean dist)))
    (t/is (m/approx-eq 0.99 (sut/variance dist))))
  (let [dist (sut/distribution :zaga {:mu 0.1 :sigma 2.0 :nu 0.3})]
    (t/is (m/approx-eq 0.07 (sut/mean dist)))
    (t/is (m/approx-eq 0.0301 (sut/variance dist)))))

(t/deftest zinbi
  (t/are [mu sigma nu vd d vp p vq q]
      (let [dist (sut/distribution :zinbi {:mu mu :sigma sigma :nu nu})]
        (and (m/approx-eq d (sut/pdf dist vd) 6)
             (m/approx-eq p (sut/cdf dist vp) 6)
             (m/approx-eq q (sut/icdf dist vq) 6)))
    1.0 1.0 0.3 3.0 0.04375 3.0 0.95625 0.9 2.0
    2.0 0.000009 0.5 3.0 0.09022352 3.0 0.9285617 0.9 3.0
    2.0 0.000009 0.5 0.0 0.5676676 0.0 0.5676676 0.5 0.0))

(t/deftest zinbi-mv
  (let [dist (sut/distribution :zinbi)]
    (t/is (m/approx-eq 0.7 (sut/mean dist)))
    (t/is (m/approx-eq 1.61 (sut/variance dist))))
  (let [dist (sut/distribution :zinbi {:mu 2.0 :sigma 0.000009 :nu 0.5})]
    (t/is (m/approx-eq 1.0 (sut/mean dist)))
    (t/is (m/approx-eq 2.0 (sut/variance dist)))))

(t/deftest zanbi
  (t/are [mu sigma nu vd d vp p vq q]
      (let [dist (sut/distribution :zanbi {:mu mu :sigma sigma :nu nu})]
        (and (m/approx-eq d (sut/pdf dist vd) 6)
             (m/approx-eq p (sut/cdf dist vp) 6)
             (m/approx-eq q (sut/icdf dist vq) 6)))
    1.0 1.0 0.3 3.0 0.0875 3.0 0.9125 0.99 7.0
    2.0 0.000009 0.5 3.0 0.1043451 3.0 0.9173804 0.99 5.0
    2.0 0.000009 0.5 0.0 0.5 0.0 0.5 0.9 3.0))

(t/deftest zanbi-mv
  (let [dist (sut/distribution :zanbi)]
    (t/is (m/approx-eq 1.4 (sut/mean dist)))
    (t/is (m/approx-eq 2.24 (sut/variance dist))))
  (let [dist (sut/distribution :zanbi {:mu 2.0 :sigma 0.000009 :nu 0.5})]
    (t/is (m/approx-eq 1.156521 (sut/mean dist)))
    (t/is (m/approx-eq 2.132043 (sut/variance dist)))))

(t/deftest zip
  (t/are [mu sigma vd d vp p vq q]
      (let [dist (sut/distribution :zip {:mu mu :sigma sigma})]
        (and (m/approx-eq d (sut/pdf dist vd) 6)
             (m/approx-eq p (sut/cdf dist vp) 6)
             (m/approx-eq q (sut/icdf dist vq) 6)))
    5.0 0.1 3.0 0.1263365 3.0 0.3385233 0.7 6.0
    5.0 0.1 0.0 0.1060642 0.0 0.1060642 0.5 5.0
    1.0 0.2 2.0 0.1471518 2.0 0.9357589 0.99 4.0))

(t/deftest zip-mv
  (let [dist (sut/distribution :zip)]
    (t/is (m/approx-eq 4.5 (sut/mean dist)))
    (t/is (m/approx-eq 6.75 (sut/variance dist))))
  (let [dist (sut/distribution :zip {:mu 1.0 :sigma 0.2})]
    (t/is (m/approx-eq 0.8 (sut/mean dist)))
    (t/is (m/approx-eq 0.96 (sut/variance dist)))))

(t/deftest zip2
  (t/are [mu sigma vd d vp p vq q]
      (let [dist (sut/distribution :zip2 {:mu mu :sigma sigma})]
        (and (m/approx-eq d (sut/pdf dist vd) 6)
             (m/approx-eq p (sut/cdf dist vp) 6)
             (m/approx-eq q (sut/icdf dist vq) 6)))
    5.0 0.1 3.0 0.0994321 3.0 0.2759344 0.7 6.0
    5.0 0.1 0.0 0.1034793 0.0 0.1034793 0.5 5.0
    1.0 0.2 2.0 0.1790654 2.0 0.8947741 0.99 4.0))

(t/deftest zip2-mv
  (let [dist (sut/distribution :zip2)]
    (t/is (m/approx-eq 5.0 (sut/mean dist)))
    (t/is (m/approx-eq 7.777777 (sut/variance dist))))
  (let [dist (sut/distribution :zip2 {:mu 1.0 :sigma 0.2})]
    (t/is (m/approx-eq 1.0(sut/mean dist)))
    (t/is (m/approx-eq 1.25 (sut/variance dist)))))

