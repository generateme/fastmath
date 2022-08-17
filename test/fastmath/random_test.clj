(ns fastmath.random-test
  (:require [fastmath.random :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]))

(t/deftest zaga
  (t/are [mu sigma nu lower-tail? d p q v]
      (let [zaga-dist (sut/distribution :zaga {:mu mu
                                               :sigma sigma
                                               :nu nu
                                               :lower-tail? lower-tail?})]
        (and (m/approx-eq d (sut/pdf zaga-dist v) 6)
             (m/approx-eq p (sut/cdf zaga-dist v) 6)
             (m/approx-eq q (sut/icdf zaga-dist v) 6)))
    0.1 2.0 0.3 true 0.1169786 0.9669273 0.001805674 0.5
    0.1 2.0 0.3 true 1.063233 0.8205746 0.0 0.1
    0.1 2.0 0.3 true 0.0276922 0.9911922 0.2159042 0.9
    0.1 2.0 0.3 true 0.3 0.3 0.0 0.0))
