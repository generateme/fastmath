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
        (and (m/delta-eq d (sut/pdf dist v))
             (m/delta-eq p (sut/cdf dist v))
             (m/delta-eq q (sut/icdf dist v))))
    0.1 2.0 0.3 true 0.1169786 0.9669273 0.001805674 0.5
    0.1 2.0 0.3 true 1.063233 0.8205746 0.0 0.1
    0.1 2.0 0.3 true 0.0276922 0.9911922 0.2159042 0.9
    0.1 2.0 0.3 true 0.3 0.3 0.0 0.0))

(t/deftest zaga-mv
  (let [dist (sut/distribution :zaga)]
    (t/is (m/delta-eq 0.9 (sut/mean dist)))
    (t/is (m/delta-eq 0.99 (sut/variance dist))))
  (let [dist (sut/distribution :zaga {:mu 0.1 :sigma 2.0 :nu 0.3})]
    (t/is (m/delta-eq 0.07 (sut/mean dist)))
    (t/is (m/delta-eq 0.0301 (sut/variance dist)))))

(t/deftest zinbi
  (t/are [mu sigma nu vd d vp p vq q]
      (let [dist (sut/distribution :zinbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    1.0 1.0 0.3 3.0 0.04375 3.0 0.95625 0.9 2.0
    2.0 0.000009 0.5 3.0 0.09022352 3.0 0.9285617 0.9 3.0
    2.0 0.000009 0.5 0.0 0.5676676 0.0 0.5676676 0.5 0.0))

(t/deftest zinbi-mv
  (let [dist (sut/distribution :zinbi)]
    (t/is (m/delta-eq 0.7 (sut/mean dist)))
    (t/is (m/delta-eq 1.61 (sut/variance dist))))
  (let [dist (sut/distribution :zinbi {:mu 2.0 :sigma 0.000009 :nu 0.5})]
    (t/is (m/delta-eq 1.0 (sut/mean dist)))
    (t/is (m/delta-eq 2.0 (sut/variance dist) 1.0e-4))))

(t/deftest zanbi
  (t/are [mu sigma nu vd d vp p vq q]
      (let [dist (sut/distribution :zanbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    1.0 1.0 0.3 3.0 0.0875 3.0 0.9125 0.99 7.0
    2.0 0.000009 0.5 3.0 0.1043451 3.0 0.9173804 0.99 5.0
    2.0 0.000009 0.5 0.0 0.5 0.0 0.5 0.9 3.0))

(t/deftest zanbi-mv
  (let [dist (sut/distribution :zanbi)]
    (t/is (m/delta-eq 1.4 (sut/mean dist)))
    (t/is (m/delta-eq 2.24 (sut/variance dist))))
  (let [dist (sut/distribution :zanbi {:mu 2.0 :sigma 0.000009 :nu 0.5})]
    (t/is (m/delta-eq 1.156521 (sut/mean dist)))
    (t/is (m/delta-eq 2.132043 (sut/variance dist)))))

(t/deftest zip
  (t/are [mu sigma vd d vp p vq q]
      (let [dist (sut/distribution :zip {:mu mu :sigma sigma})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    5.0 0.1 3.0 0.1263365 3.0 0.3385233 0.7 6.0
    5.0 0.1 0.0 0.1060642 0.0 0.1060642 0.5 5.0
    1.0 0.2 2.0 0.1471518 2.0 0.9357589 0.99 4.0))

(t/deftest zip-mv
  (let [dist (sut/distribution :zip)]
    (t/is (m/delta-eq 4.5 (sut/mean dist)))
    (t/is (m/delta-eq 6.75 (sut/variance dist))))
  (let [dist (sut/distribution :zip {:mu 1.0 :sigma 0.2})]
    (t/is (m/delta-eq 0.8 (sut/mean dist)))
    (t/is (m/delta-eq 0.96 (sut/variance dist)))))

(t/deftest zip2
  (t/are [mu sigma vd d vp p vq q]
      (let [dist (sut/distribution :zip2 {:mu mu :sigma sigma})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    5.0 0.1 3.0 0.0994321 3.0 0.2759344 0.7 6.0
    5.0 0.1 0.0 0.1034793 0.0 0.1034793 0.5 5.0
    1.0 0.2 2.0 0.1790654 2.0 0.8947741 0.99 4.0))

(t/deftest zip2-mv
  (let [dist (sut/distribution :zip2)]
    (t/is (m/delta-eq 5.0 (sut/mean dist)))
    (t/is (m/delta-eq 7.777777 (sut/variance dist))))
  (let [dist (sut/distribution :zip2 {:mu 1.0 :sigma 0.2})]
    (t/is (m/delta-eq 1.0 (sut/mean dist)))
    (t/is (m/delta-eq 1.25 (sut/variance dist)))))

;; reference values from R's BiasedUrn package:
;; dFNCHypergeo/pFNCHypergeo/qFNCHypergeo, dWNCHypergeo/pWNCHypergeo/qWNCHypergeo
;; (m1=ns, m2=nf, n=n, odds=omega)

(t/deftest fishers-noncentral-hypergeometric
  (t/are [ns nf n omega vd d vp p vq q]
      (let [dist (sut/distribution :fishers-noncentral-hypergeometric {:ns ns :nf nf :n n :omega omega})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    5 5 5 1.0       2   0.396825396825397  2   0.5                0.4                2
    10 7 8 2.0      5   0.343123905854892  5   0.512157167866174  0.5                5
    6 4 5 0.3       2   0.493575294911239  2   0.658100393214985  0.4                2
    500 500 200 2.0 120 0.0324794208221728 120 0.13791080212904   0.121671091056712  120))

(t/deftest fishers-noncentral-hypergeometric-mv
  (let [dist (sut/distribution :fishers-noncentral-hypergeometric)]
    (t/is (m/delta-eq 2.5 (sut/mean dist)))
    (t/is (m/delta-eq 0.694444444444444 (sut/variance dist))))
  (let [dist (sut/distribution :fishers-noncentral-hypergeometric {:ns 10 :nf 7 :n 8 :omega 2.0})]
    (t/is (m/delta-eq 5.44988815405563 (sut/mean dist)))
    (t/is (m/delta-eq 1.0448045002312 (sut/variance dist))))
  (let [dist (sut/distribution :fishers-noncentral-hypergeometric {:ns 6 :nf 4 :n 5 :omega 0.3})]
    (t/is (m/delta-eq 2.2244615916158 (sut/mean dist)))
    (t/is (m/delta-eq 0.59996825497418 (sut/variance dist)))))

(t/deftest wallenius-noncentral-hypergeometric
  (t/are [ns nf n omega vd d vp p vq q]
      (let [dist (sut/distribution :wallenius-noncentral-hypergeometric {:ns ns :nf nf :n n :omega omega})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))
             (m/delta-eq q (sut/icdf dist vq))))
    5 5 5 1.0       2   0.396825396825397  2   0.5                0.4                2
    10 7 8 2.0      5   0.308500769208805  5   0.431312085410887  0.5                6
    6 4 5 0.3       2   0.510523324225684  2   0.770350792583098  0.4                2
    500 500 200 2.0 120 0.016612516577956  120 0.0585082904208895 0.0502020321319115 120))

(t/deftest wallenius-noncentral-hypergeometric-mv
  (let [dist (sut/distribution :wallenius-noncentral-hypergeometric)]
    (t/is (m/delta-eq 2.5 (sut/mean dist)))
    (t/is (m/delta-eq 0.694444444444444 (sut/variance dist))))
  (let [dist (sut/distribution :wallenius-noncentral-hypergeometric {:ns 10 :nf 7 :n 8 :omega 2.0})]
    (t/is (m/delta-eq 5.64705945107076 (sut/mean dist)))
    (t/is (m/delta-eq 1.02461762013875 (sut/variance dist))))
  (let [dist (sut/distribution :wallenius-noncentral-hypergeometric {:ns 6 :nf 4 :n 5 :omega 0.3})]
    (t/is (m/delta-eq 1.99249407636023 (sut/mean dist)))
    (t/is (m/delta-eq 0.558373230211051 (sut/variance dist)))))

;; reference values for the truncated distribution derived analytically (truncated normal
;; and truncated exponential formulas) and cross-checked with R's dnorm/pnorm/qnorm and
;; numerical integration for the exponential mean/variance

(t/deftest truncated
  (t/testing "no truncation (default), matches the base normal distribution"
    (let [dist (sut/distribution :truncated)
          nd (sut/distribution :normal)]
      (t/is (m/delta-eq (sut/pdf nd 0.5) (sut/pdf dist 0.5)))
      (t/is (m/delta-eq (sut/cdf nd 0.5) (sut/cdf dist 0.5)))
      (t/is (m/delta-eq (sut/mean nd) (sut/mean dist) 1.0e-9))
      (t/is (m/delta-eq (sut/variance nd) (sut/variance dist)))))
  (t/testing "symmetric truncated standard normal on [-1, 1]"
    (let [dist (sut/distribution :truncated {:distr (sut/distribution :normal) :left -1.0 :right 1.0})]
      (t/is (m/delta-eq 0.5843685672568165 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.5157034505719383 (sut/pdf dist 0.5)))
      (t/is (m/delta-eq 0.5 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 1.0 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 0.0 (sut/icdf dist 0.5)))
      (t/is (m/delta-eq 0.0 (sut/mean dist) 1.0e-9))
      (t/is (m/delta-eq 0.291125094772793 (sut/variance dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -2.0)) "outside the left bound")
      (t/is (m/delta-eq 0.0 (sut/pdf dist 2.0)) "outside the right bound")))
  (t/testing "asymmetric truncated normal(mu=2, sd=1.5) on [1, 4]"
    (let [dist (sut/distribution :truncated {:distr (sut/distribution :normal {:mu 2.0 :sd 1.5}) :left 1.0 :right 4.0})]
      (t/is (m/delta-eq 0.405246141837459 (sut/pdf dist 2.0)))
      (t/is (m/delta-eq 0.324495743156983 (sut/pdf dist 3.0)))
      (t/is (m/delta-eq 0.377127654158999 (sut/cdf dist 2.0)))
      (t/is (m/delta-eq 0.754255308317997 (sut/cdf dist 3.0)))
      (t/is (m/delta-eq 2.30529906429291 (sut/icdf dist 0.5)))
      (t/is (m/delta-eq 2.35526166552572 (sut/mean dist)))
      (t/is (m/delta-eq 0.643966213749691 (sut/variance dist)))))
  (t/testing "truncated exponential(mean=1.0) on [0.5, 3.0]"
    (let [dist (sut/distribution :truncated {:distr (sut/distribution :exponential) :left 0.5 :right 3.0})]
      (t/is (m/delta-eq 0.660769961056685 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.243083684016409 (sut/pdf dist 2.0)))
      (t/is (m/delta-eq 0.428655528777167 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 0.846341805817443 (sut/cdf dist 2.0)))
      (t/is (m/delta-eq 1.1142574462674 (sut/icdf dist 0.5)))
      (t/is (m/delta-eq 1.27643627541537 (sut/mean dist)))
      (t/is (m/delta-eq 0.391109949588273 (sut/variance dist)))))
  (t/testing "one-sided truncation, left=0.5 only (right defaults to distr's upper-bound)"
    (let [dist (sut/distribution :truncated {:distr (sut/distribution :normal) :left 0.5})]
      (t/is (m/delta-eq 0.784250517840677 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.485782979320519 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 1.01829551596028 (sut/icdf dist 0.5)))
      (t/is (m/delta-eq 1.14107777036806 (sut/mean dist)))
      (t/is (m/delta-eq 0.268480407155879 (sut/variance dist))))))

