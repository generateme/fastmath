(ns fastmath.random-test
  (:require [fastmath.random :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]))

;; reference values from R's gamlss.dist package: dZAGA/pZAGA/qZAGA

(t/deftest zaga
  (t/testing "registered under both :zero-adjusted-gamma and :zaga keys"
    (let [d1 (sut/distribution :zero-adjusted-gamma {:mu 0.1 :sigma 2.0 :nu 0.3})
          d2 (sut/distribution :zaga {:mu 0.1 :sigma 2.0 :nu 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 0.5) (sut/pdf d2 0.5)))
      (t/is (m/delta-eq (sut/cdf d1 0.5) (sut/cdf d2 0.5)))))
  (t/testing "defaults match R's ZAGA() defaults (mu=1, sigma=1, nu=0.1)"
    (let [dist (sut/distribution :zaga)]
      (t/is (m/delta-eq 0.1 (sut/pdf dist 0)))))
  (t/testing "point mass at zero"
    (let [dist (sut/distribution :zaga {:mu 0.1 :sigma 2.0 :nu 0.3})]
      (t/is (m/delta-eq 0.3 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.3 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/icdf dist 0.0)))))
  (t/are [mu sigma nu vd d vp p]
      (let [dist (sut/distribution :zaga {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 1.0 0.1  0   0.100000000000 0   0.100000000000
    1.0 1.0 0.1  0.5 0.545877593741 0.5 0.454122406259
    1.0 1.0 0.1  1.5 0.200817144134 1.5 0.799182855866
    1.0 1.0 0.1  3   0.044808361531 3   0.955191638469
    0.1 2.0 0.3  0   0.300000000000 0   0.300000000000
    0.1 2.0 0.3  0.1 1.063232753083 0.1 0.820574561312
    0.1 2.0 0.3  0.5 0.116978584401 0.5 0.966927309199
    0.1 2.0 0.3  0.9 0.027692202630 0.9 0.991192246011
    5.0 0.5 0.05 0   0.050000000000 0   0.050000000000
    5.0 0.5 0.05 2   0.104749297443 2   0.124872812868
    5.0 0.5 0.05 5   0.148478779258 5   0.588203385652
    5.0 0.5 0.05 10  0.021755869628 10  0.959738893608
    2.0 3.0 0.4  0   0.400000000000 0   0.400000000000
    2.0 3.0 0.4  0.5 0.091964118685 0.5 0.824321852246
    2.0 3.0 0.4  2   0.024675318648 2   0.890994765299
    2.0 3.0 0.4  6   0.007441236418 6   0.943624795578
    10.0 0.2 0.2 0   0.200000000000 0   0.200000000000
    10.0 0.2 0.2 7   0.055091431932 7   0.242541043083
    10.0 0.2 0.2 10  0.159045902936 10  0.621281225155
    10.0 0.2 0.2 14  0.023208684422 14  0.974100712181))

(t/deftest zaga-icdf
  (t/are [mu sigma nu p vq]
      (let [dist (sut/distribution :zaga {:mu mu :sigma sigma :nu nu})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 1.0 0.1  0.01 0.000000000000
    1.0 1.0 0.1  0.05 0.000000000000
    1.0 1.0 0.1  0.15 0.057158413840
    1.0 1.0 0.1  0.3  0.251314428281
    1.0 1.0 0.1  0.5  0.587786664902
    1.0 1.0 0.1  0.7  1.098612288668
    1.0 1.0 0.1  0.85 1.791759469228
    1.0 1.0 0.1  0.99 4.499809670330
    0.1 2.0 0.3  0.01 0.000000000000
    0.1 2.0 0.3  0.05 0.000000000000
    0.1 2.0 0.3  0.15 0.000000000000
    0.1 2.0 0.3  0.3  0.000000000000
    0.1 2.0 0.3  0.5  0.001805673727
    0.1 2.0 0.3  0.7  0.030576892183
    0.1 2.0 0.3  0.85 0.132010251083
    0.1 2.0 0.3  0.99 0.859778766620
    5.0 0.5 0.05 0.01 0.000000000000
    5.0 0.5 0.05 0.05 0.000000000000
    5.0 0.5 0.05 0.15 2.223001024427
    5.0 0.5 0.05 0.3  3.245169924313
    5.0 0.5 0.05 0.5  4.434974150602
    5.0 0.5 0.05 0.7  5.826839996839
    5.0 0.5 0.05 0.85 7.407982152703
    5.0 0.5 0.05 0.99 12.468841812780
    2.0 3.0 0.4  0.01 0.000000000000
    2.0 3.0 0.4  0.05 0.000000000000
    2.0 3.0 0.4  0.15 0.000000000000
    2.0 3.0 0.4  0.3  0.000000000000
    2.0 3.0 0.4  0.5  0.000001093750
    2.0 3.0 0.4  0.7  0.021551484888
    2.0 3.0 0.4  0.85 0.863767481455
    2.0 3.0 0.4  0.99 23.701010264968
    10.0 0.2 0.2 0.01 0.000000000000
    10.0 0.2 0.2 0.05 0.000000000000
    10.0 0.2 0.2 0.15 0.000000000000
    10.0 0.2 0.2 0.3  7.756955470277
    10.0 0.2 0.2 0.5  9.248126434537
    10.0 0.2 0.2 0.7  10.512893291374
    10.0 0.2 0.2 0.85 11.733804793793
    10.0 0.2 0.2 0.99 15.007783664068))

(t/deftest zaga-mv
  (t/are [mu sigma nu mean-v var-v]
      (let [dist (sut/distribution :zaga {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    1.0 1.0 0.1  0.900000000000 0.990000000000
    0.1 2.0 0.3  0.070000000000 0.030100000000
    5.0 0.5 0.05 4.750000000000 7.125000000000
    2.0 3.0 0.4  1.200000000000 22.560000000000
    10.0 0.2 0.2 8.000000000000 19.200000000000))

;; reference values from R's gamlss.dist package: dZAIG/pZAIG/qZAIG

(t/deftest zaig
  (t/testing "registered under both :zero-adjusted-inverse-gaussian and :zaig keys"
    (let [d1 (sut/distribution :zero-adjusted-inverse-gaussian {:mu 0.1 :sigma 2.0 :nu 0.3})
          d2 (sut/distribution :zaig {:mu 0.1 :sigma 2.0 :nu 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 0.5) (sut/pdf d2 0.5)))
      (t/is (m/delta-eq (sut/cdf d1 0.5) (sut/cdf d2 0.5)))))
  (t/testing "defaults match R's ZAIG() defaults (mu=1, sigma=1, nu=0.1)"
    (let [dist (sut/distribution :zaig)]
      (t/is (m/delta-eq 0.1 (sut/pdf dist 0)))))
  (t/testing "point mass at zero"
    (let [dist (sut/distribution :zaig {:mu 0.1 :sigma 2.0 :nu 0.3})]
      (t/is (m/delta-eq 0.3 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.3 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/icdf dist 0.0)))))
  (t/are [mu sigma nu vd d vp p]
      (let [dist (sut/distribution :zaig {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0  1.0 0.10 0    0.100000000000 0    0.100000000000
    1.0  1.0 0.10 0.5  0.790904321042 0.5  0.428477993356
    1.0  1.0 0.10 1.5  0.179814404820 1.5  0.829691193700
    1.0  1.0 0.10 3    0.035476522173 3    0.957869128669
    0.1  2.0 0.30 0    0.300000000000 0    0.300000000000
    0.1  2.0 0.30 0.1  4.415481913535 0.1  0.731314203032
    0.1  2.0 0.30 0.5  0.007233444874 0.5  0.999510274913
    0.1  2.0 0.30 0.9  0.000022553760 0.9  0.999998379322
    5.0  0.5 0.05 0    0.050000000000 0    0.050000000000
    5.0  0.5 0.05 2    0.186970315032 2    0.350426697595
    5.0  0.5 0.05 5    0.067796716414 5    0.698248051302
    5.0  0.5 0.05 10   0.019624778804 10   0.885570045379
    2.0  3.0 0.40 0    0.400000000000 0    0.400000000000
    2.0  3.0 0.40 0.5  0.212002825875 0.5  0.803400264033
    2.0  3.0 0.40 2    0.028209479177 2    0.913675866580
    2.0  3.0 0.40 6    0.005231523796 6    0.960883065951
    10.0 0.2 0.20 0    0.200000000000 0    0.200000000000
    10.0 0.2 0.20 7    0.073371208038 7    0.506362192214
    10.0 0.2 0.20 10   0.050462650440 10   0.692930517751
    10.0 0.2 0.20 14   0.026408025750 14   0.842395474165))

;; note: uses a looser tolerance than the default 1.0e-6 because the underlying
;; inverse-gaussian icdf (SSJ's iterative root finder) differs from R's own
;; algorithm by up to ~1.5e-5 at these points.
(t/deftest zaig-icdf
  (t/are [mu sigma nu p vq]
      (let [dist (sut/distribution :zaig {:mu mu :sigma sigma :nu nu})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    1.0  1.0 0.10 0.01 0.000000000000
    1.0  1.0 0.10 0.05 0.000000000000
    1.0  1.0 0.10 0.15 0.190712515192
    1.0  1.0 0.10 0.3  0.353009176719
    1.0  1.0 0.10 0.5  0.597434553041
    1.0  1.0 0.10 0.7  0.996411820433
    1.0  1.0 0.10 0.85 1.621863305366
    1.0  1.0 0.10 0.99 4.840194437883
    0.1  2.0 0.30 0.01 0.000000000000
    0.1  2.0 0.30 0.05 0.000000000000
    0.1  2.0 0.30 0.15 0.000000000000
    0.1  2.0 0.30 0.3  0.000000000000
    0.1  2.0 0.30 0.5  0.059755334015
    0.1  2.0 0.30 0.7  0.093280653371
    0.1  2.0 0.30 0.85 0.134987600441
    0.1  2.0 0.30 0.99 0.301382301195
    5.0  0.5 0.05 0.01 0.000000000000
    5.0  0.5 0.05 0.05 0.000000000000
    5.0  0.5 0.05 0.15 1.044072955011
    5.0  0.5 0.05 0.3  1.741212051827
    5.0  0.5 0.05 0.5  2.938658382073
    5.0  0.5 0.05 0.7  5.025949006062
    5.0  0.5 0.05 0.85 8.470609980487
    5.0  0.5 0.05 0.99 27.425661235617
    2.0  3.0 0.40 0.01 0.000000000000
    2.0  3.0 0.40 0.05 0.000000000000
    2.0  3.0 0.40 0.15 0.000000000000
    2.0  3.0 0.40 0.3  0.000000000000
    2.0  3.0 0.40 0.5  0.055674779428
    2.0  3.0 0.40 0.7  0.216327917661
    2.0  3.0 0.40 0.85 0.812168548537
    2.0  3.0 0.40 0.99 24.124346197443
    10.0 0.2 0.20 0.01 0.000000000000
    10.0 0.2 0.20 0.05 0.000000000000
    10.0 0.2 0.20 0.15 0.000000000000
    10.0 0.2 0.20 0.3  4.282002338354
    10.0 0.2 0.20 0.5  6.913593872595
    10.0 0.2 0.20 0.7  10.141590907012
    10.0 0.2 0.20 0.85 14.295183958309
    10.0 0.2 0.20 0.99 30.989151459456))

(t/deftest zaig-mv
  (t/are [mu sigma nu mean-v var-v]
      (let [dist (sut/distribution :zaig {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    1.0  1.0 0.10 0.900000000000 0.990000000000
    0.1  2.0 0.30 0.070000000000 0.004900000000
    5.0  0.5 0.05 4.750000000000 30.875000000000
    2.0  3.0 0.40 1.200000000000 44.160000000000
    10.0 0.2 0.20 8.000000000000 48.000000000000))

;; reference values from R's EnvStats package: dgevd/pgevd/qgevd
;; note: EnvStats parameterizes GEVD by (location, scale, shape=kappa) where
;; fastmath's xi = -kappa (i.e. call R with shape = -xi).

(t/deftest gev
  (t/testing "registered under both :generalized-extreme-value and :gev keys"
    (let [d1 (sut/distribution :generalized-extreme-value {:mu 2.0 :sigma 1.5 :xi 0.3})
          d2 (sut/distribution :gev {:mu 2.0 :sigma 1.5 :xi 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 2.0) (sut/pdf d2 2.0)))
      (t/is (m/delta-eq (sut/cdf d1 2.0) (sut/cdf d2 2.0)))))
  (t/testing "defaults match R's GEVD() defaults (location=0, scale=1, shape=0, i.e. standard Gumbel)"
    (let [dist (sut/distribution :gev)]
      (t/is (m/delta-eq 0.367879441171 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.367879441171 (sut/cdf dist 0)))))
  (t/testing "xi > 0 (Frechet): pdf/cdf are exactly 0 below the lower bound mu - sigma/xi"
    (let [dist (sut/distribution :gev {:mu 2.0 :sigma 1.5 :xi 0.3})]
      (t/is (m/delta-eq -3.0 (sut/lower-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -3.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -3.0)))))
  (t/testing "xi < 0 (Weibull-type): pdf is 0 and cdf is 1 at/above the upper bound mu - sigma/xi"
    (let [dist (sut/distribution :gev {:mu -1.0 :sigma 0.8 :xi -0.4})]
      (t/is (m/delta-eq 1.0 (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 1.0 (sut/cdf dist 1.0)))))
  (t/testing "pdf at +/- infinity is exactly 0.0, never NaN, across xi < 0, xi = 0 and xi > 0"
    (doseq [xi [-2.0 -0.9 -0.5 -0.1 0.0 0.1 0.5 0.9 2.0]]
      (let [dist (sut/distribution :gev {:mu 0.0 :sigma 1.0 :xi xi})]
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)) (str "xi=" xi))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)) (str "xi=" xi)))))
  (t/are [mu sigma xi vd d vp p]
      (let [dist (sut/distribution :gev {:mu mu :sigma sigma :xi xi})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.0  1.0 0.00  -2  0.004566281420 -2  0.000617978989
    0.0  1.0 0.00  -1  0.179374078734 -1  0.065988035845
    0.0  1.0 0.00   0  0.367879441171  0  0.367879441171
    0.0  1.0 0.00   1  0.254646380044  1  0.692200627555
    0.0  1.0 0.00   2  0.118204951593  2  0.873423018493
    0.0  1.0 0.00   5  0.006692699678  5  0.993284702068
    2.0  1.5 0.30  -3  0.000000000000 -3  0.000000000000
    2.0  1.5 0.30  -2  0.000000000000 -2  0.000000000000
    2.0  1.5 0.30   0  0.025199773253  0  0.004131838237
    2.0  1.5 0.30   2  0.245252960781  2  0.367879441171
    2.0  1.5 0.30   5  0.070588728537  5  0.811608418650
    2.0  1.5 0.30  10  0.010179377994 10  0.959467661832
    -1.0 0.8 -0.40 -3  0.012351349788 -3  0.003493489277
    -1.0 0.8 -0.40 -1  0.459849301464 -1  0.367879441171
    -1.0 0.8 -0.40  0  0.370332542002  0  0.837966885579
    -1.0 0.8 -0.40  0.5 0.151442692887 0.5 0.969233234476
    -1.0 0.8 -0.40  0.9 0.013967614543 0.9 0.999441139227
    -1.0 0.8 -0.40  1  0.000000000000  1  1.000000000000
    5.0  2.0 0.15   0  0.000000001979  0  0.000000000108
    5.0  2.0 0.15   5  0.183939720586  5  0.367879441171
    5.0  2.0 0.15  10  0.038608034275 10  0.887213970332
    5.0  2.0 0.15  20  0.001535905817 20  0.993450908007
    1.0  0.5 -0.20 -2  0.000000000000 -2  0.000000000000
    1.0  0.5 -0.20  0  0.035465181275  0  0.004615938837
    1.0  0.5 -0.20  1  0.735758882343  1  0.367879441171
    1.0  0.5 -0.20  2  0.239808326453  2  0.925186444647
    1.0  0.5 -0.20  3  0.003198976164  3  0.999680051195
    1.0  0.5 -0.20  3.5 0.000000000000 3.5 1.000000000000))

(t/deftest gev-icdf
  (t/are [mu sigma xi p vq]
      (let [dist (sut/distribution :gev {:mu mu :sigma sigma :xi xi})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.0  1.0 0.00  0.01 -1.527179625808
    0.0  1.0 0.00  0.10 -0.834032445248
    0.0  1.0 0.00  0.25 -0.326634259978
    0.0  1.0 0.00  0.50  0.366512920582
    0.0  1.0 0.00  0.75  1.245899323707
    0.0  1.0 0.00  0.90  2.250367327312
    0.0  1.0 0.00  0.99  4.600149226777
    2.0  1.5 0.30  0.01  0.162250711969
    2.0  1.5 0.30  0.10  0.893187297636
    2.0  1.5 0.30  0.25  1.533288591217
    2.0  1.5 0.30  0.50  2.581132923157
    2.0  1.5 0.30  0.75  4.266012902591
    2.0  1.5 0.30  0.90  6.821247103910
    2.0  1.5 0.30  0.99 16.875397900026
    -1.0 0.8 -0.40 0.01 -2.684073358344
    -1.0 0.8 -0.40 0.10 -1.792005510741
    -1.0 0.8 -0.40 0.25 -1.279146166729
    -1.0 0.8 -0.40 0.50 -0.727269801205
    -1.0 0.8 -0.40 0.75 -0.215052700899
    -1.0 0.8 -0.40 0.90  0.186980147054
    -1.0 0.8 -0.40 0.99  0.682384107070
    5.0  2.0 0.15  0.01  2.270231202708
    5.0  2.0 0.15  0.10  3.432058170960
    5.0  2.0 0.15  0.25  4.362476776871
    5.0  2.0 0.15  0.50  5.753549988120
    5.0  2.0 0.15  0.75  7.739846981261
    5.0  2.0 0.15  0.90 10.353557714348
    5.0  2.0 0.15  0.99 18.250135481386
    1.0  0.5 -0.20 0.01  0.106958702753
    1.0  0.5 -0.20 0.10  0.546185987395
    1.0  0.5 -0.20 0.25  0.831230288873
    1.0  0.5 -0.20 0.50  1.176701024671
    1.0  0.5 -0.20 0.75  1.551400582390
    1.0  0.5 -0.20 0.90  1.906046725762
    1.0  0.5 -0.20 0.99  2.503732131701))

(t/deftest gev-mv
  (t/testing "mean and variance analytic formulas (mean finite for xi<1, variance finite for xi<0.5)"
    (t/are [mu sigma xi mean-v var-v]
        (let [dist (sut/distribution :gev {:mu mu :sigma sigma :xi xi})]
          (and (m/delta-eq mean-v (sut/mean dist))
               (m/delta-eq var-v (sut/variance dist))))
      0.0  1.0 0.00   0.577215664902  1.644934066848
      2.0  1.5 0.30   3.490276663238 13.330297428573
      -1.0 0.8 -0.40 -0.774527635006  0.576586756520
      5.0  2.0 0.15   6.499783159313 10.744047586264
      1.0  0.5 -0.20  1.204578144001  0.276437362394))
  (t/testing "mean is +Inf for xi >= 1, variance is +Inf for xi >= 0.5"
    (let [dist (sut/distribution :gev {:mu 0.0 :sigma 1.0 :xi 1.0})]
      (t/is (Double/isInfinite (sut/mean dist))))
    (let [dist (sut/distribution :gev {:mu 0.0 :sigma 1.0 :xi 0.5})]
      (t/is (Double/isInfinite (sut/variance dist))))))

;; reference values from R's glogis package: dglogis/pglogis/qglogis

(t/deftest generalized-logistic
  (t/testing "defaults match R's glogis() defaults (location=0, scale=1, shape=1, i.e. standard logistic)"
    (let [dist (sut/distribution :generalized-logistic)]
      (t/is (m/delta-eq 0.25 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.5 (sut/cdf dist 0)))))
  (t/testing "alpha=1 matches fastmath's own standard logistic distribution"
    (let [gl (sut/distribution :generalized-logistic {:mu 2.0 :sigma 1.5 :alpha 1.0})
          lg (sut/distribution :logistic {:mu 2.0 :s 1.5})]
      (t/is (m/delta-eq (sut/pdf gl 1.0) (sut/pdf lg 1.0)))
      (t/is (m/delta-eq (sut/cdf gl 1.0) (sut/cdf lg 1.0)))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :generalized-logistic {:mu 2.0 :sigma 1.5 :alpha 2.5})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))))
  (t/are [mu sigma alpha vd d vp p]
      (let [dist (sut/distribution :generalized-logistic {:mu mu :sigma sigma :alpha alpha})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.0  1.0 1.0 -5  0.006648056671 -5  0.006692850924
    0.0  1.0 1.0 -2  0.104993585404 -2  0.119202922022
    0.0  1.0 1.0 -1  0.196611933241 -1  0.268941421370
    0.0  1.0 1.0  0  0.250000000000  0  0.500000000000
    0.0  1.0 1.0  1  0.196611933241  1  0.731058578630
    0.0  1.0 1.0  2  0.104993585404  2  0.880797077978
    0.0  1.0 1.0  5  0.006648056671  5  0.993307149076
    2.0  1.5 2.5 -3  0.000354362143 -3  0.000220202193
    2.0  1.5 2.5  0  0.026216231233  0  0.019876052855
    2.0  1.5 2.5  2  0.147313912747  2  0.176776695297
    2.0  1.5 2.5  4  0.193713203280  4  0.557158058005
    2.0  1.5 2.5  6  0.091542104044  6  0.845405030772
    2.0  1.5 2.5 10  0.007912077619 10  0.988031368693
    -1.0 0.8 0.4 -4  0.107992788128 -4  0.221065070147
    -1.0 0.8 0.4 -2  0.213131081929 -2  0.548388318522
    -1.0 0.8 0.4 -1  0.189464570814 -1  0.757858283255
    -1.0 0.8 0.4  0  0.100675994453  0  0.904139485351
    -1.0 0.8 0.4  1  0.036750891006  1  0.968936797461
    -1.0 0.8 0.4  3  0.003337448559  3  0.997317465049
    5.0  3.0 5.0 -2  0.000008201752 -2  0.000005398255
    5.0  3.0 5.0  2  0.001714310189  2  0.001406981798
    5.0  3.0 5.0  5  0.026041666667  5  0.031250000000
    5.0  3.0 5.0  8  0.093598164911  8  0.208814613459
    5.0  3.0 5.0 12  0.092751679053 12  0.629538581927
    5.0  3.0 5.0 20  0.010786429452 20  0.966980699920
    1.0  0.5 0.1 -2  0.109463823196 -2  0.548675784430
    1.0  0.5 0.1  0  0.142408053720  0  0.808404440027
    1.0  0.5 0.1  1  0.093303299154  1  0.933032991537
    1.0  0.5 0.1  2  0.023539892954  2  0.987387412757
    1.0  0.5 0.1  3  0.003590718946  3  0.998186653312
    1.0  0.5 0.1  5  0.000067067778  5  0.999966459925))

(t/deftest generalized-logistic-icdf
  (t/are [mu sigma alpha p vq]
      (let [dist (sut/distribution :generalized-logistic {:mu mu :sigma sigma :alpha alpha})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.0  1.0 1.0 0.01  -4.595119850135
    0.0  1.0 1.0 0.10  -2.197224577336
    0.0  1.0 1.0 0.25  -1.098612288668
    0.0  1.0 1.0 0.50   0.000000000000
    0.0  1.0 1.0 0.75   1.098612288668
    0.0  1.0 1.0 0.90   2.197224577336
    0.0  1.0 1.0 0.99   4.595119850135
    2.0  1.5 2.5 0.01  -0.504267252222
    2.0  1.5 2.5 0.10   1.379962754749
    2.0  1.5 2.5 0.25   2.449427284796
    2.0  1.5 2.5 0.50   3.711459868428
    2.0  1.5 2.5 0.75   5.156152943194
    2.0  1.5 2.5 0.90   6.718267927343
    2.0  1.5 2.5 0.99  10.271643827128
    -1.0 0.8 0.4 0.01 -10.210332371936
    -1.0 0.8 0.4 0.10  -5.602636355407
    -1.0 0.8 0.4 0.25  -3.747189763588
    -1.0 0.8 0.4 0.50  -2.230672133143
    -1.0 0.8 0.4 0.75  -1.041163352565
    -1.0 0.8 0.4 0.90  -0.040410578299
    -1.0 0.8 0.4 0.99   1.937015416585
    5.0  3.0 5.0 0.01   3.759925509497
    5.0  3.0 5.0 0.10   6.608978076227
    5.0  3.0 5.0 0.25   8.422919736856
    5.0  3.0 5.0 0.50  10.717506464418
    5.0  3.0 5.0 0.75  13.479293293229
    5.0  3.0 5.0 0.90  16.547752060556
    5.0  3.0 5.0 0.99  23.625745811830
    1.0  0.5 0.1 0.01 -22.025850929940
    1.0  0.5 0.1 0.10 -10.512925464920
    1.0  0.5 0.1 0.25  -5.931471328762
    1.0  0.5 0.1 0.50  -2.465247382976
    1.0  0.5 0.1 0.75  -0.409429721740
    1.0  0.5 0.1 0.90   0.687573327264
    1.0  0.5 0.1 0.99   2.123445809029))

(t/deftest generalized-logistic-mv
  (t/are [mu sigma alpha mean-v var-v]
      (let [dist (sut/distribution :generalized-logistic {:mu mu :sigma sigma :alpha alpha})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.0  1.0 1.0  0.000000000000  3.289868133696
    2.0  1.5 2.5  3.920558458320  4.804406601634
    -1.0 0.8 0.4 -2.587335103747  5.708986020722
    5.0  3.0 5.0 11.250000000000 16.796313203268
    1.0  0.5 0.1 -3.923269637755 25.769558304410))

;; reference values from R's evd package: dgpd/pgpd/qgpd
;; note: for pdf, R's dgpd treats the support as the OPEN interval (mu, ...) and
;; returns 0 exactly at x=mu, whereas fastmath (and the standard textbook
;; definition, e.g. Coles 2001) uses the CLOSED interval [mu, ...) with
;; pdf(mu) = 1/sigma; reference pdf points below therefore avoid x=mu exactly
;; (using mu + 0.001 instead) while cdf points at/around mu agree either way.

(t/deftest gpd
  (t/testing "registered under both :generalized-pareto and :gpd keys"
    (let [d1 (sut/distribution :generalized-pareto {:mu 2.0 :sigma 1.5 :xi 0.3})
          d2 (sut/distribution :gpd {:mu 2.0 :sigma 1.5 :xi 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 3.0) (sut/pdf d2 3.0)))
      (t/is (m/delta-eq (sut/cdf d1 3.0) (sut/cdf d2 3.0)))))
  (t/testing "defaults match R's gpd() defaults (loc=0, scale=1, shape=0, i.e. standard exponential)"
    (let [dist (sut/distribution :gpd)]
      (t/is (m/delta-eq 1.0 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0)))))
  (t/testing "pdf at the lower bound mu equals 1/sigma (closed-interval convention)"
    (doseq [xi [-0.4 0.0 0.3]]
      (let [dist (sut/distribution :gpd {:mu 2.0 :sigma 1.5 :xi xi})]
        (t/is (m/delta-eq (/ 1.0 1.5) (sut/pdf dist 2.0))))))
  (t/testing "pdf/cdf are exactly 0 below the lower bound mu, for any xi"
    (doseq [xi [-0.4 0.0 0.3]]
      (let [dist (sut/distribution :gpd {:mu 2.0 :sigma 1.5 :xi xi})]
        (t/is (m/delta-eq 2.0 (sut/lower-bound dist)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist 1.0)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist 1.0))))))
  (t/testing "xi < 0: pdf is 0 and cdf is 1 at/above the finite upper bound mu - sigma/xi"
    (let [dist (sut/distribution :gpd {:mu -1.0 :sigma 0.8 :xi -0.4})]
      (t/is (m/delta-eq 1.0 (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 1.0 (sut/cdf dist 1.0)))))
  (t/testing "xi >= 0: upper bound is +Inf"
    (doseq [xi [0.0 0.3]]
      (let [dist (sut/distribution :gpd {:mu 2.0 :sigma 1.5 :xi xi})]
        (t/is (Double/isInfinite (sut/upper-bound dist))))))
  (t/testing "pdf at +/- infinity is exactly 0.0, never NaN"
    (doseq [xi [-2.0 -0.9 -0.5 -0.1 0.0 0.1 0.5 0.9 2.0]]
      (let [dist (sut/distribution :gpd {:mu 0.0 :sigma 1.0 :xi xi})]
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)) (str "xi=" xi))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)) (str "xi=" xi)))))
  (t/are [mu sigma xi vd d vp p]
      (let [dist (sut/distribution :gpd {:mu mu :sigma sigma :xi xi})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.0  1.0 0.00  -1     0.000000000000 -1     0.000000000000
    0.0  1.0 0.00   0.001 0.999000499833  0.001 0.000999500167
    0.0  1.0 0.00   0.5   0.606530659713  0.5   0.393469340287
    0.0  1.0 0.00   1     0.367879441171  1     0.632120558829
    0.0  1.0 0.00   2     0.135335283237  2     0.864664716763
    0.0  1.0 0.00   5     0.006737946999  5     0.993262053001
    2.0  1.5 0.30   1     0.000000000000  1     0.000000000000
    2.0  1.5 0.30   2.001 0.666089196907  2.001 0.000666377880
    2.0  1.5 0.30   3     0.302545019573  3     0.455418964768
    2.0  1.5 0.30   5     0.086973874241  5     0.791262701822
    2.0  1.5 0.30  10     0.010609401858 10     0.958623332753
    2.0  1.5 0.30  20     0.000895278463 20     0.993822578607
    -1.0 0.8 -0.40 -2     0.000000000000 -2     0.000000000000
    -1.0 0.8 -0.40 -0.999 1.249062617197 -0.999 0.001249531289
    -1.0 0.8 -0.40  0     0.441941738242  0     0.823223304703
    -1.0 0.8 -0.40  0.5   0.156250000000  0.5   0.968750000000
    -1.0 0.8 -0.40  0.9   0.013975424859  0.9   0.999440983006
    -1.0 0.8 -0.40  1     0.000000000000  1     1.000000000000
    5.0  2.0 0.15   4     0.000000000000  4     0.000000000000
    5.0  2.0 0.15   5.001 0.499712593415  5.001 0.000499856281
    5.0  2.0 0.15   7     0.171245786686  7     0.606134690622
    5.0  2.0 0.15  10     0.043516035101 10     0.880330903472
    5.0  2.0 0.15  20     0.001546030915 20     0.993429368610
    5.0  2.0 0.15  50     0.000006092631 50     0.999946689476
    1.0  0.5 -0.20  0     0.000000000000  0     0.000000000000
    1.0  0.5 -0.20  1.001 1.996801919488  1.001 0.001998400640
    1.0  0.5 -0.20  2     0.259200000000  2     0.922240000000
    1.0  0.5 -0.20  2.5   0.051200000000  2.5   0.989760000000
    1.0  0.5 -0.20  2.9   0.006635520000  2.9   0.999203737600
    1.0  0.5 -0.20  3     0.003200000000  3     0.999680000000))

(t/deftest gpd-icdf
  (t/are [mu sigma xi p vq]
      (let [dist (sut/distribution :gpd {:mu mu :sigma sigma :xi xi})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.0  1.0 0.00  0.01  0.010050335854
    0.0  1.0 0.00  0.10  0.105360515658
    0.0  1.0 0.00  0.25  0.287682072452
    0.0  1.0 0.00  0.50  0.693147180560
    0.0  1.0 0.00  0.75  1.386294361120
    0.0  1.0 0.00  0.90  2.302585092994
    0.0  1.0 0.00  0.99  4.605170185988
    2.0  1.5 0.30  0.01  2.015098253720
    2.0  1.5 0.30  0.10  2.160564987141
    2.0  1.5 0.30  0.25  2.450691787847
    2.0  1.5 0.30  0.50  3.155722066725
    2.0  1.5 0.30  0.75  4.578582832552
    2.0  1.5 0.30  0.90  6.976311574844
    2.0  1.5 0.30  0.99 16.905358527675
    -1.0 0.8 -0.40 0.01 -0.991975871162
    -1.0 0.8 -0.40 0.10 -0.917463031028
    -1.0 0.8 -0.40 0.25 -0.782602457966
    -1.0 0.8 -0.40 0.50 -0.515716566510
    -1.0 0.8 -0.40 0.75 -0.148698354997
    -1.0 0.8 -0.40 0.90  0.203785658893
    -1.0 0.8 -0.40 0.99  0.683021361508
    5.0  2.0 0.15  0.01  5.020115830711
    5.0  2.0 0.15  0.10  5.212394963748
    5.0  2.0 0.15  0.25  5.587958800624
    5.0  2.0 0.15  0.50  6.460926294238
    5.0  2.0 0.15  0.75  8.081925511266
    5.0  2.0 0.15  0.90 10.500500594970
    5.0  2.0 0.15  0.99 18.270164199585
    1.0  0.5 -0.20 0.01  1.005020120846
    1.0  0.5 -0.20 0.10  1.052129094098
    1.0  0.5 -0.20 0.25  1.139781221763
    1.0  0.5 -0.20 0.50  1.323623591760
    1.0  0.5 -0.20 0.75  1.605354291862
    1.0  0.5 -0.20 0.90  1.922606638800
    1.0  0.5 -0.20 0.99  2.504732073616))

(t/deftest gpd-mv
  (t/testing "mean and variance analytic formulas (mean finite for xi<1, variance finite for xi<0.5)"
    (t/are [mu sigma xi mean-v var-v]
        (let [dist (sut/distribution :gpd {:mu mu :sigma sigma :xi xi})]
          (and (m/delta-eq mean-v (sut/mean dist))
               (m/delta-eq var-v (sut/variance dist))))
      0.0  1.0 0.00   1.000000000000  1.000000000000
      2.0  1.5 0.30   4.142857142857 11.479591836735
      -1.0 0.8 -0.40 -0.428571428571  0.181405895692
      5.0  2.0 0.15   7.352941176471  7.909045971330
      1.0  0.5 -0.20  1.416666666667  0.124007936508))
  (t/testing "mean is +Inf for xi >= 1, variance is +Inf for xi >= 0.5"
    (let [dist (sut/distribution :gpd {:mu 0.0 :sigma 1.0 :xi 1.0})]
      (t/is (Double/isInfinite (sut/mean dist))))
    (let [dist (sut/distribution :gpd {:mu 0.0 :sigma 1.0 :xi 0.5})]
      (t/is (Double/isInfinite (sut/variance dist))))))

;; reference values from R's reliaR package: dgen.exp/pgen.exp/qgen.exp
;; note: reliaR's dgen.exp/pgen.exp require x > 0 strictly (error otherwise), so
;; boundary behavior at/below x=0 is checked separately against the analytic
;; pdf(0) formula (0^(alpha-1) convention) rather than against R.

(t/deftest generalized-exponential
  (t/testing "registered under both :generalized-exponential and :ge keys"
    (let [d1 (sut/distribution :generalized-exponential {:alpha 2.5 :lambda 1.3})
          d2 (sut/distribution :ge {:alpha 2.5 :lambda 1.3})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "defaults match R's reliaR::gen.exp defaults (alpha=1, lambda=1, i.e. standard exponential rate 1)"
    (let [dist (sut/distribution :generalized-exponential)]
      (t/is (m/delta-eq 1.0 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0)))))
  (t/testing "support is [0, +Inf) for any alpha, lambda"
    (let [dist (sut/distribution :generalized-exponential {:alpha 2.5 :lambda 1.3})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))))
  (t/testing "pdf/cdf are exactly 0 for negative x"
    (doseq [alpha [0.4 1.0 2.5]]
      (let [dist (sut/distribution :generalized-exponential {:alpha alpha :lambda 2.0})]
        (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0))))))
  (t/testing "pdf(0) follows the 0^(alpha-1) convention: +Inf for alpha<1, lambda for alpha=1, 0 for alpha>1"
    (let [dist (sut/distribution :generalized-exponential {:alpha 0.4 :lambda 2.0})]
      (t/is (Double/isInfinite (sut/pdf dist 0.0))))
    (let [dist (sut/distribution :generalized-exponential {:alpha 1.0 :lambda 2.0})]
      (t/is (m/delta-eq 2.0 (sut/pdf dist 0.0))))
    (let [dist (sut/distribution :generalized-exponential {:alpha 2.5 :lambda 2.0})]
      (t/is (m/delta-eq 0.0 (sut/pdf dist 0.0)))))
  (t/testing "pdf is never NaN at +/- infinity"
    (doseq [alpha [0.4 1.0 2.5]]
      (let [dist (sut/distribution :generalized-exponential {:alpha alpha :lambda 2.0})]
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf))))))
  (t/are [alpha lambda vd d vp p]
      (let [dist (sut/distribution :generalized-exponential {:alpha alpha :lambda lambda})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 1.0 0.001 0.999000499833 0.001 0.000999500167
    1.0 1.0 0.5   0.606530659713 0.5   0.393469340287
    1.0 1.0 1     0.367879441171 1     0.632120558829
    1.0 1.0 2     0.135335283237 2     0.864664716763
    1.0 1.0 5     0.006737946999 5     0.993262053001
    2.5 1.3 0.001 0.000151988390 0.001 0.000000060835
    2.5 1.3 0.3   0.403829070404 0.3   0.059267295567
    2.5 1.3 0.7   0.604165379667 0.7   0.275931165210
    2.5 1.3 1.5   0.367309304041 1.5   0.681351704464
    2.5 1.3 3     0.063798896631 3     0.950160876055
    2.5 1.3 6     0.001330820337 6     0.998975977311
    0.4 2.0 0.001 33.255534611768 0.001 0.083222030820
    0.4 2.0 0.1   1.824889667912 0.1   0.505044507267
    0.4 2.0 0.5   0.387539126342 0.5   0.832376798262
    0.4 2.0 1     0.118138752508 1     0.943493896543
    0.4 2.0 3     0.001985956824 3     0.999007760848
    5.0 0.5 0.001 0.000000000156 0.001 0.000000000000
    5.0 0.5 0.5   0.004661232103 0.5   0.000529563356
    5.0 0.5 2     0.146840274691 2     0.100925190275
    5.0 0.5 4     0.189121975400 4     0.483324364147
    5.0 0.5 8     0.042525513929 8     0.911715550327
    5.0 0.5 15    0.001379654439 15    0.997237635481
    0.8 3.0 0.001 7.649170881158 0.001 0.009575820150
    0.8 3.0 0.05  3.063951030307 0.05  0.206605081144
    0.8 3.0 0.2   1.544411313722 0.2   0.529037323561
    0.8 3.0 0.5   0.563248270154 0.5   0.817106394201
    0.8 3.0 1     0.120715658753 1     0.959967984309
    0.8 3.0 2     0.005951958840 2     0.998016506234))

(t/deftest generalized-exponential-icdf
  (t/are [alpha lambda p vq]
      (let [dist (sut/distribution :generalized-exponential {:alpha alpha :lambda lambda})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 1.0 0.01  0.010050335854
    1.0 1.0 0.10  0.105360515658
    1.0 1.0 0.25  0.287682072452
    1.0 1.0 0.50  0.693147180560
    1.0 1.0 0.75  1.386294361120
    1.0 1.0 0.90  2.302585092994
    1.0 1.0 0.99  4.605170185988
    2.5 1.3 0.01  0.132735825318
    2.5 1.3 0.10  0.390519902844
    2.5 1.3 0.25  0.657027641778
    2.5 1.3 0.50  1.090947782956
    2.5 1.3 0.75  1.707057531623
    2.5 1.3 0.90  2.452043198327
    2.5 1.3 0.99  4.244961040328
    0.4 2.0 0.01  0.000005000025
    0.4 2.0 0.10  0.001583644113
    0.4 2.0 0.25  0.015874349157
    0.4 2.0 0.50  0.097263892485
    0.4 2.0 0.75  0.333875495212
    0.4 2.0 0.90  0.731444033135
    0.4 2.0 0.99  1.848197555183
    5.0 0.5 0.01  1.015351747393
    5.0 0.5 0.10  1.993686088016
    5.0 0.5 0.25  2.836464235685
    5.0 0.5 0.50  4.088929848502
    5.0 0.5 0.75  5.767935024467
    5.0 0.5 0.90  7.740645579967
    5.0 0.5 0.99 12.421184008895
    0.8 3.0 0.01  0.001055762742
    0.8 3.0 0.10  0.019292388445
    0.8 3.0 0.25  0.064842594990
    0.8 3.0 0.50  0.181833415318
    0.8 3.0 0.75  0.399058259260
    0.8 3.0 0.90  0.697450497298
    0.8 3.0 0.99  1.461093519752))

(t/deftest generalized-exponential-mv
  (t/are [alpha lambda mean-v var-v]
      (let [dist (sut/distribution :generalized-exponential {:alpha alpha :lambda lambda})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    1.0 1.0 1.000000000000 1.000000000000
    2.5 1.3 1.292594081190 0.777855805176
    0.4 2.0 0.257915560158 0.154894369080
    5.0 0.5 4.566666666667 5.854444444444
    0.8 3.0 0.287402366065 0.100884436594))

;; reference values from R's gamlss.dist package: dGG/pGG/qGG
;; note: for nu=0 the GG distribution degenerates exactly to log-normal with
;; scale=log(mu), shape=sigma (fastmath delegates to [[sut/log-normal]] directly
;; in that case), verified separately below rather than against dGG/pGG (whose
;; nu=0 numerics agree with dlnorm/plnorm, as cross-checked during development).

(t/deftest generalized-gamma
  (t/testing "registered under both :generalized-gamma and :gg keys"
    (let [d1 (sut/distribution :generalized-gamma {:mu 2.0 :sigma 0.7 :nu 0.8})
          d2 (sut/distribution :gg {:mu 2.0 :sigma 0.7 :nu 0.8})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "defaults match R's GG() defaults (mu=1, sigma=0.5, nu=1)"
    (let [dist (sut/distribution :gg)]
      (t/is (m/delta-eq 0.7814672592526585 (sut/pdf dist 1)))))
  (t/testing "support is (0, +Inf); pdf/cdf are exactly 0 at/below 0"
    (let [dist (sut/distribution :gg {:mu 2.0 :sigma 0.7 :nu 0.8})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))))
  (t/testing "nu=0 degenerates exactly to log-normal(scale=log(mu), shape=sigma)"
    (let [mu 2.0 sigma 0.6
          gg (sut/distribution :gg {:mu mu :sigma sigma :nu 0.0})
          ln (sut/distribution :log-normal {:scale (m/log mu) :shape sigma})]
      (doseq [x [0.5 1.0 2.0 3.0 5.0]]
        (t/is (m/delta-eq (sut/pdf gg x) (sut/pdf ln x)))
        (t/is (m/delta-eq (sut/cdf gg x) (sut/cdf ln x))))))
  (t/testing "pdf is never NaN, and cdf correctly saturates to 0/1, at extreme/infinite inputs"
    (doseq [nu [-3.0 -0.5 0.5 3.0]
            sigma [0.1 1.0 3.0]]
      (let [dist (sut/distribution :gg {:mu 1.0 :sigma sigma :nu nu})]
        (doseq [x [##-Inf ##Inf 1e-300 1e300]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf nu=" nu " sigma=" sigma " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf nu=" nu " sigma=" sigma " x=" x)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [mu sigma nu vd d vp p]
      (let [dist (sut/distribution :gg {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 0.5 1.0   0.001 0.000000042496 0.001 0.000000000011
    1.0 0.5 1.0   0.5   0.721788177262 0.5   0.142876539501
    1.0 0.5 1.0   1     0.781467259253 1     0.566529879633
    1.0 0.5 1.0   2     0.114504576991 2     0.957619888008
    1.0 0.5 1.0   5     0.000010992819 5     0.999996796280
    2.0 0.7 0.8   0.001 0.000050722994 0.001 0.000000019918
    2.0 0.7 0.8   0.5   0.273963704330 0.5   0.070457652213
    2.0 0.7 0.8   1     0.368158698418 1     0.238428280606
    2.0 0.7 0.8   2     0.277631472323 2     0.574505642912
    2.0 0.7 0.8   5     0.036568694524 5     0.952121352909
    2.0 0.7 0.8  10     0.000783129618 10    0.998975671193
    2.0 0.7 -0.8  0.001 0.000000000000 0.001 0.000000000000
    2.0 0.7 -0.8  0.5   0.058631373542 0.5   0.004756692437
    2.0 0.7 -0.8  1     0.306270670808 1     0.102135882944
    2.0 0.7 -0.8  2     0.277631472323 2     0.425494357088
    2.0 0.7 -0.8  5     0.056219099320 5     0.833895138380
    2.0 0.7 -0.8 10     0.009207271954 10    0.954831988721
    3.0 0.3 2.0   0.001 0.000000000000 0.001 0.000000000000
    3.0 0.3 2.0   1     0.034077895526 1     0.006669351815
    3.0 0.3 2.0   2     0.317473233977 2     0.163798300988
    3.0 0.3 2.0   3     0.430222476297 3     0.579816206505
    3.0 0.3 2.0   5     0.031598996538 5     0.987171289248
    3.0 0.3 2.0   8     0.000001591968 8     0.999999647482
    1.5 1.2 -1.5  0.001 0.000000000000 0.001 0.000000000000
    1.5 1.2 -1.5  0.3   0.080059593038 0.3   0.003982853524
    1.5 1.2 -1.5  1     0.245946571215 1     0.170931051574
    1.5 1.2 -1.5  2     0.128715988826 2     0.350899588165
    1.5 1.2 -1.5  5     0.039128954921 5     0.560664803529
    1.5 1.2 -1.5 10     0.014666979112 10    0.678818597631))

(t/deftest generalized-gamma-icdf
  (t/are [mu sigma nu p vq]
      (let [dist (sut/distribution :gg {:mu mu :sigma sigma :nu nu})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 0.5 1.0   0.01 0.205812171586
    1.0 0.5 1.0   0.10 0.436192390706
    1.0 0.5 1.0   0.25 0.633830052975
    1.0 0.5 1.0   0.50 0.918015187213
    1.0 0.5 1.0   0.75 1.277356871281
    1.0 0.5 1.0   0.90 1.670195767064
    1.0 0.5 1.0   0.99 2.511279378708
    2.0 0.7 0.8   0.01 0.198666338543
    2.0 0.7 0.8   0.10 0.601166215640
    2.0 0.7 0.8   0.25 1.031395380778
    2.0 0.7 0.8   0.50 1.747447486222
    2.0 0.7 0.8   0.75 2.777766844807
    2.0 0.7 0.8   0.90 4.028180516814
    2.0 0.7 0.8   0.99 7.036288410659
    2.0 0.7 -0.8  0.01 0.568481529828
    2.0 0.7 -0.8  0.10 0.993004157411
    2.0 0.7 -0.8  0.25 1.440005667674
    2.0 0.7 -0.8  0.50 2.289053051115
    2.0 0.7 -0.8  0.75 3.878241142581
    2.0 0.7 -0.8  0.90 6.653733852527
    2.0 0.7 -0.8  0.99 20.134261442295
    3.0 0.3 2.0   0.01 1.083124231504
    3.0 0.3 2.0   0.10 1.770639520570
    3.0 0.3 2.0   0.25 2.242812999084
    3.0 0.3 2.0   0.50 2.818881976455
    3.0 0.3 2.0   0.75 3.439691768877
    3.0 0.3 2.0   0.90 4.029705161373
    3.0 0.3 2.0   0.99 5.099701693443
    1.5 1.2 -1.5  0.01 0.355689755987
    1.5 1.2 -1.5  0.10 0.731057061507
    1.5 1.2 -1.5  0.25 1.361121752393
    1.5 1.2 -1.5  0.50 3.728749563788
    1.5 1.2 -1.5  0.75 17.267907831352))

(t/deftest generalized-gamma-mv
  (t/testing "mean and variance analytic formulas (mean finite for nu>-1/sigma^2, variance finite for nu>-1/(2*sigma^2))"
    (t/are [mu sigma nu mean-v var-v]
        (let [dist (sut/distribution :gg {:mu mu :sigma sigma :nu nu})]
          (and (m/delta-eq mean-v (sut/mean dist))
               (m/delta-eq var-v (sut/variance dist))))
      1.0 0.5 1.0  1.000000000000 0.250000000000
      2.0 0.7 0.8  2.092912547657 2.163023510526
      2.0 0.7 -0.8 3.468827658056 27.848704949821
      3.0 0.3 2.0  2.868664629287 0.770763244680))
  (t/testing "mean and variance are +Inf when nu is too negative relative to sigma"
    (let [dist (sut/distribution :gg {:mu 1.5 :sigma 1.2 :nu -1.5})]
      (t/is (Double/isInfinite (sut/mean dist)))
      (t/is (Double/isInfinite (sut/variance dist))))))

;; reference values from R's gnorm package: dgnorm/pgnorm/qgnorm

(t/deftest generalized-normal
  (t/testing "registered under both :generalized-normal and :gnd keys"
    (let [d1 (sut/distribution :generalized-normal {:mu 1.0 :alpha 2.0 :beta 1.5})
          d2 (sut/distribution :gnd {:mu 1.0 :alpha 2.0 :beta 1.5})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "defaults match R's gnorm() defaults (mu=0, alpha=1, beta=1, i.e. standard Laplace)"
    (let [dist (sut/distribution :gnd)]
      (t/is (m/delta-eq 0.5 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.5 (sut/cdf dist 0)))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :gnd {:mu 1.0 :alpha 2.0 :beta 1.5})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))))
  (t/testing "symmetric about mu: pdf(mu-d) == pdf(mu+d), cdf(mu)=0.5"
    (doseq [beta [0.5 1.0 2.0 5.0]]
      (let [dist (sut/distribution :gnd {:mu 2.0 :alpha 1.3 :beta beta})]
        (t/is (m/delta-eq (sut/pdf dist 2.7) (sut/pdf dist 1.3)))
        (t/is (m/delta-eq 0.5 (sut/cdf dist 2.0))))))
  (t/testing "beta=2 matches normal(mu, sd=alpha/sqrt(2))"
    (let [gn (sut/distribution :gnd {:mu 0.0 :alpha 1.0 :beta 2.0})
          no (sut/distribution :normal {:mu 0.0 :sd (/ 1.0 (Math/sqrt 2.0))})]
      (doseq [x [-2.0 -0.5 0.0 0.5 2.0]]
        (t/is (m/delta-eq (sut/pdf gn x) (sut/pdf no x)))
        (t/is (m/delta-eq (sut/cdf gn x) (sut/cdf no x))))))
  (t/testing "beta=1 matches laplace(mu, scale=alpha)"
    (let [gn (sut/distribution :gnd {:mu 0.0 :alpha 1.0 :beta 1.0})
          lp (sut/distribution :laplace {:mu 0.0 :scale 1.0})]
      (doseq [x [-2.0 -0.5 0.0 0.5 2.0]]
        (t/is (m/delta-eq (sut/pdf gn x) (sut/pdf lp x)))
        (t/is (m/delta-eq (sut/cdf gn x) (sut/cdf lp x))))))
  (t/testing "pdf/cdf are never NaN at +/- infinity and other extreme inputs"
    (doseq [beta [0.1 0.5 1.0 2.0 5.0]]
      (let [dist (sut/distribution :gnd {:mu 0.0 :alpha 1.0 :beta beta})]
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf)))
        (t/is (not (Double/isNaN (sut/pdf dist 1e300))))
        (t/is (not (Double/isNaN (sut/cdf dist 1e300)))))))
  (t/are [mu alpha beta vd d vp p]
      (let [dist (sut/distribution :gnd {:mu mu :alpha alpha :beta beta})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.0  1.0 1.0 -5   0.003368973500 -5   0.003368973500
    0.0  1.0 1.0 -1   0.183939720586 -1   0.183939720586
    0.0  1.0 1.0  0   0.500000000000  0   0.500000000000
    0.0  1.0 1.0  1   0.183939720586  1   0.816060279414
    0.0  1.0 1.0  2   0.067667641618  2   0.932332358382
    1.0  2.0 1.5 -6   0.000396878664 -6   0.000270673074
    1.0  2.0 1.5 -2   0.044108761833 -2   0.042432043134
    1.0  2.0 1.5  0   0.194459197630  0   0.258250670761
    1.0  2.0 1.5  3   0.101877972681  3   0.887591235992
    1.0  2.0 1.5  6   0.005317103369  6   0.995809984618
    -2.0 0.5 0.7 -4   0.056430881565 -4   0.069571703267
    -2.0 0.5 0.7 -3   0.155637188365 -3   0.165680834041
    -2.0 0.5 0.7 -1   0.155637188365 -1   0.834319165959
    -2.0 0.5 0.7  0   0.056430881565  0   0.930428296733
    -2.0 0.5 0.7  1   0.023733517493  1   0.967852678294
    0.0  1.0 2.0 -4   0.000000063491 -4   0.000000007709
    0.0  1.0 2.0 -1   0.207553748710 -1   0.078649603525
    0.0  1.0 2.0  1   0.207553748710  1   0.921350396475
    0.0  1.0 2.0  2   0.010333492677  2   0.997661132509
    3.0  1.5 5.0  1.5 0.133555494462  1.5 0.026190215640
    3.0  1.5 5.0  2.5 0.361550545709  2.5 0.318603623378
    3.0  1.5 5.0  3.5 0.361550545709  3.5 0.681396376622
    3.0  1.5 5.0  4.5 0.133555494462  4.5 0.973809784360))

(t/deftest generalized-normal-icdf
  (t/are [mu alpha beta p vq]
      (let [dist (sut/distribution :gnd {:mu mu :alpha alpha :beta beta})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.0  1.0 1.0 0.01 -3.912023005428
    0.0  1.0 1.0 0.10 -1.609437912434
    0.0  1.0 1.0 0.25 -0.693147180560
    0.0  1.0 1.0 0.50  0.000000000000
    0.0  1.0 1.0 0.75  0.693147180560
    0.0  1.0 1.0 0.90  1.609437912434
    0.0  1.0 1.0 0.99  3.912023005428
    1.0  2.0 1.5 0.01 -3.293377493660
    1.0  2.0 1.5 0.10 -1.127793614205
    1.0  2.0 1.5 0.25 -0.042916931074
    1.0  2.0 1.5 0.50  1.000000000000
    1.0  2.0 1.5 0.75  2.042916931074
    1.0  2.0 1.5 0.90  3.127793614205
    1.0  2.0 1.5 0.99  5.293377493660
    -2.0 0.5 0.7 0.01 -6.675320871856
    -2.0 0.5 0.7 0.10 -3.564309791145
    -2.0 0.5 0.7 0.25 -2.582273226038
    -2.0 0.5 0.7 0.50 -2.000000000000
    -2.0 0.5 0.7 0.75 -1.417726773962
    -2.0 0.5 0.7 0.90 -0.435690208855
    -2.0 0.5 0.7 0.99  2.675320871856
    0.0  1.0 2.0 0.01 -1.644976357133
    0.0  1.0 2.0 0.10 -0.906193802437
    0.0  1.0 2.0 0.25 -0.476936276204
    0.0  1.0 2.0 0.50  0.000000000000
    0.0  1.0 2.0 0.75  0.476936276204
    0.0  1.0 2.0 0.90  0.906193802437
    0.0  1.0 2.0 0.99  1.644976357133
    3.0  1.5 5.0 0.01  1.337104467619
    3.0  1.5 5.0 0.10  1.851154137657
    3.0  1.5 5.0 0.25  2.308997603013
    3.0  1.5 5.0 0.50  3.000000000000
    3.0  1.5 5.0 0.75  3.691002396987
    3.0  1.5 5.0 0.90  4.148845862343
    3.0  1.5 5.0 0.99  4.662895532381))

(t/deftest generalized-normal-mv
  (t/are [mu alpha beta mean-v var-v]
      (let [dist (sut/distribution :gnd {:mu mu :alpha alpha :beta beta})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.0  1.0 1.0 0.000000000000 2.000000000000
    1.0  2.0 1.5 1.000000000000 2.953952446487
    -2.0 0.5 0.7 -2.000000000000 2.451238183054
    0.0  1.0 2.0 0.000000000000 0.500000000000
    3.0  1.5 5.0 3.000000000000 0.729862040625))

;; reference values from R's GeneralizedHyperbolic package: dgig/pgig/gigMean/gigVar
;; note: for cdf, R's pgig() DEFAULTS to a spline-interpolation approximation over a
;; 501-point grid which was found (during development of this test) to be inaccurate
;; by up to ~2% at some points (e.g. chi=2,psi=3,lambda=1.5, x=1.5 and x=3, and
;; chi=1,psi=1,lambda=-0.5, x=0.5) - so cdf reference values here are instead
;; R's dgig() integrated directly via integrate(dgig, 0, x, rel.tol=1e-12), which
;; was cross-checked against the closed-form inverse-Gaussian cdf for the
;; lambda=-0.5 case (chi=1,psi=1) and agrees to ~1e-9. icdf reference values are
;; similarly obtained via uniroot() on that same directly-integrated cdf, rather
;; than R's qgig() (which inherits the same spline-based inaccuracy).
;; also note: fastmath's own cdf/icdf are themselves interpolated from a
;; table built once at construction via integrate-pdf with :monotone
;; interpolation, for speed (see the implementation comment in
;; distributions.clj), so cdf/icdf tests below use a looser tolerance than
;; the default 1e-6; empirically the interpolation error is below ~5e-5
;; absolute (icdf) / ~4e-6 absolute (cdf) across the tested parameters.

(t/deftest gig
  (t/testing "registered under both :generalized-inverse-gaussian and :gig keys"
    (let [d1 (sut/distribution :generalized-inverse-gaussian {:chi 2.0 :psi 3.0 :lambda 1.5})
          d2 (sut/distribution :gig {:chi 2.0 :psi 3.0 :lambda 1.5})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "support is (0, +Inf); pdf/cdf are exactly 0 at/below 0"
    (let [dist (sut/distribution :gig {:chi 2.0 :psi 3.0 :lambda 1.5})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))))
  (t/testing "lambda=-1/2 degenerates exactly to the ordinary inverse-gaussian distribution"
    (let [chi 1.0 psi 1.0
          gig (sut/distribution :gig {:chi chi :psi psi :lambda -0.5})
          ig (sut/distribution :inverse-gaussian {:mu (Math/sqrt (/ chi psi)) :lambda chi})]
      (doseq [x [0.1 0.5 1.0 2.0 5.0]]
        (t/is (m/delta-eq (sut/pdf gig x) (sut/pdf ig x)))
        (t/is (m/delta-eq (sut/cdf gig x) (sut/cdf ig x) 1.0e-5)))))
  (t/testing "pdf/cdf are never NaN or throw, at extreme/infinite inputs and far into the tail"
    (doseq [lambda [-2.0 -0.5 0.0 0.5 2.0]
            chi [0.1 1.0 10.0]
            psi [0.1 1.0 10.0]]
      (let [dist (sut/distribution :gig {:chi chi :psi psi :lambda lambda})]
        (doseq [x [##-Inf ##Inf 1e-300 1e300 1e10]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf chi=" chi " psi=" psi " lambda=" lambda " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf chi=" chi " psi=" psi " lambda=" lambda " x=" x)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist 1e300))))))
  (t/are [chi psi lambda vd d vp p]
      (let [dist (sut/distribution :gig {:chi chi :psi psi :lambda lambda})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp) 1.0e-5)))
    2.0 3.0 1.5  0.1 0.000086009808 0.1 0.000000705532
    2.0 3.0 1.5  0.5 0.314638704849 0.5 0.044873356169
    2.0 3.0 1.5  1   0.571348069130 1   0.290408039823
    2.0 3.0 1.5  1.2 0.547753871888 1.2 0.403025631899
    2.0 3.0 1.5  3   0.095963836090 3   0.925339956260
    2.0 3.0 1.5  6   0.001781065303 6   0.998729558937
    1.0 1.0 -0.5 0.1 0.219794800319 0.1 0.004076111321
    1.0 1.0 -0.5 0.5 0.878782578935 0.5 0.364975548173
    1.0 1.0 -0.5 1   0.398942280401 1   0.668102001223
    1.0 1.0 -0.5 2   0.109847822367 2   0.885475425986
    1.0 1.0 -0.5 5   0.007204168934 5   0.990115297400
    1.0 1.0 -0.5 10  0.000219794800 10  0.999649585463
    5.0 0.5 2    0.5 0.000294290652 0.5 0.000019431554
    5.0 0.5 2    2   0.034401734094 2   0.020929835182
    5.0 0.5 2    5   0.086004335236 5   0.222232014744
    5.0 0.5 2   10   0.063278453570 10  0.619023070784
    5.0 0.5 2   20   0.011771626074 20  0.942299600685
    5.0 0.5 2   40   0.000168864147 40  0.999252755438
    0.3 4.0 -2   0.05 10.351780861468 0.05 0.235895202831
    0.3 4.0 -2   0.1  5.247317675776 0.1  0.632142803417
    0.3 4.0 -2   0.3  0.354119858839 0.3  0.957912454010
    0.3 4.0 -2   0.6  0.031193005857 0.6  0.994507482719
    0.3 4.0 -2   1    0.003345837471 1    0.999235674064
    0.3 4.0 -2   2    0.000061009572 2    0.999981433574
    3.0 3.0 0    0.2  0.029486242407 0.2  0.000723179922
    3.0 3.0 0    0.7  0.844127613626 0.7  0.259872878334
    3.0 3.0 0    1    0.716577125198 1    0.500000000000
    3.0 3.0 0    1.5  0.372047217490 1.5  0.768166382782
    3.0 3.0 0    3    0.032326056067 3    0.980534932652
    3.0 3.0 0    6    0.000230553213 6    0.999857438267))

(t/deftest gig-icdf
  (t/are [chi psi lambda p vq]
      (let [dist (sut/distribution :gig {:chi chi :psi psi :lambda lambda})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    2.0 3.0 1.5  0.01 0.344846583619
    2.0 3.0 1.5  0.10 0.642835598903
    2.0 3.0 1.5  0.25 0.929120116594
    2.0 3.0 1.5  0.50 1.384879376059
    2.0 3.0 1.5  0.75 2.020173627054
    2.0 3.0 1.5  0.90 2.770897033190
    2.0 3.0 1.5  0.99 4.510605185006
    1.0 1.0 -0.5 0.01 0.119841240594
    1.0 1.0 -0.5 0.10 0.237624708727
    1.0 1.0 -0.5 0.25 0.379723027460
    1.0 1.0 -0.5 0.50 0.675841305697
    1.0 1.0 -0.5 0.75 1.244059755876
    1.0 1.0 -0.5 0.90 2.143033912958
    1.0 1.0 -0.5 0.99 4.984094843321
    5.0 0.5 2    0.01  1.617814255002
    5.0 0.5 2    0.10  3.466982775327
    5.0 0.5 2    0.25  5.320709672986
    5.0 0.5 2    0.50  8.297721842386
    5.0 0.5 2    0.75 12.428758104200
    5.0 0.5 2    0.90 17.262613747342
    5.0 0.5 2    0.99 28.306349666411
    0.3 4.0 -2   0.01 0.021818256102
    0.3 4.0 -2   0.10 0.036392182172
    0.3 4.0 -2   0.25 0.051364903002
    0.3 4.0 -2   0.50 0.079008251998
    0.3 4.0 -2   0.75 0.128228385945
    0.3 4.0 -2   0.90 0.207012097005
    0.3 4.0 -2   0.99 0.499552767205))

(t/deftest gig-mv
  (t/are [chi psi lambda mean-v var-v]
      (let [dist (sut/distribution :gig {:chi chi :psi psi :lambda lambda})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    2.0 3.0 1.5   1.579795897113 0.803904751986
    1.0 1.0 -0.5  1.000000000000 1.000000000000
    5.0 0.5 2     9.547009312686 33.418724935719
    0.3 4.0 -2    0.107726655315 0.009531640077
    3.0 3.0 0     1.155929879761 0.434446032916))

;; reference values from R's GeneralizedHyperbolic package: dghyp/pghyp/qghyp/ghypMean/ghypVar
;; (mu, delta, alpha, beta, lambda) parameterization, cross-checked live in nREPL.
;; pdf/mean/variance are closed-form (default tolerance); cdf/icdf go through
;; the same integrate-pdf/:monotone interpolation table as [[gig]] - empirically
;; below ~1e-6 absolute error against R here, but an explicit looser tolerance
;; is still used below, consistent with the other interpolation-backed distributions.

(t/deftest gh
  (t/testing "registered under both :generalized-hyperbolic and :gh keys"
    (let [d1 (sut/distribution :generalized-hyperbolic {:mu 1.0 :delta 2.0 :alpha 1.5 :beta -0.5 :lambda 0.5})
          d2 (sut/distribution :gh {:mu 1.0 :delta 2.0 :alpha 1.5 :beta -0.5 :lambda 0.5})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "defaults match R's dghyp defaults (mu=0, delta=1, alpha=1, beta=0, lambda=1)"
    (let [dist (sut/distribution :gh nil)]
      (t/is (m/delta-eq 0.2715716633 (sut/pdf dist 0.5)))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :gh {:mu 0.5 :delta 1.2 :alpha 2.0 :beta 0.7 :lambda 1.0})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (neg? (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (pos? (sut/upper-bound dist)))))
  (t/testing "lambda=-1/2 degenerates exactly to the normal-inverse-gaussian distribution"
    (let [gh (sut/distribution :gh {:mu 1.0 :delta 2.0 :alpha 1.5 :beta -0.5 :lambda -0.5})
          nig (sut/distribution :normal-inverse-gaussian {:alpha 1.5 :beta -0.5 :mu 1.0 :delta 2.0})]
      (doseq [x [-3.0 -1.0 0.0 2.0 4.0]]
        (t/is (m/delta-eq (sut/pdf gh x) (sut/pdf nig x)))
        (t/is (m/delta-eq (sut/cdf gh x) (sut/cdf nig x))))
      (doseq [p [1.0e-6 0.1 0.5 0.9 (- 1.0 1.0e-6)]]
        (t/is (m/delta-eq (sut/icdf gh p) (sut/icdf nig p) 1.0e-4)))
      (t/is (m/delta-eq (sut/mean gh) (sut/mean nig)))
      (t/is (m/delta-eq (sut/variance gh) (sut/variance nig)))
      ;; :normal-inverse-gaussian used to be backed directly by the SSJ
      ;; NormalInverseGaussianDist class, whose cdf is not implemented; sample
      ;; (going through inverseF/icdf) therefore always threw. It is now
      ;; reimplemented as the lambda=-0.5 special case of generalized-hyperbolic
      ;; (numerically-integrated cdf/icdf), so sampling works.
      (t/is (number? (sut/sample nig)))))
  (t/testing "pdf/cdf are never NaN or throw, at extreme/infinite inputs across a range of lambda"
    (doseq [lambda [-3.0 -1.0 -0.5 0.0 0.5 1.0 2.5 5.0]]
      (let [dist (sut/distribution :gh {:mu 0.5 :delta 1.2 :alpha 2.0 :beta 0.7 :lambda lambda})]
        (doseq [x [##-Inf ##Inf 1e-300 1e300 -1e300 1e150 -1e150 1e10]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf lambda=" lambda " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf lambda=" lambda " x=" x)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf))))))
  (t/are [mu delta alpha beta lambda vd d vp p]
      (let [dist (sut/distribution :gh {:mu mu :delta delta :alpha alpha :beta beta :lambda lambda})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp) 1.0e-5)))
    0.0 1.0 1.0 0.0  1.0  -5   0.005069491583 -5 0.005143606904
    0.0 1.0 1.0 0.0  1.0  -1   0.201955319875 -1 0.234335935300
    0.0 1.0 1.0 0.0  1.0   0   0.305594801587  0 0.500000000000
    0.0 1.0 1.0 0.0  1.0   0.5 0.271571663331  0.5 0.646864175300
    0.0 1.0 1.0 0.0  1.0   1   0.201955319875  1  0.765664064700
    0.0 1.0 1.0 0.0  1.0   3   0.035162577802  3  0.963624141900
    0.0 1.0 1.0 0.0  1.0   5   0.005069491583  5  0.994856393100
    0.0 1.0 2.0 0.5  1.0  -5   0.000009730700 -5  0.000003945680
    0.0 1.0 2.0 0.5  1.0  -1   0.114107667510 -1  0.054110589220
    0.0 1.0 2.0 0.5  1.0   0   0.430767964314   0 0.322719773800
    0.0 1.0 2.0 0.5  1.0   0.5 0.436811434412 0.5 0.547947468500
    0.0 1.0 2.0 0.5  1.0   1   0.310176799081   1 0.736934062500
    0.0 1.0 2.0 0.5  1.0   3   0.025559632908   3 0.982069751600
    0.0 1.0 2.0 0.5  1.0   5   0.001444163921   5 0.999016699200
    1.0 2.0 1.5 -0.5 0.5  -5   0.004662186166  -5 0.004638468949
    1.0 2.0 1.5 -0.5 0.5  -1   0.176259386636  -1 0.211239220300
    1.0 2.0 1.5 -0.5 0.5   0   0.290527658585   0 0.447918115000
    1.0 2.0 1.5 -0.5 0.5   0.5 0.305401605012 0.5 0.599038336600
    1.0 2.0 1.5 -0.5 0.5   1   0.264581233073   1 0.743805885300
    1.0 2.0 1.5 -0.5 0.5   3   0.023854114013   3 0.986592126900
    1.0 2.0 1.5 -0.5 0.5   5   0.000598413111   5 0.999694418700
    -1.0 0.5 3.0 1.0 -0.5 -5   0.000000018542  -5 0.000000004276
    -1.0 0.5 3.0 1.0 -0.5 -1   1.089541771746  -1 0.349179969500
    -1.0 0.5 3.0 1.0 -0.5  0   0.125975408510   0 0.955614305100
    -1.0 0.5 3.0 1.0 -0.5  0.5 0.029979255678 0.5 0.988820173600
    -1.0 0.5 3.0 1.0 -0.5  1   0.007734209888   1 0.996981433000
    -1.0 0.5 3.0 1.0 -0.5  3   0.000055273799   3 0.999976208000
    -1.0 0.5 3.0 1.0 -0.5  5   0.000000566621   5 0.999999745400
    2.0 1.5 1.2 0.8 2.0   -5   0.000000863944  -5 0.000000464208
    2.0 1.5 1.2 0.8 2.0   -1   0.001067570091  -1 0.000625881361
    2.0 1.5 1.2 0.8 2.0    0   0.005270851056   0 0.003288045743
    2.0 1.5 1.2 0.8 2.0    0.5 0.010979348165 0.5 0.007195481132
    2.0 1.5 1.2 0.8 2.0    1   0.021417154778   1 0.015056206170
    2.0 1.5 1.2 0.8 2.0    3   0.106079862053   3 0.138640670300
    2.0 1.5 1.2 0.8 2.0    5   0.129720887446   5 0.389570227000))

(t/deftest gh-icdf
  (t/are [mu delta alpha beta lambda p vq]
      (let [dist (sut/distribution :gh {:mu mu :delta delta :alpha alpha :beta beta :lambda lambda})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    0.0 1.0 1.0 0.0  1.0  0.05 -2.669835381000
    0.0 1.0 1.0 0.0  1.0  0.25 -0.924464417500
    0.0 1.0 1.0 0.0  1.0  0.50  0.000000000000
    0.0 1.0 1.0 0.0  1.0  0.75  0.924464417500
    0.0 1.0 1.0 0.0  1.0  0.95  2.669835381000
    0.0 1.0 2.0 0.5  1.0  0.05 -1.037332207000
    0.0 1.0 2.0 0.5  1.0  0.25 -0.178270415100
    0.0 1.0 2.0 0.5  1.0  0.50  0.392238684700
    0.0 1.0 2.0 0.5  1.0  0.75  1.042966006000
    0.0 1.0 2.0 0.5  1.0  0.95  2.271338780000
    1.0 2.0 1.5 -0.5 0.5  0.05 -2.589666168000
    1.0 2.0 1.5 -0.5 0.5  0.25 -0.794663677600
    1.0 2.0 1.5 -0.5 0.5  0.50  0.175684797400
    1.0 2.0 1.5 -0.5 0.5  0.75  1.023551004000
    1.0 2.0 1.5 -0.5 0.5  0.95  2.222961341000
    -1.0 0.5 3.0 1.0 -0.5 0.05 -1.459721739000
    -1.0 0.5 3.0 1.0 -0.5 0.25 -1.097380091000
    -1.0 0.5 3.0 1.0 -0.5 0.50 -0.865648523700
    -1.0 0.5 3.0 1.0 -0.5 0.75 -0.597948455600
    -1.0 0.5 3.0 1.0 -0.5 0.95 -0.041865152320
    2.0 1.5 1.2 0.8 2.0   0.05  1.928518803000
    2.0 1.5 1.2 0.8 2.0   0.25  3.933697962000
    2.0 1.5 1.2 0.8 2.0   0.50  5.884314119000
    2.0 1.5 1.2 0.8 2.0   0.75  8.510485776000
    2.0 1.5 1.2 0.8 2.0   0.95 13.710063820000))

(t/deftest gh-mv
  (t/are [mu delta alpha beta lambda mean-v var-v]
      (let [dist (sut/distribution :gh {:mu mu :delta delta :alpha alpha :beta beta :lambda lambda})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.0 1.0 1.0 0.0  1.0   0.000000000000  2.699483936000
    0.0 1.0 2.0 0.5  1.0   0.475739211500  1.045544872000
    1.0 2.0 1.5 -0.5 0.5   0.042893218810  2.215990258000
    -1.0 0.5 3.0 1.0 -0.5 -0.823223304700  0.198873782200
    2.0 1.5 1.2 0.8 2.0    6.598150838000 13.993602450000))

;; reference values from R's `GeneralizedHyperbolic` package: dnig/pnig/qnig
;; (x/q/p, mu, delta, alpha, beta) - the same (mu, delta, alpha, beta)
;; parameterization used here. mean/variance reference values from the
;; closed-form gamma=sqrt(alpha^2-beta^2); mean=mu+delta*beta/gamma;
;; variance=delta*alpha^2/gamma^3.
;;
;; pdf is closed-form (no numerical integration involved) and matches R very
;; tightly; cdf/icdf go through fastmath's own numerical integration
;; (`integrate-pdf`, via the shared generalized-hyperbolic-core machinery) and
;; so carry a looser, but still small, tolerance - especially for the
;; alpha=5/beta=4.9 combo below, deliberately chosen close to the alpha=beta
;; boundary (heaviest skew) to stress that numerical path.

(t/deftest normal-inverse-gaussian
  (t/testing "defaults are alpha=1, beta=0, mu=0, delta=1 (standard/symmetric NIG)"
    (let [dist (sut/distribution :normal-inverse-gaussian nil)]
      (t/is (m/delta-eq 0.5208038299916704 (sut/pdf dist 0.0)))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :normal-inverse-gaussian {:mu 1.0 :delta 2.0 :alpha 1.5 :beta -0.5})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (neg? (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (pos? (sut/upper-bound dist)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, at extreme/infinite inputs across a range of alpha/beta"
    (doseq [[alpha beta] [[1.0 0.0] [1.5 -0.5] [3.0 1.0] [5.0 4.9] [0.2 0.15] [50.0 -49.5]]]
      (let [dist (sut/distribution :normal-inverse-gaussian {:alpha alpha :beta beta})]
        (doseq [x [##-Inf ##Inf 1e-300 1e300 -1e300 1e150 -1e150 1e10]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf alpha=" alpha " beta=" beta " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf alpha=" alpha " beta=" beta " x=" x)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf))))))
  (t/are [mu delta alpha beta x vd vp]
      (let [dist (sut/distribution :normal-inverse-gaussian {:mu mu :delta delta :alpha alpha :beta beta})]
        (and (m/delta-eq vd (sut/pdf dist x) 1.0e-9)
             (m/delta-eq vp (sut/cdf dist x) 1.0e-4)))
    0.0 1.0 1.0 0.0 -5     0.000614810161509777 0.000492467344744954
    0.0 1.0 1.0 0.0 -2     0.0398684291217511    0.0272228569059907
    0.0 1.0 1.0 0.0 -1     0.192235012744407     0.124034777509506
    0.0 1.0 1.0 0.0  0     0.52080382999167      0.499999999984612
    0.0 1.0 1.0 0.0  0.5   0.383145915640741     0.735169095857493
    0.0 1.0 1.0 0.0  1     0.192235012744407     0.875965222490494
    0.0 1.0 1.0 0.0  2     0.0398684291217511    0.972777143094009
    0.0 1.0 1.0 0.0  5     0.000614810161509777  0.999507532655255
    1.0 2.0 1.5 -0.5 -5    0.00164418727983464   0.001435230081675
    1.0 2.0 1.5 -0.5 -2    0.0516691035686478    0.0453911360227543
    1.0 2.0 1.5 -0.5 -1    0.147017410121953     0.137731165281906
    1.0 2.0 1.5 -0.5  0    0.314285019535692     0.366554083220295
    1.0 2.0 1.5 -0.5  0.5  0.361872442613151     0.538425873179007
    1.0 2.0 1.5 -0.5  1    0.324389499299296     0.713825618030112
    1.0 2.0 1.5 -0.5  2    0.115618997355346     0.93373440051395
    1.0 2.0 1.5 -0.5  5    0.000304317260197264  0.999858816563824
    -2.0 0.5 3.0 1.0 -5    1.5119056517774e-06   3.41243551818395e-07
    -2.0 0.5 3.0 1.0 -2    1.08954177174606      0.349179970606016
    -2.0 0.5 3.0 1.0 -1    0.125975408509954     0.955614305046348
    -2.0 0.5 3.0 1.0  0    0.00773420988794332   0.996981433359096
    -2.0 0.5 3.0 1.0  0.5  0.00212325089400032   0.999142358843805
    -2.0 0.5 3.0 1.0  1    0.000609946272971403  0.999747056249221
    -2.0 0.5 3.0 1.0  2    5.52737986202059e-05  0.999976208054767
    -2.0 0.5 3.0 1.0  5    6.1308584355612e-08   0.999999972085304
    0.0 1.0 5.0 4.9 -5     4.12054200199846e-23  4.08280463315396e-24
    0.0 1.0 5.0 4.9 -2     5.76236025241476e-10  5.75875450676322e-11
    0.0 1.0 5.0 4.9 -1     9.53556958076887e-6   1.01846159623679e-06
    0.0 1.0 5.0 4.9  0     0.0174106388248045    0.00297961443641972
    0.0 1.0 5.0 4.9  0.5   0.0939581892941343    0.0284544189868986
    0.0 1.0 5.0 4.9  1     0.171962029561146     0.0969236465084977
    0.0 1.0 5.0 4.9  2     0.187401149900837     0.287800733735118
    0.0 1.0 5.0 4.9  5     0.0785930166433385    0.670432529348062))

(t/deftest normal-inverse-gaussian-icdf
  (t/are [mu delta alpha beta p vq]
      (let [dist (sut/distribution :normal-inverse-gaussian {:mu mu :delta delta :alpha alpha :beta beta})]
        ;; the alpha=5/beta=4.9 combo's far-left tail (p=0.001) is the one row
        ;; needing the full tolerance below; every other row matches much more
        ;; tightly (~1e-5 or better).
        (m/delta-eq vq (sut/icdf dist p) 2.0e-3))
    0.0 1.0 1.0 0.0  0.001  -4.43810263029802
    0.0 1.0 1.0 0.0  0.05   -1.59135402609246
    0.0 1.0 1.0 0.0  0.25   -0.539591454418318
    0.0 1.0 1.0 0.0  0.5     0
    0.0 1.0 1.0 0.0  0.75    0.539591454418318
    0.0 1.0 1.0 0.0  0.95    1.59135402609246
    0.0 1.0 1.0 0.0  0.999   4.43810263029804
    1.0 2.0 1.5 -0.5 0.001  -5.3158077717406
    1.0 2.0 1.5 -0.5 0.05   -1.91491843937776
    1.0 2.0 1.5 -0.5 0.25   -0.417440482954573
    1.0 2.0 1.5 -0.5 0.5     0.393298744917891
    1.0 2.0 1.5 -0.5 0.75    1.11504331133029
    1.0 2.0 1.5 -0.5 0.95    2.15835367302664
    1.0 2.0 1.5 -0.5 0.999   4.08839851952259
    -2.0 0.5 3.0 1.0 0.001  -3.25802450264978
    -2.0 0.5 3.0 1.0 0.05   -2.45974433253412
    -2.0 0.5 3.0 1.0 0.25   -2.09738180379992
    -2.0 0.5 3.0 1.0 0.5    -1.86565374422419
    -2.0 0.5 3.0 1.0 0.75   -1.59794944503895
    -2.0 0.5 3.0 1.0 0.95   -1.04185634102317
    -2.0 0.5 3.0 1.0 0.999   0.438092163214403
    0.0 1.0 5.0 4.9  0.001  -0.172918689201281
    0.0 1.0 5.0 4.9  0.05    0.6923425449097
    0.0 1.0 5.0 4.9  0.25    1.80158384370533
    0.0 1.0 5.0 4.9  0.5     3.33217890755459
    0.0 1.0 5.0 4.9  0.75    6.19652646142522
    0.0 1.0 5.0 4.9  0.95   14.5916871332332
    0.0 1.0 5.0 4.9  0.999  41.7451002030271))

(t/deftest normal-inverse-gaussian-mv
  (t/are [mu delta alpha beta mean-v var-v]
      (let [dist (sut/distribution :normal-inverse-gaussian {:mu mu :delta delta :alpha alpha :beta beta})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.0 1.0 1.0 0.0    0                1
    1.0 2.0 1.5 -0.5   0.292893218813453 1.59099025766973
    -2.0 0.5 3.0 1.0  -1.82322330470336  0.198873782208716
    0.0 1.0 5.0 4.9    4.92468529477015 25.3797428095763))

;; reference values from R's bayesmeta package: dhalflogistic/phalflogistic/
;; qhalflogistic/ehalflogistic/vhalflogistic (scale parameterization).

(t/deftest half-logistic
  (t/testing "defaults are scale=1"
    (let [dist (sut/distribution :half-logistic nil)]
      (t/is (m/delta-eq 0.470007424403 (sut/pdf dist 0.5)))))
  (t/testing "support is [0, +Inf); pdf/cdf are exactly 0 below 0"
    (let [dist (sut/distribution :half-logistic {:scale 2.0})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))))
  (t/testing "matches the alpha=1 special case of generalized-half-logistic (lambda=1/scale)"
    (doseq [scale [0.5 1.0 2.5]]
      (let [hl (sut/distribution :half-logistic {:scale scale})
            ghl (sut/distribution :ghl {:alpha 1.0 :lambda (/ 1.0 scale)})]
        (doseq [x [0.0 0.3 1.0 2.5 5.0]]
          (t/is (m/delta-eq (sut/pdf hl x) (sut/pdf ghl x)))
          (t/is (m/delta-eq (sut/cdf hl x) (sut/cdf ghl x))))
        (t/is (m/delta-eq (sut/mean hl) (sut/mean ghl) 1.0e-6))
        (t/is (m/delta-eq (sut/variance hl) (sut/variance ghl) 1.0e-6)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme/infinite/negative-zero inputs"
    (doseq [scale [0.01 0.5 1.0 2.0 100.0]]
      (let [dist (sut/distribution :half-logistic {:scale scale})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf scale=" scale " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf scale=" scale " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf scale=" scale " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [scale vd d vp p]
      (let [dist (sut/distribution :half-logistic {:scale scale})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 0    0.500000000000 0   0.000000000000
    1.0 0.2  0.495033145424 0.2 0.099667994625
    1.0 0.5  0.470007424403 0.5 0.244918662404
    1.0 1    0.393223866483 1   0.462117157260
    1.0 2    0.209987170807 2   0.761594155956
    1.0 5    0.013296113342 5   0.986614298151
    1.0 10   0.000090791615 10  0.999909204263
    0.5 0    1.000000000000 0   0.000000000000
    0.5 0.2  0.961042982966 0.2 0.197375320225
    0.5 0.5  0.786447732966 0.5 0.462117157260
    0.5 1    0.419974341614 1   0.761594155956
    0.5 2    0.070650824853 2   0.964027580076
    0.5 5    0.000181583231 5   0.999909204263
    2.0 0    0.250000000000 0   0.000000000000
    2.0 0.2  0.249376040193 0.2 0.049958374958
    2.0 0.5  0.246134082738 0.5 0.124353001772
    2.0 1    0.235003712202 1   0.244918662404
    2.0 2    0.196611933241 2   0.462117157260
    2.0 5    0.070103716545 5   0.848283639958
    2.0 10   0.006648056671 10  0.986614298151
    1.5 0    0.333333333333 0   0.000000000000
    1.5 0.2  0.331856230397 0.2 0.066568076502
    1.5 0.5  0.324242881340 0.5 0.165140412925
    1.5 1    0.298876519868 1   0.321512737532
    1.5 2    0.220121346204 2   0.582782945348
    1.5 5    0.044344965549 5   0.931109608668
    3.0 0    0.166666666667 0   0.000000000000
    3.0 0.2  0.166481618569 0.2 0.033320993139
    3.0 0.5  0.165514596617 0.5 0.083140966434
    3.0 1    0.162121440670 1   0.165140412925
    3.0 2    0.149438259934 2   0.321512737532
    3.0 5    0.089086474930 5   0.682261790238))

(t/deftest half-logistic-icdf
  (t/are [scale p vq]
      (let [dist (sut/distribution :half-logistic {:scale scale})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 0.05 0.100083458557
    1.0 0.25 0.510825623766
    1.0 0.50 1.098612288668
    1.0 0.75 1.945910149055
    1.0 0.95 3.663561646130
    0.5 0.05 0.050041729278
    0.5 0.25 0.255412811883
    0.5 0.50 0.549306144334
    0.5 0.75 0.972955074528
    0.5 0.95 1.831780823065
    2.0 0.05 0.200166917114
    2.0 0.25 1.021651247532
    2.0 0.50 2.197224577336
    2.0 0.75 3.891820298111
    2.0 0.95 7.327123292259
    1.5 0.05 0.150125187835
    1.5 0.25 0.766238435649
    1.5 0.50 1.647918433002
    1.5 0.75 2.918865223583
    1.5 0.95 5.495342469194
    3.0 0.05 0.300250375671
    3.0 0.25 1.532476871298
    3.0 0.50 3.295836866004
    3.0 0.75 5.837730447166
    3.0 0.95 10.990684938389))

(t/deftest half-logistic-mv
  (t/are [scale mean-v var-v]
      (let [dist (sut/distribution :half-logistic {:scale scale})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    1.0 1.386294361120 1.368056078024
    0.5 0.693147180560 0.342014019506
    2.0 2.772588722240 5.472224312095
    1.5 2.079441541680 3.078126175553
    3.0 4.158883083360 12.312504702213))

;; reference values from R's base `stats` package: df/pf/qf(x, df1, df2, ncp).
;; mean/variance closed-form formulas cross-checked against R's own
;; high-precision numerical integration (rel.tol=1e-10) of x*df(...) /
;; x^2*df(...) over [0, Inf).

(t/deftest f-noncentral
  (t/testing "defaults are df1=1, df2=1, ncp=1"
    (let [dist (sut/distribution :f-noncentral nil)]
      (t/is (m/delta-eq 0.2499108 (sut/pdf dist 0.5) 1.0e-6))))
  (t/testing "support is [0, +Inf); pdf/cdf are exactly 0 below 0"
    (let [dist (sut/distribution :f-noncentral {:df1 3.0 :df2 5.0 :ncp 2.0})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))))
  (t/testing "pdf(0) follows the central-F boundary convention based on df1: 0 (df1>2), e^(-ncp/2) (df1=2), +Inf (df1<2)"
    (t/is (m/delta-eq 0.0 (sut/pdf (sut/distribution :f-noncentral {:df1 3.0 :df2 5.0 :ncp 2.0}) 0.0)))
    (t/is (m/delta-eq (Math/exp -5.0) (sut/pdf (sut/distribution :f-noncentral {:df1 2.0 :df2 8.0 :ncp 10.0}) 0.0)))
    (t/is (Double/isInfinite (sut/pdf (sut/distribution :f-noncentral {:df1 1.0 :df2 10.0 :ncp 5.0}) 0.0))))
  (t/testing "ncp=0 reduces exactly to the central f distribution with the same df1/df2"
    (let [fnc (sut/distribution :f-noncentral {:df1 5.0 :df2 20.0 :ncp 0.0})
          fc (sut/distribution :f {:numerator-degrees-of-freedom 5.0 :denominator-degrees-of-freedom 20.0})]
      (doseq [x [0.1 0.5 1.0 2.0 5.0 10.0]]
        (t/is (m/delta-eq (sut/pdf fnc x) (sut/pdf fc x)))
        (t/is (m/delta-eq (sut/cdf fnc x) (sut/cdf fc x))))
      (t/is (m/delta-eq (sut/mean fnc) (sut/mean fc)))
      (t/is (m/delta-eq (sut/variance fnc) (sut/variance fc)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme/infinite/negative-zero inputs"
    (doseq [df1 [0.5 1.0 2.0 3.0 10.0]
            df2 [1.0 3.0 5.0 20.0]
            ncp [0.0 1.0 5.0 20.0]]
      (let [dist (sut/distribution :f-noncentral {:df1 df1 :df2 df2 :ncp ncp})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf df1=" df1 " df2=" df2 " ncp=" ncp " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf df1=" df1 " df2=" df2 " ncp=" ncp " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf df1=" df1 " df2=" df2 " ncp=" ncp " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [df1 df2 ncp vd d vp p]
      (let [dist (sut/distribution :f-noncentral {:df1 df1 :df2 df2 :ncp ncp})]
        (and (m/delta-eq d (sut/pdf dist vd) 1.0e-8)
             (m/delta-eq p (sut/cdf dist vp) 1.0e-8)))
    3.0 5.0 2.0   0    0.000000000000 0   0.000000000000
    3.0 5.0 2.0   0.5  0.383289609397 0.5 0.156955143502
    3.0 5.0 2.0   1    0.327982758665 1   0.337743799128
    3.0 5.0 2.0   2    0.186995588358 2   0.590539309716
    3.0 5.0 2.0   5    0.040915612366 5   0.867755079741
    3.0 5.0 2.0   10   0.007595516466 10  0.960691337583
    3.0 5.0 2.0   20   0.001023647824 20  0.990637883942
    1.0 10.0 5.0  0.5  0.092929175169 0.5 0.061377192704
    1.0 10.0 5.0  1    0.095556881152 1   0.108334208577
    1.0 10.0 5.0  2    0.099172816782 2   0.206363879018
    1.0 10.0 5.0  5    0.077958962955 5   0.478890464179
    1.0 10.0 5.0  10   0.036389788937 10  0.755426066631
    1.0 10.0 5.0  20   0.007700257160 20  0.937222483692
    5.0 20.0 0.0  0    0.000000000000 0   0.000000000000
    5.0 20.0 0.0  0.5  0.718990267178 0.5 0.227395614209
    5.0 20.0 0.0  1    0.544878125182 1   0.556974815315
    5.0 20.0 0.0  2    0.157789751008 2   0.877492755318
    5.0 20.0 0.0  5    0.003925075958 5   0.996069580076
    5.0 20.0 0.0  10   0.000044342889 10  0.999934478169
    2.0 8.0 10.0  0    0.006737946999 0   0.000000000000
    2.0 8.0 10.0  0.5  0.027803804656 0.5 0.008339016296
    2.0 8.0 10.0  1    0.052264530299 1   0.028357493084
    2.0 8.0 10.0  2    0.091174725736 2   0.101813065751
    2.0 8.0 10.0  5    0.097776334101 5   0.414514982035
    2.0 8.0 10.0  10   0.041384947150 10  0.751135577845
    10.0 15.0 3.0 0    0.000000000000 0   0.000000000000
    10.0 15.0 3.0 0.5  0.377148023147 0.5 0.062209324523
    10.0 15.0 3.0 1    0.616860389219 1   0.338816177006
    10.0 15.0 3.0 2    0.257424311786 2   0.782741442022
    10.0 15.0 3.0 5    0.008159776484 5   0.991004988066
    10.0 15.0 3.0 10   0.000141802314 10  0.999754510194))

(t/deftest f-noncentral-icdf
  (t/are [df1 df2 ncp p vq]
      (let [dist (sut/distribution :f-noncentral {:df1 df1 :df2 df2 :ncp ncp})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    3.0 5.0 2.0   0.05 0.209560327500
    3.0 5.0 2.0   0.25 0.747184884000
    3.0 5.0 2.0   0.50 1.573603203000
    3.0 5.0 2.0   0.75 3.166403204000
    3.0 5.0 2.0   0.95 8.812091486000
    1.0 10.0 5.0  0.05 0.378059265200
    1.0 10.0 5.0  0.25 2.441450528000
    1.0 10.0 5.0  0.50 5.275645834000
    1.0 10.0 5.0  0.75 9.852668498000
    1.0 10.0 5.0  0.95 21.892046830000
    5.0 20.0 0.0  0.05 0.219388142000
    5.0 20.0 0.0  0.25 0.531356426300
    5.0 20.0 0.0  0.50 0.900376483300
    5.0 20.0 0.0  0.75 1.449952290000
    5.0 20.0 0.0  0.95 2.710889837000
    2.0 8.0 10.0  0.05 1.357415885000
    2.0 8.0 10.0  0.25 3.432295562000
    2.0 8.0 10.0  0.50 5.929307819000
    2.0 8.0 10.0  0.75 9.972631607000
    2.0 8.0 10.0  0.95 21.238083830000
    10.0 15.0 3.0 0.05 0.465815679900
    10.0 15.0 3.0 0.25 0.856630186000
    10.0 15.0 3.0 0.50 1.276270549000
    10.0 15.0 3.0 0.75 1.881429336000
    10.0 15.0 3.0 0.95 3.278755496000))

(t/deftest f-noncentral-mv
  (t/are [df1 df2 ncp mean-v var-v]
      (let [dist (sut/distribution :f-noncentral {:df1 df1 :df2 df2 :ncp ncp})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    3.0 5.0 2.0   2.777777778000 28.395061730000
    1.0 10.0 5.0  7.500000000000 64.583333330000
    5.0 20.0 0.0  1.111111111000  0.709876543200
    2.0 8.0 10.0  8.000000000000 61.333333330000
    10.0 15.0 3.0 1.500000000000  0.912587412600)
  (t/testing "mean/variance are +Inf outside their validity ranges (df2<=2 / df2<=4)"
    (t/is (Double/isInfinite (sut/mean (sut/distribution :f-noncentral {:df1 3.0 :df2 2.0 :ncp 1.0}))))
    (t/is (Double/isInfinite (sut/variance (sut/distribution :f-noncentral {:df1 3.0 :df2 4.0 :ncp 1.0}))))
    (t/is (m/delta-eq 3.0 (sut/mean (sut/distribution :f-noncentral {:df1 1.0 :df2 3.0 :ncp 0.0}))))))

;; reference values from R's base `stats` package: dt/pt/qt(x, df, ncp).
;; mean/variance closed-form formulas (mean = ncp*sqrt(df/2)*Gamma((df-1)/2)/Gamma(df/2),
;; variance = df*(1+ncp^2)/(df-2) - ncp^2*df/2*(Gamma((df-1)/2)/Gamma(df/2))^2)
;; cross-checked against R's own high-precision `integrate()` (rel.tol=1e-10)
;; of x*dt(...) / x^2*dt(...) over (-Inf, Inf).

(t/deftest t-noncentral
  (t/testing "defaults are df=1, ncp=1"
    (let [dist (sut/distribution :t-noncentral nil)]
      (t/is (m/delta-eq 0.193064705260 (sut/pdf dist 0.0) 1.0e-6))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :t-noncentral {:df 5.0 :ncp -3.0})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (neg? (sut/lower-bound dist)))
      (t/is (pos? (sut/upper-bound dist)))))
  (t/testing "ncp=0 reduces exactly to the central t distribution with the same df"
    (let [tnc (sut/distribution :t-noncentral {:df 7.0 :ncp 0.0})
          tc (sut/distribution :t {:degrees-of-freedom 7.0})]
      (doseq [x [-3.0 -1.0 -0.1 0.0 0.1 1.0 3.0]]
        (t/is (m/delta-eq (sut/pdf tnc x) (sut/pdf tc x)))
        (t/is (m/delta-eq (sut/cdf tnc x) (sut/cdf tc x))))
      (t/is (m/delta-eq (sut/mean tnc) (sut/mean tc)))
      (t/is (m/delta-eq (sut/variance tnc) (sut/variance tc)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme/infinite/negative-zero inputs"
    (doseq [df [0.5 1.0 2.0 5.0 30.0]
            ncp [-20.0 -1.0 0.0 1.0 20.0]]
      (let [dist (sut/distribution :t-noncentral {:df df :ncp ncp})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf df=" df " ncp=" ncp " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf df=" df " ncp=" ncp " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf df=" df " ncp=" ncp " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [df ncp vd d vp p]
      (let [dist (sut/distribution :t-noncentral {:df df :ncp ncp})]
        (and (m/delta-eq d (sut/pdf dist vd) 1.0e-8)
             (m/delta-eq p (sut/cdf dist vp) 1.0e-8)))
    3.0 2.0    -2   0.000845332062 -2  0.000687152986
    3.0 2.0    -1   0.005278237809 -1  0.003005247595
    3.0 2.0    0    0.049742834812 0   0.022750131948
    3.0 2.0    0.5  0.129665387995 0.5 0.065478985565
    3.0 2.0    1    0.236548732514 1   0.157349433970
    3.0 2.0    2    0.287783357806 2   0.443075782218
    3.0 2.0    5    0.053038652175 5   0.889052067438
    1.0 1.5    -2   0.005306351608 -2  0.011315052843
    1.0 1.5    -1   0.016849979226 -1  0.020857766993
    1.0 1.5    0    0.103340089934 0   0.066807201269
    1.0 1.5    0.5  0.213035559311 0.5 0.145889613765
    1.0 1.5    1    0.257949130174 1   0.267986599354
    1.0 1.5    2    0.176264596771 2   0.491019901518
    1.0 1.5    5    0.044146571506 5   0.763973691912
    10.0 0.0   -2   0.061145766321 -2  0.036694017385
    10.0 0.0   -1   0.230361989229 -1  0.170446566151
    10.0 0.0   0    0.389108383966 0   0.500000000000
    10.0 0.0   0.5  0.339695136352 0.5 0.686053197129
    10.0 0.0   1    0.230361989229 1   0.829553433849
    10.0 0.0   2    0.061145766321 2   0.963305982615
    10.0 0.0   5    0.000396001056 5   0.999731333199
    5.0 -3.0   -2   0.237736447074 -2  0.825197284442
    5.0 -3.0   -1   0.064128630253 -1  0.974632761821
    5.0 -3.0   0    0.004217049403 0   0.998650101968
    5.0 -3.0   0.5  0.000862309322 0.5 0.999706554724
    5.0 -3.0   1    0.000190326894 1   0.999926631614
    5.0 -3.0   2    0.000014965160 2   0.999992200364
    5.0 -3.0   5    0.000000134561 5   0.999999857991
    20.0 4.0   -1   0.000001824240 -1  0.000000412176
    20.0 4.0   0    0.000132168446 0   0.000031671242
    20.0 4.0   0.5  0.000912178865 0.5 0.000236768793
    20.0 4.0   1    0.005007705315 1   0.001464367582
    20.0 4.0   2    0.063767650558 2   0.026819856758
    20.0 4.0   5    0.221154356653 5   0.768687208932))

(t/deftest t-noncentral-icdf
  (t/are [df ncp p vq]
      (let [dist (sut/distribution :t-noncentral {:df df :ncp ncp})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    3.0 2.0  0.05 0.366968979429
    3.0 2.0  0.25 1.349720333326
    3.0 2.0  0.50 2.203826658048
    3.0 2.0  0.75 3.437835721256
    3.0 2.0  0.95 6.852347517156
    1.0 1.5  0.05 -0.195195157360
    1.0 1.5  0.25 0.930288396007
    1.0 1.5  0.50 2.051650979949
    1.0 1.5  0.75 4.700539855855
    1.0 1.5  0.95 24.368899333029
    10.0 0.0 0.05 -1.812461122812
    10.0 0.0 0.25 -0.699812061312
    10.0 0.0 0.50 0.000000000000
    10.0 0.0 0.75 0.699812061312
    10.0 0.0 0.95 1.812461122812
    5.0 -3.0 0.05 -7.099888375581
    5.0 -3.0 0.25 -4.373411292391
    5.0 -3.0 0.50 -3.183251821560
    5.0 -3.0 0.75 -2.293281791745
    5.0 -3.0 0.95 -1.287820536981
    20.0 4.0 0.05 2.279291596400
    20.0 4.0 0.25 3.283428062499
    20.0 4.0 0.50 4.055370185591
    20.0 4.0 0.75 4.917795541343
    20.0 4.0 0.95 6.385579639500))

(t/deftest t-noncentral-mv
  (t/are [df ncp mean-v var-v]
      (let [dist (sut/distribution :t-noncentral {:df df :ncp ncp})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    3.0 2.0   2.763953195771 7.360562731589
    10.0 0.0  0.000000000000 1.250000000000
    5.0 -3.0 -3.568248232306 3.934271219315
    20.0 4.0  4.158243910847 1.597896466793)
  (t/testing "mean is NaN for df<=1 and variance is NaN for df<=2 (moments genuinely do not exist)"
    (t/is (Double/isNaN (sut/mean (sut/distribution :t-noncentral {:df 1.0 :ncp 1.5}))))
    (t/is (Double/isNaN (sut/variance (sut/distribution :t-noncentral {:df 2.0 :ncp 1.5}))))
    (t/is (m/delta-eq 2.658680776358 (sut/mean (sut/distribution :t-noncentral {:df 2.0 :ncp 1.5})) 1.0e-6))))

;; reference values from R's base `stats` package: dbeta/pbeta/qbeta(x, shape1, shape2, ncp).
;; mean/variance cross-checked against R's own high-precision `integrate()`
;; (default tolerance) of x*dbeta(...) / x^2*dbeta(...) over [0,1].

(t/deftest beta-noncentral
  (t/testing "defaults are alpha=2, beta=2, ncp=1"
    (let [dist (sut/distribution :beta-noncentral nil)]
      (t/is (m/delta-eq 1.472420230494 (sut/pdf dist 0.5)))))
  (t/testing "support is [0, 1]; pdf/cdf are exactly 0 below 0 and 1 above 1"
    (let [dist (sut/distribution :beta-noncentral {:alpha 2.0 :beta 3.0 :ncp 4.0})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (m/delta-eq 1.0 (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 2.0)))
      (t/is (m/delta-eq 1.0 (sut/cdf dist 2.0)))))
  (t/testing "pdf(0)/pdf(1) follow a boundary convention based on alpha/beta: 0 (shape>1), the finite k=0-term limit (shape=1), +Inf (shape<1)"
    (t/is (m/delta-eq 0.0 (sut/pdf (sut/distribution :beta-noncentral {:alpha 2.0 :beta 3.0 :ncp 4.0}) 0.0)))
    (t/is (m/delta-eq 1.103638323514 (sut/pdf (sut/distribution :beta-noncentral {:alpha 1.0 :beta 3.0 :ncp 2.0}) 0.0)))
    (t/is (Double/isInfinite (sut/pdf (sut/distribution :beta-noncentral {:alpha 0.5 :beta 3.0 :ncp 2.0}) 0.0)))
    (t/is (m/delta-eq 0.0 (sut/pdf (sut/distribution :beta-noncentral {:alpha 2.0 :beta 3.0 :ncp 4.0}) 1.0)))
    (t/is (m/delta-eq 4.0 (sut/pdf (sut/distribution :beta-noncentral {:alpha 3.0 :beta 1.0 :ncp 2.0}) 1.0)))
    (t/is (Double/isInfinite (sut/pdf (sut/distribution :beta-noncentral {:alpha 3.0 :beta 0.5 :ncp 2.0}) 1.0))))
  (t/testing "ncp=0 reduces exactly to the central beta distribution with the same alpha/beta"
    (let [bnc (sut/distribution :beta-noncentral {:alpha 2.0 :beta 5.0 :ncp 0.0})
          bc (sut/distribution :beta {:alpha 2.0 :beta 5.0})]
      (doseq [x [0.0 0.1 0.3 0.5 0.7 0.9 1.0]]
        (t/is (m/delta-eq (sut/pdf bnc x) (sut/pdf bc x)))
        (t/is (m/delta-eq (sut/cdf bnc x) (sut/cdf bc x))))
      (t/is (m/delta-eq (sut/mean bnc) (sut/mean bc)))
      (t/is (m/delta-eq (sut/variance bnc) (sut/variance bc)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme shape/ncp and boundary/extreme x/p inputs"
    (doseq [alpha [0.1 0.5 1.0 2.0 50.0]
            beta [0.1 0.5 1.0 2.0 50.0]
            ncp [0.0 1.0 50.0 1000.0]]
      (let [dist (sut/distribution :beta-noncentral {:alpha alpha :beta beta :ncp ncp})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0 1.0 0.5]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf alpha=" alpha " beta=" beta " ncp=" ncp " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf alpha=" alpha " beta=" beta " ncp=" ncp " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf alpha=" alpha " beta=" beta " ncp=" ncp " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [alpha beta ncp vd d vp p]
      (let [dist (sut/distribution :beta-noncentral {:alpha alpha :beta beta :ncp ncp})]
        (and (m/delta-eq d (sut/pdf dist vd) 1.0e-8)
             (m/delta-eq p (sut/cdf dist vp) 1.0e-8)))
    2.0 3.0 4.0   0.1  0.212138642384 0.1  0.009772800708
    2.0 3.0 4.0   0.25 0.707958877268 0.25 0.076374141601
    2.0 3.0 4.0   0.5  1.678449950345 0.5  0.379375673481
    2.0 3.0 4.0   0.75 1.540611568274 0.75 0.823614142480
    2.0 3.0 4.0   0.9  0.491896711372 0.9  0.981627060849
    0.5 3.0 2.0   0.1  1.601682094812 0.1  0.253190173894
    0.5 3.0 2.0   0.25 1.374199180289 0.25 0.472827848131
    0.5 3.0 2.0   0.5  1.018594613391 0.5  0.777348520663
    0.5 3.0 2.0   0.75 0.420220981988 0.75 0.960316914740
    0.5 3.0 2.0   0.9  0.089674619779 0.9  0.996864733995
    3.0 0.5 2.0   0.1  0.004084317671 0.1  0.000130527469
    3.0 0.5 2.0   0.25 0.033270832372 0.25 0.002488852093
    3.0 0.5 2.0   0.5  0.217292952493 0.5  0.028683471946
    3.0 0.5 2.0   0.75 0.919587280009 0.75 0.151148834227
    3.0 0.5 2.0   0.9  2.481935662004 0.9  0.380342361878
    1.0 3.0 2.0   0.1  1.299337679930 0.1  0.120590593900
    1.0 3.0 2.0   0.25 1.471763170430 0.25 0.330518198241
    1.0 3.0 2.0   0.5  1.317308776563 0.5  0.691824033672
    1.0 3.0 2.0   0.75 0.608057837955 0.75 0.941177703958
    1.0 3.0 2.0   0.9  0.136716409678 0.9  0.995180909363
    5.0 2.0 10.0  0.1  0.000036243101 0.1  0.000000660793
    5.0 2.0 10.0  0.25 0.003208139770 0.25 0.000130622011
    5.0 2.0 10.0  0.5  0.169941598714 0.5  0.012184491666
    5.0 2.0 10.0  0.75 2.018421397558 0.75 0.216714718152
    5.0 2.0 10.0  0.9  4.148574183856 0.9  0.698393063281))

(t/deftest beta-noncentral-icdf
  (t/are [alpha beta ncp p vq]
      (let [dist (sut/distribution :beta-noncentral {:alpha alpha :beta beta :ncp ncp})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    2.0 3.0 4.0   0.05 0.208015681566
    2.0 3.0 4.0   0.25 0.415956394387
    2.0 3.0 4.0   0.50 0.568448264639
    2.0 3.0 4.0   0.75 0.705084938124
    2.0 3.0 4.0   0.95 0.853154799732
    0.5 3.0 2.0   0.05 0.005165264200
    0.5 3.0 2.0   0.25 0.098012371829
    0.5 3.0 2.0   0.50 0.269918817515
    0.5 3.0 2.0   0.75 0.473805150890
    0.5 3.0 2.0   0.95 0.727049134606
    3.0 0.5 2.0   0.05 0.577131450158
    3.0 0.5 2.0   0.25 0.832998260427
    3.0 0.5 2.0   0.50 0.939908079371
    3.0 0.5 2.0   0.75 0.986341615162
    3.0 0.5 2.0   0.95 0.999468372706
    1.0 3.0 2.0   0.05 0.043472035627
    1.0 3.0 2.0   0.25 0.194542695754
    1.0 3.0 2.0   0.50 0.364180939423
    1.0 3.0 2.0   0.75 0.545746826024
    1.0 3.0 2.0   0.95 0.765160934244
    5.0 2.0 10.0  0.05 0.610441529101
    5.0 2.0 10.0  0.25 0.765571094091
    5.0 2.0 10.0  0.50 0.849890853766
    5.0 2.0 10.0  0.75 0.912440716430
    5.0 2.0 10.0  0.95 0.967144512924))

(t/deftest beta-noncentral-mv
  (t/are [alpha beta ncp mean-v var-v]
      (let [dist (sut/distribution :beta-noncentral {:alpha alpha :beta beta :ncp ncp})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    2.0 3.0 4.0   0.554504387282 0.038472237096
    0.5 3.0 2.0   0.303394525742 0.053640516453
    3.0 0.5 2.0   0.883899074534 0.020201816126
    1.0 3.0 2.0   0.378170058914 0.049926346088
    5.0 2.0 10.0  0.827452193839 0.012571458948))

;; reference values from R's `circular` package / direct numerical integration
;; of the closed-form pdf `f(x) = e^(kappa*cos(x-mu))/(2*pi*I_0(kappa))`
;; (dvonmises/pvonmises match `integrate()` of that pdf from `mu-pi` to `x`
;; exactly, confirming both give the same [mu-pi,mu+pi] principal-branch
;; convention used here). icdf reference values via R's `uniroot` on that same
;; integral; variance reference values via R's `integrate()` (rel.tol=1e-12)
;; of `(x-mu)^2 * pdf(x)` over `[mu-pi, mu+pi]` - an independent numerical
;; method from fastmath's own gk-quadrature-based variance.
;;
;; Both R's plain `besselI(kappa, 0)` and a naive port of fastmath's own first
;; implementation overflow/lose all precision once kappa gets into the
;; hundreds: `besselI(1000, 0)` is `Inf` in R (needs `expon.scaled=TRUE`), and
;; a fixed-resolution numerical scheme (fixed quadrature step count / fixed
;; integration window) silently loses virtually all the probability mass once
;; kappa grows enough that the density's peak (width ~1/sqrt(kappa)) becomes
;; narrower than the scheme's resolution - both were hit and fixed during
;; development (see `von-mises-log-I0` and the kappa-scaled `steps`/
;; `half-width` in `distr/von-mises`); the kappa=100/1000 rows below and the
;; huge-kappa stress sweep exercise exactly that fix.

(t/deftest von-mises
  (t/testing "defaults are mu=0, kappa=1"
    (let [dist (sut/distribution :von-mises nil)]
      (t/is (m/delta-eq 0.3417104886234632 (sut/pdf dist 0.0)))))
  (t/testing "support is [mu-pi, mu+pi]; pdf/cdf are exactly 0/0/1 outside it"
    (let [dist (sut/distribution :von-mises {:mu 0.5 :kappa 3.0})]
      (t/is (m/delta-eq (- 0.5 m/PI) (sut/lower-bound dist)))
      (t/is (m/delta-eq (+ 0.5 m/PI) (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist (+ 0.5 m/PI 0.1))))
      (t/is (m/delta-eq 0.0 (sut/pdf dist (- 0.5 m/PI 0.1))))
      (t/is (m/delta-eq 0.0 (sut/cdf dist (- 0.5 m/PI 0.1))))
      (t/is (m/delta-eq 1.0 (sut/cdf dist (+ 0.5 m/PI 0.1))))
      (t/is (m/delta-eq 0.0 (sut/cdf dist (- 0.5 m/PI))))
      (t/is (m/delta-eq 1.0 (sut/cdf dist (+ 0.5 m/PI))))))
  (t/testing "kappa=0 reduces exactly to the uniform distribution on the circle [mu-pi, mu+pi]"
    (let [dist (sut/distribution :von-mises {:mu 0.5 :kappa 0.0})]
      (doseq [x [-2.0 -0.5 0.5 1.0 3.0]]
        (t/is (m/delta-eq (/ 1.0 m/TWO_PI) (sut/pdf dist x))))
      (t/is (m/delta-eq 0.5 (sut/mean dist)))
      (t/is (m/delta-eq (/ (m/sq m/TWO_PI) 12.0) (sut/variance dist)))))
  (t/testing "pdf/cdf/icdf/mean/variance are never NaN or throw, across extreme mu/kappa (including a very narrow, huge-kappa peak) and extreme x/p inputs"
    (doseq [mu [-1000.0 -3.0 0.0 3.0 1000.0]
            kappa [0.0 0.001 1.0 100.0 1.0e6 1.0e10 1.0e15]]
      (let [dist (sut/distribution :von-mises {:mu mu :kappa kappa})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0 mu (+ mu m/PI) (- mu m/PI)]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf mu=" mu " kappa=" kappa " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf mu=" mu " kappa=" kappa " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf mu=" mu " kappa=" kappa " p=" p)))
        (t/is (not (Double/isNaN (sut/mean dist))) (str "mean mu=" mu " kappa=" kappa))
        (t/is (not (Double/isNaN (sut/variance dist))) (str "variance mu=" mu " kappa=" kappa))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [mu kappa vx d vp p]
      (let [dist (sut/distribution :von-mises {:mu mu :kappa kappa})]
        (and (m/delta-eq d (sut/pdf dist vx) 1.0e-6)
             (m/delta-eq p (sut/cdf dist vp) 1.0e-6)))
    0.0 2.0    -3   0.009639793410 -3   0.001346862289
    0.0 2.0    -1.5 0.080427734601 -1.5 0.042833151659
    0.0 2.0    0    0.515885412019 0    0.500000000000
    0.0 2.0    0.5  0.403852533352 0.5  0.738192214419
    0.0 2.0    1.5  0.080427734601 1.5  0.957166848341
    0.0 2.0    3    0.009639793410 3    0.998653137711
    1.5 5.0    -1.5 0.000041387928 -1.5 0.000005668660
    1.5 5.0    0    0.008321832185 0    0.001725673014
    1.5 5.0    1.5  0.867136528542 1.5  0.500000000000
    1.5 5.0    2    0.470177012529 2    0.858662481809
    1.5 5.0    3    0.008321832185 3    0.998274326986
    1.5 5.0    4.5  0.000041387928 4.5  0.999994331340
    -1.0 0.5   -4   0.091225297646 -4   0.012873844040
    -1.0 0.5   -2.5 0.155042162215 -2.5 0.183858039433
    -1.0 0.5   -1   0.246738357394 -1   0.500000000000
    -1.0 0.5   -0.5 0.232088734384 -0.5 0.620877026816
    -1.0 0.5   0.5  0.155042162215 0.5  0.816141960567
    -1.0 0.5   2    0.091225297646 2    0.987126155960
    0.5 20.0   -1   0.000000015037 -1   0.000000000753
    0.5 20.0   0.5  1.772715417786 0.5  0.500000000000
    0.5 20.0   1    0.153226776312 1    0.986033746281
    0.5 20.0   2    0.000000015037 2    0.999999999247
    2.0 0.1    -1   0.143793828005 -1   0.020346576035
    2.0 0.1    0.5  0.159884790023 0.5  0.245385905267
    2.0 0.1    2    0.175454504089 2    0.500000000000
    2.0 0.1    2.5  0.173319728348 2.5  0.587367590508
    2.0 0.1    3.5  0.159884790023 3.5  0.754614094733
    2.0 0.1    5    0.143793828005 5    0.979653423965))

(t/deftest von-mises-icdf
  (t/are [mu kappa p vq]
      (let [dist (sut/distribution :von-mises {:mu mu :kappa kappa})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    0.0 2.0    0.05 -1.4179661935
    0.0 2.0    0.25 -0.5296631837
    0.0 2.0    0.50 0.0000000000
    0.0 2.0    0.75 0.5296631837
    0.0 2.0    0.95 1.4179661935
    1.5 5.0    0.05 0.7216414421
    1.5 5.0    0.25 1.1882983440
    1.5 5.0    0.50 1.5000000000
    1.5 5.0    0.75 1.8117016560
    1.5 5.0    0.95 2.2783585579
    -1.0 0.5   0.05 -3.6038014382
    -1.0 0.5   0.25 -2.1124469363
    -1.0 0.5   0.50 -1.0000000000
    -1.0 0.5   0.75 0.1124469363
    -1.0 0.5   0.95 1.6038014382
    0.5 20.0   0.05 0.1276420628
    0.5 20.0   0.25 0.3480598525
    0.5 20.0   0.50 0.5000000000
    0.5 20.0   0.75 0.6519401475
    0.5 20.0   0.95 0.8723579372
    2.0 0.1    0.05 -0.7942200842
    2.0 0.1    0.25 0.5288174455
    2.0 0.1    0.50 2.0000000000
    2.0 0.1    0.75 3.4711825545
    2.0 0.1    0.95 4.7942200842
    0.0 100.0  0.05 -0.1648795459
    0.0 100.0  0.25 -0.0675466545
    0.0 100.0  0.50 0.0000000000
    0.0 100.0  0.75 0.0675466545
    0.0 100.0  0.95 0.1648795459
    0.0 1000.0 0.05 -0.0520272142
    0.0 1000.0 0.25 -0.0213323110
    0.0 1000.0 0.50 0.0000000000
    0.0 1000.0 0.75 0.0213323110
    0.0 1000.0 0.95 0.0520272142))

(t/deftest von-mises-mv
  (t/testing "mean is always exactly mu, by symmetry"
    (doseq [mu [-2.0 0.0 3.0] kappa [0.0 0.5 5.0 100.0]]
      (t/is (m/delta-eq mu (sut/mean (sut/distribution :von-mises {:mu mu :kappa kappa}))))))
  (t/are [mu kappa var-v]
      (m/delta-eq var-v (sut/variance (sut/distribution :von-mises {:mu mu :kappa kappa})) 1.0e-6)
    0.0 2.0   0.764461879811
    1.5 5.0   0.227230162815
    -1.0 0.5  2.348803343669
    0.5 20.0  0.051323846750
    2.0 0.1   3.091356460618
    0.0 100.0 0.010050550607
    0.0 1000.0 0.001000500543)
  (t/testing "variance approaches the wrapped-normal limit 1/kappa as kappa grows"
    (t/is (m/delta-eq (/ 1.0 1.0e6) (sut/variance (sut/distribution :von-mises {:mu 0.0 :kappa 1.0e6})) 1.0e-9))))

;; reference values for the generalized (exponentiated) half-logistic
;; distribution: F(x) = ((1-e^(-lambda*x))/(1+e^(-lambda*x)))^alpha, x>=0.
;; No R package implements this distribution directly, so references were
;; computed independently in R from the closed-form pdf/cdf/icdf formulas
;; (cross-checked live in nREPL), and mean/variance from R's own
;; high-precision `integrate()` over the same pdf (rel.tol=1e-12) - an
;; independent numeric method from fastmath's own gk-quadrature-based
;; mean/variance, so agreement is a genuine cross-check, not circular.

(t/deftest ghl
  (t/testing "registered under both :generalized-half-logistic and :ghl keys"
    (let [d1 (sut/distribution :generalized-half-logistic {:alpha 2.0 :lambda 1.5})
          d2 (sut/distribution :ghl {:alpha 2.0 :lambda 1.5})]
      (t/is (m/delta-eq (sut/pdf d1 1.0) (sut/pdf d2 1.0)))
      (t/is (m/delta-eq (sut/cdf d1 1.0) (sut/cdf d2 1.0)))))
  (t/testing "defaults are alpha=1, lambda=1"
    (let [dist (sut/distribution :ghl nil)]
      (t/is (m/delta-eq 0.470007424403 (sut/pdf dist 0.5)))))
  (t/testing "support is [0, +Inf); pdf/cdf are exactly 0 below 0"
    (let [dist (sut/distribution :ghl {:alpha 2.0 :lambda 1.5})]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))))
  (t/testing "alpha=1 reduces exactly to the ordinary half-logistic distribution"
    (let [dist (sut/distribution :ghl {:alpha 1.0 :lambda 1.3})]
      (doseq [x [0.0 0.3 1.0 2.5 5.0]]
        (t/is (m/delta-eq (sut/pdf dist x) (* 1.3 0.5 (m/sech (* 0.5 1.3 x)) (m/sech (* 0.5 1.3 x)))))
        (t/is (m/delta-eq (sut/cdf dist x) (Math/tanh (* 0.5 1.3 x)))))))
  (t/testing "pdf(0) follows the 0^(alpha-1) convention: 0 (alpha>1), lambda/2 (alpha=1), +Inf (alpha<1)"
    (t/is (m/delta-eq 0.0 (sut/pdf (sut/distribution :ghl {:alpha 2.0 :lambda 1.0}) 0.0)))
    (t/is (m/delta-eq 0.5 (sut/pdf (sut/distribution :ghl {:alpha 1.0 :lambda 1.0}) 0.0)))
    (t/is (Double/isInfinite (sut/pdf (sut/distribution :ghl {:alpha 0.5 :lambda 1.0}) 0.0))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme/infinite/negative-zero inputs"
    (doseq [alpha [0.1 0.5 1.0 2.0 5.0 20.0]
            lambda [0.1 1.0 5.0]]
      (let [dist (sut/distribution :ghl {:alpha alpha :lambda lambda})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf alpha=" alpha " lambda=" lambda " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf alpha=" alpha " lambda=" lambda " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf alpha=" alpha " lambda=" lambda " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [alpha lambda vd d vp p]
      (let [dist (sut/distribution :ghl {:alpha alpha :lambda lambda})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 1.0 0    0.500000000000 0   0.000000000000
    1.0 1.0 0.2  0.495033145424 0.2 0.099667994625
    1.0 1.0 0.5  0.470007424403 0.5 0.244918662404
    1.0 1.0 1    0.393223866483 1   0.462117157260
    1.0 1.0 2    0.209987170807 2   0.761594155956
    1.0 1.0 5    0.013296113342 5   0.986614298151
    1.0 1.0 10   0.000090791615 10  0.999909204263
    2.0 1.0 0    0.000000000000 0   0.000000000000
    2.0 1.0 0.2  0.098677921755 0.2 0.009933709153
    2.0 1.0 0.5  0.230227179409 0.5 0.059985151194
    2.0 1.0 1    0.363430990692 1   0.213552267034
    2.0 1.0 2    0.319850004225 2   0.580025658386
    2.0 1.0 5    0.026236271065 5   0.973407773317
    2.0 1.0 10   0.000181566744 10  0.999818416769
    0.5 2.0 0.2  1.081599287710 0.2 0.444269423014
    0.5 2.0 0.5  0.578447332474 0.5 0.679791995584
    0.5 2.0 1    0.240619578027 1   0.872693620898
    0.5 2.0 2    0.035978455144 2   0.981849061758
    0.5 2.0 5    0.000090795737 5   0.999954601101
    3.0 0.5 0    0.000000000000 0   0.000000000000
    3.0 0.5 0.2  0.001867207511 0.2 0.000124688072
    3.0 0.5 0.5  0.011418407992 0.5 0.001922953665
    3.0 0.5 1    0.042290199622 1   0.014691482994
    3.0 0.5 2    0.125960772209 2   0.098686166568
    3.0 0.5 5    0.151336776754 5   0.610412296576
    3.0 0.5 10   0.019413810122 10  0.960378027086
    1.5 1.5 0.2  0.424465936273 0.2 0.057448218432
    1.5 1.5 0.5  0.586972768600 0.5 0.214523346138
    1.5 1.5 1    0.534888466504 1   0.506189787776
    1.5 1.5 2    0.193413368915 2   0.861151528479
    1.5 1.5 5    0.002484754168 5   0.998342122520
    1.5 1.5 10   0.000001376559 10  0.999999082293))

(t/deftest ghl-icdf
  (t/are [alpha lambda p vq]
      (let [dist (sut/distribution :ghl {:alpha alpha :lambda lambda})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 1.0 0.05 0.100083458557
    1.0 1.0 0.25 0.510825623766
    1.0 1.0 0.50 1.098612288668
    1.0 1.0 0.75 1.945910149055
    1.0 1.0 0.95 3.663561646130
    2.0 1.0 0.05 0.454899072012
    2.0 1.0 0.25 1.098612288668
    2.0 1.0 0.50 1.762747174039
    2.0 1.0 0.75 2.633915793850
    2.0 1.0 0.95 4.356544420602
    0.5 2.0 0.05 0.002500005208
    0.5 2.0 0.25 0.062581571477
    0.5 2.0 0.50 0.255412811883
    0.5 2.0 0.75 0.636482837906
    0.5 2.0 0.95 1.485535855866
    3.0 0.5 0.05 1.546296919419
    3.0 0.5 0.25 2.965402772811
    3.0 0.5 0.50 4.325414454657
    3.0 0.5 0.75 6.076849374933
    3.0 0.5 0.95 9.523958157970
    1.5 1.5 0.05 0.182084729578
    1.5 1.5 0.25 0.559873121896
    1.5 1.5 0.50 0.988467590937
    1.5 1.5 0.75 1.565046849143
    1.5 1.5 0.95 2.712603317056))

(t/deftest ghl-mv
  (t/are [alpha lambda mean-v var-v]
      (let [dist (sut/distribution :ghl {:alpha alpha :lambda lambda})]
        (and (m/delta-eq mean-v (sut/mean dist) 1.0e-6)
             (m/delta-eq var-v (sut/variance dist) 1.0e-6)))
    1.0 1.0 1.386294361120 1.368056078024
    2.0 1.0 2.000000000000 1.545177444480
    0.5 2.0 0.438824586260 0.266326000185
    3.0 0.5 4.772588722240 6.381869423136
    1.5 1.5 1.157370995080 0.661200406444))

;; reference values from R's gamlss.dist package: dZINBI/pZINBI/qZINBI

(t/deftest zinbi
  (t/testing "registered under both :zero-inflated-negative-binomial and :zinbi keys"
    (let [d1 (sut/distribution :zero-inflated-negative-binomial {:mu 3.0 :sigma 0.5 :nu 0.1})
          d2 (sut/distribution :zinbi {:mu 3.0 :sigma 0.5 :nu 0.1})]
      (t/is (m/delta-eq (sut/pdf d1 2) (sut/pdf d2 2)))
      (t/is (m/delta-eq (sut/cdf d1 2) (sut/cdf d2 2)))))
  (t/testing "defaults match R's ZINBI() defaults (mu=1, sigma=1, nu=0.3)"
    (let [dist (sut/distribution :zinbi)]
      (t/is (m/delta-eq 0.650000000000 (sut/pdf dist 0)))))
  (t/are [mu sigma nu vd d vp p]
      (let [dist (sut/distribution :zinbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 1.0     0.3  0 0.650000000000 0  0.650000000000
    1.0 1.0     0.3  1 0.175000000000 1  0.825000000000
    1.0 1.0     0.3  3 0.043750000000 3  0.956250000000
    1.0 1.0     0.3  8 0.001367187500 8  0.998632812500
    2.0 0.000009 0.5 0 0.567667641618 0  0.567667641618
    2.0 0.000009 0.5 1 0.135335283237 1  0.703002924855
    2.0 0.000009 0.5 3 0.090223522158 3  0.928561730249
    2.0 0.000009 0.5 8 0.000429635820 8  0.999881276336
    3.0 0.5     0.1  0 0.244000000000 0  0.244000000000
    3.0 0.5     0.1  2 0.155520000000 2  0.572320000000
    3.0 0.5     0.1  5 0.067184640000 5  0.857232640000
    3.0 0.5     0.1  10 0.009577842278 10 0.982368063078
    5.0 2.0     0.4  0 0.580906806747 0  0.580906806747
    5.0 2.0     0.4  2 0.056066159116 2  0.719203332565
    5.0 2.0     0.4  7 0.019446141948 7  0.865589900197
    5.0 2.0     0.4  15 0.006256415931 15 0.950134533725
    0.5 0.8     0.05 0 0.673826128188 0  0.673826128188
    0.5 0.8     0.05 1 0.222795045782 1  0.896621173970
    0.5 0.8     0.05 3 0.022165833636 3  0.990399700893
    0.5 0.8     0.05 6 0.000600795873 6  0.999748084684))

(t/deftest zinbi-icdf
  (t/are [mu sigma nu p vq]
      (let [dist (sut/distribution :zinbi {:mu mu :sigma sigma :nu nu})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 1.0     0.3  0.01 0
    1.0 1.0     0.3  0.15 0
    1.0 1.0     0.3  0.30 0
    1.0 1.0     0.3  0.50 0
    1.0 1.0     0.3  0.70 1
    1.0 1.0     0.3  0.85 2
    1.0 1.0     0.3  0.95 3
    1.0 1.0     0.3  0.99 6
    2.0 0.000009 0.5 0.01 0
    2.0 0.000009 0.5 0.15 0
    2.0 0.000009 0.5 0.30 0
    2.0 0.000009 0.5 0.50 0
    2.0 0.000009 0.5 0.70 1
    2.0 0.000009 0.5 0.85 3
    2.0 0.000009 0.5 0.95 4
    2.0 0.000009 0.5 0.99 5
    3.0 0.5     0.1  0.01 0
    3.0 0.5     0.1  0.15 0
    3.0 0.5     0.1  0.30 1
    3.0 0.5     0.1  0.50 2
    3.0 0.5     0.1  0.70 4
    3.0 0.5     0.1  0.85 5
    3.0 0.5     0.1  0.95 8
    3.0 0.5     0.1  0.99 12
    5.0 2.0     0.4  0.01 0
    5.0 2.0     0.4  0.15 0
    5.0 2.0     0.4  0.30 0
    5.0 2.0     0.4  0.50 0
    5.0 2.0     0.4  0.70 2
    5.0 2.0     0.4  0.85 7
    5.0 2.0     0.4  0.95 15
    5.0 2.0     0.4  0.99 30
    0.5 0.8     0.05 0.01 0
    0.5 0.8     0.05 0.15 0
    0.5 0.8     0.05 0.30 0
    0.5 0.8     0.05 0.50 0
    0.5 0.8     0.05 0.70 1
    0.5 0.8     0.05 0.85 1
    0.5 0.8     0.05 0.95 2
    0.5 0.8     0.05 0.99 3))

(t/deftest zinbi-mv
  (t/are [mu sigma nu mean-v var-v]
      (let [dist (sut/distribution :zinbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist) 1.0e-4)))
    1.0 1.0     0.3  0.700000000000 1.610000000000
    2.0 0.000009 0.5 1.000000000000 2.000000000000
    3.0 0.5     0.1  2.700000000000 7.560000000000
    5.0 2.0     0.4  3.000000000000 39.000000000000
    0.5 0.8     0.05 0.475000000000 0.676875000000))

;; reference values from R's gamlss.dist package: dZANBI/pZANBI/qZANBI

(t/deftest zanbi
  (t/testing "registered under both :zero-adjusted-negative-binomial and :zanbi keys"
    (let [d1 (sut/distribution :zero-adjusted-negative-binomial {:mu 3.0 :sigma 0.5 :nu 0.1})
          d2 (sut/distribution :zanbi {:mu 3.0 :sigma 0.5 :nu 0.1})]
      (t/is (m/delta-eq (sut/pdf d1 2) (sut/pdf d2 2)))
      (t/is (m/delta-eq (sut/cdf d1 2) (sut/cdf d2 2)))))
  (t/testing "defaults match R's ZANBI() defaults (mu=1, sigma=1, nu=0.3)"
    (let [dist (sut/distribution :zanbi)]
      (t/is (m/delta-eq 0.300000000000 (sut/pdf dist 0)))))
  (t/are [mu sigma nu vd d vp p]
      (let [dist (sut/distribution :zanbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 1.0     0.3  0 0.300000000000 0  0.300000000000
    1.0 1.0     0.3  1 0.350000000000 1  0.650000000000
    1.0 1.0     0.3  3 0.087500000000 3  0.912500000000
    1.0 1.0     0.3  8 0.002734375000 8  0.997265625000
    2.0 0.000009 0.5 0 0.500000000000 0  0.500000000000
    2.0 0.000009 0.5 1 0.156517642750 1  0.656517642750
    2.0 0.000009 0.5 3 0.104345095166 3  0.917380380666
    2.0 0.000009 0.5 8 0.000496881406 8  0.999862693988
    3.0 0.5     0.1  0 0.100000000000 0  0.100000000000
    3.0 0.5     0.1  2 0.185142857143 2  0.490857142857
    3.0 0.5     0.1  5 0.079981714286 5  0.830038857143
    3.0 0.5     0.1  10 0.011402193189 10 0.979009598903
    5.0 2.0     0.4  0 0.400000000000 0  0.400000000000
    5.0 2.0     0.4  2 0.080267816350 2  0.597993946996
    5.0 2.0     0.4  7 0.027840311789 7  0.807570103309
    5.0 2.0     0.4  15 0.008957075941 15 0.928609482934
    0.5 0.8     0.05 0 0.050000000000 0  0.050000000000
    0.5 0.8     0.05 1 0.648903274554 1  0.698903274554
    0.5 0.8     0.05 3 0.064559254356 3  0.972038581445
    0.5 0.8     0.05 6 0.001749852239 6  0.999266282278))

(t/deftest zanbi-icdf
  (t/are [mu sigma nu p vq]
      (let [dist (sut/distribution :zanbi {:mu mu :sigma sigma :nu nu})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 1.0     0.3  0.01 0
    1.0 1.0     0.3  0.15 0
    1.0 1.0     0.3  0.30 0
    1.0 1.0     0.3  0.50 1
    1.0 1.0     0.3  0.70 2
    1.0 1.0     0.3  0.85 3
    1.0 1.0     0.3  0.95 4
    1.0 1.0     0.3  0.99 7
    2.0 0.000009 0.5 0.01 0
    2.0 0.000009 0.5 0.15 0
    2.0 0.000009 0.5 0.30 0
    2.0 0.000009 0.5 0.50 0
    2.0 0.000009 0.5 0.70 2
    2.0 0.000009 0.5 0.85 3
    2.0 0.000009 0.5 0.95 4
    2.0 0.000009 0.5 0.99 5
    3.0 0.5     0.1  0.01 0
    3.0 0.5     0.1  0.15 1
    3.0 0.5     0.1  0.30 1
    3.0 0.5     0.1  0.50 3
    3.0 0.5     0.1  0.70 4
    3.0 0.5     0.1  0.85 6
    3.0 0.5     0.1  0.95 8
    3.0 0.5     0.1  0.99 12
    5.0 2.0     0.4  0.01 0
    5.0 2.0     0.4  0.15 0
    5.0 2.0     0.4  0.30 0
    5.0 2.0     0.4  0.50 1
    5.0 2.0     0.4  0.70 4
    5.0 2.0     0.4  0.85 9
    5.0 2.0     0.4  0.95 19
    5.0 2.0     0.4  0.99 33
    0.5 0.8     0.05 0.01 0
    0.5 0.8     0.05 0.15 1
    0.5 0.8     0.05 0.30 1
    0.5 0.8     0.05 0.50 1
    0.5 0.8     0.05 0.70 2
    0.5 0.8     0.05 0.85 2
    0.5 0.8     0.05 0.95 3
    0.5 0.8     0.05 0.99 4))

(t/deftest zanbi-mv
  (t/are [mu sigma nu mean-v var-v]
      (let [dist (sut/distribution :zanbi {:mu mu :sigma sigma :nu nu})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist) 1.0e-4)))
    1.0 1.0     0.3  1.400000000000 2.240000000000
    2.0 0.000009 0.5 1.156517642750 2.132019870258
    3.0 0.5     0.1  3.214285714286 7.346938775510
    5.0 2.0     0.4  4.294987437107 50.272881908802
    0.5 0.8     0.05 1.383464584375 0.714608454092))

;; reference values from R's gamlss.dist package: dZIP/pZIP/qZIP

(t/deftest zip
  (t/testing "registered under both :zero-inflated-poisson and :zip keys"
    (let [d1 (sut/distribution :zero-inflated-poisson {:mu 3.0 :sigma 0.3})
          d2 (sut/distribution :zip {:mu 3.0 :sigma 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 2) (sut/pdf d2 2)))
      (t/is (m/delta-eq (sut/cdf d1 2) (sut/cdf d2 2)))))
  (t/testing "defaults match R's ZIP() defaults (mu=5, sigma=0.1)"
    (let [dist (sut/distribution :zip)]
      (t/is (m/delta-eq 0.106064152299 (sut/pdf dist 0)))))
  (t/are [mu sigma vd d vp p]
      (let [dist (sut/distribution :zip {:mu mu :sigma sigma})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    5.0 0.1  0  0.106064152299 0  0.106064152299
    5.0 0.1  3  0.126336506233 3  0.338523323768
    5.0 0.1  8  0.058750235413 8  0.938715728750
    5.0 0.1  15 0.000141520891 15 0.999937892582
    1.0 0.2  0  0.494303552937 0  0.494303552937
    1.0 0.2  1  0.294303552937 1  0.788607105874
    1.0 0.2  2  0.147151776469 2  0.935758882343
    1.0 0.2  5  0.002452529608 5  0.999524652146
    10.0 0.05 0  0.050043129933 0  0.050043129933
    10.0 0.05 5  0.035941611062 5  0.113731664735
    10.0 0.05 10 0.118854533935 10 0.603887762683
    10.0 0.05 20 0.001772777248 20 0.998491152371
    0.5 0.4  0  0.763918395828 0  0.763918395828
    0.5 0.4  1  0.181959197914 1  0.945877593741
    0.5 0.4  3  0.007581633246 3  0.998949026466
    0.5 0.4  6  0.000007897535 6  0.999999398572
    3.0 0.3  0  0.334850947858 0  0.334850947858
    3.0 0.3  2  0.156829265359 2  0.596233056789
    3.0 0.3  4  0.117621949019 4  0.870684271167
    3.0 0.3  9  0.001890352752 9  0.999228258309))

(t/deftest zip-icdf
  (t/are [mu sigma p vq]
      (let [dist (sut/distribution :zip {:mu mu :sigma sigma})]
        (m/delta-eq vq (sut/icdf dist p)))
    5.0 0.1  0.01 0
    5.0 0.1  0.15 2
    5.0 0.1  0.30 3
    5.0 0.1  0.50 5
    5.0 0.1  0.70 6
    5.0 0.1  0.85 7
    5.0 0.1  0.95 9
    5.0 0.1  0.99 11
    1.0 0.2  0.01 0
    1.0 0.2  0.15 0
    1.0 0.2  0.30 0
    1.0 0.2  0.50 1
    1.0 0.2  0.70 1
    1.0 0.2  0.85 2
    1.0 0.2  0.95 3
    1.0 0.2  0.99 4
    10.0 0.05 0.01 0
    10.0 0.05 0.15 6
    10.0 0.05 0.30 8
    10.0 0.05 0.50 10
    10.0 0.05 0.70 11
    10.0 0.05 0.85 13
    10.0 0.05 0.95 15
    10.0 0.05 0.99 18
    0.5 0.4  0.01 0
    0.5 0.4  0.15 0
    0.5 0.4  0.30 0
    0.5 0.4  0.50 0
    0.5 0.4  0.70 0
    0.5 0.4  0.85 1
    0.5 0.4  0.95 2
    0.5 0.4  0.99 2
    3.0 0.3  0.01 0
    3.0 0.3  0.15 0
    3.0 0.3  0.30 0
    3.0 0.3  0.50 2
    3.0 0.3  0.70 3
    3.0 0.3  0.85 4
    3.0 0.3  0.95 6
    3.0 0.3  0.99 7))

(t/deftest zip-mv
  (t/are [mu sigma mean-v var-v]
      (let [dist (sut/distribution :zip {:mu mu :sigma sigma})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    5.0 0.1  4.500000000000 6.750000000000
    1.0 0.2  0.800000000000 0.960000000000
    10.0 0.05 9.500000000000 14.250000000000
    0.5 0.4  0.300000000000 0.360000000000
    3.0 0.3  2.100000000000 3.990000000000))

;; reference values from R's gamlss.dist package: dZIP2/pZIP2/qZIP2

(t/deftest zip2
  (t/testing "registered under both :zero-inflated-poisson2 and :zip2 keys"
    (let [d1 (sut/distribution :zero-inflated-poisson2 {:mu 3.0 :sigma 0.3})
          d2 (sut/distribution :zip2 {:mu 3.0 :sigma 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 2) (sut/pdf d2 2)))
      (t/is (m/delta-eq (sut/cdf d1 2) (sut/cdf d2 2)))))
  (t/testing "defaults match R's ZIP2() defaults (mu=5, sigma=0.1), mu is the mean directly"
    (let [dist (sut/distribution :zip2)]
      (t/is (m/delta-eq 0.103479328126 (sut/pdf dist 0)))
      (t/is (m/delta-eq 5.0 (sut/mean dist)))))
  (t/are [mu sigma vd d vp p]
      (let [dist (sut/distribution :zip2 {:mu mu :sigma sigma})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    5.0 0.1  0  0.103479328126 0  0.103479328126
    5.0 0.1  3  0.099432102353 3  0.275934366446
    5.0 0.1  8  0.078305960952 8  0.900623991274
    5.0 0.1  15 0.000394373932 15 0.999798947175
    1.0 0.2  0  0.429203837488 0  0.429203837488
    1.0 0.2  1  0.286504796860 1  0.715708634348
    1.0 0.2  2  0.179065498038 2  0.894774132386
    1.0 0.2  5  0.005828955014 5  0.998529531640
    10.0 0.05 0  0.050025480195 0  0.050025480195
    10.0 0.05 5  0.027441223478 5  0.097134283377
    10.0 0.05 10 0.117274518194 10 0.541613051414
    10.0 0.05 20 0.002921500873 20 0.997275586344
    0.5 0.4  0  0.660758925104 0  0.660758925104
    0.5 0.4  1  0.217299104254 1  0.878058029358
    0.5 0.4  3  0.025150359289 3  0.993749682085
    0.5 0.4  6  0.000121288384 6  0.999983905019
    3.0 0.3  0  0.309634650713 0  0.309634650713
    3.0 0.3  2  0.088481486141 2  0.439407497053
    3.0 0.3  4  0.135430846134 4  0.701240466246
    3.0 0.3  9  0.012950360981 9  0.991152427303))

(t/deftest zip2-icdf
  (t/are [mu sigma p vq]
      (let [dist (sut/distribution :zip2 {:mu mu :sigma sigma})]
        (m/delta-eq vq (sut/icdf dist p)))
    5.0 0.1  0.01 0
    5.0 0.1  0.15 2
    5.0 0.1  0.30 4
    5.0 0.1  0.50 5
    5.0 0.1  0.70 6
    5.0 0.1  0.85 8
    5.0 0.1  0.95 10
    5.0 0.1  0.99 12
    1.0 0.2  0.01 0
    1.0 0.2  0.15 0
    1.0 0.2  0.30 0
    1.0 0.2  0.50 1
    1.0 0.2  0.70 1
    1.0 0.2  0.85 2
    1.0 0.2  0.95 3
    1.0 0.2  0.99 4
    10.0 0.05 0.01 0
    10.0 0.05 0.15 7
    10.0 0.05 0.30 8
    10.0 0.05 0.50 10
    10.0 0.05 0.70 12
    10.0 0.05 0.85 14
    10.0 0.05 0.95 16
    10.0 0.05 0.99 19
    0.5 0.4  0.01 0
    0.5 0.4  0.15 0
    0.5 0.4  0.30 0
    0.5 0.4  0.50 0
    0.5 0.4  0.70 1
    0.5 0.4  0.85 1
    0.5 0.4  0.95 2
    0.5 0.4  0.99 3
    3.0 0.3  0.01 0
    3.0 0.3  0.15 0
    3.0 0.3  0.30 0
    3.0 0.3  0.50 3
    3.0 0.3  0.70 4
    3.0 0.3  0.85 6
    3.0 0.3  0.95 7
    3.0 0.3  0.99 9))

(t/deftest zip2-mv
  (t/are [mu sigma mean-v var-v]
      (let [dist (sut/distribution :zip2 {:mu mu :sigma sigma})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    5.0 0.1  5.000000000000 7.777777777778
    1.0 0.2  1.000000000000 1.250000000000
    10.0 0.05 10.000000000000 15.263157894737
    0.5 0.4  0.500000000000 0.666666666667
    3.0 0.3  3.000000000000 6.857142857143))

;; reference values from R's gamlss.dist package: dZAP/pZAP/qZAP

(t/deftest zap
  (t/testing "registered under both :zero-adjusted-poisson and :zap keys"
    (let [d1 (sut/distribution :zero-adjusted-poisson {:mu 3.0 :sigma 0.3})
          d2 (sut/distribution :zap {:mu 3.0 :sigma 0.3})]
      (t/is (m/delta-eq (sut/pdf d1 2) (sut/pdf d2 2)))
      (t/is (m/delta-eq (sut/cdf d1 2) (sut/cdf d2 2)))))
  (t/testing "defaults match R's ZAP() defaults (mu=5, sigma=0.1)"
    (let [dist (sut/distribution :zap)]
      (t/is (m/delta-eq 0.1 (sut/pdf dist 0)))))
  (t/are [mu sigma vd d vp p]
      (let [dist (sut/distribution :zap {:mu mu :sigma sigma})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    5.0 0.1   0  0.100000000000 0  0.100000000000
    5.0 0.1   3  0.127193529493 3  0.334036094267
    5.0 0.1   8  0.059148776736 8  0.938299997403
    5.0 0.1   15 0.000142480920 15 0.999937471267
    1.0 0.2   0  0.200000000000 0  0.200000000000
    1.0 0.2   1  0.465581365495 1  0.665581365495
    1.0 0.2   2  0.232790682748 2  0.898372048243
    1.0 0.2   5  0.003879844712 5  0.999248010767
    10.0 0.05 0  0.050000000000 0  0.050000000000
    10.0 0.05 5  0.035943242883 5  0.113691426388
    10.0 0.05 10 0.118859930168 10 0.603869778399
    10.0 0.05 20 0.001772857736 20 0.998491083867
    0.5 0.4   0  0.400000000000 0  0.400000000000
    0.5 0.4   1  0.462448224761 1  0.862448224761
    0.5 0.4   3  0.019268676032 3  0.997328956983
    0.5 0.4   6  0.000020071538 6  0.999998471475
    3.0 0.3   0  0.300000000000 0  0.300000000000
    3.0 0.3   2  0.165046443947 2  0.575077406579
    3.0 0.3   4  0.123784832961 4  0.863908683487
    3.0 0.3   9  0.001989399101 9  0.999187822366))

(t/deftest zap-icdf
  (t/are [mu sigma p vq]
      (let [dist (sut/distribution :zap {:mu mu :sigma sigma})]
        (m/delta-eq vq (sut/icdf dist p)))
    5.0 0.1   0.01 0
    5.0 0.1   0.15 2
    5.0 0.1   0.30 3
    5.0 0.1   0.50 5
    5.0 0.1   0.70 6
    5.0 0.1   0.85 7
    5.0 0.1   0.95 9
    5.0 0.1   0.99 11
    1.0 0.2   0.01 0
    1.0 0.2   0.15 0
    1.0 0.2   0.30 1
    1.0 0.2   0.50 1
    1.0 0.2   0.70 2
    1.0 0.2   0.85 2
    1.0 0.2   0.95 3
    1.0 0.2   0.99 4
    10.0 0.05 0.01 0
    10.0 0.05 0.15 6
    10.0 0.05 0.30 8
    10.0 0.05 0.50 10
    10.0 0.05 0.70 11
    10.0 0.05 0.85 13
    10.0 0.05 0.95 15
    10.0 0.05 0.99 18
    0.5 0.4   0.01 0
    0.5 0.4   0.15 0
    0.5 0.4   0.30 0
    0.5 0.4   0.50 1
    0.5 0.4   0.70 1
    0.5 0.4   0.85 1
    0.5 0.4   0.95 2
    0.5 0.4   0.99 3
    3.0 0.3   0.01 0
    3.0 0.3   0.15 0
    3.0 0.3   0.30 0
    3.0 0.3   0.50 2
    3.0 0.3   0.70 3
    3.0 0.3   0.85 4
    3.0 0.3   0.95 6
    3.0 0.3   0.99 7))

(t/deftest zap-mv
  (t/are [mu sigma mean-v var-v]
      (let [dist (sut/distribution :zap {:mu mu :sigma sigma})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    5.0 0.1   4.530526447078 6.657488794794
    1.0 0.2   1.265581365495 0.929466538302
    10.0 0.05 9.500431318915 14.246549262647
    0.5 0.4   0.762448224761 0.562345041700
    3.0 0.3   2.210030962632 3.955886994736))

;; reference values from R's gamlss.dist package: dZIBB/pZIBB/qZIBB

(t/deftest zero-inflated-beta-binomial
  (t/testing "registered under both :zero-inflated-beta-binomial and :zibb keys"
    (let [d1 (sut/distribution :zero-inflated-beta-binomial {:mu 0.3 :sigma 0.4 :nu 0.2 :bd 10})
          d2 (sut/distribution :zibb {:mu 0.3 :sigma 0.4 :nu 0.2 :bd 10})]
      (t/is (m/delta-eq (sut/pdf d1 3) (sut/pdf d2 3)))
      (t/is (m/delta-eq (sut/cdf d1 3) (sut/cdf d2 3)))))
  (t/are [mu sigma nu bd vd d vp p]
      (let [dist (sut/distribution :zibb {:mu mu :sigma sigma :nu nu :bd bd})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.5 0.5 0.1  1  0 0.550000000000 0  0.550000000000
    0.5 0.5 0.1  1  1 0.450000000000 1  1.000000000000
    0.3 0.4 0.2  10 0 0.383777573529 0  0.383777573529
    0.3 0.4 0.2  10 3 0.086792986425 3  0.702347285068
    0.3 0.4 0.2  10 7 0.041996606335 7  0.931617647059
    0.3 0.4 0.2  10 10 0.012821691176 10 1.000000000000
    0.7 0.2 0.05 20 0 0.050459550939 0  0.050459550939
    0.7 0.2 0.05 20 5 0.014090068093 5  0.085382416301
    0.7 0.2 0.05 20 13 0.070450340463 13 0.436240030225
    0.7 0.2 0.05 20 20 0.059282071191 20 1.000000000000
    0.1 1.5 0.3  15 0 0.831050984000 0  0.831050984000
    0.1 1.5 0.3  15 1 0.036373355069 1  0.867424339069
    0.1 1.5 0.3  15 7 0.007580091280 7  0.938283482226
    0.1 1.5 0.3  15 15 0.012956199430 15 1.000000000000
    0.6 0.9 1.0e-4 5 0 0.156750970552 0 0.156750970552
    0.6 0.9 1.0e-4 5 2 0.113698285078 2 0.387937483544
    0.6 0.9 1.0e-4 5 4 0.157428394724 4 0.669400371080
    0.6 0.9 1.0e-4 5 5 0.330599628920 5 1.000000000000))

(t/deftest zero-inflated-beta-binomial-icdf
  (t/are [mu sigma nu bd p vq]
      (let [dist (sut/distribution :zibb {:mu mu :sigma sigma :nu nu :bd bd})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.5 0.5 0.1  1  0.01 0
    0.5 0.5 0.1  1  0.50 0
    0.5 0.5 0.1  1  0.70 1
    0.5 0.5 0.1  1  0.99 1
    0.3 0.4 0.2  10 0.01 0
    0.3 0.4 0.2  10 0.30 0
    0.3 0.4 0.2  10 0.50 1
    0.3 0.4 0.2  10 0.70 3
    0.3 0.4 0.2  10 0.90 7
    0.3 0.4 0.2  10 0.99 10
    0.7 0.2 0.05 20 0.01 0
    0.7 0.2 0.05 20 0.10 6
    0.7 0.2 0.05 20 0.50 14
    0.7 0.2 0.05 20 0.90 19
    0.7 0.2 0.05 20 0.99 20
    0.1 1.5 0.3  15 0.01 0
    0.1 1.5 0.3  15 0.70 0
    0.1 1.5 0.3  15 0.90 3
    0.1 1.5 0.3  15 0.99 15
    0.6 0.9 1.0e-4 5 0.10 0
    0.6 0.9 1.0e-4 5 0.30 2
    0.6 0.9 1.0e-4 5 0.70 5
    0.6 0.9 1.0e-4 5 0.99 5))

(t/deftest zero-inflated-beta-binomial-mv
  (t/are [mu sigma nu bd mean-v var-v]
      (let [dist (sut/distribution :zibb {:mu mu :sigma sigma :nu nu :bd bd})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.5 0.5 0.1    1  0.450000000000 0.247500000000
    0.3 0.4 0.2    10 2.400000000000 7.440000000000
    0.7 0.2 0.05   20 13.300000000000 25.935000000000
    0.1 1.5 0.3    15 1.050000000000 9.355500000000
    0.6 0.9 1.0e-4 5  2.999700000000 3.474236752105))

;; reference values from R's gamlss.dist package: dZABB/pZABB/qZABB
;; note: R's ZABB default sigma is 0.1, unlike ZIBB/BB whose default sigma is 0.5

(t/deftest zero-adjusted-beta-binomial
  (t/testing "registered under both :zero-adjusted-beta-binomial and :zabb keys"
    (let [d1 (sut/distribution :zero-adjusted-beta-binomial {:mu 0.3 :sigma 0.4 :nu 0.2 :bd 10})
          d2 (sut/distribution :zabb {:mu 0.3 :sigma 0.4 :nu 0.2 :bd 10})]
      (t/is (m/delta-eq (sut/pdf d1 3) (sut/pdf d2 3)))
      (t/is (m/delta-eq (sut/cdf d1 3) (sut/cdf d2 3)))))
  (t/testing "defaults match R's ZABB() defaults (mu=0.5, sigma=0.1, nu=0.1, bd=1)"
    (let [dist (sut/distribution :zabb)]
      (t/is (m/delta-eq 0.1 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.9 (sut/pdf dist 1)))))
  (t/are [mu sigma nu bd vd d vp p]
      (let [dist (sut/distribution :zabb {:mu mu :sigma sigma :nu nu :bd bd})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.5 0.1 0.1  1  0 0.100000000000 0  0.100000000000
    0.5 0.1 0.1  1  1 0.900000000000 1  1.000000000000
    0.3 0.4 0.2  10 0 0.200000000000 0  0.200000000000
    0.3 0.4 0.2  10 3 0.112677478387 3  0.613577562716
    0.3 0.4 0.2  10 7 0.054521360510 7  0.911223804907
    0.3 0.4 0.2  10 10 0.016645536580 10 1.000000000000
    0.7 0.2 0.05 20 0 0.050000000000 0  0.050000000000
    0.7 0.2 0.05 20 5 0.014096887290 5  0.084939767049
    0.7 0.2 0.05 20 13 0.070484436452 13 0.435967186215
    0.7 0.2 0.05 20 20 0.059310762051 20 1.000000000000
    0.1 1.5 0.3  15 0 0.300000000000 0  0.300000000000
    0.1 1.5 0.3  15 1 0.150704331702 1  0.450704331702
    0.1 1.5 0.3  15 7 0.031406302456 7  0.744292310991
    0.1 1.5 0.3  15 15 0.053680925853 15 1.000000000000
    0.6 0.9 1.0e-4 5 0 0.000100000000 0 0.000100000000
    0.6 0.9 1.0e-4 5 2 0.134820096175 2 0.274234195556
    0.6 0.9 1.0e-4 5 4 0.186673979319 4 0.607984643430
    0.6 0.9 1.0e-4 5 5 0.392015356570 5 1.000000000000))

(t/deftest zero-adjusted-beta-binomial-icdf
  (t/are [mu sigma nu bd p vq]
      (let [dist (sut/distribution :zabb {:mu mu :sigma sigma :nu nu :bd bd})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.5 0.1 0.1  1  0.01 0
    0.5 0.1 0.1  1  0.30 1
    0.5 0.1 0.1  1  0.99 1
    0.3 0.4 0.2  10 0.01 0
    0.3 0.4 0.2  10 0.30 1
    0.3 0.4 0.2  10 0.50 2
    0.3 0.4 0.2  10 0.70 4
    0.3 0.4 0.2  10 0.90 7
    0.3 0.4 0.2  10 0.99 10
    0.7 0.2 0.05 20 0.01 0
    0.7 0.2 0.05 20 0.10 6
    0.7 0.2 0.05 20 0.50 14
    0.7 0.2 0.05 20 0.90 19
    0.7 0.2 0.05 20 0.99 20
    0.1 1.5 0.3  15 0.01 0
    0.1 1.5 0.3  15 0.30 0
    0.1 1.5 0.3  15 0.50 2
    0.1 1.5 0.3  15 0.90 13
    0.1 1.5 0.3  15 0.99 15
    0.6 0.9 1.0e-4 5 0.01 1
    0.6 0.9 1.0e-4 5 0.30 3
    0.6 0.9 1.0e-4 5 0.70 5
    0.6 0.9 1.0e-4 5 0.99 5))

(t/deftest zero-adjusted-beta-binomial-mv
  (t/are [mu sigma nu bd mean-v var-v]
      (let [dist (sut/distribution :zabb {:mu mu :sigma sigma :nu nu :bd bd})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.5 0.1 0.1    1  0.900000000000 0.090000000000
    0.3 0.4 0.2    10 3.115758072936 7.428721032082
    0.7 0.2 0.05   20 13.306436826889 25.861900582067
    0.1 1.5 0.3    15 4.350424864279 24.404035148482
    0.6 0.9 1.0e-4 5  3.556956397523 2.137511470557))

;; reference values from R's gamlss.dist package: dZIBI/pZIBI/qZIBI

(t/deftest zero-inflated-binomial
  (t/testing "registered under both :zero-inflated-binomial and :zibi keys"
    (let [d1 (sut/distribution :zero-inflated-binomial {:mu 0.3 :sigma 0.4 :bd 10})
          d2 (sut/distribution :zibi {:mu 0.3 :sigma 0.4 :bd 10})]
      (t/is (m/delta-eq (sut/pdf d1 3) (sut/pdf d2 3)))
      (t/is (m/delta-eq (sut/cdf d1 3) (sut/cdf d2 3)))))
  (t/testing "defaults match R's ZIBI() defaults (mu=0.5, sigma=0.1, bd=1)"
    (let [dist (sut/distribution :zibi)]
      (t/is (m/delta-eq 0.550000000000 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.450000000000 (sut/pdf dist 1)))))
  (t/are [mu sigma bd vd d vp p]
      (let [dist (sut/distribution :zibi {:mu mu :sigma sigma :bd bd})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.5 0.1     1  0 0.550000000000 0  0.550000000000
    0.5 0.1     1  1 0.450000000000 1  1.000000000000
    0.3 0.4     10 0 0.416948514940 0  0.416948514940
    0.3 0.4     10 3 0.160096759200 3  0.789766431040
    0.3 0.4     10 7 0.005401015200 7  0.999045768160
    0.3 0.4     10 10 0.000003542940 10 1.000000000000
    0.7 0.05    20 0 0.050000000033 0  0.050000000033
    0.7 0.05    20 5 0.000035520280 5  0.050040793021
    0.7 0.05    20 13 0.156048885956 13 0.422390678409
    0.7 0.05    20 20 0.000758026530 20 1.000000000000
    0.1 0.3     15 0 0.444123792466 0  0.444123792466
    0.1 0.3     15 1 0.240206320777 1  0.684330113243
    0.1 0.3     15 7 0.000193903955 7  0.999976462578
    0.1 0.3     15 15 0.000000000000 15 1.000000000000
    0.6 1.0e-4  5  0 0.010338976000 0  0.010338976000
    0.6 1.0e-4  5  2 0.230376960000 2  0.317508256000
    0.6 1.0e-4  5  4 0.259174080000 4  0.922247776000
    0.6 1.0e-4  5  5 0.077752224000 5  1.000000000000))

(t/deftest zero-inflated-binomial-icdf
  (t/are [mu sigma bd p vq]
      (let [dist (sut/distribution :zibi {:mu mu :sigma sigma :bd bd})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.5 0.1    1  0.01 0
    0.5 0.1    1  0.50 0
    0.5 0.1    1  0.70 1
    0.5 0.1    1  0.99 1
    0.3 0.4    10 0.01 0
    0.3 0.4    10 0.30 0
    0.3 0.4    10 0.50 2
    0.3 0.4    10 0.70 3
    0.3 0.4    10 0.90 4
    0.3 0.4    10 0.99 6
    0.7 0.05   20 0.01 0
    0.7 0.05   20 0.10 11
    0.7 0.05   20 0.50 14
    0.7 0.05   20 0.90 17
    0.7 0.05   20 0.99 18
    0.1 0.3    15 0.01 0
    0.1 0.3    15 0.35 0
    0.1 0.3    15 0.90 3
    0.1 0.3    15 0.99 4
    0.6 1.0e-4 5  0.10 2
    0.6 1.0e-4 5  0.30 2
    0.6 1.0e-4 5  0.70 4
    0.6 1.0e-4 5  0.99 5))

(t/deftest zero-inflated-binomial-mv
  (t/are [mu sigma bd mean-v var-v]
      (let [dist (sut/distribution :zibi {:mu mu :sigma sigma :bd bd})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.5 0.1    1  0.450000000000 0.247500000000
    0.3 0.4    10 1.800000000000 3.420000000000
    0.7 0.05   20 13.300000000000 13.300000000000
    0.1 0.3    15 1.050000000000 1.417500000000
    0.6 1.0e-4 5  2.999700000000 1.200779910000))

;; reference values from R's gamlss.dist package: dZABI/pZABI/qZABI

(t/deftest zero-adjusted-binomial
  (t/testing "registered under both :zero-adjusted-binomial and :zabi keys"
    (let [d1 (sut/distribution :zero-adjusted-binomial {:mu 0.3 :sigma 0.4 :bd 10})
          d2 (sut/distribution :zabi {:mu 0.3 :sigma 0.4 :bd 10})]
      (t/is (m/delta-eq (sut/pdf d1 3) (sut/pdf d2 3)))
      (t/is (m/delta-eq (sut/cdf d1 3) (sut/cdf d2 3)))))
  (t/testing "defaults match R's ZABI() defaults (mu=0.5, sigma=0.1, bd=1)"
    (let [dist (sut/distribution :zabi)]
      (t/is (m/delta-eq 0.100000000000 (sut/pdf dist 0)))
      (t/is (m/delta-eq 0.900000000000 (sut/pdf dist 1)))))
  (t/are [mu sigma bd vd d vp p]
      (let [dist (sut/distribution :zabi {:mu mu :sigma sigma :bd bd})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.5 0.1    1  0 0.100000000000 0  0.100000000000
    0.5 0.1    1  1 0.900000000000 1  1.000000000000
    0.3 0.4    10 0 0.400000000000 0  0.400000000000
    0.3 0.4    10 3 0.164750554593 3  0.783655226668
    0.3 0.4    10 7 0.005558015378 7  0.999018029936
    0.3 0.4    10 10 0.000003645928 10 1.000000000000
    0.7 0.05   20 0 0.050000000000 0  0.050000000000
    0.7 0.05   20 5 0.000035520280 5  0.050040792988
    0.7 0.05   20 13 0.156048885962 13 0.422390678389
    0.7 0.05   20 20 0.000758026530 20 1.000000000000
    0.1 0.3    15 0 0.300000000000 0  0.300000000000
    0.1 0.3    15 1 0.302485377617 1  0.602485377617
    0.1 0.3    15 7 0.000244178050 7  0.999970359956
    0.1 0.3    15 15 0.000000000000 15 1.000000000000
    0.6 1.0e-4 5  0 0.000100000000 0  0.000100000000
    0.6 1.0e-4 5  2 0.232760426770 2  0.310447235694
    0.6 1.0e-4 5  4 0.261855480116 4  0.921443355965
    0.6 1.0e-4 5  5 0.078556644035 5  1.000000000000))

(t/deftest zero-adjusted-binomial-icdf
  (t/are [mu sigma bd p vq]
      (let [dist (sut/distribution :zabi {:mu mu :sigma sigma :bd bd})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.5 0.1    1  0.01 0
    0.5 0.1    1  0.05 0
    0.5 0.1    1  0.50 1
    0.5 0.1    1  0.99 1
    0.3 0.4    10 0.01 0
    0.3 0.4    10 0.30 0
    0.3 0.4    10 0.50 2
    0.3 0.4    10 0.70 3
    0.3 0.4    10 0.90 4
    0.3 0.4    10 0.99 6
    0.7 0.05   20 0.01 0
    0.7 0.05   20 0.10 11
    0.7 0.05   20 0.50 14
    0.7 0.05   20 0.90 17
    0.7 0.05   20 0.99 18
    0.1 0.3    15 0.01 0
    0.1 0.3    15 0.35 1
    0.1 0.3    15 0.90 3
    0.1 0.3    15 0.99 5
    0.6 1.0e-4 5  0.10 2
    0.6 1.0e-4 5  0.30 2
    0.6 1.0e-4 5  0.70 4
    0.6 1.0e-4 5  0.99 5))

(t/deftest zero-adjusted-binomial-mv
  (t/are [mu sigma bd mean-v var-v]
      (let [dist (sut/distribution :zabi {:mu mu :sigma sigma :bd bd})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.5 0.1    1  0.900000000000 0.090000000000
    0.3 0.4    10 1.852323555764 3.422494601089
    0.7 0.05   20 13.300000000464 13.299999994296
    0.1 0.3    15 1.322236839855 1.425058154982
    0.6 1.0e-4 5  3.030734723569 1.119145095487))

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

;; reference values for the mixture distribution derived analytically as weighted
;; combinations of the component pdf/cdf, and cross-checked with R's dnorm/pnorm/dpois/
;; ppois/dexp/pexp and numerical root finding (uniroot) for the icdf

(t/deftest mixture
  (t/testing "no mixing (default), matches the base normal distribution"
    (let [dist (sut/distribution :mixture)
          nd (sut/distribution :normal)]
      (t/is (m/delta-eq (sut/pdf nd 0.7) (sut/pdf dist 0.7)))
      (t/is (m/delta-eq (sut/cdf nd 0.7) (sut/cdf dist 0.7)))
      (t/is (m/delta-eq (sut/mean nd) (sut/mean dist)))
      (t/is (m/delta-eq (sut/variance nd) (sut/variance dist)))))
  (t/testing "two normals, equal weights: N(0,1) and N(3,1)"
    (let [dist (sut/distribution :mixture {:distrs [(sut/distribution :normal)
                                                     (sut/distribution :normal {:mu 3.0})]
                                            :weights [0.5 0.5]})]
      (t/is (m/delta-eq 0.201687064406685 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.129517595665892 (sut/pdf dist 1.5)))
      (t/is (m/delta-eq 0.201687064406685 (sut/pdf dist 3.0)))
      (t/is (m/delta-eq 0.250674949015815 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.5 (sut/cdf dist 1.5)))
      (t/is (m/delta-eq 0.749325050984185 (sut/cdf dist 3.0)))
      (t/is (m/delta-eq 1.5 (sut/icdf dist 0.5) 1.0e-4))
      (t/is (m/delta-eq 1.5 (sut/mean dist)))
      (t/is (m/delta-eq 3.25 (sut/variance dist)))))
  (t/testing "two normals, unequal weights: N(0,1) w=0.3 and N(3,1) w=0.7"
    (let [dist (sut/distribution :mixture {:distrs [(sut/distribution :normal)
                                                     (sut/distribution :normal {:mu 3.0})]
                                            :weights [0.3 0.7]})]
      (t/is (m/delta-eq 0.110384893914975 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.185576797117357 (sut/pdf dist 2.0)))
      (t/is (m/delta-eq 0.268328516184288 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 0.404233638167566 (sut/cdf dist 2.0)))
      (t/is (m/delta-eq 1.28155122988272 (sut/icdf dist 0.3) 1.0e-4))
      (t/is (m/delta-eq 2.1 (sut/mean dist)))
      (t/is (m/delta-eq 2.89 (sut/variance dist)))))
  (t/testing "three heterogeneous components: N(0,1) w=0.2, N(2,0.5) w=0.3, Exp(mean=1) w=0.5"
    (let [dist (sut/distribution :mixture {:distrs [(sut/distribution :normal)
                                                     (sut/distribution :normal {:mu 2.0 :sd 0.5})
                                                     (sut/distribution :exponential)]
                                            :weights [0.2 0.3 0.5]})]
      (t/is (m/delta-eq 0.0704139573845688 (sut/pdf dist -0.5)))
      (t/is (m/delta-eq 0.264728445397463 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.189730594122149 (sut/pdf dist 2.5)))
      (t/is (m/delta-eq 0.061707593740669 (sut/cdf dist -0.5)))
      (t/is (m/delta-eq 0.491154268212441 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 0.910118991443458 (sut/cdf dist 2.5)))
      (t/is (m/delta-eq 1.03364756884511 (sut/icdf dist 0.5) 1.0e-4))
      (t/is (m/delta-eq 1.1 (sut/mean dist)))
      (t/is (m/delta-eq 1.265 (sut/variance dist)))))
  (t/testing "mixed continuous/discrete components: N(0,1) and Poisson(3), equal weights"
    (let [dist (sut/distribution :mixture {:distrs [(sut/distribution :normal)
                                                     (sut/distribution :poisson {:p 3.0})]
                                            :weights [0.5 0.5]})]
      (t/is (sut/continuous? dist) "continuous? is true when at least one component is continuous")
      (t/is (m/delta-eq 0.195665964811368 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.520246509769999 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 1.5 (sut/mean dist))))))

;; reference values pinned from the deterministic KDE-integration output itself (pdf/cdf/icdf/mean/variance
;; do not depend on rng, only sample does); a symmetric, known dataset is used so mean/variance/cdf-at-mean
;; also have simple closed-form expectations to check against independently

(t/deftest continuous-distribution
  (t/testing "registered under both :continuous-distribution and :kde keys"
    (let [data [1.0 2.0 2.0 3.0 3.0 3.0 4.0 4.0 5.0]
          opts {:data data :kde :gaussian :bandwidth 0.5 :steps 2000}
          d1 (sut/distribution :continuous-distribution opts)
          d2 (sut/distribution :kde opts)]
      (t/is (m/delta-eq (sut/pdf d1 3.0) (sut/pdf d2 3.0)))
      (t/is (m/delta-eq (sut/cdf d1 3.0) (sut/cdf d2 3.0)))))
  (t/testing "symmetric dataset around 3.0: mean/variance match sample statistics exactly"
    (let [data [1.0 2.0 2.0 3.0 3.0 3.0 4.0 4.0 5.0]
          d (sut/distribution :continuous-distribution {:data data :kde :gaussian :bandwidth 0.5 :steps 2000})]
      (t/is (m/delta-eq 3.0 (sut/mean d)))
      (t/is (m/delta-eq 1.5 (sut/variance d)))
      (t/testing "density is symmetric around the mean"
        (t/is (m/delta-eq (sut/pdf d 1.0) (sut/pdf d 5.0)))
        (t/is (m/delta-eq (sut/pdf d 2.0) (sut/pdf d 4.0))))
      (t/testing "cdf at the (symmetric) mean is exactly 0.5"
        (t/is (m/delta-eq 0.5 (sut/cdf d 3.0))))
      (t/testing "cdf is monotonically increasing and bounded in [0, 1]"
        (t/is (< (sut/cdf d 1.0) (sut/cdf d 2.0) (sut/cdf d 3.0) (sut/cdf d 4.0) (sut/cdf d 5.0)))
        (t/is (m/delta-eq 0.0 (sut/cdf d -100.0)))
        (t/is (m/delta-eq 1.0 (sut/cdf d 100.0))))
      (t/testing "icdf inverts cdf"
        (t/is (m/delta-eq 3.0 (sut/icdf d (sut/cdf d 3.0)) 1.0e-3))
        (t/is (m/delta-eq (sut/lower-bound d) (sut/icdf d 0.0) 0.1) "icdf(0.0) is near the lower bound of the (KDE-padded) support")
        (t/is (m/delta-eq 0.5 (sut/cdf d (sut/icdf d 0.5)) 1.0e-6) "icdf(0.5) round-trips through cdf"))
      (t/testing "lower/upper bounds bracket the data range with room for the kernel tails"
        (t/is (< (sut/lower-bound d) 1.0))
        (t/is (> (sut/upper-bound d) 5.0)))
      (t/testing "exact pinned values from the deterministic KDE/integration pipeline"
        (t/is (m/delta-eq 0.11273904805708312 (sut/pdf d 1.0)))
        (t/is (m/delta-eq 0.31401297060190670 (sut/pdf d 3.0)))
        (t/is (m/delta-eq 0.06062201459400837 (sut/cdf d 1.0)))
        (t/is (m/delta-eq 0.93937798540599220 (sut/cdf d 5.0)))
        (t/is (m/delta-eq 1.30918444437901240 (sut/icdf d 0.1)))
        (t/is (m/delta-eq 4.69081555562098450 (sut/icdf d 0.9))))))
  (t/testing "default data [-1 0 1] is symmetric around 0"
    (let [d (sut/distribution :continuous-distribution)]
      (t/is (m/delta-eq 0.0 (sut/mean d)))
      (t/is (m/delta-eq 0.5 (sut/cdf d 0.0) 1.0e-3))))
  (t/testing "sample draws fall within the (padded) support"
    (let [data [1.0 2.0 2.0 3.0 3.0 3.0 4.0 4.0 5.0]
          d (sut/distribution :continuous-distribution {:data data :kde :gaussian :bandwidth 0.5 :steps 2000})]
      (dotimes [_ 20]
        (let [s (sut/sample d)]
          (t/is (and (> s (sut/lower-bound d)) (< s (sut/upper-bound d)))))))))

;; reference values from R's dnbinom/pnbinom/qnbinom (size=r, prob=p)

(t/deftest negative-binomial
  (t/testing "defaults match r=20, p=0.5"
    (let [dist (sut/distribution :negative-binomial)]
      (t/is (m/delta-eq 9.536743164062499E-7 (sut/pdf dist 0)))))
  (t/are [r p vd d vp p2]
      (let [dist (sut/distribution :negative-binomial {:r r :p p})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p2 (sut/cdf dist vp))))
    20  0.5 0  0.000000953674316 0  0.000000953674316
    20  0.5 10 0.018654400482774 10 0.049368573352694
    20  0.5 20 0.062685343809790 20 0.562685343809790
    20  0.5 30 0.016743659701021 30 0.940539773720282
    5   0.3 0  0.002430000000000 0  0.002430000000000
    5   0.3 5  0.051459672600000 5  0.150268332600000
    5   0.3 10 0.068710126992507 10 0.484508940773157
    5   0.3 20 0.020603304706786 20 0.909528081445864
    1   0.9 0  0.900000000000000 0  0.900000000000000
    1   0.9 1  0.090000000000000 1  0.990000000000000
    1   0.9 2  0.009000000000000 2  0.999000000000000
    1   0.9 5  0.000009000000000 5  0.999999000000000
    3.5 0.7 0  0.286974389101188 0  0.286974389101188
    3.5 0.7 1  0.301323108556247 1  0.588297497657435
    3.5 0.7 3  0.111866204051507 3  0.903556799984409
    3.5 0.7 8  0.001672023403819 8  0.998980430640122
    50  0.8 0  0.000014272476927 0  0.000014272476927
    50  0.8 5  0.014443792322110 5  0.024539000702158
    50  0.8 12 0.101840881495356 12 0.525415515532731
    50  0.8 20 0.017305206002223 20 0.969691915248056))

(t/deftest negative-binomial-icdf
  (t/are [r p p2 vq]
      (let [dist (sut/distribution :negative-binomial {:r r :p p})]
        (m/delta-eq vq (sut/icdf dist p2)))
    20  0.5 0.01 8
    20  0.5 0.15 14
    20  0.5 0.30 16
    20  0.5 0.50 19
    20  0.5 0.70 23
    20  0.5 0.85 27
    20  0.5 0.95 31
    20  0.5 0.99 37
    5   0.3 0.01 1
    5   0.3 0.15 5
    5   0.3 0.30 8
    5   0.3 0.50 11
    5   0.3 0.70 14
    5   0.3 0.85 18
    5   0.3 0.95 23
    5   0.3 0.99 30
    1   0.9 0.01 0
    1   0.9 0.15 0
    1   0.9 0.30 0
    1   0.9 0.50 0
    1   0.9 0.70 0
    1   0.9 0.85 0
    1   0.9 0.95 1
    1   0.9 0.99 1
    3.5 0.7 0.01 0
    3.5 0.7 0.15 0
    3.5 0.7 0.30 1
    3.5 0.7 0.50 1
    3.5 0.7 0.70 2
    3.5 0.7 0.85 3
    3.5 0.7 0.95 4
    3.5 0.7 0.99 6
    50  0.8 0.01 4
    50  0.8 0.15 8
    50  0.8 0.30 10
    50  0.8 0.50 12
    50  0.8 0.70 14
    50  0.8 0.85 17
    50  0.8 0.95 19
    50  0.8 0.99 23))

(t/deftest negative-binomial-mv
  (t/are [r p mean-v var-v]
      (let [dist (sut/distribution :negative-binomial {:r r :p p})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    20  0.5 20.000000000000 40.000000000000
    5   0.3 11.666666666667 38.888888888889
    1   0.9 0.111111111111  0.123456790123
    3.5 0.7 1.500000000000  2.142857142857
    50  0.8 12.500000000000 15.625000000000))

;; reference values from R's extraDistr package: dhcauchy/phcauchy/qhcauchy

(t/deftest half-cauchy
  (t/testing "defaults match mu=0, scale=1"
    (let [dist (sut/distribution :half-cauchy)]
      (t/is (m/delta-eq 0.636619772367581 (sut/pdf dist 0.0)))))
  (t/testing "density and cdf are zero below mu, icdf(0) is mu"
    (let [dist (sut/distribution :half-cauchy {:mu 2.0 :scale 3.0})]
      (t/is (m/delta-eq 0.0 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 1.0)))
      (t/is (m/delta-eq 2.0 (sut/icdf dist 0.0)))))
  (t/testing "mean and variance are undefined (NaN), as for the parent Cauchy distribution"
    (let [dist (sut/distribution :half-cauchy)]
      (t/is (Double/isNaN (sut/mean dist)))
      (t/is (Double/isNaN (sut/variance dist)))))
  (t/are [mu scale vx d vp p]
      (let [dist (sut/distribution :half-cauchy {:mu mu :scale scale})]
        (and (m/delta-eq d (sut/pdf dist vx))
             (m/delta-eq p (sut/cdf dist vp))))
    0.0 1.0 0   0.636619772367581 0   0.000000000000000
    0.0 1.0 0.5 0.509295817894065 0.5 0.295167235300867
    0.0 1.0 1   0.318309886183791 1   0.500000000000000
    0.0 1.0 2   0.127323954473516 2   0.704832764699133
    0.0 1.0 5   0.024485375860292 5   0.874334083621998
    2.0 3.0 2   0.212206590789194 2   0.000000000000000
    2.0 3.0 3   0.190985931710274 3   0.204832764699133
    2.0 3.0 5   0.106103295394597 5   0.500000000000000
    2.0 3.0 10  0.026162456398668 10  0.771599497560184
    0.0 0.5 0   1.273239544735163 0   0.000000000000000
    0.0 0.5 0.25 1.018591635788130 0.25 0.295167235300867
    0.0 0.5 1   0.254647908947033 1   0.704832764699133
    0.0 0.5 3   0.034411879587437 3   0.894863086577493))

(t/deftest half-cauchy-icdf
  (t/are [mu scale p vq]
      (let [dist (sut/distribution :half-cauchy {:mu mu :scale scale})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.0 1.0 0.01 0.015709255323665
    0.0 1.0 0.25 0.414213562373095
    0.0 1.0 0.50 1.000000000000000
    0.0 1.0 0.75 2.414213562373095
    0.0 1.0 0.99 63.656741162871697
    2.0 3.0 0.01 2.047127765970995
    2.0 3.0 0.25 3.242640687119285
    2.0 3.0 0.50 5.000000000000000
    2.0 3.0 0.75 9.242640687119284
    2.0 3.0 0.99 192.970223488615090
    0.0 0.5 0.01 0.007854627661832
    0.0 0.5 0.25 0.207106781186548
    0.0 0.5 0.50 0.500000000000000
    0.0 0.5 0.75 1.207106781186547
    0.0 0.5 0.99 31.828370581435848))

;; reference values derived analytically from the empirical pmf formula: p(v) = sum of weights of
;; occurrences of v, normalized to sum to 1.0; mean/variance are the usual weighted moments

(t/deftest integer-discrete-distribution
  (t/testing "registered under :integer-discrete-distribution, :integer-discrete-distribution key, and :integer-discrete alias"
    (let [d1 (sut/distribution :integer-discrete-distribution {:data [1 2 2 3 3 3]})
          d2 (sut/distribution :integer-discrete {:data [1 2 2 3 3 3]})]
      (t/is (= (sut/pdf d1 2) (sut/pdf d2 2)))))
  (t/testing "repeated values accumulate probability mass (implicit equal weights)"
    (let [dist (sut/distribution :integer-discrete-distribution {:data [1 2 2 3 3 3]})]
      (t/is (= [1.0 3.0] [(sut/lower-bound dist) (sut/upper-bound dist)]))
      (t/is (m/delta-eq (/ 1.0 6) (sut/pdf dist 1)))
      (t/is (m/delta-eq (/ 1.0 3) (sut/pdf dist 2)))
      (t/is (m/delta-eq 0.5 (sut/pdf dist 3)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 4)) "no mass outside the support")
      (t/is (= [0.0 (/ 1.0 6) 0.5 1.0 1.0] (mapv #(sut/cdf dist %) [0 1 2 3 4])))
      (t/is (= [1 1 2 2 3 3] (mapv #(sut/icdf dist %) [0.0 (/ 1.0 6) 0.17 0.5 0.51 1.0])))
      (t/is (m/delta-eq (/ 7.0 3) (sut/mean dist)))
      (t/is (m/delta-eq (/ 5.0 9) (sut/variance dist)))))
  (t/testing "explicit probabilities (need not sum to 1.0, normalized internally)"
    (let [dist (sut/distribution :integer-discrete-distribution {:data [10 20 30] :probabilities [1 2 3]})]
      (t/is (m/delta-eq (/ 1.0 6) (sut/pdf dist 10)))
      (t/is (m/delta-eq (/ 1.0 3) (sut/pdf dist 20)))
      (t/is (m/delta-eq 0.5 (sut/pdf dist 30)))
      (t/is (m/delta-eq 23.333333333333332 (sut/mean dist)))
      (t/is (m/delta-eq 55.555555555555660 (sut/variance dist)))))
  (t/testing "degenerate default distribution concentrated on 1"
    (let [dist (sut/distribution :integer-discrete-distribution)]
      (t/is (m/delta-eq 1.0 (sut/pdf dist 1)))
      (t/is (m/delta-eq 1.0 (sut/mean dist)))
      (t/is (m/delta-eq 0.0 (sut/variance dist)))))
  (t/testing "sample always returns a value from the support"
    (let [dist (sut/distribution :integer-discrete-distribution {:data [1 2 2 3 3 3]})]
      (dotimes [_ 20]
        (t/is (contains? #{1 2 3} (sut/sample dist)))))))

(t/deftest real-discrete-distribution
  (t/testing "registered under :real-discrete-distribution key and :real-discrete alias"
    (let [d1 (sut/distribution :real-discrete-distribution {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]})
          d2 (sut/distribution :real-discrete {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]})]
      (t/is (= (sut/pdf d1 5.0) (sut/pdf d2 5.0)))))
  (t/testing "explicit probabilities over a real-valued support"
    (let [dist (sut/distribution :real-discrete-distribution {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]})]
      (t/is (= [0.0 10.0] [(sut/lower-bound dist) (sut/upper-bound dist)]))
      (t/is (m/delta-eq 0.5 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.3 (sut/pdf dist 5.0)))
      (t/is (m/delta-eq 0.2 (sut/pdf dist 10.0)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 3.0)) "no mass outside the support")
      (t/is (= [0.0 0.5 0.5 0.8 1.0 1.0] (mapv #(sut/cdf dist %) [-1.0 0.0 2.0 5.0 10.0 11.0])))
      (t/is (= [0.0 0.0 0.0 0.0 5.0 5.0 10.0 10.0]
               (mapv #(sut/icdf dist %) [0.0 0.4 0.49 0.5 0.79 0.8 0.99 1.0])))
      (t/is (m/delta-eq 3.5 (sut/mean dist)))
      (t/is (m/delta-eq 15.25 (sut/variance dist)))))
  (t/testing "degenerate default distribution concentrated on 1.0"
    (let [dist (sut/distribution :real-discrete-distribution)]
      (t/is (m/delta-eq 1.0 (sut/pdf dist 1.0)))
      (t/is (m/delta-eq 1.0 (sut/mean dist)))
      (t/is (m/delta-eq 0.0 (sut/variance dist)))))
  (t/testing "sample always returns a value from the support"
    (let [dist (sut/distribution :real-discrete-distribution {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]})]
      (dotimes [_ 20]
        (t/is (contains? #{0.0 5.0 10.0} (sut/sample dist)))))))

;; reference values computed independently from the Kolmogorov distribution's series formula
;; cdf(x) = 1 - 2*sum_{k>=1} (-1)^(k-1)*exp(-2*k^2*x^2), pdf = d/dx of that; icdf values cross-checked
;; against the well-known asymptotic two-sided KS-test critical constants (0.90 -> 1.22385,
;; 0.95 -> 1.35810, 0.99 -> 1.62762); mean = sqrt(pi/2)*ln(2) is a known closed form for this distribution

(t/deftest kolmogorov
  (t/testing "parameter-free, no rng needed"
    (let [dist (sut/distribution :kolmogorov)]
      (t/is (m/delta-eq 0.0 (sut/lower-bound dist)))
      (t/is (= ##Inf (sut/upper-bound dist)))))
  (t/testing "pdf/cdf are zero at and below 0"
    (let [dist (sut/distribution :kolmogorov)]
      (t/is (m/delta-eq 0.0 (sut/pdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.0)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist -1.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist -1.0)))))
  (t/testing "mean matches the known closed form sqrt(pi/2)*ln(2)"
    (let [dist (sut/distribution :kolmogorov)]
      (t/is (m/delta-eq (* (m/sqrt (/ m/PI 2.0)) (m/log 2.0)) (sut/mean dist)))))
  (t/are [vx d vp p]
      (let [dist (sut/distribution :kolmogorov)]
        (and (m/delta-eq d (sut/pdf dist vx))
             (m/delta-eq p (sut/cdf dist vp))))
    0.2    0.000000000153242 0.2    0.000000000000505
    0.5    0.639582850940456 0.5    0.036054756335125
    0.8    1.627024345636592 0.8    0.455857588425802
    1.0    1.071948558356942 1.0    0.730000328322645
    1.3581 0.271605259750462 1.3581 0.950000369568333
    1.5    0.133307227419880 1.5    0.977782037383475
    2.0    0.005367402045630 2.0    0.999329074744220))

(t/deftest kolmogorov-icdf
  (t/testing "matches the well-known two-sided KS-test asymptotic critical constants"
    (let [dist (sut/distribution :kolmogorov)]
      (t/is (m/delta-eq 1.22384787022 (sut/icdf dist 0.90) 1.0e-6))
      (t/is (m/delta-eq 1.35809863932 (sut/icdf dist 0.95) 1.0e-6))
      (t/is (m/delta-eq 1.62762361152 (sut/icdf dist 0.99) 1.0e-6))))
  (t/are [p vq]
      (let [dist (sut/distribution :kolmogorov)]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-6))
    0.10 0.571173265106
    0.25 0.676447691503
    0.50 0.827573555190
    0.75 1.019184720254
    0.90 1.223847870217
    0.95 1.358098639323
    0.99 1.627623611519))

;; reference values derived analytically from the reciprocal (log-uniform) formulas:
;; pdf(x) = 1/(x*ln(b/a)), cdf(x) = ln(x/a)/ln(b/a), icdf(p) = a*(b/a)^p,
;; mean = (b-a)/ln(b/a), var = (b^2-a^2)/(2*ln(b/a)) - mean^2

(t/deftest reciprocal
  (t/testing "defaults match a=1, b=10"
    (let [dist (sut/distribution :reciprocal)]
      (t/is (m/delta-eq 0.434294481903 (sut/pdf dist 1.0)))))
  (t/testing "density and cdf are zero outside [a, b]"
    (let [dist (sut/distribution :reciprocal {:a 1.0 :b 10.0})]
      (t/is (m/delta-eq 0.0 (sut/pdf dist 0.5)))
      (t/is (m/delta-eq 0.0 (sut/pdf dist 11.0)))
      (t/is (m/delta-eq 0.0 (sut/cdf dist 0.5)))
      (t/is (m/delta-eq 1.0 (sut/cdf dist 11.0)))))
  (t/are [a b vx d vp p]
      (let [dist (sut/distribution :reciprocal {:a a :b b})]
        (and (m/delta-eq d (sut/pdf dist vx))
             (m/delta-eq p (sut/cdf dist vp))))
    1.0 10.0  1.0   0.434294481903 1.0   0.000000000000
    1.0 10.0  3.25  0.133629071355 3.25  0.511883360979
    1.0 10.0  5.5   0.078962633073 5.5   0.740362689494
    1.0 10.0  7.75  0.056037997665 7.75  0.889301702506
    1.0 10.0  10.0  0.043429448190 10.0  1.000000000000
    0.5 5.0   0.5   0.868588963807 0.5   0.000000000000
    0.5 5.0   1.625 0.267258142710 1.625 0.511883360979
    0.5 5.0   2.75  0.157925266147 2.75  0.740362689494
    0.5 5.0   3.875 0.112075995330 3.875 0.889301702506
    0.5 5.0   5.0   0.086858896381 5.0   1.000000000000
    2.0 200.0 2.0   0.108573620476 2.0   0.000000000000
    2.0 200.0 51.5  0.004216451281 51.5  0.705388616689
    2.0 200.0 101.0 0.002149972683 101.0 0.851645689059
    2.0 200.0 150.5 0.001442838810 150.5 0.938253252133
    2.0 200.0 200.0 0.001085736205 200.0 1.000000000000))

(t/deftest reciprocal-icdf
  (t/are [a b p vq]
      (let [dist (sut/distribution :reciprocal {:a a :b b})]
        (m/delta-eq vq (sut/icdf dist p)))
    1.0 10.0  0.1  1.258925411794
    1.0 10.0  0.25 1.778279410039
    1.0 10.0  0.5  3.162277660168
    1.0 10.0  0.75 5.623413251903
    1.0 10.0  0.9  7.943282347243
    0.5 5.0   0.1  0.629462705897
    0.5 5.0   0.25 0.889139705019
    0.5 5.0   0.5  1.581138830084
    0.5 5.0   0.75 2.811706625952
    0.5 5.0   0.9  3.971641173621
    2.0 200.0 0.1  3.169786384922
    2.0 200.0 0.25 6.324555320337
    2.0 200.0 0.5  20.000000000000
    2.0 200.0 0.75 63.245553203368
    2.0 200.0 0.9  126.191468896039))

(t/deftest reciprocal-mv
  (t/are [a b mean-v var-v]
      (let [dist (sut/distribution :reciprocal {:a a :b b})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    1.0 10.0  3.908650337129  6.220029396270
    0.5 5.0   1.954325168565  1.555007349068
    2.0 200.0 42.995153708422 2493.927282139787))

;; reference values from R's `gamlss.dist` package: dexGAUS/pexGAUS/qexGAUS(x, mu, sigma, nu),
;; using fastmath's (mu, sigma, tau) parameter names for `nu`. pdf/cdf are exact closed-form
;; expressions (erf/erfc-based), so agreement is essentially to full double precision;
;; mean = mu+tau, variance = sigma^2+tau^2.

(t/deftest ex-gaussian
  (t/testing "defaults are mu=0, sigma=1, tau=1"
    (let [dist (sut/distribution :ex-gaussian nil)]
      (t/is (m/delta-eq 0.261578291865 (sut/pdf dist 0.0)))))
  (t/testing "support is the whole real line"
    (let [dist (sut/distribution :ex-gaussian {:mu 0.0 :sigma 1.0 :tau 1.0})]
      (t/is (Double/isInfinite (sut/lower-bound dist)))
      (t/is (Double/isInfinite (sut/upper-bound dist)))
      (t/is (neg? (sut/lower-bound dist)))
      (t/is (pos? (sut/upper-bound dist)))))
  (t/testing "exgaus is the same distribution under a renamed (nu instead of tau) parameter"
    (let [d1 (sut/distribution :ex-gaussian {:mu 1.0 :sigma 2.0 :tau 3.0})
          d2 (sut/distribution :exgaus {:mu 1.0 :sigma 2.0 :nu 3.0})]
      (doseq [x [-3.0 0.0 1.5 5.0]]
        (t/is (m/delta-eq (sut/pdf d1 x) (sut/pdf d2 x)))
        (t/is (m/delta-eq (sut/cdf d1 x) (sut/cdf d2 x))))
      (t/is (m/delta-eq (sut/mean d1) (sut/mean d2)))
      (t/is (m/delta-eq (sut/variance d1) (sut/variance d2)))))
  (t/testing "as tau -> 0 the exponential component vanishes and it approaches normal(mu, sigma)"
    (let [eg (sut/distribution :ex-gaussian {:mu 1.0 :sigma 2.0 :tau 1.0e-6})
          nd (sut/distribution :normal {:mu 1.0 :sd 2.0})]
      (doseq [x [-3.0 0.0 1.0 3.0 5.0]]
        (t/is (m/delta-eq (sut/pdf eg x) (sut/pdf nd x) 1.0e-4))
        (t/is (m/delta-eq (sut/cdf eg x) (sut/cdf nd x) 1.0e-4)))))
  (t/testing "pdf/cdf/icdf are never NaN or throw, across extreme/infinite/negative-zero inputs"
    (doseq [mu [-100.0 0.0 100.0]
            sigma [0.1 1.0 50.0]
            tau [0.1 1.0 50.0]]
      (let [dist (sut/distribution :ex-gaussian {:mu mu :sigma sigma :tau tau})]
        (doseq [x [##-Inf ##Inf -1e300 1e300 1e-300 0.0 -0.0]]
          (t/is (not (Double/isNaN (sut/pdf dist x))) (str "pdf mu=" mu " sigma=" sigma " tau=" tau " x=" x))
          (t/is (not (Double/isNaN (sut/cdf dist x))) (str "cdf mu=" mu " sigma=" sigma " tau=" tau " x=" x)))
        (doseq [p [0.0 1e-10 0.5 (- 1.0 1e-10) 1.0]]
          (t/is (not (Double/isNaN (sut/icdf dist p))) (str "icdf mu=" mu " sigma=" sigma " tau=" tau " p=" p)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##-Inf)))
        (t/is (m/delta-eq 0.0 (sut/pdf dist ##Inf)))
        (t/is (m/delta-eq 0.0 (sut/cdf dist ##-Inf)))
        (t/is (m/delta-eq 1.0 (sut/cdf dist ##Inf))))))
  (t/are [mu sigma tau vd d vp p]
      (let [dist (sut/distribution :ex-gaussian {:mu mu :sigma sigma :tau tau})]
        (and (m/delta-eq d (sut/pdf dist vd) 1.0e-8)
             (m/delta-eq p (sut/cdf dist vp) 1.0e-8)))
    0.0 1.0 1.0    -5 0.000000241410 -5 0.000000045242
    0.0 1.0 1.0    -1 0.101959017701 -1 0.056696236231
    0.0 1.0 1.0    0  0.261578291865 0  0.238421708135
    0.0 1.0 1.0    1  0.303265329856 1  0.538079416212
    0.0 1.0 1.0    2  0.187729387930 2  0.789520480122
    0.0 1.0 1.0    5  0.011108644703 5  0.988891068646
    0.0 1.0 1.0    10 0.000074851830 10 0.999925148170
    2.0 1.5 3.0    -5 0.000000464113 -5 0.000000138288
    2.0 1.5 3.0    -1 0.006375705254 -1 0.003623016186
    2.0 1.5 3.0    0  0.024554783065 0  0.017546870532
    2.0 1.5 3.0    1  0.064139095438 1  0.060075251233
    2.0 1.5 3.0    2  0.116539611573 2  0.150381165280
    2.0 1.5 3.0    5  0.129670878276 5  0.588237233223
    2.0 1.5 3.0    10 0.026245004137 10 0.921264939377
    -1.0 0.5 2.0   -5 0.000000000000 -5 0.000000000000
    -1.0 0.5 2.0   -1 0.207016051474 -1 0.085967897052
    -1.0 0.5 2.0   0  0.300357814900 0  0.376534238251
    -1.0 0.5 2.0   1  0.189761814383 1  0.620444699991
    -1.0 0.5 2.0   2  0.115106535360 2  0.769786928293
    -1.0 0.5 2.0   5  0.025683739784 5  0.948632520433
    -1.0 0.5 2.0   10 0.002108249745 10 0.995783500510
    0.0 2.0 0.5    -5 0.005273809488 -5 0.003572760582
    0.0 2.0 0.5    -1 0.149677461728 -1 0.233698807862
    0.0 2.0 0.5    0  0.188821282604 0  0.405589358698
    0.0 2.0 0.5    1  0.187698537373 1  0.597613192587
    0.0 2.0 0.5    2  0.147403870521 2  0.767642810808
    0.0 2.0 0.5    5  0.018082743012 5  0.984748963168
    0.0 2.0 0.5    10 0.000010338802 10 0.999994543948
    5.0 1.0 10.0   -5 0.000000000000 -5 0.000000000000
    5.0 1.0 10.0   -1 0.000000000097 -1 0.000000000015
    5.0 1.0 10.0   0  0.000000028140 0  0.000000005251
    5.0 1.0 10.0   1  0.000003097185 1  0.000000699390
    5.0 1.0 10.0   2  0.000131267471 2  0.000037223320
    5.0 1.0 10.0   5  0.046247878529 5  0.037521214712
    5.0 1.0 10.0   10 0.060957061520 10 0.390429098148))

(t/deftest ex-gaussian-icdf
  (t/are [mu sigma tau p vq]
      (let [dist (sut/distribution :ex-gaussian {:mu mu :sigma sigma :tau tau})]
        (m/delta-eq vq (sut/icdf dist p) 1.0e-4))
    0.0 1.0 1.0  0.05 -1.068877781234
    0.0 1.0 1.0  0.25 0.043767710283
    0.0 1.0 1.0  0.50 0.875798343698
    0.0 1.0 1.0  0.75 1.803377906257
    0.0 1.0 1.0  0.95 3.494166385603
    2.0 1.5 3.0  0.05 0.832225187865
    2.0 1.5 3.0  0.25 2.748544759117
    2.0 1.5 3.0  0.50 4.365199740911
    2.0 1.5 3.0  0.75 6.531432946076
    2.0 1.5 3.0  0.95 11.362196819588
    -1.0 0.5 2.0 0.05 -1.203882363153
    -1.0 0.5 2.0 0.25 -0.404231151309
    -1.0 0.5 2.0 0.50 0.448206117510
    -1.0 0.5 2.0 0.75 1.835088719744
    -1.0 0.5 2.0 0.95 5.053964547108
    0.0 2.0 0.5  0.05 -2.874692992742
    0.0 2.0 0.5  0.25 -0.893079333408
    0.0 2.0 0.5  0.50 0.491001789388
    0.0 2.0 0.5  0.75 1.882744322136
    0.0 2.0 0.5  0.95 3.904850418096
    5.0 1.0 10.0 0.05 5.246901474998
    5.0 1.0 10.0 0.25 7.926142579730
    5.0 1.0 10.0 0.50 11.981471805599
    5.0 1.0 10.0 0.75 18.912943611199
    5.0 1.0 10.0 0.95 35.007322735523))

(t/deftest ex-gaussian-mv
  (t/are [mu sigma tau mean-v var-v]
      (let [dist (sut/distribution :ex-gaussian {:mu mu :sigma sigma :tau tau})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.0 1.0 1.0   1.0  2.0
    2.0 1.5 3.0   5.0  11.25
    -1.0 0.5 2.0  1.0  4.25
    0.0 2.0 0.5   0.5  4.25
    5.0 1.0 10.0  15.0 101.0))

;; reference values from R's extraDistr package: dbbinom/pbbinom (qbbinom is absent from
;; extraDistr, so icdf reference values were derived by manual search over the pbbinom cdf table);
;; mean = n*alpha/(alpha+beta), var = n*alpha*beta*(alpha+beta+n) / ((alpha+beta)^2*(alpha+beta+1))

(t/deftest beta-binomial
  (t/testing "defaults match alpha=0.5, beta=0.5, n=10"
    (let [dist (sut/distribution :beta-binomial)]
      (t/is (m/delta-eq 0.176197052001953 (sut/pdf dist 0)))))
  (t/are [alpha beta n vd d vp p]
      (let [dist (sut/distribution :beta-binomial {:alpha alpha :beta beta :n n})]
        (and (m/delta-eq d (sut/pdf dist vd))
             (m/delta-eq p (sut/cdf dist vp))))
    0.5 0.5 10 0  0.176197052001953 0  0.176197052001952
    0.5 0.5 10 3  0.065460205078125 3  0.408035278320311
    0.5 0.5 10 5  0.060562133789063 5  0.530281066894529
    0.5 0.5 10 8  0.073642730712891 8  0.731067657470700
    0.5 0.5 10 10 0.176197052001953 10 1.000000000000000
    2.0 3.0 20 0  0.021739130434783 0  0.021739130434783
    2.0 3.0 20 5  0.076792772444946 5  0.328063241106720
    2.0 3.0 20 10 0.068322981366460 10 0.704968944099381
    2.0 3.0 20 15 0.031620553359684 15 0.940711462450596
    2.0 3.0 20 20 0.001976284584980 20 1.000000000000000
    5.0 1.0 15 0  0.000064499484004 0  0.000064499484004
    5.0 1.0 15 5  0.008126934984520 5  0.016253869969040
    5.0 1.0 15 10 0.064563983488132 10 0.193691950464396
    5.0 1.0 15 15 0.250000000000000 15 1.000000000000000))

(t/deftest beta-binomial-icdf
  (t/are [alpha beta n p vq]
      (let [dist (sut/distribution :beta-binomial {:alpha alpha :beta beta :n n})]
        (m/delta-eq vq (sut/icdf dist p)))
    0.5 0.5 10 0.10 0
    0.5 0.5 10 0.25 1
    0.5 0.5 10 0.50 5
    0.5 0.5 10 0.75 9
    0.5 0.5 10 0.90 10
    0.5 0.5 10 0.99 10
    2.0 3.0 20 0.10 2
    2.0 3.0 20 0.25 4
    2.0 3.0 20 0.50 8
    2.0 3.0 20 0.75 11
    2.0 3.0 20 0.90 14
    2.0 3.0 20 0.99 18
    5.0 1.0 15 0.10 9
    5.0 1.0 15 0.25 11
    5.0 1.0 15 0.50 13
    ;; p=0.75 lands exactly on the cdf(14) boundary (R: 0.750000000000001, fastmath:
    ;; 0.749999999999997) - a floating-point tie broken differently by each
    ;; implementation's incomplete-beta computation, not a real discrepancy
    5.0 1.0 15 0.90 15
    5.0 1.0 15 0.99 15))

(t/deftest beta-binomial-mv
  (t/are [alpha beta n mean-v var-v]
      (let [dist (sut/distribution :beta-binomial {:alpha alpha :beta beta :n n})]
        (and (m/delta-eq mean-v (sut/mean dist))
             (m/delta-eq var-v (sut/variance dist))))
    0.5 0.5 10 5.0  13.75
    2.0 3.0 20 8.0  20.0
    5.0 1.0 15 12.5 6.25))

;; reference values derived analytically from the Dirichlet pdf/mean/covariance formulas:
;; pdf(x) = (1/B(alpha)) * prod(xi^(alpha_i - 1)), mean_i = alpha_i / A,
;; var_i = alpha_i*(A - alpha_i) / (A^2*(A+1)), cov_ij = -alpha_i*alpha_j / (A^2*(A+1)),
;; where A = sum(alpha). For 2 dimensions this reduces exactly to the Beta(alpha0, alpha1) distribution.

(t/deftest dirichlet
  (t/testing "default alpha=[1 1] is uniform on the simplex (Beta(1,1))"
    (let [dist (sut/distribution :dirichlet)]
      (t/is (m/delta-eq 1.0 (sut/pdf dist [0.3 0.7])))
      (t/is (m/delta-eq 1.0 (sut/pdf dist [0.9 0.1])))
      (t/is (= [0.5 0.5] (sut/means dist)))
      (t/is (m/delta-eq 0.083333333333333 (get-in (sut/covariance dist) [0 0])))
      (t/is (m/delta-eq -0.083333333333333 (get-in (sut/covariance dist) [0 1])))))
  (t/testing "alpha=[2 3] (2d, equals Beta(2,3))"
    (let [dist (sut/distribution :dirichlet {:alpha [2 3]})]
      (t/is (m/delta-eq 1.5 (sut/pdf dist [0.5 0.5])))
      (t/is (= [0.4 0.6] (sut/means dist)))
      (t/is (m/delta-eq 0.04 (get-in (sut/covariance dist) [0 0])))
      (t/is (m/delta-eq 0.04 (get-in (sut/covariance dist) [1 1])))
      (t/is (m/delta-eq -0.04 (get-in (sut/covariance dist) [0 1])))))
  (t/testing "alpha=[1 1 1] (3d symmetric, uniform on the 2-simplex)"
    (let [dist (sut/distribution :dirichlet {:alpha [1 1 1]})]
      (t/is (m/delta-eq 2.0 (sut/pdf dist [(/ 1.0 3) (/ 1.0 3) (/ 1.0 3)])))
      (t/is (= [(/ 1.0 3) (/ 1.0 3) (/ 1.0 3)] (sut/means dist)))
      (t/is (m/delta-eq 0.055555555555556 (get-in (sut/covariance dist) [0 0])))
      (t/is (m/delta-eq -0.027777777777778 (get-in (sut/covariance dist) [0 1])))))
  (t/testing "alpha=[2 2 2] (3d symmetric, concentrated)"
    (let [dist (sut/distribution :dirichlet {:alpha [2 2 2]})]
      (t/is (m/delta-eq 4.444444444444440 (sut/pdf dist [(/ 1.0 3) (/ 1.0 3) (/ 1.0 3)])))
      (t/is (m/delta-eq 0.031746031746032 (get-in (sut/covariance dist) [0 0])))
      (t/is (m/delta-eq -0.015873015873016 (get-in (sut/covariance dist) [0 1])))))
  (t/testing "integer alpha=4 means symmetric 4d Dirichlet with concentration 1.0 each"
    (let [dist (sut/distribution :dirichlet {:alpha 4})]
      (t/is (= 4 (sut/dimensions dist)))
      (t/is (m/delta-eq 6.0 (sut/pdf dist [0.25 0.25 0.25 0.25])))
      (t/is (= [0.25 0.25 0.25 0.25] (sut/means dist)))
      (t/is (m/delta-eq 0.0375 (get-in (sut/covariance dist) [0 0])))
      (t/is (m/delta-eq -0.0125 (get-in (sut/covariance dist) [0 1])))))
  (t/testing "sample returns a probability vector summing to 1.0"
    (let [dist (sut/distribution :dirichlet {:alpha [2 3 1]})
          s (sut/sample dist)]
      (t/is (= 3 (count s)))
      (t/is (m/delta-eq 1.0 (reduce + s) 1.0e-9))
      (t/is (every? #(and (>= % 0.0) (<= % 1.0)) s)))))

;; ---------------------------------------------------------------------------
;; Custom RNG / `set-seed!` reproducibility, exercised generically across
;; *every* distribution registered under the `distribution` multimethod
;; (i.e. every key in `(methods sut/distribution)`, canonical names and
;; gamlss-style aliases alike).
;;
;; Two independent ways of controlling a distribution's randomness are
;; supported throughout `fastmath.random`/`fastmath.random.distributions`:
;;
;; 1. "external rng"  - an `rng` object is created explicitly (via [[sut/rng]])
;;    and passed in under the `:rng` key when constructing the distribution.
;;    Reseeding *that rng object itself* (via [[sut/set-seed!]]) must make
;;    the distribution reproduce its earlier sample sequence, since sampling
;;    ultimately draws from that same, shared, mutable generator.
;; 2. "internal rng" - no `:rng` is supplied; the distribution creates its own,
;;    private generator under the hood. Reseeding *the distribution itself*
;;    (`(sut/set-seed! dist seed)`, dispatching to `prot/set-seed!`, which every
;;    distribution implements - directly, or via the shared `->distribution`/
;;    `ssj-continuous`/`multinomial`/`dirichlet` helpers in `distributions.clj`)
;;    must likewise reproduce the earlier sequence.
;;
;; A few `:data`-driven distributions have no sensible parameter-free default
;; (they'd otherwise default to a single point/category, making sampling
;; trivially seed-invariant); `distribution-overrides` supplies just enough of
;; a non-degenerate parameter map for those so that the "different seeds give
;; different samples" sanity check below is meaningful, without otherwise
;; duplicating the per-distribution parameter tables already exercised by the
;; tests above.

(def ^:private distribution-overrides
  "Per-`:key` parameter overrides (merged with the caller-supplied `:rng`) used
  only by the generic reproducibility tests below, for `:data`-driven
  distributions whose parameter-free defaults are a single point/category."
  {:empirical {:data (vec (map double (range 1 101))) :bin-count 20}
   :enumerated-real {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]}
   :enumerated-int {:data [1 2 3] :probabilities [0.2 0.3 0.5]}
   :real-discrete-distribution {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]}
   :real-discrete {:data [0.0 5.0 10.0] :probabilities [0.5 0.3 0.2]}
   :integer-discrete-distribution {:data [1 2 2 3 3 3]}
   :integer-discrete {:data [1 2 2 3 3 3]}
   :categorical-distribution {:data [:a :b :c] :probabilities [0.2 0.3 0.5]}
   :categorical {:data [:a :b :c] :probabilities [0.2 0.3 0.5]}})

(def ^:private all-distribution-keys
  (keys (methods sut/distribution)))

(t/deftest custom-rng-reseed-external
  (t/testing "every distribution, constructed with an externally-supplied, non-default RNG (:isaac): reseeding that same rng object reproduces the sample sequence"
    (doseq [k all-distribution-keys]
      (let [params (get distribution-overrides k)
            seed 424242
            rng (sut/rng :isaac seed)
            dist (sut/distribution k (assoc params :rng rng))
            xs1 (vec (sut/->seq dist 20))]
        (sut/set-seed! rng seed)
        (t/is (= xs1 (vec (sut/->seq dist 20)))
              (str "reseeding the external rng did not reproduce samples for " k))))))

(t/deftest custom-rng-reseed-internal
  (t/testing "every distribution, constructed with no :rng (its own private generator): reseeding the distribution itself reproduces the sample sequence"
    (doseq [k all-distribution-keys]
      (let [params (get distribution-overrides k)
            seed 13579
            dist (sut/distribution k params)]
        (sut/set-seed! dist seed)
        (let [xs1 (vec (sut/->seq dist 20))]
          (sut/set-seed! dist seed)
          (t/is (= xs1 (vec (sut/->seq dist 20)))
                (str "reseeding the distribution's own rng did not reproduce samples for " k)))))))

(t/deftest custom-rng-different-seeds-differ
  (t/testing "sanity check: two distributions (same params) seeded differently do NOT produce the same sample sequence - guards against a reseed test that trivially 'passes' because sampling secretly ignores the rng"
    (doseq [k (remove #{:constant} all-distribution-keys)]
      (let [params (get distribution-overrides k)
            d1 (sut/distribution k (assoc params :rng (sut/rng :isaac 1)))
            d2 (sut/distribution k (assoc params :rng (sut/rng :isaac 2)))]
        (t/is (not= (vec (sut/->seq d1 20)) (vec (sut/->seq d2 20)))
              (str "different seeds produced identical samples for " k)))))
  (t/testing ":constant is legitimately seed-invariant by design (degenerate/Dirac distribution)"
    (let [d1 (sut/distribution :constant {:value 5.0 :rng (sut/rng :isaac 1)})
          d2 (sut/distribution :constant {:value 5.0 :rng (sut/rng :isaac 2)})]
      (t/is (= (vec (sut/->seq d1 5)) (vec (sut/->seq d2 5)) [5.0 5.0 5.0 5.0 5.0])))))

(t/deftest custom-rng-multiple-algorithms
  (t/testing "reseed-reproducibility holds across different RNG algorithms, not just :isaac, for a representative sample of distribution kinds (scalar continuous, scalar discrete, multivariate, data-driven, meta/combinator)"
    (doseq [rng-name [:jdk :mersenne :well19937c :well44497b]
            k [:normal :poisson :multi-normal :dirichlet :continuous-distribution :truncated :mixture]]
      (let [seed 777
            rng (sut/rng rng-name seed)
            dist (sut/distribution k {:rng rng})
            xs1 (vec (sut/->seq dist 15))]
        (sut/set-seed! rng seed)
        (t/is (= xs1 (vec (sut/->seq dist 15)))
              (str rng-name " reseed did not reproduce samples for " k))))))

;; ---------------------------------------------------------------------------
;; Golden-value regression tests: unlike the reseed tests above (which only
;; check that the *same* rng/seed reproduces *itself* - a tautology w.r.t.
;; whatever the current implementation happens to do, so it would still
;; "pass" even if some future refactor silently changed how many draws a
;; `sample` call consumes, their order, etc., as long as it changed
;; consistently both times), these pin the *exact* sample sequence produced
;; by a fixed rng algorithm + seed against a literal, hardcoded expected
;; value, captured once from the current implementation. A change to any of
;; these numbers is not necessarily a bug, but it does mean the RNG
;; consumption pattern for that distribution changed, and is worth a second
;; look.
;;
;; Covers *every* distribution key registered under `distribution` (all of
;; `all-distribution-keys`, canonical names and aliases alike - an alias
;; necessarily gets the same golden values as its canonical name, since it
;; dispatches to the exact same underlying function), not just a selection.
;; Values were captured directly via `(vec (sut/->seq dist 5))` with an
;; `:isaac`-seeded rng (seed `42`; merged with `distribution-overrides` for
;; the handful of `:data`-driven distributions that need it - same overrides
;; used by the reseed tests above), and are compared with plain `=` (exact,
;; bit-for-bit equality of primitive doubles/longs/vectors - safe here since
;; Clojure's double literal reader round-trips the shortest-representation
;; decimal strings printed by `->seq` exactly).

(def ^:private golden-values-isaac-seed-42
  "Exact `(vec (sut/->seq dist 5))` output, for every distribution key
  registered under `distribution` (`all-distribution-keys`), constructed
  with `(sut/rng :isaac 42)` under `:rng` (merged with `distribution-overrides`
  for the handful of :data-driven distributions that need it). Captured once
  from the current implementation; see `custom-rng-golden-values` below."
  {:anderson-darling
   [0.5059259287441009
    1.062802526792838
    1.4240409420706273
    0.564601413510125
    0.4550188823911598]
   :anderson-darling-quick
   [0.5059259287441009
    1.062802526792838
    1.4240409420706273
    0.564601413510125
    0.4550188823911598]
   :bb [2 8 9 3 1]
   :bernoulli [0 1 1 0 0]
   :beta
   [0.3817005118702333
    0.7303121967085704
    0.3212486429861899
    0.8944112669509152
    0.5490390312564815]
   :beta-binomial [2 8 9 3 1]
   :beta-noncentral
   [0.4454898824073701
    0.6880047901609451
    0.7590337119179655
    0.49344896619352496
    0.3864108560723533]
   :beta-symmetrical
   [0.3886854365418788
    0.6375713394242244
    0.7157867999774323
    0.43572726379469356
    0.3321728921161551]
   :binomial [9 11 12 9 9]
   :categorical [:b :c :c :b :b]
   :categorical-distribution [:b :c :c :b :b]
   :cauchy
   [-0.5671188023683934
    0.7320755164354125
    1.409490103729508
    -0.3106624066718934
    -0.9526731380984537]
   :chi
   [0.43410352920398315
    1.0389020210141662
    1.2918313137226247
    0.5303373067991106
    0.3288267458897465]
   :chi-squared
   [0.1884458736135823
    1.0793174092684488
    1.6688281430965568
    0.28125765898294713
    0.10812702864649432]
   :chi-squared-noncentral
   [0.4849011797716233
    2.382648677786544
    3.46748114606517
    0.7076347914229134
    0.284313678064495]
   :constant [0.0 0.0 0.0 0.0 0.0]
   :continuous-distribution
   [-0.5652642775322254
    0.5469059942799166
    0.812824344124532
    -0.4013080493841728
    -0.7643898871122008]
   :cramer-von-mises
   [0.11152151655219623
    0.20623605654953592
    0.24477031794585336
    0.12416196389156708
    0.09993737654754145]
   :dirichlet
   [[0.20441252563844983 0.7955874743615502]
    [0.23841041387096973 0.7615895861290303]
    [0.14648471154904386 0.8535152884509561]
    [0.5357933214156829 0.4642066785843172]
    [0.7534637912934224 0.24653620870657755]]
   :empirical
   [33.84240172746576
    70.594195545867
    80.98267480926177
    41.272110724817814
    26.740954989960905]
   :enumerated-int [2 3 3 2 2]
   :enumerated-real [0.0 5.0 10.0 0.0 0.0]
   :erlang
   [1.1956135155793854
    2.4446285559739276
    3.018408560257909
    1.3883017273476357
    0.9822547591945716]
   :ex-gaussian
   [0.1746423776720124
    2.0328143366077898
    0.3420685194177471
    1.8572703020061199
    2.912234346435775]
   :exgaus
   [0.1746423776720124
    2.0328143366077898
    0.3420685194177471
    1.8572703020061199
    2.912234346435775]
   :exponential
   [1.036294129336298
    0.4022993872419831
    0.60716886287668
    1.3096348408057512
    0.7240005967687287]
   :f
   [0.3393064249228097
    3.886432890717298
    9.845075789006323
    0.5424055412007495
    0.1835957672227515]
   :f-noncentral
   [0.8258953766419942
    8.513653885357316
    21.26712060211158
    1.2856234278464367
    0.4619340458850594]
   :fatigue-life
   [0.6564713358972849
    1.6850093586439723
    2.294297294417468
    0.7849772550022391
    0.5275771207005328]
   :fishers-noncentral-hypergeometric [2 3 3 2 2]
   :folded-normal
   [0.43410352920397494
    1.0389020210141566
    1.2918313137226245
    0.5303373067991218
    0.3288267458897569]
   :frechet
   [0.9163559193846709
    2.8166326638009735
    4.573037867872715
    1.103705642857667
    0.7375135222313798]
   :gamma
   [2.336631390950336
    5.570934838902748
    3.199803738098394
    6.975203012896114
    0.7998836569119738]
   :ge
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :generalized-exponential
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :generalized-extreme-value
   [-0.08735043149701814
    1.0355420806047824
    1.5201777253922815
    0.09867328445706379
    -0.3044708557072609]
   :generalized-gamma
   [0.731044980701355
    1.1924281892500581
    1.3867818740419124
    0.8074659145799746
    0.642761706427969]
   :generalized-half-logistic
   [0.6986724366718834
    1.7391167913653134
    2.217298652131714
    0.8571313242364567
    0.5273150689722524]
   :generalized-hyperbolic
   [-0.24272025759956362
    -0.8282255176557877
    -0.317679539345054
    -0.20575065178854876
    -0.3902434520482083]
   :generalized-inverse-gaussian
   [1.475617600429787
    3.219031336335554
    4.093338584000922
    1.7260194513542255
    1.2084640629391679]
   :generalized-logistic
   [-0.6821270291168965
    0.852778606743133
    1.4088495967909693
    -0.3883194882769483
    -1.057887547396456]
   :generalized-normal
   [-0.3981318507046895
    0.5146652979584333
    0.9343754359812362
    -0.2128914960312695
    -0.6627601619575108]
   :generalized-pareto
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :geometric [0 1 2 0 0]
   :gev
   [-0.08735043149701814
    1.0355420806047824
    1.5201777253922815
    0.09867328445706379
    -0.3044708557072609]
   :gg
   [0.731044980701355
    1.1924281892500581
    1.3867818740419124
    0.8074659145799746
    0.642761706427969]
   :gh
   [-0.24272025759956362
    -0.8282255176557877
    -0.317679539345054
    -0.20575065178854876
    -0.3902434520482083]
   :ghl
   [0.6986724366718834
    1.7391167913653134
    2.217298652131714
    0.8571313242364567
    0.5273150689722524]
   :gig
   [1.475617600429787
    3.219031336335554
    4.093338584000922
    1.7260194513542255
    1.2084640629391679]
   :gnd
   [-0.3981318507046895
    0.5146652979584333
    0.9343754359812362
    -0.2128914960312695
    -0.6627601619575108]
   :gpd
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :gumbel
   [0.8252991370059637
    3.0710841612095647
    4.040355450784563
    1.1973465689141276
    0.3910582885854782]
   :half-cauchy
   [0.5825001501544679
    1.9714037868158543
    3.1376863751383484
    0.7364818674943616
    0.42848076645177974]
   :half-logistic
   [0.6986724366718834
    1.7391167913653134
    2.217298652131714
    0.8571313242364567
    0.5273150689722524]
   :half-normal
   [0.43410352920398304
    1.0389020210141664
    1.2918313137226247
    0.5303373067991105
    0.32882674588974625]
   :hyperbolic-secant
   [-0.3440457725471681
    0.43210304197028
    0.727965609304286
    -0.194723311989008
    -0.5395412589492902]
   :hypergeometric [12 14 14 12 11]
   :hypoexponential
   [0.40915200214740227
    1.2078124785181783
    1.6275226165402286
    0.5177191883123478
    0.29801979512081445]
   :hypoexponential-equal
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :integer-discrete [2 3 3 2 2]
   :integer-discrete-distribution [2 3 3 2 2]
   :inverse-gamma
   [0.9163559193846706
    2.816632663800974
    4.573037867872715
    1.103705642857667
    0.7375135222313798]
   :inverse-gaussian
   [0.46757348548088035
    1.0883959603858404
    1.4648492596701388
    0.546133292430545
    0.38725988930610317]
   :johnson-sb
   [0.39556249420070944
    0.6289488539312967
    0.7015090848507907
    0.43962300097622964
    0.3428969024719658]
   :johnson-sl
   [0.65443075653889
    1.6950462506180797
    2.3501857150329712
    0.7845129292281703
    0.5218312069474558]
   :johnson-su
   [-0.4368073926717723
    0.5525459234670957
    0.9623437131396185
    -0.2450816612112171
    -0.6972486330517929]
   :kde
   [-0.5652642775322254
    0.5469059942799166
    0.812824344124532
    -0.4013080493841728
    -0.7643898871122008]
   :kolmogorov
   [0.7280269575530616
    0.9740595251495796
    1.0769678404433802
    0.7686521693539642
    0.681184313756874]
   :kolmogorov-smirnov
   [0.6678933685970441
    0.8505748468104958
    0.90179221571917
    0.7020609575307257
    0.6288566770260979]
   :kolmogorov-smirnov+
   [0.33578673719408814
    0.7011496936209916
    0.80358443143834
    0.40412191506145145
    0.25771335405219586]
   :kolmogorov-smirnov-quick
   [0.6678933685970441
    0.8505748468104958
    0.90179221571917
    0.7020609575307257
    0.6288566770260979]
   :laplace
   [-0.39813185070499796
    0.514665297959185
    0.93437543598128
    -0.21289149603034835
    -0.6627601619582715]
   :levy
   [1.0793874784394621
    6.789851910498397
    16.166292345314098
    1.436725912774866
    0.7806406794730248]
   :log-logistic
   [0.5055405484912872
    2.3461568506199835
    4.091246113141351
    0.678195633093517
    0.3471884553751189]
   :log-normal
   [1.7638051896776283
    5.602312706564407
    4.240595699592065
    0.7629541600690262
    2.5794542874870054]
   :logarithmic [1 1 2 1 1]
   :logistic
   [-0.6821270291168965
    0.8527786067431333
    1.408849596790969
    -0.38831948827694823
    -1.057887547396456]
   :mixture
   [-0.1998108701214906
    -0.6304138764984518
    -0.2889830525969167
    -0.10785633972323627
    -0.24976602353894722]
   :multi-normal
   [[-0.43252648520466763 0.7231794958020387]
    [0.44470375456214967 -1.2705573280489608]
    [-0.05242213991714392 1.080815091524896]
    [0.7171157560619612 -0.2154858707394434]
    [-1.466248077438233 -0.6801959988403614]]
   :multinomial [[4 6] [6 4] [4 6] [8 2] [5 5]]
   :nakagami
   [0.6396499059235821
    1.0990052225868794
    1.2757439459975797
    0.7195270587777487
    0.545911893163992]
   :nbi [0 1 2 0 0]
   :nbii [0 1 2 0 0]
   :negative-binomial [17 23 25 18 16]
   :normal
   [-0.43252648520466763
    0.7231794958020387
    0.44470375456214967
    -1.2705573280489608
    -0.05242213991714392]
   :normal-inverse-gaussian
   [-0.13662934418594996
    -0.4658810569572944
    -0.17983473202192113
    -0.12126086153080746
    -0.22196819750629498]
   :pareto
   [2.9780806959686137
    1.4262289623712692
    1.244424307984781
    2.47450079476125
    3.8802801029762137]
   :pascal [17 23 25 18 16]
   :pearson-6
   [0.5055405484912872
    2.3461568506199835
    4.091246113141351
    0.678195633093517
    0.3471884553751189]
   :poisson [0 1 0 0 0]
   :power
   [0.5794710840016852
    0.8373468180037419
    0.8964287096241061
    0.6357058400403849
    0.5076547587211173]
   :rayleigh
   [0.9046015721278035
    1.5542280904160306
    1.8041743909839898
    1.017564925017903
    0.772036003204204]
   :real-discrete [0.0 5.0 10.0 0.0 0.0]
   :real-discrete-distribution [0.0 5.0 10.0 0.0 0.0]
   :reciprocal
   [2.166639902550565
    5.025157679656955
    6.361864741914581
    2.5358403912059075
    1.8101449542753494]
   :t
   [-0.5671188023675264
    0.7320755164353214
    1.4094901037295366
    -0.310662406672404
    -0.952673138273561]
   :t-noncentral
   [0.6894996659378079
    2.6663962730648274
    4.255383102902819
    0.9330290683997322
    0.415047791194716]
   :triangular
   [-0.18050413400177778
    0.22688900359779074
    0.3732375752142444
    -0.10097617933510628
    -0.2820677552133546]
   :truncated
   [-0.43252648520466763
    0.7231794958020387
    0.44470375456214967
    -1.2705573280489608
    -0.05242213991714392]
   :uniform-int [1442193049 2058843292 1859088621 1382605351 1735690423]
   :uniform-real
   [0.33578673719408814
    0.7011496936209916
    0.80358443143834
    0.40412191506145145
    0.25771335405219586]
   :von-mises
   [-0.5004524794387464
    0.6266695622324583
    1.043566719881831
    -0.28435373978310785
    -0.7793157246280061]
   :wallenius-noncentral-hypergeometric [2 3 3 2 2]
   :watson-g
   [0.455004946833131
    0.5673901837747928
    0.6109047354469331
    0.4745544739829435
    0.43202630695806465]
   :watson-u
   [0.05576075827609812
    0.10311802827476796
    0.12238515897292668
    0.06208098194578354
    0.049968688273770726]
   :weibull
   [0.40915200214804676
    1.2078124785191304
    1.6275226165412253
    0.5177191883133454
    0.2980197951217608]
   :zabb [1.0 1.0 1.0 1.0 1.0]
   :zabi [1.0 1.0 1.0 1.0 1.0]
   :zaga
   [0.5164680434699724
    0.2730450462344595
    0.697200678970517
    0.12139884042455018
    3.6523833357778965]
   :zaig
   [1.0883959603858404
    0.546133292430545
    0.767968506909367
    1.2640084543381096
    0.40021219275560665]
   :zanbi [1.0 2.0 2.0 1.0 0.0]
   :zap [4.0 6.0 7.0 4.0 3.0]
   :zero-adjusted-beta-binomial [1.0 1.0 1.0 1.0 1.0]
   :zero-adjusted-binomial [1.0 1.0 1.0 1.0 1.0]
   :zero-adjusted-gamma
   [0.5164680434699724
    0.2730450462344595
    0.697200678970517
    0.12139884042455018
    3.6523833357778965]
   :zero-adjusted-inverse-gaussian
   [1.0883959603858404
    0.546133292430545
    0.767968506909367
    1.2640084543381096
    0.40021219275560665]
   :zero-adjusted-negative-binomial [1.0 2.0 2.0 1.0 0.0]
   :zero-adjusted-poisson [4.0 6.0 7.0 4.0 3.0]
   :zero-inflated-beta-binomial [1 0 1 1 0]
   :zero-inflated-binomial [1 0 1 1 0]
   :zero-inflated-negative-binomial [1 0 0 4 1]
   :zero-inflated-poisson [8 2 6 2 6]
   :zero-inflated-poisson2 [8 2 6 2 6]
   :zibb [1 0 1 1 0]
   :zibi [1 0 1 1 0]
   :zinbi [1 0 0 4 1]
   :zip [8 2 6 2 6]
   :zip2 [8 2 6 2 6]
   :zipf [1 1 1 1 1]})

(t/deftest custom-rng-golden-values
  (t/testing "every distribution key registered under `distribution` matches its pinned golden sample sequence (:isaac, seed=42)"
    (doseq [[k expected] golden-values-isaac-seed-42]
      (let [params (get distribution-overrides k)
            dist (sut/distribution k (assoc params :rng (sut/rng :isaac 42)))]
        (t/is (= expected (vec (sut/->seq dist (count expected))))
              (str "golden value mismatch for " k)))))
  (t/testing "every registered distribution key has a golden-value entry (and vice versa) - keeps this table honest as new distributions are added"
    (t/is (= (set all-distribution-keys) (set (keys golden-values-isaac-seed-42)))))
  (t/testing "other RNG algorithms (not :isaac), to make sure the pinning isn't an :isaac-only artifact"
    (t/is (= [4.704771863704436 3.5741221790904723 1.9613231073484556 1.2166069925080893 0.7683138172400454]
             (vec (sut/->seq (sut/distribution :gamma {:shape 2.0 :scale 1.5 :rng (sut/rng :mersenne 100)}) 5))))
    (t/is (= [1 4 2 3 5]
             (vec (sut/->seq (sut/distribution :poisson {:p 4.0 :rng (sut/rng :well19937c 7)}) 5))))))

