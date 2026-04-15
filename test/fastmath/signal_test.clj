(ns fastmath.signal-test
  (:require [fastmath.signal :as sut]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [clojure.test :as t]))

;; periodogram
;; comparing with python/scipy

;; >>> scipy.signal.welch(x=np.asarray([1,2,3,10,-1,-2,3,1]),window=[0,0.75,0.75,0],detrend=False,noverlap=2,scaling="density")
;; (array([0.  , 0.25, 0.5 ]), array([17.83333333, 42.33333333, 24.5       ]))
;; >>> scipy.signal.welch(x=np.asarray([1,2,3,10,-1,-2,3,1]),window=[0,0.75,0.75,0],detrend=False,noverlap=2,scaling="spectrum")
;;(array([0.  , 0.25, 0.5 ]), array([8.91666667, 21.16666667, 12.25      ]))

(t/deftest periodogram
  (let [sig [1,2,3,10,-1,-2,3,1]]
    (t/testing "Welsh periodogram (comparison with scipy)"
      (t/is (v/edelta-eq [17.83333333, 42.33333333, 24.5]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]}))))
      (t/is (v/edelta-eq [8.91666667, 21.16666667, 12.25]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :method :power}))))
      (t/is (v/edelta-eq [15.0 , 15.6, 15.0]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :average :umedian}))))
      (t/is (v/edelta-eq [7.5, 7.8, 7.5]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :method :power
                                                          :average :umedian})))))))
(t/deftest resampling
  (t/testing "even length signal"
    (let [sig [1 2 3 10 -1 2]]
      (t/are [cnt res] (v/edelta-eq res (sut/resample sig cnt))
        1 [2.83333333]
        2 [0.16666667 5.5]
        3 [0.16666667 5.16666667 3.16666667]
        4 [2.83333333 1.32136721 8.16666667 -0.98803387]
        5 [2.83333333 0.27137889 7.59165872 4.03785588 -0.56756016]
        6 sig
        7 [1.0 2.00611241 1.50782978 8.71029852 5.90272448 -1.74568137 2.45204951]
        8 [1.0 1.90587372 1.32136721 5.3937861 10.0 1.45139187 -0.98803387 2.58228164]
        9 [1.0 1.77534578 1.5233184  3.0 9.43579021 7.16147413 -1.0 0.03886401 2.56520747]
        10 [1.0 1.64706163 1.75457671 1.7937036 7.02512757 10. 3.47132473 -1.76009924 0.91563766 2.48600068])))
  
  (t/testing "odd length signal"
    (let [sig [1 2 3 10 -1]]
      (t/are [cnt res] (v/edelta-eq res (sut/resample sig cnt))
        1 [3.0]
        2 [-0.68328157 6.68328157]
        3 [-0.68328157 4.40470422 5.27857735]
        4 [1.0 0.81218754 8.36656315 1.82124931]
        5 sig
        6 [1.0 2.7968157 0.64602959 8.36656315 7.35397041 -2.16337885]
        7 [1.0 3.2183382 0.34969088 4.51567976 10.22038208 4.25635504 -2.56044598]
        8 [1.0 3.40706996 0.81218754 1.87941562 8.36656315 9.32953114 1.82124931 -2.61601671]
        9 [1.0 3.46356369 1.4337983 0.64602959 5.41295943 10.08828259 7.35397041 0.12347689 -2.5220809 ]
        10 [1.0 3.4472136 2.0 0.31671843 3.0 8.36656315 10.0 5.23606798 -1.0 -2.36656315]))))

(t/deftest hilbert
  (t/testing "if real values in hilber transform are same as signal"
    (let [data (map #(repeatedly % r/grand) (range 1 101))]
      (doseq [d data]
        (t/is (v/delta-eq d (map first (sut/hilbert d))))))))

(t/deftest convolution
  (t/testing "full mode"
    (t/are [in1 in2 res] (and (v/edelta-eq res (sut/convolve in1 in2))
                              (v/edelta-eq res (sut/fft-convolve in1 in2)))
      [-1] [2] [-2]
      [-1 1] [2] [-2 2]
      [2] [-1 1] [-2 2]
      [-1 1] [2 3 4 5 6 7] [-2 -1 -1 -1 -1 -1 7])))

(t/deftest correlation
  (t/testing "full mode"
    (t/are [in1 in2 res] (and (v/edelta-eq res (sut/correlate in1 in2))
                              (v/edelta-eq res (sut/fft-correlate in1 in2)))
      [-1] [2] [-2]
      [-1 1] [2] [-2 2]
      [2] [-1 1] [2 -2]
      [-1 1] [2 3 4 5 6 7] [-7 1 1 1 1 1 2])))
