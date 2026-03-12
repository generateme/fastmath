(ns fastmath.signal.pad-test
  (:require [fastmath.signal.pad :as sut]
            [fastmath.signal :as sig]
            [clojure.test :as t]

            [fastmath.vector :as v]))

(def signal-even (double-array [-1 2 0 10]))
(def signal-odd (double-array [-1 2 0 5 10]))

(def pad-methods {:zero sut/zero
                :edge sut/edge
                :periodic sut/periodic
                :symmetric sut/symmetric
                :antisymmetric sut/antisymmetric
                :reflect sut/reflect
                :antireflect sut/antireflect
                :linear sut/linear})

(defn test-paddings
  "test if result, internal padding and api return the same"
  [signal mkey side res]
  (= res
     (seq ((pad-methods mkey) signal 8 side))
     (seq (sig/pad signal 8 mkey side))))

(t/deftest padding

  (t/testing "where to put original singal in padding buffer"
    (t/is (= 90 (sut/pad-position :left 100 10)))
    (t/is (zero? (sut/pad-position :right 100 10)))
    (t/are [pos N] (= pos (sut/pad-position :both N 10))
      42 93
      42 94
      43 95
      43 96
      44 97
      44 98
      45 99
      45 100
      46 101
      46 102
      47 103
      47 104
      48 105))

  (let [mkey :zero]
    (t/testing "zero - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat [0.0 0.0 0.0 0.0] signal-even)
        :right (concat signal-even [0.0 0.0 0.0 0.0])
        :both  (concat [0.0 0.0] signal-even [0.0 0.0])))

    (t/testing "zero - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [0.0 0.0 0.0] signal-odd)
        :right (concat signal-odd [0.0 0.0 0.0])
        :both  (concat [0.0 0.0] signal-odd [0.0]))))

  (let [mkey :edge]
    (t/testing "edge - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat [-1.0 -1.0 -1.0 -1.0] signal-even)
        :right (concat signal-even [10.0 10.0 10.0 10.0])
        :both  (concat [-1.0 -1.0] signal-even [10.0 10.0])))

    (t/testing "edge - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [-1.0 -1.0 -1.0] signal-odd)
        :right (concat signal-odd [10.0 10.0 10.0])
        :both  (concat [-1.0 -1.0] signal-odd [10.0]))))

  (let [mkey :periodic]
    (t/testing "periodic - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat signal-even signal-even)
        :right (concat signal-even signal-even)
        :both  (concat [0.0 10.0] signal-even [-1.0 2.0])))

    (t/testing "periodic - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [0.0 5.0 10.0] signal-odd)
        :right (concat signal-odd [-1.0 2.0 0.0])
        :both  (concat [5.0 10.0] signal-odd [-1.0]))))

  (let [mkey :symmetric]
    (t/testing "symmetric - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat (reverse signal-even) signal-even)
        :right (concat signal-even (reverse signal-even))
        :both  (concat [2.0 -1.0] signal-even [10.0 0.0])))

    (t/testing "symmetric - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [0.0 2.0 -1.0] signal-odd)
        :right (concat signal-odd [10.0 5.0 0.0])
        :both  (concat [2.0 -1.0] signal-odd [10.0]))))

  (let [mkey :reflect]
    (t/testing "reflect - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat [0.0 10.0 0.0 2.0] signal-even)
        :right (concat signal-even [0.0 2.0 -1.0 2.0])
        :both  (concat [0.0 2.0] signal-even [0.0 2.0])))

    (t/testing "reflect - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [5.0 0.0 2.0] signal-odd)
        :right (concat signal-odd [5.0 0.0 2.0])
        :both  (concat [0.0 2.0] signal-odd [5.0]))))

  (let [mkey :antisymmetric]
    (t/testing "antisymmetric - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat (v/sub (reverse signal-even)) signal-even)
        :right (concat signal-even (v/sub (reverse signal-even)))
        :both  (concat [-2.0 1.0] signal-even [-10.0 0.0])))

    (t/testing "antisymmetric - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [0.0 -2.0 1.0] signal-odd)
        :right (concat signal-odd [-10.0 -5.0 0.0])
        :both  (concat [-2.0 1.0] signal-odd [-10.0]))))

  (let [mkey :antireflect]
    (t/testing "antireflect - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat [-22.0 -12.0 -2.0 -4.0]  signal-even)
        :right (concat signal-even [20.0 18.0 21.0 24.0])
        :both  (concat [-2.0 -4.0] signal-even [20.0 18.0])))

    (t/testing "antireflect - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [-7.0 -2.0 -4.0] signal-odd)
        :right (concat signal-odd [15.0 20.0 18.0])
        :both  (concat [-2.0 -4.0] signal-odd [15.0]))))

  (let [mkey :linear]
    (t/testing "linear - even signal"
      (t/are [side res] (test-paddings signal-even mkey side res)
        :left  (concat [-13.0 -10.0 -7.0 -4.0] signal-even)
        :right (concat signal-even [20.0 30.0 40.0 50.0])
        :both  (concat [-7.0 -4.0] signal-even [20.0 30.0])))

    (t/testing "linear - odd signal"
      (t/are [side res] (test-paddings signal-odd mkey side res)
        :left  (concat [-10.0 -7.0 -4.0] signal-odd)
        :right (concat signal-odd [15.0 20.0 25.0])
        :both  (concat [-7.0 -4.0] signal-odd [15.0])))))

