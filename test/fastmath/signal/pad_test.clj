(ns fastmath.signal.pad-test
  (:require [fastmath.signal.pad :as sut]
            [fastmath.signal :as sig]
            [clojure.test :as t]
            [clojure.java.io :as io]
            [clojure.edn :as edn]
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


;; Reference values from `pywt.pad`, generated by
;; `utils/fastmath/dev/generate_pad_reference.py` into
;; `test/resources/signal/pad_reference.edn`: every method, side, signal
;; lengths 1, 2, 3, 5, 6 and target lengths up to `3*len+7` (repeated extension).

(def ^:private pad-reference
  (delay (edn/read-string (slurp (io/resource "signal/pad_reference.edn")))))

(t/deftest padding-vs-pywavelets
  (doseq [{:keys [method side N signal result]} @pad-reference]
    (t/is (= result (vec ((pad-methods method) (double-array signal) N side)))
          (str method " " side " N=" N " signal=" signal))))

(t/deftest padding-guards
  (t/testing "pad-position rejects unknown side"
    (t/is (thrown? clojure.lang.ExceptionInfo (sut/pad-position :middle 10 5))))
  (doseq [[mkey f] pad-methods]
    (t/testing (str mkey " guards")
      (t/is (thrown? clojure.lang.ExceptionInfo (f (double-array 0) 4 :both)) "empty signal")
      (t/is (thrown? clojure.lang.ExceptionInfo (f (double-array [1 2 3]) 2 :both)) "N smaller than length")
      (t/is (thrown? clojure.lang.ExceptionInfo (f (double-array [1 2 3]) 5 :middle)) "unknown side"))
    (t/testing (str mkey " N equal to length returns a copy")
      (let [s (double-array [1 2 3])
            res (f s 3 :both)]
        (t/is (= [1.0 2.0 3.0] (vec res)))
        (t/is (not (identical? s res)))))
    (t/testing (str mkey " one sample, no padding (longer paddings in the pywt reference)")
      (t/is (= [2.5] (vec (f (double-array [2.5]) 1 :left)))))))

(t/deftest out-len-pos-placement
  (t/are [side start] (= [5 start (+ start 5)] (rest (sut/out-len-pos (double-array 5 1.0) 12 side)))
    :left 7
    :right 0
    :both 4))
