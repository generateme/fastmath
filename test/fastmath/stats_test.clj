(ns fastmath.stats-test
  (:require [fastmath.stats :as sut]
            [clojure.test :as t]))

(def one {:size 1, :step 0.0, :samples 1, :min 1, :max 1, :bins '([1.0 1])})
(def one2 {:size 1, :step 1.0, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 2])})
(def two {:size 2, :step 0.5, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 1] [1.5 1])})

(t/deftest histogram-tests
  (t/are [in method res] (= res (sut/histogram in method))
    [1 1 1 1] :sqrt  {:size 1, :step 0.0, :samples 4, :min 1.0, :max 1.0, :bins '([1.0 4])}
    [1 1] :sturges {:size 1, :step 0.0, :samples 2, :min 1.0, :max 1.0, :bins '([1.0 2])}
    [1] :rice one
    [1] :doane one
    [1 2] :sqrt one2
    [1 2] :sturges two
    [1 2] :rice two
    [1 2] :doane one2
    [1 2] :scott one2
    [1 2] :freedman-diaconis one2))
