(ns fastmath.interpolation.step
  (:require [fastmath.core :as m])
  (:import [org.apache.commons.math3.analysis.function StepFunction]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(defn step-after
  "Step function."
  [xs ys]
  (let [^StepFunction sf (StepFunction. (m/seq->double-array xs)
                                        (m/seq->double-array ys))]
    (fn ^double [^double x] (.value sf x))))

(defn step-before
  "Step function."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        l (dec (alength xs))]
    (fn ^double [^double x]
      (let [b (java.util.Arrays/binarySearch xs x)
            i (if (neg? b) (m/constrain (dec (- b)) 0 l) b)]
        (Array/aget ys i)))))

(defn step
  "Step function."
  ([xs ys] (step xs ys nil))
  ([xs ys {:keys [^double point]
           :or {point 0.5}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         l (dec (alength xs))
         lst (dec (- l))]
     (fn ^double [^double x]
       (let [b (java.util.Arrays/binarySearch ^doubles xs x) 
             i (cond
                 (== b -1) 0
                 (< b lst) l
                 (< b -1) (let [b-2 (- (- b) 2)
                                b-1 (dec (- b))
                                x1 (Array/aget xs b-2)
                                x2 (Array/aget xs b-1)
                                diff (+ x1 (* point (- x2 x1)))]
                            (if (< x diff) b-2 b-1))
                 :else b)]
         (Array/aget ys i))))))

(m/unuse-primitive-operators)
