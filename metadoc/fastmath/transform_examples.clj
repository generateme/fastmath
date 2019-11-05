(ns fastmath.transform-examples
  (:require [metadoc.examples :refer :all]
            [fastmath.transform :refer :all]
            [fastmath.vector :as v]
            [fastmath.core :as m]))


(add-examples wavelets-list
  (example "List of wavelets" (sort wavelets-list)))

(add-examples transformer
  (example "Usage" (transformer :packet :discrete-mayer)))

(add-examples forward-1d
  (example-session "Usage"
    (seq (forward-1d (transformer :packet :haar-orthogonal) [-1 8 7 6]))
    (seq (forward-1d (transformer :fast :haar-orthogonal) [-1 8 7 6]))
    (seq (forward-1d (transformer :decomposed-fast :haar-orthogonal) [10 -1 8 7 6]))
    (seq (forward-1d (transformer :orthogonal :sine) [0 -1 8 7]))
    (seq (forward-1d (transformer :standard :hadamard) [-1 8 7 6]))
    (seq (forward-1d (transformer :standard :dft) [-1 8 7 6]))))

(add-examples reverse-1d
  (example-session "Usage"
    (let [t (transformer :packet :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))
    (let [t (transformer :fast :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))
    (let [t (transformer :decomposed-fast :haar-orthogonal)]
      (seq (reverse-1d t (forward-1d t [10 -1 8 7 6]))))
    (let [t (transformer :orthogonal :sine)]
      (seq (reverse-1d t (forward-1d t [0 8 7 6]))))
    (let [t (transformer :standard :hadamard)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))
    (let [t (transformer :standard :dft)]
      (seq (reverse-1d t (forward-1d t [-1 8 7 6]))))))

(add-examples forward-2d
  (example-session "Usage"
    (m/double-double-array->seq (forward-2d (transformer :packet :daubechies-6) [[-1 8] [7 6]]))
    (m/double-double-array->seq (forward-2d (transformer :fast :daubechies-6) [[-1 8] [7 6]]))))

(add-examples reverse-2d
  (example-session "Usage"
    (let [t (transformer :packet :daubechies-20)]
      (m/double-double-array->seq (reverse-2d t (forward-2d t [[-1 8] [7 6]]))))
    (let [t (transformer :fast :daubechies-20)]
      (m/double-double-array->seq (reverse-2d t (forward-2d t [[-1 8] [7 6]]))))))

(add-examples compress
  (example "Compress 1d - explicite usage"
    (let [t (transformer :fast :symlet-5)]
      (seq (reverse-1d t (compress (forward-1d t [1 2 3 4]) 0.3)))))
  (example "Compress 1d - implicit usage"
    (let [t (transformer :fast :symlet-5)]
      (seq (compress t [1 2 3 4] 0.3))))
  (example "Compress 2d - explicit usage"
    (let [t (transformer :fast :symlet-5)]
      (m/double-double-array->seq (reverse-2d t (compress (forward-2d t [[1 2] [3 4]]) 0.5)))))
  (example "Compress 2d - implicit usage"
    (let [t (transformer :fast :symlet-5)]
      (m/double-double-array->seq (compress t [[1 2] [3 4]] 0.5)))))

(add-examples denoise
  (example "Denoise signal - implicit usage"
    (let [t (transformer :packet :haar)]
      (v/approx (seq (denoise t [1 2 3 4 5 6.5 7.5 8] true)))))
  (example "Denoise signal - explicit usage"
    (let [t (transformer :packet :haar)]
      (v/approx (vec (reverse-1d t (denoise (forward-1d t [1 2 3 4 5 6.5 7.5 8]) true)))))))
