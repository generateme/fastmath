(ns fastmath.interpolation.common
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat])
  (:import [fastmath.java Array]
           [org.apache.commons.math3.linear LUDecomposition SingularValueDecomposition
            RealMatrix RealVector SingularMatrixException MatrixUtils]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;;                                v    n=4        
;; ------0---------1--------------2------------3----------
;;  -1   0   -2    1       -3     2   -4       3     -5
;;   0   0    0    1        1     2    2       2      2 

(defn binary-search-id
  "Find array index after binary search snapping to the lower value. The last segment snaps to the previous segment."
  ^long [^long bsres ^long n]
  (if (m/neg? bsres)
    (cond
      (m/== -1 bsres) 0
      (m/> bsres (m/- n)) (m/- (m/- bsres) 2)
      :else (m/- n 2))
    (if (m/>= bsres (m/- n 2)) (m/- n 2) bsres)))

(defn binary-search
  ^long [^doubles arr ^double v]
  (-> (java.util.Arrays/binarySearch arr v)
      (binary-search-id (alength arr))))

;;

;; validation of the data given to interpolators

(defn invalid-input!
  "Throws `ex-info` describing invalid input data of an interpolation method.

  Parameters:

  - `method` (keyword): name of the interpolation method, used in the message
  - `property` (keyword): the violated property, for example `:length-mismatch`, `:too-few-points`, `:not-finite` or `:not-increasing`
  - `message` (string): description of the problem
  - `data` (map): additional data stored in the exception

  Never returns. The `ex-data` of the exception contains `:type` (`:invalid-input`), `:method`, `:property` and the keys of `data`.

  See also [[check-axis!]], [[validate-xy!]]."
  [method property message data]
  (throw (ex-info (str (name method) ": " message)
                  (assoc data :type :invalid-input :method method :property property))))

(defn check-axis!
  "Checks the coordinates of one axis of an interpolator and throws `ex-info` when they are invalid.

  Parameters:

  - `method` (keyword): name of the interpolation method, used in the message
  - `axis` (keyword): name of the axis, for example `:xs`
  - `arr` (double array): the coordinates
  - `opts` (map, optional keys):
      - `:min-points` (long, default: 1): minimal number of points
      - `:sorted?` (boolean, default: true): require strictly increasing coordinates

  Every coordinate has to be finite (not `NaN`, not infinite), also when `:sorted?` is false.

  Returns `nil`. Throws `ex-info` with `:property` set to `:too-few-points`, `:not-finite` or `:not-increasing`.

  See also [[validate-xy!]], [[invalid-input!]]."
  [method axis ^doubles arr {:keys [min-points sorted?] :or {min-points 1 sorted? true}}]
  (let [n (alength arr)
        min-points (long min-points)]
    (when (m/< n min-points)
      (invalid-input! method :too-few-points
                      (str "at least " min-points " points are required in " (name axis) ", got " n)
                      {:axis axis :min-points min-points :count n}))
    (dotimes [i n]
      (let [a (Array/aget arr i)]
        (when-not (m/valid-double? a)
          (invalid-input! method :not-finite
                          (str (name axis) " has to be finite, got " a " at index " i)
                          {:axis axis :index i :value a}))
        (when (and sorted? (m/pos? i) (m/<= a (Array/aget arr (m/dec i))))
          (invalid-input! method :not-increasing
                          (str (name axis) " has to be strictly increasing, got " a " after "
                               (Array/aget arr (m/dec i)) " at index " i)
                          {:axis axis :index i :value a}))))))

(defn validate-xy!
  "Checks the data of a one dimensional interpolator and throws `ex-info` when it is invalid.

  Parameters:

  - `method` (keyword): name of the interpolation method, used in the message
  - `xs` (double array): the nodes
  - `ys` (double array): the values; `NaN` and infinite values are allowed
  - `opts` (map, optional keys): `:min-points` (default: 1) and `:sorted?` (default: true), see [[check-axis!]]

  The number of values has to be equal to the number of nodes.

  Returns `nil`. Throws `ex-info` with `:property` set to `:length-mismatch`, `:too-few-points`, `:not-finite` or `:not-increasing`.

  See also [[check-axis!]], [[invalid-input!]]."
  [method ^doubles xs ^doubles ys opts]
  (when-not (m/== (alength xs) (alength ys))
    (invalid-input! method :length-mismatch
                    (str "xs and ys have different lengths: " (alength xs) " and " (alength ys))
                    {:xs-count (alength xs) :ys-count (alength ys)}))
  (check-axis! method :xs xs opts))

(defn solve-LU-or-SVD
  "Solves the linear system `A x = b` with LU decomposition, falls back to SVD (minimum norm least squares) when the matrix is singular.

  Parameters:

  - `A` (`RealMatrix`): square matrix
  - `b` (`RealVector`): right hand side

  Returns the solution as a vector of doubles."
  [^RealMatrix A ^RealVector b]
  (try
    (-> A LUDecomposition. .getSolver (.solve b) v/vec->array seq vec)
    (catch SingularMatrixException _
      (-> A SingularValueDecomposition. .getSolver (.solve b) v/vec->array seq vec))))

;; rbf / kriging

(defn pi-matrix-ys
  "Adds polynomial terms"
  [xss Phi ys polynomial-terms]
  (let [pterms (map polynomial-terms xss)
        P (mat/rows->RealMatrix pterms)
        pos (count xss)
        ptsize (count (first pterms))
        size (+ ptsize pos)
        ^RealMatrix target (MatrixUtils/createRealMatrix size size)]
    (.setSubMatrix target (mat/mat->array2d Phi) 0 0)
    (.setSubMatrix target (mat/mat->array2d P) 0 pos)
    (.setSubMatrix target (mat/mat->array2d (mat/transpose P)) pos 0)
    [target (concat ys (repeat ptsize 0.0))]))
