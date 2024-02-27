(ns fastmath.calculus.finite
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec3]
           [fastmath.java Array]
           [org.apache.commons.math3.linear Array2DRowRealMatrix LUDecomposition ArrayRealVector]))



;; https://github.com/JuliaMath/Richardson.jl/blob/master/src/Richardson.jl
(defn extrapolate
  ([g] (extrapolate g nil))
  ([g {:keys [^double contract ^double power ^double init-h
              ^double rel ^double abs ^long max-evals ^double tol]
       :or {contract 0.5 power 1.0 init-h 0.5
            abs m/MACHINE-EPSILON rel (m/sqrt (m/ulp init-h))
            max-evals Integer/MAX_VALUE tol 2.0}}]
   (let [invcontract (m/pow (m// contract) power)
         nf (fn ^double [^double x0 ^double h] (cond
                                                (m/valid-double? x0) (g x0 h)
                                                (m/pos-inf? x0) (g 0.0 (m// 1.0 h))
                                                (m/neg-inf? x0) (g 0.0 (m// -1.0 h))
                                                :else ##NaN))]
     (fn ^double [^double x0]
       (loop [i (long 1)
              h init-h
              ;; keep in 3-tuple of doubles to avoid destructuring
              ^Vec3 f0+err+invc (Vec3. (nf x0 init-h) ##Inf invcontract) ;; curr val, err, contraction den
              neville (list f0+err+invc)]
         (if (m/>= i max-evals)
           (.x f0+err+invc)
           (let [nh (m/* h contract)
                 nf0 (nf x0 nh)
                 res (reductions (fn [^Vec3 curr ^Vec3 ni+err+invc]
                                   (let [ni (.x ni+err+invc)
                                         ni+1 (.x curr)
                                         c (.z curr)
                                         n (m/+ ni+1 (m// (m/- ni+1 ni) (m/dec c)))]
                                     (Vec3. n (m/abs (m/- n ni)) (m/* c invcontract))))
                                 (Vec3. nf0 ##Inf invcontract) neville)
                 ^Vec3 minimal (reduce (fn [^Vec3 a ^Vec3 b]
                                         (if (m/< (.y a) (.y b)) a b)) res)
                 minimal' (.y minimal)
                 err (.y f0+err+invc)
                 ^Vec3 minimal (if (or (m/invalid-double? minimal')
                                       (m/< err minimal')) f0+err+invc minimal)]
             (if (or (m/invalid-double? minimal')
                     (m/> minimal' (m/* tol err))
                     (m/<= err (m/max (m/* rel (m/abs (.x minimal))) abs)))
               (.x minimal)
               (recur (m/inc i)
                      nh
                      minimal
                      res)))))))))

(defn fd-coeffs
  [^long n offsets]
  (let [noffsets (count offsets)
        matrix (->> (range 1 noffsets)
                    (reduce (fn [buff ^long i]
                              (conj buff (map (fn [^long j]
                                                (m/fpow j i)) offsets)))
                            [(repeat noffsets 1)])
                    (m/seq->double-double-array)
                    (Array2DRowRealMatrix.))
        rhs (->> (range noffsets)
                 (map (fn [^long id]
                        (if (m/== id n) (m/factorial n) 0)))
                 (m/seq->double-array)
                 (ArrayRealVector.))]
    [offsets (-> (LUDecomposition. matrix)
                 (.getSolver)
                 ^ArrayRealVector (.solve rhs)
                 (.getDataRef)
                 (m/double-array->seq))]))

(defn fd-coeffs-for-accuracy
  ([^long n] (fd-coeffs-for-accuracy n 2))
  ([^long n ^long acc] (fd-coeffs-for-accuracy n acc :central))
  ([^long n ^long acc kind]
   (let [num-central (-> (m/inc n)
                         (m// 2.0)
                         (m/floor)
                         (m/* 2.0)
                         (m/dec)
                         (m/+ acc)
                         (unchecked-long))
         num-side (m/quot num-central 2)
         num-coeffs (if (m/even? n) (m/inc num-central) num-central)]
     (fd-coeffs n (case kind
                    :forward (range num-coeffs)
                    :backward (range (m/inc (m/- num-coeffs)) 1)
                    (range (m/- num-side) (m/inc num-side)))))))

(defmacro produce-symbols
  [n acc kind]
  (let [[offsets coeffs] (fd-coeffs-for-accuracy n acc kind)
        x (with-meta (symbol "x") {:tag 'double})
        h (with-meta (symbol "h") {:tag 'double})
        hn (with-meta (symbol "hn") {:tag 'double})
        fname (gensym (str "local-" n "-" acc "-" (name kind)))]
    `(fn ~fname
       ([~x] (~fname ~x ~h ~hn))
       ([~x ~h] (~fname ~x ~h (m/fpow ~h ~n)))
       ([~x ~h ~hn]        
        (m// (m/+ ~@(map (fn [id coeff]
                           `(m/* ~coeff (double (~'f (m/+ ~x (m/* ~id ~h)))))) offsets coeffs)) ~hn)))))

(defmacro produce-case-for-methods
  [n acc]
  `(case ~'method
     :forward (produce-symbols ~n ~acc :forward)
     :backward (produce-symbols ~n ~acc :backward)
     (produce-symbols ~n ~acc :central)))

(defmacro produce-cases-for-accuracy
  [n]
  (let [xx (with-meta (symbol "xx") {:tag 'double})
        coeff (with-meta (symbol "coeff") {:tag 'double})
        id (with-meta (symbol "id") {:tag 'long})
        h (with-meta (symbol "h") {:tag 'double})
        hn (with-meta (symbol "hn") {:tag 'double})
        fname (gensym (str "local-" n))]
    `(case (unchecked-int ~'acc)
       1 (produce-case-for-methods ~n 1)
       2 (produce-case-for-methods ~n 2)
       3 (produce-case-for-methods ~n 3)
       4 (produce-case-for-methods ~n 4)
       5 (produce-case-for-methods ~n 5)
       6 (produce-case-for-methods ~n 6)
       (let [[offsets# coeffs#] (fd-coeffs-for-accuracy ~n ~'acc ~'method)]
         (fn ~fname
           ([~xx] (~fname ~xx ~h ~hn))
           ([~xx ~h] (~fname ~xx ~h (m/fpow ~h ~n)))
           ([~xx ~h ~hn] (m// (v/sum (map (fn [~id  ~coeff]
                                            (m/* ~coeff (double (~'f (m/+ ~xx (m/* ~id ~h)))))) offsets# coeffs#)) ~hn)))))))

(defn derivative
  ([f] (derivative f 1))
  ([f ^long n] (derivative f n nil))
  ([f ^long n {:keys [^int acc ^double h method extrapolate?]
               :or {acc 2 h 0.0 method :central}}]
   (if (m/zero? n)
     f
     (let [h (if (m/zero? h) (m/pow m/MACHINE-EPSILON (m// (m/+ n 2.0))) h)
           hn (m/fpow h n)
           result (case n
                    1 (produce-cases-for-accuracy 1)
                    2 (produce-cases-for-accuracy 2)
                    3 (produce-cases-for-accuracy 3)
                    4 (produce-cases-for-accuracy 4)
                    5 (produce-cases-for-accuracy 5)
                    6 (produce-cases-for-accuracy 6)
                    (let [[offsets coeffs] (fd-coeffs-for-accuracy n acc method)]
                      (fn local-rest
                        ([^double x] (local-rest x h hn))
                        ([^double x ^double h] (local-rest x h (m/fpow h n)))
                        ([^double x ^double h ^double hn] (m// (v/sum (map (fn [^long id ^double coeff]
                                                                             (m/* coeff ^double (f (m/+ x (m/* id h))))) offsets coeffs))
                                                               hn)))))]
       (if extrapolate?
         (let [eoptions (if (map? extrapolate?) extrapolate? {})]
           (extrapolate result (if (:power eoptions)
                                 eoptions 
                                 (assoc eoptions :power (if (= method :central) 2 1)))))
         result)))))

(defn gradient
  ([f] (gradient f nil))
  ([f {:keys [^double h ^int acc]
       :or {h 1.0e-6 acc 2}}]
   (case (unchecked-int acc)
     2 (let [h2 (m// (m/* 2.0 h))]
         (fn local-gradient-2
           ([v] (local-gradient-2 v h h2))
           ([v ^double h] (local-gradient-2 v h (m// (m/* 2.0 h))))
           ([v ^double h ^double h2]
            (let [v (vec v)]
              (map (fn [^long id]
                     (let [^double vid (v id)
                           ^double x1 (f (assoc v id (m/+ vid h)))
                           ^double x2 (f (assoc v id (m/- vid h)))]
                       (m/* (m/- x1 x2) h2))) (range (count v)))))))
     4 (fn local-gradient-4
         ([v] (local-gradient-4 v h))
         ([v ^double h]
          (let [v (vec v)]
            (map (fn [^long id]
                   (let [^double vid (v id)
                         ^double x1 (f (assoc v id (m/- vid h h)))
                         ^double x2 (f (assoc v id (m/- vid h)))
                         ^double x3 (f (assoc v id (m/+ vid h)))
                         ^double x4 (f (assoc v id (m/+ vid h h)))]
                     (m// (m/+ (m/* 0.08333333333333333 x1)
                               (m/* -0.6666666666666666 x2)
                               (m/* 0.6666666666666666 x3)
                               (m/* -0.08333333333333333 x4)) h))) (range (count v)))))))))

(defn hessian
  "Creates function returning Hessian matrix for mulitvariate function `f` and given `:h` step (default: `5.0e-3`)."
  ([f] (hessian f nil))
  ([f {:keys [^double h ]
       :or {h 5.0e-3}}]
   (let [h4 (m// (m/* 4.0 h h))
         h2 (m// (m/* h h))]
     (fn local-hessian
       ([v] (local-hessian v h h2 h4))
       ([v ^double h] (let [h2 (m// (m/* h h))
                            h4 (m// h2 4.0)]
                        (local-hessian v h h2 h4)))
       ([v ^double h ^double h2 ^double h4]
        (let [v (vec v)
              cv (count v)
              r (range cv)
              buff (double-array (m/* cv cv))
              fv-2 (m/* -2.0 ^double (f v))]
          (doseq [^long i r
                  ^long j r
                  :when (m/<= i j)]
            (if (m/== i j)
              (let [^double vi (v i)
                    ^double x1 (f (assoc v i (m/+ vi h)))
                    ^double x2 (f (assoc v i (m/- vi h)))]
                (Array/set2d buff cv i i (m/* (m/+ x1 fv-2 x2) h2)))
              (let [^double vi (v i)
                    ^double vj (v j)
                    ^double x1 (f (assoc v i (m/+ vi h) j (m/+ vj h)))
                    ^double x2 (f (assoc v i (m/+ vi h) j (m/- vj h)))
                    ^double x3 (f (assoc v i (m/- vi h) j (m/+ vj h)))
                    ^double x4 (f (assoc v i (m/- vi h) j (m/- vj h)))
                    val (m/* (m/- (m/+ x1 x4) x2 x3) h4)]
                (Array/set2d buff cv i j val)
                (Array/set2d buff cv j i val))))
          (partition cv buff)))))))

