(ns fastmath.complex-test
  (:require [fastmath.complex :as c]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is are testing]]))

;; expected reults verified with Wolfram Alpha and R

;; R code
;;
;; p<-c(-2.2, -0.8, 0.0, 0.8, 2.2)
;; d<-expand.grid(a=p,b=p)
;; f <- function(x,a) { atan(complex(real=x[1],imaginary=x[2])) }
;; res<-apply(d,1,f); print(cbind(d,Re(res), Im(res)), row.names=F)
;;

(deftest complex-operations
  (is (= (v/vec2 1.0 1.0) (c/complex 1.0 1.0)))
  (is (= (c/complex 1.0 1.0) (c/add c/I c/ONE)))
  (is (= (m/sqrt 2.0) (c/abs (c/complex 1.0 1.0))))
  (is (= (c/complex 1.0 1.0) (c/sub c/ONE c/I-)))
  (is (= -90.0 (m/degrees (c/arg c/I-))))
  (is (= (c/complex 0.0 -1.0) (c/conjugate c/I)))
  (is (= (c/complex 0.44 0.08) (c/div (c/complex 1 2) (c/complex 3 4))))
  (is (= (c/complex 0.12 -0.16) (c/reciprocal (c/complex 3 4))))
  (is (= (c/complex -5.0 10.0) (c/mult (c/complex 1 2) (c/complex 3 4))))
  (is (= (c/complex -1.0 2.0) (c/neg (c/complex 1.0 -2.0))))
  (is (= (c/complex -3.0 4.0) (c/sq (c/complex 1.0 2.0))))
  (is (= (c/complex 25.0 0.0) (c/sq (c/complex 5.0 0.0))))
  (is (= (c/complex (m/sqrt 2.0) 0.0) (c/sqrt (c/complex 2.0 0.0))))
  (is (c/delta-eq (c/complex 1.0 2.0) (c/sqrt (c/complex -3.0 4.0))))
  (is (c/delta-eq (c/complex 4.08 -2.94) (c/sqrt1z (c/complex 3 4)) 1.0e-3)))

(defn- rel-eq?
  "Relative equality of two doubles (infinities and zeros must be equal)."
  [^double expected ^double actual ^double tol]
  (or (== expected actual)
      (<= (m/abs (- expected actual)) (* tol (m/abs expected)))))

(defn- c-rel-eq? [[er ei] ^fastmath.vector.Vec2 z tol]
  (and (rel-eq? er (.x z) tol) (rel-eq? ei (.y z) tol)))

;; reference values: mpmath (60 digits)
(deftest robust-ranges
  (testing "abs does not overflow or underflow"
    (is (rel-eq? 5e200 (c/abs (c/complex 3e200 4e200)) 1e-15))
    (is (rel-eq? 5e-200 (c/abs (c/complex 3e-200 4e-200)) 1e-15))
    (is (rel-eq? 5e-310 (c/abs (c/complex 3e-310 4e-310)) 1e-3) "subnormal: the inputs have few significant bits")
    (is (= 0.0 (c/abs c/ZERO)))
    (is (= ##Inf (c/abs (c/complex ##Inf ##NaN))) "an infinite part wins over NaN")
    (is (= ##Inf (c/abs (c/complex 1.0 ##-Inf))))
    (is (m/nan? (c/abs (c/complex ##NaN 1.0)))))
  (testing "division is exact for representable quotients far from the range limits"
    (is (c-rel-eq? [0.44 0.08] (c/div (c/complex 1e200 2e200) (c/complex 3e200 4e200)) 1e-15))
    (is (c-rel-eq? [-0.2 0.4] (c/div (c/complex 1e-200 2e-200) (c/complex 3e-200 -4e-200)) 1e-15))
    (is (c-rel-eq? [0.0 1.0] (c/div (c/complex 1e200 1e200) (c/complex 1e200 -1e200)) 1e-15) "(1+i)/(1-i) = i at a huge scale")
    (is (c-rel-eq? [0.0 1.0] (c/div (c/complex 1e-200 1e-200) (c/complex 1e-200 -1e-200)) 1e-15) "and at a tiny scale")
    (is (c-rel-eq? [1.0 0.0] (c/div (c/complex 1e300 -2e300) (c/complex 1e300 -2e300)) 1e-15) "z / z")
    (is (= (c/complex 0.44 0.08) (c/div (c/complex 1 2) (c/complex 3 4))) "the exact case of the first tests")
    (is (c-rel-eq? [2.0 0.0] (c/div (c/complex 6.0 0.0) (c/complex 3.0 0.0)) 0.0) "real over real")
    (is (c-rel-eq? [0.0 -2.0] (c/div (c/complex 6.0 0.0) (c/complex 0.0 3.0)) 1e-15) "real over imaginary")
    (is (c/nan? (c/div (c/complex 1.0 1.0) c/ZERO)) "division by zero gives NaN parts")
    (is (= [##Inf 0.0] ((juxt #(.x ^fastmath.vector.Vec2 %) #(.y ^fastmath.vector.Vec2 %)) (c/div (c/complex 1e300 1e300) (c/complex 1e-300 1e-300)))) "an unrepresentable quotient overflows to infinity"))
  (testing "reciprocal"
    (is (c-rel-eq? [5e-201 -5e-201] (c/reciprocal (c/complex 1e200 1e200)) 1e-15))
    (is (c-rel-eq? [5e199 -5e199] (c/reciprocal (c/complex 1e-200 1e-200)) 1e-15))
    (is (= (c/complex 0.12 -0.16) (c/reciprocal (c/complex 3 4))))
    (is (= (c/complex 0.0 -0.25) (c/reciprocal (c/complex 0.0 4.0))) "pure imaginary")
    (is (= (c/complex 0.25 0.0) (c/reciprocal (c/complex 4.0 0.0))) "real")
    (is (= (c/complex ##Inf ##-Inf) (c/reciprocal c/ZERO)) "kept: zero gives infinity")
    (is (= (c/complex 0.0 -0.0) (c/reciprocal (c/complex ##Inf 0.0))) "kept: infinity gives zero")
    (is (c/nan? (c/reciprocal (c/complex ##NaN 1.0)))))
  (testing "sqrt is accurate over the whole range"
    (is (c-rel-eq? [1.09868411346781e+100 4.550898605622274e+99] (c/sqrt (c/complex 1e200 1e200)) 1e-15))
    (is (c-rel-eq? [1.09868411346781e-100 4.550898605622273e-101] (c/sqrt (c/complex 1e-200 1e-200)) 1e-15))
    (is (c-rel-eq? [1e150 2e150] (c/sqrt (c/complex -3e300 4e300)) 1e-15) "|z|^2 would overflow")
    (is (c-rel-eq? [5e-26 1e-5] (c/sqrt (c/complex -1e-10 1e-30)) 1e-15) "no cancellation near the negative real axis")
    (is (c-rel-eq? [1e-5 5e-26] (c/sqrt (c/complex 1e-10 1e-30)) 1e-15) "no cancellation near the positive real axis")
    (is (c-rel-eq? [1.4325088230154573e154 5.933645827121221e153] (c/sqrt (c/complex 1.7e308 1.7e308)) 1e-14) "near the largest double: |z| itself overflows, reference mpmath")
    (is (= (c/complex 0.0 2.0) (c/sqrt (c/complex -4.0 0.0))) "the negative real axis")
    (is (= (c/complex 0.0 -2.0) (c/sqrt (c/complex -4.0 -0.0))) "the other side of the branch cut")
    (is (= (c/complex 2.0 0.0) (c/sqrt (c/complex 4.0 0.0))))
    (is (= (c/complex 0.0 0.0) (c/sqrt c/ZERO)))
    (is (= (c/complex 1.0 2.0) (c/sqrt (c/complex -3.0 4.0))) "exact square")
    (is (= (c/complex 1.0 -2.0) (c/sqrt (c/complex -3.0 -4.0))) "conjugate symmetry")
    (is (= ##Inf (.x (c/sqrt (c/complex ##Inf 1.0)))))
    (is (c/nan? (c/sqrt (c/complex ##NaN 1.0))))
    (doseq [z [(c/complex 3.0 -4.0) (c/complex -0.3 7.0) (c/complex 1e-5 -1e5) (c/complex -1e150 1e-150)]
            :let [r (c/sqrt z)]]
      (is (>= (.x r) 0.0) "the real part is not negative")
      ;; the square loses the small part by cancellation: compared to the size of z
      (is (<= (c/abs (c/sub z (c/sq r))) (* 1e-14 (c/abs z))) "the square of the root"))))

(deftest complex-trig
  (is (c/delta-eq (c/complex 2.03 -3.05) (c/cos  (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 3.17 1.96)  (c/sin  (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex -0.64 1.07) (c/cosh (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex -0.49 1.4)  (c/sinh (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 0.03 1.01)  (c/tan  (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 1.17 -0.24) (c/tanh (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 0.15 0.23)  (c/sec  (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 0.23 -0.14) (c/csc  (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 1.14 -1.53) (c/acos (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 0.43 1.53)  (c/asin (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex 1.34 0.4)   (c/atan (c/complex 1 2)) 1.0e-2)))

(deftest log-exp
  (is (c/delta-eq c/ZERO (c/add (c/exp (c/complex 0 m/PI)) c/ONE)))
  (is (c/delta-eq (c/complex -1.13 2.47) (c/exp (c/complex 1 2)) 1.0e-2))
  (is (= c/ONE (c/log (c/complex m/E 0.0))))
  (is (c/delta-eq (c/complex 0.8 1.11) (c/log (c/complex 1 2)) 1.0e-2))
  (is (c/delta-eq (c/complex -0.22 0.1) (c/pow (c/complex 1 2) (c/complex 1 2)) 1.0e-2))
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/log (c/complex a b)))
    -2.2 -2.2 1.1350310 -2.3561945
    0.0 -2.2 0.7884574 -1.5707963
    2.2 -2.2 1.1350310 -0.7853982
    -2.2  0.0 0.7884574  3.1415927
    2.2  0.0 0.7884574  0.0000000
    -2.2  2.2 1.1350310  2.3561945
    0.0  2.2 0.7884574  1.5707963
    2.2  2.2 1.1350310  0.7853982))

(deftest power-branches
  (is (c/nan? (c/pow c/ZERO c/ZERO)))
  (is (c/zero? (c/pow c/ZERO (c/complex 1 1))))
  ;; (-0)^1 = -0; the earlier expectation (Inf, Inf) encoded a bug (audit CX-3)
  (is (c/zero? (c/pow (c/complex -0.0 0.0) c/ONE)))
  (let [w (c/pow c/ZERO (c/complex -1.0 0.0))]
    (is (= ##Inf (c/re w)))
    (is (m/nan? (c/im w))))
  (is (= c/ONE (c/pow (c/complex 1 1) c/ZERO)))

  ;; real powers
  (is (c/delta-eq (c/complex 4.0 0.0) (c/pow c/TWO c/TWO)))
  (are [zr ra rb] (c/delta-eq (c/complex ra rb)
                              (c/pow (c/complex -2.0 0.0) (c/complex zr 0.0)) 1.0e-4)
    -3.3 -0.0597 0.0821
    -2.2 0.1760 -0.1279
    -1.1 -0.4437 0.1442
    -0.5 0.0 -0.7071
    0.5 0.0 1.4142
    1.1 -2.0386 -0.6624
    2.2 3.7173 2.7008
    3.3 -5.7892 -7.9681)

  (are [r i zi ra rb] (c/delta-eq (c/complex ra rb)
                                  (c/pow (c/complex r i) (c/complex 0.0 zi)) 1.0e-4)
    2.0 0.0 -3.3 -0.6568167 -0.7540503
    2.0 0.0 -2.2 0.0459 -0.9989
    2.0 0.0 -1.1 0.7231378 -0.6907038
    2.0 0.0 -0.5 0.9405421 -0.3396771
    2.0 0.0 0.5  0.9405421 0.3396771
    2.0 0.0 1.1 0.7231378 0.6907038
    2.0 0.0 2.2 0.0458564 0.998948
    2.0 0.0 3.3 -0.6568167 0.7540503

    -2.0 0.0 -3.3  -20887.4855 -23979.6184
    -2.0 0.0 -2.2  46.0287 -1002.6997
    -2.0 0.0 -1.1 22.91052 -21.88295
    -2.0 0.0 -0.5 4.524457 -1.634009
    -2.0 0.0 0.5 0.1955195 0.0706119
    -2.0 0.0 1.1 0.0228248 0.02180107
    -2.0 0.0 2.2 4.56849e-05 9.952104e-04

    -2.0 -0.0 3.3 -20887.4855 23979.6184
    -2.0 -0.0 2.2 46.0287 1002.6997
    -2.0 -0.0 1.1 22.91052 21.88295
    -2.0 -0.0 0.5 4.524457 1.634009
    -2.0 -0.0 -0.5 0.1955195 -0.0706119
    -2.0 -0.0 -1.1 0.0228248 -0.02180107
    -2.0 -0.0 -2.2 4.56849e-05 -9.952104e-04)

  (are [r i zi ra rb] (c/delta-eq (c/complex ra rb)
                                  (c/pow (c/complex r i) (c/complex 1.0 zi)) 1.0e-4)

    2.0 0.0 -3.3  -1.313633 -1.508101
    2.0 0.0 -2.2  0.091713 -1.997896
    2.0 0.0 -1.1 1.446276 -1.381408
    2.0 0.0 -0.5 1.881084 -0.679354
    2.0 0.0 0.5  1.881084 0.679354
    2.0 0.0 1.1 1.446276 1.381408
    2.0 0.0 2.2  0.091713 1.997896
    2.0 0.0 3.3  -1.313633 1.508101

    -2.0 0.0 -2.2 -92.0573 2005.3994
    -2.0 0.0 -1.1  -45.82105 43.7659
    -2.0 0.0 -0.5 -9.048913 3.268018
    -2.0 0.0 0.5 -0.391039 -0.1412239
    -2.0 0.0 1.1 -0.04564961 -0.04360215
    -2.0 0.0 2.2  -9.137e-05 -1.990421e-03

    -2.0 -0.0 2.2 -92.0573 -2005.3994
    -2.0 -0.0 1.1  -45.82105 -43.7659
    -2.0 -0.0 0.5 -9.048913 -3.268018
    -2.0 -0.0 -0.5 -0.391039 0.1412239
    -2.0 -0.0 -1.1 -0.04564961 0.04360215
    -2.0 -0.0 -2.2 -9.137e-05 1.990421e-03)

  (are [i zr ra rb] (c/delta-eq (c/complex ra rb)
                                (c/pow (c/complex 0.0 i) (c/complex zr 0.0)) 1.0e-4)
    -2.2 -2.2 -0.1678326 -0.0545321
    -2.2 2.2 -5.389348 1.751105
    2.2 2.2 -5.389348 -1.751105
    2.2 -2.2 -0.1678326 0.0545321)
  
  (are [i zi ra rb] (c/delta-eq (c/complex ra rb)
                                (c/pow (c/complex 0.0 i) (c/complex 0.0 zi)) 1.0e-4)
    -2.2 -2.2  -0.00514733 -0.03114103
    -2.2 2.2 -5.16666 31.25798
    2.2 2.2 -0.00514733 0.03114103
    2.2 -2.2  -5.16666 -31.25798)
  
  (are [i zr zi ra rb] (c/delta-eq (c/complex ra rb)
                                   (c/pow (c/complex 0.0 i) (c/complex zr zi)) 1.0e-4)
    -2.2 2.2 -2.2  0.082272 0.1588163
    -2.2 2.2 2.2 -26.8911 -177.5075
    2.2 2.2 2.2  0.082272 -0.1588163
    2.2 2.2 -2.2  -26.8911 177.5075
    -2.2 -2.2 -2.2  -0.000834296 0.005507175
    -2.2 -2.2 2.2 2.571699 -4.96436
    2.2 -2.2 2.2  -0.000834296 -0.005507175
    2.2 -2.2 -2.2  2.571699 4.96436)
  
  (are [r i zr ra rb] (c/delta-eq (c/complex ra rb)
                                  (c/pow (c/complex r i) (c/complex zr 0.0)) 1.0e-4)
    -2.2 2.2 -2.2 0.03737523 0.07335302 
    -2.2 2.2 2.2   5.51454 -10.8229
    2.2 2.2 2.2  -1.90018 11.99728 
    2.2 2.2 -2.2 -0.01287863 -0.08131245 
    -2.2 -2.2 -2.2  0.03737523 -0.07335302
    -2.2 -2.2 2.2  5.51454 10.8229
    2.2 -2.2 2.2  -1.90018 -11.99728
    2.2 -2.2 -2.2  -0.01287863 0.08131245)
  
  (are [r i zi ra rb] (c/delta-eq (c/complex ra rb)
                                  (c/pow (c/complex r i) (c/complex 0.0 zi)) 1.0e-4)
    -2.2 2.2 -2.2  -142.5533 -107.1431 
    -2.2 2.2 2.2 -0.004482655 0.003369165 
    2.2 2.2 2.2  -0.1420199 0.1067422
    2.2 2.2 -2.2  -4.49949 -3.381818
    -2.2 -2.2 -2.2 -0.004482655 -0.003369165
    -2.2 -2.2 2.2  -142.5533 107.1431
    2.2 -2.2 2.2  -4.49949 3.381818
    2.2 -2.2 -2.2  -0.1420199 -0.1067422)
  
  (are [r i ri zi ra rb] (c/delta-eq (c/complex ra rb)
                                     (c/pow (c/complex r i) (c/complex ri zi)) 1.0e-4)
    2.2 2.2 2.2 2.2     -1.010753 -1.906683          
    2.2 2.2 2.2 -2.2    49.12247 -47.55557           
    2.2 2.2 -2.2 2.2    0.01050849 0.01017329        
    2.2 2.2 -2.2 -2.2   -0.2170367 0.4094177         
    2.2 -2.2 2.2 2.2    49.12247 47.55557            
    2.2 -2.2 2.2 -2.2   -1.010753 1.906683           
    2.2 -2.2 -2.2 2.2   -0.2170367 -0.4094177        
    2.2 -2.2 -2.2 -2.2  0.01050849 -0.01017329       
    -2.2 2.2 2.2 2.2    0.01174435 0.06709475        
    -2.2 2.2 2.2 -2.2   -1945.716 951.995            
    -2.2 2.2 -2.2 2.2   -4.146787e-04 -2.028929e-04  
    -2.2 2.2 -2.2 -2.2  2.53131 -14.46121            
    -2.2 -2.2 2.2 2.2   -1945.716 -951.995           
    -2.2 -2.2 2.2 -2.2  0.01174435 -0.06709475       
    -2.2 -2.2 -2.2 2.2  2.53131 14.46121             
    -2.2 -2.2 -2.2 -2.2 -4.146787e-04 2.028929e-04))

(deftest trigonometric
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/sin (c/complex a b)))
    -2.2 -2.2 -3.6931375 2.623011
    0.0 -2.2  0.0000000 -4.457105
    2.2 -2.2  3.6931375 2.623011
    -2.2  0.0 -0.8084964 0.000000
    0.0  0.0  0.0000000 0.000000
    2.2  0.0  0.8084964 0.000000
    -2.2  2.2 -3.6931375 -2.623011
    0.0  2.2  0.0000000 4.457105
    2.2  2.2  3.6931375 -2.623011)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/cos (c/complex a b)))
    -2.2 -2.2 -2.6882192 -3.603554
    0.0 -2.2  4.5679083  0.000000
    2.2 -2.2 -2.6882192  3.603554
    -2.2  0.0 -0.5885011  0.000000
    0.0  0.0  1.0000000  0.000000
    2.2  0.0 -0.5885011  0.000000
    -2.2  2.2 -2.6882192  3.603554
    0.0  2.2  4.5679083  0.000000
    2.2  2.2 -2.6882192 -3.603554)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/tan (c/complex a b)))
    -2.2 -2.2  0.02354038 -1.0072990
    0.0 -2.2  0.00000000 -0.9757431
    2.2 -2.2 -0.02354038 -1.0072990
    -2.2  0.0  1.37382306  0.0000000
    0.0  0.0  0.00000000  0.0000000
    2.2  0.0 -1.37382306  0.0000000
    -2.2  2.2  0.02354038  1.0072990
    0.0  2.2  0.00000000  0.9757431
    2.2  2.2 -0.02354038 1.0072990)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/sec (c/complex a b)))
    -2.2 -2.2 -0.1330004  0.1782868
    0.0 -2.2  0.2189186  0.0000000
    2.2 -2.2 -0.1330004 -0.1782868
    -2.2  0.0 -1.6992321  0.0000000
    0.0  0.0  1.0000000  0.0000000
    2.2  0.0 -1.6992321  0.0000000
    -2.2  2.2 -0.1330004 -0.1782868
    0.0  2.2  0.2189186  0.0000000
    2.2  2.2 -0.1330004  0.1782868)
  (is (= (c/complex ##Inf -0.0) (c/csc c/ZERO)))
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/csc (c/complex a b)))
    -2.2 -2.2 -0.1799823 -0.1278305
    0.0 -2.2  0.0000000  0.2243609
    2.2 -2.2  0.1799823 -0.1278305
    -2.2  0.0 -1.2368639  0.0000000
    2.2  0.0  1.2368639  0.0000000
    -2.2  2.2 -0.1799823  0.1278305
    0.0  2.2  0.0000000 -0.2243609
    2.2  2.2  0.1799823  0.1278305)
  (is (= (c/complex ##Inf 0.0) (c/cot c/ZERO)))
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/cot (c/complex a b)))
    -2.2 -2.2  0.0231878  0.992212
    0.0 -2.2  0.0000000  1.024860
    2.2 -2.2 -0.0231878  0.992212
    -2.2  0.0  0.7278958  0.000000
    2.2  0.0 -0.7278958  0.000000
    -2.2  2.2  0.0231878 -0.992212
    0.0  2.2  0.0000000 -1.024860
    2.2  2.2 -0.0231878 -0.992212))

(deftest hyperbolic
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/sinh (c/complex a b)))
    -2.2 -2.2  2.623011 -3.6931375
    0.0 -2.2  0.000000 -0.8084964
    2.2 -2.2 -2.623011 -3.6931375
    -2.2  0.0 -4.457105  0.0000000
    0.0  0.0  0.000000  0.0000000
    2.2  0.0  4.457105  0.0000000
    -2.2  2.2  2.623011  3.6931375
    0.0  2.2  0.000000  0.8084964
    2.2  2.2 -2.623011  3.6931375)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/cosh (c/complex a b)))
    -2.2 -2.2 -2.6882192  3.603554
    0.0 -2.2 -0.5885011  0.000000
    2.2 -2.2 -2.6882192 -3.603554
    -2.2  0.0  4.5679083  0.000000
    0.0  0.0  1.0000000  0.000000
    2.2  0.0  4.5679083  0.000000
    -2.2  2.2 -2.6882192 -3.603554
    0.0  2.2 -0.5885011  0.000000
    2.2  2.2 -2.6882192  3.603554)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/tanh (c/complex a b)))
    -2.2 -2.2 -1.0072990  0.02354038
    0.0 -2.2  0.0000000  1.37382306
    2.2 -2.2  1.0072990  0.02354038
    -2.2  0.0 -0.9757431  0.00000000
    0.0  0.0  0.0000000  0.00000000
    2.2  0.0  0.9757431  0.00000000
    -2.2  2.2 -1.0072990 -0.02354038
    0.0  2.2  0.0000000 -1.37382306
    2.2  2.2  1.0072990 -0.02354038)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/sech (c/complex a b)))
    -2.2 -2.2 -0.1330004 -0.1782868
    0.0 -2.2 -1.6992321  0.0000000
    2.2 -2.2 -0.1330004  0.1782868
    -2.2  0.0  0.2189186  0.0000000
    0.0  0.0  1.0000000  0.0000000
    2.2  0.0  0.2189186  0.0000000
    -2.2  2.2 -0.1330004  0.1782868
    0.0  2.2 -1.6992321  0.0000000
    2.2  2.2 -0.1330004 -0.1782868)
  (is (= (c/complex ##Inf -0.0) (c/csch c/ZERO)))
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/csch (c/complex a b)))
    -2.2 -2.2  0.1278305  0.1799823
    0.0 -2.2  0.0000000  1.2368639
    2.2 -2.2 -0.1278305  0.1799823
    -2.2  0.0 -0.2243609  0.0000000
    2.2  0.0  0.2243609  0.0000000
    -2.2  2.2  0.1278305 -0.1799823
    0.0  2.2  0.0000000 -1.2368639
    2.2  2.2 -0.1278305 -0.1799823)
  (is (= (c/complex ##Inf 0.0) (c/coth c/ZERO)))
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/coth (c/complex a b)))
    -2.2 -2.2 -0.992212 -0.0231878
    0.0 -2.2  0.000000 -0.7278958
    2.2 -2.2  0.992212 -0.0231878
    -2.2  0.0 -1.024860  0.0000000
    2.2  0.0  1.024860  0.0000000
    -2.2  2.2 -0.992212  0.0231878
    0.0  2.2  0.000000  0.7278958
    2.2  2.2  0.992212  0.0231878))

;;

(deftest inverse-trigonometric
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/asin (c/complex a b)))
    -2.2 -2.2 -0.7596289 -1.8291748
    -0.8 -2.2 -0.3222810 -1.5779921
    0.0 -2.2  0.0000000 -1.5296605
    0.8 -2.2  0.3222810 -1.5779921
    2.2 -2.2  0.7596289 -1.8291748
    -2.2 -0.8 -1.1893629 -1.5083878
    -0.8 -0.8 -0.6096435 -0.8641474
    0.0 -0.8  0.0000000 -0.7326683
    0.8 -0.8  0.6096435 -0.8641474
    2.2 -0.8  1.1893629 -1.5083878
    -2.2  0.0 -1.5707963  1.4254169
    -0.8  0.0 -0.9272952  0.0000000
    0.0  0.0  0.0000000  0.0000000
    0.8  0.0  0.9272952  0.0000000
    2.2  -0.0  1.5707963 -1.4254169 ;; -0.0
    -2.2  0.8 -1.1893629  1.5083878
    -0.8  0.8 -0.6096435  0.8641474
    0.0  0.8  0.0000000  0.7326683
    0.8  0.8  0.6096435  0.8641474
    2.2  0.8  1.1893629  1.5083878
    -2.2  2.2 -0.7596289  1.8291748
    -0.8  2.2 -0.3222810  1.5779921
    0.0  2.2  0.0000000  1.5296605
    0.8  2.2  0.3222810  1.5779921
    2.2  2.2  0.7596289  1.8291748)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acos (c/complex a b)))
    -2.2 -2.2 2.3304252  1.8291748
    -0.8 -2.2 1.8930773  1.5779921
    0.0 -2.2 1.5707963  1.5296605
    0.8 -2.2 1.2485153  1.5779921
    2.2 -2.2 0.8111675  1.8291748
    -2.2 -0.8 2.7601592  1.5083878
    -0.8 -0.8 2.1804398  0.8641474
    0.0 -0.8 1.5707963  0.7326683
    0.8 -0.8 0.9611528  0.8641474
    2.2 -0.8 0.3814334  1.5083878
    -2.2  0.0 3.1415927 -1.4254169
    -0.8  0.0 2.4980915  0.0000000
    0.0  0.0 1.5707963  0.0000000
    0.8  0.0 0.6435011  0.0000000
    2.2  -0.0 0.0000000  1.4254169 ;; -0.0
    -2.2  0.8 2.7601592 -1.5083878
    -0.8  0.8 2.1804398 -0.8641474
    0.0  0.8 1.5707963 -0.7326683
    0.8  0.8 0.9611528 -0.8641474
    2.2  0.8 0.3814334 -1.5083878
    -2.2  2.2 2.3304252 -1.8291748
    -0.8  2.2 1.8930773 -1.5779921
    0.0  2.2 1.5707963 -1.5296605
    0.8  2.2 1.2485153 -1.5779921
    2.2  2.2 0.8111675 -1.8291748)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/atan (c/complex a b)))
    -2.2 -2.2 -1.3362152 -0.2189998
    -0.8 -2.2 -1.3992844 -0.4136396
    -0.0 -2.2 -1.5707963 -0.4904146
    0.8 -2.2  1.3992844 -0.4136396
    2.2 -2.2  1.3362152 -0.2189998
    -2.2 -0.8 -1.1826016 -0.1260617
    -0.8 -0.8 -0.8720210 -0.4353744
    0.0 -0.8  0.0000000 -1.0986123
    0.8 -0.8  0.8720210 -0.4353744
    2.2 -0.8  1.1826016 -0.1260617
    -2.2  0.0 -1.1441688  0.0000000
    -0.8  0.0 -0.6747409  0.0000000
    0.0  0.0  0.0000000  0.0000000
    0.8  0.0  0.6747409  0.0000000
    2.2  0.0  1.1441688  0.0000000
    -2.2  0.8 -1.1826016  0.1260617
    -0.8  0.8 -0.8720210  0.4353744
    0.0  0.8  0.0000000  1.0986123
    0.8  0.8  0.8720210  0.4353744
    2.2  0.8  1.1826016  0.1260617
    -2.2  2.2 -1.3362152  0.2189998
    -0.8  2.2 -1.3992844  0.4136396
    0.0  2.2  1.5707963  0.4904146
    0.8  2.2  1.3992844  0.4136396
    2.2  2.2  1.3362152  0.2189998)
  (is (= (c/complex 0.0 ##Inf) (c/asec c/ZERO)) "the limit along the real axis, as Octave; was (Inf, Inf)")
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/asec (c/complex a b)))
    -2.2 -2.2 1.793986 -0.2309936
    -0.8 -2.2 1.706513 -0.3948458
    0.0 -2.2 1.570796 -0.4401912
    0.8 -2.2 1.435079 -0.3948458
    2.2 -2.2 1.347607 -0.2309936
    -2.2 -0.8 1.978474 -0.1583552
    -0.8 -0.8 2.101441 -0.6727468
    0.0 -0.8 1.570796 -1.0475930
    0.8 -0.8 1.040151 -0.6727468
    2.2 -0.8 1.163119 -0.1583552
    -2.2  0.0 2.042658  0.0000000
    -0.8  -0.0 3.141593 -0.6931472
    0.8  0.0 0.000000  0.6931472
    2.2  0.0 1.098934  0.0000000
    -2.2  0.8 1.978474  0.1583552
    -0.8  0.8 2.101441  0.6727468
    0.0  0.8 1.570796  1.0475930
    0.8  0.8 1.040151  0.6727468
    2.2  0.8 1.163119  0.1583552
    -2.2  2.2 1.793986  0.2309936
    -0.8  2.2 1.706513  0.3948458
    0.0  2.2 1.570796  0.4401912
    0.8  2.2 1.435079  0.3948458
    2.2  2.2 1.347607  0.2309936)
  (is (= (c/complex m/HALF_PI ##-Inf) (c/acsc c/ZERO)) "the limit along the real axis, as Octave; was (Inf, Inf)")
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acsc (c/complex a b)))
    -2.2 -2.2 -0.2231896  0.2309936
    -0.8 -2.2 -0.1357170  0.3948458
    0.0 -2.2  0.0000000  0.4401912
    0.8 -2.2  0.1357170  0.3948458
    2.2 -2.2  0.2231896  0.2309936
    -2.2 -0.8 -0.4076776  0.1583552
    -0.8 -0.8 -0.5306450  0.6727468
    0.0 -0.8  0.0000000  1.0475930
    0.8 -0.8  0.5306450  0.6727468
    2.2 -0.8  0.4076776  0.1583552
    -2.2  0.0 -0.4718618  0.0000000
    -0.8  -0.0 -1.5707963  0.6931472
    0.8  0.0  1.5707963 -0.6931472
    2.2  0.0  0.4718618  0.0000000
    -2.2  0.8 -0.4076776 -0.1583552
    -0.8  0.8 -0.5306450 -0.6727468
    0.0  0.8  0.0000000 -1.0475930
    0.8  0.8  0.5306450 -0.6727468
    2.2  0.8  0.4076776 -0.1583552
    -2.2  2.2 -0.2231896 -0.2309936
    -0.8  2.2 -0.1357170 -0.3948458
    0.0  2.2  0.0000000 -0.4401912
    0.8  2.2  0.1357170 -0.3948458
    2.2  2.2  0.2231896 -0.2309936)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acot (c/complex a b)))
    -2.2 -2.2 -0.2345811  0.2189998
    -0.8 -2.2 -0.1715120  0.4136396
    0.0 -2.2  0.0000000  0.4904146
    0.8 -2.2  0.1715120  0.4136396
    2.2 -2.2  0.2345811  0.2189998
    -2.2 -0.8 -0.3881947  0.1260617
    -0.8 -0.8 -0.6987753  0.4353744
    0.0 -0.8  1.5707963  1.0986123
    0.8 -0.8  0.6987753  0.4353744
    2.2 -0.8  0.3881947  0.1260617
    -2.2  0.0 -0.4266275  0.0000000
    -0.8  0.0 -0.8960554  0.0000000
    0.0  0.0  1.5707963  0.0000000
    0.8  0.0  0.8960554  0.0000000
    2.2  0.0  0.4266275  0.0000000
    -2.2  0.8 -0.3881947 -0.1260617
    -0.8  0.8 -0.6987753 -0.4353744
    -0.0  0.8 -1.5707963 -1.0986123
    0.8  0.8  0.6987753 -0.4353744
    2.2  0.8  0.3881947 -0.1260617
    -2.2  2.2 -0.2345811 -0.2189998
    -0.8  2.2 -0.1715120 -0.4136396
    0.0  2.2  0.0000000 -0.4904146
    0.8  2.2  0.1715120 -0.4136396
    2.2  2.2  0.2345811 -0.2189998))

(deftest inverse-hyperbolic
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/asinh (c/complex a b)))
    -2.2 -2.2 -1.8291748 -0.7596289
    -0.8 -2.2 -1.5083878 -1.1893629
    -0.0 -2.2 -1.4254169 -1.5707963
    0.8 -2.2  1.5083878 -1.1893629
    2.2 -2.2  1.8291748 -0.7596289
    -2.2 -0.8 -1.5779921 -0.3222810
    -0.8 -0.8 -0.8641474 -0.6096435
    0.0 -0.8  0.0000000 -0.9272952
    0.8 -0.8  0.8641474 -0.6096435
    2.2 -0.8  1.5779921 -0.3222810
    -2.2  0.0 -1.5296605  0.0000000
    -0.8  0.0 -0.7326683  0.0000000
    0.0  0.0  0.0000000  0.0000000
    0.8  0.0  0.7326683  0.0000000
    2.2  0.0  1.5296605  0.0000000
    -2.2  0.8 -1.5779921  0.3222810
    -0.8  0.8 -0.8641474  0.6096435
    0.0  0.8  0.0000000  0.9272952
    0.8  0.8  0.8641474  0.6096435
    2.2  0.8  1.5779921  0.3222810
    -2.2  2.2 -1.8291748  0.7596289
    -0.8  2.2 -1.5083878  1.1893629
    0.0  2.2  1.4254169  1.5707963
    0.8  2.2  1.5083878  1.1893629
    2.2  2.2  1.8291748  0.7596289)
  ;; reference values: Julia ComplexF64 (C99 branch cuts), checked against mpmath off the real axis
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acosh (c/complex a b)))
    -2.2 -2.2 1.8291748 -2.3304252
    -0.8 -2.2 1.5779921 -1.8930773
    0.0 -2.2 1.5296605 -1.5707963
    0.8 -2.2 1.5779921 -1.2485153
    2.2 -2.2 1.8291748 -0.8111675
    -2.2 -0.8 1.5083878 -2.7601592
    -0.8 -0.8 0.8641474 -2.1804398
    0.0 -0.8 0.7326683 -1.5707963
    0.8 -0.8 0.8641474 -0.9611528
    2.2 -0.8 1.5083878 -0.3814334
    -2.2 -0.0 1.4254169 -3.1415927
    -0.8 -0.0 0.0000000 -2.4980915
    0.8 -0.0 0.0000000 -0.6435011
    2.2 -0.0 1.4254169 -0.0000000
    -2.2 0.0 1.4254169 3.1415927
    -0.8 0.0 0.0000000 2.4980915
    0.8 0.0 0.0000000 0.6435011
    2.2 0.0 1.4254169 0.0000000
    -2.2 0.8 1.5083878 2.7601592
    -0.8 0.8 0.8641474 2.1804398
    0.0 0.8 0.7326683 1.5707963
    0.8 0.8 0.8641474 0.9611528
    2.2 0.8 1.5083878 0.3814334
    -2.2 2.2 1.8291748 2.3304252
    -0.8 2.2 1.5779921 1.8930773
    0.0 2.2 1.5296605 1.5707963
    0.8 2.2 1.5779921 1.2485153
    2.2 2.2 1.8291748 0.8111675
    0.0 0.0 0.0000000 1.5707963)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/atanh (c/complex a b)))
    -2.2 -2.2 -0.2189998 -1.3362152
    -0.8 -2.2 -0.1260617 -1.1826016
    0.0 -2.2  0.0000000 -1.1441688
    0.8 -2.2  0.1260617 -1.1826016
    2.2 -2.2  0.2189998 -1.3362152
    -2.2 -0.8 -0.4136396 -1.3992844
    -0.8 -0.8 -0.4353744 -0.8720210
    0.0 -0.8  0.0000000 -0.6747409
    0.8 -0.8  0.4353744 -0.8720210
    2.2 -0.8  0.4136396 -1.3992844
    -2.2  0.0 -0.4904146  1.5707963
    -0.8  0.0 -1.0986123  0.0000000
    0.0  0.0  0.0000000  0.0000000
    0.8  0.0  1.0986123  0.0000000
    2.2  -0.0  0.4904146 -1.5707963
    -2.2  0.8 -0.4136396  1.3992844
    -0.8  0.8 -0.4353744  0.8720210
    0.0  0.8  0.0000000  0.6747409
    0.8  0.8  0.4353744  0.8720210
    2.2  0.8  0.4136396  1.3992844
    -2.2  2.2 -0.2189998  1.3362152
    -0.8  2.2 -0.1260617  1.1826016
    0.0  2.2  0.0000000  1.1441688
    0.8  2.2  0.1260617  1.1826016
    2.2  2.2  0.2189998  1.3362152)
  ;; reference values: Julia ComplexF64 (C99 branch cuts), checked against mpmath off the real axis
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/asech (c/complex a b)))
    -2.2 -2.2 0.2309936 1.7939860
    -0.8 -2.2 0.3948458 1.7065133
    0.0 -2.2 0.4401912 1.5707963
    0.8 -2.2 0.3948458 1.4350793
    2.2 -2.2 0.2309936 1.3476067
    -2.2 -0.8 0.1583552 1.9784739
    -0.8 -0.8 0.6727468 2.1014414
    0.0 -0.8 1.0475930 1.5707963
    0.8 -0.8 0.6727468 1.0401513
    2.2 -0.8 0.1583552 1.1631187
    -2.2 -0.0 0.0000000 2.0426582
    -0.8 -0.0 0.6931472 3.1415927
    0.8 -0.0 0.6931472 0.0000000
    2.2 -0.0 0.0000000 1.0989345
    -2.2 0.0 0.0000000 -2.0426582
    -0.8 0.0 0.6931472 -3.1415927
    0.8 0.0 0.6931472 -0.0000000
    2.2 0.0 0.0000000 -1.0989345
    -2.2 0.8 0.1583552 -1.9784739
    -0.8 0.8 0.6727468 -2.1014414
    0.0 0.8 1.0475930 -1.5707963
    0.8 0.8 0.6727468 -1.0401513
    2.2 0.8 0.1583552 -1.1631187
    -2.2 2.2 0.2309936 -1.7939860
    -0.8 2.2 0.3948458 -1.7065133
    0.0 2.2 0.4401912 -1.5707963
    0.8 2.2 0.3948458 -1.4350793
    2.2 2.2 0.2309936 -1.3476067)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acsch (c/complex a b)))
    -2.2 -2.2 -0.2309936  0.2231896
    -0.8 -2.2 -0.1583552  0.4076776
    0.0 -2.2  0.0000000  0.4718618
    0.8 -2.2  0.1583552  0.4076776
    2.2 -2.2  0.2309936  0.2231896
    -2.2 -0.8 -0.3948458  0.1357170
    -0.8 -0.8 -0.6727468  0.5306450
    0.0 -0.8  0.6931472  1.5707963
    0.8 -0.8  0.6727468  0.5306450
    2.2 -0.8  0.3948458  0.1357170
    -2.2  0.0 -0.4401912  0.0000000
    -0.8  0.0 -1.0475930  0.0000000
    0.8  0.0  1.0475930  0.0000000
    2.2  0.0  0.4401912  0.0000000
    -2.2  0.8 -0.3948458 -0.1357170
    -0.8  0.8 -0.6727468 -0.5306450
    -0.0  0.8 -0.6931472 -1.5707963
    0.8  0.8  0.6727468 -0.5306450
    2.2  0.8  0.3948458 -0.1357170
    -2.2  2.2 -0.2309936 -0.2231896
    -0.8  2.2 -0.1583552 -0.4076776
    0.0  2.2  0.0000000 -0.4718618
    0.8  2.2  0.1583552 -0.4076776
    2.2  2.2  0.2309936 -0.2231896)
  (are [a b ra rb] (c/delta-eq (c/complex ra rb)
                               (c/acoth (c/complex a b)))
    -2.2 -2.2 -0.2189998  0.2345811
    -0.8 -2.2 -0.1260617  0.3881947
    0.0 -2.2  0.0000000  0.4266275
    0.8 -2.2  0.1260617  0.3881947
    2.2 -2.2  0.2189998  0.2345811
    -2.2 -0.8 -0.4136396  0.1715120
    -0.8 -0.8 -0.4353744  0.6987753
    0.0 -0.8  0.0000000  0.8960554
    0.8 -0.8  0.4353744  0.6987753
    2.2 -0.8  0.4136396  0.1715120
    -2.2  0.0 -0.4904146  0.0000000
    -0.8  -0.0 -1.0986123  1.5707963
    0.8  0.0  1.0986123 -1.5707963
    2.2  0.0  0.4904146  0.0000000
    -2.2  0.8 -0.4136396 -0.1715120
    -0.8  0.8 -0.4353744 -0.6987753
    0.0  0.8  0.0000000 -0.8960554
    0.8  0.8  0.4353744 -0.6987753
    2.2  0.8  0.4136396 -0.1715120
    -2.2  2.2 -0.2189998 -0.2345811
    -0.8  2.2 -0.1260617 -0.3881947
    0.0  2.2  0.0000000 -0.4266275
    0.8  2.2  0.1260617 -0.3881947
    2.2  2.2  0.2189998 -0.2345811))

;; ---------------------------------------------------------------------------
;; audit 2026-10-27: reference data, branch cuts, ranges, new functions

(def ^:private reference (delay (edn/read-string (slurp (io/resource "complex/complex_reference.edn")))))
(def ^:private cuts-reference (delay (edn/read-string (slurp (io/resource "complex/complex_cuts_reference.edn")))))

(def ^:private functions
  {:sqrt c/sqrt :exp c/exp :log c/log
   :sin c/sin :cos c/cos :tan c/tan :sec c/sec :csc c/csc :cot c/cot
   :sinh c/sinh :cosh c/cosh :tanh c/tanh :sech c/sech :csch c/csch :coth c/coth
   :asin c/asin :acos c/acos :atan c/atan :asec c/asec :acsc c/acsc :acot c/acot
   :asinh c/asinh :acosh c/acosh :atanh c/atanh :asech c/asech :acsch c/acsch :acoth c/acoth
   :log1p c/log1p :expm1 c/expm1 :cbrt c/cbrt})

(defn- same-bits?
  "Equal doubles, the sign of a zero and NaN included."
  [^double x ^double y]
  (== (Double/doubleToLongBits x) (Double/doubleToLongBits y)))

(defn- normwise-error
  "Error of `z` against the reference parts `[er ei]`, relative to the magnitude of the reference.

  Infinite reference parts must be matched exactly (also by sign), NaN gives an infinite error."
  ^double [[^double er ^double ei] ^fastmath.vector.Vec2 z]
  (let [r (.x z)
        i (.y z)]
    (cond
      (or (m/nan? r) (m/nan? i)) ##Inf
      (or (and (m/inf? er) (not= er r))
          (and (m/inf? ei) (not= ei i))
          (and (m/inf? r) (not (m/inf? er)))
          (and (m/inf? i) (not (m/inf? ei)))) ##Inf
      :else (let [fr (if (m/inf? er) 0.0 er)
                  fi (if (m/inf? ei) 0.0 ei)
                  dr (if (m/inf? er) 0.0 (- r er))
                  di (if (m/inf? ei) 0.0 (- i ei))]
              (/ (m/hypot dr di) (m/max (m/hypot fr fi) 1.0e-320))))))

;; maximum normwise errors measured against mpmath on 1424 points per function (magnitudes 1e-300 ... 1e300,
;; both signs, around 1 and around the overflow of cosh); the tolerances are set above them:
;; sqrt 2.2e-16, exp 1.1e-15 (log of the magnitude times the argument), log 1.9e-16, sin cos tan sinh cosh tanh sec csc cot
;; sech csch coth 3.8e-16 ... 5.5e-16, asin acos atan asinh acosh atanh 2.2e-16 ... 2.8e-16, log1p 2.0e-16, expm1 4.3e-16
;; asec, asech: 5.0e-11 at 1e-10 from the branch point 1, acot, acoth: 4.2e-12 at 1e-10 from +-i (points of this grid); the reciprocal 1/z
;; is rounded before the inverse function, close to a branch point the relative error grows as eps/distance (random points:
;; asec, asech 8e-10 at the distance 1.4e-8, acot 3e-7 at 2e-11, acoth 1e-6 at 3e-12; acsc, acsch eps/sqrt(distance), 3e-11);
;; acsc, acsch 1.3e-14 on this grid
(def ^:private tolerances
  {:sqrt 1.0e-15 :exp 3.0e-15 :log 1.0e-15
   :sin 1.0e-15 :cos 1.0e-15 :tan 1.0e-15 :sec 1.0e-15 :csc 1.0e-15 :cot 1.0e-15
   :sinh 1.0e-15 :cosh 1.0e-15 :tanh 1.0e-15 :sech 1.0e-15 :csch 1.0e-15 :coth 1.0e-15
   :asin 1.0e-15 :acos 1.0e-15 :atan 1.0e-15 :asinh 1.0e-15 :acosh 1.0e-15 :atanh 1.0e-15
   :asec 1.0e-9 :acsc 1.0e-13 :acot 1.0e-10 :asech 1.0e-9 :acsch 1.0e-13 :acoth 1.0e-10
   :log1p 1.0e-15 :expm1 1.0e-15 :cbrt 1.0e-14})

(deftest reference-values
  (let [rf @reference
        points (:points rf)]
    (doseq [[f tol] tolerances
            :let [fv (functions f)
                  bad (filter (fn [[[a b] e]] (and e (> (normwise-error e (fv (c/complex a b))) tol)))
                              (map vector points (get rf f)))]]
      (is (empty? bad) (str f ": " (count bad) " points above " tol ", the first " (first bad))))))

(deftest reference-pow
  (let [rf @reference
        bad (filter (fn [[[z w] e]] (and e (> (normwise-error e (c/pow (apply c/complex z) (apply c/complex w))) 1.0e-12)))
                    (map vector (:pow-args rf) (:pow rf)))]
    (is (empty? bad) (str (count bad) " points above 1e-12, the first " (first bad)))))

(deftest reference-nth-root
  (let [rf @reference
        bad (filter (fn [[[a b] e]] (and e (> (normwise-error e (c/nth-root (c/complex a b) 5)) 1.0e-13)))
                    (map vector (:points rf) (:nth-root-5 rf)))]
    (is (empty? bad) (str (count bad) " points above 1e-13, the first " (first bad)))))

(deftest branch-cuts
  ;; Julia ComplexF64, the sign of a zero part selects the side of the cut, values and signs of the zeros are compared
  (doseq [[f rows] @cuts-reference
          :let [fv (functions f)]
          [a b er ei] rows
          :let [z (fv (c/complex a b))
                msg (str f "(" a " " b ") expected [" er " " ei "] got [" (c/re z) " " (c/im z) "]")]]
    (is (< (normwise-error [er ei] z) 1.0e-9) msg)
    (when (== 0.0 er) (is (same-bits? er (c/re z)) (str "sign of the real zero: " msg)))
    (when (== 0.0 ei) (is (same-bits? ei (c/im z)) (str "sign of the imaginary zero: " msg)))))

(deftest abs-norm-proj
  (is (= 25.0 (c/abs-sq (c/complex 3 4))))
  (is (= 25.0 (c/norm (c/complex 3 4))))
  (is (= c/ZERO (c/proj c/ZERO)))
  (is (= (c/complex 1.0 2.0) (c/proj (c/complex 1.0 2.0))))
  (are [z r i] (let [p (c/proj z)]
                 (and (= ##Inf (c/re p)) (same-bits? r (c/im p)) (same-bits? i (c/im p))))
    (c/complex ##Inf 3.0) 0.0 0.0
    (c/complex 3.0 ##-Inf) -0.0 -0.0
    (c/complex ##-Inf ##NaN) 0.0 0.0)
  (is (c/nan? (c/proj (c/complex ##NaN 1.0))) "NaN stays")
  (is (= (c/complex 0.0 1.0) (c/flip (c/complex 1.0 0.0)))))

(deftest cis-expm1-log1p
  (is (= c/ONE (c/cis 0.0)))
  (is (c/delta-eq c/I (c/cis m/HALF_PI) 1.0e-15))
  (is (c/delta-eq (c/complex -1.0 0.0) (c/cis m/PI) 1.0e-15))
  (testing "expm1 and log1p do not lose the small part"
    (is (c-rel-eq? [1.0e-10 1.0e-10] (c/expm1 (c/complex 1.0e-10 1.0e-10)) 1.0e-9))
    (is (c-rel-eq? [1.5e-20 -1.0e-10] (c/log1p (c/complex 1.0e-20 -1.0e-10)) 1.0e-9) "2a + a^2 + b^2 = 3e-20")
    (is (c-rel-eq? [-5.0e-21 1.0e-10] (c/expm1 (c/complex 0.0 1.0e-10)) 1.0e-9) "real part -(1 - cos(b)) = -b^2/2"))
  (is (c/delta-eq c/ZERO (c/expm1 c/ZERO)))
  (is (= (c/complex ##-Inf 0.0) (c/log1p (c/complex -1.0 0.0)))))

(deftest roots
  (is (= (c/complex 2.0 0.0) (c/cbrt (c/complex 8.0 0.0))) "exact for perfect cubes")
  (is (c/delta-eq (c/complex 1.0 0.0) (c/pow (c/cbrt (c/complex 1.0 0.0)) (c/complex 3.0 0.0))))
  (is (c/delta-eq (c/complex 1.25992105 0.0) (c/cbrt (c/complex 2.0 0.0)) 1.0e-8))
  (is (c/delta-eq (c/complex 0.5 (/ (m/sqrt 3.0) 2.0)) (c/cbrt (c/complex -1.0 0.0)) 1.0e-15) "principal root of -1")
  (is (c/delta-eq (c/complex 0.5 (- (/ (m/sqrt 3.0) 2.0))) (c/cbrt (c/complex -1.0 -0.0)) 1.0e-15) "the other side of the cut")
  (doseq [n [1 2 3 5 8 -2 -3]
          z [(c/complex 3.0 -4.0) (c/complex -0.3 7.0) (c/complex 1e-5 -1e5) (c/complex -2.0 0.0)]]
    (let [r (c/nth-root z n)
          back (c/pow r (c/complex (double n) 0.0))]
      (is (<= (c/abs (c/sub z back)) (* 1.0e-13 (c/abs z))) (str "root " n " of " z))
      (is (<= (m/abs (c/arg r)) (+ (/ m/PI (m/abs n)) 1.0e-15)) "principal root")))
  (is (= c/ZERO (c/nth-root c/ZERO 3)))
  (is (c/delta-eq (c/complex 1.0 2.0) (c/nth-root (c/complex 1.0 2.0) 1) 1.0e-15))
  (is (thrown-with-msg? clojure.lang.ExceptionInfo #"non zero n" (c/nth-root c/ONE 0))))

(deftest variadic-arithmetic
  (let [a (c/complex 1.0 2.0) b (c/complex -3.0 0.5) d (c/complex 0.25 -4.0)]
    (is (= c/ZERO (c/add)))
    (is (= c/ONE (c/mult)))
    (is (= a (c/add a)))
    (is (= a (c/mult a)))
    (is (= (c/neg a) (c/sub a)))
    (is (= (c/add (c/add a b) d) (c/add a b d)))
    (is (= (c/add (c/add (c/add a b) d) a) (c/add a b d a)))
    (is (= (c/sub (c/sub a b) d) (c/sub a b d)))
    (is (= (c/mult (c/mult a b) d) (c/mult a b d)))
    (is (= (c/mult (c/mult (c/mult a b) d) a) (c/mult a b d a)))))

(deftest integer-powers
  (is (= (c/complex -11.0 -2.0) (c/pow (c/complex 1 2) (c/complex 3 0))) "exact")
  (is (= (c/complex -3.0 4.0) (c/pow (c/complex 1 2) (c/complex 2 0))))
  (is (= (c/complex 1.0 0.0) (c/pow (c/complex 0 1) (c/complex 4 0))))
  (is (= (c/complex 0.0 -1.0) (c/pow (c/complex 0 1) (c/complex -1 0))))
  (is (c/delta-eq (c/reciprocal (c/complex -11.0 -2.0)) (c/pow (c/complex 1 2) (c/complex -3 0)) 1.0e-15))
  (is (= (c/complex 16.0 0.0) (c/pow (c/complex -2.0 0.0) (c/complex 4 0))))
  (is (= (c/complex -8.0 0.0) (c/pow (c/complex -2.0 0.0) (c/complex 3 0))))
  (testing "the sign of the zero imaginary part"
    (is (same-bits? -0.0 (c/im (c/pow (c/complex -2.0 0.0) (c/complex 4 0)))) "d/dy of z^4 at -2 is negative, as Julia")
    (is (same-bits? 0.0 (c/im (c/pow (c/complex -2.0 -0.0) (c/complex 4 0)))))
    (is (same-bits? 0.0 (c/im (c/pow (c/complex -2.0 0.0) (c/complex 3 0))))))
  (testing "no NaN from 0 * Inf"
    (is (= (c/complex ##Inf 0.0) (c/pow (c/complex ##Inf 0.0) (c/complex 2 0))))
    (is (= (c/complex 0.0 0.0) (c/pow (c/complex ##Inf 0.0) (c/complex -2 0)))))
  (testing "large exponents use the logarithm"
    (is (c-rel-eq? [1.0e10 0.0] (c/pow (c/complex 10.0 0.0) (c/complex 10 0)) 1.0e-14))))

(deftest pow-zero-base
  (is (c/zero? (c/pow c/ZERO (c/complex 2.0 0.0))))
  (is (c/zero? (c/pow (c/complex -0.0 0.0) (c/complex 2.0 0.0))) "was (Inf, Inf)")
  (is (c/zero? (c/pow c/ZERO (c/complex 0.5 3.0))))
  (let [w (c/pow c/ZERO (c/complex -2.0 1.0))]
    (is (= ##Inf (c/re w)) "complex infinity for a negative real part of the exponent")
    (is (m/nan? (c/im w))))
  (is (c/nan? (c/pow c/ZERO (c/complex 0.0 1.0))) "pure imaginary exponent")
  (is (c/nan? (c/pow c/ZERO c/ZERO)))
  (is (= c/ONE (c/pow (c/complex 5.0 -1.0) c/ZERO))))


(deftest constants-and-utilities
  (testing "constants"
    (is (= (c/complex 0.0 -1.0) c/-I c/I-))
    (is (= (c/complex m/PI 0.0) c/PI))
    (is (= (c/complex 2.0 0.0) c/TWO))
    (is (= c/ONE (c/mult c/I c/-I)) "i * (-i) = 1")
    (is (= (c/neg c/ONE) (c/mult c/I c/I))))
  (testing "predicates"
    (is (c/real? (c/complex 3.0 0.0)))
    (is (c/real? (c/complex 3.0 -0.0)))
    (is (not (c/real? (c/complex 3.0 1.0))))
    (is (c/imaginary? (c/complex 0.0 3.0)))
    (is (c/imaginary? (c/complex -0.0 3.0)))
    (is (not (c/imaginary? (c/complex 1.0 3.0))))
    (is (c/inf? (c/complex 1.0 ##-Inf)))
    (is (c/inf? (c/complex ##Inf ##NaN)))
    (is (not (c/inf? (c/complex 1.0 2.0))))
    (is (c/invalid? (c/complex ##NaN 1.0)))
    (is (c/invalid? (c/complex 1.0 ##Inf)))
    (is (c/valid? (c/complex 1.0 2.0))))
  (testing "ensure-complex"
    (is (= (c/complex 3.0 0.0) (c/ensure-complex 3)))
    (is (= (c/complex (/ 1.0 3.0) 0.0) (c/ensure-complex 1/3)))
    (is (= (c/complex 1.0 2.0) (c/ensure-complex (c/complex 1.0 2.0))))
    (is (= (c/complex 1.0 2.0) (c/ensure-complex (org.apache.commons.math3.complex.Complex. 1.0 2.0))))
    (is (= "x" (c/ensure-complex "x")) "other values are returned"))
  (testing "scaling and real parts"
    (is (= (c/complex 3.0 2.0) (c/adds (c/complex 1.0 2.0) 2.0)))
    (is (= (c/complex 2.0 -4.0) (c/scale (c/complex 1.0 -2.0) 2.0)))
    (is (= (c/complex -2.0 1.0) (c/mult-I (c/complex 1.0 2.0))))
    (is (= (c/complex 2.0 -1.0) (c/mult-I- (c/complex 1.0 2.0))))
    (is (= (c/mult c/I (c/complex 1.0 2.0)) (c/mult-I (c/complex 1.0 2.0))))
    (is (= (c/mult c/I- (c/complex 1.0 2.0)) (c/mult-I- (c/complex 1.0 2.0)))))
  (testing "csgn"
    (is (= 1.0 (c/csgn (c/complex 2.0 -5.0))))
    (is (= -1.0 (c/csgn (c/complex -2.0 5.0))))
    (is (= 1.0 (c/csgn (c/complex 0.0 -5.0))))
    (is (= 0.0 (c/csgn c/ZERO)))
    (is (= 0.0 (c/csgn 0.0 -0.0)))
    (is (= -1.0 (c/csgn -3.0 0.0))))
  (testing "muladd"
    (is (= (c/complex 3.0 7.0) (c/muladd (c/complex 1 1) (c/complex 2 2) (c/complex 3 3))))
    (is (= (c/add (c/complex 3 3) (c/mult (c/complex 1.5 -2) (c/complex 2 0.5))) (c/muladd (c/complex 1.5 -2) (c/complex 2 0.5) (c/complex 3 3)))))
  (testing "normalize"
    (is (c/delta-eq (c/complex 0.6 0.8) (c/normalize (c/complex 3.0 4.0)) 1.0e-15))
    (is (= c/ZERO (c/normalize c/ZERO)))
    (is (c/delta-eq (c/complex 0.6 0.8) (c/normalize (c/complex 3.0e200 4.0e200)) 1.0e-15) "no overflow (was 0)")
    (is (c/delta-eq (c/complex 0.6 0.8) (c/normalize (c/complex 3.0e-200 4.0e-200)) 1.0e-15)))
  (testing "logb"
    (is (c/delta-eq (c/complex 3.0 0.0) (c/logb 8 2) 1.0e-14) "numbers")
    (is (c/delta-eq (c/complex 3.0 0.0) (c/logb (c/complex 8.0 0.0) (c/complex 2.0 0.0)) 1.0e-14))
    (is (c/delta-eq (c/complex 1.0 0.0) (c/logb (c/complex 0.0 1.0) (c/complex 0.0 1.0)) 1.0e-14))
    (is (c/delta-eq (c/complex 0.5 0.0) (c/logb 2 4) 1.0e-14)))
  (testing "arg and abs-sq"
    (is (= m/PI (c/arg (c/complex -1.0 0.0))))
    (is (= (- m/PI) (c/arg (c/complex -1.0 -0.0))))
    (is (= m/HALF_PI (c/arg c/I)))
    (is (= 25.0 (c/abs-sq (c/complex 3.0 -4.0))))))


(deftest poles-at-zero
  ;; the limits along the real axis (Octave: asec(0) = 0+Inf i, acsc(0) = pi/2-Inf i, asech(0) = Inf, acsch(0) = Inf,
  ;; acot(0) = pi/2, acoth(0) = +-pi/2 i), the signed zeros select the side
  (are [f x y re im]
      (let [r (f (c/complex x y))]
        (and (= re (c/re r)) (same-bits? im (c/im r))))
    c/asec 0.0 0.0 0.0 ##Inf
    c/asec 0.0 -0.0 0.0 ##-Inf
    c/asec -0.0 0.0 m/PI ##Inf
    c/asec -0.0 -0.0 m/PI ##-Inf
    c/acsc 0.0 0.0 m/HALF_PI ##-Inf
    c/acsc 0.0 -0.0 m/HALF_PI ##Inf
    c/acsc -0.0 0.0 (- m/HALF_PI) ##-Inf
    c/asech 0.0 0.0 ##Inf -0.0
    c/asech 0.0 -0.0 ##Inf 0.0
    c/asech -0.0 0.0 ##Inf (- m/PI)
    c/acsch 0.0 0.0 ##Inf -0.0
    c/acsch 0.0 -0.0 ##Inf 0.0
    c/acsch -0.0 0.0 ##-Inf -0.0
    c/acot 0.0 0.0 m/HALF_PI -0.0
    c/acot -0.0 0.0 (- m/HALF_PI) -0.0
    c/acoth 0.0 0.0 0.0 (- m/HALF_PI)
    c/acoth 0.0 -0.0 0.0 m/HALF_PI)
  (testing "csc, cot, csch, coth keep the infinite real part and a zero imaginary part"
    (is (= (c/complex ##Inf -0.0) (c/csc c/ZERO)))
    (is (= (c/complex ##Inf 0.0) (c/cot c/ZERO)))
    (is (= (c/complex ##Inf -0.0) (c/csch c/ZERO)))
    (is (= (c/complex ##Inf 0.0) (c/coth c/ZERO)))))


(deftest log-near-unit-circle
  ;; |z|^2 - 1 is computed in double-double: the real part keeps its digits at any distance from the unit circle
  ;; reference: mpmath, 60 digits
  (are [a b re im] (c-rel-eq? [re im] (c/log (c/complex a b)) 1.0e-14)
    1.0 1.0e-10 5.0000000000000005e-21 1.0e-10
    1.0 -1.0e-10 5.0000000000000005e-21 -1.0e-10
    1.000000000001 0.0 1.000088900581841e-12 0.0
    0.999999999999 1.0e-6 -4.999778782796285e-13 1.0000000000006666e-6
    1.0e-5 0.99999999995 -4.135768541567362e-18 1.5707863267948965
    -0.9999999999999 1.0e-7 -9.503109451873063e-14 3.1415925535897933
    0.9999999 1.0e-4 -9.50000039723646e-8 1.0000000966666756e-4)
  (testing "the point on the unit circle given by a rounded pair is not exactly on it"
    (is (c-rel-eq? [0.0 1.0471975511965976] (c/complex 0.0 (c/im (c/log (c/complex 0.5 0.8660254037844386)))) 1.0e-15))
    (is (<= (m/abs (c/re (c/log (c/complex 0.5 0.8660254037844386)))) 1.0e-16) "-4.3e-17 for the double nearest to sqrt(3)/2")))

(deftest infinite-real-part-hyperbolic-tangent
  ;; numpy and Julia: tanh(+-Inf + Inf i) = +-1 + 0i, tan(Inf + Inf i) = 0 + i (C99)
  (are [x y re im] (let [w (c/tanh (c/complex x y))]
                     (and (= re (c/re w)) (same-bits? im (c/im w))))
    ##Inf ##Inf 1.0 0.0
    ##-Inf ##Inf -1.0 0.0
    ##Inf ##-Inf 1.0 -0.0
    ##Inf ##NaN 1.0 0.0)
  (let [w (c/tan (c/complex ##Inf ##Inf))]
    (is (same-bits? 0.0 (c/re w)))
    (is (= 1.0 (c/im w))))
  (let [w (c/tan (c/complex ##Inf ##-Inf))]
    (is (= -1.0 (c/im w))))
  (testing "a finite real part and an infinite imaginary part stay NaN"
    (is (c/nan? (c/tanh (c/complex 1.0 ##Inf))))
    (is (c/nan? (c/tanh (c/complex 30.0 ##Inf)))))
  (testing "a finite imaginary part: the sign of sin(2y)"
    (is (same-bits? 0.0 (c/im (c/tanh (c/complex ##Inf 1.0)))))
    (is (same-bits? -0.0 (c/im (c/tanh (c/complex ##Inf -1.0)))))
    (is (same-bits? -0.0 (c/im (c/tanh (c/complex ##Inf 2.0)))))))

(deftest nan-with-a-zero-part
  ;; C99 G.6.2.2 (numpy and Julia agree): only casinh keeps the zero imaginary part
  (is (m/nan? (c/re (c/asinh (c/complex ##NaN 0.0)))))
  (is (same-bits? 0.0 (c/im (c/asinh (c/complex ##NaN 0.0)))))
  (is (same-bits? -0.0 (c/im (c/asinh (c/complex ##NaN -0.0)))))
  (doseq [f [c/atanh c/acos c/asin]]
    (is (c/nan? (f (c/complex ##NaN 0.0))))
    (is (m/nan? (c/im (f (c/complex ##NaN 0.0)))))))

(deftest nth-root-of-zero
  (is (= c/ZERO (c/nth-root c/ZERO 2)))
  (let [w (c/nth-root c/ZERO -2)]
    (is (= ##Inf (c/re w)))
    (is (m/nan? (c/im w))))
  (is (= c/ZERO (c/nth-root (c/complex ##Inf 1.0) -3)) "the root of infinity is zero for a negative n")
  (is (= ##Inf (c/re (c/nth-root (c/complex ##Inf 1.0) 3)))))


(deftest infinite-and-nan-parts-c99
  ;; Test phase 2026-10-27, expected values: C99 Annex G, Julia ComplexF64 (numpy and Octave where noted)
  (testing "sqrt with an infinite imaginary part, whatever the real part (BR-5)"
    (are [x y re im] (let [w (c/sqrt (c/complex x y))]
                       (and (= re (c/re w)) (= im (c/im w))))
      -1.0 ##Inf ##Inf ##Inf
      -1.0 ##-Inf ##Inf ##-Inf
      1.0e300 ##Inf ##Inf ##Inf
      ##-Inf ##Inf ##Inf ##Inf
      ##Inf ##-Inf ##Inf ##-Inf)
    (is (= ##Inf (c/re (c/sqrt (c/complex ##NaN ##Inf)))) "even with a NaN real part")
    (is (= ##Inf (c/im (c/sqrt (c/complex ##NaN ##Inf))))))
  (testing "exp with an infinite real part and a non finite imaginary part (BR-6)"
    (are [x y re im] (let [w (c/exp (c/complex x y))]
                       (and (= re (c/re w)) (if (m/nan? im) (m/nan? (c/im w)) (same-bits? im (c/im w)))))
      ##Inf ##Inf ##Inf ##NaN
      ##Inf ##-Inf ##Inf ##NaN
      ##Inf ##NaN ##Inf ##NaN
      ##-Inf ##Inf -0.0 0.0
      ##-Inf ##-Inf -0.0 -0.0
      ##-Inf ##NaN -0.0 0.0)
    (is (c/nan? (c/exp (c/complex 1.0 ##Inf))) "a finite real part: NaN"))
  (testing "sin, cos, sinh, cosh with an infinite real part and a zero imaginary part (BR-7)"
    (are [f x y im] (let [w (f (c/complex x y))]
                      (and (m/nan? (c/re w)) (same-bits? im (c/im w))))
      c/sin ##Inf 0.0 0.0
      c/sin ##-Inf -0.0 -0.0
      c/sin ##NaN 0.0 0.0
      c/cos ##Inf 0.0 -0.0
      c/cos ##-Inf 0.0 0.0
      c/cos ##Inf -0.0 0.0)
    (is (= (c/complex ##Inf 0.0) (c/cosh (c/complex ##Inf 0.0))))
    (is (= (c/complex ##-Inf 0.0) (c/sinh (c/complex ##-Inf 0.0)))))
  (testing "infinite parts of sin, cos, sinh, cosh"
    (let [w (c/sin (c/complex ##Inf ##-Inf))] (is (m/nan? (c/re w))) (is (= ##-Inf (c/im w))))
    (let [w (c/cos (c/complex ##-Inf ##Inf))] (is (= ##Inf (c/re w))) (is (m/nan? (c/im w))))
    (let [w (c/sinh (c/complex ##-Inf ##Inf))] (is (= ##-Inf (c/re w))) (is (m/nan? (c/im w))))
    (let [w (c/cosh (c/complex ##-Inf ##NaN))] (is (= ##Inf (c/re w))) (is (m/nan? (c/im w))))
    (let [w (c/sinh (c/complex -0.0 ##NaN))] (is (same-bits? -0.0 (c/re w))) (is (m/nan? (c/im w))) "csinh(+-0 + iNaN) = +-0 + iNaN"))
  (testing "atanh and atan with an infinite imaginary part or real part (BR-3)"
    (are [f x y re im] (let [w (f (c/complex x y))]
                         (and (same-bits? re (c/re w)) (same-bits? im (c/im w))))
      c/atanh -1.0e300 ##-Inf -0.0 (- m/HALF_PI)
      c/atanh -1.0 ##Inf -0.0 m/HALF_PI
      c/atanh 1.0 ##-Inf 0.0 (- m/HALF_PI)
      c/atanh 0.5 ##Inf 0.0 m/HALF_PI
      c/atan ##-Inf 1.0 (- m/HALF_PI) 0.0
      c/atan ##-Inf -1.0 (- m/HALF_PI) -0.0
      c/atan ##Inf 1.0e300 m/HALF_PI 0.0))
  (testing "NaN branches of asinh, asin and acos (BR-4)"
    (are [f x y re im] (let [w (f (c/complex x y))]
                         (and (if (m/nan? re) (m/nan? (c/re w)) (= re (c/re w)))
                              (if (m/nan? im) (m/nan? (c/im w)) (= im (c/im w)))))
      c/asinh ##-Inf ##NaN ##-Inf ##NaN
      c/asinh ##Inf ##NaN ##Inf ##NaN
      c/asin ##NaN ##-Inf ##NaN ##-Inf
      c/asin ##NaN ##Inf ##NaN ##Inf
      c/acos 0.0 ##NaN m/HALF_PI ##NaN
      c/acos -0.0 ##NaN m/HALF_PI ##NaN)))
