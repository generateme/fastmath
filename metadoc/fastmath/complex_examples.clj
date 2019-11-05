(ns fastmath.complex-examples
  (:require [fastmath.complex :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.core :as m]))

(add-examples complex (example "New complex number." (complex 2 -1)))
(add-examples abs (example "Abs" (abs (complex 1 -3))))
(add-examples add (example "Sum" (add I ONE)))
(add-examples sub (example "Subtract" (sub ONE I-)))
(add-examples arg (example "Argument" (m/degrees (arg I-))))
(add-examples conjugate (example "Conjugate" (conjugate I)))
(add-examples div (example "Divide" (div (complex 1 2) (complex 3 4))))
(add-examples reciprocal
  (example "Reciprocal of real" (reciprocal TWO))
  (example "Reciprocal of complex" (reciprocal (complex 0 2))))
(add-examples mult (example "Multiply" (mult (complex 1 2) (complex 3 4))))
(add-examples neg (example "Negate." (neg (complex 1 2))))

(add-examples sq (example "Square." (sq (complex 1 2)))
  (example "\\\\(i^2\\\\)" (sq I)))
(add-examples sqrt (example "Square root of real." (sqrt (complex 2 0)))
  (example "Square root of complex." (sqrt (complex 2 2))))
(add-examples sqrt1z (example "Example 1" (sqrt1z (complex 2 3))))

(add-examples cos (example "cos(z)" (cos (complex 2 -1))))
(add-examples sin (example "sin(z)" (sin (complex 2 -1))))
(add-examples cosh (example "cosh(z)" (cosh (complex 2 -1))))
(add-examples sinh (example "sinh(z)" (sinh (complex 2 -1))))
(add-examples tan (example "tan(z)" (tan (complex 2 -1))))
(add-examples tanh (example "tanh(z)" (tanh (complex 2 -1))))
(add-examples sec (example "sec(z)" (sec (complex 2 -1))))
(add-examples csc (example "csc(z)" (csc (complex 2 -1))))
(add-examples acos (example "acos(z)" (acos (complex 2 -1))))
(add-examples asin (example "asin(z)" (asin (complex 2 -1))))
(add-examples atan (example "atan(z)" (atan (complex 2 -1))))

(add-examples exp (example "exp(z)" (exp (complex 2 -1)))
  (example "\\\\(e^{i\\pi}+1\\\\)" (add (exp (complex 0 m/PI)) ONE)))
(add-examples log (example "log(z)" (log (complex 2 -1)))
  (example "log(e)" (log (complex m/E 0))))

(add-examples pow (example "\\\\(\\sqrt{2}\\\\)" (pow TWO (complex 0.5 0.0)))
  (example "Complex power" (pow (complex 1 2) (complex 3 4))))

(def fn-list `(atan asin acos log exp csc sec tanh tan sinh sin cosh cos sqrt sq sqrt1z reciprocal))

(defmacro ^:private add-image-examples
  []
  `(do
     ~@(for [x fn-list]
         `(add-examples ~x
            (example-image ~(str "Plot of " (name x)) ~(str "images/c/" (name x) ".jpg"))))))

(add-image-examples)
