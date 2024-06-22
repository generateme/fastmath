(ns fastmath.protocols.polynomials)

(defprotocol PolynomialProto
  (add [p1 p2])
  (negate [p])
  (scale [p v])
  (mult [p1 p2])
  (derivative [p] [p order])
  (coeffs [p])
  (degree [p])
  (evaluate [p x]))
