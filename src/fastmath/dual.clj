(ns fastmath.dual)

(defprotocol DualProto
  (get-value [d])
  (get-partials [d])
  (get-tag [d])
  (as-dual [d]))

(defrecord Dual [value partials tag]
  DualProto
  (get-value [_] value)
  (get-partials [_] partials)
  (get-tag [_] tag)
  (as-dual [d] d))

(extend-type Number
  DualProto
  (get-value [d] d)
  (get-partials [_] [])
  (get-tag [_] -1)
  (as-dual [d] (->Dual d [] -1)))


(def TestTag -3)
(def OuterTestTag -2)

(def PARTIALS [7 3])
(def PRIMAL 2)
(def FDNUM (->Dual PRIMAL PARTIALS TestTag))

(def PARTIALS2 [8 7])
(def PRIMAL2 9)
(def FDNUM2 (->Dual PRIMAL2 PARTIALS2 TestTag))

(def PARTIALS3 [7 10])
(def PRIMAL3 5)
(def FDNUM3 (->Dual PRIMAL3 PARTIALS3 TestTag))

(def M-PARTIALS [5 6 9])
(def NESTED-PARTIALS [(->Dual 7 [0 0 0] TestTag)
                      (->Dual 3 [0 0 0] TestTag)])
(def NESTED-FDNUM (->Dual (->Dual PRIMAL M-PARTIALS TestTag) NESTED-PARTIALS TestTag))

(def M-PARTIALS2 [6 4 7])
(def NESTED-PARTIALS2 [(->Dual 8 [0 0 0] TestTag)
                       (->Dual 7 [0 0 0] TestTag)])
(def NESTED-FDNUM (->Dual (->Dual PRIMAL2 M-PARTIALS2 TestTag) NESTED-PARTIALS2 TestTag))

