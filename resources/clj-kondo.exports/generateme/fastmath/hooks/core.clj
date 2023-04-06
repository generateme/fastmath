(ns hooks.core)

;; fastmath.core helpers
(defmacro proxy->declare [s & _] `(declare ~s))
(defmacro variadic-proxy->declare [_ n & _] `(declare ~n))
