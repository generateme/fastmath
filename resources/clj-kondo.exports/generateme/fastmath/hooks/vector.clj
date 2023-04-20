(ns hooks.vector)

(defmacro primitive-ops
  [ops]
  `(declare ~@(for [o ops]
                `~(symbol (name o)))))
