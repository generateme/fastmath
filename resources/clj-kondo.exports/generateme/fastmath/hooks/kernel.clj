(ns hooks.kernel)

(defmacro emit-simple
  [name form]
  `(let [~'xs 0] [~name ~form]))

(defmacro emit-beta
  [name form]
  `(let [~'xs 0 ~'beta 0] [~name ~form]))
