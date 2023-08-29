(ns hooks.quaternion)

(defmacro gen-from-complex
  [op]
  `(defn ~op [z#] z#))
