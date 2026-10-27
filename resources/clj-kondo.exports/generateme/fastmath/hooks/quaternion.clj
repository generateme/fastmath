(ns hooks.quaternion)

(defmacro gen-from-complex
  [op _description]
  `(defn ~op [z#] z#))
