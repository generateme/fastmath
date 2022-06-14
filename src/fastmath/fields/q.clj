(ns fastmath.fields.q
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn q-ode
  ([] {:type :regular
       :config (fn [] {:q-ode01 (r/drand -0.5 0.5)
                      :q-ode02 (r/drand -0.5 0.5)
                      :q-ode03 (r/drand -0.5 0.5)
                      :q-ode04 (r/drand -0.5 0.5)
                      :q-ode05 (r/drand -0.5 0.5)
                      :q-ode06 (r/drand -0.5 0.5)
                      :q-ode07 (r/drand -0.5 0.5)
                      :q-ode08 (r/drand -0.5 0.5)
                      :q-ode09 (r/drand -0.5 0.5)
                      :q-ode10 (r/drand -0.5 0.5)
                      :q-ode11 (r/drand -0.5 0.5)
                      :q-ode12 (r/drand -0.5 0.5)})})
  ([^double amount {:keys [^double q-ode01 ^double q-ode02 ^double q-ode03 ^double q-ode04
                           ^double q-ode05 ^double q-ode06 ^double q-ode07 ^double q-ode08
                           ^double q-ode09 ^double q-ode10 ^double q-ode11 ^double q-ode12]}]
   (let [aq-ode02 (* amount q-ode02)
         aq-ode11 (* amount q-ode11)]
     (fn [^Vec2 v]
       (let [xx (* (.x v) (.x v))
             yy (* (.y v) (.y v))
             xy (* (.x v) (.y v))]
         (Vec2. (+ q-ode01 (* aq-ode02 (.x v)) (* q-ode03 xx)
                   (* q-ode04 xy) (* q-ode05 (.y v)) (* q-ode06 yy))
                (+ q-ode07 (* q-ode08 (.x v)) (* q-ode09 xx)
                   (* q-ode10 xy) (* aq-ode11 (.y v)) (* q-ode12 yy))))))))
