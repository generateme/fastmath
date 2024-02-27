(ns fastmath.fields.k
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]
           [fastmath.fields.utils Tetrad]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn kaleidoscope
  ([] {:type :regular
       :config (fn [] {:pull (r/drand -1.0 1.0)
                      :x (r/drand -1.0 1.0)
                      :y (r/drand -1.0 1.0)
                      :line-up (r/drand -1.0 1.0)
                      :rotate (u/sdrand 0.25 1.5)})})
  ([^double amount {:keys [^double pull ^double line-up ^double rotate ^double x ^double y]}]
   (let [w (* 0.5253219888177296 rotate) ;; cos(45)
         qe (+ pull line-up)] 
     (fn [^Vec2 v]
       (let [xx (* w (.x v))
             yy (* 0.8509035245341184 (.y v))]
         (Vec2. (* amount (+ (- xx yy) x))
                (let [a (+ xx yy)]
                  (* amount (if (pos? (.y v))
                              (+ a qe y)
                              (- a qe))))))))))

(defn- complex-inverse
  [^Tetrad m]
  (Tetrad. (.d m) (c/scale (.b m) -1.0) (c/scale (.c m) -1.0) (.a m)))

(defn- kleingroup-jorgensen
  [traceA traceB]
  (let [b (c/scale (c/mult traceA traceB) -1.0)
        c (c/add (c/mult traceA traceA) (c/mult traceB traceB))
        bsq (c/mult b b)
        ac4 (c/scale c 4.0)
        traceAB (-> (c/scale b -1.0)
                    (c/sub (c/sqrt (c/sub bsq ac4)))
                    (c/div c/TWO))]
    [(Tetrad. (c/sub traceA (c/div traceB traceAB))
              (c/div traceA (c/mult traceAB traceAB))
              traceA
              (c/div traceB traceAB))
     (Tetrad. (c/sub traceB (c/div traceA traceAB))
              (c/div (c/scale traceB -1.0) (c/mult traceAB traceAB))
              (c/scale traceB -1.0)
              (c/div traceA traceAB))]))

(defn- kleingroup-riley
  [c _]
  [(Tetrad. c/ONE c/ZERO c c/ONE)
   (Tetrad. c/ONE c/TWO c/ZERO c/ONE)])

(defn- kleingroup-modifiedriley
  [c b1]
  [(Tetrad. c/ONE c/ZERO c c/ONE)
   (Tetrad. c/ONE b1 c/ZERO c/ONE)])

(defn- kleingroup-maskit
  [mu _]
  [(Tetrad. (c/mult (c/scale mu -1.0) c/I) c/-I c/-I c/ZERO)
   (Tetrad. c/ONE c/TWO c/ZERO c/ONE)])

(defn- kleingroup-modifiedmaskit
  [mu b1]
  [(Tetrad. (c/mult (c/scale mu -1.0) c/I) c/-I c/-I c/ZERO)
   (Tetrad. c/ONE b1 c/ZERO c/ONE)])

(defn- kleingroup-maskitleysmodified
  [mu b1]
  (let [b1 (c/complex (* 2.0 (m/cos (/ m/PI (c/re b1)))) (c/im b1))]
    [(Tetrad. (c/mult (c/scale mu -1.0) c/I) c/-I c/-I c/ZERO)
     (Tetrad. c/ONE b1 c/ZERO c/ONE)]))

(def ^:private ITWO (Vec2. 0.0 2.0))
(def ^:pricate FOUR (Vec2. 4.0 0.0))
(def ^:private IFOUR (Vec2. 0.0 4.0))

(defn- kleingroup-grandma
  [traceA traceB]
  (let [b (c/scale (c/mult traceA traceB) -1.0)
        c (c/add (c/mult traceA traceA) (c/mult traceB traceB))
        bsq (c/mult b b)
        ac4 (c/scale c 4.0)
        traceAB (-> (c/scale b -1.0)
                    (c/sub (c/sqrt (c/sub bsq ac4)))
                    (c/div c/TWO))
        z0 (c/div (-> (c/sub traceAB c/TWO) (c/mult traceB))
                  (-> (c/mult traceB traceAB)
                      (c/sub (c/mult traceA c/TWO))
                      (c/add (c/mult traceAB ITWO))))
        ]
    [(Tetrad. (c/div traceA c/TWO)
              (c/div (-> (c/mult traceA traceAB)
                         (c/sub (c/mult traceB c/TWO))
                         (c/add IFOUR))
                     (-> (c/mult traceAB c/TWO)
                         (c/add FOUR)
                         (c/mult z0)))
              (c/div (-> (c/mult traceA traceAB)
                         (c/sub (c/mult traceB c/TWO))
                         (c/sub IFOUR)
                         (c/mult z0))
                     (-> (c/mult traceAB c/TWO)
                         (c/sub FOUR)))
              (c/div traceA c/TWO))
     (Tetrad. (c/div (c/sub traceB ITWO) c/TWO)
              (c/div traceB c/TWO)
              (c/div traceB c/TWO)
              (c/div (c/add traceB ITWO) c/TWO))]))

(defn kleingroup
  ([] {:type :random
       :config (fn [] {:a-re (r/drand -0.7 0.7)
                      :a-im (r/randval 0.0 (r/drand -0.7 0.7))
                      :b-re (r/drand -0.7 0.7)
                      :b-im (r/randval 0.0 (r/drand -0.7 0.7))
                      :recipe (r/irand 7)})})
  ([^double amount {:keys [^double a-re ^double a-im ^double b-re ^double b-im
                           ^long recipe]}]
   (let [a (c/complex a-re a-im)
         b (c/complex b-re b-im)
         recipe ([kleingroup-grandma
                  kleingroup-maskit
                  kleingroup-jorgensen
                  kleingroup-riley
                  kleingroup-modifiedriley
                  kleingroup-modifiedmaskit
                  kleingroup-maskitleysmodified] (mod recipe 7))
         [mat-a mat-b] (recipe a b)
         mat-inv-a (complex-inverse mat-a)
         mat-inv-b (complex-inverse mat-b)
         all-matrices [mat-a mat-inv-a mat-b mat-inv-b]]
     (fn [^Vec2 v]
       (let [^Tetrad mat (rand-nth all-matrices)
             win (v/div v amount)]
         (v/mult (-> (c/mult win (.a mat))
                     (c/add (.b mat))
                     (c/div (-> (c/mult win (.c mat))
                                (c/add (.d mat))))) amount))))))

(m/unuse-primitive-operators)
