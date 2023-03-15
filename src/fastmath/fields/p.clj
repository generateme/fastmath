(ns fastmath.fields.p
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.fields.utils :as u]
            [fastmath.complex :as c])
  (:import [fastmath.vector Vec2]
           [fastmath.fields.utils Pair]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn pdj
  "PDJ"
  ([] {:type :regular
       :config (fn [] {:a (r/drand -6.0 6.0)
                      :b (r/drand -6.0 6.0)
                      :c (r/drand -6.0 6.0)
                      :d (r/drand -6.0 6.0)})})
  ([^double amount {:keys [^double a ^double b ^double c ^double d]}]
   (fn [^Vec2 v]
     (Vec2. (* amount (- (m/sin (* a (.y v))) (m/cos (* b (.x v)))))
            (* amount (- (m/sin (* c (.x v))) (m/cos (* d (.y v)))))))))

(defn ptransform
  ([] {:type :regular
       :config (fn [] {:rotate (r/drand m/TWO_PI)
                      :power (u/sdrand 0.5 2.0)
                      :move (r/drand -1.0 1.0)
                      :split (r/drand -1.0 1.0)
                      :log? (r/brand)})})
  ([^double amount {:keys [^double rotate ^double power ^double move ^double split log?]}]
   (fn [^Vec2 v]
     (let [rho (+ (/ (if log?
                       (m/log (v/mag v))
                       (v/mag v)) power) move)
           theta (+ (v/heading v) rotate)
           rho (if (pos? (.x v))
                 (+ rho split)
                 (- rho split))
           arho (* amount (if log? (m/exp rho) rho))]
       (Vec2. (* arho (m/cos theta))
              (* arho (m/sin theta)))))))

(defn panorama1
  "Panorama1"
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [aux (/ (m/sqrt (inc (v/magsq v))))
           nv (v/mult v aux)
           aux (v/mag nv)]
       (Vec2. (* amount m/M_1_PI (v/heading nv))
              (* amount (- aux 0.5)))))))

(defn panorama2
  "Panorama2"
  ([] {:type :regular})
  ([^double amount _]
   (fn [v]
     (let [aux (/ (inc (m/sqrt (v/magsq v))))
           nv (v/mult v aux)
           aux (v/mag nv)]
       (Vec2. (* amount m/M_1_PI (v/heading nv))
              (* amount (- aux 0.5)))))))
(defn parabola
  "Parabola fn"
  ([] {:type :random
       :config (fn [] {:width (u/sdrand 0.5 2.0)
                      :height (u/sdrand 0.5 2.0)})})
  ([^double amount {:keys [^double width ^double height]}]
   (fn [v]
     (let [r (v/mag v)
           sr (m/sin r)
           cr (m/cos r)]
       (Vec2. (* amount height sr sr (r/drand))
              (* amount width cr (r/drand)))))))

(defn parallel
  ([] {:type :random
       :config (fn [] {:x1width (u/sdrand 0.2 4.0)
                      :x1tilesize (u/sdrand 0.2 2.0)
                      :x1mod1 (u/sdrand 0.1 2.0)
                      :x1mod2 (u/sdrand 0.1 2.0)
                      :x1height (u/sdrand 0.2 4.0)
                      :x1move (r/drand -1.0 1.0)
                      :x2width (u/sdrand 0.2 4.0)
                      :x2tilesize (u/sdrand 0.2 2.0)
                      :x2mod1 (u/sdrand 0.1 2.0)
                      :x2mod2 (u/sdrand 0.1 2.0)
                      :x2height (u/sdrand 0.2 4.0)
                      :x2move (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double x1width ^double x1tilesize ^double x1mod1
                           ^double x1mod2 ^double x1height ^double x1move
                           ^double x2width ^double x2tilesize ^double x2mod1
                           ^double x2mod2 ^double x2height ^double x2move]}]
   (let [xr1 (* x1mod2 x1mod1)
         xr2 (* x2mod2 x2mod1)]
     (fn [^Vec2 v]
       (r/randval
        (let [x1 (r/randval x1width (- x1width))
              amove (* amount x1move)]
          (Vec2. (* x1tilesize (+ (.x v) (m/round (* x1 (m/log (r/drand))))))
                 (cond
                   (> (.y v) x1mod1) (+ (* x1height (- (mod (+ (.y v) x1mod1) xr1) x1mod1)) amove)
                   (< (.y v) (- x1mod1)) (+ (* x1height (- x1mod1 (mod (- x1mod1 (.y v)) xr1)) amove))
                   :else (+ (* x1height (.y v)) amove))))
        (let [x2 (r/randval x2width (- x2width))
              amove (* amount x2move)]
          (Vec2. (* x2tilesize (+ (.x v) (m/round (* x2 (m/log (r/drand))))))
                 (cond
                   (> (.y v) x2mod1) (+ (* x2height (- (mod (+ (.y v) x2mod1) xr2) x2mod1)) amove)
                   (< (.y v) (- x2mod1)) (+ (* x2height (- x2mod1 (mod (- x2mod1 (.y v)) xr2)) amove))
                   :else (+ (* x2height (.y v)) amove)))))))))

(defn perspective
  "Perspective"
  ([] {:type :regular
       :config (fn [] {:angle (u/sdrand 0.2 m/PI)
                      :dist (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double angle ^double dist]}]
   (let [ang (* m/HALF_PI angle)
         vsin (m/sin ang)
         vfcos (* dist (m/cos ang))]
     (fn [^Vec2 v]
       (let [t (/ amount (- dist (* (.y v) vsin)))]
         (Vec2. (* t dist (.x v))
                (* t vfcos (.y v))))))))

(defn petal
  "Petal"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [a (m/cos (.x v))
           bx (m/pow (* (m/cos (.x v)) (m/cos (.y v))) 3.0)
           by (m/pow (* (m/sin (.x v)) (m/cos (.y v))) 3.0)]
       (Vec2. (* amount a bx)
              (* amount a by))))))

(defn phoenix-julia
  "Phoenix julia"
  ([] {:type :random
       :config (fn [] {:power (u/sirand 1 10)
                      :dist (r/drand -2.0 2.0)
                      :x-distort (r/drand -2.0 2.0)
                      :y-distort (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double power ^double dist ^double x-distort ^double y-distort]}]
   (let [inv-n (/ dist power)
         inv2pi-n (/ m/TWO_PI power)
         c-n (* 0.5 inv-n)]
     (fn [^Vec2 v]
       (let [pre-x (* (.x v) x-distort)
             pre-y (* (.y v) y-distort)
             a (+ (* (m/atan2 pre-y pre-x) inv-n)
                  (* (r/lrand) inv2pi-n))
             sina (m/sin a)
             cosa (m/cos a)
             r (* amount (m/pow (v/magsq v) c-n))]
         (Vec2. (* r cosa)
                (* r sina)))))))

(defn pie
  "pie from jwildfire"
  ([] {:type :pattern
       :config (fn [] {:slices (r/randval (u/sirand 1 10) (u/sdrand 0.01 10.0))
                      :rotation (r/drand m/TWO_PI)
                      :thickness (r/drand -1.0 1.0)})})
  ([^double amount {:keys [^double slices ^double rotation ^double thickness]}]
   (fn [_]
     (let [sl (m/round (+ 0.5 (r/drand slices)))
           a (-> (r/drand thickness)
                 (+ sl)
                 (* m/TWO_PI)
                 (/ slices)
                 (+ rotation))
           r (* amount (r/drand))]
       (Vec2. (* r (m/cos a))
              (* r (m/sin a)))))))

(defn- pixelflow-hash
  ^double [^long a]
  (as-> (unchecked-int a) a
    (bit-xor (bit-xor a 61)
             (bit-shift-right a 16))
    (+ a (bit-shift-left a 3))
    (bit-xor a (bit-shift-right a 4))
    (unchecked-int (* a 0x27d4eb2d))
    (bit-xor a (bit-shift-right a 15))
    (/ (double (bit-and a 0xffffffff)) Integer/MAX_VALUE)))

(defn pixelflow
  "Pixel Flow"
  ([] {:type :random
       :config (fn [] {:angle (r/drand m/TWO_PI)
                      :len (u/sdrand 0.01 3.0)
                      :width (u/sdrand 0.1 10.0)
                      :seed (r/lrand)})})
  ([^double amount {:keys [^double angle ^double len ^double width ^long seed]}]
   (let [vin (Vec2. (m/cos angle) (m/sin angle))
         seed2 (/ seed 2)
         -seed (- seed)]
     (fn [^Vec2 v]
       (let [blockx (int (m/floor (* width (.x v))))
             blockx (+ blockx (- 2.0 (* 4.0 (pixelflow-hash (inc (* seed blockx))))))
             blocky (int (m/floor (* width (.y v))))
             blocky (+ blocky (- 2.0 (* 4.0 (pixelflow-hash (inc (* seed blocky))))))
             flen (* 0.5 (+ (pixelflow-hash (+ blocky (* blockx -seed)))
                            (pixelflow-hash (+ blockx (* blocky seed2)))))]
         (v/add v (v/mult vin (* amount flen (m/sq (m/sq (r/drand))) len))))))))

(defn plusrecip
  ([] {:type :regular
       :config (fn [] {:ar (u/sdrand 0.1 5.0)
                      :ai (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double ar ^double ai orig?]}]
   (let [a (c/complex ar ai)]
     (fn [^Vec2 z]
       (let [aa (m/sqrt (+ (v/magsq z) m/EPSILON))
             k (-> (c/sq z)
                   (c/sub a)
                   (c/sqrt)
                   (c/add z))
             k (if (< (m/sqrt (+ (v/magsq (c/sq k)) m/EPSILON)) aa)
                 (-> (c/conjugate k)
                     (c/mult (c/scale a (/ -1.0 aa))))
                 k)
             k (if (neg? (c/re k)) (c/neg k) k)]
         (c/scale k amount))))))

(defn polar2
  "Polar2"
  ([] {:type :regular})
  ([^double amount _]
   (let [p2v (/ amount m/PI)
         p2v2 (* 0.5 p2v)]
     (fn [^Vec2 v] (Vec2. (* p2v (v/heading v)) (* p2v2 (m/log (v/magsq v))))))))

(defn polar
  "Polar"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [ny (dec (v/mag v))]
       (Vec2. (* amount (v/heading v) m/M_1_PI)
              (* amount ny))))))

(def ^:private rznegodd [-0.08333333333333333,
                       0.008333333333333333,
                       -0.003968253968253968,
                       0.004166666666666667,
                       -0.007575757575757576,
                       0.021092796092796094,
                       -0.08333333333333333,
                       0.4432598039215686,
                       -3.0539543302701198,
                       26.456212121212122,
                       -281.46014492753625,
                       3607.5105463980462,
                       -54827.583333333336,
                       974936.8238505747,
                       -2.005269579668808e+7,
                       4.723848677216299e+8,
                       -1.2635724795916667e+10,
                       3.8087931125245367e+11,
                       -1.2850850499305083e+13,
                       4.824144835485017e+14,
                       -2.0040310656516253e+16,
                       9.167743603195331e+17,
                       -4.5979888343656505e+19,
                       2.518047192145109e+21,
                       -1.500173349215393e+23,
                       9.689957887463594e+24,
                       -6.764588237929281e+26,
                       5.089065946866229e+28,
                       -4.114728879255798e+30,
                       3.5666582095375554e+32,
                       -3.306608987657758e+34,
                       3.2715634236478713e+36,
                       -3.4473782558278054e+38,
                       3.861427983270526e+40,
                       -4.589297443245433e+42,
                       5.777538634277042e+44,
                       -7.691985875950713e+46,
                       1.0813635449971654e+49,
                       -1.6029364522008965e+51,
                       2.501947904156046e+53,
                       -4.106705233581021e+55,
                       7.079877440849458e+57,
                       -1.280454688793951e+60,
                       2.4267340392333522e+62,
                       -4.8143218874045757e+64,
                       9.987557417572751e+66,
                       -2.1645634868435182e+69,
                       4.8962327039620546e+71,
                       -1.1549023923963518e+74,
                       2.83822495706937e+76,
                       -7.26120088036067e+78,
                       1.932351423341981e+81,
                       -5.345016042528861e+83,
                       1.535602884642242e+86,
                       -4.5789872682265786e+88,
                       1.4162025212194806e+91,
                       -4.540065229609264e+93,
                       1.5076656758807855e+96,
                       -5.183094914826456e+98,
                       1.843564742725653e+101,
                       -6.780555475309095e+103,
                       2.57733267027546e+106,
                       -1.01191128757046e+109,
                       4.101634616154228e+111,
                       -1.715524453403202e+114,
                       7.400342570526908e+116,
                       -3.290922535705444e+119,
                       1.507983153416477e+122,
                       -7.116987918825454e+124,
                       3.458042914157776e+127,
                       -1.729090760667674e+130,
                       8.893699169503295e+132,
                       -4.7038470619636e+135,
                       2.55719382310602e+138,
                       -1.428406750044352e+141,
                       8.195215221831376e+143,
                       -4.827648542272736e+146,
                       2.918961237477031e+149,
                       -1.81089321625689e+152,
                       1.152357722002117e+155,
                       -7.519231195198174e+157,
                       5.029401657641103e+160,
                       -3.447342044447767e+163,
                       2.420745864586851e+166,
                       -1.740946592037767e+169,
                       1.281948986348224e+172,
                       -9.662412110856088e+174,
                       7.452691030430086e+177,
                       -5.880839331167436e+180,
                       4.74627186549076e+183,
                       -3.916913259477282e+186,
                       3.304507144322603e+189,
                       -2.849289055099457e+192,
                       2.510332934507758e+195,
                       -2.259390199547524e+198,
                       2.07691380042876e+201,
                       -1.949473217492725e+204,
                       1.868073147126591e+207,
                       -1.827075266281457e+210,
                       1.823538632259567e+213])

(def ^:private rzposint [-0.50000000000000000000000000000,
                       0.0, 
                       1.64493406684822643647241516665, 
                       1.20205690315959428539973816151,
                       1.08232323371113819151600369654,
                       1.03692775514336992633136548646,
                       1.01734306198444913971451792979,
                       1.00834927738192282683979754985,
                       1.00407735619794433937868523851,
                       1.00200839282608221441785276923,
                       1.00099457512781808533714595890,
                       1.00049418860411946455870228253,
                       1.00024608655330804829863799805,
                       1.00012271334757848914675183653,
                       1.00006124813505870482925854511,
                       1.00003058823630702049355172851,
                       1.00001528225940865187173257149,
                       1.00000763719763789976227360029,
                       1.00000381729326499983985646164,
                       1.00000190821271655393892565696,
                       1.00000095396203387279611315204,
                       1.00000047693298678780646311672,
                       1.00000023845050272773299000365,
                       1.00000011921992596531107306779,
                       1.00000005960818905125947961244,
                       1.00000002980350351465228018606,
                       1.00000001490155482836504123466,
                       1.00000000745071178983542949198,
                       1.00000000372533402478845705482,
                       1.00000000186265972351304900640,
                       1.00000000093132743241966818287,
                       1.00000000046566290650337840730,
                       1.00000000023283118336765054920,
                       1.00000000011641550172700519776])

(defn- polylogarithm-Z
  ^double [^long N]
  (cond
    (> N 33) 1.0
    (< N -200) 1.0e250
    (neg? N) (let [M (- N)]
               (if (even? M)
                 0.0
                 (rznegodd (/ (dec M) 2))))
    :else (rzposint N)))

(defn polylogarithm
  ([] {:type :regular
       :config (fn [] {:n (r/irand 1 21)
                      :zpow (u/sdrand 0.1 3.0)})})
  ([^double amount {:keys [^long n ^double zpow]}]
   (let [N (min 19 n)
         N- (dec N)
         cpow (c/complex zpow 0.0)
         HSTerm (c/complex (if (< N 1)
                             0.0
                             (reduce m/fast+ (map (fn [^double i]
                                                    (/ (inc i))) (range N)))) 0.0)
         series (mapv #(c/complex (m/pow % (- N)) 0.0) (range 1 20))]
     (fn [^Vec2 v]
       (let [z (-> v (c/pow cpow))
             zmag2 (v/magsq z)]
         (if (> zmag2 250000.0)
           (v/mult v amount)
           (if (< zmag2 0.07)
             (-> (->> (iterate #(c/mult % z) z)
                      (map c/mult series)
                      (reduce c/add))
                 (v/mult amount))
             (let [z (c/log z)
                   ^Pair zlin (loop [i (long 0)
                                     lin c/ZERO
                                     zpl1 c/ZERO
                                     z' z]
                                (if (< i 20)
                                  (let [rf (/ 1.0 (m/factorial20 i))
                                        zee (* rf (polylogarithm-Z (- N i)) )
                                        nlin (c/add lin (c/scale z' zee))
                                        nz (c/mult z' z)]
                                    (if (== i N-)
                                      (recur (inc i) nlin (c/scale z' rf) nz)
                                      (recur (inc i) nlin zpl1 nz)))
                                  (Pair. zpl1 lin)))]
               (-> (c/neg z)
                   (c/log)
                   (c/neg)
                   (c/add HSTerm)
                   (c/mult (.a zlin))
                   (c/add (.b zlin))
                   (v/mult amount))))))))))


(defn popcorn2
  "popcorn2 from apophysis"
  ([] {:type :regular
       :config (fn [] {:x (r/drand -1.5 1.5)
                      :y (r/drand -1.5 1.5)
                      :c (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double x ^double y ^double c]}]
   (fn [^Vec2 v]
     (let [xx (->> (.y v)
                   (* c)
                   (m/tan)
                   (m/sin)
                   (* x)
                   (+ (.x v))
                   (* amount))
           yy (->> (.x v)
                   (* c)
                   (m/tan)
                   (m/sin)
                   (* y)
                   (+ (.y v))
                   (* amount))]
       (Vec2. xx yy)))))

(defn popcorn
  ([] {:type :regular
       :config (fn [] {:coeff20 (r/drand -1.5 1.5)
                      :coeff21 (r/drand -1.5 1.5)})})
  ([^double amount {:keys [^double coeff20 ^double coeff21]}]
   (popcorn2 amount {:x coeff20 :y coeff21 :c 3.0})))

(defn powblock
  "PowBlock"
  ([] {:type :random
       :config (fn [] {:numerator (r/drand -20.0 20.0)
                      :denominator (r/drand -20.0 20.0)
                      :root (r/drand -6.0 6.0)
                      :correctn (r/drand -2.0 2.0)
                      :correctd (r/drand -2.0 2.0)})})
  ([^double amount {:keys [^double numerator ^double denominator ^double root
                           ^double correctn ^double correctd]}]
   (let [power (/ (* denominator correctn) (+ m/EPSILON (m/abs correctd)))
         power (if (< (m/abs power) m/EPSILON) m/EPSILON power)
         power (/ (* 0.5 numerator) power) 
         deneps (/ (if (< (m/abs denominator) m/EPSILON) m/EPSILON denominator))]
     (fn [^Vec2 v]
       (let [theta (v/heading v)
             r2 (* amount (m/pow (v/magsq v) power))
             ran (+ (* numerator (+ (* theta deneps) (* root m/TWO_PI (m/floor (r/drand denominator)) deneps))))]
         (Vec2. (* r2 (m/cos ran))
                (* r2 (m/sin ran))))))))
(defn power
  "Power"
  ([] {:type :regular})
  ([^double amount _]
   (fn [^Vec2 v]
     (let [theta (v/heading v)
           sa (m/sin theta)
           ca (m/cos theta)
           pow (* amount (m/pow (v/mag v) sa))]
       (Vec2. (* pow ca) (* pow sa))))))

(defn pressure-wave
  "Pressure Wave"
  ([] {:type :regular
       :config (fn [] {:x-freq (r/drand -6.0 6.0)
                      :y-freq (r/drand -6.0 6.0)})})
  ([^double amount {:keys [^double x-freq ^double y-freq]}]
   (let [[^double pwx ^double ipwx] (if (zero? x-freq)
                                      [1.0 1.0]
                                      (let [pwx (* x-freq m/TWO_PI)]
                                        [pwx (/ pwx)]))
         [^double pwy ^double ipwy] (if (zero? y-freq)
                                      [1.0 1.0]
                                      (let [pwy (* y-freq m/TWO_PI)]
                                        [pwy (/ pwy)]))]
     (fn [^Vec2 v]
       (Vec2. (* amount (+ (.x v) (* ipwx (m/sin (* pwx (.x v))))))
              (* amount (+ (.y v) (* ipwy (m/sin (* pwy (.y v)))))))))))

(defn projective
  ([] {:type :regular
       :config (fn [] {:A  (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :B  (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :C  (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :A1 (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :B1 (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :C1 (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :A2 (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :B2 (r/randval 0.3 0.0 (r/drand -1.5 1.5))
                      :C2 (r/randval 0.3 0.0 (r/drand -1.5 1.5))})})
  ([^double amount {:keys [^double A ^double B ^double C
                           ^double A1 ^double B1 ^double C1
                           ^double A2 ^double B2 ^double C2]}]
   (fn [^Vec2 v]
     (let [aU (/ amount (+ (* A (.x v)) (* B (.y v)) C))]
       (Vec2. (* aU (+ (* A1 (.x v)) (* B1 (.y v)) C1))
              (* aU (+ (* A2 (.x v)) (* B2 (.y v)) C2)))))))

(defn pulse
  ([] {:type :regular
       :config (fn [] {:freqx (r/randval 0.3 (u/sdrand 5.0 50.0) (u/sdrand 0.1 10.0))
                      :freqy (r/randval 0.3 (u/sdrand 5.0 50.0) (u/sdrand 0.1 10.0))
                      :scalex (u/sdrand 0.5 1.5)
                      :scaley (u/sdrand 0.5 1.5)})})
  ([^double amount {:keys [^double freqx ^double freqy ^double scalex ^double scaley]}]
   (let [freq (Vec2. freqx freqy)
         scale (Vec2. scalex scaley)]
     (fn [^Vec2 v]
       (-> (v/emult v freq)
           (v/sin)
           (v/emult scale)
           (v/add v)
           (v/mult amount))))))

;;

(defn perlin
  "Perlin noise"
  ([] {:type :regular
       :config (fn [] {:seed (r/irand)
                      :scale (u/sdrand 0.1 1.5)
                      :octaves (r/irand 1 6)})})
  ([amount {:keys [seed octaves scale]}]
   (let [n (r/fbm-noise {:seed seed :octaves octaves :gain 6})]
     (u/make-noise-variation amount scale 2.0 n))))

(defn perlin2
  "Perlin noise"
  ([] {:type :regular
       :config (fn [] {:seed (r/irand)
                      :scale (u/sdrand 0.1 1.5)
                      :octaves (r/irand 1 6)})})
  ([^double amount {:keys [^int seed ^int octaves ^double scale]}]
   (let [n (r/fbm-noise {:seed seed :octaves octaves})]
     (u/make-noise-variation2 amount scale n))))

(defn plusrecip2
  ([] {:type :regular
       :config (fn [] {:ar (u/sdrand 0.1 5.0)
                      :ai (r/drand -5.0 5.0)})})
  ([^double amount {:keys [^double ar ^double ai]}]
   (let [a (c/complex ar ai)]
     (fn [^Vec2 z]
       (let [aa (m/sqrt (+ (v/magsq z) m/EPSILON))
             k (-> (c/sq z)
                   (c/sub a)
                   (c/sqrt)
                   (c/add z))
             k (if (< (m/sqrt (+ (v/magsq (c/sq k)) m/EPSILON)) aa)
                 (-> (c/conjugate k)
                     (c/scale (/ -1.0 aa)))
                 k)]
         (c/scale k amount))))))
