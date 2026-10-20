(ns fastmath.kernel.registry-test
  (:require [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.kernel.vector :as vk]
            [fastmath.vector :as v]))

;; The registry (`k/kernel`, `k/rbf`) and the kernel combinators of `fastmath.kernel`.

(defn- ex-of
  "The `ExceptionInfo` thrown by `f`, nil when nothing is thrown; other exceptions are not caught."
  [f]
  (try (f) nil (catch clojure.lang.ExceptionInfo e e)))

(defn- close?
  [a b tol]
  (<= (m/abs (- (double a) (double b))) (* (double tol) (m/max 1.0 (m/abs (double b))))))

(def ^:private vector-keys (vec (sort (keys (methods k/kernel)))))
(def ^:private rbf-keys (vec (sort (keys (methods k/rbf)))))

;; documented minimal maps: the kernels which have no default for some parameter
(def ^:private vector-minimal {:geometric {:n 4}})

(def ^:private rbf-minimal
  {:truncated-power {:k 2.0}
   :gaussians-laguerre {:dimension 1.0 :degree 1}
   :poisson {:d 3.0}
   :matern {:order 3}
   :generalized-multiquadratic {:beta 0.5}
   :radial-powers {:beta 3.0}
   :thin-plate-splines {:beta 1.0}
   :shifted-surface-splines {:s 2 :beta 2.0}
   :wendland {:s 3.0 :k 1}
   :gneiting {:s 3.0 :l 4.0}
   :wu {:l 2.0 :k 1}
   :whittaker {:alpha 0.0 :k 2 :beta 1.0}})

(def ^:private vector-presets
  {:matern-12 {:order 1} :matern-32 {:order 3} :matern-52 {:order 5}
   :triangular {:n 1} :circular {:n 2} :spherical {:n 3}})

(def ^:private rbf-presets
  {:truncated-power-1 {:k 1.0} :truncated-power-2 {:k 2.0} :truncated-power-3 {:k 3.0}
   :truncated-power-half {:k 0.5} :truncated-power-third {:k m/THIRD}
   :gaussians-laguerre-11 {:dimension 1.0 :degree 1.0} :gaussians-laguerre-12 {:dimension 1.0 :degree 2.0}
   :gaussians-laguerre-21 {:dimension 2.0 :degree 1.0} :gaussians-laguerre-22 {:dimension 2.0 :degree 2.0}
   :poisson-2 {:d 2.0} :poisson-3 {:d 3.0} :poisson-4 {:d 4.0}
   :matern-c0 {:order 1.0} :matern-c2 {:order 3.0} :matern-c4 {:order 5.0}
   :multiquadratic {:beta 0.5 :negate? false} :inverse-multiquadratic {:beta -0.5 :negate? false}
   :radial-powers-3 {:beta 3.0 :negate? false} :thin-plate {:beta 1.0 :negate? false}
   :wu-10 {:l 1.0 :k 0.0} :wu-11 {:l 1.0 :k 1.0} :wu-20 {:l 2.0 :k 0.0} :wu-21 {:l 2.0 :k 1.0}
   :wu-22 {:l 2.0 :k 2.0} :wu-30 {:l 3.0 :k 0.0} :wu-31 {:l 3.0 :k 1.0} :wu-32 {:l 3.0 :k 2.0}
   :wu-33 {:l 3.0 :k 3.0}})

(defn- different
  "A value which differs from the preset value `x`."
  [x]
  (cond (boolean? x) (not x)
        (number? x) (+ 1.0 (double x))
        :else :other))

(defn- minimal-map
  [kind kw]
  (get (if (= kind :rbf) rbf-minimal vector-minimal) kw))

(deftest registry-coverage-test
  (testing "the minimal maps and presets name registered keywords"
    (is (every? (set vector-keys) (concat (keys vector-minimal) (keys vector-presets))))
    (is (every? (set rbf-keys) (concat (keys rbf-minimal) (keys rbf-presets)))))
  (testing "the registry holds 36 vector kernel names and 42 rbf names"
    (is (= 36 (count vector-keys)))
    (is (= 42 (count rbf-keys)))
    (is (= vector-keys (vec (:vector k/kernel-list))))
    (is (= rbf-keys (vec (:rbf k/kernel-list))))))

(deftest minimal-map-test
  (testing "every vector kernel keyword constructs with its documented minimal map and evaluates"
    (doseq [kw vector-keys
            :let [kf (k/kernel kw (minimal-map :vector kw))
                  val (kf [0.2 0.5] [0.4 0.1])]]
      (is (fn? kf) (str kw))
      (is (Double/isFinite (double val)) (str kw " " val))))
  (testing "every vector kernel gives the same value for 1-D numbers and one element vectors"
    (doseq [kw vector-keys
            :let [kf (k/kernel kw (minimal-map :vector kw))]]
      (is (= (kf 0.3 0.8) (kf [0.3] [0.8])) (str kw))))
  (testing "every rbf keyword constructs with its documented minimal map and evaluates"
    (doseq [kw rbf-keys
            :let [rf (k/rbf kw (minimal-map :rbf kw))
                  val (rf 0.5)]]
      (is (fn? rf) (str kw))
      (is (Double/isFinite (double val)) (str kw " " val)))))

(deftest non-map-params-test
  (testing "a number or any other non-map second argument throws ex-info naming the keyword and the value"
    (doseq [[multi kws] [[k/kernel vector-keys] [k/rbf rbf-keys]]
            kw kws
            bad [2.0 1 "a" [1.0] :a '(1) true]
            :let [e (ex-of #(multi kw bad))]]
      (is (some? e) (str kw " " bad))
      (is (= :invalid-input (:type (ex-data e))) (str kw " " bad))
      (is (= bad (:value (ex-data e))))
      (is (str/includes? (ex-message e) (str kw)))
      (is (str/includes? (ex-message e) (pr-str bad)))))
  (testing "the callers which passed a number (gp default kernel) and presets (ClassCastException before)"
    (is (some? (ex-of #(k/kernel :gaussian 1.0))))
    (is (some? (ex-of #(k/kernel :matern-52 1.5))))
    (is (some? (ex-of #(k/rbf :gaussian 2.0)))))
  (testing "nil and an empty map are accepted and equal to no argument"
    (doseq [kw vector-keys
            :let [base (minimal-map :vector kw)
                  x [0.2 0.5] y [0.4 0.1]
                  a ((k/kernel kw base) x y)]]
      (is (= a ((k/kernel kw (or base {})) x y)) (str kw))
      (when-not base
        (is (= a ((k/kernel kw nil) x y)) (str kw))
        (is (= a ((k/kernel kw) x y)) (str kw))))
    (is (= ((k/rbf :gaussian) 0.7) ((k/rbf :gaussian nil) 0.7) ((k/rbf :gaussian {}) 0.7)))
    (is (= ((k/rbf :linear) 0.7) ((k/rbf :linear nil) 0.7) ((k/rbf :linear {}) 0.7))))
  (testing "a map is used, scale applies"
    (is (close? ((k/kernel :gaussian {:sigma 2.0}) [0.0] [3.0]) (Math/exp (- (/ 9.0 8.0))) 1e-15))
    (is (close? ((k/rbf :gaussian {:scale 2.0}) 3.0) (Math/exp -2.25) 1e-15))))

(deftest missing-key-test
  (testing "a missing required parameter throws ex-info naming the key (vector kernels)"
    (doseq [p [nil {} {:r 1.0}]]
      (let [e (ex-of #(k/kernel :geometric p))]
        (is (some? e))
        (is (= :n (:option (ex-data e))))
        (is (str/includes? (ex-message e) ":n"))))
    (is (= :n (:option (ex-data (ex-of #(k/kernel :geometric))))))
    (is (some? (ex-of #(k/kernel :matern {:order 2})))))
  (testing "the registry passes the exception of the rbf constructor through unchanged"
    (doseq [[kw key] [[:truncated-power :k] [:poisson :d] [:matern :order] [:generalized-multiquadratic :beta]
                      [:wendland :s] [:wu :l]]
            p [nil {} {:scale 2.0}]
            :let [e (ex-of #(k/rbf kw p))]]
      (is (some? e) (str kw " " p))
      (is (= :invalid-input (:type (ex-data e))) (str kw " " p))
      (is (= key (:property (ex-data e))) (str kw " " p " " (ex-data e)))
      (is (str/includes? (ex-message e) (name key))))))

(deftest preset-test
  (testing "the acceptance condition: matern-32 with order 5 throws, with order 3 passes"
    (let [e (ex-of #(k/kernel :matern-32 {:order 5}))]
      (is (some? e))
      (is (= {:type :invalid-input :kind :kernel :name :matern-32 :option :order :value 5 :expected 3} (ex-data e)))
      (is (str/includes? (ex-message e) ":matern-32")))
    (is (fn? (k/kernel :matern-32 {:order 3})))
    (is (= ((k/kernel :matern-32 {:order 3}) [0.0] [1.0]) ((k/kernel :matern-32) [0.0] [1.0]))))
  (testing "the same value as a double or other keys pass; other keys are used"
    (is (= ((k/kernel :matern-32 {:order 3.0}) [0.0] [1.0]) ((k/kernel :matern-32) [0.0] [1.0])))
    (is (= ((k/kernel :matern-32 {:theta 2.0}) [0.0] [1.0]) ((k/kernel :matern {:order 3 :theta 2.0}) [0.0] [1.0])))
    (is (= ((k/kernel :triangular {:r 2.0}) [0.0] [1.0]) ((k/kernel :geometric {:n 1 :r 2.0}) [0.0] [1.0])))
    (is (= ((k/kernel :circular {:n 2}) [0.0] [0.5]) ((k/kernel :circular) [0.0] [0.5]))))
  (testing "a nil value of a preset key is a missing key"
    (is (= ((k/kernel :matern-32 {:order nil}) [0.0] [1.0]) ((k/kernel :matern-32) [0.0] [1.0]))))
  (testing "every vector preset: a conflicting value throws, the same value passes"
    (doseq [[kw preset] vector-presets
            [key value] preset
            :let [e (ex-of #(k/kernel kw {key (different value)}))]]
      (is (some? e) (str kw))
      (is (= key (:option (ex-data e))))
      (is (= (different value) (:value (ex-data e))))
      (is (= value (:expected (ex-data e))))
      (is (fn? (k/kernel kw {key value})))))
  (testing "every rbf preset: a conflicting value throws, the same value passes"
    (doseq [[kw preset] rbf-presets
            [key value] preset
            :let [e (ex-of #(k/rbf kw {key (different value)}))]]
      (is (some? e) (str kw " " key))
      (is (= key (:option (ex-data e))) (str kw))
      (is (= :rbf (:kind (ex-data e))))
      (is (= value (:expected (ex-data e))))
      (is (fn? (k/rbf kw {key value})) (str kw " " key))
      (is (fn? (k/rbf kw preset)) (str kw))
      (is (fn? (k/rbf kw (assoc preset :scale 2.0))) (str kw))))
  (testing "a preset with a user scale"
    (is (close? ((k/rbf :gaussians-laguerre-11 {:scale 2.0}) 1.0) ((k/rbf :gaussians-laguerre-11) 0.5) 1e-15))))

(deftest truncated-power-3-test
  (testing "the preset k is the number 3.0 (a vector before: ClassCastException)"
    (is (close? ((k/rbf :truncated-power-3) 0.5) 0.125 1e-15))
    (is (close? ((k/rbf :truncated-power-3) 0.0) 1.0 1e-15))
    (is (close? ((k/rbf :truncated-power-3 {:scale 2.0}) 1.0) 0.125 1e-15))
    (is (== 0.0 ((k/rbf :truncated-power-3) 1.5)))
    (is (close? ((k/rbf :truncated-power-3) 0.3) (Math/pow 0.7 3.0) 1e-15))
    (doseq [r [0.0 0.2 0.5 0.9 1.0 1.2]]
      (is (close? ((k/rbf :truncated-power-3) r) ((k/rbf :truncated-power {:k 3.0}) r) 1e-15))
      (is (close? ((k/rbf :truncated-power-2) r) ((k/rbf :truncated-power {:k 2.0}) r) 1e-15)))))

(def ^:private g (k/kernel :gaussian))
(def ^:private l (k/kernel :laplacian))
(def ^:private c (k/kernel :cauchy))
(def ^:private p (k/kernel :polynomial))

(def ^:private points
  [[[0.0] [1.0]] [[0.3 0.4] [1.2 -0.7]] [[1.0 2.0 3.0] [1.5 2.5 3.5]] [[0.0 0.0] [0.0 0.0]] [[5.0] [-3.0]]])

(deftest wmean-test
  (testing "weighted mean: (g + 3 l) / 4 for weights [1 3] (a function object before)"
    (doseq [[x y] points
            :let [expected (/ (+ (g x y) (* 3.0 (l x y))) 4.0)
                  w (k/wmean [g l] [1 3])]]
      (is (number? (w x y)))
      (is (close? (w x y) expected 1e-15) (str x y))))
  (testing "the weights are normalised: only the ratios matter"
    (doseq [[x y] points]
      (is (close? ((k/wmean [g l] [2 6]) x y) ((k/wmean [g l] [1 3]) x y) 1e-15))
      (is (close? ((k/wmean [g l] [0.25 0.75]) x y) ((k/wmean [g l] [1 3]) x y) 1e-15))))
  (testing "default weights give the arithmetic mean; one kernel"
    (doseq [[x y] points]
      (is (close? ((k/wmean [g l c]) x y) (/ (+ (g x y) (l x y) (c x y)) 3.0) 1e-15))
      (is (close? ((k/wmean [g]) x y) (g x y) 1e-15))
      (is (close? ((k/wmean [g l] [1 0]) x y) (g x y) 1e-15))))
  (testing "negative weights with a non-zero sum are a plain linear combination over their sum"
    (is (close? ((k/wmean [g l] [3 -1]) [0.0] [1.0]) (/ (- (* 3.0 (g [0.0] [1.0])) (l [0.0] [1.0])) 2.0) 1e-15)))
  (testing "wrong arguments throw ex-info"
    (doseq [f [#(k/wmean [] [])
               #(k/wmean [])
               #(k/wmean [g l] [1])
               #(k/wmean [g l] [1 2 3])
               #(k/wmean [g l] [1 -1])
               #(k/wmean [g l] [0 0])
               #(k/wmean [g l] [1 ##NaN])
               #(k/wmean [g l] [1 ##Inf])
               #(k/wmean [g 1.0] [1 1])]
            :let [e (ex-of f)]]
      (is (some? e))
      (is (= :invalid-input (:type (ex-data e)))))))

(deftest sum-and-product-test
  (testing "sum of 1 … 5 kernels, values exact at 1e-15"
    (doseq [[x y] points
            :let [kernels [g l c p g]]]
      (doseq [n (range 1 6)
              :let [ks (take n kernels)
                    vals (map #(% x y) ks)]]
        (is (close? ((apply k/sum ks) x y) (reduce + vals) 1e-15) (str n))
        (is (close? ((apply k/product ks) x y) (reduce * vals) 1e-15) (str n))
        (is (= ((apply k/product ks) x y) ((apply k/mult ks) x y)) (str n)))))
  (testing "hand values"
    (let [a (fn [_ _] 1.5) b (fn [_ _] -0.25) d (fn [_ _] 4.0)]
      (is (== 1.5 ((k/sum a) 0 0)))
      (is (== 1.25 ((k/sum a b) 0 0)))
      (is (== 5.25 ((k/sum a b d) 0 0)))
      (is (== 6.5 ((k/sum a b d a b) 0 0)))
      (is (== -0.375 ((k/product a b) 0 0)))
      (is (== -1.5 ((k/product a b d) 0 0)))
      (is (== 0.5625 ((k/product a b d a b) 0 0)))))
  (testing "one kernel is returned as it is"
    (is (identical? g (k/sum g)))
    (is (identical? g (k/product g)))
    (is (identical? g (k/mult g))))
  (testing "results are doubles"
    (is (double? ((k/sum g l) [0.0] [1.0])))
    (is (double? ((k/product g l c) [0.0] [1.0])))
    (is (double? ((k/sum g l c p) [0.0] [1.0]))))
  (testing "at least one kernel, only functions"
    (is (thrown? clojure.lang.ArityException (apply k/sum [])))
    (is (thrown? clojure.lang.ArityException (apply k/product [])))
    (doseq [f [#(k/sum g 1.0) #(k/sum 1.0) #(k/product g l 2.0) #(k/product nil) #(k/mult g "a")]
            :let [e (ex-of f)]]
      (is (some? e))
      (is (= :invalid-input (:type (ex-data e)))))))

(deftest exp-scale-test
  (testing "exp and scale"
    (doseq [[x y] points]
      (is (close? ((k/exp g) x y) (Math/exp (g x y)) 1e-15))
      (is (close? ((k/exp g 2.5) x y) (Math/exp (* 2.5 (g x y))) 1e-15))
      (is (close? ((k/exp p -0.1) x y) (Math/exp (* -0.1 (p x y))) 1e-15))
      (is (close? ((k/scale g 3.0) x y) (* 3.0 (g x y)) 1e-15))
      (is (close? ((k/scale l 0.0) x y) 0.0 1e-15))))
  (testing "combinations"
    (is (close? ((k/sum (k/scale g 2.0) (k/scale l 3.0)) [0.0] [1.0]) (+ (* 2.0 (g [0.0] [1.0])) (* 3.0 (l [0.0] [1.0]))) 1e-15))))

(deftest rbf->kernel-test
  (testing "re-exported from fastmath.kernel.vector"
    (is (identical? vk/rbf->kernel @#'k/rbf->kernel))
    (is (close? ((k/rbf->kernel (k/rbf :gaussian)) [0.0 0.0] [3.0 4.0]) (Math/exp -25.0) 1e-15))
    (is (close? ((k/rbf->kernel (k/rbf :gaussian) v/dist-abs) [0.0 0.0] [3.0 4.0]) (Math/exp -49.0) 1e-15))
    (is (close? ((k/rbf->kernel (k/rbf :truncated-power-3)) [0.0] [0.5]) 0.125 1e-15))
    (is (= '([rbf-kernel] [rbf-kernel distance]) (:arglists (meta #'k/rbf->kernel))))))

(deftest documentation-test
  (testing "kernel-list has its docstring on the var (on the map value before)"
    (is (some? (:doc (meta #'k/kernel-list))))
    (is (str/includes? (:doc (meta #'k/kernel-list)) ":vector"))
    (is (nil? (:doc (meta k/kernel-list))))
    (is (= #{:vector :rbf :kde :window} (set (keys k/kernel-list))))
    (is (every? seq (vals k/kernel-list)))
    (is (some? (:doc (meta #'k/rbf->kernel)))))
  (testing "every public var of fastmath.kernel.vector and the combinators has a docstring in the house style"
    (doseq [var (concat (vals (ns-publics 'fastmath.kernel.vector))
                        (map #(ns-resolve 'fastmath.kernel %) '[exp scale sum product mult wmean]))
            :let [doc (:doc (meta var))
                  nm (str (:name (meta var)))]]
      (is (some? doc) nm)
      (is (not (str/includes? doc "\"")) (str nm " double quote"))
      (is (not (str/includes? doc "\\")) (str nm " backslash"))
      (is (not (str/includes? doc "```")) (str nm " code block"))
      (is (not (re-find #"(?m)^#" doc)) (str nm " heading"))
      (is (str/includes? doc "Returns") (str nm " Returns"))
      (is (str/includes? doc "See also") (str nm " See also"))
      (is (str/includes? doc "Parameters") (str nm " Parameters"))))
  (testing "the parameter blocks are bullet lists"
    (doseq [var (vals (ns-publics 'fastmath.kernel.vector))
            :let [doc (:doc (meta var))]
            :when (str/includes? doc "Parameters:")
            :let [after (second (str/split doc #"Parameters:\n\n" 2))]]
      (is (str/starts-with? (or after "") "  - ") (str (:name (meta var)))))))

(deftest window-registry-test
  (let [names (set (:window k/kernel-list))]
    (testing "the registered windows: 61 of the earlier version without :low-sidelobe plus the 15 added in the audit"
      (is (= 76 (count names)))
      (is (not (contains? names :low-sidelobe)))
      (is (every? names [:cosine :kaiser :general-gaussian :flat-top-hp :flat-top-srs :flat-top-bk :hft70 :hft90d
                         :hft95 :hft116d :hft144d :hft169d :hft196d :hft223d :hft248d])))
    (testing "every registered window is a function of the same name in fastmath.kernel.window"
      (doseq [nm names
              :let [v (ns-resolve 'fastmath.kernel.window (symbol (name nm)))]]
        (is (some? v) (str nm))))
    (testing "every window: N = 0 gives (), N = 1 gives (1.0) (symmetric and periodic), N = 8 gives 8 finite numbers"
      (doseq [nm names
              symmetric? [true false]
              :let [opts {:symmetric? symmetric?}
                    default-1 (k/window nm 1 opts)]]
        (is (= () (k/window nm 0 opts)) (str nm))
        (is (= 1 (count default-1)) (str nm))
        (let [w (k/window nm 8 opts)]
          (is (= 8 (count w)) (str nm))
          (is (every? #(Double/isFinite (double %)) w) (str nm)))))
    (testing "the periodic window is the symmetric window of N + 1 samples without the last one"
      (doseq [nm [:hann :tukey :kaiser :cosine :hft95 :dolph-chebyshev]]
        (is (= (take 9 (k/window nm 10)) (k/window nm 9 {:symmetric? false})) (str nm))))
    (testing "nil :symmetric? is the default"
      (is (= (k/window :hann 7) (k/window :hann 7 {:symmetric? nil}))))))

(deftest window-validation-test
  (testing "unknown window: the message lists the valid keys"
    (let [e (try (k/window :low-sidelobe 8) (catch clojure.lang.ExceptionInfo e e))]
      (is (= :invalid-input (:type (ex-data e))))
      (is (= :low-sidelobe (:value (ex-data e))))
      (is (str/includes? (ex-message e) ":hann"))))
  (testing "N: an integer >= 0"
    (doseq [bad [-1 3.5 nil "5" ##NaN ##Inf]]
      (is (thrown? clojure.lang.ExceptionInfo (k/window :hann bad)) (pr-str bad))
      (is (thrown? clojure.lang.ExceptionInfo (k/window (fn [n _] (repeat n 1.0)) bad)) (pr-str bad)))
    (is (= 5 (count (k/window :hann 5.0)))))
  (testing "a custom window function must return N coefficients, in both modes"
    (is (= [1 2 3 4] (k/window (fn [n _] (range 1 (inc n))) 4)))
    (is (thrown? clojure.lang.ExceptionInfo (k/window (fn [_ _] [1.0 2.0 3.0]) 4)))
    (is (thrown? clojure.lang.ExceptionInfo (k/window (fn [_ _] [1.0 2.0 3.0]) 4 {:symmetric? false})))
    (is (= (range 1 5) (k/window (fn [n _] (range 1 (inc n))) 4 {:symmetric? false})))
    (is (thrown? clojure.lang.ExceptionInfo (k/window (fn [_ _] 5) 4))))
  (testing "a custom window function may return an array"
    (let [f (fn [n _] (double-array (range 1 (inc n))))]
      (is (= 4 (alength ^doubles (k/window f 4))))
      (is (= [1.0 2.0 3.0] (vec (k/window f 3 {:symmetric? false}))))
      (is (thrown? clojure.lang.ExceptionInfo (k/window (fn [_ _] (double-array 3)) 4))))))
