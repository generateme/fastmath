(ns fastmath.signal-test
  (:require [fastmath.signal :as sut]
            [fastmath.signal.iir :as iir]
            [fastmath.transform :as trans]
            [fastmath.kernel :as k]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [clojure.java.io :as io]
            [clojure.edn :as edn]
            [clojure.test :as t]))

;; Reference data: `utils/fastmath/dev/generate_signal_reference.py` (scipy 1.18.1, numpy 2.5.2),
;; Python ports of the original effect sources and scipy spectra, filters and resampling.

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "signal/signal_reference.edn")))))

(defn- scaled-error
  "`max |a - b| / max(1, max |b|)`, infinite when lengths differ."
  ^double [a b]
  (if-not (== (count a) (count b))
    ##Inf
    (/ (double (reduce max 0.0 (map (fn [^double p ^double q] (m/abs (- p q))) a b)))
       (max 1.0 (double (reduce max 0.0 (map (fn [^double q] (m/abs q)) b)))))))

(defn- throws?
  [f]
  (try (f) false (catch clojure.lang.ExceptionInfo _ true)))

;; effects

(defn- run-effect
  [name params input]
  (let [e (if (= name :gain) (sut/effect :gain (:gain params)) (sut/effect name params))]
    (vec (sut/apply-effects-raw input e))))

(t/deftest effects-vs-original-sources
  (let [{:keys [x effects]} @reference
        x (double-array x)
        impulse (double-array (cons 1.0 (repeat (dec (alength x)) 0.0)))]
    (doseq [[name params input expected] effects]
      (t/testing (str name " " params)
        ;; 1e-12 (spec WP13); measured max 1.2e-13 (`:dj-eq`), most effects are exact
        (t/is (< (scaled-error (run-effect name params (if (= input :x) x impulse)) expected) 1.0e-12))))))

(t/deftest effects-construction
  (t/testing "every effect except :iir is constructible with defaults and processes a signal"
    (let [x (double-array (map (fn [^long i] (m/sin (* 0.1 i))) (range 200)))]
      (doseq [e (remove #{:iir} sut/effects-list)]
        (t/is (every? m/valid-double? (sut/apply-effects-raw x (sut/effect e))) (str e)))))
  (t/testing ":iir"
    (t/is (throws? #(sut/effect :iir)))
    (t/is (throws? #(sut/effect :iir {})))
    (let [f (iir/butterworth :lowpass 3 {:cutoff 0.1})
          xs (map (fn [^long i] (m/sin (* 0.3 i))) (range 50))
          by-effect (sut/apply-effects xs (sut/effect :iir f))]
      (sut/reset-IIR! f)
      (t/is (= (vec by-effect) (sut/filter-signal f xs)))))
  (t/testing ":gain accepts a number and a map"
    (let [x (double-array [1.0 -2.0 0.5])]
      (t/is (= (vec (sut/apply-effects-raw x (sut/effect :gain 0.5)))
               (vec (sut/apply-effects-raw x (sut/effect :gain {:gain 0.5})))
               [0.5 -1.0 0.25]))
      (t/is (= [1.0 -2.0 0.5] (vec (sut/apply-effects-raw x (sut/effect :gain)))))))
  (t/testing "invalid parameters"
    (t/is (throws? #(sut/effect :echo {:delay 0.0})))
    (t/is (throws? #(sut/effect :echo {:delay 1.0e-5})))
    (t/is (sut/effect :echo {:delay (/ 44100.0)}))
    (t/is (throws? #(sut/effect :svf1 {:Q 1.07})))
    (t/is (sut/effect :svf1 {:Q 1.065}))
    (t/is (throws? #(sut/effect :svf1 {:oversamples 0})))
    (t/is (throws? #(sut/effect :svf1 {:kind :peak})))
    (t/is (throws? #(sut/effect :svf2 {:kind :allpass})))
    (t/is (throws? #(sut/effect :svf2 {:tan :approx})))
    (t/is (throws? #(sut/effect :svf3 {:kind :band})))
    (t/is (throws? #(sut/effect :clipping {:method :cubic})))
    (t/is (throws? #(sut/effect :divider {:denom 0})))
    (t/is (throws? #(sut/effect :divider {:denom 1.5})) "not an integer")
    (t/is (sut/effect :divider {:denom 3.0}))
    (t/is (throws? #(sut/effect :echo {:delay 1000.0})) "more than 1e7 samples")
    (t/is (sut/effect :echo {:delay (/ 1.0e7 44100.0)}) "exactly 1e7 samples")
    (doseq [e [:svf2 :svf3] q [0 0.0 -1.0]]
      (t/is (throws? #(sut/effect e {:Q q})) (str e " Q " q)))
    (t/is (throws? #(sut/effect :biquad-hs {:gain 36.0})) "shelf slope 1.5 too large for 36 dB")
    (t/is (throws? #(sut/effect :biquad-lp {:fc 30000.0})))))

(t/deftest compose-order
  (t/testing "the last effect is applied first"
    (let [e (sut/compose-effects (sut/effect :gain 10.0) (sut/effect :clipping))]
      (t/is (= [10.0] (vec (sut/apply-effects-raw (double-array [1.5]) e)))))))

(t/deftest effects-boundaries
  (let [x (double-array [0.3 -0.7 0.0 0.9 -0.2])]
    (t/testing "decimator without bit reduction and rate reduction is identity"
      (t/is (= (vec x) (vec (sut/apply-effects-raw x (sut/effect :decimator {:bits 31.0 :fs 44100.0}))))))
    (t/testing "decimator of zero gives -step/2 as the original"
      (t/is (m/delta-eq (* -0.5 (m/pow 0.5 (- 4.0 0.999)))
                        (first (sut/apply-effects-raw (double-array 1) (sut/effect :decimator {:bits 4.0 :fs 44100.0})))
                        1.0e-15)))
    (t/testing "one sample echo"
      (t/is (v/edelta-eq [0.3 -0.55 -0.275 0.7625 0.18125]
                         (vec (sut/apply-effects-raw x (sut/effect :echo {:delay (/ 44100.0)})))
                         1.0e-15))))
  (t/testing "reset each 2 samples"
    (let [x (double-array [1.0 0.0 1.0 0.0])
          e (sut/effect :simple-lowpass)
          res (sut/apply-effects-raw x e 2)]
      (t/is (= (aget res 0) (aget res 2)))
      (t/is (= (aget res 1) (aget res 3))))))

;; spectra

(def ^:private spectrum-methods [:psd :power :energy :magnitude :amplitude :asd :phase :real :imag])

(t/deftest spectrum-vs-scipy
  (doseq [{:keys [x one-sided? freqs] :as e} (:spectra @reference)
          :let [opts {:fs 2.5 :spectrum (if one-sided? :single-sided :double-sided)}
                txs (trans/fft x (select-keys opts [:spectrum]))]]
    (t/testing (str "N=" (count x) " one-sided? " one-sided?)
      ;; measured max 6.5e-16
      (doseq [method spectrum-methods]
        (t/is (< (scaled-error (:spectrum (sut/spectrum x (assoc opts :method method))) (e method)) 1.0e-14) (str method)))
      (t/is (< (scaled-error (:freqs (sut/spectrum x opts)) freqs) 1.0e-14))
      (t/is (== (count x) (:N (sut/spectrum x opts))))
      (t/testing "normalization of fft does not change values"
        (t/is (< (scaled-error (:spectrum (sut/spectrum x (assoc opts :method :psd :norm :ortho))) (:psd e)) 1.0e-14))
        (t/is (< (scaled-error (:spectrum (sut/spectrum x (assoc opts :method :psd :norm :forward))) (:psd e)) 1.0e-14)))
      (t/testing "frequency domain, metadata or explicit layout"
        (t/is (< (scaled-error (:spectrum (sut/spectrum txs {:domain :frequency :method :psd :fs 2.5})) (:psd e)) 1.0e-14))
        (t/is (< (scaled-error (:spectrum (sut/spectrum (with-meta txs nil)
                                                        {:domain :frequency :method :psd :fs 2.5
                                                         :kind (if one-sided? :real :complex) :n (count x)}))
                               (:psd e)) 1.0e-14))
        (when one-sided?
          (t/is (< (scaled-error (:spectrum (sut/spectrum (with-meta txs nil)
                                                          {:domain :frequency :method :psd :fs 2.5
                                                           :nyquist? (even? (count x))}))
                                 (:psd e)) 1.0e-14))))
      (t/testing "fft-* functions"
        (t/is (< (scaled-error (sut/fft-energy txs) (:energy e)) 1.0e-14))
        (t/is (< (scaled-error (sut/fft-magnitude txs) (:magnitude e)) 1.0e-14))
        (t/is (< (scaled-error (sut/fft-amplitude txs) (:amplitude e)) 1.0e-14))
        (t/is (< (scaled-error (sut/fft-power txs) (:power e)) 1.0e-14))))))

(t/deftest spectrum-boundaries
  (t/testing "DC of a constant signal (S-2)"
    (t/is (= [1.0 0.0 0.0] (:spectrum (sut/spectrum [1 1 1 1]))))
    (t/is (= [4.0 0.0 0.0] (:spectrum (sut/spectrum [1 1 1 1] {:method :magnitude})))))
  (t/testing "one sample"
    (t/is (= {:N 1 :spectrum [2.0] :freqs [0.0]} (sut/spectrum [2.0])))
    (t/is (= [4.0] (:spectrum (sut/spectrum [2.0] {:method :psd})))))
  (t/testing "sine amplitude, one-sided"
    (let [xs (map (fn [^long i] (* 3.0 (m/sin (* m/TWO_PI 2.0 (/ i 16.0))))) (range 16))]
      (t/is (m/delta-eq 3.0 ((:spectrum (sut/spectrum xs)) 2) 1.0e-14))))
  (t/testing "decibels are not floored: 10 log10 for power, 20 log10 for amplitude, -Infinity for 0"
    (t/is (= [-120.0 ##-Inf] (:spectrum (sut/spectrum [1e-6 1e-6] {:method :power :db? true}))))
    (t/is (= [-120.0 ##-Inf] (:spectrum (sut/spectrum [1e-6 1e-6] {:method :amplitude :db? true})))))
  (t/testing "decibels: 10 log10 for power, 20 log10 for amplitude"
    (t/is (m/delta-eq 20.0 (first (:spectrum (sut/spectrum [10 10] {:method :amplitude :db? true}))) 1.0e-14))
    (t/is (m/delta-eq 20.0 (first (:spectrum (sut/spectrum [10 10] {:method :power :db? true}))) 1.0e-14)))
  (t/testing "errors"
    (t/is (throws? #(sut/spectrum [])))
    (t/is (throws? #(sut/spectrum [1 2 3] {:method :db})))
    (t/is (throws? #(sut/spectrum [1 2 3] {:domain :time-domain})))
    (t/is (throws? #(sut/fft-energy (trans/fft [1 2 3 4]) {:n 7})))
    (t/is (throws? #(sut/fft-energy [] {:kind :complex})))
    (t/is (throws? #(sut/fft-energy (trans/fft [1 2 3 4]) {:kind :half})))
    (t/is (throws? #(sut/spectrum (map identity (trans/fft [1 2 3])) {:domain :frequency})) "one-sided without metadata, :n and :nyquist?")
    (t/is (throws? #(sut/fft-power (vec (trans/fft [1 2 3 4])))))
    (t/is (== 3 (:N (sut/spectrum (map identity (trans/fft [1 2 3])) {:domain :frequency :nyquist? false}))))
    (t/is (== 3 (:N (sut/spectrum (map identity (trans/fft [1 2 3])) {:domain :frequency :n 3})))))
  (t/testing "fft-frequencies"
    (t/is (= [0.0 2.0 4.0 6.0] (sut/fft-frequencies 8.0 4)))
    (t/is (= [0.0] (sut/fft-frequencies 8.0 1)))))

(t/deftest periodogram
  (let [sig [1,2,3,10,-1,-2,3,1]]
    (t/testing "Welch periodogram (comparison with scipy)"
      ;; >>> scipy.signal.welch(x=np.asarray([1,2,3,10,-1,-2,3,1]),window=[0,0.75,0.75,0],detrend=False,noverlap=2,scaling="density")
      ;; (array([0.  , 0.25, 0.5 ]), array([17.83333333, 42.33333333, 24.5       ]))
      (t/is (v/edelta-eq [17.83333333, 42.33333333, 24.5]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]}))))
      (t/is (v/edelta-eq [8.91666667, 21.16666667, 12.25]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :method :power}))))
      (t/is (v/edelta-eq [15.0 , 15.6, 15.0]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :average :umedian}))))
      (t/is (v/edelta-eq [7.5, 7.8, 7.5]
                         (:spectrum (sut/periodogram sig {:window [0 0.75 0.75 0]
                                                          :method :power
                                                          :average :umedian}))))))
  (t/testing "6 and 7 segments, mean and median, density and spectrum"
    (doseq [{:keys [x window freqs] :as e} (:welch @reference)
            [average prefix] [[:mean "mean"] [:umedian "median"]]
            [method scaling] [[:psd "density"] [:power "spectrum"]]
            :let [res (sut/periodogram x {:window window :fs 3.0 :method method :average average})]]
      ;; measured max 7.7e-16
      (t/is (< (scaled-error (:spectrum res) (e (keyword (str prefix "-" scaling)))) 1.0e-14)
            (str (count x) " " average " " method))
      (t/is (< (scaled-error (:freqs res) freqs) 1.0e-14))))
  (t/testing "one segment: median = mean"
    (let [xs [1.0 2.0 -1.0 0.5]]
      (t/is (= (:spectrum (sut/periodogram xs {:window [1 1 1 1]}))
               (:spectrum (sut/periodogram xs {:window [1 1 1 1] :average :umedian}))))))
  (t/testing "errors"
    (t/is (throws? #(sut/periodogram [1 2 3] {:window [1 1 1 1]})))
    (t/is (throws? #(sut/periodogram [1 2 3 4] {:window [1 1] :method :complex})))
    (t/is (throws? #(sut/periodogram [1 2 3 4] {:window [1 1] :average :mode})))))

(t/deftest stft
  (t/testing "odd window vs scipy.signal.spectrogram"
    (let [{:keys [x window psd power freqs times]} (:spectrogram @reference)
          opts {:window window :overlap (/ 2.0 3.0) :fs 2.0}
          r-psd (sut/stft x (assoc opts :method :psd))
          r-power (sut/stft x (assoc opts :method :power))]
      ;; measured max 5.2e-16
      (t/is (< (scaled-error (flatten (:spectrum r-psd)) (flatten psd)) 1.0e-14))
      (t/is (< (scaled-error (flatten (:spectrum r-power)) (flatten power)) 1.0e-14))
      (t/is (= (count psd) (count (:spectrum r-psd))))
      (t/is (< (scaled-error (:freqs r-psd) freqs) 1.0e-14))
      (t/is (< (scaled-error (:times r-psd) times) 1.0e-14))))
  (t/testing "signal of window length gives one segment, partial segment dropped"
    (t/is (= 1 (count (:spectrum (sut/stft [1 2 3] {:window [1 1 1]})))))
    (t/is (= 2 (count (:spectrum (sut/stft [1 2 3 4 5 6 7] {:window [1 1 1] :overlap 0.0})))))
    (t/is (= [0.75] (:times (sut/stft [1 2 3] {:window [1 1 1] :fs 2.0})))))
  (t/testing "errors"
    (t/is (throws? #(sut/stft [1 2] {:window [1 1 1]})))
    (t/is (throws? #(sut/stft [1 2] {:window []})))
    (t/is (throws? #(sut/stft [1 2 3] {:window [1 1] :method :db})))
    (t/is (throws? #(sut/stft [])))
    (doseq [o [1.0 1.5 -0.1 nil]]
      (t/is (throws? #(sut/stft [1 2 3 4] {:window [1 1] :overlap o})) (str "overlap " o))))
  (t/testing "default window: Gaussian of length min(N, 256)"
    (t/is (= 3 (count (:freqs (sut/stft (range 4))))) "N=4, window 4")
    (t/is (= 129 (count (:freqs (sut/stft (range 1000))))) "window 256")
    (t/is (= [[2.0]] (:spectrum (sut/stft [2.0] {:method :amplitude}))) "one sample, window [1.0]"))
  (t/testing "default window is (k/window :gaussian (min N 256)): the Gaussian with the default :alpha 1.0, unchanged by the window audit"
    (let [x (map #(Math/sin (* 0.7 %)) (range 40))
          explicit (fn [n] (sut/stft x {:window (vec (k/window :gaussian n))}))]
      (t/is (= (explicit 40) (sut/stft x)))
      (t/is (= (sut/stft (range 300) {:window (vec (k/window :gaussian 256))}) (sut/stft (range 300)))
            "window of 256 samples for N > 256"))))

;; convolution and filtering

(t/deftest convolution
  (t/testing "full mode"
    (t/are [in1 in2 res] (and (v/edelta-eq res (sut/convolve in1 in2))
                              (v/edelta-eq res (sut/fft-convolve in1 in2)))
      [-1] [2] [-2]
      [-1 1] [2] [-2 2]
      [2] [-1 1] [-2 2]
      [-1 1] [2 3 4 5 6 7] [-2 -1 -1 -1 -1 -1 7]))
  (t/testing "all modes vs scipy"
    (let [{:keys [a b] :as e} (:convolve @reference)]
      (doseq [mode [:full :same :valid]
              [nm direct fft] [["convolve" sut/convolve sut/fft-convolve]
                               ["correlate" sut/correlate sut/fft-correlate]]
              [order p q] [["ab" a b] ["ba" b a]]
              :let [expected (e (keyword (str nm "-" order "-" (name mode))))]]
        ;; measured max 2.5e-16
        (t/is (< (scaled-error (direct p q mode) expected) 1.0e-14) (str nm " " order " " mode))
        (t/is (< (scaled-error (fft p q mode) expected) 1.0e-14) (str "fft " nm " " order " " mode)))))
  (t/testing "aliases, vectors, errors"
    (t/is (= (sut/convolve [1 2 3] [1 1] :full) (sut/convolve [1 2 3] [1 1] :all) [1.0 3.0 5.0 3.0]))
    (t/is (= (sut/convolve [1 2 3] [1 1] :same) (sut/convolve [1 2 3] [1 1] :first) [1.0 3.0 5.0]))
    (t/is (vector? (sut/fft-convolve [1 2 3] [1 1] :valid)))
    (t/is (throws? #(sut/convolve [] [1])))
    (t/is (throws? #(sut/fft-convolve [1] [])))
    (t/is (throws? #(sut/convolve [1 2] [1] :middle)))))

(t/deftest correlation
  (t/testing "full mode"
    (t/are [in1 in2 res] (and (v/edelta-eq res (sut/correlate in1 in2))
                              (v/edelta-eq res (sut/fft-correlate in1 in2)))
      [-1] [2] [-2]
      [-1 1] [2] [-2 2]
      [2] [-1 1] [2 -2]
      [-1 1] [2 3 4 5 6 7] [-7 1 1 1 1 1 2])))

(t/deftest filtering
  (t/testing "FIR is centred convolution, direct and fft give the same result"
    (let [fir [0.25 0.5 0.25]
          short (vec (range 10))
          long (mapv (fn [^long i] (m/sin (* 0.01 i))) (range 3000))]
      (t/is (= (sut/convolve short fir :same) (sut/filter-signal fir short)))
      (t/is (< (scaled-error (sut/filter-signal fir long) (sut/convolve long fir :same)) 1.0e-14))
      (t/is (vector? (sut/filter-filter-signal fir short)))
      (t/is (= 10 (count (sut/filter-filter-signal fir short))))))
  (t/testing "IIR keeps state until reset"
    (let [f (iir/butterworth :lowpass 2 {:cutoff 0.2})
          xs [1.0 0.0 0.0 0.0]
          a (sut/filter-signal f xs)
          b (sut/filter-signal f xs)]
      (t/is (vector? a))
      (t/is (not= a b))
      (sut/reset-IIR! f)
      (t/is (= a (sut/filter-signal f xs))))))

;; smoothing

(t/deftest savgol
  (let [{:keys [x default] :as e} (:savgol @reference)]
    (t/testing "vs scipy.signal.savgol_filter(mode='constant')"
      ;; measured max 2.7e-15 (third derivative)
      (t/is (< (scaled-error ((sut/savgol-filter) x) default) 1.0e-14))
      (doseq [d (range 4)]
        (t/is (< (scaled-error ((sut/savgol-filter 7 3 d 0.5) x) (e (keyword (str "deriv-" d)))) 1.0e-14) (str d)))))
  (t/testing "boundaries"
    (t/is (= [1.0 2.0 3.0] ((sut/savgol-filter 1 0) [1 2 3])))
    (t/is (v/edelta-eq [1.2 2.0 2.0 1.8 1.4] ((sut/moving-average-filter 5) [1 2 3 4 0]) 1.0e-15))
    (t/is (throws? #(sut/savgol-filter 4)))
    (t/is (throws? #(sut/savgol-filter 5 5)))
    (t/is (throws? #(sut/savgol-filter 5 2 3)))
    (t/is (throws? #(sut/savgol-filter 5 2 -1)))
    (t/is (throws? #(sut/savgol-filter 5 2 1 0.0)))
    (t/is (= [1.0 2.0 3.0] ((sut/moving-average-filter 1) [1 2 3])) "length 1 is identity")
    (t/is (v/edelta-eq [2.0 3.0 2.0] ((sut/moving-average-filter 3) [3 3 3]) 1.0e-15))))

(t/deftest kernel-smoothing
  (let [f (sut/kernel-smoothing-filter (k/kernel :gaussian) 5)
        res (f (repeat 9 2.0))]
    (t/is (vector? res))
    (t/is (m/delta-eq 2.0 (res 4) 1.0e-15))
    (t/is (< (res 0) 2.0) "zero padding at edges")
    (t/is (throws? #(sut/kernel-smoothing-filter (k/kernel :gaussian) 4)))))

;; transforms

(t/deftest resampling
  (t/testing "even length signal"
    (let [sig [1 2 3 10 -1 2]]
      (t/are [cnt res] (v/edelta-eq res (sut/resample sig cnt))
        1 [2.83333333]
        2 [0.16666667 5.5]
        3 [0.16666667 5.16666667 3.16666667]
        4 [2.83333333 1.32136721 8.16666667 -0.98803387]
        5 [2.83333333 0.27137889 7.59165872 4.03785588 -0.56756016]
        6 sig
        7 [1.0 2.00611241 1.50782978 8.71029852 5.90272448 -1.74568137 2.45204951]
        8 [1.0 1.90587372 1.32136721 5.3937861 10.0 1.45139187 -0.98803387 2.58228164]
        9 [1.0 1.77534578 1.5233184  3.0 9.43579021 7.16147413 -1.0 0.03886401 2.56520747]
        10 [1.0 1.64706163 1.75457671 1.7937036 7.02512757 10. 3.47132473 -1.76009924 0.91563766 2.48600068])))
  
  (t/testing "odd length signal"
    (let [sig [1 2 3 10 -1]]
      (t/are [cnt res] (v/edelta-eq res (sut/resample sig cnt))
        1 [3.0]
        2 [-0.68328157 6.68328157]
        3 [-0.68328157 4.40470422 5.27857735]
        4 [1.0 0.81218754 8.36656315 1.82124931]
        5 sig
        6 [1.0 2.7968157 0.64602959 8.36656315 7.35397041 -2.16337885]
        7 [1.0 3.2183382 0.34969088 4.51567976 10.22038208 4.25635504 -2.56044598]
        8 [1.0 3.40706996 0.81218754 1.87941562 8.36656315 9.32953114 1.82124931 -2.61601671]
        9 [1.0 3.46356369 1.4337983 0.64602959 5.41295943 10.08828259 7.35397041 0.12347689 -2.5220809 ]
        10 [1.0 3.4472136 2.0 0.31671843 3.0 8.36656315 10.0 5.23606798 -1.0 -2.36656315])))
  (t/testing "vs scipy.signal.resample"
    (doseq [{:keys [x n y]} (:resample @reference)]
      ;; measured max 8.9e-16
      (t/is (< (scaled-error (sut/resample x n) y) 1.0e-14) (str (count x) "->" n))))
  (t/testing "boundaries"
    (t/is (= [1.0 2.0] (sut/resample [1 2] 2)))
    (t/is (= [5.0 5.0 5.0] (sut/resample [5] 3)))
    (t/is (throws? #(sut/resample [] 3)))
    (t/is (throws? #(sut/resample [1 2] 0)))))

(t/deftest hilbert
  (t/testing "if real values in hilber transform are same as signal"
    (let [data (map #(repeatedly % r/grand) (range 1 101))]
      (doseq [d data]
        (t/is (v/delta-eq d (map first (sut/hilbert d)))))))
  (t/testing "vs scipy.signal.hilbert"
    (doseq [{:keys [x real imag]} (:hilbert @reference)
            :let [h (sut/hilbert x)]]
      ;; measured max 5.6e-16
      (t/is (< (scaled-error (map first h) real) 1.0e-14))
      (t/is (< (scaled-error (map second h) imag) 1.0e-14))))
  (t/testing "instantaneous frequency of a sine"
    (let [xs (map (fn [^long i] (m/cos (* m/TWO_PI 0.125 i))) (range 64))
          f (sut/instantaneous-frequency xs 2.0)]
      (t/is (= 63 (count f)))
      ;; periodic sine at an exact bin, analytic signal is exact: measured max 2.2e-16 * fs
      (t/is (every? #(m/delta-eq 0.25 % 1.0e-12) f)))))

(t/deftest sample-waveform
  (let [w (sut/waveform :sine {:f 1.0})]
    (t/is (= 4 (count (sut/sample-waveform w 4.0))))
    (t/is (v/edelta-eq [0.0 1.0 0.0 -1.0] (sut/sample-waveform w 4.0) 1.0e-15))
    (t/is (= 10 (count (sut/sample-waveform w 4.0 2.5))))
    (t/is (= [] (sut/sample-waveform w 4.0 0.0)))
    (t/is (= [0.0 1.0 2.0] (sut/sample-waveform identity 1.0 2.4)) "t < time")
    (t/is (= 3 (count (sut/sample-waveform identity 10.0 0.3))) "fs time = 3.0000000000000004")
    (t/is (= 1 (count (sut/sample-waveform identity 4.0 0.1))) "partial step")
    (t/is (throws? #(sut/sample-waveform w 0.0)))))

(t/deftest file-roundtrip
  (let [f (java.io.File/createTempFile "fastmath-signal" ".raw")]
    (try
      (sut/save-signal [0.0 1.0 -1.0 2.0 0.5] f)
      (let [res (vec (sut/load-signal f))]
        (t/is (= 5 (count res)))
        ;; 16 bit quantization: step 2/65535
        (t/is (v/edelta-eq [0.0 1.0 -1.0 1.0 0.5] res (/ 2.0 65535.0))))
      (finally (.delete f)))))

(t/deftest oscillators-removed
  (t/is (not-any? (ns-publics 'fastmath.signal)
                  '[oscillator oscillators oscillators-sum oscillator-gain oscillator->signal signal->oscillator])))
