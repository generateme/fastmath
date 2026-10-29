(ns fastmath.signal.pad
  "Signal padding (boundary extension).

  Every padding function takes a signal as a `double-array`, a target length `N` and a `side`, and returns a new `double-array` of length `N`. The original signal is placed according to `side`:

  * `:left` - padding added before the signal (signal ends at the last position)
  * `:right` - padding added after the signal (signal starts at position 0)
  * `:both` - padding split between both ends; when the padding length is odd, the extra sample goes to the left

  Available methods (names of the corresponding PyWavelets / NumPy modes in brackets):

  * [[zero]] - zeros (`zero`, `constant` with 0)
  * [[edge]] - repeated edge value (`constant` in PyWavelets, `edge` in NumPy)
  * [[periodic]] - periodic copies (`periodic`, `wrap`)
  * [[symmetric]] - mirrored copies including the edge sample (`symmetric`)
  * [[reflect]] - mirrored copies excluding the edge sample (`reflect`)
  * [[antisymmetric]] - mirrored and negated copies (`antisymmetric`)
  * [[antireflect]] - point reflection around the edge sample (`antireflect`, NumPy `reflect` with `reflect_type='odd'`)
  * [[linear]] - linear extrapolation of the two edge samples (`smooth`)

  The padding can be longer than the signal, in which case the extension pattern repeats. A one-sample signal is padded with copies of that sample in every method except [[zero]] and [[antisymmetric]] (alternating sign), as in PyWavelets.

  Every method throws `ex-info` when the signal is empty, when `N` is smaller than the signal length or when `side` is unknown."
  (:require [fastmath.core :as m])
  (:import [java.util Arrays]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn pad-position
  "Returns the index at which the original signal starts in the padded one.

  * `side` - `:left`, `:right` or `:both` (for `:both` an odd padding length puts the extra sample on the left)
  * `N` - target (padded) length
  * `len` - signal length

  Throws `ex-info` for an unknown `side`."
  ^long [side ^long N ^long len]
  (case side
    :left (m/- N len)
    :right 0
    :both (m/round (m// (m/- N len) 2.0))
    (throw (ex-info "Unknown padding side, use :left, :right or :both." {:side side}))))

(defn out-len-pos
  "Returns a zero-padded copy of `signal` together with its placement: `[out len start end]`.

  * `out` - new `double-array` of length `N` with `signal` copied into it, zeros elsewhere
  * `len` - signal length
  * `start` - index of the first signal sample in `out` (see [[pad-position]])
  * `end` - index just after the last signal sample

  This is the building block of all padding methods. Throws `ex-info` when `signal` is empty, `N` is smaller than the signal length or `side` is unknown."
  [^doubles signal ^long N side]
  (let [len (alength signal)]
    (when (m/zero? len)
      (throw (ex-info "Signal to pad is empty." {:length len})))
    (when (m/< N len)
      (throw (ex-info "Target length is smaller than the signal length." {:N N :length len})))
    (let [out (double-array N)
          start (pad-position side N len)
          end (m/+ start len)]
      (System/arraycopy signal 0 out start len)
      [out len start end])))

(defn zero
  "Returns zero-padded signal.

  0 ... 0 | s1 ... sn | 0 ... 0"
  [^doubles signal ^long N side]
  (first (out-len-pos signal N side)))

(defn edge
  "Returns signal padded with edge values.

  s1 ... s1 | s1 ... sn | sn ... sn"
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]
    (when (m/< end N) (Arrays/fill out end N (Array/aget signal (m/dec len))))
    (when-not (m/zero? start) (Arrays/fill out 0 start (Array/aget signal 0)))
    out))

(defn periodic
  "Returns signal padded with copies of the original signal.

  ... sn s1 ... sn | s1 ... sn | s1 ... sn s1 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]
    (when (m/< end N) (let [diff (m/long-sub N end)
                            copies (m/long-quot diff len)
                            remainder (m/long-mod diff len)]
                        (dotimes [id copies]
                          (System/arraycopy signal 0 out (m/+ end (m/* id len)) len))
                        (System/arraycopy signal 0 out (m/+ end (m/* copies len)) remainder)))
    (when-not (m/zero? start) (let [copies (m/long-quot start len)
                                    remainder (m/long-mod start len)]
                                (dotimes [id copies]                                  
                                  (System/arraycopy signal 0 out (m/- start (m/* (m/inc id) len)) len))
                                (System/arraycopy signal (m/- len remainder) out 0 remainder)))
    out))

(defn- reverse-signal
  [^doubles signal]
  (let [len (alength signal)
        nsignal (double-array len)]
    (dotimes [id len]
      (Array/aset nsignal id (Array/aget signal (m/- len id 1))))
    nsignal))

(defn symmetric
  "Returns signal padded with mirrored copies of the original signal.

  ... sn sn ... s1 | s1 ... sn | sn ... s1 s1 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
        ^doubles mirror (reverse-signal signal)]

    ;; right end
    (when (m/< end N) (let [diff (m/long-sub N end)
                            copies (m/long-quot diff len)
                            remainder (m/long-mod diff len)]
                        (dotimes [id copies]
                          (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/+ end (m/* id len)) len))
                        (System/arraycopy (if (m/even? copies) mirror signal) 0 out (m/+ end (m/* copies len)) remainder)))

    ;; left end
    (when-not (m/zero? start) (let [copies (m/long-quot start len)
                                    remainder (m/long-mod start len)]
                                (dotimes [id copies]
                                  (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (m/inc id) len)) len))
                                (System/arraycopy (if (m/even? copies) mirror signal) (m/- len remainder) out 0 remainder)))
    out))

(defn reflect
  "Returns signal padded with reflected copies of the original signal, the edge sample is not repeated.

  ... s3 s2 | s1 s2 ... sn-1 sn | sn-1 sn-2 ...

  A one-sample signal is padded with its value."
  [^doubles signal ^long N side]
  (if (m/one? (alength signal))
    (edge signal N side)
    (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
          ^doubles mirror (reverse-signal signal)
          len- (m/long-dec len)]

      ;; right end
      (when (m/< end N) (let [diff (m/long-sub N end)
                              copies (m/long-quot diff len-)
                              remainder (m/long-mod diff len-)]
                          (dotimes [id copies]
                            (System/arraycopy (if (m/even? id) mirror signal) 1 out (m/+ end (m/* id len-)) len-))
                          (System/arraycopy (if (m/even? copies) mirror signal) 1 out (m/+ end (m/* copies len-)) remainder)))

      ;; left end
      (when-not (m/zero? start) (let [copies (m/long-quot start len-)
                                      remainder (m/long-mod start len-)]
                                  (dotimes [id copies]
                                    (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (m/inc id) len-)) len-))
                                  (System/arraycopy (if (m/even? copies) mirror signal) (m/- len- remainder) out 0 remainder)))
      out)))

(defn- negate-signal
  [^doubles signal]
  (let [len (alength signal)
        nsignal (double-array len)]
    (dotimes [id len]
      (Array/aset nsignal id (m/- (Array/aget signal id))))
    nsignal))

(defn antisymmetric
  "Returns signal padded with mirrored and negated copies of the original signal.

  ... sn -sn ... -s1 | s1 ... sn | -sn ... -s1 s1 ...."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
        ^doubles mirror (negate-signal (reverse-signal signal))]

    ;; right end
    (when (m/< end N) (let [diff (m/long-sub N end)
                            copies (m/long-quot diff len)
                            remainder (m/long-mod diff len)]
                        (dotimes [id copies]
                          (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/+ end (m/* id len)) len))
                        (System/arraycopy (if (m/even? copies) mirror signal) 0 out (m/+ end (m/* copies len)) remainder)))

    ;; left end
    (when-not (m/zero? start) (let [copies (m/long-quot start len)
                                    remainder (m/long-mod start len)]
                                (dotimes [id copies]
                                  (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (m/inc id) len)) len))
                                (System/arraycopy (if (m/even? copies) mirror signal) (m/- len remainder) out 0 remainder)))
    
    out))

(defn linear
  "Returns signal padded with a linear extrapolation of the two edge samples on each side.

  ... s1-2d1 s1-d1 | s1 s2 ... sn-1 sn | sn+dn sn+2dn ...

  where `d1 = s2 - s1` and `dn = sn - sn-1`. A one-sample signal is padded with its value."
  [^doubles signal ^long N side]
  (if (m/one? (alength signal))
    (edge signal N side)
    (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]

      ;; right end
      (when (m/< end N) (let [lst (Array/aget signal (m/dec len))
                              diff (m/- lst (Array/aget signal (m/- len 2)))]
                          (dotimes [id (m/long-sub N end)]
                            (Array/aset out (m/+ end id) (m/+ lst (m/* (m/inc id) diff))))))

      ;; left end
      (when-not (m/zero? start) (let [frst (Array/aget signal 0)
                                      diff (m/- (Array/aget signal 1) frst)]
                                  (dotimes [id start]
                                    (Array/aset out (m/- start id 1) (m/- frst (m/* (m/inc id) diff))))))
      
      out)))

(defn antireflect
  "Returns signal reflected anti-symmetrically (point reflection) around the edge samples.

  ... 2s1-s3 2s1-s2 | s1 s2 ... sn-1 sn | 2sn-sn-1 2sn-sn-2 ...

  Longer paddings continue the pattern around the newly created edge sample. A one-sample signal is padded with its value."
  [^doubles signal ^long N side]
  (if (m/one? (alength signal))
    (edge signal N side)
    (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]

      ;; right end
      (when (m/< end N) (let [end- (m/dec end)]
                          (loop [pos end-
                                 id (long 1)
                                 curr (m/* 2.0 (Array/aget out end-))]
                            (let [npos (m/inc pos)]
                              (when (m/< npos N)
                                (if (m/== id len)
                                  (recur pos 1 (m/* 2.0 (Array/aget out pos)))
                                  (do
                                    (Array/aset out npos (m/- curr (Array/aget out (m/- npos id id))))
                                    (recur npos (m/inc id) curr))))))))

      ;; left end
      (when-not (m/zero? start) (loop [pos start
                                       id (long 1)
                                       curr (m/* 2.0 (Array/aget out start))]
                                  (let [npos (m/dec pos)]
                                    (when (m/not-neg? npos)
                                      (if (m/== id len)
                                        (recur pos 1 (m/* 2.0 (Array/aget out pos)))
                                        (do
                                          (Array/aset out npos (m/- curr (Array/aget out (m/+ npos id id))))
                                          (recur npos (m/inc id) curr)))))))
      
      out)))
