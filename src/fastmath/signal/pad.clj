(ns fastmath.signal.pad
  (:require [fastmath.core :as m])
  (:import [java.util Arrays]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn pad-position
  "Calculate position of original signal in the padded one based on `side` method (`:left`, `:right`, `:both`)."
  ^long [side ^long N ^long len]
  (case side
    :left (m/- N len)
    :right 0
    (m/round-even (m// (m/- N len) 2.0))))

(defn out-len-pos
  "Returns new, zero-padded signal with information about position." 
  [^doubles signal ^long N side]
  (let [len (alength signal)
        out (double-array N)
        start (pad-position side N len)
        end (m/+ start len)]
    (System/arraycopy signal 0 out start len)
    [out len start end]))

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
    (when (< end N) (Arrays/fill out end N (Array/aget signal (m/dec len))))
    (when-not (m/zero? start) (Arrays/fill out 0 start (Array/aget signal 0)))
    out))

(defn periodic
  "Returns signal padded with copies of the original signal.

  ... sn s1 ... sn | s1 ... sn | s1 ... sn s1 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]
    (when (< end N) (let [diff (m/- N end)
                          copies (m/quot diff len)
                          remainder (m/mod diff len)]
                      (dotimes [id copies]
                        (System/arraycopy signal 0 out (m/+ end (m/* id len)) len))
                      (System/arraycopy signal 0 out (m/+ end (m/* copies len)) remainder)))
    (when-not (m/zero? start) (let [copies (m/quot start len)
                                    remainder (m/mod start len)]
                                (dotimes [id copies]                                  
                                  (System/arraycopy signal 0 out (m/- start (m/* (inc id) len)) len))
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
  "Returns signal padded with mirrored copies of the orignal signal.

  ... sn sn ... s1 | s1 ... sn | sn ... s1 s1 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
        ^doubles mirror (reverse-signal signal)]

    ;; right end
    (when (< end N) (let [diff (m/- N end)
                          copies (m/quot diff len)
                          remainder (m/mod diff len)]
                      (dotimes [id copies]
                        (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/+ end (m/* id len)) len))
                      (System/arraycopy (if (m/even? copies) mirror signal) 0 out (m/+ end (m/* copies len)) remainder)))

    ;; left end
    (when-not (m/zero? start) (let [copies (m/quot start len)
                                    remainder (m/mod start len)]
                                (dotimes [id copies]
                                  (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (inc id) len)) len))
                                (System/arraycopy (if (m/even? copies) mirror signal) (m/- len remainder) out 0 remainder)))
    out))

(defn reflect
  "Returns signal padded with reflected copies of the orignal signal.

  ... sn-1 sn ... s1 s2 | s1 ... sn | sn-1 ... s2 s1 s2 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
        ^doubles mirror (reverse-signal signal)
        len- (m/dec len)]

    ;; right end
    (when (< end N) (let [diff (m/- N end)
                          copies (m/quot diff len-)
                          remainder (m/mod diff len-)]
                      (dotimes [id copies]
                        (System/arraycopy (if (m/even? id) mirror signal) 1 out (m/+ end (m/* id len-)) len-))
                      (System/arraycopy (if (m/even? copies) mirror signal) 1 out (m/+ end (m/* copies len-)) remainder)))

    ;; left end
    (when-not (m/zero? start) (let [copies (m/quot start len-)
                                    remainder (m/mod start len-)]
                                (dotimes [id copies]
                                  (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (inc id) len-)) len-))
                                (System/arraycopy (if (m/even? copies) mirror signal) (m/- len- remainder) out 0 remainder)))
    out))

(defn- negate-signal
  [^doubles signal]
  (let [len (alength signal)
        nsignal (double-array len)]
    (dotimes [id len]
      (Array/aset nsignal id (m/- (Array/aget signal id))))
    nsignal))

(defn antisymmetric
  "Returns signal padded with mirrored and negated copies of the orignal signal.

  ... sn -sn ... -s1 | s1 ... sn | -sn ... -s1 s1 ...."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)
        ^doubles mirror (negate-signal (reverse-signal signal))]

    ;; right end
    (when (< end N) (let [diff (m/- N end)
                          copies (m/quot diff len)
                          remainder (m/mod diff len)]
                      (dotimes [id copies]
                        (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/+ end (m/* id len)) len))
                      (System/arraycopy (if (m/even? copies) mirror signal) 0 out (m/+ end (m/* copies len)) remainder)))

    ;; left end
    (when-not (m/zero? start) (let [copies (m/quot start len)
                                    remainder (m/mod start len)]
                                (dotimes [id copies]
                                  (System/arraycopy (if (m/even? id) mirror signal) 0 out (m/- start (m/* (inc id) len)) len))
                                (System/arraycopy (if (m/even? copies) mirror signal) (m/- len remainder) out 0 remainder)))
    
    out))

(defn linear
  "Returns signal padded with linear extension of the edges."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]

    ;; right end
    (when (< end N) (let [lst (Array/aget signal (m/dec len))
                          diff (m/- lst (Array/aget signal (m/- len 2)))]
                      (dotimes [id (m/- N end)]
                        (Array/aset out (m/+ end id) (m/+ lst (m/* (m/inc id) diff))))))

    ;; left end
    (when-not (m/zero? start) (let [frst (Array/aget signal 0)
                                    diff (m/- (Array/aget signal 1) frst)]
                                (dotimes [id start]
                                  (Array/aset out (m/- start id 1) (m/- frst (m/* (m/inc id) diff))))))
    
    out))

(defn antireflect
  "Returns signal reflected anti-symmetrically around the edges.

  ... 2*s1-sn ... 2*s1-s2 | s1 ... sn | 2*sn-sn-1 ... 2*sn-s1 ..."
  [^doubles signal ^long N side]
  (let [[^doubles out ^long len ^long start ^long end] (out-len-pos signal N side)]

    ;; right end
    (when (< end N) (let [end- (m/dec end)]
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
    
    out))
