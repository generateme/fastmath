(ns fastmath.protocols.wavelets)

(defprotocol TransformProto
  "Transformer functions."
  (forward-1d [t xs] [t xs options] "Forward transform of sequence or array.")
  (reverse-1d [t xs] [t xs options] "Reverse transform of sequence or array.")
  (forward-2d [t xss] [t xss options] "Forward transform of sequence of sequences.")
  (reverse-2d [t xss] [t xss options] "Reverse transform of sequence of sequences."))

(defprotocol WaveletProto
  "Wavelet filter bank information."
  (wavelet-name [w] "Return name of the wavelet")
  (coeffs-size [w] "Length of the filters")
  (phi [w] [w kind] "Scaling, low-pass, filter (decomposition or reconstruction)")
  (psi [w] [w kind] "Wavelet, high-pass, filter (decomposition or reconstruction)"))
