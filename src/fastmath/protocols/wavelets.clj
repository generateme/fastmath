(ns fastmath.protocols.wavelets)

(defprotocol TransformProto
  "Transformer functions."
  (forward-1d [t xs] [t xs options] "Forward transform of sequence or array.")
  (reverse-1d [t xs] [t xs options] "Reverse transform of sequence or array.")
  (forward-2d [t xss] [t xss options] "Forward transform of sequence of sequences.")
  (reverse-2d [t xss] [t xss options] "Reverse transform of sequence of sequences."))

(defprotocol WaveletProto
  "Wavelet information and operations"
  (wavelet-name [w] "Return name of the wavelet")
  (wavelet-forward [w signal length] "Forward transform, length - how much data to transform")
  (wavelet-reverse [w signal length] "Reverse transform, length - how much data to transform")
  (coeffs-size [w] "Length of the coefficients")
  (phi [w] [w kind] "Scaling, low-pass, coefficients (deconstruction or reconstruction)")
  (psi [w] [w kind] "Wavelet, high-pass, coefficients (deconstruction or reconstruction)"))
