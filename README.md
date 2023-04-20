[![Clojars Project](https://img.shields.io/clojars/v/generateme/fastmath.svg)](https://clojars.org/generateme/fastmath)

# fastmath

Math library.

## Documentation

[Documentation with examples](https://generateme.github.io/fastmath/index.html)

### Clerk notebooks - towards 2.2.0

[Index](https://generateme.github.io/fastmath/notebooks/)

* [fastmath.core](https://generateme.github.io/fastmath/notebooks/notebooks/core.html)
* [fastmath.random](https://generateme.github.io/fastmath/notebooks/notebooks/random.html)
* [fastmath.stats](https://generateme.github.io/fastmath/notebooks/notebooks/stats.html)
    * [fastmath.stats.bootstrap](https://generateme.github.io/fastmath/notebooks/notebooks/bootstrap.html)


### Previous version

Based on SMILE 1.5.3 and including: regression and classification.

```clojure
[generateme/fastmath "1.5.3"]
```

[Documentation with examples](https://generateme.github.io/fastmath/1.5/index.html)


## Installation

```clojure
[generateme/fastmath "2.1.8"]
```

Snapshot:

```clojure
[generateme/fastmath "2.1.9-SNAPSHOT"]
```

### MKL

**Important Note**

Fastmath relies on [SMILE](https://haifengl.github.io/) 2.6.0 which relies on BLAS/LAPACK via MKL and/or OpenBlas. MKL (preferred) and OpenBlas are included as dependencies in `fastmath`. This leads to addional 1GB of jar files. I can't assure that all functionalities will work on OpenBlas (see: [#15](https://github.com/generateme/fastmath/issues/15#issuecomment-1090323385)) but 99% should. If you need `fastmath` to be lighter, please exclude MKL from your path.

#### lein / project.clj

```clojure
[generateme/fastmath "2.1.8" :exclusions [com.github.haifengl/smile-mkl]]
```

#### deps.edn

```clojure
{:deps {generateme/fastmath {:mvn/version "2.1.8"
                             :exclusions [com.github.haifengl/smile-mkl]}}}
```

If you don't need certain interpolation or clustering methods you can exclude OpenBlas as well (be warned that other things can break):

```clojure
:exclusions [com.github.haifengl/smile-mkl org.bytedeco/openblas]
```

#### MKL Exception

When MKL is not available `fastmath` (SMILE actually, [here](https://github.com/haifengl/smile/blob/master/base/src/main/java/smile/math/blas/BLAS.java#L58) and [here](https://github.com/haifengl/smile/blob/master/base/src/main/java/smile/math/blas/LAPACK.java#L60)) will throw two exceptions with full stack traces about lack of MKL. You can safely ignore them.

```
[main] DEBUG smile.math.blas.LAPACK - Failed to create MKL instance:
java.lang.ClassNotFoundException: smile.math.blas.mkl.MKL
[...]

[main] DEBUG smile.math.blas.BLAS - Failed to create MKL instance:
java.lang.ClassNotFoundException: smile.math.blas.mkl.MKL
[...]
```

## Content

### [PrimitiveMath](https://github.com/ztellman/primitive-math)

Code adopted from Zach Tellmans' library.

```
[* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? bool-and bool-or bool-xor bool-not << >> >>> not==]
```

### Math functions

* Trigonometric functions
* Power: log, ln, logb, exp, pow, sqrt
* Rounding functions: round, floor, ceil, trunc, frac, approx + other
* Normalizations: norm  wrap, constrain
* Interpolations: lerp, cos-interpolation, smooth-interpolation, quad-interpolation, smoothstep
* Special functions: erf, beta, gamma + other
* Distance: dist, hypot
* Sign: sgn, signum, abs
* Other: gcd, lcm

Most of them backed by [Jafama FastMath 2.3.1](https://github.com/jeffhain/jafama) or [Apache Commons Math 3.6.1](http://commons.apache.org/proper/commons-math/index.html)

### Vector operations protocol and implementations

* 2d (`Vec2`), 3d (`Vec3`) and 4d (`Vec4`) vector types.
* ArrayVector for fixed length long vectors (fixed sized double-array)
* Clojure vector, double array, Number (as 1d vector)

### Matrix operations protocol and implementations

* 2d (`Mat2x2`), 3d (`Mat3x3`) and 4d (`Mat4x4`) matrix types.

With typical basic matrix operations

### Complex number functions

* primitive operations: mult, div, add, sub
* abs, arg, conjugate, reciprocal, neg
* atan, asin, acos, csc, sec, tanh, tan, sinh, sin, cosh, cos
* log, exp, pow
* sqrt, sq, sqrt1z

### Random numbers

* Collection of random number generators
* Random generator functions for each primitive type (drand - double, lrand - long, frand - float, irand - int)
* Additional RNG functions: brand - true/false, grand - gaussian distributed double
* Random sequences: from distribution, halton, sobol, R2, sphere, uniform

### Distributions

* Collection of distributions (60+)

### Noise

* 4 noise types: value, gradient, simplex, discrete
* 3 noise blends: fbm, ridgedmulti, billow
* Ready to use fbm functions: noise (perlin), vnoise (value noise), simplex
* Warp noise

### Statistics

* Descriptive statistics: size, min, max, mode, mean, median, percentiles, kurtosis, skewness, IQR, LAV, UAV and other
* Correlations
* t-test
* ACF/PACF
* histogram
* Bootstrap
* Confidence intervals

#### Bootstrap

Bootstrap functions and confidence intervals

### Interpolations

1d, 2d interpolations

### Easings

Several easing functions (in, out, in-out)

### Transforms

* Wavelets: 1d, 2d (haar, biorthogonal, symlet, coiflet, daubechies, legendre)
* 1d Fast Sine, Cosine and Hadamard

### Vector fields

Great collection (100+) of R^2->R^2 functions.

### Gaussian Processes

Gaussian Processes

### Clustering

Various clustering algorithms including K-Means++, DBSCAN, CLARANS, DENCLUE, MEC, Spectral, Deterministic Annealing

### Optimization

L-BFGS-B, Gradient, Nelder-Mead, Simplex, Powell, BOBYQA, CMAES, BayesianOptimizer

### Grids

Hexagonal, squared, triangular, rhomboidal grid functions

### Kernels

Collection of various kernels (density, RBF, correlation)

### Signal

* Signal (audio) processing filters and oscillators
* Smoothing filters: Savitzky-Golay, moving average, kernel smoothing

### EFloat

Floating point operations with error bounds

### Other

Plenty of constant values

Almost all functions optimized to work with `double` and `long` primitives

## Supporting libraries

* [Apache Commons Math 3.6.1](http://commons.apache.org/proper/commons-math/index.html) - Apache 2.0 Licence
* [SMILE 2.5.0](http://haifengl.github.io/smile/) - Apache 2.0 Licence
* [Jafama FastMath 2.3.1](https://github.com/jeffhain/jafama) - Apache 2.0 Licence
* [PrimitiveMath](https://github.com/ztellman/primitive-math) - MIT Licence
* [JWave](https://github.com/cscheiblich/JWave/) - MIT Licence
* [SSJ](https://github.com/umontreal-simul/ssj) - Apache 2.0 Licence

## Alternatives

Since this library is only JVM version, you can check following Clojure/ClojureScript libraries as replacement

* [PrimitiveMath](https://github.com/ztellman/primitive-math) - for primitive operators
* [Kixi stats](https://github.com/MastodonC/kixi.stats) - for pure clj(s) statistics/distributions (tends to be 2-10x slower)
* [thi.ng](http://thi.ng/) - for vectors, general math, noise, complex numbers, transforms (fourier)
* [vectorz-clj](https://github.com/mikera/vectorz-clj) - fast vector operations
* [Incanter](https://github.com/incanter/incanter) - statistics/distributions/probability

## Java

Java classes are compiled for java 1.8

## How To Help

If you see place of improvement, I'm accepting PRs.

## Licence

* The Unlicence - up to 1.5.2
* MIT License - from 1.5.3

Copyright (c) 2023 generateme
