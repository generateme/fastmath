# fastmath

Fast and primitive based math library.

Originally it was a part of generative art/glich [Clojure2d](https://github.com/Clojure2D/clojure2d) library.

## Installation

```clojure
[generateme/fastmath "1.1.2"]
```

## Documentation

[HERE](https://generateme.github.io/fastmath/index.html)

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
* Other: gcd

Most of them backed by [Jafama FastMath 2.3.1](https://github.com/jeffhain/jafama) or [Apache Commons Math 3.6.1](http://commons.apache.org/proper/commons-math/index.html)

### Vector operations protocol and implementations

* 2d (`Vec2`), 3d (`Vec3`) and 4d (`Vec4`) vector types.
* ArrayVector for fixed length long vectors (fixed sized double-array)
* Clojure vector

With following groups of functions:

* Basic linear operations: add, mult, div, sub, dot, cross, hadamard product
* mag, magsq, heading, angle-between, limit, normalize
* rotations (2d,3d)
* translations (2d, 3d)
* various distances
* centroid, interpolations
* and other

### Complex number functions

* primitive operations: mult, div, add, sub
* abs, arg, conjugate, reciprocal, neg
* atan, asin, acos, csc, sec, tanh, tan, sinh, sin, cosh, cos
* log, exp, pow
* sqrt, sq, sqrt1z

### Random numbers

* Collection of random number generators
* Collection of distributions (count: 28)
* Random generator functions for each primitive type (drand - double, lrand - long, frand - float, irand - int)
* Additional RNG functions: brand - true/false, grand - gaussian distributed double
* Random sequences: from distribution, halton, sobol, sphere, uniform

### Noise

* 4 noise types: value, gradient, simplex, discrete
* 3 noise blends: fbm, ridgedmulti, billow
* Ready to use fbm functions: noise (perlin), vnoise (value noise), simplex

### Statistics

* Descriptive statistics: size, min, max, mode, mean, median, percentiles, kurtosis, skewness, IQR, LAV, UAV and other
* Correlations

### Interpolations

1d, 2d interpolations

### Easings

Several easing functions (in, out, in-out)

### Transforms

* Wavelets: 1d, 2d (haar, biorthogonal, symlet, coiflet, daubechies, legendre)
* 1d Fast Sine, Cosine and Hadamard

### Vector fields

Great collection (100+) of R^2->R^2 functions.

### Clustering

SMILE bindings for clustering

### Other

Plenty of constant values

Almost all functions optimized to work on `double` and `long` primitives

## Supporting libraries

* [Apache Commons Math 3.6.1](http://commons.apache.org/proper/commons-math/index.html) - Apache 2.0 Licence
* [SMILE 1.5.1](http://haifengl.github.io/smile/) - Apache 2.0 Licence
* [Jafama FastMath 2.3.1](https://github.com/jeffhain/jafama) - Apache 2.0 Licence
* [PrimitiveMath](https://github.com/ztellman/primitive-math) - MIT Licence
* [JWave](https://github.com/cscheiblich/JWave/) - MIT Licence

## Alternatives

Since this library is only JVM version, you can check following Clojure/ClojureScript libraries as replacement

* [PrimitiveMath](https://github.com/ztellman/primitive-math) - for primitive operators
* [Kixi stats](https://github.com/MastodonC/kixi.stats) - for pure clj(s) statistics/distributions (tends to be 2-10x slower)
* [thi.ng](http://thi.ng/) - for vectors, general math, noise, complex numbers, transforms (fourier)
* [vectorz-clj](https://github.com/mikera/vectorz-clj) - fast vector operations
* [Incanter](https://github.com/incanter/incanter) - statistics/distributions/probability

## Java

Java classes are compiled for java 1.8

## TODO

* Change images to plots with axes
* More tests
* Move vector fields

## How To Help

If you see place of improvement, I'm accepting PRs.

## Licence

The Unlicence

