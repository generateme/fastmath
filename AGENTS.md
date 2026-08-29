# AGENTS.md — fastmath

## General info

Fastmath is a general math library.
Mathematics, statistics, machine learning, numerical methods/computing.

### Architecture

- **Foundation:** `fastmath.core` (alias `m`) wraps `fastmath.java.PrimitiveMath` (Java); busiest dependency in the codebase, used by nearly every other namespace.
- **Abstraction layer:** `fastmath.protocols` (`RNGProto`, `VectorProto`, `DistributionProto`, `GridProto`, ...), `fastmath.protocols.matrix` (`MatrixProto`), `fastmath.protocols.polynomials`, `fastmath.protocols.wavelets` — decouple implementations from consumers.
- **Multi-file namespace packages** (dir of files aggregated by a top-level `.clj`): `fields/` (`a.clj`…`z.clj`, one file per letter), `interpolation/`, `calculus/`, `signal/`, `kernel/`, `special/`, `optimization/` (`lbfgsb.clj`, `bo.clj`), `ml/` (`regression.clj` + `regression/terms.clj`, `regression/contrast.clj`, `clustering.clj`), `dual/` (`partials.clj`), `transform/` (`wavelets.clj`).
- **Java interop:** `src/fastmath/java/` — `PrimitiveMath.java`, `Array.java`, `Monotone.java`, `R2.java`, `noise/` (Perlin/value/gradient noise backing `fastmath.random`).
- **Vendored dependencies (bundled in-repo):** `LBFGSBJava/` (L-BFGS-B optimizer, used only by `fastmath.optimization.lbfgsb`); `fastmath.core.matrix/` (`clojure.core.matrix` protocol-extension shim).
- **Docs/notebooks** (`clay/`, `notebooks/`, `metadoc/`, `docs/`) are downstream consumers of the public API (Clay/Codox-generated), not part of the library.

### Namespaces

| Namespace  | Alias | Content |
|:-----------|:--|:--|
| `fastmath.core` | m | basic math (arithmetic, trigonometric/hyperbolic, log/power), bitwise, floating point ops, combinatorics, lerp, gcm/lcm |
| `fastmath.vector` | v | 2,3,4d and general vectors, vector ops, rotations, distances, elementwise ops |
| `fastmath.matrix` | mat | 2,3,4d and general matrices, matrix ops, rotations, metrics, decompositions, elemenentwise ops |
| `fastmath.random` | r | RNGs, random and low-discrepancy sequences, noise, probabilistic distributions |
| `fastmath.stats` | stats | descriptive statistics, quantiles, extents, data transformations, correlation, similarity, contingency tables, classification metrics (`fastmath.stats.binary`), effect size, statistical tests, acf/pacf, histograms |
| `fastmath.stats.bootstrap` | boot | bootstrap data with confidence intervals |
| `fastmath.complex` | cplx | complex number ops |
| `fastmath.quaternion` | quat | quaternion ops |
| `fastmath.special` | spec | special functions: bessel, elliptic, jacobi, gamma, beta, erf, airy, zeta, integrals (Si, Ci, li/Li, Ei, Ein, En), hypergeometric, lambert, minkowski, harmonic, owen's t |
| `fastmath.calculus` | calc | integration: vegas/vegas+, h-cubature, quadrature, romberg/trapezoid/simson/midpoint; finite differences derivatives |
| `fastmath.distance` | dist | distance, metric functions |
| `fastmath.easings` | ease | easing functions |
| `fastmath.efloat` | efloat | floating point ops with error tracking |
| `fastmath.grid` | grid | euclidean (2d) grid cells: square, hexagonal, rhombus, triangle and coordinate ops |
| `fastmath.interpolation` | interp | 1d and nd interpolation functions: linear, cubic, monotone, sprague, step, barycentric, neville, divided difference, loess, akima, b-spline, microsphere projection, bicubic, bilinear, shepard, rbf, kriging, gaussian processes; isotonic: pava, cir; extrapolation |
| `fastmath.kernel` | k, kernel | Kernel functions: rbf (`fastmath.kernel.rbf`), covariance (`fastmath.kernel.vector`), KDE (`fastmath.kernel.density`), window (`fastmath.kernel.vector`), variograms (`fastmath.kernel.variogram`) with estimation and fitting |
| `fastmath.ml.regression` | reg, regr | LM, GLM |
| `fastmath.ml.regression.contrast` | contrast | categorical encoding functions |
| `fastmath.ml.clustering` | clust | kmeans++, fuzzy-kmeans, dbscan |
| `fastmath.optimization` | opt | brent, bobyqa, powell, nelder-mead, simplex, cmaes, gradient, l-bfgs-b, bayesian optimisation, linear problem solver |
| `fastmath.polynomials` | poly | real and complex polynomial evaluation and ops, orthogonal polynomials: bernstein, laguerre, chebyshev, legendre, gegenbauer, hermite, jacobi, bessel, meixner-pollaczek, ince beams |
| `fastmath.solver` | solver | Root finder: brent, bisection, illinois, muller, muller2, pegasus, regula-falsi, ridders, secant, quatdratic, cubic |
| `fastmath.signal` | signal processing, audio filters, convolution, correlation, iir/fir filters, smoothing, spectrum, stft, periodogram, waveforms, chirp, padding |
| `fastmath.transform` | wavelets: dwt, wpt (`fastmath.transform.wavelets`); dft/fft, cosine, sine, hadamard, hartley |

## Build & Test Commands

- **Run all tests:** `lein test`
- **Run a single namespace:** `lein test fastmath.core-test`
- **Run a single var:** `lein test :only fastmath.core-test/my-test`
- **Lint (Eastwood):** `lein with-profile eastwood eastwood`
- **Build JAR:** `lein jar`

## Code Style

### Namespace & Requires

- Use `fastmath.*` namespace hierarchy; one `:require` / `:import` entry per line.

### Performance (mandatory in computational namespaces)

Every source file starts with:
```clojure
(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
```

### Naming & Formatting

- All names use `kebab-case`; predicate functions end with `?`.
- Private helpers: `defn ^:private foo`.
- Numeric constants use `{:const true}` metadata.
- Type-hint aggressively: `^double`, `^long`, `^doubles`, `^Vec2`, etc. on params and return types.
- 2-space indentation; no trailing commas.

### Skills

Available skill for this project

- `docstring` - Load this skill for function or var documentation.
- `obsidian` - for general topic memory store and retrieve in Obsidian
- `codebase-memory` - knowledge-graph tools (search_graph, trace_path, get_architecture, ...) for structural code exploration; project indexed as `home-ts-clojure-fastmath`
- `clojure-eval` - evaluate Clojure against a running nREPL via `clj-nrepl-eval`, to verify edits compile and behave as expected
