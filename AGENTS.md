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
| `fastmath.interpolation` | interp | 1d interpolation functions: linear, cubic, monotone, pchip, akima (and makima), sprague, step, barycentric, polynomial, neville, divided difference, b-spline; smoothing: loess, cubic smoothing spline; 2d grid: bilinear, bicubic, cubic-2d; scattered data: microsphere projection, shepard, rbf, kriging, gaussian processes (`fastmath.interpolation.gp`); isotonic: pava, cir; extrapolation |
| `fastmath.kernel` | k, kernel | Kernel functions: rbf (`fastmath.kernel.rbf`), vector/covariance kernels (`fastmath.kernel.vector`), KDE (`fastmath.kernel.density`), window (`fastmath.kernel.window`), kernel combinators (`sum`, `product`, `wmean`); the variograms are not re-exported here |
| `fastmath.kernel.variogram` | vgram | semivariogram models, empirical variogram (classical, Cressie, Dowd, Genton, quantile and M-robust estimators), model fitting, variogram to covariance |
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

Use:

- `fastmath.core` (alias `m`) for every primitive math operations.
- `fastmath.java.Array` (import it first) for `Array/aget` and `Array/aset`.
- `^{:private true :const true :tag 'double}` metadata for constant, private, double values (use `'long` for long values, etc.)

### Naming & Formatting

- All names use `kebab-case`; predicate functions end with `?`.
- Private helpers: `defn ^:private foo`.
- Numeric constants use `{:const true}` metadata.
- Type-hint aggressively: `^double`, `^long`, `^doubles`, `^Vec2`, etc. on params and return types.
- 2-space indentation; no trailing commas.

### Auditing

Audit is a proces of verification of numerical correctness and completeness of given namespace or namespaces, writing docstings, writing tests, filling gaps and implementing low hanging fruits.

### Documentation of development

Use `obisdian` vault to store all information about the development process: findings, spec, progress log, overview and maybe other related files documenting the process.

### Writing tests

- For a reference values in tests you can use (with caution):
  - `Rscript` - for any R package
  - `/home/ts/.juliaup/bin/julia` - for any Julia package
  - `uv` Python package manager run in `/home/ts/penv` environment - for any Python package 
  - `octave` - for any Octave/Matlab package
- Use `m/delta-eq` or `v/delta-eq`  with the default tolerance for anything closed-form; use an explicit looser tolerance (with a one-line comment stating why and the empirically observed error) when needed
- when creating EDN reference test file, always keep script in the `utils/fastmath/dev` folder for reproducibility
- more in `fastmath-testing` skill

## Skills

Available skill for this project

- `docstring` - load this skill for function or var documentation.
- `obsidian` - for general topic memory store and retrieve in Obsidian
- `clojure-eval` - evaluate Clojure against a running nREPL via `clj-nrepl-eval`, to verify edits compile and behave as expected
- `fastmath-testing` - instructions about writing tests and fixing bugs

## Gotchas

- never commit to the fastmath repo, it requires signing
- add to the CHANGELOG.md only the most important changes or additions, concise as much as possible (one liners), only things which differ from previous release
