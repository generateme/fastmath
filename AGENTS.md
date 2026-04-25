# AGENTS.md — fastmath

## General info

Fastmath is general math library.
Use common math, statistics, machine learning, numerical computing and related knowledge.

## Build & Test Commands
- **Run all tests:** `lein test`
- **Run a single namespace:** `lein test fastmath.core-test`
- **Run a single var:** `lein test :only fastmath.core-test/my-test`
- **Lint (Eastwood):** `lein with-profile eastwood eastwood`
- **Build JAR:** `lein jar`

## Code Style

### Namespace & Requires
- Use `fastmath.*` namespace hierarchy; one `:require` / `:import` entry per line.
- Standard aliases: `[fastmath.core :as m]`, `[fastmath.vector :as v]`, `[fastmath.random :as r]`.

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
