### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.matrix` namespace. 

### Please write a short paragraph for `clay/vector_matrix.clj` book. Do not copy content of the `clay/vector_matrix.clj` file, write only asked chapter. Put text as a Clojure comments. Put actual Clojure code without comment block and markdown wrapper.  Follow the convention of other chapters.

- Write a paragraph for the "#### Matrix Norms and Condition Number" chapter in the Matrix section.
- Describe all mentioned functions, group them in points and describe them in brief
- When possible add latex formulas
- Illustrate with examples in the `utls/example-note` call.

Write only asked chapter. Refer to previous chapters as examples of chapter structure.



```clojure
;; #### Matrix Norms and Condition Number

;; Functions to compute various matrix norms and the condition number. Matrix norms measure the "size" or "magnitude" of a matrix, and the condition number quantifies how sensitive the solution of a linear system is to changes in the input data. These concepts are vital for analyzing the properties of matrices, assessing the numerical stability of algorithms (like solving linear systems or computing inverses), and understanding the behaviour of transformations represented by matrices.

;; ::: {.callout-tip title="Defined functions"}
;; * `norm`
;; * `condition`
;; :::

;; *   **`norm`**: Computes a specific matrix norm of `A`. The type of norm is specified by the `norm-type` argument. Supported types include:
;;     *   `1`: The L1 norm (maximum absolute column sum). $\Vert A \Vert_1 = \max_{j} \sum_{i=1}^m |a_{ij}|$.
;;     *   `:inf`: The L-infinity norm (maximum absolute row sum). $\Vert A \Vert_\infty = \max_{i} \sum_{j=1}^n |a_{ij}|$.
;;     *   `2`: The spectral norm (largest singular value). $\Vert A \Vert_2 = \sigma_{\max}(A)$.
;;     *   `:max`: The maximum absolute value norm. $\Vert A \Vert_{\max} = \max_{i,j} |a_{ij}|$.
;;     *   `:frobenius`: The Frobenius norm. $\Vert A \Vert_F = \sqrt{\sum_{i=1}^m \sum_{j=1}^n |a_{ij}|^2} = \sqrt{\operatorname{tr}(A^T A)}$. This is a special case of the generalized L$_{p,q}$ norm with $p=2, q=2$.
;;     *   `[p, q]`: The generalized L$_{p,q}$ norm. A specific implementation is provided for $[2,2]$ (Frobenius) and $[p,p]$ (entrywise p-norm).
;;     *   `[p]`: The Schatten p-norm, which is the L-p norm of the vector of singular values. $\Vert A \Vert_p = (\sum_{i=1}^{\min(m,n)} \sigma_i^p)^{1/p}$. Includes `:nuclear` or `:trace` norm for $p=1$.
;;
;;     If no `norm-type` is provided, it defaults to the L1 norm (`1`).
;;
;; *   **`condition`**: Computes the condition number of a matrix `A` with respect to a given norm `norm-type`. It is defined as $\operatorname{cond}(A) = \Vert A \Vert \Vert A^{-1} \Vert$. A large condition number indicates that the matrix is close to being singular, and linear systems involving the matrix may be ill-conditioned, meaning small changes in the input can lead to large changes in the output solution.
;;
;;     Defaults to the L2 norm (`2`) for calculation.

(utls/examples-note
  (mat/norm M2x2) ;; Default (L1)
  (mat/norm M3x3 :inf)
  (mat/norm M4x4 2) ;; Spectral norm
  (mat/norm RealMat :frobenius)
  (mat/norm M2x2 :max)
  (mat/norm M3x3 [2 2]) ;; Frobenius
  (mat/norm RealMat [3 3]) ;; Entrywise L3
  (mat/norm M4x4 [1]) ;; Nuclear norm
  (mat/condition M2x2) ;; Default (L2)
  (mat/condition M3x3 :inf)
  (mat/condition M4x4 1))
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-flash-preview-04-17 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (440 444) (445 478) (781 3493))) -->
<!-- End: -->
