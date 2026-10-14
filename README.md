# linalg

Linear algebra for [sysl](https://github.com/sysl-lang/sysl): vectors, matrices, and elimination over
any element type that behaves like a number.

```
sh/sysl/linalg/
    matrix.sysl     matrices, and the two products whose results are different types
    vector.sysl     vectors, and the half of a vector space whose product yields a scalar
    scalar.sysl     what a cell has to be — and the one promise no signature can carry
    gauss.sysl      elimination with partial pivoting: solve, det, rank, inverse
    bareiss.sysl    fraction-free elimination, for rings whose division cannot be trusted
    factor.sysl     Householder QR, least squares through it, and Cholesky — at real and f32
    eigen.sysl      symmetric eigenvalues and eigenvectors, by cyclic Jacobi — at real and f32
    svd.sysl        the singular value decomposition, by one-sided Jacobi, and the rank,
                    pseudo-inverse and condition number read off it — at real and f32
    fail.sysl       what a decomposition refuses, and why
    tests.sysl      what all of it claims, run by `sysl test .`
    tests_factor.sysl   the factorizations, against numpy and against their definitions
    tests_eigen.sysl    the eigenproblem, against numpy and against A v = λ v
    tests_svd.sysl      the SVD, against numpy and against A = U diag(σ) Vt and Moore-Penrose
    tests_readme.sysl   every example in this file, printing what this file says it prints
package.hocon       who this package is, and what it needs of the machine
```

The module is **`sh.sysl.linalg`**. It depends on nothing outside the standard library —
`sysl.buf`, `sysl.math` and `sysl.math.complex` are the whole of it.

## Using it

```
dependencies {
  linalg { git = "github.com/sysl-lang/linalg", version = "0.1.1" }
}
```

```
import sh.sysl.linalg.{Matrix, Vector, solve}

var a = Matrix.of(2, 2, [2.0, 1.0, 1.0, 3.0])
var b = Vector.of([5.0, 10.0])

solve(a, b) match
    Ok(x) -> print(f"x = $x")
    Err(e) -> print(f"no solution: $e")
```

```output
x = (1, 3)
```

The factorizations take any width `sysl.math.Float` covers — `real` and `f32` — since a reflector's
length and a Cholesky pivot are square roots, and a rank decision is judged against the width's own
`epsilon`:

```
import sh.sysl.linalg.{Matrix, Vector, qr, lstsq, cholesky}

// Fit a line y = c0 + c1 t through four points, which no line passes through exactly.
var a = Matrix.of(4, 2, [1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0])
var y = Vector.of([1.0, 2.9, 5.1, 7.0])

lstsq(a, y) match
    Ok(c) -> print(f"c = ${c}%.6s")
    Err(e) -> print(f"no unique fit: $e")      // RankDeficient, ShapeMismatch

qr(a) match
    Ok((q, r)) -> print(f"Q = ${q}%.6s\nR = ${r}%.6s")
    Err(e) -> print(f"$e")

cholesky(Matrix.of(2, 2, [4.0, 2.0, 2.0, 3.0])) match
    Ok(l) -> print(f"L = ${l}%.6s")
    Err(e) -> print(f"$e")                     // NotSymmetric, NotPositiveDefinite, NotSquare
```

```output
c = (0.97, 2.02)
Q = (-0.5, 0.67082); (-0.5, 0.223607); (-0.5, -0.223607); (-0.5, -0.67082)
R = (-2, -3); (0, -2.23607)
L = (2, 0); (1, 1.41421)
```

A symmetric matrix's eigenvalues come back ascending, with the eigenvectors as the columns of an
orthogonal matrix in the same order; `eigvalsh` answers the values alone:

```
import sh.sysl.linalg.{Matrix, eigh, eigvalsh}

var a = Matrix.of(3, 3, [2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0])

eigh(a) match
    Ok((values, vectors)) -> print(f"λ = ${values}%.6s\nV = ${vectors}%.6s")
    Err(e) -> print(f"$e")                     // NotSymmetric, NotSquare, NotConverged

eigvalsh(a) match
    Ok(values) -> print(f"λ = ${values}%.6s")
    Err(e) -> print(f"$e")
```

```output
λ = (1, 1, 4)
V = (0.707107, -0.408248, 0.57735); (-0.707107, -0.408248, 0.57735); (0, 0.816497, 0.57735)
λ = (1, 1, 4)
```

Any matrix, tall or wide, has a singular value decomposition `A = U diag(σ) Vt` — the thin one, with
`min(m, n)` values in descending order. `singular_values` answers them alone, and the rank, the
pseudo-inverse and the condition number are read off them:

```
import sh.sysl.linalg.{Matrix, svd, singular_values, svd_rank, pinv, cond}

var a = Matrix.of(2, 3, [3.0, 2.0, 2.0, 2.0, 3.0, -2.0])

svd(a) match
    Ok((u, s, vt)) -> print(f"σ = ${s}%.6s\nU = ${u}%.6s\nVt = ${vt}%.6s")
    Err(e) -> print(f"$e")                     // NotConverged

singular_values(a) match
    Ok(s) -> print(f"σ = ${s}%.6s")
    Err(e) -> print(f"$e")

// The second row is twice the first: rank 1, and no finite condition number.
var b = Matrix.of(2, 2, [1.0, 2.0, 2.0, 4.0])

svd_rank(b) match
    Ok(r) -> print(f"rank = $r")
    Err(e) -> print(f"$e")

cond(b) match
    Ok(k) -> print(f"cond = $k")
    Err(e) -> print(f"$e")

pinv(b) match
    Ok(p) -> print(f"pinv = $p")
    Err(e) -> print(f"$e")
```

```output
σ = (5, 3)
U = (0.707107, 0.707107); (0.707107, -0.707107)
Vt = (0.707107, 0.707107, 0); (0.235702, -0.235702, 0.942809)
σ = (5, 3)
rank = 1
cond = inf
pinv = (0.04, 0.08); (0.08, 0.16)
```

A sign of a column of `U`, paired with the same row of `Vt`, is not normalized — another
implementation may print both negated. Every block above, with the output under it, is a test in
`tests_readme.sysl`.

## One body per operation, not one per element type

`solve`, `det`, `rank` and `inverse` are written once against a bound and run at `real`, at `f32` and
over `Complex[F]` unchanged. That is the claim the test suite is built to make: every section after
the first runs the *same* routines at a different element type, so one that had quietly specialised to
the reals would pass the first section and fail the complex one.

**Two bounds, and the difference between them is a safeguard rather than a taxonomy.**

`Scalar` is what a cell has to be to live in a matrix at all: it adds, subtracts, multiplies, divides,
negates, has both identities, compares, prints, and has a **magnitude**. That last one is not an
operator, and it is what pivoting needs — Gaussian elimination chooses the largest remaining cell in a
column, and the complexes have no `<` at all, because no order on the plane respects arithmetic. An
associated `Size` is what lets a `Complex[F]` measure in `F` while an integer measures in itself.

`Field` adds the one promise no signature can carry: that `/` answers the quotient it was asked for.
`int` satisfies `Scalar` and is deliberately **not** a `Field`, so `solve` on a `Matrix[int]` is a
**compile error** naming the promise `int` cannot make — rather than a routine that truncates every
multiplier it forms and returns a plausible wrong answer.

For that case there is `bareiss.sysl`, whose every intermediate entry is a minor of the original
matrix, so every division it performs comes out exact. It asks for no `Field`, no `Magnitude`, no
tolerance and no ordering.

## Where it came from

This was `sysl.math.matrix` in the standard library. It moved out because it was a **leaf**: nothing
else in `library/sysl` imported it, and it was the only *domain* in a library otherwise made of the
platform (`io`, `fs`, `process`, `time`) and of things all code touches (`slices`, `text`, `buf`).
`sysl.math.complex` stays where it is — a complex *type* is what standard libraries across languages
consistently carry, and the algorithms over it are what they consistently leave to packages.

The extraction found one thing worth recording: `vector.sysl` and `tests.sysl` called `.sqrt()` while
importing nothing that brought `sysl.math.Float` into scope, and compiled anyway — they were reaching
it through the standard library's own scope, as a module living under `sysl.math`. Both files declare
that import now, which is what they always needed.

## Licence

ISC. See `LICENSE`.
