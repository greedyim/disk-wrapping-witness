# disk-wrapping-witness

A certified quantitative lower bound for the maximum volume wrappable by a disk.

This repository accompanies the preprint **“A Certified Lower Bound for the Maximum Volume Wrappable by a Disk.”** It contains a 17-vertex, 30-face nonconvex polyhedral witness together with two independent exact-arithmetic verification programs.

## Main result

For the closed disk \(D_R\) of radius \(R\), write

\[
W(R)=\sup\{\operatorname{Vol}(B): B \text{ is wrappable by } D_R\}
\]

in the 1-Lipschitz wrapping formalization used by Nandakumar. The certified construction proves

\[
W(R) \ge C R^3,
\qquad
C=0.17984349416291615\ldots,
\]

where \(C\) is represented exactly by a rational number in the certificate.

More concretely, the source disk maps continuously and 1-Lipschitz onto an embedded nonconvex polyhedral 2-sphere with 17 vertices and 30 triangular faces.

## Why this is stronger than Nandakumar's Question 1

Nandakumar's Question 1 asks whether \(D_\pi\) can wrap *some* nonconvex solid with volume greater than \(4\pi/3\). Under the stated 1-Lipschitz formalization, that qualitative question has a short soft affirmative answer and does not require the polyhedral certificate.

Take an oblate spheroid with equatorial radius \(23/20\) and polar radius \(4/5\). Its volume is greater than that of the unit ball, while a Cauchy-Schwarz estimate shows that its north-south meridian length is strictly less than \(\pi\). A sufficiently small inward rotational notch preserves both strict inequalities and makes the body nonconvex. Any surface of revolution whose generating north-south meridian has length \(L\le R\) is wrapped by \(D_R\) via a direct radial 1-Lipschitz map.

Accordingly, the point of this repository is the **substantially stronger quantitative bound** \(C=0.179843494\ldots\), not merely existence of a nonconvex improvement over the sphere.

For comparison, the classical Mylar balloon is the sharp convex body-of-revolution benchmark

\[
V_{\mathrm{Mylar}}=0.152315527\ldots R^3,
\]

so the certified coefficient is about 18.07% larger. No claim is made that Mylar is optimal among all convex bodies wrappable by a disk.

## Preprint

- [`paper/preprint.tex`](paper/preprint.tex)

Compile from the repository root so the figure paths resolve:

```bash
pdflatex paper/preprint.tex
pdflatex paper/preprint.tex
```

The compiled PDF is distributed alongside the current source package and can be placed at `paper/preprint.pdf`.

## Exact verification

The repository contains two independent exact verifiers. They check:

1. the quotient is a connected orientable triangulated 2-sphere with \((V,E,F)=(17,45,30)\);
2. all 30 affine triangle maps are strictly 1-Lipschitz;
3. all 435 unordered target-triangle pairs intersect exactly in their prescribed common simplex;
4. an exact face-plane witness proves nonconvexity; and
5. the enclosed volume equals the exact rational constant used in the paper.

The important checks use explicit `require(...)/raise` logic rather than Python `assert`, so they remain active under `python -O`.

```bash
# Verifier A: Python 3 + SymPy
python certificate/exact_certificate.py

# Verifier B: Python 3 standard library only
python certificate/independent_verify.py
```

Both report

```text
volume decimal = 0.17984349416291615
```

and the exact rational value

```text
359686988325832312742347877865122384142690317017149
/ 2000000000000000000000000000000000000000000000000000
```

## Status

Research preprint / external-review candidate. The current claim is deliberately quantitative: an explicit, reproducible lower bound for disk-wrapping volume. It does not determine the optimum and does not solve Nandakumar's general conjecture.

Feedback, independent verification, and pointers to prior equivalent constructions or bounds are very welcome.
