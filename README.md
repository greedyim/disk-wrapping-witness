# disk-wrapping-witness


A small explicit, exactly certified nonconvex wrapping of a disk.

This repository accompanies the preprint **“A Certified Nonconvex Wrapping of a Disk with Volume Greater than the Sphere.”** It contains a 17-vertex, 30-face polyhedral witness together with two independent exact-arithmetic verification programs.

## Main result

For every source radius \(R>0\), the construction gives a continuous 1-Lipschitz map from the closed disk \(D_R\subset\mathbb R^2\) onto an embedded nonconvex polyhedral 2-sphere bounding a solid of volume

\[
V = C R^3,
\qquad
C = 0.17984349416291615\ldots
\]

where \(C\) is represented exactly by a rational number in the certificate.

For \(R=\pi\), this exceeds \(4\pi/3\), the volume of the unit ball. The construction is intended as a **small explicit and exactly certified witness** for the disk-wrapping question discussed by Nandakumar. We do **not** claim priority for abstract existence alone; related volume-increasing isometric-deformation results of Bleecker and Pak may imply the general existence phenomenon indirectly.

For comparison, the classical Mylar balloon gives the sharp convex axisymmetric benchmark

\[
V_{\mathrm{Mylar}} = 0.152315527\ldots R^3,
\]

so the certified nonconvex witness has about 18.07% larger volume.

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

Both report the decimal volume

```text
0.17984349416291615
```

and the exact rational value

```text
359686988325832312742347877865122384142690317017149
/ 2000000000000000000000000000000000000000000000000000
```

## Status

This is a research preprint / external-review candidate. The claim is intentionally narrower than a general solution of the disk-wrapping maximization problem: the contribution is the explicit witness, its quantitative lower bound, and a reproducible exact certificate.

Feedback, independent verification, and pointers to prior equivalent constructions are very welcome.
