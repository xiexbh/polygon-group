# Supplementary computational materials

Baohua Xie — accompanying the revised polygon-group manuscript (2026).

This folder is self-contained and ready to upload to a repository. No repository URL or DOI is asserted before publication. The paper cites these files as *Supplementary computational materials*.

## Run

Use Mathematica 14.0 or a compatible Wolfram kernel. No optional packages or network access are required.

```text
math -script section6_verification.wl
math -script side_verification.wl
math -script boundary_cells.wl
```

Alternatively evaluate `Get["/absolute/path/to/the/file.wl"]` in Mathematica. Each script is standalone and clears `Global` symbols; use a separate kernel if needed. A failed identity, unresolved elimination, or time limit terminates with a failure. A timeout is never accepted as a proof. Runtime depends on the computer; the per-check time limit is explicitly set in each script.

## What is verified

| File | Successful checks | Role in the manuscript |
| --- | ---: | --- |
| `section6_verification.wl` | 126 | Normalization and symmetry identities; optional checks of analytic scalar estimates; local row determinants and witnesses; mirror-circle identities; complete triangle inequalities on the closed ball; exclusion of additional relatively open ridges. |
| `side_verification.wl` | 30 | Sufficient inequalities on entire sides (Proposition 6.15); parameter-free local implications and the two remote F-pair exclusions. |
| `boundary_cells.wl` | 16 | Strictness off the triangle cutting lines; lower-dimensional intersections of discarded pairs; absence of additional ideal strata (Proposition 6.16). |

All parameters and variables are real. The parameter interval is `3/5 <= c < 1`; since `c = cos(2*pi/p) > 3/5` for every integer `p >= 7`, the tests cover all required p at once. Polynomial identities use exact arithmetic. Universal inequalities are verified by asking `Reduce[..., Reals]` for counterexamples and requiring `False`. Denominators are checked before clearing them. No floating-point sampling or numerical optimization is used.

The main script uses `N <= 0` for the closed-ball region and no-ridge tests. This is stronger than the open-ball version in the preceding revision. The side script uses single-bisector coordinates `(t,x,y)`; these auxiliary coordinates differ from the two-bisector coordinates `(t,h)` in the main script. Its normalized ball polynomial is affine in c; positivity at c=0 reduces the local implication to the stronger condition at c=1. That reduction is checked explicitly.

## Scope of the certificate

The 172 checks certify the exact algebraic assertions encoded in these files. They do not constitute a formal proof of the entire paper. The real half-space arguments, triangle and hexagon descriptions, side-pairing deductions, Poincare theorem application, and convex-cocompactness argument must be read in the manuscript. In particular, an absence-of-ridge test does not by itself prove that a closed intersection is empty; the boundary script records the stronger statements actually used later.

The three logs record successful runs on Mathematica 14.0. `verification-summary.json` records scope and counts; `SHA256.json` records the delivered files. When editing a script, rerun it and replace its log and hashes before publishing.
