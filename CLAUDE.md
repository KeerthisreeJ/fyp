# Project context for Claude

23CSE498 Final Year Project, Phase II. Amrita Vishwa Vidyapeetham, Coimbatore.
Guide: Dr. T. Ramraj. Owner: P. M. Radha Krishna (+ 4 teammates).

**The question:** can a descriptor computed from a MOF's connectivity alone — in
seconds, no molecular simulation — carry enough structural information to decide
which frameworks are worth the hours of simulation? The pipeline is
CIF → bond perception → metal-oxo decomposition → labelled quotient graph →
supercell → multifractal spectrum → Δα.

---

## THE OPEN QUESTION — RESOLVED, NEGATIVELY (2026-09-14)

**Δα does not converge with supercell size, at any q range.** Measured on
HKUST-1 with `2_python/code_12_converged_band.py`:

```
diam= 8   Δα=1.3445     diam=24   Δα=2.2065
diam=12   Δα=1.6666     diam=32   Δα=2.4297   (n=8, from an earlier run)
diam=16   Δα=1.9104
```

That is Δα ≈ 0.78·ln(D) − 0.28, fitting every point to within 0.01. Clean
logarithmic divergence, no plateau. Linearity R² stayed 0.78–0.92, so the
power-law fits are sound — it is the *width* that will not settle.

**Why:** Δα is taken over q ∈ [−10, 10]. At |q| = 10 the partition function
Σ p_i^q is dominated by a single box. The smallest possible box is one vertex,
p = 1/N, so α_max inherits a ln(N) term and drifts forever. At extreme q this is
extreme-value statistics of one box, not a scaling measurement.

**The test that decides everything** was `q_range_convergence()` in
`code_12_converged_band.py`, run on 2026-09-14 over q caps {10, 8, 5, 4, 3, 2}
on all four demonstration structures (HKUST-1, MOF5, ZIF-8, UiO-66 — pcu, pcu,
sod and fcu nets respectively, so this is not one topology's quirk):

```
HKUST-1   |q|<=10..2   Δα drift (largest two sizes) = 0.101 -> 0.140   NO CAP CONVERGES
MOF5      |q|<=10..2   Δα drift (largest two sizes) = 0.103 -> 0.138   NO CAP CONVERGES
ZIF-8     |q|<=10..2   Δα drift (largest two sizes) = 0.142 -> 0.194   NO CAP CONVERGES
UiO-66    |q|<=10..2   Δα drift (largest two sizes) = 0.057 -> 0.064   NO CAP CONVERGES
```

Narrowing the q range makes the drift *smaller in absolute terms* (UiO-66 goes
from 0.101 down toward 0.057) but it never drops below the 0.05 convergence
tolerance, and the slope of Δα against ln(diameter) stays 0.7–0.9 — essentially
unchanged — at every cap tested, on every net. There is no q window, narrow or
wide, over which Δα stops tracking supercell size.

**Consequence — this is now a genuine negative result about the published
method (Xiao et al. 2021), not an open question:**

- Every Δα in this project — the 77-framework hMOF band
  (`9_dataset/mof_multifractal_descriptors_v1.csv`), the 61-framework CoRE MOF
  band (`9_dataset/mof_multifractal_descriptors_real_v2.csv`), and a teammate's
  3 sub-bands — was computed at diameter 8–12 and is therefore a readout of the
  supercell rule, not a property of the framework, **at every q range tested**.
- This also explains the r = +0.62 size correlation: that was not a confound
  contaminating the descriptor, it *was* the descriptor.
- **Neither existing band should be presented as a structural finding.** Do
  not recompute them "at the converged cap" — there isn't one. See
  `9_dataset/README.md` and the site's Results page for the writeup; both now
  carry this withdrawal the same way the pore-size and four-structure claims
  were withdrawn — documented, not deleted.
- What would actually fix this: a box-covering measure that is not dominated by
  single-vertex boxes at extreme q (e.g. a fixed, framework-size-independent q
  range chosen a priori, or a different partition-function normalisation), or
  reporting a per-diameter Δα(D) curve — the slope ~0.78 ln(D) itself, rather
  than a single number — as the descriptor. Neither has been implemented.

---

## Defects already found and fixed (do not reintroduce)

1. **τ(q) fitted with a free intercept** instead of through the origin. Flattened
   every spectrum's apex. 0/77 → 77/77 valid. τ must be `Σxy/Σx²`.
2. **Spectrum computed on the induced subgraph of influential nodes.** Metal
   clusters never bond to each other — they connect *through* linkers — so that
   subgraph has zero edges. Box-cover the FULL graph, read the measure only at
   influential nodes.
3. **Averaging order.** Build Z per trial, take logs, then average ln Z. Averaging
   p_i(r) first biases the spectrum narrow by >10×.
4. **Periodicity discarded.** Must use the labelled quotient graph (edges carry
   lattice translation t). Without it UiO-66's 12-connected Zr₆ cluster reports
   as 6-connected or degenerates to a star.
5. **Tercile classes narrower than the error bar.** A teammate's run split Δα into
   Narrow/Medium/Wide at 0.032 wide against ±0.03–0.05 error bars. Her own
   isomorphism check listed hMOF-94 and hMOF-96 as the *same graph*; they landed
   in different classes. Never classify when class width < error bar —
   `assign_classes()` refuses.
6. **nx.diameter() and all-pairs shortest paths** are both O(N²)-ish in pure
   Python and stall above ~3000 vertices. Use `fast_diameter()` (double sweep)
   and the BFS box covering in `code_12`.

## Known traps

- **MOFX-DB provenance.** `database=` is advisory on that endpoint and loses to
  the `gases[]` filter. It has twice returned hMOF records for a CoRE MOF 2019
  request. Filter each returned record on its own `database` field, or work from
  CIFs on disk (`code_12.analyse_directory` does the latter deliberately).
- **Sample is topologically degenerate.** ~48 of 76 hMOFs share
  256 nodes / 384 edges / diameter 12 — a 4×4×4 pcu net, 64 six-connected Zn₄O
  nodes and 192 two-connected linkers. They differ only in linker chemistry,
  which coarse-graining deletes by design. A band over these is narrow for a
  trivial reason. Real sub-band work needs deliberate spread across *nets*
  (pcu, tbo, sod, fcu, rht), not random sampling.
- **Coarse-grain BEFORE replicating.** Building an atomic supercell then
  coarse-graining caps you at MAX_ATOMS and hence diameter ~12. The quotient
  graph replicates combinatorially with no atoms — diameter 35 is cheap.
- **`metal_oxo_deep_dive` is not in git.** It is the local folder containing
  `.git`; the repo root holds `1_website/`, `2_python/` etc. directly. Several
  files still describe a layout with that prefix — those strings are wrong.

## Conventions

- **The website is built, not hand-edited.** Edit `1_website/_build/fragments/*.html`
  or `_build/content_source.html`, then run `python 1_website/_build/build_site.py`.
  Never edit `1_website/index.html` and friends directly.
- **Tests:** `cd 1_website/tests && node test_*.js`. They boot the real pages in
  jsdom and assert on computed output, because this site's failure mode is a
  silently blank panel, not an error.
- **CI** (`.github/workflows/pages.yml`) runs the suites on push. The deploy job
  is manual-only (`workflow_dispatch`) because GitHub Pages is not enabled.
- **Two independent implementations** (Python in `2_python/`, JavaScript engine in
  `1_website/`) were written from the same spec, not from each other. Where they
  disagree, one is wrong — that disagreement has caught real bugs. Keep both.

## House style for this project

Claims are stated with their evidence and their limits. Two claims have already
been withdrawn (the pore-size correlation, and a four-structure real-MOF gap)
and both withdrawals are documented on the site rather than quietly deleted.
Maintain that. Do not describe a motivation as a result. If a number is not
converged, say so.

## Outstanding work

1. ~~Run `q_range_convergence()`~~ — done, 2026-09-14: no q range converges,
   on any of the four demonstration nets. See the resolved section above.
2. ~~Write up the negative result~~ — done: `9_dataset/README.md` and the
   site's Results page (`results.html`) now carry a withdrawal notice for both
   bands. Nothing was recomputed, because there is no converged cap to
   recompute at.
3. Site prose: the overview page (`index.html`) needs rewriting in plainer,
   more formal academic English — current text is over-stylised for a panel.
4. Fix the stale `metal_oxo_deep_dive/` path strings in README and 6 other files.
5. README claims the Python side needs "NumPy. Nothing else" — it needs networkx.
6. **New, from the convergence result:** either implement a size-independent
   Δα (fixed a-priori q range, or reporting the Δα(D) slope itself) and
   recompute both bands under that definition, or retire Δα as this project's
   headline descriptor and say so on the Results and Overview pages.
