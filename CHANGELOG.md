# Changelog

## 2.0.0 — 2026-11-03

A rewrite of the numerics with the same model and a backward-compatible
interface. Tested on GNU Octave 8.4 and written for MATLAB R2016b or later
(no toolboxes required); v1.0 had only been tested on MATLAB R2019b.

### Fixed

- **Solver tolerance was never applied.** `cpdem_Forward` passed `1e-5` to
  `ode45` as a bare fourth argument. MATLAB ignores a non-struct there, so
  the solver ran at its default `RelTol = 1e-3`; Octave passes it on to the
  derivative function and fails. Because the right-hand side behaves as
  `1/σ` at low conductivity (where `m̄ ∝ σ`), the v1.0 forward model carried
  solver noise of a few percent and was not monotonic in conductivity, and
  `cpdem_InvertConduct`/`cpdem_SolveEPAR` consequently converged to the
  wrong values even on exact synthetic data. The system is now integrated in
  `ln σ` with the moduli scaled by the background bulk modulus, at
  `RelTol = 1e-10`, `AbsTol = 1e-13` (configurable through `cpdem_Options`);
  round trips agree to about 1e-9.
- **Spheres (EPAR = 1) returned NaN** from `P`, `Q`, `mbar` and `R`
  (0/0 in the closed forms), and `P` used `epar < 1` where `Q` used
  `epar <= 1`. A shared `cpdem_Shape` now evaluates θ, f and L with a
  13th-order series about the sphere for |EPAR − 1| < 0.1, where the closed
  forms lose accuracy, and the sphere limits reduce exactly to the classical
  formulae.
- Invalid `modFlag`, `normFlag` and `paramFlag` values printed a message
  and then crashed on an undefined variable; they now raise errors.
- `cpdem_InvertConduct` did not preallocate its output; the cost functions
  called `clearvars` on every evaluation.
- Example scripts began with `clear all` and referred to variables they
  never defined, so they could not run as written.

### Added

- `cpdem_Inverse`: direct inversion of a bulk or shear modulus for
  conductivity by integrating the cross-property equations with the modulus
  as the independent variable (no search; about 60× faster than v1.0 in
  Octave and exact to solver tolerance). `cpdem_InvertConduct` now wraps it.
- `cpdem_CalibrateEPAR`: calibrates α_K and α_μ together, in log(EPAR) by
  default, with l1/l2, relative and weighted misfits, per-sample EPARs and
  fit statistics (RMS, MAPE, R²). `cpdem_SolveEPAR` now wraps it. One ODE
  integration per trial EPAR models the whole data set (about 15× faster
  than v1.0 for 24 samples in Octave).
- `cpdem_Forward` accepts arrays of measurements of any shape (one
  integration when the parameters are shared), per-sample parameters without
  a flag, per-sample EPARs, an options struct, and returns the implied
  porosity as a third output. Measurements outside (σ₁, σ₂) give NaN with a
  warning (`CPDEM:outOfRange`).
- EPAR may be a function handle of porosity, e.g. the power law of Cilli
  and Chapman (2020), in every forward, inverse and single-property routine.
- `cpdem_ElecDEM`, `cpdem_ElasDEM`, `cpdem_ImpliedPorosity`: the
  single-property DEMs and the porosity behind a conductivity.
- `cpdem_ForwardConduct`: cross-property DEM between two conductivity-like
  properties (thermal from electrical; Cilli and Chapman 2021,
  arXiv:2102.12434), invertible by swapping the roles of the properties.
- Thermal–electrical model of Cilli and Chapman (arXiv:2102.12434) in
  full: `cpdem_ForwardConduct` (equations 5–6, exact disk/sphere/needle
  cases 7–9 with EPAR 0 and Inf accepted), `cpdem_CalibrateConduct`
  (Section 3: aspect ratio over [0, ∞] and bounded effective inclusion
  conductivities by l2 misfit, found by a profile search over aspect ratio
  that is robust to the strong parameter trade-off of this problem),
  `cpdem_CondCondHS` (thermal–electrical Hashin–Shtrikman
  bounds by Carcione's stratagem) and `cpdem_PublishedCuGraphite` (the
  paper's constituent values and solved parameters).
- Anisotropic inclusion conductivities (`sig2 = [axis plane]`,
  `cpdem_CondIncrement`) in the conductivity DEM, implied porosity and the
  thermal–electrical model, for inclusions such as graphite flakes.
- `cpdem_AlignedCondDEM`, `cpdem_AlignedForwardConduct`: transversely
  isotropic DEM for aligned spheroids (affine stretch of the host) and its
  cross-property form, tested against the exact laminate, fibre and sphere
  limits.
- `cpdem_ModelVpVs`: Vp, Vs and Vp/Vs from conductivity with measured,
  implied-porosity or Gardner densities (v1.0 had only Gardner, in a script).
- `cpdem_ElecElasHS`: input validation and the porosity bounds as extra
  outputs.
- `cpdem_PublishedEPAR`: the EPAR table of Cilli and Chapman (2021).
- `cpdem_Options`, `cpdem_Setup`, `cpdem_Version`, `Functions/Contents.m`.
- Test suite (`Tests/cpdem_RunTests`, 10 test files) and GitHub Actions
  workflow running the tests and examples on Octave and the Python tests.
- Python reference port (`python/cpdem.py`) with tests that cross-check it
  against values written by the Octave toolbox.
- Six self-contained example scripts with a synthetic data generator and a
  runner (`cpdem_RunExamples`), and `docs/Theory.md`.

### Changed

- `cpdem_Deriv` has a new signature (scaled state `[K/K₁; μ/K₁; φ]`,
  independent variable `ln(σ/σ₁)`); the old `cpdem_SolveEPARCost` and
  `cpdem_InvertConductCost` subfunctions were removed.
- `cpdem_InvertConduct` clips its result to `[lb, ub]` and returns the
  nearest bound for moduli the model cannot reach (what the v1.0 search did
  in effect); `normFlag` and `paramFlag` are accepted but no longer needed.
- `P.m`, `Q.m` are thin wrappers around `cpdem_PQ`; `mbar` and `R` are
  vectorised over conductivity.

## 1.0 — 2021-02-06

Initial release accompanying Cilli and Chapman (2021); see
`ReleaseNotes_CP_DEM.pdf`.
