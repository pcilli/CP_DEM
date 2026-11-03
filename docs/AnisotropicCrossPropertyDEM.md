# Anisotropic cross-property DEM (aligned spheroids)

This note documents the toolbox's implementation of the manuscript
*Anisotropic cross-property differential effective medium modelling*
(P. Cilli, draft, 2026), the audit of that manuscript that was carried out
before implementing it, and how the new functions relate to the published
isotropic models. The extension lives in its own set of files and does not
alter any existing function.

## 1. What the manuscript proposes

Aligned spheroidal inclusions (semi-axes a, a, c; aspect ratio α = c/a;
symmetry axis x₃) make the effective medium transversely isotropic (TI):
two conductivities (σ∥, σ⊥) and five stiffnesses (C₁₁, C₁₃, C₃₃, C₄₄, C₆₆,
with C₁₂ = C₁₁ − 2C₆₆). The manuscript keeps the DEM analytically closed by
two *closures*, both exact at φ = 0 and in the spherical limit:

- **Electrical (eqs 7–12).** Each component uses the depolarisation
  factor of the spheroid in an *isotropic* host, L along the axis (eq. 3)
  and (1 − L)/2 across it:
  (1 − φ) dσ∥/dφ = (σ₂ − σ∥) σ∥/(σ∥ + L(σ₂ − σ∥)),
  (1 − φ) dσ⊥/dφ = (σ₂ − σ⊥) σ⊥/(σ⊥ + ½(1 − L)(σ₂ − σ⊥)).
- **Elastic (eqs 20–39).** (1 − φ) dC/dφ = (C₂ − C)[I + S_cl C⁻¹(C₂ − C)]⁻¹ in
  engineering-shear Voigt notation, with S_cl the spheroidal Eshelby tensor
  of an isotropic host whose Poisson ratio comes from the Voigt–Reuss–Hill
  average of the current TI stiffness (eqs 27–32); only M₁₁, M₁₃, M₃₃, M₄₄,
  M₆₆ of M = (C₂ − C)A are integrated.
- **Cross-property (eqs 40–52).** Dividing by the electrical equation of
  the measured component eliminates φ: dC_IJ/dσ_p = M_IJ/((σ₂ − σ_p) m̄_p),
  dσ_q/dσ_p = (σ₂ − σ_q) m̄_q/((σ₂ − σ_p) m̄_p). The inverse (section 2.7)
  uses a monotonic stiffness as the independent variable (eq. 63).
- **Partial alignment (eqs 53–57).** A heuristic interpolation in the
  Hermans order parameter S between the aligned factors and the randomly
  oriented ones (electrical: m̄_rand = (2m̄⊥ + m̄∥)/3; elastic: between the
  sphere and the spheroid Eshelby tensors).

## 2. Audit

Everything below was checked independently in Python before the MATLAB
code was written, and the same checks are in `Tests/test_anisotropic_paper.m`
and `python/test_cpdem_aniso.py`.

| Claim | Result |
|---|---|
| Eq. 3 (depolarisation factor, prolate and oblate forms) | Identical to Osborn's factor used by the toolbox to 1e-13 for α from 0.05 to 50 |
| Limits of section 2.1.2 (needle, sphere, disc) | Correct |
| Sphere Eshelby tensor, eqs 25–26 (engineering-shear convention, S₄₄ = 2S₂₃₂₃) | Correct; the spheroidal tensor reduces to it as α → 1 |
| Spheroidal Eshelby tensor (Mura 1987 / Tandon & Weng 1984 forms used here) | The orientation average of the dilute strain-concentration tensor built from it reproduces Berryman's P and Q (cpdem_PQ) to 1e-15 for fluid and solid inclusions, α from 0.05 to 30 |
| VRH reduction, eqs 27–32 | Matches the general Voigt and Reuss definitions to 1e-16; returns the exact moduli for an isotropic stiffness |
| Dry-inclusion closed form, eqs 68–70 | Reproduced to 1e-8 by the ODE |
| Spherical limit "collapses to Berryman/Bruggeman" (sections 2.3.7, 2.8.1) | Exact: K, μ and σ agree with cpdem_ElasDEM and cpdem_ElecDEM to 1e-9 |
| S = 0 electrical branch equals the randomly oriented isotropic DEM (eq. 16, section 2.5) | Exact to 1e-9; the S = 0 elastic branch is the spherical reference, as the manuscript states |
| Section 3 end points (φ = 0.3, S = 1): conductive and pore-like needles, spheres, discs, all five stiffnesses | Every quoted value reproduced to the printed precision (e.g. needles σ∥ 3.536, σ⊥ 1.739, C₃₃ 54.15, C₁₁ 32.39, C₄₄ 14.90, C₆₆ 10.06 GPa; discs σ⊥ 3.173, σ∥ 1.430, C₁₁ 52.54, C₃₃ 9.88) |
| Chain rule (eqs 40–52) and inverse (eq. 63) | Forward from σ∥ or σ⊥ reproduces the φ-parametrised DEMs to 1e-7; the inverse from C₁₁, C₃₃, C₄₄ or C₆₆ round-trips to 1e-6 |

The manuscript is therefore internally correct: its equations follow from
its stated closures, its limits are right, and its numbers are
reproducible. Three remarks for the author:

1. **The electrical closure is unnecessary.** For a coaxial TI host the
   exact one-particle solution of section 4.1 is available in closed form:
   stretching the coordinates by 1/√σᵢ makes the host isotropic and turns
   the spheroid into a spheroid of aspect ratio α' = α√(σ⊥/σ∥), so the
   exact "Hill tensor" is just the depolarisation factors of α'
   (`cpdem_AlignedIncrement`, option `Closure = 'exact'`). At φ = 0.3 the
   closure differs from the exact solution by 2–3 % for conductive
   inclusions and by up to 10 % for resistive aligned discs (σ∥ 0.255
   against 0.285 S/m in the pore-like disc case of section 3).
2. **Near the sphere the closed-form Eshelby components cancel
   catastrophically** (the 1/(α² − 1) terms); at α = 1 ± 10⁻⁶ naive
   evaluation is wrong by an order of magnitude. The toolbox uses an
   eighth-order series about the sphere for |α − 1| < 0.02, derived
   symbolically, and Berryman's θ (itself series-protected) for g.
3. **M = (C₂ − C)A is not exactly symmetric** when S_cl is evaluated for a
   different (isotropic) host than the TI C, so M₁₃ ≠ M₃₁ away from the
   sphere. Integrating M₁₃ (as the manuscript does) is a choice the text
   should state. The manuscript's reference list also gives the
   thermal–electrical paper as published (Materials Letters: X, 12,
   100101), which the README of the main branch still cites as a preprint.

## 3. Functions

| Function | Equations | Purpose |
|---|---|---|
| `cpdem_EshelbySpheroid(nu, alpha)` | 24–26 | 6×6 Eshelby tensor of a spheroid in an isotropic host (engineering Voigt) |
| `cpdem_VRH(C)`, `cpdem_TIStiffness(c5)` | 17–18, 27–32 | TI stiffness assembly and Voigt–Reuss–Hill reduction |
| `cpdem_AnisoElecDEM(phi, sig1, sig2, alpha, S, ...)` | 7–12, 54–56 | TI electrical DEM versus inclusion fraction |
| `cpdem_AnisoElasDEM(phi, k1, mu1, k2, alpha, S, ...)` | 20–39, 57 | TI elastic DEM versus inclusion fraction |
| `cpdem_AnisoForward(k1, mu1, sig1, k2, sig2, sig0, alpha, S, which, ...)` | 40–52 | Stiffnesses and the other conductivity from one measured component |
| `cpdem_AnisoInverse(cMeas, which, k1, mu1, sig1, k2, sig2, alpha, S, ...)` | 63 | Conductivity tensor from one TI stiffness |
| `cpdem_TIVelocities(C11, C13, C33, C44, C66, rho)` | — | Phase velocities and Thomsen parameters |

All accept the usual `cpdem_Options` fields (as a struct or name/value
pairs) plus `Closure = 'isotropic-host'` (the manuscript, default) or
`'exact'` (electrical part exact, S = 1 only). The Python mirror is
`python/cpdem_aniso.py`, cross-checked against Octave through
`python/reference_aniso.json` (written by `Tests/cpdem_WriteAnisoReference.m`).
The example `ExampleScripts/cpdem_AnisotropicCrossProperty_Example.m`
reproduces the end-point table and the trajectories of Figures 1–6.

```matlab
% Section 3 setup of the manuscript: aligned needles, conductive fluid
[sigPar, sigPerp] = cpdem_AnisoElecDEM(0.3, 1, 10, 10);          % 3.536, 1.739 S/m
el = cpdem_AnisoElasDEM(0.3, 40e9, 30e9, 2.2e9, 10);              % el.C33 = 54.15 GPa
f  = cpdem_AnisoForward(40e9, 30e9, 1, 2.2e9, 10, [2 3], 10, 1, 'par');   % stiffnesses from sigma_par
g  = cpdem_AnisoInverse(f.C33, 'C33', 40e9, 30e9, 1, 2.2e9, 10, 10);     % and back
```

## 4. Scope and limitations

The model is the manuscript's closure-based approximation: exact at φ = 0,
exact for spheres, approximate once the medium becomes strongly
anisotropic; the partial-alignment branch is a trend-level surrogate, not
an orientation-distribution average. The fully anisotropic-host framework
of section 4.1 (exact Hill and Eshelby tensors in the evolving TI medium)
is not implemented, except that the electrical half of it is available
through the `'exact'` closure. Inclusions are fluid (μ₂ = 0) with isotropic
conductivity, and the matrix is isotropic, as in the manuscript.
