# CP_DEM — Cross-Property DEM Toolbox

[![tests](https://github.com/pcilli/CP_DEM/actions/workflows/tests.yml/badge.svg)](https://github.com/pcilli/CP_DEM/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Paper](https://img.shields.io/badge/DOI-10.1093%2Fgji%2Fggab046-orange.svg)](https://doi.org/10.1093/gji/ggab046)

MATLAB / GNU Octave toolbox (with a Python reference port) for the
**cross-property differential effective medium (DEM) model** of
Cilli and Chapman (2021), which links the elastic moduli and the electrical
conductivity of a porous or composite material through a single geometric
parameter, the **equivalent pore aspect ratio (EPAR)**, with no reference
to porosity. Written by Phil Cilli. Version 1.0 was released on
6 February 2021; version 2.0 is a full rewrite of the numerics with the
same model and a backward-compatible interface.

## What it does

- **Forward model** bulk and shear moduli from electrical conductivity
  (`cpdem_Forward`), for one measurement or a whole data set in a single
  ODE integration, with the porosity the model implies as a by-product.
- **Invert** a bulk or shear modulus for conductivity directly, without a
  search (`cpdem_Inverse`), exact to solver tolerance.
- **Model Vp, Vs and Vp/Vs** from resistivity (`cpdem_ModelVpVs`), with the
  density taken from measurements, from the model's implied porosity, or
  from the Gardner relation.
- **Calibrate** the bulk and shear modulus EPARs to joint
  elastic–electrical data sets, or per sample (`cpdem_CalibrateEPAR`), with
  fit statistics.
- **Bound** the moduli from a conductivity with the electrical–elastic
  Hashin–Shtrikman bounds of Carcione et al. (2007) (`cpdem_ElecElasHS`).
- Evaluate the **single-property DEMs** the model is built from
  (`cpdem_ElecDEM`, `cpdem_ElasDEM`, `cpdem_ImpliedPorosity`), including a
  porosity-dependent aspect ratio such as the power law of Cilli and
  Chapman (2020).
- Link **thermal and electrical conductivity** with the same construction
  (`cpdem_ForwardConduct`; Cilli and Chapman, 2021, arXiv:2102.12434),
  with that paper's three-parameter calibration (`cpdem_CalibrateConduct`),
  its Hashin–Shtrikman bounds (`cpdem_CondCondHS`) and its copper–graphite
  parameters (`cpdem_PublishedCuGraphite`), plus extensions for inclusion
  materials that are themselves anisotropic and for **aligned inclusions**
  (transversely isotropic DEM and its cross-property form,
  `cpdem_AlignedCondDEM`, `cpdem_AlignedForwardConduct`).

Everything runs on MATLAB (R2016b or later, no toolboxes needed) and on
GNU Octave (tested on 8.4), is covered by a test suite that runs in
continuous integration, and is cross-checked against an independent Python
implementation.

## What changed in version 2.0

Version 1.0 was written by hand in January 2021. The rewrite fixes several
problems in it and adds the features above; see [CHANGELOG.md](CHANGELOG.md)
for the full list. The most important points:

- **Accuracy.** v1.0 passed a tolerance to `ode45` as a bare fourth argument,
  which MATLAB silently ignores (and which breaks Octave). The solver
  therefore ran at its default `RelTol = 1e-3`, and because the cross-property
  equations have a `1/σ` singularity at low conductivity the forward model
  carried noise of a few percent and was not even monotonic in conductivity;
  the v1.0 inversions converged to the wrong values on exact synthetic data.
  v2.0 integrates in log conductivity with scaled moduli at `RelTol = 1e-10`
  and round-trips to about 1e-9.
- **Spheres.** EPAR = 1 gave NaN; the shape functions now use a series
  about the sphere (and `P`, `Q` and `m̄` reduce exactly to the classical
  sphere formulae).
- **Speed.** A data set is modelled with one integration instead of one per
  sample, and conductivity is inverted directly instead of by a bounded
  search. Measured in Octave, calibrating both EPARs to 24 samples takes
  about 2 s (v1.0: about half a minute), and inverting them for conductivity
  about 50 ms (v1.0: several seconds).
- **Interface.** The v1.0 functions keep their signatures; `paramFlag`
  is no longer needed (per-sample parameters are detected), invalid flags
  raise errors instead of printing a message, and options are passed in a
  struct (`cpdem_Options`).

## Installation

Clone or download the repository and, in MATLAB or Octave,

```matlab
cd CP_DEM
cpdem_Setup          % adds Functions/, Tests/ and ExampleScripts/ to the path
cpdem_RunTests       % optional: runs the test suite (about 30 s)
```

The Python port needs only `numpy` and `scipy`:

```bash
pip install numpy scipy pytest
pytest -q python
```

## Quick start

```matlab
% Constituent properties (SI units): quartz host, brine inclusions
k1 = 37e9;  mu1 = 44e9;  sig1 = 1e-4;     % mineral: K, mu (Pa), conductivity (S/m)
k2 = 2.25e9; mu2 = 0;    sig2 = 5;        % fluid
sig0 = 1./[2 5 10 20 50];                 % measured conductivities (S/m)

% Published calibration for the sandstones of Han et al. (2011)
e = cpdem_PublishedEPAR('all');           % eparK = 16.4, eparMu = 12.8

% Moduli from conductivity (bulk modulus with alpha_K, shear modulus with alpha_mu)
[K, ~, phi] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2, sig0, e.eparK);
[~, mu]     = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2, sig0, e.eparMu);

% Conductivity from a measured bulk modulus
sigK = cpdem_Inverse(K, 'K', k1,mu1,sig1,k2,mu2,sig2, e.eparK);

% Vp/Vs versus resistivity with a density from the implied porosity
out = cpdem_ModelVpVs(k1,mu1,sig1,k2,mu2,sig2, sig0, e.eparK, e.eparMu, ...
                      struct('method','implied','rho1',2650,'rho2',1050));

% Calibrate both EPARs to your own joint data (K0, mu0, sig0)
cal = cpdem_CalibrateEPAR(K0, mu0, sig0, k1,mu1,sig1,k2,mu2,sig2, 'Bounds', [1 100]);
```

The mineral conductivity `sig1` must be greater than zero: a DEM with a
perfectly insulating host never conducts, however much fluid is added.

## The model in brief

The electrical DEM of a host (phase 1) with randomly oriented spheroidal
inclusions (phase 2) of aspect ratio α, and the elastic DEM of Berryman
(1992), are

$$
(1-\phi)\frac{d\sigma}{d\phi} = (\sigma_2-\sigma)\,\bar m(\sigma,\sigma_2,\alpha),\qquad
(1-\phi)\frac{dK}{d\phi} = (K_2-K)\,P,\qquad
(1-\phi)\frac{d\mu}{d\phi} = (\mu_2-\mu)\,Q ,
$$

with $P$ and $Q$ the polarisation factors of Berryman (1980) and
$\bar m = \tfrac{1}{3}\sigma\sum_{j=1}^{3}\big[\sigma + L_j(\sigma_2-\sigma)\big]^{-1}$
built from the depolarisation factors $L_j$ of the spheroid (Osborn 1945).
Eliminating porosity by the chain rule gives the cross-property DEM
(Cilli and Chapman 2021, equations 18–19; equations 20–21 are the same
with Berryman's $R = \bar m/3\sigma$):

$$
\frac{dK}{d\sigma} = \frac{(K_2-K)\,P}{(\sigma_2-\sigma)\,\bar m},\qquad
\frac{d\mu}{d\sigma} = \frac{(\mu_2-\mu)\,Q}{(\sigma_2-\sigma)\,\bar m},\qquad
K(\sigma_1)=K_1,\ \mu(\sigma_1)=\mu_1 .
$$

The relation depends only on the end-member properties and on α. For real
rocks a different EPAR is needed for each modulus (the bulk modulus EPAR
$\alpha_K$ and the shear modulus EPAR $\alpha_\mu$); prolate ("needle")
EPARs fit sandstones best. The toolbox integrates the system in
$u=\ln(\sigma/\sigma_1)$ with the moduli scaled by $K_1$, carrying porosity
as a third state so that it can be reported and so that α may be a function
of porosity. The theory, conventions, limits and numerical scheme are set
out in [docs/Theory.md](docs/Theory.md).

## What is from the papers and what is not

The core of the toolbox is exactly the model of Cilli and Chapman (2021);
the table says where each piece comes from. Features marked *extension*
are additions of this toolbox that are consistent with the theory but are
not described in the papers; none of them is used unless you ask for it,
and every default reproduces the paper's choice.

| Piece | Source |
|---|---|
| Cross-property equations `dK/dσ`, `dμ/dσ` (`cpdem_Forward`, `cpdem_Deriv`) | Cilli and Chapman (2021), equations 18–19 (m̄ form, default) and 20–21 (R form, `Form = 'R'`), identical |
| $\bar m$ (`mbar`), $R$ (`R`) | Cilli and Chapman (2021); Berryman (1995); depolarisation factors from Osborn (1945) |
| $P$, $Q$, θ, f (`cpdem_PQ`, `cpdem_Shape`) | Berryman (1980, 1992), as used in the paper |
| Separate bulk and shear modulus EPARs, prolate pores, published values | Cilli and Chapman (2021) |
| Calibration by minimising the misfit in one modulus (`cpdem_CalibrateEPAR`, `cpdem_SolveEPAR`) | Cilli and Chapman (2021); the log-space bounded search is a numerical detail |
| Conductivity from a modulus (`cpdem_Inverse`, `cpdem_InvertConduct`), $\sigma_K$ and $\sigma_\mu$ | Cilli and Chapman (2021); solving the same equations with the modulus as the independent variable is a numerical detail |
| Vp/Vs from resistivity with the Gardner (1974) density (`cpdem_ModelVpVs`, default) | Cilli and Chapman (2021) |
| Electrical–elastic Hashin–Shtrikman bounds (`cpdem_ElecElasHS`) | Carcione et al. (2007), as used in the paper |
| Single-property DEMs (`cpdem_ElecDEM`, `cpdem_ElasDEM`) | The two schemes the paper combines; `cpdem_ElecDEM` also evaluates the grain-inclusion form of Cilli and Chapman (2020) and Sen et al. (1981), see its help |
| Porosity-dependent aspect ratio (`epar = @(phi) ...`) | The power-law construction of Cilli and Chapman (2020), with the law supplied by you |
| Thermal from electrical conductivity and the reverse (`cpdem_ForwardConduct`), exact disk/sphere/needle cases | Cilli and Chapman (2021), arXiv:2102.12434, equations 1–6 and 7–9 |
| Three-parameter calibration: aspect ratio over [0, ∞] and effective inclusion conductivities within the flake's principal values, l2 misfit (`cpdem_CalibrateConduct`) | arXiv:2102.12434, Section 3 |
| Thermal–electrical Hashin–Shtrikman bounds (`cpdem_CondCondHS`), published Cu–graphite parameters (`cpdem_PublishedCuGraphite`) | arXiv:2102.12434, Sections 3–4; Carcione et al. (2007) |
| Anisotropic inclusion conductivities, `sig2 = [axis plane]` (`cpdem_CondIncrement`) | *extension*: orientation average of the per-axis concentration factors of an anisotropic ellipsoid, for the graphite-flake case the preprint discusses with bounding isotropic values |
| Aligned inclusions, transversely isotropic DEM and cross-property (`cpdem_AlignedCondDEM`, `cpdem_AlignedForwardConduct`) | *extension*: classical affine-stretch treatment of aligned spheroids in an anisotropic matrix, checked against the exact laminate and fibre limits; to be reconciled with the author's forthcoming anisotropic cross-property DEM |
| Implied porosity output and `cpdem_ImpliedPorosity` | *extension*: a diagnostic; porosity is a dummy variable of the model and the DEM porosity need not equal the true porosity |
| Density from the implied porosity (`dens.method = 'implied'`) | *extension*: not in the paper, which used Gardner's relation |
| Integration in log conductivity, scaled state, tolerances, one integration per data set | numerical scheme only; the equations are unchanged (see `docs/Theory.md`) |

Han et al. (2022) have since assessed the applicability of the
cross-property DEM to reservoir sandstones independently; users modelling
real data should read both that paper and the discussion of clay-bearing
rocks in Cilli and Chapman (2021).

## Functions

| Function | Purpose |
|---|---|
| `cpdem_Forward` | Bulk modulus, shear modulus and implied porosity from conductivity |
| `cpdem_Inverse` | Conductivity from a bulk or shear modulus (direct integration) |
| `cpdem_ModelVpVs` | Vp, Vs and Vp/Vs from conductivity, with three density options |
| `cpdem_CalibrateEPAR` | Calibrate $\alpha_K$ and $\alpha_\mu$ (whole data set or per sample) with fit statistics |
| `cpdem_ForwardConduct` | Cross-property DEM between two conductivity-like properties (thermal from electrical and the reverse); inclusions may be anisotropic |
| `cpdem_CalibrateConduct` | Aspect ratio and effective inclusion conductivities from thermal–electrical data (Section 3 of the preprint) |
| `cpdem_CondCondHS` | Hashin–Shtrikman bounds between two conductivity-like properties |
| `cpdem_PublishedCuGraphite` | Constituent values and solved parameters of the preprint's copper–graphite example |
| `cpdem_AlignedCondDEM`, `cpdem_AlignedForwardConduct` | Transversely isotropic DEM for aligned spheroids and its cross-property form |
| `cpdem_ElecDEM`, `cpdem_ElasDEM` | Single-property electrical and elastic DEMs versus porosity |
| `cpdem_ImpliedPorosity` | Porosity the electrical DEM implies for a conductivity |
| `cpdem_ElecElasHS` | Electrical–elastic Hashin–Shtrikman bounds and the porosity bounds behind them |
| `cpdem_PublishedEPAR` | The EPAR calibrations published by Cilli and Chapman (2021) |
| `cpdem_Options` | Solver tolerances, the written form of the equations (`'mbar'` or `'R'`) and calibration options |
| `cpdem_SolveEPAR`, `cpdem_InvertConduct` | v1.0 interfaces, kept for backward compatibility |

Subfunctions in `Functions/Subfunctions` evaluate the geometry
(`cpdem_Shape`: θ, f and L; `cpdem_PQ`, `P`, `Q`; `mbar`, `R`), the
right-hand side (`cpdem_Deriv`) and the multi-target integration
(`cpdem_Integrate`). All functions have help text (`help cpdem_Forward`);
`help Contents` in the Functions folder lists them.

Parameters may be scalars or per-sample arrays, measurements may be arrays
of any shape, and EPAR may be a scalar, a per-sample array, or a function
handle `@(phi) ...` for a porosity-dependent aspect ratio.

## Published EPARs

Cilli and Chapman (2021) calibrated the model to the brine-saturated
reservoir sandstones of Han et al. (2011) at 60 MPa differential pressure
(`cpdem_PublishedEPAR`):

| Data set | $\alpha_K$ | $\alpha_\mu$ |
|---|---|---|
| All samples | 16.4 | 12.8 |
| Clean sandstones (clay < 3.5 %) | 16.6 | 12.3 |
| Clay-bearing sandstones (clay > 3.5 %) | 16.3 | 12.9 |

EPAR is an effective parameter rather than a literal pore shape; it depends
on the constituent properties used in the calibration and, weakly, on
pressure (Cilli and Chapman, EGU 2022). If you have joint data, calibrate
your own EPARs with `cpdem_CalibrateEPAR` first.

## Examples

The scripts in `ExampleScripts` run out of the box (they use a synthetic
data set generated from the model, clearly labelled as such, in place of
laboratory data that are not redistributed here):

| Script | Shows |
|---|---|
| `cpdem_ModelModuli_Example` | Moduli from resistivity for the three published calibrations, Hashin–Shtrikman bounds, per-sample parameters |
| `cpdem_ModelVpVs_Example` | Vp, Vs and Vp/Vs from resistivity with implied, measured and Gardner densities |
| `cpdem_SolveEPAR_Example` | Calibration of $\alpha_K$ and $\alpha_\mu$, options, per-sample EPARs, the shape of the objective |
| `cpdem_SolveConduct_Example` | Direct inversion of moduli for conductivity, unreachable measurements |
| `cpdem_PorosityAndArchie_Example` | The single-property DEMs, implied porosity, $\bar m$ as a generalised cementation exponent, EPAR(φ) |
| `cpdem_ThermalConduct_Example` | Reproduces the trends and bounds of Figure 1 of the thermal–electrical preprint and its three-parameter calibration (on synthetic data at the paper's volume fractions) |
| `cpdem_AnisotropicConduct_Example` | Graphite-like anisotropic flakes, randomly oriented and aligned; the transversely isotropic cross-property relation |

`cpdem_RunExamples` runs them all (`cpdem_RunExamples(true)` also saves the
figures to `ExampleScripts/output`); set the environment variable
`CPDEM_NOPLOT=1` to skip plotting on a machine without graphics.

## Tests and continuous integration

`Tests/cpdem_RunTests` runs the suite on MATLAB or Octave. It checks the
sphere and Archie limits of the geometry functions, that the cross-property
model equals the composition of the two single-property DEMs, forward /
inverse round trips, EPAR recovery from synthetic data, containment within
the Hashin–Shtrikman bounds, the thermal–electrical model, input validation,
and agreement with v1.0 to the accuracy v1.0 had. The GitHub Actions
workflow runs the suite and the examples on Octave and the Python tests on
every push. `Tests/cpdem_WriteReference` writes the values the Python tests
compare against.

## Python port

`python/cpdem.py` is a function-for-function port on numpy and scipy
(`forward`, `inverse`, `elec_dem`, `elas_dem`, `implied_porosity`,
`forward_conduct`, `elec_elas_hs`, `calibrate_epar`, `model_vpvs`). It
exists mainly as an independent check of the MATLAB code, which it agrees
with to about 1e-8, but it is complete and usable on its own.

## Limitations

The model is derived for a composite that is infinite, linearly elastic
and isotropic, with two isotropic phases, randomly oriented spheroidal
inclusions so dilute (at each DEM step) that they do not interact, and a
static (DC) electric or elastic field applied far from any inclusion. It
has no account of surface or double-layer conduction (Waxman and Smits
1968), which can matter in clay-bearing rocks, as discussed by Cilli and
Chapman (2021). A DEM with a perfectly insulating host never conducts, so
the mineral conductivity must be finite. The aspect ratio is an effective
parameter: when it is given as a function of porosity it must stay
positive at zero porosity, where every DEM path starts.

## Citing

If you use the toolbox, please cite the paper that describes the model
(and, optionally, the software, see [CITATION.cff](CITATION.cff)):

```bibtex
@article{CilliChapman2021,
  author  = {Cilli, Phillip A. and Chapman, Mark},
  title   = {Linking elastic and electrical properties of rocks using cross-property {DEM}},
  journal = {Geophysical Journal International},
  volume  = {225},
  number  = {3},
  pages   = {1812--1823},
  year    = {2021},
  doi     = {10.1093/gji/ggab046}
}
```

## References

- Berryman, J. G. (1980). Long-wavelength propagation in composite elastic media II. Ellipsoidal inclusions. *J. Acoust. Soc. Am.*, 68(6), 1820–1831. doi:10.1121/1.385172
- Berryman, J. G. (1992). Single-scattering approximations for coefficients in Biot's equations of poroelasticity. *J. Acoust. Soc. Am.*, 91(2), 551–571. doi:10.1121/1.402518
- Berryman, J. G. (1995). Mixture theories for rock properties. In *Rock Physics and Phase Relations: A Handbook of Physical Constants*, AGU, 205–228. doi:10.1029/RF003p0205
- Carcione, J. M., Ursin, B., and Nordskag, J. I. (2007). Cross-property relations between electrical conductivity and the seismic velocity of rocks. *Geophysics*, 72(5), E193–E204. doi:10.1190/1.2762224
- Cilli, P. A., and Chapman, M. (2020). The power-law relation between inclusion aspect ratio and porosity: Implications for electrical and elastic modeling. *J. Geophys. Res. Solid Earth*, 125, e2019JB019187. doi:10.1029/2019JB019187
- Cilli, P. A., and Chapman, M. (2021). Linking elastic and electrical properties of rocks using cross-property DEM. *Geophys. J. Int.*, 225(3), 1812–1823. doi:10.1093/gji/ggab046
- Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal conductivity through cross-property inclusion modelling. arXiv:2102.12434
- Mazloum, A., Kováčik, J., Emmer, Š., and Sevostianov, I. (2016). Copper–graphite composites: thermal expansion, thermal and electrical conductivities, and cross-property connections. *J. Mater. Sci.*, 51(17), 7977–7990.
- Norris, A. N. (1985). A differential scheme for the effective moduli of composites. *Mech. Mater.*, 4(1), 1–16.
- Cilli, P. A., and Chapman, M. (2022). On the pressure (in)dependence of cross-property relations. EGU General Assembly 2022, EGU22-9251.
- Eshelby, J. D. (1957). The determination of the elastic field of an ellipsoidal inclusion, and related problems. *Proc. R. Soc. Lond. A*, 241(1226), 376–396.
- Gardner, G. H. F., Gardner, L. W., and Gregory, A. R. (1974). Formation velocity and density — The diagnostic basics for stratigraphic traps. *Geophysics*, 39(6), 770–780.
- Han, T., et al. (2022). Applicability of cross-property differential effective medium model to the joint elastic–electrical properties of reservoir sandstones. *Geophys. Prospect.* doi:10.1111/1365-2478.13233
- Han, T., Best, A. I., Sothcott, J., and MacGregor, L. M. (2011). Joint elastic-electrical properties of reservoir sandstones and their relationships with petrophysical parameters. *Geophys. Prospect.*, 59(3), 518–535. doi:10.1111/j.1365-2478.2010.00940.x
- Mendelson, K. S., and Cohen, M. H. (1982). The effect of grain anisotropy on the electrical properties of sedimentary rocks. *Geophysics*, 47(2), 257–263.
- Osborn, J. A. (1945). Demagnetizing factors of the general ellipsoid. *Phys. Rev.*, 67(11–12), 351–357. doi:10.1103/PhysRev.67.351
- Sen, P. N., Scala, C., and Cohen, M. H. (1981). A self-similar model for sedimentary rocks with application to the dielectric constant of fused glass beads. *Geophysics*, 46(5), 781–795.
- Torquato, S. (2002). *Random Heterogeneous Materials: Microstructure and Macroscopic Properties*. Springer.
- Waxman, M. H., and Smits, L. J. M. (1968). Electrical conductivities in oil-bearing shaly sands. *SPE J.*, 8(2), 107–122.

## Support

Please raise issues, discussion points or suggestions for new features on
this repository's issue tracker, or by email (found on the author's
website).

## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## Disclaimer

Although the author thinks everything in this toolbox is coded correctly,
they take no responsibility for the performance and accuracy of its scripts
and functions, and do not provide any warranties. Ultimately, the toolbox
is used solely at the user's risk.
