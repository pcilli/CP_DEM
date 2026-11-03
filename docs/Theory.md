# Theory and numerical notes

This note collects the equations implemented in the Cross-Property DEM
Toolbox, the conventions behind them, the limits used in the tests, and the
numerical scheme of version 2.0. The model itself is that of Cilli and
Chapman (2021); nothing here changes it.

## 1. Single-property DEMs

Consider a host (phase 1) into which randomly oriented spheroidal inclusions
(phase 2) of aspect ratio α are added incrementally, so that at every step
the inclusions are dilute in the current effective medium (the differential
effective medium, DEM, scheme). With φ the inclusion volume fraction
(porosity for a rock), the effective properties obey

$$
(1-\phi)\frac{d\sigma}{d\phi} = (\sigma_2-\sigma)\,\bar m(\sigma,\sigma_2,\alpha),
\qquad \sigma(0)=\sigma_1,
$$

for a conductivity-type property (electrical conductivity σ, but equally
thermal conductivity, permittivity or permeability, since each obeys
Laplace's equation), and

$$
(1-\phi)\frac{dK}{d\phi} = (K_2-K)\,P(K,\mu,K_2,\mu_2,\alpha),\qquad
(1-\phi)\frac{d\mu}{d\phi} = (\mu_2-\mu)\,Q(K,\mu,K_2,\mu_2,\alpha),
\qquad K(0)=K_1,\ \mu(0)=\mu_1,
$$

for the elastic moduli (Berryman 1992).

### Aspect ratio convention

Throughout the toolbox **EPAR** is the ratio of the length of the unique
(symmetry) axis of the spheroid to the length of the two equal axes:
EPAR < 1 is an oblate spheroid (disk, penny), EPAR = 1 a sphere, EPAR > 1 a
prolate spheroid (needle). Berryman (1980) instead writes α = a/c with
a > b = c for prolate and a = b < c for oblate spheroids so that his α is in
(0, 1] for both; Osborn (1945) uses m = (long axis)/(short axis) ≥ 1.
`cpdem_Shape` converts from EPAR to each.

### Geometric functions

For the elastic P and Q, Berryman (1980) defines, with α = EPAR,

$$
\theta = \begin{cases}
\dfrac{\alpha}{(1-\alpha^2)^{3/2}}\left[\arccos\alpha - \alpha\sqrt{1-\alpha^2}\right], & \alpha<1,\\[2ex]
\dfrac{\alpha}{(\alpha^2-1)^{3/2}}\left[\alpha\sqrt{\alpha^2-1} - \operatorname{arccosh}\alpha\right], & \alpha>1,
\end{cases}
\qquad
f = \frac{\alpha^2}{1-\alpha^2}\,(3\theta-2),
$$

and the nine functions $F_1 \ldots F_9$ of $A=\mu_2/\mu-1$,
$B=\tfrac13(K_2/K-\mu_2/\mu)$ and $R=3\mu/(3K+4\mu)$ listed in the appendix
of Berryman (1980) and in Mavko et al. (2009), from which

$$
P = \frac{F_1}{F_2},\qquad
Q = \frac15\left[\frac{2}{F_3}+\frac{1}{F_4}+\frac{F_4F_5+F_6F_7-F_8F_9}{F_2F_4}\right].
$$

For the electrical $\bar m$, the depolarisation factor along the unique
axis of the spheroid is (Osborn 1945; m = 1/α for oblate, m = α for prolate)

$$
L = \begin{cases}
\dfrac{m^2}{m^2-1}\left[1-\dfrac{\arcsin\!\big(\sqrt{m^2-1}/m\big)}{\sqrt{m^2-1}}\right], & \text{oblate},\\[2ex]
\dfrac{1}{m^2-1}\left[\dfrac{m}{2\sqrt{m^2-1}}\ln\dfrac{m+\sqrt{m^2-1}}{m-\sqrt{m^2-1}}-1\right], & \text{prolate},
\end{cases}
$$

with the two degenerate axes sharing $L_2=L_3=(1-L)/2$ so that
$\sum_j L_j = 1$. One can verify that $L = 1-\theta$ identically (Berryman's
θ is the depolarisation factor of the degenerate axes summed, $2L_2$); the
toolbox keeps both expressions and tests the identity.

Then (Cilli and Chapman 2021)

$$
\bar m(\sigma,\sigma_2,\alpha) = \frac{\sigma}{3}\sum_{j=1}^{3}\frac{1}{\sigma+L_j(\sigma_2-\sigma)}
= \frac{\sigma}{3}\left[\frac{1}{\sigma-L(\sigma-\sigma_2)}+\frac{4}{\sigma+\sigma_2+L(\sigma-\sigma_2)}\right],
$$

and Berryman's (1995) $R = \bar m/(3\sigma)$, in terms of which the
electrical DEM reads $(1-\phi)\,d\sigma/d\phi = 3\sigma(\sigma_2-\sigma)R$.
The two forms are identical; the toolbox provides both (`mbar`, `R`).

### Limits used in the tests

- **Spheres** (α → 1): θ → 2/3, f → −2/5, L → 1/3;
  $P \to (K+\tfrac43\mu)/(K_2+\tfrac43\mu)$,
  $Q \to (\mu+\zeta)/(\mu_2+\zeta)$ with $\zeta=\tfrac{\mu}{6}\,\tfrac{9K+8\mu}{K+2\mu}$;
  $\bar m \to 3\sigma/(2\sigma+\sigma_2)$, so the electrical DEM becomes
  Bruggeman's asymmetric model with the closed form
  $\dfrac{\sigma-\sigma_2}{\sigma_1-\sigma_2}\left(\dfrac{\sigma_1}{\sigma}\right)^{1/3} = 1-\phi$.
- **Insulating inclusions** ($\sigma_2 = 0$): $\bar m = \tfrac13\sum_j (1-L_j)^{-1}$
  is independent of σ, the DEM integrates to $\sigma = \sigma_1(1-\phi)^{\bar m}$,
  and with the fluid as host this is Archie's law with the aspect-ratio
  dependent cementation exponent of Sen et al. (1981) and Mendelson and
  Cohen (1982): 3/2 for spheres, larger for disks, 1.66 for needles of
  EPAR 16.4. For a conducting fluid $\bar m$ depends on σ and the apparent
  exponent varies with porosity; letting α itself vary with porosity (Cilli
  and Chapman 2020) generalises Archie's law further.

## 2. Cross-property DEM

Porosity appears in both single-property DEMs only through the common factor
$d\phi/(1-\phi)$. Dividing the elastic equations by the electrical one
eliminates it (Cilli and Chapman 2021, equations 18 and 19; equations 20
and 21 are the same written with R):

$$
\frac{dK}{d\sigma} = \frac{(K_2-K)\,P}{(\sigma_2-\sigma)\,\bar m},\qquad
\frac{d\mu}{d\sigma} = \frac{(\mu_2-\mu)\,Q}{(\sigma_2-\sigma)\,\bar m},\qquad
K(\sigma_1)=K_1,\ \mu(\sigma_1)=\mu_1 .
$$

Given the end-member properties and α, the moduli are a function of
conductivity alone. Porosity has become a dummy variable, but it can be
carried along,

$$
\frac{d\phi}{d\sigma} = \frac{1-\phi}{(\sigma_2-\sigma)\,\bar m},\qquad \phi(\sigma_1)=0,
$$

which the toolbox does for three reasons: it reports the porosity the model
implies for each measurement (`cpdem_ImpliedPorosity`, third output of
`cpdem_Forward`), it gives a consistent bulk density
$\rho=(1-\phi)\rho_1+\phi\rho_2$ for velocities (`cpdem_ModelVpVs`), and it
allows α to be a function of porosity along the path (`cpdem_EvalEPAR`).

The same construction links any two DEM-modelled properties. For two
conductivity-type properties, e.g. thermal conductivity λ and electrical
conductivity σ (Cilli and Chapman 2021, arXiv:2102.12434),

$$
\frac{d\lambda}{d\sigma} = \frac{(\lambda_2-\lambda)\,\bar m(\lambda,\lambda_2,\alpha)}{(\sigma_2-\sigma)\,\bar m(\sigma,\sigma_2,\alpha)},
$$

implemented in `cpdem_ForwardConduct`.

### Two EPARs

Theoretical inclusion modelling uses one aspect ratio. For real rocks Cilli
and Chapman (2021) found that the aspect ratio minimising the misfit in the
bulk modulus ($\alpha_K$) differs from that minimising the misfit in the
shear modulus ($\alpha_\mu$), and that only the modulus whose EPAR was used
is predicted accurately. The toolbox therefore treats them separately:
`cpdem_CalibrateEPAR` calibrates both, `cpdem_ModelVpVs` uses $\alpha_K$ for
K and $\alpha_\mu$ for μ, and `cpdem_Inverse` returns $\sigma_K$ from the
bulk modulus and $\sigma_\mu$ from the shear modulus. The published values
for the sandstones of Han et al. (2011) are in `cpdem_PublishedEPAR`.

### Direct inversion

Because K (and μ) is a strictly monotonic function of σ along a DEM path,
the relation can be inverted without a search by taking the modulus as the
independent variable:

$$
\frac{d\sigma}{dK} = \frac{(\sigma_2-\sigma)\,\bar m}{(K_2-K)\,P},\qquad
\frac{d\mu}{dK} = \frac{(\mu_2-\mu)\,Q}{(K_2-K)\,P},\qquad
\frac{d\phi}{dK} = \frac{1-\phi}{(K_2-K)\,P},
$$

integrated from $K=K_1$ to the measured modulus (`cpdem_Inverse`; the shear
modulus case is analogous). A measured modulus outside the interval between
the inclusion and host moduli has no solution and returns NaN.

### Hashin–Shtrikman bounds

Carcione et al. (2007) obtain porosity bounds by inverting the
Hashin–Shtrikman conductivity bounds for the measured σ and then evaluate the
elastic Hashin–Shtrikman bounds at those porosities. With an inviscid fluid
the lower shear bound is zero and the lower bulk bound is the Reuss average.
Inverting the *upper* conductivity bound gives the *minimum* porosity, at
which the *upper* elastic bounds are evaluated, and vice versa; any
cross-property DEM prediction must lie between them (`cpdem_ElecElasHS`,
checked in the tests for oblate, spherical and prolate EPARs).

## 3. Numerical scheme (version 2.0)

### Why v1.0 was noisy

At low conductivity $\sigma \ll \sigma_2$, $\bar m \approx \tfrac{\sigma}{3\sigma_2}\sum_j L_j^{-1}$
is proportional to σ, so the right-hand sides of the cross-property
equations behave as $1/\sigma$: the moduli change linearly in $\ln\sigma$. A
fixed relative tolerance in σ then spends its error budget unevenly, and
v1.0's call `ode45(f, [sig1; sig0], y0, 1e-5)` did not set a tolerance at
all (MATLAB ignores a non-struct fourth argument; Octave passes it to the
derivative function and fails). Run at the default `RelTol = 1e-3` in Octave,
the v1.0 forward model showed non-monotonic scatter of up to about 3 % in
the moduli between neighbouring conductivities, and the bounded searches of
`cpdem_InvertConduct` and `cpdem_SolveEPAR` converged to wrong values on
exact synthetic data (e.g. 0.0158 S/m for a true 0.02 S/m).

### Scaled formulation

Version 2.0 integrates every path in scaled variables:

- independent variable $u=\ln(\sigma/\sigma_1)$ (for modulus-driven paths,
  $x=(m-m_1)/K_1$; for porosity-driven paths, φ itself);
- state $[K/K_1,\ \mu/K_1,\ \phi]$ (or $[\ln(\sigma/\sigma_1),\ \cdot,\ \phi]$
  for the inverse), all O(1), so that `AbsTol` is meaningful;
- `ode45` with `RelTol = 1e-10`, `AbsTol = 1e-13` by default
  (`cpdem_Options`).

All targets (e.g. every measured conductivity of a data set) are passed to
the solver as its output grid, so a data set with shared parameters costs
one integration (`cpdem_Integrate`), which also makes the modelled values
mutually consistent (no step-sequence noise between samples).

Measured accuracy in Octave (`test_forward_composition`, `test_inverse`):
the cross-property model agrees with the composition of the two
single-property DEMs to about 1e-9, forward/inverse round trips to about
1e-9, and the independent Python port (scipy's DOP853) to about 1e-8.

### Near the sphere

The closed forms for θ, f and L suffer catastrophic cancellation as α → 1
and are 0/0 at α = 1. In double precision f is only accurate to about
1e-8 at |α − 1| = 1e-3 and to about 1e-12 at |α − 1| = 1e-2. For
|α − 1| < 0.1 the toolbox therefore uses the Taylor series about the sphere,
$\theta = \tfrac23 + \tfrac{4}{15}e - \tfrac{6}{35}e^2 + \tfrac{32}{315}e^3 - \ldots$,
$f = -\tfrac25 - \tfrac{12}{35}e + \tfrac{2}{15}e^2 - \tfrac{32}{1155}e^3 - \ldots$,
$L = 1-\theta$, with $e=\alpha-1$, to 13th order (exact rational
coefficients are listed in `cpdem_Shape.m`); its truncation error is below
1e-16 on that interval and the two branches agree to 1e-14 at the switch.

### EPAR as a function of porosity

With α = α(φ), the DEM equations are integrated with α evaluated at the
current porosity state. Every path starts at φ = 0, so the law must be
positive there: a pure power law $a\phi^b$ (Cilli and Chapman 2020) is
zero at φ = 0 and must be given a floor, `@(phi) max(a*phi.^b, 1e-3)`.
Solver trial states can lie marginally outside [0, 1); the law is evaluated
at the nearest physical porosity.

## 4. Assumptions and limitations

The model is derived for a composite that is infinitely extended, linearly
elastic and isotropic, with two isotropic phases, spheroidal inclusions
with at least two equal semi-axes that are randomly oriented, a static
(DC) electrical or elastic source infinitely far from any inclusion, and
inclusions so dilute at each DEM step that their fields do not interact. It
has no account of surface or double-layer conduction (Waxman and Smits
1968), which can be significant in clay-bearing rocks (discussed by Cilli
and Chapman 2021). A DEM with a perfectly insulating host never conducts, so
the host conductivity must be finite. EPAR is an effective parameter, not a
literal pore shape; it depends on the constituent properties used in the
calibration and varies weakly with pressure (Cilli and Chapman, EGU 2022).

## References

See the list in the [README](../README.md).
