% Cross-Property DEM Toolbox (CP_DEM)
% Version 2.0.0
%
% MATLAB/GNU Octave implementation of the cross-property differential
% effective medium (DEM) model of Cilli and Chapman (2021), which links the
% elastic moduli and electrical conductivity of a two-phase composite
% (e.g. a fluid-saturated rock) through the equivalent pore aspect ratio
% (EPAR) without reference to porosity.
%
% Cross-property modelling (elastic <-> electrical)
%   cpdem_Forward         - Elastic moduli (and implied porosity) from conductivity
%   cpdem_Inverse         - Conductivity from a bulk or shear modulus (direct ODE)
%   cpdem_ModelVpVs       - Vp, Vs and Vp/Vs from conductivity
%   cpdem_CalibrateEPAR   - Calibrate alpha_K and alpha_mu to joint data sets
%   cpdem_ForwardConduct  - Cross-property DEM between two conductivity-like
%                           properties (e.g. thermal from electrical)
%   cpdem_CalibrateConduct - Aspect ratio and effective inclusion conductivities
%                           from thermal-electrical data (arXiv:2102.12434, Sec. 3)
%   cpdem_CondCondHS      - Hashin-Shtrikman bounds between two conductivities
%   cpdem_PublishedCuGraphite - Parameters of the copper-graphite example
%   cpdem_AlignedCondDEM  - Transversely isotropic DEM for aligned spheroids
%   cpdem_AlignedForwardConduct - Its cross-property form
%
% Single-property DEM and diagnostics
%   cpdem_ElecDEM         - Electrical DEM: conductivity vs porosity
%   cpdem_ElasDEM         - Elastic DEM (Berryman 1992): moduli vs porosity
%   cpdem_ImpliedPorosity - Porosity the electrical DEM implies for a conductivity
%   cpdem_ElecElasHS      - Electrical-elastic Hashin-Shtrikman bounds (Carcione et al. 2007)
%
% v1.0 interfaces (kept for backward compatibility)
%   cpdem_SolveEPAR       - Solve for one EPAR (wraps cpdem_CalibrateEPAR)
%   cpdem_InvertConduct   - Invert for conductivity (wraps cpdem_Inverse)
%
% Options and utilities
%   cpdem_Options         - Solver and calibration options
%   cpdem_Version         - Toolbox version string
%
% Subfunctions (Functions/Subfunctions)
%   cpdem_Shape           - theta, f (Berryman 1980) and depolarisation L (Osborn 1945)
%   cpdem_PQ, P, Q        - Berryman's (1980) P and Q polarisation factors
%   mbar, R               - Electrical geometric factors (Cilli & Chapman 2021; Berryman 1995)
%   cpdem_Deriv           - Right-hand side of the cross-property equations
%   cpdem_Integrate       - Multi-target ODE integration helper
%   cpdem_EvalEPAR        - Constant or porosity-dependent EPAR
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
