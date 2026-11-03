function C = cpdem_TIStiffness(c5)
% Build the 6x6 Voigt stiffness matrix of a transversely isotropic medium
%
%   C = cpdem_TIStiffness([C11 C13 C33 C44 C66])
%
% Symmetry axis along x3; C12 = C11 - 2 C66. Voigt order 11, 22, 33, 23,
% 13, 12 with engineering shear strains, the single convention used by the
% anisotropic cross-property functions (Cilli, "Anisotropic cross-property
% differential effective medium modelling", eqs 17-18).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

C11 = c5(1); C13 = c5(2); C33 = c5(3); C44 = c5(4); C66 = c5(5);
C12 = C11 - 2*C66;
C = [C11 C12 C13 0   0   0
     C12 C11 C13 0   0   0
     C13 C13 C33 0   0   0
     0   0   0   C44 0   0
     0   0   0   0   C44 0
     0   0   0   0   0   C66];

end
