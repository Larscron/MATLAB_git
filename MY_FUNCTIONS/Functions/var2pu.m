function [puMatrix] = var2pu(varMatrix,S_base)
%VAR2PU converts form var units to pu units
%   this function takes in a susceptance matric in var and the unit base
%   and returns a susceptance matrix in pu
puMatrix=varMatrix/S_base;
end

