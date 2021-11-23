function [ GE ] = J_ind_det( J, nbus, nmeas, alpha, GE_loc )
%J_IND_DET: [ GE ] = J_ind_det( J, nbus, nmeas, alpha, GE_loc )
%   Detailed explanation goes here
% J index detection
N=nbus;                     % number of state variables
m=nmeas;                    % number of measurements
dof=m-N;                    % number of degrees of freedom
chi2=chi2cdf(J,dof);        % prcentage chnace J should be this or lower
if chi2>(1-alpha)
    disp("A GE IS SUSPECTED")
    GE=1;
else
    disp("NO GE IS SUSPECTED")
    GE=0;
    if (exist ('GE_loc','var'))
        fprintf('However, there is a GE at measurement\n%i\n',GE_loc)
    end
end

end

