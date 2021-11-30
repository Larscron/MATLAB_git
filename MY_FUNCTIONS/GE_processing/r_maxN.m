function [ GE_loc ] = r_maxN( del_z, Om )
%R_MAXN Summary of this function goes here
%   Detailed explanation goes here
r_N=abs(del_z)./(sqrt(diag(Om)));
k = find(r_N == (max(r_N)));
beta=3;
if r_N(k)>beta
    fprintf('The r_max^N test predicts a GE at measurement:\n%i\n',k)
    GE_loc = k; % if a GE is found GE_loc returns its location
else
    fprintf('The r_max^N test DOES NOT predicts any GEs\n')
    GE_loc = 0; % if no a GE is found GE_loc returns 0
end

end

