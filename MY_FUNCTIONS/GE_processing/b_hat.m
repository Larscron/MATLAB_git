function [ measMtx, GE_loc ] = b_hat( del_z, Om, stdev, measMtx )
%B_HAT Summary of this function goes here
%   Detailed explanation goes here
r_N=abs(del_z)./(sqrt(diag(Om)));
k = find(r_N == (max(r_N)));
b_hat=(stdev./(sqrt(diag(Om)))).*r_N;
c=4;
if b_hat(k)>c
    GE_loc = k; % if a GE is found GE_loc returns its location
    measMtx(GE_loc,:)=[];
else
    GE_loc = 0; % if no a GE is found GE_loc returns 0
end

end
%{
function [ measMtx, GE_loc ] = b_hat( del_z, Om, stdev, measMtx )
%B_HAT Summary of this function goes here
%   Detailed explanation goes here
r_N=abs(del_z)./(sqrt(diag(Om)));
k = find(r_N == (max(r_N)));
b_hat=(stdev./(sqrt(diag(Om)))).*r_N;
c=4;
if b_hat(k)>c
    %fprintf('The b^hat test predicts a GE at measurement:\n%i\n',k)
    GE_loc = k; % if a GE is found GE_loc returns its location
    measMtx(GE_loc,:)=[];
    %fprintf('We have deleted measuremetn %i\n',GE_loc)
else
    %fprintf('The b^hat test DOES NOT predicts any GEs\n')
    GE_loc = 0; % if no a GE is found GE_loc returns 0
end

end
%}
