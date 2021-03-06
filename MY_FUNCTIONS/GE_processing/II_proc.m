function [ z_corr, GE_loc_II, II ] = II_proc( K, r, z, stdev, Om)
%II_proc( K, r, z, stdev) returns a corrected measurments vector and the
%suspected GE location

I=eye(length(K));                           % Identity matrix of smae dimensions as K
II=sqrt(diag(I)-diag(K))./sqrt(diag(K));    % Innovation Index calculation

CME=r.*sqrt(1+(1./(II.^2)));
CME_N=CME./stdev;
CNE=CME./sqrt(diag(Om));

[val, k] = max(abs(CME_N));                 % locate the largest normal
beta=3;
if abs(val)>beta
    GE_loc_II = k; % if a GE is found GE_loc returns its location
    %fprintf('The II test predicts a GE at measurement:\n%i\n',GE_loc_II)
    z(GE_loc_II)=z(GE_loc_II)-CNE(GE_loc_II)*stdev(GE_loc_II);
    z_corr=z;
    %fprintf('And this measurement has now been corrected\n')
else
    %fprintf('The II test DOES NOT predicts any GEs\n')
    z_corr=z;
    GE_loc_II = 0; % if no a GE is found GE_loc returns 0
end
end

