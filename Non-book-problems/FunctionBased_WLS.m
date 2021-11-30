%% Weighted Least Square Calculations
clear all; close all; clc

nbus=14;
pctSTdev=1; % should be made an optional function argument

measMtx=extr_meas_mtx(nbus);
lineMtx=extr_line_mtx(nbus);
b_bus_sh=extr_bus_sh_vec(nbus);             % this on is empty for this example

% adding normal errors and a randome GE to the measurements
measMtx=add_stdev(measMtx,pctSTdev,3);              % adds a column with the sandard deviation of 1% of each measurment. The measuremnts are located in the 3rd column 
stdev=measMtx(:,6);
measMtx(:,3)=add_nor_error(measMtx(:,3),pctSTdev);  % adds a normal error of proportional to the percentage standard deviation to all measurments
[measMtx(:,3),GE_loc_true]=add_GE(measMtx(:,3),stdev);

for i=1:20
% All  WLS state estimation is done inside this function 
WLSres=WLS_function( nbus, measMtx, lineMtx, b_bus_sh);

H=WLSres.H; G=WLSres.G; R=WLSres.R; J=WLSres.J;
del_z=WLSres.del_z; del_x_hat=WLSres.del_x_hat;

stdev=WLSres.STdev;

nmeas=WLSres.Nmeas;

% This is where GE analysis will ocure
K=H*inv(G)*transpose(H)*inv(R);    % The projection matrix
I=eye(length(K));                   % Identity matrix of smae dimensions as K

S=I-K;                              % Residual sensitivity matrix
Om=S*R;                             % Covariance matrix

del_z_hat=H*del_x_hat;
r=del_z-del_z_hat; % residuals

% J index detection
%GE = J_ind_det(J, nbus, nmeas, 0.05, GE_loc_true);

% r_max^N
%GE_loc_r=r_maxN(del_z,Om);

% b^hat test
GE_loc_b=b_hat(del_z,Om,stdev);

if GE_loc_b
    measMtx(GE_loc_b,:)=[];
    fprintf('We have deleted measuremetn %i\n',GE_loc_b)
else
    break
end

end

%{
% Innovation Index
II=sqrt(diag(I)-diag(K))./sqrt(diag(K));
% II=zeros(nmeas,1); % for i=1:nmeas %     II(i)=sqrt(K(i,i))/sqrt(I(i,i)-K(i,i)); % end
e_D=r;  % residuals are the detectable errors
e_U=(1./II).*e_D;
CME=r.*sqrt(1+(1./(II.^2)));
CME_N=CME./stdev;
% k = find(CME_N == (max(abs(CME_N))));
[val, k] = max(abs(CME_N));
beta=3;
if abs(val)>beta
    fprintf('The II test predicts a GE at measurement:\n%i\n',k)
    GE_loc_II = k; % if a GE is found GE_loc returns its location
else
    fprintf('The II test DOES NOT predicts any GEs\n')
    GE_loc_II = 0; % if no a GE is found GE_loc returns 0
end
%}