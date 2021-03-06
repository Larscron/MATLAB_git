%% Weighted Least Square Calculations
clear all; close all; clc

nbus=14;
pctSTdev=0.5; % should be made an optional function argument

% MATPower Extraction
lineMtx=extr_line_mtx(nbus);
b_bus_sh=extr_bus_sh_vec(nbus);             % this on is empty for this example
measMtx_base=extr_meas_mtx(nbus);

measMtx=measMtx_base;

% adding normal errors and a randome GE to the measurements
measMtx=add_stdev(measMtx,pctSTdev,3);              % adds a column with the sandard deviation of 1% of each measurment. The measuremnts are located in the 3rd column 
stdev=measMtx(:,6);
measMtx(:,3)=add_nor_error(measMtx(:,3),pctSTdev);  % adds a normal error of proportional to the percentage standard deviation to all measurments
[measMtx(:,3),GE_loc_true]=add_GE(measMtx(:,3),stdev);

% WLS state estimation and GE Processing
%WLSres=WLS_function( nbus, measMtx, lineMtx, b_bus_sh);


for it=1:1000 % will loop unitl break or 1000 loops
% WLS state estimation and GE Processing
WLSres=WLS_function( nbus, measMtx, lineMtx, b_bus_sh);

H=WLSres.H; G=WLSres.G; R=WLSres.R; J=WLSres.J;
del_z=WLSres.del_z; del_x_hat=WLSres.del_x_hat;

nmeas=WLSres.Nmeas; % not really used....

% GE analysis
K=H*inv(G)*transpose(H)*inv(R);    % The projection matrix
I=eye(length(K));                   % Identity matrix of smae dimensions as K
S=I-K;                              % Residual sensitivity matrix
Om=S*R;                             % Covariance matrix
del_z_hat=H*del_x_hat;
r=del_z-del_z_hat;                  % residuals

% Innovation Index
[measMtx(:,3), GE_loc_guess, II] = II_proc( K, r, measMtx(:,3), stdev, Om);

% b^hat test
%stdev=measMtx(:,6);
%[measMtx, GE_loc_guess]=b_hat(del_z,Om,stdev,measMtx);

% r_max^N
%[measMtx, GE_loc_guess]=r_maxN(del_z,Om,measMtx);

if GE_loc_guess==GE_loc_true
    GE_found=1;
end
if GE_loc_guess==0
    break
end
end
%}