%% Weighted Least Square Calculations
clear all; close all; clc

nbus=14;
pctSTdev=0.5; % should be made an optional function argument
Nrun=10000;

% MATPower Extraction
lineMtx=extr_line_mtx(nbus);
b_bus_sh=extr_bus_sh_vec(nbus);             % this on is empty for this example
measMtx_base=extr_meas_mtx(nbus);

errorLog=zeros(Nrun,2);
for Run=1:Nrun
    clearvars -except nbus pctSTdev lineMtx b_bus_sh Nrun errorLog Run measMtx_base
measMtx=measMtx_base;

% adding normal errors and a randome GE to the measurements
measMtx=add_stdev(measMtx,pctSTdev,3);              % adds a column with the sandard deviation of 1% of each measurment. The measuremnts are located in the 3rd column 
stdev=measMtx(:,6);
measMtx(:,3)=add_nor_error(measMtx(:,3),pctSTdev);  % adds a normal error of proportional to the percentage standard deviation to all measurments
[measMtx(:,3),GE_loc_true]=add_GE(measMtx(:,3),stdev);

% WLS state estimation and GE Processing
GE_found=0;
for it=1:1000 % will loop unitl break or 1000 loops
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

% logging errors
if (GE_found && it==2)
    error=0;
elseif (~GE_found && it==1)
    error=1;
elseif (GE_found && it>2)
    error=2;
elseif (~GE_found && it>=2)
    error=3;
else
    fprintf('How did you even get to this poin?!?!?')
    break
end

errorLog(Run,1)=error;
errorLog(Run,2)=GE_loc_true;

fprintf('Run nr %i done\n',Run)
%{
J index detection
GE = J_ind_det(J, nbus, nmeas, 0.05, GE_loc_true);

r_max^N
GE_loc_r=r_maxN(del_z,Om);

b^hat test
GE_loc_b=b_hat(del_z,Om,stdev);

if GE_loc_b
    measMtx(GE_loc_b,:)=[];
    fprintf('We have deleted measuremetn %i\n',GE_loc_b)
else
    break
end
%}
end

%% display the resutls
errors=find(errorLog(:,1)~=0);
Nerrors=length(errors);
noGEerror=find(errorLog(:,1)==1); % no GE was found
Nno=length(noGEerror);
tmGEerror=find(errorLog(:,1)==2); % too many GE were found, including the true velue
Ntm=length(tmGEerror);
tmnoGEerror=find(errorLog(:,1)==3); % too many GE were found, not including the true velue
Ntmno=length(tmnoGEerror);

pctSuc=100*(1-Nerrors/Nrun);

fprintf('\n\nAfter %i runs there were %i errors (success rate = %.1f%% )\n',Run,Nerrors,pctSuc)
fprintf('Error types and occurrence:\n')
fprintf('- %i No GE was found\n',Nno)
fprintf('- %i Too may GEs were found (including the true value)\n',Ntm)
fprintf('- %i Too may GEs were found (NOT including the true value)\n',Ntmno)



%}