function [measMtx] = extr_meas_mtx(busNum)
%EXTR_MEAS_MTX Summary of this function goes here
%   Detailed explanation goes here
switch busNum
    case 14   
        mpc=loadcase(case14); % this is where we declare which case we want to look at
        res=runpf(mpc);
        
    otherwise

measMtx = inputArg1;

end

