function [bus_sh_vec] = extr_bus_sh_vec(busNum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mpopt = mpoption('verbose', 0, 'out.all', 0);
switch busNum
    case 30   
        mpc=loadcase(case30); % this is where we declare which case we want to look at
        res=runpf(mpc,mpopt);
    case 14   
        mpc=loadcase(case14); % this is where we declare which case we want to look at
        res=runpf(mpc,mpopt);
    case 9
        mpc=loadcase(case9); % this is where we declare which case we want to look at
        res=runpf(mpc,mpopt);
    case 5
        mpc=loadcase(case5); % this is where we declare which case we want to look at
        res=runpf(mpc,mpopt);
    otherwise
        disp ('We do not have this matrix saved')
        quit
end

bus_sh_vec = res.bus(:,6)./100;

end

