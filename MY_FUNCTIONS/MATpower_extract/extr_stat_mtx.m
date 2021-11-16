function [statMtx] = extr_stat_mtx(busNum)

switch busNum
    case 30   
        mpc=loadcase(case30); % this is where we declare which case we want to look at
        res=runpf(mpc);
    case 14   
        mpc=loadcase(case14); % this is where we declare which case we want to look at
        res=runpf(mpc);
    case 9
        mpc=loadcase(case9); % this is where we declare which case we want to look at
        res=runpf(mpc);
    case 5
        mpc=loadcase(case5); % this is where we declare which case we want to look at
        res=runpf(mpc);
        
    otherwise
        disp ('We do not have this matrix saved')
        quit
end

V=res.bus(:,8);
thet_deg=res.bus(:,9);
thet_rad=thet_deg.*(pi/180);

statMtx=[V,thet_rad];
    
end

