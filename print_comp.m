function [] = print_comp(V_mag,V_ang)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nbus=length(V_mag);
switch nbus
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

V_ang_deg = V_ang*180/pi;

disp('---------My Result--------   -----MATPOWER Result------');
disp('| Bus |    V   |  Angle  |   | Bus |    V   |  Angle  |');
disp('| No  |   pu   |  Degree |   | No  |   pu   |  Degree |');
disp('--------------------------   --------------------------');
for m = 1:nbus
    fprintf('%4g', m); fprintf('  %8.4f', V_mag(m)); fprintf('   %8.4f', V_ang_deg(m)); fprintf('  |');
    fprintf('%4g', m); fprintf('  %8.4f', res.bus(m,8)); fprintf('   %8.4f', res.bus(m,9)); fprintf('\n');
end
disp('-------------------------------------------------------');
end

