function [] = printStateVar(Vmagnitude,Vangle)
%PRINTSTATEVAR Prints the voltage magnitude and angle (in degrees) of the
%state variables

Vangle_deg = Vangle*180/pi;
nbus=length(Vmagnitude);

disp('--------------------------');
disp('| Bus |    V   |  Angle  | ');
disp('| No  |   pu   |  Degree | ');
disp('--------------------------');
for m = 1:nbus
    fprintf('%4g', m); fprintf('  %8.4f', Vmagnitude(m)); fprintf('   %8.4f', Vangle_deg(m)); fprintf('\n');
end
disp('--------------------------');
end

