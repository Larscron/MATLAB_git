function [lineMtx] = extr_line_mtx(busNum)
%EXTR_LINE_MTX Summary of this function goes here
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

fbus=res.branch(:,1);   % from bus
tbus=res.branch(:,2);   % to bus
R   =res.branch(:,3);
X   =res.branch(:,4);
B   =res.branch(:,5); % the matpower line shunt is given in B and not B/2!!
tap =res.branch(:,9); % taps that are not used are set to 0!!!!
% rat A, B, and C, are only relevant when we look at individual lines
% i think the angle column, colu,n 10, refere to pahse shifting transformers

% B an tap need som more processing:
B_2 = B.*0.5;

A=ones(length(tap),1);
for i=1:length(tap)
    if tap(i)~=0
        A(i)=tap(i);
    end
end

%combining branch vectors to make the branch data matr
%        |  From |  To    | R    |  X   |   B/2  |  X'mer   |
%        |  Bus  |  Bus   | pu   |  pu  |   pu   |  TAP (a) |
lineMtx=[   fbus,   tbus,   R,      X,      B_2,    A];

end

