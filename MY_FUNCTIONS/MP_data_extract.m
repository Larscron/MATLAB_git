%% matpower data extraction
% this script is used to extract and re-organize matpower data outputs
clear all; clc;
mpc=loadcase(case14); % this is where we declare which case we want to look at
res=runpf(mpc);
% starting with power flow data

%% Measurement data

% we still have not finished the meas matrix, so we will start with that
busnum  =res.bus(:,1);
bustype =res.bus(:,2); % 1=load bus, 2=gen bus, 3=ref bus
P_dis   =res.bus(:,3); % P power distributed / sent out form bus
Q_dis   =res.bus(:,4); % Q power distributed / sent out form bus
V_bus   =res.bus(:,8); % Voltage at bus
%there is more data, but is is not relevant for now

% bus voltage data
V_bus_meas_type = ones(length(busnum),1);   % this is not the bus type!!
V_bus_from = busnum;                        % the from bus is simply the bus
V_bus_to   = zeros(length(busnum),1);       % there is no to bus, so it is set to 0
% combine each of the vectors to a new measurement matrix
Vbus = [V_bus_meas_type, V_bus, V_bus_from, V_bus_to];
% I might want to do something to eliminate the measurements form the bus
% types that do not normally return measuremetns
to_be_removed = find(bustype == 1);         % returns a list of the type 1 buses
removed=0;
for i=1:length(to_be_removed)
    Vbus((to_be_removed(i)-removed),:)=[]; % removes each voltage that would not be measured
    removed=removed+1;
end
% net power injection at each bus (net=gen-dis)
% the power generation data is stored in the .gen field 
gennum  =res.gen(:,1);
P_gen    =res.gen(:,2);
Q_gen    =res.gen(:,3);

%net bus power injection
P_net=zeros(length(busnum),1);
Q_net=zeros(length(busnum),1);
gen_pos =1; gennum_0=[gennum;0];
% go through all buses and subtract distributed power from generated power
for i=1:length(P_dis)
    if gennum_0(gen_pos)==i % checks if there is andy generated power
        P_net(i)=P_gen(gen_pos)-P_dis(i);
        Q_net(i)=Q_gen(gen_pos)-Q_dis(i);
        gen_pos=gen_pos+1;
    else                    % if no generated power, net=-dis
        P_net(i)=-P_dis(i);
        Q_net(i)=-Q_dis(i);
    end
end
% meas types
P_bus_meas_type = ones(length(busnum),1)*2;
Q_bus_meas_type = ones(length(busnum),1)*3;
P_bus_from = busnum;                        % form (just the bus number)
Q_bus_from = busnum;
P_bus_to   = zeros(length(busnum),1);       % ther is not to bus so it is set to 0
Q_bus_to   = zeros(length(busnum),1);
% combine each of the vectors to a new measurement matrix
Pnet = [P_bus_meas_type, P_net, P_bus_from, P_bus_to];
Qnet = [Q_bus_meas_type, Q_net, Q_bus_from, P_bus_to];

% power flow data
fbus=res.branch(:,1);   % from bus
tbus=res.branch(:,2);   % to bus
P_ft=res.branch(:,14);  % the P power flow from the from bus to the to bus           from --> to
Q_ft=res.branch(:,15);  % the Q power flow from the from to the to bus               from --> to
P_tf=res.branch(:,16);  % the P power flow from the to bus to the from bus             to --> from
Q_tf=res.branch(:,17);  % the Q power flow from the to bus to the from bus             to --> from


P_flow_type = ones(length(fbus),1)*4;
Q_flow_type = ones(length(fbus),1)*5;
% combine vectors to form a new measurement matri
% there are two lines to acount include forwards and backwards flow
Pflow=[ P_flow_type, P_ft, fbus, tbus;  
        P_flow_type, P_tf, tbus, fbus];
Qflow=[ Q_flow_type, Q_ft, fbus, tbus;
        Q_flow_type, Q_tf, tbus, fbus];

%combine all measurement matrices
measMtx=[Vbus; Pnet;Qnet;Pflow;Qflow];

% get total number of measurments
numMeas=size(measMtx);
numMeas=numMeas(1);

% add a measurment number to each measurement
measurment = (1:numMeas)';
measMtx=[measurment, measMtx];

% add standard deviation to each measurement
measMtx=add_stdev(measMtx,1,3); % adds a column with the sandard deviation of 1% of each measurment. The measuremnts are located in the 3rd column 

% in the end, we need to arange the meassurement data as:
%         |Msnt |Type | Value | From | To | Rii | 
    
%% line data
% now we need to reformat the matpower data to fit the branch data matrix below
%         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |

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

%% bus data?
% the current structure of my code does not require any bus data exept the
% bus shunt b_bus_sh
busnum=res.bus(:,1); % this has been declared before, but it i doen again for redundancy. 
b_bus_sh=res.bus(:,6);
b_bus_sh_Mtx=[busnum,b_bus_sh];

% fprintf('| Type |  Q_flow  | from |  to  |        | Type |  Q_flow  | from |  to  |\n')
% for bra=1:(2*length(fbus))
%     
%     fprintf('%4i %12.4f %5i %6i           ',Pflow(bra,1),Pflow(bra,2),Pflow(bra,3),Pflow(bra,4))
%     fprintf('%4i %12.4f %5i %6i\n',Qflow(bra,1),Qflow(bra,2),Qflow(bra,3),Qflow(bra,4))
% end



% in the end, we need to arange the meassurement data as:
%         |Msnt |Type | Value | From | To | Rii | 

clearvars -except measMtx lineMtx b_bus_sh

