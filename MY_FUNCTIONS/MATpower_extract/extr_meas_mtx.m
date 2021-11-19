function [measMtx] = extr_meas_mtx(busNum)
%EXTR_MEAS_MTX Summary of this function goes here
%   Detailed explanation goes here
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
gennum   =res.gen(:,1); % this list when gennerator we are talking about
P_gen    =res.gen(:,2);
Q_gen    =res.gen(:,3);

% making a P and Q gen matrix of the same dimentions as the P and Q dis matices:
P_gen_temp=zeros(length(busnum),1);
Q_gen_temp=zeros(length(busnum),1);
for i=1:length(gennum)
    P_gen_temp(gennum(i))=P_gen_temp(gennum(i))+P_gen(i);
    Q_gen_temp(gennum(i))=Q_gen_temp(gennum(i))+Q_gen(i);
end
P_gen=P_gen_temp;
Q_gen=Q_gen_temp;

%net bus power injection
P_net = P_gen-P_dis;
Q_net = Q_gen-Q_dis;

% meas types
P_bus_meas_type = ones(length(busnum),1)*2;
Q_bus_meas_type = ones(length(busnum),1)*3;
P_bus_from = busnum;                        % form (just the bus number)
Q_bus_from = busnum;
P_bus_to   = zeros(length(busnum),1);       % ther is not to bus so it is set to 0
Q_bus_to   = zeros(length(busnum),1);
% combine each of the vectors to a new measurement matrix
Pnet = [P_bus_meas_type, P_net, P_bus_from, P_bus_to];
Qnet = [Q_bus_meas_type, Q_net, Q_bus_from, Q_bus_to];

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


% scale down power measurements
P_meas=find(measMtx(:,2)~=1);
measMtx(P_meas,3)=measMtx(P_meas,3)./100; % this might have ruined the compatibility with test_h_x_vector

fprintf('test this')

end

