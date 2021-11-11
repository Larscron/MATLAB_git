%% matpower data extraction
% this script is used to extract and re-organize matpower data outputs
clear all; clc;
mpc=loadcase(case14); % this is where we declare which case we want to look at
res=runpf(mpc);
% starting with power flow data


%% branch data
res.branch
fbus=res.branch(:,1);
tbus=res.branch(:,2);
r   =res.branch(:,3);
x   =res.branch(:,4);
bsh =res.branch(:,5);
% rat A, B, and C, are only relevant when we look at individual lines
tap =res.branch(:,9); % taps that are not used are set to 0!!!!
% i think the angle column, colu,n 10, refere to pahse shifting transformers
P_ft=res.branch(:,14); % the P power flow from the from bus to the to bus           from --> to
Q_ft=res.branch(:,15); % the Q power flow from the from to the to bus               from --> to
P_tf=res.branch(:,16); % the P power flow from the to bus to the from bus             to --> from
Q_tf=res.branch(:,17); % the Q power flow from the to bus to the from bus             to --> from



P_flow_type = ones(length(fbus),1)*4;
Q_flow_type = ones(length(fbus),1)*5;

Pflow=[ P_flow_type, P_ft, fbus, tbus;
        P_flow_type, P_tf, tbus, fbus];
Qflow=[ Q_flow_type, Q_ft, fbus, tbus;
        Q_flow_type, Q_tf, tbus, fbus];
    
% fprintf('| Type |  Q_flow  | from |  to  |        | Type |  Q_flow  | from |  to  |\n')
% for bra=1:(2*length(fbus))
%     
%     fprintf('%4i %12.4f %5i %6i           ',Pflow(bra,1),Pflow(bra,2),Pflow(bra,3),Pflow(bra,4))
%     fprintf('%4i %12.4f %5i %6i\n',Qflow(bra,1),Qflow(bra,2),Qflow(bra,3),Qflow(bra,4))
% end



% in the end, we need to arange the meassurement data as:
%         |Msnt |Type | Value | From | To | Rii | 

%% Bus and gen data

% we still have not finished the meas matrix, so we will start with that
busnum  =res.bus(:,1);
bustype =res.bus(:,2); % 1=load bus, 2=gen bus, 3=ref bus
Pdis    =res.bus(:,3);
Qdis    =res.bus(:,4);
V_bus    =res.bus(:,8);
%there is more data, but is is not relevant for now

gennum  =res.gen(:,1);
Pgen    =res.gen(:,2);
Qgen    =res.gen(:,3);

%net bus power injection
P_net=zeros(length(busnum),1);
Q_net=zeros(length(busnum),1);

gen_pos =1; gennum_0=[gennum;0];
for i=1:length(Pdis)
    if gennum_0(gen_pos)==i
        P_net(i)=Pgen(gen_pos)-Pdis(i);
        Q_net(i)=Qgen(gen_pos)-Qdis(i);
        gen_pos=gen_pos+1;
    else
        P_net(i)=-Pdis(i);
        Q_net(i)=-Qdis(i);
    end
end

% arrange in meas format
% meas types
P_bus_meas_type = ones(length(busnum),1)*2;
Q_bus_meas_type = ones(length(busnum),1)*3;
%form (just the bus number)
P_bus_from = busnum;
Q_bus_from = busnum;
P_bus_to   = zeros(length(busnum),1);
Q_bus_to   = zeros(length(busnum),1);

Pnet = [P_bus_meas_type, P_net, P_bus_from, P_bus_to];
Qnet = [Q_bus_meas_type, Q_net, Q_bus_from, P_bus_to];

% fprintf('| Type |  P_net  | from |  to  |        | Type |  Q_net  | from |  to  |\n')
% for bra=1:length(busnum)
%     
%     fprintf('%4i %12.4f %5i %6i           ',Pnet(bra,1),Pnet(bra,2),Pnet(bra,3),Pnet(bra,4))
%     fprintf('%4i %12.4f %5i %6i\n',Qnet(bra,1),Qnet(bra,2),Qnet(bra,3),Qnet(bra,4))
% end

% voltage meassurements
V_bus_meas_type = ones(length(busnum),1);
V_bus_from = busnum;
V_bus_to   = zeros(length(busnum),1);

Vbus = [V_bus_meas_type, V_bus, V_bus_from, V_bus_to];

% fprintf('| Type |  V_bus  | from |  to  |\n')
% for bus=1:length(busnum)
%     fprintf('%4i %12.4f %5i %6i\n',Vbus(bus,1),Vbus(bus,2),Vbus(bus,3),Vbus(bus,4))
% end

%combining all meassurements
measMtx=[Vbus; Pnet;Qnet;Pflow;Qflow];

nummes=size(measMtx);
nummes=nummes(1);

% fprintf('| Type |  V_bus  | from |  to  |\n')
% for bus=1:nummes
%     fprintf('%4i %12.4f %5i %6i\n',measMtx(bus,1),measMtx(bus,2),measMtx(bus,3),measMtx(bus,4))
% end
% add measurment number number
measNum=[1:nummes]';
measMtx=[measNum, measMtx];

% add standard deviation to each measurement
measMtx=add_stdev(measMtx,1,3)

% in the end, we need to arange the meassurement data as:
%         |Msnt |Type | Value | From | To | Rii | 

fprintf('just ran the old one\n')
