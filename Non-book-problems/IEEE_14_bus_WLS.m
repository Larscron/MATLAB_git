%% WLS section
clear all; close all; clc;

% there is something strange going on with the equations for h(5) and H52 and potentially h(3) and H42 as well 

% declare the number of buses
busCase=14;

% run MP_data_extract; 

% get our case data form the function files
lineMtx=extr_line_mtx(busCase);
%busMtx=busData(busCase);
measMtx=extr_meas_mtx(busCase);

% Measurement Error matrix
if busCase==6
    stdev=abs(measMtx(:,3))./100;
    stdev_sqr=stdev.^2;
    Ri=diag(stdev_sqr); % declaring the R matrix
else
    stdev=measMtx(:,6);
    stdev_sqr=stdev.^2;
    Ri = diag(stdev_sqr);
end

% from and two 
fb = lineMtx(:,1);                  % From bus number...
tb = lineMtx(:,2);                  % To bus number...
fb_mes = measMtx(:,4);
tb_mes = measMtx(:,5);

nbus = max(max(fb),max(tb));        % no. of buses...
nbra = length(fb);                  % no. of branch

z = measMtx(:,3);                    % measurements data
V = ones(nbus,1);                   % we use a flat start even though we have some (1) V measurements
thet = zeros(nbus,1);               % Initialize the bus angles..
E = [thet(2:end); V];               % State Vector..

y=smally(lineMtx);
g=real(y);
b=imag(y);
bsh=lineShMtx(lineMtx);
b_bus_sh=zeros(nbus,1);             % this on is empty for this example
a=aMtx(lineMtx);                    % the A matrix procedure would be dependent on which side the taps are located

% bus and measurement type
%busType=busMtx(:,2); % not used...
mesType=measMtx(:,2);

vi  = find(mesType == 1); % Index of voltage magnitude measurements..
pin = find(mesType == 2); % Index of real power injection measurements..
qin = find(mesType == 3); % Index of reactive power injection measurements..
pf  = find(mesType == 4); % Index of real powerflow measurements..
qf  = find(mesType == 5); % Index of reactive powerflow measurements..

nvi = length(vi); % Number of Voltage measurements..
npi = length(pin); % Number of Real Power Injection measurements..
nqi = length(qin); % Number of Reactive Power Injection measurements..
npf = length(pf); % Number of Real Power Flow measurements..
nqf = length(qf); % Number of Reactive Power Flow measurements..


it = 1;
tol = 1e-5;
maxerror = 1;


while ((maxerror > tol) && it < 10 )
   %Measurement Function, h
    h1 = V(fb_mes(vi),1);   % this might not be compatible with larger sections
    h2 = zeros(npi,1);
    h3 = zeros(nqi,1);
    h4 = zeros(npf,1);
    h5 = zeros(nqf,1);
    
    for i = 1:npi           % for the real power injection measurements
        k = fb_mes(pin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h2(i) = h2(i) + (a(k,m)*V(k))^2*g(k,m)              -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
        end
    end
    
    for i = 1:nqi           % for the reactive power injection measurements
        k = fb_mes(qin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h3(i) = h3(i) -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
        end
        h3(i) = h3(i) - V(k)^2*b_bus_sh(k);
    end
    
    for i = 1:npf           % for the real power flow measurements
        k = fb_mes(pf(i));
        m = tb_mes(pf(i));
        thet_km=thet(k)-thet(m);
        h4(i) = (a(k,m)*V(k))^2*g(k,m)              -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
    end
    
    for i = 1:nqf           % for the reactive power flow measurements
        k = fb_mes(qf(i));
        m = tb_mes(qf(i));
        thet_km=thet(k)-thet(m);
        h5(i) = -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
    end
    
    h = [h1; h2; h3; h4; h5];
   
    % Residue..
    r = z - h;
    
    % Jacobian..
    % H11 - Derivative of V with respect to angles.. All Zeros
    H11 = zeros(nvi,nbus-1);

    % H12 - Derivative of V with respect to V.. allways 1 for Vk/Vk
    H12 = zeros(nvi,nbus);
    for k = 1:nvi
        for m = 1:nbus
            if m == k
                H12(k,m) = 1;                                               %ERROR FOUND HERE
            end
        end
    end
    
    % H21 - Derivative of Real Power Injections with Angles..
    H21 = zeros(npi,nbus-1);
    for i = 1:npi % step through all power injection measurements
        k = fb_mes(pin(i)); % which bus is this power injection measurement at
        for var = 1:(nbus-1) % step through all variable locations of the H (thet(2) to thet thet(14))
            if var+1 == k % check if power injection bus match thet (Pk/thet_k)
                for m = 1:nbus % step through all busses connected to the m bus, this alos include bus #1
                    thet_km=thet(k)-thet(m);
                    H21(i,var) = H21(i,var) +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km) -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);
                end
            else % Pk/thet_m
                m=var+1;
                thet_km=thet(k)-thet(m);
                H21(i,var) = -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);
            end
        end
    end
    
    % H22 - Derivative of Real Power Injections with V..
    H22 = zeros(npi,nbus);
    for i = 1:npi
        k = fb_mes(pin(i)); % which bus is this power injection measurement at
        for var = 1:(nbus) % step through all variable locations of the H (V(1) to thet V(14))
            if var == k % check if power injection bus match the variable location (Pk/Vk)
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    H22(i,var) = H22(i,var)     +(2*a(k,m)^2*V(k)*g(k,m)-(a(k,m)*a(m,k)*V(m))*g(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(m))*b(k,m)*sin(thet_km));
                end
            else % Pk/Vm
                m=var;
                thet_km=thet(k)-thet(m);
                H22(i,var) =                                              -(a(k,m)*a(m,k)*V(k))*g(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(k))*b(k,m)*sin(thet_km);
            end
        end
    end
    
    % H31 - Derivative of Reactive Power Injections with Angles..
    H31 = zeros(nqi,nbus-1);
    for i = 1:nqi                   % step through all Q injections measurements
        k = fb_mes(qin(i));         % which bus is this Q injection measurement at
        for var = 1:(nbus-1)        % step through all variable locations of the H (thet(2) to thet thet(14))
            if var+1 == k           % check if power injection bus match thet
                for m = 1:nbus      % step through all busses connected to the m bus, this alos include bus 1 Qk/thet_k
                    thet_km=thet(k)-thet(m);
                    H31(i,var) = H31(i,var)-(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)-(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km);
                end
            else % Qk/thet_m
                m=var+1;
                thet_km=thet(k)-thet(m);
                H31(i,var) = +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km) +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km);
            end
        end
    end
    
    % H32 - Derivative of Reactive Power Injections with V..
    H32 = zeros(nqi,nbus);
    for i = 1:nqi
        k = fb_mes(qin(i));         % which bus is this Q injection measurement at
        for var = 1:(nbus)          % step through all variable locations of the H (V(1) to thet V(14))
            if var == k             % check if Q injection bus match the variable location (Qk/Vk)
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    H32(i,var) = H32(i,var) - 2*a(k,m)^2*V(k)*(b(k,m)+bsh(k,m))+(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km);
                end
                H32(i,var) = H32(i,var) -2*V(k)*b_bus_sh(k);
            else % Qk/Vm
                m=var;
                thet_km=thet(k)-thet(m);
                H32(i,var) = +(a(k,m)*a(m,k)*V(k))*b(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(k))*g(k,m)*sin(thet_km);
            end
        end
    end
    
    % H41 - Derivative of Real Power Flows with Angles..
    H41 = zeros(npf,nbus-1);
    for i = 1:npf
        k = fb_mes(pf(i));          % form bus
        m = tb_mes(pf(i));          % to bus
        for var = 1:(nbus-1)        % we step through all the thet variables
            if var+1 == k           % P_km/thet_k
                thet_km=thet(k)-thet(m);
                H41(i,var) = +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km) -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);
            else if var+1 == m      % P_km/thet_m
                thet_km=thet(k)-thet(m);
                H41(i,var) = -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km) +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);
                else
                    H41(i,var) = 0;
                end
            end
        end
    end
    
    % H42 - Derivative of Real Power Flows with V..
    H42 = zeros(npf,nbus);
    for i = 1:npf
        k = fb_mes(pf(i));            % form bus
        m = tb_mes(pf(i));            % to bus
        for var = 1:nbus            % we step through all the V variables
            if var == k             % P_km/V_k
                thet_km=thet(k)-thet(m);
                H42(i,var) = 2*a(k,m)^2*V(k)*g(k,m)-(a(k,m)*a(m,k)*V(m))*g(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(m))*b(k,m)*sin(thet_km);
            else if var == m        % P_km/V_m
                thet_km=thet(k)-thet(m);
                H42(i,var) = -(a(k,m)*a(m,k)*V(k))*g(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(k))*b(k,m)*sin(thet_km);
                else
                    H42(i,var) = 0;
                end
            end
        end
    end
    
    % H51 - Derivative of Reactive Power Flows with Angles..
    H51 = zeros(nqf,nbus-1);
    for i = 1:nqf
        k = fb_mes(qf(i));            % form bus
        m = tb_mes(qf(i));            % to bus
        for var = 1:(nbus-1)
            if var+1 == k           % Q_km/thet_k
                thet_km=thet(k)-thet(m);
                H51(i,var) = -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)-(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km);
            else if var+1 == m      % Q_km/thet_m
                thet_km=thet(k)-thet(m);
                H51(i,var) = +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km) +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km);
                else
                    H51(i,var) = 0;
                end
            end
        end
    end
    
    % H52 - Derivative of Reactive Power Flows with V..
    H52 = zeros(nqf,nbus);                                                  % THERE IS NO THET DECLARATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:nqf
        k = fb_mes(qf(i));
        m = tb_mes(qf(i));
        for var = 1:nbus
            if var == k             % Q_km/V_k
                % i do not understand why the commented out equation is not
                % correct, but it is the active equation tha makes the
                % system converge to the correct values
                H52(i,var) = -2*a(k,m)^2*V(k)*(b(k,m)+bsh(k,m)) +(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km);
                %H52(i,var) = -2*a(k,m)*V(k)*(b(k,m)+bsh(k,m)) +(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km);
            else if var == m        % Q_km/V_m
                H52(i,var) =                                    +(a(k,m)*a(m,k)*V(k))*b(k,m)*cos(thet_km)-(a(k,m)*a(m,k)*V(k))*g(k,m)*sin(thet_km);
                else
                    H52(i,var) = 0;
                end
            end
        end
    end
    
    
    % Measurement Jacobian, H..
    H =[H11 H12; 
        H21 H22;
        H31 H32;
        H41 H42;
        H51 H52];
    
    % Gain Matrix, Gm..
    Gm = H'*inv(Ri)*H;
    
    %Objective Function..
    J = sum(inv(Ri)*r.^2);
    
    % State Vector..
    dE = inv(Gm)*(H'*inv(Ri)*r);
    E = E + dE;
    thet(2:end) = E(1:nbus-1);
    V = E(nbus:end);
    it = it + 1;
    maxerror = max(abs(dE))
end

fprintf('convergent occurred!!\nWE ARE DONE!!!\n')

printStateVar(V,thet)


        
