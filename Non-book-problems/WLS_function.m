function [ WLSres ] = WLS_function(nbus, measMtx, lineMtx, b_bus_sh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%pctSTdev=1;

fb = lineMtx(:,1);                      % From bus number...
tb = lineMtx(:,2);                      % To bus number...
fb_mes = measMtx(:,4);                  % Form bus, measurement
tb_mes = measMtx(:,5);                  % To bus, measurement

z = measMtx(:,3);                       % measurements data
stdev=measMtx(:,6);
stdev2=stdev.^2;                        % squaring the standard deviation to get Rii
R=diag(stdev2);                         % adding the square of the standard deviation of the measurements

V = ones(nbus,1);                       % we use a flat start even though we have some (1) V measurements
thet = zeros(nbus,1);                   % Initialize the bus angles.
x = [thet(2:end); V];                   % State Vector..

y=smally(lineMtx);
g=real(y);
b=imag(y);
bsh=lineShMtx(lineMtx);
%b_bus_sh=extr_bus_sh_vec(nbus);             % this on is empty for this example
a=aMtx(lineMtx);                    % the A matrix procedure would be dependent on which side the taps are located

% bus and measurement type
%busType=busMtx(:,2); % not used...
nmeas=length(measMtx(:,2));
mesType=measMtx(:,2);
vi  = find(mesType == 1); % Index of voltage magnitude measurements..
vi_bus = fb_mes(vi);
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
tol = 1e-8;
maxerror = 1;


while ((maxerror > tol) && it < 10 )

%Measurement Function, h
    h1 = V(fb_mes(vi),1);   % Voltage values
    h2 = zeros(npi,1);      % P injection
    h3 = zeros(nqi,1);      % Q injection
    h4 = zeros(npf,1);      % P flow
    h5 = zeros(nqf,1);      % Q flow
     
    for i = 1:npi           % for the real power injection measurements
        k = fb_mes(pin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h2(i) = h2(i) +     (a(k,m)*V(k))^2*g(k,m)      -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)      -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
        end
    end     %h2
    
    for i = 1:nqi           % for the reactive power injection measurements
        k = fb_mes(qin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h3(i) = h3(i) -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
        end
        h3(i) = h3(i) - V(k)^2*b_bus_sh(k);
    end %h3
    
    for i = 1:npf           % for the real power flow measurements
        k = fb_mes(pf(i));
        m = tb_mes(pf(i));
        thet_km=thet(k)-thet(m);
        h4(i) =                 (a(k,m)*V(k))^2*g(k,m)      -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)      -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
    end          %h4
    
    for i = 1:nqf           % for the reactive power flow measurements
        k = fb_mes(qf(i));
        m = tb_mes(qf(i));
        thet_km=thet(k)-thet(m);
        h5(i) = -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
    end      %h5
    
    %combining all the h functions
    h = [h1; h2; h3; h4; h5];

    % Jacobian calculations
    % H11 - Derivative of V with respect to angles.. All Zeros              Vk/?k and Vk/?m
    H11 = zeros(nvi,nbus-1);
    
    % H12 - Derivative of V with respect to V.. allways 1 for               Vk/Vk and Vk/Vm 
    H12 = zeros(nvi,nbus);
    for k = 1:nvi
        H12(k,fb_mes(vi(k))) = 1; % consider expanding for more clarity
    end
    
    % H21 - Derivative of Real Power Injections with Angles..               Pk/?k and Pk/?m
    H21 = zeros(npi,nbus-1);
    for i = 1:npi % step through all power injection measurements
        k = fb_mes(pin(i)); % which bus is this power injection measurement at
        for var = 1:(nbus-1) % step through all variable locations of the H (?_2 to thet ?_14)
            if var+1 == k % check if power injection bus match thet (Pk/?_k)
                for m = 1:nbus % step through all busses connected to the m bus, this alos include bus #1
                    thet_km=thet(k)-thet(m);
                    H21(i,var) = H21(i,var) +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km) -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);% Pk/?_k
                end
            else % Pk/?_m
                m=var+1;
                thet_km=thet(k)-thet(m);          
                H21(i,var) =                -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km);% Pk/?_m
            end
        end
    end
    
    % H22 - Derivative of Real Power Injections with V..                    Pk/Vk and Pk/Vm
    H22 = zeros(npi,nbus);
    for i = 1:npi
        k = fb_mes(pin(i)); % which bus is this power injection measurement at
        for var = 1:(nbus) % step through all variable locations of the H (V(1) to thet V(14))
            if var == k % check if power injection bus match the variable location (Pk/Vk)
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);                 
                    H22(i,var) = H22(i,var)     +(2*a(k,m)^2*V(k)*g(k,m)-(a(k,m)*a(m,k)*V(m))*g(k,m)*cos(thet_km)       -(a(k,m)*a(m,k)*V(m))*b(k,m)*sin(thet_km));% Pk/Vk
                end
            else % Pk/Vm
                m=var;
                thet_km=thet(k)-thet(m);
                H22(i,var) =                                            -(a(k,m)*a(m,k)*V(k))*g(k,m)*cos(thet_km)       -(a(k,m)*a(m,k)*V(k))*b(k,m)*sin(thet_km); % Pk/Vm
            end
        end
    end
    
    % H31 - Derivative of Reactive Power Injections with Angles..           Qk/?k and Qk/?m
    H31 = zeros(nqi,nbus-1);
    for i = 1:nqi                   % step through all Q injections measurements
        k = fb_mes(qin(i));         % which bus is this Q injection measurement at
        for var = 1:(nbus-1)        % step through all variable locations of the H (thet(2) to thet thet(14))
            if var+1 == k           % check if power injection bus match thet
                for m = 1:nbus      % step through all busses connected to the m bus, this alos include bus 1 Qk/?_k
                    thet_km=thet(k)-thet(m);                  
                    H31(i,var) = H31(i,var)             -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km); %Qk/?_k
                end            
            else % Qk/thet_m
                m=var+1;
                thet_km=thet(k)-thet(m);              
                H31(i,var) =                            +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km); %Qk/?_m       
            end
        end
    end
    
    % H32 - Derivative of Reactive Power Injections with V..                Qk/Vk and Qk/Vm
    H32 = zeros(nqi,nbus);
    for i = 1:nqi
        k = fb_mes(qin(i));         % which bus is this Q injection measurement at
        for var = 1:(nbus)          % step through all variable locations of the H (V(1) to thet V(14))
            if var == k             % check if Q injection bus match the variable location (Qk/Vk)
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);                    
                    H32(i,var) = H32(i,var) - 2*a(k,m)^2*V(k)*(b(k,m)+bsh(k,m))+(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)        -(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km); % Qk/Vk
                end          
                H32(i,var) = H32(i,var) -2*V(k)*b_bus_sh(k);
            else % Qk/Vm
                m=var;
                thet_km=thet(k)-thet(m);               
                H32(i,var) =                            +(a(k,m)*a(m,k)*V(k))*b(k,m)*cos(thet_km)       -(a(k,m)*a(m,k)*V(k))*g(k,m)*sin(thet_km); % Qk/Vm)                
            end
        end
    end
    
    % H41 - Derivative of Real Power Flows with Angles..                    Pkm/?k and Pkm/?m
    H41 = zeros(npf,nbus-1);
    for i = 1:npf
        k = fb_mes(pf(i));          % form bus
        m = tb_mes(pf(i));          % to bus
        thet_km=thet(k)-thet(m);
        for var = 1:(nbus-1)        % we step through all the thet variables
            if var+1 == k           % P_km/thet_k               
                H41(i,var) =            +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)      -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km); % Pkm/?k
            else
                if var+1 == m      % P_km/thet_m                  
                    H41(i,var) =        -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)      +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km); % Pkm/?k
                end
            end
        end
    end
    
    % H42 - Derivative of Real Power Flows with V..                         Pkm/Vk and Pkm/Vm
    H42 = zeros(npf,nbus);
    for i = 1:npf
        k = fb_mes(pf(i));            % form bus
        m = tb_mes(pf(i));            % to bus
        thet_km=thet(k)-thet(m);
        for var = 1:nbus              % we step through all the V variables
            if var == k               % P_km/V_k               
                H42(i,var) = 2*a(k,m)^2*V(k)*g(k,m)     -(a(k,m)*a(m,k)*V(m))*g(k,m)*cos(thet_km)           -(a(k,m)*a(m,k)*V(m))*b(k,m)*sin(thet_km); % Pkm/Vk
            else
                if var == m        % P_km/V_m                 
                    H42(i,var) =    -(a(k,m)*a(m,k)*V(k))*g(k,m)*cos(thet_km)           -(a(k,m)*a(m,k)*V(k))*b(k,m)*sin(thet_km); % Pkm/Vm
                end
            end
        end
    end
    
    % H51 - Derivative of Reactive Power Flows with Angles..                Qkm/?k and Qkm/?m
    H51 = zeros(nqf,nbus-1);
    for i = 1:nqf
        k = fb_mes(qf(i));            % form bus
        m = tb_mes(qf(i));            % to bus
        thet_km=thet(k)-thet(m);
        for var = 1:(nbus-1)
            if var+1 == k           % Q_km/thet_k               
                H51(i,var) =                    -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km); % Qkm/?k
            else
                if var+1 == m      % Q_km/thet_m
                    H51(i,var) =                +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km); % Qkm/?m
                end
            end
        end
    end
    
    % H52 - Derivative of Reactive Power Flows with V..                     Qkm/Vk and Qkm/Vm
    H52 = zeros(nqf,nbus);
    for i = 1:nqf
        k = fb_mes(qf(i));
        m = tb_mes(qf(i));
        thet_km=thet(k)-thet(m);
        for var = 1:nbus
            if var == k             % Q_km/V_k               
                H52(i,var) = -2*a(k,m)^2*V(k)*(b(k,m)+bsh(k,m))     +(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)       -(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km); % Qkm/Vk

            else
                if var == m        % Q_km/V_m                   
                    H52(i,var) =                                    +(a(k,m)*a(m,k)*V(k))*b(k,m)*cos(thet_km)       -(a(k,m)*a(m,k)*V(k))*g(k,m)*sin(thet_km); % Qkm/Vm
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
    Gm = H'*inv(R)*H;
    
    % Residue..
    del_z = z - h;
    
    % Objective Function..
    J = sum(inv(R)*del_z.^2);
    
    % The projection matrix
    K=H*inv(Gm)*transpose(H)*inv(R);
    % Residual vector?
    
    % State Vector..
    del_x_hat = inv(Gm)*(H'*inv(R)*del_z);
    x = x + del_x_hat;
    thet(2:end) = x(1:nbus-1);
    V = x(nbus:end);
    it = it + 1;
    maxerror = max(abs(del_x_hat));
end
it=it-1;

% save data to output structure
WLSres=struct();
WLSres.H=H;
WLSres.G=Gm;
WLSres.R=R;
WLSres.J=J;

WLSres.del_z=del_z; WLSres.del_x_hat=del_x_hat;
WLSres.Nmeas=nmeas; % could probably be moved outside of function
%WLSres.STdev=stdev;
%WLSres.GE_loc=GE_loc; %just so that I will know for testing

if 1
    fprintf('%i Iterations were completed\n',it)
    if maxerror < tol
        %printStateVar(V,thet)
        print_comp(V,thet)
        fprintf('THERE WAS CONVERGENCE!!!\n')
        fprintf('%i Iterations were completed\n',it)
        fprintf('The velue of J was:\n%f\n',J)
        %fprintf('The maximum remaining error is:\n%.14f\n',max(abs(del_x_hat)))

    else
        fprintf('there was no convergence :(\n')
    end
end

end

