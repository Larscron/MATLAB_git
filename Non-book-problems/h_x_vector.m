%% h() calculations
clear all; close all; clc

nbus=5;
%nbus=14;

measMtx=extr_meas_mtx(nbus);
lineMtx=extr_line_mtx(nbus);


fb = lineMtx(:,1);                  % From bus number...
tb = lineMtx(:,2);                  % To bus number...
fb_mes = measMtx(:,4);              % Form bus, measurement
tb_mes = measMtx(:,5);              % To bus, measurement

z = measMtx(:,3);                    % measurements data
V = ones(nbus,1);                   % we use a flat start even though we have some (1) V measurements
thet = zeros(nbus,1);               % Initialize the bus angles.
E = [thet(2:end); V];               % State Vector..

y=smally(lineMtx);
g=real(y);
b=imag(y);
bsh=lineShMtx(lineMtx);
b_bus_sh=extr_bus_sh_vec(nbus);             % this on is empty for this example
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

%updating the stat variables with the finall MATPower values
stateMtx=extr_stat_mtx(nbus);
V=stateMtx(:,1);
thet=stateMtx(:,2);
%Measurement Function, h
    h1 = V(fb_mes(vi),1);   % Voltage values
    h2 = zeros(npi,1);      % P injection
    h3 = zeros(nqi,1);      % Q injection
    h4 = zeros(npf,1);      % P flow
    h5 = zeros(nqf,1);      % Q flow
    
    P_inj=zeros(nbus,3);
    Q_inj=zeros(nbus,3);
    P_flo=zeros(nbus,3);
    Q_flo=zeros(nbus,3);
    
    for i = 1:npi           % for the real power injection measurements
        k = fb_mes(pin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h2(i) = h2(i) +     (a(k,m)*V(k))^2*g(k,m)      -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)      -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
        end
        P_inj(i,1)=h2(i);   % P injection
        P_inj(i,2)=i;       % bus
    end
    
    for i = 1:nqi           % for the reactive power injection measurements
        k = fb_mes(qin(i));
        for m = 1:nbus
            thet_km=thet(k)-thet(m);
            h3(i) = h3(i) -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
        end
        h3(i) = h3(i) - V(k)^2*b_bus_sh(k);
        
        Q_inj(i,1)=h2(i);   % Q injection
        Q_inj(i,2)=i;       % bus
    end
    
    for i = 1:npf           % for the real power flow measurements
        k = fb_mes(pf(i));
        m = tb_mes(pf(i));
        thet_km=thet(k)-thet(m);
        h4(i) =                 (a(k,m)*V(k))^2*g(k,m)      -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)      -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
        
        P_flo(i,1)=h4(i);   % P injection
        P_flo(i,2)=k;       % from bus
        P_flo(i,3)=m;       % to bus
    end
    
    for i = 1:nqf           % for the reactive power flow measurements
        k = fb_mes(qf(i));
        m = tb_mes(qf(i));
        thet_km=thet(k)-thet(m);
        
        h5(i) = -(a(k,m)*V(k))^2*(b(k,m)+bsh(k,m))   +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);

        Q_flo(i,1)=h5(i);   % P injection
        Q_flo(i,2)=k;       % from bus
        Q_flo(i,3)=m;       % to bus
    end
    
% multiplying the power calculations with the base (100 MVA)
h2=h2.*100;
h3=h3.*100;
h4=h4.*100;
h5=h5.*100;
    
h = [h1; h2; h3; h4; h5];

comp=h-measMtx(:,3);
error=0;
for i=1:length(comp)
    if abs(comp(i)) >= 0.00001
        fprintf('We have an error at equation %i\n',i)
        
        switch measMtx(i,2)
            case 1
                disp('This is a V problem')
            case 2
                fprintf('This is a P injection problem \nCheck h2 for %i\n', measMtx(i,4)')
                
            case 3
                fprintf('This is a Q injection problem \nCheck h3 for bus %i\n', measMtx(i,4))
            case 4
                disp('This is a P flow problem')
            case 5
                disp('This is a Q flow problem')
        end
        error=error+1;
    end
end

if error==0
    disp('NO ERRORS!!!')
else
    fprintf('the error count was: %i\n',error)
end