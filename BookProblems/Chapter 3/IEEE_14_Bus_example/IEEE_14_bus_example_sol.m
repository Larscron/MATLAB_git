% Power System Simulation Lab
% Ashish Kumar 16MT000947

%% data
clear;clc
% IEEE 14 bus system
%         Bus No  Vm  Theta  Pg    Qg     P1     Q1     Qmax    Qmin  Type
busdata= [ 1    1.060   0     0    0     0       0       0       0     1;
           2    1.045   0   40.0 45.41  21.7  12.7      -40     50     2;
           3    1.010   0     0  25.28  94.2  19.0       0      40     2;
           4    1.000   0     0    0    47.8   -3.9      0       0     3;
           5    1.000   0     0    0     7.6     1.6     0       0     3;
           6    1.070   0     0  13.62  11.2   7.5    -6.0     24.0    3;
           7    1.000   0     0    0      0      0      0       0      3;
           8    1.090   0     0  18.24    0      0     -6.0    24.0    3;
           9    1.000   0     0    0      29.5   16.6    0       0     3;
           10   1.000   0     0    0      9.0     5.8    0       0     3;
           11   1.000   0     0    0      3.5     1.8    0       0     3;
           12   1.000   0     0    0      6.1     1.6    0       0     3;
           13   1.000   0     0    0      13.5    5.8    0       0     3;
           14   1.000   0     0    0      14.9    5.0    0       0     3];
       
% Susceptance
Shuntdata=[0 0 0 0 0 0 0 0 sqrt(-1)*0.190 0 0 0 0 0]' ;
% Line Data for IEEE 14 bus system
%         Line No   From to    R          X        B/2      T    
linedata= [ 1        1   2    0.01938   0.05917  0.02640    1;
            2        2   3    0.04699   0.19797  0.02190    1;
            3        2   4    0.05811   0.17632  0.01870    1;
            4        1   5    0.05403   0.22304  0.02460    1;
            5        2   5    0.05695   0.17388  0.01700    1;
            6        3   4    0.06701   0.17103  0.01730    1;
            7        4   5    0.01335   0.04211  0.00640    1;
            8        5   6    0.00000   0.25202  0.00000   0.932;
            9        4   7    0.00000   0.20912  0.00000   0.978;
            10       7   8    0.00000   0.17615  0.00000    1;
            11       4   9    0.00000   0.55618  0.00000    0.969;
            12       7   9    0.00000   0.11001  0.00000    1;
            13       9   10   0.03181   0.08450  0.00000    1;
            14       6   11   0.09498   0.19890  0.00000    1;
            15       6   12   0.12291   0.25581  0.00000    1;
            16       6   13   0.06615   0.13027  0.00000    1;
            17       9   14   0.12711   0.27038  0.00000    1;
            18      10   11   0.08205   0.19207  0.00000    1;
            19      12   13   0.22092   0.19988  0.00000    1;
            20      13   14   0.17093   0.34802  0.00000    1];
        
% Ybus program
fb = linedata(:,2);
tb = linedata(:,3);
R =  linedata(:,4);
X =  linedata(:,5);
B = linedata(:,6);
T =  linedata(:,7);
z =  R+sqrt(-1)*X;
y = 1./z;
b= sqrt(-1)*B;
nbus = max(max(fb),max(tb)); % no of bus
nbranch = length(fb); % no of branch
Y = zeros(nbus,nbus); % initialise Ybus
%formation of offdiagonal element
for k=1: nbranch
    Y(fb(k),tb(k))= Y(fb(k),tb(k))-y(k)/T(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
%formation of diagonal element
for m=1 : nbus
    for n = 1:nbranch
        if fb(n)== m
            Y(m,m)= Y(m,m)+(y(n)/(T(k)^2))+b(n);
        elseif tb(n)== m
            Y(m,m)= Y(m,m)+y(n)+b(n);
        end
    end 
end
%entry of shuntdata in IEEE 14 Bus system
 for i = 1:14
      Y(i,i)= Y(i,i)+ Shuntdata(i);  
 end
 
 % Newton-Raphson Load Flow Analysis..
%clear;clc
%Data14;
Y;          % Calling ybusppg.m to get Y-Bus Matrix..
busd = busdata;      % Calling busdatas..
BMva = 100;                 % Base MVA..
bus = busdata(:,1);            % Bus Number..
type = busd(:,10);           % Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,2);              % Specified Voltage..
del = busd(:,3);            % Voltage Angle..
Pg = busd(:,4)/BMva;        % PGi..
Qg = busd(:,5)/BMva;        % QGi..
Pl = busd(:,6)/BMva;        % PLi..
Ql = busd(:,7)/BMva;        % QLi..
Qmin = busd(:,8)/BMva;      % Minimum Reactive Power Limit..
Qmax = busd(:,9)/BMva;     % Maximum Reactive Power Limit..
P = Pg - Pl;                % Pi = PGi - PLi..
Q = Qg - Ql;                % Qi = QGi - QLi..
Psp = P;                    % P Specified..
Qsp = Q;                    % Q Specified..
G = real(Y);                % Conductance matrix..
B = imag(Y);                % Susceptance matrix..
pv = find(type == 2 | type == 1);   % PV Buses..
pq = find(type == 3);               % PQ Buses..
npv = length(pv);                   % No. of PV buses..
npq = length(pq);                   % No. of PQ buses..
Tol = 1;  
Iter = 1;

%%
while (Tol > 0.00001)   % Iteration starting..
    
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculate P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
    % Checking Q-limit violations..
    if Iter <= 7 && Iter > 2    % Only checked up to 7th iterations..
        for n = 2:nbus
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];     % Jacobian Matrix..
    
    X = inv(J)*M;           % Correction Vector
    dTh = X(1:nbus-1);      % Change in Voltage Angle..
    dV = X(nbus:end);       % Change in Voltage Magnitude..
    
    % Updating State Vectors..
    del(2:nbus) = dTh + del(2:nbus);    % Voltage Angle..
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end
    
    Iter = Iter + 1;
    Tol = max(abs(M));                  % Tolerance..
    
end
P
Q
G
B
V
del
%Loss Calculation
Vs=V(fb);
Vr=V(tb);
dels=del(fb);
delr=del(tb);
G=real(y);
for i=1:20
    
Pl(i)=G(i)*((Vs(i)^2)+(Vr(i)^2)-2*Vs(i)*Vr(i)*cos(dels(i)-delr(i)));
    
end
sum(Pl)
%% citation:
% Ashish Kumar (2021). Newton Raphson Load Flow for IEEE 14 Bus system (https://www.mathworks.com/matlabcentral/fileexchange/63717-newton-raphson-load-flow-for-ieee-14-bus-system), MATLAB Central File Exchange. Retrieved September 14, 2021.