clear all; close all; clc;

%       Bus No  Vm      Theta   P_gen   Q_gen   P_1oad  Q_1oad  Qmax    Qmin    shunt(pu)   Type
busdata= [ 1    1.060   0       0       0       0       0       0       0       0           1;
           2    1.045   0       40.0    45.41   21.7    12.7    -40     50      0           2;
           3    1.010   0       0       25.28   94.2    19.0    0       40      0           2;
           4    1.000   0       0       0       47.8    -3.9    0       0       0           3;
           5    1.000   0       0       0       7.6     1.6     0       0       0           3;
           6    1.070   0       0       13.62   11.2    7.5     -6.0    24.0    0           3;
           7    1.000   0       0       0       0       0       0       0       0           3;
           8    1.090   0       0       18.24   0       0       -6.0    24.0    0           3;
           9    1.000   0       0       0       29.5    16.6    0       0       0.190       3;
           10   1.000   0       0       0       9.0     5.8     0       0       0           3;
           11   1.000   0       0       0       3.5     1.8     0       0       0           3;
           12   1.000   0       0       0       6.1     1.6     0       0       0           3;
           13   1.000   0       0       0       13.5    5.8     0       0       0           3;
           14   1.000   0       0    	0       14.9    5.0     0       0       0           3];
                                    % sometimes Pg and Pl are combined,

                                            % same for Qg and Ql
Base_Mva=100;                               % The power base is 100 MVA
V=busdata(:,2);                             % Bus voltages
thet=busdata(:,3);                          % Bus angles
Pg=busdata(:,4)/Base_Mva;                   % Real power generated at bus
Qg=busdata(:,5)/Base_Mva;                   % Reactive power generated bus
Pl=busdata(:,6)/Base_Mva;                   % Real power consumed at bus
Ql=busdata(:,7)/Base_Mva;                   % Reactive power consumed bus
P_sp= Pg - Pl;                              % Calculating net real power at bus
Q_sp= Qg - Ql;                              % Calculating net reactive power at bus
Qmax = busdata(:,8)/Base_Mva;               % Minimum Reactive Power Limit..
Qmin = busdata(:,9)/Base_Mva;               % Maximum Reactive Power Limit..

b_sh_bus=-busdata(:,10);                    % i do not know why the example wanted b_shunt to be negative? Probably because it is capacitive and not inductive!
bus_type=busdata(:,11);                     % All of the bus types

pv = find(bus_type == 2 | bus_type == 1);   % List of the PV buses, including the V-theta/refference bus
pq = find(bus_type == 3);                   % List of the PQ Buses..
npv = length(pv);                           % No. of PV buses..
npq = length(pq);                           % No. of PQ buses.

% Line Data for IEEE 14 bus system
%        Line No   From to  R           X           B/2         T    
linedata= [ 1       1   2   0.01938     0.05917     0.02640     1;
            2       2   3   0.04699     0.19797     0.02190     1;
            3       2   4   0.05811     0.17632     0.01870     1;
            4       1   5   0.05403     0.22304     0.02460     1;
            5       2   5   0.05695     0.17388     0.01700     1;
            6       3   4   0.06701     0.17103     0.01730     1;
            7       4   5   0.01335     0.04211     0.00640     1;
            8       5   6   0.00000     0.25202     0.00000     0.932;
            9       4   7   0.00000     0.20912     0.00000     0.978;
            10      7   8   0.00000     0.17615     0.00000     1;
            11      4   9   0.00000     0.55618     0.00000     0.969;
            12      7   9   0.00000     0.11001     0.00000     1;
            13      9   10  0.03181     0.08450     0.00000     1;
            14      6   11  0.09498  	0.19890     0.00000     1;
            15      6   12  0.12291     0.25581     0.00000     1;
            16      6   13  0.06615     0.13027     0.00000     1;
            17      9   14  0.12711     0.27038     0.00000     1;
            18      10  11  0.08205     0.19207     0.00000     1;
            19      12  13  0.22092     0.19988     0.00000     1;
            20      13  14  0.17093     0.34802     0.00000     1];

fb = linedata(:,2);         %from bus
tb = linedata(:,3);         %to bus
R =  linedata(:,4);         %branch resistance
X =  linedata(:,5);         %branch reluctance
bra_sh = linedata(:,6);     %branch shunt
T =  linedata(:,7);         %tap values
a = 1./T;                   %the relative transformer value                 % could be that this should be changed
%a=T;
z =  R+sqrt(-1)*X;          %branh impedance
y = 1./z;                   %branch admitance


nbus=length(busdata(:,1));
nbra = length(fb); % no of branch
nmes=nbra*4+nbus*3;
nvar=nbus*2;

% constructing the constant matrixes
y_matrix=zeros(nbus,nbus);
b_sh_matrix=zeros(nbus,nbus);
A=ones(nbus,nbus); % if we want to express a in a matrix, not used here
for t=1:nbra % convert parameters to matrices
    y_matrix(fb(t),tb(t))= y(t); %does not includ transformer pahse shift because it does not exist in this problem
    y_matrix(tb(t),fb(t))=y_matrix(fb(t),tb(t));
        
    b_sh_matrix(fb(t),tb(t))= bra_sh(t); %does not includ transformer pahse shift because it does not exist in this problem
    b_sh_matrix(tb(t),fb(t))=b_sh_matrix(fb(t),tb(t));
    
    A(fb(t),tb(t))=a(t);
end

% all the variables are now matrices
y=y_matrix;
g=real(y);
b=imag(y);
b_sh=b_sh_matrix;
a=A;

% measurements

%declaring the initial state variables
X=[thet;V]; %flat start for theta and meassured values for the voltage
X(1)=[]; %removing theta_1 from the state variables, bacause it is the refference voltage

% calculating the net power of each bus through the power flow equations
it = 1;         % we are starting at itteration 1
tol = 1e-5;     % we declare the 
error = 1;      % we set the error to an arbitrary value that will be greater that the tollerance

while (error>tol)
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    for k=1:nbus
        for m=1:nbus
            thet_km=thet(k)-thet(m);
            thet_mk=-thet_km;   
            P(k)=P(k)+( (a(k,m)*V(k))^2*g(k,m)              -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)   );
            Q(k)=Q(k)+(-(a(k,m)*V(k))^2*(b(k,m)+b_sh(k,m))  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)   );
        end
        Q(k)=Q(k)-(V(k)^2*b_sh_bus(k));
    end

    %difference form spesified power
    diffP=P_sp-P;
    diffQ=Q_sp-Q;
    
    % Checking Q-limit violations..
    if it <= 7 && it > 2    % Only checked up to 7th iterations..
        for n = 2:nbus      % we are going through all but the 1st bus
            if bus_type(n) == 2
                QG = Q(n)+Ql(n);    %retrevign the generated reactive power form the calculated net Q value
                if QG < Qmin(n)     % if the generated reactive power is less than the treshold
                    V(n) = V(n) + 0.01;        % we increase the bus voltage -> more reactive power
                elseif QG > Qmax(n) % if the generated reactive power is more than the treshold
                    V(n) = V(n) - 0.01;        % we reduse the bus voltage -> less reactive power
                end
            end
         end
    end

    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if bus_type(i) == 3
            dQ(k,1) = diffQ(i);
            k = k+1;
        end
    end

    dP = diffP(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    %constructing the jacobian from the partial derivatives of P and Q
% Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for f = 1:(nbus-1)% f for function
        k = f+1;
        for v = 1:(nbus-1)% v for variable
            m = v+1;
            if m == k % if we are taking the derivative with respect to the theta of the bus, thet_k
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    J1(f,v) = J1(f,v) + (   (a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)   );
                end
            else        % if we are taking the partial derivative with respect to thet_m, there is no need to sum
                J1(f,v) =               (  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)   );
            end
        end
    end
    
    J2 = zeros(nbus-1,npq);
    for f = 1:(nbus-1)% f for function
        k = f+1;
        for v = 1:npq% v for variable
            m = pq(v);
            if m == k % if we are taking the derivative with respect to the voltage of the bus, V_k
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    J2(f,v) = J2(f,v) +( 2*a(k,m)^2*V(k)*g(k,m)     -(a(k,m)*a(m,k)*V(m))*g(k,m)*cos(thet_km)   -(a(k,m)*a(m,k)*V(m))*b(k,m)*sin(thet_km)   );
                end
            else        % if we are taking the partial derivative with respect to V_m, there is no need to sum
                J2(f,v) = (                                         -(a(k,m)*a(m,k)*V(k))*g(k,m)*cos(thet_km)   -(a(k,m)*a(m,k)*V(k))*b(k,m)*sin(thet_km)   );
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for f = 1:npq% f for function
        k = pq(f);
        for v = 1:(nbus-1)% v for variable
            m = v+1;
            if m == k % if we are taking the derivative with respect to the theta of the bus, thet_k
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    J3(f,v) = J3(f,v) + (  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)  );
                end
            else        % if we are taking the partial derivative with respect to thet_m, there is no need to sum
                J3(f,v) =               (  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km)  +(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)   );
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for f = 1:npq   % f for function
        k = pq(f);
        for v = 1:npq% v for variable
            m = pq(v);
            if m == k % if we are taking the derivative with respect to the voltage of the bus, V_k
                for m = 1:nbus
                    thet_km=thet(k)-thet(m);
                    J4(f,v) = J4(f,v) +( -2*a(k,m)^2*V(k)*(b(k,m)+b_sh(k,m))    +(a(k,m)*a(m,k)*V(m))*b(k,m)*cos(thet_km)   -(a(k,m)*a(m,k)*V(m))*g(k,m)*sin(thet_km)   );
                end
                J4(f,v) = J4(f,v) -(2*V(k)*b_sh_bus(k));
            else        % if we are taking the partial derivative with respect to V_m, there is no need to sum
                J4(f,v) = (                                                     +(a(k,m)*a(m,k)*V(k))*b(k,m)*cos(thet_km)   -(a(k,m)*a(m,k)*V(k))*g(k,m)*sin(thet_km)   );
            end
        end
    end

    J=[J1,J2;J3,J4];
    
    dX = inv(J)*M;           % Correction vector for state variables
    dTheta = dX(1:nbus-1);   % Change in Voltage Angle..
    dV = dX(nbus:end);       % Change in Voltage Magnitude..
    
    % updating the state variables 
    thet(2:nbus)=dTheta+thet(2:nbus);       % updating all but the 1st voltage angle
    
    k = 1;
    for i = 2:nbus 
        if bus_type(i) == 3                 % updating all of the PQ voltage magnitudes
            V(i) = dV(k)+V(i);
            k = k+1;
        end
    end
    
    it = it + 1;                            % increasing itteration count
    error = max(abs(M));                    % the current maximum error is derived form the missmatch vector
end
    
fprintf 'we are done\n'

printStateVar(V,thet)
    
    
%% Power flow calculations
% this section will be used

Pflow=zeros(nbus,nbus);
Qflow=zeros(nbus,nbus);

for k=1:nbus
    for m=1:nbus
        thet_km=thet(k)-thet(m);
        
        Pflow(k,m)= (a(k,m)*V(k))^2*g(k,m)              -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*sin(thet_km);
        

        Qflow(k,m)=-(a(k,m)*V(k))^2*(b(k,m)+b_sh(k,m))  +(a(k,m)*a(m,k)*V(k)*V(m))*b(k,m)*cos(thet_km)  -(a(k,m)*a(m,k)*V(k)*V(m))*g(k,m)*sin(thet_km);
        
        if ((Qflow(k,m) ~= 0) && (k<m) && (0))
            fprintf('the Q flow from %g to %g was:\n',k,m)
            fprintf('%f\n',Qflow(k,m))
        end

    end
end