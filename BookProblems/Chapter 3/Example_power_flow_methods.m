%% settining up the system parameters
clc; close all; clear all
N_bus=4;

% declaring impedance matrix
Z=zeros(N_bus);
Z(1,1)=1j;
Z(1,2)=0.25j; Z(2,1)=Z(1,2);
Z(1,3)=0.4j; Z(3,1)=Z(1,3);
Z(2,2)=0.8j;
Z(2,4)=0.4j; Z(4,2)=Z(2,4);
Z(3,4)=0.5j; Z(4,3)=Z(3,4);

% ampping the admitances
y=z2y(Z);
% creating the admitance matrix:
Y=y2Y(y)

%% Gauss-Seidel power flow solution
clc; close all; clear all
N_bus=3;
V=ones(1,3);
V(1)=1.02+0j; V(3)=1.03

P_sch=[0,-2,1.5]; %incerting the scheduled powers
P=P_sch;% first power guess is the same as the scheduled ones

Q_sch=[0,-0.5,0];
Q=Q_sch;% first power guess is the same as the scheduled ones

Z=zeros(N_bus);
z(1,2)=.02+.06j; z(2,1)=z(1,2);
z(1,3)=.0059+.0235j; z(3,1)=z(1,3);
z(2,3)=.0055+.0183j; z(3,2)=z(2,3);

% ampping the admitances
y=z2y(z)
% creating the admitance matrix:
Y=y2Y(y)

for p=1:6 %decides how many iterations
    
    %PQ bus calculations
    k=2;
    sum=0;
    for l=1:N_bus
        if k~=l
            sum=sum+Y(k,l)*V(l);
        end
    end                         
    V(k)=1/Y(k,k)*((P(k)-Q(k)*j)/(conj(V(k)))-sum);

    %PV bus calculations
    k=3;
    %Q is not given and thus has to be calculated
    sum=0;
    for l=1:N_bus
        if k~=l
            sum=sum+Y(k,l)*V(l);
        end
    end
    Q(k)=-imag(conj(V(k))*(V(k)*Y(k,k)+sum));
    %now that Q is nown we can move on to calculating the voltage
    sum=0;
    for l=1:N_bus
        if k~=l
            sum=sum+Y(k,l)*V(l);
        end
    end                         
    V_tmp=1/Y(k,k)*((P(k)-Q(k)*j)/(conj(V(k)))-sum); %Temporary V(k) value
    V_imag=imag(V_tmp); %the magnitude of V(k) is already deffined, so we will only keep the imaginary part of V(k), because it is the smalles one
    V_real=sqrt(abs(V(k))^2-V_imag^2);

    V(3)=V_real+V_imag*1i;
    mag=abs(V(k));
    %ang=180/pi*angle(V(k))

    fprintf('for iteration #%i we get:\nV_2 = ',p)
    polPrint(V(2))
    fprintf('V_3 = ')
    polPrint(V(3))
end

%% Newton-Raphson technique
clc; close all; clear all;
N_bus=4;
%we start with a completed admittance matrix
Y=[15-55j,  -5+15j, -10+40j; 
   -5+15j,  20-65j, -15+50j;
  -10+40j, -15+50j, 25-90j];
%declaring symbolic and numeric admittance matrices
Y_mag_sym=sym('Y_mag',[N_bus,N_bus]);
Y_ang_sym=sym('Y_ang',[N_bus,N_bus]);
Y_mag_num=abs(Y);
Y_ang_num=angle(Y);

%declaring symbolic and numeric voltage magnitudes and angles
V_mag_sym=sym('V_mag',[3,1]);
V_ang_sym=sym('V_ang',[3,1]);
V_mag_num=[1;1;1.03]; %these matrices include initial guesses and known values
V_ang_num=[0;0;0];

%declaring the variable matrix numerically and numerically
X_sym=[V_ang_sym(2);V_ang_sym(3);V_mag_sym(2)];
X_num=[V_ang_num(2);V_ang_num(3);V_mag_num(2)];

%declaring the specified bus-powers
P_2_sp=-2; P_3_sp=1.5; Q_2_sp=-0.5;
C=[P_2_sp;
   P_3_sp;
   Q_2_sp];


%the power equations are defined
F=zeros(3,1); %declares an empty matrix/vector to hold the power functions 
F=sym(F); %defines the function matrix as symbolic
for n=1:3 
    %P_2
    F(1)=F(1)+V_mag_sym(2)*V_mag_sym(n)*Y_mag_sym(2,n)*cos(Y_ang_sym(2,n)-V_ang_sym(2)+V_ang_sym(n));
    %P_3
    F(2)=F(2)+V_mag_sym(3)*V_mag_sym(n)*Y_mag_sym(3,n)*cos(Y_ang_sym(3,n)-V_ang_sym(3)+V_ang_sym(n));
    %Q_2
    F(3)=F(3)-V_mag_sym(2)*V_mag_sym(n)*Y_mag_sym(2,n)*sin(Y_ang_sym(2,n)-V_ang_sym(2)+V_ang_sym(n));
end

% taking the jacobian of the symbolic power functions
J=jacobian(F,X_sym);%we now get the jacobian of the functions F with respect to X

% we now have all the symbolic functions we will need, but before starting
% the itterative solution loop we want to substitute in the constant that
% will not change with sucessive itterations

for m=1:3
    for n=1:3
        %subs in the value of the admittance magnitudes and angles
        F=subs(F,Y_mag_sym(m,n),Y_mag_num(m,n));
        F=subs(F,Y_ang_sym(m,n),Y_ang_num(m,n));
        J=subs(J,Y_mag_sym(m,n),Y_mag_num(m,n));
        J=subs(J,Y_ang_sym(m,n),Y_ang_num(m,n));
    end
end
% all constants have now been substituted into F and J
  
fprintf('| iteration |    V1    |         V2         |         V3         |\n')

for t=1:4 % will be replaced by a while loop, but itterations are simpler for now
    %fprintf('This is itteration #%i\n', t)
    
    F_it=F; %declare the F and J itteration matrices form the original matrices
    J_it=J;
    for m=1:3
        %subs in the most recent value of the voltage magnitudes and angles
        F_it=subs(F_it,V_mag_sym(m),V_mag_num(m));
        F_it=subs(F_it,V_ang_sym(m),V_ang_num(m));
        J_it=subs(J_it,V_mag_sym(m),V_mag_num(m)); 
        J_it=subs(J_it,V_ang_sym(m),V_ang_num(m));
    end
    F_it=eval(F_it); %evaluate and get the numerical value of F and J
    J_it=eval(J_it);
    
    
    P_change=C-F_it;
    
    X_change=J_it\P_change; %apparently this is better than individually taking the inverse of J_it
    %J_it_inv=inv(J_it);
    %X_change=J_it_inv*P_change;
    
    X_num=X_num+X_change; %finding the new variable values
    
    V_ang_num(2)=X_num(1); %using the new variables to update the voltage magnitueds and angles
    V_ang_num(3)=X_num(2);
    V_mag_num(2)=X_num(3);
    
    V_ang_num_deg=V_mag_num*180/pi;
    
    fprintf('%7d%11d %d%16d %d %13d %d\n',t,V_mag_num(1),V_ang_num(1),V_mag_num(2),V_ang_num(2),V_mag_num(3),V_ang_num(3))
%     F_it
%     J_it
%     
%     V_mag_num'
%     V_ang_num'*180/pi
end





