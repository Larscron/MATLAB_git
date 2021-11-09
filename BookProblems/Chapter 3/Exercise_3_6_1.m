%% Exercises 3.1, #1
clc; close all; clear all

N_bus=14;%total number of buses
N_PQ=9;%number of PQ buses
N_PV=4;%number of PV buses
N_eq=2*N_PQ+N_PV;%number of equations to solve (Newton-Raphson technique)
S_base=100;%units of MVA

% branch shunts
bra_b_sh=zeros(N_bus);
bra_b_sh(1,2)=5.28; bra_b_sh(1,5)=4.92; bra_b_sh(2,3)=4.38; bra_b_sh(2,4)=3.40; bra_b_sh(2,5)=3.46; bra_b_sh(3,4)=1.28;
bra_b_sh=bra_b_sh/2;                                                       %The sunt is split to have half on each side of each line!!!!!!!!!!
bra_b_sh=bra_b_sh+transpose(bra_b_sh);%fill cells below the diagonal 
%converting to per unit:
bra_b_sh=bra_b_sh/S_base; %is this correct?

% bus shunt
bus_b_sh=zeros(1,14);
bus_b_sh(9)=19; %should this just be summed with the Q at bus 9?
%converting to per unit:
bus_b_sh=bus_b_sh/S_base; %is this correct?

% declaring the branch admittance matrix
bra_y=zeros(N_bus);%the impedances are converted directly to admittances through y=(1/z)
bra_y(1,2)=1/(1.938+5.917j); bra_y(1,5)=1/(5.403+22.304j);
bra_y(2,3)=1/(4.699+19.797j);bra_y(2,4)=1/(5.811+17.632j);bra_y(2,5)=1/(5.695+17.388j);
bra_y(3,4)=1/(6.701+17.103j);
bra_y(4,5)=1/(1.335+4.211j);bra_y(4,7)=1/(20.912j);bra_y(4,9)=1/(55.618j);
bra_y(5,6)=1/(25.202j);
bra_y(6,11)=1/(9.498+19.890j);bra_y(6,12)=1/(12.291+25.581j);bra_y(6,13)=1/(6.615+13.027j);
bra_y(7,8)=1/(17.615j);bra_y(7,9)=1/(11.001j);
bra_y(9,10)=1/(3.181+8.45j);bra_y(9,14)=1/(12.711+27.038j);
bra_y(10,11)=1/(8.205+19.207j);
bra_y(12,13)=1/(22.092+19.988j);
bra_y(13,14)=1/(17.093+34.802j);
bra_y=bra_y+transpose(bra_y);%fill cells below the diagonal 
%converting to per unit:
bra_y=bra_y*S_base; %multiplying with S_base rather than dividing each impedance term by S_base

%declaring the turns ratio matrix
V_1=132e3;
V_2=33e3;
a_tran=V_2/V_1;

a=ones(N_bus);
% a(4,7)=a_tran*0.978; %there are only 3 transformers, and I have assumed that the tap value is assosiated with the low-voltage side
% %    a(7,4)=1/a(4,7); %This is not the case! the transformer is only
% %    signifficant when moving grom 4 to 7, not the other way around
% a(4,9)=a_tran*0.969;
% a(5,6)=a_tran*0.932;

a(4,7)=1*0.978; %there are only 3 transformers, and I have assumed that the tap value is assosiated with the low-voltage side
%    a(7,4)=1/a(4,7); %This is not the case! the transformer is only
%    signifficant when moving grom 4 to 7, not the other way around
a(4,9)=1*0.969;
a(5,6)=1*0.932;


% The admittance matrix
Y=zeros(N_bus); %creates an empty matrix Y
for k=1:N_bus
    for l=1:N_bus
        if(k==l)% checks if we are on a diagonal
            sum=0;
            for t=1:N_bus
                sum=sum+a(k,t)^2*bra_y(k,t)+bra_b_sh(k,t)*1j;
            end         
            Y(k,k)=bus_b_sh(k)*1j+sum;
            
        else % if we are not on a diagonal
            Y(k,l)=-a(k,l)*bra_y(k,l);% there are no transformer phase angles, so this is not included in the formula
        end
    end 
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of admittance matrix setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%apparent power vector
S_bus=zeros(1,N_bus);
S_bus(1)=0;
S_bus(2)=18.3;
S_bus(3)=-94.2;
S_bus(4)=-47.8+3.9j;
S_bus(5)=-7.6-1.6j;
S_bus(6)=-11.2;
S_bus(7)=0;
S_bus(8)=0;
S_bus(9)=-29.5-16.6j+19j;
S_bus(10)=-9-5.8j;
S_bus(11)=-3.5-1.8j;
S_bus(12)=-6.1-1.6j;
S_bus(13)=-13.5-5.8j;
S_bus(14)=-14.9-5j;
%converting to per unit:
S_bus=S_bus/S_base;


%nominal voltage vector
V_nom=zeros(1,N_bus);
V_nom(1:5)=V_1;
V_nom(6:14)=V_2;


% Newton-Raphson technique

%declaring symbolic and numeric admittance matrices
Y_mag_sym=sym('Y_mag',[N_bus,N_bus]);
Y_ang_sym=sym('Y_ang',[N_bus,N_bus]);
Y_mag_num=abs(Y);
Y_ang_num=angle(Y);

%declaring symbolic and numeric voltage magnitudes and angles
V_mag_sym=sym('V_mag',[N_bus,1]);
V_ang_sym=sym('V_ang',[N_bus,1]);
V_mag_num=ones(N_bus,1);
V_mag_num(1)=1.060; V_mag_num(2)=1.045; V_mag_num(3)=1.010; V_mag_num(6)=1.070; V_mag_num(8)=1.090; %input the know V magnitudes
V_ang_num=zeros(N_bus,1);%no known V angles (except V_1=0), so initialized with zeros

%declaring the variable matrix numerically and numerically
X_sym=[V_ang_sym(2);V_ang_sym(3);V_ang_sym(4);V_ang_sym(5);V_ang_sym(6);
    V_ang_sym(7);V_ang_sym(8);V_ang_sym(9);V_ang_sym(10);V_ang_sym(11);
    V_ang_sym(12);V_ang_sym(13);V_ang_sym(14);
    V_mag_sym(4);V_mag_sym(5);V_mag_sym(7);V_mag_sym(9);V_mag_sym(10);
    V_mag_sym(11);V_mag_sym(12);V_mag_sym(13);V_mag_sym(14);];

X_num=[V_ang_num(2);V_ang_num(3);V_ang_num(4);V_ang_num(5);V_ang_num(6);
    V_ang_num(7);V_ang_num(8);V_ang_num(9);V_ang_num(10);V_ang_num(11);
    V_ang_num(12);V_ang_num(13);V_ang_num(14);
    V_mag_num(4);V_mag_num(5);V_mag_num(7);V_mag_num(9);V_mag_num(10);
    V_mag_num(11);V_mag_num(12);V_mag_num(13);V_mag_num(14);];

%declaring the specified bus-powers
P_sp=real(S_bus);
Q_sp=imag(S_bus);

C=[P_sp(2);P_sp(3);P_sp(4);P_sp(5);P_sp(6);P_sp(7);P_sp(8);P_sp(9);P_sp(10);
    P_sp(11);P_sp(12);P_sp(13);P_sp(14);
    Q_sp(4);Q_sp(5);Q_sp(7);Q_sp(9);Q_sp(10);Q_sp(11);Q_sp(12);Q_sp(13);Q_sp(14);];

%the power equations are defined
F=zeros(N_eq,1); %declares an empty matrix/vector to hold the power functions 
F_max=zeros(N_eq,1);
F=sym(F); %defines the function matrix as symbolic
F_sym=F;
t=0;
for n=2:N_bus
    for m=1:N_bus
        %real power functions
        F(n-1)=F(n-1)+V_mag_sym(n)*V_mag_sym(m)*Y_mag_num(n,m)*cos(Y_ang_num(n,m)-V_ang_sym(n)+V_ang_sym(m));
        
        if(n==4||n==5||n==7||n==9||n==10||n==11||n==12||n==13||n==14) %reactive power, but only for PQ buses
            F(N_bus+t)=F(N_bus+t)-V_mag_sym(n)*V_mag_sym(m)*Y_mag_num(n,m)*sin(Y_ang_num(n,m)-V_ang_sym(n)+V_ang_sym(m));
        end
    end
    if(n==4||n==5||n==7||n==9||n==10||n==11||n==12||n==13||n==14) %reactive power, but only for PQ buse
        t=t+1;
    end
end
%substituting in all the known values
F=subs(F,V_mag_sym(1),V_mag_num(1));
F=subs(F,V_mag_sym(2),V_mag_num(2));
F=subs(F,V_mag_sym(3),V_mag_num(3));
F=subs(F,V_mag_sym(6),V_mag_num(6));
F=subs(F,V_mag_sym(8),V_mag_num(8));
F=subs(F,V_ang_sym(1),0);

% thaking the jacobina:
J=jacobian(F,X_sym);
conv=1;
t=0;

while conv>1e-5
    J_it=J; %declare functions and jacobian that will be used for iterations
    F_it=F;
    
    for n=1:N_eq %sub in numeric values of variables
        J_it=subs(J_it,X_sym(n),X_num(n));
        F_it=subs(F_it,X_sym(n),X_num(n));
    end

    % Evaluate jacobian and power functions
    J_it=eval(J_it);
    F_it=eval(F_it);
    
    P_change=C-F_it;
    
    conv=max(abs(P_change));
    if conv>1e-5
        X_change=J_it\P_change; %apparently this is better than individually taking the inverse of J_it
        X_num=X_num+X_change; %finding the new variable values
        
        fprintf('we are on iteration: %i\n',t)
        fprintf('the constant values are:\n')
        disp(C')
        fprintf('the calculated values are:\n')
        disp(F_it')
        fprintf('the state variable errores are are:\n')
        disp(X_change')
        t=t+1;
    else
        disp('We are done!!')
        disp('the finall X values were:')
        disp(X_sym.')
        disp(X_num')
        fprintf('this took %i itterations\n',t)
    end
    
end 

% re-introducing the stat variables into the correct places
V_ang_num(2)=X_num(1);
V_ang_num(3)=X_num(2);
V_ang_num(4)=X_num(3);
V_ang_num(5)=X_num(4);
V_ang_num(6)=X_num(5);
V_ang_num(7)=X_num(6);
V_ang_num(8)=X_num(7);
V_ang_num(9)=X_num(8);
V_ang_num(10)=X_num(9);
V_ang_num(11)=X_num(10);
V_ang_num(12)=X_num(11);
V_ang_num(13)=X_num(12);
V_ang_num(14)=X_num(13);
V_mag_num(4)=X_num(14);
V_mag_num(5)=X_num(15);
V_mag_num(7)=X_num(16);
V_mag_num(9)=X_num(17);
V_mag_num(10)=X_num(18);
V_mag_num(11)=X_num(19);
V_mag_num(12)=X_num(20);
V_mag_num(13)=X_num(21);
V_mag_num(14)=X_num(22);





