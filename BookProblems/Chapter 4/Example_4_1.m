clc; clear all; close all

%setting up constants
N_bus=3;
N_var=N_bus*2-1;
a=ones(N_bus);%there are no transformers in this system
phi=0;% should phi thecnically also be a matrix?????????????????????????????????????????????????????????????????

g=zeros(N_bus);
g(1,2)=4.99913;
g(2,3)=1.13502;
g=g+g';%adds the conductance for the other way

b=zeros(N_bus);
b(1,2)=-15.26308;
b(2,3)=-4.78186;
b=b+b';

b_sh=zeros(N_bus);
b_sh(1,2)=0.02640;
b_sh(2,3)=0.02190;
b_sh=b_sh+b_sh';

y=g+b*1j;
Y=y2Y(y,b_sh);
G=real(Y);  B=imag(Y);

P_sent=zeros(N_bus); P_sent_st=zeros(N_bus);
P_sent(1,2)=0.8190;  P_sent_st(1,2)=0.005460;
P_sent(2,1)=-0.8070; P_sent_st(2,1)=0.005380;
P_sent(2,3)=0.9900;  P_sent_st(2,3)=0.006600;
P_sent(3,2)=-0.9420; P_sent_st(3,2)=0.006280;
P_sent_sym=sym('P_sent',[N_bus,N_bus]);
P_sym=sym('P',[N_bus,N_bus]);

Q_sent=zeros(N_bus); Q_sent_st=zeros(N_bus);
Q_sent(1,2)=-0.0110; Q_sent_st(1,2)=0.000073;
Q_sent(2,1)=-0.0120; Q_sent_st(2,1)=0.00080;
Q_sent(2,3)=0.3490;  Q_sent_st(2,3)=0.002327;
Q_sent(3,2)=-0.1900; Q_sent_st(3,2)=0.001267;
Q_sym=sym('Q',[N_bus,N_bus]);

P_bus_num=ones(1,N_bus);P_bus_st=zeros(1,N_bus);
P_bus_num(2)=0.1830;    P_bus_st(2)=0.001220;
P_bus_sym=sym('P',[N_bus,1]); %since P and Q bus contains equations it seemed better to arange them vertically

Q_bus_num=ones(1,N_bus);Q_bus_st=zeros(1,N_bus);
Q_bus_num(2)=0.3380;    Q_bus_st(2)=0.002253;
Q_bus_sym=sym('Q',[N_bus,1]);

V_num=ones(1,N_bus); V_st=zeros(1,N_bus);
V_num(1)=1.0600;        V_st(1)=0.007067;
V_num(3)=0.9450;        V_st(3)=0.006300;
V_sym=sym('V',[1,N_bus]);

thet_bus_num=zeros(1,N_bus);%we initially assume that the voltage phase angle between busse is 0 rad.
thet_bus_sym=sym('thet',[1,N_bus]);
thet_bus_sym(1)=0;%we establish the voltage angle at bus 1 to be 0

thet_sym=sym('thet',[N_bus,N_bus]); %this is the voltage angle differences between the busses (not all of them are relevant)
for k=1:N_bus                                                               % making thet matrix
    for l=1:N_bus
        thet_sym(k,l)=thet_bus_sym(k)-thet_bus_sym(l);
    end
end

%Generalized ? model power equations
for k=1:N_bus                                                               % branch calculations
    for l=1:N_bus
        %numeric
        %P(k,l)=a(k,l)^2*V(k)^2*g(k,l)-a(k,l)*V(k)*V(l)*g(k,l)*cos(thet(k,l)+phi)-a(k,l)*V(k)*V(l)*b(k,l)*sin(thet(k,l)+phi);
        %Q(k,l)=-a(k,l)^2*V(k)^2*(b(k,l)+b_sh(k,l))+a(k,l)*V(k)*V(l)*b(k,l)*cos(thet(k,l)+phi)-a(k,l)*V(k)*V(l)*g(k,l)*sin(thet(k,l)+phi);
        %symbolic
        P_sym(k,l)=a(k,l)^2*V_sym(k)^2*g(k,l)-a(k,l)*V_sym(k)*V_sym(l)*g(k,l)*cos(thet_sym(k,l)+phi)-a(k,l)*V_sym(k)*V_sym(l)*b(k,l)*sin(thet_sym(k,l)+phi); %does not contain error vector!
        Q_sym(k,l)=-a(k,l)^2*V_sym(k)^2*(b(k,l)+b_sh(k,l))+a(k,l)*V_sym(k)*V_sym(l)*b(k,l)*cos(thet_sym(k,l)+phi)-a(k,l)*V_sym(k)*V_sym(l)*g(k,l)*sin(thet_sym(k,l)+phi);
    end
end

for k=1:N_bus                                                               % bus calculations
    sum_P=0;
    sum_Q=0;
    for l=1:N_bus
        sum_P=sum_P+V_sym(l)*(G(k,l)*cos(thet_sym(k,l))+B(k,l)*sin(thet_sym(k,l))); %does not contain error vector!
        sum_Q=sum_Q+V_sym(l)*(G(k,l)*sin(thet_sym(k,l))-B(k,l)*cos(thet_sym(k,l)));
    %Q_2
    end
    P_bus_sym(k)=V_sym(k)*sum_P;
    Q_bus_sym(k)=Q_sym(k)*sum_P;
end
%
%the chosen variables X:
X_sym=sym('X',[1,N_var]);
X_sym(1)=thet_bus_sym(2);
X_sym(2)=thet_bus_sym(3);
X_sym(3)=V_sym(1);
X_sym(4)=V_sym(2);
X_sym(5)=V_sym(3);
X_num=[thet_bus_num(2),thet_bus_num(3),V_num(1),V_num(2),V_num(3)];

%the chosen functions F:
F_sym=sym('F',[12,1]);      z_mes=zeros(12,1);
F_sym(1)=P_sym(1,2);        z_mes(1)=P_sent(1,2);
F_sym(2)=P_sym(2,1);        z_mes(2)=P_sent(2,1);
F_sym(3)=P_sym(2,3);        z_mes(3)=P_sent(2,3);
F_sym(4)=P_sym(3,2);        z_mes(4)=P_sent(3,2);
F_sym(5)=Q_sym(1,2);        z_mes(5)=Q_sent(1,2);
F_sym(6)=Q_sym(2,1);        z_mes(6)=Q_sent(2,1);
F_sym(7)=Q_sym(2,3);        z_mes(7)=Q_sent(2,3);
F_sym(8)=Q_sym(3,2);        z_mes(8)=Q_sent(3,2);
F_sym(9)=P_bus_sym(2);      z_mes(9)=P_bus_num(2);
F_sym(10)=Q_bus_sym(2);     z_mes(10)=Q_bus_num(2);
F_sym(11)=V_sym(1);         z_mes(11)=V_num(1);
F_sym(12)=V_sym(3);         z_mes(12)=V_num(3);

R=diag([P_sent_st(1,2)^2,P_sent_st(2,1)^2,P_sent_st(2,3)^2,P_sent_st(3,2)^2, Q_sent_st(1,2)^2,Q_sent_st(2,1)^2,Q_sent_st(2,3)^2,Q_sent_st(3,2)^2, P_bus_st(2)^2,Q_bus_st(2)^2,V_st(1)^2,V_st(1)^2]);
%to use the jacobina, thet_sym bust be reorganized in terms of bus voltage angles

H=jacobian(F_sym,X_sym);%we now get the jacobian of the functions F with respect to X

conv=1;

while conv>0.001
    H_it=H;% declare the jacobian that will be solved for this itteration
    F_it=F_sym;

    H_it=subs(H_it,thet_bus_sym(2),X_num(1));                                   % H_it=subs(H_it,thet_bus_sym(2),thet_bus_num(2));
    H_it=subs(H_it,thet_bus_sym(3),X_num(2));                                   % H_it=subs(H_it,thet_bus_sym(3),thet_bus_num(3));
    H_it=subs(H_it,V_sym(1),X_num(3));                                          % H_it=subs(H_it,V_sym(1),V_num(1));
    H_it=subs(H_it,V_sym(2),X_num(4));                                          % H_it=subs(H_it,V_sym(2),V_num(2));
    H_it=subs(H_it,V_sym(3),X_num(5));                                          % H_it=subs(H_it,V_sym(3),V_num(3));
    F_it=subs(F_it,thet_bus_sym(2),X_num(1));                                   % F_it=subs(F_it,thet_bus_sym(2),thet_bus_num(2));
    F_it=subs(F_it,thet_bus_sym(3),X_num(2));                                   % F_it=subs(F_it,thet_bus_sym(3),thet_bus_num(3));
    F_it=subs(F_it,V_sym(1),X_num(3));                                          % F_it=subs(F_it,V_sym(1),V_num(1));
    F_it=subs(F_it,V_sym(2),X_num(4));                                          % F_it=subs(F_it,V_sym(2),V_num(2));
    F_it=subs(F_it,V_sym(3),X_num(5));                                          % F_it=subs(F_it,V_sym(3),V_num(3));

    H_it=eval(H_it); %evaluate and get the numerical value of H_it
    F_it=eval(F_it);

    %Calculate the gain matrix
    Gain=H_it'*inv(R)*H_it;

    %calculate the correction matrix
    X_cor=inv(Gain)*H_it'*inv(R)*[z_mes-F_it];

    %declare new X values
    X_num=X_num+X_cor'

    %thest for convergence
    fprintf('the deviation from known values are:\n')
    conv=max(abs(X_cor));

    %if no convergance the loop is repeted:
    % we could feed the new variable values inot the respective theta and V
    % matrices, but this is not necessary
   
    %X_num=[thet_bus_num(2),thet_bus_num(3),V_num(1),V_num(2),V_num(3)];
end

fprintf('final X values were:\n')
X_num
