%% Chapter 3, Exercise 2
clc; close all; clear all
N_bus=5;
S_base=1e6;
V_base=12.66e3;
Z_base=V_base^2/S_base;

P=[0,1000,900,1200,800]*1000;
Q=[0,600,400,800,600]*1000;
S=(P+Q.*j)/S_base;

R=[0,.322,.493,.366,.3811];
X=[0,.27,.2511,.1864,.1941];
Z=(R+X.*j)/Z_base;

V=ones(1,5);
V_next=ones(1,5);
I_sent=zeros(1,N_bus);
V_dif=1;
n=0;

while V_dif>0.00001
    n=n+1;
    I_dem=conj(S./V); %we have no shunts
    fprintf('for iteration #%i\n',n)
    I_dem;

    for k=N_bus-1:-1:1 %the backwards sweep to find the current output of each bus
        I_sent(k)=sum(I_dem(k+1:N_bus));
    end
    I_sent;

    for k=2:N_bus %forwards sweep to 
        V_next(k)=V(k-1)-Z(k)*I_sent(k-1);
    end
    V_dif=max(abs(V-V_next))
    V=V_next;
    I_sent;
    V;
end

V_print='the nodal complex voltages ended up beeing:\n';
for n=1:N_bus
    str=num2str(V(n));
    V_print=append(V_print,str,'\n');
end
fprintf(V_print)
    
