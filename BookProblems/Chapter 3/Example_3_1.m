%% Exampe 3.1
close all; clear all
P_2_sp=-.48;
conv=.001; %the convergen criteria is 0.001
V=[1,1]; %assuming both voltages start as 1PU
thet=zeros(2); %assume no intial phas lagg between busses
z_12=0.16+0.96j; %given impedance between the two lines
b_sh_12=.02j; %given shunt asmitance
y_12=1/z_12; %get the admitanc from the impedance

it=0; %we start on iteration 0

% creating the admitance matrix:
Y=zeros(2); %creates an empty matrix Y
for m=1:2
    for n=1:2
        if(m==n)% checks if we are on a diagonal
            Y(m,m)=y_12+b_sh_12;
        else % if we are not on a diagonal
            Y(m,n)=-y_12;
        end
    end 
end
%seperate the impedance matrix into real and imaginary part
G=real(Y);
B=imag(Y);

% calculating the power at 2
P=zeros(1,2);
k=2;
sum=0;
for l=1:k
    sum=sum+V(l)*(G(k,l)*cos(thet(k,l))+B(k,l)*sin(thet(k,l)));
end
P(k)=V(k)*sum;

P_2_cor=P_2_sp-P(2)

if(abs(P_2_cor)<conv)
    disp('we are done')
else 
    disp('another round is needed')
    it=it+1
end 

fprintf('in this case we need only one jacobian matrix cell \nyay:)\n')

J=V(2)*((-V(2)*B(2,2)-V(1)*(G(2,1)*sin(thet(2,1))-B(2,1)*cos(thet(k,l)))+(G(2,2)*sin(thet(2,2))+B(2,2)*cos(thet(2,2)))))

thet_cor=P_2_cor/J
thet(2)=thet(2)+thet_cor
% all good
%%

sum=0;
k=2;
for l=1:k
    sum=sum+V(l)*(G(k,l)*cos(thet(k,l))+B(k,l)*sin(thet(k,l)));
end
P(k)=V(k)*sum;
P
P_2_cor=P_2_sp-P(2)

if(abs(P_2_cor)<conv)
    disp('we are done')
else 
    disp('another round is needed')
    it=it+1
end 

%rouble here...
J=V(2)*((-V(2)*B(2,2)-V(1)*(G(2,1)*sin(thet(2,1))-B(2,1)*cos(thet(k,l)))+(G(2,2)*sin(thet(2,2))+B(2,2)*cos(thet(2,2)))))

%%
thet_cor=P_2_cor/J;
thet(2)=thet(2)+thet_cor

%% testing area
close all; clear all
