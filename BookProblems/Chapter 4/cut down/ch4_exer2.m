% Chapter 4, exercise 1
clc; clear all; close all;

%bus data       bus     P_bus(pu)   Q_bus(pu)   V_bus(pu)
bus_data=[      1       1.4902      -0.0368     1.0600;
                2       0.1830      0.11199     1.0450;
                3       -0.9421     0.0221      1.0100;
                4       -0.4779     0.0390      1.0250;
                5       -0.0760     -0.0158     1.0280;
                6       -0.11120    -0.1399     1.0700]; 
P_bus=bus_data(:,2);
Q_bus=bus_data(:,3);
V=bus_data(:,4);
bus_sh=zeros(length(P_bus),1);% just a placeholder in for this problem
thet=zeros(length(V),1);%start theta with 0 values

%               from(k) to(l)   P_kl(pu)    Q_kl(pu)    P_lk(pu)    Q_lk(pu)
branch_flow=[   1       2       1.0310      -0.0694     -1.0127     0.0669;
                1       5       0.4592      0.0325      -0.4489     -0.0436;
                2       3       0.6248      0.0474      -0.6078 	-0.0220;
                2       4       0.3428      -0.0057     -0.3366     -0.0118;
                2       5       0.2281      0.0112      -0.2253     -0.0400;
                3       4       -0.3343     0.0440      0.3418      -0.0381;
                4       5       -0.4831     0.0889      0.4862      -0.0793;
                5       6       0.1120      0.1470      -0.1120     -0.1399];
from=branch_flow(:,1);
to=branch_flow(:,2);
P_kl=branch_flow(:,3);
Q_kl=branch_flow(:,4);
P_lk=branch_flow(:,5);
Q_lk=branch_flow(:,6);

%               from(k) to(l)   r(pu)       x(pu)       b_sh(pu)    tap(pu)
branch_param=[  1       2       0.0194      0.0592      0.0264      1;
                1       5       0.0540      0.2230      0.0246      1;
                2       3       0.0470      0.1980      0.0219      1;
                2       4       0.0581      0.1763      0.0170      1;
                2       5       0.0570      0.1379      0.0173      1;
                3       4       0.0670      0.1710      0.0064      1;
                4       5       0.0134      0.0421      0.0000      1;
                5       6       0.0000      0.2520      0.0000      1.0730];

r=branch_param(:,3);
x=branch_param(:,4);
g=r./(r.^2+x.^2);% Conductance calculation
b=-x./(r.^2+x.^2 );% Susceptance calculation
y=g+1i*b;
% b_sh=branch_param(:,5);
b_sh=0.5*branch_param(:,5); % we multiply with 0.5 to get half the shunt values (B/2)

%a=1./branch_param(:,6);        
a=branch_param(:,6); %I'm not sure which one of these should be implemented, but neither method results in convergence
%a(8)=1;

N_bus=length(P_bus);
N_bra=length(from);
N_mes=N_bra*4+N_bus*3;
N_var=N_bus*2;

phi=zeros(N_bus,N_bus);% just a placeholder in for this problem
b_sh_bus=zeros(1,N_bus); % there are no bus shunts in this problem

y_ma=zeros(N_bus,N_bus);
b_sh_ma=zeros(N_bus,N_bus);
A=ones(N_bus,N_bus); % if we want to express a in a matrix, not used here

for t=1:N_bra
    y_ma(from(t),to(t))= y(t); %does not includ transformer pahse shift because it does not exist in this problem
    y_ma(to(t),from(t))=y_ma(from(t),to(t));
    
    
    b_sh_ma(from(t),to(t))= b_sh(t); %does not includ transformer pahse shift because it does not exist in this problem
    b_sh_ma(to(t),from(t))=b_sh_ma(from(t),to(t));
    
    A(from(t),to(t))=a(t);
end

%building the Y matrix from y_ma, B_sh_ma, and A (eq 3.40-3.41)
Y=zeros(N_bus,N_bus);
for k=1:N_bus
    for l=1:N_bus
        if k==l
            for m=1:N_bus
                Y(k,k)=Y(k,k)+A(k,m)^2*y_ma(k,m)+1i*b_sh_ma(k,m);
            end
            Y(k,k)=Y(k,k)+1i*b_sh_bus(k); % there is no bus shunt so this is just symbolic here
        else
            Y(k,l)=-A(k,l)*exp(-1i*phi(k,l))*y_ma(k,l);
        end
    end
end
          
% adding shuntdata to diagonals,     has no effect in this problem
%  for i = 1:N_bus
%       Y(i,i)= Y(i,i)+ bus_sh(i);  
%  end
 
G=real(Y);
B=imag(Y);

%declaring the initial state variables
X=[thet;V]; %flat start for theta and meassured values for the voltage
X(1)=[]; %removing theta_1 from the state variables, bacause it is the refference voltage


%order the measurements:
z_mea=zeros(N_mes,1);
pos=1;
Q_pos=2*N_bra; %the displacement of the Q measurements/equations  
for bra=1:N_bra
    z_mea(pos)=P_kl(bra);         %P_kl
    z_mea(pos+Q_pos)=Q_kl(bra);   %Q_kl
    pos=pos+1;                  
    z_mea(pos)=P_lk(bra);         %P_lk
    z_mea(pos+Q_pos)=Q_lk(bra);   %Q_lk
    pos=pos+1;
end
pos=pos+Q_pos;
for bus=1:N_bus
    z_mea(pos)=P_bus(bus);        %P_bus
    z_mea(pos+N_bus)=Q_bus(bus);  %Q_bus
    z_mea(pos+2*N_bus)=V(bus);    %V_bus
    pos=pos+1;
end

%standard deviations of each measurement:
stdev=abs(z_mea)./100;
stdev_sqr=stdev.^2;
R=diag(stdev_sqr); % declaring the R matrix
conv=1e-4;


% start the iterative solution process
for it=1:15
    %delta z:
    %fill the formula vector h()  
        h=zeros(N_mes,1);                                                           % declaring and clean the function vector h(x)
        P_bus_t=zeros(N_bus,1);                                                     % this is to test for faulty caclulations
        Q_bus_t=zeros(N_bus,1);                                                     % this is to test for faulty caclulations      
        pos=1;
        for t=1:N_bra% for the brach measurements
            thet_kl=thet(from(t))-thet(to(t));% yes, these are unessesary, but they make the equations nicer
            thet_lk=-thet_kl;
            phi_kl=phi(from(t),to(t));
            V_k=V(from(t)); V_l=V(to(t));
            
            %P_kl=          a_kl^2 V_k^2 g_kl           -a_kl V_k V_l g_kl  cos(?_kl+?)                         -a_kl*V_k       *V_l     *b_kl*sin(?_kl+?)
            h(pos)=         a(t)^2*V_k^2*g(t)    -a(t)*V_k*V_l*g(t)*cos(thet_kl+phi_kl)      -a(t)*V_k*V_l*b(t)*sin(thet_kl+phi_kl); % eq (3.32) looks good
            P_bus_t(from(t))=P_bus_t(from(t))+h(pos);                         
            
            %Q_kl=          -a_kl^2 V_k^2 (b_kl+b_kl^sh )           +a_kl V_k V_l b_kl*cos(?_kl+?)                      -a_kl*V_k*V_l*g_kl*sin(?_kl+?)
            h(pos+Q_pos)=  -a(t)^2*V_k^2*(b(t)+b_sh(t))     +a(t)*V_k*V_l*b(t)*cos(thet_kl+phi_kl)  -a(t)*V_k*V_l*g(t)*sin(thet_kl+phi_kl); % eq (3.33) looks good
            Q_bus_t(from(t))=Q_bus_t(from(t))+h(pos+Q_pos);
            pos=pos+1;
            
            %P_lk=          V_l^2 g_kl                  -a_kl V_k V_l g_kl*cos(?_lk-?)                          -a_kl*V_k*V_l*b_kl*sin*(?_lk-?)
            h(pos)=         V_l^2*g(t)             -a(t)*V_k*V_l*g(t)*cos(thet_lk-phi_kl)      -a(t)*V_k*V_l*b(t)*sin(thet_lk-phi_kl); % eq (3.34) looks good
            P_bus_t(to(t))=P_bus_t(to(t))+h(pos);
            
            %Q_lk=          
            h(pos+Q_pos)=  -V_l^2*(b(t)+b_sh(t))              +a(t)*V_k*V_l*b(t)*cos(thet_lk-phi_kl)  -a(t)*V_k*V_l*g(t)*sin(thet_lk-phi_kl); % eq (3.35) looks good
%             Q_bus_t(to(t))=Q_bus_t(to(t))+h(pos+Q_pos);
            pos=pos+1;
        end
        pos=pos+Q_pos;
        for k=1:N_bus %for the bus measurements
            for l=1:N_bus
                % test the if condition here
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 if k~=l
                    %P_bus_net
                    h(pos)=h(pos)               +V(k)*V(l)*(G(k,l)*cos(thet(k)-thet(l))+B(k,l)*sin(thet(k)-thet(l)));  % confirmed to be correct (3.46)
                    %Q_bus_net
                    h(pos+N_bus)=h(pos+N_bus)   +V(k)*V(l)*(G(k,l)*sin(thet(k)-thet(l))-B(k,l)*cos(thet(k)-thet(l)));  % confirmed to be correct (3.47)
%                 end
            end
            h(pos+2*N_bus)=V(k); % V_bus
            pos=pos+1;
        end
        %h(44)=Q_bus_t(6);% to heck if this helps the convergance
    
    
    % Jacobina matrix:
        %declare and clean the Jacobian matrix, it will be mostly 0s
        H=zeros(N_mes,N_var);
        pos=1;
        for t=1:N_bra% this loop will step through each branch
            thet_kl=thet(from(t))-thet(to(t));% yes, these are unessesary, but they make the equations nicer
            thet_lk=-thet_kl;
            phi_kl=phi(from(t),to(t));
            V_k=V(from(t)); V_l=V(to(t));
            %P_kl partial derivatives
                    %h(pos)=         a(t)^2*V_k^2*g(t)    -a(t)*V_k*V_l*g(t)*cos(thet_kl+phi(t))      -a(t)*V_k*V_l*b(t)*sin(thet_kl+phi(t));
                %Partial derrivative of P_kl/thet_k
                H(pos,from(t))= a(t)*V_k*V_l*g(t)*sin(thet_kl+phi_kl)   -a(t)*V_k*V_l*b(t)*cos(thet_kl+phi_kl); % eq (4.21) looks good
                %Partial derrivative of P_kl/thet_l
                H(pos,to(t))=  -a(t)*V_k*V_l*g(t)*sin(thet_kl+phi_kl)   +a(t)*V_k*V_l*b(t)*cos(thet_kl+phi_kl); % eq (4.21) looks good
                %Partial derrivative of P_kl/V_k
                H(pos,(from(t)+N_bus))= 2*a(t)^2*V_k*g(t)    -a(t)*V_l*g(t)*cos(thet_kl+phi_kl)     -a(t)*V_l*b(t)*sin(thet_kl+phi_kl); % eq (4.21) looks good
                %Partial derrivative of P_kl/V_l
                H(pos,(to(t)+N_bus))=                               -a(t)*V_k*g(t)*cos(thet_kl+phi_kl)   -a(t)*V_k*b(t)*sin(thet_kl+phi_kl); % eq (4.21) looks good

            %Q_kl partial derivatives
                    % h(pos+Q_pos)= -a(t)^2*V_k^2*(b(t)+b_sh(t))     +a(t)*V_k*V_l*b(t)*cos(thet_kl+phi(t))  -a(t)*V_k*V_l*g(t)*sin(thet_kl+phi(t));
                %Partial derrivative of Q_kl/thet_k
                H(pos+Q_pos,from(t))=   -a(t)*V_k*V_l*b(t)*sin(thet_kl+phi_kl)   -a(t)*V_k*V_l*g(t)*cos(thet_kl+phi_kl); % eq (4.23) looks good
                %Partial derrivative of Q_kl/thet_l
                H(pos+Q_pos,to(t))=      a(t)*V_k*V_l*b(t)*sin(thet_kl+phi_kl)   +a(t)*V_k*V_l*g(t)*cos(thet_kl+phi_kl); % eq (4.23) looks good
                %Partial derrivative of Q_kl/V_k
                H(pos+Q_pos,(from(t)+N_bus))= 2*a(t)^2*V_k*(b(t)+b_sh(t))    +a(t)*V_l*b(t)*cos(thet_kl+phi_kl)     -a(t)*V_l*g(t)*sin(thet_kl+phi_kl); % eq (4.23) looks good
                %Partial derrivative of Q_kl/V_l
                H(pos+Q_pos,(to(t)+N_bus))=                                         +a(t)*V_k*b(t)*cos(thet_kl+phi_kl)   -a(t)*V_k*g(t)*sin(thet_kl+phi_kl); % eq (4.23) looks good      
            pos=pos+1;

            %P_lk partial derivatives
                    %h(pos)=         V_l^2*g(t)             -a(t)*V_k*V_l*g(t)*cos(thet_lk-phi(t))      -a(t)*V_k*V_l*b(t)*sin(thet_lk-phi(t));
                %Partial derrivative of P_lk/thet_k
                H(pos,from(t))=-a(t)*V_k*V_l*g(t)*sin(thet_lk-phi_kl)   +a(t)*V_k*V_l*b(t)*cos(thet_lk-phi_kl); % eq (4.22) looks good 
                %Partial derrivative of P_lk/thet_l
                H(pos,to(t))=   a(t)*V_k*V_l*g(t)*sin(thet_lk-phi_kl)   -a(t)*V_k*V_l*b(t)*cos(thet_lk-phi_kl); % eq (4.22) looks good
                %Partial derrivative of P_lk/V_k
                H(pos,(from(t)+N_bus))=                             -a(t)*V_l*g(t)*cos(thet_lk-phi_kl)     -a(t)*V_l*b(t)*sin(thet_lk-phi_kl); % eq (4.22) looks good
                %Partial derrivative of P_lk/V_l
                H(pos,(to(t)+N_bus))=   2*V_l*g(t)             -a(t)*V_k*g(t)*cos(thet_lk-phi_kl)   -a(t)*V_k*b(t)*sin(thet_lk-phi_kl); % eq (4.22) looks good

            %Q_lk partial derivatives
                    %h(pos+2*N_bra-1)= -V_l^2*(b(t)+b_sh(t))              +a(t)*V_k*V_l*b(t)*cos(thet_lk-phi(t))  -a(t)*V_k*V_l*g(t)*sin(thet_lk-phi(t));
                %Partial derrivative of Q_lk/thet_k
                H(pos+Q_pos,from(t))=   a(t)*V_k*V_l*b(t)*sin(thet_lk-phi_kl)   +a(t)*V_k*V_l*g(t)*cos(thet_lk-phi_kl); % eq (4.24) looks good
                %Partial derrivative of Q_lk/thet_l
                H(pos+Q_pos,to(t))=    -a(t)*V_k*V_l*b(t)*sin(thet_lk-phi_kl)   -a(t)*V_k*V_l*g(t)*cos(thet_lk-phi_kl); % eq (4.24) looks good
                %Partial derrivative of Q_lk/V_k
                H(pos+Q_pos,(from(t)+N_bus))=                              +a(t)*V_l*b(t)*cos(thet_lk-phi_kl)      -a(t)*V_l*g(t)*sin(thet_lk-phi_kl); % eq (4.24) looks good
                %Partial derrivative of Q_lk/V_l
                H(pos+Q_pos,(to(t)+N_bus))=  -2*V_l*(b(t)+b_sh(t))    +a(t)*V_k*b(t)*cos(thet_lk-phi_kl)    -a(t)*V_k*g(t)*sin(thet_lk-phi_kl); % eq (4.24) looks good
            pos=pos+1;
        end
        pos=pos+2*N_bra;    
        for k=1:N_bus %steps through all buses
            for l=1:N_bus %steps through the buses connected to the k bus
                if k==l % if the k bus is connected to itself
                    for n=1:N_bus % step through all the buses
                        %P_k/thet_k
                        H(pos,k)=H(pos,k)               +V(k)*V(n)*(-G(k,n)*sin(thet(k)-thet(n))+B(k,n)*cos(thet(k)-thet(n)));% eq (4.27) looka good !!!!!!!!!!!!!!!!!!! shold they be inside or outside of this one?
                        %P_k/V_k
                        H(pos,k+N_bus)=H(pos,k+N_bus)+V(n)*(G(k,n)*cos(thet(k)-thet(n))+B(k,n)*sin(thet(k)-thet(n)));% eq (4.27) looka good !!!!!!!!!!!!!!!!!!! shold they be inside or outside of this one?
                        
                        %Q_k/thet_k 
                        H(pos+N_bus,k)=H(pos+N_bus,k)   +V(k)*V(n)*(G(k,n)*cos(thet(k)-thet(n))+B(k,n)*sin(thet(k)-thet(n)));% eq (4.28) looks good !!!!!!!!!!!!!!!!!!! shold they be inside or outside of this one?
                        %Q_k/V_k
                        H(pos+N_bus,k+N_bus)=H(pos+N_bus,k+N_bus)+V(n)*(G(k,n)*sin(thet(k)-thet(n))-B(k,n)*cos(thet(k)-thet(n)));% eq (4.28) looks good !!!!!!!!!!!!!!!!!!! shold they be inside or outside of this one?  
                    end
                    %add the las part of the P_k and Q_k with respect to
                    %V_k 
                    %And, according to eq 3.64-71 with respect to thet_k 
                    %P_k/thet_k
                    H(pos,k)=H(pos,k)               -V(k)^2*B(k,k);% eq (4.27) and eq (3.64) are different, this is 3.64
                    %Q_k/thet_k
                    H(pos+N_bus,k)=H(pos+N_bus,k)   -V(k)^2*G(k,k);% eq (4.28) and eq (3.66) are different, this is 3.66
                    
                    %P_k/V_k
                    H(pos,k+N_bus)=         H(pos,k+N_bus)          +V(k)*G(k,k);% eq (4.27) and 3.65 are different, this is 3.65!!!!!!!!!!!!!!!
                    %Q_k/V_k
                    H(pos+N_bus,k+N_bus)=   H(pos+N_bus,k+N_bus)    -V(k)*B(k,k);% eq (4.28) and 3.67 are different, this is 3.67!!!!!!!!!!!!!!!
                    
                else % if the k bus is not connected with itself, no sum is needed, so just find the partials
                   %partial derivative of P_k/thet_l
                   H(pos,l)=       V(k)*V(l)*(G(k,l)*sin(thet(k)-thet(l))  -B(k,l)*cos(thet(k)-thet(l)));% eq (4.27) looka good 
                   %partial derivative of P_k/V_l
                   H(pos,l+N_bus)= V(k)     *(G(k,l)*cos(thet(k)-thet(l))  +B(k,l)*sin(thet(k)-thet(l)));% eq (4.27) looka good
                   
                   %partial derivative of Q_k/thet_l
                   H(pos+N_bus,l)=         -V(k)*V(l)*(G(k,l)*cos(thet(k)-thet(l)) +B(k,l)*sin(thet(k)-thet(l)));% eq (4.28) looks good 
                   %partial derivative of Q_k/V_l
                   H(pos+N_bus,l+N_bus)=    V(k)*(G(k,l)*sin(thet(k)-thet(l))      -B(k,l)*cos(thet(k)-thet(l)));% eq (4.28) looks good 
                end       
            end
            H(pos+2*N_bus,k+N_bus)=1; % the derivatie of the voltage measurements are 1
            pos=pos+1;
        end
        
        H(:,1)=[];%removing the 1st column of H thar correspond to thet_1, the refference angle:
    
    % gain matrix   
        Gain=transpose(H)*inv(R)*H;
    % z_delta 
        z_del=z_mea-h;
    % the state variable correction vector
        X_del=(inv(Gain))*transpose(H)*inv(R)*z_del;
        
    % checking for convergance    
    if max(abs(X_del))>conv
        fprintf('iteration # %i compete, no convergence yet \n',it)
        %z_del'
        %X_del'
        max(abs(X_del))
%         if it>1
%             diff=H-old; diff;
%             old=H;
%             if it>2
%                 fprintf('the diffeence of difference is:\n')
%                 diffdiff=diff-olddiff;diffdiff
%                 olddiff=diff;
%             else
%                 fprintf('nothing to compare yet\n')
%                 olddiff=diff;
%             end
%         else
%             fprintf('nothing to compare yet\n')
%             old=H;
%         end
        
        X=X+X_del;
        thet(2:6)=X(1:5);
        V=X(6:11);
    else
        fprintf('We are done!!')
        it=1e9;
    end
end
%% P and Q bus calculations
P_bus=zeros(N_bus,1);
Q_bus=zeros(N_bus,1);
for bra=1:N_bra
    P_bus(from(bra))=P_bus(from(bra))+P_kl(bra);
    P_bus(to(bra))=P_bus(to(bra))+P_lk(bra);  
    
    Q_bus(from(bra))=Q_bus(from(bra))+Q_kl(bra);
    Q_bus(to(bra))=Q_bus(to(bra))+Q_lk(bra);      
end

%% P abd Q partail test calculations

P_part=zeros(N_bus,N_bus*2); % the N_bus*2 component is because there are two state varaibles for each bus
Q_part=zeros(N_bus,N_bus*2);
for bus=1:N_bus % we need to step through all buses
    for bra=1:N_bra % then we must check each branch at each bus
        thet_kl=thet(from(t))-thet(to(t));% yes, these are unessesary, but they make the equations nicer
        thet_lk=-thet_kl;
        phi_kl=phi(from(bra),to(bra));
        if bus==from(bra) % this is only done when thre are braches leaving this bus
            %P_k/thet_k
            P_part(bus,bus)=        a(bra)*V(from(bra))*V(to(bra))                  *(g(bra)*sin(thet_kl+phi_kl)-b(bra)*cos(thet_kl+phi_kl))    +P_part(bus,bus);
            %P_k/thet_l
            P_part(bus,to(bra))=    a(bra)*V(from(bra))*V(to(bra))                 *(-g(bra)*sin(thet_kl+phi_kl)+b(bra)*cos(thet_kl+phi_kl)); 
            %P_k/V_k
            P_part(bus,bus+N_bus)=2*a(bra)^2*V(from(bra))*g(bra)-a(bra)*V(to(bra))  *(g(bra)*cos(thet_kl+phi_kl)+b(bra)*sin(thet_kl+phi_kl))    +P_part(bus,bus+N_bus);
            %P_k/V_l
            P_part(bus,to(bra)+N_bus)=                          -a(bra)*V(from(bra))*(g(bra)*cos(thet_kl+phi_kl)+b(bra)*sin(thet_kl+phi_kl));
        
            % These are the Q partials
            %Q_k/thet_k
            Q_part(bus,bus)=    -a(bra)*V(from(bra))*V(to(bra))                     *(b(bra)*sin(thet_kl+phi_kl)+g(bra)*cos(thet_kl+phi_kl))    +Q_part(bus,bus);
            %Q_k/thet_l
            Q_part(bus,to(bra))= a(bra)*V(from(bra))*V(to(bra))                     *(b(bra)*sin(thet_kl+phi_kl)+g(bra)*cos(thet_kl+phi_kl));
            %Q_k/V_k !!! error in book????
            Q_part(bus,bus+N_bus)= -2*a(bra)^2*(b(bra)+b_sh(bra))+a(bra)*V(to(bra)) *(b(bra)*cos(thet_kl+phi_kl)-g(bra)*sin(thet_kl+phi_kl))    +Q_part(bus,bus+N_bus);
            %Q_k/V_l
            Q_part(bus,to(bra)+N_bus)=                  a(bra)*V(from(bra))         *(b(bra)*cos(thet_kl+phi_kl)-g(bra)*sin(thet_kl+phi_kl));
        end %end of if
    end %end of branch
    
    % incert the shunt for Q_k
    Q_part(bus,bus+N_bus)= -2*V(bus)*b_sh_bus(bus)  +   Q_part(bus,bus+N_bus); %there are no bus shunts, so this has no effect in this problem 
end %end of bus