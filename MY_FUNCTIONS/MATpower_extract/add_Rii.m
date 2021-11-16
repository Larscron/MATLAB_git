function [ mtxWithRii ] = add_Rii(mtx)
%ADD_ERROR Summary of this function goes here
%   this adds a standard deviation of 1% of the measurement value
mtx_dim=size(mtx);
nmes=mtx_dim(1);
Rii=zeros(nmes,1);
mes_type=mtx(:,2);
for i=1:nmes
    switch mes_type(i)
        
        case 1 %voltage meas
            Rii(i)=9e-4;
        case 2 %Pinj
            Rii(i)=1e-4;
        case 3 %Qimj
            Rii(i)=1e-4;
        case 4 %Pflow
            Rii(i)=64e-6;
        case 5 %Qflow
            Rii(i)=64e-6;
        otherwise
            fprintf('SOMETHING HAS GONE TERRIBLY WRONG!!!\n')
            quit
    end
end
    
mtxWithRii=[mtx,Rii]; % the standard deviation of each meassurement is appended to the end of the matrix
fprintf('test')
end

