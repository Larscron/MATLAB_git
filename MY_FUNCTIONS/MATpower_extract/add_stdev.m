function [ mtxWithStdev ] = add_stdev( mtx, stdevPct, measCol )
%ADD_ERROR Summary of this function goes here
%   this adds a standard deviation of 1% of the measurement value
meas=mtx(:,measCol);
stdev=abs((stdevPct/100).*meas); % creates a vector with the numerical stdev values. Abs is used becaus stdev is not signed
zeroPoints=find(stdev<1e-4);

mes_type=mtx(:,2); % hard coded where the measurements are
for i=1:length(zeroPoints)
    stdev(zeroPoints(i))=0.001;
    switch mes_type(i)
        case 1 %voltage meas
            stdev(zeroPoints(i))=3e-2;
        case 2 %Pinj
            stdev(zeroPoints(i))=1e-2;
        case 3 %Qimj
            stdev(zeroPoints(i))=1e-2;
        case 4 %Pflow
            stdev(zeroPoints(i))=8e-3;
        case 5 %Qflow
            stdev(zeroPoints(i))=8e-3;
        otherwise
            fprintf('SOMETHING HAS GONE TERRIBLY WRONG!!!\n')
            quit
    end
end
    
mtxWithStdev=[mtx,stdev]; % the standard deviation of each meassurement is appended to the end of the matrix

end

