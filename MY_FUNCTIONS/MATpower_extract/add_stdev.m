function [ mtxWithStdev ] = add_stdev( mtx, stdevPct, measCol )
%ADD_ERROR Summary of this function goes here
%   this adds a standard deviation of 1% of the measurement value
meas=mtx(:,measCol);
stdev=abs((stdevPct/100).*meas); % creates a vector with the numerical stdev values. Abs is used becaus stdev is not signed

mtxWithStdev=[mtx,stdev]; % the standard deviation of each meassurement is appended to the end of the matrix

end
