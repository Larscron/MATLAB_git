function [ measurement_with_error ] = add_nor_error(measurement,percentage_error)
%ADD_NOR_ERROR adds a normally distributed percentage error

%   This function adds normally distributed percentage error to the input measurement

pct_error=percentage_error.*0.01; % turn into percentage
pct_meas=measurement.*pct_error;
pct_meas=abs(pct_meas); % the normrnd only takes positive values, so we need the absolute value of the percentage

measurement_with_error=normrnd(measurement,pct_meas);
end

