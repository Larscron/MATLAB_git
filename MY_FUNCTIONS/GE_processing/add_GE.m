function [ meas_with_error, error_location ] = add_GE( meas, STdev )
%ADD_GE adds a GE to a random measurement within the given function
%   
loc = randi(length(meas));
switch randi(2)
    case 1
        %add + GE
        meas(loc)=meas(loc)+10*STdev(loc);
    case 2
        %add - GE
        meas(loc)=meas(loc)-10*STdev(loc);
end

meas_with_error = meas;
error_location = loc;
end

