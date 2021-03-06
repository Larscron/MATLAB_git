function [busDataMtx] = busData(nbuses)
%LINEDATA Summary of this function goes here
%   Detailed explanation goes here
    switch nbuses
        case 14
            %          |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi |  Qmin | Qmax |
            busDataMtx=[ 1     1    1.060   0       0     0     0     0       0       0;
                         2     2    1.045   0      40   42.4  21.7   12.7    -40     50;
                         3     2    1.010   0       0   23.4  94.2   19.0     0      40;
                         4     3    1.0     0       0     0   47.8   -3.9     0       0;
                         5     3    1.0     0       0     0    7.6    1.6     0       0;
                         6     2    1.070   0       0   12.2  11.2    7.5    -6      24;
                         7     3    1.0     0       0     0    0.0    0.0     0       0;
                         8     2    1.090   0       0   17.4   0.0    0.0    -6      24;
                         9     3    1.0     0       0     0   29.5   16.6     0       0;
                         10    3    1.0     0       0     0    9.0    5.8     0       0;
                         11    3    1.0     0       0     0    3.5    1.8     0       0;
                         12    3    1.0     0       0     0    6.1    1.6     0       0;
                         13    3    1.0     0       0     0   13.5    5.8     0       0;
                         14    3    1.0     0       0     0   14.9    5.0     0       0 ];
        otherwise
            disp ('We do not have this matrix saved')
            busDataMtx = 0;
    end
end

