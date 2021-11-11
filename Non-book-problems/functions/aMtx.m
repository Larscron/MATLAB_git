function [aMtx] = aMtx(lineMtx)
%SMALLY takes in the line data of a bus system and returns the small y
%matrix
%   the line data must be in the order:
%   |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%   |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |

fb = lineMtx(:,1);                  % From bus number...
tb = lineMtx(:,2);                  % To bus number...
a= 1./lineMtx(:,6);                 % of some reasonn the taps must be expressed as 1/a
%a= lineMtx(:,6);

nbus = max(max(fb),max(tb));        % no. of buses...
nbra = length(fb);                  % no. of branch

aMtx=ones(nbus,nbus);

for t=1:nbra % convert parameters to matrices
    aMtx(fb(t),tb(t))= a(t); %does not includ transformer pahse shift because it does not exist in this problem
end

end