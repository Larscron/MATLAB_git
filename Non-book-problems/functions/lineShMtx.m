function [lineShMtx] = lineShMtx(lineMtx)
%SMALLY takes in the line data of a bus system and returns the small y
%matrix
%   the line data must be in the order:
%   |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%   |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |

fb = lineMtx(:,1);                  % From bus number...
tb = lineMtx(:,2);                  % To bus number...
bsh= lineMtx(:,5);

nbus = max(max(fb),max(tb));        % no. of buses...
nbra = length(fb);                  % no. of branch

bsh_mtx=zeros(nbus,nbus);

for t=1:nbra % convert parameters to matrices
    bsh_mtx(fb(t),tb(t))= bsh(t); %does not includ transformer pahse shift because it does not exist in this problem
    bsh_mtx(tb(t),fb(t))= bsh_mtx(fb(t),tb(t));
end

lineShMtx=bsh_mtx;

end