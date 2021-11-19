function [smallyMtx] = smally(lineMtx)
%SMALLY takes in the line data of a bus system and returns the small y
%matrix
%   the line data must be in the order:
%   |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%   |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |

fb = lineMtx(:,1);                  % From bus number...
tb = lineMtx(:,2);                  % To bus number...
r=lineMtx(:,3);
x=lineMtx(:,4);
z = r + 1i*x;                       % Z matrix...
y=1./z;

nbus = max(max(fb),max(tb));        % no. of buses...
nbra = length(fb);                  % no. of branch

y_mtx=zeros(nbus,nbus);

for t=1:nbra % convert parameters to matrices
    y_mtx(fb(t),tb(t))= y(t); %does not includ transformer pahse shift because it does not exist in this problem
    y_mtx(tb(t),fb(t))= y_mtx(fb(t),tb(t));
end

smallyMtx=y_mtx;

end

