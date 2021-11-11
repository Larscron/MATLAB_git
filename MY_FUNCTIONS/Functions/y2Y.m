function [Y_matrix] = y2Y(y_matrix,b_sh)
    %y2Y  transforms and admittance-location matrix to a admittance matrix
    %   Detailed explanation goes here
                                                            %to be added if (y_matrix)==y_matrix
    dim=length(y_matrix);% since this is a square matrix we can use the length to give the dimentions of the new matris
    %fprintf('the number of inputs are %i\n',nargin)
    
    if 2>nargin
        b_sh=zeros(dim);
        %fprintf('there was no shunt input\n')
    end
    %creating the admitance matrix:
    Y=zeros(dim); %creates an empty matrix Y
    for m=1:dim
        for n=1:dim
            if(m==n)% checks if we are on a diagonal
                for i=1:dim                                                 %i don't know it if is an issue to have thsi be i, since it could be comfused with current
                    Y(m,m)=Y(m,m)+y_matrix(m,i)+b_sh(m,i)*1j;
                end
            else % if we are not on a diagonal
                Y(m,n)=-y_matrix(m,n);
            end
        end 
    end
    Y_matrix=Y;
end

