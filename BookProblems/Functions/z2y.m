function [y_matrix] = z2y(z_matrix)
    %Z2Y Summary of this function goes here
    %   Detailed explanation goes here
    
    dim=length(z_matrix);    
    y=zeros(dim);
    for m=1:dim
        for n=1:dim
            if z_matrix(m,n)==0
                y(m,n)=0;
            else
                y(m,n)=1/z_matrix(m,n);
            end
        end
    end
    y_matrix=y;
end

