function [err] = MSE(X,Y)
%MSE Mean square error computation for two matrices
%  
   
    D = abs(X-Y).^2;
    err = sum(D(:))/numel(X);
    
end

