function [err] = totalError(X,Y)

    D = abs(X-Y);
    err = sum(D(:));

end