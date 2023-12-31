function y = mrdivide(X,Y)

if (isa(Y,'sdpvar'))
    if Y.dim(1)*Y.dim(2) == 1
        y = X*Y^-1;
        return
    else
        error('Division of matrix variables not possible (maybe you meant ./).')
    end
end

try
    
    % Quick exit on scalar division
    if length(Y)==1
        if isinf(Y)
            y = zeros(X.dim);
        else
            y = X;
            y.basis = y.basis/Y;
        end
        y.conicinfo = [0 0];
        return
    end

    lmi_variables = getvariables(X);
    nv = length(lmi_variables);
    y  = X;
    n = X.dim(1);
    m = X.dim(2);    
    if (n==1) & (m==1) % SPECIAL CODE FOR FREAKY MATLAB BUG        
        y.basis = sparse(reshape(reshape(full(X.basis(:,1)),n,m)/Y,n*m,1));    
        for i = 1:nv
            temp = reshape(full(X.basis(:,i+1)),n,m)/Y;
            y.basis(:,i+1) = sparse(temp(:));   
        end;
    else
        y.basis = reshape(reshape(X.basis(:,1),n,m)/Y,n*m,1);    
        for i = 1:nv
            temp = reshape(X.basis(:,i+1),n,m)/Y;
            y.basis(:,i+1) = temp(:);   
        end;   
    end
    
    y.dim(1) = size(temp,1);
    y.dim(2) = size(temp,2);    
catch
    error(lasterr);
end

