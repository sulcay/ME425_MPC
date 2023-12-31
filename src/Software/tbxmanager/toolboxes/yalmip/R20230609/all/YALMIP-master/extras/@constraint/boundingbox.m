function varargout = boundingbox(F,ops,x,fun_eval)

%BOUNDINGBOX Computes bounding box of a constraint
%
% If only one output is requested, only the symbolic model is returned
%  B = boundingbox(F)
%
% If three outputs are requested, numerical values are returned too
%  [B,L,U] = boundingbox(F)
%
% A second argument can be used to specificy solver settings
%  B = boundingbox(F,sdpsettings('solver','cplex'))
%
% A third argument can (should!) be used to obtain a bounding box in a
% particular set of variables, i.e. the bounding box of a projection.
% Unless you specify which variables you are computing the bounding box
% w.r.t, it will be very hard for you to understand which variables the
% values in L and U relate to
%  [B,L,U] = boundingbox(F,[],x)
% B will now be the box [L <= x <= U] (infinite bounds not included)
%
% See also VERTEX, CHEBYBALL, POLYTOPE

if nargin < 2
    ops = sdpsettings('verbose',0);
end

if nargin < 3
    [B,L,U,sols,vals] = boundingbox(lmi(F),ops);
elseif nargin < 4
    [B,L,U,sols,vals] = boundingbox(lmi(F),ops,x);
else
    [B,L,U,sols,vals] = boundingbox(lmi(F),ops,x,fun_eval);
end
switch nargout
    case 0
        B
    case 1
        varargout{1} = B;
    case 2
        varargout{1} = B;
        varargout{2} = L;
    case 3
        varargout{1} = B;
        varargout{2} = L;
        varargout{3} = U;
    case 5
        varargout{1} = B;
        varargout{2} = L;
        varargout{3} = U; 
        varargout{4} = sols; 
        varargout{5} = vals; 
    otherwise
end
