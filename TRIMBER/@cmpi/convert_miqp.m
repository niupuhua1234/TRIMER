function [mip] = convert_miqp(mip)
% CONVERT_MIQP  Prepare a MIQP for solution
%
%   [MIP] = CONVERT_MIQP(MIP)
%
%   Converts a general MIQP problem to standard form:
%       min sense*( x'*Q*x + obj*x )
%       subject to
%           A*x (<=/=/>=) b
%           lb <= x <= ub
%
%   General MIQP problems are defined by three fields:
%       Q     Matrix with standard quadratic weights for x'*Q*x.
%
%       Qd    If Qd(i,j) = k, then the objective entry is 
%             k*(x(i) - x(j))^2.  A difference variable d is added such
%             that d = (x(i) - x(j)), and the entry Q(d,d) = k.
%
%       Qc.w  If Qc.w(i) and Qc.c(i) are nonzero (Qc.w and Qc.c are 
%       Qc.c  vectors the same size as x), then the objective entry is
%             Qc.w(i)*(x(i) - Qc.c(i))^2, i.e., the weighted least-square
%             distance between x(i) and a constant c(i).

has_Q  = ~isempty(cmpi.check_field('Q',mip));
has_Qd = ~isempty(cmpi.check_field('Qd',mip));
has_Qc = ~isempty(cmpi.check_field('Qc',mip));

assert(sum([has_Q has_Qd has_Qc]) <= 1, ...
       'only one of Q, Qd, and Qc can be specified');

[m,n] = size(mip.A);
 
if has_Q
    % move all nonzero elements in Q to the lower half
    mip.Q  = tril(mip.Q,-1) + triu(mip.Q)';
elseif has_Qd
    % move all nonzero elements in Q to the lower half
    Qd = mip.Qd;
    Qd = tril(Qd,-1) + triu(Qd)';
    
    mip.Q = spalloc(n,n,nnz(Qd));
    [I,J,w] = find(Qd);
    Nadd = length(w);
    mip = add_column(mip,Nadd,'c');
    mip = add_row(mip,Nadd);
    for i = 1 : Nadd
        mip.A(m+i,[I(i) J(i) n+i]) = [1 -1 1];
        mip.Q(n+i,n+i) = w(i);
        lb = min(mip.lb([I(i) J(i)]));
        ub = max(mip.ub([I(i) J(i)]));
        mip.lb(n+i) = -(ub - lb);
        mip.ub(n+i) = ub - lb;
    end
elseif has_Qc
    mip.Q = spalloc(n,n,nnz(mip.Qc.w));
    Nadd = length(find(mip.Qc.w));
    mip = add_column(mip,Nadd,'c');
    mip = add_row(mip,Nadd);
    rowidx = m;
    colidx = 0;
    for i = 1 : length(mip.Qc.w)
        if mip.Qc.w(i) ~= 0
            rowidx = rowidx + 1;
            colidx = colidx + 1;
            mip.A(rowidx,[n+colidx i]) = [1 -1];
            mip.b(rowidx) = -mip.Qc.c(i);
            mip.Q(n+colidx,n+colidx) = mip.Qc.w(i);
            bound = mip.ub(i) - mip.lb(i);
            mip.ub(n+colidx) =  bound;
            mip.lb(n+colidx) = -bound;
        end
    end
end

% make Q symmetric
mip.Q = tril(mip.Q) + tril(mip.Q,-1)';

% ensure that Q is sparse
mip.Q = sparse(mip.Q);
