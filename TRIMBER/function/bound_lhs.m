function [lb,ub] = bound_lhs(mip,rows)
% BOUND_MIP  Find bounds on LHS of a MIP
%
%   [LB,UB] = BOUND_LHS(MIP,ROWS)
%
%   LB(i) = min(A(rows(i),:)*x) for any lb <= x <= ub.
%   UB(i) = max(A(rows(i),:)*x) for any lb <= x <= ub.
%
%   If ROWS is not given, the bounds for all rows are computed.

if nargin < 2
    rows = 1 : size(mip.A,1);
end

N = length(rows);
lb = zeros(N,1);
ub = zeros(N,1);

for i = 1 : N
    lbs = mip.A(rows(i),:) .* mip.lb(:)';
    ubs = mip.A(rows(i),:) .* mip.ub(:)';
    lb(i) = sum(min(lbs,ubs));
    ub(i) = sum(max(lbs,ubs));
end