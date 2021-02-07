function [tf] = verify_sol(mip,sol,tol)
% VERIFY_SOL  Verify a solution is feasible
%
%   [TF] = VERIFY_SOL(MIP,SOL,TOL)
%
%   Verifies that the solution structure SOL satisfies the constraints
%   in MIP to within tolerance TOL.  The default tolerance is 1e-8.

if nargin < 3,  tol = 1e-8; end

switch class(sol)
    case 'double'
        x = sol(:);
    case 'struct'
        x = sol.x(:);
    otherwise
        error(['invalid class for argument sol: ' class(sol)]);
end

resid = mip.A * x - mip.b(:);

leq = (mip.ctypes == '<');
geq = (mip.ctypes == '>');
eq  = (mip.ctypes == '=');

tf = all(resid(leq) <= tol) ...
        && all(resid(geq) >= -tol) ...
        && all( abs(resid(eq)) <= tol );

tf = tf && ( all(x <= mip.ub(:) + tol) ...
                && all(x >= mip.lb(:) - tol) );

