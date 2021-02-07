function [tf] = near(x,y,tol)
% NEAR  Test if two values are close to each other
%
%   [TF] = NEAR(X,Y,TOL)
%
%   Tests if |X - Y| <= TOL.  If TOL is not given, the default is 1e-5.
%   If Y is not given, the default is 0.  X and Y can be single numbers,
%   vectors, or matrices.

if nargin < 3
    tol = 1e-5;
end

if nargin < 2
    y = 0;
end

tf = all(abs(x(:) - y(:)) <= tol);
