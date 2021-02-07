function [x] = expand_to(x,N,dim)
% EXPAND_TO  Expand a vector to length N
%
%   [X] = EXPAND_TO(X,N,DIM)
%   [X] = EXPAND_TO(X,[M N])
%
%   Expands the vector X to length N by filling with zeros.  If a vector 
%   sizes is given, X is expanded as a M by N matrix.
%
%   If X is empty, a vector of zeros is created.  If DIM = 1 (default), 
%   the result is a column vector.  If DIM = 2, the a row vector is
%   returned.

if nargin < 4 || isempty(dim)
    dim = 1;
end
    
assert(nargin >= 2, 'at least two arguments required');

if length(N) == 2
    dims = N;
elseif dim == 1
    dims = [N,1];
else
    dims = [1,N];
end

if any(size(x) == 0)
    x = zeros(dims);
elseif ~all(size(x) == dims)
    x(dims(1),dims(2)) = 0;
end
