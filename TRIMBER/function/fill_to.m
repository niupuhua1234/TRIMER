function [filled] = fill_to(x,N,default,dim)
% FILL_TO  Fill a short to empty vector
%
%   [FILLED] = FILL_TO(X,N,DEFAULT,DIM)
%
%   If X has only one element, replicate it to have N elements.  If X is
%   empty, use DEFAULT (DEFAULT = 0 if not specified).  Vector orientation
%   is determined by DIM:  1 (default) -> column, 2 -> row.

if nargin < 4 || isempty(dim)
    dim = 1;
end
assert(dim == 1 || dim == 2, 'DIM must be either 1 or 2.');

if nargin < 3 || isempty(default)
    default = 0;
end

if isempty(x)
    filled = default;
else
    filled = x;
end

if length(x) < N
    if dim == 1
        filled = repmat(filled,N,1);
    else
        filled = repmat(filled,1,N);
    end
end

    