function [names] = array2names(fmt,array,dim)
% ARRAY2NAMES  Create a cell of names from an array of numbers
%
%   [NAMES] = ARRAY2NAMES(FMT,ARRAY,DIM)
%
%   Using a printf format string FMT that accepts a single integer, 
%   create a cell of names that are sequentially numbered by the values
%   in array.  If DIM = 1, the result is a column cell.  Otherwise, the
%   result is a row cell.
%
%   Example:
%       array2names('var%i',1:3)
%       ans = 
%           'var1'
%           'var2'
%           'var3'

if nargin < 3 || isempty(dim)
    dim = 1;
end

names = arrayfun(@(x) sprintf(fmt,x),array,'Uniform',false);
if dim == 1
    names = names(:);
end

