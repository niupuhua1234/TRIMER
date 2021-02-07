function [ziped] = cellzipn(f,varargin)
% CELLZIP  Zip an unlimited number of cell arrays by a function
%
%   [ZIPED] = CELLZIP(F,...)
%
%   Computes ZIPED{i} = F(...), iterating element by element for each
%   input cell.

if isempty(varargin)
    ziped = [];
else
    N = length(varargin{1});
    if argmax(size(varargin),1) == 1
        ziped = cell(N,1);
    else
        ziped = cell(1,N);
    end
    
    extract = @(c,i) map(@(x) x{i},c);
    for i = 1 : N
        ziped{i} = f(extract(varargin,i));
    end
end

