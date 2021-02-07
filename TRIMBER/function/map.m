function [mapped] = map(f,C)
% MAP  Generate a new list by applying a function
%
%   [MAPPED] = MAP(F,C) applies function handle F to each element in cell
%   C to create a new list MAPPED from the return values, i.e., 
%   MAPPED{i} = F(C{i})
%
%   MAP also works on arrays using ARRAYFUN.  For cells and arrays with
%   uniform return values, ARRAYFUN and CELLFUN are faster.

if isa(C,'cell')
    mapped = cellfun(f,C,'Uniform',false);
elseif isa(C,'double') || length(C) > 1
    mapped = arrayfun(f,C,'Uniform',false);
else
    mapped = cellfun(f,{C},'Uniform',false);
end


