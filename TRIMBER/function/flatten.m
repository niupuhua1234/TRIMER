function [list] = flatten(lists,depth)
% FLATTEN  Flatten a cell of cells into a single cell.
%
%   [LIST] = FLATTEN(LISTS) flattens a cell of cells into a single cell.
%   If no element of LISTS is a cell, LIST = LISTS.
%
%   [LIST] = FLATTEN(LISTS,DEPTH) flattens to a depth of DEPTH.  The
%   default depth is 1.
%
%   >> A = FLATTEN({{1,2},3,{4}})
%   A = 
%       [1]  [2]  [3]  [4]

if nargin < 2
    depth = 1;
end

if depth > 1
    lists = flatten(lists,depth-1);
end

if ~any(cellfun(@iscell,lists))
    list = lists;
else
    %list = [lists{:}];
    % hack to avoid error with horzcat and expr objects
    n = 0;
    for i = 1 : length(lists)
        if isa(lists{i},'cell')
            n = n + length(lists{i});
        else
            n = n + 1;
        end
    end
    list = cell(1,n);
    idx = 0;
    for i = 1 : length(lists)
        if isa(lists{i},'cell')
            for j = 1 : length(lists{i})
                idx = idx + 1;
                list{idx} = lists{i}{j};
            end
        else
            idx = idx + 1;
            list{idx} = lists{i};
        end
    end
end
