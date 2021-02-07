function [names,idxs,logic] = convert_ids(all_names,ids,type)
% CONVERT_IDS  Create name, indices, and logical indices from an array
%
%   [NAMES,IDXS,LOGIC] = CONVERT_IDS(ALL_NAMES,IDS)
%   [IDS] = CONVERT_IDS(ALL_NAMES,IDS,TYPE)
%
%   Creates a cell of names (NAMES), linear indices (IDXS), and logical
%   indices (LOGIC) for the values IDS in the cell ALL_NAMES.  IDS can be 
%   any of the forms previously mentioned.
%
%   If three inputs are given, then only a single value is returned.  The
%   type of id returned is determined by type:
%       'name'    -> NAMES
%       'index'   -> IDXS
%       'logical' -> LOGIC
%
%   Example:
%       all_names = {'a','b','c','d'};
%       ids = [2 4];
%       [n,i,l] = convert_ids(all_names,ids)
%           n = {'b','d'}
%           i = [2 4]
%           l = [0 1 0 1]
%
%       idx = convert_ids(all_names,'c','index')
%           idx = 3

if isa(ids,'logical')
    logic = ids;
    names = all_names(logic);
    idxs  = find(logic);
elseif isa(ids,'double')
    idxs  = ids;
    logic = ismember(1:length(all_names),idxs);
    names = all_names(idxs);
else
    names = ids;
    logic = ismember(all_names,names);
    [~,idxs] = ismember(names,all_names);
end

if nargin == 3
    assert(nargout == 1,'only one output allowed when TYPE is given');
    switch validatestring(type,{'name','index','logical'})
        case 'name'
            names = names;
        case 'index'
            names = idxs;
        case 'logical'
            names = logic;
    end
end
