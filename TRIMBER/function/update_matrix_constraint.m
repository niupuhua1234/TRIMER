function [trimer] = update_matrix_constraint(trimer,linalgs,row_names,matching)
%  
% UPDATE_MATRIX_CONSTRAINTupdate the constrain for given constrain names
%   Inputs
%   TRIMER   TRIMER model structure.
%   linalgs  constraint object  
%   rownames  constraint names
%
%   Outputs
%   TRIMER   TRIMER model with new constraints added.
%
%   Parameters
%   matching     select all the rownames that contain the input
%                name(default) or just select the one that match the name.
if nargin < 3
    row_names = {}; 
elseif isa( row_names,'char')
    row_names = {row_names};
end
if nargin < 4, matching = true;end 
% check that a TRIMER model was given (and convert if COBRA)
linalgs = assert_cell(linalgs);
for i = 1 : length(linalgs)
    
    var= reshape([linalgs{i}.vars{:}],[],1);
    coefs=[linalgs{i}.coefs{:}];

   
    if isempty(row_names)
       roff=ones(length(trimer.rownames),1);
    else
        if matching && ischar(row_names{i})
             roff=startsWith(trimer.rownames,row_names{i});
        else
            [roff,~]=ismember(trimer.rownames,row_names{i});
        end
    end
                        
    [loc,~] = ismember(trimer.varnames,var);
    
    trimer.A(roff,loc) =  coefs;
    trimer.ctypes(roff) = linalgs{i}.op;
    trimer.b(roff) = linalgs{i}.rhs;
end



