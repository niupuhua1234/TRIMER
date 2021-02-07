function [trimer] = add_matrix_constraint(trimer,linalgs,row_names)
%   ADD_MATRIX_CONSTRAINT   add the constrain with  given constrain names  to model
%
%    [TRIMER,SOL] = ADD_MATRIX_CONSTRAINT (TRIMER,LINALG,ROW_NAMES)
%
%   Inputs
%   TRIMER   TRIMER model structure.
%   linalgs  constraint object  
%   rownames  constraint names
%
%   Outputs
%   TRIMER   TRIMER model with new constraints added.


if nargin < 3 || isempty(row_names), row_names = []; end

% check that a TRIMER model was given (and convert if COBRA)
linalgs = assert_cell(linalgs);
for i = 1 : length(linalgs)
    
    Nrows = size(linalgs{1}.rhs,1);
    Nprev = size(trimer.A,1);
    if isempty(row_names)
        trimer = add_row(trimer,Nrows);
    else        
        trimer = add_row(trimer,Nrows,[],[],row_names{i});
    end        
    roff = Nprev +  [1:Nrows];
    var= reshape([linalgs{i}.vars{:}],[],1);
    [tf,loc] = ismember(var,trimer.varnames);
    if any(~tf)
        trimer = add_column(trimer,var(tf));
        [~,loc] = ismember(linalgs{i}.vars,trimer.varnames);
    end
    coefs=[linalgs{i}.coefs{:}];
    trimer.A(roff,loc) =  coefs;
    trimer.ctypes(roff) = linalgs{i}.op;
    trimer.b(roff) = linalgs{i}.rhs;
end
