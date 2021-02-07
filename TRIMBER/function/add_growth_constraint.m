function [trimer] = add_growth_constraint(trimer,val)
% ADD_GROWTH_CONSTRAINT  Add minimum growth constraint to a model.
%
%   [TRIMER,SOL] = ADD_GROWTH_CONSTRAINT(TRIMER,VAL,...params...)
%  
%   Adds a constraint that requires a minimum flux through the objective
%   reaction.
%
%   Inputs
%   TRIMER   TRIMER model structure.
%   VAL     Constraining value (see the 'valtype' parameter for details).
%
%   Outputs
%   TRIMER   TRIMER model with growth constraint added.
%   SOL     CMPI solution object from the FBA calculation.
%
%   Parameters
%   'ctype'     Character indicating the type of constraint to add:
%                   '>'  -->  v_obj >= VAL  (default)
%                   '='  -->  v_obj  = VAL

if nargin < 2
    error('two inputs required');
end

p = inputParser;
p.addParamValue('ctype','>');
p.addParamValue('valtype','frac');
p.parse(varargin{:});

value = val;
trimer = add_row(trimer,1);
trimer.A(end,:) = trimer.obj';
trimer.b(end) = value;
trimer.ctypes(end) = p.Results.ctype;
trimer.rownames{end} = 'GROWTH';

