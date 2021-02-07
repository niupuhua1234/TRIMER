function set_solver_option(option,val)
% SET_SOLVER_OPTION  Set a default solver option
%
%   SET_SOLVER_OPTION(OPTION,VAL) sets an option in the default option 
%   structure.
%
%   SET_SOLVER_OPTION(STRUCT) sets the default option structure to the
%   structure STRUCT.
%
%   For a description of the solver options, see the documentation
%   for SOLVE_MIP.

if nargin == 1
    cmpi.set_option(option);
else
    cmpi.set_option(option,val);
end
