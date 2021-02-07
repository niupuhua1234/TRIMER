function [val] = check_field(field,structure)
% CHECK_FIELD  Check if a field exists in a structure.
%
%   [VAL] = CHECK_FIELD(FIELD,STRUCTURE)
%   
%   Checks if FIELD is a defined field in STRUCTURE.  If defined, 
%   returns the value in STRUCTURE.FIELD.  If not definied, returns
%   the empty array.

if isfield(structure,field)
    val = structure.(field);
else
    val = [];
end