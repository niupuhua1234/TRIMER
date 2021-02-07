function [c,was_cell] = assert_cell(str)
% ASSERT_CELL  Assert that variable is a cell array.
%
%   [C,WAS_CELL] = ASSERT_CELL(STR)
%   
%   Checks that STR is a cell array.  If it is not (i.e., if it is a
%   single string), then converts it to a cell array of length one.
%
%   WAS_CELL is true if STR was originially a cell.

was_cell = isa(str,'cell');

if isempty(str) && ~isa(str,'char')
    c = {};
elseif ~isa(str,'cell')
    c = {str};
else
    c = str;
end
