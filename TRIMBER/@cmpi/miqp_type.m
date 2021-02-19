function [qtype] = miqp_type(mip)
% MIQP_TYPE  Return the quadratic type of a MIQP problem
%
%   [QTYPE] = MIQP_TYPE(MIP)
%
%   Return value:
%       'Q'   Standard form quadratic
%       'Qd'  Quadratic difference (x_i - x_j)^2
%       'Qc'  Quadratic constant (x_i - c)^2
%       []    Not quadratic

if isfield(mip,'Q') && ~isempty(mip.Q) && nnz(mip.Q) > 0
    qtype = 'Q';
elseif isfield(mip,'Qd') && ~isempty(mip.Qd)
    qtype = 'Qd';
elseif isfield(mip,'Qc') && ~isempty(mip.Qc)
    qtype = 'Qc';
else
    qtype = [];
end
