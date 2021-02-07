function [mip] = bounds_to_constraints(mip,varargin)
% BOUNDS_TO_CONSTRAINTS  Convert LB and UB to inequalities
%
%   Adds inequalities x_i <= UB(i) and x_i >= LB(i), effectively moving
%   these constraints from UB and LB and into A.
%
%   Parameters
%   'bounds'    If 'all' (default), all bounds are re-written.  If 'equal',
%               only bounds where UB(i) == LB(i) are moved.
%   'equal_tol' Tolerance for determining if UB(i) == LB(i).  Default is
%               1e-10.
%   'max_bound' Replacement upper and lower bounds.  Default is 1e10 for 
%               upper boundsand -1e10 for lower bounds.
%
%   This function is used as a work-around for the GLPK solver where 
%   bounds of the form UB(i) == LB(i) are not allowed.

p = inputParser;
p.addParamValue('bounds','all');
p.addParamValue('max_bound',1e10);
p.addParamValue('equal_tol',1e-10);
p.parse(varargin{:});

max_bound = p.Results.max_bound;

if strcmpi(p.Results.bounds,'all')
    idxs = 1:length(mip.ub);
elseif strcmpi(p.Results.bounds,'equal')
    idxs = find(abs(mip.ub - mip.lb) <= p.Results.equal_tol);
end

N = length(idxs);
roff = size(mip.A,1);
mip = add_row(mip,2*N);

for i = 1 : N
    mip.A(roff+i,idxs(i)) = 1;
    mip.ctypes(roff+i) = '<';
    mip.b(roff+i) = mip.ub(idxs(i));
    
    mip.A(roff+i+N,idxs(i)) = 1;
    mip.ctypes(roff+i+N) = '>';
    mip.b(roff+i) = mip.lb(idxs(i));
end

mip.ub(idxs) =  max_bound;
mip.lb(idxs) = -max_bound;
