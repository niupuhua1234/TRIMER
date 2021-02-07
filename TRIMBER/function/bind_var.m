function [trimer] = bind_var(trimer,vars,inds,varargin)
% BIND_VAR  Bind variables to a indicator variable
%
%   [TRIMER] = BIND_VAR(TRIMER,VARS,INDS,...params...)
%
%   For each variable v in VARS and corresponding indicator I in INDS,
%   adds constraints such that v=0 if I=0.
%
%   LB*I<=v<=UB*I.
%    
%   The bounds used when adding these rules are determined by the
%   parameters:
%       (default)  LB = min(TRIMER.lb), i.e. the lowest lower bound in the
%                  entire model.  UB = max(TRIMER.ub), the largest upper
%                  bound in the model.
%       'tight'    If true, the upper and lower bounds for v in trimber are used.
%                  This may be more numerically stable, but can add
%                  complications if the variable bounds are changed later,
%                  as these changes will not be reflected in the binding
%                  constraints.
%       'lb','ub'  Number specifying LB and UB; these are used for every
%                  variable.

p = inputParser;
p.addParamValue('tight',false);
p.addParamValue('lb',[]);
p.addParamValue('ub',[]);
p.parse(varargin{:});

tight = p.Results.tight;
default_lb = p.Results.lb;
default_ub = p.Results.ub;
if isempty(default_lb)
    default_lb = min(trimer.lb);
end
if isempty(default_ub)
    default_ub = max(trimer.ub);
end

assert(length(vars) == length(inds), ...
       'VARS and INDS must have the same length.');

% make sure we have names
[vars,var_idxs] = convert_ids(trimer.varnames,vars);
[inds,ind_idxs] = convert_ids(trimer.varnames,inds);

N = length(var_idxs);
A = zeros(2*N,size(trimer.A,2));
ctypes = repmat(' ',2*N,1);
for i = 1 : N
    if tight
        A(  i,[var_idxs(i) ind_idxs(i)]) = [1 -trimer.ub(var_idxs(i))];
        A(i+N,[var_idxs(i) ind_idxs(i)]) = [1 -trimer.lb(var_idxs(i))];
    else
        A(  i,[var_idxs(i) ind_idxs(i)]) = [1 -default_ub];
        A(i+N,[var_idxs(i) ind_idxs(i)]) = [1 -default_lb];
    end
    ctypes(  i) = '<';
    ctypes(i+N) = '>';
end

trimer = add_row(trimer,A,ctypes);

