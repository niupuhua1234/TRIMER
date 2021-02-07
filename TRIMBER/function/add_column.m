function [trimer,varname] = add_column(trimer,name,vartype,lb,ub,obj,A)
% ADD_COLUMN  Add a column to a TRIMER model structure
%
%   [TRIMER,VARNAME] = ADD_COLUMN(TRIMER,N)
%   [TRIMER,VARNAME] = ADD_COLUMN(TRIMER,NAME,VARTYPE,LB,UB,OBJ,A)
%
%   Add a column to an existing structure, updating the corresponding
%   vectors.  The following default values are used:
%       NAME     'VARi', where i is the column index
%       VARTYPE  'b'
%       LB       0
%       UB       1
%       OBJ      0
%       A        column of zeros
%
%   If called as ADD_COLUMN(TRIMER,N), N columns are added with the default
%   values.  More than one name can be given; if only a single value is
%   given for the other arguments, it will be repeated for all new 
%   columns.
%
%   The return value VARNAME is the name(s) of the column(s) added.

[m,n] = size(trimer.A);
loc = n+1;

if nargin < 2 || isempty(name)
    name = {};
elseif isa(name,'double')
    name = zeros(1,name);  % holder until filled after calculating N
elseif isa(name,'char')
    name = {name};
end

if nargin < 3 || isempty(vartype), vartype = 'b'; end

if nargin < 4 || isempty(lb), lb = 0; end

if nargin < 5 || isempty(ub), ub = 1; end

if nargin < 6 || isempty(obj), obj = 0; end

if nargin < 7 || isempty(A), A = [];  end

N = max([length(name), ...
         length(lb),   ...
         length(ub),   ...
         length(obj),  ...
         size(A,2),    ...
         length(vartype)]);
         
if length(name) < N || isa(name,'double')
    % TODO:  add to names instead of replace
    name = array2names('VAR%i',loc:loc+N-1);
end
vartype = fill_to(vartype,N);
lb = fill_to(lb,N);
ub = fill_to(ub,N);
obj = fill_to(obj,N);
A = expand_to(A,[m N]);

locs = loc : loc + N-1;

trimer.varnames(locs) = name;
trimer.vartypes(locs) = vartype;
trimer.lb(locs) = lb;
trimer.ub(locs) = ub;
trimer.obj(locs) = obj;
trimer.A(:,locs) = A;

% expand the Q, Qd, and Qc fields
if isfield(trimer,'Q') && ~isempty(trimer.Q)
    trimer.Q(n+N,n+N) = 0;
end
if isfield(trimer,'Qd') && ~isempty(trimer.Qd)
    trimer.Qd(n+N,n+N) = 0;
end
if isfield(trimer,'Qc') && ~isempty(trimer.Qc)
    trimer.Qc.w(n+N) = 0;
    trimer.Qc.c(n+N) = 0;
end 


trimer = check_trimer(trimer);

if nargout > 1
    varname = name;
end



