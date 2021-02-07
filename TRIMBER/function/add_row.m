function [trimer] = add_row(trimer,A,ctype,b,name,ind,indtype)
% ADD_ROW  Add a row to a TRIMER model structure
%
%   [TRIMER] = ADD_ROW(TRIMER,N)
%   [TRIMER] = ADD_ROW(TRIMER,A,CTYPE,B,NAME,IND,INDTYPE)
%
%   Add a row to an existing structure, updating the corresponding 
%   vectors.  The following default values are used:
%       A        row of zeros
%       CTYPE    '='
%       B        0
%       NAME     'ROWi', where i is the row index
%       IND      0
%       INDTYPE  ' '
%
%   If called as ADD_ROW(TRIMER,N), N rows are added with the default
%   values.  If only a single value is given for each argument, it will
%   be repeated for all new rows.

[m,n] = size(trimer.A);
loc = m+1;

if nargin < 2 || isempty(A)
    A = [];
elseif length(A) == 1
    if nargin == 2 && A == 0
        return
    else
        A = zeros(A,n);
    end
end

if nargin < 3 || isempty(ctype), ctype = '='; end

if nargin < 4 || isempty(b), b = 0; end

if nargin < 5 || isempty(name)
    name = {};
elseif isa(name,'char')
    name = {name};
end


N = max([length(name),  ...
         length(ctype), ...
         length(b),     ...
         size(A,1)]);

A = expand_to(A,[N n]);
ctype = fill_to(ctype,N);
b = fill_to(b,N);

if length(name) < N
    if length(name) == 1
        name = array2names([name{1},'%i'],loc:loc+N-1);
    else
        name = array2names('ROW%i',loc:loc+N-1);
    end
end

locs = loc + (0:N-1);

trimer.A(locs,:) = A;
trimer.ctypes(locs) = ctype;
trimer.b(locs) = b;
trimer.rownames(locs) = name;

trimer = check_trimer(trimer);
