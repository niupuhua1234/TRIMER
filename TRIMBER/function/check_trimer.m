function [tiger] = check_tiger(tiger)
% CHECK_TIGER  Check size and orientation of TIGER fields
%
%   [TIGER] = CHECK_TIGER(TIGER)
%
%   Checks that fields in a TIGER model structure are the correct size 
%   and orientation.  Returns the validated TIGER model.

col_fields = {'b','lb','ub','obj','varnames','rownames', ...
              'vartypes','ctypes'};

[m,n] = size(tiger.A);

m_fields = {'b','rownames','ctypes'};
n_fields = {'lb','ub','obj','varnames','vartypes'};

for i = 1 : length(col_fields)
    if size(tiger.(col_fields{i}),1) ~= length(tiger.(col_fields{i}))
        tiger.(col_fields{i}) = tiger.(col_fields{i})';
    end
end

for i = 1 : length(m_fields)
    assert(length(tiger.(m_fields{i})) == m, ...
           'field %s is wrong size\n', m_fields{i});
end

for i = 1 : length(n_fields)
    assert(length(tiger.(n_fields{i})) == n, ...
           'field %s is wrong size\n', n_fields{i});
end
