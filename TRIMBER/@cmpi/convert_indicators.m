function [mip] = convert_indicators(mip)
% CONVERT_INDICATORS  Convert indicators to MILP constraints
%
%   [MIP] = CONVERT_INDICATORS(MIP)
%
%   Converts the indicators indentified by the 'ind' and 'indtypes' fields
%   to MILP constraints.  The constraints generated depend on the 
%   indicator types:
%       'p'  Ax (<=/=/>=) b  => I = 1  becomes  Ax + s (<=/=/>=) b
%                                               I + s >= IND_EPS
%
%       'b'  Ax (<=/=/>=) b <=> I = 1  becomes  Ax + s (<=/=/>=) b
%                                               I + s >= IND_EPS
%                                               s_lb*(1-I) <= s
%                                               s <= s_ub*(1-I)
%
%   The function adds a slack variable s to each constraint that is tied 
%   to an indicator and either 1 or 3 constraints per indicator, depending
%   on if the indicator was a positive ('p') or bound ('b') indicator.
%
%   The CMPI setting IND_EPS is used to construct one of the constraints
%   for each indicator.

IND_EPS = cmpi.get_ind_eps();
% temporary fix
if isempty(IND_EPS)
    IND_EPS = 1e-8;
end

% assert that ind and indtypes fields are present
mip.ind = cmpi.check_field('ind',mip);
mip.indtypes = cmpi.check_field('indtypes',mip);
if isempty(mip.ind)
    mip.ind = zeros(size(mip.b));
end
if isempty(mip.indtypes)
    mip.indtypes = repmat(' ',size(mip.ctypes));
end

% change Ax > b to -Ax < -b and Ax >= b to -Ax <= -b
to_flip = mip.ind > 0 & (mip.ctypes == 'g' | mip.ctypes == '>');
mip.A(to_flip,:) = -mip.A(to_flip,:);
mip.b(to_flip) = -mip.b(to_flip);
mip.ctypes(to_flip & mip.ctypes == 'g') = 'l';
mip.ctypes(to_flip & mip.ctypes == '>') = '<';

rows = find(mip.ind);
inds = mip.ind(rows);
indtypes = mip.indtypes(rows);
[lb,ub] = bound_lhs(mip,rows);

% add extra rows for the <=> indicators
roff = size(mip.A,1);
b_rows = rows(indtypes == 'b');
if ~isempty(b_rows)
    mip = add_row(mip,mip.A(b_rows,:),'<',mip.b(b_rows));
end

for i = 1 : length(rows)
    row = rows(i);
    ind = inds(i);
    indtype = indtypes(i);
    
    if mip.ctypes(row) == '<'
        mip.b(row) = mip.b(row) + IND_EPS;
    end
    
    mip.A(row,ind) = mip.b(row) - lb(i);
    
    if indtype == 'b'
        roff = roff + 1;
        if mip.ctypes(row) == 'l'
            mip.b(roff) = mip.b(roff) - IND_EPS;
        end
        alpha = ub(i) - mip.b(roff);
        mip.A(roff,ind) = alpha;
        mip.b(roff) = mip.b(roff) + alpha;
        mip.ctypes(roff) = '<';
    end
    
    mip.ctypes(row) = '>';
    
    mip.ind(row) = 0;
    mip.indtypes(row) = ' ';
end
mip.ind='';
mip.indtypes='';

