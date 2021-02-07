function [inf_rows,sol] = find_infeasible_rows(mip,varargin)

p = inputParser;
p.addParamValue('rows',1:size(mip.A,1));
p.addParamValue('obj','count');
p.addParamValue('eps',1e-5);
p.parse(varargin{:});

rows = convert_ids(mip.varnames,p.Results.rows,'index');

assert(ismember(upper(p.Results.obj),{'COUNT','ABS'}), ...
       'parameter ''obj'' must be either ''count'' or ''abs''');
count_obj = strcmpi(p.Results.obj,'count');

nrows = length(rows);

slack_ub = zeros(nrows,1);
slack_lb = zeros(nrows,1);
for i = 1 : nrows
    row_up = mip.A(rows(i),:) .* mip.ub';
    row_dn = mip.A(rows(i),:) .* mip.lb';
    slack_lb(i) = mip.b(rows(i)) - sum(max(row_up,row_dn));
    slack_ub(i) = mip.b(rows(i)) - sum(min(row_up,row_dn));
end

coff = size(mip.A,2);
[mip,slacks] = add_column(mip,[],'c',slack_lb,slack_ub);
for i = 1 : nrows
    mip.A(rows(i),coff+i) = 1;
end

if count_obj
    inds = map(@(x) ['I_',x],slacks);
    mip = add_column(mip,inds,'b');
    mip = bind_var(mip,slacks,inds,'tight',true);
    
    mip.obj(:) = 0;
    mip = set_fieldval(mip,'obj',inds,1);
else
    mip.obj(:) = 0;
    slack_locs = convert_ids(mip,slacks,'index');
    coefs = ones(1,length(slack_locs));
    mip.Q = sparse(slack_locs,slack_locs,coefs, ...
                   size(mip.A,1),size(mip.A,2));
end

mip.sense = 1;
sol = cmpi.solve_mip(mip);
if isempty(sol.x)
    warning('Model could not be made feasible.  Flag: %i\n',sol.flag);
    inf_rows = [];
    return
end

if count_obj
    var_locs = convert_ids(mip.varnames,inds,'index');
else
    var_locs = slack_locs;
end

inf_rows = rows(abs(sol.x(var_locs)) > p.Results.eps);
