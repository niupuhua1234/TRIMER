function show_mip(mip,varargin)
% SHOW_MIP  Show equations for a MIP structure
%
%   SHOW_MIP(MIP,...params...)
%
%   Print algebraic equations for a MIP structure.
%
%   Inputs
%   MIP      CMPI model structure
%
%   Parameters
%   'rowidxs'  Indexes of rows (constraints) to show. (default = all rows)
%   'varidxs'  Indexes of variables to show. (default = all rows)
%   'rownames' Cell array of names for each row.  Empty values are 
%              replaced with 'ROW_x'.  May be shorter than the number of 
%              rows; extra 'ROW_x' names are added automatically.
%   'varnames' Cell array of name for each column.  Format is the same as
%              for 'rownames'.
%   'showvars' If true, show bounds on each variable.  (default = true)

[nrows,nvars] = size(mip.A);

default_varnames = arrayfun(@(x) ['x(' num2str(x) ')'], 1:nvars, ...
                            'Uniform', false);
default_rownames = arrayfun(@(x) ['ROW_' num2str(x)], 1:nrows, ...
                            'Uniform', false);

p = inputParser;

p.addParamValue('rowidxs',1 : nrows);
p.addParamValue('varidxs',1 : nvars);
p.addParamValue('rownames',default_rownames);
p.addParamValue('varnames',default_varnames);
p.addParamValue('showvars',true);

p.parse(varargin{:});

showvars = p.Results.showvars;

% default ranges and names
if isfield(mip,'varnames')
    varnames = mip.varnames;
else
    varnames = p.Results.varnames;
end
if isfield(mip,'rownames')
    rownames = mip.rownames;
else
    rownames = p.Results.rownames;
end

varidxs = p.Results.varidxs;
if isa(varidxs,'logical')
    varidxs = find(varidxs);
end

rowidxs = p.Results.rowidxs;
if isa(rowidxs,'logical')
    rowidxs = find(rowidxs);
end                     

% expand the names if an incomplete list was given
varnames = zip_names(default_varnames,varnames);
rownames = zip_names(default_rownames,rownames);

% show optimization direction
if isfield(mip,'sense')
    sense = mip.sense;
else
    sense = 1;
end
fprintf('\n\nOptimization sense:  %s', ...
        iif(sense == 1,'minimize','maximize'));

% show quadratic
qptype = cmpi.miqp_type(mip);
if ~isempty(qptype)
    fprintf('\n\n----- Quadratic Objective -----\n');
    pb = printbuffer();
    pb.start_wrap;
    switch qptype
        case 'Q'
            fprintf('Q:  ');
            [I,J] = find(mip.Q);
            for i = 1 : length(I)
                if I(i) == J(i)
                    rxn_name = [mip.varnames{I(i)} '^2'];
                else
                    rxn_name = [mip.varanems{I(i)} '*' ...
                                mip.varnames{J(i)}];
                end
                pb.printf(make_coef_name_pair(mip.Q(I(i),J(i)), ...
                                              rxn_name, ...
                                              i == 1));
            end
        case 'Qd'
            fprintf('Qd:  ');
            [I,J] = find(mip.Qd);
            for i = 1 : length(I)
                rxn_name = sprintf('(%s - %s)^2',mip.varnames{I(i)}, ...
                                                 mip.varnames{J(i)});
                pb.printf(make_coef_name_pair(mip.Qd(I(i),J(i)), ...
                                              rxn_name, ...
                                              i == 1));
            end
        case 'Qc'
            I = find(mip.Qc.w);
            for i = 1 : length(I)
                rxn_name = sprintf('(%s - %s)^2',mip.varnames{I(i)}, ...
                                                 num2str(mip.Qc.c(I(i))));
                pb.printf(make_coef_name_pair(mip.Qc.w(I(i)), ...
                                              rxn_name, ...
                                              i == 1));
            end
            
    end
    pb.stop_wrap;
end
    
% show objective
fprintf('\n\n----- Objective -----\n');
pb = printbuffer;
pb.start_wrap;
pb.printf('obj:  ');
show_coef_list(mip.obj);
pb.stop_wrap;

% show constraints
fprintf('\n\n----- Constraints -----\n');
for r = 1 : length(rowidxs)
    row = rowidxs(r);
    pb = printbuffer;
    pb.start_wrap;
    pb.printf('%s:  ',rownames{row});
    show_coef_list(mip.A(row,:));
    switch mip.ctypes(row)
        case '<'
            fprintf(' <= ');
        case '>'
            fprintf(' >= ');
        case '='
            fprintf(' = ');
        case 'l'
            fprintf(' < ');
        case 'g'
            fprintf(' > ');
    end
    fprintf('%g\n',mip.b(row));
    pb.stop_wrap;
end

% show indicators
if isfield(mip,'ind') && any(mip.ind)
    fprintf('\n\n----- Indicators -----\n');
    rows = find(mip.ind);
    inds = mip.ind(rows);
    types = mip.indtypes(rows);
    for i = 1 : length(inds)
        if types(i) == 'p'
            fprintf('   %s  => %s\n',rownames{rows(i)},varnames{inds(i)});
        else
            fprintf('   %s <=> %s\n',rownames{rows(i)},varnames{inds(i)});
        end
    end
end

% show bindings
if isfield(mip,'bounds') && any(mip.bounds.var)
    vars = convert_ids(mip.varnames,mip.bounds.var);
    inds = convert_ids(mip.varnames,mip.bounds.ind);
    fprintf('\n\n----- Binding Constraints -----\n');
    for i = 1 : length(vars)
        if mip.bounds.type(i) == 'b'
            fprintf('   %s <= %s <= %s\n', ...
                    make_coef_name_pair(mip.lb(mip.bounds.var(i)), ...
                                        inds{i},true), ...
                    vars{i}, ...
                    make_coef_name_pair(mip.ub(mip.bounds.var(i)), ...
                                        inds{i},true));
        else
            fprintf('   %s <= %s\n', ...
                    vars{i}, ...
                    make_coef_name_pair(mip.ub(mip.bounds.var(i)), ...
                                        inds{i},true));
        end
    end
end

if showvars
    % show variable bounds
    maxlength = max(cellfun(@(x) length(x), varnames));
    fmt = ['%' num2str(maxlength + 2) 's'];
    fprintf('\n\n----- Variable bounds -----\n');
    for i = 1 : length(varidxs)
        fprintf(fmt,[varnames{varidxs(i)} ':']);
        fprintf('  %s',mip.vartypes(varidxs(i)));
        fprintf('  [%g,%g]\n', [mip.lb(varidxs(i)) mip.ub(varidxs(i))]);
    end
end


function [zipped] = zip_names(default,given)
    zipped = default;
    for j = 1 : length(given)
        if ~isempty(given{j})
            zipped{j} = given{j};
        end
    end
end
    
function show_coef_list(constraint)
    nonzeros = find(constraint(varidxs));
    for j = 1 : length(nonzeros)
        col = nonzeros(j);
        coef = constraint(col);
        pb.printf(make_coef_name_pair(coef,varnames{col},j==1));
    end
end

function [str] = make_coef_name_pair(coef,name,first)
    if coef < 0 && first
        sign = '-';
    elseif coef < 0 && ~first
        sign = ' - ';
    elseif first
        sign = '';
    else
        sign = ' + ';
    end
    
    if abs(coef) == 1
        numstr = '';
    else
        numstr = [num2str(abs(coef)) '*'];
    end
    
    str = sprintf('%s%s%s',sign,numstr,name);
end      

end
    