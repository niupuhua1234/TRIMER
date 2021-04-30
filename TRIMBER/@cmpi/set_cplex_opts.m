function [opts] = set_cplex_opts(options)
% SET_CPLEX_OPTS  Convert a CMPI options structure into CPLEX opts

opts = cplexoptimset('cplex');

if isfield(options,'MaxTime')
    opts.timelimit = options.MaxTime;
end
if isfield(options,'MaxIter')
    opts.mip.limits.solutions = options.MaxIter;
end
if isfield(options,'MaxNodes')
    opts.mip.limits.nodes = options.MaxNodes;
end
if isfield(options,'FeasTol')
    opts.simplex.tolerances.feasibility = options.FeasTol;
end
if isfield(options,'IntFeasTol')
    opts.mip.tolerances.integrality = options.IntFeasTol;
end
if isfield(options,'OptTol')
    opts.mip.tolerances.mipgap = options.OptTol;
end
if isfield(options,'AbsOptTol')
    opts.mip.tolerances.absmipgap= options.AbsOptTol;
end
if isfield(options,'Emphasis')
    opts.emphasis.mip=options.Emphasis;
end
if isfield(options,'Display') && strcmpi(options.Display,'on')
    opts.display = 'iter';
end
    
