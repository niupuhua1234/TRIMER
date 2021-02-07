function [opts] = set_glpk_opts(options)
% SET_GLPK_OPTS  Convert a CMPI options structure into GLPK opts

opts = [];

if isfield(options,'MaxTime')
    opts.tmlim = options.MaxTime;
end
if isfield(options,'MaxIter')
    opts.itlim = options.MaxIter;
end
if isfield(options,'MaxNodes')
    opts.x = options.MaxNodes;
end
if isfield(options,'FeasTol')
    opts.x = options.FeasTol;
end
if isfield(options,'IntFeasTol')
    opts.tolint = options.IntFeasTol;
end
if isfield(options,'OptTol')
    opts.tolobj = options.OptTol;
end
if isfield(options,'AbsOptTol')
    opts.mipgap = options.AbsOptTol;
end

if isfield(options,'Display') && strcmpi(options.Display,'on')
    opts.msglev = 3;  % full output; use 2 for normal output
end
