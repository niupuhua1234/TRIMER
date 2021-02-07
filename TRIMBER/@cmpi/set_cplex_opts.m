function [opts,cplex] = set_cplex_opts(options,cplex)
% SET_CPLEX_OPTS  Convert a CMPI options structure into CPLEX opts

if nargin > 1
    opts = [];
    
    if isfield(options,'MaxTime')
        cplex.Param.timelimit.Cur = options.MaxTime;
    end
    if isfield(options,'MaxIter')
        cplex.mip.limits.solutions.Cur = options.MaxIter;
    end
    if isfield(options,'MaxNodes')
        cplex.Param.mip.limits.nodes.Cur = options.MaxNodes;
    end
    if isfield(options,'FeasTol')
        cplex.Param.simplex.tolerances.feasibility.Cur = options.FeasTol;
    end
    if isfield(options,'IntFeasTol')
        cplex.Param.mip.tolerances.integrality.Cur = options.IntFeasTol;
    end
    if isfield(options,'OptTol')
        cplex.Param.mip.tolerances.mipgap.Cur = options.OptTol;
    end
    if isfield(options,'AbsOptTol')
        cplex.Param.mip.tolerances.absmipgap.Cur = options.AbsOptTol;
    end
else
    cplex = [];
    
    opts = cplexoptimset;
    setifdef('MaxTime','MaxTime');
    setifdef('MaxIter','MaxIter');
    setifdef('MaxNodes','MaxNodes');
    setifdef('FeasTol','EpRHS');
    setifdef('IntFeasTol','TolXInteger');
    setifdef('OptTol','EpGap');
    setifdef('Display','Display');
    if isfield(options,'Display') && strcmpi(options.Display,'on')
        opts.Display = 'iter';
    end
end
    
function setifdef(cmpifield,cplexfield)
    if isfield(options,cmpifield)
        opts.(cplexfield) = options.(cmpifield);
    end
end

end