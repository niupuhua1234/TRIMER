function [sol] = solve_mip(mip)
% SOLVE_MIP  Solve a Mixed Integer Programming problem
%
%   [SOL] = SOLVE_MIP(MIP,SOLVER)
%
%   Solves the Mixed Integer Linear Programming problem
%       min sense*obj*x
%       subject to
%           A*x (<=/=/>=) b
%           lb <= x <= ub
%
%   or the Mixed Integer Quadratic Programming problem
%       min sense*(x'*Q*x + obj*x)
%       subject to
%           A*x (<=/=/>=) b
%           lb <= x <= ub
%
%   Inputs
%   MIP     Problem structure.  Fields include:
%               obj       Objective coefficients
%               sense     Direction of optimization:
%                              1    Minimization
%                             -1    Maximization
%               A         Coefficient matrix for constraints.
%               b         Right-hand side of constraints.
%               ctypes    Char array denoting the type of each constraint:
%                             '<'   a*x <= b
%                             '='   a*x == b
%                             '>'   a*x >= b
%               lb        Lower bound for each variable.
%               ub        Upper bound for each variable.
%               vartypes  Char array denoting the type of each variable:
%                             'C'   Continuous
%                             'B'   Binary
%                             'I'   Integer
%               options   Solver options (see below).  If not given, the
%                         options in the global variable CMPI_OPTIONS are
%                         used (if defined).
%                             MaxTime     Maximum solution time (seconds)
%                             MaxIter     Maximum simplex iterations
%                             MaxNodes    Maximum number of nodes
%                             Display     Turn reporting to the screen
%                                         'on' or 'off' (default)
%                             FeasTol     Feasibility tolerance
%                             IntFeasTol  Integer feasibility tolerance
%                             OptTol      Optimality tolerance
%                             AbsOptTol   Absolute optimality tolerance
%               Q         Quadratic objective matrix.  If given, the
%                         problem is solved as a MIQP.  See CONVERT_MIQP
%                         for details on Q and related fields.
%               prepared  True if MIP has been preprocessed by PREPARE_MIP.
%
%   Outputs
%   SOL     Solution structure with fields:
%               x       Optimal vector
%               val     Objective values (sense*c*x)
%               flag    Exit flag.  Possible values are:
%                           1   Not started
%                           2   Optimal
%                           3   Infeasible
%                           4   Infeasible or unbounded
%                           5   Unbounded
%                           6   Objective worse than user cutoff
%                           7   Iteration limit reached
%                           8   Node limit reached
%                           9   Time limit reached
%                           10  Solution limit reached
%                           11  User interruption
%                           12  Numerical difficulties
%                           13  Suboptimal solution
%               output  Other solver-specific output

if ~issparse(mip.A)
    mip.A = sparse(mip.A);
end

if ~isfield(mip,'sense') || isempty(mip.sense)
    mip.sense = 1;
end
if ~isfield(mip,'qp') || isempty(mip.qp)
    mip.qp = 0;
end
if ~isfield(mip,'options') || isempty(mip.options)
    mip.options = cmpi.get_options();
end
if mip.qp
    mip = cmpi.convert_miqp(mip);
end

sol = cmpi.run_solver(mip);
