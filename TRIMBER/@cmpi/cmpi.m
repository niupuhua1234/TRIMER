classdef cmpi
% CMPI  Common Mathematical Programming Interface
%
%   CMPI defines a common interface for mathematical programming solvers.
%   Its methods help construct and solve LP/QP/MILP/MIQP problems.
%
%   CMPI represents a mathematical program as a structure with various
%   fields.  A description of these fields can be found in the 
%   documentation for the SOLVE_MIP function (help solve_mip).
    

properties (Constant)
    % Initial (default) values
    init_SOLVER = 'cplex';
    init_IND_EPS = 1e-8;
end

methods (Static)
    function [ctype] = make_ctype(leq,geq,eq)
        % MAKE_CTYPE  Make contraint type identifiers
        %
        %   [CTYPE] = MAKE_CTYPE(LEQ,GEQ,EQ)
        %
        %   Returns a char array defining the constraint types for
        %   LEQ '<=' inequalities, GEQ '<=' inequalities, and 
        %   EQ equalities.
        
        if nargin < 3,   eq = 0; end
        if nargin < 2,  geq = 0; end

        ctype = [repmat('<', 1, leq) ...
                 repmat('>', 1, geq) ...
                 repmat('=', 1, eq)];
    end

    function [vartype] = make_vartype(c,b,i,s,n)
        % MAKE_VARTYPE  Make variable type identifiers
        %
        %   [VARTYPE] = MAKE_VARTYPE(C,B,I,S,N)
        %
        %   Returns a char array defining the variable types for
        %   C continuous, B binary, I integer, S semi-continuous,
        %   and N semi-integer-continuous variables.
        
        if nargin < 5,  n = 0; end
        if nargin < 4,  s = 0; end
        if nargin < 3,  i = 0; end
        if nargin < 2,  b = 0; end

        vartype = [repmat('C', 1, c) ...
                   repmat('B', 1, b) ...
                   repmat('I', 1, i) ...
                   repmat('S', 1, s) ...
                   repmat('N', 1, n)];
    end

    function [sense] = make_sense(s)
        % MAKE_SENSE  Convert optimization sense to an integer
        %
        %   Converts the strings 'min', 'minimize', 'max', and 'maximize'
        %   into the corresponding optimization directions (1 and -1).
        
        switch lower(s)
            case {'min','minimize'}
                sense =  1;
            case {'max','maximize'}
                sense = -1;
        end
    end

    function [tf] = is_acceptable_exit(sol)
        % IS_ACCEPTABLE_EXIT  Determine if a solution was found
        %
        %   Determines if a CMPI solution structure contains an error
        %   code that returns a valid, feasible solution.  Acceptable
        %   exits are:
        %       Optimal
        %       Iteration limit reached
        %       Node limit reached
        %       Time limit reached
        %       Solution limit reached
        
        tf = ( (sol.flag == 2) || (sol.flag == 7) ...
              || (sol.flag == 8) || (sol.flag == 9) ...
              || (sol.flag == 10) ) && ~isempty(sol.x);
    end

    function set_option(option,val)
        % SET_OPTION  Set a default solver option
        %
        %   SET_OPTION(OPTION,VAL) sets an option in the default option 
        %   structure.
        %
        %   SET_OPTION(STRUCT) sets the default option structure to the
        %   structure STRUCT.
        %
        %   For a description of the solver options, see the documentation
        %   for SOLVE_MIP.
        
        global CMPI_OPTIONS
        
        if nargin == 1 && isa(option,'struct')
            CMPI_OPTIONS = option;
        elseif nargin == 1 && isempty(option)
            CMPI_OPTIONS = [];
        else
            CMPI_OPTIONS.(option) = val;
        end
    end

    function clear_options()
        % CLEAR_OPTIONS  Clear the default solver options
        global CMPI_OPTIONS
        CMPI_OPTIONS = [];
    end

    function [opts] = get_options()
        % GET_OPTIONS  Return the default solver options
        cmpi.assert_init();
        global CMPI_OPTIONS
        opts = CMPI_OPTIONS;
    end

    function set_solver(solver)
        % SET_SOLVER  Set the default solver name
        global CMPI_SOLVER
        CMPI_SOLVER = solver;
    end

    function [solver] = get_solver()
        % GET_SOLVER  Get the default solver name
        cmpi.assert_init();
        global CMPI_SOLVER
        solver = CMPI_SOLVER;
    end

    function [ind_eps] = get_ind_eps()
        % GET_IND_EPS  Get epsilon for indicator constraints
        cmpi.assert_init();
        global CMPI_IND_EPS
        ind_eps = CMPI_IND_EPS;
    end
    
    function init()
        % INIT  Initialze CMPI
        %
        %   Creates the global variables and default settings for the
        %   solver name and options.
        
        global CMPI_OPTIONS
        global CMPI_SOLVER
        global CMPI_IND_EPS
        global CMPI_INITIALIZED

        CMPI_SOLVER = cmpi.init_SOLVER;
        CMPI_IND_EPS = cmpi.init_IND_EPS;
        
        CMPI_INITIALIZED = true;
        
        cmpi.set_option('Display','off');
    end

    function assert_init()
        % ASSERT_INIT  Assert that CMPI has been initialized
        global CMPI_INITIALIZED
        
        if ~CMPI_INITIALIZED
            cmpi.init();
        end
    end
    
    % external static methods
    [sol] = solve_mip(mip)
    [tf] = verify_sol(milp,sol,tol)
    show_mip(mip,rowidxs,colidxs,rownames,colnames,showvars)
    [opts,cplex] = set_cplex_opts(options,cplex)
    [flag] = get_cplex_flag(status)
    [flag] = get_glpk_flag(cmpi_flag)
    [opts] = set_glpk_opts(options)
    [mip] = convert_miqp(mip)
    [qtype] = miqp_type(mip)
    [inf_rows,sol] = find_infeasible_rows(mip,varargin)
    [mip] = bounds_to_constraints(mip,varargin)
    [sol] = run_solver(mip)
end

end % class
        
        
    
