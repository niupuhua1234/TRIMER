function [sol] = run_solver(mip)

solver = cmpi.get_solver();
if isempty(solver)
    error('No solver selected.  Use set_solver() to choose a solver');
end

switch solver                   
    case 'cplex'        
        opts = cmpi.set_cplex_opts(mip.options);
        
        mip.b = mip.b(:);
        mip.vartypes = upper(mip.vartypes);

        Aineq = [ mip.A(mip.ctypes == '<',:); 
                 -mip.A(mip.ctypes == '>',:)];
        bineq = [ mip.b(mip.ctypes == '<');
                 -mip.b(mip.ctypes == '>')];
             
        Aeq = mip.A(mip.ctypes == '=',:);
        beq = mip.b(mip.ctypes == '=');
        if mip.qp
            [sol.x,sol.val,~,sol.output] = ...
                cplexmiqp(mip.sense.*(2*mip.Q), ...
                          mip.sense.*mip.obj(:), ...
                          Aineq, bineq, ...
                          Aeq, beq, ...
                          [], [], [], ...
                          mip.lb(:), mip.ub(:), ...
                          mip.vartypes(:)', ...
                          [], ...
                          opts);
        else
            [sol.x,sol.val,~,sol.output] = ...
                cplexmilp(mip.sense*mip.obj(:), ...
                          Aineq, bineq, ...
                          Aeq, beq, ...
                          [], [], [], ...
                          mip.lb(:), mip.ub(:), ...
                          mip.vartypes(:)', ...
                          [], ...
                          opts);
        end

        sol.val = mip.sense*sol.val;
                   
        sol.flag = cmpi.get_cplex_flag(sol.output.cplexstatus);
        
    case 'glpk'
        % move equality bounds to A (problem with GLPK?)
        mip = cmpi.bounds_to_constraints(mip,'bounds','equal');
        
        opts = cmpi.set_glpk_opts(mip.options);
        
        mip.b = mip.b(:);
        if mip.qp
            error('GLPK does not support MIQP problems.');
        end
        
        ctypes = mip.ctypes;
        ctypes(ctypes == '<') = 'U';
        ctypes(ctypes == '=') = 'S';
        ctypes(ctypes == '>') = 'L';
        
        if ~isempty(opts)
            [sol.x,sol.val,glpk_flag,sol.output] = ...
                glpk(mip.obj(:),mip.A,mip.b(:),mip.lb(:),mip.ub(:), ...
                     ctypes,upper(mip.vartypes),mip.sense,opts);
        else
            [sol.x,sol.val,glpk_flag,sol.output] = ...
                glpk(mip.obj(:),mip.A,mip.b(:),mip.lb(:),mip.ub(:), ...
                     ctypes,upper(mip.vartypes),mip.sense);
        end
         
        sol.flag = cmpi.get_glpk_flag(glpk_flag);
        
    otherwise
        error('Unrecognized solver: %s',solver);
end
