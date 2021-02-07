function [sol] = fba(trimer)
% FBA  Run Flux Balance Analysis on a TRIMER model.
%   [SOL] = FBA(TRIMER,FLUXNORM,OBJ_FRAC)
%           maximize obj*v
%           s.t. Sv = 0
%                lb <= v <= ub
%
%   Unless specified, OBJ_FRAC defaults to 1.0.
% standard FBA
trimer.sense = -1;
sol = cmpi.solve_mip(trimer);

        


