function [sol] = moma(trimer,flux_vals,flux_ids)
% MOMA  Minimization of Metabolic Adjustment
%
%   [SOL] = MOMA(TRIMER,FLUX_VALS)
%   [SOL] = MOMA(TRIMER,FLUX_VALS,FLUX_IDS)
%
%   algorithm:
%       minimize (v - FLUX_VALS)^2
%       s.t.  Sv = 0
%             lb <= v <= ub
%
%   If only a subset of fluxes are specified, the indices can be given in
%   the vector FLUX_IDS.  The corresponding objective is then:
%       minimize (v(FLUX_IDS) - FLUX_VALS)^2
N = size(trimer.A,2);

if nargin < 3
    flux_ids = 1 : size(trimer.S,2);
end

flux_idxs = convert_ids(trimer.varnames,flux_ids,'index');

trimer.Qc.w = zeros(N,1);
trimer.Qc.c = zeros(N,1);

trimer.Qc.w(flux_idxs) = 1;
trimer.Qc.c(flux_idxs) = flux_vals;

trimer.obj(:) = 0;
trimer.qp= true;
sol = cmpi.solve_mip(trimer);
