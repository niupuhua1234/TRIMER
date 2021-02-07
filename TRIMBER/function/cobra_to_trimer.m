function [trimer] = cobra_to_trimer(cobra)
% COBRA_TO_TIGER  Convert a COBRA model to a TRIMER model
%
%   [TIGER] = COBRA_TO_TIGER(COBRA,...params...)
%
%   Convert a COBRA model structure to a TRIMER model structure.
%

trimer =cobra;
[m,n] = size(trimer.S);

trimer.varnames = cobra.rxns(:);
trimer.rownames = cobra.mets(:);
trimer.A = trimer.S;
% some Cobra models do not have the RHS vector b; add it
if ~isfield(trimer,'b')
    trimer.b = zeros(size(m,1));
end

trimer.obj = cobra.c(:);
trimer.ctypes = repmat('=',m,1);
trimer.vartypes = repmat('c',n,1);

trimer.grRules = cobra.grRules(:);
trimer.genes = cobra.genes(:);




    
