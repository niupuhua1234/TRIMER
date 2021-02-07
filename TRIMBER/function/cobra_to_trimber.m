function [timber] = cobra_to_timber(cobra)
% COBRA_TO_TIGER  Convert a COBRA model to a TIMBER model
%
%   [TIGER] = COBRA_TO_TIGER(COBRA,...params...)
%
%   Convert a COBRA model structure to a TIMBER model structure.
%

timber =cobra;
[m,n] = size(timber.S);

timber.varnames = cobra.rxns(:);
timber.rownames = cobra.mets(:);
timber.A = timber.S;
% some Cobra models do not have the RHS vector b; add it
if ~isfield(timber,'b')
    timber.b = zeros(size(m,1));
end

timber.obj = cobra.c(:);
timber.ctypes = repmat('=',m,1);
timber.vartypes = repmat('c',n,1);

timber.grRules = cobra.grRules(:);
timber.genes = cobra.genes(:);




    
