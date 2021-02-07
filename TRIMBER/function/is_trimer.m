function [tf] = is_trimer(model)
% IS_TRIMER  Returns true if a structure is a TRIMER model.

fields = {'A','b','vartypes','ctypes','rownames','varnames','obj'};
tf = isstruct(model) && all(isfield(model,fields));
