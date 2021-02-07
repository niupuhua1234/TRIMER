function model=change_obj(model,value,rxnname)
if nargin <3 || isempty(rxnname),rxnname=[1:length(model.lb)]; end        

if isfield(model,'varnames')
    bool = convert_ids(model.varnames,rxnname,'index');
else
    bool = convert_ids(model.rxns,rxnname,'index');
end
if ~any(bool)
    warning('no reaction found');
end
value(bool==0)='';
bool(bool==0)='';

if isfield(model,'obj')
    model.obj(bool)=value;
else
    model.c(bool)=value;
end
end
