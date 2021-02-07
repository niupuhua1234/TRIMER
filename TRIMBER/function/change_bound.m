function model=change_bound(model,value,type,rxnname)
%default is to change the lower and upper bound at the same time 
if nargin <3 || isempty(type), type = 'b'; end
if nargin <4 || isempty(rxnname),rxnname=[1:length(model.lb)]; end        

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
switch lower(type)
    case 'l'
        model.lb(bool)=value;
    case 'u'
        model.ub(bool)=value;
    case 'b'
         model.lb(bool)=value;
         model.ub(bool)=value;
end

end

