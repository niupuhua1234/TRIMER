function [infeasible] = find_infeasible_constrain(trimer,lb_est,ub_est,rxn_affected,varargin)
% FIND_INFEASIBLE_RULES Find a minimal set of rules which cannot be satisfied when finding a
%   feasible solution.  Reports which rules are not feasible and which 
%   side of the rule is not satisfiable.
%
%   INPUTS
%   lb_est/ub_est  estimated bounds for TF knocks out
%   rxn_affected   reactions affected for each TFs knocks out
%
%   OUPUTS
%   Iinfeasible indices of reactions which are over constrained by
%               estimated bounds
%
%   PARAMETERS
%   'obj_frac'  Define a fraction of the objective value that must be
%               obtained when the rules are added for the trimer to be
%               declared feasible.  Default is 0.01.  It is also possible
%               add the objective constraint to TIGER before calling this
%               function (see ADD_GROWTH_CONSTRAINT).

p = inputParser;
p.addParamValue('display',true);
p.addParamValue('obj_frac',0.01);
p.addParamValue('growth_pos',find(trimer.c));
p.parse(varargin{:});
obj_frac = p.Results.obj_frac;
growth_pos=p.Results.growth_pos;

if ~iscell(rxn_affected)
    rxn_affected={rxn_affected};
end

trimer =check_trimer(trimer);

sol=fba(trimer);
obj_abs=obj_frac*sol.x(growth_pos);

RXN_PRE = 'RXN__';
rxns  = find(~cellfun(@isempty,trimer.grRules));
rxn_names = map(@(x) [RXN_PRE x],trimer.varnames);
trimer =add_column(trimer,rxn_names,'b',0,1);
loc = convert_ids(trimer.varnames,rxn_names,'index');


trimer.options=cmpi.get_options();
trimer.options.IntFeasTol=power(10, -round(log10(max(trimer.ub)))-6);

infeasible=[];

for ci = 1:length(ub_est)
    
    lbg=lb_est{ci}; ubg=ub_est{ci};
     if obj_abs ~= 0
        trimerR = add_growth_constraint(trimer,obj_abs,'pos',growth_pos,'valtype','abs');
        rxn_a=intersect(rxn_affected{ci},rxns);

        rxn_k=1:length(trimer.rxns);rxn_k(rxn_a)='';
        trimerR=change_bound(trimerR, lbg(rxn_k),'l',trimer.rxns(rxn_k));
        trimerR=change_bound(trimerR, ubg(rxn_k),'u',trimer.rxns(rxn_k));

        lnrl=T_linalg({{eye(length(rxn_a)),trimer.rxns(rxn_a)},{ -diag(trimer.ub(rxn_a) -ubg(rxn_a) ) , rxn_names(rxn_a) }},'<',ubg(rxn_a) );
        lnru=T_linalg({{eye(length(rxn_a)),trimer.rxns(rxn_a)},{ -diag(trimer.lb(rxn_a)-lbg(rxn_a)) ,  rxn_names(rxn_a ) }},'>',lbg(rxn_a) );  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        trimerR=add_matrix_constraint(trimerR,{lnrl,lnru},{'upperTol','lowerTol'});
        trimerR.obj(:)=0;
        trimerR.obj(loc(rxn_a))=1;
        sol = cmpi.solve_mip(trimerR);


        if isempty(sol.x)
            warning('Model cannot be made feasible.');
        else
            state = find(sol.x(loc));
            infeasible = [infeasible(:);state];
        end
    end
end

