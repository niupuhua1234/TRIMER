function [lbg,ubg,vm]=regulatory_bound(trimer, bnumstobekoed,rxn_affected,rxn_prob,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        REGULATORY_BOUND      bound estimation for each reaction flux in  metabolic model
%        INPUTS:
%        regulatory network   - Format - cell array of regulators and matching target genes
%        rxn_affected   - Reaction affected for each knock out
%        rxn_prob      - Prob associated with reactions affected for each knock out
%        ko_tf         - TF list to be knocked out(default is to  knock all TF in regulatory network one by one )
%
%        PARAMETER:                   
%        thresh         - threshhold for computation stability 
%
%        OUTPUT:
%        lbg,ubg          -  bound estimation  after knock out of all the regulators in the regulatory model;
%        vm               - maximum flux for each reaction 



p = inputParser;
p.addParameter('thresh',10^(-6));
p.parse(varargin{:});
thresh=p.Results.thresh;

if ~iscell(bnumstobekoed)
    bnumstobekoed={bnumstobekoed};
end
if ~iscell(rxn_affected)
    rxn_affected={rxn_affected};
end 
if ~iscell(rxn_prob)
    rxn_prob={rxn_prob};
end 
if ~isfield(trimer,'c')
    grwthpos = find(trimer.obj == 1);
else
    grwthpos = find(trimer.c == 1);
end
[nMets, nRxns] = size(trimer.S);
[rxnpos,genelist] = find(trimer.rxnGeneMat);   % genelist - corresponding genes number.
                                               %rxnpos - metabolic reaction number ,there maybe many reactions under the same gene.

vm = zeros(nRxns,1);                           % initial all Vmax = 0
pos_history=false(nRxns,1);                    % vm estimation history

lbg =map(@(x)trimer.lb(1:nRxns),[1:length(bnumstobekoed)] );
ubg =map(@(x)trimer.ub(1:nRxns),[1:length(bnumstobekoed)] );

sol=fba(trimer); 
if isempty(sol.x)
    warning('Error: no feasible solution ,WT bounds will be returned ')    
    return;
elseif sol.x(grwthpos)<=0
    warning('Error: no feasible solution ,WT bounds will be returned ')    
    return;
else
    v0=sol.x;        %flux extimate under original matrix
end

statbar = statusbar(length(bnumstobekoed),true);
statbar.start('Doing bound estimations');

for ci = 1:length(bnumstobekoed)
%% check if its a metabolic or regulatory gene or both 
    if any(ismember(trimer.genes,bnumstobekoed{ci}))
        temppos = rxnpos(ismember(genelist,find(ismember(trimer.genes,bnumstobekoed{ci}))));
        %temppos =states_check(trimer,{bnumstobekoed{ci}}, temppos ); 
        for jj = 1:length(temppos)
            if trimer.rev(temppos(jj))
                lbg{ci}(temppos) = -thresh;
                ubg{ci}(temppos) = thresh;                
            else% all irrversible reaction has lower bound ==0 , reversible can have lower bound either being == 0 or < 0.
                lbg{ci}(temppos) = -thresh; 
                %ubg{ci}(temppos) = thresh;      
            end
        end
        
    end

    temprxnpos=  rxn_affected{ci};temprxnprob= rxn_prob{ci};    
%% Constrain flux through the reactions associated with these genes        
        for m = 1:length(temprxnpos)                            
                if  ~(temprxnprob(m) < 1),continue;end    % prob=1 means  no change in the original bounds    
                if   (temprxnprob(m)  ~= 0)                % no need for vm estimation,  ust set vm =0 and 1/vm = Inf               
                    if ~pos_history(temprxnpos(m))  % done to save time , if estimated already use it
                        [v11,v12]=T_fva(trimer,'vars',temprxnpos(m),'frac',v0(grwthpos),'status',false);
                        %vm:  maximum flux  vm=0 for v=0
                        if  v0(temprxnpos(m)) < 0     
                            vm(temprxnpos(m)) =  min([v11,v12,v0(temprxnpos(m))]);% make sure vm is -
                        elseif v0(temprxnpos(m)) > 0
                            vm(temprxnpos(m)) =  max([v11,v12,v0(temprxnpos(m))]);% make sure vm is +
                        else
                            vm(temprxnpos(m)) =  max([abs(v11),abs(v12),abs(v0(temprxnpos(m)))]);
                        end
                        pos_history(temprxnpos(m))=true;
                    end
                end

                xx =  vm(temprxnpos(m))*  temprxnprob(m) ;%vm *prob
                if  v0(temprxnpos(m)) < 0                                    
                    tem = max([trimer.lb(temprxnpos(m)),-abs(xx),lbg{ci}(temprxnpos(m))]);  %make sure we arent violating the original bounds trimer.lb ; also get the lowest value if there were multiple modifications for the rxn                    
                    if tem> -thresh,tem=0;end
                    lbg{ci}(temprxnpos(m)) = tem;               % prevents the solver from crashing      
                    %ubg{ci}(temprxnpos(m)) =thresh;               % prevents the solver from crashing     
                end
                if v0(temprxnpos(m)) > 0
                    tem = min([trimer.ub(temprxnpos(m)),abs(xx),ubg{ci}(temprxnpos(m))]);
                    if tem< thresh,tem=0;end
                    ubg{ci}(temprxnpos(m)) = tem;      
                    %lbg{ci}(temprxnpos(m)) =-thresh;       
                end
                                                
        end
           
    statbar.update(ci);
    clear tempgenepos tempgeneprobs temprxnpos k
end
end