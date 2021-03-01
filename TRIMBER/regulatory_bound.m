function [lbg,ubg,rxn_affcted,vm]=regulatory_bound(trimer,regulator,regulated,probtfgene,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        REGULATORY_BOUND      bound estimation for each reaction flux in  metabolic model
%        INPUTS:
%        Model                -  metabolic model obtained from COBRA toolbox through readcbmodel command
%        regulatory network   - format - cell array of regulators and matching target genes
%
%        Parameter:           
%        bnumstobekoed          -    TF list to be knocked out(default is to  knock all TF in regulatory network one by one )
%        method                   - BN : joint probability  CN: conditional  probability
%        thresh                  - threshhold for computation stability 
%       OUTPUT:
%       lbg,ubg                -  bound estimation  after knock out of all the regulators in the regulatory model;
%       rxn_affcted             - reaction affected for each knock out 
%              vm               - maximum flux for each reaction 
%  


p = inputParser;
p.addParameter('bnumstobekoed',unique(regulator));
p.addParameter('method','naive');
p.addParameter('thresh',10^(-6));

p.parse(varargin{:});
bnumstobekoed= p.Results.bnumstobekoed;
method=p.Results.method;
thresh=p.Results.thresh;
if ~iscell(bnumstobekoed)
    bnumstobekoed={bnumstobekoed};
end
   
grwthpos = find(trimer.c == 1);
[nMets, nRxns] = size(trimer.S);
tfnames=unique(regulator);
[rxnpos,genelist] = find(trimer.rxnGeneMat);   % genelist - corresponding genes number.
                                               %rxnpos - metabolic reaction number ,there maybe many reactions under the same gene.

vm = zeros(nRxns,1);                           % initial all Vmax = 0
pos_history=false(nRxns,1);                    % vm estimation history

lbg =map(@(x)trimer.lb,[1:length(bnumstobekoed)] );
ubg =map(@(x)trimer.ub,[1:length(bnumstobekoed)] );
rxn_affcted=cell(length(bnumstobekoed),1);


sol=fba(trimer); 
if sol.flag~=2
    warning('Error: no feasible solution ,WT bounds will be returned ')    
    return;
else
    v0=sol.x;        %flux extimate under original matrix
end

statbar = statusbar(length(bnumstobekoed),true);
statbar.start('Doing bound estimations');

for ci = 1:length(bnumstobekoed)

    % check if its a metabolic or regulatory gene or both
    if any(ismember(trimer.genes,bnumstobekoed{ci}))
        temppos = rxnpos(ismember(genelist,find(ismember(trimer.genes,bnumstobekoed{ci}))));
        temppos =states_check(trimer,{bnumstobekoed{ci}}, temppos ); 
        for jj = 1:length(temppos)
            if trimer.rev(temppos(jj))
                lbg{ci}(temppos) = -thresh;
                ubg{ci}(temppos) = thresh;                
            else
                lbg{ci}(temppos) = -thresh;
            end
        end
        
    end
    if any(ismember(tfnames,bnumstobekoed{ci}))
        [temprxnpos,temprxnprob]=rxn_probvector(trimer,bnumstobekoed{ci},regulator,regulated,probtfgene);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrain flux through the reactions associated with these genes        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rxn_affcted{ci}=   temprxnpos;
        for m = 1:length(temprxnpos)
                rxn_prob =temprxnprob(m);
                switch method
                    case {'general'}
                        state_list=dec2bin(0:2^sum( kgenepos)-1,sum( kgenepos));
                        state_list= arrayfun(@(x) str2double(x) ,state_list);  
                        state_index=false(length(state_list),1);
                        for i =1:length(state_list)
                             x(tempgenepos(kgenepos))=logical(state_list(i,:));
                             if (eval(trimer.rules{temprxnpos(m)}))
                                  state_index(i)=true;
                             end
                             x = true(size(trimer.genes));
                        end
                        valid_state=state_list(state_index,:);
                         %rxn_prob =tempgeneprobs   (kgenepos) ;
                end
                                
                if  ~(rxn_prob < 1),continue;end    % prob=1 means  no change in the original bounds    
                if   (rxn_prob  ~= 0)                % no need for vm estimation,  ust set vm =0 and 1/vm = Inf               
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

                xx =  vm(temprxnpos(m))*  rxn_prob ;%vm *prob
                if  v0(temprxnpos(m)) < 0                                    
                    tem = max([trimer.lb(temprxnpos(m)),-abs(xx),lbg{ci}(temprxnpos(m))]);  %make sure we arent violating the original bounds trimer.lb ; also get the lowest value if there were multiple modifications for the rxn
                    lbg{ci}(temprxnpos(m)) = min([tem,-thresh]);               % prevents the solver from crashing                           
                end
                if v0(temprxnpos(m)) > 0
                    tem = min([trimer.ub(temprxnpos(m)),abs(xx),ubg{ci}(temprxnpos(m))]);
                    ubg{ci}(temprxnpos(m)) = max([tem,thresh]);      
                end
                                                
        end
    
       
    end
    statbar.update(ci);
    clear tempgenepos tempgeneprobs temprxnpos k
end
end