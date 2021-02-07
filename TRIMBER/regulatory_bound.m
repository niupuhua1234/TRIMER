function [lbg,ubg,rxn_affcted,vm]=regulatory_bound(timber,regulator,regulated,probtfgene,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        REGULATORY_BOUND      bound estimation for each reaction flux in  metabolic model
%        INPUTS:
%        Model                -  metabolic model obtained from COBRA toolbox through readcbmodel command
%        regulatory network   - format - cell array of regulators and matching target genes
%
%        Parameter:           
%        bnumstobekoed          -    TF list to be knocked out(default is to  knock all TF in regulatory network one by one )
%        mode                   - BN : joint probability  CN: conditional  probability
%
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


grwthpos = find(timber.c == 1);
[nMets, nRxns] = size(timber.S);
[~,posgenelist] = ismember(regulated,timber.genes);  %3704 个interaction 中的regulated genes 在model.genes 这个表 中的位置.可能包含0值，即找不到对应。
tfnames=unique(regulator);
[rxnpos,genelist] = find(timber.rxnGeneMat);   % genelist - corresponding genes number.
                                               %rxnpos - metabolic reaction number ,there maybe many reactions under the same gene.

vm = zeros(nRxns,1);                           % initial all Vmax = 0
pos_history=false(nRxns,1);                    % vm estimation history
sol=fba(timber); 
if sol.flag~=2
    warning('error ');
    return;
else
v0=sol.x;        %flux extimate under original matrix
end

lbg =map(@(x)timber.lb,[1:length(bnumstobekoed)] );
ubg =map(@(x)timber.ub,[1:length(bnumstobekoed)] );
rxn_affcted=cell(length(bnumstobekoed),1);


statbar = statusbar(length(bnumstobekoed),true);
statbar.start('Doing bound estimations');

for ci = 1:length(bnumstobekoed)

    % check if its a metabolic or regulatory gene or both
    if any(strcmpi(timber.genes,bnumstobekoed(ci)))
        temppos = rxnpos(genelist == find(strcmp(timber.genes,bnumstobekoed(ci))));
        for jj = 1:length(temppos)
            if timber.rev(temppos(jj))
                lbg{ci}(temppos) = -thresh;
                ubg{ci}(temppos) = thresh;
                
            else
                lbg{ci}(temppos) = -thresh;
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section is for tf-gene relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(ismember(tfnames,bnumstobekoed(ci)))
        
        tfstate = logical(zeros(size(tfnames)));                   %初始化 tfstate
        tfstate( ismember(tfnames,bnumstobekoed(ci))  ) = 1;       %被KO的 TF 状态标记为1

        k = ismember(regulator,tfnames(tfstate))  ;                % 218 个interaction 中，  被KO的TF 参与的interaction 的序号 
        %tempgene = regulated(k);                                  % 这些interaction 中  受控制的基因的名称 ,没有重复.
        tempgeneprobs = probtfgene(k);                             % 这些interaction 中  受该TF 控制的基因的 conditional probility
        
        tempgenepos = posgenelist(k);                              % 这些interaction 中  受调控的基因在 timber.gene 的序号.可能为0 没有对应。
        temprxnpos = rxnpos( ismember(genelist,tempgenepos)  );    % 这些gene 对应 反应的 序号  mapping from gene to  reaction
        
        tempgeneprobs(tempgenepos == 0)  = '';    %delete prob 
        tempgenepos(tempgenepos == 0)  = '';      %delete postion number equal 0 , which means that the tempgene is not found in timber.gene
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section is for gene-protein-reaction relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %change the state of gene 
        x = true(size(timber.genes));
        x(tempgenepos) = false;    %[geneInd,~] =ismember(timber.genes,tempgene) ;  x(geneInd) = false;  
        
        %Figure out if any of the reaction states is changed
        constrainRxn = false(length(temprxnpos),1);
        for j = 1:length(temprxnpos)
            if (~eval(timber.rules{temprxnpos(j)}))
                constrainRxn(j) = true;
            end
        end
        temprxnpos(~constrainRxn)='';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrain flux through the reactions associated with these genes        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temprxnpos has the rxns that are going to  be affected by this tf
% kgenepos  are the indicator for one rxn that will be associated by these target gene
% we loop around all the rxn. 
        rxn_affcted{ci}=   temprxnpos;
        for m = 1:length(temprxnpos)
                %if~constrainRxn(m),continue;end % if any of the reaction states is changed 
                kgenepos = ismember(  tempgenepos ,  genelist(ismember(rxnpos,temprxnpos(m))));%  标记所有受影响的反应中  影响当前单个反应  的所有基因为1.
                kgene=timber.genes(tempgenepos ( kgenepos ));
                x = true(size(timber.genes));
                
                switch method
                    case {'general'}
                        state_list=dec2bin(0:2^sum( kgenepos)-1,sum( kgenepos));
                        state_list= arrayfun(@(x) str2double(x) ,state_list);  
                        state_index=false(length(state_list),1);
                        for i =1:length(state_list)
                             x(tempgenepos(kgenepos))=logical(state_list(i,:));
                             if (eval(timber.rules{temprxnpos(m)}))
                                  state_index(i)=true;
                             end
                             x = true(size(timber.genes));
                        end
                        valid_state=state_list(state_index,:);
                         %rxn_prob =tempgeneprobs   (kgenepos) ;
                    case{'naive'}
                        rxn_prob =min(tempgeneprobs   (kgenepos) );
                end
                                
                if  ~(rxn_prob < 1),continue;end    % prob=1 means  no change in the original bounds    
                if   (rxn_prob  ~= 0)                % no need for vm estimation,  ust set vm =0 and 1/vm = Inf               
                    if ~pos_history(temprxnpos(m))  % done to save time , if estimated already use it
                        [v11,v12]=T_fva(timber,'vars',temprxnpos(m),'frac',v0(grwthpos),'status',false);
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
                    tem = max([timber.lb(temprxnpos(m)),-abs(xx),lbg{ci}(temprxnpos(m))]);  %make sure we arent violating the original bounds timber.lb ; also get the lowest value if there were multiple modifications for the rxn
                    lbg{ci}(temprxnpos(m)) = min([tem,-thresh]);               % prevents the solver from crashing                           
                end
                if v0(temprxnpos(m)) > 0
                    tem = min([timber.ub(temprxnpos(m)),abs(xx),ubg{ci}(temprxnpos(m))]);
                    ubg{ci}(temprxnpos(m)) = max([tem,thresh]);      
                end
                                                
        end
    
       
    end
    statbar.update(ci);
    clear tempgenepos tempgeneprobs temprxnpos k
end
end