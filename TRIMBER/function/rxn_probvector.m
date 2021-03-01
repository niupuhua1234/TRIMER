function  [temprxnpos, temprxnprob]=rxn_probvector(trimer,bnumstobekoed,regulator,regulated,probtfgene,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%        RXN_PROBVECTOR             Find affected reaction list and the probabilities vector by each TF
%        INPUTS:
%        Model                -  metabolic model obtained from COBRA toolbox through readcbmodel command
%        regulatory network   - format - cell array of regulators and matching target genes
%
%        Parameter:           
%        bnumstobekoed          -    TF list to be knocked out(default is to  knock all TF in regulatory network one by one )
%
%       OUTPUT:
%       temprxnpos            - reaction affected for each knock out 
%       temprxnprob               - probabilities for temprxnpos  


[~,posgenelist] = ismember(regulated,trimer.genes);  %find the position of each target genes in the gene list,"model.gene"
tfnames=unique(regulator);
[rxnpos,genelist] = find(trimer.rxnGeneMat);         % genelist - corresponding genes number.
                                                     % rxnpos - metabolic reaction number ,there maybe many reactions under the same gene.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section is for tf-gene relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if any(ismember(tfnames,bnumstobekoed))        
        tfstate = logical(zeros(size(tfnames)));                   %初始化 tfstate
        tfstate( ismember(tfnames,bnumstobekoed)  ) = 1;       %被KO的 TF 状态标记为1

        k = ismember(regulator,tfnames(tfstate))  ;                % 218 个interaction 中，  被KO的TF 参与的interaction 的序号 
        %tempgene = regulated(k);                                  % 这些interaction 中  受控制的基因的名称 ,没有重复.
        tempgeneprobs = probtfgene(k);                             % 这些interaction 中  受该TF 控制的基因的 conditional probility
        
        tempgenepos = posgenelist(k);                              % 这些interaction 中  受调控的基因在 trimer.gene 的序号.可能为0 没有对应。
        temprxnpos = rxnpos( ismember(genelist,tempgenepos)  );    % 这些gene 对应 反应的 序号  mapping from gene to  reaction
        
        tempgeneprobs(tempgenepos == 0)  = '';    %delete prob 
        tempgenepos(tempgenepos == 0)  = '';      %delete postion number equal 0 , which means that the tempgene is not found in trimer.gene  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section is for gene-protein relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        temprxnpos=states_check(trimer,tempgenepos,temprxnpos);   % reaction affected with state change for  each TF
        temprxnprob=zeros(1,length(temprxnpos));
        for m = 1:length(temprxnpos)
            kgenepos = ismember(  tempgenepos ,  genelist(ismember(rxnpos,temprxnpos(m))));%  标记所有受影响的反应中  影响当前单个反应  的所有基因为1.
            %kgene=trimer.genes(tempgenepos ( kgenepos ));
            temprxnprob(m) =min(tempgeneprobs(kgenepos) );
        end
else
     temprxnpos=[];
     temprxnprob=[];
end


