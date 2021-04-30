function  [rxn_affected,rxn_prob]=rxn_probvector(trimer,bnumstobekoed,regulator,regulated,probtfgene)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%        RXN_PROBVECTOR             Find affected reaction list and the probabilities vector by each TF
%        INPUTS:
%       bnumstobekoed         -tf to be knocked out 
%        regulatory network   - format - cell array of regulators and matching target genes
%
%        Parameter:           
%        bnumstobekoed          -    TF list to be knocked out(default is to  knock all TF in regulatory network one by one )
%
%       OUTPUT:
%       ko_rxn.rxn_affected            - reaction affected for each knock out 
%       ko_rxn.rxn_prob               - probabilities for temprxnpos  
if ~iscell(bnumstobekoed)
    bnumstobekoed={bnumstobekoed};
end

[~,posgenelist] = ismember(regulated,trimer.genes);  %find the position of each target genes in the gene list,"model.gene"
tfnames=unique(regulator);
[rxnpos,genelist] = find(trimer.rxnGeneMat);         % genelist - corresponding genes number.
                                                     % rxnpos - metabolic reaction number ,there maybe many reactions under the same gene.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section is for tf-gene relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
rxn_affected=cell(length(bnumstobekoed),1);
rxn_prob=cell(length(bnumstobekoed),1);
for ci = 1:length(bnumstobekoed)
    if any(ismember(tfnames,bnumstobekoed{ci}))        
            tfstate = logical(zeros(size(tfnames)));                
            tfstate( ismember(tfnames,bnumstobekoed{ci})  ) = 1;      

            k = ismember(regulator,tfnames(tfstate))  ;                % find all the interaction of the KO-TFs 
            %tempgene = regulated(k);                                  % %find the targets gene of these interaction pair
            tempgeneprobs = probtfgene(k);                             % conditional probility for these interation pairs

            tempgenepos = posgenelist(k);                                   % find the index of these target genes in trimer.gene.There may not be match for some target genes
            temprxnpos = unique(rxnpos( ismember(genelist,tempgenepos)  )); % find the rxns index for these genes. There maybe gene which is controlling by multiple reactions.

            tempgeneprobs(tempgenepos == 0)  = '';    %delete prob 
            tempgenepos(tempgenepos == 0)  = '';      %delete postion number equal 0 , which means that the tempgene is not found in trimer.gene  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this section is for gene-protein relationship
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            temprxnpos=states_check(trimer,tempgenepos,temprxnpos);   % reaction affected with state change for  each TF
            temprxnprob=zeros(1,length(temprxnpos));
            for m = 1:length(temprxnpos)
                kgenepos = ismember(  tempgenepos ,  genelist(ismember(rxnpos,temprxnpos(m))));%  find all the genes which control the current reaction. 
                %kgene=trimer.genes(tempgenepos ( kgenepos ));
                temprxnprob(m) =min(tempgeneprobs(kgenepos) );
            end
            rxn_affected{ci}=temprxnpos;
            rxn_prob{ci}=temprxnprob;
    else
         rxn_affected{ci}=[];
         rxn_prob{ci}=[];
    end
end



