function [rxn_affected,rxn_prob]=prob_estimation(trimer,ko_tf,expression,expressionid,regulator,targets,varargin)        
%      PROB_ESTIMATION          Eestimate the conditinonal probability of  gene when the TF is off in the way of PROM
%        INPUTS:
%        expression           - Binarizized gene expression data
%        regulatory network   - format     cell array of regulators and matching target genes
%                               example:   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets ={'GeneA';'GeneB';'GeneC'}
%        ko_tf                -list of TF to be knocked out 
%
%       PARAMETER:           
%       litevidence          -       The flag vector for  high confidence interactions have the same length as the regulator/target arrays
%                                   high confidence interactions (not necessarily based on literature) should be flagged as one in litevidence 
%                                   array and the rest should be set as 0.
%       prob_prior            -       Prob_prior for flagged interaction. it should be set values between 0 and 1 for those interactions with litevidence.
%                                   ( other values in the array would be ignored)
%       OUTPUT:
%       rxn_affected            - reaction affected for each knock out 
%       rxn_prob               - probabilities for temprxnpos  

p = inputParser;
p.addParameter('litevidence',[]);
p.addParameter('prob_prior',[]);
p.parse(varargin{:});

prob_prior = p.Results.prob_prior;
litevidence = p.Results.litevidence;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
probtfgene=zeros(length(targets),1);
disp('finding probabilities')

data = expression;

% u could set those interactions that u think have strong literature evidence to have predefined  probabilities
if ~isempty(litevidence)
    assert(length(litevidence)==length(prob_prior),'the lag vector must have the same length as probability vector');
    if max(litevidence)>1
         probtfgene(litevidence)= prob_prior;
         litevidence_binary=zeros(length(targets),1);
         litevidence_binary(litevidence)=1;
         litevidence=litevidence_binary;
    else
        probtfgene(litevidence) = prob_prior(litevidence); 
    end
else
    litevidence=zeros(length(targets),1);
end                        
for  i = 1:length(targets)
    if ~litevidence(i)
        k = find(ismember(expressionid,targets(i)));
        l = find(ismember(expressionid,regulator(i)));        
        if  ~isempty(k) && ~isempty(l)
            tec = data(k,:); tec1 = data(l,:);           
            probtfgene(i)= sum(tec(tec1 == 0))/length(tec(tec1==0));
        else
            probtfgene=1;        
        end
    end
end

probtfgene(isnan(probtfgene))=1;
[rxn_affected,rxn_prob]=rxn_probvector(trimer,ko_tf,regulator,targets,probtfgene); 
