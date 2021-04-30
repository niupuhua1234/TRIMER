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
p.addParameter('litevidence',zeros(length(targets),1));
p.addParameter('prob_prior',[]);
p.parse(varargin{:});

prob_prior = p.Results.prob_prior;
litevidence = logical(p.Results.litevidence);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lost_xn = false(size(targets));
elm_xn = false(size(targets));
probtfgene=zeros(length(targets),1);
disp('finding probabilities')

data = expression;

% u could set those interactions that u think have strong literature evidence to have predefined  probabilities
if ~isempty(litevidence)
    probtfgene(litevidence) = prob_prior(litevidence); 
end                        
for  i = 1:length(targets)
    k = ismember(expressionid,targets(i)); % [m,n] =ismember  (A,B)   : 判断A中元素是否属于B 。 如果属于 ，m=1 ,  n=在B中的位置.否则 m=0.
    l = ismember(expressionid,regulator(i));        
    tec = data(k,:); tec1 = data(l,:);       
    if ~litevidence(i)
        probtfgene(i)= sum(tec(tec1 == 0))/length(tec(tec1==0));
    end
    % this formula also gives the same answer  - (sum(~tec1(tec == 1))/length(tec1(tec==1))) * (sum(tec)/sum(~tec1))     
end

probtfgene(isnan(probtfgene))=1;
[rxn_affected,rxn_prob]=rxn_probvector(trimer,ko_tf,regulator,targets,probtfgene); 
