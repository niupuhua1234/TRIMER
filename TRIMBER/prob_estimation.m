function [probtfgene,elm_xn]=prob_estimation(expression,expressionid,regulator,regulated,varargin)
% The PROM algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
%        INPUTS:
%
%        Gene expression data - rows - genes,columns - conditions; no need to normalize or impute
%        Expressionid         - an array of identifiers for each row/gene should be included
%        regulatory network   - format - cell array of regulators and matching target genes
%                               example:   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets ={'GeneA';'GeneB';'GeneC'}
%
%      Parameter:           
%      litevidence          -       The flag vector for  high confidence interactions have the same length as the regulator/target arrays
%                                   high confidence interactions (not necessarily based on literature) should be flagged as one in litevidence 
%                                   array and the rest should be set as 0.
%      prob_prior            -       Prob_prior for flagged interaction. it should be set values between 0 and 1 for those interactions with litevidence.
%                                   ( other values in the array would be ignored)
%      thresh_value         -        threshhold for binarization (default   0.33, (0.2 - 0.4) works for most cases)  
%                                    eg.[0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1]
%      OUTPUT:
%      probtfgene           -     the prob of gene when TF are off

p = inputParser;
p.addParameter('thresh_value',0.33);
p.addParameter('litevidence',zeros(length(regulated),1));
p.addParameter('prob_prior',[]);
p.parse(varargin{:});

prob_prior = p.Results.prob_prior;
litevidence = logical(p.Results.litevidence);
thresh_value = p.Results.thresh_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
lost_xn = false(size(regulated));
elm_xn = false(size(regulated));
probtfgene=zeros(length(regulated),1);
disp('finding probabilities')

data = expression;
if max(data)~=1
    data = knnimpute(data);
    data = quantilenorm(data);
    data1 = data;
    datathresh = quantile(data(:),thresh_value);
    if datathresh < 0
        data(data>=datathresh) = 1;
        data(data < datathresh) = 0;
    else
        data(data < datathresh) = 0;
        data(data>=datathresh) = 1;
    end
else
    data1=data;
end

% u could set those interactions that u think have strong literature evidence to have predefined  probabilities
if ~isempty(litevidence)
    probtfgene(litevidence) = prob_prior(litevidence); 
end                        
for  i = 1:length(regulated)
    k = find(ismember(expressionid,regulated(i))); % [m,n] =ismember  (A,B)   : 判断A中元素是否属于B 。 如果属于 ，m=1 ,  n=在B中的位置.否则 m=0.
    l = find(ismember(expressionid,regulator(i)));
    if ~isempty(k) && ~isempty(l)                 % chech whther gene k and  TF l are included in the microarray data.  (bnumsinexpsn = expressionid)
        
        te = data1(k,:); te1 = data1(l,:);       %te: gene data  ; tel :TF data
        tec = data(k,:); tec1 = data(l,:);
        
        % dependence test for gene data when TF is on and off 
        try kstest2(te(tec1 == 1),te(tec1== 0)); 
            if  (kstest2(te(tec1 == 1),te(tec1== 0)) == 1)
                if ~litevidence(i)
                    probtfgene(i)= sum(tec(tec1 == 0))/length(tec(tec1==0));
                end
                % this formula also gives the same answer  - (sum(~tec1(tec == 1))/length(tec1(tec==1))) * (sum(tec)/sum(~tec1))     
            else
                probtfgene(i) = 1;    %  TF have no effect on gene 认为是无关的 概率为1
                elm_xn(i) = 1;        % if it has strong evidence, consider setting this to zero  
            end   
        % cant be estimated from microarray :te(tec1 == 1)or te(tec1== 0) is empty  ; 
        catch ERRLG   
             probtfgene(i) = 1;       % TF have no effect on gene  无法估计    概率为1 
             lost_xn(i) = 1;          % if it has strong evidence, consider setting this to zero       
             elm_xn(i) = 1;
        end 
    else
        probtfgene(i) = 1;
        elm_xn(i) = 1;
    end
    
end

probtfgene(isnan(probtfgene))=1;
toc

%% check if there is problem with binarizartion
%usually there ll be some interactions that'd be missed, but if most of it (like >75% ) are missed then its time to change the threshold
if (sum(lost_xn) > 0.75*length(probtfgene))  
   datathreshflag = 1;
   disp('change binarization threshold')
else
   datathreshflag = 0;
end
                                   

toc