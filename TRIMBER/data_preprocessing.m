function [data,regulator,target]=data_preprocessing(expression,expressionid,regulator,target,varargin) 
%      DATA_PREPROCESSIN    perform quantile normalization over gene
%      expression data and select interaction that are tested to be depedent by 
%      Kolmogorov-Smirnov test.
%      INPUTS:
%
%      expression           - unnormalized gene expression data: rows -> genes,columns -> conditions
%      expressionid         - an array of identifiers for each row/gene should be included
%      regulatory network   - format - cell array of regulators and matching target genes
%                             example:   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets ={'GeneA';'GeneB';'GeneC'}
%
%      PARAMETER:           
%      thresh_value         -  threshhold for binarization (default   0.33, (0.2 - 0.4) works for most cases)  
%                              eg.[0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1]
%      OUTPUT:
%      regulatory network   - filtered regulatory network

p = inputParser;
p.addParameter('thresh_value',0.33);
p.parse(varargin{:});
thresh_value = p.Results.thresh_value;

if ~all(ismember(target,expressionid)) || ~all(ismember(regulator,expressionid)),
    warning('not enough gene expression for gene in interaction list  ');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = expression;
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

lost_xn= zeros(length(target),1);
elm_xn=zeros(length(target),1);

for  i = 1:length(target)
    k = find(ismember(expressionid,target(i))); % [m,n] =ismember  (A,B)   : 判断A中元素是否属于B 。 如果属于 ，m=1 ,  n=在B中的位置.否则 m=0.
    l = find(ismember(expressionid,regulator(i)));
    if ~isempty(k) && ~isempty(l)
        te = data1(k,:); te1 = data1(l,:);       %te: gene data  ; tel :TF data
        tec = data(k,:); tec1 = data(l,:);
        % dependence test for gene data when TF is on and off 
        try kstest2(te(tec1 == 1),te(tec1== 0)); 
            if  ~(kstest2(te(tec1 == 1),te(tec1== 0)) == 1) 
                elm_xn(i) = 1;        % if it has strong evidence, consider setting this to zero  
            end   
        % cant be estimated from microarray :te(tec1 == 1)or te(tec1== 0) is empty  ; 
        catch ERRLG   
             lost_xn(i) = 1;          % if it has strong evidence, consider setting this to zero       
             elm_xn(i) = 1;
        end     
    else
        lost_xn(i) = 1;          % if it has strong evidence, consider setting this to zero       
        elm_xn(i) = 1;
    end
    
end

%% check if there is problem with binarizartion
%usually there ll be some interactions that'd be missed, but if most of it (like >75% ) are missed then its time to change the threshold
if (sum(lost_xn) > 0.75*length(target))  
   error('change binarization threshold');
end
regulator=regulator(~elm_xn);
target=target(~elm_xn);

