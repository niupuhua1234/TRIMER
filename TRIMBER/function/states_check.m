function  temprxnpos=states_check(timber,tempgenepos,temprxnpos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        STATES_CHECK             Find affected reactions with states change for each TF 
%        INPUTS:
%        Model                -  metabolic model obtained from COBRA toolbox through readcbmodel command
%        tempgenepos          - gene affected for each knock out 
%        temprxnpos          -reaction affected for each knock out 
%
%       OUTPUT:
%        temprxnpos         - reaction affected with state change for  each TF
        
%change the state of gene 
x = true(size(timber.genes));
if iscell(tempgenepos) || ischar(tempgenepos)
    [geneInd,~] =ismember(timber.genes,tempgenepos) ; 
    x(geneInd) = false; 
else
    x(tempgenepos) = false;
end

%Figure out if any of the reaction states is changed
constrainRxn = false(length(temprxnpos),1);
for j = 1:length(temprxnpos)
    if (~eval(timber.rules{temprxnpos(j)}))
        constrainRxn(j) = true;
    end
end
temprxnpos(~constrainRxn)='';


