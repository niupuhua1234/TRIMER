

load('ecoli_prom_growth_data.mat');
load('E.colidataset.mat')
model.rxns=replace(model.rxns,'_e','[e]');
trimer=cobra_to_trimer(model);


%% PROM running 
growth_pos=[1005,2270,2271,2272];
source={'EX_glc[e]','EX_succ[e]','EX_nh4[e]','EX_o2[e]' , 'EX_co2[e]', 'EX_pi[e]' ,'EX_so4[e]','EX_h[e]' ,'EX_h2o[e]'};
ko_tf={'fnr'	'soxS'	'crp'	'lysR'	'fucR'	'malI'	'phoB'	'cpxR'	'trpB'	'trpD' 	'trpE'	'paaX'	'trpA'  'tnaA'	'tnaB'	'dhaR'};
%ko_tf={'appY'   'arcA'   'fnr'  'oxyR'  'soxS'};
%ko_tf={'fnr'	'soxS'	'crp'	'lysR'	'fucR'	'malI'	'phoB'	'cpxR'	'creB'  'trpB'	'trpD' 	'trpE'	'paaX'	'trpA'  'tnaA'	'trpL'  'tnaC'  'tnaB'	'dhaR'}
[~,loc]=ismember(ko_tf,expressionname);
ko_tf=expressionid(loc,1)';

set_solver('glpk');
method='PROM';
datathreshvals=0.33;
if ~strcmp(method,'TRIMBER')
    [prior,interaction]=xlsread('Output_data/bn_learn_interaction_marginal.xlsx');
    targets=interaction(:,2);regulator=interaction(:,1);
    prior    =[prior;zeros(length(ko_tf(~ismember(ko_tf,regulator))),1)];
    targets=[targets;ko_tf(~ismember(ko_tf,regulator))'];
    regulator=[regulator;ko_tf(~ismember(ko_tf,regulator))'];
    
    [probtfgene,elm_xn]              =prob_estimation(expression,expressionid,regulator,targets,'thresh_value',datathreshvals,...
                                                   'litevidence',ones(length(regulator),1),'prob_prior',prior ); 
else
    targets=[targets;ko_tf(~ismember(ko_tf,regulator))'];
    regulator=[regulator;ko_tf(~ismember(ko_tf,regulator))'];
    [probtfgene,elm_xn]              =prob_estimation(expression,expressionid,regulator,targets,'thresh_value',datathreshvals);
end

%trimer=change_bound(trimer, [ -10     0   -10   -10   -15   -15   -10   -10   -55],'l',source);  
[lb_est,ub_est,rxn_affected,vmax]  =regulatory_bound(trimer,regulator,targets,probtfgene,'bnumstobekoed',ko_tf);     
[f,v,status1]       =ko_prediction(trimer,ko_tf,lb_est,ub_est,rxn_affected,vmax,'growth_pos',growth_pos,'method','sfba');   

result=[0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		            0	0	0	0.039333333	0	0.038	                    	0.04	0.040333333];
%0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		0.038333333 0	0	0	0.039333333	0	0.038	0.039333333 0.039666667	0.04	0.040333333
%WT     appY   arcA fnr  arcA\fnr oxyR soxS
%0.71   0.476 0.377 0.41 0.301   0.481 0.465
%0.485 0.636 0.686 0.635 0.648 0.637 0.724

