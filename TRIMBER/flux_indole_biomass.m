

%load('Ecoli_dataset_PROM.mat');
load('Ecoli_dataset_EcoMac.mat')
model.rxns=replace(model.rxns,'(e)','[e]');
trimer=cobra_to_trimer(model);


options.Display='off';options.MaxTime=100;
cmpi.set_solver('glpk');
cmpi.set_option(options);

%% PROM running 
growth_pos=[1005,2270,2271,2272];
source={'EX_glc[e]','EX_succ[e]','EX_nh4[e]','EX_o2[e]' , 'EX_co2[e]', 'EX_pi[e]' ,'EX_so4[e]','EX_h[e]' ,'EX_h2o[e]'};
%ko_tf={'fnr'	'soxS'	'crp'	'lysR'	'fucR'	'malI'	'phoB'	'cpxR'	'trpB'	'trpD' 	'trpE'	'paaX'	'trpA'  'tnaA'	'tnaB'	'dhaR'};
%ko_tf={'arcA'  ,  'fnr' ,{ 'arcA' ,  'fnr' } ,'appY' , 'oxyR' , 'soxS'};
ko_tf={'fnr'	'soxS'	'crp'	'lysR'	'fucR'	'malI'	'phoB'	'cpxR'	'creB'  'trpB'	'trpD' 	'trpE'	'paaX'	'trpA'  'tnaA'	'trpL'  'tnaC'  'tnaB'	'dhaR'};
ko_tf= map(@(x)replace(x,expressionname,expressionid),ko_tf);
ko_tf_unique=unique(flatten(ko_tf));
ko_tf_target=ko_tf_unique(~ismember(ko_tf_unique,regulator));


method='TRIMER';
datathreshvals=0.33;
if  strcmp(method,'TRIMER')
    [probtfgene,interaction]=xlsread('Output_data/bn_learn_interaction.csv');
    targets=interaction(:,2);regulator=interaction(:,1);
    targets=[targets;ko_tf_target'];
    regulator=[regulator;ko_tf_target'];
    probtfgene=[probtfgene;zeros(length(ko_tf_target),1)];
else
    targets=[targets;ko_tf_target'];
    regulator=[regulator;ko_tf_target'];
    [probtfgene,elm_xn]              =prob_estimation(expression,expressionid,regulator,targets,'thresh_value',datathreshvals);
end
source={'EX_glc[e]','EX_o2[e]' };
%trimer=change_bound(trimer, [ -15     0   -10   -9   -15   -15   -10   -10   -55],'l',source);  
%trimer=change_bound(trimer, [ -15.0     0],'l',source);  
[lb_est,ub_est,rxn_affected,vmax]  =regulatory_bound(trimer,regulator,targets,probtfgene,'bnumstobekoed',ko_tf,'thresh',1e-6);     
[fff,v,status1]       =ko_prediction(trimer,lb_est,ub_est,rxn_affected,vmax,'growth_pos',growth_pos,'method','sfba');   

%result=[0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		            0	0	0	0.039333333	0	0.038	                    	0.04	0.040333333];
result=[0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		0.038333333 0	0	0	0.039333333	0	0.038	0.039333333 0.039666667	0.04	0.040333333];
%WT     appY   arcA fnr  arcA\fnr oxyR soxS  
%result=[0.636 0.686 0.635 0.648 0.637 0.724];
%result=[0.69,0.63,  0.65,  0.64,0.64,  0.72,0.38,0.41,0.3,0.48,0.48,0.46];
