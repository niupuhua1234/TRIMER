load('Ecoli_model_EcoMac.mat')
load('EcoMac_data.mat')
model.rxns=replace(model.rxns,'(e)','[e]');
trimer=cobra_to_trimer(model);


options.Display='off';options.MaxTime=100;
%options.Emphasis = 1;
cmpi.set_solver('cplex');
cmpi.set_option(options);

%% PROM running 
growth_pos=[1005,2270,2271,2272];
ko_tf={'fnr'	'soxS'	'crp'	'lysR'	'fucR'	'malI'	'phoB'	'cpxR'	'creB'  'trpB'	'trpD' 	'trpE'	'paaX'	'trpA'  'tnaA'	'trpL'  'tnaC'  'tnaB'	'dhaR'};
ko_tf= map(@(x)replace(x,expressionname,expressionid),ko_tf);
ko_tf_unique=unique(flatten(ko_tf));

source={'EX_glc[e]','EX_o2[e]' };
trimer=change_bound(trimer, [- 9.5,-13.0],'l',source);  
algorithm='TRIME';
method='CN';
datathreshvals=0.33;
[expression,regulator,targets]=data_preprocessing(expression,expressionid,regulator,targets,'thresh_value',datathreshvals); 
R_path=['"C:\Program Files\R\R-4.0.3\bin\Rscript.exe"  '];
Rfun_path=[pwd,'\BN_Module'];
Rmodel_path=[pwd,'\learned_network.bif'];
if  strcmp(algorithm,'TRIMER')
         [rxn_affected_ko,rxn_prob_ko]  =prob_estimation_R(trimer,ko_tf,regulator,targets,R_path,Rmodel_path,'Rfun_path',Rfun_path,'mode','BN');
else
         [rxn_affected_ko,rxn_prob_ko]   =prob_estimation(trimer,ko_tf,expression,expressionid,regulator,targets);   
end
[lb_est,ub_est,vmax]  =regulatory_bound(trimer,ko_tf,rxn_affected_ko,rxn_prob_ko,'thresh',1e-6); 
[f,v,status1]       =ko_prediction(trimer,lb_est,ub_est,rxn_affected_ko,vmax,'growth_pos',growth_pos,'method','ROOM');   

result=[0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		0.038333333 0	0	0	0.039333333	0	0.038	0.039333333 0.039666667	0.04	0.040333333];
