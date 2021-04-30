load('Ecoli_model_EcoMac.mat')
load('EcoMac_data.mat')
model.rxns=replace(model.rxns,'(e)','[e]');
trimer=cobra_to_trimer(model);

options.Display='off';options.MaxTime=60;
cmpi.set_solver('cplex');
cmpi.set_option(options);

%% PROM running 
growth_pos=[1005,2270,2271,2272];
ko_tf={'arcA'  ,  'fnr' ,{ 'arcA' ,  'fnr' } ,'appY' , 'oxyR' , 'soxS'};
%ko_tf={'trpR'  ,  'trpL' , 'tyrR' ,'crp', 'cpxR' , 'soxR','zraR','lrp'};
ko_tf= map(@(x)replace(x,expressionname,expressionid),ko_tf);
ko_tf_unique=unique(flatten(ko_tf));
ko_tf_target=ko_tf_unique(~ismember(ko_tf_unique,regulator));
targets=[targets;ko_tf_target']; regulator=[regulator;ko_tf_target'];

source={'EX_glc[e]','EX_o2[e]' };
%trimer=change_bound(trimer,  [- 8.5,-14.6],'l',source);    
trimer=change_bound(trimer, [ -21.2     0],'l',source);  
algorithm='TRIMER';
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
%cmpi.set_solver('cplex');
[lb_est,ub_est,vmax]  =regulatory_bound(trimer,ko_tf,rxn_affected_ko,rxn_prob_ko,'thresh',1e-6); 
[ff,v,status1]       =ko_prediction(trimer,lb_est,ub_est,rxn_affected_ko,vmax,'growth_pos',growth_pos,'method','sfba');   

result=[0.69,0.63,  0.65,  0.64,0.64,  0.72,0.38,0.41,0.3,0.48,0.48,0.46];
corr([0.7077;f(:,1);0.4914;ff(:,1)],[0.7077;result(1:6)';0.4914;result(7:12)'],'type','Pearson')