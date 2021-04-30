load('Ecoli_dataset_EcoMac.mat')
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

if method=='BN'   
        [prob_joint,  regulator] =xlsread('Output_data/bn_learn_interaction_joint.csv');   
        ko_tf_target=ko_tf_unique(~ismember(ko_tf_unique,regulator));
        rxn_prob=map(@(x) prob_joint(prob_joint(:,x)~=1,x),1:length(regulator))';          rxn_affected=map(@(x) find(prob_joint(:,x)~=1),1:length(regulator))';
        [rxn_affected2,rxn_prob2]=rxn_probvector(trimer,ko_tf_target',ko_tf_target',ko_tf_target',zeros(length(ko_tf_target),1)); 
        rxn_affected= [ rxn_affected;rxn_affected2];            rxn_prob=[ rxn_prob; rxn_prob2];       regulator=[regulator,ko_tf_target];    
        [~,kpos]=ismember(ko_tf,regulator);
        rxn_affected_ko=rxn_affected(kpos);                        rxn_prob_ko=rxn_prob(kpos);
elseif method=='CN'
    if  strcmp(algorithm,'TRIMER')
        [probtfgene,interaction]=xlsread('Output_data/bn_learn_interaction.csv'); 
        targets=interaction(:,2);regulator=interaction(:,1);  
    else
        [probtfgene,~]              =prob_estimation(expression,expressionid,regulator,targets,'thresh_value',datathreshvals);
    end    
    ko_tf_target=ko_tf_unique(~ismember(ko_tf_unique,regulator));
    targets=[targets;ko_tf_target'];
    regulator=[regulator;ko_tf_target'];
    probtfgene=[probtfgene;zeros(length(ko_tf_target),1)];
    [rxn_affected_ko,rxn_prob_ko]=rxn_probvector(trimer,ko_tf,regulator,targets,probtfgene); 
end
%cmpi.set_solver('glpk');
[lb_est,ub_est,vmax]  =regulatory_bound(trimer,ko_tf,rxn_affected_ko,rxn_prob_ko,'thresh',1e-6); 
[f,v,status1]       =ko_prediction(trimer,lb_est,ub_est,rxn_affected_ko,vmax,'growth_pos',growth_pos,'method','ROOM');   

result=[0.042666667	0.038666667	0.039666667	0.04	0.039	0.040333333	0.039	0.039333333		0.038333333 0	0	0	0.039333333	0	0.038	0.039333333 0.039666667	0.04	0.040333333];
