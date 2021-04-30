

load('Ecoli_model_EcoMac.mat')
load('EcoMac_data.mat')

trimer=cobra_to_trimer(model);
options.Display='off';options.MaxTime=1000;
options.IntFeasTol=power(10, -round(log10(max(trimer.ub)))-9);
cmpi.set_solver('cplex');
cmpi.set_option(options);

%tiger=cmpi.convert_indicators(tiger);
[~,ex_name]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'M13:M162');
lb=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'Q13:Q162');
%[lb_ub,ex_name]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'M13:P162');
ex_name=replace(ex_name,'(e)','[e]');
bool= ismember(ex_name(:,1),trimer.rxns);ex_name=ex_name(:,1);
trimer=change_bound(trimer,lb(bool,1),'l',ex_name(bool,1));
%tiger=change_bound(tiger,lb_ub(:,1),'l',ex_name(:,1));
%tiger=change_bound(tiger,lb_ub(:,2),'u',ex_name(:,1));
source={'EX_succ[e]','EX_nh4[e]','EX_o2[e]' , 'EX_co2[e]', 'EX_pi[e]' ,'EX_so4[e]','EX_h[e]' ,'EX_h2o[e]'};
[~,label,~]=xlsread('asap_data', 'B3:P127');
label=replace(label,{'+','-'},{'1','0'});
label= arrayfun(@(x) str2double(x) ,label);  


[WT_lb_ub,WT_source1,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B14:K15');
[carbon1_lb_ub,carbon_source1,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B18:K69');
[carbon2_lb_ub,carbon_source2,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B72:K88');
[nitrogen_lb_ub,nitrogen_source,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B91:K146');



ko_lb_ub=[WT_lb_ub;carbon1_lb_ub;carbon2_lb_ub;nitrogen_lb_ub];
ko_source= [WT_source1;carbon_source1;carbon_source2;nitrogen_source];
ko_source=replace(ko_source,'(e)','[e]');
ko_source=map(@(x) regexp(x,'; ','split'),ko_source);
num_cond=size(ko_lb_ub,1);
all_source=[map(@(x) ['EX_',x],unique(flatten(ko_source))),source];
trimer=change_bound(trimer,zeros(length(all_source),1),'l',all_source);

orig_lb=trimer.lb;orig_ub=trimer.ub;


%% PROM running 
growth_pos=[1005,2270,2271,2272];

ko_tf={'crp'};%{'tdcR','crp','malT','glpR','gntR','xylR','asnC','rbsR','ilvY','glnG','rhaS','cpxR','cytR','soxR','melR'};
ko_tf= map(@(x)replace(x,expressionname,expressionid),ko_tf);
grRateKO=cell(num_cond,1);



algorithm='TRIMER';
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

conditon_Set=[3,4,9,11,28,35,38,39,45,48,52,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,71,72,73,77,78,79,81,84,93,94,95,97,98,101,103,104,105,106,107,108,109,112,113,114,115,118,119,120,121,122,123,124,125,126,127];
grRateKO=zeros(127,15);
for i =3:127
    %i=conditon_Set(k);
    fprintf('fba in condition %i \n',i)
    ko=ko_source{i};
    for j=1:length(ko)       
        trimer=change_bound(trimer,-ko_lb_ub(i,:),'l',[{['EX_',ko{j}]},source]');  
    end
    %tiger.lb(1:2382)=trimer.lb;
    %mip=cmpi.convert_indicators(tiger);
    sol=fba(trimer)
    if ~isempty( sol.x)
        if  sol.x(1005)>0
            [lb_est,ub_est,vmax]  =regulatory_bound(trimer,ko_tf,rxn_affected_ko,rxn_prob_ko,'thresh',1e-6); 
            [f,v,status1]       =ko_prediction(mip,lb_est,ub_est,rxn_affected_ko,vmax,'growth_pos',growth_pos,'method','sfba');     
            grRateKOO(i-2,2)=f(:,1);
        else
            grRateKOO(i-2,2)=zeros(length(ko_tf),1);
        end;
    end
    trimer.lb=orig_lb; trimer.ub=orig_ub;
end

sa=zeros(length(0.01:0.01:0.3),1);j=1
for i=0.01:0.01:0.3
test=grRateKO;
test(test<=i)=0;
test(test> i)=1;
sa(j)=length(find(test == label))/(125*15);j=j+1;
end
out=cell(125,15);

a=grRateKO;
b=grRateKOO;
a(a<=0.15)=0;
a(a>0.15)=1;
b(b<=0.15)=0;
b(b>0.15)=1;
for i=1:125
    for j=1:15
    out{i,j}=[num2str(label(i,j)),'/',num2str(a(i,j)),'/',num2str(b(i,j))];
    end
end
out=replace(out,{'0','1'},{'-','+'});

