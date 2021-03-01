

load('ecoli_prom_growth_data.mat');
%load('E.colidataset.mat')
model.rxns=replace(model.rxns,'_e','[e]');
model.rxns=replace(model.rxns,'__','-');
model.rxns=replace(model.rxns,{'(',')'},{'[',']'});
trimber=cobra_to_trimer(model);
options.Display='off';
cmpi.set_option(options);
cmpi.set_solver('cplex');

[~,ex_name]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'M13:M162');
lb=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'Q13:Q162');
%[lb_ub,ex_name]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'M13:P162');
bool= ismember(ex_name(:,1),trimber.rxns);ex_name=ex_name(:,1);
trimber=change_bound(trimber,lb,'l',ex_name);
%trimber=change_bound(trimber,lb_ub(:,1),'l',ex_name(:,1));
%trimber=change_bound(trimber,lb_ub(:,2),'u',ex_name(:,1));
source={'EX_succ[e]','EX_nh4[e]','EX_o2[e]' , 'EX_co2[e]', 'EX_pi[e]' ,'EX_so4[e]','EX_h[e]' ,'EX_h2o[e]'};
[~,label,~]=xlsread('asap_data', 'B4:P128');
label=replace(label,{'+','-'},{'1','0'});
label= arrayfun(@(x) str2double(x) ,label);  


[WT_lb_ub,WT_source1,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B14:K15');
[carbon1_lb_ub,carbon_source1,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B18:K69');
[carbon2_lb_ub,carbon_source2,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B72:K88');
[nitrogen_lb_ub,nitrogen_source,~]=xlsread('41586_2004_BFnature02456_MOESM3_ESM.xls', 'Parameters', 'B91:K146');
[~,rules]=xlsread('environ_tf.xlsx', 'MET_NOGENE', 'A1:D216');
index=map(@(x)['x(',num2str(x),')'],1:length(ex_name));
rules(:,3:4)=replace(rules(:,3:4),replace(ex_name,'EX_',''),index');
rec_states={'FBP','TKT2','TALA','PGI','AGDC','PPM2','MNNH','ALTRH','GUI1','MANAO','GUI2','TAGURr','ME2','ME1' ,'GLCptspp','PYK','PFK' ,'LDH_D' , 'LDH_D2', 'SUCCt2_2pp', 'SUCCt2_3pp'};
rec_bool={'true','true','true','true','false','true','true','true','true','true','true',   'true','true','true','true'    ,'true','true','true','true',   'true','true'};                  
cell_states={'"dipyridyl"','"heat shock"','"high NAD"','"high osmolarity"','"LBMedia"','"Rich Medium"','"Oxidative Stress"','"Salicylate"','"Stress"','"Stringent"','pH4','pH7','Growth'};
cell_bool={'true'         ,'true'        ,'true'      ,'true'             ,'true'     ,'false'        ,'true'              ,'true'        ,'true'    ,'false'      ,'true' ,'true' ,'true' };
rules(:,3:4)=replace(rules(:,3:4),cell_states,cell_bool);
rules(:,3:4)=replace(rules(:,3:4),rec_states,rec_bool);
rules1=rules(find(~cellfun(@isempty,rules(:,3))),[1,2,3]);rules2=rules(find(~cellfun(@isempty,rules(:,4))),[1,2,4]);

ko_lb_ub=[WT_lb_ub;carbon1_lb_ub;carbon2_lb_ub;nitrogen_lb_ub];
ko_source= [WT_source1;carbon_source1;carbon_source2;nitrogen_source];
ko_source=replace(ko_source,'(e)','[e]');
ko_source=map(@(x) regexp(x,'; ','split'),ko_source);
num_cond=size(ko_lb_ub,1);

init_bound={};
for i =1:size(ko_lb_ub,1)
    init_bound=[init_bound{:},ko_source{i}];
end
init_bound=map(@(x) ['EX_',x],init_bound);
init_bound= unique({init_bound{:},source{:}});
trimber=change_bound(trimber,zeros(length(init_bound),1),'l',init_bound);

orig_lb=trimber.lb;orig_ub=trimber.ub;


%% PROM running 
%growth_pos=[926,2257,2266,2267];
growth_pos=[1005,2270,2271,2272];
%growth_pos=[1005,2270,2271,2272,1429 2269];
%growth_pos=[1005,887,1438,1439,1440,2271,2272,2269];

ko_tf={'tdcR','crp','malT','glpR','gntR','xylR','asnC','rbsR','ilvY','glnG','rhaS','cpxR','cytR','soxR','melR'};
%ko_tf={'b3119','b3357','b3418','b3423','b3438','b3569','b3743','b3753','b3773','b3868','b3905','b3912','b3934','b4063','b4118'};
[~,loc]=ismember(ko_tf,expressionname);
ko_tf=expressionid(loc,1)';

datathreshvals=0.33;
grRatio=cell(num_cond,1);grRateKO=cell(num_cond,1);grRateWT=cell(num_cond,1);
[~,POS]=ismember(rec_states,trimber.rxns);
trimber.lb(POS)=0;
[probtfgene,elm_xn]=prob_estimation(expression,expressionid,regulator,targets);
ko=ko_source{1};
trimber=change_bound(trimber,-ko_lb_ub(1,:),'l',[{['EX_',ko{1}]},source]); 
sol=fba(trimber); 
 trimber.lb=orig_lb; trimber.ub=orig_ub;
for i =52:127
    fprintf('fba in condition %i \n',i)
    ko=ko_source{i};
    for j=1:length(ko)       
        trimber=change_bound(trimber,-ko_lb_ub(i,:),'l',[{['EX_',ko{j}]},source]);  
    end
    x=true(length(ex_name),1);ex_lb=zeros(length(ex_name),1);
    [bool,pos]=ismember(ex_name,trimber.varnames);
    for j =1:length(pos)
        if pos(j)~=0,ex_lb(j)=trimber.lb(pos(j));else ex_lb(j)=-10; end
    end
    x(find(ex_lb>=0))=false;
    x(~ismember(ex_name,init_bound))=true;
    
    ge_en=zeros(length(rules1(:,3)),1);
    for j =1:length(rules1(:,3))
        ge_en(j)=eval(rules1{j,3});
    end
    gene_en1=rules1(find(~ge_en),1); 
    rules2(:,3)=replace(rules2(:,3),rules1(find(~ge_en),2),'false');
    rules2(:,3)=replace(rules2(:,3),rules1(find(ge_en),2),'true');
    ge_en=zeros(length(rules2(:,3)),1);
    for j =1:length(rules2(:,3))
        ge_en(j)=eval(rules2{j,3});
    end
    gene_en2=rules2(find(~ge_en),1); 
    gene_en=[gene_en1;gene_en2];
    regulator_tmp=[regulator; gene_en(~ismember(gene_en,regulator))]; 
    targets_tmp=[targets; gene_en(~ismember(gene_en,regulator))];
    probtfgene_tmp=[probtfgene;zeros(length(gene_en(~ismember(gene_en,regulator))),1)];
    [lb_est,ub_est,rxn_affected,vmax]  =regulatory_bound(trimber,regulator_tmp,targets_tmp,probtfgene_tmp,'bnumstobekoed',map(@(x) {x,gene_en1{:},gene_en2{:}} ,ko_tf));     
    [grRateKO{i},v,status1]       =ko_prediction(trimber,ko_tf,lb_est,ub_est,rxn_affected,vmax,'growth_pos',growth_pos,'method','sfba');   
        trimber.lb=orig_lb; trimber.ub=orig_ub;
end
