function [rxn_affected,rxn_prob]=prob_estimation_R(trimer,ko_tf,regulator,targets,R_path,Rmodel_path,varargin)
%        Interaface function which first preposses the required data and then call R script for BN learning.
%       INPUTS:
%       ko_tf                -list of TF to be knocked out 
%       regulatory network    - format - cell array of regulators and matching target genes
%                                example:   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets ={'GeneA';'GeneB';'GeneC'}
%       R_path                - The path for Rscript.exe in your operation system 
%       Rmodel_path           -  path for saved bayesian network in bif format
%
%       OUTPUT:
%       rxn_affected            - reactions affected for each knock-out TF 
%       rxn_prob               -  probabilities for all reactions affeced by each KO TF  
%
%       PARAMETER:           
%       Rfun_path              -   path for R function (defaut:  current  work directory)
%       mode                  -    CN:conditional vector for regualtory bounds 
%                                 BN:conditional table for regulatory bounds(default)
p = inputParser;
p.addParameter('Rfun_path',[]);
p.addParameter('mode','CN');
p.parse(varargin{:});
Rfun_path = p.Results.Rfun_path;
mode= p.Results.mode;

if isempty(Rfun_path)
    Rfun_path=pwd;
end

rules=replace(trimer.rules,{'x(','1)','2)','3)','4)','5)','6)','7)','8)','9)','0)'},...
{'x[','1]','2]','3]','4]','5]','6]','7]','8]','9]','0]'});
rxnGeneMat=full(trimer.rxnGeneMat);
rxns      = trimer.rxns;
genes      = trimer.genes;
save([Rfun_path,'\info_R.mat'],'rules','rxnGeneMat','rxns','genes','regulator','targets','ko_tf');
%% runing R script for probability inference 
flag=system([R_path,' ',...                      %R path          
             Rfun_path,'\bn_inference.R',' ',... %function path
             Rfun_path,'\info_R.mat',' ',...     % data path
             Rmodel_path,' ',...                 % BN path
             mode,' ',...                        % mode
             Rfun_path]);                        % output path
assert(flag ==0,'R script running fail!')
system(['rm ',Rfun_path,'\info_R.mat']);

[~,~,prob_joint]=xlsread([Rfun_path,'\bn_learn_interaction.csv'] );
system(['rm ',Rfun_path,'\bn_learn_interaction.csv']);
prob_joint=cell2mat(prob_joint(2:end,2:end));
rxn_prob_ko=map(@(x) prob_joint(prob_joint(:,x)~=1,x),1:size(prob_joint,2))';
rxn_affected_ko=map(@(x) find(prob_joint(:,x)~=1),1:size(prob_joint,2))';
