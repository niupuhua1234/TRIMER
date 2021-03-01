function [f,v,status] =  PROM_Ecoil(trimer,lb_est,ub_est,rxn_affected,vm,varargin)

% The  algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
%      INPUTS:
%      Model                -  metabolic trimer obtained from COBRA toolbox through readcbmodel command
%      lb_est\ub_est        -  bound estimation for each reaction flux
%      rxn_affected         -  reaction affected for each gene knock out 
%      vm                   -   maximum of each reaction flux 
%
%      Parameter:                                          
%      growth_pos           -       the flux of reaction of interest(default biomass),
%      kappa                -       the strength of regulation on  constraints 
%                                   eg.[0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,1,5,10];
%      method               -       sFBA, FBA,ROOM ,MOMA 
%      delta,epsilon        -       parameter for ROOM      (default      ,0,1,0,001)
%
%      OUTPUT:
%        f                  -growth rate    after knock out of all the regulators in the regulatory trimer;
%        v                  -flux response  after knock out of all the regulators in the regulatory trimer;
%       status              - solver status 


%===========================================================
%% INPUT HANDLING
%===========================================================
p = inputParser;
p.addParameter('kappa',1);
p.addParameter('growth_pos',find(trimer.c));
p.addParameter('method','sfba');
p.addParameter('delta',0.3);
p.addParameter('epsilon',0.001);

p.parse(varargin{:});
growth_pos= p.Results.growth_pos;
kappa=p.Results.kappa;
method=p.Results.method;
delta=p.Results.delta;
epsilon=p.Results.epsilon;
if ~iscell(bnumstobekoed)
    bnumstobekoed={bnumstobekoed};
end
%===========================================================
%% SOME BASIC INITIALIZATION
%===========================================================
if ~is_trimer(trimer)
    trimer=cobra_to_trimer(trimer); 
end
grwthpos = find(trimer.obj == 1);
[nMets, nRxns] = size(trimer.S);                                            
sol=fba(trimer);   f0=sol.val;v0=sol.x;%flux extimate under original matrix
%trimer=add_growth_constraint(trimer, 1);


%===========================================================
%% trimer setup
%===========================================================
switch method
    case {'fba','FBA'}
        trimerF=trimer;
    case {'sFBA','SFBA','sfba'}
        var_beta=map(@(x) ['Beta_' x],trimer.rxns);
        var_alpha=map(@(x) ['Aplha_' x],trimer.rxns);
        trimerS= add_column(trimer,var_alpha,'c',0,0);%lower bound for v       ,¦Á>0    ,¦Â>0,  upper bound for v       ,¦Á    ,¦Â
        trimerS= add_column(trimerS,var_beta,'c',0,0); 
        lna=T_linalg({{eye(nRxns),trimerS.rxns},{eye(nRxns),var_alpha}},'>',trimer.lb);
        lnb=T_linalg({{eye(nRxns),trimerS.rxns},{-eye(nRxns),var_beta}},'<',trimer.ub);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        trimerS=add_matrix_constraint(trimerS,{lna,lnb},{'Alpha_','Beta_'});
        
        trimerS=change_bound(trimerS,max(trimer.ub),'u',trimerS.rxns);
        trimerS=change_bound(trimerS,min(trimer.lb),'l',trimerS.rxns);
        
               
   case {'room','ROOM'}
        % paramater for ROOM   
        var_bin=map(@(x) ['Room_' x],trimer.rxns);
        trimerR=add_column(trimer, var_bin,'b',0,1);
        
        % Eliminate almost-zero fluxes
        fluxWT=v0(1:nRxns);fluxWT(abs(fluxWT)<epsilon) = 0;
        % generate auxiliary variables
        WT_upperTol = fluxWT + delta*abs(fluxWT) + epsilon;  
        WT_lowerTol = fluxWT - delta*abs(fluxWT) - epsilon; 

        lnrl=T_linalg({{eye(nRxns),trimer.rxns},{ -diag(trimer.ub -WT_upperTol ) ,var_bin}},'<',WT_upperTol);
        lnru=T_linalg({{eye(nRxns),trimer.rxns},{-diag(trimer.lb-WT_lowerTol ) ,var_bin}},'>',WT_lowerTol);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        lnub=T_linalg({eye(nRxns),trimer.rxns},'<',trimer.ub);
        lnlb=T_linalg({eye(nRxns),trimer.rxns},'>',trimer.lb);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        trimerR.options=cmpi.get_options();
        trimerR.options.IntFeasTol=power(10, -round(log10(max(trimer.ub)))-1);
        
        
        trimerR=add_matrix_constraint(trimerR,{lnrl,lnru,lnub,lnlb},{'upperTol','lowerTol','UB','LB'});
        trimerR=change_obj(trimerR,[zeros(nRxns,1);ones(nRxns,1)]);
    case {'moma','MOOA'}        
        % Eliminate almost-zero fluxes
        fluxWT=v0(1:nRxns);fluxWT(abs(fluxWT)<1e-6) = 0;
        trimerM=trimer;
end


%===========================================================
%%  FBA  prediction 
%===========================================================

mthresh = 10^(-3);
weights_alpha=zeros(nRxns,1); 
weights_beta=zeros(nRxns,1); 
disp('FBA prediction');
%hw = waitbar(0,'FBA prediction');
statbar = statusbar(length(lb_est),true);
statbar.start('Doing FBA');
for ci = 1:length(lb_est)
    
    ub_beta = zeros(nRxns,1); 
    ub_alpha =zeros(nRxns,1);    
    lb_bin =ones(nRxns,1);
    
    temprxnpos=rxn_affected{ci};
    lbg=lb_est{ci}; ubg=ub_est{ci};

    switch method
        case {'sFBA','SFBA','sfba'}

            ub_alpha(temprxnpos(v0(temprxnpos)<0))=max(trimer.ub);
            ub_beta(temprxnpos(v0(temprxnpos)>0))=max(trimer.ub);

            vv=abs(vm(temprxnpos)); vv(vv<mthresh)=mthresh;vv=(kappa*(-1)*abs(f0))./ vv;
            weights_alpha(temprxnpos(v0(temprxnpos)<0)) = vv(v0(temprxnpos)<0);
            weights_beta(temprxnpos(v0(temprxnpos)>0)) = vv(v0(temprxnpos)>0); % new weights based on kappa, normalized with growth rate

            lnpl=T_linalg({{eye(nRxns),trimer.rxns},{eye(nRxns),var_alpha}},'>',lbg);
            lnpu=T_linalg({{eye(nRxns),trimer.rxns},{-eye(nRxns),var_beta}},'<',ubg);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
            trimerS=update_matrix_constraint(trimerS,{lnpl,lnpu},{'Alpha_','Beta_'});
           
            trimerS=change_obj(trimerS, weights_alpha, var_alpha);
            trimerS=change_obj(trimerS, weights_beta, var_beta);
            trimerS=change_bound(trimerS,ub_alpha,'u' ,var_alpha);
            trimerS=change_bound(trimerS,ub_beta, 'u',var_beta);
            sol=fba(trimerS);

       case {'fba','FBA'}
            trimerF=change_bound(trimerF,lbg,'l');
            trimerF=change_bound(trimerF,ubg,'u');
            sol=fba(trimerF);
        case{'moma','MOMA'}
            trimerM=change_bound(trimerM,lbg,'l');
            trimerM=change_bound(trimerM,ubg,'u');
            sol=moma(trimerM,fluxWT);
        case{'room','ROOM'}

            lnrl=T_linalg({{eye(nRxns),trimer.rxns},{ -diag(trimer.ub -WT_upperTol )  ,var_bin}},'<',WT_upperTol);
            lnru=T_linalg({{eye(nRxns),trimer.rxns},{-diag(trimer.lb-WT_lowerTol ) ,var_bin}},'>',WT_lowerTol);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0            
            lnub=T_linalg({eye(nRxns),trimer.rxns},'<',ubg);
            lnlb=T_linalg({eye(nRxns),trimer.rxns},'>',lbg);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
            
            trimerR=update_matrix_constraint(trimerR,{lnrl,lnru,lnub,lnlb},{'upperTol','lowerTol','UB','LB'});
            sol=   cmpi.solve_mip(trimerR);
    end
    if  ((sol.flag ~= 2) || any(sol.x(grwthpos) < 0)||isempty(sol.x))
        disp(' problem in'); disp(ci);
        v00(ci,1:nRxns)=zeros(nRxns,1);f00(ci)=0;status(ci)=sol.flag;
    else
        v00(ci,1:nRxns)=sol.x(1:nRxns); f00(ci)=sol.val;status(ci)=sol.flag;   
    end         % if that doesnt work,  display a warning       
        
    statbar.update(ci);
end

f=v00(:,growth_pos)  ;    v = v00;
end


