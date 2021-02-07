function [f,v,status] =  PROM_Ecoil(trimber,bnumstobekoed,lb_est,ub_est,rxn_affected,vm,varargin)

% The PROM algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
%      INPUTS:
%      Model                -  metabolic trimber obtained from COBRA toolbox through readcbmodel command
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
%        f                  -growth rate    after knock out of all the regulators in the regulatory trimber;
%        v                  -flux response  after knock out of all the regulators in the regulatory trimber;
%       status              - solver status 
%  
%       EXAMPLE:
%       load mtbpromdata
%       [f_ko,v_ko] = prom(trimber,lb_est,ub_est,rxn_affected,vm);

%===========================================================
%% INPUT HANDLING
%===========================================================
p = inputParser;
p.addParameter('kappa',1);
p.addParameter('growth_pos',find(trimber.c));
p.addParameter('method','sfba');
p.addParameter('delta',0.1);
p.addParameter('epsilon',0.001);

p.parse(varargin{:});
growth_pos= p.Results.growth_pos;
kappa=p.Results.kappa;
method=p.Results.method;
delta=p.Results.delta;
epsilon=p.Results.epsilon;
%===========================================================
%% SOME BASIC INITIALIZATION
%===========================================================
if ~is_trimber(trimber)
    trimber=cobra_to_trimber(trimber); 
end
grwthpos = find(trimber.obj == 1);
[nMets, nRxns] = size(trimber.S);                                            
sol=fba(trimber);   f0=sol.val;v0=sol.x;%flux extimate under original matrix
%trimber=add_growth_constraint(trimber, 1);


%===========================================================
%% trimber setup
%===========================================================
switch method
    case {'fba','FBA'}
        trimberF=trimber;
    case {'sFBA','SFBA','sfba'}
        var_beta=map(@(x) ['Beta_' x],trimber.rxns);
        var_alpha=map(@(x) ['Aplha_' x],trimber.rxns);
        trimberS= add_column(trimber,var_alpha,'c',0,0);%lower bound for v       ,¦Á>0    ,¦Â>0,  upper bound for v       ,¦Á    ,¦Â
        trimberS= add_column(trimberS,var_beta,'c',0,0); 
        lna=T_linalg({{eye(nRxns),trimberS.rxns},{eye(nRxns),var_alpha}},'>',trimber.lb);
        lnb=T_linalg({{eye(nRxns),trimberS.rxns},{-eye(nRxns),var_beta}},'<',trimber.ub);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        trimberS=add_matrix_constraint(trimberS,{lna,lnb},{'Alpha_','Beta_'});
        
        trimberS=change_bound(trimberS,max(trimber.ub),'u',trimberS.rxns);
        trimberS=change_bound(trimberS,min(trimber.lb),'l',trimberS.rxns);
        
               
   case {'room','ROOM'}
        % paramater for ROOM   
        var_bin=map(@(x) ['Room_' x],trimber.rxns);
        trimberR=add_column(trimber, var_bin,'b',0,1);
        
        % Eliminate almost-zero fluxes
        fluxWT=v0(1:nRxns);fluxWT(abs(fluxWT)<1e-6) = 0;
        % generate auxiliary variables
        WT_upperTol = fluxWT + delta*abs(fluxWT) + epsilon;
        WT_lowerTol = fluxWT - delta*abs(fluxWT) - epsilon;

        lnrl=T_linalg({{eye(nRxns),trimber.rxns},{ -diag(trimber.ub -WT_upperTol ) ,var_bin}},'<',WT_upperTol);
        lnru=T_linalg({{eye(nRxns),trimber.rxns},{-diag(trimber.lb-WT_lowerTol ) ,var_bin}},'>',WT_lowerTol);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        lnub=T_linalg({eye(nRxns),trimber.rxns},'<',trimber.ub);
        lnlb=T_linalg({eye(nRxns),trimber.rxns},'>',trimber.lb);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        
        trimberR=add_matrix_constraint(trimberR,{lnrl,lnru,lnub,lnlb},{'upperTol','lowerTol','UB','LB'});
        trimberR=change_obj(trimberR,[zeros(nRxns,1);-ones(nRxns,1)]);
    case {'moma','MOOA'}        
        % Eliminate almost-zero fluxes
        fluxWT=v0(1:nRxns);fluxWT(abs(fluxWT)<1e-6) = 0;
        trimberM=trimber;
end


%===========================================================
%%  FBA  prediction 
%===========================================================

mthresh = 10^(-3);
weights_alpha=zeros(nRxns,1); 
weights_beta=zeros(nRxns,1); 
disp('FBA prediction');
%hw = waitbar(0,'FBA prediction');
statbar = statusbar(length(bnumstobekoed),true);
statbar.start('Doing FBA');
for ci = 1:length(bnumstobekoed)
    
    ub_beta = zeros(nRxns,1); 
    ub_alpha =zeros(nRxns,1);    
    lb_bin =ones(nRxns,1);
    
    temprxnpos=rxn_affected{ci};
    lbg=lb_est{ci}; ubg=ub_est{ci};

    switch method
        case {'sFBA','SFBA','sfba'}

            ub_alpha(temprxnpos(v0(temprxnpos)<0))=max(trimber.ub);
            ub_beta(temprxnpos(v0(temprxnpos)>0))=max(trimber.ub);

            vv=abs(vm(temprxnpos)); vv(vv<mthresh)=mthresh;vv=(kappa*(-1)*abs(f0))./ vv;
            weights_alpha(temprxnpos(v0(temprxnpos)<0)) = vv(v0(temprxnpos)<0);
            weights_beta(temprxnpos(v0(temprxnpos)>0)) = vv(v0(temprxnpos)>0); % new weights based on kappa, normalized with growth rate

            lnpl=T_linalg({{eye(nRxns),trimber.rxns},{eye(nRxns),var_alpha}},'>',lbg);
            lnpu=T_linalg({{eye(nRxns),trimber.rxns},{-eye(nRxns),var_beta}},'<',ubg);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
            trimberS=update_matrix_constraint(trimberS,{lnpl,lnpu},{'Alpha_','Beta_'});
           
            trimberS=change_obj(trimberS, weights_alpha, var_alpha);
            trimberS=change_obj(trimberS, weights_beta, var_beta);
            trimberS=change_bound(trimberS,ub_alpha,'u' ,var_alpha);
            trimberS=change_bound(trimberS,ub_beta, 'u',var_beta);
            sol=fba(trimberS);

       case {'fba','FBA'}
            trimberF=change_bound(trimberF,lbg,'l');
            trimberF=change_bound(trimberF,ubg,'u');
            sol=fba(trimberF);
        case{'moma','MOMA'}
            trimberM=change_bound(trimberM,lbg,'l');
            trimberM=change_bound(trimberM,ubg,'u');
            sol=moma(trimberM,fluxWT);
        case{'room','ROOM'}
            lb_bin(temprxnpos) =1;

            lnrl=T_linalg({{eye(nRxns),trimber.rxns},{ -diag(trimber.ub -WT_upperTol )  ,var_bin}},'<',WT_upperTol);
            lnru=T_linalg({{eye(nRxns),trimber.rxns},{-diag(trimber.lb-WT_lowerTol ) ,var_bin}},'>',WT_lowerTol);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0            
            lnub=T_linalg({eye(nRxns),trimber.rxns},'<',ubg);
            lnlb=T_linalg({eye(nRxns),trimber.rxns},'>',lbg);  % constraint for A*v=0, ¦Á+v>0, v-¦Â<0
        
            trimberR=add_matrix_constraint(trimberR,{lnrl,lnru,lnub,lnlb},{'upperTol','lowerTol','UB','LB'});
            trimberR=change_bound(trimberR,lb_bin,'l', var_bin);    
            sol=   fba(trimberR);
    end
    if  ((sol.flag ~= 2) || any(sol.x(grwthpos) < 0))
        disp(' problem in'); disp(ci);
        v00(ci,1:nRxns)=zeros(nRxns,1);f00(ci)=0;status(ci)=sol.flag
    else
        v00(ci,1:nRxns)=sol.x(1:nRxns); f00(ci)=sol.val;status(ci)=sol.flag;   
    end         % if that doesnt work,  display a warning       
        
    statbar.update(ci);
end

f=v00(:,growth_pos)  ;    v = v00;
end


