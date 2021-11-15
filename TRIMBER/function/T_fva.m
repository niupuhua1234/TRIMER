function [minflux,maxflux] = T_fva(trimer,varargin)
% FVA  Flux Variability Analysis
%
%   [MINFLUX,MAXFLUX] = FVA(trimer,...params...)
%
%   Calculates the minimum and maximum allowable flux through a reaction
%   given a minimum fraction of the objective.
%
%   Parameters
%   'vars'    Indices of variables for the variability calculation.
%             Default is all reactions in the S matrix.
%   'frac'    Fraction of the objective fraction that must be satisfied by
%             each flux distribution.  Default is 1.0.
%   'status'  If true (default), a status bar is displayed.
%   'valtype'  growth constraint type
p = inputParser;
p.addParamValue('vars',1:size(trimer.S,2));
p.addParamValue('frac',1.0);
p.addParamValue('valtype','abs');
p.addParamValue('status',false);
p.parse(varargin{:});

frac = p.Results.frac;
status = p.Results.status;
valtype=p.Results.valtype;

trimer.sense = -1;
vars = convert_ids(trimer.varnames,p.Results.vars,'index');
nvars = length(vars);
minflux = zeros(nvars,1);
maxflux = zeros(nvars,1);
if ~isfield(trimer,'c')
    grwpos=find(trimer.obj);
else
    grwpos=find(trimer.c);
end

if strcmp(valtype,'frac')
    sol = cmpi.solve_mip(trimer);
    frac=sol.x(grwpos);
end

trimer.lb(grwpos)   = frac;
statbar = statusbar(nvars,status);
statbar.start('Flux Variability status');
for i = 1 : nvars
    trimer.obj(:)=0;
    trimer.obj(grwpos)=1; 
    

    trimer.obj(vars(i)) = -1;        
    sol = cmpi.solve_mip(trimer);
    if isempty(sol.x)
        minflux(i) =trimer.lb(vars(i));
    else
        minflux(i) = sol.x(vars(i));
    end
    
    trimer.obj(vars(i)) = 1;    
    sol = cmpi.solve_mip(trimer);
    if isempty(sol.x)
        maxflux(i)=trimer.ub(vars(i));
    else
        maxflux(i) = sol.x(vars(i));
    end
    
    statbar.update(i);    
end

        