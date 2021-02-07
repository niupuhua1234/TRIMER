function [minflux,maxflux] = T_fva(trimber,varargin)
% FVA  Flux Variability Analysis
%
%   [MINFLUX,MAXFLUX] = FVA(TRIMBER,...params...)
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
p.addParamValue('vars',1:size(trimber.S,2));
p.addParamValue('frac',1.0);
p.addParamValue('valtype','abs');
p.addParamValue('status',false);
p.parse(varargin{:});

frac = p.Results.frac;
status = p.Results.status;
valtype=p.Results.valtype;

trimber.sense = -1;
vars = convert_ids(trimber.varnames,p.Results.vars,'index');
nvars = length(vars);
minflux = zeros(nvars,1);
maxflux = zeros(nvars,1);
grwpos=find(trimber.c);

if strcmp(valtype,'frac')
    sol = cmpi.solve_mip(trimber);
    frac=sol.x(grwpos);
end

trimber.lb(grwpos)   = frac;
statbar = statusbar(nvars,status);
statbar.start('Flux Variability status');
for i = 1 : nvars
    trimber.obj = trimber.c;    
    
    trimber.obj(vars(i)) = -1;        
    sol = cmpi.solve_mip(trimber);
    minflux(i) = sol.x(vars(i));
    
    trimber.obj(vars(i)) = 1;    
    sol = cmpi.solve_mip(trimber);
    maxflux(i) = sol.x(vars(i));
    
    statbar.update(i);    
end

        