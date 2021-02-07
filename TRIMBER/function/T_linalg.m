classdef T_linalg
    
properties
    vars
    coefs
    constant
    op
    rhs
end

properties (Dependent)
    terms
    
    is_simple
    is_equation
    is_expression
end

methods (Static)
    function [cs] = get_coefs(terms)
        if ~isa(terms{1},'cell')
            cs = terms(1);
        else
            cs = cellfun(@(x) x{1},terms, 'UniformOutput' ,false);
        end
    end
    
    function [vs] = get_vars(terms)
        if ~isa(terms{1},'cell')
            vs = terms(2);
        else
            vs = map(@(x) x{2},terms);
        end
    end
end

methods
    function [obj] = T_linalg(terms,op,rhs)
        obj.vars = T_linalg.get_vars(terms);
        obj.coefs = T_linalg.get_coefs(terms);
        obj.constant = 0;
        obj.op = op;
        obj.rhs = rhs;
    end
    
    function [ts] = get.terms(obj)
        N = length(obj.vars);
        ts = cell(1,N);
        for i = 1 : N
            ts{i} = {obj.coefs(i),obj.vars{i}};
        end
    end
end

end
            
    