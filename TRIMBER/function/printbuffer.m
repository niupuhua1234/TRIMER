classdef printbuffer < handle
    
properties
    wrap = false;
    width = 80;
    curr = 0;
    indent = '   '
end

methods
    function printf(obj,fmt,varargin)
        str = sprintf(fmt,varargin{:});
        l = length(str);
        if obj.wrap
            if l + obj.curr > obj.width
                obj.newline;
                obj.output(obj.indent);
                obj.curr = l + length(obj.indent);
            else
                obj.curr = obj.curr + l;
            end
        end
        obj.output(str);
    end
    
    function start_wrap(obj)
        obj.wrap = true;
        obj.curr = 0;
    end
    
    function stop_wrap(obj)
        obj.wrap = false;
    end
    
    function newline(obj,n)
        if nargin < 2
            n = 1;
        end
        obj.output(repmat('\n',1,n));
    end
    
    function output(obj,str)
        fprintf(str);
    end
    
    function show(obj,str)
        obj.output(str);
        obj.newline();
    end
end

end
    