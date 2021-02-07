classdef statusbar < handle
    
properties
    N
    width = 80;
    n = 0;
    
    show_N = true;
    show_percent = true;
    show_time = true;
    
    n_before_estimate = 3;
    
    start_tic = [];
    
    margin = '  ';
    barchar = '=';
    
    reprint = false;
    display = true;
    
    update_every = 1;
    last_update = 0;
end

methods
    function [obj] = statusbar(Ntotal,display)
        obj.N = Ntotal;
        
        if nargin == 2
            obj.display = display;
        end
    end
    
    function start(obj,msg)
        if ~obj.display
            return;
        end
        
        obj.start_tic = tic;

        if nargin < 2 || isempty(msg)
            fprintf('\n');
        else
            fprintf('\n%s:\n',msg);
        end
        
        obj.update(0);
    end
    
    function update(obj,n)
        if ~obj.display ...
              || ((n - obj.last_update) < obj.update_every && n > 0)
            return;
        end

        obj.last_update = n;
        
        bar = obj.getbar(n);
        
        if n > 0
            if ~obj.reprint
                fprintf(repmat('\b',1,length(bar)));
            else
                fprintf('\n');
            end
        end
        fprintf('%s',bar);
        
        if n >= obj.N
            fprintf('\n');
        end
    end
    
    function increment(obj,n)
        obj.update(obj.last_update + n)
    end
    
    function [bar] = getbar(obj,n)
        trailer = '';
        frac = n / obj.N;
        
        if obj.show_N
            Nlen = length(sprintf('%i',obj.N));
            trailer = sprintf(' %*i/%*i',Nlen,n,Nlen,obj.N);
        end
        if obj.show_percent
            trailer = sprintf('%s (%5.1f%%)',trailer,frac*100);
        end
        
        if obj.show_time
            timebar = obj.get_timebar(n);
        else
            timebar = '';
        end
        
        barwidth = obj.width - 2*length(obj.margin) - length(trailer) ...
                      - length(timebar) - 2;
        filled = ceil(frac*barwidth);
        if filled == barwidth && n == obj.N
            indbar = repmat(obj.barchar,1,barwidth);
        elseif filled == 0
            indbar = ['>' repmat(' ',1,barwidth - 1)];
        else
            indbar = [repmat(obj.barchar,1,filled-1) '>' ...
                      repmat(' ',1,barwidth - filled)];
        end
        
        if ~isempty(timebar)
            timebar = [timebar ' '];
        end
        bar = [obj.margin timebar '[' indbar ']' trailer obj.margin];
    end
    
    function [bar] = get_timebar(obj,n)
        if n < obj.n_before_estimate
            timestr = '--:--:--';
        else
            elapsed = toc(obj.start_tic);
            estimate = (obj.N - n) * elapsed / n + (n < obj.N);
            timestr = statusbar.make_timestr(estimate);
        end
        
        bar = sprintf('E%s R%s', ...
                      statusbar.make_timestr(toc(obj.start_tic)), ...
                      timestr);
    end
end

methods (Static)
    function test(pause_length)
        if nargin == 0
            pause_length = 0.1;
        end
        
        N = 100;
        s = statusbar(N);
        s.start('Testing statusbar');
        for i = 1 : N
            s.update(i);
            pause(pause_length);
        end
    end
    
    function [timestr] = make_timestr(seconds)
        timestr = datestr(datevec(num2str(seconds),'SS'),'HH:MM:SS');
    end
end

end
        