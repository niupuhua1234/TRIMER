

load('Ecoli_dataset_PROM.mat');

expression = knnimpute(expression);
expression = quantilenorm(expression);

fd=fopen('./Output_data/regulator_PROM.txt','w');
 for j=1:length(regulator)
     if ismember(regulator(j),expressionid)
        fprintf(fd,'%s\n',regulator{j}); 
     end
 end
 fclose(fd);

fd=fopen('./Output_data/target_PROM.txt','w');
 for j=1:length(targets)
     if ismember(regulator(j),expressionid)
        fprintf(fd,'%s\n',targets{j}); 
     end
 end
 fclose(fd);
 

 fd=fopen('./Output_data/expressionname_PROM.txt','w');
 for j=1:length(expressionname)
        fprintf(fd,'%s\n',expressionname{j}); 
 end
 fclose(fd);
 
 fd=fopen('./Output_data/expressionid_PROM.txt','w');
 for j=1:length(expressionid)
        fprintf(fd,'%s\n',expressionid{j}); 
 end
 fclose(fd);
 
fd=fopen('./Output_data/expression_norm_PROM.txt','w');

for i=1:size(expression,1)
    for j=1:size(expression,2)
        fprintf(fd,'%g\t',expression(i ,j)); 
    end
    fprintf(fd,'\n');
end
fclose(fd);


