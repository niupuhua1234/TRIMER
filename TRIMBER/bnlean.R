a=getwd() 
setwd('C:/Users/sweet/Desktop/ts/csbl-tiger-454756cefb77/prom')
library(bnlearn)
library(parallel)
library(dplyr)
library(purrr)

nNodes=detectCores()
cpucl <- makeCluster(nNodes)
stopCluster(cpucl)

print('data  loading \n')
expression.data =read.table('./Output_data/expression_norm_EcoMac.txt')
expression.id =unlist(read.table('./Output_data/expressionid_EcoMac.txt'))
expression.name =unlist(read.table('./output_data/expressionname_EcoMac.txt'))

regulators =unlist(read.table('./Output_data/regulator_EcoMac.txt'))
targets =unlist(read.table('./Output_data/target_EcoMac.txt'))

#filter the interaction with target=TF
interaction.id =data.frame('from'=regulators,'to'=targets)
interaction.id<-interaction.id[-which(interaction.id[,1]==interaction.id[,2]),]
#interaction.edge<- data.frame(lapply(interaction.id,function (x) expression.id[x]))
interaction.tf_tf<-filter(interaction.id,from %in% regulators & to %in% regulators)
interaction.tf_target<-filter(interaction.edge,from %in% regulators & to %in% targets)
interaction.tf_tf.directed<-setdiff(interaction.tf_tf,data.frame(from =interaction.tf_tf[,2],to= interaction.tf_tf[,1]))
interaction.tf_tf.undirected<-intersect(interaction.tf_tf,data.frame(from =interaction.tf_tf[,2],to= interaction.tf_tf[,1]))
#########################
print('binirization \n')
Threshold=0.33
datathresh=quantile(as.matrix(expression.data),Threshold)
expression=data.frame(matrix(data=0, nrow =ncol(expression.data),ncol = nrow(expression.data), 
                             dimnames = list(NULL,expression.id)))
expression[t(expression.data)>= datathresh]=1
expression[t(expression.data)<  datathresh]=0
print('data arrangement\n')
all_gene=c(regulators,targets)
gene.index=match(all_gene,expression.id)
tf.index=match(regulators,expression.id)
expression.est=data.frame(lapply(expression[,gene.index],
                                 function (x) factor(x,levels = c('0','1'))))
expression.tf =data.frame(lapply(expression[,tf.index],
                                 function (x) factor(x,levels = c('0','1'))))
#############################
# Structure constraint
#gene can not be parent of regulator
#gene can not be parent of gene

TFnum=length(regulators)
targetnum=length(targets)


#black.list=data.frame('from'= unlist(map(targets,~unname(rep(.x,TFnum)))),'to'  = unlist(rep(regulators,targetnum)))
#black.list=rbind(black.list,set2blacklist(targets))


#black.list=data.frame(set2blacklist(regulators))
#black.list=data.frame(set2blacklist(all_gene))
#black.list=setdiff(black.list,interaction.edge)
data.frame("from"=c(interaction.edge$from,interaction.edge$to),"to"=c(interaction.edge$to,interaction.edge$from))


#black.list=setdiff(black.list,interaction.tf_tf.undirected)
#white.list=interaction.tf_tf.directed


print('structrue learning ')
#blacklisted arcs are never present in the learned graph.
#arcs can only be whitelisted in a single direction== fixed 
# constraint learning is not proper for general BN ,the learned edge is only around 6.
#edege  learn by  bic socre  is too small when restrict the seaching space to exprimental 
#verified  interaction list.  without the restric above , the learn network of TFs are too dense 570 edges.

bn.learn=tabu(expression.tf, blacklist =black.list,tabu = 100,score = 'bic',debug = TRUE)
#bn.est= empty.graph(c(regulators,targets))
#arcs(bn.est)=rbind(arcs(bn.learn),as.matrix(interaction.tf_target))
bn.est=bn.fit(bn.est,data = expression.est)

dsep.list=function(model,nodes_list,nodes)
  {Filter(function (x) !dsep(bn =model ,x = x,y = nodes),nodes_list)}



