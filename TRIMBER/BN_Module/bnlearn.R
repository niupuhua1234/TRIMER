a=getwd() 
setwd('C:/Users/sweet/Desktop/ts/trimer')
library(bnlearn)
library(parallel)
library(dplyr)
library(purrr)
library(R.matlab)
nNodes=detectCores()
cpucl <- makeCluster(nNodes)
stopCluster(cpucl)
source("util.R")
binarization<-function(expression.data,expression.id,Threshold=0.33)
  {
  print('binirization \n')
  datathresh=quantile(as.matrix(expression.data),Threshold)
  expression=data.frame(matrix(data=0, nrow =nrow(expression.data),ncol = ncol(expression.data), 
                               dimnames = list(NULL,expression.id)))
  expression[expression.data>= datathresh]=1
  expression[expression.data<  datathresh]=0
  expression=data.frame(lapply(expression,function (x) factor(x,levels = c('0','1'))))
  return (expression)
  }
data.preprocess<-function(regulators,targets)
  {
  # edge:all directed edge without self loop , tf_tf: edge between TF
  # tf_taget :edge  between TF target 
  names(regulators)=NULL
  names(targets)=NULL
  interaction.id=data.frame('from'=regulators,'to'=targets)
  regulators=unique(regulators)
  targets=unique(targets)
  interaction.edge=interaction.id[-which(interaction.id[,1]==interaction.id[,2]),]
  interaction.tf_tf  =filter(interaction.edge,from %in% regulators & to %in% regulators)
  # filtering the self-loop edge
  interaction=list('raw'=interaction.id,
                   'edge'=  interaction.edge,
                   'tf_tf'=interaction.tf_tf,
                   'tf_target'=dplyr::filter(interaction.edge,from %in% regulators & to %in% targets),
                   'tf_tf_directed'=dplyr::setdiff(interaction.tf_tf,data.frame(from =interaction.tf_tf[,2],to= interaction.tf_tf[,1])),
                   'tf_tf_undirected'=dplyr::intersect(interaction.tf_tf,data.frame(from =interaction.tf_tf[,2],to= interaction.tf_tf[,1])))
  return (interaction)
  }



print('data  loading \n')
iaf1260<- readMat("./Output_data/Ecoli_model_EcoMac.mat")
expression.id=unlist(iaf1260$expressionid)
expression.name=unlist(iaf1260$expressionname)
expression.data =read.table('./Output_data/expression_norm_EcoMac.txt')
expression.data=as.data.frame(t(expression.data))
names(expression.data)=expression.id


model.genes= unlist(iaf1260$model[['genes',1,1]])
model.rxns=  unlist(iaf1260$model[['rxns',1,1]])
model.rxnGeneMat=as.data.frame(t(iaf1260$model[['rxnGeneMat',1,1]]),row.names = model.genes)
colnames(model.rxnGeneMat)=model.rxns
model.rules=iaf1260$model[['rulesforR',1,1]]
model.rules=map(model.rules,~unlist(.x))
names(model.rules)=model.rxns

#########################
print('data arrangement\n')
interaction=data.preprocess(unlist(iaf1260$regulator),unlist(iaf1260$targets))
regulators=unique(unlist(iaf1260$regulator))
targets=unique(unlist(iaf1260$targets))
targets=targets[-which(targets %in% regulators)]

expression=binarization(expression.data,expression.id,0.33)
expression.est=expression[,c(regulators,targets)]
colnames(expression.est)=c(regulators,targets)
#############################
print('structrue learning ')
#blacklisted arcs are never present in the learned graph.
#arcs can only be whitelisted in a single direction== fixed 
bn.learn=bnlearning(expression.est,regulators,targets,interaction$edge,space='specified',algorithm='tabu',rep=1)
bn.est=bn.fit(bn.learn,data = expression.est)

#selection option1
#pa.ch=dsep.search(bn.learn,targets,regulators)
#selection option2
pa.ch=unlist(map2(interaction$raw$from,interaction$raw$to,
           function (x,y){!dsep(bn =bn.learn ,x ,y)}))
pa.ch=interaction$raw[pa.ch,]
#kstest option1 
#interaction.ind=ks_test(expression.data,expression,interaction$raw,alpha=0.05)
#kstest option2 import  the list  return by kstest2 in MATLAB
interaction.ind=read.table('./Output_data/bn_learn_interaction_kstest.csv',sep = ',',col.names=c('from','to') )
pa.ch=intersect(pa.ch,interaction.ind)

#inference option1
#prob.list=appro.query(bn.est,pa.ch)
#interaction.prob=data.frame("from"=pa.ch$from,"to"=pa.ch$to,"prob"=unlist(prob.list))
#filename=sprintf('./Output_data/bn_learn_interaction.csv')
#write.csv(interaction.prob,filename)  # saving
#inference option2
prob.list=Joint.query(bn.est,pa.ch,model.rxnGeneMat,model.rules)
filename=sprintf('./Output_data/bn_learn_interaction_joint.csv')
write.csv(prob.list,filename)  # saving









