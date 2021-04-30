#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
a=getwd() 
setwd('C:/Users/sweet/Desktop/ts/trimer')
library(bnlearn)
library(parallel)
library(purrr)


source("util.R")

json_data <-fromJSON(paste(readLines('bn_saved_model.json'), collapse=""))
nNodes=detectCores()
cpucl <- makeCluster(nNodes)
stopCluster(cpucl)

print('data  loading \n')
expression.data =read.table('./Output_data/expression_norm_EcoMac.txt')
expression.id =unlist(read.table('./Output_data/expressionid_EcoMac.txt'))
expression.name =unlist(read.table('./output_data/expressionname_EcoMac.txt'))


regulators =unlist(read.table('./output_data/regulator_EcoMac_coexpression.txt'))
targets =unlist(read.table('./Output_data/target_EcoMac_coexpression.txt'))


all_gene=c(regulators,targets)
all.index=match(all_gene,expression.id)
TFnum=length(regulators)
targetnum=length(targets)

###################
# model generation
##################
true.model=empty.graph(nodes = c(regulators,targets))
#arcs.tf<-arcs(random.graph(regulators,method = "melancon", max.degree = 6))
#new_edge=arcs(true.model)
#for (reg in regulators){
#num_edge=sample(targetnum/2,1)
#new_edge=rbind(new_edge,as.matrix(data.frame(from=rep(reg,num_edge),to=targets[sample(targetnum, num_edge)])))}

new_edge =read.table('bn_saved_model_edge.txt')

arcs(true.model)<-new_edge
edge.list=arcs(true.model)

#########################
random.dist=list()
for (rdv in nodes(true.model)){
  temp.parent=true.model$nodes[[rdv]]$parents
  temp.prob=round(runif(2^length(temp.parent)), 3) #2:number  of decimals
  temp.dist<-rbind(temp.prob,1-temp.prob)
  
  dim(temp.dist)=rep(2,length(temp.parent)+1)
  temp.name=rep(list(c('0','1')),length(temp.parent)+1)
  names(temp.name)=c(rdv,temp.parent)
  dimnames(temp.dist)=temp.name
  print(temp.dist)
  
  random.dist<-c(random.dist,list(temp.dist))}
names(random.dist)=nodes(true.model)
#################################
true.est=custom.fit(true.model,dist = random.dist)
save(random.dist,edge.list,file='./simulated_data/True_model.Rdata')

#############################
print('data generation \n')
size.list=c(100,200,400,800,1600)
data.sim = rbn(true.est,2000)
names(colnames(data.sim ))=NULL
data.list=list()
bn.graph=list()


for (kk in seq(10)){
  ##########################
  for (size in size.list){
    sub_data=data.sim[sample(2000,size),]#
    filename=sprintf('./simulated_data/bn_simulated_data_%d_%d.csv',size,kk)
    write.csv(sub_data,filename)  # saving
    data.list<-c(data.list,list(sub_data))}
  #############################
  # Structure constraint
  #gene can not be parent of regulator
  #gene can not be parent of gene
  ################################
  print('structrue learning ')
  black.list=data.frame('from'= unlist(map(targets,~unname(rep(.x,TFnum)))),
                        'to'  = unlist(rep(regulators,targetnum)))
  black.list=rbind(black.list,set2blacklist(targets))
  
  for (i in seq_along(size.list)){
    bn.learn=tabu(data.list[[i]], blacklist =black.list,tabu = 100,score = 'aic',debug = TRUE)
    bn.graph=c(bn.graph,   bn.learn)    
    diff.edge=dplyr::setdiff(as.data.frame(arcs(bn.learn)),as.data.frame(arcs(true.model)))
    difo.edge=dplyr::setdiff(as.data.frame(arcs(true.model)),as.data.frame(arcs(bn.learn)))
    filename=sprintf('./simulated_data/bn_learn_%d_%d.pdf',size.list[i],kk)
    pdf(filename) 
    bn.graph=graphviz.plot(bn.learn,highlight = list(nodes = unique(c(diff.edge[,1],diff.edge[,2])),
                          arcs =  diff.edge,col = "red",lwd = 1, lty = "solid"),layout = 'fdp')
    dev.off()
  ##################################################
  # infer the conditioning probabilty  for all model
  ##################################################
    bn.est =bn.fit(bn.learn,data = data.list[[i]])
    pa.ch=dsep.search(bn.learn,regulators,targets)
    prob.list=appro.query(bn.est,pa.ch$from,pa.ch$to)
    interaction.prob=data.frame("from"=pa.ch$from,"to"=pa.ch$to,"prob"=unlist(prob.list))
    filename=sprintf('./simulated_data/bn_learn_interaction_%d_%d.csv',size.list[i],kk)
    write.csv(interaction.prob,filename)  # saving
  }
}






  
  

