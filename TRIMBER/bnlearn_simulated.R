#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
a=getwd() 
setwd('D:/я╦ювобть/PROM_Ecoil')
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

#regulators =unlist(read.table('./Output_data/regulator_EcoMac.txt'))
#targets =unlist(read.table('./Output_data/target_EcoMac.txt'))
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
arcs.tf<-arcs(random.graph(regulators,method = "melancon", max.degree = 6))

new_edge=arcs(true.model)
for (reg in regulators){
num_edge=sample(targetnum/2,1)
new_edge=rbind(new_edge,
               as.matrix(data.frame(from=rep(reg,num_edge),to=targets[sample(targetnum, num_edge)])))
}
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
#data=matrix(sample(2,10000*50,replace =TRUE)-1,nrow=10000,ncol=50,dimnames=list(NULL,c(regulators,targets)))
#data =data.frame(map(as.data.frame(data),~factor(.x,levels = c('0','1'))))
#############################
print('data generation \n')
size.list=c(100,200,400,800,1600)
data.sim = rbn(true.est,2000)
names(colnames(data.sim ))=NULL
data.list=list()

for (kk in range(0,10,1)){
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
  
  
  for (i in seq(5)){
    bn.learn=tabu(data.list[[i]], blacklist =black.list,tabu = 100,score = 'bic',debug = TRUE)
    diff.edge=setdiff(as.data.frame(arcs(bn.learn)),as.data.frame(arcs(true.model)))
    
  
    filename=sprintf('./simulated_data/bn_learn_%d_%d.pdf',size.list[i],kk)
    pdf(filename) 
    graphviz.plot(bn.learn,highlight = list(nodes = unique(c(diff.edge[,1],diff.edge[,2])),
                          arcs =  diff.edge,col = "darkblue",lwd = 3, lty = "dashed"))
    dev.off()
  
  ##################################################
  # infer the conditioning probabilty  for all model
  ##################################################
    bn.est =bn.fit(bn.learn,data = data.list[[i]])
    junction = compile(as.grain(bn.est))
    pa.ch=ancestors.search(bn.learn,targets)
    #exact.query(pa.ch$from,pa.ch$to)
    prob.list=appro.query(bn.est,pa.ch$from,pa.ch$to)
    interaction.prob=data.frame("from"=pa.ch$from,"to"=pa.ch$to,"prob"=unlist(prob.list))
    filename=sprintf('./simulated_data/bn_learn_interaction_%d_%d.csv',size.list[i],kk)
    write.csv(interaction.prob,filename)  # saving
    
  }
}




appro.query<-function(bn,pa.list,ch.list,method="ls"){
  if (method=="ls"){
    str1=paste("(", pa.list, " == '",as.character(0), "')", sep = "")
    str2 = paste("(", ch.list, " == '",as.character(1), "')", sep = "")
    prob=vector("list",length(str1))
    for(k in seq_along(str1)){
      prob[[k]]=cpquery(bn.est, eval(parse(text=str2[[k]])), eval(parse(text=str1[[k]])),method = "ls")}}
  else if (method=="lw"){
    str1=as.list(factor(rep(0,length(pa.list)),levels = c('0','1')))
    names(str1)=pa.list
    str2 = paste("(", ch.list, " == '",as.character(1), "')", sep = "")
    prob=vector("list",length(str2))
    for(k in seq_along(str1)){
      prob[[k]]=cpquery(bn.est, eval(parse(text=str2[[k]])),str1[k],method = "lw")}}
  return(prob)}
  

exact.query<-function(junc,pa.list,ch.list){
  query.grain=gRain::querygrain
  prob=map2(pa.list,ch.list,function(x,y){
  query.grain(junc,nodes = c(y,x), type = "conditional")})
  prob=map(prob,~.x[2,1])
  return (prob)}


ancestors.search<-function(bn,node.list) {
  ancestor=bnlearn::ancestors
  ancestor.list=map(node.list,function(x){ancestor(bn,x)})
  child.list=map2(ancestor.list,node.list,function(x,y){rep(y,length(x))})  
  ances_chil.list=data.frame("from" =unlist(ancestor.list), "to" =unlist(child.list))
  return( ances_chil.list)
}
dsep.search<-function(model,query.list,nodes_list=model$nodes){
  pairx.list=vector("list" ,length(query.list))
  sep=function(query){
    result=unlist(Filter(function (x) !dsep(bn =model ,x = x,y = query),setdiff(nodes_list,query)))
    print(sprintf('query %s',query))
    names(result)=NULL
    return(result)}
  pairx.list=map(query.list,sep)
  names(pairx.list)=NULL
  pairy.list=map2(pairx.list,query.list,function(x,y){rep(y,length(x))}) 
  pair.list =data.frame("from" =unlist(pairx.list), "to" =unlist(pairy.list))
  return (pair.list)
}

  
  

