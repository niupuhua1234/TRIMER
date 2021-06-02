setwd('C:/Users/sweet/Desktop/ts/trimer')
interaction.preprocess<-function(regulators,targets)
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
bnlearning<-function(expression,regulators,targets,interaction,space='specific',
                     algorithm,rep,start=vector("list",0),debug=TRUE)
{
  #general : no prior knowledege about network structure
  #specific: using interaction list as candidate edges
  #rep : number of  model to be learned for model averaging 
  #space: seaching space (default all possbile arcs wihout arcs from target to regualtor)
  White.list=data.frame("from"=c("b4063","b4004"),"to"=c("b2601","b3281"))
  TFnum=length(regulators)
  targetnum=length(targets)
  all_gene=c(regulators,targets)
  if (space=='general'){
    # Structure constraint
    #gene can not be parent of regulator
    #gene can not be parent of gene
    black.list=data.frame('from'= unlist(map(targets,~unname(rep(.x,TFnum)))),'to'  = unlist(rep(regulators,targetnum)))
    black.list=rbind(black.list,set2blacklist(targets))}
  else{
    black.list=data.frame(set2blacklist(all_gene))
    interaction.edge=data.frame("from"=c(interaction$from,interaction$to),
                                "to"=c(interaction$to,interaction$from))
    black.list=setdiff(black.list,interaction.edge)
  }
  if (algorithm  == "tabu") {
    if (length(start) == 0){
      start=empty.graph(nodes = all_gene,1)}
    if (rep==1){
      bn.learn=tabu(expression, blacklist =black.list,tabu = 100,score = 'bic',start = start ,debug = debug)}
    else{
      boot = boot.strength(data = expression, R =rep,algorithm = "hc",
                           algorithm.args = list(score = "bic", blacklist =black.list,whitelist =White.list,tabu = 100,start=start,debug = debug))
      bn   = cextend(averaged.network(boot))
    }
  }
  if (algorithm  == "aracne") {
    bn.learn=aracne(expression,debug = debug)
    bn.learn=pdag2dag(bn,names(expression))
  }
  if (algorithm  == "chow_liu") {
    bn.learn=chow.liu(expression,debug = debug)
    bn.learn=pdag2dag(bn,names(expression))
  }
  
  return (bn.learn)}
instant_pkgs <- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  print(pkgs_miss)
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed.\n")
  }
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded.\n")
  }
}
instant_pkgs(c("bnlearn","dplyr","purrr","R.matlab"))
print('data  loading \n')
iaf1260<- readMat("./Output_data/Ecoli_model_EcoMac.mat")
expression.data =readMat('./Output_data/EcoMac_data_binary.mat')
print('data arrangement\n')
expression.data =expression.data$expression
expression.data=as.data.frame(t(expression.data))
expression.id=unlist(iaf1260$expressionid)
expression.name=unlist(iaf1260$expressionname)
names(expression.data)=expression.id
expression=data.frame(lapply(expression.data,function (x) factor(x,levels = c('0','1'))))
#############################
interaction.reg=unlist(iaf1260$regulator)
interaction.tar=unlist(iaf1260$targets)
regulators=unique(interaction.reg)
targets=unique(interaction.tar)
targets=targets[-which(targets %in% regulators)]
#########################
interaction=interaction.preprocess(interaction.reg,interaction.tar)
expression.est=expression[,c(regulators,targets,'b3281')]
#############################
print('structrue learning ')
#blacklisted arcs are never present in the learned graph.
#arcs can only be whitelisted in a single direction== fixed 

#soxR - aroF activator  b4063 b2601
#ZraR - aroE activator  b4004 b3281

White.list=data.frame("from"=c("b4063","b4004"),"to"=c("b2601","b3281"))

#for (i in seq(10)){
#  start= empty.graph(nodes=c(regulators,targets,'b3281'),1)
#  arcs(start)=arcs(random.graph(nodes =c(regulators), method = "ic-dag", max.degree=10))
#  bn.learn=bnlearning(expression.est,regulators,targets,interaction$edge,space='specified',algorithm='tabu',rep=1,start= start)
#  bn.est=bn.fit(bn.learn,data = expression.est)  
#  filename=sprintf('learned_network_en_%d.bif',i)
#  write.bif(file=  filename,bn.est)
#}



bn.learn=bn.net(read.bif('learned_network_1.bif'))
num.perturb=ceiling(dim(arcs(bn.learn))[1]*0.05)
for (i in seq(10))
{
  bn.new=bn.learn
  arc.per=arcs(bn.learn)
  order=sample(dim(arc.per)[1],num.perturb)
  count=1
  for (j in seq(num.perturb)){
    bn.old=bn.new
    bn.new=reverse.arc(bn.old, arc.per[ order[j],1], arc.per[ order[j],2],check.cycles = FALSE)
    while (!acyclic(bn.new)){
      new_order=sample(dim(arc.per)[1],1)
      while (new_order %in% order){
        new_order=sample(dim(arc.per)[1],1)}
      order=c(order,new_order)
      bn.new=reverse.arc(bn.old, arc.per[new_order,1], arc.per[new_order,2],check.cycles = FALSE)}
  }
  bn.est=bn.fit(bn.new,data = expression.est)  
  filename=sprintf('learned_network_%d_per.bif',i)
  write.bif(file=  filename,bn.est)
}

  
  








