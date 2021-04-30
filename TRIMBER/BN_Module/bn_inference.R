args = commandArgs(trailingOnly=TRUE)
ks_test<-function(expression.data,expression,interaction,alpha=0.05)
{ 
  #expression.data :raw data with continous value 
  #expression     :binarized data 
  index=rep(0, length(interaction$from))
  for  (i in seq_along(targets)){
    te=as.numeric(paste(expression.data[,interaction$to[[i]]]))  
    tec1=as.numeric(paste(expression[,interaction$from[[i]]]))
    if ((length(te[tec1 == 1])>0) & (length(te[tec1 == 0])>0)){
      ind=ks.test(te[tec1 == 1],te[tec1 == 0])
      print(ind$p.value)
      if (ind$p.value> alpha){index[[i]]=1}
    }
    else{index[[i]]=1}
  }
  return(index)
}

dsep.search<-function(model,query.list,nodes_list){
  pairx.list=vector("list" ,length(query.list))
  sep=function(query){
    result=unlist(Filter(function (x) !dsep(bn =model ,x = x,y = query),dplyr::setdiff(nodes_list,query)))
    print(sprintf('query %s',query))
    names(result)=NULL
    return(result)}
  pairx.list=purrr::map(query.list,function(x) {sep(x)})
  names(pairx.list)=NULL
  pairy.list=purrr::map2(pairx.list,query.list,function(x,y){rep(y,length(x))}) 
  pair.list =data.frame("from" =unlist(pairx.list), "to" =unlist(pairy.list))
  return (pair.list)
}
ancestors.search<-function(bn,node.list) {
  ancestor=bnlearn::ancestors
  ancestor.list=purrr::map(node.list,function(x){ancestor(bn,x)})
  child.list=purrr::map2(ancestor.list,node.list,function(x,y){rep(y,length(x))})  
  ances_chil.list=data.frame("from" =unlist(ancestor.list), "to" =unlist(child.list))
  return( ances_chil.list)
}

Joint.query <-function (bn,ko.tf,pa.ch,rxnGeneMat,rules,max_gene=10,mode='BN')
{
  # a prob vector indicating the prob associated with each rxn wiil be returned.
  # prob is set to be 1 for rxn that is not affected 
  # max_gene: the max number of genes for prob table computation 
  # if max_gene is bigger than 10(default), prob vector as in PROM will be computed.
  # BN :prob table CN:prob_vector
  reg.list=ko.tf
  genes=row.names(rxnGeneMat)
  rxns=colnames(rxnGeneMat) 
  
  prob.list=as.data.frame(matrix(1,nrow=length(rxns),ncol = length(reg.list)))
  names(prob.list)=reg.list
  for (j in seq_along(reg.list)){
    reg=unlist(reg.list[[j]])
    prob=rep(1,length(rxns))
    names(prob)= rxns
    #find genes regulated and rxns affected for each all TF
    tar=pa.ch[which(pa.ch$from %in% reg),'to']
    if (all(!(tar %in% genes))){next}
    tar=tar[tar %in% genes]
    rxn.index=purrr::keep(rxns,~any(as.logical(rxnGeneMat[tar,.x])))
    #Figure out if any of the reaction states is changed
    rxn.index=states_check(rxn.index,tar,genes,rules)
    ##############
    evidence=paste("(", reg, " == '",as.character(0), "')", sep = "",collapse = " & ")
    for (i in rxn.index){
      prob[[i]]=0
      gene.ko =dplyr::intersect(genes[which(rxnGeneMat[,i]==1)],tar)
      gene.num=sum(length( gene.ko))
      print(sprintf('TF:%s -> reaction :%s regulated by %d genes',paste(reg,sep = "",collapse = " & ") ,i,gene.num))
      #compute the minimum  of  prob vector if gene.num >max_gene or model =='CN'
      if (gene.num>10 || mode=='CN')
      {tmp.pach=data.frame('from'=rep(reg,gene.num),'to'=gene.ko)
      prob[[i]]=max(unlist(appro.query(bn,tmp.pach)))
      next}
      #compute the summation of valid states' probs in prob table .
      state.list = binaryLogic::as.binary(1:(2^gene.num-1),n=gene.num)
      state.list=states_filter(state.list,gene.ko,genes,rules[[i]])
      for (state in state.list){
        gene.on=gene.ko[state] 
        gene.off=gene.ko[!state]
        event.on=purrr::map(gene.on,~paste("(",.x,"==",as.character(1),")"))
        event.off=purrr::map(gene.off,~paste("(",.x,"==",as.character(0),")"))
        event=c(event.on,event.off)
        str= paste(c(event.on,event.off),sep = "",collapse = " & ")
        cmd=paste("cpquery(bn.est, ", str, ", ", evidence, ")", sep = "")
        prob[[i]]= prob[[i]]+eval(parse(text = cmd))
      }
    }
    prob.list[,j]=prob
  }
  return(prob.list)
}


states_check<-function(rxns,tar,genes,rules,mode=FALSE){
  #Mode=TRUE for  unchanged rxns selection
  #Mode=FALSE for changed   rxns
  x=rep(TRUE,length(genes))
  names(x)=genes
  x[tar]=FALSE
  if(mode){
    constrainRxn=unlist(purrr::map(rxns,~eval(parse(text=rules[[.x]]))))
    rxns=rxns[constrainRxn]}
  else{
    constrainRxn=unlist(purrr::map(rxns,~!eval(parse(text=rules[[.x]]))))
    rxns=rxns[constrainRxn]}
  return(rxns)
}

states_filter<-function(state.list,tar,genes,bool_rule){
  
  x=rep(TRUE,length(genes))
  s=rep(FALSE,length(state.list))
  names(x)=genes
  for (k in seq_along(state.list)){
    x[tar]=FALSE
    x[tar[state.list[[k]]]]=TRUE 
    s[[k]]=eval(parse(text=bool_rule))
  }
  state.list=state.list[s]
  return(state.list)
}

appro.query<-function(bn,pa.ch){
  pa.list=pa.ch$from
  ch.list=pa.ch$to
  str1=paste("(", pa.list, " == '",as.character(0), "')", sep = "")
  str2 = paste("(", ch.list, " == '",as.character(1), "')", sep = "")
  prob=vector("list",length(str1))
  for(k in seq_along(str1)){print(sprintf('interaction no:%d',k))
    cmd=paste("cpquery(bn, ", str2[[k]], ", ", str1[[k]],")", sep = "")
    prob[[k]]=eval(parse(text=cmd))}
  interaction.prob=data.frame("from"=pa.ch$from,"to"=pa.ch$to,"prob"=unlist(prob))
  return(interaction.prob)}

exact.query<-function(bn,pa.ch){
  pa.list=pa.ch$from
  ch.list=pa.ch$to
  junc = compile(as.grain(bn))
  query.grain=gRain::querygrain
  prob=purrr::map2(pa.list,ch.list,function(x,y){
    query.grain(junc,nodes = c(y,x), type = "conditional")})
  prob=purrr::map(prob,~.x[2,1])
  interaction.prob=data.frame("from"=pa.ch$from,"to"=pa.ch$to,"prob"=unlist(prob))
  
  return (interaction.prob)}

instant_pkgs <- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  print(pkgs_miss)
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss,quiet = TRUE)
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed.\n")
  }
}

instant_pkgs(c("bnlearn", "parallel","dplyr","purrr","R.matlab"))
library(bnlearn)
print('data  loading \n')
model<-R.matlab::readMat(args[1])
print('model  loading \n')
bn.est<-read.bif(args[2])
bn.learn<-bn.net(bn.est)
######################
model.genes= unlist(model$genes)
model.rxns=  unlist(model$rxns)
model.rxnGeneMat=as.data.frame(t(model$rxnGeneMat),row.names = model.genes)
colnames(model.rxnGeneMat)=model.rxns
model.rules=model$rules
model.rules=purrr::map(model.rules,~unlist(.x))
names(model.rules)=model.rxns
ko.tf=lapply(model$ko.tf,function (x) unlist(x))
#######################
interaction.reg=unlist(model$regulator)
interaction.tar=unlist(model$targets)
interaction=which(interaction.reg %in% unique(unlist(ko.tf)))
interaction=data.frame('from'=interaction.reg[interaction],'to'=interaction.tar[interaction])
#######################
#selection
pa.ch=unlist(purrr::map2(interaction$from,interaction$to,
                  function (x,y){!dsep(bn =  bn.learn ,x ,y)}))
pa.ch=interaction[pa.ch,]


prob.list=Joint.query(bn.est,ko.tf,pa.ch,model.rxnGeneMat,model.rules,mode = args[3])
filename=sprintf('%s/bn_learn_interaction.csv',args[4])
write.csv(prob.list,filename)  # saving








