appro.query<-function(bn,pa.ch){
  pa.list=pa.ch$from
  ch.list=pa.ch$to
  str1=paste("(", pa.list, " == '",as.character(0), "')", sep = "")
  str2 = paste("(", ch.list, " == '",as.character(1), "')", sep = "")
  prob=vector("list",length(str1))
  for(k in seq_along(str1)){print(sprintf('interaction no:%d',k))
    cmd=paste("cpquery(bn, ", str2[[k]], ", ", str1[[k]],")", sep = "")
    prob[[k]]=eval(parse(text=cmd))}
  return(prob)}

exact.query<-function(bn,pa.ch){
  pa.list=pa.ch$from
  ch.list=pa.ch$to
  junc = compile(as.grain(bn))
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
dsep.search<-function(model,query.list,nodes_list){
  setdiff=dplyr::setdiff
  pairx.list=vector("list" ,length(query.list))
  sep=function(query){
    result=unlist(Filter(function (x) !dsep(bn =model ,x = x,y = query),setdiff(nodes_list,query)))
    print(sprintf('query %s',query))
    names(result)=NULL
    return(result)}
  pairx.list=map(query.list,function(x) {sep(x)})
  names(pairx.list)=NULL
  pairy.list=map2(pairx.list,query.list,function(x,y){rep(y,length(x))}) 
  pair.list =data.frame("from" =unlist(pairx.list), "to" =unlist(pairy.list))
  return (pair.list)
}

ks_test<-function(expression.data,expression,interaction,alpha=0.05)
  {
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

Joint.query <-function (bn,pa.ch,rxnGeneMat,rules,max_gene=10)
{
  # a prob vector indicating the prob associated with each rxn wiil be returned.
  # prob is set to be 1 for rxn that is not affected 
  # max_gene: the max number of genes for prob table computation 
  # if max_gene is bigger than 10(default), prob vector as in PROM will be computed.
  reg.list=unique(pa.ch$from)
  genes=row.names(rxnGeneMat)
  rxns=colnames(rxnGeneMat) 

  prob.list=as.data.frame(matrix(1,nrow=length(rxns),ncol = length(reg.list),
                                 dimnames = list(NULL,reg.list)))
  names(prob.list)=reg.list
  for (reg in reg.list){
    prob=rep(1,length(rxns))
    names(prob)= rxns
    #find genes regulated and rxns affected for each TF
    tar=pa.ch[which(pa.ch$from==reg),'to']
    if (all(!(tar %in% genes))){next}
    tar=tar[tar %in% genes]
    rxn.index=keep(rxns,~any(as.logical(rxnGeneMat[tar,.x])))
    #Figure out if any of the reaction states is changed
    rxn.index=states_check(rxn.index,tar,genes,rules)
    ##############
    evidence=paste("(", reg, " == '",as.character(0), "')", sep = "")
    for (i in rxn.index){
      prob[[i]]=0
      gene.ko =intersect(genes[which(rxnGeneMat[,i]==1)],tar)
      gene.num=sum(length( gene.ko))
      print(sprintf('TF:%s -> reaction :%s regulated by %d genes',reg ,i,gene.num))
      #compute the minimum  of  prob vector if gene.num >max_gene
      if (gene.num>10){
        tmp.pach=pa.ch[which(reg == pa.ch$from),]
        prob[[i]]=max(unlist(appro.query(bn,tmp.pach)))
        next}
      #compute the summation of valid states' probs in prob table .
      state.list = binaryLogic::as.binary(1:(2^gene.num-1),n=gene.num)
      state.list=states_filter(state.list,gene.ko,genes,rules[[i]])
      for (state in state.list){
        gene.on=gene.ko[state] 
        gene.off=gene.ko[!state]
        event.on=map(gene.on,~paste("(",.x,"==",as.character(1),")"))
        event.off=map(gene.off,~paste("(",.x,"==",as.character(0),")"))
        event=c(event.on,event.off)
        str= paste(c(event.on,event.off),sep = "",collapse = " & ")
        cmd=paste("cpquery(bn.est, ", str, ", ", evidence, ")", sep = "")
        prob[[i]]= prob[[i]]+eval(parse(text = cmd))
        }
    }
    prob.list[,reg]=prob
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
    constrainRxn=unlist(map(rxns,~eval(parse(text=rules[[.x]]))))
    rxns=rxns[constrainRxn]}
  else{
    constrainRxn=unlist(map(rxns,~!eval(parse(text=rules[[.x]]))))
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


bnlearning<-function(expression,regulators,targets,interaction,space,algorithm,rep){
  TFnum=length(regulators)
  targetnum=length(targets)
  if (space=='general'){
    # Structure constraint
    #gene can not be parent of regulator
    #gene can not be parent of gene
    black.list=data.frame('from'= unlist(map(targets,~unname(rep(.x,TFnum)))),'to'  = unlist(rep(regulators,targetnum)))
    black.list=rbind(black.list,set2blacklist(targets))}
  else{
    all_gene=c(regulators,targets)
    black.list=data.frame(set2blacklist(all_gene))
    interaction.edge=data.frame("from"=c(interaction$from,interaction$to),
                                "to"=c(interaction$to,interaction$from))
    black.list=setdiff(black.list,interaction.edge)
  }
  if (algorithm  == "tabu") {
    if (rep==1){
      bn.learn=tabu(expression, blacklist =black.list,tabu = 100,score = 'bic',debug = TRUE)}
    else{
      boot = boot.strength(data = expression, R =rep,algorithm = "hc",
                           algorithm.args = list(score = "bic", blacklist =black.list,tabu = 100,debug = TRUE))
      bn   = cextend(averaged.network(boot))
      }
  }
  if (algorithm  == "aracne") {
    bn.learn=aracne(expression,debug = TRUE)
    bn.learn=pdag2dag(bn,names(expression))
  }
  if (algorithm  == "chow_liu") {
    bn.learn=chow.liu(expression,debug = TRUE)
    bn.learn=pdag2dag(bn,names(expression))
  }
  
  return (bn.learn)}

