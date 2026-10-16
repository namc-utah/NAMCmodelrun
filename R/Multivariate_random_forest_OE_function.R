multivariate_OE_calc=function(bugnew,prednew,ranfor.mod,Pc){
  datanew=cbind.data.frame(bugnew,prednew)
  testpredictions=randomForestSRC::predict.rfsrc(ranfor.mod,datanew)
  testresults=data.frame(sapply(testpredictions$regrOutput, function(x) x$predicted))
  names(testresults) <- names(testpredictions$regrOutput)
  #testresults=data.frame(sapply(testpredictions$regrOutput, function(x) x$predicted))
  #test=as.data.frame(testpredictions$regrOutput)
  #testresults<-sapply(test, function(x) x$predicted)
  Os=list()
  Es=list()
  BC=list()
  for(i in 1:nrow(testresults)) {
    cur.prd<-testresults[i,]; #vector of predicted probs for current sample;
    spdyn<-names(cur.prd)[cur.prd>=Pc];  #subset of taxa with Pi>=Pcutoff;
    cur.prd<-cur.prd[spdyn]
    cur.obs<-datanew[i,spdyn]; #vector of observed P/A for those taxa;
    Os[[i]]<-sum(cur.obs); #observed richness (O);
    Es[[i]]<-sum(cur.prd); #Expected richness (E);
    BC[[i]]<-sum(abs(cur.obs-cur.prd))/ (sum(cur.obs)+sum(cur.prd)); #BC value;
  }
  OE=as.data.frame(cbind(unlist(Os),unlist(Es),unlist(BC)))

  names(OE)<-c('O','E','BC')
  OE$OoverE=OE$O/OE$E
  row.names(OE)=row.names(bugnew)
  return(OE)
}
