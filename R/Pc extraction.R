#get Pc for a client
#1. Source the appropriate code, in this case, RF 4.2, to get important objects
#2. Run the Run_models.R code to ensure that you have the appropriate bugnew
#and to get the model results saved to the database.
#4. Ctrl+A and run this code.
Cal.OOB=F
Eprd<-list()

  #first convert bug matrix to P/A (1/0);
  temp.pa<-bugnew;
  temp.pa[temp.pa>0]<-1;
  #rm(bugnew);

  #1. - initial definitions;
  names(grps.final)<-row.names(bugcal.pa);
  nsite.cal<-length(grps.final); #number of calibration sites;
  npreds<-length(preds.final); #number of predictor variables;
  grpsiz<-table(grps.final); #tabulate sites per group;
  ngrps<-length(grpsiz);  #number of groups;

  #2. Alignment of new predictor and bug data with model data;
  #2a) Align the rows (samples) of the new bug data to the new predictor data;
  temp.pa<-temp.pa[row.names(prednew),];
  #2b)reshape bugnew columns (taxa) to match those in bugcal.pa, and be in same order;
  # New bug data might have fewer or more columns (taxa) than calibration bug data;
  # create a new empty site x taxa matrix with bugcal.pa columns and bugnew rows, fill it with zeros;
  nsite.new<-dim(temp.pa)[[1]];
  ntaxa<-dim(bugcal.pa)[[2]];
  bugnew.pa<-matrix(rep(0,nsite.new*ntaxa),nrow=nsite.new,dimnames=list(rownames(temp.pa),names(bugcal.pa)));
  #loop through columns of new matrix and fill with columns of the original test data matrix;
  col.match<-match(dimnames(bugnew.pa)[[2]],dimnames(temp.pa)[[2]]);
  for(kcol in 1:ntaxa) if(!is.na(col.match[kcol]))bugnew.pa[,kcol]<-temp.pa[,col.match[kcol]];
  ################;

  ## STEP 3. -- Use RF to predict the group (cluster) membership for all new sites. ;
  # Does not use RIVPACS assumption of weighting the membership probabilities by Calibration group size, as a prior probability;
  # Also, RF predictions do not have an outlier test, unlike DFA predictions;
  # Predicted probs are outputted as a matrix, sites are rows, columns are groups;

  #If Cal.OOB is true, do OOB predictions, appropriate ONLY for CAL data;
  # If it is false, do a new prediction;
  if(Cal.OOB==TRUE) grpprobs<-ranfor.mod$votes else grpprobs<-predict(ranfor.mod,newdata=prednew[,preds.final],type='prob');

  ############;
  #STEP 4 -- Compute predicted occurrence probability for each modeled taxon at each new sample;
  # "modeled OTU's" consist of all taxa that were found at >=1 calibration sample;
  #To do this, first calculate the occurrence freqs of all modeled taxa in the Calibration sample groups;
  grpocc<-apply(bugcal.pa,2,function(x){tapply(x,grps.final,function(y){sum(y)/length(y)})});

  #finally, compute the matrix of predicted occurrence (capture) probabilities, for all new samples and all modeled taxa;
  #This is the matrix-algebra form of the RIVPACS combining formula (e.g., Clarke et al. 2003, Eq. 4)
  site.pred.dfa<-grpprobs%*%grpocc;

  #######################;

  # STEP 5. Compute O, E, O/E and BC for all samples. ;
  # Also compute O/E and BC for the null model;

  #5.1 loop over all samples. Compute and store  O, predicted E, predicted BC for each sample. ;
  #temporary data frame to hold nonnull results for all samples. ;
  nsit.new<-dim(prednew)[[1]];
  OE.stats<-data.frame(OBS=rep(NA,nsit.new), E.prd=rep(NA,nsit.new),BC.prd=rep(NA,nsit.new),row.names=row.names(prednew));
  for(i in 1:nsit.new) {;
    #i<-1;
    cur.prd<-site.pred.dfa[i,]; #vector of predicted probs for current sample;
    #spdyn<-names(cur.prd)[cur.prd>=Pc];  #subset of taxa with Pi>=Pcutoff;
    #cur.prd<-cur.prd[spdyn]; #vector of Pi for subset of included taxa;
    #cur.obs<-bugnew.pa[i,spdyn]; #vector of observed P/A for those taxa;
    #OE.stats$OBS[i]<-sum(cur.obs); #observed richness (O);
    #OE.stats$E.prd[i]<-sum(cur.prd); #Expected richness (E);
    #OE.stats$BC.prd[i]<-sum(abs(cur.obs-cur.prd))/ (OE.stats$OBS[i]+OE.stats$E.prd[i]); #BC value;
  Eprd[[i]]<-cur.prd
    }; #finish sample loop;




pc_list<-list()
for(z in 1:length(Eprd)){
  finalPc<-Eprd[[z]]
  pp<-data.frame(Taxon=names(finalPc),
                 Pc=unname(finalPc))
  sampID<-row.names(bugnew)[z]
  pp$sampleId<-sampID
  pc_list[[z]]<-pp
  message(paste(z, ' of ',length(Eprd)))
}
#
library(dplyr)
dt_list<-lapply(pc_list,data.table::as.data.table)
 long_dt<-data.table::rbindlist(dt_list)
 Pcs<-as.data.frame(data.table::dcast(long_dt,sampleId~Taxon,value.var = 'Pc'))
 row.names(Pcs)<-Pcs$sampleId;Pcs<-Pcs[,-1]
 Os<-bugnew.pa[,order(colnames(bugnew.pa))]
 #write.csv(do.call(rbind,pc_list),'C://Users//andrew.caudillo//Box//NAMC//OutgoingData//Data exports//AREMP//AREMPOEs_Pcs.csv')
 #write.csv(Os,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Observed_PA.csv')
 #write.csv(Pcs,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Expected.csv')
  write.csv(Os,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Prob_O.csv')
