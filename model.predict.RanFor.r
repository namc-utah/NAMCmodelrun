# R code to make predictions of O/E for a new set of sites,;
#   based on a Random Forest predictive model:

# J. Van Sickle, 02/25/10;
###########################;
#STEP 1 -- # Obtain the needed R objects;

# Next, create the 2 objects needed to define the data for new samples/sites at which predictions are desired;
#   "prednew" = data frame containing predictor variables (columns) at all new sites/samples (rows);
#                at which predictions are desired (e.g., 'test' sites);
#                See predcal for desired format;
#                Columns must include all predictors named in preds.final;
#                Do not include a site/sample in this data frame if it has missing data for 1 or more of the final predictors;
#   "bugnew.pa" = corresponding data frame of sites (rows) by taxa (columns) of observed
#               presence/absence (coded 1/0) at new samples/sites;
#               Columns must include all taxa (columns )in bugall.pa;
# bugnew and prednew must have site/sample ID's as their row names;

#Sample usage:
    # OE.assess<-model.predict.RF(bugcal.pa,grps.final,preds.final,rf.mod,prednew=pred.vld,bugnew.pa=bug.vld.pa,Pc);

## See Step 7 in model.build.RanFor.r for more details on making RanFor predictions;

#The function prints out basic statistics of O/E to the screen;
#The function result (OE.assess) is a list with 2 components:
# "OE.scores" -- A data frame for all samples, containing O, E, O/E and BC from the predictive and null models, as well as outlier flags;
# "Capture.Probs" - Matrix of model-predicted capture (occurrence) probabilties for all taxa in all samples;
###################;

#STEP 2 -- Source or compile the following code;  --
# All arguments are defined above, except Pc = the cutoff value;
# of predicted capture probability, for a taxon to be included in O/E calculations;
# A typical value is Pc=0.5. ;
# Set Pc equal to a very small positive number, such as 1.E-14, to have all reference taxa included in O/E;

#' OE predict version 0
#'
#' @param bugcal.pa 
#' @param grps.final 
#' @param preds.final 
#' @param ranfor.mod 
#' @param prednew 
#' @param bugnew.pa 
#' @param Pc 
#'
#' @return O, E, OoverE, BC, null 
#' @export
#'
#' @examples
model.predict.RanFor<-function(bugcal.pa,grps.final,preds.final,ranfor.mod, prednew,bugnew.pa,Pc=0.5) {;

#1. - initial definitions;
names(grps.final)<-row.names(bugcal.pa);
nsite.cal<-length(grps.final); #number of calibration sites;
npreds<-length(preds.final); #number of predictor variables;
grpsiz<-table(grps.final); #tabulate sites per group;
ngrps<-length(grpsiz);  #number of groups;

#2. Alignment of new predictor and bug data;
  # Need to align the rows (samples) of the new bug data to the new predictor data;
  #Also, keep only those taxa (columns) that are in the model calibration data set, and put them in same order;
  #     This taxon alignment is required for BC calculation. ;
 bugnew.pa<-bugnew.pa[row.names(prednew),names(bugcal.pa)];
 ################;

## STEP 3. -- Use RF to predict the group (cluster) membership for all new sites. ;
# Does not use RIVPACS assumption of weighting the membership probabilities by Calibration group size, as a prior probability;
# Also, RF predictions do not have an outlier test, unlike DFA predictions;
# Predicted probs are outputted as a matrix, sites are rows, columns are groups;

grpprobs<-predict(ranfor.mod,newdata=prednew[,preds.final],type='prob');

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
   spdyn<-names(cur.prd)[cur.prd>=Pc];  #subset of taxa with Pi>=Pcutoff;
   cur.prd<-cur.prd[spdyn]; #vector of Pi for subset of included taxa;
   cur.obs<-bugnew.pa[i,spdyn]; #vector of observed P/A for those taxa;
   OE.stats$OBS[i]<-sum(cur.obs); #observed richness (O);
   OE.stats$E.prd[i]<-sum(cur.prd); #Expected richness (E);
   OE.stats$BC.prd[i]<-sum(abs(cur.obs-cur.prd))/ (OE.stats$OBS[i]+OE.stats$E.prd[i]); #BC value;
         }; #finish sample loop;

#5.2 - Compute Expected richness (E) and BC for null model using taxa >= Pc.
# Note that the set of taxa included in the null model is fixed for all samples;
pnull<-apply(bugcal.pa,2,sum)/dim(bugcal.pa)[[1]];  #null model predicted occurrnece probabilities, all taxa;
nulltax<-names(pnull[pnull>=Pc]); #subset of taxa with Pnull >= Pc;
Enull<-sum(pnull[nulltax]);
Obsnull<-apply(bugnew.pa[,nulltax],1,sum); #vector of Observed richness, new samples, under null model;
BC.null<-apply(bugnew.pa[,nulltax],1,function(x)sum(abs(x-pnull[nulltax])))/(Obsnull+Enull); #vector of null-model BC;

#5.3 - Final data frame contains values of O, E, O/E, Onull, Enull, Onull/Enull, BC.prd and BC.null, for all samples;
#Also includes outlier flags;

OE.final<-data.frame(O=OE.stats$OBS,E=OE.stats$E.prd,
                      OoverE=OE.stats$OBS/OE.stats$E.prd, 
   Onull=Obsnull,Enull=rep(Enull,length(Obsnull)),OoverE.null=Obsnull/Enull,
   BC= OE.stats$BC.prd,BC.null=BC.null,
          row.names=row.names(bugnew.pa));
###########;
#print some summary statistics of O/E to wrap up;
print(' ',quote=F)
print('Statistics of O/E for new samples',quote=F);
print(' ',quote=F)
  c1<-mean(OE.final$OoverE); c2<-mean(OE.final$OoverE.null);
  s1<-sqrt(var(OE.final$OoverE)); s2<-sqrt(var(OE.final$OoverE.null));
print(' Mean(O/E) and SD(O/E), from predictive model: ',quote=F)
print( c(c1,s1),digits=3);
print(' ',quote=F)
print(' Mean(O/E) and SD(O/E), from null model: ',quote=F)
print( c(c2,s2),digits=3);
#print outlier count;
print(' ',quote=F)
print('RanFor predictions complete', quote=F);
#return a list with OE.final and the matrix of predicted probs as elements;
list(OE.scores=OE.final,Capture.Probs=site.pred.dfa); #return data frame as final object;
 }; #end of function;



