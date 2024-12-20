# R code to make predictions of O/E for a new set of sites,;
#     based on a single 'final' predictive model:

#Version 3, June 25, 2007 -- Includes BC index;

#Version 4 - Restructured to use a final model from model.build.v4.r as its input;

#version 4.1 - add code to
   # a)align the test bug and predictor data, by sample;
   #b) convert bug data matrix to P/A (1/0), in case it contains abundances;

###########################;
#STEP 1 -- # Obtain the needed R objects;
 # first, obtain the five R objects which define a single model (see model.build.v4.r);
#   "bugcal.pa"= Matrix of observed presence/absence (1/0) at calibration sites, for all taxa.
#
#   "grps.final" = Corresponding vector identifying the group ID's of calibration sites in bugcal.pa;
#
#   "preds.final" = Vector with names of the chosen predictor variables, all of which must be available;
#
#   "grpmns" = Matrix of group (cluster) means at calibration sites;
#
#   "covpinv" = Inverse of pooled covariance matrix among final predictor variables at calibration sites;
#
# Next, create the 2 objects needed to define the data for new samples/sites at which predictions are desired;
#   "prednew" = data frame containing predictor variables (columns) at all new sites/samples (rows);
#                at which predictions are desired (e.g., 'test' sites);
#                Columns must include all predictors named in preds.final;
#                Do not include a site/sample in this data frame if it has missing data for 1 or more of the final predictors;
#   "bugnew" = Corresponding observed species matrix of sites (rows) by taxa (columns).
#              Matrix entries are either observed abundance (count), or else presence/absence (coded as 1/0);
#   bugnew and prednew must have site/sample ID's as their row names;

#Sample usage:
    # OE.assess<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.vld,bugnew=bug.vld.pa,Pc)

#The function prints out basic statistics of O/E to the screen;
#The function result (OE.assess) is a list with 3 components:
# "OE.scores" -- A data frame for all samples, containing O, E, O/E and BC from the predictive and null models, as well as outlier flags;
# "Capture.Probs" - Matrix of model-predicted capture (occurrence) probabilties for all taxa in all samples;
#  "Group.Occurrence.Probs" = Matrix of predicted probabilities of occurrence for each sample in each sample/site group;
###################;

#STEP 2 -- Source or compile the following code;  --
# All arguments are defined above, except Pc = the cutoff value;
# of predicted capture probability, for a taxon to be included in O/E calculations;
# A typical value is Pc=0.5. ;
# Set Pc equal to a very small positive number, such as 1.E-14, to have all reference taxa included in O/E;

#'  OE predict version 4.1
#'
#' @param bugcal.pa
#' @param grps.final
#' @param preds.final
#' @param grpmns
#' @param covpinv
#' @param prednew
#' @param bugnew
#' @param Pc
#'
#' @return O, E, OoverE, BC, null
#' @export
#'
#' @examples
model.predict.v4.1<-function(bugcal.pa,grps.final,preds.final,grpmns,covpinv,prednew,bugnew,Pc=0.5) {

#first convert bug matrix to P/A (1/0);
   temp.pa<-bugnew;
   temp.pa[temp.pa>0]<-1;
   rm(bugnew);

#1. - initial definitions;
names(grps.final)<-row.names(bugcal.pa);
nsite.cal<-length(grps.final); #number of calibration sites;
npreds<-length(preds.final); #number of predictor variables;
grpsiz<-table(grps.final); #tabulate sites per group;
ngrps<-length(grpsiz);  #number of groups;

#2. Alignment of new predictor and bug data with model data;
   #2a) Align the rows (samples) of the new bug data to the new predictor data;
    temp.pa<-temp.pa[row.names(prednew),];
   #reshape bugnew.pa columns (taxa) to match those in bugcal.pa, and be in same order;
   # New bug data might have fewer or more columns, and it might have extraneous columns;
   # create a new empty site x taxa matrix with bugcal.pa columns and bugnew rows, fill it with zeros;
   nsite.new<-dim(temp.pa)[[1]];
   ntaxa<-dim(bugcal.pa)[[2]];
   bugnew.pa<-matrix(rep(0,nsite.new*ntaxa),nrow=nsite.new,dimnames=list(rownames(temp.pa),names(bugcal.pa)));
   #loop through columns of new matrix and fill with columns of the original test data matrix;
   col.match<-match(dimnames(bugnew.pa)[[2]],dimnames(temp.pa)[[2]]);
   for(kcol in 1:ntaxa) if(!is.na(col.match[kcol]))bugnew.pa[,kcol]<-temp.pa[,col.match[kcol]];
 ################;

## 3. -- predict the group (cluster) membership for all new sites. ;
# Follow RIVPACS assumption of weighting ;
# the membership probabilities by Calibration group size, as a prior probability;
# Also, flag any outlier sites, using the chi-squared statistic;
#Store probs in matrix, sites are rows, columns are groups;
#Uses mahalanobis function, where new vector is taken as the 'center', mu,;
#and matrix of means is taken as the 'data matrix', x;

dmat<-as.matrix(prednew[,preds.final]);
dmat<-apply(dmat,2,as.numeric);row.names(dmat)<-row.names(prednew)#matrix of predictor data for new samples, include only the predictor variables;

#3.1 -- compute the critical chi-squared values for flagging outlier samples;
# df = the MINIMUM of (a)(number of groups-1), and (b) number of predictor variables;
# will flag each sample at P-value =.05 and also P-value =.01 level;
dff<-min(c(npreds,(ngrps-1)));
crit.01<-qchisq(0.99,df=dff);
crit.05<-qchisq(0.95,df=dff);

# 3.2 - construct empty matrix for predicted membership probabilities;
nsit.new<-dim(dmat)[[1]]; #number of new samples;
grpprobs<-matrix(rep(0,nsit.new*ngrps),nrow=nsit.new,
            dimnames=list(dimnames(dmat)[[1]],dimnames(grpmns)[[1]]));

#Also construct empty data.frame for outlier flag;
# include minimum (squared)distance;
# Each sample is either a PASS (denote by 0) or FAIL (denote by 1) for the outlier test;
 outlier.flag<-data.frame(outlier.05=rep(0,nsit.new),outlier.01=rep(0,nsit.new),dismin=rep(0,nsit.new),row.names=dimnames(dmat)[[1]]);

# 3.3 - Loop over ALL samples, compute vector of group membership probs for each sample and set its outlier flags;
#execute the following code piece as a single block;
##;
for(i in 1:nsit.new){;
#vector of squared Mahal. dist from current sample to each group mean;
dist<-mahalanobis(grpmns,dmat[i,],covpinv,inverted=T); #vector of distances;
grpprobs[i,]<-grpsiz*exp(-0.5*dist); # see Clarke et al. (2000);
grpprobs[i,]<-grpprobs[i,]/sum(grpprobs[i,]);
#check for outlier;
outlier.flag$dismin[i]<-min(dist); #save minimum distance;
if(outlier.flag$dismin[i]>crit.05)outlier.flag[i,'outlier.05']<-1;
if(outlier.flag$dismin[i]>crit.01)outlier.flag[i,'outlier.01']<-1;
        }; #finish sample loop;
#### sample membership probabilities complete for new samples;

############;
#STEP 4 -- Compute predicted occurrence probability for each modeled taxon at each new sample;
# "modeled taxa" consist of all taxa that were found at >=1 calibration sample;
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
#this addition will be the workaround for sites with only 1 sample
#note the simple ifelse.
#if bugnew.pa is only 1 row of data,
#apply will not work because bugnew.pa
#technically has no dimensions,
#thus apply cannot pull the needed data.
Obsnull<-if (is.null(dim(bugnew.pa[,nulltax]))){
  sum(bugnew.pa[,nulltax])
}else{ rowSums(bugnew.pa[,nulltax])} #vector of Observed richness, new samples, under null model;
#this addition will be the workaround for sites with only 1 sample
#note the simple ifelse.
#if bugnew.pa is only 1 row of data,
#apply will not work because bugnew.pa
#technically has no dimensions,
#thus apply cannot pull the needed data.
BC.null<-if(is.null(dim(bugnew.pa[,nulltax]))){
  sum(abs(bugnew.pa[,nulltax]-pnull[nulltax]))/(Obsnull+Enull)
}else{
  apply(bugnew.pa[,nulltax],1,function(x)sum(abs(x-pnull[nulltax])))/(Obsnull+Enull)
} #vector of null-model BC;


#5.3 - Final data frame contains values of O, E, O/E, Onull, Enull, Onull/Enull, BC.prd and BC.null, for all samples;
#Also includes outlier flags;

OE.final<-data.frame(O=OE.stats$OBS,E=OE.stats$E.prd,
                      OoverE=OE.stats$OBS/OE.stats$E.prd,
   Onull=Obsnull,Enull=rep(Enull,length(Obsnull)),OoverE.null=Obsnull/Enull,
   BC= OE.stats$BC.prd,BC.null=BC.null,
     outlier.05=outlier.flag$outlier.05,outlier.01=outlier.flag$outlier.01,
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
print('Counts of OK (=0) and outlier (=1) samples, assessed at P=0.01 level of chi-square',quote=F);
print(table(OE.final[,'outlier.01']));
print(' ',quote=F)
print('Predictions complete', quote=F);
#return a list with OE.final and the matrix of predicted probs as elements;
list(OE.scores=OE.final,Capture.Probs=site.pred.dfa,Group.Occurrence.Probs=grpprobs); #return data frame as final object;
OE.final<-OE.final[row.names(OE.final)==row.names(temp.pa),]
}; #end of function;



