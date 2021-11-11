   # User must supply a sample-by-taxa matrix of taxon presence/absence (coded as 1 or 0), for all new samples;
   # User must also supply a corresponding file of predictor data for those same samples;
   # These 2 files should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
   # Notes on format --
   #   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
   #   B) Predictor data set -- Must include columns with the same names, units, etc.,
   #        as the model's predictor variables. All other columns will be ignored;
   #        Column order does not matter;
   #        Predictions and calculations of O/E will be made only for those samples that have;
   #        complete data for all model predictors.;
   #   C)  Sample-by-taxa matrix.  Sample ID's (row names) must match those of predictor data.
   #        OTU names (column names) must include all taxa for which the model can make predictions. Other taxa (columns) are ignored.
   #        The model can make predictions for all taxa in the calibration data (bugcal.pa) that were present at 1 or more calibration sites;
   #        To see a list of these taxa, do the following:

PIBO_model<-function (bugsOTU_matrix,test_preds,rf_model){
### what is this associated with??
         names(bugall)[colSums(bugall)>0];
#########;

#### Trip is this needed???         
bug.otu <- read.csv("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\OTUbugs1.csv",header=T);
# merges test bug data with a null set of all OTU bug names in order for VanSickle code to work properly.
bugsOTU_matrix.pa <- merge(x=bug.otu, y=bugsOTU_matrix, all=TRUE)         
# converts na data produced from columns above into 0s;
bugsOTU_matrix.pa[is.na(bugsOTU_matrix.pa)]<- 0 

# Reorder sampleid in ascending order
bugsOTU_matrix.pa = bugsOTU_matrix.pa[order(bugsOTU_matrix.pa$sampleId), ];
test_preds = test_preds[order(test_preds$sampleId), ];

# makes predictions for test data;
OE<-model.predict.RanFor(bugall, grps.final, preds.final, ranfor.mod=rf.mod, prednew=test_preds, bugnew.pa=bugsOTU_matrix.pa, Pc=0.5);

###is this needed????
# Append Sampleid to results
OE$OE.scores = data.frame(bugsOTU_matrix.pa$sampleId, OE$OE.scores)
names(OE$OE.scores) = c('sampleId',names(OE$OE.scores)[2:ncol(OE$OE.scores)] )
OE$Capture.Probs = data.frame(bugsOTU_matrix.pa$sampleId, OE$Capture.Probs)
names(OE$Capture.Probs) = c('sampleId',names(OE$Capture.Probs)[2:ncol(OE$Capture.Probs)] )

return(OE)
}
