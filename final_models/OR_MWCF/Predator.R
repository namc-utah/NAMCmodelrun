############################
# From JVS original code:
# The 2 files (bugs and predictors) should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
# Notes on format --
#   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
#   B) Predictor data set -- Must include columns with the same names, units, etc.,
#        as the model's predictor variables. All other columns will be ignored;
#        Column order does not matter;
#        Predictions and calculations of O/E will be made only for those samples that have;
#        complete data for all model predictors.;
#   C)  Sample-by-taxa matrix.  Sample ID's (row names) must match those of predictor data.
#        Model predictions will include only those taxa (columns) in this file that are also
#        in the calibration data (bugcal.pa);
#        To see a list of the calibration taxa, do the following:
#                                                         names(bugcal.pa)[colSums(bugcal.pa)>0];
#  Code modified by TWA (January 2019) in order to utilize NAMC file structure and output needs
##########;

OR_MWCF_model<- function(bugsOTU_matrix,test_preds,rf_model){
#make predictions for test data;
OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,
                                   prednew=test_preds,bugnew=bugsOTU_matrix,Pc=0.5);
return(OE)
}


#############

# 5.2  WCCP model

#############

OR_WCCP_model<- function (bugsOTU_matrix,test_preds,rf_model){
#make predictions for test data;
OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,
                                   prednew=test_preds.F,bugnew=bugsOTU_matrix,Pc=0.5);
return(OE)
}
