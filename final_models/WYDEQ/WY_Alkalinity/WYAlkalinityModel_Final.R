##This is the streamlined version of the alkalinity model for Wyoming.  Use this to predict alkalinity values for inclusion in the predictor file used in the Wyoming OE model.  

##The input file is generated using the 'pts2scat_20180208.r' code, which pulls all StreamCAt data from the repository and creates the file that can be fed into this code.  Once the StreamCat file is generated, place it into Z:\buglab\OE_Modeling\NAMC_Supported_OEmodels\Wyoming\Alkalinity\StreamCat_input folder.
##The pts2cat script has been updated (May 2019, TWA)

setwd("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\Wyoming\\Alkalinity\\")
load(paste0("AlkalinityModel.RData"))
library(randomForest)

# Now, bring in the validation data:

test.data = read.csv(paste0(".\\StreamCat_input\\AlkPred2110.csv"))
row.names(test.data) = test.data$sampleid

test.data11 = test.data[,c("COMID","CompStrgthWs","ElevWs","PctConif2006Ws","PctShrb2011Ws","PctSilicicWs","PermWs","Pestic97Ws","Precip08Ws","SandWs","Tmax8110Ws","Tmean08Ws")]
 

# Predict alkalinity for the test points:
predict.alkalinity.from.11preds.for.val=predict(m1,newdata=test.data11,type="response")

predict.alkalinity.from.11preds.for.val

predicted.alkalinity.val=as.data.frame(predict.alkalinity.from.11preds.for.val)
dim(predicted.alkalinity.val)
head(predicted.alkalinity.val)

write.csv(predicted.alkalinity.val,file=".\\WYAlkPred2110_ModelPredictions_trip.csv")
