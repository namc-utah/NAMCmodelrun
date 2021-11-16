# NAMCmodelrun

This repository contains the R code and accompanying files for NAMC models. 

models are of three types
1. OE
2. MMI
3. WQ

models have 3 main inputs needed and seperate functions depending on model type (model function R file)
1. new bug data (OE and MMIs only)
2. new predictors 
3. reference data and random forest models used to build the model

process sample models and process box model functions 
1. get needed predictors stored in the database
2. get bug data stored in the database by calling bug functions in the bug functions R file
3. call model functions from the model function R file
4. call the model appicability function
5. save the model results in the database
