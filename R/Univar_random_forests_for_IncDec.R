#Traits as a function of the environment
trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')
names(trait_table)

traits=names(trait_table)[2:ncol(trait_table)]
#How well can we predict a trait using environmental factors?

#PA of a trait:
PA_sites_prob=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//PA_all_traits_probabilistic.csv')
PA_sites_prob=PA_sites_prob[PA_sites_prob$sampleId %in% sites_w_no_bugs==F,]
PA_sites_prob=PA_sites_prob[PA_sites_prob$sampleId!='sampleId',]
PA_sites_prob[traits]=lapply(PA_sites_prob[traits], as.factor)
library(randomForest)
library(ggplot2)

# predictors (everything except response)
response_varsPA=names(PA_sites_prob)[1:25]
predictor_vars <- names(PA_sites_prob)[c(27:40,43,47)]
#oob_results_list <- list()
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'Prob_PA')
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = PA_sites_prob,
    importance = TRUE
  )

  # Extract the final OOB error
  # final_oob <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
  #
  # # Store the current iteration as a small 1-row data frame in the list
  # oob_results_list[[resp]] <- data.frame(
  #   Response_Variable = resp,
  #   OOB_Error = final_oob,
  #   stringsAsFactors = FALSE
  # )

  # print(paste(resp, final_oob, sep=' '))
# }

# 3. Combine all rows into a single full data frame after the loop
# final_oob_df <- do.call(rbind, oob_results_list)
#
# # Reset row names for clean indexing
# rownames(final_oob_df) <- NULL
# print(paste(resp,rf_model$err.rate[nrow(rf_model$err.rate), "OOB"],sep=' '))


  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = PA_sites_prob,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          strip.text = element_text(size = 24),
          strip.background = element_blank(),
          legend.text = element_text(size=20)
        )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
  }




Pct_sites_prob=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
Pct_sites_prob=Pct_sites_prob[Pct_sites_prob$sampleId %in% sites_w_no_bugs==F,]
Pct_sites_prob=Pct_sites_prob[Pct_sites_prob$sampleId!='sampleId',]

library(randomForest)
library(ggplot2)

# predictors (everything except response)
response_varsPA=names(Pct_sites_prob)[1:25]
predictor_vars <- names(Pct_sites_prob)[c(27:40, 43, 47)]
reg_results_list <- list()
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  #base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Prob_Pct")
  #outdir <- file.path(base_dir, resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = Pct_sites_prob,
    importance = TRUE
  )
#   # Extract final R-squared and calculate final RMSE
#   final_rsq  <- rf_model$rsq[length(rf_model$rsq)]
#   final_rmse <- sqrt(rf_model$mse[length(rf_model$mse)])
#
#   # Store the metrics in a 1-row data frame
#   reg_results_list[[resp]] <- data.frame(
#     Response_Variable = resp,
#     R_Squared         = final_rsq,
#     RMSE              = final_rmse,
#     stringsAsFactors  = FALSE
#   )
#
#   # Print console update
#   print(paste(resp, "R2:", round(final_rsq, 4), "RMSE:", round(final_rmse, 4), sep='  '))
# }

# 3. Combine all rows into a single data frame after the loop
#final_reg_df <- do.call(rbind, reg_results_list)

# Reset row names for clean indexing
#rownames(final_reg_df) <- NULL

  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = PA_sites_prob,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}
#getting indivs into wide format for % traits
Pct_indiv_sites_prob=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Prob_abundances.csv')
Pct_indiv_sites_prob=Pct_indiv_sites_prob[Pct_indiv_sites_prob$sample_id %in% c(sites_w_no_bugs,failed_sites$sampleId)==F,]

trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')
#Pct_indiv_sites_prob=Pct_indiv_sites_prob[,names(Pct_indiv_sites_prob)!='sample_Id']
Pct_indiv_sites_prob=Pct_indiv_sites_prob[Pct_indiv_sites_prob$OTU %in% trait_table$OTU,]

trait_table=trait_table[!duplicated(trait_table$OTU),]
#join the abundances to the trait table by OTU
abun_traits=Pct_indiv_sites_prob %>%
  left_join(trait_table,by='OTU')
#keep OTU and the trait names so we can summarize by sampleID
#and OTU to get the abundances and trait states we need.
trait_cols_abun=names(trait_table[,2:ncol(trait_table)])
#for each sample,
#find the taxa whose trait state is 1, sum them divide bthe total,
#and multiply by 100.
pct_indivs_w_Trait_summary <- abun_traits %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    across(
      all_of(trait_cols_abun),
      ~ 100 * sum(split_count * .x, na.rm = TRUE) / sum(split_count, na.rm = TRUE),
      .names = "Pct_{.col}"
    )
  )
#environmental_for_indivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
#environmental_for_indivs=
#environmental_for_indivs=environmental_for_indivs[,names(environmental_for_indivs)%in% c(predictor_vars,'sampleId')]
#names(environmental_for_indivs)[names(environmental_for_indivs)=='sampleId']<-'sample_id'
#pct_indivs_w_Trait_summary=plyr::join(pct_indivs_w_Trait_summary,environmental_for_indivs)

library(randomForest)
library(ggplot2)
#names(pct_indivs_w_Trait_summary)[2:26]<-gsub("^Pct_", "", names(pct_indivs_w_Trait_summary)[2:26])
#pct_indivs_w_Trait_summary=pct_indivs_w_Trait_summary[,names(pct_indivs_w_Trait_summary)!='sample_id']
# predictors (everything except response)
response_varsPA=names(pct_indivs_w_Trait_summary)[1:25]
predictor_vars <- names(pct_indivs_w_Trait_summary)[c(26:39,42,46)]

pct_indivs_w_Trait_summary

for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Prob_PctIndiv")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = pct_indivs_w_Trait_summary,
    importance = TRUE
  )

  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = pct_indivs_w_Trait_summary,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}


#Reference PA

PA_sites_R=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//PA_all_traits_ref.csv')
#PA_sites_R[!complete.cases(PA_sites_R),]
PA_sites_R=PA_sites_R[PA_sites_R$sampleId %in% c('HAWK-33','WE-131','WE-209')==F,]
PA_sites_R=PA_sites_R[,names(PA_sites_R)!='sampleId',]
PA_sites_R[traits]=lapply(PA_sites_R[traits], as.factor)




library(randomForest)
library(ggplot2)

# predictors (everything except response)
response_varsPA=names(PA_sites_R)[1:25]
predictor_vars <- names(PA_sites_R)[c(26:39,42,46)]

for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'Ref_PA_trait')
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = PA_sites_R,
    importance = TRUE
  )
  # Extract the final OOB error
   final_oob <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
  #
  # # Store the current iteration as a small 1-row data frame in the list
   oob_results_list[[resp]] <- data.frame(
     Response_Variable = resp,
     OOB_Error = final_oob,
     stringsAsFactors = FALSE
   )

   #print(paste(resp, final_oob, sep=' '))
   #}

  # 3. Combine all rows into a single full data frame after the loop
  # final_oob_df <- do.call(rbind, oob_results_list)
  #
  # # Reset row names for clean indexing
 #  rownames(final_oob_df) <- NULL
  # print(paste(resp,rf_model$err.rate[nrow(rf_model$err.rate), "OOB"],sep=' '))
#}
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = PA_sites_R,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}


Pct_sites_R=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_ref.csv')
Pct_sites_R=Pct_sites_R[Pct_sites_R$sampleId %in% c('HAWK-33','WE-131','WE-209')==F,]
Pct_sites_R=Pct_sites_R[,names(Pct_sites_R)!='sampleId',]

library(randomForest)
library(ggplot2)

# predictors (everything except response)

response_varsPA=names(Pct_sites_R)[1:25]
predictor_vars <- names(Pct_sites_R)[c(26:39,42,46)]

names(PA_sites_R)

for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Ref_Pct")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = Pct_sites_R,
    importance = TRUE
  )
  # # #   # Extract final R-squared and calculate final RMSE
  #    final_rsq  <- rf_model$rsq[length(rf_model$rsq)]
  #    final_rmse <- sqrt(rf_model$mse[length(rf_model$mse)])
  # # #
  # # #   # Store the metrics in a 1-row data frame
  #    reg_results_list[[resp]] <- data.frame(
  #      Response_Variable = resp,
  #      R_Squared         = final_rsq,
  #     RMSE              = final_rmse,
  #      stringsAsFactors  = FALSE
  #    )
  # # #
  # # #   # Print console update
  #    print(paste(resp, "R2:", round(final_rsq, 4), "RMSE:", round(final_rmse, 4), sep='  '))
  #  }
  # #
  # # # 3. Combine all rows into a single data frame after the loop
  # final_reg_df <- do.call(rbind, reg_results_list)
  # #
  # # # Reset row names for clean indexing
  # rownames(final_reg_df) <- NULL
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = Pct_sites_R,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}
 ## % I and %D at sites (prob only has indivs)

P_comm=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
row.names(P_comm)<-P_comm$X;P_comm=P_comm[,-c(1,2)]
P_comm=P_comm[row.names(P_comm) %in% c(failed_sites$sampleId, sites_w_no_bugs)==F,]
P_comm$Site=row.names(P_comm)
P_comm_long=P_comm %>%
  pivot_longer(
    cols=-Site,
    names_to = "taxon",
    values_to = "Present"
  ) %>%
  filter(Present == 1)

P_comm_long <- P_comm_long %>%
  dplyr::left_join(
    P_Results[, c("taxon", "Regional_Response")],
    by = "taxon"
  )

P_comm_long <- P_comm_long %>%
  dplyr::mutate(
    is_I = Regional_Response == "Increaser",
    is_D = Regional_Response == "Decreaser"
  )

##
PA_P_response <- P_comm_long %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(

    TotalTaxa = dplyr::n(),

    I_taxa = sum(is_I, na.rm = TRUE),
    D_taxa = sum(is_D, na.rm = TRUE),

    Percent_I = 100 * I_taxa / TotalTaxa,
    Percent_D = 100 * D_taxa / TotalTaxa,

    Presence_I = as.integer(I_taxa > 0),
    Presence_D = as.integer(D_taxa > 0),

    .groups = "drop"
  )
##


PA_P_response

Prob_abun=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Prob_abundances_OG.csv')
#Prob_abun=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Prob_abundances.csv")

names(Prob_abun)[2]<-'taxon'
names(Prob_abun)[3]<-'Site'
Prob_abun=Prob_abun[,-1]
Prob_abun=Prob_abun[Prob_abun$Site %in% c(failed_sites$sampleId,sites_w_no_bugs)==F,]
Prob_abun_wresp=Prob_abun %>%
  dplyr::left_join(
    P_Results[, c("taxon", "Regional_Response")],
    by = "taxon"
  )
Prob_abun_wresp = Prob_abun_wresp %>%
  dplyr::mutate(
    is_I = Regional_Response == "Increaser",
    is_D = Regional_Response == "Decreaser"
  )



##
Prob_abun_wresp2=Prob_abun_wresp %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(

    TotalIndividuals = sum(split_count, na.rm = TRUE),

    I_individuals = sum(split_count * is_I, na.rm = TRUE),
    D_individuals = sum(split_count * is_D, na.rm = TRUE),

    Percent_I = 100 * I_individuals / TotalIndividuals,
    Percent_D = 100 * D_individuals / TotalIndividuals,

    Presence_I = as.integer(I_individuals > 0),
    Presence_D = as.integer(D_individuals > 0),

    .groups = "drop"
  )
##
## % indivs that are increasers
environmental_for_indivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
environmental_for_indivs=environmental_for_indivs[,names(environmental_for_indivs)%in% c(predictor_vars,'sampleId')]
names(environmental_for_indivs)[1]<-'Site'

PA_P_response_wtraits=as.data.frame(plyr::join(PA_P_response,environmental_for_indivs))
PA_P_response_wtraits=PA_P_response_wtraits[PA_P_response_wtraits$Site %in% c(failed_sites, sites_w_no_bugs)==F,]
I_P_response=PA_P_response_wtraits[,c('Percent_I',predictor_vars)]

predictor_vars=names(I_P_response)[2:ncol(I_P_response)]
response_varsPA=names(I_P_response)[1]
predictor_vars=c('Precip8110',
                 'MAST_mean08091314',
                 'WsAreaSqKm',
                 'MSST_mean08091314',
                 'LONG')
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Increase_top5")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = I_P_response,
    importance = TRUE
  )

 # rf_model$importance
  sqrt(rf_model$mse[length(rf_model$mse)])
  rf_model$rsq[length(rf_model$rsq)]
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = I_P_response,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}



environmental_for_indivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
environmental_for_indivs=environmental_for_indivs[,names(environmental_for_indivs)%in% c(predictor_vars,'sampleId')]
names(environmental_for_indivs)[1]<-'Site'

abun_P_response_wtraits=as.data.frame(plyr::join(Prob_abun_wresp2,environmental_for_indivs))
abun_P_response_wtraits=abun_P_response_wtraits[abun_P_response_wtraits$Site %in% c(failed_sites, sites_w_no_bugs)==F,]
I_abun_response=abun_P_response_wtraits[,c('Percent_I',predictor_vars)]

predictor_vars=names(I_P_response)[2:ncol(I_P_response)]
response_varsPA=names(I_P_response)[1]
predictor_vars=c('LONG','WsAreaSqKm','MSST_mean08091314',
                 'Precip8110','MWST_mean08091314')
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Indiv_Increase_top5")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = I_abun_response,
    importance = TRUE
  )
  #varImpPlot(rf_model)
   sqrt(rf_model$mse[length(rf_model$mse)])
   rf_model$rsq[length(rf_model$rsq)]
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = I_abun_response,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}


## % Decreasers
PA_P_response_wtraits=as.data.frame(plyr::join(PA_P_response,environmental_for_indivs))
PA_P_response_wtraits=PA_P_response_wtraits[PA_P_response_wtraits$Site %in% c(failed_sites, sites_w_no_bugs)==F,]
D_P_response=PA_P_response_wtraits[,c('Percent_D',predictor_vars)]

predictor_vars=names(D_P_response)[2:ncol(D_P_response)]
response_varsPA=names(D_P_response)[1]
predictor_vars=c('Precip8110','MAST_mean08091314',
                 'ElevCat','LONG','WsAreaSqKm')
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Decrease_top5")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = D_P_response,
    importance = TRUE
  )
  #varImpPlot(rf_model)
   #sqrt(rf_model$mse[length(rf_model$mse)])
   #rf_model$rsq[length(rf_model$rsq)]
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = D_P_response,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}


abun_P_response_wtraits=as.data.frame(plyr::join(Prob_abun_wresp2,environmental_for_indivs))
abun_P_response_wtraits=abun_P_response_wtraits[abun_P_response_wtraits$Site %in% c(failed_sites, sites_w_no_bugs)==F,]
D_abun_response=abun_P_response_wtraits[,c('Percent_D',predictor_vars)]

predictor_vars=names(D_abun_response)[2:ncol(D_abun_response)]
response_varsPA=names(D_abun_response)[1]
predictor_vars=c('Precip8110','Tmean8110Ws','MSST_mean08091314',
                 'LONG','MAST_mean08091314')
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Indiv_Decrease_top5")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))

  rf_model <- randomForest(
    form,
    data = D_abun_response,
    importance = TRUE
  )
  #varImpPlot(rf_model)
  sqrt(rf_model$mse[length(rf_model$mse)])
  rf_model$rsq[length(rf_model$rsq)]
  # -----------------------------
  # 3. Variable importance plot
  # -----------------------------
  vi <- data.frame(
    variable = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]
  )

  vi_plot <- ggplot(vi, aes(x = reorder(variable, importance), y = importance)) +
    geom_point(size=5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          panel.grid.major.y = element_line(color = "gray5", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(
    filename = file.path(outdir, "VarImp.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )
  predictor_vars <- rownames(importance(rf_model))
  #savp(10,8,paste0(outdir,'//VarImp.png'))

  # -----------------------------
  # 4. Partial dependence plots
  # -----------------------------
  for (v in predictor_vars) {

    pd <- randomForest::partialPlot(
      rf_model,
      pred.data = D_abun_response,
      x.var = paste(v),
      plot = FALSE
    )

    df <- data.frame(x = pd$x, y = pd$y)

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = v,
        x = v,
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size=20)
      )
    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, ".png")),
      plot = p,
      width = 10,
      height = 8
    )
    #savp(10,8,paste0(outdir,'//PDP.png'))

  }
}


#P/A of a trait state at sites as fxn of the environment
#top 5 only!

# -----------------------------
# Store OOB results
# -----------------------------

response_varsPA=names(PA_sites_prob)[1:25]
predictor_vars <- names(PA_sites_prob)[c(27:40,43,47)]
oob_results_list <- list()

for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results", "Prob_PA_top5s")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = TRUE)

  # -----------------------------
  # 2. Fit full RF model
  # -----------------------------
  full_form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + "))
  )

  rf_model_full <- randomForest(
    full_form,
    data = PA_sites_prob,
    importance = TRUE
  )

  # -----------------------------
  # 3. Extract top 5 predictors
  # -----------------------------
  vi_full <- data.frame(
    variable = rownames(importance(rf_model_full)),
    importance = importance(rf_model_full)[, "MeanDecreaseGini"]
  )

  vi_full <- vi_full[order(vi_full$importance, decreasing = TRUE), ]

  top5_predictors <- vi_full$variable[1:5]

  print(paste(resp, "Top 5 predictors:"))
  print(top5_predictors)

  # -----------------------------
  # 4. Fit reduced RF model
  # -----------------------------
  reduced_form <- as.formula(
    paste(resp, "~", paste(top5_predictors, collapse = " + "))
  )

  rf_model_top5 <- randomForest(
    reduced_form,
    data = PA_sites_prob,
    importance = TRUE
  )

  # -----------------------------
  # 5. Extract OOB error
  # -----------------------------
  final_oob <- rf_model_top5$err.rate[
    nrow(rf_model_top5$err.rate),
    "OOB"
  ]

  oob_results_list[[resp]] <- data.frame(
    Response_Variable = resp,
    OOB_Error = final_oob,
    Top5_Predictors = paste(top5_predictors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  print(paste(resp, "OOB =", round(final_oob, 4)))

  # -----------------------------
  # 6. Variable importance plot
  # -----------------------------
  vi_top5 <- data.frame(
    variable = rownames(importance(rf_model_top5)),
    importance = importance(rf_model_top5)[, "MeanDecreaseGini"]
  )

  vi_plot <- ggplot(
    vi_top5,
    aes(x = reorder(variable, importance), y = importance)
  ) +
    geom_point(size = 5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Top 5 Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      panel.grid.major.y = element_line(
        color = "gray50",
        linetype = "dotted"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(
    filename = file.path(outdir, "VarImp_top5.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )

  # -----------------------------
  # 7. PDPs for TOP 5 MODEL ONLY
  # -----------------------------
  thetop5=row.names(rf_model_top5$importance)
  for (v in 1:length(thetop5)) {

    pd <- randomForest::partialPlot(
      rf_model_top5,
      pred.data = PA_sites_prob,
      x.var = thetop5[v],
      plot = F
    )

    df <- data.frame(
      x = pd$x,
      y = pd$y
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = thetop5[v],
        x = thetop5[v],
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size = 20)
      )

    ggsave(
      filename = file.path(outdir, paste0("PDP_", thetop5[v], "_top5.png")),
      plot = p,
      width = 10,
      height = 8
    )
  }
}

# -----------------------------
# 8. Combine OOB results
# -----------------------------
final_oob_df <- do.call(rbind, oob_results_list)

rownames(final_oob_df) <- NULL


response_varsPA=names(PA_sites_R)[1:25]
predictor_vars <- names(PA_sites_R)[c(26:39,42,46)]
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results", "Ref_PA_top5s")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = TRUE)

  # -----------------------------
  # 2. Fit full RF model
  # -----------------------------
  full_form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + "))
  )

  rf_model_full <- randomForest(
    full_form,
    data = PA_sites_R,
    importance = TRUE
  )

  # -----------------------------
  # 3. Extract top 5 predictors
  # -----------------------------
  vi_full <- data.frame(
    variable = rownames(importance(rf_model_full)),
    importance = importance(rf_model_full)[, "MeanDecreaseGini"]
  )

  vi_full <- vi_full[order(vi_full$importance, decreasing = TRUE), ]

  top5_predictors <- vi_full$variable[1:5]

  print(paste(resp, "Top 5 predictors:"))
  print(top5_predictors)

  # -----------------------------
  # 4. Fit reduced RF model
  # -----------------------------
  reduced_form <- as.formula(
    paste(resp, "~", paste(top5_predictors, collapse = " + "))
  )

  rf_model_top5 <- randomForest(
    reduced_form,
    data = PA_sites_R,
    importance = TRUE
  )

  # -----------------------------
  # 5. Extract OOB error
  # -----------------------------
  final_oob <- rf_model_top5$err.rate[
    nrow(rf_model_top5$err.rate),
    "OOB"
  ]

  oob_results_list[[resp]] <- data.frame(
    Response_Variable = resp,
    OOB_Error = final_oob,
    Top5_Predictors = paste(top5_predictors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  print(paste(resp, "OOB =", round(final_oob, 4)))

  # -----------------------------
  # 6. Variable importance plot
  # -----------------------------
  vi_top5 <- data.frame(
    variable = rownames(importance(rf_model_top5)),
    importance = importance(rf_model_top5)[, "MeanDecreaseGini"]
  )

  vi_plot <- ggplot(
    vi_top5,
    aes(x = reorder(variable, importance), y = importance)
  ) +
    geom_point(size = 5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Top 5 Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      panel.grid.major.y = element_line(
        color = "gray50",
        linetype = "dotted"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(
    filename = file.path(outdir, "VarImp_top5.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )

  # -----------------------------
  # 7. PDPs for TOP 5 MODEL ONLY
  # -----------------------------
  thetop5=row.names(rf_model_top5$importance)
  for (v in 1:length(thetop5)) {

    pd <- randomForest::partialPlot(
      rf_model_top5,
      pred.data = PA_sites_R,
      x.var = thetop5[v],
      plot = F
    )

    df <- data.frame(
      x = pd$x,
      y = pd$y
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = thetop5[v],
        x = thetop5[v],
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size = 20)
      )

    ggsave(
      filename = file.path(outdir, paste0("PDP_", thetop5[v], "_top5.png")),
      plot = p,
      width = 10,
      height = 8
    )
  }
}

# -----------------------------
# 8. Combine OOB results
# -----------------------------
final_oob_df <- do.call(rbind, oob_results_list)

rownames(final_oob_df) <- NULL

clipr::write_clip(final_oob_df)


##now for regression top5s


# -----------------------------
# Store model metrics
# -----------------------------
model_results_list <- list()
response_varsPA=names(Pct_sites_prob)[1:25]
predictor_vars <- names(Pct_sites_prob)[c(27:40, 43, 47)]
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(
    path.expand("~"),
    "New_for_SFS",
    "RF_results",
    "Prob_Pct_top5s"
  )

  outdir <- file.path(base_dir, resp)

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = TRUE
  )

  # -----------------------------
  # 2. Fit full RF model
  # -----------------------------
  full_form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + "))
  )

  rf_model_full <- randomForest(
    full_form,
    data = Pct_sites_prob,
    importance = TRUE
  )

  # -----------------------------
  # 3. Extract top 5 predictors
  # -----------------------------
  vi_full <- data.frame(
    variable = rownames(importance(rf_model_full)),
    importance = importance(rf_model_full)[, "%IncMSE"]
  )

  vi_full <- vi_full[
    order(vi_full$importance, decreasing = TRUE),
  ]

  top5_predictors <- vi_full$variable[1:5]

  print(paste(resp, "Top 5 predictors:"))
  print(top5_predictors)

  # -----------------------------
  # 4. Fit reduced RF model
  # -----------------------------
  reduced_form <- as.formula(
    paste(resp, "~", paste(top5_predictors, collapse = " + "))
  )

  rf_model_top5 <- randomForest(
    reduced_form,
    data = Pct_sites_prob,
    importance = TRUE
  )

  # -----------------------------
  # 5. Extract RMSE + pseudo R2
  # -----------------------------

  # Final OOB MSE
  final_mse <- tail(rf_model_top5$mse, 1)

  # RMSE
  final_rmse <- sqrt(final_mse)

  # Pseudo R2
  final_r2 <- tail(rf_model_top5$rsq, 1)

  model_results_list[[resp]] <- data.frame(
    Response_Variable = resp,
    RMSE = final_rmse,
    Pseudo_R2 = final_r2,
    Top5_Predictors = paste(top5_predictors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  print(
    paste(
      resp,
      "RMSE =",
      round(final_rmse, 4),
      "| Pseudo R2 =",
      round(final_r2, 4)
    )
  )

  # -----------------------------
  # 6. Variable importance plot
  # -----------------------------
  vi_top5 <- data.frame(
    variable = rownames(importance(rf_model_top5)),
    importance = importance(rf_model_top5)[, "%IncMSE"]
  )

  vi_plot <- ggplot(
    vi_top5,
    aes(x = reorder(variable, importance), y = importance)
  ) +
    geom_point(size = 5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Top 5 Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      panel.grid.major.y = element_line(
        color = "gray50",
        linetype = "dotted"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(
    filename = file.path(outdir, "VarImp_top5.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )

  # -----------------------------
  # 7. PDPs for TOP 5 MODEL ONLY
  # -----------------------------
  thetop5 <- rownames(rf_model_top5$importance)

  for (v in 1:length(thetop5)) {

    pd <- randomForest::partialPlot(
      rf_model_top5,
      pred.data = Pct_sites_prob,
      x.var = thetop5[v],
      plot = FALSE
    )

    df <- data.frame(
      x = pd$x,
      y = pd$y
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = thetop5[v],
        x = thetop5[v],
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size = 20)
      )

    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, "_top5.png")),
      plot = p,
      width = 10,
      height = 8
    )
  }
}

# -----------------------------
# 8. Combine model metrics
# -----------------------------
final_model_metrics_df <- do.call(
  rbind,
  model_results_list
)

rownames(final_model_metrics_df) <- NULL

print(final_model_metrics_df)

clipr::write_clip(final_model_metrics_df
                  )


model_results_list <- list()
response_varsPA=names(Pct_sites_R)[1:25]
predictor_vars <- names(Pct_sites_R)[c(26:39, 42, 46)]
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(
    path.expand("~"),
    "New_for_SFS",
    "RF_results",
    "Ref_Pct_top5s"
  )

  outdir <- file.path(base_dir, resp)

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = TRUE
  )

  # -----------------------------
  # 2. Fit full RF model
  # -----------------------------
  full_form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + "))
  )

  rf_model_full <- randomForest(
    full_form,
    data = Pct_sites_R,
    importance = TRUE
  )

  # -----------------------------
  # 3. Extract top 5 predictors
  # -----------------------------
  vi_full <- data.frame(
    variable = rownames(importance(rf_model_full)),
    importance = importance(rf_model_full)[, "%IncMSE"]
  )

  vi_full <- vi_full[
    order(vi_full$importance, decreasing = TRUE),
  ]

  top5_predictors <- vi_full$variable[1:5]

  print(paste(resp, "Top 5 predictors:"))
  print(top5_predictors)

  # -----------------------------
  # 4. Fit reduced RF model
  # -----------------------------
  reduced_form <- as.formula(
    paste(resp, "~", paste(top5_predictors, collapse = " + "))
  )

  rf_model_top5 <- randomForest(
    reduced_form,
    data = Pct_sites_R,
    importance = TRUE
  )

  # -----------------------------
  # 5. Extract RMSE + pseudo R2
  # -----------------------------

  # Final OOB MSE
  final_mse <- tail(rf_model_top5$mse, 1)

  # RMSE
  final_rmse <- sqrt(final_mse)

  # Pseudo R2
  final_r2 <- tail(rf_model_top5$rsq, 1)

  model_results_list[[resp]] <- data.frame(
    Response_Variable = resp,
    RMSE = final_rmse,
    Pseudo_R2 = final_r2,
    Top5_Predictors = paste(top5_predictors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  print(
    paste(
      resp,
      "RMSE =",
      round(final_rmse, 4),
      "| Pseudo R2 =",
      round(final_r2, 4)
    )
  )

  # -----------------------------
  # 6. Variable importance plot
  # -----------------------------
  vi_top5 <- data.frame(
    variable = rownames(importance(rf_model_top5)),
    importance = importance(rf_model_top5)[, "%IncMSE"]
  )

  vi_plot <- ggplot(
    vi_top5,
    aes(x = reorder(variable, importance), y = importance)
  ) +
    geom_point(size = 5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Top 5 Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      panel.grid.major.y = element_line(
        color = "gray50",
        linetype = "dotted"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(
    filename = file.path(outdir, "VarImp_top5.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )

  # -----------------------------
  # 7. PDPs for TOP 5 MODEL ONLY
  # -----------------------------
  thetop5 <- rownames(rf_model_top5$importance)

  for (v in 1:length(thetop5)) {

    pd <- randomForest::partialPlot(
      rf_model_top5,
      pred.data = Pct_sites_R,
      x.var = thetop5[v],
      plot = FALSE
    )

    df <- data.frame(
      x = pd$x,
      y = pd$y
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = thetop5[v],
        x = thetop5[v],
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size = 20)
      )

    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, "_top5.png")),
      plot = p,
      width = 10,
      height = 8
    )
  }
}

# -----------------------------
# 8. Combine model metrics
# -----------------------------
final_model_metrics_df <- do.call(
  rbind,
  model_results_list
)

rownames(final_model_metrics_df) <- NULL

print(final_model_metrics_df)

clipr::write_clip(final_model_metrics_df
)



#Pct indivs with a trait state top5
Pct_indivis_w_env=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Prob_Pct_indivs_wtrait_and_envs.csv')

model_results_list <- list()
response_varsPA=names(Pct_indivis_w_env)[1:25]
predictor_vars <- names(Pct_sites_R)[c(26:39, 42, 46)]
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(
    path.expand("~"),
    "New_for_SFS",
    "RF_results",
    "Pct_Indivs_w_Trait_State_top5"
  )

  outdir <- file.path(base_dir, resp)

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = TRUE
  )

  # -----------------------------
  # 2. Fit full RF model
  # -----------------------------
  full_form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + "))
  )

  rf_model_full <- randomForest(
    full_form,
    data = Pct_indivis_w_env,
    importance = TRUE
  )

  # -----------------------------
  # 3. Extract top 5 predictors
  # -----------------------------
  vi_full <- data.frame(
    variable = rownames(importance(rf_model_full)),
    importance = importance(rf_model_full)[, "%IncMSE"]
  )

  vi_full <- vi_full[
    order(vi_full$importance, decreasing = TRUE),
  ]

  top5_predictors <- vi_full$variable[1:5]

  print(paste(resp, "Top 5 predictors:"))
  print(top5_predictors)

  # -----------------------------
  # 4. Fit reduced RF model
  # -----------------------------
  reduced_form <- as.formula(
    paste(resp, "~", paste(top5_predictors, collapse = " + "))
  )

  rf_model_top5 <- randomForest(
    reduced_form,
    data = Pct_indivis_w_env,
    importance = TRUE
  )

  # -----------------------------
  # 5. Extract RMSE + pseudo R2
  # -----------------------------

  # Final OOB MSE
  final_mse <- tail(rf_model_top5$mse, 1)

  # RMSE
  final_rmse <- sqrt(final_mse)

  # Pseudo R2
  final_r2 <- tail(rf_model_top5$rsq, 1)

  model_results_list[[resp]] <- data.frame(
    Response_Variable = resp,
    RMSE = final_rmse,
    Pseudo_R2 = final_r2,
    Top5_Predictors = paste(top5_predictors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  print(
    paste(
      resp,
      "RMSE =",
      round(final_rmse, 4),
      "| Pseudo R2 =",
      round(final_r2, 4)
    )
  )

  # -----------------------------
  # 6. Variable importance plot
  # -----------------------------
  vi_top5 <- data.frame(
    variable = rownames(importance(rf_model_top5)),
    importance = importance(rf_model_top5)[, "%IncMSE"]
  )

  vi_plot <- ggplot(
    vi_top5,
    aes(x = reorder(variable, importance), y = importance)
  ) +
    geom_point(size = 5) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      title = paste("Top 5 Variable Importance -", resp),
      x = NULL,
      y = "Mean Decrease Gini"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      panel.grid.major.y = element_line(
        color = "gray50",
        linetype = "dotted"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(
    filename = file.path(outdir, "VarImp_top5.png"),
    plot = vi_plot,
    width = 10,
    height = 8
  )

  # -----------------------------
  # 7. PDPs for TOP 5 MODEL ONLY
  # -----------------------------
  thetop5 <- rownames(rf_model_top5$importance)

  for (v in 1:length(thetop5)) {

    pd <- randomForest::partialPlot(
      rf_model_top5,
      pred.data = Pct_sites_R,
      x.var = thetop5[v],
      plot = FALSE
    )

    df <- data.frame(
      x = pd$x,
      y = pd$y
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) +
      theme_classic(base_size = 16) +
      labs(
        title = thetop5[v],
        x = thetop5[v],
        y = "Partial Effect"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        legend.text = element_text(size = 20)
      )

    ggsave(
      filename = file.path(outdir, paste0("PDP_", v, "_top5.png")),
      plot = p,
      width = 10,
      height = 8
    )
  }
}

# -----------------------------
# 8. Combine model metrics
# -----------------------------
final_model_metrics_df <- do.call(
  rbind,
  model_results_list
)

rownames(final_model_metrics_df) <- NULL

print(final_model_metrics_df)

clipr::write_clip(final_model_metrics_df
)

