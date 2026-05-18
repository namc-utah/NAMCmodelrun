#site status as function of traits
trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')
sites_w_no_bugs=c(173029,
                  184847,
                  185060,
                  187112,
                  187223,
                  187301,
                  187304,
                  187318,
                  187780,
                  190728,
                  190768
)
PCT_sites=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//site_status_Pct_Traitstates.csv')
PCT_sites=PCT_sites[PCT_sites$sampleId %in% sites_w_no_bugs==F,]
PCT_sites=PCT_sites[,names(PCT_sites)!='sampleId']

PCT_sites_reg=PCT_sites
PCT_sites_reg$Status=ifelse(PCT_sites_reg$Status=='Probabilistic',as.numeric(1),as.numeric(0))
set.seed(42)
form2 <- as.formula(
  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
sitestatus_PCT_RF_reg=randomForest(Status~better_names,data=PCT_sites_reg,ntree=500)
sitestatus_PCT_RF_reg=randomForest(form2,data=PCT_sites_reg,ntree=500)
clipr::write_clip(sitestatus_PCT_RF_reg$importance)

sqrt(rf_model$mse[length(rf_model$mse)])
rf_model$rsq[length(rf_model$rsq)]

respvar="Status"

predictor_vars <- names(PCT_sites_reg)[names(PCT_sites_reg) %in% names(trait_table)[2:ncol(trait_table)]]
predictor_vars=c('Thermal_pref_Cold_cool_eurythermal_0_15_C',
                 'Occurance_drift_rare',
                 'Habit_prim_Clinger',
                 'Crawl_rate_high',
                 'Survive_desiccation_yes') #top5


for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'SiteStatus_Regression_top5')
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(resp, "~", paste(predictor_vars, collapse = " + ")))
  set.seed(42)
  rf_model <- randomForest(
    form,
    data = PCT_sites_reg,
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
      pred.data = PCT_sites_reg,
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




#Classification now
set.seed(42)
traits=names(PCT_sites)[names(PCT_sites) %in% names(trait_table)[2:ncol(trait_table)]]
PCT_sitesC=PCT_sites
PCT_sitesC$Status=as.factor(PCT_sitesC$Status)
sitestatus_PCT_RF_class=randomForest(Status~.,data=PCT_sitesC,ntree=500)
form2 <- as.formula(
  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
sitestatus_PCT_RF_class=randomForest(form2,data=PCT_sitesC,ntree=500)


clipr::write_clip(sitestatus_PCT_RF_class$importance)
varImpPlot(sitestatus_PCT_RF_class)

respvar=as.factor("Status")
predictor_vars=traits
predictor_vars=c('Habit_prim_Clinger',
                 'Thermal_pref_Cold_cool_eurythermal_0_15_C',
                 'Emerge_synch_abbrev_Poorly',
                 'Crawl_rate_high',
                 'Habit_prim_Burrower') #top 5 of this model
for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'SiteStatus_Classification_Top5')
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(as.factor(resp), "~", paste(predictor_vars, collapse = " + ")))
  set.seed(42)
  rf_model <- randomForest(
    form,
    data = PCT_sitesC,
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
      pred.data = PCT_sitesC,
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

##  now with PA of trait states

PA_sites=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//site_status_PA_Traitstates.csv')
PA_sites=PA_sites[PA_sites$sampleId %in% sites_w_no_bugs==F,]
PA_sites=PA_sites[,names(PA_sites)!='sampleId']

PA_sites_reg=PA_sites
PA_sites_reg$status=ifelse(PA_sites_reg$status=='Reference',as.numeric(0),as.numeric(1))
set.seed(42)
form2 <- as.formula(
  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
sitestatus_PA_RF_reg=randomForest(status~.,data=PA_sites_reg,ntree=500)
sitestatus_PA_RF_reg=randomForest(form2,data=PA_sites_reg,ntree=500)
sitestatus_PA_RF_reg
clipr::write_clip(sitestatus_PA_RF_reg$importance)
sqrt(rf_model$mse[length(rf_model$mse)])
rf_model$rsq[length(rf_model$rsq)]

respvar="status"
predictor_vars <- names(PA_sites_reg)[names(PA_sites_reg) %in% names(trait_table)[2:ncol(trait_table)]]
predictor_vars=c('Crawl_rate_high',
                 'Feed_prim_abbrev_PR',
                 'Survive_desiccation_yes',
                 'Max_body_size_abbrev_Large',
                 'Feed_prim_abbrev_CF') #top5
for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'SiteStatus_PA_Regression_top5')
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
    data = PA_sites_reg,
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
      pred.data = PA_sites_reg,
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
## now classification
PA_sites_class=PA_sites
PA_sites_class$status=as.factor(PA_sites_class$status)
set.seed(42)
form2 <- as.formula(
  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
sitestatus_PA_RF_class=randomForest(status~.,data=PA_sites_class,ntree=500)
sitestatus_PA_RF_class=randomForest(form2,data=PA_sites_class,ntree=500)
sitestatus_PA_RF_class
clipr::write_clip(sitestatus_PA_RF_class$importance)

respvar=as.factor("status")
predictor_vars <- names(PA_sites_class)[names(PA_sites_class) %in% names(trait_table)[2:ncol(trait_table)]]
predictor_vars = c('Habit_prim_Clinger',
                   'Habit_prim_Burrower',
                   'Armoring_good',
                   'Develop_slow_season',
                   'Feed_prim_abbrev_SC')
for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'SiteStatus_PA_classification_top5')
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
    data = PA_sites_class,
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
      pred.data = PA_sites_class,
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

#
# for(i in 1:length(traits)){
#   pd <- partialPlot(
#     ID_trait_mod2,
#     pred.data = PCT_sites,
#     x.var = traits[i],
#     plot = FALSE
#   )
#   df <- data.frame(
#     x = pd$x,
#     y = pd$y
#   )
#   library(ggplot2)
#
#   ggplot(df, aes(x = x, y = y)) +
#     geom_line(linewidth = 1.5) +
#     theme_classic(base_size = 18) +
#     labs(
#       title = paste("Partial Dependence on", traits[i]),
#       x = traits[i],
#       y = "Partial Effect"
#     ) +
#     theme(
#       plot.title = element_text(size = 24, hjust = 0.5),
#       axis.title = element_text(size = 20),
#       axis.text = element_text(size = 18)
#     )
#   partialPlot(ID_trait_mod2,x.var = traits[i],pred.data=PCT_sites,
#               main=paste("Partial Dependence on", traits[i]),
#               xlab='',
#               ylab='Partial Effect')
#   savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//','Top5PDP_for',names(P_for_RF)[i],'_260515.png'))
# }
#

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

for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'Prob_PA')
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
    data = PA_sites_prob,
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
  print(paste(resp,rf_model$rsq[length(rf_model$rsq)], sep='  '))

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

abun_traits=Pct_indiv_sites_prob %>%
  left_join(trait_table,by='OTU')
trait_cols_abun=names(trait_table[,2:ncol(trait_table)])

pct_indivs_w_Trait_summary <- abun_traits %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    across(
      all_of(trait_cols_abun),
      ~ 100 * sum(split_count * .x, na.rm = TRUE) / sum(split_count, na.rm = TRUE),
      .names = "Pct_{.col}"
    )
  )
environmental_for_indivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
environmental_for_indivs=environmental_for_indivs[,names(environmental_for_indivs)%in% c(predictor_vars,'sampleId')]
names(environmental_for_indivs)[names(environmental_for_indivs)=='sampleId']<-'sample_id'
pct_indivs_w_Trait_summary=plyr::join(pct_indivs_w_Trait_summary,environmental_for_indivs)

library(randomForest)
library(ggplot2)
names(pct_indivs_w_Trait_summary)[2:26]<-gsub("^Pct_", "", names(pct_indivs_w_Trait_summary)[2:26])
pct_indivs_w_Trait_summary=pct_indivs_w_Trait_summary[,names(pct_indivs_w_Trait_summary)!='sample_id']
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

colSums(PA_sites_R[predictor_vars])


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


Pct_sites_R=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
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
    data = PA_sites_prob,
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
 ## % I and %D at sites (prob only has indivs)

P_comm=ProbOs ## read this in as csv and filter out fails too
Probs=ProbOs[row.names(ProbOs) %in% sites_w_no_bugs==F,]
P_comm$Site=row.names(ProbOs)
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

environmental_for_indivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
environmental_for_indivs=environmental_for_indivs[,names(environmental_for_indivs)%in% c(predictor_vars,'sampleId')]
names(environmental_for_indivs)[1]<-'Site'

PA_P_response_wtraits=as.data.frame(plyr::join(PA_P_response,environmental_for_indivs))
PA_P_response_wtraits=PA_P_response_wtraits[PA_P_response_wtraits$Site %in% c(failed_sites, sites_w_no_bugs)==F,]
I_P_response=PA_P_response_wtraits[,c('Percent_I',predictor_vars)]

predictor_vars=names(I_P_response)[2:ncol(I_P_response)]
response_varsPA=names(I_P_response)[1]
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Increase")
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
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Indiv_Increase")
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
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Decrease")
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
for (resp in response_varsPA) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results","Pct_Indiv_Decrease")
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



