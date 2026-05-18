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
