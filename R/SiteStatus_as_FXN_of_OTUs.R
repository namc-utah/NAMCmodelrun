## Site status as a fxn of OTUs
site_PA_mat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//sites_PA_matrix.csv')
site_PA_mat=site_PA_mat[site_PA_mat$sampleId %in% c('HAWK-33','WE-131','WE-209')==F,]
row.names(site_PA_mat)=site_PA_mat$sampleId;site_PA_mat=site_PA_mat[,-1]
site_PA_mat=site_PA_mat[,names(site_PA_mat) %in% c('Farula','Stactobiella')==F,]

trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')

respvar='status'
predictor_vars=names(site_PA_mat)[1:(ncol(site_PA_mat)-1)]

predictor_vars=c('Ephemerella',
                        'Epeorus',
                        'RHYACOPHILA',
                        'Serratella',
                        'Micrasema')
##regression
library(randomForest)
site_PA_mat_reg=site_PA_mat
site_PA_mat_reg$status=ifelse(site_PA_mat_reg$status=='Probabilistic',as.numeric(1),as.numeric(0))
set.seed(42)
form2 <- as.formula(
  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
site_PA_mat_reg_RF=randomForest(status~.,data=site_PA_mat_reg,ntree=500)
site_PA_mat_reg_RF=randomForest(form2,data=site_PA_mat_reg,ntree=500)
clipr::write_clip(site_PA_mat_reg_RF$importance)
site_PA_mat_reg_RF
sqrt(site_PA_mat_reg_RF$mse[length(site_PA_mat_reg_RF$mse)])
site_PA_mat_reg_RF$rsq[length(site_PA_mat_reg_RF$rsq)]



for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "RF_results",'SiteStatus_FXN_ofOTUS_Regression_top5')
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
    data = site_PA_mat_reg,
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
      pred.data = site_PA_mat_reg,
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

site_PA_mat_class=site_PA_mat
#site_PA_mat_class$status=ifelse(site_PA_mat_class$status=='Probabilistic',as.factor(site_PA_mat_class$status),as.factor(site_PA_mat_class$status))
set.seed(42)
#form2 <- as.formula(
#  paste(respvar, "~", paste(predictor_vars, collapse = " + ")))
site_PA_mat_class_RF=randomForest(as.factor(status)~.,data=site_PA_mat_class,ntree=500)

clipr::write_clip(site_PA_mat_class_RF$importance)

sqrt(site_PA_mat_class_RF$mse[length(site_PA_mat_class_RF$mse)])
site_PA_mat_class_RF$rsq[length(site_PA_mat_class_RF$rsq)]


site_PA_mat_class$status=as.factor(site_PA_mat_class$status)

predictor_vars=names(site_PA_mat)[1:(ncol(site_PA_mat)-1)]

respvar='status'
predictor_vars=c('Ephemerella',
                 'Epeorus',
                 'RHYACOPHILA',
                 'Serratella',
                 'Micrasema')
for (resp in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "SiteStatus_FXN_ofOTUS_Classification")
  outdir <- file.path(base_dir, resp)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(as.factor(respvar), "~", paste(predictor_vars, collapse = " + ")))
  set.seed(42)
  rf_model <- randomForest(
    form,
    data = site_PA_mat_class,
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
      pred.data = site_PA_mat_class,
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
