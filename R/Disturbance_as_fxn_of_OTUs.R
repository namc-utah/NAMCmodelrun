#Disturbance as fxn of OTUs

disturbs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Disturbance_by_site.csv')

trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')

site_taxa_matrix=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//sites_PA_matrix.csv")
site_taxa_matrix=site_taxa_matrix[site_taxa_matrix$sampleId %in% disturbs$sampleId,]

PA_prob_w_d=plyr::join(disturbs,site_taxa_matrix,by='sampleId','left')
PA_prob_w_d=PA_prob_w_d[!duplicated(PA_prob_w_d$sampleId),]
PA_prob_w_d=PA_prob_w_d[,names(PA_prob_w_d) !='sampleId']

PA_prob_w_d$Disturbance =as.factor(PA_prob_w_d$Disturbance)

PA_prob_w_d=PA_prob_w_d[,names(PA_prob_w_d) %in% c('Disturbance',trait_table$OTU)]
#taxa which do not occur at any site, R or P.
PA_prob_w_d=PA_prob_w_d[,names(PA_prob_w_d) %in% c('Farula','Stactobiella')==F]
predictor_vars=names(PA_prob_w_d)[2:ncol(PA_prob_w_d)]

respvar=names(PA_prob_w_d)[1]
predictor_vars=c('Ephemerella',
                 'Epeorus',
                 'RHYACOPHILA',
                 'Rhithrogena',
                 'Micrasema')
for (j in respvar) {

  # -----------------------------
  # 1. Create output folder
  # -----------------------------
  base_dir <- file.path(path.expand("~"), "New_for_SFS", "Disturbance_all_sites_OTUs_top5")
  outdir <- file.path(base_dir, j)
  dir.create(outdir, recursive = TRUE, showWarnings = T)
  #outdir <- paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//RF_results//",resp)
  #dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # 2. Build formula + model
  # -----------------------------
  form <- as.formula(
    paste(j, "~", paste(predictor_vars, collapse = " + ")))
  set.seed(42)
  rf_model <- randomForest(
    form,
    data = PA_prob_w_d,
    importance = TRUE)
  rf_model
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
      title = paste("Variable Importance -", j),
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
      pred.data = PA_prob_w_d,
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


rf_model
rf_mod_imp=as.data.frame(rf_model$importance)
rf_mod_imp[order(rf_mod_imp$MeanDecreaseGini, decreasing = T),]
