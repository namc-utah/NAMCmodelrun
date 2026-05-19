### Ordination for traits in environmental space
#site x taxa
#taxa x trait
#site x environment

#site by taxa
#this is the "L matrix"
s_x_taxa=PA

#taxa x trait
#This is the "Q matrix"
t_x_t=trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')
t_x_t=t_x_t[t_x_t$OTU %in% c('Farula','Stactobiella')==F,]
t_x_t=t_x_t[!duplicated(t_x_t),]
row.names(t_x_t)=t_x_t$OTU;t_x_t=t_x_t[,-1]



#site x environment
#this is the "R matrix"
s_x_e=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//ordination_sites_envs.csv')
row.names(s_x_e)=s_x_e$sampleId
s_x_e=s_x_e[,names(s_x_e)!='sampleId']
just_ref=s_x_taxa[s_x_taxa$status=='Reference',]
just_prob=s_x_taxa[s_x_taxa$status=='Probabilistic',]

just_ref_env=s_x_e[s_x_e$sampleId %in% row.names(just_ref),]

just_prob_envs=s_x_e[s_x_e$sampleId %in% row.names(just_prob),]

# check alignments
all(rownames(s_x_e) == rownames(s_x_taxa))
#rows match up, so just reassign to a new variable
R_matrix=s_x_e

Q_matrix <-t_x_t[colnames(s_x_taxa[,-ncol(s_x_taxa)]), , drop = FALSE]

L_matrix=s_x_taxa[,-ncol(s_x_taxa)]

better_Trait_names=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//better_trait_names.csv')
better_Trait_names$Trait==names(t_x_t)
names(Q_matrix)<-better_Trait_names$Abbrev

better_env_names=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//better_env_names.csv')
names(R_matrix)==better_env_names$EnvVar
names(R_matrix)<-better_env_names$Abbrev
all(rownames(R_matrix) == rownames(L_matrix))
all(colnames(L_matrix) == rownames(Q_matrix))




library(ade4)

# 1. Correspondence Analysis on the community data (L)
ord_L <- dudi.coa(L_matrix, scannf = FALSE, nf = 2)

# 2. PCA on continuous environmental data (R) weighted by sites
ord_R <- dudi.pca(R_matrix, row.w = ord_L$lw, scannf = FALSE, nf = 2)

# 3. PCA (or Hill-Smith) on binary trait data (Q) weighted by species
ord_Q <- dudi.pca(Q_matrix, row.w = ord_L$cw, scannf = FALSE, nf = 2)

rlq_result <- rlq(dudiR = ord_R, dudiL = ord_L, dudiQ = ord_Q, scannf = FALSE, nf = 2)

# Plot the multi-panel summary graph
plot(rlq_result)


s.arrow(rlq_result$l1, boxes = FALSE)

# Plot the trait coefficients (groupings for binary traits)
s.label(rlq_result$c1, add.plot = TRUE, clabel = 1.2)


#fourth_test <- fourthcorner(R_matrix, L_matrix, Q_matrix, modeltype = 6, nrepet = 999)

# Plot the significance matrix
#plot(fourth_test, alpha = 0.05, stat = "D2")




library(ggplot2)
library(ggrepel) # Prevents text labels from overlapping
Fish_trait_names=c("Emerge_synch_abbrev_Poorly",
                   "Feed_prim_abbrev_CG",
                   "Female_disp_abbrev_High",
                   "Habit_prim_Burrower",
                   "Habit_prim_Climber",
                   "Habit_prim_Sprawler",
                   "Habit_prim_Swimmer",
                   "Rheophily_abbrev_depo",
                   "Lifespan_very_short")
# 1. Extract environmental continuous variables (vectors/arrows)
trait_coords <- rlq_result$co
trait_coords$OG=better_Trait_names$Trait
trait_coords=trait_coords[trait_coords$OG %in% Fish_trait_names,]
colnames(trait_coords) <- c("Axis1", "Axis2")
trait_coords$Variable <- rownames(trait_coords)

# 2. Extract binary trait variables (points/labels)
env_coords <- rlq_result$li
colnames(env_coords) <- c("Axis1", "Axis2")
env_coords$Trait <- rownames(env_coords)

# Initialize the plot space
ggplot() +
  # Add subtle origin crosshairs
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.4) +

  # 1. Plot Continuous Environmental Variables as Arrows
  geom_segment(data = trait_coords,
               aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.8) +

  # Label the environmental arrows (moved further away using nudging and higher padding)
  # geom_text_repel(data = trait_coords,
  #                 aes(x = Axis1, y = Axis2, label = Variable),
  #                 color = "black", fontface = "bold", size = 4.5,
  #                 box.padding = 1.2,        # Higher values push text further from points
  #                 point.padding = 0.8,      # Keeps text clear of the arrow tips
  #                 ) +

  # 2. Plot Binary Traits as Points
  geom_point(data = env_coords,
             aes(x = Axis1, y = Axis2),
             color = "dodgerblue4", size = 3, alpha = 0.8) +

  # Label the trait points safely without overlapping
  geom_text_repel(data = env_coords,
                  aes(x = Axis1, y = Axis2, label = Trait),
                  color = "dodgerblue4", fontface = "italic", size = 6,
                  box.padding = 0.5, segment.color = "gray70",
                  nudge_x = env_coords$Axis1 * 0.15, # Pushes text slightly outward radially
                                   nudge_y = env_coords$Axis2 * 0.5,
                  min.segment.length = Inf) +

  # 3. Polish the Layout and Styling
  labs(
    x = "RLQ Axis 1",
    y = "RLQ Axis 2"
  ) +
  theme_classic() + # 1. BASE THEME GOES FIRST
  theme(            # 2. CUSTOM OVERRIDES GO SECOND
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25), # Now successfully active!
    axis.title.x = element_text(margin=ggplot2::margin(t=20)),
    axis.title.y  = element_text(margin=ggplot2::margin(r=20)),
    legend.key.size = unit(1, "lines"),
    axis.ticks = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.13)
  )
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//traits_in_envspace_all_nolabs.png')
