#number of reference sites / prob sites
#for creating frequency
n_ref= 656
n_p=349
library(tidyverse)
#read in data
CI_plot_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//CI_plot_dat.csv')
#create frequency
#Fo / all sites
CI_plot_dat$frequency=ifelse(CI_plot_dat$status=='Probabilistic',
                              CI_plot_dat$Fo/n_p,
                              CI_plot_dat$Fo/n_ref)
#force name to be factor for plotting
CI_plot_dat$taxon_proper<-as.factor(CI_plot_dat$taxon_proper)
CI_plot_dat$Group<-as.factor(CI_plot_dat$Group)

#plot saving fxn
savp = function(W,H,fn) {
  dev.copy(dev=png,file=fn,wi=W,he=H,un="in",res=650)
  dev.off()
}

#this will create positions that ggplot will plot taxa on
taxa_positions <- CI_plot_dat %>%
  distinct(taxon_proper, .keep_all = TRUE) %>%
  arrange(Group, taxon_proper) %>%      # Ascending by Group then Taxon
  dplyr::mutate(ypos = row_number()) %>%
  select(taxon_proper, ypos, Group)

#join the positions back to the main data frame
CI_plot_dat <- CI_plot_dat %>%
  left_join(taxa_positions, by = "taxon_proper")
#drop extra group that arises via join
CI_plot_dat<-CI_plot_dat[,names(CI_plot_dat)!='Group.y']
names(CI_plot_dat)[1]<-'Group'
# Compute Order rectangles
#this will help color code taxa
order_rects <- taxa_positions %>%
  group_by(Group) %>%
  dplyr::summarize(
    ymin = min(ypos),
    ymax = max(ypos),
    n_taxa = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Ensure rectangle visible for single-taxon groups
    ymin_pad = ifelse(n_taxa == 1, ymin - 0.5, ymin - 0.2),
    ymax_pad = ifelse(n_taxa == 1, ymax + 0.5, ymax + 0.2),
    # Center label
    y_label = (ymin_pad + ymax_pad) / 2
  )

#add padding to the rectangles so we can
#show groups with only 1 record against larger sizes
order_rects <- order_rects %>%
  dplyr::mutate(
    # Compute padded rectangle coordinates
    ymin_pad = ifelse(ymax - ymin < 1, ymin - 0.5, ymin - 0.2),
    ymax_pad = ifelse(ymax - ymin < 1, ymax + 0.5, ymax + 0.2),
    # Label in center
    y_label = (ymin_pad + ymax_pad) / 2
  )
# Define fill colors
#which will be used for the rectangles
fill_colors <- c(
  'Amphipoda'='firebrick4',
  'Annelida'='red',
  'Arachnida'='firebrick2',
  'Bivalvia'='indianred3',
  'Coleoptera'='pink',
  'Decapoda'='orange4',
  'Diptera'='orange',
  'Ephemeroptera'='dodgerblue4',
  'Gastropoda'= 'blue',
  'Hemiptera'='dodgerblue',
  'Hirudinea'='darkslategray4',
  'Isopoda'='darkslategray1',
  'Lepidoptera'='violetred4',
  'Megaloptera'='violetred',
  'Non-Annelid worms'='violet',
  'Odonata'='khaki2',
  'Platyhelminthes'='lemonchiffon4',
  'Plecoptera'= 'grey',
  'Trichoptera'='black'
)
#plot the data
ggplot(CI_plot_dat, aes(x = frequency, y = ypos, color = status)) +
  # 1. Draw Group rectangles first
  geom_rect(
    data = order_rects,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_pad, ymax = ymax_pad, fill = Group),
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  #Points and error bars
  geom_point(size=1) +
  geom_errorbar(aes(y = ypos, xmin =LL_n, xmax = UL_n), width = 0.2) +
  # Order labels (omitted because it is just too messy)
  # geom_text(
  #   data = order_rects,
  #   # Corrected x position to use a numeric value
  #   aes(x = -0.02, y = y_label, label = Group),
  #   hjust = 1,
  #   size = 2,
  #   inherit.aes = FALSE
  # ) +
  # Axis scales and colors
  scale_y_continuous(
    breaks = taxa_positions$ypos,
    labels = taxa_positions$taxon_proper,
    trans = "reverse" # ascending from top to bottom,
    #so Amphipoda is plotted first
  ) +
  #change point color based on status
  #used these deep colors instead for visibility
  scale_color_manual(values = c('Reference' = 'blue', 'Probabilistic' = 'red')) +
  #fill the rectangles
  scale_fill_manual(values = fill_colors) +
  # graph settings
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = 'Frequency\n(proportion of sites present)',
    y = 'Taxon',
    color = 'Site status',
    fill='Order/Class/Phylum'
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 2, hjust = 1),
    axis.ticks = element_line(size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
    #plot.margin = margin(5, 5, 5, 20))+
  #move legend for readability
  theme(legend.position = 'left')
#save the file
savp(10,10,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//test.png')

small_plot_names=CI_plot_dat$taxon_proper[CI_plot_dat$frequency>=0.1]
small_plot<-CI_plot_dat[CI_plot_dat$taxon_proper %in% small_plot_names,]
taxa_positions2 <- small_plot %>%
  distinct(taxon_proper, .keep_all = TRUE) %>%
  arrange(Group, taxon_proper) %>%      # Ascending by Group then Taxon
  dplyr::mutate(ypos = row_number()) %>%
  select(taxon_proper, ypos, Group)

#join the positions back to the main data frame
small_plot<- small_plot %>%
  left_join(taxa_positions2, by = "taxon_proper")
#drop extra group that arises via join
small_plot<-small_plot[,names(small_plot)%in%c('Group.y','ypos.x')==F]
names(small_plot)[1]<-'Group';names(small_plot)[ncol(small_plot)]<-'ypos'
# Compute Order rectangles
#this will help color code taxa
order_rects2 <- taxa_positions2 %>%
  group_by(Group) %>%
  dplyr::summarize(
    ymin = min(ypos),
    ymax = max(ypos),
    n_taxa = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Ensure rectangle visible for single-taxon groups
    ymin_pad = ifelse(n_taxa == 1, ymin - 0.5, ymin - 0.2),
    ymax_pad = ifelse(n_taxa == 1, ymax + 0.5, ymax + 0.2),
    # Center label
    y_label = (ymin_pad + ymax_pad) / 2
  )

#add padding to the rectangles so we can
#show groups with only 1 record against larger sizes
order_rects2 <- order_rects2 %>%
  dplyr::mutate(
    # Compute padded rectangle coordinates
    ymin_pad = ifelse(ymax - ymin < 1, ymin - 0.5, ymin - 0.2),
    ymax_pad = ifelse(ymax - ymin < 1, ymax + 0.5, ymax + 0.2),
    # Label in center
    y_label = (ymin_pad + ymax_pad) / 2
  )


small_plot$LL_n[small_plot$status=='Probabilistic']<-NA
small_plot$UL_n[small_plot$status=='Probabilistic']<-NA
Probs=small_plot[small_plot$status=='Probabilistic',]
ggplot(small_plot[small_plot$status=='Reference',], aes(x = frequency, y = ypos, color = status)) +
  # 1. Draw Group rectangles first
  geom_rect(
    data = order_rects2,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_pad, ymax = ymax_pad, fill = Group),
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  #Points and error bars
  geom_point(size=2) +
  geom_point(data=Probs,
             aes(y=ypos,x=frequency,
                 col=status),
             size=2)+
  geom_errorbar(aes(y = ypos, xmin =LL_n, xmax = UL_n), width = 0.2) +
  # Order labels (omitted because it is just too messy)
  # geom_text(
  #   data = order_rects,
  #   # Corrected x position to use a numeric value
  #   aes(x = -0.02, y = y_label, label = Group),
  #   hjust = 1,
  #   size = 2,
  #   inherit.aes = FALSE
  # ) +
  # Axis scales and colors
  scale_y_continuous(
    breaks = taxa_positions2$ypos,
    labels = taxa_positions2$taxon_proper,
    trans = "reverse",
    expand = expansion(mult = c(0.05, 0.05))) +
  #change point color based on status
  #used these deep colors instead for visibility
  scale_color_manual(values = c('Reference' = 'blue', 'Probabilistic' = 'red')) +
  #fill the rectangles
  scale_fill_manual(values = fill_colors) +
  # graph settings
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = 'Frequency\n(proportion of sites present)',
    y = 'OTU',
    color = 'Site status',
    fill='Taxonomic group'
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 5, hjust = 1,color='black'),
    axis.ticks = element_line(size = 0.5),
    panel.grid.major = element_line(color='grey50'))+
    #panel.grid.minor = element_blank())+
  #plot.margin = margin(5, 5, 5, 20))+
  #move legend for readability
  theme(legend.position = 'left')
savp(10,10,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//small_CI_plot.png')


