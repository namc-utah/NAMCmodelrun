#"pivot table"

piv_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//all_ratios_251121.csv')
piv_dat$Inc=ifelse(piv_dat$Presponse=='Increaser',1,0)
piv_dat$Dec=ifelse(piv_dat$Presponse=='Decreaser',1,0)


piv_clean=piv_dat%>%
  filter(!is.na(Group)) %>%           # remove empty groups (if any)
  mutate(
    Inc = ifelse(is.na(Inc), 0, Inc), # convert NA to 0
    Dec = ifelse(is.na(Dec), 0, Dec)
  )
piv2=piv_clean %>%
  group_by (taxon, Group) %>%
  dplyr::summarise(Fo=sum(pFo),
                  Fe=sum(pFe),
                  ratio=Fo/Fe,
                  response=ifelse(ratio > 0.8 & ratio < 1.2,"Neutral",
                                  ifelse( ratio < 0.8, 'Decreaser','Increaser')),
                  .groups='drop')

piv2$Inc=ifelse(piv2$response=='Increaser',1,0)
piv2$Dec=ifelse(piv2$response=='Decreaser',1,0)

piv2 = piv2 %>% filter(!is.na(Group)) %>%           # remove empty groups (if any)
  mutate(
    Inc = ifelse(is.na(Inc), 0, Inc), # convert NA to 0
    Dec = ifelse(is.na(Dec), 0, Dec)
  )


X=piv_clean[piv_clean$Ecoregion!='Other',]
X=X[X$Fo>0,]
O=piv_clean[piv_clean$Ecoregion=='Other',]
O=O[O$Fo>0,]


oth_summ=as.data.frame(O %>%
  group_by(Group) %>%
  dplyr::summarise(
    Total_Taxa = n_distinct(taxon),
    Num_Increasers = sum(Inc),
    Num_Decreasers = sum(Dec),
    Prop_Increaser = Num_Increasers / Total_Taxa,
    Prop_Decreaser = Num_Decreasers / Total_Taxa,
    .groups = "drop"
  ))
oth_total_tab=O %>%
  dplyr::summarise(
    Group = "Total",
    Total_Taxa = n_distinct(taxon),
    Num_Increasers = sum(Inc),
    Num_Decreasers = sum(Dec),
    Prop_Increaser = Num_Increasers / Total_Taxa,
    Prop_Decreaser = Num_Decreasers / Total_Taxa
  )


Other_final_tab=bind_rows(oth_total_tab,oth_summ)
Other_final_tab[,5:6]=Other_final_tab[,5:6]*100
clipr::write_clip(Other_final_tab)

xer_summ=as.data.frame(X %>%
                         group_by(Group) %>%
                         dplyr::summarise(
                           Total_Taxa = n_distinct(taxon),
                           Num_Increasers = sum(Inc),
                           Num_Decreasers = sum(Dec),
                           Prop_Increaser = Num_Increasers / Total_Taxa,
                           Prop_Decreaser = Num_Decreasers / Total_Taxa,
                           .groups = "drop"
                         ))
xer_total_tab=X %>%
  dplyr::summarise(
    Group = "Total",
    Total_Taxa = n_distinct(taxon),
    Num_Increasers = sum(Inc),
    Num_Decreasers = sum(Dec),
    Prop_Increaser = Num_Increasers / Total_Taxa,
    Prop_Decreaser = Num_Decreasers / Total_Taxa
  )
xer_final_tab=bind_rows(xer_total_tab,xer_summ)
xer_final_tab[,5:6]=xer_final_tab[,5:6]*100

clipr::write_clip(xer_final_tab)
all_summ=as.data.frame(piv2 %>%
                         group_by(Group) %>%
                         dplyr::summarise(
                           Total_Taxa = n_distinct(taxon),
                           Num_Increasers = sum(Inc),
                           Num_Decreasers = sum(Dec),
                           Prop_Increaser = Num_Increasers / Total_Taxa,
                           Prop_Decreaser = Num_Decreasers / Total_Taxa,
                           .groups = "drop"
                         ))
all_total_tab=piv2 %>%
  dplyr::summarise(
    Group = "Total",
    Total_Taxa = n_distinct(taxon),
    Num_Increasers = sum(Inc),
    Num_Decreasers = sum(Dec),
    Prop_Increaser = Num_Increasers / Total_Taxa,
    Prop_Decreaser = Num_Decreasers / Total_Taxa
  )
all_final_tab=bind_rows(all_total_tab,all_summ)
all_final_tab[,5:6]=all_final_tab[,5:6]*100
all_final_tab

clipr::write_clip(all_final_tab)
