#Predictor inspection
query<- ("SELECT
    r.sample_id AS sampleId,
    r.translation_id AS translationId,
    z.taxonomy_id AS taxonomyId,
    t.scientific_name AS scientificName,
    z.split_count,
    p.abbreviation,
    sp.predictor_value,
    (
    SELECT DISTINCT COALESCE((
        WITH RECURSIVE TaxonomyHierarchy AS (/* create the CTE, essentially a for loop*/
            SELECT taxonomy_id, scientific_name, parent_id /* this is the base query of the CTE. Checking for matches {*/
            FROM taxonomy
            WHERE taxonomy_id = z.taxonomy_id  /* }*/
            UNION ALL /* combines recursive and base queries*/
            SELECT t.taxonomy_id, t.scientific_name, t.parent_id /* this is the looping part {*/
            FROM taxonomy t
            JOIN TaxonomyHierarchy th ON t.taxonomy_id = th.parent_id
        ) /* }*/ /*leave the recursive section*/
        SELECT g.scientific_name /* now let's take the scientific name*/
        FROM TaxonomyHierarchy th /* from the big nasty CTE we just made*/
        LEFT JOIN translation_taxa f ON th.taxonomy_id = f.taxonomy_id /* join f and th, by taxId, to start getting the OTU name*/
        LEFT JOIN translation_groups g ON f.translation_group_id = g.translation_group_id /*join g and f. Now we can get the OTU name*/
        WHERE g.translation_id = 2 /* filter*/
        LIMIT 1 /* make sure only one row is returned*/
    ), NULL)
) As OTU, /*name the column*/
  /*the rest is all pretty standard SQL and joins*/
  z.split_count AS splitCount
FROM
rarefactions r
JOIN
rarefied_organisms z ON r.rarefaction_id = z.rarefaction_id
JOIN
taxonomy t ON z.taxonomy_id = t.taxonomy_id
join samples s on s.sample_id = r.sample_id
join site_predictors sp on sp.site_id = s.site_id
join predictors p on sp.predictor_id = p.predictor_id
join model_predictors mp on mp.predictor_id = sp.predictor_id
WHERE
r.sample_id in (select sample_id from samples where sample_id
                /* change to get a group of samples via a project, or box, for example
                or just change to sample_id =...*/
                  /*this huge query is pulling all of AREMP*/

                  in (select sample_id from samples where box_id in (select box_id from boxes where customer_id = 4396)) )
AND r.translation_id = 2
and mp.model_id = 25;")
#connect to the database
mydb=DBI::dbConnect(RSQLite::SQLite(),'C://NAMC_S3//LegacyDatabases//instar.sqlite')
#now execute the query and save as a data frame
ben_dat<-DBI::dbGetQuery(mydb,query)
library(tidyr)
library(dplyr)
ben_dat<-ben_dat[is.na(ben_dat$OTU)==F,]
# First, separate taxon data
taxon_wide <- ben_dat %>%
  group_by(sampleId, OTU) %>%
  summarize(split_count = sum(split_count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = OTU, values_from = split_count, values_fill = 0)

# 2. Wide-format predictor values — ensures only one column per predictor_id
predictor_wide <- ben_dat %>%
  select(sampleId, abbreviation, predictor_value) %>%
  distinct() %>%
  pivot_wider(names_from = abbreviation, values_from = predictor_value)

# 3. Merge: join both by sample_id (one row per sample)
final_df <- left_join(taxon_wide, predictor_wide, by = "sampleId")
final_df<-as.data.frame(final_df)
rownames(final_df) <- final_df$sampleId
final_df <- final_df %>% select(-sampleId)

if(0){
library(randomForestSRC)

# -------------------------------
# 1. Simulate data
set.seed(123)
n_sites <- 100
n_taxa <- 5

# Environmental gradient (e.g., % urbanization)
env_data <- data.frame(
  Urbanization = runif(n_sites, 0, 100),    # 0–100% urban
  Elevation = runif(n_sites, 100, 500)
)

# Simulate taxon response to disturbance (% Urban)
# Taxon1 decreases, Taxon5 increases, others random
bio_data <- data.frame(
  Taxon1 = rbinom(n_sites, 1, prob = 1 - env_data$Urbanization/150),  # decreasing
  Taxon2 = rbinom(n_sites, 1, 0.5),                                   # neutral
  Taxon3 = rbinom(n_sites, 1, runif(n_sites, 0.3, 0.7)),              # random
  Taxon4 = rbinom(n_sites, 1, prob = 0.7 - env_data$Urbanization/300),# slightly decreasing
  Taxon5 = rbinom(n_sites, 1, prob = env_data$Urbanization/100)       # increasing
)

# Combine into a modeling dataset
data <- cbind(env_data, bio_data)

# -------------------------------
# 2. Fit MRF model
taxa_names <- paste0("Taxon", 1:5)
mrf_model <- rfsrc(
  formula = as.formula(paste("cbind(", paste(taxa_names, collapse = ", "), ") ~ .")),
  data = data,
  ntree = 300
)

# -------------------------------
# 3. Create a prediction gradient (% Urban from 0 to 100)
grad_data <- data.frame(
  Urbanization = seq(0, 100, length.out = 50),
  Elevation = mean(env_data$Elevation)  # hold Elevation constant
)

pred <- predict(mrf_model, newdata = grad_data)
pred_probs <- sapply(pred$regrOutput, function(x) x$predicted)
pred_probs <- as.data.frame(pred_probs)
pred_probs$Urbanization <- grad_data$Urbanization

# -------------------------------
# 4. Plot taxon response curves
library(ggplot2)
library(tidyr)

plot_data <- pivot_longer(pred_probs, cols = starts_with("Taxon"),
                          names_to = "Taxon", values_to = "Probability")

ggplot(plot_data, aes(x = Urbanization, y = Probability, color = Taxon)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(title = "Taxon Response Across Urbanization Gradient",
       y = "Predicted Probability of Occurrence") +
  theme(legend.position = "bottom")
}
