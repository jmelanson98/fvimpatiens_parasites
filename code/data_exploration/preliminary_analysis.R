## UBC Parasite Analysis Code

#load packages
library(dplyr)
library(stringr)
library(tidyverse)
library(vegan)
library(labdsv)
library(nlme)
library(gtools)
library(lme4)


#load dataframes
specimenData = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/allsamplepoints.csv", sep = ","))
sampleEffort = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022sampledata.csv", sep = ","))
vegData = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022vegetationdata.csv", sep = ","))
parasiteScores = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022parasitedata.csv", sep = ","), header = T)

#correct headers
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite", "sub2")
colnames(sampleEffort) = sampleEffort[1,]
sampleEffort = sampleEffort[-1,]
colnames(vegData) = vegData[1,]
vegData = vegData[-1,]
colnames(parasiteScores) = parasiteScores[1,]
parasiteScores = parasiteScores[-1,]

#join dataframes to get specimenData, GPS/subsite data, and parasite scores into a single dataframe
#leave sample effort separate

#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,-2]
samplePoints = samplePoints[,-2] #run this twice to get rid of both columns lol

#join sample points to specimen data
allspecimens = left_join(specimenData, samplePoints, by = "sample_pt")

#join parasite scores
parasiteScores$barcode_id[parasiteScores$barcode_id == "HR1_14_01"] = "HR1_14A_01"
parasiteScores = mutate(parasiteScores, across(-c(barcode_id), as.numeric))
allspecimens = left_join(allspecimens, parasiteScores, by = "barcode_id")

#clean up: remove males / bees without final_id / add periods
allspecimens %>% filter(final_id != "") %>%
  filter(final_id != "B. huntii") %>%
  mutate(period = case_when(round < 3 ~ "1",
                          round > 2 & round < 5 ~ "2",
                          round > 4 & round < 7 ~ "3",
                          round > 6 & round < 9 ~ "4",
                          round > 8 ~ "5")) %>%
  relocate(period, .after = round) -> cleanedspecimens
cleanedspecimens = cleanedspecimens[str_detect(cleanedspecimens$notes, "released|male|lost", negate = T),]

#join sample points to sample effort data frame
sampleEffort = left_join(sampleEffort, samplePoints, by = "sample_pt")
sampleEffort$date = paste(sampleEffort$day, sampleEffort$month, sampleEffort$year)
sampleEffort$date = gsub(" ", "", sampleEffort$date, fixed = TRUE)
sampleEffort$date <- as.POSIXlt(sampleEffort$date, format = "%d%b%y")
sampleEffort$julian_date = sampleEffort$date$yday

#####CODE FOR IF YOU WANT TO LOOK AT COMMUNITY COMP / FLORAl COMP AT THE INDIVIDUAL SAMPLING EVENT LEVEL

#add Julian dates to specimen dataframe
datesdf = sampleEffort[,colnames(sampleEffort) %in% c("sample_id", "julian_date")]
specimens = left_join(cleanedspecimens, datesdf, by = "sample_id")

#add bee richness to subsite sample effort data frame
richness_df = group_by(specimens, sample_id) %>%
  count(final_id) %>%
  count()
colnames(richness_df) = c("sample_id", "bombus_richness")
specimens = left_join(specimens, richness_df, by = "sample_id")
specimens$bombus_richness[is.na(specimens$bombus_richness)] <- 0

#add bee diversity (Shannon Index) to specimen data frame
community_df_long = group_by(specimens, sample_id) %>%
  count(final_id)
community_df_wide = as.data.frame(pivot_wider(community_df_long, names_from = "final_id", values_from = "n"))
community_df_wide[is.na(community_df_wide)] <- 0
rownames(community_df_wide) = community_df_wide[,1]
community_df = community_df_wide[,-1]
community_df_wide$H <- vegan::diversity(community_df)
colnames(community_df_wide) = c("sample_id", "impatiens_abundance", "rufocinctus_abundance", "californicus_abundance", "vos_abundance", "flavifrons_abundance", "mixtus_abundance", "melanopygus_abundance", "sitkensis_abundance", "shannon_diversity")

specimens = left_join(specimens, community_df_wide, by = "sample_id")

#add bee abundance to subsite sample effort data frame
specimens$bombus_abundance = specimens$impatiens_abundance + specimens$melanopygus_abundance + specimens$mixtus_abundance + specimens$flavifrons_abundance + specimens$rufocinctus_abundance + specimens$californicus_abundance + specimens$vos_abundance

# #plot impatiens abundance as a function of julian date
# library(ggplot2)
# colnames(sampleEffort)[18] = "impatiens_abundance"
# subsiteEffort$date = paste(subsiteEffort$day, subsiteEffort$month, subsiteEffort$year)
# subsiteEffort$date = gsub(" ", "", subsiteEffort$date, fixed = TRUE)
# subsiteEffort$date <- as.POSIXlt(subsiteEffort$date, format = "%d%b%y")
# subsiteEffort$julian_date = subsiteEffort$date$yday
# ggplot(subsiteEffort,aes(x = `B. impatiens`,y = bombus_abundance), fill = julian_date) + 
#   geom_point(aes(color=julian_date)) +
#   theme_bw()
# correlation <- cor(subsiteEffort$`B. impatiens`, subsiteEffort$bombus_abundance, method = 'pearson')
# 

#create a list of all the flowers visited by bees
beeflowers = unique(specimens$active_flower)
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu", negate = T)]

#prepare veg data frame
colnames(vegData)[3] = "sample_pt"
vegData = left_join(vegData, samplePoints, by = "sample_pt")
vegData -> subsiteVeg

#calculate floral richness at each sample_id
columnlist = c("sample_id", beeflowers)
floral_df_wide = subsiteVeg[,colnames(subsiteVeg) %in% columnlist]
floral_df_wide = mutate(floral_df_wide, across(-c(sample_id), as.numeric))
floral_df_wide %>%
 mutate_at(vars(-1), ~ +(.x > 0)) -> floral_point_binary
floral_point_binary$point_floral_richness = rowSums(floral_point_binary[-1], na.rm = TRUE)
point_richness = floral_point_binary[,c(1,ncol(floral_point_binary))]
specimens = left_join(specimens, point_richness, by = "sample_id")

#calculate floral abundance at each sample_id
floral_df_wide %>%
  pivot_longer(!sample_id, names_to = "flower", values_to = "floral_abundance") -> floral_df_long
floral_df_long$floral_abundance = 10^(floral_df_long$floral_abundance - 1)
abundance_aggregate = aggregate(floral_abundance ~ sample_id, floral_df_long, sum)
abundance_aggregate$floral_abundance = log10(abundance_aggregate$floral_abundance) #don't round this log value...yet
specimens = left_join(specimens, abundance_aggregate, by = "sample_id")

#calculate floral diversity at each sample_id
floral_df_wide %>%
  mutate(across(-c(sample_id), function(x) 10^x)) %>%
  mutate(across(-c(sample_id), ~replace(., is.na(.), 0))) -> floral_df_wide_exponential 
rownames(floral_df_wide_exponential) = floral_df_wide_exponential[,1]
floral_df_wide_exponential = floral_df_wide_exponential[,-1]
floral_df_wide_exponential$floral_diversity <- vegan::diversity(floral_df_wide_exponential)
floral_df_wide_exponential$sample_id = rownames(floral_df_wide_exponential)
floral_df = floral_df_wide_exponential[,colnames(floral_df_wide_exponential) %in% c("sample_id", "floral_diversity")]
specimens = left_join(specimens, floral_df, by = "sample_id")


#filter out only the bees that were screened for parasites
screenedbees = filter(specimens, gut_plate != "")
screenedbees = screenedbees[str_detect(screenedbees$gut_plate, "bcparasite_extra", negate = T),]

#filter out bees where apidae DNA didn't amplify
screenedbees = filter(screenedbees, apidae == 1)

#create parasite richness column
screenedbees$parasite_richness = screenedbees$apicystis + screenedbees$ascosphaera + screenedbees$nbombii + screenedbees$nceranae + screenedbees$cbombii + screenedbees$cexpoeki + screenedbees$crithidiaspp

#add a column for queen/worker
screenedbees$caste = ifelse(str_detect(screenedbees$notes, "queen"), "queen", "worker")

### START ANALYZING DATA!

#standardize richness and abundance values
screenedbees$bombus_richness_standardized = as.vector(scale(screenedbees$bombus_richness, center = TRUE, scale = TRUE))
screenedbees$bombus_abundance_standardized = as.vector(scale(screenedbees$bombus_abundance, center = TRUE, scale = TRUE))
screenedbees$floral_richness_standardized = as.vector(scale(screenedbees$point_floral_richness, center = TRUE, scale = TRUE))
screenedbees$floral_abundance_standardized = as.vector(scale(screenedbees$floral_abundance, center = TRUE, scale = TRUE))
screenedbees$impatiens_abundance_standardized = as.vector(scale(screenedbees$impatiens_abundance, center = TRUE, scale = TRUE))
screenedbees$julian_date_standardized = as.vector(scale(screenedbees$julian_date, center = TRUE, scale = TRUE))


#any parasite
screenedbees$anyparasite = screenedbees$apicystis + screenedbees$ascosphaera + 
  screenedbees$cbombii + screenedbees$cexpoeki + screenedbees$crithidiaspp +
  screenedbees$nbombii + screenedbees$nceranae
screenedbees$anyparasite[screenedbees$anyparasite > 0] <- 1

# a little data exploration
library(ggplot2)

ggplot(screenedbees,aes(x = final_id,y = parasite_richness)) +
  geom_violin() +
  theme_bw()

ggplot(screenedbees,aes(x = caste,y = parasite_richness)) +
  geom_violin() +
  theme_bw()

#make barplot of parasite richness for different castes
# Calculate total count of individuals in each caste
total_count <- screenedbees %>%
  group_by(caste) %>%
  summarize(total_count = n())

# Generate the summary table with counts of individuals in each caste and richness group
summary_table <- screenedbees %>%
  group_by(caste, parasite_richness) %>%
  summarize(n = n())

# Merge the total counts with the summary table
summary_table <- merge(summary_table, total_count, by = "caste")

# Calculate proportion of individuals in each richness group for each caste
summary_table <- summary_table %>%
  mutate(proportion = n / total_count)

# Plot the bar plot
ggplot(summary_table, aes(x = parasite_richness, y = proportion, fill = caste)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Parasite Richness", y = "Proportion", title = "Proportion of Individuals in Each Richness Group by Caste")




#make a heatmap of parasites x species
screenedbees %>%
  group_by(final_id) %>%
  summarise(apicystis = sum(apicystis),
            ascosphaera = sum(ascosphaera),
            nbombii = sum(nbombii),
            nceranae = sum(nceranae),
            cbombii = sum(cbombii),
            cexpoeki = sum(cexpoeki),
            crithidiaspp = sum(crithidiaspp)) -> parasite_sums

#normalize parasite prevalence by dividing by the number of each species screened
screenedbees %>%
  group_by(final_id) %>%
  summarize(n=n()) -> species_sums
parasite_sums$apicystis = parasite_sums$apicystis/species_sums$n
parasite_sums$ascosphaera = parasite_sums$ascosphaera/species_sums$n
parasite_sums$nbombii = parasite_sums$nbombii/species_sums$n
parasite_sums$nceranae = parasite_sums$nceranae/species_sums$n
parasite_sums$cbombii = parasite_sums$cbombii/species_sums$n
parasite_sums$cexpoeki = parasite_sums$cexpoeki/species_sums$n
parasite_sums$crithidiaspp = parasite_sums$crithidiaspp/species_sums$n

#heatmap!
parasite_sums = as.data.frame(parasite_sums)
rownames(parasite_sums) = parasite_sums[,1]
parasite_sums = parasite_sums[,-1]
parasite_sums = as.matrix(parasite_sums)
colnames(parasite_sums) = c("Apicystis spp.", "Ascosphaera spp.", "Nosema bombi", "Nosema ceranae", "Crithidia bombi", "Crithidia expoeki", "Crithidia spp.")

heatmap(parasite_sums)
levelplot(parasite_sums,
          xlab = list(label = "Bumble bee species", cex = 1.2), 
          ylab = list(label = "Parasite species/genera", cex = 1.2),
          scales= list(x= list(rot=90, cex = 1.2), y = list(cex = 1.2)))


scale=list(x=list(rot=90))


#normalize each parasite by dividing by number of times that parasite was detected overall
apicystis_sum = sum(screenedbees$apicystis)
ascosphaera_sum = sum(screenedbees$ascosphaera)
nbombii_sum = sum(screenedbees$nbombii)
nceranae_sum = sum(screenedbees$nceranae)
cbombii_sum = sum(screenedbees$cbombii)
cexpoeki_sum = sum(screenedbees$cexpoeki)
crithidiaspp_sum = sum(screenedbees$crithidiaspp)

parasite_sums = as.data.frame(parasite_sums)
parasite_sums$apicystis = parasite_sums$apicystis/apicystis_sum
parasite_sums$ascosphaera = parasite_sums$ascosphaera/ascosphaera_sum
parasite_sums$nbombii = parasite_sums$nbombii/nbombii_sum
parasite_sums$nceranae = parasite_sums$nceranae/nceranae_sum
parasite_sums$cbombii = parasite_sums$cbombii/cbombii_sum
parasite_sums$cexpoeki = parasite_sums$cexpoeki/cexpoeki_sum
parasite_sums$crithidiaspp = parasite_sums$crithidiaspp/crithidiaspp_sum

parasite_sums = as.matrix(parasite_sums)
heatmap(parasite_sums)

##make boxplots like Lauren's--average prevalence of each parasite at each subsite x period combo
fvimp_brmsdf <- read.csv("/Users/jenna1/Documents/UBC/Bombus Project/fvimpatiens_parasites/data/fvimp_brmsdf.csv", sep = ",", header = T, row.names = 1)
fvimp_brmsdf %>% mutate(period = case_when(round < 3 ~ "1",
                                           round > 2 & round < 5 ~ "2",
                                           round > 4 & round < 7 ~ "3",
                                           round > 6 & round < 9 ~ "4",
                                           round > 8 ~ "5")) -> fvimp_brmsdf

library(cowplot)
fvimp_brmsdf %>%
  filter(caste != "queen") %>%
  filter(apidae == 1) %>%
  filter(subsite != "") %>%
  group_by(subsite, period, final_id) %>%
  summarize(
    crith_prop = sum(crithidiaspp)/n(),
    ce_prop = sum(cexpoeki)/n(),
    cb_prop = sum(cbombii)/n(),
    nb_prop = sum(nbombii)/n(),
    asco_prop = sum(ascosphaera)/n(),
    api_prop = sum(apicystis)/n(),
    ncer_prop = sum(nceranae)/n()
  ) -> prevalence_df
truncated = prevalence_df[,3:ncol(prevalence_df)]

library(tidyr)
data_long <- gather(truncated, key = "Variable", value = "Proportion", -final_id)
data_long$Variable <- factor(data_long$Variable, labels = c("Apicystis spp.", "Ascosphaera spp.", "Crithidia bombi", "Crithidia expoeki", "Crithidia spp", "Vairimorpha bombi", "Vairimorpha ceranae"))


ggplot(data_long, aes(x = Proportion, y = Variable, fill = Variable)) +
  geom_boxplot() +
  facet_grid(.~ final_id, scales = "free_y") +
  labs(x = "Rate of Infection", y = "Parasite Group") +
  theme_minimal() +
  theme(legend.position = "none", panel.spacing = unit(1, "lines"))



library(ggpubr)
library(broom)
library(betareg)

# Calculate sample size for each species
data_long <- data_long %>%
  group_by(final_id) %>%
  mutate(species_with_n = paste0(final_id, "\n(n=", n(), ")"))

# Plot with updated y-axis labels
p <- ggplot(data_long, aes(x = Proportion, y = species_with_n, fill = final_id)) +
  geom_boxplot() +
  facet_grid(.~ Variable, scales = "free_y") +
  labs(x = "Rate", y = "Bumble bee species") +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(5, "lines"),
        text = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) # Adjust text size for y-axis labels if needed


ggsave("boxplot_parasiteprev.png", plot = p, width = 40, height = 8, dpi = 500)





# Load the multcomp package
library(multcomp)
screenedbees = fvimp_brmsdf %>% filter(apidae ==1 & caste == "worker")
# Fit a generalized linear model (GLM) to count data
screenedbees$final_id = as.factor(screenedbees$final_id)
model <- glm(nbombii ~ final_id, data = screenedbees, family = binomial)
summary(model)

# Conduct post-hoc tests for pairwise comparisons using Tukey's method
# Adjust the 'method' argument as needed (e.g., "Tukey", "Bonferroni", "Holm", etc.)
comparison <- glht(model, linfct = mcp(final_id = "Tukey"))

# Summarize the results
summary(comparison)



screenedbees %>% 
  group_by(caste) %>% 
  summarize(ApicystisSpp = sum(apicystis)/n(),
            AscosphaeraSpp = sum(ascosphaera)/n(),
            CrithidiaBombi = sum(cbombii)/n(),
            CrithidiaExpoeki = sum(cexpoeki)/n(),
            NosemaBombi = sum(nbombii)/n(),
            NosemaCeranae = sum(nceranae)/n()) -> caste_proportions
caste_long = gather(caste_proportions, key = "Parasite", value = "Proportion", -caste)


p <- ggplot(caste_long, aes(x = Parasite, y = Proportion, fill = caste)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Parasite Group", y = "Rate") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1))
ggsave("barplot_parasiteprev.png", plot = p, dpi = 500)







#run multinomial glm for parasite richness
screenedbees$site = as.factor(screenedbees$site)
screenedbees$subsite = as.factor(screenedbees$subsite)

trials = rep(7,length(screenedbees$parasite_richness))
glm.richness <- glmer(cbind(parasite_richness, 7-parasite_richness) ~ bombus_richness_standardized +
                 floral_richness_standardized + floral_abundance_standardized + 
                 final_id + bombus_abundance_standardized +
                 julian_date_standardized*impatiens_abundance_standardized + julian_date_standardized^2 +
                 (1|site/subsite) , 
               data = screenedbees, 
               family = binomial(),
               weights = trials)

#plot partial regressions for significant (or almost significant) factors
plot_model(glm.richness, type = "pred", terms = "final_id")
plot_model(glm.richness, type = "pred", terms = "caste")
plot_model(glm.richness, type = "pred", terms = "julian_date_standardized")
plot_model(glm.richness, type = "pred", terms = "bombus_richness_standardized")
plot_model(glm.richness, type = "pred", terms = "floral_abundance_standardized [all]")
plot_model(glm.richness, type = "int", terms = "julian_date_standardized:impatiens_abundance_standardized")

#run binomial glms for individual parasites

#works!
glm.nb <- glmer(nbombii ~ bombus_richness_standardized + caste +
                 floral_richness_standardized + floral_abundance_standardized + 
                 final_id + bombus_abundance_standardized +
                 impatiens_abundance_standardized + julian_date_standardized + julian_date_standardized^2 +
                 (1|site) , 
               data = screenedbees, 
               family = binomial())
plot_model(glm.nb, type = "pred", terms = "final_id")
plot_model(glm.nb, type = "pred", terms = "caste")

#works!
glm.api <- glmer(apicystis ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + caste + bombus_abundance_standardized +
                  julian_date_standardized*impatiens_abundance_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())
plot_model(glm.api, type = "pred", terms = "floral_abundance_standardized")
plot_model(glm.api, type = "pred", terms = "julian_date_standardized")
plot_model(glm.api, type = "int", terms = "julian_date_standardized:impatiens_abundance_standardized")



#singular fit:
glm.nc <- glmer(nceranae ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + bombus_abundance_standardized +
                  impatiens_abundance_standardized + julian_date_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())
#doesn't converge:
glm.cb <- glmer(cbombii ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + bombus_abundance_standardized +
                  impatiens_abundance_standardized + julian_date_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())

#large eigenvalue, nearly unidentifiable
glm.ce <- glmer(cexpoeki ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + bombus_abundance_standardized +
                  impatiens_abundance_standardized + julian_date_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())

#singular fit:
glm.crith <- glmer(crithidiaspp ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + bombus_abundance_standardized +
                  impatiens_abundance_standardized*julian_date_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())


#is singular
glm.asco <- glmer(ascosphaera ~ bombus_richness_standardized +
                  floral_richness_standardized + floral_abundance_standardized + 
                  final_id + bombus_abundance_standardized +
                  impatiens_abundance_standardized*julian_date_standardized + julian_date_standardized^2 +
                  (1|site) , 
                data = screenedbees, 
                family = binomial())


vars = data.frame(
  screenedbees$bombus_abundance_standardized,
  screenedbees$bombus_richness_standardized,
  screenedbees$floral_abundance_standardized,
  screenedbees$floral_richness_standardized,
  screenedbees$julian_date_standardized,
  screenedbees$impatiens_abundance_standardized
  )
corr_matrix = cor(vars, use = "pairwise.complete.obs")


my_data <- data.frame(
  var1 = c(1, 2, 3, 4, 5),
  var2 = c(2, 4, 6, 8, 10),
  var3 = c(3, 6, 9, 12, 15)
)









####CODE FOR IF YOU WANT TO LOOK AT COMMUNITY COMP / FLORAL COMPparasite_sums$apicystis = parasite_sums$apicystis/species_sums$n AT THE SUBSITE X TIMEPOINT LEVEL

#add Julian dates to specimen dataframe
datesdf = sampleEffort[,colnames(sampleEffort) %in% c("sample_id", "julian_date")]
specimens = left_join(cleanedspecimens, datesdf, by = "sample_id")

#calculate sampling effort at group level
sampleEffort$round = as.numeric(sampleEffort$round)
sampleEffort %>% mutate(period = case_when(round < 3 ~ "1",
                            round > 2 & round < 5 ~ "2",
                            round > 4 & round < 7 ~ "3",
                            round > 6 & round < 9 ~ "4",
                            round > 8 ~ "5")) %>%
  relocate(period, .after = round) -> sampleEffort

#reduce sample effort to just include points in a subsite
sampleEffort %>%
  filter(subsite != "") -> subsiteEffort

#calculate sampling effort
subsiteEffort$group = paste(sampleEffort$subsite, sampleEffort$period)
subsiteEffort %>% 
  group_by(group) %>%
  summarize(sampling_minutes=n()) -> groupeffort
groupeffort = groupeffort[colnames(groupeffort) %in% c("group", "sampling_minutes")]

#filter specimens dataframe down to just subsite specimens
specimens %>%
  filter(subsite != "") -> subsiteSpecimens

#add sampling effort to specimens dataframe
subsiteSpecimens$group = paste(subsiteSpecimens$subsite, subsiteSpecimens$period)
subsiteSpecimens = left_join(subsiteSpecimens, groupeffort, by = "group")
subsiteSpecimens %>%
  relocate(group, .after = round) %>%
  relocate(sampling_minutes, .after = round) -> subsiteSpecimens
subsiteSpecimens$sampling_minutes = subsiteSpecimens$sampling_minutes*5

#calculate species richness at the group level

#calculate richness of entire site, then divide by sampling minutes
richnessdf = group_by(subsiteSpecimens, group) %>%
  count(final_id) %>%
  count(name = "bombus_richness")
subsiteSpecimens = left_join(subsiteSpecimens, richnessdf, by = "group")
subsiteSpecimens$bombus_richness = subsiteSpecimens$bombus_richness/subsiteSpecimens$sampling_minutes

#add bee diversity (Shannon Index) and bombus abundance at the group level
#bees are POOLED for sampling points and then we take diversity of total pool
community_df_long = group_by(subsiteSpecimens, group) %>%
  count(final_id)
community_df_wide = as.data.frame(pivot_wider(community_df_long, names_from = "final_id", values_from = "n"))
community_df_wide[is.na(community_df_wide)] <- 0
rownames(community_df_wide) = community_df_wide[,1]
community_df = community_df_wide[,-1]
community_df_wide$H <- diversity(community_df)
colnames(community_df_wide) = c("group", "californicus_abundance", "flavifrons_abunance", "impatiens_abundance", "mixtus_abundance", "rufocinctus_abundance", "vosnesenskii_abundance", "melanopygus_abundance", "shannon_diversity")
community_df_wide$bombus_abundance = community_df_wide$californicus_abundance + community_df_wide$flavifrons_abunance + community_df_wide$impatiens_abundance + community_df_wide$mixtus_abundance + community_df_wide$rufocinctus_abundance + community_df_wide$vosnesenskii_abundance + community_df_wide$melanopygus_abundance

#add diversity and abundance to specimens dataframe, then normalize by dividing by sampling effort
subsiteSpecimens = left_join(subsiteSpecimens, community_df_wide, by = "group")
subsiteSpecimens$shannon_diversity = subsiteSpecimens$shannon_diversity/subsiteSpecimens$sampling_minutes
subsiteSpecimens$bombus_abundance = subsiteSpecimens$bombus_abundance/subsiteSpecimens$sampling_minutes

#calculate average diversity (averaged across points that fall into group)
# make sure to start with subsiteEffort dataframe, not specimens
#
#??? is this necessary?
#
#
#
#
#
#
#
#
#


#create a list of all the flowers visited by bees
beeflowers = unique(specimens$active_flower)
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu", negate = T)]

#create a floral community dataframe at the group level
columnlist = c("group", beeflowers)
group_to_point = subsiteEffort[,colnames(subsiteEffort) %in% c("group", "sample_id")]
vegData = left_join(vegData, group_to_point, by = "sample_id")
vegData %>% filter(!is.nan(group)) -> subsiteVeg
floral_df_wide = subsiteVeg[,colnames(subsiteVeg) %in% columnlist]
floral_df_wide = mutate(floral_df_wide, across(-c(group), as.numeric))

#calculate floral richness
#all flowers are POOLED and we take richness of entire pool
floral_df_wide %>%
  na.replace(0) %>%
  group_by(group) %>% 
  summarise(across(where(is.numeric), max)) %>%
  mutate_at(vars(-1), ~ +(.x > 0)) -> floral_grouped_binary
floral_grouped_binary$floral_richness = rowSums(floral_grouped_binary[-1])
floral_df_richness = floral_grouped_binary[,c(1,ncol(floral_grouped_binary))]
subsiteSpecimens = left_join(subsiteSpecimens, floral_df_richness, by = 'group')
subsiteSpecimens$floral_richness = subsiteSpecimens$floral_richness/subsiteSpecimens$sampling_minutes

##calculate average floral richness per point and sum
#columnlist = c("sample_id", beeflowers)
#floral_df_wide = subsiteVeg[,colnames(subsiteVeg) %in% columnlist]
#floral_df_wide = mutate(floral_df_wide, across(-c(sample_id), as.numeric))
#floral_df_wide %>%
#  mutate_at(vars(-1), ~ +(.x > 0)) -> floral_point_binary
#floral_point_binary$point_floral_richness = rowSums(floral_point_binary[-1], na.rm = TRUE)
#point_richness = floral_point_binary[,c(1,ncol(floral_point_binary))]
#point_richness = left_join(point_richness, group_to_point, by = "sample_id")
#average_richness = aggregate(point_floral_richness ~ group, point_richness, mean)
#colnames(average_richness) = c("group", "average_floral_richness")
#specimens = left_join(specimens, average_richness, by = 'group')


#calculate floral abundance
#create a floral community dataframe at the sampling point level (so we can average over points)
columnlist = c("group", beeflowers)
floral_df_wide = subsiteVeg[,colnames(subsiteVeg) %in% columnlist]
floral_df_wide = mutate(floral_df_wide, across(-c(group), as.numeric))

#sum across all sampling points in group then divide by sampling effort
floral_df_wide %>%
  pivot_longer(!group, names_to = "flower", values_to = "floral_abundance") -> floral_df_long
floral_df_long$floral_abundance = 10^(floral_df_long$floral_abundance - 1)
abundance_aggregate = aggregate(floral_abundance ~ group, floral_df_long, sum)
subsiteSpecimens = left_join(subsiteSpecimens, abundance_aggregate, by = "group")
subsiteSpecimens$floral_abundance = subsiteSpecimens$floral_abundance/(subsiteSpecimens$sampling_minutes/5)
subsiteSpecimens$floral_abundance = log10(subsiteSpecimens$floral_abundance)

# #calculate floral abundance by getting abundance of each point (log val) and averaging across points
# #create a floral community dataframe at the sampling point level (so we can average over points)
# columnlist = c("sample_id", beeflowers)
# floral_df_wide = specimens[,colnames(specimens) %in% columnlist]
# floral_df_wide = mutate(floral_df_wide, across(-c(sample_id), as.numeric))
# 
# #we calculate abundance at each sampling point and average across points
# floral_df_wide %>%
#   pivot_longer(!sample_id, names_to = "flower", values_to = "floral_abundance") -> floral_df_long
# floral_df_long$floral_abundance = 10^(floral_df_long$floral_abundance - 1)
# abundance_aggregate = aggregate(floral_abundance ~ sample_id, floral_df_long, sum)
# abundance_aggregate$floral_abundance = log10(abundance_aggregate$floral_abundance) #don't round this log value...yet
# 
# #join 'group' variable to abundance_aggregate so you can group the sample points accordingly
# group_to_point = subsiteEffort[,colnames(subsiteEffort) %in% c("group", "sample_id")]
# abundance_aggregate = left_join(abundance_aggregate, group_to_point, by = "sample_id")
# abundance_aggregate %>% 
#   group_by(group) %>%
#   summarise(across(where(is.numeric), mean)) -> average_abundance
# colnames(average_abundance) = c("group", "average_floral_abundance")
# 
# #now calculate total floral abundance by summing across all sampling points in group
# abundance_aggregate$floral_abundance = 10^abundance_aggregate$floral_abundance
# total_abundance = aggregate(floral_abundance ~ group, abundance_aggregate, sum)
# total_abundance$floral_abundance = log10(total_abundance$floral_abundance)
# colnames(total_abundance) = c("group", "total_floral_abundance")
# abundance_df = left_join(average_abundance, total_abundance, by = "group")
# specimens = left_join(specimens, abundance_df, by = "group")

#calculate floral diversity
#POOL floral communities in group and calculate overall diversity
floral_df_wide %>%
  mutate(across(-c(group), function(x) 10^x)) %>%
  mutate(across(-c(group), ~replace(., is.na(.), 0))) -> floral_df_wide_exponential 

floral_df_aggregated <- as.data.frame(floral_df_wide_exponential %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), sum)))
floral_df_aggregated = floral_df_aggregated[1:120,]
rownames(floral_df_aggregated) = floral_df_aggregated[,1]
floral_df_aggregated = floral_df_aggregated[,-1]
floral_df_aggregated$floral_diversity <- diversity(floral_df_aggregated)
floral_df_aggregated$group = rownames(floral_df_aggregated)
floral_df = floral_df_aggregated[,colnames(floral_df_aggregated) %in% c("group", "floral_diversity")]

#add floral diversity to specimens dataframe and normalize by dividing by sampling minutes
subsiteSpecimens = left_join(subsiteSpecimens, floral_df, by = "group")
subsiteSpecimens$floral_diversity = subsiteSpecimens$floral_diversity/subsiteSpecimens$sampling_minutes


#calculate floral diversity
#calculate average floral diversity across points
#
#
#
# is this necessary???
#
#
#
#
#

#create a df which includes just the bees that have been screened for parasites
screenedbees = filter(subsiteSpecimens, gut_plate != "")
screenedbees = screenedbees[str_detect(screenedbees$gut_plate, "bcparasite_extra", negate = T),]


#####END HERE FOR GROUP LEVEL
