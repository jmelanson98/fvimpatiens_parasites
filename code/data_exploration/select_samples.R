yearOne = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022_samplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)

##################################################################################
###### How should we group the bees (# of timepoints)
##################################################################################

#if we use three timepoints and 4 subsites per site, summarize:
yearOneJoined = left_join(yearOne, samplePoints,  by = "sample_pt") %>%
  mutate(period = case_when(round < 4 ~ "1",
                            round > 3 & round < 7 ~ "2",
                            round > 6 ~ "3")) %>%
  relocate(period, .after = round) 
yearOneJoined = yearOneJoined[str_detect(yearOneJoined$notes, "released|male|lost", negate = T),]
yearOneJoined %>% filter(!is.na(subsite) & subsite != "") %>%
                  filter(final_id != "B. huntii" & final_id != "B. melanopygus" & final_id != "") %>%
                  group_by(subsite, final_id, period) %>%
                  summarise(n=n()) -> subset4summary
subset4summary %>% filter(n >= 5) -> threeroundgood

#with five timepoints and four subsites, summarize:
twoRoundGrouping = left_join(yearOne, samplePoints,  by = "sample_pt") %>%
  mutate(period = case_when(round < 3 ~ "1",
                            round > 2 & round < 5 ~ "2",
                            round > 4 & round < 7 ~ "3",
                            round > 6 & round < 9 ~ "4",
                            round > 8 ~ "5")) %>%
  relocate(period, .after = round) 
twoRoundGrouping = twoRoundGrouping[str_detect(yearOneJoined$notes, "released|male|lost", negate = T),]
twoRoundGrouping %>% filter(!is.na(subsite) & subsite != "") %>%
  filter(final_id != "B. huntii" & final_id != "B. melanopygus" & final_id != "") %>%
  group_by(subsite, final_id, period) %>%
  summarise(n=n()) -> tworoundsummary
tworoundsummary %>% filter(n >= 5) -> tworoundgood

#make histograms of how many points *work* (>=5 individuals) for each species
tworoundgood %>% group_by(final_id, period) %>%
  summarise(n=n()) -> tworoundforplot
threeroundgood %>% group_by(final_id, period) %>%
  summarise(n=n()) -> threeroundforplot

ggplot(tworoundforplot,aes(x = subsite,y = n, fill=final_id)) + 
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  labs(y = "number of subsite w/ 5+ bees", x = "timepoint", fill = "species", title = "five timepoints of two rounds each") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(threeroundforplot,aes(x = period,y = n, fill=final_id)) + 
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  labs(y = "number of subsite w/ 5+ bees", x = "timepoint", fill = "species", title = "three timepoints (3, 3, 4 rounds respectively)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# now plot how these are distributed spatially
tworoundgood %>% group_by(final_id, subsite) %>%
  summarise(n=n()) -> tworoundspatial

ggplot(tworoundspatial,aes(x = subsite,y = n, fill=final_id)) + 
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  labs(y = "number of timepoints w/ 5+ bees", x = "subsite", fill = "species", title = "five timepoints") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

threeroundgood %>% group_by(final_id, subsite) %>%
  summarise(n=n()) -> threeroundspatial

ggplot(threeroundspatial,aes(x = subsite,y = n, fill=final_id)) + 
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  labs(y = "number of timepoints w/ 5+ bees", x = "subsite", fill = "species", title = "three timepoints") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



##########################################################################################
#####Decided to use 5 timepoints, with each site split into 4 pieces (e.g., tworoundgood)
#Now must take a random subset of 5 bees from each site x species x time combo
##########################################################################################

twoRoundGrouping %>% filter(subsite != "") %>%
  group_by(final_id, subsite, period) %>%
  slice_sample(n=5) %>%
  filter(n() >= 5) -> parasitebees
write.csv(parasitebees, "parasitebees.csv")

parasitebees %>% group_by(final_id, subsite, period) %>% count() %>% View()



##########################################################################################
##### Create a csv subset of specimen data, which includes all gutted bees and associated info
##########################################################################################

#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,-2]
samplePoints = samplePoints[,-2] #run this twice to get rid of both columns lol

guttedbees = left_join(yearOne, samplePoints,  by = "sample_pt") %>%
  filter(gut_plate != "")
guttedbees = guttedbees[str_detect(guttedbees$gut_plate, "bcparasite_extra", negate = T),]
write.csv(guttedbees, "bcparasite.csv")
write.csv(t(guttedbees), "transposebees.csv")
