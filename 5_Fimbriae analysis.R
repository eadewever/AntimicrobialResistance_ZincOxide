#This script aims at comparing the fimbriae expressed in the pre- vs post-withdrawal groups

# ctrl + alt + R to run the whole script

#Set working directory
setwd("indicate working directory")


#Install necessary packages, doing it if and only if not already installed
packages <- c("ggplot2", "tidyr", "dplyr")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

#Loading packages
library(ggplot2)
library(tidyr)
library(dplyr)


#Creating a folder to save the graphs
if (!dir.exists("Graphs")) {
  dir.create("Graphs")}


#Closing all previous plots
if (dev.cur() != 1) {
  dev.off()
}






#This script has been designed to work with a .csv file organised as follows:
#One row = one isolate/bacteria
#Column 1 = Farm = ID of the farm from which the isolate comes from
#Column 2 = Visit = sampling time point (V1, V2 or V4 in our case)
#Column 3 = Group = pre-withdrawal or post-withdrawal (respectively, received or did not receive zinc oxide in our case)
#Column 4 = Sample = ID of the sample within each sampling time point (e.g. W01, W02, etc. in our case)
#Column 5 = Agar = to distinguish isolates collected on Unselective or Selective agar
#Column 6 = Pick = individual number of the isolate within each sample (e.g. P1, P2, etc in our case)

#Columns 7 to 20 = AST values (zone diameter in millimeters)
#Antibiotics should be in the following order: Tetracycline, Ampicillin, Amoxicillin_Clavulanate,
#Apramycin, Gentamicin, Streptomycin, Spectinomycin, Sulfamethoxazole_Trimethoprim,
#Florfenicol, Chloramphenicol, Ceftiofur, Cefotaxime, Enrofloxacin, Ciprofloxacin

#Columns 21 to 29 = results of PCR analysis (if applicable)
#Columns 21 to 23 = toxins (Yes if toxin detected by PCR, No otherwise) (e.g. STa, LTa, Stx2e)
#Columns 24 to 29 = fimbriae (Yes if fimbriae detected by PCR, No otherwise) (e.g. F4, F18, etc)


#Loading the Antimicrobial Susceptibility (AST) data
AST_data <- read.csv("add_filename.csv",
                     header = TRUE,
                     sep = ",",
                     na.strings = "") #If empty values = analysis not done, then add NA


#Filtering AST_data to remove isolates with missing info about AMR profiles (i.e. at least one NA in columns of ATB)
AST_data_filtered <- AST_data[!apply(is.na(AST_data[, 7:20]), 1, any), ]


#Retrieving the list of toxins, fimbriae and visits
toxins <- colnames(AST_data_filtered)[21:22] #Only STa and LTa as no Stx2e positive isolates
fimbriae <- colnames(AST_data_filtered)[24:25] #Only F4 and F18 as other fimbriae not detected
visits <- unique(AST_data_filtered$Visit)


#Adding a column to AST_data_filtered with one unique ID per isolate (removing Farm ID for confidentiality reasons)
AST_data_filtered$IsolateID <- paste(AST_data_filtered$Visit,
                                     AST_data_filtered$Sample,
                                     AST_data_filtered$Pick,
                                     sep = "_")


#Creating a new dataframe with only toxins and fimbriae
data_ToxFim <- AST_data_filtered %>%
  #Filtering to keep isolates expressing either LTa or STa (removing isolates without toxins)
  filter(if_any(all_of(toxins), ~ . == "Yes")) %>%
  #Selecting only isolateID, toxins and fimbriae (i.e. not keeping AST values)
  select(IsolateID, Visit, all_of(toxins), all_of(fimbriae)) %>%
  #Adding a column to identify isolates expressing no fimbriae (Yes if does not express any fimbriae, no otherwise)
  mutate(NoFimbriae = if_else(F4 == "No" & F18 == "No", "Yes", "No"))


#Creating a long version of data_ToxFim
data_ToxFim_long <- data_ToxFim %>%
  #Transforming data based on toxins
  pivot_longer(cols = toxins,
               names_to = "Toxin", values_to = "PresenceToxin") %>%
  #Filtering to keep only toxins expressed (removes combination isolate + toxin if isolate does not express toxin)
  #Because we are only looking at ETEC isolates = must express at least one toxin
  filter(PresenceToxin != "No") %>%
  
  #Transforming data based on fimbriae (F4, F18, or NoFimbriae when no fimbriae detected through PCR)
  pivot_longer(cols = c(fimbriae, "NoFimbriae"),
               names_to = "Fimbriae", values_to = "PresenceFimbriae") %>%
  #Filtering to keep only rows if PresenceFimbriae = Yes
  #As backup since all the isolates should be assigned either F4, F18 or NoFimbriae
  filter(PresenceFimbriae != "No") %>%
  
  #Removing columns PresenceToxin and PresenceFimbriae (because all = "Yes" and do not need these columns afterwards)
  select(-PresenceToxin, -PresenceFimbriae)


#As of now, we are only considering STa- and LTa-positive isolates = ETEC
#For each visit, we are going to calculate how many isolate express each fimbriae (or no fimbriae)
#And determine the percentage by dividing by the total number of STa- and LTa-positive isolates
data_ToxFim_percent <- data_ToxFim_long %>%
  #Grouping by visit to calculate number of STa- and LTa-positive isolates
  group_by(Visit) %>%
  mutate(Total_isolates = n_distinct(IsolateID)) %>%
  
  #Grouping by visit, fimbriae and toxin to calculate nb of isolates for each fimbriae and toxin (per visit)
  group_by(Visit, Fimbriae, Toxin) %>%
  summarise(Total_isolates = first(Total_isolates),
            Nb_isolates = n(),
            .groups = "drop") %>%
  #Calculating corresponding percentage (dividing by total number of STa- and LTa-positive isolates)
  mutate(Percentage_isolates = 100 * Nb_isolates / Total_isolates)





#Now, creating different graphs, separated by fimbriae or visit depending on the chosen focus

#Creating the graph, separating by fimbriae
fimbriae_barplot_perFim <- ggplot(data_ToxFim_percent,
                           aes(x = Visit, y = Percentage_isolates, fill = Toxin)) +
  
  geom_bar(stat = "identity", position = "stack") +
  
  theme_bw() +
  
  facet_wrap(~ Fimbriae,
             labeller = as_labeller(c("F4" = "F4",
                                      "F18" = "F18",
                                      "NoFimbriae" = "No fimbriae"))) +
  
  scale_fill_manual(values = c("LTa" = "powderblue", "STa" = "cadetblue")) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 5), # forcing R to display from 0% to 100% on the y axis
                     expand = c(0, 1)) +
  
  labs(x = "Visit", y = "Percentage of toxin-positive isolates (%)", fill = "Toxin") +
  
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14))


#Saving the barplot in the "Graphs" folder
png("Graphs/fimbriae_barplot_perFim.png", width = 3000, height = 2000, res = 300)
print(fimbriae_barplot_perFim)
dev.off()





#Creating the graph, separating by visit
fimbriae_barplot_perVisit <- ggplot(data_ToxFim_percent,
                                  aes(x = Fimbriae, y = Percentage_isolates, fill = Toxin)) +
  
  geom_bar(stat = "identity", position = "stack") +
  
  theme_bw() +
  
  facet_wrap(~ Visit) +
  
  scale_fill_manual(values = c("LTa" = "powderblue", "STa" = "cadetblue")) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 5), # forcing R to display from 0% to 100% on the y axis
                     expand = c(0, 1)) +
  
  scale_x_discrete(labels = c("F4" = "F4",
                              "F18" = "F18",
                              "NoFimbriae" = "No fimbriae")) +
  
  labs(x = "Fimbriae", y = "Percentage of toxin-positive isolates (%)", fill = "Toxin") +
  
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14))


#Saving the barplot in the "Graphs" folder
png("Graphs/fimbriae_barplot_perVisit.png", width = 3000, height = 2000, res = 300)
print(fimbriae_barplot_perVisit)
dev.off()





#Creating the graph, separating by toxin
fimbriae_barplot_perToxin <- ggplot(data_ToxFim_percent,
                                    aes(x = Toxin, y = Percentage_isolates, fill = Fimbriae)) +
  
  geom_bar(stat = "identity", position = "stack") +
  
  theme_bw() +
  
  facet_wrap(~ Visit) +
  
  scale_fill_manual(values = c("F4" = "peachpuff",
                               "F18" = "sandybrown",
                               "NoFimbriae" = "darkgrey")) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 5), # forcing R to display from 0% to 100% on the y axis
                     expand = c(0, 1)) +
  
  labs(x = "Toxin", y = "Percentage of toxin-positive isolates (%)", fill = "Fimbriae") +
  
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14))


#Saving the barplot in the "Graphs" folder
png("Graphs/fimbriae_barplot_perToxin.png", width = 3000, height = 2000, res = 300)
print(fimbriae_barplot_perToxin)
dev.off()





#Closing all plots
if (dev.cur() != 1) {
  dev.off()
}