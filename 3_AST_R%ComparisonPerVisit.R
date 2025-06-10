#This script aims at comparing the resistance profiles in the pre- vs post-withdrawal groups by plotting the resistance prevalence across visits as boxplots

# ctrl + alt + R to run the whole script

#Set working directory
setwd("indicate working directory")


#Install necessary packages, doing it if and only if not already installed
packages <- c("ggplot2", "tidyr", "dplyr", "pheatmap", "tibble", "purrr", "writexl", "rstatix", "ggpubr", "FSA")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

#Loading packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(tibble)
library(purrr)
library(writexl)
library(rstatix)
library(ggpubr)
library(FSA)


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


#Filtering AST_data to remove isolates with missing AST values (i.e. at least one NA in antibiotics columns = columns 7 to 20)
AST_data_filtered <- AST_data[!apply(is.na(AST_data[, 7:20]), 1, any), ]


#Adding a column with the sample ID (i.e. Farn_V1_W01)
AST_data_filtered$SampleID <- paste(AST_data_filtered$Farm,
                                     AST_data_filtered$Visit,
                                     AST_data_filtered$Sample,
                                     sep = "_")

#Adding a column with the full isolate ID (i.e. Farm_V1_W01_P1)
AST_data_filtered$IsolateID <- paste(AST_data_filtered$Farm,
                                     AST_data_filtered$Visit,
                                     AST_data_filtered$Sample,
                                     AST_data_filtered$Pick,
                                     sep = "_")


#Retrieving the list of antibiotics, visits and groups (pre-/post-withdrawal)
antibiotics <- colnames(AST_data_filtered)[7:20]
visits <- unique(AST_data_filtered$Visit)
groups <- unique(AST_data_filtered$Group)


#Creating a table with, for each antibiotic, the threshold zone diameters defining resistant and susceptible isolates
#Epidemiological breakpoint values according to EUCAST 2024 and VARSS 2024 (values listed in the order of the antibiotics aforementioned)
#Diameters defining resistant isolates (resistant if <= list_R_diameter_adapted value)
list_R_diameter <- c(11,13,18,11,16,10,17,10,18,16,17,16,16,21)
#Diameters defining susceptible isolates (susceptible if >= list_S_diameter_adapted value)
list_S_diameter <- c(15,14,19,15,17,15,21,14,19,17,21,20,23,25)
#If list_R_diameter != list_S_diameter - 1 = an "intermediary" category exists, i.e. isolate considered as susceptible, increased exposure
table_RS_ATB <- data.frame(Antibiotic = antibiotics, R = list_R_diameter, S = list_S_diameter)


#Creating a table replacing zone diameter values by R (Resistant), I (Susceptible, Increased exposure) or S (Susceptible, Standard regimen) for each isolate and each antibiotic
RIS_data <- AST_data_filtered[, c(1:6, ncol(AST_data_filtered) - 1, ncol(AST_data_filtered))]

for (antibiotic in antibiotics) {
  for (isolate in RIS_data$IsolateID) {
    
    zone_diameter <- AST_data_filtered[AST_data_filtered$IsolateID == isolate, antibiotic]
    R_threshold <- table_RS_ATB[table_RS_ATB$Antibiotic == antibiotic, "R"][[1]]
    S_threshold <- table_RS_ATB[table_RS_ATB$Antibiotic == antibiotic, "S"][[1]]
      
    if (zone_diameter <= R_threshold) {
      RIS_data[RIS_data$IsolateID == isolate, antibiotic] <- "R"
    }
    else if (zone_diameter >= S_threshold) {
      RIS_data[RIS_data$IsolateID == isolate, antibiotic] <- "S"
    }
    else {
      RIS_data[RIS_data$IsolateID == isolate, antibiotic] <- "I"
    }
  }
}


#Replacing I by S in a new data frame
#We want to compare percentage of resistant vs susceptible bacteria, thus I isolates considered as S
RS_data <- RIS_data

RS_data[, antibiotics] <- lapply(RS_data[, antibiotics], function(x) replace(x, x == "I", "S"))


#Transforming the data set so that each antibiotic is displayed on rows and not columns
RS_data_long <- pivot_longer(
  RS_data,
  cols = c(all_of(antibiotics)),
  names_to = "Antibiotics",         
  values_to = "AMR_profile"            
)










#Transforming the data set to have the following columns = antibiotics, sample ID, visit (V1, V2, V4) and agar type (selective or unselective)
#And for each category indicated by the columns above, indicate number of R bacteria and %, total number of bacteria
RS_data_grouped_agar <- RS_data_long %>%
  group_by(Antibiotics, SampleID, Visit, Agar) %>%
  summarise(Number_R = sum(AMR_profile == "R"),
            Number_Total = n(),
            Percentage_R = 100 * Number_R / Number_Total,
            .groups = "drop")


#Computing a few statistic parameters for each antibiotic R rates (median and mean)
RS_data_stats_agar <- RS_data_grouped_agar %>%
  group_by(Antibiotics, Visit, Agar) %>%
  summarise(#Data not normally distributed = indicate median + IQR
            Median = median(Percentage_R),
            Min = min(Percentage_R),
            Max = max(Percentage_R),
            Q1 = quantile(Percentage_R, probs = c(0.25)),
            Q3 = quantile(Percentage_R, probs = c(0.75)),
            IQR = Q3 - Q1,
            .groups = "drop")





#We consider all isolates from one sample as linked but independent from the isolates from another sample
#Thus, we can perform a Kruskal-Wallis test (data not normally distributed) to compare the resistance percentages between the V1, V2 and V4 groups

#Renaming column names of the dataframe to give 3 letter code for each antibiotic
RS_data_grouped_agar_codeATB <- RS_data_grouped_agar %>%
  mutate(Antibiotics = recode(Antibiotics,
                              "Tetracycline" = "TET",
                              "Ampicillin" = "AMP",
                              "Amoxicillin_Clavulanate" = "AMC",
                              "Apramycin" = "APR",
                              "Gentamicin" = "GEN",
                              "Streptomycin" = "STR",
                              "Spectinomycin" = "SPE",
                              "Sulfamethoxazole_Trimethoprim" = "SUT",
                              "Florfenicol" = "FLO",
                              "Chloramphenicol" = "CHL",
                              "Ceftiofur" = "CFT",
                              "Cefotaxime" = "CTX",
                              "Enrofloxacin" = "ENR",
                              "Ciprofloxacin" = "CIP"))


#In order to display the antibiotics in the same order as on the heatmap (cf. clustering of 2_AST_Pheatmap4)
display_ATB <- c("CFT", "CIP", "CTX", "ENR", "SUT", "AMP", "TET",
                 "STR", "FLO", "CHL", "SPE", "APR", "AMC", "GEN")

RS_data_grouped_agar_codeATB$Antibiotics <- factor(RS_data_grouped_agar_codeATB$Antibiotics,
                                                   levels = display_ATB)


# Filtering to display only picks from unselective agar
RS_data_filtered_unselective <- RS_data_grouped_agar_codeATB %>%
  filter(Agar == "Unselective")


#Statistical tests to determine whether there is a significant difference between groups
#Kruskal-Wallis to determine whether significant difference between groups
#Dunn's test (BH correction) if Kruskal-Wallis significant 

#Taking into account only picks from unselective UTI
RS_results_unselective <- RS_data_filtered_unselective %>%
  group_by(Antibiotics) %>%
  do({
    #Kruskal-Wallis
    kruskal_test <- kruskal.test(Percentage_R ~ Visit, data = .)
    
    #Dunn's test to compare values between visits if Kruskal_Wallis significant (i.e. p-value < 0.05)
    if (kruskal_test$p.value < 0.05) {
      # BH method for adjusting value
      dunn_test <- dunnTest(Percentage_R ~ Visit, data = ., method = "bh")
      dunn_results <- dunn_test$res
    } else {
      dunn_results <- data.frame(Comparison = "N/A", Z = NA, P.unadj = NA, P.adj = NA)
    }
    
    # Add significance (i.e. ns, *, **, ***) to the Dunn's results using add_significance
    dunn_results <- dunn_results %>%
      add_significance("P.adj")
    
    # Table with all the results
    data.frame(Antibiotics = unique(.$Antibiotics),
               kruskal_p_value = kruskal_test$p.value,
               dunn_results)
  })


#Print results and saving as an excel file
print(RS_results_unselective)
write_xlsx(RS_results_unselective, "RS_results_unselective.xlsx")


#Third step = plotting the data as boxplots with one boxplot per antibiotic
RS_boxplots_unselective <- ggplot(RS_data_filtered_unselective,
                                  aes(x = Visit, y = Percentage_R, fill = Visit)) +
  
  geom_boxplot() +
  
  theme_bw() +
  
  scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
  
  facet_wrap(~ Antibiotics, ncol = 5, scales = "free_y") + # separating the boxplots with one per antibiotic
  
  scale_y_continuous(limits = c(0, 120),
                     breaks = seq(0, 100, by = 10)) + # forcing R to display from 0% to 100% on the y axis
  
  stat_summary(fun = mean, geom = "point",
               shape = 4, size = 2, color = "black", fill = "black") + #Showing the mean as a black cross
  
  labs(title = "Resistance Percentage by Antibiotic (unselective agar)",
       x = NULL,
       y = "Percentage of resistant isolates (%)") +
  
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text.x = element_blank(), # Removing x labels as already indicated by fill = Group
        axis.ticks.x = element_blank(), # removing ticks on x axis
        #Changing font size
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  
  #Changing legend overall size
  guides(fill = guide_legend(keywidth = 2.5, keyheight = 3))
        

#Displaying the boxplots
RS_boxplots_unselective


#Saving the boxplots in the "Graphs" folder
ggsave("Graphs/RS_boxplots_unselective.png", plot = RS_boxplots_unselective,
       width = 22, height = 15, dpi = 800)










#Doing the same as previously, this time by looking at all the isolates (from selective and unselective agar)

#Transforming the data set to have the following columns = antibiotics, sample ID, visit (V1, V2, V4)
#And for each category indicated by the columns above, indicate number of R bacteria and %, total number of bacteria
RS_data_grouped_all <- RS_data_long %>%
  group_by(Antibiotics, SampleID, Visit) %>%
  summarise(Number_R = sum(AMR_profile == "R"),
            Number_Total = n(),
            Percentage_R = 100 * Number_R / Number_Total,
            .groups = "drop")


#Computing a few statistic parameters for each antibiotic R rates (median and mean)
RS_data_stats_all <- RS_data_grouped_all %>%
  group_by(Antibiotics, Visit) %>%
  summarise(#In case data not normally distributed = indicate median + IQR
    Median = median(Percentage_R),
    Min = min(Percentage_R),
    Max = max(Percentage_R),
    Q1 = quantile(Percentage_R, probs = c(0.25)),
    Q3 = quantile(Percentage_R, probs = c(0.75)),
    IQR = Q3 - Q1,
    .groups = "drop")


#Renaming column names of the dataframe to give 3 letter code for each antibiotic
RS_data_grouped_all_codeATB <- RS_data_grouped_all %>%
  mutate(Antibiotics = recode(Antibiotics,
                              "Tetracycline" = "TET",
                              "Ampicillin" = "AMP",
                              "Amoxicillin_Clavulanate" = "AMC",
                              "Apramycin" = "APR",
                              "Gentamicin" = "GEN",
                              "Streptomycin" = "STR",
                              "Spectinomycin" = "SPE",
                              "Sulfamethoxazole_Trimethoprim" = "SUT",
                              "Florfenicol" = "FLO",
                              "Chloramphenicol" = "CHL",
                              "Ceftiofur" = "CFT",
                              "Cefotaxime" = "CTX",
                              "Enrofloxacin" = "ENR",
                              "Ciprofloxacin" = "CIP"))


#In order to display the antibiotics in the same order as on the heatmap (cf. 2_AST_Pheatmap4)
display_ATB <- c("CFT", "CIP", "CTX", "ENR", "SUT", "AMP", "TET",
                 "STR", "FLO", "CHL", "SPE", "APR", "AMC", "GEN")

RS_data_grouped_all_codeATB$Antibiotics <- factor(RS_data_grouped_all_codeATB$Antibiotics,
                                                  levels = display_ATB)


#Statistical tests to determine whether there is a significant difference between groups
#Kruskal-Wallis to determine whether significant difference between groups
#Dunn's test (BH correction) if Kruskal-Wallis significant 

#Taking into account all picks (from unselective and selective combined)
RS_results_all <- RS_data_grouped_all_codeATB %>%
  group_by(Antibiotics) %>%
  do({
    #Kruskal-Wallis
    kruskal_test <- kruskal.test(Percentage_R ~ Visit, data = .)
    
    #Dunn's test to compare values between visits if Kruskal-Wallis significant (i.e. p-value < 0.05)
    if (kruskal_test$p.value < 0.05) {
      # BH method for adjusting value
      dunn_test <- dunnTest(Percentage_R ~ Visit, data = ., method = "bh")
      dunn_results <- dunn_test$res
    } else {
      dunn_results <- data.frame(Comparison = "N/A", Z = NA, P.unadj = NA, P.adj = NA)
    }
    
    # Add significance to the Dunn's results using add_significance
    dunn_results <- dunn_results %>%
      rename(p.adj = P.adj) %>%  # The column has to be named "p.adj"
      add_significance("p.adj")
    
    # Table with all the results
    data.frame(Antibiotics = unique(.$Antibiotics),
               kruskal_p_value = kruskal_test$p.value,
               dunn_results)
  })

#Print results and saving as an excel file
print(RS_results_all)
write_xlsx(RS_results_all, "RS_results_all.xlsx")


#Third step = plotting the data as boxplots with one boxplot per antibiotic
RS_boxplots_all <- ggplot(RS_data_grouped_all_codeATB,
                          aes(x = Visit, y = Percentage_R, fill = Visit)) +
  
  geom_boxplot() +
  
  theme_bw() +
  
  scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
  
  facet_wrap(~ Antibiotics, ncol = 5, scales = "free_y") + # separating the boxplots with one per antibiotic
  
  scale_y_continuous(limits = c(0, 120),
                     breaks = seq(0, 100, by = 10)) + # forcing R to display from 0% to 100% on the y axis
  
  stat_summary(fun = mean, geom = "point",
               shape = 4, size = 2, color = "black", fill = "black") + #Showing the mean as a black cross
  
  labs(title = "Resistance Percentage by Antibiotic",
       x = NULL,
       y = "Percentage of resistant isolates (%)") +
  
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text.x = element_blank(), # Removing x labels as already indicated by fill = Group
        axis.ticks.x = element_blank(), # removing ticks on x axis
        #Changing font size
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  
  #Changing legend overall size
  guides(fill = guide_legend(keywidth = 2.5, keyheight = 3))


#Displaying the boxplots
RS_boxplots_all


#Saving the boxplots in the "Graphs" folder
ggsave("Graphs/RS_boxplots_all.png", plot = RS_boxplots_all,
       width = 22, height = 15, dpi = 300)