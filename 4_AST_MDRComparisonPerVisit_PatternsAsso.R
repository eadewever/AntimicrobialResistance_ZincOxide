#This script aims at comparing multi-drug resistance prevalence in the pre- vs post-withdrawal groups by plotting the MDR prevalence across visits as boxplots
#It also looks at frequent resistance combinations, both at the level of the antibiotics and the level of the classes

# ctrl + alt + R to run the whole script

#Set working directory
setwd("indicate working directory")


#Install necessary packages, doing it if and only if not already installed
packages <- c("ggplot2", "tidyr", "dplyr", "pheatmap", "tibble", "FSA", "writexl")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

#Loading packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(tibble)
library(FSA)
library(writexl)


#Creating a folder to save the graphs
if (!dir.exists("Graphs")) {
  dir.create("Graphs")}


#Closing all plots (just in case one remained open from previous run)
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


#Adding a column with the sample ID (i.e. V1_W01)
AST_data_filtered$SampleID <- paste(AST_data_filtered$Visit,
                                    AST_data_filtered$Sample,
                                    sep = "_")

#Adding a column with the full isolate ID (i.e. V1_W01_P1)
AST_data_filtered$IsolateID <- paste(AST_data_filtered$Visit,
                                     AST_data_filtered$Sample,
                                     AST_data_filtered$Pick,
                                     sep = "_")


#Retrieving the list of antibiotics, toxins, fimbriae, visits and groups (pre-/post-withdrawal)
antibiotics <- colnames(AST_data_filtered)[7:20]
toxins <- colnames(AST_data_filtered)[21:22] #Only STa and LTa as no Stx2e positive isolates
fimbriae <- colnames(AST_data_filtered)[24:25] #Only F4 and F18 as other fimbriae not detected by PCR
visits <- unique(AST_data_filtered$Visit)
groups <- unique(AST_data_filtered$Group)


#Creating lists for each class of antibiotics studied
classes <- list(
  Tetracyclines = c("Tetracycline"),
  Penicillins = c("Ampicillin", "Amoxicillin_Clavulanate"),
  Aminoglycosides = c("Apramycin", "Gentamicin", "Streptomycin", "Spectinomycin"),
  Sulphonamides = c("Sulfamethoxazole_Trimethoprim"),
  Amphenicols = c("Florfenicol", "Chloramphenicol"),
  Cephalosporins3G = c("Ceftiofur", "Cefotaxime"),
  Quinolones = c("Enrofloxacin", "Ciprofloxacin")
)


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


#Creating a table indicating to how many antibiotics an isolate is resistant
MDR_data <- RS_data[, c(1:8)]

MDR_data$Nb_ATB_R <- apply(RS_data[, antibiotics], 1,
                              function(x) sum(x == "R"))


#Transforming the data et to have the following columns = visit (V1, V2, V4) and agar type (i.e. separate picks from unselective or selective)
#and the number of antibiotics each isolate is resistant to = Nb_ATB_R
#And for each category indicated by the columns above, indicate number of R bacteria and %, total number of bacteria
MDR_data_groupedNbR <- MDR_data %>%
  group_by(Visit, Nb_ATB_R, Agar) %>%
  summarise(Number_Isolates = n(), .groups = "drop") %>% # Counting the number of isolates for each visit, agar type and level of resistance
  group_by(Visit, Agar) %>%  # Grouping by visit and agar type (goal = that per agar type and per visit, sum proportions = 100%)
  mutate(
    Total_Isolates_visit_agar = sum(Number_Isolates),  # Number of isolates by visit and agar type
    Proportion = Number_Isolates / Total_Isolates_visit_agar,
    Percentage = Proportion * 100,
    # Calculating to display the confidence interval 95%
    Standard_error = sqrt((Proportion * (1 - Proportion)) / Total_Isolates_visit_agar),
    CI95 = 1.96 * Standard_error,
    ymin = (Proportion - CI95) * 100, #Lower limit CI95
    ymax = (Proportion + CI95) * 100, #Upper limit CI95
    ymin = pmax(ymin, 0),
    ymax = pmin(ymax, 100)
  ) %>%
  ungroup()
  

#Plotting the data as barplots (x axis = level of resistance, 2 graphs = selective agar vs unselective agar)
ymax <- ceiling(max(max(MDR_data_groupedNbR$Percentage), max(MDR_data_groupedNbR$ymax)) / 5) * 5
xmin = min(MDR_data_groupedNbR$Nb_ATB_R)
xmax = max(MDR_data_groupedNbR$Nb_ATB_R)

NbR_ATB_barplots <- ggplot(MDR_data_groupedNbR,
                      aes(x = Nb_ATB_R, fill = Visit)) +
  
  #Plotting picks from selective and unselective agar separately
  facet_wrap(~ Agar, ncol = 1,
             labeller = as_labeller(c("Selective" = "Selective agar",
                                      "Unselective" = "Unselective agar"))) +
  
  geom_bar(aes(y = Percentage),
           width = 0.5,
           stat = "identity", position = position_dodge(width = 0.45)) +
  
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),  # CI95
    position = position_dodge(width = 0.45), 
    width = 0.1, 
    color = "black") +
  
  theme_bw() +
  
  scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
  
  scale_y_continuous(limits = c(0, ymax),
                     breaks = seq(0, ymax, by = 5)) + # forcing R to display from 0% to 100% on the y axis
  
  scale_x_continuous(breaks = seq(xmin, xmax, by = 1)) +
  
  labs(title = "Distribution of Isolates by Number of Antibiotic Resistances",
       x = "Number of antibiotic resistances",
       y = "Percentage of resistant isolates (%)") +
  
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        
        #Changing font size
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 22)) +
  
  #Changing legend overall size
  guides(fill = guide_legend(keywidth = 2, keyheight = 2.5))      


#Saving the barplots in the "Graphs" folder
ggsave("Graphs/NbR_ATB_barplots.png", plot = NbR_ATB_barplots,
       width = 15, height = 15, dpi = 300)





#Adding a boxplot to see whether differences in terms of median level of resistance between V1, V2 and V4 isolates
median_ATB_R_visits <- MDR_data %>%
  group_by(Visit, SampleID, Agar) %>%
  summarise(Median_ATB_R = median(Nb_ATB_R),
            .groups = "drop") #in order to remove group structure

#Statistical tests to determine whether there is a significant difference between groups
#Kruskal-Wallis to determine whether significant difference between groups
kruskal.test(Median_ATB_R ~ Visit, data = median_ATB_R_visits, 
             subset = Agar == "Unselective")

kruskal.test(Median_ATB_R ~ Visit, data = median_ATB_R_visits, 
             subset = Agar == "Selective")

#Dunn test to calculate pairwise comparisons between groups (post-hoc test after Kruskal-Wallis)
dunnTest(Median_ATB_R ~ Visit,
         data = median_ATB_R_visits[median_ATB_R_visits$Agar == "Unselective", ],
         method = "bh") #Benjamini-Hochberg method

dunnTest(Median_ATB_R ~ Visit,
         data = median_ATB_R_visits[median_ATB_R_visits$Agar == "Selective", ],
         method = "bh") #Benjamini-Hochberg method


#Computing a few statistic parameters for each visit and agar type
median_ATB_R_visits_stats <- median_ATB_R_visits %>%
  group_by(Visit, Agar) %>%
  summarise(
    Median = median(Median_ATB_R),
    Min = min(Median_ATB_R),
    Max = max(Median_ATB_R),
    Q1 = quantile(Median_ATB_R, probs = c(0.25)),
    Q3 = quantile(Median_ATB_R, probs = c(0.75)),
    IQR = Q3 - Q1,
    .groups = "drop")


#Adding the boxplot to display the data
Median_ATB_R_boxplot <- ggplot(median_ATB_R_visits,
                          aes(x = Visit, y = Median_ATB_R, fill = Visit)) +
  
  geom_boxplot() +
  
  theme_bw() +
  
  scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
  
  facet_wrap(~ Agar, ncol = 1,
             labeller = as_labeller(c("Selective" = "Selective agar",
                                      "Unselective" = "Unselective agar"))) +
  
  scale_y_continuous(limits = c(0, 14),
                     breaks = seq(0, 14, by = 1)) + # forcing R to display y axis from 0 to 14 = max nb of resistances
  
  stat_summary(fun = mean, geom = "point",
               shape = 4, size = 2, color = "black", fill = "black") + #Showing the mean as a black cross
  
  labs(title = "Median level of resistance",
       x = NULL,
       y = "Median level of resistance") +
  
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text.x = element_blank(), # Removing x labels as already indicated by fill = Group
        axis.ticks.x = element_blank(), # removing ticks on x axis
        #Changing font size
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 22)) +
  
  #Changing legend overall size
  guides(fill = guide_legend(keywidth = 2, keyheight = 2.5))


#Saving the boxplots in the "Graphs" folder
ggsave("Graphs/Median_ATB_R_boxplot.png", plot = Median_ATB_R_boxplot,
       width = 10, height = 15, dpi = 300)










#Adding to MDR_data an column indicating whether an isolate is resistant to at least one antibiotic in each class
for (class in names(classes)) {
  
  antibiotics <- classes[[class]]
  
  resistance_class <- character(nrow(RS_data))
  
  for (i in seq_along(RS_data$IsolateID)) {
    isolat <- RS_data$IsolateID[i]
    isolate_R_class <- RS_data[RS_data$IsolateID == isolat, antibiotics, drop = FALSE]
    
    resistance_class[i] <- if (any(isolate_R_class == "R", na.rm = TRUE)) "Yes" else "No"
  }
  
  MDR_data[[paste0("R_", class)]] <- resistance_class
}


#And adding another column indicating whether an isolate is MDR or not
#In that case, isolate classified as MDR if resistant to at least one antibiotic in 3 or more classes
MDR_data$MDR <- apply(MDR_data[, 9:15], 1, function(x) {
  if (sum(x == "Yes", na.rm = TRUE) >= 3) "Yes" else "No"
})


#Transforming the dataset to have the following columns = sample ID, visit (V1, V2, V4) and agar type (selective or unselective)
#And for each category indicated by the columns above, indicate number of R bacteria and %, total number of bacteria
MDR_data_groupedMDR <- MDR_data %>%
  group_by(SampleID, Visit, Agar) %>%
  summarise(Number_MDR = sum(MDR == "Yes"),
            Number_Total = n(),
            Percentage_MDR = 100 * Number_MDR / Number_Total,
            .groups = "drop")


#Computing a few statistic parameters for each antibiotic R rates (median and mean) and separating picks from selective and unselective agar
MDR_data_groupedMDR_stats <- MDR_data_groupedMDR %>%
  group_by(Visit, Agar) %>%
  summarise(Mean = mean(Percentage_MDR),
            Median = median(Percentage_MDR),
            Min = min(Percentage_MDR),
            Q1 = quantile(Percentage_MDR, probs = c(0.25)),
            Q3 = quantile(Percentage_MDR, probs = c(0.75)),
            Max = max(Percentage_MDR),
            .groups = "drop")


#Statistical tests to determine whether there is a significant difference between groups
#Kruskal-Wallis to determine whether significant difference between groups
kruskal.test(Percentage_MDR ~ Visit, data = MDR_data_groupedMDR, 
             subset = Agar == "Unselective")

kruskal.test(Percentage_MDR ~ Visit, data = MDR_data_groupedMDR, 
             subset = Agar == "Selective")

#Dunn test to calculate pairwise comparisons between groups (post-hoc test after Kruskal-Wallis)
dunnTest(Percentage_MDR ~ Visit,
         data = MDR_data_groupedMDR[MDR_data_groupedMDR$Agar == "Unselective", ],
         method = "bh") #Benjamini-Hochberg method

dunnTest(Percentage_MDR ~ Visit,
         data = MDR_data_groupedMDR[MDR_data_groupedMDR$Agar == "Selective", ],
         method = "bh") #Benjamini-Hochberg method


#Plotting the data as boxplots with one boxplot for picks from selective and one for picks from unselective
MDR_boxplots <- ggplot(MDR_data_groupedMDR,
                       aes(x = Visit, y = Percentage_MDR, fill = Visit)) +
  
  geom_boxplot() +
  
  theme_bw() +
  
  # Plotting picks from selective and unselective agar side by side
  facet_wrap(~ Agar, ncol = 2,
             labeller = as_labeller(c("Selective" = "Selective agar",
                                      "Unselective" = "Unselective agar"))) +
  
  scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
  
  scale_y_continuous(limits = c(0, 110),
                     breaks = seq(0, 100, by = 10)) + # forcing R to display from 0% to 100% on the y axis
  
  stat_summary(fun = mean, geom = "point",
               shape = 4, size = 2, color = "black", fill = "black") + #Showing the mean as a black cross
  
  labs(title = "Percentage of MDR isolates",
       x = NULL,
       y = "Percentage of MDR isolates (%)") +
  
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text.x = element_blank(), # Removing x labels as already indicated by fill = Group
        axis.ticks.x = element_blank(), # removing ticks on x axis
        #Changing font size
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 18)) +
  
  #Changing legend overall size
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 2))


#Saving the boxplots in the "Graphs" folder
ggsave("Graphs/MDR_boxplots.png", plot = MDR_boxplots,
       width = 12, height = 7, dpi = 300)










# Frequent patterns of phenotypic resistance can be identified (i.e. associations/combinations of antibiotics the isolates are resistant to)
# We can examine the overall proportion of resistant isolates for each pattern
RS_data_codeATB <- RS_data

#Renaming the columns with 3 letter codes
colnames(RS_data_codeATB)[9:22] <- c("TET", "AMP", "AMC", "APR", "GEN", "STR", "SPE",
                                     "SUT", "FLO", "CHL", "CFT", "CTX", "ENR", "CIP")


#Generating all the antibiotics combinations possible of 1 to 11 antibiotics (nb max of resistances observed in our dataset)
combinations_ATB <- lapply(1:11, function(k) {
  combn(colnames(RS_data_codeATB)[9:22], k, simplify = FALSE)
})

all_combinations_ATB <- unlist(combinations_ATB, recursive = FALSE)

#Giving names to the associations (i.e. asso1, asso2, etc)
names(all_combinations_ATB) <- paste0("asso", seq_along(all_combinations_ATB))

#Creating a list of the combinations
ATBasso_list <- all_combinations_ATB

#Sorting the association from longest to shortest (because we want to find longest association an isolate is resistant to)
ATBasso_list_sorted <- ATBasso_list[order(-sapply(ATBasso_list, length))]


#Creating a list to indicate association to which the isolate is resistant (and avoid assigning an isolate twice)
assigned_asso <- rep(NA, nrow(RS_data_codeATB))

#Loop for each association
for (asso_name in names(ATBasso_list_sorted)) {
  #Retrievng list of ATB of that association
  atb_group <- ATBasso_list_sorted[[asso_name]]
  #Checking antibiotics name are in the raw dataframe (e.g. RS_data_codeATB)
  atb_group <- atb_group[atb_group %in% colnames(RS_data_codeATB)]
  
  if (length(atb_group) > 0) {
    #Identifying the isolates resistant to all the ATB of that association
    resistant_all <- apply(RS_data_codeATB[, atb_group, drop = FALSE], 1, function(row) all(row == "R"))
   
    #If isolate not yet assigned to an association (is.na(assigned_asso))
    #and identified as resistant to the association considered (resistant_all = TRUE)
    newly_assigned <- which(is.na(assigned_asso) & resistant_all)
    #If isolate resistant to association, adding association in the list
    assigned_asso[newly_assigned] <- asso_name
  }
}


#Adding the patterns to the whole dataframe (indicating Unassigned if no association = isolate susceptible to all the antibiotics tested)
RS_data_codeATB$Name_ATB_asso <- ifelse(is.na(assigned_asso), "Unassigned",
                                        assigned_asso)

ATB_asso_counts <- RS_data_codeATB %>%
  group_by(Name_ATB_asso, Visit, Agar) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Visit, values_from = n, values_fill = 0) %>%
  
  mutate(
    #Total = total nb of picks with that asso (from V1 + V2 + V4)
    Total = V1 + V2 + V4,
    
    #Adding list of the ATB in that association
    ATB_asso = sapply(as.character(Name_ATB_asso), function(name) {
      if (name %in% names(ATBasso_list)) {
        paste(ATBasso_list[[name]], collapse = ", ")
      } else {
        NA  #If unassigned
      }
    })
  ) %>%
  
  group_by(Agar) %>%
  mutate(Prop_V1 = 100 * V1 / sum(V1),
         Prop_V2 = 100 * V2 / sum(V2),
         Prop_V4 = 100 *V4 / sum(V4)) %>%   #Calculating the proportion by agar type for each visit
  ungroup() %>%
  
  filter(Total > 0) %>% #Keeping only associations for which there is at least one isolate
  arrange(desc(Total)) %>%
  
  #Forcing R to display columns in a specific order
  select(Name_ATB_asso, Agar, V1, Prop_V1, V2, Prop_V2, V4, Prop_V4, Total, ATB_asso)


#Verifying that all the isolates have been assigned
sum(ATB_asso_counts$Total) == nrow(RS_data_codeATB)

#Saving results as an excel file
write_xlsx(ATB_asso_counts, "ATB_asso_counts.xlsx")










#In order to look at patterns of MDR (i.e. for an MDR isolate, resistance to which classes is the most frequently found?)
#First of all, displaying data of MDR data on a heatmap
#Converting  data as a matrix (required for pheatmap)
MDR_data_4heatmap <- as.matrix(MDR_data[, 10:16])

#Converting Yes as 1 (resistant to at least one antibiotic of the class) and No as 0 (susceptible to all antibiotics) and all values as numeric
MDR_data_4heatmap <- ifelse(MDR_data_4heatmap == "Yes", 1,
                            ifelse(MDR_data_4heatmap == "No", 0, NA))

MDR_data_4heatmap <- matrix(as.numeric(MDR_data_4heatmap),
                            nrow = nrow(MDR_data_4heatmap))


#Renaming the rows of with IsolateID (without Farm ID for confidentiality reasons)
rownames(MDR_data_4heatmap) <- AST_data_filtered$IsolateID


#Adding another annotation based on whether the isolates are from pre- (V1) or post- (V2 and V4) sampling
annotation_visits_row <- MDR_data %>%
  {rownames(.) <- NULL; .} %>%
  select(IsolateID, Visit) %>%
  mutate(Visit = as.factor(Visit)) %>%
  column_to_rownames(var = "IsolateID")


#Adding a way to annotate graph to indicate whether an isolate is MDR
annotation_MDR_row <- data.frame(MDR_data$MDR)

#Renaming the rows of annotation_MDR_row with IsolateID (so that they match with data_antibiotics_normalized)
rownames(annotation_MDR_row) <- MDR_data$IsolateID

#And changing name of MDR column
colnames(annotation_MDR_row) <- c("MDR")

#Adding annotation of toxins and fimbriae
data_toxins <- AST_data_filtered[, toxins]
annotation_toxins_row <- data.frame(data_toxins)
rownames(annotation_toxins_row) <- AST_data_filtered$IsolateID

data_fimbriae <- AST_data_filtered[, fimbriae]
annotation_fimbriae_row <- data.frame(data_fimbriae)
rownames(annotation_fimbriae_row) <- AST_data_filtered$IsolateID

#And changing NA as text for fimbriae (analysis not done) so that it can be understood by heatmap
annotation_fimbriae_row[is.na(annotation_fimbriae_row)] <- "NA"


#Defining colours for the resistance to each class (red if resistant, green if susceptible)
#Note, two values = 0 and 1, the way the palette below is defined, Rstudio colours 0 in green and 1 in red (resistant to at least one antibiotic of the class)
palette_classes_colours <- c("#006600", "#990000")


#Defining colours for the visit time points
annotation_visits_colours <- list(
  Visit = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")
)

#Defining colours for to indicate whether an isolate is MDR (black) or not (white)
annotation_MDR_colours <- list(
  MDR = c("Yes" = "black", "No" = "white")
)


#Defining colours for the toxins
annotation_toxins_colours <- list()
for (toxin in toxins) {
  if (toxin == "Stx2e") {
    annotation_toxins_colours[[toxin]] <- c("Yes" = "lightsteelblue", "No" = "white")
  } else if (toxin == "LTa") {
    annotation_toxins_colours[[toxin]] <- c("Yes" = "powderblue", "No" = "white")
  } else if (toxin == "STa") {
    annotation_toxins_colours[[toxin]] <- c("Yes" = "cadetblue", "No" = "white")
  }
}


#Defining colours for the fimbriae
annotation_fimbriae_colours <- list()
for (fim in fimbriae) {
  if (fim == "F4") {
    annotation_fimbriae_colours[[fim]] <- c("Yes" = "peachpuff", "No" = "white", "NA" = "#E0E0E0")
  } else if (fim == "F18") {
    annotation_fimbriae_colours[[fim]] <- c("Yes" = "sandybrown", "No" = "white", "NA" = "#E0E0E0")
  } 
}


#Combining the annotations without (combined), or with the fimbriae (all)
annotation_row_combined <- cbind(annotation_MDR_row, annotation_toxins_row, annotation_visits_row)
annotation_row_all <- cbind(annotation_MDR_row, annotation_fimbriae_row, annotation_toxins_row, annotation_visits_row)


#Combining the colours without (combined), or with the fimbriae (all)
annotation_colours_combined <- c(annotation_MDR_colours, annotation_toxins_colours, annotation_visits_colours)
annotation_colours_all <- c(annotation_MDR_colours, annotation_fimbriae_colours, annotation_toxins_colours, annotation_visits_colours)


#Renaming column names of the dataframe
colnames(MDR_data_4heatmap) <- c("Tetracyclines", "Penicillins", "Aminoglycosides", "Sulphonamides",
                                 "Amphenicols", "Cephalosporins3G", "Quinolones")


#Creating the heatmap with toxins and visit and saving it in the graphs folder
filename <- "Graphs/Classes_heatmap_ToxVisit.png"
print(paste("Saving:", filename))
png(filename, width = 3000, height = 5500, res = 400)

pheatmap(MDR_data_4heatmap,
         cluster_rows = TRUE, #clustering the isolates
         cluster_cols = TRUE, #clustering the classes
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = "ward.D2",
        
         annotation_row = annotation_row_combined,
         annotation_colors = annotation_colours_combined,
         color = palette_classes_colours,
         scale = "none",
         fontsize_row = 4,
         fontsize_col = 10,
         fontsize = 10,
         #main = "Isolates - MDR Profiles",
         legend_breaks = c(0, 1),  #Antibiotics legend breaks
         legend_labels = c("Susceptible", "Resistant"),
         angle_col = 45) #rotating antibiotics and toxins names

dev.off()


#Creating the heatmap with toxins, fimbriae and visit and saving it in the graphs folder
filename <- "Graphs/Classes_heatmap_ToxFimVisit.png"
print(paste("Saving:", filename))
png(filename, width = 3000, height = 5500, res = 400)

pheatmap(MDR_data_4heatmap,
         cluster_rows = TRUE, #clustering the isolates
         cluster_cols = TRUE, #clustering the classes
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = "ward.D2",
         
         annotation_row = annotation_row_all,
         annotation_colors = annotation_colours_all,
         color = palette_classes_colours,
         scale = "none",
         fontsize_row = 4,
         fontsize_col = 10,
         fontsize = 10,
         #main = "Isolates - MDR Profiles",
         legend_breaks = c(0, 1),  #Antibiotics legend breaks
         legend_labels = c("Susceptible", "Resistant"),
         angle_col = 45) #rotating antibiotics and toxins names

dev.off()


#Generating all the classes combinations/associations possible of 3 (nb min to define MDR) to 6 (nb max observed in our dataset) classes
combinations_classes <- lapply(3:6, function(k) {
  combn(names(classes), k, simplify = FALSE)
})

all_combinations_classes <- unlist(combinations_classes, recursive = FALSE)

#Giving names to the associations (i.e. asso1, asso2, etc)
names(all_combinations_classes) <- paste0("asso", seq_along(all_combinations_classes))

#Creating a list of the combinations
Classesasso_list <- all_combinations_classes

#Sorting the association from longest to shortest (because we want to find longest association an isolate is resistant to)
Classesasso_list_sorted <- Classesasso_list[order(-sapply(Classesasso_list, length))]


#Renaming column names of MDR_data
colnames(MDR_data)[10:16] <- c("Tetracyclines", "Penicillins", "Aminoglycosides", 
                               "Sulphonamides", "Amphenicols", "Cephalosporins3G", 
                               "Quinolones")


#Creating a list to indicate association to which the isolate is resistant (and avoid assigning an isolate twice)
assigned_classes_asso <- rep(NA, nrow(MDR_data))

#Loop for each association
for (asso_name in names(Classesasso_list_sorted)) {
  #Retrievng list of classes of that association
  class_group <- Classesasso_list_sorted[[asso_name]]
  #Checking antibiotics name are in the raw dataframe (e.g. RS_data_codeATB)
  class_group <- class_group[class_group %in% colnames(MDR_data)]
  
  if (length(class_group) > 0) {
    #Identifying the isolates resistant to all the classes of that association
    resistant_all <- apply(MDR_data[, class_group, drop = FALSE], 1, function(row) all(row == "Yes"))
    
    #If isolate not yet assigned to an association (is.na(assigned_asso))
    #and identified as resistant to the association considered (resistant_all = TRUE)
    newly_assigned <- which(is.na(assigned_classes_asso) & resistant_all)
    #If isolate resistant to association, adding association in the list
    assigned_classes_asso[newly_assigned] <- asso_name
  }
}


#Adding the patterns to the whole dataframe (indicating Unassigned if no association)
MDR_data$Name_Classe_asso <- ifelse(is.na(assigned_classes_asso), "Unassigned",
                                    assigned_classes_asso)

Classes_asso_counts <- MDR_data %>%
  group_by(Name_Classe_asso, Visit, Agar) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Visit, values_from = n, values_fill = 0) %>%
  
  mutate(
    #Total = total nb of picks with that asso (from V1 + V2 + V4)
    Total = V1 + V2 + V4,
    
    #Adding list of the ATB in that association
    Classe_asso = sapply(as.character(Name_Classe_asso), function(name) {
      if (name %in% names(Classesasso_list)) {
        paste(Classesasso_list[[name]], collapse = ", ")
      } else {
        NA  #If unassigned
      }
    })
  ) %>%
  
  group_by(Agar) %>%
  mutate(Prop_V1 = 100 * V1 / sum(V1),
         Prop_V2 = 100 * V2 / sum(V2),
         Prop_V4 = 100 *V4 / sum(V4)) %>%   #Calculating the proportion by agar type for each visit
  ungroup() %>%
  
  filter(Total > 0) %>% #Keeping only associations for which at least one isolate has it
  arrange(desc(Total)) %>%
  
  #Forcing R to display columns in a specific order
  select(Name_Classe_asso, Agar, V1, Prop_V1, V2, Prop_V2, V4, Prop_V4, Total, Classe_asso)


#Verifying that all the isolates have been assigned
sum(Classes_asso_counts$Total) == nrow(MDR_data)


#Saving as an excel file
write_xlsx(Classes_asso_counts, "Classes_asso_counts.xlsx")


#Closing all plots (just in case)
if (dev.cur() != 1) {
  dev.off()
}