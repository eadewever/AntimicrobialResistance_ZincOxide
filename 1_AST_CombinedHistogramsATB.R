#This script aims at displaying the distribution of the isolates based on their Antimicrobial Susceptibility values (i.e. zone diameter in millimeters)

# ctrl + alt + R to run the whole script

#Set working directory
setwd("indicate working directory")


#Install necessary packages
packages <- c("ggplot2", "tidyr", "dplyr", "tibble")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

#Loading packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

#Creating a folder to save the graphs
if (!dir.exists("Graphs")) {
  dir.create("Graphs")}





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


#Retrieving the list of antibiotics, toxins, visits and samples
antibiotics <- colnames(AST_data_filtered)[7:20]
toxins <- colnames(AST_data_filtered)[21:23]
visits <- unique(AST_data_filtered$Visit)
samples <- unique(AST_data_filtered$Sample)


#Creating a table with, for each antibiotic, the threshold zone diameters defining resistant and susceptible isolates
#Epidemiological breakpoint values according to EUCAST 2024 and VARSS 2024 (values listed in the order of the antibiotics aforementioned)
#Diameters defining resistant isolates (resistant if <= list_R_diameter_adapted value)
list_R_diameter <- c(11,13,18,11,16,10,17,10,18,16,17,16,16,21)
#Diameters defining susceptible isolates (susceptible if >= list_S_diameter_adapted value)
list_S_diameter <- c(15,14,19,15,17,15,21,14,19,17,21,20,23,25)
#If list_R_diameter != list_S_diameter - 1 = an "intermediary" category exists, i.e. isolate considered as susceptible, increased exposure
table_RS_ATB <- data.frame(Antibiotic = antibiotics, R = list_R_diameter, S = list_S_diameter)


#Transforming the data set so that each antibiotic is displayed on rows and not columns
#Here, we keep the information of the toxins so that we can refer to AST_data_filtered to look at the toxins expressed
AST_data_filtered <- AST_data_filtered %>%
  mutate(across(c(all_of(antibiotics), all_of(toxins)), as.character))

AST_data_long <- pivot_longer(
  AST_data_filtered,
  cols = c(all_of(antibiotics), all_of(toxins)),
  names_to = "Variable",         # Antibiotic or toxin
  values_to = "Value"            # Zone diameter for antibiotic, Yes/No for toxin
) %>%
  mutate(Type = case_when(
    Variable %in% antibiotics ~ "Antibiotic",
    Variable %in% toxins ~ "Toxin",
    TRUE ~ "Other"
  ))

AST_data_long <- AST_data_long %>%
  mutate(
    Zone_diameter = case_when(
      Type == "Antibiotic" ~ suppressWarnings(as.numeric(Value)), #suppressWarnings because Rstudio gives warning as Value column contains numerical values for ATB and Yes/No for toxins
      TRUE ~ NA_real_ 
    )
  )


#Forcing R to display the picks in the order of the file and not alphabetical order
AST_data_long$Pick <- factor(AST_data_long$Pick, 
                             levels = c(paste0("P", 1:10), "AmpR_P1", "AmpR_P2", "AprR_P1", "AprR_P2"))





#Displaying the distributions of the picks = f(zone diameter in mm)

#Slightly changing table_RS_ATB as breakpoint values for a nicer display on the graphs
#Diameters defining resistant isolates 
list_R_diameter_adapted <- c(11.5, 13.5, 18.5, 11.5, 16.5, 10.5, 17.5,
                             10.5, 18.5, 16.5, 17.5, 16.5, 16.5, 21.5)
#Diameters defining susceptible isolates
list_S_diameter_adapted <- c(14.5, 13.5, 18.5, 14.5, 16.5, 14.5, 20.5,
                             13.5, 18.5, 16.5, 20.5, 19.5, 22.5, 24.5)
#Combining the lists
table_RS_ATB_adapted <- data.frame(Antibiotic = antibiotics, R = list_R_diameter_adapted, S = list_S_diameter_adapted)


#For loop for each visit and each antibiotic
for (antibiotic in antibiotics) {

  #Filtering the data to get the desired visit and antibiotic
  #And calculating the proportion of picks rather than the absolute number of picks
  proportion_data <- AST_data_long %>%
    filter(Type == "Antibiotic", Variable == antibiotic) %>%
    count(Visit, Zone_diameter, name = "Nb_picks") %>%  # Counting the number of picks per visit and per zone diameter
    group_by(Visit) %>%
    mutate(
      Nb_picks_visit = sum(Nb_picks),  # Total number of picks per visit
      Proportion = Nb_picks / Nb_picks_visit,  # For each visit, proportion per zone diameter
      Percentage = Proportion * 100, # Corresponding percentage
      # Calculating the confidence interval at 95% (can be displayed on the graphs)
      Standard_error = sqrt((Proportion * (1 - Proportion)) / Nb_picks_visit),
      CI95 = 1.96 * Standard_error,
      ymin = (Proportion - CI95) * 100,
      ymax = (Proportion + CI95) * 100,
      ymin = pmax(ymin, 0),
      ymax = pmin(ymax, 100)
    ) %>%
    ungroup()
  
  #Max value of the x axis = max of Zone_Diameter or S value of table_RS_ATB_adapted (if all picks susceptible)
  #To avoid cases when all picks are resistant, i.e. sensitive area not displayed
  max_x_limit <- max(max(proportion_data$Zone_diameter, na.rm = TRUE),
                     table_RS_ATB_adapted[table_RS_ATB_adapted$Antibiotic == antibiotic, "S"])
  
  #Creating the graph, displaying V1, V2 and V4 isolates (fill = Visit)
  p <- ggplot(proportion_data, aes(x = Zone_diameter, fill = Visit)) +
    
    geom_bar(aes(y = Percentage), stat = "identity", position = position_dodge(width = 0.9)) +
    
    #Indicating colour for each visit and displaying the legend
    scale_fill_manual(values = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")) +
    labs(fill = "Visit") +

    theme_bw() +

    labs(title = antibiotic,
         x = "Zone diameter (mm)",
         y = "Percentage of picks (%)") +

    #Limits x and y axis
    scale_x_continuous(limits = c(-0.5, max_x_limit +1),
                       breaks = seq(0, max_x_limit + 1, by = 2),
                       expand = c(0, 1)) +

    scale_y_continuous(breaks = seq(0, 100, by = 5),
                       expand = c(0, 3)) +

    #Rectangles to indicate the Susceptible and Resistant zones/diameters
    geom_rect(aes(xmin = table_RS_ATB_adapted[table_RS_ATB_adapted$Antibiotic == antibiotic, "S"], xmax = Inf,
                  ymin = -Inf, ymax = Inf),
              fill = NA,
              color = "darkgreen",
              linetype = "solid",
              linewidth = 1.1) +
      
    geom_rect(aes(xmin = -Inf, xmax = table_RS_ATB_adapted[table_RS_ATB_adapted$Antibiotic == antibiotic, "R"],
                  ymin = -Inf, ymax = Inf),
              fill = NA,
              color = "darkred",
              linetype = "dashed",
              linewidth = 1.1) +
      
    #Annotating the Resistant and Susceptible zones
    annotate("text", x = table_RS_ATB_adapted[table_RS_ATB_adapted$Antibiotic == antibiotic, "R"] - 3,
              y = -3,
              label = "Resistant", color = "darkred",
              size = 5, fontface = "bold") +
      
    annotate("text", x = table_RS_ATB_adapted[table_RS_ATB_adapted$Antibiotic == antibiotic, "S"] + 3,
              y = -3,
              label = "Susceptible", color = "darkgreen",
              size = 5, fontface = "bold") +
      
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          panel.grid.minor.y = element_blank(), #Removing minor grid lines y axis
          #Changing font size
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16))

  #Creating filename and printing it
  filename <- paste0("Graphs/", antibiotic, "_histogram.png")
  print(paste("Saving:", filename))

  #Saving the file
  ggsave(filename,
          plot = p,
          width = 10, height = 6, dpi = 200)
}