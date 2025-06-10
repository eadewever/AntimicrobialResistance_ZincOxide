#This script aims at displaying the phenotypic resistance profiles of additional isolates using heatmaps (i.e. for additional isolates not considered in the main analysis)

# ctrl + alt + R to run the whole script

#Set working directory
setwd("indicate working directory")


#Install necessary packages, doing it if and only if not already installed
packages <- c("ggplot2", "tidyr", "dplyr", "pheatmap", "tibble")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

#Loading packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(tibble)

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


#Retrieving the list of antibiotics, toxins, fimbriae and visits
antibiotics <- colnames(AST_data_filtered)[7:20]
toxins <- colnames(AST_data_filtered)[21:22] #Only STa and LTa as no Stx2e positive isolates
fimbriae <- colnames(AST_data_filtered)[24:25] #Only F4 and F18 as other fimbriae not detected
visits <- unique(AST_data_filtered$Visit)


#Adding a column to AST_data_filtered with one unique ID per isolate (removing Farm ID for confidentiality reasons)
AST_data_filtered$IsolateID <- paste(AST_data_filtered$Visit,
                                     AST_data_filtered$Sample,
                                     AST_data_filtered$Pick,
                                     sep = "_")


#Creating a table with, for each antibiotic, the threshold zone diameters defining resistant and susceptible isolates
#Epidemiological breakpoint values according to EUCAST 2024 and VARSS 2024 (values listed in the order of the antibiotics aforementioned)
#Diameters defining resistant isolates (resistant if <= list_R_diameter_adapted value)
list_R_diameter <- c(11,13,18,11,16,10,17,10,18,16,17,16,16,21)
#Diameters defining susceptible isolates (susceptible if >= list_S_diameter_adapted value)
list_S_diameter <- c(15,14,19,15,17,15,21,14,19,17,21,20,23,25)
#If list_R_diameter != list_S_diameter - 1 = an "intermediary" category exists, i.e. isolate considered as susceptible, increased exposure
table_RS_ATB <- data.frame(Antibiotic = antibiotics, R = list_R_diameter, S = list_S_diameter)


#Creating separate tables, for the antibiotics, toxins and fimbriae
data_antibiotics <- AST_data_filtered[, antibiotics]
data_toxins <- AST_data_filtered[, toxins]
data_fimbriae <- AST_data_filtered[, fimbriae]


#Convert antibiotics data as numeric values
data_antibiotics <- data.frame(lapply(data_antibiotics,
                                      function(x) as.numeric(as.character(x))))


#Normalisation of the zone diameter values
#We divide by resistant breakpoint values
#i.e. if normalisation > 1 = pick susceptible (S or I), vs if <= 1 = pick resistant
ref_values_R <- data.frame(
  Antibiotic = antibiotics,
  Reference_Zone_R = list_R_diameter)


#The order in which the antibiotics are listed should be the same in the two dataframe but we realign them in case
Reference_Zone_R <- ref_values_R$Reference_Zone_R[match(colnames(data_antibiotics),
                                                        ref_values_R$Antibiotic)]


#We normalise the values for the zone diameter
#normalise based on columns = 2 ; and by dividing by reference values = FUN = "/"
data_antibiotics_normalised <- sweep(data_antibiotics, 2, Reference_Zone_R, FUN = "/")


#Converting antibiotics data as a matrix (required for pheatmap)
data_antibiotics_normalised <- as.matrix(data_antibiotics_normalised)


#Renaming column names of the dataframe to give 3 letter code for each antibiotic
colnames(data_antibiotics_normalised) <- c("TET", "AMP", "AMC", "APR", "GEN", "STR", "SPE",
                                           "SUT", "FLO", "CHL", "CFT", "CTX", "ENR", "CIP")


#Renaming the rows of data_antibiotics_normalised with IsolateID
rownames(data_antibiotics_normalised) <- AST_data_filtered$IsolateID


#Adding a way to annotate graph with results for the toxins and fimbriae analysis
annotation_toxins_row <- data.frame(data_toxins)
annotation_fimbriae_row <- data.frame(data_fimbriae)

#And changing NA as text for fimbriae (i.e. PCR not done) so that it can be understood by heatmap
annotation_fimbriae_row[is.na(annotation_fimbriae_row)] <- "NA"


#Renaming the rows with IsolateID (so that they match with data_antibiotics_normalised)
rownames(annotation_toxins_row) <- AST_data_filtered$IsolateID
rownames(annotation_fimbriae_row) <- AST_data_filtered$IsolateID


#Adding another annotation based on whether the isolates are from pre- (V1) or post- (V2 and V4) sampling
annotation_visits_row <- AST_data_filtered %>%
  select(Visit, IsolateID) %>%
  mutate(Visit = as.factor(Visit)) %>%
  column_to_rownames(var = "IsolateID") %>%
  select(Visit)


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


#Defining colours for the visit time points
annotation_visits_colours <- list(
  Visit = c("V1" = "lightcoral", "V2" = "darkolivegreen3", "V4" = "mediumturquoise")
)


#Combining the annotations for the toxins and visit (combined), or toxins, fimbriae and visit (all)
annotation_row_combined <- cbind(annotation_toxins_row, annotation_visits_row)
annotation_row_all <- cbind(annotation_fimbriae_row, annotation_toxins_row, annotation_visits_row)


#Combining the colours for the toxins and visits annotations (combined), or toxins, fimbriae and visit (all)
annotation_colors_combined <- c(annotation_toxins_colours, annotation_visits_colours)
annotation_colors_all <- c(annotation_fimbriae_colours, annotation_toxins_colours, annotation_visits_colours)






#Creating one heatmap per sample to see the phenotypic resistance diversity within each sample
#Extracting the group_ID for each group abovementioned to create one heatmap per group
#group_ID = Sample_ID if additional isolates isolated from one sample
#group_ID != Sample_ID if one isolate subcultured and new isolates obtained from that single isolate (should all be the same bacteria = same phenotypic resistance profile)
groups_ID <- sub("(_[^_]+)$", "", rownames(data_antibiotics_normalised))

#Creating lists of all the isolates ID per sample
samples_isolateslist <- split(rownames(data_antibiotics_normalised), groups_ID)

#Separating the whole matrix into sub-matrices with one by groups
samples_submatrices <- lapply(samples_isolateslist,
                              function(rows) data_antibiotics_normalised[rows, , drop = FALSE])


#Generating one heatmap per group
for (group_name in names(samples_submatrices)) {
  
  # Creating filename and printing it
  filename <- paste0("Graphs/", group_name, "_heatmap.png")
  print(paste("Saving:", filename))
  
  # Opening the PNG device to save the plot
  png(filename, width = 4000, height = 2000, res = 400)
  
  #Creating annotation colors for normalised zones
  max_val_palette <- max(samples_submatrices[[group_name]], na.rm = TRUE)
  
  n_low <- 1000
  n_high <- (max_val_palette - 1) * n_low
  
  palette_antibiotics_colours <- c(
    colorRampPalette(c("#990000", "#FFE5E5"))(n_low), #colouring resistant isolates in shades of red
    colorRampPalette(c("#E5FFE5", "#006600"))(n_high) #colouring sensitive isolates in shades of green
  )
  
  #Creating breaks to associate with thresholds of palette_antibiotic_colours and say red for values between 0 and 1, green for values > 1
  #Otherwise #990000 (i.e. red) going to be used for lowest normalised AST value (problem if lowest value != 0 in the dataset)
  breaks <- c(
    seq(0, 1, length.out = n_low),
    seq(1 + .Machine$double.eps, max_val_palette, length.out = n_high)
  )
  
  # Creating the heatmap
  pheatmap(samples_submatrices[[group_name]],
           cluster_rows = FALSE,  # Not clustering the isolates
           cluster_cols = FALSE,  # Not clustering the antibiotics
           annotation_row = annotation_row_all, #using annotations with toxins, fimbriae and visit (change all to combine to remove fimbriae)
           annotation_colors = annotation_colors_all, #using colours for toxins, fimbriae and visit (change all to combine to remove fimbriae)
           color = palette_antibiotics_colours,
           breaks = breaks,
           scale = "none",
           cellheight = 25,
           cellwidth = 30,
           fontsize_row = 8,
           fontsize_col = 15,
           fontsize = 10,
           main = paste(group_name, "heatmap"),
           legend_breaks = c(0, 1, 1.5),  # Antibiotics legend breaks
           legend_labels = c("Resistant (≤ 1)", "1.0", "Susceptible (> 1)"),
           angle_col = 45)  # Rotating antibiotics and toxins names
  
  # Closing the device (and saving the plot)
  dev.off()
}


#Closing all plots (just in case)
if (dev.cur() != 1) {
  dev.off()
}