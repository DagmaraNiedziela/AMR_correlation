library(tidyverse)
library(readxl)
#detach(package:plyr)    
#library(dplyr)

# Manipulate gene file ####
AMR_data <- read_excel("AMRFinder_all_results_short.xlsx", sheet = 1)
AMR_data

# Clean up the data to get only the columns I need - cfs_name, gene name, subclass, percent coverage and identity
AMR_data <- AMR_data %>% select(cfs_name, gene_symbol, percent_coverage_of_reference_sequence,percent_identity_to_reference_sequence,subclass)
View(AMR_data)
colnames(AMR_data) <- c("cfs_name", "gene_symbol", "percent_coverage","percent_identity","resistance")

typeof(AMR_data$resistance)
AMR_data$resistance <- factor(AMR_data$resistance)
levels(AMR_data$resistance) 
View(AMR_data)

#[1] "AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE" ??? gsub? Quinolone
#[2] "BETA-LACTAM"                                               
##[3] "CEPHALOSPORIN"                                             
##[4] "CHLORAMPHENICOL"                                           
#[5] "CHLORAMPHENICOL/FLORFENICOL"                              
##[6] "COLISTIN"                                                  
##[7] "EFFLUX"                                                    
##[8] "FOSFOMYCIN"                                                
##[9] "FOSMIDOMYCIN"                                              
#[10] "GENTAMICIN"                                                
#[11] "GENTAMICIN/KANAMYCIN/TOBRAMYCIN"    gsub!                       
##[12] "HYGROMYCIN"                                                
##[13] "KANAMYCIN"                                                 
##[14] "MACROLIDE"                                                 
##[15] "NITROFURANTOIN"                                            
#[16] "QUINOLONE"                                                 
#[17] "QUINOLONE/TRICLOSAN"     gsub!                                  
##[18] "STREPTOMYCIN"                                              
##[19] "STREPTOTHRICIN"                                            
#[20] "SULFONAMIDE"                                               
#[21] "TETRACYCLINE"                                              
#[22] "TRIMETHOPRIM"   

#"CEPHALOSPORIN", - amox-clav
AMR_data <- AMR_data %>% filter(!resistance %in% c("CHLORAMPHENICOL","COLISTIN","EFFLUX","FOSFOMYCIN","FOSMIDOMYCIN","HYGROMYCIN","KANAMYCIN","MACROLIDE","NITROFURANTOIN","STREPTOMYCIN", "STREPTOTHRICIN"))
AMR_data %>% filter(resistance == "AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE")

[1] "AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE"
[2] "BETA-LACTAM"                                               
[3] "CHLORAMPHENICOL/FLORFENICOL"                               
[4] "GENTAMICIN"                                                
[5] "GENTAMICIN/KANAMYCIN/TOBRAMYCIN"                           
[6] "QUINOLONE"                                                 
[7] "QUINOLONE/TRICLOSAN"                                       
[8] "SULFONAMIDE"                                               
[9] "TETRACYCLINE"                                              
[10] "TRIMETHOPRIM"  

# Now change all categories into the ones I need - with gsub 
?gsub
AMR_data$resistance <- gsub("AMPICILLIN/CHLORAMPHENICOL/QUINOLONE/RIFAMPIN/TETRACYCLINE","Quinolone",AMR_data$resistance)
AMR_data$resistance <- gsub("QUINOLONE/TRICLOSAN","Quinolone",AMR_data$resistance)
AMR_data$resistance <- gsub("QUINOLONE","Quinolone",AMR_data$resistance)

AMR_data$resistance <- gsub("CEPHALOSPORIN","Amoxicillin-clavulanic_acid",AMR_data$resistance) 
AMR_data$resistance <- gsub("BETA-LACTAM","Amoxicillin",AMR_data$resistance) 
AMR_data$resistance <- gsub("TETRACYCLINE","Tetracycline",AMR_data$resistance) 
AMR_data$resistance <- gsub("TRIMETHOPRIM","Trimetophrim-sulfonamide",AMR_data$resistance) 
AMR_data$resistance <- gsub("SULFONAMIDE","Trimetophrim-sulfonamide",AMR_data$resistance) 
AMR_data$resistance <- gsub("CHLORAMPHENICOL/FLORFENICOL","Florfenicol",AMR_data$resistance) 
AMR_data$resistance <- gsub("GENTAMICIN","Gentamicin",AMR_data$resistance) 
AMR_data$resistance <- gsub("GENTAMICIN/KANAMYCIN/TOBRAMYCIN","Gentamicin",AMR_data$resistance) 

# Change factor names to match 
#library(plyr)

#disk_data_long$Antibiotic <- mapvalues(disk_data_long$Antibiotic, from = c("amox_clav_average","amoxicillin_average","cefepime_average","cefotaxime_average","cefoxitin_average","cephalotin_average","colistin_average","florfenicol_average","gentamicin_average","tetracycline_average","trim_sulf_average"), 
                                       to = c("Amoxicillin-clavulanic_acid","Amoxicillin","Cefepime","Cefotaxime","Cefoxitin","Cephalothin","Colistin","Florfenicol","Gentamicin","Tetracycline","Trimetophrim-sulfonamide"))

#detach(package:plyr)

AMR_data_sorted <- AMR_data %>% arrange(resistance, gene_symbol) %>% group_by(resistance,cfs_name) %>% mutate(Genotype = paste0(gene_symbol, collapse = "_"))
View(AMR_data_sorted)
AMR_data_sorted %>% arrange(cfs_name) 

AMR_data_duplicates <- AMR_data_sorted %>% distinct(Genotype, .keep_all = TRUE) %>% arrange(cfs_name)
AMR_data_duplicates 
colnames(AMR_data_duplicates) 

# Copy all beta lactams to cover "Amoxicillin-clavulanic_acid","Amoxicillin","Cephalothin" 
AMR_data_betalactam <- AMR_data_duplicates %>% filter(resistance == "Amoxicillin")
AMR_data_betalactam2 <- AMR_data_betalactam

AMR_data_betalactam$resistance <- gsub("Amoxicillin","Amoxicillin-clavulanic_acid",AMR_data_betalactam$resistance) 

AMR_data_betalactam2$resistance <- gsub("Amoxicillin","Cephalothin",AMR_data_betalactam2$resistance) 

AMR_data_duplicates <- rbind(AMR_data_duplicates,AMR_data_betalactam,AMR_data_betalactam2)
View(AMR_data_duplicates)

# Copy all quinolones to cover "Ciprofloxacin","Nalixidic_acid","Marbofloxacin" 
AMR_data_quinolone <- AMR_data_duplicates %>% filter(resistance == "Quinolone")
AMR_data_quinolone2 <- AMR_data_quinolone
AMR_data_quinolone3 <- AMR_data_quinolone

AMR_data_quinolone$resistance <- gsub("Quinolone","Ciprofloxacin",AMR_data_quinolone$resistance) 

AMR_data_quinolone2$resistance <- gsub("Quinolone","Nalixidic_acid",AMR_data_quinolone2$resistance)

AMR_data_quinolone3$resistance <- gsub("Quinolone","Marbofloxacin",AMR_data_quinolone3$resistance) 

AMR_data_duplicates <- rbind(AMR_data_duplicates,AMR_data_quinolone,AMR_data_quinolone2,AMR_data_quinolone3)
View(AMR_data_duplicates) 

# Genotype summary #### 
View(AMR_data)

AMR_data_summary <- AMR_data %>% group_by(resistance, gene_symbol) %>% 
  summarise(isolate_count = n(), isolate_percent = isolate_count/143*100)
View(AMR_data_summary)

write.csv(AMR_data_summary, "Table1.csv")

# NOTE!!!! 
# Daniel's file contains 147 samples --- remove ones with bad assemblies 

# AMR_data_duplicates <- AMR_data_duplicates %>% filter(!Sample %in% c("CFS3241", "CFS3306", "CFS3348", "CFS3367", "Macrolide"))

# Add zone of inhibition ####

# Deleted non-E coli - isolate numbers 93, 125, 140 
# Need to figure out what CFS numbers were submitted to ENA - the new CFS numbers file has completely different numbers and they dont match disk diffusion 
disk_data <- read_excel("VTQ project database 2020-3-28.xlsx", sheet = "AST_replicates3")

View(disk_data)  

library(janitor)
disk_data <- disk_data %>% clean_names() 
colnames(disk_data) 

disk_data <- disk_data %>% select(cfs_number, contains("average"))
disk_data 
colnames(disk_data) 
[1] "cfs_number"             "gentamicin_average"     "amox_clav_average"      "cephalotin_average"    
[5] "cefoxitin_average"      "cefotaxime_average"     "cefepime_average"       "amoxicillin_average"   
[9] "tetracycline_average"   "colistin_average"       "florfenicol_average"    "trim_sulf_average"     
[13] "marbofloxacin_average"  "ciprofloxacin_average"  "nalidixic_acid_average"
nrow(disk_data) #149

#Filter for the samples we don't have 
disk_data <- disk_data %>% filter(!cfs_number %in% c("CFS3241", "CFS3306", "CFS3348", "CFS3367"))
nrow(disk_data) #143

disk_data <- disk_data[-c(144,145),]
# There were two blank rows with NAs only 

# Long format by antibiotic 
disk_data_long <- disk_data %>% pivot_longer(cols = -cfs_number, names_to = "Antibiotic", values_to = "Zone_of_inhibition")
View(disk_data_long) 
typeof(disk_data_long$Antibiotic)
disk_data_long$Antibiotic <- factor(disk_data_long$Antibiotic)
levels(disk_data_long$Antibiotic) 

[1] "amox_clav_average"      "amoxicillin_average"    "cefepime_average"       "cefotaxime_average"    
[5] "cefoxitin_average"      "cephalotin_average"     "ciprofloxacin_average"  "colistin_average"      
[9] "florfenicol_average"    "gentamicin_average"     "marbofloxacin_average"  "nalidixic_acid_average"
[13] "tetracycline_average"   "trim_sulf_average"      

colnames(AMR_data_duplicates)
colnames(disk_data_long)
levels(AMR_data_duplicates$resistance)
AMR_data_duplicates$resistance <- factor(AMR_data_duplicates$resistance)
levels(AMR_data_duplicates$resistance)
AMR_data_duplicates$resistance <- gsub("Gentamicin/KANAMYCIN/TOBRAMYCIN","Gentamicin",AMR_data_duplicates$resistance) 
# after fixing above, these are the levels 
[1] "Amoxicillin"                 "Amoxicillin-clavulanic_acid" "Cephalothin"                
[4] "Florfenicol"                 "Gentamicin"                  "Quinolone"                  
[7] "Tetracycline"                "Trimetophrim-sulfonamide" 
# Above add 3 different quinolones 
[1] "Amoxicillin"                 "Amoxicillin-clavulanic_acid" "Cephalothin"                
[4] "Ciprofloxacin"               "Florfenicol"                 "Gentamicin"                 
[7] "Marbofloxacin"               "Nalixidic_acid"              "Quinolone"                  
[10] "Tetracycline"                "Trimetophrim-sulfonamide"  

# Change factor names to match ####
library(plyr)
 
disk_data_long$Antibiotic <- mapvalues(disk_data_long$Antibiotic, from = c("amox_clav_average","amoxicillin_average","cefepime_average","cefotaxime_average","cefoxitin_average","cephalotin_average","colistin_average","ciprofloxacin_average","florfenicol_average","gentamicin_average","marbofloxacin_average","nalidixic_acid_average","tetracycline_average","trim_sulf_average"), 
          to = c("Amoxicillin-clavulanic_acid","Amoxicillin","Cefepime","Cefotaxime","Cefoxitin","Cephalothin","Colistin","Ciprofloxacin","Florfenicol","Gentamicin","Marbofloxacin","Nalixidic_acid","Tetracycline","Trimetophrim-sulfonamide"))
detach(package:plyr)
colnames(disk_data_long)
colnames(disk_data_long) <- c("cfs_name", "resistance", "Zone_of_inhibition") 

#disk_data_long <- disk_data_long %>% filter(!Sample %in% c("CFS3241", "CFS3306", "CFS3348", "CFS3367"))

View(disk_data_long)

# Join both tables ####

correlation <- full_join(AMR_data_duplicates, disk_data_long, by=c("cfs_name", "resistance"))
#levels(disk_data_long$Resistance)
#levels(AMR_data_duplicates$Resistance) 
# Warning message:
# Column `Resistance` joining factors with different levels, coercing to character vector 

correlation
View(correlation)
nrow(correlation) #2057 

#correlation_genes <- correlation[complete.cases(correlation),] #Don't do this, just make NA "No genes"
View(correlation_genes)
#nrow(correlation_genes) #279 

correlation_genes <- correlation %>% mutate(Genotype = replace_na(Genotype, "No genes"))

write.csv(correlation_genes, "correlation_genes_new_data.csv")

# Round zones of inhibition 
correlation_genes$Zone_of_inhibition <- round(correlation_genes$Zone_of_inhibition)

# Plot data per gene #### 

# ** Gentamicin ####

# New plot with antibiotics removed 
correlation_genes_gentamicin <- subset(correlation_genes, resistance == "Gentamicin")
#correlation_genes_gentamicin1$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_gentamicin$Genotype)
#correlation_genes_gentamicin1 <- correlation_genes_gentamicin %>% 
#  filter(str_detect(Genotype, 'aac\\(3\\)-IIe|aac\\(3\\)-IVa|ant\\(2\'\'\\)-Ia|No genes')) 
#nrow(correlation_genes_gentamicin1) #111
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
#correlation_genes_gentamicin2 <- correlation_genes_gentamicin %>% 
#  filter(!str_detect(Genotype, 'aac\\(3\\)-IIa|aac\\(3\\)-IVa|ant\\(2\'\'\\)-Ia|No genes')) 
#nrow(correlation_genes_gentamicin2) #32 

#correlation_genes_gentamicin2 <- correlation_genes_gentamicin2 %>% 
#  mutate(Genotype = paste("No genes"))
#correlation_genes_gentamicin3 <- rbind(correlation_genes_gentamicin1,correlation_genes_gentamicin2)

#nrow(correlation_genes_gentamicin3) #143 

# There is no need to fix gentamicin's gene names because there are only the necessary 3 

# Text under x axis 
library(grid)
text_resistant <- textGrob("R<", gp=gpar(fontsize=13, fontface="bold"))
text_susceptible <- textGrob(">S", gp=gpar(fontsize=13, fontface="bold"))

gentamicin_DD <- ggplot(data = correlation_genes_gentamicin) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 15, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
 # theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=15,xmax=15,ymin=-0.2,ymax=-0.2) + 
 # annotation_custom(text_resistant,xmin=12,xmax=12,ymin=-0.2,ymax=-0.2) + 
 # coord_cartesian(clip = "off")

gentamicin_DD
ggsave("Gentamicin_DD.jpeg") 


# Plot without position_dodge and count instead of percentages - I don't like this 

ggplot(data = correlation_genes_gentamicin) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 15, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(count), fill = Genotype, group = Genotype), width = 1) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(count), group = Genotype, label = stat(count)), stat = "count", position=position_stack(0.5)) + 
  xlab("Zone of inhibition (cm)") + ylab("Number of isolates with a genotype") + 
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=15,xmax=15,ymin=-7,ymax=-7) + 
  annotation_custom(text_resistant,xmin=12,xmax=12,ymin=-7,ymax=-7) + 
  coord_cartesian(clip = "off") 

# Note! Had to lower ymin and ymax in annotation_custom to put text below graph 

ggsave("Gentamicin_no_dodge.jpeg") 

?ggsave # Pulls up the help section 

# ** Cephalothin ####
# New plot with antibiotics removed 
correlation_genes_cephalothin <- subset(correlation_genes, resistance == "Cephalothin")
correlation_genes_cephalothin$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_gentamicin$Genotype)
correlation_genes_cephalothin1 <- correlation_genes_cephalothin %>% 
  filter(str_detect(Genotype, 'blaTEM-1|No genes')) 
nrow(correlation_genes_cephalothin1) #54
#blaTEM-1A	blaTEM-1B	blaTEM-1C	blaTEM-1D
# In new one only bla_TEM-1
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_genes_cephalothin2 <- correlation_genes_cephalothin %>% 
  filter(!str_detect(Genotype, 'blaTEM-1|No genes')) 
nrow(correlation_genes_cephalothin2) #89 

correlation_genes_cephalothin2 <- correlation_genes_cephalothin2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_cephalothin3 <- rbind(correlation_genes_cephalothin1,correlation_genes_cephalothin2)

nrow(correlation_genes_cephalothin3) #143 

# Get rid of blaTEM135 with gsub 
correlation_genes_cephalothin3$Genotype <- gsub("blaEC_blaTEM-135","No genes",correlation_genes_cephalothin3$Genotype)

cephalotin_DD <- ggplot(data = correlation_genes_cephalothin3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
 # theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
#  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off")

cephalotin_DD
ggsave("Cephalothin_DD.jpeg") 

# ** Cefoxitin - NOT included ####

# New plot with antibiotics removed 
correlation_genes_cefoxitin <- subset(correlation_genes, Resistance == "Cefoxitin")
correlation_genes_cefoxitin1$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_cefoxitin$Genotype)
correlation_genes_cefoxitin1 <- correlation_genes_cefoxitin %>% 
  filter(str_detect(Genotype, 'blaCMY-2|No genes')) 

#blaCMY-2
nrow(correlation_genes_cefoxitin1) #81
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_genes_cefoxitin2 <- correlation_genes_cefoxitin %>% 
  filter(!str_detect(Genotype, 'blaCMY-2|No genes')) 
nrow(correlation_genes_cefoxitin2) #62

correlation_genes_cefoxitin2 <- correlation_genes_cefoxitin2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_cefoxitin3 <- rbind(correlation_genes_cefoxitin1,correlation_genes_cefoxitin2)

nrow(correlation_genes_cefoxitin3) #143 

ggplot(data = correlation_genes_cefoxitin3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefoxitin_new.jpeg") 




#old plot
plotC <- ggplot(data = subset(correlation_genes, Resistance == "Cefoxitin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefoxitin_rects_DN.jpeg") 
# Dodged numbers escape to the right 

# ** Cefotaxime - NOT included ####
# New plot with antibiotics removed 
correlation_genes_cefotaxime <- subset(correlation_genes, Resistance == "Cefotaxime")
correlation_genes_cefotaxime1$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_cefotaxime$Genotype)
correlation_genes_cefotaxime1 <- correlation_genes_cefotaxime %>% 
  filter(str_detect(Genotype, 'blaCMY-2|blaCTX-M-1_|blaTEM-106|No genes')) 
#blaCMY-2|blaCTX-M-1|blaTEM-106

nrow(correlation_genes_cefotaxime1) #84
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_genes_cefotaxime2 <- correlation_genes_cefotaxime %>% 
  filter(!str_detect(Genotype, 'blaCMY-2|blaCTX-M-1_|blaTEM-106|No genes')) 
nrow(correlation_genes_cefotaxime2) #59 

correlation_genes_cefotaxime2 <- correlation_genes_cefotaxime2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_cefotaxime3 <- rbind(correlation_genes_cefotaxime1,correlation_genes_cefotaxime2)

nrow(correlation_genes_cefotaxime3) #143 

ggplot(data = correlation_genes_cefotaxime3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 22, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 26, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=26,xmax=26,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=22,xmax=22,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefotaxime_new.jpeg") 

#old plot
plotD <- ggplot(data = subset(correlation_genes, Resistance == "Cefotaxime")) + 
  annotate(geom = "rect", xmin = -Inf, xmax =22, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 26, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=26,xmax=26,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=22,xmax=22,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefotaxime_rects_DN.jpeg") 

# ** Cefepime - NOT included ####

# New plot with antibiotics removed 
correlation_genes_cefepime <- subset(correlation_genes, Resistance == "Cefepime")
correlation_genes_cefepime1$Genotype
typeof(correlation_genes_cefepime$Genotype)
correlation_genes_cefepime1 <- correlation_genes_cefepime %>% 
  filter(str_detect(Genotype, 'blaCTX-M-1_|blaTEM-106|blaOXA-1|No genes'))
#blaCTX-M-1 blaTEM-106 blaOXA-1

nrow(correlation_genes_cefepime1) #88 

correlation_genes_cefepime2 <- correlation_genes_cefepime %>% 
  filter(!str_detect(Genotype, 'blaCTX-M-1_|blaTEM-106|blaOXA-1|No genes')) 
nrow(correlation_genes_cefepime2) #55

correlation_genes_cefepime2 <- correlation_genes_cefepime2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_cefepime3 <- rbind(correlation_genes_cefepime1,correlation_genes_cefepime2)

nrow(correlation_genes_cefepime3) #143 

ggplot(data = correlation_genes_cefepime3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "green", alpha = 0.1) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefepime_new.jpeg") 

#old plot
plotE <- ggplot(data = subset(correlation_genes, Resistance == "Cefepime")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Cefepime_rects_DN.jpeg") 


# ** Amoxicillin #### 

# blaCTX-M-1	blaOXA-1	blaOXA-2	blaTEM	blaTEM-1	blaTEM-135	blaTEM-30	blaTEM-33


correlation_genes_amoxicillin <- subset(correlation_genes, resistance == "Amoxicillin")
correlation_genes_amoxicillin$Genotype
typeof(correlation_genes_amoxicillin$Genotype)
correlation_genes_amoxicillin1 <- correlation_genes_amoxicillin %>% 
  filter(str_detect(Genotype, 'blaCTX-M-1|blaOXA-1|blaOXA-2|blaTEM|blaTEM-30|blaTEM-33|blaTEM-135|blaTEM-1|No genes'))

nrow(correlation_genes_amoxicillin1) #65 

correlation_genes_amoxicillin2 <- correlation_genes_amoxicillin %>% 
  filter(!str_detect(Genotype, 'blaCTX-M-1|blaOXA-1|blaOXA-2|blaTEM|blaTEM-30|blaTEM-33|blaTEM-135|blaTEM-1|No genes')) 
nrow(correlation_genes_amoxicillin2) #78

correlation_genes_amoxicillin2 <- correlation_genes_amoxicillin2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_amoxicillin3 <- rbind(correlation_genes_amoxicillin1,correlation_genes_amoxicillin2)

nrow(correlation_genes_amoxicillin3) #143 

library(ggrepel)
amoxicillin_DD <- ggplot(data = correlation_genes_amoxicillin3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 16, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 17, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "green", alpha = 0.1) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
  
#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=17,xmax=17,ymin=-0.2,ymax=-0.2) + 
 # annotation_custom(text_resistant,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off")

amoxicillin_DD
ggsave("Amoxicillin_DD.jpeg")
# Saving 11.5 x 7.58 in image


# ** Amoxicillin clavulanic acid #### 

# ampC_C-42T	blaCMY-2	blaOXA-1	blaOXA-2	blaTEM-30

# New plot with antibiotics removed 
correlation_genes_amoxicillin_clavulanic_acid <- subset(correlation_genes, resistance == "Amoxicillin-clavulanic_acid")
correlation_genes_amoxicillin_clavulanic_acid$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_amoxicillin_clavulanic_acid$Genotype)
correlation_genes_amoxicillin_clavulanic_acid1 <- correlation_genes_amoxicillin_clavulanic_acid %>% 
  filter(str_detect(Genotype, 'blaCMY-2|blaOXA-1|blaOXA-2|blaTEM-30|ampC_C-42T|No genes')) 

nrow(correlation_genes_amoxicillin_clavulanic_acid1) #11
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_genes_amoxicillin_clavulanic_acid2 <- correlation_genes_amoxicillin_clavulanic_acid %>% 
  filter(!str_detect(Genotype, 'blaCMY-2|blaOXA-1|blaOXA-2|blaTEM-30|ampC_C-42T|No genes')) 
nrow(correlation_genes_amoxicillin_clavulanic_acid2) #137 
# I think here lies the issue, or previously in adding amox clav twice???? 

correlation_genes_amoxicillin_clavulanic_acid2 <- correlation_genes_amoxicillin_clavulanic_acid2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_amoxicillin_clavulanic_acid3 <- rbind(correlation_genes_amoxicillin_clavulanic_acid1,correlation_genes_amoxicillin_clavulanic_acid2)

nrow(correlation_genes_amoxicillin_clavulanic_acid3) #148 --- why not 143?????  
View(correlation_genes_amoxicillin_clavulanic_acid3) 
complete.cases(correlation_genes_amoxicillin_clavulanic_acid3)

correlation_genes_amoxicillin_clavulanic_acid3 <- correlation_genes_amoxicillin_clavulanic_acid3 %>% distinct(cfs_name, .keep_all = TRUE)

amox_clav_DD <- ggplot(data = correlation_genes_amoxicillin_clavulanic_acid3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 13, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # theme(panel.background=element_rect(fill="white", colour="black")) +
 # annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
#  annotation_custom(text_resistant,xmin=13,xmax=13,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off")
# Here only the R and S text works, why??? 
amox_clav_DD

ggsave("Amoxicillin-clavulanic_acid_DD.jpeg") 

# ** Tetracycline #### 

# Can play with colors here: https://www.w3schools.com/colors/colors_picker.asp 
# Change color in "fill" command of geom_rect
# Alpha in geom_rect code is transparency of the color, sometimes if it is set too low the color disappears 
# The transparency is important to see the grid lines underneath 

tetracycline_DD <- ggplot(data = subset(correlation_genes, resistance == "Tetracycline")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
 # theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=19,xmax=19,ymin=-0.2,ymax=-0.2) + 
 # annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
 # coord_cartesian(clip = "off")

tetracycline_DD
ggsave("Tetracycline_DD.jpeg") 

# Plot without dodged columns  

ggplot(data = subset(correlation_genes, Resistance == "Tetracycline")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(count), fill = Genotype, group = Genotype), width = 1) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(count), group = Genotype, label = stat(count)), stat = "count", position=position_stack(0.5)) + 
  xlab("Zone of inhibition (cm)") + ylab("Number of isolates with a genotype") + 
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=19,xmax=19,ymin=-6,ymax=-6) + 
  annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-6,ymax=-6) + 
  coord_cartesian(clip = "off")

ggsave("Tetracycline_no_dodge.jpeg") 


# ** Colistin - NOT #### 
# There are only two positive isolates for that resistance gene, and susceptible isolates exist 
# The cutoffs for susceptible and resistant are correct 
# Which means that the S and R data from the master file are not correct 
# The 22 seem to be susceptible but the dodge made it look bad 

ggplot(data = subset(correlation_genes, Resistance == "Colistin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 11, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 14, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(text_susceptible,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
  annotation_custom(text_resistant,xmin=11,xmax=11,ymin=-0.2,ymax=-0.2) + 
  coord_cartesian(clip = "off")

ggsave("Colistin_new.jpeg") 

# ** Florfenicol #### 

# New plot, filtered for genes 
correlation_genes_florfenicol <- subset(correlation_genes, resistance == "Florfenicol")
correlation_genes_florfenicol$Genotype
#typeof(correlation_genes_florfenicol$Genotype)
#correlation_genes_florfenicol1 <- correlation_genes_florfenicol %>% 
#  filter(str_detect(Genotype, 'floR|No genes')) 
#nrow(correlation_genes_florfenicol1) #119

#correlation_genes_florfenicol2 <- correlation_genes_florfenicol %>% 
#  filter(!str_detect(Genotype, 'floR|No genes')) 
#nrow(correlation_genes_florfenicol2) #24

#correlation_genes_florfenicol2 <- correlation_genes_florfenicol2 %>% 
 # mutate(Genotype = paste("No genes"))
#correlation_genes_florfenicol3 <- rbind(correlation_genes_florfenicol1,correlation_genes_florfenicol2)

#nrow(correlation_genes_florfenicol3) #143 

florfenicol_DD <- ggplot(data = correlation_genes_florfenicol) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=19,xmax=19,ymin=-0.2,ymax=-0.2) + 
 # annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off") 

florfenicol_DD
ggsave("Florfenicol_DD.jpeg") 

# ** Trimethoprim-Sulfamethoxazole #### 

trim_sulph_DD <- ggplot(data = subset(correlation_genes, resistance == "Trimetophrim-sulfonamide")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 10, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 16, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
 # theme(plot.margin = unit(c(1,1,2,1), "lines")) +
 # annotation_custom(text_susceptible,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
 # annotation_custom(text_resistant,xmin=10,xmax=10,ymin=-0.2,ymax=-0.2) + 
 # coord_cartesian(clip = "off")

trim_sulph_DD
ggsave("Trimetophrim-sulfonamide_DD.jpeg") 

# ** Marbofloxacin #### 

marbo_DD <- ggplot(data = subset(correlation_genes, resistance == "Marbofloxacin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 20, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
# theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# annotation_custom(text_susceptible,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
# annotation_custom(text_resistant,xmin=10,xmax=10,ymin=-0.2,ymax=-0.2) + 
# coord_cartesian(clip = "off")

marbo_DD
ggsave("Marbofloxacin_DD.jpeg") 

# ** Nalixidic acid #### 

nalix_DD <- ggplot(data = subset(correlation_genes, resistance == "Nalixidic_acid")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 13, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
# theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# annotation_custom(text_susceptible,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
# annotation_custom(text_resistant,xmin=10,xmax=10,ymin=-0.2,ymax=-0.2) + 
# coord_cartesian(clip = "off")

nalix_DD
ggsave("Nalixidic_acid_DD.jpeg") 

# ** Ciprofloxacin #### 

cipro_DD <- ggplot(data = subset(correlation_genes, resistance == "Ciprofloxacin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 15, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()
# theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# annotation_custom(text_susceptible,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
# annotation_custom(text_resistant,xmin=10,xmax=10,ymin=-0.2,ymax=-0.2) + 
# coord_cartesian(clip = "off")

cipro_DD
ggsave("Ciprofloxacin_DD.jpeg") 


# Arrange the plots together #### 
library(ggpubr)

plot <- ggarrange(amoxicillin_DD,amox_clav_DD,cephalotin_DD,
                  florfenicol_DD,gentamicin_DD,tetracycline_DD,
                  trim_sulph_DD,
                  labels = c("A", "B", "C", "D", "E", "F", "G"), 
                  nrow = 4, ncol = 2,align = "v") 

plot
ggsave("Combined_plot_DD.jpeg",scale = 2)
# Potentially make No genes always the same colour? 

plot1 <- ggarrange(amoxicillin_DD,amox_clav_DD,cephalotin_DD,
                  labels = c("A", "B", "C"), 
                  nrow = 3, ncol = 1,align = "v") 
plot1
#plot1 + theme(axis.title = element_text(size=16), axis.text.y = element_text(size=16),
 #         axis.text.x = element_text(size=16))
ggsave("Combined_plot_beta_lactams.jpeg",width = 7, height=11, units = "in")
ggsave("Combined_plot_beta_lactams2.jpeg")
ggsave("Combined_plot_beta_lactams3.jpeg",width = 8, height=11, units = "in")
# Original Saving 7.31 x 10.8 in image - make higher 

# For this one use Cairo - nope, changes nothing 
library(Cairo)
Cairo(file="Combined_plot_beta_lactams.png", 
      type="png",
      units="in", 
      width=7, 
      height=10, 
      pointsize=6, 
      dpi=600)
plot1;
dev.off() 
?Cairo

plot2 <- ggarrange(florfenicol_DD,gentamicin_DD,tetracycline_DD,
                   trim_sulph_DD,
                   labels = c("A", "B", "C", "D"), 
                   nrow = 2, ncol = 2,align = "v") 

plot2
ggsave("Combined_plot_others.jpeg",scale =2) 
ggsave("Combined_plot_others2.jpeg",scale =1.5) 
ggsave("Combined_plot_others3.jpeg") 
# NEW - Saving 14.2 x 9.53 in image
# Saving 13.6 x 7.92 in image
# Gentamicin, Tetracycline, Colistin, Florfenicol, Trim_sulf

plot3 <- ggarrange(cipro_DD,nalix_DD,marbo_DD,
                   labels = c("A", "B", "C"), 
                   nrow = 3, ncol = 1,align = "v",common.legend = TRUE,
                   legend = "right") 
plot3
#plot1 + theme(axis.title = element_text(size=16), axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=16))
ggsave("Combined_plot_quinolones.jpeg",width = 9, height=11, units = "in")
ggsave("Combined_plot_quinolones2.jpeg") # Saving 10.4 x 11.7 in image
# Can do a white plot bakcground instead of the grey 
# Also get the plot area to be a bit higher so the numbers are visible 


plot4 <- ggarrange(ggarrange(amoxicillin_DD,amox_clav_DD,cephalotin_DD,
                             labels = c("A", "B", "C"), 
                             nrow = 3, align = "v"),
                   ggarrange(cipro_DD,nalix_DD,marbo_DD,
                             labels = c("D", "E", "F"), 
                             nrow = 3, align = "v",common.legend = TRUE,
                             legend = "right"),ncol = 2, widths = c(0.9,1))
plot4 
ggsave("Combined_DD_betalactams_and_quinolones.jpeg")


# NEW! Add MIC ####
library(tidyverse)
library(readxl)

# Deleted non-E coli - isolate numbers 93, 125, 140 
# Need to figure out what CFS numbers were submitted to ENA - the new CFS numbers file has completely different numbers and they dont match disk diffusion 
MIC_data <- read_excel("MIC.xlsx", sheet = 1)
# hidden columns!!!!
View(MIC_data)  
MIC_data <- MIC_data[,c(1:7)]

#library(janitor)
#MIC_data <- MIC_data %>% clean_names() 
colnames(MIC_data) 

MIC_data 
colnames(MIC_data) 
[1] "cfs_name"                      "Gentamicin"                   
[3] "Cephalotin"                    "Tetracycline"                 
[5] "Marbofloxacin"                 "Ciprofloxacin"                
[7] "Trimethoprim-sulfamethoxazole"
nrow(MIC_data) #104 

#Filter for the samples we don't have 
MIC_data <- MIC_data %>% filter(!cfs_name %in% c("CFS3241", "CFS3306", "CFS3348", "CFS3367"))
nrow(MIC_data) #104

colnames(MIC_data) <- c("cfs_name","Gentamicin","Cephalothin","Tetracycline",
                        "Marbofloxacin","Ciprofloxacin","Trimetophrim-sulfonamide")

# Long format by antibiotic 
MIC_data_long <- MIC_data %>% pivot_longer(cols = -cfs_name, names_to = "Antibiotic", values_to = "MIC")
View(MIC_data_long) 
typeof(MIC_data_long$Antibiotic)
MIC_data_long$Antibiotic <- factor(MIC_data_long$Antibiotic)
levels(MIC_data_long$Antibiotic) 
[1] "Cephalotin"                    "Ciprofloxacin"                
[3] "Gentamicin"                    "Marbofloxacin"                
[5] "Tetracycline"                  "Trimetophrim-sulfamethoxazole"
# "Trimethoprim-sulfonamide"
Trimetophrim-sulfonamide

colnames(MIC_data_long)
levels(AMR_data_duplicates$resistance)
View(AMR_data_duplicates)

colnames(MIC_data_long) <- c("cfs_name", "resistance", "MIC") 

#MIC_data_long <- MIC_data_long %>% filter(!Sample %in% c("CFS3241", "CFS3306", "CFS3348", "CFS3367"))

View(MIC_data_long)

# Join both tables ####

correlation_MIC <- full_join(AMR_data_duplicates, MIC_data_long, by=c("cfs_name", "resistance"))
#levels(MIC_data_long$Resistance)
#levels(AMR_data_duplicates$Resistance) 
# Warning message:
# Column `Resistance` joining factors with different levels, coercing to character vector 

correlation_MIC
View(correlation_MIC)
nrow(correlation_MIC) #1076 

correlation_genes_MIC <- correlation_MIC %>% mutate(Genotype = replace_na(Genotype, "No genes"))

#write.csv(correlation_genes, "correlation_genes_new_data.csv")

# Round MIC? 
#correlation_genes_MIC$MIC <- round(correlation_genes$Zone_of_inhibition)

# From MIC 
Trimethophrim-sulfonamide 
# From AMR_gene data
Trimetophrim-sulfonamide
# OK I think its done - changed in MIC 

# Plot data per gene #### 

# ** Gentamicin ####

# New plot with antibiotics removed 
correlation_gentamicin_MIC <- subset(correlation_genes_MIC, resistance == "Gentamicin")

# There is no need to fix gentamicin's gene names because there are only the necessary 3 

# Text under x axis 
library(grid)
text_resistant <- textGrob("R<", gp=gpar(fontsize=13, fontface="bold"))
text_susceptible <- textGrob(">S", gp=gpar(fontsize=13, fontface="bold"))

gentamicin_MIC <- ggplot(data = correlation_gentamicin_MIC) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 4, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 16, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()

gentamicin_MIC
ggsave("Gentamicin_MIC.jpeg") 

# ** Cephalothin ####
# New plot with antibiotics removed 
# Cephalothin was misspelled in MIC file - Cephalotin - now have to rerun all 

correlation_cephalothin_MIC <- subset(correlation_genes_MIC, resistance == "Cephalothin")
correlation_cephalothin_MIC$Genotype # check if it worked and our genes are there 

correlation_cephalothin_MIC1 <- correlation_cephalothin_MIC %>% 
  filter(str_detect(Genotype, 'blaTEM-1|No genes')) 
nrow(correlation_genes_cephalothin1) #54
#blaTEM-1A	blaTEM-1B	blaTEM-1C	blaTEM-1D
# In new one only bla_TEM-1
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_cephalothin_MIC2 <- correlation_cephalothin_MIC %>% 
  filter(!str_detect(Genotype, 'blaTEM-1|No genes')) 
nrow(correlation_cephalothin_MIC2) #89 

correlation_cephalothin_MIC2 <- correlation_cephalothin_MIC2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_cephalothin_MIC3 <- rbind(correlation_cephalothin_MIC1,correlation_cephalothin_MIC2)

nrow(correlation_cephalothin_MIC3) #165!!!! 
# Some rows of NAs 
View(correlation_cephalothin_MIC3)

# Get rid of blaTEM135 with gsub 
correlation_cephalothin_MIC3$Genotype <- gsub("blaEC_blaTEM-135","No genes",correlation_cephalothin_MIC3$Genotype)

cephalotin_MIC <- ggplot(data = correlation_cephalothin_MIC3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 8, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 32, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()

cephalotin_MIC
ggsave("Cephalothin_MIC.jpeg") 

# ** Amoxicillin - NOT #### 

# blaCTX-M-1	blaOXA-1	blaOXA-2	blaTEM	blaTEM-1	blaTEM-135	blaTEM-30	blaTEM-33


correlation_genes_amoxicillin <- subset(correlation_genes, resistance == "Amoxicillin")
correlation_genes_amoxicillin$Genotype
typeof(correlation_genes_amoxicillin$Genotype)
correlation_genes_amoxicillin1 <- correlation_genes_amoxicillin %>% 
  filter(str_detect(Genotype, 'blaCTX-M-1|blaOXA-1|blaOXA-2|blaTEM|blaTEM-30|blaTEM-33|blaTEM-135|blaTEM-1|No genes'))

nrow(correlation_genes_amoxicillin1) #65 

correlation_genes_amoxicillin2 <- correlation_genes_amoxicillin %>% 
  filter(!str_detect(Genotype, 'blaCTX-M-1|blaOXA-1|blaOXA-2|blaTEM|blaTEM-30|blaTEM-33|blaTEM-135|blaTEM-1|No genes')) 
nrow(correlation_genes_amoxicillin2) #78

correlation_genes_amoxicillin2 <- correlation_genes_amoxicillin2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_amoxicillin3 <- rbind(correlation_genes_amoxicillin1,correlation_genes_amoxicillin2)

nrow(correlation_genes_amoxicillin3) #143 

amoxicillin_DD <- ggplot(data = correlation_genes_amoxicillin3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 16, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 17, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "green", alpha = 0.1) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_bw()

#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# annotation_custom(text_susceptible,xmin=17,xmax=17,ymin=-0.2,ymax=-0.2) + 
# annotation_custom(text_resistant,xmin=16,xmax=16,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off")

amoxicillin_DD
ggsave("Amoxicillin_DD.jpeg")


# ** Amoxicillin clavulanic acid - NOT #### 

# ampC_C-42T	blaCMY-2	blaOXA-1	blaOXA-2	blaTEM-30

# New plot with antibiotics removed 
correlation_genes_amoxicillin_clavulanic_acid <- subset(correlation_genes, resistance == "Amoxicillin-clavulanic_acid")
correlation_genes_amoxicillin_clavulanic_acid$Genotype # check if it worked and our genes are there 
typeof(correlation_genes_amoxicillin_clavulanic_acid$Genotype)
correlation_genes_amoxicillin_clavulanic_acid1 <- correlation_genes_amoxicillin_clavulanic_acid %>% 
  filter(str_detect(Genotype, 'blaCMY-2|blaOXA-1|blaOXA-2|blaTEM-30|ampC_C-42T|No genes')) 

nrow(correlation_genes_amoxicillin_clavulanic_acid1) #11
# | signs stands for "or" 
# \\ to escape double brackets 
# \ to escape single quotations
# str_detect will take the characters and keep all phrases that contain these 
correlation_genes_amoxicillin_clavulanic_acid2 <- correlation_genes_amoxicillin_clavulanic_acid %>% 
  filter(!str_detect(Genotype, 'blaCMY-2|blaOXA-1|blaOXA-2|blaTEM-30|ampC_C-42T|No genes')) 
nrow(correlation_genes_amoxicillin_clavulanic_acid2) #137 
# I think here lies the issue, or previously in adding amox clav twice???? 

correlation_genes_amoxicillin_clavulanic_acid2 <- correlation_genes_amoxicillin_clavulanic_acid2 %>% 
  mutate(Genotype = paste("No genes"))
correlation_genes_amoxicillin_clavulanic_acid3 <- rbind(correlation_genes_amoxicillin_clavulanic_acid1,correlation_genes_amoxicillin_clavulanic_acid2)

nrow(correlation_genes_amoxicillin_clavulanic_acid3) #148 --- why not 143?????  
View(correlation_genes_amoxicillin_clavulanic_acid3) 
complete.cases(correlation_genes_amoxicillin_clavulanic_acid3)

correlation_genes_amoxicillin_clavulanic_acid3 <- correlation_genes_amoxicillin_clavulanic_acid3 %>% distinct(cfs_name, .keep_all = TRUE)

amox_clav_DD <- ggplot(data = correlation_genes_amoxicillin_clavulanic_acid3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 13, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_bw()
#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# theme(panel.background=element_rect(fill="white", colour="black")) +
# annotation_custom(text_susceptible,xmin=18,xmax=18,ymin=-0.2,ymax=-0.2) + 
#  annotation_custom(text_resistant,xmin=13,xmax=13,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off")
# Here only the R and S text works, why??? 
amox_clav_DD

ggsave("Amoxicillin-clavulanic_acid_DD.jpeg") 

# ** Tetracycline #### 

# Can play with colors here: https://www.w3schools.com/colors/colors_picker.asp 
# Change color in "fill" command of geom_rect
# Alpha in geom_rect code is transparency of the color, sometimes if it is set too low the color disappears 
# The transparency is important to see the grid lines underneath 

tetracycline_MIC <- ggplot(data = subset(correlation_genes_MIC, resistance == "Tetracycline")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 4, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 16, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()

tetracycline_MIC
ggsave("Tetracycline_MIC.jpeg") 

# ** Florfenicol - NOT #### 

# New plot, filtered for genes 
correlation_genes_florfenicol <- subset(correlation_genes, resistance == "Florfenicol")
correlation_genes_florfenicol$Genotype
#typeof(correlation_genes_florfenicol$Genotype)
#correlation_genes_florfenicol1 <- correlation_genes_florfenicol %>% 
#  filter(str_detect(Genotype, 'floR|No genes')) 
#nrow(correlation_genes_florfenicol1) #119

#correlation_genes_florfenicol2 <- correlation_genes_florfenicol %>% 
#  filter(!str_detect(Genotype, 'floR|No genes')) 
#nrow(correlation_genes_florfenicol2) #24

#correlation_genes_florfenicol2 <- correlation_genes_florfenicol2 %>% 
# mutate(Genotype = paste("No genes"))
#correlation_genes_florfenicol3 <- rbind(correlation_genes_florfenicol1,correlation_genes_florfenicol2)

#nrow(correlation_genes_florfenicol3) #143 

florfenicol_DD <- ggplot(data = correlation_genes_florfenicol) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 14, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.2) + 
  geom_bar(mapping = aes(x = Zone_of_inhibition, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = Zone_of_inhibition, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_bw()
#  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
# annotation_custom(text_susceptible,xmin=19,xmax=19,ymin=-0.2,ymax=-0.2) + 
# annotation_custom(text_resistant,xmin=14,xmax=14,ymin=-0.2,ymax=-0.2) + 
#  coord_cartesian(clip = "off") 

florfenicol_DD
ggsave("Florfenicol_DD.jpeg") 

# ** Trimethoprim-Sulfamethoxazole #### 
# This one also misspelled, this time in antibiotics! Just rename the column
View(subset(correlation_genes_MIC, resistance == "Trimetophrim-sulfonamide"))

trim_sulph_MIC <- ggplot(data = subset(correlation_genes_MIC, resistance == "Trimetophrim-sulfonamide")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 38, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 76, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent,expand = c(0.05, 0.05)) +
  scale_x_continuous(breaks = c(19,38,57,76), labels = c("1/19","2/38","3/57","4/76"))+
  theme_bw()

trim_sulph_MIC
ggsave("Trimetophrim-sulfonamide_MIC.jpeg") 

# ** Marbofloxacin #### 

marbo_MIC <- ggplot(data = subset(correlation_genes_MIC, resistance == "Marbofloxacin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent, expand = c(0.05, 0.05)) +
  theme_bw()

marbo_MIC
ggsave("Marbofloxacin_MIC.jpeg") 

# ** Nalixidic acid - NOT #### 

View(subset(correlation_genes_MIC, resistance == "Nalixidic_acid"))
nalix_MIC <- ggplot(data = subset(correlation_genes_MIC, resistance == "Nalixidic_acid")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 13, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 19, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("Zone of inhibition (cm)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_bw()

nalix_MIC
ggsave("Nalixidic_acid_MIC.jpeg") 

# ** Ciprofloxacin #### 

cipro_MIC <- ggplot(data = subset(correlation_genes_MIC, resistance == "Ciprofloxacin")) + 
  annotate(geom = "rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#00CED1", alpha = 0.1) + 
  annotate(geom = "rect", xmin = 4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#DC143C", alpha = 0.2) + 
  geom_bar(mapping = aes(x = MIC, y = stat(prop), fill = Genotype, group = Genotype), width = 1, position = position_dodge(preserve = "single")) +
  geom_text_repel(mapping = aes(x = MIC, y = stat(prop), group = Genotype, label = stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + 
  xlab("MIC (μg/mL)") + ylab("Percentage of isolates with a genotype") + 
  scale_y_continuous(labels = scales::percent,expand = c(0.05, 0.05)) +
  theme_bw()

cipro_MIC
ggsave("Ciprofloxacin_MIC.jpeg") 


# Arrange the plots together #### 
library(ggpubr)

plot_MIC <- ggarrange(cephalotin_MIC,gentamicin_MIC,tetracycline_MIC,
                  trim_sulph_MIC,cipro_MIC,marbo_MIC,
                  labels = c("A", "B", "C", "D", "E", "F"), 
                  nrow = 3, ncol = 2,align = "v") 

plot_MIC
ggsave("Combined_plot_MIC.jpeg")
# Saving 16.6 x 11.7 in image
# Potentially make No genes always the same colour? 

# Not done #### 
plot1_MIC <- ggarrange(amoxicillin_DD,amox_clav_DD,cephalotin_DD,
                   labels = c("A", "B", "C"), 
                   nrow = 3, ncol = 1,align = "v") 
plot1_MIC
#plot1 + theme(axis.title = element_text(size=16), axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=16))
ggsave("Combined_plot_beta_lactams.jpeg",width = 7, height=11, units = "in")
# Original Saving 7.31 x 10.8 in image - make higher 

plot2_MIC <- ggarrange(florfenicol_DD,gentamicin_DD,tetracycline_DD,
                   trim_sulph_DD,
                   labels = c("A", "B", "C", "D"), 
                   nrow = 2, ncol = 2,align = "v") 

plot2_MIC
ggsave("Combined_plot_others.jpeg",scale =2) 
ggsave("Combined_plot_others2.jpeg",scale =1.5) 
ggsave("Combined_plot_others3.jpeg",scale = 1.1) 
# Saving 13.6 x 7.92 in image
# Gentamicin, Tetracycline, Colistin, Florfenicol, Trim_sulf

plot3_MIC <- ggarrange(cipro_DD,nalix_DD,marbo_DD,
                   labels = c("A", "B", "C"), 
                   nrow = 3, ncol = 1,align = "v",common.legend = TRUE,
                   legend = "right") 
plot3_MIC
#plot1 + theme(axis.title = element_text(size=16), axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=16))
ggsave("Combined_plot_quinolones_MIC.jpeg",width = 9, height=11, units = "in")
# Can do a white plot bakcground instead of the grey 
# Also get the plot area to be a bit higher so the numbers are visible 