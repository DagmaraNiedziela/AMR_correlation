library(tidyverse)
library(readxl)
#detach(package:plyr)    
#library(dplyr)

# Manipulate gene file ####
AMR_data2 <- read_excel("AMRFinder_all_results_short.xlsx", sheet = 1)
AMR_data2

# Clean up the data to get only the columns I need - cfs_name, gene name, subclass, percent coverage and identity
AMR_data2 <- AMR_data2 %>% select(cfs_name, gene_symbol, percent_coverage_of_reference_sequence,percent_identity_to_reference_sequence,class)
View(AMR_data2)
colnames(AMR_data2) <- c("cfs_name", "gene_symbol", "percent_coverage","percent_identity","resistance")

typeof(AMR_data2$resistance)
AMR_data2$resistance <- factor(AMR_data2$resistance)
levels(AMR_data2$resistance) 
View(AMR_data2)

[1] "AMINOGLYCOSIDE"      "BETA-LACTAM"         "COLISTIN"            "EFFLUX"             
[5] "FOSFOMYCIN"          "FOSMIDOMYCIN"        "MACROLIDE"           "MULTIDRUG"          
[9] "NITROFURAN"          "PHENICOL"            "QUINOLONE"           "QUINOLONE/TRICLOSAN"
[13] "STREPTOTHRICIN"      "SULFONAMIDE"         "TETRACYCLINE"        "TRIMETHOPRIM"  

AMR_data$resistance <- gsub("QUINOLONE/TRICLOSAN","Quinolone",AMR_data$resistance)

# Heatmap sorted alphabetically 
ggplot(AMR_data2, aes(x = gene_symbol, y = cfs_name, fill= resistance)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black')#+
#  theme(panel.background=element_rect(fill="white", colour="black"))

# In Chiara's code there are some aesthetics I could use 

ggsave("heatmap_sorted_by_name2.jpeg", scale = 2)
?geom_tile

# Try to add NA columns by pivoting it twice? ####
colnames(AMR_data2)
AMR_data2_long <- AMR_data2 %>% pivot_wider(id_cols =gene_symbol, names_from = cfs_name,values_from=percent_identity,values_fn = max)

AMR_data3 <- AMR_data2_long %>% pivot_longer(cols = -gene_symbol, names_to= "cfs_name", values_to= "percent_identity")
AMR_data3 <- left_join(AMR_data3, AMR_data2, by = c("cfs_name","gene_symbol"))

AMR_data3$cfs_name <- factor(AMR_data3$cfs_name)
levels(AMR_data3$cfs_name)

ggplot(AMR_data3, aes(x = gene_symbol, y = cfs_name, fill= resistance)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_discrete(na.value = 'white',na.translate=FALSE) +
  scale_color_discrete(na.value = 'black')+
  scale_y_discrete(limits = rev(levels(AMR_data3$cfs_name))) +
  theme(panel.background=element_rect(fill="white", colour="black"))

?scale_fill_manual

# In Chiara's code there are some aesthetics I could use 

ggsave("heatmap_sorted_by_name_NAs.jpeg", scale = 2) 

# Could also do one where height of the image is larger 
# Saving 16 x 14.4 in image 

ggsave("heatmap_sorted_by_name_NAs_2.jpeg", width = 14.4, height = 18, units = "in") 
?geom_tile

# Add phylogroup ####
# Heatmap with phylogroup as annotation of cfs numbers? 

phylogroup <- read_excel("phylogroup2.xlsx", sheet =1)
View(phylogroup)

AMR_data4 <- full_join(AMR_data3,phylogroup, by = "cfs_name")

AMR_data4 <- AMR_data4 %>% arrange(phylogroup,cfs_name)
typeof(AMR_data4$cfs_name)
AMR_data4$cfs_name <- factor(AMR_data4$cfs_name)
AMR_data4$phylogroup <- factor(AMR_data4$phylogroup)
levels(AMR_data4$phylogroup)
# "A"  "B1" "B2" "C"  "D"  "E"  "G" 
#levels(AMR_data4$phylogroup) <- c("A","B1","B2","D")
AMR_data4$cfs_name <- reorder(AMR_data4$cfs_name, AMR_data4$phylogroup)
levels(AMR_data4$cfs_name)
#levels(AMR_data2$resistance) 

library(forcats)
AMR_data4$phylogroup2 <- as.numeric(AMR_data4$phylogroup)
View(AMR_data4)

AMR_data4 <- AMR_data4 %>%
  mutate(cfs_name2 = fct_reorder(cfs_name, phylogroup2)) 
levels(AMR_data4$cfs_name2)

B <- ggplot(AMR_data4, aes(x = gene_symbol, y = cfs_name2, fill= resistance)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_discrete(na.value = 'white',na.translate=FALSE) +
  scale_color_discrete(na.value = 'black')+
  scale_y_discrete(limits = rev(levels(AMR_data4$cfs_name2))) +
  theme(panel.background=element_rect(fill="white", colour="black")) 

B

ggsave("heatmap_sorted_by_phylogroup_NAs.jpeg")
ggsave("heatmap_sorted_by_phylogroup_NAs_2.jpeg", width = 14.4, height = 18, units = "in") 

ggplot(AMR_data4, aes(x = 1, y = cfs_name2, fill= phylogroup)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_blank())+
  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black') +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  scale_y_discrete(limits = rev(levels(AMR_data4$cfs_name2)))

ggsave("phylogroup.jpeg", width = 3.34, height = 15, units = "in") 
# Saving 3.34 x 11.7 in image

A <- ggplot(AMR_data4, aes(x = 1, y = cfs_name2, fill= phylogroup)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, colour = "white")) +
  theme(axis.text.y = element_blank())+
  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black') +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  scale_y_discrete(limits = rev(levels(AMR_data4$cfs_name2)))
A

# Arrange with phylogroup ####
library(ggpubr)

ggarrange(
          ggarrange(A,NULL, nrow = 2, heights = c(1,0.05)),
          B,ncol=2,widths = c(0.15,1))
?ggarrange
# It looks crappy on the screen but good when saved 
# If there are any issues on another computer, play with heights and widths 

ggsave("heatmap_sorted_by_phylogroup_AB_newgroups.jpeg", width = 18.1, height = 15, units = "in", dpi = 1200) 
ggsave("heatmap_sorted_by_phylogroup_AB_newgroups_lower_res.jpeg", width = 18.1, height = 15, units = "in") 
# Saving 18.1 x 14.4 in image 
# Saving 14.5 x 11.7 in image

