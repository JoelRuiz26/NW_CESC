#Joel Ruiz
#24-Jun-2024
#Obtain HPVclade of TCGA-CESC project from literature (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7459114/)
  #Input HPV_genotypes.tsv from a supplementary PDF of doi: 10.1038/s41598-020-71300-7. PMID: 32868873; PMCID: PMC7459114.
  #Output: 3 tsv with information of HPV type, clade and their TCGAIDsample

###Load data and packages
setwd("/home/jjruiz/0_HPV_Distribution")

library(dplyr)
library(vroom)
library(viridis)
library(stringr)
library(ggplot2)
library(gg.gap)
library(plotrix)

HPV_IDsample <- vroom('0_HPV_genotypes.tsv') 
#`TCGA Case ID` Total number of reads (h…¹ Included in TCGA-CES…² `TCGA final genotypes` `TCGA HPV genotypes`
# 1 TCGA           Case                       ID                     other)                 study?              
#2 TCGA-BI-A0VR   185589496                  TRUE                   HPV16                  HPV16               
#3 TCGA-BI-A0VS   213513044                  TRUE                   HPV16                  HPV16               
#4 TCGA-BI-A20A   156561650                  TRUE                   HPV16                  HPV16               
#5 TCGA-C5-A0TN   116784232                  TRUE                   HPV16                  HPV16               

###Clean table (Select only HPV positive and remove multiple-infection)
TCGA_HPV <- HPV_IDsample %>% 
  select(1,6) %>% #Select only col TCGAID and HPV_type
  slice(2:305) %>% #Cause the first line was wrong, i delete this
  filter(HPV_estimated!="negative") %>%  #Select only the HPV_positive
  mutate(HPV_true = str_replace(HPV_estimated, "^HPV\\d+\\s+", "") %>% 
      str_trim()
      ) %>% #Create a col that exclude the error in conversion step from pdf
  #^exclude HPV\\d+ two digits \\ plus space(s) and after delete any space
  mutate(tipos_HPV = str_extract_all(HPV_true,"HPV\\d+") %>% #Create new column
           #that extract everything after "HPV\\d+ folowed by two digits
           lapply(paste, collapse=", ") %>% #and put in the other culumn the HPV
           #type. In the case of multiple HPV, it separates with ","
           as.character()
         ) %>% 
  arrange(HPV_true) #ArrangeHPV number
#Seleccionar solo los que contengan un solo tipo de HPV y cortar la tabla
TCGA_HPV_mono <- TCGA_HPV %>% 
  filter(!str_detect(tipos_HPV, ",")) %>% 
  select(1,4)
#head(TCGA_HPV_mono) #275   2
#  `TCGA Case ID` tipos_HPV   
#1 TCGA-HM-A4S6   HPV16    
#2 TCGA-ZJ-AAXD   HPV16    
#3 TCGA-EA-A43B   HPV16 

#Distribution of HPV
HPV_Distrib <- TCGA_HPV_mono %>% 
  group_by(tipos_HPV) %>% 
  summarise(frec = n()) %>%
  mutate(percent = frec / nrow(TCGA_HPV_mono) * 100)%>%
  arrange(desc(frec))
#head(HPV_Distrib) #15 different kind of HPV   
#1 HPV16            166       60.4 
#2 HPV18             38       13.8 
#3 HPV45             20       7.27
#4 HPV33              8       2.91
#5 HPV52              8       2.91
#6 HPV31              7       2.55

# First 2 digits in %
HPV_Distrib$percent <- sprintf("%.2f", HPV_Distrib$percent)
write.table(x=HPV_Distrib, file = "1_1_HPV_Distrib.tsv",append=FALSE, quote=FALSE,
            sep="\t", row.names = FALSE, col.names = TRUE)

#Anexar los clados filogenéticos: Referencia: Rader JS, Tsaih SW, Fullin D,
  #Murray MW, Iden M, Zimmermann MT, Flister MJ. Genetic variations in human 
  #papillomavirus and cervical cancer outcomes. Int J Cancer.
  #2019 May 1;144(9):2206-2214. doi: 10.1002/ijc.32038. Epub 2019 Jan 4. 
  #PMID: 30515767; PMCID: PMC6450540.
  #LINk https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6450540/
  
TCGA_HPV_Clades <- TCGA_HPV_mono %>% 
  mutate(
    Clado_filogenetico = case_when(
      tipos_HPV %in% c("HPV16", "HPV31", "HPV33", "HPV35", "HPV52", "HPV58", "HPV67") ~ "A9",
      tipos_HPV %in% c("HPV18", "HPV39", "HPV45", "HPV59", "HPV68", "HPV70", "HPV85", "HPV97") ~ "A7",
      TRUE ~ "otro"
    )
  )
#head(TCGA_HPV_Clades)
#`TCGA Case ID` tipos_HPV Clado_filogenetico
#1 TCGA-HM-A4S6   HPV16     A9                
#2 TCGA-ZJ-AAXD   HPV16     A9                
#3 TCGA-EA-A43B   HPV16     A9                
#4 TCGA-MA-AA3W   HPV16     A9                
#5 TCGA-VS-A9UQ   HPV16     A9   

write.table(x=TCGA_HPV_Clades, file = "1_2_HPV_Clades.tsv",append=FALSE, quote=FALSE,
            sep="\t", row.names = FALSE, col.names = TRUE)

#Clades distribution
samples_caldes <- nrow(TCGA_HPV_Clades)
# Crear el dataframe con la frecuencia y porcentaje
HPV_Clades_dist <- TCGA_HPV_Clades %>% 
  group_by(Clado_filogenetico) %>% 
  summarise(frec = n()) %>%
  mutate(percent = frec / samples_caldes * 100)%>%
  arrange(desc(frec))

# FIrst 2 digit in percent
HPV_Clades_dist$percent <- sprintf("%.2f", HPV_Clades_dist$percent)

write.table(x=HPV_Clades_dist, file = "1_3_HPV_Clade_Distrib.tsv",append=FALSE, quote=FALSE,
            sep="\t", row.names = FALSE, col.names = TRUE)


###Graph HPV distribution in the samples
Graph_distribution <- ggplot(HPV_Distrib, aes(x = reorder(tipos_HPV, -frec), y = frec, fill = as.factor(-frec))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = viridis_pal(option = "E")(length(unique(HPV_Distrib$tipos_HPV)))) +  # Usar colores viridis
  theme_minimal() +
  labs(title = "Distribución de los tipos de HPV en cáncer cervicouterino",
       subtitle = "275 tumores de pacientes del TCGA", 
       x = "",  # Eliminar nombre del eje x
       y = "Número de pacientes") +
  theme(axis.line = element_line(color = "black"),  
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
        axis.text.y = element_text(color = "black"),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = NA),  
        legend.position = "none",  
        plot.title = element_text(margin = margin(t = 10, b = 5), size = 20, color = "black"),  
        plot.subtitle = element_text(margin = margin(b = 25), size = 15, color = "grey20"),  
        plot.margin = margin(t = 20, r = 20, b = 10, l = 20)) +  # Ajustar los márgenes
  geom_text(aes(label = paste("n =", frec, "\n(", percent, "%)")), vjust = -0.5, color = "grey40", size = 3.5, hjust = 0.5) +  
  scale_y_continuous(expand = c(0, 0))  # Eliminar espacio en el límite inferior

#VIew
print(Graph_distribution)

# SAve graph in png and PDF
ggsave("1_4Graph_HPV_Distrib_Bar.png", plot = Graph_distribution, width = 20, height = 12, units = "in", dpi = 300)
ggsave("1_4Graph_HPV_Distrib_Bar.pdf", plot = Graph_distribution, width = 20, height = 12, units = "in", dpi = 300)

