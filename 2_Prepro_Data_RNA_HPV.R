# Script: Preprocess TCGA-CESC RNA-Seq data adapted from https://github.com/CSB-IG/SGCCA/blob/main/prepro-mRNA.R
#                                       and adapted from https://github.com/CSB-IG/SGCCA/blob/main/prepro-mRNA.R
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
# Author: Joel Ruiz Hernandez
# Contact: jjrhernandez@gmail.com
# Description: This script preprocesses RNA-Seq data from the TCGA-CESC project. Only tumor samples HPV high risk positive
#              Steps include quality control (QC), bias detection, correction, and normalization
#Input: Metadata and unstranded counts
#Output: QC_report graphs and Final counts_df, metadata and annotated matrix

## Set working directory --- ---
setwd("~/2_Prepro_Data_TCGA/A7_A9/")

#load(file = "~/2_Prepro_Data_TCGA/A7_A9/1_Prepro_Data_RNA_HPV.RData")

## Packages required --- ---
library(tidyselect)
library(biomaRt)
library(NOISeq)
library(edgeR)
library(EDASeq)
library(ggplot2)
library(vroom)
library(tidyverse)
#Seed for reproducibility
set.seed(123)

### 1) Input data --- ---
  #Load data
metadata <- vroom(file = "~/1_Get_Data_TCGA/1_2_Metadata.tsv") #[1] 278   8
unstranded_counts <- vroom(file = "~/1_Get_Data_TCGA/1_1_unstranded_counts.tsv") #[1] 60660   279 (1 col add for gene names)

  # Define factors correctly: Tumor samples (A7 and A9)
factors <- metadata %>% 
  dplyr::select(specimenID,HPV_clade, sample_type, HPV_type) %>% 
  dplyr::filter(!HPV_clade=="otro", !HPV_type=="HPV70", !sample_type=="Solid Tissue Normal") %>% #HPV 70 It's A7 but not high risk
  as.data.frame()   #[1] 268   4  (268 tumor samples with HPV A7 and A9 high risk)

  #Count matrix (samples only included in factors)
counts <- as.data.frame(unstranded_counts) %>% 
  dplyr::select(gene, factors$specimenID)#[1] 60660   269 (col -1, rownames)

  #Clean data
#keep only transcript id not version numbers
counts <- counts %>% mutate(gene = str_remove(gene, "\\..*$")) 
#Delete cero counts
counts <- counts[rowSums(counts[,-1]) != 0, ] 
#%>%  column_to_rownames(var = "gene") #[1] 56385   269


  #Generate mart object
    #Annnotate GC content, length & biotype per transcript
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
    # Retrieve annotations for genes in your counts matrix
myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position", "end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values = counts$gene,  # annotate the genes in the count matrix 
                 mart = mart)
    # Rename the column for clarity
myannot <- myannot %>% dplyr::rename(feature = ensembl_gene_id)  #[1] 56081     7
  # Add length column
myannot$length <- abs(myannot$end_position - myannot$start_position) #[1] 56081     8
    #filter transcripts without annotation
myannot <-  myannot[myannot$hgnc_id!="" & myannot$hgnc_symbol != "", ] #[1]  38610     8
    #Duplicates
myannot <-  myannot[!duplicated(myannot$feature), ]  #[1] 38607   8
    #synchronization
counts <- counts[counts$gene%in%myannot$feature,] #[1] 38607   269
matched_samples <- factors$specimenID[factors$specimenID %in% colnames(counts)]
factors <- factors %>% filter(specimenID %in% matched_samples)
counts <- counts %>% dplyr::select(gene, all_of(matched_samples))

all(counts$gene == myannot$feature) #TRUE
all(colnames(counts[-1])== factors$specimenID) #TRUE
  
  #Remove duplicate genes
  #Check if there are duplicated genes
  #Delete duplicates function (if were necessary)
  #del_dupl <- function(counts) {
  # Identify duplicated genes
  #  repeated_values <- counts %>%
  #    group_by(gene) %>%
  #    filter(n() > 1) %>%
  #    distinct(gene) %>%
  #    pull(gene)
  # Get the duplicated rows
  #  repeated_rows <- counts[counts$gene %in% repeated_values, ]
  # Sort and calculate the median for duplicated values
  #  repeated_rows <- repeated_rows[order(repeated_rows$gene), ]
  #  repeated_rows <- repeated_rows %>%
  #    group_by(gene) %>%
  #    summarize(across(everything(), median, na.rm = TRUE))
  # Remove duplicates and add the rows with the calculated median
  #  counts <- counts %>% filter(!gene %in% repeated_rows$gene)
  #  counts <- bind_rows(counts, repeated_rows)
  # Convert selected columns to integers
  #  counts <- counts %>% mutate(across(-gene, as.integer))
  #  return(counts)
  #}
  # Apply the function to delete duplicate genes
  #counts <- del_dupl(counts)

  #Rows in factors must be equal to number of columns in data.
rownames(counts) <- NULL
counts <- counts %>%
  column_to_rownames(var = "gene") #only once
nrow(factors) == ncol(counts) #TRUE
all(factors$specimenID == colnames(counts)) #TRUE

  #duplicated_symbols <- myannot$hgnc_symbol[duplicated(myannot$hgnc_symbol) | duplicated(myannot$hgnc_symbol, fromLast = TRUE)]
  #print(duplicated_symbols)

### 2)Converting data into a NOISeq object--- ---

  # Create named vectors for annotations
mylength <- setNames(myannot$length, myannot$feature)
mygc <- setNames(myannot$percentage_gene_gc_content, myannot$feature)
mybiotype <- setNames(myannot$gene_biotype, myannot$feature)

  #Create the NOISeq object
noiseqData_beforeNormal <- NOISeq::readData(data = counts,
                                            factors = factors,           # Variables indicating the experimental group for each sample
                                            gc = mygc,                   # %GC in myannot
                                            biotype = mybiotype,         # Biotype
                                            length = mylength)          # Gene length

### 3)Quality control of count data--- ---

  # 3.1)Biotype detection

    # 3.1.1) Biodetection plot: 
    #Detecting an enrichment in the sample of any other biotype
mybiodetection_HPV <- dat(noiseqData_beforeNormal, type = "biodetection", factor = "HPV_clade" )
png("1_Biodetection_HPV_persample.png")
par(mfrow = c(1, 2))
explo.plot(mybiodetection_HPV,
           plottype = "persample", #type of plot
           samples = c(1,2))
dev.off()  

png("1_Biodetection_HPV_comparison.png", width = 900, height = 350) # Ajusta el width
par(mfrow = c(1, 2))
explo.plot(mybiodetection_HPV,
           plottype = "comparison",
           samples = c(1,2))
dev.off()
  
  #[1] "Percentage of protein_coding biotype in each sample:"
  #A7      A9 
  #51.6564 49.9909 
  #[1] "Confidence interval at 95% for the difference of percentages: A7 - A9"
  #[1] 0.9578 2.3732
  #[1] "The percentage of this biotype is significantly DIFFERENT for these two samples (p-value = 3.811e-06 )."
  
  # 3.2) Sequencing depth & Expression Quantification

    # 3.2.1) Count distribution per sample  
mycountsbio_HPV = dat(noiseqData_beforeNormal, type = "countsbio", factor = "HPV_clade")
png("2_Counts_distribution_HPV_clade.png",width = 500, height = 520)
explo.plot(mycountsbio_HPV,
           plottype = "boxplot", #type of plot
           samples = 1:2)
dev.off()

  # 3.2.2) Sensitivity ploT: check for low count genes
    #Low counts may introduce noise in the data that makes more difficult to extract the relevant information,

png("3_LowCounts_HPV.png", width=500, height=520)
par(cex.axis=16,    # Aumentar tamaño de los ejes
    cex.lab=16,     # Aumentar tamaño de las etiquetas
    cex.main=24,      # Aumentar tamaño del título
    mar=c(5, 5, 4, 2) + 4)  # Ajustar márgenes
explo.plot(mycountsbio_HPV,
           plottype = "barplot", # tipo de gráfico
           samples = 1:2) 
dev.off()

  # Plot global distribution of CPM 
png("4_lowCountThres.png", width=800, height=600)
hist_values <- hist(rowMeans(cpm(counts, log=TRUE)), 
                    ylab="Número de genes", 
                    xlab="Media de log CPM", 
                    col="lightblue", 
                    border="darkblue", 
                    main="Distribución global de CPM",
                    breaks=30, 
                    xlim=c(min(rowMeans(cpm(counts, log=TRUE))), max(rowMeans(cpm(counts, log=TRUE)))), 
                    ylim=c(0, 12000),  # Ajustar límite superior del eje Y a 12000
                    xaxs="i", 
                    yaxs="i", 
                    cex.axis=1.5,  # Aumentar tamaño de los ejes
                    cex.lab=1.5,   # Aumentar tamaño de las etiquetas
                    cex.main=2)    # Aumentar tamaño del título
abline(v=0, col="red", lwd=2)  # Línea vertical en 0 con mayor grosor
box()  # Agrega un marco al gráfico
dev.off()


  # 3.3) Sequencing bias detection
    # 3.3.1) Length bias
    #Relationship between the gene length and the expression values
    #If the model p-value is significant and R2 value is high (more than 70%), the
    #expression depends on the feature length and the curve shows the type of dependence.
mylengthbias_HPV = dat(noiseqData_beforeNormal, factor = "HPV_clade", type = "lengthbias")
png("5_Length_bias_HPV.png", width=1300, height=550)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mylengthbias_HPV, samples = c(1,2), toplot = "global")
dev.off()

  # 3.3.2) GC content bias
  # Relationship between the gene GC content and the expression values
  # If the model p-value is signifficant and R2 value is high (more than 70%), the expression will depend on
  #the feature GC content and the curve will show the type of dependence.
png("6_GC_HPV.png", width=1300, height=550)
myGCbias_HPV = dat(noiseqData_beforeNormal, factor = "HPV_clade", type = "GCbias")
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(myGCbias_HPV, samples = 1:2, toplot = "global")
dev.off()

  # 3.3.3) RNA composition
    #Diagnostic of data

    #each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
    #by computing M values=log2(counts = countsr). 
    #Confidence intervals (CI) for the M median is computed by bootstrapping.
    #If the median of M values for each comparison is not in the CI, the deviation
    # of the sample is significant, therefore, normalization is needed 
    #"cd" means "Cumulative Distribution.

mycd <- dat(noiseqData_beforeNormal, type = "cd", norm = FALSE) # Slow operation
    #[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."
    # Verify diagnostic
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])
    #FAILED PASSED
    #264      3 

png("7_Mvalues_mycd.png", width=1000, height=600)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mycd, samples = c(250:261))
dev.off()

  # 3.4) PCA
  #used to visualize if the experimental samples are clustered according to the experimental desig
myPCA_HPV = dat(noiseqData_beforeNormal, type = "PCA", norm = F,logtransf = F)
png("8_PCA_Ori_hpv.png")
explo.plot(myPCA_HPV, plottype= "scores", factor = "HPV_clade")
dev.off()

  #3.5) Quality Control report
QCreport(noiseqData_beforeNormal, samples = NULL, factor = "HPV_clade", norm = FALSE)

### 4) Solve biases: Normalization, Low-count fltering & Batch efect correction --- ---
  # 4.1)Normalization
  #correct expression values based on gen length and expression values
myTMM = NOISeq::tmm(assayData(noiseqData_beforeNormal)$exprs, long = 1000, lc = 1) #lc = 1 <- apply standar length

  #4.2) Low-count fltering

count_filtered <- filtered.data(counts, factor = "HPV_clade",
                                norm = FALSE, depth = NULL, method = 1, cpm = 1, p.adj = "fdr")
  #CPM=1 <- 11351 features are to be kept for differential expression analysis with filtering method 1
 
  #Compare distribution
png("9_lowCountThres_after_HPV.png", width=800, height=600)
hist_values <- hist(rowMeans(cpm(count_filtered, log=TRUE)), 
                    ylab="Número de genes", 
                    xlab="Media de log CPM", 
                    col="lightblue", 
                    border="darkblue", 
                    main="Distribución global de CPM",
                    breaks=30, 
                    xlim=c(min(rowMeans(cpm(counts, log=TRUE))), max(rowMeans(cpm(counts, log=TRUE)))), 
                    ylim=c(0, 12000),  # Ajustar límite superior del eje Y a 12000
                    xaxs="i", 
                    yaxs="i", 
                    cex.axis=1.5,  # Aumentar tamaño de los ejes
                    cex.lab=1.5,   # Aumentar tamaño de las etiquetas
                    cex.main=2)    # Aumentar tamaño del título
abline(v=0, col="red", lwd=2)  # Línea vertical en 0 con mayor grosor
box()  # Agrega un marco al gráfico
dev.off()

#Synchronization
myannot_after <- myannot[myannot$feature%in%rownames(count_filtered),]

#all names must match
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(count_filtered),
  featureData=data.frame(myannot_after,row.names=myannot_after$feature),
  phenoData=data.frame(factors,row.names=factors$specimenID))

#order for less bias
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 

fullfullTMM <-NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)

#norm.counts <- betweenLaneNormalization(normCounts(lFull),
#                                        which = "median", offset = FALSE)

noiseqData_after <-  NOISeq::readData(data = fullfullTMM, factors=factors)

#Diagnostic Test
mycd_after <- NOISeq::dat(noiseqData_after,type="cd",norm=TRUE)
table(mycd_after@dat$DiagnosticTest[,  "Diagnostic Test"])
#With CPM=1
#FAILED PASSED 
#9    258 

#4.3) Batch effect correction
ffTMMARSyn_HPV=ARSyNseq(noiseqData_after, factor = "HPV_clade", batch = F,
                           norm = "n",  logtransf = T)

myPCA_after_norm = dat(ffTMMARSyn_HPV, type = "PCA", norm = T,logtransf = T)
png("10_myPCA_HPV_norm.png")
explo.plot(myPCA_after_norm, plottype = "scores", factor = "HPV_clade")
dev.off()
#save
saveRDS(myPCA_after_norm, "myPCA_HPV_norm.rds")
print(explo.plot(myPCA_after_norm, plottype = "scores", factor = "HPV_clade"))

##Same PCA but make step by step (optional)
#Extract the expression matrix
norm_ARSyn <- exprs(ffTMMARSyn_HPV)# [1] 11245   271
# Build the transposed matrix for PCA
pca_norm_ARSyn <- norm_ARSyn %>% 
  t()  # Transpose the matrix
# Perform PCA using `prcomp`
pca_norm_ARSyn <- prcomp(pca_norm_ARSyn, retx = TRUE, center = TRUE, scale. = FALSE)
#Convert PCA results to a data frame
pca_norm_ARSyn_df <- pca_norm_ARSyn$x %>% as.data.frame()

#Calculate variance explained by each principal component (PC)
variance_table_norm_ARSyn <- data.frame(
  PC = 1:length(pca_norm_ARSyn$sdev),  # Principal component number
  Variance_Percentage = pca_norm_ARSyn$sdev^2 / sum(pca_norm_ARSyn$sdev^2) * 100,  # Percentage of variance explained
  cumulative_percentage = cumsum(pca_norm_ARSyn$sdev^2 / sum(pca_norm_ARSyn$sdev^2) * 100))  # Cumulative variance

#Elbow_plot
# Calculate the differences between consecutive variance percentages
variance_diff <- diff(variance_table_norm_ARSyn$Variance_Percentage)
# Calculate the second difference to find the "elbow" (largest change in the slope)
second_diff <- diff(variance_diff)
# Find the index of the maximum second difference (elbow point)
elbow_index <- which.max(second_diff) + 1  # Adding 1 because diff reduces the length by 1
# Print the number of components where the elbow occurs
cat("The elbow point is at PC", elbow_index, "with", 
    round(variance_table_norm_ARSyn$Variance_Percentage[elbow_index], 2), "% variance explained.\n")
#Highlight the elbow point in the elbow plot
elbow_plot <- ggplot(variance_table_norm_ARSyn, aes(x = PC, y = Variance_Percentage)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  geom_vline(xintercept = elbow_index, linetype = "dashed", color = "green") +
  labs(title = "Elbow Plot of PCA",
       x = "Principal Component (PC)",
       y = "Variance Percentage (%)") +
  theme_minimal()
# Save the plot as a PNG file
ggsave("11_elbow_plot_HPV.png", elbow_plot, width = 6, height = 4)

# Plot PC2 vs PC1 colored by sample type, without specimen names
PC2_PC1_norm_ARSyn_HPV <- pca_norm_ARSyn_df %>% 
  ggplot() +
  aes(x = PC2, y = PC1, colour = factors$HPV_clade) +
  geom_point() +
  labs(title = "PCA Scatterplot by sample type",
       subtitle = "PC2 vs PC1", 
       x = paste("PC2 (", sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[2]), "%)"),
       y = paste("PC1 (", sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[1]), "%)")) +
  theme_minimal() +
  theme(legend.title = element_blank())

print(PC2_PC1_norm_ARSyn_HPV)

# Save the plot as a PNG file
ggsave("12_PC1_PC2_norm_ARSyn_HPV.png", PC2_PC1_norm_ARSyn_HPV , width = 6, height = 4)

#FINALY QUALITY CHECK
#Synchronization
counts_after_norm <- exprs(ffTMMARSyn_HPV)
myannot_after_norm <- myannot_after[myannot_after$feature%in%rownames(counts_after_norm),]

#By sample
noiseqData_after_norm = NOISeq::readData(data = counts_after_norm , gc = myannot_after_norm[,1:2],
                                         biotype = myannot_after_norm[,c(1,3)],
                                         factors = factors,
                                         length=myannot_after_norm[,c(1,8)])

mycountsbio_after_norm = dat(noiseqData_after_norm, type = "countsbio", factor = "HPV_clade",
                             norm=T)

png("13_CountsFinal_HPV.png",width = 500, height = 520)
explo.plot(mycountsbio_after_norm, plottype = "boxplot",samples=1:2)
dev.off()

#GCcontent per sample
myGCcontent_after <- dat(noiseqData_after_norm, k = 0, type = "GCbias", 
                         factor = "HPV_clade",norm=T)
#Residual standard error: 413.8 on 47 degrees of freedom
#Multiple R-squared:  0.5393,	Adjusted R-squared:  0.4608 
#F-statistic: 6.876 on 8 and 47 DF,  p-value: 5.935e-06
png("14_GCbiasFinal_HPV.png",width=1300, height=550)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(myGCcontent_after, samples = NULL, toplot = "global")
dev.off() 

##lenBias per sample
mylenBias_sample <- dat(noiseqData_after_norm, k = 0, type = "lengthbias", 
                        factor = "HPV_clade",norm=T)
#Residual standard error: 364.9 on 47 degrees of freedom
#Multiple R-squared:  0.1274,	Adjusted R-squared:  -0.02117 
#F-statistic: 0.8575 on 8 and 47 DF,  p-value: 0.5583

png("15_lengthbiasFinal_HPV.png", width=1300, height=550)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mylenBias_sample, samples = c(1,2), toplot = "global")
dev.off()

# Sensitivity ploT: check for low count genes
#Low counts may introduce noise in the data that makes more dificult to extract the relevant information,
png("16_LowCounts_HPV_norm.png", width=500, height=520)
par(cex.axis=8,    # Aumentar tamaño de los ejes
    cex.lab=8,     # Aumentar tamaño de las etiquetas
    cex.main=12,      # Aumentar tamaño del título
    mar=c(5, 5, 4, 2) + 4)  # Ajustar márgenes
explo.plot(mycountsbio_after_norm,
           plottype = "barplot", # tipo de gráfico
           samples = NULL) 
dev.off()


png("15_Mvalues_mycd_after_norm.png", width=1000, height=600)
par(cex.axis=1.5,   
    cex.lab=1.5,     
    cex.main=2)
explo.plot(mycd_after, samples = c(250:261))
dev.off()

#Quality Control report of normalized data
pdf("QCreport_after_HPV.pdf")
QCreport(noiseqData_after, samples = NULL, factor = "HPV_clade", norm = T)
dev.off()

###Save data
  #Final counts as TSV
Final_A7_A9 <- counts_after_norm %>% 
  as.data.frame() #[1] 11351   268
vroom::vroom_write(Final_A7_A9, file = "~/2_Prepro_Data_TCGA/A7_A9/Final_A7_A9.tsv", delim = "\t")

  #Factors
Factors_A7_A9 <- factors %>% as.data.frame() #[1] 268   4
vroom::vroom_write(Factors_A7_A9, file = "~/2_Prepro_Data_TCGA/A7_A9/Factors_A7_A9.tsv", delim = "\t")

#Annotated matrix
Myannot_A7_A9 <- myannot_after_norm %>% as.data.frame() #[1] 11351     8
vroom::vroom_write(Myannot_A7_A9, file = "~/2_Prepro_Data_TCGA/A7_A9/Myannot_A7_A9.tsv", delim = "\t")

#For AracneMUlticore:
Final_A7_A9_NW <- Final_A7_A9 %>% mutate(gene = rownames(Final_A7_A9), .before = 1) #[1] 11351   269 (+1)
# Save as TSV
vroom::vroom_write(Final_A7_A9_NW, file = "~/2_Prepro_Data_TCGA/A7_A9/Final_A7_A9_NW.tsv", delim = "\t")

   #A7
Factors_A7 <- Factors_A7_A9 %>% filter(HPV_clade=="A7") #[1] 66  4
Final_A7_NW <- Final_A7_A9 %>% select(Factors_A7$specimenID) %>% 
  mutate(gene = rownames(Final_A7_A9), .before = 1) #[1] 11351    67 (+1)
  # Save as TSV
vroom::vroom_write(Final_A7_NW, file = "~/2_Prepro_Data_TCGA/A7_A9/Final_A7_NW.tsv", delim = "\t")

  #A9
Factors_A9 <- Factors_A7_A9 %>% filter(HPV_clade=="A9") #[1] 202   4
Final_A9_NW <- Final_A7_A9 %>% select(Factors_A9$specimenID) %>% 
  mutate(gene = rownames(Final_A7_A9), .before = 1) #[1] 11351   203 (+1)
  # Save as TSV
vroom::vroom_write(Final_A9_NW, file = "~/2_Prepro_Data_TCGA/A7_A9/Final_A9_NW.tsv", delim = "\t")


#Save_image
save.image(file = "~/2_Prepro_Data_TCGA/A7_A9/1_Prepro_Data_RNA_HPV.RData")
