# Beijing Hospital Cohort Processing
rt <- read.csv("TCR_CDR3_length_merged.csv")
clin <- readxl::read_excel("Timepoint_Distinguished_Basic_Information.xlsx")

# Merge and filter data
dat1 <- merge(rt, clin, by.x = "phenotype", by.y = "sampleID2")
dat1 <- subset(dat1, timepoint == "1")
dat1 <- dat1[, -5]

# Thoracic Hospital Cohort Processing
rt2 <- read.csv("72TCR_CDR3_length_merged.csv")
dat2 <- rbind(dat1, rt2)

# Output merged data
write.csv(dat2, "TCR.120.72_Length_Merged_nt.csv", row.names = FALSE)
clin <- readxl::read_excel("Integrated_Basic_Information.xlsx")

# Process frequency data
freq <- dat2[, -3]
freq$CDR3_Length_bp <- as.factor(freq$CDR3_Length_bp)

# Merge with clinical data and adjust frequency scale
freq <- merge(freq, clin, by.x = "phenotype", by.y = "sampleID2")
freq$Frequency <- freq$Frequency * 100

# Convert length to numeric for segmentation
freq$CDR3_Length_bp <- as.numeric(freq$CDR3_Length_bp)

# Function to categorize CDR3 lengths into segments
get_segment <- function(length) {
  if (length <= 18) {
    return("<=18")
  } else if (length >= 52) {
    return(">=52")
  } else {
    lower_bound <- floor((length - 1) / 3) * 3 + 1
    upper_bound <- lower_bound + 2
    return(paste(lower_bound, upper_bound, sep = "~"))
  }
}

# Add segment category column
freq$Segment <- sapply(freq$CDR3_Length_bp, get_segment)

# Calculate correlation between frequency and age for each segment
outTab <- data.frame()
Segment_aa <- unique(freq$Segment)

for (n in Segment_aa) {
  freq_aa <- subset(freq, Segment == n)
  corT <- cor.test(freq_aa$Frequency, freq_aa$Age, method = "spearman")
  outTab <- rbind(outTab, cbind(Segment = n, cor = corT$estimate, pvalue = corT$p.value))
}

outTab$cor <- as.numeric(outTab$cor)
outTab$pvalue <- as.numeric(outTab$pvalue)

# Create significance label for p-values
outTab$Pvalue <- ifelse(outTab$pvalue < 0.05, "<0.05", ">=0.05")

# Order segments for plotting
outTab$Segment <- factor(outTab$Segment,
                         levels = c("<=18", "19~21", "22~24", "25~27", "28~30", "31~33", "34~36",
                                    "37~39", "40~42", "43~45", "46~48", "49~51", ">=52"),
                         ordered = TRUE)

# Plot correlation by segment
library(ggplot2)
library(ggeasy)

ggplot(outTab, aes(x = Segment, y = cor, fill = Pvalue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "TRB CDR3nt length (bp)", y = "Correlation with Age") +
  scale_fill_manual(values = c("<0.05" = "#E64B35CC", ">=0.05" = "#8491B4CC")) +
  theme_minimal() +
  easy_rotate_x_labels(angle = 45, side = "right")

ggsave(filename = "TRB_Segment_Age_Correlation.pdf", width = 6, height = 3.5)

# Calculate mean CDR3 length for each sample
rt <- read.csv("TCR.120.72_Length_Merged_nt.csv")
samples <- unique(rt$phenotype)

outTab <- data.frame()
for (n in samples) {
  rt2 <- subset(rt, phenotype == n)
  rt2$lengh_sum <- rt2$CDR3_Length_bp * rt2$Number
  meanLength2 <- sum(rt2$lengh_sum) / sum(rt2$Number)
  outTab <- rbind(outTab, cbind(sample = n, meanLength = meanLength2))
}

outTab$meanLength <- as.numeric(outTab$meanLength)

# Merge with clinical data for analysis
clin <- readxl::read_excel("Re_Aged_Grouped_Basic_Information.xlsx")
dat <- merge(outTab, clin, by.x = "sample", by.y = "sampleID2")

# Plot mean CDR3 length vs age
library(dplyr)
library(ggplot2)

p <- dat %>%
  ggplot(aes(Age, meanLength)) +
  labs(x = "Age", y = "Mean CDR3nt length (bp)") +
  geom_point(aes(color = cohort), size = 1.5, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#8491B4FF", fill = "#cbc9e2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(ggpubr)
p <- p + stat_cor(method = "spearman") + ggtitle(label = "BCR")
ggsave(filename = "All_Samples_Mean_CDR3_Length.pdf", width = 6, height = 4)

# Boxplot by age group
p2 <- ggboxplot(dat, x = "Group2", y = "meanLength", color = "timepoint", 
                ylab = "Mean CDR3nt length",
                xlab = "",
                add = "jitter",
                legend.title = "timepoint",
                palette = "npg",
                width = 0.6,
                order = c("30y_Group", "45y_Group", "60y_Group", "75y_Group"))

p2 <- p2 + rotate_x_text(45)
p2 <- p2 + stat_compare_means(aes(group = timepoint),
                              method = "wilcox.test",
                              symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                              label = "p.signif")

ggsave("Boxplot_Age_Group.pdf", width = 6, height = 4)

# Compare by gender
p1 <- dat %>%
  ggplot(aes(Age, meanLength)) +
  labs(x = "Age", y = "Mean CDR3nt length") +
  geom_point(size = 1.5, alpha = 0.3, color = "#4DBBD5FF") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#8491B4FF", fill = "#cbc9e2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ Gender, ncol = 2)

p1 <- p1 + stat_cor(method = "spearman") + ggtitle(label = "")
ggsave("Gender_Comparison.pdf", width = 5.5, height = 3)

# Boxplot comparing genders within age groups
p2 <- ggboxplot(dat, x = "Group2", y = "meanLength", color = "Gender", 
                ylab = "Mean CDR3nt length",
                xlab = "",
                add = "jitter",
                legend.title = "Gender",
                palette = "npg",
                width = 0.6,
                order = c("30y_Group", "45y_Group", "60y_Group", "75y_Group"))

p2 <- p2 + rotate_x_text(45)
p2 <- p2 + stat_compare_means(aes(group = Gender),
                              method = "wilcox.test",
                              symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                              label = "p.signif")

ggsave(filename = "Gender_Comparison_by_Age_Group.pdf", width = 7, height = 4.5)

# Combine plots
library(patchwork)
(p + p1) / p2
ggsave("Combined_Analysis_Plot.png", width = 9, height = 8)

# Beijing Hospital Cohort Processing
rt <- read.csv("TCR_CDR3_length_aa_merged.csv")
clin <- readxl::read_excel("Timepoint_Distinguished_Basic_Information.xlsx")

# Merge and filter data
dat1 <- merge(rt, clin, by.x = "phenotype", by.y = "sampleID2")
dat1 <- subset(dat1, timepoint == "1")
dat1 <- dat1[, -5]

# Thoracic Hospital Cohort Processing
rt2 <- read.csv("72TCR_CDR3_length_aa_merged.csv")
dat2 <- rbind(dat1, rt2)
clin <- readxl::read_excel("Re_Aged_Grouped_Basic_Information.xlsx")

# Process frequency data
freq <- dat2[, -3]
freq$CDR3.length.aa <- as.factor(freq$CDR3.length.aa)

# Merge with clinical data
dat <- merge(freq, clin, by.x = "phenotype", by.y = "sampleID2")
unique(dat$Group)

# Split data by age group for correlation analysis
time1 <- subset(dat, Group == "Young")[, c(1, 2, 3)]
time2 <- subset(dat, Group == "Middle")[, c(1, 2, 3)]
time3 <- subset(dat, Group == "Aged")[, c(1, 2, 3)]
time4 <- subset(dat, Group == "Old")[, c(1, 2, 3)]

# Convert to wide format for correlation analysis
library(tidyr)
spread_tmp1 <- spread(time1, phenotype, Frequency)
spread_tmp2 <- spread(time2, phenotype, Frequency)
spread_tmp3 <- spread(time3, phenotype, Frequency)
spread_tmp4 <- spread(time4, phenotype, Frequency)

# Replace NA with 0
spread_tmp1[is.na(spread_tmp1)] <- 0
spread_tmp2[is.na(spread_tmp2)] <- 0
spread_tmp3[is.na(spread_tmp3)] <- 0
spread_tmp4[is.na(spread_tmp4)] <- 0

# Format matrices for correlation calculation
rownames(spread_tmp1) <- spread_tmp1[, 1]
spread_tmp1 <- spread_tmp1[, -1]

rownames(spread_tmp2) <- spread_tmp2[, 1]
spread_tmp2 <- spread_tmp2[, -1]

rownames(spread_tmp3) <- spread_tmp3[, 1]
spread_tmp3 <- spread_tmp3[, -1]

rownames(spread_tmp4) <- spread_tmp4[, 1]
spread_tmp4 <- spread_tmp4[, -1]

# Calculate correlation matrices
library(Hmisc)
res1 <- rcorr(as.matrix(spread_tmp1))  # Default Pearson correlation
res2 <- rcorr(as.matrix(spread_tmp2))
res3 <- rcorr(as.matrix(spread_tmp3))
res4 <- rcorr(as.matrix(spread_tmp4))

# Heatmap visualization
library(pheatmap)
col <- colorRampPalette(c("blue", "white", "red"))(20)

pheatmap(res1$r, scale = "none", show_rownames = FALSE, show_colnames = FALSE, color = col, border = FALSE)
pheatmap(res2$r, scale = "none", show_rownames = FALSE, show_colnames = FALSE, color = col, border = FALSE)

# Convert correlation matrices to long format
library(reshape2)
melted_data1 <- melt(res1$r)
melted_data2 <- melt(res2$r)
melted_data3 <- melt(res3$r)
melted_data4 <- melt(res4$r)

# Filter correlation pairs by threshold
melted_data1_select <- subset(melted_data1, value >= 0.8)
num1_1 <- length(melted_data1_select$Var1)

melted_data2_select <- subset(melted_data2, value >= 0.8)
num2_1 <- length(melted_data2_select$Var1)

melted_data3_select <- subset(melted_data3, value >= 0.8)
num3_1 <- length(melted_data3_select$Var1)

melted_data4_select <- subset(melted_data4, value >= 0.8)
num4_1 <- length(melted_data4_select$Var1)

# Filter anti-correlated pairs
melted_data1_select2 <- subset(melted_data1, value < 0.8)
num1_2 <- length(melted_data1_select2$Var1)

melted_data2_select2 <- subset(melted_data2, value < 0.8)
num2_2 <- length(melted_data2_select2$Var1)

melted_data3_select2 <- subset(melted_data3, value < 0.8)
num3_2 <- length(melted_data3_select2$Var1)

melted_data4_select2 <- subset(melted_data4, value < 0.8)
num4_2 <- length(melted_data4_select2$Var1)

# Prepare data for stacked bar plot
values <- c(num1_1, num1_2, num2_1, num2_2, num3_1, num3_2, num4_1, num4_2)
groups <- c("30y_Group", "30y_Group", "45y_Group", "45y_Group", "60y_Group", "60y_Group", "75y_Group", "75y_Group")
labels <- c("cor>=0.8", "cor<0.8", "cor>=0.8", "cor<0.8", "cor>=0.8", "cor<0.8", "cor>=0.8", "cor<0.8")

data_df <- data.frame(Group = groups, Value = values, Label = labels)

# Create stacked bar plot
library(ggplot2)
library(scales)
library(ggsci)

p <- ggplot(data = data_df, aes(x = Group, y = Value, fill = Label)) + 
  theme_bw() +
  geom_bar(stat = 'identity', position = 'fill', width = 0.5) +
  labs(x = '', y = 'Relative Abundance', fill = NULL) +
  scale_fill_manual(values = c("cor<0.8" = "#8491B4CC", "cor>=0.8" = "#F39B7FCC")) +
  scale_y_continuous(labels = percent) +
  guides(fill = guide_legend(ncol = 1, bycol = TRUE, override.aes = list(size = 5))) +
  theme(
    axis.title.y = element_text(face = 'bold', color = 'black', size = 14),
    axis.title.x = element_text(face = 'bold', color = 'black', size = 14, vjust = -1.2),
    axis.text.y = element_text(face = 'bold', color = 'black', size = 10),
    axis.text.x = element_text(face = 'bold', color = 'black', size = 12, angle = 45, vjust = 0.5),
    panel.grid = element_blank(),
    legend.position = 'right',
    legend.key.height = unit(0.6, 'cm'),
    legend.text = element_text(face = 'bold', color = 'black', size = 10)
  )

# Chi-square test for contingency table
library(dplyr)

contingency_table <- data_df %>%
  group_by(Group, Label) %>%
  summarise(Frequency = sum(Value)) %>%
  pivot_wider(names_from = Label, values_from = Frequency, values_fill = 0)

chisq_test_result <- chisq.test(contingency_table[, -1])
print(chisq_test_result)

p <- p + ggtitle(label = paste0("TRB p < 0.001"))
ggsave("TRB_Correlation_Barplot.pdf", width = 5, height = 3)
# Heatmap for All Samples
spread_tmp3 <- spread(rt[, -2], phenotype, Frequency)
spread_tmp3[is.na(spread_tmp3)] <- 0  # Replace NA with 0 (frequency = 0)
rownames(spread_tmp3) <- spread_tmp3[, 1]
spread_tmp3 <- spread_tmp3[, -1]

# Calculate correlation matrix
res3 <- rcorr(as.matrix(spread_tmp3))

# Basic heatmap visualization
pheatmap(res3$r, scale = "none", show_rownames = TRUE, show_colnames = FALSE, color = col)
heatmap(x = res3$r, col = col, symm = TRUE, scale = "none", xlab = NULL, ylab = NULL)  # symm = TRUE for symmetric matrix

# Create annotation file
ann <- data.frame(Group = clin[, 5])
rownames(ann) <- clin$sampleID2  # Match row names with sample IDs

# Enhanced heatmap with annotations
col <- colorRampPalette(c("blue", "white", "red"))(20)  # Custom color palette
pheatmap(res3$r, 
         scale = "none", 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         color = col,
         annotation_row = ann,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_colors = list(timepoint = c("1" = "green", "2" = "yellow")), 
         border = FALSE)

# ComplexHeatmap visualization for better clustering
library(ComplexHeatmap)
mat <- as.data.frame(res3$r)
mat$samples <- rownames(mat)

# Prepare annotation data
ann2 <- ann
ann2$sample <- rownames(ann2)
ann3 <- merge(ann2, mat, by.x = "sample", by.y = "samples")

# Subset matrix to paired samples only
mat <- mat[ann3$sample, ann3$sample]

# Define row annotations
ha <- rowAnnotation(
  Group = ann3$Group,
  col = list( 
    Group = c("30y_Group" = "#E64B35CC", "45y_Group" = "#4DBBD5CC",
              "60y_Group" = "#00A087CC", "75y_Group" = "#3C5488CC")
  )
)

# Define column annotations
ha2 <- HeatmapAnnotation(
  Group = ann3$Group,
  col = list( 
    Group = c("30y_Group" = "#E64B35CC", "45y_Group" = "#4DBBD5CC",
              "60y_Group" = "#00A087CC", "75y_Group" = "#3C5488CC")
  )
)

# Plot and save ComplexHeatmap
Heatmap(as.matrix(mat), 
        name = "Cor", 
        column_title = "TRB",
        show_row_names = FALSE,
        show_column_names = FALSE,
        right_annotation = ha,
        top_annotation = ha2)

pdf("TRB_ComplexHeatmap.pdf", width = 7, height = 5.5)
Heatmap(as.matrix(mat), 
        name = "Cor", 
        column_title = "TRB",
        show_row_names = FALSE,
        show_column_names = FALSE,
        right_annotation = ha,
        top_annotation = ha2)
dev.off()


# Amino Acid Frequency Analysis
rt <- read.csv("all_clonotype_statistics_merged.csv")
rt <- rt[, c(1, 5, 2, 3)]  # Subset relevant columns

# Calculate amino acid counts per clone and sample
samples <- unique(rt$phenotype)
outTab <- data.frame()

for (n in samples) {
  rt_sample <- subset(rt, phenotype == n)
  
  # Get unique amino acids from CDR3 sequences
  characters_to_count <- unique(unlist(strsplit(rt_sample$cdr3aa, "")))
  sequences <- rt_sample$cdr3aa
  
  # Initialize result dataframe
  result_df <- data.frame(matrix(nrow = length(sequences), ncol = length(characters_to_count)))
  colnames(result_df) <- characters_to_count
  result_df$Sequence <- sequences
  
  # Count amino acid occurrences in each sequence
  for (i in 1:length(sequences)) {
    sequence <- sequences[i]
    for (j in 1:length(characters_to_count)) {
      char <- characters_to_count[j]
      count <- nchar(gsub(paste0("[^", char, "]"), "", sequence))
      result_df[i, j] <- count
      result_df$sample <- n
    }
  }
  
  outTab <- rbind(outTab, result_df)
  print(which(n == samples))  # Progress tracker
}

save(outTab, file = "aa_count_per_clone.Rdata")
write.csv(outTab, "aa_count_per_clone.csv", row.names = FALSE)


# Calculate total amino acid frequencies per sample
load("aa_count_per_clone.Rdata")
clin <- readxl::read_excel("re_aged_grouped_basic_information.xlsx")
clin <- as.data.frame(clin)

sampless <- unique(outTab$sample)
outTab2 <- data.frame()  # Amino acid counts per clone
outTab3 <- data.frame()  # Total amino acid counts per sample

for (m in sampless) {
  rt_sample2 <- subset(rt, phenotype == m)
  outTab_sample <- subset(outTab, sample == m)
  
  # Merge with clone counts
  dat_new <- merge(outTab_sample, rt_sample2, by.x = "Sequence", by.y = "cdr3aa")
  dat_new <- unique(dat_new)
  
  # Calculate weighted amino acid counts (count * clone frequency)
  result <- dat_new %>%
    mutate(
      C_total = C * count,
      P_total = P * count,
      R_total = R * count,
      S_total = S * count,
      N_total = N * count,
      Y_total = Y * count,
      E_total = E * count,
      Q_total = Q * count,
      F_total = F * count,
      A_total = A * count,
      I_total = I * count,
      L_total = L * count,
      W_total = W * count,
      V_total = V * count,
      G_total = G * count,
      T_total = T * count,
      D_total = D * count,
      H_total = H * count,
      K_total = K * count,
      M_total = M * count
    ) %>%
    select(phenotype, Sequence, ends_with("_total"))
  
  outTab2 <- rbind(outTab2, result)
  
  # Sum counts per sample
  result2 <- colSums(result[, -c(1, 2)])
  result2 <- as.data.frame(t(result2))
  result2$sample <- m
  outTab3 <- rbind(outTab3, result2)
  
  print(which(m == sampless))  # Progress tracker
}

write.csv(outTab2, "per_clone_aa_counts.csv", row.names = FALSE)
write.csv(outTab3, "sample_aa_summary.csv", row.names = FALSE)


# Convert counts to percentages
rt2 <- read.csv("sample_aa_summary.csv", check.names = FALSE)
rownames(rt2) <- rt2$sample
rt2 <- rt2[, -which(colnames(rt2) %in% c("sample", "_total", "*_total"))]  # Remove non-amino acid columns

# Calculate percentages
percentage_data <- rt2 / rowSums(rt2) * 100
write.csv(percentage_data, "aa_percentage_results.csv", row.names = TRUE)

# Merge with clinical data
percentage_data$sample <- rownames(percentage_data)
dat <- merge(percentage_data, clin, by.x = "sample", by.y = "sampleID2")


# Correlation between age and amino acid frequencies
library(ggplot2)
library(ggpubr)

theme_set(ggpubr::theme_pubr() + theme(legend.position = "top"))
genes <- colnames(rt2)

# Generate correlation plots for each amino acid
for (i in genes) {
  dat_new <- data.frame(Age = dat$Age, Frequency = dat[, i])
  
  ggscatter(dat_new, x = "Age", y = "Frequency",
            add = "reg.line", conf.int = TRUE,    
            add.params = list(fill = "lightgray")) +
    stat_cor(method = "spearman") +
    geom_point(color = "#8DA0CB") +
    ggtitle(label = i) +
    ylab("Frequency (%)") +
    xlab("Age")
  
  ggsave(filename = paste0("AA_", i, "_age_correlation.pdf"), width = 5, height = 4)
}

# Summarize correlation results
outTab <- data.frame()
for (n in genes) {
  cor_T <- cor.test(as.numeric(dat[, n]), dat$Age, method = "spearman")
  outTab <- rbind(outTab, data.frame(aa = n, cor = cor_T$estimate, pvalue = cor_T$p.value))
}

# Clean up amino acid names
outTab$aa <- sapply(strsplit(outTab$aa, "_"), "[", 1)
write.csv(outTab, "TCR_AA_correlation_test.csv", row.names = FALSE)


# Visualize average amino acid frequencies
dat2_average <- colMeans(dat[, genes, drop = FALSE])
dat2_average <- as.data.frame(dat2_average)
dat2_average$aa <- rownames(dat2_average)
dat2_average$aa <- sapply(strsplit(dat2_average$aa, "_"), "[", 1)

# Categorize amino acids by polarity
library(dplyr)
dat2_average <- dat2_average %>%
  mutate(category = case_when(
    aa %in% c("I", "V", "L", "F", "C", "M", "A", "W") ~ "HYDROPHOBIC",
    aa %in% c("G", "T", "S", "Y", "P", "H") ~ "NEUTRAL",
    aa %in% c("N", "D", "Q", "E", "K", "R") ~ "HYDROPHILIC"
  ))

dat2_average$percent <- round(dat2_average$dat2_average, 1)

# Plot average frequencies
library(ggsci)
p5 <- ggdotchart(dat2_average, x = "aa", y = "percent",
                 color = "category",
                 dot.size = 6,
                 palette = pal_npg("nrc", alpha = 0.8)(3),
                 add = "segments",
                 sorting = "descending",
                 group = "category",
                 label = "percent",
                 font.label = list(color = "white", face = "bold", size = 8, vjust = 0.5),
                 xlab = "",
                 ylab = "Average frequency (%)",
                 x.text.col = TRUE,
                 add.params = list(color = "category", size = 1),
                 ggtheme = theme(axis.text.x = element_text(angle = 45, hjust = 1)),
                 title = "TRB",
                 rotate = TRUE)

ggpar(p5, legend.title = "Category")
ggsave(filename = "TRB_AA_average_frequency.pdf", width = 6, height = 5)


# V-Gene Usage Analysis
rt <- read.csv("all_TCR_V_gene_usage_merged.csv")
data <- rt[, -3]  # Remove unnecessary column

# Check total frequency per sample (should sum to 100)
bb <- subset(data, phenotype == "GKZL0220205130")
sum(bb$Frequency)

# Reshape data for correlation analysis
library(tidyr)
spread_tmp <- spread(data, phenotype, Frequency)
spread_tmp[is.na(spread_tmp)] <- 0  # Replace NA with 0
write.table(spread_tmp, file = "V_gene_matrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Prepare matrix for correlation
rownames(spread_tmp) <- spread_tmp[, 1]
spread_tmp <- spread_tmp[, -1]
genes <- rownames(spread_tmp)
spread_tmp <- as.data.frame(t(spread_tmp))

# Merge with clinical data
spread_tmp$sample <- rownames(spread_tmp)
clin <- readxl::read_excel("re_aged_grouped_basic_information.xlsx")
clin <- as.data.frame(clin)
dat <- merge(spread_tmp, clin, by.x = "sample", by.y = "sampleID2")

# Correlation between V-gene usage and age
outTab <- data.frame()
for (n in genes) {
  cor_T <- cor.test(as.numeric(dat[, n]), dat$Age, method = "spearman")
  outTab <- rbind(outTab, data.frame(gene = n, cor = cor_T$estimate, pvalue = cor_T$p.value))
}

outTab$cor <- as.numeric(outTab$cor)
outTab$pvalue <- as.numeric(outTab$pvalue)
write.csv(outTab, "TCR_V_gene_correlation_test.csv", row.names = FALSE)

# Extract significant V-genes
outTab_sig <- subset(outTab, pvalue < 0.05)
gene_sig <- outTab_sig$gene
gene_sig_matrix <- spread_tmp[, gene_sig]
write.csv(gene_sig_matrix, "TCR_V_gene_significant_matrix.csv", row.names = FALSE)


# Visualize V-gene correlations
outTab <- na.omit(outTab)
outTab$p_significance <- outTab$pvalue < 0.05  # Flag significant correlations

# Horizontal bar plot
ggplot(outTab, aes(x = cor, y = reorder(gene, -cor), fill = p_significance)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("FALSE" = "#8491B4CC", "TRUE" = "#E64B35CC")) +
  ylab("TRB V gene") +
  xlab("Correlation with Age") +
  labs(title = "") +
  theme_light() +
  guides(fill = guide_legend(title = "P<0.05")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))

ggsave(filename = "TCR_V_gene_correlation_plot.pdf", width = 5, height = 11)

# Generate individual correlation plots for each V-gene
for (i in genes) {
  dat_new <- data.frame(Age = dat$Age, Frequency = dat[, i])
  
  ggscatter(dat_new, x = "Age", y = "Frequency",
            add = "reg.line", conf.int = TRUE,    
            add.params = list(fill = "lightgray")) +
    stat_cor(method = "spearman") +
    geom_point(color = "#8DA0CB") +
    ggtitle(label = i) +
    ylab("Frequency (%)") +
    xlab("Age")
  
  ggsave(filename = paste0("V_gene_", i, "_age_correlation.pdf"), width = 5, height = 4)
}


# Advanced Immun repertoire Analysis with immunarch
library(immunarch)
library(ggsci)

# Load data (10x format)
clin <- readxl::read_excel("re_aged_grouped_basic_information.xlsx")
immdata_10x <- repLoad("TCR_clone_files", .mode = "single")  # Adjust path as needed

# Basic repertoire statistics
repExplore(immdata_10x$data, "clones") %>% vis()  # Number of clones per sample
cloneNum <- as.data.frame(repExplore(immdata_10x$data, "clones"))

# Correlation between clone count and age
dat <- merge(cloneNum, clin, by.x = "Sample", by.y = "sampleID2")
ggscatter(dat, x = "Age", y = "Clones",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray")) +
  stat_cor(method = "spearman") +
  geom_point(color = "#8DA0CB") +
  ylab("Clones") +
  xlab("Age")

ggsave(filename = "TRB_clone_count_age_correlation.pdf", width = 5, height = 4)


# CDR3 Length Distribution Analysis
# Nucleotide length distribution
repExplore(immdata_10x$data, "lens", "nt") %>% vis()
ggsave(filename = "CDR3nt_length_distribution.pdf", width = 20, height = 7)

repExplore(immdata_10x$data, "lens", "nt") %>% vis(.by = "Group", .meta = immdata_10x$meta) + scale_fill_npg()
ggsave(filename = "CDR3nt_length_by_group.pdf", width = 10, height = 7)

# Amino acid length distribution
repExplore(immdata_10x$data, "lens", "aa") %>% vis()
ggsave(filename = "CDR3aa_length_distribution.pdf", width = 20, height = 7)

repExplore(immdata_10x$data, "lens", "aa") %>% vis(.by = "Group", .meta = immdata_10x$meta) + scale_fill_npg()
ggsave(filename = "CDR3aa_length_by_group.pdf", width = 10, height = 7)


# Clonality Analysis
# Homeostasis of clonotype abundance
repClonality(immdata_10x$data, "homeo") %>% vis() + scale_fill_npg()
ggsave(filename = "clonality_homeostasis.pdf", width = 20, height = 7)

repClonality(immdata_10x$data, "homeo") %>% vis(.by = "Group", .meta = immdata_10x$meta) + scale_fill_npg()
ggsave(filename = "clonality_homeostasis_by_group.pdf", width = 10, height = 7)

# Clonal proportion analysis
repClonality(immdata_10x$data, "clonal.prop") %>% vis()
ggsave(filename = "clonal_proportion.pdf", width = 20, height = 7)

repClonality(immdata_10x$data, "clonal.prop") %>% vis(.by = "Group", .meta = immdata_10x$meta) + scale_fill_npg()
ggsave(filename = "clonal_proportion_by_group.pdf", width = 6, height = 5)

# Top clonotype analysis
repClonality(immdata_10x$data, "top", .head = c(10, 50, 100, 500, 1000, 1e+05)) %>% vis() + scale_fill_npg()
ggsave(filename = "top_clonotypes.pdf", width = 20, height = 7)

repClonality(immdata_10x$data, "top") %>% vis(.by = "Group", .meta = immdata_10x$meta) + scale_fill_npg()
ggsave(filename = "top_clonotypes_by_group.pdf", width = 6, height = 5)


# Repertoire Overlap Analysis
# Public clonotypes between samples
repOverlap(immdata_10x$data, "public") %>% vis()
ggsave(filename = "public_clonotypes_heatmap.pdf", width = 25, height = 25)

# Jaccard index for repertoire similarity
repOverlap(immdata_10x$data, "jaccard") %>% vis()
ggsave(filename = "jaccard_similarity.pdf", width = 25, height = 25)

# MDS visualization of overlap
imm_ov1 <- repOverlap(immdata_10x$data, .method = "public", .verbose = FALSE)
repOverlapAnalysis(imm_ov1, "mds") %>% vis()
ggsave(filename = "repertoire_overlap_mds.pdf", width = 8, height = 6)

# Public clonotype table (amino acid + V gene)
pr.aav <- pubRep(immdata_10x$data, "aa+v", .verbose = TRUE)
write.csv(pr.aav, "public_clonotypes_aa_v.csv", row.names = FALSE)

# Group-specific clonotypes
pr1 <- pubRepFilter(pr.aav, immdata_10x$meta, c(Group = "30y_Group"))
pr2 <- pubRepFilter(pr.aav, immdata_10x$meta, c(Group = "45y_Group"))
pr3 <- pubRepFilter(pr.aav, immdata_10x$meta, c(Group = "60y_Group"))
pr4 <- pubRepFilter(pr.aav, immdata_10x$meta, c(Group = "75y_Group"))

# Shared clonotypes across all groups
seq.list <- list(pr1$CDR3.aa, pr2$CDR3.aa, pr3$CDR3.aa, pr4$CDR3.aa)
all_intersections <- Reduce(intersect, seq.list)
write.csv(data.frame(shared_CDR3aa = all_intersections), "shared_clonotypes_all_groups.csv", row.names = FALSE)


# Gene Usage Analysis
# V-gene usage (IGHV)
geneUsage(immdata_10x$data, "hs.ighv", .type = "segment") %>% vis("hist", .by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "IGHV_gene_usage_by_group.pdf", width = 20, height = 7)

# J-gene usage (IGHJ)
geneUsage(immdata_10x$data, "hs.ighj", .type = "segment") %>% vis("hist", .by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "IGHJ_gene_usage_by_group.pdf", width = 20, height = 7)

# D-gene usage (IGHD)
geneUsage(immdata_10x$data, "hs.ighd", .type = "segment") %>% vis("hist", .by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "IGHD_gene_usage_by_group.pdf", width = 20, height = 7)

# Gene usage divergence (Jensen-Shannon)
imm_gu <- geneUsage(immdata_10x$data, "hs.ighv", .norm = TRUE)
imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = FALSE)
p1 <- vis(imm_gu_js, .title = "IGHV JS Divergence", .leg.title = "JS")
ggsave(filename = "IGHV_JS_divergence.pdf", width = 10, height = 8)


# Spectratyping (length distribution by gene)
p1 <- vis(spectratype(immdata_10x$data[[1]], .quant = "id", .col = "aa"))
p2 <- vis(spectratype(immdata_10x$data[[1]], .quant = "id", .col = "aa+v"))
(p1 + p2)
ggsave(filename = "spectratype_aa_and_aa_v.pdf", width = 12, height = 6)


# Diversity Analysis
# Chao1 diversity index
repDiversity(immdata_10x$data) %>% vis(.by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "chao1_diversity_by_group.pdf", width = 6, height = 5)

# D50 diversity index
repDiversity(immdata_10x$data, .method = "d50") %>% vis(.by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "d50_diversity_by_group.pdf", width = 6, height = 5)

# Diversity by sequence feature (nucleotide/amino acid/V gene)
repDiversity(immdata_10x$data, .method = "div", .col = "nt") %>% vis(.by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "nucleotide_diversity_by_group.pdf", width = 6, height = 5)

repDiversity(immdata_10x$data, .method = "div", .col = "aa+v") %>% vis(.by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "aa_v_diversity_by_group.pdf", width = 6, height = 5)


# Clonotype Tracking
# Track top 5 clonotypes from first sample
tc1 <- trackClonotypes(immdata_10x$data, list(1, 5), .col = "nt")
p1 <- vis(tc1, .by = "Group", .meta = immdata_10x$meta)
ggsave(filename = "top_clonotype_tracking.pdf", width = 13, height = 9)

# Track specific clonotypes by amino acid sequence
target <- c("CASSLEETQYF", "CASSDSSGGANEQFF", "CASSDSSGSTDTQYF")
tc <- trackClonotypes(immdata_10x$data, target, .col = "aa")
p1 <- vis(tc, .plot = "smooth")
p2 <- vis(tc, .plot = "area")
p3 <- vis(tc, .plot = "line")
(p1 + p2 + p3)
ggsave(filename = "specific_clonotype_tracking.pdf", width = 25, height = 9)


# Kmer and Motif Analysis
# Calculate 5-mers from first sample
kmers <- getKmers(immdata_10x$data[[1]], 5)
write.csv(kmers, "5mers_first_sample.csv", row.names = FALSE)

# Visualize top kmers
p1 <- vis(kmers, .head = 10, .position = "stack")
p2 <- vis(kmers, .head = 10, .position = "dodge", .log = TRUE)
(p1 + p2)
ggsave(filename = "kmer_visualization.pdf", width = 15, height = 6)

# Motif analysis (sequence logo)
kp <- kmer_profile(kmers, "self")
p1 <- vis(kp)  # Text logo
p2 <- vis(kp, .plot = "seq")  # Sequence logo
(p1 + p2)
ggsave(filename = "sequence_motif_logos.pdf", width = 10, height = 5)


# XGBoost Regression Model for Age Prediction Using TCR/BCR Repertoire Data
# Based on tidymodels framework
# References:
# https://www.tidymodels.org/find/parsnip/
# https://parsnip.tidymodels.org/reference/boost_tree.html
# https://parsnip.tidymodels.org/reference/details_boost_tree_xgboost.html
# Evaluation metrics: https://cran.r-project.org/web/packages/yardstick/vignettes/metric-types.html

library(tidymodels)
library(xgboost)
library(iml)
library(fastshap)
library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(forcats)
library(ggbeeswarm)

set.seed(1234)  # Ensure reproducibility

# 1. Data Loading and Preprocessing
# Load training data (BJH cohort)
rt_bjh <- read.csv("BJH.TCR.BCR.select.freq.csv", row.names = 1, check.names = TRUE) %>%
  mutate(sample = rownames(.))
clin <- read_excel("重新年龄分组整合basic_infomation.xlsx") %>% select(2, 4) %>% rename(sampleID2 = 1, Age = 2)
dat_bjh <- merge(rt_bjh, clin, by.x = "sample", by.y = "sampleID2") %>%
  column_to_rownames("sample") %>% select(-sample)
traindata <- dat_bjh

# Load test data (XKH cohort)
rt_xkh <- read.csv("XKH.TCR.BCR.select.freq.csv", row.names = 1, check.names = TRUE) %>%
  mutate(sample = rownames(.))
dat_xkh <- merge(rt_xkh, clin, by.x = "sample", by.y = "sampleID2") %>%
  column_to_rownames("sample") %>% select(-sample)
testdata <- dat_xkh

# Load COVID test data
rt_covid <- read.csv("新冠.TCR.BCR.select.freq.csv", row.names = 1, check.names = TRUE) %>%
  mutate(sample = rownames(.))
dat_covid <- merge(rt_covid, clin, by.x = "sample", by.y = "sampleID2") %>%
  column_to_rownames("sample") %>% select(-sample)
testdata_COVID <- dat_covid

# Keep common features across all datasets
common_genes <- Reduce(intersect, list(colnames(traindata), colnames(testdata), colnames(testdata_COVID)))
traindata <- traindata[, common_genes]
testdata <- testdata[, common_genes]
testdata_COVID <- testdata_COVID[, common_genes]

# 2. Data Preprocessing with Recipes
datarecipe <- recipe(Age ~ ., data = traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2 <- bake(datarecipe, new_data = NULL) %>% select(Age, everything())
testdata2 <- bake(datarecipe, new_data = testdata) %>% select(Age, everything())

# 3. Model Definition and Training
# Define XGBoost model with hyperparameters to tune
model_xgboost <- boost_tree(
  mode = "regression",
  engine = "xgboost",
  mtry = tune(),
  trees = 1000,
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  stop_iter = 25
) %>% set_args(validation = 0.2)

# Create workflow
wk_xgboost <- workflow() %>%
  add_model(model_xgboost) %>%
  add_formula(Age ~ .)

# 5-fold cross-validation
folds <- vfold_cv(traindata2, v = 5)

# Define hyperparameter search space
hpset_xgboost <- parameters(
  mtry(range = c(2, 8)),
  min_n(range = c(5, 10)),
  tree_depth(range = c(1, 3)),
  learn_rate(range = c(-3, -1)),
  loss_reduction(range = c(-3, 0)),
  sample_prop(range = c(0.8, 1))
)

# Random search with 5 combinations
hpgrid_xgboost <- grid_random(hpset_xgboost, size = 5)

# Tune hyperparameters
tune_xgboost <- wk_xgboost %>%
  tune_grid(resamples = folds,
            grid = hpgrid_xgboost,
            metrics = metric_set(rmse, rsq, mae),
            control = control_grid(save_pred = TRUE, verbose = TRUE))

# Select optimal hyperparameters based on RMSE
hpbest_xgboost <- tune_xgboost %>% select_best(metric = "rmse")

# Train final model with optimal hyperparameters
final_xgboost <- wk_xgboost %>%
  finalize_workflow(hpbest_xgboost) %>%
  fit(traindata2)

# 4. Model Interpretation
# Feature importance
final_xgboost2 <- extract_fit_engine(final_xgboost)
importance_matrix <- xgb.importance(model = final_xgboost2)
xgb.plot.importance(importance_matrix, measure = "Cover")

# SHAP values visualization
xgb.plot.shap(data = as.matrix(traindata2[,-1]), 
              model = final_xgboost2, top_n = 5)

# 5. Prediction and Evaluation
# Predict on training set
predtrain_xgboost <- final_xgboost %>%
  predict(new_data = traindata2) %>%
  bind_cols(traindata2 %>% select(Age)) %>%
  mutate(dataset = "train")

# Predict on test sets (XKH and COVID cohorts)
predtest_xgboost <- final_xgboost %>%
  predict(new_data = testdata2) %>%
  bind_cols(testdata2 %>% select(Age)) %>%
  mutate(dataset = "test", model = "xgboost")

predtest_xgboost_COVID <- final_xgboost %>%
  predict(new_data = testdata_COVID) %>%
  bind_cols(testdata_COVID %>% select(Age)) %>%
  mutate(dataset = "test_COVID", model = "xgboost")

# Correlation analysis
cor.test(predtrain_xgboost$Age, predtrain_xgboost$.pred, method = "spearman")
cor.test(predtest_xgboost$Age, predtest_xgboost$.pred, method = "spearman")
cor.test(predtest_xgboost_COVID$Age, predtest_xgboost_COVID$.pred, method = "spearman")

# 6. Visualization of Predictions
# Training set predictions
predtrain_xgboost %>%
  ggplot(aes(x = Age, y = .pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1.2) +
  labs(x = "Actual Age", y = "Predicted Age", title = "Training Set Predictions") +
  theme_bw()

# Test set (XKH) predictions
predtest_xgboost %>%
  ggplot(aes(x = Age, y = .pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1.2) +
  labs(x = "Actual Age", y = "Predicted Age", title = "XKH Test Set Predictions") +
  theme_bw()

# COVID test set predictions
predtest_xgboost_COVID %>%
  ggplot(aes(x = Age, y = .pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1.2) +
  labs(x = "Actual Age", y = "Predicted Age", title = "COVID Test Set Predictions") +
  theme_bw()

# 7. Comprehensive Evaluation
# Combine evaluation metrics
eval_xgboost <- predtrain_xgboost %>%
  bind_rows(predtest_xgboost) %>%
  group_by(dataset) %>%
  metrics(truth = Age, estimate = .pred) %>%
  mutate(model = "xgboost")

# Cross-validation results for optimal hyperparameters
eval_best_cv5_xgboost <- tune_xgboost %>%
  collect_predictions() %>%
  inner_join(hpbest_xgboost[, 1:6]) %>%
  group_by(id) %>%
  metrics(truth = Age, estimate = .pred) %>%
  mutate(model = "xgboost")

# Save results
save(final_xgboost, predtest_xgboost, eval_xgboost, eval_best_cv5_xgboost, 
     file = "evalresult_xgboost.RData")

# 8. Advanced Model Interpretation with IML
traindatax <- traindata2[, -1]
predictor_model <- Predictor$new(final_xgboost, data = traindatax, y = traindata2$Age)

# Permutation feature importance
imp_model <- FeatureImp$new(predictor_model, loss = "rmse")
imp_model$plot() + theme_bw()

# Partial dependence plots
effs_model <- FeatureEffects$new(predictor_model, method = "pdp")
effs_model$plot()

# SHAP analysis with fastshap
shap <- explain(final_xgboost, 
                X = as.data.frame(traindatax),
                nsim = 10,
                adjust = TRUE,
                pred_wrapper = function(model, newdata) {predict(model, newdata) %>% pull(1)})

# SHAP summary plot
data1 <- shap %>%
  as.data.frame() %>%
  mutate(id = 1:n()) %>%
  pivot_longer(cols = -(ncol(traindatax)+1), values_to = "shap")
shapimp <- data1 %>%
  group_by(name) %>%
  summarise(shap.abs.mean = mean(abs(shap))) %>%
  arrange(shap.abs.mean) %>%
  mutate(name = as_factor(name))
data2 <- traindatax %>%
  mutate(id = 1:n()) %>%
  pivot_longer(cols = -(ncol(traindatax)+1))

data1 %>%
  left_join(data2) %>%
  rename("feature" = "name") %>%
  group_by(feature) %>%
  mutate(value = (value - min(value)) / (max(value) - min(value)),
         feature = factor(feature, levels = levels(shapimp$name))) %>%
  ggplot(aes(x = shap, y = feature, color = value)) +
  geom_quasirandom(width = 0.2) +
  scale_color_gradient(low = "red", high = "blue", 
                       breaks = c(0, 1), labels = c("Low", "High")) +
  labs(x = "SHAP value", color = "Feature value") +
  theme_bw() +
  theme(legend.title = element_text(angle = -90))



