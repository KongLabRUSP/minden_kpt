# |---------------------------------------------------------------------------------|
# | Project:  KPT                                                                   |
# | Script:   DNA vs. RNA                                                           |
# | Author:   Audrey Minden, Davit Sargsyan                                         |
# | Created:  11/02/2018                                                            |
# | Modified:                                                                       | 
# |---------------------------------------------------------------------------------|
# sink(file = "tmp/log_minden_ngs_v1.txt")
# Reference:
# https://rstudio-pubs-static.s3.amazonaws.com/222140_9ae622e197da48ebbe93d7044d8bdc36.html  

# Header----
require(data.table)
require(ggplot2)
require(ggdendro)
require(DESeq2)
require(DEGseq)

# Read data----
# # Format the input file----
# rna <- fread("data/minden_rna.csv",
#              header = FALSE)
# rna[[1]] <- gsub(pattern = "\"",
#                  replacement = "",
#                  x = rna[[1]])
# 
# rna[[3]] <- gsub(pattern = "\"",
#                  replacement = "",
#                  x = rna[[3]])
# rna
# write.csv(rna,
#           file = "data/minden_rna.csv",
#           row.names = FALSE,
#           col.names = FALSE)
# rm(rna)
# gc()

rna <- fread("data/minden_rna.csv")
rna

# Remove genes with all zeros----
ndx.keep <- rowSums(rna[, -1]) > 0
sum(ndx.keep)
# 19,046 our of 26,364 genes kept
rna <- rna[ndx.keep, ]
rna

# Remove low counts----
row.counts <- rowSums(rna[, -1])
# Keep only the once with at least 100 counts per row----
rna <- rna[row.counts >= 100, ]
# 11,018, down from 26,364

# Remove duplicate genes----
genes.rm <- rna$gene[duplicated(rna$gene)]
genes.rm
rna[gene %in% genes.rm, ]
rna <- rna[!(gene %in% genes.rm), ]
# Down to 11,016 genes

# Normalize data to Fragments per million (FPM)----
mat <- data.frame(sample = colnames(rna)[-1],
                  trt = colnames(rna)[-1],
                  repl = factor(rep(1, ncol(rna) - 1)))
mat

dtm <- as.matrix(rna[, -1, with = FALSE])
rownames(dtm) <- rna$gene
head(dtm)

dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ trt)
dds <- estimateSizeFactors(dds)
dds

# Fragments per million (FPM) normalization----
rna.fpm <- data.table(gene = rna$gene,
                      fpm(dds,
                          robust = FALSE))
colnames(rna.fpm)[-1] <- paste(colnames(rna.fpm)[-1],
                               "fpm",
                               sep = "_")
rna.fpm

# DEGSeq----
DEGexp(geneExpMatrix1 = rna,
       geneCol1 = which(colnames(rna) == "gene"), 
       expCol1 = which(colnames(rna) == "K2"), 
       groupLabel1 = "RA",
       
       geneExpMatrix2 = rna,
       geneCol2 = which(colnames(rna) == "gene"), 
       expCol2 = which(colnames(rna) == "D3"),
       groupLabel2 = "C",
       
       foldChange = 2,
       qValue = 0.01,
       thresholdKind = 5, 
       rawCount = TRUE,
       normalMethod = "none",
       method = "MARS",
       outputDir = "tmp")

k2_d3 <- fread("tmp/output_score.txt")
k2_d3
k2_d3[k2_d3$`Signature(q-value(Storey et al. 2003) < 0.01)`,]

# # CHECKPOINT: raw vs normalized values----
# tmp <- k2_d3[, c("GeneNames",
#                 "value1",
#                 "value2")]
# colnames(tmp)[1] <- "gene"
# tmp <- merge(rna, tmp, by = "gene")
# tmp
# plot(tmp$C ~ tmp$value2)
# # Sum of square errors
# sum((tmp$C - tmp$value2)^2,
#     na.rm = TRUE)
# # 0, i.e. raw counts are printed

# Write as CSV----
write.csv(k2_d3,
          file = "tmp/rna_diff.csv",
          row.names = FALSE)


# Merge all data----
dt1 <- merge(rna,
             rna.fpm,
             by = "gene")
dt1

tmp <- k2_d3[, c("GeneNames",
                "log2(Fold_change) normalized",
                "q-value(Storey et al. 2003)")]
colnames(tmp) <- c("gene",
                   "log2(K2-D3)",
                   "q-val(K2-D3)")
dt1 <- merge(dt1,
             tmp,
             by = "gene")
dt1

# Save combined data file----
write.csv(dt1,
          file = "tmp/dt1.csv")
# Clean memory
rm(list = ls()[ls() != "dt1"])
gc()

# # Select genes----
# hist(log2(dt1$D3_fpm + 1), 100)
# qnt <- quantile(x = log2(dt1$D3_fpm + 1),
#                 prob = c(0.005, 0.995))
# qnt
# # 7.217905
# qnt <- quantile(x =dt1$D3_fpm,
#                 prob = 0.90)
# qnt
#   # 0.5%     99.5% 
#   # 1.451462 10.586535
# 
# # Select genes with high expression in Control group (90st quantile)
# # significant differences in RA groups vs Control, and no significant 
# # differences in non-RA groups vs. control, all at alpa = 0.01.
# dt2 <- dt1[dt1$`q-val(K2-D3)` <= 0.05 &
#              log2(dt1$D3_fpm + 1) >= qnt, ]
# dt2

dt2 <- dt1[abs(dt1$`log2(K2-D3)`) >= 4 &
             dt1$`q-val(K2-D3)` <= 0.01, ]

# Heatmap of FPM----
# Distances----
tmp <- as.matrix(dt2[, c("D3_fpm",
                         "K2_fpm")])
rownames(tmp) <- dt2$gene
geneDist <- dist(tmp[, -1])

# Make dendrogram data----
dhc <- as.dendrogram(hclust(d = geneDist),
                     horiz = TRUE)
ddata <- dendro_data(dhc, 
                     type = "rectangle")

# Segment data----
dtp1 <- segment(ddata)
# Hitmap data----
dtp2 <- melt.data.table(dt2,
                        id.vars = "gene",
                        measure.vars = c("D3_fpm",
                                         "K2_fpm"),
                        variable.name = "Treatment",
                        value.name = "RNA")
dtp2$gene <- factor(dtp2$gene,
                    levels = ddata$labels$label)

dtp2

offset.size <- 2
p1 <- ggplot(data = dtp2) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_segment(data = dtp1,
               aes(x = -log(y) -2.5,
                   y = x, 
                   xend = -log(yend) -2.5,
                   yend = xend),
               size = 1) +
  geom_tile(aes(x =  as.numeric(Treatment),
                y = gene, 
                fill = log2(RNA + 1)),
            color = "white") +
  geom_text(data = dtp2[Treatment == "D3_fpm", ],
            aes(x = rep(2.75,
                        nlevels(gene)),
                y = gene,
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(gene))[as.numeric(gene)] + 
                  offset.size,
                label = unique(gene)),
            hjust = 0) +
  geom_text(data = dtp2[gene == levels(dtp2$gene)[1], ],
            aes(x = 1:nlevels(Treatment),
                y = rep(-offset.size,
                        nlevels(Treatment)),
                angle = 0,
                label = gsub(x = levels(Treatment),
                             pattern = "_fpm",
                             replacement = "")),
            hjust = 1) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "log2(FPM)") +
  scale_x_continuous("",
                     breaks = as.numeric(dt2$Treatment),
                     labels = unique(dt2$Treatment)) +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("RNA Normalized Counts (FPM) in Differentially Expressed Genes") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p1

tiff(filename = "tmp/heatmap_fpm.tiff",
     height = 10,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Keep RNA data for the selected genes----
rna <- dt2
rm(list = ls()[ls() != "rna"])
gc()

# DNA----
dna <- fread("data/minden_dna.csv")
dna

# In DNA data, keep only the genes selected from RNA expressions----
dna <- dna[gene %in% rna$gene, ]
dna

length(unique(dna$gene))
# 153

length(unique(rna$gene))
# 160

# Separate genes with meaningfull (>10%) differences----
dna$dMethyl <- 100*(dna$`7.KPT9274` - dna$`6.Control`)
hist(dna$dMethyl, 100)

plot(dna$dMethyl ~ dna$`6.Control`)

gene.keep <- unique(dna$gene[abs(dna$dMethyl) >= 10])
gene.keep
# 29 genes

dna <- dna[gene %in% gene.keep, ]
dna[, distRank := rev(rank(distance)),
    by = gene]

dna

# Transform to Long format----
dt1 <- dna[, c("gene",
               "CpG",
               "feature",
               "distance",
               "distRank",
               "dMethyl"), ]
dt1

dt1$feature[substr(dt1$feature, 1, 4) == "Exon"] <- "Exon"
dt1$feature[substr(dt1$feature, 1, 6) == "Intron"] <- "Intron"
dt1$feature[substr(dt1$feature, 1, 8) == "Promoter"] <- "Promoter"
dt1$feature[substr(dt1$feature, 1, 4) == "Down"] <- "Downstream"
dt1$feature <- factor(dt1$feature)
summary(dt1$CpG)

dt1$reg <- "5 to 10"
dt1$reg[dt1$CpG > 10] <- "11 to 20"
dt1$reg[dt1$CpG > 20] <- ">20"
dt1$reg <- factor(dt1$reg,
                      levels = c("5 to 10",
                                 "11 to 20",
                                 ">20"))
summary(dt1)

# Merge DNA with RNA----
dt1 <- merge(dt1,
             rna,
             by = c("gene"))
dt1

# Starburst plot----
for (i in 1:length(unique(dt1$Treatment))) {
  cmp <- unique(dt1$Treatment)[i]
  
  
  tmp <- dt1[dt1$Treatment == cmp, ]
  tmp
  
  g1 <- dt1[dMethyl >= 10 & 
              `log2(K2-D3)` < 0 &
              feature == "Promoter", ]
  g1
  length(unique(g1$gene))
  # 5 genes
  
  g2 <- dt1[dMethyl <= -10 & 
              `log2(K2-D3)` > 0 &
              feature == "Promoter"]
  g2
  length(unique(g2$gene))
  # 0 genes
  
  p1 <- ggplot(data = tmp,
               aes(x = DNA,
                   y = RNA,
                   fill = feature)) +
    geom_point(alpha = 0.7,
               size = 2,
               shape = 21) +
    geom_text(data = unique(tmp[gene %in% unique(g1$gene) &
                                  feature == "Promoter",
                                c("gene",
                                  "feature",
                                  "RNA")]),
              aes(x = 40,
                  y = RNA,
                  label = gene),
              color = "blue",
              size = 2) +
    geom_text(data = unique(tmp[gene %in% unique(g2$gene) &
                                  feature == "Promoter", 
                                c("gene",
                                  "feature",
                                  "RNA")]),
              aes(x = -40,
                  y = RNA,
                  label = gene),
              color = "blue",
              size = 2) +
    geom_hline(yintercept = c(-0.5, 0.5),
               linetype = "dashed") +
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed") +
    scale_x_continuous("DNA Methylation Difference(%)",
                       breaks = seq(-50, 30, 10)) +
    scale_y_continuous("RNA Expression Difference (log2)",
                       breaks = seq(-5, 10, 1)) +
    ggtitle(cmp) +
    # scale_fill_manual("Region",
    #                   values = c("Promoter" = "green",
    #                              "5' UTR" = "white",
    #                              "Body" = "blue",
    #                              "3' UTR" = "grey",
    #                              "Downstream" = "red")) +
    theme(plot.title = element_text(hjust = 0.5))
  p1
  tiff(filename = paste("tmp/",
                        cmp,
                        ".tiff",
                        sep = ""),
       height = 10,
       width = 10,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

# Isolate genes----
for (i in 1:length(unique(dna$gene))) {
  gX <- unique(dt1$gene)[i]
  dna.gX <- dt1[dt1$gene %in% gX, ]
  dna.gX$y0 <- 0
  
  dna.gX$Treatment <- paste(dna.gX$Treatment,
                            " (RNA = ",
                            round(dna.gX$RNA, 3),
                            ")",
                            sep = "")
  
  p1 <- ggplot(dna.gX,
               aes(x = distRank,
                   y = DNA)) +
    facet_wrap(.~ Treatment,
               scales = "free_y",
               ncol = 1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = -Inf,
                  ymax = -10),
              fill = "red",
              alpha = 0.1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = 10,
                  ymax = Inf),
              fill = "green",
              alpha = 0.1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(10,
                              -10),
               linetype = "dashed") +
    geom_segment(aes(x = distRank,
                     y = y0,
                     xend = distRank,
                     yend = DNA)) + 
    geom_point(aes(x = distRank,
                   y = DNA,
                   fill = feature,
                   size = reg),
               shape = 21) +
    # geom_rect(aes(xmin = -Inf,
    #               xmax = Inf,
    #               ymin = -10,
    #               ymax = 10),
    #           fill = "white",
    #           alpha = 0.1) +
    ggtitle(paste("Gene:",
                  gX)) +
    scale_x_continuous("Distance from TSS",
                       breaks = dna.gX$distRank,
                       labels = dna.gX$distanceToTSS) +
    scale_y_continuous("% Methylation") +
    scale_fill_manual("Region",
                      values = c("Distal Intergenic" = "purple",
                                 "Exon" = "blue",
                                 "Intron" = "white",
                                 "Promoter" = "brown",
                                 "3' UTR" = "black",
                                 "5' UTR" = "yellow",
                                 "Downstream" = "orange")) +
    scale_size_manual("Number of CpG-s",
                      values = c("5 to 10" = 5,
                                 "11 to 20" = 6,
                                 ">20" = 7)) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position = "top",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  p1
  tiff(filename = paste("tmp/",
                        gX,
                        ".tiff",
                        sep = ""),
       height = 8,
       width = 8,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

# sessionInfo()
# sink()