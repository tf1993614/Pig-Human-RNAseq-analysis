###### Use kallisto counts for human oocyte do DEG analysis ######

# load necessary packages
library(tidyverse)
library(DRIMSeq)
library(stageR)
library(AnnotationHub)
library(ggsci)
library(clusterProfiler)
library(GenomicFeatures)
library(limma)
library(edgeR)
library(VennDiagram)
library(tximport)
library(ggrepel)
library(data.table)
library(rgl)
library(org.Hs.eg.db)
library(ggrepel)

# prepare colors for visualization
colors = c(pal_lancet("lanonc")(9),
           pal_npg("nrc")(10)[c(1:7,9,10)],
           pal_nejm("default")(8))

# prepare annotation files for tximport function
ah = AnnotationHub()
ah = ah %>% subset(species == "Homo sapiens" & rdataclass == "EnsDb")
ensdb = ah[["AH109336"]]

tx2gene = transcripts(ensdb) %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(tx_id_version, gene_id)

geneID2symbol = genes(ensdb) %>% mcols() %>%
  as.data.frame() %>% dplyr::select(gene_id, gene_name, description)

rownames(geneID2symbol) = NULL


# get the location of Kallisto quantification files
files = list.files("kallisto_counts/", pattern = ".tsv", full.names = T)
filename = basename(files) %>% str_remove("\\.tsv")
names(files) = filename

# since using the cDNA fastq files from genecode database for estabilshing
# the kallisto index, the quantification result would not only contain the
# Ensembl transcript id but also transcript ids in other format. We need to
# only keep Ensebl transcript id for downstream application

kallisto_files = map(files, ~ fread(.x))
kallisto_files = map(kallisto_files, function(x){
  data = x
  data$target_id = map_chr(data$target_id, function(x){
    res = str_split(x, "\\|") %>% unlist()
    return(res[1])
  })
  return(data)
})

# write the files with transcript id correction into disk
walk(names(kallisto_files), ~ write_tsv(kallisto_files[[.x]],file = str_c(.x, ".tsv")))

# get the location of files after tidying up the transcript id
files = list.files("kallisto_counts_after_tiding_transcript_id/", pattern = ".tsv", full.names = T)
filename = basename(files) %>% str_remove("\\.tsv")
names(files) = filename

# import kallisto counts files using tximport function
txi = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene,
                         countsFromAbundance = "lengthScaledTPM")

# create dgelist object
dge = DGEList(txi$counts[,c(1:10)] %>% as.data.frame()) %>% calcNormFactors()

dge$samples$group = c("GV", "GV", "GV", "GV", "GV",
                      "MII", "MII", "MII",
                      "MII", "MII")

# add gene attributes into DGElist
geneGR = genes(ensdb)
mcols(geneGR) = mcols(geneGR)[,c("gene_id", "gene_name", "gene_biotype",
                                 "description", "entrezid")]
dge$genes = geneGR[rownames(dge),]

# filter out low-counts genes
condition = filterByExpr(dge, group = dge$samples$group)
filter_dge = dge[condition, keep.lib.size = F] %>% calcNormFactors()

# compare gene expression distribution before/after trimming
par(mfrow = c(1,2))
plotDensities(dge %>% cpm(log=T), legend = FALSE, main = "before trimming")
plotDensities(filter_dge %>% cpm(log = T), legend = FALSE, main = "after trimming")

# pearson correlation to show the samples correlation
cor_mat = cor(filter_dge %>% cpm(log = T))

# tidy up column and row names of this person matrix
colnames(cor_mat) = c(str_c("GV_", as.character(seq_len(5))), str_c("MII_", as.character(seq_len(5))))
rownames(cor_mat) = colnames(cor_mat)

# prepare side annotation file
annotation_df = tibble(sample = colnames(dge$counts),
                       type = rep(c("GV", "MII"), each = 5))

annotation_df = annotation_df %>% column_to_rownames("sample")
rownames(annotation_df) = rownames(cor_mat)

# prepare side annotation color
annotation_color = list(type = c(GV = colors[1], MII = colors[3]))

png("correlation heatmap.png", res = 300, width = 2000, height = 2000)
pheatmap::pheatmap(1-cor_mat, annotation_col = annotation_df,
                   annotation_colors = annotation_color,
                   color = rev(colorRampPalette(c("white", colors[7]))(50)[4:50]))
dev.off()

# PCA analysis and graphing
p1 = PCA_analysis(filter_dge, save = F,
                  plot = T,
                  colors = c("#4DBBD5FF", "#00A087FF"),
                  fileName = "pca", show.legend = F,
                  point.size = 3,
                  reName_samples = rownames(cor_mat))

saveRDS(p1, "human_pca.rds")

# DGE analysis by limma-voom pipeline
model_matrix <- model.matrix(~ 0 + group, filter_dge$samples)

contr_matrix <- makeContrasts(MII_stageVsGv_stage = groupMII - groupGV,
                              levels = colnames(model_matrix))

v = voom(filter_dge, model_matrix, plot = T)

top_table_MG = v %>% lmFit(model_matrix) %>% contrasts.fit(contr_matrix) %>%
  eBayes() %>% topTable(coef = "MII_stageVsGv_stage", number = Inf) %>%
  mutate(DE = case_when(adj.P.Val < 0.01 & logFC > 1.5 ~ "Up",
                        adj.P.Val < 0.01 & logFC < -1.5~ "Down",
                        T ~ "No difference"))

# write DGE analysis into disk
write.csv(top_table_MG %>% dplyr::select(-ID.entrezid),
          file = "human_MII_Vs_GV.csv", row.names = F)


# visualize DGE analysis result
p1 = Volcano.plot(top_table_MG, title = "DEG_MII_Vs_GV", save = T,
                            pvalue.cutoff = 0.01, log2fc.cutoff = 1.5)

saveRDS(p1, "human_volcano.rds")


# GO/KEGG analysis based on only up-regulated or down-regulated genes
up_genes = top_table_MG %>% subset(logFC > 1.5 & adj.P.Val < 0.01) %>%
  subset(!is.na(ID.entrezid)) %>% .$ID.entrezid %>% unlist() %>% unique() %>% as.character()

down_genes = top_table_MG %>% subset(logFC < -1.5 & adj.P.Val < 0.01) %>%
  subset(!is.na(ID.entrezid)) %>% .$ID.entrezid %>% unlist() %>% unique() %>% as.character()

# all genes remained after QC as background for GO enrichment analysis
un_genes = top_table_MG %>% .$ID.entrezid %>% unlist() %>% unique() %>%
  na.omit() %>%
  as.character()

# do GO enrichment analysis
GO_DGE = lapply(list(down_genes, up_genes), function(x) enrichGO(
  x,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = un_genes,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = T,
  pool = T))

names(GO_DGE) = c("down", "up")

# Visualization of GO enrichment result
p2 = GO_DGE[["down"]]@result %>% enrich.plot(str.wrap.size = 35, save = F,
                                        text.location = 16.5,
                                        fraction.factor = 20,
                                        bar.width = 0.2,
                                        file.name = "GO_down_regulated_genes",
                                        show.size = 10, show.legend = T,
                                        bar.font.size = 8)

saveRDS(p2, file = "human_GO_down.rds")

p3 = GO_DGE[["up"]]@result %>% enrich.plot(str.wrap.size = 35, save = F,
                                          text.location = 6.5,
                                          fraction.factor = 10,
                                          file.name = "GO_up_regulated_genes",
                                          show.size = 6,
                                          show.legend = T,
                                          bar.font.size = 8)
saveRDS(p3, file = "human_GO_up.rds")

# export GO enrichment for up-regulated genes for pathway overlap
# analysis between human and pig dataset
saveRDS(GO_DGE[["up"]], "human_GO_DGE.rds")
