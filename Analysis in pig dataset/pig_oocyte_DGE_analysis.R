###### Use kallisto counts for pig oocyte do DEG analysis ######

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
library(org.Ss.eg.db)
library(ggrepel)
library(patchwork)
library(RColorBrewer)

# prepare colors for visualization
colors = c(pal_lancet("lanonc")(9),
           pal_npg("nrc")(10)[c(1:7,9,10)],
           pal_nejm("default")(8))

# prepare annotation files for tximport function
ah = AnnotationHub()
ah = ah %>% subset(species == "Sus scrofa" & rdataclass == "EnsDb")
ensdb = ah[["AH104988"]]

tx2gene = transcripts(ensdb) %>% mcols() %>%
  as.data.frame() %>%
  dplyr::select(tx_id_version, gene_id)

geneID2symbol = genes(ensdb) %>% mcols() %>%
  as.data.frame() %>% dplyr::select(gene_id, gene_name, description, gene_biotype)

rownames(geneID2symbol) = NULL

# get the location of Kallisto quantification files
files = list.files("kallisto_quant/", pattern = ".tsv", full.names = T)
filename = basename(files) %>% str_remove("\\.tsv")
names(files) = filename

# import kallisto counts files using tximport function
txi = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene,
                         countsFromAbundance = "lengthScaledTPM")

# create dgelist object
# remove '319' (GV) sample since it only contains 13M sequencing reads which is only
# half of other sequencing samples
dge = DGEList(txi$counts %>% as.data.frame()) %>% calcNormFactors()


# Since Damage samples and MII samples are very similar
# combine them as one big group to do downstream analysis
dge$samples$group = c(rep(c("GV", "Mix", "Mix"), 2), "Mix", "Mix",
                      rep(c("GV", "Mix", "Mix"), 2))

# add gene attributes to DGElist
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

## PCA analysis on clean data

# have a copy of filter_dge object as filter_dge2
# change the column names, row nmaes, group column in samples array of filter_dge2
filter_dge2 = filter_dge
colnames(filter_dge2$counts) =  c("GV_1", "MII_1", "Damage_1",
                                             "GV_2", "MII_2", "Damage_2",
                                             "MII_3", "Damage_3",
                                             "GV_3", "MII_4", "Damage_4",
                                             "GV_4", "MII_5", "Damage_5")

rownames(filter_dge2$samples) = c("GV_1", "MII_1", "Damage_1",
                                  "GV_2", "MII_2", "Damage_2",
                                  "MII_3", "Damage_3",
                                  "GV_3", "MII_4", "Damage_4",
                                  "GV_4", "MII_5", "Damage_5")

filter_dge2$samples$group = c("GV", "MII", "Damage",
                              "GV", "MII", "Damage",
                              "MII", "Damage",
                              "GV", "MII", "Damage",
                              "GV", "MII", "Damage")

# create a PCA_analysis plot
p1 = PCA_analysis(filter_dge2, save = F, show.legend = F, plot = T)

saveRDS(p1, "pig_pca.rds")

# pearson correlation to show the samples correlation
cor_mat = cor(filter_dge %>% cpm(log = T))

# prepare side annotation table
annotation_df = tibble(sample = colnames(dge$counts),
                       type = rep(c(rep(c("GV", "MII", "Damage"), 2),
                                    "MII", "Damage",
                                    rep(c("GV", "MII", "Damage"), 2))))

annotation_df = annotation_df %>% column_to_rownames("sample")

# prepare side annotation colors
annotation_color = list(type = c(GV = colors[1],  MII = colors[3],
                                 Damage = colors[10]

                                 ))

# draw a heatmap of pearson distance martix
png("correlation heatmap.png", res = 300, width = 2000, height = 2000)
pheatmap::pheatmap(1-cor_mat, annotation_col = annotation_df,
                   annotation_colors = annotation_color,
                   color = colorRampPalette(c("white", colors[7]))(50)[4:50])
dev.off()

# DGE analysis by limma-voom pipeline
model_matrix <- model.matrix(~ 0 + group, filter_dge$samples)

contr_matrix <- makeContrasts(Mix_stageVsGv_stage = groupMix - groupGV,
                              levels = colnames(model_matrix))

v = voom(filter_dge, model_matrix, plot = T)

top_table_MG = v %>% lmFit(model_matrix) %>% contrasts.fit(contr_matrix) %>%
  eBayes() %>% topTable(coef = "Mix_stageVsGv_stage", number = Inf) %>%
  mutate(DE = case_when(logFC > 2 & adj.P.Val < 0.01 ~ "Up",
                        logFC < -2 & adj.P.Val < 0.01 ~ "Down",
                        TRUE ~ "No difference"
                        ))

# tidy up column names
colnames(top_table_MG) = colnames(top_table_MG) %>% str_remove("ID\\.")

# write DE analysis result to disk
write.csv(top_table_MG %>% dplyr::select(-entrezid), row.names = F,
          file = "DE result(Mix vs GV).csv")

saveRDS(top_table_MG, "top_tabe_MG.rds")


# volcano plot to show DE genes
p2 = Volcano.plot(top_table_MG, title = "DEG_Mix_Vs_GV", save = T,
             pvalue.cutoff = 0.01, log2fc.cutoff = 2)

saveRDS(p2, "pig_volcano.rds")

# counts the DE genes
top_table_MG %>% dplyr::count(DE)


# GO/KEGG analysis based on only up-regulated genes or only
# down-regulated genes
up_genes = top_table_MG %>% subset(logFC > 2 & adj.P.Val < 0.01) %>%
  subset(!is.na(entrezid)) %>% .$entrezid %>% unlist() %>% unique() %>% as.character()

down_genes = top_table_MG %>% subset(logFC < -2 & adj.P.Val < 0.01) %>%
  subset(!is.na(entrezid)) %>% .$entrezid %>% unlist() %>% unique() %>% as.character()

# all genes remained after QC as background for GO enrichment analysis
un_genes = top_table_MG %>% .$entrezid %>% unlist() %>% unique() %>%
  na.omit() %>%
  as.character()

# do GO enrichment analysis
GO_DGE = lapply(list(down_genes, up_genes), function(x) enrichGO(
  x,
  org.Ss.eg.db,
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

# visualize GO enrichment results
p3 = GO_DGE[["down"]]@result  %>%
  enrich.plot(str.wrap.size = 35, text.location = 8,
              save = F, fraction.factor = 12,
              file.name = "GO_down_regulated_genes",
              show.size = 10,
              bar.font.size = 8,
              show.legend = F)

saveRDS(p3, file = "pig_GO_down.rds")

p4 = GO_DGE[["up"]]@result  %>%
  enrich.plot(str.wrap.size = 35, save = F, fraction.factor = 15,
              text.location = 7.5,
              file.name = "GO_up_regulated_genes",
              show.size = 10,
              bar.font.size = 8,
              show.legend = F)

saveRDS(p4, file = "pig_GO_up.rds")


# use fisher exact test check whether there is significant
# overlapping regarding the number of DEGs on a specific pathway
# between pig and human dataset

# identifical pathways enriched for DEGS in both pig and human datasets
pathways = c("GO:0007059", "GO:0098687", "GO:0022626", "GO:0003735")

# read GO enrcihment results of human dataset in
human_GO_up  = readRDS("../human_oocyte/human_GO_DGE.rds")
human_GO_dataset = human_GO_up@result %>% subset(ID %in% pathways)
human_geneset = human_GO_up@geneSets %>% subset(names(.) %in% pathways)

# Since two candidate pathways enriched for up-regulated genes, the other
# two for down-regulated genes in pig dataset, we retrieve theme from each
# GO result dataframe and combine them together to form a new dataframe
pig_dataset = purrr::reduce(list(GO_DGE[["up"]]@result %>% subset(ID %in% pathways[1:2]),
                                 GO_DGE[["down"]]@result %>% subset(ID %in% pathways[3:4])),
                            rbind)

pig_geneset = GO_DGE[["up"]]@geneSets %>% subset(names(.) %in% pathways)

# do fisher exact test to check pathway overlap
Overlap_res = pathway_overlap(pig_dataset, pig_geneset, "pig",
                              human_GO_dataset, human_geneset, "human")

# visualize result by pheatmap
overlap_pathway_plot = Overlap_res[[2]] %>%
  left_join(pig_dataset |> dplyr::select(ID, Description), by = c("GOID_1" = "ID")) %>%
  mutate(description_1 = str_to_sentence(Description)) %>%
  left_join(pig_dataset |> dplyr::select(ID, Description), by = c("GOID_2" = "ID")) %>%
  mutate(description_2 = str_to_sentence(Description.y)) %>%
  ggplot(aes(description_1, description_2)) +
  geom_tile(aes(fill = -log10(p.value)), colour = "black", size = 0.8) +
  theme_classic() +
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                        "RdYlBu")))(10),
                       name = "-Log10(p.value)") +
  geom_text(aes(label = round(-log10(p.value),2))) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text = element_text(face = "bold", colour = "black"),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

saveRDS(overlap_pathway_plot, "pathway overlap heatmap.rds")
