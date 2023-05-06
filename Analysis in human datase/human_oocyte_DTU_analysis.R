###### Use kallisto counts for human oocyte do DTU analysis ######

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
library(ggsci)
library(tximport)
library(IsoformSwitchAnalyzeR)
library(org.Hs.eg.db)


# get the path for kallisto counts files
files = list.files("kallisto_counts_after_tiding_transcript_id/", pattern = ".tsv", full.names = T)
files = files[1:10]
filename = basename(files) %>% str_remove("\\.tsv")
names(files) = filename

# prepare sample attribute data frame
samples = tibble(sampleID = filename,
                 condition = c(rep("GV", 5), rep("MII",5)))


# import Kallisto counts matrix
kallistoQuant = importIsoformExpression(sampleVector = files,
                                        addIsofomIdAsColumn = T)

# create aSwitchlist object
# annotation GTF file can be downloaded and cDNA fastq file
# can be downloaed from genecode database uisng the following link
# https://www.gencodegenes.org/human/
aSwitchList = importRdata(isoformCountMatrix = kallistoQuant$counts,
                          isoformRepExpression = kallistoQuant$abundance,
                          designMatrix = samples,
                          isoformExonAnnoation = "./annotation file/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz",
                          isoformNtFasta = "./annotation file/gencode.v41.transcripts.fa.gz")

# filter out single isoform genes or non-expressed isoforms
trim_switchlist = preFilter(switchAnalyzeRlist = aSwitchList,
                            geneExpressionCutoff = 1,
                            isoformExpressionCutoff = 0,
                            removeSingleIsoformGenes = T)

# identifying isoform switches
dtu_result = isoformSwitchTestDEXSeq(switchAnalyzeRlist = trim_switchlist,
                                     reduceToSwitchingGenes = T)

# have a quick glance of isoform switching events
extractSwitchSummary(dtu_result, dIFcutoff = 0.1, alpha = 0.01)

# prepare colors for visulization
color = c("#ED000099","#00468B99","#ADB6B699")

# draw a volcano plot to show isoform switcing events
p1 = dtu_result$isoformFeatures %>%
  mutate(DTU = case_when(isoform_switch_q_value < 0.01 & dIF > 0.1 ~ "Up",
                         isoform_switch_q_value < 0.01 & dIF < -0.1 ~ "Down",
                         TRUE ~ "No difference")) %>%
  ggplot(aes(dIF, -log10(isoform_switch_q_value), col = DTU )) +
  geom_point(position = "jitter") + geom_hline(yintercept = -log10(0.01), linetype = "dashed", size =0.8) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", size = 0.8) +
  theme_classic() +  scale_color_manual(values = c(color[2], color[3], color[1])) +
  theme(axis.title.x.bottom = element_text(size =13, face = "bold"),
        axis.text.x.bottom = element_text(face = "bold"),
        axis.title.y.left = element_text(size =13, face = "bold"),
        axis.text.y.left = element_text(face = "bold"),
        axis.line.x.bottom = element_line(size = 0.8),
        axis.line.y.left = element_line(size = 0.8),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, vjust = 1.5,
                                  face = "bold", size = 15),
        legend.position = "top") +
  ggtitle("DTU_MII_Vs_GV") +
  labs(y = "-log10(FDR)") +
  scale_y_continuous(limits = c(0,25))

saveRDS(p1, "human_DTU_MII_Vs_GV.rds")

# read DGE analysis results in
DEG_toptable = data.table::fread("./human_MII_Vs_GV.csv") %>%
  mutate(DE = case_when(adj.P.Val < 0.01 & logFC > 1.5 ~ "Up",
                        adj.P.Val < 0.01 & logFC < -1.5~ "Down",
                        T ~ "No difference"))

# get the gene name of those differential expression genes
DEG = DEG_toptable %>%
  subset(DE != "No difference") %>% .$ID.gene_name %>%
  na.omit() %>% unique()

# get the gene name of genes with isoform switching events
gDTU = dtu_result$isoformFeatures %>%
  subset(isoform_switch_q_value < 0.01 & abs(dIF) > 0.1) %>%
  .$gene_name %>% unique() %>% na.omit()

# check the intersection between DEG and genes with DTU events
dplyr::intersect(DEG, gDTU)


# prepare gene symbol to entrezid conversion table
gene_info = AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c("ENTREZID", "SYMBOL"))

# get those genes only identified by DTU analysis
gDTU_only = tibble(gene_name = dplyr::setdiff(gDTU, DEG)) %>%
  left_join(gene_info, by = c("gene_name" = "SYMBOL")) %>%
  .$ENTREZID %>% unique() %>% na.omit()

# all genes after QC remained for DGE analysis as background
un_genes = tibble(gene_name = DEG_toptable$ID.gene_name) %>%
  left_join(gene_info, by = c("gene_name" = "SYMBOL")) %>%
  .$ENTREZID %>% unique() %>% na.omit()

# use those genes only identified by DTU analysis to do
# GO enrichment analysis
GO_gDTU = enrichGO(
  gDTU_only,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = un_genes,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T)

# visualize GO enrichment result
p2= enrich.plot(GO_gDTU@result, fraction.factor = 20,
            show.size = 10,
            text.location = 10, str.wrap.size = 25, save = F,
            show.legend = F)

saveRDS(p2, "human_DTU_GO.rds")

# summary the biotype of genes with isoform switch events
biotype_summary = dtu_result$isoformFeatures %>%
  subset(isoform_switch_q_value < 0.01 & abs(dIF) > 0.1) %>%
  group_by(gene_biotype) %>% tally()

biotype_summary = biotype_summary %>% mutate(percentage = n/sum(n))

# bar plot of biotype summary
p3= biotype_summary %>% Barplot(x = "gene_biotype", y = "percentage")
saveRDS(p3, "human_DTU_biotype_summary.rds")


# extract sequences for external tools (CPC2, pfam, SignalIP and Iupred2A)
# The use of those external tools can follow tutorial of IsoformSwitchAnalyzeR
# package in this link
# https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
dtu_result = extractSequence(dtu_result, pathToOutput = "./files for external tools/",
                             removeLongAAseq = T, alsoSplitFastaFile = T,
                             alpha = 0.01, dIFcutoff = 0.1, writeToFile = T)

# add CPC2 analysis
dtu_result = analyzeCPC2(switchAnalyzeRlist = dtu_result,
                         pathToCPC2resultFile = "./output from external tools/CPC2/cpc2_integrated.txt",
                         codingCutoff = 0.725,
                         removeNoncodinORFs = T)

# add Pfam analysis
dtu_result = analyzePFAM(switchAnalyzeRlist = dtu_result,
                         pathToPFAMresultFile = "./output from external tools/pfam/human_pfam.txt",
                         showProgress = F)

# add SignalP analysis
dtu_result = analyzeSignalP(switchAnalyzeRlist = dtu_result,
                            pathToSignalPresultFile = list.files("./output from external tools/signalIP/", pattern = "txt", full.names = T))

# add Iupred2A analysis
dtu_result = analyzeIUPred2A(switchAnalyzeRlist = dtu_result,
                             pathToIUPred2AresultFile = list.files("./output from external tools/IUPred2A/", pattern = "result", full.names = T))


# predict alternative splicing
dtu_result = analyzeAlternativeSplicing(switchAnalyzeRlist = dtu_result)

# predict switch consequence
dtu_result = analyzeSwitchConsequences(dtu_result, dIFcutoff = 0.1)


# extract predictive splicing events
splicing_summary = extractSplicingSummary(dtu_result, returnResult = T, plot = F)
splicing_summary = splicing_summary %>% mutate(splicingResult = str_extract(splicingResult, "less|more"))


consequence_summary = extractConsequenceSummary(dtu_result, returnResult = T, plot = F)
consequence_summary = consequence_summary %>%
  mutate(switchConsequence = str_extract(switchConsequence, "gain|loss|switch|insensitive|sensitive|Noncoding|coding|shorter|longer")) %>%
  mutate(switchConsequence = factor(switchConsequence, levels = c("gain", "loss", "switch","longer", "shorter",
                                                                  "coding", "Noncoding", "sensitive", "insensitive")))

# visualize splicing events and functional consequence results
p4 = stackBar(splicing_summary, x = "AStype", y = "nrIsoWithConsequences", fill = "splicingResult",
              colors = c(color[2], color[1]), legend = T)

p5 = stackBar(consequence_summary, x = "featureCompared", y = "nrIsoWithConsequences", fill = "switchConsequence",
              colors = pal_npg("nrc")(10), legend = T)

# viusalize differential splicing events and fuctional consequence results
p6 = extractSplicingEnrichment(dtu_result, returnResult = F, returnSummary = F, countGenes = F)

p7 = extractConsequenceEnrichment(dtu_result, returnResult = F, countGenes = F)


saveRDS(p4, "human_splicing_bar.rds")
saveRDS(p5, "human_consequence_bar.rds")
saveRDS(p6, "human_splicing_enrichment.rds")
saveRDS(p7, "human_consequence_enrichment.rds")
