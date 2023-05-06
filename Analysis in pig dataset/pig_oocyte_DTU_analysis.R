###### Use kallisto counts for pig oocyte do DTU analysis ######

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
library(ggsci)
library(tximport)
library(IsoformSwitchAnalyzeR)
library(ggrepel)
library(VennDiagram)
library(org.Ss.eg.db)

# get the path for kallisto counts files
files = list.files("kallisto_quant/", pattern = ".tsv", full.names = T)
filename = basename(files) %>% str_remove("\\.tsv")
names(files) = filename

# prepare sample attribute data frame
samples = tibble(sampleID = filename,
                 condition = c(rep(c("GV", "Mix", "Mix"), 2), "Mix", "Mix",
                          rep(c("GV", "Mix", "Mix"), 2)))

# import data from kallisto
kallistoQuant = importIsoformExpression(sampleVector = files,
                                        addIsofomIdAsColumn = T)


#create aSwitchlist object
# annotation GTF file can be downloaded and cDNA fastq file
# can be downloaed from Ensembl database uisng the following link
# https://asia.ensembl.org/info/data/ftp/index.html
aSwitchList = importRdata(isoformCountMatrix = kallistoQuant$counts,
                          isoformRepExpression = kallistoQuant$abundance,
                          designMatrix = samples,
                          comparisonsToMake = tibble(condition_1 = c("GV"),
                                                     condition_2 = c("Mix")),
                          isoformExonAnnoation = "./DTU analysis/Sus_scrofa.Sscrofa11.1.108.gtf.gz",
                          isoformNtFasta = "./DTU analysis/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz")

# filter out single isoform genes or non-expressed isoforms
trim_switchlist = preFilter(switchAnalyzeRlist = aSwitchList,
                            geneExpressionCutoff = 1,
                            isoformExpressionCutoff = 0,
                            IFcutoff = 0.05,
                            removeSingleIsoformGenes = T)

# identifying isoform switches
dtu_result = isoformSwitchTestDEXSeq(switchAnalyzeRlist = trim_switchlist,
                                     reduceToSwitchingGenes = T)

# have a quick glance of isoform switching events
extractSwitchSummary(dtu_result, alpha = 0.01, dIFcutoff = 0.1)

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
  theme_classic() + scale_color_manual(values = c(color[2], color[3], color[1])) +
  theme(axis.title.x.bottom = element_text(size =13, face = "bold"),
        axis.text.x.bottom = element_text(face = "bold", colour = "black"),
        axis.title.y.left = element_text(size =13, face = "bold"),
        axis.text.y.left = element_text(face = "bold", colour = "black"),
        axis.line.x.bottom = element_line(size = 0.8),
        axis.line.y.left = element_line(size = 0.8),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, vjust = 1.5,
                                  face = "bold", size = 15),
        legend.position = "top") +
  ggtitle("DTU_Mix_Vs_GV") +
  labs(y = "-Log10(FDR)") +
  scale_y_continuous(limits = c(0,25))

ggsave(p1, filename = "p2.png", width = 2000, height = 2000, units = "px")
saveRDS(p1, "DTU_Mix_Damage_volcano.rds")

# summary the number of biotype of genes with isoform switching (also known as DTU) events (gDTU)
DTU_summary = dtu_result$isoformFeatures %>%
  subset(isoform_switch_q_value < 0.01 & abs(dIF) > 0.1) %>%
  group_by(iso_biotype) %>% tally()


# read DGE analysis results in
DEG_toptable = data.table::fread("./Mix sample result/DE result(Mix vs GV).csv")

# get the gene name of those differential expression genes
DEG = DEG_toptable %>% subset(DE != "No difference") %>% .$gene_name %>%
  na.omit()

# get the gene name of genes with isoform switching events
gDTU = dtu_result$isoformFeatures %>%
  subset(isoform_switch_q_value < 0.01 & abs(dIF) > 0.1)  %>%
  .$gene_name %>% unique() %>% na.omit()

# check the intersection between DEG and genes with DTU events
dplyr::intersect(DEG, gDTU)

# prepare gene symbol to entrezid conversion table
gene_info = AnnotationDbi::select(org.Ss.eg.db,
                                  keys = keys(org.Ss.eg.db),
                                  columns = c("ENTREZID", "SYMBOL"))

# use gDTU genes that doesn't ovaerlap with DEG to do GO analysis
gDTU_only = tibble(gene_name = dplyr::setdiff(gDTU, DEG)) %>%
  left_join(gene_info, by = c("gene_name" = "SYMBOL")) %>%
  .$ENTREZID %>% unique() %>% na.omit()

# all genes after QC remained for DGE analysis as background
un_genes = tibble(gene_name = DEG_toptable$gene_name) %>%
  left_join(gene_info, by = c("gene_name" = "SYMBOL")) %>%
  .$ENTREZID %>% unique() %>% na.omit()

# use those genes only identified by DTU analysis to do
# GO enrichment analysis
GO_gDTU = map(list(gDTU_only, gDTU_and_DEG, all_gDTU_DEG), ~ enrichGO(
  .x,
  org.Ss.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = un_genes,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T))

names(GO_gDTU) = c("only gDTU", "gDTU&DEG", "all")

# viusalize GO enrichment result
p2 = enrich.plot(GO_gDTU[[1]]@result, fraction.factor = 15,
            text.location = 10, str.wrap.size = 32, save = F,
            bar.font.size = 6,
            show.legend = F)
saveRDS(p2, "pig_DTU_GO.rds")

# extract sequences for external tools (CPC2, pfam, SignalIP and Iupred2A)
# The use of those external tools can follow tutorial of IsoformSwitchAnalyzeR
# package in this link
# https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
dtu_result = extractSequence(dtu_result, pathToOutput = "./files for external tools/",
              removeLongAAseq = T, alsoSplitFastaFile = T,
              alpha = 0.01, dIFcutoff = 0.1, writeToFile = T)

# add CPC2 analysis
dtu_result = analyzeCPC2(switchAnalyzeRlist = dtu_result,
                         pathToCPC2resultFile = "./output from external tools/cpc2/cpc2_integrated_result.txt",
                         codingCutoff = 0.725,
                         removeNoncodinORFs = T)

# add Pfam analysis
dtu_result = analyzePFAM(switchAnalyzeRlist = dtu_result,
                         pathToPFAMresultFile = "./output from external tools/pfam/pig_pfam.txt",
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
splicing_summary = extractSplicingSummary(dtu_result, returnResult = T)
splicing_summary = splicing_summary %>% mutate(splicingResult = str_extract(splicingResult, "less|more"))

# extract predictive consequences on protein expression
consequence_summary = extractConsequenceSummary(dtu_result, returnResult = T)
consequence_summary = consequence_summary %>%
  mutate(switchConsequence = str_extract(switchConsequence, "gain|loss|switch|insensitive|sensitive|Noncoding|coding|shorter|longer")) %>%
  mutate(switchConsequence = factor(switchConsequence, levels = c("gain", "loss", "switch","longer", "shorter",
                                                              "coding", "Noncoding", "sensitive", "insensitive")))

# visualize splicing events and functional consequence results
p3 = stackBar(splicing_summary, x = "AStype", y = "nrIsoWithConsequences", fill = "splicingResult",
         colors = c(color[2], color[1]), legend = F)

p4 = stackBar(consequence_summary, x = "featureCompared", y = "nrIsoWithConsequences", fill = "switchConsequence",
         colors = pal_npg("nrc")(10), legend = T)

# viusalize differential splicing events and fuctional consequence results
p5 = extractSplicingEnrichment(dtu_result, returnResult = F, returnSummary = F, countGenes = F)

p6 = extractConsequenceEnrichment(dtu_result, returnResult = F, countGenes = F)


saveRDS(p3, "pig_splicing_bar.rds")
saveRDS(p4, "pig_consequence_bar.rds")
saveRDS(p5, "pig_splicing_enrichment.rds")
saveRDS(p6, "pig_consequence_enrichment.rds")

# export switching, AS, switching consequence result
switching = dtu_result$isoformFeatures %>% dplyr::select(1:9, gene_value_1, gene_value_2,
                                                         iso_value_1, iso_value_2,
                                                         IF1, IF2, dIF, isoform_switch_q_value:IR)
AS = dtu_result$AlternativeSplicingAnalysis
switching_consequence = dtu_result$switchConsequence

writexl::write_xlsx(list(isoform_switch = switching, alternative_splicing = AS,
                         switching_consequence = switching_consequence), path = "./Mix sample result/DTU_analysis.xlsx")
