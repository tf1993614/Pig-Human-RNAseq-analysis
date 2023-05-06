###### customize functions for analysis and graphing ######

# load necessary packages
library(tidyverse)
library(ggrepel)
library(RcolorBrewer)


# create a function to draw volcano plot showing DGE analysis result
Volcano.plot = function(data,
                        colors = c("#00468B99","#ADB6B699", "#ED000099"),
                        title = NULL,
                        save = F,
                        log2fc.cutoff = 1,
                        pvalue.cutoff = 0.05,
                        format = ".png",
                        fig.wdith = 3000,
                        fig.height = 2000,
                        show.legend = "top"
){
  p = data %>% mutate(DE = case_when(adj.P.Val < pvalue.cutoff & logFC > log2fc.cutoff ~ "Up",
                                     adj.P.Val < pvalue.cutoff & logFC < -log2fc.cutoff ~ "Down",
                                     TRUE ~ "No difference")) %>%
    ggplot(aes(logFC, -log10(adj.P.Val), col = DE)) +
    geom_point() +
    geom_hline(yintercept = -log10(pvalue.cutoff), linetype = "dashed", size =0.8) +
    geom_vline(xintercept = c(-log2fc.cutoff, log2fc.cutoff), linetype = "dashed", size = 0.8) +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(text = element_text(colour = "black"),
          axis.title.x.bottom = element_text(size =13, face = "bold"),
          axis.text.x.bottom = element_text(face = "bold",colour = "black"),
          axis.title.y.left = element_text(size =13, face = "bold"),
          axis.text.y.left = element_text(face = "bold", colour = "black"),
          axis.line.x.bottom = element_line(size = 0.8),
          axis.line.y.left = element_line(size = 0.8),
          legend.title = element_blank(),
          legend.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 3,
                                        face = "bold", size = 15)) +
    xlim(c(-10,10))+
    ggtitle(title) +
    labs(y = "-Log10(FDR)", x = "Log2(Fold Change)")
  if(! is.null(show.legend)){
    p = p + theme(legend.position = show.legend)
  }
  if(save){
    ggsave(p, width = fig.wdith, height = fig.height, filename = str_c(title, format),
           units = "px")
  }
  return(p)
}



# create a function to find the number of up/down-regulated genes
# in the clusterProfile enrichment results
FindUpDown = function(data, up.list, down.list){
  geneID.ls = data$geneID %>% str_split("/")
  names(geneID.ls) = data$ID

  result = map_dfr(names(geneID.ls), function(x){
    dt = geneID.ls[[x]]
    numOfUp = sum(dt %in% up.list)
    numOfDown = sum(dt %in% down.list)
    output = tibble(ID = x, Up = numOfUp, Down = numOfDown)
  })

  return(result)
}



# create a function to find the total number of genes involved in the GO pathways
GO_summary = function(species, features = NULL){

  x = switch(species, human = "org.Hs.eg.db",
                   mouse = "org.Mm.eg.db",
                   pig = "org.Ss.eg.db")

  print(x)

  if(!require(x, character.only = T)){
    message("Downloading ", x, " package")
    BiocManager::install(x)
  }

  if(x == "org.Hs.eg.db"){
    library(org.Hs.eg.db)
    result = AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                                 columns = c("GOALL")) %>%
      subset(GOALL %in% features) %>%
      group_by(GOALL) %>% tally()
  }
  else if(x == "org.Mm.eg.db"){
    library(org.Mm.eg.db)
    result = AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db),
    columns = c("GOALL")) %>% subset(GOALL %in% features) %>%
      group_by(GOALL) %>% tally()
  }
  else if(x == "org.Ss.eg.db"){
    library(org.Ss.eg.db)
    result = AnnotationDbi::select(org.Ss.eg.db, keys = keys(org.Ss.eg.db),
    columns = c("GOALL")) %>% subset(GOALL %in% features) %>%
      group_by(GOALL) %>% tally()
  }
  return(result)
}



# create a function to visualize enrichment results done by
# clusterProfile package
enrich.plot = function(data,
                       colors = c("#00468B99","#ADB6B699", "#ED000099"),
                       type = "GO", # GO or KEGG
                       subtype = "UpDown", # All or UpDown
                       upGeneList = NULL, # required for when subtype is "All"
                       downGeneList = NULL, # required for when subtype is "All"
                       order.by = "p.adjust", # p.adjut or count
                       show.size = 5,
                       bar.width = NULL,
                       bar.font.size = 5,
                       point.size = 10,
                       line.size = 5,
                       show.legend = T,
                       fraction.factor = 3,
                       text.location = 3,
                       str.wrap.size = 10,
                       save = F,
                       format = ".png",
                       fig.height = 2000,
                       fig.width = 2000,
                       file.name = "enrichplot"){

  data = data %>% subset(p.adjust < 0.05) %>% arrange(p.adjust)


  if(type == "GO" & subtype == "UpDown"){

    if(order.by == "p.adjust"){
      new.data = data %>% arrange(p.adjust) %>% group_by(ONTOLOGY) %>%
        slice_head(n = show.size) %>%
        mutate(Description = factor(Description, levels = Description))
    }

    if(order.by == "Count"){
      new.data = data %>% arrange(p.adjust) %>% group_by(ONTOLOGY) %>%
        slice_head(n = show.size) %>% arrange(Count) %>%
        mutate(Description = factor(Description, levels = Description))
    }

    p = new.data %>% ggplot(aes(y = Description, fill = ONTOLOGY)) +
      geom_bar(aes(x = Count), stat = "identity", width = bar.width) +
      geom_path(aes(x = -log10(p.adjust) * fraction.factor, group = ONTOLOGY,
                    size = line.size), show.legend = F) +
      geom_point(aes(x = -log10(p.adjust) * fraction.factor,
                     size = point.size), show.legend = F) +
      scale_x_continuous(sec.axis = sec_axis( ~ . / fraction.factor,
                                              name = "-Log10(FDR)")) +
      scale_y_discrete(label = function(x) str_wrap(x, str.wrap.size))+
      geom_text(aes(label = str_c(str_remove(GeneRatio, "/\\d+")," (",
                                  str_remove(BgRatio, "/\\d+"), ")"),
                                  x = text.location, y = Description
      ), size = bar.font.size, show.legend = F) +
      theme_classic() +
      theme(axis.title.x.bottom = element_text(size =13, face = "bold", colour = "black"),
            axis.title.x.top = element_text(size =13, face = "bold", colour = "black"),
            axis.text.x.bottom = element_text(face = "bold", colour = "black"),
            axis.text.x.top = element_text(face = "bold", colour = "black"),
            axis.text.y.left = element_text(face = "bold", colour = "black"),
            axis.line.x.bottom = element_line(size = 0.8),
            axis.line.y.left = element_line(size = 0.8),
            axis.line.x.top = element_line(size = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(face = "bold", size =11),
            axis.title.y.left = element_blank(),
            panel.grid.major.x = element_line(linewidth = 0.8, linetype = "dashed"),
      ) +
      labs(x = "Number of DEGs involved in the pathway")
    if(! show.legend){
      p = p + theme(legend.position = "none")
    }
  }
  else if( (type == "GO" | type == "KEGG") & subtype == "All"){
    if(order.by == "p.adjust"){
      new.data = data %>% arrange(p.adjust) %>% head(n = show.size) %>%
        mutate(Description = factor(Description, levels = Description))
    }

    if(order.by == "Count"){
      new.data = data %>% arrange(desc(Count)) %>% head(show.size) %>%
        mutate(Description = factor(Description, levels = Description))
    }

    summary_table = FindUpDown(new.data, up.list = upGeneList,
                               down.list = downGeneList)

    new.data = new.data %>% left_join(summary_table, by = "ID") %>%
      pivot_longer(c("Up", "Down"), names_to = "UpDown", values_to = "N")

    p = new.data %>% ggplot(aes(y = Description)) +
      geom_bar(aes(x = N, fill = UpDown), stat = "identity") +
      geom_path(aes(x = -log10(p.adjust) * fraction.factor, group = 1,
                    size = line.size), show.legend = F) +
      geom_point(aes(x = -log10(p.adjust) * fraction.factor,
                     size = point.size), show.legend = F) +
      scale_x_continuous(sec.axis = sec_axis( ~ . / fraction.factor,
                                              name = "-Log10(FDR)")) +
      scale_y_discrete(label = function(x) str_wrap(x, str.wrap.size))+
      geom_text(aes(label = str_remove(GeneRatio, "/\\d+"), x = text.location, y = Description,
      ), show.legend = F) +
      scale_fill_manual(values = colors)+
      theme_classic() +
      theme(axis.title.x.bottom = element_text(size =13, face = "bold", colour = "black"),
            axis.title.x.top = element_text(size =13, face = "bold", colour = "black"),
            axis.text.x.bottom = element_text(face = "bold", colour = "black"),
            axis.text.x.top = element_text(face = "bold", colour = "black"),
            axis.text.y.left = element_text(face = "bold", colour = "black"),
            axis.line.x.bottom = element_line(size = 0.8),
            axis.line.y.left = element_line(size = 0.8),
            axis.line.x.top = element_line(size = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(face = "bold", size =11),
            axis.title.y.left = element_blank(),
            panel.grid.major.x = element_line(linewidth = 0.8, linetype = "dashed"),
      ) +
      labs(x = "The number of DE genes involved in the pathway")
  }
  else if(type == "KEGG" & subtype == "UpDown"){
    p = data %>%
      mutate(Description = factor(Description, levels = Description)) %>%
      ggplot(aes(y = Description)) +
      geom_bar(aes(x = Count), stat = "identity" ,fill = sample(colors, 1)) +
      geom_path(aes(x = -log10(p.adjust) * fraction.factor,
                    size = line.size, group = 1), show.legend = F, colour = "black") +
      geom_point(aes(x = -log10(p.adjust) * fraction.factor,
                     size = point.size), show.legend = F) +
      scale_x_continuous(sec.axis = sec_axis( ~ . / fraction.factor,
                                              name = "-Log10(FDR)")) +
      scale_y_discrete(label = function(x) str_wrap(x, str.wrap.size))+
      geom_text(aes(label = str_remove(GeneRatio, "/\\d+"), x = text.location, y = Description,
      ), show.legend = F) +
      theme_classic() +
      theme(axis.title.x.bottom = element_text(size =13, face = "bold", colour = "black"),
            axis.title.x.top = element_text(size =13, face = "bold", colour = "black"),
            axis.text.x.bottom = element_text(face = "bold", colour = "black"),
            axis.text.x.top = element_text(face = "bold", colour = "black"),
            axis.text.y.left = element_text(face = "bold", colour = "black"),
            axis.line.x.bottom = element_line(size = 0.8),
            axis.line.y.left = element_line(size = 0.8),
            axis.line.x.top = element_line(size = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(face = "bold", size =11),
            axis.title.y.left = element_blank(),
            panel.grid.major.x = element_line(linewidth = 0.8, linetype = "dashed"),
      ) +
      labs(x = "The number of DE genes involved in the pathway")
  }



  if(!is.null(colors) & subtype != "All"){
    p = p + scale_fill_manual(values = colors)
  }

  if(save){
    ggsave(p, width = fig.width, height = fig.height, units = "px",
           filename = str_c(file.name, format))
  }

  return(p)

}



# create a function to do PCA analysis and then graphing
PCA_analysis = function(DGElistObject,
                        transpose = TRUE,
                        features = NULL,
                        scale = FALSE,
                        plot = FALSE,
                        group = NULL,
                        label = FALSE,
                        point.size = 4,
                        text.size = 6,
                        reName_samples = NULL,
                        colors = c("#E64B35FF","#4DBBD5FF", "#00A087FF"),
                        save = FALSE,
                        format = ".png",
                        width = 2000,
                        height = 2000,
                        fileName = "pca",
                        show.legend = TRUE){

  if(is.null(features)){
    features = rownames(DGElistObject)
  }

  if(transpose){
    pca_data = DGElistObject  %>% cpm(log = T) %>%
      subset(rownames(.) %in% features) %>%
      t() %>% prcomp(scale. = scale, center = T)
  }
  else{
    pca_data = DGElistObject %>% cpm(log = T) %>%
      subset(rownames(.) %in% features) %>%
      prcomp(scale. = scale, center = T)
  }

  pca_summary = summary(pca_data)$importance %>% as.data.frame()


  samples_attr = DGElistObject$samples %>% rownames_to_column("sampleName") %>%
    dplyr::select(sampleName, group)


  if(transpose){
    pca_data = pca_data$x %>% as.data.frame() %>%
      rownames_to_column("sampleName") %>%
      left_join(samples_attr, by = "sampleName")

    if(! is.null(reName_samples)){
      pca_data[["sampleName"]] = reName_samples
    }

    pca_f =  pca_data %>%
      ggplot(aes(PC1,PC2, color = group)) + geom_point(size = point.size) +
      labs(x = str_c("PC1 (", as.character(scales::percent(pca_summary[2,1])), ")"),
           y = str_c("PC2 (", as.character(scales::percent(pca_summary[2,2])), ")")) +
      scale_color_manual(values = colors) +
      theme_classic() +
      theme(axis.title.x.bottom = element_text(size =13, face = "bold"),
            axis.text.x.bottom = element_text(face = "bold", colour = "black"),
            axis.title.y.left = element_text(size =13, face = "bold"),
            axis.text.y.left = element_text(face = "bold", colour = "black"),
            axis.line.x.bottom = element_line(size = 0.8),
            axis.line.y.left = element_line(size = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, vjust = 3,
                                      face = "bold", size = 15),
            panel.grid.major = element_line(colour = "grey", linewidth = 0.5,
                                            linetype = 2),
            panel.grid.minor = element_line(colour = "grey", linewidth = 0.5,
                                            linetype = 2))
    if(! show.legend){
      pca_f = pca_f + theme(legend.position = "none")
    }
  }
  else{
    pca_data = pca_data$x %>% as.data.frame() %>% rownames_to_column("sampleName")

    pca_f =  pca_data %>%
      ggplot(aes(PC1,PC2)) + geom_point() +
      labs(x = str_c("PC1 (", as.character(scales::percent(pca_summary[2,1])), ")"),
           y = str_c("PC2 (", as.character(scales::percent(pca_summary[2,2])), ")"))
  }

  if(label){
    pca_f = pca_f + geom_label_repel(aes(label = sampleName), max.overlaps = 20,
                     size = text.size)
  }

  if(plot){
    return(pca_f)
  }
  else{
    return(pca_data)
  }

  if(save){
    ggsave(pca_f, width = width, height = height, units = "px",
           filename = str_c(fileName, format))
  }
}



# create a function based on fisher exact test to checck
# whether there is an overlap of DEGs on specific GO pathway
# between two datasets
pathway_overlap = function(dataset1, # enrichment result dataframe produced by clusterProfile package
                           geneSet1, # a named list of genes involved in GO pathway
                           species1, # species name
                           dataset2,
                           geneSet2,
                           species2){

  vector1 = dataset1[["geneID"]] %>% str_split("/")
  names(vector1) = dataset1[["ID"]]

  vector2 = dataset2[["geneID"]] %>% str_split("/")
  names(vector2) = dataset2[["ID"]]

  species_1_GO_db = switch(species1, pig = org.Ss.eg.db,
                           human = org.Hs.eg.db,
                           mouse = org.Mm.eg.db)

  species_2_GO_db = switch(species2, pig = org.Ss.eg.db,
                           human = org.Hs.eg.db,
                           mouse = org.Mm.eg.db)

  species_1_GO_all = geneSet1 %>% subset(names(.) %in% dataset1[["ID"]])
  species_2_GO_all = geneSet2 %>% subset(names(.) %in% dataset2[["ID"]])

  species_1_GO_all = map_dfr(names(species_1_GO_all), function(x){
    GOID = x
    SYMBOL = mapIds(species_1_GO_db, keys = species_1_GO_all[[x]],
                    column = "SYMBOL", keytype = "ENTREZID")

    data = data.frame(GOID = GOID, SYMBOL = SYMBOL)
  })

  species_2_GO_all = map_dfr(names(species_2_GO_all), function(x){
    GOID = x
    SYMBOL = mapIds(species_2_GO_db, keys = species_2_GO_all[[x]],
                    column = "SYMBOL", keytype = "ENTREZID")

    data = data.frame(GOID = GOID, SYMBOL = SYMBOL)
  })

  IDs = unique(c(dataset1[["ID"]], dataset2[["ID"]]))

  statistic_table = map_dfr(IDs, function(x){


    gene1 = species_1_GO_all %>% dplyr::filter(GOID == x) %>% .$SYMBOL %>% str_to_upper()

    data = map_dfr(IDs, function(ID){
      gene2 = species_2_GO_all %>% dplyr::filter(GOID == ID) %>% .$SYMBOL %>% str_to_upper()
      geneAll = union(gene1, gene2)
      InA = length(setdiff(str_to_upper(vector1[[x]]), str_to_upper(vector2[[ID]])))
      InB = length(setdiff(str_to_upper(vector2[[x]]), str_to_upper(vector1[[ID]])))
      intersect_size = length(intersect(str_to_upper(vector1[[x]]), str_to_upper(vector2[[ID]])))
      notAnotB = length(geneAll) - length(union(InA, InB))
      data = data.frame(GOID_1 = x, GOID_2 = ID, InA = InA, InB = InB,
                        intersect = intersect_size, notAnotB = notAnotB)
    })

    return(data)
  })


  fisher_res = map_dfr(IDs, function(x){

      data = statistic_table %>% subset(GOID_1 == x)

      data2 = map(data[["GOID_2"]],function(y){

        data3 = data %>% subset(GOID_2 == y)
        mat = matrix(c(data3[["notAnotB"]], data3[["InB"]],
                     data3[["InA"]],
                     data3[["intersect"]]),
                   nrow = 2)

       fisher = fisher.test(mat, alternative = "greater") %>% broom::tidy() %>%
        mutate(GOID_1 = x, GOID_2 = y) %>% dplyr::select(GOID_1, GOID_2, everything())
      })


       return(data2)
  })

  return(list(statistic_table, fisher_res))

}



# create stack bar plot function
stackBar = function(data, x , y, fill = NULL, colors = NULL, legend = T){
  x = sym(x)
  y = sym(y)
  fill = sym(fill)

  p = data %>% ggplot(aes(x = !! x, y = !! y, fill = !! fill)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(
      axis.text.x.bottom = element_text(face = "bold", colour = "black"),
      axis.title.y.left = element_text(size =13, face = "bold"),
      axis.text.y.left = element_text(face = "bold", colour = "black"),
      axis.line.x.bottom = element_line(size = 0.8),
      axis.line.y.left = element_line(size = 0.8),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, vjust = 3,
                                face = "bold", size = 15),
      panel.grid.major.y  = element_line(colour = "grey", linewidth = 0.5,
                                         linetype = 2),
      axis.title.x.bottom = element_blank()) +
    labs(y = "Number of significant isoforms")

  if(!is.null(colors)){
    p = p + scale_fill_manual(values = colors, name = "Type")
  }

  if(! legend){
    p = p + theme(legend.position = "none")
  }

  return(p)
}


# create a bar plot function
Barplot = function(data, x = "marker_gene", y = "n",
                          color = NULL,
                         y_title = "Proportion",
                         RcolorBrewer.palette = "Set1"
                         ){
  x = sym(x)
  y = sym(y)

  data = data %>% arrange(desc(!! y)) %>% mutate(!! x := factor(!! x, levels = !! x))

  p = data %>% ggplot(aes(!! x, !! y, fill = !! x)) +
    geom_bar(stat = "identity") +
    geom_text(aes(x = !! x , y = !! y + 0.05, label = scales::percent(!! y, accuracy = 0.1))) +
    theme_classic()+
    bar_theme() + labs(y = y_title) +
    scale_fill_brewer(palette = RcolorBrewer.palette) +
    scale_y_continuous(labels =  scales::percent, limits = c(0,1)) +
    theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none")

  return(p)
}
