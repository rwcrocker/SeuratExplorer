
### Ligand-Receptor Pipeline  ###

##  --------------  ##
##  Technical Fxns  ##

#Function to read the csv "pairmates" column as a vector
csv_to_vector = function(string){
  vector = unlist(strsplit(string, split=", "))
  return(vector)
}

##  -------------  ##
##  Analysis Fxns  ##

# Determine effect of condition for each cluster
Find_Condition_DEGs=function(seurat_obj, condition, reference.cond, experimental.cond, padj.cut=0.05){
  require("Seurat", "tidyverse", "ggplot2", "cowplot")

  DEGs = list()
  for (clust in unique(seurat_obj@active.ident)){
    DEGs[[clust]] = FindMarkers(object = seurat_obj, ident.1 = experimental.cond, ident.2 = reference.cond, group.by = condition, subset.ident = clust, logfc.threshold = 0.1)
    DEGs[[clust]]$gene = rownames(DEGs[[clust]])
    DEGs[[clust]] = DEGs[[clust]]%>%
      filter(p_val_adj < padj.cut)%>%
      mutate(cluster = clust)
  }
  DEGs = bind_rows(DEGs)
  return(DEGs)
}


#Subset markers to LR members with CellTalk database and generate a LR table
Analyze_LR = function(marker_df, db){
  req_packages = c("Seurat", "tidyverse", "ggplot2", "cowplot")
  lapply(req_packages, require, character.only = TRUE)
  
  #Annotate Ls and Rs in marker set
  marker_df$is_ligand=NA
  marker_df$is_receptor=NA
  for (num in 1:length(marker_df$gene)){
    marker_gene = marker_df$gene[num]
    marker_df$is_ligand[num] = marker_gene %in% db$ligand_gene_symbol
    marker_df$is_receptor[num] = marker_gene %in% db$receptor_gene_symbol
  }
  
  marker_df = marker_df %>% filter(is_ligand == TRUE | is_receptor == TRUE) # subset to LRs only
  
  #Find all potential and real pairmates for each gene
  marker_df$possible_pairmates = NA
  marker_df$pairmates = NA
  for (num in 1:length(marker_df$gene)){
    marker_gene = marker_df$gene[num]
    subset_db = db %>% filter(ligand_gene_symbol == marker_gene | receptor_gene_symbol == marker_gene)
    
    potential_pairmates = unique(c(subset_db$ligand_gene_symbol, subset_db$receptor_gene_symbol))
    potential_pairmates = potential_pairmates[which(!(potential_pairmates == marker_gene))]
    potential_pairmates_string = paste(potential_pairmates, collapse = ", ")
    
    marker_df$possible_pairmates[num] = potential_pairmates_string
    
    pairmates = potential_pairmates[which(potential_pairmates %in% marker_df$gene)]
    pairmates_string = NA
    if(length(pairmates) > 0){
      pairmates_string = paste(pairmates, collapse = ", ")
    }
    marker_df$pairmates[num] = pairmates_string
  }
  
  #subset to those w/ >= one real pair
  marker_df = marker_df %>% filter(!(is.na(pairmates)))
  
  #generate another df with only ligands
  lig_markers = marker_df%>%filter(is_ligand == TRUE)
  
  #Generate LR interaction_df from marker_df
  interaction = c()
  ligand = c()
  ligand_cluster=c()
  ligand_marker_log2fc=c()
  ligand_marker_fdr=c()
  receptor = c()
  receptor_cluster=c()
  receptor_marker_log2fc=c()
  receptor_marker_fdr=c()
  for (num in 1:length(lig_markers$gene)){
    pairmates = csv_to_vector(lig_markers$pairmates[num])
    for (pairmate in pairmates){
      pairmate_clusters = as.vector(unique((marker_df%>%filter(gene == pairmate))$cluster))
      for (pairmate_cluster in pairmate_clusters){
        ligand = c(ligand, lig_markers$gene[num])
        ligand_cluster = c(ligand_cluster, as.vector(lig_markers$cluster[num]))
        ligand_marker_fdr = c(ligand_marker_fdr, lig_markers$p_val_adj[num])
        ligand_marker_log2fc = c(ligand_marker_log2fc, lig_markers$avg_log2FC[num])
        
        receptor = c(receptor, pairmate)
        receptor_cluster = c(receptor_cluster, pairmate_cluster)
        receptor_marker_fdr = c(receptor_marker_fdr, marker_df$p_val_adj[which(marker_df$gene == pairmate & marker_df$cluster == pairmate_cluster)])
        receptor_marker_log2fc = c(receptor_marker_log2fc, marker_df$avg_log2FC[which(marker_df$gene == pairmate & marker_df$cluster == pairmate_cluster)])
        
        interaction = c(interaction, paste0(lig_markers$gene[num], "_",lig_markers$cluster[num], "_", pairmate, "_", pairmate_cluster))
      }
    }
  }
  
  LR_table = data.frame(interaction, ligand, ligand_cluster, ligand_marker_log2fc, ligand_marker_fdr,
                        receptor, receptor_cluster, receptor_marker_log2fc, receptor_marker_fdr)
  
  return(LR_table)
  
}

Crossreference_LR = function(LR.table, DEG.table, db){
  req_packages = c("Seurat", "tidyverse", "ggplot2", "cowplot", "plyr")
  lapply(req_packages, require, character.only = TRUE)
  
  LR.table$receptor_is_condition_DE=NA
  LR.table$ligand_is_condition_DE=NA
  
  LR.table$ligand_condition_log2fc =NA
  LR.table$receptor_condition_log2fc=NA
  LR.table$ligand_condition_fdr=NA
  LR.table$receptor_condition_fdr=NA
  
  for (num in 1:length(LR.table$interaction)){
    Lig = LR.table$ligand[num]
    Lig_clust = LR.table$ligand_cluster[num]
    Rec = LR.table$receptor[num]
    Rec_clust = LR.table$receptor_cluster[num]
    
    Lig_DEG = DEG.table %>% filter(gene == Lig & cluster == Lig_clust)
    LR.table$ligand_is_condition_DE[num] = as.logical(dim(DEG.table %>% filter(gene == Lig & cluster == Lig_clust))[1])
    if (LR.table$ligand_is_condition_DE[num]){
      LR.table$ligand_condition_log2fc[num] = Lig_DEG$avg_log2FC
      LR.table$ligand_condition_fdr[num] = Lig_DEG$p_val_adj
    }
    
    Rec_DEG = DEG.table %>% filter(gene == Rec & cluster == Rec_clust)
    LR.table$receptor_is_condition_DE[num] = as.logical(dim(DEG.table %>% filter(gene == Rec & cluster == Rec_clust))[1])
    if(LR.table$receptor_is_condition_DE[num]){
      LR.table$receptor_condition_log2fc[num] = Rec_DEG$avg_log2FC
      LR.table$receptor_condition_fdr[num] = Rec_DEG$p_val_adj
    }
  }
  
  LR.table = LR.table %>% filter(ligand_is_condition_DE == TRUE | receptor_is_condition_DE == TRUE)
  
  
  #Check to make sure there are not DEG-based LR pairs that are not signature genes...
  #Use the same function on DEGs instead of signatures
  #Rename stats to be consistent with the LR.table
  backcross = Analyze_LR(DEG.table, db = db)
  
  if(!(all(backcross$interaction %in% LR.table$interaction))){
    backcross = backcross %>%
      filter(!(interaction %in% LR.table$interaction))%>%
      mutate(ligand_is_condition_DE = TRUE)%>%
      mutate(receptor_is_condition_DE = TRUE)%>%
      dplyr::rename(ligand_condition_log2fc = ligand_marker_log2fc)%>%
      dplyr::rename(ligand_condition_fdr = ligand_marker_fdr)%>%
      dplyr::rename(receptor_condition_log2fc = receptor_marker_log2fc)%>%
      dplyr::rename(receptor_condition_fdr = receptor_marker_fdr)
    LR.table = rbind.fill(LR.table, backcross)
  }
  else{
    print("All DE LR pairs are also signature genes",)
  }

  LR.table = LR.table %>%
    mutate(both_are_condition_DE = ( receptor_is_condition_DE & ligand_is_condition_DE ) )

  return(LR.table)
}
