library(liana)
library(CellChat)
library(dplyr)
library(Seurat)
library(magrittr)
library(openxlsx)
library(SCpubr)

#Was actually Spleen and Lymph together
#Receptors for T cells
#ligands for DCs

###
#This produces outputs for Spleen data. To analyse Lymph Node data, change `Spleen` on lines 177 and 179 to `LN`
###

##Load Data ##
seu_T <- readRDS('seu_Tcell.Rds')
seu_myeloid <- readRDS('seu_DC.Rds')

##adjust annotations ##
seu_myeloid$annotation<-as.character(seu_myeloid$SCT2_snn_res.0.2)
seu_myeloid$annotation[seu_myeloid$annotation == "3"]<-"migDCs"
seu_myeloid$annotation[seu_myeloid$annotation == "1"]<-"moDCs"
seu_myeloid$annotation[seu_myeloid$annotation == "2"]<-"cDC1"
seu_myeloid$annotation[seu_myeloid$annotation %in% c("0","5")]<-"cDC2"
seu_myeloid$annotation[seu_myeloid$annotation == "4"]<-"pDC"
seu_myeloid$annotation[seu_myeloid$annotation == "6"]<-"unknown"
Idents(seu_myeloid)<-"annotation"
seu_myeloid<-seu_myeloid[,!seu_myeloid$annotation == "unknown"]
do_DimPlot(seu_myeloid)

seu_T$annotation<-seu_T$SCT2_snn_res.0.2
seu_T$annotation<-ifelse(seu_T$annotation == "0", paste0("TPEX"),paste0("TEX"))
do_DimPlot(seu_T, group.by = "annotation", reduction = "umap3")
Idents(seu_T)<-"annotation"



## Prepare database##

#Use cellchat database of ligands/receptors
cellchatDB <- CellChat::CellChatDB.mouse
ligands <- select(cellchatDB$interaction, ligand)
receptors <- select(cellchatDB$interaction, receptor)
receptors_complex <- cellchatDB$complex
receptors <- mutate(receptors, is_complex = receptor %in% rownames(receptors_complex)) %>% tibble::rownames_to_column(var = 'interaction')
ligands <- mutate(ligands, is_complex = ligand %in% rownames(receptors_complex)) %>% tibble::rownames_to_column(var = 'interaction')

#set up expression tables
##simple ligand simple receptor
simple_ligands_idx <- which(ligands$is_complex == F)
simple_receptors_idx <- which(receptors$is_complex == F)
simple_simple_idx <- intersect(simple_ligands_idx, simple_receptors_idx)
simple_simple_df <- cbind(ligands, receptors)[simple_simple_idx,c(2,5,1)]


##complex expressions (receptor in complex or ligand in complex)
##complex receptors
complex_receptor_vector <- receptors_complex[filter(receptors, is_complex ==T)$receptor,] %>% distinct() %>% tibble::rownames_to_column(var = 'name') %>% tidyr::pivot_longer(cols = !name, values_to = 'gene', names_to = 'subunit')
complex_receptor_vector <- filter(complex_receptor_vector, gene != "")
#match order for receptors
idx_complex_receptor_to_ligand <- which(cellchatDB$interaction$receptor %in% complex_receptor_vector$name)
complex_interaction_ligands <- cellchatDB$interaction[idx_complex_receptor_to_ligand,'ligand', drop = F]
#is the ligand complex too? (ie. both receptor and ligand are complex)
complex_ligands_for_complex_receptors <- receptors_complex[complex_interaction_ligands$ligand,] %>% na.omit() %>% distinct()

simple_ligands <- filter(complex_interaction_ligands, !ligand %in% rownames(complex_ligands_for_complex_receptors))

#map complex receptor subunits with simple ligands
complex_receptor_simple_ligand <- cellchatDB$interaction %>% filter(receptor %in% complex_receptor_vector$name & ligand %in% simple_ligands$ligand)

complex_receptor_simple_ligand_pairings <- apply( complex_receptor_simple_ligand, 1, function(x) {
  receptor = x[4]
  ligand = x[3]
  
  subunits_receptor = filter(complex_receptor_vector, name == receptor) %>% select(gene) %>% unlist()
  return(data.frame(ligand = ligand, receptor = subunits_receptor, interaction = x[1]))
})
complex_receptor_simple_ligand_df <- do.call(rbind, complex_receptor_simple_ligand_pairings) %>% distinct


##complex ligands
complex_ligand_vector <- receptors_complex[filter(ligands, is_complex ==T)$ligand,] %>% distinct() %>% tibble::rownames_to_column(var = 'name') %>% tidyr::pivot_longer(cols = !name, values_to = 'gene', names_to = 'subunit')
complex_ligand_vector <- filter(complex_ligand_vector, gene != "")
#match order for ligands
idx_complex_ligand_to_receptor <- which(cellchatDB$interaction$ligand %in% complex_ligand_vector$name)
complex_interaction_receptors <- cellchatDB$interaction[idx_complex_ligand_to_receptor,'receptor', drop = F]
#is the receptor complex too?
complex_receptors_for_complex_ligands <- receptors_complex[complex_interaction_receptors$receptor,] %>% na.omit() %>% distinct()
simple_receptors <- filter(complex_interaction_receptors, !receptor %in% rownames(complex_receptors_for_complex_ligands))
#map complex ligand subunits with simple receptors
complex_ligand_simple_receptor <- cellchatDB$interaction %>% filter(ligand %in% complex_ligand_vector$name & receptor %in% simple_receptors$receptor)

complex_ligand_simple_receptor_pairings <- apply( complex_ligand_simple_receptor, 1, function(x) {
  receptor = x[4]
  ligand = x[3]
  
  subunits_ligand = filter(complex_ligand_vector, name == ligand) %>% select(gene) %>% unlist()
  return(data.frame(ligand = subunits_ligand, receptor = receptor, interaction = x[1]))
})

complex_ligand_simple_receptor_df <- do.call(rbind, complex_ligand_simple_receptor_pairings) %>% distinct



#complex ligand and complex receptor
complex_to_complex_interactions <- cellchatDB$interaction %>% filter(ligand %in% rownames(complex_ligands_for_complex_receptors) & receptor %in% rownames(complex_receptors_for_complex_ligands))

complex_complex_list <- apply(complex_to_complex_interactions, 1, function(x) {
  ligand = x[3]
  receptor = x[4]
  subunits_receptors <- receptors_complex[ligand,] %>% {.[which(. != "")]} %>% unlist() 
  subunits_ligands <- receptors_complex[receptor,] %>% {.[which(. != "")]} %>% unlist()
  out_grid <- expand.grid(subunits_ligands, subunits_receptors)
  colnames(out_grid) <- c('ligand', 'receptor')
  out_grid$interaction = x[1]
  return(out_grid)
}, simplify = F)

complex_complex_df <- do.call(rbind, complex_complex_list) %>% distinct 

#put all together. all of these need to have interaction name in a column
all_interactions <- rbind(simple_simple_df, complex_ligand_simple_receptor_df, complex_receptor_simple_ligand_df, complex_complex_df) %>% distinct()



##Get expressions from Seurat object ##
##Reminder - ligands from  Myeloid, receptors from T

#ligand expressions
ligand_expr <- Seurat::AverageExpression(seu_myeloid, features = unique(all_interactions$ligand), group.by = 'annotation', assays = 'RNA', use.scale = T)

#receptor expressions
receptor_expr <- Seurat::AverageExpression(seu_T, features = unique(all_interactions$receptor), group.by = 'annotation', assays = 'RNA', use.scale = T)

#put expression in right orders
ligand_expr_tables <- left_join(all_interactions  %>% select(c(interaction, ligand)), as.data.frame(ligand_expr$RNA) %>% tibble::rownames_to_column('gene'), join_by(ligand == gene))
receptor_expr_tables <- left_join(all_interactions  %>% select(c(interaction, receptor)), as.data.frame(receptor_expr$RNA) %>% tibble::rownames_to_column('gene'), join_by(receptor == gene))

#percent cells expressing tables
##ligands

  clusters <- levels(as.factor(seu_myeloid$annotation))
  subset_cluster_ligands <- lapply(clusters, function(cluster) {
    subs <- subset(seu_myeloid, annotation == cluster, features = all_interactions$ligand)
    return(subs)
  })
  names(subset_cluster_ligands) <- clusters


pct_express_ligands <-  lapply(subset_cluster_ligands, function(subs) {
    n_cells <- dim(subs)[2]
    n_expressing <- rowSums(as.matrix(subs@assays$RNA@counts > 0))
    pct.expressing <- n_expressing / n_cells
    return(pct.expressing)
  }) %>% data.frame()

##receptors

  clusters <- levels(as.factor(seu_T$annotation))
  subset_cluster_receptors <- lapply(clusters, function(cluster) {
    subs <- subset(seu_T, annotation == cluster, features = all_interactions$receptor)
    return(subs)
  }) 
  names(subset_cluster_receptors) <- clusters


pct_express_receptors <- lapply(subset_cluster_receptors, function(subs) {
    n_cells <- dim(subs)[2]
    n_expressing <- rowSums(as.matrix(subs@assays$RNA@counts > 0))
    pct.expressing <- n_expressing / n_cells
    return(pct.expressing)
  }) %>% data.frame()

#put percent cells expressing in right orders
ligand_pct_tables <- left_join(select(all_interactions, 'ligand'), as.data.frame(pct_express_ligands) %>% tibble::rownames_to_column('gene'), join_by(ligand == gene))
receptor_pct_tables <- left_join(select(all_interactions, 'receptor'), as.data.frame(pct_express_receptors) %>% tibble::rownames_to_column('gene'), join_by(receptor == gene))

##Plot interactions ##
#plot expression correlations first
receptor_expr_for_bind <- receptor_expr_tables %>% as.data.frame() %>% rename_with(
  ~ paste0("receptor.", .x, recycle0 = TRUE) , !(receptor | interaction )) 
ligand_expr_for_bind <- ligand_expr_tables %>% as.data.frame() %>% rename_with(
  ~ paste0("ligand.", .x, recycle0 = TRUE), !(ligand | interaction)) 

all_expr_plot <- cbind(receptor_expr_for_bind, ligand_expr_for_bind) %>% tidyr::pivot_longer(cols = starts_with('receptor.'), names_to = 'receptor_cluster', values_to = 'receptor_expr') %>% tidyr::pivot_longer(cols = starts_with('ligand.'), names_to = 'ligand_cluster', values_to = 'ligand_expr')

scatter_plots <- ggplot(all_expr_plot, aes(x = ligand_expr, y = receptor_expr)) + geom_point() + facet_grid(rows = vars(receptor_cluster), cols = vars(ligand_cluster)) + coord_cartesian(xlim = c(1,10), ylim = c(1,10)) + geom_rug(col=rgb(.5,0,0,alpha=.2))
scatter_plots

#plot the average pct of cells expressing ligand/receptor
#Using spleen here. change to LN to do Lymph Node
receptor_pct_for_bind <- receptor_pct_tables %>% as.data.frame() %>% rename_with(
  ~ paste0("receptor.", .x, recycle0 = TRUE) , !receptor) 
ligand_pct_for_bind <- ligand_pct_tables %>% as.data.frame() %>% rename_with(
  ~ paste0("ligand.", .x, recycle0 = TRUE), !ligand) 

all_pct_plot <- cbind(receptor_pct_for_bind, ligand_pct_for_bind) %>% tidyr::pivot_longer(cols = starts_with('receptor.'), names_to = 'receptor_cluster', values_to = 'receptor_pct') %>% tidyr::pivot_longer(cols = starts_with('ligand.'), names_to = 'ligand_cluster', values_to = 'ligand_pct') %>% distinct()
#distinct() used because data is longer than it needs to be. features from different lists may be input twice
pct_plots <- ggplot(all_pct_plot, aes(x = ligand_pct, y = receptor_pct)) + geom_point() + facet_grid(rows = vars(receptor_cluster), cols = vars(ligand_cluster)) + coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + geom_rug(col=rgb(.5,0,0,alpha=.2))
pct_plots

##Identify putative interactions ##
#Set thresholds for expression quantity and percentage of cells expressing ligand/receptor
pct_thresh = 0.1
expr_thresh = 0.5

#identify interactions passing both thresholds
expr_pass_thresh <- filter(all_expr_plot, receptor_expr > expr_thresh & ligand_expr > expr_thresh)
pct_pass_thresh <- filter(all_pct_plot, receptor_pct > pct_thresh & ligand_pct > pct_thresh)

pass_both <- inner_join(x = expr_pass_thresh, y = pct_pass_thresh)
pass_both$receptor_cluster<-gsub("receptor.", "",pass_both$receptor_cluster)
pass_both$ligand_cluster<-gsub("ligand.", "",pass_both$ligand_cluster)

pre_pass<-inner_join(all_expr_plot,all_pct_plot) %>% unique()
pre_pass$receptor_cluster<-gsub("receptor.", "",pre_pass$receptor_cluster)
pre_pass$ligand_cluster<-gsub("ligand.", "",pre_pass$ligand_cluster)

expr_plots_pass_both_thresh <- ggplot(pass_both, aes(x = ligand_expr, y = receptor_expr)) + geom_point() + facet_grid(rows = vars(receptor_cluster), cols = vars(ligand_cluster)) + coord_cartesian(xlim = c(0,10), ylim = c(0,10))
expr_plots_pass_both_thresh

pass_both_plusZ <- filter(pre_pass, receptor_expr > expr_thresh & ligand_expr > expr_thresh & receptor_pct > pct_thresh & ligand_pct > pct_thresh)
pass_both_plusZ<-rbind(pass_both,pre_pass[pre_pass$receptor == "Tgfbr2" & pre_pass$ligand == "Tgfb1",]) %>% unique()
pass_both_plusZ<-rbind(pass_both_plusZ,pre_pass[pre_pass$receptor == "Ccr5" & pre_pass$ligand == "Ccl5",]) %>% unique()
pass_both_plusZ<-rbind(pass_both_plusZ,pre_pass[pre_pass$receptor == "Sell" & pre_pass$ligand == "Selplg",]) %>% unique()

pass_both_plusZ<-pass_both_plusZ %>% group_by(interaction) %>% mutate(ligand_Z = (ligand_expr - mean(ligand_expr))/sd(ligand_expr),receptor_Z = (receptor_expr - mean(receptor_expr))/sd(receptor_expr)) %>% replace(is.na(.),0)




#bubble plot
df.plot <- pass_both_plusZ  %>% filter(!grepl("H2-",ligand)) 
MHC_pb <- pass_both_plusZ[grepl("H2-",pass_both_plusZ$ligand),] %>% group_by(ligand_cluster,receptor_cluster) %>% 
  summarise(ligand_Z = mean(ligand_Z),receptor_Z = mean(receptor_Z), ligand_pct = mean(ligand_pct), receptor_pct = mean(receptor_pct)) %>% 
  mutate(ligand = "MHC", receptor = "MHC_receptor") %>% mutate(interaction = "MHC_MHC-receptor")

df.plotfull<-full_join(df.plot,MHC_pb) %>% group_by(ligand_cluster, receptor_cluster) 

Tcell_dotZ2<-df.plot %>%
  ggplot() +
  aes(x=receptor_cluster,y=paste(ligand,"->",receptor), size = ligand_Z, color = receptor_Z) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_point() +
  facet_grid(.~ligand_cluster, scales = "free_x") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, size=10, vjust = .4,hjust = 1), axis.text.y = element_text(size=10)) +
  labs(x="Target cells", y="Ligand/Receptor complex", color = "Receptor Expression (Z-score)", size = "Ligand Expression (Z-score)")
Tcell_dotZ2

#save interactions to excel
pass_both_list <- pass_both %>% group_by(ligand_cluster, receptor_cluster) %T>% { group_keys(.) ->> group_names} %>% group_split()

# wb <- openxlsx::createWorkbook()
# 
# Map(function(data, nameofsheet){     
#   
#   addWorksheet(wb, nameofsheet)
#   writeData(wb, nameofsheet, data)
#   
# }, pass_both_list, apply(group_names, 1, paste, collapse = "_"))
# 
# saveWorkbook(wb, file = 'bespoke_interactions.xlsx')