library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(SingleR)
library(scater)
library(ExperimentHub)
library(AnnotationHub)
library(AnnotationDbi)
library(celldex)
library(reshape2)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(tidyr)
source('../Rscript/stat.r')

dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

dir.create(str_c('Results_',dtVar))
workpath <- "../3.combination.DeDoublets/"
file_name = 'combined'
# file_name = 'DeDoublets_C2_D21'
rdsfile <- str_c(workpath,"/",file_name,"_t-SNE_30PCA_0.6Resolution/",file_name,"_t-SNE_30PCA_0.6Resolution.rds")
print(file.exists(rdsfile))

immune.combined <- readRDS(rdsfile)
DefaultAssay(immune.combined) <- "RNA"
sample_list <- c('G1','G2','G3','G4')
immune.combined$sample_label <- factor(immune.combined$sample_label, levels=sample_list)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4) + NoLegend()
Idents(immune.combined) = 'sample_label'
Idents(immune.combined) = factor(Idents(immune.combined), levels = c('G1','G2','G3','G4'))
p1 <- DimPlot(immune.combined, reduction = "tsne", label = FALSE, pt.size = 0.01) + 
      ggplot2::scale_color_manual(values = c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_bw() + ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Total.DimPlot_samples.',dtVar,'.pdf'), p1, width=5, height=3)


Idents(immune.combined) <- "seurat_clusters"
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, pt.size = 0.01) + 
      theme_bw() + ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Total.DimPlot_clusters.',dtVar,'.pdf'), p1, width=4.5, height=3)


Idents(immune.combined) <- factor(Idents(immune.combined), levels=c(12,1,3,17,5,0,11,6,19,2,4,8,16,9,21,13,7,10,20,18,14,15))
marker <- c('Cd19','Cd79a','Cd3e','Cd4','Cd8a','Nkg7','Klrb1c','Ccr2','Ly6g','Fcgr3','Cxcr2','Arg2','Il1b','Lyz2','Cd14','Csf1r','Adgre1','Spic','Mrc1','Sirpa','Cd86',
            'Itgax','Siglech','Bst2','Clec9a','Xcr1','Cst3','Fcer1g','Sdc1','Tnfrsf17','Kit','Cd34','Mpo','Gata1','Klf1','Epor','Tfrc','Hba-a1','Hbb-bt')
p1 <- DotPlot(immune.combined, features = unique(c(marker)), dot.scale = 3, col.min=-0.5, cols=c('white','red','red4')) + RotatedAxis() +
      theme_bw() + theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('Results_',dtVar, '/Total.DotPlot_clusters.',dtVar,'.pdf'), p1, width=10, height=4)

# table(immune.combined$seurat_clusters)
table(immune.combined$seurat_clusters)
immune.combined$celltype = dplyr::case_when(
  immune.combined$seurat_clusters %in% c(12,1,3,17,5,0,11) ~ "B cell",
  immune.combined$seurat_clusters %in% c(6,19,2) ~ "CD4 T cell",
  immune.combined$seurat_clusters %in% c(4) ~ "CD8 T cell",
  immune.combined$seurat_clusters %in% c(8,16) ~ "NK T cell",
  immune.combined$seurat_clusters %in% c(9) ~ "NK cell",
  immune.combined$seurat_clusters %in% c(13) ~ "MonoNeutro",
  immune.combined$seurat_clusters %in% c(7) ~ "Macrophage",
  immune.combined$seurat_clusters %in% c(10,20) ~ "DC",
  immune.combined$seurat_clusters %in% c(18) ~ "Plasma cell",
  immune.combined$seurat_clusters %in% c(14) ~ "Progenitor cell",
  immune.combined$seurat_clusters %in% c(15) ~ "Erythrocyte")
table(immune.combined$celltype)

marjor_cluster <- c('B cell','CD4 T cell','CD8 T cell','NK T cell','NK cell', 'MonoNeutro','Macrophage','DC',"Plasma cell", 'Progenitor cell', 'Erythrocyte')
immune.combined$celltype <- factor(immune.combined$celltype, levels = marjor_cluster)

Idents(immune.combined) = 'celltype'
celltype_colors <- c('#7E6148E5', '#B09C85E5', pal_nejm(alpha = 0.2)(8)[2:8], 'grey', pal_npg(alpha=0.2)(1))
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4, repel=TRUE, cols=celltype_colors)+
            theme_bw() + theme(legend.position="bottom")+ 
            ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Total.Anno1_SeuratClusters.',dtVar,'.pdf'), p1, width=4.5, height=5.5)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = FALSE, label.size = 4, repel=FALSE, cols=celltype_colors)+
            theme_bw() + theme(legend.position="none", axis.text=element_blank())+ NoLegend()+ylab('')+xlab('')+ 
            ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Total.Anno1_SeuratClusters.',dtVar,'.no.pdf'), p1, width=3, height=3)

p1 <- DotPlot(immune.combined, features = unique(c(marker)), dot.scale = 3, col.min=-0.5, cols=c('white','red','red4')) + RotatedAxis() +
      theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) + ylab('Cell type')
ggsave(str_c('Results_',dtVar, '/Total.DotPlot_celltype.',dtVar,'.pdf'), p1, width=9, height=4)


### Calculating Cell numbers and ratio
Idents(immune.combined) = 'celltype'
rds = immune.combined
rds2 = immune.combined
colors = celltype_colors
sample_labels = c('G1','G2','G3','G4')
sample_groups = sample_labels
celltypes = marjor_cluster
group_colors = c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")
name = str_c('Results_',dtVar, '/Total.', dtVar)
cellnumber_summary(rds, rds2, sample_labels, sample_labels, marjor_cluster, colors, group_colors, name)

### cellcycle
library(Hmisc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- c('Cenpu', capitalize(tolower(cc.genes$s.genes)))
g2m.genes <- c('Pimreg','Jpt1', capitalize(tolower(cc.genes$g2m.genes)))
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# s.genes_old <- c('Cenpu', capitalize(tolower(cc.genes$s.genes)))
# g2m.genes_old <- c('Pimreg','Jpt1', capitalize(tolower(cc.genes$g2m.genes)))

# biomart_mouse2human = read.table('/Users/xingyun/Downloads/mart_export.mouse2human.txt', sep='\t', head=TRUE)
# biomart_mouse2human = unique(biomart_mouse2human[, c('Gene.name','Human.gene.name')])

# biomart_s_genes = biomart_mouse2human[which(biomart_mouse2human$Human.gene.name %in% cc.genes$s.genes), 'Gene.name']
# biomart_g2m_genes = biomart_mouse2human[which(biomart_mouse2human$Human.gene.name %in% cc.genes$g2m.genes), 'Gene.name']

# intersect(s.genes,biomart_mouse2human$Human.gene.name )
# setdiff(s.genes,biomart_mouse2human$Human.gene.name )
# ## ncbi-genes: MLF1IP > Cenpu (mouse); POLD3 > Pold3 (mouse); ATAD2 > Atad2 (mouse)
# intersect(s.genes_old, biomart_s_genes)
# setdiff(s.genes_old, biomart_s_genes)
# ## "Cenpu"  "Mlf1ip" "Pold3"  "Atad2"
# intersect(g2m.genes,biomart_mouse2human$Human.gene.name )
# setdiff(g2m.genes,biomart_mouse2human$Human.gene.name )
# ## ncbi-genes: FAM64A > Pimreg (mouse); HN1 > Jpt1 (mouse)
# intersect(g2m.genes_old, biomart_g2m_genes )
# setdiff(g2m.genes_old, biomart_g2m_genes )
# ## "Pimreg" "Jpt1"   "Fam64a" "Hn1"
Idents(immune.combined) = 'celltype'
saveRDS(immune.combined, file = str_c(workpath,"/",file_name,"_t-SNE_30PCA_0.6Resolution/",file_name,"_t-SNE_30PCA_0.6Resolution.AnnoManual.rds"))







### Roe and Proportion
source('../Rscript/stat.r')
df = as.data.frame(table(immune.combined@meta.data[, c("celltype","sample_label")]))
df_wide_total <- df %>%
      pivot_wider(
      names_from = sample_label,
      values_from = Freq
      ) %>%
      arrange(celltype)
total_cell = as.data.frame(df_wide_total); rownames(total_cell) = total_cell$celltype
total_cell = total_cell[, 2:5]
total_cells = apply(total_cell, 2, sum); total = sum(total_cells)
total_clusters = apply(total_cell, 1, sum)
sample_groups_new = colnames(total_cell)
chisq_total_result2 = calculate_roe(total_clusters, total_cell, total_cells, total, sample_groups_new)
labels = c("c_ur", "c_urctx", "c_ctx",  "ctx_ur", "ctx_urctx", "ur_urctx")

chisq_total_result2_sym = add_symnum_label(chisq_total_result2, labels)
colnames(chisq_total_result2_sym) = gsub('CellNum','real_vs_expected',colnames(chisq_total_result2_sym))
total_cells = apply(total_cell,2,sum)
chisq_total_result2_melt = melt(chisq_total_result2[,c('cl', colnames(chisq_total_result2)[2:5])])
chisq_total_result2_melt$variable = factor(chisq_total_result2_melt$variable, levels=c('G1','G2','G3','G4'))
p1 <- ggplot(data=chisq_total_result2_melt, aes(x=variable, y=value, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+geom_hline(yintercept=1, linetype='dashed', col = 'red')+
      theme(axis.text.x=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7))+
      xlab('') + ylab('Roe (in total)') + facet_grid(.~cl)
ggsave(str_c('Results_',dtVar, '/Total.Roe.',dtVar,'.pdf'), p1, width=9, height=3)



for(s in names(total_cells)){
    total_cell[, str_c(gsub('CellNum','',s), ' in Total')] = total_cell[,s]/total_cells[s]
}
total_cell_in_Total = total_cell[, grep('Total', colnames(total_cell))]
total_cell_in_Total$cl = rownames(total_cell_in_Total)
prob_Total_result3_melt = melt(total_cell_in_Total, id.vars='cl')
prob_Total_result3_melt$variable = gsub(" in Total","",prob_Total_result3_melt$variable)
prob_Total_result3_melt$variable = factor(prob_Total_result3_melt$variable, levels=c('G1','G2','G3','G4'))
p1_1 <- ggplot(data=prob_Total_result3_melt[prob_Total_result3_melt$cl == 'B cell', ], aes(x=variable, y=value * 100, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+theme(axis.text.x=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7), legend.position='none')+
      xlab('') + ylab('Proportion in Total (%)') + facet_grid(.~cl)
p1_2 <- ggplot(data=prob_Total_result3_melt[prob_Total_result3_melt$cl != 'B cell', ], aes(x=variable, y=value * 100, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+theme(axis.text.x=element_blank(), axis.title.y=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7))+
      xlab('') + facet_grid(.~cl)
g = plot_grid(p1_1, p1_2, nrow=1, rel_widths = c(0.25,1.8))
ggsave(str_c('Results_',dtVar, '/Total.Proportion.',dtVar,'.pdf'), g, width=9, height=3)

write.xlsx(list(total_cell, chisq_total_result2_sym), str_c('Results_',dtVar, '/Total.Roe_intotal.',dtVar,'.xlsx'), rowNames=TRUE)




### Myeloid marker genes

# 1. Macrophage
immunosuppressive_genes = c('Marco', 'Il10', 'Socs3', 'Cd163', 'Trem2', 'Lgmn','C1qa','C1qb','Il18','Lap3','S100a8') # 'Mertk', 
proinflammatory_markers = c('Itgax','Icosl','Irf8','Ikzf1','Il16','Tlr4')

pbmc_macro <- subset(immune.combined, idents=c('Macrophage'))#))
DefaultAssay(pbmc_macro) <- "RNA"
pbmc_macro$sample_label = factor(pbmc_macro$sample_label, levels=c('G1','G2','G3','G4'))
Idents(pbmc_macro) = 'sample_label'

p3 <- DotPlot(subset(pbmc_macro, idents=c('G2','G3','G4')), features = unique(c(immunosuppressive_genes, proinflammatory_markers)), 
      dot.scale = 3, , cols = c("blue", "red", "grey") ) + RotatedAxis() +theme_bw()+
      theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5), legend.position='bottom', axis.title.y=element_blank())
ggsave(str_c('Results_',dtVar, '/Mye_Macro.3samples_selectedmarkers',dtVar,'.pdf'), p3, width=7, height=2.5)


# 2. DC
immuno_infla_markers = c('Cd63','Bcl2a1d','Bcl7c','Bcl2a1a','Ccr3','Fcgr3','Irf8','Csf2rb','Vcam1')
pbmc_DC <- subset(immune.combined, idents=c('DC'))#))
DefaultAssay(pbmc_DC) <- "RNA"
pbmc_DC$sample_label = factor(pbmc_DC$sample_label, levels=c('G1','G2','G3','G4'))
Idents(pbmc_DC) = 'sample_label'

p3 <- DotPlot(subset(pbmc_DC, idents=c('G2','G3','G4')), features =immuno_infla_markers, 
      dot.scale = 3, , cols = c("blue", "red", "grey") ) + RotatedAxis() +theme_bw()+
      theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5), legend.position='bottom', axis.title.y=element_blank())
ggsave(str_c('Results_',dtVar, '/Mye_DC.3samples_selectedmarkers',dtVar,'.pdf'), p3, width=7, height=2.5)


# 3. Mono/Neutro
immuno_infla_markers = c('Ptgs2','Atf3','Arg2','Cebpb','Cebpz','Tnfrsf1b','Tnfaip3','Ifitm3','Clec4a1','Ifitm6','Irf8')
pbmc_mononeu <- subset(immune.combined, idents=c('MonoNeutro'))#))
DefaultAssay(pbmc_mononeu) <- "RNA"
pbmc_mononeu$sample_label = factor(pbmc_mononeu$sample_label, levels=c('G1','G2','G3','G4'))
Idents(pbmc_mononeu) = 'sample_label'

p3 <- DotPlot(subset(pbmc_mononeu, idents=c('G2','G3','G4')), features =immuno_infla_markers, 
      dot.scale = 3, , cols = c("blue", "red", "grey") ) + RotatedAxis() +theme_bw()+
      theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5), legend.position='bottom', axis.title.y=element_blank())
ggsave(str_c('Results_',dtVar, '/Mye_MonoNeu.3samples_selectedmarkers',dtVar,'.pdf'), p3, width=7, height=2.5)
