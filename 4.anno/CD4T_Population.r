library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(SingleR)
library(scater)
library(ExperimentHub)
library(AnnotationHub)
library(AnnotationDbi)
library(celldex)
library(reshape2)
library(ggplot2)
library(ggsci)
library(tidyr)
library(openxlsx)
source('../Rscript/stat.r')
dtVar <- Sys.Date() 
dtVar <- as.Date(dtVar, tz="UTC")

workpath <- "../3.combination.DeDoublets/"
dir.create(str_c('Results_',dtVar))
file_name = 'combined'
# file_name = 'DeDoublets_C2_D21'
rdsfile <- str_c(workpath,"/",file_name,"_t-SNE_30PCA_0.6Resolution/",file_name,"_t-SNE_30PCA_0.6Resolution.AnnoManual.rds")
print(file.exists(rdsfile))

immune.combined <- readRDS(rdsfile)
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- 'celltype'
sample_list <- c('G1','G2','G3','G4')
immune.combined$sample_label <- factor(immune.combined$sample_label, levels=sample_list)
p1 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, label.size = 4) + NoLegend()
ggsave(str_c('Results_',dtVar, '/Tsubset_DimPlot_clusters.',dtVar,'.check.pdf'), p1, width=5.5, height=5.5)


### T细胞亚群细分
Idents(immune.combined) <- immune.combined$celltype
pbmc_T <- subset(immune.combined, idents=c('CD4 T cell'))#))
DefaultAssay(pbmc_T) <- "integrated"
# Run the standard workflow for visualization and clustering
pbmc_T <- ScaleData(pbmc_T, verbose = FALSE)
pbmc_T <- RunPCA(pbmc_T, npcs = 30, verbose = FALSE)
pbmc_T <- RunTSNE(pbmc_T, reduction = "pca", dims = 1:30)
pbmc_T <- FindNeighbors(pbmc_T, reduction = "pca", dims = 1:30)
pbmc_T <- FindClusters(pbmc_T, resolution = 0.6)
p2 <- DimPlot(pbmc_T, reduction = "tsne", label = TRUE, repel = TRUE) + theme_bw()
ggsave(str_c('Results_',dtVar, '/Tsubset_DimPlot_clusters.',dtVar,'.pdf'), p2, width=5.5, height=4.5)

Idents(pbmc_T) = 'sample_label'
Idents(pbmc_T) = factor(Idents(pbmc_T), levels = c('G1','G2','G3','G4'))
p1 <- DimPlot(pbmc_T, reduction = "tsne", label = FALSE, pt.size = 0.01) + 
      ggplot2::scale_color_manual(values = c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_bw() + ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Tsubset_DimPlot_samples.',dtVar,'.pdf'), p1, width=5, height=3)


DefaultAssay(pbmc_T) <- "RNA"
Idents(pbmc_T) <- factor(pbmc_T$seurat_clusters, levels=c(3,9,5,4,2,0,1,8,6,7))
      markers <- c('Cd3e','Cd4','Cd8a','Cd44','Foxp3','Il2ra','Ctla4','Ikzf2',
                        'Lag3','Icos','Tnfrsf4','Tnfrsf9','Tnfrsf18','Klrg1','Il10','Gzmb',
                        'Sell','Tcf7','Bcl2','Il7r','Il6ra','Ccr6','Cd24a','Cd69','Xcl1','Nr4a1')
      markers = intersect(markers, rownames(pbmc_T))
p1 <- DotPlot(pbmc_T, features = markers, dot.scale = 3, col.min=-0.5, cols=c('white','red','red4')) + RotatedAxis() +theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave(str_c('Results_',dtVar, '/Tsubset.DotPlot_clusters.',dtVar,'.pdf'), p1, width=7, height=4)
## conserved markers
# for (i in names(table(pbmc_T$seurat_clusters))) {
#       cluster_i.markers <- FindConservedMarkers(pbmc_T, ident.1 = i, grouping.var = "sample_label", verbose = FALSE)
#       write.xlsx(cluster_i.markers, str_c('t_c',i,'.xlsx'),rowNames=TRUE)
# }

table(pbmc_T$seurat_clusters)
pbmc_T$celltype = dplyr::case_when(
  pbmc_T$seurat_clusters %in% c(3,9) ~ "CD4 Tem",
  pbmc_T$seurat_clusters %in% c(5) ~ "Effector Treg",
  pbmc_T$seurat_clusters %in% c(4) ~ "Central Treg",
  pbmc_T$seurat_clusters %in% c(0,1,2,8) ~ "CD4 Naive",
  pbmc_T$seurat_clusters %in% c(6) ~ "Activated CD4 T",
  pbmc_T$seurat_clusters %in% c(7) ~ "Unknown")
table(pbmc_T$celltype)

marjor_cluster <- c("CD4 Naive","CD4 Tem","Effector Treg","Central Treg","Activated CD4 T", "Unknown")
pbmc_T$celltype <- factor(pbmc_T$celltype, levels = marjor_cluster)

Idents(pbmc_T) = 'celltype'
p1 <- DimPlot(pbmc_T, reduction = "tsne", label = FALSE, label.size = 4, repel=TRUE, cols=c(pal_nejm(alpha = 0.4)(8), '#7E6148E5', '#B09C85E5','grey'))+theme_bw() +
      theme(legend.position="right") + ggplot2::theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(str_c('Results_',dtVar, '/Tsubset_Anno1_SeuratClusters.',dtVar,'.pdf'), p1, width=4.8, height=3)
p1 <- DimPlot(pbmc_T, reduction = "tsne", label = FALSE, label.size = 4, repel=FALSE, cols=c(pal_nejm(alpha = 0.3)(8), '#7E6148E5', '#B09C85E5','grey'))+theme_bw() +
      theme(legend.position="bottom", axis.text=element_blank())+ NoLegend()+ylab('')+xlab('')
ggsave(str_c('Results_',dtVar, '/Tsubset_Anno1_SeuratClusters.',dtVar,'.no.pdf'), p1, width=3, height=3)
p1 <- DotPlot(pbmc_T, features = markers, dot.scale = 3, col.min=-0.5, cols=c('white','red','red4')) + RotatedAxis() +
      theme_bw()+theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5))+ ylab('Cell type')
ggsave(str_c('Results_',dtVar, '/Tsubset.DotPlot_celltype.',dtVar,'.pdf'), p1, width=7, height=3.8)

### 保存文件
saveRDS(pbmc_T, file = str_c(workpath,"/",file_name,"_t-SNE_30PCA_0.6Resolution/",file_name,"_t-SNE_30PCA_0.6Resolution.AnnoManual.Tsubset.rds"))




markers.to.plot = c("Foxp3",'Ctla4','Ikzf2','Icos','Klrg1','Lag3','Ccl5','Gimap7','Pglyrp1','Gimap3','Pdcd1','Cd200r1')#,'Il10','Fgl2','Mmp9','Il10','Cxcr3','Capg','Gzmb','Ikzf2')
pbmc_effetcorTreg = subset(pbmc_T, idents=c("Effector Treg"))
pbmc_effetcorTreg$sample_label = factor(pbmc_effetcorTreg$sample_label, levels=c('G1','G2','G3','G4'))
Idents(pbmc_effetcorTreg) = "sample_label"
p2 = VlnPlot(pbmc_effetcorTreg, features = markers.to.plot, ncol=6, cols = c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) 
ggsave(str_c('Results_',dtVar, '/Tsubset.VlnPlot_celltype.',dtVar,'.pdf'), p2, width=8, height=6)


# effetor_treg_data = FetchData(pbmc_effetcorTreg, c("Foxp3",'Ctla4','Ikzf2','Icos','Klrg1','Lag3','Ccl5','Gimap7','Pglyrp1','Gimap3','Pdcd1','Cd200r1', 'sample_label'))
# pvalue_c = c()
# for(gene in  c("Foxp3",'Ctla4','Ikzf2','Icos','Klrg1','Lag3','Ccl5','Gimap7','Pglyrp1','Gimap3','Pdcd1','Cd200r1')){
#   pvalue_c = c(pvalue_c, wilcox.test(effetor_treg_data[effetor_treg_data$sample_label == 'D21_CTX', gene], 
#               effetor_treg_data[effetor_treg_data$sample_label == 'D21_UR_CTX', gene])$p.value)
# }
# cbind( c("Foxp3",'Ctla4','Ikzf2','Icos','Klrg1','Lag3','Ccl5','Gimap7','Pglyrp1','Gimap3','Pdcd1','Cd200r1'),pvalue_c)




##### Roe
df = as.data.frame(table(pbmc_T@meta.data[, c("celltype","sample_label")]))
df <- aggregate(df$Freq, list(df$celltype, df$sample_label), sum)
colnames(df) <- c('celltype', 'sample_label', 'Freq')

df_wide_T <- df %>%
      pivot_wider(
      names_from = sample_label,
      values_from = Freq
      ) %>%
      arrange(celltype) %>%
      mutate(across(-celltype, 
                  .fns = list(ratio = ~ .x / sum(.x)),
                  .names = "{.col}_ratio"))

      
Tsubset = as.data.frame(df_wide_T); rownames(Tsubset) = Tsubset$celltype
Tsubset = Tsubset[, 2:5]
total_T_cells = apply(Tsubset,2,sum)
total_T_clusters = apply(Tsubset, 1, sum); total_t_cells = sum(total_T_clusters)
sample_groups_new = colnames(Tsubset)
total_cells = apply(total_cell,2,sum)

chisq_T_result3 = calculate_roe(total_T_clusters, Tsubset, total_T_cells, total_t_cells, sample_groups_new)
labels = c("c_ur", "c_urctx", "c_ctx",  "ctx_ur", "ctx_urctx", "ur_urctx")
chisq_T_result3_sym = add_symnum_label(chisq_T_result3, labels); colnames(chisq_T_result3_sym) = gsub('CellNum','real_vs_expected',colnames(chisq_T_result3_sym))
chisq_T_result3_melt = melt(chisq_T_result3[,c('cl', colnames(chisq_T_result3)[2:5])])

chisq_T_result3_melt$variable = as.character(chisq_T_result3_melt$variable)
chisq_T_result3_melt$variable = factor(chisq_T_result3_melt$variable, levels=c('G1','G2','G3','G4'))
p1 <- ggplot(data=chisq_T_result3_melt, aes(x=variable, y=value, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+geom_hline(yintercept=1, linetype='dashed', col = 'red')+theme(axis.text.x=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7))+xlab('')+
      facet_grid(.~cl)+ylab('Roe (in CD4T)')
ggsave(str_c('Results_',dtVar, '/Tsubset.Roe.',dtVar,'.pdf'), p1, width=7, height=3)

total_cells = apply(total_cell,2,sum)
for(s in names(total_cells)){
    Tsubset[, str_c(gsub('CellNum','',s), ' in Total')] = Tsubset[,s]/total_cells[s]
    Tsubset[, str_c(gsub('CellNum','',s), ' in CD4T')] = Tsubset[,s]/total_T_cells[s]
}
Tsubset_in_CD4T = Tsubset[, grep('CD4T', colnames(Tsubset))]
Tsubset_in_CD4T$cl = rownames(Tsubset_in_CD4T)
prob_T_result3_melt = melt(Tsubset_in_CD4T, id.vars='cl')
prob_T_result3_melt$variable = gsub(" in CD4T","",prob_T_result3_melt$variable)
prob_T_result3_melt$variable = factor(prob_T_result3_melt$variable, levels=c('G1','G2','G3','G4'))
p1_1 <- ggplot(data=prob_T_result3_melt[prob_T_result3_melt$cl == 'CD4 Naive', ], aes(x=variable, y=value * 100, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+theme(axis.text.x=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7), legend.position='none')+xlab('')+
      facet_grid(.~cl)+ylab('Proportion in CD4T (%)')
p1_2 <- ggplot(data=prob_T_result3_melt[prob_T_result3_melt$cl != 'CD4 Naive', ], aes(x=variable, y=value * 100, color=variable)) +
      geom_point() + geom_bar(stat="identity", fill='white', linewidth=0.7, width=0.7)+ scale_color_manual(values= c("#999999","#FFBE7A","#8ECFC9","#FE7C6D")) + 
      theme_classic()+theme(axis.text.x=element_blank(), axis.title.y=element_blank(), strip.background=element_blank(), strip.text.x=element_text(size=7))+xlab('')+
      facet_grid(.~cl)
g = plot_grid(p1_1, p1_2, nrow=1, rel_widths = c(0.4,1.8))
ggsave(str_c('Results_',dtVar, '/Tsubset.Proportion.',dtVar,'.pdf'), g, width=9, height=3)

write.xlsx(list(Tsubset, chisq_T_result3_sym, chisq_T_result3), str_c('Results_',dtVar, '/Tsubset.Roe_inCD4T.',dtVar,'.xlsx'), rowNames=TRUE)
