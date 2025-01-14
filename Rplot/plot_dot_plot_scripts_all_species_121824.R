library(Matrix)
library(ggplot2)
args <- commandArgs(trailingOnly=T)


input_1_users_provide_configure_fl <- as.character(args[1])

input_2_users_provided_target_gene_ID_list_fl <- as.character(args[2])

input_3_zscore_all_organs_dir <- '../data/input_3_zscore_all_organs_dir'

input_4_accprop_all_organs_dir <-  '../data/input_4_accprop_all_organs_dir'

input_output_file <- as.character(args[3])


################
##read the files
source(input_1_users_provide_configure_fl)
#ipt_target_species <- 'Osativa'
#ipt_target_organ <- 'semroot'
#plot_width <- 10
#plot_height <- 10

# prefix <- paste0(ipt_target_species,'_',ipt_target_organ)

input_2_users_provided_target_gene_ID_list_dt <- read.delim(input_2_users_provided_target_gene_ID_list_fl,header = F)
colnames(input_2_users_provided_target_gene_ID_list_dt) <- c('geneID')

ipt_prepare_plot_rds <- readRDS(paste0(input_3_zscore_all_organs_dir,'/',ipt_target_species,'/opt_prepare_',ipt_target_organ,'_plot.rds'))

ipt_cellprop_mtx_dt <- read.delim(paste0(input_4_accprop_all_organs_dir,'/',ipt_target_species,'/opt_',ipt_target_organ,'_gene_acc_cellprop_mtx.txt'),header= T,row.names = 1)
# head(ipt_cellprop_mtx_dt)


##process for the Accprop
ipt_target_gene_dt <- input_2_users_provided_target_gene_ID_list_dt
shared_IDs <- intersect(rownames(ipt_cellprop_mtx_dt),ipt_target_gene_dt$geneID)
ipt_cellprop_mtx_target_dt <- ipt_cellprop_mtx_dt[shared_IDs,]
ipt_cellprop_mtx_target_dt <- as(as.matrix(ipt_cellprop_mtx_target_dt), "sparseMatrix")
ia <- as.data.frame(summary(ipt_cellprop_mtx_target_dt))
ia$i <- rownames(ipt_cellprop_mtx_target_dt)[as.numeric(ia$i)]
ia$j <- colnames(ipt_cellprop_mtx_target_dt)[as.numeric(ia$j)]
colnames(ia) <- c('geneID','celltype','Accprop')
opt_propAcc_dt <- ia
opt_propAcc_dt$geneID_celltype <- paste0(opt_propAcc_dt$geneID,':',opt_propAcc_dt$celltype)

##process for the zscore
shared_IDs <- intersect(rownames(ipt_prepare_plot_rds$zscore),ipt_target_gene_dt$geneID)
ipt_prepare_plot_dt <- ipt_prepare_plot_rds$zscore[shared_IDs,]
ipt_prepare_plot_dt <- as(as.matrix(ipt_prepare_plot_dt), "sparseMatrix")
ia <- as.data.frame(summary(ipt_prepare_plot_dt))
ia$i <- rownames(ipt_prepare_plot_dt)[as.numeric(ia$i)]
ia$j <- colnames(ipt_prepare_plot_dt)[as.numeric(ia$j)]
colnames(ia) <- c('geneID','celltype','zscore')
opt_zscore_dt <- ia
opt_zscore_dt$geneID_celltype <- paste0(opt_zscore_dt$geneID,':',opt_zscore_dt$celltype)

##merge accprop and zscore
merged_dt <- merge(opt_zscore_dt,opt_propAcc_dt,by.x = 'geneID_celltype',by.y = 'geneID_celltype')
merged_dt <- merged_dt[complete.cases(merged_dt), ]

merged_target_dt <- merged_dt[c('geneID.x','celltype.x','zscore','Accprop')]
colnames(merged_target_dt) <- c('geneID','celltype','zscore','Accprop')

p <- ggplot(merged_target_dt, aes(y=celltype, x = geneID, 
                                  color = zscore, size = Accprop)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  #scale_colour_gradient2(high = "gray0",  
  #                       mid = 'gray90',
  #                       low = "gray100",
  #                       midpoint = 0.5) +
  scale_colour_gradient2(midpoint = 0.5) +  ##original is 0.5
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text.x = element_text(colour = "black", size=12,angle = 45, vjust = 1,hjust =1),
        axis.text.y = element_text(colour = "black", size=12),
        axis.line.x = element_line(color="black", size = 1), 
        axis.line.y = element_line(color="black", size = 1)) +
  ggtitle("Z Score of Markers - Across Annotations")
# tiff(input_output_file, width = plot_width,height = plot_height)
# p
# dev.off()

ggsave(input_output_file, p, width = plot_width, height = plot_height)











