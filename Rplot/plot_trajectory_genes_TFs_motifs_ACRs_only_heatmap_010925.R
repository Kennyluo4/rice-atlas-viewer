##this pipeline will plot the target trajectories

##step01 libraries
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(dplyr)

##step02 libraries
library(scales)
library(phytools)

# load arguments
arg <- commandArgs(T)
print(arg)

input_0_users_provide_species_name <- as.character(arg[1])
##Os Zm

input_1_users_select_cate <- as.character(arg[2])
##Gene TF Mt ACR

input_2_users_select_path_name <- as.character(arg[3])
##Eseedling_SAMPhloem

input_3_users_provide_ID_fl <- as.character(arg[4])

# input_sorted_pt_dir <- as.character(arg[5])
input_sorted_pt_dir <-  '../data/input_sorted_pt_dir'
##it contains all the sorted pt files used for the plotting

# input_annot_TF_dir <- as.character(arg[6])
input_annot_TF_dir <- '../data/input_annot_TF_dir'
##it contains all species TF annotation to show their super family information

# input_gene_TF_motif_impute_sparse_rds_dir <- as.character(arg[7])
input_gene_TF_motif_impute_sparse_rds_dir <- '../data/input_gene_TF_motif_impute_sparse_rds_dir'
##it contains the impute sparse rds fl for all genes as well as motifs

# input_traject_meta_fl_dir <- as.character(arg[8])
input_traject_meta_fl_dir <- '../data/input_traject_meta_fl_dir'

input_output_dir <- as.character(arg[5])


###########################################################
##step01 plot the pt heatmap of show the order of target ID
ipt_target_ID_dt <- read.delim(input_3_users_provide_ID_fl,header = F)
head(ipt_target_ID_dt)
target_ID_list <- unique(ipt_target_ID_dt$V1)

ipt_pt_dt <- read.delim(paste0(input_sorted_pt_dir,'/',input_2_users_select_path_name,'/',input_2_users_select_path_name,'.',input_1_users_select_cate,'_pt_sorted.txt'))
ipt_pt_dt <- as.matrix(ipt_pt_dt)

##add the annotation to the gene name
if (input_1_users_select_cate == 'TF'){
  
  ipt_TF_annot_dt <- read.delim(paste0(input_annot_TF_dir,'/Os_TF_list.txt'))
  ipt_TF_annot_dt <- ipt_TF_annot_dt[c('Gene_ID','Family')]
  ipt_TF_annot_dt <- unique(ipt_TF_annot_dt)
  dim(ipt_TF_annot_dt)

  ipt_TF_annot_dt <- ipt_TF_annot_dt %>% distinct(Gene_ID, .keep_all = TRUE)
  rownames(ipt_TF_annot_dt) <- ipt_TF_annot_dt$Gene_ID
  ipt_TF_annot_dt <- ipt_TF_annot_dt[c('Family')]
  
  ##keep the same order as the merge will change the order
  ori_gene_order <- rownames(ipt_pt_dt)
  ipt_pt_merged_dt <- merge(ipt_pt_dt,ipt_TF_annot_dt,by = 'row.names', all.x=TRUE)
  ipt_pt_merged_dt$Row.names <- factor(ipt_pt_merged_dt$Row.names, levels = ori_gene_order)  
  ipt_pt_merged_dt <- ipt_pt_merged_dt[order(ipt_pt_merged_dt$Row.names), ]
  
  dim(ipt_pt_merged_dt)
  dim(ipt_pt_dt)
  ipt_pt_merged_dt$name_TFfam <- paste0(ipt_pt_merged_dt$Row.names,';',ipt_pt_merged_dt$Family)
  head(ipt_pt_merged_dt)
  rownames(ipt_pt_merged_dt) <- ipt_pt_merged_dt$name_TFfam
  ipt_pt_merged_dt <- ipt_pt_merged_dt[,1:(ncol(ipt_pt_merged_dt)-2)]
  ipt_pt_merged_dt <- ipt_pt_merged_dt[,-1]
  
  ##for the target ID list
  target_ID_list_dt <- as.data.frame(target_ID_list)
  rownames(target_ID_list_dt) <- target_ID_list_dt$target_ID_list
  mergedt_target_ID_dt <- merge(target_ID_list_dt,ipt_TF_annot_dt,by ='row.names',all.x=TRUE)
  mergedt_target_ID_dt$name_TFfam <- paste0(mergedt_target_ID_dt$Row.names,';',mergedt_target_ID_dt$Family)
  
  target_ID_list <- mergedt_target_ID_dt$name_TFfam
  ipt_pt_dt <- ipt_pt_merged_dt
}


##we will check which TF we will plot
if (input_1_users_select_cate == 'Mt'){
  ##for the Mt we will modify the ID
  all_ID_motifnm_list <- rownames(ipt_pt_dt)
  all_ID_motifnm_dt <- as.data.frame(all_ID_motifnm_list)
  all_ID_motifnm_dt$ID <- gsub('_.+','',all_ID_motifnm_dt$all_ID_motifnm_list)
  all_ID_motifnm_dt$IDconcise <- gsub('\\..+','',all_ID_motifnm_dt$ID)
  all_ID_motifnm_dt$motif <- gsub('.+_','',all_ID_motifnm_dt$all_ID_motifnm_list)
  head(all_ID_motifnm_dt)
  
  all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$motif %in% target_ID_list,]
  if (nrow(all_ID_motifnm_dt_flt) == 0){
    
    all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$IDconcise %in% target_ID_list,]
    
    if (nrow(all_ID_motifnm_dt_flt) == 0){
      all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$ID %in% target_ID_list,]
      
      if (nrow(all_ID_motifnm_dt_flt) == 0){
        
        all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$all_ID_motifnm_list %in% target_ID_list,]
        
      }
      
    }
    
  }
  
  target_rownames <- all_ID_motifnm_dt_flt$all_ID_motifnm_list

}else{
  
  if (input_1_users_select_cate == 'ACR'){
    ##for the ACR we will modify the Chr information
    all_ID_motifnm_list <- rownames(ipt_pt_dt)
    all_ID_motifnm_list <- gsub('chrChr','Chr',all_ID_motifnm_list)
    target_rownames <- intersect(all_ID_motifnm_list,target_ID_list)
  }else{
    ##for the TF and gene we will keep the same
    all_ID_motifnm_list <- rownames(ipt_pt_dt)
    target_rownames <- intersect(all_ID_motifnm_list,target_ID_list)
  }
}

##assign to a new z
z <- ipt_pt_dt
head(z)

##change chrChr to Chr
rownames(z) <- gsub('chrChr','Chr',rownames(z))
head(z)

##add the rank information in the z
rownames(z) <- paste0(rownames(z) ,' (rank ',seq_len(nrow(z)),')')

##add the rank information to the target_rownames
z_dt <- as.data.frame(rownames(z))
colnames(z_dt) <- c('name')
z_dt$ori_name <- gsub(' \\(rank.+\\)','',z_dt$name)
z_dt_target <- z_dt[z_dt$ori_name %in% target_rownames,]
final_target_rownames <- z_dt_target$name

ha = rowAnnotation(foo = anno_mark(at = which(rownames(z) %in% final_target_rownames),
                                   labels = rownames(z)[rownames(z)%in%final_target_rownames]))

if (input_1_users_select_cate == 'Mt'){
  cate_col = colorRampPalette(c(rev(brewer.pal(8,'PiYG'))))(100)
}
if (input_1_users_select_cate == 'TF'){
  cate_col = colorRampPalette(c('grey80','grey75',brewer.pal(8,'BuPu')[2:8]))(100)
}
if (input_1_users_select_cate == 'Gene'){
  cate_col = viridis(100)
}
if (input_1_users_select_cate == 'ACR'){
  cate_col = colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
}

# heatmap
par(mar = c(7, 7, 7, 7))
# pdf(paste0(input_output_dir,'/',input_2_users_select_path_name,".trajectory",input_1_users_select_cate,".pdf"), width=8, height=8)
tiff(paste0(input_output_dir,'/',input_2_users_select_path_name,".trajectory.",input_1_users_select_cate,"_heatmap.tiff"), width =24, height = 30,res = 200, units = 'cm')

par(cex = 1.5)

##remove the zero for the column sum
#column_sums <- colSums(z)
# Identify columns with sum equal to 0
#columns_to_drop <- names(column_sums[column_sums == 0])
# Drop those columns
#z <- z[, !(names(z) %in% columns_to_drop)]
#dim(z)
p <- Heatmap(z, name = paste0(input_1_users_select_cate,' ',"enrichment"), cluster_rows = F,cluster_columns=F, 
        #cluster_columns=col.dendro, 
        #top_annotation = ha.col,
        col = cate_col,
        #col=colorRamp2(as.numeric(n.range),c("white","#fabd73","#f7aa4d","darkorange","firebrick3")),
        #col=colorRamp2(as.numeric(n.range),c("#440154","#3b528b","#21918c","#5ec962","#fde725")),
        use_raster=T,
        show_row_names = F,
        show_column_names = F,
        right_annotation = ha,
        row_title = paste0(input_1_users_select_cate,' count = ',nrow(z)),
        column_title = 'Pseudotime (begin to end)')
print(p)
dev.off()


















