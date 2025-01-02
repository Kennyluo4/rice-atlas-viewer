library(Matrix)
library(pheatmap)
args <- commandArgs(trailingOnly=T)

input_1_users_provide_configure_fl <- as.character(args[1])

input_2_users_provide_motif_ID_list_fl <- as.character(args[2])

# input_3_required_raw_data_dir <- as.character(args[3])
input_3_required_raw_data_dir <- 'Rplot/data/input_3_required_raw_data_dir'

# input_4_motif_enrichment_all_organ_dir <- as.character(args[4])
input_4_motif_enrichment_all_organ_dir <- 'Rplot/data/input_4_motif_enrichment_all_organ_dir'

input_output_dir <- as.character(args[3])

################
##read the files
source(input_1_users_provide_configure_fl)
#ipt_target_species <- 'Osativa'
#ipt_target_organ <- 'semroot'
#plot_width <- 10
#plot_height <- 10
#FDR <- 0.01

prefix <- paste0(ipt_target_species,'_',ipt_target_organ)

ipt_target_motif_ID_dt <- read.delim(input_2_users_provide_motif_ID_list_fl,header = F)
colnames(ipt_target_motif_ID_dt) <- 'motifname'
ipt_target_motif_name_list <- ipt_target_motif_ID_dt$motifname

meta_name_dt <- read.delim(paste0(input_3_required_raw_data_dir,'/opt_motif_common_name.txt'),header = F)

TFmotif_dt <- read.delim(paste0(input_4_motif_enrichment_all_organ_dir,'/',ipt_target_species,'/','opt_',ipt_target_organ,'_motif_enrichment_cluster.txt'),row.names = 1)

TFmotif_dt <- merge(TFmotif_dt,meta_name_dt,by.x='motif',by.y= 'V1')
head(TFmotif_dt)


##processing
TFmotif_dt <- TFmotif_dt[TFmotif_dt$FDR < FDR,]
TFmotif_dt$log10FC <- log10(TFmotif_dt$fc)
TFmotif_dt <- TFmotif_dt[c('V2','organct','log10FC')]##we will use the beta other than the fc
colnames(TFmotif_dt) <- c('motif','organct','enrichment_score')

dim(TFmotif_dt)
write.table(TFmotif_dt,paste0(input_output_dir,'/opt_',prefix,'_enrichment_score.txt'),quote = F,sep = '\t')
head(TFmotif_dt)

TFmotif_dt <- read.table(paste0(input_output_dir,'/opt_',prefix,'_enrichment_score.txt'),stringsAsFactors = T)
dim(TFmotif_dt)
head(TFmotif_dt)
##transfer to matrix
TFmotif_mtx <- sparseMatrix(i=as.numeric(TFmotif_dt$motif),
                            j=as.numeric(TFmotif_dt$organct),
                            x=as.numeric(TFmotif_dt$enrichment_score),
                            dimnames=list(levels(TFmotif_dt$motif), levels(TFmotif_dt$organct)))
mat <- as(TFmotif_mtx, "dgCMatrix")
dim(mat)

shared_motifID <- intersect(rownames(mat),ipt_target_motif_name_list)
mat <- mat[shared_motifID,]


##plot the heatmap
##order the enrichment by order of values
signal_enrichment_matrix <- mat
z <- t(as.matrix(scale(t(as.matrix(signal_enrichment_matrix)))))


o.order <- colnames(z)
col.clust <- hclust(as.dist(1-cor(z)))
col.o <- col.clust$order
col.dendro <- as.dendrogram(col.clust)

signal_enrichment_matrix_reorder <- signal_enrichment_matrix[,col.o]
signal_enrichment_matrix_reorder <- signal_enrichment_matrix_reorder[order(apply(signal_enrichment_matrix_reorder, 1, which.max), decreasing=F),]
head(signal_enrichment_matrix_reorder)
signal_enrichment_matrix_reorder <- t(signal_enrichment_matrix_reorder)
message(plot_width)
message(plot_height)
# pdf(paste0(input_output_dir,'/',prefix,'_FDR',FDR,'_heat_map_motif_enrich.pdf'),width = plot_width, height = plot_height)
tiff(paste0(input_output_dir,'/',prefix,'_FDR',FDR,'_heat_map_motif_enrich.tiff'), width = plot_width, height = plot_height,res = 250, units = 'cm')
signal_enrichment_matrix_reorder <- as.matrix(signal_enrichment_matrix_reorder)
p <- pheatmap(signal_enrichment_matrix_reorder,
              #scale="row",
              #gaps_row = c(6,12,18,24,30,36,42,48,54,60),
              #gaps_col = c(3,5,8,11,12,13,14,15,16,17),
              cluster_cols = F,
              cluster_rows = F,
              border_color = "black",
              #border_color='black',
              fontsize_col = 8,
              fontsize_row = 8,
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              #annotation_col = annotdf_col,
              #color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow"))(50),
              na_col = "white",
              #color=colorRampPalette(c("white","#df5757","#df5757","#df5757"))(100),
              #color=colorRampPalette(c("white","#ed9f9f","#ed8787","#fa2525"))(100),
              #color=colorRampPalette(c("#440154","#3b528b","#21918c","#5ec962","#fde725"))(100),
              #color=colorRampPalette(c("#fcfdbf","#fc8961","#b73779","#51127c","#000004"))(100),
              color=colorRampPalette(c("white","#fc8961","#b73779","#51127c","#000004"))(100),
              
              #color=colorRampPalette(c("#df5757", "lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",
              #                           "lightgoldenrodyellow",'lightgoldenrodyellow','lightgoldenrodyellow',"lightgoldenrodyellow","lightgoldenrodyellow","lightgoldenrodyellow",'white','white','white'))(100),
              #color=colorRampPalette(viridis)(100),
              
              ##second candidate
              #color = viridis(n = 256, alpha = 1, 
              #                     begin = 0, end = 1, option = "viridis"),
              
              #annotation_row = annotdf,
              #annotation_colors = mycolors,
              
              #annotation_col = annotdf_col,
              #annotation_colors = mycolors_col,
              #color=colorRampPalette(c("lightgoldenrodyellow", "firebrick3"))(2)
              #color = inferno(length(mat_breaks) - 1),
              #breaks = mat_breaks,
              
              ##open it when necessary
              #display_numbers = matrix(ifelse(signal_enrichment_matrix < 0.05, "*", ""), nrow(signal_enrichment_matrix)),
              fontsize_number = 25
) 
print(p)
dev.off()


