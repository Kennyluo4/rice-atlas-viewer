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
input_0_users_provide_species_name <- "Osativa"
#Os Zm

input_1_users_select_cate <- as.character(arg[2])
# input_1_users_select_cate <- "Gene"
##Gene TF Mt ACR

input_2_users_select_path_name <- as.character(arg[3])
# input_2_users_select_path_name <- "bud_ProcamPhloem"
##Eseedling_SAMPhloem

input_3_users_provide_ID <- as.character(arg[4])
# input_3_users_provide_ID <- "LOC_Os04g12080"

# input_sorted_pt_dir <- as.character(arg[5])
input_sorted_pt_dir <-  'Rplot/data/input_sorted_pt_dir'
##it contains all the sorted pt files used for the plotting

# input_annot_TF_dir <- as.character(arg[6])
input_annot_TF_dir <- 'Rplot/data/input_annot_TF_dir'
##it contains all species TF annotation to show their super family information

# input_gene_TF_motif_impute_sparse_rds_dir <- as.character(arg[7])
input_gene_TF_motif_impute_sparse_rds_dir <- 'Rplot/data/input_gene_TF_motif_impute_sparse_rds_dir'
##it contains the impute sparse rds fl for all genes as well as motifs

# input_traject_meta_fl_dir <- as.character(arg[8])
input_traject_meta_fl_dir <- 'Rplot/data/input_traject_meta_fl_dir'

input_output_dir <- as.character(arg[5])
input_output_dir <- "Rplot/temp"


#################################################
##step02 plot the UMAP information for target IDs
##the step02 is not available for the ACR
if (input_1_users_select_cate != 'ACR'){

  if (input_1_users_select_cate == 'TF'){
    all_gene_impute_acc_sparse_rds_fl <- paste0(input_gene_TF_motif_impute_sparse_rds_dir,'/Gene/',input_2_users_select_path_name,'_Gene.rds')
  }
  if (input_1_users_select_cate == 'Gene'){
    all_gene_impute_acc_sparse_rds_fl <- paste0(input_gene_TF_motif_impute_sparse_rds_dir,'/Gene/',input_2_users_select_path_name,'_Gene.rds')
  }
  if (input_1_users_select_cate == 'Mt'){
    all_gene_impute_acc_sparse_rds_fl <- paste0(input_gene_TF_motif_impute_sparse_rds_dir,'/Mt/',input_2_users_select_path_name,'_Mt.rds')
  }
  
  meta_fl <- paste0(input_traject_meta_fl_dir,'/',input_2_users_select_path_name,'.trajectory.txt')
  
  markers <- input_3_users_provide_ID
  
  
  plot.act.scores    <- function(meta_fl,output_dir,
                                 act_fl=all_gene_impute_acc_sparse_rds_fl, 
                                 info_fl=NULL, 
                                 top=NULL,
                                 logT=F,
                                 marker.dist=NULL,
                                 outname="markerActivityScores", 
                                 lim=0.95)
    {
    
    
    acts <- readRDS(act_fl)
    
    ##read the meta file
    df <- read.delim(meta_fl,row.names = 1)
    
    ##read the marker
    #info <- read.delim(info_fl,header=F)
    
    ##corresponding the row and col of motifs, acc, and meta dt
    intersect_cells <- intersect(rownames(df),colnames(acts))
    df <- df[intersect_cells,]
    acts <- acts[,intersect_cells]
    
    # prep data
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,which(rownames(df) %in% colnames(acts))]
    
    # reorder rows
    #rownames(info) <- info$V1
    #info <- info[order(info$type),]
    #info.genes <- rownames(info)
    
    info.genes <- c(info_fl)
    
    act.genes <- rownames(acts)
    rd.cells <- rownames(df)
    
    #get the common genes
  
    ##we will check the motif name
    if (input_1_users_select_cate == 'Mt'){
      ##for the Mt we will modify the ID
      all_ID_motifnm_list <- rownames(acts)
      all_ID_motifnm_dt <- as.data.frame(all_ID_motifnm_list)
      all_ID_motifnm_dt$ID <- gsub('_.+','',all_ID_motifnm_dt$all_ID_motifnm_list)
      all_ID_motifnm_dt$IDconcise <- gsub('\\..+','',all_ID_motifnm_dt$ID)
      all_ID_motifnm_dt$motif <- gsub('.+_','',all_ID_motifnm_dt$all_ID_motifnm_list)
      head(all_ID_motifnm_dt)
    
      all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$motif %in% info.genes,]
      if (nrow(all_ID_motifnm_dt_flt) == 0){
        
        all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$IDconcise %in% info.genes,]
        
        if (nrow(all_ID_motifnm_dt_flt) == 0){
          all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$ID %in% info.genes,]
          
          if (nrow(all_ID_motifnm_dt_flt) == 0){
            
            all_ID_motifnm_dt_flt <- all_ID_motifnm_dt[all_ID_motifnm_dt$all_ID_motifnm_list %in% info.genes,]
            
          }
          
        }
        
      }
      
      target_rownames <- all_ID_motifnm_dt_flt$all_ID_motifnm_list
      common <- target_rownames
      
    ##if the target is not the Mt we will directly get the common gene
    }else{

      common <- intersect(info.genes, act.genes)
      #info <- info[which(rownames(info) %in% common),]
    }
      
      
    info.ordered <- common
    
    sub.scores <- as.matrix(acts[info.ordered,])
    colnames(sub.scores) <- info.ordered
    sub.scores <- t(sub.scores)
    
    gids <- info.ordered
    
    # setup plot size
    #nrows <- ceiling(length(gids)/3)
    #totals <- nrows*3
    #ratio <- nrows/3
    
    
    # params
    #png(file=paste0(output_dir,'/',outname), width=6, height=ratio*6, units="in", res=300, type="cairo")
    
    #pdf(file=paste0(output_dir,'/',outname), width=6, height=ratio*6)
    
   
    
    
    #if (length(gids) >= 6){
    #layout(matrix(c(1:totals), ncol=3, byrow=T))
    #}else{
    #  layout(matrix(c(1:length(gids)), ncol=length(gids), byrow=T))
    #}
    
   
    #par(mar=c(2,2,1,1))
    
    # adjust cluster IDs
    message("begin plotting pre-defined markers...")
    for (i in 1:length(gids)){
      
      # copy meta data
      gene.index <- which(rownames(sub.scores) == gids[i])
      acv <- sub.scores[gene.index,]
      
      acv <- rescale(acv, c(-1, 1))
      
      # set up plot cols/sizes
      orderRow <- order(acv, decreasing=F)
      #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
      #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
      #cols <- inferno(100)
      #cols <- plasma(100)
      cols <- colorRampPalette(rev(c(brewer.pal(11, "RdYlGn")[2:11])), bias=0.7)(100)
      #cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
      acv <- as.numeric(acv[orderRow])
      if(logT==T){
        acv <- log2(acv+1)
      }
      
      ##Note: important here, we will allow the df2 to have the same order as the acv that has acc value from small to large
      df2 <- df[orderRow,]
      ##change na to be -1 other than 0
      acv[is.na(acv)] <- -1
      acv[is.infinite(acv)] <- -1
      #upper.lim <- quantile(acv, lim)
      #acv[acv > upper.lim] <- upper.lim
      if(!is.null(marker.dist)){
        message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
        colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
      }else{
        min.acv <- min(acv) - (1e-6*min(acv))
        max.acv <- max(acv) + (1e-6*max(acv))
        message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
        if(min.acv == max.acv){
          next
        }
        
        ##Note: for each color it corresponds to one range of color: like color 1 (-1,-0.98]  and color 2 (-0.98,-0.96]
        ##Note: it allows the acv to be cut into different range, and allows each range to correspond one color
        ##Note: here we allow the acc value within different range 
        colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
      }
      colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
      colvec[is.na(colvec) & acv == -1] <- cols[1]
      ##change 0 to be -1
      #sizes <- rescale(acv, c(0.25, 0.3))
      
      
      par(mar = c(7, 7, 7, 7))
      # pdf(paste0(input_output_dir,'/',input_2_users_select_path_name,".trajectory",input_1_users_select_cate,".pdf"), width=8, height=8)
      fig_name <- paste0(output_dir,'/',outname)
      message('Output figure: ',fig_name)
      tiff(fig_name, width =24, height = 30,res = 200, units = 'cm')
      # pdf(fig_name)
      par(cex = 1.5)
      

      # plot
      plot(df2[[umap1_colnm]], df2[[umap2_colnm]], col=colvec,
           #main=paste(info$name[i],info$type[i],sep="-"),
           main=paste(gids[i]),
           #main=info$name[i],
           cex.main=1,
           xlab="", ylab="", bty="n",
           xaxt="n", yaxt="n", pch=16, cex=0.25)
      
      ##Note: we will generate a bar that shows exact 100 colors other than the range exactly corresponding to the real plotting.
      add.color.bar(3, cols, title='low -> high',subtitle='', lims=NULL,prompt=F,x=min(df2[[umap1_colnm]]+8),y=min(df2[[umap2_colnm]])+1,cex=0.15)
      
      # turn device off
      dev.off()
    }
    
   
    
  }
  
  lim <- 0.98
  umap1_colnm <- 'umapPT_1'
  umap2_colnm <- 'umapPT_2'
  plot.act.scores(meta_fl,input_output_dir,
                  act_fl=all_gene_impute_acc_sparse_rds_fl,
                  info_fl=markers,
                  logT=F,
                  lim=lim,
                  marker.dist=NULL,
                  outname=paste0(input_2_users_select_path_name,".trajectory.",input_1_users_select_cate,"_UMAP.tiff"))
  
 
}

















