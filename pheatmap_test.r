#设置操作目录
getwd()
setwd("/home/zj/铁代谢")
dir.create("all")
dir.create("all/2_26")
save_folders <- "/home/zj/铁代谢/all/2_26"

#导入汇总文件
merge_salmon_file <- read.csv("merge_salmon_file.csv",header = T)
head(merge_salmon_file)
merge_salmon_file <- merge_salmon_file[-c(50:55)]
head(merge_salmon_file)

#TPM基因表达热图

#数据处理
TPM_salmon_file <- merge_salmon_file

rownames(TPM_salmon_file) <- TPM_salmon_file$Name

TPM_salmon_file <- TPM_salmon_file[,-c(1)]


JBC_gene_list <- read.csv("JBC_gene_list.csv",header = T)
head(JBC_gene_list)

JBC_gene_list <- JBC_gene_list[c("JBC.gene.list","Aliases")]

colnames(JBC_gene_list) <- c("gene","gene_id")


TF_gene_list <- read.csv("TF_GENE_LIST_S288C.csv",header = T)
head(TF_gene_list)



TF_all_gene_filted <- TPM_salmon_file[c(JBC_gene_list$gene_id),]

TF_all_gene_filted <- na.omit(TF_all_gene_filted)
dim(TF_all_gene_filted)
head(TF_all_gene_filted)

TF_all_gene_filted$AFT1_0h <- apply(TF_all_gene_filted[,1:3], 1, mean,na.rm=T)
TF_all_gene_filted$AFT1_4h <- apply(TF_all_gene_filted[,4:6], 1, mean,na.rm=T)
TF_all_gene_filted$DOT1_0h <- apply(TF_all_gene_filted[,7:9], 1, mean,na.rm=T)
TF_all_gene_filted$DOT1_4h <- apply(TF_all_gene_filted[,10:12], 1, mean,na.rm=T)
TF_all_gene_filted$GCN5_0h <- apply(TF_all_gene_filted[,13:15], 1, mean,na.rm=T)
TF_all_gene_filted$GCN5_4h <- apply(TF_all_gene_filted[,16:18], 1, mean,na.rm=T)
TF_all_gene_filted$RTT109_0h <- apply(TF_all_gene_filted[,19:21], 1, mean,na.rm=T)
TF_all_gene_filted$RTT109_4h <- apply(TF_all_gene_filted[,22:24], 1, mean,na.rm=T)
TF_all_gene_filted$SAS2_0h <- apply(TF_all_gene_filted[,25:27], 1, mean,na.rm=T)
TF_all_gene_filted$SAS2_4h <- apply(TF_all_gene_filted[,28:30], 1, mean,na.rm=T)
TF_all_gene_filted$SET1_0h <- apply(TF_all_gene_filted[,31:33], 1, mean,na.rm=T)
TF_all_gene_filted$SET1_4h <- apply(TF_all_gene_filted[,34:36], 1, mean,na.rm=T)
TF_all_gene_filted$SET2_0h <- apply(TF_all_gene_filted[,37:39], 1, mean,na.rm=T)
TF_all_gene_filted$SET2_4h <- apply(TF_all_gene_filted[,40:42], 1, mean,na.rm=T)
TF_all_gene_filted$WT_0h <- apply(TF_all_gene_filted[,43:45], 1, mean,na.rm=T)
TF_all_gene_filted$WT_4h <- apply(TF_all_gene_filted[,46:48], 1, mean,na.rm=T)
TF_all_gene_filted$YNG2_0h <- apply(TF_all_gene_filted[,49:51], 1, mean,na.rm=T)
TF_all_gene_filted$YNG2_4h <- apply(TF_all_gene_filted[,52:54], 1, mean,na.rm=T)

TF_all_gene_filted <- TF_all_gene_filted[,-c(1:54)]


#热图绘制
#准备注释文件
ann_col <- data.frame(Group = factor(c(rep("AFT1",2),rep("DOT1",2),rep("GCN5",2),rep("RTT109",2),rep("SAS2",2)
                                       ,rep("SET1",2),rep("SET2",2),rep("WT",2),rep("YNG2",2))),
                      Time =factor(c(rep(c("0h","4h"),9))))
rownames(ann_col) <- colnames(TF_all_gene_filted)

TF_all_gene_filted$gene_id <- rownames(TF_all_gene_filted)
head(TF_all_gene_filted)

head(JBC_gene_list)


deal_TF_all_gene_filted <- merge(as.data.frame(TF_all_gene_filted),
                                 as.data.frame(JBC_gene_list),by='gene_id',sort=FALSE)

head(deal_TF_all_gene_filted)

deal_TF_all_gene_filted <- deal_TF_all_gene_filted[,-c(1)]

rownames(deal_TF_all_gene_filted) <- deal_TF_all_gene_filted$gene

dim(deal_TF_all_gene_filted)

deal_TF_all_gene_filted <- deal_TF_all_gene_filted[,-c(19)]
deal_TF_all_gene_filted


deal_TF_all_gene_filted_1 <- deal_TF_all_gene_filted[c("WT_0h","WT_4h","AFT1_0h",
                                                       "AFT1_4h","GCN5_0h","GCN5_4h","RTT109_0h",
                                                       "RTT109_4h","SAS2_0h","SAS2_4h","YNG2_0h",
                                                       "YNG2_4h","SET1_0h","SET1_4h","SET2_0h",
                                                       "SET2_4h","DOT1_0h","DOT1_4h")]



deal_TF_all_gene_filted_1$WT_4h <- log2(deal_TF_all_gene_filted_1$WT_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$AFT1_0h <- log2(deal_TF_all_gene_filted_1$AFT1_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$AFT1_4h <- log2(deal_TF_all_gene_filted_1$AFT1_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$GCN5_0h <- log2(deal_TF_all_gene_filted_1$GCN5_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$GCN5_4h <- log2(deal_TF_all_gene_filted_1$GCN5_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$RTT109_0h <- log2(deal_TF_all_gene_filted_1$RTT109_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$RTT109_4h <- log2(deal_TF_all_gene_filted_1$RTT109_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SAS2_0h <- log2(deal_TF_all_gene_filted_1$SAS2_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SAS2_4h <- log2(deal_TF_all_gene_filted_1$SAS2_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$YNG2_0h <- log2(deal_TF_all_gene_filted_1$YNG2_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$YNG2_4h <- log2(deal_TF_all_gene_filted_1$YNG2_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SET1_0h <- log2(deal_TF_all_gene_filted_1$SET1_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SET1_4h <- log2(deal_TF_all_gene_filted_1$SET1_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SET2_0h <- log2(deal_TF_all_gene_filted_1$SET2_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$SET2_4h <- log2(deal_TF_all_gene_filted_1$SET2_4h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$DOT1_0h <- log2(deal_TF_all_gene_filted_1$DOT1_0h/deal_TF_all_gene_filted_1$WT_0h)
deal_TF_all_gene_filted_1$DOT1_4h <- log2(deal_TF_all_gene_filted_1$DOT1_4h/deal_TF_all_gene_filted_1$WT_0h)

head(deal_TF_all_gene_filted_1)

deal_TF_all_gene_filted_1 <- deal_TF_all_gene_filted_1[,-c(1)]

library(pheatmap)


bk <- c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))

pheatmap(deal_TF_all_gene_filted_1,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,   
         scale = "none",
         clustering_method = "average",
         cluster_row = F,
         cluster_col = F,
         border_size = NA,
         border_color = "white",
         main = "TF_gene_TPM_pheatmap_JBC_fe_gene_list_NO_cluster",
         treeheight_row = 30,
         treeheight_col = 30,
         cellwidth = 15, 
         cellheight = 10,
         show_rownames = T,
         display_numbers = T,
         fontsize = 5,
         fontsize_row = 5,
         fontsize_col = 7,
         angle_col = 45,
         annotation_col = ann_col,
         filename = file.path(save_folders, "TF_gene_TPM_pheatmap_JBC_NO_cluster_numbers_log2.pdf"))

pheatmap(deal_TF_all_gene_filted_1,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk, 
         scale = "none",
         clustering_method = "average",
         cluster_row = T,
         cluster_col = F,
         border_size = NA,
         border_color = "white",
         main = "TF_gene_TPM_pheatmap_JBC_fe_gene_list_NO_cluster",
         treeheight_row = 30,
         treeheight_col = 30,
         cellwidth = 15, 
         cellheight = 10,
         show_rownames = T,
         display_numbers = T,
         fontsize = 5,
         fontsize_row = 5,
         fontsize_col = 7,
         angle_col = 45,
         annotation_col = ann_col,
         filename = file.path(save_folders, "TF_gene_TPM_pheatmap_JBC_cluster_numbers_log2.pdf"))
