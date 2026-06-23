
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas")
}
library(RColorBrewer)
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidydr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(survcomp)
library(customLayout)
#library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('z:/projects/codes/mg_base.R')


#TCGA############
dir.create('results/00.pre_data')
tcga.pheno.all=read.delim('origin_datas/TCGA/TCGA.LIHC.sampleMap_LIHC_clinicalMatrix')
colnames(tcga.pheno.all)
table(tcga.pheno.all$X_primary_disease)
table(tcga.pheno.all$histological_type)
tcga.pheno.all=tcga.pheno.all %>%subset(histological_type=='Hepatocellular Carcinoma')

tcga.pheno=data.frame(Samples=tcga.pheno.all$sampleID,
                      Age=tcga.pheno.all$age_at_initial_pathologic_diagnosis,
                      Gender=tcga.pheno.all$gender,
                      tcga.pheno.all[,c('pathologic_T','pathologic_N','pathologic_M','pathologic_stage')],
                      neoplasm_histologic_grade=tcga.pheno.all$neoplasm_histologic_grade,
                      radiation_therapy=tcga.pheno.all$radiation_therapy )
table(tcga.pheno$pathologic_T)
tcga.pheno$pathologic_T[tcga.pheno$pathologic_T %in% c('','[Discrepancy]','TX')]=NA
tcga.pheno$pathologic_T=gsub('[abc]','',tcga.pheno$pathologic_T)

table(tcga.pheno$pathologic_N)
tcga.pheno$pathologic_N[tcga.pheno$pathologic_N %in% c('','NX')]=NA

table(tcga.pheno$pathologic_M)
tcga.pheno$pathologic_M[tcga.pheno$pathologic_M %in% c('','MX')]=NA

table(tcga.pheno$pathologic_stage)
tcga.pheno$pathologic_stage[tcga.pheno$pathologic_stage %in% c('','[Discrepancy]')]=NA
tcga.pheno$pathologic_stage=gsub('[ABC]','',tcga.pheno$pathologic_stage)

table(tcga.pheno$neoplasm_histologic_grade)
tcga.pheno$neoplasm_histologic_grade[tcga.pheno$neoplasm_histologic_grade %in% c('')]=NA

table(tcga.pheno$radiation_therapy)
tcga.pheno$radiation_therapy[tcga.pheno$radiation_therapy %in% c('','[Discrepancy]')]=NA


tcga.survival=read.delim('origin_datas/TCGA/survival_LIHC_survival.txt')
head(tcga.survival)
tcga.survival=tcga.survival[,c(1,3:10)]
colnames(tcga.survival)[1]='Samples'
tcga.survival$OS.time
tcga.survival=tcga.survival %>% drop_na(OS.time) %>%subset(OS.time>30)
dim(tcga.survival)
#407

tcga.cli=merge(tcga.pheno,tcga.survival,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
dim(tcga.cli)
#399
head(tcga.cli)


genecode=read.delim('z:/users/project/public/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
head(genecode)
PCD_SYMBOL=genecode[which(genecode$TYPE=='protein_coding'),'SYMBOL']


tcga.mat=read.delim('origin_datas/TCGA/TCGA.LIHC.sampleMap_HiSeqV2.gz',row.names = 1,check.names = F)
tcga.mat[1:5,1:5]
dim(tcga.mat)
tcga.mat=tcga.mat[rownames(tcga.mat) %in% PCD_SYMBOL,]
range(tcga.mat)
table(substr(colnames(tcga.mat),14,15))

tcga.sample.T=colnames(tcga.mat)[which(as.numeric(substr(colnames(tcga.mat),14,15))==1)]
tcga.sample.N=colnames(tcga.mat)[which(as.numeric(substr(colnames(tcga.mat),14,15))==11)]样本
tcga.type=data.frame(Samples=c(tcga.sample.T,tcga.sample.N),Type=rep(c('HCC','Adjacent'),c(length(tcga.sample.T),length(tcga.sample.N))))
rownames(tcga.type)=tcga.type$Samples
table(tcga.type$Type)

tcga.sample.T=intersect(tcga.sample.T,tcga.cli$Samples)
tcga.cli=tcga.cli[tcga.sample.T,]
tcga.exp=tcga.mat[,tcga.sample.T]
dim(tcga.cli);dim(tcga.exp)
# [1] 335  16
# [1] 16788   335

save(tcga.mat,tcga.exp,tcga.cli,tcga.type,file = 'results/00.pre_data/TCGA_data.RData')



#ICGC############
icgc.mat=read.delim('origin_datas/ICGC/HCCDB18_mRNA_level3.txt',row.names = 2,check.names = F)
icgc.mat[1:5,1:5]
icgc.mat=icgc.mat[rownames(icgc.mat) %in% PCD_SYMBOL,-1]
range(icgc.mat)
dim(icgc.mat)

icgc.patient=read.delim('origin_datas/ICGC/HCCDB18.patient.txt',check.names = F,row.names = 1)
icgc.patient=t(icgc.patient)
icgc.patient=as.data.frame(icgc.patient)
head(icgc.patient)
dim(icgc.patient)

icgc.sample=read.delim('origin_datas/ICGC/HCCDB18.sample.txt',check.names = F,row.names = 1)
icgc.sample=t(icgc.sample)
icgc.sample=as.data.frame(icgc.sample)
head(icgc.sample)
icgc.sample$SAMPLE_ID=rownames(icgc.sample)
dim(icgc.sample)

icgc.cli=merge(icgc.patient,icgc.sample,by='PATIENT',all = FALSE)
head(icgc.cli)
rownames(icgc.cli)=icgc.cli$SAMPLE_ID
table(icgc.cli$TYPE)

icgc.sample.T=intersect(colnames(icgc.mat),icgc.cli$SAMPLE_ID[icgc.cli$TYPE=='HCC'])
icgc.sample.N=intersect(colnames(icgc.mat),icgc.cli$SAMPLE_ID[icgc.cli$TYPE=='Adjacent'])

icgc.type=data.frame(Samples=c(icgc.sample.T,icgc.sample.N),Type=rep(c('HCC','Adjacent'),c(length(icgc.sample.T),length(icgc.sample.N))))
rownames(icgc.type)=icgc.type$Samples
table(icgc.type$Type)


icgc.cli=icgc.cli[icgc.sample.T,]
icgc.exp=icgc.mat[,icgc.sample.T]
dim(icgc.cli);dim(icgc.exp)
# [1] 212  22
# [1] 17055   212

save(icgc.mat,icgc.exp,icgc.cli,icgc.type,file = 'results/00.pre_data/ICGC_data.RData')


#01.TCMSP###################
dir.create('results/01.TCMSP')


library(readxl)
library(dplyr)
file_path <- "results/01.TCMSP/_TCMSP.xlsx"
sheet_names <- excel_sheets(file_path)
herb_TCMSP <- lapply(sheet_names, function(s) {
  df=read_excel(file_path, sheet = s) %>% mutate(herb=s)
}) %>%  dplyr::bind_rows()
 

head(herb_TCMSP)
herb_TCMSP_filter=herb_TCMSP %>% subset(`OB(%)`>30 & DL>0.18)
herb_TCMSP_filter=herb_TCMSP_filter[!duplicated(herb_TCMSP_filter$`Mol ID`),]
dim(herb_TCMSP_filter)
table(herb_TCMSP_filter$herb)
unique(herb_TCMSP_filter$`Mol ID`)
write.csv(herb_TCMSP_filter,'results/01.TCMSP/Table_1.csv')


file_path2 <- "results/01.TCMSP/_TCMSP_targets.xlsx"
target_TCMSP <- lapply(sheet_names, function(s) {
  df=read_excel(file_path2, sheet = s) 
}) %>%  dplyr::bind_rows()
head(target_TCMSP)
target_TCMSP_filter=target_TCMSP[target_TCMSP$`Mol ID` %in% unique(herb_TCMSP_filter$`Mol ID`),] %>% 
  subset(`Target name` !='')
dim(target_TCMSP_filter)
# write.xlsx(target_TCMSP_filter,file = 'results/111.xlsx')

herb_targets=read.xlsx('results/01.TCMSP/converted_gene_symbol.xlsx')
herb_targets=herb_targets %>% drop_na(Gene_Symbol)
herb_targets=herb_targets$Gene_Symbol
length(herb_targets)
#257

##network#########
herb_targets2=read.xlsx('results/01.TCMSP/converted_gene_symbol.xlsx')
herb_targets2=herb_targets2 %>% drop_na(Gene_Symbol)
head(herb_targets2)

head(herb_TCMSP_filter)
TCMSP_filter=merge(herb_TCMSP_filter[,c(1,13)],target_TCMSP_filter,by='Mol ID')
head(TCMSP_filter)
write.csv(TCMSP_filter,'results/01.TCMSP/TCMSP_filter.csv')
TCMSP_filter=as.data.frame(TCMSP_filter)

node.net=merge(TCMSP_filter,herb_targets2,by.x='Target name',by.y='gene_full_name')
node.net=node.net[,c('herb','Mol ID','Gene_Symbol')]
head(node.net)

length(unique(node.net$herb))#9
length(unique(node.net$`Mol ID`))#143
length(unique(node.net$Gene_Symbol))#255


node.table <- rbind(
  data.frame(node1 = node.net$herb, node2 = node.net$`Mol ID`),
  data.frame(node1 = node.net$Gene_Symbol, node2 = node.net$`Mol ID`)
)
node.table=node.table %>% distinct()

node.table.anno=data.frame(node=c(unique(node.net$herb),
                                   unique(node.net$`Mol ID`),
                                   unique(node.net$Gene_Symbol)),
                           type=rep(c('Herb','Molecule','Target'),c(9,143,255)))
write.table(node.table,'results/01.TCMSP/node.table.txt',quote = F,sep = '\t',row.names = F)
write.table(node.table.anno,'results/01.TCMSP/node.table.anno.txt',quote = F,sep = '\t',row.names = F)


library(clusterProfiler)
library(org.Hs.eg.db)
target.entrez_id = mapIds(x = org.Hs.eg.db,
                          keys = herb_targets,
                          keytype = "SYMBOL",
                          column = "ENTREZID")
target.entrez_id = na.omit(target.entrez_id)

target.enrichKEGG=enrichKEGG(target.entrez_id,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.25,
                              pAdjustMethod = "BH",
                              minGSSize = 1,
                              maxGSSize = 500)
kegg.plot=enrichplot::dotplot(target.enrichKEGG)+ggtitle('KEGG')+
  theme(text = element_text(size = 12,color='black'),
        axis.text.x = element_text(size = 12,color='black'),
        axis.text.y = element_text(size = 12,color='black'))
kegg.plot

target.enrichKEGG.res=target.enrichKEGG@result
write.xlsx(target.enrichKEGG.res,'results/01.TCMSP/KEGG_enrichment.xlsx',overwrite = T)

target.enrichGO=enrichGO(target.entrez_id,
                       OrgDb = "org.Hs.eg.db", 
                       ont="ALL", pvalueCutoff=0.05)

target.enrichGO.res=target.enrichGO@result
write.xlsx(target.enrichGO.res,'results/01.TCMSP/GO_enrichment.xlsx',overwrite = T)
target.enrichGO.res=target.enrichGO.res %>% group_by(ONTOLOGY) %>% slice_max(n =5, order_by = Count)
target.enrichGO.res$Description <- factor(target.enrichGO.res$Description, 
                                        levels = target.enrichGO.res$Description[order(target.enrichGO.res$ONTOLOGY,decreasing = T)])

GP.plot=ggplot(target.enrichGO.res, aes(x = Description, y = Count, 
                              fill = ONTOLOGY, group = Description)) +
  geom_bar(stat="identity", position="dodge", colour=aes(ONTOLOGY))+
  scale_fill_manual(values =c("#80B1D3","#FDB462","#91C133"))+
  xlab('')+ggtitle('GO enrichment')+theme_bw()+coord_flip()+
  scale_x_discrete(labels=function(y)str_wrap(y,width = 40))+
  theme(text = element_text(size = 12,color='black'),
        axis.text.x = element_text(size = 12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        legend.position = 'top')
GP.plot

enrichment.plot=mg_merge_plot(kegg.plot,GP.plot,labels = c('B','C'))
ggsave('results/01.TCMSP/enrichmnet_plot.pdf',enrichment.plot,height = 8,width = 15)



dir.create('results/02.DEGs')
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","grey"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+
    theme_bw()+
    theme(
      legend.title = element_blank(),
      legend.position = leg.pos,
      text = element_text(size=12)
    )+
    ylab(ylab)+
    xlab(xlab)+
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)padj<0.05
  return(p)
}

tcga.limma=mg_limma_DEG(exp = tcga.mat[,tcga.type$Samples],group = tcga.type$Type,
                        ulab = 'HCC',dlab = 'Adjacent')
tcga.limma$Summary
head(tcga.limma$DEG)
tcga.DEGs=tcga.limma$DEG %>% subset(abs(logFC)>1 & adj.P.Val<0.05)
write.csv(tcga.DEGs,'results/02.DEGs/tcga.DEGs.csv')

icgc.limma=mg_limma_DEG(exp = icgc.mat[,icgc.type$Samples],group = icgc.type$Type,
                        ulab = 'HCC',dlab = 'Adjacent')
icgc.limma$Summary
head(icgc.limma$DEG)
icgc.DEGs=icgc.limma$DEG %>% subset(abs(logFC)>1 & adj.P.Val<0.05)
write.csv(icgc.DEGs,'results/02.DEGs/icgc.DEGs.csv')

up.DEGs=intersect(rownames(tcga.DEGs[tcga.DEGs$logFC>0,]),
                  rownames(icgc.DEGs[icgc.DEGs$logFC>0,]))

dn.DEGs=intersect(rownames(tcga.DEGs[tcga.DEGs$logFC<0,]),
                  rownames(icgc.DEGs[icgc.DEGs$logFC<0,]))

all.DEGs=c(up.DEGs,dn.DEGs)

targets_DEGs=Reduce(intersect,list(herb_targets,all.DEGs))
length(targets_DEGs)
# 30
fig2a=my_volcano(tcga.limma)+ggtitle('TCGA')
fig2b=my_volcano(icgc.limma)+ggtitle('ICGC')

library(eulerr)
v1=list(TCGA.upDEGs=rownames(tcga.DEGs[tcga.DEGs$logFC>0,]),
        ICGC.upDEGs=rownames(icgc.DEGs[icgc.DEGs$logFC>0,]))
v2=list(TCGA.dnDEGs=rownames(tcga.DEGs[tcga.DEGs$logFC<0,]),
        ICGC.dnDEGs=rownames(icgc.DEGs[icgc.DEGs$logFC<0,]))
fig2c=plot(venn(v1),labels = list(col = "gray20", font = 2), 
                     edges = list(col="gray60", lex=1),
                     fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
                     quantities = list(cex=.8, col='gray20'))
fig2c
fig2d=plot(venn(v2),labels = list(col = "gray20", font = 2), 
           edges = list(col="gray60", lex=1),
           fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
           quantities = list(cex=.8, col='gray20'))
fig2d

figure2=mg_merge_plot(fig2a,fig2b,fig2c,fig2d,nrow=2,ncol=2,heights = c(1.2,1),
                      labels = LETTERS[1:4])
ggsave('results/02.DEGs/Figure2.pdf',figure2,height = 9,width = 12)

#03.PPI
dir.create('results/03.PPI')
# write.table(targets_DEGs,'results/03.PPI/targets_DEGs.txt',quote = F,row.names = F,col.names = F)

ppi.path='results/03.PPI/cytoHubba/'
ppi.f=list.files(ppi.path)
cytoHubba_res <- lapply(ppi.f, function(f) {
  df=read.csv(paste0(ppi.path,f),skip = 1)
  g=df$Name
}) 
names(cytoHubba_res)=str_split_fixed(ppi.f,'_',2)[,1]
names(cytoHubba_res)[c(1,3,4,7)]
hub.genes=Reduce(intersect,cytoHubba_res[c(1,3,4,7)])
hub.genes

library(ggvenn)
cytoHubba_res_venn=cytoHubba_res[c(1,3,4,7)]
pdf('results/03.PPI/cytoHubba_vennplot.pdf',height = 8,width = 8)
ggvenn(cytoHubba_res_venn,fill_alpha = .7,stroke_linetype = "solid",
       set_name_size = 5,fill_color = brewer.pal(4,"Set3"),
       text_size=3) 
dev.off()

#04.
dir.create('results/04.Hubgene_expr')

hub.gene.df1=data.frame(t(tcga.mat[hub.genes,tcga.type$Samples]),
                        Type=tcga.type$Type)
head(hub.gene.df1)
hub.gene.df1=melt(hub.gene.df1)
fig4a=hub.gene.df1 %>% 
  ggplot(aes(x=variable,y=value,fill=Type))+
  scale_fill_manual(values = c('skyblue','pink'))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif')+
  geom_boxplot(size=1,outlier.fill="white",outlier.color="white")+
  theme_bw()+labs(y='mRNA Expression',x='',title = 'TCGA')+
  theme(legend.position = 'top',
        text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))
fig4a

hub.gene.df2=data.frame(t(icgc.mat[hub.genes,icgc.type$Samples]),
                        Type=icgc.type$Type)
head(hub.gene.df2)
hub.gene.df2=melt(hub.gene.df2)
fig4b=hub.gene.df2 %>% 
  ggplot(aes(x=variable,y=value,fill=Type))+
  scale_fill_manual(values = c('skyblue','pink'))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif')+
  geom_boxplot(size=1,outlier.fill="white",outlier.color="white")+
  theme_bw()+labs(y='mRNA Expression',x='',title = 'ICGC')+
  theme(legend.position = 'top',
        text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))
fig4b


table(tcga.cli$pathologic_stage)
# Stage IV 
table(tcga.cli$neoplasm_histologic_grade)
# Stage IV 
colnames(tcga.cli)
hub.gene.df3=data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),
                       tcga.cli[,4:9])
head(hub.gene.df3)
hub.gene.df3=melt(hub.gene.df3)

fig4c=hub.gene.df3 %>% drop_na(pathologic_stage) %>% 
  ggplot(aes(x=variable,y=value,fill=pathologic_stage))+
  stat_compare_means(method = 'kruskal.test',label = 'p.signif')+
  geom_boxplot(size=1,outlier.fill="white",outlier.color="white")+
  theme_bw()+labs(y='mRNA Expression',x='pathologic_stage')+
  theme(legend.position = 'top',
        text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))
fig4c

fig4d=hub.gene.df3 %>% drop_na(neoplasm_histologic_grade) %>% 
  ggplot(aes(x=variable,y=value,fill=neoplasm_histologic_grade))+
  stat_compare_means(method = 'kruskal.test',label = 'p.signif')+
  geom_boxplot(size=1,outlier.fill="white",outlier.color="white")+
  theme_bw()+labs(y='mRNA Expression',x='neoplasm_histologic_grade')+
  theme(legend.position = 'top',
        text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))
fig4d


colnames(icgc.cli)
hub.gene.df4=data.frame(t(icgc.exp[hub.genes,icgc.cli$SAMPLE_ID]),
                        pathologic_stage=icgc.cli$TNM_STAGE_T1)
head(hub.gene.df4)
hub.gene.df4=melt(hub.gene.df4)
fig4e=hub.gene.df4 %>% drop_na(pathologic_stage) %>% 
  ggplot(aes(x=variable,y=value,fill=pathologic_stage))+
  stat_compare_means(method = 'kruskal.test',label = 'p.signif')+
  geom_boxplot(size=1,outlier.fill="white",outlier.color="white")+
  theme_bw()+labs(y='mRNA Expression',x='pathologic_stage')+
  theme(legend.position = 'top',
        text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'))
fig4e


hub.genes.cox=cox_batch(dat = tcga.exp[hub.genes,tcga.cli$Samples],
                        time = tcga.cli$OS.time,event = tcga.cli$OS)
hub.genes.cox

tcga_cox_datas=data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),
                          OS.time=tcga.cli$OS.time,
                          OS=tcga.cli$OS)
head(tcga_cox_datas)
library(ezcox)
uni_ezcox <- ezcox(data = tcga_cox_datas, covariates = hub.genes,
                   time = 'OS.time',status = 'OS')
uni_ezcox


fig4f=show_forest(tcga_cox_datas, covariates = hub.genes,time = 'OS.time',status = 'OS',
               add_caption = T)
fig4f

fig4=mg_merge_plot(fig4a,fig4c,fig4d,fig4f,nrow=2,ncol=2,labels = LETTERS[1:4])
ggsave('results/04.Hubgene_expr/Figure4.pdf',fig4,height = 9,width = 12)


dir.create('results/05.Immune_cor')
library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(magrittr)
library(reshape)
library(psych)
library(devtools)
# install_github('Hy4m/linkET',force = TRUE) 
library(linkET)

##estimate########
tcga.estimate=read.delim('results/TCGA_ESTIMATE_score.txt',row.names = 1,check.names = F)
head(tcga.estimate)
tcga.estimate.df=data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),
                            tcga.estimate[tcga.cli$Samples,])

table1=tcga.estimate.df[,hub.genes]
table2=tcga.estimate.df[,colnames(tcga.estimate)]


cor.df=list()
for (i in 1:length(hub.genes)) {
  # i=1
  cr=psych::corr.test(x=table1[,hub.genes[i]],
                      y=table2
                      ,method = 'spearman')
  cor.df[[i]]=t(rbind(cr$r,cr$p))
  colnames(cor.df[[i]])=c('r','p.value')
  cor.df[[i]]=data.frame(gene=hub.genes[i],cell=colnames(table2),cor.df[[i]])
  cor.df[[i]]
  
}
df=do.call(rbind, cor.df)
head(df)
range(df$r)
mantel <- df %>%
  mutate(correlation = cut(r, breaks = c(-1, 0, 1),
                           labels = c("Neg", "Pos")),
         pval = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                    labels = c("< 0.01", "< 0.05", ">= 0.05"),
                    right = FALSE, include.lowest = TRUE),
         cor= cut(abs(r), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                  labels = c("<0.2", "0.2<= r <0.4","0.4<= r <0.6","0.6<= r <0.8")))
mantel


pdf('results/05.Immune_cor/estimate_cor.pdf',height = 7,width = 7)
pal = c("orange","lightslateblue","lightgray")
qcorrplot(correlate(table2), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(data = mantel,aes(colour = pval,size=cor, linetype = correlation),curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
  scale_size_manual(values = c(0.5,1,1.5,2)) +
  scale_colour_manual(values = pal)+ scale_linetype_manual(values = c("dotted", "solid"))
dev.off()

##cibersort##########
tcga.cibersort=read.delim('results/TCGA_CIBERSORT_Results.txt',row.names = 1,check.names = F)
head(tcga.cibersort)
tcga.cibersort.df=data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),
                            tcga.cibersort[tcga.cli$Samples,1:22],check.names = F)
head(tcga.cibersort.df)


cor_res <- Hmisc::rcorr(as.matrix(tcga.cibersort.df),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0
library(corrplot)
pdf('results/05.Immune_cor/cibersort_cor.pdf',height = 5,width = 7)
corrplot(as.matrix(cor_res$r[hub.genes,colnames(tcga.cibersort)[1:22]]),
         p.mat = as.matrix(cor_res$P[hub.genes,colnames(tcga.cibersort)[1:22]]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('blue', 'white','red'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[2],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()

#06.GSEA
dir.create('results/06.GSEA')
cutoff.value=sapply(hub.genes, function(s) {
  cutoff<-surv_cutpoint(tcga_cox_datas,time="OS.time",event="OS",variables=s)
  cutoff$cutpoint$cutpoint
})
cutoff.value

getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  

  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  return(geneList)
}

hall.gmt=read.gmt('z:/users/project/public/gmt/h.all.v7.5.1.entrez.gmt')
source('z:/users/project/R_packages/GseaVis/R/dotplotGsea.R')
p.gsea=list()
for (i in 1:5) {
  tcga_cox_datas$group=ifelse(tcga_cox_datas[,hub.genes[i]]>unname(cutoff.value[hub.genes[i]]),'High','Low')
  
  tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.cli$Samples],
                          group=tcga_cox_datas$group
                          ,ulab='High',dlab = 'Low')
  set.seed(i)
  tcga.gsea<-GSEA(tcga.geneList,TERM2GENE = hall.gmt,seed=T)
  tcga.gsea.res=tcga.gsea@result
  write.csv(tcga.gsea.res,paste0('results/06.GSEA/',hub.genes[i],'_GSEA_res.csv'))
  tcga.gsea.res$group=ifelse(tcga.gsea.res$NES>0,"Activated","Suppressed")
  # topPathways = tcga.gsea.res%>% group_by(group)%>% slice_max(n =3,order_by = abs(NES))
  
  
  gsea.dotplot=dotplotGsea(data = tcga.gsea,topn = 5,order.by = 'NES')
  p.gsea[[i]]=gsea.dotplot$plot+theme(text = element_text(size=12))+
    scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))+ggtitle(hub.genes[i])
  # dev.off()
  
}

gsea_plot=mg_merge_plot(p.gsea,labels = LETTERS[1:5],nrow = 3,ncol=2,align = 'hv')
ggsave('results/06.GSEA/GSEA_plot.pdf',gsea_plot,height = 15,width = 20)



library(progeny)
tcga.pathway.activ=progeny(as.matrix(tcga.exp),scale = T)
dim(tcga.pathway.activ)
range(tcga.pathway.activ)
tcga.pathway.activ=as.data.frame(tcga.pathway.activ)
head(tcga.pathway.activ)

tcga.pathway.df=data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),
                            tcga.pathway.activ[tcga.cli$Samples,],
                           check.names = F)

table1=tcga.pathway.df[,hub.genes]
table2=tcga.pathway.df[,colnames(tcga.pathway.activ)]


cor.df=list()
for (i in 1:length(hub.genes)) {
  # i=1
  cr=psych::corr.test(x=table1[,hub.genes[i]],
                      y=table2
                      ,method = 'spearman')
  cor.df[[i]]=t(rbind(cr$r,cr$p))
  colnames(cor.df[[i]])=c('r','p.value')
  cor.df[[i]]=data.frame(gene=hub.genes[i],cell=colnames(table2),cor.df[[i]])
  cor.df[[i]]
  
}
df=do.call(rbind, cor.df)
head(df)
range(df$r)
mantel <- df %>%
  mutate(correlation = cut(r, breaks = c(-1, 0, 1),
                           labels = c("Neg", "Pos")),
         pval = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                    labels = c("< 0.01", "< 0.05", ">= 0.05"),
                    right = FALSE, include.lowest = TRUE),
         cor= cut(abs(r), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                  labels = c("<0.2", "0.2<= r <0.4","0.4<= r <0.6","0.6<= r <0.8")))
mantel


pdf('results/06.GSEA/pathway_cor.pdf',height = 15,width =15)
pal = c("orange","lightslateblue","lightgray")
qcorrplot(correlate(table2), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(data = mantel,aes(colour = pval,size=cor, linetype = correlation),curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
  scale_size_manual(values = c(0.5,1,1.5,2)) +
  scale_colour_manual(values = pal)+ scale_linetype_manual(values = c("dotted", "solid"))
dev.off()

save.image(file = 'project.RData')

#07.
dir.create('results/07.scRNA')
library(Seurat)
f_name=list.files('origin_datas/scRNA_GSE166635/')
datalist=list()
for (i in 1:length(f_name)){
  dir.10x = paste0("origin_datas/scRNA_GSE166635/",f_name[i])
  list.files(dir.10x)
  sce <- Read10X(data.dir = dir.10x)
  datalist[[i]]=CreateSeuratObject(counts = sce, project = f_name[i], min.cells = 3, min.features = 200)
  datalist[[i]]$Samples=f_name[i]
}
rm(sce)
names(datalist)=f_name


sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")

raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  24832
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
ggsave('results/07.scRNA/pearplot_befor.pdf',pearplot_befor,height = 5,width = 10)


sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<8000  & percent.mt<25)
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#20806
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/07.scRNA/pearplot_after.pdf',pearplot_after,height = 5,width = 10)



# options(future.globals.maxSize = 1024^5) 
sce = SCTransform(sce, vars.to.regress = "percent.mt", verbose = F) 

sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)

library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples", max.iter.harmony=50, lambda=0.5,assay.use = "SCT")

pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(size = 12))
pca.plot
ggsave('results/07.scRNA/pca.plot.pdf',pca.plot,height = 5,width = 5)




sce <- RunTSNE(sce, dims=1:20, reduction="harmony")
DimPlot(sce,group.by='Samples',reduction="tsne",label = F,pt.size = .5)

library(clustree) 
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")

sce <- FindClusters(object = sce,resolution = 0.2)
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))

seurat_clusters_tsne=DimPlot(sce,group.by='seurat_clusters',reduction="tsne",
                             label = T,pt.size = 1)+
  tidydr::theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(size = 12),legend.position = 'none')
seurat_clusters_tsne

ggsave('results/07.scRNA/seurat_clusters_tsne.pdf',seurat_clusters_tsne,height = 5,width = 5)

Idents(sce)<-'seurat_clusters'
All_DEmarkers <- FindAllMarkers(sce, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.csv(All_DEmarkers,"results/07.scRNA/All_cluster_DEG.csv")

# Hepatocytes  2
VlnPlot(sce,features = c('ANGPTL3','LBP','FGA'),pt.size = 0,group.by = 'seurat_clusters')
#Epithelial cells 2,4
VlnPlot(sce,features = c('EPCAM','KRT8','KRT18','KRT19'),pt.size = 0,group.by = 'seurat_clusters')
#Fibroblasts 7
VlnPlot(sce,features = c('ACTA2','COL1A2','COL3A1'),pt.size = 0,group.by = 'seurat_clusters')
#Endothelial cell	 8
VlnPlot(sce,features = c('PLVAP','AQP1','ACKR1','VWF'),pt.size = 0,group.by = 'seurat_clusters')

# NK/T cell  1,3,10
VlnPlot(sce,features = c('CD3D','TRAC','CD2','GZMK','NKG7','GZMA','CST7'),pt.size = 0,group.by = 'seurat_clusters')
#B cell   6
VlnPlot(sce,features = c('CD79A','IGHG1','JCHAIN'),pt.size = 0,group.by = 'seurat_clusters')
#Macrophage 0,5
VlnPlot(sce,features = c('C1QB','AIF1','CD68'),pt.size = 0,group.by = 'seurat_clusters')
#mast 9
VlnPlot(sce,features = c('TPSB2','TPSAB1','CPA3'),pt.size = 0,group.by = 'seurat_clusters')
# Proliferating  5,10
VlnPlot(sce,features = c('TOP2A','MKI67'),pt.size = 0,group.by = 'seurat_clusters')

genes = c('ANGPTL3','LBP','FGA',
          'EPCAM','KRT8','KRT18','KRT19',
          'ACTA2','COL1A2','COL3A1',
          'PLVAP','AQP1','ACKR1','VWF',
          'CD3D','TRAC','CD2'	,'IL7R','GZMK','NKG7','GZMA','CST7',
          'CD79A','IGHG1','JCHAIN','C1QB','AIF1','CD68',
          'TPSB2','TPSAB1','CPA3','TOP2A','MKI67'
)
DotPlot(sce, features=genes,group.by = 'seurat_clusters')+coord_flip()+
  scale_color_gradientn(colors=brewer.pal(5,"YlOrRd"))


marker <- data.frame(cluster = 0:10,cell = 0:10)
marker[marker$cluster %in% c(0),2] <- 'Macrophage'
marker[marker$cluster %in% c(1),2] <- 'Natural killer cell'
marker[marker$cluster %in% c(2),2] <- 'Malignant hepatocytes'
marker[marker$cluster %in% c(3),2] <- 'T cell'
marker[marker$cluster %in% c(4),2] <- 'Epithelial cell'
marker[marker$cluster %in% c(5),2] <- ' Proliferating macrophage'
marker[marker$cluster %in% c(6),2] <- 'B cells'
marker[marker$cluster %in% c(7),2] <- 'Fibroblast'
marker[marker$cluster %in% c(8),2] <- 'Endothelial cells'
marker[marker$cluster %in% c(9),2] <- 'Mast cell'
marker[marker$cluster %in% c(10),2] <- 'Proliferating T cell'
marker

sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols=brewer.pal(12,"Paired")
cell_type_tsne=DimPlot(sce,group.by='cell_type',reduction="tsne",label = F,pt.size = .5,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(size = 12),
        legend.position = 'right')
cell_type_tsne
ggsave('results/07.scRNA/cell_type_tsne.pdf',cell_type_tsne,height = 5,width = 5.5)

saveRDS(sce,file = 'results/07.scRNA/sce.rds')


marker.genes2=c('CD3D','TRAC','CD2','CST7','IL7R','GZMK','NKG7','GZMA',
                'PLVAP','AQP1','ACKR1','VWF',
                'TPSB2','TPSAB1','CPA3',
                'CD79A','IGHG1','JCHAIN',
                'ACTA2','COL1A2','COL3A1',
                'C1QB','AIF1','CD68',
                'ANGPTL3','LBP','FGA',
                'KRT8','KRT18','KRT19','EPCAM',
                'TOP2A','MKI67'
)
dotplot_gene_marker=DotPlot(sce, features=marker.genes2,group.by = 'cell_type')
dotplot_data <- dotplot_gene_marker$data
head(dotplot_data)
heatmap_data <- dotplot_data %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled)
head(heatmap_data)

heatmap_data <- column_to_rownames(heatmap_data, var ="id")


heatmap_data_sorted <- heatmap_data#[order(rowMeans(heatmap_data)), ]
# 4. 
marker.heatmap=pheatmap(t(heatmap_data),name = 'Average Expr',
                        border_color = NA,
                        breaks = c(-3, 0, 3),legend_breaks =  seq(-3,3,1),
                        # cellwidth = 16,
                        # cellheight = 15,
                        cluster_rows = F,
                        cluster_cols = T,
                        color = rev(brewer.pal(11,"RdYlBu")),fontsize = 12)
library(ggplotify)
marker.heatmap = as.ggplot(marker.heatmap)
marker.heatmap


# Hepatocytes  2
VlnPlot(sce,features = c('ANGPTL3','LBP','FGA'),pt.size = 0,group.by = 'cell_type')+NoLegend()+xlab('')
#Epithelial cells 2,4
VlnPlot(sce,features = c('EPCAM','KRT8','KRT18','KRT19'),pt.size = 0,group.by = 'cell_type')
#Fibroblasts 7
VlnPlot(sce,features = c('ACTA2','COL1A2','COL3A1'),pt.size = 0,group.by = 'cell_type')
#Endothelial cell	 8
VlnPlot(sce,features = c('PLVAP','AQP1','ACKR1','VWF'),pt.size = 0,group.by = 'cell_type')

# NK/T cell  1,3,10
VlnPlot(sce,features = c('CD3D','TRAC','CD2','GZMK','NKG7','GZMA','CST7'),pt.size = 0,group.by = 'cell_type')
#B cell   6
VlnPlot(sce,features = c('CD79A','IGHG1','JCHAIN'),pt.size = 0,group.by = 'cell_type')
#Macrophage 0,5
VlnPlot(sce,features = c('C1QB','AIF1','CD68'),pt.size = 0,group.by = 'cell_type')
#mast 9
VlnPlot(sce,features = c('TPSB2','TPSAB1','CPA3'),pt.size = 0,group.by = 'cell_type')
# Proliferating  5,10
VlnPlot(sce,features = c('TOP2A','MKI67'),pt.size = 0,group.by = 'cell_type')



cell.prop2=table(sce$cell_type,sce$Samples)
head(cell.prop2)
cell.prop2=as.data.frame(cell.prop2)
df_wide <- pivot_wider(cell.prop2, names_from = Var2, values_from = Freq, values_fill = 0)
head(df_wide)
df_wide=as.data.frame(df_wide)
df_wide$total=as.numeric(apply(df_wide[,-1],1,sum))
df_wide$Percentage <- df_wide$total / sum(df_wide$total) * 100
head(df_wide)
colnames(df_wide)[1]='Cell_Type'

dat22 <- df_wide %>%
  mutate(perc = total/sum(total),
         y = cumsum(total) - 0.5*total,
         label = paste0(round(perc*100,2),"%"))
library(ggrepel)
cell_prop_fig1=ggplot() +
  geom_bar( data = dat22,
            aes(x = 2, y = total, fill =fct_reorder(Cell_Type , y, .desc = T)),
            stat = "identity", width = 1, color = "white")+ coord_polar("y") + 
  geom_text_repel(data = dat22,
                  aes(x = 2.45, y = as.numeric(y), label = label), 
                  size = 3.5,direction = "both", 
                  point.padding = 0.5, box.padding = 1,
                  nudge_x = 0.8,    
                  nudge_y = 0.5,     
                  segment.size = 0.5, 
                  segment.color = "grey50",
                  max.overlaps = 50, 
                  min.segment.length = 0.2, 
                  xlim = c(2.2, 3.5))+ 
  theme_void() +  scale_fill_manual(values =rev(my.cols[1:11]))+labs(fill='Cell type')+
  theme(legend.text = element_text(size=12),legend.title  = element_text(size=14),
        legend.position = 'right')
cell_prop_fig1




library(paletteer)

color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))


p1 <- VlnPlot(sce,features=hub.genes,#c('ESR1','MMP9',"FOS",'CDK1','CDKN2A'),
              group.by="cell_type",
              stack=T,cols=color,flip=T
)


hub.vlnplot=p1+NoLegend()+xlab('')
hub.vlnplot

pdf('results/07.scRNA/Figure7.pdf',height = 12,width = 10)
mg_merge_plot(mg_merge_plot(seurat_clusters_tsne,cell_type_tsne,
                            labels = c('A','B'),widths = c(.6,1,2)),
              mg_merge_plot(marker.heatmap,
                            mg_merge_plot(cell_prop_fig1,hub.vlnplot,
                                          nrow=2,labels = c('D','E')),
                            labels = c('C',''),widths = c(1.2,1)),
              nrow=2,heights = c(1,2))
dev.off()
