
















#######################  Main script for all RNA-seq related analysis 

cts = read.csv('MMRNHep/data/RSEMCountTable.tab',sep = '\t',header = TRUE)
cts = unique(cts)
row.names(cts) = cts$Gene
cts$Gene = NULL
cts$X = NULL
sample = colnames(cts)
coldata = as.data.frame(sample)
coldata$experiment = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 1 )
coldata$tissue = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 2 )
coldata$diet = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 3 )
coldata$group = paste(coldata$experiment, coldata$tissue, coldata$diet,sep = '_')
row.names(coldata) = coldata$sample

head(cts,n=2)
head(coldata,n=2)

#######################  Load DESeq2 object

dds = DESeqDataSetFromMatrix(countData = cts,colData = coldata,design= ~ group)

#######################  PCA plot - Figure 2C (write to output/Figures)

rld = rlog(dds, blind=FALSE)
pcaData = plotPCA(rld, intgroup=c('group'),returnData=TRUE,ntop=99999999999999)
pcaData$group = factor(pcaData$group,levels=level_order,ordered = TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))
g = ggplot(pcaData, aes(PC1, PC2,color=group) )
g = g + geom_point(size=4)
g = g + ggtitle("Principal component analysis")
g = g + xlab(paste0("PC1: ",percentVar[1],"% variance"))
g = g + ylab(paste0("PC2: ",percentVar[2],"% variance"))
g = g + coord_fixed() 
g = g + theme 
g = g + ggtitle('')
pdf('MMRNHep/output/Figures/Figure2C.pdf')
print(g)
dev.off()

#######################  Perform differential gene expression

comps = c('nonDEN_Liver_WD_vs_nonDEN_Liver_CD',
          'DEN_Liver_CD_vs_nonDEN_Liver_CD',
          'DEN_AdjLiver_WD_vs_nonDEN_Liver_CD',
          'DEN_Tumour_WD_vs_nonDEN_Liver_CD',
          'DEN_AdjLiver_WD_vs_nonDEN_Liver_WD',
          'DEN_Liver_CD_vs_nonDEN_Liver_WD',
          'DEN_Tumour_WD_vs_nonDEN_Liver_WD',
          'DEN_AdjLiver_WD_vs_DEN_Liver_CD',
          'DEN_Tumour_WD_vs_DEN_AdjLiver_WD',
          'DEN_Tumour_WD_vs_DEN_Liver_CD')      # manually specified to ensure correct baseling for each DGE

counter = 0
for (comp in comps){
  counter = counter + 1
  
  #perform DESeq2
  ref = strsplit(comp,'_vs_')[[1]][2]
  dds = DESeqDataSetFromMatrix(countData = cts,colData = coldata,design= ~ group)
  dds$group = relevel(dds$group,ref = ref)
  dds = DESeq(dds)
  res = lfcShrink(dds,coef=paste('group',comp,sep='_'),type='apeglm')
  
  if(counter == 1) {
    genes = row.names(res)
    resTable = as.data.frame(genes)
  }
  
  #add in all the columns
  player = paste('contrast',counter,'lg2BaseMean',comp,sep='_')
  resTable[player] = res$baseMean
  
  player = paste('contrast',counter,'lg10p',comp,sep='_')
  resTable[player] = res$pvalue
  
  player = paste('contrast',counter,'padj',comp,sep='_')
  resTable[player] = res$padj
  
  player = paste('contrast',counter,'logFC',comp,sep='_')
  resTable[player] = res$log2FoldChange
}
# table with DEGs for each contrast
write.csv(resTable,'MMRNHep/output/resLFC.csv')


#######################  get DEGs against nonDEN_liver_CD  - Figure 2D (write to output/Figures)

getDEGs = function(df,contrast){
  df = df[,grepl(contrast, colnames(df))]
  df = subset(df,df[,3] <= 0.05)
  df = subset(df,abs(df[,4]) >= 1)
  return(df)
}

resTable = read.csv('MMRNHep/output/resLFC.csv')
row.names(resTable) = resTable$genes
resTable$genes = NULL
resTable$X = NULL
getDEGs(resTable,'nonDEN_Liver_WD_vs_nonDEN_Liver_CD')
DEGOverlap = list(nonDEN_Liver_WD = rownames(getDEGs(resTable,'nonDEN_Liver_WD_vs_nonDEN_Liver_CD')),
                  DEN_Liver_CD = rownames(getDEGs(resTable,'DEN_Liver_CD_vs_nonDEN_Liver_CD')),
                  DEN_AdjLiver_WD = rownames(getDEGs(resTable,'DEN_AdjLiver_WD_vs_nonDEN_Liver_CD')),
                  DEN_Tumour_WD = rownames(getDEGs(resTable,'DEN_Tumour_WD_vs_nonDEN_Liver_CD')))
colors = scales::viridis_pal()(5)
colors = colors[1:4]
colors = rev(colors)
g = ggvenn(DEGOverlap,stroke_size = 0.5,set_name_size = 3,fill_color = colors)
pdf('MMRNHep/output/Figures/Figure2D.pdf')
print(g)
dev.off()










