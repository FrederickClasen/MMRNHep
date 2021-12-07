
















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

#######################  get enriched GO Biological Processes  - Figure 2E (write to output/Figures)

combinedDotPlot = function(resLFC,homology,contrasts){
  gseDat = list()
  for (contrast in contrasts){
    search = paste('contrast_',contrast,'_',sep='')
    resLFC$entrez = homologyMap$NCBI.gene..formerly.Entrezgene..ID[match(resLFC$genes,homologyMap$Gene.stable.ID)]
    df = resLFC[ , grepl( search , names( resLFC ) ) ]
    parts = str_split(colnames(df)[1],'lg2BaseMean_')
    title = parts[[1]][2]
    print(title)
    df = data.frame(entrez = resLFC$entrez,lfc = df[,4],padj=df[,3])
    df = df[!is.na(df$entrez),]
    df = subset(df,abs(df$lfc) > 1 & df$padj <= 0.05)
    de = unique(df$entrez)
    gseDat[[title]] = de
  }
  ck = compareCluster(geneCluster = gseDat, 
                      fun = "enrichGO",
                      pAdjustMethod = "BH",
                      OrgDb = "org.Mm.eg.db", 
                      ont = "bp", 
                      readable = TRUE, 
                      qvalueCutoff  = 0.05)
  return(ck)
}

resLFC = read.csv('MMRNHep/output/resLFC.csv',header = TRUE,sep = ',')
resLFC$X = NULL
homologyMap = read.csv('MMRNHep/data/homologyMap.tab',header = TRUE, sep = '\t')
contrasts = c(1,2,3,4,5,6,7,8,9,10)
p = combinedDotPlot(resLFC,homology,contrasts)




rxnKO = read.csv('output/reactionKO.csv',sep=',',header = T)
rowAnnot = rxnKO[,1:2]
row.names(rowAnnot) = rowAnnot$X
rowAnnot$X = NULL
row.names(rxnKO) = rxnKO$X
rxnKO$X = NULL
rxnKO$biomass_SUBSYSTEM = NULL
pheatmap(rxnKO,
         cluster_cols = F,
         cluster_rows = T,
         annotation_row = rowAnnot,
         fontsize_row = 1,
         fontsize_col = 10)








