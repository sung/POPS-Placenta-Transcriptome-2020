source("lib/local.R") # boiler plates

my.RData<-file.path('RData',paste0(my.gtex.ver,'.li.enriched.RData'))
if(!file.exists(my.RData)){
    message("see placenta_enriched_transcript.R")
}else{
    message(paste("loading ",my.RData))
    load(my.RData)
}

li.enriched[["Placenta.T3.rZ"]][,.N,gene_biotype]

dt.pt.specific<-list()
##########################################
# Heat-map of pt-specific protein-coding #
##########################################
dt.pt.specific[['pc']]<-li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding",-"Placenta.T3.pA"][order(-Placenta.T3.rZ)]
setnames(dt.pt.specific[['pc']],"Placenta.T3.rZ","Placenta")
dt.pt.specific[['pc']] # n=71

mat.top.pt.a<-as.matrix(log10(dt.pt.specific[['pc']][Placenta>100,5:54]+0.001)) # n=43
rownames(mat.top.pt.a)<-dt.pt.specific[['pc']][Placenta>100]$hgnc_symbol
pm.a<-pheatmap::pheatmap(mat.top.pt.a, cluster_rows=F, cluster_cols=T, fontsize_col=13, cutree_cols=2)

#########################################################
# Heat-map of pt-specific long non-coding (FPKM>2) n=48 #
#########################################################
dt.pt.specific[['lincRNA']]<-li.enriched[["Placenta.T3.rZ"]][gene_biotype=="lincRNA",-"Placenta.T3.pA"][order(-Placenta.T3.rZ)]
setnames(dt.pt.specific[['lincRNA']],"Placenta.T3.rZ","Placenta")
dt.pt.specific[['lincRNA']]

mat.top.pt.b<-as.matrix(log10(dt.pt.specific[['lincRNA']][Placenta>10,5:54]+0.001)) # n=26
rownames(mat.top.pt.b)<-dt.pt.specific[['lincRNA']][Placenta>10]$hgnc_symbol
pm.b<-pheatmap::pheatmap(mat.top.pt.b, cluster_rows=F, cluster_cols=T, fontsize_col=13, cutree_cols=2)

#################################
## ERV: endogenous retro virus ##
## retrotransposon             ##
#################################
dt.foo<-fread('data/Tissue-specific/ERV.tpm.csv')[,-"Placenta.T3.pA"] # n=14 ERV
dt.foo[,.N,gene_biotype]
setnames(dt.foo,"Placenta.T3.rZ","Placenta")
dt.target.tpm<-rbind(
                        dt.foo[gene_biotype=="protein_coding"],
                        dt.foo[gene_biotype=="lincRNA"],
                        dt.foo[gene_biotype=="processed_pseudogene"],
                        dt.foo[gene_biotype=="processed_transcript"],
                        dt.foo[gene_biotype=="transcribed_processed_pseudogene"]
                )
dt.target.tpm[,`Placenta-enriched?`:=ifelse(Placenta>1 & Tau.rZ>.99,'Yes',"No")]
dt.target.tpm[,`Placenta-enriched?(relaxed)`:=ifelse(Enriched=="Enriched",'Yes','No')]
dt.target.tpm[,.(ensembl_gene_id,hgnc_symbol,gene_biotype,Placenta,Tau.rZ,Enriched,`Placenta-enriched?`,`Placenta-enriched?(relaxed)`)]
dt.target.tpm[,c("ensembl_gene_id","hgnc_symbol","gene_biotype","Placenta","Tau.rZ","Enriched")]
colnames(dt.target.tpm)

mat.target<-as.matrix(log10(dt.target.tpm[,6:56]+0.001))
rownames(mat.target)<-dt.target.tpm$hgnc_symbol

annotation_row=data.frame(dt.target.tpm[,.(`Gene type`=gsub("_"," ",gene_biotype),`Placenta-enriched?`,`Placenta-enriched?(relaxed)`)])
colnames(annotation_row)=c("Gene type","Placenta-enriched?","Placenta-enriched?(relaxed)")
rownames(annotation_row)=dt.target.tpm$hgnc_symbol

my.Col=c(cbPalette2[1:3],cbPalette2[6:7])
names(my.Col)=gsub("_"," ",dt.target.tpm[,.N,gene_biotype]$gene_biotype)
ann_colors=list(`Gene type`=my.Col,
                    `Placenta-enriched?`=c(`Yes`="grey40",`No`="grey90"),
                    `Placenta-enriched?(relaxed)`=c(`Yes`="grey40",`No`="grey90"))

pm.c<-pheatmap::pheatmap(mat.target, annotation_row=annotation_row, annotation_colors=ann_colors, cluster_rows=F, cluster_cols=T,fontsize_col=8) #, cutree_cols=2, fontsize_col=12)

##############################
## coding + non-coding + ERV #
##############################
my.filename="Figures/SI.Fig3a.pt.specific.rev1.tiff"
tiff(filename=my.filename, width=9, height=12, units="in", res=300, compression='lzw')
cowplot::plot_grid(NULL,pm.a$gtable,nrow=2, labels=c("a",""), rel_heights=c(.2,5),label_size=27)
dev.off()

my.filename="Figures/SI.Fig3b.pt.specific.rev1.tiff"
tiff(filename=my.filename, width=9, height=10, units="in", res=300, compression='lzw')
cowplot::plot_grid(NULL,pm.b$gtable,nrow=2, labels=c("b",""), rel_heights=c(.2,5),label_size=27)
dev.off()

my.filename="Figures/SI.Fig3c.pt.specific.rev1.tiff"
tiff(filename=my.filename, width=9, height=5, units="in", res=300, compression='lzw')
cowplot::plot_grid(NULL,pm.c$gtable,nrow=2, labels=c("c",""), rel_heights=c(.2,5),label_size=27)
dev.off()

