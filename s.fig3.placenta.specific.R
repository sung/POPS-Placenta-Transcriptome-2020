source("lib/local.R") # boiler plates
load("RData/dl.gtex.deseq.RData")

    dt.gtex.deseq<-rbindlist(dl.gtex.deseq)
    dt.gtex.deseq[,.N,Tissue] # 20 GTEx + 1 Placenta
    my.good.ensg<-dt.gtex.deseq[,.N,ensembl_gene_id][N==nrow(dt.gtex.deseq[,.N,Tissue])]$ensembl_gene_id 
    length(my.good.ensg) # n=56192

    dt.gtex.tpm<-dt.gtex.deseq[ensembl_gene_id %in% my.good.ensg,.(ensembl_gene_id,hgnc_symbol,chromosome_name,gene_biotype,baseMean,meanFpkm,Tissue)][,TPM:=meanFpkm/sum(meanFpkm)*10^6,Tissue]
    dt.gtex.tpm[,rank:= frank(-TPM,ties.method='dense'),ensembl_gene_id]
    dt.gtex.tpm[,.(sumCnt=sum(baseMean),sumTPM=sum(TPM),sumFPKM=sum(meanFpkm),.N),Tissue]

## Lists of tissue-specific transcripts
    dt.tissue.specific<-rbind(
        fread("Data/Tissue-specific/Tissue.specific.protein_coding.genes.list.TPM.csv"),
        fread("Data/Tissue-specific/Tissue.specific.lincRNA.genes.list.TPM.csv")
        )

#################################################
# Heat-map of pt-specific protein-coding (n=55)
#################################################
    dt.pt.specific.pc<-dt.tissue.specific[Which_Tissue=="Placenta" & gene_biotype.y=="protein_coding"][order(-Placenta)]

    dt.pt.specific.pc[!grepl("HIST",hgnc_symbol.y)] # n=76
    dt.top.pt.a<-dt.pt.specific.pc[!grepl("HIST",hgnc_symbol.y) & Placenta>40,3:23] # n=47
    mat.top.pt.a<-log10(as.matrix(dt.top.pt.a)+0.001)
    colnames(mat.top.pt.a)<-gsub("_"," ",colnames(mat.top.pt.a))
    rownames(mat.top.pt.a)<-dt.pt.specific.pc[!grepl("HIST",hgnc_symbol.y) & Placenta>40,hgnc_symbol.y]

	pm.a<-pheatmap::pheatmap(mat.top.pt.a, cluster_rows=F, cluster_cols=T, fontsize_col=15, cutree_cols=2)

#########################################################
# Heat-map of pt-specific long non-coding (FPKM>2) n=48 #
#########################################################
    dt.pt.specific.lincRNA<-dt.tissue.specific[Which_Tissue=="Placenta" & gene_biotype.y=="lincRNA"][order(-Placenta)]
    dt.pt.specific.lincRNA[!hgnc_symbol.y %in% c("RMRP","AL356488.2")] # n=78
    dt.top.pt.b<-dt.pt.specific.lincRNA[!hgnc_symbol.y %in% c("RMRP","AL356488.2") & Placenta>4,3:23] # n=36
    mat.top.pt.b<-log10(as.matrix(dt.top.pt.b)+0.001)
    colnames(mat.top.pt.b)<-gsub("_"," ",colnames(mat.top.pt.b))
    rownames(mat.top.pt.b)<-dt.pt.specific.lincRNA[!hgnc_symbol.y %in% c("RMRP","AL356488.2") & Placenta>4,hgnc_symbol.y]

	pm.b<-pheatmap::pheatmap(mat.top.pt.b, cluster_rows=F, cluster_cols=T, fontsize_col=15, cutree_cols=2)

#################################
## ERV: endogenous retro virus ##
## retrotransposon             ##
#################################
    dt.target<-fread("Data/ERV.genes.csv") # n=23 ERV genes

    dt.foo<-merge(dt.gtex.tpm[ensembl_gene_id %in% dt.target$ensembl_gene_id,-c("hgnc_symbol","gene_biotype")], dt.target[,.(ensembl_gene_id,description)],all.x=T)
    dt.tpm.foo<-dcast.data.table(merge(dt.foo,dt.ensg[,.(ensembl_gene_id,hgnc_symbol,gene_biotype)]), ensembl_gene_id+hgnc_symbol+gene_biotype+description~Tissue,value.var='TPM')[order(gene_biotype,-Placenta)] # n=17

    dt.target.tpm<-rbind(
                         dt.tpm.foo[gene_biotype=="protein_coding"],
                         dt.tpm.foo[gene_biotype=="lincRNA"],
                         dt.tpm.foo[gene_biotype=="processed_pseudogene"],
                         dt.tpm.foo[gene_biotype=="processed_transcript"],
                         dt.tpm.foo[gene_biotype=="transcribed_processed_pseudogene"]
                    )
    dt.target.tpm[,c("ensembl_gene_id","hgnc_symbol","gene_biotype")]

    m.tpm<-as.matrix(dt.target.tpm[,-c("ensembl_gene_id","hgnc_symbol","gene_biotype","description")])
    mat.target<-log10(m.tpm+0.001)
    colnames(mat.target)<-gsub("_"," ",colnames(mat.target))
    rownames(mat.target)<-dt.target.tpm$hgnc_symbol

    dt.target.tpm[,`Placenta-specific?`:=ifelse(ensembl_gene_id %in% dt.foo[Placenta>1 & FC>100 & Tau>0.99]$ensembl_gene_id,'Yes',"No")]
    dt.target.tpm[,`Placenta-specific?(relaxed)`:=ifelse(ensembl_gene_id %in% dt.foo[Placenta>1 & FC>10 & Tau>0.9]$ensembl_gene_id,'Yes','No')]

    annotation_row=data.frame(dt.target.tpm[,.(`Gene type`=gsub("_"," ",gene_biotype),`Placenta-specific?`,`Placenta-specific?(relaxed)`)])
    colnames(annotation_row)=c("Gene type","Placenta-specific?","Placenta-specific?(relaxed)")
    rownames(annotation_row)=dt.target.tpm$hgnc_symbol

    my.Col=c(cbPalette2[1:3],cbPalette2[6:7])
    names(my.Col)=gsub("_"," ",dt.target.tpm[,.N,gene_biotype]$gene_biotype)
    ann_colors=list(`Gene type`=my.Col,
                    `Placenta-specific?`=c(`Yes`="grey40",`No`="grey100"),
                    `Placenta-specific?(relaxed)`=c(`Yes`="grey40",`No`="grey100"))

	pm.c<-pheatmap::pheatmap(mat.target, annotation_row=annotation_row, annotation_colors=ann_colors, cluster_rows=F, cluster_cols=T, cutree_cols=2, fontsize_col=15)

##############################
## coding + non-coding + ERV #
##############################
    p.ab<-cowplot::plot_grid(NULL,pm.a$gtable,NULL,pm.b$gtable, nrow=1, labels=c("a","","b",""), rel_widths=c(0.2,5,0.2,5),label_size=27, align="h")
    p.cd<-cowplot::plot_grid(pm.c$gtable,NULL, nrow=1, labels=c("c",""), label_size=27,rel_widths=c(5,1),align="h")
	my.filename="Figures/SI.Fig3.pt.specific.tiff"
    tiff(filename=my.filename, width=12, height=14, units="in", res=300, compression='lzw')
    cowplot::plot_grid(p.ab, p.cd, nrow=2, labels=c("","c"), rel_heights=c(2,1.2),label_size=27)
    dev.off()
