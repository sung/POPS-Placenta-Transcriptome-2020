system.time(source("libs/local.R"))

init_meta_gene() # df.ensg, gr.ensg, dt.ensg

my.RData<-file.path('RData',paste0(my.gtex.ver,'.dt.tpm.tau.RData'))
if(file.exists(my.RData)){
        message(paste("loading",my.RData,"..."))
        system.time(load(my.RData))
}else{
    init_norm_tissue() # load 'li.rpkm.tissue'
    # calculate TPM
    mat.tpm<-sweep(li.rpkm.tissue[["TMM"]],2,colSums(li.rpkm.tissue[["TMM"]]),"/")*1e+6
    dim(mat.tpm)
    apply(mat.tpm,2,sum)
    # calculate tau (omitting placenta polyA data and omitting placenta ribo-zero data)
    tau<-list()
    tau[["Placenta.T3.rZ"]]<-apply(mat.tpm[,-which(colnames(mat.tpm)=="Placenta.T3.pA")],1, fTau)
    tau[["Placenta.T3.pA"]]<-apply(mat.tpm[,-which(colnames(mat.tpm)=="Placenta.T3.rZ")],1, fTau)
    # calculate mean TPM of GTEx tissues
    gtex.mean.tpm<-apply(mat.tpm[,!grepl("Placenta",colnames(mat.tpm))],1,mean)
    gtex.median.tpm<-apply(mat.tpm[,!grepl("Placenta",colnames(mat.tpm))],1,median)
    gtex.max.tpm<-apply(mat.tpm[,!grepl("Placenta",colnames(mat.tpm))],1,max)

    dt.tpm.tau<-data.table(
        df.ensg[rownames(mat.tpm),],
        mat.tpm,
        `GTEx_meanTPM`=gtex.mean.tpm,
        `GTEx_medianTPM`=gtex.median.tpm,
        `GTEx_maxTPM`=gtex.max.tpm,
        `Tau.rZ`=tau[["Placenta.T3.rZ"]],
        `Tau.pA`=tau[["Placenta.T3.pA"]]
        )
    save(dt.tpm.tau,file=my.RData)
}
dt.tpm.tau[Placenta.T3.rZ>1 & Tau.rZ>0.99 & Placenta.T3.rZ/GTEx_meanTPM>100,.N,gene_biotype]

####################################
# Find tissue-enriched transcripts #
####################################
my.RData<-file.path('RData',paste0(my.gtex.ver,'.li.enriched.RData'))
if(file.exists(my.RData)){
    load(my.RData) # load li.enriched
else{
    li.enriched<-list()
    for(my.tissue in colnames(mat.tpm)){
        if(my.tissue=="Placenta.T3.pA"){
            my.pt=my.tissue
            not.my.pt="Placenta.T3.rZ"
        }else{
            my.pt="Placenta.T3.rZ"
            not.my.pt="Placenta.T3.pA"
        }
        other.tissues<-!colnames(mat.tpm) %in% c(my.tissue,not.my.pt)
        other.mean.tpm<-apply(mat.tpm[,other.tissues],1,mean)

        # filter 1: TPM(placenta) > 1
        ft1 <- mat.tpm[,my.tissue] > 1
        # filter 2: Tau>0.99
        ft2 <- tau[[my.pt]]>.99
        # filter 3: Placenta> mean(GTEx)*100
        ft3 <- mat.tpm[,my.tissue] > other.mean.tpm * 100
        # filter 4: protein-coding or lincRNA
        ft4<-df.ensg[rownames(mat.tpm),"gene_biotype"] %in% c("protein_coding","lincRNA")

        my.genes<-rownames(mat.tpm)[ft1 & ft2 & ft3 & ft4]

        if(length(my.genes)>1){
            li.enriched[[my.tissue]]<-data.table(
                Which_Tissue=my.tissue,
                cbind(
                    df.ensg[my.genes,c("ensembl_gene_id","hgnc_symbol","gene_biotype")],
                    mat.tpm[my.genes,],
                    `Other_meanTPM`=other.mean.tpm[my.genes],
                    `Tau`=tau[[my.pt]][my.genes]
                    )
            )
        }else if(length(my.genes)==1){
            li.enriched[[my.tissue]]<-data.table(
                Which_Tissue=my.tissue,
                cbind(
                    df.ensg[my.genes,c("ensembl_gene_id","hgnc_symbol","gene_biotype")],
                    t(data.frame(mat.tpm[my.genes,])),
                    `Other_meanTPM`=other.mean.tpm[my.genes],
                    `Tau`=tau[[my.pt]][my.genes]
                    )
            )
        }
    }
    save(li.enriched,file=my.RData)
    # SI table of tissue-enriched transcript #
    fwrite(rbindlist(li.enriched), file.path('data/Tissue-specific',paste0(my.gtex.ver,'.tissue-enriched.csv')))
}

##
## ERV
##
dt.tpm.tau[grepl('ERV',hgnc_symbol) & Placenta.T3.rZ>1 & Tau.rZ>0.9 & Placenta.T3.rZ/GTEx_meanTPM>10,.(hgnc_symbol,gene_biotype,Placenta.T3.rZ/GTEx_meanTPM,Tau.rZ)]
dt.erv<-merge(
    dt.tpm.tau[grepl('ERV',hgnc_symbol)],
    dt.tpm.tau[grepl('ERV',hgnc_symbol) & Placenta.T3.rZ>1 & Tau.rZ>0.9 & Placenta.T3.rZ/GTEx_meanTPM>10,.(ensembl_gene_id,`Enriched`="Enriched")]
    ,all.x=T
    )[order(gene_biotype,-Placenta.T3.rZ)]
fwrite(dt.erv,file=file.path('data/Tissue-specific/ERV.tpm.csv'))

######################################################
## Chi-sq test for the NO. of tissue-enriched genes ##
######################################################
print(foo1<-rbindlist(li.enriched)[Which_Tissue!="Placenta.T3.pA",.N,.(Which_Tissue,gene_biotype)][,.(total=sum(N)),gene_biotype])

print(foo2<-rbindlist(li.enriched)[Which_Tissue!="Placenta.T3.pA",.(`observed`=.N),.(Which_Tissue,gene_biotype)] )

print(bar<-merge(
      merge(foo1, foo2),
    merge(foo1, foo2)[,.N,gene_biotype]
    )[,.(observed,expected=total/N,ratio=observed/(total/N)),.(gene_biotype,Which_Tissue)]
)

# protein-coding
x<-bar[gene_biotype=="protein_coding"]$observed
names(x)<-bar[gene_biotype=="protein_coding"]$Which_Tissue
print(y<-chisq.test(x))

names(y)
data.frame(y$observed,y$expected,y$observed/y$expected)
y$p.value
y$statistic
y$parameter
y$method

options(scipen=0)
formatC(y$p.value, format = "e", digits = 2)

# lincRNA
x<-bar[gene_biotype=="lincRNA"]$observed
names(x)<-bar[gene_biotype=="lincRNA"]$Which_Tissue
print(y<-chisq.test(x,rescale.p=T))

data.frame(y$observed,y$expected,y$observed/y$expected)
y$p.value
y$statistic
y$parameter
y$method

# via lapply
lapply(split(bar, bar$gene_biotype), function(i){
        i[,.(Which_Tissue,observed,expected)]
        x<-i$observed
        names(x)<-i$Which_Tissue
        print(y<-chisq.test(x))
        y$p.value
    } 
)

dcast(rbindlist(li.enriched)[Which_Tissue!="Placenta.T3.pA",.N,.(Which_Tissue,gene_biotype)], Which_Tissue~gene_biotype,fill=0)[order(protein_coding)]
dcast(rbindlist(li.enriched)[Which_Tissue!="Placenta.T3.pA",.N,.(Which_Tissue,gene_biotype)], Which_Tissue~gene_biotype,fill=0)[order(lincRNA)]

##
## Several groups of closely related placenta-enriched transcripts
#
if(TRUE){
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding" & grepl("PSG",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding" & grepl("LGA",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding" & grepl("MAGE",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding" & grepl("SLC",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding" & grepl("ERV",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
    li.enriched[["Placenta.T3.rZ"]][gene_biotype=="lincRNA" & grepl("ERV",hgnc_symbol),.(ensembl_gene_id,hgnc_symbol,Placenta.T3.rZ,Other_meanTPM,Tau)][order(-Placenta.T3.rZ)]
}

##
## Check some genes
##
if(TRUE){
    sort(mat.tpm[dt.ensg[hgnc_symbol=="EXPH5"]$ensembl_gene_id,])
    tau[["Placenta.T3.rZ"]][dt.ensg[hgnc_symbol=="EXPH5"]$ensembl_gene_id]
    sort(mat.tpm[dt.ensg[hgnc_symbol=="EXPH5"]$ensembl_gene_id,])

    sort(mat.tpm[dt.ensg[hgnc_symbol=="MAGEA4"]$ensembl_gene_id,])
    tau[["Placenta.T3.rZ"]][dt.ensg[hgnc_symbol=="MAGEA4"]$ensembl_gene_id]

    foo2<-counts(dds,normalized=T)["ENSG00000147381",]
    sort(sapply(split(colData(dds), colData(dds)$Tissue), function(j) mean(foo2[rownames(j)])))
}

##
## Ribo-zero vs polyA+
##
if(TRUE){
    li.enriched[["Placenta.T3.rZ"]][,.N,gene_biotype]
    li.enriched[["Placenta.T3.pA"]][,.N,gene_biotype]

    merge(
        li.enriched[["Placenta.T3.rZ"]][,.(ensembl_gene_id,hgnc_symbol,gene_biotype,Placenta.T3.rZ,Tau)],
        li.enriched[["Placenta.T3.pA"]][,.(ensembl_gene_id,Placenta.T3.pA,Tau)],by="ensembl_gene_id"
        )[,.N,gene_biotype]
}

##
## current (GTEx v8.p2) vs. previous (GTEx v6p)
##
if(TRUE){
    pt.pc.v1<-fread('data/Tissue-specific/Tissue.specific.protein_coding.genes.list.TPM.csv')[Which_Tissue=="Placenta" & !grepl('HIST',hgnc_symbol.y)]
    pt.linc.v1<-fread('data/Tissue-specific/Tissue.specific.lincRNA.genes.list.TPM.csv')[Which_Tissue=="Placenta"]

    merge(
        li.enriched[["Placenta.T3.rZ"]][gene_biotype=="protein_coding",.(ensembl_gene_id,hgnc_symbol,gene_biotype,Placenta.T3.rZ,Tau)],
        pt.pc.v1[,.(ensembl_gene_id,hgnc_symbol.y,Me_TPM,Tau)]
        ,by="ensembl_gene_id",all=T
        )[!is.na(Tau.x) & is.na(Tau.y)][order(-Me_TPM)]

    merge(
        li.enriched[["Placenta.T3.rZ"]][gene_biotype=="lincRNA",.(ensembl_gene_id,hgnc_symbol,gene_biotype,Placenta.T3.rZ,Tau)],
        pt.linc.v1[,.(ensembl_gene_id,hgnc_symbol.y,Me_TPM,Tau)]
        ,by="ensembl_gene_id",all=T
        )[!is.na(Tau.x) & is.na(Tau.y)][order(-Me_TPM)]
}
