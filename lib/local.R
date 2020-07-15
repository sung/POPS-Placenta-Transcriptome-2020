library(scales)
library(data.table)
library(rtracklayer)
library(DESeq2)
library(BiocParallel)
library(edgeR)
library(RColorBrewer) # for brewer.pal
library(ggplot2)
library(RColorBrewer) # for brewer.pal
library(R.utils) # to read .gz file via fread
library(ggrepel)
library(ggthemes)
library(ggsci)

source("lib/theme_publish.R") # defines theme_Publication()
source("lib/graphic.R")

time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
my.gtex.ver="v8.p2" # v8.p2 (gencode v26; Ensembl v88)

##########################################
## Meta info of RNA-Seq GTEx & placenta ##
##########################################
if(TRUE){
    message("loading meta info...")
    my.pt3.pa.meta<-"data/Meta/placenta.T3.pA.meta.csv"
    my.pt3.rz.meta<-"data/Meta/placenta.T3.rZ.meta.csv"

    ## 1. read GTEx sample meta-info 
    if(!file.exists("data/Meta/GTEx.v8.p2.meta.ALL.csv")){
        makeGTExMeta(my.gtex.ver)
    }
    dt.gtex.meta<-fread(file.path("data/Meta",paste0("GTEx.",my.gtex.ver,".meta.ALL.csv")))
    dim(dt.gtex.meta) # n=13,070 (v8.p2)
    dt.gtex.meta[,.N,SAMPID] # n=13,070 (v8.p2)

    ## 2. read Placenta T3 polyA+ #
    dt.pt3.pa.meta<-fread(my.pt3.pa.meta)[!SampleName %in%c("88")][,SampleName:=paste0('pA',SampleName)]
    nrow(dt.pt3.pa.meta) # n=59
    dt.pt3.pa.meta[,.N,Condition]

    ## 3. Placenta T3 riboZero #
    dt.pt3.rz.meta<-fread(my.pt3.rz.meta)[!SampleName %in% c("80C","84C","88")] 
    nrow(dt.pt3.rz.meta) # n=321
    dt.pt3.rz.meta[Condition==0,.N,CRN] # n=152
    dt.pt3.rz.meta[Condition==0,.N,CRN][N>1]  # n=19

    ## merge the meta info
    dt.pt.meta<-rbind(
        dt.pt3.pa.meta[Condition==0,.(SampleName,CntFile,Tissue="Placenta.T3.pA",Sex,Source="POPS")],
        dt.pt3.rz.meta[Condition==0,.(SampleName,CntFile,Tissue="Placenta.T3.rZ",Sex,Source="POPS")]
    )
    dt.meta<-rbind(
        dt.gtex.meta[,.(SampleName=SAMPID,CntFile="",Tissue=SMTSD,Sex=SEX,Source="GTEx")],
        dt.pt.meta
    )
    dt.meta[,.N,Source]
}


init_meta_gene<-function(){
    my.RData=file.path('RData',"dt.ensg.RData")
    if(file.exists(my.RData)){
        message(paste("loading",my.RData,"..."))
        system.time(attach(my.RData))
    }else{
        my.gtf="data/Homo_sapiens.GRCh38.88.collapsed.from.gencode.v26.gtf"

        gr.ensg<-rtracklayer::import.gff2(my.gtf, format="gtf", feature.type="gene") # isa 'GRanges' 
        mcols(gr.ensg)$gene_id<-substr(mcols(gr.ensg)$gene_id,1,15)

        gr.exon<-rtracklayer::import.gff2(my.gtf, format="gtf", feature.type="exon") # isa 'GRanges'
        mcols(gr.exon)$gene_id<-substr(mcols(gr.exon)$gene_id,1,15)

        gr.target.list<-split(gr.exon, mcols(gr.exon)[["gene_id"]]) #
        system.time(exonic.gene.sizes <- parallel::mclapply(gr.target.list,function(x){sum(width(reduce(x)))},mc.cores=32)) # isa 'list' 
        dt.ensg<-
            data.table(data.frame(chromosome_name=as.character(seqnames(gr.ensg)),
                    strand=strand(gr.ensg),
                    ensembl_gene_id=as.character(mcols(gr.ensg)$gene_id),
                    hgnc_symbol=as.character(mcols(gr.ensg)$gene_name),
                    gene_biotype=as.character(mcols(gr.ensg)$gene_type),
                    length=unlist(exonic.gene.sizes)[mcols(gr.ensg)$gene_id],
                    tag=as.character(mcols(gr.ensg)$tag),
                    stringsAsFactors=F
                    ))
        dt.ensg[tag!="PAR" | is.na(tag)]

        # TODO:  <23-09-20, yourname> # PAR e.g. ENSG00000228572 may have X and Y loci assgined with <PAR> tag on Y
        # RPKM would not be correct
        df.ensg<-data.frame(rbind(
                    dt.ensg[tag=="PAR",.(chromosome_name="X",ensembl_gene_id,hgnc_symbol,gene_biotype,`length`=length/2)],
                    dt.ensg[!ensembl_gene_id %in% dt.ensg[tag=="PAR"]$ensembl_gene_id,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,length)]
                    ))
        rownames(df.ensg)<-df.ensg$ensembl_gene_id
        save(gr.ensg,gr.target.list,dt.ensg,df.ensg,file=my.RData)
    }
}

# Suppl information from: https://academic.oup.com/bib/article/18/2/205/2562739
# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues
fTau <- function(x){
    if(all(!is.na(x))){
        if(min(x, na.rm=TRUE) >= 0){
            if(max(x)!=0){
                x <- (1-(x/max(x)))
                res <- sum(x, na.rm=TRUE)
                res <- res/(length(x)-1)
            }else{
                res <- 0
            }
 		}else{
            res <- NA
            #print("Expression values have to be positive!")
 		} 
 	}else{
        res <- NA
        #print("No data for this gene avalable.")
    } 
    return(res)
}

############################
## Parsing GTEx Meta Data ##
## 'v6p' or 'v8.p2'       ##
############################
makeGTExMeta<-function(my.gtex.ver="v8.p2"){
	gtex.dir<-file.path("data", paste("GTEx",my.gtex.ver,sep='.'))
    stopifnot(file.exists(gtex.dir))

    list.files(gtex.dir,'Annotations_SubjectPhenotypesDS.txt')
    list.files(gtex.dir,'Annotations_SampleAttributesDS.txt')
    list.files(gtex.dir,'gene_reads')

    # gene-count file
    dt.gtex.cnt<-fread(file.path(gtex.dir, list.files(gtex.dir,'gene_reads')[1]), na.strings="")
    dt.gtex.cnt[,Name:=substr(Name,1,15)] # ENSG00000223972.4 => ENSG00000223972
    sample.in.cnt<-colnames(dt.gtex.cnt)[-c(1,2)]
    length(sample.in.cnt) # n=17,382 samples
    dim(dt.gtex.cnt) #  56,200 x 17,384 (v8.p2)
    dim(dt.gtex.cnt[,.N,Name])
    dt.gtex.cnt[,.N,Name][N>1]
    dt.gtex.cnt[1:5,1:5]
    dt.gtex.cnt<-dt.gtex.cnt[,lapply(.SD,sum),.(Name,Description)] # aggregate cnt for the same gene 
    dim(dt.gtex.cnt) #  56,156 x 17,384 (v8.p2)
    dim(dt.gtex.cnt[,.N,Name])

    # sample meta-info
	dt.gtex.subject<-fread(file.path(gtex.dir,list.files(gtex.dir,'Annotations_SubjectPhenotypesDS.txt')[1]), na.strings="", select=1:3) # patient (subject) meta-info 
    # run meta-info
	dt.gtex.meta<-fread(file.path(gtex.dir,list.files(gtex.dir,'Annotations_SampleAttributesDS.txt')[1]), na.strings="") # rna-seq run meta-info
    dim(dt.gtex.subject) # 980 x 3 (v8.p2)
    dim(dt.gtex.meta) #  22,951 x 63 (v8.p2)

    dt.gtex.meta[,.N,SAMPID %in% sample.in.cnt] # samples where cnt not available
    table(sample.in.cnt %in% dt.gtex.meta$SAMPID) # samples where run meta-info not available

	# http://stackoverflow.com/questions/18154556/split-text-string-in-a-data-table-columns
	# A-B-C-D => A-B
	# data.table > v1.9.6+
	dt.gtex.meta[,c("foo","bar"):= tstrsplit(SAMPID, "-", fixed=TRUE, keep=c(1,2))]
	dt.gtex.meta[,SUBJID:=paste(foo,bar,sep="-")]
	dt.gtex.meta[,c("foo","bar"):=list(NULL,NULL)]
	dt.gtex.meta<-merge(dt.gtex.meta, dt.gtex.subject) # add GENDER|SEX

    if("GENDER"%in%colnames(dt.gtex.meta)){
	    dt.gtex.meta[,SEX:=ifelse(GENDER==1,"M","F")]
    	dt.gtex.meta[,GENDER:=NULL]
    }else{
	    dt.gtex.meta[,SEX:=ifelse(SEX==1,"M","F")]
    }

	# SMRIN: RIN
	# SMMAPRT: mapping ratio
	# SMEXNCRT: mapping ratio to exon
	# SMMNCPB: coverage per base (only v6p; NA for v8.p2)
    dt.gtex.meta[SAMPID %in% sample.in.cnt,.N,SMTS][order(N)] # no-filtering: n=30 tissues (v8.p2) 
    dt.gtex.meta[SAMPID %in% sample.in.cnt,.N,SMTSD][order(N)] # no-filtering: n=54 sub-tissues (v8.p2) 
    dt.gtex.meta[SAMPID %in% sample.in.cnt,.N,SMTSD][order(SMTSD)] # no-filtering: n=54 sub-tissues (v8.p2) 
    dt.gtex.meta[SAMPID %in% sample.in.cnt,quantile(SMRIN)] 
    dt.gtex.meta[SAMPID %in% sample.in.cnt & SMTS=="Testis",.N,.(SMEXNCRT<0.7)]
    dt.gtex.meta[SAMPID %in% sample.in.cnt & SMTSD=="Kidney - Medulla"]

    if(my.gtex.ver=="v6p"){
	    dt.gtex.meta[!is.na(SMRIN) & SMMAPRT>0.9 & SMEXNCRT>0.8 & SMMNCPB>20]
	    dt.gtex.sample<-dt.gtex.meta[!is.na(SMRIN) & SMMAPRT>0.9 & SMEXNCRT>0.8 & SMMNCPB>20, .(SAMPID,SUBJID,SMTS,SMTSD,SMRIN,SEX,AGE)] # n=6611 (v6p); 
    }else if(my.gtex.ver=="v8.p2"){
        # NB, SMMNCPB not avilable in v8.p2
	    dt.gtex.sample<-dt.gtex.meta[SMRIN>=6 & SMMAPRT>=0.9 & SMEXNCRT>=0.75 & SAMPID %in% sample.in.cnt, .(SAMPID,SUBJID,SMTS,SMTSD,SMRIN,SMMAPRT,SMEXNCRT,SEX,AGE)]
    }
	#dt.gtex.sample[,SMTS:=gsub(" ","_",SMTS)]
    dt.gtex.sample[,.N,SMTS][order(N)] # n=30 tissues (v8.p2) 
    dt.gtex.sample[,.N,"SMTS,SMTSD"][order(N)] # n=53 sub-tissues (v8.p2) 

    ###########
    ## Filter #
    ###########
	# 1. tissue of less than 20 sample per sub-tissue 
	dt.gtex.sample[,.N,"SMTS,SMTSD"][order(SMTS,SMTSD)][N<20]
	dt.gtex.sample[,.N,"SMTS,SMTSD"][order(SMTS,SMTSD)][N>=20][order(N)]#[,sum(N)] # 49 sub-tissues
	less.than.x.SMTSD<-dt.gtex.sample[,.N,"SMTS,SMTSD"][order(SMTS,SMTSD)][N<20,unique(SMTSD)] # n=4 sub-tissues (v8.p2)
    #  "Bladder"             "Cervix - Ectocervix" "Cervix - Endocervix" "Fallopian Tube"

    # No. of samples before-after filtering
    dt.foo<-
        merge(dt.gtex.meta[SAMPID %in% sample.in.cnt,.N,"SMTS,SMTSD"],
        dt.gtex.sample[,.N,SMTSD]
        ,by="SMTSD",all.x=T
    )
    fwrite(dt.foo,file="data/Tissue-specific/sample.number.before.after.filtering.txt")

    ###################
	# Write meta file #
    ###################
	dt.gtex.final<-
        dt.gtex.sample[!SMTSD %in% less.than.x.SMTSD] # n=13070
	dt.gtex.final[,.N,"SMTS,SMTSD"][order(SMTS,SMTSD)]
	dt.gtex.final
	write.csv(dt.gtex.final,file=paste0("data/Meta/GTEx.",my.gtex.ver,".meta.ALL.csv"),row.names=F, quote=F)

	##############################
	# per sample gene-count file #
	##############################
    if(FALSE){
        parallel::mclapply(dt.gtex.final$SAMPID, function(i){
                            write.table(dt.gtex.cnt[,.(Name,get(i))],
                                        sep="\t", quote=FALSE, 
                                        col.names=FALSE, 
                                        row.names=FALSE, 
                                        file=gzfile(paste0("data/GTEx.",my.gtex.ver,"/CNT/",i,".txt.gz")))}
        ,mc.cores=32)
        }
} # end of makeMeta


##############
# run DESeq2 #
##############
runDESeq<-function(mat.cnt, df.meta){
    register(MulticoreParam(32))
    message("DESeqDataSetFromMatrix...")
    deseq.design=formula(~Tissue)
    dds <- DESeqDataSetFromMatrix(mat.cnt, df.meta, deseq.design)
    #dds <- DESeqDataSetFromMatrix(mat.cnt, df.meta, formula(~1))
    message("estimateSizeFactors...")
    dds <- estimateSizeFactors(dds)
    rowRanges(dds) <- gr.target.list[rownames(dds)] # sorted by the original rownames
    return(dds)
}

##############################
# DESeq2 across all tissues ##
##############################
init_dds<-function(){
    my.RData<-file.path("RData",paste0(my.gtex.ver,".dds.RData"))
    if(file.exists(my.RData)){
        message(paste("loading",my.RData,"..."))
        system.time(attach(my.RData))
    }else{
        ###########################
        ## merge gene count file ##
        ###########################
        # 1. GTEx
        gtex.dir<-file.path("data", paste("GTEx",my.gtex.ver,sep='.'))
        stopifnot(file.exists(gtex.dir))
        dt.gtex.cnt<-fread(file.path(gtex.dir, list.files(gtex.dir,'gene_reads')[1]), na.strings="")
        dt.gtex.cnt[,Name:=substr(Name,1,15)] # ENSG00000223972.4 => ENSG00000223972
        dim(dt.gtex.cnt) #  56,200 x 17,384 (v8.p2)
        dt.gtex.cnt[,.N,Name]
        dt.gtex.cnt<-dt.gtex.cnt[,lapply(.SD,sum),.(Name,Description)] # aggregate cnt for the same gene 
        dim(dt.gtex.cnt) #  56,156 x 17,384 (v8.p2)
        setnames(dt.gtex.cnt,"Name","ensembl_gene_id")
        dt.gtex.cnt[1:5,1:5]

        # 2. placenta
        foo<-apply(dt.pt.meta, 1, function(j){
                            #fread(j[2],col.names=c("ensembl_gene_id","cnt"))[,`:=`(SampleName=j[1])]
                            fread(j[2],col.names=c("ensembl_gene_id","cnt"))[,.(`cnt`=sum(cnt)),ensembl_gene_id][,`:=`(SampleName=j[1])]
                        }
                    ) # apply fread for each row (1); isa 'list'
        dt.pt.cnt<-dcast.data.table(
            rbindlist(foo),
            ensembl_gene_id~SampleName,value.var="cnt")[order(ensembl_gene_id)]
        dim(dt.pt.cnt) # 56,156 x 201
        dt.pt.cnt[1:5,1:5]
        dt.pt.cnt[,.N,ensembl_gene_id]

        this.col<-c("ensembl_gene_id",dt.gtex.meta$SAMPID)
        dim(dt.gtex.cnt[,..this.col])
        dt.cnt.all<-merge(dt.gtex.cnt[,..this.col], dt.pt.cnt,by="ensembl_gene_id")[order(ensembl_gene_id)]
        mat.cnt.all<-as.matrix(dt.cnt.all[,-"ensembl_gene_id"])
        rownames(mat.cnt.all)<-dt.cnt.all$ensembl_gene_id
        stopifnot(ncol(mat.cnt.all)==nrow(dt.meta))
        dim(mat.cnt.all) # 56,156 x 13,270

        ########################################
        # get a list of common genes (ensg id) #
        ########################################
        # filter 1
        filter1<-rowSums(mat.cnt.all)==0
        table(filter1) # TRUE: 169; FALSE: 55,987
        zero.ensg<-rownames(mat.cnt.all[filter1,])
        length(zero.ensg) # n=169 (v8.p2)
        non.zero.ensg<-rownames(mat.cnt.all[!filter1,])
        length(non.zero.ensg) # n=55,987 (v8.p2)
        #stopifnot(length(non.zero.ensg)==nrow(mat.cnt.all)-table(filter1)["TRUE"])

        # filter 2: polyA(-) non-polyadenylated transcript from Yang et al
        filter2<-grepl("^HIST",mcols(gr.ensg)$gene_name) | mcols(gr.ensg)$gene_name %in% c("RMRP","AL356488.2")
        table(filter2) # TRUE: 90; FALSE: 56110
        non.pA.ensg<-mcols(gr.ensg)[filter2,]$gene_id
        length(non.pA.ensg) # n=90
        table(non.zero.ensg %in% non.pA.ensg) # 90 TRUE; 55,897 FALSE

        # filter 3: non protein-coding transcript
        filter3<-mcols(gr.ensg)$gene_biotype!="protein_coding"
        table(filter3)
        non.pc.ensg<-mcols(gr.ensg)[filter3,]$gene_id
        length(non.pc.ensg) # n=38,405 

        table(
        non.zero.ensg[!(non.zero.ensg %in% non.pA.ensg)] %in% non.pc.ensg
        )

        # legacy codes
        #filter.gene<-mcols(gr.ensg)$gene_id %in% non.zero.ensg & !grepl("^HIST",mcols(gr.ensg)$gene_name) & !mcols(gr.ensg)$gene_name %in% c("RMRP","AL356488.2")
        #table(filter.gene) # 259  F; 55,941 T
        #my.common.ensg<-sort(gr.ensg[filter.gene]$gene_id)
        #length(my.common.ensg) # n=55,897
        #my.common.ensg[table(my.common.ensg)>1]
        #gr.ensg[gr.ensg$gene_id=="ENSG00000002586",]

        filter.gene<-rownames(df.ensg) %in% non.zero.ensg & !grepl("^HIST",df.ensg$hgnc_symbol) & !df.ensg$hgnc_symbol %in% c("RMRP","AL356488.2")
        table(filter.gene) # 259  F; 55,897 T
        my.common.ensg<-df.ensg[filter.gene,]$ensembl_gene_id
        length(my.common.ensg) # n=55,897


        filter.gene<-mcols(gr.ensg)$gene_id %in% my.common.ensg & mcols(gr.ensg)$gene_biotype=="protein_coding"
        my.common.ensg.pc<-sort(gr.ensg[filter.gene]$gene_id)
        length(my.common.ensg.pc) # n=19,173

        ################################
        ## protein-coding only (mRNA) ##
        ################################
        if(protein.coding){
            mat.cnt<-mat.cnt.all[my.common.ensg.pc,] # isa 'matrix'
            dim(mat.cnt) # 19,156 x 4,527 | 19,203 x 4,527
            df.meta<-data.frame(dt.meta, row.names=dt.meta$SampleName)[colnames(mat.cnt),] # same order of samples from mat.cnt
            df.meta[c(1:ncol(df.meta))]=lapply(df.meta[c(1:ncol(df.meta))], as.factor) 
            dds.pc<-runDESeq(mat.cnt, df.meta)
            message("saving dds...")
            save(dds.pc, file=my.RData)
            attach(dds.pc)
        }else{
            mat.cnt<-mat.cnt.all[my.common.ensg,] # isa 'matrix'
            dim(mat.cnt) # 55,897 x 13,270

            df.meta<-data.frame(dt.meta, row.names=dt.meta$SampleName)[colnames(mat.cnt),] # same order of samples from mat.cnt
            df.meta[c(1:ncol(df.meta))]=lapply(df.meta[c(1:ncol(df.meta))], as.factor) 

            save(mat.cnt,df.meta,file=file.path("RData",paste0(my.gtex.ver,".mat.cnt.RData")))

            dds<-runDESeq(mat.cnt, df.meta)
            message("saving dds...")
            save(dds, file=my.RData)
            attach(dds)
        }
    }
} # end of init_dds()

##
## FPKM & TPM from dds all-tissue vs. dds a-tissue
##
init_norm_tissue<-function(){
    my.RData<-file.path("RData",paste0(my.gtex.ver,".li.rpkm.tissue.RData")) # 78M 
    if(file.exists(my.RData)){
        message(paste("loading",my.RData,"..."))
        system.time(attach(my.RData))
    }else{
        init_dds() # load 'dds.pc' or 'dds'
        message("setting edgeR object...")

        system.time(d<-DGEList(counts(dds),
                            samples=colData(dds)[,c("SampleName","Tissue","Sex","Source")],
                            group=colData(dds)$Tissue))
        system.time(d<-calcNormFactors(d,method="TMM")) # TMM is the default 

        #li.norm<-list()
        #message("DESeq normalized count..")
        #system.time(li.norm[["DESeq.ncnt"]]<-counts(dds,normalized=T))
        #message("DESeq CPM...")
        #system.time(li.norm[["DESeq"]]<-fpm(dds)) # use the size factors to normalize 
        #                                                # rather than taking the column sums of the raw counts.
        #message("TMM normalized count...")
        #e.lib<-d$samples$lib.size * d$samples$norm.factors # effective lib size
        #system.time(li.norm[["TMM.ncnt"]]<-sweep(d$counts,2,e.lib,"/")*mean(e.lib))
        #
        #message("aggregating normalized values by tissues...")
        #system.time(
        #    li.norm.tissue<-lapply(li.norm,function(x){
        #                        simplify2array(
        #                                        lapply(split(colData(dds), colData(dds)$Tissue), function(j)
        #                                                rowMeans(x[,rownames(j)]))
        #                            )}
        #                )
        #)
        #system.time(li.norm.tissue[["TMM"]]<- cpmByGroup(d,log=F)) # a function of 'edgeR'
        li.rpkm.tissue<-list()
        d$genes<-df.ensg[rownames(d),]
        system.time(li.rpkm.tissue[["TMM"]]<- rpkmByGroup(d,log=F)) # a function of 'edgeR'
        save(li.rpkm.tissue,file=my.RData)
    }
} # end of init_norm_tissue()
