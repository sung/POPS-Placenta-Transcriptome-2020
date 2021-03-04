#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# circRNA prediced by CIRI2
# https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbx014/3058729

TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
ENS_VER=90 # 82 (for reconstruction) 90 (exon/intron/intergenic info only)
source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R")

my.RData<-"~/results/RNA-Seq/Placentome/RData/circRNA.RData"
if(file.exists(my.RData)){
    load(my.RData)
}else{
    #########################
    # FG CIRI2 result files #
    #########################
    my.cmd<-"for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-916*.Homo_sapiens.SE125.v2/CIRI2/D*/SLX-916*.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done"
    fg.meta<-system(my.cmd,intern=T)

    ########################
    # JD CIRI2 result file #
    ########################
    my.cmd<-c(
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-9792.Homo_sapiens.v1/CIRI2/D*/SLX-9792.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10281.Homo_sapiens.v1/CIRI2/D*/SLX-10281.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10283.Homo_sapiens.v1/CIRI2/D*/SLX-10283.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10284.Homo_sapiens.v1/CIRI2/D*/SLX-10284.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10285.Homo_sapiens.v1/CIRI2/D*/SLX-10285.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10287.Homo_sapiens.v1/CIRI2/D*/SLX-10287.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
            "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10402.Homo_sapiens.v1/CIRI2/D*/SLX-10402.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done"
            )
    jd.meta<-unlist(
                sapply(my.cmd, system, intern=T, USE.NAMES=F)
                )

    dt.fg.jd<-data.table(t(simplify2array(strsplit(c(fg.meta,jd.meta), " ")))) # n=324 (114 FG + 210 JD)

    dt.pops<-fread("~/results/RNA-Seq/Placentome/Meta/POPS.stringtie.GRCh38.82.gtf.csv")[,-"StringTieFile"] # see bin/R/Placentome/local.R
                                                                                         # final list of samples for trasncriptome reconstruction (input for StringTie)
                                                                                         # n=295 (152 controls; 91 PET; 52 SGA)
    outLiers=list(`Failed`=c("91"), `IGFBP1`=c("84C","88","80C"), `SGA.PET`=c("77","71","79","19P","02P","72P"))
    dt.pops[SampleName %in% unlist(outLiers)] # should be 0 records

    # 12-pairs excluded (and identified by JD) either batch-effect or decidual contamination (30/APR/2018) 
	outLiers=c("PET08","PET16","PET20","PET60","PET64","PET75","PET76","PET77","PET78","PET79","PET80","PET84") # 12 pairs (n=24)
    dt.pops[Pair %in% outLiers]

    # read all the CIRI2 result files for the samples used in placenta re-construction
    dl.pops.ciri<-apply(merge(dt.fg.jd, dt.pops, by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n= 295
                   1, function(i){fread(i[3])[,-"junction_reads_ID"][,`:=`(SLX=i[1],Barcode=i[2])]})
    dt.pops.ciri<-rbindlist(dl.pops.ciri)
    dt.pops.ciri[,.N,.(SLX,Barcode)] # No of circRNA per sample 
    write.csv(dt.pops.ciri, file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.CIRI2.txt.gz"), row.names=F, quote=F)

    # read all the CIRI2 result files for the samples used in all the placenta samples 
    dl.pops.ciri.all<-apply(dt.fg.jd, # n= 324
                   1, function(i){fread(i[3])[,-"junction_reads_ID"][,`:=`(SLX=i[1],Barcode=i[2])]})
    dt.pops.ciri.all<-rbindlist(dl.pops.ciri.all)
    dt.pops.ciri.all[,.N,.(SLX,Barcode)] # No of circRNA per sample 
    write.csv(dt.pops.ciri.all, file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.CIRI2.all.txt.gz"), row.names=F)

    #############
    # circRNA   #
    # No filter #
    #############
    dt.circRNA<-dt.pops.ciri[,
                             .(`sum.junc.reads`=sum(`#junction_reads`),
                               `sum.non.junc.read`=sum(`#non_junction_reads`),
                               `avg.junc.ratio`=sum(`#junction_reads`)*2/(sum(`#junction_reads`)*2+sum(`#non_junction_reads`)),
                               `freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops),
                               `cnt.sample`=length(unique(paste(SLX,Barcode,sep=".")))
                               ),
                    .(circRNA_ID,chr,circRNA_start,circRNA_end,circRNA_type,gene_id,strand)
                    ][order(-cnt.sample)]
    # remove ',' from gene_id
    #####################
    # circRNA by Cohort #
    #####################
    dt.circRNA.cohort<-merge(dt.pops.ciri, dt.pops,by.x=c("SLX","Barcode"),by.y=c("Library","BarCode"))[,
                            .(`sum.junc.reads`=sum(`#junction_reads`),
                            `sum.non.junc.read`=sum(`#non_junction_reads`),
                            `avg.junc.ratio`=sum(`#junction_reads`)*2/(sum(`#junction_reads`)*2+sum(`#non_junction_reads`)),
                            `freq`=ifelse(Cohort=="PET",
                                          length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops[Cohort=="PET"]),
                                          ifelse(Cohort=="SGA",
                                                 length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops[Cohort=="SGA"]),
                                                length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops[Cohort=="CTL"])
                                            )
                                    ),
                            `cnt.sample`=length(unique(paste(SLX,Barcode,sep=".")))
                            ),
                .(Cohort,circRNA_ID,chr,circRNA_start,circRNA_end,circRNA_type,gene_id,strand)
                ][order(Cohort,-cnt.sample)]
    dt.circRNA.cohort[,.(length(unique(circRNA_ID)),.N),Cohort]

    if(TRUE){
        # this way, circRNA1 may be redundant
        # circRNA1, GeneA|GeneB => circRNA1, GeneA 
        #                          circRNA1, GeneB
        dt.circRNA.clean<-dt.circRNA[,lapply(.SD,function(i)strsplit(i,",")[[1]]),by=.(circRNA_ID,chr,circRNA_start,circRNA_end,circRNA_type,strand,sum.junc.reads,sum.non.junc.read,avg.junc.ratio,freq,cnt.sample),.SDcol="gene_id"] # split gene name by "," 
        dt.circRNA.clean[gene_id=="n/a",gene_id:=NA]
        setnames(dt.circRNA.clean,"gene_id","ensembl_gene_id") # replace 'n/a' with NA
    }else{
        # below does not deal with multiple gene names (i.e. gene_id may have GeneA,GeneB seperated by comma)
        dt.circRNA.clean[,ensembl_gene_id:=lapply(.SD,function(i) ifelse(i=="n/a",NA,gsub(",$","",i))),.SDcol=c("gene_id")] # remove last ','
        dt.circRNA.clean[,c("gene_id","N"):=NULL] # remove two columns
    }

    gr.circRNA<-makeGRangesFromDataFrame(dt.circRNA,keep.extra.column=T)
    genome(seqinfo(gr.circRNA))<-"GRCh38"
    save(dt.fg.jd, dt.pops, dt.pops.ciri, dt.pops.ciri.all, dt.circRNA, gr.circRNA, dt.circRNA.cohort, dt.circRNA.clean, file=my.RData)

}# end of RData

#################################
## circRNA from PolyA+ RNA-Seq ##
## & circRNA-PT30              ##
## see 'circRNA.polyA.open.data.R'
#################################
if(TRUE){
    load("~/results/RNA-Seq/Placentome/RData/dl.circRNA.polyA.RData")
    names(dl.circRNA.polyA)
    dt.circRNA.PA<-dl.circRNA.polyA[["FGmRNA2014"]]

    dt.circRNA.PA.clean<-dt.circRNA.PA[,lapply(.SD,function(i)strsplit(i,",")[[1]]),by=.(circRNA_ID,chr,circRNA_start,circRNA_end,circRNA_type,strand,sum.junc.reads,sum.non.junc.read,avg.junc.ratio,freq,cnt.sample),.SDcol="gene_id"] # split gene name by "," 
    dt.circRNA.PA.clean[gene_id=="n/a",gene_id:=NA]
    setnames(dt.circRNA.PA.clean,"gene_id","ensembl_gene_id") # replace 'n/a' with NA

    dt.circRNA.PA.clean[,avg.junction.reads:=sum.junc.reads/cnt.sample]

    #write.csv(merge(dt.circRNA.PA.clean, dt.ensg[,.(ensembl_gene_id,hgnc_symbol,`ensembl.strand`=strand,gene_biotype)] ,all.x=T)[order(-avg.junc.ratio)], 
    #            file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA_PolyA+.GRCh38.anno.csv.gz"), row.names=F, quote=F)


    dt.circRNA[,avg.junction.reads:=sum.junc.reads/cnt.sample]
    dt.circRNA[,gene_id:=gsub(",$","",gene_id)]
    dt.circRNA[,gene_id:=gsub(",","|",gene_id)] 
    dt.circRNA[,`in_polyA`:=ifelse(circRNA_ID %in% dt.circRNA.PA$circRNA_ID,TRUE,FALSE)]
    dt.circRNA[freq>=0.3,.N,`in_polyA`]
    #write.csv(dt.circRNA[freq>=0.3], file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA_PT30.GRCh38.csv.gz"), row.names=F, quote=F)

    dt.circRNA.clean[,avg.junction.reads:=sum.junc.reads/cnt.sample]
    dt.circRNA.clean[,`in_polyA`:=ifelse(circRNA_ID %in% dt.circRNA.PA$circRNA_ID,TRUE,FALSE)]
    dt.circRNA.clean[,.N,`in_polyA`]

    # below to share from ShinyPlacentome
    if(F){
        load("~/results/RNA-Seq/Boy.Girl.FG.JD.circRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
        dt.foo<-dt.deseq.anno[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID,.(circRNA_ID,meanCount,meanFpm)][order(-meanCount)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)] # n=3399; n=3279
        write.csv(
            foo<-merge(
                merge(dt.circRNA.clean[freq>=.3 & is_rt==F,-c("is_rt","is_known")], 
                    dt.ensg[,.(ensembl_gene_id,hgnc_symbol,`ensembl.strand`=strand,gene_biotype)] ,all.x=T)[order(-avg.junc.ratio)], 
                merge(dt.foo[,.(circRNA_ID,`norm.bsj.read`=meanCount,RANK)], dt.circbase.hg38[,.(circRNA_ID,`circbase_ID`=name)], all.x=T),
                all.x=T,
                by="circRNA_ID"
                )[order(RANK)]
            ,file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA_PT30.GRCh38.anno.csv.gz"), row.names=F, quote=F
        )
   }
}

#########################
## readthrough circRNA ##
#########################
if(TRUE){
    setkeyv(dt.ensg,c("chromosome_name","strand","start_position","end_position"))
    dt.foo<-foverlaps(dt.circRNA[circRNA_type %in% c("intergenic_region","intron")], dt.ensg, by.x=c("chr","strand","circRNA_start","circRNA_end"),nomatch=0L)

    dt.rt.circ<-dt.foo[!grepl("MIR",hgnc_symbol),.(.N,`genes`=paste(hgnc_symbol,collapse=":")),.(circRNA_ID,circRNA_type,sum.junc.reads,sum.non.junc.read,avg.junc.ratio,freq,avg.junction.reads)][N>1][order(-N,circRNA_ID)]
    write.csv(dt.rt.circ, gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/read-through-circRNA.csv.gz"), quote=F, row.names=F)

    dt.rt.circ[,.(.N,length(unique(circRNA_ID))),circRNA_type]
    dt.rt.circ[circRNA_ID %in% dt.circRNA[in_polyA==FALSE]$circRNA_ID & freq>=.3,.N,circRNA_type]
    
    dt.rt.circ.POPS30<-dt.rt.circ[circRNA_ID %in% dt.circRNA[in_polyA==FALSE]$circRNA_ID & freq>=.3][order(circRNA_type,circRNA_ID)]
    write.csv(dt.rt.circ.POPS30,file="~/results/RNA-Seq/Placentome/circRNA/CSV/read-through-circRNA.POPS30.csv", quote=F, row.names=F)

    ## now using method2 to flag read-through circRNA
    dt.circRNA[,`is_rt`:=ifelse(circRNA_ID %in% dt.rt.circ$circRNA_ID,TRUE,FALSE)]
    dt.circRNA.clean[,`is_rt`:=ifelse(circRNA_ID %in% dt.rt.circ$circRNA_ID,TRUE,FALSE)]

    dt.rt.circ.POPS30[,.N,circRNA_type]
    dt.circRNA[circRNA_ID=="X:7514882|7516290"]
}

###############################
# previously reported circRNA #
###############################
if(TRUE){
    ################################################################ 
    ## run LifeOver tool to convert hg19 to grch38 coordinate     ##
    ## read ~/data/Annotation/circBase/Homo_sapiens/GRCh37/README ##
    ################################################################ 
    if(TRUE){
        suppl_excel <- "/home/ssg29/data/Annotation/circBase/Homo_sapiens/GRCh37/Maass.et.al/109_2017_1582_MOESM3_ESM.xlsx"
        sheet_names <- readxl::excel_sheets(suppl_excel)
        # Read GRCh38 BED files
        my.bed.files<-tools::list_files_with_exts(file.path("~/data/Annotation/circBase/Homo_sapiens/GRCh37/Maass.et.al") , "GRCh38.bed")
        # make a bed file list
        foo<-data.table(`tissue`=sheet_names[order(sheet_names)], `bed.file`=my.bed.files[order(my.bed.files)])
        # import bed to GRangeList
        gl.maass.hg38<-lapply(foo$bed.file, rtracklayer::import.bed)
        names(gl.maass.hg38)<-foo$tissue
        # GRangeList to data.table
        dt.maass.hg38<-rbindlist(lapply(names(gl.maass.hg38), function(i) data.table(`tissue`=i,data.frame(gl.maass.hg38[[i]]))))
        dt.maass.hg38[,.N,seqnames]
        dt.maass.hg38[,seqnames:=substr(seqnames,4,5)] # remove 'chr' prefix: eg chrX => X 
        dt.maass.hg38[seqnames=="M",seqnames:="MT"] # M => MT
        dt.maass.hg38[,circRNA_ID:=paste0(seqnames,':',start,'|',end)] 
        dt.maass.hg38[,Source:="Maass"] 

        dt.maass.hg38[seqnames=="X" & start==140783175 & end==140784659] # CDR1-AS
        dt.maass.hg38[seqnames=="X" & start==7514882 & end==7516290] # novel sponge? 

        #chr19 cluster
        dt.confident<-dt.circRNA[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID & freq>=0.3]
        merge(
            dt.confident[chr=="19" & circRNA_start>=42729778 & circRNA_end<=43262138], # n=63
            dt.maass.hg38, 
            by.x=c("chr","circRNA_start","circRNA_end"),by.y=c("seqnames","start","end")
        )
        dt.maass.hg38[seqnames=="19" & start>=42729778 & end<=43262138]
        dt.circRNA[circRNA_ID=="19:43175870|43194594"]
    }

    ###############################
    ## 2. import circBase GRCh38 ##
    ###############################
    if(TRUE){
        dt.circbase.hg38<-data.table(`tissue`="circBase",data.frame(rtracklayer::import.bed("~/data/Annotation/circBase/Homo_sapiens/GRCh38/hsa_hg38.merged.uniq.bed")))
        dt.circbase.hg38[,.N,seqnames]
        dt.circbase.hg38<-dt.circbase.hg38[seqnames %in% my.ucsc.chr] # my.ucsc.chr defined in Annotation.R
        dt.circbase.hg38[,seqnames:=substr(seqnames,4,5)] # remove 'chr' prefix: eg chrX => X 
        dt.circbase.hg38[seqnames=="M",seqnames:="MT"] # M => MT
        dt.circbase.hg38[,circRNA_ID:=paste0(seqnames,':',start,'|',end)] 
        dt.circbase.hg38[,Source:="circBase"] 

        dt.circbase.hg38[seqnames=="X" & start==140783175 & end==140784659] # CDR1-AS
        # hsa_hg19_Memczak2013.bed:chrX  139865339       139866824       hsa_circ_0001946
        # hsa_hg19_Rybak2015.bed:chrX    139865339       139866824       hsa_circ_0001946
        dt.circbase.hg38[seqnames=="X" & start==7514882 & end==7516290] # novel sponge? 
        #hsa_hg19_Rybak2015.bed:chrX    7432922 7434331 hsa_circ_0140572

        #chr19 cluster (n=0)
        merge(
            dt.confident[chr=="19" & circRNA_start>=42729778 & circRNA_end<=43262138], # n=63
            dt.circbase.hg38,
        )
        dt.circbase.hg38[seqnames=="19" & start>=42729778 & end<=43262138]
    }

    ###################
    ## 3. MiOncoCirc ##
    ###################
    if(TRUE){
        # from the cell journal web site
        dt.MioncoCirc.hg38<-data.table(readxl::read_excel("~/data/Annotation/MiOncoCirc/1-s2.0-S0092867418316350-mmc2.xlsx"))
        dt.MioncoCirc.hg38[,chr:=substr(chr,4,5)]
        dt.MioncoCirc.hg38[chr=="M",chr:="MT"] # M => MT
        dt.MioncoCirc.hg38[,start:=start+1] #  this is BED format
        dt.MioncoCirc.hg38[,circRNA_ID:=paste0(chr,':',start,'|',end)] 
        dt.MioncoCirc.hg38[,Source:="MiOncoCirc"] 
        dt.MioncoCirc.hg38[,.N,`Found in circbase`]
        dt.MioncoCirc.hg38[,.N,circRNA_ID]
        dt.MioncoCirc.hg38[circRNA_ID %in% c("X:140783175|140784659", "X:7514882|7516290")]

        #chr19 cluster (n=0)
        merge(
            dt.confident[chr=="19" & circRNA_start>=42729778 & circRNA_end<=43262138], # n=63
            dt.MioncoCirc.hg38, 
            by.x=c("chr","circRNA_start","circRNA_end"),by.y=c("chr","start","end")
        )
        dt.MioncoCirc.hg38[chr=="19" & start>=42729778 & end<=43262138]
    }
} # end of known and reportec circRNA

############################################
## combined Maass + circBase + MiOncoCirc ##
############################################
if(TRUE){
    dt.known.circRNA<-rbind(
                        unique(dt.maass.hg38[!tissue %in% c('ADA-SCID','LCL','WAS'),.(Source,chr=seqnames,start,end,strand,circRNA_ID)]),
                        unique(dt.circbase.hg38[,.(Source,chr=seqnames,start,end,strand,circRNA_ID)]),
                        unique(dt.MioncoCirc.hg38[,.(Source,chr,start,end,strand='*',circRNA_ID)])
                        )
    dt.known.circRNA[,.N,Source]
    dt.known.circRNA<-dt.known.circRNA[chr%in%my.ens.chr]
    dt.known.circRNA[,.N,chr]
    dt.known.circRNA[,.N,circRNA_ID]

    dt.circRNA[,`is_known`:=ifelse(circRNA_ID %in% dt.known.circRNA$circRNA_ID,TRUE,FALSE)]
    dt.circRNA.clean[,`is_known`:=ifelse(circRNA_ID %in% dt.known.circRNA$circRNA_ID,TRUE,FALSE)]

    dt.circRNA[`in_polyA`==FALSE & freq>=.3 & is_rt==FALSE,.N,is_known]
    ##################################
    # POPS placenta-specific circRNA #
    ##################################
    dt.novel.circRNA<-merge(
                            dt.circRNA.clean[is_known==FALSE & in_polyA==FALSE & is_rt==FALSE]
                            ,dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)],by="ensembl_gene_id",all.x=T
                            )
    write.csv(dt.novel.circRNA, file=gzfile("~/results/RNA-Seq/Placentome/circRNA/CSV/novel.circRNA_PT.txt.gz"), row.names=F, quote=F)
    dt.novel.circRNA[freq==1][order(chr,circRNA_start)]

    cnt.novel.circRNA<-lapply(c(0,seq(0.1,1,by=0.2),1), function(i) dt.novel.circRNA[freq>=i,length(unique(circRNA_ID))])
    names(cnt.novel.circRNA)<-c(0,seq(0.1,1,by=0.2),1)

    dt.novel.circRNA.multi.hosts<-merge(
            merge(
                merge(
                    merge(
                        dt.novel.circRNA[freq>=.3,.(`N_30`=length(unique(circRNA_ID))),.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)][N_30>=5]
                        ,dt.novel.circRNA[freq>=.5,.(`N_50`=length(unique(circRNA_ID))),ensembl_gene_id],by="ensembl_gene_id"
                        )
                    ,dt.novel.circRNA[freq>=.7,.(`N_70`=length(unique(circRNA_ID))),ensembl_gene_id],by="ensembl_gene_id",all.x=T
                    )
                ,dt.novel.circRNA[freq>=.9,.(`N_90`=length(unique(circRNA_ID))),ensembl_gene_id],by="ensembl_gene_id",all.x=T
                )
            ,dt.novel.circRNA[freq>=1,.(`N_100`=length(unique(circRNA_ID))),ensembl_gene_id],by="ensembl_gene_id",all.x=T
            )[!is.na(ensembl_gene_id)][order(-N_30)]

    write.csv(dt.novel.circRNA.multi.hosts, file="~/results/RNA-Seq/Placentome/circRNA/CSV/gene.hosting.more.than.five.novel.circRNA.csv", quote=F,row.names=F)
} # end of setting is_known


########################
## BED File of cirRNA ##
########################
if(TRUE){
    dt.circRNA # 51517
    gr.foo.circRNA<-makeGRangesFromDataFrame(dt.circRNA[in_polyA==FALSE & is_rt==FALSE,-c("in_polyA","is_rt","is_known")],keep.extra.column=T) #49545
    genome(seqinfo(gr.foo.circRNA))<-"GRCh38" # e.g. X, MT
    names(mcols(gr.foo.circRNA))<-c(names(mcols(gr.foo.circRNA))[1:6], "score", "cnt.sample","avg.junction.reads") # freq=>score
    #names(mcols(gr.foo.circRNA))<-c(c("name",names(mcols(gr.foo.circRNA))[-1])[1:6], "score", "cnt.sample","avg.junction.reads") # freq=>score
    gr.foo.circRNA$score<-round(gr.foo.circRNA$score*1000)
    gr.foo.circRNA[gr.foo.circRNA$score>=300,] # POPS30 (i.e. observed > 30% of samples)


    rtracklayer::export.bed(gr.foo.circRNA,con="~/results/RNA-Seq/Placentome/circRNA/circRNA.CIRI2.GRCh38.bed")
    rtracklayer::export.gff(gr.foo.circRNA,con="~/results/RNA-Seq/Placentome/circRNA/circRNA.CIRI2.GRCh38.gtf")
    rtracklayer::export.gff(gr.foo.circRNA[gr.foo.circRNA$score>=300,],con="~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS30.CIRI2.GRCh38.gtf")
    rtracklayer::export.bed(gr.foo.circRNA[gr.foo.circRNA$score>=300,],con="~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS30.CIRI2.GRCh38.bed")

    #UCSC-style
    gr.foo<-gr.foo.circRNA[gr.foo.circRNA$score>=300,]
    GenomeInfoDb::seqlevelsStyle(gr.foo) = "UCSC"  # e.g. X => chrX
    rtracklayer::export.gff(gr.foo, con="~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS30.CIRI2.GRCh38.UCSC.gtf")
    rtracklayer::export.gff(gr.foo, con="~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS30.CIRI2.GRCh38.UCSC.gff3")
    rtracklayer::export.bed(gr.foo, con="~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS30.CIRI2.GRCh38.UCSC.bed")
}


############################
## Excel file for sharing ##
############################
if(FALSE){
    lapply(dl.circRNA.polyA, nrow)

    dt.circRNA.PA[,.(length(unique(circRNA_ID)),.N),chr][order(N)]
    dt.circRNA.PA[freq>=.3 & chr==19 & circRNA_start>=42729778 & circRNA_end<=43262138,.N,circRNA_ID] # n=5

    foo<-list()
    ## sample freq & No of circ
    foo[["No of circRNA"]]=data.table(
        `min freq`=seq(0,1,by=.1),
        `ribo-`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA[freq>=i])),
        `polyA+`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA.PA[freq>=i])),
        ## overlap
        `common`=sapply(seq(0,1,by=.1), function(i) 
            nrow(merge(dt.circRNA[freq>=i], dt.circRNA.PA[freq>=i], by=c("circRNA_ID","chr","circRNA_start","circRNA_end","strand")))
            )
        )

    ## chr19
    foo[["chr19"]]=data.table(
        `min freq`=seq(0,1,by=.1),
        `ribo-`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA[chr=="19" & freq>=i])),
        `polyA+`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA.PA[chr=="19" & freq>=i])),
        ## overlap
        `common`=sapply(seq(0,1,by=.1), function(i) 
            nrow(merge(dt.circRNA[chr=="19" & freq>=i], dt.circRNA.PA[chr=="19" & freq>=i], by=c("circRNA_ID","chr","circRNA_start","circRNA_end","strand")))
            )
        )

    ## chr19 cluster (chr19:42729778-43262138, -strand)
    my.cl=c(42729778,43262138)
    foo[["chr19 cluster"]]<-data.table(
        `min freq`=seq(0,1,by=.1),
        `ribo-`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA[chr=="19" & circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=i])),
        `polyA+`=sapply(seq(0,1,by=.1), function(i) nrow(dt.circRNA.PA[chr=="19" & circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=i])),
        ## overlap
        `common`=sapply(seq(0,1,by=.1), function(i) 
            nrow(merge(dt.circRNA[chr=="19"& circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=i], dt.circRNA.PA[chr=="19"& circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=i], by=c("circRNA_ID","chr","circRNA_start","circRNA_end","strand")))
            )
        )

    ## chr19 cluster (chr19:42729778-43262138, -strand)
    my.cl=c(42729778,43262138)
    foo[["chr19 cluster (POPS30 riboZ)"]]<-data.table(
        `min sample`=c(seq(1,5),seq(10,60,by=5)),
        `min freq`=round(c(seq(1,5),seq(10,60,by=5))/60,2),
        `ribo-(POPS30)`=rep(nrow(dt.circRNA[chr=="19" & circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=.3]),length(c(seq(1,5),seq(10,60,by=5)))),
        `polyA+`=sapply(c(seq(1,5),seq(10,60,by=5)), function(i) nrow(dt.circRNA.PA[chr=="19" & circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & cnt.sample>=i])),
        `common`=sapply(c(seq(1,5),seq(10,60,by=5)), function(i) 
            nrow(merge(dt.circRNA[chr=="19"& circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & freq>=.3], dt.circRNA.PA[chr=="19"& circRNA_start>=my.cl[1] & circRNA_end<=my.cl[2] & cnt.sample>=i], by=c("circRNA_ID","chr","circRNA_start","circRNA_end","strand")))
            )
        )

    foo
    openxlsx::write.xlsx(foo, file ="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.polyA.riboZ.xlsx", creator="Sung Gong")
}


#############################################
## 2 circRNA in chrX that could be sponges ##
#############################################
dt.circRNA[circRNA_ID %in% c("X:140783175|140784659", "X:7514882|7516290")]
dt.circRNA.PA[circRNA_ID %in% c("X:140783175|140784659", "X:7514882|7516290")]

############################
## No. of circRNA vs. FPKM #
############################
## Question: correlation between the abundance of circRNA and the expression level of their hosting genes
if(TRUE){
    #load("~/results/RNA-Seq/Boy.Girl.FG.JD.circRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData") # circRNA
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData") # by total RNA-Seq quantified by featureCount

    dt.circRNA.total<-merge(dt.circRNA.clean[freq>=.3 & in_polyA==FALSE & is_rt==FALSE],dt.deseq.anno)
    # avg.junc.reads vs FPKM 
    cor.test(dt.circRNA.total[!(`in_polyA`),sum.junc.reads/cnt.sample], dt.circRNA.total[!(`in_polyA`),meanFpkm], method="pearson")
    # avg. non-junc reads vs FPKM
    cor.test(dt.circRNA.total[!(`in_polyA`),sum.non.junc.read/cnt.sample], dt.circRNA.total[!(`in_polyA`),meanFpkm], method="pearson")
    # avg.junc.ratio (i.e. circular over linear ratio) vs FPKM
    cor.test(dt.circRNA.total[!(`in_polyA`),avg.junc.ratio], dt.circRNA.total[!(`in_polyA`),meanFpkm], method="pearson")

    dt.foo<-merge(dt.circRNA.clean[freq>=.3,.N,ensembl_gene_id], dt.deseq.anno)
    # number of circRNAs and the FPKM of the hosting genes
    cor.test(dt.foo$N, dt.foo$meanFpkm, method="pearson")

    dt.foo<-merge(dt.circRNA.clean[freq>=.3,.N,ensembl_gene_id], dt.deseq.anno)
    p1<-ggplot(dt.foo, aes(log10(meanFpkm),N)) + 
        geom_point(alpha=.5,size=3) + 
        theme_Publication()
    print(p1)

    dt.circRNA.total[,cor(avg.junc.ratio,meanFpkm,method="pearson")]

    summary(mgcv::gam(formula=y ~ s(x,bs="cs"),data=dt.circRNA.total[meanFpkm>0,.(y=log(avg.junction.reads),x=log(meanFpkm))]))
    foo<-dt.circRNA.total[,cor.test(avg.junction.reads,meanFpkm,method="pearson")]
    foo$p.value; foo$estimate
    p.cor1<-ggplot(dt.circRNA.total, aes(log(meanFpkm),log(avg.junction.reads) )) +
        geom_point(shape=1,size=4) +
        xlab("log(FPKM)") + ylab("log(Number of back-spliced reads)") +
        #geom_smooth(se=F,size=2.5) +
        #geom_line(stat="smooth",method = "gam", formula = y ~ s(x,bs="cs"),col='blue', size =2.5,alpha = 0.7) +
        ggtitle(paste("Pearson coefficient=",round(foo$estimate,3),"(P<",round(foo$p.value,4),")")) +
        theme_Publication() + theme(plot.title = element_text(face = "bold",size = rel(1.1), hjust = 0.5))


    summary(mgcv::gam(formula=y ~ s(x,bs="cs"),data=dt.circRNA.total[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID & meanFpkm>0,.(x=avg.junc.ratio,y=log(meanFpkm))]))
    foo<-dt.circRNA.total[,cor.test(avg.junc.ratio,meanFpkm,method="pearson")]
    foo$p.value; foo$estimate
    p.cor2<-ggplot(dt.circRNA.total,aes(log(meanFpkm),avg.junc.ratio)) +
        geom_point(shape=1,size=4) +
        xlab("log(FPKM)") + ylab("Circular junction read ratio") +
        #geom_smooth(se=F,size=2.5) +
        #geom_line(stat="smooth",method = "gam", formula = y ~ s(x,bs="cs"),col='blue', size =2.5,alpha = 0.7) +
        #ggtitle(paste("Pearson coefficient=",round(foo$estimate,3),"(P<",round(foo$p.value,4),")")) +
        ggtitle(paste("Pearson coefficient=",round(foo$estimate,3),"(P<",1.5e-16,")")) +
        theme_Publication() + theme(plot.title = element_text(face = "bold",size = rel(1.1), hjust = 0.5))

    file.name<-file.path("~/results/RNA-Seq/Placentome/circRNA/avg.junc.ratio.FPKM.correlation.tiff")
    tiff(filename=file.name,width=10, height=6,units="in",res=300, compression = 'lzw') #A4 size 
    cowplot::plot_grid(p.cor1, p.cor2,labels="auto",label_size=20,align="v")
    dev.off()

    ##
    ##
    p.a<-cowplot::plot_grid(NULL,p.frac, NULL, nrow=1, labels=c("a","",""), label_size=27,rel_widths=c(1,7,2))
    p.bc<-cowplot::plot_grid(p.cor1, p.cor2, labels=c("b","c"),label_size=27,align="v")

    file.name<-file.path("~/results/RNA-Seq/Placentome/circRNA/bsj.location.avg.junc.ratio.FPKM.correlation.tiff")
    tiff(filename=file.name,width=10, height=9,units="in",res=300, compression = 'lzw') #A4 size 
    cowplot::plot_grid(p.a, p.bc, nrow=2, labels=c("a",""),label_size=27, rel_heights=c(2,3))
    dev.off()
}

##########################
## Exploratory Analysis ##
##########################
if(TRUE){
    my.file.name<-"~/results/RNA-Seq/Placentome/circRNA/circRNA.number.CIRI2"
    pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title="circRNA Number") # A4 size

    ##########################################
    # 1. NO. of circRNA by junction and freq #
    ##########################################
    i<-0; dl.dummy<-list()
    for(MIN.JNC.CNT in seq(2,10)){ # 2~10
        for(MIN.FREQ in seq(0,1,by=0.1)){
            i<-i+1
            cnt<-dt.pops.ciri[`#junction_reads`>=MIN.JNC.CNT,.(`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops)),circRNA_ID][freq>=MIN.FREQ,.N]
            dl.dummy[[i]]<-data.table(`min junction count`=MIN.JNC.CNT,`min freq`=MIN.FREQ,`Number of circRNA`=cnt)
        }
    }

    if(TRUE){
        p<-ggplot(rbindlist(dl.dummy), aes(`min freq`,`Number of circRNA`,group=`min junction count`)) + 
            geom_point(aes(color=factor(`min junction count`)),size=5,alpha=.7) + 
            geom_line(linetype="dotted",size=.8) +
            scale_y_continuous(breaks=seq(0,2e4,by=1000)) +
            scale_x_continuous(name="Minimum Sample Frequency",breaks=seq(0.1,1,by=0.1)) +
            scale_color_discrete(name="Minimum\nJunction Read") +
            theme_Publication() +
            theme(legend.position=c(.9,.9), axis.text.x = element_text(angle=45,vjust=0.2))
        print(p)

        p<-ggplot(rbindlist(dl.dummy)[`min junction count`==2], aes(`min freq`,`Number of circRNA`)) + 
            geom_point(size=5,alpha=.7) + 
            geom_line(linetype="dotted",size=.8) +
            scale_y_continuous(trans = log10_trans(),
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))) +
            scale_x_continuous(name="Sample Frequency",breaks=seq(0,1,by=0.1),labels=paste0(">=",seq(0,100,by=10),"%")) +
            scale_color_discrete(name="Minimum\nJunction Read") +
            theme_Publication() +
            theme(axis.text.x = element_text(angle=45,hjust=1))
        print(p)

        p<-ggplot(rbindlist(dl.dummy)[`min junction count`==2], aes(`min freq`,`Number of circRNA`)) + 
            geom_point(size=5,alpha=.7) + 
            geom_line(linetype="dotted",size=.8) +
            scale_y_continuous(breaks=seq(0,2e4,by=1000)) +
            scale_x_continuous(name="Minimum Sample Frequency",breaks=c(0,1/100,seq(0.1,1,by=0.1),1)) +
            scale_color_discrete(name="Minimum\nJunction Read") +
            theme_Publication() +
            theme(legend.position=c(.9,.9), axis.text.x = element_text(angle=45,vjust=0.2))
        print(p)
    }

    ###################
    # 2. circRNA_type #
    ###################
    if(TRUE){
        dl.circRNA.cnt<-lapply(c(0,seq(0.1,1,by=0.2),1), function(i) dt.circRNA[in_polyA==F & is_rt==F & freq>=i,.(freq=as.character(i),.N,`total`=nrow(dt.circRNA[in_polyA==F & is_rt==F & freq>=i]),`pct`=round(.N/nrow(dt.circRNA[in_polyA==F & is_rt==F & freq>=i])*100,2)),circRNA_type])
        dt.circRNA.cnt<-rbindlist(dl.circRNA.cnt)
        dt.circRNA.cnt$`circRNA_type`<-factor(dt.circRNA.cnt$`circRNA_type`, levels=c("intergenic_region","intron","exon")) #character to factor  
        p.frac<-ggplot(dt.circRNA.cnt, aes(freq,pct)) +
            geom_bar(aes(fill=circRNA_type),stat="identity",alpha=.8) +
            scale_fill_grey(name="Genomic location of\nback-splice junction", labels=c('intergenic','intron','exon')) +
            scale_x_discrete(name="Sample frequency",breaks=c(0,seq(0.1,1,by=0.2),1),labels=c(">0%",paste0(">=",seq(10,100,by=20),"%"),"100%")) +
            xlab("Sample frequency") + ylab("Proportion (%)") +
            theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))
        file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/TranscriptNumber/circRNA.number.by.bsj.position.tiff")
        tiff(file.name, width=6, height=5,units="in",res=300, compression = 'lzw') #A4 size
        print(p.frac)
        dev.off()

        dcast.data.table(rbindlist(dl.circRNA.cnt), circRNA_type~freq, value.var=c("N"))
        write.csv(dcast.data.table(rbindlist(dl.circRNA.cnt), circRNA_type~freq, value.var=c("N")), file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.cnt.by.type.csv",row.names=F,quote=F)
    }

    ################################
    # 3. Hosting Genes of circRNA ##
    ################################
    if(TRUE){
        # first clean the data
        MIN.FREQ=0.3 # 50% out of 324 samples (at least 162 samples)
        # print
        dt.circRNA.clean[freq>=MIN.FREQ & in_polyA==F & is_rt==F,length(unique(circRNA_ID))] # No. of unique circRNA

        # multi-hosting genes
        dt.circRNA.ensg<-merge(dt.circRNA.clean[`in_polyA`==F & is_rt==F], dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,strand)],by="ensembl_gene_id")
        dt.circRNA.ensg[freq>=MIN.FREQ,.N,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)][N>=10][order(-N)]
        my.ensg<-dt.circRNA.ensg[freq>=MIN.FREQ,.N,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)][N>=10,ensembl_gene_id]

        dt.multiC<-merge(
                merge(
                    merge(
                        merge(
                            dt.circRNA.ensg[freq>=.3,.(`N_30`=.N),.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)][N_30>=10],
                            dt.circRNA.ensg[ensembl_gene_id %in% my.ensg & freq>=.5,.(`N_50`=.N),ensembl_gene_id],by="ensembl_gene_id"
                            ),
                        dt.circRNA.ensg[ensembl_gene_id %in% my.ensg & freq>=.70,.(`N_70`=.N),ensembl_gene_id],by="ensembl_gene_id"
                        ),
                    dt.circRNA.ensg[ensembl_gene_id %in% my.ensg & freq>=.90,.(`N_90`=.N),ensembl_gene_id],by="ensembl_gene_id"
                    ),
                dt.circRNA.ensg[ensembl_gene_id %in% my.ensg & freq>=1,.(`N_100`=.N),ensembl_gene_id],by="ensembl_gene_id"
                )[order(-N_30)]
        write.csv(dt.multiC
            ,file=paste0("~/results/RNA-Seq/Placentome/circRNA/CSV/gene.hosting.at.least.10.circRNAs.at.",MIN.FREQ,".freq.-PA.csv"), row.names=F, quote=F)

        # yet another way
        if(FALSE){
            dt.confident<-dt.circRNA.clean[in_polyA==F & is_rt==F & freq>=MIN.FREQ]
            dt.multi.hosts<-merge(
                                merge(dt.confident, dt.confident[!is.na(ensembl_gene_id),.(`N_circRNA`=.N),.(ensembl_gene_id)][N_circRNA>1])
                                ,dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,strand)])[order(-N_circRNA)]
            write.csv(dt.multi.hosts, file=paste0("~/results/RNA-Seq/Placentome/circRNA/CSV/multi.hosts.freq.",MIN.FREQ,".csv"), row.names=F, quote=F)

            dt.multi.hosts[N_circRNA>=10 & freq>=.5,.N,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,N_circRNA)][order(-N_circRNA)]
            dt.multi.hosts[N_circRNA>=10 & freq>=.9,.N,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,N_circRNA)][order(-N_circRNA)]

            merge(dt.novel.circRNA.multi.hosts, dt.multi.hosts[N_circRNA>=10,.N,.(ensembl_gene_id,hgnc_symbol,gene_biotype)], all.x=T)[order(-N)]
            merge(dt.multi.hosts[N_circRNA>=10,.N,.(ensembl_gene_id,hgnc_symbol,gene_biotype)], dt.novel.circRNA.multi.hosts, all=T)[order(-N)]
            dt.multi.hosts[hgnc_symbol %in% c("MYO9A","PSG8","CYP19A1","AC022872.1"),.N,.(ensembl_gene_id,hgnc_symbol,gene_biotype)][order(-N)]
        }

    }

    dev.off()

    dt.confident[,`length`:=circRNA_end-circRNA_start+1]
    lapply(split(dt.confident, dt.confident$circRNA_type), function(i) i[,summary(length)])

    #######
    # 5. 
    #######
    dt.circRNA[freq>=0.3,.(.N,`No_circRNA`=length(unique(circRNA_ID))),chr][order(chr)] # by chr
    dt.circRNA[freq>=0.9,.(.N,`No_circRNA`=length(unique(circRNA_ID))),chr][order(chr)] # by chr
    dt.circRNA[,.N,chr]
}

##########
## Venn ##
##########
if(TRUE){
    # NB. ADA-SCID, LCL, WAS: not tissues
    library(VennDiagram)
    venn.list=list()
    venn.list[["POPS30"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.3,paste(chr,circRNA_start,circRNA_end,strand,sep=":")]
    venn.list[["POPS90"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.9,paste(chr,circRNA_start,circRNA_end,strand,sep=":")]
    venn.list[["Maass et al"]]=dt.maass.hg38[!tissue %in% c('ADA-SCID','LCL','WAS'),paste(seqnames,start,end,strand,sep=":")]
    venn.list[["circBase"]]=dt.circbase.hg38[,paste(seqnames,start,end,strand,sep=":")]

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/Venn/circRNA/circRNA-PT30.PT90.Maass.circBase.venn.a.-PA.tiff")
    venn.diagram(
        x=venn.list,
        filename = file.name,
        height=1500,
        width=1700,
        main="a",
        main.pos=c(0.02,0.98),
        main.cex = 1.1,
        main.fontfamily = "sans",
        main.fontface = "bold",
        col = "black",
        fill = c(cbPalette[1],cbPalette[2],"white","white"),
        alpha = 0.50,
        cex = 0.95,
        fontfamily = "sans",
        cat.cex = 0.85,
        cat.pos=c(-10,8,0,0),
        cat.fontfamily = "sans",
        cat.fontface = "bold",
        margin=0.05
    )

    venn.list=list()
    venn.list[["POPS30"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.3,paste(chr,circRNA_start,circRNA_end,strand,sep=":")]
    venn.list[["POPS90"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.9,paste(chr,circRNA_start,circRNA_end,strand,sep=":")]
    venn.list[["Placenta\n(Maass et al)"]]=dt.maass.hg38[tissue=="placenta",paste(seqnames,start,end,strand,sep=":")]
    venn.list[["Decidua\n(Maass et al)"]]=dt.maass.hg38[tissue=="decidua",paste(seqnames,start,end,strand,sep=":")]

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/Venn/circRNA/circRNA-PT30.PT90.Maass.placenta.decidua.venn.b.-PA.tiff")
    venn.diagram(
        x=venn.list,
        filename = file.name,
        height=1500,
        width=1700,
        main="b",
        main.pos=c(0.02,0.98),
        main.cex = 1.1,
        main.fontfamily = "sans",
        main.fontface = "bold",
        col = "black",
        fill = c(cbPalette[1],cbPalette[2],"white","white"),
        alpha = 0.50,
        cex = 0.95,
        fontfamily = "sans",
        cat.cex = 0.85,
        #cat.pos=c(-10,8,20,10),
        cat.pos=c(-10,8,0,-10),
        cat.fontfamily = "sans",
        cat.fontface = "bold",
        margin=0.05
    )

    venn.list=list()
    venn.list[["POPS30"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.3,paste(chr,circRNA_start,circRNA_end,sep=":")]
    venn.list[["MiOncoCirc"]]=dt.MioncoCirc.hg38[,paste(chr,start,end,sep=":")]
    venn.list[["Maass et al"]]=dt.maass.hg38[!tissue %in% c('ADA-SCID','LCL','WAS'),paste(seqnames,start,end,sep=":")]
    venn.list[["circBase"]]=dt.circbase.hg38[,paste(seqnames,start,end,sep=":")]

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/Venn/circRNA/circRNA-PT30.Maass.circBase.MiOncoCirc.venn.c.-PA.tiff")
    VennDiagram::venn.diagram(
        x=venn.list,
        filename = file.name,
        height=1500,
        width=1700,
        main="c",
        main.pos=c(0.02,0.98),
        main.cex = 1.1,
        main.fontfamily = "sans",
        main.fontface = "bold",
        col = "black",
        fill = c(cbPalette[1],"white","white","white"),
        alpha = 0.50,
        cex = 0.8,
        fontfamily = "sans",
        cat.cex = 0.8,
        cat.pos=c(-10,8,0,0),
        cat.fontfamily = "sans",
        cat.fontface = "bold"
    )

    venn.list=list()
    venn.list[["POPS90"]]=dt.circRNA[in_polyA==F & is_rt==F & freq>=0.9,paste(chr,circRNA_start,circRNA_end,sep=":")]
    venn.list[["MiOncoCirc"]]=dt.MioncoCirc.hg38[,paste(chr,start,end,sep=":")]
    venn.list[["Maass et al"]]=dt.maass.hg38[!tissue %in% c('ADA-SCID','LCL','WAS'),paste(seqnames,start,end,sep=":")]
    venn.list[["circBase"]]=dt.circbase.hg38[,paste(seqnames,start,end,sep=":")]

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/Venn/circRNA/circRNA-PT90.Maass.circBase.MiOncoCirc.venn.d.-PA..tiff")
    VennDiagram::venn.diagram(
        x=venn.list,
        filename = file.name,
        height=1500,
        width=1700,
        main="d",
        main.pos=c(0.02,0.98),
        main.cex = 1.1,
        main.fontfamily = "sans",
        main.fontface = "bold",
        col = "black",
        fill = c(cbPalette[2],"white","white","white"),
        alpha = 0.50,
        cex = 0.8,
        fontfamily = "sans",
        cat.cex = 0.8,
        cat.pos=c(-10,8,0,0),
        cat.fontfamily = "sans",
        cat.fontface = "bold"
    )
}

############################
## Translation of circRNA ##
############################
if(TRUE){
    #gr.circPOPS30<-gr.circRNA[mcols(gr.circRNA)$freq>=0.3]  #n=3452
    dt.circRNA[`in_polyA`==FALSE & freq>=.3 & is_rt==FALSE]
    gr.circPOPS30<-makeGRangesFromDataFrame(dt.circRNA[`in_polyA`==FALSE & freq>=.3 & is_rt==FALSE,-c("in_polyA","is_rt","is_known")], keep.extra.column=T) # n=3279
    seqlevelsStyle(gr.circPOPS30) = "UCSC" 
    names(gr.circPOPS30)<-gr.circPOPS30$circRNA_ID
    genome(seqinfo(gr.circPOPS30))<-"hg38"

    # prepare fasta seq
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38 # chr1, chr2 .. chrM, chrX, chrY
    #library(BSgenome.Hsapiens.NCBI.GRCh38)
    #genome <- BSgenome.Hsapiens.NCBI.GRCh38 # 1, 2 .. MT, X, Y

    # fasta in a single file 
    # but this is just from the start-end of BSJ
    bs.circ<-BSgenome::getSeq(genome, gr.circPOPS30) # isa 'DNAStringSet'
    Biostrings::writeXStringSet(bs.circ, "~/results/RNA-Seq/Placentome/circRNA/circRNA.POPS.30.fa")

    # fasta files per circRNA
    gl.circPOPS30<-split(gr.circPOPS30, names(gr.circPOPS30))
    sapply(names(gl.circPOPS30), function(i){Biostrings::writeXStringSet(getSeq(genome,gl.circPOPS30[[i]]), file.path("~/results/RNA-Seq/Placentome/circRNA/FASTA",paste(i,"fa",sep=".")))})

    ###############################
    ## last 100 nt + first 100 nt #
    ###############################
    bs.circ.200nt<-DNAStringSet(
                        lapply(bs.circ, function(i){
                                BString(
                                    paste0(
                                        # last 100 nt
                                        toString(Biostrings::subseq(i, start=ifelse(length(i)<100,length(i)+1,ifelse(length(i)<200, 101, length(i)-100+1)) , end=length(i))), 
                                        # first 100 nt
                                        toString(Biostrings::subseq(i, start=1, end=ifelse(length(i)<100,length(i),100))) 
                                        )
                                    )
                            }) # isa list 
                        )

    ########################
    ## 3-frame translation #
    ########################
    if(TRUE){
        # USE THIS
        li.circ.3frame <- lapply(names(bs.circ.200nt), function(i){
                                    t<-AAStringSet(lapply(1:3, function(pos) translate(subseq(bs.circ.200nt[[i]], start=pos))))
                                    names(t)<-paste(i,seq(1:3),sep="_")
                                    return(t) # isa 'AAStringSet'
                                }
                            ) # isa 'list'
        names(li.circ.3frame)<-names(bs.circ.200nt)
    }

    sapply(names(li.circ.3frame), function(i){Biostrings::writeXStringSet(li.circ.3frame[[i]], file.path("~/results/RNA-Seq/Placentome/circRNA/POPSTrCircRNA",paste(i,"tr.fa.gz",sep=".")), compress=T)})
} # end of translation of circRNA

########################
## Single-exon circRNA #
########################
if(TRUE){
    dt.sExon<-merge(
                    dt.circRNA[`in_polyA`==FALSE & freq>=0.3 & is_rt==FALSE], # n=3279
                    dt.exon, 
                    by.x=c("chr","circRNA_start","circRNA_end","strand"),
                    by.y=c("chromosome_name","start_position","end_position","strand"))[,-"gene_id"]
    dt.sExon[,.N,.(circRNA_ID)] # n=236
    
    gr.sExon<-gr.circRNA[gr.circRNA$circRNA_ID %in% dt.sExon[,.N,.(circRNA_ID)]$circRNA_ID] # n=236
    seqlevelsStyle(gr.sExon) = "UCSC" 
    names(gr.sExon)<-gr.sExon$circRNA_ID
    genome(seqinfo(gr.sExon))<-"hg38"

    # prepare fasta seq
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38 # chr1, chr2 .. chrM, chrX, chrY

    bs.sExon<-BSgenome::getSeq(genome, gr.sExon) # isa 'DNAStringSet'
    ########################
    ## 3-frame translation #
    ########################
    li.sExon.3frame <- lapply(names(bs.sExon), function(i){
                                t<-AAStringSet(lapply(1:3, function(pos) translate(subseq(bs.sExon[[i]], start=pos))))
                                names(t)<-paste(i,seq(1:3),sep="_")
                                return(t) # isa 'AAStringSet'
                            }
                        ) # isa 'list'
    names(li.sExon.3frame)<-names(bs.sExon)
    length(li.sExon.3frame) # n=236

    sapply(names(li.sExon.3frame), function(i){Biostrings::writeXStringSet(li.sExon.3frame[[i]], file.path("~/results/RNA-Seq/Placentome/circRNA/POPS30circRNASingleExon3FrameTr",paste(i,"tr.fa.gz",sep=".")), compress=T)})

    #li.sExon.3frame[["5:68226290|68227009"]]
    dt.foo<-merge(dt.sExon, dt.ensg[,.(ensembl_gene_id,gene_biotype)], by="ensembl_gene_id")[,.(circRNA_ID,chr,circRNA_start,circRNA_end,strand,circRNA_type,sum.junc.reads,sum.non.junc.read,avg.junc.ratio,freq,cnt.sample,ensembl_gene_id,ensembl_transcript_id,ensembl_exon_id,exon_number,hgnc_symbol,gene_biotype)]
    write.csv(dt.foo[,-"circRNA_type"], file="~/results/RNA-Seq/Placentome/circRNA/POPS30circRNASingleExon.csv", quote=F,row.names=F)

}


################
## rt-circRNAs #
################
# POPS
# method1
dt.query<-dt.circRNA[freq>=0.3 & circRNA_type=="exon"]
setkeyv(dt.query,c("chr","strand","circRNA_start","circRNA_end"))
#setkeyv(dt.ensg,c("chromosome_name","strand","start_position","end_position"))
#dt.overlap<-foverlaps(dt.query, dt.ensg, nomatch=0L)
#dt.foo<-dt.overlap[gene_biotype=="protein_coding",.(.N,length(unique(ensembl_gene_id))),circRNA_ID]
setkeyv(dt.exon,c("chromosome_name","strand","start_position","end_position"))
dt.overlap<-foverlaps(dt.query, dt.exon, nomatch=0L)
dt.foo<-dt.overlap[,.(.N,length(unique(ensembl_gene_id))),circRNA_ID]
dt.overlap[circRNA_ID %in% dt.foo[V2>1,circRNA_ID]]
dt.overlap[circRNA_ID %in% dt.foo[V2>1,circRNA_ID],.N,circRNA_ID][order(-N)][1:50] 
dt.overlap[circRNA_ID=="15:89871475|89889140"]
my.circRNA_ID<-dt.foo[V2>2,circRNA_ID]


# method 2
dt.foo<-merge(dt.circRNA.clean[freq>=0.3 & circRNA_type=="exon"], dt.exon, by.x=c("ensembl_gene_id","chr","circRNA_start","strand"), by.y=c("ensembl_gene_id","chromosome_name","start_position","strand") )

dt.bar<-merge(dt.foo, dt.exon, by.x=c("chr","circRNA_end","strand"), by.y=c("chromosome_name","end_position","strand") )[as.character(ensembl_gene_id.x)!=as.character(ensembl_gene_id.y)]

merge(dt.foo, dt.exon, by.x=c("chr","circRNA_end","strand"), by.y=c("chromosome_name","end_position","strand") )[as.character(ensembl_gene_id.x)!=as.character(ensembl_gene_id.y),.N,.(circRNA_ID,ensembl_gene_id.x,hgnc_symbol.x,ensembl_gene_id.y,hgnc_symbol.y)]
merge(dt.foo, dt.exon, by.x=c("chr","circRNA_end","strand"), by.y=c("chromosome_name","end_position","strand") )[as.character(ensembl_gene_id.x)!=as.character(ensembl_gene_id.y),.N,.(circRNA_ID)]
my.circRNA_ID<-merge(dt.foo, dt.exon, by.x=c("chr","circRNA_end","strand"), by.y=c("chromosome_name","end_position","strand") )[as.character(ensembl_gene_id.x)!=as.character(ensembl_gene_id.y),.N,.(circRNA_ID)]$circRNA_ID


# method 3
dt.circRNA.clean[freq>=0.3,.(.N,length(unique(ensembl_gene_id))),circRNA_ID][V2>1]
dt.circRNA[circRNA_ID %in% dt.circRNA.clean[freq>=0.3,.(.N,length(unique(ensembl_gene_id))),circRNA_ID][V2>1,circRNA_ID]]
dt.circRNA.clean[circRNA_ID %in% dt.circRNA.clean[freq>=0.3,.(.N,length(unique(ensembl_gene_id))),circRNA_ID][V2>1,circRNA_ID]]

my.circRNA_ID<-dt.circRNA.clean[freq>=0.3,.(.N,length(unique(ensembl_gene_id))),circRNA_ID][V2>1,circRNA_ID]


##
## chrX miR sponge (X:7514882|7516290)
## parsing miranda result 
if(TRUE){
    library(Biostrings)
    bs.circX<-Biostrings::readDNAStringSet("~/results/RNA-Seq/Placentome/circRNA/FASTA/X:7514882|7516290.fa")

    my.cmd1<-"grep '^>hsa-miR-5584-5p' ~/results/RNA-Seq/Placentome/circRNA/miranda/X:7514882\\|7516290.miranda.txt | sed 's/^>//g'"
    my.cmd2<-"grep '^>hsa-miR-7113-5p' ~/results/RNA-Seq/Placentome/circRNA/miranda/X:7514882\\|7516290.miranda.txt | sed 's/^>//g'"
    foo<-system(my.cmd1,intern=T)
    bar<-system(my.cmd2,intern=T)

    dt.foo<-data.table(t(simplify2array(strsplit(c(foo,bar), "\t")))) 
    setnames(dt.foo, c("V1","V2","V3","V4","V7","V8","V9"), c("miRNA","ciRNA","score","energy","alnLen","pid1","pid2"))

    dt.foo[, c("m.start","m.end") := tstrsplit(V5, " ", fixed=TRUE)]
    dt.foo[, c("c.start","c.end") := tstrsplit(V6, " ", fixed=TRUE)]
    dt.spg<-dt.foo[,.(miRNA,score=as.integer(score),energy=as.numeric(energy),
              m.start=as.numeric(m.start),m.end=as.numeric(m.end),
              c.start=as.numeric(c.start),c.end=as.numeric(c.end),
              alnLen=as.integer(alnLen),pid1,pid2)][order(miRNA,c.start)]

    ir.circ<-IRanges(start=1,end=width(bs.circX))
    ir.m1<-IRanges(start=dt.spg[miRNA=="hsa-miR-5584-5p",c.start],end=dt.spg[miRNA=="hsa-miR-5584-5p",c.end])
    ir.m2<-IRanges(start=dt.spg[miRNA=="hsa-miR-7113-5p",c.start],end=dt.spg[miRNA=="hsa-miR-7113-5p",c.end])

    ##
    ##  ggbio
    ##
    library(ggbio)
    g.circ<-autoplot(gr.circX,) # autoplot(ir.circ)
    g.m1<-autoplot(gr.m1) # autoplot(ir.m1)
    g.m2<-autoplot(gr.m2) # autoplot(ir.m2)
    tracks(g.circ,g.m1,g.m2, heights=c(2,1,1)) + theme_bw()

    g.circ<-autoplot(ir.circ, layout="circle")
    g.m1<-autoplot(ir.m1, layout="circle", fill="grey80")
    g.m2<-autoplot(ir.m2, layout="circle", fill="yellow")

    #

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/circRNA/chrX.sponge.miR.binding.circ.tiff")
    tiff(filename=file.name,width=20, height=20,units="cm",res=300, compression = 'lzw') #A4 size 

    ggbio(radius=100) + 
        circle(ir.m1, layout="circle", geom="segment", size=2,col=ggsci::pal_d3("category10")(10)[1]) +
        circle(ir.m2, layout="circle", geom="segment", size=2,col=ggsci::pal_d3("category10")(10)[2]) + 
        circle(ir.circ, layout="circle", geom="segment",size=2) +
        circle(ir.circ, layout="circle", geom="scale", size=3, scale.n=8, scale.type="sci") 

    dev.off()

        
    ###########
    ##  Gviz ##
    ###########
    library(Gviz)
    library(BSgenome.Hsapiens.UCSC.hg38)
    #library(BSgenome.Hsapiens.NCBI.GRCh38)
    #options(ucscChromosomeNames=FALSE)
    seqinfo(Hsapiens)
    #gv.strack <- SequenceTrack(Hsapiens, cex=1.3)

    ## GRanges
    gr.circX<-gr.circRNA[gr.circRNA$circRNA_ID=="X:7514882|7516290",] # NCBI style (e.g. X)
    GenomeInfoDb::seqlevelsStyle(gr.circX) = "UCSC"
    gr.m1<-makeGRangesFromDataFrame(dt.spg[miRNA=="hsa-miR-5584-5p",.(`seqnames`="chrX",`start`=c.start+7514881,`end`=c.end+7514881,score,energy)],keep.extra.columns=T)
    gr.m2<-makeGRangesFromDataFrame(dt.spg[miRNA=="hsa-miR-7113-5p",.(`seqnames`="chrX",`start`=c.start+7514881,`end`=c.end+7514881,score,energy)],keep.extra.columns=T)

    #
    gv.axis <-GenomeAxisTrack(cex=2,add53=T) # littleTicks=T)
    gv.circ<-AnnotationTrack(gr.circX, fill="grey50", id="chrX:7514882-7516290", col.line="black",name="", featureAnnotation="id", fontsize=17, cex.group=1)
    gv.m1<-AnnotationTrack(gr.m1, name="", fill="grey10")
    gv.m2<-AnnotationTrack(gr.m2, name="",fill="cornflowerblue")

    gv.m3<-DataTrack(gr.m1, data=gr.m1$energy, type="heatmap")
    gv.m4<-DataTrack(gr.m2, data=gr.m2$energy, type="heatmap")

    #plotTracks(list(gv.axis, gv.circ, gv.m1,gv.m2, gv.m3,gv.m4), sizes=c(1.3,.3,.5,.5,.5,.5))

    #plotTracks(list(gv.axis, gv.circ, gv.m1,gv.m2 ), sizes=c(1.3,.3,.2,.2))

    plotTracks(list(gv.circ, gv.m1, gv.m2), sizes=c(3,1,1), panel.only=T)

    ##################
    ## Get Sequences # 
    ##################
    Biostrings::subseq(bs.circX, start=1, end=10)
    li.bs.circX<-lapply(split(dt.spg, dt.spg$miRNA), function(i){
                    lapply(unname(split(i,i$c.start)), function(j){
                                    Biostrings::subseq(bs.circX, start=j$c.start, end=j$c.end)
                                })
                }) # list of list of Biostrings 
    names(li.bs.circX)
    do.call(c,li.bs.circX[["hsa-miR-7113-5p"]])
    split(dt.spg, dt.spg$miRNA)[[1]]

    Biostrings::writeXStringSet(do.call(c,li.bs.circX[["hsa-miR-5584-5p"]]), "~/results/RNA-Seq/Placentome/circRNA/hsa-miR-5584-5p.vs.X:7514882|7516290.fa")
    Biostrings::writeXStringSet(do.call(c,li.bs.circX[["hsa-miR-7113-5p"]]), "~/results/RNA-Seq/Placentome/circRNA/hsa-miR-7113-5p.vs.X:7514882|7516290.fa")
}

############################
# chr19 cluster & abudance #
############################
if(TRUE){
    my.bed=file.path("~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.CIRI2.GRCh38.gtf.gz")
    gr.target<-rtracklayer::import(my.bed) # isa 'GRanges'
    dt.target<-data.table(data.frame(gr.target))[,c("seqnames","start","end","strand","circRNA_ID")]
    setnames(dt.target, c("chromosome_name","start_position","end_position","strand","circRNA_ID"))

    load("~/results/RNA-Seq/Boy.Girl.FG.JD.circRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData") # freq>=30% (n=3452)
    #dt.foo<-dt.mean[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID][order(-meanCount)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)] # -polyA+
    dt.foo<-dt.mean[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID & !circRNA_ID %in% dt.rt.circ$circRNA_ID][order(-meanCount)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)] # -polyA+
    dt.circ<-merge(dt.foo, dt.target)[order(RANK)] # n=3399 | n=3322

    this.col<-c("circRNA_ID","meanCount","chromosome_name","start_position","end_position","strand","RANK","INDEX_PCT")
    dt.circ[,..this.col]                                          # n=3452, 3399
    dt.circ[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID,..this.col] # n=3399 (53 were also in polyA+)

    # C19
    dt.circ[chromosome_name==19 & start_position>=40000000 & end_position<=45000000 & INDEX_PCT<=5,..this.col] # n=79, 20
    dt.circ[chromosome_name==19 & start_position>=40000000 & end_position<=45000000 & INDEX_PCT<=5,][,.(min(start_position),max(end_position))]
    dt.circ[chromosome_name==19 & start_position>=42729778 & end_position<=43262138,..this.col]                # n=79,63
    dt.circ[chromosome_name==19 & start_position>=42729778 & end_position<=43262138 & INDEX_PCT<=5,..this.col] # n=25,20
    dt.circ[chromosome_name==19 & start_position>=42729778 & end_position<=43262138 & INDEX_PCT<=1,..this.col] # n=16,12

    # C17
    dt.circ[chromosome_name==17][order(RANK)][1:20]
    dt.ensg[chromosome_name==17 & start_position>=63872577 & end_position<=63895850][order(hgnc_symbol)]

    setkeyv(dt.ensg,c("chromosome_name","strand","start_position","end_position"))
    foverlaps(dt.circ, dt.ensg, by.x=c("chromosome_name","strand","start_position","end_position"),nomatch=0L)[gene_biotype=="protein_coding"]

    #
    dt.circ[circRNA_ID %in% c("X:140783175|140784659", "X:7514882|7516290")]
    #
    merge(dt.circ[RANK<=10,..this.col], dt.known.circRNA, by.x=c("chromosome_name","start_position","end_position"),by.y=c("seqnames","start","end"))[,.N,.(circRNA_ID,RANK)][order(RANK)]
    merge(dt.confident, dt.known.circRNA, by.x=c("chr","circRNA_start","circRNA_end"),by.y=c("seqnames","start","end"))[]

    ## 
    merge(
          merge(dt.foo[,.(circRNA_ID,RANK,INDEX_PCT)],dt.circRNA.clean[,.(circRNA_ID,circRNA_type,chr,circRNA_start,circRNA_end,ensembl_gene_id)]),
          dt.ensg[,.(ensembl_gene_id,hgnc_symbol)],by="ensembl_gene_id"
    )[order(RANK)]

    [order(RANK)][chr=="19" & circRNA_start>=42729778 & circRNA_end<=43262138]


}

#################################################
## Selection of circRNAs to validate by RT-PCR ##
## 2020-01-03
#################################################
if(TRUE){
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.circRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    dt.foo<-dt.deseq.anno[!circRNA_ID %in% dt.circRNA.PA$circRNA_ID,.(circRNA_ID,meanCount,meanFpm)][order(-meanCount)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)] # n=3399; n=3279

    ##
    ## 0. for SI table 
    ##
    dt.circRNA[freq>=.3]
    dt.circRNA[freq>=.3,.N, is_rt]
    dt.circRNA[freq>=.3 & is_rt==F,.N,in_polyA]
    dt.circRNA[freq>=.3 & is_rt==F & in_polyA==F]
    dt.circRNA[freq>=.9 & is_rt==F & in_polyA==F]

    dt.circRNA.SI.table<-merge(
            merge(dt.circRNA.clean[freq>=.3 & in_polyA==F & is_rt==F,-c("in_polyA","is_rt","is_known")], dt.ensg[,.(ensembl_gene_id,hgnc_symbol,gene_biotype)], all.x=T),
            merge(dt.foo[,.(circRNA_ID,`norm.bsj.read`=meanCount,RANK)], dt.circbase.hg38[,.(circRNA_ID,`circbase_ID`=name)], all.x=T),
            by="circRNA_ID"
            )[order(RANK)]
    write.csv(dt.circRNA.SI.table, ,file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.POPS30.with.cirbaseID.csv", quote=F,row.names=F)
    dt.circRNA.SI.table[,.N,circRNA_ID]

    ######################
    ## I. Known circRNA ## 
    ######################
    # 1. hsa_circ_0001946: aka. CDR1as from Memczak 2013 Nature (rank 892)
    # 2. hsa_circ_0000268: aka. hsa_circ_000002 from Memczak 2013 Nature; ZRANB1 exon1 (rank 320)
    # 3. hsa_circ_0001445: aka. hsa_circ_000003 from Memczak 2013 Nature (rank 22)
    # 4. hsa_circ_0000284: aka. hsa_circ_000016 from Memczak 2013 Nature (highly expressed in the placenta 11:33286413|33287511; rank7)
    # 5. hsa_circ_0140572: chrX putative sponge (X:7514882|7516290; rank 816)
    # 6. hsa_circ_0005868: DEG in PE (20:34470048|34481206; rank 1319)
    # 7. hsa_circ_0007444: highly expressed in the placenta (5:95755396|95763620; rank4)
    my.known.circ=c("hsa_circ_0001946","hsa_circ_0000268","hsa_circ_0001445","hsa_circ_0000284","hsa_circ_0140572","hsa_circ_0005868","hsa_circ_0007444")
    dt.known.target<-merge(dt.circRNA, dt.circbase.hg38[name %in% my.known.circ,.(circRNA_ID,`circbase_ID`=name)])


    ###############################
    ## II. Novel in the placenta ##
    ###############################
    dt.novel.target<-rbind(
    # 8-9. highly expressed novel (5:87238815|87239762;  15:51236859|51242950)
    dt.circRNA[circRNA_type=="exon" & avg.junction.reads>20 & freq>=0.7 &  in_polyA==FALSE & is_rt==FALSE & is_known==FALSE & chr!=19][order(-avg.junction.reads)],
    # 10. DEG in PE (19:53687886|53694145)
    dt.circRNA[circRNA_ID=="19:53687886|53694145"]
    # 11. rank 1 (PSG3,PSG8)
    #dt.circRNA[circRNA_ID=="19:42740321|42764281"],
    # 12. rank 2 (PSG3,PSG8,PSG1)
    #dt.circRNA[circRNA_ID=="19:42740321|42878278"],
    # 13. rank 3 (PSG3,PSG8,PSG1,PSG6,PSG7)
    #dt.circRNA[circRNA_ID=="19:42740321|42935769"]
    )[,circbase_ID:='NULL']

    dt.final.target<-merge(merge(rbind(dt.known.target,dt.novel.target), dt.foo), dt.ensg[,.(ensembl_gene_id,hgnc_symbol)],by.x="gene_id",by.y="ensembl_gene_id",all.x=T)[order(-is_known,RANK)]

    # 3prime_100nt + 5prime_100nt
    data.frame(bs.circ.200nt["X:7514882|7516290"])
    foo<-as.data.frame(bs.circ.200nt[dt.final.target$circRNA_ID])
    dt.seq<-data.table(circRNA_ID=rownames(foo), three_prime_100nt=subseq(foo$x,1,100), five_prime_100nt=subseq(foo$x,101,200))

    merge(dt.final.target, dt.seq)[order(-is_known,RANK)]
    write.csv(merge(dt.final.target, dt.seq)[order(-is_known,RANK)],
              file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.validation.candidate.seq.csv", quote=F,row.names=F)


    ## find the best samples to validate the 10 circRNAs
    dt.pops.ciri[circRNA_ID %in% dt.final.target$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169"),.N,circRNA_ID] # n=10
    dt.pops.ciri[circRNA_ID %in% dt.final.target$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169"),.(circRNA_ID,`#junction_reads`,`#non_junction_reads`,junction_reads_ratio,SLX,Barcode)]

    dt.final.samples<-dt.pops.ciri[circRNA_ID %in% dt.final.target$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169"),
                                   .(.N,
                                     mean_nj=mean(`#non_junction_reads`),
                                     median_nj=median(`#non_junction_reads`),
                                     min_nj=min(`#non_junction_reads`),
                                     max_nj=max(`#non_junction_reads`),
                                     mean_bsj=mean(`#junction_reads`),
                                     median_bsj=median(`#junction_reads`),
                                     min_bsj=min(`#junction_reads`),
                                     max_bsj=max(`#junction_reads`))
                                   ,.(SLX,Barcode)][N==10][order(-min_bsj,-median_bsj)]
    write.csv(merge(dt.final.samples,dt.pops[,.('SLX'=Library,'Barcode'=BarCode,SampleName,CRN)])[order(-min_bsj,-median_bsj)],
              file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.validation.candidate.samples.csv", quote=F,row.names=F)

    dt.final.samples.circRNA<-dcast.data.table(
                    merge(
                    dt.pops.ciri[circRNA_ID %in% dt.final.target$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169")]
                    ,dt.pops.ciri[circRNA_ID %in% dt.final.target$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169"),.(.N),.(SLX,Barcode)][N==10]
                    ),
                    SLX+Barcode~circRNA_ID,value.var='#junction_reads'
    )
    write.csv(dt.final.samples.circRNA,file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.validation.candidate.samples.bsj.csv", quote=F,row.names=F)

    ################
    ## -ve samples #
    ################
    dt.foo<-merge(
        dt.pops[Library %in% c("SLX-9168","SLX-9169"),.('SLX'=Library,'Barcode'=BarCode,SampleName,CRN)], # n=107 samples
        #dt.pops.ciri[circRNA_ID %in% dt.final.target[freq!=1]$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169"),.N,circRNA_ID] # n=5 circRNAs
        dt.pops.ciri[circRNA_ID %in% dt.final.target[freq!=1]$circRNA_ID & SLX %in% c("SLX-9168","SLX-9169")],
        all.x=T
        )
    dt.foo[is.na(circRNA_ID)] # all samples have the target circRNA
    dt.foo[,.N,.(SLX,Barcode)] # n=107 samples

    write.csv(
    dcast.data.table(dt.foo,SLX+Barcode+SampleName+CRN~circRNA_ID,value.var="#junction_reads")[is.na(`X:140783175|140784659`) & is.na(`X:7514882|7516290`) & is.na(`20:34470048|34481206`)]
        ,file="~/results/RNA-Seq/Placentome/circRNA/CSV/circRNA.validation.candidate.negative.samples.csv", quote=F,row.names=F)


}
