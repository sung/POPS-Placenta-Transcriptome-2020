##
## Size of small RNAs
## USE THIS
if(TRUE){
    library(DESeq2)
    my.ucsc.chr<-paste0("chr",c(seq(1:22),"X","Y","M"))
    my.ens.chr<-c(seq(1:22),"X","Y","MT")

    # 1. known mature miRNA
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.mature.miRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    dt.mature.miRNA.fpkm<-merge(dt.mean[,.(`type`="mature miRNA",`ID`=mirbase_id,meanCount,meanFpm,meanFpkm)],
                        as.data.table(unlist(rowRanges(dds)))[seqnames %in% my.ucsc.chr,.(.N,width=mean(width)),Name], # n=2588
                        by.x="ID",
                        by.y="Name")

    # 2. known pre-miRNA
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.pre.miRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    dt.pre.miRNA.fpkm<-merge(dt.mean[,.(`type`="precursor miRNA",`ID`=mirbase_id,meanCount,meanFpm,meanFpkm)],
                        as.data.table(unlist(rowRanges(dds)))[seqnames %in% my.ucsc.chr,.(.N,width=mean(width)),Name],
                        by.x="ID",
                        by.y="Name")

    # 4 known piRNA (piRBase v1.0)
    #load("~/results/RNA-Seq/Boy.Girl.FG.JD.piRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.piRNA.v1-miR.overlap30.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    dt.piRNA.v1.fpkm<-merge(dt.mean[,.(`type`="piRNA",`ID`=pirbase_id,meanCount,meanFpm,meanFpkm)],
                        as.data.table(unlist(rowRanges(dds)))[seqnames %in% my.ucsc.chr,.(.N,width=mean(width)),name],
                        by.x="ID",
                        by.y="name")

    dt.small.RNA.fpkm<-rbind(
                        dt.mature.miRNA.fpkm,
                        dt.piRNA.v1.fpkm,
                        dt.pre.miRNA.fpkm
                        #dt.novel.smallRNA.fpkm,
                        #dt.novel.miRNA.fpkm
                        )
    dt.small.RNA.fpkm[,.(median(width),.N),type]

    # boxplot
    p.known.small.size<-ggplot(dt.small.RNA.fpkm, aes(type,width)) + 
        geom_boxplot(outlier.shape=NA, fill="grey60",width=.5) + 
        #ylim(0,120) +
        coord_cartesian(ylim=c(0,120)) +
        ylab("Size (bp)") + 
        xlab("Types of known small-RNA") +
        theme_Publication() + 
        theme(legend.position="",axis.text.x=element_text(angle=20, hjust=1))
    print(p.known.small.size)

    p.size3<-ggplot(dt.small.RNA.fpkm, aes(width,meanFpkm/10^6)) + 
        geom_point(aes(col=type),size=2,alpha=0.8) +
        geom_segment(data=dt.small.RNA.fpkm,aes(x=width,xend=width,y=0,yend=meanFpkm/10^6,col=type),size=1.2) +
        ggsci::scale_color_jco() +
        xlab("Size (bp) of small RNA") + 
        ylab("Relative Abundance") +
        facet_wrap(~type,nrow=5) +
        theme_Publication() +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank()
        )
    print(p.size3)
}

if(TRUE){
    # novel small-RNA
    load("~/results/RNA-Seq/Placentome/RData/novel.smallRNA.v2.RData") # load 'dt.novel.smallRNA'
    gr.novel.smallRNA<-makeGRangesFromDataFrame(data.frame(dt.novel.smallRNA[seqnames!="MT"]),keep.extra.columns=T) # n=29358 (299 from 'MT' NOT INCLUDED)
    genome(seqinfo(gr.novel.smallRNA))<-"GRCh38"

    # novel miRNA
    load("~/results/RNA-Seq/Boy.Girl.FG.JD.novel.miRNA.2.3.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
    mat.count<-counts(dds)
    novel.miRNA.freq<-apply(mat.count, 1, function(i) sum(i>0)/nrow(dt.samples))
    novel.miRNA.freq[1:3]

    dt.novel.miRNA<-data.table(
                               seqnames=sapply(strsplit(names(novel.miRNA.freq), ":"), function(i) i[1]),
                               start=sapply(strsplit(names(novel.miRNA.freq), ":"), function(i)  as.numeric(i[2])),
                               end=sapply(strsplit(names(novel.miRNA.freq), ":"), function(i)  as.numeric(i[3])),
                               width=sapply(strsplit(names(novel.miRNA.freq), ":"), function(i)  as.numeric(i[3])-as.numeric(i[2])),
                               strand=sapply(strsplit(names(novel.miRNA.freq), ":"), function(i)  i[4]),
                               freq=unname(novel.miRNA.freq)
    )

    i<-0; dl.small.size<-list()
    for(MIN.FREQ in seq(0,1,by=0.3)){
        i<-i+1
        dl.small.size[[i]]=cbind(`Type`="Novel small-RNA",`Sample Frequency`=paste0(">", MIN.FREQ*100,"%"), data.table(data.frame(reduce(gr.novel.smallRNA[gr.novel.smallRNA$freq>MIN.FREQ,]))))

        i<-i+1
        dl.small.size[[i]]=
            cbind(`Type`="Novel miRNA precursor",`Sample Frequency`=paste0(">", MIN.FREQ*100,"%"), dt.novel.miRNA[freq>=MIN.FREQ,-"freq"])
        #    cbind(`Type`="Novel precursor miRNA",`Sample Frequency`=paste0(">", MIN.FREQ*100,"%"), data.table(data.frame(reduce(gr.novel.miRNA[gr.novel.miRNA$freq>=MIN.FREQ,]))))
    }

    rbindlist(dl.small.size)[,.(.N,median(width)),.(Type,`Sample Frequency`)]

    p.novel.small.size<-ggplot(rbindlist(dl.small.size), aes(`Sample Frequency`, width)) +
        geom_boxplot(aes(fill=Type),outlier.shape=NA,width=.7) +
        scale_fill_manual(values=c(cbPalette2[6],cbPalette2[3])) +
        coord_cartesian(ylim=c(0,120)) +
        ylab("Size (bp)") +
        theme_Publication() +
        theme(legend.position="top", axis.text.x=element_text(angle=20, hjust=1))
    print(p.novel.small.size)

}

#################################################
# see bin/count.unmapped.reads.from.mirdeep2.sh #
# see bin/count.mapped.reads.from.mirdeep2.sh #
#################################################
my.RData="~/results/RNA-Seq/Placentome/RData/dt.mapped.smallRNA.RData" # 71MB
if(file.exists(my.RData)){
    load(my.RData)
}else{
    dt.small.rna<-fread("~/results/RNA-Seq/Boy.Girl.FG.JD.mature.miRNA.GRCh38/Meta/meta.ALL.BR.csv")[,-"HtseqFile"] # n=328
    bad.samples=c("91","80C","84C","88") # decidual contamination
    MIN.SCORE=10 # at least this coverage for each reads 

    #########################
    ## 1. Mapped small-RNA ##
    #########################
    # see bin/count.mapped.reads.from.mirdeep2.sh
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.smallRNA<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V4>=MIN.SCORE,.(`size.read`=mean(V2),`size.mapped`=mean(V3),`depth`=mean(V4)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})

    # Q: size of mapped region same as size of read?
    # A: yes 
    sapply(dl.mapped.smallRNA, function(i) nrow(i[size.read!=size.mapped]))

    dt.mapped.smallRNA<-rbindlist(dl.mapped.smallRNA)

    ##############################
    ## 2. Reads Mapped to miRNAs #
    ##############################
    # see bin/count.mapped.reads.from.mirdeep2.sh
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.miRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.miRNA<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.miRNA<-rbindlist(dl.mapped.miRNA) # n=1,406,130

    ##############################
    ## 3. Reads Mapped to piRNAs #
    ##############################
    # see bin/count.mapped.reads.from.mirdeep2.sh
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.piRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.piRNA<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.piRNA<-rbindlist(dl.mapped.piRNA) # n=395,599

    #############################
    ## 4. Reads Mapped to tRNAs #
    #############################
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.tRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.tRNA<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.tRNA<-rbindlist(dl.mapped.tRNA) # n=86,041

    ###############################
    ## 5. Reads Mapped to sncRNAs #
    ###############################
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.sncRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.sncRNA<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.sncRNA<-rbindlist(dl.mapped.sncRNA) # n=1,117,540

    ###############################################
    ## 6-a. Reads Mapped to other exonic (Ens 82) #
    ## -miRNA, -piRNA
    ###############################################
    # see bin/count.mapped.reads.from.mirdeep2.sh
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.exon.82/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.exon<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.exon<-rbindlist(dl.mapped.exon) # n=1,588,569

    ###############################################
    ## 6-b. Reads Mapped to other exonic (Ens 82) #
    ## -miRNA, -piRNA, -tRNA, -sncRNA
    ###############################################
    # see bin/count.mapped.reads.from.mirdeep2.sh
    my.cmd.mapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountMappedRead.exon.82.rev1/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.mapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.mapped.exon.rev1<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328-4
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})
    dt.mapped.exon.rev1<-rbindlist(dl.mapped.exon.rev1) # n=448,388

    ################################### 
    ## reads after mapping with followings 
    ## 1. mirBase (v21) mature miR only
    ################################### 
    ## remaining reads
    ## see bin/count.unmapped.reads.from.mirdeep2.sh
    my.cmd.unmapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountUnmappedRead.miRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.unmapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.unmapped.smallRNA.miR<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328 
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})

    #dl.unmapped.smallRNA.miR[[1]]
    dt.unmapped.smallRNA.miR<-rbindlist(dl.unmapped.smallRNA.miR)

    ################################### 
    ## reads after mapping with followings 
    ## 1. mirBase (v21) mature miR only
    ## 2. pirBase (v1)
    ################################### 
    ## remaining reads
    ## see bin/count.unmapped.reads.from.mirdeep2.sh
    my.cmd.unmapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountUnmappedRead.miRNA.piRNA/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )
    foo.v2<-sapply(my.cmd.unmapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v2<-data.table(t(simplify2array(strsplit(foo.v2, " ")))) # n=330

    dl.unmapped.smallRNA.miR.piR<-apply(merge(dt.fg.jd.v2, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328 
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})

    #dl.unmapped.smallRNA.miR.piR[[1]]
    dt.unmapped.smallRNA.miR.piR<-rbindlist(dl.unmapped.smallRNA.miR.piR)

    ################################### 
    ## reads unmapped with following ##
    ## 1. mirBase (v22) 
    ## 2. novel mirBase
    ## 3. pirBase (v2)
    ## 4. exonic region of GRCh38.94
    ################################### 
    ## remaining reads
    ## see bin/count.unmapped.reads.from.mirdeep2.sh
    my.cmd.unmapped<-c(
            "for i in `ls /home/ssg29/results/RNA-Seq/Placentome/smallRNA/CountUnmappedRead/*csv.gz`; do SLX=`echo $i | cut -d'/' -f 9 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f 9 | cut -d'.' -f4`; printf \"$SLX $BARCODE $i\n\" ; done"
            )

    foo.v3<-sapply(my.cmd.unmapped, system, intern=T, USE.NAMES=F)
    dt.fg.jd.v3<-data.table(t(simplify2array(strsplit(foo.v3, " ")))) # n=330

    dl.unmapped.smallRNA<-apply(merge(dt.fg.jd.v3, dt.small.rna[!SampleName %in% bad.samples], by.x=c("V1","V2"), by.y=c("Library","BarCode")), # n=328 
                    1, function(i){fread(cmd=paste("zcat",i[3]))[V3>=MIN.SCORE,.(`size.mapped`=mean(V2),`depth`=mean(V3)),V1][,`:=`(SLX=i[1],Barcode=i[2],SampleName=i[4],Source=i[8],CRN=as.integer(i[9]))]})

    dl.unmapped.smallRNA[[1]]
    dt.unmapped.smallRNA<-rbindlist(dl.unmapped.smallRNA)

    ##
    ## Save it
    save(dt.mapped.smallRNA, dt.mapped.miRNA, dt.mapped.piRNA, dt.mapped.tRNA, dt.mapped.sncRNA, dt.mapped.exon, dt.mapped.exon.rev1, 
         dt.unmapped.smallRNA.miR, dt.unmapped.smallRNA.miR.piR, dt.unmapped.smallRNA, file=my.RData)
}


if(TRUE){
    dt.small.mapped2<-rbind(
        dt.mapped.smallRNA[,.(`Read type`="Mapped reads",V1,size.mapped,depth)],
        dt.mapped.miRNA[,.(`Read type`="+miRNA",V1,size.mapped,depth)],
        dt.mapped.piRNA[,.(`Read type`="+piRNA",V1,size.mapped,depth)],
        dt.unmapped.smallRNA.miR[,.(`Read type`="-miRNA",V1,size.mapped,depth)],
        dt.unmapped.smallRNA.miR.piR[,.(`Read type`="-miRNA-piRNA",V1,size.mapped,depth)],
        dt.unmapped.smallRNA[,.(`Read type`="-miRNA-piRNA-exon",V1,size.mapped,depth)]
        )
    dt.small.mapped2$`Read type`=factor(dt.small.mapped2$`Read type`, levels=c("Mapped reads","+miRNA","+piRNA","-miRNA","-miRNA-piRNA","-miRNA-piRNA-exon"))

    p.abc.v2<-ggplot(dt.small.mapped2[`Read type` %in% c("Mapped reads","-miRNA","-miRNA-piRNA-exon"),.(.N,count=sum(depth)),.(`Read type`,size=round(size.mapped))], aes(size, count)) +
        geom_point(aes(col=`Read type`),size=5,alpha=.8) +
        geom_line(aes(col=`Read type`),size=1,alpha=.8) +
        scale_y_continuous(trans = log10_trans(),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x))) +
        ylab("Read count") + 
        xlab("Size (bp)") +
        #scale_color_manual(values=c(cbPalette2[1],cbPalette2[2],cbPalette2[6],cbPalette2[3])) +
        scale_color_manual(values=c(cbPalette2[1],cbPalette2[2],cbPalette2[3])) +
        theme_Publication() +
        theme(legend.position=c(.65,.85))
    print(p.abc.v2)

    p.abc.v3<-ggplot(dt.small.mapped2[`Read type` %in% c("Mapped","+miRNA","+piRNA","-miRNA-piRNA"),.(.N,count=sum(depth)),.(`Read type`,size=round(size.mapped))], aes(size, count)) +
        geom_point(aes(col=`Read type`),size=5,alpha=.85) +
        geom_line(aes(col=`Read type`),size=1,alpha=.85) +
        scale_y_continuous(trans = log10_trans(),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x))) +
        ylab("Read count") + 
        xlab("Size (bp)") +
        scale_color_manual(values=c(cbPalette2[1],cbPalette2[4],cbPalette2[5],cbPalette2[6])) +
        theme_Publication() +
        theme(legend.position=c(.65,.85))
    print(p.abc.v3)

    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/smallRNA/mapped.unmapped.read.size.tiff")
    tiff(filename=file.name,width=15, height=7,units="in",res=300, compression = 'lzw') #A4 size (2400px x 1500px)
    cowplot::plot_grid(p.abc.v2+coord_cartesian(ylim=c(10^2,10^9)), p.abc.v3+coord_cartesian(ylim=c(10^2,10^9)), nrow=1, labels="auto",label_size=20)
    dev.off()

    p.known.novel<-cowplot::plot_grid(p.known.small.size, p.novel.small.size, labels=c("b","c"), label_size=27, ncol=1)
    file.name<-file.path("~/results/RNA-Seq/Placentome/Figures/smallRNA/mapped.unmapped.smallRNA.read.size.tiff")
    tiff(filename=file.name,width=10, height=8,units="in",res=300, compression = 'lzw')
    cowplot::plot_grid(p.abc.v2, p.known.novel, labels=c("a",""), label_size=27, axis="bt",ncol=2, rel_widths=c(1,1))
    dev.off()

}

