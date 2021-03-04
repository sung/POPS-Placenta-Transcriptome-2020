minFpkm=0.1
TR_PREFIX='GRCh38' 
ENS_VER=82 # re-construction is based on GRCh38.82 which is used for StringTie
mySource="Placentome" 
myCohort="POPS" 
myProject=paste(mySource,myCohort,sep=".") 

source("lib/cufflink.R")
source("lib/local.R")
loadData() # dt.frag, dt.novel.miRNA, dt.pops.ciri
library(scales)

		cat("making class_code ...\n")
		dummy<-sapply(names(rev(sampleFreq)), function(i) table(dt.frag[fpkm.iqr >= minFpkm & evi.ratio >=sampleFreq[i],class_code])) # isa 'matrix'
		dummy<-reshape2::melt(dummy) # isa 'data.frame'
		colnames(dummy)<-c("Code","freq","cnt")
		#   Var1  Var2 value
		#1     = >=90% 31503
		#2     . >=90%     0
		#3     c >=90%   301
		#4     e >=90%     0
		#5     i >=90%   751
		dummy<-data.table(merge(dummy, dt.cuffClass, by.x="Code",by.y="class_code"))
		dummy2<-dummy[,list(cnt=sum(cnt)),list(freq,class)][order(freq,class)]
        dummy2$freq<-factor(dummy2$freq, levels=rev(names(sampleFreq))) #character to factor  
        
		p.1 <-ggplot(dummy2, aes(x=freq, y=cnt, fill=class)) + 
			geom_bar(position="dodge", stat="identity", width=0.5) +  
			scale_fill_manual(name="Type",values=cbPalette) +
            xlab("Sample frequency") + ylab("Number") +
            scale_y_continuous(trans = log10_trans(),
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))) +
			theme_Publication()	+
			theme(axis.text.x=element_text(angle=45, hjust=1),legend.position=c(.76,.8),
                  legend.title=element_text(size=rel(.8)),
                  legend.text=element_text(size=rel(.8)),
                  legend.key.size= unit(0.5, "cm"),
                  legend.background=element_rect(fill=alpha("white",0.6)))
		print(p.1)


		#######################
		## Novel Transcripts ##
		## by Sample Freq    ##
		## j, i, u only      ##
		#######################
		cat("counting novel transcript...\n")
		dt.frag[,class_code:=as.factor(class_code)]
		dummy<-sapply(seq(0, 1, by=.1), function(i) table(dt.frag[ fpkm.iqr>=minFpkm & evi.ratio>=i,class_code]) )
		#paste0(">=",seq(0, 1, by=.1)*100,"%") # isa matrix

		colnames(dummy)<-c(">0%",paste0(">=",seq(0.1, 0.9, by=.1)*100,"%"),"100%") # isa matrix
		dummy<-reshape2::melt(dummy)
		colnames(dummy)<-c("Code","freq","cnt")
		dummy<-data.table(merge(dummy, dt.cuffClass, by.x="Code",by.y="class_code"))
		dt.novel<-dummy[Code %in% c("j","i","u")]
		dt.novel$Code<-droplevels(dt.novel$Code)
		dt.novel$desc<-factor(dt.novel$desc, levels=dt.cuffClass[class_code %in% c("j","i","u"),desc]) #character to factor

		# all bar-chart (log10 y-axis)
		p.2<-ggplot(dt.novel[cnt!=0], aes(x=freq, y=cnt, fill=desc)) +  
			geom_bar(stat="identity", position="dodge") +
			scale_fill_brewer(name="Class", palette="Set2") +
            xlab("Sample frequency") + ylab("Number") +
            scale_y_continuous(trans = log10_trans(),
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1),legend.position=c(.7,.8),
                  legend.title=element_text(size=rel(.8)),
                  legend.text=element_text(size=rel(.8)),
                  legend.key.size= unit(0.5, "cm"),
                  legend.background=element_rect(fill=alpha("white",0.6)))
		print(p.2)

        ###############
        # Novel miRNA #
        # n=5819 (from gr.novel.miRNA or dt.novel.miRNA)
        # n=4351 (from reduced(gr.novel.miRNA) - same as novel.miRNA.GRCh38.merged.bed)
        ###############
        library(DESeq2)
        load("RData/novel.miRNA.2.3.GRCh38.deseq.ALL.RData")
        mat.count<-counts(dds)
        dim(mat.count)
        novel.miRNA.freq<-apply(mat.count, 1, function(i) sum(i>0)/nrow(dt.samples))
        novel.miRNA.freq[1:3]

        i<-0; dl.dummy.mi<-list()
        for(MIN.FREQ in seq(0,1,by=0.1)){
            i<-i+1
            ##cnt=dt.novel.miRNA[freq>=MIN.FREQ,.N] # use below
            #cnt=length(reduce(gr.novel.miRNA[gr.novel.miRNA$freq>=MIN.FREQ,]))
            cnt=sum(novel.miRNA.freq>=MIN.FREQ)
            dl.dummy.mi[[i]]<-data.table(`min freq`=MIN.FREQ,cnt)
        }

        ####################
        ## Novel smallRNA ##
        # n=18611 novel smallRNA assembled from 293 samples
        ####################
        library('GenomicRanges')
        load("RData/novel.smallRNA.v2.RData") # load  dt.novel.smallRNA (n=29657) 
        gr.novel.smallRNA<-GenomicRanges::makeGRangesFromDataFrame(data.frame(dt.novel.smallRNA[seqnames!="MT"]),keep.extra.columns=T) # n=29358
        genome(seqinfo(gr.novel.smallRNA))<-"GRCh38"
        
        i<-0; dl.dummy.smallRNA<-list()
        for(MIN.FREQ in seq(0,1,by=0.1)){
            i<-i+1
            #cnt=length(reduce(gr.novel.smallRNA[gr.novel.smallRNA$freq>=MIN.FREQ,]))
            dl.dummy.smallRNA[[i]]<-data.table(`min freq`=MIN.FREQ,cnt)
        }

        dt.dummy.small<-rbind(
            cbind(`RNA-type`="novel miRNA",rbindlist(dl.dummy.mi)),
            cbind(`RNA-type`="novel small-RNA",rbindlist(dl.dummy.smallRNA))
            #cbind(`RNA-type`="circRNA",dt.dummy.circ[`min junction count`==2,.(`min freq`,cnt)])
            )
        #dt.dummy.small$`RNA-type`<-factor(dt.dummy.small$`RNA-type`, levels=c("circRNA","novel small-RNA","novel miRNA")) #character to factor  
        dt.dummy.small$`RNA-type`<-factor(dt.dummy.small$`RNA-type`, levels=c("novel miRNA","novel small-RNA")) #character to factor  

        p.3<-ggplot(dt.dummy.small, aes(`min freq`,cnt)) + 
            geom_bar(aes(fill=`RNA-type`),position="dodge",stat="identity",size=5,alpha=.7) + 
			scale_fill_brewer(palette="Set1") +
            #ggsci::scale_fill_d3() +
            scale_y_continuous(name="Number",
                            trans = log10_trans(),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x))) +
            #scale_y_log10(name="Number of novel miRNA",
            #              breaks=c(0,10,10^2,10^3,10^4),
            #              #limits=c(0,10^4),
            #              labels=trans_format("log10",math_format(10^.x))) +
            scale_x_continuous(name="Sample frequency",breaks=seq(0,1,by=0.1),labels=c(">0%",paste0(">=",seq(10,90,by=10),"%"),"100%")) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle=45,hjust=1),legend.position=c(.8,.8), 
                  legend.title=element_text(size=rel(.8)),
                  legend.text=element_text(size=rel(.8)),
                  legend.key.size= unit(0.5, "cm"),
                  legend.background=element_rect(fill=alpha("white",0.6)))
        print(p.3)

        ############
        ## circRNA #
        ############
        load("RData/dl.circRNA.polyA.RData")
        dt.circRNA.PA<-dl.circRNA.polyA[["FGmRNA2014"]]

        i<-0; dl.dummy<-list()
        for(MIN.JNC.CNT in seq(2,10)){ # 2~10
            for(MIN.FREQ in seq(0,1,by=0.1)){
                i<-i+1
                #cnt<-dt.pops.ciri[`#junction_reads`>=MIN.JNC.CNT,.(`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops)),circRNA_ID][freq>=MIN.FREQ,.N]
                cnt<-dt.pops.ciri[`#junction_reads`>=MIN.JNC.CNT & !circRNA_ID %in% dt.circRNA.PA$circRNA_ID,.(`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops)),circRNA_ID][freq>=MIN.FREQ,.N]
                dl.dummy[[i]]<-data.table(`min junction count`=MIN.JNC.CNT,`min freq`=MIN.FREQ,cnt)
            }
        }

        #write.csv(rbindlist(dl.dummy), file="~/results/RNA-Seq/Placentome/Figures/TranscriptNumber/circRNA.number.by.sample.freq.csv", quote=F,row.names=F)
        dt.dummy.circ<-rbindlist(dl.dummy)
		dt.dummy.circ<-dt.dummy.circ[,`min junction count`:=as.factor(`min junction count`)]

        dt.dummy.circ[`min junction count`==2]

        p.4.a<-ggplot(dt.dummy.circ[`min junction count` %in% c(2,5,10)], aes(`min freq`,cnt)) + 
            geom_point(stat="identity",size=4,alpha=.7) + 
            geom_line(aes(col=`min junction count`),size=1.2,alpha=.7,linetype="longdash") +
            ggsci::scale_color_d3(name="Back-spliced junction (≥)") +
            scale_y_continuous(name="Number",trans = log10_trans(),
                                breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x))) +
            scale_x_continuous(name="Sample frequency",breaks=seq(0,1,by=0.1),labels=c(">0%",paste0(">=",seq(10,90,by=10),"%"),"100%")) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle=45,hjust=1),legend.position=c(.74,.8),  
                  legend.title=element_text(size=rel(.8)),
                  legend.text=element_text(size=rel(.8)),
                  legend.key.size= unit(0.5, "cm"),
                  legend.background=element_rect(fill=alpha("white",0.6)))
        print(p.4.a)

        file.name<-file.path("Figures/transcript.number.cuffcompare.a.b.c.d.tiff")
        tiff(file.name, width=12, height=10,units="in",res=300, compression = 'lzw') #A4 size
        cowplot::plot_grid(p.1,p.2,p.3,p.4.a, labels="auto", label_size=20)
        dev.off()

#############################
## Assessment of assembler ##
#############################
if(TRUE){
    ## with -R option to adjust sensitivity
    # splicing-pattern ~= Transcript level?
    # splice-junction ~= Intron-chain?
    dt.foo.R<-fread("data/assessment.cuffcompare.gffcompare.taco.str-merge.csv")
    dt.bar.R<-fread("data/assessment.cuffcompare.by.FPKM.sample.Freq.csv")
    ## with -Q option to adjust precision 
    # splicing-pattern ~= Transcript level?
    # splice-junction ~= Intron-chain?
    dt.foo.Q<-fread("data/assessment.cuffcompare.gffcompare.taco.str-merge.Q.csv")
    dt.bar.Q<-fread("data/assessment.cuffcompare.by.FPKM.sample.Freq.Q.csv")
	dt.bar.Q[,`minimum FPKM`:=as.factor(`minimum FPKM`)]

    ## Both R+Q
    dt.foo<-rbind(
                  cbind(`Adjust`="Adjusted Sensitivity", dt.foo.R),
                  cbind(`Adjust`="Adjusted Precision", dt.foo.Q))
    dt.bar<-rbind(
                  cbind(`Adjust`="Adjusted Sensitivity", dt.bar.R),
                  cbind(`Adjust`="Adjusted Precision", dt.bar.Q))
    setnames(dt.bar, "minimum FPKM","FPKM (≥)")

    p5<-ggplot(dt.bar[Level %in% c("Base","Intron-chain","Transcript")], aes(Sensitivity,Precision)) + 
        geom_point(aes(col=`FPKM (≥)`),size=5,alpha=0.7) +
        geom_point(data=dt.foo[Level %in% c("Base","Intron-chain","Transcript")], aes(shape=`Meta-Assembler`),size=5,alpha=0.7) +
        geom_line(aes(linetype=`Sample Frequency`),size=0.8) +
        ggsci::scale_color_d3() +
        facet_grid(Adjust~Level) +
        theme_Publication()

	my.file.name<- file.path("Figures/tr.number.POPS.cuffcompare.gffcompare.TACO.STR-merge.tiff") # per transcript (TCON_)
	tiff(filename=my.file.name,width=12,height=15,units="in",res=300, compression = 'lzw')
    p.1234<-cowplot::plot_grid(p.1,p.2,p.3,p.4.a, nrow=2,labels="auto",label_size=27)
    cowplot::plot_grid(p.1234,p5,labels=c("","e"),nrow=2,label_size=27,rel_heights=c(1.2,1))
    dev.off()


}


