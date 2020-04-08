source("lib/local.R")
load("RData/dl.manhattan.RData") # loads dl.manhattan (for miRNA, piRNA and circRNA)

lapply(dl.manhattan,nrow)

#################
## A. chr-wide ##
#################
dt.manhattan<-unique(rbindlist(dl.manhattan))
dt.manhattan[,`score`:=ifelse(Type=="circRNA",meanCount/10^3,meanFpkm/10^6)]
gr.manhattan<-GenomicRanges::makeGRangesFromDataFrame(dt.manhattan,keep.extra.columns=T)
seqinfo(gr.manhattan)
gr.manhattan <- keepSeqlevels(gr.manhattan,my.ucsc.chr)
seqinfo(gr.manhattan)
seqlengths(gr.manhattan)
GenomeInfoDb::seqlevelsStyle(gr.manhattan) = "Ensembl"  # e.g. chrX => X
gr.manhattan$Type<-factor(gr.manhattan$Type,levels=c("circRNA","miRNA","piRNA"))

p.man<-ggbio::plotGrandLinear(gr.manhattan, facets=Type ~., aes(y=score), 
                              color=c("gray1", "gray50"), size=1.8, alpha=.9) + 
    xlab("Chromosomes") + 
    theme_Publication() + 
    theme(axis.text.x = element_text(size=rel(.7)), legend.position="")
class(p.man) # isa 'ggbio'

# Let's hack the code of ggbio
foo<-ggplot_build(p.man@ggplot) # isa 'list'
bar<-cbind(data.table(foo$data[[1]]), ID=dt.manhattan$ID, RANK=dt.manhattan$RANK, INDEX_PCT=dt.manhattan$INDEX_PCT)
bar[,Type:=ifelse(PANEL==1,"circRNA",ifelse(PANEL==2,"miRNA","piRNA"))]; bar$Type<-factor(bar$Type)
# Let's put things together as a ggplot object
p.all<-p.man@ggplot + facet_wrap(~Type,nrow=3,scales="free_y",strip.position="right") +
    geom_point(data=bar[Type %in% c("circRNA","miRNA") & INDEX_PCT<=1], aes(x,y),col="orange2",size=2) +
    geom_point(data=bar[Type=="piRNA" & INDEX_PCT<=.1], aes(x,y),col="orange2",size=2) +
    ggrepel::geom_text_repel(data=bar[Type=="circRNA" & RANK<=3], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="miRNA" & RANK<=4], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="piRNA" & RANK==1], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) 

# NB. this version of ggbio (ggbio_1.26.1) does not work with 'auptoplot' if ggplot > 2.2.1 used
# hence, ggplot2 downgraded to 2.2.1 to use ::autoplot
#devtools::install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org")
#####################################
## B. chr19 only (miRNA & circRNA) ##
#####################################
gr.manhattan.chr19<-gr.manhattan[seqnames(gr.manhattan)=="19" & gr.manhattan$Type %in% c("miRNA"),]
seqlevels(gr.manhattan.chr19)<-"19"
gr.manhattan.chr19$Type<-droplevels(gr.manhattan.chr19$Type)

p.man.chr19<-ggbio::autoplot(gr.manhattan.chr19, geom="point", aes(y=score), color="gray1",size=2.3, alpha=.9, spaceline = T) + 
            xlab("Chromosome 19") + 
            facet_grid(Type~.) + 
            theme_Publication()
class(p.man.chr19) # isa 'ggbio'

# Let's hack the code of ggbio
foo<-ggplot_build(p.man.chr19@ggplot)
bar<-cbind(data.table(foo$data[[1]]), ID=gr.manhattan.chr19$ID, RANK=gr.manhattan.chr19$RANK, INDEX_PCT=gr.manhattan.chr19$INDEX_PCT)
# text annotation for RANK <=10
p.chr19<-p.man.chr19@ggplot + facet_wrap(~Type,nrow=2,scales="free_y",strip.position="right") +
    geom_point(data=bar[INDEX_PCT<=1], aes(x,y),col="orange2",size=2.5) +
    ggrepel::geom_text_repel(data=bar[RANK<=6], aes(x,y,label=ID), direction="y", nudge_x=-10^8/7, segment.size=0.2, size=3.7, segment.color = "grey50") +
    ggrepel::geom_text_repel(data=bar[RANK>7 & RANK<=10], aes(x,y,label=ID), direction="y", nudge_y=-.5,nudge_x=-10^8/4, segment.size=0.2, size=3.7, segment.color = "grey50") 

##########################
## C. chrM only (piRNA) ##
##########################
gr.manhattan.chrMT<-gr.manhattan[seqnames(gr.manhattan)=="MT" & gr.manhattan$Type=="piRNA",]
seqlevels(gr.manhattan.chrMT)<-"MT"
gr.manhattan.chrMT$Type<-droplevels(gr.manhattan.chrMT$Type)
p.man.chrMT<-ggbio::autoplot(gr.manhattan.chrMT, geom="point", aes(y=score), color="gray1",size=2.5, alpha=.9, spaceline = T) + ylim(0,1.5) + 
            xlab("Mitochondrial Chromosome") +
            facet_grid(Type~.) + 
            theme_Publication()

foo<-ggplot_build(p.man.chrMT@ggplot)
bar<-cbind(data.table(foo$data[[1]]), ID=gr.manhattan.chrMT$ID, RANK=gr.manhattan.chrMT$RANK, INDEX_PCT=gr.manhattan.chrMT$INDEX_PCT)
p.chrMT<-p.man.chrMT@ggplot + 
    geom_point(data=bar[INDEX_PCT<=.1], aes(x,y),col="orange2",size=3) +
    ggrepel::geom_text_repel(data=bar[RANK<=10], aes(x,y,label=ID), direction="y", nudge_x=-10^4/3, hjust=1, segment.size=0.2, size=3.7, segment.color = "grey50")

###########################
## D. chr14 only (piRNA) ##
###########################
gr.manhattan.chr14<-gr.manhattan[seqnames(gr.manhattan)=="14" & gr.manhattan$Type=="piRNA",]
seqlevels(gr.manhattan.chr14)<-"14"
gr.manhattan.chr14$Type<-droplevels(gr.manhattan.chr14$Type)
p.man.chr14<-ggbio::autoplot(gr.manhattan.chr14, geom="point", aes(y=score), color="gray45",size=2.5, alpha=.95, spaceline = T) + 
            xlab("Chromosome 14") +
            facet_grid(Type~.) + 
            theme_Publication()

foo<-ggplot_build(p.man.chr14@ggplot)
bar<-cbind(data.table(foo$data[[1]]), ID=gr.manhattan.chr14$ID, RANK=gr.manhattan.chr14$RANK, INDEX_PCT=gr.manhattan.chr14$INDEX_PCT)
p.chr14<-p.man.chr14@ggplot + 
    geom_point(data=bar[INDEX_PCT<=.1], aes(x,y),col="orange2",size=3) +
    scale_y_continuous(breaks=c(0,0.5,1)) +
    ggrepel::geom_text_repel(data=bar[RANK<=10], aes(x,y,label=ID), direction="y", hjust=1, nudge_x=-10^7,segment.size=0.2, segment.color = "grey50", size=3.7)

##########################
## E. chr1 only (piRNA) ##
##########################
gr.manhattan.chr1<-gr.manhattan[seqnames(gr.manhattan)=="1" & gr.manhattan$Type=="piRNA",]
seqlevels(gr.manhattan.chr1)<-"1"
gr.manhattan.chr1$Type<-droplevels(gr.manhattan.chr1$Type)

p.man.chr1<-ggbio::autoplot(gr.manhattan.chr1, geom="point", aes(y=score), color="gray1",size=2.5, alpha=.9, spaceline = T) + 
            xlab("Chromosome 1") +
            facet_grid(Type~.) + 
            theme_Publication()

foo<-ggplot_build(p.man.chr1@ggplot)
bar<-cbind(data.table(foo$data[[1]]), ID=gr.manhattan.chr1$ID, RANK=gr.manhattan.chr1$RANK, INDEX_PCT=gr.manhattan.chr1$INDEX_PCT)
p.chr1<-p.man.chr1@ggplot + 
    geom_point(data=bar[INDEX_PCT<=.1], aes(x,y),col="orange2",size=3) +
    ggrepel::geom_text_repel(data=bar[order(RANK)][1], aes(x,y,label=ID), direction="y", hjust=0, nudge_x=-10^8/2,segment.size=0.2, segment.color = "grey50", size=3.5) +
    ggrepel::geom_text_repel(data=bar[order(RANK)][2:6], aes(x,y,label=ID), direction="y", hjust=0, nudge_x=10^8/2,segment.size=0.2, segment.color = "grey50", size=3.5) +
    ggrepel::geom_text_repel(data=bar[order(RANK)][7], aes(x,y,label=ID), direction="y", hjust=0, nudge_x=-10^8/3,segment.size=0.2, segment.color = "grey50", size=3.5) 
print(p.chr1)

###########################
# combine A,B,C,D,E above #
###########################
file.name<-file.path("Figures/Fig2.manhattan.tiff")
tiff(filename=file.name,width=10, height=14,units="in",res=300, compression = 'lzw')
cowplot::plot_grid(p.all, p.chr19, p.chrMT, p.chr14, p.chr1, 
                   labels="auto",label_size=20,ncol=1, align="v", 
                   rel_heights=c(2.5,1,1,1,1))
dev.off()
