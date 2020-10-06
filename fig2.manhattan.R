source("lib/local.R")
load("RData/dl.manhattan.RData") # loads dl.manhattan (for miRNA, piRNA, circRNA and sncRNA)

## NB, no chrM (or MT) included
my.ucsc.chr<-paste0("chr",c(seq(1:22),"X","Y","M")) # 1,2 ... 22, X, Y
lapply(dl.manhattan,nrow)

#################
## A. chr-wide ##
#################
dt.manhattan<-unique(rbindlist(dl.manhattan))
dt.manhattan[,.N,Type]
dt.manhattan[,`score`:=ifelse(Type=="circRNA",meanCount/10^3,meanFpkm/10^6)]
gr.manhattan<-GenomicRanges::makeGRangesFromDataFrame(dt.manhattan,keep.extra.columns=T)
seqinfo(gr.manhattan)
gr.manhattan <- keepSeqlevels(gr.manhattan,my.ucsc.chr)
seqinfo(gr.manhattan)
seqlengths(gr.manhattan)
GenomeInfoDb::seqlevelsStyle(gr.manhattan) = "Ensembl"  # e.g. chrX => X
gr.manhattan$Type<-factor(gr.manhattan$Type,levels=names(dl.manhattan))

p.man<-ggbio::plotGrandLinear(gr.manhattan, facets=Type ~., aes(y=score), scales="free_y" ,
                              color=c("gray1", "gray50"), size=1.8, alpha=.9) + 
    xlab("Chromosomes") + 
    theme_Publication() + 
    theme(axis.text.x = element_text(size=rel(.7)), legend.position="")
class(p.man) # isa 'ggbio'

# Let's hack the code of ggbio
foo<-ggplot_build(p.man@ggplot) # isa 'list'
bar<-cbind(data.table(foo$data[[1]]), ID=dt.manhattan$ID, RANK=dt.manhattan$RANK, INDEX_PCT=dt.manhattan$INDEX_PCT)
bar[,Type:=ifelse(PANEL==1,"circRNA",ifelse(PANEL==2,"miRNA",ifelse(PANEL==3,"piRNA","sncRNA")))]; bar$Type<-factor(bar$Type)
# Let's put things together as a ggplot object
p.all<-p.man@ggplot + 
    #facet_wrap(~Type,nrow=length(dl.manhattan),scales="free_y",strip.position="right") + # NOT NEEDED
    geom_point(data=bar[Type %in% c("circRNA","miRNA","sncRNA") & INDEX_PCT<=1], aes(x,y),col="orange2",size=2) +
    geom_point(data=bar[Type=="piRNA" & INDEX_PCT<=.1], aes(x,y),col="orange2",size=2) +
    ggrepel::geom_text_repel(data=bar[Type=="circRNA" & RANK<=3], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="miRNA" & RANK<=4], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="piRNA" & RANK==1], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="sncRNA" & RANK<=3], aes(x,y,label=ID), hjust=0, segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar[Type=="sncRNA" & group==11][1:3], aes(x,y,label=ID), hjust=0,nudge_x=-10^6,segment.size=0.2, segment.color = "grey50", size=3.7) +
    ylab("") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

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
p.chr19<-p.man.chr19@ggplot + 
    #facet_wrap(~Type,nrow=2,scales="free_y",strip.position="right") + # NOT NEEDED
    geom_point(data=bar[INDEX_PCT<=1], aes(x,y),col="orange2",size=2.5) +
    ggrepel::geom_text_repel(data=bar[RANK<=6], aes(x,y,label=ID), direction="y", nudge_x=-10^8/7, segment.size=0.2, size=3.7, segment.color = "grey50") +
    ggrepel::geom_text_repel(data=bar[RANK>7 & RANK<=10], aes(x,y,label=ID), direction="y", nudge_y=-.5,nudge_x=-10^8/4, segment.size=0.2, size=3.7, segment.color = "grey50") +
    ylab("") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

###########
## C19MC ##
###########
gr.manhattan.c19mc<-gr.manhattan[seqnames(gr.manhattan)=="19" & gr.manhattan$Type %in% c("miRNA") & start(gr.manhattan)>5e7 & end(gr.manhattan)<5.6e7]
seqlevels(gr.manhattan.c19mc)<-"19"
gr.manhattan.c19mc$Type<-droplevels(gr.manhattan.c19mc$Type)

p.man.c19mc<-ggbio::autoplot(gr.manhattan.c19mc, geom="point", aes(y=score), color="gray1",size=2.3, alpha=.9, spaceline = T) + 
            xlab("Chromosome 19") + 
            facet_grid(Type~.) + 
            theme_Publication()
class(p.man.c19mc) # isa 'ggbio'

# Let's hack the code of ggbio
foo.c19mc<-ggplot_build(p.man.c19mc@ggplot)
bar.c19mc<-cbind(data.table(foo.c19mc$data[[1]]), ID=gr.manhattan.c19mc$ID, RANK=gr.manhattan.c19mc$RANK, INDEX_PCT=gr.manhattan.c19mc$INDEX_PCT)
# text annotation for RANK <=10
p.c19mc<-p.man.c19mc@ggplot + 
    #facet_wrap(~Type,nrow=2,scales="free_y",strip.position="right") + # NOT NEEDED
    geom_point(data=bar.c19mc[INDEX_PCT<=1], aes(x,y),col="orange2",size=2.5) + # n=12
    geom_point(data=bar.c19mc[INDEX_PCT>1 & INDEX_PCT<=5], aes(x,y),col="cornflowerblue",size=2.5, alpha=.8) +
    ggrepel::geom_text_repel(data=bar.c19mc[INDEX_PCT<=1 & RANK<=7,.(x=mean(x),y=mean(y)),ID], aes(x,y,label=ID), direction="y", 
                             nudge_x=bar.c19mc[INDEX_PCT<=1 & RANK<=7,.(x=mean(x),y=mean(y)),ID]$x-5.4e7, hjust=1, 
                             segment.size=0.2, size=3.7, segment.color = "grey50") +
    ggrepel::geom_text_repel(data=bar.c19mc[INDEX_PCT<=1 & RANK>=8 & RANK<=13], aes(x,y,label=ID), direction="y", 
                             nudge_x=5.4e7-bar.c19mc[INDEX_PCT<=1 & RANK>=8 & RANK<=13]$x, hjust=0, 
                             segment.size=0.2, size=3.7, segment.color = "grey50") +
    ggrepel::geom_text_repel(data=bar.c19mc[INDEX_PCT<=1 & RANK==16], aes(x,y,label=ID), direction="y", 
                             nudge_x=5.2e7-bar.c19mc[INDEX_PCT<=1 & RANK==16]$x, hjust=0, 
                             segment.size=0.2, size=3.7, segment.color = "grey50") +
    ggrepel::geom_text_repel(data=bar.c19mc[INDEX_PCT<=1 & RANK==23], aes(x,y,label=ID), direction="y", 
                             nudge_x=bar.c19mc[INDEX_PCT<=1 & RANK==23]$x-5.4e7, hjust=1, 
                             segment.size=0.2, size=3.7, segment.color = "grey50") +
    ylab("") + xlab("C19MC") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))


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
    ggrepel::geom_text_repel(data=bar[RANK<=10], aes(x,y,label=ID), direction="y", nudge_x=-10^4/3, hjust=1, segment.size=0.2, size=3.7, segment.color = "grey50") +
    ylab("") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

####################################
## D-1. chr14 only (piRNA+sncRNA) ##
####################################
gr.manhattan.chr14<-gr.manhattan[seqnames(gr.manhattan)=="14" & gr.manhattan$Type%in%c("piRNA","sncRNA"),]
seqlevels(gr.manhattan.chr14)<-"14"
gr.manhattan.chr14$Type<-droplevels(gr.manhattan.chr14$Type)
p.man.chr14<-ggbio::autoplot(gr.manhattan.chr14, geom="point", aes(y=score), color="gray45",size=2.2, alpha=.95, spaceline = T) + 
            xlab("Chromosome 14") +
            facet_grid(Type~.) + 
            theme_Publication()

foo.14<-ggplot_build(p.man.chr14@ggplot)
bar.chr14<-cbind(data.table(foo.14$data[[1]]), ID=gr.manhattan.chr14$ID, RANK=gr.manhattan.chr14$RANK, INDEX_PCT=gr.manhattan.chr14$INDEX_PCT)
bar.chr14[,Type:=ifelse(PANEL==1,"piRNA","sncRNA")]
p.chr14<-p.man.chr14@ggplot + 
    facet_wrap(~Type,nrow=2,scales="free_y",strip.position="right") +
    geom_point(data=bar.chr14[PANEL==1 & INDEX_PCT<=.1], aes(x,y),col="orange2",size=3) +# piRNA
    geom_point(data=bar.chr14[PANEL==2 & INDEX_PCT<=1], aes(x,y),col="orange2",size=3) + # sncRNA
    scale_y_continuous(breaks=c(0,0.5,1)) +
    ggrepel::geom_text_repel(data=bar.chr14[PANEL==1 & RANK<=10], aes(x,y,label=ID), direction="y", hjust=1, nudge_x=-10^7,segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar.chr14[PANEL==2 & RANK<=15], aes(x,y,label=ID), direction="y", hjust=1, nudge_x=-10^7,segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar.chr14[PANEL==2 & RANK>15 & INDEX_PCT<1], aes(x,y,label=ID), direction="y", hjust=1, nudge_x=-(10^7)*2.2,segment.size=0.2, segment.color = "grey50", size=3.7) +
    ylab("") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

#####################################
## D-2. C14MC (miRNA+piRNA+sncRNA) ##
#####################################
gr.manhattan.c14mc<-gr.manhattan[seqnames(gr.manhattan)=="14" & gr.manhattan$Type%in%c("miRNA","piRNA","sncRNA") 
                                 & start(gr.manhattan)>1e8+.8e6 &end(gr.manhattan)<1e8+1.2e6,]
seqlevels(gr.manhattan.c14mc)<-"14"
gr.manhattan.c14mc$Type<-droplevels(gr.manhattan.c14mc$Type)
p.man.c14mc<-ggbio::autoplot(gr.manhattan.c14mc, geom="point", aes(y=score,label=ID), color="gray45",size=2.2, alpha=.95, spaceline = T) + 
            xlab("Chromosome 14") +
            facet_grid(Type~.) + 
            theme_Publication()

foo.c14mc<-ggplot_build(p.man.c14mc@ggplot)
bar.c14mc<-cbind(data.table(foo.c14mc$data[[1]]), ID=gr.manhattan.c14mc$ID, RANK=gr.manhattan.c14mc$RANK, INDEX_PCT=gr.manhattan.c14mc$INDEX_PCT)
bar.c14mc[,Type:=ifelse(PANEL==1,"miRNA",ifelse(PANEL==2,"piRNA","sncRNA"))]
p.c14mc<-p.man.c14mc@ggplot + 
    #facet_wrap(~Type,nrow=3,scales="free_y",strip.position="right") +
    geom_point(data=bar.c14mc[Type=="miRNA" & INDEX_PCT<=5], aes(x,y),col="cornflowerblue",size=3,alpha=.8) +# miRNA
    geom_point(data=bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1], aes(x,y),col="orange2",size=3) +# piRNA
    geom_point(data=bar.c14mc[Type=="sncRNA" & INDEX_PCT<=1], aes(x,y),col="orange2",size=3) + # sncRNA
    scale_y_continuous(breaks=c(0,0.5,1)) +
    ggrepel::geom_text_repel(data=bar.c14mc[Type=="miRNA" & INDEX_PCT<=5],aes(x,y,label=ID),direction="x",
                             nudge_y=0.15-bar.c14mc[Type=="miRNA" & INDEX_PCT<=5]$y,
                             angle= 90,vjust= 0,segment.size=0.2, segment.color = "grey50", size=3.7) +

    ggrepel::geom_text_repel(data=bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & x<100900000], aes(x,y,label=ID), direction="y", hjust=1, 
                             nudge_x=bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & x<100900000]$x-100930000,
                             segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & RANK<=7], aes(x,y,label=ID), direction="y", hjust=0, 
                             nudge_x=100970000-bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & RANK<=7]$x,
                             segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & RANK>=18], aes(x,y,label=ID), direction="y", hjust=0, 
                             nudge_x=101060000-bar.c14mc[Type=="piRNA" & INDEX_PCT<=.1 & RANK>=18]$x, nudge_y=0.2,
                             segment.size=0.2, segment.color = "grey50", size=3.7) +

    ggrepel::geom_text_repel(data=bar.c14mc[Type=="sncRNA" & INDEX_PCT<=1 & x<100960000], aes(x,y,label=ID), direction="y", hjust=1, 
                             nudge_x=100900000-bar.c14mc[Type=="sncRNA" & INDEX_PCT<=1 & x<100960000]$x, nudge_y=0.2,
                             segment.size=0.2, segment.color = "grey50", size=3.7) +
    ggrepel::geom_text_repel(data=bar.c14mc[Type=="sncRNA" & INDEX_PCT<=1 & x>100960000], aes(x,y,label=ID), direction="y", hjust=0, 
                             nudge_x=101000000-bar.c14mc[Type=="sncRNA" & INDEX_PCT<=1 & x>100960000]$x, nudge_y=0.2,
                             segment.size=0.2, segment.color = "grey50", size=3.7) +
    ylab("") + xlab("C14MC") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

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
    ggrepel::geom_text_repel(data=bar[order(RANK)][7], aes(x,y,label=ID), direction="y", hjust=0, nudge_x=-10^8/3,segment.size=0.2, segment.color = "grey50", size=3.5) +
    ylab("") +
    theme(plot.margin = unit(c(8,5.1,5,8), "pt"))

###########################
# combine A,B,C,D,E above #
###########################
file.name<-file.path("Figures/Fig2.manhattan.rev1.tiff")
tiff(filename=file.name,width=10, height=15,units="in",res=300, compression = 'lzw')
cowplot::plot_grid(p.all, p.chr19, p.chrMT, p.chr14, p.chr1, 
                   labels="auto",label_size=20,ncol=1, align="v", 
                   rel_heights=c(3.3,1,1,1.7,1.1))
dev.off()

###################
## C19MC & C14MC ##
###################
file.name<-file.path("Figures/Fig2.manhattan.rev1.mc.tiff")
tiff(filename=file.name,width=10, height=14,units="in",res=300, compression = 'lzw')
cowplot::plot_grid(p.all, p.c19mc, p.c14mc, p.chrMT,
                   labels="auto",label_size=20,ncol=1, align="v", 
                   rel_heights=c(3.5,1,3,1))
dev.off()
