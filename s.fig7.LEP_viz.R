#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# first created 8/Aug/2018
# last modified 8/Aug/2018

minFpkm=0.1
TR_PREFIX='GRCh38' # GRCh37|GRCh38
ENS_VER=82 # re-construction is based on GRCh38.82 which is used for StringTie
mySource="Placentome" 
myCohort="POPS" 
library(Gviz)
myGene="LEP" # LEP

    #################
    # Reference LEP #
    #################
    my.file=file.path("data/Homo_sapiens.GRCh38.82.LEP.gtf")
    txdb.lep<-GenomicFeatures::makeTxDbFromGFF(my.file) 
    gv.lep<-GeneRegionTrack(txdb.lep, name="", fontcolor.group="black", showId=T, fontsize=15, cex.group=1)
    plotTracks(gv.lep)

    ######################
    # Re-constructed LEP #
    ######################
    my.file=file.path('data/POPS.LEP.GRCh38.82.reconstructed.tr.gtf')
    gr.pops.lep<-rtracklayer::import(my.file)
    txdb.pops.lep<-GenomicFeatures::makeTxDbFromGFF(my.file) 
   
    gv.pops.lep<-GeneRegionTrack(txdb.pops.lep,name="Re-constructed LEP transcripts", col.title="black",fontcolor.group="black", fill="cornflowerblue", showId=T, fontsize=17, cex.group=1)
    plotTracks(gv.pops.lep)


    #######################
    # Coverage (bedgraph) #
    #######################
    my.file=file.path('data/POPS.GRCh38.82.LEP.tr.500.bedgraph')
    gv.lep.cov<-DataTrack(my.file, type="histogram", col.title="black",name="Coverage",fontsize=20)
    plotTracks(gv.lep.cov)
    plotTracks(gv.lep.cov, from = 128239000,to=128258000)


    #############
    ## Sequence #
    #############
    library(BSgenome.Hsapiens.UCSC.hg38)
    #library(BSgenome.Hsapiens.NCBI.GRCh38)
    seqinfo(Hsapiens)
    gv.strack <- SequenceTrack(Hsapiens, cex=1.3)

    ################################### 
    ## put individual tracks together #
    ################################### 
    # entire LEP loci
    gv.axis<-GenomeAxisTrack(add53=T, cex=1.3, littleTicks=T)
    # highlight 1st/3rd exon
    gv.hl1<-HighlightTrack(trackList=list(gv.lep, gv.pops.lep, gv.lep.cov), start=c(128241180, 128257600), end=c(128241380,128257680))
    # highlight 1st/2nd/3rd exon
    gv.hl2<-HighlightTrack(trackList=list(gv.lep, gv.pops.lep, gv.lep.cov), start=c(128241180, 128254350, 128257600), end=c(128241380, 128254410, 128257680))
    # highlight 2nd exon
    gv.hl3<-HighlightTrack(trackList=list(gv.lep, gv.pops.lep, gv.lep.cov), start=c(128254350), end=c(128254410))

    file.name<-file.path("Figures/LEP.ref.rpt.cov.a.tiff")
    tiff(file.name, width=12, height=4,units="in",res=300, compression = 'lzw') #A4 size
    # highlight 1st exon and the 3' of the 3rd exon
    plotTracks(list(gv.axis, gv.hl1),
               main="a                                                                                                                             ")
    dev.off()

    file.name<-file.path("Figures/LEP.ref.rpt.cov.a2.tiff")
    tiff(file.name, width=12, height=4,units="in",res=300, compression = 'lzw') #A4 size
    # highlight 1st exon and the 3' of the 3rd exon
    plotTracks(list(gv.axis, gv.hl2),
               main="a                                                                                                                             ")
    dev.off()

    file.name<-file.path("Figures/LEP.ref.rpt.cov.a.4x1.tiff")
    tiff(file.name, width=12, height=3,units="in",res=300, compression = 'lzw') #A4 size
    # highlight 1st exon and the 3' of the 3rd exon
    plotTracks(list(gv.axis, gv.hl3),
               main="a                                                                                                                             ")
    dev.off()


    # 1st exon 
    gv.axis<-GenomeAxisTrack(fontsize=10)
    file.name<-file.path("Figures/LEP.ref.rpt.cov.1st.exon.b.tiff")
    tiff(file.name, width=6, height=4,units="in",res=300, compression = 'lzw') #A4 size
    plotTracks(c(gv.axis, gv.lep, gv.pops.lep, gv.lep.cov),from=128241180, to=128241380, extend.left=50, fontsize=12.5,
               main="b                                                            ")
    dev.off()

    # 3rd exon
    gv.axis<-GenomeAxisTrack(fontsize=10)
    file.name<-file.path("Figures/LEP.ref.rpt.cov.3rd.exon.c.tiff")
    tiff(file.name, width=6, height=4,units="in",res=300, compression = 'lzw') #A4 size
    plotTracks(c(gv.axis, gv.lep, gv.pops.lep, gv.lep.cov),from=128257500, to=128257700, just.group="right", extend.right=50, fontsize=12.5,
               main="c                                                            ")
    dev.off()

    # the first 3-bases from the 3rd exon 
    gv.axis<-GenomeAxisTrack(fontsize=13)
    gv.lep<-GeneRegionTrack(txdb.lep, name="", fontcolor.group="black", showId=T, fontsize=15, cex.group=0.5)
    gv.pops.lep<-GeneRegionTrack(txdb.pops.lep,name="Re-constructed LEP transcripts", col.title="black",fontcolor.group="black", fill="cornflowerblue", showId=T, fontsize=13, cex.group=0.5)
    my.file=file.path("data/POPS.GRCh38.82.LEP.tr.500.bedgraph")
    file.name<-file.path("Figures/LEP.ref.rpt.cov.Gln.b.3x1.tiff")
    tiff(file.name, width=12, height=4,units="in",res=300, compression = 'lzw') #A4 size
    plotTracks(c(gv.axis, gv.lep, gv.pops.lep, gv.lep.cov, gv.strack),from=128254400, to=128254410,
               main="b                                                                                                                             ")
    dev.off()

    # A: 12 * 4
    # B/C: 6 * 4
    # D: 12 * 3
    # A-(B+C)-D: 12 * (4+4+3)
    # A2-(B+C)-D2: 12 * (4+4+4)
