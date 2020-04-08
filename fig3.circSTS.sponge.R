library(GenomicRanges)
library(Biostrings)
library(ggbio)

## chrX miR sponge (X:7514882|7516290)
bs.circX<-Biostrings::readDNAStringSet("Data/circSTS/X:7514882|7516290.fa")
my.cmd1<-"grep '^>hsa-miR-5584-5p' Data/circSTS/X:7514882\\|7516290.miranda.txt | sed 's/^>//g'"
my.cmd2<-"grep '^>hsa-miR-7113-5p' Data/circSTS/X:7514882\\|7516290.miranda.txt | sed 's/^>//g'"
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

ir.circ<-IRanges(start=1,end=width(gr.circX))
ir.m1<-IRanges(start=dt.spg[miRNA=="hsa-miR-5584-5p",c.start],end=dt.spg[miRNA=="hsa-miR-5584-5p",c.end])
ir.m2<-IRanges(start=dt.spg[miRNA=="hsa-miR-7113-5p",c.start],end=dt.spg[miRNA=="hsa-miR-7113-5p",c.end])

##  ggbio circle bit
g.circ<-autoplot(ir.circ, layout="circle")
g.m1<-autoplot(ir.m1, layout="circle", fill="grey80")
g.m2<-autoplot(ir.m2, layout="circle", fill="yellow")

# save the figure
file.name<-file.path("Figures/Fig3.circSTS.sponge.tiff")
tiff(filename=file.name,width=20, height=20,units="cm",res=300, compression = 'lzw') #A4 size 

ggbio(radius=100) + 
    circle(ir.m1, layout="circle", geom="segment", size=2,col=ggsci::pal_d3("category10")(10)[1]) +
    circle(ir.m2, layout="circle", geom="segment", size=2,col=ggsci::pal_d3("category10")(10)[2]) + 
    circle(ir.circ, layout="circle", geom="segment",size=2) +
    circle(ir.circ, layout="circle", geom="scale", size=3, scale.n=8, scale.type="sci") 

dev.off()
