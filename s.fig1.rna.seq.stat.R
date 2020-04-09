source("lib/local.R") # boiler plates

dt.pops<-fread("Data/Meta/295.samples.csv") # sample meta info (n=295)
dt.read.cnt<-fread("Data/Meta/read.cnt.csv") # rna-seq read number
dt.pops.read.cnt<-merge(dt.pops,dt.read.cnt)

dt.foo<-rbind(dt.pops.read.cnt[,list(Read=sum(Read)),"SampleName,Library,Cohort"][,`Category`:="Sequenced Read"],
              dt.pops.read.cnt[Category!="Unmapped",list(Read=sum(Read)),"SampleName,Library,Cohort"][,`Category`:="Total Mapped Read"])

p.read2<-ggplot(dt.foo, aes(Cohort,Read/10^6)) + 
    geom_boxplot(aes(fill=Category),width=.5,size=1,alpha=0.7) +
    scale_fill_manual(values=cbPalette) +
    xlab("Patient Group") +
    ylab("No. of reads (in million)") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="none")
print(p.read2)

p.read3<-ggplot(dt.foo, aes(Library,Read/10^6)) + 
    geom_boxplot(aes(fill=Category),width=.5,size=1,alpha=0.7) +
    scale_fill_manual(values=cbPalette) +
    ylab("") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="top")
print(p.read3)

# print them out
jpeg(filename=file.path(top.dir,"Figures/SI.Fig1.rna.seq.read.stat.jpeg"), width=10, height=6, units="in", res=300)
cowplot::plot_grid(p.read2, p.read3, labels=c("b","c"), rel_widths=c(2,5),label_size=27, align="h",axis="tb")
dev.off()
