source("lib/local.R") # boiler plates

load("~/results/RNA-Seq/Boy.Girl.FG.JD.mature.miRNA.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
dt.mirdeep<-dt.deseq.anno

load("~/results/RNA-Seq/Boy.Girl.FG.JD.mature.miRNA.sRNAbench.genome.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
dt.sRNAbench.g<-dt.deseq.anno

load("~/results/RNA-Seq/Boy.Girl.FG.JD.mature.miRNA.sRNAbench.lib.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
dt.sRNAbench.l<-dt.deseq.anno

my.cols<-c("mirbase_id","baseMean","meanFpm")
dt.foo<-rbind(
    cbind(`tool`="mirDeep2",dt.mirdeep[,..my.cols]),
    cbind(`tool`="sRNAbench(genome)",dt.sRNAbench.g[,..my.cols]),
    cbind(`tool`="sRNAbench(library)",dt.sRNAbench.l[,..my.cols])
    )
dt.foo[,.N,mirbase_id][N>3]
dt.foo[mirbase_id=="hsa-miR-103b"]

dt.bar<-dcast.data.table(dt.foo, mirbase_id~tool,value.var="meanFpm",fun.aggregate=mean)

## 
## sRNAbehcn(genome) vs sRNAbehcn(library)
##
summary(lm(data=dt.bar[,.(`sRNAbench(genome)`,`sRNAbench(library)`)]))$r.squared
print(dt.bar.label<-dt.bar[`sRNAbench(genome)`>10 & `sRNAbench(library)`>10 & (`sRNAbench(genome)`/`sRNAbench(library)`>3 | `sRNAbench(genome)`/`sRNAbench(library)` < 1/3)])
ggplot(dt.bar[`sRNAbench(genome)`>10 & `sRNAbench(library)`>10], aes(`sRNAbench(genome)`,`sRNAbench(library)`)) +
#ggplot(dt.bar[mirbase_id %in% dt.foo[,.(.N,mean(baseMean)),mirbase_id][V2>1]$mirbase_id], aes(`sRNAbench(genome)`,`sRNAbench(library)`)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    ggrepel::geom_text_repel(data=dt.bar.label,aes(label=mirbase_id)) +
    theme_Publication()

## 
## sRNAbehcn(genome) vs miRDeep2
##
summary(lm(data=dt.bar[,.(`sRNAbench(genome)`,`mirDeep2`)]))$r.squared
print(dt.bar.label<-dt.bar[`sRNAbench(genome)`>10 & `mirDeep2`>10 & (`sRNAbench(genome)`/`mirDeep2`>3 | `sRNAbench(genome)`/`mirDeep2` < 1/3)][order(mirbase_id)])
ggplot(dt.bar[`sRNAbench(genome)`>5 & `mirDeep2`>5], aes(`sRNAbench(genome)`,`mirDeep2`)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    ggrepel::geom_text_repel(data=dt.bar.label, aes(label=mirbase_id)) +
    theme_Publication()

## 
## sRNAbehcn(library) vs miRDeep2
##
summary(lm(data=dt.bar[,.(`sRNAbench(library)`,`mirDeep2`)]))$r.squared
print(dt.bar.label<-dt.bar[`sRNAbench(library)`>5 & `mirDeep2`>5 & (`sRNAbench(library)`/`mirDeep2`>3 | `sRNAbench(library)`/`mirDeep2` < 1/3)][order(mirbase_id)])
ggplot(dt.bar[`sRNAbench(library)`>1 & `mirDeep2`>1], aes(`sRNAbench(library)`,`mirDeep2`)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    ggrepel::geom_text_repel(data=dt.bar.label, aes(label=mirbase_id)) +
    theme_Publication()
