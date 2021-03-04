source("lib/local.R")
load("RData/dl.deseq.RData") # 12M (dl.deseq.anno, dl.abundance, dl.cum.fpkm)
load("RData/dl.manhattan.RData") # loads dl.manhattan (for miRNA, piRNA, circRNA and sncRNA)
load("RData/dt.ensg.RData") 

# chr19:53,641,443-53,780,750
dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750][order(RANK)]
dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750 & INDEX_PCT<=5][order(RANK)]
dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750 & INDEX_PCT<=5,.N,ID]
my.target.mir<-dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750,.N,ID]$ID
length(my.target.mir)
c19mc.top5.pct<-dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750 & INDEX_PCT<=5,.N,ID]$ID
length(c19mc.top5.pct)

####################
## miRTarbase v8.0 #
####################
dt.mtb<-data.table(readxl::read_excel("~/results/RNA-Seq/Placentome/ObsGynaeNetworkAnalysisData/miRTarBase/miRTarBase.8.0/miRTarBase_SE_WR_hsa.xlsx"))
dt.mtb[,.N,miRNA][order(-N)]
dt.mtb[,.(`miRTarBase ID`,miRNA,`Target Gene`,`Support Type`,`References (PMID)`)]
merge(dt.mtb, dt.ensg[,.(hgnc_symbol,gene_biotype)],by.x='Target Gene',by.y="hgnc_symbol")[,.N,gene_biotype]
dt.mtb[,No_evi:=lengths(regmatches(Experiments, gregexpr("//", Experiments)))+1] # count how many evi
dt.mtb<-dt.mtb[order(`miRTarBase ID`,miRNA,`Target Gene`,-No_evi)] # order by No of evi
dt.mtb.f<-dt.mtb[No_evi>2 & `Support Type`=="Functional MTI",.SD[c(1)],.(`miRTarBase ID`,miRNA,`Target Gene`)][,-c("Species (miRNA)","Target Gene (Entrez ID)","Species (Target Gene)","Support Type")] # get the entry having the best evi

dt.mtb.f[,.N,miRNA][order(-N)]
dt.mtb.f[miRNA=="hsa-miR-34a-5p",-"Experiments"]
dt.mtb.f[miRNA %in% my.target.mir,.N,`miRTarBase ID`]
dt.mtb.f[miRNA %in% my.target.mir,.N,miRNA]
dt.mtb[miRNA %in% c19mc.top5.pct,.N,miRNA]
dt.mtb.f[miRNA %in% c19mc.top5.pct,.N,miRNA]

fwrite(dt.mtb.f, file=gzfile("~/tmp/filtered.miRTarBase_SE_WR_hsa.csv.gz"))
write.csv(dt.mtb.f, file=gzfile("~/tmp/filtered.miRTarBase_SE_WR_hsa.csv.gz"),quote=F)

##
dt.mir<-merge(
    #dl.manhattan[["miRNA"]][seqnames=="chr19" & start>=53641443 & end<=53780750][,.N,.(ID,meanCount,meanFpkm,RANK,INDEX_PCT)],
    dl.manhattan[["miRNA"]][,.N,.(`miRNA`=ID,meanCount,meanFpkm,RANK,INDEX_PCT)],
    dt.mtb.f,
    by="miRNA")

## POPS total RNA-Seq
#dt.target<-dl.deseq.anno[["total-RNA-Seq"]][gene_biotype=="protein_coding"][order(-meanFpkm)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)]
dt.target<-dl.deseq.anno[["total-RNA-Seq"]][hgnc_symbol %in% dt.mtb.f$`Target Gene`][order(-meanFpkm)][,`:=`(RANK=seq_len(.N),INDEX_PCT=(seq_len(.N)/.N)*100)]

#
dt.mir.target<-merge(
                     dt.mir,
                     dt.target[,.(hgnc_symbol,ensembl_gene_id,baseMean,meanFpm,meanFpkm,gene_biotype,RANK,INDEX_PCT)],
                     by.x="Target Gene",by.y="hgnc_symbol")[order(RANK.x)]

##
dt.mir.target[miRNA %in% c19mc.top5.pct,-"Experiments"]
dt.mir.target[miRNA %in% c19mc.top5.pct,.N,miRNA][order(-N)]
dt.mir.target[miRNA %in% c19mc.top5.pct,.N,.(`Target Gene`,RANK.y,meanFpkm.y)][order(RANK.y,meanFpkm.y,-N)]
dt.mir.target[miRNA=="hsa-miR-519a-3p"][order(RANK.y)][,-"Experiments"]
dt.mir.target[`Target Gene`=="CDKN1A"][order(RANK.y)][,-"Experiments"]


p1<-ggplot(dt.mir.target, aes(meanFpkm.x,meanFpkm.y)) + 
    #geom_point(data=dt.mir.target[miRNA %in% my.target.mir & !miRNA %in% c19mc.top5.pct], shape=1,col="grey",size=2) + 
    geom_point(data=dt.mir.target[miRNA %in% c19mc.top5.pct], shape=1,col="grey10",size=5) + 
    geom_smooth(data=dt.mir.target[miRNA %in% c19mc.top5.pct], size=2, alpha=.5, method="lm") +
    scale_x_continuous(trans = log10_trans(),
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = log10_trans(),
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
    annotate("text",x=10^5.8,y=10^2,label="Pearson=-0.5\n(P=0.002)",size=5) +
    #annotate("rect",xmin=10^5.5,xmax=10^6,ymin=10^1.7,ymax=10^2.5,alpha=.2) +
    xlab("Expression level of miRNA (RPKM)") +
    ylab("Expression level of mRNA (RPKM)") +
    theme_Publication()
print(p1)

cor.test(log2(dt.mir.target[miRNA %in% c19mc.top5.pct]$meanFpkm.x), log2(dt.mir.target[miRNA %in% c19mc.top5.pct]$meanFpkm.y), method="spearman")
cor.test(log2(dt.mir.target[miRNA %in% c19mc.top5.pct]$meanFpkm.x), log2(dt.mir.target[miRNA %in% c19mc.top5.pct]$meanFpkm.y),method="pearson")


##
##
library(biomaRt)
my.target.ensg<-as.character(dt.mir.target[miRNA %in% c19mc.top5.pct,.N,ensembl_gene_id]$ensembl_gene_id)
length(my.target.ensg)
myMart= biomaRt::useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
dt.target.desc<-data.table(
                    biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"), filters = "ensembl_gene_id", values = my.target.ensg, mart = myMart)
                           )

foo<-merge(
    dt.mir.target[miRNA %in% c19mc.top5.pct][,.(miRNA,`Target Gene`,FPKM_miRNA=meanFpkm.x,FPKM_target=meanFpkm.y,`miRTarBase ID`,Experiments,`References (PMID)`)],
    dt.target.desc
    ,by.x="Target Gene",by.y="hgnc_symbol"
    )[,.(miRNA,`Target Gene`,ensembl_gene_id,description,FPKM_miRNA,FPKM_target,`miRTarBase ID`,Experiments,`References (PMID)`)]
fwrite(foo,file="data/c19mc.miRNAs.target.miRTarBase.csv")

##
## GO
if(TRUE){
	library(goseq) # for goseq
	library(GO.db) # for GOBPPARENTS
	library(GOstats) # for GOGraph

    dt.gP<-fread("data/gProfiler_hsapiens_03-11-2020_15-21-45__intersections.csv")
    dt.gP[,.N,source]
    dt.gP.go<-dt.gP[grepl("GO:",source)]
    
    dt.gP.go[["depth"]]<-as.numeric(sapply(dt.gP.go$term_id, function(i) {
                if(is.na(Ontology(i))){
                    0
                }else{
                    if(Ontology(i)=="BP"){
                        graph::DFS(GOGraph(i,GOBPPARENTS), "all")[i]
                    }else if(Ontology(i)=="MF"){
                        graph::DFS(GOGraph(i,GOMFPARENTS), "all")[i]
                    }else{
                        graph::DFS(GOGraph(i,GOCCPARENTS), "all")[i]
                    }
                }
            } # end of function(i)
        )# end of sapply
    )
    dt.gP.go[intersection_size>=10,.N,.(source,depth)][order(source,depth)]
    dt.foo<-dt.gP.go[depth==6 & intersection_size>=10][,term:=ifelse(nchar(term_name)>80,paste0(substr(term_name,1,80),"..."),term_name)]
    p.gop<-ggplot(dt.foo, aes(term,-log10(adjusted_p_value))) +
        geom_bar(aes(fill=source),stat="identity",width=.7, alpha=.6) +
        coord_flip() +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey30", size=1) +
        #ggsci::scale_fill_d3(alpha=.5) + 
        scale_fill_manual(values=c(`GO:BP`=ggsci::pal_d3()(3)[1],
                                `GO:CC`=ggsci::pal_d3()(3)[2],
                                `GO:MF`=ggsci::pal_d3()(3)[3])) +
        xlab("") + ylab("-log10(adj. P)") + 
        scale_x_discrete(limits=dt.foo[order(source,-adjusted_p_value)]$term) +
        theme_Publication()
    print(p.gop)

    file.name<-file.path('Figures/SI.Fig.c19mc.target.go.jpg')
    jpeg(file=file.name, width=12, height=8,units="in",res=150) 
    cowplot::plot_grid(p1, 
                       p.gop+
                           theme(axis.text=element_text(size=rel(.9)),
                                 legend.position="top",
                                 legend.title=element_blank(),
                                 legend.text=element_text(size=rel(.5)),
                                 legend.key.width=unit(.3,"cm"), 
                                 legend.key.height=unit(.5,"cm")), 
                       labels="auto",label_size=25,nrow=1,rel_widths=c(1,1.65))
    dev.off()
}

if(FALSE){
	cat("doing GO analysis...\n")
	library(goseq) # for goseq
	library(GO.db) # for GOBPPARENTS
	library(GOstats) # for GOGraph

    myGenes<-as.integer(
               #as.character(dt.target$ensembl_gene_id) 
                as.character(dl.deseq.anno[["total-RNA-Seq"]][gene_biotype=="protein_coding"]$ensembl_gene_id)
                        %in%
                   as.character(dt.mir.target[miRNA %in% my.target.mir,.N,ensembl_gene_id]$ensembl_gene_id)
               )
    #names(myGenes)<-as.character(dt.target$ensembl_gene_id)
    names(myGenes)<-as.character(dl.deseq.anno[["total-RNA-Seq"]][gene_biotype=="protein_coding"]$ensembl_gene_id)
    table(myGenes)

	pwf=goseq::nullp(myGenes,'hg19',"ensGene") # isa data.frame  (check available ID via supportedGeneIDs())
    GO.wall=goseq::goseq(pwf,'hg19',"ensGene") # isa data.frame

    dt.go<-data.table(GO.wall)
    dt.go.sure<-rbind(
        dt.go[,`padj` := p.adjust(over_represented_pvalue, method="BH")][padj<0.05][,type:="over"],
        dt.go[,`padj` := p.adjust(under_represented_pvalue, method="BH")][padj<0.05][,type:="under"]
        )
    dt.go.sure[["depth"]]<-as.numeric(sapply(dt.go.sure$category, function(i) {
                if(Ontology(i)=="BP"){
                    graph::DFS(GOGraph(i,GOBPPARENTS), "all")[i]
                }else if(Ontology(i)=="MF"){
                    graph::DFS(GOGraph(i,GOMFPARENTS), "all")[i]
                }else{
                    graph::DFS(GOGraph(i,GOCCPARENTS), "all")[i]
                }
            } # end of function(i)
        )# end of sapply
    )
    dt.go.sure[type=="over",.N,ontology]
    dt.go.sure[type=="over" & ontology=="BP",.N,depth]
    dt.go.sure[type=="over" & depth==7 & ontology=="MF"]

    p.go.over<-ggplot(dt.go.sure[type=="over" & depth==8 & ontology=="BP"][order(padj)], aes(term,-log10(padj))) +
        geom_bar(stat="identity",width=.7,fill="grey") +
        coord_flip() +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey30", size=1) +
        xlab("") + ylab("-log10(adjusted p-value)") + #ggtitle("over") +
        #facet_grid(~ontology) +
        theme_Publication() +
        #theme(axis.text = element_text(size=rel(1.2)))
        theme(axis.text = element_text(size=rel(1)))
    print(p.go.over)
}
fwrite(dl.deseq.anno[["total-RNA-Seq"]][gene_biotype=="protein_coding",.(ensembl_gene_id)], file="~/tmp/ensg_protein_coding.txt", col.names=F)
