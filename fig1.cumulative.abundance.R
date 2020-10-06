source("lib/local.R")
load("RData/dl.deseq.RData") # 12M (dl.deseq.anno, dl.abundance, dl.cum.fpkm)
#load("RData/dl.gtex.deseq.RData") # 66M # TODO:  <30-09-20, Sung> # error while loading - rebuild
load("RData/dt.mapped.smallRNA.RData") # 80M 

############################################
## Pie chart: a) totla-RNA, b) small-RNA   #
## Percentage of Mapped Reads by RNA types #
############################################
# Total-RNA 
dt.biotype<-data.table(
    rbind(
        cbind("mRNA","protein_coding"),
        cbind("lincRNA","lincRNA"),
        cbind("pseudogene",c("processed_pseudogene","unprocessed_pseudogene","transcribed_unprocessed_pseudogene","transcribed_processed_pseudogene","pseudogene")),
        cbind("Mt rRNA","Mt_rRNA"),
        cbind("misc RNA","misc_RNA")
        )
    )
setnames(dt.biotype, c("RNA type","gene_biotype"))
sumT<-dl.deseq.anno[["total-RNA-Seq"]][,sum(meanCount)]
dt.foo<-merge(dt.biotype,  dl.deseq.anno[["total-RNA-Seq"]], all.y=T) # n=60619
#dt.foo[`RNA type`=="misc RNA"]
#dt.foo[is.na(`RNA type`),.N,gene_biotype][order(-N)]

dt.total<-dt.foo[,.(`totalCount`=sum(meanCount),`percentage`=round(sum(meanCount)/sumT*100,1)),.(`RNA type`)][order(percentage)]
dt.total[is.na(`RNA type`),`RNA type`:="others"]

my.levels<-c("mRNA","lincRNA","pseudogene","Mt rRNA","misc RNA","others")
dt.total<-rbindlist(lapply(rev(my.levels), function(i) dt.total[`RNA type`==i]))
dt.total$`RNA type`<-factor(dt.total$`RNA type`, levels=my.levels)

# https://stackoverflow.com/questions/44436214/ggplot-pie-chart-labeling/44438500
# https://stackoverflow.com/questions/46277894/ggplot2-pie-plot-with-geom-text-repel
dt.total[,pos:=cumsum(percentage) - percentage/2]

p.pie.total<-ggplot(dt.total, aes(x=2, y=percentage, fill=`RNA type`)) + 
    xlim(0.5, 2.5) +
    geom_bar(stat="identity", width=1,alpha=.8) +
    coord_polar("y", start=0) + 
    geom_text(data=dt.total[`RNA type`%in%c("mRNA")],
                aes(label=paste0(`RNA type`," (",percentage, "%)")), 
                position = position_stack(vjust = .66),
                size=6,col="grey90") +
    ggrepel::geom_text_repel(data=dt.total[!`RNA type`%in%c("mRNA")], 
                                aes(x=2.5,y=pos,label=paste0(`RNA type`,"(",percentage, "%)")), 
                                size=6, nudge_x = .9, segment.size = .4, show.legend = T, segment.color = "grey10") +
    labs(x = NULL, y = NULL, fill = NULL) +
    annotate("text", x=.5,y = 85, label = "Total RNA-Seq",size=8) +
    scale_fill_manual(values=c(cbPalette2[1:3],"grey30","grey50","grey70")) +
    theme_void() +
    #theme(legend.position="")
    theme(legend.position="", plot.margin = unit(c(.5,0,-1.7,0), "cm"))

##############
## Small-RNA #
##############
dt.mapped.smallRNA[,.(count=sum(depth))]
print(dt.small.v1<-rbind(
    dt.mapped.miRNA[,.(`RNA type`="miRNA",count=sum(depth))],
    dt.mapped.piRNA[,.(`RNA type`="piRNA",count=sum(depth))],
    dt.mapped.exon[,.(`RNA type`="other exonic",count=sum(depth))],
    #dt.unmapped.smallRNA[,.(`RNA type`="remaining",count=sum(depth))]
    data.table(`RNA type`="remaining",count=dt.unmapped.smallRNA.miR.piR[,sum(depth)] - dt.mapped.exon[,sum(depth)])
    ))
dt.small.v1[,sum(count)]

print(dt.small.v2<-rbind(
    dt.mapped.miRNA[,.(`RNA type`="miRNA",count=sum(depth))],
    dt.mapped.piRNA[,.(`RNA type`="piRNA",count=sum(depth))],
    dt.mapped.tRNA[,.(`RNA type`="tRNA",count=sum(depth))],
    dt.mapped.sncRNA[,.(`RNA type`="sncRNA",count=sum(depth))],
    dt.mapped.exon.rev1[,.(`RNA type`="other exonic",count=sum(depth))],
    data.table(`RNA type`="remaining",count=dt.unmapped.smallRNA.miR.piR[,sum(depth)] - dt.mapped.tRNA[,sum(depth)] - dt.mapped.sncRNA[,sum(depth)] - dt.mapped.exon.rev1[,sum(depth)])
    ))
dt.small.v2[,sum(count)]

print(dt.small<-rbind(
    dt.mapped.miRNA[,.(`RNA type`="miRNA",count=sum(depth))],
    dt.mapped.piRNA[,.(`RNA type`="piRNA",count=sum(depth))],
    dt.mapped.tRNA[,.(`RNA type`="tRNA",count=sum(depth))],
    dt.mapped.sncRNA[,.(`RNA type`="sncRNA",count=sum(depth))],
    dt.mapped.exon.rev1[,.(`RNA type`="other exonic",count=sum(depth))],
    data.table(`RNA type`="remaining",count=dt.mapped.smallRNA[,sum(depth)]- dt.mapped.miRNA[,sum(depth)] - dt.mapped.piRNA[,sum(depth)] - dt.mapped.tRNA[,sum(depth)] - dt.mapped.sncRNA[,sum(depth)] - dt.mapped.exon.rev1[,sum(depth)])
    ))
dt.small[,sum(count)]

sumS<-dt.small[,sum(count)]
dt.small[,`percentage`:=round(count/sumS*100,1)]

my.levels<-c("miRNA","piRNA","sncRNA","tRNA","other exonic","remaining")
dt.small<-rbindlist(lapply(rev(my.levels), function(i) dt.small[`RNA type`==i]))
dt.small$`RNA type`<-factor(dt.small$`RNA type`, levels=my.levels)
dt.small

# https://stackoverflow.com/questions/44436214/ggplot-pie-chart-labeling/44438500
# https://stackoverflow.com/questions/46277894/ggplot2-pie-plot-with-geom-text-repel
dt.small[,pos:=cumsum(percentage) - percentage/2]
dt.small$`RNA type`

p.pie.small<-ggplot(dt.small, aes(x=2, y=percentage, fill=`RNA type`)) + 
    xlim(0.5, 2.5) +
    geom_bar(stat="identity", width=1,alpha=.85) +
    coord_polar("y", start=0) + 
    geom_text(data=dt.small[`RNA type`=="miRNA"],aes(label=paste0(`RNA type`," (",percentage, "%)")), 
                position = position_stack(vjust = 0.55),size=6) +
    ggrepel::geom_text_repel(data=dt.small[!`RNA type`=="miRNA"], 
                                aes(x=2.5,y=pos,label=paste0(`RNA type`,"(",percentage, "%)")), 
                                size=6, nudge_x = .4, segment.size = .4, show.legend = T, segment.color = "grey10") +
    labs(x = NULL, y = NULL, fill = NULL) +
    annotate("text", x=.5,y = 85, label = "Small RNA-Seq",size=8) +
    #scale_fill_manual(values=c(cbPalette[4],cbPalette[5],cbPalette[1],"grey80")) +
    scale_fill_manual(values=c(cbPalette[4],cbPalette[5],cbPalette[6],cbPalette[7],cbPalette[1],"grey80")) +
    theme_void() +
    theme(legend.position="", plot.margin = unit(c(.5,0,-1.7,0), "cm"))

##################################
## Figure: Transcript abundance ##
##################################
dl.abundance<-lapply(dl.abundance, function(i) i[,meanTPM:=meanFpkm/sum(meanFpkm)*10^6]) # add TPM
dt.abundance<-rbindlist(dl.abundance) # dl.abundance[-1]: "pre-miRNA"
dt.abundance[`RNA Type`=="protein_coding",`RNA Type`:="mRNA"]
dt.abundance$`RNA Type`<-factor(dt.abundance$`RNA Type`, levels=c("mRNA","lincRNA","pseudogene","miRNA","piRNA","sncRNA","pre-miRNA","piRNA-miRNA"))

fpkm.cutoff<-c(0,10^-3, 10^-2 ,10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)
dl.abundance.summary<-list()
for(i in names(dl.abundance)){
    dl.abundance.summary[[i]]<-data.table(`RNA Type`=i,
                                fpkm.cutoff,
                                `Count.FPKM`=sapply(fpkm.cutoff, function(j) dl.abundance[[i]][meanFpkm<=j,.N]), # max this level
                                `Count.TPM`=sapply(fpkm.cutoff, function(j) dl.abundance[[i]][meanTPM<=j,.N]), # max this level
                                `Total`=dl.abundance[[i]][,.N]
                        )
}
dt.abundance.summary<-rbindlist(dl.abundance.summary)
dt.abundance.summary[`RNA Type`=="protein_coding",`RNA Type`:="mRNA"]
#dt.abundance.summary$`RNA Type`<-factor(dt.abundance.summary$`RNA Type`, levels=c("mRNA","lincRNA","pseudogene","miRNA","piRNA","pre-miRNA"))
dt.abundance.summary$`RNA Type`<-factor(dt.abundance.summary$`RNA Type`, levels=c("mRNA","lincRNA","pseudogene","miRNA","piRNA","sncRNA","pre-miRNA","piRNA-miRNA"))

# Figure 1. Transcript Abundance
p.a<-ggplot(dt.abundance.summary[`RNA Type`!="pre-miRNA"], aes(fpkm.cutoff+0.0001, Count.FPKM/Total*100)) +
    geom_line(aes(col=`RNA Type`),size=2.7,alpha=.85) +
    geom_point(size=4.5,alpha=.85,) +
    scale_x_continuous(trans = log10_trans(),
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values=cbPalette2) +
    xlab("Maximum abundance level (FPKM)") + 
    ylab("Percentage of transcripts (%)") + 
    theme_Publication() +
    theme(legend.position=c(0.8,0.25), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size = rel(1.1)), legend.background=element_rect(fill=alpha("white",0.7)))

p.b<-ggplot(dt.abundance[`RNA Type`!="pre-miRNA"], aes(meanFpkm+0.0001)) +
    geom_line(aes(col=`RNA Type`),stat="density", size=2.7, alpha=.85) +
    scale_x_continuous(trans = log10_trans(),
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(name="RNA types", values=cbPalette2) +
    xlab("Abundance level (FPKM)") + 
    ylab("Density") + 
    theme_Publication() +
    theme(legend.position="none") 

#######################################################
## Figure: Transcript complexity (tissue comparision) #
#######################################################
dt.cum.fpkm<-rbindlist(dl.cum.fpkm) # FPKM>0.1 only
dt.cum.fpkm[,`:=`(INDEX=seq_len(.N),COUNT=.N,INDEX_PCT=(seq_len(.N)/.N)*100),`RNA Type`] # index 
dt.cum.fpkm[`RNA Type`=="protein_coding",`RNA Type`:="mRNA"]
dt.cum.fpkm$`RNA Type`=factor(dt.cum.fpkm$`RNA Type`,levels=c("mRNA","lincRNA","pseudogene","miRNA","piRNA","sncRNA","pre-miRNA"))
dt.cum.fpkm[,.(.N,max(COUNT)),`RNA Type`]

p.c<-ggplot(dt.cum.fpkm[`RNA Type`!="pre-miRNA"], aes(index.pct,cumsum.fpkm.pct,group=`RNA Type`)) +
    geom_point(aes(col=`RNA Type`),size=3,alpha=.8) +
    geom_hline(yintercept=50,linetype="longdash",size=1.1,alpha=.8) +
    scale_x_continuous(trans = log10_trans(),
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(name="RNA types",values=cbPalette2) +
    xlab("Percentage of transcripts (%)") + 
    ylab("Sum of transcript abundance (%)") +
    theme_Publication() +
    theme(legend.position="none") 

###########################
## comparision with GTEx ##
###########################
if(my.gtex.ver=="v8.p2"){
    load("RData/v8.p2.dt.tpm.tau.RData")
    dt.tpm.tau
    colnames(dt.tpm.tau)
    dt.cum.fpkm.gtex<-melt.data.table(dt.tpm.tau[gene_biotype=="protein_coding" & !chromosome_name %in% c("MT","Y"),c(2,3,6:56)],
                    id.vars=c("ensembl_gene_id","hgnc_symbol"),variable.name="Tissue",value="meanFpkm")[Tissue!="Placenta.T3.pA"][order(Tissue,-meanFpkm)][,cumsum.fpkm:=cumsum(meanFpkm),Tissue] # meanFPKM is actuall TPM 
    dt.cum.fpkm.gtex[Tissue=="Placenta.T3.rZ",Tissue:="Placenta"]
    dt.cum.fpkm.gtex[,.N,Tissue] # n=19,115 gene across 50 tissues
    dt.cum.fpkm.gtex$Tissue<-as.character(dt.cum.fpkm.gtex$Tissue)

}else{
    # v6.p
    my.tissue<-c("Blood","Pancreas","Liver","Muscle","Pituitary","Brain","Placenta")
    dt.gtex.deseq<-rbindlist(dl.gtex.deseq)

    # common ensembl_gene_id across 21 tissues (including placenta)
    my.ensg<-dt.gtex.deseq[!chromosome_name %in% c("MT","Y") & gene_biotype=="protein_coding",.N,ensembl_gene_id][N==21,ensembl_gene_id] 
    #length(my.ensg) # n=19,751 genes (-MT, -Y) common across 21 tissues

    dt.cum.fpkm.gtex<-dt.gtex.deseq[meanFpkm>0.1 & ensembl_gene_id %in% my.ensg, 
                                .(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype,meanFpkm,Tissue)][order(Tissue,-meanFpkm)][,cumsum.fpkm:=cumsum(meanFpkm),Tissue]
}

    dt.cum.fpkm.gtex[,`:=`(cumsum.fpkm.pct=cumsum.fpkm/sum(meanFpkm)*100),Tissue] # index 
    dt.cum.fpkm.gtex[,`:=`(INDEX=seq_len(.N),COUNT=.N,INDEX_PCT=(seq_len(.N)/.N)*100),Tissue] # index 
    dt.cum.fpkm.gtex[,`Source`:=ifelse(Tissue=='Placenta',"POP Study","GTEx")]

    # GTEx TA50
    print(dt.ta50<-dt.cum.fpkm.gtex[cumsum.fpkm.pct>50,.(min(INDEX),min(COUNT),min(INDEX_PCT),min(cumsum.fpkm.pct)),Tissue][order(V1)]) 
    print(top.tissues<-c(as.character(dt.ta50[V1<=dt.ta50[Tissue=="Placenta"]$V1]$Tissue)))
    print(bottom.tissues<-c(dt.ta50[.N-1]$Tissue,dt.ta50[.N]$Tissue))
    print(my.tissue<-c(top.tissues,bottom.tissues))

    dt.cum.fpkm.gtex[,Tissues:=ifelse(Tissue %in% my.tissue, Tissue, "Other somatic")]
    dt.cum.fpkm.gtex$Tissues<-factor(dt.cum.fpkm.gtex$Tissues, levels=c(top.tissues,"Other somatic",bottom.tissues))
    dt.cum.fpkm.gtex[,.N,Tissues]
    dt.cum.fpkm.gtex[INDEX==1,1:6][order(cumsum.fpkm.pct)]

    p.d<-ggplot(dt.cum.fpkm.gtex, aes(INDEX_PCT,cumsum.fpkm.pct)) +
        geom_point(aes(col=Tissues),size=2.5,alpha=.7) +
        geom_line(data=dt.cum.fpkm.gtex[Tissue %in% my.tissue],aes(col=Tissues),size=1,alpha=.7) +
        geom_hline(yintercept=50,linetype="longdash",size=1.1,alpha=.8) +
        scale_x_continuous(trans = log10_trans(),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x))) +
        xlab("Percentage of protein coding transcripts (%)") + 
        ylab("Sum of transcript abundance (%)") +
        theme_Publication() +
        theme(legend.position=c(0.84,0.22), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size = rel(1.1)), legend.background=element_rect(fill=alpha("white",0.7))) +
        scale_color_manual(values=c(ggsci::pal_nejm()(6),"black","darkgrey","#FFDC91FF", "#cc8b00"))

    dt.foo<-dt.cum.fpkm.gtex[INDEX_PCT<1,.(.N,max(cumsum.fpkm), max(cumsum.fpkm.pct)),.(Tissue,Source)][order(-V3)]
    p.f<-ggplot(dt.foo,aes(Tissue,V3,fill=Source)) +
        geom_bar(stat="identity",alpha=.8) +
        scale_x_discrete(limits=dt.foo$Tissue, labels=gsub("_"," ",dt.foo$Tissue)) +
        scale_fill_manual(name="Sources",values=c(cbPalette[1],cbPalette2[1])) +
        xlab("") +
        ylab("Sum of transcript abundance (%)") +
        coord_flip() +
        theme_Publication() +
        #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(.7)), legend.position=c(0.75,0.85), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size = rel(1.1)), legend.background=element_rect(fill=alpha("white",0.7)))
        theme(axis.text.y = element_text(size=rel(.75)), legend.position=c(0.75,0.85), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size = rel(1.1)), legend.background=element_rect(fill=alpha("white",0.7)))
#############################
# A4: 8.27 × 11.69 inches   #
# combine all figures above #
#############################
file.name<-file.path("Figures/Fig1.transcript.pie.abundance.complexity.rev1.tiff")
cp.all.pie<-cowplot::plot_grid(p.pie.total, p.pie.small, labels="auto",label_size=25,nrow=1)
cp.all.top<-cowplot::plot_grid(p.a, p.b, p.c, labels=c("c","d","e"), label_size=25, align="h", nrow=1)
cp.all.bottom<-cowplot::plot_grid(p.d, p.f, labels=c("f","g"),rel_widths=c(1.2,1), label_size=25, align="h", nrow=1)
cp.all<-cowplot::plot_grid(cp.all.pie, cp.all.top, cp.all.bottom, rel_heights=c(1.1,1,1.4),nrow=3)
tiff(file.name, width=16, height=21,units="in",res=300, compression = 'lzw') #A4 size
print(cp.all)
dev.off()
