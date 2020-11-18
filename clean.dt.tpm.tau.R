# this is to clean up data to be used in ShinyPlacentome website
system.time(source("libs/local.R"))

my.RData<-file.path('RData',paste0(my.gtex.ver,'.dt.tpm.tau.RData'))
message(paste("loading",my.RData,"..."))
system.time(load(my.RData))

my.gene.types=c("protein_coding","lincRNA","processed_pseudogene")
dt.foo<-dt.tpm.tau[gene_biotype%in%my.gene.types,-c("Placenta.T3.pA","length","Tau.pA","GTEx_medianTPM","GTEx_maxTPM")]

data.table::setnames(dt.foo,"GTEx_meanTPM","meanTPMGTEx")
data.table::setnames(dt.foo,"Tau.rZ","Tau")
data.table::setnames(dt.foo,"Placenta.T3.rZ","Placenta")

my.ids<-c("ensembl_gene_id","chromosome_name","hgnc_symbol","gene_biotype")
my.non.ids<-colnames(dt.foo)[!colnames(dt.foo) %in% my.ids]
dt.gtex.pt.tpm.tau<-dt.foo[,lapply(.SD,round,3),my.ids,.SDcol=my.non.ids]
dt.gtex.pt.tpm.tau[1]

save(dt.gtex.pt.tpm.tau,file="RData/dt.gtex.pt.tpm.tau.RData")

##
## Placenta-enriched 
##
dt.gtex.pt.tpm.tau[Placenta>1 & Tau>.99 & Placenta > meanTPMGTEx*100,.N,gene_biotype] # Placenta-enriched transcripts
dcast.data.table(
                dt.gtex.pt.tpm.tau[Placenta>1 & Tau>.99 & Placenta > meanTPMGTEx*100,.N,.(chromosome_name,gene_biotype)]
                ,chromosome_name~gene_biotype)[order(protein_coding)]

melt.data.table(dt.gtex.pt.tpm.tau,id.vars=c(my.ids,"meanTPMGTEx"),value.name="TPM",variable.name="Tissue")[hgnc_symbol=="PSG7"]
