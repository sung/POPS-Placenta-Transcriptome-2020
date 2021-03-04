minFpkm=0.1
TR_PREFIX='GRCh38' # GRCh37|GRCh38
ENS_VER=82 # re-construction is based on GRCh38.82 which is used for StringTie
#ENS_VER=94 # re-construction is based on GRCh38.94 which is used for StringTie
mySource="Placentome" # "FG"
myCohort="POPS" # "CTLM", "CTLF", "POPS", "CTL", "PET", "SGA"
myProject=paste(mySource,myCohort,sep=".") # e.g. Placentome.CTL
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"
