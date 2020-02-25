############
## colour ##
############
hmcol.RdBu <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
hmcol.GnBu <- colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
hmcol<-hmcol.RdBu

# A colorblind-friendly palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# 12 colours
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#662506", "#997a00", "#800080", "#000000")
# The palette with black:
cbPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

my.col=list(
	`Chromosome Types`=c(`Autosomes`=cbPalette[1], `X-chromosome`=cbPalette[2]),
	#`Condition`=c(`Case`="#D95F02", `Control`="#1B9E77"), # orange/green (RnBeads default)
	`Sex`=c(`M`="#D95F02", `F`="#1B9E77"), # orange/green (RnBeads default) for RNA-Seq data
	`Gender`=c(`Male`="#D95F02", `Female`="#1B9E77"), # orange/green (RnBeads default) for RoadMap via ggplot
	#`Gender`=c(`Male`="#69b4ff", `Female`="#ff69b4"), # blue/pink 
	`Expression?`=c(`female<male`="#D95F02", `female>male`="#1B9E77"), # orange/green (RnBeads default) for RoadMap via ggplot
	`Source`=c(`Plasma`="#D95F02", `Placenta`="#1B9E77"), # orange/green (RnBeads default)
    `Outcomes`=c(`CTL`=cbPalette[1],`PE`=cbPalette[2],`SGA`=cbPalette[3])
)
my.col.Condition=c(`SGA`="#F8766D", `AGA`="#00BFC4")
#my.col.Sex=c(`Boy`="#F8766D", `Girl`="#00BFC4")

############
# graphics #
############
#par(cex=0.5, cex.main=0.7, cex.lab=0.7)
#par(oma=c(0.2,0.2,0.2,0.2)) # outter margin area
#par(mar=c(5,5,2,2)) # margin area. the default is c(5, 4, 4, 2) + 0.1
myDotCol<-"#00000020" # transparent grey 

# FDR color
myCol <- c(`FDR<=0.1`="darkgreen", `0.1<FDR<=0.2`="orange", `0.2<FDR<=0.3`="purple", `0.3<FDR<=0.4`="red", `0.4<FDR<=0.6`="red4", `0.6<FDR`="snow4")

blankPlot <- ggplot()+geom_blank(aes(1,1))+
	theme(
		plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank()
		)
