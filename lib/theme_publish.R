# code from https://rpubs.com/Koundy/71792

#theme_Publication <- function(base_size=14, base_family="helvetica") {
theme_Publication <- function(base_size=14, base_family="") {
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",size = rel(1.3), hjust = 0.5),
			text = element_text(),
			panel.background = element_rect(colour = NA),
			plot.background = element_rect(colour = NA),
			#panel.border = element_rect(colour = NA),
			axis.title = element_text(face = "bold",size = rel(1.2)),
			axis.title.y = element_text(angle=90,vjust =2),
			axis.title.x = element_text(vjust = -0.2),
			axis.text = element_text(size=rel(1.1)), 
			axis.line = element_line(colour="black"),
			axis.ticks = element_line(),
			panel.grid.major = element_line(colour="#f0f0f0"),
			panel.grid.minor = element_blank(),
			legend.key = element_rect(colour = NA),
			#legend.background = element_rect(colour = 'black', fill='grey',linetype='dashed'),
			legend.position = "right",
			#legend.direction = "horizontal",
			legend.key.size= unit(0.7, "cm"),
			legend.spacing= unit(0.1, "mm"),
			legend.title = element_text(face="bold.italic",size=rel(1)),
			legend.text = element_text(size = rel(0.9),family = "sans"),
			plot.margin=unit(c(10,5,5,5),"mm"),
			strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
			strip.text = element_text(face="bold",size=rel(1.1))
          ) + if(packageVersion("ggplot2")<=2.1){theme(legend.margin = unit(0.1, "mm"))}else{theme(legend.spacing = unit(0.1, "mm"))}
	   )
      
}

scale_fill_Publication <- function(...){
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


# code from https://gist.github.com/adamwespiser/3077780
theme_publish <- function(base_size = 12) {
  structure(list(
    axis.line =         element_blank(),
    axis.text.x =       theme_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
    axis.text.y =       theme_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
    axis.ticks =        theme_segment(colour = "black", size = 0.2),
    axis.title.x =      theme_text(size = base_size, vjust = 1),
    axis.title.y =      theme_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
    
    legend.background = theme_rect(colour=NA), 
    legend.key =        theme_rect(colour = "grey80"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       theme_text(size = base_size * 0.8,family = "sans"),
    legend.title =      theme_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
    
    panel.background =  theme_rect(fill = "white", colour = NA), 
    panel.border =      theme_rect(fill = NA, colour="grey50"), 
    panel.grid.major =  element_blank(),
    panel.grid.minor =  element_blank(),
    panel.margin =      unit(0.25, "lines"),
    
    strip.background =  theme_rect(fill = "grey80", colour = "grey50"), 
    strip.text.x =      theme_text(size = base_size * 0.8),
    strip.text.y =      theme_text(size = base_size * 0.8, angle = -90),
    
    plot.background =   theme_rect(colour = NA),
    plot.title =        theme_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  ), class = "options")
}
