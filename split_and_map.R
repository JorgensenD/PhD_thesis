#Code to split phylogenies generated in BEAST and produce animated movement maps.

# The code presented uses SPDF format maps, these are now deprecated in favour of the SF format and this code may not work with newer versions of the packages loaded or with newer R format maps.

require(bezier)
require(matlib)
require(ggplot2)
require(ggtree)
require(plyr)
require(dplyr)
require(spatialEco)
require(lubridate)
require(forcats)
require(randomcoloR)
require(RColorBrewer)
require(foreach)
require(doParallel)
require(parallel)
require(ape)
require(ggpubr)
library(grid)
library(gridExtra)
require(gtable)
require(lubridate)
require(OutbreakTools)
require(devtools)
require(phylobase)
require(stringr)
require(lemon)
require(treeio)
#install_github("GuangchuangYu/ggtree")

# MAP DATA USED IN THE ANALYSIS IS ADMIN 2 LEVEL DATA FROM THE WORLD HEALTH ORGANIZATION - REMOVED FOR REPRODUCTION HERE


# BEZIER CURVE DRAWING ----------------------------------------------------

# A function to draw a Bezier curve between any two linked locations where the names of these locations are available on the linked map.
plotbezier <- function(loc1,loc2,ID,height,mapdata,name_column){
#get x and y of each polygon and add to a dataframe
coords <- list()
print(paste0(loc1,loc2))
for(i in 1:length(mapdata@polygons)){
 name <- mapdata@data[i,name_column]
 #get midpoint of each column
 coord <- mapdata@polygons[[i]]@labpt
 #add name of the region to the coordinates
 coords[[name]] <- coord
}
if(loc2==loc1){
  # control points for internal transmission (4 to give loop)
  controlx1 <- coords[[loc1]][1]+sample(c(-1,1),1)*height
  controly1 <- coords[[loc1]][2]
  controlx2 <- coords[[loc1]][1]
  controly2 <- coords[[loc1]][2]+sample(c(-1,1),1)*height
  x <- c(coords[[loc1]][1], controlx1, controlx2, coords[[loc2]][1])
  y <- c(coords[[loc1]][2], controly1, controly2, coords[[loc2]][2])
} else{
# control point calculations for lines linking two locations
m <- (coords[[loc2]][2] - coords[[loc1]][2])/(coords[[loc2]][1] - coords[[loc1]][1])
b <-  coords[[loc2]][2] - m*coords[[loc2]][1]
#print the equation
#if (b>0){bprint <- paste("+",b, sep="")} else {bprint <- paste(b)}
#paste("y=",m,"x",bprint, sep="")
#intercept +- height to give parallel line (ifelse so different for each direction)
ifelse(coords[[loc2]][1]-coords[[loc1]][1]<0, bplus <- b - height, bplus <- b + height) 
#if (bplus>0){bplusprint <- paste("+",bplus, sep="")} else {bplusprint <- paste(bplus)}
#paste("parallel line: y=",m,"x",bplusprint, sep="")
#find the midpoint between the locations
midx <- (coords[[loc2]][1] + coords[[loc1]][1])/2
midy <- (coords[[loc2]][2] + coords[[loc1]][2])/2
midm <- -1/m
midb <- midy-midm*midx
#if (midb>0){midbprint <- paste("+",midb, sep="")} else {midbprint <- paste(midb)}
#paste("perpendicular line: y=",midm,"x",midbprint, sep="")
# find intercept of perpendicular and parallel lines
A <- matrix(c(1,1,-midm,-m),2,2)
B <- c(midb,bplus)
#showEqn(A,b)
controlyx <- solve(A,B)
# standardise the height of the arc so it doesnt depend on the angle of the line
d <- sqrt((controlyx[2]-midx)^2+(controlyx[1]-midy)^2)
t <- height/d
controlx <- ((1-t)*midx+t*controlyx[2])
controly <- ((1-t)*midy+t*controlyx[1])

#set up 3 control points for the curve (start, middle and end)
x <- c(coords[[loc1]][1], controlx, coords[[loc2]][1])
y <- c(coords[[loc1]][2], controly, coords[[loc2]][2])
}
coordsp <- cbind(x,y)
curve <- pointsOnBezier(n=20,method='evenly_spaced',p=coordsp)

curvegg <- data.frame(curve$points)

curvegg$ID <- as.numeric(ID)
curvegg$START <- loc1
Curve <- list(curvegg)
names(Curve)[[1]] <- paste(loc1, loc2, sep="_")
return(Curve)
}

# TREE INPUT --------------------------------------------------------------

## Import tree
# Load in your BEAST consensus tree and set the most recent sample date.
# We assume throughout that the discrete trait in the loaded phylogeny is named 'state'.

directory_interest <- "FILEPATH"
tree <- read.beast(paste0(directory_interest, "MJ_AFP_ca.trees"))
mrsd=2022.4932
mrsd=as.Date(date_decimal(mrsd))

## Lookup function to pull the ancestor location of each node
ancestor <- function(node){
  if(!is.na(node)){
    return(tree@phylo[["edge"]][tree@phylo[["edge"]][,2]==node,][1])
  }
}
tree@data$ancestor <- mapply(ancestor, tree@data$node)



ancestor_loc <- tree@data[,c("node","state")]
names(ancestor_loc) <- c("ancestor","ancestor_loc")
tree@data <- as_tibble(join(tree@data,ancestor_loc, by="ancestor"))

#add branch length/2 to node age to get the branching point
tree@data$branching <- as.numeric(tree@data$height) + as.numeric(tree@data$length)/2
tree@data$branching <- date_decimal(decimal_date(mrsd)-tree@data$branching)


# PLOT CURVES ONLY FOR LOCATIONS WHICH ARE LINKED IN THE TREE -------------

locationmatrix <- unique(na.omit(tree@data[,c("ancestor_loc","state")]))
ID <- c(1:nrow(locationmatrix))
locationmatrix <- cbind(locationmatrix, ID)
locationmatrix$ancestor_loc <- as.character(locationmatrix$ancestor_loc)

locationmatrix$state <- as.character(locationmatrix$state)


## run through the matrix of pairs with mapply - need to specify the name of your map and the number of the column where the region names will be found
ptm <- proc.time()
out <- mapply(plotbezier,locationmatrix[,1], locationmatrix[,2], locationmatrix[,3],  MoreArgs = list(height=2, mapdata=MTT3_map, name_column=1))
proc.time() - ptm


# Code to generate animated maps ------------------------------------------

# custom palette for the regions included
getPalette <- c("#7735c6", "#26a0d0", "#9b0058", "#007d32", "#dc241f", "#ffd329", "#0019a8", "#04d690", "#ed9993","#ef7b10")
names(getPalette) <- c("SOUTH-CORRIDOR", "EAST-PAK", "CENTRE-PAK", "KARACHI", "CENTRAL-CORRIDOR", "EAST-AFG", "NORTH-CORRIDOR", "WEST-AFG", "WEST-PAK", "NORTH-AFG")
newpal <- getPalette

## Generate background maps to plot arrows over

mo <- seq(from=as.Date(floor_date(min(tree@data$branching), "month")), to=mrsd,by='months')


# Centroid points for each location
coords <- list()
for(i in 1:length(MTT3_map@polygons)){
  name <- MTT3_map@data[i,1]
  #get midpoint of each column
  coord <- MTT3_map@polygons[[i]]@labpt
  #add name of the region to the coordinates
  coords[[name]] <- coord
}
centroids <- data.frame(do.call(rbind,coords))
centroids$group <- rownames(centroids)



plot_bg <- list()

for(i in seq_along(mo)){
  date <- mo[i]
  
  #Join the map data and proportion unvacciated data by the id (row name)
  p <- ggplot()+
    geom_polygon(data=model_mapf, aes(y=lat, x=long, group=group, fill=Regions), alpha=0.2,color="grey80", show.legend = F)+
    scale_fill_manual(values = getPalette)+
    coord_map()+
    scale_x_continuous()+  
    geom_path(data=MTT3_fort, aes(y=lat,x=long,group=group), size=.7, colour="white")+
    geom_path(data=adm0AfPkf, aes(y=lat, x=long, group=group), size=.7, colour="black")+
    coord_cartesian(xlim = c(60, 78), ylim = c(22, 40)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin=unit(c(0,0,0,0), "cm"))+
    annotate("text",x=75, y=26, label=format(as.Date(date), "%Y-%m"), size=6 , fontface=2)
  # insets <- ggdraw() +
  #   draw_plot(p + coord_cartesian(xlim = c(60, 80), ylim = c(20, 40))) +
  #   draw_plot(p + coord_cartesian(xlim = c(30, 36), ylim = c(-18, -13)), x= 0.05, y= 0.05)
  
  plot_bg[[i]] <- p
  names(plot_bg)[[i]] <- paste(mo[i])
  
}


# FOREGROUND LINE PLOTTING ------------------------------------------------

linedata <- as.data.frame(tree@data)
linedata <- linedata[!is.na(linedata$ancestor_loc),]
linedata$branching <- as.Date(linedata$branching)
frames <- seq(from=min(as.Date(linedata$branching)), to=max(as.Date(mrsd)),by='days')
linedata$drawing <- 0

foreach(j=1:length(frames), .packages=c("ggplot2", "dplyr", "plyr", "lubridate", "RColorBrewer")) %do% {

  linedata[linedata$branching <= frames[j],]$drawing <- difftime(frames[j] ,linedata[linedata$branching <= frames[j],]$branching , units = c("days"))+1
  #counts over 30 do not contribute to the plots
  current_lines <- linedata[linedata$drawing > 0 & linedata$drawing <= 30,]
  if(nrow(current_lines)>0){
    current_lines$ref <- paste(current_lines$ancestor_loc, ".", current_lines$ancestor_loc, "_", current_lines$state, sep="")
    ## look up the line in the list of lines and reference drawing in the number of points of the line to draw
    current_lines <- current_lines[order(-current_lines$drawing),]
    rownames(current_lines) <- NULL
  pltlsts <- list()
  for (k in 1:nrow(current_lines)){
  if(current_lines[k,]$drawing <=10){
  line <- out[[current_lines[k,]$ref]][1:current_lines[k,]$drawing,]
  } else if (current_lines[k,]$drawing >10 & current_lines[k,]$drawing <=20) {
  line <- out[[current_lines[k,]$ref]][(current_lines[k,]$drawing-10):current_lines[k,]$drawing,]
  } else {
    line <- out[[current_lines[k,]$ref]][(current_lines[k,]$drawing-10):20,]
  }
  line$ID2 <- as.numeric(k)
  pltlsts[[k]] <- line
  }
  pltlst <- do.call(rbind, pltlsts)
  pltlst$size <- 1.3
  bg_list <- pltlst
  bg_list$START <- NA
  bg_list$size <- 1.5

  bg_list$ID2 <- bg_list$ID2-0.00001
  tst <- rbind(pltlst, bg_list)
  tst <- tst[order(tst$ID2),]
  tst$ID2 <- ordered(factor(tst$ID2))
  tst$size <- factor(tst$size)
  } else {pltlst <- data.frame()}

  bg_month <- as.Date(floor_date(frames[j], "month"))
  if(nrow(pltlst)>0){
   plot <-  plot_bg[[paste(bg_month)]]+
      geom_path(data=tst, aes(y=X2, x=X1, group=ID2, colour=factor(START), size=size), lineend="round", show.legend = F)+
      geom_point(data = centroids, aes(y=X2, x=X1, group = group, fill = group), shape = 21, color = "white", size = 5, show.legend = F) +
      scale_color_manual(values = newpal, drop=FALSE, na.value="white")+
      scale_size_manual(values=c(1.3,2))

  } else {
    plot <- plot_bg[[paste(bg_month)]]+
      geom_point(data = centroids, aes(y=X2, x=X1, group = group, fill = group), shape = 21, color = "white", size = 5, show.legend = F)
  }
  # plots will be daily
  ggsave(paste0("PkAf_Anim_", j, ".png"), plot = print(plot), height = 6.78, width = 6, units = "in", device="png", dpi=100, path = paste0(directory_interest, "anim/"))
  gc()
  rm(plot)
}




# Plotting maps and circos of supported movements from SpreaD3 ------------

## Arrow map of supported movements
# manually edit the txt from spread3 to show only supported results and add in the values (probability,numbers of movements and HPDs) from tracer
transitions <- read.csv(paste0(directory_interest, "BF_support.csv"))
colnames(transitions) <- c("From", "To", "Bayes_Factor", "Posterior_Probability", "Median_Transitions", "Lower_HPD", "Upper_HPD")

## subset curves to this list
transitions <- left_join(transitions, locationmatrix, by = c("From" = "ancestor_loc", "To" = "state"))

# make one data frame
out2 <- do.call("rbind", out)
out2 <- left_join(transitions[, c("ID", "Median_Transitions")], out2 , by="ID")


# Single background plot fro the mapping
bg_blank <- ggplot()+
  geom_polygon(data=model_mapf, aes(y=lat, x=long, group=group, fill=Regions), alpha=0.2,color="grey80", linewidth=.7, show.legend = F)+
  scale_fill_manual(values = getPalette)+
  coord_map()+
  scale_x_continuous()+  
  geom_path(data=MTT3_fort, aes(y=lat,x=long,group=group), linewidth=.7, colour="white")+
  geom_path(data=adm0AfPkf, aes(y=lat, x=long, group=group), linewidth=.7, colour="black")+
  coord_cartesian(xlim = c(60, 78), ylim = c(22, 40)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"))



arrows <- 
  bg_blank + 
  geom_path(data=out2[out2$Median_Transitions>0,], aes(y=X2, x=X1, group=ID,size=Median_Transitions, color = factor(START)),lineend="butt", linejoin="mitre")+
  geom_point(data = centroids, aes(y=X2, x=X1, group = group, fill = group), shape = 21, color = "white", size = 5, show.legend = F) +
  scale_size(range=c(1,4))+
  scale_color_manual(values = getPalette)+
  guides(color="none", size=guide_legend(title = "Number of\nTransitions", override.aes = list(alpha = 1)))+
  theme(legend.position=c(0.84,0.2),
        legend.background = element_rect(colour="black"),
        legend.key=element_blank())

  ggsave(paste0(directory_interest, "Supported_mvmt_",latestdate,".png"), width = 6, height = 7)
  ggsave(paste0(directory_interest, "Supported_mvmt_",latestdate,".svg"), width = 6, height = 7)
  
# circos plot of same info
require(circlize)
 
png(paste0(directory_interest, "Supported_mvmt_circ_",latestdate,".png"), units="in", width=5, height=5, res=350)
chordDiagram(transitions[,c("From", "To", "Median_Transitions")], grid.col = newpal, transparency = 0.6)
dev.off()
circos.clear()

svg(paste0(directory_interest, "Supported_mvmt_circ.svg"), width=5, height=5)
chordDiagram(transitions[,c("From", "To", "Median_Transitions")], grid.col = newpal, transparency = 0.6)
dev.off()
circos.clear()

# Table of supported values
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"

require(formattable)
plottab <- subset(transitions, select = -c(ID,Lower_HPD, Upper_HPD))
colnames(plottab) <- c("From", "To", "Bayes Factor Support", "Posterior Probability", "Median Transitions")
plottab <- plottab[,c("From", "To", "Bayes Factor Support", "Posterior Probability", "Median Transitions")]
plottab$`Bayes Factor Support` <- round(plottab$`Bayes Factor Support`, 2)
plottab$`Posterior Probability` <- round(plottab$`Posterior Probability`, 2)
plottab$`(95% HPD)` = paste0(" (", transitions$Lower_HPD, "-", transitions$Upper_HPD, ")")
formattable(plottab,
            align = c("l", "l", "c", "c", "c", "c"),
            list(
              `Median Transitions` = color_tile(customGreen0, customGreen)
            ))



  
## tree subplotting in R
#need to plot tree in full first as the subtree functions just zoom the original plot
tr <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none")+
  geom_rootpoint()


tr2 <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  scale_y_reverse()+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_rootpoint()

# extract subtree at any node where the ancestor location and current location do not match 
# collapse any subtree within this where there is a change to a different location and label this with a symbol in the colour of the location it travels to

## FIND NODES WHERE THE ANCESTOR DOES NOT MATCH THE CURRENT LOCATION
treelvls <- tree@data[order(tree@data$branching),]#
treelvls <- unique(treelvls$state)
subtree <- tree@data[tree@data$state!=tree@data$ancestor_loc & !is.na(tree@data$state)& !is.na(tree@data$ancestor),]
subtree$state <- factor(subtree$state, levels=treelvls)
subtree <- subtree[order(subtree$state, subtree$branching),]## SORT THIS BY LOCATION SO IT IS ALREADY DONE BEFORE PLOTTING TO LIST
subtree$state <- as.character(subtree$state)
\

off <- function(offspring){
  if(tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$parent,]$colour &
     tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$colour!=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$parent,]$colour
     ########################################## THIS NEEDS TO CHECK FURTHER BACK IN TIME ############################################
  ){  
    return(offspring) 
  }
}

'%!in%' <- Negate('%in%')
