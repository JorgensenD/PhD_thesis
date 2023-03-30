
# read in nexus format alongside BEAST format previously loaded
tre <- read.nexus(paste0(directory_interest, "MJ_AFP_ca.trees"))

# Drop tip function to collapse onward transmissions --------
custom.drop.tip <- function (phy, tip, tree, trim.internal = TRUE, subtree = TRUE, root.edge = 0, 
                             rooted = is.rooted(phy), collapse.singles = FALSE, interactive = FALSE) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)  # total tips
  if (interactive) {
    cat("Left-click close to the tips you want to drop; right-click when finished...\n")
    xy <- locator()
    nToDrop <- length(xy$x)
    tip <- integer(nToDrop)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    for (i in 1:nToDrop) {
      d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
      tip[i] <- which.min(d)
    }
  }
  else {
    if (is.character(tip)) 
      tip <- which(phy$tip.label %in% tip) # convert to tip numbers
  }
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) 
    return(phy)
  if (length(tip) == Ntip) {
    if (Nnode(phy) < 3 || trim.internal) {
      warning("drop all tips of the tree: returning NULL")
      return(NULL)
    }
  }
  og_tree <- phy
  wbl <- !is.null(phy$edge.length)
  if (length(tip) == Ntip - 1 && trim.internal) { ## only works if there is only one remaining
    i <- which(phy$edge[, 2] == (1:Ntip)[-tip])   ## drop edges ending with the dropped tips
    res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i, 
                                                                             2]], Nnode = 1L)
    class(res) <- "phylo"
    if (wbl) 
      res$edge.length <- phy$edge.length[i]
    if (!is.null(phy$node.label)) 
      res$node.label <- phy$node.label[phy$edge[i, 1] - 
                                         Ntip]
    return(res)
  }
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  if (subtree) {
    trim.internal <- TRUE
    tr <- reorder(phy, "postorder")
    N <- .C(node_depth, as.integer(Ntip), as.integer(tr$edge[, 
                                                             1]), as.integer(tr$edge[, 2]), as.integer(Nedge), 
            double(Ntip + Nnode), 1L)[[5]]
  }
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  if (trim.internal) {
    ints <- edge2 > Ntip
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep
      if (!sum(sel)) 
        break
      keep[sel] <- FALSE
    }
    if (subtree) {
      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
      keep[subt] <- TRUE
    }
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          degree <- tabulate(edge1[keep])
          if (degree[NEWROOT] > 1) 
            break
        }
        keep[j] <- FALSE
        if (length(j) > root.edge) 
          j <- 1:root.edge
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge)) 
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  if (!root.edge) 
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  if (subtree) {
    i <- which(tip %in% oldNo.ofNewTips)
    if (length(i)) {
      # phy$tip.label[tip[i]] <- "[1_tip]"
      phy$tip.label[tip[i]] <- paste("collapse",sapply( strsplit( phy$tip.label[tip[i]], '\\_' ), function(x){ tail(x,1)[1] }), sep = "_")
      tip <- tip[-i]
    }
  }
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  if (length(tip)) 
    phy$tip.label <- phy$tip.label[-tip]
  if (subtree || !trim.internal) {
    node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
    length(node2tip)
    new.tip.label <- if (!length(node2tip)) {
      character(0)
    } else if (subtree) {
      paste("collapse", lapply(node2tip, function(x) tree@data[tree@data$node==MRCA(tree, og_tree$tip.label[offspring(og_tree, x, type = "tips")]),]$state), sep = "_")
    } else {
      if (is.null(phy$node.label)) 
        rep("NA", length(node2tip))
      else phy$node.label[node2tip - Ntip]
    }
    phy$tip.label <- c(phy$tip.label, new.tip.label)
  }
  if(any(grepl("collapse", phy$tip.label))){
    # half any collapsed branches
    edgenos <- which(phy$edge[,2] %in% which(grepl("collapse", phy$tip.label)))
    phy$edge.length[edgenos] <- phy$edge.length[edgenos]/2
  }
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[which(newNb > 0) - 
                                       Ntip]
  if (collapse.singles) 
    phy <- collapse.singles(phy)
  phy
}




# Year bands for tree plots -----------------------------------------------

## stripes for years rather than seasons
low_season <- matrix(c("2006-01-01","2007-01-01","2008-01-01","2009-01-01","2010-01-01","2011-01-01","2012-01-01","2013-01-01","2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01","2019-01-01", "2020-01-01", "2021-01-01", "2022-01-01", "2023-01-01"),ncol=2,byrow=TRUE)

colnames(low_season) <- c("Start","End")
low_season <- as.data.frame(low_season)
low_season$Start <- as.Date(low_season$Start)
low_season$End <- as.Date(low_season$End)




require(phylobase)
require(stringr)
require(tidytree)
require(gginnards)
tstcl <- ggplot_build(tr) ## need to build to get the coordinates of each node
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")
treeroot <- rootnode(tre)
rm(nodenum)
subtreelist <- list()
treedata <- list()



# Tree splitting function - could be optimised and cleaned up -------------


subtrefunc <- function(i,subtree, tree){ 
  nodenum <<- as.numeric(subtree[i,]$node)
  offspring <- NULL
  xmin <- decimal_date(min(tree@data$branching))
  xmax <- decimal_date(mrsd)
  if(nrow(tstcl[["data"]][[2]][tstcl[["data"]][[2]]$parent==nodenum,])>0){
    #droptips from different locations
    offspring <- tidytree::offspring(tre,subtree[i,]$node)  # find all offspring of the node of interest
    
    ## get tips
    tips <- offspring[offspring<=nTips(tre)]
    parentlist <- list()
    for(k in 1:length(tips)){
      parent <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[k],]$parent
      ancest <- parent
      while(treeroot%!in%parent){
        ancest <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==ancest,]$parent
        parent <- c(parent, ancest)
      }
      parentlist[[k]] <- parent
    }
    tocollapse <- unlist(lapply(as.list(offspring), off))
    
    ## use parentlist instead of working out the tips
    
    collapsetips <- list()
    for(j in 1:length(tips)){
      if(any(parentlist[[j]]%in%tocollapse | tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[j],]$colour!=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[j],]$parent,]$colour)){
        collapsetips[[j]] <- tips[j]  
      }
    } 
    collapsetips <- unlist(collapsetips)
    if(!is.null(collapsetips)){
    collapsenames <- list()
      for(r in 1:length(collapsetips)){
      ## get names for these tips
       collapsenames[[r]] <- tre[["tip.label"]][collapsetips[r]]
      }
    collapsenames <- unlist(collapsenames)
# Custom function above implemented here
    droptree <-  custom.drop.tip(extract.clade(tre,nodenum), tip=collapsenames, tree=tree, trim.internal = T, collapse.singles = F, subtree = T)
    
    }else{
      droptree <- extract.clade(tre,nodenum)
    }
    # blank plot if no tree
    if(is.null(droptree)){
        subtreelist[[i]] <- ggplot()+
        theme_bw() +
          coord_cartesian(xlim=c(xmin-1.5, xmax+.2), ylim=c(0,2))+ #clip off v. important to allow the plots to overlap
          theme(legend.position = "none",
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.background=element_rect(fill = "transparent"),
                panel.border = element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_rect(fill = "transparent"),
                plot.title = element_blank(),
                legend.spacing = unit(0,"cm"),
                plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
        tipshapes <- NA
        tipnames <- NA
        alltips <- 0
        edges <- as.data.frame(tre[["edge"]])
        rootdate <- as.character(as.Date(subtree[i,]$branching))
          #as.character(as.Date(tree@data[tree@data$node==as.numeric(edges[edges[,2]==nodenum,][1]),]$branching))
        tipdateout <- NA
      }else if (nTips(droptree)>1){   ###  ADD A LESS THAN 5 TO THIS ONE
    edges <- as.data.frame(tre[["edge"]])
    ## can we get this from the branching date somehow? 
    rootdate <- decimal_date(as.Date(subtree[i,]$branching))
    current_root <- decimal_date(mrsd) - as.numeric(tree@data[tree@data$node==nodenum,]$height)
    droptree$root.edge <- current_root - rootdate
    rootdate <- as.character(date_decimal(rootdate))
    ## extract mrsd from the tip dates
    tipnames <- droptree$tip.label
    tipdates <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
    tipdate <- max(tipdates, na.rm = T)

    tipshapes <-     sapply( strsplit( tipnames, '\\_' ), function(x){  head(x,2)[2]})
     tipshapes[which(tipshapes %!in% c("AFP", "ES"))] <- NA
    ## extract AFP and ES and set all the other tips to NA
      
      length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
    alltips <- length(tipshapes[!is.na(tipshapes)])
    subtreelist[[i]] <- ggtree(droptree, mrsd=tipdate,colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour) +
      scale_y_reverse()+
      coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ #clip off v. important to allow the plots to overlap
      ###ADD YLIMS A BIT WIDER THAN THE TREEE
      geom_tippoint(aes(subset = !grepl("collapse", label), shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2) +
      scale_shape_manual(values=tipsh)+
      geom_tiplab(aes(subset = grepl("collapse", label),label = "\u2B9E", color = c(sapply( strsplit( droptree$tip.label, '\\_' ), function(x){ tail(x,1)[1] }), rep("NA", Nnode2(droptree)-Ntip(droptree)))), show.legend = F, offset = -.02, vjust = 0.4) +
      scale_color_manual(values = getPalette, drop=FALSE, na.value="white") +
      
      geom_rootedge(colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
      theme_tree2() +
      geom_rootpoint(shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, position = position_nudge(x=-droptree$root.edge), size=2)+
      theme_bw()+
      theme(legend.position = "none",
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_rect(fill = "transparent"),
            panel.border = element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_rect(fill = "transparent"),
            plot.title = element_blank(),
            legend.spacing = unit(0,"cm"),
            plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
    tipdateout <- as.character(tipdate)
    
    }else{
      tipnames <- droptree[["tip.label"]][1]
      # extract the date from this string
      tipdate <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
      tipshapes <- read.table(text = as.character(tipnames), sep = "_")$V2
      length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
      #tipdate <- as.POSIXct(max(tipdates))
      edges <- as.data.frame(tre[["edge"]])
      rootdate <- as.character(as.Date(subtree[i,]$branching))
      
      alltips <- 1
      subtreelist[[i]] <- ggplot() +
        geom_line(aes(x=c(decimal_date(as.Date(rootdate)),decimal_date(tipdate)), y=c(1,1)), colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
        geom_point(aes(x=decimal_date(as.Date(rootdate)),y=1),shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, size=2)+
        geom_point(aes(x=decimal_date(tipdate),y=1, shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2)+
        scale_shape_manual(values=tipsh)+
        coord_cartesian(xlim= c(xmin-1.5, xmax+.2), clip = 'off')+
        theme_bw()+
        theme(legend.position = "none",
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill = "transparent"),
              panel.border = element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_rect(fill = "transparent"),
              plot.title = element_blank(),
              legend.spacing = unit(0,"cm"),
              plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
      tipdateout <- as.character(tipdate)
    }
    #nudge the point by the length of the root edge
  }else{
    # if a single tip is left need to plot it as a point on an empty ggplot axis then add a line left of it the length of the distance to the mrca and then another point in the colour of the mrca
    # therefore need the date of the node from it's name
    nodename <- tre[["tip.label"]][nodenum]
    tipnames <- nodename
    # extract the date from this string
    tipdate <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
    tipshapes <- read.table(text = as.character(nodename), sep = "_")$V2
    tipshapes <- tipshapes[!is.na(tipshapes)]
    #tipdate <- as.POSIXct(max(tipdates))
    
    edges <- as.data.frame(tre[["edge"]])
    ## can we get this from the branching date somehow? 
    rootdate <- decimal_date(as.Date(subtree[i,]$branching))
    current_root <- decimal_date(mrsd) - as.numeric(tree@data[tree@data$node==nodenum,]$height)
    root.edge <- current_root - rootdate
    rootdate <- as.character(date_decimal(rootdate))
    
    alltips <- 1
    subtreelist[[i]] <- ggplot() +
      geom_line(aes(x=c(decimal_date(as.Date(rootdate)),decimal_date(tipdate)), y=c(1,1)), colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
      geom_point(aes(x=decimal_date(as.Date(rootdate)),y=1), shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, size=2)+
      geom_point(aes(x=decimal_date(tipdate),y=1, shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2)+
      scale_shape_manual(values=tipsh)+
      coord_cartesian(xlim= c(xmin-1.5, xmax+.2), clip = 'off')+
      theme_bw()+
      theme(legend.position = "none",
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_rect(fill = "transparent"),
            panel.border = element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_rect(fill = "transparent"),
            plot.title = element_blank(),
            legend.spacing = unit(0,"cm"),
            plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
    
      tipdateout <- as.character(tipdate)
    #need to set the x axis based on tr
    #from the lowest date in tree@data$branching to the MRSD
  }

  subtreelist[[i]] <- print(subtreelist[[i]])
  ancestor_loc <-  tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour
  location <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour
  tipshap <- as.character(table(factor(tipshapes, levels=c("AFP", "ENV"))))
  out <- c(ancestor_loc,location,rootdate,tipdateout,alltips,tipshap, toString(tipnames))
  print(i)
  
  return(list(subtreelist[[i]],out))
  }


#apply this fuction
treedat <- lapply(1:nrow(subtree),subtrefunc, subtree=subtree, tree=tree)




## split the lists, first 9 elements to the ggplot and the rest to downstreamtree data frame
subtreelist <- list()
treedata <- list()
for(i in 1:length(treedat)){
  subtreelist[[i]] <- append_layers(treedat[[i]][[1]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")
  class(subtreelist[[i]]) <- class(treedat[[i]][[1]])
  treedata[[i]] <- treedat[[i]][[2]]
  }

downstreamtree <- data.frame(treedata)
downstreamtree <- data.table::transpose(downstreamtree)


# Root LTL ----------------------------------------------------------------
# The root subtree cannot be edited inside of the function as implemented so must be generated separately
# Need to respecifiy some parameter which are set in the full function 

##get root subtree
xmin <- decimal_date(min(tree@data$branching))
xmax <- decimal_date(mrsd)
## parent list of any node in the tree gives the root as the last entry
treroot <- last(ancestors(phylo4(tre), 1))
nodenum <- treroot
tiplist <- tidytree::offspring(tre,treroot)
# find all nodes to collapse
tocollapse <- unlist(lapply(as.list(tiplist), off))
# find all of the tips which are ancestors of these
collapse <- as.vector(unlist(descendants(phylo4(tre), tocollapse, type = c("tips"))))
droptree <- custom.drop.tip(tre, tip=collapse, tree=tree)
tipnames <- droptree$tip.label
tipdates <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
tipdate <- max(tipdates, na.rm = T)

tipshapes <-     sapply( strsplit( tipnames, '\\_' ), function(x){  head(x,2)[2]})
tipshapes[which(tipshapes %!in% c("AFP", "ES"))] <- NA

edges <- as.data.frame(tre[["edge"]])
rootlng <- 1
droptree$root.edge <- tre[["edge.length"]][rootlng]

length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
roottips <- length(tipshapes[!is.na(tipshapes)])

rootdate <- as.character(as.Date(tree@data[tree@data$node==treroot,]$branching))

treelist <- list()
treelist[[1]] <- ggplot()
treelist[[2]] <- ggtree(droptree, mrsd=tipdate,colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour) +
  scale_y_reverse()+
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ #clip off v. important to allow the plots to overlap
  ###ADD YLIMS A BIT WIDER THAN THE TREEE
  geom_tippoint(aes(subset = !grepl("collapse", label), shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2) +
  scale_shape_manual(values=tipsh)+
  geom_tiplab(aes(subset = grepl("collapse", label),label = "\u2B9E", color = c(sapply( strsplit( droptree$tip.label, '\\_' ), function(x){ tail(x,1)[1] }), rep("NA", Nnode2(droptree)-Ntip(droptree)))), show.legend = F, offset = -.02, vjust = 0.4) +
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white") +
  
  geom_rootedge(colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
  theme_tree2() +
  geom_rootpoint(shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, position = position_nudge(x=-droptree$root.edge), size=2)+
  theme_bw()+
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill = "transparent"),
        plot.title = element_blank(),
        legend.spacing = unit(0,"cm"),
        plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
treelist[[2]] <- append_layers(treelist[[2]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

## Add the root tree in with the others
treelist <- append(treelist, subtreelist)


# Generate axis for each subtree plot -------------------------------------

treaxis <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ 
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none")+
  geom_rootpoint()

# extract only the axis
gTable <- ggplot_gtable(ggplot_build(treaxis))
grid.newpage()
groblst <- list()
groblst[[1]] <- as_ggplot(gTable$grobs[[7]])
treetst <- append(treelist, groblst)


# Set tipshapes for two data sources --------------------------------------

tipls <-  as.character(table(factor(tipshapes, levels=c("AFP", "ES"))))


#add root
downstreamtree[nrow(downstreamtree)+1,] <- c(NA, tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==treroot,]$colour, rootdate ,max(tipdates),nTips(droptree),tipls, toString(tipnames))


# Extract locations from the branch colours -------------------------------

# convert the colours in the list in to locations
locatefunc <- function(x){
  if(!is.na(x)){
return(names(getPalette[getPalette==x]))}
  else(return(NA))}
downstreamtree[,1] <- unlist(lapply(downstreamtree[,1], locatefunc)) # for some reason wont let me apply to both columns together, maybe thinks i am passing two elements to the function?
downstreamtree[,2] <- unlist(lapply(downstreamtree[,2], locatefunc))


# Unused plots of imports/exports etc. ------------------------------------

#number of importations per area
imp <- as.data.frame(table(downstreamtree$V2))

sum1 <- ggplot(downstreamtree, aes(factor(downstreamtree$V2), fill=downstreamtree$V2)) +
  geom_bar(stat="count", position = "dodge") + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        #axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  xlab("Region")+
  ylab(paste("Number of importations"))


#number of exportations per region
downstreamtree$V1 <- factor(downstreamtree$V1, levels=sort(unique(downstreamtree$V2)))
dta <- as.data.frame(table(downstreamtree$V1))
sum5 <- ggplot(dta, aes(x=dta$Var1,y=dta$Freq, fill=dta$Var1)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        #axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  xlab("Region")+
  ylab(paste("Number of exportations"))

#proportion of these importations which come from each of the other locations
sum2 <- ggplot(subset(downstreamtree, !is.na(downstreamtree$V1)), aes(factor(V2), fill=V1)) +
  geom_bar(stat="count", position = "fill", na.rm = TRUE) + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title.x=element_blank())+
  guides(fill = guide_legend(title="Importation from",title.position="top", title.hjust = 0.5))+
  xlab("Region")+
  ylab(paste("Proportion of importations \n from each region")
)

#stripchart of the survival time of lineages in each area
downstreamtree$V3 <- as.Date(downstreamtree$V3)
downstreamtree$V4 <- as.Date(downstreamtree$V4)
downstreamtree$diff <- as.numeric(downstreamtree$V4-downstreamtree$V3)
#root
rootlin <- as.numeric(as.Date(max(tipdates))-as.Date(min(tree@data$branching)))

sum3 <- ggplot(downstreamtree, aes(x=diff)) + 
  geom_histogram(aes(fill=V2), binwidth = 20)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  scale_x_continuous(expand = c(0,0))+
  #annotation_logticks(sides="b")+
  guides(fill = guide_legend(title="Region",title.position="top"))+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Lineage survival (days)")

#similar representation of the number of points per tree in each area
downstreamtree$V5 <- as.numeric(downstreamtree$V5)

sum4 <- ggplot(downstreamtree, aes(x=V5)) + 
  geom_histogram(aes(fill=V2), binwidth = 2)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  #scale_y_log10()+
  #annotation_logticks(sides="l")+
  guides(fill = guide_legend(title="Region",title.position="top", title.hjust = 0.5, nrow = 2))+
  theme(legend.position = "bottom",
        legend.justification="center",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Number of tips")
  #annotate("text", x=1,y=120,label=paste("root subtree tips:", roottips, sep=" ") )
downstreamtree$V6 <- as.numeric(downstreamtree$V6)


# AFP only from the full tree
AFP_tips <- ggplot(downstreamtree, aes(x=V6)) + 
  geom_histogram(aes(fill=V2), binwidth = 2)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  #scale_y_log10()+
  #annotation_logticks(sides="l")+
  guides(fill = guide_legend(title="Region",title.position="top", title.hjust = 0.5, nrow = 2))+
  theme(legend.position = "bottom",
        legend.justification="center",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Number of AFP tips")

grid_arrange_shared_legend(sum4, AFP_tips, sum3, arrangeGrob( sum1,sum5,sum2, ncol=1, nrow=3, heights = c(1/4, 1/4, 1/2)), ncol=4, nrow=1, position = "bottom")




polioclusters <- downstreamtree
names(polioclusters) <- c("root_location", "location", "root_date", "max_tipdate", "number_tips","AFP","ENV", "cluster_survival") 
# save to file
#save(polioclusters, file="C:/Users/dnj13/Documents/DTA_10_cluster(corr)/cluster_detail.Rdata")


# Simple movement matrix --------------------------------------------------

require(reshape2)
polio_matrix <- polioclusters[!is.na(polioclusters$root_location),]
polio_matrix <- table(polio_matrix$root_location, polio_matrix$location)
diag(polio_matrix) <- NA
melted_pm <- melt(polio_matrix)


matrix <- ggplot(data = melted_pm, aes(x=Var2, y=Var1)) + 
  geom_tile(aes(fill=value))+
  geom_text(aes(label = value)) +
  theme( legend.position = "none",
         axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
  scale_fill_gradient(low = "white", high = "red")+
  ylab("location exported from")+
  xlab("location imported to")
matrix
## save the datasets for pierre to use
#save(polioclusters, file = "C:/Users/dnj13/Documents/Pak_Afg_Phylo/DTA/DTA_10_PLOTS/polioclusters.RData")
#write.csv(polioclusters, file = "C:/Users/dnj13/Documents/Pak_Afg_Phylo/DTA/DTA_10_PLOTS/polioclusters.csv",row.names=FALSE)


downstreamtree$ID <- 2:(nrow(downstreamtree)+1)
downstreamtree[downstreamtree$ID==nrow(downstreamtree)+1,]$ID <- 1
downstreamtree <- downstreamtree[order(downstreamtree$ID),]
rownames(downstreamtree) <- NULL

require(egg)
tiplist <- as.list(as.numeric(downstreamtree$V5))

alltips <- unlist(tiplist)

##rescale with tanh to ensure the plots are not dominated by large or small subtrees but there is still some scaling with size.
alltips <- 1-2/((exp(1)^(.1*alltips)) +1)
alltips[length(alltips)+1] <- 0.4


spacing <- alltips/sum(alltips)

clusternames <- c("East Pakistan"="EAST-PAK", "South Corridor"="SOUTH-CORRIDOR", "North Corridor"="NORTH-CORRIDOR", "West Afghanistan"="WEST-AFG", "Central Pakistan"="CENTRE-PAK", "Central Corridor"="CENTRAL-CORRIDOR", 
            "Karachi"="KARACHI", "East Afghanistan"="EAST-AFG", "West Pakistan"="WEST-PAK", "North Afghanistan"="NORTH-AFG", "Other"="OTHER")
cluster_shrt <- unique(downstreamtree[downstreamtree$V2!="NORTH-CORRIDOR+EAST-PAK",]$V2)




# Plot cluster blocks -----------------------------------------------------

plttrees <- treetst
plttrees[[1]] <- ggplot() + 
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), ylim=c(0,2), clip = 'off')+ 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill = "transparent"),
        plot.title = element_blank(),
        legend.spacing = unit(0,"cm"),
        plot.margin=unit(c(-0.15,0,-0.15,0), "lines"))
plttrees[[1]] <- append_layers(plttrees[[1]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

cluster_block <- list()
for (i in 1:length(cluster_shrt)){
plttrees_tips <- c(0.4,alltips)
plttrees_spacing <- plttrees_tips[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(plttrees_tips))]#/sum(plttrees_tips[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(plttrees_tips))])
arrange <- ggarrange(plots = plttrees[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(treetst))], heights = plttrees_spacing, ncol=1,padding = unit(0, "line"))
cluster_block[[i]] <- annotate_figure(arrange, top = grobTree( rectGrob(gp=gpar(fill="grey"), height = unit(3, "npc")), textGrob(names(which(clusternames==cluster_shrt[i])), gp=gpar(fontsize=14, col="black"), vjust = 0.8)))
}
## 4npc

tiff("cluster_blocks.tiff",width = 12, height = 60, units = 'in', res = 100)
grid.arrange(cluster_block[[1]],cluster_block[[2]],cluster_block[[3]],cluster_block[[4]],cluster_block[[5]],cluster_block[[6]],cluster_block[[7]],cluster_block[[8]],cluster_block[[9]],cluster_block[[10]],cluster_block[[11]], nrow=11, ncol=1, heights=c(14,11,10,3,12,11,12,1,1,1,2))
dev.off()

png("cluster_examples.png",width = 12, height = 14, units = 'in', res = 200)
grid.arrange(cluster_block[[1]],cluster_block[[2]],cluster_block[[6]], nrow=3, ncol=1, heights=c(12,9,10))
dev.off()


png("cluster_block1.png" ,width = 12, height = 14   , units = 'in', res = 200)
cluster_block[[1]]
dev.off()
png("cluster_block2.png" ,width = 12, height = 11   , units = 'in', res = 200)
cluster_block[[2]]
dev.off()
png("cluster_block3.png" ,width = 12, height = 10, units = 'in', res = 200)
cluster_block[[3]]
dev.off()
png("cluster_block4.png" ,width = 12, height = 3, units = 'in', res = 200)
cluster_block[[4]]
dev.off()
png("cluster_block5.png" ,width = 12, height = 12, units = 'in', res = 200)
cluster_block[[5]]
dev.off()
png("cluster_block6.png" ,width = 12, height = 11, units = 'in', res = 200)
cluster_block[[6]]
dev.off()
png("cluster_block7.png" ,width = 12, height = 12 , units = 'in', res = 200)
cluster_block[[7]]
dev.off()
png("cluster_block8.png" ,width = 12, height = 0.9 , units = 'in', res = 200)
cluster_block[[8]]
dev.off()
png("cluster_block9.png" ,width = 12, height = 0.9 , units = 'in', res = 200)
cluster_block[[9]]
dev.off()
png("cluster_block10.png",width = 12, height = 0.9 , units = 'in', res = 200)
cluster_block[[10]]
dev.off()
png("cluster_block11.png",width = 12, height = 0.9 , units = 'in', res = 200)
cluster_block[[11]]
dev.off()


# Full tree with tip shapes -----------------------------------------------

tipshapes_n <- read.table(text = as.character(tre$tip.label), sep = "_")$V2
length(tipshapes_n) <- length(tre$tip.label)+tre$Nnode
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")

tr2 <- ggtree(tree, mrsd=mrsd, aes(colour=state), lineend = "square") +
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions")+
  geom_tippoint(aes( fill=state, shape = tipshapes_n), colour="black")+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions")+
  scale_shape_manual(values = tipsh, name = "Tip shape")+
  coord_cartesian(xlim= c(xmin, xmax+.2), clip = 'off')+
  geom_vline(xintercept = decimal_date(earliest_date), size = 1, color = "darkorange")+
  theme_tree2() +
  theme(legend.position=c(.15, .45),
        axis.text.x=element_text(size=11),
        legend.key.width=unit(.6, "cm"),
        legend.spacing.x = unit(.3, 'cm'),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", 
                                         colour ="black"),
        legend.margin= margin(6,6,6,10) ,
        plot.margin=unit(c(0.1,0.1,-0.5,0.4), "cm"))+
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(override.aes = list(shape = NA, linewidth=5)))

# tree without cases
tr2 <- append_layers(tr2, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")


earliest_date = as.Date("2012-01-07")


tr2 

# Need to import case data - not included here - linelist data, columns are yearmonth and region (epiblock)

case_bars <- ggplot() +
  geom_bar(WT1_AFP_region, mapping = aes(x = YM, fill = epiblock), stat = "count", show.legend = F) +
  scale_y_continuous("Reported cases", breaks = seq(0,80, by = 10), expand = c(0,0)) +
  theme_cowplot() +
  coord_cartesian(xlim= c(as.Date(date_decimal(xmin)), as.Date(date_decimal(xmax+.2)))) +
  geom_vline(xintercept = earliest_date, size = 1, color = "darkorange") +
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white") +
  labs(x = "Date") +
  theme(plot.margin=unit(c(-0.2,0.1,0.1,0.4), "cm"))
case_bars <- append_layers(case_bars, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

# remove some themeing
tr3 <- tr2 +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line())

# plot together with plot grid
png("./PkAf_20220902/2022_tree_cases.png",width = 8, height = 8.5 , units = 'in', res = 200)
plot_grid(tr3, case_bars, ncol = 1, align = "v", axis = "b",rel_heights = c(8,2))
dev.off()

tree_case <- plot_grid(tr3, case_bars, ncol = 1, align = "v", axis = "b",rel_heights = c(8,2))
ggsave("./PkAf_20220902/2022_tree_cases.svg", tree_case, device = "svg", width = 8, height = 8.5, units = "in")


## case facet plots
case_bars <- ggplot() +
  geom_bar(WT1_AFP_region, mapping = aes(x = YM, fill = epiblock), stat = "count", show.legend = F) +
  scale_y_continuous("Reported cases", breaks = seq(0,80, by = 10), expand = c(0,0)) +
  theme_classic()+
  coord_cartesian(xlim= c(as.Date(date_decimal(xmin)), as.Date(date_decimal(xmax+.2))))+
  geom_vline(xintercept = earliest_date, size = 1, color = "darkorange")+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  labs(x = "Date")+
  facet_wrap(~epiblock, ncol = 1)

case_bars <- append_layers(case_bars, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

case_bars

ggsave("case_region.png", height = 9, width = 8)

