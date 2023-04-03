# Data subsets and starting trees - want these in their own folders
require(treedater)
require(dplyr)

# derived from sarscov2Rutils 
make_starting_trees = function(tree , treeoutfn = 'startTrees.nwk', ncpu = 1, ntres = 1){
  tr <- tree
  tr <- di2multi( tr, tol = 1e-5 ) 
  dates <- sapply( strsplit( tr$tip.label, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
  names(dates) <- tr$tip.label
  # resolve polytomies randomly 
  tres <- lapply( 1:ntres, function(i) { 
    tr = unroot( multi2di( tr )  )  
    tr$edge.length <- pmax( 1/29e3/5, tr$edge.length  ) #ensures that edge is >0, makes treedater fit a little nicer 
    tr
  })
  tds <- lapply( tres, function(tr){
    dater( unroot(tr), dates[tr$tip.label], s= 906, omega0 = .01, numStartConditions=0, meanRateLimits = c(0.005,0.015), ncpu = ncpu )
  })
  tds
  outtrees = lapply( tds, function(x){
    class(x) <- 'phylo'
    x
  })
  class( outtrees ) <- 'multiPhylo' 
  write.tree( outtrees 
              , file = treeoutfn 
  )
  invisible( tds )
}

# number of starting trees required
ntres = 8

## All ----
seqs <- read.FASTA(paste0("./Meta_matched/longlabs_seq_AFGPAKMWIMOZ_",latestdate,".fasta"))
meta <- read.csv(paste0("./Meta_matched/longlabs_metadata_AFGPAKMWIMOZ_", latestdate,".csv"))

# Save all
if(file.exists(paste0("./all_",latestdate))){cat("folder exists")} else {dir.create(paste0("./all_",latestdate))}
# use ape write not seqinr
write.FASTA(seqs, paste0("./all_",latestdate,"/longlabs_seq_AFGPAKMWIMOZ_",latestdate,"_all.fasta"))

# Starting trees
fn = paste0("./all_",latestdate,"/longlabs_seq_AFGPAKMWIMOZ_",latestdate,"_all.fasta")

tree <- read.tree(paste0(fn, ".treefile"))
make_starting_trees(tree, treeoutfn = paste0("./all_",latestdate,"/startTrees.nwk"), ncpu =1, ntres = ntres)

## AFP only ----
AFP_seqs <- seqs[grep("_AFP_", names(seqs))]
AFP_seqs_PkAf <- AFP_seqs[grep("_PAKISTAN_|_AFGHANISTAN_", names(AFP_seqs))]

dir.create(paste0("./AFP_PkAf_",latestdate))
if(file.exists(paste0("./AFP_",latestdate))){cat("folder exists")} else {dir.create(paste0("./AFP_",latestdate))}
write.FASTA(AFP_seqs, paste0("./AFP_PkAf_",latestdate,"/longlabs_seq_AFGPAKMWIMOZ_",latestdate,"_AFP.fasta"))

# Starting trees
fn = paste0("./AFP_PkAf_",latestdate,"/longlabs_seq_AFGPAKMWIMOZ_",latestdate,"_AFP.fasta")
system( paste0( 'iqtree -nt AUTO -redo -m HKY -s ', fn ), intern=FALSE)
tree <- read.tree(paste0(fn, ".treefile"))

make_starting_trees(tree, treeoutfn = paste0("./AFP_PkAf_",latestdate,"/startTrees.nwk"), ncpu =1, ntres = 8)
