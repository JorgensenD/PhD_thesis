## xml format functions ----
.seq_format <- function( d ){
  y = '<sequence> 
  <taxon idref="NAME"/>
  SEQUENCE
  </sequence>'
  seqd = lapply(  rownames(d) , function(sid){
    y1 = gsub( pattern = 'NAME', replace=sid, y )
    seq = paste( collapse='', as.character( d[sid, ] )[1,]  ) 
    y2 = gsub( pattern ='SEQUENCE', replace= seq, y1 )
    y2 
  })
  paste( seqd, collapse = '\n' )
}

.taxon_format <- function( d ){
  y = '<taxon id="NAME">
  <date value="DATE" direction="forwards" units="years"/>
		</taxon>'
  seqd = lapply(  rownames(d) , function(sid){
    y1 = gsub( pattern = 'NAME', replace=sid, y )
    date = sapply( strsplit( sid, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
    y2 = gsub( pattern ='DATE', replace= date, y1 )
    y2 
  })
  paste( seqd, collapse = '\n' )
}

## Fill templates ----


require(ape)


# Paste data and trees into premade beast template with the other parameters set
directory_interest <- paste0("./AFP_PkAf_",latestdate, "/")
fasta_name <- "longlabs_seq_AFGPAKMWIMOZ_20220902_AFP.fasta"
xml_name <- paste0(directory_interest, "TREEGEN_TEMPLATE.xml")
## Load template - copy to a new folder
file.copy("TREEGEN_TEMPLATE.xml", xml_name)
file.copy("qsub_anaconda_array_resub_beast1.pbs", paste0(directory_interest, "qsub_anaconda_array_resub_beast1.pbs"))



## load trees
trees <- read.tree(paste0(directory_interest, 'startTrees.nwk'))
ntres <- length(trees)

# format xml tags
d = read.dna(paste0(directory_interest, fasta_name), format = 'fasta')
taxon_date <- .taxon_format(d)
seq_data <- .seq_format(d)


## read in skeleton
x = readLines( xml_name ) 
xmlofn = gsub( xml_name, pattern='TEMPLATE', replacement='FINAL' )

## paste data
for ( k in 1:ntres ){
  xk0 = gsub( x , pattern = 'STARTTREE', replacement = write.tree( trees[[k]] )  )
  xk1 = gsub( xk0, pattern = 'SEQUENCES', replacement = seq_data )
  xk2  = gsub( xk1, pattern='TAXONLIST', replacement= taxon_date )

  if ( !grepl( pattern = '\\.xml$', xmlofn )  )
    writeLines( xk2, con =  paste0( xmlofn, '.', k, '.xml' )  )
  else 
    writeLines( xk2, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.',k,'\\.xml'), xmlofn ) )
}



