# PhD_thesis

Presented here are code snippets related to my PhD thesis:  
"Tools and methods for rapid phylogeographic analysis of pathogen genomic surveillance data with application to poliovirus"  
The full thesis is currently not publicly avaialable due to the use of unreleased sequence data.

# Local Transmission Lineages
R code to produce and summarise 'local transmission lineages' and plot a map of movements over time from a discrete location labelled phylogeny generated with BEAST is provided in three files, these should be run in order and depend on outputs of previous files:  
1) ```split_and_map.R```
2) ```LTL.R```
3) ```Cluster_duration_code.R```

The first file generates the maps presented in the thesis together with an animated map over time; The second file produces full sets of local transmission lineages and the combined tree and case figure; The third file genreates the plots of time to sampling and cluster type presented in the thesis.  

# Discrete trait analysis sampling correction
The correction for discrete trait analysis to account for sampling rate differences by location is presented in ```probin.R```. This codes a version of the ```ace``` function from the ```ape``` R package adapted to incorporate a diagonal multiplier to the Q matrix, used to supply sampling rates for each location. The function uses a relative rate setup where the rate from location 1 to 2 is expressed relative to the rate from location 2 to 1. This is carried out to aid with model fitting on the heavily ridged likelihood surface output by likelihood maximisation. The code here uses an eigenvector-based solution, as is standard in the original ```ace``` function. The function can also use matrix multiplication, although this is not reproduced here. This function runs on pre-simulated/inferred phylogenies and a phylogeny must be supplied by the user.

BEAST 1.10.5 XML files implementing the correction at the tips and over the full phylogeny are included. These are labelled ```MJ_tipcorrection.xml``` and ```MJ_wholetreecorrection.xml``` respectively. The sequence data and taxon assignmnents which include the location and date of the sequences have been removed from these files but operators and priors are retained for reference. The tree likelihood and substitution models are of particular importance as this is where the correction is implemented.

# Markov jump DTA
A beast template for Markov Jump DTA with the proposed correction to account for sampling differences between locations is provided. The numbered R code files work forward from attaching metadata to sequence data, selecting a sample to analyse, producing a tree generation script for BEAST and producing a Markov Jump script for BEAST. A set of functions are included to automatically generate sections of the xml from the metadata-attached sequences.

# SARS-CoV-2 model
SARS-CoV-2 analyses presented in the thesis are derived from code avaialble [in an existing repository](https://github.com/JorgensenD/sarscov2Rutils). This includes a simple tutorial to run a Saudi Arabia model.


