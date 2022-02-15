#; ARG Pipeline
#' Filename: run_arg_inference.R
#' Written by: Lauren Castro
#' Created:   August, 2018
#' Edited on: November 4, 2019

#  Description                   ----
#' How the script works
#' Step 0: Source in functions and load needed libraries. Other libraries may be in other scripts
#' Step 1: Read Command Line Inputs and Default Parameters 
#' Step 2: Define paths for experiment
#' Step 3: Set up sampling scheme based on patient of interest 
#' Step 4: Run the ARG simulation
#' Step 5: Get distance matrix of the binary tree for sample of residues 
#' Step 6: Calculate average distance matrix 
#' Step 7: Calculate average tree
#' Step 8: Get metrics and save score 


#### Step 0 ####
rm(list=ls())
source('../scripts/arg_functions.R')
source('../scripts/decompose_functions.R')
source('../scripts/inference_functions.R')

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(CollessLike))



#### Step 1 ####

args         <- commandArgs(trailingOnly = TRUE) 
patient      <- as.numeric(as.character(args[1]))   
recombo.rate <- as.numeric(as.character(args[2]))
num.cores    <- as.numeric(as.character(args[3])) # Used in Step 5

# Default setting if no arguments are supplied 

if(length(args)==0){
 print("No arguments supplied.")
  patient      <- 1 
  recombo.rate <- 0.01
  num.cores   <- 4
} 

patient <- paste0("p", patient)

## Sites to decompose the binary tree -- chosing 100 sites from max of 700 
site.min <- 1 
site.max <- 700 
sites    <-cbind(sample(seq(site.min, site.max), size = 100, replace = FALSE)) 


#### Step 2: Setting the experiment path  ####
if(dir.exists(path = "../output/") == FALSE) 
dir.create(path = "../output/", recursive = TRUE)
experiment.path = "../output/"
        
### Step 3: Get the sampling scheme from the empirical tree ####
rep.times <- 1
tr        <- read.tree(paste0("../data/reduced_tree_",patient, ".tre"))

## Just for this
tips            <- as.numeric(c(sub(".*?\\.(\\w{1})([0-9]{1,2})(\\w{1})(\\d{3}).*", "\\4", tr$tip.label)))
sampling.scheme <- data.frame(table(tips))
sampling.sizes  <- as.numeric(as.character(sampling.scheme$Freq))
sampling.times  <- as.numeric(as.character(sampling.scheme$tips))
final.time      <- sampling.times[length(sampling.times)]*30 # In days 
rm(tr)

## Read in scores so these will be distributed to every core when running in parallel
y           <- readRDS(paste0("../data/", patient, "_tn93_reduced.rds"))
tree.scores <- read.csv(paste0("../data/metrics", patient, "_reduced.csv"))


# Setting up sim parameters 
recombo.rate        <- recombo.rate
selection.strength  <- runif(n=1, min=1, max= 100)
selection.frequency <- runif(n = 1, min = 14, max = 362)
latent.pool         <- runif(n = 1, min = 10^1, max = 10^4)
activate.rate       <- latent.pool*(log(2)/(44*30)) ## based on 44 mo half life
  
if(dir.exists(path = paste0(experiment.path, "args/", patient, "/")) == FALSE) 
  dir.create(path = paste0(experiment.path, "args/", patient, "/"), recursive=TRUE)
  
  
##### Step 4 Runs the ARG and saves ####
#' Input: Simluation Parameters 
#' Output: ARG object (named trial here). ARG objects contains 
#' (1) list that contains one or more arg entries depending on its size. Function get.arg
#' in fast_decomposition.R will collapse the list to return the full arg. 
#' (2) track.lineages: data.frame with summary of num of lineages through time 
#' (3-5) nodes, valid, time: all used for bookkeeping the ARG in between sampling moments 


if(dir.exists(path = paste0(experiment.path, "args/", patient, "/")) == FALSE) 
    dir.create(path = paste0(experiment.path, "args/", patient, "/"), recursive=TRUE)
  
  
trial <- run_trial_sensitivity(activate.rate = activate.rate, 
                                  rep.times = 1, 
                                  experiment.path = experiment.path,
                                  selection.strength = selection.strength,
                                  selection.frequency = selection.frequency,
                                  recombo.rate = recombo.rate,
                                  final.time = final.time,
                                  sample.sizes = rev(sampling.sizes[1:length(sampling.sizes)-1]),
                                  sample.times = sampling.times[1:length(sampling.times)-1],
                                  starting.k = sampling.sizes[length(sampling.sizes)],
                                  latent.pool = latent.pool)

 
file.extension <- paste0(round(selection.strength, digits =3), "_", round(selection.frequency, digits = 3), "_", 
                        recombo.rate, "_", round(activate.rate, digits = 3))


saveRDS(object = trial, file = paste0(experiment.path,"args/", patient, "/", file.extension, ".rds"))

#### Step 5: Get distance matrix of the binary tree of each residue ####
#' Input: ARG object (the actual ARG will be extracted within get.residue.dmat)
#' Output: object dat contains is a list where each element contains the distance matrix of the binary tree
#' at the decomposed site

if(dir.exists(path = paste0(experiment.path, "distance_matrices/", patient, "/")) == FALSE)
  dir.create(path = paste0(experiment.path, "distance_matrices/", patient, "/"), recursive = TRUE)

dm.path <- paste0(experiment.path, "distance_matrices/", patient, "/")
dat     <- mclapply(sites,mc.cores = num.cores,function(x) get.residue.dmat(arg = trial, res.ind = x))
rm(trial)

names(dat) <- sites
file.name  <- paste0(dm.path, "dm_", file.extension, ".RData")
save(x = dat, file = file.name)


#### Step 6: Calculate average distance matrix and make it ready for metrics ####
#' Input: The list of 100 distance matrices (dat)
#' Output: An average distance matrix (avg.dm.df)

avg.dm              <- apply(simplify2array(dat), 1:2, mean)
avg.dm.df           <- data.frame(avg.dm)
time_extension      <- rep(rev(sampling.times), rev(sampling.sizes))
colnames(avg.dm.df) <- paste0(colnames(avg.dm.df), "_", time_extension)
rownames(avg.dm.df) <- colnames(avg.dm.df)
rm(dat)
if(dir.exists(path = paste0(experiment.path, "avg_dm/", patient, "/")) == FALSE)
  dir.create(path = paste0(experiment.path, "avg_dm/", patient, "/"), recursive = TRUE)
write.csv(x = avg.dm.df, file = paste0(experiment.path, "avg_dm/",
                                       patient, "/avgDM_", file.extension, ".csv"), row.names = T)

#### Step 7: Calculate average tree ####
# Input: The average distance matrix (in matrix format)
# Output: Average rooted tree 
rooted.tree <- get_avg_tree_modified(avg.dm = avg.dm, sample.sizes = rev(sampling.sizes))
if(dir.exists(path = paste0(experiment.path, "trees/", patient, "/")) == FALSE)
  dir.create(path = paste0(experiment.path, "trees/", patient, "/"), recursive = TRUE)
write.tree(phy = rooted.tree, file = paste0(experiment.path, "trees/", patient, "/", file.extension, ".tre"))

#### Step 8: Get Metrics and Save Score ####
#' Input: The average distance matrix (in data.frame format)
#' Output: Scores for EI ratio, sackin, mlt, distance matrix similarity, and CV 

distance.results <- find_pairwise_objective(y=y, s = avg.dm.df)
distance.results <- data.frame(matrix(unlist(distance.results), nrow=1, byrow=T),stringsAsFactors=FALSE)
colnames(distance.results)=c("objective", "lambda")

tree.metrics <- get_tree_metrics(tree = rooted.tree, sample.times = sampling.times)
# Combine for score
score           <- data.frame(matrix(nrow = length(sampling.times),ncol = 11))
colnames(score) <- c("args", "sackin", "mlt", "e.i.ratio", "dm", "lambda",
                  "cv", "raw.sackin", "raw.mlt", "raw.e.i.ratio", "raw.cv")

score$args       <- file.extension
score$sackin     <-(tree.metrics$sackin - tree.scores$sackin[1])^2
score$raw.sackin <- tree.metrics$sackin
score$mlt        <- (tree.metrics$mean.lin.time - tree.scores$mean.lin.time[1])^2
score$raw.mlt    <- tree.metrics$mean.lin.time
score$e.i.ratio  <-(tree.metrics$e.i.ratio - tree.scores$e.i.ratio[1])^2
score$raw.e.i.ratio <- tree.metrics$e.i.ratio
score$dm            <-distance.results$objective
score$lambda        <-distance.results$lambda

## find coefficient of variation of simulated
cv              <- find_cv(avg.dm.df)
score$raw.cv    <- cv$CV
colnames(cv)[2] <- "CV_sim"
sim.score       <- tree.scores %>% 
  left_join(cv, by = "time1") %>% 
  mutate(score = abs(CV-CV_sim)) %>%
  summarize(total = sum(score)) %>% 
  pull(total)
score$cv <- sim.score

if(dir.exists(path = paste0(experiment.path, "scores/", patient, "/")) == FALSE)
  dir.create(path = paste0(experiment.path, "scores/", patient, "/"), recursive = TRUE)
write.csv(x=score, file = paste0("../output/scores/", patient,"/",
                                 patient, "_", file.extension, ".csv"), row.names = F)
