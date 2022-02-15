# ARG_simul-decomp

Code from Castro et al "Recombination smooths the time-signal disrupted by latency in within-host HIV phylogenies" for ancestral recombination graph (ARG) simulation and break down into binary trees. ARGs are simulated according to a specific sampling scheme. The code is written to import a sampling scheme provided by one of nine serially sampled HIV patients (see Data).  The ARG is then decomposed into a series of binary trees based on random break points. Using the average distance between tips over the population of decomposed trees, we reconstruct a single distance matrix which is the basis of an average tree. In our application, we compare the average tree and average distance matrix to those of the nine serially sampled HIV patients  [^1].  

## Required R Packages and version used in deveopment 
* CollessLike v1.0
* ape v5.4.1
* igraph v1.2.6
* data.table v1.12.6
* purrr v0.3.4
* magrittr v2.0.1
* tidyr v1.0.0
* string v1.4.0
* tibble v3.1.4
* dplyr v1.0.7
* plyr v1.8.6
* parallel

### Directory Setup
The code assumes the following level 1 directory structure:

```
Project
 |
 +-- data
 |
+-- output
 |
+-- scripts
 |  |  
 |  +-- arg_functions.R
 |  +-- decompose_functions.R
 |  +-- inference_functions.R 
 |  +-- run_arg_inference.R
 |   
 ```
The following output directory structure will be automatically generated with level 3 subfolders for each patient (not shown): 
```
+-- output
 |  | 
 +  |-- args
 +  |-- avg_dm
 +  |-- distance_matrices
 +  |-- scores
 +  |-- trees
 ```

### Data  
This folder contains three file types for each of the nine serially sampled HIV patients that we analyzed:
1. `*csv`: the six calculated metrics used in the scoring algorithm
2. `*rds`: an R object of the distance matrix of the samples
3. `*tre`: a generated binary tree 

## Workflow 
The code and workflow assume an ARG is being simulated, decomposed and scored against a known set of serially sampled HIV data. We have indicated where changes would need to be made for more general use of the ARG simulator. 

Use `run_arg_inference.R`to run the full pipeline composed of:
1. Simulating an ARG based on a specific sampling scheme and input parameters. 
    *The code is set up to import a sampling scheme from an empirical set of serially sampled data. However, a custom sampling scheme can be specified with a series of sampling times (in days) and number of clones sampled at each time.* 
2. Decomposing the ARG into a series of distances matrices based on specific residues
3. Calculating an average distance matrix
4. Calculating an average binary tree
5. Scoring the simulated ARG based on its similarity to longitudinally sampled data
   *If using the ARG simulator to generate non-patient specific ARGs, this step can be commented out.* 

**Input**

The arguments of `run_arg_inference.R` are:
1. patient id: This refers to the integer id of the nine serially sampled HIV patients. *(default = 1)*
2. recombination rate: There are no upper and lower restrictions on this value  *(default = 0.01)*
3. number of cores to use in the decomposition step (step 2 above).  *(default = 4)*

An example of running the pipeline from the command line with arguments is: 
```r
Rscript run_arg_inference.R 2 0.001 6
```

**Output**

As the code is written, in order there are five outputs from `run_arg_inference.R` that will be saved into the level 3 subfolders of `output/`.
1. Data pertaining to the ARG (saved in args/ as an R object). 
    * A list that contains (1) a list of ARG elements (2) a data frame that tracks the number of lineages through time  (3) A bookkeeping number of the number of nodes in the ARG at the time of each sampling event (4) A bookkeeping boolean indicator of whether to advance the ARG (will always be saved in `TRUE`) (5) A bookkeeping entry of the "time" in the simulation (will always be saved as `NULL`)
2. a list of distance matricies corresponding to the binary break down of a specific residue (saved in distance_matrices/ in R Data format)
3. a mean distance matrix (saved in avg_dm/ as a `*csv`)
4. a rooted tree based on the mean distance matrix (saved in trees/ as a `*tre`)
5. a scores compared to the patient's empirical data (saved in scores/ as a `*csv`)

All output filenames follow the structure **patientID_bottleneckStrength_bottleneckFrequency_recombinationRate_activationRate**.  An example is 
`p2_13.63_43.192_0.001_1.452.csv`

#### Description of other scripts  
* `arg_functions.R`: This script contains the functions for simulating the ARG
* `decompose_functions.R`: This script contains the functions for decomposing the ARG into a series of distance matrices based on a specific residue. 
* `inference_functions.R`: This script contains the remaining functions needed for calculating an average distance matrix, average binary tree, and evaluating how well a simulated ARG matches serially sampled HIV data. 


[^1]: Shankarappa RA, Margolick JB, Gange SJ, Rodrigo AG, Upchurch D, Farzadegan H, Gupta P, Rinaldo CR, Learn GH, He XI, Huang XL. Consistent viral evolutionary changes associated with the progression of human immunodeficiency virus type 1 infection. Journal of virology. 1999 Dec 1;73(12):10489-502.


