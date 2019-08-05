# Jointly inferring the dynamics of population size and sampling intensity from molecular sequences: BESP analyses

### KV Parag, L du Plessis and OG Pybus

## Summary

This repository contains the raw data, configuration files and scripts necessary to reproduce the BESP analyses presented in [https://www.biorxiv.org/content/10.1101/686378v1.abstract](https://www.biorxiv.org/content/10.1101/686378v1.abstract). Some of the scripts may need some adjustment depending on the local setup. 

## Abstract

_Estimating past population dynamics from molecular sequences that have been sampled through time is an important problem in infectious disease epidemiology, molecular ecology and macroevolution. Popular solutions, such as the skyline and skygrid approaches, infer past effective population sizes from the coalescent event times of phylogenies reconstructed from longitudinally sampled sequences. These methods assume that sequence sampling times are uninformative about population size dynamics. Recent work has started to question this assumption by exploring how sample time information can aid coalescent inference. Here we develop, investigate, and implement a new skyline method, termed the epoch sampling skyline plot (ESP), to jointly estimate the dynamics of population size and sampling rate through time. The ESP is inspired by real-world data collection practices. It comprises a flexible model in which the sequence sampling rate is proportional to the population size within an epoch, but can change discontinuously between epochs. We show that the ESP achieves accurate estimates under several realistic sampling protocols and we prove analytically that it can at least double the best precision attainable by standard approaches. We generalise the ESP to incorporate phylogenetic uncertainty in a new Bayesian package (BESP) in BEAST2. We further examine two well-studied empirical datasets from virus epidemiology and macroevolution. We find that the BESP improves upon previous coalescent estimators and generates new and biologically-useful insights into the sampling protocols underpinning these empirical datasets. Sequence sampling times  comprise a rich source of information for coalescent inference  that will become increasingly important as sequence collection intensifies and becomes more formalised._


## Prerequisites

### R


### Python
- numpy
- yaml
- BioPython
- TreeTime

### BEAST2
BEAST v2.5.0 with the BESP package installed. For now we just provide a `.jar` file on this repository until the package is ready to be released. 



## Bayesian simulation study

Simulate 100 replicate trees under a density-defined (linear preferential sampling) protocol and a frequency-defined (constant number of samples per epoch) protocol with 8 and 24 sampling epochs. Trees are simulated under five demographic scenarios (1) constant-size, (2) bottleneck, (3) boom-bust, (4) cyclical boom-bust and (5) logistic growth and decline. 

XML files are then created for each simulated tree and run in BEAST2 (only for frequency-defined sampling with 24 epochs). Finally, analyses are checked for convergence and summary statistics are computed in R.

To make the simulations run faster reduce the number of replicates, or increase the lower bound on the population size used to simulate coalescent trees (`SimUtils.R`, line 78). Keep in mind that if the lower bound is too high the demographic function for each simulated tree will turn into a constant-size coalescent once this bound is reached.

1. **Simulate trees:** 
	- Input: 
		- _N/A_
	- Scripts: 
		- `scripts/palettes.R`
		- `scripts/SimUtils.R`
	- Workflow: 
		- `workflows:simulate_trees.Rmd`	
	- Output: 
		- `results/simulations/`: Stored in `.RData`, `.csv`, `.trees` and `.json` files. _(Only the `.trees` and `.json` files are uploaded on the repository)._
2. **Create XML files:** 
	- Input: 
		- `templates/`: Template BEAST2 XML files for fixed tree analyses with the BESP with placeholders for data and parameters. 
		- `results/simulation_results/linear_epoch_24/config/:` Configuration files for the five different demographic scenarios.
	- Scripts: 
 		- `scripts/beastutils.py`
		- `scripts/MakeBEASTXML.py`
	- Workflow:
		- Run: `python MakeBEASTXML.py -i ../results/simulation_results/linear_epoch_24/config/` from `scripts/`
	- Output:
		- `results/simulation_results/linear_epoch_24/`: XML files and bash script for running XML files.
3. **Run simulation inferences:**
	- Input: 
		- `results/simulation_results/linear_epoch_24/`
	- Workflow:
		- Run the bash script in the directory of each of the demographic scenarios and provide the BEAST2 `.jar` file as the first argument.
		- Or use GNU parallel to speed up the process if you have many cores, e.g. run: `ls *.xml | parallel --delay 1 --jobs 75% --results outdir -I% --max-args 1 java -jar ../../../../BEAST2_BESP.jar -overwrite -seed 127 % &` from `results/simulation_results/linear_epoch_24/bottleneck/` 
	- Output:  
		- `results/simulation_results/linear_epoch_24/`: 
4. **Check convergence and compute summary statistics:**



## Case study 1: Seasonal Human Influenza

1. **Raw data:**
	- `data/H3N2/HA.fas` and `data/H3N2/HA.csv`: Full HA dataset from Rambaut _et al._ (2008) (687 sequences). 
	- `data/H3N2/HA_trim.fas` and `data/H3N2/HA_trim.csv`: HA dataset from Rambaut _et al._ (2008) with the first incomplete season removed (637 sequences).
2. **Initial trees:**
	- Workflow: (Run from `results/initial_trees/H3N2_trim/`)
		1. RAxML tree with 10 random starts

			```
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -T 2 -# 10 -p 12345 -s ../../H3N2/HA_trim.fas -n HA_trim.GTR+G
			
			```

		2. Rapid bootstrapping with 100 replicates (unnecessary if you only want a dated tree)

			```
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -# 100 -p 12345 -x 12345 -s ../../H3N2/HA_trim.fas -n HA_trim.GTR+G.RAPIDBS
			
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f b -z RAxML_bootstrap.HA_trim.GTR+G.RAPIDBS -t RAxML_bestTree.HA_trim.GTR+G -n HA_trim.GTR+G.tree
				
			```

		3. TreeTime dated tree
	
			```
			treetime --aln ../../H3N2/HA_trim.fas --tree RAxML_bipartitionsBranchLabels.HA_trim.GTR+G.tree --dates ../../H3N2/HA_trim.csv --coalescent skyline --confidence --outdir treetime_HA
			
			``` 
3. **Run analyses:**
4. **Create MCC tree:**
5. **Check convergence and analyse results:**


## Case study 2: Steppe Bison

1. **Raw data:**
2. **Run analyses:**
3. **Check convergence and analyse results:**
