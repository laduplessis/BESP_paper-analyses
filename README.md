# Jointly inferring the dynamics of population size and sampling intensity from molecular sequences: BESP analyses

### KV Parag, L du Plessis and OG Pybus

## Summary

This repository contains the raw data, configuration files and scripts necessary to reproduce the BESP analyses presented in [https://www.biorxiv.org/content/10.1101/686378v1.abstract](https://www.biorxiv.org/content/10.1101/686378v1.abstract). Some of the scripts may need some adjustment depending on the local setup. 

### Abstract

_Estimating past population dynamics from molecular sequences that have been sampled longitudinally through time is an important problem in infectious disease epidemiology, molecular ecology and macroevolution. Popular solutions, such as the skyline and skygrid methods, infer past effective population sizes from the coalescent event times of phylogenies reconstructed from sampled sequences, but assume that sequence sampling times are uninformative about population size changes. Recent work has started to question this assumption by exploring how sampling time information can aid coalescent inference. Here we develop, investigate, and implement a new skyline method, termed the epoch sampling skyline plot (ESP), to jointly estimate the dynamics of population size and sampling rate through time. The ESP is inspired by real-world data collection practices and comprises a flexible model in which the sequence sampling rate is proportional to the population size within an epoch, but can change discontinuously between epochs. We show that the ESP is accurate under several realistic sampling protocols and we prove analytically that it can at least double the best precision achievable by standard approaches. We generalise the ESP to incorporate phylogenetic uncertainty in a new Bayesian package (BESP) in BEAST2. We re-examine two well-studied empirical datasets from virus epidemiology and molecular evolution and find that the BESP improves upon previous coalescent estimators and generates new, biologically-useful insights into the sampling protocols underpinning these datasets. Sequence sampling times provide a rich source of information for coalescent inference that will become increasingly important as sequence collection intensifies and becomes more formalised._


## Dependencies

### BEAST2
BEAST v2.6.0 with the besp package installed. Alternatively, use the `.jar` file included with this repository.

### Python
- numpy
- yaml
- BioPython
- TreeTime

### CRAN R-packages
- ape
- coda
- phylodyn
- lubridate
- TreeSim

### Other R-packages
These packages are developed by Louis du Plessis and not available on CRAN yet. These packages can be installed from Github using `devtools::install_github()`.

- [beastio](https://github.com/laduplessis/beastio): 
	- _R-package with functions for pre- and post-processing of BEAST and BEAST2 files (good for automating the functionality in Tracer or Logcombiner e.g. checking convergence in hundreds of replicates from a simulation study or combining many chains)._
- [rskylinetools](https://github.com/laduplessis/rskylinetools): 
	- _R-package implementing function for gridding skylines (so this work doesn't need to be painstakingly done in Tracer)._



## Bayesian Simulation Study

Simulate 100 replicate trees under a density-defined (linear preferential sampling) protocol and a frequency-defined (constant number of samples per epoch) protocol with 8 and 24 sampling epochs. Trees are simulated under five demographic scenarios (1) constant-size, (2) bottleneck, (3) boom-bust, (4) cyclical boom-bust and (5) logistic growth and decline. 

XML files are then created for each simulated tree and run in BEAST2. Finally, analyses are checked for convergence and summary statistics are computed in R. The results report contains some plots from example analyses, but plots for every replicate are stored in `results/`. Runs where the automatic check for convergence failed can be checked in Tracer (the ESS calculation is slightly different in coda to Tracer).

To make the simulations run faster reduce the number of replicates, or increase the lower bound on the population size used to simulate coalescent trees (`SimUtils.R`, line 78). Keep in mind that if the lower bound is too high the demographic function for each simulated tree will turn into a constant-size coalescent once this bound is reached.

To save space we only provide XML files and output figures for the 24-epoch BESP results and only figures on a log-scale. The output workflow for the single-epoch BESP (linear preferential sampling) are provided for comparison (but are not presented in the paper).

1. **Simulate trees:** 
	- _Input:_ N/A
	- _Scripts:_ 
		- `scripts/palettes.R`
		- `scripts/SimUtils.R`
	- _Workflow:_  `workflows:simulate_trees.Rmd`	
	- _Output:_ 
		- `results/simulations/`: Stored in `.RData`, `.csv`, `.trees` and `.json` files. _(Only the `.trees` and `.json` files are uploaded on the repository)._
2. **Create XML files:** 
	- _Input:_
		- `templates/`: Template BEAST2 XML files for fixed tree analyses with the BESP with placeholders for data and parameters. 
		- `results/simulation_results/linear_epoch_24/config/:` Configuration files for the five different demographic scenarios.
	- _Scripts:_ 
 		- `scripts/beastutils.py`
		- `scripts/MakeBEASTXML.py`
	- _Workflow:_ Run from `scripts/`
	
		```bash
		python MakeBEASTXML.py -i ../results/simulation_results/linear_epoch_24/config/
		```
	- _Output:_
		- `results/simulation_results/linear_epoch_24/`: XML files and bash script for running XML files.
3. **Run simulation inferences:**
	- _Input:_ `results/simulation_results/linear_epoch_24/`
	- _Workflow:_ Run the bash script in the directory of each of the demographic scenarios and provide the BEAST2 `.jar` file as the first argument or use GNU parallel to speed up the process if you have many cores, e.g. run from the directory with the XML files: 

		```bash
		ls *.xml | parallel --delay 1 --jobs 75% --results outdir -I% --max-args 1 java -jar ../../../../BEAST2_BESP.jar -overwrite -seed 127 % &` from `results/simulation_results/linear_epoch_24/bottleneck/
		```
	- _Output:_  
		- `results/simulation_results/linear_epoch_24/`: 
4. **Check convergence and compute summary statistics:**
	- _Input:_ `results/simulation_results/linear_epoch_24/`
	- _Scripts:_
		- `scripts/palettes.R`
		- `scripts/SimUtils.R`
		- `scripts/SimResultUtils.R`
	- _Workflow:_ `workflows/simulation_results_linear_epoch_24/`
	- _Output:_ 
		- `results/simulation_results/linear_epoch_24/`: Summary statistic are stored in `.csv` files (population size statistics are calculated from the present to the most ancient sample and from the present to the TMRCA). For each replicate 4 pdf figures are produced showing the results in real/log space and the results in real/log space cut off at the most ancient sample.


## Case study 1: Seasonal Human Influenza

1. **Raw data:**
	- `data/H3N2/HA.fas` and `data/H3N2/HA.csv`: Full HA dataset from Rambaut _et al._ (2008) (687 sequences). 
	- `data/H3N2/HA_trim.fas` and `data/H3N2/HA_trim.csv`: HA dataset from Rambaut _et al._ (2008) with the first incomplete season removed (637 sequences).
2. **Initial trees:**
	- _Input:_ `data/H3N2/`
	- _Workflow:_ Run from `results/initial_trees/H3N2_trim/`:
		1. RAxML tree with 10 random starts

			```bash
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -T 2 -# 10 -p 12345 -s ../../../data/H3N2/HA_trim.fas -n HA_trim.GTR+
			```

		2. Rapid bootstrapping with 100 replicates (unnecessary if you only want a dated tree)

			```bash
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -# 100 -p 12345 -x 12345 -s ../../../data/H3N2/HA_trim.fas -n HA_trim.GTR+G.RAPIDBS
			
			raxmlHPC-PTHREADS-AVX -m GTRGAMMA -f b -z RAxML_bootstrap.HA_trim.GTR+G.RAPIDBS -t RAxML_bestTree.HA_trim.GTR+G -n HA_trim.GTR+G.tree
			```

		3. TreeTime dated tree
	
			```bash
			treetime --aln ../../../data/H3N2/HA_trim.fas --tree RAxML_bipartitionsBranchLabels.HA_trim.GTR+G.tree --dates ../../../data/H3N2/HA_trim.csv --coalescent skyline --confidence --outdir treetime_HA
			``` 
3. **Run analyses:**
	- _Input:_
		- `results/H3N2/HA_trim/HA.HKY+G+F.BESP40_12.xml`: BESP with 40 segments, 12 epochs and initial tree calculated above.
		- `results/H3N2/HA_trim/HA.HKY+G+F.BESP40_1.xml`: BESP with 40 segments, 1 epoch and initial tree calculated above.
		- `results/H3N2/HA_trim/HA.HKY+G+F.BSP40.xml`: BSP with 40 segments and initial tree calculated above.
	- _Workflow:_ Run 7 chains for each XML file, starting from different seeds (this is best done on a remote server since it will take a long time to run):
 
		```bash
		mkdir output
		for SEED in 125 126 127 128 129 130 131 
		do 
			nohup java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED -threads 2 HA.HKY+G+F.BSP40.xml > HA.HKY+G+F.BSP40.${SEED}.out &
			nohup java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED -threads 2 HA.HKY+G+F.BESP40_1.xml > HA.HKY+G+F.BESP40_1.${SEED}.out &
			nohup java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED -threads 2 HA.HKY+G+F.BESP40_12.xml > HA.HKY+G+F.BESP40_12.${SEED}.out &
		done
		
		```		
	- _Output:_
		- `results/H3N2/HA_trim/output/`
5. **Check convergence and analyse results:**
	- _Input:_ `results/H3N2/HA_trim/output/`
	- _Workflow:_ `workflows/HA_trim_results.Rmd`
	- _Output:_ `results/H3N2/HA_trim/figures/`

## Case study 2: Steppe Bison

1. **Raw data:**
	- `data/Bison/bison_2013.fasta` and `data/Bison/bison_2013_taxa.csv`: Bison dataset from Gill _et al._ (2013) (152 sequences).
2. **Run analyses:**
	- _Input:_
		- `results/Bison/Bison.HKY.BESP20_12.xml`: BESP with 20 segments and 12 epochs.
		- `results/Bison/Bison.HKY.BESP20_1.xml`: BESP with 20 segments and 1 epoch.
		- `results/Bison/Bison.HKY.BSP20.xml`: BSP with 20 segments.
	- _Workflow:_ Run 3 chains for each XML file, starting from different seeds:

		```bash
		mkdir output
		for SEED in 125 126 127
		do 
			java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED Bison.HKY.BSP20.xml > Bison.HKY.BSP20.${SEED}.out&
			java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED Bison.HKY.BESP20_1.xml > Bison.HKY.BESP20_1.${SEED}.out
			java -jar ../../../BEAST2_BESP.jar -overwrite -seed $SEED Bison.HKY.BESP20_12.xml > Bison.HKY.BESP20_12.${SEED}.out
		done		
		```

	- _Output:_
		- `results/Bison/output/`
3. **Check convergence and analyse results:**
 	- _Input:_ `results/Bison/output/`
	- _Workflow:_ `workflows/bison_results.Rmd`
	- _Output:_ `results/Bison/figures/`
