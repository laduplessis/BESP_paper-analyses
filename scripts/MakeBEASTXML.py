import os, sys, io, yaml
import numpy as np
from fnmatch import fnmatch
from subprocess import call
from optparse import OptionParser
from Bio import Phylo, Nexus

from beastutils import *

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

parser.add_option("-b","--beast",
                  dest = "beast",
                  default = "",
                  metavar = "path",
                  help = "Path to BEAST2 jar file [required]")

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

parser.add_option("-s","--seed",
                  dest = "seed",
                  default = "127",
                  metavar = "integer",
                  help = "Seed or comma separated list of seeds to use [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"


################################################################################################################################  


	

################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
	if (fnmatch(filename,config)):

		sys.stdout.write(filename+"\t"+config+"\n")

		# Load config file
		configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
		pars 	   = yaml.load(configfile)
	
		# Set BEAST specific parameters	
		seeds      = list(map(int, options.seed.split(',')))				
		basename   = pars["name"] if options.name == '' else pars["name"]+"_"+options.name
		outputpath = os.path.abspath(pars["outputpath"] if options.outputpath == '' else options.outputpath)
		template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()
		treefile   = open(pars["trees"], 'r')

		# Output scripts
		if (not os.path.exists(outputpath)):
			os.makedirs(outputpath)
		scriptfile = open(outputpath+"/"+basename+".sh",'w')

		i = 0		
		for tree in treefile:
				
			# Set parameters not in the config file
			formatPars(pars, tree)

			# Create XML file			
			pars["name"] = basename+".T"+str(i)			
			makeXMLFile(pars, template, outputfile=pars["name"], outputpath=outputpath)

			# Write command to scripts
			for seed in seeds:
				cmd ="java -jar -Xms2G -Xmx4G $1 -overwrite -seed %d %s" % (seed, pars["name"]+".xml") 
				scriptfile.write("%s\n" % (cmd))

			i += 1
		#
		treefile.close()		
		scriptfile.close()
	#
#
