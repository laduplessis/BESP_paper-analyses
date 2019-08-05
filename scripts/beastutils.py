import os, sys, io
import numpy as np
from subprocess import call
from Bio import Phylo, Nexus

def writeDateTrait(dates,  output):

	traits = []	
	for i in dates.keys():
		traits.append('\n\t\t\t\t\t%s=%.10f' %(i, dates[i]))
	output.write(','.join(traits)+'\n')
#

def writeAlignment(newicktree, output):

	tree   = Nexus.Trees.Tree(newicktree)	
	taxa   = tree.get_taxa()

	for i in range(0,len(taxa)):
		output.write('\n\t\t\t\t\t<sequence spec="Sequence" taxon="%d" value="?" />' % (i+1))
#


# Read sampling times from tree
# 
# If forwards = True  gives times as time from tMRCA
# If forwards = False gives times as time from most recent sample
#
def getSamplingTimes(newicktree, forwards=False):

	# Sampling times from tree
	tree    = Nexus.Trees.Tree(newicktree)
	leaves  = tree.get_terminals()
	heights = np.zeros(len(leaves))	
	for i in range(0,len(leaves)):
		heights[i] = tree.sum_branchlength(node=leaves[i])

	# Forwards or backwards
	if (forwards == True):
		treetimes = heights
	else:
		treetimes = max(heights) - heights

	# Associate with label
	times = dict()
	for i in range(0,len(treetimes)):
		label = tree.get_taxa(node_id=leaves[i])[0]
		times[label] = treetimes[i]		
	#	

	return times	
#


def makeXMLFile(pars, template, outputfile="", outputpath=""):

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.makedirs(outputpath)

	if (outputfile == ""):
		outputfile = pars["name"]

	sys.stdout.write(outputfile+"...\n")

	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	outfile = open(outputpath+"/"+outputfile+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars, tree):

	# Set inputtree, sampling times and dummy alignment
	output_align = io.StringIO()
	output_dates = io.StringIO()
	output_model = io.StringIO()


	# Add default date direction if not defined
	if ("dateTrait" not in pars.keys()):
		pars["dateTrait"] = "date"
	samplingTimes = getSamplingTimes(tree, forwards=(pars["dateTrait"] != "date-backward"))
	
	writeDateTrait(samplingTimes, output_dates)
	writeAlignment(tree, output_align)

	pars["tree"]  = tree
	pars["dates"] = output_dates.getvalue()
	pars["taxa"]  = output_align.getvalue()
#


def get_beast_call(beast, filename, seed):	
	return ["java", "-jar", beast, "-seed", str(seed), "-overwrite", filename]