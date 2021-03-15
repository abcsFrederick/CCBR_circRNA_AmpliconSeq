import shutil
import sys
import os
import pandas as pd
import yaml
import glob

CONFIGFILE = str(workflow.overwrite_configfiles[0])

def check_existence(filename):
	"""Checks if file exists on filesystem
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	if not os.path.exists(filename):
		sys.exit("File: {} does not exists!".format(filename))
	return True


def check_readaccess(filename):
	"""Checks permissions to see if user can read a file
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	check_existence(filename)
	if not os.access(filename,os.R_OK):
		sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))
	return True


def check_writeaccess(filename):
	"""Checks permissions to see if user can write to a file
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	check_existence(filename)
	if not os.access(filename,os.W_OK):
		sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))
	return True


#
MEMORY="100"
# get working dir from config
WORKDIR = config["workdir"]

# get resources folder
try:
	RESOURCESDIR = config["resourcesdir"]
except KeyError:
	RESOURCESDIR = join(WORKDIR,"resources")

# get scripts folder
try:
	SCRIPTSDIR = config["scriptsdir"]
except KeyError:
	SCRIPTSDIR = join(WORKDIR,"scripts")

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)

if not os.path.exists(join(WORKDIR,"fastqs")):
	os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(join(WORKDIR,"results")):
	os.mkdir(join(WORKDIR,"results"))
for f in ["samples", "tools", "cluster"]:
	check_readaccess(config[f])

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["R2"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["PEorSE"]="PE"

for sample in SAMPLES:
	R1file=SAMPLESDF["path_to_R1_fastq"][sample]
	R2file=SAMPLESDF["path_to_R2_fastq"][sample]
	# print(sample,R1file,R2file)
	check_readaccess(R1file)
	R1filenewname=join(WORKDIR,"fastqs",sample+".R1.fastq.gz")
	if not os.path.exists(R1filenewname):
		os.symlink(R1file,R1filenewname)
	SAMPLESDF.loc[[sample],"R1"]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
		if not os.path.exists(R2filenewname):
			os.symlink(R2file,R2filenewname)
		SAMPLESDF.loc[[sample],"R2"]=R2filenewname
	else:
		SAMPLESDF.loc[[sample],"PEorSE"]="SE"
# print(SAMPLESDF)
# sys.exit()


