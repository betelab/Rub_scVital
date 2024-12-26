import sys
import os
import random
import time
import scanpy as sc
import pandas as pd
import numpy as np
import pdb
import pickle
from scib_metrics.benchmark import Benchmarker

from matplotlib import pyplot as plt

import FuncVizBench as vb

#print(snakemake.input)
#print(snakemake.params)

inAdataFile = snakemake.input["inAdata"]
name = inAdataFile.split("~")[1].split("/")[0]

dataName = snakemake.params[0]["filename"] #inAdataFile.split("~")[1].split("/")[0]
batchName = snakemake.params[0]["batchName"] #inAdataFile.split("/")[4].split("_")[9].split("~")[1]
labelName = snakemake.params[0]["cellName"] #inAdataFile.split("/")[4].split("_")[10].split("~")[1]
dataDir = ("/").join(inAdataFile.split("/")[:-1])

#outBenchFile = snakemake.output["outBench"]
outAdataFile = snakemake.output["outAdata"]

latNameSet = set()
clutNameSet = set()

if(labelName != ""):
	adata = sc.read(inAdataFile)    #set adata  
	
	for infoDir in os.listdir(os.path.join(dataDir)):##in directory with latents and clusters
	    if(os.path.isdir(os.path.join(dataDir,infoDir))): 
	        for file in os.listdir(os.path.join(dataDir,"latents")):
	            if(".csv" in file):
	                name = file.split(".")[0]
	                adata.obsm[f"X_{name}"] = pd.read_table(os.path.join(dataDir,"latents",file), sep=",",header=None).to_numpy()
	                latNameSet.add(f"X_{name}")
	
	        for file in os.listdir(os.path.join(dataDir,"clusters")):
	            if(".csv" in file):
	                name = file.split(".")[0]
	                adata.obs[name] = pd.read_table(os.path.join(dataDir,"clusters",file), sep=",",index_col=0,dtype=str).iloc[:,0].to_numpy()
	                clutNameSet.add(f"{name}")
	
	latents = list(latNameSet)
	batches = adata.obs[batchName].cat.categories.values

	bm = Benchmarker(adata,
					 batch_key = batchName, 
					 label_key = "overlapLabel",
					 embedding_obsm_keys = latents,
					 n_jobs=len(latents))
	bm.benchmark()
	benchDF = bm.get_results(min_max_scale=False)
	benchDF.to_csv(os.path.join(dataDir,"figures","bench.csv"))
	bm.plot_results_table(min_max_scale=False, show=False, save_dir=os.path.join(dataDir,"figures"))
	
	adata.write(outAdataFile)





