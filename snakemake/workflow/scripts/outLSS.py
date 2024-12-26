import sys
import os
import random
import time
import scanpy as sc
import pandas as pd
import numpy as np
import pdb
import pickle

from matplotlib import pyplot as plt

import FuncVizBench as vb

#sys.path.insert(0, 'scripts/scVitalPackage/src')
import scVital as sv

#print(snakemake.input)
#print(snakemake.params)

inAdataFile = snakemake.input["inAdata"]
name = inAdataFile.split("~")[1].split("/")[0]

dataName = snakemake.params[0]["filename"] #inAdataFile.split("~")[1].split("/")[0]
batchName = snakemake.params[0]["batchName"] #inAdataFile.split("/")[4].split("_")[9].split("~")[1]
labelName = snakemake.params[0]["cellName"] #inAdataFile.split("/")[4].split("_")[10].split("~")[1]
dataDir = ("/").join(inAdataFile.split("/")[:-1])

#outBenchFile = snakemake.output["outBench"]
outLssFile = snakemake.output["outLSS"]
outAdataFile = snakemake.output["outAdata"]


latNameSet = set()
clutNameSet = set()

adata = sc.read_h5ad(inAdataFile)    #set adata  

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
aucDf = pd.DataFrame(np.zeros(len(latents)),columns=[dataName],index=latents).T	
batches = adata.obs[batchName].cat.categories.values

if(labelName != ""):

	for latent in latents:
		clustDist, lssAUC, totalDist, allCellTypes, ctPairs = sv.lss.calcPairsLSS(adata, latent=latent, batchName=batchName, cellTypeLabel=labelName)
		
		batchDict, annoToColorDict = sv.lss.plotGraphLSS(adata, labelName, batchName, clustDist, name=dataName, 
             ctColors=plt.get_cmap('tab20').colors, btColors=None, shapes="ospx^><.....", 
             prog="neato", wLab=False, qCut = 0.28, plot=False, save=os.path.join(dataDir,"figures"))
		
		sv.lss.plotHeatLSS(adata, clustDist, latent, allCellTypes, ctPairs, plot=False, save=os.path.join(dataDir,"figures"))
		aucDf[latent] = lssAUC

aucDf.to_csv(outLssFile)
adata.write(outAdataFile)  
