#!/usr/bin/env python
#evaluate the DL model on simulated data
import sys
import os
import time
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
from scipy.sparse import isspmatrix

import pdb

import FuncVizBench as vb

inAdataFile = snakemake.input["inAdata"]

normLatentOutFile = snakemake.output["normalLatent"]
normClustersOutFile=snakemake.output["normalCluster"]

paramDict = snakemake.params["paramVals"]

#print(paramDict)

batchLabel = paramDict["batchName"]

try:
	cellLabel = paramDict["cellName"]
except:
	cellLabel = ""

#
numOutLayer = 40#paramDict["lastLayer"]

try:
	res = float(paramDict["res"])
except:
	res= 0.3
#nNeighborsUMAP = paramDict["nNeighborsUMAP"]
#inters = paramDict["inters"]
#cellsPerIter = paramDict["cellsPerIter"]
#nNeighborsKbet = paramDict["nNeighborsKbet"]
#	outAdataFile: {outAdataFile} \n

allVars = f"inAdataFile: {inAdataFile} \n\
	batchLabel: {batchLabel} \n\
	cellLabel: {cellLabel}\n\
	is NA: {pd.isna(cellLabel)}\n\
	res: {res}"
	#numOutLayer: {numOutLayer} nNeighborsUMAP: {nNeighborsUMAP} inters: {inters} cellsPerIter: {cellsPerIter}, nNeighborsKbet {nNeighborsKbet} \n\
#print(allVars)

outClustStatDir = os.sep.join(normClustersOutFile.split("/")[:-2]+["figures"])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir

if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

adata = sc.read(inAdataFile)

clustKey = "normal"
latentRep = "X_pca"

startTrain = time.process_time() 
sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
endTrain = time.process_time()
scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", clustKey])
scale.to_csv(f"{outClustStatDir}scale_normal.csv")

vb.testClustAndStats(adata, umapKey = "norm", neighborsKey=clustKey, pcaRep=latentRep, 
					cellLabel=cellLabel, batchLabel=batchLabel, res=res,
					numOutLayer=numOutLayer, outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 

normLatent = adata.obsm[latentRep].copy()

pd.DataFrame(normLatent).to_csv(normLatentOutFile,header=False, index=False)
pd.DataFrame(adata.obs[clustKey]).to_csv(normClustersOutFile)






