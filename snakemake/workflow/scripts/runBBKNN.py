#!/usr/bin/env python
#evaluate the DL model on simulated data
import sys
import os
import time
import scanpy as sc
import pandas as pd
import numpy as np
import pdb
import bbknn

import FuncVizBench as vb

inAdataFile = snakemake.input["inAdata"]

bbknnLatentOutFile = snakemake.output["bbknnLatent"]
bbknnClustersOutFile = snakemake.output["bbknnCluster"]

paramDict = snakemake.params["paramVals"]

batchLabel = paramDict["batchName"]

try:
	cellLabel = paramDict["cellName"]
except:
	cellLabel = ""

try:
	res = float(paramDict["res"])
except:
	res= 0.3
#nNeighborsUMAP = paramDict["nNeighborsUMAP"]
#inters = paramDict["inters"]
#cellsPerIter = paramDict["cellsPerIter"]
#nNeighborsKbet = paramDict["nNeighborsKbet"]
#	outAdataFile: {outAdataFile} \n\

allVars = f"inAdataFile: {inAdataFile} \n\
	batchLabel: {batchLabel} \n\
	cellLabel: {cellLabel}\n\
	is NA: {pd.isna(cellLabel)}\n\
	res: {res}"
	#numOutLayer: {numOutLayer} nNeighborsUMAP: {nNeighborsUMAP} inters: {inters} cellsPerIter: {cellsPerIter}, nNeighborsKbet {nNeighborsKbet} \n\
#print(allVars)

outClustStatDir = os.sep.join(bbknnClustersOutFile.split("/")[:-2]+["figures"])+os.sep
sc.settings.figdir = outClustStatDir

if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

adata = sc.read(inAdataFile)

clustKey = "BBKNN"

sc.tl.pca(adata, svd_solver='arpack', n_comps=100)

vb.testClustAndStats(adata, 
					umapKey = None, 
					neighborsKey=clustKey, 
					pcaRep="X_pca", 
					cellLabel=cellLabel,
					batchLabel=batchLabel, 
					res=res*0.25,
					numOutLayer=40, 
					outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 
bbknnLatent= adata.obsm["X_umap"].copy()

pd.DataFrame(bbknnLatent).to_csv(bbknnLatentOutFile,header=False, index=False)
pd.DataFrame(adata.obs[clustKey]).to_csv(bbknnClustersOutFile)




























