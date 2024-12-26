#!/usr/bin/env python
#evaluate the DL model on simulated data
import sys
import os
import time
import scanpy as sc
import pandas as pd
import numpy as np

import harmonypy

import FuncVizBench as vb


inAdataFile = snakemake.input["inAdata"]

harmonyLatentOutFile = snakemake.output["harmonyLatent"]
harmonyClusterOutFile = snakemake.output["harmonyCluster"]

paramDict = snakemake.params["paramVals"]

batchLabel = paramDict["batchName"]

try:
	cellLabel = paramDict["cellName"]
except:
	cellLabel = ""

#numOutLayer = paramDict["lastLayer"]

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

outClustStatDir = os.sep.join(harmonyClusterOutFile.split("/")[:-2]+["figures"])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir

if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

adata = sc.read(inAdataFile)

clustKey = "Harmony"
latentRep = "X_pca_harmony"

sc.tl.pca(adata, svd_solver='arpack', n_comps=100)

startTrain = time.process_time() 
sc.external.pp.harmony_integrate(adata, key=batchLabel)  #has to be str
endTrain = time.process_time()
scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", clustKey])
scale.to_csv(f"{outClustStatDir}scale_Harmony.csv")

vb.testClustAndStats(adata, umapKey="harm", neighborsKey=clustKey, pcaRep=latentRep, 
					cellLabel=cellLabel, batchLabel=batchLabel, res=res*0.25,
					numOutLayer=40, outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 

harmonyLatent= adata.obsm[latentRep].copy()

pd.DataFrame(harmonyLatent).to_csv(harmonyLatentOutFile,header=False, index=False)
pd.DataFrame(adata.obs[clustKey]).to_csv(harmonyClusterOutFile)







