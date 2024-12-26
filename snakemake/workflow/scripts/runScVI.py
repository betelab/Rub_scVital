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

import scvi
from scvi.model.utils import mde

import pdb
import torch

import FuncVizBench as vb

inAdataFile = snakemake.input["inAdata"]

scViLatentOutFile = snakemake.output["scViLatent"]
scViClustersOutFile=snakemake.output["scViClusters"]

paramDict = snakemake.params["paramVals"]

batchLabel = paramDict["batchName"]
lastLayer = int(paramDict["lastLayer"])
epoch = paramDict["numEpoch"]

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
#outAdataFile: {outAdataFile} \n\

allVars = f"inAdataFile: {inAdataFile} \n\
	batchLabel: {batchLabel} \n\
	cellLabel: {cellLabel}\n\
	is NA: {pd.isna(cellLabel)}\n\
	res: {res}"
	#numOutLayer: {numOutLayer} nNeighborsUMAP: {nNeighborsUMAP} inters: {inters} cellsPerIter: {cellsPerIter}, nNeighborsKbet {nNeighborsKbet} \n\
#print(allVars)

outClustStatDir = os.sep.join(scViClustersOutFile.split("/")[:-2]+["figures"])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir
if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

adata = sc.read(inAdataFile)

clustKey = "scVI"
latentRep = "X_scVI"

startTrain = time.process_time() 
numOutLayer = lastLayer

scvi.model.SCVI.setup_anndata(adata, batch_key=batchLabel, layer="counts")
scviVAE = scvi.model.SCVI(adata, n_layers=2, n_latent=numOutLayer, gene_likelihood="nb")
scviVAE.train()#max_epochs=10, early_stopping=True)#int(epoch)
scViLatent = scviVAE.get_latent_representation()
adata.obsm[latentRep] = scViLatent
endTrain = time.process_time()
scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", "scVI"])
scale.to_csv(f"{outClustStatDir}scale_scVI.csv")

vb.testClustAndStats(adata, umapKey = clustKey, neighborsKey=clustKey, pcaRep=latentRep, 
					cellLabel=cellLabel, batchLabel=batchLabel, res=res,
					numOutLayer=numOutLayer, outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 


pd.DataFrame(scViLatent).to_csv(scViLatentOutFile,header=False, index=False)
pd.DataFrame(adata.obs[clustKey]).to_csv(scViClustersOutFile)






