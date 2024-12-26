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

import matplotlib.pyplot as plt
import matplotlib.image as img

import pdb
import torch

#sys.path.insert(0, 'scripts/scVitalPackage/src')
import scVital as sv

import FuncVizBench as vb

inAdataFile = snakemake.input["inAdata"]

scVitalModelOutFile= snakemake.output["outModel"]
scVitalLatentOutFile= snakemake.output["scVitalLatent"]
scVitalClusterOutFile=snakemake.output["scVitalCluster"]
outAdataFile = snakemake.output["outAdata"]
outCombFigFile = snakemake.output["outCombFig"]

paramDict = snakemake.params["paramVals"]
batchLabel = paramDict["batchName"]
batchSize= int(paramDict['batchSize'])
numEpoch = int(paramDict["numEpoch"])#int64
learningRate = float(paramDict["learningRate"])#np.float64(
lastLayer = int(paramDict["lastLayer"])
inLayerDims = [int(l) for l in paramDict["inLayerDims"].split("-")]
inDiscLayer = int(paramDict["inDiscLayer"])
reconCoef = float(paramDict["reconCoef"]) 
klCoef = float(paramDict["klCoef"])
discCoef = float(paramDict["discCoef"])

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

allVars = f"inAdataFile: {inAdataFile} \n\
	outAdataFile: {outAdataFile} \n\
	batchLabel: {batchLabel} \n\
	cellLabel: {cellLabel}\n\
	is NA: {pd.isna(cellLabel)}\n\
	res: {res}"
	# nNeighborsUMAP: {nNeighborsUMAP} inters: {inters} cellsPerIter: {cellsPerIter}, nNeighborsKbet {nNeighborsKbet} \n\
#print(allVars)

outClustStatDir = os.sep.join(scVitalClusterOutFile.split("/")[:-2]+["figures"])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir
if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

clustKey = "scVital"
latentRep = "X_scVital"

adata = sc.read(inAdataFile)


scVitalModel = sv.makeScVital(adata, 
							  batchLabel = batchLabel, 
							  miniBatchSize = batchSize, 
							  numEpoch = numEpoch, 
							  learningRate = learningRate,
							  hid1 = inLayerDims[0], 
							  hid2 = inLayerDims[1], 
							  latentSize = lastLayer, 
							  discHid = inDiscLayer, 
							  reconCoef = reconCoef, 
							  klCoef = klCoef, 
							  discCoef = discCoef, 
							  train = False, 
							  seed = 18, 
							  verbose = True)
print(scVitalModel)

startTrain = time.process_time() 

scVitalModel.runTrainScVital()

endTrain = time.process_time()
print(f"{(endTrain-startTrain)/60} mins")
scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", "scVital"])
scale.to_csv(f"{outClustStatDir}scale_scVital.csv")

adata = scVitalModel.getAdata()


vb.testClustAndStats(adata, umapKey = "scVital", neighborsKey=clustKey, pcaRep=latentRep, 
					cellLabel=cellLabel, batchLabel=batchLabel, res=res,
					numOutLayer=lastLayer, outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 


adata.write(outAdataFile)
scVitalLatent= adata.obsm[latentRep].copy()

pd.DataFrame(scVitalLatent).to_csv(scVitalLatentOutFile,header=False, index=False)
pd.DataFrame(adata.obs[clustKey]).to_csv(scVitalClusterOutFile)

outClustStatDir = os.sep.join(outCombFigFile.split(os.sep)[:-1])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir
if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)
	
scVitalModel.plotLoss(outCombFigFile)

scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", "scVital"])
scale.to_csv(f"{outClustStatDir}scale_scVital.csv")

scVitalModel.saveModel(scVitalModelOutFile)

