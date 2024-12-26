#!/usr/bin/env python
import sys
import os
import random
import time
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
from scipy.sparse import isspmatrix

import tensorflow as tf2
import tensorflow.compat.v1 as tf

import scipy.io
from sklearn.decomposition import PCA
import pdb

import pandas as pd
import scanpy as sc


import scipy.sparse
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from scipy import stats 
from scipy import * 
import datetime 
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

from scDREAMER.scDREAMER.model import scDREAMER as model

import FuncVizBench as vb


inAdataFile = snakemake.input["inAdata"]

#outAdataFile = snakemake.output["outAdata"]
scDreamLatentOutFile = snakemake.output["scDreamLatent"]
clustersOutFile = snakemake.output["scDreamCluster"]

paramDict = snakemake.params["paramVals"]
name = paramDict["filename"]
epoch = paramDict["numEpoch"]
learningRate = paramDict["learningRate"]
batchLabel = paramDict["batchName"]
res = 0.3
try:
	cellLabel = paramDict["cellName"]
except:
	cellLabel = ""

lastLayer = int(paramDict["lastLayer"])

#	outAdataFile: {outAdataFile} \n\

allVars = f"inAdataFile: {inAdataFile} \n\
	batchLabel: {batchLabel} \n\
	epoch: {epoch} \n\
	cellLabel: {cellLabel}\n\
	is NA: {pd.isna(cellLabel)}\n\
	"
#print(allVars)

outClustStatDir = os.sep.join(clustersOutFile.split("/")[:-2]+["figures"])+os.sep
#print(outClustStatDir)
sc.settings.figdir = outClustStatDir
if not os.path.exists(outClustStatDir):
	os.mkdir(outClustStatDir)

adata = sc.read(inAdataFile)
if (adata.raw != None):
	adata = adata.raw.to_adata()
#print(adata)
adata.write(f"{name}scDreamRaw.h5ad")
#print(name)
#pdb.set_trace()
dreamerOutFile = f"{name}latent_matrix_300.csv"

print(dreamerOutFile , (not os.path.exists(dreamerOutFile)))

clustKey = "scDREAMER"
latentRep = "X_scDREAM"

if (not os.path.exists(dreamerOutFile)):
	startTrain = time.process_time() 
	
	np.random.seed(0)
	tf.set_random_seed(0)
	random.seed(0)
	tf2.random.set_seed(0)
	tf2.keras.utils.set_random_seed(0)
	
	run_config = tf.ConfigProto()
	
	run_config.gpu_options.per_process_gpu_memory_fraction = 0.333
	run_config.gpu_options.allow_growth = True
	
	with tf.Session(config = run_config) as sess:
		dreamer = model(
			sess,
			#epoch = 10,
			dataset_name = f"{name}scDreamRaw.h5ad",
			batch = batchLabel,
			cell_type = cellLabel,
			name = name,
			lr_ae = 0.0001,
			lr_dis = 0.00001,
			num_layers = 2,
			#X_dim = adata.shape[1],
			z_dim = lastLayer
			)
	
		dreamer.train_cluster()
	
	endTrain = time.process_time()
	scale = pd.DataFrame(np.array([[adata.X.size, endTrain-startTrain]]), columns=["Size", "scDREAMER"])
	scale.to_csv(f"{outClustStatDir}scale_scDREAMER.csv")
	
os.rename(dreamerOutFile, scDreamLatentOutFile)
#print(scDreamLatentOutFile)

scDREAMLatent = pd.read_table(scDreamLatentOutFile,sep=",",header=None).to_numpy()
#print(f"\t {scDREAMLatent.shape}\n")
#print(adata)
adata.obsm[latentRep] = scDREAMLatent

numOutLayer = adata.obsm[latentRep].shape[1]

vb.testClustAndStats(adata, umapKey = clustKey, neighborsKey=clustKey, pcaRep=latentRep, 
					cellLabel=cellLabel, batchLabel=batchLabel, res=res,
					numOutLayer=numOutLayer, outClustStatDir=outClustStatDir)
#					nNeighborsKbet=nNeighborsKbet, inters=inters, cellsPerIter=cellsPerIter, outUMAPFile=outUMAPFile, 

adata.obs[clustKey].to_csv(clustersOutFile)

os.remove(f"{name}scDreamRaw.h5ad")

