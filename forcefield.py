import sys
import subprocess
import numpy as N
from sgdml.train import GDMLTrain
from sgdml.predict import GDMLPredict
from sgdml.utils import io
#This module contains all routines for training GDML and sGDML models.
#extendxyz=open("trajectories.xyz","r")
subprocess.call(["cat trajectories.*.ani > trajectories.xyz"],shell=True)
subprocess.call(["sgdml_dataset_from_xyz.py trajectories.xyz"],shell=True)
#Force field reconstruction
dataset = N.load('trajectories.npz')
n_train = 200

gdml_train = GDMLTrain()
task = gdml_train.create_task(dataset, n_train,\
        valid_dataset=dataset, n_valid=1000,\
        sig=10, lam=1e-15)

try:
        model = gdml_train.train(task)
except Exception, err:
        sys.exit(err)
else:
        N.savez_compressed('FFftraj.npz', **model)
        
#Force field query
#model = N.load('FFftraj.npz')
gdml = GDMLPredict(model)

r,_ = io.read_xyz('target.xyz') 
e,f = gdml.predict(r)

print r.shape 
print e.shape 
print f.shape 