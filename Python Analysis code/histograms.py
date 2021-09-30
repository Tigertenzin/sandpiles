import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
# import scipy as sp
# from scipy.sparse import linalg as ln
from IPython.display import HTML
from random import choices
import math
import itertools
from os import path
from decimal import Decimal 
from scipy.optimize import curve_fit
from scipy import stats
import os.path
from os import path

try:
    plt.style.use('classic')
except:
    pass


def readSizes(seed, filetype):
    """
    Read the file of avalanche sizes for given seed. 
    
    Input: 
        - seed for run (contained in filename)

    Return: 
        - sizeS_Raw: (type:List) of each avalanche's size, S. 
        - params: List of parameters' values
    """
    
    with open("sp" + str(seed) + '_' + filetype + '.txt', "r") as f:
        header=[x for x in next(f)]
        params = [float(x) for x in next(f).split()]
        array = [[int(x) for x in line.split() if '\x00' not in x ] for line in f]

    sizeS_raw = [x[0] for x in array if len(x) >= 1]
    params[1] = int(params[1])
    
    return sizeS_raw, params
    
    
def plotSizes_raw(seed, filetype):
    """
    Plot log-log size histogram (log-binned)""
    """
    sizeS, params = readSizes(seed,filetype)
    sizeS_Log = np.asarray(sizeS)
    hist, bin_edges = np.histogram(sizeS_Log,bins=np.logspace(0,25,50,base=2),density=True)
    for i in range(len(hist)-1,-1,-1):      # Delete any 0 counts. 
        if hist[i]==0:
            hist = np.delete(hist,i)
            bin_edges = np.delete(bin_edges,i)   
    return bin_edges[:-1], hist, params


def writeSizes(seed, filetype): 
    """
    Write the histogram of avalanche sizes to a file (for quicker plotting later on...)
    """
    ## Create histogram from raw data
    x,y,params = plotSizes_raw(seed, filetype)
    z = np.column_stack((x, y))
    
    ## Write params to the file, commented with '#'
    filename = "Histo_fortran/sp_hist" + str(seed) + "_" + filetype + ".txt"
    listToStr = ','.join([str(elem) for elem in params])
    
    ## write histogram array to file
    np.savetxt(filename, z, delimiter=',', header=listToStr, comments='#')

# ===============================================================
filetypes = ['s', 'a', 'd']
heightRange = np.round(2**(np.unique(np.arange(1,11.5, 0.5))))
lambdaRange = np.logspace(-6, -1, 25).astype(float)

setRange = np.zeros((len(heightRange), len(lambdaRange),2))
for j in range(len(heightRange)): 
    for k in range(len(lambdaRange)): 
        setRange[j,k,0] = heightRange[j]
        setRange[j,k,1] = lambdaRange[k]
        
fixed_SetRange = setRange.transpose(2,0,1).reshape(2,-1)
print(np.shape(fixed_SetRange))

# ===============================================================
for filetype in filetypes: 
    for i in range(len(fixed_SetRange[0,:])): 
        seed = i + 100
        
        filename = "sp" + str(seed) + "_" + filetype + ".txt"
        if path.exists(filename): 
            print(str(seed)+filetype, end=' ')
            writeSizes(seed, filetype)
        else: 
            print(str(seed)+filetype+"DNE", end=' ')