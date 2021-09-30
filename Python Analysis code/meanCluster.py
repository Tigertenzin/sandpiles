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
    
def calc_chi_disp(seed_list, filetype):
    chi_s = np.zeros_like(seed_list);  lambda_s = np.zeros_like(seed_list).astype('float')
    for s in range(len(seed_list)):
        print(seed_list[s], end=filetype+' ')
        sizeS, params = readSizes(seed_list[s],filetype)
        lambda_s[s] = params[2]
        chi_s[s] = np.sum(np.power(sizeS, 2)) / np.sum(np.power(sizeS, 1))
    return lambda_s, chi_s

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
for hc in heightRange: 
    for filetype in filetypes: 
        indexRange = []
        for j in range(len(fixed_SetRange[0,:])): 
            if fixed_SetRange[0,j] == hc:
                filename = "sp" + str(j+100) + "_" + filetype + ".txt"
                if path.exists(filename):
                    indexRange += [j+100] 
                    
        print(filetype, hc, indexRange[0], indexRange[-1], end=' ')
        
        lambda_x, chi_x = calc_chi_disp(indexRange, filetype)
    
        results = np.stack((lambda_x, chi_x))
    
        filename = "mean_cluster/chi_"+filetype+str(hc)+".txt"
        np.savetxt(filename, results, delimiter=',')
    
    
    