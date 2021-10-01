from os import path
import numpy as np
import matplotlib as mpl

def readRG(seed, filetype):
    """
    Read the file of avalanche sizes for given seed. 
    
    Input: 
        - seed for run (contained in filename)
    
    Return: 
        - sizeS: (type:List) of each avalanche's size, S. 
        - rg: (type: List) of each avalanche's radius of gyration, r_G
        - params: List of parameters' values
    """
    # Read in from size file
    with open("sp" + str(seed) + '_' + filetype + '.txt', "r") as f:
        header=[x for x in next(f)]
        params = [float(x) for x in next(f).split()]
        array_list = [[float(x) for x in line.split() if '\x00' not in x ] for line in f]
    sizeS = np.asarray([x[0] for x in array_list if len(x) >= 1])
    
    # Read in from rg_s file
    with open("sp" + str(seed) + '_rg' + filetype + '.txt', "r") as f:
        header=[x for x in next(f)]
        params = [float(x) for x in next(f).split()]
        array_list = [[float(x) for x in line.split() if '\x00' not in x ] for line in f]
    rg = np.asarray([x[0] for x in array_list if len(x) >= 1])
    
    # Get the dissipation 
    params[1] = int(params[1])
    
    # Weird discrepencies between length of the two files. So get whichever list is smallest 
    if len(rg) > len(sizeS): 
        return sizeS, rg[:len(sizeS)], params
    else: 
        return sizeS[:len(rg)], rg, params


def calc_xi_disp(seed_list, filetype):
    """ 
    Compute the correlation length from the radius of Gyration
    
    Input: 
        - seed for run (contained in filename)
    
    Return: 
        - lambda_s: (type:List) of the dissipations, lambda
        - xi_s: the correlation length (squared)
    """
    xi_s = np.zeros_like(seed_list);  lambda_s = np.zeros_like(seed_list).astype('float')
    for s in range(len(seed_list)):
        sizeS, rg, params = readRG(seed_list[s],filetype)
        lambda_s[s] = params[2]
        xi_s[s] = np.sum(sizeS * rg * rg) / np.sum(sizeS)
    return lambda_s, np.sqrt(xi_s)


# ===============================================================
filetypes = ['s', 'a']
heightRange = np.round(2**(np.unique(np.arange(8,17, 0.5))))
lambdaRange = 10**np.arange(-6, -1, 0.25).astype(float)

setRange = np.zeros((len(heightRange), len(lambdaRange),2))
for j in range(len(heightRange)): 
    for k in range(len(lambdaRange)): 
        setRange[j,k,0] = heightRange[j]
        setRange[j,k,1] = lambdaRange[k]

fixed_SetRange = setRange.transpose(2,0,1).reshape(2,-1)

# ===============================================================
for hc in heightRange: 
    for filetype in filetypes: 
        indexRange = []
        for j in range(len(fixed_SetRange[0,:])): 
            if fixed_SetRange[0,j] == hc:
                filename = "sp" + str(j+1000) + "_rg" + filetype + ".txt"
                if path.exists(filename):
                    indexRange += [j+1000] 
                    
        print(filetype, hc, indexRange[0], indexRange[-1], end=' ')
        
        lambda_x, xi_x = calc_xi_disp(indexRange, filetype)
    
        results = np.stack((lambda_x, xi_x))
    
        filename = "Histo_fortran/xi_"+filetype+str(hc)+".txt"
        np.savetxt(filename, results, delimiter=',')
    
    
    