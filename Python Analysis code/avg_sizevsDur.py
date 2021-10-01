from os import path
import numpy as np
import matplotlib as mpl

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
    
def findSizeVsDur(seed):
    size, params = readSizes(seed,'s')
    duration, params = readSizes(seed,'d')
        
    # Create arrays for durations and corresponding average size
    size_avg = np.arange(np.min(duration),np.max(duration)+1) * 0
    size_counts = np.arange(np.min(duration),np.max(duration)+1) * 0
    duration_each = np.arange(np.min(duration),np.max(duration)+1)

    end_iter = min(len(size), len(duration))
    for i in range(end_iter):      # count up the avergae size for each given duration. 
        size_avg[duration[i]-1] += size[i]
        size_counts[duration[i]-1] += 1

    for i in range(len(size_avg)-1,-1,-1):      # Delete any 0 counts. 
        if size_avg[i]==0:
            size_avg = np.delete(size_avg,i)
            size_counts = np.delete(size_counts,i)
            duration_each = np.delete(duration_each,i)
    size_avg = np.divide(size_avg, size_counts)
    return size_avg, duration_each, params

def findSizeAVsDur(seed):
    size, params = readSizes(seed,'a')
    duration, params = readSizes(seed,'d')
        
    # Create arrays for durations and corresponding average size
    size_avg = np.arange(np.min(duration),np.max(duration)+1) * 0
    size_counts = np.arange(np.min(duration),np.max(duration)+1) * 0
    duration_each = np.arange(np.min(duration),np.max(duration)+1)

    end_iter = min(len(size), len(duration))
    for i in range(end_iter):      # count up the avergae size for each given duration. 
        size_avg[duration[i]-1] += size[i]
        size_counts[duration[i]-1] += 1

    for i in range(len(size_avg)-1,-1,-1):      # Delete any 0 counts. 
        if size_avg[i]==0:
            size_avg = np.delete(size_avg,i)
            size_counts = np.delete(size_counts,i)
            duration_each = np.delete(duration_each,i)
    size_avg = np.divide(size_avg, size_counts)
    return size_avg, duration_each, params

def plotSizeVsDur(seed):
    print(str(seed)+"s", end=" ")
    y, x, params = findSizeVsDur(seed)
    z = np.column_stack((x, y))
    ## Write params to the file, commented with '#'
    filename = "Histo_fortran/dur_size" + str(seed) + ".txt"
    np.savetxt(filename, z, delimiter=',')
    
    print(str(seed)+"a", end=" ")
    ya, x, params = findSizeAVsDur(seed)
    za = np.column_stack((x, ya))
    ## Write params to the file, commented with '#'
    filename = "Histo_fortran/dur_area" + str(seed) + ".txt"
    np.savetxt(filename, za, delimiter=',')


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
print("starting")
for i in range(len(fixed_SetRange[0,:])): 
    seed = i + 1000
    filename = "sp" + str(seed) + "_" + "s" + ".txt"
    if path.exists(filename): 
        plotSizeVsDur(seed)