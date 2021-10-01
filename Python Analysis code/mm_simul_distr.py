 # Load necessary packages
# %pylab notebook
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from os import path

try:
    plt.style.use('classic')
except:
    pass

# ===========================================================================================================

def driving_m(lattice, L, hc, lambd):
    """
    Randomly choose sites on the lattice to add a height-unit to 
    until one site overcomes the threhsold, hc
    """
    while np.any(lattice < hc): 
        ##Choose a ranom site (assuming 2d lattice)
        site_x = int(L*random.random())
        site_y = int(L*random.random())
        lattice[site_x, site_y] += 1
        
        ##If incrimented site height >= hc, proceed to the toppling rules: 
        if lattice[site_x, site_y] >= hc:
            lattice, [siz, dur, area, rgs, rga] = relax_m(lattice, L, hc, lambd)
            break
            
    return lattice, [siz, dur, area, rgs, rga]


def relax_m(lattice, L, hc, lambd):
    """
    Relaxes the (active) lattice to a stable configuration. 
    Returns: avalanche size, area, duration, and corresponding radii of gyration. 
    """
    siz = 0
    dur = 0
    ##Temp Lattice with extra sites for the open boundaries (ease of toppling code). 
    latticeOBC = np.zeros((L+2, L+2))
    latticeOBC[1:-1, 1:-1] = lattice
    ##Lattice-clone to keep track of sites that fired at-least once (for computing area). 
    latticeTriggers = np.zeros_like(lattice)
    
    ## Do topplings
    while(np.any(lattice >= hc)):
        dur += 1
        ##Loop over lattice, checking for sites >= hc that need to be toppled. 
        for x in range( 1, L+1 ):
            for y in range( 1, L+1 ):
                if lattice[x-1,y-1] >= hc:
                    ##Draw final configuartion and set this configuration to the (temp) lattice. 
                    final_configs = np.random.multinomial(hc, [(1-lambd)/4.]*4 + [lambd], size=1)
                    final_config = final_configs[0]
                    latticeOBC[x,y] -= hc
                    latticeOBC[x+1,y] += final_config[0]
                    latticeOBC[x-1,y] += final_config[1]
                    latticeOBC[x,y+1] += final_config[2]
                    latticeOBC[x,y-1] += final_config[3]
                    ##Also record the site that triggered
                    latticeTriggers[x-1,y-1] += 1
                    
        ##Reset the boundaries to be 0 (so they never get too big)
        latticeOBC[:,0] = 0; latticeOBC[:,-1] = 0
        latticeOBC[0,:] = 0; latticeOBC[-1,:] = 0
        ## Set the new (temp) lattice onto the real (recorded) lattice. 
        lattice = latticeOBC[1:-1, 1:-1]
    siz = np.sum(latticeTriggers) #Size = total number of sites that toppled
    area = np.count_nonzero(latticeTriggers) #Area = total number of sites that toppled at least once. 
    
    # Compute the radius of Gyration (s/a)
    rsx=0; rsy=0; rax=0; ray = 0; rgs = 0; rga = 0
    for x in range(0, L):
        for y in range(0,L):
            rsx += x * latticeTriggers[x,y] / siz
            rsy += y * latticeTriggers[x,y] / siz
            rax += x * latticeTriggers[x,y] / (latticeTriggers[x,y] - 0.001) / area
            ray += y * latticeTriggers[x,y] / (latticeTriggers[x,y] - 0.001) / area
    for x in range(0, L):
        for y in range(0,L):
            rgs += latticeTriggers[x,y]/siz * ( (x-rsx)**2 + (y - rsy)**2)
            rga += latticeTriggers[x,y] / (latticeTriggers[x,y] - 0.001) / area * ( (x-rax)**2 + (y - ray)**2)            
        
    return lattice, [siz, dur, area, rgs, rga]

# ===========================================================================================================

def runSystem_m(lattice, L, hc, lambd,timer):
    
    lattice = lattice.astype(int)
    sizList = []
    durList = []
    areList = []
    rgsList = []
    rgaList = []
    HList = []
    counter = 0
    while counter < timer: 
        lattice, results = driving_m(lattice, L, hc, lambd)
        
        if counter%10==0:
            print('reached driving: '+str(counter+1), results)

        if results[0]>0:
            sizList += [results[0]]
            durList += [results[1]]
            areList += [results[2]]
            rgsList += [results[3]]
            rgaList += [results[4]]
            counter += 1        
            
    return sizList,durList,areList, rgsList, rgaList
    
# ===========================================================================================================

## Running a bunch of simulations... 

params = sys.argv

L = 512
hc = int(params[1])
lamb = float(params[2])
timer = 5000000
lattice = (hc-1)*np.random.rand(L,L)

sizList0 , durList0 , areList0 , rgsList0 , rgaList0 = runSystem_m(lattice, L, hc, lamb, timer)

## Write the observables text file
jobname = 'mm_' + str(hc) + "_"+ str(lamb)
listToStr = ','.join([str(elem) for elem in [L,hc,lamb,timer]])

filename = "results/"+jobname+"_s.txt"
np.savetxt(filename, sizList0, delimiter=',', header=listToStr, comments='#')

filename = "results/"+jobname+"_d.txt"
np.savetxt(filename, durList0, delimiter=',', header=listToStr, comments='#')

filename = "results/"+jobname+"_a.txt"
np.savetxt(filename, areList0, delimiter=',', header=listToStr, comments='#')

filename = "results/"+jobname+"_rs.txt"
np.savetxt(filename, rgsList0, delimiter=',', header=listToStr, comments='#')

filename = "results/"+jobname+"_ra.txt"
np.savetxt(filename, rgaList0, delimiter=',', header=listToStr, comments='#')
