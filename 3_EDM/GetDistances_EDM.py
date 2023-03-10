import sys, os, re
import numpy as np
from numpy import linalg as LA
import mrcfile
import matplotlib
from matplotlib import rc
#matplotlib.rc('text', usetex = True)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from pylab import imshow, show, loadtxt, axes
import scipy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator

# =============================================================================
# Python EDM Distance Calculator (run via 'python GetDistances_EDM.py') 
# Author:    E. Seitz @ Columbia University - Frank Lab - 2019-2020
# Contact:   evan.e.seitz@gmail.com
# =============================================================================

def calc_dist(vol1, vol2):
    dist = LA.norm(vol1 - vol2)
    #return np.divide(dist, np.shape(vol1)[0]**3)
    return dist
    
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

pyDir = os.getcwd() #python file location
parDir = os.path.abspath(os.path.join(pyDir, os.pardir))
dataDir = os.path.join(pyDir, 'SS2_MRC_Pristine') #location of all MRCs (i.e., entire state space)

MRC_paths0 = []
for root, dirs, files in os.walk(dataDir):
    for file in sorted(files):
        if not file.startswith('.'): #ignore hidden files
            if file.endswith(".mrc"):
                MRC_paths0.append(os.path.join(root, file))               
MRC_paths = natural_sort(MRC_paths0)

m = len(MRC_paths) #total number of states
states = range(0,m)
D = np.zeros(shape=(m,m), dtype=float)

for i in states:
    for j in states:
        if i < j:
            print('Index [%s, %s]:' % (i,j))
            vol1 = mrcfile.open(MRC_paths[i])
            vol2 = mrcfile.open(MRC_paths[j])
            D[i,j] = float(calc_dist(vol1.data, vol2.data))
            print('\t%s' % D[i,j])
            vol1.close()
            vol2.close()
            
D = D + D.T

if 1: #save to file
    np.save('EDM_SS2_dist.npy', D)
        
if 1:
    plt.imshow(D, origin='lower', interpolation='nearest', cmap='jet')
    plt.colorbar()
    plt.tight_layout()
    plt.show()