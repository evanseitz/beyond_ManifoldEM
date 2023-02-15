import sys, os
sys.dont_write_bytecode = True
import numpy as np
from numpy import linalg as LA
from itertools import permutations, combinations
from scipy.spatial.distance import pdist, cdist, squareform
import pandas as pd
import matplotlib
from matplotlib import rc
matplotlib.rc('text', usetex = True)
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from pylab import imshow, show, loadtxt, axes
import scipy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import cm
import GaussianBandwidth

pyDir = os.path.dirname(os.path.abspath(__file__)) #python file location

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

#################
# generate data:

Dist0 = np.load('EDM_SS2_dist.npy')
m = 400 #number of states

Dist = Dist0[0:m,0:m]*(250**3)

# =============================================================================
# EDM plot used in Figure 8
# =============================================================================
if 1: #plot distances of state_01_01 to all others
    fig, ax = plt.subplots()
    ax.scatter(np.linspace(1,m,m), Dist[0,:], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
    pointsX1 = []
    pointsY1 = []
    pointsX2 = []
    pointsY2 = []
    for i in range(0,400):
        if i%20 == 0:
            pointsX1.append(i+1)
            pointsY1.append(Dist[0,:][i])
        if i < 20:
            pointsX2.append(i+1)
            pointsY2.append(Dist[0,:][i])           
        
    ax.plot(pointsX1, pointsY1, c='red', linewidth=1.5, zorder=-1, label='$\mathrm{CM_1}$')
    ax.plot(pointsX2, pointsY2, c='blue', linewidth=1.5, zorder=-1, label='$\mathrm{CM_2}$')
    ax.legend(loc='lower right', fontsize=14, frameon=False)
    ax.set_xlabel(r'States', fontsize=20, labelpad=12)
    ax.set_ylabel(r'$D_{1,k}$', fontsize=22, labelpad=10)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    #plt.show()
    fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig8_EDM_dists.pdf'), dpi=600)
    plt.clf()

# ========================================================================
# Generate CM2 ground-truth indexing:
# ======================================================================== 
if 1: 
    # CM1:
    CM1_idx = np.arange(0,m)
    
    # CM2:
    binsActual = []
    Idx = 0
    tau = 1
    Nx = 20
    Ny = 20
    CM2_idx = np.ndarray(shape=(Nx, Nx*tau), dtype=int)  
    for s in range(0,Nx):
        state_list = []
        for r in range(0,Ny):
            state_list.append(np.arange(Idx,Idx+tau))
            Idx+=(Nx*tau)
        binsActual.append([item for sublist in state_list for item in sublist])
        Idx-=(Nx*Ny-tau)
    
    for i in range(Nx):
        for j in range(Ny*tau):
            CM2_idx[i,j] = binsActual[i][j]     


# =============================================================================
# # generate optimal kernel:
# =============================================================================
if 0: #automated approach
    logEps = np.arange(-100, 100.2, 0.2) #may need to widen range if a hyperbolic tangent is not captured
    a0 = 1*(np.random.rand(4,1)-.5)
    popt, logSumWij, resnorm, R2 = GaussianBandwidth.op(Dist, logEps, a0)
    
    def fun(xx, aa0, aa1, aa2, aa3): #fit tanh()
        F = aa3 + aa2 * np.tanh(aa0 * xx + aa1)
        return F
    
    plt.scatter(logEps, logSumWij, s=1, c='C0', edgecolor='#1f77b4', zorder=.1, label='data')
    plt.plot(logEps, fun(logEps, popt[0], popt[1], popt[2], popt[3]), c='C1', linewidth=.5, zorder=.2, label='tanh(x)')
    plt.plot(logEps, popt[0]*popt[2]*(logEps+popt[1]/popt[0]) + popt[3], c='C2', linewidth=.5, zorder=.3, label='slope')
    plt.axvline(x=-(popt[1] / popt[0]), c='C3', linewidth=.5, label='epsilon')
    
    plt.legend(loc='best')
    plt.xlabel(r'$\mathrm{ln \ \epsilon}$', fontsize=16)
    plt.ylabel(r'$\mathrm{ln \ \sum_{i,j} \ A_{i,j}}$', fontsize=18, rotation=90)
    plt.ylim(np.amin(fun(logEps, popt[0], popt[1], popt[2], popt[3]))-1, np.amax(fun(logEps, popt[0], popt[1], popt[2], popt[3]))+1)

    slope = popt[0]*popt[2] #slope of tanh
    ln_eps = -(popt[1] / popt[0]) #x-axis line through center of tanh
    eps = np.exp(ln_eps)
    print('Coefficient of Determination: %s' % R2)
    print('Slope: %s' % slope) 
    print('ln(epsilon): %s; epsilon: %s' % (ln_eps, eps))
    if 0: #plot Gaussian Bandwidth (graph should be a hyperbolic tangent)
        plt.show()
    plt.clf()
    BigEps = False
    
else:
    eps = 2000000
    BigEps = False


# =========================================================================
# Generate optimal Gaussian kernel for Similarity Matrix (A)
# =========================================================================
A = np.exp(-1.*((Dist**2.) / (2.*eps))) #similarity matrix

if 0: #plot similarity matrix A
    imshow(A, cmap='jet', origin='lower')
    plt.title(r'Gaussian Kernel, $\mathit{\epsilon}$=%s' % eps, fontsize=20)
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    if 0:
        rowSums = np.sum(A, axis=1)
        print('minimum affinity:', np.where(rowSums == np.amin(rowSums)))
        plt.scatter(np.linspace(1,m,m), rowSums)
        plt.xlim(1, m)
        plt.ylim(np.amin(rowSums), np.amax(rowSums))
        plt.show()
    if 0: #similarity of state_01_01 to all others
        plt.scatter(np.linspace(1,m,m), A[0,:])
        plt.show()

# =========================================================================
# Construction of Markov Transition Matrix (M):        
# =========================================================================
D = np.ndarray(shape=(m,m), dtype=float) #diagonal matrix D
for i in range(0,m):
    for j in range(0,m):
        if i == j:
            D[i,j] = np.sum(A[i], axis=0)
        else:
            D[i,j] = 0
            
Dinv = np.linalg.inv(D)
M = np.matmul(A, Dinv) #Markov matrix via normalization of A to right stochastic form
if 0: #check M is right (row) stochastic
    print(np.sum(M, axis=0)) #should be all 1's
    
# =========================================================================
# cite: 'Systematic determination of order parameters for chain dynamics...
#       ...using diffusion maps'; PNAS, Ferguson (2010), SI
# =========================================================================
Dhalf = scipy.linalg.sqrtm(D)
Dinv_half = scipy.linalg.sqrtm(Dinv)
if 0: #sanity-check for D^-1/2 
    print(Dinv_half.dot(D).dot(Dinv_half)) #should be Identity matrix
Ms = Dinv_half.dot(M).dot(Dhalf) #M is adjoint to symmetric matrix Ms
    
def is_symmetric(A):
    # for positive semidefinite, need to additionally check that all...
    # ...eigenvalues are (significantly) non-negative
    if np.allclose(A, A.T, rtol=1e-08, atol=1e-08):
        print('Matrix is symmetric')    
    else:
        print('Matrix is not symmetric')
is_symmetric(Ms) #check if matrix is symmetric

# =========================================================================
# Eigendecomposition
# =========================================================================
U, sdiag, vh = np.linalg.svd(Ms) #vh = U.T
sdiag = sdiag**(2.) #eigenvalues given by s**2    


# =========================================================================
# Analysis of diffusion map
# =========================================================================
s = 20
lw = .5
enum = np.arange(1,(m)+1)

if 0: #orthogonality check
    print(np.linalg.det(U)) #should be +- 1
    print(np.sum(U[:,0]*U[:,0])) #should be 1
    print(np.sum(U[:,1]*U[:,1])) #should be 1
    print(np.sum(U[:,1]*U[:,2])) #should be ~0
    print(np.sum(U[:,1]*U[:,3])) #should be ~0
    print(np.sum(U[:,2]*U[:,3])) #should be ~0

if 0: #plot eignevalue spectrum
    x = range(1,len(sdiag[1:])+1)
    plt.scatter(x, sdiag[1:])
    plt.title('Eigenvalue Spectrum, $\mathit{\epsilon}$=%s' % eps)
    plt.xlabel(r'$\mathrm{\Psi}$')
    plt.ylabel(r'$\mathrm{\lambda}$', rotation=0)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlim([0,15])
    plt.ylim(sdiag[15], sdiag[1])
    plt.locator_params(nbins=15)
    plt.axhline(y=0, color='k', alpha=.5, linestyle='--', linewidth=1)
    plt.show()
    
if 0: #eigenfunction analysis
    v1,v2 = 1,2 #note to avoid steady-state
    plt.subplot(2,1,1)
    plt.scatter(enum, U[:,v1-1], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
    plt.ylabel(r'$\psi_{%s}$' % (v1-1))
    plt.gca().set_box_aspect(1)
    plt.subplot(2,1,2)
    plt.scatter(enum, U[:,v2-1], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
    plt.xlabel(r'$\mathrm{CM_{1} \ Index}$')
    plt.ylabel(r'$\psi_{%s}$' % (v2-1))
    plt.gca().set_box_aspect(1)
    plt.show()
    plt.clf()
    
    
# =============================================================================
# All plots used in Figure 7
# =============================================================================
if 1: #ordered array 2D manifold subspaces
    dimRows = 2#4
    dimCols = 10#6

    # hardcoded coordinates for figure subplots
    V1_top = [1,1,1,1,1,1]
    V2_top = [2,3,4,5,6,7]
    V1_bot = [2,2,2,2,2,2]
    V2_bot = [3,4,5,6,7,8]
    V_idx = 0
            
    widths = []
    for i in range(dimCols):
        widths.append(1)
    widths.append(.2)
    
    fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
    gs = GridSpec(dimRows, dimCols+1, left=0.05, right=0.95, wspace=0.3, hspace=0.3, width_ratios=widths)
    
    idx_col1 = 0  
    idx_col2 = 5
    for v1 in range(1,6):#(2*int(dimCols/2)+1)):

        # CM1 side:
        ax = fig.add_subplot(gs[0, idx_col1])
        ax.scatter(list(CM1_idx.flatten()), U[:,v1], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
        ax.set_xlabel(r'$\mathrm{CM}_{1}$', fontsize=14, labelpad=2.5)
        ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

        # CM2 side:
        ax = fig.add_subplot(gs[1, idx_col1])
        ax.scatter(list(CM2_idx.flatten()), U[:,v1], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
        ax.set_xlabel(r'$\mathrm{CM}_{2}$', fontsize=14, labelpad=2.5)
        ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Cartesian-product side:
        ax = fig.add_subplot(gs[0, idx_col2])
        ax.scatter(U[:,V1_top[V_idx]], U[:,V2_top[V_idx]], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
        ax.set_xlabel(r'$\Psi_{%s}$' % V1_top[V_idx], fontsize=18, labelpad=2.5)
        ax.set_ylabel(r'$\Psi_{%s}$' % V2_top[V_idx], fontsize=18, labelpad=2.5)
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

        ax = fig.add_subplot(gs[1, idx_col2])
        ax.scatter(U[:,V1_bot[V_idx]], U[:,V2_bot[V_idx]], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
        ax.set_xlabel(r'$\Psi_{%s}$' % V1_bot[V_idx], fontsize=18, labelpad=2.5)
        ax.set_ylabel(r'$\Psi_{%s}$' % V2_bot[V_idx], fontsize=18, labelpad=2.5)
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        
        idx_col1 += 1
        idx_col2 += 1
        V_idx += 1

    ax = fig.add_subplot(gs[:, dimCols])
    
    sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cm.tab20).reversed())
    sm.set_clim(vmin=0, vmax=U.shape[0])
    cbar = fig.colorbar(sm, cax=ax)
    cbar.ax.tick_params(labelsize=16)

    fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig7.pdf'), dpi=600)
    plt.clf()
    
    