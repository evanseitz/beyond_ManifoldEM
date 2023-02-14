import sys, os
import numpy as np
from numpy import linalg as LA
from itertools import permutations, combinations
from scipy.spatial.distance import pdist, cdist, squareform
import pandas as pd
import matplotlib
from matplotlib import rc
matplotlib.rc('text', usetex = True)
#matplotlib.use('Qt5Agg')
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from pylab import imshow, show, loadtxt, axes
import scipy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
#external scripts:
import fergusonE

pyDir = os.path.dirname(os.path.abspath(__file__)) #python file location

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

#################
# generate data:

Dist = np.load('RMSD_Dist_All.npy')
m = 400 #number of PDBs

if 1: #plot distances of state_01_01 to all others
    plt.scatter(np.linspace(1,m,m), Dist[0,:], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
    point1 = [1, Dist[0,:][0]]
    point2 = [20, Dist[0,:][19]]
    x_values1 = [point1[0], point2[0]]
    y_values1 = [point1[1], point2[1]]
    point3 = [1, Dist[0,:][0]]
    point4 = [381, Dist[0,:][380]]
    x_values2 = [point3[0], point4[0]]
    y_values2 = [point3[1], point4[1]]
    plt.plot(x_values2, y_values2, c='red', linewidth=1.5, zorder=-1, label='$\mathrm{CM_1}$')
    plt.plot(x_values1, y_values1, c='blue', linewidth=1.5, zorder=-1, label='$\mathrm{CM_2}$')

    plt.legend(loc='lower right')
    plt.xlabel(r'States', fontsize=14, labelpad=7.5)
    plt.ylabel(r'$D_{1,k}$', fontsize=15, labelpad=10)
    plt.show()


#####################################
# perform search for optimal epsilon:

'''if 0:
    logEps = np.arange(-10,10.2,0.2)
    a0 = 1*(np.random.rand(4,1)-.5)
    popt, logSumWij, resnorm = fergusonE.op(RMSD,logEps,a0)

    def fun(xx, aa0, aa1, aa2, aa3): #fit tanh()
        F = aa3 + aa2 * np.tanh(aa0 * xx + aa1)
        return F
    
    plt.scatter(logEps, logSumWij, s=1, c='#1f77b4', edgecolor='#1f77b4', zorder=.1)
    plt.plot(logEps, fun(logEps, popt[0], popt[1], popt[2],popt[3]), c='red', linewidth=.5, zorder=.2)    
    #ax.axvline(-(popt[1] / popt[0]), color='lightgray', alpha=.5, linewidth=1, linestyle='-', zorder=0)
    print(-(popt[1] / popt[0]))
    plt.xlabel(r'$\mathrm{ln \ \epsilon}$', fontsize=20)
    plt.ylabel(r'$\mathrm{ln \ \sum_{i,j} \ A_{i,j}}$', fontsize=20, rotation=90)
    plt.show()'''

###########################
# generate optimal kernel:
    
#use 0.1 for cosines, 10000 for rectilinear
if 0:
    eps = 1000
    BigEps = True
else:
    eps = .12
    BigEps = False
#A = np.exp(-(RMSD**2 / 2*eps0)) #distance matrix


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
# High-dimensional counter rotation (if unaligned from preferred eigenbasis)
# =========================================================================  
'''dim = 20 #analyze 'dim'-dimensional manifold subspace with 'dim'('dim'-1)/2 rotation operators
if 0:
    theta_total = int(dim*(dim-1)/2) #total number of rotation operators
    U_init = U[:,1:dim+1]
    thetas = np.zeros(shape=(theta_total,1), dtype=float)
    thetas[0] = -19*np.pi/180 #Rij operator
    #thetas[theta_total-1] = 45*np.pi/180 #Rij operator
    thetas[((dim-1)+(dim-2)+(dim-3)+1)-1] = 45*np.pi/180 #Rij operator
    R = Generate_Nd_Rot.genNdRotations(dim, thetas)
    U_rot = np.matmul(R, U_init.T)
    U = U_rot.T
    rotated = True
else:
    rotated = False'''

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
    #print(np.amin(U[:,v1-1]),np.amax(U[:,v1-1]))
    #print(np.amin(U[:,v2-1]),np.amax(U[:,v2-1]))
    #plt.xlabel(r'$\mathrm{CM_{1} \ Index}$')
    plt.ylabel(r'$\psi_{%s}$' % (v1-1))
    plt.gca().set_box_aspect(1)
    plt.subplot(2,1,2)
    plt.scatter(enum, U[:,v2-1], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
    plt.xlabel(r'$\mathrm{CM_{1} \ Index}$')
    plt.ylabel(r'$\psi_{%s}$' % (v2-1))
    plt.gca().set_box_aspect(1)
    plt.show()
    plt.clf()
            
if 1: #view an organized array of (1D) eigenfunctions  
    fig = plt.figure()
    dimRowsTotal = 4
    dimCols = 6
    idx = 1
    fS = 12 #font size
    for v1 in range(1,2*dimCols+1):
        plt.subplot(dimRowsTotal, dimCols, idx)
        plt.scatter(enum, U[:,v1], s=s, c='white', edgecolor='k', linewidths=lw, zorder=1)
        #plt.scatter(enum, U[:,v1][Y_idx], s=s, c='white', edgecolor='k', linewidths=lw, zorder=1)
        #plt.plot(enum, U[:,v1], color='lightgray', linewidth=1, zorder=-1)
        plt.xlabel(r'$\mathrm{CM_{1}}$', fontsize=fS-.5)
        plt.ylabel(r'$\Psi_{%s}$' % (v1), fontsize=fS+2, labelpad=2.5)
        
        if BigEps is True:
            #plt.xlim(-1,1)
            if np.abs(np.amax(U[:,v1])) >= np.abs(np.amin(U[:,v1])) :
                plt.ylim(-1*np.amax(U[:,v1]),np.amax(U[:,v1]))
            else:
                plt.ylim(np.amin(U[:,v1]),-1*np.amin(U[:,v1]))
        
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.rc('font', size=6)
        if 1:
            frame = plt.gca()
            frame.axes.xaxis.set_ticklabels([])
            frame.axes.yaxis.set_ticklabels([])
            plt.gca().set_xticks([])
            plt.xticks([])
            plt.gca().set_yticks([])
            plt.yticks([])
        else:
            plt.tick_params(axis="x", labelsize=6)
            plt.tick_params(axis="y", labelsize=6)
        plt.gca().set_box_aspect(1)
        idx += 1
    #plt.tight_layout()
    #plt.show()

#if 1: #ordered array 2D manifold subspaces    
    dimRows = 2
    #dimCols = 6
    idx = 2*dimCols+1#1
    for v1 in range(1,dimRows+1):
        for v2 in range(v1+1, v1+dimCols+1):
            #print(v1,v2)
            plt.subplot(dimRowsTotal, dimCols, idx)
            if 0:#if groundTruth is True:
                color=iter(cm.tab20(np.linspace(1, 0, np.shape(CM_idx)[0])))
                for b in range(np.shape(CM_idx)[0]):
                    c=next(color)
                    plt.scatter(U[:,v1][CM_idx[b]], U[:,v2][CM_idx[b]], color=c, s=15, edgecolor='k', linewidths=.1, zorder=1)
                #plt.scatter(U[:,v1], U[:,v2], c=np.arange(1,m+1), cmap=cmap, s=s, linewidths=lw, edgecolor='k') #cmap: 'nipy_spectral', 'gist_rainbow'
            else:
                plt.scatter(U[:,v1], U[:,v2], c='white', s=20, linewidths=.5, edgecolor='k', zorder=0)
            plt.xlabel(r'$\Psi_{%s}$' % (v1), fontsize=fS+2, labelpad=2.5)
            plt.ylabel(r'$\Psi_{%s}$' % (v2), fontsize=fS+2, labelpad=2.5)
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 
            plt.rc('font', size=6)
            if 1:
                frame = plt.gca()
                frame.axes.xaxis.set_ticklabels([])
                frame.axes.yaxis.set_ticklabels([])
                plt.gca().set_xticks([])
                plt.xticks([])
                plt.gca().set_yticks([])
                plt.yticks([])
            else:
                plt.tick_params(axis="x", labelsize=6)
                plt.tick_params(axis="y", labelsize=6) 

            #plt.xlim(np.amin(U[:,v1])*1.1, np.amax(U[:,v1])*1.1)
            #plt.ylim(np.amin(U[:,v2])*1.1, np.amax(U[:,v2])*1.1)
            plt.gca().set_box_aspect(1)

            idx += 1 
    plt.tight_layout()
    plt.subplots_adjust(left=0.01, right=0.45, bottom=0.5, top=0.99, wspace=0.26, hspace=0.23)
    plt.show()
    #fig = plt.gcf()
    #fig.savefig(os.path.join(outDir,'Fig_PD_%s.png' % PD), dpi=200)
    #plt.clf()

'''A = np.exp(-1.*((RMSD**2.) / (2.*eps0))) #similarity matrix

alpha = 1
#alpha = 1.0: Laplace-Beltrami operator
#alpha = 0.5: Fokker-Planck diffusion
#alpha = 0.0: graph Laplacian normalization

if 0:
    imshow(A, cmap='jet', origin='lower')
    plt.title(r'Gaussian Kernel, $\mathit{\epsilon}$=%s' % eps0, fontsize=20)
    plt.colorbar()
    plt.tight_layout()
    show()
    

###########################################
# normalized graph laplacian construction:

D = np.ndarray(shape=(m,m), dtype=float) #diagonal matrix
for i in range(0,m):
    for j in range(0,m):
        if i == j:
            D[i,j] = np.sum(A[i], axis=0)
        else:
            D[i,j] = 0    

if 0:
    print(D.diagonal())
    imshow(D, cmap='jet', origin='upper')
    plt.title(r'$\mathit{M}, \mathit{\alpha}$=%s' % alpha, fontsize=20)
    plt.colorbar()
    plt.tight_layout()
    show()
    
D_inv = scipy.linalg.fractional_matrix_power(D, -1)
M = np.matmul(A, D_inv)

if 0:
    print(np.sum(M, axis=0)) #should be all 1's
    imshow(M, cmap='jet', origin='lower')
    plt.title('Markov Transition Matrix', fontsize=20)
    plt.colorbar()
    plt.tight_layout()
    show()

#####################
# SVD decompisition:
    
def tidyUp(D,EV):
    order = np.argsort(D)[::-1]
    D = np.sort(D)[::-1]
    EV = EV[:,order]
    sqrtD = np.sqrt(D)
    S = np.diag(sqrtD)
    invS = np.diag(1./sqrtD)
    return (D,EV,S,invS)

D,U = np.linalg.eigh(np.matmul(M,M.T))
D,U,S,invS = tidyUp(D,U)
V = np.matmul(M.T,np.matmul(U,invS))
sdiag = np.diag(S)

if 0:
    if 0:
        print('U:', U)
        print('S:', sdiag)
        print('V:', V)
        print(sdiag[0], sdiag[1])
        print(U[0], U[1])
    if 1: #eignevalue spectrum
        x = range(1,len(sdiag[1:])+1)
        print(sdiag[1:])
        plt.scatter(x, sdiag[1:])
        plt.title('Eigenvalue Spectrum, $\mathit{\epsilon}$=%s' % eps0, fontsize=20)
        plt.xlabel(r'$\mathrm{\Psi}$')
        plt.ylabel(r'$\mathrm{\lambda}$', rotation=0)
        plt.xlim([0,15])
        plt.locator_params(nbins=15)
        plt.axhline(y=0, color='k', alpha=.5, linestyle='--', linewidth=1)
        plt.show()
    
if 1:
    if 1: #2d diffusion map
        #plt.scatter(U[:,1]*sdiag[1], U[:,2]*sdiag[2], c=U[:,1]*sdiag[1])
        plt.scatter(U[:,1], U[:,2], c=U[:,1])
        plt.title(r'Diffusion Map, $\mathit{\epsilon}$=%s' % eps0, fontsize=20)
        plt.xlabel(r'$\psi_1$')
        plt.ylabel(r'$\psi_2$')
        plt.colorbar()
        plt.show()
        
    if 1: #eigenfunction analysis
        plt.scatter(np.arange(1,m+1), U[:,1], c=U[:,1], cmap='gist_rainbow', s=15, edgecolor='k', linewidths=.1, zorder=1)
        plt.title(r'Diffusion Map, $\mathit{\epsilon}$=%s' % eps0, fontsize=20)
        plt.xlabel(r'$CM_1 Index$')
        plt.ylabel(r'$\psi_1$')
        plt.colorbar()
        plt.show()
        
    if 1: #ordered array 2D manifold subspaces    
        dimRows = 4
        dimCols = 6
        idx = 1
        for v1 in range(1,dimRows+1):
            for v2 in range(v1+1, v1+dimCols+1):
                plt.subplot(dimRows, dimCols, idx)
               # color=iter(cm.tab20(np.linspace(1, 0, np.shape(CM_idx)[0])))
                #for b in range(np.shape(CM_idx)[0]):
                    #c=next(color)
                plt.scatter(U[:,v1], U[:,v2], c=np.arange(1,m+1), cmap='gist_rainbow', s=15, edgecolor='k', linewidths=.1, zorder=1)
                #plt.scatter(U[:,v1], U[:,v2], c=np.arange(1,m+1), cmap=cmap, s=s, linewidths=lw, edgecolor='k') #cmap: 'nipy_spectral', 'gist_rainbow'
                plt.xlabel(r'$\Psi_{%s}$' % v1, fontsize=6, labelpad=2.5)
                plt.ylabel(r'$\Psi_{%s}$' % v2, fontsize=6, labelpad=2.5)
                plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 
                plt.rc('font', size=6)
                if 1:
                    frame = plt.gca()
                    frame.axes.xaxis.set_ticklabels([])
                    frame.axes.yaxis.set_ticklabels([])
                    plt.gca().set_xticks([])
                    plt.xticks([])
                    plt.gca().set_yticks([])
                    plt.yticks([])
                else:
                    plt.tick_params(axis="x", labelsize=6)
                    plt.tick_params(axis="y", labelsize=6) 
    
                plt.xlim(np.amin(U[:,v1])*1.1, np.amax(U[:,v1])*1.1)
                plt.ylim(np.amin(U[:,v2])*1.1, np.amax(U[:,v2])*1.1)
                idx += 1 
        plt.tight_layout()
        plt.subplots_adjust(left=0.02, right=0.99, bottom=0.05, top=0.99, wspace=0.26, hspace=0.23)
        plt.show()
        #fig = plt.gcf()
        #fig.savefig(os.path.join(outDir,'Fig_PD_%s.png' % PD), dpi=200)

    if 0: #3d diffusion map
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(U[:,2], U[:,1], U[:,3], c=U[:,1], cmap='gist_rainbow', s=15)#, edgecolor='k', linewidths=.1, zorder=1)
        ax.view_init(elev=15, azim=-60)
        #plt.title('3D Embedding')
        ax.set_xlabel(r'$\psi_2$')
        ax.set_ylabel(r'$\psi_1$')
        ax.set_zlabel(r'$\psi_3$')
        plt.show()'''
        
    