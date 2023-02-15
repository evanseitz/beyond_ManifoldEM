import sys, os
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar
from matplotlib import rc
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import imshow, show
import scipy
import GaussianBandwidth

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

# =============================================================================
# Embed distance files for each PD via diffusion maps framework
# =============================================================================
# SETUP: Before running this script, distance matrices must first be created...
#   ...via the script '1_Dist_Batch'. Next, make sure the input file path matches...
#   ...the one created for your dataset via the 'Dist' variable below.
# RUNNING: To run a series of PDs at once: first edit '2_DM_Batch.sh'...
#   ...for the total number of PDs requested; e.g., {1...5} for 5 PDs...
#   ...then start batch processing via 'sh 2_DM_Batch.sh'
# =============================================================================
# Author:    E. Seitz @ Columbia University - Frank Lab - 2020-2021
# Contact:   evan.e.seitz@gmail.com
# =============================================================================

def op(pyDir, PD):
    parDir1 = os.path.abspath(os.path.join(pyDir, os.pardir))
    parDir2 = os.path.abspath(os.path.join(parDir1, os.pardir))
    
    # =========================================================================
    # User parameters:
    # =========================================================================
    groundTruth = True #use GT indices for visualizations; see '0_Data_Inputs/GroundTruth_Indices'
    viewCM = 1 #{1,2, etc.}; if using ground-truth, CM reference frame to use for color map indices
    SSn = 2 #options: {1 or 2} in 'SS2_PD_Pristine' for 1 or 2 degrees of freedom
    
    # =========================================================================
    # Prepare directories and import files:
    # =========================================================================    
    outDir = os.path.join(pyDir, 'Data_Manifolds')
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    dataDir = os.path.join(pyDir, 'Data_Distances')
    Dist = np.load(os.path.join(dataDir, 'PD_%s_SS%s_dist.npy' % (PD,SSn)))
    
    if groundTruth is True:
        if SSn == 2:
            if viewCM == 1: #view in reference frame of CM1
                CM_idx = np.load(os.path.join(pyDir, 'GroundTruth_Indices/CM1_Indices.npy'), allow_pickle=True)
            elif viewCM == 2: #view in reference frame of CM2
                CM_idx = np.load(os.path.join(pyDir, 'GroundTruth_Indices/CM2_Indices.npy'), allow_pickle=True)
    
            CM1_idx = np.load(os.path.join(pyDir, 'GroundTruth_Indices/CM1_Indices.npy'), allow_pickle=True)
            CM2_idx = np.load(os.path.join(pyDir, 'GroundTruth_Indices/CM2_Indices.npy'), allow_pickle=True)
        elif SSn == 1:
            CM_idx = np.arange(20)
        
    # =========================================================================
    # Distances matrix analysis
    # =========================================================================
    m = np.shape(Dist)[0]
    
    if 0: #plot Distance matrix
        imshow(Dist, cmap='jet', origin='lower', interpolation='nearest')
        plt.title('Distances', fontsize=20)
        plt.xlabel('State')
        plt.ylabel('State')
        plt.colorbar()
        plt.tight_layout()
        plt.show()
        
    if 0: #plot distances of state_01_01 to all others
        plt.scatter(np.linspace(1,m,m), Dist[0,:])
        plt.show()
    
    # =========================================================================
    # Method to estimate optimal epsilon; see Ferguson SI (2010)
    # =========================================================================
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
        if 0: #save Gaussian Bandwidth plot to file
            np.save('GaussianBandwidth_PD%s.npy' % PD, [logEps, logSumWij])
        plt.clf()

    else: #manually input different epsilon (trial and error) if manifolds generated with above methods not converged
        if SSn == 1:    
            eps = .5
        elif SSn == 2:
            eps = 1e4
    
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
    np.save(os.path.join(outDir, 'PD_%s_SS%s_val.npy' % (PD, SSn)), sdiag)
    np.save(os.path.join(outDir, 'PD_%s_SS%s_vec.npy' % (PD, SSn)), U)
        
    # =========================================================================
    # Analysis of diffusion map
    # =========================================================================
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
        plt.ylim(-sdiag[1]/8., sdiag[1]+sdiag[1]/8.)
        plt.locator_params(nbins=15)
        plt.axhline(y=0, color='k', alpha=.5, linestyle='--', linewidth=1)
        plt.show()
    
    
    # =============================================================================
    # All plots used in Figure 11 and Figure 14
    # =============================================================================
    if SSn == 2: #ordered array 2D manifold subspaces
        dimRows = 4
        dimCols = 8
                
        widths = []
        for i in range(dimCols):
            widths.append(1)
        widths.append(.25)
        
        fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
        gs = GridSpec(dimRows, dimCols+1, left=0.05, right=0.95, wspace=0.3, hspace=0.3, width_ratios=widths)
        
        idx_row = 0
        for v1 in range(1,dimRows+1):
            idx_col = 0
            for v2 in range(v1+1, v1+dimCols+1):
                ax = fig.add_subplot(gs[idx_row, idx_col])

                color=iter(cm.tab20(np.linspace(0, 1, np.shape(CM_idx)[0])))
                for b in range(np.shape(CM_idx)[0]):
                    c=next(color)
                    ax.scatter(U[:,v1][CM_idx[b]], U[:,v2][CM_idx[b]], color=c, s=16, edgecolor='k', linewidths=.1, zorder=1)

                ax.set_xlabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
                ax.set_ylabel(r'$\Psi_{%s}$' % v2, fontsize=18, labelpad=2.5)

                ax.axes.xaxis.set_ticklabels([])
                ax.axes.yaxis.set_ticklabels([])
                ax.set_xticks([])
                ax.set_yticks([])
 
                ax.set_xlim(np.amin(U[:,v1])*1.1, np.amax(U[:,v1])*1.1)
                ax.set_ylim(np.amin(U[:,v2])*1.1, np.amax(U[:,v2])*1.1)
                
                idx_col += 1
            idx_row += 1
            
        ax = fig.add_subplot(gs[:, dimCols])
        
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cm.tab20).reversed())
        sm.set_clim(vmin=0, vmax=U.shape[0])
        cbar = fig.colorbar(sm, cax=ax)
        cbar.ax.tick_params(labelsize=16)

        fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig11_14_PD_%s.pdf' % PD), dpi=600)
        plt.clf()
        

    # =============================================================================
    # All plots used in Figure 9
    # =============================================================================
    if SSn == 2: #ordered array 2D manifold subspaces
        dimRows = 3
        dimCols = 8
                
        widths = []
        for i in range(dimCols):
            widths.append(1)
        widths.append(.25)
        
        fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
        gs = GridSpec(dimRows, dimCols+1, left=0.05, right=0.95, wspace=0.3, hspace=0.3, width_ratios=widths)
        
        idx_row1 = 0
        idx_col1 = 0
        idx_row2 = 0
        idx_col2 = int(dimCols/2)
        for v1 in range(1,(dimRows*int(dimCols/2)+1)):

            # CM1 side:
            ax = fig.add_subplot(gs[idx_row1, idx_col1])
            ax.scatter(list(CM1_idx.flatten()), U[:,v1], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
            ax.set_xlabel(r'$\mathrm{CM}_{1}$', fontsize=14, labelpad=2.5)
            ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

            # CM2 side:
            ax = fig.add_subplot(gs[idx_row2, idx_col2])
            ax.scatter(list(CM2_idx.flatten()), U[:,v1], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
            ax.set_xlabel(r'$\mathrm{CM}_{2}$', fontsize=14, labelpad=2.5)
            ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            
            if idx_col1 == 3:#(dimCols-1):
                idx_row1 += 1
                idx_col1 = 0
                idx_row2 += 1
                idx_col2 = int(dimCols/2)
            else:
                idx_col1 += 1
                idx_col2 += 1

        ax = fig.add_subplot(gs[:, dimCols])
        
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cm.tab20).reversed())
        sm.set_clim(vmin=0, vmax=U.shape[0])
        cbar = fig.colorbar(sm, cax=ax)
        cbar.ax.tick_params(labelsize=16)

        fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig9_PD_%s.pdf' % PD), dpi=600)
        plt.clf()
        
        
    # =============================================================================
    # All plots used in Figure 11 and Figure 14
    # =============================================================================
    if SSn == 1: #ordered array 2D manifold subspaces
        dimRows = 3
        dimCols = 8
        
        # hardcoded coordinates for figure subplots
        V1 = [1,1,1,1,2,2,2,2,3,3,3,3]
        V2 = [2,3,4,5,3,4,5,6,4,5,6,7]
        V_idx = 0
                
        widths = []
        for i in range(dimCols):
            widths.append(1)
        widths.append(.25)
        
        fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
        gs = GridSpec(dimRows, dimCols+1, left=0.05, right=0.95, wspace=0.3, hspace=0.3, width_ratios=widths)
        
        idx_row1 = 0
        idx_col1 = 0
        idx_row2 = 0
        idx_col2 = int(dimCols/2)
        for v1 in range(1,(dimRows*int(dimCols/2)+1)):

            # Indexed side:
            ax = fig.add_subplot(gs[idx_row1, idx_col1])
            ax.scatter(CM_idx, U[:,v1], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
            ax.plot(CM_idx, U[:,v1], c='k', linewidth=1, zorder=-5, alpha=.3)

            ax.set_xlabel(r'$\mathrm{CM}_{1}$', fontsize=14, labelpad=2.5)
            ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=18, labelpad=2.5)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

            # Cartesian-product side:
            ax = fig.add_subplot(gs[idx_row2, idx_col2])
            ax.scatter(U[:,V1[V_idx]], U[:,V2[V_idx]], c=np.arange(U.shape[0]), cmap='tab20', s=16, linewidths=.1, edgecolor='k', zorder=1)
            ax.plot(U[:,V1[V_idx]], U[:,V2[V_idx]], c='k', linewidth=1, zorder=-5, alpha=.3)
            ax.set_xlabel(r'$\Psi_{%s}$' % V1[V_idx], fontsize=18, labelpad=2.5)
            ax.set_ylabel(r'$\Psi_{%s}$' % V2[V_idx], fontsize=18, labelpad=2.5)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            
            if idx_col1 == 3:#(dimCols-1):
                idx_row1 += 1
                idx_col1 = 0
                idx_row2 += 1
                idx_col2 = int(dimCols/2)
            else:
                idx_col1 += 1
                idx_col2 += 1
            V_idx += 1

        ax = fig.add_subplot(gs[:, dimCols])
        
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cm.tab20).reversed())
        sm.set_clim(vmin=0, vmax=U.shape[0])
        cbar = fig.colorbar(sm, cax=ax)
        cbar.ax.tick_params(labelsize=16)

        fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig10_PD_%s.pdf' % PD), dpi=600)
        plt.clf()
        
        
if __name__ == '__main__':
    path1 = os.path.splitext(sys.argv[0])[0]
    path2, tail = os.path.split(path1)
    op(path2, sys.argv[1])