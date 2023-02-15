import sys, os, re
sys.dont_write_bytecode = True
import numpy as np
from numpy import linalg as LA
import matplotlib
from matplotlib import rc
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from pylab import imshow, show, loadtxt, axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import cm
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import scipy
from scipy import linalg
import Generate_Nd_Rot


pyDir = os.path.dirname(os.path.abspath(__file__)) #python file location

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

# =============================================================================
# User parameters for rendering plots for Figures 2, 4 and 8
# =============================================================================

SSn = 2 #options : {1 or 2} for Figure 2 (top) or (bottom), respectively

Nx=50
Ny=50

if 1: #Figure 2
    x = np.linspace(0,1.1,Nx) #nondegenerate
    bc = 'rect'
    rotated = False
else: #Figure 4 (top)
    x = np.linspace(0,1,Nx)
    bc = 'square'
    if 0: #counter rotate degenerate eigenspace (Figure 4 LHS/RHS before and after, respectively)
        rotated = True
    else: #no counter rotation performed
        rotated = False

y = np.linspace(0,1,Ny)
xx, yy = np.meshgrid(x, y)
pts = np.vstack([xx.ravel(), yy.ravel()])

if 1: #Figure 2 (LHS) or Figure 4 (top) depending on above choice
    BigEps = False
    figname = 'smallEps'
    eps = .00005 #or use 0.0003 for {Nx=20, Ny=20}
else: #Figure 2 (RHS)
    BigEps = True
    figname = 'largeEps'
    eps = 10


# ========================================================================
# Generate CM2 ground-truth indexing:
# ======================================================================== 
if 0:
    binsActual = []
    Idx = 0
    tau = 1
    Y_idx = np.ndarray(shape=(Nx, Nx*tau), dtype=int)  
    for s in range(0,Nx):
        state_list = []
        for r in range(0,Ny):
            state_list.append(np.arange(Idx,Idx+tau))
            Idx+=(Nx*tau)
        binsActual.append([item for sublist in state_list for item in sublist])
        Idx-=(Nx*Ny-tau)
    
    for i in range(Nx):
        for j in range(Ny*tau):
            Y_idx[i,j] = binsActual[i][j]    
              

if 0: #plot latent space
    plt.scatter(pts[0,:],pts[1,:], s=1)
    plt.axis('equal')
    plt.show()
    plt.clf()
        
    
if 0: #affine transformation
    # ========================================================================
    # general affine transform in 2D:
    #   xs = a1*x + a2*y + b1
    #   ys = a3*x + a4*y + b2
    # ========================================================================
    if 1:
        ptsCP = np.copy(pts)
        for i in range(Nx*Ny):
            pts[0,i] = ptsCP[0,i] + .25*ptsCP[1,i]
    elif 0:
        A = np.array([[1,1.0001],[0,1]]) #shear
        #A = np.array([[.5,0],[0,1]]) #expansion/contraction
        pts = np.matmul(A, pts)
    else:
        print('...')
    if 1:
        plt.scatter(pts[0,:],pts[1,:], s=5)
        plt.axis('equal')
        plt.show()
    
if 1:
    if SSn == 1:
        m = Nx #1 doF
    elif SSn == 2:
        m = Nx*Ny #2 doF
    
    # Compute distances:
    Dist = np.zeros((m,m)) #distance matrix
    p = 2. #Minkowski distance metric: p1=Manhattan, p2=Euclidean, etc.
    for i in range(0,m):
        for j in range(0,m):
            Dist[i,j] = np.sqrt((pts[0,i]-pts[0,j])**2 + (pts[1,i]-pts[1,j])**2)
    
    if 0: #plot distance matrix
        imshow(Dist, cmap='jet', origin='lower', interpolation='nearest')
        plt.title('Distances', fontsize=20)
        plt.colorbar()
        plt.tight_layout()
        plt.show()

    # =============================================================================
    # LS plot used in Figure 8
    # =============================================================================
    if Nx == 20 and Ny == 20: #plot distances of state_01_01 to all others
        """
            Note: specifically for this plot, must set Nx=20 and Ny=20 at top of this scipt
            Change back to Nx=50 and Ny=50 for all other plots
        """
        fig, ax = plt.subplots()
        ax.scatter(np.linspace(1,m,m), Dist[0,:], s=20, c='white', edgecolor='k', linewidths=.5, zorder=1)
        point1 = [1, Dist[0,:][0]]
        point2 = [20, Dist[0,:][19]]
        x_values1 = [point1[0], point2[0]]
        y_values1 = [point1[1], point2[1]]
        ax.plot(x_values1, y_values1, c='red', linewidth=1.5, zorder=-1, label='$X$')
        point3 = [1, Dist[0,:][0]]
        point4 = [381, Dist[0,:][380]]
        x_values2 = [point3[0], point4[0]]
        y_values2 = [point3[1], point4[1]]
        ax.plot(x_values2, y_values2, c='blue', linewidth=1.5, zorder=-1, label='$Y$')
        plt.legend(loc='lower right', fontsize=14, frameon=False)
        plt.xlabel(r'States', fontsize=20, labelpad=12)
        plt.ylabel(r'$D_{1,k}$', fontsize=22, labelpad=10)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        #plt.show()
        fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig8_LS_dists.pdf'), dpi=600)
        plt.clf()
            
        
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
    Dhalf = linalg.sqrtm(D)
    Dinv_half = linalg.sqrtm(Dinv)
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
    dim = 20 #analyze 'dim'-dimensional manifold subspace with 'dim'('dim'-1)/2 rotation operators
    if rotated is True:
        theta_total = int(dim*(dim-1)/2) #total number of rotation operators
        U_init = U[:,1:dim+1]
        thetas = np.zeros(shape=(theta_total,1), dtype=float)
        thetas[0] = 45*np.pi/180#-19*np.pi/180 #Rij operator
        #thetas[theta_total-1] = 45*np.pi/180 #Rij operator
        thetas[((dim-1)+(dim-2)+(dim-3)+1)-1] = 45*np.pi/180 #Rij operator
        R = Generate_Nd_Rot.genNdRotations(dim, thetas)
        U_rot = np.matmul(R, U_init.T)
        U = U_rot.T

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
    # All plots used in Figure 2 (repeat for large epsilon regime, and vice versa above)
    # =============================================================================
    if 1: #ordered array 2D manifold subspaces
    
        dimRows = 4
        dimCols = 5
    
        labpadX = 3.25
        labpadY = 3
        labsizeX = 18
        labsizeY = 22
    
        # hardcoded coordinates for figure subplots
        V1_top = [1,1,1,1,1]
        V2_top = [2,3,4,5,6]
        V1_bot = [2,2,2,2,2]
        V2_bot = [3,4,5,6,7]
        V_idx = 0
                
        widths = []
        for i in range(dimCols):
            widths.append(1)
        
        fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
        gs = GridSpec(dimRows, dimCols, left=0.05, right=0.95, wspace=0.3, hspace=0.3, width_ratios=widths)
        
        idx_col1 = 0  
        for v1 in range(1,6):#(2*int(dimCols/2)+1)):
    
            # CM1 side:
            ax = fig.add_subplot(gs[0, idx_col1])
            if rotated is False:
                ax.scatter(enum, U[:,v1], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            elif rotated is True:
                ax.scatter(enum, U[:,v1-1], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            if SSn == 1:
                ax.plot(enum, U[:,v1], c='k', linewidth=1, zorder=-5, alpha=.3)
            ax.set_xlabel(r'$X$', fontsize=labsizeX, labelpad=labpadX)
            ax.set_ylabel(r'$\Psi_{%s}$' % v1, fontsize=labsizeY, labelpad=labpadY)            
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            
            if BigEps is True: #don't hide Legendre-like form
                ylims = ax.get_ylim()
                if np.abs(ylims[0]) > np.abs(ylims[1]):
                    ax.set_ylim(-1.*np.abs(ylims[0]), np.abs(ylims[0]))
                else:
                    ax.set_ylim(-1.*np.abs(ylims[1]), np.abs(ylims[1]))
    
            # CM2 side:
            ax = fig.add_subplot(gs[1, idx_col1])
            if rotated is False:
                ax.scatter(enum, U[:,v1+dimCols], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            elif rotated is True:
                ax.scatter(enum, U[:,v1+dimCols-1], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            if SSn == 1:
                ax.plot(enum, U[:,v1+dimCols], c='k', linewidth=1, zorder=-5, alpha=.3)
            ax.set_xlabel(r'$X$', fontsize=labsizeX, labelpad=labpadX)
            ax.set_ylabel(r'$\Psi_{%s}$' % (v1+dimCols), fontsize=labsizeY, labelpad=labpadY)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
    
            if BigEps is True: #don't hide Legendre-like form
                ylims = ax.get_ylim()
                if np.abs(ylims[0]) > np.abs(ylims[1]):
                    ax.set_ylim(-1.*np.abs(ylims[0]), np.abs(ylims[0]))
                else:
                    ax.set_ylim(-1.*np.abs(ylims[1]), np.abs(ylims[1]))
            
            # Cartesian-product side:
            ax = fig.add_subplot(gs[2, idx_col1])
            ax.scatter(U[:,V1_top[V_idx]], U[:,V2_top[V_idx]], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            ax.set_xlabel(r'$\Psi_{%s}$' % V1_top[V_idx], fontsize=labsizeY, labelpad=labpadY)
            ax.set_ylabel(r'$\Psi_{%s}$' % V2_top[V_idx], fontsize=labsizeY, labelpad=labpadY)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
    
            ax = fig.add_subplot(gs[3, idx_col1])
            ax.scatter(U[:,V1_bot[V_idx]], U[:,V2_bot[V_idx]], s=s, c='white', linewidths=lw, edgecolor='k', zorder=1)
            ax.set_xlabel(r'$\Psi_{%s}$' % V1_bot[V_idx], fontsize=labsizeY, labelpad=labpadY)
            ax.set_ylabel(r'$\Psi_{%s}$' % V2_bot[V_idx], fontsize=labsizeY, labelpad=labpadY)
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            
            idx_col1 += 1
            V_idx += 1
    
        if bc == 'rect':
            fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig2_%s_SS%s.pdf' % (figname,SSn)), dpi=600)
        elif bc == 'square':
            if rotated is True:
                fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig4_%s_SS%s_rot.pdf' % (figname,SSn)), dpi=600)
            else:
                fig.savefig(os.path.join(pyDir,'_figure_assets_/Fig4_%s_SS%s.pdf' % (figname,SSn)), dpi=600)

        plt.clf()

