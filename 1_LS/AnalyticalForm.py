import sys, os, re
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
from matplotlib import rc
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from pylab import imshow, show, axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams



plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
})

'''states_gt = 20 #ground-truth number of states along each degree of freedom (assumed symmetric)
tau = 1
ss_gt = 400
CM1_idx = np.ndarray(shape=(states_gt, states_gt*tau), dtype=int)  
CM2_idx = np.ndarray(shape=(states_gt, states_gt*tau), dtype=int)  

# ========================================================================
# Generate CM1 ground-truth indexing:
# ========================================================================
idx=0
shift=0
Idx=0
for i in range(ss_gt):
    if Idx*tau <= i < (Idx+states_gt)*tau:
        CM1_idx[shift, idx-Idx*tau] = i
        idx+=1
        if idx%(states_gt*tau) == 0:
            Idx+=states_gt
            shift+=1
            
# ========================================================================
# Generate CM2 ground-truth indexing:
# ======================================================================== 
binsActual = []
Idx = 0
for s in range(0,states_gt):
    state_list = []
    for r in range(0,states_gt):
        state_list.append(np.arange(Idx,Idx+tau))
        Idx+=(states_gt*tau)
    binsActual.append([item for sublist in state_list for item in sublist])
    Idx-=(ss_gt-tau)

for i in range(states_gt):
    for j in range(states_gt*tau):
        CM2_idx[i,j] = binsActual[i][j]'''
        
        
# ========================================================================
# Generate analytical manifold:
# ========================================================================



#if 1:
    #a, b = 3,1
    #z = np.cos(a*np.pi*xx)*np.cos(b*np.pi*yy)
    
    
# =============================================================================
# Bottom of Figure 4
# =============================================================================
if 1:
    x = np.linspace(0,1,20)
    xx, yy = np.meshgrid(x, x)

    dimCols = 7
    dimRows = 2
    fig = plt.figure(figsize=[dimCols*2,dimRows*2], constrained_layout=True)
    gs = GridSpec(dimRows, dimCols, left=0.05, right=0.85, wspace=0.1, hspace=0.1)
    
    
    #Ord = [1,8,2,9,3,10,4,11,5,12,6,13,7,14]
    idx_col1 = 0
    for a in np.arange(0,90+1,15):
        ang = -a*np.pi/180
        z1 = (np.cos(ang))*(np.cos(1*np.pi*xx)) + (np.sin(ang))*(np.cos(1*np.pi*yy))
        z2 = (-np.sin(ang))*(np.cos(1*np.pi*xx)) + (np.cos(ang))*(np.cos(1*np.pi*yy))
        
        rcParams['axes.titlepad'] = 10
        ax = fig.add_subplot(gs[0, idx_col1])
        ax.set_title(r'$R_{1,2}(%s^\circ)$' % a, fontsize=18)
        ax.scatter(xx, z1, c='white', s=20, linewidths=.5, edgecolor='k', zorder=0)
        ax.set_box_aspect(1)
        ax.set_xlim(-.354,1.354)
        ax.set_ylim(-1.55,1.55)
        ax.axes.xaxis.set_ticklabels([])
        if idx_col1 != 0:
            ax.axes.yaxis.set_ticklabels([])
        else:
            ax.set_ylabel(r'$\Psi_1$', fontsize=18, labelpad=3)
        
        ax = fig.add_subplot(gs[1, idx_col1])            
        ax.scatter(xx, z2, c='white', s=20, linewidths=.5, edgecolor='k', zorder=0)
        ax.set_box_aspect(1)
        ax.set_xlim(-.354,1.354)
        ax.set_ylim(-1.55,1.55)
        if idx_col1 != 0:
            ax.axes.yaxis.set_ticklabels([])
        else:
            ax.set_ylabel(r'$\Psi_2$', fontsize=18, labelpad=3)
        ax.set_xlabel(r'$X$', fontsize=18, labelpad=3)
        
        idx_col1 += 1
        
    plt.tight_layout()
    fig.savefig('_figure_assets_/Fig4_bot.pdf', dpi=600)
    #plt.show()
    plt.clf()
    plt.close()
    

# =============================================================================
# Bottom of Figure 5
# =============================================================================
x = np.linspace(0,1,50)
xx, yy = np.meshgrid(x, x)
if 1:
    angle = 45 #angle used in DiffusionMaps_LS.py to counter rotate manifold
    if 1:
        ang1 = (angle)*np.pi/180
        psi = (1*np.cos(ang1))*(np.cos(1*np.pi*xx)) + (np.sin(ang1))*(np.cos(1*np.pi*yy)) #PSI 1
        psi_label = '1'
    elif 0:
         ang2 = (angle)*np.pi/180
         psi = (-np.sin(ang2))*(np.cos(1*np.pi*xx)) + (np.cos(ang2))*(np.cos(1*np.pi*yy)) #PSI 2
         psi_label = '2'
    else: #used for Table of Contents
        psi = -np.cos(1*np.pi*xx)*np.cos(2*np.pi*yy) - np.cos(2*np.pi*xx)*np.cos(np.pi*yy) #PSI 6
        psi_label = '6'

    if 1: #external view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=4, fontsize=18)
        ax.set_ylabel('$Y$', labelpad=4, fontsize=18)
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=4, rotation=0, fontsize=18)
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        ax.zaxis.set_tick_params(labelsize=10)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=1)
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('k')
        ax.yaxis.pane.set_edgecolor('k')
        ax.zaxis.pane.set_edgecolor('k')
        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.xaxis._axinfo["grid"]['linewidth'] = .1
        ax.yaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.yaxis._axinfo["grid"]['linewidth'] = .1
        ax.zaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.zaxis._axinfo["grid"]['linewidth'] = .1
        #ax.grid(False)
        #ax.view_init(elev=30, azim=45) #external view for Psi6
        ax.view_init(elev=35, azim=-45)
        ax.auto_scale_xyz
        #plt.show()
        fig.savefig('_figure_assets_/Fig5B1_psi%s.pdf' % psi_label, dpi=600)
        plt.clf()
        plt.close()
        
        
    if 1: # X-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=8, fontsize=18)
        ax.set_ylabel('$Y$', labelpad=-3, fontsize=18)
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=23, rotation=0, fontsize=18)
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        ax.zaxis.set_tick_params(labelsize=10)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=8)
        ax.set_yticks([])
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('k')
        ax.yaxis.pane.set_edgecolor('k')
        ax.zaxis.pane.set_edgecolor('k')
        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.xaxis._axinfo["grid"]['linewidth'] = .1
        ax.yaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.yaxis._axinfo["grid"]['linewidth'] = .1
        ax.zaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.zaxis._axinfo["grid"]['linewidth'] = .1
        #ax.grid(False)
        if psi_label == '6':
            ax.view_init(elev=0, azim=-90.1)
            for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
                axis.set_ticklabels([])
                axis._axinfo['axisline']['linewidth'] = 1
                axis._axinfo['axisline']['color'] = (0, 0, 0)
                axis._axinfo['grid']['linewidth'] = 0.5
                axis._axinfo['grid']['linestyle'] = "-"
                axis._axinfo['grid']['color'] = (0, 0, 0)
                axis._axinfo['tick']['inward_factor'] = 0.0
                axis._axinfo['tick']['outward_factor'] = 0.0
                axis.set_pane_color((0.95, 0.95, 0.95))

            
        else:
            ax.view_init(elev=0, azim=-90.1)
        ax.auto_scale_xyz
        #plt.show()
        fig.savefig('_figure_assets_/Fig5B2_psi%s.pdf' % psi_label, dpi=600)
        plt.clf()
        plt.close()
        
    if 1: # Y-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=-3, fontsize=18)
        ax.set_ylabel('$Y$', labelpad=8, fontsize=18)      
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=5, rotation=0, fontsize=18)       
        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        ax.zaxis.set_tick_params(labelsize=10)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=0)
        ax.set_xticks([])
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('k')
        ax.yaxis.pane.set_edgecolor('k')
        ax.zaxis.pane.set_edgecolor('k')
        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.xaxis._axinfo["grid"]['linewidth'] = .1
        ax.yaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.yaxis._axinfo["grid"]['linewidth'] = .1
        ax.zaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.zaxis._axinfo["grid"]['linewidth'] = .1
        #ax.grid(False)
        ax.view_init(elev=0, azim=0)
        ax.auto_scale_xyz
        #plt.show()
        fig.savefig('_figure_assets_/Fig5C_psi%s.pdf' % psi_label, dpi=600)
        plt.clf()
        plt.close()

      
        
      
        
      
        
      
        
      
        
      
        
# =============================================================================
# Figure 12
# =============================================================================
if 0: #FIGURE IN MS: comparing analytical to Figure 3 and Figure 4
    #a1 = .90
    #a2 = .99
    #z1 = (a1)*np.cos(1*np.pi*xx) + (1-a1)*np.cos(2*np.pi*yy)
    
    #Psi1: -181, {1,1}, -X, Y
    #Psi2: -260, {2,1}, X, Y
    #Psi3: 187, {2,1}, X, -Y
    #Psi4: 5, {3,1}, X, -Y
    #Psi5: -92, {4,2}, X, -Y
    #Psi6:  179, {4,2}, X, Y

    th1 = (-181)*np.pi/180
    k,j = 1, 1
    z1 = (np.cos(th1))*np.cos(k*np.pi*(xx)) + (np.sin(th1))*np.cos(j*np.pi*yy)
    
    th2 = (-260)*np.pi/180
    k,j = 2, 1
    z2 = (np.cos(th2))*np.cos(k*np.pi*(xx)) + (np.sin(th2))*np.cos(j*np.pi*yy)
    
    th3 = (187)*np.pi/180
    k,j = 2, 1
    z3 = (np.cos(th3))*np.cos(k*np.pi*(xx)) + (np.sin(th3))*np.cos(j*np.pi*yy)
    
    th5 = (-92)*np.pi/180
    k,j = 4, 2
    z5 = (np.cos(th5))*np.cos(k*np.pi*(xx)) + (np.sin(th5))*np.cos(j*np.pi*yy)
    
    th6 = (179)*np.pi/180
    k,j = 4, 2
    z6 = (np.cos(th6))*np.cos(k*np.pi*(xx)) + (np.sin(th6))*np.cos(j*np.pi*yy)

    psi_label = '6'
    
    if 1: # X-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, z6, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=.5)
        ax.set_ylabel('$Y$', labelpad=-10)
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=5, rotation=0)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.zaxis.set_tick_params(labelsize=6)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=3)
        ax.set_yticks([])
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('k')
        ax.yaxis.pane.set_edgecolor('k')
        ax.zaxis.pane.set_edgecolor('k')
        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.xaxis._axinfo["grid"]['linewidth'] = .1
        ax.yaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.yaxis._axinfo["grid"]['linewidth'] = .1
        ax.zaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.zaxis._axinfo["grid"]['linewidth'] = .1
        #ax.grid(False)
        ax.view_init(elev=0, azim=-90.1)
        #ax.view_init(elev=0, azim=89.9)

        ax.auto_scale_xyz
        plt.show()
        
    if 1: # Y-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, z6, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=-10)
        ax.set_ylabel('$Y$', labelpad=.5)        
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=-2, rotation=0)        
        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.zaxis.set_tick_params(labelsize=6)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=0)
        ax.set_xticks([])
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('k')
        ax.yaxis.pane.set_edgecolor('k')
        ax.zaxis.pane.set_edgecolor('k')
        ax.xaxis.pane.set_alpha(1)
        ax.yaxis.pane.set_alpha(1)
        ax.zaxis.pane.set_alpha(1)
        ax.xaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.xaxis._axinfo["grid"]['linewidth'] = .1
        ax.yaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.yaxis._axinfo["grid"]['linewidth'] = .1
        ax.zaxis._axinfo["grid"].update({"linewidth":1, "color" : "black"})
        ax.zaxis._axinfo["grid"]['linewidth'] = .1
        #ax.grid(False)
        ax.view_init(elev=0, azim=0)
        #ax.view_init(elev=0, azim=180.1)
        ax.auto_scale_xyz
        plt.show()
        
    if 1: # composite
        fS = 16
        s = 30
        lw = .5
        #plt.scatter(z1, z2, s=s, c='white', edgecolor='k', linewidths=lw, zorder=1) #PSI 1, PSI 2
        #plt.scatter(z1[::-1], z3[::-1], s=s, c='white', edgecolor='k', linewidths=lw, zorder=1) #PSI 1, PSI 3
        #plt.scatter(z2, z3, s=s, c='white', edgecolor='k', linewidths=lw, zorder=1) #PSI 2, PSI 3

        
        plt.scatter(z2, z5, s=s, c='white', edgecolor='k', linewidths=lw, zorder=1)

        plt.xlabel(r'$\Psi_{%s}$' % (2), fontsize=fS+2, labelpad=2.5)
        plt.ylabel(r'$\Psi_{%s}$' % (5), fontsize=fS+2, labelpad=2.5)
        

        
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
    
        #fig = plt.figure()
        #ax = plt.axes()
        #ax.scatter(z1, z2, c='white', s=20, linewidths=.5, edgecolor='k')
        #ax.set_xlabel('$X$', labelpad=.5)
        #ax.set_ylabel('$Y$', labelpad=-10)
        #ax.zaxis.set_rotate_label(False)
        #ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=5, rotation=0)
        #ax.xaxis.set_tick_params(labelsize=6)
        #ax.yaxis.set_tick_params(labelsize=6)
        #ax.zaxis.set_tick_params(labelsize=6)
        #ax.tick_params(axis='x', which='major', pad=-1)
        #ax.tick_params(axis='y', which='major', pad=-1)
        #ax.tick_params(axis='z', which='major', pad=3)
        #ax.set_yticks([])

        #ax.auto_scale_xy
        plt.show()
    
    
    if 0:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(z1, z2, z1.T)
        ax.set_xlabel('$X$')
        ax.set_ylabel('$Y$')
        ax.set_zlabel('$Z$')
        ax.auto_scale_xyz
        plt.show()
        
    if 0:
        plt.scatter(z1,z2, s=1)
        plt.xlabel('$Z_{1}$')
        plt.xlabel('$Z_{2}$')
        plt.show()
        
    if 0:
        plt.scatter(xx,z)
        #plt.scatter(yy,z)
        plt.show()
        plt.scatter(z,z.T)
        plt.show()
        
        
    
if 0: #UNKNOWN
    A, B = 1, 1
    z1 = A*np.cos(1*np.pi*xx/1) + B*np.cos(1*np.pi*yy/1)*-2
    #z2 = A*np.cos(2*np.pi*xx) + B*np.cos(1*np.pi*yy)
    z2 = A*np.cos(2*np.pi*xx/1) + B*np.cos(2*np.pi*yy/1)
    
    if 1:
        plt.suptitle('$Z_1=%s\mathrm{cos}(1\pi x) + %s\mathrm{cos}(1\pi y)$, $\ Z_2=%s\mathrm{cos}(2\pi x) + %s\mathrm{cos}(2\pi y)$' % (A,B,A,B), y=.81)
        plt.subplot(1,5,1)
        plt.scatter(xx,z1, c='C0', s=15, edgecolor='k', linewidths=.1)
        plt.xlabel('$X$')
        plt.ylabel('$Z_1$')
        plt.gca().set_box_aspect(1)
        
        plt.subplot(1,5,2)
        plt.scatter(xx,z2, c='C0', s=15, edgecolor='k', linewidths=.1)
        plt.xlabel('$X$')
        plt.ylabel('$Z_2$')
        plt.gca().set_box_aspect(1)
    
        plt.subplot(1,5,3)
        plt.scatter(yy,z1, c='C1', s=15, edgecolor='k', linewidths=.1)
        plt.xlabel('$Y$')
        plt.ylabel('$Z_1$')
        plt.gca().set_box_aspect(1)
    
        plt.subplot(1,5,4)
        plt.scatter(yy,z2, c='C1', s=15, edgecolor='k', linewidths=.1)
        plt.xlabel('$Y$')
        plt.ylabel('$Z_2$')
        plt.gca().set_box_aspect(1)
    
        if 1:
            plt.subplot(1,5,5)
            plt.scatter(z1,z2, c='C2', s=15, edgecolor='k', linewidths=.1) #parabola
            plt.xlabel('$Z_1$')
            plt.ylabel('$Z_2$')
            plt.gca().set_box_aspect(1)
        else:
            plt.subplot(1,5,5)
            plt.scatter(z1.T,z2.T, c='C2', s=15, edgecolor='k', linewidths=.1) #parabola
            plt.xlabel('$Z_1^T$')
            plt.ylabel('$Z_2^T$')
            plt.gca().set_box_aspect(1)
        
        plt.tight_layout()
        plt.subplots_adjust(left=.05, bottom=0.3, right=.95, top=1, wspace=.25, hspace=.2)
        
    if 1:
        plt.scatter(z1,z2, c='C2', s=15, edgecolor='k', linewidths=.1) #parabola
        plt.xlabel('$Z_1$')
        plt.ylabel('$Z_2$')
        plt.gca().set_box_aspect(1)

    plt.show()
    


