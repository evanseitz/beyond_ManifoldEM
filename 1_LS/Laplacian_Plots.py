import sys, os, re
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
from matplotlib import rc
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from pylab import imshow, show, axes
from mpl_toolkits.mplot3d import Axes3D

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

x = np.linspace(0,1,20)
xx, yy = np.meshgrid(x, x)


# =============================================================================
# Bottom of Figure 4
# =============================================================================
if 0:
    a, b = 3,1
    z = np.cos(a*np.pi*xx)*np.cos(b*np.pi*yy) #test for z1, z2

    if 1: #ROTATIONS PLOT FOR MS FIGURE
        Ord = [1,8,2,9,3,10,4,11,5,12,6,13,7,14]
        idx = 0
        for a in np.arange(0,90+1,15):
            ang = -a*np.pi/180
            z1 = (np.cos(ang))*(np.cos(1*np.pi*xx)) + (np.sin(ang))*(np.cos(1*np.pi*yy))
            z2 = (-np.sin(ang))*(np.cos(1*np.pi*xx)) + (np.cos(ang))*(np.cos(1*np.pi*yy))
            plt.subplot(2,7,Ord[idx])
            plt.title(r'$R_{1,2}(%s^\circ)$' % a, fontsize=16)
            idx += 1
            plt.scatter(xx, z1, c='white', s=20, linewidths=.5, edgecolor='k', zorder=0)

            plt.gca().set_box_aspect(1)
            plt.xlabel(r'$X$', fontsize=12)
            plt.ylabel(r'$\Psi_1$', fontsize=14)
            plt.subplot(2,7,Ord[idx])
            idx += 1
            plt.scatter(xx, z2, c='white', s=20, linewidths=.5, edgecolor='k', zorder=0)
            plt.gca().set_box_aspect(1)
            plt.xlabel(r'$X$', fontsize=12)
            plt.ylabel(r'$\Psi_2$', fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(left=0.094, right=0.913, bottom=0.031, top=0.426, wspace=0.488, hspace=0.0)

        plt.show()
        plt.clf()
        
        
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
    
    
    '''if 1:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(z1, z2, z1.T)
        ax.set_xlabel('$X$')
        ax.set_ylabel('$Y$')
        ax.set_zlabel('$Z$')
        ax.auto_scale_xyz
        plt.show()
        
    if 1:
        plt.scatter(z1,z2, s=1)
        plt.xlabel('$Z_{1}$')
        plt.xlabel('$Z_{2}$')
        plt.show()
        
    if 0:
        plt.scatter(xx,z)
        #plt.scatter(yy,z)
        plt.show()
        plt.scatter(z,z.T)
        plt.show()'''
        
        
    
if 0:
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
    

# =============================================================================
# Bottom of Figure 5
# =============================================================================
if 0: #FIGURE FOR MS (REFERENCE FRAMES)
    psi = -np.cos(1*np.pi*xx)*np.cos(2*np.pi*yy) - np.cos(2*np.pi*xx)*np.cos(np.pi*yy) #PSI 6
    #ang1 = (250)*np.pi/180
    #psi = (1*np.cos(ang1))*(np.cos(1*np.pi*xx)) + (np.sin(ang1))*(np.cos(1*np.pi*yy)) #PSI 1
    ang2 = (250)*np.pi/180
    psi = (-np.sin(ang2))*(np.cos(1*np.pi*xx)) + (np.cos(ang2))*(np.cos(1*np.pi*yy)) #PSI 2
    psi_label='2'

    if 1: #external view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
        ax.set_xlabel('$X$', labelpad=-2)
        ax.set_ylabel('$Y$', labelpad=-2)
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=-2, rotation=0)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.zaxis.set_tick_params(labelsize=6)
        ax.tick_params(axis='x', which='major', pad=-1)
        ax.tick_params(axis='y', which='major', pad=-1)
        ax.tick_params(axis='z', which='major', pad=0)
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
        ax.view_init(elev=30, azim=-45)
        ax.auto_scale_xyz
        plt.show()
        
    if 1: # X-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
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
        ax.auto_scale_xyz
        plt.show()
        
    if 1: # Y-view
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(xx, yy, psi, c='white', s=20, linewidths=.5, edgecolor='k')
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
        ax.auto_scale_xyz
        plt.show()
