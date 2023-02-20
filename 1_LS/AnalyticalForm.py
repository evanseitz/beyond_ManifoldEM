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
if 1: #FIGURE IN MS: comparing analytical to Figure 3 and Figure 4
    x = np.linspace(0,1,20)
    xx, yy = np.meshgrid(x, x)

    th1 = (-177)*np.pi/180
    k,j = 1, 1
    z1 = (np.cos(th1))*np.cos(k*np.pi*(xx)) + (np.sin(th1))*np.cos(j*np.pi*yy)
    
    th2 = (-100)*np.pi/180
    k,j = 2, 1
    z2 = (np.cos(th2))*np.cos(k*np.pi*(xx)) + (np.sin(th2))*np.cos(j*np.pi*yy)
    
    th3 = (168)*np.pi/180
    k,j = 2, 1
    z3 = (np.cos(th3))*np.cos(k*np.pi*(xx)) + (np.sin(th3))*np.cos(j*np.pi*yy)
    
    th4 = (-182)*np.pi/180
    k,j = 3, 1
    z4 = (np.cos(th4))*np.cos(k*np.pi*(xx)) + (np.sin(th4))*np.cos(j*np.pi*yy)
    
    th5 = (-86)*np.pi/180
    k,j = 4, 2
    z5 = (np.cos(th5))*np.cos(k*np.pi*(xx)) + (np.sin(th5))*np.cos(j*np.pi*yy)
    
    th6 = (-177)*np.pi/180 #183
    k,j = 4, 2
    z6 = (np.cos(th6))*np.cos(k*np.pi*(xx)) + (np.sin(th6))*np.cos(j*np.pi*yy)
        
    fig_idx1 = [1,2,3,4,5,6]
    for f_idx in fig_idx1:
        psi_label = str(f_idx)
        if f_idx == 1:
            Z = z1
            rev = True
        elif f_idx == 2:
            Z = z2
            rev = True
        elif f_idx == 3:
            Z = z3
            rev = True
        elif f_idx == 4:
            Z = z4
            rev = False
        elif f_idx == 5:
            Z = z5
            rev = False
        elif f_idx == 6:
            Z = z6
            rev = False#True
        
        if 1: #X-view
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            if rev is True:
                ax.scatter(xx, yy[::-1], Z, c='white', s=20, linewidths=.5, edgecolor='k')
            else:
                ax.scatter(xx, yy, Z, c='white', s=20, linewidths=.5, edgecolor='k')
            #ax.set_xlabel('$X$', labelpad=lpadx, fontsize=18)
            #ax.set_ylabel('$Y$', labelpad=lpady, fontsize=18)  
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=0, rotation=0, fontsize=18)
            if 1:
                frame = plt.gca()
                frame.axes.xaxis.set_ticklabels([])
                frame.axes.yaxis.set_ticklabels([])
                frame.axes.zaxis.set_ticklabels([])
            else:
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
            ax.view_init(elev=0, azim=-90.1)
            #ax.view_init(elev=0, azim=89.9)
            ax.auto_scale_xyz
            #plt.show()
            fig.savefig('_figure_assets_/Fig12_psi%s_X.pdf' % f_idx, dpi=600)
            plt.clf()
            
        if 1: #Y-view
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.scatter(xx, yy, Z, c='white', s=20, linewidths=.5, edgecolor='k')
            #ax.set_xlabel('$X$', labelpad=lpadx, fontsize=18)
            #ax.set_ylabel('$Y$', labelpad=lpady, fontsize=18)        
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel('$\Psi_{%s}$' % psi_label, labelpad=-5, rotation=0, fontsize=18)    
            if 1:
                frame = plt.gca()
                frame.axes.xaxis.set_ticklabels([])
                frame.axes.yaxis.set_ticklabels([])
                frame.axes.zaxis.set_ticklabels([])
            else:
                ax.xaxis.set_tick_params(labelsize=10)
                ax.yaxis.set_tick_params(labelsize=10)
                ax.zaxis.set_tick_params(labelsize=10)
                ax.tick_params(axis='x', which='major', pad=-1)
                ax.tick_params(axis='y', which='major', pad=-1)
                ax.tick_params(axis='z', which='major', pad=3)
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
            #plt.show()
            fig.savefig('_figure_assets_/Fig12_psi%s_Y.pdf' % f_idx, dpi=600)
            plt.clf()
            
    fig_idx2a = [1,1,2,2]
    fig_idx2b = [2,3,3,5]
    for f_idx in range(4):
        if f_idx == 0:
            Z1, Z2 = z1, z2
        elif f_idx == 1:
            Z1, Z2 = z1, z3
        elif f_idx == 2:
            Z1, Z2 = z2, z3
        elif f_idx == 3:
            Z1, Z2 = z2, z5
        
        if 1: #composite
            s = 30
            lw = .5
            plt.scatter(Z1, Z2, s=s, c='white', edgecolor='k', linewidths=lw, zorder=1)
            plt.xlabel(r'$\Psi_{%s}$' % (fig_idx2a[f_idx]), fontsize=20, labelpad=3)
            plt.ylabel(r'$\Psi_{%s}$' % (fig_idx2b[f_idx]), fontsize=20, labelpad=3)
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
            #plt.show()
            fig.savefig('_figure_assets_/Fig12_psi%s_psi%s.pdf' % (fig_idx2a[f_idx], fig_idx2b[f_idx]), dpi=600)
            plt.clf()
        
        
