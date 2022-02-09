#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################
############################################################

# Author: Viktor Belay 
# Creation date: 02/06/2022

# Purpose of script: Find positions of water molecules in the human TMEM175 pore + plot over time

# Version: 1.0.0

# Funtionalty:
    
    # Inputs:
        # Required: path to topology
        # Required: path to trajectory
        # Optional: plot (True or False). Will plot results with no input
    # Outputs:
        
# Notes: Currently configured only for TMEM175
# Future plans: Configure script to work with any channel

############################################################
############################################################

from MDAnalysis import *
from matplotlib import pyplot as pl
import numpy as np

def get_water_positions(topology,traj,plot=True):
    
    ############################################################
    # Define MDA universe
    
    u = Universe(topology,traj) # Define MDA universe
    
    ############################################################
    
    
    ############################################################
    # Designate waters 
    
    w = u.select_atoms('resname TIP3',updating=True)
    
    ############################################################


    ############################################################
    # Calculate water positions over the length of the trajectory
    
    w_pos = []
    
    for ts in u.trajectory:
        pos.append(w.positions)
        
    w_pos = np.array(w_pos)
    
    
    ############################################################
    
    
    ############################################################
    # Create an array of only z water positions
    
    z_w_pos = []
    
    for i in list(range(0,len(w_pos))):
        
        zz_w_pos = []
        z_w_pos.append(zz_w_pos)
        
        for j in list(range(0,len(w_pos[0]))):
            
            zz_w_pos.append(w_pos[i][j][2])
            
    z_w_pos = np.array(z_w_pos)
    
    
    ############################################################
    
    
    ############################################################
    # Plot stuff
    
    if plot == True:
    
        for i in list(range(0,500)):
            #print(np.array([i]*len(zzzpos[i][(zzzpos[i]>30)&(zzzpos[i]<90)])))
            pl.plot(np.array([i]*len(zzzpos_closed[i][(zzzpos_closed[i]>30)&(zzzpos_closed[i]<90)])),zzzpos_closed[i][(zzzpos_closed[i]>30)&(zzzpos_closed[i]<90)],'b.')
            pl.ylim(52,58)
            pl.xlabel('time (ns)')
            pl.ylabel('water position (~55 = isoleucine constriction)')
            
    else:
        
        pass
    
    ############################################################
    
if __name__ == '__main__':
        
    get_water_positions(topology,traj,plot=True)