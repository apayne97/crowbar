#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################
############################################################

# Author: Viktor Belay 
# Creation date: 02/06/2022

# Purpose of script: Calculate angles of kinks in 
# transmembrane helices which form the human TMEM175 permeation pore.

# Version: 1.0.0

# Funtionalty:
    
    # Inputs:
        # Required: path to topology
        # Required: path to trajectory
        # Optional: protomer (a or b or none)
        # Optional: tm helicers (1 or 7 or none)
    # Outputs:
        
# Notes: Currently configured only for TMEM175
# Future plans: Configure script to work with any protein 

############################################################
############################################################

# Import necessary modules/libraries

from MDAnalysis import *
import numpy as np
from numpy.linalg import norm


def get_pore_angles(topology,traj,protomer='both',tm='both'):
    
    u = Universe(topology,traj) # define MDA universe for top/traj
    
   ############################################################
    
    # Define angle calculation functions
    
    def get_tm1_a_angle(universe):
        
        A = universe.select_atoms('resid 34 and backbone and segid PROA').center_of_geometry()
        B = universe.select_atoms('resid 48 and backbone and segid PROA').center_of_geometry()
        C = universe.select_atoms('resid 56 and backbone and segid PROA').center_of_geometry()

        AB = B - A
        BC = B - C
        angle = np.arccos((np.dot(AB,BC))/(norm(AB)*norm(BC)))
        angle1=np.rad2deg(angle)
        return angle1
        
    def get_tm1_b_angle(universe):
        
        A = universe.select_atoms('resid 34 and backbone and segid PROB').center_of_geometry()
        B = universe.select_atoms('resid 48 and backbone and segid PROB').center_of_geometry()
        C = universe.select_atoms('resid 56 and backbone and segid PROB').center_of_geometry()

        AB = B - A
        BC = B - C
        angle = np.arccos((np.dot(AB,BC))/(norm(AB)*norm(BC)))
        angle1=np.rad2deg(angle)
        return angle1


    def get_tm7_a_angle(universe):
        
        A = universe.select_atoms('resid 260 and backbone and segid PROA').center_of_geometry()
        B = universe.select_atoms('resid 270 and backbone and segid PROA').center_of_geometry()
        C = universe.select_atoms('resid 283 and backbone and segid PROA').center_of_geometry()
        AB = B - A
        BC = B - C    
        angle = np.arccos((np.dot(AB,BC))/(norm(AB)*norm(BC)))
        angle1=np.rad2deg(angle)
        return angle1
    
    def get_tm7_b_angle(universe):
        
        A = universe.select_atoms('resid 260 and backbone and segid PROB').center_of_geometry()
        B = universe.select_atoms('resid 270 and backbone and segid PROB').center_of_geometry()
        C = universe.select_atoms('resid 283 and backbone and segid PROB').center_of_geometry()
        AB = B - A
        BC = B - C    
        angle = np.arccos((np.dot(AB,BC))/(norm(AB)*norm(BC)))
        angle1=np.rad2deg(angle)
        return angle1
    
    
    ############################################################
    # Create main data array
    global data
    data = [] # all data generated below (time, angles) will be appended here
    #data = np.array(data)
    
    
    ############################################################
    
    
    
    ############################################################
    
    ############################################################
    
    # Calculate trajectory time
    
    t = []
    
    for ts in u.trajectory:
        t.append(u.trajectory.time)
        
    real_t = list(range(0,len(t)))
    real_t = np.array(real_t)
    
    data.append(real_t); print('time appended')
    
    ############################################################
    
    ############################################################
    
    # Define set of if statements which compute angles based on input from user
    
    if protomer == 'both':
        
        if tm == 'both':
        
            tm1_a_angle = []
            tm1_b_angle = []
            tm7_a_angle = []
            tm7_b_angle = []
        
            for ts in u.trajectory:
            
                tm1_a_angle.append(get_tm1_a_angle(u))
                tm1_b_angle.append(get_tm1_b_angle(u))
                tm7_a_angle.append(get_tm7_a_angle(u))
                tm7_b_angle.append(get_tm7_b_angle(u))
            
            tm1_a_angle = np.array(tm1_a_angle)
            tm1_b_angle = np.array(tm1_b_angle)
            tm7_a_angle = np.array(tm7_a_angle)
            tm7_b_angle = np.array(tm7_b_angle)
            
            data.append(tm1_a_angle)
            data.append(tm1_b_angle)
            data.append(tm7_a_angle)
            data.append(tm7_b_angle)
            print('angles appended')
            
            data = np.array(data)
            return data

        

        
        elif tm == '1':
            
            tm1_a_angle = []
            tm1_b_angle = []
            
            for ts in u.trajectory:
            
                tm1_a_angle.append(get_tm1_a_angle(u))
                tm1_b_angle.append(get_tm1_b_angle(u))
                
            tm1_a_angle = np.array(tm1_a_angle)
            tm1_b_angle = np.array(tm1_b_angle)
            
            data.append(tm1_a_angle)
            data.append(tm1_b_angle)
            
            return data
            
        elif tm == '7':
            
            tm7_a_angle = []
            tm7_b_angle = []
            
            for ts in u.trajectory:
            
                tm7_a_angle.append(get_tm7_a_angle(u))
                tm7_b_angle.append(get_tm7_b_angle(u))
                
            tm7_a_angle = np.array(tm7_a_angle)
            tm7_b_angle = np.array(tm7_b_angle)
            
            data.append(tm7_a_angle)
            data.append(tm7_b_angle)
            
            return data
            
        else:
            
            print('Error. invalid tm input. leave blank or enter "1" or "7". exiting.')
        
        
    elif protomer==('a' or 'A'):
        
        if tm == 'both':
        
            tm1_a_angle = []

            tm7_a_angle = []
            
            for ts in u.trajectory:
            
                tm1_a_angle.append(get_tm1_a_angle(u))

                tm7_a_angle.append(get_tm7_a_angle(u))

            
            tm1_a_angle = np.array(tm1_a_angle)

            tm7_a_angle = np.array(tm7_a_angle)
            
            data.append(tm1_a_angle)
            data.append(tm7_a_angle)
            
            return data
        
        elif tm == '1':
            
            tm1_a_angle = []


            
            for ts in u.trajectory:
            
                tm1_a_angle.append(get_tm1_a_angle(u))

            tm1_a_angle = np.array(tm1_a_angle)
            
            data.append(tm1_a_angle)
            
            return data
            
        elif tm == '7':

            tm7_a_angle = []
            
            for ts in u.trajectory:

                tm7_a_angle.append(get_tm7_a_angle(u))

        

            tm7_a_angle = np.array(tm7_a_angle)
            
            data.append(tm7_a_angle)
            
            return data
            
            
        else: 
            
            print('error. invalid tm input. leave blank or enter "1" or "7. exiting. ')
            

            



        
    elif protomer == ('b' or 'B'):
        
                
        if tm == 'both':
        
            tm1_b_angle = []

            tm7_b_angle = []
            
            for ts in u.trajectory:
            
                tm1_b_angle.append(get_tm1_a_angle(u))

                tm7_b_angle.append(get_tm7_a_angle(u))

            
            tm1_b_angle = np.array(tm1_b_angle)

            tm7_b_angle = np.array(tm7_b_angle)
            
            data.append(tm1_b_angle)
            data.append(tm7_b_angle)
            
            return data
            
        elif tm == '1':
            
            tm1_b_angle = []


            
            for ts in u.trajectory:
            
                tm1_b_angle.append(get_tm1_b_angle(u))

            tm1_b_angle = np.array(tm1_b_angle)
            
            data.append(tm1_b_angle)
            
            #return data
            
        elif tm == '7':

            tm7_b_angle = []
            
            for ts in u.trajectory:

                tm7_b_angle.append(get_tm7_b_angle(u))

        

            tm7_b_angle = np.array(tm7_b_angle)
            data.append(tm7_b_angle)
            
            return data
            
            
        else: 
            
            print('invalid tm input. leave blank or enter "1" or "7. exiting. ')
            
    else:
        print('invalid protomer input. leave blank or use "a" or "b"')
    
        
if __name__ == '__main__':
    
     get_pore_angles(topology,traj,protomer='both',tm='both')


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
