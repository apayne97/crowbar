#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################
############################################################

# Author: Viktor Belay 
# Creation date: 02/06/2022

# Purpose of script: Save CHAP data as numpy arrays and plot things in a pretty way :3
# Also, gets output files

# Version: 0.1.0

# Funtionalty:

# Future plans: 

############################################################
############################################################

import numpy as np
from matplotlib import pyplot as pl
import matplotlib.mlab as mlab
import json
from scipy.stats import norm
import mdtraj as md
import os
import sys
import subprocess
import getpass

def do_chap(topology,trajectory,output_path,solvent='15'):
    from scp import SCPClient, SCPException
    import paramiko
    
    def dcd_to_xtc(topology, trajectory):
        
        dcd_traj = md.load(trajectory,top=topology)
        dcd_traj.save_xtc(output_path+'/equilibrated.xtc')
        
    host = 'hitegpu2'
    username = 'hiter'
    password = getpass.getpass("Enter password for "+username+'@'+host+":")
        
    port = 22 
        
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
    ssh.connect(host,port,username=username,password=password); x=1
        
    if x == 1:
            
        print('Connection to '+username+'@'+host+' successful')
        
        if trajectory.endswith('.dcd'):
            print('Converting DCD to XTC...')
            
            dcd_to_xtc(topology, trajectory)
            
        elif trajectory.endswith('.xtc') == False:
            
            print('Please select a trajecory in either xtc or dcd format. Exiting.')
            
        else:
            
            print('Input trajectory is xtc, no conversion necessary.')
            
            
        ############################################################    
        # Create a temp dir in hitegpu2 for CHAP execution
        stdin,stdout,stderr=ssh.exec_command('cd belayv ; mkdir vb_temp_chap_dir ;')
        outlines1 = stdout.readlines()
        outlines2 = stderr.readlines()
        resp1=''.join(outlines1)
        resp2=''.join(outlines2)
        print(resp1)
        print(resp2)
        
        ############################################################     
            
        
        ############################################################
        # Create SCP client and copy traj and topology to hitegpu
        
        traj_path = output_path+'/'+'equilibrated.xtc'
        
        scp=SCPClient(ssh.get_transport())
        scp.put(traj_path,remote_path='/home/hiter/belayv/vb_temp_chap_dir/.')
        scp.put(topology,remote_path="/home/hiter/belayv/vb_temp_chap_dir/.")
        scp.close()
        print('Trajectory and topology have been temporarily copied to hitegpu2.')

        
        
        
        ############################################################



        ############################################################
        # Run CHAP execution
        
        print('Starting CHAP execution...')
        
        stdin,stdout,stderr=ssh.exec_command('cd belayv;cd vb_temp_chap_dir;source /programs/sbgrid.shrc;/programs/x86_64-linux/chap/0.9.1/bin.capsules/chap -f equilibrated.xtc -s step5_input.pdb -sel-pathway 1 -sel-solvent 15 -hydrophob-fallback -0.56846473')
        outlines3 = stdout.readlines()
        outlines4=stderr.readlines()
        resp3=''.join(outlines3)
        resp4=''.join(outlines4)
        print(resp3)
        print(resp4)
        print('CHAP execution successful.')
        
        
        ############################################################
        
        
        ############################################################
        # Copy output data from hitegpu2 back to original output dir
        

        scp = SCPClient(ssh.get_transport())
        os.chdir(output_path)
        
        scp.get('/home/hiter/belayv/vb_temp_chap_dir/output.json')
        scp.get('/home/hiter/belayv/vb_temp_chap_dir/output.pdb')
        scp.get('/home/hiter/belayv/vb_temp_chap_dir/output.mtl')
        
        scp.close()
        
        print('Data transferred from hitegpu2 successfully.')
        
        
        ############################################################



        ############################################################
        # Delete temp folder from hitehpu2
        
        stdin,stdout,stderr = ssh.exec_command('cd belayv;rm -r vb_temp_chap_dir')
        outlines5=stdout.readline()
        resp5=''.join(outlines5)
        
        print(resp5)
        
        ############################################################

    ###
    # Close SSH clinet
    ssh.close()
    print('SSH connection to '+username+'@'+host+' successfully closed.')
    
    
    ###


def load_data(data):
    
    with open(data) as data_file:
        chap_data=json.load(data_file)
        
    return chap_data

    
        
def get_pore_radius_profile(chap_dat,plot='True'):
    
    pass

def get_min_pore_radius(chap_data):
    dat=[]
    min_radius = np.array(chap_data['pathwayScalarTimeSeries']['minRadius'])*10 # Min pore (nm) radians * 10 = min pore radius in A 
    t = np.array(chap_data['pathwayScalarTimeSeries']['t'])
    
    dat.append(t)
    dat.append(min_radius)
    
    dat=np.array(dat)
    
    return dat

def plot_min_pore_radius(chap_data,hist=True,timeseries=True):
    min_radius = np.array(chap_data['pathwayScalarTimeSeries']['minRadius'])*10 # Min pore (nm) radians * 10 = min pore radius in A 
    t = np.array(chap_data['pathwayScalarTimeSeries']['t'])
    if hist == True:
        
        if timeseries==True:
            
            pl.figure()
            
            mean,std=norm.fit(min_radius[min_radius>=0])
            pl.hist(min_radius[min_radius>=0],bins=(int(len(min_radius)/4)),alpha=0.9)
            pl.xlabel('Minimum pore radius (A)')
            pl.ylabel('Frequency')
            pl.title('mean= '+str(mean)+' SD= '+str(std))
            
            pl.figure()
            
            pl.plot(t,min_radius)
            pl.xlabel('time (ns)')
            pl.ylabel('minimum pore radius (A)')
            
            pl.ylim(0,max(min_radius)+2)
            
        elif plot_time_series==False:
            
            pl.figure()
            
            mean,std=norm.fit(min_radius[min_radius>=0])
            pl.hist(min_radius[min_radius>=0],bins=(int(len(min_radius)/4)),alpha=0.9)
            pl.xlabel('Minimum pore radius (A)')
            pl.ylabel('Frequency')
            pl.title('mean= '+str(mean)+' SD= '+str(std))
            
    elif plot_hist==False:
        
        if plot_time_series==True:
            
            pl.figure()
            
            pl.plot(t,min_radius)
            pl.xlabel('time (ns)')
            pl.ylabel('minimum pore radius (A)')
            
            pl.ylim(0,max(min_radius)+2)
            
        elif plot_time_series==False:
            
            pass
        

def get_open_probability(chap_data,radius=1.33):
        
        min_radius = np.array(chap_data['pathwayScalarTimeSeries']['minRadius'])*10

        open_probability = (len(min_radius[min_radius>=radius]))/(len(min_radius))

        return open_probability
        
def get_pore_solvent_density(chap_data):
    
    dat = []
    z = np.array(chap_data['pathwayProfile']['s'])*10
    mean_solvent_density = np.array(chap_data['pathwayProfile']['densityMean'])
    
    dat.append(z)
    dat.append(mean_solvent_density)
    
    dat = np.array(dat)

    return dat
    
def plot_get_pore_solvent_desnity(chap_data):
    
    z = np.array(chap_data['pathwayProfile']['s'])*10
    mean_solvent_density = np.array(chap_data['pathwayProfile']['densityMean'])
    pl.figure()
    pl.plot(z,mean_solvent_density,'k')
    pl.vlines(3,0,100,'k',linestyles='--',alpha=0.4)
    pl.vlines(-4,0,100,'k',linestyles='--',alpha=0.4)
    pl.axvspan(z[500], z[580], alpha=0.3, color='orange')
    pl.xlim(-4,6)
    pl.ylim(0,100)

    pl.xlabel('z (nm)')
    pl.ylabel('Mean water density (nm^-3)')
    pl.title('mean water density BB open clone00')



    
    
    
    
    
    
        
        
    
    
    
    
    
    
    
    

