# crowbar
Analysis scripts and tools for simulating and analyzing proteins, with special focus on handling membranes and multiple chains.

# `iterate.py` and the sys_dict
With the idea that several trajectories are going to be analyzed all at once, I have written the `iterate` functions to use the 'sys-dict'.
This is a dictionary of trajectory dictionaries. It could have been its own object but that seemed like too much at the time.
The iterate scripts take the `sys_dict` as an argument and return it with a new item added to each traj dict. 
They can also take a plotting function and make a plot for a particular item in the traj dictionaries.

# Subdivision of scripts
munge.py - not transferred yet but I've used this to 'munge' existing trajectories on our HPC into single long trajectories
load.py - loads the sys_dict. currently loads the trajectories as mdtraj trajectories.
calc.py - functions for calculating features from the trajectory data
convert.py - functions for converting one datatype to another - usually having to do with taking information from several chains or residues and returning a new desired dataframe
plot.py - functions for plotting data! includes some useful wrappers as well!
iterate.py - discussed above, used to iterate a function through the sys_dict




