"""
Contains functions for calculating values from trajectories. Mainly uses MDTraj and OpenMM.
"""

import mdtraj as md
import matplotlib.pyplot as plt
import simtk.unit as unit
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from simtk.openmm.app import PDBxFile, CharmmPsfFile
from math import pi
import numpy as np
import pandas as pd

VERSION = '0.0.1'


## RMSD and RMSF Calculations

def get_total_rmsf(t):
    idx = t.topology.select(f'protein and name CA')
    subset = t.atom_slice(idx)

    rmsf = md.rmsf(target=t,
                   reference=None,
                   frame=0,
                   atom_indices=idx)
    total_rmsf = rmsf.sum()
    return total_rmsf

def get_rmsf_by_chain(t):
    chain_dict = {'Chain A': 0, 'Chain B': 1}
    chain_rmsf_dict = {}

    ## make figure larger for easy reading
    plt.figure(figsize=(8, 6), dpi=80)
    plt.ylim([0, 10])
    for chain, chainid in chain_dict.items():
        idx = t.topology.select(f'protein and chainid {chainid} and name CA')
        subset = t.atom_slice(idx)

        rmsf = md.rmsf(target=t,
                       reference=None,
                       frame=0,
                       atom_indices=idx)

        resids = [id for id in range(30, subset.n_residues + 30)]

        ##Save as angstroms
        chain_rmsf_dict[chain] = rmsf * 10

## Dihedral angles
def construct_chi12_df(t, selection):
    """
    Returns a dataframe of the chi1 and chi2 angles of the desired selection.

    Uses MDTraj selection math, i.e. 'protein and (residue 46 or residue 271)'

    :param t:
    :param selection:
    :return:
    """

    ## First get the chi1 and chi2 dihedrals for the given selection
    df = pd.DataFrame()
    resSlice = t.atom_slice(t.topology.select(selection))
    chi1s = md.compute_chi1(resSlice)
    chi2s = md.compute_chi2(resSlice)
    #reslist = [f'{str(residue)}_{residue.chain.index}' for residue in resSlice.topology.residues]
    reslist = [f'{str(residue)}' for residue in resSlice.topology.residues]

    ## convert the output of the compute method into a dataframe
    chi1df = pd.DataFrame(chi1s[1] * 180 / 3.14, columns = res)
    chi2df = pd.DataFrame(chi2s[1] * 180 / 3.14, columns = reslist)

    ## combine the dataframes into one
    combined_df = pd.DataFrame()
    for i in range(len(reslist)):
        res = reslist[i]
        chi1s = chi1df[res]
        chi2s = chi2df[res]
        newdf = pd.DataFrame({'Time': chi1df.index, 'Res': res, 'Chi1': chi1s, 'Chi2': chi2s})
        combined_df = combined_df.append(newdf)
    return combined_df

def construct_chi2_df(t, selection):
    """
    Returns a dataframe of the chi1 and chi2 angles of the desired selection.

    Uses MDTraj selection math, i.e. 'protein and (residue 46 or residue 271)'

    :param t:
    :param selection:
    :return:
    """

    # First get the chi2 dihedrals for the given selection
    df = pd.DataFrame()
    resSlice = t.atom_slice(t.topology.select(selection))
    chi2s = md.compute_chi2(resSlice)
    reslist = [f'{str(residue)}_Chain {residue.chain.index}' for residue in resSlice.topology.residues]

    ## convert the output of the compute method into a dataframe
    chi2df = pd.DataFrame(chi2s[1] * 180 / np.pi, columns=reslist)
    return chi2df


def get_dist_df_from_idx_dict(t, names, idx_dict, scheme='closest', return_long_df=True):
    """
    Given an idx_dict as generated from convert.get_mdtraj_idx_dict_from_dist_array(),
    return a distance df.

    :param t:
    :param names:
    :param idx_dict:
    :param scheme:
    :param return_long_df:
    :return:
    """
    dfs = {}
    for chain, idx_array in idx_dict.items():
        contacts, residue_pairs = md.compute_contacts(t, contacts=idx_array, scheme=scheme)
        contacts_angstroms = contacts * 10
        df = pd.DataFrame(contacts_angstroms, columns=names)

        if return_long_df:
            df['Time'] = t.time
            newdf = pd.melt(df, id_vars='Time')
            newdf.columns = ['Time', 'Label', 'Minimum Heavy Atom Distance (Å)']
            newdf['Chain'] = chain
        else:
            newdf = pd.DataFrame(df, columns=names)
            newdf['Time'] = t.time

        dfs[chain] = newdf
    if return_long_df:
        return_obj = pd.concat(dfs.values())
    else:
        print('return chain dict')
        # return_df = pd.merge(dfs['ChainA'],dfs['ChainB'], on="Time")
        return_obj = dfs
    return return_obj