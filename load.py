import os.path
import yaml
import pandas as pd
import mdtraj as md
from simtk.openmm.app import CharmmPsfFile

VERSION = '0.2.1'

def load_sysdf(simulation_table_path="simulations.yaml"):
    """
    loads the simulations yaml file and returns a transformed dataframe
    :param simulation_table_path:
    :return:
    """
    with open(simulation_table_path, 'r') as file: sim_dict = yaml.safe_load(file)
    sysdf = pd.DataFrame(sim_dict).T
    return sysdf


def import_traj_from_munged(traj_path, psf_path, pdb_path):
    """
    Assumes trajectory strided to 1ns / frame and aligned with md.superpose().

    Returns an mdtraj trajectory object and a pdb file
    """

    ## load dcd into mdtraj object
    t = md.load_dcd(traj_path,
                    top=psf_path,
                    stride=1)

    pdb = md.load(pdb_path)

    ## get topology from psf file using openmm
    psf = CharmmPsfFile(psf_path)

    ## convert topology to mdtraj topology and save to mdtraj object
    t.topology = t.topology.from_openmm(psf.topology)
    pdb.topology = pdb.topology.from_openmm(psf.topology)

    return (t, pdb)

def load_test():
    crow_path = '/Users/alexpayne/Scientific_Projects/crowbar'
    traj_path = crow_path + '/testing/step7_1.dcd'
    pdb_path = crow_path + '/testing/step5_input.pdb'
    psf_path = crow_path + '/testing/step5_input.psf'
    t, pdb = import_traj_from_munged(traj_path, psf_path, pdb_path)
    return (t, pdb)


def import_systems_from_munged(sys_names, traj_prefixes, sim_yaml, prefix_dict, traj_file_extension, n_clones, munged_path):
    """
    Assumes trajectory strided to 1ns / frame and aligned with md.superpose().

    Returns the sys_dict. Trajectories are mdtraj trajectories.

    :param sys_names:
    :param traj_prefixes:
    :param sim_yaml:
    :param prefix_dict:
    :param traj_file_extension:
    :param n_clones:
    :param munged_path:
    :return sys_dict:
    """


    ## Make sure everything is as we expect
    assert type(sys_names) == list
    assert type(n_clones) == int
    assert type(traj_prefixes) == list

    ## Load the dictionary of simulations
    with open(sim_yaml, 'r') as file:
        sim_dict = yaml.safe_load(file)
    print(sim_dict)

    ## Load a dictionary of what the file name prefixes mean
    with open(prefix_dict, 'r') as file:
        prefix_dict = yaml.safe_load(file)
    print(prefix_dict)

    sys_dict = {}
    for system in sys_names:
        sys_info = sim_dict[system]
        # clone_dict = {}
        for idx in range(n_clones):
            for traj_prefix in traj_prefixes:
                clone_info = {}
                clone_info.update(sys_info)
                clone = f'{traj_prefix}_clone{idx:02d}'

                traj_path = f'{munged_path}{system}/{clone}.{traj_file_extension}'
                psf_path = f'{munged_path}{system}/step5_input.psf'
                pdb_path = f'{munged_path}{system}/step5_input.pdb'

                assert os.path.exists(traj_path), f'{traj_path} does not exist!'
                assert os.path.exists(psf_path), f'{psf_path} does not exist!'

                clone_info['Length'] = prefix_dict[traj_prefix]['Length']

                full_name = f'{system}_{clone}'
                clone_info['Title'] = full_name

                ## Save mdtraj trajectory object of the loaded pdb path
                pdb_path = f'{munged_path}{system}/step5_input.pdb'

                print(f'Loading {full_name} from {traj_path}')

                traj, pdb = import_traj_from_munged(traj_path, psf_path, pdb_path)
                clone_info['Input PDB'] = pdb
                clone_info['traj'] = traj

                ## this values is used to combine replicates (clones) of the same system
                clone_info['Sys'] = f'{clone_info["State"]} {clone_info["Equilibration"]}'

                sys_dict[full_name] = clone_info

    return sys_dict

def load_dist_dict(dist_yaml='/Users/alexpayne/Scientific_Projects/crowbar/testing/distances.yaml'):
    print(f'Getting dist_dict from {dist_yaml}')

    with open(dist_yaml) as f:
        dist_dict = yaml.safe_load(f)

    return dist_dict
