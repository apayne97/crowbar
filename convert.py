"""
Contains useful scripts for converting between various types of data.
i.e. includes methods to convert a dictionary of dataframes by chain into a dictionary of smaller dfs given a selection.
Doesn't require any unusual packages.
"""
import pandas as pd


VERSION = '0.0.1'


def pymol_to_mdtraj(pymol_string):
    """
    Simple function to convert a pymol-like residue string into mdtraj-like
    """
    x = pymol_string.replace('+', ' or residue ')
    y = x.replace('-', ' to ')
    final = f'residue {y}'
    return final

def get_chain_df_dict_from_large_chain_df_dict(chain_df_dict, single_selection):
    """
    Given a dictionary of dataframes by chain, this selects a subset of each dataframe
    based on the selection of interest.
    """
    new_chain_df_dict = {}
    chainids = range(len(chain_df_dict))
    chainletters = list(chain_df_dict.keys())

    dist0id = single_selection[0][0]
    dist0name = single_selection[0][1]
    dist1id = single_selection[1][0]
    dist1name = single_selection[1][1]
    time = chain_df_dict[chainletters[0]]['Time']

    for chainid in chainids:

        if chainid == 1:
            if dist0name == dist1name:
                continue
            dist0id = chainids[dist0id - 1]
            dist1id = chainids[dist1id - 1]

        dist0chainletter = chainletters[dist0id]
        dist1chainletter = chainletters[dist1id]

        if dist0name == dist1name and chainid == 0:
            series0name = f'{dist0chainletter}_{dist0name}'
            series1name = f'{dist1chainletter}_{dist1name}'
        else:
            series0name = f'{dist0name}'
            series1name = f'{dist1name}'

        series0 = chain_df_dict[dist0chainletter][dist0name]
        series1 = chain_df_dict[dist1chainletter][dist1name]

        new_chain_df_dict[chainletters[chainid]] = pd.DataFrame({'Time': time,
                                                                 series0name: series0,
                                                                 series1name: series1,
                                                                 })
    return new_chain_df_dict


def get_mdtraj_idx_array(t, dist_array):
    """
    Given a numpy array that looks like this:
    (
    ((XX, chainid), (YY, chainid)),
    ((XX, chainid), (YY, chainid)),
    )
    """
    dist_label_dict = {}

    ## Not super necessary but nice to get the names:
    protein = t.atom_slice(t.topology.select('chainid 0 or chainid 1'))
    chainA = t.atom_slice(t.topology.select('chainid 0'))
    n_residues = chainA.topology.n_residues

    ## this will returned
    names = []

    ## these will be return as a dict
    chainA_idx_list = []
    chainB_idx_list = []

    ## for each pair of residues in the dist_array
    for res1_set, res2_set in dist_array:

        ## get the naive res numbers
        res1_number = res1_set[0]
        res2_number = res2_set[0]

        ## for chain A, chain B
        for chain in [0, 1]:

            ## set res1
            res1_chainid = chain

            ## if it is zero, then we can take the chainid's as is
            if res1_chainid == 0:
                res2_chainid = res2_set[1]

                ## have to account for missing residues
                res1id = res1_number - 30

            ## otherwise, we have to flip the chainid
            elif res1_chainid == 1:
                res2_chainid = int(not res2_set[1])

                ## we also have to add n_residues to res1 to give it chain 0's number
                res1id = res1_number - 30 + n_residues

            ## now we have set res2's chainid, so we can set the resid
            if res2_chainid == 0:
                res2id = res2_number - 30
            elif res2_chainid == 1:
                res2id = res2_number - 30 + n_residues

            print('RES NUMBER, CHAIN:', res1_number, res1_chainid, res2_number, res2_chainid)

            print('RES ID:', res1id, res2id)

            ## get the residues!
            res1 = protein.topology.residue(res1id)
            res2 = protein.topology.residue(res2id)

            if res1_chainid == 0:
                ## this is a nice check just to make sure the chain index is the one we think we are using
                ## we will save the name based on the first residue being in chain 0
                string = f'{res1}_{res1.chain.index} to {res2}_{res2.chain.index}'
                names.append(string)

                chainA_idx_list.append([res1id, res2id])


            elif res1_chainid == 1:
                chainB_idx_list.append([res1id, res2id])

    ## now we combine the two lists

    idx_dict = {'ChainA': chainA_idx_list, 'ChainB': chainB_idx_list}

    return names, idx_dict


def get_dihedral_bin_probabilities_from_df(df, bin_boundaries = [-180, 0, 120, 180]):
    cols = list(df.columns)
    print(cols)
    bin_list = [(col, bin_max) for col in cols for bin_max in bin_boundaries]

    total_bins = len(cols) * (len(bin_boundaries) - 1)

    df_dict = {}

    for i in range(len(cols)):
        col_dict = {}
        for j in range(len(bin_boundaries) - 1):
            col = cols[i]
            lower = bin_boundaries[j]
            upper = bin_boundaries[j + 1]
            new_df = df[(df[col] > lower) & (df[col] < upper)]
            #return_df = new_df.loc[:, col]
            col_dict[f'{lower}>{col}<{upper}'] = new_df
        df_dict[col] = col_dict
        #print(df_dict)

    state_list = [f'{r1}_{r2}' for r1 in [-60, 60, 160] for r2 in [-60, 60, 160]]
    print(state_list)
    # for col in cols

    chaina = cols[0]
    chainb = cols[1]

    state_names = []
    counts = []
    for chaina_tuple in df_dict[chaina].items():
        for chainb_tuple in df_dict[chainb].items():
            chaina_index = list(chaina_tuple[1].index)
            chainb_index = list(chainb_tuple[1].index)
            index_list = [value for value in chaina_index if value in chainb_index]
            state_names.append(f'{chaina_tuple[0]}_{chainb_tuple[0]}')
            counts.append(len(index_list))

    prob_df = pd.DataFrame({'State Name': state_names, 'Counts': counts})

    return prob_df

