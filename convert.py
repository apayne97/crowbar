"""
Contains useful scripts for converting between various types of data.
i.e. includes methods to convert a dictionary of dataframes by chain into a dictionary of smaller dfs given a selection.
Doesn't require any unusual packages.
"""
import pandas as pd


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

def get_dihedral_bin_probabilities(df, bin_boundaries = [-180, 0, 120, 180]):
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

    prob_df = pd.DataFrame('State Name' = state_names, 'Counts' = counts)

    return prob_df

