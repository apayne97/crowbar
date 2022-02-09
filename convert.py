"""
Contains useful scripts for converting between various types of data.
i.e. includes methods to convert a dictionary of dataframes by chain into a dictionary of smaller dfs given a selection.
Doesn't require any unusual packages.
"""
import pandas as pd
from pymbar import timeseries
import numpy as np


VERSION = '0.2.5'


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


def get_mdtraj_idx_dict_from_dist_array(t, dist_array):
    """
    Given a numpy array that looks like this:
    (
    ((XX, chainid), (YY, chainid)),
    ((XX, chainid), (YY, chainid)),
    )

    This function will return a set of dist_names and mdtraj selection strings like this:
    dist_names = ['ILE46_0 to ILE271_0', 'ILE46_0 to ILE271_1']
    mdtraj_dict = {'ChainA': [[16, 241], [16, 688]], 'ChainB': [[463, 688], [463, 241]]}
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


def get_state_df_from_dims_and_bin_bounds(dims, bin_bounds):
    import itertools

    ## bins include left and right outer bounds
    n_bins = len(bin_bounds) - 1
    n_dims = len(dims)

    ## the number of states scales as the number of bins to the power of the number of dimensions
    n_states = n_bins ** n_dims

    ## Get bins from bin_bounds
    bin_list = [(bin_bounds[i], bin_bounds[i + 1]) for i in range(n_bins)]

    ## This is just a list of the new state ids
    state_idx = [i for i in range(n_states)]

    ## This gets all the combinations (where order matters) of the bins for each dimension
    state_bins = list(itertools.product(bin_list, repeat=n_dims))

    ## This needs to be true!
    assert len(state_idx) == len(state_bins)

    ## Love this dict comprehension
    ## this unzips the list of tuples from the state_bins into two lists
    dim_bin_dict = {dim: list(zip(*state_bins))[i] for i, dim in enumerate(dims)}

    ## Probably not necessary but I think its kinda nice to have?
    dim_bin_dict['State'] = state_idx

    ## DataFrames are nice so lets return that
    df = pd.DataFrame(dim_bin_dict)

    return df

def get_state_from_snapshot(snapshot, state_df):
    """
    For each snapshot in time, captures which state the snapshot is in.

    Notably, does not require a specific number of states or dimensions.

    But is probably very slow for many dimensions / states.

    State bins are assumed to be [lower, upper)

    :param snapshot:
    :param state_df:
    :return:
    """

    ## start with no state label

    state_label = None
    n_state = len(state_df)

    ## Now we iterate through each state and see if it is in that state
    for n in range(n_state):
        ## get our state definitions from the state_df
        state_def = list(state_df.iloc[n])
        state_idx = state_def[-1]
        state_bins = state_def[:-1]

        ## Assume true, then mark false and break if any of the dimensions is not in the state definition
        ## I think this is faster because it will quit faster
        this_label = True
        for i, dim_bin in enumerate(state_bins):
            lower, upper = dim_bin
            data = snapshot[i]
            if not data >= lower or not data < upper:
                this_label = False
                break

        ## if all necessary conditions have still been met, then call this snapshot this state and quit
        if this_label == True:
            state_label = state_idx
            break

    return state_label

def get_state_timeseries_from_df(df, state_df):
    """
    Assumes 'State' is last column in state_df.

    Also assumes that it needs to index the df based on dimensions in state_df.

    :param df:
    :param state_df:
    :return:
    """

    dims = list(state_df.columns[:-1])

    data_df = df[dims]

    result = data_df.apply(get_state_from_snapshot, axis=1, state_df=state_df)

    return result

def compress_state_space(tseries, state_df):
    """
    Converts result from get_state_timeseries_from_df into a dataframe and then compresses into the state space of the provided state_df

    :param tseries:
    :param state_df:
    :return:
    """

    tseries_df = pd.DataFrame({'All State': tseries})

    two_state_series = tseries_df.apply(get_state_from_snapshot, axis=1, state_df=state_df)

    return two_state_series

def get_uncorrelated_tseries(tseries):
    """

    :param tseries:
    :return:
    """
    try:
        idx = timeseries.subsampleCorrelatedData(tseries)
        uncorr = tseries[idx]
    except timeseries.ParameterError:
        print("Failed to find autocorrelation time")
        uncorr = tseries

    return uncorr


def get_state_prob_from_tseries(tseries, state_dict):
    """
    Based on a dictionary of state names and idx, convert a tseries to a dataframe of probabilities associated with that state.


    :param state_dict:
    :return:
    """

    probs = []
    for name, idx in state_dict.items():
        prob = len(tseries[tseries == idx]) / len(tseries)
        probs.append(prob)

    df = pd.DataFrame({'State': list(state_dict.keys()), 'Probability': probs})

    return df

def bootstrap_error_bars(tseries, n=10000, state_idx=1):
    """
    Assumes binary data.
    """

    lower_bound = int(n * 0.025)
    upper_bound = int(n * 0.975)

    probs = []
    for i in range(n):
        sample = tseries.sample(frac=1, replace=True, axis=0)
        prob = len(sample[sample == state_idx]) / len(sample)
        probs.append(prob)

    probs.sort()

    lower = probs[lower_bound]
    upper = probs[upper_bound]

    mean = np.mean(probs)
    print(lower, mean, upper)
    print(mean-lower, mean, upper-mean)

    return (mean-lower, mean, upper-mean)

def get_binary_state_prob_from_tseries(tseries):
    """
    Based on a dictionary of state names and idx, convert a tseries to a dataframe of probabilities associated with that state.


    :param state_dict:
    :return:
    """

    lower, mean, upper = bootstrap_error_bars(tseries)

    df = pd.DataFrame({'State': 'Open', 'Probability': mean, 'Lower Bound': lower, 'Upper Bound': upper},index=[0])

    return df






def get_dihedral_bin_probabilities_from_df(sys_name, df):
    cols = list(df.columns)[:-1] ## This drops the 'sys name' column
    print(cols)
    bin_list = [(col, bin_max) for col in cols for bin_max in bin_boundaries]

    total_bins = len(cols) * (len(bin_boundaries) - 1)

    total_len = len(df)

    df_dict = {}

    for i in range(len(cols)):
        col_dict = {}
        for j in range(len(bin_boundaries) - 1):
            col = cols[i]
            lower = bin_boundaries[j]
            upper = bin_boundaries[j + 1]
            # print(upper, lower)
            # print(df[col])
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
    probabilities = []
    for chaina_tuple in df_dict[chaina].items():
        for chainb_tuple in df_dict[chainb].items():
            chaina_index = list(chaina_tuple[1].index)
            chainb_index = list(chainb_tuple[1].index)
            index_list = [value for value in chaina_index if value in chainb_index]
            state_names.append(f'{chaina_tuple[0]}_{chainb_tuple[0]}')

            count = len(index_list)
            prob = count / total_len

            counts.append(count)
            probabilities.append(prob)

    state_names = ['Closed-Closed',
                   'Closed-Intermediate',
                   'Closed-Open',
                   'Closed-Intermediate',
                   'Intermediate-Intermediate',
                   'Intermediate-Open',
                   'Closed-Open',
                   'Intermediate-Open',
                   'Open-Open'
                   ]

    prob_df = pd.DataFrame({'State Name': state_names, 'Counts': counts, 'Probability': probabilities, 'Sys Name': sys_name})

    return prob_df

def get_mean_from_long_dist_df(sys_name, long_df, data_name):
    """
    Collapses data split by chain into a mean.
    """

    ## get list of distance names
    dist_names = list(set(long_df['Label']))

    means = []

    for dist_name in dist_names:
        concat_df = long_df[long_df['Label'] == dist_name]
        mean = concat_df[data_name].mean()
        means.append(mean)

    df = pd.DataFrame({'Dist Name': dist_names, 'Mean Val (Å)': means, 'Sys Name': sys_name})

    return df

def get_replicate_df(sys_dict, df_name):
    """
    Concatenates similar dfs in a sys_dict.

    :param sys_dict:
    :param df_name:
    :return:
    """

    df_list = []
    idx_list = []
    for sys, info in sys_dict.items():
        idx_list.append(info['CloneIDX'])
        df = info[df_name]
        sys_name = info['Sys']
        df['Sys Name'] = sys_name
        df_list.append(df)
    full_df = pd.concat(df_list)
    full_df['Clone ID'] = idx_list
    return full_df