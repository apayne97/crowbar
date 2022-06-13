"""
Contains functions for iterating functions from plot, calc, and convert over a set of trajectories organized by a
'sys_dict'.

"""
import random

VERSION = '0.1.0'

def sys_dict_data_iterator(func, sys_dict, data_name, *args, **kwargs):
    for name, info in sys_dict.items():
        print(name)
        t = info['traj']
        sys_dict[name][data_name] = func(t, *args, **kwargs)
    return sys_dict


def sys_dict_data_convert_iterator(func, sys_dict, input_data_dict, output_data_name, *args, **kwargs):
    for name, info in sys_dict.items():
        print(name)
        data_dict = {}
        for param, value in input_data_dict.items():
            data_dict[param] = info[value]
        sys_dict[name][output_data_name] = func(*args, **data_dict, **kwargs)
    return sys_dict


def sys_dict_plot_iterator(func, sys_dict, data_name, selector=False, pdf_dir=False, pdf_tag=False,
                           update_layout_kwargs=False, update_xaxes_kwargs=False, update_yaxes_kwargs=False,
                           update_traces_kwargs=False, show_fig=False, **kwargs):
    for name, info in sys_dict.items():
        print(name)
        data = info[data_name]

        if selector:
            data = data[data['Label'].isin(selector)]

        if info.get('Plot Title'):
            title = info['Plot Title']
        else:
            title = f'{info["Title"]}, {info["State"]} Equilibrated with {info["Equilibration"]} and run for {info["Length"]}ns'
        fig = func(data, **kwargs)
        fig.update_layout(title=title)
        if update_layout_kwargs:
            fig.update_layout(**update_layout_kwargs)

        if update_xaxes_kwargs:
            fig.update_xaxes(**update_xaxes_kwargs)

        if update_yaxes_kwargs:
            fig.update_yaxes(**update_yaxes_kwargs)

        if update_traces_kwargs:
            fig.update_traces(**update_traces_kwargs)

        if show_fig:
            print('Showing figure...')
            fig.show()

        if pdf_dir and pdf_tag:
            filepath = f'{pdf_dir}/{name}_{pdf_tag}.pdf'
            print('Saving figure...')
            fig.write_image(filepath)
    return sys_dict

def return_sample_item(sys_dict:dict, item):
    """
    For grabbing a trajectory to test something real quick.

    :param sys_dict:
    :return:
    """
    sys, info = random.choice(list(sys_dict.items()))
    item = info.get(item)

    print(f'{info["Title"]}')
    return item


def return_sample_sys(sys_dict: dict):
    """
    For grabbing a trajectory to test something real quick.

    :param sys_dict:
    :return:
    """
    sys, info = random.choice(list(sys_dict.items()))
    print(f'{info["Title"]}')
    return info

def print_sys_dict(sys_dict):
    for sys, info in sys_dict.items():
        print_str = []
        for item in ['Sys', 'Length', 'CloneIDX', 'Plot Title']:
            print_str.append(f'{item}: {info[item]}')
        print('\t'.join(print_str))

def get_sys_names(sys_dict):
    sys_names = set([info['Sys'] for info in sys_dict.values()])
    sys_list = list(sys_names)
    sys_list.sort(reverse=True)
    return sys_list