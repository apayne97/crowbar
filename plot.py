"""
Contains functions for plotting data!

The standard is for these functions to return a figure.
"""

import plotly as pt
import plotly.express as px
import plotly.graph_objects as go
from functools import wraps

VERSION = '0.2.0'

def make_square(func):

    ## thanks to this link <https://stackoverflow.com/questions/1782843/python-decorator-handling-docstrings> for pointing me to this
    @wraps(func)
    def wrapper(*args, **kwargs):

        fig = func(*args, **kwargs)
        fig.update_layout(height=600, width=600)
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
        )
        return fig
    return wrapper

def dist_hist_wrapper(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        fig = func(*args, **kwargs)
        fig.update_layout(height=800, width=800, barmode='overlay')
        fig.update_yaxes(range=[0, 0.5])
        fig.update_traces(opacity=0.5, xbins={'size': 0.25})
        #     'opacity': 0.5,
        #     'xbins': {'size': 0.25}
        return fig
    return wrapper

def remove_silly_annotations(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        fig = func(*args, **kwargs)
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        return fig
    return wrapper

def plot_tseries(tseries):
    fig = px.line(tseries)
    return fig

def plot_total_rmsd(total_rmsd_df):
    """
    Meant to be iterated for each trajectory, plots both Chains and Total separately.

    :param total_rmsd_df:
    :return:
    """
    fig = px.line(total_rmsd_df, y='RMSD (Å)', x='Time (ns)', color='Label')
    return fig


@remove_silly_annotations
def plot_combined_df_lines(combined_df, value_name, label="Label", time="Time (ns)", facet="System"):
    """
    Uses combined long_df.

    :param combined_df:
    :param value_name:
    :param label:
    :param time:
    :param facet:
    :return:
    """

    fig = px.line(combined_df,
                  y=value_name,
                  x=time,
                  color=label,
                  facet_col=facet,)
    fig.update_yaxes(range=[0,2])

    return fig

@remove_silly_annotations
def plot_rmsd_histograms(combined_df):
    fig = px.histogram(combined_df, x="RMSD (Å)", color="Label", facet_col="System", barmode="overlay", nbins=50)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis1=dict(title="Frequency"))
    # fig.update_xaxes(range=[0,2])
    # fig.update_yaxes(range=[0, 100])
    return fig

@make_square
def plot_dihedrals(df):
    #combined_df = construct_chi12_df(t, selection)
    fig = px.scatter(df, x='Chi1', y='Chi2', facet_col='Res', color='Time', facet_col_wrap=2)
    fig.update_yaxes(range=[-185, 185])
    fig.update_xaxes(range=[-185, 185])
    return fig

@make_square
def plot_dihedral_by_chain(df, resname, chainids = [0, 1], plot_type = 'scatter'):
    x = f'{resname}_Chain {chainids[0]}'
    y = f'{resname}_Chain {chainids[1]}'
    if plot_type == 'scatter':
        fig = px.scatter(df,
                         x=x,
                         y=y,
                         color=df.index,
                         )
        fig.update_yaxes(range=[-185, 185],
                         nticks=20
        )
        fig.update_xaxes(range=[-185, 185],
                         nticks=20)
    elif plot_type == 'density_contour':
        fig = px.density_contour(df,
                                 x=x,
                                 y=y,
                            )
        
        # bin_dict = {'start': -180,'end':180,'size':15}
        # fig.update_traces(xbins=bin_dict, ybins=bin_dict)
        # fig.update_traces(contours_coloring='fill', colorscale='Viridis')   

    elif plot_type == 'density_heatmap':
        fig = px.density_heatmap(df,
                                 x=x,
                                 y=y,
                                 nbinsx=72,
                                 nbinsy=72,
                                 marginal_x="histogram",
                                 marginal_y="histogram",
                                 # range_x=[-180, 180],
                                 # range_y=[-180, 180],
                                 )
    # fig.update_layout(labels={  # replaces default labels by column name
    #                          "index": "Time (ns)"
    #                      })
    return fig

@dist_hist_wrapper
def plot_dist_from_long_df(long_df):
    fig = px.histogram(long_df,
                       x='Minimum Heavy Atom Distance (Å)',
                       facet_col='Label',
                       color='Chain',
                      histnorm='probability'
                      )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig

def plot_dihedral_by_chain_histogram(df, resname, chainids = [0, 1]):
    x = f'{resname}_Chain {chainids[0]}'
    y = f'{resname}_Chain {chainids[1]}'
    fig = px.density_heatmap(df,
                             x=x,
                             y=y,
                             nbinsx=72,
                             nbinsy=72,
                             marginal_x="histogram",
                             marginal_y="histogram",
                             range_x=[-180, 180],
                             range_y=[-180, 180],
                             )
    return fig

@remove_silly_annotations
@make_square
def plot_replicates_with_error(replicate_df, data_name, sys_list = ['Open Default', 'Open 10x', 'Closed Default', 'Closed 10x']):
    """
    Assumes you have bootstrapped the 'Upper Bound' and 'Lower Bound' of the data_name you would like to plot.
    Splits up the 'Clone' on the x-axis and facets by 'System'.

    :param replicate_df:
    :param data_name:
    :return:
    """
    fig = px.scatter(replicate_df, x='Clone', y=data_name, facet_col='System', text='N Samples',
                        category_orders={  # replaces default order by column name,
                            "System": sys_list
                        },
                   error_y="Upper Bound", error_y_minus="Lower Bound"
                        )
    fig.update_traces(textposition="bottom center")
    return fig

@remove_silly_annotations
def plot_systems_with_error(replicate_df, data_name, sys_list = ['Open Default', 'Open 10x', 'Closed Default', 'Closed 10x']):
    """
    Assumes you have bootstrapped the 'Upper Bound' and 'Lower Bound' of the data_name you would like to plot.
    'System' is plotted on the x axis, the y axis is the 'data_name' data.

    :param replicate_df:
    :param data_name:
    :return:
    """
    fig = px.scatter(replicate_df, x='System', y=data_name, text='N Samples',
                        category_orders={  # replaces default order by column name,
                            "System": sys_list
                        },
                   error_y="Upper Bound", error_y_minus="Lower Bound"
                        )
    fig.update_traces(textposition="bottom center")
    return fig

def plot_combined_dihedral_plots(df, resname, chainids = [0, 1]):
    x = f'{resname}_Chain {chainids[0]}'
    y = f'{resname}_Chain {chainids[1]}'
    fig = px.density_heatmap(df,
                             x=x,
                             y=y,
                             nbinsx=72,
                             nbinsy=72,
    #                          marginal_x="histogram",
    #                          marginal_y="histogram",
                                   color_continuous_scale='YlGnBu',
                             range_x=[-180, 180],
                             range_y=[-180, 180],
                             facet_col='Sys Name',
                                  category_orders={ # replaces default order by column name
    #                     "Dist Name": ["ILE46_0 to ILE271_1", "ILE46_0 to ILE271_0", "ILE46_0 to ILE46_1", "ILE271_0 to ILE271_1"],
                        "Sys Name": ["Open CHARMM-GUI", "Open CGUI 10x + 100ns bb", "Closed CHARMM-GUI", "Closed CGUI 10x + 100ns bb"]
                    }
                             )
    fig.update_layout(height=400, width=1200)
    fig.update_xaxes(range=[-135, 185], nticks=20)
    fig.update_yaxes(range=[-135, 185], nticks=20)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # fig.update_yaxes(
    #     scaleanchor="x",
    #     scaleratio=1,
    # )
    return fig

