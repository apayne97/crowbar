"""
Contains functions for plotting data!

The standard is for these functions to return a figure.
"""

import plotly as pt
import plotly.express as px
import plotly.graph_objects as go

VERSION = '0.2.0'

def make_square(func):
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
    def wrapper(*args, **kwargs):
        fig = func(*args, **kwargs)
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        return fig
    return wrapper

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
                         color=df.index)
        fig.update_yaxes(range=[-185, 185], nticks=20)
        fig.update_xaxes(range=[-185, 185], nticks=20)
    elif plot_type == 'contour':
        fig = px.density_heatmap(df,
                                 x=x,
                                 y=y,
                                 marginal_x="histogram",
                                 marginal_y="histogram")
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
def plot_distance_replicate_df(replicate_df):
    fig = px.strip(replicate_df, x='Dist Name', y='Mean Val (Å)', color='Sys Name',
                        category_orders={ # replaces default order by column name
                    "Dist Name": ["ILE46_0 to ILE271_1", "ILE46_0 to ILE271_0", "ILE46_0 to ILE46_1", "ILE271_0 to ILE271_1"],
                    "Sys Name": ["Open CHARMM-GUI", "Open CGUI 10x + 100ns bb", "Closed CHARMM-GUI", "Closed CGUI 10x + 100ns bb"]
                })
    fig.update_xaxes(title='Residue-Residue Distance')
    fig.update_yaxes(title='Mean CA-CA Distance (Å)')
    return fig

@remove_silly_annotations
def plot_dihedral_prob_replicate_df(replicate_df):
    fig = px.strip(replicate_df, x='Sys Name', y='Probability', facet_col='State', color='Sys Name',
                        category_orders={  # replaces default order by column name,
                            "Sys Name": ["Open CHARMM-GUI", "Open CGUI 10x + 100ns bb", "Closed CHARMM-GUI",
                                         "Closed CGUI 10x + 100ns bb"]
                        },
                        )
    return fig

@remove_silly_annotations
def plot_dihedral_prob_replicate_df_with_error(replicate_df):
    fig = px.scatter(replicate_df, x='Clone ID', y='Probability', facet_col='Sys Name',
                        category_orders={  # replaces default order by column name,
                            "Sys Name": ["Open CHARMM-GUI", "Open CGUI 10x + 100ns bb", "Closed CHARMM-GUI",
                                         "Closed CGUI 10x + 100ns bb"]
                        },
                   error_y="Upper Bound", error_y_minus="Lower Bound"
                        )
    return fig

@remove_silly_annotations
def plot_dihedral_prob_replicate_df_with_error_combined(replicate_df):
    fig = px.scatter(replicate_df, x='Sys Name', y='Probability',
                        category_orders={  # replaces default order by column name,
                            "Sys Name": ["Open CHARMM-GUI", "Open CGUI 10x + 100ns bb", "Closed CHARMM-GUI",
                                         "Closed CGUI 10x + 100ns bb"]
                        },
                   error_y="Upper Bound", error_y_minus="Lower Bound"
                        )
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