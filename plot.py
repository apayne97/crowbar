"""
Contains functions for plotting data!

The standard is for these functions to return a figure.
"""

import plotly as pt
import plotly.express as px
import plotly.graph_objects as go

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

    # fig.update_layout(
    #     xaxis=dict(
    #         tickmode='linear',
    #         tick0=-180,
    #         dtick=25,
    #     ),
    #     yaxis=dict(
    #         tickmode='linear',
    #         tick0=-180,
    #         dtick=25
    #     )
    # )
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

def plot_state_probabilities(prob_dict):
    pass