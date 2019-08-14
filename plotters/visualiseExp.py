'''
Created on 7 Oct 2016

@author: shsymoen
'''

from plotters.generalPlotters import colorbar_index, f_save
from plotters.comparison import (
    create_hor_validation_lines,
    plot_errBar_exp
)
from dataProcessing.utils import (
    conf_int_dist_chi2,
    convert_smiles_sort,
    inchi_2_smiles,
    popup_error,
)


def annotate_tsne_plot(ax, X_embed, cmap, cmap_labels, colors, ann_texts):
    '''Puts annotations next to the datapoints in the t-SNE representation
    of the data

    Parameters
    ----------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    X_embed : ndarray
        t-SNE embedded data
    cmap : matplotlib cmap
    cmap_labels : lst
        List of labels used for the colormap
    colors : lst, ndarray
        List of the datapoints with a color encoding value
    ann_texts : lst, ndarray
        List of annotation text for plotting

    Returns
    -------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    from sklearn.preprocessing import MinMaxScaler

    # Scale the data in the same way as the t-SNE representation
    mm_sc = MinMaxScaler()
    X_embed = mm_sc.fit_transform(X_embed)
    X, Y = X_embed[:, 0], X_embed[:, 1]
    # Add annotations to the datapoints
    for c, xi, yi, txt in zip(colors, X, Y, ann_texts):
        cmap_index = cmap_labels.index(c)
        ax.annotate(txt, (xi, yi), color=cmap(cmap_index))

    return ax


def dbscan_color_encoding(X_embed,):
    '''Scales the t-SNE input data with a standardscaler prior to performing
    the DBSCAN clustering algorithm and returning the result

    Parameters
    ----------
    X_embed : ndarray
        Embedding of the training data (X) in low-dimensional space.

    Returns
    -------
    db : ndarray
        List of which point belongs to which cluster
    '''
    from sklearn.cluster import DBSCAN
    from sklearn.preprocessing import StandardScaler

    # DBSCAN for cluster based coloring (first scale the data)
    std_sc = StandardScaler()
    X_embed_sc = std_sc.fit_transform(X_embed)
    db = DBSCAN(eps=0.3, min_samples=3).fit_predict(X_embed_sc)

    return db


def ellips_sub(
    x,
    y,
):
    '''A subroutine for the calculation of the eigenvalues and
    eigenvectors

    Parameters
    ----------
    x : ndarray
    y : ndarray

    Returns
    -------
    v : ndarray
    w : ndarray
    xy : ndarray
    '''
    import numpy as np
    import pandas as pd

    # Put x, y in a pandas DataFrame
    subdata = pd.DataFrame({'x': x, 'y': y})

    # Calculate the mean center
    xy = subdata.mean().values

    # Calculate the covariance matrix
    cov_mat = subdata.cov()

    # Compute the eigenvalues and eigenvectors
    # of the covariance matrix
    v, w = np.linalg.eig(cov_mat)

    return v, w, xy


def ellips_parameters(
    v,
    w,
    xy,
    conf_int_dist=5.991,
    verbose=False
):
    '''Creates the parameters necessary to draw the
    confidence interval Mahalanobis distance.

    http://www.visiondummy.com/2014/04/
    draw-error-ellipse-representing-covariance-matrix/

    Parameters
    ----------
    x : list, ndarray
        x-values of the data
    y : list, ndrray
        y-values of the data
    conf_int_dist : float, optional (default: 5.991
        Confidence interval distance measure
    verbose : boolean, optional (default: False)
        print some usefull information to stdout

    Returns
    -------
    ellips_paras : dict
        xy : ndarray
            The coordinates of the Ellips center
            [x, y]
        width : float
            The width of the Ellips
        height : float
            The height of the Ellips
        angle : float
            The angle in degrees counter-clockwise relative
            to the y-axis
    '''
    import math
    import numpy as np

    # Get the index of the min, max eigenvalue
    id_max = np.argmax(v)
    id_min = np.argmin(v)

    # Calculate the angle for the ellips in radians
    angle_rad = np.arctan2(w[id_max][1], w[id_max][0])
    # angle counter clockwise starting vertical
    # convert from radians to degrees
    angle = 90 - math.degrees(angle_rad)

    # Calculate the height
    height = 2 * np.sqrt(conf_int_dist * v[id_max])

    # Calculate the width
    width = 2 * np.sqrt(conf_int_dist * v[id_min])

    # create the output dict
    ellips_paras = {
        'xy': xy,
        'angle': angle,
        'height': height,
        'width': width
    }

    # print some information if verbose is True
    if verbose:
        # print('cov_matrix\n', cov_mat)
        print('eigenvalues\n', v, '\neigenvectors\n', w, '\n')
        print('center:', xy)
        print('height:', height)
        print('width:', width)
        print('angle:', angle)

    return ellips_paras


def create_ellips_data(data, conf_int_dists, groupby='temperature'):
    '''Creates the ellips data DataFrame that can be used to plot
    the Mahalanobis distances. It groups the data by groupby and
    creates a number of ellipses for the number of confidance
    levels in conf_lvls

    Parameters
    ----------
    data : DataFrame
        a column 'groupby' element
        a column x and y
    conf_int_dists : list
        List of confidence interval distances
    groupby : string, optional (default: temperature)
        Element to groupby

    Returns
    -------
    ellips_data : DataFrame
    '''
    import pandas as pd

    # Initialisations
    ellips_data = {}

    # Group the data based on the same groupby element
    grouped = data.groupby(groupby)

    # Iterate over all groups
    for group_name in grouped.groups:
        # Select x and y column
        subdata = grouped.get_group(group_name)[['x', 'y']]

        # If the grouped data only contains 1 value, skip MD creation
        if len(subdata) > 1:

            # Compute the eigenvalues and eigenvectors
            v, w, xy = ellips_sub(
                x=subdata['x'].values,
                y=subdata['y'].values,
            )

            # Calculate the ellips parameters
            ellips_paras = {}

            for conf_int_dist in conf_int_dists:
                ellips_paras[conf_int_dist] = ellips_parameters(
                    v=v,
                    w=w,
                    xy=xy,
                    conf_int_dist=conf_int_dist,
                )
                # Add to ellips_data as DataFrame
                ellips_data[group_name] = pd.DataFrame(ellips_paras)

    # Create a DataFrame with the input data for the ellipses
    ellips_data = pd.concat(
        ellips_data,
        axis=1,
        names=['Temperature', 'conf_int']
    )

    return ellips_data


def ellipses_plotter(ax, ellips_data, conf_lvls):
    '''Creates the matplotlib Ellips patches with the
    xy, width, height, and angle information from the
    ellips_data DataFrame

    Parameters
    ----------
    ax : matplotlib axes
    ellips_data : DataFrame
        A dataframe that contains the xy, width, height and angle
        information for the different ellipses.
        The DataFrame contains a multiindex column
    conf_lvls : list
        A list of the confidence levels

    Returns
    -------
    ax : matplotlib.axes.Axes object

    Creates
    -------
    Ellipses on a given axes
    '''
    from matplotlib.patches import Ellipse  # Ellips patches
    import matplotlib.lines as mlines
    import numpy as np
    import pandas as pd
    idx = pd.IndexSlice  # MultiIndex slicing

    # Groubpy element
    temperatures = ellips_data.columns.get_level_values(0).unique()

    # Conf_intervals
    conf_int_dists = ellips_data.columns.get_level_values(1).unique()

    # Initialisations
    line_handles = []
    colors_single = ['Green', 'Orange', 'Red']
    ellipses = np.array([
        [
            Ellipse(
                xy=ellips_data.loc['xy', idx[T, conf_int_dist]],
                width=ellips_data.loc['width', idx[T, conf_int_dist]],
                height=ellips_data.loc['height', idx[T, conf_int_dist]],
                angle=ellips_data.loc['angle', idx[T, conf_int_dist]]
            )
            for conf_int_dist
            in conf_int_dists
        ] for T in temperatures
    ])
    ellipses = ellipses.ravel()

    # Add the Ellips patches to the figure
    colors = len(temperatures) * colors_single

    # Iteration
    for ellips, color in zip(ellipses, colors):
        ax.add_artist(ellips)
        ellips.set_facecolor("None")
        ellips.set_edgecolor(color)
        ellips.set_alpha(0.5)
        ellips.set_linestyle('dashed')
        line_handles.append(
            mlines.Line2D([], [], color=color, ls='dashed')
        )

    # Create labels
    lbls = ['{} %'.format(conf_lvl * 100) for conf_lvl in conf_lvls]

    # Add legend
    ax.legend(handles=line_handles[:3], loc='best', labels=lbls)


def create_scatterplot_with_annotations_colorbar(
    ax,
    x,
    y,
    c,
    cmap,
    cmap_labels,
    cbar_label,
    ann_text,
    alpha=1,
):
    '''Creates a scatterplot with annotated points and adds a colorbar

    Parameters
    ----------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    x : ndarray
    y : ndarray
    c : lst, ndarray
    cmap : matplotlib cmap
    cmap_labels : lst
        List of labels used for the colormap
    cbar_label : string
    ann_text : lst, ndarray
        List of annotation text for plotting
    alpha : float
        The opacity level of the scatterpoints and annotated text

    Returns
    -------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    from sklearn import preprocessing

    # Create an encoding based on c to color the scatterplot points
    le = preprocessing.LabelEncoder()
    cc = le.fit_transform(c)

    ax.scatter(x=x, y=y, c=cc, cmap=cmap, alpha=alpha)
    colorbar_index(
        ncolors=cmap.N,
        cmap=cmap,
        tick_labels=cmap_labels,
        cbar_label=cbar_label
    )
    for ci, xi, yi, txt in zip(c, x, y, ann_text):
        cmap_index = cmap_labels.index(ci)
        ax.annotate(txt, (xi, yi), color=cmap(cmap_index), alpha=alpha)

    return ax


def centered_temperature_profiles(temperatureProfiles):
    """Creates a figure showing the centered temperature profiles
    of the experiment

    Parameters
    ----------
    temperatureProfiles : lst, ndarray
        The temperatureprofiles during injections

    Returns
    -------
    f0 : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object

    Creates
    -------
    A figure with the mean centered temperatures
    """
    import matplotlib.pyplot as plt

    f0, ax = plt.subplots(figsize=(6.5, 6.5))
    m_temps = temperatureProfiles.transpose().mean()
    centered_temps = temperatureProfiles.apply(lambda x: x - m_temps)
    centered_temps.plot(ax=ax, colormap='jet')
    ax.set(
        title='Centered temperature profiles (injections)',
        xlabel='Axial positions (m)',
        ylabel=r'Delta Temperature (K)'
    )
    ax.legend(loc='best', fontsize=7.5, ncol=3)
    f0.tight_layout()
    plt.draw()

    return f0, ax


def mass_balance_difference(sums):
    """Creates a figure showing the mass balance differences with 100 %

    Parameters
    ----------
    sums : DataFrame
        Contains the sums (C4-, C5+ and Total)
        only Total is used

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object

    Creates
    -------
    Figure that displays the mass balance difference with 100 wt.%
    """
    import matplotlib.pyplot as plt

    f, ax = plt.subplots(1)
    sums = sums.transpose()
    sums["SUM_diff"] = sums["Total"].apply(lambda x: 100-x)
    sums["SUM_diff"].plot(kind='bar', ax=ax)
    # Validation lines
    ylines = [5, 10, 20]
    create_hor_validation_lines(ax, ylines)
    ax.set(title='Mass balance (wt.%)', ylabel='100 - sum of yields (wt.%)')
    f.tight_layout()
    plt.draw()

    return f, ax


def fid_vs_tcd(store):
    """Figure relative difference FID vs TCD

    Creates a figure showing the relative difference between
    the FID and TCD channels from the RGA.

    Parameters
    ----------
    store : HDF5 store
        opened via pandas HDFStore

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    idx = pd.IndexSlice

    f, ax = plt.subplots(1)
    fid_vs_tcd_components = [
        ind
        for
        ind in store['processed/yields/RGA_TCD'].loc[
            idx['raw_yields', :]
        ].index
        if
        (ind in store['processed/yields/RGA_FID'].loc[
            idx['raw_yields', :]
        ].index)
        ]

    average_tcd_fid = (
        store[
            'processed/yields/RGA_FID'
        ].loc[idx['raw_yields'], :].loc[fid_vs_tcd_components, :]
        + store[
            'processed/yields/RGA_TCD'
        ].loc[idx['raw_yields'], :].loc[fid_vs_tcd_components, :]) \
        / 2

    fid_vs_tcd = 100 * (
        store[
            'processed/yields/RGA_FID'
        ].loc[idx['raw_yields'], :].loc[fid_vs_tcd_components, :]
        - store[
            'processed/yields/RGA_TCD'
        ].loc[idx['raw_yields'], :].loc[fid_vs_tcd_components, :]) \
        / average_tcd_fid

    # Convert InChI to Smiles with RDKit
    fid_vs_tcd.index = [
        inchi_2_smiles(inchi_name)
        for inchi_name
        in fid_vs_tcd.index
    ]

    fid_vs_tcd.plot(kind='bar', ax=ax)
    # Validation lines
    ylines = [1, 3, 5]
    create_hor_validation_lines(ax, ylines)
    ax.legend(loc='best')
    ax.set(
        title="Relative difference between RGA's TCD and FID detectors",
        ylabel="(FID-TCD)/mean(FID,TCD) %"
    )
    f.tight_layout()
    plt.draw()

    return f, ax


def elemental_balances(store):
    """Creates a figure showing the elemental balance differences with 100 %

    Parameters
    ----------
    store : HDF5 store
        opened via pandas HDFStore

    Returns
    -------
    f3 : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """

    ######################################
    # Dilutions
    ######################################
    # HC_inflow should contain summed inflow
    # Still BUG
    ######################################
    import matplotlib.pyplot as plt
    import pandas as pd
    idx = pd.IndexSlice

    f3, ax = plt.subplots(5, sharex=True)
    elements = ['C', 'H', 'O', 'N', 'S']

    balance_sum = store['processed/yields/sums']
    inflow_balances = store['processed/yields/HC_inflow']

    for element, axx in zip(elements, ax):

        norm_element_balance = 100 * (
            inflow_balances.loc[
                idx['{}_balance'.format(element), :]
            ].sum(axis=0)
            - balance_sum.loc[
                idx['{}_balance_norm'.format(element), :]
            ].loc['Total', :]
        ) \
            / inflow_balances.loc[
                idx['{}_balance'.format(element), :]
            ].sum(axis=0)
        norm_element_balance.plot(kind='bar', ax=axx)
        ylines = [5, 7.5, 10]
        create_hor_validation_lines(axx, ylines)
        axx.set(
            title='{} balance after mass normalisation'.format(element),
            ylabel='(in - out) / in %',
            xlabel="injection tag"
        )
    f3.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.2)
    plt.draw()

    return f3, ax


def hydrogen_balance_after_carbon_norm(store):
    '''Creates a figure showing the hydrogen balance after
    carbon normalisation. And after taking into account
    the oxygen balance based on CO and CO2 yields.

    Parameters
    ----------
    store : HDF5 store
        opened through pandas HDFStore

    Returns
    -------
    f4 : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    idx = pd.IndexSlice

    ######################################
    # Dilutions
    ######################################
    # HC_inflow should contain summed inflow
    # Still BUG
    ######################################

    try:
        # Only applicable for components that are composed
        # out of solely C and H
        f4, ax = plt.subplots(1)

        H_in_HC = store[
            'processed/yields/HC_inflow'
        ].loc[idx['H_balance', :]].sum(axis=0)

        # Calculate H in from water
        # Calculate CO and CO2
        CO_inchi = '1S/CO/c1-2'
        CO2_inchi = '1S/CO2/c2-1-3'

        CO2 = store[
            'processed/yields/total'
        ].loc[idx['C_norm_yields', :]].loc[CO2_inchi, :]
        CO = store[
            'processed/yields/total'
        ].loc[idx['C_norm_yields', :]].loc[CO_inchi, :]

        H_norm_out = store[
            'processed/yields/sums'
        ].loc[idx['H_balance_norm', :]].loc['Total', :]
        H_C_norm_out = store[
            'processed/yields/sums'
        ].loc[idx['H_balance_after_C_norm', :]].loc['Total', :]

        # [C, H, O, N, S, MW]
        # O from component = Weight component times weight of O in component
        O_CO2 = CO2 * 2 * 15.999 / 44.009
        O_CO = CO * 1 * 15.999 / 28.01
        O_H2O = O_CO + O_CO2
        H_in_H2O = O_H2O * 2 * 1.008 / 15.999
        H_in_HC_H2O = H_in_HC + H_in_H2O
        Cnorm_h_balance = pd.DataFrame({
            'mass normalisation': 100 * (H_in_HC - H_norm_out)
            / H_in_HC,
            'C normalisation H2O excl': 100 * (H_in_HC - H_C_norm_out)
            / H_in_HC,
            'C normalisation H2O incl': 100 * (H_in_HC_H2O - H_C_norm_out)
            / H_in_HC_H2O
        })

        Cnorm_h_balance.plot(kind='bar', ax=ax)
        ax.set(
            title='H balance after C normalisation',
            ylabel='(in - out) / in %'
        )
        ylines = [5, 7.5, 10]
        create_hor_validation_lines(ax, ylines)
        f4.tight_layout()
        plt.draw()
    except:
        popup_error('CO or CO2 not found in dataset')

    return f4, ax


def hydrogen_vs_ethene(store, temperaturesKelvin, cmap, cmap_labels):
    """Creates a scatterplot showing hydrogen as a function of ethene.
    Ethene is taken from the TCD channel.

    Parameters
    ----------
    store : HDF5 store
    temperaturesKelvin : ndarray, list
        Experimental temperatures in Kelvin
    cmap : matplotlib cmap
        discrete colormap (matplotlib)
    cmap_labels : ndarray, lst
        List of labels used for the colormap

    Returns
    -------
    f5 : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    idx = pd.IndexSlice

    f5, ax = plt.subplots(1)

    yieldsPlotter = store[
        'processed/yields/total'
    ].loc[idx['mol_perc_wet', :]].loc[[
            r'1S/H2/h1H',  # H2
            r'1S/CO/c1-2',  # CO
            r'1S/CO2/c2-1-3',  # CO2
            r'1S/C2H4/c1-2/h1-2H2'  # C2H4
    ], :].transpose()

    create_scatterplot_with_annotations_colorbar(
        ax,
        x=yieldsPlotter[r'1S/C2H4/c1-2/h1-2H2'],
        y=yieldsPlotter[r'1S/H2/h1H'],
        c=temperaturesKelvin,
        cmap=cmap,
        cmap_labels=cmap_labels,
        cbar_label='Temperature (K)',
        ann_text=yieldsPlotter.index
    )

    ax.set(
        title='H2 vs. C2H4',
        xlabel='C2H4 normalised mol% (wet)',
        ylabel='H2 normalised mol% (wet)'
    )
    f5.tight_layout()
    plt.draw()

    return f5, ax


def hydrogen_vs_sum_ethene_propene(
    store,
    temperaturesKelvin,
    cmap,
    cmap_labels
):
    """Creates a scatterplot showing hydrogen as a function
    of the sum ethene and propene.
    Both are taken from the FID channel.

    Parameters
    ----------
    store : HDF5 store
    temperaturesKelvin : ndarray
        Experimental temperatures in Kelvin
    cmap : matplotlib cmap
        discrete colormap (matplotlib)
    cmap_labels : ndarray, lst
        List of labels used for the colormap

    Returns
    -------
    f6 : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    idx = pd.IndexSlice

    f6, ax = plt.subplots(1)

    yieldsPlotter = store[
        'processed/yields/RGA_TCD'
    ].loc[idx['mol_perc_wet', :]].loc[[
        r'1S/H2/h1H',  # H2
        r'1S/CO/c1-2',  # CO
        r'1S/CO2/c2-1-3',  # CO2
        r'1S/C2H4/c1-2/h1-2H2'  # C2H4
    ], :].transpose()
    C2H4andC3H6 = store[
        'processed/yields/total'
    ].loc[idx['mol_perc_wet', :]].loc[r'1S/C2H4/c1-2/h1-2H2', :] \
        + store[
            'processed/yields/total'
        ].loc[
            idx['mol_perc_wet', :]
        ].loc['1S/C3H6/c1-3-2/h3H,1H2,2H3', :]
    C2H4andC3H6 = C2H4andC3H6.transpose()

    create_scatterplot_with_annotations_colorbar(
        ax,
        x=C2H4andC3H6,
        y=yieldsPlotter[r'1S/H2/h1H'],
        c=temperaturesKelvin,
        cmap=cmap,
        cmap_labels=cmap_labels,
        cbar_label='Temperature (K)',
        ann_text=C2H4andC3H6.index
     )

    # Set figure parameters (title, x-, y-lbl)
    ax.set(
        title=r'H2 vs C2H4 + C3H6',
        xlabel='C2H4+C3H6 normalised mol% (wet)',
        ylabel='H2 normalised mol% (wet)'
    )
    f6.tight_layout()
    plt.draw()

    return f6, ax


def ethene_propene_sum(store):
    """Creates a figures showing the ethene, propene and C4- yield

    Parameters
    ----------
    store : HDF5 store

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    idx = pd.IndexSlice

    f, ax = plt.subplots(1)

    df_C2H4_C3H6_SUM = pd.concat([
        store[
            'processed/yields/total'
        ].loc[
            idx['C_norm_yields', :]
        ].loc[r'1S/C2H4/c1-2/h1-2H2', :],
        store[
            'processed/yields/total'
        ].loc[
            idx['C_norm_yields', :]
        ].loc[r'1S/C3H6/c1-3-2/h3H,1H2,2H3', :],
        store[
            'processed/yields/sums'
        ].loc[idx['C_norm_yields', :]].loc['C4', :],
    ], axis=1)

    # Change InChI names to Smiles
    df_C2H4_C3H6_SUM.columns = [
        inchi_2_smiles(inchi_name)
        for
        inchi_name
        in
        df_C2H4_C3H6_SUM.columns
    ]

    df_C2H4_C3H6_SUM.plot(ax=ax)
    ax.legend(loc='best')
    ax.set(
        title='C2H4, C3H6 yield and C4- fraction (wt.%)',
        xlabel='Injection tag',
        ylabel="wt.% (carbon normalised yields)"
    )
    f.tight_layout()

    plt.draw()

    return f, ax


def PCA_yields(
    df_results_transpose,
    temperaturesKelvin,
    cmap,
    cmap_labels,
    cbar_label,
    mol
):
    """Create a figure of the 1st 2 principal components of the yield data

    Parameters
    ----------
    df_results_transpose : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    temperaturesKelvin : ndarray, lst
        The setpoint temperatures for every injection
    cmap : discrete colormap
    cmap_labels : ndarray, lst
        labels that correspond with the discrete colormap
    cbar_label : string
    mol: boolean
        Distinction between mol.% wet and carbon norm wt.% dry

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : matplotlib axes
    pca_results : DataFrame
        with the injection tags as index
        and component names (InChI's) as column names
        values are PCA1 and PCA2 results

    Creates
    -------
    Figure
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    wt_mol = 'normalised mol.% wet' if mol else 'carbon normalised wt.% dry'

    f, ax = plt.subplots(1)
    X = df_results_transpose
    ax, pca_results = create_pca_figure(
        X,
        ax,
        temperaturesKelvin,
        cmap,
        cmap_labels,
        cbar_label,
#         md_dist=False,
#         loading_plot=True,
    )
    ax.set(
        title='PCA of yields ({})'.format(wt_mol)
    )
    f.tight_layout()
    plt.draw()

    pca_results = pd.DataFrame(
        data=pca_results,
        columns=['x', 'y']
    )

    return f, ax, pca_results


def PCA_temperatures(
    temperatureProfiles,
    temperaturesKelvin,
    cmap,
    cmap_labels,
    cbar_label
):
    """Create a figure of the 1st 2 principal components
    of the temperature profiles

    Parameters
    ----------
    temperatureProfiles : DataFrame
        with the injection tags as index
        and axial reactor positions as column names
        values are carbon normalised yields
    temperaturesKelvin : ndarray, lst
        The setpoint temperatures for every injection
    cmap : matplotlib cmap
    cmap_labels : lst
        List of labels used for the colormap
    cbar_label : string

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    """
    import matplotlib.pyplot as plt

    f, ax = plt.subplots(1)
    X = temperatureProfiles.transpose()
    X.columns = ['{} cm'.format(col) for col in X.columns]
    create_pca_figure(
        X,
        ax,
        temperaturesKelvin,
        cmap,
        cmap_labels,
        cbar_label,
#         md_dist=False,
#         loading_plot=True,
    )
    ax.set(title='PCA of temperature profiles')
    f.tight_layout()
    plt.draw()

    return f, ax


def create_pca_figure(
        X,
        ax,
        c,
        cmap,
        cmap_labels,
        cbar_label,
        md_dist=True,
        loading_plot=False,
    ):
    """Decompose X in its 1st two PC and make a scatterplot
    in this new axis set

    Parameters
    ----------
    X : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    c : lst
        List with different values for the colors
    cmap : matplotlib cmap
        discrete colormap
    cmap_labels : ndarray, lst
        labels that correspond with the discrete colormap
    cbar_label : String
    md_dist : boolean
        Flag to plot Mahalanobis distances or not
    loading_plot : boolean
        Flag to plot the loadings onto the score plot or not

    Returns
    -------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    pca_results : sklearn PCA instance
        The object that contains all information on the PCA results
    """
    from sklearn.decomposition import PCA
    import pandas as pd
   
    # Set opacity for scatterpoints and text 
    alpha=1

    # Create PCA instance with 2 components
    pca = PCA(n_components=2)

    # Fit the PCA instance with the data
    # Transform the input data to the new coordinates
    pca_results = pca.fit_transform(X)

    # Loading plot
    def loading_plotter(ax, loading_coeff, labels=None):
        n = loading_coeff.shape[0]
        if labels is None:
            labels = ["Var"+str(i+1) for i in range(n)]
        for i in range(n):
            ax.arrow(0, 0, loading_coeff[i,0], loading_coeff[i,1],color = 'r',alpha = 0.5)
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', facecolor='wheat')
            ax.text(
                loading_coeff[i, 0] * 1.15,
                loading_coeff[i, 1] * 1.15,
                labels[i],
                color='g',
                ha='center',
                va='center',
                bbox=props,
            )

    if loading_plot:
        import numpy as np
        from sklearn.preprocessing import MinMaxScaler
        pca_results_sc = MinMaxScaler(feature_range=(-1, 1)).fit_transform(pca_results)
        labels = [inchi_2_smiles(inchi) for inchi in X.columns]
        loading_coeff_all = np.transpose(pca.components_[0:2, :])
        # Create a DataFrame from the loadings
        df_loadings_all = pd.DataFrame(
            loading_coeff_all,
            index=labels,
            columns=['x', 'y']
        )
        # Add the magnitude of the loading vectors to the DataFrame for sorting
        # and take the 5 most import contributions
        df_loadings_all['norm'] = [
            np.linalg.norm(vector) for vector in loading_coeff_all
        ]
        df_loading = df_loadings_all.sort_values(
            by='norm', ascending=False
        ).iloc[:5]
        loading_plotter(ax, df_loading[['x', 'y']].values, labels=df_loading.index)
        alpha = 0.3
#         ax.get_legend().remove()
    else:
        pca_results_sc = pca_results

    # Create a scatterplot with annotated points
    create_scatterplot_with_annotations_colorbar(
        ax,
        x=pca_results_sc[:, 0],
        y=pca_results_sc[:, 1],
        c=c,
        cmap=cmap,
        cmap_labels=cmap_labels,
        cbar_label=cbar_label,
        ann_text=X.index,
        alpha=alpha
    )

    # Mahalanobis distances
    if md_dist:
        data = pd.DataFrame(
            {'x': pca_results[:, 0], 'y': pca_results[:, 1], 'temperature': c}
        )

        # Confidence level
        conf_lvls = [0.7, 0.8, 0.95]

        # Get confidence interval distances
        conf_int_dists = conf_int_dist_chi2(
            conf_lvls,
            2  # Number of degrees
        )
        # Create the ellips DataFrame
        ellips_data = create_ellips_data(
            data,
            conf_int_dists,
            groupby='temperature'
        )

        # Add the ellipses to the axes
        ellipses_plotter(ax, ellips_data, conf_lvls)

    ax.axhline(y=0, ls='dashed', color='grey')
    ax.axvline(x=0, ls='dashed', color='grey')
    ax.set(
        xlabel='PC 1 ({:.2%})'.format(pca.explained_variance_ratio_[0]),
        ylabel='PC 2 ({:.2%})'.format(pca.explained_variance_ratio_[1])
    )

    return ax, pca_results


def corr_scatter_matrix(store, no_corr=20, no_scatter=10):
    """
    Reorganize the data for the correlation matrix plotting and for the
    scatter matrix graph and make the graphs

    Parameters
    ----------
    store : HDF5 store
        opened via Pandas HDFStore
    no_corr : int, optional (default: 20)
        Maximal number of components to add in the correlation matrix
    no_scatter : int, optional (default: 10)
        Maximal number of components to add in the scatter matrix

    Returns
    -------
    f_corr : matplotlib.figure.Figure object
        Correlation graph figure
    ax_corr : Axes object
        ax is a single matplotlib.axes.Axes object
        Correlation graph axes
    f_scatter : matplotlib.figure.Figure object
        Scatter matrix graph figure
    sm : Axes object
        sm is a matplotlib.axes.Axes object
        Scatter matrix axes
    """
    df_results_tr_corr, df_results_tr_scatter = reorganize_corr_scatter_matrix(
        store,
        no_corr=20,
        no_scatter=10
    )

    # Plot correlation matrix
    print('Graphing correlation matrix')
    f_corr, ax_corr = plot_corr_matrix(df_results_tr_corr)

    # plot scatter matrix
    print('Graphing scatter matrix')
    f_scatter, sm = plot_scatter_matrix(df_results_tr_scatter)

    return f_corr, ax_corr, f_scatter, sm


def plot_scatter_matrix(df_results_transpose_scatter, hideticks=True):
    '''
    Create an overview figure with every scatterplot possible
    in a matrix

    Parameters
    ----------
    df_results_transpose_scatter : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    hideticks : boolean, optional (default: True)
        Switch to hide or show the tick labels in every small figure

    Returns
    -------
    f : matplotlib.figure.Figure object
    sm : Axes object
        sm is a matplotlib.axes.Axes object
    '''
    from pandas.plotting import scatter_matrix
    import matplotlib.pyplot as plt

    wt_mol = 'normalised mol.% wet'
    width = 18
    golden_ratio = 1.6
    f, sm = plt.subplots(figsize=(width, width/golden_ratio))

    sm = scatter_matrix(
        df_results_transpose_scatter,
        diagonal='kde',
        ax=sm
    )
    # Change label rotation
    [s.xaxis.label.set_rotation(45) for s in sm.reshape(-1)]
    [s.yaxis.label.set_rotation(0) for s in sm.reshape(-1)]

    # May need to offset label when rotating to prevent overlap of figure
    [s.get_yaxis().set_label_coords(-0.8, 0.5) for s in sm.reshape(-1)]

    # Hide all ticks
    if hideticks:
        [s.set_xticks(()) for s in sm.reshape(-1)]
        [s.set_yticks(()) for s in sm.reshape(-1)]

    # Give overall title
    plt.suptitle(
        'Scatter matrix with density function ({})'.format(wt_mol),
        fontsize=14
    )
    plt.draw()

    return f, sm


def plot_tsne(ax, X_embed, color_encoding=None, title=None):
    '''Creates the 2D figure of the t-SNE representation of X_embed

    Parameters
    ----------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    X_embed : ndarray
        2D t-SNE representation of initial data
        Embedding of the training data in low-dimensional space.
    color_encoding : lst, ndarray, optional
        Encoding list for the scatter point colors
    title : string, optional
        Title for the figure if necessary

    Returns
    -------
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    from sklearn.preprocessing import MinMaxScaler

    # Scale the data to make it easier for representation
    mm_sc = MinMaxScaler()
    X_embed = mm_sc.fit_transform(X_embed)

    # Plot the data
    X, Y = X_embed[:, 0], X_embed[:, 1]
    if color_encoding is not None:
        ax.scatter(x=X, y=Y, c=color_encoding)
    else:
        ax.scatter(x=X, y=Y)

    # Remove x-, and y-ticks
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    # Add a figure title
    if title is not None:
        ax.set(title=title,)

    return ax


def plot_tsne_annot_dbscan(
    df_results_transpose,
    cmap,
    cmap_labels,
    annot_colors,
    ann_texts,
    title=None
):
    '''Creates the 2D figure of the t-SNE representation of the data and
    annotates the datapoints based on the ann_texts, annotations are colored
    based on the list of annot_colors added as a parameter. The color of the
    datapoints is based the clustering algorithm DBSCAN.

    Parameters
    ----------
    df_results_transpose : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    cmap : matplotlib cmap
    cmap_labels : lst
        List of labels used for the colormap
    annot_colors : lst, ndarray
        List of the datapoints with a color encoding value
    ann_texts : lst, ndarray
        List of annotation text for plotting
    title : string, optional
        Title for the figure if necessary

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    import matplotlib.pyplot as plt
    f, ax = plt.subplots()

    X_embed = tsne_2D(df_results_transpose)
    db = dbscan_color_encoding(X_embed)
    plot_tsne(ax, X_embed, color_encoding=db, title=title)
    annotate_tsne_plot(ax, X_embed, cmap, cmap_labels, annot_colors, ann_texts)

    f.tight_layout()
    plt.draw()

    return f, ax


def plot_corr_matrix(dataframe):
    ''' Creates a graph of the correlation matrix of a given dataframe

    Parameters
    ----------
    dataframe : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields

    Returns
    -------
    f : matplotlib.figure.Figure object
    ax : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    # Create correlation matrix
    df_corr = dataframe.corr()

    # Save as excel file
    df_corr.to_excel('correlation_matrix.xlsx')

    # Generate a mask for the upper triangle
    mask = np.zeros_like(df_corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    f, ax = plt.subplots(figsize=(7, 7))
    sns.heatmap(
        df_corr,
        mask=mask,
        cmap=cmap,
        center=0,
        square=True,
        linewidths=.5,
        cbar_kws={"shrink": .5, 'label': 'Correlation coefficient'},
        annot=True,
        fmt=".2f",
        annot_kws={"size": 5.5},
        ax=ax,
    )
    # Change orientation of x- and y-ticks
    ax.tick_params(axis='x', labelrotation=90,)
    ax.tick_params(axis='y', labelrotation=0,)

    # Create the title
    wt_mol = 'normalised mol.% wet'
    ax.set(
        title='Correlation matrix for yields ({})'.format(wt_mol)
    )

    # Restore default Matplotlib style
    plt.style.use('default')

    try:
        f.tight_layout()
    except:
        pass

    return f, ax


def subplots_compounds_vs_temperature(
    pivotted_avg,
    pivotted_std,
    title,
    mol
):
    '''Creates subplots of different compounds versus
    temperature with errorbars

    Parameters
    ----------
    pivotted_avg : DataFrame
    pivotted_std : DataFrame
    title : string
    mol : boolean
        Distinction between mol.% wet and carbon norm wt.% dry

    Returns
    -------
    f_major : matplotlib.figure.Figure object
    axes : Axes object
        ax is a single matplotlib.axes.Axes object
    '''
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    import numpy as np

    number = int(np.ceil(np.sqrt(len(pivotted_avg.columns))))
    f_major, axes = plt.subplots(
        nrows=number,
        ncols=number,
        figsize=(16, 10),
        sharex='col'
    )
    for inchiName, ax in zip(pivotted_avg.columns, np.ravel(axes)):
        exp_name = 'BSSC Experiment_exp'
        x = pivotted_avg[inchiName].unstack('exp_sim')[exp_name].index.tolist()
        y = pivotted_avg[inchiName].unstack('exp_sim')[exp_name].tolist()
        errbars = pivotted_std[inchiName].tolist()

        smileName = inchi_2_smiles(inchiName)

        plot_errBar_exp(ax, x, y, errbars, 'exp')

        wt_mol = 'mol.% wet' if mol else 'wt.% dry'
        subtitle = '{} {}'.format(smileName, wt_mol)
        ax.text(
            .5,
            .85,
            subtitle,
            horizontalalignment='center',
            transform=ax.transAxes
        )

        if number > 5:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    # Set window title for the overall figure
    f_major.canvas.set_window_title(title)
    # Adjust label
    for ax in axes[-1, :]:
        ax.set(xlabel='Temperature (K)')
        [tick.set_rotation(45) for tick in ax.get_xticklabels()]

    f_major.tight_layout()
    plt.draw()

    return f_major, axes


def tsne_2D(X):
    '''Creates the 2D t-distributed Stochastic Neighbor Embedding
    representation of the input X data that contains multiple dimensions.

    Parameters
    ----------
    X : ndarray, DataFrame
        Multidimensional input data.
        If the metric is precomputed X must be a square distance matrix.
        Otherwise it contains a sample per row. If the method is exact, X may
        be a sparse matrix of type csr, csc, or coo.

    Returns
    -------
    X_embedded : ndarray
        Embedding of the training data (X) in low-dimensional space.
    '''
    from sklearn.manifold import TSNE

    X_embedded = TSNE(
        n_components=2,
        # init='pca',
        perplexity=4
    ).fit_transform(X)

    return X_embedded


def reorganize_corr_scatter_matrix(
    store,
    no_corr=20,
    no_scatter=10
):
    """
    Reorganise the data for the correlation matrix plotting and for the
    scatter matrix graph

    Parameters
    ----------
    store : HDF5 store
        opened via Pandas HDFStore
    no_corr : int, optional (default: 20)
        Maximal number of components to add in the correlation matrix
    no_scatter : int, optional (default: 10)
        Maximal number of components to add in the scatter matrix

    Returns
    -------
    df_results_tr_corr : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    df_results_tr_scatter : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    """
    import pandas as pd
    idx = pd.IndexSlice  # MultiIndex slicing

    df_results = store['processed/yields/total'].loc[idx['mol_perc_wet', :]]
    rs_sum_sort = df_results.sum(axis=1).sort_values(ascending=False)

    # Take only the 20 (corr), and 10 (scatter) most abundant components
    mask_ind_corr = rs_sum_sort.head(no_corr).index
    mask_ind_scatter = rs_sum_sort.head(no_scatter).index
    df_results_corr_tr = df_results.loc[mask_ind_corr].transpose()
    df_results_scatter_tr = df_results.loc[mask_ind_scatter].transpose()

    # Sort components based on their carbon number and convert InChI to SMILES
    df_results_tr_corr = convert_smiles_sort(df_results_corr_tr)
    df_results_tr_scatter = convert_smiles_sort(df_results_scatter_tr)

    return df_results_tr_corr, df_results_tr_scatter


def graph_experiment(
        store,
        temperatureProfiles,
        sums,
        temperaturesKelvin,
        cmap,
        cmap_labels,
        df_results_transpose,
        cbar_label,
        sim_exp,
        mol
):
    """
    Module plotting all validation graphs for an experiment

    Parameters
    ----------
    store : HDF5 store
        opened via pandas HDFStore
    temperatureProfiles : DataFrame
        with the injection tags as index
        and axial reactor positions as column names
        values are carbon normalised yields
    sums : DataFrame
        Contains the sums (C4-, C5+ and Total)
        columns : injection tags
        index : [C4-, C5+, Total]
    temperaturesKelvin : ndarray, lst
        The setpoint temperatures for every injection
    cmap : matplotlib cmap
    cmap_labels : lst
        List of labels used for the colormap
    df_results_transpose : DataFrame
        with the injection tags as index
        and component names (Smiles) as column names
        values are carbon normalised yields
    cbar_label : string
    sim_exp : DataFrame
    mol : boolean
        Distinction between mol.% wet and carbon norm wt.% dry

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    # Create folder to save the figures
    fig_loc = r"Figures_experiment"
    if not os.path.exists(fig_loc):
        os.makedirs(fig_loc)

    # Create figure to show the differences with 100 %
    # of the raw GC results data
    print('Graphing how closed the overall mass balance is.')
    f_mass, _ = mass_balance_difference(sums)
    f_save(
        f_mass,
        r"{}/mass_balance.tif".format(fig_loc)
    )
 
    # Create figure to show elemental balance differences from normalised data
    print('Graphing elemental balances')
    f_el, _ = elemental_balances(store)
    f_save(
        f_el,
        r"{}/elemental_balances.tif".format(fig_loc)
    )
 
    # Create H balance from carbon normalised data
    print('Graphing H balance from carbon normalised data')
    f_h, _ = hydrogen_balance_after_carbon_norm(store)
    f_save(
        f_h,
        r"{}/hydrogen_bal_after_c_norm.tif".format(fig_loc)
    )
 
    # Create figure to show differences between RGA and FID data
    print('Graphing relative difference between RGA and FID data.')
    f_fid_tcd, _ = fid_vs_tcd(store)
    f_save(
        f_fid_tcd,
        r"{}/fid_vs_tcd.tif".format(fig_loc)
    )
 
    # Ethylene, propylene and sum of C4-
    print('Graphing ethylene, propylene and the sum of the C4- fraction')
    f_eth_prop_sum, _ = ethene_propene_sum(store)
    f_save(
        f_eth_prop_sum,
        r"{}/ethene_propene_sum.tif".format(fig_loc)
    )
 
    # Scatter plot displaying hydrogen as a function of ethene
    print('Graphing hydrogen as a function of ethene')
    f_h2_c2h4, _ = hydrogen_vs_ethene(
        store,
        temperaturesKelvin,
        cmap,
        cmap_labels
    )
    f_save(
        f_h2_c2h4,
        r"{}/H2_C2H4.tif".format(fig_loc)
    )
 
    # Hydrogen versus sum of ethylene and propylene
    print('Graphing hydrogen as a function of the sum of ethene and propene')
    f_h2_sum, _ = hydrogen_vs_sum_ethene_propene(
        store,
        temperaturesKelvin,
        cmap,
        cmap_labels
    )
    f_save(
        f_h2_sum,
        r"{}/H2_sum.tif".format(fig_loc)
    )
 
    # t-SNE of yield_data
    print('Graphing t-SNE embedding of results')
    f_tsne, _ = plot_tsne_annot_dbscan(
        df_results_transpose,
        cmap,
        cmap_labels,
        temperaturesKelvin,  # color_encoding
        df_results_transpose.index,  # ann_texts
        title='t-SNE visualisation of the carbon normalised yield data'
    )
    f_save(
        f_tsne,
        r"{}/t-SNE.tif".format(fig_loc)
    )

    # PCA of results to detect outliers
    print('Graphing principal component analysis of results')
    f_pca_yields, _, _ = PCA_yields(
        df_results_transpose,
        temperaturesKelvin,
        cmap,
        cmap_labels,
        cbar_label,
        mol
    )
    f_save(
        f_pca_yields,
        r"{}/pca_yields.tif".format(fig_loc)
    )

    # PCA of temperature profiles
    print('Graphing principal component analysis of temperature profiles')
    f_pca_temperatures, _ = PCA_temperatures(
        temperatureProfiles,
        temperaturesKelvin,
        cmap,
        cmap_labels,
        cbar_label
    )
    f_save(
        f_pca_temperatures,
        r"{}/pca_temperatures.tif".format(fig_loc)
    )

    # Create temperature profiles during every injection
    print("Graphing temperature profiles.")
    f_c_T, _ = centered_temperature_profiles(temperatureProfiles)
    f_save(
        f_c_T,
        r"{}/centered_temperatures.tif".format(fig_loc)
    )
 
    # Correlation matrix and scatter matrix
    print('Rearranging data for scatter and correlation matrix')
    f_corr, _, f_scatter, _ = corr_scatter_matrix(
        store,
        no_corr=20,
        no_scatter=10
    )
    f_save(
        f_corr,
        r"{}/correlation_matrix.tif".format(fig_loc)
    )
    f_save(
        f_scatter,
        r"{}/scatter_matrix.tif".format(fig_loc)
    )
 
    # Yields vs Temperature
    print('Graphing major yields versus Temperature')
    pivotted_avg = sim_exp.pivot_table(
        columns=['exp_sim', 'Temperature'],
        aggfunc=np.mean
    ).unstack('InChI')
    pivotted_std = sim_exp.pivot_table(
        columns=['exp_sim', 'Temperature'],
        aggfunc=np.std
    ).unstack('InChI')
    mask_major = pivotted_avg.max() > 1.0
    mask_minor = -mask_major
    pivotted_avg_major = pivotted_avg.loc[:, mask_major]
    pivotted_std_major = pivotted_std.loc[:, mask_major]
    pivotted_avg_minor = pivotted_avg.loc[:, mask_minor]
    pivotted_std_minor = pivotted_std.loc[:, mask_minor]
    pivotted_avg_major.loc[:, 'Sum'] = pivotted_avg_major.sum(axis=1)
    pivotted_std_major.loc[:, 'Sum'] = 0
    pivotted_avg_minor.loc[:, 'Sum'] = pivotted_avg_minor.sum(axis=1)
    pivotted_std_minor.loc[:, 'Sum'] = 0
 
    # Plot major products (> 1 mol.% wet/wt. dry)
    f_major, _ = subplots_compounds_vs_temperature(
        pivotted_avg_major,
        pivotted_std_major,
        'Major products vs Temperature',
        mol
    )
    f_save(
        f_major,
        r"{}/major_components.tif".format(fig_loc)
    )
#
#     # Plot minor products (< 1 mol.% wet/wt. dry)
#     print('Graphing minor yields versus Temperature')
#     subplots_compounds_vs_temperature(
#         pivotted_avg_minor,
#         pivotted_std_minor,
#         'Minor products vs Temperature',
#         mol
#     )

    plt.show()

    return None
