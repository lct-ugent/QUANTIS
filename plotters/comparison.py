'''
Created on 7 Oct 2016

@author: shsymoen
'''

from plotters import generalPlotters
from sklearn.metrics import r2_score


def create_hor_validation_lines(ax, ylines):

    sorted(ylines)
    colors = ['green', 'orange', 'red']
    lowerlimit, upperlimit = ax.get_xlim()
    pos = xlim_to_01(lowerlimit, upperlimit, 95)
    ax.axhline(y=0)

    for yline, color in zip(ylines, colors):
        ax.axhline(y=yline, color=color, ls='dashed')
        ax.axhline(y=-yline, color=color, ls='dashed')
        ax.text(
            pos,
            yline,
            '{}%'.format(yline),
            va='bottom',
            ha='right',
            color=color,
            bbox=dict(facecolor="white", edgecolor='none', pad=1)
        )
        ax.text(
            pos,
            -yline,
            '{}%'.format(-yline),
            va='bottom',
            ha='right',
            color=color,
            bbox=dict(facecolor="white", edgecolor='none', pad=1)
        )


def xlim_to_01(lowerlimit, upperlimit, percentage):
    ''' Calculates the position for a text element '''
    return lowerlimit + percentage / 100 * (upperlimit - lowerlimit)


def xpos_to_ypos(xpos, xp, fp):
    ''' Calculates the linear interpolated value '''
    import numpy as np

    return np.interp(xpos, xp, fp)


def add_intervals_parity_plot(ax, highest, lowest):
    ''' Puts interval lines in a parity scatter plot '''
    ax.set_xlim((lowest, highest))
    ax.set_ylim((lowest, highest))
    bbox_props = dict(facecolor="white", edgecolor='none', pad=1, alpha=0.5)
    # Get ax limits
    xlow = ax.get_xlim()[0]
    xhigh = ax.get_xlim()[1]
    ylow = ax.get_ylim()[0]
    yhigh = ax.get_ylim()[1]

    # Put middle line
    ax.plot((0, highest), (0, highest), c=".3", alpha=0.5)

    # 5% ranges
    percentage = 95
    xposfix = xlim_to_01(xlow, xhigh, percentage)
    yposfix = xlim_to_01(ylow, yhigh, percentage)
    xp = [xlow, xhigh]
    yp = [ylow, yhigh]
    ax.plot(
        (0, highest*percentage/100),
        (0, highest),
        ls="dotted",
        c=".2",
        alpha=0.5
    )
    ax.plot(
        (0, highest),
        (0, highest*percentage/100),
        ls="dotted",
        c=".2",
        alpha=0.5
    )
    ax.text(
        xposfix,
        xpos_to_ypos(xposfix, xp, [ylow, highest*percentage/100]),
        '-5%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )
    ax.text(
        xpos_to_ypos(yposfix, yp, [xlow, highest*percentage/100]),
        yposfix,
        '+5%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )

    # 10% ranges
    percentage = 90
    # xpos = xlim_to_01(xlow, xhigh, percentage)
    # ypos = xlim_to_01(ylow, yhigh, percentage)
    ax.plot(
        (0, highest*percentage/100),
        (0, highest),
        ls="-.",
        c=".4",
        alpha=0.5
    )
    ax.plot(
        (0, highest),
        (0, highest*percentage/100),
        ls="-.",
        c=".4",
        alpha=0.5
    )
    ax.text(
        xposfix,
        xpos_to_ypos(xposfix, xp, [ylow, highest*percentage/100]),
        '-10%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )
    ax.text(
        xpos_to_ypos(yposfix*0.95, yp, [xlow, highest*percentage/100]),
        yposfix*0.95,
        '+10%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )

    # 20% ranges
    percentage = 80
    # xpos = xlim_to_01(xlow, xhigh, percentage)
    # ypos = xlim_to_01(ylow, yhigh, percentage)
    ax.plot(
        (0, highest*percentage/100),
        (0, highest),
        ls="--",
        c=".6",
        alpha=0.5
    )
    ax.plot(
        (0, highest),
        (0, highest*percentage/100),
        ls="--",
        c=".6",
        alpha=0.5
    )
    ax.text(
        xposfix,
        xpos_to_ypos(xposfix, xp, [ylow, highest*percentage/100]),
        '-20%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )
    ax.text(
        xpos_to_ypos(yposfix, yp, [xlow, highest*percentage/100]),
        yposfix,
        '+20%',
        alpha=0.5,
        va='center',
        ha='center',
        bbox=bbox_props
    )


def create_parity_plot(
    ax,
    experiment,
    simulation_names,
    compounds,
    grouped_sim_exp,
    compound_Smiles,
    compound_iupac,
    wt_mol
):
    '''Creates a parity plot for a specific
    component for different simulations
    '''

    no = len(simulation_names)
    colors, markers = generalPlotters.colors_markers(no)
    single = False

    # Check if multiple compounds or single
    if not isinstance(compounds, list):
        compounds = [compounds]

    # Check if single simulation or multiple
    if not isinstance(simulation_names, list):
        simulation_names = [simulation_names]
        single = True

    for compound in compounds:
        for color, marker, simulation_name in zip(
            colors,
            markers,
            simulation_names
        ):
            fc = color if single else 'white'
            ax.scatter(x=grouped_sim_exp.get_group(experiment)[compound],
                       y=grouped_sim_exp.get_group(simulation_name)[compound],
                       marker=marker,
                       s=20,
                       facecolor=fc,
                       edgecolor=color,
                       label=simulation_name.replace('_sim', '')
                       )

    # Get limits of plot
    if ax.get_xlim()[0] < ax.get_ylim()[0]:
        lowest = ax.get_xlim()[0]
    else:
        lowest = ax.get_ylim()[0]
    if lowest < 0:
        lowest = 0
    if ax.get_xlim()[1] > ax.get_ylim()[1]:
        highest = ax.get_xlim()[1]
    else:
        highest = ax.get_ylim()[1]

    if (highest > 100) & (compound != 'Temperature'):
        highest = 100

    # Add interval lines
    add_intervals_parity_plot(ax, highest, lowest)

    # Set title and labels
    if not single:
        ylabel = 'Simulation {}'.format(wt_mol)
    else:
        ylabel = '{} {}'.format(simulation_names[0], wt_mol)
    ax.set(title='{} {}'.format(compound_Smiles, wt_mol),
           xlabel='{} {}'.format(experiment.replace('_exp', ''), wt_mol),
           ylabel=ylabel
           )
    if not single:
        handles, labels = ax.get_legend_handles_labels()
        number = len(simulation_names)
        ax.legend(handles=handles[:number],
                  labels=labels[:number],
                  loc='best',
                  scatterpoints=1,
                  fontsize=10,
                  frameon=False
                  )


def create_plot_component_vs_conversion(
    ax,
    exps,
    sims,
    compound,
    exps_std,
    name_experiment,
    names_simulations,
    inchi2smiles_total,
    inchi2iupac,
    wt_mol
):
    ''' Creates a plot for a component vs the conversion '''
    import numpy as np

    no = len(names_simulations)
    colors, _ = generalPlotters.colors_markers(no)

    y_to_plot = exps[compound].tolist()
    errbars_to_plot = exps_std[compound].tolist()
    xerrbars = exps_std['Conversion'].tolist()
    conv = exps['Conversion'].tolist()

    if np.isnan(y_to_plot).sum() < len(y_to_plot):
        ax.errorbar(
            conv,
            y_to_plot,
            fmt='o',
            yerr=errbars_to_plot,
            xerr=xerrbars,
            label=name_experiment[0].replace('_exp', '')
        )

        # iterate over all simulations
        for simulation, color in zip(names_simulations, colors):
            y_sim = sims.loc[simulation].set_index('Conversion').sort_index()
            # Don't show simulations without values
            if y_sim[compound].isnull().sum() == 0:
                r2 = r2_score(y_to_plot, y_sim[compound])
                lbl_legend = '{} - R^2: {:.2%}'.format(
                    simulation.replace('_sim', ''), r2
                )
                sims.loc[simulation].plot(
                    x='Conversion',
                    y=compound,
                    ax=ax,
                    label=lbl_legend,
                    color=color,
                )

        # Changes for the figure
        ax.autoscale(tight=False)
        compound_Smiles = inchi2smiles_total.loc[compound]
        try:
            compound_iupac = inchi2iupac[compound]
        except:
            compound_iupac = compound_Smiles
        ax.set(xlabel='Conversion {}'.format(wt_mol),
               ylabel='{} {}'.format(compound_iupac, wt_mol),
               title=compound_Smiles)
        ax.legend(loc='best', framealpha=0.5, frameon=False, fontsize=10)
        ax.margins(0.1)


def plot_errBar_exp(ax, x, y, errbars, label=None):
    ''' Creates an errorbar plot '''
    ax.errorbar(x, y, fmt='o', yerr=errbars, label=label)


def create_comparison_figures(
    sim_exp,
    inchi2iupac,
    feedName_InChI,
    mol
):
    ''' Creates all comparison figures '''
    from rdkit import Chem
    import matplotlib.pyplot as plt
    import os
    import pandas as pd
    import numpy as np

    wt_mol = 'mol.% wet' if mol else 'wt.% dry'

    # Take only components that are available
    # in experiments and all simulations
    sim_exp_reduced = sim_exp.dropna(axis=1)

    # Group results based on simulation/experiment label
    grouped_sim_exp = sim_exp_reduced.groupby('exp_sim')

    # Use pivot on the comparison file with aggregation function mean and stdev
    ######################################
    # Enhancement
    ######################################
    # Pivot really necessay?
    ######################################
    pivotted_avg = sim_exp_reduced.pivot_table(
        columns=['exp_sim', 'Temperature'],
        aggfunc=np.mean
    ).unstack('InChI')
    pivotted_std = sim_exp_reduced.pivot_table(
        columns=['exp_sim', 'Temperature'],
        aggfunc=np.std
    ).unstack('InChI')

    # Get InChI's to Smiles converter from rdkit
    compounds = pivotted_avg.columns.dropna()
    try:
        compounds.remove('Conversion')
    except:
        pass

    smiles = []
    compound_list = []
    for compound in compounds:
        if compound != 'No_INCHI':
            try:
                m = Chem.MolFromInchi("InChI="+compound.replace(' ', ''))
                smiles.append(Chem.MolToSmiles(m))
                compound_list.append(compound)
            except:
                pass

    inchi2smiles_total = pd.Series(data=smiles, index=compound_list)

    # Number needed for subplots (amount of rows and columns)
    number = int(np.ceil(np.sqrt(len(sim_exp_reduced.columns)-3.0)))

    # Conversion vs Temperature
    f1, ax = plt.subplots()

    names_simulations = [sim for sim in grouped_sim_exp.groups if 'sim' in sim]
    names_simulations.sort()
    name_experiment = [exp for exp in grouped_sim_exp.groups if 'exp' in exp]

    if 'Conversion' in pivotted_avg.columns:
        for experiment in name_experiment:
            # Collect errorbars on standard deviations
            errbars_to_plot = pivotted_std.loc[
                :, 'Conversion'
            ].unstack('exp_sim')[experiment].tolist()
            # Collect averages
            exp_to_ploty = pivotted_avg.loc[
                :, 'Conversion'
            ].unstack('exp_sim')[experiment].tolist()
            exp_to_plotx = pivotted_avg.loc[
                :, 'Conversion'
            ].unstack('exp_sim')[experiment].index.tolist()
            lbl = experiment.replace('_exp', '')

            plot_errBar_exp(
                ax,
                exp_to_plotx,
                exp_to_ploty,
                errbars_to_plot,
                label=lbl,
            )

        no = len(names_simulations)
        colors, _ = generalPlotters.colors_markers(no)
     
        for simulation, color in zip(names_simulations, colors):
            sim_y = pivotted_avg.loc[
                :, 'Conversion'
            ].unstack('exp_sim')[simulation]
            r2 = r2_score(exp_to_ploty, sim_y)
            lbl_legend = '{} - R^2: {:.2%}'.format(
                simulation.replace('_sim', ''), r2
            )
            sim_y.plot(ax=ax, label=lbl_legend, color=color)

        # Adapt figure parameters
        ax.autoscale(tight=False)
        ax.set(
            xlabel='Temperature (K)',
            ylabel='Conversion {}'.format(wt_mol),
            title='Conversion\n{} {}\nvs. Temperature'.format(
                inchi2smiles_total.loc[feedName_InChI[0]],
                wt_mol
                )
               )
        ax.legend(loc='best', framealpha=0.5, frameon=False, fontsize=10)
        ax.margins(0.1)
        generalPlotters.setFigLinesBW(f1)
        f1.tight_layout()
        print("\nSaving figure conversion versus temperature.")
        f1.savefig(
            r"Comparison\Figures\Conversion_vs_temperature.tif",
            format='tif',
            dpi=200
        )
        print("Figure conversion versus temperature saved.")

        plt.close()

        # Yield vs conversion
        if not os.path.exists(r"Comparison\Figures\Component_vs_conversion"):
            os.makedirs(r'Comparison\Figures\Component_vs_conversion')

        f3 = plt.figure(figsize=(20, 20))

        exps = pivotted_avg.loc[name_experiment]
        exps_std = pivotted_std.loc[name_experiment]
        sims = pivotted_avg.loc[names_simulations]

        print("\nCreating plots components vs conversion...")
        for i, compound in enumerate(exps.columns):
            if compound != 'Conversion':  # and i < 25:
                try:
                    compound_iupac = inchi2iupac[compound]
                except:
                    compound_iupac = inchi2smiles_total.loc[
                        compound
                    ].replace('/', 'slash').replace("\\", "backslash")
                # All figures as subplots in one big figure
                ax = f3.add_subplot(number, number, i+1)
                create_plot_component_vs_conversion(
                    ax,
                    exps,
                    sims,
                    compound,
                    exps_std,
                    name_experiment,
                    names_simulations,
                    inchi2smiles_total,
                    inchi2iupac,
                    wt_mol
                    )

                # All figures stored separately
                f5, ax3 = plt.subplots()
                create_plot_component_vs_conversion(
                    ax3,
                    exps,
                    sims,
                    compound,
                    exps_std,
                    name_experiment,
                    names_simulations,
                    inchi2smiles_total,
                    inchi2iupac,
                    wt_mol
                    )
                generalPlotters.setFigLinesBW(f5)
                fig_name = os.path.join(
                    r"Comparison\Figures\Component_vs_conversion",
                    r"{}.tif".format(compound_iupac),
                )
                f5.savefig(
                    fig_name,
                    format='tif',
                    dpi=200
                )
                plt.close(f5)

        plt.tight_layout(h_pad=3.0)
        generalPlotters.setFigLinesBW(f3)
        print("\nSaving figure 'components versus conversion'.")
        fig_name = os.path.join(
            r"Comparison\Figures\Component_vs_conversion",
            r"components_vs_conversion.tif"
        )
        f3.savefig(
            fig_name,
            format='tif',
            dpi=500
        )
        print("Figure components versus conversion saved.")

        plt.close()

    else:
        print(
            "No conversion was calculated so the conversion vs temperature "
            "figure is not created."
        )

    # Parity plots
    experiment = name_experiment[0]
    compounds = grouped_sim_exp.get_group(experiment).columns

    number = int(np.ceil(np.sqrt(len(compounds))))

    # Parity plots all simulations together
    if not os.path.exists(r"Comparison\Figures\Parity_plots"):
        os.makedirs(r'Comparison\Figures\Parity_plots')

    if len(names_simulations) > 1:
        f2 = plt.figure(figsize=(20, 20))
        if not os.path.exists(
            r"Comparison\Figures\Parity_plots\all_simulations"
        ):
            os.makedirs(r"Comparison\Figures\Parity_plots\all_simulations")

        for i, compound in enumerate(compounds):
            if (compound != 'exp_sim') and (compound != 'Temperature'):
                if (compound != 'Conversion'):
                    compound_Smiles = inchi2smiles_total.loc[compound]
                    try:
                        compound_iupac = inchi2iupac[compound]
                    except:
                        compound_iupac = compound_Smiles.replace(
                            '/', 'slash'
                        ).replace("\\", "backslash")
                else:
                    compound_iupac = compound
                    compound_Smiles = compound

                # Create big overall figure with all parity plots
                ax = f2.add_subplot(number, number, i+1)
                create_parity_plot(
                    ax,
                    experiment,
                    names_simulations,
                    compound,
                    grouped_sim_exp,
                    compound_Smiles,
                    compound_iupac,
                    wt_mol
                )

                # Create and save separate figures
                f4, ax2 = plt.subplots()
                create_parity_plot(
                    ax2,
                    experiment,
                    names_simulations,
                    compound,
                    grouped_sim_exp,
                    compound_Smiles,
                    compound_iupac,
                    wt_mol
                )
                fig_name = os.path.join(
                    r"Comparison\Figures\Parity_plots\all_simulations",
                    r"{}_parity_plot.tif".format(compound_iupac)

                )
                f4.savefig(
                    fig_name,
                    format='tif',
                    dpi=200
                )

                # Close the created figure
                plt.close(f4)

        f2.tight_layout(h_pad=3.0)
        print("Saving figure parity plots.")
        fig_name = os.path.join(
            r"Comparison\Figures\Parity_plots\all_simulations",
            r"parity_plots.tif"
        )
        f2.savefig(
            fig_name,
            format='tif',
            dpi=500
        )

        # Close the subplots figure
        plt.close(f2)

    # Parity plots all simulations separate
    # Iterate over all simulations
    for simulation_name in names_simulations:
        print("\nProcessing parity plots of {}".format(simulation_name))
        if not os.path.exists(
            r"Comparison\Figures\Parity_plots\{}".format(
                simulation_name.replace("CHEMKIN_", "").replace("_sim", "")
            )
        ):
            os.makedirs(
                r"Comparison\Figures\Parity_plots\{}".format(
                    simulation_name.replace("CHEMKIN_", "").replace("_sim", "")
                )
            )

        f2 = plt.figure(figsize=(20, 20))

        # Iterate over all compounds to create separate parity plots
        for i, compound in enumerate(compounds):
            # Get IUPAC names from the InChI name of the component
            if (compound != 'exp_sim') and (compound != 'Temperature'):
                if (compound != 'Conversion'):
                    compound_Smiles = inchi2smiles_total.loc[compound]
                    try:
                        compound_iupac = inchi2iupac[compound]
                    except:
                        compound_iupac = compound_Smiles.replace(
                            '/', 'slash'
                        ).replace("\\", "backslash")
                else:
                    compound_iupac = compound
                    compound_Smiles = compound

                # Create big overall figure with all parity plots
                ax = f2.add_subplot(number, number, i+1)
                create_parity_plot(
                    ax,
                    experiment,
                    simulation_name,
                    compound,
                    grouped_sim_exp,
                    compound_Smiles,
                    compound_iupac,
                    wt_mol
                )

                # Create and save separate figures
                f4, ax2 = plt.subplots()
                create_parity_plot(
                    ax2,
                    experiment,
                    simulation_name,
                    compound,
                    grouped_sim_exp,
                    compound_Smiles,
                    compound_iupac,
                    wt_mol
                )
                fig_name = os.path.join(
                    r"Comparison\Figures\Parity_plots",
                    "{}\{}_parity_plot.tif".format(
                        simulation_name.replace("CHEMKIN_", "").replace(
                            "_sim", ""),
                        compound_iupac
                    )
                )
                f4.savefig(
                    fig_name,
                    format='tif',
                    dpi=200
                )

                plt.close(f4)

        f2.tight_layout(h_pad=3.0)
        print("Saving figure parity plots.")
        f2.savefig(
            r'Comparison\Figures\Parity_plots\{}\parity_plots.tif'.format(
                simulation_name.replace("CHEMKIN_", "").replace("_sim", "")),
            format='tif',
            dpi=500
        )
        print("Figure parity plots saved.")

    # Try all compounds

    f, ax3 = plt.subplots()
    compounds = [
        comp
        for comp
        in compounds if (comp != 'exp_sim' and comp != 'Temperature')
    ]
    create_parity_plot(
        ax3,
        experiment,
        names_simulations,
        compounds,
        grouped_sim_exp,
        "All compounds",
        compound_iupac,
        wt_mol
    )
    f.tight_layout()
    try:
        fig_name = os.path.join(
            r"Comparison\Figures\Parity_plots\all_simulations",
            "all_compound.tif"
        )

        f.savefig(
            fig_name,
            format='tif',
            dpi=200
        )
    except:
        os.mkdir(r'Comparison\Figures\Parity_plots\all_simulations')
        fig_name = os.path.join(
            r"Comparison\Figures\Parity_plots\all_simulations",
            "all_compound.tif"
        )
        f.savefig(
            fig_name,
            format='tif',
            dpi=200
        )

    plt.close()
