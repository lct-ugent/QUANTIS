'''
Created on 13 Oct 2016

@author: shsymoen
'''

import os
import pandas as pd
from dataProcessing import utils


def read_coilsim1D_sensitivityOverview(
    pathName,
    coilsimConverter,
    mol,
    store,
    v39=0
):
    '''Reads the COILSIM1D sensitivityOverview and
    adds it to the HDF5 store
    '''

    # Save experimental names
    col_indices = store['processing/conditions'].index

    wt_mass = 'mol' if mol else 'mass'
    mol_dry = 'dry_wt' if mol == 0 else 'wet_mol'

    if not v39:
        sens_Overview = 'SensitivityOverview_{}.csv'.format(wt_mass)
        hdfkey = 'processed/simulations/COILSIM1D_{}'.format(mol_dry)
    else:
        sens_Overview = 'SensitivityOverview_v39_{}.csv'.format(wt_mass)
        hdfkey = 'processed/simulations/COILSIM1D_v39_{}'.format(mol_dry)

    coilsimSummary = os.path.join(pathName, sens_Overview)

    if os.path.isfile(coilsimSummary):

        # Read COILSIM1D's summary file 'sensitivityOverview'
        df_sim_coilsim = pd.read_csv(
            coilsimSummary,
            skipfooter=24,
            engine='python',
            index_col='Component name'
        )
        # Clean up the dataframe
        df_sim_coilsim.columns = [
            col.strip() for col in df_sim_coilsim.columns
        ]
        df_sim_coilsim.index = [ind.strip() for ind in df_sim_coilsim.index]
        if 'Component Number' in df_sim_coilsim.columns:
            df_sim_coilsim.drop('Component Number', axis=1, inplace=True)
        df_sim_coilsim.dropna(axis=1, how='all', inplace=True)

        # Change index of COILSIM1D converter file
        coilsimConverterInChI = coilsimConverter.reset_index()
        coilsimConverterInChI.set_index(keys='name COILSIM', inplace=True)

        # Append InChI names to COILSIM1D DataFrame as index
        df_sim_coilsim.index = coilsimConverterInChI.loc[
            df_sim_coilsim.index, 'InChI'
        ].tolist()
        df_sim_coilsim.index.name = 'InChI'

        # Store COILSIM1D simulation results to HDF5 file
        df_sim_coilsim.columns = col_indices

        # Store simulation results to HDF5 file
        utils.storing_2_HDF5_store(
            df_sim_coilsim.reset_index(),
            store,
            hdfkey,
            data_column=False
        )

        print(
            "COILSIM1D sensitivityOverview successfully added "
            "to the HFD5 file container."
        )
        
    else:
        print(
            "COILSIM1D {} couldn't be found in {}".format(
                sens_Overview,
                pathName
            )
        )


def save_coilsim1D_sensitivityOverview_2_summary(
    sim_exp,
    temperaturesKelvin,
    mol,
    store
):
    '''
    Reads COILSIM1D sensitivityOverview if it exists from the HDF5 store
    and adds the results to the complete overview of
    simulations and experiments.
    '''

    # Change mol switch to string format
    mol_dry = 'dry_wt' if mol == 0 else 'wet_mol'

    # Check if COILSIM1D simulation exists in HDF5 store
    sens_overview_v38 = 'processed/simulations/COILSIM1D_{}'.format(mol_dry)
    sens_overview_v39 = 'processed/simulations/COILSIM1D_v39_{}'.format(
        mol_dry
    )

    if sens_overview_v38 in store:
        # Read simulation results from HDF5 file
        df_sim_coilsim = store[sens_overview_v38].set_index('InChI')
        # Combine experimental results and COILSIM1D simulations
        sim_exp = utils.add2sim_exp(
            df_sim_coilsim,
            'COILSIM1D_CRACKSIM_sim',
            temperaturesKelvin,
            sim_exp
        )
        print("\nCOILSIM1D v3.8 results added to overview experiments")

    if sens_overview_v39 in store:
        # Read simulation results from HDF5 file
        df_sim_coilsim = store[sens_overview_v39].set_index('InChI')
        # Combine experimental results and COILSIM1D simulations
        sim_exp = utils.add2sim_exp(
            df_sim_coilsim,
            'COILSIM1D_CRACKSIM_39_sim',
            temperaturesKelvin,
            sim_exp
        )
        print("\nCOILSIM1D v3.9 results added to overview experiments")

    return sim_exp


def read_CHEMKINpar_summary(
    pathName,
    summaryName,
    chemkinConverters,
    mol,
    store
):
    '''
    Reads CHEMKINPar summary and adds it to the HDF5 file
    '''

    wt_mol = 'mol' if mol else 'mass'

    # Check if CHEMKIN summary file exists.
    chemkinSummary = os.path.join(pathName, 'Summary_{}_{}.xlsx'.format(
        wt_mol,
        summaryName
    ))

    if os.path.isfile(chemkinSummary):

        # Read DataFrame
        df_sim_chemkin = pd.read_excel(chemkinSummary, index_col=0)
        # Clean up DataFrame
        df_sim_chemkin.index = [
            ind.strip().replace('Mass_percentage_', '').replace(
                'Mole_percentage_', ''
            )
            for ind in df_sim_chemkin.index
            ]
        df_sim_chemkin.dropna(axis=1, how='all', inplace=True)

        # Change index of CHEMKIN converter file
        inchi2smiles_chemkin = chemkinConverters[summaryName].reset_index()
        inchi2smiles_chemkin.set_index(keys='simulation names', inplace=True)

        # Append InChI names to CHEMKIN DataFrame as index
        df_sim_chemkin.index = inchi2smiles_chemkin.loc[
            df_sim_chemkin.index, 'InChI'
        ].tolist()
        df_sim_chemkin = utils.duplicate_remover(df_sim_chemkin)
        df_sim_chemkin.index.name = 'InChI'
        # Add CHEMKIN InchI's to InchI converter

        # Save experimental names
        col_indices = store['processing/conditions'].index

        # Store CHEMKIN simulations to HDF5 file
        mol_dry = 'dry_wt' if mol == 0 else 'wet_mol'
        df_sim_chemkin.columns = col_indices
        utils.storing_2_HDF5_store(
            df_sim_chemkin.reset_index(),
            store,
            'processed/simulations/CHEMKIN_{}_{}'.format(summaryName, mol_dry),
            data_column=False
        )
        print('CHEMKIN simulation: CHEMKIN_{}_{} added to HDF5 file'.format(
            summaryName,
            mol_dry)
        )


def read_CHEMKINpar_summaries(
    pathName,
    chemkin_names,
    chemkinConverters,
    mol,
    store
):
    '''
    Reads CHEMKINPar summaries and adds them to the HDF5 file
    '''
    for summaryName in chemkin_names:
        read_CHEMKINpar_summary(
            pathName,
            summaryName,
            chemkinConverters,
            mol,
            store
        )


def save_CHEMKINpar_result_2_summary(
    summaryName,
    sim_exp,
    temperaturesKelvin,
    mol,
    store
):
    '''
    Reads CHEMKINpar summary if it exists from HDF5 file
    and adds the results to the complete comparison overview of
    simulations and experiments.
    '''
    wt_mol = 'wet_mol' if mol else 'dry_wt'

    simulation_name = 'processed/simulations/CHEMKIN_{}_{}'.format(
        summaryName,
        wt_mol
    )

    if simulation_name in store:
        # Read simulation results from HDF5 file
        df_sim_chemkin = store[simulation_name].set_index('InChI')

        # Combine experimental results and CHEMKIN simulations
        sim_exp = utils.add2sim_exp(df_sim_chemkin, 'CHEMKIN_{}_sim'.format(
            summaryName), temperaturesKelvin, sim_exp
        )
        print('\nCHEMKIN_{} compared with experiments.'.format(summaryName))
    else:
        print('\nCHEMKIN_{} not found in HDF5 file'.format(summaryName))

    return sim_exp


def save_CHEMKINpar_results_2_summary(
    chemkin_names,
    sim_exp,
    temperaturesKelvin,
    mol,
    store
):
    '''
    Iterate over all CHEMKIN summaries
    files that are given as input and available in the HDF5 file

    And adds the results to the overall
    sim_exp file.
    '''

    # Iterate over all CHEMKIN simulations
    # that are provided as input and available in HDF5 file
    for summaryName in chemkin_names:
        sim_exp = save_CHEMKINpar_result_2_summary(
            summaryName,
            sim_exp,
            temperaturesKelvin,
            mol,
            store
        )

    return sim_exp
