'''
Created on 12 Oct 2016

@author: shsymoen
'''

from dataProcessing.GC_processor_HDF5 import (
    copy2database_folder,
    process_experimental_results,
    repack,
)

from dataProcessing.utils import (
    add2sim_exp,
    check_drive_available,
    check_file_available,
    duplicate_remover,
    popup_error,
    popup_warning,
)

from plotters import (
    comparison,
    generalPlotters,
    visualiseExp,
)
from readers import simulation_output_files as sim_out_f
from writers import inputFileCreator


def runCOILSIM1DParallellizer(path_coilsim, fileDir, store, mol, v39=0):
    '''
    Copies necessary files and runs COILSIM1DParallellizer

    A distinction can be made between COILSIM1D v3.9
    (mainly in the name of the sensitivity overview file)

    Adds sensitivity results to HDF5 file
    '''
    import os
    import shutil
    import subprocess

    # Copy COILSIM1DParallellizer from N: Drive
    coilsimPar_src = r'N:\Project\918-Experimental set-ups\CoilsimParallelizer'
    if not os.path.exists(path_coilsim):
        shutil.copytree(coilsimPar_src, path_coilsim)

    # Copy input files for COILSIM1DParallellizer
    cwd = os.getcwd()
    src_files = os.listdir(fileDir)
    dst = os.path.join(
        path_coilsim,
        'Example files\simulationInputFiles\COILSIM1D'
    )
    for file_name in src_files:
        full_file_name = os.path.join(fileDir, file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy2(full_file_name, dst)

    # Run COILSIM1DParallellizer
    os.chdir(path_coilsim)
    print('COILSIM1DParallellizer started.')
    subprocess.call('run_COILSIM1DParallelizer_example.bat')
    print('COILSIM1DParallellizer stopped.')
    os.chdir(cwd)

    # Mol percentage or dry mass
    mol_mass = 'mol' if mol else 'mass'

    # Copy simulation results
    if not v39:
        src = os.path.join(
            path_coilsim,
            r'Example files\simulationSummaries\SensitivityOverview.csv'
        )
        dst = r'simulationSummaries\SensitivityOverview_{}.csv'.format(
            mol_mass
        )
    else:
        src = os.path.join(
            path_coilsim,
            r"Example files/simulationSummaries/"
            "Sensitivity_Overvieuw_Furnace.csv"
        )
        dst = r'simulationSummaries\SensitivityOverview_v39_{}.csv'.format(
            mol_mass
        )

    # Check if directory exists, otherwise create it
    # and copy simulation results
    if not os.path.exists('simulationSummaries'):
        os.makedirs('simulationSummaries')
    shutil.copy2(src, dst)


def runCHEMKINParallellizer(
    chemkin_name,
    number_of_experiments,
    path_chemkin,
    fileDir,
    mol
):
    '''Copies necessary files and runs CHEMKINParallellizer
    from folder on desktop
    '''
    import numpy as np
    import os
    import re
    import shutil
    import subprocess
    import time

    chemistry_input_location = r"N:/Project/918-Experimental set-ups/" \
        "CHEMKIN networks/{}.inp".format(chemkin_name)
    chemkinPar_src = r"N:/Project/918-Experimental set-ups" \
        "/ChemkinParallellizer"

    tags = [
        'jar_location',
        'concentrations',
        'chemistry_input',
        'total_number_of_experiments'
    ]

    if not os.path.exists(path_chemkin):
        shutil.copytree(chemkinPar_src, path_chemkin)

    # Check if results directory exists, if not create it
    if not os.path.exists('simulationSummaries'):
        os.makedirs('simulationSummaries')

    mol_mass = 'mol' if mol else 'mass'

    tag_texts = [
        path_chemkin,
        mol_mass,
        'TEMP_CHEMISTRY_INPUT',
        number_of_experiments
    ]

    # Copy input files for CHEMKINParallellizer
    cwd = os.getcwd()
    src_file = 'reactor_input.csv'
    dst = path_chemkin
    full_file_name = os.path.join(fileDir, src_file)
    shutil.copy2(full_file_name, dst)

    # Run CHEMKINParallellizer
    os.chdir(path_chemkin)
    # Get all folder names that are present at this point
    directories_old = np.array([
        dirr
        for dirr
        in os.listdir(path_chemkin)
        if 'output_' in dirr
    ])
    f = open('INPUT_auto.xml', 'r')
    g = f.read()

    for tag, tag_text in zip(tags, tag_texts):
        g = re.sub(
            r'<{}>.*</{}>'.format(tag, tag),
            r'<{}>{}</{}>'.format(tag, tag_text, tag),
            g
        )
    g = g.replace('TEMP_CHEMISTRY_INPUT', chemistry_input_location)
    f.close()
    f = open('INPUT_auto.xml', 'w')
    f.write(g)
    f.close()

    time.sleep(2)
    print('CHEMKINParallellizer started.')
    subprocess.call('run.bat')
    print('CHEMKINParallellizer stopped.')

    # Check which folder is created after the simulations
    directories_new = np.array([
        dirr
        for dirr
        in os.listdir(path_chemkin)
        if 'output_' in dirr
    ])
    result_directory = np.setdiff1d(directories_new, directories_old)
    # If create summary failed, no result directory is present
    # so check this and try createsummary once more.
    while len(result_directory) == 0:
        print("Joining results CHEMKIN in summary")
        subprocess.call('java -jar CreateSummary_06_01.jar output')
        directories_new = np.array([
            dirr
            for dirr
            in os.listdir(path_chemkin)
            if 'output_' in dirr
        ])
        result_directory = np.setdiff1d(directories_new, directories_old)

    os.chdir(result_directory[0])
    result_directory_full_path = os.getcwd()
    summary_name = 'Summary.xlsx' if mol else 'Summary_without_h2o.xlsx'
    print(summary_name)

    # Copy simulation results
    src = os.path.join(result_directory_full_path, summary_name)
    mol_mass = 'mol' if mol else 'mass'
    dest = os.path.join(
        cwd,
        r'simulationSummaries\Summary_{}_{}.xlsx'.format(
            mol_mass,
            chemkin_name
        )
    )
    print('Summary_{}_{}.xlsx copied'.format(mol_mass, chemkin_name))

    shutil.copy2(src, dest)
    os.chdir(cwd)


def check_necessary_drives_available():
    '''Verifies if the necessary drives are accessible throughout the system
    for online usage of the program. A boolean flag states the online
    useability of the program depending on the operating system.

    Parameters
    ----------
    None

    Returns
    -------
    flag : boolean
        A check flag if all necessary drives are accessible for online usage
    '''
    import sys

    a, b = True, True

    if sys.platform == 'linux':
        if not check_drive_available(r'/home/public/n_schijf'):
            popup_warning(
                "/home/public/n_schijf drive not accessible. "
                "Make sure you are connected to the UGent and LCT network "
                "(VPN if outside university) for online usage."
            )
            a = False
        # Check if N:\Project\918-Experimental set-ups
        if not check_drive_available(
            r'/home/public/n_schijf/Project/918-Experimental set-ups'
        ):
            popup_warning(
                r"'918-Experimental set-ups' project folder not accessible. "
                r"Make sure you have access rights to "
                r"'N:\Project\918-Experimental set-ups' (contact Georges"
                r" Verenghen)"
            )
            b = False
    elif 'win' in sys.platform:
        # Check if N:\ drive is available
        if not check_drive_available(r'N:'):
            popup_warning(
                "N:\ drive not accessible. "
                "Make sure you are connected to the UGent and LCT network "
                "(VPN if outside university) for online usage."
            )
            a = False
        # Check if N:\Project\918-Experimental set-ups
        if not check_drive_available(r'N:\Project\918-Experimental set-ups'):
            popup_warning(
                r"'918-Experimental set-ups' project folder not accessible. "
                r"Make sure you have access rights to "
                r"'N:\Project\918-Experimental set-ups' (contact Georges"
                r" Verenghen)"
            )
            b = False

    return a & b


def run_visualise(
        process_exp_results,
        write_chemkin,
        write_coilsim,
        visualise,
        compare,
        results_file_name,
        chemkin_names,
        mol,
        fid_sum=1,
        write_coilsim_v39=0
):
    """ Main program to visualise experimental data """
    import numpy as np
    import os
    import pandas as pd
    import sys
    idx = pd.IndexSlice

    online = check_necessary_drives_available()

    # No availability check needed (should exist in any case)
    # offline usage of COILSIM1D and/or CHEMKIN doesn't work anyway
    if write_coilsim_v39:
        coilsim_converter_location = r"N:/Project/918-Experimental set-ups/" \
            "InChI_converters/COILSIM1D_InChI_v39.xlsx"
    else:
        if online:
            coilsim_converter_location = r"N:/Project/" \
                "918-Experimental set-ups/InChI_converters/COILSIM1D_InChI.xlsx"
        else:
            coilsim_converter_location = r"InChI_converters/COILSIM1D_InChI.xlsx" 

    if 'win' in sys.platform:
        path_coilsim = r'C:\Users\{}\Desktop\COILSIM1DParallelizer'.format(
            os.getlogin()
        )
        path_chemkin = r'C:\Users\{}\Desktop\CHEMKINParallellizer'.format(
            os.getlogin()
        )
    elif sys.platform == 'linux':
        # Not implemented at the moment
        # COILSIM1D and CHEMKIN do not run on linux anyway
        pass

    # Read CHEMKIN converter files
    if (
        write_chemkin == 'True' or
        write_chemkin == '1' or
        write_chemkin == 1 or
        compare == 'True' or
        compare == '1' or
        compare == 1
    ):

        chemkinConverters = {}
        if chemkin_names != ['']:
            for chemkin_name in chemkin_names:
                if online:
                    loc = r'N:\Project\918-Experimental set-ups\InChI_converters'
                else:
                    loc = r'InChI_converters'
                chemkin_converter_location = os.path.join(
                    loc, r'ids_{}.csv'.format(chemkin_name)
                )
                if not check_file_available(chemkin_converter_location):
                    popup_error(
                        "Converter file for network {} not available in {}. "
                        "Program terminated.".format(chemkin_name, loc)
                        )
                    sys.exit()
                chemkinConverter_sub = pd.read_csv(
                    chemkin_converter_location,
                    header=0,
                    index_col='InChI'
                )
                chemkinConverters[chemkin_name] = duplicate_remover(
                    chemkinConverter_sub
                )

    yields_wt_mol = 'mol.% wet' if mol == 1 else 'wt.%'
    print(
        "The experimental data will be processed: {}\n"
        "The yields will be visualised: {}\n"
        "CHEMKIN input files will be written: {}\n"
        "COILSIM1D input files will be written: {}\n"
        "An overview comparing experiments with simulations will "
        "be created: {}\n"
        "COILSIM1D for version 3.9 will be written: {}\n"
        "The yields will be in {}".format(
            process_exp_results,
            visualise,
            write_chemkin,
            write_coilsim,
            compare,
            write_coilsim_v39,
            yields_wt_mol,
        )
    )

    # Process experimental results if necessary
    if (
        process_exp_results == 'True' or
        process_exp_results == '1' or
        process_exp_results == 1
    ):
        filename = results_file_name.split(" ")
        h5filename = '{}.h5'.format(
            process_experimental_results(filename, fid_sum)
        )
    else:
        h5filename = results_file_name

    # Read all results
    print("\nReading HDF5 results file")

    if not check_file_available(h5filename):
        popup_error(
            "{} not available in {}. Program terminated.".format(
                h5filename,
                os.getcwd()
            )
        )
        sys.exit()

    store = pd.HDFStore(h5filename)
    try:
        geometry = store['processed/geometry'].transpose()

        # Reading temperature profiles
        temperatureProfiles = store['processed/profiles/temperature']
        # Reading pressure profiles
        pressures = store['processed/profiles/pressure']
    except KeyError:  # JSR reactor does not have T,p-profiles
        temperatureProfiles = pd.DataFrame(
            np.zeros(
                (2, len(store['processing/conditions']))
            )
        )

    # Read general library
    # Check if available
    if online:
        lib_path = r"N:/Project/918-Experimental set-ups/GC_library/" \
            "GCxGC_Response_copy.xlsx"
    else:
        lib_path = r"{}/GCxGC_Response_copy.xlsx".format(
            store['processing/gen_configs'].set_index('sleutels').loc[
                'gcxgc_lib_loc'
            ].values[0]
        )
    if not check_file_available(lib_path):
        popup_error(
            "General library not available: {}. Program terminated.".format(
                lib_path
            )
        )
        sys.exit()

    gc_lib = pd.read_excel(lib_path, index_col=0)
    gc_lib = duplicate_remover(gc_lib)
    gc_lib_inchi = gc_lib.reset_index().set_index('InChI')
    gc_lib_inchi = duplicate_remover(gc_lib_inchi)

    panelname = 'mol_perc_wet' if mol == 1 else 'C_norm_yields'
    df_results = store['processed/yields/total'].loc[idx[panelname, :]]
    df_results.index.name = 'InChI'

    # Get processing conditions for all datapoints(=injections)
    datapoints = store['processing/conditions']  # .set_index('index')
    gen_configs = store['processing/gen_configs'].set_index('sleutels')
    ########################################
    # Diluents
    ########################################
    # Read it from the store
    # Save it in the store?
    ########################################
    mask_feedrates = [
        col_name
        for col_name
        in datapoints.columns
        if (('Feed_' in col_name) and ('_FeedRate' in col_name))
    ]
    feedFlow = datapoints.loc[df_results.columns, mask_feedrates]

    # Find feednames in GCxGC-library and convert to InChI
    ########################################
    # Enhancement
    ########################################
    # The converted names should be included in the HDF5_processor
    ########################################
    mask_feedrates_name = [
        col.replace('Feed_', '').replace('_FeedRate', '')
        for col
        in mask_feedrates
    ]
    feedName_InChI = gc_lib.loc[mask_feedrates_name, 'InChI']

    # Convert experimental temperatures from degrees C to Kelvin
    temp_degrees = datapoints.loc[df_results.columns, 'Temperature']
    temperaturesKelvin = [
        float(temp) + 273.15
        for temp
        in temp_degrees
    ]

    # Read sums from raw yield results
    sums = store['processed/yields/sums'].loc[idx['raw_yields', :]]

    df_results_transpose = df_results.transpose()

    # # Create InChI to Smiles converter (first time)
    # inchi2smiles_total = gc_lib[['InChI', 'Smiles']].reset_index(drop=True)
    # inchi2smiles_total.set_index(keys='InChI', inplace=True)
    # inchi2smiles_total = duplicate_remover(inchi2smiles_total)

    sim_exp = add2sim_exp(
        df_results,
        'BSSC Experiment_exp',
        temperaturesKelvin,
        pd.DataFrame()  # empty dataframe
    )

    print("HDF5 results file successfully read")

    # Visualise results
    if visualise == 'True' or visualise == '1' or visualise == 1:
        print("\nStarted creating graphs...")

        # Create discrete colormap for the temperatures
        cmap, cmap_labels = generalPlotters.create_discrete_colormap_column(
            temperaturesKelvin,
            'jet'
        )

        # Visualise the experimental data with some validation graphs
        visualiseExp.graph_experiment(
            store,
            temperatureProfiles,
            sums,
            temperaturesKelvin,
            cmap,
            cmap_labels,
            df_results_transpose,
            'Temperature (K)',
            sim_exp,
            mol
        )

    # Writing CHEMKIN parallellizer input files
    if write_chemkin == 'True' or write_chemkin == '1' or \
            write_chemkin == 1:
        print('\nWriting CHEMKIN parallelizer input file')

        ########################################
        # Enhancement
        ########################################
        # Use RDKit to compute the list of MW's
        # Can be changed in lower level in inputFileCreator
        ########################################
        MW = gc_lib.loc[mask_feedrates_name, 'Molecular Weight']

        inputFileCreator.storeMultipleCHEMKINpar_input(
            chemkin_names,
            geometry,
            chemkinConverters,
            feedName_InChI,
            MW,
            temperatureProfiles,
            pressures,
            feedFlow
            )
        print('\nCHEMKIN parallellizer input files successfully written.')

        print('Starting CHEMKIN simulations.')
        for chemkin_name in chemkin_names:
            fileDir = r'simulationInputFiles\{}'.format(chemkin_name)
            number_of_experiments = len(temperatureProfiles.transpose())
            runCHEMKINParallellizer(
                chemkin_name,
                number_of_experiments,
                path_chemkin,
                fileDir,
                mol
            )
        print('CHEMKIN simulations successfully finished.')

        # Read CHEMKIN simulation results and copy to HDF5 store
        src = 'simulationSummaries'
        sim_out_f.read_CHEMKINpar_summaries(
            src,
            chemkin_names,
            chemkinConverters,
            mol,
            store
        )
        print('\nReady for comparison between COILSIM1D and experiments.')

    # Writing COILSIM1D parallellizer input files
    run_coilsim = (
        write_coilsim == 'True'
        or write_coilsim == '1'
        or write_coilsim == 1
    )
    run_coilsim_v39 = (
        write_coilsim_v39 == 'True'
        or write_coilsim_v39 == '1'
        or write_coilsim_v39 == 1
    )
    if (run_coilsim or run_coilsim_v39):

        v39 = int(run_coilsim_v39)

        print('\nWriting COILSIM1D parallellizer input files.')

        # Read converter files
        coilsimConverter = pd.read_excel(
            coilsim_converter_location,
            header=0,
            index_col='InChI'
        )
        coilsimConverter = duplicate_remover(coilsimConverter)

        ########################################
        # diluent
        ########################################
        # No other diluent possible in COILSIM1D
        # so at the moment kept at water
        ########################################
        diluent = 'H2O'
        mask_not_diluent = [
            name
            for name
            in feedName_InChI.index
            if name != diluent
        ]

        ########################################
        # Enhancement
        ########################################
        # It should be possible to include multiple feeds
        # At the moment only one is possible
        ########################################
        if len(feedName_InChI[mask_not_diluent].values) > 1:
            popup_error(
                "At the moment it is not possible to have more "
                "than one HC for COILSIM1D"
            )
            sys.exit()

        inputData = temperatureProfiles, pressures, feedFlow.transpose(), \
            geometry, coilsimConverter, \
            feedName_InChI[mask_not_diluent].values[0]
        fileDir = r'simulationInputFiles\COILSIM1D'

        reactor = gen_configs.loc['Reactor', :].values[0]
        inputFileCreator.storeCOILSIM1Dpar_input(
            fileDir,
            inputData,
            reactor,
            mol,
            v39
        )
        print(
            '\nCOILSIM1D parallellizer input files are successfully written.'
        )

        # Run COILSIM1D simulations
        print('\nStarting COILSIM1D.')
        runCOILSIM1DParallellizer(path_coilsim, fileDir, store, mol, v39)
        print('COILSIM1D simulations finished.')

        # Read COILSIM1D sensitivityOverview and copy to HDF5 store
        src = 'simulationSummaries'
        ########################################
        # Enhancement
        ########################################
        # store input necessary?
        ########################################
        sim_out_f.read_coilsim1D_sensitivityOverview(
            src,
            coilsimConverter,
            mol,
            store,
            v39
        )
        print('\nReady for comparison between COILSIM1D and experiments.')

    # Writing comparison overview file
    if compare == 'True' or compare == '1' or compare == 1:

        if not os.path.exists('Comparison'):
            os.makedirs('Comparison')
        if not os.path.exists(r"Comparison\Figures"):
            os.makedirs(r'Comparison\Figures')

        # Read converter files
        coilsimConverter = pd.read_excel(
            coilsim_converter_location,
            header=0,
            index_col='InChI'
        )
        coilsimConverter = duplicate_remover(coilsimConverter)

        # Sum trans- and cis-2-butene to 2-butene
        sim_exp['1S/C4H8/c1-3-4-2/h3-4H,1-2H3'] = \
            sim_exp['1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'] \
            + sim_exp['1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-']
        print('Trans- and cis-2-butene summed to 2-butene for comparison.')


        inchi2iupac = coilsimConverter['IUPAC']

        # Add COILSIM1D simulation results
        # if sensitivityOverview exists.
        sim_exp = sim_out_f.save_coilsim1D_sensitivityOverview_2_summary(
            sim_exp,
            temperaturesKelvin,
            mol,
            store
            )

        # Add CHEMKIN simulation results
        if chemkin_names != ['']:
            sim_exp = sim_out_f.save_CHEMKINpar_results_2_summary(
                chemkin_names,
                sim_exp,
                temperaturesKelvin,
                mol,
                store
            )

        print("\nWriting comparison file")
        try:
            diluent = 'H2O'
            mask_not_diluent = [
                name for name in feedName_InChI.index if name != diluent
            ]
            sim_exp['Conversion'] = 100 - sim_exp.loc[
                :, feedName_InChI[mask_not_diluent]
            ]
        except:
            print("Conversion could not be calculated.")
        sim_exp.columns.name = 'InChI'
        sim_exp.to_csv(r"Comparison\comparison_exp_sim.csv")
        print("Comparison file is written.")

        # Start making comparison figures
        comparison.create_comparison_figures(
            sim_exp,
            inchi2iupac,
            feedName_InChI,
            mol
        )
        print(
            "\nComparison figures are successfully created "
            "and stored in the 'Comparison' folder"
        )

    # Closes the HDF5 file container when exiting the program.
    store.close()

    # Reduce size of HDF5 file by repacking and compressing
    if write_chemkin or write_coilsim or process_exp_results:
        if h5filename.endswith('.h5'):
            h5filename = h5filename.replace('.h5', '')

        repack(h5filename)

    # Copy to designated folder
    if online:
        if 'win' in sys.platform:
            copy2database_folder(h5filename, gen_configs)
        elif sys.platform == 'linux':
            pass
