'''
Created on 22 Sep 2017

@author: shsymoen
'''

import sys
# Add some more folders to the Python path
sys.path.append('..')
sys.path.append('.')

# Local imports
from dataProcessing.GC_processor_HDF5 import (
    calculate_areas_2_raw,
    calculate_elemental_balance,
    calculate_mass_2_mol,
    calculate_norm_2_wet,
    calculate_normalised,
    processing_chrom_gcxgc_fid_scd_ncd,
    processing_chrom_loa,
    processing_chrom_rga,
    processing_chroms_2_raw,
    processing_elemental_balances,
    processing_log_file,
    processing_norm_2_wet,
    processing_norm_dry_molar,
    processing_norm_wet_molar,
    process_start_end_log_file,
    processing_raw_2_carbon_normalized,
    processing_raw_2_norm,
    read_gcxgc_lib,
    reading_agilent_data,
    reading_cal_factor,
    reading_hyperchrom_data,
    reading_input_conditions,
    reading_loggingfile,
    reading_reactor_config,
    reading_storing_GCxGC_fid_data,
    reading_storing_LOA_data,
    reading_storing_RGA_data,
    reading_storing_cal_factors,
    store_processing_conditions_hdf,
    sum_generator,
    yield_data,
)

from dataProcessing.utils import (
    inchi_mol_wt
)


def check_processed_yields(exp_yields, sim_yields):
    '''Asserts the processed yields are equal between two HDF5 files

    Parameters
    ----------
    exp_yields : DataFrame
    sim_yields : DataFrame

    Returns
    -------
    None
        If the processed yields are identical
        Throws assertion error if not identical
    '''
    import pandas as pd

    for sheet in exp_yields.index.get_level_values(0).unique():
        print(sheet)  # Easier for bug detection

        # Get expected yield
        exp_yield = exp_yields.loc[sheet]

        # Extract processed yields
        sim_yield = sim_yields.total[sheet]
        sim_yield.index.name = 'molecules'

        # Assert expected and processed are equal
        pd.testing.assert_frame_equal(
            sim_yield.sort_index(),
            exp_yield.sort_index()
        )


def create_check_processed_yields(exp_name):
    '''Creates the yields for a certain experiment and compares
    the newly created yields against a reference.

    Parameters
    ----------
    exp_name : string
        The name of the folder in which all the data can be found.

    Returns
    -------
    None
        If all yields are identical with the reference, otherwise
        an assertion error is raised.
    '''
    import os
    import pandas as pd

    # HDF5 file name
    h5filename = r"{}_sim".format(exp_name)
    h5filename_exp = r"{}_exp".format(exp_name)

    # Read the expected yields from the HDF5 file
    h5filepath = os.path.join(
        r'testFiles_dataProcessing/'
        'yield_processor_test/{}'.format(exp_name),
        '{}.h5'.format(h5filename_exp)
    )
    exp_yields = pd.read_hdf(
        h5filepath,
        "processed/yields/total"
    )

    # Read in the conditions and gen_configs
    filename = "conditions.csv gen_configs.csv".split(' ')
    filename = [
        r'testFiles_dataProcessing/yield_processor_test/{}/{}'.format(
            exp_name,
            file
        )
        for file in filename
    ]
    conditions, gen_configs = reading_input_conditions(filename)

    store_processing_conditions_hdf(conditions, gen_configs, h5filename)

    # Reading and storing RGA calibration factors
    CFs = reading_storing_cal_factors(h5filename, gen_configs)

    # Read in the RGA data
    rga_channels = reading_storing_RGA_data(h5filename)

    # Reading GCxGC library (NOT STORED in HDF5 file)
    gcxgc_lib = read_gcxgc_lib(gen_configs)

    # Read in the GCxGC data
    reading_storing_GCxGC_fid_data(
        h5filename,
        conditions,
        gen_configs,
        detector_switch='FID'
    )

    # Read in the LOA data if applicable
    if gen_configs['Process_LOA'] == 'TRUE':
        reading_storing_LOA_data(
            h5filename,
            conditions,
            gen_configs
        )

    # Create yield data instance
    exp_data = yield_data(
        rga_channels,
        conditions,
        CFs,
        gen_configs,
        gcxgc_lib
    )

    # Extract area and volume data from raw chromatograms
    processing_chrom_rga(h5filename, rga_channels, conditions)
    processing_chrom_gcxgc_fid_scd_ncd(h5filename, conditions)
    processing_chrom_loa(h5filename, conditions)

    # Process the chromatograms to raw yields
    fid_sum = True
    processing_chroms_2_raw(
        exp_data,
        rga_channels,
        gcxgc_lib,
        fid_sum,
        h5filename
    )

    # Normalise raw yields
    processing_raw_2_norm(exp_data, fid_sum)
    # Convert the normalised yields to wet yields
    processing_norm_2_wet(exp_data, fid_sum)
    # Calculate elemental balances
    processing_elemental_balances(exp_data, gcxgc_lib, fid_sum)
    # Normalise yields based on carbon balance
    processing_raw_2_carbon_normalized(exp_data, fid_sum)
    # Convert normalised dry results to molar dry
    processing_norm_dry_molar(exp_data, fid_sum)
    # Convert normalised wet results to molar wet
    processing_norm_wet_molar(exp_data, fid_sum)

    # Check all processed yields with true data
    check_processed_yields(exp_yields, exp_data)

    # Remove the created HDF5 file
    os.remove('{}.h5'.format(h5filename))
    # Remove the created GCxGC_fid_index.csv file
    os.remove('GCxGC_fid_index.csv')


class TestGC_processor(object):

    def test_reading_agilent_data(self):
        """ Automatic pytest for reading_agilent_data """

        import numpy as np
        import pandas as pd

        # Chromatogram with data
        sim_output = reading_agilent_data(
            r'testFiles_dataProcessing',
            r'fid_agilent.csv'
        )
        exp_output = pd.read_csv(
            r'testFiles_dataProcessing/exp_out_agilent_reduced.csv',
            index_col=0
        )
        exp_output.index = [ind.strip() for ind in exp_output.index]

        pd.testing.assert_frame_equal(sim_output, exp_output)

        # Empty chromatogram
        sim_output = reading_agilent_data(
            r'testFiles_dataProcessing',
            r'fid_agilent_empty.csv'
        )
        column_names = ['Area_perc', 'Ret_Time', 'Area']
        exp_output = pd.DataFrame(columns=column_names)

        # Change to numeric
        exp_output['Area_perc'] = exp_output['Area_perc'].astype(np.int64)
        exp_output['Area'] = exp_output['Area'].astype(np.int64)
        exp_output['Ret_Time'] = exp_output['Ret_Time'].astype(np.float64)

        pd.testing.assert_frame_equal(sim_output, exp_output)

    def test_reading_hyperchrom_data(self):
        """ Automatic pytest for reading_hyperchrom_data """
        import pandas as pd

        # Chromatogram with data
        sim_output = reading_hyperchrom_data(
            r'testFiles_dataProcessing',
            r'propane_800_FID_Run1.XLS'
        )
        exp_output = pd.read_csv(
            r'testFiles_dataProcessing/exp_out_hyperchrom.csv',
            index_col=0
        )
        exp_output.index = [ind.strip() for ind in exp_output.index]

        pd.testing.assert_frame_equal(sim_output, exp_output)

        # Empty chromatogram
        sim_output = reading_hyperchrom_data(
            r'testFiles_dataProcessing',
            r'empty_chrom.XLS'
        )
        exp_output = pd.DataFrame(columns=['Area_perc', 'Ret_Time', 'Area'])
        for col in exp_output.columns:
            exp_output[col] = pd.to_numeric(exp_output[col])

        pd.testing.assert_frame_equal(sim_output, exp_output)

    def test_reading_cal_factor(self):
        """ Automatic pytest for reading_cal_factor """
        # INCOMPLETE
        # Add additional response files as input

        import pandas as pd

        sim_output = reading_cal_factor(
            r'testFiles_dataProcessing',
            r'LOA_Response_MTHF_12082016.xlsx'
        )
        exp_output = pd.read_csv(
            r'testFiles_dataProcessing/output_reading_cal_factor.csv',
            index_col=0
        )

        pd.testing.assert_frame_equal(sim_output, exp_output)

    def test_reading_input_conditions(self):
        """ Automatic pytest for reading_input_conditions """
        # UNFINISHED
        # Assert all necessary input in the gen_configs file
        # Check conditions file

        filename = "conditions_propane.csv gen_configs_propane.csv".split(' ')
        filename = [
            r'testFiles_dataProcessing/{}'.format(file) for file in filename
        ]
        conditions, gen_configs = reading_input_conditions(filename)

        assert 'Reactor' in gen_configs


    def test_reading_reactor_config(self):
        """ Automatic pytest for reading_reactor_config """
        # UNFINISHED
        # Not only BSSC metal should be tested
        import pandas as pd
        import numpy as np

        loggingFileNames_exp = np.array([
            'Coriflows',
            'Pressures',
            'Temperatures_Elicra',
            'Temperatures_Nabertherm',
            'Temperatures_Samplebox'
        ])

        decimal_exp = '.'
        geometry_exp = pd.DataFrame(
            data=[1.49, 0.006, 0.002],
            index=['length (m)', 'diameter (m)', 'thickness (m)'],
        )

        temperature_axial_exp = pd.DataFrame(
            data=[0.19, 0.38, 0.48, 0.67, 0.86, 1.05, 1.24, 1.43],
            index=[
                "Inlet_Temp_Reactor.Zone_1",
                "Inlet_Temp_Reactor.Zone_2",
                "Inlet_Temp_Reactor.Zone_3",
                "Inlet_Temp_Reactor.Zone_4",
                "Inlet_Temp_Reactor.Zone_5",
                "Inlet_Temp_Reactor.Zone_6",
                "Inlet_Temp_Reactor.Zone_7",
                "Inlet_Temp_Reactor.Zone_8",
            ],
            columns=['axial_position[m]']
        )
        temperature_axial_exp.index.name = 'Temperature_name'

        temperature_names_exp = temperature_axial_exp.index.values

        pressure_axial_exp = pd.DataFrame(
            data=[0, 1.49],
            index=[
                "PI-05 Elicra",
                "PI-07 Outlet Reactor",
            ],
            columns=['axial_position[m]']
        )
        pressure_axial_exp.index.name = 'Pressure_name'

        pressure_names_exp = pressure_axial_exp.index.values

        keys_exp = [
            "loggingFileNames",
            "decimal",
            "geometry",
            "log_name_temperature",
            "log_name_pressure",
            "axial_temperatures",
            "axial_pressures",
            "names_temperature",
            "names_pressure",
        ]

        reactor_config = reading_reactor_config(
            r'testFiles_dataProcessing', r'BSSC_metal.xlsx'
        )

        # Check if everything is there
        for key in keys_exp:
            assert key in reactor_config

        # Check if decimal character is the same
        assert reactor_config['decimal'] == decimal_exp
        assert reactor_config['log_name_temperature'] == "Temperatures_Elicra"
        assert reactor_config['log_name_pressure'] == "Pressures"
        assert reactor_config['divide_by'] == 1
        assert reactor_config['barg'] == 1
        np.testing.assert_array_equal(
            pressure_names_exp,
            reactor_config['names_pressure'],
        )
        pd.testing.assert_frame_equal(
            temperature_axial_exp,
            reactor_config['axial_temperatures']
        )
        pd.testing.assert_frame_equal(
            pressure_axial_exp,
            reactor_config['axial_pressures']
        )
        np.testing.assert_array_equal(
            loggingFileNames_exp,
            reactor_config['loggingFileNames'],
        )
        np.testing.assert_array_equal(
            temperature_names_exp,
            reactor_config['names_temperature'],
        )
        pd.testing.assert_frame_equal(
            geometry_exp,
            reactor_config['geometry']
        )

    def test_sum_generator(self):
        """ Automatic pytest for sum_generator """

        import pandas as pd

        df_fid = pd.DataFrame(
            {'800-1': [1.1, 2], '800-2': [2, 3]},
            index=['CH4', 'C3H6']
        )
        df_tcd = pd.DataFrame(
            {'800-1': [5.1, 7, 45], '800-2': [2, 8, 45]},
            index=['CH4', 'C2H4', 'N2']
        )

        df_gcxgc = pd.DataFrame(
            {'800-1': [10, 12], '800-2': [0.5, 0.5]},
            index=['benzene', 'toluene']
        )

        df_loa = pd.DataFrame(
            {'800-1': [10.1, 4.0, 13.1], '800-2': [0.5, 0.0, 76]},
            index=['CH4', 'H2O', '1S/CH4CO/c1-2/h2H,1H3']
        )

        list_1 = {
            'FID': df_fid,
            'TCD': df_tcd,
            'GCxGC': df_gcxgc
        }

        # First test FID first
        sums_exp_1 = pd.DataFrame(
            {'800-1': [10.1, 22, 32.1], '800-2': [13.0, 1, 14]},
            index=['C4', 'C5', 'Total']
        )
        total_exp_1 = pd.DataFrame(
            {'800-1': [1.1, 2, 7, 10, 12.0], '800-2': [2, 3, 8, 0.5, 0.5]},
            index=['CH4', 'C3H6', 'C2H4', 'benzene', 'toluene']
        )
        sums, total = sum_generator(
            list_1,
            internal_standards=['N2'],
            fid_sum=1
        )
        pd.testing.assert_frame_equal(sums_exp_1, sums)
        pd.testing.assert_frame_equal(total_exp_1, total)

        # Second test TCD first
        sums_exp_2 = sums_exp_1.copy()
        total_exp_2 = total_exp_1.copy()
        sums_exp_2.loc['C4', '800-1'] = 14.1
        sums_exp_2.loc['Total', '800-1'] = 36.1
        total_exp_2.loc['CH4', '800-1'] = 5.1
        sums, total = sum_generator(
            list_1,
            internal_standards=['N2'],
            fid_sum=0
        )
        pd.testing.assert_frame_equal(
            sums_exp_2.sort_index(),
            sums.sort_index()
        )
        pd.testing.assert_frame_equal(
            total_exp_2.sort_index(),
            total.sort_index()
        )

        # Third test add LOA
        list_2 = list_1.copy()
        list_2['LOA'] = df_loa
        total_exp_3 = total_exp_1.copy()
        total_exp_3.loc['H2O', :] = [4.0, 0.0]
        sums, total = sum_generator(
            list_2,
            internal_standards=['N2'],
            fid_sum=1
        )
        pd.testing.assert_frame_equal(
            sums_exp_2.sort_index(),
            sums.sort_index()
        )
        pd.testing.assert_frame_equal(
            total_exp_3.sort_index(),
            total.sort_index()
        )

        # Fourth test add SCD
        list_3 = list_1.copy()
        list_3['GCxGC-SCD'] = df_loa.loc[['CH4', 'H2O'], :]
        sums_exp_3 = sums_exp_1.copy()
        sums_exp_3.loc['C5', '800-1'] = 26.0
        sums_exp_3.loc['Total', '800-1'] = 36.1
        sums, total = sum_generator(
            list_3,
            internal_standards=['N2'],
            fid_sum=1
        )
        pd.testing.assert_frame_equal(
            sums_exp_3.sort_index(),
            sums.sort_index()
        )
        pd.testing.assert_frame_equal(
            total_exp_3.sort_index(),
            total.sort_index()
        )

    def test_calculate_mass_2_mol(self):
        """ Automatic pytest for calculate_mass_2_mol """
        import pandas as pd
        import numpy as np

        excelfile = pd.ExcelFile(r'testFiles_dataProcessing/mass_2_mol.xlsx')

        norm_mass = excelfile.parse('norm_mass')
        mol_perc = excelfile.parse('mol_perc')
        sub = excelfile.parse('sub')
        # Remove sum from sheet (has no index)
        sub.drop(np.nan, inplace=True)

        # Compute molecular weights
        # mol_wts = norm_mass.index.map(lambda x: inchi_mol_wt(x))
        mol_wts = [28.05, 16.04, 42.08, 28.01]
        # Put them in a series
        mws = pd.Series(
            index=norm_mass.index,
            data=mol_wts,
            name='MWs'
        )
        sim = calculate_mass_2_mol(norm_mass, mws)
        # drop nitrogen in simulated DataFrame
        sim.drop('1S/N2/c1-2', inplace=True)
        pd.testing.assert_frame_equal(sim, sub)

    def test_calculate_normalised(self):
        """ Automatic pytest for calculate_normalised """
        import pandas as pd

        # Series
        sheet = pd.Series([2, 3, 5])
        sums = 10
        sim = calculate_normalised(sheet, sums)
        exp = pd.Series([20.0, 30, 50])
        pd.testing.assert_series_equal(sim, exp)

        # DataFrame
        sheet = pd.DataFrame([[1, 3, 4, 6, 6], [1, 1, 1, 1, 6]]).transpose()
        sums = sheet.sum()
        exp = pd.DataFrame(
            [[5, 15, 20, 30, 30.0], [10.0, 10, 10, 10, 60]]
        ).transpose()
        sim = calculate_normalised(sheet, sums)
        pd.testing.assert_frame_equal(sim, exp)

        # Different normalising factor
        sim = calculate_normalised(sheet, sums, normalising_factor=80)
        exp = pd.DataFrame(
            [[4, 12, 16, 24, 24.0], [8.0, 8, 8, 8, 48]]
        ).transpose()

    def test_calculate_norm_2_wet(self):
        """ Automatic pytest for calculate_norm_2_wet """
        import pandas as pd

        # Create the input
        norm_yields = pd.DataFrame(
            [[5, 10.0], [15, 10], [20, 10], [30, 10], [30.0, 60.0]],
            columns=["800-1", "800-2"],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )
        HC_feed = pd.Series(
            [200, 200],
            index=["800-1", "800-2"]
        )
        total_feedrate = pd.Series(
            [400, 400],
            index=["800-1", "800-2"]
        )

        # Create the expected output
        exp = pd.DataFrame(
            [[2.5, 5.0], [7.5, 5], [10, 5], [15, 5], [15.0, 30.0]],
            columns=["800-1", "800-2"],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )

        # Do the calculations with created input
        sim = calculate_norm_2_wet(norm_yields, HC_feed, total_feedrate)

        # Assert calc results equal expected
        pd.testing.assert_frame_equal(sim, exp)

    def test_calculate_elemental_balance(self):
        """ Automatic pytest for calculate_norm_2_wet """
        import pandas as pd

        # Create inputs
        mass_yields = pd.DataFrame(
            [[5, 10.0], [15, 10], [20, 10], [30, 10], [30.0, 60.0]],
            columns=["800-1", "800-2"],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )
        mass_wt = pd.Series(
            [10, 5, 2, 10, 30],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )
        atoms_in_comp = pd.Series(
            [1, 2, 3, 0, 2],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )
        mass_atom = 12.0

        # Calculate the output
        sim = calculate_elemental_balance(
            mass_yields,
            mass_wt,
            atoms_in_comp,
            mass_atom
        )

        # Create expected output
        exp = pd.DataFrame(
            [[6, 12.0], [72, 48], [360, 180], [0, 0], [24.0, 48.0]],
            columns=["800-1", "800-2"],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )

        # Assert calculation equal to expected
        pd.testing.assert_frame_equal(sim, exp)

        # Different mass of atom
        mass_atom = 6.0

        # Calculate the output
        sim = calculate_elemental_balance(
            mass_yields,
            mass_wt,
            atoms_in_comp,
            mass_atom
        )

        # Create expected output
        exp = pd.DataFrame(
            [[3, 6.0], [36, 24], [180, 90], [0, 0], [12.0, 24.0]],
            columns=["800-1", "800-2"],
            index=["methane", "ethene", "propene", "hydrogen", "ethane"]
        )

        # Assert calculation equal to expected
        pd.testing.assert_frame_equal(sim, exp)

    def test_calculate_areas_2_raw(self):
        """ Automatic pytest for calculate_areas_2_raw """

        import pandas as pd

        # Test RGA TCD channel
        raw_areas_dct = {
            '800-1': [1671119, 1163770, 1774434],
            '800-2': [3228769, 2256780, 3430734]
            }
        component_names = ['C2H4', 'CH4', 'N2']
        raw_areas = pd.DataFrame(raw_areas_dct, index=component_names)
        channel_RFs = pd.Series(
            {'C2H4': 1.158, 'CH4': 1.000, 'N2': 1.485},
            name='Response'
        )
        tcd_constant = pd.Series(
            {'800-1': 1.459622e-07, '800-2': 7.549412e-08}
        )
        expected_output = pd.DataFrame(
            {
                '800-1': [28.245956, 16.986641, 38.461538],
                '800-2': [28.226608, 17.037363, 38.461538]
            }, index=component_names
        )

        pd.testing.assert_frame_equal(
            calculate_areas_2_raw(raw_areas, channel_RFs, tcd_constant),
            expected_output
            )

        # Test RGA FID channel (the same for all other detectors)
        raw_areas_dct = {
            '800-1': [1177721, 663345, 819952],
            '800-2': [2327519, 1335464, 1701984]
            }
        component_names = ['C2H4', 'CH4', 'C3H6']
        raw_areas = pd.DataFrame(raw_areas_dct, index=component_names)
        channel_RFs = pd.Series(
            {'C2H4': 0.962940, 'CH4': 1.000, 'C3H6': 0.835000},
            name='Response'
            )
        is_remote_area = pd.Series({'800-1': 1163770, '800-2': 2256780})
        is_remote_response = pd.Series({'800-1': 1, '800-2': 1})
        is_channel_area = pd.Series({'800-1': 663345, '800-2': 1335464})
        is_channel_RF = pd.Series({'800-1': 1, '800-2': 1})
        expected_output = pd.DataFrame(
            {
                '800-1': [29.040862, 17.532464, 16.986641],
                '800-2': [28.593180, 18.130598, 17.037363],
                }, index=['C2H4', 'C3H6', 'CH4']
            )
        sim_output = calculate_areas_2_raw(
            raw_areas,
            channel_RFs,
            tcd_constant,
            is_remote_area,
            is_channel_area,
            is_channel_RF,
            is_remote_response
        )

        pd.testing.assert_frame_equal(
            sim_output.sort_index(),
            expected_output.sort_index()
        )


    def test_processing_yields(self):
        """Automatic PyTest for processing yields"""
        exp_names = [
            'propane',
        ]

        for exp_name in exp_names:
            print('\n\nFolder being examined: {}...'.format(exp_name))
            create_check_processed_yields(exp_name)
            print('\n\n{} folder examined.'.format(exp_name))
