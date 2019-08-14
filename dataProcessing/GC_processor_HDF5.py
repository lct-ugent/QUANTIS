'''
Created on 12 Dec 2016

@author: shsymoen
'''
from dataProcessing.utils import (
    add_excel_formatting,
    check_diff_index,
    check_file_available,
    diff_index,
    duplicate_remover,
    inchi_2_smiles,
    popup_entry,
    popup_error,
    storing_2_HDF5_store,
)

GC_processor_HDF5_version = 'v0.3.0'


class yield_data():

    def __init__(self, rga_channels, datapoints, CFs, gen_configs, gcxgc_lib):
        '''
        Constructor
        '''
        import sys
        ######################################
        # Can all the processing modules be implemented in this class?
        ######################################

        # Input files
        ######################################
        self.rga_channels = rga_channels
        self.datapoints = datapoints
        self.gen_configs = gen_configs

        # Initialisations
        self.CFs = CFs
        self.gcxgc_lib = gcxgc_lib
        self.rga_tcd = {}
        self.rga_fid = {}
        self.loa = {}
        self.gcxgc_fid = {}
        self.sums = {}
        self.total = {}

        self.fid_raw_data = {}
        self.tcd_raw_data = {}
        self.loa_raw_data = {}
        self.gcxgc_fid_raw_data = {}
        self.gcxgc_scd_raw_data = {}
        self.gcxgc_ncd_raw_data = {}

        # Initialisations internal standards
        self.primary_is_response = {}
        self.primary_is_area = {}
        self.sec_is_area_tcd = {}
        self.sec_is_response = {}
        self.sec_is_response_tcd = {}
        self.sec_is_area_fid = {}
        self.tert_is_area = {}
        self.tert_is_response_tcd = {}
        self.tert_is_area_gcxgc = {}
        self.tert_is_response = {}
        self.tert_in_tcd = {}
        self.quat_is_area_fid = {}
        self.quat_is_response = {}
        self.quat_is_response_fid = {}
        self.quat_is_area_loa = {}
        self.quat_in_tcd = {}

        # Detectors
        self.detectors = {
            'FID': self.rga_fid,
            'TCD': self.rga_tcd,
            'GCxGC': self.gcxgc_fid,
            }
        # Add additional detectors if available
        if 'LOA' in self.datapoints:
            self.detectors['LOA'] = self.loa
        if 'SCD' in self.datapoints:
            self.detectors['SCD'] = self.gcxgc_scd
        if 'NCD' in self.datapoints:
            self.detectors['NCD'] = self.gcxgc_ncd

        self.is_inchi = None

        # Feeds
        ######################################
        # experimental names of feeds
        HC_feed_names_iupac = [
            gen_configs[feed_tag]
            for feed_tag
            in gen_configs.index
            if 'Feed_' in feed_tag
        ]
        dil_names_iupac = [
            gen_configs[dil_tag]
            for dil_tag
            in gen_configs.index
            if 'Diluent_' in dil_tag
        ]
        # mask for conditions.csv
        mask_feedrates = [
            "Feed_{}_FeedRate".format(feedname)
            for feedname
            in HC_feed_names_iupac
        ]
        mask_dilrates = [
            "Feed_{}_FeedRate".format(dilname)
            for dilname
            in dil_names_iupac
        ]

        # Get the InChI names from the GCxGC library
        try:
            import numpy as np
            HC_feed_names = gcxgc_lib.loc[HC_feed_names_iupac, 'InChI'].values
            HC_feed_names = np.unique(HC_feed_names)
        except:
            print(
                "Given feednames in 'gen_configs'"
                " are not available in GCxGC library"
            )
            sys.exit(0)
        try:
            dil_feed_names = gcxgc_lib.loc[dil_names_iupac, 'InChI'].values
        except:
            print(
                "Given diluent names in 'gen_configs'"
                " are not available in GCxGC library"
            )
            sys.exit(0)

        # Get the feedrates from datapoints
        try:
            HC_feed_rates = self.datapoints.loc[:, mask_feedrates]
            HC_feed_rates.columns = HC_feed_names
        except:
            print(
                "Mismatch feed names in gen_configs.csv and in conditions.csv"
            )
            sys.exit(0)

        try:
            dil_feed_rates = self.datapoints.loc[:, mask_dilrates]
            dil_feed_rates.columns = dil_feed_names
        except:
            print(
                "Mismatch diluent names in gen_configs.csv"
                "and in conditions.csv"
            )
            sys.exit(0)

        self.diluents = feed(dil_feed_names, dil_feed_rates, dil_names_iupac)
        self.HC_feeds = feed(HC_feed_names, HC_feed_rates, HC_feed_names_iupac)

    def get_detectors_sheet(self, sheetname):
        ''' Returns a dict with detectors and associated DataFrames for a
        given sheetname

        Extended descript

        Parameters
        ----------
        self :
        sheetname : String
            Name of the sheetname

        Returns
        -------
        detectors_sheet : Dict
            A list of the detectors with associated DataFrames
        '''
        detectors_sheet = {
            detector: self.detectors[detector][sheetname]
            for detector
            in self.detectors
            }

        return detectors_sheet

    def get_CFs_inchi(self):
        ''' Returns a dict with the CFs DataFrames with the InChI's as index

        Extended description

        Parameters
        ----------
        self :

        Returns
        -------
        cal_factors_inchi : Dict
        '''
        ######################################
        # Enhancement
        ######################################
        # cal_factors_inchi should be stored
        # and created during init
        ######################################

        # initialise
        cal_factors_inchi = {}

        # RGA and LOA calibration factors
        for detector in self.CFs:
            cal_factors_inchi[detector] = self.CFs[detector].set_index('InChI')

        # GCxGC calibration factors
        # take only components present in index of raw_yields
        gcxgc_lib_inchi = self.gcxgc_lib.set_index('InChI')  # .loc[
#             self.gcxgc_fid['raw_yields'].index, :
#         ]
        gcxgc_lib_inchi = duplicate_remover(
            gcxgc_lib_inchi,
            first_last='first'
            )
        # rename column molecular weight to MW
        converter = {
            '#{}'.format(el): el for el in 'CHONS'
            }
        converter['Molecular Weight'] = 'MW'
        gcxgc_lib_inchi.rename(
            columns=converter,
            inplace=True
            )
        # Add to cal_factors_inchi
        cal_factors_inchi['GCxGC'] = gcxgc_lib_inchi

        return cal_factors_inchi

    def get_fractions(self):
        ''' Returns the feed fractions '''
        import pandas as pd

        flowrates_seperate = pd.concat(
            [self.diluents.get_flows(), self.HC_feeds.get_flows()], axis=1
        )
        total_flow = self.HC_feeds.get_flow_total() \
            + self.diluents.get_flow_total()

        return flowrates_seperate.divide(total_flow, axis=0)

    def get_feed_total_flow(self):
        '''Returns the total feedflow rate. Diluent and HCs.

        Parameters
        ----------
        self

        Returns
        -------
        total_flow : float
            Total flow rate of the inflow
        '''
        total_flow = self.HC_feeds.get_flow_total() \
            + self.diluents.get_flow_total()
        return total_flow

    def store_yield_data_to_hdf5(self, h5filename):
        ''' Converts the yield data and stores it to the HDF5 file '''

        import pandas as pd

        subgroups = {
            'RGA_TCD': self.rga_tcd,
            'RGA_FID': self.rga_fid,
            'GCxGC_FID': self.gcxgc_fid,
            'sums': self.sums,
            'total': self.total,
            'HC_inflow': self.HC_feeds.inflow,
            'Diluent_inflow': self.diluents.inflow
            }

        # Add LOA data if it is not empty
        if self.loa:
            subgroups['LOA'] = self.loa

        for subgroupname in subgroups:
            # Convert to pandas MultiIndexed DataFrame
            # Necessary to sort for slicing afterwards
            subgrp_panel = pd.concat(subgroups[subgroupname]).sort_index()
            subgrp_panel.index.names = ['sheets', 'molecules']

            # Let's hope these weren't necessary
#             subgrp_panel.dropna(axis=0, how='all', inplace=True)
#             subgrp_panel.dropna(axis=1, how='all', inplace=True)

            # Store to HDF5 file
            storing_2_HDF5_store(
                subgrp_panel,
                '{}.h5'.format(h5filename),
                'processed/yields/{}'.format(subgroupname),
            )
            print("\t{} appended to HDF5 file.".format(subgroupname))

        print(
            'Yield_data successfully appended to HDF5 {}.h5'.format(
                h5filename
            )
        )

    def read_prim_IS_data(self, inj_name):
        '''Read all the necessary information of the primary Internal Standard
        for the TCD channel
        (multiplying might be a bit faster rather than dividing in next step)

        Parameters
        ----------
        self : yield_data instance
        inj_name : String
            Name of the injection tag
        '''
        prim_is_name = self.datapoints.loc[
            inj_name, 'PrimaryInternalStandard_Name'
        ]
        self.primary_is_response[inj_name] = 1 / self.CFs['TCD'].loc[
            prim_is_name, 'Response'
        ]

        # Read primaryIS area
        self.primary_is_area[inj_name] = self.tcd_raw_data[
            inj_name
        ][prim_is_name]

    def read_sec_IS_data(self, inj_name):
        '''Read all the necessary information of the secondary Internal Standard
        for the FID channel

        Parameters
        ----------
        self : yield_data instance
        inj_name : String
            Name of the injection tag
        '''
        # Retrieve the name of the secondary IS
        sec_is_name = self.datapoints.loc[
            inj_name,
            'SecondaryInternalStandard_Name'
        ]
        # Read secondaryIS area TCD
        self.sec_is_area_tcd[inj_name] = self.tcd_raw_data[
            inj_name
        ][sec_is_name]
        # Read secondaryIS area FID
        self.sec_is_area_fid[inj_name] = self.fid_raw_data[
            inj_name
        ][sec_is_name]
        # Read the response factor of the secondary IS (FID)
        self.sec_is_response[inj_name] = self.CFs['FID'].loc[
            sec_is_name,
            'Response'
        ]
        # Read the response factor of the secondary IS (TCD)
        self.sec_is_response_tcd[inj_name] = self.CFs['TCD'].loc[
            sec_is_name,
            'Response'
        ]

    def read_tert_IS_data(self, inj_name):
        '''Read all the necessary information of the tertiary Internal
        Standard for the GCxGC channel

        Parameters
        ----------
        self : yield_data instance
        inj_name : String
            Name of the injection tag

        Returns
        -------
        tert_in_tcd : boolean
            True if the reference IS is in the TCD channel
        '''
        tert_is_name = self.datapoints.loc[
            inj_name,
            'TertiaryInternalStandard_Name'
        ]
        ######################################
        # Enhancement
        ######################################
        # Change to InChI before comparison
        ######################################
        tert_in_tcd = True

        if tert_is_name == 'methane':
            self.tert_is_area[inj_name] = self.tcd_raw_data[
                inj_name
            ]['CH4']
            self.tert_is_response_tcd[inj_name] = self.CFs['TCD'].loc[
                'CH4',
                'Response'
            ]
        elif tert_is_name == 'CH4':
            self.tert_is_area[inj_name] = self.tcd_raw_data[
                inj_name
            ][tert_is_name]
            self.tert_is_response_tcd[inj_name] = self.CFs['TCD'].loc[
                'CH4',
                'Response'
            ]
            tert_is_name = 'methane'
        elif tert_is_name == 'ethene' or tert_is_name == 'ethylene':
            self.tert_is_area[inj_name] = self.tcd_raw_data[
                inj_name
            ]['C2H4']
            self.tert_is_response_tcd[inj_name] = self.CFs['TCD'].loc[
                'C2H4',
                'Response'
            ]
        elif tert_is_name == 'propylene':
            self.tert_is_area[inj_name] = self.fid_raw_data[
                inj_name
            ]['C3H6']
            self.tert_is_response_tcd[inj_name] = self.CFs['FID'].loc[
                'C3H6',
                'Response'
            ]
            tert_in_tcd = False
        elif tert_is_name == 'propene':
            self.tert_is_area[inj_name] = self.fid_raw_data[
                inj_name
            ]['C3H6']
            self.tert_is_response_tcd[inj_name] = self.CFs['FID'].loc[
                'C3H6',
                'Response'
            ]
            tert_in_tcd = False
        elif tert_is_name == '1,3-butadiene':
            self.tert_is_area[inj_name] = self.fid_raw_data[
                inj_name
            ]['1,3-C4H6']
            self.tert_is_response_tcd[inj_name] = self.CFs['FID'].loc[
                '1,3-C4H6',
                'Response'
            ]
            tert_in_tcd = False

        # Read tertiaryIS area gcxgc
        self.tert_is_area_gcxgc[inj_name] = self.gcxgc_fid_raw_data[
            inj_name
        ][tert_is_name]
        # Read tert is response factor
        self.tert_is_response[inj_name] = self.gcxgc_lib.loc[
            tert_is_name,
            'Response Factor'
        ]

        return tert_in_tcd

    def read_quat_IS_data(self, inj_name):
        '''Read all the necessary information of the Quaternary Internal
        Standard for the LOA channel

        Parameters
        ----------
        self : yield_data instance
        inj_name : String
            Name of the injection tag

        Returns
        -------
        quat_is_in_tcd : boolean
            True if the reference IS is in the TCD channel
        '''
        quat_is_name = self.datapoints.loc[
            inj_name,
            'QuaternaryInternalStandard_Name'
        ]
        self.quat_is_area_loa[inj_name] = self.loa_raw_data[
            inj_name
        ][quat_is_name]
        self.quat_is_response[inj_name] = self.CFs['LOA'].loc[
            quat_is_name,
            'Response'
        ]
        # Check which reference detector should be chosen
        quat_is_in_tcd = quat_is_name in self.CFs['TCD'].index

        # TCD as reference even though name says FID
        if quat_is_in_tcd:
            self.quat_is_area_fid[inj_name] = self.tcd_raw_data[
                inj_name
            ][quat_is_name]
            self.quat_is_response_fid[inj_name] = self.CFs['TCD'].loc[
                quat_is_name,
                'Response'
            ]
        # FID as reference
        else:
            self.quat_is_area_fid[inj_name] = self.fid_raw_data[
                inj_name
            ][quat_is_name]
            self.quat_is_response_fid[inj_name] = self.CFs['FID'].loc[
                quat_is_name,
                'Response'
            ]

        return quat_is_in_tcd

    def read_IS_data(self):
        '''Read all the necessary information of the available Internal Standards

        Parameters
        ----------
        self : yield_data instance
        '''
        # Iterate over all injections
        import pandas as pd

        for inj_name in self.datapoints.index:
            # Read primary IS data (RGA TCD)
            self.read_prim_IS_data(inj_name)
            # Read secondary IS data (RGA FID)
            self.read_sec_IS_data(inj_name)
            # Read tertiary IS data (GCxGC FID)
            if 'TertiaryInternalStandard_Name' in self.datapoints.columns:
                self.tert_in_tcd[inj_name] = self.read_tert_IS_data(inj_name)
            # Read quaternary IS data (LOA)
            if 'QuaternaryInternalStandard_Name' in self.datapoints.columns:
                self.quat_in_tcd[inj_name] = self.read_quat_IS_data(inj_name)

        # Transform IS data to Pandas Series
        self.primary_is_response = pd.Series(self.primary_is_response)
        self.primary_is_area = pd.Series(self.primary_is_area)

        self.sec_is_area_tcd = pd.Series(self.sec_is_area_tcd)
        self.sec_is_area_fid = pd.Series(self.sec_is_area_fid)
        self.sec_is_response = pd.Series(self.sec_is_response)
        self.sec_is_response_tcd = pd.Series(self.sec_is_response_tcd)

        self.tert_is_area = pd.Series(self.tert_is_area)
        self.tert_is_area_gcxgc = pd.Series(self.tert_is_area_gcxgc)
        self.tert_is_response = pd.Series(self.tert_is_response)
        self.tert_is_response_tcd = pd.Series(self.tert_is_response_tcd)

        self.quat_is_area_loa = pd.Series(self.quat_is_area_loa)
        self.quat_is_area_fid = pd.Series(self.quat_is_area_fid)
        self.quat_is_response = pd.Series(self.quat_is_response)
        self.quat_is_response_fid = pd.Series(self.quat_is_response_fid)


class feed():

    def __init__(self, feed_names, feed_rates, feed_names_iupac):
        '''
        Constructor
        '''

        # Variables that compose a feed
        self.feed_names = feed_names
        self.feed_names_iupac = dict(zip(self.feed_names, feed_names_iupac))
        self.feed_rates = feed_rates
        # Elemental inflows per feed as a percentage
        self.inflow = {}

    def get_flow_total(self):
        ''' Returns the sum of all HCs '''
        return self.feed_rates.sum(axis=1)

    def get_flows(self):
        ''' Returns the feed rate of all HCs independently'''
        return self.feed_rates

    def get_fractions(self):
        ''' Returns the fractions for the feeds '''
        flows = self.get_flows().copy()
        total_flow = self.get_flow_total().copy()
        returner = flows.divide(total_flow, axis=0) * 100
        return returner

    def get_names_inchi(self):
        ''' Returns the names of the feeds '''
        return self.feed_names

    def get_names_smiles(self):
        ''' Returns the names of the feeds in Smiles format '''

        smiles = []
        for name in self.get_names_inchi():
            smiles.append(inchi_2_smiles(name))

        return smiles

    def get_total_elemental_flow(self, element):
        ''' Returns the total elemental inflow percentage wise'''
        elemental_inflow = self.inflow[element].transpose()
        if len(elemental_inflow.columns) == 1:
            returner = elemental_inflow
        else:
            returner = elemental_inflow.sum(axis=1)

        return returner


def parse_XML_input_file(filename, config_names):
    '''
    Converts Ruben Debruycker's XML input xmlbeans
    to pandas DataFrame: datapoints and pandas Series: configs
    '''

    # imports
    import xml.etree.ElementTree as ET
    import pandas as pd

    # Store general configurations separately
    configs = {config_name: None for config_name in config_names}

    # Store all experimental datapoints
    datapoints = {}

    # Open XML file and parse to elementTree
    tree = ET.parse(filename)
    root = tree.getroot()
    # Iterate over all
    for child in root:
        # Store all configs
        if child.tag in configs:
            configs[child.tag] = child.text
        # Store all datapoints
        if child.tag == 'Datapoint':
            # Iterate over everything in datapoint
            for grandchild in child:

                # Find Datapoints name
                if grandchild.tag == 'Name':
                    name = grandchild.text
                    datapoint_info = {}

                # Filter out all internal standard information first
                elif 'Standard' in grandchild.tag:
                    for grgrandchild in grandchild:
                        datapoint_info['{}_{}'.format(
                            grandchild.tag,
                            grgrandchild.tag
                        )] = grgrandchild.text

                # Filter out all feed information
                elif grandchild.tag == 'Feed':
                    for grgrandchild in grandchild:
                        if grgrandchild.tag == 'Name':
                            feedname = grgrandchild.text
                    datapoint_info[
                        '{}_{}_{}'.format(
                            grandchild.tag,
                            feedname,
                            grgrandchild.tag
                        )
                    ] = grgrandchild.text

                # All other data
                else:
                    datapoint_info[grandchild.tag] = grandchild.text

                datapoints[name] = datapoint_info

    conditions = pd.DataFrame(datapoints).transpose().convert_objects()
    configs = pd.Series(configs)

    return conditions, configs


def gcxgc_remover(x):
    ''' Can be used to remove path included in GCxGC filename '''
    import re
    x = re.sub(r'^.*\GCxGC\(.*)$', r'\1', x)
    return x


def rga_remover(x):
    ''' Can be used to remove path included in RGA filename '''
    import re
    x = re.sub(r'^.*\RGA\(.*)$', r'\1', x)
    return x


def reading_input_conditions(filename):
    ''' Reading input conditions for GC processing '''
    ######################################
    # Enhancement
    ######################################
    # There should be a check whether all
    # necessary tags are included in the
    # input files
    ######################################

    ######################################
    # Enhancement
    ######################################
    # Check if delimited correctly
    # ',' instead of ';'
    ######################################

    import pandas as pd

    # Check which file is given as input xml or csv
    if filename[0].endswith('.csv'):
        # Datapoints
        datapoints = pd.read_csv(filename[0], index_col=0)
        # Clean up columns and indices
        datapoints.columns = [col.strip() for col in datapoints.columns]
        datapoints.index = [ind.strip() for ind in datapoints.index]

        # General configurations
        gen_configs = pd.read_csv(
            filename[1],
            index_col=0,
            header=None,
            squeeze=True
        )
        # Clean up indices
        gen_configs.index = [ind.strip() for ind in gen_configs.index]

    elif filename[0].endswith('.xml'):
        config_names = [
            'PROCESS_LOA',
            'RF_GCxGC',
            'RF-RGA-TCD',
            'RF-RGA-FID',
            'RF_locations',
            'EXPERIMENTAL_name',
            'RGA_location',
            'GCxGC_location',
            'Experimental_filename',
            'LoggingFile_location',
            'Reactor'
         ]

        filename = filename[0]
        datapoints, gen_configs = parse_XML_input_file(filename, config_names)
#         datapoints['GCxGC'] = datapoints['GCxGC'].apply(gcxgc_remover)
#         datapoints[['RGA-FID', 'RGA-TCD-L', 'RGA-TCD-R']] = datapoints[
#             [
#                 'RGA-FID',
#                 'RGA-TCD-L',
#                 'RGA-TCD-R'
#             ]
#         ].applymap(rga_remover)

        # Create input files
        datapoints.to_csv('conditions.csv')
        gen_configs.to_csv('gen_configs.csv')

    gen_configs.dropna(inplace=True)
    gen_configs.index.name = 'sleutels'

    print("\nInput conditions successfully read")
    return datapoints, gen_configs


def store_processing_conditions_hdf(datapoints, gen_configs, h5filename=None):
    '''Stores datapoints and gen_configs to HDF5 file container

    Parameters
    ----------
    datapoints : dct
    gen_configs : dct
    h5filename : str, optional (default None)

    Returns
    -------
    None
    '''
    import pandas as pd
    import os

    sublist = list(datapoints.groupby('Date').groups.keys())
    sublist.sort()
    date_to_append = sublist[0].replace('/', '_')
    experimental_name = gen_configs['EXPERIMENTAL_name']

    if h5filename == None:
        h5filename = '{}_{}.h5'.format(experimental_name, date_to_append)
    else:
        h5filename = '{}.h5'.format(h5filename)

    # If HDF5 file already exists delete folders
    if os.path.isfile(h5filename):
        store = pd.HDFStore(h5filename)
        if 'processing/conditions' in store:
            store.remove('processing/conditions')
        if 'processing/gen_configs' in store:
            store.remove('processing/gen_configs')

    # Store datapoints
    storing_2_HDF5_store(
        datapoints,
        h5filename,
        'processing/conditions',
    )

    # Store general configurations
    storing_2_HDF5_store(
        gen_configs.reset_index(),
        h5filename,
        'processing/gen_configs',
        data_column=False
    )

    # Close HDF5 file
    try:
        store.close()
    except:
        pass


def reading_agilent_data(filelocation, filename):
    '''
    Reads the csv summary created by the OpenLab software

    Parameters
    ----------
    filelocation : The path location where the file can be found
    filename : The name of the file

    Returns
    -------
    data : DataFrame that contains the raw data
        columns : Area_perc, Ret_Time, Area
        index : components (human readable name)

    Examples
    --------
    reading_agilent_data(
    r'unit_tests\testFiles_dataProcessing',
    r'fid_agilent.csv'
    )
    >>> assert (check pandas assert equal...)

    reading_hyperchrom_data(
    r'unit_tests\testFiles_dataProcessing',
    r'fid_agilent_empty.csv'
    )
    >>> None

    '''
    import pandas as pd

    # Read the csv file with pandas
    print(filename)
    df = pd.read_csv(
        '{}/{}'.format(filelocation, filename),
        header=None,
        index_col=0,
        names=['Ret_Time', 'Type', 'Width', 'Area', 'Area_perc', 'Name'],
        encoding='utf-16',
        delimiter=','
    )

    # Drop unnecessary columns
    df.drop(['Type', 'Width'], axis=1, inplace=True)

    # Drop rows with zeros in 'Area' column
    df.query('Area != 0', inplace=True)

    # Drop rows with no name
    df.query('Name != "?"', inplace=True)

    # Drop rows with C6+ as name
    df.query('Name != "C6+"', inplace=True)

    # Set 'Name' column as the new index
    df.set_index('Name', inplace=True, drop=True)

    # Clean up the indices
    df.index = [ind.strip() for ind in df.index]

    # Put the columns in the same order as Hyperchrom
    df = df[['Area_perc', 'Ret_Time', 'Area']]

    # Convert all columns to numeric objects
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])

    return df


def reading_hyperchrom_data(filelocation, filename):
    '''
    Reads the excel summary created by the HyperChrom software

    Parameters
    ----------
    filelocation : The path location where the file can be found
    filename : The name of the file

    Returns
    -------
    data : DataFrame that contains the raw data
        columns :
        index : components (human readable name)

    Examples
    --------
    reading_hyperchrom_data(
    r'unit_tests\testFiles_dataProcessing',
    r'propane_800_FID_Run1.XLS'
    )
    >>> assert (check pandas assert equal...)

    reading_hyperchrom_data(
    r'unit_tests\testFiles_dataProcessing',
    r'empty_chrom.XLS'
    )
    >>> None

    '''
    import pandas as pd
    import xlrd

    # Check if the files are available before opening
    if not check_file_available(filelocation + '/' + filename):
        popup_error(
            '{} could not be found in {}'.format(filename, filelocation)
        )
        sys.exit(0)
    
    # Open the excel file
    wb = xlrd.open_workbook('{}/{}'.format(filelocation, filename))

    # Find number of rows to skip and the column number of the index
    sheet = wb.sheet_by_index(0)
    # Get total number of columns
    cols = sheet.ncols
    # Initialize
    rows_to_skip = 0
    # col_ind = 0
    # Go over all cells
    for i in range(cols):
        elements = [str(el).strip() for el in sheet.col_values(i)]
        if 'Peak Number #' in elements:
            rows_to_skip = elements.index('Peak Number #')
            # col_ind = i
        elif 'Component Name' in elements:
            rows_to_skip = elements.index('Component Name')
            # col_ind = i

    # Read the excel file with pandas
    df = pd.read_excel(
        '{}/{}'.format(filelocation, filename),
        skiprows=rows_to_skip
    )

    # Clean up the column names
    df.columns = [
        col.strip().replace(
            '#', ''
        ).replace(' ', '_').replace('%', 'perc').replace('.', '_')
        for col in df.columns
    ]

    # Drop rows with empty lines
    df.dropna(axis=0, how='any', inplace=True)

    # Set Peak_number/Component_Name column as the new index
    # and give the index a name
    try:
        df.set_index('Peak_Number_', inplace=True)
    except:
        df.set_index('Component_Name', inplace=True)
    df.index.name = 'Name'

    # Clean up the indices
    df.index = [ind.strip() for ind in df.index]

    # Drop column BC
    try:
        df.drop('BC', axis=1, inplace=True)
    except ValueError:
        df.drop('Ident__Number', axis=1, inplace=True)

    # Convert all columns to numeric objects
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])

    return df


def reading_cal_factor(file_location, file_name):
    '''
    Reads the excel files with the calibration factors for the RGA and LOA

    Parameters
    ----------
    file_location : String that contains the folder path location
    file_name : String that contains the name of the file

    Returns
    -------
    CF : DataFrame with calibration factors

    Examples
    --------
    reading_cal_factor(
    r'unit_tests/testFiles_dataProcessing',
    r'LOA_Response_MTHF_12082016.xlsx'
    )
    >>> assert equal
    '''
    import os
    import pandas as pd
    import numpy as np

    # File location
    file_location_cf = os.path.join(file_location, file_name)

    # Read the calibration file
    CF = pd.read_excel(file_location_cf, index_col=0)

    # Clean up the indices and columns
    CF.index = [ind.strip() for ind in CF.index]
    CF.columns = [col.strip() for col in CF.columns]

    # Create an extra column with the molecular weight
    CF['MW'] = np.dot(
        CF[['C', 'H', 'O', 'N', 'S']],
        [12.011, 1.008, 15.999, 14.007, 32.065]
    )

    return CF


def reading_reactor_config(file_location, file_name):
    '''Reads the reactor configuration file for reading the loggingfiles

    Extended description

    Parameters
    ----------
    file_location : String
        Pathname where the reactor configuration file can be found
    file_name : String
        Name of the reactor configuration file

    Returns
    -------
    reactor_config : dict
        loggingFileNames : ndarray
            A list of the loggingfile names
        decimal : String
            The character that is used for a decimal
        geometry : DataFrame
            reactor geometry
        log_name_temperature : String
            Name of the loggingfile that holds the temperatures
        log_name_pressure : String
            Name of the loggingfile that holds the pressures
        divide_by : int
            Badly logged BSSC quartz divide by 10
        axial_temperatures : DataFrame
            Loggingfile tags for temperature with axial position
        axial_pressures : DataFrame
            Loggingfile tags for pressure with axial position
        names_temperature : ndarray
            A list of the temperature tags
        names_pressure : ndarray
            A list of the pressure tags

    Examples
    --------
    >>> reading_reactor_config(r'')
    '''
    import pandas as pd
    import os

    # Initialisation
    reactor_config = {}

    full_path_name = os.path.join(file_location, file_name)

    reactor_config['loggingFileNames'] = pd.read_excel(
        full_path_name,
        sheet_name='loggingFileNames',
        index_col=0,
        header=None
    ).index.values

    config = pd.read_excel(
        full_path_name,
        sheet_name='config',
        index_col=0,
    )

    # Read in the reactor geometry data
    reactor_config['geometry'] = config.loc[
        ['length (m)', 'diameter (m)', 'thickness (m)'],
        'values'
    ].to_frame()

    # Badly logged temperatures BSSC
    try:
        reactor_config['divide_by'] = config.loc[
            'divide_by',
            'values'
        ]
    except:
        reactor_config['divide_by'] = 1
    # Are the pressures logging in barg or bara
    try:
        reactor_config['barg'] = config.loc[
            'barg',
            'values'
        ]
    except:
        reactor_config['barg'] = 0

    # Get the name of the loggingFile with temperatures
    reactor_config['log_name_temperature'] = config.loc[
        'Temperatures',
        'values'
    ]
    # Get the name of the loggingFile with pressures
    reactor_config['log_name_pressure'] = config.loc[
        'Pressures',
        'values'
    ]

    # Change the name of the column of the geometry df
    reactor_config['geometry'].columns = [0]
    # Change dtype to float or int
    reactor_config['geometry'][0] = pd.to_numeric(
        reactor_config['geometry'][0],
        errors='coerce'
    )

    reactor_config['decimal'] = config.loc['decimal', 'values']

    sheetnames = ['axial_temperatures', 'axial_pressures']
    names = ['names_temperature', 'names_pressure']

    for sheetname, name in zip(sheetnames, names):
        reactor_config[sheetname] = pd.read_excel(
            full_path_name,
            sheet_name=sheetname,
            index_col=0
        )
        reactor_config[name] = reactor_config[sheetname].index.values

    return reactor_config


def reading_storing_RGA_data(h5filename):
    '''Reading all chromatographic RGA data and store it in HDF5 file

    Parameters
    ----------
    h5filename : String
        HDF5 file name without ".h5"

    Returns
    -------
    rga_channels : MultiIndexed DataFrame
        They are also stored to the HDF5 file
    '''

    import numpy as np
    import os
    import pandas as pd
    import sys

    # Open the HDF5 store
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Get general configuration file and conditions file
    gen_configs = store['processing/gen_configs'].set_index('sleutels')
    conditions = store['processing/conditions']
    # Get some values from the general configurations
    reactor = gen_configs.loc['Reactor'].values[0]
    rga_gen_location = r'{}'.format(gen_configs.loc['RGA_location'].values[0])
    xls_csv = '.XLS' if reactor != 'IHEB' else '.csv'
    # Get the folder name for the raw chromatographic data
    exp_filename = gen_configs.loc["Experimental_filename"].values[0]
    rga_location = os.path.join(rga_gen_location, exp_filename)
    # Get the names of the RGA files
    fid_names = conditions['RGA-FID'].values
    tcdl_names = conditions['RGA-TCD-L'].values
    if reactor not in ['JSR']:
        tcdr_names = conditions['RGA-TCD-R'].values
        rga_channel_names = ['channel_A', 'channel_B', 'channel_C']
        chrom_names_in_conditions = np.concatenate(
            (fid_names, tcdl_names, tcdr_names)
        ).astype(str)
    else:
        rga_channel_names = ['channel_A', 'channel_B']
        chrom_names_in_conditions = np.concatenate(
            (fid_names, tcdl_names)
        ).astype(str)

    # If the RGA of the IHEB needs to be processed. The structure of the files
    # is different. This module transforms the structure to the regular one.
    if reactor in ['IHEB']:
        restructure_agilent_RGA_data(rga_location)

    # Check if any RGA files are available in the HDF5 store
    # and if they correspond with the ones from the conditions file
    if 'raw/RGA' in store:
        chrom_in_HDF5 = []
        rga_channels = []

        # Iterate over the available channels and store the results
        for rga_channel_name in rga_channel_names:
            rga_channel = store['raw/RGA/{}'.format(rga_channel_name)]
            rga_channels.append(rga_channel)
            chrom_in_HDF5.append(
                rga_channel.index.get_level_values(0).unique()
            )
        # Get a complete ndarray of the RGA files present
        chrom_in_HDF5_flat = np.concatenate(chrom_in_HDF5).astype(str)

        arrays_equal = np.array_equal(
            np.sort(np.char.replace(chrom_names_in_conditions, '.XLS', '')),
            np.sort(chrom_in_HDF5_flat)
        )

        if not arrays_equal:
            popup_error(
                "The RGA chromatograms in the HDF5 file are not the same "
                "as the ones in the conditions file.\n"
                "Try processing the raw chromatograms again. "
                "by deleting the HDF5 file."
            )
            store.close()
            sys.exit()

        print('\nRGA data checked\n')

    # The RGA files need to be read and processed
    else:
        # Initialize RGA results
        rga = {}
        # Iterate over all filenames in the conditions file
        for filename in chrom_names_in_conditions:
            if reactor not in ['IHEB']:
                rga_df = reading_hyperchrom_data(rga_location, filename)
            else:
                rga_df = reading_agilent_data(rga_location, filename)

            # Make sure the file is not empty
            if rga_df.empty:
                popup_error("Chromatogram: '{}' is empty".format(filename))
                store.close()
                sys.exit()
            else:
                rga[filename] = rga_df

        # Create MultiIndexed DataFrames (sorted for slice selecting)
        rga_chA_dct = {
            k.replace(xls_csv, ''): i
            for k, i
            in rga.items()
            if k in tcdl_names
        }
        rga_chA = pd.concat(rga_chA_dct).sort_index()
        rga_chC_dct = {
            k.replace(xls_csv, ''): i
            for k, i
            in rga.items()
            if k in fid_names
        }
        rga_chC = pd.concat(rga_chC_dct).sort_index()

        # The JSR Reactor has only 2 channels for the RGA (no Channel C)
        if reactor not in ['JSR']:
            rga_chB_dct = {
                k.replace(xls_csv, ''): i
                for k, i
                in rga.items()
                if k in tcdr_names
            }
            rga_chB = pd.concat(rga_chB_dct).sort_index()

            # Combine all RGA results
            rga_channels = [rga_chA, rga_chB, rga_chC]
            rga_channelNames = ['channel_A', 'channel_B', 'channel_C']
        else:
            rga_channels = [rga_chA, rga_chC]
            rga_channelNames = ['channel_A', 'channel_B']

        # Iterate over all channels and store them to the HDF5 file
        for channel, channelName in zip(rga_channels, rga_channelNames):
            # Save to HDF5 file
            storing_2_HDF5_store(
                channel,
                store,
                'raw/RGA/{}'.format(channelName),
            )

        print("\nRGA data stored to {}.h5".format(h5filename))

    # Close the HDF5 file
    store.close()

    return rga_channels


def reading_gcxgc_blob_table(
    gcxgc_location,
    gcxgc_filename,
    store,
):
    '''Reads GC Image blob tables and returns a DataFrame

    Parameters
    ----------
    gcxgc_location : String
        Folder name in which the GCxGC are stored
    gcxgc_filename : String
        Name of the GCxGC blob table
    type : String
        csv or xls

    Returns
    -------
    gcxgc : DataFrame
    '''
    import os
    import pandas as pd
    import sys

    # Get file name and path name
    gcxgc_filename_and_path = os.path.join(
        gcxgc_location,
        gcxgc_filename
    )

    # Read GCxGC blob table first check if present
    if not check_file_available(gcxgc_filename_and_path):
        popup_error(
            "{} GCxGC file not found in {}.\n"
            "Program terminated.".format(
                gcxgc_filename,
                gcxgc_location
            )
        )
        store.close()
        sys.exit()

    gcxgc = pd.read_csv(
        '{}'.format(gcxgc_filename_and_path)
    )

    # Reformat column names
    gcxgc.columns = [
        col.strip().replace(' ', '_').replace('(', '').replace(')', '')
        for col in gcxgc.columns
    ]

    # Combine Compound name and group name
    gcxgc['Name'] = gcxgc['Compound_Name'].combine_first(
        gcxgc['Group_Name']
    )

    # Disregard not named blobs
    gcxgc = gcxgc[pd.notnull(gcxgc["Name"])]

    # Reformat names don't sum same named blobs,
    # because SCD needs to be subtracted still
    # I don't think SCD will be completely combined.
    gcxgc['Name'] = [name.strip() for name in gcxgc['Name']]

    return gcxgc


def reading_storing_LOA_data(h5filename, datapoints, gen_configs):
    '''
    Reading all LOA files and stores them in the HDF5 file

    Parameters
    ----------
    h5filename : String
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection
    gen_configs : DataFrame

    Returns
    -------
    loa_channel : MultiIndexed DataFrame with all loa data

    Examples
    --------
    reading_storing_LOA_data(h5filename, datapoints, gen_configs)
    >>> assert df equal multiIndex

    '''

    import os
    import pandas as pd
    import sys

    loa_gen_location = r'{}'.format(gen_configs['LOA_location'])

    exp_filename = gen_configs['Experimental_filename']
    loa_location = os.path.join(loa_gen_location, exp_filename)

    _, _, loa_filenames = next(os.walk(loa_location), (None, None, []))
    loa = {}

    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Check if the raw LOA data is already stored in the HDF5 file
    if 'raw/LOA' not in store:

        # Find all excel files in the given folder
        loa_filenames = [
            loa_file
            for loa_file
            in loa_filenames
            if loa_file.endswith('.XLS')
        ]

        # Check if there are LOA files present in the given folder.
        if not len(loa_filenames):
            print('No LOA files found.')

        # Not all of them present
        else:
            # Iterate over all files
            for filename in loa_filenames:
                loa_df = reading_hyperchrom_data(loa_location, filename)

                # Don't save empty DataFrames
                if not loa_df.empty:
                    loa[filename] = loa_df

        # Create MultiIndexed DataFrames (sorted for slice selecting)
        loa = pd.concat(
            {
                k.replace('.XLS', ''): i
                for k, i
                in loa.items()
            }
        ).sort_index()

        # Save to hdf5 file
        storing_2_HDF5_store(
            loa,
            store,
            'raw/LOA',
        )

        print("\nLOA data stored to {}.h5".format(h5filename))

    # The key 'raw/LOA' exists -> LOA data is processed -> check if correct
    else:
        # Check if all LOA files are present in the HDF5 file
        loa = store['raw/LOA']
        chrom_in_HDF5 = loa.index.get_level_values(0).unique().tolist()

        chrom_in_datapoints = list(
            store[
                'processing/conditions'
                ][['LOA']].as_matrix().flatten()
        )

        # Check if all chromatograms can be found
        chrom_not_found = [
            chrom
            for chrom
            in chrom_in_datapoints
            if chrom.replace('.XLS', '') not in chrom_in_HDF5
        ]
        if not len(chrom_not_found) == 0:
            popup_error(
                "Following chromatograms\n{}\nare not found"
                " in folder:\n{}\n".format(
                    chrom_not_found,
                    loa_location
                )
            )
            sys.exit()

        print('\nLOA data checked\n')

    # Close the HDF5 file
    store.close()

    return loa


def reading_storing_GCxGC_scd_data(h5filename, datapoints, gen_configs):
    ''' Reading all GCxGC SCD blob tables and stores them in the HDF5 file '''
    ######################################
    # SCD_NCD_detector
    ######################################
    # This should remain the same I think
    ######################################
    reading_storing_GCxGC_fid_data(
        h5filename,
        datapoints,
        gen_configs,
        detector_switch='SCD'
    )


def reading_storing_GCxGC_ncd_data(h5filename, datapoints, gen_configs):
    ''' Reading all GCxGC NCD blob tables and stores them in the HDF5 file '''
    ######################################
    # SCD_NCD_detector
    ######################################
    # This should remain the same I think
    # Check in 'reading_storing_GCxGC_fid_data'
    # if gcxgc at the end looks similar
    # for FID and SCD blob table
    ######################################
    reading_storing_GCxGC_fid_data(
        h5filename,
        datapoints,
        gen_configs,
        detector_switch='NCD'
    )


def reading_storing_GCxGC_fid_data(
    h5filename,
    datapoints,
    gen_configs,
    detector_switch='FID'
):
    '''Reading all GCxGC FID blob tables and store them in HDF5 file

    Parameters
    ----------
    h5filename : String
        HDF5 file name without ".h5"
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection
    gen_configs : dct
    detector_switch : String
        FID, SCD, or NCD

    Returns
    -------
    None
    '''
    import os
    import pandas as pd
    import sys

    # Get reactor switch (JSR has 1D C5+)
    reactor = gen_configs['Reactor']

    # Get proper file location
    gcxgc_gen_location = gen_configs[
        'GCxGC_{}_location'.format(detector_switch)
    ]
    exp_filename = gen_configs['Experimental_filename']
    gcxgc_location = os.path.join(gcxgc_gen_location, exp_filename)

    # Open HDF5 file
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Check if GCxGC data is present in HDF5 file
    gcxgc_in_store = 'raw/GCxGC_{}'.format(detector_switch) in store

    # Initiate dictionary to store the files together
    # Read existing data
    if gcxgc_in_store:
        df = pd.read_hdf(store, key='raw/GCxGC_{}'.format(detector_switch))
        grouped_df = df.groupby(level=0)
        combined_gcxgc_data = {
            name: grouped_df.get_group(name).xs(name)
            for name
            in grouped_df.groups
        }
        # Delete GCxGC data from store to append an entirely new one
        store.remove('raw/GCxGC_{}'.format(detector_switch))
    else:
        combined_gcxgc_data = {}

    # Iterate over all GCxGC blob tables that are given in datapoints.
    for gcxgc_filename in datapoints.groupby('GCxGC').groups.keys():

        # Get a shorter filename
        gcxgc_filename_no_csv = gcxgc_filename.replace('.csv', '').replace(
            '.XLS', ''
        )
        print(gcxgc_filename)

        # Create a flag to check if single GCxGC chrom in HDF5 store
        gcxgc_chrom_in_store = False
        if gcxgc_in_store:
            gcxgc_chrom_in_store = gcxgc_filename_no_csv in combined_gcxgc_data

        # Check if GCxGC file is already present in HDF5 file
        # If not, read it in, clean it and store it to HDF5 file
        if not gcxgc_chrom_in_store:

            if reactor == 'JSR':
                gcxgc = reading_hyperchrom_data(
                    gcxgc_location,
                    gcxgc_filename,
                )
            else:
                gcxgc = reading_gcxgc_blob_table(
                    gcxgc_location,
                    gcxgc_filename,
                    store,
                )

            # Store files to dict
            combined_gcxgc_data[gcxgc_filename_no_csv] = gcxgc

        # Transform dict to MultiIndexed DataFrame via concatenating
        combined_gcxgc_data_df = pd.concat(combined_gcxgc_data)

        # Store files to HDF5 file
        storing_2_HDF5_store(
            combined_gcxgc_data_df,
            store,
            'raw/GCxGC_{}'.format(detector_switch),
        )

    # Close HDF5 file
    store.close()
    print("\nGCxGC data stored to {}.h5".format(h5filename))

    return None


def reading_storing_cal_factors(h5filename, gen_configs):
    '''
    Reads the calibration factors for the RGA and stores them to the HDF5 file

    Parameters
    ----------
    h5filename : String, HDF5 file name
    gen_configs : DataFrame with general configurations

    Returns
    -------
    CFs : dict
        values are DataFrames that contain the calibration factors for every
        detector
    '''

    import os
    import pandas as pd
    import sys

    # Check if LOA has to be processed as well
    if gen_configs['Process_LOA'] == 'TRUE':
        names = ['RF-RGA-TCD', 'RF-RGA-FID', 'RF-LOA']
    else:
        names = ['RF-RGA-TCD', 'RF-RGA-FID']

    CFs = {}

    # Open the HDF5 store
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Replace existing calibration factors with new ones
    if 'processing/CF' in store:
        del store['processing/CF']

    # Iterate over all detectors
    for name in names:

        # File locations
        file_location = gen_configs['RF_locations']
        file_name = gen_configs[name]

        file_location_cf = os.path.join(
            file_location,
            file_name
            )

        # File calibration factors
        if not check_file_available(file_location_cf):
            popup_error(
                "Calibration factors not found: {} in {}\n"
                "Program terminated.".format(
                    file_name,
                    file_location
                )
            )
            store.close()
            sys.exit()

        # Read the CF from the excel file
        CF = reading_cal_factor(file_location, file_name)

        # Store them to the dict
        CFs[name.split('-')[-1]] = CF

        # Store the calibration factors to the HDF5 file
        storing_2_HDF5_store(
            CF,
            store,
            'processing/CF/{}'.format(name.split('-')[-1]),
        )

    # Close HDF5 file
    store.close()

    print(
        "Calibration factors successfully read and "
        "stored to {}.h5".format(h5filename)
    )
    return CFs


def reading_loggingfile(
    log_location,
    log_file_full_name,
    decimal='.',
    old=False,
    jsr=False
):
    '''Reads in the loggingfile from the log_location and returns the DataFrame

    More description

    Parameters
    ----------
    log_location : String
        Location of the logging file
    log_file_full_name : String
        Full name of the logging_file
    decimal : String
        decimal character that is used in the loggingfile
    old : boolean, optional
        The structure of the loggingfile is different for
        the old BSSC reactor
    jsr : boolean, optional
        The structure of the loggingfile is different for
        the JSR reactor

    Returns
    -------
    data : DataFrame
        Cleaned up DataFrame with the content of that loggingfile

    Examples
    --------
    >>> reading_loggingfiles(r'', r'')

    '''
    import pandas as pd
    import os

    # Get full path name
    full_path_name = os.path.join(log_location, log_file_full_name)

    # Read in the csv loggingfile
    if old:  # older formatted loggingfiles
        data = pd.read_csv(
            full_path_name,
            index_col=24,
            sep='\t',
            parse_dates=True,
            dayfirst=True
            )
    elif jsr:  # JSR formatted loggingfiles
        data = pd.read_table(
            full_path_name,
            index_col=0
        )
        data.dropna(inplace=True, axis=1)
    else:  # PLC logging files
        data = pd.read_csv(
            full_path_name,
            index_col=0,
            sep=';',
            decimal=decimal,
            parse_dates={'DateTime': [0, 1]},
            dayfirst=True
            )

    # Clean up the column and index names
    data.columns = [col.strip() for col in data.columns]
    data.index.name = 'DateTime'

    return data


def reading_storing_loggingfiles(datapoints, gen_configs, h5filename):
    '''Reading and storing the loggingfiles

    Reading and storing the logging files from the days the
    experiment was performed

    Parameters
    ----------
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection
    gen_configs : Series
        general configuration file
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None

    Examples
    --------

    '''
    import os
    import pandas as pd
    import re
    import sys

    dates = datapoints['Date'].unique()
    loggingFile_location = gen_configs['LoggingFile_location']
    reactor = gen_configs['Reactor']

    # Open HDF5 file
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Reactor configuration
    try:
        reactor_config = reading_reactor_config(
            r'N:\Project\918-Experimental set-ups\Reactor_configurations',
            r'{}.xlsx'.format(reactor)
        )
    except FileNotFoundError:  # offline version of the software
        try:
            reactor_config = reading_reactor_config(
                gen_configs['reactor_config_loc'],
                r'{}.xlsx'.format(reactor)
            )
        except FileNotFoundError:
            sys.exit()

    # Initialisations
    old = False
    jsr = False

    # Get the loggingFileNames for the appropriate setup
    decimal = reactor_config['decimal']
    loggingFileNames = reactor_config['loggingFileNames']

    if reactor == 'BSSC_metal_old':
        # BSSC before 2015
        # Older logging files need to be processed in a different way
        # Maybe better in submodule
        old = True
    elif reactor == 'JSR':
        jsr = True

    # Check if all loggingfiles are available
    for date in dates:
        # Get the DD/MM/YYYY numbers
        date_split = date.split('/')

        # Get loggingFileNames for old reactor
        # Search in loggingfile folder for loggingfiles that start with
        # the experimental date
        if old:
            date_old = '{}{}{}'.format(
                date_split[2],
                date_split[1],
                date_split[0],
            )
            loggingFileNames = [
                file
                for file
                in os.listdir(loggingFile_location)
                if
                (os.path.isfile(os.path.join(loggingFile_location, file)))
                and
                (re.sub(r'.*([0-9]{8}).*', r'\1', file) in date_old)
            ]
        elif jsr:
            loggingFileNames = [
                '{}.xls'.format(
                    gen_configs['Experimental_filename']
                )
            ]

        # Iterate over all loggingFileNames
        for loggingFileName in loggingFileNames:

            if not old and not jsr:
                # Get the full file name
                logging_file_full_name = '{}_{}_{}_{}.csv'.format(
                    date_split[2].lstrip('0'),  # Remove starting zeros
                    date_split[1].lstrip('0'),
                    date_split[0].lstrip('0'),
                    loggingFileName
                )
            else:
                logging_file_full_name = loggingFileName

            # Get the full name: path and file name
            full_file_path_name = os.path.join(
                loggingFile_location,
                logging_file_full_name
            )

            # Check if file is available
            if not check_file_available(full_file_path_name):
                popup_error(
                    "Loggingfile {} not found in {}.\n"
                    "Program terminated.".format(
                        logging_file_full_name,
                        loggingFile_location
                    )
                )
                store.close()
                sys.exit()

            # Change name for storing the data if old BSSC
            if old:
                loggingFileName = 'log_temperatures'
            elif jsr:
                loggingFileName = 'log_file_JSR'

            # Check if loggingfile is not already
            # appended to the HDF5 file
            log_key = 'raw/logging/{}_{}'.format(
                loggingFileName,
                date.replace('/', '_')
                )

            if log_key not in store:
                # Read in the loggingfile
                log_file = reading_loggingfile(
                    loggingFile_location,
                    logging_file_full_name,
                    decimal,
                    old,
                    jsr
                    )

                # Store it to HDF5 store
                storing_2_HDF5_store(
                    log_file,
                    store,
                    log_key
                    )

    # Close HDF5 file
    store.close()

    print(
        "\nLoggingFiles successfully appended "
        "to HDF5 file: {}.h5".format(h5filename)
    )

    return None


def restructure_agilent_RGA_data(rga_location):
    '''Restructures the RGA data obtained with the ChemStation software
    to be compliant to the structure of the Hyperchrom software

    Parameters
    ----------
    rga_loctation : String
        Full folder path of the RGA data files

    Returns
    -------
    None
        However it restructures and renames the data
    '''
    import os
    import shutil

    # Get a list of the folders/injections in the rga_location folder
    injection_names = [f for f in os.listdir(rga_location) if f.endswith('.D')]

    # Iterate over all folders
    for injection_name in injection_names:
        injection_number = injection_name.replace('.D', '')
        # Iterate over all channels
        for i, j in zip([1, 2, 3], ['C', 'B', 'A']):
            try:
                # Get the source and destination location for copying
                src_sub = os.path.join(rga_location, injection_name)
                src = os.path.join(src_sub, 'REPORT0{}.CSV'.format(i))
                dst = os.path.join(
                    rga_location,
                    'Dfch{}{}.csv'.format(j, injection_number)
                )
                shutil.copy(src, dst)

            # Error handling
            except FileNotFoundError:
                print(
                    "WARNING!! REPORT0{}.CSV not found in {}".format(
                        i,
                        src_sub
                    )
                )

    return None


def fid_tcd_merge(tcd_location, fid_location, rga_location):
    '''Merges the excel chromatograms from the separate FID and TCD channel
    to a combined RGA folder for JSR setup.

    Parameters
    ----------
    tcd_location : String
        Location of TCD files
    fid_location : String
        Location of FID files
    rga_location : String
        Target location for the merged files

    Returns
    -------
    folder_name : String
        Returns the merge folder name

    Examples
    --------
    >>> tcd_location = r"TCD/20170614-Al alloy 1st CC"
    >>> fid_location = r"FID/20170614-Al alloy 1st CC"
    >>> rga_location = r"RGA"
    >>> fid_tcd_merge(tcd_location, fid_location, rga_location)
    r"20170614-Al alloy 1st CC"
    '''
    import glob
    import os
    import shutil
    import sys

    # Make sure TCD and FID base folders have the same name
    base_folder_tcd = os.path.basename(os.path.normpath(tcd_location))
    base_folder_fid = os.path.basename(os.path.normpath(fid_location))

    if base_folder_fid != base_folder_tcd:
        print(
            "TCD and FID folder don't share the same name.\n"
            "Program terminated."
        )
        sys.exit(0)

    # Create the merged folder
    rga_location = os.path.join(rga_location, base_folder_fid)
    try:
        os.makedirs(rga_location)
    except FileExistsError:
        pass

    # Get current working directory (CWD)
    cwd = os.getcwd()
    rga_location = os.path.join(cwd, rga_location)

    # Get all TCD file names
    try:
        os.chdir(tcd_location)
    except FileNotFoundError:
        print(
            "TCD file location not found.\n"
            "Program terminated"
        )
        sys.exit(0)

    tcd_files = glob.glob('*.XLS')
    # Move and replace one-by-one
    for tcd_file in tcd_files:
        print(tcd_file)
        dst = os.path.join(rga_location, tcd_file)
        print(dst)
        shutil.move(tcd_file, dst)

    # Move back to starting working directory
    os.chdir(cwd)

    # Get all FID file names
    try:
        os.chdir(fid_location)
    except FileNotFoundError:
        print(
            "FID file location not found.\n"
            "Program terminated"
        )
        sys.exit(0)
    fid_files = glob.glob('*.XLS')
    # Move and replace one-by-one
    for fid_file in fid_files:
        dst = os.path.join(rga_location, fid_file.replace('DfchA', 'DfchB'))
        shutil.move(fid_file, dst)

    # Move back to starting working directory
    os.chdir(cwd)

    folder_name = base_folder_fid

    return folder_name


def read_gcxgc_lib(gen_configs):
    ''' Read GCxGC library '''

    import os
    import pandas as pd

    print('\nReading GCxGC library')
    gcxgc_filename = gen_configs['RF_GCxGC']

    ######################################
    # Enhancement
    ######################################
    # Should be taken from gen_configs
    ######################################
    if "gcxgc_lib_loc" in gen_configs:
        gcxgc_filelocation = gen_configs["gcxgc_lib_loc"]
    else:
        gcxgc_filelocation = r'N:\Project\918-Experimental set-ups\GC_library'

    gcxgc_lib = pd.read_excel(
        os.path.join(gcxgc_filelocation, gcxgc_filename),
        index_col=0
    )
    gcxgc_lib.index.name = 'Name'
    gcxgc_lib.index = [name.strip() for name in gcxgc_lib.index]
    gcxgc_lib.columns = [name.strip() for name in gcxgc_lib.columns]
    gcxgc_lib = duplicate_remover(gcxgc_lib)

    print('GCxGC library successfully read.\n')
    return gcxgc_lib


def processing_chrom_loa(h5filename, datapoints):
    '''Extracts the area data from the raw chromatograms and returns
    stores the DataFrame to the HDF5 file

    Parameters
    ----------
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None
    '''
    import pandas as pd

    store = pd.HDFStore('{}.h5'.format(h5filename))

    if '/raw/LOA' in store.keys():
        loa = pd.read_hdf(
            store,
            'raw/LOA',
            columns=['Area']
        ).unstack(level=0).fillna(0)

        loa.columns = loa.columns.droplevel(0)

        # Store to the HDF5 store
        storing_2_HDF5_store(
            loa,
            store,
            'processing/raw/LOA',
        )
        print('\nArea data from LOA detector successfully extracted.\n')

        # Extract only the chroms that are in the datapoints and
        # give the correct column names
        loa_datapoints = datapoints['LOA'].map(
            lambda x: x.replace('.XLS', '')
        )
        loa_raw_data = pd.read_hdf(
            store,
            'processing/raw/LOA',
            columns=loa_datapoints
        )

        loa_datapoints = {y: x for x, y in loa_datapoints.items()}
        # Change column names from rga_fid_raw
        loa_raw_data = loa_raw_data.rename(columns=loa_datapoints)

        # store the LOA data consistent with datapoints to the HDF5 store
        storing_2_HDF5_store(
            loa_raw_data,
            store,
            'processing/LOA',
        )

    store.close()
    return None


def reading_chrom_HDF(h5filename, keyname, chromname, blob=False):
    '''Read the chromatographic data for a given key from a given HDF5 store

    Parameters
    ----------
    keyname : String
        The name of the key inside the HDF5 store
    chromname : String
        The name of the chromatogram without suffix
    h5filename : String
        The name of the HDF5 file without ending '.h5'
    blob : boolean
        Flag to distinguish blob table or Hyperchromdata

    Returns
    -------
    gc_data : DataFrame
    '''
    import pandas as pd
    idx = pd.IndexSlice

    # Blob tables are read with Volume and Name columns only
    if blob:
        gc_data = pd.read_hdf(
            '{}.h5'.format(h5filename),
            keyname,
            columns=['Volume', 'Name']
        ).loc[idx[chromname, :]]

    else:
        gc_data = pd.read_hdf(
            '{}.h5'.format(h5filename),
            keyname,
        ).loc[idx[chromname, :]]['Area']

    return gc_data


def processing_chrom_gcxgc_fid_scd_ncd(h5filename, datapoints):
    '''Extracts the volume data from the raw chromatograms and returns
    stores the DataFrame to the HDF5 file

    Parameters
    ----------
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None
    '''
    import pandas as pd
    idx = pd.IndexSlice

    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Read GCxGC raw data from HDF5 file
    if '/raw/GCxGC_FID' in store.keys():

        try:  # GCxGC GCImage blob table
            gcxgc = pd.read_hdf(
                store,
                'raw/GCxGC_FID',
                columns=['Volume', 'Name']
            )

            # Initialise
            combined = {}

            # Iterate over all chromatograms
            for i in gcxgc.index.get_level_values(level=0):
                combined[i] = gcxgc.loc[idx[i, :]].groupby(by='Name').sum()

            # Recombine to DataFrame
            gcxgc = pd.concat(combined).unstack(level=0).fillna(0)

        except:  # JSR hyperchrom data
            gcxgc = pd.read_hdf(
                store,
                'raw/GCxGC_FID',
                columns=['Area']
            ).unstack(level=0).fillna(0)

        # Drop multiIndex
        gcxgc.columns = gcxgc.columns.droplevel(level=0)
        # this should not be necessary if InChI conversion prior to
        # finding reference internal standard
        converter = {
            'CH4': 'methane',
        }
        gcxgc = gcxgc.rename(index=converter)

        # Store to the HDF5 store
        storing_2_HDF5_store(
            gcxgc,
            store,
            'processing/raw/GCxGC_FID',
        )

        # Extract only the chroms that are in the datapoints and
        # give the correct column names
        gcxgc_names = datapoints['GCxGC'].map(
            lambda x: x.replace('.csv', '').replace('.XLS', '')
        )

        gcxgc_raw_data = {}
        # Iterate over all datapoints
        for inj_name in gcxgc_names.index:
            gcxgc_raw_data[inj_name] = gcxgc[
                gcxgc_names[inj_name]
            ]

        # Recombine in DataFrame
        gcxgc_raw_data = pd.concat(gcxgc_raw_data, axis=1)

        # store the GCxGC data consistent with datapoints to the HDF5 store
        storing_2_HDF5_store(
            gcxgc_raw_data,
            store,
            'processing/GCxGC_FID',
        )

    # Read GCxGC SCD raw data from HDF5 file
    if '/raw/GCxGC_SCD' in store.keys():
        gcxgc_scd = pd.read_hdf(
            store,
            'raw/GCxGC_SCD',
            columns=['Volume', 'Name']
        )
        combi_scd = {}
        # Iterate over all chroms
        for i in gcxgc_scd.get_level_values(level=0):
            combi_scd[i] = gcxgc_scd.loc[idx[i, :]].groupby(by='Name').sum()
        gcxgc_scd = pd.concat(combi_scd).unstack(level=0).fillna(0)
        gcxgc_scd.columns = gcxgc_scd.columns.droplevel(level=0)
        # Store to the HDF5 store
        storing_2_HDF5_store(
            gcxgc_scd,
            store,
            'processing/raw/GCxGC_SCD',
        )

        # Extract only the chroms that are in the datapoints and
        # give the correct column names
        gcxgc_scd_names = datapoints['GCxGC_SCD'].map(
            lambda x: x.replace('.csv', '').replace('.XLS', '')
        )

        gcxgc_scd_data = {}
        # Iterate over all datapoints
        for inj_name in gcxgc_scd_names.index:
            gcxgc_scd_data[inj_name] = gcxgc_scd[
                gcxgc_scd_names[inj_name]
            ]

        # Recombine in DataFrame
        gcxgc_scd_data = pd.concat(gcxgc_scd_data, axis=1)

        # store the GCxGC-SCD data consistent with datapoints to the HDF5 store
        storing_2_HDF5_store(
            gcxgc_scd_data,
            store,
            'processing/GCxGC_SCD',
        )

    # Read GCxGC NCD raw data from HDF5 file
    if '/raw/GCxGC_NCD' in store.keys():
        gcxgc_ncd = pd.read_hdf(
            store,
            'raw/GCxGC_NCD',
            columns=['Volume', 'Name']
        )
        combi_ncd = {}
        # Iterate over all chroms
        for i in gcxgc_ncd.get_level_values(level=0):
            combi_ncd[i] = gcxgc_ncd.loc[idx[i, :]].groupby(by='Name').sum()
        gcxgc_ncd = pd.concat(combi_ncd).unstack(level=0).fillna(0)
        gcxgc_ncd.columns = gcxgc_ncd.columns.droplevel(level=0)
        # Store to the HDF5 store
        storing_2_HDF5_store(
            gcxgc_ncd,
            store,
            'processing/raw/GCxGC_NCD',
        )

        # Extract only the chroms that are in the datapoints and
        # give the correct column names
        gcxgc_ncd_names = datapoints['GCxGC_NCD'].map(
            lambda x: x.replace('.csv', '').replace('.XLS', '')
        )

        gcxgc_ncd_data = {}
        # Iterate over all datapoints
        for inj_name in gcxgc_ncd_names.index:
            gcxgc_ncd_data[inj_name] = gcxgc_ncd[
                gcxgc_ncd_names[inj_name]
            ]

        # Recombine in DataFrame
        gcxgc_ncd_data = pd.concat(gcxgc_ncd_data, axis=1)

        # store the GCxGC-NCD data consistent with datapoints to the HDF5 store
        storing_2_HDF5_store(
            gcxgc_ncd_data,
            store,
            'processing/GCxGC_NCD',
        )

    store.close()
    print(
        "\nGCxGC volume data extracted and stored to HDF5 file, "
        "if present."
    )
    return None


def processing_chrom_rga(h5filename, rga_channels, datapoints):
    '''Extracts the area data from the raw chromatograms and returns
    stores the DataFrame to the HDF5 file

    Parameters
    ----------
    h5filename : String
        The name of the HDF5 file without ending '.h5'
    rga_channels : MultiIndexed DataFrame
        They are also stored to the HDF5 file
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection

    Returns
    -------
    None
    '''
    import pandas as pd
    import sys

    try:
        # Deconvolute rga_channels
        rga_channelA, rga_channelB, rga_channelC = rga_channels

        rga_tcd_1 = rga_channelA['Area'].unstack(level=0).fillna(0)
        rga_tcd_2 = rga_channelB['Area'].unstack(level=0).fillna(0)
        if rga_tcd_2.columns.shape == rga_tcd_1.columns.shape:
            rga_tcd_2.columns = rga_tcd_1.columns
        else:
            popup_error(
                'TCD channels do not have equal length.\n'
                'Check if none of the chromatograms is empty.\n'
                'Program terminated.'
            )
            sys.exit(0)

        rga_tcd = pd.concat([rga_tcd_1, rga_tcd_2])

    except:  # JSR reactor has only 2 channels
        rga_channelA, rga_channelC = rga_channels
        rga_tcd = rga_channelA['Area'].unstack(level=0).fillna(0)

    # FID channel
    rga_fid = rga_channelC['Area'].unstack(level=0).fillna(0)

    store = pd.HDFStore('{}.h5'.format(h5filename))

    # store the extracted RGA areas to the HDF5 store
    storing_2_HDF5_store(
        rga_tcd,
        store,
        'processing/raw/RGA/TCD',
    )
    storing_2_HDF5_store(
        rga_fid,
        store,
        'processing/raw/RGA/FID',
    )

    # Extract only the chroms that are in the datapoints and
    # give the correct column names
    rga_fid_names = datapoints['RGA-FID'].map(
        lambda x: x.replace('.XLS', '').replace('.csv', '')
    )
    rga_fid_raw = pd.read_hdf(
        store,
        'processing/raw/RGA/FID',
        columns=rga_fid_names
    )
    rga_fid_names = {y: x for x, y in rga_fid_names.items()}
    # Change column names from rga_fid_raw
    fid_raw_data = rga_fid_raw.rename(columns=rga_fid_names)

    # TCD data
    rga_tcd_names = datapoints['RGA-TCD-L'].map(
        lambda x: x.replace('.XLS', '').replace('B', 'A').replace(
            '.csv', ''
        )
    )
    rga_tcd_raw = pd.read_hdf(
        store,
        'processing/raw/RGA/TCD',
        columns=rga_tcd_names
    )
    rga_tcd_names = {y: x for x, y in rga_tcd_names.items()}
    # Change column names from rga_fid_raw
    tcd_raw_data = rga_tcd_raw.rename(columns=rga_tcd_names)

    # store the RGA data consistent with datapoints to the HDF5 store
    storing_2_HDF5_store(
        tcd_raw_data,
        store,
        'processing/RGA/TCD',
    )
    storing_2_HDF5_store(
        fid_raw_data,
        store,
        'processing/RGA/FID',
    )

    store.close()

    print("Areas extracted from RGA raw chromatograms.")
    return None


def processing_log_file(df, mask, times):
    '''Processes the logging file by returning only the relevant data

    Extended description

    Parameters
    ----------
    df : DataFrame
        LoggingFile dataframe (timeSeries index)
    mask : ndarray, list
        Mask for the column names
    times : list of pandas Timestamps
        Mask for the injection times

    Returns
    -------
    df : DataFrame
        Extracted T, p-profile from logging file
    '''

    return df.resample('min').mean().loc[times, mask]


def processing_yields(
    datapoints,
    rga_channels,
    CFs,
    gen_configs,
    gcxgc_lib,
    fid_sum,
    h5filename
):
    ''' Process the experimental data to yields '''
    # Create a new class that contains all experimental yield data
    exp_data = yield_data(
        rga_channels,
        datapoints,
        CFs,
        gen_configs,
        gcxgc_lib
    )

    # Extract area and volume data from raw chromatograms
    processing_chrom_rga(h5filename, rga_channels, datapoints)
    processing_chrom_gcxgc_fid_scd_ncd(h5filename, datapoints)
    processing_chrom_loa(h5filename, datapoints)

    # Convert chromatograms from RGA, GCxGC, LOA to raw yields
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

    # Store all yield data to HDF5 file container
    exp_data.store_yield_data_to_hdf5(h5filename)

    # Add attributes to HDF5 file
    add_attributes(h5filename, exp_data)

    print('\nYield processing successfully finished.')
    print('Appending results to HDF5 file.')


def calculate_norm_2_wet(norm_yields, HC_feedrate, total_feedrate):
    '''Calculate the wet yields from the normalised dry yields

    Extended description

    Parameters
    ----------
    norm_yields : DataFrame
    HC_feedrate : Series
    total_feedrate : Series

    Returns
    -------
    wet_yields : DataFrame
    '''
    wet_yields = norm_yields.divide(
        total_feedrate
    ).multiply(HC_feedrate)

    return wet_yields


def calculate_mass_2_mol(mass_yields, mass_wt):
    '''Calculates the molar yields from a DataFrame with mass yields

    Extended description

    Parameters
    ----------
    mass_yields : DataFrame
    mass_wt : Series

    Returns
    -------
    mol_yields : DataFrame
    '''
    mol_yields = mass_yields.divide(mass_wt, axis=0)
    return mol_yields


def calculate_elemental_balance(
        mass_yields,
        mass_wt,
        atoms_in_comp,
        mass_atom
        ):
    '''Calculates the elemental effluent per component

    Extended description

    Parameters
    ----------
    mass_yields : DataFrame
    mass_wt : Series
    atoms_in_comp : Series
    mass_atom : Float
        The atomic mass of the element

    Returns
    -------
    elemental_balance : DataFrame
    '''
    elemental_balance = mass_yields.divide(
        mass_wt, axis=0
        ).multiply(
            atoms_in_comp, axis=0
        ) * mass_atom
    return elemental_balance


def processing_norm_2_wet(exp_data, fid_sum):
    '''Convert normalised dry yields to wet yields

    Extended description

    Parameters
    ----------
    exp_data : yield_data instance
        Holds all the information concerning the yield_data
    fid_sum : boolean
        Trigger for computing the sums with FID or TCD first

    Returns
    -------
    None
    '''

    # Get total HC inflow
    HC_feedrate = exp_data.HC_feeds.get_flow_total()

    # Calculate total mass influx's reciprocal
    total_feedrate = exp_data.get_feed_total_flow()

    # List of detectors
    detectors = exp_data.detectors

    # Convert dry yields to wet yields
    for detector in detectors:
        detectors[detector]['norm_yields_wet'] = calculate_norm_2_wet(
            detectors[detector]['norm_yields'],
            HC_feedrate,
            total_feedrate
            )

    # Calculate sums removing IS but no diluents
    is_inchi_name = []
    if exp_data.is_inchi not in exp_data.diluents.get_names_inchi():
        is_inchi_name = [exp_data.is_inchi]

    # Get list of detectors at given sheet for summation
    list_of_channels = exp_data.get_detectors_sheet('norm_yields_wet')

    # Compute the sums
    exp_data.sums[
        'norm_yields_wet'
    ], exp_data.total['norm_yields_wet'] = sum_generator(
        list_of_channels,
        is_inchi_name,
        fid_sum
        )

    print('Dry normalised yields converted to wet yields.')
    return None


def processing_elemental_balances(exp_data, gcxgc_lib, fid_sum):
    ''' Calculates the elemental balances for all detectors '''

    # imports
    import numpy as np

    # Calibration factors
    cfs = exp_data.get_CFs_inchi()

    # All elements and their atomic masses
    elements = 'CHONS'
    atomic_mass_elements = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'S': 32.065,
        'O': 15.999
    }

    # Get the HC_feedflows
    # use copy() → without modifying the original DF
    HC_feedflows = exp_data.HC_feeds.get_flows().transpose().copy()
    HC_inchis = exp_data.HC_feeds.get_names_inchi()
    total_HC_flow = exp_data.HC_feeds.get_flow_total()

    # Get the diluent_feedflows
    # use copy() → without modifying the original DF
    diluents_feedflows = exp_data.diluents.get_flows().transpose().copy()
    diluent_inchis = exp_data.diluents.get_names_inchi()
    total_dil_flow = exp_data.diluents.get_flow_total()

    # Initialisations
    ######################################
    for element in elements:

        # Inflows
        # Get elemental inflow rates for the HC feeds
        exp_data.HC_feeds.inflow[
            '{}_balance'.format(element)
            ] = HC_feedflows.divide(
                    cfs['GCxGC'].loc[HC_inchis, 'MW'], axis=0
                ).multiply(
                    cfs['GCxGC'].loc[HC_inchis, element], axis=0
                ).divide(
                    total_HC_flow
                ) * atomic_mass_elements[element] * 100

        # Get elemental inflow rates for the diluents
        exp_data.diluents.inflow[
            '{}_balance'.format(element)
            ] = diluents_feedflows.divide(
                    cfs['GCxGC'].loc[diluent_inchis, 'MW'], axis=0
                ).multiply(
                    cfs['GCxGC'].loc[diluent_inchis, element], axis=0
                ).divide(
                    total_dil_flow
                ) * atomic_mass_elements[element] * 100

        # Compute elemental balance DataFrame for all detectors
        for detector in exp_data.detectors:
            # Elemental balances on raw yield data
            exp_data.detectors[detector][
                '{}_balance'.format(element)
                ] = calculate_elemental_balance(
                    mass_yields=exp_data.detectors[detector]['raw_yields'],
                    mass_wt=cfs[detector]['MW'],
                    atoms_in_comp=cfs[detector][element],
                    mass_atom=atomic_mass_elements[element]
                    )
            exp_data.detectors[detector][
                '{}_balance'.format(element)
                ].dropna(how='any', inplace=True)

            # Elemental balances on normalised data
            exp_data.detectors[detector][
                '{}_balance_norm'.format(element)
                ] = calculate_elemental_balance(
                    mass_yields=exp_data.detectors[detector]['norm_yields'],
                    mass_wt=cfs[detector]['MW'],
                    atoms_in_comp=cfs[detector][element],
                    mass_atom=atomic_mass_elements[element]
                    )
            exp_data.detectors[detector][
                '{}_balance_norm'.format(element)
                ].dropna(how='any', inplace=True)

        # Sums
        # Calculate the sums for the elemental balances on raw yield data
        diluents = np.append(
            exp_data.diluents.get_names_inchi(),
            exp_data.is_inchi
        )

        # Get a list of the employed detectors for elemental balance sheet
        list_of_channels = exp_data.get_detectors_sheet(
            '{}_balance'.format(element)
            )

        exp_data.sums[
            '{}_balance'.format(element)
            ], exp_data.total['{}_balance'.format(element)] = sum_generator(
            list_of_channels,
            diluents,
            fid_sum
            )

        # Normalised sums
        # Get a list of the employed detectors for elemental balance sheet
        list_of_channels = exp_data.get_detectors_sheet(
            '{}_balance_norm'.format(element)
            )

        # Calculate the sums for the elemental balances
        # on normalised yield data
        exp_data.sums['{}_balance_norm'.format(element)], \
            exp_data.total['{}_balance_norm'.format(element)] = sum_generator(
                list_of_channels,
                diluents,
                fid_sum
            )

    print('Elemental balances calculated.')
    return None


def processing_norm_dry_molar(exp_data, fid_sum):
    '''Convert normalised weight percent yields to molar yields

    Extended description

    Parameters
    ----------
    exp_data : yield_data instance
    fid_sum : boolean

    Returns
    -------
    None
    '''

    # imports
    import numpy as np

    # Initialisations
    mol_subs = {}
    # Find diluents
    diluents = np.append(
        exp_data.diluents.get_names_inchi(),
        exp_data.is_inchi
        )

    # Calibration factors
    cfs = exp_data.get_CFs_inchi()

    # Sub dfs with mass to mols converted
    for detector in exp_data.detectors:
        mol_subs[detector] = calculate_mass_2_mol(
            mass_yields=exp_data.detectors[detector]['norm_yields'],
            mass_wt=cfs[detector]['MW']
            )
        # Remove empty lines
        mol_subs[detector].dropna(how='any', inplace=True)

    # Temporary sums
    mol_sub_sums, _ = sum_generator(
        mol_subs,
        diluents,
        fid_sum
        )

    # Iterate over the detectors to normalise the computed mol yields
    for detector in exp_data.detectors:
        exp_data.detectors[detector]['mol_perc'] = calculate_normalised(
            sheet=mol_subs[detector],
            sums=mol_sub_sums.loc['Total', :]
            )

    # get the list of channels to compute the sums
    list_of_channels = exp_data.get_detectors_sheet('mol_perc')
    # Sums mol percentage check for 100 %
    exp_data.sums['mol_perc'], exp_data.total['mol_perc'] = sum_generator(
        list_of_channels,
        diluents,
        fid_sum
        )

    print('Dry molar yields calculated.')
    return None


def processing_norm_wet_molar(exp_data, fid_sum):
    '''Convert normalised weight percent yields to wet molar yields

    Extended description

    Parameters
    ----------
    exp_data : yield_data instance
        Holds all the yield data
    fid_sum : boolean

    Returns
    -------
    None
    '''
    # Calibration factors
    cfs = exp_data.get_CFs_inchi()

    # Initialisations
    mol_subs = {}
    dil_inchis = exp_data.diluents.get_names_inchi()

    # Calculate relative molar influx of diluent
    dil_feedrate_mol = exp_data.get_fractions().divide(
        cfs['GCxGC'].loc[dil_inchis, 'MW']
        ).sum(axis=1) * 100

    # Calculate sub normalised yields wet
    for detector in exp_data.detectors:
        mol_subs[detector] = calculate_mass_2_mol(
            mass_yields=exp_data.detectors[detector]['norm_yields_wet'],
            mass_wt=cfs[detector]['MW']
            )
        # Remove empty lines
        mol_subs[detector].dropna(how='any', inplace=True)

    # Calculate sums removing IS but no diluents
    # find IS name
    is_inchi_name = []
    if exp_data.is_inchi not in exp_data.diluents.get_names_inchi():
        is_inchi_name = [exp_data.is_inchi]
    # Calculate actual sum
    mol_sub_wet_sums, _ = sum_generator(
        mol_subs,
        is_inchi_name,
        fid_sum
        )

    # Add amount of diluent
    ######################################
    # Enhancement
    ######################################
    # Should check if diluent is measured
    # with one of the detectors
    # otherwise this is not necessary
    ######################################
    mol_sub_wet_sums.loc['Total', :] = mol_sub_wet_sums.loc['Total', :] \
        + dil_feedrate_mol

    # Normalise molar sub DataFrames
    for detector in exp_data.detectors:
        exp_data.detectors[detector]['mol_perc_wet'] = calculate_normalised(
            sheet=mol_subs[detector],
            sums=mol_sub_wet_sums.loc['Total', :]
            )

    # Get the list of the detectors for calculating the sums
    list_of_channels = exp_data.get_detectors_sheet('mol_perc_wet')

    # Sum wet molar percentage
    exp_data.sums[
        'mol_perc_wet'
    ], exp_data.total['mol_perc_wet'] = sum_generator(
        list_of_channels,
        is_inchi_name,
        fid_sum
        )

    print('Wet molar yields calculated.')
    return None


def processing_chroms_2_raw(
    exp_data,
    rga_channels,
    gcxgc_lib,
    fid_sum,
    h5filename
):
    '''
    Convert chromatograms to raw yields (without normalisation)

    Parameters
    ----------
    exp_data :
    rga_channels :
    gcxgc_lib :
    fid_sum :
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None
    '''
    import pandas as pd
    import numpy as np
    import sys

    # Open HDF5 store
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Read RGA FID cleaned data
    exp_data.fid_raw_data = pd.read_hdf(
        store,
        'processing/RGA/FID',
    )

    # Read RGA TCD cleaned data
    exp_data.tcd_raw_data = pd.read_hdf(
        store,
        'processing/RGA/TCD',
    )

    # Read LOA data if available
    cond1 = '/processing/LOA' in store
    cond2 = exp_data.gen_configs["Process_LOA"] == 'TRUE'
    if cond1 and cond2:
        exp_data.loa_raw_data = pd.read_hdf(
            store,
            'processing/LOA',
        )

    # Read GCxGC data if available
    if '/processing/GCxGC_FID' in store:
        exp_data.gcxgc_fid_raw_data = pd.read_hdf(
            store,
            'processing/GCxGC_FID'
        )

    # Read GCxGC-SCD data if available
    cond1 = '/processing/GCxGC_SCD' in store
    cond2 = exp_data.gen_configs['Process_SCD'] == 'TRUE'
    if cond1 and cond2:
        exp_data.gcxgc_scd_raw_data = pd.read_hdf(
            store,
            'processing/GCxGC_SCD'
        )

    # Read GCxGC-NCD data if available
    cond1 = '/processing/GCxGC_NCD' in store
    cond2 = exp_data.gen_configs['Process_NCD'] == 'TRUE'
    if cond1 and cond2:
        exp_data.gcxgc_ncd_raw_data = pd.read_hdf(
            store,
            'processing/GCxGC_NCD'
        )

    # Read all internal standard data
    exp_data.read_IS_data()

    # Create converters for area/volume to raw yields
    primary_is_feedrate = exp_data.datapoints.loc[
        :,
        'PrimaryInternalStandard_FeedRate'
    ]

    # Get the total hydrocarbon inflow
    HC_feed_flow = exp_data.HC_feeds.get_flow_total()

    # Calculate constants
    tcd_area_to_yield = exp_data.primary_is_response.divide(
        exp_data.primary_is_area
    ).multiply(
        primary_is_feedrate
    ).divide(
       HC_feed_flow
    )

    # Convert RGA-TCD results to raw yields
    ######################################
    # Enhancement
    ######################################
    # Maybe this is not really necessary
    ######################################
    tcd_raw_data = duplicate_remover(exp_data.tcd_raw_data, first_last='first')
    # Remove duplicates
    tcd_area_to_yield = duplicate_remover(
        tcd_area_to_yield,
        first_last='first'
    )

    # Use the calculator for areas to raw yields
    exp_data.rga_tcd['raw_yields'] = calculate_areas_2_raw(
        tcd_raw_data,
        exp_data.CFs['TCD']['Response'],
        tcd_area_to_yield
        )

    # Transform indices to InChI
    exp_data.rga_tcd['raw_yields'].index = exp_data.CFs['TCD'].loc[
        exp_data.rga_tcd['raw_yields'].index, 'InChI'
        ].tolist()

    # Convert RGA-FID results to raw yields
    exp_data.rga_fid['raw_yields'] = calculate_areas_2_raw(
        raw_areas=exp_data.fid_raw_data,
        channel_RFs=exp_data.CFs['FID']['Response'],
        tcd_constant=tcd_area_to_yield,
        is_remote_area=exp_data.sec_is_area_tcd,
        is_channel_area=exp_data.sec_is_area_fid,
        is_channel_RF=exp_data.sec_is_response,
        is_remote_response=exp_data.sec_is_response_tcd,
    )

    # Add column Catoms for sorting on carbon number
    # Transform indices to InChI
    exp_data.rga_fid['raw_yields'].index = exp_data.CFs['FID'].loc[
        exp_data.rga_fid['raw_yields'].index, 'InChI'
        ].tolist()

    # An extra factor is necessary to convert TCD IS to FID
    # when using a reference from the FID
    tcd_area_to_yield_to_fid = tcd_area_to_yield.multiply(
            exp_data.sec_is_area_tcd
        ).multiply(
            exp_data.sec_is_response_tcd
        ).divide(
            exp_data.sec_is_area_fid
        ).divide(
            exp_data.sec_is_response
        )

    if 'LOA' in exp_data.datapoints:
        # A check which detector is chosen
        # because the tcd_area_to_yield is
        # multiplied with TCD->FID factors
        tcd_area_to_yield_loa = {}
        for inj_name in exp_data.datapoints.index:
            if exp_data.quat_in_tcd[inj_name]:  # Reference from TCD
                tcd_area_to_yield_loa[inj_name] = tcd_area_to_yield[inj_name]
            else:  # Reference detector chosen from FID
                tcd_area_to_yield_loa[inj_name] = tcd_area_to_yield_to_fid[
                    inj_name
                ]
        tcd_area_to_yield_loa = pd.Series(tcd_area_to_yield_loa)

        # Calculate LOA to percentages
        exp_data.loa['raw_yields'] = calculate_areas_2_raw(
            raw_areas=exp_data.loa_raw_data,
            channel_RFs=exp_data.CFs['LOA']['Response'],
            tcd_constant=tcd_area_to_yield_loa,
            is_remote_area=exp_data.quat_is_area_fid,
            is_channel_area=exp_data.quat_is_area_loa,
            is_channel_RF=exp_data.quat_is_response,
            is_remote_response=exp_data.quat_is_response_fid,
            )

        # Add column Catoms for sorting on carbon number
        # Transform indices to InChI
        exp_data.loa['raw_yields'].index = exp_data.CFs['LOA'].loc[
            exp_data.loa['raw_yields'].index,
            'InChI'
        ].tolist()

    if 'GCxGC' in exp_data.datapoints:
        # Convert GCxGC-FID results to raw yields
        if 'Factor_GCxGC' not in exp_data.datapoints:
            exp_data.datapoints['Factor_GCxGC'] = 1

        # A check which detector is chosen
        # because the tcd_area_to_yield is
        # multiplied with TCD->FID factors
        tcd_area_to_yield_gcxgc = {}
        for inj_name in exp_data.datapoints.index:
            if exp_data.tert_in_tcd[inj_name]:  # Reference from TCD
                tcd_area_to_yield_gcxgc[inj_name] = tcd_area_to_yield[inj_name]
            else:  # Reference detector chosen from FID
                tcd_area_to_yield_gcxgc[inj_name] = tcd_area_to_yield_to_fid[
                    inj_name
                ]
        tcd_area_to_yield_gcxgc = pd.Series(tcd_area_to_yield_gcxgc)

        # Create an output excel with all GCxGC_FID
        # components that are present in the chromatograms
        gcxgc_index = exp_data.gcxgc_fid_raw_data.index.unique().tolist()
        np.savetxt('GCxGC_fid_index.csv', gcxgc_index, delimiter=";", fmt="%s")
        # Calculate GCxGC data to percentages
        exp_data.gcxgc_fid['raw_yields'] = calculate_areas_2_raw(
            raw_areas=exp_data.gcxgc_fid_raw_data,
            channel_RFs=gcxgc_lib.loc[
                exp_data.gcxgc_fid_raw_data.index,
                'Response Factor'
            ],
            tcd_constant=tcd_area_to_yield_gcxgc,
            is_remote_area=exp_data.tert_is_area,
            is_channel_area=exp_data.tert_is_area_gcxgc,
            is_channel_RF=exp_data.tert_is_response,
            is_remote_response=exp_data.tert_is_response_tcd,
            ) * exp_data.datapoints['Factor_GCxGC']

        # Drop nans
        exp_data.gcxgc_fid['raw_yields'] = exp_data.gcxgc_fid[
            'raw_yields'
        ].dropna()
        # Sum blobs with the same name together
        exp_data.gcxgc_fid['raw_yields'] = exp_data.gcxgc_fid[
            'raw_yields'
            ].groupby(level=0).sum()

        # Check if all indices can be found in the library
        # If not create csv with not found components and stop program
        if not check_diff_index(exp_data.gcxgc_fid['raw_yields'], gcxgc_lib):
            # Get a list of not found components
            not_found = diff_index(exp_data.gcxgc_fid['raw_yields'], gcxgc_lib)
            # Write the list to a csv file
            np.savetxt(
                'not_found_components.csv',
                not_found,
                delimiter=';',
                fmt='%s'
            )
            # Show an error message
            popup_error(
                "Some components cannot be found in the GCxGC library."
                "Check the not_found_components file"
                )
            # Stop the program
            sys.exit()

        # Transform component names to InChI's
        exp_data.gcxgc_fid['raw_yields']['InChI'] = gcxgc_lib.loc[
            exp_data.gcxgc_fid['raw_yields'].index,
            'InChI'
            ].tolist()
        # Over different chromatograms the name of certain components
        # e.g. ethene-ethylene is altered
        # With following command they are summed together (based on InChI)
        exp_data.gcxgc_fid['raw_yields'] = exp_data.gcxgc_fid[
            'raw_yields'
            ].groupby('InChI').sum()

        # For the JSR, 1,3-butadiene should be taken from the C5+
        # dirty hack by copying the value to the the RGA-FID
        if exp_data.gen_configs['Reactor'] == 'JSR':
            exp_data.rga_fid['raw_yields'].loc[
                '1S/C4H6/c1-3-4-2/h3-4H,1-2H2'
            ] = exp_data.gcxgc_fid['raw_yields'].loc[
                '1S/C4H6/c1-3-4-2/h3-4H,1-2H2'
            ]

    # There is no GCxGC data available
    else:
        # Create empty dataframe
        exp_data.gcxgc_fid['raw_yields'] = pd.DataFrame(
            columns=exp_data.rga_tcd['raw_yields'].columns
        )

    # Get InChI name of primary internal standard to reject it in the sum
    # Unless it is used as a diluent
    prim_is_name = exp_data.datapoints['PrimaryInternalStandard_Name'][0]
    exp_data.is_inchi = gcxgc_lib.loc[prim_is_name, 'InChI']
    diluents = np.append(
        exp_data.diluents.get_names_inchi(),
        exp_data.is_inchi
    )

    # Calculate the sums
    list_of_channels = {
        'FID': exp_data.rga_fid['raw_yields'],
        'TCD': exp_data.rga_tcd['raw_yields'],
        'GCxGC': exp_data.gcxgc_fid['raw_yields']
        }
    ######################################
    # Enhancement
    ######################################
    # This should be different
    ######################################
    # Add LOA if available
    if 'LOA' in exp_data.datapoints:
        list_of_channels['LOA'] = exp_data.loa['raw_yields']
    # Add SCD if available
    if 'SCD' in exp_data.datapoints:
        list_of_channels['SCD'] = None  # Data not available yet
    # Add NCD if available
    if 'NCD' in exp_data.datapoints:
        list_of_channels['NCD'] = None  # Data not available yet

    exp_data.sums['raw_yields'], exp_data.total['raw_yields'] = sum_generator(
        list_of_channels,
        diluents,
        fid_sum
    )

    store.close()
    print('Raw yield results calculated.')
    return None


def calculate_areas_2_raw(
        raw_areas,
        channel_RFs,
        tcd_constant,
        is_remote_area=1,
        is_channel_area=1,
        is_channel_RF=1,
        is_remote_response=1
        ):
    '''
    Transforms the integrated areas to raw mass yields

    Parameters
    ----------
    raw_areas : DataFrame
        Contains the areas of the detected components.
        Indices : component names
        columns : Datapoint tags
    channel_RFs : Series
        Contains all response factors for the channel that is being processed
        Indices : Component names
    tcd_constant : Series
        Contains the partially processed TCD data
        Indices : Datapoint tags
    is_remote_area : Series
        Holds the area of the IS on the other channel
        Indices : Datapoint tags
    is_channel_area : Series
        Holds the area of the IS in the channel that is being processed
        Indices : Datapoint tags
    is_channel_RF : Series
        Holds the response factor of the IS on processing channel
        Indices : Datapoint tags
    is_remote_response : Series
        Holds the response factor of the IS on the remote channel
        Indices : Datapoint tags

    Returns
    -------
    raw_yields : DataFrame
        Contains the mass percentage raw yields
    '''
    raw_yields = raw_areas.multiply(
        channel_RFs,
        axis=0
    ).multiply(
        tcd_constant
    ).multiply(
        is_remote_area
    ).divide(
        is_channel_area
    ).divide(
        is_channel_RF
    ).multiply(
        is_remote_response
    ) * 100

    return raw_yields.dropna(how='all')


def calculate_normalised(sheet, sums, normalising_factor=100):
    '''Normalises the datapoints to 100% with given sum

    Extended description

    Parameters
    ----------
    sheet : DataFrame
        Detector yielddata
    sums : Series
        The sums of the datapoints

    Returns
    -------
    normalised : DataFrame
        Normalised yield data
    '''
    normalised = sheet / sums * normalising_factor
    return normalised


def processing_raw_2_norm(exp_data, fid_sum):
    '''
    Normalise raw yields mathematically to 100 %

    Parameters
    ---------
    exp_data: yield_data instance
    fid_sum: boolean

    Returns
    -------
    None

    '''

    import numpy as np

    # Normalising to 100 % for every detector
    for detector in exp_data.detectors:
        exp_data.detectors[detector]['norm_yields'] = calculate_normalised(
            sheet=exp_data.detectors[detector]['raw_yields'],
            sums=exp_data.sums['raw_yields'].loc['Total', :]
            )

    # Calculate sum to test normalisation
    diluents = np.append(
        exp_data.diluents.get_names_inchi(),
        exp_data.is_inchi
    )

    # Calculate the sums
    list_of_channels = exp_data.get_detectors_sheet('norm_yields')

    # Calculate the sums
    exp_data.sums[
        'norm_yields'
    ], exp_data.total['norm_yields'] = sum_generator(
        list_of_channels,
        diluents,
        fid_sum
    )

    print('Yields normalised.')
    return None


def processing_raw_2_carbon_normalized(exp_data, fid_sum):
    '''
    Carbon normalises the raw yields

    Parameters
    ----------
    exp_data: instance of yield_data
    fid_sum: boolean
        Should FID be taken as first values or TCD

    Returns
    -------
    None
    '''
    # imports
    import numpy as np

    # Calibration factors
    cfs = exp_data.get_CFs_inchi()

    # All elements and their atomic masses
    elements = 'CHONS'
    atomic_mass_elements = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'S': 32.065,
        'O': 15.999
    }

    # Only normalise based on C
    norm_el = 'C'

    # Normalise first raw carbon balance

    try:  # 1 feed
        HC_total_el_inflow = exp_data.HC_feeds.get_total_elemental_flow(
            '{}_balance'.format(norm_el)
            ).iloc[:, 0]  # iloc changes df to Series
    except:  # Multiple feed
        HC_total_el_inflow = exp_data.HC_feeds.get_total_elemental_flow(
            '{}_balance'.format(norm_el)
            )

    # Calculate the normalised el balance for all detectors
    for detector in exp_data.detectors:
        exp_data.detectors[detector][
            'norm_{}_balance'.format(norm_el)
            ] = calculate_normalised(
                sheet=exp_data.detectors[detector][
                    '{}_balance'.format(norm_el)
                    ],
                sums=exp_data.sums[
                    '{}_balance'.format(norm_el)
                    ].loc['Total', :],
                normalising_factor=HC_total_el_inflow
                )

    # Calculate sums
    diluents = np.append(
        exp_data.diluents.get_names_inchi(), exp_data.is_inchi
        )

    # Create a list of the channels
    list_of_channels = exp_data.get_detectors_sheet(
        'norm_{}_balance'.format(norm_el)
        )

    # Calculate the sums
    exp_data.sums[
        'norm_{}_balance'.format(norm_el)
        ], exp_data.total['norm_{}_balance'.format(norm_el)] = sum_generator(
        list_of_channels,
        diluents,
        fid_sum
        )

    # Normalise yields based on closed carbon balance
    for detector in exp_data.detectors:
        exp_data.detectors[detector][
            '{}_norm_yields'.format(norm_el)
            ] = calculate_normalised(
                sheet=exp_data.detectors[detector]['raw_yields'],
                sums=exp_data.sums[
                    '{}_balance'.format(norm_el)
                ].loc['Total', :],
                normalising_factor=HC_total_el_inflow
                )
        # Drop empty lines
        exp_data.detectors[detector]['{}_norm_yields'.format(norm_el)].dropna(
            how='any', inplace=True
            )

    # Get a list of the used detectors for el_norm_yields sheet
    list_of_channels = exp_data.get_detectors_sheet(
        '{}_norm_yields'.format(norm_el)
        )
    # Calculate sums
    exp_data.sums[
        '{}_norm_yields'.format(norm_el)
        ], exp_data.total['{}_norm_yields'.format(norm_el)] = sum_generator(
        list_of_channels,
        diluents,
        fid_sum
        )

    # Calculate elemental balances after carbon normalisation
    for element in elements:
        sheetname = '{}_balance_after_{}_norm'.format(element, norm_el)

        for detector in exp_data.detectors:
            # Calculate the elemental balance
            exp_data.detectors[detector][
                sheetname
                ] = calculate_elemental_balance(
                    mass_yields=exp_data.detectors[detector][
                        '{}_norm_yields'.format(norm_el)
                        ],
                    mass_wt=cfs[detector]['MW'],
                    atoms_in_comp=cfs[detector][element],
                    mass_atom=atomic_mass_elements[element]
                    )
            exp_data.detectors[detector][sheetname].dropna(
                how='any',
                inplace=True
            )

        # Get the list of used channels
        list_of_channels = exp_data.get_detectors_sheet(sheetname)
        # Calculate the sums for the elemental balances on raw yield data
        exp_data.sums[sheetname], exp_data.total[sheetname] = sum_generator(
            list_of_channels,
            diluents,
            fid_sum
            )

    print('Yields are carbon normalised.')

    return None


def sum_generator(list_dfs, internal_standards, fid_sum=1):
    '''
    Module calculates the total sum of the inserted list of dataframes

    Parameters
    ----------
    list_dfs : dict
        {'FID': df_fid, 'TCD': df_tcd, 'GCxGC': df_gcxgc}
        from FID channel > TCD channel
    internal_standards : lst
        list of internal standard names
    fid_sum: boolean

    Returns
    -------
        sums: DataFrame
            C4, C5, Total
        total: DataFrame
            summed yield overview
    '''
    import pandas as pd
    import numpy as np

    # Create a list with the names
    if fid_sum:
        c4_names = ['FID', 'TCD']  # FID > TCD channel (RGA)
    else:
        c4_names = ['TCD', 'FID']  # TCD < FID channel
    # If LOA is used implement it in the c4- names
    if 'LOA' in list_dfs:
        # Remove CO/CH4 peak
        try:
            list_dfs['LOA'].drop('1S/CH4CO/c1-2/h2H,1H3', inplace=True)
        except:
            pass
        c4_names.append('LOA')
    # Take all components for total sum
    total_names = c4_names.copy()
    total_names.append('GCxGC')
    # If SCD is available add to total sums
    if 'GCxGC-SCD' in list_dfs:
        total_names.append('GCxGC-SCD')
    # If NCD is available add to total sums
    if 'GCxGC-NCD' in list_dfs:
        total_names.append('GCxGC-NCD')

    # Create a list of DataFrames with the names for C4- and total
    c4_list = [list_dfs[name] for name in c4_names]
    total_list = [list_dfs[name] for name in total_names]

    # Make a DataFrame from the list of DataFrames
    c4 = pd.concat(c4_list)
    total = pd.concat(total_list)

    # Remove measured diluents from the sum
    for internal_standard in internal_standards:
        if internal_standard in c4.index:
            c4.drop(internal_standard, inplace=True)
    c4.replace([np.inf, -np.inf], np.nan, inplace=True)
    c4.dropna(inplace=True)
    c4_sumsub = duplicate_remover(c4, first_last='first')  # .sum()
    c4_sum = c4_sumsub.sum()

    # Remove measured diluents from the sum
    for internal_standard in internal_standards:
        if internal_standard in total.index:
            total.drop(internal_standard, inplace=True)
    total = duplicate_remover(total, first_last='first')
    total.replace([np.inf, -np.inf], np.nan, inplace=True)
    total.dropna(inplace=True)
    total_sum = total.sum()

    sums = pd.concat([c4_sum, total_sum-c4_sum, total_sum], axis=1).transpose()
    sums.index = ['C4', 'C5', 'Total']

    return sums, total


def processing_loggingfiles_2_temperature_pressure_profiles(
    datapoints,
    gen_configs,
    h5filename
):
    '''Extracts the temperature and pressure profiles
    during specific injections

    Extended description

    Parameters
    ----------
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection
    gen_configs : Series
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None
    '''
    import pandas as pd

    # Get the name of the reactor and thus reactor configuration file
    reactor = gen_configs['Reactor']
    # Reactor configuration
    try:
        reactor_config = reading_reactor_config(
            r'N:\Project\918-Experimental set-ups\Reactor_configurations',
            r'{}.xlsx'.format(reactor)
        )
    except FileNotFoundError:  # offline version of the software
        try:
            reactor_config = reading_reactor_config(
                gen_configs['reactor_config_loc'],
                r'{}.xlsx'.format(reactor)
            )
        except FileNotFoundError:  # location not given in gen_configs
            sys.exit()

    # Initialisations
    profile_names = ['temperature', 'pressure']
    profiles = {
        'temperature': [],
        'pressure': []
    }

    # Geometry
    geometry = reactor_config['geometry']
    # Get the logging file names (T and p)
    loggingFileName = {
        'temperature': reactor_config['log_name_temperature'],
        'pressure': reactor_config['log_name_pressure']
    }
    # Get the tags for the temperatures and pressures
    mask = {
        'temperature': reactor_config['names_temperature'],
        'pressure': reactor_config['names_pressure']
    }

    # Group the datapoints based on the date
    datapoints_per_date = datapoints.groupby('Date')

    for group in datapoints_per_date.groups:
        # Get a list of the injection times as a Pandas Timestamp
        time = (group + ' ' + datapoints_per_date.get_group(group).loc[
                        :, 'Time'
            ])
        time = pd.to_datetime(time, dayfirst=True)

        # Iterate over all profiles (T and p)
        for profile_name in profile_names:

            # Old BSSC does not have pressure log file
            if not pd.isnull(loggingFileName[profile_name]):
                # Read the logging file from the HDF5 store
                log = pd.read_hdf(
                    '{}.h5'.format(h5filename),
                    'raw/logging/{}_{}'.format(
                        loggingFileName[profile_name],
                        group.replace('/', '_')
                        ),
                    )

                # Extract the profiles during injection times
                # from the loggingfile
                profile = processing_log_file(
                    log,
                    mask[profile_name],
                    time
                )

                # Change pandas datetime index to the experimental labels
                profile.index = datapoints_per_date.get_group(group).index

                # Add the profile to the profiles dict
                profiles[profile_name].append(profile)

    # Make some changes on the profiles
    for profile_name in profile_names:
        # BSSC old has no pressures logged
        if pd.isnull(loggingFileName[profile_name]):
            # Create own DataFrame with zeros
            profiles[profile_name] = pd.DataFrame(
                0,
                index=['inlet', 'outlet'],
                columns=datapoints.index
                )
        else:
            # Concatenate all different days
            profiles[profile_name] = pd.concat(
                profiles[profile_name]
            ).transpose()

        # Add axial positions to the profiles
        profiles[profile_name].index = reactor_config['axial_{}s'.format(
            profile_name
            )].loc[
                profiles[profile_name].index,
                'axial_position[m]'
                ].values

    # Add 1 to the pressure profiles if the measurements are
    # recorded in the logging files in barg instead of bara
    # Fix for BSSC old (Add 1.7 bara)
    profiles['pressure'] = profiles['pressure'] + reactor_config['barg']

    # Divide by 10 if the measurements on the BSSC quartz reactor
    # are performed before March 2017, due to bad logging of the T
    profiles['temperature'] = profiles['temperature'] \
        / reactor_config['divide_by']

    # The inlet and outlet values for the metal reactor are
    # not logged so they are manually added
    if reactor == 'BSSC_metal':
        profiles['temperature'].loc[0, :] = 400
        profiles['temperature'].loc[1.49, :] = 300

    # Check if all inlet values are higher than the outlet values
    # if not make the inlet pressure == outlet + 0.01
    # This is sometimes a problem for the BSSC quartz reactor
    p_prof_tp = profiles['pressure'].transpose()
    cond = p_prof_tp.iloc[:, 0] > p_prof_tp.iloc[:, -1]
    profiles['pressure'] = p_prof_tp.where(
        cond,
        p_prof_tp.iloc[:, -1] + 0.01,
        axis=0
        ).transpose()

    # Save the relevant data to the HDF5 store
    # open the store
    store = pd.HDFStore('{}.h5'.format(h5filename))

    # Save geometry to HDF5 file
    storing_2_HDF5_store(
        geometry,
        store,
        key='processed/{}'.format('geometry'),
        )

    # Save profiles to HDF5 file
    for profile_name in profile_names:
        # Sort the index and columns
        profiles[profile_name].sort_index(inplace=True)
        profiles[profile_name].sort_index(inplace=True, axis=1)

        # Save the profiles to the HDF5 file
        storing_2_HDF5_store(
            profiles[profile_name],
            store,
            key='processed/profiles/{}'.format(profile_name),
            )

    # Close the HDF5 store
    store.close()

    print(
        "\nTemperature and pressure profiles extracted and stored "
        "in HDF5 file: {}.h5".format(h5filename)
    )
    return None


def processing_loggingfiles_extract_cracking(h5filename, reactor):
    '''Extract the cracking data for all the loggingfiles

    Parameters
    ----------
    h5filename : string
        The name of the HDF5 file without ending '.h5'
    reactor : string
        The name of the reactor/setup to read the appropriate reactor config

    Returns
    -------
    None
        Adds the extracted data to the HDF5 store
    '''
    import pandas as pd

    # Check if the old BSSC reactor is used
    old = True if reactor == 'BSSC_metal_old' else False

    # Open the store
    store = pd.HDFStore('{}.h5'.format(h5filename))

    logging_file_names = [
        log for log in store.keys() if 'raw/logging' in log
    ]

    start_time = store.root._v_attrs['Start_time']
    end_time = store.root._v_attrs['End_time']

    for log_file_name in logging_file_names:
        df_log = pd.read_hdf(
            store,
            log_file_name
        )

        # Extract the relevant cracking time data
        df_extracted = process_start_end_log_file(
            df_log,
            start_time,
            end_time,
            old=old
        )

        # Save all data to HDF5 file
        storing_2_HDF5_store(
            df_extracted,
            store,
            key='processed/log_file_cracking/{}'.format(log_file_name.replace('/raw/logging/', '')),
        )

    # Close the HDF5 store
    store.close()

    print(
        "\nCracking data extracted from all logging files "
        "and stored in the HDF5 file: {}.h5".format(h5filename)
    )
    return None


def process_start_end_log_file(df_log, time_start, time_end, old=False):
    '''Read in the DataFrame of the logginfile and
    extract the data inbetween the time_start and time_end
    Change the index to seconds beginning from time_start
    and use the timestep from the loggingfile itself

    Parameters
    ----------
    df_log : DataFrame
        DataFrame from the loggingfile
    time_start : datetime, string
        Start time of the experiment
        String format : YYYY-MM-DD hh:mm
    time_end : datetime, string
        End time of the experiment
        String format : YYYY-MM-DD hh:mm
    old : boolean, optional
        Old logging file structure

    Returns
    -------
    df_sub : DataFrame
        Modified loggingfile data.
        Data extracted between start and end time
        of the experiment.
    '''
    import numpy as np
    import pandas as pd

    # Convert start and end time to Pandas Timestamp if necessary
    try:
        if not isinstance(time_start, pd.Timestamp):
            time_start = pd.to_datetime(time_start, dayfirst=True)
        if not isinstance(time_end, pd.Timestamp):
            time_end = pd.to_datetime(time_end, dayfirst=True)
    except ValueError:
        print(
            "Either start or end time is incorrectly written.\n"
            "The string format should be YYYY-MM-DD hh:mm"
        )
        return

    # Downsample old logging file structure to every 5 seconds
    if old:
        df_log = df_log.resample('5S').mean()

    # Get lower and upper limits for the extraction
    ll = df_log.index > time_start
    ul = df_log.index <= time_end
    cond = ll & ul

    # Extract given condition
    df_sub = df_log.loc[cond, :]
    df_sub.reset_index(drop=True, inplace=True)

    # Create new time (in seconds) index based on the
    # timestep found in the loggingfile
    if old:
        timestep = 5
    else:
        timestep = int(np.rint(
            (df_log.index[-1] - df_log.index[0]).total_seconds() \
            / (len(df_log.index) - 1)
        ))
    df_sub.index = np.arange(
        start=0,
        stop=len(df_sub.index) * timestep,
        step=timestep
    )
    df_sub.index.name = 'Time (s)'

    return df_sub


def processing_loggingfiles_jsr(
    datapoints,
    gen_configs,
    h5filename
):
    '''Extracts the temperature and pressure profiles
    during specific injections

    Extended description

    Parameters
    ----------
    datapoints : DataFrame
        Contains the relevant information per datapoint/injection
    gen_configs : Series
    h5filename : string
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None
    '''
    print("***No JSR loggingfiles processing implemented yet.***")


def repack(h5filename):
    '''
    Uses the utility ptrepack from PyTables (in the command line)
    to compress the data

    Parameters
    ----------
    h5filename : String
        The name of the HDF5 file without ending '.h5'

    Returns
    -------
    None

    '''
    import os
    import subprocess
    import sys

    # Start the ptrepack utility in the command line
    process_line = (
        r"ptrepack --chunkshape=auto --complevel=9 --complib=zlib "
        r"--overwrite {}.h5 temp.h5".format(h5filename)
    )
    if 'win' in sys.platform:
        subprocess.call(process_line)
    elif sys.platform == 'linux':
        subprocess.call(process_line, shell=True)

    # Rename it to original name and remove old file
    os.remove('{}.h5'.format(h5filename))
    os.rename('temp.h5', '{}.h5'.format(h5filename))

    print(
        "\n{}.h5 successfully repacked to "
        "a smaller form factor.".format(h5filename)
    )
    return None


def check_add_attribute(tag_name, text, store, exp_data, selection_list=None):
    '''
    Adds an attribute with given tag_name
    If the tag_name was not provided in the gen_configs file
    a pop-up window will ask to provide it
    '''
    if tag_name in exp_data.gen_configs:
        if selection_list is not None:
            if exp_data.gen_configs[tag_name] not in selection_list:
                exp_data.gen_configs.drop(tag_name, inplace=True)
                check_add_attribute(
                    tag_name,
                    text,
                    store,
                    exp_data,
                    selection_list
                )
            else:
                store.root._v_attrs[tag_name] = exp_data.gen_configs[tag_name]
        else:
            store.root._v_attrs[tag_name] = exp_data.gen_configs[tag_name]
    else:
        # It should be possible to write the selection more elegantly
        if selection_list is None:
            inserted_attr = popup_entry(text)
        else:
            inserted_attr = None
            while inserted_attr not in selection_list:
                inserted_attr = popup_entry(text)
        exp_data.gen_configs[tag_name] = inserted_attr
        store.root._v_attrs[tag_name] = exp_data.gen_configs[tag_name]


def add_attributes(h5filename, exp_data):
    ''' Add attributes/metadata of experimental data to HDF5 file container '''

    import pandas as pd
    import time

    store = pd.HDFStore('{}.h5'.format(h5filename))
    # program version
    store.root._v_attrs[
        'GC_processor_HDF5_version'
    ] = GC_processor_HDF5_version
    # Feeds
    store.root._v_attrs['Feeds'] = exp_data.HC_feeds.get_names_smiles()
    store.root._v_attrs[
        'Feed_rate_total_flow_gh'
    ] = exp_data.HC_feeds.get_flow_total().unique()
    feeds_dct = {}
    feeds_flows = exp_data.HC_feeds.get_flows()
    for feed in feeds_flows.columns:
        feeds_dct[feed] = feeds_flows[feed].unique()
    store.root._v_attrs['Feeds_flows_gh'] = feeds_dct
    # Diluents
    store.root._v_attrs['Diluents'] = exp_data.diluents.get_names_smiles()
    store.root._v_attrs[
        'Diluent_total_flow_gh'
    ] = exp_data.diluents.get_flow_total().unique()
    diluents_dct = {}
    diluents_flows = exp_data.diluents.get_flows()
    for dil in diluents_flows.columns:
        diluents_dct[dil] = diluents_flows[dil].unique()
    store.root._v_attrs['Diluents_flows_gh'] = diluents_dct
    # Temperatures
    store.root._v_attrs['Temperatures_C'] = list(
        exp_data.datapoints.groupby('Temperature').groups.keys()
    )
    # LOA detector
    cond = exp_data.gen_configs['Process_LOA'] == 'TRUE'
    store.root._v_attrs['Detector_LOA'] = 'yes' if cond else 'no'
    # SCD detector
    cond = exp_data.gen_configs['Process_SCD'] == 'TRUE'
    store.root._v_attrs['Detector_SCD'] = 'yes' if cond else 'no'
    # NCD detector
    cond = exp_data.gen_configs['Process_NCD'] == 'TRUE'
    store.root._v_attrs['Detector_NCD'] = 'yes' if cond else 'no'
    # Processing date
    store.root._v_attrs['Date_processing'] = time.strftime("%x")
    # Reactor
    store.root._v_attrs['Reactor'] = exp_data.gen_configs['Reactor']
    if exp_data.gen_configs['Reactor'] == 'Pilot':
        check_add_attribute(
            'Pilot_type_experiment',
            "Was it a 'coke' or 'yield' study",
            store,
            exp_data,
            ['coke', 'yield']
        )
    # Dates
    store.root._v_attrs['Dates'] = list(
        exp_data.datapoints.groupby('Date').groups.keys()
    )
    # Check and add attributes
    tag_names = {
        # Researcher who processed the experiments
        'AAA_Researcher_exp': ['Please add the experimental researcher', None],
        # Researcher who processed the experiments
        'AAA_Researcher_proc': ['Please add your name as processor', None],
        # Confidentiality tag
        'Confidential': [
            'Is the data confidential (type TOP, yes, or no)',
            ['yes', 'no', 'TOP']
        ],
        # Start time of the experiment
        'Start_time': [
            'Provide the starting time of the (cracking experiment)', None
        ],
        # End time of the experiment
        'End_time': [
            'Provide the ending time of the (cracking experiment)', None
        ],
        # Published data or not
        'Published': [
            'Is the data been published? (type yes, or no)',
            ['yes', 'no']
        ],
        # Experimental remarks
        'Remarks': ['Are there any experimental remarks/comments?', None],
        # Reliability of the data
        'Reliability_1_to_5': [
            'How reliable is the data set? (type: 1-5, with 1 being very bad)',
            ['{}'.format(i+1) for i in range(5)]
        ],
        # Mass spectrum available
        'Mass_spectrum_available': [
            'Is a mass spectrum available? (yes, no)',
            ['yes', 'no']
        ],
        # Pressures
        'Pressures_bara': ['Please insert the pressure (bara)', None],
        # Campaign or kinetic study
        'Campaign_kinetic': [
            "'campaign' or 'kinetic' study",
            ['campaign', 'kinetic']
        ]
    }

    for tag_name, txt in tag_names.items():
            check_add_attribute(
                tag_name,
                txt[0],
                store,
                exp_data,
                selection_list=txt[1]
            )
    # Date for confidentiality
    if exp_data.gen_configs['Confidential'] in ['TOP', 'yes']:
        check_add_attribute(
            'Confidentiality_date',
            'Please type in the confidentiality date (DD/MM/YY)',
            store,
            exp_data,
            None
            )
    # Additional article information
    if exp_data.gen_configs['Published'] == 'yes':
        check_add_attribute(
            'doi',
            'Please enter the doi of the article',
            store,
            exp_data,
            None
            )

    ######################################
    # Attributes
    ######################################
    # Add printing functionality
    ######################################

    print('\nAttributes added.')
    store.close()


def create_excel_summary(
    h5file,
    keys,
    file_name,
    foldername='excelSummaries',
    lib=None
):
    '''Creates summary excel file from all detectors

    Extended description

    Parameters
    ----------
    h5file : String
    keys : List
    file_name : String
    foldername : String
    lib : DataFrame or None

    Returns
    -------
    None
    '''
    from functools import reduce
    import numpy as np
    import pandas as pd
    import sys
    idx = pd.IndexSlice

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    if 'win' in sys.platform:
        writer = pd.ExcelWriter(
            r"{}\{}.xlsx".format(foldername, file_name),
            engine='xlsxwriter'
        )
    elif sys.platform == 'linux':
        writer = pd.ExcelWriter(
            r"{}/{}.xlsx".format(foldername, file_name),
            engine='xlsxwriter'
        )

    # Initialisations
    dfs = {}

    # Read in the DataFrames for given keys
    for key in keys:
        df = pd.read_hdf(h5file, key=key)
        # Convert InChI to Smiles
        if 'sums' not in key:
            unique_inchis = df.index.levels[1]
            inchi_converter = {
                inchi: inchi_2_smiles(inchi) for inchi in unique_inchis
                }
            df.index = df.index.set_levels(
                unique_inchis.map(inchi_converter.get),
                level='molecules'
                )
            # Add column with number of carbon atoms
            # Take just the sum of C's in Smiles name
            df['Catoms'] = df.index.get_level_values('molecules').map(
                lambda x: str(x).count('C') + str(x).count('c')
                )  # RDKIT has sometimes small c's in Smiles name

        dfs[key] = df

    # get the union of sheetnames (should be identical but better to check it)
    sheets = reduce(
        np.union1d,
        [df.index.get_level_values('sheets') for df in dfs.values()]
        )

    # Iterate over all sheets
    for sheetname in sheets:

        # Initialise
        df_list = []

        # Iterate over the keys
        for key in keys:
            # sums df should be treated differently
            if 'sums' not in key:
                # Sort the df based on the Catoms column
                # Delete Catoms column afterwards
                small_df = dfs[key].loc[idx[sheetname, :]].sort_values(
                    'Catoms', axis=0
                ).drop(labels='Catoms', axis=1)
            else:
                small_df = dfs[key].loc[idx[sheetname, :]]
            # Add empty line with detector description
            empty_name = pd.DataFrame(
                index=[key.split('/')[-1]],
                columns=small_df.columns
                )
            df_list.append(empty_name)
            # Add small df to df list
            df_list.append(small_df)
            # Add empty line
            empty_line = pd.DataFrame(
                index=[None],
                columns=small_df.columns
                )
            df_list.append(empty_line)

        # Create big DataFrame per sheet
        joined_data = pd.concat(df_list[:-1], axis=0)

        # Write to Pandas ExcelWriter
        joined_data.to_excel(writer, sheet_name=sheetname)

    # Add Excel formatting
    writer = add_excel_formatting(writer)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
    return None


def create_excel_profile(h5file, key, file_name, foldername='excelSummaries'):
    '''Creates an excel of the axial profiles

    Extended description

    Parameters
    ----------
    h5file : String
    key : String
    file_name : String
    foldername : String

    Returns
    -------
    None
    '''
    import pandas as pd
    import sys

    # Create a Pandas Excel Writer
    if 'win' in sys.platform:
        writer = pd.ExcelWriter(
            r'{}\{}.xlsx'.format(foldername, file_name),
            engine='xlsxwriter',
        )
    elif sys.platform == 'linux':
        writer = pd.ExcelWriter(
            r'{}/{}.xlsx'.format(foldername, file_name),
            engine='xlsxwriter',
        )

    # Dump profile to writer
    pd.read_hdf(h5file, key).to_excel(writer)

    # Apply number formatting to the Excel spreadsheet
    writer = add_excel_formatting(writer)

    # Close the Pandas Excel writer
    writer.save()

    return None


def create_excel(
    h5file,
    key,
    file_name,
    foldername='excelSummaries',
    lib=None
):
    ''' Creates excel files from the HDF5 file '''
    import pandas as pd
    import sys
    idx = pd.IndexSlice

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    if 'win' in sys.platform:
        writer = pd.ExcelWriter(
            r"{}\{}.xlsx".format(foldername, file_name),
            engine='xlsxwriter'
        )
    elif sys.platform == 'linux':
        writer = pd.ExcelWriter(
            r"{}/{}.xlsx".format(foldername, file_name),
            engine='xlsxwriter'
        )

    # Read data from HDF5 file
    df = pd.read_hdf(h5file, key=key)
    sheets = df.index.get_level_values('sheets').unique()

    # Convert InChI's to Smiles
    if lib is not None:
        unique_inchis = df.index.levels[1]
        inchi_converter = {
            inchi: inchi_2_smiles(inchi_name=inchi)
            for inchi
            in unique_inchis
        }
        df.index = df.index.set_levels(
            unique_inchis.map(inchi_converter.get),
            level='molecules'
        )
        # Sort Smiles based on carbon number
        df['Catoms'] = df.index.get_level_values('molecules').map(
            lambda x: str(x).count('C') + str(x).count('c')
        )  # RDKIT has sometimes small c's in smiles name

    # Write to excel writer
    for sheetname in sheets:
        # Write to excelWriter
        if lib is not None:
            small_df = df.loc[idx[sheetname, :]].sort_values(
                'Catoms',
                axis=0
            ).drop(labels='Catoms', axis=1)
        else:  # This will be the sums file
            small_df = df.loc[idx[sheetname, :]]
        small_df.to_excel(writer, sheet_name=sheetname)

    # Apply number formatting to the Excel spreadsheet
    writer = add_excel_formatting(writer)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


def create_excel_files(h5filename, lib):
    ''' Creates all possible excels from HDF5 file '''
    import os
    import pandas as pd

    print("Creating excel files")

    # Check if folder exists otherwise create it
    foldername = 'excelSummaries'
    if not os.path.exists(foldername):
        os.makedirs(foldername)

    # Change library
    lib = lib.set_index('InChI')
    lib = duplicate_remover(lib)

    # Get all group names in HDF5 file
    h5file = '{}.h5'.format(h5filename)
    store = pd.HDFStore(h5file)
    yieldnames = [
        name.replace('/processed/yields/', '')
        for name
        in store.keys()
        if '/processed/yields' in name
    ]
    store.close()

    # Yields
    for name in yieldnames:
        if name in ['sums']:
            create_excel(
                h5file,
                'processed/yields/{}'.format(name),
                '{}_{}'.format(name, h5filename),
                foldername=foldername,
                lib=None
            )
        else:
            ######################################
            # Enhancement
            ######################################
            # Maybe this is not really necessary
            # -> Creates a lot of files..
            # Only total and summary should be
            # sufficient, I guess
            ######################################
            create_excel(
                h5file,
                'processed/yields/{}'.format(name),
                '{}_{}'.format(name, h5filename),
                foldername,
                lib
            )
    # Create yields summary
    names = ['RGA_TCD', 'RGA_FID', 'GCxGC_FID', 'sums']
    store = pd.HDFStore(h5file)
    if 'processed/yields/LOA' in store:
        # Insert LOA before the second element in the list
        names.insert(2, 'LOA')
    store.close()
    keys = ['processed/yields/{}'.format(name) for name in names]
    create_excel_summary(
        h5file,
        keys,
        '{}_{}'.format('overview', h5filename),
        foldername,
        lib
    )

    # T,p-profiles
    for profile in ['temperature', 'pressure']:
        try:
            create_excel_profile(
                h5file,
                key='processed/profiles/{}'.format(profile),
                file_name=r'{}s_{}'.format(profile, h5filename),
                foldername=foldername
            )
        except KeyError:
            print("No temperature and pressure profiles found.")

    print(
        '\nExcel files successfully written and saved in {}.'.format(
            os.path.join(
                os.getcwd(), foldername
            )
        )
    )


def copy2database_folder(h5filename, gen_configs):
    """Copies the formed HDF5 folder to the central database
    folder of the setup
    """
    import shutil

    reactor = gen_configs.loc['Reactor', :].values[0]
    folder_location = None

    cond1 = reactor == 'BSSC_metal'
    cond2 = reactor == 'BSSC_metal_old'
    cond3 = reactor == 'BSSC_quartz'
    if cond1 or cond2 or cond3:
        folder_location = r'N:\Setup\BSSC\Database'

    elif reactor == 'pilot' or reactor == 'pilot - IMPROOF':
        folder_location = r'N:\Setup\Pilot\Database'

    elif reactor == 'JSR':
        folder_location = r'N:\Setup\JSR\Database'

    elif reactor == 'IHEB':
        folder_location = r'N:\Setup\IHEB\Database'

    if h5filename.endswith('.h5'):
        h5filename = h5filename.replace('.h5', '')

    try:
        shutil.copy('{}.h5'.format(h5filename), folder_location)
        print(
            '\n{}.h5 successfully copied to {}.'.format(
                h5filename, folder_location
            )
        )
    except(shutil.SameFileError):
        print('\n{}.h5 not copied.'.format(h5filename))


def process_experimental_results(filename, fid_sum=1):
    '''Starts the experimental results processor and produces a HDF5 file

    Parameters
    ----------
    filename : String
        The names of the input csv files conditions and gen_configs
        separated by a space
    fid_sum : int
        1: Take the yield sum with a preference of the FID detector
        0: Take the yield sum with a preference of the FID detector

    Returns
    -------
    h5filename : String
        The name of the created HDF5 file storage
    '''

    # Reading input conditions
    datapoints, gen_configs = reading_input_conditions(filename)

    # Store input conditions
    store_processing_conditions_hdf(datapoints, gen_configs)

    # Make HDF5 file name
    sublist = list(datapoints.groupby('Date').groups.keys())
    sublist.sort()
    date_to_append = sublist[0].replace('/', '_')
    experimental_name = gen_configs['EXPERIMENTAL_name']
    h5filename = '{}_{}'.format(experimental_name, date_to_append)

    # Reading and storing RGA data to HDF5 file Container
    rga_channels = reading_storing_RGA_data(h5filename)

    # Reading and storing GCxGC FID data to HDF5 file Container
    if 'GCxGC' in datapoints.columns:
        reading_storing_GCxGC_fid_data(h5filename, datapoints, gen_configs)

    if gen_configs['Process_SCD'] == 'TRUE':
        # Reading and storing GCxGC FID data to HDF5 file Container
        reading_storing_GCxGC_scd_data(h5filename, datapoints, gen_configs)

    if gen_configs['Process_NCD'] == 'TRUE':
        # Reading and storing GCxGC FID data to HDF5 file Container
        reading_storing_GCxGC_ncd_data(h5filename, datapoints, gen_configs)

    if gen_configs['Process_LOA'] == 'TRUE':
        # Reading and storing GCxGC FID data to HDF5 file Container
        reading_storing_LOA_data(h5filename, datapoints, gen_configs)

    # Reading and storing RGA calibration factors
    CFs = reading_storing_cal_factors(h5filename, gen_configs)

    # Reading GCxGC library (NOT STORED in HDF5 file)
    gcxgc_lib = read_gcxgc_lib(gen_configs)

    ######################################
    # SCD_NCD_detector
    ######################################
    # Look for a way to read in the
    # response factors and name to InChI
    ######################################

    # Processing the experimental data
    ######################################
    # SCD_NCD_detector
    ######################################
    # The main changes will be in "processing_yields"
    ######################################
    processing_yields(
        datapoints,
        rga_channels,
        CFs,
        gen_configs,
        gcxgc_lib,
        fid_sum,
        h5filename
    )

    # Processing logging files
    reading_storing_loggingfiles(datapoints, gen_configs, h5filename)

    # Extracting temperatureprofiles during injections
    if gen_configs["Reactor"] == "JSR":
        processing_loggingfiles_jsr(
            datapoints,
            gen_configs,
            h5filename,
        )
    else:
        processing_loggingfiles_2_temperature_pressure_profiles(
            datapoints,
            gen_configs,
            h5filename
        )

        # Extract the cracking experiment data for all logging files
        processing_loggingfiles_extract_cracking(
            h5filename,
            gen_configs["Reactor"]
        )

    # Reduce size of HDF5 file by repacking and compressing
    repack(h5filename)

    # Create excel files
    ######################################
    # SCD_NCD_detector
    ######################################
    # Check whether an excel summary is
    # created for the GCxGC-SCD/NCD data
    # at the end in 'create_excel_files'
    ######################################
    create_excel_files(h5filename, gcxgc_lib)

#     # Copy to designated folder
#     copy2database_folder(h5filename, gen_configs)

    return h5filename


if __name__ == '__main__':
    import sys

    print(sys.argv[1:])
    filename = sys.argv[1:]

    process_experimental_results(filename, fid_sum=1)
