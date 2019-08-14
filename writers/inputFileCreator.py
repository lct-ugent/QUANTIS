'''
Created on 11 Oct 2016

@author: shsymoen
'''

from dataProcessing.utils import popup_error
import numpy as np
import os
import pandas as pd
import xlsxwriter


def write_chemkinParInputFile(fileDir, inputData):
    ''' Write CHEMKIN parallelizer input files'''

    ########################################
    # Enhancement
    ########################################
    # Remove MW, reactants -> store reactants together with flows
    # Use RDKit to compute MW
    ########################################

    # Collect data
    reactor_dimensions, reactants, MW, temperatureProfiles, \
        axialTempPos, pressures, flows = inputData

    # Open excel workbook
    wbk = xlsxwriter.Workbook(fileDir)
    sheet = wbk.add_worksheet("reactor_input")

    # Write standard labels
    labels = ["REACTOR", "SPCS", "MW", "RUN"]
    sheet.write_column(0, 0, labels)

    # Write reactor dimensions
    sheet.write(0, 1, reactor_dimensions[0])
    sheet.write(0, 2, reactor_dimensions[1])

    # Write species and molecular weight
    for i, specie in enumerate(reactants):
        sheet.write(1, i + 1, specie)
        sheet.write(2, i + 1, MW[i])
        sheet.write(3, i + 1, specie)

    # Write more labels
    axialPresPos = [0, axialTempPos[-1]]
    numberOfSpecies = len(reactants)
    numberOfAxialPos = len(temperatureProfiles[0])
    sheet.write(3, 1+numberOfSpecies, 'TEMP')
    sheet.write_row(3, 1+numberOfSpecies+1, axialTempPos)  # convert to m
    sheet.write(3, 1 + numberOfSpecies + 1 + numberOfAxialPos, 'PRESS')
    sheet.write_row(
        3, 1 + numberOfSpecies + 1 + numberOfAxialPos + 1, axialPresPos
    )
    sheet.write_row(
        3, 1 + numberOfSpecies + 1 + numberOfAxialPos + 3, ['ATOL', 'RTOL']
    )

    # Write inputFlows, temperature profiles
    # and inlet and outlet pressure and tolerances
    sheet.write_column(4, 0, range(1, len(temperatureProfiles)+1))
    for j, temperatureProfile in enumerate(temperatureProfiles):
        sheet.write_row(4+j, 1, flows[j])
        sheet.write_row(4+j, 1 + numberOfSpecies + 1, temperatureProfile)
        sheet.write_row(
            4+j, 1 + numberOfSpecies + 1 + numberOfAxialPos + 1, pressures[j]
        )
        sheet.write_row(
            4+j,
            1 + numberOfSpecies + 1 + numberOfAxialPos + 3,
            ['1E-09', '0.000001']
        )

    wbk.close()


def excel2csv(fileName):
    ''' Transforms excel file to csv file '''

    df = pd.read_excel(fileName, header=None)
    fileDir = os.path.dirname(fileName)
    fileDirName = os.path.join(fileDir, 'reactor_input.csv')
    df.to_csv(fileDirName, header=None, index=None)


def write_coilsimParInputFiles(fileDir, inputData, reactor, mol, v39):
    ''' Write COILSIM1D parallelizer input files'''

    import pandas as pd
    import numpy as np
    import sys

    coilsimFeedName, temperatureProfiles, pressures, \
        feedFlow, geometry = inputData

    if np.nan is coilsimFeedName:
        popup_error(
            "Not all components found in the COILSIM1D translator network.\n\n"
            "Program terminated unsuccessfully!\ncheckflag for writing"
            "COILSIM1DPar input files should be switched off."
            )
        sys.exit()

    writer_conditions = pd.ExcelWriter(
        os.path.join(fileDir, 'conditions.xlsx')
    )

    # Last temperature point can't be run with COILSIM1D
    # because the heat flux drop is too big for BSSC reactor.
    if reactor == 'BSSC_metal':
        temperatureProfiles.iloc[:-1, :].to_excel(
            writer_conditions,
            sheet_name='temperatureProfiles'
        )
    else:
        temperatureProfiles.to_excel(
            writer_conditions,
            sheet_name='temperatureProfiles'
        )
    pressures.to_excel(writer_conditions, sheet_name='pressureProfiles')
    feedFlow.index = ['{} flowrate'.format(coilsimFeedName), 'water flowrate']
    feedFlow.to_excel(writer_conditions, sheet_name='feedFlow')

    pd.DataFrame(
        data=[[0, 0, 1, mol, v39]],
        columns=['calculateProfit', 'EkviCalc', 'Plotting', 'mol', 'v3.9']
    ).to_excel(
        os.path.join(fileDir, 'parameters.xlsx'),
        sheet_name='parameters'
    )
    # Saving conditions file for COILSIM1D parallellizer
    writer_conditions.save()

    # Saving geometry file
    geometry.to_excel(
        os.path.join(fileDir, 'geometry.xlsx'),
        sheet_name='geometry',
        index=None
    )


def storeMultipleCHEMKINpar_input(chemkin_names,
                                  geometry,
                                  chemkinConverters,
                                  feedName_InChI,
                                  MW,
                                  temperatureProfiles,
                                  pressures,
                                  feedFlow):
    '''
    Iterates over all CHEMKIN names and creates input files
    from the temperature and pressureprofiles of the experiment
    '''
    for chemkin_name in chemkin_names:

        storeSingleCHEMKINpar_input(
            chemkin_name,
            geometry,
            chemkinConverters,
            feedName_InChI,
            MW,
            temperatureProfiles,
            pressures,
            feedFlow
            )


def storeSingleCHEMKINpar_input(chemkin_name,
                                geometry,
                                chemkinConverters,
                                feedName_InChI,
                                MW,
                                temperatureProfiles,
                                pressures,
                                feedFlow
                                ):
    '''Creates input files for the CHEMKIN parallellizer
    for a given chemistry network name.
    The network should include the inflow components

    Parameters
    ----------
    chemkin_name

    Returns
    -------
    CHEMKINParallellizer input files to use with the automatic PFR plugin
    '''
    import os
    import sys

    # Create output folder
    if not os.path.exists(r'simulationInputFiles\{}'.format(chemkin_name)):
        os.makedirs(r'simulationInputFiles\{}'.format(chemkin_name))

    # Rearrange the data
    CP_reactorcond = geometry.values.tolist()[0][:-1]

    CP_reactants = chemkinConverters[chemkin_name].loc[
        feedName_InChI,
        'simulation names'
    ].tolist()
    print(feedName_InChI)
    if np.nan in CP_reactants:
        popup_error(
            "Not all components found in the CHEMKIN network "
            "{}.\n\nProgram terminated unsuccessfully!"
            "\n\nOther network should be chosen or\ncheckflag for writing "
            "CHEMKIN input files should be switched off.".format(chemkin_name)
            )
        sys.exit()

    CP_MW = MW
    CP_Tprofiles = temperatureProfiles.transpose().as_matrix().tolist()
    CP_axialpos = temperatureProfiles.index.tolist()
    CP_pressureprofiles = pressures.transpose().as_matrix().tolist()
    CP_flows = feedFlow.as_matrix().tolist()

    inputData = [
        CP_reactorcond,
        CP_reactants,
        CP_MW,
        CP_Tprofiles,
        CP_axialpos,
        CP_pressureprofiles,
        CP_flows
    ]

    # Store the files in the designated folder
    fileName = r'simulationInputFiles\{}\reactor_input.xlsx'.format(
        chemkin_name
    )
    write_chemkinParInputFile(fileName, inputData)
    excel2csv(fileName)
    os.remove(fileName)


def storeCOILSIM1Dpar_input(fileDir, inputData, reactor, mol, v39=0):
    '''
    Check if the folder exists, rearrange data and store
    input files for the COILSIM1Dparallellizer
    '''
    import os

    # Deconvolute inputData
    TProfiles, pressures, feedFlow, geometry, coilsimConverter, \
        feedName_InChI = inputData

    coilsimFeedName = None

    if not os.path.exists(fileDir):
        os.makedirs(fileDir)

    if feedName_InChI in coilsimConverter.index:
        coilsimFeedName = coilsimConverter.loc[feedName_InChI, 'name COILSIM']
    # If component not in translator just try feedName_InChI
    else:
        # Some imports
        import glob
        # Remember current working directory
        cwd = os.getcwd()
        os.chdir(r'N:\Project\918-Experimental set-ups\Feed analyses')
        feed_analyses = [
            analysis.replace('.xlsx', '') for analysis in glob.glob('*.xlsx')
        ]
        if feedName_InChI in feed_analyses:
            coilsimFeedName = feedName_InChI
        # Change directory back
        os.chdir(cwd)

    inputData = coilsimFeedName, TProfiles, pressures, feedFlow, geometry

    write_coilsimParInputFiles(fileDir, inputData, reactor, mol, v39)
