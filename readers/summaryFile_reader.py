'''
Created on 11 Oct 2016

@author: shsymoen
'''

import pandas as pd


def get_excelResultsFile(results_file_name):
    ''' Reading excel file created by GC_processor '''

    # Results sheet
    xl = pd.ExcelFile(results_file_name)
    results = pd.read_excel(xl, sheetname='Results', index_col=0, header=None)
    results.columns = [
        'Exp_{:0>3d}'.format(i) for i in range(len(results.columns))
    ]

    # Find indices
    index_rga_tcd = results.index.get_loc('RGA TCD Results')
    index_rga_fid = results.index.get_loc('RGA FID Results')
    index_gcxgc = results.index.get_loc('GCxGC Results')
    index_sum = results.index.get_loc("C4-")
    try:
        index_loa = results.index.get_loc('LOA Results')
    except:
        index_loa = index_sum
    index_geometry = results.index.get_loc("Reactor length (m)")
    index_temperatureProfile = results.index.get_loc(
        "Axial position reactor (mm)"
    )
    index_pressures = results.index.get_loc("CIP (abar)")

    indices = [
               (0, index_rga_tcd - 1),
               (index_rga_tcd, index_rga_fid - 1),
               (index_rga_fid, index_gcxgc - 1),
               (index_gcxgc, index_loa - 1),
               (index_loa, index_sum - 1),
               (index_sum, index_sum + 3)
               ]

    data = {}
    for sheetname in xl.sheet_names:
        df_data = pd.read_excel(
            xl,
            sheetname=sheetname,
            index_col=0,
            header=None
        )
        data[sheetname] = read_sheet_results(df_data, indices, sheetname)

    geometry = results.iloc[
        index_geometry: index_temperatureProfile - 1, :1
    ].transpose()
    temperatureProfiles = results.iloc[
        index_temperatureProfile + 1: index_pressures - 1, :
    ]
    temperatureProfiles.index = [
        int(pos) / 1000 for pos in temperatureProfiles.index
    ]
    pressures = results.iloc[index_pressures: index_pressures + 2]
    pressures.index = [0, 1.43]

    # Extract feed flows from conditions
    geometry.columns = ['length (m)', 'diameter (m)']
    geometry['thickness (m)'] = 0.002

    return temperatureProfiles, pressures, geometry, data


def read_sheet_results(df_sheet, indices, sheetname):
    '''
    Read sheet results of GC_processor results.xlsx file
    '''

    dict_names = [
        'inputConditions',
        'RGA_TCD',
        'RGA_FID',
        'GCxGC',
        'LOA',
        'sums'
    ]
    data = {}
    for i in range(len(indices)):
        k, e = indices[i]
        data[dict_names[i]] = df_sheet.iloc[k:e, :]
        if sheetname != 'CHO':
            data[dict_names[i]].columns = [
                'Exp_{:0>3d}'.format(i)
                for i
                in range(len(data[dict_names[i]].columns))
            ]
        else:
            if dict_names[i] == 'inputConditions':
                new_header = data[dict_names[i]].loc['Name', :]
                data[dict_names[i]] = data[dict_names[i]].iloc[1:, :]
            data[dict_names[i]].columns = new_header

    return data
