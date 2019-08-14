'''
Created on 12 Oct 2016

@author: shsymoen
'''

import pandas as pd


def convert_smiles_sort(df):
    '''Converts the column names (InChI) of a dataframe to SMILES and sorts
    them based on the number of C-atoms

    Parameters
    ----------
    df : DataFrame
        column : InChI names
        index : injection tags

    Returns
    -------
    sort_df : DataFrame
        Sorted based on C-atoms
    '''
    df_copy = df.copy()
    # Convert InChI names to Smiles
    df_copy.columns = [
        inchi_2_smiles(inchi_name)
        for inchi_name
        in df_copy.columns
    ]
    # Sort components based on their carbon number
    df_copy.loc['Catoms', :] = df_copy.columns.map(
        lambda x: str(x).count('C') + str(x).count('c')
    )  # RDKIT has sometimes small c's in smiles name
    sort_df = df_copy.sort_values('Catoms', axis=1).drop(labels="Catoms", axis=0)

    return sort_df


def duplicate_remover(dataframe, first_last='first'):
    '''Removes the duplicated lines in a dataframe based on the index

    Extended description module

    Parameters
    ----------
    dataframe : DataFrame
        A dataframe with duplicated indices
    first_last : String
        Describes which indices to keep (the first occurrence or the last)

    Returns
    -------
    frame : DataFrame
        the dataframe that was given as input is returned with the removal
        of the duplicated indices.

    Examples
    --------
    >>> duplicate_remover(dataframe, first_last='first')

    '''
    return dataframe[~dataframe.index.duplicated(keep=first_last)]


def add2sim_exp(df, label, temperaturesKelvin, sim_exp):
    """
    Add simulation to overview DataFrame

    Parameters
    ----------
    df : DataFrame
        Data to be added to the overview
    label : String
        A label that describes the set of data
        that is added to the sim_exp file
    temperaturesKelvin : list
        Describe
    sim_exp : DataFrame
        Pre-existing data where the new data will be appended to

    Returns
    -------
    df : DataFrame
        Newly created dataframe with new data concatenated to the old one

    Examples
    --------
    To be added
    """

    df_transpose = df.transpose()
    df_transpose['exp_sim'] = label
    df_transpose['Temperature'] = temperaturesKelvin
    df_transpose.reset_index(drop=True, inplace=True)
    return pd.concat([sim_exp, df_transpose], ignore_index=True, sort=True)


def deconv_data(data):
    """Transforms an overall DataFrame to smaller DataFrames

    DEPRECATED

    Parameters
    ----------
    data : DataFrame

    Returns
    -------
    feedflow :
    temperaturesKelvin : List
    reactants_results :
    feedName_InChI : String

    """
    feedFlow = data['Results']['inputConditions'][1:-1]
    temperaturesCelsius = data['Results']['inputConditions'].loc[
        'Temperature', :
    ].tolist()
    temperaturesKelvin = [T + 273.15 for T in temperaturesCelsius]
    reactants_results = [
        reactant.replace(' g/hr', '')
        for reactant in feedFlow.index
    ]
    feedName_InChI = data['CHO']['inputConditions'].ix[
        reactants_results, 'InChI'].tolist()

    return feedFlow, temperaturesKelvin, reactants_results, feedName_InChI


def data2results(data, sheetname='Mol% wet'):
    '''Convert data structure to an easier to handle results dataframe

    DEPRECATED

    '''

    # Get all experimental data
    df_results = pd.concat([data[sheetname]['RGA_TCD'],
                            data[sheetname]['RGA_FID'],
                            data[sheetname]['GCxGC'],
                            data[sheetname]['LOA']
                            ])

    if 'wet' in sheetname and 'H2O' not in df_results.index:
        water = data[sheetname]['inputConditions']
        water = water.loc['H2O g/hr', :]
        water.name = 'H2O'
        df_results = df_results.append(water)

    # Clean up results dataframe
    df_results.dropna(inplace=True)

    # Create results converter file
    sheet_converter = 'CHO'
    results_converter = pd.concat(
        [data[sheet_converter]['RGA_TCD'],
         data[sheet_converter]['RGA_FID'],
         data[sheet_converter]['GCxGC'],
         data[sheet_converter]['LOA']]
    )

    if 'wet' in sheetname and 'H2O' not in results_converter.index:
        water = data[sheet_converter]['inputConditions']
        water = water.loc['H2O', :]
        water.name = 'H2O'
        results_converter = results_converter.append(water)

    # Remove duplicates
    results_converter.dropna(inplace=True)
    results_converter.drop('Name', axis=0, inplace=True)

    # Remove duplicate indices take components
    # from TCD > RGA > GCxGC (depending on order in GC_processor summary file)
    results_converter.drop_duplicates(inplace=True, keep='first')

    # Change index of results dataframe to InChI's
    df_results.index = results_converter.loc[
        df_results.index, 'InChI'
    ].tolist()

    # Remove duplicate indices (based on InChI)
    # take components from TCD > RGA > GCxGC
    df_results = duplicate_remover(df_results)
    df_results.index.name = 'InChI'

    return df_results, results_converter


def check_drive_available(pathname):
    """
    Checks if a certain pathname is available.

    Parameters
    ----------
    pathname : String
        A string of the path location that needs to be tested

    Returns
    -------
    flag : boolean

    Examples
    --------
    >>> check_drive_available(r'C:\')
    True
    >>> check_drive_available(r'bar')
    False
    """
    import os

    return os.path.isdir(pathname)


def check_file_available(pathname):
    """
    Checks if a certain file is available.

    Parameters
    ----------
    pathname : String
        A string of the file location that needs to be tested

    Returns
    -------
    flag : boolean

    Examples
    --------
    >>> check_file_available(r'C:\Windows\explorer.exe')
    True
    >>> check_file_available(r'bar')
    False
    """
    import os

    return os.path.exists(pathname)


def popup_error(text):
    """Creates a popup error window with the inserted text

    The tkinter python package is employed to show popup error messages.

    Parameters
    ----------
    text : String
        The text that has to be shown in the error window

    Returns
    -------
    None

    Creates
    -------
    Popup error window with inserted text
    """
    import tkinter as tk
    from tkinter.messagebox import showerror

    root = tk.Tk()
    root.withdraw()
    showerror(message=text)

    return None


def popup_warning(text):
    """Creates a popup error window with the inserted text

    The tkinter python package is employed to show popup warning messages.

    Parameters
    ----------
    text : String
        The text that has to be shown in the error window

    Returns
    -------
    None

    Creates
    -------
    Popup error window with inserted text
    """
    import tkinter as tk
    from tkinter.messagebox import showwarning

    root = tk.Tk()
    root.withdraw()
    showwarning(message=text)

    return None


def popup_entry(text):
    """
    Creates a popup window to ask for additional input
    ---------
    Returns: string with additional input
    """
    import tkinter as tk
    from tkinter import simpledialog

    root = tk.Tk()
    root.withdraw()
    answer = simpledialog.askstring('Additional input', text, parent=root)
    root.destroy()

    return answer


def inchi_2_smiles(inchi_name):
    '''Returns the smiles name of the given InChI name

    More description

    Parameters
    ----------
    inchi_name : String
        InChI identifier

    Returns
    -------
    smiles : String
        Smiles representation of the molecule that was
        given as input with its InChI identifier.
        If the Smiles representation can not be found,
        the given input is returned.

    Examples
    --------
    >>> inchi_2_smiles("1S/C2H4/c1-2/h1-2H2")
    C=C
    >>> inchi_2_smiles("not_inchi")
    not_inchi
    >>> inchi_2_smiles("1S/C3H6/c1-3-2/h3H,1H2,2H3")
    C=CC
    '''
    from rdkit import Chem

    with suppress_stdout_stderr():
        m = Chem.MolFromInchi('InChI={}'.format(inchi_name))
    try:
        smiles = Chem.MolToSmiles(m, isomericSmiles=True)
    except:
        smiles = inchi_name

    return smiles


def smiles_2_inchi(smiles_name):
    '''Returns the InChI name of the given Smiles name

    More description

    Parameters
    ----------
    smiles_name : String
        Smiles identifier

    Returns
    -------
    inchi : String
        InChI representation of the molecule that was
        given as input with its Smiles identifier.
        If the InChI representation can not be found,
        the given input is returned.

    Examples
    --------
    >>> smiles_2_inchi("C=C")
    1S/C2H4/c1-2/h1-2H2
    >>> smiles_2_inchi("not_inchi")
    not_inchi
    >>> smiles_2_inchi("C=CC")
    1S/C3H6/c1-3-2/h3H,1H2,2H3
    '''
    from rdkit import Chem

    with suppress_stdout_stderr():
        m = Chem.MolFromSmiles(smiles_name)
    try:
        inchi = Chem.MolToInchi(m)
    except:
        inchi = smiles_name

    return inchi


def inchi_mol_wt(inchi_name):
    ''' Returns the molecular weight of the molecule

    Parameters
    ----------
    inchi_name : String
        InChI identifier of the component

    Returns
    -------
    mw : float
        Molecular weight of the molecule

    Examples
    --------
    >>> inchi_2_smiles("1S/C2H4/c1-2/h1-2H2")
    28.05
    '''
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors

    MolWt = lambda *x, **y: _rdMolDescriptors._CalcMolWt(*x, **y)

    mw = 0

    with suppress_stdout_stderr():
        m = Chem.MolFromInchi('InChI={}'.format(inchi_name))

    try:
        mw = MolWt(m)
    except:
        pass

    return mw


def diff_index(df1, df2):
    '''Computes an index with the differences between the index df1 and df2

    Extended description

    Parameters
    ----------
    df1 : DataFrame
    df2 : DataFrame

    Returns
    -------
    difference : Index
    '''

    return df1.index.difference(df2.index)


def check_diff_index(df1, df2):
    '''Returns a boolean if df1 has indices not in df2

    Extended description

    Parameters
    ----------
    df1 : DataFrame
    df2 : DataFrame

    Returns
    -------
    different : boolean
    '''

    diff_ind = diff_index(df1, df2)
    return diff_ind.empty


def conf_int_dist_chi2(conf_int, degrees):
    '''Calculates the confidence interval distance
    for a given confidence interval.

    Parameters
    ----------
    conf_int : float
        The confidence level
    degrees : int
        Number of degrees

    Returns
    -------
    conf_int_dst : float
        The distance for the given confidence level

    Examples
    --------
    >>> conf_int_dist_chi2(0.95, 2)
    5.9914645471079799
    '''
    import scipy.stats as stats

    conf_int_dst = stats.chi2.ppf(
        conf_int,
        degrees,
    )

    return conf_int_dst


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''

    def __init__(self):
        import os
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        import os
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        import os
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


def add_excel_formatting(writer, nuumber_formatting='0.0#'):
    '''Adds formatting to a given excel spreadsheet

    Adds number and cell formatting to a Pandas Excel writer object
    for every sheet available in the writer.
    It further reduces the zoom of the sheet.

    Parameters
    ----------
    writer : Pandas Excel writer object
    number_formatting : String

    Returns
    -------
    writer : Pandas Excel writer object
    '''
    # Get access to the workbook
    workbook = writer.book

    # Add formatting
    format1 = workbook.add_format({
            'align': 'center',
            'num_format': '0.0#'
    })

    # Apply formatting to the excel sheet
    for sheetname in writer.sheets:
        worksheet = writer.sheets[sheetname]
        worksheet.set_column('A:XFD', 10, format1)
        # Reduce the zoom of the excel file
        worksheet.set_zoom(90)

    return writer


def storing_2_HDF5_store(df, store, key, data_column=True):
    '''Stores the DataFrame to the HDF5 store

    Extended description

    Parameters
    ----------
    df : DataFrame
        The DataFrame that needs to be added to the HDF5 store
    store : HDF5 store
        Opened via pd.HDFStore(filename)
    key : String
        The key to which the DataFrame should be stored

    Returns
    -------
    None

    Examples
    --------
    >>> storing_2_HDF5_store(df, store, key)

    '''
    df.to_hdf(
        store,
        key,
        format='t',
        mode='a',
        encoding='utf-8',
        data_column=data_column
    )

    return None


if __name__ == '__main__':
    import doctest
    doctest.testmod()
