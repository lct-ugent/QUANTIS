'''
Created on 21 Sep 2017

@author: shsymoen
'''

import sys

# Add some more folders to the Python path
sys.path.append('..')
sys.path.append('.')

# Local imports
from dataProcessing.utils import (
    add2sim_exp,
    check_diff_index,
    check_drive_available,
    check_file_available,
    conf_int_dist_chi2,
    diff_index,
    duplicate_remover,
    inchi_2_smiles,
    inchi_mol_wt,
)


class TestUtils(object):

    def test_check_drive_available(self):
        """ Automatic pytest for check_drive_available """
        import sys
        if sys.platform == "linux":
            assert check_drive_available(r"/home")
        elif 'win' in sys.platform:
            assert check_drive_available(r"C:\\")
        assert check_drive_available(r"bar") is False

    def test_conf_int_dist_chi2(self):
        """ Automatic pytest for check_drive_available """
        import pytest

        conf_int_sim = conf_int_dist_chi2(0.95, 2)
        conf_int_exp = 5.9914645471079799
        assert pytest.approx(conf_int_exp) == conf_int_sim

    def test_check_file_available(self):
        """ Automatic pytest for check_file_available """

        if sys.platform == "linux":
            import os
            open(r"tmp_python.txt", 'x')
            assert check_file_available("tmp_python.txt")
            os.remove("tmp_python.txt")
        elif 'win' in sys.platform:
            assert check_file_available(r'C:\Windows\explorer.exe')
        assert check_file_available(r'bar') is False

    def test_inchi_2_smiles(self):
        """ Automatic pytest for inchi_2_smiles """
        assert inchi_2_smiles("1S/C2H4/c1-2/h1-2H2") == "C=C"
        assert inchi_2_smiles("not_inchi") == "not_inchi"
        assert inchi_2_smiles("1S/C3H6/c1-3-2/h3H,1H2,2H3") == "C=CC"

    def test_inchi_mol_wt(self):
        """ Automatic pytest for inchi_mol_wt """
        import numpy as np

        np.testing.assert_almost_equal(
            inchi_mol_wt("1S/C2H4/c1-2/h1-2H2"),
            28.054
        )
        assert inchi_mol_wt("not_inchi") == 0

    def test_popup_entry(self):
        """ Automatic pytest for inchi_2_smiles """
        # not implemented yet
        assert 1

    def test_duplicate_remover(self):
        """ Automatic pytest for duplicate_remover """
        import pandas as pd

        # Create the input DataFrames
        dct = {
            i+1: [
                int('{}{}'.format(j+1, i+1)) for j in range(3)
            ]
            for i in range(3)
        }
        df_input = pd.DataFrame(dct, index=['name_1', 'name_2', 'name_1'])
        df_input3 = pd.DataFrame(dct, index=['name_1', 'name_2', 'name_3'])

        # Expected output DataFrames
        df_output3 = df_input3.copy()
        df_output_first = df_input.iloc[:2, :]
        df_output_last = df_input.iloc[1:, :]

        # Perform the tests
        pd.testing.assert_frame_equal(
            df_output3,
            duplicate_remover(df_input3)
        )
        pd.testing.assert_frame_equal(
            df_output_first,
            duplicate_remover(df_input, 'first')
        )
        pd.testing.assert_frame_equal(
            df_output_last,
            duplicate_remover(df_input, 'last')
        )

    def test_add2sim_exp(self):
        """ Automatic pytest for add2sim_exp """

        import pandas as pd

        # Input parameters for add2sim_exp module
        label = 'COILSIM1D_sim'
        temperaturesKelvin = [1073.15, 1098.15, 1123.15, 1123.15]

        input_df = pd.read_csv(
            r'testFiles_dataProcessing/add2_sim_df.csv',
            index_col=0
        )
        input_sim_exp = pd.read_csv(
            r'testFiles_dataProcessing/add2_sim_simexp_input.csv', index_col=0
        )

        # Read in the expected output of add2sim_exp
        expected_output = pd.read_csv(
            r'testFiles_dataProcessing/add2_sim_out.csv', index_col=0
        )

        # Get the output of the add2sim_exp module
        sim_output = add2sim_exp(
            input_df,
            label,
            temperaturesKelvin,
            input_sim_exp
        )

        # Perform the test
        pd.testing.assert_frame_equal(sim_output, expected_output)

    def test_diff_index(self):
        """ Automatic pytest for diff_index """
        import numpy as np
        import pandas as pd

        # Create first df
        df1 = pd.DataFrame(np.random.rand(3, 3), index=[1, 3, 5])
        # Create second df
        df2 = pd.DataFrame(np.random.rand(3, 3), index=[1, 4, 5])

        sim_out = diff_index(df1, df2)
        exp_out = pd.Index([3])

        pd.testing.assert_index_equal(sim_out, exp_out)

    def test_check_diff_index(self):
        """ Automatic pytest for check_diff_index """
        import numpy as np
        import pandas as pd

        # Create first df
        df1 = pd.DataFrame(np.random.rand(3, 3), index=[1, 3, 5])
        # Create second df
        df2 = pd.DataFrame(np.random.rand(3, 3), index=[1, 4, 5])

        assert check_diff_index(df1, df2) is False
        assert check_diff_index(df1, df1)
