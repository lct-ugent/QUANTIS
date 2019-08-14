'''
Created on 11 Aug 2017

@author: shsymoen
'''


class metadata_from_hdf5():
    '''
    Metadata class
    '''

    def __init__(self, store):
        ''' Constructor '''

        self.tags = store.root._v_attrs._f_list()
        self.dct = {
            k: store.root._v_attrs[k]
            for k in self.tags
            }
        self.__repr__()

    def __str__(self):
        """
        Shows the metadata from a certain HDF5 file
        """

        returner = ''
        for tag in self.tags:
            tagr = tag.replace(
                    '_gh',
                    ' (g/h)'
                ).replace(
                    '_bara',
                    ' (bara)'
                ).replace(
                    '_C',
                    ' (°C)'
                ).replace(
                    '_1_to_5',
                    ' (1-5)'
                ).replace(
                    'AAA_',
                    ''
                )
            returner += '{}:\n\t{}\n\n'.format(tagr, self.dct[tag])
        return returner[:-1]

    def __repr__(self):
        """ Shows a representation of the metadata """
        return "{}".format(self.tags)


class configurations_from_hdf5():
    '''
    Configurations class
    '''

    def __init__(self, store):
        ''' Constructor '''

        self.conditions = store['processing/conditions']
        self.gen_configs = store['processing/gen_configs']

    def write_conditions(self, file_name='conditions.csv'):
        '''Writes the conditions csv from a DHF5 store

        Parameters
        ----------
        self : Configurations_from_hdf5

        Returns
        -------
        None

        Creates
        -------
        csv file
        '''
        self.conditions.to_csv(file_name)
        return None

    def write_gen_configs(self, file_name='gen_configs.csv'):
        '''Writes the gen_configs csv from a DHF5 store

        Parameters
        ----------
        self : Configurations_from_hdf5

        Returns
        -------
        None

        Creates
        -------
        csv file
        '''
        self.gen_configs.to_csv(
            file_name,
            header=None,
            index=None
            )
        return None
