import sys
sys.path.append('.')
sys.path.append('..')

from dataProcessing.GC_processor_HDF5 import (
    processing_chrom_gcxgc_fid_scd_ncd,
    processing_chrom_loa,
    processing_chrom_rga,
    processing_chroms_2_raw,
    fid_tcd_merge,
    read_gcxgc_lib,
    reading_input_conditions,
    reading_storing_loggingfiles,
    reading_storing_GCxGC_fid_data,
    reading_storing_RGA_data,
    reading_storing_cal_factors,
    yield_data,
)


if __name__ == '__main__':

    import os

    # HDF5 file name
    h5filename = r"JSR_Al_alloy_1CC_14_06_2017"

    # Read the conditions and general configurations file
    filename = "conditions_JSR.csv gen_configs_JSR.csv".split(' ')
#     filename = [
#         r'testFiles_dataProcessing/yield_processor_test/JSR/{}'.format(file)
#         for file
#         in filename
#     ]
    conditions, gen_configs = reading_input_conditions(filename)

    # RGA subdivided locations (should be added in gen_configs)
    tcd_location = \
        r'TCD/20170614-Al alloy 1st CC'
    fid_location = \
        r'FID/20170614-Al alloy 1st CC'
    rga_location = \
        r'RGA'

    foldername = fid_tcd_merge(tcd_location, fid_location, rga_location)

    # Reading and storing RGA calibration factors
    CFs = reading_storing_cal_factors(h5filename, gen_configs)

    # Reading GCxGC library (NOT STORED in HDF5 file)
    gcxgc_lib = read_gcxgc_lib(gen_configs)

    # Read JSR RGA data
    rga_channels = reading_storing_RGA_data(h5filename, gen_configs)

    # Read C5+ JSR data
    reading_storing_GCxGC_fid_data(
        h5filename,
        conditions,
        gen_configs,
        detector_switch='FID'
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
    fid_sum = False
    processing_chroms_2_raw(
        exp_data,
        rga_channels,
        gcxgc_lib,
        fid_sum,
        h5filename
    )

    reading_storing_loggingfiles(conditions, gen_configs, h5filename)

    # Remove the created HDF5 file
    os.remove('{}.h5'.format(h5filename))
    print('\n{}.h5 removed.'.format(h5filename))
