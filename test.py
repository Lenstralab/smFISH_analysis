#!/usr/bin/env python3

import sys, os, yaml
import smfish
from smfish.FISH_pipeline import FISH_pipeline, get_exp_info, get_filelist
import tllab_common
from tllab_common.misc import getConfig
from tllab_common.wimread import imread

fname = os.path.realpath(__file__)
test_files = os.path.join(os.path.dirname(fname), 'test_files')

# This file defines tests to be run to assert the correct working of our scripts
# after updates. Add a test below as a function, name starting with 'test', and 
# optionally using 'assert'.
#
# Place extra files used for these tests in the folder test_files, add imports 
# above this text.
#
# Then navigate to the directory containing this file and run ./test.py directly
# from the terminal. If you see red text then something is wrong and you need to
# fix the code before committing to gitlab.
#
#wp@tl20200124


def make_test_FISH_pipeline(parameter_file):
    def pipeline_fun(tmp_path):
        outputfolder = str(tmp_path)
        analysisfolder = os.path.join(outputfolder, 'analysis')

        parameters = getConfig(os.path.join(test_files, parameter_file))
        parameters['outputfolder'] = outputfolder
        parameters['analysisfolder'] = analysisfolder
        tmp_parameter_file = os.path.join(outputfolder, 'params.yml')

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        with open(tmp_parameter_file, 'w') as f:
            yaml.dump(parameters, f)

        parameters, rawdistributions = FISH_pipeline(tmp_parameter_file)

        files = ('_max.tif', '_max_nucleus_mask.tif', '_max_cell_mask.tif', '_loc_results_cy3.tif',  '_mask_cell+nuc+spots.tif', '_TS.tif')

        for file in range(0, parameters['lenfileListIn']):
            tmp_paramsfile = get_exp_info(parameters['fileListIn'][file], parameters)
            tmp_lfn = get_filelist(parameters['fileListIn'][file], parameters)

            assert os.path.exists(tmp_paramsfile['pathOut']), 'Output folder has not been generated'
            assert len(os.listdir(tmp_paramsfile['pathOut']))>4, 'There aren''t enough files in the output folder'
            for fn in range(0, len(tmp_lfn)):
                for file in files:
                    assert os.path.exists(tmp_paramsfile['pathOut'] + tmp_paramsfile['date'] + "_" + tmp_lfn[fn] + file),\
                            'File {} has not been generated'.format(tmp_paramsfile['pathOut'] + tmp_paramsfile['date'] + "_" + tmp_lfn[fn] + file)

        assert rawdistributions, 'rawdistributions not generated'
    return pipeline_fun

test_demo  = make_test_FISH_pipeline('FISH_pipeline_test.yml')

## ----- This part runs the tests -----
    
if __name__ == '__main__':
    if len(sys.argv)<2:
        py = ['3.8']
    else:
        py = sys.argv[1:]

    for p in py:
        print('Testing using python {}'.format(p))
        print('Loaded smfish code from {}'.format(smfish.__file__))
        print('Loaded tllab_common code from {}'.format(tllab_common.__file__))
        os.system('python{} -m pytest -n=12 -p no:warnings --verbose {}'.format(p, fname))
        print('')

        imread.kill_vm()
