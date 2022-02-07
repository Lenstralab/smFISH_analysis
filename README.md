# smFISH
Lenstra labs code for single molecule FISH as used by Brouwer et al. This script uses Python 3.

## Installation
Do this the first time only:

    git clone https://github.com/Lenstralab/smFISH_analysis.git
    pip install -e smFISH_analysis --user

This will install the smfish package in your personal path in editable mode.

## Running the pipeline
From a terminal:

    FISH_pipeline /path/to/parameter_file.yml

Note that because of the pip step in the installation the FISH_pipeline script is in your path
and accessible from anywhere. So you don't have to specify the path. If this is not what you
want, you can still call the script by path, just notice it's in the subfolder 'smfish':

    /path/to/repo/smfish/FISH_pipeline.py /path/to/parameter_file.yml

Or from an ipython window:

    %run /path/to/repo/smfish/FISH_pipeline.py /path/to/parameter_file.yml

