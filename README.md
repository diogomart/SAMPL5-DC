# SAMPL5-DC

Estimates the following:
1. free energy of hydration
2. free energy of solvation in cyclohexane
3. cyclohexane/water partition coefficient

## Requirements / dependencies
* Works on Linux
* Python 
* openbabel (v2.3.2 available through apt-get)
* python-openbabel (available through apt-get)
* MSMS (http://mgltools.scripps.edu/downloads)

## Setup
1. Download the script (empiricalSES.py)
2. Install the dependencies
3. Edit line X in empiricalSES.py to match the path to MSMS binary

## Usage
python /path/to/empiricalSES.py molecule.pdb

